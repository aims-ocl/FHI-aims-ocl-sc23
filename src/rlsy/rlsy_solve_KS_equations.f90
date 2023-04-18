module rlsy_solve_KS_equations
!!
!! Solve KS equations perhaps? Yes no maybe!
!!
use rlsy_constants, only: r8,i8,rl_huge,rl_hugeint,rl_tol,rl_sqtol,rl_twopi,&
    rl_exitcode_symmetry,rl_exitcode_param,rl_exitcode_memory,rl_iou
use rlsy_memtracker, only: rl_memtracker
use rlsy_timer, only: rl_timer_eigenvalue
use rlsy_helpers, only: tochar,rl_sqnorm,rl_chop,tochar,rl_mom_real,rl_mom_complex
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_distancetable, only: rl_distancetable
use rlsy_extended_cluster, only: rl_extended_cluster
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_verletlist, only: rl_verletbox
use rlsy_spacegroup, only: rl_spacegroup,rl_spacegroup_operation
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list
use rlsy_kpointmesh, only: rl_kpoint_mesh
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc,&
    rl_kspace_eigenproblem_prob,rl_kspace_eigenproblem_prob_real,&
    rl_kspace_eigenproblem_prob_complex,rl_kspace_eigenproblem_kpoint_real,&
    rl_kspace_eigenproblem_kpoint_complex
use rlsy_realspace_matrix, only: rl_realspace_matrix,rl_realspace_matrix_notdistributed,&
    rl_realspace_matrix_mediumdistributed,rl_realspace_matrix_fulldistributed
use elsi, only: elsi_ev_real,elsi_ev_complex,elsi_compute_mu_and_occ,&
    elsi_dm_real,elsi_dm_complex,elsi_compute_entropy,&
    elsi_get_mu,elsi_get_entropy,elsi_get_eval,elsi_get_evec_real,elsi_get_evec_complex

implicit none
private

public :: rl_solve_KS

contains

!> Takes all the fancy objects and solves the KS equations on a k-point mesh
subroutine rl_solve_KS(rmtx,KS,p,basis,kmesh,sym,mw,mem,verbosity,tmr)
    !> distributed matrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> eigenproblem handle
    class(rl_kspace_eigenproblem), intent(inout) :: KS
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> k-point mesh
    type(rl_kpoint_mesh), intent(in) :: kmesh
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timing information
    type(rl_timer_eigenvalue), intent(inout) :: tmr

    real(r8) :: timer

    timer=mpi_wtime()
    if ( verbosity .gt. 0 ) then
        write(rl_iou,*) ''
        write(rl_iou,*) 'SOLVING KS EQUATIONS'
    endif

    ! This routine takes the realspace Hamiltonian and overlap, stored in rmtx,
    ! fourier-transforms it to k-space and solves the Kohn-Sham eigenvalue problem.
    ! While I have a bunch of auxiliary information availabe, I calculate the
    ! chemical potential, occupation numbers and the realspace density matrix.
    !
    ! It is slightly confusing, since depending on how we have distributed things,
    ! it will flow slightly differently. I've tried to keep everything as coherent
    ! and similar as I could, but it went so-so.
    !
    ! I send this along to another routine. This is a little ugly, but it
    ! avoids tons of select type() everywhere, also the flow for single-proc
    ! and multiproc is so different it makes sense to have them as separate
    ! routines. But it's neat to only expose the parent routine, so that from
    ! a developer point of view, there is no doubt which routine is supposed
    ! to be called when solving the KS equations.
    select type(KS)
    type is(rl_kspace_eigenproblem_singleproc)
        call solve_KS_singleproc(KS,rmtx,basis,p,sym,kmesh,mw,mem,verbosity,tmr)
    type is(rl_kspace_eigenproblem_multiproc)
        call solve_KS_multiproc(KS,rmtx,basis,p,sym,kmesh,mw,mem,verbosity,tmr)
    end select

    ! And keep track of the total time
    tmr%total=tmr%total+mpi_wtime()-timer
end subroutine

subroutine solve_KS_multiproc(KS,rmtx,basis,p,sym,kmesh,mw,mem,verbosity,tmr)
    !> KS solver handle. The multiproc version.
    type(rl_kspace_eigenproblem_multiproc), intent(inout) :: KS
    !> distributed matrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> k-point mesh
    type(rl_kpoint_mesh), intent(in) :: kmesh
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timing information
    type(rl_timer_eigenvalue), intent(inout) :: tmr

    ! 50MB? That should not be too bad. Should be decided on-the-fly.
    real(r8), dimension(8) :: mellantid
    real(r8), parameter :: safemem=50*1024*1024
    real(r8) :: timer,t0,t1

    !ftr: block
        type(rl_mom_real), dimension(:,:), allocatable :: hbuf,obuf
        complex(r8), dimension(:,:,:), allocatable :: commbuf_ovl,commbuf_ham
        real(r8), dimension(3) :: kvector,rvector
        real(r8) :: kdotr,f0
        complex(r8) :: expikr
        integer :: clo,chi,coff,l
        integer :: ispin,ipair,jpair,iop,ii,jj
        integer :: b1,b2,l1,l2
        integer :: ia1,ia2,ja1,ja2,is1,is2,js1,js2
        integer :: iter,niter,nproblem,ncol
    !egvsol: block
        real(r8) :: bandstructure_energy
        integer :: i
    !invft: block
        integer, parameter :: safenumentries=10000000 ! or some other arbitrary number
        complex(r8), dimension(:,:), allocatable :: buf_rot
        complex(r8), dimension(:,:), allocatable :: buf_dm
        !complex(r8) :: expikr
        real(r8), dimension(:,:,:,:), allocatable :: buf_comm
        real(r8), dimension(:,:), allocatable :: buf_lv
        !real(r8), dimension(3) :: Rvector,kvector
        !real(r8), dimension(3) :: Rvector,kvector
        !real(r8) :: kdotr
        integer, dimension(:,:), allocatable :: buf_ia
        !integer :: ikp,ipair,iop,ikpoint
        integer :: ikp,ikpoint
        !integer :: a1,a2,s1,s2,b1,b2,nb1,nb2,ob1,ob2,l1,l2,g1,g2
        integer :: a1,a2,s1,s2,nb1,nb2,ob1,ob2,g1,g2
        integer :: n_pair_per_iter
        !integer :: iter,niter,plo,phi,poff,lnp,pair_ctr,ispin
        integer :: plo,phi,poff,lnp,pair_ctr
        integer :: j,k

    ! Start timers
    timer=mpi_wtime()
    t0=timer
    t1=t0
    mellantid=0.0_r8

    ! First we have to Fourier-transform to reciprocal space. Depending on the distribution
    ! of the realspace Hamiltonian this is various levels of annoying, and perhaps
    ! a little confusing. I have tried to explain as much as possible, but I will add
    ! more as more questions come.
    !ftr: block

        ! Reset the buffer to nothing. This is where we are going to store
        ! the Fourier-transformed Hamiltonian and Overlap. At the end of
        ! this block that is. Please not that 'KS%buf' holds matrices that
        ! are ScaLAPACK distributed, so the straight arrays here are only
        ! subsets of the whole thing.
        KS%buf%eigenvalue=0.0_r8
        select type(b=>KS%buf)
        type is(rl_kspace_eigenproblem_prob_real)
            b%hamiltonian=0.0_r8
            b%overlap=0.0_r8
            b%eigenvector=0.0_r8
            b%densitymatrix=0.0_r8
        type is(rl_kspace_eigenproblem_prob_complex)
            b%hamiltonian=0.0_r8
            b%overlap=0.0_r8
            b%eigenvector=0.0_r8
            b%densitymatrix=0.0_r8
        end select

        ! Decide on the number of times we have to iterate over all pairs and stuff
        ! to sort out all the communication. The logic here is that I calculate
        ! some number of elements of the resulting k-space quantities, then
        ! communicate the whole chunk. This is then repeated enough times so that
        ! everyone has everything they need. This approach is reasonably
        ! conservative with memory: I pick some arbitrary upper limit of to the
        ! amount of data to communicate at once. If there is more memory available,
        ! there will be fewer calls to communication and much faster. For small
        ! calculations one can usually afford to do everything in a single shot,
        ! that adds only minimal overhead, but it's not particularly scalable.
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! Here we don't need any communication. Allocate dummy arrays only.
            call mem%allocate(commbuf_ovl,[1,1,1],persistent=.false.,scalable=.false.)
            call mem%allocate(commbuf_ham,[1,1,1],persistent=.false.,scalable=.false.)
            commbuf_ovl=-rl_huge
            commbuf_ham=-rl_huge
            ! Only one iteration, and the column index things mean nothing, and
            ! the number of problems is irrelevant
            niter=1
            nproblem=-rl_hugeint
            if ( verbosity .gt. 0 ) then
                write(rl_iou,*) '... no communication needed for Fourier transform'
            endif
        type is(rl_realspace_matrix_mediumdistributed)
            ! So, in this case every group of ranks that work with this problem
            ! has a copy of the Hamiltonian. So the number of problems that we
            ! have to communicate over is just a single one!
            nproblem=1
            ! The smallest unit I will bother with is a number of whole columns
            ! at the time, otherwise it's just too tedious. How many columns could
            ! I sensibly fit into memory at the same time?
            f0=(KS%n_basis*32)        ! bytes per column, kinda
            f0=safemem/f0             ! How many columns can I store in memory?
            ncol=floor(f0)            ! Now as an integer.
            ncol=max(1,ncol)          ! At least one, for sure. Too annoying otherwise.
            ncol=min(ncol,KS%n_basis) ! And no more than what is actually necessary.
            ! Decide on how many times I have to iterate to have everything communicated.
            niter=0
            do
                if ( niter*ncol .ge. KS%n_basis ) exit
                niter=niter+1
            enddo
            ! Now make some space for the communication buffers. I think I will keep
            ! them complex, at least for a little while. It should be a good idea and
            ! cover the case where every k-point decides for itself wether it should
            ! be real or complex, without way too many if-clauses.
            call mem%allocate(commbuf_ovl,[KS%n_basis,ncol,1],persistent=.false.,scalable=.false.,supress_error=.true.)
            call mem%allocate(commbuf_ham,[KS%n_basis,ncol,1],persistent=.false.,scalable=.false.,supress_error=.true.)
            if ( verbosity .gt. 0 ) then
                write(rl_iou,*) '... semilocal communication needed for Fourier transform'
            endif
        type is(rl_realspace_matrix_fulldistributed)
            ! Way more annoying.
            if ( mw%talk ) write(*,*) 'FIXME fulldistribution FFT'
            call mw%destroy()
            stop
            if ( verbosity .gt. 0 ) then
                write(rl_iou,*) '... global communication needed for Fourier transforms'
            endif
        end select

        ! Make space for rotated buffers
        allocate(hbuf(p%n_species,p%n_species))
        allocate(obuf(p%n_species,p%n_species))
        do is1=1,p%n_species
        do is2=1,p%n_species
            call mem%allocate(hbuf(is1,is2)%m,[basis%species(is1)%n_basis,basis%species(is2)%n_basis],persistent=.false.,scalable=.false.)
            call mem%allocate(obuf(is1,is2)%m,[basis%species(is1)%n_basis,basis%species(is2)%n_basis],persistent=.false.,scalable=.false.)
            hbuf(is1,is2)%m=0.0_r8
            obuf(is1,is2)%m=0.0_r8
        enddo
        enddo

        ! Fetch the actual k-vector? Note factor 2*pi here.
        kvector=kmesh%ip( KS%index_irreducible_kpoint )%r*rl_twopi
        ! Which spin-channel are we working on now?
        ispin=KS%index_spin

        ! Make note of timing for initialization phase
        t1=mpi_wtime(); mellantid(1)=mellantid(1)+t1-t0; t0=t1

        ! Now go over all pairs. We will do it in a slightly convoluted way
        ! to throw people off. I will first loop over irreducible pairs, and
        ! for each pair I will unfold the irreducible pairs to the set of all
        ! pairs. In the case of no symmetry that is quite straightforward,
        ! but in the symmetric case we have to rotate a little back-and-forth.
        iterloop: do iter=1,niter
            ! Here decide on which columns to treat.
            select type(m=>rmtx)
            type is(rl_realspace_matrix_notdistributed)
                ! Do nothing!
            type is(rl_realspace_matrix_mediumdistributed)
                ! Decide on which columns to treat this iteration.
                clo=(iter-1)*ncol+1
                chi=min(iter*ncol,KS%n_basis)
                coff=(iter-1)*ncol
                ! Reset communication buffer
                commbuf_ovl=0.0_r8
                commbuf_ham=0.0_r8
            type is(rl_realspace_matrix_fulldistributed)
                ! Do something!
                if ( mw%talk ) write(*,*) 'FIXME fulldistribution FFT'
                call mw%destroy()
                stop
            end select
            !@TODO Insert problem loop here? That should cover everything.
            irrpairloop: do ipair=1,rmtx%n_irr_pair
                ! The two atoms involved in the irreducible pair
                ia1=rmtx%irr_pair(ipair)%a1
                ia2=rmtx%irr_pair(ipair)%a2
                ! Species of the irreducible pair
                is1=p%species(ia1)
                is2=p%species(ia2)
                ! Now go over all the pairs that this can unfold to:
                symunfoldloop: do ii=1,rmtx%irr_pair(ipair)%n_unfold_pair
                    jpair=rmtx%irr_pair(ipair)%unfold_index(ii)    ! Index of unfolded pair
                    iop=rmtx%irr_pair(ipair)%unfold_operation(ii)  ! Operation that unfolds. Negative index means operation+transpose
                    ja1=rmtx%full_pair(jpair)%a1           ! first atom of transformed pair
                    ja2=rmtx%full_pair(jpair)%a2           ! second atom of transformed pair

                    ! Quick exit test to see if this pair is currently relevant!
                    select type(m=>rmtx)
                    type is(rl_realspace_matrix_notdistributed)
                        ! Check if jpair pair has relevant entries for the local matrix
                        if ( m%full_pair(jpair)%n_ft_element .eq. 0 ) cycle symunfoldloop
                    type is(rl_realspace_matrix_mediumdistributed)
                        ! Equivalent check. Hmm. Enough to check column index, I think.
                        l=0
                        do b2=1,rmtx%full_pair(jpair)%nb2
                            l2=b2+basis%offset(ja2)
                            if ( l2 .ge. clo .and. l2 .le. chi ) l=l+1
                        enddo
                        if ( l .eq. 0 ) cycle symunfoldloop
                    type is(rl_realspace_matrix_fulldistributed)
                        ! Way more annoying.
                        if ( mw%talk ) write(*,*) 'FIXME fulldistribution FFT'
                        call mw%destroy()
                        stop
                    end select

                    ! If we made it here this pair is relevant. Fetch some intermediate information.
                    js1=p%species(ja1)                     ! first species of transformed pair
                    js2=p%species(ja2)                     ! second species of transformed pair
                    Rvector=rmtx%full_pair(jpair)%lv       ! Lattice vector for Fourier transform
                    kdotr=dot_product(rvector,kvector)     ! well, k dot R
                    expikr=cmplx(cos(kdotr),sin(kdotr),r8) ! phase factor. Plus or minus, who knows.

                    ! Now, here is for the slightly tricky part. I have to rotate the block
                    ! from the irreducible representation to here. Negative index means transpose
                    ! I rotate it into hbuf and obuf. That's why I allocate n_species x n_species work
                    ! buffers in the beginning, so that I always have a matrix of the appropriate size to store
                    ! intermediate results in, without have to use lots of allocations or half-empty arrays.
                    ! In the case of trivial symmetry, I do it slightly faster.
                    select case(sym%n_operation)
                    case(1)
                        ! This means there are no symmetry operations except
                        ! for the trivial. So either it's a straight copy or
                        ! a transpose, with not need for rotations.
                        select case(iop)
                        case(1)
                            ! Just a copy
                            obuf(js1,js2)%m=rmtx%irr_pair(ipair)%overlap
                            hbuf(js1,js2)%m=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin)
                        case(-1)
                            ! Transpose
                            obuf(js1,js2)%m=transpose(rmtx%irr_pair(ipair)%overlap)
                            hbuf(js1,js2)%m=transpose(rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin))
                        end select
                    case default
                        ! This is the case with one or more non-trivial symmetry operations.
                        ! Then we have to rotate the block associated with the irreducible
                        ! pair here.
                        if ( iop .lt. 0 ) then
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(-iop),&
                                forward=.true.,transposition=.true.,&
                                original_block=rmtx%irr_pair(ipair)%overlap,&
                                rotated_block=obuf(js1,js2)%m )
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(-iop),&
                                forward=.true.,transposition=.true.,&
                                original_block=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin),&
                                rotated_block=hbuf(js1,js2)%m )
                        else
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(iop),&
                                forward=.true.,transposition=.false.,&
                                original_block=rmtx%irr_pair(ipair)%overlap,&
                                rotated_block=obuf(js1,js2)%m )
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(iop),&
                                forward=.true.,transposition=.false.,&
                                original_block=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin),&
                                rotated_block=hbuf(js1,js2)%m )
                        endif
                    end select

                    ! At this point, the realspace Hamiltonian/Overlap is held in hbuf/obuf respectively.
                    ! What is left now is to take the realspace values, multiply by expikr and store in
                    ! the correct place. What the correct place is, that varies.
                    select type(m=>rmtx)
                    type is(rl_realspace_matrix_notdistributed)
                        ! This is by far the easiest: No communication at the end, just straight storage.
                        select type(b=>KS%buf)
                        type is(rl_kspace_eigenproblem_prob_real)
                            do jj=1,rmtx%full_pair(jpair)%n_ft_element
                                b1=rmtx%full_pair(jpair)%ft_element(1,jj)
                                b2=rmtx%full_pair(jpair)%ft_element(2,jj)
                                l1=rmtx%full_pair(jpair)%ft_element(3,jj)
                                l2=rmtx%full_pair(jpair)%ft_element(4,jj)
                                b%overlap(l1,l2)=b%overlap(l1,l2)+obuf(js1,js2)%m(b1,b2)*real(expikr)
                                b%hamiltonian(l1,l2)=b%hamiltonian(l1,l2)+hbuf(js1,js2)%m(b1,b2)*real(expikr)
                            enddo
                        type is(rl_kspace_eigenproblem_prob_complex)
                            do jj=1,rmtx%full_pair(jpair)%n_ft_element
                                b1=rmtx%full_pair(jpair)%ft_element(1,jj)
                                b2=rmtx%full_pair(jpair)%ft_element(2,jj)
                                l1=rmtx%full_pair(jpair)%ft_element(3,jj)
                                l2=rmtx%full_pair(jpair)%ft_element(4,jj)
                                b%overlap(l1,l2)=b%overlap(l1,l2)+obuf(js1,js2)%m(b1,b2)*expikr
                                b%hamiltonian(l1,l2)=b%hamiltonian(l1,l2)+hbuf(js1,js2)%m(b1,b2)*expikr
                            enddo
                        end select
                    type is(rl_realspace_matrix_mediumdistributed)
                        ! Store the relevant rows for this iteration.
                        do b2=1,rmtx%full_pair(jpair)%nb2
                            l2=b2+basis%offset(ja2)
                            if ( l2 .lt. clo ) cycle
                            if ( l2 .gt. chi ) cycle
                            do b1=1,rmtx%full_pair(jpair)%nb1
                                l1=b1+basis%offset(ja1)
                                commbuf_ovl(l1,l2-coff,1)=commbuf_ovl(l1,l2-coff,1)+obuf(js1,js2)%m(b1,b2)*expikr
                                commbuf_ham(l1,l2-coff,1)=commbuf_ham(l1,l2-coff,1)+hbuf(js1,js2)%m(b1,b2)*expikr
                            enddo
                        enddo
                    type is(rl_realspace_matrix_fulldistributed)
                        ! Way more annoying.
                        if ( mw%talk ) write(*,*) 'FIXME fulldistribution FFT'
                        call mw%destroy()
                        stop
                    end select
                enddo symunfoldloop
            enddo irrpairloop

            ! Now sort out communication to make sure the right rank has the right thing.
            select type(m=>rmtx)
            type is(rl_realspace_matrix_notdistributed)
                ! Do nothing, we already have everything!
            type is(rl_realspace_matrix_mediumdistributed)
                ! Now we have to accumulate over local ranks
                call KS%ml%allreduce('sum',commbuf_ovl)
                call KS%ml%allreduce('sum',commbuf_ham)
                ! And then eventually store things in the right place.
                do b2=clo,chi
                do b1=1,KS%n_basis
                    if ( KS%bl%rank_from_global_indices(KS%n_basis,KS%n_basis,KS%buf%blocksize,b1,b2) .ne. KS%ml%r ) cycle
                    l1=KS%bl%local_row_from_global(KS%buf%blocksize,b1)
                    l2=KS%bl%local_col_from_global(KS%buf%blocksize,b2)
                    select type(b=>KS%buf)
                    type is(rl_kspace_eigenproblem_prob_real)
                        b%overlap(l1,l2)=real(commbuf_ovl(b1,b2-coff,1),r8)
                        b%hamiltonian(l1,l2)=real(commbuf_ham(b1,b2-coff,1),r8)
                    type is(rl_kspace_eigenproblem_prob_complex)
                        b%overlap(l1,l2)=commbuf_ovl(b1,b2-coff,1)
                        b%hamiltonian(l1,l2)=commbuf_ham(b1,b2-coff,1)
                    end select
                enddo
                enddo
            type is(rl_realspace_matrix_fulldistributed)
                ! Now we have to sum over the global communicator
            end select
        enddo iterloop

        ! And a little cleanup at the end, now that we are done with the transform!
        do is1=1,p%n_species
        do is2=1,p%n_species
            call mem%deallocate(obuf(is1,is2)%m,persistent=.false.,scalable=.false.)
            call mem%deallocate(hbuf(is1,is2)%m,persistent=.false.,scalable=.false.)
        enddo
        enddo
        deallocate(obuf)
        deallocate(hbuf)
        call mem%deallocate(commbuf_ovl,persistent=.false.,scalable=.false.)
        call mem%deallocate(commbuf_ham,persistent=.false.,scalable=.false.)

        ! Timings for Fourier transform
        t1=mpi_wtime(); mellantid(2)=mellantid(2)+t1-t0; t0=t1
    !end block ftr

    if ( verbosity .gt. 0 ) then
        write(rl_iou,*) '... Fourier-transformed Hamiltonian and Overlap'
    endif

    ! Now, if we made it here KS%buf%overlap and KS%buf%hamiltonian should hold the
    ! exactly that, so solve the KS equations! ELSI makes this rather simple.
    !egvsol: block

        ! ELSI is pretty neat. This will fix everything in a single call, almost!
        ! I don't know if I actually need the eigenvalues, eigenvectors or anything else really?
        ! time will tell what I actually keep in the end.
        select type(b=>KS%buf)
        type is(rl_kspace_eigenproblem_prob_real)
            call elsi_dm_real(KS%eh,b%hamiltonian,b%overlap,b%densitymatrix,bandstructure_energy)
            call elsi_get_evec_real(KS%eh,b%eigenvector)
        type is(rl_kspace_eigenproblem_prob_complex)
            call elsi_dm_complex(KS%eh,b%hamiltonian,b%overlap,b%densitymatrix,bandstructure_energy)
            call elsi_get_evec_complex(KS%eh,b%eigenvector)
        end select

        ! Time for eigenvalue solution
        t1=mpi_wtime(); mellantid(3)=mellantid(3)+t1-t0; t0=t1

        ! Eigenvalues might come handy
        call elsi_get_eval(KS%eh,KS%buf%eigenvalue)
        ! Fetch chemical potential and entropy?
        call elsi_get_mu(KS%eh,KS%chemical_potential)
        call elsi_get_entropy(KS%eh,KS%entropy)

        ! Time for occupation numbers
        t1=mpi_wtime(); mellantid(4)=mellantid(4)+t1-t0; t0=t1
    !end block egvsol

    if ( verbosity .gt. 0 ) then
        write(rl_iou,*) '... solved eigenvalue problem'
        write(rl_iou,*) '... building density matrix'
    endif

    ! Inverse-transform the density matrix to realspace. Basically the same as the Fourier transform
    ! above, but backwards, with yet another annoying communication pattern.
    !invft: block

        ! First things first, make sure we zero out the realspace
        ! density matrix before the transform. Always a good idea.
        do ipair=1,rmtx%n_irr_pair
            rmtx%irr_pair(ipair)%densitymatrix=0.0_r8
        enddo

        ! Largest number of basis functions per species, need to know this
        ! to decide on communucation scheme.
        l=0
        do s1=1,p%n_species
            l=max(l,basis%species(s1)%n_basis)
        enddo

        ! Then decide how communication is going to happen, since everything
        ! is distributed all over the place.
        select type(rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! We can do all pairs in one go! I assume. Probably. Maybe not.
            n_pair_per_iter=safenumentries/real(KS%n_spin*l**2,r8)
            n_pair_per_iter=min(n_pair_per_iter,rmtx%n_irr_pair)
            niter=0
            do
                if ( niter*n_pair_per_iter .ge. rmtx%n_irr_pair ) then
                    exit
                else
                    niter=niter+1
                endif
            enddo
        type is(rl_realspace_matrix_mediumdistributed)
            ! We will likely have to do a few pairs at the time.
            n_pair_per_iter=safenumentries/real(KS%n_spin*l**2,r8)
            n_pair_per_iter=min(n_pair_per_iter,rmtx%n_irr_pair_global)
            niter=0
            do
                if ( niter*n_pair_per_iter .ge. rmtx%n_irr_pair_global ) then
                    exit
                else
                    niter=niter+1
                endif
            enddo
        type is(rl_realspace_matrix_fulldistributed)
            ! We will likely have to do a few pairs at the time.
            n_pair_per_iter=safenumentries/real(KS%n_spin*l**2,r8)
            n_pair_per_iter=min(n_pair_per_iter,rmtx%n_irr_pair_global)
            niter=0
            do
                if ( niter*n_pair_per_iter .ge. rmtx%n_irr_pair_global ) then
                    exit
                else
                    niter=niter+1
                endif
            enddo
        end select



        ! Make space for communication buffers
        call mem%allocate(buf_comm,[l,l,KS%n_spin,n_pair_per_iter],persistent=.false.,scalable=.false.)
        call mem%allocate(buf_lv,[3,n_pair_per_iter],persistent=.false.,scalable=.false.)
        call mem%allocate(buf_ia,[5,n_pair_per_iter],persistent=.false.,scalable=.false.)
        buf_comm=0.0_r8
        buf_lv=0.0_r8
        buf_ia=0

        ! Make some space for the rotation matrix and the transformed density matrix.
        call mem%allocate(buf_rot,[KS%buf%n_row_local,KS%buf%n_col_local],persistent=.true.,scalable=.true.)
        call mem%allocate(buf_dm,[KS%buf%n_row_local,KS%buf%n_col_local],persistent=.true.,scalable=.true.)
        buf_rot=0.0_r8
        buf_dm=0.0_r8

        iterloop2: do iter=1,niter
            ! The iterations are iterations over the pairs, in case I can't
            ! do all the pairs at once. First reset the communication buffer
            buf_comm=0.0_r8
            buf_lv=0.0_r8
            buf_ia=0
            ! Fetch information about the currently relevant pairs
            select type(m=>rmtx)
            type is(rl_realspace_matrix_notdistributed)
                ! No comm needed, we already have all the pairs
                plo=(iter-1)*n_pair_per_iter+1
                phi=min(iter*n_pair_per_iter,m%n_irr_pair)
                poff=(iter-1)*n_pair_per_iter
                pair_ctr=0
                do ipair=plo,phi
                    pair_ctr=pair_ctr+1
                    a1=m%irr_pair(ipair)%a1
                    a2=m%irr_pair(ipair)%a2
                    buf_lv(:,pair_ctr)=m%irr_pair(ipair)%lv
                    buf_ia(:,pair_ctr)=[m%irr_pair(ipair)%nb1,m%irr_pair(ipair)%nb2,basis%offset(a1),basis%offset(a2),ipair]
                enddo
            type is(rl_realspace_matrix_mediumdistributed)
                ! We have some pairs per rank, collect them
                plo=(iter-1)*n_pair_per_iter+1
                phi=min(iter*n_pair_per_iter,m%n_irr_pair_global)
                poff=(iter-1)*n_pair_per_iter
                pair_ctr=0
                do ipair=1,m%n_irr_pair
                    j=ipair+m%irr_offset
                    if ( j .lt. plo ) cycle
                    if ( j .gt. phi ) cycle
                    j=j-poff
                    a1=m%irr_pair(ipair)%a1
                    a2=m%irr_pair(ipair)%a2
                    buf_lv(:,j)=m%irr_pair(ipair)%lv
                    buf_ia(:,j)=[m%irr_pair(ipair)%nb1,m%irr_pair(ipair)%nb2,basis%offset(a1),basis%offset(a2),j+poff]
                    pair_ctr=pair_ctr+1
                enddo
                call KS%ml%allreduce('sum',pair_ctr)
                call KS%ml%allreduce('sum',buf_lv)
                call KS%ml%allreduce('sum',buf_ia)
            type is(rl_realspace_matrix_fulldistributed)
                ! We will likely have to do a few pairs at the time.
                ! to come to a sensible conclusion I need to know how
                ! large a pair buffer maximally is:
                if ( mw%talk ) write(*,*) 'continue here multiproc'
                call mw%destroy()
                stop
            end select

            kploop: do ikp=1,kmesh%ip( KS%index_irreducible_kpoint )%n_full_point
                ! Full index to k-point
                ikpoint=kmesh%ip( KS%index_irreducible_kpoint )%index_full_point(ikp)
                ! Which operation to transform the density matrix with
                iop=kmesh%ip( KS%index_irreducible_kpoint )%operation_full_point(ikp)
                ! Fetch the actual k-vector? Note factor 2*pi here.
                kvector=kmesh%ap( ikpoint )%r*rl_twopi
                ! Which spin-channel are we working on now?
                ispin=KS%index_spin

                ! Now either transform the density matrix, to this k-point, or not.
                select case(sym%n_operation)
                case(1)
                    ! This is the trivial case, don't need to do anything
                    select type(b=>KS%buf)
                    type is(rl_kspace_eigenproblem_prob_real)
                        buf_dm=b%densitymatrix
                    type is(rl_kspace_eigenproblem_prob_complex)
                        buf_dm=b%densitymatrix
                    end select
                case default
                    ! we have to rotate and stuff. Frustrating. How about no for now.
                    if ( mw%talk ) write(*,*) 'continue here rotate kspace multiproc'
                    call mw%destroy()
                    stop
                end select

                pairloop: do i=1,pair_ctr
                    Rvector=buf_lv(:,i)                    ! Lattice vector for Fourier transform
                    kdotr=dot_product(rvector,kvector)     ! well, k dot R
                    expikr=cmplx(cos(kdotr),sin(kdotr),r8) ! phase factor. Plus or minus, who knows.
                    expikr=expikr*kmesh%ap(ikpoint)%integration_weight ! to normalize the transform

                    nb1=buf_ia(1,i)
                    nb2=buf_ia(2,i)
                    ob1=buf_ia(3,i)
                    ob2=buf_ia(4,i)

                    ! Start storing things
                    do b1=1,nb1 ! index 1 in pair block
                    do b2=1,nb2 ! index 2 in pair block
                        g1=b1+ob1 ! global index 1 in density matrix
                        g2=b2+ob2 ! global index 2 in density matrix
                        ! Skip if we don't have anything on this rank.
                        j=KS%bl%rank_from_global_indices(KS%n_basis,KS%n_basis,KS%buf%blocksize,g1,g2)
                        if ( j .ne. KS%ml%r ) cycle
                        ! Get the local indices
                        l1=KS%bl%local_row_from_global(KS%buf%blocksize,g1)
                        l2=KS%bl%local_col_from_global(KS%buf%blocksize,g2)
                        ! Accumulate transform
                        buf_comm(b1,b2,ispin,i)=buf_comm(b1,b2,ispin,i)+real(buf_dm(l1,l2)*expikr,r8)
                    enddo
                    enddo
                enddo pairloop
            enddo kploop

            ! Add together the transform over (global) ranks:
            call mw%allreduce('sum',buf_comm)
            ! Store in the right place:
            select type(m=>rmtx)
            type is(rl_realspace_matrix_notdistributed)
                ! Store everywhere, no problem.
                pair_ctr=0
                do ipair=plo,phi
                    pair_ctr=pair_ctr+1
                    nb1=m%irr_pair(ipair)%nb1
                    nb2=m%irr_pair(ipair)%nb2
                    m%irr_pair(ipair)%densitymatrix=buf_comm(1:nb1,1:nb2,:,pair_ctr)
                enddo
            type is(rl_realspace_matrix_mediumdistributed)
                ! Store maybe not averywhere.
                do i=1,pair_ctr
                    ipair=buf_ia(5,i)        ! global index to pair
                    ipair=ipair-m%irr_offset ! local index of pair
                    if ( ipair .lt. 1 ) cycle
                    if ( ipair .gt. m%n_irr_pair ) cycle
                    nb1=m%irr_pair(ipair)%nb1
                    nb2=m%irr_pair(ipair)%nb2
                    m%irr_pair(ipair)%densitymatrix=buf_comm(1:nb1,1:nb2,:,i)
                enddo
            type is(rl_realspace_matrix_fulldistributed)
                ! We will likely have to do a few pairs at the time.
                ! to come to a sensible conclusion I need to know how
                ! large a pair buffer maximally is:
                if ( mw%talk ) write(*,*) 'continue here multiproc'
                call mw%destroy()
                stop
            end select
        enddo iterloop2

        ! Note time to get density matrix
        t1=mpi_wtime(); mellantid(5)=mellantid(5)+t1-t0; t0=t1

        ! Some cleanup at the end
        call mem%deallocate(buf_comm,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_lv,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_ia,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_rot,persistent=.true.,scalable=.true.)
        call mem%deallocate(buf_dm,persistent=.true.,scalable=.true.)
    !end block invft

    if ( verbosity .gt. 0 ) then
        write(rl_iou,*) '... calculated realspace density matrix'
    endif

    ! Check that I'm not stupid with memory.
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

    ! Make note of the total time:
    mellantid(6)=mpi_wtime()-timer
    ! Store the timers in the persistent thing for analysis later.
    tmr%init          =tmr%init          +mellantid(1)
    tmr%forward_ft    =tmr%forward_ft    +mellantid(2)
    tmr%solve         =tmr%solve         +mellantid(3)
    tmr%occnumbers    =tmr%occnumbers    +mellantid(4)
    tmr%densitymatrix =tmr%densitymatrix +mellantid(5)
    tmr%idle          =tmr%idle          +mellantid(6)-sum(mellantid(1:5))
end subroutine

subroutine solve_KS_singleproc(KS,rmtx,basis,p,sym,kmesh,mw,mem,verbosity,tmr)
    !> KS solver handle. The multiproc version.
    type(rl_kspace_eigenproblem_singleproc), intent(inout) :: KS
    !> distributed matrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> k-point mesh
    type(rl_kpoint_mesh), intent(in) :: kmesh
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timing information
    type(rl_timer_eigenvalue), intent(inout) :: tmr

    ! helpers for symmetry transform
    type(rl_mom_real), dimension(:,:,:), allocatable :: hbuf
    type(rl_mom_real), dimension(:,:), allocatable :: obuf
    type(rl_mom_complex), dimension(:), allocatable :: dum_ovl,dum_ham
    ! buffers for fourier transform
    complex(r8), dimension(:,:,:), allocatable :: buf_ham_c
    complex(r8), dimension(:,:), allocatable :: buf_ovl_c
    real(r8), dimension(:,:,:), allocatable :: buf_ham_r
    real(r8), dimension(:,:), allocatable :: buf_ovl_r
    ! timers
    real(r8) :: timer,t0,t1
    real(r8), dimension(8) :: mellantid
    ! counters
    integer :: ikp
    !init: block
        integer :: is1,is2,ispin
        integer :: i
        !fouriertransform: block
            complex(r8) :: expikr
            real(r8), dimension(3) :: kvector,rvector
            real(r8) :: kdotr
            !integer :: ispin,ipair,jpair,iop,ii,jj
            integer :: ipair,jpair,iop,ii,jj
            integer :: i1,i2,j1,j2
            !integer :: ia1,ia2,ja1,ja2,is1,is2,js1,js2
            integer :: ia1,ia2,ja1,ja2,js1,js2
            integer :: ctr
        !egvsolve: block
            real(r8) :: f0,f1
            !integer :: ispin,i,j
            integer :: j
    !cleanup: block
        !integer :: is1,is2,ispin
    !muandocc: block
        real(r8), dimension(:,:,:), allocatable :: buf_egv,buf_occ
        real(r8), dimension(:), allocatable :: buf_kwt
        !real(r8) :: f0
        !integer :: ispin
        !integer :: i,j,k
        integer :: k
    !densitymatrix: block
        type(rl_mom_real), dimension(:,:), allocatable :: BM
        complex(r8), dimension(:,:), allocatable :: wegv,wegw,wdm,wrot !,wdd
        !complex(r8) :: expikr,tracectr,c0
        complex(r8) :: tracectr,c0
        real(r8), dimension(:), allocatable :: wocc
        real(r8), dimension(3) :: qv,qw,v0,v1
        !real(r8), dimension(3) :: kvector,rvector
        !real(r8) :: kdotr,f0,elcount0,elcount1
        real(r8) :: elcount0,elcount1
        !integer :: a1,a2,nb1,nb2,s1,s2,i1,i2,j1,j2
        integer :: a1,a2,nb1,nb2,s1,s2
        !integer :: i,j,ii,jj
        !integer :: ispin,ikp,ikpoint,jkpoint,istatemax,iop
        integer :: ikpoint,jkpoint,istatemax
        !integer :: ipair

    ! Start timers
    timer=mpi_wtime()
    t0=timer
    t1=t0
    mellantid=0.0_r8

    !init: block

        ! Make sure the Hamiltonian is distributed the way I think it is.
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! Do nothing, all is fine.
        type is(rl_realspace_matrix_mediumdistributed)
            call rl_stop_gracefully(['Distributed realspace Hamiltonian in single-proc mode, inefficient and strange.'],rl_exitcode_param,mw%comm)
        type is(rl_realspace_matrix_fulldistributed)
            call rl_stop_gracefully(['Distributed realspace Hamiltonian in single-proc mode, inefficient and strange.'],rl_exitcode_param,mw%comm)
        end select

        ! Make space for rotated buffers
        allocate(hbuf(p%n_species,p%n_species,KS%n_spin))
        allocate(obuf(p%n_species,p%n_species))
        do is1=1,p%n_species
        do is2=1,p%n_species
            do ispin=1,KS%n_spin
                call mem%allocate(hbuf(is1,is2,ispin)%m,[basis%species(is1)%n_basis,basis%species(is2)%n_basis],persistent=.false.,scalable=.false.)
                hbuf(is1,is2,ispin)%m=0.0_r8
            enddo
            call mem%allocate(obuf(is1,is2)%m,[basis%species(is1)%n_basis,basis%species(is2)%n_basis],persistent=.false.,scalable=.false.)
            obuf(is1,is2)%m=0.0_r8
        enddo
        enddo

        ! Make space for the Hamiltonian and overlap buffers.
        if ( KS%n_kpoint .gt. 0 ) then
            select type(k=>KS%kpoint(1))
            type is(rl_kspace_eigenproblem_kpoint_real)
                call mem%allocate(buf_ham_r,[KS%n_basis,KS%n_basis,KS%n_spin],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ovl_r,[KS%n_basis,KS%n_basis],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ham_c,[1,1,1],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ovl_c,[1,1],persistent=.false.,scalable=.false.)
                buf_ham_r=0.0_r8
                buf_ovl_r=0.0_r8
                buf_ham_c=0.0_r8
                buf_ovl_c=0.0_r8
            type is(rl_kspace_eigenproblem_kpoint_complex)
                call mem%allocate(buf_ham_r,[1,1,1],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ovl_r,[1,1],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ham_c,[KS%n_basis,KS%n_basis,KS%n_spin],persistent=.false.,scalable=.false.)
                call mem%allocate(buf_ovl_c,[KS%n_basis,KS%n_basis],persistent=.false.,scalable=.false.)
                buf_ham_r=0.0_r8
                buf_ovl_r=0.0_r8
                buf_ham_c=0.0_r8
                buf_ovl_c=0.0_r8
            end select
        else
            ! no k-points on this rank, allocate dummies
            call mem%allocate(buf_ham_r,[1,1,1],persistent=.false.,scalable=.false.)
            call mem%allocate(buf_ovl_r,[1,1],persistent=.false.,scalable=.false.)
            call mem%allocate(buf_ham_c,[1,1,1],persistent=.false.,scalable=.false.)
            call mem%allocate(buf_ovl_c,[1,1],persistent=.false.,scalable=.false.)
            buf_ham_r=0.0_r8
            buf_ovl_r=0.0_r8
            buf_ham_c=0.0_r8
            buf_ovl_c=0.0_r8
        endif

        ! Make note of timing for initialization phase
        t1=mpi_wtime(); mellantid(1)=mellantid(1)+t1-t0; t0=t1

        !@TODO If the Fourier transform becomes a bottleneck, I can create temporary space and
        ! unfold the realspace matrix only once. That wastes a bit of memory, but it is decidedly
        ! faster. Maybe that is worth thinking about.
allocate(dum_ovl(kmesh%n_irr_kpoint))
allocate(dum_ham(kmesh%n_irr_kpoint))
do i=1,kmesh%n_irr_kpoint
    allocate(dum_ovl(i)%m(KS%n_basis,KS%n_basis))
    allocate(dum_ham(i)%m(KS%n_basis,KS%n_basis))
    dum_ovl(i)%m=0.0_r8
    dum_ham(i)%m=0.0_r8
enddo
    !end block init

    ! Loop over k-points
    kploop: do ikp=1,KS%n_kpoint

        ! So, we have to Fourier transform the Hamiltonian and overlap at every k-point.
        !fouriertransform: block

            ! Fetch the actual k-vector? Note factor 2*pi here.
            kvector=kmesh%ip( KS%kpoint(ikp)%irreducible_index )%r*rl_twopi
            ! reset buffers
            buf_ham_r=0.0_r8
            buf_ovl_r=0.0_r8
            buf_ham_c=0.0_r8
            buf_ovl_c=0.0_r8
            ctr=0

            ! Go over all pairs for the transform!
            irrpairloop: do ipair=1,rmtx%n_irr_pair
                ! The two atoms involved in the irreducible pair
                ia1=rmtx%irr_pair(ipair)%a1
                ia2=rmtx%irr_pair(ipair)%a2
                ! Species of the irreducible pair
                is1=p%species(ia1)
                is2=p%species(ia2)
                ! Now go over all the pairs that this can unfold to:
                symunfoldloop: do ii=1,rmtx%irr_pair(ipair)%n_unfold_pair
                    jpair=rmtx%irr_pair(ipair)%unfold_index(ii)    ! Index of unfolded pair
                    iop=rmtx%irr_pair(ipair)%unfold_operation(ii)  ! Operation that unfolds. Negative index means operation+transpose
                    ja1=rmtx%full_pair(jpair)%a1           ! first atom of transformed pair
                    ja2=rmtx%full_pair(jpair)%a2           ! second atom of transformed pair
                    js1=p%species(ja1)                     ! first species of transformed pair
                    js2=p%species(ja2)                     ! second species of transformed pair
                    Rvector=rmtx%full_pair(jpair)%lv       ! Lattice vector for Fourier transform
                    kdotr=dot_product(Rvector,kvector)     ! well, k dot R
                    expikr=cmplx(cos(kdotr),sin(kdotr),r8) ! phase factor. Plus or minus, who knows.

                    ! Now, here is for the slightly tricky part. I have to rotate the block
                    ! from the irreducible representation to here. Negative index means transpose
                    ! I rotate it into hbuf and obuf. That's why I allocate n_species x n_species work
                    ! buffers in the beginning, so that I always have a matrix of the appropriate size to store
                    ! intermediate results in, without have to use lots of allocations or half-empty arrays.
                    ! Anyway, if there is no symmetry, there is no need for actual rotations and it can
                    ! be resolved very quickly.
                    select case(sym%n_operation)
                    case(1)
                        ! This means there are no symmetry operations except
                        ! for the trivial. So either it's a straight copy or
                        ! a transpose, with not need for multiplications.
                        select case(iop)
                        case(1)
                            ! Just a copy
                            obuf(js1,js2)%m=rmtx%irr_pair(ipair)%overlap
                            do ispin=1,KS%n_spin
                                hbuf(js1,js2,ispin)%m=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin)
                            enddo
                        case(-1)
                            ! Transpose
                            obuf(js1,js2)%m=transpose(rmtx%irr_pair(ipair)%overlap)
                            do ispin=1,KS%n_spin
                                hbuf(js1,js2,ispin)%m=transpose(rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin))
                            enddo
                        end select
                    case default
                        ! This is the case with one or more non-trivial symmetry operations.
                        ! Then we have to rotate the block associated with the irreducible
                        ! pair here.
                        if ( iop .lt. 0 ) then
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(-iop),&
                                forward=.true.,transposition=.true.,&
                                original_block=rmtx%irr_pair(ipair)%overlap,&
                                rotated_block=obuf(js1,js2)%m )
                            do ispin=1,KS%n_spin
                                call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(-iop),&
                                    forward=.true.,transposition=.true.,&
                                    original_block=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin),&
                                    rotated_block=hbuf(js1,js2,ispin)%m )
                            enddo
                        else
                            call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(iop),&
                                forward=.true.,transposition=.false.,&
                                original_block=rmtx%irr_pair(ipair)%overlap,&
                                rotated_block=obuf(js1,js2)%m )
                            do ispin=1,KS%n_spin
                                call rmtx%irr_pair(ipair)%rotate( basis,p,sym%op(iop),&
                                    forward=.true.,transposition=.false.,&
                                    original_block=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin),&
                                    rotated_block=hbuf(js1,js2,ispin)%m )
                            enddo
                        endif
                    end select
                    ! At this point, the realspace Hamiltonian/Overlap is held in hbuf/obuf respectively.
                    ! What is left now is to take the realspace values, multiply by expikr and store in
                    ! the correct place.
                    i1=basis%offset(ja1)+1
                    i2=basis%offset(ja1)+rmtx%full_pair(jpair)%nb1
                    j1=basis%offset(ja2)+1
                    j2=basis%offset(ja2)+rmtx%full_pair(jpair)%nb2
                    select type(k=>KS%kpoint(ikp))
                    type is(rl_kspace_eigenproblem_kpoint_real)
                        do ispin=1,KS%n_spin
                            buf_ham_r(i1:i2,j1:j2,ispin)=buf_ham_r(i1:i2,j1:j2,ispin)+hbuf(js1,js2,ispin)%m*real(expikr)
                        enddo
                        buf_ovl_r(i1:i2,j1:j2)=buf_ovl_r(i1:i2,j1:j2)+obuf(js1,js2)%m*real(expikr)
                        ctr=ctr+1
                    type is(rl_kspace_eigenproblem_kpoint_complex)
                        do ispin=1,KS%n_spin
                            buf_ham_c(i1:i2,j1:j2,ispin)=buf_ham_c(i1:i2,j1:j2,ispin)+hbuf(js1,js2,ispin)%m*expikr
                        enddo
                        buf_ovl_c(i1:i2,j1:j2)=buf_ovl_c(i1:i2,j1:j2)+obuf(js1,js2)%m*expikr
                        ctr=ctr+1
                    end select
                    ! And that should be it for this pair!
                enddo symunfoldloop
            enddo irrpairloop


            ! Sanity check. Should never trigger, ideally.
            if ( ctr .ne. rmtx%n_full_pair ) then
                call rl_stop_gracefully(['Indexing is strange in Fourier transform'],rl_exitcode_symmetry,mw%comm)
            endif
        !end block fouriertransform
        t1=mpi_wtime(); mellantid(2)=mellantid(2)+t1-t0; t0=t1

        ! Now solve the KS equations in reciprocal space
        !egvsolve: block

            ! First I have a slight sanity check here: I actually calculate the full Hamiltonian and
            ! overlap, so that I can check that they really are Hermitian. They always should be, and
            ! this is a good place to catch any bugs in the symmetry stuff, since any mistake in the
            ! transformation stuff will almost inevitably break the Hermiticity of the Fourier-transformed
            ! things.
            f0=0.0_r8
            f1=0.0_r8
            select type(k=>KS%kpoint(ikp))
            type is(rl_kspace_eigenproblem_kpoint_real)
                do i=1,KS%n_basis
                do j=i+1,KS%n_basis
                    f0=f0+abs(buf_ovl_r(i,j)-buf_ovl_r(j,i))
                enddo
                enddo
                do ispin=1,KS%n_spin
                    do i=1,KS%n_basis
                    do j=i+1,KS%n_basis
                        f1=f1+abs(buf_ham_r(i,j,ispin)-buf_ham_r(j,i,ispin))
                    enddo
                    enddo
                enddo
            type is(rl_kspace_eigenproblem_kpoint_complex)
                do i=1,KS%n_basis
                    f0=f0+abs(aimag(buf_ovl_c(i,i)))
                    do j=i+1,KS%n_basis
                        f0=f0+abs(buf_ovl_c(i,j)-conjg(buf_ovl_c(j,i)))
                    enddo
                enddo
                do ispin=1,KS%n_spin
                    do i=1,KS%n_basis
                        f1=f1+abs(aimag(buf_ham_c(i,i,ispin)))
                        do j=i+1,KS%n_basis
                            f1=f1+abs(buf_ham_c(i,j,ispin)-conjg(buf_ham_c(j,i,ispin)))
                        enddo
                    enddo
                enddo
            end select
            f0=f0/(KS%n_basis**2)
            f1=f1/(KS%n_basis**2)
            if ( f0+f1 .gt. 1E-14_r8 ) then
                write(*,*) 'nonherm',f0,f1
                !call rl_stop_gracefully(['Non-Hermitian overlap and Hamiltonian, should never happen'],rl_exitcode_symmetry,mw%comm)
            endif

            ! Now that I am sure that things are Hermitian, we can solve the eigenvalue problem
            select type(k=>KS%kpoint(ikp))
            type is(rl_kspace_eigenproblem_kpoint_real)
                do ispin=1,KS%n_spin
                    call elsi_ev_real(KS%eh,buf_ham_r(:,:,ispin),buf_ovl_r,k%eigenvalue(:,ispin),k%eigenvector(:,:,ispin))
                enddo
            type is(rl_kspace_eigenproblem_kpoint_complex)
dum_ovl( KS%kpoint(ikp)%irreducible_index )%m=buf_ovl_c
dum_ham( KS%kpoint(ikp)%irreducible_index )%m=buf_ham_c(:,:,ispin)
                do ispin=1,KS%n_spin
                    call elsi_ev_complex(KS%eh,buf_ham_c(:,:,ispin),buf_ovl_c,k%eigenvalue(:,ispin),k%eigenvector(:,:,ispin))
                enddo
            end select
        !end block egvsolve
        t1=mpi_wtime(); mellantid(3)=mellantid(3)+t1-t0; t0=t1
    enddo kploop

    ! And a little cleanup at the end, now that we are done with the transform!
    !cleanup: block

        do is1=1,p%n_species
        do is2=1,p%n_species
            do ispin=1,KS%n_spin
                call mem%deallocate(hbuf(is1,is2,ispin)%m,persistent=.false.,scalable=.false.)
            enddo
            call mem%deallocate(obuf(is1,is2)%m,persistent=.false.,scalable=.false.)
        enddo
        enddo
        deallocate(obuf)
        deallocate(hbuf)
        call mem%deallocate(buf_ham_r,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_ovl_r,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_ham_c,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_ovl_c,persistent=.false.,scalable=.false.)
    !end block cleanup

    ! Calculate the chemical potential and occupation numbers
    !muandocc: block

do i=1,kmesh%n_irr_kpoint
    call mw%allreduce('sum',dum_ovl(i)%m)
    call mw%allreduce('sum',dum_ham(i)%m)
enddo

        ! It makes sense to use ELSI to calculate the chemical potential and
        ! occupation numbers, even though it is not well documented. The code looks
        ! perfectly fine, so I'm pretty sure it's safe. First I have to collect
        ! KS eigenvalues from all the ranks, as well as create some dummy arrays for
        ! that interface.
        call mem%allocate(buf_egv,[KS%n_state,KS%n_spin,kmesh%n_irr_kpoint],persistent=.false.,scalable=.false.)
        call mem%allocate(buf_occ,[KS%n_state,KS%n_spin,kmesh%n_irr_kpoint],persistent=.false.,scalable=.false.)
        call mem%allocate(buf_kwt,kmesh%n_irr_kpoint,persistent=.false.,scalable=.false.)
        buf_egv=0.0_r8
        buf_occ=0.0_r8
        buf_kwt=0.0_r8
        ! Fetch eigenvalues from all ranks
        do i=1,KS%n_kpoint
            j=KS%kpoint(i)%irreducible_index
            do k=1,KS%n_state
                buf_egv(k,:,j)=KS%kpoint(i)%eigenvalue(k,:)
            enddo
        enddo
        call mw%allreduce('sum',buf_egv)
        ! Fetch k-point weights. Not sure how ELSI expects it.
        do i=1,kmesh%n_irr_kpoint
            buf_kwt(i)=kmesh%ip(i)%integration_weight !*kmesh%n_full_kpoint
        enddo
        ! Get Fermi level and occupation numbers
        call elsi_compute_mu_and_occ(ks%eh,KS%n_electron,KS%n_state,KS%n_spin,kmesh%n_irr_kpoint,buf_kwt,buf_egv,buf_occ,KS%chemical_potential)
        ! Get the electronic entropy
        call elsi_compute_entropy(ks%eh,KS%n_state,KS%n_spin,kmesh%n_irr_kpoint,buf_kwt,buf_egv,buf_occ,KS%chemical_potential,KS%entropy)

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... got chemical potential and entropy'
        endif

        ! Keep the occupation numbers!
        f0=0.0_r8
        do i=1,KS%n_kpoint
            j=KS%kpoint(i)%irreducible_index
            KS%kpoint(i)%occupation(1:KS%n_state,:)=buf_occ(:,:,j)
            ! Count electrons via occupation numbers?
            !f0=f0+sum(KS%kpoint(i)%occupation)*kmesh%ip(j)%integration_weight
        enddo
        !call mw%allreduce('sum',f0)
        !if ( mw%talk ) write(*,*) 'Count electrons via occupation numbers',f0,abs(f0-KS%n_electron)

        call mem%deallocate(buf_egv,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_occ,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_kwt,persistent=.false.,scalable=.false.)
    !end block muandocc

    ! Note time to get occupation numbers
    t1=mpi_wtime(); mellantid(4)=mellantid(4)+t1-t0; t0=t1

    ! Calculate the density matrix
    !densitymatrix: block

        ! Some work-matrices
        call mem%allocate(wrot,[KS%n_basis,KS%n_basis],persistent=.false.,scalable=.false.)
        call mem%allocate(wdm ,[KS%n_basis,KS%n_basis],persistent=.false.,scalable=.false.)
        call mem%allocate(wegv,[KS%n_basis,KS%n_state],persistent=.false.,scalable=.false.)
        call mem%allocate(wegw,[KS%n_basis,KS%n_state],persistent=.false.,scalable=.false.)
        call mem%allocate(wocc,KS%n_state,persistent=.false.,scalable=.false.)
!allocate(wdd(KS%n_basis,KS%n_basis))
        wrot=0.0_r8
        wegv=0.0_r8
        wegw=0.0_r8
        wocc=0.0_r8
        wdm =0.0_r8

        ! Block-matrix helper
        allocate(BM(p%n_atom,p%n_atom))
        do a1=1,p%n_atom
        do a2=1,p%n_atom
            s1=p%species(a1)
            s2=p%species(a2)
            nb1=basis%species(s1)%n_basis
            nb2=basis%species(s2)%n_basis
            call mem%allocate(BM(a1,a2)%m,[nb1,nb2],persistent=.false.,scalable=.false.)
            BM(a1,a2)%m=0.0_r8
        enddo
        enddo

        ! Make sure the current realspace DM is zeroed properly.
        do ipair=1,rmtx%n_irr_pair
            rmtx%irr_pair(ipair)%densitymatrix=0.0_r8
        enddo

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... building density matrix'
        endif

        tracectr=0.0_r8
        ! Start fiddling with the density matrix.
        do ispin=1,KS%n_spin ! Not sure about where to have spin-loop.
        do ikp=1,KS%n_kpoint
            ! Index of actual k-point
            ikpoint=KS%kpoint(ikp)%irreducible_index
            ! Grab occupation numbers?
            wocc=KS%kpoint(ikp)%occupation(1:KS%n_state,ispin)
            ! Grab eigenvectors
            wegw=0.0_r8
            wegv=0.0_r8
            select type(k=>KS%kpoint(ikp))
            type is(rl_kspace_eigenproblem_kpoint_real)
                wegw=k%eigenvector(:,1:KS%n_state,ispin)
            type is(rl_kspace_eigenproblem_kpoint_complex)
                wegw=k%eigenvector(:,1:KS%n_state,ispin)
            end select

            ! Sqrt and slice the occupation numbers
            wocc=rl_chop(sqrt(wocc),1E-15_r8)
            ! Multipliy in occupation numbers to the eigenvectors, and count number of relevant states
            istatemax=0
            do i=1,KS%n_state
                if ( wocc(i) .gt. 0.0_r8 ) then
                    wegv(:,i)=wegw(:,i)*wocc(i)
                    istatemax=i
                else
                    exit
                endif
            enddo

            ! Now go over all equivalent k-points:
            do ii=1,kmesh%ip(ikpoint)%n_full_point
                ! Identify the k-point we are unfolding to
                jkpoint=kmesh%ip(ikpoint)%index_full_point(ii)
                iop=kmesh%ip(ikpoint)%operation_full_point(ii)
                kvector=kmesh%ap(jkpoint)%r
                ! Now transform the eigenvectors to this k-point. This can be various
                ! levels of complicated.
                wegw=0.0_r8
                if ( iop .eq. 1 ) then
                    ! do nothing, this is just the identity operation
                    wegw(:,1:istatemax)=wegv(:,1:istatemax)
                elseif ( iop .eq. -1 ) then
                    ! this is pure time-reversal. Just conjugate the eigenvectors.
                    wegw(:,1:istatemax)=conjg(wegv(:,1:istatemax))
                elseif ( iop .lt. -1 ) then
                    ! this is a non-trivial operation, followed by conjugate
                    call basis%kspace_rotation_matrix(sym%op(-iop),p,kvector,wrot)
                    ! rotate the eigenvectors
                    call zgemm('N','N', KS%n_basis, istatemax, KS%n_basis, (1.0_r8,0.0_r8), wrot, KS%n_basis, wegv(:,1:istatemax), KS%n_basis, (0.0_r8,0.0_r8), wegw(:,1:istatemax), KS%n_basis)
                    ! and conjugate them
                    wegw(:,1:istatemax)=conjg(wegw(:,1:istatemax))
                elseif ( iop .gt. 1 ) then
                    ! this is a normal rotation with a non-trivial operation.
                    call basis%kspace_rotation_matrix(sym%op(iop),p,kvector,wrot)
                    ! rotate eigenvectors
                    call zgemm('N','N', KS%n_basis, istatemax, KS%n_basis, (1.0_r8,0.0_r8), wrot, KS%n_basis, wegv(:,1:istatemax), KS%n_basis, (0.0_r8,0.0_r8), wegw(:,1:istatemax), KS%n_basis)
                endif
                ! Do the magic outer product guy to get the density matrix
                wdm=0.0_r8
                call zherk('U', 'N', KS%n_basis, istatemax, (1.0_r8,0.0_r8), wegw(:,1:istatemax), KS%n_basis, (0.0_r8,0.0_r8), wdm, KS%n_basis)
                ! Ok, good. now wdm holds the k-space density matrix at this k-point.
                ! fill it out a little.
                do i=1,KS%n_basis
                do j=i+1,KS%n_basis
                    wdm(j,i)=conjg(wdm(i,j))
                enddo
                enddo

                ! ! ! Store the trace, to check I got things right.
                ! do i=1,KS%n_state !basis
                !     tracectr=tracectr+wdm(i,i)*kmesh%ap(jkpoint)%integration_weight
                !     ! c0=dot_product(wegw(:,i),wegw(:,i))*wocc(i)**2
                !     ! write(*,*) i,dot_product(wegw(:,i),wegw(:,i)),dot_product(wegv(:,i),wegv(:,i))
                !     ! if ( abs(c0) .gt. rl_sqtol ) then
                !     ! tracectr=tracectr+wdd(i,i)*kmesh%ap(jkpoint)%integration_weight/c0
                !     ! endif
                !     !write(*,*) i,dot_product(wegw(:,i),wegw(:,i))
                ! enddo

                ! Inverse-transform to realspace, to the irreducible pairs.
                do ipair=1,rmtx%n_irr_pair
                    rvector=rmtx%irr_pair(ipair)%lv
                    kdotr=-rl_twopi*dot_product(rvector,kvector) ! Note negative sign here, as compared to the forward transform above
                    expikr=cmplx(cos(kdotr),sin(kdotr),r8)       ! Get the phase factor
                    expikr=expikr*kmesh%ap(jkpoint)%integration_weight ! This normalizes by number of k-points
                    a1=rmtx%irr_pair(ipair)%a1
                    a2=rmtx%irr_pair(ipair)%a2
                    nb1=rmtx%irr_pair(ipair)%nb1
                    nb2=rmtx%irr_pair(ipair)%nb2
                    i1=basis%offset(a1)+1
                    i2=basis%offset(a1)+nb1
                    j1=basis%offset(a2)+1
                    j2=basis%offset(a2)+nb2
                    rmtx%irr_pair(ipair)%densitymatrix(:,:,ispin)=rmtx%irr_pair(ipair)%densitymatrix(:,:,ispin)+real(wdm(i1:i2,j1:j2)*expikr)
                enddo
            enddo
        enddo
        enddo
        ! Note time to get density matrix, pre-communication
        t1=mpi_wtime(); mellantid(5)=mellantid(5)+t1-t0; t0=t1

        ! Add it up over ranks. I know what the realspace matrix distribution is, so no
        ! real need to worry about that. I will do this really stupidly for now, will
        ! fix later once I know it works.
        f0=0.0_r8
        do ipair=1,rmtx%n_irr_pair
            call mw%allreduce('sum',rmtx%irr_pair(ipair)%densitymatrix)
        enddo

        ! Some cleanup
        call mem%deallocate(wrot,persistent=.false.,scalable=.false.)
        call mem%deallocate(wdm ,persistent=.false.,scalable=.false.)
        call mem%deallocate(wegv,persistent=.false.,scalable=.false.)
        call mem%deallocate(wegw,persistent=.false.,scalable=.false.)
        call mem%deallocate(wocc,persistent=.false.,scalable=.false.)
        do a1=1,p%n_atom
        do a2=1,p%n_atom
            call mem%deallocate(BM(a1,a2)%m,persistent=.false.,scalable=.false.)
        enddo
        enddo
        deallocate(BM)
    !end block densitymatrix

    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

    ! Make note of the total time:
    mellantid(6)=mpi_wtime()-timer
    ! Store the timers in the persistent thing for analysis later.
    tmr%init          =tmr%init          +mellantid(1)
    tmr%forward_ft    =tmr%forward_ft    +mellantid(2)
    tmr%solve         =tmr%solve         +mellantid(3)
    tmr%occnumbers    =tmr%occnumbers    +mellantid(4)
    tmr%densitymatrix =tmr%densitymatrix +mellantid(5)
    tmr%idle          =tmr%idle          +mellantid(6)-sum(mellantid(1:5))
end subroutine

end module
