module rlsy_calculate_density
!!
!! Uses the realspace density matrix to calculate the electron density and it's gradients.
!!
use rlsy_constants, only: r8,rl_iou,rl_exitcode_param,rl_exitcode_memory,rl_hugeint,rl_pi
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_helpers, only: tochar,rl_mom_real, norm2
use rlsy_sorting, only: rl_qsort
use rlsy_integration_grid, only: rl_integration_grid,rl_integration_batch
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_extended_cluster, only: rl_extended_cluster,rl_extended_cluster_hashedspace
use rlsy_realspace_matrix, only: rl_realspace_matrix,rl_realspace_matrix_notdistributed,&
    rl_realspace_matrix_mediumdistributed,rl_realspace_matrix_fulldistributed
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc,&
    rl_kspace_eigenproblem_prob,rl_kspace_eigenproblem_prob_real,&
    rl_kspace_eigenproblem_prob_complex,rl_kspace_eigenproblem_kpoint_real,&
    rl_kspace_eigenproblem_kpoint_complex
use rlsy_verletlist, only: rl_verletbox
use rlsy_electron_density, only: rl_electron_density
use rlsy_timer, only: rl_timer_density
use rlsy_pair_locator, only: rl_pair_locator

implicit none
private

public :: rl_get_density_from_densitymatrix

!> information about an active set of basis functions. All entries are subject to change.
type rl_active_set
    !> how many orders of derivative should be calculated?
    integer :: derivative_level=-rl_hugeint
    !> number of spin channels
    integer :: n_spin=-rl_hugeint

    !> helper arrays to enumerate the currently relevant atoms
    integer :: n_active_atom=-rl_hugeint
    integer, dimension(:), allocatable :: active_atom

    !> basis function index things
    integer :: n_nonzero_basis=-rl_hugeint
    integer, dimension(:), allocatable :: active_atom_offset
    integer, dimension(:), allocatable :: active_atom_nbasis
    integer, dimension(:,:), allocatable :: basis_ind

    !> how many integration points in the batch, and tabulated coordinates.
    integer :: n_integration_point=-rl_hugeint
    real(r8), dimension(:,:), allocatable :: x  !< x/norm(r)
    real(r8), dimension(:,:), allocatable :: y  !< y/norm(r)
    real(r8), dimension(:,:), allocatable :: z  !< z/norm(r)
    real(r8), dimension(:,:), allocatable :: r  !< norm(r)
    real(r8), dimension(:,:), allocatable :: ir !< 1/norm(r)

    !> helper arrays that hold intermediate basis function values
    real(r8), dimension(:,:), allocatable :: f          !< radial function
    real(r8), dimension(:,:), allocatable :: df         !< radial derivative
    real(r8), dimension(:,:), allocatable :: ddf        !< radial second derivative
    real(r8), dimension(:,:), allocatable :: Ylm        !< Ylm spherical harmonic
    real(r8), dimension(:,:), allocatable :: dYlmdx     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: dYlmdy     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: dYlmdz     !< Ylm gradient
    real(r8), dimension(:,:), allocatable :: ddYlmdxdx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdx  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdxdy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdy  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdxdz  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdydz  !< Ylm Hessian
    real(r8), dimension(:,:), allocatable :: ddYlmdzdz  !< Ylm Hessian

    !> flat arrays that hold basis function values and derivatives
    real(r8), dimension(:), allocatable :: basis_V          !< basis fn
    real(r8), dimension(:), allocatable :: basis_dVdx       !< basis gradient
    real(r8), dimension(:), allocatable :: basis_dVdy       !< basis gradient
    real(r8), dimension(:), allocatable :: basis_dVdz       !< basis gradient
    real(r8), dimension(:), allocatable :: basis_ddVdxdx    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdydx    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdzdx    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdxdy    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdydy    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdzdy    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdxdz    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdydz    !< basis Hessian
    real(r8), dimension(:), allocatable :: basis_ddVdzdz    !< basis Hessian

    real(r8), dimension(:,:,:), allocatable :: dm_D
    real(r8), dimension(:,:), allocatable :: dm_DtimesPhi
    real(r8), dimension(:,:), allocatable :: dm_Phi
    real(r8), dimension(:,:), allocatable :: dm_Phi_gx
    real(r8), dimension(:,:), allocatable :: dm_Phi_gy
    real(r8), dimension(:,:), allocatable :: dm_Phi_gz

    real(r8), dimension(:,:), allocatable :: dm_PhiT
    real(r8), dimension(:,:), allocatable :: dm_PhiT_gx
    real(r8), dimension(:,:), allocatable :: dm_PhiT_gy
    real(r8), dimension(:,:), allocatable :: dm_PhiT_gz

    !> density matrix with everything except the relevant entries zeroed out
    real(r8), dimension(:,:), allocatable :: pruned_density_matrix

    !> density
    real(r8), dimension(:,:), allocatable :: rho
    real(r8), dimension(:,:,:), allocatable :: grad_rho
    real(r8), dimension(:,:,:), allocatable :: hess_rho

    contains
        !> initialize temporary storage
        procedure :: setup
        !> reset the work arrays in the active set to nothing.
        procedure :: reset
        !> update list of active basis functions
        procedure :: update_active_basis_functions
        !> tabulate coordinates
        procedure :: tabulate_coordinates
        !> tabulate basis functions
        procedure :: tabulate_basis_functions
        !> fetch pruned density matrix
        procedure :: fetch_pruned_density_matrix
        !> evaluate density
        procedure :: evaluate_density
end type

contains

!> Calulcate the density
subroutine rl_get_density_from_densitymatrix(density,grid,rmtx,KS,p,basis,ec,sym,mw,mem,verbosity,tmr)
    !> electron density
    type(rl_electron_density), intent(inout) :: density
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> distributed realspace matrices
    class(rl_realspace_matrix), allocatable, intent(in) :: rmtx
    !> Kohn-Sham solution thingy handle
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timer
    type(rl_timer_density), intent(inout) :: tmr

    type(rl_active_set) :: active
    type(rl_pair_locator) :: pl
    type(rl_mom_real), dimension(:,:), allocatable :: dm
    real(r8), dimension(KS%n_spin) :: acc_vol
    real(r8), dimension(KS%n_spin) :: acc_charge
    real(r8), dimension(10) :: mellantid
    real(r8) :: timer,t0,t1

    !calcdens: block
    integer :: ibatch,ispin,i,j
    integer :: np
    !cleanup: block
    real(r8) :: scalefactor
    !integer :: i,j

    ! Do some stuff. Maybe allocate temporary storage and things?
    !init: block

        ! Start timers
        timer=mpi_wtime()
        t0=timer
        t1=timer
        mellantid=0.0_r8
        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) ''
            write(rl_iou,*) 'CALCULATING ELECTRON DENSITY'
        endif

        ! Create pair locator? I need a way to know how to take
        ! two generic atoms in the extended cluster and know which
        ! pair that corresponds to. That is indeed a little messy
        ! to get done in a fast enough way, I think. Could become
        ! a bottleneck I ave to revisit.
        call pl%generate(rmtx,KS,ec,irreducible=.false.,mw=mw,mem=mem)

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... built pair locator'
        endif

        ! Create the temporary space for basis functions and stuff like that.
        !@TODO Obviously te derivative level should be input or guesstimated.
        call active%setup(ec,basis,grid,mem,KS%n_spin,derivative_level=1)

        ! Set density to nothing, for now
        select case(active%derivative_level)
        case(0)
            density%irr_rho=0.0_r8
        case(1)
            density%irr_rho=0.0_r8
            density%irr_grad_rho=0.0_r8
        end select
        ! Set total charge counters to zero
        acc_charge=0.0_r8
        acc_vol=0.0_r8

        ! Note how long initialization took. Includes allocations and the pair locator.
        t1=mpi_wtime(); mellantid(1)=t1-t0; t0=t1
        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... set up temporary storage'
        endif

        ! Fetch the full density density matrix to all ranks. Danger zone.
        call fetch_full_density_matrix(rmtx,KS,p,basis,sym,dm,mw,mem)

        ! Note how long grabbing densitymatrix took
        t1=mpi_wtime(); mellantid(2)=t1-t0; t0=t1
        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... fetched full density matrix'
        endif
    !end block init

    ! begin operation: 'actual operation'
    !calcdens: block

        do ibatch=1,grid%n_irr_batch

            ! First thing, reset temporary storage.
            call active%reset(ec,basis,grid,grid%irr_batch(ibatch),mem)
            t1=mpi_wtime(); mellantid(3)=mellantid(3)+t1-t0; t0=t1

            ! Tabulate some coordinates
            call active%tabulate_coordinates(ec,grid%irr_batch(ibatch))
            t1=mpi_wtime(); mellantid(4)=mellantid(4)+t1-t0; t0=t1

            ! Grab the list of active basis functions
            call active%update_active_basis_functions(ec,p,basis)
            t1=mpi_wtime(); mellantid(5)=mellantid(5)+t1-t0; t0=t1

            ! Evaluate basis functions?
            call active%tabulate_basis_functions(basis)
            t1=mpi_wtime(); mellantid(6)=mellantid(6)+t1-t0; t0=t1

            ! Fetch density matrix.
            call active%fetch_pruned_density_matrix(ec,rmtx,basis,pl,dm)
            t1=mpi_wtime(); mellantid(7)=mellantid(7)+t1-t0; t0=t1

            ! Calculate actual density
            call active%evaluate_density()

            ! Accumulate total charge
            np=grid%irr_batch(ibatch)%n_point
            do ispin=1,KS%n_spin
                acc_charge(ispin)=acc_charge(ispin)+sum(active%rho(1:np,ispin)*grid%irr_batch(ibatch)%partition_function*grid%irr_batch(ibatch)%integration_weight)
                acc_vol(ispin)=acc_vol(ispin)+sum(grid%irr_batch(ibatch)%partition_function*grid%irr_batch(ibatch)%integration_weight)
            enddo

            ! Store it in the right place
            !@TODO Actually, why not store it in batches directly? No indexing needed, and never gets lost.
            do i=1,grid%irr_batch(ibatch)%n_point
                j=grid%irr_batch(ibatch)%semilocal_irr_offset+i
                select case(active%derivative_level)
                case(0)
                    density%irr_rho(j,:)=active%rho(i,:)
                case(1)
                    density%irr_rho(j,:)=active%rho(i,:)
                    density%irr_grad_rho(:,j,:)=active%grad_rho(:,i,:)
                end select
            enddo
            t1=mpi_wtime(); mellantid(8)=mellantid(8)+t1-t0; t0=t1
        enddo
    !end block calcdens

    ! Free space and communicate
    !cleanup: block

        ! Clean temporary storage neatly
        do j=1,size(dm,2)
        do i=1,size(dm,1)
            call mem%deallocate(dm(i,j)%m,persistent=.true.,scalable=.false.)
        enddo
        enddo

        ! Add together integrated charge
        call mw%allreduce('sum',acc_charge)
        call mw%allreduce('sum',acc_vol)

        ! Not sure about this, but does probably not hurt either.
        scalefactor=density%n_electron/acc_charge(1)


        ! Here I should communicate the density properly.
        select case(active%derivative_level)
        case(0)
            call grid%ml%allreduce('sum',density%irr_rho)
            if ( abs(scalefactor-1.0_r8) .gt. 1E-13_r8 ) then
                density%irr_rho=density%irr_rho*scalefactor
            endif
        case(1)
            call grid%ml%allreduce('sum',density%irr_rho)
            call grid%ml%allreduce('sum',density%irr_grad_rho)
            if ( abs(scalefactor-1.0_r8) .gt. 1E-13_r8 ) then
                density%irr_rho=density%irr_rho*scalefactor
                density%irr_grad_rho=density%irr_grad_rho*scalefactor
            endif
        end select



    !end block cleanup

    ! And finally, sort out the timings
    mellantid(10)=mpi_wtime()-timer                ! total timer
    mellantid(9)=mellantid(10)-sum(mellantid(1:8)) ! time idle
    tmr%init        = tmr%init +        mellantid(1)
    tmr%idle        = tmr%idle +        mellantid(9)
    tmr%total       = tmr%total +       mellantid(10)
    tmr%fetchdm     = tmr%fetchdm +     mellantid(2)
    tmr%reset       = tmr%reset +       mellantid(3)
    tmr%tab_coord   = tmr%tab_coord +   mellantid(4)
    tmr%prune_basis = tmr%prune_basis + mellantid(5)
    tmr%tab_basis   = tmr%tab_basis +   mellantid(6)
    tmr%prune_dm    = tmr%prune_dm +    mellantid(7)
    tmr%matrixop    = tmr%matrixop +    mellantid(8)

    if ( verbosity .gt. 0 ) then
        t1=mpi_wtime()
        write(*,*) '... calculated density (',tochar(t1-timer),'s)'
        write(*,*) '  electron count:',acc_charge
        write(*,*) '   proper target:',KS%n_electron
        write(rl_iou,'(1X,A,F15.6)') '  timer: init        ',mellantid(1)
        write(rl_iou,'(1X,A,F15.6)') '  timer: fullDM      ',mellantid(2)
        write(rl_iou,'(1X,A,F15.6)') '  timer: reset       ',mellantid(3)
        write(rl_iou,'(1X,A,F15.6)') '  timer: prune atom  ',mellantid(4)
        write(rl_iou,'(1X,A,F15.6)') '  timer: prune basis ',mellantid(5)
        write(rl_iou,'(1X,A,F15.6)') '  timer: eval basis  ',mellantid(6)
        write(rl_iou,'(1X,A,F15.6)') '  timer: prune DM    ',mellantid(7)
        write(rl_iou,'(1X,A,F15.6)') '  timer: matrixop    ',mellantid(8)
        t0=t1
    endif

    ! Check for stupidity
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

end subroutine

!@TODO make a destructor

!> Fetch the full realspace density matrix into a buffer.
subroutine fetch_full_density_matrix(rmtx,KS,p,basis,sym,dm,mw,mem)
    !> distributed realspace matrices
    class(rl_realspace_matrix), allocatable, intent(in) :: rmtx
    !> Kohn-Sham solution thingy handle
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> full realspace density matrix
    type(rl_mom_real), dimension(:,:), allocatable, intent(out) :: dm
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    integer, parameter :: safenumentries=1000000 ! or some other arbitrary number
    type(rl_mom_real), dimension(:,:), allocatable :: rbuf
    real(r8), dimension(:,:,:), allocatable :: buf_comm
    integer, dimension(:,:), allocatable :: buf_ia
    integer :: n_pair_per_iter,niter,npair
    !init: block
    integer :: ispin,ipar,ipair,jpair,iop
    integer :: l,s1,s2,nb1,nb2

    !fixdm: block
    !integer :: plo,phi,poff,a1,a2,s1,s2,nb1,nb2
    integer :: plo,phi,poff,a1,a2
    !integer :: iter,ipair,jpair,iop,ispin,i
    integer :: iter,i
    !finalize: block
    !    integer :: s1,s2

    !init: block

        ! First things first: look for an easy exit. In case the Hamiltonian is not
        ! distributed at all, and there is no symmetry, there is little use bothering
        ! with being parallel or other fancy stuff. I will just copy the data.
        if ( sym%n_operation .eq. 1 ) then
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! Yup, this is the trivial case. Just do a copy.
            ! Make a little space
            allocate(dm(m%n_full_pair,KS%n_spin))
            do ispin=1,KS%n_spin
            do ipair=1,m%n_full_pair
                call mem%allocate(dm(ipair,ispin)%m,[m%full_pair(ipair)%nb1,m%full_pair(ipair)%nb2],persistent=.true.,scalable=.false.)
                dm(ipair,ispin)%m=0.0_r8
            enddo
            enddo
            ! Populate it.
            do ispin=1,KS%n_spin
            do ipair=1,rmtx%n_full_pair
                jpair=rmtx%full_pair(ipair)%index_unique
                iop=rmtx%full_pair(ipair)%index_operation
                select case(iop)
                case(1)
                    dm(ipair,ispin)%m=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin)
                case(-1)
                    dm(ipair,ispin)%m=transpose(rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin))
                end select
            enddo
            enddo
            ! And we are done and can return early!
            return
        end select
        endif

        ! So, if we make it here, it was not the trivial case and life is a little
        ! annoying. First get the global number of pairs.
        select type(rmtx)
        type is(rl_realspace_matrix_notdistributed)
            npair=rmtx%n_full_pair
        type is(rl_realspace_matrix_mediumdistributed)
            npair=rmtx%n_full_pair_global
        type is(rl_realspace_matrix_fulldistributed)
            npair=rmtx%n_full_pair_global
        end select

        ! Largest number of basis functions per species, need to know this
        ! to decide on communication scheme.
        l=0
        do s1=1,p%n_species
            l=max(l,basis%species(s1)%n_basis)
        enddo

        ! Then decide how communication is going to happen, since everything
        ! is distributed all over the place.
        n_pair_per_iter=int(anint(safenumentries/real(KS%n_spin*l**2,r8)))
        n_pair_per_iter=max(n_pair_per_iter,1)
        n_pair_per_iter=min(n_pair_per_iter,npair)
        niter=0
        do
            if ( niter*n_pair_per_iter .ge. npair ) then
                exit
            else
                niter=niter+1
            endif
        enddo

        ! And make a little space for rotated buffers.
        allocate(rbuf(p%n_species,p%n_species))
        do s1=1,p%n_species
        do s2=1,p%n_species
            call mem%allocate(rbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
            rbuf(s1,s2)%m=0.0_r8
        enddo
        enddo

        ! And space for communication buffers
        call mem%allocate(buf_comm,[l*l,KS%n_spin,n_pair_per_iter],persistent=.false.,scalable=.false.)
        call mem%allocate(buf_ia,[2,n_pair_per_iter],persistent=.false.,scalable=.false.)
        buf_comm=0.0_r8
        buf_ia=0

        ! And, most importantly, space for the density-matrix container
        allocate(dm(npair,KS%n_spin))
    !end block init

    !fixdm: block

        iterloop: do iter=1,niter
            ! Reset communication buffers
            buf_comm=0.0_r8
            buf_ia=0
            ! Set relevant pair indices for this iteration
            plo=(iter-1)*n_pair_per_iter+1
            phi=min(iter*n_pair_per_iter,npair)
            poff=(iter-1)*n_pair_per_iter

            ! Now it depends a bit on distribution how we will collect things.
            select type(m=>rmtx)
            type is(rl_realspace_matrix_notdistributed)
                do ipair=plo,phi
                    if ( mod(ipair,mw%n) .ne. mw%r ) cycle
                    jpair=m%full_pair(ipair)%index_unique
                    iop=m%full_pair(ipair)%index_operation
                    a1=m%full_pair(ipair)%a1
                    a2=m%full_pair(ipair)%a2
                    s1=p%species(a1)
                    s2=p%species(a2)
                    nb1=m%full_pair(ipair)%nb1
                    nb2=m%full_pair(ipair)%nb2

                    do ispin=1,KS%n_spin
                        select case(sym%n_operation)
                        case(1)
                            if ( iop .eq. 1 ) then
                                rbuf(s1,s2)%m=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin)
                            else
                                rbuf(s1,s2)%m=transpose(rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin))
                            endif
                        case default
                            if ( iop .gt. 0 ) then
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(iop),&
                                    forward=.true.,transposition=.false.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            else
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(-iop),&
                                    forward=.true.,transposition=.true.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            endif
                        end select
                        buf_comm(1:nb1*nb2,ispin,ipair-poff)=reshape(rbuf(s1,s2)%m,[nb1*nb2])
                    enddo
                    buf_ia(:,ipair-poff)=[nb1,nb2]
                enddo
                ! Reduce over all ranks
                call mw%allreduce('sum',buf_comm)
                call mw%allreduce('sum',buf_ia)
            type is(rl_realspace_matrix_mediumdistributed)
                ! So in this case the full thing is distributed. I have to think
                ! a little about how to be smart about this. Hmmm. Think think.
                do i=1,m%n_full_pair
                    ipair=i+m%full_offset
                    if ( ipair .lt. plo ) cycle
                    if ( ipair .gt. phi ) cycle
                    jpair=m%full_pair(i)%index_unique
                    iop=m%full_pair(i)%index_operation
                    nb1=m%full_pair(i)%nb1
                    nb2=m%full_pair(i)%nb2
                    a1=m%full_pair(i)%a1
                    a2=m%full_pair(i)%a2
                    s1=p%species(a1)
                    s2=p%species(a2)
                    do ispin=1,KS%n_spin
                        select case(sym%n_operation)
                        case(1)
                            if ( iop .eq. 1 ) then
                                rbuf(s1,s2)%m=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin)
                            else
                                rbuf(s1,s2)%m=transpose(rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin))
                            endif
                        case default
                            if ( iop .gt. 0 ) then
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(iop),&
                                    forward=.true.,transposition=.false.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            else
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(-iop),&
                                    forward=.true.,transposition=.true.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            endif
                        end select
                        buf_comm(1:nb1*nb2,ispin,ipair-poff)=reshape(rbuf(s1,s2)%m,[nb1*nb2])
                    enddo
                    buf_ia(:,ipair-poff)=[nb1,nb2]
                enddo
                ! And reduce over ranks, but only the local ranks so we don't get annoying
                ! everyone-to-everyone communication, if possible.
                select type(KS)
                type is(rl_kspace_eigenproblem_multiproc)
                    call KS%ml%allreduce('sum',buf_comm)
                    call KS%ml%allreduce('sum',buf_ia)
                class default
                    call rl_stop_gracefully(['Distributed Hamiltonian and multiproc ELSI at the same time, strange'],rl_exitcode_param,mw%comm)
                end select
            type is(rl_realspace_matrix_fulldistributed)
                ! Same thing as above, only difference is that we reduce with the world communicator
                ! at the end.
                do i=1,m%n_full_pair
                    ipair=i+m%full_offset
                    if ( ipair .lt. plo ) cycle
                    if ( ipair .gt. phi ) cycle
                    jpair=m%full_pair(i)%index_unique
                    iop=m%full_pair(i)%index_operation
                    nb1=m%full_pair(i)%nb1
                    nb2=m%full_pair(i)%nb2
                    a1=m%full_pair(i)%a1
                    a2=m%full_pair(i)%a2
                    s1=p%species(a1)
                    s2=p%species(a2)
                    do ispin=1,KS%n_spin
                        select case(sym%n_operation)
                        case(1)
                            if ( iop .eq. 1 ) then
                                rbuf(s1,s2)%m=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin)
                            else
                                rbuf(s1,s2)%m=transpose(rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin))
                            endif
                        case default
                            if ( iop .gt. 0 ) then
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(iop),&
                                    forward=.true.,transposition=.false.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            else
                                call rmtx%irr_pair(jpair)%rotate(basis,p,sym%op(-iop),&
                                    forward=.true.,transposition=.true.,&
                                    original_block=rmtx%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                                    rotated_block=rbuf(s1,s2)%m)
                            endif
                        end select
                        buf_comm(1:nb1*nb2,ispin,ipair-poff)=reshape(rbuf(s1,s2)%m,[nb1*nb2])
                    enddo
                    buf_ia(:,ipair-poff)=[nb1,nb2]
                enddo
                ! And reduce over ranks,
                call mw%allreduce('sum',buf_comm)
                call mw%allreduce('sum',buf_ia)
            end select

            ! Now we can store things
            do ipair=plo,phi
                nb1=buf_ia(1,ipair-poff)
                nb2=buf_ia(2,ipair-poff)
                do ispin=1,KS%n_spin
                    call mem%allocate(dm(ipair,ispin)%m,[nb1,nb2],persistent=.true.,scalable=.false.)
                    dm(ipair,ispin)%m=reshape( buf_comm(1:nb1*nb2,ispin,ipair-poff),[nb1,nb2])
                enddo
            enddo
        enddo iterloop
    !end block fixdm

    !finalize: block
        ! Cleanup temporary things
        do s1=1,p%n_species
        do s2=1,p%n_species
            call mem%deallocate(rbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
        enddo
        enddo
        deallocate(rbuf)
        call mem%deallocate(buf_comm,persistent=.false.,scalable=.false.)
        call mem%deallocate(buf_ia,persistent=.false.,scalable=.false.)
    !end block finalize
end subroutine

!> evaluate the density
subroutine evaluate_density(active)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active

    integer :: i,ispin,nb,n,m

    select case(active%derivative_level)
    case(0)
        n=size(active%dm_D,1)
        m=size(active%dm_Phi,2)
        do ispin=1,active%n_spin
            !call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            call dgemm('N','N',n,m,n,1.0_r8,active%dm_D(:,:,ispin),n,active%dm_Phi,n,0.0_r8,active%dm_DtimesPhi,n)
            !active%dm_DtimesPhi=matmul(active%dm_D(:,:,ispin),active%dm_Phi)
            do i=1,active%n_integration_point
                active%rho(i,ispin)=dot_product(active%dm_Phi(:,i),active%dm_DtimesPhi(:,i))
            enddo
        enddo
        deallocate(active%dm_D)
        deallocate(active%dm_DtimesPhi)
        deallocate(active%dm_Phi)
        deallocate(active%dm_PhiT)
    case(1)
        n=size(active%dm_D,1)
        m=size(active%dm_Phi,2)
        do ispin=1,active%n_spin
            ! density
            !@TODOdsymm?
            call dgemm('N','N',n,m,n,1.0_r8,active%dm_D(:,:,ispin),n,active%dm_Phi,n,0.0_r8,active%dm_DtimesPhi,n)
            do i=1,active%n_integration_point
                active%rho(i,ispin)=dot_product(active%dm_Phi(:,i),active%dm_DtimesPhi(:,i))
            enddo

            ! gradient of density.
            !@TODO Could be one huge dsymm.
            call dgemm('N','N',n,m,n,1.0_r8,active%dm_D(:,:,ispin),n,active%dm_Phi_gx,n,0.0_r8,active%dm_DtimesPhi,n)
            do i=1,active%n_integration_point
                active%grad_rho(1,i,ispin)=2*dot_product(active%dm_Phi(:,i),active%dm_DtimesPhi(:,i))
            enddo
            call dgemm('N','N',n,m,n,1.0_r8,active%dm_D(:,:,ispin),n,active%dm_Phi_gy,n,0.0_r8,active%dm_DtimesPhi,n)
            do i=1,active%n_integration_point
                active%grad_rho(2,i,ispin)=2*dot_product(active%dm_Phi(:,i),active%dm_DtimesPhi(:,i))
            enddo
            call dgemm('N','N',n,m,n,1.0_r8,active%dm_D(:,:,ispin),n,active%dm_Phi_gz,n,0.0_r8,active%dm_DtimesPhi,n)
            do i=1,active%n_integration_point
                active%grad_rho(3,i,ispin)=2*dot_product(active%dm_Phi(:,i),active%dm_DtimesPhi(:,i))
            enddo
        enddo
        deallocate(active%dm_D)
        deallocate(active%dm_DtimesPhi)
        deallocate(active%dm_Phi)
        deallocate(active%dm_Phi_gx)
        deallocate(active%dm_Phi_gy)
        deallocate(active%dm_Phi_gz)
        deallocate(active%dm_PhiT)
        deallocate(active%dm_PhiT_gx)
        deallocate(active%dm_PhiT_gy)
        deallocate(active%dm_PhiT_gz)
    end select
end subroutine

!> get the pruned density matrix. This is actually a bottleneck now, have to think.
subroutine fetch_pruned_density_matrix(active,ec,rmtx,basis,pl,dm)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> matrix
    type(rl_realspace_matrix), intent(in) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> pair locator
    type(rl_pair_locator), intent(in) :: pl
    !> realspace density matrix, sorted into pairs
    type(rl_mom_real), dimension(:,:), intent(in) :: dm

    integer :: i,j,ie,je,ipair,ispin
    integer :: b1,b2,l1,l2,o1,o2
    integer :: l

    active%pruned_density_matrix=0.0_r8
    do i=1,active%n_active_atom
        ie=active%active_atom(i)
        do j=i,active%n_active_atom
            je=active%active_atom(j)
            ipair=pl%locate(ie,je)
            if ( ipair .lt. 0 ) then
                !write(*,*) 'NEGATIVE IPAIR'
                cycle
            endif

            do l1=1,active%active_atom_nbasis(i)
            do l2=1,active%active_atom_nbasis(j)
                b1=active%active_atom_offset(i)+l1
                b2=active%active_atom_offset(j)+l2
                l=(b2-1)*active%n_nonzero_basis+b1
                ! Don't forget ispin here.
                !active%pruned_density_matrix(l,1)=dm(ipair,1)%m(l1,l2)
                do ispin=1,active%n_spin
                    active%dm_D(b1,b2,ispin)=dm(ipair,ispin)%m(l1,l2)
                enddo
            enddo
            enddo
        enddo
    enddo
    ! Then fill out the transpose, to be on the safe side.
    do ispin=1,active%n_spin
        do b1=1,active%n_nonzero_basis
        do b2=b1+1,active%n_nonzero_basis
            active%dm_D(b2,b1,ispin)=active%dm_D(b1,b2,ispin)
            ! i=(b2-1)*active%n_nonzero_basis+b1
            ! j=(b1-1)*active%n_nonzero_basis+b2
            ! active%pruned_density_matrix(j,1)=active%pruned_density_matrix(i,1)
        enddo
        enddo
    enddo
end subroutine

!> tabulate the coordinates of all atoms with respect to all points
subroutine tabulate_coordinates(active,ec,batch)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> batch information
    class(rl_integration_batch), intent(in) :: batch

    real(r8), dimension(3) :: v,w
    real(r8) :: f0,f1,tol
    integer :: np,na,i,j,k

    ! Some shorthand
    np=batch%n_point
    na=active%n_active_atom
    tol=1E-10_r8

    ! calculate coordinates
    do i=1,na
        k=active%active_atom(i)
        w=ec%cartesian_coordinate(:,k)
        do j=1,np
            v=batch%folded_coordinate(:,j)-w
            f0=norm2(v)
            if ( f0 .gt. tol ) then
                f1=1.0_r8/f0
                v=v*f1
                active%x(j,i)=v(1)
                active%y(j,i)=v(2)
                active%z(j,i)=v(3)
                active%r(j,i)=f0
                active%ir(j,i)=f1
            else
                active%x(j,i)=0.0_r8
                active%y(j,i)=0.0_r8
                active%z(j,i)=0.0_r8
                active%r(j,i)=f0
                active%ir(j,i)=0.0_r8
            endif
        enddo
    enddo
end subroutine

!> convert the list of active atoms to a list of active basis functions
subroutine update_active_basis_functions(active,ec,p,basis)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set information
    type(rl_lcao_basis_set), intent(in) :: basis


    integer :: i,j,l
    integer :: is,iu,ie,nb

    ! This could probably be made faster by first generating a possible pool for
    ! all the batches or something. Don't think this is a bottleneck though.
    active%basis_ind=0
    l=0
    do i=1,active%n_active_atom
        ie=active%active_atom(i)
        iu=ec%index_unit_cell(ie)
        is=p%species(iu)
        nb=basis%species(is)%n_basis
        active%active_atom_offset(i)=l
        active%active_atom_nbasis(i)=nb
        do j=1,nb
            l=l+1
            active%basis_ind(:,l)=[basis%species(is)%basis_fn(j),basis%species(is)%basis_ang(j)]
        enddo
    enddo
    active%n_nonzero_basis=l

    ! Equipped with this information, I can make some temporary space:
    select case(active%derivative_level)
    case(0)
        allocate(active%dm_D(active%n_nonzero_basis,active%n_nonzero_basis,active%n_spin))
        allocate(active%dm_DtimesPhi(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_Phi(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_PhiT(active%n_integration_point,active%n_nonzero_basis))
    case(1)
        allocate(active%dm_D(active%n_nonzero_basis,active%n_nonzero_basis,active%n_spin))
        allocate(active%dm_DtimesPhi(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_Phi(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_Phi_gx(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_Phi_gy(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_Phi_gz(active%n_nonzero_basis,active%n_integration_point))
        allocate(active%dm_PhiT(active%n_integration_point,active%n_nonzero_basis))
        allocate(active%dm_PhiT_gx(active%n_integration_point,active%n_nonzero_basis))
        allocate(active%dm_PhiT_gy(active%n_integration_point,active%n_nonzero_basis))
        allocate(active%dm_PhiT_gz(active%n_integration_point,active%n_nonzero_basis))
    end select
end subroutine

!> tabulate the values of the wave functions over the batch
subroutine tabulate_basis_functions(active,basis)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> basis set information
    type(rl_lcao_basis_set), intent(in) :: basis

    integer :: i,j,i1,i2,ir,ia,np

    np=active%n_integration_point
    ! Do it slightly different depending on how many derivatives we need.
    select case(active%derivative_level)
    case(0)
        ! No derivatives, just function values
        do i=1,active%n_active_atom
            ! Evaluate radial functions
            do ir=1,basis%n_radial
                call basis%radial(ir)%val( active%r(1:np,i),active%f(1:np,ir) )
            enddo
            ! Evaluate spherical harmonics
            do ia=1,basis%n_angular
                call basis%Ylm( active%x(1:np,i),active%y(1:np,i),active%z(1:np,i),ia,active%Ylm(1:np,ia) )
            enddo
            ! Combine basis functions to the full ones
            do j=active%active_atom_offset(i)+1,active%active_atom_offset(i)+active%active_atom_nbasis(i)
                ir = active%basis_ind(1,j) ! radial index of basis fn j
                ia = active%basis_ind(2,j) ! angular index of basis fn j
                ! where to store values in flat array
                i1=(j-1)*np+1
                i2=j*np
                !active%basis_v(i1:i2)=active%f(1:np,ir) * active%Ylm(1:np,ia)
                active%dm_PhiT(:,j)=active%f(1:np,ir) * active%Ylm(1:np,ia)
            enddo
        enddo
        active%dm_Phi=transpose(active%dm_PhiT)
    case(1)
        ! Get values + gradient
        do i=1,active%n_active_atom

            ! Evaluate radial functions + first derivative
            do ir=1,basis%n_radial
                call basis%radial(ir)%val_der( active%r(1:np,i),active%f(1:np,ir),active%df(1:np,ir) )
            enddo
            ! Evaluate spherical harmonics + gradient
            do ia=1,basis%n_angular
                call basis%Ylm_grad( active%x(1:np,i),active%y(1:np,i),active%z(1:np,i),active%ir(1:np,i),ia,&
                active%Ylm(1:np,ia),active%dYlmdx(1:np,ia),active%dYlmdy(1:np,ia),active%dYlmdz(1:np,ia) )
            enddo

            ! Loop over the basis functions on this atom.
            do j=active%active_atom_offset(i)+1,active%active_atom_offset(i)+active%active_atom_nbasis(i)
                ir = active%basis_ind(1,j) ! radial index of basis fn j
                ia = active%basis_ind(2,j) ! angular index of basis fn j
                ! where to store values in flat array
                i1=(j-1)*np+1
                i2=j*np
                ! Store function values
                active%dm_PhiT(:,j)=active%f(1:np,ir) * active%Ylm(1:np,ia)
                ! Store gradient
                active%dm_PhiT_gx(:,j)= active%x(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdx(1:np,ia)
                active%dm_PhiT_gy(:,j)= active%y(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdy(1:np,ia)
                active%dm_PhiT_gz(:,j)= active%z(1:np,i)*active%df(1:np,ir)*active%Ylm(1:np,ia) + active%f(1:np,ir)*active%dYlmdz(1:np,ia)
                ! active%dm_PhiT_gx(:,j)= active%dYlmdx(1:np,ia)
                ! active%dm_PhiT_gy(:,j)= active%dYlmdy(1:np,ia)
                ! active%dm_PhiT_gz(:,j)= active%dYlmdz(1:np,ia)

            enddo
        enddo
        active%dm_Phi=transpose(active%dm_PhiT)
        active%dm_Phi_gx=transpose(active%dm_PhiT_gx)
        active%dm_Phi_gy=transpose(active%dm_PhiT_gy)
        active%dm_Phi_gz=transpose(active%dm_PhiT_gz)
    case(2)
        ! Get values + gradient + hessian
    end select
end subroutine

!> reset the active set for work on a new batch
subroutine reset(active,ec,basis,grid,batch,mem)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> basis set information
    type(rl_lcao_basis_set), intent(in) :: basis
    !> grid information
    type(rl_integration_grid), intent(in) :: grid
    !> batch information
    class(rl_integration_batch), intent(in) :: batch
    !> memory tracker
    type(rl_memtracker), intent(in) :: mem

    ! Zero all basis functions
    active%n_nonzero_basis    = 0
    active%active_atom_offset = 0
    active%active_atom_nbasis = 0
    active%basis_ind          = 0

    ! Make a note of the number of integration points
    active%n_integration_point=batch%n_point

    ! set the active atoms to the ones determined by the batch?
    active%n_active_atom=batch%n_relevant_ext_atom
    active%active_atom = 0
    active%active_atom(1:active%n_active_atom)=batch%relevant_ext_atom(1:active%n_active_atom)

    ! space for coordinates
    active%x  = 0.0_r8
    active%y  = 0.0_r8
    active%z  = 0.0_r8
    active%r  = 0.0_r8
    active%ir = 0.0_r8
    select case(active%derivative_level)
    case(0)
        active%rho     = 0.0_r8
        ! space for basis function values
        active%basis_V = 0.0_r8
        active%f       = 0.0_r8
        active%Ylm     = 0.0_r8
    case(1)
        active%rho     = 0.0_r8
        active%grad_rho= 0.0_r8
        ! space for basis function values
        active%basis_V    = 0.0_r8
        active%basis_dVdx = 0.0_r8
        active%basis_dVdy = 0.0_r8
        active%basis_dVdz = 0.0_r8
        ! Intermediate things for derivatives
        active%f      = 0.0_r8
        active%df     = 0.0_r8
        active%Ylm    = 0.0_r8
        active%dYlmdx = 0.0_r8
        active%dYlmdy = 0.0_r8
        active%dYlmdz = 0.0_r8
    case(2)
        ! space for basis function values
        active%basis_V       = 0.0_r8
        active%basis_dVdx    = 0.0_r8
        active%basis_dVdy    = 0.0_r8
        active%basis_dVdz    = 0.0_r8
        active%basis_ddVdxdx = 0.0_r8
        active%basis_ddVdydx = 0.0_r8
        active%basis_ddVdzdx = 0.0_r8
        active%basis_ddVdxdy = 0.0_r8
        active%basis_ddVdydy = 0.0_r8
        active%basis_ddVdzdy = 0.0_r8
        active%basis_ddVdxdz = 0.0_r8
        active%basis_ddVdydz = 0.0_r8
        active%basis_ddVdzdz = 0.0_r8
        ! Intermediate things for derivatives
        active%f         = 0.0_r8
        active%df        = 0.0_r8
        active%ddf       = 0.0_r8
        active%Ylm       = 0.0_r8
        active%dYlmdx    = 0.0_r8
        active%dYlmdy    = 0.0_r8
        active%dYlmdz    = 0.0_r8
        active%ddYlmdxdx = 0.0_r8
        active%ddYlmdydx = 0.0_r8
        active%ddYlmdzdx = 0.0_r8
        active%ddYlmdxdy = 0.0_r8
        active%ddYlmdydy = 0.0_r8
        active%ddYlmdzdy = 0.0_r8
        active%ddYlmdxdz = 0.0_r8
        active%ddYlmdydz = 0.0_r8
        active%ddYlmdzdz = 0.0_r8
    end select
end subroutine

!> make some temporary space
subroutine setup(active,ec,basis,grid,mem,n_spin,derivative_level)
    !> subset of basis functions and atoms
    class(rl_active_set), intent(inout) :: active
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> basis set information
    type(rl_lcao_basis_set), intent(in) :: basis
    !> grid information
    type(rl_integration_grid), intent(in) :: grid
    !> memory tracker
    type(rl_memtracker), intent(in) :: mem
    !> number of spin channels
    integer, intent(in) :: n_spin
    !> how many derivatives are to be taken
    integer, intent(in) :: derivative_level

    ! Set some constants
    active%derivative_level=derivative_level
    active%n_spin = n_spin

    ! Space for tracking relevant atoms and basis functions
    allocate(active%active_atom(        grid%max_n_atom_per_batch) )
    allocate(active%active_atom_offset( grid%max_n_atom_per_batch) )
    allocate(active%active_atom_nbasis( grid%max_n_atom_per_batch) )
    allocate(active%basis_ind(        2,grid%max_n_basis_per_batch))

    ! Space for coordinates
    allocate(active%x (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(active%y (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(active%z (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(active%r (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(active%ir(grid%max_n_point_per_batch,grid%max_n_atom_per_batch))


    ! space for pruned density matrix @TODO don't forget nspin here, eventually.
    allocate(active%pruned_density_matrix(grid%max_n_basis_per_batch**2,active%n_spin))

    select case(active%derivative_level)
    case(0)
        ! Space for density
        allocate(active%rho(grid%max_n_point_per_batch,active%n_spin))
        ! space for basis function values
        allocate(active%basis_V      (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        ! Intermediate things for derivatives
        allocate(active%f      ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%Ylm    ( grid%max_n_point_per_batch,basis%n_angular ))
    case(1)
        ! Space for density
        allocate(active%rho(grid%max_n_point_per_batch,active%n_spin))
        ! Space for gradient
        allocate(active%grad_rho(3,grid%max_n_point_per_batch,active%n_spin))
        ! space for basis function values
        allocate(active%basis_V      (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdx   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdy   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdz   (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        ! Intermediate things for derivatives
        allocate(active%f      ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%df     ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%Ylm    ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdx ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdy ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdz ( grid%max_n_point_per_batch,basis%n_angular ))
    case(2)
        ! space for basis function values
        allocate(active%basis_V       (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdx    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdy    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_dVdz    (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdxdx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdydx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdzdx (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdxdy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdydy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdzdy (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdxdz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdydz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        allocate(active%basis_ddVdzdz (grid%max_n_point_per_batch*grid%max_n_basis_per_batch))
        ! Intermediate things for derivatives
        allocate(active%f         ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%df        ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%ddf       ( grid%max_n_point_per_batch,basis%n_radial  ))
        allocate(active%Ylm       ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdx    ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdy    ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%dYlmdz    ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdxdx ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdydx ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdzdx ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdxdy ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdydy ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdzdy ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdxdz ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdydz ( grid%max_n_point_per_batch,basis%n_angular ))
        allocate(active%ddYlmdzdz ( grid%max_n_point_per_batch,basis%n_angular ))
    end select
end subroutine

end module
