module rlsy_calculate_overlap
!!
!! Calculates the realspace overlap matrix
!!
!! @TODO insert ascii drawing
!! @TODO move integration helper
!! @TODO flatten arrays to avoid reallocation.
!!
use rlsy_constants, only: r8,rl_iou,rl_exitcode_param,rl_exitcode_memory,rl_hugeint,rl_pi,rl_exitcode_symmetry
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_helpers, only: tochar,rl_mom_real,norm2
!use rlsy_sorting, only: rl_qsort
use rlsy_integration_grid, only: rl_integration_grid,rl_integration_batch
use rlsy_crystalstructure, only: rl_crystalstructure
!use rlsy_spacegroup, only: rl_spacegroup
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
use rlsy_timer, only: rl_timer_overlap
use rlsy_pair_locator, only: rl_pair_locator

implicit none
private

public :: rl_calculate_overlapmatrix

type integrationhelper
    !> how many pairs are there, globally
    integer :: n_pair_global=-rl_hugeint
    !> I need to be able to quickly look up what the block dimensions are, per pair
    integer, dimension(:), allocatable :: nbasis_atom1
    integer, dimension(:), allocatable :: nbasis_atom2
    !> which way sould I update overlap?
    integer :: update_mode=-rl_hugeint

    !> helper arrays to enumerate the currently relevant atoms
    integer :: n_active_atom=-rl_hugeint
    integer, dimension(:), allocatable :: active_atom

    !> basis function index things
    integer :: n_nonzero_basis=-rl_hugeint
    integer, dimension(:), allocatable :: active_atom_offset
    integer, dimension(:), allocatable :: active_atom_nbasis
    integer, dimension(:,:), allocatable :: basis_ind

    ! Space for coordinate things
    integer :: n_integration_point=-rl_hugeint
    real(r8), dimension(:,:), allocatable :: x  !< x/norm(r)
    real(r8), dimension(:,:), allocatable :: y  !< y/norm(r)
    real(r8), dimension(:,:), allocatable :: z  !< z/norm(r)
    real(r8), dimension(:,:), allocatable :: r  !< norm(r)
    !real(r8), dimension(:,:), allocatable :: ir !< 1/norm(r)
    ! Sqrt of integration weight
    real(r8), dimension(:), allocatable :: wt
    ! intermediate basis function values
    real(r8), dimension(:,:), allocatable :: f   !< radial basis function values
    real(r8), dimension(:,:), allocatable :: Ylm !< spherical harmonic values
    real(r8), dimension(:,:), allocatable :: basis_times_wt !< basis function values times weight
    real(r8), dimension(:,:), allocatable :: basis_times_wt_tr !< basis function values times weight
end type

contains

!> Calulcate the density
subroutine rl_calculate_overlapmatrix(grid,rmtx,KS,p,basis,ec,mw,mem,verbosity,tmr)
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> distributed realspace matrices. After tis routine, te overlap will be updated.
    class(rl_realspace_matrix), allocatable, intent(inout) :: rmtx
    !> Kohn-Sham solution thingy handle.
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> timer
    type(rl_timer_overlap), intent(inout) :: tmr

    ! helper tings needed for te integration
    type(rl_mom_real), dimension(:), allocatable :: ovl_buf
    type(rl_pair_locator) :: ploc
    type(integrationhelper) :: ih
    real(r8), dimension(8) :: mellantid
    real(r8) :: timer,t0,t1
    !init: block
    integer :: a1,nb1,nb2,ipair
    !getovl: block
    integer :: ibatch
    !store: block
    !integer :: ipair,i
    integer :: i

    ! Do some stuff. Maybe allocate temporary storage and things?
    !init: block

        ! Start timers
        timer=mpi_wtime()
        t0=timer
        t1=timer
        mellantid=0.0_r8

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CALCULATING OVERLAP'
        endif

        ! First we generate the pair locator for lookup.
        call ploc%generate(rmtx,KS,ec,irreducible=.true.,mw=mw,mem=mem)
        t1=mpi_wtime(); mellantid(1)=t1-t0; t0=t1

        ! Initialize the integration helper
        call init_integration_helper(ih,p,grid,basis,rmtx)

        ! Create space for the overlap buffer. This is temporarey storage were
        ! I accumulate te overlap matrix, and towards the end I store it in
        ! the right place.
        allocate(ovl_buf(ih%n_pair_global))
        do ipair=1,ih%n_pair_global
            nb1=ih%nbasis_atom1(ipair)
            nb2=ih%nbasis_atom2(ipair)
            call mem%allocate(ovl_buf(ipair)%m,[nb1,nb2],persistent=.false.,scalable=.false.)
            ovl_buf(ipair)%m=0.0_r8
        enddo
        t1=mpi_wtime(); mellantid(2)=t1-t0; t0=t1

        if ( verbosity .gt. 0 ) then
            write(*,*) '... initialized integration helper'
        endif
    !end block init

    ! Start integrating. Instead of making this block very long, I chopped it
    ! into bite-sized pieces that each is quite easy to understand.
    !getovl: block

        do ibatch=1,grid%n_full_batch
            ! This resets all the temporary values to be relevant for this batch
            call reset_integration_helper(ih,ec,basis,grid,grid%full_batch(ibatch),mem)
            t1=mpi_wtime(); mellantid(3)=mellantid(3)+t1-t0; t0=t1

            ! First thing we need to do is tabulate the coordinates of all the atoms
            ! whose basis functions can touch the points in this batch:
            call tabulate_coordinates(ih,ec,grid%full_batch(ibatch))

            ! Also a good idea to refresh the list of which basis
            ! functions we have to evaluate
            call update_active_basis_functions(ih,ec,p,basis)
            t1=mpi_wtime(); mellantid(4)=mellantid(4)+t1-t0; t0=t1

            ! With the list of coordinates and the list of basis functions,
            ! we evaluate the basis functions
            call tabulate_basis_functions(ih,basis)
            t1=mpi_wtime(); mellantid(5)=mellantid(5)+t1-t0; t0=t1

            ! And finally, we accumulate the contributions to the overlap matrix
            call buildbigmatrix(ih,ploc,ovl_buf)
            t1=mpi_wtime(); mellantid(6)=mellantid(6)+t1-t0; t0=t1
        enddo
    !end block getovl

    ! Store the overlap matrix in the right place
    !store: block

        ! First I have to communicate the overlap.
        !@todo make faster.
        do ipair=1,size(ovl_buf)
            call mw%allreduce('sum',ovl_buf(ipair)%m)
        enddo

        ! Then store it in the right place
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! Simple enough, just store it
            do ipair=1,m%n_irr_pair
                m%irr_pair(ipair)%overlap = ovl_buf(ipair)%m
            enddo
        type is(rl_realspace_matrix_mediumdistributed)
            do i=1,m%n_irr_pair
                ipair=m%irr_offset+i
                m%irr_pair(i)%overlap = ovl_buf(ipair)%m
            enddo
        type is(rl_realspace_matrix_fulldistributed)
            do i=1,m%n_irr_pair
                ipair=m%irr_offset+i
                m%irr_pair(i)%overlap = ovl_buf(ipair)%m
            enddo
        end select

        ! And destroy the temporary storage
        do ipair=1,size(ovl_buf)
            call mem%deallocate(ovl_buf(ipair)%m,persistent=.false.,scalable=.false.)
        enddo
        deallocate(ovl_buf)

        if ( verbosity .gt. 0 ) then
            write(*,*) '... communicated and stored overlap'
        endif
    !end block store

    ! Time the communications step
    t1=mpi_wtime(); mellantid(7)=mellantid(7)+t1-t0; t0=t1

    ! And finally, sort out the timings
    mellantid(8)=mpi_wtime()-timer                ! total timer
    mellantid(7)=mellantid(8)-sum(mellantid(1:6)) ! time idle
    tmr%init        = tmr%init +        mellantid(2)
    tmr%idle        = tmr%idle +        mellantid(7)
    tmr%total       = tmr%total +       mellantid(8)
    tmr%pairlocator = tmr%pairlocator + mellantid(1)
    tmr%reset       = tmr%reset +       mellantid(3)
    tmr%tab_coord   = tmr%tab_coord +   mellantid(4)
    tmr%tab_basis   = tmr%tab_basis +   mellantid(5)
    tmr%matrixop    = tmr%matrixop +    mellantid(6)

    if ( verbosity .gt. 0 ) then
        write(*,*) '... done integration overlap'
    endif

    ! Check for stupidity
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)
end subroutine

!> tabulate the coordinates of all atoms with respect to all points
subroutine tabulate_coordinates(ih,ec,batch)
    !> subset of basis functions and atoms
    type(integrationhelper), intent(inout) :: ih
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> batch information
    class(rl_integration_batch), intent(in) :: batch

    real(r8), parameter :: radtol=1E-10_r8
    real(r8), dimension(3) :: v,w
    real(r8) :: f0,f1
    integer :: ia,ie,ip

    ! calculate coordinates
    do ia=1,ih%n_active_atom
        ie=ih%active_atom(ia)
        w=ec%cartesian_coordinate(:,ie)
        do ip=1,batch%n_point
            v=batch%folded_coordinate(:,ip)-w
            f0=norm2(v)
            if ( f0 .gt. radtol ) then
                f1=1.0_r8/f0
                v=v*f1
                ih%x(ip,ia)=v(1)
                ih%y(ip,ia)=v(2)
                ih%z(ip,ia)=v(3)
                ih%r(ip,ia)=f0
                !ih%ir(j,i)=f1
            else
                ih%x(ip,ia)=0.0_r8
                ih%y(ip,ia)=0.0_r8
                ih%z(ip,ia)=0.0_r8
                ih%r(ip,ia)=f0
                !ih%ir(j,i)=0.0_r8
            endif
        enddo
    enddo
end subroutine

!> convert the list of active atoms to a list of active basis functions
subroutine update_active_basis_functions(ih,ec,p,basis)
    !> subset of basis functions and atoms
    type(integrationhelper), intent(inout) :: ih
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
    ih%basis_ind=0
    l=0
    do i=1,ih%n_active_atom
        ie=ih%active_atom(i)
        iu=ec%index_unit_cell(ie)
        is=p%species(iu)
        nb=basis%species(is)%n_basis
        ih%active_atom_offset(i)=l
        ih%active_atom_nbasis(i)=nb
        do j=1,nb
            l=l+1
            ih%basis_ind(:,l)=[basis%species(is)%basis_fn(j),basis%species(is)%basis_ang(j)]
        enddo
    enddo
    ih%n_nonzero_basis=l
    ! Some space for the basis function values
    allocate(ih%basis_times_wt(ih%n_nonzero_basis,ih%n_integration_point))
    allocate(ih%basis_times_wt_tr(ih%n_integration_point,ih%n_nonzero_basis))
    ih%basis_times_wt=0.0_r8
    ih%basis_times_wt_tr=0.0_r8
end subroutine

!> build the huge matrices and start accumulating the overlap. Hmmm.
subroutine buildbigmatrix(ih,ploc,ovl_buf)
    !> subset of basis functions and atoms
    type(integrationhelper), intent(inout) :: ih
    !> pair locator
    type(rl_pair_locator), intent(in) :: ploc
    !> actual overlap matrix
    type(rl_mom_real), dimension(:), intent(inout) :: ovl_buf

    real(r8), dimension(:,:), allocatable :: mA
    integer :: a1,s1,np,ipair
    integer :: nbi,nbj
    integer :: ia,ja,ie,je
    integer :: i1,i2,j1,j2

    ! Here we can switch depending on the symmetry reduction, with some
    ! sensible crossover maybe.
    select case(ih%update_mode)
    case(1)
        ! Not very much symmetry, do normal update
        allocate(mA(ih%n_nonzero_basis,ih%n_nonzero_basis))
        mA=0.0_r8
        call dsyrk('U','N', ih%n_nonzero_basis, ih%n_integration_point, 1.0_r8, ih%basis_times_wt, ih%n_nonzero_basis, 0.0_r8, mA, ih%n_nonzero_basis)
        do i1=1,ih%n_nonzero_basis
            do i2=i1+1,ih%n_nonzero_basis
                mA(i2,i1)=mA(i1,i2)
            enddo
        enddo

        do ia=1,ih%n_active_atom
        do ja=1,ih%n_active_atom
            ie=ih%active_atom(ia)
            je=ih%active_atom(ja)
            ipair=ploc%locate(ie,je)
            if ( ipair .lt. 1 ) cycle
            i1=ih%active_atom_offset(ia)+1
            i2=ih%active_atom_offset(ia)+ih%active_atom_nbasis(ia)
            j1=ih%active_atom_offset(ja)+1
            j2=ih%active_atom_offset(ja)+ih%active_atom_nbasis(ja)
            ovl_buf(ipair)%m=ovl_buf(ipair)%m+mA(i1:i2,j1:j2)
        enddo
        enddo
        deallocate(mA)
    case(2)
        ! This is the incremental version of the update. Slower on the matrix
        ! operations, but way fewer matrix operations, I hope.
        do ia=1,ih%n_active_atom
            ie=ih%active_atom(ia)
            do ja=1,ih%n_active_atom
                je=ih%active_atom(ja)
                ipair=ploc%locate(ie,je)
                if ( ipair .gt. 0 ) then
                    nbi=ih%active_atom_nbasis(ia)
                    nbj=ih%active_atom_nbasis(ja)
                    i1=ih%active_atom_offset(ia)+1
                    i2=ih%active_atom_offset(ia)+nbi
                    j1=ih%active_atom_offset(ja)+1
                    j2=ih%active_atom_offset(ja)+nbj
                    ! Accumulate overlap
                    call dgemm('T','N', nbi, nbj, ih%n_integration_point, 1.0_r8, ih%basis_times_wt_tr(:,i1:i2), ih%n_integration_point, ih%basis_times_wt_tr(:,j1:j2), ih%n_integration_point, 1.0_r8, ovl_buf(ipair)%m, nbi)
                endif
            enddo
        enddo
    end select
end subroutine

!> tabulate the values of the wave functions over the batch
subroutine tabulate_basis_functions(ih,basis)
    !> subset of basis functions and atoms
    type(integrationhelper), intent(inout) :: ih
    !> basis set information
    type(rl_lcao_basis_set), intent(in) :: basis

    real(r8), dimension(:,:), allocatable :: buf
    integer :: i,j,i1,i2,ir,ia,np

    allocate(buf(ih%n_integration_point,ih%n_nonzero_basis))
    buf=0.0_r8

    np=ih%n_integration_point
    ! No derivatives, just function values
    do i=1,ih%n_active_atom
        ! Evaluate radial functions
        do ir=1,basis%n_radial
            call basis%radial(ir)%val( ih%r(1:np,i),ih%f(1:np,ir) )
        enddo
        ! Evaluate spherical harmonics
        do ia=1,basis%n_angular
            call basis%Ylm( ih%x(1:np,i),ih%y(1:np,i),ih%z(1:np,i),ia,ih%Ylm(1:np,ia) )
        enddo
        ! Combine basis functions to the full ones
        do j=ih%active_atom_offset(i)+1,ih%active_atom_offset(i)+ih%active_atom_nbasis(i)
            ir = ih%basis_ind(1,j) ! radial index of basis fn j
            ia = ih%basis_ind(2,j) ! angular index of basis fn j
            ! where to store values in flat array
            ! i1=(j-1)*np+1
            ! i2=j*np
            ! ih%sqrtwt_times_basis(i1:i2)=ih%f(1:np,ir) * ih%Ylm(1:np,ia)
            ih%basis_times_wt_tr(:,j)=ih%f(1:np,ir) * ih%Ylm(1:np,ia)
        enddo
    enddo
    ! Transpose it, and multiply in sqrt of integration weight
    ih%basis_times_wt=transpose(ih%basis_times_wt_tr)
    do i=1,np
        ih%basis_times_wt(:,i)=ih%basis_times_wt(:,i)*ih%wt(i)
    enddo
    ih%basis_times_wt_tr=transpose(ih%basis_times_wt)

    deallocate(buf)
end subroutine

!> reset the integration helper
subroutine reset_integration_helper(ih,ec,basis,grid,batch,mem)
    !> subset of basis functions and atoms
    class(integrationhelper), intent(inout) :: ih
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

    integer :: i

    ! Zero all basis functions
    ih%n_nonzero_basis    = 0
    ih%active_atom_offset = 0
    ih%active_atom_nbasis = 0
    ih%basis_ind          = 0

    ! Make a note of the number of integration points
    ih%n_integration_point=batch%n_point

    ! set the active atoms to the ones determined by the batch?
    ih%n_active_atom=batch%n_relevant_ext_atom
    ih%active_atom = 0
    ih%active_atom(1:ih%n_active_atom)=batch%relevant_ext_atom(1:ih%n_active_atom)

    ! set coordinates to noting
    ih%x  = 0.0_r8
    ih%y  = 0.0_r8
    ih%z  = 0.0_r8
    ih%r  = 0.0_r8
    ! I need the square root of the integration weights
    ih%wt = 0.0_r8
    do i=1,batch%n_point
        ih%wt(i)=sqrt(batch%integration_weight(i)*batch%partition_function(i))
    enddo

    ! space for basis function values
    ih%f       = 0.0_r8
    ih%Ylm     = 0.0_r8

    if ( allocated(ih%basis_times_wt) ) deallocate(ih%basis_times_wt)
    if ( allocated(ih%basis_times_wt_tr) ) deallocate(ih%basis_times_wt_tr)
end subroutine

!> make some space in the integration helper
subroutine init_integration_helper(ih,p,grid,basis,rmtx)
    !> integration helper
    type(integrationhelper), intent(out) :: ih
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> integration grid
    type(rl_integration_grid), intent(in) :: grid
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> realspace matrix
    class(rl_realspace_matrix), intent(in) :: rmtx

    real(r8) :: f0
    integer :: a1,ipair

    allocate(ih%active_atom(        grid%max_n_atom_per_batch) )
    allocate(ih%active_atom_offset( grid%max_n_atom_per_batch) )
    allocate(ih%active_atom_nbasis( grid%max_n_atom_per_batch) )
    allocate(ih%basis_ind(        2,grid%max_n_basis_per_batch))

    ! Space for coordinates
    allocate(ih%x (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(ih%y (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(ih%z (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(ih%r (grid%max_n_point_per_batch,grid%max_n_atom_per_batch))
    allocate(ih%wt(grid%max_n_point_per_batch))
    ih%x =0.0_r8
    ih%y =0.0_r8
    ih%z =0.0_r8
    ih%r =0.0_r8
    ih%wt=0.0_r8

    ! Space for intermediate basis function values
    allocate(ih%f   ( grid%max_n_point_per_batch,basis%n_radial  ))
    allocate(ih%Ylm ( grid%max_n_point_per_batch,basis%n_angular ))
    ih%f  =0.0_r8
    ih%Ylm=0.0_r8

    ! Start counting pair things, will need this information eventually.
    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        ! make a note of the global number of pairs
        ih%n_pair_global=m%n_irr_pair
        ! space for dummy counting array
        allocate(ih%nbasis_atom1(ih%n_pair_global))
        allocate(ih%nbasis_atom2(ih%n_pair_global))
        ih%nbasis_atom1=0
        ih%nbasis_atom2=0
        do ipair=1,rmtx%n_irr_pair
            ih%nbasis_atom1(ipair)=rmtx%irr_pair(ipair)%nb1
            ih%nbasis_atom2(ipair)=rmtx%irr_pair(ipair)%nb2
        enddo

        ! Compare the number of irreducible pairs to the total, and make
        ! some decision based on that:
        f0=m%n_irr_pair/real(m%n_full_pair,r8)
        if ( f0 .gt. 0.4_r8 ) then
            ih%update_mode=1
        else
            ih%update_mode=2
        endif
    class default
        write(*,*) 'OLLE MUST FIX ALL DISTRIBUTIONS IN OVERLAP INTEGRATION HELPER INIT'
        stop
    end select

end subroutine

end module
