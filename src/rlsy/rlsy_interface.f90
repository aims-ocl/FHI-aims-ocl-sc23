module rlsy_interface
!!
!! This module defines the handle used to deal with symmetry. And maybe a few more things.
!!
use rlsy_constants, only: r8,rl_iou,rl_exitcode_symmetry,rl_exitcode_memory,rl_sqtol,rl_hugeint
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: tochar,rl_sqnorm,rl_clean_fractional_coordinates
use rlsy_mpi_helper, only: rl_mpi_helper,mpi_wtime,rl_stop_gracefully
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_extended_cluster, only: rl_extended_cluster
use rlsy_integration_grid, only: rl_integration_grid
use rlsy_free_atom, only: rl_free_atom
use rlsy_electron_density, only: rl_electron_density
use rlsy_realspace_matrix, only: rl_realspace_matrix,rl_setup_realspace_storage
use rlsy_kpointmesh, only: rl_kpoint_mesh
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,rl_setup_kspace_eigenproblem
use rlsy_hartree_potential, only: rl_multipole_expansion
use rlsy_timer, only: rl_timer

implicit none

private
public :: rlsy_handle
public :: rlsy_h
public :: generate_rlsy_handle

! Symmetry handle, contains all information regarding the symmetry of the system.
! Definitely a work in progress, will get trimmed down to the bare minimum, eventually.
! The reason there is so much repeated information from AIMS is that I don't know what
! is actually needed, or how AIMS actually works, so I need a bunch of things for
! debugging/development.
type rlsy_handle
    !> how much to talk?
    integer :: verbosity=-rl_hugeint
    !> Handle to write to stdout. Shh, don't tell anyone.
    integer :: stdout=-rl_hugeint
    !> Are we running as a library? Always good to know.
    logical :: run_as_library=.false.
    !> contains that holds everything related to the crystal structure
    type(rl_crystalstructure) :: structure
    !> spacegroup information
    type(rl_spacegroup) :: spacegroup
    !> basis set handle (retire when I know exactly what is needed)
    type(rl_lcao_basis_set) :: basis
    !> free atom handle (retire obviously)
    type(rl_free_atom) :: free
    !> extended cluster of atoms
    type(rl_extended_cluster) :: ec
    !> integration grid, divided into batches in various ways
    type(rl_integration_grid) :: grid
    !> electron density, defined on chunks of the irreducible points
    type(rl_electron_density) :: density
    !> electron multipole expansion
    type(rl_multipole_expansion) :: multipole
    !> handle for realspace matrices
    class(rl_realspace_matrix), allocatable :: rmtx
    !> kpoint mesh
    type(rl_kpoint_mesh) :: kmesh
    !> handle for eigenvalue problem. Should rename to something less stupid.
    class(rl_kspace_eigenproblem), allocatable :: KS
    !> MPI helper
    type(rl_mpi_helper) :: mw
    !> Memory tracker
    type(rl_memtracker) :: mem
    !> timing information
    type(rl_timer) :: timer
    contains
        !> create the handle
        procedure :: generate=>generate_rlsy_handle
        !> check if we need to rebuild the handle
        procedure :: check_handle_rebuild
end type rlsy_handle

! I suppose I have no choice but to expose a global variable, even if it hurts.
type(rlsy_handle), save :: rlsy_h

contains

!> check if we need to rebuild the handle.
subroutine check_handle_rebuild(rh,frac_coords,lattice_vector,need_reinit)
    !> symmetry handle
    class(rlsy_handle), intent(inout) :: rh
    !> geometry from AIMS
    real(r8), dimension(:,:), intent(in) :: frac_coords
    real(r8), dimension(3,3), intent(in) :: lattice_vector
    !> do we need to reinitialize?
    logical, intent(out) :: need_reinit

    integer :: i
    real(r8), dimension(3) :: v0
    real(r8) :: f0

    need_reinit=.false.
    ! Check if the lattice vectors have changed?
    f0=sum(abs(lattice_vector-rh%structure%latticevectors))
    if ( f0 .gt. 1E-14_r8 ) then
        need_reinit=.true.
        return
    endif

    ! Check if coordinates of atoms have changed
    do i=1,size(frac_coords,2)
        v0=frac_coords(:,i)-rh%structure%fractional_coordinate(:,i)
        v0=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
        if ( sum(abs(v0)) .gt. 1E-14_r8 ) then
            need_reinit=.true.
            return
        endif
    enddo
end subroutine

!> Create all the information needed to play around with symmetry. This is likely far too much, and should get pruned down to the bare necessities once things converge.
subroutine generate_rlsy_handle(rh,&
    species,species_z,frac_coords,lattice_vector,&
    basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
    r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl,outer_partition_radius,&
    n_k_points_xyz,k_point_list,k_points_offset,&
    atom_radius,free_rho_spl,free_drho_dr_spl,free_pot_es_spl,free_pot_es_at_zero, &
    n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,scale_radial,l_hartree,&
    n_electrons,n_states,n_spin,&
    elsi_solver,elsi_out_level,occupation_type,solver_method,&
    elsi_elpa_n_single,occupation_width,occupation_acc,basis_threshold,&
    verbosity,use_unit,mpi_world_communicator,nosym)
    !> container that holds all symmetry information
    class(rlsy_handle), intent(out) :: rh
    !> structure information from AIMS
    integer, dimension(:), intent(in) :: species
    real(r8), dimension(:), intent(in) :: species_z
    real(r8), dimension(:,:), intent(in) :: frac_coords
    real(r8), dimension(3,3), intent(in) :: lattice_vector
    !> basis set information from AIMS
    integer, dimension(:), intent(in) :: basis_m
    integer, dimension(:), intent(in) :: basis_l
    integer, dimension(:), intent(in) :: basis_fn
    integer, dimension(:), intent(in) :: basis_atom
    real(r8), dimension(:), intent(in) :: outer_radius
    real(r8), dimension(:), intent(in) :: r_grid_min
    real(r8), dimension(:), intent(in) :: r_grid_inc
    real(r8), dimension(:,:,:), intent(in) :: basis_wave_spl
    real(r8), dimension(:,:,:), intent(in) :: basis_deriv_spl
    real(r8), dimension(:), intent(in) :: outer_partition_radius
    !> k-mesh information from AIMS
    integer, dimension(3), intent(in) :: n_k_points_xyz
    real(r8), dimension(:,:), intent(in) :: k_point_list
    real(r8), dimension(3), intent(in) :: k_points_offset
    !> grid information from AIMS
    integer, dimension(:), intent(in) :: n_radial
    integer, dimension(:,:), intent(in) :: n_angular
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:,:), intent(in) :: w_radial
    real(r8), dimension(:,:,:), intent(in) :: w_angular
    real(r8), dimension(:), intent(in) :: scale_radial
    integer, dimension(:), intent(in) :: l_hartree
    ! integer, intent(in) :: use_batch_permutation
    ! integer, intent(in) :: n_my_batches_work
    ! type(batch_of_points), dimension(:), pointer, intent(in) :: batches_work
    !> free atom
    real(r8), dimension(:), intent(in) :: atom_radius
    real(r8), dimension(:,:,:), intent(in) :: free_rho_spl
    real(r8), dimension(:,:,:), intent(in) :: free_drho_dr_spl
    real(r8), dimension(:,:,:), intent(in) :: free_pot_es_spl
    real(r8), dimension(:), intent(in) :: free_pot_es_at_zero
    !> KS solution related things
    real(r8), intent(in) :: n_electrons
    integer, intent(in) :: n_states
    integer, intent(in) :: n_spin
    !> ELSI related things
    integer, intent(in) :: elsi_solver
    integer, intent(in) :: elsi_out_level
    integer, intent(in) :: occupation_type
    integer, intent(in) :: solver_method
    integer, intent(in) :: elsi_elpa_n_single
    real(r8), intent(in) :: occupation_width
    real(r8), intent(in) :: occupation_acc
    real(r8), intent(in) :: basis_threshold
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> which unit to talk to
    integer, intent(in) :: use_unit
    !> mpi world communicator
    integer, intent(in) :: mpi_world_communicator
    !> should I switch off all symmetries except for identity? Only for debugging.
    logical, intent(in) :: nosym

    real(r8) :: timer,t0,t1

    !init: block
        integer, dimension(:), allocatable :: atomic_numbers
        integer :: i,na
    !sortallatoms: block
        real(r8), dimension(3) :: v0
        real(r8) :: cutoff
        !integer :: a1,i
        integer :: a1
    !fin: block
        character(len=500) :: opf
        real(r8) :: tomb,totmem,f0

    ! I know the argument list is quite ridiculously long, but I prefer it that way
    ! as opposed to blindly using modules and all side effects that can come with that.
    ! This way it is clean, and clear that nothing is modified, and also clear exactly
    ! what is needed to sort out all the symmetries.

    ! Set some basic things
    !init: block

        ! Start local clock
        timer=mpi_wtime()
        t0=timer

        ! Initialize all timers to zero
        call rlsy_h%timer%init()

        ! Make sure we only write to the unit that AIMS also writes to
        rl_iou=use_unit
        ! Now for an ugly kludge: when I, Olle, personally run AIMS as a library it is neat
        ! to be able to monitor progress, and print som minimal messages during runs. But I
        ! don't want to mess up the usual output file and break thing. My somewhat odd solution
        ! is to check if aims is writing to unit 6, since that is hard-coded in aims.f90. If we
        ! are not writing to that unit, my guesstimate is that we are running as a library.
        ! Someone invented precompilers to resolve these kinds of issues, but apparently that
        ! is not an option, hence I resort to this ugly kludge.
        if ( use_unit .ne. 6 ) then
            rh%stdout=6
            rh%run_as_library=.true.
        else
            ! This is a normal run.
            rh%stdout=-rl_hugeint
            rh%run_as_library=.false.
        endif

        ! Initialize the memory tracker for dynamic memory
        call rh%mem%init()

        ! First, I need handy access to the MPI communicator, so stash that
        call rh%mw%init(communicator=mpi_world_communicator)
        ! Make all but one rank quiet
        if ( rh%mw%talk .eqv. .false. ) then
            rh%verbosity=-100
        else
            rh%verbosity=verbosity
        endif

        if ( rh%verbosity .gt. 0 ) then
            write(rl_iou,*) ''
            write(rl_iou,*) 'SORTING OUT SYMMETRIES'
        endif

        ! Generate my structure object thingy with information from
        ! AIMS. What I collect and use is
        !   lattice_vector
        !   species
        !   species_z
        !   frac_coords
        ! This might get extended with magnetic things later to get the magnetic
        ! space group as well. Anyway, I take this information and generate a crystal
        ! structure object.
        na=size(species)
        call rh%mem%allocate(atomic_numbers,na,persistent=.false.,scalable=.false.)
        do i=1,na
            atomic_numbers(i)=int(anint(species_z(species(i))))
        enddo
        ! And pack the structure into my handle
        call rh%structure%generate(lattice_vector,frac_coords,atomic_numbers,species,rh%mw,rh%mem)
        call rh%mem%deallocate(atomic_numbers,persistent=.false.,scalable=.false.)

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got structure ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Generate the spacegroup and the irreducible representations
        call rh%spacegroup%generate(&
            latticevectors=rh%structure%latticevectors,&
            timereversal=.true.,&
            fractional_coordinate=rh%structure%fractional_coordinate,&
            species=rh%structure%species,symmorphic=.false.,tol=1E-5_r8,&
            verbosity=-1,mw=rh%mw,nosym=nosym,mem=rh%mem)
        !call rh%spacegroup%get_character_table(verbosity=-1,mw=rh%mw)

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got spacegroup ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Now that the structure is settled on, I have to make sure I understand the
        ! basis set. I take the full information about the basis functions and
        ! rearrange it in a more compact format. This will likely be removed in the
        ! future, but I need it for understanding.
        call rh%basis%generate(rh%structure,rh%mw,rh%mem,&
             basis_m,basis_l,basis_fn,basis_atom,outer_radius,&
             r_grid_min,r_grid_inc,basis_wave_spl,basis_deriv_spl)

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... sorted basis functions ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Likewise, I have to understand how the free atom quantities work, so I will
        ! package those as well. Will also be removed once my understanding is complete.
        ! So likely never, but you never know.
        call rh%free%generate(rh%structure,rh%mw,rh%mem,&
            free_rho_spl,free_drho_dr_spl,free_pot_es_spl,&
            free_pot_es_at_zero,r_grid_min,r_grid_inc,atom_radius)

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... sorted free atom quantities ('//tochar(t1-t0)//'s)'
            t0=t1
        endif
    !end block init

    if ( rh%mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,rh%mw%comm)

    ! Create a huge list of all the atoms in an extended scheme, always handy thing to have.
    ! There is a sensible reason why I don't use the list from AIMS,
    !sortallatoms: block

        ! This is a safe upper bound on the cutoff
        cutoff=maxval(outer_radius)*2
        ! This generates the full extended cluster in a way that preserves symmetry
        call rh%ec%generate(rh%structure,cutoff,rh%mw,rh%verbosity,nosym,rh%mem)

        ! Generate irrep projections. Later problem.
        !call as%irrep%generate(p,as%ec,as%basis,mw,verbosity)

        ! A small sanity test to make sure the cluster makes sense. Maybe remove.
        do i=1,rh%ec%n_extended_atom
            a1=rh%ec%index_unit_cell(i)
            v0=rh%ec%cartesian_coordinate(:,i)-rh%structure%cartesian_coordinate(:,a1)
            v0=matmul(rh%structure%inv_latticevectors,v0)
            v0=v0-anint(v0)
            if ( rl_sqnorm(v0) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['Clearly I do not understand symmetry'],rl_exitcode_memory,rh%mw%comm)
            endif
        enddo

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... done with extended cluster ('//tochar(t1-t0)//'s)'
            t0=t1
        endif
    !end block sortallatoms

    ! And start worrying about the grids
    !setgrids: block
        !integer :: l_hartree
        !@TODO I have not pruned the cluster now, for that I need to do the irrep thing first.

        ! Actually not sure about that pruning, maybe later problem when I project onto irreps.
        call rh%grid%generate(rh%structure,rh%spacegroup,rh%ec,rh%basis,rh%free,rh%mw,rh%mem,rh%verbosity,&
             n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,outer_partition_radius,scale_radial)

        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... done partitioning grid ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! dummy max l_hartree

        ! Initialize storage on the grid.
        call rh%density%generate(rh%grid,rh%structure,rh%free,rh%ec,n_spin,n_electrons,rh%mw)

        if ( rh%verbosity .gt. 0 ) then
           t1=mpi_wtime()
           write(rl_iou,*) '... done setting up density ('//tochar(t1-t0)//'s)'
           t0=t1
        endif

        ! ! Initialize multipole expansion.
        ! call rh%multipole%generate(l_hartree,rh%structure,rh%spacegroup,rh%grid,rh%mw,rh%mem,rh%verbosity)
        !
        ! if ( rh%verbosity .gt. 0 ) then
        !    t1=mpi_wtime()
        !    write(rl_iou,*) '... done setting up multipole expansion ('//tochar(t1-t0)//'s)'
        !    t0=t1
        ! endif

    !end block setgrids

    ! Sort the k-point mesh
    !kspace: block

        ! Get the irreducible k-mesh
        call rh%kmesh%generate(rh%structure,rh%spacegroup,rh%mw,rh%mem,rh%verbosity,n_k_points_xyz,k_point_list,k_points_offset)

!        ! Setup ELSI and the rest needed for the solver, eventually.
!        call rl_setup_kspace_eigenproblem(rh%KS,rh%structure,rh%basis,rh%kmesh,n_states,n_spin,n_electrons,rh%mw,rh%mem,rh%verbosity,&
!            elsi_solver,elsi_out_level,occupation_type,solver_method,elsi_elpa_n_single,&
!            occupation_width,occupation_acc,basis_threshold)
!
!        if ( rh%verbosity .gt. 0 ) then
!            t1=mpi_wtime()
!            write(rl_iou,*) '... generated irreducible k-mesh and storage ('//tochar(t1-t0)//'s)'
!            t0=t1
!        endif
    !end block kspace

    ! Matrix storage goes here.
    !setmatrixstorage: block
        ! This creates my kind of distributed realspace storage.
!        call rl_setup_realspace_storage(rh%rmtx,rh%structure,rh%basis,rh%ec,rh%spacegroup,rh%KS,rh%mw,rh%verbosity,rh%mem)
!
!        if ( rh%verbosity .gt. 0 ) then
!            t1=mpi_wtime()
!            write(rl_iou,*) '... initialized storage ('//tochar(t1-t0)//'s)'
!            t0=t1
!        endif

        ! Quick test to see that Fourier transforms work.
        !call rl_dummy_fourier_transform_to_test_that_things_work(rh%rmtx,rh%KS,rh%structure,rh%spacegroup,rh%basis,rh%kmesh,rh%mw,rh%mem,rh%verbosity)
    !end block setmatrixstorage

    ! Then some report
    !fin: block


        opf='(1X,A20,1X,F12.4,1X,"MiB")'
        tomb=1.0_r8/1024.0_r8**2
        totmem=0.0_r8
        if ( rh%verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) ''
            write(rl_iou,*) 'Memory allocated for symmetry things:'
            f0=rh%structure%size_in_mem()*tomb
            totmem=totmem+f0
            write(rl_iou,opf) '       structure:',f0
            f0=rh%basis%size_in_mem()*tomb
            totmem=totmem+f0
            write(rl_iou,opf) '           basis:',f0
            f0=rh%spacegroup%size_in_mem()*tomb
            totmem=totmem+f0
            write(rl_iou,opf) '      spacegroup:',f0
            f0=rh%ec%size_in_mem()*tomb
            totmem=totmem+f0
            write(rl_iou,opf) '         cluster:',f0
            !f0=rh%rmtx%size_in_mem()*tomb
            !totmem=totmem+f0
            !write(rl_iou,opf) '           H/S/D:',f0
            write(rl_iou,opf) '           total:',totmem
            write(rl_iou,opf) ''
            write(rl_iou,*) 'Initialized symmetry handler ('//tochar(t1-timer)//'s)'
        endif
    !end block fin

    ! And be sure I do temporary memory correctly
    if ( rh%mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,rh%mw%comm)
    if ( rh%mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,rh%mw%comm)


end subroutine

end module rlsy_interface
