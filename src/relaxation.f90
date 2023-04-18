!****h* FHI-aims/relaxation
!  NAME
!    relaxation
!  SYNOPSIS

module relaxation

  !  PURPOSE
  !   This module takes care of:
  !   * allocation/deallocation of relaxation variables
  !   * relaxation restarts & saving of all necessary information
  !   * relaxation with a number of different algorithms
  !
  !  CONVENTION FOR VARIABLES
  !   o mobile_atom_coords:   Index picture of thoses components of atomic coordinates which are allowed to relax
  !   o mobile_lv_components: Index picture thoses components of lattice vectors which are allowed to relax
  !   o mobile_coords:        Index picture of all the coordinates (from atoms AND lattice) which relax
  !   o n_*                   Number of elements in the corresponding arrays
  !
  !  USES
  use localorb_io
  use runtime_choices
  use dimensions
  use constants
  use geometry
  use scalapack_wrapper
  use mpi_tasks
  use generate_aims_uuid, only: write_aims_uuid
  use aims_memory_tracking, only: aims_allocate, aims_deallocate
  use distributed_hessian
  use symmetry_constrained_relaxation
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !    to the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE
  implicit none

  ! to keep all the information of a single relaxation step for future use
  real*8, dimension(:,:), allocatable :: stored_coords
  real*8 :: stored_energy
  real*8, dimension(:,:), allocatable :: stored_forces
  ! VA: similar for the lattice vectors
  real*8, dimension(:,:), allocatable :: stored_lv
  real*8, dimension(:,:), allocatable :: stored_forces_lv
  ! FlK: The step in coords should be explicitly stored
  real*8, dimension(:), allocatable   :: stored_DX
  ! FK: Need to store shift of atoms for every geometry
  real*8, dimension(:,:), allocatable :: stored_periodic_unit_cell_trans

  ! ... that was not enough: might need to remember two relaxation steps, actually
  real*8, dimension(:,:), allocatable :: stored_coords_2
  real*8 :: stored_energy_2
  real*8, dimension(:,:), allocatable :: stored_forces_2
  ! VA: similar for the lattice vectors
  real*8, dimension(:,:), allocatable :: stored_lv_2
  real*8, dimension(:,:), allocatable :: stored_forces_lv_2
  ! FK: again
  real*8, dimension(:,:), allocatable :: stored_periodic_unit_cell_trans_2

  real*8, dimension(:,:,:,:), allocatable :: stored_KS_eigenvector
  complex*16, dimension(:,:,:,:), allocatable :: stored_KS_eigenvector_complex
  real*8, dimension(:,:,:), allocatable :: stored_occ_numbers

  ! For handling of atoms whose relaxation is constrained
  logical, dimension(:), allocatable :: constrain_relaxation
  logical, dimension(:,:), allocatable :: constrain_coords

  ! VA: For handling of lattice vectors whose relaxation is constrained
  logical, dimension(:), allocatable :: constrain_lv_relaxation
  logical, dimension(:,:), allocatable :: constrain_lv_components

  ! index lists of atoms for relaxation - fixed vs. mobile
  integer :: n_mobile_atoms
  integer, dimension(:), allocatable :: mobile_atom
  integer, dimension(:), allocatable :: mobile_atom_coord
  integer :: n_mobile_atom_coords
  integer :: n_mobile_coords

  ! VA: index lists for unit cell relaxation - fixed vs. mobile
  integer :: n_mobile_lv
  integer, dimension(:), allocatable :: mobile_lv
  integer :: n_mobile_lv_components
  integer, dimension(:), allocatable :: mobile_lv_component

  ! this is for descent-directed MD
  real*8 :: relaxation_step
  real*8 :: step_decrease
  real*8 :: step_minimum
  real*8, dimension(:,:), allocatable :: velocity
  logical :: step_found
  logical :: initial_geo
  logical :: switched_to_bfgs

  ! remember minimum energy ever encountered across the whole
  ! relaxation trajectory (safeguard against protracted uphill relaxation in small steps)
  real*8, private :: minimum_energy_ever
  logical, private :: minimum_energy_ever_not_initialized

  real*8 :: step_maximum
  real*8 :: mix_start
  real*8 :: mix
  real*8 :: mix_decrease
  real*8 :: step_increase
  integer :: min_mix_steps
  integer :: positive_steps
  real*8 :: maximum_displacement ! steepest-descent-guided MD only: fix the maximum displacement of atoms to prevent a cluster from exploding
  logical :: limited_displacement
  logical :: t_extrap ! was the last step extrapolated?
  logical :: t_first_pred ! is this the first step after resetting the Hessian?

  ! For BFGS
  real*8, dimension(:,:), allocatable :: hessian
  real*8, dimension(:,:), allocatable :: hessian_work
  real*8, dimension(:), allocatable :: search_direction
  real*8 :: line_step, stored_line_step
  real*8 :: force_times_searchdir_hat, stored_force_times_searchdir_hat
  logical :: ini_relaxation, line_search_on, step_was_restricted
  logical :: initial_hessian

  ! For trust radius method
  real*8 :: trust_radius

  ! For Lapack handling of BFGS only (private!!)
  integer, dimension(:), allocatable, private :: ipiv
  integer, private :: info
  integer, private :: lwork
  real*8, dimension(:), allocatable, private :: work

  ! for localorb_io
  character*100, private :: info_str

  ! for Lapack handling of cleaning_forces
  real*8, dimension(3,3), private :: tensor_of_inertia
  real*8, dimension(3,3), private :: tensor_of_inertia_lv
  real*8, dimension(3), private :: moments_of_inertia
  real*8, dimension(3), private :: moments_of_inertia_lv
  integer, private :: info_cf
  integer, private :: lwork_cf
  real*8, dimension(:), allocatable, private :: work_cf
  real*8, dimension(:,:,:), allocatable, private :: translation_vectors
  real*8, dimension(:,:,:), allocatable, private :: rotation_vectors
  real*8, dimension(:,:,:), allocatable, private :: translation_vectors_lv
  real*8, dimension(:,:,:), allocatable, private :: rotation_vectors_lv

  character*40,private :: file_name
  character*40,private,parameter :: restart_file_name='relaxation_restart_file.FHIaims'
  logical:: eigenvectors_stored


  !MOL: symmetry-constrained relaxation
  !real*8, dimension(:), allocatable      :: SCR_coords_sr_1d, SCR_forces_sr_1d
  !real*8, dimension(:), allocatable      :: SCR_lv_sr_1d, SCR_forces_lv_sr_1d
  !real*8, dimension(:,:),allocatable    :: SCR_A_lv  ! trafo matrix lattice
  !real*8, dimension(:,:),allocatable    :: SCR_A_lv_inv
  !real*8, dimension(:,:),allocatable   :: SCR_A_coords  ! trafo matrix coords
  !real*8, dimension(:,:),allocatable   :: SCR_A_coords_inv
  !real*8, dimension(:), allocatable   :: SCR_B_lv        ! vector for coords trafo
  !real*8, dimension(:), allocatable    :: SCR_B_coords    ! vector for lattice trafo


contains
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/allocate_relaxation
  !  NAME
  !    allocate_relaxation
  !  SYNOPSIS
  subroutine allocate_relaxation( )
    !  PURPOSE
    !    Allocation command
    !  USES
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! we only allocate the variables which need to be read in here (done in read_geo.f90)
    ! everything else is checked and allocated in initialize_relaxation further down.

    !DB: 10/14/13: we need to have this flag here, when empty atoms are present.
    !    in case those empty sites won't have basis sets we should set
    !    relaxation constrains true by default. in order the array constrain_relaxation
    !    is allocated, we need a .true. here.
    if (n_empty_atoms.gt.0) use_relaxation_constraints = .true.

    if (use_relaxation_constraints) then
       call aims_allocate(constrain_relaxation,n_atoms,"constrain_relaxation")
       call aims_allocate(constrain_coords,3,n_atoms,"constrain_coords")
    end if
    if (use_relaxation_constraints .or. (relax_unit_cell.eq.2)) then
       call aims_allocate(constrain_lv_relaxation,n_periodic,"constrain_lv_relaxation")
       call aims_allocate(constrain_lv_components,3,n_periodic,"constrain_lv_components")
       constrain_lv_relaxation = .false.
       constrain_lv_components = .false.
    end if

    trust_radius = 0.d0
    eigenvectors_stored = .false.

  end subroutine allocate_relaxation
  !******
  !****s* relaxation/initialize_relaxation
  !  NAME
  !    initialize_relaxation
  !  SYNOPSIS
  subroutine initialize_relaxation(coords, lattice_vector)

    !  PURPOSE
    !    initialization of all variables pertaining to the relaxation
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(in) :: coords(3, n_atoms)
    real*8, intent(in) :: lattice_vector(3, 3)

    !  INPUTS
    !    o coords -- atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! Local variables

    integer :: i_atom
    integer :: i_periodic
    integer :: i_mobile_atom
    integer :: i_coord

    ! begin work

    write(info_str,'(A)') ""
    call localorb_info(info_str,use_unit,'(A)')
    write(info_str,'(2X,A)') &
         "Initializing relaxation algorithms."
    call localorb_info(info_str,use_unit,'(A)')

    n_mobile_lv = 0
    n_mobile_lv_components = 0
    minimum_energy_ever_not_initialized = .true.

    ! first, check for constraints on relaxation
    ! The following two arrays store the indexing between the coordinate index in the Hessian
    ! (contains only mobile coordinates) and the actual atom and coordinate that this column/row
    ! in the Hessian belongs to.
    if (.not. allocated(mobile_atom)) then
       call aims_allocate(mobile_atom,3*n_atoms,"mobile_atom")
    end if
    if (.not. allocated(mobile_atom_coord)) then
       call aims_allocate(mobile_atom_coord,3*n_atoms,"mobile_atom_coord")
    end if
    if (relax_unit_cell .ne. 0) then
       if (.not. allocated(mobile_lv)) then
          call aims_allocate(mobile_lv,3*n_periodic,"mobile_lv")
       end if
       if (.not. allocated(mobile_lv_component)) then
          call aims_allocate(mobile_lv_component,3*n_periodic,"mobile_lv_component")
       end if
    end if

    ! FlK: Some safeguards when use_symm_const_geo_relaxation is used:
    if (use_symm_const_geo_relaxation) then
      use_distributed_hessian = .false.
    end if

    if (use_relaxation_constraints) then
       n_mobile_atoms = 0
       n_mobile_atom_coords = 0
       n_mobile_lv = 0
       n_mobile_lv_components = 0
       do i_atom = 1, n_atoms, 1
          if (.not. constrain_relaxation(i_atom)) then
             ! if the entire atom _is_ constrained, nothing is mobile.
             ! else add the relevant constraints as needed ...
             n_mobile_atoms = n_mobile_atoms+1
             ! check x coordinate
             if (.not. constrain_coords(1,i_atom)) then
               n_mobile_atom_coords = n_mobile_atom_coords + 1
               mobile_atom(n_mobile_atom_coords) = i_atom
               mobile_atom_coord(n_mobile_atom_coords) = 1
             end if
             ! check y coordinate
             if (.not. constrain_coords(2,i_atom)) then
               n_mobile_atom_coords = n_mobile_atom_coords + 1
               mobile_atom(n_mobile_atom_coords) = i_atom
               mobile_atom_coord(n_mobile_atom_coords) = 2
             end if
             ! check z coordinate
             if (.not. constrain_coords(3,i_atom)) then
               n_mobile_atom_coords = n_mobile_atom_coords + 1
               mobile_atom(n_mobile_atom_coords) = i_atom
               mobile_atom_coord(n_mobile_atom_coords) = 3
             end if
          end if
       end do
       ! lattice vectors
       if (relax_unit_cell .ne. 0) then
          do i_periodic = 1, n_periodic, 1
             if (.not. constrain_lv_relaxation(i_periodic)) then
                ! if the entire lattice vector _is_ constrained, nothing is mobile.
                ! else add the relevant constraints as needed ...
                n_mobile_lv = n_mobile_lv +1
                ! check x coordinate
                if (.not. constrain_lv_components(1,i_periodic)) then
                  n_mobile_lv_components = n_mobile_lv_components + 1
                  mobile_lv(n_mobile_lv_components) = i_periodic
                  mobile_lv_component(n_mobile_lv_components) = 1
                end if
                ! check y coordinate
                if (.not. constrain_lv_components(2,i_periodic)) then
                  n_mobile_lv_components = n_mobile_lv_components + 1
                  mobile_lv(n_mobile_lv_components) = i_periodic
                  mobile_lv_component(n_mobile_lv_components) = 2
                end if
                ! check z coordinate
                if (.not. constrain_lv_components(3,i_periodic)) then
                  n_mobile_lv_components = n_mobile_lv_components + 1
                  mobile_lv(n_mobile_lv_components) = i_periodic
                  mobile_lv_component(n_mobile_lv_components) = 3
                end if
             end if
          end do
       end if

    else if ((use_symm_const_geo_relaxation) .and. (init_hess_type .eq. HESS_EXPL)) then
       call reduce_index_vectors()
    else
       ! .not. use_relaxation_constraints and .not. use_symm_const_geo_relaxation
       n_mobile_atoms = n_atoms
       n_mobile_atom_coords = 3*n_atoms
       if (relax_unit_cell .ne. 0) then
          n_mobile_lv = n_periodic
          n_mobile_lv_components = 3*n_periodic
       end if

       do i_atom = 1, n_atoms, 1
          mobile_atom(3*(i_atom-1)+1:3*(i_atom-1)+3) = i_atom
          mobile_atom_coord(3*(i_atom-1)+1) = 1
          mobile_atom_coord(3*(i_atom-1)+2) = 2
          mobile_atom_coord(3*(i_atom-1)+3) = 3
       end do
       ! lattice vectors
       if (relax_unit_cell .ne. 0) then
          do i_periodic = 1, n_periodic, 1
             mobile_lv(3*(i_periodic-1)+1:3*(i_periodic-1)+3) = i_periodic
             mobile_lv_component(3*(i_periodic-1)+1) = 1
             mobile_lv_component(3*(i_periodic-1)+2) = 2
             mobile_lv_component(3*(i_periodic-1)+3) = 3
          end do
       end if
    end if ! use_relaxation_constraint

    ! Add the number of mobile atoms and lattice vectors -> dimension of the Hessian.
    n_mobile_coords = n_mobile_atom_coords + n_mobile_lv_components

    if (n_mobile_coords.gt.0) then
       if (remove_unitary_force_components .eq. 2) then
          call initialize_cleaning_forces()
       end if

       ! Next, initialize all algorithm-specific quantities for relaxation

       if (relax_mode.eq.RELAX_BFGS_TB .or. relax_mode.eq.RELAX_TRM &
          .or. relax_mode.eq.RELAX_LTRM) then

          call initialize_bfgs(coords, lattice_vector)

       end if

    else

       write(info_str,'(2X,A)') &
            "| The whole structure is constrained to be fixed."
       call localorb_info(info_str,use_unit,'(A)' )
       write(info_str,'(2X,A)') &
            "| Switching relaxation off."
       call localorb_info(info_str,use_unit,'(A)' )
       use_geo_relaxation = .false.
    end if
    ! TARP: Moved to read_geo subroutine to initialize SCR_coords_sr to non-cell centered coordinates
    ! if (use_symm_const_geo_relaxation) then
    !    call allocate_SCR()
    !    call SCR_create_coords_matrices()
    !    call SCR_create_lv_matrices()
    ! end if

  end subroutine initialize_relaxation
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/initialize_bfgs
  !  NAME
  !    initialize_bfgs
  !  SYNOPSIS

  subroutine initialize_bfgs(coords, lattice_vector)

    !  PURPOSE
    !    Initialization for BFGS relaxation
    !
    !  USES
    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    real*8, intent(in) :: coords(3, n_atoms)
    real*8, intent(in) :: lattice_vector(3, 3)

    !  INPUTS
    !    o coords -- Atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    ! TARP: Parameters to reduce Hessian
    real*8, dimension(SCR_n_params, SCR_n_params) :: hessian_sr
    real*8, dimension(3, 3) :: hessian_store, lv_inv
    integer,dimension(3):: ipivot
    integer,dimension(3):: work_ri
    integer :: info
    integer :: i_coord
    integer :: j_coord


    ! Keep all linear algebra restricted to 'mobile' atoms (if unit cell is
    ! relaxed then this degree of freedom is already added in)

    call find_hess_dimension(n_mobile_coords)
    call init_hess_id()

    if (.not. allocated(hessian)) then
       call aims_allocate(hessian,nrow_hess,ncol_hess,"hessian")
    end if

    ! TAP: Initialize the Hessian to be the full one, if using symmetry constraints
    !      and not HESS_EXPL, reduce the hessian and index vectors
    ! Initialize hessian to be unity or something else...
    ! This means the first step will be a simple steepest descent (or...).
    call bfgs_init_hessian(coords, lattice_vector)

    if (use_symm_const_geo_relaxation .and. init_hess_type .ne. HESS_EXPL) then
      ! TARP: Reduce the Hessian and index vectors
      call SCR_reduce_hessian(hessian, hessian_sr, lattice_vector)
      call reduce_index_vectors()
      ! TARP: With reduced index vectors calculate the new Hessian dimensions
      call find_hess_dimension(n_mobile_coords)
      ! TARP: Reset hess_id to be correct for symmetry-constrained relaxation
      if (init_hess_type .eq. HESS_LINDH) then
        hess_id = 0.d0
        do i_coord = 1, SCR_n_params
           hess_id(i_coord) = i_coord
        end do
      end if
      call aims_deallocate(hessian, "hessian")
      call aims_allocate(hessian, nrow_hess, ncol_hess, "hessian")
      hessian(:,:) = hessian_sr(1:nrow_hess, 1:ncol_hess)
    end if

    if (.not. allocated(hessian_work) .and. relax_mode == RELAX_BFGS_TB) then
       call aims_allocate(hessian_work,nrow_hess,ncol_hess,"hessian_work")
    end if

    if (.not. allocated(search_direction)) then
       call aims_allocate(search_direction,n_mobile_coords,"search_direction")
    end if

    ! store coordinates, forces for all atoms
    if (.not. allocated(stored_coords)) then
       call aims_allocate(stored_coords,3,n_atoms,"stored_coords")
    end if
    if (.not. allocated(stored_forces)) then
       call aims_allocate(stored_forces,3,n_atoms,"stored_forces")
    end if
    if (.not. allocated(stored_DX)) then
        call aims_allocate(stored_DX,n_mobile_coords,"stored_DX")
    end if
    if (.not. allocated(stored_coords_2)) then
       call aims_allocate(stored_coords_2,3,n_atoms,"stored_coords_2")
    end if
    if (.not. allocated(stored_forces_2)) then
       call aims_allocate(stored_forces_2,3,n_atoms,"stored_forces_2")
    end if

    ! VA: store lattice vectors and derivatives on lattice vectors
    if (.not. allocated(stored_lv)) then
       call aims_allocate(stored_lv,3,n_periodic,"stored_lv")
    end if
    if (.not. allocated(stored_forces_lv)) then
       call aims_allocate(stored_forces_lv,3,n_periodic,"stored_forces_lv")
    end if
    if (.not. allocated(stored_lv_2)) then
       call aims_allocate(stored_lv_2,3,n_periodic,"stored_lv_2,3")
    end if
    if (.not. allocated(stored_forces_lv_2)) then
       call aims_allocate(stored_forces_lv_2,3,n_periodic,"stored_forces_lv_2")
    end if

    ! FK: Store shift of atoms
    if (.not. allocated(stored_periodic_unit_cell_trans)) then
       call aims_allocate(stored_periodic_unit_cell_trans,3,n_atoms,"stored_periodic_unit_cell_trans")
    end if
    if (.not. allocated(stored_periodic_unit_cell_trans_2)) then
       call aims_allocate(stored_periodic_unit_cell_trans_2,3,n_atoms,"stored_periodic_unit_cell_trans_2")
    end if

    ! TARP: Moving this call right after Hessian is allocated. This is done so the initial Hessian for
    !      SCR relaxation can be reduced before reduce_index_vectors() is called
    ! Initialize hessian to be unity or something else...
    ! This means the first step will be a simple steepest descent (or...).
    ! call bfgs_init_hessian(coords, lattice_vector)

    initial_hessian = .true.
    t_first_pred    = .true.

    ! initialize step along line search direction to 1
    ! (the correct value for a purely quadratic extrapolation scheme)
    line_step      = 1.d0
    line_search_on = .false.

    ! Do not assume relaxation history in first iteration
    ini_relaxation = .true.

    ! For the time being, take this as initializer for the trust radius.
    ! Please note that trust_radius is meant as upper bound of the 2-norm
    ! whereas max_atomic_move in BFGS belongs to the infinum norm.

    ! If no initial value has been provided, the current trust radius should still be zero.
    if (trust_radius == 0.d0) then
      trust_radius = max_atomic_move
    end if

    ! determine work space dimensions for bfgs lapack call
    ! lapack test call only, no work is done
    if (.not. allocated(work)) then
       call aims_allocate(work,1,"work")
    end if
    if (.not. allocated(ipiv)) then
       call aims_allocate(ipiv,n_mobile_coords,"ipiv")
    end if

    ipiv  = 0
    lwork = -1

    call dsysv('u',n_mobile_coords,1,hessian,n_mobile_coords,ipiv,&
            search_direction,n_mobile_coords,work,lwork,info)

    lwork = work(1)
    call aims_deallocate(work,"work")
    call aims_allocate(work,lwork,"work")

    if (store_eigenvectors_to_disk_in_relaxation) then

       write(file_name,'(A)') 'relaxation_store'

       if (myid < 10) then
          write(file_name,'(2A,I1)') trim(file_name), '0000', myid
       else if (myid < 100) then
          write(file_name,'(2A,I2)') trim(file_name), '000', myid
       else if (myid < 1000) then
          write(file_name,'(2A,I3)') trim(file_name), '00', myid
       else if (myid < 10000) then
          write(file_name,'(2A,I4)') trim(file_name), '0', myid
       else if (myid < 100000) then
          write(file_name,'(A,I5)') trim(file_name), myid
       end if

    end if

    call read_relaxation_data_from_restart(coords, lattice_vector)

  end subroutine initialize_bfgs
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/read_relaxation_data_from_restart
  !  NAME
  !    read_relaxation_data_from_restart
  !  SYNOPSIS

  subroutine read_relaxation_data_from_restart(coords, lattice_vector)

    !  PURPOSE
    !    read the relaxation data from a previous run from file: avoids unnecessary
    !    reinitialization of the Hessian Matrix
    !  USES

    use mpi_tasks, only: myid
    use synchronize_mpi_basic, only: sync_vector, bcast_real, bcast_logical

    !  ARGUMENTS

    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)

    !  INPUTS
    !    o coords -- Atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    logical :: relaxation_restart_file_exists
    integer :: i_atom_1, i_atom_2, i
    integer :: i_periodic_1, i_periodic_2
    integer :: ios, old_relax_mode, old_relax_unit_cell
    real*8  :: this_coord(3)
    real*8  :: this_lv(3)

    if (.not. restart_relaxations) return
    if (myid == 0) then
       inquire(FILE=restart_file_name,EXIST=relaxation_restart_file_exists)
    end if
    call bcast_logical(relaxation_restart_file_exists, 0)
    if (.not. relaxation_restart_file_exists) return
    open(file = restart_file_name, unit = 88, status = 'old', form = 'unformatted', action='read')
    if (myid == 0) then
       do i_atom_1 = 1, n_mobile_coords
          do i_atom_2 = 1, n_mobile_coords
             read(88)  hessian(i_atom_1,i_atom_2)
          end do
       end do
       do i_atom_1 = 1, n_mobile_coords
          read(88)  search_direction(i_atom_1)
       end do
       do i_atom_1 = 1, n_atoms
          read(88)  (stored_coords(i,i_atom_1), i = 1,3)
       end do
       if (relax_unit_cell .ne. 0) then
       do i_periodic_1 = 1, n_periodic
          read(88)  (stored_lv(i_periodic_2,i_periodic_1), i_periodic_2 = 1,3)
       end do
       end if
       do i_atom_1 = 1, n_atoms
          read(88)  (stored_forces(i,i_atom_1), i = 1,3)
       end do
       if (relax_unit_cell .ne. 0) then
       do i_periodic_1 = 1, n_periodic
          read(88)  (stored_forces_lv(i_periodic_2,i_periodic_1), i_periodic_2 = 1,3)
       end do
       end if
       read(88) stored_energy
       read(88) initial_hessian
       read(88) ini_relaxation
       read(88) line_step
       read(88) stored_force_times_searchdir_hat
       read(88) line_search_on
       read(88) step_was_restricted
       read(88) t_extrap
       read(88) t_first_pred
       read(88, IOSTAT=ios) trust_radius
       if (ios /= 0) trust_radius = max_atomic_move
       do i_atom_1 = 1, n_atoms
          read(88, IOSTAT=ios) (this_coord(i), i = 1,3)
          if (ios /= 0) exit ! Legacy: Do not check for old restart files.
          if (any(abs(this_coord - coords(:, i_atom_1)) > 1d-5)) then
             call localorb_info('*** Trying to restart, but atoms in  geometry.in do not match!')
             call localorb_info('*** Replace geometry.in by geometry.in.next_step before restart.')
             call aims_stop
          end if
       end do
       if (n_periodic .ne. 0) then
          do i_periodic_1 = 1, n_periodic
             read(88, IOSTAT=ios) (this_lv(i_periodic_2), i_periodic_2 = 1,3)
             if (ios /= 0) exit ! Legacy: Do not check for old restart files.
             if (any(abs(this_lv - lattice_vector(:, i_periodic_1)) > 1d-5)) then
                call localorb_info('*** Trying to restart, but lattice vectors in  geometry.in do not match!')
                call localorb_info('*** Replace geometry.in by geometry.in.next_step before restart.')
                call aims_stop
             end if
          end do
       end if
       read(88, IOSTAT=ios) old_relax_mode
       if (ios == 0) then
          if (relax_mode == RELAX_BFGS_TB .and. &
          &   old_relax_mode /= RELAX_BFGS_TB) then
             call localorb_info('*** Cannot change to BFGS method during restart.')
             call aims_stop
          end if
       end if
       read(88, IOSTAT=ios) old_relax_unit_cell
       if (ios == 0) then
          if (old_relax_unit_cell .ne. relax_unit_cell) then
             call localorb_info('*** Cannot change unit cell relaxation mode during restart.')
             call aims_stop
          end if
       end if
       close(88)
       ! complete hessian (needed for BFGS->TRM transition)
       do i_atom_1 = 1, n_mobile_coords
          do i_atom_2 = i_atom_1+1, n_mobile_coords
             hessian(i_atom_2, i_atom_1) = hessian(i_atom_1,i_atom_2)
          end do
       end do
    else
       hessian             = 0d0
       search_direction    = 0d0
       stored_coords       = 0d0
       stored_lv           = 0d0
       stored_forces       = 0d0
       stored_forces_lv    = 0d0
       stored_energy       = 0d0
       ini_relaxation      = .true.
       initial_hessian     = .true.
       line_search_on      = .false.
       step_was_restricted = .false.
       line_step           = 1d0
       step_was_restricted = .false.
       t_extrap            = .false.
       t_first_pred        = .true.
       trust_radius        = 0d0
    end if

    call sync_vector(hessian, n_mobile_coords*n_mobile_coords)
    call sync_vector(search_direction, n_mobile_coords)
    call sync_vector(stored_forces, n_atoms*3)
    call sync_vector(stored_forces_lv, n_periodic*3)
    call sync_vector(stored_coords, n_atoms*3)
    call sync_vector(stored_lv, n_periodic*3)
    call bcast_real(stored_energy, 0)
    call bcast_real(line_step, 0)
    call bcast_real(stored_force_times_searchdir_hat, 0)
    call bcast_logical(ini_relaxation, 0)
    call bcast_logical(initial_hessian, 0)
    call bcast_logical(line_search_on, 0)
    call bcast_logical(step_was_restricted, 0)
    call bcast_logical(t_extrap, 0)
    call bcast_logical(t_first_pred, 0)
    call bcast_real(trust_radius, 0)

    call bfgs_allocate_eigenvectors

  end subroutine read_relaxation_data_from_restart
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/save_relaxation_data_for_restart
  !  NAME
  !    save_relaxation_data_for_restart
  !  SYNOPSIS

  subroutine save_relaxation_data_for_restart(coords, lattice_vector)

    !  PURPOSE
    !
    !    save the relaxation data for current run from file
    !
    !  USES

    use mpi_tasks, only: myid
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)

    !  INPUTS
    !    o coords -- atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_atom_1, i_atom_2
    integer:: i_periodic_1, i_periodic_2

    if (.not. restart_relaxations) return
    if (myid == 0) then
       !            inquire(FILE=restart_file_name,EXIST=restart_file_exists)
       open(file = restart_file_name, unit = 88, status = 'replace', form = 'unformatted', action='write')
       !            open(file = restart_file_name, unit =88)
       !            open(88, file = restart_file_name)
       do i_atom_1 = 1, n_mobile_coords
          do i_atom_2 = 1, n_mobile_coords
             write(88)  hessian(i_atom_1,i_atom_2)
          end do
       end do
       do i_atom_1 = 1, n_mobile_coords
          write(88)  search_direction(i_atom_1)
       end do
       do i_atom_1 = 1, n_atoms
          write(88)  (stored_coords(i_atom_2,i_atom_1), i_atom_2 = 1,3)
       end do
       if (relax_unit_cell .ne. 0) then
          do i_periodic_1 = 1, n_periodic
             write(88)  (stored_lv(i_periodic_2,i_periodic_1), i_periodic_2 = 1,3)
          end do
       end if
       do i_atom_1 = 1, n_atoms
          write(88)  (stored_forces(i_atom_2,i_atom_1), i_atom_2 = 1,3)
       end do
       if (relax_unit_cell .ne. 0) then
          do i_periodic_1 = 1, n_periodic
             write(88)  (stored_forces_lv(i_periodic_2,i_periodic_1), i_periodic_2 = 1,3)
          end do
       end if
       write(88) stored_energy
       write(88) initial_hessian
       write(88) ini_relaxation
       write(88) line_step
       write(88) stored_force_times_searchdir_hat
       write(88) line_search_on
       write(88) step_was_restricted
       write(88) t_extrap
       write(88) t_first_pred
       write(88) trust_radius
       do i_atom_1 = 1, n_atoms
          write(88)  (coords(i_atom_2,i_atom_1), i_atom_2 = 1,3)
       end do
       if (n_periodic .ne. 0)  then
          do i_periodic_1 = 1, n_periodic
             write(88)  (lattice_vector(i_periodic_2,i_periodic_1), i_periodic_2 = 1,3)
          end do
       end if
       write(88) relax_mode
       write(88) relax_unit_cell
       close(88)
    end if

  end subroutine save_relaxation_data_for_restart
  !******
  !----------------------------------------------------------------------------
  !****s* relaxation/bfgs_init_hessian
  !  NAME
  !    bfgs_init_hessian
  !  SYNOPSIS
  subroutine bfgs_init_hessian(coords, lattice_vector)
    !  PURPOSE
    !    Initialize Hessian matrix, either:
    !      * as a diagonal matrix of magnitude init_hess,
    !      * as specified by 'hessian_block' in geometry.in,
    !      * or as a sum of a diagonal matrix and the Lindh matrix.
    !     (* for the lattice:) as the square of the reciprocal lattice
    !        to effectively implement a diagonal strain transformation in
    !        the first step
    !  USES
    use mpi_tasks
    use species_data
    use lindh
    use debug_output
    use numerical_utilities
    use bravais
    use pbc_lists, only: check_geometry_for_sanity
    implicit none

    !  ARGUMENTS

    real*8, intent(in)  :: coords(3,n_atoms)
    real*8, intent(in)  :: lattice_vector(3,3)

    !  INPUTS
    !    o coords -- atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    o hessian -- Hessian matrix as used within this module.
    !                 If no file is given, use unity scaled by init_hess.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable :: Htmp(:,:)
    real*8, allocatable :: Dtmp(:,:)

    integer :: i_block
    integer :: i_atom
    integer :: j_atom
    integer :: i
    integer :: j
    integer :: i_mobile
    integer :: j_mobile
    integer :: i_periodic
    integer :: j_periodic
    integer :: i_mobile_lv
    integer :: j_mobile_lv
    integer :: i_row
    integer :: i_col
    real*8  :: recip_lattice_vector(3,3)
    real*8  :: recip_lattice_vector_2(3,3)
    real*8  :: hessian_temp(9,9)
    real*8  :: asym

    character*150 :: info_str

    hessian = 0.d0
    hessian_temp = 0.d0

    call check_geometry_for_sanity()

    if (.not. use_distributed_hessian) then
       select case (init_hess_type)
       case (HESS_EXPL) ! Initialize Hessian matrix from hess_blocks
          ! Atom
          do i_block = 1, n_hess_blocks
             i_atom = hess_block_pos(1, i_block)
             j_atom = hess_block_pos(2, i_block)
             do i = 1, 3
                do i_mobile = 1, n_mobile_atom_coords
                   if (i_atom == mobile_atom(i_mobile) .and. &
                   &   i == mobile_atom_coord(i_mobile)) exit
                end do
                if (i_mobile > n_mobile_atom_coords) cycle ! not needed
                do j = 1, 3
                   do j_mobile = 1, n_mobile_atom_coords
                      if (j_atom == mobile_atom(j_mobile) .and. &
                      &   j == mobile_atom_coord(j_mobile)) exit
                   end do
                   if (j_mobile > n_mobile_atom_coords) cycle
                   hessian(i_mobile, j_mobile) = hess_blocks(i, j, i_block)
                   hessian(j_mobile, i_mobile) = hess_blocks(i, j, i_block)
                end do
             end do
          end do

          ! Lattice
          do i_block = 1, n_hess_blocks_lv
             i_periodic = hess_block_pos_lv(1, i_block)
             j_periodic = hess_block_pos_lv(2, i_block)
             do i = 1, 3
                do i_mobile_lv = 1, n_mobile_lv_components
                   if (i_periodic == mobile_lv(i_mobile_lv) .and. &
                   &   i == mobile_lv_component(i_mobile_lv)) exit
                end do
                if (i_mobile_lv > n_mobile_lv_components) cycle ! not needed
                do j = 1, 3
                   do j_mobile_lv = 1, n_mobile_lv_components
                      if (j_periodic == mobile_lv(j_mobile_lv) .and. &
                      &   j == mobile_lv_component(j_mobile_lv)) exit
                   end do
                   if (j_mobile_lv > n_mobile_lv_components) cycle
                   hessian(i_mobile_lv+n_mobile_atom_coords, &
                   &       j_mobile_lv+n_mobile_atom_coords) &
                           = hess_blocks_lv(i, j, i_block)
                   hessian(j_mobile_lv+n_mobile_atom_coords, &
                   &       i_mobile_lv+n_mobile_atom_coords) &
                           = hess_blocks_lv(i, j, i_block)
                end do
             end do
          end do

          ! Lattice-atom
          do i_block = 1, n_hess_blocks_lv_atom
             i_periodic = hess_block_pos_lv_atom(1, i_block)
             j_atom = hess_block_pos_lv_atom(2, i_block)
             do i = 1, 3
                do i_mobile_lv = 1, n_mobile_lv_components
                   if (i_periodic == mobile_lv(i_mobile_lv) .and. &
                   &   i == mobile_lv_component(i_mobile_lv)) exit
                end do
                if (i_mobile_lv > n_mobile_lv_components) cycle ! not needed
                do j = 1, 3
                   do j_mobile = 1, n_mobile_atom_coords
                      if (j_atom == mobile_atom(j_mobile) .and. &
                      &   j == mobile_atom_coord(j_mobile)) exit
                   end do
                   if (j_mobile > n_mobile_atom_coords) cycle
                   hessian(i_mobile_lv+n_mobile_atom_coords, j_mobile) &
                        =  hess_blocks_lv_atom(i, j, i_block)
                   hessian(j_mobile, i_mobile_lv+n_mobile_atom_coords) &
                        = hess_blocks_lv_atom(i, j, i_block)
                end do
             end do
          end do

       case (HESS_DIAG) ! Set diagonal init_hess matrix
          do i_mobile = 1, n_mobile_coords
             if (i_mobile <= n_mobile_atom_coords) then
                hessian(i_mobile, i_mobile) = hessian(i_mobile, i_mobile) &
                &                             + init_hess_diag
             else
                ! This must be a lattice vector coordinate. Here we need a
                ! potentially different initial diagonal element
                hessian(i_mobile, i_mobile) = hessian(i_mobile, i_mobile) &
                &                             + init_hess_lv_diag
             end if
          end do

       case (HESS_LINDH)
          call get_lindh_hessian(coords, lattice_vector, &
                  init_hess_lindh_thres, hessian)

          ! Add diagonal
          do i_mobile = 1, n_mobile_coords
             if (i_mobile <= n_mobile_atom_coords) then
                hessian(i_mobile, i_mobile) = hessian(i_mobile, i_mobile) &
                &                             + init_hess_lindh_diag
             else
                ! This must be a lattice vector coordinate. Here we need a
                ! potentially much larger initial diagonal element
                hessian(i_mobile, i_mobile) = hessian(i_mobile, i_mobile) &
                &                             + init_hess_lv_diag
             end if
          end do

          ! FlK: initialize hessian as in case(HESS_LATTICE) when lattice_trm is
          ! used:
          if ((relax_geo .eq. RELAX_LTRM) .and. (n_periodic > 0)) then
            call get_reciprocal_vectors( &
              n_periodic, lattice_vector, recip_lattice_vector)

            recip_lattice_vector_2 = &
              matmul(recip_lattice_vector, transpose(recip_lattice_vector))

            do j = 1, 3
              do i = 1, 3
                do i_block = 0, 2
                  hessian_temp(i+i_block*3, j+i_block*3) &
                    = recip_lattice_vector_2(i, j) &
                    * init_hess_lv_diag &
                    / frobnorm_3x3_real(recip_lattice_vector_2)
                end do
              end do
            end do

            ! Assign the components of the hessian that are allowed to move
            do j_mobile_lv = 1, n_mobile_lv_components
              do i_mobile_lv = 1, n_mobile_lv_components
                ! assign x,y,z from the mobile coordinates to match the
                ! correct line from the full initial hessian
                i = (mobile_lv(i_mobile_lv) - 1) * 3 + mobile_lv_component(i_mobile_lv)
                j = (mobile_lv(j_mobile_lv) - 1) * 3 + mobile_lv_component(j_mobile_lv)

                hessian(i_mobile_lv+n_mobile_atom_coords, j_mobile_lv+n_mobile_atom_coords) &
                  = hessian_temp(i, j)
              end do
            end do
          end if

       case (HESS_LATTICE)
         ! FlK:
         ! Initialize Hessian with the reciprocal lattice to effectivetly
         ! implement a diagonal Hessian in strain coordinates.

         ! Set diagonal hessian for atomic positions
         do i_mobile = 1, n_mobile_atom_coords
           hessian(i_mobile, i_mobile) = hessian(i_mobile, i_mobile) &
                                         + init_hess_diag
         end do

         call get_reciprocal_vectors(n_periodic, lattice_vector, recip_lattice_vector)

         ! square of reciprocal lattice
         recip_lattice_vector_2 = &
            matmul(recip_lattice_vector, transpose(recip_lattice_vector))

         ! choose
         !   H = c * lattice^-2
         ! and choose c such that |H| = |init_hess_lv_diag * unit_matrix|
         do j = 1, 3
           do i = 1, 3
             do i_block = 0, 2
               hessian_temp(i+i_block*3, j+i_block*3) &
                 = recip_lattice_vector_2(i, j) &
                 * init_hess_lv_diag &
                 / frobnorm_3x3_real(recip_lattice_vector_2)
             end do
           end do
         end do

         ! Assign the components of the hessian that are allowed to move
         do j_mobile_lv = 1, n_mobile_lv_components
           do i_mobile_lv = 1, n_mobile_lv_components
             ! assign x,y,z from the mobile coordinates to match the
             ! correct line from the full initial hessian
             i = (mobile_lv(i_mobile_lv) - 1) * 3 + mobile_lv_component(i_mobile_lv)
             j = (mobile_lv(j_mobile_lv) - 1) * 3 + mobile_lv_component(j_mobile_lv)

             hessian(i_mobile_lv+n_mobile_atom_coords, j_mobile_lv+n_mobile_atom_coords) &
                = hessian_temp(i, j)
           end do
         end do
       end select

       ! Check symmetry
       asym = maxval(abs(hessian - transpose(hessian)))
       if (asym > 1.d-10) then
          write(info_str, "('Hessian is asymmetric by',ES10.2,' eV/A^2.')") &
          & asym*hartree/bohr**2
          call localorb_info(info_str, use_unit, "(2X,'** ',A)")
       end if

    else if (is_worker) then ! use_distributed_hessian
       select case (init_hess_type)
       case (HESS_EXPL) ! Initialize Hessian matrix from file
          if (hess_in_file) then
             call read_previous_hess(hessian)
          else
             ! Atom-atom
             do i_block = 1, n_hess_blocks
                i_atom = hess_block_pos(1, i_block)
                j_atom = hess_block_pos(2, i_block)

                do i = 1, 3
                   do i_mobile = 1, n_mobile_atom_coords
                      if (i_atom == mobile_atom(i_mobile) &
                         .and. i == mobile_atom_coord(i_mobile)) then
                         exit
                      end if
                   end do

                   if (i_mobile > n_mobile_atom_coords) then
                      cycle
                   end if

                   do j = 1, 3
                      do j_mobile = 1, n_mobile_atom_coords
                         if (j_atom == mobile_atom(j_mobile) &
                            .and. j == mobile_atom_coord(j_mobile)) then
                            exit
                         end if
                      end do

                      if (j_mobile > n_mobile_atom_coords) then
                         cycle
                      end if

                      i_row = i_mobile
                      i_col = j_mobile

                      if (i_row > offset .and. i_row <= offset+blk_hess) then
                         hessian(i_row-offset, i_col) &
                            = hess_blocks(i, j, i_block)
                      end if

                      if (i_col > offset .and. i_col <= offset+blk_hess) then
                         hessian(i_col-offset, i_row) &
                            = hess_blocks(i, j, i_block)
                      end if
                   end do
                end do
             end do

             ! Lattice-lattice
             do i_block = 1, n_hess_blocks_lv
                i_periodic = hess_block_pos_lv(1, i_block)
                j_periodic = hess_block_pos_lv(2, i_block)

                do i = 1, 3
                   do i_mobile_lv = 1, n_mobile_lv_components
                      if (i_periodic == mobile_lv(i_mobile_lv) &
                         .and. i == mobile_lv_component(i_mobile_lv)) then
                         exit
                      end if
                   end do

                   if (i_mobile_lv > n_mobile_lv_components) then
                      cycle
                   end if

                   do j = 1, 3
                      do j_mobile_lv = 1, n_mobile_lv_components
                         if (j_periodic == mobile_lv(j_mobile_lv) &
                            .and. j == mobile_lv_component(j_mobile_lv)) then
                            exit
                         end if
                      end do

                      if (j_mobile_lv > n_mobile_lv_components) then
                         cycle
                      end if

                      i_row = i_mobile_lv+n_mobile_atom_coords
                      i_col = j_mobile_lv+n_mobile_atom_coords

                      if (i_row > offset .and. i_row <= offset+blk_hess) then
                         hessian(i_row-offset, i_col) &
                            = hess_blocks_lv(i, j, i_block)
                      end if

                      if (i_col > offset .and. i_col <= offset+blk_hess) then
                         hessian(i_col-offset, i_row) &
                            = hess_blocks_lv(i, j, i_block)
                      end if
                   end do
                end do
             end do

             ! Lattice-atom
             do i_block = 1, n_hess_blocks_lv_atom
                i_periodic = hess_block_pos_lv_atom(1, i_block)
                j_atom = hess_block_pos_lv_atom(2, i_block)

                do i = 1, 3
                   do i_mobile_lv = 1, n_mobile_lv_components
                      if (i_periodic == mobile_lv(i_mobile_lv) &
                         .and. i == mobile_lv_component(i_mobile_lv)) then
                         exit
                      end if
                   end do

                   if (i_mobile_lv > n_mobile_lv_components) then
                      cycle
                   end if

                   do j = 1, 3
                      do j_mobile = 1, n_mobile_atom_coords
                         if (j_atom == mobile_atom(j_mobile) &
                            .and. j == mobile_atom_coord(j_mobile)) then
                            exit
                         end if
                      end do

                      if (j_mobile > n_mobile_atom_coords) then
                         cycle
                      end if

                      i_row = i_mobile_lv+n_mobile_atom_coords
                      i_col = j_mobile

                      if (i_row > offset .and. i_row <= offset+blk_hess) then
                         hessian(i_row-offset, i_col) &
                            = hess_blocks_lv_atom(i, j, i_block)
                      end if

                      if (i_col > offset .and. i_col <= offset+blk_hess) then
                         hessian(i_col-offset, i_row) &
                            = hess_blocks_lv_atom(i, j, i_block)
                      end if
                   end do
                end do
             end do
          end if

       case (HESS_DIAG) ! Set diagonal matrix
          do i_row = 1, nrow_hess
             if (i_row+offset <= n_mobile_atom_coords) then ! Atom
                hessian(i_row, i_row+offset) = init_hess_diag
             else if (i_row+offset <= n_mobile_coords) then ! Lattice
                hessian(i_row, i_row+offset) = init_hess_lv_diag
             end if
          end do

       case (HESS_LINDH)
          call get_lindh_hessian(coords, lattice_vector, &
               init_hess_lindh_thres, hessian)

          ! Add diagonal
          do i_row = 1, nrow_hess
             if (i_row+offset <= n_mobile_atom_coords) then ! Atom
                hessian(i_row, i_row+offset) = hessian(i_row, i_row+offset) &
                   + init_hess_lindh_diag
             else if (i_row+offset <= n_mobile_coords) then ! Lattice
                hessian(i_row, i_row+offset) = hessian(i_row, i_row+offset) &
                   + init_hess_lv_diag
             end if
          end do

          ! lattice_trm
          if ((relax_geo == RELAX_LTRM) .and. (n_periodic > 0)) then
             call get_reciprocal_vectors(n_periodic, lattice_vector, &
                  recip_lattice_vector)

             recip_lattice_vector_2 = matmul(recip_lattice_vector, &
                transpose(recip_lattice_vector))

             do j = 1, 3
                do i = 1, 3
                   do i_block = 0, 2
                      hessian_temp(i+i_block*3, j+i_block*3) &
                         = recip_lattice_vector_2(i, j) * init_hess_lv_diag &
                         / frobnorm_3x3_real(recip_lattice_vector_2)
                   end do
                end do
             end do

             ! Assign the components of the hessian that are allowed to move
             do j_mobile_lv = 1, n_mobile_lv_components
                do i_mobile_lv = 1, n_mobile_lv_components
                   i = (mobile_lv(i_mobile_lv) - 1) * 3 &
                      + mobile_lv_component(i_mobile_lv)
                   j = (mobile_lv(j_mobile_lv) - 1) * 3 &
                      + mobile_lv_component(j_mobile_lv)

                   i_row = i_mobile_lv + n_mobile_atom_coords
                   i_col = j_mobile_lv + n_mobile_atom_coords

                   if ((i_row > offset) .and. (i_row <= offset+blk_hess)) then
                      hessian(i_row-offset, i_col) = hessian_temp(i, j)
                   end if
                end do
             end do
          end if

       case default
          write(info_str,'(X,A)') "*** Invalid init_hess_type."
          call localorb_info(info_str,use_unit,'(A)')
          call aims_stop()

       end select

    end if ! use_distributed_hessian

    ! Hessian read from geometry.in can be deallocated now.
    if (allocated(hess_blocks)) then
       call aims_deallocate(hess_blocks, "hess_blocks")
    end if

    if (allocated(hess_blocks_lv)) then
       call aims_deallocate(hess_blocks_lv, "hess_blocks_lv")
    end if

    if (allocated(hess_blocks_lv_atom)) then
       call aims_deallocate(hess_blocks_lv_atom, "hess_blocks_lv_atom")
    end if

    if (allocated(hess_block_pos)) then
       call aims_deallocate(hess_block_pos, "hess_block_pos")
    end if

    if (allocated(hess_block_pos_lv)) then
       call aims_deallocate(hess_block_pos_lv, "hess_block_pos_lv")
    end if

    if (allocated(hess_block_pos_lv_atom)) then
       call aims_deallocate(hess_block_pos_lv_atom, "hess_block_pos_lv_atom")
    end if

  end subroutine bfgs_init_hessian
  !******
  !----------------------------------------------------------------------------
  !****s* relaxation/initialize_cleaning_forces
  !  NAME
  !    initialize_cleaning_forces
  !  SYNOPSIS
  subroutine initialize_cleaning_forces

    !  PURPOSE
    !    initialize force cleaning by sayvetz conditions
    !
    !  USES
    implicit none
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! determine work space dimensions for lapack call
    ! in remove_translation_and_rotation
    ! lapack test call only, no work is done

    ! local variables
    real*8 :: inv_norm

    ! functions
    real*8 :: dnrm2

    ! counter
    integer :: i_atom
    integer :: i_periodic
    integer :: i_coord
    integer :: i_coord_2

    if (.not. allocated(work_cf)) then
       call aims_allocate(work_cf,1,"work_cf")
    end if
    if (.not. allocated(translation_vectors)) then
       call aims_allocate(translation_vectors,3,n_atoms,3,"translation_vectors")
    end if
    if (.not. allocated(rotation_vectors)) then
       call aims_allocate(rotation_vectors,3,n_atoms,3,"rotation_vectors")
    end if

    if (n_periodic .ne. 0) then
       if (.not. allocated(translation_vectors_lv)) then
          call aims_allocate(translation_vectors_lv,3,n_periodic,3,"translation_vectors_lv")
       end if
       if (.not. allocated(rotation_vectors_lv)) then
          call aims_allocate(rotation_vectors_lv,3,n_periodic,3,"rotation_vectors_lv")
       end if
    end if

    tensor_of_inertia = 1.d0
    lwork_cf = -1

    call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work_cf, lwork_cf, info)

    ! Now allocate work to the correct size
    lwork_cf = work_cf(1)
    call aims_deallocate(work_cf,"work_cf")
    call aims_allocate(work_cf,lwork_cf,"lwork_cf")

    ! initialize vectors of translation since they are configuration independent
    translation_vectors = 0.d0
    inv_norm = 1.d0 / dsqrt(dble(n_atoms))
    do i_atom = 1, n_atoms, 1
       do i_coord = 1, 3, 1
          translation_vectors(i_coord, i_atom, i_coord) = inv_norm
       end do
    end do

    ! Similar for lattice vectors
    if (n_periodic .ne. 0) then
       translation_vectors_lv = 0.d0
       inv_norm = 1.d0 / dsqrt(dble(n_periodic))
       do i_periodic = 1, n_periodic, 1
          do i_coord = 1, 3, 1
             translation_vectors_lv(i_coord, i_periodic, i_coord) = inv_norm
          end do
       end do
   end if

  end subroutine initialize_cleaning_forces
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/symmetrize_hessian
  !  NAME
  !    symmetrize_hessian
  !  SYNOPSIS
  !
  subroutine symmetrize_hessian()

    !  PURPOSE
    !    Overwrites the lower triangle of the Hessian matrix with values from
    !    the upper triangle.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Dec 2017.
    !  SOURCE

    ! FlK Q: why not averaging?

    integer :: i_count

    real*8, dimension(:,:), allocatable :: trans_hess

    if (.not. use_distributed_hessian) then
       do i_count = 1, ncol_hess
          hessian(i_count+1:ncol_hess,i_count) = &
             hessian(i_count,i_count+1:ncol_hess)
       end do
    else if (is_worker) then
       call aims_allocate(trans_hess, nrow_hess, ncol_hess, "trans_hess")

       call transpose_hess(hessian, trans_hess)

       if (offset > 0) then
          hessian(1:nrow_hess,1:offset) = trans_hess(1:nrow_hess,1:offset)
       end if

       do i_count = 1, nrow_hess
          hessian(i_count:nrow_hess,:) = trans_hess(i_count:nrow_hess,:)
       end do

       call aims_deallocate(trans_hess, "trans_hess")
    end if

  end subroutine symmetrize_hessian
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/write_new_geometry_file
  !  NAME
  !    write_new_geometry_file
  !  SYNOPSIS

  subroutine write_new_geometry_file(coords, lattice_vector)

    !  PURPOSE
    !    write geometry file of the updated structure for relaxation restart
    !  USES

    use runtime_choices
    use species_data
    use mpi_tasks

    !  ARGUMENTS

    real*8, intent(in) :: coords(3, n_atoms)
    real*8, intent(in) :: lattice_vector(3, 3)

    !  INPUTS
    !    o coords -- Current atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_atom, j_atom, i_mobile, j_mobile, i, j
    integer:: i_periodic, j_periodic, i_mobile_lv
    real*8 :: block_3x3(3,3)
    character(len=80) :: info_str

    if (write_restart_geometry) then

       if (myid == 0) then
          write(use_unit,'(2X,A)') 'Writing the current geometry to file "geometry.in.next_step".'

          open(88, file='geometry.in.next_step')
          write(88,'(A)') '# '
          write(88,'(A)') '# This is the geometry file that corresponds to the current relaxation step.'
          write(88,'(A)') '# If you do not want this file to be written, set the "write_restart_geometry" flag to .false.'
          call write_aims_uuid(info_str)
          write(88,'(A,2X,A)') '#', info_str
          write(88,'(A)') '# '
          call write_new_geometry_content(88, coords, lattice_vector, flag_set_vacuum_level, vacuum_z_level)

          if (hessian_to_restart_geometry) then
             write(88,'(A)') '# '
             write(88,'(A)') '# What follows is the current estimated Hessian matrix constructed by the BFGS algorithm.'
             write(88,'(A)') '# This is NOT the true Hessian matrix of the system.'
             write(88,'(2A)') '# If you do not want this information here, ', &
                'switch it off using the "hessian_to_restart_geometry" keyword.'
             write(88,'(A)') '# '

             ! Remember the current trust radius
             write(88,"('trust_radius          ',F15.10)") &
                trust_radius*bohr

             if (use_distributed_hessian) then
                write(88,'(A)') 'hessian_file'
             else
                ! atoms
                do i_atom = 1, n_atoms
                   do j_atom = i_atom, n_atoms
                      block_3x3 = 0.d0
                      do i_mobile = 1, n_mobile_atom_coords
                         if (i_atom /= mobile_atom(i_mobile)) cycle
                         i = mobile_atom_coord(i_mobile)
                         do j_mobile = 1, n_mobile_atom_coords
                            if (j_atom /= mobile_atom(j_mobile)) cycle
                            j = mobile_atom_coord(j_mobile)
                            block_3x3(i, j) = hessian(i_mobile, j_mobile)
                         end do
                      end do
                      if (any(block_3x3 /= 0.d0)) then
                         write(88, "('hessian_block         ',2I8,9(1X,F18.6))") &
                         & i_atom, j_atom, block_3x3*hartree/bohr**2
                      end if
                   end do
                end do
                ! Lattice
                do i_periodic = 1, n_periodic
                   do j_periodic = i_periodic, n_periodic
                      block_3x3 = 0.d0
                      do i_mobile = 1, n_mobile_lv_components
                         if (i_periodic /= mobile_lv(i_mobile)) cycle
                         i = mobile_lv_component(i_mobile)
                         do j_mobile = 1, n_mobile_lv_components
                         if (j_periodic /= mobile_lv(j_mobile)) cycle
                            j = mobile_lv_component(j_mobile)
                            block_3x3(i, j) = &
                            hessian(i_mobile+n_mobile_atom_coords, j_mobile+n_mobile_atom_coords)
                         end do
                      end do
                      if (any(block_3x3 /= 0.d0)) then
                         write(88, "('hessian_block_lv      ',2I8,9(1X,F18.6))") &
                         & i_periodic, j_periodic, block_3x3*hartree/bohr**2
                      end if
                   end do
                end do
                ! Lattice-atom
                do i_periodic = 1, n_periodic
                   do j_atom  = 1, n_atoms
                      block_3x3 = 0.d0
                      do i_mobile_lv = 1, n_mobile_lv_components
                         if (i_periodic /= mobile_lv(i_mobile_lv)) cycle
                         i = mobile_lv_component(i_mobile_lv)
                         do j_mobile = 1, n_mobile_atom_coords
                            if (j_atom /= mobile_atom(j_mobile)) cycle
                            j = mobile_atom_coord(j_mobile)
                            block_3x3(i, j) = &
                            hessian(i_mobile_lv+n_mobile_atom_coords, j_mobile)
                         end do
                      end do
                      if (any(block_3x3 /= 0.d0)) then
                         write(88, "('hessian_block_lv_atom ',2I8,9(1X,F18.6))") &
                         & i_periodic, j_atom, block_3x3*hartree/bohr**2
                      end if
                   end do
                end do
             end if ! use_distributed_hessian
          end if ! hessian_to_restart_geometry

          close(88)

       end if ! myid == 0

       if (hessian_to_restart_geometry .and. use_distributed_hessian) then
           call write_current_hess(hessian)
       end if

    end if ! write_restart_geometry

  end subroutine write_new_geometry_file
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/write_new_geometry_content
  !  NAME
  !    write_new_geometry_content
  !  SYNOPSIS

  subroutine write_new_geometry_content(fileno, coords, lattice_vector, flag_set_vacuum_level, vacuum_level)

    !  PURPOSE
    !    write content into geometry file of the updated structure
    !    for relaxation restart
    !  USES

    use species_data
    use mpi_tasks
    use geometry, only : periodic_unit_cell_translations, coord_basis_atoms, &
                         multipole_coords, multipole_charge, multipole_order, &
                         frac2cart, cart2frac
                  ! Do not use any more. coords, lattice_vector, etc. are already passed
                  ! as explicit arguments.
    use pseudodata, only : pp_species2species

    !  ARGUMENTS

    integer, intent(IN) :: fileno
    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)
    real*8, intent(IN) :: vacuum_level
    logical, intent(IN) :: flag_set_vacuum_level

    !  INPUTS
    !    o fileno -- file descriptor for output
    !    o coords -- atomic coordinates
    !    o lattice_vector -- Bravais lattice
    !    o flag_set_vacuum_level -- has the vacuum level been set_vacuum_level
    !    o vacuum_level -- coordinate of the used vacuum level
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    integer:: i_atom, i_periodic, i_multipole, i_coord, i_pp_atom

    real*8, dimension(3,n_atoms) :: coords_temp
    real*8, dimension(3,n_atoms) :: frac_coords_temp
    real*8, dimension(3,n_atoms) :: periodic_unit_cell_translations_cart

    character(*), parameter :: func = 'write_new_geometry_content'

    ! In the periodic case, we want the output coordinates to be close to the user input
    ! coordinates, not the internal (possibly shifted to Wigner Seitz cell) coordinate that
    ! the code used. Thus, we need to generate this coordinate. For convenience (negligible cost)
    ! we do this here for absolute length (cartesian in bohr radii) and fractional coordinates.
    if (n_periodic.gt.0) then

       ! Create unshifted Cartesian coordinates
       call frac2cart (lattice_vector,                    &
                       periodic_unit_cell_translations,   &
                       periodic_unit_cell_translations_cart)
       coords_temp(:,:) = coords(:,:)+periodic_unit_cell_translations_cart(:,:)

       ! ... and create fractional coordinates from the unshifted Cartesian coordinates.
       call cart2frac(lattice_vector, coords_temp, frac_coords_temp)

    else
       coords_temp(:,:) = coords(:,:)
    end if

    if (myid == 0) then
       if (flag_set_vacuum_level) then
          write(fileno,'(2X,A,F16.6)')'set_vacuum_level ', vacuum_level*bohr
       end if

       if ( allocated (homogeneous_field) ) then
          write(fileno,'(2X,A,3F16.6)')'homogeneous_field ',   homogeneous_field(1:3)*hartree_over_bohr
       end if

       if (n_periodic > 0) then
          do i_periodic = 1, n_periodic
             write(fileno,'(A,3F16.8)') &
             'lattice_vector', lattice_vector(:,i_periodic)*bohr

             if (use_relaxation_constraints) then
                if (constrain_lv_relaxation(i_periodic)) then
                   write(fileno,'(A)') '  constrain_relaxation .true.'
                else
                   if (constrain_lv_components (1, i_periodic)) then
                      write(fileno,'(A)') '  constrain_relaxation x'
                   end if
                   if (constrain_lv_components (2, i_periodic)) then
                      write(fileno,'(A)') '  constrain_relaxation y'
                   end if
                   if (constrain_lv_components (3, i_periodic)) then
                      write(fileno,'(A)') '  constrain_relaxation z'
                   end if
                end if
             end if
          end do
       end if

       do i_atom = 1, n_atoms

          if (coord_basis_atoms(i_atom).eq.0) then
            ! cartesian coordinates used in input
            write(fileno,'(A,3F16.8,A,A)') &
            & 'atom ', coords_temp(:,i_atom)*bohr, ' ',trim(species_name(species(i_atom)))
          else if (coord_basis_atoms(i_atom).eq.1) then
            ! fractional coordinates used in input
            write(fileno,'(A,3F16.8,A,A)') &
            & 'atom_frac ', frac_coords_temp(:,i_atom), ' ',trim(species_name(species(i_atom)))
          else
            ! no other meaningful value at the time of this writing
            write(use_unit,*) "* Warning! While writing the file 'geometry.in.next_step', "
            write(use_unit,*) "* found an atom for which the input was neither fractional nor integer coordinates."
            write(use_unit,*) "* The atom in question was atom number ", i_atom, " ."
            write(use_unit,*) "* Stopping because this behaviour indicates a programming error."
            call aims_stop_coll('', func)
          end if

          if (use_relaxation_constraints) then
             if (constrain_relaxation(i_atom)) then

                write(fileno,'(A)') '  constrain_relaxation .true.'

             else
               if (constrain_coords(1,i_atom)) then
                 write(fileno,'(A)') '  constrain_relaxation x'
               end if
               if (constrain_coords(2,i_atom)) then
                 write(fileno,'(A)') '  constrain_relaxation y'
               end if
               if (constrain_coords(3,i_atom)) then
                 write(fileno,'(A)') '  constrain_relaxation z'
               end if

             end if
          end if

          if (external_forces_on) then
             if (any(abs(external_forces(:,i_atom)) > 1.d-10)) then
                write(fileno,'(A,3F16.8)') '  external_force  ', &
                    external_forces(:,i_atom) * hartree_over_bohr
             end if
          end if

       end do

       if (n_pp_atoms.gt.0 .and. (.not. add_embedding_grids)) then
          do i_pp_atom = 1, n_pp_atoms, 1
             write(fileno,'(A,3F16.8,A,A)') &
             & 'pseudocore', pp_coords(:,i_pp_atom)*bohr, ' ', &
              trim(species_name(species(n_atoms + i_pp_atom)))
          end do
       end if

       ! not entirely complete information yet - only monopoles.
       ! Written for completeness only, for now.
       ! If this becomes an issue in production, please add missing output.
       if (n_multipoles.gt.0) then
          do i_multipole = 1, n_multipoles, 1
             write(fileno,'(7X,A,3(2X,F16.8),2X,I4,2X,F10.3)') &
               "multipole ", &
               (multipole_coords(i_coord,i_multipole)*bohr, i_coord=1,3,1), &
                multipole_order(i_multipole), multipole_charge(i_multipole)
          end do
       end if

    end if

  end subroutine write_new_geometry_content
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/clean_hessian
  !  NAME
  !    clean_hessian
  !  SYNOPSIS

  subroutine clean_hessian(hessian, mobile_ref_coords)
    !  PURPOSE
    !    Clean Hessian from translations and rotations.
    !    Effectively, each row and column is cleaned just like a force
    !    component.
    !    For the BFGS method, this is a no-op.
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, intent(inout) :: hessian(nrow_hess, ncol_hess)
    real*8, intent(in)    :: mobile_ref_coords(n_mobile_coords)

    !  INPUTS
    !    o total_mobile_forces -- total forces before cleaning
    !    o mobile_ref_coords -- corresponding geometry
    !  OUTPUTS
    !    o total_mobile_forces -- total forces after cleaning
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: i_mobile_coord

    real*8, dimension(:,:), allocatable :: trans_hess

    ! BFGS would get troubled by zero modes.
    if (relax_geo == RELAX_BFGS_TB) then
       return
    end if

    if (.not. use_distributed_hessian) then
       ! Columns
       do i_mobile_coord = 1, n_mobile_coords
          call clean_mobile_force_components(hessian(:, i_mobile_coord), mobile_ref_coords)
       end do

       ! Rows
       do i_mobile_coord = 1, n_mobile_coords
          call clean_mobile_force_components(hessian(i_mobile_coord, :), mobile_ref_coords)
       end do
    else if (is_worker) then
       ! Columns
       do i_mobile_coord = 1, nrow_hess
          call clean_mobile_force_components(hessian(i_mobile_coord, :), mobile_ref_coords)
       end do

       call aims_allocate(trans_hess, nrow_hess, ncol_hess, "trans_hess")

       call transpose_hess(hessian, trans_hess)

       hessian = trans_hess

       call aims_deallocate(trans_hess, "trans_hess")

       ! Rows
       do i_mobile_coord = 1, nrow_hess
          call clean_mobile_force_components(hessian(i_mobile_coord, :), mobile_ref_coords)
       end do
    end if

  end subroutine clean_hessian
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/clean_force_components
  !  NAME
  !    clean_force_components
  !  SYNOPSIS

  subroutine clean_force_components(total_forces, forces_lv)

    !  PURPOSE
    !    Wrapper around relaxation/clean_force_components_work()
    !    to avoid explicit dependence.
    !  USES

    use runtime_choices
    use dimensions
    use localorb_io
    implicit none

    !  ARGUMENTS

    real*8, dimension(3,n_atoms),    intent(INOUT) :: total_forces
    real*8, dimension(3,n_periodic), intent(INOUT) :: forces_lv

    !  INPUTS
    !    o total_forces - input forces straight from scf_solver
    !  OUTPUT
    !    o total_force - cleaned forces
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    logical :: actually_remove
    logical :: actually_remove_lv
    character*100 :: info_str

    ! VA: more conditions if unit cell is relaxed
    actually_remove    = (.not. use_geo_relaxation) .or. ((remove_unitary_force_components == 2 .and. &
                          (n_mobile_atom_coords > (3*n_atoms-6)) ))

    if (use_qmmm) then
       actually_remove    = ((.not. use_geo_relaxation) .or. (remove_unitary_force_components == 2))
    end if

    ! VA: at a later point (analytical stress tensor) one should also consider the case where
    !     no actual relaxation is performed -- just like in the atomic case 3 lines above
    actually_remove_lv = (remove_unitary_force_components == 2 .and. &
                          (n_mobile_lv_components > (3*n_periodic-6)) .and. &
                          relax_unit_cell > 1)

    if (actually_remove) then
       write(info_str,'(2X,A)') &
       &     "Removing unitary transformations (pure translations, rotations) from forces on atoms."
       call localorb_info(info_str)
       call output_unitary_components('Atomic forces before filtering:', total_forces)
    end if
    if (actually_remove_lv) then
       write(info_str,'(2X,A)') &
       &     "Removing unitary transformations (rotations) from forces on lattice."
       call localorb_info(info_str)
       call output_unitary_components_lv('Lattice forces before filtering:', forces_lv)
    end if

    call clean_force_components_work(coords, lattice_vector, total_forces, forces_lv)

    if (actually_remove) then
       call output_unitary_components('Atomic forces after filtering:', total_forces)
    end if
    if (actually_remove_lv) then
       call output_unitary_components_lv('Lattice forces after filtering:', forces_lv)
    end if

  end subroutine clean_force_components
  !******
  !----------------------------------------------------------------------------
  !****s* relaxation/output_unitary_components
  !  NAME
  !    output_unitary_components
  !  SYNOPSIS

  subroutine output_unitary_components(msg, total_forces)

    !  PURPOSE
    !    Output magnitude of unitary components.
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: msg
    real*8, intent(IN) :: total_forces(3, n_atoms)

    !  INPUTS
    !    o msg -- Introductory message
    !  OUTPUT
    !    o total_force -- forces to check
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: center_of_mass(3), relative_coords(3, n_atoms)
    real*8 :: net_force(3), net_torque(3)
    real*8 :: conversion
    character*100 :: info_str

    ! begin work

    ! first get center of mass and relative coords
    call get_relative_coords(coords, center_of_mass, relative_coords)

    write(info_str,'(2X,A)') trim(msg)
    call localorb_info(info_str)

    conversion = hartree / bohr
    ! get net force before filtering
    call get_net_force(total_forces, net_force)
    write(info_str,'(2X,A,3(1X, E13.6),A)') &
    "| Net force on center of mass : ", net_force(:) * conversion, " eV/A"
    call localorb_info(info_str)

    if (n_periodic.eq.0) then
       conversion = hartree
       ! get net torque on center of mass before filtering
       call get_net_torque(total_forces, relative_coords, net_torque)
       write(info_str,'(2X,A,3(1X, E13.6),A)') &
       "| Net torque on center of mass: ", net_torque(:) * conversion, " eV"
       call localorb_info(info_str)
    end if

  end subroutine output_unitary_components
  !******
  !----------------------------------------------------------------------------
  !****s* relaxation/output_unitary_components_lv
  !  NAME
  !    output_unitary_components_lv
  !  SYNOPSIS

  subroutine output_unitary_components_lv (msg, forces_lv)

    !  PURPOSE
    !    Output magnitude of unitary components.
    !  USES

    implicit none

    !  ARGUMENTS

    character*(*), intent(IN) :: msg
    real*8, intent(IN) :: forces_lv(3, n_periodic)

    !  INPUTS
    !    o msg -- Introductory message
    !  OUTPUT
    !    o total_force -- forces to check
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: center_of_mass_lv(3), relative_lv(3, n_periodic)
    real*8 :: net_force_lv(3), net_torque_lv(3)
    real*8 :: conversion
    character*100 :: info_str

    ! begin work

    ! first get center of mass and relative coords
    call get_relative_lv(lattice_vector, center_of_mass_lv, relative_lv)

    write(info_str,'(2X,A)') trim(msg)
    call localorb_info(info_str)

    conversion = hartree / bohr
    ! get net force before filtering

    ! VA: not for lattice vectors

    !call get_net_force_lv(forces_lv, net_force_lv)
    !write(info_str,'(2X,A,3(1X, E12.6),A)') &
    !"| Net force on lattice vector center of mass : ", net_force_lv(:) * conversion, " eV/A"
    !call localorb_info(info_str)

       conversion = hartree
       ! get net torque on center of mass before filtering
       call get_net_torque_lv(forces_lv, relative_lv, net_torque_lv)
       write(info_str,'(2X,A,3(1X, E13.6),A)') &
       "| Net torque relative to origin: ", net_torque_lv(:) * conversion, " eV"
       call localorb_info(info_str)

  end subroutine output_unitary_components_lv
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/clean_mobile_force_components
  !  NAME
  !    clean_mobile_force_components
  !  SYNOPSIS

  subroutine clean_mobile_force_components(total_mobile_forces, mobile_ref_coords)
    !  PURPOSE
    !    Wrapper around clean_force_components().
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(n_mobile_coords), intent(INOUT) :: total_mobile_forces
    real*8, dimension(n_mobile_coords), intent(IN) :: mobile_ref_coords

    !  INPUTS
    !    o total_mobile_forces -- total forces before cleaning
    !    o mobile_ref_coords -- reference geometry
    !  OUTPUTS
    !    o total_mobile_forces -- total forces after cleaning
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: i_mobile_coord, i_atom, i, n_dims
    real*8  :: total_forces(3, n_atoms)
    real*8  :: forces_lv(3, n_periodic)
    real*8  :: ref_coords(3, n_atoms)
    real*8  :: ref_lv(3, n_periodic)

    total_forces = 0.d0
    forces_lv    = 0.d0
    ref_coords   = coords
    ref_lv       = 0.d0

    ! lattice_vector not allocated in non-periodic systems, must never touch it!
    if (n_periodic > 0) then
       ref_lv = lattice_vector
    end if

    call X2coords(ref_coords, ref_lv, mobile_ref_coords)
    call X2coords(total_forces, forces_lv, total_mobile_forces)

    call clean_force_components_work(ref_coords, ref_lv, total_forces, forces_lv)

    call coords2X(total_forces, forces_lv, total_mobile_forces)

  end subroutine clean_mobile_force_components
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/clean_force_components_work
  !  NAME
  !    clean_force_components_work
  !  SYNOPSIS

  subroutine clean_force_components_work(coords, lv, total_forces, forces_lv)

    !  PURPOSE
    !    central coordination of the cleaning of translations and
    !    rotations from forces
    !  USES
    !
    !  ARGUMENTS
    !  INPUTS
    !    o coords -- Coordinates of atoms
    !    o total_forces - input forces straight from scf_solver
    !    o lv           - lattice vector
    !    o forces_lv    - input forces on lattice vectors
    !  OUTPUT
    !    o total_force - cleaned forces
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    implicit none

    character*100 :: info_str
    real*8, dimension(3,n_atoms),    intent(IN)    :: coords
    real*8, dimension(3,n_periodic), intent(IN)    :: lv
    real*8, dimension(3,n_atoms),    intent(INOUT) :: total_forces
    real*8, dimension(3,n_periodic), intent(INOUT) :: forces_lv

    select case (remove_unitary_force_components)

    case (0)

       ! only enforce relaxation constraint by setting forces to zero on those atoms
       if (use_relaxation_constraints) then
          call enforce_relaxation_constraints (total_forces, forces_lv)
       end if

    case (2)

       ! in this case, translation and rotation cleaning is done, but only
       ! if less than three atoms are constrained
       ! Relaxation constraints are enforced explicitly thereafter.
       ! MOL: Also clean if symmetry-constrained relaxation (has to be called with
       ! with symmetrized forces though).

       if ((.not. use_geo_relaxation) .or. (n_mobile_atom_coords.gt.(3*n_atoms-6)) .or. use_symm_const_geo_relaxation ) then
          call remove_translation_and_rotation(coords, total_forces)
       end if

       if (use_relaxation_constraints) then
          call enforce_relaxation_constraints (total_forces, forces_lv)
       end if

    case default
       write (info_str,'(1X,A,A)') "* Removal of rotations and translations: ", &
       " not set - internal inconsistency."
       call localorb_info(info_str,use_unit,'(A)')
       stop

    end select

  end subroutine clean_force_components_work
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/remove_translation_and_rotation
  !  NAME
  !    remove_translation_and_rotation
  !  SYNOPSIS

  subroutine remove_translation_and_rotation (coords, total_forces)

    !  PURPOSE
    !    removes noise in the forces that cause unphysical translation or rotation of a cluster
    !    by fulfilling the sayvetz-conditions (net force vanishes, net torque on center of mass vanishes)
    !      (see http://www.gaussian.com/g_whitepap/vib.htm)
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3, n_atoms)    :: coords
    real*8, dimension(3, n_atoms)    :: total_forces

    !  INPUTS
    !    o coords -- coordanates of atoms
    !    o total_forces -- total forces before cleaning
    !  OUTPUTS
    !    o total_forces -- total forces after cleaning
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! local variables
    real*8, dimension(3,n_atoms) :: relative_coords
    real*8 :: center_of_mass(3)
    integer :: i_coord_p1
    integer :: i_coord_p2
    real*8 :: inv_norm
    real*8 :: translation_component
    real*8 :: rotation_component
    real*8 :: norm

    ! functions
    real*8 :: ddot
    real*8 :: dnrm2

    ! counter
    integer :: i_atom
    integer :: i_coord
    integer :: i_coord_2

    ! first get center of mass and relative coords
    call get_relative_coords(coords, center_of_mass, relative_coords)

    if (n_periodic.eq.0) then
    !VA in the cluster case take out the whole rotation
       ! now, get moments of inertia
       tensor_of_inertia = 0.d0
       do i_atom = 1, n_atoms, 1
          do i_coord = 1, 3, 1
             i_coord_p1 = mod(i_coord , 3)  + 1
             i_coord_p2 = mod(i_coord + 1, 3) + 1
             ! diagonal elements
             tensor_of_inertia(i_coord, i_coord) = tensor_of_inertia(i_coord, i_coord) &
                  + (relative_coords(i_coord_p1, i_atom) * relative_coords(i_coord_p1, i_atom)) + &
                  (relative_coords(i_coord_p2, i_atom) * relative_coords(i_coord_p2, i_atom))
             ! non-diagonal elements
             do i_coord_2 = i_coord + 1, 3, 1
                tensor_of_inertia(i_coord, i_coord_2) = tensor_of_inertia(i_coord, i_coord_2) &
                     - relative_coords(i_coord, i_atom) * relative_coords(i_coord_2, i_atom)
             end do
          end do
       end do

       ! diagonalize it
       call dsyev('V', 'U', 3, tensor_of_inertia, 3, moments_of_inertia, work_cf, lwork_cf, info)

       ! generate vectors of rotation according to gaussian
       do i_coord = 1, 3, 1
          i_coord_p1 = mod(i_coord, 3) + 1
          i_coord_p2 = mod(i_coord + 1, 3) + 1
          do i_atom = 1, n_atoms, 1
             do i_coord_2 = 1, 3, 1
                rotation_vectors(i_coord_2, i_atom, i_coord) = &
                     (ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:,i_coord_p1), 1) * &
                     tensor_of_inertia(i_coord_2, i_coord_p2) - &
                     ddot(3, relative_coords(:,i_atom), 1, tensor_of_inertia(:, i_coord_p2), 1) * &
                     tensor_of_inertia(i_coord_2, i_coord_p1))
                if (i_coord .eq. 2) then
                   rotation_vectors(i_coord_2, i_atom, i_coord) = - rotation_vectors(i_coord_2, i_atom, i_coord)
                end if
             end do
          end do
          ! normalize them
          norm = dnrm2(3*n_atoms, rotation_vectors(:, :, i_coord), 1)
          if (norm .gt. 0) then
             inv_norm = 1.d0 / norm
             rotation_vectors(:, :, i_coord) = rotation_vectors(:, :, i_coord) * inv_norm
             ! project out rotation
             rotation_component = ddot(3*n_atoms, rotation_vectors(1,1,i_coord), 1, total_forces, 1)
             total_forces(:,:) = total_forces(:,:) - rotation_component * rotation_vectors(:,:,i_coord)
          end if
       end do
    end if

    ! project out translation on atoms
    do i_coord = 1, 3, 1
       translation_component = ddot(3*n_atoms, translation_vectors(1,1,i_coord), 1, total_forces, 1)
       total_forces(:,:) = total_forces(:,:) - translation_component * translation_vectors(:,:,i_coord)
    end do

  end subroutine remove_translation_and_rotation
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/remove_translation_and_rotation
  !  NAME
  !    remove_translation_and_rotation
  !  SYNOPSIS

  subroutine remove_translation_and_rotation_lv (lv, forces_lv)

    !  PURPOSE
    !    removes noise in the forces that cause unphysical translation or rotation of a cluster
    !    by fulfilling the sayvetz-conditions (net force vanishes, net torque on center of mass vanishes)
    !      (see http://www.gaussian.com/g_whitepap/vib.htm)
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3, n_periodic) :: lv
    real*8, dimension(3, n_periodic) :: forces_lv

    !  INPUTS
    !    o coords -- coordanates of atoms
    !    o total_forces -- total forces before cleaning
    !  OUTPUTS
    !    o total_forces -- total forces after cleaning
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! local variables
    real*8, dimension(3,n_periodic) :: relative_lv
    real*8 :: center_of_mass_lv(3)
    integer :: i_coord_p1
    integer :: i_coord_p2
    real*8 :: inv_norm
    real*8 :: translation_component_lv
    real*8 :: rotation_component_lv
    real*8 :: norm

    ! functions
    real*8 :: ddot
    real*8 :: dnrm2

    ! counter
    integer :: i_periodic
    integer :: i_coord
    integer :: i_coord_2

    ! first get center of mass and relative coords
    ! Do nothing - everything in origin!
    call get_relative_lv(lv, center_of_mass_lv, relative_lv)

    if ((n_periodic.ne.0) .and. (use_numerical_stress)) then
    !VA: in the periodic case take out rotations of lattice vectors only!
    !VA: currently forces on lattice vectors are only used when geometry is relaxed.
    !    this might change later when an analytical stress is available

       ! now, get moments of inertia
       tensor_of_inertia_lv = 0.d0
       do i_periodic = 1, n_periodic, 1
          do i_coord = 1, 3, 1
             i_coord_p1 = mod(i_coord , 3)  + 1
             i_coord_p2 = mod(i_coord + 1, 3) + 1
             ! diagonal elements
             tensor_of_inertia_lv(i_coord, i_coord) = tensor_of_inertia_lv(i_coord, i_coord) &
                  + (relative_lv(i_coord_p1, i_periodic) * relative_lv(i_coord_p1, i_periodic)) + &
                  (relative_lv(i_coord_p2, i_periodic) * relative_lv(i_coord_p2, i_periodic))
             ! non-diagonal elements
             do i_coord_2 = i_coord + 1, 3, 1
                tensor_of_inertia_lv(i_coord, i_coord_2) = tensor_of_inertia_lv(i_coord, i_coord_2) &
                     - relative_lv(i_coord, i_periodic) * relative_lv(i_coord_2, i_periodic)
             end do
          end do
       end do

       ! diagonalize it
       call dsyev('V', 'U', 3, tensor_of_inertia_lv, 3, moments_of_inertia_lv, work_cf, lwork_cf, info)

       ! generate vectors of rotation according to gaussian
       do i_coord = 1, 3, 1
          i_coord_p1 = mod(i_coord, 3) + 1
          i_coord_p2 = mod(i_coord + 1, 3) + 1
          do i_periodic = 1, n_periodic, 1
             do i_coord_2 = 1, 3, 1
                rotation_vectors_lv(i_coord_2, i_periodic, i_coord) = &
                     (ddot(3, relative_lv(:,i_periodic), 1, tensor_of_inertia_lv(:,i_coord_p1), 1) * &
                     tensor_of_inertia_lv(i_coord_2, i_coord_p2) - &
                     ddot(3, relative_lv(:,i_periodic), 1, tensor_of_inertia_lv(:, i_coord_p2), 1) * &
                     tensor_of_inertia_lv(i_coord_2, i_coord_p1))
                if (i_coord .eq. 2) then
                   rotation_vectors_lv(i_coord_2, i_periodic, i_coord) = - rotation_vectors_lv(i_coord_2, i_periodic, i_coord)
                end if
             end do
          end do
          ! normalize them
          norm = dnrm2(3*n_periodic, rotation_vectors_lv(:, :, i_coord), 1)
          if (norm .gt. 0) then
             inv_norm = 1.d0 / norm
             rotation_vectors_lv(:, :, i_coord) = rotation_vectors_lv(:, :, i_coord) * inv_norm
             ! project out rotation
             rotation_component_lv = ddot(3*n_periodic, rotation_vectors_lv(1,1,i_coord), 1, forces_lv, 1)
             forces_lv(:,:)     = forces_lv(:,:) - rotation_component_lv * rotation_vectors_lv(:,:,i_coord)
          end if
       end do
    end if

  end subroutine remove_translation_and_rotation_lv
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/enforce_relaxation_constraints
  !  NAME
  !    enforce_relaxation_constraints
  !  SYNOPSIS

  subroutine enforce_relaxation_constraints ( total_forces, forces_lv )

    !  PURPOSE
    !   Deletes the forces on fixed atoms
    !   Trivial way to enforce relaxation constraints: Simply set to zero the
    !   forces on atoms whose position shall not be moved.
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms) :: total_forces
    real*8, dimension(3,n_periodic) :: forces_lv

    !  INPUTS
    !    o total_forces -- total forces
    !  OUTPUTS
    !   o total_forces -- total forces after constrained forces are put to 0
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! local variables

    integer :: i_atom, i_lv

    ! begin work

    do i_atom = 1, n_atoms, 1
       if (constrain_relaxation(i_atom)) then
          total_forces(1:3,i_atom) = 0.d0
       else
         if (constrain_coords(1,i_atom)) then
           total_forces(1,i_atom) = 0.d0
         end if
         if (constrain_coords(2,i_atom)) then
           total_forces(2,i_atom) = 0.d0
         end if
         if (constrain_coords(3,i_atom)) then
           total_forces(3,i_atom) = 0.d0
         end if
       end if
    end do

    ! lattice vectors -- only if unit cell is relaxed
    if (relax_unit_cell .ne. 0)  then
       do i_lv = 1, n_periodic, 1
          if (constrain_lv_relaxation(i_lv)) then
             forces_lv(1:3,i_lv) = 0.d0
          else
            if (constrain_lv_components(1,i_lv)) then
              forces_lv(1,i_lv) = 0.d0
            end if
            if (constrain_lv_components(2,i_lv)) then
              forces_lv(2,i_lv) = 0.d0
            end if
            if (constrain_lv_components(3,i_lv)) then
             forces_lv(3,i_lv) = 0.d0
            end if
          end if
       end do
    end if

  end subroutine enforce_relaxation_constraints
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/check_maximum_displacement
  !  NAME
  !    check_maximum_displacement
  !  SYNOPSIS

  subroutine check_maximum_displacement(forces)

    !  PURPOSE
    !    simple fix to prevent clusters from exploding
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms) :: forces

    !  INPUTS
    !   o forces -- atomic forces
    !  OUTPUT
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    ! local variables
    real*8 :: max_force_component
    character*80 :: info_str

    ! counter
    integer :: i_atom
    integer :: i_coords

    max_force_component = 0.d0
    do i_atom = 1, n_atoms, 1
       do i_coords = 1, 3, 1
          max_force_component = max(max_force_component, abs(forces(i_coords, i_atom)))
       end do
    end do
    if ((relaxation_step * relaxation_step).gt.(maximum_displacement / max_force_component)) then
       relaxation_step = sqrt(maximum_displacement / max_force_component)
       limited_displacement = .true.
       write (info_str,'(2X,A,E14.6)') "Displacement of atoms too large. Decreasing relaxation step to ", &
            relaxation_step
       call localorb_info(info_str)
    else
       limited_displacement = .false.
    end if

  end subroutine check_maximum_displacement
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/get_next_relaxation_step
  !  NAME
  !   get_next_relaxation_step
  !  SYNOPSIS

  subroutine get_next_relaxation_step   &
            ( free_energy,              &
              coords,                   &
              total_forces,             &
              lv,                       &
              forces_lv,                &
              periodic_unit_cell_trans, &
              KS_eigenvector,           &
              KS_eigenvector_complex,   &
              occ_numbers,              &
              valid_geo )

    !  PURPOSE
    !    Wrapper around the relaxation routines
    !
    !  USES
    use timing_core, only: get_timestamps, get_times
    implicit none
    !  ARGUMENTS
    real*8 ::  free_energy
    real*8, dimension(3,n_atoms), intent(INOUT) :: coords
    real*8, dimension(3,n_atoms) :: total_forces
    real*8, dimension(3,n_periodic), intent(INOUT) :: lv
    real*8, dimension(3,n_periodic) :: forces_lv
    real*8, dimension(3,n_atoms), intent(INOUT) :: periodic_unit_cell_trans
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task)     :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points_task) :: occ_numbers
    logical :: valid_geo

    ! INPUTS
    ! o energy -- total energy
    ! o coords -- coordinates
    ! o lv     -- lattice vectors
    ! o total_forces -- forces on atoms
    ! o forces_lv    -- generalized 'forces' on lattice vectors (Atoms still move)
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! OUTPUT
    ! o coords -- updated coordinates
    ! o lv     -- updated lattice vectors
    ! o valid_geo -- geometry checking flag
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2008).
    ! SOURCE
    real*8 :: t_relax
    real*8 :: t_tmp

    call get_timestamps(t_relax,t_tmp)

    if (relax_geo.eq.RELAX_BFGS_TB) then

      call bfgs_descent                   &
               (free_energy,              &
                coords,                   &
                total_forces,             &
                lv,                       &
                forces_lv,                &
                periodic_unit_cell_trans, &
                KS_eigenvector,           &
                KS_eigenvector_complex,   &
                occ_numbers,              &
                valid_geo )

    else if (relax_geo.eq.RELAX_TRM .or. relax_geo.eq.RELAX_LTRM) then
       call trusted_descent               &
               (free_energy,              &
                coords,                   &
                total_forces,             &
                lv,                       &
                forces_lv,                &
                periodic_unit_cell_trans, &
                KS_eigenvector,           &
                KS_eigenvector_complex,   &
                occ_numbers,              &
                valid_geo)
    end if

    call get_times(t_relax,t_tmp)

    write(info_str,'(2X,A)') "Finished advancing geometry"
    call localorb_info(info_str,use_unit)
    write(info_str,'(2X,A,F10.3,A)') "| Time : ",t_relax," s"
    call localorb_info(info_str,use_unit)

  end subroutine get_next_relaxation_step
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/extract_lv_forces
  !  NAME
  !    extract_lv_forces
  !  SYNOPSIS

  subroutine extract_lv_forces(forces_lv, stress_tensor, atomic_forces)

    !  PURPOSE
    !  This routine calculates and cleans the forces on the
    !  lattice vectors obtained from the stress tensor.
    !
    ! STEP 1:
    !       Get  generalized derivativesi/forces on the lattice vectors,
    !       where all atoms are still allowed to move.
    !
    !
    ! STEP 2:
    !       Clean these generalized forced from it's atomic contributions.
    !       In index notation this matrix operation reads:
    !
    !  [F_lattice]_ml = [F_stress]_ml - \sum_n{  [(Fractional_coord)^T]_mn*[F_atom]_nl }
    !
    !  [F_atom]_nl            : l-th component of the force on atom n
    !  [F_lattice]_ml         : l-th component of the force on lattice vector m
    !  [Fractional_coord]_nm  : m-th component of atom n in fractional charges
    !  (Fractional_coord)^T   : transposed of Fractional_coord
    !   F_stress              : 'generalized forces' on lattice vectors derived from the stress tensor
    !                           -> atoms still 'move'
    !
    !  The transposed of the fractional cooredinates (Fractional_coord)^T can
    !  can be obtained by:
    !
    !  (Fractional_coord) = (coord)*(L)^-1  => (Fractional_coord)^T = (L)^-1T * (coord)^T
    !  .... with a atom-vector (Fractional_coord)/ (coord) of a given atom-index n.
    !
    !
    !  NOTE: these formulars refer to the above defined index pictures. within aims typically the
    !        second index refers to the atom number -- in contrast to these formulars!
    !
    !  USES
    implicit none
    !  ARGUMENTS
    real*8, dimension(3,n_periodic), intent(OUT) :: forces_lv
    real*8, dimension(3,n_periodic), intent(IN)  :: stress_tensor
    real*8, dimension(3,n_atoms),    intent(IN)  :: atomic_forces
    ! INPUTS
    !
    ! OUTPUT
    !
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2008).
    ! SOURCE

    ! local variables
    real*8, dimension(3,n_atoms) :: effective_atomic_forces
    real*8, dimension(3,n_atoms) :: transp_fractional_coords
    real*8, dimension(3,3) :: inverse_lv
    real*8, dimension(3):: work
    integer,dimension(3):: ipivot
    integer :: i_atom, j, i


    ! FlK: In lattice_trm, the
    ! lattice contribution to atomic forces is taken to account, do nothing here

    ! 0. Get from the stress tensor the generalized forces
      call get_derivatives_on_lattice (stress_tensor, forces_lv)
      forces_lv=-forces_lv

    ! fkdev
    if ((relax_geo.eq.RELAX_TRM) .and. (.not. use_symm_const_geo_relaxation)) then

      !Get index picture of the above formulars
      forces_lv= transpose(forces_lv)

      ! 1. Obtain  inverse of lattice vector matrix
      inverse_lv = lattice_vector
      call DGETRF(3, 3, inverse_lv, 3, ipivot, info )
      if (info /= 0) then
         write(use_unit,*) 'ERROR: new lattice_vector is singular!'
         stop
      end if
      call DGETRI(3, inverse_lv, 3, ipivot, work, 3, info)
      if (info /= 0) then
         write(use_unit,*) 'ERROR: new lattice_vector is singular!'
         stop
      end if

      ! 2. Obtain transposed of fractional coord matrix
      do i_atom = 1, n_atoms, 1
        transp_fractional_coords(:,i_atom) = matmul(inverse_lv, coords(:,i_atom))
      end do

      ! 3. Substract atomic part from the
      !    stress-derived forces on the lattice vectors

      effective_atomic_forces =  atomic_forces

      ! VA: for FRACTIONALLY constrained atoms, don't take into account their contribution
      !     Currently this does not have any effect, because fixed force contributions
      !     are set to 0 (do we really always want to have this?)!
      !     Currently this is done for whole atom -> the check for that is in read_geo.f90
      !     However here we still asign the forces component wise to zero
      if (use_relaxation_constraints) then
         do i_atom = 1, n_atoms,1
            do i= 1, 3, 1
               if ((coord_basis_atoms(i_atom) == 1) .and. &   ! frac
                  (constrain_coords(i,i_atom))) then
                   effective_atomic_forces (i, i_atom) = 0d0
               end if
            end do
         end do
      end if

      forces_lv = forces_lv -&
      matmul (transp_fractional_coords, transpose (effective_atomic_forces))

      ! 4. get back aims-index picture
      forces_lv= transpose(forces_lv)
    end if
  end subroutine extract_lv_forces

  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/project_lv_forces_4_shape_relaxation
  !  NAME
  !    project_lv_forces_4_shape_relaxation
  !  SYNOPSIS

  subroutine project_lv_forces_4_shape_relaxation(forces_lv)

  !  PURPOSE
  !
  !  USES
  implicit none
  !  ARGUMENTS
  real*8, dimension(3,n_periodic), intent(INOUT) :: forces_lv
  ! INPUTS
  !
  ! OUTPUT
  !
  ! AUTHOR
  !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  ! HISTORY
  !   Release version, FHI-aims (2008).
  ! SOURCE

  ! local variables
  integer                         ::  i, i_periodic
  real*8                          ::  norm_squared, force_dot_lv
  real*8, dimension(3,n_periodic) ::  forces_lv_proj

    do i_periodic = 1, n_periodic, 1
       norm_squared  = sum(lattice_vector(:,i_periodic)**2)
       force_dot_lv  = dot_product (forces_lv(:,i_periodic), lattice_vector(:,i_periodic))
       forces_lv_proj (:,i_periodic) = lattice_vector(:,i_periodic)*force_dot_lv/norm_squared
    end do

    forces_lv = forces_lv_proj

  end subroutine project_lv_forces_4_shape_relaxation
  !******
  !---------------------------------------------------------------------------------------------------------------------------------
  !**** NAME get_derivatives_on_lattice
  !          get_derivatives_on_lattice
  !  SYNOPSIS

  subroutine get_derivatives_on_lattice(stress_tensor, energy_gradient_on_lattice)

  !  PURPOSE
  !  Calculates the derivative of the (free) energy with respect
  !  to the lattice vectors once the stress tensor is given.
  !  In essence this is done by a simple matrix multiplication
  !
  !  Mind the trap!:
  !  According to [doll] we have (read twice and mind the index picture)
  !
  !         V*Stress = L^T dE/dL
  !
  !  with:
  !        (L)_nl       the matrix of lattice vector, with lattice vector number n and component l
  !         Stress      Stress tensor
  !         V           Unit cell volume
  !         E           Free energy
  !
  ! Since the index picture in aims is different: L_aims = L^T we conclude:
  !
  !        V*((L_aims)^-1 * Stress)^T = dE/dL_aims
  !
  !

  implicit none

  real*8, dimension(3,3), intent(in)        :: stress_tensor
  real*8, dimension(3,3), intent(out)       :: energy_gradient_on_lattice
  integer                                   :: i_latt, i_atom
  real*8                                    :: inv_lattice_vector(3,3)
  real*8, dimension(3)                      :: work
  integer,dimension(3)                      :: ipivot
  integer                                   :: info, i, j

  !  INPUTS
  !    none
  !  OUTPUT
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !  SOURCE

    energy_gradient_on_lattice(:,:)=0d0

    ! Get inverse matrix
    do i_latt = 1, n_periodic
       inv_lattice_vector(:,i_latt) = lattice_vector(:,i_latt)
    end do

    call DGETRF(3, 3, inv_lattice_vector, 3, ipivot, info )
    if (info /= 0) then
       write(use_unit,*) 'ERROR: lattice_vector is singular!'
       stop
    end if
    call DGETRI(3, inv_lattice_vector, 3, ipivot, work, 3, info)
    if (info /= 0) then
       write(use_unit,*) 'ERROR: lattice_vector is singular!'
       stop
    end if

    ! Get gradient of the energy with respect to lattice vectors
    energy_gradient_on_lattice  = matmul (inv_lattice_vector,stress_tensor)
    energy_gradient_on_lattice  =  energy_gradient_on_lattice * cell_volume

    ! Transpose on order to have the same index picture as lattice_vector
    energy_gradient_on_lattice  = transpose (energy_gradient_on_lattice)

  end subroutine get_derivatives_on_lattice
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/bfgs_descent
  !  NAME
  !    bfgs_descent
  !  SYNOPSIS

  subroutine bfgs_descent &
                    (energy, coords, total_forces, lv, &
                     forces_lv, periodic_unit_cell_trans, &
                     KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers, valid_geo)
    !  PURPOSE
    !    Driver routine for the BFGS-based relaxation - it contains the main decision tree.
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8 :: energy
    real*8, dimension(3,n_atoms) :: coords
    real*8, dimension(3,n_atoms) :: frac_coords_temp
    real*8, dimension(3,n_atoms) :: total_forces
    real*8, dimension(3,n_periodic) :: lv
    real*8, dimension(3,n_periodic) :: forces_lv
    real*8, dimension(3,n_atoms) :: periodic_unit_cell_trans

    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task)     :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points_task)     :: occ_numbers
    logical :: valid_geo

    ! INPUTS
    ! o energy -- total energy
    ! o coords -- coordinates
    ! o lv     -- lattice vectors
    ! o total_forces -- forces on atoms
    ! o forces_lv    -- generalized 'forces' on lattice vectors
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! OUTPUT
    ! o coords -- updated coordinates
    ! o lv     -- updated lattice vectors
    ! o valid_geo -- geometry checking flag
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2008).
    ! SOURCE

    ! Local variables

    real*8 :: force_times_searchdir
    real*8 :: line_step_predicted, line_step_change

    real*8, external :: ddot

    character(*), parameter :: func = 'bfgs_descent'

    ! counters
    integer :: i_mobile_coord, i_atom, jj
    logical :: line_search_accepted, step_accepted, negative_curvature_in_searchdir

    ! begin work

    call localorb_info("Advancing geometry using BFGS.",use_unit,'(2X,A)' )

    ! Determine whether relaxation step was successful, current
    ! coordinates, Hessian, search direction, and line step
    if (ini_relaxation) then
       ! This is the initial relaxation step - must store energy, forces etc.
       ! Store energy, coordinates, forces of present iteration for next iteration
       stored_energy       = energy
       stored_coords       = coords
       stored_lv           = lv
       stored_forces       = total_forces
       stored_forces_lv    = forces_lv
       stored_periodic_unit_cell_trans = periodic_unit_cell_trans

       line_search_on      = .false.
       step_was_restricted = .false.
       negative_curvature_in_searchdir = .false.
       t_extrap = .false.

       ! allocate data for storing the KS eigenvectors during relaxation ...
       call bfgs_allocate_eigenvectors()

       ! memorize the eigenvectors for the next step
       call bfgs_store_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

       ! first relaxation step
       call update_searchdir(hessian, search_direction, total_forces, forces_lv)

       ! remember direction*force for line search
       call get_force_times_searchdir_hat &
            (total_forces, forces_lv, search_direction, stored_force_times_searchdir_hat)

    else
       ! not initial relaxation:

       ! calculate inner product of current forces and previous search direction ...
       call get_force_times_searchdir_hat &
            (total_forces, forces_lv, search_direction, force_times_searchdir_hat)

       ! ... to predict the line step that would have been optimal in retrospect and then see if it did the trick ...

       if (abs((stored_force_times_searchdir_hat - force_times_searchdir_hat)/force_times_searchdir_hat).gt.1d-15) then
          ! prevent division by zero ...
          line_step_predicted = line_step * stored_force_times_searchdir_hat / &
               (stored_force_times_searchdir_hat - force_times_searchdir_hat)
       else
          write(info_str,'(2X,A)') 'The energy landscape appears LINEAR in the search direction. '
          call localorb_info(info_str)
          write(info_str,'(2X,A)') 'This seems higly improbable, make sure that the geometry is correct ... '
          call localorb_info(info_str)
          line_step_predicted = line_step
       end if

       ! what happens if line_step_predicted < 0, i.e. the curvature in search_dir is negative? Then the line search
       ! would eventually converge to a saddle point ... We definitely don't want that and require special treatment
       ! in this particular case:
       negative_curvature_in_searchdir = (line_step_predicted.lt.0d0)

       ! check whether they the current line search seems to be alright:
       ! EITHER if the (a posteriori) predicted line step is within a factor line_search_tol of the actual line_search
       !        AND the hessian had not been reset on the previous step (where BFGS knows nothing about curvature ... )
       ! OR     if we were already doing a line search AND the hessian is new - in that case there is nothing more to be
       !        done than updating the hessian using all the available information for a new hessian
       ! AND    - curvature is positive, otherwise line search = stupid idea
       !        - the energy went down - otherwise there's no point for a line search
       !
       ! Note: a line search acceptance does not guarantee that the actual step is accepted as well, see below

       line_search_accepted = ( &
            .not. ((abs(line_step_predicted/line_step-1d0).ge.line_search_tol) &
            .and. (.not. t_first_pred)) .or. &
            ((line_search_on) .and. initial_hessian))       &
            .and. (.not. negative_curvature_in_searchdir)   &
            .and. (energy.lt.stored_energy)                &
            .or. negative_curvature_in_searchdir

       ! status output for this thing:
       write(info_str,'(2X,A,E14.6)') 'Predicted optimal line step for last step: ', line_step_predicted
       call localorb_info(info_str)
       write(info_str,'(2X,A,E14.6)') 'Line step actually used for last step:     ', line_step
       call localorb_info(info_str)
       if (negative_curvature_in_searchdir) then
          write(info_str,'(2X,A)') 'It looks like the energy surface in the search direction has negative curvature!'
          call localorb_info(info_str)
          write(info_str,'(2X,A)') 'Accepting line search and continuing with maximally permissible line step to get out of this.'
          call localorb_info(info_str)
       end if

       ! check whether or not a single step is accepted:

       ! (1) Simplest case: energy went down and line search was appreciated:
       step_accepted = line_search_accepted

       ! (2) if there is a small energy increase and the hessian is initial and line_step = small, then it might
       !     be useful to accept toward the end of a run:
       step_accepted = step_accepted .or. ( &
            ((energy - stored_energy).lt.energy_tolerance) .and. &
            initial_hessian .and. (line_step.lt.min_line_step))

       ! (3) accept step if the curvature was negative in search_direction AND the energy went down
       !     as the other treatments fail in this situation
       step_accepted = step_accepted .or. &
            ((energy.lt.stored_energy) .and. negative_curvature_in_searchdir)

       ! (4) ensure backward compatibility: if not line_search_automatic, one should still accept the
       !     step anyway if the energy went down!
       step_accepted = step_accepted .or. &
            ((energy.lt.stored_energy) .and. (.not. line_step_reduce_automatic))

       ! (5) if the line step was restricted due to too large displacements, then also accept if energy went down
       step_accepted = step_accepted .or. (step_was_restricted .and. (energy.lt.stored_energy))

       ! (5) if there was a line search AND it failed so far, then one should at least accept either the
       !     initially rejected step, OR the line_search step, whichever was better.
       !     Otherwise, it might be that there are three subsequent single-point calculations in one search direction,
       !     which is NOT what we want in a single-point line search
       if (line_search_on .and. (.not. step_accepted)) then

          ! check on the energies of the two last steps respectively, and then accept the one that
          ! was best
          if (energy.gt.stored_energy_2) then
             energy       = stored_energy_2
             coords       = stored_coords_2
             lv           = stored_lv_2
             total_forces = stored_forces_2
             forces_lv    = stored_forces_lv_2
             periodic_unit_cell_trans = stored_periodic_unit_cell_trans_2
             line_step    = stored_line_step
             write(info_str,'(2X,A)') "Line search failed, the initial guess was better. "
             call localorb_info(info_str)
             write(info_str,'(2X,A)') "Continuing relaxation with BFGS update from initial step. "
             call localorb_info(info_str)
          else
             ! the second, line-searched, step was better but not ideal. Let's keep it anyway
             ! and notify the user of that occurrence; all things remain as they were before
             write(info_str,'(2X,A)') "Line search was far from good, but better than initial guess. "
             call localorb_info(info_str)
             write(info_str,'(2X,A)') "Accepting current step for that reason and continuing BFGS relaxation. "
             call localorb_info(info_str)
          end if

          step_accepted   = .true.
       end if

       ! check whether to continue with bfgs ... ?
       if (step_accepted) then

          ! remember eigenvectors and -values for next step in case the current one was not accepted!
          call bfgs_store_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

          ! update hessian matrix and search direction using BFGS from wikipedia
          call update_hessian_and_searchdir_BFGS(hessian, search_direction, total_forces, forces_lv)
          ! JW: should be tested before uncommenting
          ! call clean_mobile_force_components(search_direction, mob_coords)

          ! Store energy, coordinates, forces of present iteration for next iteration

          stored_energy     = energy
          stored_forces     = total_forces
          stored_forces_lv  = forces_lv
          stored_coords     = coords
          stored_lv         = lv
          stored_periodic_unit_cell_trans = periodic_unit_cell_trans

          ! check whether the search direction makes any sense compared to the current forces
          ! scalar product: force_difference * search_direction_unit_vector
          call get_force_times_searchdir_hat &
               (total_forces, forces_lv, search_direction, force_times_searchdir)

          ! if the forces point away from the current search direction, the Hessian can't be any good
          ! in that case, reinitialize the Hessian and restart with a steepest descent anyway
          if (force_times_searchdir.lt.0.d0) then

             write(info_str,'(2X,A,E14.6)') &
                  "BFGS: Forces * search direction is negative: ", force_times_searchdir
             call localorb_info( info_str )
             write(info_str,'(2X,A)') &
                  "Search direction should not point away from forces - reinitializing Hessian."
             call localorb_info( info_str )

             call update_hessian_and_searchdir_steepest_descent &
                  (hessian, search_direction, total_forces, forces_lv)

          end if

          ! remember for next time:
          line_search_on = .false.
          call get_force_times_searchdir_hat &
               (total_forces, forces_lv, search_direction, stored_force_times_searchdir_hat)

          ! If the last step was extrapolated, reset line_step to 1
          if (t_extrap) then
             line_step = 1.d0
             t_extrap = .false.
          end if

          ! attempt to do a predefined line_step mixing, possibly speeding up the process ?
          if (bfgs_extrapolation_on) then
             ! mix the calculated line step with the actually used line step: This makes sure that the experience with line
             ! searches from the previous few iterations is taken into account and that the line searches are longer (or shorter!)
             ! by default if they have always been that previously.
             line_step = line_step*(1d0-bfgs_mixing_factor) &
                  + min(line_step_predicted, bfgs_mixing_cap) * bfgs_mixing_factor
             write(info_str,'(2X,A,F15.8)') "BFGS: mixed line step with predicted value, now using : ", line_step
             call localorb_info( info_str )
          else
             ! quadratic extrapolation: line step 1 corresponds to the
             ! solution on a perfectly quadratic energy surface
             line_step = 1d0
          end if

       else
          ! Step was not accepted - revert to previous coordinates,

          ! keep information on the current situation in case the next step is rejected as well,
          ! then you can at least the less evil of two guesses and prevent infinite line searches in the
          ! same direction!
          stored_coords_2    = coords
          stored_lv_2        = lv
          stored_forces_2    = total_forces
          stored_forces_lv_2 = forces_lv
          stored_energy_2    = energy
          stored_periodic_unit_cell_trans_2 = periodic_unit_cell_trans

          ! for now: use same search direction as before but change line step according to preselected criteria
          coords       = stored_coords
          lv           = stored_lv
          total_forces = stored_forces
          forces_lv    = stored_forces_lv
          energy       = stored_energy
          periodic_unit_cell_trans = stored_periodic_unit_cell_trans

          ! also obtain scf solution from previous iteration ...
          call bfgs_load_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

          ! check why step was rejected in the first place  ...
          if ((.not. line_search_on) .and. (.not. line_search_accepted)) then         ! line search was not on & step was garbage, so search now ...

             write(info_str,'(2X,A)') &
                  "BFGS line minimization: Last trial step was not optimal, reverting to previous geometry."
             call localorb_info( info_str )

             if (line_step_reduce_automatic) then
                line_step_change = line_step_predicted/line_step
                ! keep original line step in case line search fails anyway
                stored_line_step = line_step

                if (.not. t_first_pred) then
                   line_step = line_step_predicted                ! use originally predicted line step from above
                else
                   if (line_step_predicted.lt.line_step) then
                      ! Do not extrapolated if the Hessian had just been reinitialized
                      line_step = line_step_predicted
                   end if
                end if

                ! Avoid too dangerous extrapolations:
                line_step = min(line_step,bfgs_extrapolation_cap)
                t_extrap = .true.

                write(info_str,'(2X,A,E14.6)') &
                     "Using a one-step line search with a line step of: ", line_step
                call localorb_info( info_str )

                ! this is a line search, remember for next time:
                line_search_on = .true.
             else
                line_step = line_step * line_step_reduce       ! for backward-compatibility: simply reduce the line step by constant factor
                write(info_str,'(2X,A,E14.6)') &
                     "Reducing line step to: ", line_step
                call localorb_info( info_str )
                line_search_on = .false.
             end if

          end if

          ! step was way too short, this looks like noise, restart calculation from steepest descent and
          ! reinitialize Hessian.
          if ((line_step.le.min_line_step) .and. (.not. initial_hessian)) then

             write(info_str,'(2X,A)') &
                  "BFGS: Extrapolation using current Hessian along current search direction"
             call localorb_info( info_str )
             write(info_str,'(2X,A)') &
                  "is not nearly quadratic - incorrect Hessian for current point?"
             call localorb_info( info_str )
             write(info_str,'(2X,A)') &
                  "Resetting search direction and Hessian to steepest descent."
             call localorb_info( info_str )

             call update_hessian_and_searchdir_steepest_descent &
                  (hessian, search_direction, total_forces, forces_lv)

             ! search direction changed, no line search:
             line_search_on = .false.
             call get_force_times_searchdir_hat  &
                  (total_forces, forces_lv, search_direction, stored_force_times_searchdir_hat)

          end if ! if ((line_step.le.min_line_step) .and. (.not. initial_hessian))
       end if    ! if (step_accepted)
    end if       ! if (ini_relaxation)

    ! now that we have a current search direction, make sure that the displacement remains below
    ! a threshold - Check for overly large movements and reduce line_step if needed
    ! this can only happen after the search direction was changed, i.e. when line_step.eq.1)
    ! can also happen when the line search is on

    step_was_restricted = .false.
    call restrict_max_displacement(search_direction, line_step, step_was_restricted, &
         negative_curvature_in_searchdir)

    ! -----------------------------------------------------------------------------------------------
    ! We now have a new search direction AND/OR a valid line step; must UPDATE the actual coordinates

    ! VA: Constrained relaxation + unit cell relaxation:
    if (relax_unit_cell > 0) then
       call cart2frac(lv, coords, frac_coords_temp)     ! Frac coord w.r.t old lv
       do i_mobile_coord = 1,n_mobile_lv_components     ! Get new lv
             lv(mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord)) = &
                  lv (mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord)) + &
                  line_step * search_direction(i_mobile_coord+n_mobile_atom_coords)
       end do
       ! Cart coords w.r.t. new lv -> frac coords remain const
       call frac2cart(lv, frac_coords_temp, coords)

       ! Assign coordinates if they were fixed (constrained) cartesian
       if (use_relaxation_constraints) then
          do i_atom = 1, n_atoms , 1
             if (coord_basis_atoms(i_atom) == 0) then  ! cartesian
                do jj = 1, 3
                   if (constrain_coords(jj, i_atom)) then
                      coords (jj, i_atom) = input_coords(jj, i_atom)    ! do whole atom since the atom is frac or cart fixed.
                   end if
                end do
             end if
          end do
       end if

       ! VA: Now update atoms that move -> frac/cart constrained atoms are fixed in frac/cart coords
       do i_mobile_coord = 1,n_mobile_atom_coords
             coords(mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord)) = &
                  coords(mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord)) + &
                  line_step * search_direction(i_mobile_coord)
       end do

    else  ! unit cell fix/cluster
       do i_mobile_coord = 1,n_mobile_lv_components     ! Get new lv
             lv(mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord)) = &
                  lv (mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord)) + &
                  line_step * search_direction(i_mobile_coord+n_mobile_atom_coords)
       end do
       do i_mobile_coord = 1,n_mobile_atom_coords
             coords(mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord)) = &
                  coords(mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord)) + &
                  line_step * search_direction(i_mobile_coord)
       end do
    end if

    ! relaxation exit options:
    valid_geo      = .true.
    ini_relaxation = .false.

    call save_relaxation_data_for_restart(coords, lv)

  end subroutine bfgs_descent
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/trusted_descent
  !  NAME
  !    trusted_descent
  !  SYNOPSIS

  subroutine trusted_descent &
                    (energy, coords, total_forces, lv, &
                     forces_lv, periodic_unit_cell_trans, &
                     KS_eigenvector, KS_eigenvector_complex, &
                     occ_numbers, valid_geo)

    !   (energy, coords, total_forces, lv, forces_lv, &
    !     KS_eigenvector, KS_eigenvector_complex, occ_numbers, valid_geo )

    !  REMARK: in principle it should be enough to pass the 'effective' coordinates (X) and forces (F)
    !          to the optimizer. However, in order to check that the atoms are not too close when a step
    !          is made, the whole structure information is needed.

    !  PURPOSE
    !    Driver routine for the trust-radius based relaxation - it contains the main decision tree.
    !
    !  USES
    use synchronize_mpi_basic, only: sync_real_number
    use trust_radius_method
    use bravais

    implicit none
    !  ARGUMENTS

    real*8 :: energy
    real*8, dimension(3,n_atoms)    :: coords
    real*8, dimension(3,n_atoms)    :: frac_coords_temp
    real*8, dimension(3,n_periodic) :: lv
    real*8, dimension(3,n_atoms)    :: total_forces
    real*8, dimension(3,n_periodic) :: forces_lv
    real*8, dimension(3,n_atoms)    :: periodic_unit_cell_trans
    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task)     :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points_task) :: occ_numbers
    logical :: valid_geo, constrain_frac

    ! INPUTS
    ! o energy -- total energy
    ! o X      -- Coordinates including lattice vector that the optimizer 'sees'.
    ! o F      -- Forces on coordinates including lattice vector that the optimizer 'sees'.
    ! o coords -- coordinates
    ! o lv     -- lattice vectors
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! OUTPUT
    ! o coords -- updated coordinates
    ! o lv     -- updated lattice vectors
    ! o valid_geo -- geometry checking flag
    ! o KS_eigenvector -- eigenvector, for resetting the BFGS
    ! o KS_eigenvector_complex -- same thing
    ! o occ_numbers -- occupation numbers (also for resetting)
    !
    ! AUTHOR
    !   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    ! HISTORY
    !   Release version, FHI-aims (2010).
    ! SOURCE

    ! Local variables
    ! Current iterate (last force evaluation)
    real*8, allocatable, dimension(:) :: X, F
    real*8 :: E
    real*8 :: V
    ! Last iterate (lowest energy before X, E, F)
    real*8, allocatable, dimension(:) :: X0, F0
    real*8 :: E0
    real*8 :: V0
    ! |X-X0|_2
    real*8 :: DX_norm
    ! true gain / expected gain; >1 is good, <0 is bad
    real*8 :: expected_gain, true_gain, harmonic_gain, quality_ratio
    ! next step
    real*8, allocatable, dimension(:) :: DX_step, DX_last, DX_tmp, HDX
    ! temporary coordinate difference
    real*8, allocatable, dimension(:,:) :: cleaned_hessian
    ! counters
    real*8 :: old_dist, this_dist, max_dist, max_radius, volume_change
    logical :: accept_step

    real*8 :: max_force_atom
    real*8, external :: get_max_force

    real*8 :: tmp_real
    real*8, external :: ddot

    character*150 :: info_str
    character(*), parameter :: func = 'trusted_descent'

    ! FK: If the new unit cell volume decreases percent wise by more than volume_threshold,
    !     we advance the geometry with a smaller trust radius.
    real*8, parameter :: volume_threshold = 0.3d0

    call aims_allocate(X,n_mobile_coords,"X")
    call aims_allocate(F,n_mobile_coords,"F")
    call aims_allocate(X0,n_mobile_coords,"X0")
    call aims_allocate(F0,n_mobile_coords,"F0")
    call aims_allocate(DX_step,n_mobile_coords,"DX_step")
    call aims_allocate(DX_last,n_mobile_coords,"DX_last")
    call aims_allocate(DX_tmp,n_mobile_coords,"DX_tmp")
    call aims_allocate(HDX,n_mobile_coords,"HDX")
    call aims_allocate(cleaned_hessian,nrow_hess,ncol_hess,"cleaned_hessian")

    call coords2X(coords, lv, X)
    call coords2X(total_forces, forces_lv, F)

    if (minimum_energy_ever_not_initialized) then
       minimum_energy_ever                 = energy
       minimum_energy_ever_not_initialized = .false.
    else
       if (energy < minimum_energy_ever) then
          minimum_energy_ever = energy
       end if
    end if

    E = energy
    ! MOL don't compute cell volume when symmetry-constrained relaxation
    if ((n_periodic /= 0) .and. (.not. use_symm_const_geo_relaxation)) then
       call get_cell_volume(lv, V0 )
    end if

    write(info_str, "(2X, 'Advancing geometry using trust radius method.')")
    call localorb_info(info_str)

    if (n_periodic /= 0) then
       max_radius = 0.5d0 * minval(length_of_lattice_vector)
    else
       max_radius = huge(max_radius)
    end if

    if (ini_relaxation) then
      ! allocate data for storing the KS eigenvectors during relaxation ...
      call bfgs_allocate_eigenvectors()

      stored_energy                   = energy
      stored_forces                   = total_forces
      stored_coords                   = coords
      stored_forces_lv                = forces_lv
      stored_lv                       = lv
      stored_periodic_unit_cell_trans = periodic_unit_cell_trans
      stored_DX                       = 0.0d0

      call bfgs_store_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

      if (trust_radius > max_radius) then
         trust_radius = max_radius
      end if

    else ! -> .not. ini_relaxation

       ! - Update Hessian

       ! Get X0, E0, F0, DX
       call coords2X(stored_coords, stored_lv,  X0)
       call coords2X(stored_forces, stored_forces_lv, F0)

       E0      = stored_energy
       ! DX_last = X - X0
       !> FlK: replace on-the-fly computation of DX by taking the
       ! saved DX from last step. This avoids confusion when the atomic
       ! positions have been scaled with the lattice.
       DX_last = stored_DX
       ! Mapping to the center does not work in reduced space
       if (.not. use_symm_const_geo_relaxation) then
         call map_DX_to_center_cell(DX_last)
       end if

       DX_norm = sqrt(sum(DX_last**2))

       if (DX_norm > trust_radius + 1.d-5) then
          write(info_str, "(2X,A,ES10.2,2A,ES10.2,A)") &
          & '** Last step length is', DX_norm*bohr, ' A ', &
          & 'instead of', trust_radius*bohr, ' A.'
          call localorb_info(info_str)
       end if

       ! Evaluate expected energy gain
       ! HDX = matmul(hessian, DX)
       if (is_worker) then
          call dgemv('N', nrow_hess, ncol_hess, 1.d0, hessian, nrow_hess, &
                  DX_last, 1, 0.d0, HDX(1+offset), 1)
       end if

       ! Textbook trust radius
       ! expected_gain = - dot_product(F0, DX_last) + 0.5d0 * dot_product(DX_last, HDX)
       expected_gain = ddot(n_mobile_coords, F0, 1, DX_last, 1)

       if (is_worker) then
          tmp_real = ddot(nrow_hess, DX_last(1+offset), 1, HDX(1+offset), 1)
       else
          tmp_real = 0.d0
       end if

       if (use_distributed_hessian) then
          call sync_real_number(tmp_real)
       end if

       expected_gain = - expected_gain + 0.5d0 * tmp_real

       if (expected_gain > 0.d0) then
          call localorb_info('  ** Detected expected_gain > 0: Force revert.')
          call localorb_info('  ** This should only happen for BFGS->TRM switch')
          expected_gain = - 1.d10
       end if

       ! Actually update Hessian
       call trm_BFGS_update(F-F0, DX_last, hessian)

       ! - Gains and Quality

       ! actual energy difference
       true_gain = energy - stored_energy
       write(info_str, "(2X,A,ES10.2,2A,ES10.2,A,F10.4)") &
       & '| True / expected gain:', true_gain*hartree, ' eV ', &
       & '/', expected_gain*hartree, ' eV =', true_gain / expected_gain
       call localorb_info(info_str)

       ! energy difference from forces assuming harmonic PES
       ! harmonic_gain = - dot_product(DX_last, 0.5d0*(F0+F))
       harmonic_gain = - ddot(n_mobile_coords, DX_last, 1, 0.5d0*(F0+F), 1)

       write(info_str, "(2X,A,ES10.2,2A,ES10.2,A,F10.4)") &
       & '| Harmonic / expected gain:', harmonic_gain*hartree, ' eV ', &
       & '/', expected_gain*hartree, ' eV =', harmonic_gain / expected_gain
       call localorb_info(info_str)


       ! Notify that the energy gain per atom becomes quite small.
       if (abs(true_gain*hartree / n_atoms) < 1e-5) then
         write(info_str, "(2X,A,ES10.2,2A)") &
         & '* True gain / number of atoms:', true_gain*hartree / n_atoms, ' eV '
         call localorb_info(info_str)

         ! ! Mix in the harmonic_gain to smooth out the decision
         ! if (true_gain > 0.d0) then
         !   quality_ratio = 0.1d0 * (9.d0*true_gain + harmonic_gain) / expected_gain
         ! else
         !   quality_ratio = 0.2d0 * (4.d0*true_gain + harmonic_gain) / expected_gain
         ! end if

         ! accept_step = (quality_ratio > 0.1d0)  ! Slightly more conservative
         ! write(info_str, "(2X,'| ',A)") &
         ! & 'Using harmonic gain <DX|-(F0+F1)/2> instead of DE to judge step.'
         ! call localorb_info(info_str)
       end if

       if (DX_norm > harmonic_length_scale) then
         ! Textbook TRM:
         quality_ratio = true_gain / expected_gain
         ! JW: Strictly speaking, the threshold could be even slightly
         ! positive.  See Nocedal, Wright: "Numerical Optimization".
         accept_step = (quality_ratio > 0.d0)
       else
         ! Wiggle adjustment:
         ! JW: At these length scales, the user (or the default) ensures us
         ! that the actual potential energy surface is at least qualitatively
         ! harmonic, i.e., we do not expect "hidden" features on the line
         ! between X and X0.  If the energy difference does not match the
         ! value expected from the integrated force (dot(DX, F_avg)), it is
         ! most likely the energy which got it wrong (grid errors tend to
         ! affect the energy stronger than the forces).  Do not look at the
         ! energy in this case.

         quality_ratio = harmonic_gain / expected_gain

         accept_step = (quality_ratio > 0.1d0)  ! Slightly more conservative
         write(info_str, "(2X,'| ',A)") &
         & 'Using harmonic gain <DX|-(F0+F1)/2> instead of DE to judge step.'
         call localorb_info(info_str)
       end if

       ! - Check forces-energy consistency

       if (accept_step .and. true_gain > energy_tolerance) then
          ! This only happens for
          !   * |DX| < harmonic_length_scale [~ 0.025 AA]
          !   * <DX|-(F0+F1)/2> < 0
          !   * DE > energy_tolerance [~ 0.2 meV]
          ! That is, the step is pretty short and the forces point into a
          ! direction where the energy goes up; forces and energy are not
          ! consistent anymore.
          ! I'm not completely sure if this is the right^TM criterion to stop.
          ! But at least, it is well adjustable in control.in by
          ! harmonic_length_scale and energy_tolerance. -- JW
          write(info_str, "(2X,'** ',A)") &
          & 'Inconsistency of forces<->energy above specified tolerance. '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'Specifically, we have: '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A,ES10.2,A)") &
          & 'True total energy step DE (positive means uphill) : ', true_gain*hartree, ' eV.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A,ES10.2,A)") &
          & "Current 'energy_tolerance' value for uphill steps : ", energy_tolerance*hartree, ' eV.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'Expected energy gain based on harmonic approximation: <DX | -(F0+F1)/2> /= DE:'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',14X,4(ES10.2,A))") &
          & DX_norm*bohr, ' A x', &
          & harmonic_gain*hartree/(DX_norm*bohr), ' eV/A =', &
          & harmonic_gain*hartree, ' eV /=', &
          & true_gain*hartree, ' eV.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'In principle, there should be no such inconsistencies.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'In fact, it is possible that the relaxation algorithm mis-estimated and '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'that a simple restart of the relaxation is all that is needed.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'We recommend to proceed as follows:'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ''
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '(1) Simply continue the calculation by copying the "geometry.in.next_step" file'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    produced in the last relaxation step to "geometry.in" and restarting FHI-aims.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ''
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '(2) Use the "get_relaxation_info.pl" script found in the "utilities" folder to'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '   monitor the progress of the relaxation. This will tell you if the relaxation is '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '   proceeding downwards in energy (the forces can fluctuate but the energy should generally go down).'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ''
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '(3) The tolerances of the standard relaxation algorithms can be increased to ignore'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    larger inconsistencies. To try this, increase "energy_tolerance" or'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    decrease "harmonic_length_scale".'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ''
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '(4) The initial guess for the second derivative (Hessian matrix) of the energy landscape'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    may not reflect your system well. The "init_hess" keyword can be used to specify '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    a more conservative initial guess. For example, "init_hess diag 0.97" is a rather'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    conservative choice.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ' '
          call localorb_info(info_str)
          if (use_mbd_old .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd) then
             write(info_str, "(2X,'** ',A)") &
             & '(5) You are using a variant of the many_body_dispersion method.'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & '    Unfortunately, the many_body_dispersion correction in its '
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & '    current default implementation can lead to small inconsistencies that can cause'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & "    FHI-aims' relaxation algorithm to stop with this warning. The warning is real,"
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & '    but the inconsistency may still be smaller than the numerical accuracy needed'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & '    for your production simulations. Please monitor the total energy change along'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & '    the relaxation trajectory carefully and continue the relaxation as outlined above.'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & ' '
             call localorb_info(info_str)
          end if
          write(info_str, "(2X,'** ',A)") &
          & '    Please also simply contact us for help and advice, especially if you see this'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    error message too often despite the steps outlined above. If there is a real'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & '    force-energy inconsistency, we would like to know about it.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & ' '
          call localorb_info(info_str)
          call aims_stop_coll('Numerical inconsistency of forces and energy above energy_tolerance.', func)
       end if

       ! Ensure that we did not take multiple uphill steps that took us uphill in small
       ! increments.
       if (accept_step .and. (energy-minimum_energy_ever) > aggregated_energy_tolerance) then
          write(info_str, "(2X,'** ',A)") &
          & 'The relaxation trajectory moved uphill from its lowest known energy so far '
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A,F20.5,A)") &
          & 'by more than: ', aggregated_energy_tolerance*hartree, ' eV.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & 'This likely means that the relaxation trajectory went uphill in a series of multiple small steps.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "In general, this should not happen - we therefore stop here to warn the user at this stage."
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "Please check (e.g., using the get_relaxation_info.pl utility) whether the relaxation indeed"
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "climbed uphill in multiple steps."
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "Possible remedies include: "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "(1) Just restart the relaxation using geometry.in.next_step as the new geometry.in. "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "(2) You can adjust the maximum allowed increase using the aggregated_energy_tolerance keyword. "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "(3) If you indeed observe multiple successive uphill energy steps that were accepted by "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "    the relaxation algorithm, please let us know. This should not happen and actually,"
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "    we are not aware under which circumstances it might happen. (Would be good to know.)"
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & "Sorry about this inconvenience. "
          call localorb_info(info_str)
          write(info_str, "(2X,'** ',A)") &
          & " "
          call localorb_info(info_str)
          call aims_stop_coll('Too many uphill relaxation steps, overall above aggregated_energy_tolerance.', func)
       end if

       ! - Update trust radius

       if (quality_ratio < 0.25d0) then
         ! Poor step: shorten trust radius
         trust_radius = 0.5d0 * min(trust_radius, DX_norm)
         write(info_str, "(2X,'| Reduce trust radius to ',ES10.2,' A.')") trust_radius*bohr
         call localorb_info(info_str)
        !> FlK: maybe have an intermediate criterion?
        ! else if ((quality_ratio < 0.5d0) .or. (quality_ratio > 1.75d0)) then
        !   ! Poor step: shorten trust radius a bit
        !   trust_radius = 0.75d0 * min(trust_radius, DX_norm)
        !   write(info_str, "(2X,'| Slightly reduce trust radius to ',ES10.2,' A.')") &
        !     trust_radius*bohr
        !   call localorb_info(info_str)
       else if (quality_ratio > 0.75d0 .and. DX_norm > 0.8d0*trust_radius) then
         ! Good step: show courage
         trust_radius = 2.d0 * trust_radius
         trust_radius = min(trust_radius, max_radius)
         write(info_str, "(2X,'| Enlarge trust radius to ',ES10.2,' A.')") trust_radius*bohr
         call localorb_info(info_str)
       end if

       ! - Revert or Accept (i.e. store)

       if (accept_step) then
          stored_energy                   = energy
          stored_forces                   = total_forces
          stored_coords                   = coords
          stored_forces_lv                = forces_lv
          stored_lv                       = lv
          stored_DX                       = DX_step
          stored_periodic_unit_cell_trans = periodic_unit_cell_trans

          call bfgs_store_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

       else
          write(info_str, "(2X,'| Counterproductive step -> revert!')")
          call localorb_info(info_str)
          DX_tmp = X0 - X
          DX_last = 0.
          ! TARP: Mapping to the center does not work for SCR relaxation
          if (.not. use_symm_const_geo_relaxation) then
            call map_DX_to_center_cell(DX_tmp)
          end if
          X = X + DX_tmp   ! == X0 (as near to X as possible)
          F = F0
          E = E0
          ! JW: It is important not to do "coords = stored_coords", because we
          ! must not change the constrained coordinates from within this
          ! procedure (these coordinates may have been changed by a
          ! map_to_center_cell() outside).
          if (.not. use_symm_const_geo_relaxation) then
            call ensure_constraints(coords, lv, X)
          endif

          call X2coords(coords, lv,  X)

          call bfgs_load_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)
       end if

       ! For the RI version of LVL_fast as of August, 2012, it is possible to find entirely
       ! nonsensical forces due to an unpatched ill-conditioning problem that affects the
       ! forces, not the energies. In this case, we need a graceful exit.
       if ((.not. accept_step) .and. (RI_type == RI_LVL)) then
          max_force_atom = get_max_force(total_forces)
          if (max_force_atom > 1.d0) then
             write(info_str, "(2X,'** ',A)") &
             & 'Warning. This is a relaxation involving a hybrid functional,'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'and a large force-energy inconsistency appears to have been detected.'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'In the present version of FHI-aims, this functionality is still exerimental,'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'and only relatively small basis sets (tier 1) avoid such inconsistencies.'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'We stop the run here out of precaution. If you think that this message'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'does not apply to your case, let us know - we are actively working on'
             call localorb_info(info_str)
             write(info_str, "(2X,'** ',A)") &
             & 'repairing this issue.'
             call localorb_info(info_str)

             call aims_stop_coll('LVL force inconsistency?','trusted_descent')
          end if
       end if
    end if ! -> .not. ini_relaxation

    !MOL: in symmetry-constrained relaxation atom parameters and lattice parameters
    !     are independent of each other --> delete interaction terms in upper right
    !     and lower left parts of Hessian (lv-atom interaction)
    ! TARP: Removing for now as parametric relaxation performs better with the coupling
    ! if (use_symm_const_geo_relaxation) then
    !    hessian(1:SCR_n_params_coords,SCR_n_params_coords+1:SCR_n_params) = 0.
    !    hessian(SCR_n_params_coords+1:SCR_n_params,1:SCR_n_params_coords) = 0.
    ! endif

    ! - Step

    ! Get step
    ! MOL: Don't clean in case of symm-constrained relaxation, only works with
    ! full symmetrized forces but not the symmetry-reduced forces/hessian here

    cleaned_hessian = hessian
    if (.not. use_symm_const_geo_relaxation) then
        call clean_hessian(cleaned_hessian, X)
    endif

    call trm_trusted_step(F, cleaned_hessian, trust_radius, DX_step)

    if (.not. use_symm_const_geo_relaxation) then
        call clean_mobile_force_components(DX_step, X)
    endif

    ! convert X+DX to coords
    if (.not. use_symm_const_geo_relaxation) then
      call min_atomic_dist(lv, coords, old_dist)
      ! call ensure_constraints(coords, lv, X+DX_step)
      !> FlK: integrated within X2new_coords (?)
    endif
    call X2new_coords(coords, lv, X, DX_step)
    ! MOL: WARNING: if you ever decide to scale the new coordinates to the new lattice
    ! vectors (which is already done in BFGS_textbook), please make sure to exclude the
    ! case of the symmetry-constrained by checking the flag `use_symm_const_geo_relaxation`
    ! relaxation because it works in a parameter space where the atom parameters are
    ! independent of the lattice parameters. Thank you.
    ! FlK: This is being taken care of be X2new_coords

    ! check for close encounters
    if (.not. use_symm_const_geo_relaxation) then
      call min_atomic_dist(lv, coords, this_dist)
      max_dist = min(old_dist-0.3d0/bohr, 0.69d0/bohr)

      if (this_dist < max_dist) then
         ! The smallest bond length is both smaller than H-H (0.69 A)
         ! and significantly smaller than in the last iteration.
         ! -> Safe guard:
         write(info_str, "(2X,'| Use smaller step to avoid collisions.')")
         call localorb_info(info_str)
         call trm_trusted_step(F, cleaned_hessian, 0.2d0/bohr, DX_step)
         call clean_mobile_force_components(DX_step, X)

         !> FlK: Before we can move with the new step, we have to restore the old geometry:
          if (.not. use_symm_const_geo_relaxation) then
            call ensure_constraints(coords, lv, X)
          endif
          call X2coords(coords, lv,  X)
          call X2new_coords(coords, lv, X, DX_step)
      end if
    endif

    ! FK: Check for too fast decreasing unit cell volume
    ! MOL: but not if symm-const. relaxation
    if ((n_periodic /= 0).and.(.not. use_symm_const_geo_relaxation)) then
       call get_cell_volume(lv ,V)
       volume_change = (V0-V)/V0

       if (volume_change > volume_threshold) then
          ! FK: The unit cell volume decreased percent wise more than volume_threshold.
          !     We don't want this to avoid the collapse of the unit cell.
          write(info_str, "(2X,A,F4.1,A)") '** The unit cell volume has decreased by ', 100*volume_change, ' percent.'
          call localorb_info(info_str)
          write(info_str, "(2X,A,F4.1,A)") '** This is more than the threshold of ', 100*volume_threshold, ' percent.'
          call localorb_info(info_str)
          write(info_str, "(2X,'** Use smaller step to avoid collapse of unit cell.')")
          call localorb_info(info_str)
          ! FK: Advance geometry with smaller radius
          call trm_trusted_step(F, cleaned_hessian, 0.4d0/bohr, DX_step)
          call clean_mobile_force_components(DX_step, X)

          !> FlK: Before we can move with the new step, we have to restore the old geometry:
          if (.not. use_symm_const_geo_relaxation) then
            call ensure_constraints(coords, lv, X)
          endif
          call X2coords(coords, lv,  X)
          call X2new_coords(coords, lv, X, DX_step)
       end if
    end if

    ! relaxation exit options:
    valid_geo      = .true.
    ini_relaxation = .false.

    stored_DX = DX_step

    call aims_deallocate(X,"X")
    call aims_deallocate(X0,"X0")
    call aims_deallocate(F,"F")
    call aims_deallocate(F0,"F0")
    call aims_deallocate(DX_step,"DX_step")
    call aims_deallocate(DX_last,"DX_last")
    call aims_deallocate(DX_tmp,"DX_tmp")
    call aims_deallocate(HDX,"HDX")
    call aims_deallocate(cleaned_hessian,"cleaned_hessian")

    call save_relaxation_data_for_restart(coords, lv)

  end subroutine trusted_descent
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/coords2X
  !  NAME
  !    coords2X
  !  SYNOPSIS

  subroutine coords2X(coords, lattice_vector, X)

    !  PURPOSE
    !     Convert the two sets of mobile_atom_coords(1:3, 1:n_mobile_atoms)
    !       and mobile_lv_components (1:3,1:n_mobile_lv) to  common 1d array X(1:n_mobile_atoms).
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN)      :: coords(3, n_atoms)
    real*8, intent(IN)      :: lattice_vector(3, n_periodic)
    real*8, intent(OUT)     :: X(n_mobile_coords)

    !  INPUTS
    !    o coords -- Complete set of Cartesian coordinates
    !  OUTPUTS
    !    o X -- Reduced set of (unconstrained) Cartesian coordinates
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_mobile_coord, i_atom, i_lv, i

    X = 0   !MOL: initialize in case n_mobile_atom_coords and n_mobile_lv_components are smaller

    do i_mobile_coord = 1, n_mobile_atom_coords,1
       i_atom = mobile_atom(i_mobile_coord)
       i = mobile_atom_coord(i_mobile_coord)
       X(i_mobile_coord) = coords(i, i_atom)
    end do
   ! VA: similar for lattice_vectors
    do i_mobile_coord = 1, n_mobile_lv_components,1
       i_lv = mobile_lv (i_mobile_coord)
       i = mobile_lv_component(i_mobile_coord)
       X(i_mobile_coord+n_mobile_atom_coords) = lattice_vector(i, i_lv)
    end do

  end subroutine coords2X
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/X2coords
  !  NAME
  !    X2coords
  !  SYNOPSIS

  subroutine X2coords(coords, lattice_vector, X)

    !  PURPOSE
    !     Convert the set of  X(1:n_mobile_atoms) to coords(1:3, 1:n_atoms).
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: coords(3, n_atoms)
    real*8, intent(INOUT) :: lattice_vector(3, n_periodic)
    real*8, intent(IN)    :: X(n_mobile_coords)

    !  INPUTS
    !    o coords -- Some reference geometry (needed for constraints)
    !    o X -- Reduced set of (unconstrained) Cartesian coordinates
    !  OUTPUTS
    !    o coords -- Complete set of Cartesian coordinates
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_mobile_coord, i_atom, i_lv, i_coord
    real*8  :: frac_coords_temp (3, n_atoms)

    do i_mobile_coord = 1, n_mobile_atom_coords
       i_atom = mobile_atom(i_mobile_coord)
       i_coord = mobile_atom_coord(i_mobile_coord)
       coords(i_coord, i_atom) = X(i_mobile_coord)
    end do

   ! VA: lattice_vectors
    do i_mobile_coord = 1, n_mobile_lv_components
       i_lv = mobile_lv(i_mobile_coord)
       i_coord = mobile_lv_component(i_mobile_coord)
       lattice_vector(i_coord, i_lv) = X(i_mobile_coord+n_mobile_atom_coords)
    end do

    end subroutine X2coords

    subroutine X2new_coords(coords, lattice_vector, X, DX) !, DX_last)

      !  PURPOSE
      !    Convert the set of X(1:n_mobile_atoms) to new_coords(1:3, 1:n_atoms).
      !    When atomic positions are given in fractional coordinates,
      !    changing the lattice will leave the fractional coordinates
      !    invariant, i.e., the Cartesian coordinates are scaled
      !    accordingly.

      implicit none

      real*8, intent(INOUT)         :: coords(3, n_atoms)
      real*8, intent(INOUT)         :: lattice_vector(3, n_periodic)
      real*8, intent(IN)            :: X(n_mobile_coords)
      real*8, intent(IN)            :: DX(n_mobile_coords)

      !  AUTHOR
      !    Florian Knoop (FlK)
      !  HISTORY
      !    06.08.2018 Beijing: Initial discussion with Volker
      !    08.08.: Step DX need to know about this as well
      !    10.08.: Remember to scale DX as well. Make more elegant by
      !            setting positions with it?
      !    31.08.: Push to master
      !    02.11.: Necessary fixes

      integer :: i_mobile_coord, i_atom, i_lv, i_coord
      real*8  :: frac_coords_temp (3, n_atoms)
      real*8  :: frac_coords_init (3, n_atoms)
      real*8  :: scaled_coords(3, n_atoms)
      real*8  :: scaled_input_coords(3, n_atoms)

      ! Initial frac coords
      if (.not. use_symm_const_geo_relaxation) then
        call cart2frac(lattice_vector, coords, frac_coords_init)
      end if

      ! Move atomic positions
      do i_mobile_coord = 1, n_mobile_atom_coords
        i_atom = mobile_atom(i_mobile_coord)
        i_coord = mobile_atom_coord(i_mobile_coord)
        coords(i_coord, i_atom) = X(i_mobile_coord) + DX(i_mobile_coord)
      end do

      ! Save the fractional coordinates of the new atomic positions
      if (((relax_unit_cell == 1) .and. (n_periodic>0)) .and. &
        (.not. use_symm_const_geo_relaxation) .and. (relax_geo.eq.RELAX_LTRM)) &
        then

        ! Frac coord w.r.t old lv
        call cart2frac(lattice_vector, coords, frac_coords_temp)
      end if

     ! Move lattice_vectors
      do i_mobile_coord = 1, n_mobile_lv_components
        i_lv = mobile_lv(i_mobile_coord)
        i_coord = mobile_lv_component(i_mobile_coord)
        lattice_vector(i_coord, i_lv) = &
            X(i_mobile_coord+n_mobile_atom_coords) &
            + DX(i_mobile_coord+n_mobile_atom_coords)
      end do

      ! now set the new positions by assigning the fractional coords _after_
      ! the atoms were moved (if it is actually what we want)
      if (((relax_unit_cell == 1) .and. (n_periodic>0)) .and. &
        (.not. use_symm_const_geo_relaxation)) then

        if (relax_geo .eq. RELAX_LTRM) then

          call frac2cart(lattice_vector, frac_coords_temp, scaled_coords)
          coords = scaled_coords

        end if

        ! This is from ensure_constraints:
        !(2) Overwrite those coordinates which are fixed in a kartesian basis.
        ! if this is wanted:
        if (use_relaxation_constraints) then
          do i_atom = 1, n_atoms , 1
            ! FlK: check if positions are cartesian and if there are relaxation
            ! constraints
            if (coord_basis_atoms(i_atom) == 0) then
              if (constrain_relaxation(i_atom)) then
                coords (:,i_atom) = input_coords(:,i_atom)
              else
                do i_coord = 1, 3
                  if (constrain_coords(i_coord, i_atom)) then
                    coords (i_coord, i_atom) = input_coords(i_coord, i_atom)
                  end if
                end do
              end if
            !> 28.11.2018: fix for fractionally constrained atoms:
            else if (coord_basis_atoms(i_atom) == 1) then
              ! New cartesian position of fractionally constrained atom
              call frac2cart(lattice_vector, frac_coords_init, scaled_input_coords)
              if (constrain_relaxation(i_atom)) then
                coords (:,i_atom) = scaled_input_coords(:,i_atom)
              else
                do i_coord = 1, 3
                  if (constrain_coords(i_coord, i_atom)) then
                    coords (i_coord, i_atom) = scaled_input_coords(i_coord, i_atom)
                  end if
                end do
              end if
            else
              write(info_str,'(1X,A)') &
              "* Don't know in which basis the atoms are constrained !"
              call localorb_info(info_str,use_unit,'(A)')
              call aims_stop ()
            end if
          end do
        end if
      end if
    end subroutine X2new_coords

   !******
  !------------------------------------------------------------------------------
  !****s* relaxation/map_DX_to_center_cell
  !  NAME
  !    map_DX_to_center_cell
  !  SYNOPSIS

  subroutine map_DX_to_center_cell(DX)

    !  PURPOSE
    !     Map a geometry difference vector DX(1:n_mobile_atoms) to the center cell.
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(INOUT) :: DX(n_mobile_coords)

    !  INPUTS
    !    o X -- Reduced set of (unconstrained) Cartesian coordinates
    !  OUTPUTS
    !    o coords -- Complete set of Cartesian coordinates
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: dcoord(3, n_atoms)
    real*8 :: dlv(3, n_periodic)
    real*8 :: new_lattice_vector(3, n_periodic)
    real*8 :: new_map_to_center_cell_matrix(3, n_periodic)
    real*8 :: length(3)
    integer,dimension(3):: ipivot
    integer :: i_atom, i_latt

    dcoord = 0.d0
    dlv    = 0.d0

    ! VA since we use dummy variable lattice_vector
    if (n_periodic /= 0) then
       call X2coords(dcoord, dlv, DX)

       ! in case that lattice vectors changed
       new_lattice_vector            = lattice_vector + dlv
       new_map_to_center_cell_matrix = lattice_vector

       ! Get inverse matrix
       call dgetrf(3, 3, new_map_to_center_cell_matrix, 3, ipivot, info)
       if (info /= 0) then
          write(info_str,'(X,A)') &
             "*** New lattice_vector is singular."
          call localorb_info(info_str,use_unit,'(A)')
          call aims_stop()
       end if

       call dgetri(3, new_map_to_center_cell_matrix, 3, ipivot, work, 3, info)
       if (info /= 0) then
          write(info_str,'(X,A)') &
             "*** New lattice_vector is singular."
          call localorb_info(info_str,use_unit,'(A)')
          call aims_stop()
       end if

       do i_atom = 1, n_atoms
          length           = matmul(new_map_to_center_cell_matrix, dcoord(:,i_atom))
          length           = length - nint(length)
          dcoord(:,i_atom) = matmul(lattice_vector, length)
       end do

       call coords2X(dcoord, dlv, DX)
    end if

  end subroutine map_DX_to_center_cell
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/ensure_constraints
  !  NAME
  !    ensure_constraints
  !  SYNOPSIS

  subroutine ensure_constraints(coords, lv, X)

    !  PURPOSE
    !    Ensures constrained atoms for unit cell translations.
    !    If the unit cell varies, then fractional and cartesian
    !    constraints give different results.
    !
    !  USES
    implicit none
    !  ARGUMENTS
    real*8, intent(INOUT)   :: coords(3, n_atoms)
    real*8, intent(INOUT)   :: lv(3, n_periodic)
    real*8, intent(IN)      :: X(n_mobile_coords)
    real*8                  :: frac_coords_temp (3, n_atoms)
    integer                 :: i_atom

    !  INPUTS
    !   o coords
    !   o lv
    !
    !  OUTPUTS
    !   o coords
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    if ((relax_unit_cell == 0) .or. (n_periodic==0)) return

    !(1) Constrain all coordinates fractionally
    call cart2frac(lv, coords, frac_coords_temp)     ! Frac coord w.r.t old lv
    call X2coords (coords, lv, X)                    ! Get new lv
    call frac2cart(lv, frac_coords_temp, coords)     ! Cart coords w.r.t. new lv -> frac coords remain const.

    !(2) Overwrite those coordinates which are fixed in a kartesian basis.
    do i_atom = 1, n_atoms , 1

       if (coord_basis_atoms(i_atom) == 0) then        ! cartesian
           coords (:,i_atom) = input_coords(:,i_atom)  ! do whole atom since the atom is frac or cart fixed.
       else if (coord_basis_atoms(i_atom) == 1) then
          cycle
       else
           write(info_str,'(1X,A)') &
           "* Don't know in which basis the atoms are constrained !"
           call localorb_info(info_str,use_unit,'(A)')
           call aims_stop ()
       end if
    end do

  end subroutine ensure_constraints
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/get_force_times_searchdir_hat
  !  NAME
  !    get_force_times_searchdir_hat
  !  SYNOPSIS

  subroutine get_force_times_searchdir_hat &
             (forces,                      &
              forces_lv,                   &
              search_dir,                  &
              force_times_searchdir_hat)

    !  PURPOSE
    !    calculate single vector product; required for BFGS and line searches
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms)       :: forces
    real*8, dimension(3,n_periodic)    :: forces_lv
    real*8, dimension(n_mobile_coords) :: search_dir
    real*8 :: force_times_searchdir_hat

    !  INPUTS
    !   o forces
    !   o search_dir
    !
    !  OUTPUTS
    !   o force_times_searchdir_hat
    !
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: ddot
    integer :: i_coord
    force_times_searchdir_hat = 0d0
    do i_coord = 1, n_mobile_atom_coords
       force_times_searchdir_hat = force_times_searchdir_hat+ &
          forces(mobile_atom_coord(i_coord),mobile_atom(i_coord)) &
          *search_dir(i_coord)
    end do

    ! VA similar for lattice vectors
    do i_coord = 1, n_mobile_lv_components
          force_times_searchdir_hat = force_times_searchdir_hat+ &
               forces_lv(mobile_lv_component(i_coord),mobile_lv(i_coord))&
               *search_dir(n_mobile_atom_coords + i_coord)
    end do

    ! Normalize
    force_times_searchdir_hat = force_times_searchdir_hat / &
         sqrt(ddot(n_mobile_coords,search_dir,1,search_dir,1))

  end subroutine get_force_times_searchdir_hat
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/update_hessian_and_searchdir_steepest_descent
  !  NAME
  !    update_hessian_and_searchdir_steepest_descent
  !  SYNOPSIS

  subroutine update_hessian_and_searchdir_steepest_descent &
             (hessian_int,                                 &
              searchdir,                                   &
              forces,                                      &
              forces_lv)

    !  PURPOSE
    !  single steepest descent step; resets Hessian matrix to unity
    !  simple steepest descent type geometry relaxation for initial steps and useful when the
    !  BFGS probably fails, also improves the Hessian matrix of the system
    !
    !  The simple diagonal matrix is even used to reset the Hessian when the
    !  initial Hessian is specified to be the Lindh matrix.  When the code
    !  arrives here, there has already been some trouble, most likely due to
    !  wiggles in the energy landscape on the order of sub-meV and pm.  The
    !  Lindh matrix is not adequate for these wiggles.
    !
    !  USES
    use bravais

    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms)                          :: forces
    real*8, dimension(3,n_periodic)                       :: forces_lv
    real*8, dimension(n_mobile_coords, n_mobile_coords)   :: hessian_int
    real*8, dimension(9, 9)                               :: hessian_temp
    real*8, dimension(n_mobile_coords)                    :: searchdir
    real*8  :: recip_lattice_vector(3, 3), recip_lattice_vector_2(3, 3)
    integer :: i_mobile, i_mobile_lv, j_mobile_lv, ii, jj, n_block

    !  INPUTS
    !   o forces
    !   o hessian_int -- internal hessian
    !   o searchdir
    !
    !  OUTPUT
    !   o hessian_int
    !   o searchdir
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: i_hessian
    line_step = 1.0
    hessian_int = 0.d0
    do i_hessian = 1, n_mobile_coords, 1
       hessian_int(i_hessian,i_hessian) = init_hess_diag
    end do

    ! FlK: re-initialize hessian with reciprocal lattice if requested
    ! TARP: But not if using symmetry constraints
    if((.not. use_symm_const_geo_relaxation) .and. (init_hess_type == HESS_LATTICE .or. relax_geo == RELAX_LTRM)) then
      ! Set diagonal init_hess matrix:
      do i_mobile = 1, n_mobile_atom_coords
         hessian_int(i_mobile, i_mobile) = init_hess_diag
      end do

      call get_reciprocal_vectors(n_periodic, lattice_vector, recip_lattice_vector)

      recip_lattice_vector_2 = &
             matmul(recip_lattice_vector, transpose(recip_lattice_vector))

      do jj = 1, 3
         do ii = 1, 3
            do n_block = 0, 2
               hessian_temp(ii+n_block*3, jj+n_block*3) &
                   = recip_lattice_vector_2(ii, jj) &
                   * init_hess_lv_diag &
                   / frobnorm_3x3_real(recip_lattice_vector_2)
            end do
         end do
      end do

      ! Assign the components of the hessian that are allowed to move
      do j_mobile_lv = 1, n_mobile_lv_components
         do i_mobile_lv = 1, n_mobile_lv_components
            ! assign x,y,z from the mobile coordinates to match the
            ! correct line from the full initial hessian
            ii = (mobile_lv(i_mobile_lv) - 1) * 3 + mobile_lv_component(i_mobile_lv)
            jj = (mobile_lv(j_mobile_lv) - 1) * 3 + mobile_lv_component(j_mobile_lv)
            hessian(i_mobile_lv + n_mobile_atom_coords, &
                j_mobile_lv + n_mobile_atom_coords) &
                = hessian_temp(ii, jj)
         end do
      end do
    end if

    initial_hessian = .true.
    ! FK: not hessian_int ?
    call update_searchdir(hessian, searchdir, forces, forces_lv)

  end subroutine update_hessian_and_searchdir_steepest_descent
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/update_hessian_and_searchdir_BFGS
  !  NAME
  !    update_hessian_and_searchdir_BFGS
  !  SYNOPSIS

  subroutine update_hessian_and_searchdir_BFGS &
             (hessian_int,                     &
              searchdir,                       &
              forces,                          &
              forces_lv)

    !  PURPOSE
    !    Broyden-Fletcher-Goldfarb-Shanno method, straight from Wikipedia
    !    (http://en.wikipedia.org/wiki/BFGS_method, the description there appears to be correct as of Jan 28 '08)
    !    Line Search for the minimum along a given search direction may be done,
    !    but algorithm may also be run with a simple quadratic extrapolation (no line search)
    !
    !   Note: The update runs only over those n_mobile_coords that are allowed to move
    !   during a relaxation.
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms)                          :: forces
    real*8, dimension(3,n_periodic)                       :: forces_lv
    real*8, dimension(n_mobile_coords, n_mobile_coords)   :: hessian_int
    real*8, dimension(n_mobile_coords) :: searchdir

    ! INPUTS
    ! o forces
    ! o hessian_int -- internal hessian
    ! o searchdir
    !
    !  OUTPUT
    ! o hessian_int -- intenal hessian
    ! o searchdir
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8, dimension(n_mobile_coords) :: force_difference
    real*8 :: forcediff_times_searchdir, sHs, ddot
    real*8, dimension(n_mobile_coords) :: Hs
    integer :: i_mobile_coord

    ! compute force difference between present and previous geometry
    do i_mobile_coord = 1, n_mobile_atom_coords
       force_difference(i_mobile_coord) = &
             stored_forces(mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord)) &
           - forces       (mobile_atom_coord(i_mobile_coord),mobile_atom(i_mobile_coord))
    end do
    do i_mobile_coord = 1, n_mobile_lv_components
       force_difference(i_mobile_coord+n_mobile_atom_coords) = &
              stored_forces_lv (mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord)) &
            - forces_lv        (mobile_lv_component(i_mobile_coord),mobile_lv(i_mobile_coord))
    end do
    ! include the line step factor required by wikipedia
    force_difference = force_difference/line_step

    ! scalar product: force_difference * search_direction
    forcediff_times_searchdir = ddot(n_mobile_coords,force_difference,1,searchdir,1)

    ! Hessian times search_direction: Hs
    call dsymv( 'u',n_mobile_coords,1d0,hessian_int,n_mobile_coords,searchdir,1,0.d0,Hs,1 )

    ! scalar product: search_direction * hessian * search_direction
    sHs = ddot(n_mobile_coords,searchdir,1,Hs,1)

    ! Now update the hessian ...
    call dsyr( 'u', n_mobile_coords, 1.d0/forcediff_times_searchdir, force_difference, 1, hessian_int, n_mobile_coords )
    call dsyr( 'u', n_mobile_coords, -1.d0/sHs, Hs, 1, hessian_int, n_mobile_coords )

    if (initial_hessian) then
       t_first_pred = .true.
    else
       t_first_pred = .false.
    end if
    initial_hessian = .false.

    call update_searchdir(hessian_int, searchdir, forces, forces_lv)

  end subroutine update_hessian_and_searchdir_BFGS
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/update_searchdir
  !  NAME
  !    update_searchdir
  !  SYNOPSIS

  subroutine update_searchdir(hessian_int, searchdir, forces, forces_lv)

    !  PURPOSE
    !    Calculate search direction H^-1 F from Hessian (H) and forces (F).
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(3,n_atoms), intent(IN)    :: forces
    real*8, dimension(3,n_periodic), intent(IN) :: forces_lv
    real*8, dimension(n_mobile_coords, n_mobile_coords), intent(IN) :: hessian_int
    real*8, dimension(n_mobile_coords), intent(OUT) :: searchdir

    !  INPUTS
    !   o forces
    !   o hessian_int -- internal hessian
    !  OUTPUT
    !   o searchdir
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer :: info
    character*150 :: info_str

    ! Use current Hessian matrix and current forces to predict search
    ! direction s, by inverting hessian * search_direction = - total_forces

    ! preload current forces into search_direction (after dsysv, this array
    ! will contain the actual search direction) note that search_direction is
    ! a 1D array for technical reasons but contains the search direction in
    ! the same order as a 2D array search_direction(i_coord,i_atom)
    call coords2X(forces, forces_lv,  searchdir)

    ! actually solve system of linear equations (lapack) to obtain the next
    ! search direction
    hessian_work = hessian_int
    call dsysv ( 'u', n_mobile_coords, 1, hessian_work, n_mobile_coords, &
    &           ipiv, searchdir, n_mobile_coords, work, lwork, info )

    ! check whether dgesv gave a reasonable result:
    if (info.lt.0) then
       write(info_str,'(1X,A,I3,A)') &
            "* Search direction: Argument number ",-info," to dsysv had an illegal value."
       call localorb_info(info_str,use_unit,'(A)' )
       stop
    else if (info.gt.0) then
       write(info_str,'(1X,A,I3,A)') &
            "* Search direction: Hessian factorisation ",info," is exactly zero."
       call localorb_info(info_str,use_unit,'(A)' )
       write(info_str,'(1X,A)') &
            "Your input Hessian may be singular - cannot proceed."
       call localorb_info(info_str,use_unit,'(A)' )
       stop
    end if

  end subroutine update_searchdir
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/restrict_max_displacement
  !  NAME
  !    restrict_max_displacement
  !  SYNOPSIS

  subroutine restrict_max_displacement(search_dir, line_step, restricted, negative_curvature)

    !  PURPOSE
    !    keep maximal displacement below a certain value, set flags if the line step had to be restricted
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8 :: line_step
    real*8, dimension(n_mobile_coords) :: search_dir
    logical :: restricted, negative_curvature

    ! INPUTS
    ! o line_step
    ! o search_dir
    !
    !  OUTPUTS
    !   o restricted
    !   o negative_curvature
    !   o line_step
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: max_displacement
    integer :: i_mobile_coord

    max_displacement = abs(line_step*search_direction(1))
    do i_mobile_coord = 1, n_mobile_coords
          max_displacement = max(max_displacement,abs(line_step*search_dir(i_mobile_coord)))
    end do

    write(info_str,'(2X,A,E14.6,A)')"Predicted maximal displacement in this step: ",max_displacement*bohr," A"
    call localorb_info(info_str)

    if (max_displacement.gt.max_atomic_move) then
       line_step = line_step*max_atomic_move/max_displacement
       restricted = .true.
       write(info_str,'(2X,A)') &
            "This displacement is too big!!!"
       call localorb_info( info_str )
       write(info_str,'(2X,A,E14.6,A,E14.6,A)') &
            "Reducing line step to ", line_step, ", to avoid exceeding ", &
            max_atomic_move*bohr, " A."
       call localorb_info( info_str )
    else if (negative_curvature) then
       ! want to get into region with positive curvature, explicitly move with maximally
       ! allowed line step into search direction!!!
       restricted = .true.
       line_step = line_step*max_atomic_move/max_displacement
       write(info_str,'(2X,A)') &
            "Using maximally permissible displacement to get out of region of negative energy curvature"
       call localorb_info(info_str)
       write(info_str,'(2X,A,E14.6,A,E14.6,A)') "Increasing line step to " , &
            line_step, " in order to displace maximally ", max_atomic_move*bohr," A."
       call localorb_info(info_str)
    else
       restricted = .false.
    end if

  end subroutine restrict_max_displacement
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/bfgs_allocate_eigenvectors
  !  NAME
  !    bfgs_allocate_eigenvectors
  !  SYNOPSIS

  subroutine bfgs_allocate_eigenvectors

    !  PURPOSE
    !    Allocate space for (old) eigenvectors
    !
    !  USES
    implicit none
    !  ARGUMENTS

    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE
    character(*), parameter :: func = 'bfgs_allocate_eigenvectors'

    if (n_k_points_task > 0 .and. .not. use_scalapack .and. &
    &  .not. store_eigenvectors_to_disk_in_relaxation .and. &
    &   energy_tolerance < 1e10) then
       if (real_eigenvectors) then
          if (.not. allocated(stored_KS_eigenvector)) then
             call aims_allocate(stored_KS_eigenvector,n_basis,n_states,n_spin,n_k_points_task,"+stored_KS_eigenvector")
          end if
       else
          if (.not. allocated( stored_KS_eigenvector_complex)) then
             call aims_allocate(stored_KS_eigenvector_complex,n_basis,n_states,n_spin,n_k_points_task,"+stored_KS_eigenvector_complex")
          end if
       end if
       if (.not. allocated(stored_occ_numbers)) then
          call aims_allocate(stored_occ_numbers,n_states,n_spin,n_k_points_task,"stored_occ_numbers")
       end if
    end if

  end subroutine bfgs_allocate_eigenvectors
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/bfgs_store_eigenvectors
  !  NAME
  !    bfgs_store_eigenvectors
  !  SYNOPSIS

  subroutine bfgs_store_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

    !  PURPOSE
    !    keep eigenvectors in case BFGS step fails: shortens the scf cycle
    !
    !  USES
    implicit none
    !  ARGUMENTS

    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points_task) :: occ_numbers

    !  INPUTS
    !    o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !    o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_k, i_basis, i_states, i_spin
    character*40 :: form1

    if (store_eigenvectors_to_disk_in_relaxation) then

       open(file = file_name, unit = 7, status = 'unknown')

       if (real_eigenvectors) then
          write(form1,'(A,I3,A)') '(',n_spin,'E30.18E4)'
       else
          write(form1,'(A,I3,A)') '(',2*n_spin,'E30.18E4)'
          write(form1,'(A,I3,A)') '(',2*n_spin,'E30.18E4)'
       end if

       if (real_eigenvectors) then
          ! KS_eigenvector(n_basis,n_states,n_spin,n_k_points_task)
          do i_k = 1, n_k_points_task
             do i_basis = 1, n_basis
                do i_states = 1, n_states
                   write(unit=7,fmt=form1) (KS_eigenvector(i_basis,i_states,i_spin,i_k),i_spin = 1, n_spin)
                end do
             end do
          end do
       else
          do i_k = 1, n_k_points_task
             do i_basis = 1, n_basis
                do i_states = 1, n_states
                   write(unit=7,fmt=form1) (KS_eigenvector_complex(i_basis,i_states,i_spin,i_k) ,i_spin = 1, n_spin)
                end do
             end do
          end do
       end if

       do i_k = 1, n_k_points_task
          do i_states = 1, n_states
             write(unit=7,fmt="(2F26.18)") (occ_numbers(i_states,i_spin,i_k) ,i_spin = 1, n_spin)
          end do
       end do

       close(7)

    else

       if (energy_tolerance < 1.d10) then
          if (use_scalapack) then

             call store_eigenvectors_scalapack

          else

             if (n_k_points_task > 0) then
                if (real_eigenvectors) then
                   stored_KS_eigenvector = KS_eigenvector
                else
                   stored_KS_eigenvector_complex = KS_eigenvector_complex
                end if
                stored_occ_numbers = occ_numbers
             end if
          end if
       end if
    end if

    eigenvectors_stored = .true.

  end subroutine bfgs_store_eigenvectors
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/bfgs_load_eigenvectors
  !  NAME
  !    bfgs_load_eigenvectors
  !  SYNOPSIS

  subroutine bfgs_load_eigenvectors(KS_eigenvector, KS_eigenvector_complex, occ_numbers)

    !  PURPOSE
    !    load old eigenvectors again
    !
    !  USES


    implicit none
    !  ARGUMENTS

    real*8, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis,n_states,n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8, dimension(n_states,n_spin,n_k_points_task) :: occ_numbers

    !  INPUTS
    !    o KS_eigenvector -- Kohn-Sham eigenvectors real format
    !    o KS_eigenvector_complex -- Kohn-Sham eigenvectors complex format
    !
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    integer:: i_k, i_basis, i_states,i_spin
    character*40 :: form1

    if (.not. eigenvectors_stored) return

    if (store_eigenvectors_to_disk_in_relaxation) then

       if (real_eigenvectors) then
          write(form1,'(A,I3,A)') '(',n_spin,'E30.18E4)'
       else
          write(form1,'(A,I3,A)') '(',2*n_spin,'E30.18E4)'
       end if
       open(file = file_name, unit = 7, status = 'unknown')

       if (real_eigenvectors) then

          do i_k = 1,n_k_points_task
             do i_basis = 1, n_basis
                do i_states = 1, n_states
                   read(7,*) (KS_eigenvector(i_basis,i_states,i_spin,i_k),i_spin = 1, n_spin)
                end do
             end do
          end do

       else ! complex eigenvectors

          do i_k = 1,n_k_points_task
             do i_basis = 1, n_basis
                do i_states = 1, n_states
                   read(unit=7,fmt=form1) (KS_eigenvector_complex(i_basis,i_states,i_spin,i_k),i_spin = 1, n_spin)
                end do
             end do
          end do
       end if

       do i_k = 1,n_k_points_task
          do i_states = 1, n_states
             read(unit=7,fmt=*) (occ_numbers(i_states,i_spin,i_k),i_spin = 1, n_spin)
          end do
       end do

       close(7)

    else

       if (use_scalapack) then

          call load_eigenvectors_scalapack

       else

          if (n_k_points_task > 0) then
             if (real_eigenvectors) then
                KS_eigenvector = stored_KS_eigenvector
             else
                KS_eigenvector_complex = stored_KS_eigenvector_complex
             end if
             occ_numbers = stored_occ_numbers
          end if
       end if
    end if

  end subroutine bfgs_load_eigenvectors
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/cleanup_relaxation
  !  NAME
  !    cleanup_relaxation
  !  SYNOPSIS

  subroutine cleanup_relaxation

    !  PURPOSE
    !    deallocate all relaxation related variables
    !
    !  USES
    implicit none
    !  ARGUMENTS
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE
    integer :: ierr

    if (allocated(constrain_relaxation)) then
       call aims_deallocate(constrain_relaxation,"constrain_relaxation")
    end if

    if (allocated(constrain_coords)) then
       call aims_deallocate(constrain_coords,"constrain_coords")
    end if

    if (allocated(constrain_lv_relaxation)) then
       call aims_deallocate(constrain_lv_relaxation,"constrain_lv_relaxation")
    end if

    if (allocated(constrain_lv_components)) then
       call aims_deallocate(constrain_lv_components,"constrain_lv_components")
    end if

    if (allocated(velocity)) then
       call aims_deallocate(velocity,"velocity")
    end if

    if (allocated(stored_coords)) then
       call aims_deallocate(stored_coords,"stored_coords")
    end if

    if (allocated(stored_forces)) then
       call aims_deallocate(stored_forces,"stored_forces")
    end if

    if (allocated(stored_DX)) then
        call aims_deallocate(stored_DX, 'stored_DX')
    end if

    if (allocated(stored_coords_2)) then
       call aims_deallocate(stored_coords_2,"stored_coords_2")
    end if

    if (allocated(stored_forces_2)) then
       call aims_deallocate(stored_forces_2,"stored_forces_2")
    end if

    if (allocated(stored_lv)) then
       call aims_deallocate(stored_lv,"stored_lv")
    end if

    if (allocated(stored_forces_lv)) then
       call aims_deallocate(stored_forces_lv,"stored_forces_lv")
    end if

    if (allocated(stored_lv_2)) then
       call aims_deallocate(stored_lv_2,"stored_lv_2")
    end if

    if (allocated(stored_forces_lv_2)) then
       call aims_deallocate(stored_forces_lv_2,"stored_forces_lv_2")
    end if

    if (allocated(stored_periodic_unit_cell_trans)) then
       call aims_deallocate(stored_periodic_unit_cell_trans,"stored_periodic_unit_cell_trans")
    end if

    if (allocated(stored_periodic_unit_cell_trans_2)) then
       call aims_deallocate(stored_periodic_unit_cell_trans_2,"stored_periodic_unit_cell_trans_2")
    end if

    if (allocated(stored_KS_eigenvector)) then
       call aims_deallocate(stored_KS_eigenvector,"stored_KS_eigenvector")
    end if

    if (allocated(stored_KS_eigenvector_complex)) then
       call aims_deallocate(stored_KS_eigenvector_complex,"stored_KS_eigenvector_complex")
    end if

    if (allocated(stored_occ_numbers)) then
       call aims_deallocate(stored_occ_numbers,"stored_occ_numbers")
    end if

    if (allocated(hessian)) then
       call aims_deallocate(hessian,"hessian")
    end if

    if (allocated(hessian_work)) then
       call aims_deallocate(hessian_work,"hessian_work")
    end if

    if (allocated(search_direction)) then
       call aims_deallocate(search_direction,"search_direction")
    end if

    if (allocated(ipiv)) then
       call aims_deallocate(ipiv,"ipiv")
    end if

    if (allocated(work)) then
       call aims_deallocate(work,"work")
    end if

    if (allocated(work_cf)) then
       call aims_deallocate(work_cf,"work_cf")
    end if

    if (allocated(mobile_lv)) then
       call aims_deallocate(mobile_lv,"mobile_lv")
    end if

    if (allocated(mobile_lv_component)) then
       call aims_deallocate(mobile_lv_component,"mobile_lv_component")
    end if

    if (allocated(translation_vectors)) then
       call aims_deallocate(translation_vectors,"translation_vectors")
    end if

    if (allocated(rotation_vectors)) then
       call aims_deallocate(rotation_vectors,"rotation_vectors")
    end if

    if (allocated(translation_vectors_lv)) then
       call aims_deallocate(translation_vectors_lv,"translation_vectors_lv")
    end if

    if (allocated(rotation_vectors_lv)) then
       call aims_deallocate(rotation_vectors_lv,"rotation_vectors_lv")
    end if

    if (allocated(mobile_atom)) then
       call aims_deallocate(mobile_atom,"mobile_atom")
    end if

    if (allocated(mobile_atom_coord)) then
       call aims_deallocate(mobile_atom_coord,"mobile_atom_coord")
    end if

    if (allocated(hess_id)) then
       call aims_deallocate(hess_id,"hess_id")
    end if

    if (use_distributed_hessian) then
       call MPI_Comm_free(comm_hess,ierr)

       if (is_worker) then
          call BLACS_Gridexit(ctxt_hess)
       end if
    end if

    if (use_symm_const_geo_relaxation) then
        call deallocate_SCR()
    end if

  end subroutine cleanup_relaxation
  !******
  !----------------------------------------------------------------------------
  !****s* relaxation/test_output_geometry
  !  NAME
  !  test_output_geometry
  !  SYNOPSIS

  subroutine test_output_geometry(lattice_vector, coords)

    !  PURPOSE
    !
    !  USES

    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !  OUTPUTS
    !    geometry
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE
    !  ARGUMENTS

    real*8, intent(IN) :: coords(3, n_atoms)
    real*8, intent(IN) :: lattice_vector(3, 3)

    integer       :: i_atom, i
    character*150 :: info_str

    write(info_str,'(2X,A)') "TEST geometry output in Angstroms"
    if (n_periodic /=3)  then
       write(info_str,'(2X,A)') "Non-periodic --> lattice vectors are dummy variables"
    end if
    call localorb_info(info_str,use_unit,'(A)')
    do i_atom = 1, n_periodic, 1  ! print row-by-row
          write(info_str,'(2X,A1,A,3(2X,F15.6))') "|", "lattice_vector",  &
               (lattice_vector(i,i_atom)*bohr ,i=1,3,1)
          call localorb_info(info_str,use_unit,'(A)')
    end do
    call localorb_info(info_str,use_unit,'(A)')
    do i_atom = 1, n_atoms, 1  ! print row-by-row
          write(info_str,'(2X,A1,A,3(2X,F15.6))') "|", "atom",  &
               (coords(i,i_atom)*bohr ,i=1,3,1)
          call localorb_info(info_str,use_unit,'(A)')
    end do

  end subroutine test_output_geometry
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/X2coords
  !  NAME
  !    X2coords
  !  SYNOPSIS

  subroutine test_output_X(X)

    !  PURPOSE
    !     Test output routine for  X(1:n_mobile_atoms) .
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: X(n_mobile_coords)
    integer            :: i

    write(info_str,'(2X,A)') "TEST X output in Angstroms"
    call localorb_info(info_str,use_unit,'(A)')
    do i = 1, n_mobile_coords, 1  ! print row-by-row
       write(info_str,'(2X,A1,A,1(2X,F15.6))') "|", "X",  &
             X(i)*bohr
       call localorb_info(info_str,use_unit,'(A)')
    end do

  end subroutine test_output_X
  !******
  !------------------------------------------------------------------------------
  !****s* relaxation/init_hess_id
  !  NAME
  !    init_hess_id
  !  SYNOPSIS

  subroutine init_hess_id()

    !  PURPOSE
    !    Create an auxiliary 1D array, hess_id, to help compute the Lindh
    !    Hessian in parallel.
    !  USES

    implicit none

    integer :: id
    integer :: i0
    integer :: i1
    integer :: l_count
    integer :: g_count
    integer :: i_atom
    integer :: i_periodic

    if (init_hess_type == HESS_LINDH) then
       call aims_allocate(hess_id, 3*(n_atoms+n_periodic), "hess_id")

       hess_id = 0
       g_count = 1
       l_count = 0

       if (is_worker) then
          i0 = offset+1
          i1 = min(offset+blk_hess,n_mobile_coords)
       else
          i0 = -1
          i1 = -1
       end if

       if (use_relaxation_constraints) then
          do i_atom = 1,n_atoms
             id = (i_atom-1)*3+1

             if (.not. constrain_relaxation(i_atom)) then
                ! x coordinate
                if (.not. constrain_coords(1,i_atom)) then
                   hess_id(id) = -g_count

                   if (-hess_id(id) >= i0 .and. -hess_id(id) <= i1) then
                      l_count = l_count+1
                      hess_id(id) = l_count
                   end if

                   g_count = g_count+1
                end if

                ! y coordinate
                if (.not. constrain_coords(2,i_atom)) then
                   hess_id(id+1) = -g_count

                   if (-hess_id(id+1) >= i0 .and. -hess_id(id+1) <= i1) then
                      l_count = l_count+1
                      hess_id(id+1) = l_count
                   end if

                   g_count = g_count+1
                end if

                ! z coordinate
                if (.not. constrain_coords(3,i_atom)) then
                   hess_id(id+2) = -g_count

                   if (-hess_id(id+2) >= i0 .and. -hess_id(id+2) <= i1) then
                      l_count = l_count+1
                      hess_id(id+2) = l_count
                   end if

                   g_count = g_count+1
                end if
             end if
          end do

          if (relax_unit_cell /= 0) then
             do i_periodic = 1,n_periodic
                id = n_atoms*3+(i_periodic-1)*3+1

                if (.not. constrain_lv_relaxation(i_periodic)) then
                   ! x coordinate
                   if (.not. constrain_lv_components(1,i_periodic)) then
                      hess_id(id) = -g_count

                      if (-hess_id(id) >= i0 .and. -hess_id(id) <= i1) then
                         l_count = l_count+1
                         hess_id(id) = l_count
                      end if

                      g_count = g_count+1
                   end if

                   ! y coordinate
                   if (.not. constrain_lv_components(2,i_periodic)) then
                      hess_id(id+1) = -g_count

                      if (-hess_id(id+1) >= i0 .and. -hess_id(id+1) <= i1) then
                         l_count = l_count+1
                         hess_id(id+1) = l_count
                      end if

                      g_count = g_count+1
                   end if

                   ! z coordinate
                   if (.not. constrain_lv_components(3,i_periodic)) then
                      hess_id(id+2) = -g_count

                      if (-hess_id(id+2) >= i0 .and. -hess_id(id+2) <= i1) then
                         l_count = l_count+1
                         hess_id(id+2) = l_count
                      end if

                      g_count = g_count+1
                   end if
                end if
             end do
          end if ! relax_unit_cell
       else ! use_relaxation_constraints
          if (is_worker) then
             do id = 1, 3*(n_atoms+n_periodic)
                hess_id(id) = -id
             end do

             do id = i0, i1
                hess_id(id) = id-offset
             end do
          end if
       end if ! use_relaxation_constraints
    end if ! HESS_LINDH

  end subroutine
  !******

  !MOL: change between symmetry reduced and full mobile index vectors
  subroutine reduce_index_vectors()
     use dimensions
     implicit none
     integer :: i_atom
     integer :: i_periodic
     integer :: i_count

     n_mobile_atoms = (SCR_n_params_coords + 2) / 3
     n_mobile_atom_coords = SCR_n_params_coords
     if (relax_unit_cell .ne. 0) then
        n_mobile_lv = (SCR_n_params_lv + 2) / 3
        n_mobile_lv_components = SCR_n_params_lv
     end if

     ! index arrays
     i_atom = 1
     mobile_atom = 0
     mobile_atom_coord = 0
     do i_count=1,n_mobile_atom_coords
       mobile_atom(i_count) = i_atom
       if (modulo(i_count,3) == 1) then
         mobile_atom_coord(i_count) = 1
       else if (modulo(i_count,3) == 2) then
         mobile_atom_coord(i_count) = 2
       else if (modulo(i_count,3) == 0) then
         mobile_atom_coord(i_count) = 3
         i_atom = i_atom + 1
       end if
     end do
     ! lattice vectors
     if (relax_unit_cell .ne. 0) then
       i_periodic = 1
       mobile_lv = 0
       mobile_lv_component = 0
       do i_count=1,n_mobile_lv_components
         mobile_lv(i_count) = i_periodic
         if (modulo(i_count,3) == 1) then
           mobile_lv_component(i_count) = 1
         else if (modulo(i_count,3) == 2) then
           mobile_lv_component(i_count) = 2
         else if (modulo(i_count,3) == 0) then
           mobile_lv_component(i_count) = 3
           i_periodic = i_periodic + 1
         end if
       end do
     end if

    n_mobile_coords = n_mobile_atom_coords + n_mobile_lv_components
  end subroutine

  subroutine full_index_vectors()
    use dimensions
    implicit none
    integer :: i_atom
    integer :: i_periodic

    n_mobile_atoms = n_atoms
    n_mobile_atom_coords = 3*n_atoms
    if (relax_unit_cell .ne. 0) then
       n_mobile_lv = n_periodic
       n_mobile_lv_components = 3*n_periodic
    end if

    do i_atom = 1, n_atoms, 1
       mobile_atom(3*(i_atom-1)+1:3*(i_atom-1)+3) = i_atom
       mobile_atom_coord(3*(i_atom-1)+1) = 1
       mobile_atom_coord(3*(i_atom-1)+2) = 2
       mobile_atom_coord(3*(i_atom-1)+3) = 3
    end do
    ! lattice vectors
    if (relax_unit_cell .ne. 0) then
       do i_periodic = 1, n_periodic, 1
          mobile_lv(3*(i_periodic-1)+1:3*(i_periodic-1)+3) = i_periodic
          mobile_lv_component(3*(i_periodic-1)+1) = 1
          mobile_lv_component(3*(i_periodic-1)+2) = 2
          mobile_lv_component(3*(i_periodic-1)+3) = 3
       end do
    end if

    n_mobile_coords = n_mobile_atom_coords + n_mobile_lv_components
  end subroutine

  ! FlK: math helpers
  !> Determinant of a 3x3 matrix, doubles
  function determinant_3x3_real(a) result(det)
    !> 3x3 matrix
    real*8, dimension(3,3), intent(in) :: a
    !> the determinant
    real*8 :: det
    !
    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
        + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function

  !> Trace of a 3x3 matrix, doubles
  function trace_3x3_real(a) result(trace)
    !> 3x3 matrix
    real*8, dimension(3,3), intent(in) :: a
    !> the trace
    real*8 :: trace
    !
    trace = (a(1,1) + a(2,2) + a(3,3)) / 3.d0
  end function

  !> Frobenius norm of a 3x3 matrix, doubles
  function frobnorm_3x3_real(a) result(frob)
    !> 3x3 matrix
    real*8, dimension(3,3), intent(in) :: a
    !> the Frobenius norm
    real*8 :: frob
    !
    frob = sqrt(trace_3x3_real(matmul(a, transpose(a))))
  end function
end module relaxation
