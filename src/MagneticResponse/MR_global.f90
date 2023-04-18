!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides the magnetic response work variables that need to be
!!  accessed globally, and the subroutines for initializing them and
!!  cleaning up.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module MR_global

  use aims_memory_tracking, only: aims_allocate
  use dimensions,           only: n_hamiltonian_matrix_size, n_my_batches, &
       & use_gga
  use geometry,             only: mr_atoms
  use grids,                only: batches
  use load_balancing,       only: batch_perm, permute_point_array, n_bp_integ
  use pbc_lists,            only: column_index_hamiltonian, index_hamiltonian
  use runtime_choices,      only: magnetic_response_pars, &
       & packed_matrix_format, PM_index, use_load_balancing, use_scalapack
  use scalapack_wrapper,    only: eigenvec, l_col, l_row, mxld, mxcol, my_col, &
       & my_row, n_my_rows, n_my_cols
  use dimensions,           only: n_basis, n_full_points, n_spin
  use mpi_tasks,            only: aims_stop, myid
  use physics,              only: hartree_potential, KS_eigenvector, &
       & occ_numbers, partition_tab, rho, rho_gradient
  use tools,                only: safe_deallocate
  use types,                only: dp

  implicit none

  private
  public :: initialize_MR, cleanup_MR, process_magnetic_keywords

  ! VARIABLES DEFINING THE CALCULATIONS' ENVIRONMENT

  ! Number of terms in the NMR Hamiltonian currently supported
  integer, parameter :: N_TERMS = 9
  character(*), parameter, public :: term_names(N_TERMS) = [character(28) :: &
       & 'Fermi contact', 'Paramagnetic spin-orbit', 'Spin-dipole', &
       & 'Diamagnetic spin-orbit', 'Paramagnetic shielding', &
       & 'Diamagnetic shielding', 'Paramagnetic magnetizability', &
       & 'Diamagnetic magnetizability', 'Electric field gradient']
  ! Information about which terms are calculated. Ordering is the same
  ! as in term_names. For instance, if the spin-dipole term was
  ! requested, which_terms(3) = .true., else .false..
  logical, public :: which_terms(N_TERMS)
  ! Information about which properties are being calculated. If
  ! any(which_terms(1:4))==.true., then j_is_active=.true. Similarly,
  ! shield_is_active=any(which_terms(5:6)),
  ! magnet_is_active=any(which_terms(7:8)), and
  ! efg_is_active=which_terms(9))
  ! These are set in check_consistency_of_keywords
  logical, public :: j_is_active, shield_is_active, magnet_is_active, &
       & efg_is_active
  ! Whether the shieldings and magnetizability are calculated with the
  ! dependence on the gauge origin (default: no)
  logical, public :: no_giao

  ! MAIN WORK VARIABLES

  ! The magnetic shielding tensors (3,3,size(active_nuclei))
  real(dp), allocatable, public :: shield_para_tensor(:,:,:), &
       & shield_dia_tensor(:,:,:)
  ! J-coupling tensors (3,3,size(active_nuclei),size(active_nuclei))
  real(dp), dimension(:,:,:,:), allocatable, public :: &
       & PSO_tensor, FC_tensor, spin_dipole_tensor, DSO_tensor
  ! Dipolar tensor
  real(dp), allocatable, public :: dipolar_tensor(:,:,:,:)
  ! Magnetizability tensors.
  real(dp), public :: mag_para_tensor(3,3), mag_dia_tensor(3,3)
  ! Electric field gradient
  real(dp), allocatable, public :: EFG_tensor(:,:,:)
  ! Whether to calculated the xx,yy,zz components of the spin-dipole term
  logical, public :: calc_spin_dipole_diag
  ! Whether to calculate the full tensor for the second-order response
  ! or only the diagonal elements.
  logical, public :: calc_full_tensor = .false.
  ! Indices of NMR active nuclei
  integer, allocatable, public :: active_nuclei(:)

  ! OTHER

  ! Gauge origin (default is center of mass, else mr_gauge_origin)
  real(dp), public :: gauge_origin(3)
  ! For printing output
  character(14), allocatable, public :: c_atoms(:)
  ! Highest occupied state for each spin. This is filled in
  ! initialize_MR.
  integer, public :: n_occ_states(2)
  ! Row and column indices of the highest occupied state in a 2D block
  ! cyclic array for each spin.
  integer, public :: my_max_occ_col(2), my_max_occ_row(2)
  ! Timers for the non-self-consistent parts of the first-order
  ! Hamiltonian
  real(dp), public :: timers(4,N_TERMS)
  ! Memory peaks for the whole of MR part and DFPT separately.
  integer*8, public :: max_total_memory, max_dfpt_memory
  ! Total wall time spent on pdgemms or zgemms
  integer, public :: walltime_mul
  ! Total time spent on calls to MR_core
  integer, public :: walltime_all

contains

  !!  FUNCTION
  !!
  !!  Before this subroutine is called, the magnetic response input
  !!  must be extracted from control.in and geometry.in in a raw
  !!  format. Now process that input to determine what are the
  !!  specific values of some of the essential keywords for magnetic
  !!  response. If there are any inconsistencies in the keywords, this
  !!  should be revealed in a call to check_consistency_of_keywords
  !!  before any calculations start.
  !!
  subroutine process_magnetic_keywords()
    character(*), parameter :: &
         & THIS_SUB = 'MR_global::process_magnetic_keywords::'
    character(:), allocatable :: string, word
    integer :: i_atom, j_atom
    which_terms = .false. ! Default is actually .true.; see below
    no_giao  = .false.
    string = adjustl(trim(magnetic_response_pars))
    ! Determine which terms need to be calculated.
    do while (index(string,' ') > 0)
       string = adjustl(trim(string(index(string,' '):)))
       word = string(:index(string,' '))
       if (len(word) == 0) word = string
       select case(word)
       case('J_coupling', 'J', 'j')
          which_terms(1:4) = .true.
       case('fc')
          which_terms(1) = .true.
       case('po')
          which_terms(2) = .true.
       case('sd')
          which_terms(3) = .true.
       case('do')
          which_terms(4) = .true.
       case('shielding', 's')
          which_terms(5:6) = .true.
       case('shielding_p', 's_p')
          which_terms(5) = .true.
       case('shielding_d', 's_d')
          which_terms(6) = .true.
       case('magnet', 'm')
          which_terms(7:8) = .true.
       case('magnet_p', 'm_p')
          which_terms(7) = .true.
       case('magnet_d', 'm_d')
          which_terms(8) = .true.
       case('efg')
          which_terms(9) = .true.
       case('full')
          calc_full_tensor = .true.
       case('no_giao')
          no_giao = .true.
       case default
          if (myid == 0 .and. len(trim(word)) > 0) &
               & call aims_stop('Unsupported keyword following &
               &magnetic_response: '//word, THIS_SUB)
       end select
    end do
    if (.not. any(which_terms)) which_terms = .true.
    ! Next determine which atoms in geometry.in are NMR active.
    call aims_allocate(active_nuclei, count(mr_atoms), &
         & THIS_SUB//'active_nuclei')
    j_atom = 1
    do i_atom = 1, size(mr_atoms)
       if (mr_atoms(i_atom)) then
          active_nuclei(j_atom) = i_atom
          j_atom = j_atom + 1
       end if
    end do
    j_is_active = any(which_terms(1:4))
    shield_is_active = any(which_terms(5:6))
    magnet_is_active = any(which_terms(7:8))
    efg_is_active = which_terms(9)
  end subroutine process_magnetic_keywords

  !!  FUNCTION
  !!
  !!  Allocates and initializes the main MR work variables. This
  !!  subroutine assumes we already know the number of NMR active
  !!  nuclei, i.e., active_nuclei must have been allocated and
  !!  initialized (see MR_main).
  !!
  subroutine initialize_MR()
    character(*), parameter :: THIS_SUB = 'MR_global::initialize_MR::'
    ! The MR module works with spin-polarized quantities only. If the
    ! initial nonperturbative calculation was non-spin-polarized, copy
    ! the values from the first dimension of rho/rho_gradient to the
    ! second.
    real(dp), allocatable :: rho_tmp(:), rho_gradient_tmp(:,:)
    integer :: i_row, i_col, i_counter, i_spin
    max_dfpt_memory = 0
    ! The MR module always assumes we are using Scalapack and packed
    ! matrices. If necessary, manually allocate and initialize
    ! quantities such as my_row and index_hamiltonian here and
    ! deallocate them at the end of the calculation.
    if (.not. use_scalapack) then
       call aims_allocate(eigenvec, n_basis, n_basis, n_spin, &
            & THIS_SUB//'eigenvec')
       call aims_allocate(l_row, n_basis, THIS_SUB//'l_row')
       call aims_allocate(l_col, n_basis, THIS_SUB//'n_basis')
       call aims_allocate(my_row, n_basis, THIS_SUB//'my_row')
       call aims_allocate(my_col, n_basis, THIS_SUB//'my_col')
       eigenvec = KS_eigenvector(:,:,:,1)
       mxld = n_basis
       n_my_rows = n_basis
       mxcol = n_basis
       n_my_cols = n_basis
       l_row = [(i_counter, i_counter=1,mxld)]
       l_col = l_row
       my_row = l_row
       my_col = l_col
    end if
    if (packed_matrix_format /= PM_index) then
       call aims_allocate(index_hamiltonian, 2, 1, n_basis, &
            & THIS_SUB//'index_hamiltonian')
       call aims_allocate(column_index_hamiltonian, n_basis**2, &
            & THIS_SUB//'column_index_hamiltonian')
       index_hamiltonian(:,1,1) = [1,1]
       column_index_hamiltonian(1) = 1
       i_counter = 1
       do i_row = 2, n_basis
          index_hamiltonian(1,1,i_row) = index_hamiltonian(2,1,i_row-1)+1
          index_hamiltonian(2,1,i_row) = index_hamiltonian(1,1,i_row) + &
               & i_counter
          i_counter = i_counter + 1
          column_index_hamiltonian(index_hamiltonian(1,1,i_row): &
                 & index_hamiltonian(1,1,i_row)+i_counter) = &
               & [(i_col, i_col=1,i_counter)]
       end do
    end if
    if (.not. use_load_balancing) then
       batch_perm(n_bp_integ)%n_local_matrix_size = n_hamiltonian_matrix_size
       batch_perm(n_bp_integ)%n_my_batches = n_my_batches
       batch_perm(n_bp_integ)%n_basis_local = n_basis
       batch_perm(n_bp_integ)%batches => batches
       batch_perm(n_bp_integ)%partition_tab => partition_tab
       batch_perm(n_bp_integ)%n_full_points = n_full_points
    end if
    if (n_spin == 1) then
       call aims_allocate(rho_tmp, size(rho,2), THIS_SUB//'rho_tmp')
       rho_tmp = rho(1,:)
       call safe_deallocate(rho)
       call aims_allocate(rho, 2, size(rho_tmp), THIS_SUB//'rho')
       rho = spread(rho_tmp,1,2)/2
       call safe_deallocate(rho_tmp)
       if (use_gga) then
          call aims_allocate(rho_gradient_tmp, 3, size(rho_gradient,3), &
               & THIS_SUB//'rho_gradient_tmp')
          rho_gradient_tmp = rho_gradient(:,1,:)
          call safe_deallocate(rho_gradient)
          call aims_allocate(rho_gradient, 3, 2, size(rho_gradient_tmp,2), &
               & THIS_SUB//'rho_gradient')
          rho_gradient = spread(rho_gradient_tmp,2,2)/2
          call safe_deallocate(rho_gradient_tmp)
       end if
    end if
    ! Allocate and initialize tensors for the second-order
    ! response. These are all very small. Just allocate all of them
    ! even if not used.
    if (shield_is_active) then
       call aims_allocate(shield_para_tensor, 3, 3, size(active_nuclei), &
            & THIS_SUB//'shield_para_tensor')
       shield_para_tensor = 0d0
       call aims_allocate(shield_dia_tensor, 3, 3, size(active_nuclei), &
            & THIS_SUB//'shield_dia_tensor')
       shield_dia_tensor = 0d0
    end if
    if (j_is_active) then
       call aims_allocate(FC_tensor, 3, 3, size(active_nuclei), &
            & size(active_nuclei), THIS_SUB//'FC_tensor')
       FC_tensor = 0d0
       call aims_allocate(PSO_tensor, 3, 3, size(active_nuclei), &
            & size(active_nuclei), THIS_SUB//'PSO_tensor')
       PSO_tensor = 0d0
       call aims_allocate(spin_dipole_tensor, 3, 3, size(active_nuclei), &
            & size(active_nuclei), THIS_SUB//'spin_dipole_tensor')
       spin_dipole_tensor = 0d0
       call aims_allocate(DSO_tensor, 3, 3, size(active_nuclei), &
            & size(active_nuclei), THIS_SUB//'DSO_tensor')
       DSO_tensor = 0d0
       call aims_allocate(dipolar_tensor, 3, 3, size(active_nuclei), &
            & size(active_nuclei), THIS_SUB//'dipolar_tensor')
       dipolar_tensor = 0d0
    end if
    if (efg_is_active) then
       call aims_allocate(EFG_tensor, 3, 3, size(active_nuclei), 'EFG_tensor')
       EFG_tensor = 0d0
    end if
    mag_para_tensor = 0d0
    mag_dia_tensor = 0d0
    ! Reset timers
    timers = 0d0
    walltime_mul = 0
    ! Find the number of occupied states for each spin channel
    do i_spin = 1, n_spin
       n_occ_states(i_spin) = minloc(occ_numbers(:,i_spin,1), 1, &
            & occ_numbers(:,i_spin,1) <= epsilon(1d0)) - 1
    end do
    if (n_spin == 1) n_occ_states(2) = n_occ_states(1)
    ! Find the corresponding indices in 2D block cyclic arrays.
    if (use_scalapack) then
       do i_spin = 1, 2
          my_max_occ_row(i_spin) = &
               & maxloc(my_row, 1, my_row <= n_occ_states(i_spin))
          my_max_occ_col(i_spin) = &
               & maxloc(my_col, 1, my_col <= n_occ_states(i_spin))
          if (my_row(1) > n_occ_states(i_spin)) my_max_occ_row(i_spin) = 0
          if (my_col(1) > n_occ_states(i_spin)) my_max_occ_col(i_spin) = 0
       end do
    else
       my_max_occ_row = n_occ_states
       my_max_occ_col = n_occ_states
    end if
  end subroutine initialize_MR

  subroutine cleanup_MR()
    character(*), parameter :: THIS_SUB = 'MR_global::cleanup_MR::'
    real(dp), allocatable :: rho_tmp(:), rho_gradient_tmp(:,:)
    ! If size(rho,1)==1 prior to calling MR_main(), restore
    ! rho/rho_gradient to its initial shape.
    if (n_spin == 1) then
       call aims_allocate(rho_tmp, size(rho,2), THIS_SUB//'rho_tmp')
       rho_tmp = rho(1,:)
       call safe_deallocate(rho)
       call aims_allocate(rho, 1, size(rho_tmp), THIS_SUB//'rho')
       rho(1,:) = 2*rho_tmp
       call safe_deallocate(rho_tmp)
       if (use_gga) then
          call aims_allocate(rho_gradient_tmp, 3, size(rho_gradient,3), &
              & THIS_SUB//'rho_gradient_tmp')
          rho_gradient_tmp = rho_gradient(:,1,:)
          call safe_deallocate(rho_gradient)
          call aims_allocate(rho_gradient, 3, 1, size(rho_gradient_tmp,2), &
              & THIS_SUB//'rho_gradient')
          rho_gradient(:,1,:) = 2*rho_gradient_tmp
          call safe_deallocate(rho_gradient_tmp)
       end if
    end if
    ! WPH: Since the integration hijacks the load balancing infrastructure, we
    !      need to undo the changes to avoid possible issues down the line
    if (.not. use_load_balancing) then
       batch_perm(n_bp_integ)%n_local_matrix_size = 0
       batch_perm(n_bp_integ)%n_my_batches = 0
       batch_perm(n_bp_integ)%n_basis_local = 0
       nullify(batch_perm(n_bp_integ)%batches)
       nullify(batch_perm(n_bp_integ)%partition_tab)
       batch_perm(n_bp_integ)%n_full_points = 0
    end if
    ! If arrays related to Scalapack and/or packed matrices were
    ! manually allocated in MR_initialize, deallocate them here.
    if (.not. use_scalapack) then
       mxld = 0
       n_my_rows = 0
       mxcol = 0
       n_my_cols = 0
       call safe_deallocate(eigenvec)
       call safe_deallocate(l_row)
       call safe_deallocate(l_col)
       call safe_deallocate(my_row)
       call safe_deallocate(my_col)
    end if
    if (packed_matrix_format /= PM_index) then
       call safe_deallocate(column_index_hamiltonian)
       call safe_deallocate(index_hamiltonian)
    end if
    if (allocated(c_atoms)) deallocate(c_atoms)
    call safe_deallocate(active_nuclei)
    call safe_deallocate(FC_tensor)
    call safe_deallocate(PSO_tensor)
    call safe_deallocate(spin_dipole_tensor)
    call safe_deallocate(DSO_tensor)
    call safe_deallocate(dipolar_tensor)
    call safe_deallocate(shield_para_tensor)
    call safe_deallocate(shield_dia_tensor)
    call safe_deallocate(EFG_tensor)
  end subroutine cleanup_MR
end module MR_global
