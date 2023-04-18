!****f* FHI-aims/calculate_second_variational_soc
!*  NAME
!*    calculate_second_variational_soc
!*  SYNOPSIS
subroutine calculate_second_variational_soc
!*  PURPOSE
!*    Parent subroutine for post-processing method to include spin-orbit
!*    coupling in a second-variational manner on the self-consistent Kohn-Sham
!*    eigenfunctions, as previously calculated by the scf cycle.
!*    This should be called before using any of the spin-orbit coupling
!*    functionality in FHI-aims, as it sets up much of the underlying
!*    architecture.
!*  USES
  use localorb_io, only: localorb_info, use_unit, default_unit, OL_norm
  use applicable_citations, only: cite_reference
  use dimensions, only: flag_out_dipmat, &
      gap_for_min_energy_in_soc, gap_for_min_energy_in_soc_set, &
      gap_for_saved_min_energy_in_soc, max_energy_include_in_soc, &
      max_energy_include_in_soc_set, max_energy_save_in_soc, &
      min_energy_include_in_soc, min_energy_save_in_soc, &
      gap_for_saved_min_energy_in_soc_set, max_energy_save_in_soc_set, &
      min_energy_include_in_soc_set, min_energy_save_in_soc_set, &
      n_k_points_task, n_basis, n_hamiltonian_matrix_size, n_k_points, &
      n_periodic, n_spin, n_states, n_write_soc_eigenvectors, &
      save_soc_perturbed_eigenvectors, force_occupation_projector
  use runtime_choices, only: calculate_work_function, occupation_type, out_dos,&
      out_k_points_eigenvalues, pert_dos_on, real_eigenvectors, &
      use_dipole_correction, use_scalapack, use_symmetry_reduced_k_grid, &
      write_soc_eigenvectors, use_local_index, use_load_balancing, &
      out_soc_eigenvalues
  use load_balancing, only: use_batch_permutation, batch_perm, n_bp_integ, &
      permute_point_array
  use constants, only: hartree
  use mpi_tasks, only: myid, n_tasks, mpi_comm_global, aims_stop_coll
  use synchronize_mpi_basic, only: sync_vector
  use physics, only: hartree_partition_tab, hartree_potential, partition_tab, &
      KS_eigenvector,KS_eigenvector_complex, KS_eigenvector_soc_perturbed, &
      KS_eigenvalue, KS_eigenvalue_soc_perturbed, n_electrons, &
      occ_numbers, occ_numbers_soc, chemical_potential_soc, &
      total_energy, soc_non_sc_total_energy, entropy_correction, ev_sum, &
      rho
  use timing, only: get_times, get_timestamps, tot_time_soc, tot_clock_time_soc
  use scalapack_wrapper, only: my_scalapack_id, eigenvec, eigenvec_complex, &
      mxld, mxcol, my_k_point
  use pbc_lists, only: k_weights, k_point_list
  use aims_memory_tracking, only: aims_mem_current_output, aims_allocate, &
      aims_deallocate
  use soc_utilities, only: calculate_spin_expectation, &
      find_num_core_states_from_energy, find_min_energy_from_gap, &
      find_num_high_states_from_energy, convert_wf_basis_to_compute_basis, &
      perform_soc_perturbation, write_soc_values, &
      write_soc_perturbed_eigenvectors, convert_sr_to_soc_environment, &
      revert_soc_to_sr_environment, create_sorted_sr_eigenvalues_list, &
      create_sorted_sr_eigenvalues_idxmap
  use scalapack_soc          ! We really do use almost everything here, so...
  use dimensions_soc, only: n_states_sr, sr_state_start, n_states_soc, &
      n_high_states_omit_from_soc, n_core_states_omit_from_soc, &
      n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll, &
      n_saved_states_soc, soc_saved_state_start
  use species_data, only: l_shell_max
  use force_occupation, only: get_occupation_numbers_occ_p0, adjust_force_occ_sr
  implicit none
!*  ARGUMENTS
!*    None (called from main(), uses global variables)
!*  INPUTS
!*    None (uses global variables)
!*  OUTPUTS
!*    None (uses global variables)
!*  AUTHORS
!*    William Huhn (Duke University), based off cluster code written by Matthias
!*    Gramzow
!*  NOTES
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*
!*    When calling this subroutine, it is assumed that steps 1+2 in
!*      Section III.3 from Huhn and Blum, Phys. Rev. Mater. (2017) have
!*      already been performed.
!*  TODO
!*    o Create new packing for soc_matrix instead of relying on the real-space
!*      Hamiltonian's packing (which is much less sparse than soc_matrix)
!*    o Better sorting algorithm than bubblesort...
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject
!*    the terms and conditions of the respective license agreement."
!*  SOURCE

! WPH:   This is the main wrapper function to include the spin-orbit interaction
!        in FHI-aims. It requires the non-relativistic self-consistent KS
!        eigenfunctions be calculated, and self-consistency of the electron
!        density will be lost (as the method is single-shot.)
!
!        To understand what is going on here, please read the Phys. Rev.
!        Mater. paper, where the steps are documented.

  ! Computed quantities
  real*8,dimension(:,:),      allocatable :: soc_matrix  ! The real-space SOC matrix, analogous to hamiltonian in the main code
                                                         ! In the Phys Rev Materials paper, these are the \Pi matrices
  complex*16, dimension(:,:), allocatable :: soc_ham     ! The matrix elements of the (full) SOC operator.  This can either be a
                                                         ! (2n_states,2n_states) or (mxld_soc, mcol_soc) sized matrix, depending
                                                         ! on whether LAPACK or ScaLAPACK is used
  ! The expectation values for five spin operators, in order:
  ! 1.   Projector onto spin-up component in z-axis
  ! 2.   Projector onto spin-down component in z-axis
  ! 3-5. The three components of the Pauli matrix vector (sigma_x, sigma_y, sigma_z)
  ! 4 October 2017:  When moving the code to fully support BLACS, I made the decision to temporarily remove calculation and output
  !                  of these variables, but they're still allocated in the code (and passed around in interfaces) because they may
  !                  be useful someday.
  real*8, dimension(:,:,:), allocatable   :: spin_expectation

  ! Total energy variables
  ! Known to give non-physical results, because we calculate only the correction due to sum-of-eigenvalues.  It can be argued that
  ! it should be removed entirely, but we decided to keep it in specifically so that people know that the "real" total energy doesn't
  ! include SOC effects and that we can clearly state in the output that they shouldn't use these values
  real*8  ::  ev_sum_new, total_energy_corrected_loc, total_energy_loc

  ! Parallelization variables
  integer, dimension(n_k_points_task) :: my_k_points ! Which k-points are assigned to this task (this variable is like
                                                     ! my_k_point for ScaLAPACK calculations, but for the more general
                                                     ! case where the MPI rank can have more than one k-point)

  ! Used to re-introduce the notion of spin channels back into non-spin-polarized calculations
  integer :: spin_index

  ! Used for outputting information to stdout
  character*300 :: info_str

  ! debug
  integer :: info
  integer, dimension(3) :: max_soc_matrix_ind
  real*8, dimension(3) :: max_soc_matrix_val
  real*8 :: this_perturb_matrix_max
  integer :: pauli_coup
  real*8, dimension(3) :: temp
  character l_to_str ! Function found in convert_l_str.f90
  character l_char_1, l_char_2
  character*9 :: k_point_string

  ! Counters
  integer ::  i_k_point, this_k_point, i_state, j_state, i_spin, i_point, i_basis
  integer ::  i, j, k
  integer :: i_basis_1, i_basis_2, i_fn_1, i_fn_2, i_coord, i_cell, i_index_real, i_size, i_subspace, temp_int
  integer :: saved_state_offset

  ! timings
  real*8 :: time_soc, clock_time_soc
  real*8 :: time_matrix, clock_time_matrix, tot_time_matrix = 0.0d0, tot_clock_time_matrix = 0.0d0
  real*8 :: time_ham, clock_time_ham, tot_time_ham = 0.0d0, tot_clock_time_ham = 0.0d0
  real*8 :: time_diag_soc, clock_time_diag_soc, tot_time_diag_soc = 0.0d0, tot_clock_time_diag_soc = 0.0d0
  real*8 :: time_eigenvec, clock_time_eigenvec, tot_time_eigenvec = 0.0d0, tot_clock_time_eigenvec = 0.0d0
  real*8 :: time_sync, clock_time_sync, tot_time_sync = 0.0d0, tot_clock_time_sync = 0.0d0
  real*8 :: time_post, clock_time_post, tot_time_post = 0.0d0, tot_clock_time_post = 0.0d0

  ! Energy window temporary variable
  integer :: n_core_states_omit_from_saved_soc = 0
  integer :: n_high_states_omit_from_saved_soc = 0
  real*8, dimension(:), allocatable     :: KS_eigenvalue_soc_perturbed_temp
  integer, dimension(:), allocatable    :: sr_to_soc_idxmap
  real*8, dimension(:,:), allocatable   :: spin_expectation_temp

  ! Occupation number determination
  logical :: t_out = .true.
  real*8  :: chemical_potential_spin_soc(1) ! A trivial array for interfacing with output_energy_levels

  ! Used for calculating the second-variational eigenvectors from the diagonalization output
  complex*16, dimension(:,:), allocatable :: eigenvec_soc_wf_basis

  ! Variables for self-consistency SOC (currently broken!)
  logical :: soc_converged, soc_enough_walltime_left
  real*8 :: old_total_energy, old_entropy_correction

  ! Variables for determining sparsity of soc_matrix
  real*8,       parameter :: soc_mat_rel_thresh = 1.0d-15  ! Relative threshhold for determining sparsity in soc_matrix
  real*8                  :: soc_mat_thresh                ! Actual threshhold

  ! Variables for load balancing
  integer                 :: ld_soc_matrix

  character(*), parameter :: func = 'calculate_second_variational_soc'

  ! Variable that tells whether the density is converged or not. Since this is post-processing, it is for now
  ! set to .true. always, but could become an argument that is set depending on the case in the furure
  logical :: conv
  conv=.true.

  ! Making sure I get paid before anything gets done.  Yes, even before starting
  ! the timer.
  call cite_reference("Spin_Orbit_Coupling")

  call get_timestamps( time_soc, clock_time_soc )

  ! SOC header information
  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **          STARTING SECOND VARIATIONAL SOC CALCULATION         **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **                     Stable functionality:                    **"
  call localorb_info( info_str )
  write(info_str, *)  " **                 Writes eigenvalues to stdout                 **"
  call localorb_info( info_str )
  write(info_str, *)  " **                    Recalculates band gaps                    **"
  call localorb_info( info_str )
  write(info_str, *)  " **                Regular and interpolated DOSs                 **"
  call localorb_info( info_str )
  write(info_str, *)  " **                       Band structures                        **"
  call localorb_info( info_str )
  write(info_str, *)  " **                      Mulliken analysis                       **"
  call localorb_info( info_str )
  write(info_str, *)  " **                         ELSI Support                         **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **               New/Experimental functionality:                **"
  call localorb_info( info_str )
  write(info_str, *)  " **         use_local_index and load balancing support           **"
  call localorb_info( info_str )
  write(info_str, *)  " **         Dielectric constant/absorption coefficents           **"
  call localorb_info( info_str )
  write(info_str, *)  " **             Specifying restricted energy window              **"
  call localorb_info( info_str )
  write(info_str, *)  " **                         Cube files                           **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **            Known buggy functionality (in SOC code):          **"
  call localorb_info( info_str )
  write(info_str, *)  " **   SOC-perturbed total energies are meaningless (DO NOT USE)  **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **                             Note:                            **"
  call localorb_info( info_str )
  write(info_str, *)  " **             You may see oscillations or blips in             **"
  call localorb_info( info_str )
  write(info_str, *)  " **       high-lying conduction bands or slight sub-banding      **"
  call localorb_info( info_str )
  write(info_str, *)  " **            in calculations.  This is likely due to           **"
  call localorb_info( info_str )
  write(info_str, *)  " **        to not enough empty bands being mixed into the        **"
  call localorb_info( info_str )
  write(info_str, *)  " **             eigensolver, and you should increase             **"
  call localorb_info( info_str )
  write(info_str, *)  " **      the keyword 'empty_states' until the behavior goes      **"
  call localorb_info( info_str )
  write(info_str, *)  " **                        away to fix it.                       **"
  call localorb_info( info_str )
  write(info_str, *)  " **      If you are feeling particularly paranoid, consider      **"
  call localorb_info( info_str )
  write(info_str, *)  " **        setting the calculate_all_eigenstates flag in         **"
  call localorb_info( info_str )
  write(info_str, *)  " **      control.in to calculate all empty states possible.      **"
  call localorb_info( info_str )
  write(info_str, *)  " **                    This is rarely needed.                    **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )
  write(info_str, *)  " **                      Future Developments:                    **"
  call localorb_info( info_str )
  write(info_str, *)  " **                       Self-Consistency                       **"
  call localorb_info( info_str )
  write(info_str, *)  " ******************************************************************"
  call localorb_info( info_str )

  ! TODO:  WARNING IF ATOM TOO HEAVY

  if (use_symmetry_reduced_k_grid.and.n_spin .gt. 1) then
    write(info_str, *)
    call localorb_info( info_str )
    write(info_str, *)  " ******************************************************************"
    call localorb_info( info_str )
    write(info_str, *)  " **                          NOTE                                **"
    call localorb_info( info_str )
    write(info_str, *)  " **    YOU CURRENTLY HAVE SYMMETRY REDUCED K_GRID TURNED ON      **"
    call localorb_info( info_str )
    write(info_str, *)  " **       AND ARE RUNNING A SPIN POLARIZED CALCULATION           **"
    call localorb_info( info_str )
    write(info_str, *)  " **                 WITH SPIN-ORBIT COUPLING                     **"
    call localorb_info( info_str )
    write(info_str, *)  " ******************************************************************"
    call localorb_info( info_str )
  end if

  call mpi_barrier(mpi_comm_global, info)

  ! Set up the energy window in which to perform second-variational SOC
  ! This determines the number of states in the SOC Hamiltonian matrix
  ! REQUIREMENTS:
  !   - Energy window must include Fermi level (omitted core states must be
  !     fully occupied, omitted high-lying valence must be unoccupied)
  !   - Energy window must include same number of states in both spin channels
  !   - Energy window must include same number of states for every k-point
  write(info_str, *)
  call localorb_info(info_str)
  ! Inform user which of their choices will take priority
  if (gap_for_min_energy_in_soc_set.and.(min_energy_include_in_soc_set.or.n_core_states_omit_from_soc.gt.0)) then
    write (info_str, '(2X,A)') &
         "You have specified multiple methods for determining the lower bound of the energy"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "window used in SOC; out of the choices you made, the automated selection based on"
    call localorb_info( info_str )
    write(info_str, '(2X,A)')  &
         "the gap criterion will take priority."
    call localorb_info( info_str )
  else if (min_energy_include_in_soc_set.and.n_core_states_omit_from_soc.gt.0) then
    write (info_str, '(2X,A)') &
         "You have specified multiple methods for determining the lower bound of the energy"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "made, the minimum energy eigenvalue to include will take priority."
    call localorb_info( info_str )
  end if
  if (max_energy_include_in_soc_set.and.n_high_states_omit_from_soc.gt.0) then
    write (info_str, '(2X,A)') &
         "You have specified multiple methods for determining the upper bound of the energy"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "window used in SOC; out of the choices you made, the maximum energy eigenvalue "
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "to include will take priority."
    call localorb_info( info_str )
  end if

  ! Now determine the bounds of the energy window
  if (gap_for_min_energy_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') "Requested gap for setting min energy in SOC  : ", gap_for_min_energy_in_soc, " eV"
    call localorb_info( info_str )

    call find_min_energy_from_gap( gap_for_min_energy_in_soc, &
         KS_eigenvalue, occ_numbers, min_energy_include_in_soc )
  end if
  if (min_energy_include_in_soc_set.or.gap_for_min_energy_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') "Minimum energy eigenvalue to include in SOC  : ", min_energy_include_in_soc, " eV"
    call localorb_info( info_str )

    call find_num_core_states_from_energy( min_energy_include_in_soc, &
         KS_eigenvalue, occ_numbers, n_core_states_omit_from_soc )
  end if
  if (max_energy_include_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') "Maximum energy eigenvalue to include in SOC  : ", max_energy_include_in_soc, " eV"
    call localorb_info( info_str )

    call find_num_high_states_from_energy( max_energy_include_in_soc, &
         KS_eigenvalue, occ_numbers, n_high_states_omit_from_soc )
  end if
  write (info_str, '(2X,A,I9)') "Number of core states to omit from SOC       :         ", n_core_states_omit_from_soc
  call localorb_info( info_str )
  if (mod(n_core_states_omit_from_soc,2) .eq. 1) then
    write (info_str, '(2X,A)' ) "* The number of core states omitted from SOC must be even to ensure both spin channels have&
         & the same number and characters of states (necessary but not sufficient!)  Exiting."
    call aims_stop_coll( info_str )
  end if
  if (dble(n_core_states_omit_from_soc) .ge. n_electrons) then
    write(info_str,'(1X,A)') '* Specified number of core states to omit from SOC&
         & would eliminate all occupied states from calculation, preventing&
         & determination of Fermi level.  Exiting.'
    call aims_stop_coll( info_str )
  end if
  write (info_str, '(2X,A,I9)') "Number of high-lying states to omit from SOC :         ", n_high_states_omit_from_soc
  call localorb_info( info_str )
  if (mod(n_high_states_omit_from_soc,2) .eq. 1) then
    write (info_str, '(2X,A)' ) "* The number of high-lying states omitted from SOC must be even to ensure both spin channels&
         & have the same number and characters of states (necessary but not sufficient!)  Exiting."
    call aims_stop_coll( info_str )
  end if
  if (dble(2*n_states - n_high_states_omit_from_soc) .le. n_electrons) then
    write(info_str,'(1X,A)') '* Specified number of high-lying states to omit&
         & from SOC would eliminate occupied states from calculation, preventing&
         & determination of Fermi level.  Exiting.'
    call aims_stop_coll( info_str )
  end if
  write(info_str, *)
  call localorb_info( info_str )

  ! Set up energy window for saving results of second-variational SOC
  ! This determines the number of states in the eigenvalues/eigenvector matrices
  ! This will be a subset of the previous energy window
  ! As of this writing, there isn't a way to manually set indices to include/omit states

  ! Inform user which of their choices will take priority
  if (gap_for_saved_min_energy_in_soc_set.and.(min_energy_save_in_soc_set.or.soc_saved_state_start.gt.1)) then
    write (info_str, '(2X,A)') &
         "You have specified multiple methods for determining the lower bound of the energy"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "window to save in SOC; out of the choices you made, the automated selection based"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "on the gap criterion will take priority."
    call localorb_info( info_str )
  else if (min_energy_save_in_soc_set.and.soc_saved_state_start.gt.1) then
    write (info_str, '(2X,A)') &
         "You have specified multiple methods for determining the lower bound of the energy"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "window to save in SOC; out of the choices you made, the minimum energy eigenvalue"
    call localorb_info( info_str )
    write (info_str, '(2X,A)') &
         "to include will take priority."
    call localorb_info( info_str )
  end if
  ! Now determine the bounds of the energy window
  if (gap_for_saved_min_energy_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') &
         "Requested gap for setting min energy to be saved in SOC  : ", gap_for_saved_min_energy_in_soc, " eV"
    call localorb_info( info_str )

    call find_min_energy_from_gap( gap_for_saved_min_energy_in_soc, &
         KS_eigenvalue, occ_numbers, min_energy_save_in_soc )
  end if
  if (min_energy_save_in_soc_set.or.gap_for_saved_min_energy_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') "Minimum energy eigenvalue to save in SOC     : ", min_energy_save_in_soc, " eV"
    call localorb_info( info_str )

    call find_num_core_states_from_energy( min_energy_save_in_soc, &
         KS_eigenvalue, occ_numbers, n_core_states_omit_from_saved_soc )
  end if
  if (max_energy_save_in_soc_set) then
    write (info_str, '(2X,A,F17.5,A)') "Maximum energy eigenvalue to save in SOC     : ", max_energy_save_in_soc, " eV"
    call localorb_info( info_str )

    call find_num_high_states_from_energy( max_energy_save_in_soc, &
         KS_eigenvalue, occ_numbers, n_high_states_omit_from_saved_soc )
  end if
  ! Make sure that we're not saving more states than we've calculated
  if (n_core_states_omit_from_saved_soc < n_core_states_omit_from_soc) then
    n_core_states_omit_from_saved_soc = n_core_states_omit_from_soc
  end if
  if (n_high_states_omit_from_saved_soc < n_high_states_omit_from_soc) then
    n_high_states_omit_from_saved_soc = n_high_states_omit_from_soc
  end if

  write (info_str, '(2X,A,I9)') "Index of first SOC-perturbed state to save   :         ", n_core_states_omit_from_saved_soc+1
  call localorb_info( info_str )
  if (mod(n_core_states_omit_from_soc,2) .eq. 1) then
    write (info_str, '(2X,A)' ) "* The number of core states omitted from SOC must be even to ensure both spin channels have&
         & the same number and characters of states (necessary but not sufficient!)  Exiting."
    call aims_stop_coll( info_str )
  end if
  if (dble(n_core_states_omit_from_soc) .ge. n_electrons) then
    write(info_str,'(1X,A)') '* Specified number of core states to omit from SOC&
         & would eliminate all occupied states from calculation, preventing&
         & determination of Fermi level.  Exiting.'
    call aims_stop_coll( info_str )
  end if

  write (info_str, '(2X,A,I9)') "Index of last SOC-perturbed state to save    :         ", 2*n_states - n_high_states_omit_from_saved_soc
  call localorb_info( info_str )
  if (mod(n_high_states_omit_from_soc,2) .eq. 1) then
    write (info_str, '(2X,A)' ) "* The number of high-lying states omitted from SOC must be even to ensure both spin channels&
         & have the same number and characters of states (necessary but not sufficient!)  Exiting."
    call aims_stop_coll( info_str )
  end if
  if (dble(2*n_states - n_high_states_omit_from_soc) .le. n_electrons) then
    write(info_str,'(1X,A)') '* Specified number of high-lying states to omit&
         & from SOC would eliminate occupied states from calculation, preventing&
         & determination of Fermi level.  Exiting.'
    call aims_stop_coll( info_str )
  end if
  write(info_str, *)
  call localorb_info( info_str )

  ! Set up various indexing variables (see dimensions_soc for more details)
  ! Move into its own subroutine?
  n_states_soc          = 2*n_states - n_core_states_omit_from_soc &
                         - n_high_states_omit_from_soc
  n_states_sr           = n_states_soc/2
  sr_state_start        = n_core_states_omit_from_soc/2 + 1 ! Starting index for SR states to include
  n_basis_soc_coll      = 2*n_basis
  n_basis_soc_ncoll     = 0
  n_basis_soc           = n_basis_soc_coll + n_basis_soc_ncoll
  n_saved_states_soc    = 2*n_states - n_high_states_omit_from_saved_soc - n_core_states_omit_from_saved_soc
  soc_saved_state_start = n_core_states_omit_from_saved_soc + 1

  ! Set up the various BLACS-related variables in scalapack_soc that this
  ! and subsequent SOC-related subroutines will use
  if (use_scalapack) then
    call initialize_scalapack_soc( n_states_soc, n_basis_soc )
  end if

  ! Allocate memory
  write(info_str,'(2X,A)') &
       "Allocating memory for main spin-orbit coupling matrices..."
  call localorb_info(info_str)

  if (use_scalapack) then
    call aims_allocate( soc_ham, mxld_soc, mxcol_soc, "+soc_ham" )
    call aims_allocate(eigenvec_soc_wf_basis, mxld_soc, mxcol_soc, "+eigenvec_soc_wf_basis")
  else

    call aims_allocate( soc_ham, n_states_soc, n_states_soc, "+soc_ham" )
    call aims_allocate(eigenvec_soc_wf_basis, n_states_soc, n_states_soc, "+eigenvec_soc_wf_basis")
  end if

  ! To keep the rank of the eigenvalues array consistent with the SR version, we
  ! include an additional dummy index, where the spin index should be
  ! This allows us to pass the eigenvalues array into subroutines originally
  ! designed for SR quantities
  call aims_allocate(KS_eigenvalue_soc_perturbed, n_saved_states_soc, 1, n_k_points, "+KS_eigenvalue_soc_perturbed")
  call aims_allocate(occ_numbers_soc, n_saved_states_soc, 1, n_k_points, "+occ_numbers_soc")
  call aims_allocate(spin_expectation, n_saved_states_soc, 5, n_k_points, "+spin_expectation")

  call aims_allocate(KS_eigenvalue_soc_perturbed_temp, n_states_soc, "+KS_eigenvalue_soc_perturbed_temp")
  call aims_allocate(sr_to_soc_idxmap, n_states_soc, "+sr_to_soc_idxmap")
  call aims_allocate(spin_expectation_temp, n_states_soc, 5, "+spin_expectation_temp")

  ! Calculate arrays needed for SOC-perturbed eigenfunctions, if they are needed.
  ! Otherwise, allocate dummy arrays.
  if ( save_soc_perturbed_eigenvectors ) then
    ! As SOC induces non-collinear spin polarization, there will be spin channel 
    ! mixing, and the basis elements acquire an explicit spin index.  Each basis 
    ! element may then be written as a tensor product of a computational basis
    ! element and a spinor, doubling the effective size of the basis set.  This
    ! means that the eigenvector has the general form
    ! (2*n_basis,2*n_states,n_k_points_task), where there is no longer an
    ! explicit spin index because spin is no longer a good quantum number of the
    ! system.
    ! Like the eigenvalues array, we include a dummy index in place of the spin
    ! index.
    ! Unlike for the scalar-relativistic eigenvectors (KS_eigenvector and
    ! eigenvec), the SOC-perturbed eigenvectors do not have separate version for
    ! BLACS and the full/global/LAPACK versions.  If use_scalapack is set, then
    ! the vector will be BLACS.  If use_scalapack is not set, then the vector
    ! will be full.
    if (use_scalapack) then
      call aims_allocate( KS_eigenvector_soc_perturbed, mxld_soc_vec, mxcol_soc_vec, 1, 1, "+KS_eigenvector_soc_perturbed")
    else
      call aims_allocate( KS_eigenvector_soc_perturbed, n_basis_soc, n_saved_states_soc, 1, n_k_points_task, "+KS_eigenvector_soc_perturbed")
    end if
  else
    call aims_allocate(KS_eigenvector_soc_perturbed, 1, 1, 1, 1, "+KS_eigenvector_soc_perturbed")
  end if

  ! Set up load balancing
  if (use_local_index.and.use_load_balancing) then
    call aims_allocate(soc_matrix, batch_perm(n_bp_integ)%n_local_matrix_size, 3, "+soc_matrix")
    use_batch_permutation = n_bp_integ
    ld_soc_matrix         = batch_perm(n_bp_integ)%n_local_matrix_size
  else
    call aims_allocate(soc_matrix, n_hamiltonian_matrix_size, 3,                  "+soc_matrix")
    ld_soc_matrix         = n_hamiltonian_matrix_size
  end if

  write (info_str, *)
  call localorb_info( info_str )

  KS_eigenvalue_soc_perturbed = 0.0d0
  occ_numbers_soc = 0.0d0
  spin_expectation = 0.0d0
  saved_state_offset = (soc_saved_state_start-1) - n_core_states_omit_from_soc

  soc_converged = .false.
  soc_enough_walltime_left = .true.

  ! Here is where various quantities related to parallelizing over k-points
  ! are calculated
  ! This should be rolled into a generic function at some point.
  ! Or better yet, be added as a global variable in pbc_list.
  if (n_k_points_task .gt. 0) then
    if (n_k_points .ge. n_tasks) then
      ! Undo the round-robin allocation
      i = 1
      do i_k_point = 1, n_k_points, 1
        if(myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
          my_k_points(i) = i_k_point
          i = i + 1
        end if
      ! NEEDS AN ERROR CONDITION
      end do
    ! The following case can occur when, for some reason, LAPACK is being
    ! used even though there are more processes than k-points.  This can
    ! arise when the number of k-points is less than double the number of
    ! tasks before time-inversion symmetry is applied.  The code determines
    ! that LAPACK should be used, but then reduces the number of k-points
    ! down to less than the number of tasks, so some tasks are empty.  This
    ! can also arise when ScaLAPACK is not compiled against or LAPACK is
    ! manually specified in control.in.
    else if (n_k_points .lt. n_tasks .and. .not.use_scalapack) then
      if (myid .le. n_k_points) then
        my_k_points(1) = myid ! Note that in this case, myid 0 will never have a k-point
      else
        my_k_points(1) = 0 ! This is used later to denote "this MPI rank has no k-points assigned"
      end if
    ! The ScaLAPACK case.
    else
      my_k_points(1) = my_k_point
    end if
  end if

  ! Fill in the matrix for the real-space SOC matrix elements (\Pi from the
  ! Phys Rev Mater. paper)
  call get_timestamps(time_matrix, clock_time_matrix)
  ! Step 3 in Section III.3 from Huhn and Blum, Phys. Rev. Mater. (2017)
  call integrate_soc_matrix (rho,hartree_potential,partition_tab,soc_matrix)
  ! When local indexing is not used, every process has the full real-space
  ! matrices, so we can calculate statistics trivially and output them.
  ! We can do the same when local indexing is used via MPI calls, but I haven't
  ! gotten around to it, to be honest.
  if (.not.use_local_index) then
    ! Find maximum value of SOC matrix
    temp = 0.d0
    do i_coord = 1, 3, 1
      do i_basis_1 = 1, ld_soc_matrix, 1
        if (dabs(soc_matrix(i_basis_1,i_coord)) .gt. temp(1)) then
          temp(1) = dabs(soc_matrix(i_basis_1,i_coord))
          temp(2) = i_basis_1
          temp(3) = i_coord
        end if
      end do
    end do
    ! Determine sparsity of SOC matrix
    ! Absolute threshhold is set by the maximum value times a relative threshhold
    soc_mat_thresh = temp(1)*soc_mat_rel_thresh
    temp_int = 0
    do i_coord = 1, 3, 1
      do i_basis_1 = 1, ld_soc_matrix, 1
        if (dabs(soc_matrix(i_basis_1,i_coord)) .gt. soc_mat_thresh) then
          temp_int = temp_int + 1
        end if
      end do
    end do

    ! Write out some useful information that's debug-sh but might be useful for user
    write(info_str, *)
    call localorb_info( info_str )
    write(info_str,'(2X,A,F17.5)') &
          "Largest matrix element of SOC matrix is            ", temp(1)
    call localorb_info( info_str )
    write(info_str,'(2X,A,E11.5,A,F7.5)') &
         "Sparsity of soc_matrix based on threshhold of ", soc_mat_thresh, " is ", &
         1.0d0 - dble(temp_int)/dble(ld_soc_matrix*3)
    call localorb_info( info_str )
    write(info_str,*)
    call localorb_info( info_str )
  else
    write(info_str, *)
    call localorb_info( info_str )
    write(info_str, '(2X,A)') 'Because use_local_index is enabled, statistics about soc_matrix will not be output.'
    call localorb_info( info_str )
    write(info_str,*)
    call localorb_info( info_str )
  end if

  call get_times(time_matrix, clock_time_matrix, tot_time_matrix, tot_clock_time_matrix, .true.)

  ! The work-horse loop of the code:  perform second-variational SOC at each k-point
  do i_k_point = 1, n_k_points_task
    this_k_point = my_k_points(i_k_point)
    if (this_k_point .eq. 0) then
      exit
    end if

    ! Create the SOC Hamiltonian for the current k-point
    ! Step 4 (non-periodic) and steps 4-5 (periodic) in Section III.3
    !      from Huhn and Blum, Phys. Rev. Mater. (2017)
    call get_timestamps( time_ham, clock_time_ham )
    if (use_scalapack) then
      call construct_SOC_Hamiltonian(ld_soc_matrix, soc_matrix, mxld, mxcol, eigenvec, eigenvec_complex, &
           this_k_point, mxld_soc, mxcol_soc, soc_ham)
    else
      if (real_eigenvectors) then
        call construct_SOC_Hamiltonian(ld_soc_matrix, soc_matrix, n_basis, n_states, KS_eigenvector(1,1,1,i_k_point), &
           KS_eigenvector_complex(1,1,1,1), this_k_point, n_states_soc, n_states_soc, soc_ham)
      else
        call construct_SOC_Hamiltonian(ld_soc_matrix, soc_matrix, n_basis, n_states, KS_eigenvector(1,1,1,1), &
           KS_eigenvector_complex(1,1,1,i_k_point), this_k_point, n_states_soc, n_states_soc, soc_ham)
      end if
    end if
    call get_times(time_ham, clock_time_ham, tot_time_ham, tot_clock_time_ham, .true.)

    ! Output the SOC Hamiltonian to file.  Only for the trivial case of a cluster calculation with LAPACK (i.e. one task)
    ! This has been requested multiple times, so I'm leaving it in here to stop re-implementing it every time.
    ! It should be turned into its own proper subroutine.
    if (.false..and.n_periodic.eq.0 .and. n_tasks.eq.1 .and. myid.eq.0) then 
      open(66,file="SOC_Hamiltonians.out")
      write(66,*) "     i_state     j_state        H_SR [Ha]                               V_SOC [Ha]"
      do i_state = 1, n_states_soc
        do j_state = 1, n_states_soc
          if (i_state.ne.j_state) then
            write(66,*) i_state, j_state, 0.0d0, soc_ham(i_state,j_state)
          else
            if (i_state .le. n_states) then
              write(66,*) i_state, j_state, KS_eigenvalue(i_state,1,1), soc_ham(i_state,j_state)
            else
              if (n_spin.eq.1) then
                write(66,*) i_state, j_state, KS_eigenvalue(i_state-n_states,1,1), soc_ham(i_state,j_state)
              else
                write(66,*) i_state, j_state, KS_eigenvalue(i_state-n_states,1,2), soc_ham(i_state,j_state)
              end if
            end if
          end if
        end do
      end do
      close(66)
    end if

    ! And do the second-variational step
    ! Steps 5-6 (non-periodic) and steps 6-7 (periodic) in Section III.3
    !      from Huhn and Blum, Phys. Rev. Mater. (2017)
    call get_timestamps( time_diag_soc, clock_time_diag_soc )
    if (use_scalapack) then
      call perform_soc_perturbation( mxld_soc, mxcol_soc, soc_ham, KS_eigenvalue(1,1,this_k_point),&
           KS_eigenvalue_soc_perturbed_temp, this_perturb_matrix_max, &
           mxld_soc, mxcol_soc, eigenvec_soc_wf_basis )
    else
      call perform_soc_perturbation( n_states_soc, n_states_soc, soc_ham, KS_eigenvalue(1,1,this_k_point),&
           KS_eigenvalue_soc_perturbed_temp, this_perturb_matrix_max, &
           n_states_soc, n_states_soc, eigenvec_soc_wf_basis )
    end if

    ! Eigensolver will return the eigenvalues for the entire second-variational
    ! window, so we must further restrict to eigenvalues we wish to save
    KS_eigenvalue_soc_perturbed( 1:n_saved_states_soc, 1, this_k_point ) = &
         KS_eigenvalue_soc_perturbed_temp( saved_state_offset+1: saved_state_offset + n_saved_states_soc )
    call get_times(time_diag_soc, clock_time_diag_soc, tot_time_diag_soc, tot_clock_time_diag_soc, .true.)

    ! Calculate the expectation values of spin operators for every eigenstate
    ! 7 September 2017:  I am no longer calculating or outputting these quantities, as I am unclear
    !                    how physical the resulting values are.  That being said, I would like to
    !                    return to this in the future, however, so I am not completely eliminating the code.
    !                    (Note that this code has not been updated to support BLACS.)
!      call calculate_spin_expectation( n_states_sr, eigenvec_soc_wf_basis, spin_expectation_temp )
!      spin_expectation( 1:n_saved_states_soc, :, this_k_point ) = &
!         spin_expectation_temp( saved_state_offset+1: saved_state_offset + n_saved_states_soc, : )

! check cl for orbit
!    if (myid == 0) then
!      do i_basis = 1, n_states_soc
!        write(info_str, '(I3, F10.6, F10.6)') i_basis, real(eigenvec_soc_wf_basis(i_basis, 261)),&
! aimag(eigenvec_soc_wf_basis(i_basis, 261))
!        call localorb_info(info_str)
!      end do
!   end if
! end check

    ! Store the second-variational eigenvectors, but only if we need them
    if ( save_soc_perturbed_eigenvectors ) then
      call get_timestamps( time_eigenvec, clock_time_eigenvec )

      ! The output wavefunctions from the second-variational solution is in terms
      ! of the variational basis set, e.g. the unperturbed wavefunctions.  This
      ! converts them into the more familiar (and usable) form in terms of the
      ! computational basis (NAOs, Gaussians, etc.)
      if (use_scalapack) then
        call convert_wf_basis_to_compute_basis( &
             mxld,         mxcol,         eigenvec, eigenvec_complex,&
             mxld_soc,     mxld_soc,      eigenvec_soc_wf_basis, &
             mxld_soc_vec, mxcol_soc_vec, KS_eigenvector_soc_perturbed )
      else
        if (real_eigenvectors) then
          call convert_wf_basis_to_compute_basis( &
               n_basis,      n_states,           KS_eigenvector(1,1,1,i_k_point), KS_eigenvector_complex(1,1,1,1),&
               n_states_soc, n_states_soc,       eigenvec_soc_wf_basis, &
               n_basis_soc,  n_saved_states_soc, KS_eigenvector_soc_perturbed(1,1,1,i_k_point) )
        else
          call convert_wf_basis_to_compute_basis( &
               n_basis,      n_states,           KS_eigenvector(1,1,1,1), KS_eigenvector_complex(1,1,1,i_k_point),&
               n_states_soc, n_states_soc,       eigenvec_soc_wf_basis, &
               n_basis_soc,  n_saved_states_soc, KS_eigenvector_soc_perturbed(1,1,1,i_k_point) )
        end if
      end if
      call get_times(time_eigenvec, clock_time_eigenvec, tot_time_eigenvec, tot_clock_time_eigenvec, .true.)
    end if
  end do  ! end loop over k-points

  ! Immediaetely deallocate all arrays which are no longer needed
  call aims_deallocate(soc_ham,                                                  "soc_ham" )
  call aims_deallocate(eigenvec_soc_wf_basis,                       "eigenvec_soc_wf_basis")
  call aims_deallocate(soc_matrix,                                             "soc_matrix")

  call get_timestamps( time_sync, clock_time_sync)
  ! Synchronize eigenvalues across all ranks
  if ( use_scalapack .and. my_scalapack_id.ne.0 ) then
    KS_eigenvalue_soc_perturbed = 0.0d0
  end if
  call sync_vector(KS_eigenvalue_soc_perturbed, n_saved_states_soc*n_k_points, mpi_comm_global)

!  7 September 2017:  See previous comment
!  if ( use_scalapack .and. my_scalapack_id.ne.0 ) then
!    spin_expectation = 0.0d0
!  end if
!  call sync_vector(spin_expectation, n_saved_states_soc*5*n_k_points, mpi_comm_global)

  ! Get occupation numbers
  call convert_sr_to_soc_environment ()
  ! Adjust force_occupation constraints to usage in nscf(!)-SOC
  if (.not. force_occupation_projector) then
    call get_occupation_numbers_p0(KS_eigenvalue_soc_perturbed,n_electrons,t_out,occ_numbers_soc,chemical_potential_soc)
  else
    ! GSM: This was only tested for cluster systems since force_occ_periodic is currently broken (Feb 2018)
    do i_k_point = 1, n_k_points
      call create_sorted_sr_eigenvalues_idxmap(n_states, 1, KS_eigenvalue(:,:,i_k_point), sr_to_soc_idxmap)
      call adjust_force_occ_sr(sr_to_soc_idxmap, i_k_point)
    end do
    call get_occupation_numbers_occ_p0(KS_eigenvalue_soc_perturbed,n_electrons,t_out,occ_numbers_soc,chemical_potential_soc)
  end if
  call revert_soc_to_sr_environment ()

  ! Deallocate indexmap, not needed anymore
  call aims_deallocate(sr_to_soc_idxmap, "sr_to_soc_idxmap")

  ! We're now essentially done.  From here, all there is left to do is some minor reshuffling, and then
  ! output the values.
  call get_times(time_sync, clock_time_sync, tot_time_sync, tot_clock_time_sync)
  call get_timestamps(time_post, clock_time_post)

  ! Will output results to main output file as well as stdout
  if (out_soc_eigenvalues) then
    if (myid.eq.0) then
      open(66,file="SOC_eigenvalues.dat")
      write(info_str,'(A,F12.5)') "biggest value in soc-matrix", temp(1)
      call localorb_info(info_str,66,'(A)')
      write(info_str,'(A,F5.0,A,F5.0,A,F5.0)') "at index ", temp(2)," for coordinate ",temp(3)
      call localorb_info(info_str,66,'(A)')
      write(info_str,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      call localorb_info(info_str,66,'(A)')
      write(info_str,*)"Perturbed KS eigenvalues"
      call localorb_info(info_str,66,'(A)')

      ! In retrospect, this seems like overkill.
      do i_k_point = 1, n_k_points
        call write_soc_values( KS_eigenvalue(:,:,i_k_point), KS_eigenvalue_soc_perturbed(1,1,i_k_point), &
             occ_numbers_soc(1,1,i_k_point), spin_expectation(:,:,i_k_point), i_k_point, k_point_list(i_k_point,:), &
             66 )
        write (info_str,*)
        call localorb_info(info_str, 66,'(A)')
      end do
      write(info_str,*)"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      call localorb_info(info_str,66,'(A)')

      close(66)
    end if
  end if

  write (info_str,'(2X,A)') &
    "Writing SOC-perturbed Kohn-Sham eigenvalues."
  call localorb_info(info_str)

  ! Write out SOC values on the s.c.f. k-grid
  do i_k_point = 1, min(n_k_points,out_k_points_eigenvalues)
    call write_soc_values( KS_eigenvalue(:,:,i_k_point), KS_eigenvalue_soc_perturbed(1,1,i_k_point), &
         occ_numbers_soc(1,1,i_k_point), spin_expectation(:,:,i_k_point), i_k_point, k_point_list(i_k_point,:), &
         default_unit )
    write (info_str,*)
    call localorb_info(info_str)
  end do ! i_k_point

  ! In all but the rarest cases, SOC will change the frontier orbitals, so they
  ! and the band gap need to be recalculated
  ! The only exception I can think of is a system where both frontier orbitals
  ! are s-states; however I have yet to see a system where the band gap didn't
  ! change, however slightly
  call convert_sr_to_soc_environment ()
  call find_and_output_homo_lumo_gap(KS_eigenvalue_soc_perturbed, occ_numbers_soc, temp(1), temp(2))
  call revert_soc_to_sr_environment ()

  ! If requested, output the SOC-perturbed work function
  if (n_periodic > 0) then
    if (use_dipole_correction .or.  calculate_work_function) then
      chemical_potential_spin_soc(1) = chemical_potential_soc

      call convert_sr_to_soc_environment ()
      call output_energy_levels(KS_eigenvalue_soc_perturbed, occ_numbers_soc, &
           chemical_potential_soc, chemical_potential_spin_soc, conv)
      call revert_soc_to_sr_environment ()
    end if
  end if

  ! Recalculate the new sum-of-occupied-eigenvalues contribution to total
  ! energy for SOC-perturbed eigenvalues
  ! On the non-self-consistent level of theory, this is the only change to
  ! the total energy possible
  ! We also know from testing that the resulting *total energies* are complete
  ! trash; the conduction/valence eigenvalues however look quite good and can
  ! be improved by include p1/2 orbitals, as WIEN2k does
  ! As eigenvectors can change, at minimum the Hartree term should also be
  ! re-calculated
  ev_sum_new = 0.d0
  call aims_deallocate(KS_eigenvalue_soc_perturbed_temp,             "KS_eigenvalue_soc_perturbed_temp")
  call aims_allocate( KS_eigenvalue_soc_perturbed_temp, 2*n_states, "KS_eigenvalue_soc_perturbed_temp" )
  do i_k_point = 1, n_k_points, 1
    do i_state = 1, n_saved_states_soc, 1
      ev_sum_new = ev_sum_new + k_weights(i_k_point) * occ_numbers_soc(i_state, 1, i_k_point) * &
          KS_eigenvalue_soc_perturbed(i_state, 1, i_k_point)
    end do
    ! If we have chosen not to save low-lying SOC-perturbed states, we now add
    ! in the original scalar-relativistic values for those states, which we
    ! obtain by including both spin channels in KS_eigenvaules and sorting
    ! the resulting array
    if (soc_saved_state_start .gt. 1) then
      call create_sorted_sr_eigenvalues_list( n_states, 1, KS_eigenvalue(:,:,i_k_point), &
           KS_eigenvalue_soc_perturbed_temp )
      do i_state = 1, soc_saved_state_start-1, 1
        ! We require that the occupation numbers for these omitted states are one
        ev_sum_new = ev_sum_new + k_weights(i_k_point) * &
             KS_eigenvalue_soc_perturbed_temp(i_state)
      end do
    end if
  end do
  total_energy_loc = total_energy
  total_energy_loc = total_energy_loc - ev_sum + ev_sum_new
  SOC_non_sc_total_energy = total_energy_loc
  total_energy_corrected_loc = total_energy_loc + entropy_correction

  ! If requested, output the SOC-perturbed DOS
  ! The species/atom projected DOS will be calculated later in the code as part of the Mulliken anlysis
  if(out_dos.and..not.pert_dos_on) then
    call convert_sr_to_soc_environment ()
    call output_KS_dos(KS_eigenvalue_soc_perturbed, occ_numbers_soc, chemical_potential_soc, n_electrons, &
         "KS_DOS_total.dat", "KS_DOS_total_raw.dat")
    call revert_soc_to_sr_environment ()
    write(info_str, *)
    call localorb_info( info_str )
  end if

  ! Write out the SOC-perturbed eigenvectors, if requested
  if (n_write_soc_eigenvectors .gt. 0) then
    write(info_str, '(2X)')
    call localorb_info( info_str )
    write(info_str, '(2X,A)') "Writing SOC-perturbed eigenvectors to disk, as requested."
    call localorb_info( info_str )

    do i = 1, n_write_soc_eigenvectors
      ! The k-point to output SOC eigenvetors at
      this_k_point = write_soc_eigenvectors(i)

      ! Make sure user is sane
      if (this_k_point .lt. 1 .or. this_k_point .gt. n_k_points) then
        write(info_str, '(2X,A,I9,A)') "The requested k-point for outputting SOC-perturbed eigenvectors, ", this_k_point, &
                                       " does not exist.  Skipping."
        call localorb_info( info_str )
        continue
      else
       ! It would be nice to know which task did the writing
        write(info_str, '(2X,A,I9)') "Writing out SOC-perturbed eigenvectors at k-point ", this_k_point
        call localorb_info( info_str )
      end  if

      ! Construct the output filename
      write(k_point_string, '(I9)') this_k_point
      write(info_str, '(A,A)') "SOC_eigenvectors_k_point_", adjustl(k_point_string)

      ! Sort through all k-points assigned to this MPI task, see if the requested 
      ! k-points matches, then output if it does
      if (use_scalapack) then
        if (my_k_point.eq.this_k_point .and. my_scalapack_id.eq.0) then
          call write_SOC_perturbed_eigenvectors( info_str, this_k_point,&
                                                  KS_eigenvector_soc_perturbed(1,1,1,1), &
                                                  KS_eigenvalue_soc_perturbed )
        end if
      else
        do i_k_point = 1, n_k_points_task
          if (my_k_points(i_k_point).eq.this_k_point) then
            call write_SOC_perturbed_eigenvectors( info_str, this_k_point, &
                                                    KS_eigenvector_soc_perturbed(:,:,:,i_k_point), &
                                                    KS_eigenvalue_soc_perturbed(1,1,i_k_point) )
          end if
        end do
      end if
    end do

    write(info_str, '(2X)')
    call localorb_info( info_str )
  end if

  ! out dipole
  if (flag_out_dipmat) then
    write(info_str, *)
    call localorb_info( info_str )
    write(info_str, '(2X,A)') "calculated soc trans dipole mat"
    call localorb_info( info_str )
    if (myid == 0) then
      write(info_str, '(2X,A,F17.5)') 'chemical_potential_soc', chemical_potential_soc
      call localorb_info( info_str )
    end if
    call get_dipolematrix_soc( KS_eigenvalue_soc_perturbed, occ_numbers_soc, &
                               KS_eigenvector_soc_perturbed, &
                               chemical_potential_soc, partition_tab, l_shell_max)
  end if

  ! Now deallocate arrays still lying around
  call aims_deallocate(spin_expectation,                                 "spin_expectation")
  call aims_deallocate(KS_eigenvalue_soc_perturbed_temp, "KS_eigenvalue_soc_perturbed_temp")
  call aims_deallocate(spin_expectation_temp,                       "spin_expectation_temp")

  ! And write out timing information before exiting
  call get_times(time_post, clock_time_post, tot_time_post, tot_clock_time_post)
  call get_times(time_soc, clock_time_soc, tot_time_soc, tot_clock_time_soc)
  write(info_str, *)
  call localorb_info( info_str )
  write(info_str, *)  " Spin-orbit coupling                                     :  max(cpu_time)    wall_clock(cpu1)"
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Computing matrix elements for SOC operator            :", &
       tot_time_matrix, tot_clock_time_matrix
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Constructing second-variational Hamiltonian           :", &
      tot_time_ham, tot_clock_time_ham
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Diagonalizing second-variational Hamiltonian          :", &
      tot_time_diag_soc, tot_clock_time_diag_soc
  call localorb_info( info_str )
  if ( save_soc_perturbed_eigenvectors ) then
    write(info_str, "(2X, A,F15.3,F20.3)")  "| Constructing SOC-perturbed eigenvectors               :", &
        tot_time_eigenvec, tot_clock_time_eigenvec
    call localorb_info( info_str )
  end if
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Syncing data and calculating occupations              :", &
      tot_time_sync, tot_clock_time_sync
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Writing out data, calculating band gaps, etc.         :", &
      tot_time_post, tot_clock_time_post
  call localorb_info( info_str )
  write(info_str, "(2X, A,F15.3,F20.3)")  "| Total SOC time                                        :", &
      tot_time_soc, tot_clock_time_soc
  call localorb_info( info_str )
  write(info_str, *)
  call localorb_info( info_str )
  call aims_mem_current_output ()
  write(info_str, *) "------------------------------------------------------------"
  call localorb_info( info_str )

  ! FIN

end subroutine calculate_second_variational_soc
!******
