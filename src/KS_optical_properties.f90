!***h* FHI-aims/KS_optical_properties
!  NAME
!    KS_optical_properties
!  SYNOPSIS
module KS_optical_properties
  !  PURPOSE
  !    This module, in conjunction with calculate_mommat_base, contains
  !    subroutines used to computes optical properties of a material via the
  !    Lindhard formalism (i.e. RPA approximation using single-particle (g)KS
  !    orbitals. This module should not be confused with the optical_response
  !    module, which uses TDDFT.
  !
  !    To understand the formalism used, please refer to Section 2 and Appendix
  !    A of Ambrosch-Draxl and Sofo, Comput. Phys. Commun. 175 (2006), 1-14.
  !    The Ambrosch-Draxl paper refers to the implementation of the Lindhard
  !    formalism in WIEN2k, and thus the messy implementation details are
  !    different.  However, as both codes are all-electron, the overarching
  !    design choices made are similar.
  !
  !    We calculate the components of the macroscopic dielectric tensor for a
  !    material by first calculating the interband and intraband contributions
  !    to the polarizability, and then obtaining the dielectric tensor via
  !    dielectric = vacuum + polarizability = \delta_ij + polarizablity.
  !    We have chosen to make this (seemingly) petty distinction because
  !    multiple literature references contain equations which double-count the
  !    vacuum contribution/Kronecker delta by adding it to the interband
  !    and intraband contribution to the dielectric constant in different
  !    contexts.  By working with the polarizability instead, this source of
  !    error is avoided.  [There is still the original physically arbitrary
  !    distinction of interband versus intraband, but this is a useful
  !    computational distinction.]
  !  USES
  use scalapack_wrapper, only: dlen_
  implicit none
  !  AUTHOR
  !    William Huhn and Tong Zhu (Duke University), based on earlier code by
  !    Bjoern Bieniek (formerly FHI-Berlin)
  !  HISTORY
  !    May 2017 - Created
  !    September 2017 - Updated to support ScaLAPACK
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement."
  !  TODO
  !    o Move mommat_full_{x,y,z}_* directly into subroutines instead of module
  !      variables
  !    o Change name of {x,y,z}_on to something intelligible out of context
  !    o Eliminate relics of old dielectric code throughout code
  !  SOURCE

  ! Temporary variables (should be moved into subroutines rather than here)
  real*8, allocatable :: mommat_full_x_up(:)
  real*8, allocatable :: mommat_full_x_low(:)
  real*8, allocatable :: mommat_full_y_up(:)
  real*8, allocatable :: mommat_full_y_low(:)
  real*8, allocatable :: mommat_full_z_up(:)
  real*8, allocatable :: mommat_full_z_low(:)

  ! ScaLAPACK-related variables for indexing momentum matrix
  integer :: block_size_moment
  integer :: sc_desc_moment(dlen_)
  integer :: mb_moment, nb_moment
  integer :: mxld_moment, mxcol_moment
  integer :: n_my_rows_moment, n_my_cols_moment
  integer, allocatable :: l_row_moment(:), l_col_moment(:)
  integer, allocatable :: my_row_moment(:), my_col_moment(:)

  ! Module variables
  real*8,     allocatable :: omega_plasma_sq(:)
  complex*16, allocatable :: interband_polarizability(:,:)
  complex*16, allocatable :: intraband_polarizability(:,:)
  complex*16, allocatable :: interband_polarizability_reset(:,:)
  complex*16, allocatable :: dielectric(:,:)
  real*8,     allocatable :: absorption_coeff(:,:)
  complex*16, allocatable :: interband_polarizability_imagfreq(:,:)

  ! Timings (should be placed in different module, to be honest)
  real*8 :: tot_time_real_space_mommat = 0.0d0
  real*8 :: tot_time_bloch_mommat = 0.0d0
  real*8 :: tot_time_plasma_interband = 0.0d0
  real*8 :: tot_time_sync_postprocess = 0.0d0
  real*8 :: tot_clock_time_real_space_mommat = 0.0d0
  real*8 :: tot_clock_time_bloch_mommat = 0.0d0
  real*8 :: tot_clock_time_plasma_interband = 0.0d0
  real*8 :: tot_clock_time_sync_postprocess = 0.0d0

  ! Drivers
  public :: KS_dielectric_calculation
  public :: KS_dielectric_calculation_imagfreq
  ! Utility functions
  public :: calculate_moment_matrix
  public :: calc_moment_w0
  public :: calc_moment_w0_soc
  public :: calculate_plasma_frequency_sq
  public :: calculate_plasma_frequency_sq_imagfreq
  public :: calculate_interband_polarizability
  public :: calculate_interband_polarizability_imagfreq
  public :: kramerskronig_im_to_re
  public :: calculate_drude_like_intraband_polarizability
  public :: calculate_dielectric_from_polarizability
  public :: calculate_absorption_coeff
  public :: output_optical_properties
  public :: simpson_int
  ! Memory management functions
  public :: allocate_mommat_TZ
  public :: init_KS_optical_properties
  public :: cleanup_KS_optical_properties

contains

  !!!!!!!!!!!!!!!
  ! MAIN DRIVER !
  !!!!!!!!!!!!!!!

  subroutine KS_dielectric_calculation &
       ( n_rows, n_cols, n_spin_in, KS_vec, KS_vec_complex, &
         n_states_in, KS_eigen, occ_numbers, chemical_potential, &
         spin_degeneracy_in, partition_tab, l_shell_max )
    ! PURPOSE
    !
    !   Input:
    !   n_plot_dielectric
    !   broaden_type(n_plot_dielectric), broaden_para(n_plot_dielectric)
    !
    !    Wrapper function for calculating and outputting of the marcroscopic
    !    linear dielectric function in three steps
    !    0. ep1_in, ep2_in from character to integer (x->1,y->2,z->3);
    !       allocations
    !    1. function 'calculate_mommat_p0', real space integral of
    !       atom-centered basis function i with the gradient of function j for
    !       in each cell (cell distance)
    !    2. function 'construct_overlap' 'fourier transform' of real space
    !       integrals by summing up the phases resulting from cell distances
    !       function is called once for the upper triangle and once for the
    !       lower triangle of the sparse matrix from first step
    !    3. function 'calc_moment_w0' basis transformation from atom-centered
    !       basis to KS-orbitals
    !    4. function 'calc_dielectric_function' momentum matrix elements are
    !       summed up to yield the dielectric function
    !    5. function 'output_optical_properties' output of result to file
    ! USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use calculate_mommat_base, only: get_state_minmax_K, calculate_mommat_p0
    use constants, only: hartree
    use dimensions, only: n_full_points, n_species, n_k_points, &
        n_k_points_task, calculate_perturbative_soc, n_plot_dielectric
    use dimensions_soc, only: n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll
    use localorb_io, only: localorb_info
    use mpi_tasks, only: myid, n_tasks, aims_stop_coll, mpi_comm_global
    use pbc_lists, only: k_weights
    use load_balancing, only: use_batch_permutation, n_bp_integ, batch_perm
    use physics, only: n_electrons
    use runtime_choices, only: n_omega, n_omega_reset, omega_min,  omega_max, &
        Emin, Emax, use_local_index, calc_px_dielectric, calc_py_dielectric, &
        calc_pz_dielectric, direc_1, direc_2, broaden_para, broaden_type, &
        real_eigenvectors, use_scalapack, use_load_balancing, &
        dipole_trans_mat_k, dipole_trans_mat_out, dipole_trans_mat_orbit_num
    use scalapack_wrapper, only: my_k_point, my_scalapack_id, &
        my_scalapack_comm_all
    use soc_utilities, only: convert_sr_to_soc_environment, &
        revert_soc_to_sr_environment
    use synchronize_mpi_basic, only: sync_vector_complex, sync_real_number
    use timing, only: get_times, get_timestamps, tot_time_dielectric, &
        tot_clock_time_dielectric
    implicit none
    !  ARGUMENTS
    integer,    intent(in) :: n_rows
    integer,    intent(in) :: n_cols
    integer,    intent(in) :: n_spin_in
    real*8,     intent(in) :: KS_vec(n_rows,n_cols,n_spin_in,n_k_points_task)
    complex*16, intent(in) :: &
         KS_vec_complex(n_rows,n_cols,n_spin_in,n_k_points_task)
    integer,    intent(in) :: n_states_in
    real*8,     intent(in) :: KS_eigen(n_states_in,n_spin_in,n_k_points)
    real*8,     intent(in) :: occ_numbers(n_states_in,n_spin_in,n_k_points)
    real*8,     intent(in) :: chemical_potential
    real*8,     intent(in) :: spin_degeneracy_in
    real*8,     intent(in), target :: partition_tab(n_full_points)
    integer,    intent(in) :: l_shell_max(n_species)
    !  INPUTS
    !    o n_rows - Numbers of rows for the eigenvectors.  Should be the number
    !      of basis functions for a LAPACK-style eigenvector and the appropriate
    !      mxld-like variable for a ScaLAPACK eigenvector.
    !    o n_cols - Numbers of columns for the eigenvector.  Should be the
    !      number of states for a LAPACK-style eigenvector and the appropriate
    !      mxcol-like variable for a ScaLAPACK eigenvector.
    !    o n_spin_in - Numbers of spin channels in the eigenvectors and
    !      eigenvalues.  Should be n_spin for scalar-relativistic calculations
    !      and 1 for spin-orbit-coupled calculations.
    !    o KS_vec - Real eigenvector.  for spin-orbit-coupled calculations, a
    !      dummy variable.
    !    o KS_vec_complex - Complex eigenvector; for spin-orbit-coupled
    !      calculations, *the* eigenvector.
    !    o n_states_in - Numbers of states in eigenvalues and occupation numbers
    !    o KS_eigen - Eigenvalues
    !    o occ_numbers - Occupation numbers.
    !    o chemical_potential - Fermi level.
    !    o spin_degeneracy - The degeneracy factor due to spin channels.  For
    !      scalar-relativistic calculations, this will be the spin_degeneracy
    !      global variable (itself 1.0d0 for spin-polarized calculations and
    !      2.0d0 for spin-non-polarized calculations).  For spin-orbit-coupled
    !      calculations, this is always 1.0d0.
    !    o partition_tab - Integration weights.
    !    o l_shell_max - Maximum angular momentum channels.
    !  OUTPUT
    !    o None (writes to screen and disc)
    !  AUTHOR
    !    William Huhn and Tong Zhu (Duke University), based on earlier code by
    !    Bjoern Bieniek (formerly FHI-Berlin)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    !  local variables
    character*128 :: info_str
    integer :: info
    integer :: n_state_min_in
    integer :: n_state_max_in
    real*8  :: homo_level, lumo_level
    logical :: component_is_diagonal, real_eigenvectors_save
    integer :: n_states_window
    integer :: ld_mommat

    !  counters
    integer :: i, i_k, i_k_point
    integer :: new_k_point
    integer :: coord_1, coord_2
    integer :: my_k_points(n_k_points)
    !TZ changes for output many direction dielectric function
    integer :: i_dielectric, i_omega
    ! Timings
    real*8 :: time_dielectric = 0.0d0, clock_time_dielectric = 0.0d0
    real*8 :: time_real_space_mommat = 0.0d0, clock_time_real_space_mommat = 0.0d0
    real*8 :: time_bloch_mommat = 0.0d0, clock_time_bloch_mommat = 0.0d0
    real*8 :: time_plasma_interband = 0.0d0, clock_time_plasma_interband = 0.0d0
    real*8 :: time_sync_postprocess = 0.0d0, clock_time_sync_postprocess = 0.0d0

    character*50 :: file_name_x, file_name_y, file_name_z
    integer :: current_row_num, current_row_start, current_row_end
    integer :: row_id, col_id, moment_id, n, state_low, state_high, homo_id
    integer :: row_id_real, col_id_real
    real*8 :: trans_energy

    complex*16, allocatable :: moment_TZ(:,:)
    real*8   :: omega_step

    character(*), parameter :: func = 'KS_dielectric_calculation'

    call get_timestamps( time_dielectric, clock_time_dielectric )

    ! WPH:  Local indexing and load balancing currently do not work with
    ! dielectric calculations.  The problem is the packing for the various
    ! mommat arrays, which are misnamed.  They're not real-space momentum
    ! matrices as the names imply, they're actually the real-space matrix
    ! elements for the gradient operator, which is anti-symmetric.   The
    ! handling of matrices in aims is largely set up for Hermitian matrices, so
    ! the data structures must be massaged to support anti-symmetric matrices
    ! without complete rewrites. This is why the choice was made to store both
    ! the lower and upper matrix elements, whereas usually one only stores the
    ! upper portion, and why update_full_matrix_p0X was converted into the
    ! virtually-identical update_full_matrix_mommat, where only the indexing of
    ! the lower half of the shell differs.  Whether this was a good choice, I
    ! can't say.  Seems to be me that only one half needs to be stored and the
    ! other half can be filled in as the negative, but there's probably some
    ! cancellation of terms that I'm overlooking.
    !
    ! However, local indexing and load balancing were written assuming the
    ! real-space matrix is symmetric, so we can't drop-in subroutines as easily
    ! as is done for the standard packed matrix format.  The errors for my test
    ! cases for local indexing without load balancing were minor, around a 1%
    ! difference in the interband contribution to the polarizability.  Why it's
    ! so small, I can't say, (possibly my cases are small so that each processor
    ! has most of the packed real-space Hamiltonian already) but the numbers
    ! should be exactly the same as the LAPACK cases, so this can't be blamed on
    ! processor jitter.  As for local indexing with load balancing, the numbers
    ! are completely wrong, by factors of two or more.
    !
    ! The parts that I believe need to be modified are the
    ! update_full_matrix_mommat() call in calculate_mommat_p0() (as well as the
    ! load balancing part immediately before the call) and the first half of
    ! calculate_moment_matrix() where the real-space matrix is converted into a
    ! BLACS format.

    ! Should have been caught already in read_control, but just in case
    if (use_local_index) then
      call aims_stop_coll ('* Dielectric calculations currently do not support &
                           &use_local_index keyword.  Exiting', func)
    end if
    if (use_load_balancing) then
      call aims_stop_coll ('* Dielectric calculations currently do not support &
                           &load_balancing keyword.  Exiting', func)
    end if

    if (calculate_perturbative_soc) then
      ! Force the code to traverse the complex path
      ! in all cases, since SOC-perturbed eigenvectors
      ! are complex
      real_eigenvectors_save = real_eigenvectors
      real_eigenvectors      = .false.

      ! Some light error checking
      if (mod(n_basis_soc_coll,2) .eq. 1) then
        call aims_stop_coll('The number of collinear basis functions in SOC is &
                           &odd, exiting.', func)
      end if
      if (n_basis_soc_ncoll .gt. 0) then
        call aims_stop_coll('This subroutine does not support usage of &
                            &non-collinear basis functions in SOC, exiting.', &
                            func)
      end if
    end if

    ! Here is where various quantities related to parallelizing over k-points
    ! are calculated
    ! This should be rolled into a generic function at some point.
    ! Or better yet, be added as a global variable in pbc_list.
    if (n_k_points_task .gt. 0) then
      if (n_k_points .ge. n_tasks) then
        ! Undo the round-robin allocation
        i = 1
        do i_k_point = 1, n_k_points, 1
          if (myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
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
          ! Note that in this case, myid 0 will never have a k-point
          my_k_points(1) = myid
        else
          ! This is used later to denote this MPI rank has no k-points assigned
          my_k_points(1) = 0
        end if
      ! The ScaLAPACK case.
      else
        my_k_points(1) = my_k_point
      end if
    end if

    ! -------------------------------------------------------------------------!
    ! FIRST STEP:
    !              allocate needed memory
    ! -------------------------------------------------------------------------!

    ! Determine the state indices corresponding to the user-provided energy
    ! window
    if (calculate_perturbative_soc) then
      call convert_sr_to_soc_environment ()
    end if
    write (info_str, '(2X,A)') "Recalculating HOMO and LUMO levels to &
                               &determine eigenstates needed for dielectric &
                               &calculations..."
    call localorb_info( info_str )
    write (info_str, *) ""
    call localorb_info( info_str )
    call find_and_output_homo_lumo_gap &
                ( KS_eigen, occ_numbers, homo_level, lumo_level )

    ! The energy window are defined acoording to the frequency window you want
    ! to output
    !  Emin  = HOMO - Omega_max - safety
    !  Emax = LUMO + omega_max + safety     the safety value are setting to 10.0
    !  eV
    Emin = homo_level*hartree - omega_max - 10.0
    Emax = lumo_level*hartree + omega_max + 10.0

    call get_state_minmax_K( KS_eigen, n_state_min_in, n_state_max_in )
    if (calculate_perturbative_soc) then
      call revert_soc_to_sr_environment ()
    end if
    n_states_window = n_state_max_in - n_state_min_in + 1

    write (info_str, *) ""
    call localorb_info( info_str )
    write (info_str, '(2X,A,I9)') "Index of first state to include in &
                                  &dielectric calculations : ", n_state_min_in
    call localorb_info( info_str )
    write (info_str, '(2X,A,I9)') "Index of last state to include in &
                                  &dielectric calculations  : ", n_state_max_in
    call localorb_info( info_str )
    write (info_str, *)
    call localorb_info( info_str )

    ! Initialize the ScaLAPACK environment for the moment matrix (which has a
    ! leading dimension of n_states_window)
    if (use_scalapack) then
      write(info_str, '(2X,A)') "Creating ScaLAPACK distribution for matrices&
                                 & with leading dimension n_states_window"
      call localorb_info( info_str )
      call aims_allocate( l_row_moment, n_states_window,   "l_row_moment" )
      call aims_allocate( l_col_moment, n_states_window,   "l_col_moment" )
      call aims_allocate( my_row_moment, n_states_window, "my_row_moment" )
      call aims_allocate( my_col_moment, n_states_window, "my_col_moment" )
      block_size_moment = 0
      call initialize_scalapack_descriptor_matrix( &
           n_states_window, n_states_window, block_size_moment, sc_desc_moment,&
           mb_moment, nb_moment, mxld_moment, mxcol_moment, n_my_rows_moment, &
           n_my_cols_moment, l_row_moment, l_col_moment, my_row_moment, &
           my_col_moment )
      write(info_str, *)
      call localorb_info( info_str )
    else
      call aims_allocate( l_row_moment, 1,   "l_row_moment" )
      call aims_allocate( l_col_moment, 1,   "l_col_moment" )
      call aims_allocate( my_row_moment, 1, "my_row_moment" )
      call aims_allocate( my_col_moment, 1, "my_col_moment" )
    end if

    ! Initialize load balancing, if activated
    if (use_local_index.and.use_load_balancing) then
      use_batch_permutation = n_bp_integ
    end if

    ! Now allocate memory
    write(info_str,'(2X,A)') &
         "Allocating memory for main dielectric-related matrices..."
    call localorb_info(info_str)

    call init_KS_optical_properties()
    call allocate_mommat_TZ( ld_mommat )

    ! This part is request for gaussian broadening, in the gaussian
    ! broadening, as we use the Kramers-Kronig transformation to obtain the real
    ! part, which needs a converged omega values to obtain a converged values.
    ! As a result, we let the code to calculated the imaginary part in the
    ! omega range of all the considered "possible" transitions (omega_max =
    ! Emax-Emin). But still output the results in the User requested range.
    omega_step = (omega_max - omega_min)/n_omega
    n_omega_reset = n_omega + INT((Emax-Emin - omega_max)/omega_step)

    call aims_allocate(interband_polarizability_reset, &
         n_omega_reset,n_plot_dielectric, "interband_polarizability_reset" )

    ! WPH:  moment_TZ are stored in a lower triangular format, which will cause
    ! a memory bottleneck for large calculations.  They should be stored in
    ! BLACS, and indeed, when running with use_scalapack they are originally
    ! calculated in a BLACS format then converted back to low triangular and
    ! synced across processes.
    !
    ! This is done because the calculation of the plasma frequency and interband
    ! contribution are written assuming the low triangular format.  However, the
    ! calculations would be trivial to rewrite to support BLACS; both are loops
    ! over matrix elements with no communication between processes needed.  The
    ! same loop could be done over the local matrix (with suitable indexing
    ! arrays between local and global indices of course) then sync'd across
    ! processes.
    !
    ! I haven't implemented this, since my gut instinct is that this memory
    ! bottleneck will hit after the lack of use_local_index/load_balancing
    ! support.  While moment_TZ is a gloval matrix, it is still a global matrix
    ! in a restricted energy window, which is a lower dimensionality than the
    ! basis set.  Still, this should be done.
    if (calculate_perturbative_soc) then
      if (.not.allocated( moment_TZ )) then
        call aims_allocate(moment_TZ, &
             (n_states_window*(n_states_window+1))/2, 3, "+moment_TZ")
      end if
    else
      if (.not.allocated( moment_TZ )) then
        call aims_allocate(moment_TZ, &
             (n_states_window*(n_states_window+1))/2*n_spin_in, 3, &
             "+moment_TZ")
      end if
    end if

    moment_TZ = (0.0d0, 0.0d0)
    ! -------------------------------------------------------------------------!
    ! SECOND STEP:
    !              calculate real-space momentum matrix elements needed based on
    !              tensor components requested by the user (i.e. if they request
    !              the xy component, calculate the p_x and p_y matrix elements
    !              but not p_z)
    ! -------------------------------------------------------------------------!
    call get_timestamps( time_real_space_mommat, clock_time_real_space_mommat )
    write(info_str, *)
    call localorb_info ( info_str )
    if (calc_px_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_x matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0 ( partition_tab, l_shell_max, &
                                 mommat_full_x_up, mommat_full_x_low,1)
      write(info_str, *)
      call localorb_info ( info_str )
    end if

    if (calc_py_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_y matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0 ( partition_tab, l_shell_max, &
                                 mommat_full_y_up, mommat_full_y_low,2)
      write(info_str, *)
      call localorb_info ( info_str )
    end if

    if (calc_pz_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_z matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0 ( partition_tab, l_shell_max, &
                                 mommat_full_z_up, mommat_full_z_low,3)
      write(info_str, *)
      call localorb_info ( info_str )
    end if
    call get_times(time_real_space_mommat, clock_time_real_space_mommat, &
         tot_time_real_space_mommat, tot_clock_time_real_space_mommat, .true.)

    ! -------------------------------------------------------------------------!
    ! THIRD STEP:
    !             Calculate quantities requiring summation of k-points.
    !             At this time, these would be the plasma frequency (squared)
    !             and the interband contribution to the polarizability
    ! -------------------------------------------------------------------------!

    ! WPH:  This part of the code currently uses too much memory.  For any given
    ! component of the dielectric function, we need at most two momentum
    ! matrices, however, here we compute and store all three.  This is
    ! particularly severe if only the diagonal elements of the dielectric tensor
    ! are needed, which is a very common case, either because the system is
    ! orthorhombic (or better) or because we only need the trace of the
    ! dielectric tensor.
    !
    ! Pseudocode for a better way to do this loop would be:
    !
    ! do i_k = 1,n_k_points_tasks
    !   calculate p_z matrix elements
    !   compute all zz tensor components
    !
    !   calculate p_y matrix elements
    !   compute all zy tensor components
    !   compute all yy tensor components
    !
    !   throw away p_y matrix elements
    !   calculate p_x matrix elements
    !   compute all zx tensor components
    !   compute all xx tensor components
    !
    !   throw away p_z matrix elements
    !   re-calculate p_y matrix elements
    !   compute all xy tensor components
    ! end do
    !
    ! where it is omitted that we only calculate p_i matrix elements when the
    ! associated tensor components are needed, and we throw away matrix elements
    ! earlier if they're no longer needed
    !
    ! The ordering is specifically chosen:  when using the standard unit cell
    ! conventions, the xy tensor component often turns out to be the component
    ! most likely to be zero or equivalent to xx or yy, sidestepping the
    ! recalculation of p_y.
    !
    ! In the best case (where only diagonal components are needed), this
    ! version has a third of the memory cost with no difference in performance.
    ! In the worst case (where all components are needed), this loop reduces the
    ! memory cost by a third, but increases the runtime by a third, since the
    ! p_y matrix elements are calculated twice.  I think that this is an
    ! acceptable compromise, but I would like to see how the current code
    ! performs before making these changes.

    write(info_str,'(2X,A)') "Calculating dielectric tensor elements..."
    call localorb_info ( info_str )

    omega_plasma_sq = 0.0d0
    interband_polarizability  = (0.0d0,0.0d0)
    interband_polarizability_reset = (0.0d0,0.0d0)
    do i_k = 1,n_k_points_task
      new_k_point = my_k_points(i_k)

      ! TZ:  calculate possible moment matrix for needed direction
      !------------------------------------------------------------------------!
      ! THIRD STEP, PART 1:
      !            Calculate the momentum matrix elements at this k-point (e.g.
      !            in the Bloch basis representation)
      !------------------------------------------------------------------------!
      ! WPH:  these moment matrices are O(n_state^2) and will lead to memory
      ! issues down the road.  (See early comment when allocating.)
      call get_timestamps( time_bloch_mommat, clock_time_bloch_mommat )
      if (calc_px_dielectric) then
        call calculate_moment_matrix(ld_mommat, mommat_full_x_up, &
             mommat_full_x_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
             KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
             n_states_window, moment_TZ(:,1))
      end if

      if (calc_py_dielectric) then
        call calculate_moment_matrix(ld_mommat, mommat_full_y_up, &
             mommat_full_y_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
             KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
             n_states_window, moment_TZ(:,2))
      end if

      if (calc_pz_dielectric) then
        call calculate_moment_matrix(ld_mommat, mommat_full_z_up,&
             mommat_full_z_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
             KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
             n_states_window, moment_TZ(:,3) )
      end if

      ! LC: output dipole transition matrix
      !------------------------------------------------------------------------!
      if (dipole_trans_mat_out .and. new_k_point == dipole_trans_mat_k) then
        if (calculate_perturbative_soc) then
          homo_id = int(n_electrons)
        else
          homo_id = int(n_electrons) / 2
        end if
        state_low = homo_id - dipole_trans_mat_orbit_num + 1
        state_high = homo_id + dipole_trans_mat_orbit_num
        ! check
        write(file_name_x, '(A,I1,A4)') 'dipole_transition_matrix_x_', &
             dipole_trans_mat_k,'.out'
        write(file_name_y, '(A,I1,A4)') 'dipole_transition_matrix_y_', &
             dipole_trans_mat_k,'.out'
        write(file_name_z, '(A,I1,A4)') 'dipole_transition_matrix_z_', &
             dipole_trans_mat_k,'.out'
        open(unit = 55, file = file_name_x)
        open(unit = 56, file = file_name_y)
        open(unit = 57, file = file_name_z)
        write(55, '(2X,A10,2X,A10,2X,A10,2X,A10,2X,A7,2X,A7,2X,A10,2X,A10)') &
             'E_trans', 'Re', 'Im', 'ID 1', 'ID 2', 'E1', 'E2'
        write(56, '(2X,A10,2X,A10,2X,A10,2X,A10,2X,A7,2X,A7,2X,A10,2X,A10)') &
             'E_trans', 'Re', 'Im', 'ID 1', 'ID 2', 'E1', 'E2'
        write(57, '(2X,A10,2X,A10,2X,A10,2X,A10,2X,A7,2X,A7,2X,A10,2X,A10)') &
             'E_trans', 'Re', 'Im', 'ID 1', 'ID 2', 'E1', 'E2'

        n = n_states_window

        do col_id = 1, n ! col_id denotes homo orbitals
          col_id_real = col_id + n_state_min_in
          if (col_id_real >= state_low .and. col_id_real <= homo_id) then
            current_row_num = n - col_id + 1
            current_row_end = (n + current_row_num)*(n - current_row_num + 1)/2
            current_row_start = current_row_end - current_row_num + 1

            do row_id = 1, current_row_num ! row_id denotes lumo orbitals
              row_id_real = row_id + n_state_min_in
              if (row_id_real > homo_id .and. row_id_real <= state_high) then
                moment_id = current_row_start + row_id - 1
                trans_energy = KS_eigen(row_id_real, 1, dipole_trans_mat_k) - &
                     KS_eigen(col_id_real, 1, dipole_trans_mat_k)

                write(55, '(2X,F10.5,2X,F10.5,2X,F10.5,2X,I7,2X,I7,2X,F10.5,2X,F10.5)') &
                     trans_energy * hartree, real(moment_TZ(moment_id, 1)), &
                     aimag(moment_TZ(moment_id, 1)), col_id_real, row_id_real, &
                     KS_eigen(col_id_real, 1, dipole_trans_mat_k) * hartree, &
                     KS_eigen(row_id_real, 1, dipole_trans_mat_k) * hartree

                write(56, '(2X,F10.5,2X,F10.5,2X,F10.5,2X,I7,2X,I7,2X,F10.5,2X,F10.5)') &
                     trans_energy * hartree, real(moment_TZ(moment_id, 2)), &
                     aimag(moment_TZ(moment_id, 1)), col_id_real, row_id_real, &
                     KS_eigen(col_id_real, 1, dipole_trans_mat_k) * hartree, &
                     KS_eigen(row_id_real, 1, dipole_trans_mat_k) * hartree

                write(57, '(2X,F10.5,2X,F10.5,2X,F10.5,2X,I7,2X,I7,2X,F10.5,2X,F10.5)') &
                     trans_energy * hartree, real(moment_TZ(moment_id, 3)), &
                     aimag(moment_TZ(moment_id, 1)), col_id_real, row_id_real, &
                     KS_eigen(col_id_real, 1, dipole_trans_mat_k) * hartree, &
                     KS_eigen(row_id_real, 1, dipole_trans_mat_k) * hartree
              end if
            end do ! row_id
          end if
        end do ! col_id
        close(55)
        close(56)
        close(57)
      end if !new_k_point == dipole_trans_mat_k

      call get_times(time_bloch_mommat, clock_time_bloch_mommat, &
           tot_time_bloch_mommat, tot_clock_time_bloch_mommat, .true.)

      ! -----------------------------------------------------------------------!
      ! THIRD STEP, PART 2:
      !             Calculate the plasma frequency (squared) and the interband
      !             contribution to the polarizability for all requested
      !             tensor components and broadening functions
      ! -----------------------------------------------------------------------!
      call get_timestamps( time_plasma_interband, clock_time_plasma_interband )
      do i_dielectric = 1, n_plot_dielectric
        ! WPH:  The following subroutine call re-uses the broadening parameter
        ! from the interband broadening (which may be Lorentzian) as the
        ! Gaussian broadening.  This probably shouldn't be done, but I can't
        ! quite say it's unphysical.
        call calculate_plasma_frequency_sq( n_state_min_in, n_states_window, &
             moment_TZ(:,direc_1(i_dielectric)), &
             moment_TZ(:,direc_2(i_dielectric)), n_states_in, n_spin_in, &
             KS_eigen(:,:,new_k_point), chemical_potential, spin_degeneracy_in,&
             k_weights(new_k_point), omega_plasma_sq(i_dielectric) )

        ! WPH:  The following subroutine recalculates the occupation numbers
        ! from a Fermi-Dirac distribution.  The occupation numbers should be
        ! passed in and used directly, instead.
        call calculate_interband_polarizability(n_state_min_in, &
             n_states_window, moment_TZ(:,direc_1(i_dielectric)), &
             moment_TZ(:,direc_2(i_dielectric)), n_states_in, n_spin_in, &
             KS_eigen(:,:,new_k_point), occ_numbers(:,:,new_k_point), &
             chemical_potential, spin_degeneracy_in,&
             k_weights(new_k_point), broaden_type(i_dielectric), &
             broaden_para(i_dielectric), &
             interband_polarizability(:,i_dielectric), &
             interband_polarizability_reset(:,i_dielectric))
      end do
      call get_times(time_plasma_interband, clock_time_plasma_interband, &
           tot_time_plasma_interband, tot_clock_time_plasma_interband, .true. )
    end do

    !--------------------------------------------------------------------------!
    ! THIRD STEP, PART 3:
    !             Each MPI task still has only a partial k-point sum, now sum
    !             over all MPI tasks to get the full k-point sum
    !--------------------------------------------------------------------------!
    call get_timestamps( time_sync_postprocess, clock_time_sync_postprocess )
    if (use_scalapack .and. my_scalapack_id.ne.0) then
      interband_polarizability = (0.0d0,0.0d0)
      omega_plasma_sq = 0.0d0
      interband_polarizability_reset = (0.0d0,0.0d0)
    end if
    do i_dielectric = 1, n_plot_dielectric
      call sync_vector_complex(interband_polarizability(:,i_dielectric), &
           n_omega)
      call sync_real_number (omega_plasma_sq(i_dielectric) )
      call sync_vector_complex (interband_polarizability_reset(:,i_dielectric),&
           n_omega_reset)
    end do

    !--------------------------------------------------------------------------!
    ! FOURTH STEP:
    !             Calculate all remaining quantities (intraband contribution to
    !             polarizability, dielectric tensor, absorption coefficients)
    !--------------------------------------------------------------------------!
    do i_dielectric = 1, n_plot_dielectric
      component_is_diagonal = direc_1(i_dielectric).eq.direc_2(i_dielectric)

      ! When using Gaussian broadening, the Kramers-Kronig transformation is
      ! used to obtain the real part of the interband contribution from the
      ! imaginary part
      if (broaden_type(i_dielectric) == 1) then
        write(info_str,'(2X,A,1X,I4)') &
             "Performing Kramers-Kronig transform for interband contribution &
             &to dielectric component # ", i_dielectric
        call localorb_info (info_str)
        call kramerskronig_im_to_re( &
             interband_polarizability_reset(:,i_dielectric))
        do i_omega = 1, n_omega
          interband_polarizability(i_omega,i_dielectric) = &
               interband_polarizability_reset(i_omega,i_dielectric)
        end do
      end if

      ! At this point, we have the interband contribution to the polarizability
      ! (real and complex) and the plasma frequencies
      ! We now calculate the intraband contributions to the polarizability via
      ! the Drude-like model
      call calculate_drude_like_intraband_polarizability &
        ( omega_min, omega_max, n_omega, omega_plasma_sq(i_dielectric), &
          broaden_para(i_dielectric), intraband_polarizability(1,i_dielectric) )
      ! Once we have the polarizability (interband and intraband, real and
      ! complex,) all other quantities immediately follow
      call calculate_dielectric_from_polarizability &
        ( n_omega, broaden_type(i_dielectric), component_is_diagonal, &
          interband_polarizability(1,i_dielectric), &
          intraband_polarizability(1,i_dielectric), &
          dielectric(1,i_dielectric) )
      ! Note:  the output of the following only applies when we are working in a
      ! coordinate system where the dielectric tensor is diagonal
      if (component_is_diagonal) then
        call calculate_absorption_coeff( &
             omega_min, omega_max, n_omega, dielectric(1,i_dielectric), &
             absorption_coeff(1,i_dielectric))
      end if
    end do

    !--------------------------------------------------------------------------!
    ! FIFTH AND FINAL STEP:
    !           Output quantities and deallocate memory
    !------------------------------------------------------------------------- !
    if (myid == 0) then
      do i_dielectric = 1, n_plot_dielectric
        component_is_diagonal = direc_1(i_dielectric).eq.direc_2(i_dielectric)
        call output_optical_properties(i_dielectric, component_is_diagonal, &
             broaden_type(i_dielectric), &
             broaden_para(i_dielectric), &
             interband_polarizability(1,i_dielectric), &
             intraband_polarizability(1,i_dielectric), &
             dielectric(1,i_dielectric), absorption_coeff(1,i_dielectric),&
             omega_plasma_sq(i_dielectric), direc_1(i_dielectric), &
             direc_2(i_dielectric))
      end do
    end if

    call aims_deallocate(l_row_moment, "l_row_moment")
    call aims_deallocate(l_col_moment, "l_col_moment")
    call aims_deallocate(my_row_moment, "my_row_moment")
    call aims_deallocate(my_col_moment, "my_col_moment")

    call deallocate_mommat_TZ()
    call cleanup_KS_optical_properties()
    if (allocated(moment_TZ)) call aims_deallocate( moment_TZ, "moment_TZ" )
    if (allocated(interband_polarizability_reset)) &
         call aims_deallocate(interband_polarizability_reset, &
              "interband_polarizability_reset")
    if (calculate_perturbative_soc) then
      real_eigenvectors = real_eigenvectors_save
    end if

    call get_times(time_sync_postprocess, clock_time_sync_postprocess, &
         tot_time_sync_postprocess, tot_clock_time_sync_postprocess, .true.)
    call get_times(time_dielectric, clock_time_dielectric, tot_time_dielectric,&
         tot_clock_time_dielectric, .true. )
    ! -------------------------------------------------------------------------!
    ! Fin.
    ! -------------------------------------------------------------------------!
  end subroutine KS_dielectric_calculation
  !******

  subroutine KS_dielectric_calculation_imagfreq &
       (n_rows, n_cols, n_spin_in, n_states_in, KS_eigen, KS_vec, KS_vec_complex, occ_numbers, &
        chemical_potential, partition_tab, l_shell_max)
    ! PURPOSE
    !   Input:
    !   n_plot_dielectric
    !   broaden_type(n_plot_dielectric), broaden_para(n_plot_dielectric)
    !
    !   Wrapper function for calculating and outputting of the marcroscopic
    !   linear dielectric function in three steps
    !   0. ep1_in, ep2_in from character to integer (x->1,y->2,z->3);
    !      allocations
    !   1. function 'calculate_mommat_p0', real space integral of atom centered
    !      basis function i with the gradient of function j for in each cell
    !      (cell distance)
    !   2. function 'construct_overlap' 'fourier transform' of real space
    !      integrals by summing up the phases resulting from cell distances
    !      function is called once for the upper triangle and once for the lower
    !      triangle of the sparse matrix from first step
    !   3. function 'calc_moment_p0' basis transformation from atom centered
    !      basis to KS-orbitals
    !   4. function 'calc_dielectric_function' momentum matrix elements are
    !      summed up to yield the dielectric function
    !  USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use calculate_mommat_base, only: get_state_minmax_K, calculate_mommat_p0, &
        allocate_mommat_k, calc_moment_p0, clean_mommat, clean_mommat_final, &
        mommat_full_w_up, mommat_full_w_complex_up, work_ovl_mom, &
        mommat_full_w_low, mommat_full_w_complex_low
    use dimensions, only: n_basis, n_states, n_spin, n_k_points, &
        n_k_points_task, n_full_points, n_species, n_plot_dielectric, calculate_perturbative_soc
    use dimensions_soc, only: n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll
    use gw_para, only: dielec_func_imagfreq, n_full_freq, omega_full_grid
    use localorb_io, only: localorb_info
    use mpi_tasks, only: myid, n_tasks, aims_stop_coll, mpi_comm_global
    use pbc_lists, only: k_weights
    use runtime_choices, only: broaden_type, direc_1, direc_2, &
        calc_px_dielectric, calc_py_dielectric, calc_pz_dielectric, &
        broaden_para, use_local_index, use_scalapack, &
        real_eigenvectors, use_load_balancing, Emin, Emax
    use scalapack_wrapper, only: my_k_point, my_scalapack_id, &
        my_scalapack_comm_all
    use soc_utilities, only: convert_sr_to_soc_environment, &
        revert_soc_to_sr_environment
    use synchronize_mpi_basic, only: sync_real_number, sync_vector_complex
    implicit none
    !  ARGUMENTS
    integer,    intent(in) :: n_rows
    integer,    intent(in) :: n_cols
    integer,    intent(in) :: n_spin_in
    integer,    intent(in) :: n_states_in
    real*8 ,    intent(in) :: KS_eigen(n_states_in, n_spin, n_k_points)
    real*8,     intent(in) :: KS_vec(n_rows, n_cols, n_spin, n_k_points_task)
    complex*16, intent(in) :: &
         KS_vec_complex(n_rows, n_cols, n_spin, n_k_points_task)
    real*8,     intent(in) :: occ_numbers(n_states, n_spin,n_k_points)
    real*8,     intent(in) :: chemical_potential
    real*8,     intent(in), target :: partition_tab(n_full_points)
    integer,    intent(in) :: l_shell_max(n_species)
    !  INPUTS
    !    o KS_eigen
    !    o KS_vec/KS_vec_complex
    !    o occ_numbers
    !    o chemical_potential
    !    o partition_tab
    !    o l_shell_max
    !  OUTPUT
    !    None (modified module variables)
    !  AUTHOR
    !    William Huhn and Tong Zhu (Duke University), based on earlier code by
    !    Bjoern Bieniek (formerly FHI-Berlin)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    !  local variables
    character*128 :: info_str
    integer :: info
    integer :: n_state_min_in
    integer :: n_state_max_in
    integer :: n_states_window
    logical :: component_is_diagonal, real_eigenvectors_save
    integer :: ld_mommat
    !  counters
    integer :: i, i_k, i_k_point
    integer :: new_k_point
    integer :: coord_1, coord_2
    integer :: my_k_points(n_k_points)
    !TZ changes for output many direction dielectric function
    integer :: i_dielectric, i_omega

    complex*16, allocatable :: moment_TZ(:,:)

    character(*), parameter :: func = 'KS_dielectric_calculation_imagfreq'

    write(info_str,'(2X,A)') "Calculating macroscopic dielectric function &
                             &along the imaginary frequency axis ..."
    call localorb_info ( info_str )

    ! Should have been caught already in read_control, but just in case
    if (use_local_index) then
      call aims_stop_coll ('* Dielectric calculations currently do not support &
                           &use_local_index keyword.  Exiting', func)
    end if
    if (use_load_balancing) then
      call aims_stop_coll ('* Dielectric calculations currently do not support &
                           &load_balancing keyword.  Exiting', func)
    end if

    if (calculate_perturbative_soc) then
      ! Force the code to traverse the complex path
      ! in all cases, since SOC-perturbed eigenvectors
      ! are complex
      real_eigenvectors_save = real_eigenvectors
      real_eigenvectors      = .false.

      ! Some light error checking
      if (mod(n_basis_soc_coll,2) .eq. 1) then
        call aims_stop_coll('The number of collinear basis functions in SOC is &
                           &odd, exiting.', func)
      end if
      if (n_basis_soc_ncoll .gt. 0) then
        call aims_stop_coll('This subroutine does not support usage of &
                            &non-collinear basis functions in SOC, exiting.', &
                            func)
      end if
    end if

    ! Here is where various quantities related to parallelizing over k-points
    ! are calculated
    ! This should be rolled into a generic function at some point.
    ! Or better yet, be added as a global variable in pbc_list.
    if (n_k_points_task .gt. 0) then
      if (n_k_points .ge. n_tasks) then
        ! Undo the round-robin allocation
        i = 1
        do i_k_point = 1, n_k_points, 1
          if (myid ==  MOD(i_k_point, n_tasks) .and. myid <= n_k_points) then
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
          ! Note that in this case, myid 0 will never have a k-point
          my_k_points(1) = myid
        else
          ! This is used later to denote this MPI rank has no k-points assigned
          my_k_points(1) = 0
        end if
      ! The ScaLAPACK case.
      else
        my_k_points(1) = my_k_point
      end if
    end if

    ! -------------------------------------------------------------------------!
    ! FIRST STEP:
    !              allocate needed memory
    ! -------------------------------------------------------------------------!

    ! Determine the state indices corresponding to the user-provided energy
    ! window
    if (calculate_perturbative_soc) then
      call convert_sr_to_soc_environment ()
    end if
    write (info_str, '(2X,A)') "Recalculating HOMO and LUMO levels to &
                               &determine eigenstates needed for dielectric &
                               &calculations..."
    call localorb_info( info_str )
    write (info_str, *)

    call localorb_info( info_str )
    calc_px_dielectric = .true.
    calc_py_dielectric = .true.
    calc_pz_dielectric = .true.

    if (.not.allocated(broaden_para)) then
       call aims_allocate( broaden_para, n_plot_dielectric, "broaden_para" )
    end if
    if (.not.allocated(broaden_type)) then
       call aims_allocate( broaden_type, n_plot_dielectric, "broaden_type" )
    end if
    broaden_type(:) = 2
    broaden_para(:) = 0.01

    ! -------------------------------------------------------------------------!
    ! FIRST STEP:
    !              allocate needed memory
    ! -------------------------------------------------------------------------!
!  Emin and Emax have already been set in control.in. We reset them here to be safe.
!  The idea is to convert the entire spectrum.
    Emin=-1.e10
    Emax= 1.e10
    call get_state_minmax_K(KS_eigen, n_state_min_in, n_state_max_in)
    n_states_window = n_state_max_in - n_state_min_in + 1

    write (info_str, *) ""
    call localorb_info( info_str )
    write (info_str, '(2X,A,I9)') "Index of the first state to include in &
                                  &the calculations of the dielectric function along the imaginary axis : ", n_state_min_in
    call localorb_info( info_str )
    write (info_str, '(2X,A,I9)') "Index of the last state to include in &
                                  &the calculations of the dielectric function along the imaginary axis : ", n_state_max_in
    call localorb_info( info_str )
    write (info_str, *)
    call localorb_info( info_str )

    ! Initialize the ScaLAPACK environment for the moment matrix (which has a
    ! leading dimension of n_states_window)
    if (use_scalapack) then
      write(info_str, '(2X,A)') "Creating ScaLAPACK distribution for matrices&
                                 & with leading dimension n_states_window"
      call localorb_info( info_str )
      call aims_allocate( l_row_moment, n_states_window,   "l_row_moment" )
      call aims_allocate( l_col_moment, n_states_window,   "l_col_moment" )
      call aims_allocate( my_row_moment, n_states_window, "my_row_moment" )
      call aims_allocate( my_col_moment, n_states_window, "my_col_moment" )
      block_size_moment = 0
      call initialize_scalapack_descriptor_matrix( &
           n_states_window, n_states_window, block_size_moment, sc_desc_moment,&
           mb_moment, nb_moment, mxld_moment, mxcol_moment, n_my_rows_moment, &
           n_my_cols_moment, l_row_moment, l_col_moment, my_row_moment, &
           my_col_moment )
      write(info_str, *)
      call localorb_info( info_str )
    else
      call aims_allocate( l_row_moment, 1,   "l_row_moment" )
      call aims_allocate( l_col_moment, 1,   "l_col_moment" )
      call aims_allocate( my_row_moment, 1, "my_row_moment" )
      call aims_allocate( my_col_moment, 1, "my_col_moment" )
    end if

    call allocate_mommat_TZ( ld_mommat )
!    if (.not.allocated( moment_TZ)) then
!      call aims_allocate( moment_TZ, (((n_state_max_in-n_state_min_in+1)+1)* &
!           (n_state_max_in-n_state_min_in+1)/2)*((n_spin*(n_spin+1))/2), 3, &
!           "moment_TZ" )
!    end if
    if (calculate_perturbative_soc) then
      if(.not.allocated( moment_TZ )) then
        call aims_allocate(moment_TZ, &
             (n_states_window*(n_states_window+1))/2, 3, "+moment_TZ")
      end if
    else
      if(.not.allocated( moment_TZ )) then
        call aims_allocate(moment_TZ, &
             (n_states_window*(n_states_window+1))/2*n_spin_in, 3, &
             "+moment_TZ")
      end if
    end if

    if (.not.allocated(omega_plasma_sq)) then
      call aims_allocate(omega_plasma_sq, n_plot_dielectric, "omega_plasma_sq")
    end if

    if (.not.allocated(interband_polarizability_imagfreq)) then
      call aims_allocate( interband_polarizability_imagfreq, n_full_freq, &
           n_plot_dielectric, "interband_polarizability_imagfreq" )
    end if
    if (.not.allocated(direc_1)) then
      call aims_allocate( direc_1, n_plot_dielectric, "direc_1" )
    end if
    if (.not.allocated(direc_2)) then
      call aims_allocate( direc_2, n_plot_dielectric, "direc_2" )
    end if

    direc_1(1) = 1
    direc_1(2) = 2
    direc_1(3) = 3
    direc_2(1) = 1
    direc_2(2) = 2
    direc_2(3) = 3

    moment_TZ = 0.0
    !--------------------------------------------------------------------------!
    ! SECOND STEP:
    !              calculate real-space momentum matrix elements needed based on
    !              tensor components requested by the user (i.e. if they request
    !              the xy component, calculate the p_x and p_y matrix elements
    !              but not p_z)
    !--------------------------------------------------------------------------!
    write(info_str, *)
    call localorb_info ( info_str )
    if (calc_px_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_x matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0( partition_tab, l_shell_max, &
           mommat_full_x_up, mommat_full_x_low,1)
      write(info_str, *)
      call localorb_info ( info_str )
    end if

    if (calc_py_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_y matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0(partition_tab, l_shell_max, &
           mommat_full_y_up, mommat_full_y_low,2)
      write(info_str, *)
      call localorb_info ( info_str )
    end if

    if (calc_pz_dielectric) then
      write(info_str,'(2X,A)') "Calculating p_z matrix elements..."
      call localorb_info ( info_str )
      call calculate_mommat_p0 ( partition_tab, l_shell_max, &
           mommat_full_z_up, mommat_full_z_low,3)
      write(info_str, *)
      call localorb_info ( info_str )
    end if

    !--------------------------------------------------------------------------!
    ! THIRD STEP:
    !             Calculate quantities requiring summation of k-points.
    !             At this time, these would be the plasma frequency (squared)
    !             and the interband contribution to the polarizability
    !--------------------------------------------------------------------------!
    write(info_str,'(2X,A)') "Calculating dielectric tensor elements..."
    call localorb_info ( info_str )

    omega_plasma_sq = 0.0d0
    interband_polarizability_imagfreq  = (0.0d0,0.0d0)

    do i_k = 1,n_k_points_task
      new_k_point = my_k_points(i_k)

!      if (myid ==  modulo(new_k_point, n_tasks) .and. myid <= n_k_points ) then
        ! TZ: for the dense k-grid needed, the dielectric function part are now
        ! only support the lapack case. not work for the scalapack. using
        ! scalpack would lead to some un-predicted results.
!        call allocate_mommat_k()
        ! TZ:  calculate possible moment matrix for needed direction
        !----------------------------------------------------------------------!
        ! THIRD STEP, PART 1:
        !            Calculate the momentum matrix elements at this k-point
        !            (e.g. in the Bloch basis representation)
        !----------------------------------------------------------------------!
        ! WPH:  these moment matrices are O(n_state^2) and will lead
        !       to memory issues down the road.
        if (calc_px_dielectric) then
!          call construct_overlap(mommat_full_x_up, mommat_full_w_up,&
!                                 mommat_full_w_complex_up, new_k_point,&
!                                 work_ovl_mom)
!          call construct_overlap(mommat_full_x_low, mommat_full_w_low, &
!                                 mommat_full_w_complex_low, new_k_point,&
!                                 work_ovl_mom)
!          call calc_moment_p0(moment_TZ(:,1),mommat_full_w_up, &
!                     mommat_full_w_low, mommat_full_w_complex_up, &
!                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
!                     KS_vec_complex(:,:,:,i_k) , new_k_point, 1, &
!                     n_state_min_in, n_state_max_in)
          call calculate_moment_matrix(ld_mommat, mommat_full_x_up, &
                 mommat_full_x_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
                 KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
                 n_states_window, moment_TZ(:,1))
        end if

        if (calc_py_dielectric) then
!          call construct_overlap(mommat_full_y_up, mommat_full_w_up,&
!                                 mommat_full_w_complex_up, new_k_point,&
!                                 work_ovl_mom)
!          call construct_overlap(mommat_full_y_low, mommat_full_w_low, &
!                                 mommat_full_w_complex_low, new_k_point,&
!                                 work_ovl_mom)
!          call calc_moment_p0(moment_TZ(:,2),mommat_full_w_up, &
!                     mommat_full_w_low, mommat_full_w_complex_up, &
!                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
!                     KS_vec_complex(:,:,:,i_k) , new_k_point, 2, &
!                     n_state_min_in, n_state_max_in)
          call calculate_moment_matrix(ld_mommat, mommat_full_y_up, &
                  mommat_full_y_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
                  KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
                  n_states_window, moment_TZ(:,2))
        end if

        if (calc_pz_dielectric) then
!          call construct_overlap(mommat_full_z_up, mommat_full_w_up,&
!                                 mommat_full_w_complex_up, new_k_point,&
!                                 work_ovl_mom)
!          call construct_overlap(mommat_full_z_low, mommat_full_w_low, &
!                                 mommat_full_w_complex_low, new_k_point,&
!                                 work_ovl_mom)
!          call calc_moment_p0(moment_TZ(:,3),mommat_full_w_up, &
!                     mommat_full_w_low, mommat_full_w_complex_up, &
!                     mommat_full_w_complex_low, KS_vec(:,:,:,i_k) ,&
!                     KS_vec_complex(:,:,:,i_k) , new_k_point, 3, &
!                     n_state_min_in, n_state_max_in)
          call calculate_moment_matrix(ld_mommat, mommat_full_z_up,&
                  mommat_full_z_low, n_rows, n_cols, KS_vec(:,:,:,i_k), &
                  KS_vec_complex(:,:,:,i_k), new_k_point, n_state_min_in, &
                  n_states_window, moment_TZ(:,3) )
        end if
        ! ---------------------------------------------------------------------!
        ! THIRD STEP, PART 2:
        !             Calculate the plasma frequency (squared) and the interband
        !             contribuation to the polarizability for all requested
        !             tensor components and broadening functions
        ! ---------------------------------------------------------------------!
        do i_dielectric = 1, n_plot_dielectric
          call calculate_plasma_frequency_sq_imagfreq( &
               moment_TZ(:,direc_1(i_dielectric)), &
               moment_TZ(:,direc_2(i_dielectric)), &
               omega_plasma_sq(i_dielectric), chemical_potential, &
               KS_eigen(:,:,new_k_point), k_weights(new_k_point), &
               broaden_para(i_dielectric), n_state_min_in, n_state_max_in )
          call calculate_interband_polarizability_imagfreq( &
               moment_TZ(:,direc_1(i_dielectric)), &
               moment_TZ(:,direc_1(i_dielectric)),&
               interband_polarizability_imagfreq(:,i_dielectric), &
               chemical_potential, KS_eigen(:,:,new_k_point), &
               k_weights(new_k_point), broaden_type(i_dielectric), &
               broaden_para(i_dielectric), n_state_min_in, n_state_max_in)
!              if(i_dielectric.eq.1) then
!               write(6,'(2X,A,I4, 2f16.8, I4, 2f16.8)') "new_k_point, k_weights(new_k_point),  myid, interpolar", &
!                   new_k_point, KS_eigen(1,1,new_k_point), KS_eigen(19,1,new_k_point), myid, interband_polarizability_imagfreq(1,i_dielectric)
!              endif
        end do
        call clean_mommat()
!      end if
    end do

    !--------------------------------------------------------------------------!
    ! THIRD STEP, PART 3:
    !             Each MPI task still has only a partial k-point sum, now sum
    !             over all MPI tasks to get the full k-point sum
    !--------------------------------------------------------------------------!
!    write(6,*) "use_scalapack", use_scalapack
    if (use_scalapack .and. my_scalapack_id.ne.0) then
      interband_polarizability_imagfreq = (0.0d0,0.0d0)
      omega_plasma_sq = 0.0d0
    end if
    do i_dielectric = 1, n_plot_dielectric
      if (.not. use_local_index) &
           call sync_real_number(omega_plasma_sq(i_dielectric))
      if (.not. use_local_index) &
           call sync_vector_complex( &
                interband_polarizability_imagfreq(:,i_dielectric), n_full_freq )
    end do

    write(info_str, *)
    call localorb_info ( info_str )

    ! defined in module "gw_para", the dielectric function at imaginary axis
    ! FIXME: needs to be modified to include the intraband contributions and
    ! considering the non-isotropic situation.
    dielec_func_imagfreq(:) = (interband_polarizability_imagfreq(:,1) + &
         interband_polarizability_imagfreq(:,2) + &
         interband_polarizability_imagfreq(:,3))/3.d0
    !--------------------------------------------------------------------------!
    ! FOURTH AND FINAL STEP:
    !           Output quantities and deallocate memory
    !--------------------------------------------------------------------------!

!    open(unit=8, file='dielec_imagfreq.dat',ACTION='WRITE')
!!    write(8,'(A)')"# \omega(Ha),   \epsilon(i\omega) "
!    do i_omega = 1,n_full_freq
!       WRITE (8,'(f18.8, 4f20.10)') omega_full_grid(i_omega), dielec_func_imagfreq(i_omega), interband_polarizability_imagfreq(:,1)
!    end do
!    close(unit=8)

    call deallocate_mommat_TZ ()
    if (allocated(broaden_para)) &
         call aims_deallocate(broaden_para, "broaden_para")
    if (allocated(broaden_type)) &
         call aims_deallocate(broaden_type, "broaden_type")
    if (allocated(interband_polarizability_imagfreq)) &
         call aims_deallocate( interband_polarizability_imagfreq, &
              "interband_polarizability_imagfreq" )
    if (allocated(omega_plasma_sq)) call aims_deallocate(omega_plasma_sq, &
         "omega_plasma_sq" )
    if (allocated(direc_1)) call aims_deallocate(direc_1, "direc_1")
    if (allocated(direc_2)) call aims_deallocate(direc_2, "direc_2")
    if (allocated(moment_TZ)) call aims_deallocate(moment_TZ, "moment_TZ")

    if (calculate_perturbative_soc) then
      real_eigenvectors = real_eigenvectors_save
    end if

    ! -------------------------------------------------------------------------!
    ! Fin.
    ! -------------------------------------------------------------------------!
  end subroutine KS_dielectric_calculation_imagfreq
  !******

  !!!!!!!!!!!!!!!!!!!!!!!
  ! UTILITY SUBROUTINES !
  !!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate_moment_matrix(ld_mommat, mommat_up, mommat_low, &
       n_rows, n_cols, KS_vec, KS_vec_complex, k_point_global, &
       n_state_min_in, n_states_window, moment)
    !  PURPOSE
    !    Wrapper function for calculating momentum matrix element between
    !    eigenstates for a given p_i and a specified k-point.
    !  USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use calculate_mommat_base, only: calc_moment_p0_SOC
    use dimensions, only: calculate_perturbative_soc, &
        n_hamiltonian_matrix_size, n_basis, n_centers_basis_I, n_spin
    use load_balancing,only: batch_perm, n_bp_integ, &
        init_comm_full_local_matrix
    use mpi_tasks, only: myid, aims_stop
    use runtime_choices, only: real_eigenvectors, use_scalapack, &
        use_local_index, use_load_balancing
    use scalapack_wrapper, only: construct_hamiltonian_like_matrix_scalapack, &
        construct_hamiltonian_like_matrix_zero_diag_scalapack, &
        set_sparse_local_matrix_scalapack_generic, &
        set_full_local_matrix_scalapack_generic, mxld, mxcol, sc_desc
    implicit none
    !  ARGUMENTS
    integer,    intent(in)  :: ld_mommat
    real*8,     intent(in)  :: mommat_up(ld_mommat)
    real*8,     intent(in)  :: mommat_low(ld_mommat)
    integer,    intent(in)  :: n_rows
    integer,    intent(in)  :: n_cols
    real*8,     intent(in)  :: KS_vec(:,:,:)
    complex*16, intent(in)  :: KS_vec_complex(:,:,:)
    integer,    intent(in)  :: k_point_global
    integer,    intent(in)  :: n_state_min_in
    integer,    intent(in)  :: n_states_window
    complex*16, intent(out) :: moment(:)
    !  INPUT
    !    o ld_mommat - Number of elements in mommat matrices.  Used for error
    !      checking to ensure that matrix size matches the code path traversed
    !    o mommat_up - Upper half of real-space momentum matrix element between
    !      basis functions for the current p_i operator
    !    o mommat_low - Lower half of real-space momentum matrix element between
    !      basis functions for the current p_i operator
    !    o n_rows - Numbers of rows for the eigenvectors.  Should be the number
    !      of basis functions for a LAPACK-style eigenvector and the appropriate
    !      mxld-like variable for a ScaLAPACK eigenvector.
    !    o n_cols - Numbers of columns for the eigenvector.  Should be the
    !      number of states for a a LAPACK-style eigenvector and the appropriate
    !      mxcol-like variable for a ScaLAPACK eigenvector.
    !    o KS_vec - Real eigenvector for the current k-point.
    !    o KS_vec_complex - Complex eigenvector for the current k-point.
    !    o k_point_global - The k-point index for the current k-point, i.e. what
    !      indexes the eigennvalues and occupation numbers.
    !    o n_state_min_in - The minimum state to include in the energy window
    !    o n_states_window - The number of states to include in the energy
    !      window
    !  OUTPUT
    !    o moment - The momentum matrix element between eigenstates at the
    !      current k-point for the current p_i operator.  Stored in a
    !      lower-triangular format.
    !  AUTHOR
    !    William Huhn (Duke University)
    !  HISTORY
    !    September 2017 - Created
    !  SOURCE
    integer :: i_basis_1, i_basis_2, i_index
    integer :: n_spin_save
    real*8,     allocatable :: mommat_temp(:)
    real*8,     allocatable :: mommat_full(:,:)
    complex*16, allocatable :: mommat_temp_complex(:)
    complex*16, allocatable :: mommat_full_complex(:,:)
    real*8,     allocatable :: work_ovl_mom(:,:)

    character(*), parameter :: func = 'calculate_moment_matrix'

    ! Check to ensure that reported matrix size corresponds to the current code
    ! path.  I keep going back and forth on whether this should be rigidly
    ! enforced, or if we should just pass the block of memory through with no
    ! check...
    if (use_load_balancing) then
      if (ld_mommat .ne. batch_perm(n_bp_integ)%n_local_matrix_size) then
        call aims_stop("* Input momentum matrices have the incorrect size, &
                       &exiting.", func)
      end if
    else
      if (ld_mommat .ne. n_hamiltonian_matrix_size) then
        call aims_stop("* Input momentum matrices have the incorrect size, &
                       &exiting.", func)
      end if
    end if

    ! Allocate temporary arrays for various code paths
    ! TODO:  We should be able to eliminate the mommat_temp_* matrices by
    ! performing in-place construction and accumulation for the "low" matrix.
    if (use_scalapack) then
      if (real_eigenvectors) then
        call aims_allocate(mommat_temp, mxld*mxcol, "mommat_temp")
        call aims_allocate(mommat_full, mxld,mxcol, "mommat_full")
        call aims_allocate(mommat_temp_complex, 1, "mommat_temp_complex")
        call aims_allocate(mommat_full_complex, 1,1, "mommat_full_complex")
      else
        call aims_allocate(mommat_temp, 1, "mommat_temp")
        call aims_allocate(mommat_full, 1,1, "mommat_full")
        call aims_allocate(mommat_temp_complex, mxld*mxcol, &
             "mommat_temp_complex")
        call aims_allocate(mommat_full_complex, mxld,mxcol, &
             "mommat_full_complex")
      end if
      call aims_allocate(work_ovl_mom, 1,1, "work_ovl_mom" )
    else ! LAPACK
      if (real_eigenvectors) then
        call aims_allocate(mommat_temp, (n_basis*(n_basis+1))/2, "mommat_temp")
        call aims_allocate(mommat_full, n_basis, n_basis, "mommat_full" )
        call aims_allocate(mommat_temp_complex, 1, "mommat_temp_complex")
        call aims_allocate(mommat_full_complex, 1, 1, "mommat_full_complex")
      else
        call aims_allocate(mommat_temp, 1, "mommat_temp")
        call aims_allocate(mommat_full, 1,1, "mommat_full")
        call aims_allocate(mommat_temp_complex, (n_basis*(n_basis+1))/2, &
             "mommat_temp_complex")
        call aims_allocate(mommat_full_complex, n_basis, n_basis, &
             "mommat_full_complex")
      end if
      call aims_allocate(work_ovl_mom, n_centers_basis_I, n_centers_basis_I, &
           "work_ovl_mom")
    end if

    ! Initialize to be on the safe side
    mommat_full = 0.0d0
    mommat_temp = 0.0d0
    mommat_full_complex = (0.0d0, 0.0d0)
    mommat_temp_complex = (0.0d0, 0.0d0)

    ! Convert momentum matrix from real-space basis to Bloch basis
    ! In this section, we have to dance around the fact that aims generally
    ! expects upper-triangular matrices, but we are storing the upper and lower
    ! triangles independently, both in an upper-triangular form.  Accordingly,
    ! in all code paths, the general strategy is to construct the upper triangle
    ! directly and put it immediately in the final "full" matrix, then construct
    ! the lower triangle as an upper-triangular matrix and conjugate it, adding
    ! the conjugate to the "full" matrix.
    if (use_scalapack) then
      ! The "up" portion of the matrix are already in the correct positions,
      ! insert them into the full matrix immediately
      if (use_local_index) then
        if (use_load_balancing) then
          ! THIS DOES NOT WORK!  See header for more details.
          call init_comm_full_local_matrix( &
               batch_perm(n_bp_integ)%n_basis_local, &
               batch_perm(n_bp_integ)%i_basis_local )
          call set_full_local_matrix_scalapack_generic(mommat_up, mommat_full, &
               mommat_full_complex )
          call set_full_local_matrix_scalapack_generic(mommat_low, mommat_temp,&
               mommat_temp_complex )
        else
          ! THIS (probably) DOES NOT WORK!  See header for more details.
          call set_sparse_local_matrix_scalapack_generic(mommat_up, &
               mommat_full, mommat_full_complex)
          call set_sparse_local_matrix_scalapack_generic(mommat_low, &
               mommat_temp, mommat_temp_complex)
        end if
      else
        ! construct_hamiltonian_like_matrix_scalapack loops over both spin
        ! channels, here we don't have a spin channel so hack it away
        n_spin_save = n_spin
        n_spin      = 1
        call construct_hamiltonian_like_matrix_zero_diag_scalapack(mommat_up, &
             mommat_full, mommat_full_complex)
        call construct_hamiltonian_like_matrix_scalapack(mommat_low, &
             mommat_temp, mommat_temp_complex)
        n_spin = n_spin_save
      end if ! use_load_balancing

      ! Now conjugate the "low" portion and add it to the full matrix
      ! WPH: Ordinarily, we would be doubling the diagonal here, but the
      ! gradient operator is anti-hermitian so the diagonal should be zero.
      if (real_eigenvectors) then
        call pdtran(n_basis, n_basis, &
                    1.0d0, &
                    mommat_temp, 1,1,sc_desc, &
                    1.0d0, &
                    mommat_full, 1,1,sc_desc)
      else
        call pztranc(n_basis, n_basis, &
                     (1.0d0,0.0d0), &
                     mommat_temp_complex, 1,1,sc_desc, &
                     (1.0d0,0.0d0), &
                     mommat_full_complex, 1,1,sc_desc)
      end if
    else ! .not.use_local_index
      ! Construct the upper half of the matrix
      call construct_overlap(mommat_up,  mommat_temp,  mommat_temp_complex, &
           k_point_global, work_ovl_mom)
      i_index = 0
      do i_basis_2 = 1, n_basis
        do i_basis_1 = 1, i_basis_2
          i_index=i_index+1
          if (real_eigenvectors) then
            mommat_full(i_basis_1, i_basis_2) = mommat_temp(i_index)
          else
            mommat_full_complex(i_basis_1, i_basis_2) = &
                 mommat_temp_complex(i_index)
          end if
        end do
      end do

      ! Zero out the diagonal, as this will be set by the lower half
      if (real_eigenvectors) then
        do i_basis_1 = 1, n_basis
           mommat_full(i_basis_1, i_basis_1) = 0.0d0
        end do
      else
        do i_basis_1 = 1, n_basis
           mommat_full_complex(i_basis_1, i_basis_1) = (0.0d0, 0.0d0)
        end do
      end if

      ! Now construct the lower half
      call construct_overlap(mommat_low, mommat_temp, mommat_temp_complex, &
           k_point_global, work_ovl_mom)
      i_index = 0
      do i_basis_2 = 1, n_basis
        do i_basis_1 = 1, i_basis_2
          i_index=i_index+1
          if (real_eigenvectors) then
            mommat_full(i_basis_2, i_basis_1) = mommat_temp(i_index)
          else
            mommat_full_complex(i_basis_2, i_basis_1) = &
                 conjg(mommat_temp_complex(i_index))
          end if
        end do
      end do
    end if ! use_local_index

    ! Deallocate unneeded arrays early to give next subroutine some breathing
    ! room
    call aims_deallocate(mommat_temp, "mommat_temp")
    call aims_deallocate(mommat_temp_complex, "mommat_temp_complex")
    call aims_deallocate(work_ovl_mom, "work_ovl_mom")

    ! Calculate matrix elements of the momentum matrix for the eigenvectors
    if (calculate_perturbative_soc) then
      call calc_moment_w0_soc( n_rows, n_cols, KS_vec_complex, &
                               mommat_full_complex, &
                               n_state_min_in, n_states_window, moment )
    else
      call calc_moment_w0( n_rows, n_cols, KS_vec, KS_vec_complex, &
                           mommat_full, mommat_full_complex, &
                           n_state_min_in, n_states_window, moment )
    end if

    ! Deallocate final arrays and exit
    call aims_deallocate(mommat_full, "mommat_full")
    call aims_deallocate(mommat_full_complex, "mommat_full_complex")
  end subroutine calculate_moment_matrix
  !******

  subroutine calc_moment_w0( n_rows, n_cols, KS_vec, KS_vec_complex, &
                             mommat_full,  mommat_full_complex, &
                             n_state_min_in, n_states_window, moment )
    !  PURPOSE
    !    Given momentum matrix elements in the Bloch basis for a given p_i and
    !    at a specified k-point, calculate momentum matrix elements between
    !    eigenstates
    !
    !    This subroutine is intended for scalar-relativistic eigenvectors.
    !  USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use dimensions, only: n_basis, n_states, n_spin
    use mpi_tasks, only: myid, aims_stop
    use runtime_choices, only: real_eigenvectors, use_scalapack
    use scalapack_wrapper, only: sc_desc, mxld, mxcol, my_scalapack_comm_all, &
        l_col, l_row
    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none
    !  ARGUMENTS
    integer,    intent(in)  :: n_rows
    integer,    intent(in)  :: n_cols
    real*8,     intent(in)  :: KS_vec(n_rows, n_cols, n_spin)
    complex*16, intent(in)  :: KS_vec_complex(n_rows, n_cols, n_spin)
    real*8,     intent(in)  :: mommat_full(:,:)
    complex*16, intent(in)  :: mommat_full_complex(:,:)
    integer,    intent(in)  :: n_state_min_in
    integer,    intent(in)  :: n_states_window
    complex*16, intent(out) :: &
         moment((n_states_window*(n_states_window+1))/2*n_spin)
    !  INPUT
    !    o n_rows - Numbers of rows for the eigenvectors.  Should be the number
    !      of basis functions for a LAPACK-style eigenvector and the appropriate
    !      mxld-like variable for a ScaLAPACK eigenvector.
    !    o n_cols - Numbers of columns for the eigenvector.  Should be the
    !      number of states for a LAPACK-style eigenvector and the appropriate
    !      mxcol-like variable for a ScaLAPACK eigenvector.
    !    o KS_vec - Real (scalar-relativistic) eigenvector for the current
    !      k-point.
    !    o KS_vec_complex - Complex (scalar-relativistic) eigenvector for the
    !      current k-point.
    !    o mommat_full - Momentum matrix elements for a given p_i in the Bloch
    !      basis at the current k-point (real version)
    !    o mommat_full_complex - Momentum matrix elements for a given p_i in
    !      the Bloch basis at the current k-point (complex version)
    !    o n_state_min_in - The minimum state to include in the energy window
    !    o n_states_window - The number of states to include in the energy
    !      window
    !  OUTPUT
    !    o moment - The momentum matrix element between eigenstates at the
    !      current k-point for the current p_i operator.  Stored in a
    !      lower-triangular format.
    !  AUTHOR
    !    William Huhn (Duke University)
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !    the terms and conditions of the respective license agreement."
    !  HISTORY
    !    September 2017 - Created.
    !  SOURCE
    integer :: i_spin
    integer :: i_state, j_state, i_off
    integer :: i_basis_1, i_basis_2
    integer :: i_index, i_col, i_col_global, i_row, i_row_global
    ! Intermediate arrays
    real*8,     allocatable :: temp(:,:)
    complex*16, allocatable :: temp_complex(:,:)
    real*8,     allocatable :: moment_temp(:,:)
    complex*16, allocatable :: moment_temp_complex(:,:)

    character*300           :: info_str

    character(*), parameter :: func = 'calc_moment_w0'

    ! TODO: The moment_temp_* matrices are needed to interface to the
    ! subroutines that calculate the plasma frequency and interband
    ! polarizability.  When those subroutines are updated to support ScaLAPACK,
    ! the matrices can be safely eliminated
    if (use_scalapack) then
      if (n_rows .ne. mxld .or. n_cols .ne. mxcol) then
        write(info_str,'(1X,A)') '* Provided eigenvector does not appear to be &
             &in BLACS format.  Exiting.'
        call aims_stop( info_str, func )
      end if
      if (real_eigenvectors) then
        call aims_allocate(temp, mxld, mxcol, "temp")
        call aims_allocate(moment_temp, mxld_moment, mxcol_moment, "moment")
        call aims_allocate(temp_complex, 1, 1, "temp_complex")
        call aims_allocate(moment_temp_complex, 1, 1, "moment_complex")
      else
        call aims_allocate(temp, 1, 1, "temp")
        call aims_allocate(moment_temp, 1, 1, "moment")
        call aims_allocate(temp_complex, mxld, mxcol, "temp_complex")
        call aims_allocate(moment_temp_complex, mxld_moment, mxcol_moment, &
             "moment_complex" )
      end if
    else
      if (n_rows .ne. n_basis .or. n_cols .ne. n_states) then
        write(info_str,'(1X,A)') '* Provided eigenvector does not appear to be &
             &in full format.  Exiting.'
        call aims_stop( info_str, func )
      end if
      if (real_eigenvectors) then
        call aims_allocate(temp, n_basis, n_states_window, "temp")
        call aims_allocate(moment_temp, n_states_window, n_states_window, &
             "moment")
        call aims_allocate(temp_complex, 1, 1, "temp_complex")
        call aims_allocate(moment_temp_complex, 1, 1, "moment_complex")
      else
        call aims_allocate(temp, 1, 1, "temp")
        call aims_allocate(moment_temp, 1, 1, "moment")
        call aims_allocate(temp_complex, n_basis, n_states_window, &
             "temp_complex")
        call aims_allocate(moment_temp_complex, n_states_window, &
             n_states_window, "moment_complex")
      end if
    end if

    moment = (0.0d0, 0.0d0)
    do i_spin = 1, n_spin
      if (use_scalapack) then
        if (real_eigenvectors) then
          ! Calculate momentum matrix elements between eigenstates
          call pdgemm( 'N', 'N', &
                       n_basis, n_states_window, n_basis, &
                       1.0d0, &
                       mommat_full, 1,1,sc_desc, &
                       KS_vec(1,1,i_spin), 1,n_state_min_in,sc_desc, &
                       0.0d0, &
                       temp, 1,1,sc_desc )
          call pdgemm('T', 'N', &
                       n_states_window, n_states_window, n_basis, &
                       1.0d0, &
                       KS_vec(1,1,i_spin), 1,n_state_min_in,sc_desc, &
                       temp, 1,1,sc_desc, &
                       0.0d0, &
                       moment_temp, 1,1,sc_desc_moment )

          ! Index results back into lower-triangular form expected by rest of
          ! code
          i_off = (n_states_window*(n_states_window+1))/2 * (i_spin-1)
          ! Here we must inline what should be a utility function converting
          ! BLACS to a lower-triangular matrix, due to the real-to-complex
          ! conversion
          do i_row = 1, n_my_rows_moment
            do i_col = 1, n_my_cols_moment
              i_row_global = my_row_moment(i_row)
              i_col_global = my_col_moment(i_col)
              if (i_row_global > 0 .and. i_col_global > 0 .and. &
                  (i_col_global .le. i_row_global) .and. &
                  (i_col_global .le. n_states_window) ) then
                i_index = (n_states_window * (n_states_window+1))/2 &
                     - ((n_states_window-i_col_global+1) * &
                    ((n_states_window-i_col_global+1)+1))/2
                i_index = i_index + (i_row_global - i_col_global + 1)

                moment(i_index + i_off) = dcmplx(moment_temp(i_row,i_col))
              end if
            end do ! i_col
          end do ! i_row
        else ! .not.real_eigenvectors
          ! Calculate momentum matrix elements between eigenstates
          call pzgemm( 'N', 'N', &
                       n_basis, n_states_window, n_basis, &
                       (1.0d0,0.0d0), &
                       mommat_full_complex, 1,1,sc_desc, &
                       KS_vec_complex(1,1,i_spin), 1,n_state_min_in,sc_desc, &
                       (0.0d0,0.0d0), &
                       temp_complex, 1,1,sc_desc )
          call pzgemm( 'C', 'N', &
                       n_states_window, n_states_window, n_basis, &
                       (1.0d0,0.0d0), &
                       KS_vec_complex(1,1,i_spin), 1,n_state_min_in,sc_desc, &
                       temp_complex, 1,1,sc_desc, &
                       (0.0d0,0.0d0), &
                       moment_temp_complex, 1,1,sc_desc_moment )

          ! Index results back into lower-triangular form expected by rest of
          ! code
          ! See comments in main subroutine for plan to eliminate this
          ! conversion from BLACS to lower-triangular
          i_off = (n_states_window*(n_states_window+1))/2 * (i_spin-1)
          call convert_local_to_global_low_tri_matrix_complex( &
               n_my_rows_moment, n_my_cols_moment, my_row_moment, &
               my_col_moment, mxld_moment, mxcol_moment, moment_temp_complex, &
               n_states_window, moment(1+i_off) )
        end if
      else ! .not.use_scalapack
        if (real_eigenvectors) then
          ! Calculate momentum matrix elements between eigenstates
          call dgemm( 'N', 'N', &
                      n_basis, n_states_window, n_basis, &
                      1.0d0, &
                      mommat_full, n_basis, &
                      KS_vec(1,n_state_min_in,i_spin), n_basis, &
                      0.0d0, &
                      temp, n_basis )
          call dgemm( 'T', 'N', &
                      n_states_window, n_states_window, n_basis, &
                      1.0d0, &
                      KS_vec(1,n_state_min_in,i_spin), n_basis, &
                      temp, n_basis, &
                      0.0d0, &
                      moment_temp, n_states_window )

          ! Index results back into lower-triangular form expected by rest of
          ! code
          i_off = (n_states_window*(n_states_window+1))/2 * (i_spin-1)
          i_index = i_off
          do i_state = 1, n_states_window
            do j_state = i_state, n_states_window
              i_index = i_index + 1
              moment(i_index) = dcmplx(moment_temp(j_state,i_state))
            end do
          end do
        else ! .not.real_eigenvectors
          ! Calculate momentum matrix elements between eigenstates
          call zgemm( 'N', 'N', &
                      n_basis, n_states_window, n_basis, &
                      (1.0d0,0.0d0), &
                      mommat_full_complex, n_basis, &
                      KS_vec_complex(1,n_state_min_in,i_spin), n_basis, &
                      (0.0d0,0.0d0), &
                      temp_complex, n_basis )
          call zgemm( 'C', 'N', &
                      n_states_window, n_states_window, n_basis, &
                      (1.0d0,0.0d0), &
                      KS_vec_complex(1,n_state_min_in,i_spin), n_basis, &
                      temp_complex, n_basis, &
                      (0.0d0,0.0d0), &
                      moment_temp_complex, n_states_window )

          ! Index results back into lower-triangular form expected by rest of
          ! code
          i_off = (n_states_window*(n_states_window+1))/2 * (i_spin-1)
          i_index = i_off
          do i_state = 1, n_states_window
            do j_state = i_state, n_states_window
              i_index = i_index + 1
              moment(i_index) = moment_temp_complex(j_state,i_state)
            end do
          end do
        end if ! real_eigenvectors
      end if ! use_scalapack
    end do ! i_spin

    ! See comments in main subroutine for plan to eliminate this conversion from
    ! BLACS to lower-triangular
    if (use_scalapack) then
      call sync_vector_complex(moment, &
           ((n_states_window*(n_states_window+1))/2)*n_spin, &
           my_scalapack_comm_all )
    end if

    call aims_deallocate( temp, "temp" )
    call aims_deallocate( temp_complex, "temp_complex" )
    call aims_deallocate( moment_temp, "moment_temp" )
    call aims_deallocate( moment_temp_complex, "moment_temp_complex" )
  end subroutine calc_moment_w0
  !******

  subroutine calc_moment_w0_soc(n_rows, n_cols, KS_vec_soc, &
        mommat_full_complex, n_state_min_in, n_states_window, moment)
    !  PURPOSE
    !    Given momentum matrix elements in the Bloch basis for a given p_i and
    !    at a specied k-point, calculate momentum matrix elements between
    !    eigenstates
    !
    !    This subroutine is inteneded for SOC-perturbed eigenvectors.
    !  USES
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use dimensions, only: n_basis
    use dimensions_soc, only: n_basis_soc, n_states_soc
    use mpi_tasks, only: myid, aims_stop
    use runtime_choices, only: use_scalapack
    use scalapack_soc, only: sc_desc_soc_vec, mxld_soc_vec, mxcol_soc_vec, &
        l_col_soc_vec, l_row_soc_vec
    use scalapack_wrapper, only: sc_desc, mxld, mxcol, my_scalapack_comm_all
    use synchronize_mpi_basic, only: sync_vector_complex
    implicit none
    !  ARGUMENTS
    integer,    intent(in)  :: n_rows
    integer,    intent(in)  :: n_cols
    complex*16, intent(in)  :: KS_vec_soc(n_rows, n_cols)
    complex*16, intent(in)  :: mommat_full_complex(:,:)
    integer,    intent(in)  :: n_state_min_in
    integer,    intent(in)  :: n_states_window
    complex*16, intent(out) :: moment((n_states_window*(n_states_window+1))/2)
    !  INPUT
    !    o n_rows - Numbers of rows for the eigenvectors.  Should be the number
    !      of basis functions for a LAPACK-style eigenvector and the appropriate
    !      mxld-like variable for a ScaLAPACK eigenvector.
    !    o n_cols - Numbers of columns for the eigenvector.  Should be the
    !      number of states for a a LAPACK-style eigenvector and the appropriate
    !      mxcol-like variable for a ScaLAPACK eigenvector.
    !    o KS_vec_soc - Complex (spin-orbit-coupled) eigenvector for the current
    !      k-point.
    !    o mommat_full_complex - Momentum matrix elements for a given p_i in the
    !      Bloch basis at the current k-point (always complex for SOC)
    !    o n_state_min_in - The minimum state to include in the energy window
    !    o n_states_window - The number of states to include in the energy
    !      window
    !  OUTPUT
    !    o moment - The momentum matrix element between eigenstates at the
    !      current k-point for the current p_i operator.  Stored in a
    !      lower-triangular format.
    !  AUTHOR
    !    William Huhn (Duke University)
    !  COPYRIGHT
    !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !    the terms and conditions of the respective license agreement."
    !  HISTORY
    !    September 2017 - Created.
    !  SOURCE
    integer :: i_state, j_state, i_index, i_spin, i_off, basis_off

    ! Intermediate arrays
    complex*16, allocatable :: temp_complex(:,:)
    complex*16, allocatable :: moment_temp_complex(:,:)

    character*300 :: info_str

    character(*), parameter :: func = 'calc_moment_w0_soc'

    ! TODO: The moment_temp_* matrices are needed to interface to the
    ! subroutines that calculate the plasma frequency and interband
    ! polarizability.  When those subroutines are updated to support
    ! LAPACK/ScaLAPACK, the matrices can be safely eliminated
    if (use_scalapack) then
      if (n_rows .ne. mxld_soc_vec .or. n_cols .ne. mxcol_soc_vec) then
        write(info_str,'(1X,A)') '* Provided eigenvector does not appear to be &
                                 &in BLACS format.  Exiting.'
        call aims_stop( info_str, func )
      end if
      call aims_allocate(temp_complex, mxld_soc_vec, mxcol_soc_vec, &
           "temp_complex")
      call aims_allocate(moment_temp_complex, mxld_moment, mxcol_moment, &
           "moment_complex")
    else
      if (n_rows .ne. n_basis_soc.or. n_cols .ne. n_states_soc) then
        write(info_str,'(1X,A)') '* Provided eigenvector does not appear to be &
             &in full format.  Exiting.'
        call aims_stop(info_str, func)
      end if
      call aims_allocate(temp_complex, n_basis_soc, n_states_window, &
           "temp_complex" )
      call aims_allocate(moment_temp_complex, n_states_window, n_states_window,&
           "moment_complex" )
    end if

    moment = (0.0d0,0.0d0)
    temp_complex = (0.0d0,0.0d0)
    moment_temp_complex = (0.0d0,0.0d0)

    if (use_scalapack) then
      ! Calculate momentum matrix elements between eigenstates
      ! Only the spin-up basis elements coupled to one another
      call pzgemm( 'N', 'N', &
                   n_basis, n_states_window, n_basis, &
                   (1.0d0,0.0d0), &
                   mommat_full_complex, 1,1,sc_desc, &
                   KS_vec_soc, 1,n_state_min_in,sc_desc_soc_vec, &
                   (0.0d0,0.0d0), &
                   temp_complex, 1,1,sc_desc_soc_vec )
      ! Only the spin-down basis elements coupled to one another
      call pzgemm( 'N', 'N', &
                   n_basis, n_states_window, n_basis, &
                   (1.0d0,0.0d0), &
                   mommat_full_complex, 1,1,sc_desc, &
                   KS_vec_soc, n_basis+1, n_state_min_in, sc_desc_soc_vec, &
                   (1.0d0,0.0d0), &
                   temp_complex, n_basis+1,1,sc_desc_soc_vec)
      ! Now do the entire thing
      call pzgemm( 'C', 'N', &
                   n_states_window, n_states_window, n_basis_soc, &
                   (1.0d0,0.0d0), &
                   KS_vec_soc, 1,n_state_min_in,sc_desc_soc_vec, &
                   temp_complex, 1,1,sc_desc_soc_vec, &
                   (0.0d0,0.0d0), &
                   moment_temp_complex, 1,1,sc_desc_moment )

      ! Index results back into lower-triangular form expected by rest of code
      ! See comments in main subroutine for plan to eliminate this conversion
      ! from BLACS to lower-triangular
      call convert_local_to_global_low_tri_matrix_complex( &
           n_my_rows_moment, n_my_cols_moment, my_row_moment, my_col_moment, &
           mxld_moment, mxcol_moment, moment_temp_complex, &
           n_states_window, moment )
    else ! .not.use_scalapack
      ! Calculate momentum matrix elements between eigenstates
      ! Only the spin-up basis elements coupled to one another
      call zgemm( 'N', 'N', &
                  n_basis, n_states_window, n_basis, &
                  (1.0d0,0.0d0), &
                  mommat_full_complex, n_basis, &
                  KS_vec_soc(1,n_state_min_in), n_basis_soc, &
                  (0.0d0,0.0d0), &
                  temp_complex(1,1), n_basis_soc )
      ! Only the spin-down basis elements coupled to one another
      call zgemm( 'N', 'N', &
                  n_basis, n_states_window, n_basis, &
                  (1.0d0,0.0d0), &
                  mommat_full_complex, n_basis, &
                  KS_vec_soc(n_basis+1,n_state_min_in), n_basis_soc, &
                  (1.0d0,0.0d0), &
                  temp_complex(n_basis+1,1), n_basis_soc )
      ! Now do the entire thing
      call zgemm( 'C', 'N', &
                  n_states_window, n_states_window, n_basis_soc, &
                  (1.0d0,0.0d0), &
                  KS_vec_soc(1,n_state_min_in), n_basis_soc, &
                  temp_complex, n_basis_soc, &
                  (0.0d0,0.0d0), &
                  moment_temp_complex, n_states_window )

      ! Index results back into lower-triangular form expected by rest of code
      i_index = 0
      do i_state = 1, n_states_window
        do j_state = i_state, n_states_window
          i_index = i_index + 1
          moment(i_index) = moment_temp_complex(j_state,i_state)
        end do
      end do
    end if ! use_scalapack

    ! See comments in main subroutine for plan to eliminate this conversion from
    ! BLACS to lower-triangular
    if (use_scalapack) then
      call sync_vector_complex( moment, &
           ((n_states_window*(n_states_window+1))/2), my_scalapack_comm_all )
    end if

    call aims_deallocate(temp_complex, "temp_complex")
    call aims_deallocate(moment_temp_complex, "moment_temp_complex")
  end subroutine calc_moment_w0_soc
  !******

  subroutine calculate_plasma_frequency_sq(n_state_min_in, n_states_window, &
       moment_xi, moment_xj, n_states_in, n_spin_in, KS_eigenvalue, &
       chemical_potential, spin_degeneracy_in, k_weight, omega_pl_sq)
    !  PURPOSE
    !    Calculate the contribution of a particular k-point to the plasma
    !    frequency (squared) using Eq. 21 from Ambrosch-Draxl and Sofo,
    !    Comput. Phys. Commun. 175 (2006), 1-14.
    !
    !    The plasma frequency (squared) ultimately enters into the Drude-like
    !    fit for the intraband contribution to the polarizability.
    !  USES
    use constants, only: hartree, pi
    use geometry, only: cell_volume
    use dimensions, only: calculate_perturbative_soc
    use runtime_choices, only: occupation_width
    implicit none
    ! ARGUMENTS
    integer,    intent(in)    :: n_state_min_in
    integer,    intent(in)    :: n_states_window
    complex*16, intent(in)    :: &
         moment_xi((n_states_window*(n_states_window+1))/2 * n_spin_in)
    complex*16, intent(in)    :: &
         moment_xj((n_states_window*(n_states_window+1))/2 * n_spin_in)
    integer,    intent(in)    :: n_states_in
    integer,    intent(in)    :: n_spin_in
    real*8,     intent(in)    :: KS_eigenvalue(n_states_in, n_spin_in)
    real*8,     intent(in)    :: chemical_potential
    real*8,     intent(in)    :: spin_degeneracy_in
    real*8,     intent(in)    :: k_weight
    real*8,     intent(inout) :: omega_pl_sq
    !  INPUTS
    !    o n_state_min_in - Index for minimum state to include in summation
    !    o n_states_window - Total number of states (at current k-point) to
    !      include in the summation
    !    o moment_xi - i^th component of "momentum" matrix
    !    o moment_xj - j^th component of "momentum" matrix
    !    o n_states_in - Number of states in eigenvalue array (not to be
    !      confused with n_states_window!)
    !    o n_spin_in - Number of spin channels in eigenvalue array
    !    o KS_eigenvalue - Eigenvalues at current k_point
    !    o chemical_potential - Fermi level/chemical potential of electrons
    !    o spin_degeneracy_in - Degeneracy for a given spin channel
    !    o k_weight - Weight of current k_point
    !    o omegapl_sq - Current value for plasma frequency (squared)
    !  OUTPUT
    !    o omegapl_sq - New value for plasma frequency (squared)
    !  AUTHOR
    !    Tong Zhu and William Huhn (Duke University)
    !  NOTE
    !    WPH: Currently written with support only for lower-triangular moment
    !    matrices.  See main subroutine for plan to convert to BLACS support.
    !  HISTORY
    !    May 2017  - Created
    !    Sept 2017 - Updated to support SOC
    !  SOURCE
    real*8     :: width
    real*8     :: norm_gaussian
    real*8     :: prefactor
    real*8     :: gaussian
    complex*16 :: moment_product

    integer    :: m_state, n_state, i_index
    integer    :: n_state_max_in
    integer    :: i_spin

    character(*), parameter :: func = 'calculate_plasma_frequency_sq'

    ! WPH: To model the Dirac delta in Eq. 21 of Ambrosch-Draxl, we use a
    !      Gaussian function.  It is *critical* that we use the same broadening
    !      for this Gaussian as was used when determining the chemical
    !      potential/ occupation numbers.  If we use a larger broadening and
    !      our system is insulating, there is the possibility that the chemical
    !      potential will be close enough to the occupied bands to cause a
    !      spurious leakage of the Gaussian into the occupied bands, yielding an
    !      artificially non-zero value for the plasma frequency (which should be
    !      zero for an insulator on this level of theory.)
    width           = occupation_width
    norm_gaussian   = sqrt(1.0/(2.0*pi))*(1.0/width)
    prefactor       = (0.5*spin_degeneracy_in)*4.0*pi**2*(k_weight/cell_volume)
    n_state_max_in  = n_state_min_in + n_states_window - 1

    i_index = 0
    do i_spin = 1, n_spin_in
      do n_state = n_state_min_in, n_state_max_in
        gaussian = norm_gaussian * exp(-0.5*((KS_eigenvalue(n_state,i_spin) - &
                                 chemical_potential)**2/(width**2)))
        do m_state = n_state, n_state_max_in
          i_index = i_index + 1
          if (n_state .eq. m_state) then
            moment_product = moment_xi(i_index)*conjg(moment_xj(i_index))
            ! WPH:  I've added an explicit conversion to real here.  I need to
            ! go back and look at the the equations to verify if it is
            ! justified.  (This was implicit in the old code.)
            omega_pl_sq = omega_pl_sq + prefactor*gaussian*dble(moment_product)
          end if
        end do ! m_state
      end do ! n_state
    end do ! i_spin
  end subroutine calculate_plasma_frequency_sq
  !******

  subroutine calculate_plasma_frequency_sq_imagfreq( dipelementxi, &
       dipelementxj, omegapl_sq, chemical_potential, KS_eigen, k_weight, &
       widthone_in, n_state_min_in, n_state_max_in)
    !  PURPOSE
    !    Calculate the plasma frequency squared which enters into the Drude-like
    !    fit for the intraband contribution to the polarizability
    !
    !    WPH: This subroutine is identical to the real frequency case; it's been
    !    forked off as a temporary measure to allow me to rewrite the real
    !    frequency branch for ScaLAPACK and keep the imaginary frequency code
    !    as-is; after that branch is finished, we'll port the changes over to
    !    the imaginary frequency branch and unfork this subroutine
    !  USES
    use constants, only: hartree, pi
    use dimensions, only: n_states, n_spin
    use geometry, only: cell_volume
    implicit none
    !  ARGUMENTS
    integer,    intent(in)    :: n_state_min_in
    integer,    intent(in)    :: n_state_max_in
    complex*16, intent(in)    :: &
         dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                     *(n_state_max_in-n_state_min_in+1)/2)&
                     *((n_spin*(n_spin+1))/2))
    complex*16, intent(in)    :: &
         dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                      *(n_state_max_in-n_state_min_in+1)/2)&
                      *((n_spin*(n_spin+1))/2))
    real*8,     intent(inout) :: omegapl_sq
    real*8,     intent(in)    :: chemical_potential
    real*8,     intent(in)    :: KS_eigen(n_states, n_spin)
    real*8,     intent(in)    :: k_weight
    real*8,     intent(in)    :: widthone_in
    !  INPUTS
    !    o dipelementxi -- i component of momentummatrix
    !    o dipelementxj -- j component of momentummatrix
    !    o omegapl_sq -- plasma frequency**2
    !    o chemical_potential -- \epsilon_F
    !    o KS_eigen -- eigenvalues at current k_point
    !    o k_weight -- weight of current k_point
    !    o widthone_in -- broadening of Gaussian/Lorentzian
    !    o n_state_min_in -- Minimum state concidered
    !    o n_state_max_in -- Maximum state concidered
    !  OUTPUT
    !    o die_el -- dielectric function
    !    o omegapl_sq -- plasma frequency**2
    !  AUTHOR
    !    Tong Zhu and William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    real*8     :: width
    real*8     :: norm_gauss
    complex*16 :: dipmult

    integer    :: m_state, n_state, num
    integer    :: i_spin , j_spin
    integer    :: occmax

    character(*), parameter :: func = 'calculate_plasma_frequency_sq_imagfreq'

    if (n_spin.gt.1) then
      occmax = 1
    else
      occmax = 2
    end if

    width=widthone_in/hartree   ! In hartree now
    norm_gauss=sqrt(1.0/(2.0*pi))*(1.0/width)

    num=0
    do i_spin = 1, n_spin
      do n_state = n_state_min_in, n_state_max_in
        do m_state = n_state, n_state_max_in
          num = num + 1
          if (n_state == m_state) then
            dipmult = dipelementxi(num) * conjg(dipelementxj(num))
            omegapl_sq = omegapl_sq + norm_gauss*exp(-0.5*&
                 ((KS_eigen(n_state,i_spin) - chemical_potential)**2/&
                 ((width)**2)))*dipmult*(0.5*occmax)*4.0*pi**2* &
                 (k_weight/cell_volume)
          end if
        end do ! m_state
      end do ! n_state
    end do ! i_spin
  end subroutine calculate_plasma_frequency_sq_imagfreq
  !******

  subroutine calculate_interband_polarizability( n_state_min_in, &
             n_states_window, moment_xi, moment_xj, n_states_in, n_spin_in, &
             KS_eigenvalue, occ_numbers, chemical_potential, spin_degeneracy_in, &
             k_weight, broaden_method, width_in, inter_polar,inter_polar_reset )
    !  PURPOSE
    !    Calculate the contribution of a particular k-point to the interband
    !    contribution to polarizability.
    !  USES
    use constants, only: hartree, pi
    use geometry, only: cell_volume
    use runtime_choices, only: n_omega, n_omega_reset, omega_min, omega_max, &
        Emax, Emin
    implicit none
    ! ARGUMENTS
    integer,    intent(in)    :: n_state_min_in
    integer,    intent(in)    :: n_states_window
    complex*16, intent(in)    :: &
         moment_xi( ((n_states_window*(n_states_window+1))/2) * n_spin_in )
    complex*16, intent(in)    :: &
         moment_xj( ((n_states_window*(n_states_window+1))/2) * n_spin_in )
    integer,    intent(in)    :: n_states_in
    integer,    intent(in)    :: n_spin_in
    real*8,     intent(in)    :: KS_eigenvalue(n_states_in, n_spin_in)
    real*8,     intent(in)    :: occ_numbers(n_states_in, n_spin_in)
    real*8,     intent(in)    :: chemical_potential
    real*8,     intent(in)    :: spin_degeneracy_in
    real*8,     intent(in)    :: k_weight
    integer,    intent(in)    :: broaden_method
    real*8,     intent(in)    :: width_in
    complex*16, intent(inout) :: inter_polar(n_omega)
    complex*16, intent(inout) :: inter_polar_reset(n_omega_reset)
    ! INPUTS
    !   o n_state_min_in - Index for minimum state to include in summation
    !   o n_states_window - Total number of states (at current k-point) to
    !     include in the summation
    !   o moment_xi - i^th component of "momentum" matrix
    !   o moment_xj - j^th component of "momentum" matrix
    !   o n_states_in - Number of states in eigenvalue array (not to be confused
    !     with n_states_window!)
    !   o n_spin_in - Number of spin channels in eigenvalue array
    !   o KS_eigen - Eigenvalues at current k_point
    !   o occ_numbers - Occupation number at current k_point
    !   o spin_degeneracy_in - Degeneracy for a given spin channel
    !   o k_weight - Weight of current k_point
    !   o broaden_method - Choice of Lorentz-vs-Gaussian broadening
    !   o width_in  - Broadening to be used for Gaussian or Lorentz
    !   o inter_polar - Current value for interband contributions to
    !     polarizabilty
    ! OUTPUT
    !   o inter_polar - New value for interband contribution to polarizability
    ! AUTHOR
    !   Tong Zhu and William Huhn (Duke University)
    ! NOTE
    !   WPH: Currently written with support only for lower-triangular moment
    !   matrices.  See main subroutine for plan to convert to BLACS support.
    ! HISTORY
    !   May 2017  - Created
    !   Sept 2017 - Updated to support SOC
    ! SOURCE
    real*8     :: omega
    real*8     :: gauss
    real*8     :: width
    real*8     :: norm_lorentz
    real*8     :: norm_gauss
    complex*16 :: sigma
    complex*16 :: tau
    integer    :: n_state_max_in

    complex*16 :: moment_product
    real*8     :: prefactor
    real*8     :: diff_eigen
    real*8     :: diff_occ
    real*8     :: omega_step
    ! Counters
    integer    :: m_state, n_state, i_index, i_omega, i_spin
    real*8     :: diff_re_part, diff_im_part

    character(*), parameter :: func = 'calculate_interband_polarizability'

    n_state_max_in = n_state_min_in + n_states_window - 1
    width = width_in/hartree   ! Convert from eV to Hartree
    norm_lorentz = (1.0d0/(width*pi))
    norm_gauss = sqrt(1.0d0/(2.d00*pi))*(1.0d0/width)

    ! Calculate the interband contribution to the polarizability.  For Gaussian
    ! broadening, only the imaginary part will be calculated here, and the real
    ! part will be calculated later via a Kramers-Kronig transformation.  For
    ! Lorentzian broadening, the real and imaginary components will be
    ! calculated.
    if (broaden_method .eq. 1) then  ! Gaussian
      ! No spin degeneracy here, as spin degeneracy is accounted for in the
      ! occupation numbers
      prefactor = pi*4.0d0*pi*k_weight/cell_volume
      diff_re_part = 0.0d0 ! will calculate using KK transformation later
      ! as the demand for the KK transformation, we need the Im(omega) close to
      ! zero to make KK transformation converges. As aresult, for lorentzian
      ! broadening we calculate all the possible phonon ranges according to the
      ! dipole matrix window we defined before (omega_max = Emax - Emin), so
      ! that we contain all the "possible" transitions.
      ! However, as inthe output options, we still want output the omega range
      ! the user request, so we didn't change the step of the omega axis.

      omega_step = (omega_max - omega_min)/n_omega
      do i_omega = 1, n_omega_reset
        omega = (omega_min+(i_omega-1)*omega_step)/hartree


        i_index = 0
        do i_spin = 1, n_spin_in
          do n_state = n_state_min_in, n_state_max_in
            do m_state = n_state, n_state_max_in
              i_index = i_index + 1
              if (KS_eigenvalue(n_state,i_spin) .le. chemical_potential .and. &
                   KS_eigenvalue(m_state,i_spin) .ge. chemical_potential) then
                diff_eigen = KS_eigenvalue(m_state,i_spin) - &
                    KS_eigenvalue(n_state,i_spin)
                moment_product = moment_xi(i_index)* conjg(moment_xj(i_index))

                diff_occ = occ_numbers(n_state,i_spin) - occ_numbers(m_state,i_spin)
                gauss = norm_gauss*exp( -0.5d0*(diff_eigen-omega)**2/width**2 )

                ! WPH: Added real type conversion here, which was implicit
                ! before.
                diff_im_part = prefactor * dble(moment_product) * diff_occ &
                     * gauss / (diff_eigen**2)
                inter_polar_reset(i_omega) = inter_polar_reset(i_omega) &
                     + cmplx(diff_re_part, diff_im_part)
              end if
            end do ! m_state
          end do ! n_state
        end do ! i_spin
      end do ! i_omega
    else if (broaden_method .eq. 2) then ! Lorentzian
      prefactor = spin_degeneracy_in* 4.0d0 * pi * k_weight / cell_volume
      tau = cmplx(0.0d0,width)

      do i_omega = 1, n_omega
        omega = (omega_min+(i_omega-1)*((omega_max-omega_min)/n_omega))/hartree

        i_index = 0
        do i_spin = 1, n_spin_in
          do n_state = n_state_min_in, n_state_max_in
            do m_state = n_state, n_state_max_in
              i_index = i_index + 1
              if (KS_eigenvalue(n_state,i_spin) .le. chemical_potential .and. &
                   KS_eigenvalue(m_state,i_spin) .ge. chemical_potential) then
                diff_eigen = KS_eigenvalue(m_state,i_spin) &
                      - KS_eigenvalue(n_state,i_spin)
                moment_product = moment_xi(i_index)* conjg(moment_xj(i_index))

                sigma = (1.d0 / (diff_eigen + width)) * &
                          (( moment_product) /(omega-diff_eigen+tau)+&
                          conjg( moment_product)/(omega+diff_eigen+tau))
                sigma = (0.0d0, 1.0d0) * sigma

                diff_re_part = -1.0d0 * prefactor * aimag(sigma / (omega+tau))
                diff_im_part = prefactor * dble(sigma / (omega+tau))
                inter_polar(i_omega) = inter_polar(i_omega) + &
                     cmplx(diff_re_part, diff_im_part)
              end if
            end do ! m_state
          end do ! n_state
        end do ! i_spin
      end do ! i_omega
    end if ! broaden_method
  end subroutine calculate_interband_polarizability
  !******

  subroutine calculate_interband_polarizability_imagfreq(dipelementxi, &
       dipelementxj, inter_polar, chemical_potential, KS_eigen, k_weight, &
       broaden_method, widthone_in, n_state_min_in, n_state_max_in)
    !  PURPOSE
    !    Sum up momentum matrix elements to generate interband contributions
    !    to polarizability.
    !  USES
    use constants, only: hartree, pi
    use dimensions, only: n_states, n_spin
    use geometry, only: cell_volume
    use gw_para, only: n_full_freq, freq_grid_type, n_freq, omega_grid, &
        omegamax, womega, womega_full, omega_full_grid, allocate_gw, &
        deallocate_gw
    use mpi_tasks, only: aims_stop, myid
    use runtime_choices, only: n_omega, omega_min, omega_max
    implicit none
    !  ARGUMENTS
    integer,    intent(in)    :: n_state_min_in
    integer,    intent(in)    :: n_state_max_in
    complex*16, intent(in)    :: &
         dipelementxi((((n_state_max_in-n_state_min_in+1)+1)&
                     *(n_state_max_in-n_state_min_in+1)/2)&
                     *((n_spin*(n_spin+1))/2))
    complex*16, intent(in)    :: &
         dipelementxj ((((n_state_max_in-n_state_min_in+1)+1)&
                      *(n_state_max_in-n_state_min_in+1)/2)&
                      *((n_spin*(n_spin+1))/2))
    complex*16, intent(inout) :: inter_polar(n_full_freq)
    real*8,     intent(in)    :: chemical_potential
    real*8,     intent(in)    :: KS_eigen(n_states, n_spin)
    real*8,     intent(in)    :: k_weight
    real*8,     intent(in)    :: widthone_in
    integer,    intent(in)    :: broaden_method
    !  INPUTS
    !    o dipelementxi -- i component of momentummatrix
    !    o dipelementxj -- j component of momentummatrix
    !    o die_el -- dielectric function
    !    o chemical_potential -- \epsilon_F
    !    o KS_eigen -- eigenvalues at current k_point
    !    o k_weight -- weight of current k_point
    !    o widthone_in -- broadening of Gaussian/Lorentzian
    !    o broaden_method -- broadening type: 1 for Gaussian; 2 for Lorentzian
    !    o n_state_min_in -- Minimum state concidered
    !    o n_state_max_in -- Maximum state concidered
    !  OUTPUT
    !    o die_el -- dielectric function
    !  AUTHOR
    !    Tong Zhu and William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    real*8     :: omega
    real*8     :: fermione
    real*8     :: fermitwo
    real*8     :: scaling
    real*8     :: width
    real*8     :: width_onset
    real*8     :: norm_gauss_onset
    real*8     :: norm_lorentz
    real*8     :: norm_gauss
    complex*16 :: sigma
    complex*16 :: tau

    complex*16 :: dipmult
    real*8     :: ohm
    real*8     :: dfermi

    integer    :: n_state
    integer    :: m_state
    integer    :: num
    integer    :: omegaind
    integer    :: i_spin , j_spin
    integer    :: occmax

    real*8 :: diff_re_part, diff_im_part

    character(*), parameter :: func = &
         'calculate_interband_polarizability_imagfreq'

    if (n_spin.gt.1) then
      occmax = 1
    else
      occmax = 2
    end if

    width = widthone_in/hartree   ! In hartree now
    width_onset = 0.00001/hartree ! small broadening restriced for metals like
                                  ! materials in the first transition matrix
                                  ! component
    norm_lorentz = (1.0/(width*pi))
    norm_gauss = sqrt(1.0/(2.0*pi))*(1.0/width)
    norm_gauss_onset = sqrt(1.0/(2.0*pi))*(1.0/width_onset)

    call deallocate_gw() ! NEW: for atom_bsse
    call allocate_gw()
    !print * ,'freq_grid_type:' ,freq_grid_type
    if (freq_grid_type.eq.0) then
      call tf_ini(n_freq,n_full_freq, omegamax,omegamax, &
           omega_grid,omega_full_grid, &
           womega,womega_full,.true.)
    else if (freq_grid_type.eq.1) then
      n_freq=n_full_freq
      call tf_ini_trans(n_freq,n_full_freq, omegamax,omegamax, &
           omega_grid,omega_full_grid, &
           womega,womega_full,.true.)
    else
      call aims_stop("Unsupported freq_grid_type", func)
    end if

    ! Calculate the interband contribution to the polarizability.  For Gaussian
    ! broadening, only the imaginary part will be calculated here, and the real
    ! part will be calculated later via a Kramers-Kronig transformation.  For
    ! Lorentzian broadening, the real and imaginary components will be
    ! calculated.
    if (broaden_method == 1) then  ! Gaussian
      diff_re_part = 0.0d0 ! will calculate using KK transformation later
      do omegaind = 0, n_omega-1
        omega=(omega_min+omegaind*((omega_max-omega_min)/n_omega))/hartree

        num=0
        do i_spin = 1, n_spin
          do j_spin = i_spin, n_spin
            do n_state = n_state_min_in, n_state_max_in
              fermitwo=(1.0/(1.0 + exp((KS_eigen(n_state,i_spin) - &
                   chemical_potential)/(0.01/hartree))))
              do m_state = n_state, n_state_max_in
                num = num + 1
                if (KS_eigen(n_state,j_spin) <= chemical_potential .and. &
                     KS_eigen(m_state,i_spin) >=chemical_potential) then
                  fermione = (1.0/(1.0 + exp((KS_eigen(m_state,1) - &
                       chemical_potential)/(0.01/hartree))))
                  dfermi = fermitwo - fermione
                  ohm = KS_eigen(m_state,j_spin)-KS_eigen(n_state,i_spin)
                  dipmult = dipelementxi(num)* conjg(dipelementxj(num))
                  scaling=norm_gauss*exp(-0.5*(((ohm-omega)**2)/((width)**2)))

                  if ( ohm <= 0.2/hartree.and. omega < ohm ) then
                    scaling = norm_gauss_onset*exp(-0.5*&
                         (((ohm-omega)**2)/((width_onset)**2)))
                  end if
                  diff_im_part = (occmax*pi*4.0*pi*k_weight/cell_volume)*&
                       dipmult*dfermi*scaling/((omega**2))
                  inter_polar(omegaind+1) = inter_polar(omegaind+1) + &
                       cmplx(diff_re_part, diff_im_part)
                end if
              end do ! m_state
            end do ! n_state
          end do ! j_spin
        end do ! i_spin
      end do ! omegaind
      if (omega_min == 0.0d0) then
        inter_polar(1) = (0.0d0,0.0d0)
      end if
    else if (broaden_method == 2) then ! Lorentzian
      tau = dcmplx(0d0,width)
      do omegaind = 0, n_full_freq - 1
        omega=omega_full_grid(omegaind+1)
        num=0
        do i_spin = 1, n_spin
          do n_state = n_state_min_in, n_state_max_in
            do m_state = n_state, n_state_max_in
              num = num + 1
              if (KS_eigen(n_state,i_spin) <= chemical_potential .and. &
                  KS_eigen(m_state,i_spin) >= chemical_potential) then

                ohm = KS_eigen(m_state,i_spin) - KS_eigen(n_state,i_spin)
                dipmult = dipelementxi(num)* conjg(dipelementxj(num))
                sigma = 1.d0/ohm**2 * occmax* dipmult* 2.d0 * &
                     ohm/(omega*omega+ohm*ohm)
                sigma =  (k_weight / cell_volume) * sigma

                diff_re_part = (4.0*pi)*dble(sigma)
                diff_im_part = 4.0*pi*aimag(sigma)
                inter_polar(omegaind+1) = inter_polar(omegaind+1) &
                     + cmplx(diff_re_part, diff_im_part)
    !            if(myid.eq.0 .and. omegaind .eq. 0) then
    !               write(6,'(2I4,6f18.8)') m_state, n_state, ohm, dipmult, inter_polar(omegaind+1)
    !            endif
              end if
            end do ! m_state
          end do ! n_state
        end do ! i_spin
      end do ! omegaind
    end if ! broaden_method
  end subroutine calculate_interband_polarizability_imagfreq
  !******

  subroutine kramerskronig_im_to_re(inter_polar)
    !  PURPOSE
    !   Do Kramer-Kroning transformation to obtain real part of the interband
    !   contribution to the polarizability via the imaginary part
    !  USES
    use constants, only: hartree, pi
    use runtime_choices, only: omega_min, omega_max, n_omega, n_omega_reset
    implicit none
    !  ARGUMENTS
    complex*16, intent(inout) :: inter_polar(n_omega_reset)
    !  INPUTS
    !    o inter_polar (imaginary part) - interband contribution to
    !      polarizability
    !  OUTPUT
    !    o inter_polar (real part) - interband contribution to polarizability
    !  AUTHOR
    !    Tong Zhu and William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  TODO
    !    Use better quadrature.
    !    Allow the re -> im transform via flag (and rename)
    !  SOURCE
    real*8  :: omega(n_omega_reset)
    real*8  :: omega_step
    integer :: omegaind
    integer :: ip
    real*8  :: myomega
    real*8  :: myomegap
    real*8  :: acc
    real*8  :: re_inter_polar, im_inter_polar, im_inter_polarp

    character(*), parameter :: func = 'kramerskronig_im_to_re'

    !set the frequency grid
    omega = 0.0
    omega_step = (omega_max - omega_min)/n_omega
    do omegaind = 1, n_omega_reset
      omega(omegaind) =(omega_min+(omegaind-1)*omega_step)/hartree
    end do
    omega_step = omega_step/hartree

    ! Perform KK transform
    do omegaind = 1, n_omega_reset
      myomega = omega(omegaind)
      im_inter_polar = aimag(inter_polar(omegaind))
      acc = 0.0
      do ip = 1, n_omega_reset
        if (ip == omegaind) cycle
        im_inter_polarp = aimag(inter_polar(ip))
        myomegap = omega(ip)
        acc = acc + myomegap/(myomegap**2 - myomega**2)*im_inter_polarp
      end do
      re_inter_polar = 2.0/pi*acc*omega_step
      inter_polar(omegaind) = cmplx(re_inter_polar, im_inter_polar)
    end do
  end subroutine kramerskronig_im_to_re
  !******

  subroutine calculate_drude_like_intraband_polarizability &
       (omega_min, omega_max, n_omega, omega_pl_sq, widthone_in, &
        intra_polar)
    !  PURPOSE
    !   Calculate intraband contribution to a tensor component of the
    !   macroscopic polarizability based on a Drude-like model fitted
    !   to the calculated plasma frequency
    !   (See Ambrosch-Draxl and Sofo, Comp Phys Comm 175, 1-14 (2006))
    !   Note we do not add the "1" (actually Kronecker delta) to the real
    !   part like the reference did, as we are calculating the polarizability
    !   here (i.e. the dielectric function minus the identity tensor)
    !  USES
    use constants, only : hartree
    implicit none
    !  ARGUMENTS
    real*8,     intent(in)  :: omega_min
    real*8,     intent(in)  :: omega_max
    integer,    intent(in)  :: n_omega
    real*8,     intent(in)  :: omega_pl_sq
    real*8,     intent(in)  :: widthone_in
    complex*16, intent(out) :: intra_polar(n_omega)
    !  AUTHOR
    !    Tong Zhu and William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer :: i_omega
    real*8  :: diff_omega, omega
    real*8  :: broadening
    real*8  :: re_intra_polar, im_intra_polar

    character(*), parameter :: func = &
         'calculate_drude_like_intraband_polarizability'

    broadening=widthone_in/hartree
    diff_omega = (omega_max-omega_min)/n_omega
    do i_omega = 1, n_omega
      omega = (omega_min + (i_omega-1) * diff_omega) / hartree
      re_intra_polar = -1.0d0 * omega_pl_sq / (omega**2 + broadening**2)
      im_intra_polar =  &
           broadening * omega_pl_sq / ( omega* (omega**2 + broadening**2) )
      intra_polar = cmplx(re_intra_polar,im_intra_polar)
    end do
    ! Why do we do this?
    if (omega_min == 0.0d0) then
      im_intra_polar = cmplx(real(intra_polar(1)),0.0d0)
    end if
  end subroutine calculate_drude_like_intraband_polarizability
  !******

  subroutine calculate_dielectric_from_polarizability &
       (n_omega, broaden_method, component_is_diagonal, &
        inter_polar, intra_polar, eps )
    !  PURPOSE
    !   Calculate the tensor component of the macroscopic dielectric tensor
    !   from the interband and intraband contributions to the macroscopic
    !   polarizability.  That is, sum up the contributions and then add the
    !   vacuum contribution (i.e. the Kronecker delta)
    !  USES
    use constants, only: hartree
    implicit none
    !  ARGUMENTS
    integer,    intent(in)  :: n_omega
    integer,    intent(in)  :: broaden_method
    logical,    intent(in)  :: component_is_diagonal
    complex*16, intent(in)  :: inter_polar(n_omega)
    complex*16, intent(in)  :: intra_polar(n_omega)
    complex*16, intent(out) :: eps(n_omega)
    !  AUTHOR
    !    William Huhn and Tong Zhu (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer    :: i_omega
    complex*16 :: vacuum_contribution

    character(*), parameter :: func = 'calculate_dielectric_from_polarizability'

    if (component_is_diagonal) then
      vacuum_contribution = (1.0d0,0.0d0)
    else
      vacuum_contribution = (0.0d0,0.0d0)
    end if

    eps = vacuum_contribution + intra_polar + inter_polar
  end subroutine calculate_dielectric_from_polarizability
  !******

  subroutine calculate_absorption_coeff &
       (omega_min, omega_max, n_omega, eps, absorption)
    !  PURPOSE
    !    Calculate the components of the absorption coefficient from the
    !    macroscopic dielectric tensor.  Note this formula assumes that the
    !    dielectric tensor is diagonal.
    !    (See Ambrosch-Draxal and Sofo, Comp Phys Comm 175, 1-14 (2006) or any
    !    graduate-level EM textbook)
    !  USES
    use constants, only: hartree, pi
    implicit none
    !  ARGUMENTS
    real*8,     intent(in)  :: omega_min
    real*8,     intent(in)  :: omega_max
    integer,    intent(in)  :: n_omega
    complex*16, intent(in)  :: eps(n_omega)
    real*8,     intent(out) :: absorption(n_omega)
    !  AUTHOR
    !    William Huhn and Tong Zhu (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer :: i_omega
    real*8  :: diff_omega, omega, magnitude
    real*8  :: evtocmm1

    character(*), parameter :: func = 'calculate_absorption_coeff'

    evtocmm1= 10000.0d0/1.23981d0  !Unit conversion: ev to cm-1

    diff_omega = (omega_max-omega_min)/n_omega
    do i_omega = 1, n_omega
      omega = omega_min + (i_omega-1) * diff_omega
      magnitude = abs(eps(i_omega))

      absorption(i_omega) = evtocmm1* 4.0d0 * pi * &
           omega * sqrt( (magnitude-real(eps(i_omega)))/2.0d0)
    end do
  end subroutine calculate_absorption_coeff
  !******

  subroutine output_optical_properties (i_direction,component_is_diagonal, &
       broaden_method, broaden_width, inter_polar, intra_polar, eps, &
       absorption, omegapl_sq, direction1, direction2 )
    ! PURPOSE
    !   Writes dielectric function \epsilon_ep1_ep2 to file.
    ! USES
    use constants, only: hartree
    use dimensions, only: calculate_perturbative_soc
    use localorb_io, only: localorb_info
    use runtime_choices, only: n_omega, omega_min, omega_max
    implicit none
    !  ARGUMENTS
    character*120          :: info_str
    integer,    intent(in) :: i_direction
    logical,    intent(in) :: component_is_diagonal
    complex*16, intent(in) :: inter_polar(n_omega)
    complex*16, intent(in) :: intra_polar(n_omega)
    complex*16, intent(in) :: eps(n_omega)
    real*8,     intent(in) :: absorption(n_omega)
    real*8,     intent(in) :: omegapl_sq
    integer,    intent(in) :: broaden_method
    real*8,     intent(in) :: broaden_width
    integer,    intent(in) :: direction1
    integer,    intent(in) :: direction2
    !  AUTHOR
    !    William Huhn and Tong Zhu (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer           :: i_omega
    real*8            :: omega, diff_omega
    CHARACTER(len=100) :: fmt
    CHARACTER(len=100) :: name
    character*30      :: ep1_in
    character*30      :: ep2_in
    character*30      :: broaden_method_name
    character*30      :: width_name

    character(*), parameter :: func = 'output_optical_properties'

    diff_omega = (omega_max-omega_min)/n_omega

    if (direction1 == 1) then
      ep1_in  = "x"
    else if (direction1 == 2) then
      ep1_in = "y"
    else if (direction1 ==3) then
      ep1_in = "z"
    end if

    if (direction2 == 1) then
      ep2_in  = "x"
    else if (direction2 == 2) then
      ep2_in = "y"
    else if (direction2 ==3) then
      ep2_in = "z"
    end if

    if (broaden_method .eq. 1) then
      broaden_method_name = "Gaussian"
    else if (broaden_method .eq. 2) then
      broaden_method_name = "Lorentzian"
    else
      broaden_method_name = "Unknown" ! Not worthy exiting aims over.
    end if

    ! write output information in FHI-aims output
    write(info_str, *)
    call localorb_info ( info_str )
    write(info_str,'(2X,A,1X,A,A,A)') &
         "Summary for dielectric component # ", trim(ep1_in), '_', trim(ep2_in)
    call localorb_info ( info_str )
    write(info_str, '(2X,A,1X,ES17.5,A,1x,A,A,A)') &
         "| Broadening parameters    :", broaden_width, " eV", "(", &
         trim(broaden_method_name), ")"
    call localorb_info ( info_str )
    write(info_str,'(2X,A,1X,ES17.5,A)') &
         "| Plasma frequency         :",  sqrt(abs(omegapl_sq))*hartree, " eV"
    call localorb_info ( info_str )
    write(info_str,'(2X,A,1X,ES17.5,A)') &
         "| Minimum omega value      :",  omega_min, " eV"
    call localorb_info ( info_str )
    write(info_str,'(2X,A,1X,ES17.5)') &
         "| Epsilon at minimum omega :",  real(eps(1))
    call localorb_info ( info_str )
    write(info_str, *)
    call localorb_info ( info_str )

    if (calculate_perturbative_soc) then
      write (width_name,'(f8.4)') broaden_width
      write (unit=name,fmt='(A,A,A,A,A,A,A,A,A)') &
           'dielectric_function_soc_', trim(broaden_method_name), '_', &
           trim(adjustl(width_name)), '_', trim(ep1_in), '_', trim(ep2_in), &
           '.out'
    else
      write (width_name,'(f8.4)') broaden_width
      write (unit=name,fmt='(A,A,A,A,A,A,A,A,A)') &
           'dielectric_function_', trim(broaden_method_name), '_', &
            trim(adjustl(width_name)), '_', trim(ep1_in), '_', trim(ep2_in), &
            '.out'
    end if
    name = trim(name)

    ! Write header for output dielectric file
    open(unit=8, file=name,ACTION='WRITE')
    write (8,'(A,A,A,ES17.4,A)') "# ", trim(broaden_method_name), &
         " broadening used, width ", broaden_width, " eV"
    write(8,423) sqrt(abs(omegapl_sq))*hartree
    423 FORMAT ('# Calculated plasma frequency [eV]: ',ES17.4)
    write(8,522) trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in), &
                 trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in), &
                 trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in), &
                 trim(ep1_in), trim(ep2_in),trim(ep1_in), trim(ep2_in)
    522 FORMAT ('# omega [eV],  &
                &Re(Vacuum)_{',A1,A1,'}, &
                &Im(Vacuum)_{',A1,A1,'}, &
                &Re(Polarizability)^{inter}_{',A1,A1,'}, &
                &Im(Polarizability)^{inter}_{',A1,A1,'}, &
                &Re(Polarizability)^{intra}_{',A1,A1,'}, &
                &Im(Polarizability)^{intra}_{',A1,A1,'}, &
                &Re(Dielectric)_{',A1,A1,'}, &
                &Im(Dielectric)_{',A1,A1,'}')

    fmt = '(9ES19.6E3)'
    ! Write values for dielectric component at each omega to output file
    do i_omega = 1, n_omega
      omega = omega_min + (i_omega-1) * diff_omega
      if (direction1 .eq. direction2) then
        write(8,fmt) omega, 1.0d0, 0.0d0, real(inter_polar(i_omega)), &
             aimag(inter_polar(i_omega)), real(intra_polar(i_omega)), &
             aimag(intra_polar(i_omega)), real(eps(i_omega)), &
             aimag(eps(i_omega))
      else
        write(8,fmt) omega, 0.0d0, 0.0d0, real(inter_polar(i_omega)), &
             aimag(inter_polar(i_omega)), real(intra_polar(i_omega)), &
             aimag(intra_polar(i_omega)), real(eps(i_omega)), &
             aimag(eps(i_omega))
      end if
    end do
    close(unit=8)

    ! Write header for output absorption file
    if (component_is_diagonal) then
      if (calculate_perturbative_soc) then
        write (width_name,'(f8.4)') broaden_width
        write (unit=name,fmt='(A,A,A,A,A,A,A,A,A)') &
             'absorption_soc_', trim(broaden_method_name), '_', &
             trim(adjustl(width_name)), '_', trim(ep1_in), '_', trim(ep2_in), &
             '.out'
      else
        write (width_name,'(f8.4)') broaden_width
        write (unit=name,fmt='(A,A,A,A,A,A,A,A,A)') &
             'absorption_',trim(broaden_method_name), '_', &
             trim(adjustl(width_name)), '_', trim(ep1_in), '_', trim(ep2_in), &
             '.out'
      end if
      fmt = '(9ES17.4E3)'
      open(unit=8, file=name,ACTION='WRITE')
      write (8,'(A,A,A,ES17.4,A)') "# ", trim(broaden_method_name), &
           " broadening used, width ", broaden_width, " eV"
      write(8,122) trim(ep1_in)
      122 FORMAT ('# \omega [eV], \alpha_{',A1,'}')

      ! Write values for absorption coefficient at each omega to output file
      do i_omega = 1,n_omega
        omega = omega_min + (i_omega-1) * diff_omega
        if (abs(absorption(i_omega))>1.0E-30) then
          write(8,fmt) omega, absorption(i_omega)
        else
          write(8,fmt) omega, 0.0
        end if
      end do
      close(unit=8)
    end if
  end subroutine output_optical_properties
  !******

  subroutine simpson_int(simp_delta,simp_input,simp_res)
    !  PURPOSE
    !    Perform Simpson integration.  Currently not used.
    !  USES
    use runtime_choices, only: n_omega
    implicit none
    !  ARGUMENTS
    real*8, intent(in)  :: simp_delta
    real*8, intent(in)  :: simp_input(n_omega)
    real*8, intent(out) :: simp_res(n_omega)
    !  AUTHOR
    !    Tong Zhu (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer :: ii
    real*8 :: coef1, coef2, coef3

    character(*), parameter :: func = 'simpson_int'

    coef1 = 9.0/ 24.0
    coef2 = 28.0/ 24.0
    coef3 = 23.0/ 24.0

    simp_res(1) =  coef1*simp_input(1)
    simp_res(2) = simp_res(1) + coef2*simp_input(2)
    simp_res(3) = simp_res(2) + coef3*simp_input(3)

    do ii = 4, n_omega-3
      simp_res(ii) = simp_res(ii-1)+simp_input(ii)
    end do

    simp_res(n_omega-2) = simp_res(n_omega-3) + coef3*simp_input(n_omega-2)
    simp_res(n_omega-1) = simp_res(n_omega-2) + coef2*simp_input(n_omega-1)
    simp_res(n_omega) = simp_res(n_omega-1) + coef1*simp_input(n_omega)

    simp_res = simp_res*simp_delta
  end subroutine simpson_int
  !******

  !!!!!!!!!!!!!!!!!!!!!
  ! MEMORY MANAGEMENT !
  !!!!!!!!!!!!!!!!!!!!!

  subroutine allocate_mommat_TZ( ld_mommat )
    !  PURPOSE
    !    allocation of momentum matrix
    !  USES
    use aims_memory_tracking, only: aims_allocate
    use dimensions, only: n_hamiltonian_matrix_size
    use load_balancing, only: batch_perm, use_batch_permutation, n_bp_integ
    use runtime_choices, only : calc_px_dielectric, calc_py_dielectric, &
        calc_pz_dielectric, use_load_balancing
    implicit none
    !  ARGUMENTS
    integer, intent(out) :: ld_mommat
    !  INPUTS
    !    None
    !  OUTPUT
    !    o ld_mommat - The sizes of the matrix allocated
    !  AUTHOR
    !    Tong Zhu (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE

    character(*), parameter :: func = 'allocate_mommat_TZ'

    if (use_batch_permutation .eq. n_bp_integ) then
      ld_mommat = batch_perm(n_bp_integ)%n_local_matrix_size
    else
      ld_mommat = n_hamiltonian_matrix_size
    end if

    if (calc_px_dielectric) then
      if (.not.allocated(mommat_full_x_up)) then
        call aims_allocate(mommat_full_x_up, ld_mommat, "+mommat_full_x_up")
      end if
      if (.not.allocated(mommat_full_x_low)) then
        call aims_allocate(mommat_full_x_low, ld_mommat, "+mommat_full_x_low")
      end if
    end if
    if (calc_py_dielectric) then
      if (.not.allocated(mommat_full_y_up)) then
        call aims_allocate(mommat_full_y_up, ld_mommat, "+mommat_full_y_up")
      end if
      if (.not.allocated(mommat_full_y_low)) then
        call aims_allocate(mommat_full_y_low, ld_mommat, "+mommat_full_y_low")
      end if
    end if
    if (calc_pz_dielectric) then
      if (.not.allocated(mommat_full_z_up)) then
        call aims_allocate(mommat_full_z_up, ld_mommat, "+mommat_full_z_up")
      end if
      if (.not.allocated(mommat_full_z_low)) then
        call aims_allocate(mommat_full_z_low, ld_mommat, "+mommat_full_z_low")
      end if
    end if
  end subroutine allocate_mommat_TZ
  !******

  subroutine deallocate_mommat_TZ
    !  PURPOSE
    !    deallocation of momentum matrix
    !  USES
    use aims_memory_tracking, only: aims_deallocate
    implicit none
    !  ARGUMENTS
    !    None
    !  INPUTS
    !    None
    !  OUTPUT
    !    None
    !  AUTHOR
    !    William Huhn (Duke University)
    !  HISTORY
    !    Created in September 2017
    !  SOURCE

    character(*), parameter :: func = 'deallocate_mommat_TZ'

    if (allocated(mommat_full_x_up )) &
      call aims_deallocate(mommat_full_x_up,"mommat_full_x_up")
    if (allocated(mommat_full_x_low)) &
      call aims_deallocate(mommat_full_x_low, "mommat_full_x_low")

    if (allocated(mommat_full_y_up )) &
      call aims_deallocate(mommat_full_y_up, "mommat_full_y_up")
    if (allocated(mommat_full_y_low)) &
      call aims_deallocate(mommat_full_y_low,"mommat_full_y_low")

    if (allocated(mommat_full_z_up )) &
      call aims_deallocate(mommat_full_z_up, "mommat_full_z_up")
    if (allocated(mommat_full_z_low)) &
      call aims_deallocate(mommat_full_z_low, "mommat_full_z_low")
  end subroutine deallocate_mommat_TZ
  !******

  subroutine init_KS_optical_properties
    !  PURPOSE
    !    Allocate module variables for KS_optical_properties
    !  USES
    use aims_memory_tracking, only: aims_allocate
    use dimensions, only: n_plot_dielectric
    use runtime_choices, only: n_omega
    implicit none
    !  ARGUMENTS
    !    None
    !  INPUTS
    !    None
    !  OUTPUT
    !    None
    !  AUTHOR
    !    William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer :: info

    character(*), parameter :: func = 'init_KS_optical_properties'

    if (.not.allocated(omega_plasma_sq)) then
      call aims_allocate(omega_plasma_sq, n_plot_dielectric, "+omega_plasma_sq")
    end if

    if (.not.allocated(interband_polarizability)) then
      call aims_allocate( interband_polarizability, n_omega,n_plot_dielectric, &
           "+interband_polarizability")
    end if

    if (.not.allocated(intraband_polarizability)) then
      call aims_allocate(intraband_polarizability, n_omega,n_plot_dielectric, &
           "+intraband_polarizability")
    end if

    if (.not.allocated(dielectric)) then
      call aims_allocate(dielectric, n_omega,n_plot_dielectric, "+dielectric" )
    end if

    if (.not.allocated(absorption_coeff)) then
      call aims_allocate(absorption_coeff, n_omega,n_plot_dielectric, &
           "+absorption_coeff")
    end if
  end subroutine init_KS_optical_properties
  !******

  subroutine cleanup_KS_optical_properties
    !  PURPOSE
    !    Deallocate module variables for KS_optical_properties
    !  USES
    use aims_memory_tracking, only: aims_deallocate
    implicit none
    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUT
    !    none
    !  AUTHOR
    !    William Huhn (Duke University)
    !  HISTORY
    !    Created in May 2017
    !  SOURCE
    integer :: info

    character(*), parameter :: func = 'cleanup_KS_optical_properties'

    if (allocated(omega_plasma_sq)) &
         call aims_deallocate(omega_plasma_sq, "omega_plasma_sq")
    if (allocated(interband_polarizability)) &
         call aims_deallocate(interband_polarizability, &
                              "interband_polarizability" )
    if (allocated(intraband_polarizability)) &
         call aims_deallocate(intraband_polarizability, &
                              "intraband_polariazbility" )
    if (allocated(dielectric)) call aims_deallocate( dielectric, "dielectric" )
    if (allocated(absorption_coeff)) call aims_deallocate( absorption_coeff, &
         "absorption_coeff" )
  end subroutine cleanup_KS_optical_properties
  !******

end module KS_optical_properties
!******
