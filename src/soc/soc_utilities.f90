!****h* FHI-aims/soc_utilities
!*  NAME
!*    soc_utilities
!*  SYNOPSIS
module soc_utilities
!*  PURPOSE
!*    This module contains various subroutines used by the SOC implementation in
!*    FHI-aims which have been developed specifically for SOC, i.e., have no
!*    analogues in the main FHI-aims code base.  Forks of core FHI-aims
!*    functionality (such as batch integration for SOC, construction of the
!*    Hamiltonian for SOC, etc.) should remain in their own source files for
!*    consistency with main code body.
!*  USES
  implicit none
!*  AUTHOR
!*    William Huhn (Duke University)
!*  NOTES
!*    The implementation of second-variational SOC in FHI-aims is published in
!*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
!*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
!*  HISTORY
!*    July 2017 - Created from combining subroutines in scattered source files.
!*  COPYRIGHT
!*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!*    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!*    the terms and conditions of the respective license agreement."
!*  SOURCE

  private

  ! Stored SR values, used when changing between SR/SOC environments
  ! We store these as module variables instead of arguments to subroutines
  ! to avoid changing the interface if we ever need to modify additional
  ! module variables
  integer :: n_basis_sr_save
  integer :: n_states_sr_save
  integer :: n_spin_sr_save
  real*8 :: n_electrons_sr_save
  real*8 :: spin_degeneracy_sr_save
  logical :: real_eigenvectors_sr_save

  public ::  convert_sr_to_soc_environment
  public ::  revert_soc_to_sr_environment
  public ::  convert_wf_basis_to_compute_basis
  public ::  calculate_spin_expectation
  public ::  find_min_energy_from_gap
  public ::  find_num_core_states_from_energy
  public ::  find_num_high_states_from_energy
  public ::  create_sorted_sr_eigenvalues_list
  public ::  create_sorted_sr_eigenvalues_idxmap
  public ::  perform_soc_perturbation
  public ::  write_soc_values
  public ::  write_soc_perturbed_eigenvectors

contains

  !****f* soc_utilities/convert_sr_to_soc_environment
  !*  NAME
  !*    convert_sr_to_soc_environment
  !*  SYNOPSIS
  subroutine convert_sr_to_soc_environment( )
  !*  PURPOSE
  !*    Swap out various indexing parameters from the original
  !*    (scalar-relativistic) values to the spin-orbit coupled values,
  !*    converting the working environment to a spin-orbit-coupled environment.
  !*    The original values are stored later by revert_soc_to_sr_environment().
  !*    This allows many subroutines written for scalar-relativistic quantities
  !*    to be reused by spin-orbit-coupled quantities, as the underlying math is
  !*    often the same and only the indices and spin-related quantities differ.
  !*  USES
    use dimensions, only: n_basis, n_states, n_spin, spin_degeneracy
    use physics, only: n_electrons
    use runtime_choices, only: real_eigenvectors
    use dimensions_soc, only: n_basis_soc, n_saved_states_soc, &
        soc_saved_state_start
    implicit none
  !*  ARGUMENTS
  !*    o None (modifies module variables)
  !*  INPUTS
  !*    o None (modifies module variables)
  !*  OUTPUTS
  !*    o None (modifies module variables)
  !*  AUTHORS
  !*    William Huhn
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    n_basis_sr_save           = n_basis
    n_states_sr_save          = n_states
    n_spin_sr_save            = n_spin
    n_electrons_sr_save       = n_electrons
    spin_degeneracy_sr_save   = spin_degeneracy
    real_eigenvectors_sr_save = real_eigenvectors

    n_basis                   = n_basis_soc
    n_states                  = n_saved_states_soc
    n_spin                    = 1
    n_electrons               = n_electrons - soc_saved_state_start + 1
    spin_degeneracy           = 1.0d0
    real_eigenvectors         = .false.
  end subroutine convert_sr_to_soc_environment
  !******

  !****f* soc_utilities/revert_soc_to_sr_environment
  !*  NAME
  !*    revert_soc_to_sr_environment
  !*  SYNOPSIS
  subroutine revert_soc_to_sr_environment( )
  !*  PURPOSE
  !*    Reverse operation of convert_sr_to_soc_environment():  restore the
  !*    stored original/scalar-relativistic values for indexing parameters,
  !*    reverting back to the original scalar-relativistic environment.
  !*  USES
    use dimensions, only: n_basis, n_states, n_spin, spin_degeneracy
    use physics, only: n_electrons
    use runtime_choices, only: real_eigenvectors
    implicit none
  !*  ARGUMENTS
  !*    o None (modifies module variables)
  !*  INPUTS
  !*    o None (modifies module variables)
  !*  OUTPUTS
  !*    o None (modifies module variables)
  !*  AUTHORS
  !*    William Huhn
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE
    n_basis                 = n_basis_sr_save
    n_states                = n_states_sr_save
    n_spin                  = n_spin_sr_save
    n_electrons             = n_electrons_sr_save
    spin_degeneracy         = spin_degeneracy_sr_save
    real_eigenvectors       = real_eigenvectors_sr_save
  end subroutine revert_soc_to_sr_environment
  !******

  !****f* soc_utilities/convert_wf_basis_to_compute_basis
  !*  NAME
  !*    convert_wf_basis_to_compute_basis
  !*  SYNOPSIS
  subroutine convert_wf_basis_to_compute_basis( &
       n_eigenvec_rows, n_eigenvec_cols, eigenvec, eigenvec_complex, &
       n_wf_bas_rows, n_wf_bas_cols, eigenvec_soc_wf_basis, &
       n_eigenvec_soc_rows, n_eigenvec_soc_cols, eigenvec_soc )
  !*  PURPOSE
  !*    The result of the second-variational step will be eigenvectors in the
  !*    variational basis, e.g. as a linear sum of the unperturbed KS
  !*    eigenvectors.  This functional converts them into the computational
  !*    basis, e.g. as a linear sum of computational basis elements (NAOs,
  !*    Gaussians, etc.) at a given k-point.
  !*  USES
    use dimensions, only: n_basis, n_states, n_spin
    use runtime_choices, only: real_eigenvectors, use_scalapack
    use mpi_tasks, only: aims_stop
    use dimensions_soc, only: n_states_sr, sr_state_start, n_states_soc,&
        n_basis_soc, n_basis_soc_coll, n_basis_soc_ncoll, n_saved_states_soc, &
        soc_saved_state_start, n_core_states_omit_from_soc
    use synchronize_mpi_basic, only: sync_vector_complex
    use scalapack_wrapper, only: sc_desc, my_scalapack_comm_all
    use scalapack_soc
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: n_eigenvec_rows
    integer, intent(in) :: n_eigenvec_cols
    real*8, intent(in) :: eigenvec(n_eigenvec_rows, n_eigenvec_cols, n_spin)
    complex*16, intent(in) :: eigenvec_complex(n_eigenvec_rows, n_eigenvec_cols, n_spin)
    integer, intent(in) :: n_wf_bas_rows
    integer, intent(in) :: n_wf_bas_cols
    complex*16, intent(in) :: eigenvec_soc_wf_basis(n_wf_bas_rows,n_wf_bas_cols)
    integer, intent(in) :: n_eigenvec_soc_rows
    integer, intent(in) :: n_eigenvec_soc_cols
    complex*16, intent(out) :: eigenvec_soc(n_eigenvec_soc_rows, n_eigenvec_soc_cols)
  !*  INPUTS
  !*    o n_rows                - Number of rows for scalar-relativistic eigenvectors
  !*    o n_cols                - Number of columns for scalar-relativsitic eigenvectors
  !*    o eigenvec              - Scalar-relativistic KS eigenvectors in real case
  !*    o eigenvec_complex      - Scalar-relativistic KS eigenvectors in complex case
  !*    o n_wf_bas_rows         - Number of rows for SOC eigenvectors in SR eigenvector basis
  !*    o n_wf_bas_rows         - Number of columns for SOC eigenvectors in SR eigenvector basis
  !*    o eigenvec_soc_wf_basis - The wavefunctions output by the second-variational step, which are
  !*                              in the basis set of unperturbed scalar-relativistic eigenvectors.
  !*                              First index is the basis element (which are unperturbed
  !*                              wavefunctions), second index is the SOC-perturbed state index.
  !*    o n_wf_bas_rows         - Number of rows for SOC eigenvectors in computational basis
  !*    o n_wf_bas_rows         - Number of columns for SOC eigenvectors in computational basis
  !*  OUTPUTS
  !*    o eigenvec_soc          - the SOC-perturbed wavefunctions in terms of the computational basis.
  !*                              First index is the basis element (which are computational basis
  !*                              elements), second index is the SOC-perturbed state index
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  NOTES
  !*    This subroutine applies for both LAPACK and ScaLAPACK branches of the code.
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! Iterators
    integer :: i_basis, i_state, j_state, i_spin, spin_index, basis_offset, state_offset

    complex*16, dimension(:,:,:), allocatable :: temp_eigenvec_complex

    character*300                  :: info_str
    character(*), parameter :: func = 'convert_wf_basis_to_compute_basis'

    ! WPH: Create a complex copy of the eigenvectors, if they are not already complex
    !      This shouldn't be necessary.  There should be a way to play with real/imaginary
    !      parts individually to work with real eigenvectors and eliminate this matrix copy.
    if (real_eigenvectors) then
      call aims_allocate(temp_eigenvec_complex, n_eigenvec_rows, n_eigenvec_cols, n_spin, "temp_eigenvec_complex")
      temp_eigenvec_complex = dcmplx(eigenvec)
    else
      call aims_allocate(temp_eigenvec_complex, 1, 1, 1, "temp_cmplx_eienvec")
    end if

    ! Code is written assuming that one is using all possible scalar-relativistic basis
    ! functions and spinors
    if (n_basis_soc .ne. 2*n_basis) then
      write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
      call aims_stop(info_str, func)
    end if
    if (mod(n_basis_soc_coll,2) .eq. 1) then
      call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
    end if
    if (n_basis_soc_ncoll .gt. 0) then
      call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                     & in SOC, exiting.', func)
    end if
    if (soc_saved_state_start .le. n_core_states_omit_from_soc) then
      call aims_stop('More SOC states are requested to be saved than were calculated, this&
                     & indicates something went wrong.  Exiting.', func)
    end if

    ! Specifies which states should be saved out of all the states that were
    ! calculated.
    state_offset = (soc_saved_state_start - 1) - n_core_states_omit_from_soc

    eigenvec_soc = (0.0d0, 0.0d0)
    do i_spin = 1, 2
      if (n_spin .eq. 1) then
        spin_index = 1
      else
        spin_index = i_spin
      end if
      if (i_spin .eq. 1) then
        basis_offset = 0
      else
        basis_offset = n_basis_soc_coll/2
      end if

      if (use_scalapack) then
        if (real_eigenvectors) then
          call pzgemm('N', 'N', &
                      n_basis_soc_coll/2, n_saved_states_soc, n_states_sr, &
                      (1.d0,0.d0), &
                      temp_eigenvec_complex(:,:,spin_index), 1,sr_state_start,sc_desc, &
                      eigenvec_soc_wf_basis, 1+(i_spin-1)*n_states_sr,1+state_offset,sc_desc_soc,&
                      (1.d0,0.d0), &
                      eigenvec_soc, basis_offset+1,1,sc_desc_soc_vec )
        else
          call pzgemm('N', 'N', &
                      n_basis_soc_coll/2, n_saved_states_soc, n_states_sr, &
                      (1.d0,0.d0), &
                      eigenvec_complex(:,:,spin_index), 1,sr_state_start,sc_desc, &
                      eigenvec_soc_wf_basis, 1+(i_spin-1)*n_states_sr,1+state_offset,sc_desc_soc,&
                      (1.d0,0.d0), &
                      eigenvec_soc, basis_offset+1,1,sc_desc_soc_vec )
        end if
      else
        if (real_eigenvectors) then
          call zgemm("N", "N", &
                     n_basis_soc_coll/2, n_saved_states_soc, n_states_sr, &
                     (1.0d0,0.0d0), &
                     temp_eigenvec_complex(1,sr_state_start,spin_index), n_eigenvec_rows, &
                     eigenvec_soc_wf_basis(1+(i_spin-1)*n_states_sr,1+state_offset), n_wf_bas_rows, &
                     (1.0d0,0.0d0), &
                     eigenvec_soc(basis_offset+1,1), n_eigenvec_soc_rows)
        else
          call zgemm("N", "N", &
                     n_basis_soc_coll/2, n_saved_states_soc, n_states_sr, &
                     (1.0d0,0.0d0), &
                     eigenvec_complex(1,sr_state_start,spin_index), n_eigenvec_rows, &
                     eigenvec_soc_wf_basis(1+(i_spin-1)*n_states_sr,1+state_offset), n_wf_bas_rows, &
                     (1.0d0,0.0d0), &
                     eigenvec_soc(basis_offset+1,1), n_eigenvec_soc_rows)
        end if
      end if
    end do

    call aims_deallocate(temp_eigenvec_complex, "temp_eigenvec_complex")
  end subroutine convert_wf_basis_to_compute_basis
  !******

  !****f* soc_utilities/calculate_spin_expectation
  !*  NAME
  !*    calculate_spin_expectation
  !*  SYNOPSIS
  subroutine calculate_spin_expectation( n_states_sr , KS_eigenvector_non_collinear_wf_basis, spin_expectation )
  !*  PURPOSE
  !*    Calculates the expectation values of non-collinear eigenvectors for a set of five spin operators:
  !*    1)   Projector onto spin-up (z-axis)
  !*    2)   Projector onto spin-down (z-axis)
  !*    3-5) Pauli spin matrices (sigma_x, sigma_y, sigma_z)
  !*  USES
    implicit none
  !*  ARGUMENTS
    integer, intent(in) :: n_states_sr
    complex*16, intent(in) :: KS_eigenvector_non_collinear_wf_basis(2*n_states_sr, 2*n_states_sr)
    real*8, intent(out) :: spin_expectation(2*n_states_sr, 5)
  !*  INPUTS
  !*    o n_states_sr                            - the number of SR states used.  (This is n_states in most of aims.)
  !*    o KS_eigenvector_non_collinear_wf_basis: - The eigenvectors expressed in a non-collinear expansion.
  !*                                               The first number_of_states rows are spin-up basis elements,
  !*                                               and the last number_of_states rows are spin-down basis elements.
  !*  OUTPUTS
  !*    o spin_expectation - The expectation values of spin operators.
  !*  NOTES
  !*    This subroutine is currently not used, but may be in the future.
  !*  TO-DO
  !*    o Update this subroutine to support BLACS, currently assumes that every process has the full matrix
  !*  AUTHORS
  !*    William Huhn (Duke University)
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    integer :: i_state, j_state
    complex*16 :: up_star_up, up_star_down, down_star_down

    spin_expectation = 0.0d0

    ! Calculate the expectation values of spin operators for every eigenstate
    do j_state = 1,2*n_states_sr
      ! Sum up the needed combinations of coefficients for the SOC-perturbed eigenvectors in
      ! the basis of unperturbed eigenstates
      ! In this basis, the overlap matrix is unity (since unperturbed eigenstates are orthnormal),
      ! making the sum trivial
      up_star_up     = (0.0d0,0.0d0)
      up_star_down   = (0.0d0,0.0d0)
      down_star_down = (0.0d0,0.0d0)
      do i_state = 1, n_states_sr
        up_star_up     = up_star_up     + dconjg(KS_eigenvector_non_collinear_wf_basis(i_state, j_state)) * &
                              KS_eigenvector_non_collinear_wf_basis(i_state, j_state)
        up_star_down   = up_star_down   + dconjg(KS_eigenvector_non_collinear_wf_basis(i_state, j_state)) * &
                              KS_eigenvector_non_collinear_wf_basis(i_state + n_states_sr, j_state)
        down_star_down = down_star_down + dconjg(KS_eigenvector_non_collinear_wf_basis(i_state + n_states_sr, j_state)) * &
                              KS_eigenvector_non_collinear_wf_basis(i_state + n_states_sr, j_state)
      end do
      ! Expectation value of spin-up (z-axis) projector
      spin_expectation(j_state, 1) = real( up_star_up )
      ! Expectation value of spin-down (z-axis) projector
      spin_expectation(j_state, 2) = real( down_star_down )
      ! Expectation value of sigma_x operator
      spin_expectation(j_state, 3) = 2.0d0 *  real( up_star_down )
                                 ! = real( up_star_down + dconjg(up_star_down) )
      ! Expectation value of sigma_x operator
      spin_expectation(j_state, 4) = 2.0d0 * aimag( up_star_down )
                                 ! = real( (0.0d0,-1.0d0)*up_star_down + (0.0d0,1.0d0)*dconjg(up_star_down) )
      ! Expectation value of sigma_z operator
      spin_expectation(j_state, 5) = real( up_star_up - down_star_down )
    end do
  end subroutine calculate_spin_expectation
  !******

  !****f* soc_utilities/find_num_core_states_from_energy
  !*  NAME
  !*    find_num_core_states_from_energy
  !*  SYNOPSIS
  subroutine find_min_energy_from_gap( gap, KS_eigenvalue, &
       occ_numbers, min_energy )
  !*  PURPOSE
  !*    Determines what the minimum energy threshhold for core states should be
  !*    based on a user-provided separation between energy eigenvalues.
  !*
  !*    Note that we are creating indices for SOC-perturbed quantities but
  !*    approximating them from scalar-relativistic quantities (as, when this
  !*    function is called, we don't have the SOC-perturbed quantities yet).
  !*  USES
    use localorb_io, only: localorb_info
    use dimensions, only: n_states, n_spin, n_k_points
    use constants, only: hartree
    use physics, only: n_electrons
    implicit none
  !*  ARGUMENTS
    real*8, intent(in)  :: gap
    real*8, intent(in)  :: KS_eigenvalue(n_states, n_spin, n_k_points)
    real*8, intent(in)  :: occ_numbers(n_states, n_spin, n_k_points)
    real*8, intent(out) :: min_energy
  !*  INPUTS
  !*    o gap            - Threshhold energy for gaps between states
  !*    o KS_eigenvalue  - The scalar-relativistic KS eigenvalues across all k-points
  !*    o occ_numbers    - The scalar-relativistic occupation numbers across all k-points
  !*  OUTPUTS
  !*    o min_energy     - Minimum energy for energy window
  !*  AUTHORS
  !*    William Huhn
  !*  HISTORY
  !*    August 2017, created
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! local variables
    integer :: i_spin, i_state, i_k_point
    integer :: n_core_states, n_core_states_last_spin, n_core_states_last_k
    real*8  :: eigenvalue, last_eigenvalue
    real*8  :: upper_eigenvalue_gap, lower_eigenvalue_gap
    real*8  :: failsafe_min_energy

    character*300                  :: info_str

    character(*), parameter :: func = 'find_min_energy_from_gap'

    ! In this algorithm, we do not stop aims when something goes wrong, because
    ! it may be due to the algorithm (i.e. me) and not PEBKAC.  Instead, we
    ! set a minimum energy that's sufficiently low to lie below every
    ! eigenvalue, essentially including all states in the energy window.
    failsafe_min_energy = KS_eigenvalue(1,1,1)*hartree - 100000

    ! Ideally, what should be done when the algorithm fails is we should repeat
    ! the algorithm, but stopping before encountering the set of bands which
    ! were selected by the previous application of the algorithm and which
    ! caused the algorithm to fail.

    if (gap .le. 0.0d0) then
      write(info_str,'(2X,A)') "A zero or negative value was specified for&
           & the gap between states.&
           & Setting minimum energy sufficiently low to include all states."
      call localorb_info( info_str )
      min_energy = failsafe_min_energy
      return
    end if

    ! Although the desired quantity is the minimum energy threshhold, it is
    ! more convenient to work with the state index for the highest-lying state
    ! below this threshhold (whereas the minimum energy threshhold can lie
    ! anywhere in the gap)

    do i_k_point = 1, n_k_points
      ! Save number of core states predicted from previous k-point
      if (i_k_point .gt. 1) then
        n_core_states_last_k = n_core_states
      end if

      do i_spin = 1, n_spin
        ! Save number of core states predicted from previous spin channel
        if (i_spin .gt. 1) then
          n_core_states_last_spin = n_core_states
        end if

        ! Find number of core states for this k-point and spin channel
        ! It is assumed that the eigenvalues are already sorted.
        last_eigenvalue = KS_eigenvalue(1, i_spin, i_k_point)*hartree
        do i_state = 2, n_states
          eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_point)*hartree
          if (eigenvalue - last_eigenvalue .gt. gap) then
            n_core_states        = i_state - 1
            upper_eigenvalue_gap = eigenvalue
            lower_eigenvalue_gap = last_eigenvalue
          end if
          last_eigenvalue = eigenvalue
        end do
        if (i_k_point .eq. 1 .and. i_spin .eq. 1) then
          ! Now set the minimum energy to lie at 0.5*gap above the eigenvalue
          ! for the predicted core state at the first k-point and in the first
          ! spin channel.  This is the value that will be returned.
          ! We will continue iterating over all k-points and spin channels to
          ! ensure that this value holds for all k-points and spin channels,
          ! as it should as long as the provided gap actually does pick out
          ! core states (i.e. flat and fixed across k-points and spin channels)
          min_energy = KS_eigenvalue(n_core_states,1,1)*hartree + 0.5*gap
        end if

        if (n_core_states .ge. n_states) then
          write(info_str,'(1X,A)') '* Predicted number of core states using gap&
               & criterion is greater than or equal to the total number of&
               & states.&
               & Setting minimum energy sufficiently low to include all states.'
          call localorb_info( info_str )
          min_energy = failsafe_min_energy
          return
        end if

        ! If core is spin-polarized (i.e. predicted number of core states differs
        ! between spin channels), exit
        if (i_spin .gt. 1 .and. n_core_states_last_spin .ne. n_core_states) then
          write(info_str,'(1X,A)') '* Predicted number of core states using gap&
               & criterion had a different number of states in the core between&
               & the spin channels.&
               & Setting minimum energy sufficiently low to include all states.'
          call localorb_info( info_str )
          min_energy = failsafe_min_energy
          return
        end if

        ! This is necessary but not sufficient; the user may still remove too
        ! many core states for determination of Fermi level in metallic systems.
        if (dble(2*n_core_states) .ge. n_electrons) then
          write(info_str,'(1X,A)') '* Predicted number of core states using gap&
               & criterion would eliminate all occupied states from calculation,&
               & preventing determination of Fermi level.&
               & Setting minimum energy sufficiently low to include all states.'
          call localorb_info( info_str )
          min_energy = failsafe_min_energy
          return
        end if

        ! Make sure that the minimum energy (obtained for the first k-point and
        ! spin channel) lies between the boundary eigenvalues for all k-points
        ! and spin channels, i.e. the bands aren't too curved.
        if (min_energy .lt. lower_eigenvalue_gap .or. &
             min_energy .gt. upper_eigenvalue_gap) then
          write(info_str,'(1X,A)') '* Predicted minimum energy for energy window&
               & using gap criterion does not hold for all k-points or spin&
               & channels.  This is likely due to your gap criterion being too&
               & small and/or the relevant bands having too large a curvage.&
               & Setting minimum energy sufficiently low to include all states.'
          call localorb_info( info_str )
          min_energy = failsafe_min_energy
          return
        end if
      end do

      ! If prediction of number of core states differs between k-points, exit
      if (i_k_point .gt. 1 .and. n_core_states_last_k .ne. n_core_states) then
        write(info_str,'(1X,A)') '* Predicted number of core states using gap&
             & criterion had a differing number of states for different k-points&
             &. Setting minimum energy sufficiently low to include all states.'
         call localorb_info( info_str )
         min_energy = failsafe_min_energy
         return
      end if
    end do
  end subroutine find_min_energy_from_gap
  !******

  !****f* soc_utilities/find_num_core_states_from_energy
  !*  NAME
  !*    find_num_core_states_from_energy
  subroutine find_num_core_states_from_energy( min_energy, KS_eigenvalue, &
       occ_numbers, n_core_states )
  !*  PURPOSE
  !*    Find the number of states with energy less than a specified threshhold.
  !*    In principle, trivial, but error checking needs to be done to make sure
  !*    the user-provided energy makes sense.
  !*
  !*    Note that we are creating indices for SOC-perturbed quantities but
  !*    approximating them from scalar-relativistic quantities (as, when this
  !*    function is called, we don't have the SOC-perturbed quantities yet).
  !*  USES
    use mpi_tasks, only: aims_stop
    use dimensions, only: n_states, n_spin, n_k_points
    use constants, only: hartree
    use physics, only: n_electrons
    implicit none
  !*  ARGUMENTS
    real*8,  intent(in)  :: min_energy
    real*8,  intent(in)  :: KS_eigenvalue(n_states, n_spin, n_k_points)
    real*8,  intent(in)  :: occ_numbers(n_states, n_spin, n_k_points)
    integer, intent(out) :: n_core_states
  !*  INPUTS
  !*    o min_energy     - Threshhold energy for core states
  !*    o KS_eigenvalue  - The scalar-relativistic KS eigenvalues across all k-points
  !*    o occ_numbers    - The scalar-relativistic occupation numbers across all k-points
  !*  OUTPUTS
  !*    o n_core_states  - Number of core states for SOC-perturbed quantities
  !*  AUTHORS
  !*    William Huhn
  !*  HISTORY
  !*    August 2017, created
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! local variables
    integer :: i_spin, i_state, i_k_point
    integer :: n_core_states_last_spin, n_core_states_last_k
    real*8  :: eigenvalue, last_eigenvalue

    character*300                  :: info_str

    character(*), parameter :: func = 'find_num_core_states_from_energy'

    do i_k_point = 1, n_k_points
      ! Save number of core states predicted from previous k-point
      if (i_k_point .gt. 1) then
        n_core_states_last_k = n_core_states
      end if

      do i_spin = 1, n_spin
        ! Save number of core states predicted from previous spin channel
        if (i_spin .gt. 1) then
          n_core_states_last_spin = n_core_states
        end if

        ! Find number of core states for this k-point and spin channel
        ! It is assumed that the eigenvalues are already sorted.
        last_eigenvalue = KS_eigenvalue(1, i_spin, i_k_point)*hartree
        if (last_eigenvalue .ge. min_energy) then
          n_core_states = 0
        else
          do i_state = 2, n_states
            eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_point)*hartree
            if (last_eigenvalue .lt. min_energy .and. eigenvalue .ge. min_energy) then
              n_core_states = i_state - 1
              exit
            end if
            last_eigenvalue = eigenvalue
          end do
        end if
        if (n_core_states .ge. n_states) then
          write(info_str,'(1X,A)') '* Predicted number of core states is greater&
               & than or equal to the total number of states.  Exiting.'
          call aims_stop(info_str, func)
        end if

        ! If core is spin-polarized (i.e. predicted number of core states differs
        ! between spin channels), exit
        if (i_spin .gt. 1 .and. n_core_states_last_spin .ne. n_core_states) then
          write(info_str,'(1X,A)') '* Specified energy threshhold for core states &
               &had a differing numbers of states in the core betweem the spin&
               &channels.  Exiting.'
          call aims_stop(info_str, func)
        end if

        ! This is necessary but not sufficient; the user may still remove too
        ! many core states for determination of Fermi level in metallic systems.
        if (dble(2*n_core_states) .ge. n_electrons) then
          write(info_str,'(1X,A)') '* Predicted number of core states would&
               & eliminate all occupied states from calculation, preventing&
               & determination of Fermi level.  Exiting.'
          call aims_stop(info_str, func)
        end if
      end do

      ! If prediction of number of core states differs between k-points, exit
      if (i_k_point .gt. 1 .and. n_core_states_last_k .ne. n_core_states) then
        write(info_str,'(1X,A)') '* Specified energy threshhold for core states &
             &had a differing number of states for different k-points.  Exiting.'
        call aims_stop(info_str, func)
      end if
    end do

    ! Convert from SR to SOC
    n_core_states = 2*n_core_states

  end subroutine find_num_core_states_from_energy
  !******

  !****f* soc_utilities/find_num_high_states_from_energy
  !*  NAME
  !*    find_num_high_states_from_energy
  subroutine find_num_high_states_from_energy( max_energy, KS_eigenvalue, &
       occ_numbers, n_high_states )
  !*  SYNOPSIS
  !*    call find_num_high_states_from_energy( max_energy, KS_eigenvalue, &
  !*         n_high_states )
  !*  FUNCTION
  !*    Find the number of states with energy greater than a specified
  !*    threshhold.  In principle, trivial, but error checking needs to be done
  !*    to make sure the user-provided energy makes sense.
  !*
  !*    Note that we are creating indices for SOC-perturbed quantities but
  !*    approximating them from scalar-relativistic quantities (as, when this
  !*    function is called, we don't have the SOC-perturbed quantities yet).
  !*  USES
    use mpi_tasks, only: aims_stop
    use dimensions, only: n_states, n_spin, n_k_points
    use constants, only: hartree
    use physics, only: n_electrons
    implicit none
  !*  ARGUMENTS
    real*8,  intent(in)  :: max_energy
    real*8,  intent(in)  :: KS_eigenvalue(n_states, n_spin, n_k_points)
    real*8,  intent(in)  :: occ_numbers(n_states, n_spin, n_k_points)
    integer, intent(out) :: n_high_states
  !*  INPUTS
  !*    o max_energy     - Threshhold energy for high-lying states
  !*    o KS_eigenvalue  - The scalar-relativistic KS eigenvalues across all k-points
  !*    o occ_numbers    - The scalar-relativistic occupation numbers across all k-points
  !*  OUTPUTS
  !*    o n_high_states  - Number of high-lying states for SOC-perturbed quantities
  !*  AUTHORS
  !*    William Huhn
  !*  HISTORY
  !*    August 2017, created
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! local variables
    integer :: i_spin, i_state, i_k_point
    integer :: n_high_states_last_spin, n_high_states_last_k
    real*8  :: eigenvalue, last_eigenvalue

    character*300                  :: info_str

    character(*), parameter :: func = 'find_num_high_states_from_energy'

    do i_k_point = 1, n_k_points
      ! Save number of high-lying states predicted from previous k-point
      if (i_k_point .gt. 1) then
        n_high_states_last_k = n_high_states
      end if

      do i_spin = 1, n_spin
        ! Save number of high-lying states predicted from previous spin channel
        if (i_spin .gt. 1) then
          n_high_states_last_spin = n_high_states
        end if

        ! Find number of high-lying states for this k-point and spin channel
        ! It is assumed that the eigenvalues are already sorted.
        last_eigenvalue = KS_eigenvalue(n_states, i_spin, i_k_point)*hartree
        if (last_eigenvalue .le. max_energy) then
          n_high_states = 0
        else
          do i_state = n_states - 1, 1, -1
            eigenvalue = KS_eigenvalue(i_state, i_spin, i_k_point)*hartree
            if (last_eigenvalue .gt. max_energy .and. eigenvalue .le. max_energy) then
              n_high_states = n_states - i_state
              exit
            end if
            last_eigenvalue = eigenvalue
          end do
        end if
        if (n_high_states .ge. n_states) then
          write(info_str,'(1X,A)') '* Predicted number of high-lying states is greater&
               & than or equal to the total number of states.  Exiting.'
          call aims_stop(info_str, func)
        end if

        ! If differing number of states in spin channels, exit
        if (i_spin .gt. 1 .and. n_high_states_last_spin .ne. n_high_states) then
          write(info_str,'(1X,A)') '* Specified energy threshhold for high-lying states &
               &had a differing number of states in spin channel.  Exiting.'
          call aims_stop(info_str, func)
        end if

        ! This is necessary but not sufficient; the user may still remove too
        ! many high-lying states for determination of Fermi level in metallic
        ! systems.
        if (dble(2*n_states - 2*n_high_states) .le. n_electrons) then
          write(info_str,'(1X,A)') '* Predicted number of high-lying states would&
               & eliminate occupied states from calculation, preventing&
               & determination of Fermi level.  Exiting.'
          call aims_stop(info_str, func)
        end if
      end do

      ! If prediction of number of high-lying states differs between k-points, exit
      if (i_k_point .gt. 1 .and. n_high_states_last_k .ne. n_high_states) then
        write(info_str,'(1X,A)') '* Specified energy threshhold for high-lying states &
             &had a differing number of states for different k-points.  Exiting.'
        call aims_stop(info_str, func)
      end if
    end do

    ! Convert from SR to SOC
    n_high_states = 2*n_high_states

  end subroutine find_num_high_states_from_energy
  !******

  !****f* soc_utilities/create_sorted_sr_eigenvalues_list
  !*  NAME
  !*    create_sorted_sr_eigenvalues_list
  !*  SYNOPSIS
  subroutine create_sorted_sr_eigenvalues_list( n_states_sr, sr_state_start, &
       KS_eigenvalue, sorted_eigenvalues )
  !*  PURPOSE
  !*    Create a list of sorted SR eigenvalues within a window of eigenvalues,
  !*    merging the spin channels together to form one continuous array.
  !*  USES
    use localorb_io
    use dimensions, only: n_states, n_spin
    implicit none
  !*  ARGUMENTS
    integer, intent(in)  :: n_states_sr
    integer, intent(in)  :: sr_state_start
    real*8,  intent(in)  :: KS_eigenvalue(n_states, n_spin)
    real*8,  intent(out) :: sorted_eigenvalues(2*n_states_sr)
  !*  INPUTS
  !*    o n_states_sr        - Number of SR states to include (does not include spin index)
  !*    o sr_state_start     - Starting index for SR states to include
  !*    o KS_eigenvalue      - The KS eigenvalues for the current k-point.
  !*  OUTPUTS
  !*    o sorted_eigenvalues - The list of sorted eigenvalues over both spin channels
  !*  AUTHORS
  !*    William Huhn
  !*  HISTORY
  !*    August 2017, created
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! local variables
    integer       :: i_state, j_state
    character*300 :: info_str
    real*8        :: temp

    character(*), parameter :: func = 'create_sorted_sr_eigenvalues_list'

    ! Create a mostly ordered list of eigenvalues including both spin
    ! channels (exactly ordered for spin-non-polarized calculations)
    do i_state = 1, n_states_sr, 1
      sorted_eigenvalues(2*i_state-1) = KS_eigenvalue(i_state+sr_state_start-1,1)
    end do
    do i_state = 1, n_states_sr, 1
      if (n_spin.eq.2) then
        sorted_eigenvalues(2*i_state) = KS_eigenvalue(i_state+sr_state_start-1,2)
      else
        sorted_eigenvalues(2*i_state) = KS_eigenvalue(i_state+sr_state_start-1,1)
      end if
    end do

    ! Yes, it's bubblesort on a mostly-ordered list.  Oh well.
    do i_state = 1, 2*n_states_sr, 1
      do j_state = i_state + 1, 2*n_states_sr, 1
        if (sorted_eigenvalues(i_state) > &
             sorted_eigenvalues(j_state)) then
          temp = sorted_eigenvalues(i_state)
          sorted_eigenvalues(i_state) = sorted_eigenvalues(j_state)
          sorted_eigenvalues(j_state) = temp
        end if
      end do
    end do

  end subroutine create_sorted_sr_eigenvalues_list
  !******

  !****f* soc_utilities/create_sorted_sr_eigenvalues_idxmap
  !*  NAME
  !*    create_sorted_sr_eigenvalues_idxmap
  !*  SYNOPSIS
  subroutine create_sorted_sr_eigenvalues_idxmap(n_states_sr, sr_state_start, &
       KS_eigenvalue, sr_to_soc_idxmap )
  !*  PURPOSE
  !*    Create a list of mapping indices by sorting the SR eigenvalues within a
  !*    window of eigenvalues,  merging the spin channels together to form one
  !*    continuous array.
  !*  USES
    use localorb_io
    implicit none
  !*  ARGUMENTS
    integer, intent(in)  :: n_states_sr
    integer, intent(in)  :: sr_state_start
    real*8,  intent(in)  :: KS_eigenvalue(n_states_sr_save,n_spin_sr_save)
    integer, intent(out) :: sr_to_soc_idxmap(n_states_sr)
  !*  INPUTS
  !*    o n_states_sr        - Number of SR states to include (does not include spin index)
  !*    o sr_state_start     - Starting index for SR states to include
  !*    o KS_eigenvalue      - The KS eigenvalues for the current k-point.
  !*  OUTPUTS
  !*    o sr_to_soc_idxmap - The list of indices mapping the soc-KS_eigenvals to sr
  !*  AUTHORS
  !*    Georg Michelitsch
  !*  HISTORY
  !*    Feb 2018
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! local variables
    integer       :: i_state, j_state, n_states_local
    character*300 :: info_str
    real*8        :: temp
    integer       :: tempidx
    real*8        :: sorted_eigenvalues(n_states_sr)

    character(*), parameter :: func = 'create_sorted_sr_eigenvalues_idxmap'

    ! adjust to the local environment and initialize arrays
    n_states_local = n_states_sr / 2
    sr_to_soc_idxmap = 0

    ! Create a mostly ordered list of eigenvalues including both spin
    ! channels (exactly ordered for spin-non-polarized calculations)
    do i_state = 1, n_states_local, 1
      sorted_eigenvalues(2*i_state-1) = KS_eigenvalue(i_state+sr_state_start-1,1)
      sr_to_soc_idxmap(2*i_state-1) = i_state+sr_state_start-1
    end do
    do i_state = 1, n_states_local, 1
      if (n_spin_sr_save.eq.2) then
        sorted_eigenvalues(2*i_state) = KS_eigenvalue(i_state+sr_state_start-1,2)
        sr_to_soc_idxmap(2*i_state) = n_states_local + i_state+sr_state_start-1
      else
        sorted_eigenvalues(2*i_state) = KS_eigenvalue(i_state+sr_state_start-1,1)
        sr_to_soc_idxmap(2*i_state) = n_states_local + i_state+sr_state_start-1
      end if
    end do

    ! This is *exactly* the same bubble-sort that Will implemented
    do i_state = 1, n_states_sr, 1
      do j_state = i_state + 1, n_states_sr, 1
        if (sorted_eigenvalues(i_state) > &
             sorted_eigenvalues(j_state)) then
          temp = sorted_eigenvalues(i_state)
          tempidx = sr_to_soc_idxmap(i_state)
          sorted_eigenvalues(i_state) = sorted_eigenvalues(j_state)
          sr_to_soc_idxmap(i_state) = sr_to_soc_idxmap(j_state)
          sorted_eigenvalues(j_state) = temp
          sr_to_soc_idxmap(j_state) = tempidx
        end if
      end do
    end do

  end subroutine create_sorted_sr_eigenvalues_idxmap
  !******


  !****f* soc_utilities/perform_soc_perturbation
  !*  NAME
  !*    perform_soc_pertubation
  !*  SYNOPSIS
  subroutine perform_soc_perturbation( n_rows, n_cols, perturb_hamil, &
                                       KS_eigenvalue, KS_eigenvalue_perturb, &
                                       perturb_matrix_max, n_wf_bas_rows, &
                                       n_wf_bas_cols, eigenvec_soc_wf_basis )
  !*  PURPOSE
  !*    Perform second-variational perturbations by adding the eigenvalues to
  !*    the diagonal of perturb_hamil, generating the perturbed SOC Hamiltonian
  !*    (called perturb_mat) and then diagonalize the Hamiltonian using ELSI,
  !*    ELPA (2013 version), or LAPACK.
  !*
  !*    This function unites both the LAPACK and ScaLAPACK code paths; the code
  !*    will interpret how to index perturb_hamil based on use_scalapack
  !*  USES
    use localorb_io, only: OL_norm, use_unit, localorb_info
    use dimensions, only: n_states, n_spin
    use runtime_choices, only: use_elpa, use_elsi, use_scalapack
    use applicable_citations, only: cite_reference
    use scalapack_wrapper, only: mpi_comm_cols, mpi_comm_rows, my_blacs_ctxt, &
        my_scalapack_comm_all, my_scalapack_comm_work, solver_method_used
    use mpi_tasks, only: aims_stop
    use aims_memory_tracking, only: aims_allocate, aims_deallocate
    use physics, only: n_electrons
    use elsi_wrapper, only: aims_elsi_stdevp
    use scalapack_soc
    use dimensions_soc, only: n_states_sr, n_states_soc, sr_state_start
    implicit none
  !*  ARGUMENTS
    integer, intent(in)  :: n_rows
    integer, intent(in)  :: n_cols
    complex*16, intent(in)  :: perturb_hamil(n_rows, n_cols)
    real*8, intent(in)  :: KS_eigenvalue(n_states, n_spin)
    real*8, intent(out) :: KS_eigenvalue_perturb(n_states_soc)
    real*8, intent(out) :: perturb_matrix_max
    integer, intent(in)  :: n_wf_bas_rows
    integer, intent(in)  :: n_wf_bas_cols
    complex*16, intent(out) :: eigenvec_soc_wf_basis(n_wf_bas_rows, n_wf_bas_cols)
  !*  INPUTS
  !*    o n_rows          - Number of rows for perturb_hamil.   In the LAPACK case, should be n_states_soc.  In the ScaLAPACK case,
  !*                        should be mxld_soc
  !*    o n_cols          - Number of columns for perturb_hamil.   In the LAPACK case, should be n_states_soc.  In the ScaLAPACK case,
  !*                        should be mxcol_soc
  !*    o perturb_hamil   - The matrix elements of the SOC opertor over the unperturbed eigenvectors, <Psi|V_SOC|Phi>.
  !*                        Whether this is a full/global matrix or local/2D block cyclic matrix depends on use_scalapack
  !*    o KS_eigenvalue   - The unperturbed KS eigenvalues for a single k-point.  Both perturb_hamil and KS_eigenvalue should
  !*                        correspond to the same k-point.
  !*    o n_wf_bas_rows   - Number of rows for eigenvector_soc_wf_bas.  In the LAPACK case, should be n_states_soc.  In the ScaLAPACK
  !*                        case, should be mxld_soc.
  !*    o n_wf_bas_cols   - Number of columns for eigenvector_soc_wf_bas.   In the LAPACK case, should be n_states_soc.  In the
  !*                        ScaLAPACK case, should be mxcol_soc.
  !*  OUTPUTS
  !*    o KS_eigenvalue_perturb    - The SOC-perturbed KS eigenvalues.
  !*    o perturb_matrix_max       - The maximum absolute value for elements of the perturbation matrix
  !*    o eigenvector_soc_wf_basis - The coefficients of the SOC-perturbed eigenvectors in the basis set of the unperturbed
  !*                                 eigenvectors.
  !*  AUTHORS
  !*    William Huhn
  !*  HISTORY
  !*    November 2016, Refactor of perform_qd_perturbation
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  NOTES
  !*    The implementation of second-variational SOC in FHI-aims is published in
  !*      Huhn and Blum, Phys. Rev. Materials 1, 033803 (2017)
  !*      https://dx.doi.org/10.1103/PhysRevMaterials.1.033803
  !*
  !*    This subroutine implements steps 5-6 (non-periodic) and steps 6-7 &
  !*    (periodic) in Section III.3 from said paper.
  !*  SOURCE

  ! local variables
    integer ::  temp_int_1, temp_int_2, i_state, j_state
    integer ::  i, j, k, i_index_real
    real*8  ::  temp, dummy_1, dummy_2, temp_ev_diff
    integer :: up_index, down_index
    character*300 :: info_str

 ! Variables for LAPACK
    ! The perturbed Hamiltonian in the LAPACK case
    complex*16, dimension(:), allocatable :: perturbation_matrix
    integer, dimension(:), allocatable    :: ifail
    real*8, dimension(:), allocatable     :: work
    integer, dimension(:), allocatable    :: iwork
    complex*16, dimension(:), allocatable :: zwork
    integer :: n_found
    integer :: info = 0

  ! Varibles for ScaLAPACK/ELPA
    ! The local matrix of the perturbed Hamiltonian in the ScaLAPACK case
    complex*16, dimension(:,:), allocatable :: perturb_mat

   ! WPH:  Once the perturbation's matrix elements are calculated, perturbation
   !       theory doesn't care whether the system is periodic or cluster;  it
   !       only "sees" an unperturbed eigenvector basis.  So (thankfully) the
   !       periodic and cluster codes are identical

    character(*), parameter :: func = 'perform_soc_perturbation'

    if (.not.use_scalapack) then
      ! Allocate work arrays needed for LAPACK eigensolver
      if (.not.allocated(work)) then
        call aims_allocate(work,  7*n_states_soc, "work")
      end if
      if (.not.allocated(iwork)) then
        call aims_allocate(iwork, 5*n_states_soc, "iwork")
      end if
      if (.not.allocated(ifail)) then
        call aims_allocate(ifail, n_states_soc, "ifail")
      end if
      if (.not.allocated(zwork)) then
        call aims_allocate(zwork, 2*n_states_soc, "zwork")
      end if
   else
      if (.not.allocated(work)) then
        call aims_allocate(work,  1, "work")
      end if
      if (.not.allocated(iwork)) then
        call aims_allocate(iwork, 1, "iwork")
      end if
      if (.not.allocated(ifail)) then
        call aims_allocate(ifail, 1, "ifail")
      end if
      if (.not.allocated(zwork)) then
        call aims_allocate(zwork, 1, "zwork")
      end if
   end if

    perturb_matrix_max = 0.0d0
    KS_eigenvalue_perturb = 0.0d0

    ! These are indices for the spin index of the eigenvectors, used to handle
    ! both "spin none" and "spin collinear" runs in the same framework.  They
    ! have the side benefit of making the code a little easier to convert back
    ! to physical equations.
    if (n_spin.eq.1) then
      up_index = 1
      down_index = 1
    else
      up_index = 1
      down_index = 2
    end if

    if ( use_scalapack .and. (use_elpa.or.use_elsi) ) then
      if ( (n_rows.ne.mxld_soc) .or. (n_cols.ne.mxcol_soc) ) then
       call aims_stop("Array indexing error in ScaLAPACK branch of &
                      &perform_soc_perturbation!  Exiting.")
      end if

      call aims_allocate(perturb_mat, mxld_soc, mxcol_soc, "perturb_mat")
      call aims_allocate(perturbation_matrix, 1, "perturbation_matrix")

      perturb_mat = (0.0d0, 0.0d0)

      ! Generate the SOC-perturbed Hamiltonian by adding eigenvalues to diagonal
      ! ELPA requires this be the unpacked, full matrix (though in dense 2D
      ! block cyclic form)

      ! Step 5 (non-periodic) and step 6 (periodic) from Section III.3
      ! in Huhn and Blum, Phys. Rev. Mater. (2017)
      do i_state = 1, n_rows, 1
        do j_state = 1, n_cols, 1
          if (my_row_soc(i_state) .gt. 0 .and. my_col_soc(j_state) .gt. 0) then
            if (my_row_soc(i_state) .eq. my_col_soc(j_state) .and. &
                my_row_soc(i_state) .le. n_states_sr) then
              ! Diagonal element for spin-up state
              perturb_mat(i_state, j_state) = perturb_hamil(i_state, j_state) &
                   + KS_eigenvalue( (sr_state_start-1) + my_row_soc(i_state), up_index )
            else if (my_row_soc(i_state) .eq. my_col_soc(j_state) .and. &
                     my_row_soc(i_state) .gt. n_states_sr) then
              ! Diagonal element for spin-down state
              perturb_mat(i_state, j_state) = perturb_hamil(i_state, j_state) &
                   + KS_eigenvalue( (sr_state_start-1) + my_row_soc(i_state) - n_states_sr, down_index )
            else
              ! Off-diagonal element
              perturb_mat(i_state, j_state) = perturb_hamil(i_state, j_state)
            end if

            if (abs(perturb_mat(i_state, j_state)) .gt. perturb_matrix_max) then
              perturb_matrix_max = abs(perturb_mat(i_state, j_state))
            end if
          end if
        end do
      end do
      call set_full_square_matrix_complex(n_states_soc, nb_soc, sc_desc_soc, &
                                          l_row_soc, l_col_soc, &
                                          mxld_soc, mxcol_soc, perturb_mat )

      ! Now diagonalize the Hamiltonian using the chosen solver
      !
      ! Step 6 (non-periodic) and step 7 (periodic) from Section III.3
      ! in Huhn and Blum, Phys. Rev. Mater. (2017)
      if (use_elsi .or. use_elpa) then
        call aims_elsi_stdevp(n_states_soc, mxld_soc, mxcol_soc, n_states_soc, &
                              perturb_mat, KS_eigenvalue_perturb, &
                              eigenvec_soc_wf_basis, my_blacs_ctxt, mb_soc, &
                              my_scalapack_comm_work)
      else
        call aims_stop("You have reached a branch in the SOC code that should &
                       &not be possible.  Exiting.", func)
      end if
    else ! LAPACK
      if ( (n_rows.ne.n_states_soc) .or. (n_cols.ne.n_states_soc) ) then
        call aims_stop("Array indexing error in LAPACK branch of &
                       &perform_soc_perturbation!  Exiting.", func)
      end if

      call aims_allocate(perturb_mat, 1,1, "perturb_mat")
      call aims_allocate(perturbation_matrix, n_states_soc*(n_states_soc+1)/2, &
           "perturbation_matrix")

      perturbation_matrix = (0.0d0,0.0d0)

      ! Pack perturbation matrix into upper triangular form that LAPACK expects
      ! while adding SR eigenvalues to diagonal
      ! The first n_states_sr states in perturbation_matrix should correspond to
      ! spin up, and the second n_states_sr states should correspond to spin down
      do i_state = 1, n_states_soc, 1
        do j_state = i_state, n_states_soc, 1
          i_index_real = i_state + (j_state-1)*j_state/2
          if (i_state .eq. j_state .and. i_state .le. n_states_sr) then
            ! Diagonal element for spin-up state
            perturbation_matrix(i_index_real) = perturb_hamil(i_state,i_state) + &
                 KS_eigenvalue((sr_state_start-1)+i_state,up_index)
          else if (i_state .eq. j_state .and. i_state .gt. n_states_sr) then
            ! Diagonal element for spin-down state
            perturbation_matrix(i_index_real) = perturb_hamil(i_state,i_state) + &
                 KS_eigenvalue((sr_state_start-1)+i_state-n_states_sr,down_index)
          else
            ! Off-diagonal element
            perturbation_matrix(i_index_real) = perturb_hamil(i_state,j_state)
          end if

          if (abs(perturbation_matrix(i_index_real)) .gt. perturb_matrix_max) then
            perturb_matrix_max = abs(perturbation_matrix(i_index_real))
          end if
        end do
      end do

      call zhpevx('V', 'A', 'U', &
                  n_states_soc, perturbation_matrix, &
                  dummy_1, dummy_2, temp_int_1, temp_int_2, 1.0d-12, &
                  n_found, KS_eigenvalue_perturb, &
                  eigenvec_soc_wf_basis, n_states_soc, &
                  zwork, work, iwork, ifail, info )
    end if ! use_scalapack ...

    ! Output failed LAPACK diagonalizations if requested
    ! In practice, this reported a number of errors for noise, so I don't use it.
    ! I am still loathe to remove it.
    ! CC: FIXME It should be either commented out or removed, 
    !           since the IF(.false.) is easily overseen
    !if (.not.use_scalapack .and. .false. .and. info .ne. 0) then
    if (.not.use_scalapack .and. info .ne. 0) then
  !    open(88,file=file_name)
      write(88,'(A)')        "********************************************************"
      write(88,'(A)')        "**                  WARNINGS                          **"
      write(88,'(A)')        "********************************************************"
      write(88,*)"In the following subspaces LAPACK failed to diagonalize the perturbation matrix"
      write(88,*)
      write(88,'(A,I5,A,I5)')"SOC-perturbed Hamiltonian with dimension ",&
           2*n_states
      write(88,*)"----------------------------------------------------------------"
      write(88,*)
      write(88,'(2X,A,I0)') "WARNING:  ERROR in diagonalizing PMatrix in SOC-perturbed Hamiltonian"
      write(88,'(2X,A,I0)') "ZHPEVX returns with info = ",info
      write(88,'(2X,A)')"Perturbation matrix elements: "
      write(88,*)"--------------------------------------------------"
      do i_state = 1, n_states_soc, 1
        do j_state = i_state, n_states_soc, 1
          i_index_real = i_state + (j_state-1)*j_state/2
          write(88,*) i_state, j_state, perturbation_matrix(i_index_real)
        end do
      end do
    end if

    ! Allocatable arrays that are tracked
    if (allocated(perturbation_matrix)) &
         call aims_deallocate(perturbation_matrix, "perturbation_matrix")
    if (allocated(work)) call aims_deallocate(work, "work")
    if (allocated(iwork)) call aims_deallocate(iwork, "iwork")
    if (allocated(ifail)) call aims_deallocate(ifail, "ifail")
    if (allocated(zwork)) call aims_deallocate(zwork, "zwork")
    if (allocated(perturb_mat)) call aims_deallocate(perturb_mat, "perturb_mat")

  end subroutine perform_soc_perturbation
  !******

  !****f* soc_utilities/write_soc_values
  !*  NAME
  !*    write_soc_values
  !*  SYNOPSIS
  subroutine write_soc_values(KS_eigenvalue, KS_eigenvalue_soc_perturbed, occ_numbers_soc, &
                              spin_expectation, this_k_point, k_point_comps, unit )
  !*  PURPOSE
  !*    Writes a comparison of spin-orbit-coupled eigenvalues and occupation
  !*    numbers at a given k-point to a user-specified unit.
  !*  USES
    use dimensions, only: n_states, n_spin, save_soc_perturbed_eigenvectors
    use constants, only: hartree
    use localorb_io, only: localorb_info
    use mpi_tasks, only: aims_stop
    use dimensions_soc, only: n_states_sr, sr_state_start, n_saved_states_soc,&
        soc_saved_state_start
    implicit none
  !*  ARGUMENTS
    real*8,  intent(in) :: KS_eigenvalue(n_states, n_spin)
    real*8,  intent(in) :: KS_eigenvalue_soc_perturbed(n_saved_states_soc)
    real*8,  intent(in) :: occ_numbers_soc(n_saved_states_soc)
    real*8,  intent(in) :: spin_expectation(n_saved_states_soc, 5)
    integer, intent(in) :: this_k_point
    real*8,  intent(in) :: k_point_comps(3)
    integer, intent(in) :: unit
  !*  INPUTS
  !*    o KS_eigenvalue               - The scalar-relativistic eigenvalues at the current k-point
  !*    o KS_eigenvalue_soc_perturbed - The spin-orbit-coupled eigenvalues at the current k-point
  !*    o occ_numbers_soc             - The spin-orbit-coupled occupation numbers at the current k-point
  !*    o spin_expection              - Expectation values for Pauli operators at current k-point
  !*    o this_k_point                - The "global" k-point index for the current k-point
  !*    o k_point_comps               - Reciprocal lattice coordinates for current k-point (no idea why
  !*                                    I chose this name for this variable, to be honest)
  !*    o unit                        - The unit to write output to.
  !*  OUTPUTS
  !*    o None (writes to screen or file)
  !*  AUTHORS
  !*    William Huhn
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    ! Counters
    integer :: i_state

    character*300                  :: info_str

    ! Used for sorting levels by energy
    real*8,  dimension(2*n_states)         :: reordered_unperturb_energies_temp
    real*8,  dimension(n_saved_states_soc) :: reordered_unperturb_energies
    real*8,  dimension(n_saved_states_soc) :: energy_spacings

    character(*), parameter :: func = 'write_soc_values'

    ! Create an ordered list of unperturbed eigenvalues including both spin
    ! channels
    call create_sorted_sr_eigenvalues_list( n_states, 1, &
         KS_eigenvalue, reordered_unperturb_energies_temp )
    reordered_unperturb_energies(1:n_saved_states_soc) = &
         reordered_unperturb_energies_temp(soc_saved_state_start:soc_saved_state_start + n_saved_states_soc - 1)

    ! Calculate the spacings between (SOC-perturbed) energy levels
    energy_spacings = 0.0d0
    if (n_saved_states_soc.gt.1) then
      do i_state = 2, n_saved_states_soc, 1
        energy_spacings(i_state) = KS_eigenvalue_soc_perturbed(i_state) - &
             KS_eigenvalue_soc_perturbed(i_state-1)
      end do
    end if

    ! And output the eigenvalues at k-points requested
    write(info_str,'(2X,A,I8,A,3(2X,F10.6), A)')   'K-point:', this_k_point, ' at',k_point_comps(1), &
      k_point_comps(2), k_point_comps(3), ' (in units of recip. lattice)'
    call localorb_info(info_str, unit)
    write(info_str, *)
    call localorb_info(info_str, unit)
!    write(info_str,'(2X,A, 4X,A, 4X,A, 4X,A, 4X,A, 4X,A, 3X,A, 2X,A, 2X,A, 2X,A)') &
!         "State", "Occupation", "Unperturbed Eigenvalue [eV]", "Eigenvalue [eV]",&
!         "Level Spacing [eV]", "spin-up", "spin-down", "<sigma_x>", "<sigma_y>", "<sigma_z>"
    write(info_str,'(2X,A, 4X,A, 4X,A, 4X,A, 4X,A)') &
         "State", "Occupation", "Unperturbed Eigenvalue [eV]", "Eigenvalue [eV]",&
         "Level Spacing [eV]"
    call localorb_info(info_str,unit)
    do i_state = 1, n_saved_states_soc, 1
    ! 7 September 2017:  I am no longer calculating or outputting these quantities, as I am unclear
    !                    how physical the resulting values are.  That being said, I would like to
    !                    return to this in the future, however, so I am not completely eliminating the code.
    !                    (The evaluation of the spin expectation valueschas not been updated to support BLACS.)
!        write(info_str,'(2X,I5, 6X,F8.5, 17X,F14.6, 5X,F14.6, 7X,F14.6, 5X,F7.4, 4X,F7.4, 4X,F7.4, 4X,F7.4, 4X,F7.4)') &
!             i_state+soc_saved_state_start-1, occ_numbers_soc(i_state), &
!             (reordered_unperturb_energies(i_state)*hartree), &
!             (KS_eigenvalue_soc_perturbed(i_state)*hartree), &
!              energy_spacings(i_state)*hartree, &
!              spin_expectation(i_state,1), &
!              spin_expectation(i_state,2), &
!              spin_expectation(i_state,3), &
!              spin_expectation(i_state,4), &
!              spin_expectation(i_state,5)
        write(info_str,'(2X,I5, 6X,F8.5, 17X,F14.6, 5X,F14.6, 7X,F14.6)') &
             i_state+soc_saved_state_start-1, occ_numbers_soc(i_state), &
             (reordered_unperturb_energies(i_state)*hartree), &
             (KS_eigenvalue_soc_perturbed(i_state)*hartree), &
              energy_spacings(i_state)*hartree

      call localorb_info(info_str,unit)
    end do
  end subroutine write_soc_values
  !******

  !****f* soc_utilities/write_soc_perturbed_eigenvectors
  !*  NAME
  !*    write_soc_perturbed_eigenvectors
  !*  SYNOPSIS
  subroutine write_soc_perturbed_eigenvectors( file_name, k_point, KS_eigenvector_soc_perturbed, &
                                               KS_eigenvalue_soc_perturbed )
  !*  PURPOSE
  !*    Writes out SOC-perturbed eigenvectors at a specified k-point to a file.
  !*  USES
    use constants, only: hartree
    use mpi_tasks, only: aims_stop
    use dimensions, only: n_basis
    use dimensions_soc, only: n_saved_states_soc, n_basis_soc, &
        n_basis_soc_coll, n_basis_soc_ncoll
    use basis, only: basis_fn, basis_atom, basisfn_type, basisfn_n, basis_l, &
        basis_m
    implicit none
  !*  ARGUMENTS
    character(len=*), intent(in) :: file_name
    integer, intent(in) :: k_point
    complex*16, intent(in) :: KS_eigenvector_soc_perturbed(n_basis_soc,n_saved_states_soc)
    real*8, intent(in) :: KS_eigenvalue_soc_perturbed(n_saved_states_soc)
  !*  INPUTS
  !*    o file_name                    - Name of file where spin-orbit-coupled eigenvector will be output
  !*    o k_point                      - The index for the k-point (only used for output)
  !*    o KS_eigenvector_soc_perturbed - The spin-orbit-coupled eigenvectors at the current k-point
  !*    o KS_eigenvalue_soc_perturbed  - The spin-orbit-coupled eigenvalues at the current k-point
  !*  OUTPUTS
  !*    o None (writes to file)
  !*  AUTHORS
  !*    William Huhn
  !*  COPYRIGHT
  !*    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !*    e.V. Please note that any use of the "FHI-aims-Software" is subject
  !*    the terms and conditions of the respective license agreement."
  !*  SOURCE

    integer :: i_basis, i_state, i_spin, basis_offset ! iterators
    integer :: i_fn
    character l_char
    character l_to_str
    character*300 :: info_str
    character*3 :: spin_str

    character(*), parameter :: func = 'write_SOC_perturbed_eigenvectors'

    ! Because we are using the basis set indexing arrays from the scalar-relativistic code,
    ! we must use the same scalar-relativistic basis set
    if (n_basis_soc .ne. 2*n_basis) then
      write(info_str,'(1X,A,A)') '* Incorrect number of SOC-perturbed basis elements specified.  Exiting.'
      call aims_stop(info_str, func)
    end if
    if (mod(n_basis_soc_coll,2) .eq. 1) then
      call aims_stop('The number of collinear basis functions in SOC is odd, exiting.', func)
    end if

    if (n_basis_soc_ncoll .gt. 0) then
      call aims_stop('This subroutine does not support usage of non-collinear basis functions&
                     & in SOC, exiting.', func)
    end if

    open (50, FILE=file_name)

    do i_state = 1, n_saved_states_soc
      write(50,'(A,I9,A,F15.8,A,F15.8,A,I9,A)') "Eigenstate ", i_state, " with energy ", KS_eigenvalue_soc_perturbed(i_state), &
           " Ha = ", KS_eigenvalue_soc_perturbed(i_state)*hartree, " eV at k-point ", k_point, " : "
      write(50, '(A)') "Basis#  Atom# type      n l   m s     coefficient"

      do i_basis = 1, n_basis_soc_coll/2, 1
        do i_spin = 1, 2, 1
          if (i_spin .eq. 1) then
            basis_offset = 0
            spin_str = " up"
          else
            basis_offset = n_basis_soc_coll/2
            spin_str = " dn"
          end if
          i_fn = basis_fn(i_basis)
          l_char = l_to_str(basis_l(i_basis))

          write(50, &
               '(2X, I5,1X,I5,1X,A8,I3,1X,A1,1X,I3,A3,3X,F19.15,A,F19.15,A)') &
               i_basis, basis_atom(i_basis), basisfn_type(i_fn), &
               basisfn_n(i_fn), l_char, basis_m(i_basis), spin_str, &
               real(KS_eigenvector_soc_perturbed(basis_offset+i_basis, i_state)), &
               " + ", aimag(KS_eigenvector_soc_perturbed(basis_offset+i_basis, i_state)),&
               " * i"
        end do
      end do

      write(50,*)
    end do

    close(50)
  end subroutine write_soc_perturbed_eigenvectors
  !******
end module soc_utilities
!******
