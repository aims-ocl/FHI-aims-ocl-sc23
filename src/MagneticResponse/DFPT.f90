!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Computes the first-order wavefunctions for some perturbation using
!!  density functional perturbation theory (DFPT) (also called
!!  coupled-perturbed Kohn-Sham, CPKS). The main procedure is do_DFPT,
!!  for which the primary input variables are the non-self-consistent
!!  part of the perturbing Hamiltonian and the first-order overlap matrix.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module DFPT

  use aims_memory_tracking,  only: aims_allocate
  use calculate_fock_matrix_p0, only: evaluate_exchange_matr_realspace_p0
  use DFPT_pulay_mixing,     only: cleanup_pulay_mixing, pulay_mix
  use dimensions,            only: n_basis, n_full_points, n_loc_prodbas, &
       & n_spin, n_states_k, use_hartree_fock, use_gga
  use hartree_fock,          only: ovlp_3fn
  use hartree_fock_p0,       only: hf_exchange_matr_real
  use integrands
  use integration,           only: integrate
  use load_balancing,        only: batch_perm, n_bp_integ, &
       & permute_point_array, permute_point_array_back
  use localorb_io,           only: localorb_multi, output_priority
  use MR_global,             only: my_max_occ_row, my_max_occ_col, &
       & n_occ_states, walltime_mul
  use tools,                 only: mul, safe_deallocate, start_wall, &
       & stop_wall, str
  use physics,               only: delta_v_hartree_deriv_l0_at_zero, &
       & delta_v_hartree_part_at_zero, en_density_embed, &
       & free_hartree_superpos, free_rho_superpos, hartree_partition_tab, &
       & KS_eigenvalue, l_hartree_max_far_distance, multipole_moments, &
       & multipole_radius_sq, outer_potential_radius, partition_tab
  use prodbas,               only: basis_nghr
  use runtime_choices,       only: dfpt_accuracy_n1, dfpt_iter_limit, &
       & dfpt_linear_mix_param, dfpt_pulay_steps, flag_xc, hybrid_coeff, &
       & output_level, RI_type, use_load_balancing
  use scalapack_wrapper,     only: eigenvec, my_row, my_col, mxld, mxcol, &
       & n_my_rows
  use synchronize_mpi_basic, only: sync_real_number, sync_vector
  use timing,                only: start_timer, stop_timer
  use types,                 only: dp

  implicit none

  private
  public :: do_DFPT, cleanup_DFPT, reset_DFPT_timers, timestamps_H1, &
       & timestamps_Ha1, timestamps_rho1, walltime_Ha1_update

  ! First-order matrix elements between eigenfunctions. In CPHF
  ! terminology, this is sometimes called the U matrix. Global
  ! dimensions: (n_basis, n_basis)
  real(dp), allocatable :: Psi0_H1_Psi0(:,:)
  ! Same as Psi0_H1_Psi0, but contains only the non-self-consistent part.
  real(dp), allocatable :: Psi0_H1_Psi0_nsc(:,:,:,:)
  ! Matrix elements of the self-consistent part of the first-order
  ! Hamiltonian between basis functions.
  ! Global dimensions: (n_basis,n_basis,n_spin,n_dirs)
  real(dp), allocatable :: H1_sc(:,:,:,:)
  ! First-order electron density (2,n_full_points,n_dirs)
  real(dp), allocatable :: rho1(:,:,:)
  ! First-order electron density gradient (3,2,n_full_points,n_dirs)
  real(dp), allocatable :: rho1_gradient(:,:,:,:)
  ! First-order Hartree potential (n_full_points,n_dirs)
  real(dp), allocatable :: hartree_pot1(:,:)
  ! H1|Psi0>, where Psi0 are the unperturbed eigenstates.
  ! dimensions: (n_basis,n_basis)
  real(dp), allocatable :: H1_Psi0(:,:)
  ! Matrix elements of the first order Hartree potential between basis
  ! functions. Global dimensions: (n_basis,n_basis,n_dirs)
  real(dp), allocatable :: H1_Ha(:,:,:)
  ! First-order density matrix. Only calculated if there is more than
  ! 1 DFPT step. Global dimensions: (n_basis,n_basis,n_spin,n_dirs)
  real(dp), allocatable, public :: density_matrix1(:,:,:,:)
  ! First-order density matrix from the previous iteration. Global
  ! dimensions: (n_basis,n_basis,n_spin,n_dirs)
  real(dp), allocatable :: density_matrix1_prev(:,:,:,:)
  ! Timers for the self-consistent part of the first-order
  ! Hamiltonian, first-order Hartree potential, and for the update of
  ! first-order density (including first-order density gradient).
  real(dp) :: timestamps_H1(4), timestamps_Ha1(4), timestamps_rho1(4)
  ! A separate timer is for updating the first-order Hartree potential
  ! since it is wrapped around routines that perform syncing across
  ! CPUs and thus individual CPU times are not meaningful.
  integer :: walltime_Ha1_update

contains

  !!  COMMENT
  !!
  !!  This is run automatically with first call to do_DFPT.
  !!
  subroutine allocate_DFPT(n_dirs, imaginary_pert, calc_H1_Ha)
    integer, intent(in) :: n_dirs
    logical, intent(in) :: imaginary_pert
    logical, intent(in) :: calc_H1_Ha ! Same as in do_DFPT
    character(*), parameter :: THIS_SUB = 'DFPT::allocate_DFPT::'
    call aims_allocate(Psi0_H1_Psi0, mxld, mxcol, THIS_SUB//'Psi0_H1_Psi0')
    call aims_allocate(Psi0_H1_Psi0_nsc, mxld, mxcol, n_spin, n_dirs, &
         & THIS_SUB//'Psi0_H1_Psi0_nsc')
    call aims_allocate(H1_sc, mxld, mxcol, n_spin, n_dirs, THIS_SUB//'H1_sc')
    call aims_allocate(H1_Psi0, mxld, mxcol, THIS_SUB//'H1_Psi0')
    ! If only a sum-over-states and not the full DFPT cycle is
    ! requested, return here. First-order quantities such as rho1 are
    ! only required for a self-consistent calculation (dfpt_iter_limit
    ! > 1) and are zero otherwise.
    if (dfpt_iter_limit == 1) return
    ! Even for dfpt_iter_limit > 1, if the perturbation is purely
    ! imaginary, there is no contribution from first-order density of
    ! the Hartree potential.
    if (.not. imaginary_pert) then
       call aims_allocate(rho1, 2, batch_perm(n_bp_integ)%n_full_points, &
            & n_dirs, THIS_SUB//'rho1')
       rho1 = 0d0
       call aims_allocate(hartree_pot1, batch_perm(n_bp_integ)%n_full_points, &
            & n_dirs, THIS_SUB//'hartree_pot1')
       if (use_gga) call aims_allocate(rho1_gradient, 3, 2, &
            & batch_perm(n_bp_integ)%n_full_points, n_dirs, &
            & THIS_SUB//'rho1_gradient')
    end if
    ! If the perturbation is imaginary but a portion of exact exchange
    ! is included, the first-order density matrix is computed
    ! self-consistently even though the first-order density vanishes.
    if (.not. imaginary_pert .or. use_hartree_fock) then
       call aims_allocate(density_matrix1, mxld, mxcol, n_spin, n_dirs, &
            & THIS_SUB//'density_matrix1')
       call aims_allocate(density_matrix1_prev, mxld, mxcol, n_spin, n_dirs, &
            & THIS_SUB//'density_matrix1_prev')
    end if
    ! Allocate H1_Ha only if first-order Hartree potential was
    ! explicitly requested.
    if (calc_H1_Ha) then
       call aims_allocate(H1_Ha, mxld, mxcol, n_dirs, THIS_SUB//'H1_Ha')
       H1_Ha = 0d0
    end if
  end subroutine allocate_DFPT

  !!  FUNCTION
  !!
  !!  We cannot reset timers in allocate_DFPT because that subroutine
  !!  could be called repeatedly.
  !!
  subroutine reset_DFPT_timers()
    timestamps_H1 = 0d0
    timestamps_Ha1 = 0d0
    timestamps_rho1 = 0d0
    walltime_Ha1_update = 0d0
  end subroutine reset_DFPT_timers

  !!  FUNCTION
  !!
  !!  Performs a DFPT cycle. The first call to this procedure
  !!  automatically allocates the work variables. Afterwards, it is
  !!  the caller's responsibility to call cleanup_DFPT().
  !!
  !!  The first-order wavefunction coefficient are calculated as
  !!  follows:
  !!
  !!  for i -> occupied, j -> virtual:
  !!    C1 = C0*U,
  !!    U_ij = U0_ij/(epsilon_j-epsilon_i)
  !!    U0 = C0*H1*C0 - C0*S1*C0*epsilon,
  !!  for i -> occupied, j -> occupied:
  !!    C1 = C0*U,
  !!    U_ij = -1/2*C0*S1*C0,
  !!
  !!  where S1 - overlap1, H1 - H1_nsc+H1_sc, U - Psi0_H1_Psi0.
  !!
  subroutine do_DFPT(H1_nsc, overlap1, imaginary_pert, n_dirs, calc_H1_Ha, &
       & include_overlap1, spin_operator, Psi1)
    ! Matrix elements between the basis functions of the
    ! non-self-consistent (no dependence on rho1) part of the
    ! first-order Hamiltonian. Global dimensions:
    ! (n_basis,n_basis,n_dirs*n_spin). See n_dirs below.
    real(dp), intent(in) :: H1_nsc(:,:,:)
    ! First-order overlap matrix. If include_overlap1=.false., this
    ! may be left unallocated.
    real(dp), intent(in out) :: overlap1(:,:,:)
    ! Whether the perturbation is purely imaginary. If so, there is no
    ! first-order density and, in case of LDA/GGA, the DFPT cycle
    ! reduces to a single step. With HF, there is still no first-order
    ! density but the first-order density matrix needs to be computed
    ! self-consistently. With a hybrid functional, the LDA/GGA part is
    ! ignored during the self-consistency cycle and only the
    ! contribution from exact exchanged is calculated.
    logical, intent(in) :: imaginary_pert
    ! Number of perturbation directions simultaneously
    ! processed. Current implementation allows max 3.
    integer, intent(in) :: n_dirs
    ! Whether to calculate the first-order Hartree potential
    logical, intent(in) :: calc_H1_Ha
    ! Whether to include the first-order overlap matrix in the
    ! calculations.
    logical, intent(in) :: include_overlap1
    ! This is relevant if n_spin==1 and the perturbation contains the
    ! Sz operator. In that case both spin channels are explicitly
    ! required when computing the first-order xc response, even though
    ! the zero-order density is non-spin-polarized. It is applied by
    ! manually setting rho1(2,:,:) = -rho1(1,:,:) and likewise for the
    ! first-order density gradient.
    logical, intent(in) :: spin_operator
    ! First-order wavefunctions. Global dimensions:
    ! (n_basis,<n_occ_states>,2,n_dirs)
    real(dp), intent(out) :: Psi1(:,:,:,:)
    ! RMT change in the first-order density matrix
    real(dp) :: density_matrix1_change
    ! Other
    real(dp) :: beta
    character(128) :: info_str
    integer :: i_spin, iter, i_row, i_col, i_dir, i_min

    if (.not. allocated(Psi0_H1_Psi0)) &
         & call allocate_DFPT(n_dirs, imaginary_pert, calc_H1_Ha)
    main: do iter = 1, dfpt_iter_limit
       ! STEP 1 - In the first iteration, find the non-self-consistent
       !          parts of the U_ij coefficients.
       H1_sc = 0d0
       if (iter == 1) then
          do i_dir = 1, n_dirs
             do i_spin = 1, n_spin
                ! -C0*S1*C0
                if (include_overlap1) then
                   call start_wall(walltime_mul)
                   call mul(overlap1(:,:,i_dir), eigenvec(:,:,i_spin), &
                        & H1_psi0, N=n_occ_states(i_spin))
                   call mul(eigenvec(:,:,i_spin), H1_psi0, &
                        & Psi0_H1_Psi0_nsc(:,:,i_spin,i_dir), transa='t', &
                        & alpha=-1d0)
                   call stop_wall(walltime_mul)
                   ! U0 = -1/2*C0*S1*C0
                   i_row = my_max_occ_row(i_spin)
                   Psi0_H1_Psi0_nsc(:i_row,:,i_spin,i_dir) = &
                        & Psi0_H1_Psi0_nsc(:i_row,:,i_spin,i_dir)/2
                   ! U0 = -C0*S1*C0*epsilon
                   do i_row = my_max_occ_row(i_spin)+1, n_my_rows
                      Psi0_H1_Psi0_nsc(i_row,:,i_spin,i_dir) = &
                           & Psi0_H1_Psi0_nsc(i_row,:,i_spin,i_dir)* &
                           & KS_eigenvalue(my_col,i_spin,1)
                   end do
                else
                   Psi0_H1_Psi0_nsc(:,:,i_spin,i_dir) = 0d0
                end if
                ! U0 += C0*H1_nsc*C0
                call start_wall(walltime_mul)
                call mul(H1_nsc(:,:,(i_dir-1)*n_spin+i_spin), &
                     & eigenvec(:,:,i_spin), H1_Psi0, N=n_occ_states(i_spin))
                call mul(eigenvec(:,:,i_spin), H1_psi0, &
                     & Psi0_H1_Psi0_nsc(:,:,i_spin,i_dir), transa='t', &
                     & M=n_states_k(1)-n_occ_states(i_spin), &
                     & N=n_occ_states(i_spin), jA=n_occ_states(i_spin)+1, &
                     & iC=n_occ_states(i_spin)+1, beta=1d0)
                call stop_wall(walltime_mul)
             end do
          end do
          ! STEP 2 - In subsequent iterations, compute the matrix
          !          elements of the self-consistent part of the
          !          first-order Hamiltonian (H1_sc). This only takes
          !          place if the perturbation is not purely imaginary
          !          or if we include a part of exact exchange (else
          !          there is only one iteration).
       else
          ! Unless we are doing pure Hartree-Fock, calculate the
          ! first-order response of the xc-kernel.
          if (.not. imaginary_pert .and. flag_xc /= 0) then
             ! The output is not a density matrix. We use
             ! density_matrix1_aux here in order to save memory from
             ! allocating another real array.
             INT_XC_KERNEL%n_dimensions = n_dirs
             call integrate(INT_XC_KERNEL, H1_sc, timer=timestamps_H1, &
                  & rho1=rho1, rho1_gradient=rho1_gradient)
          end if
          ! With HF or hybrids, calculate the first-order response of
          ! exact exchange.
          if (use_hartree_fock) call add_exact_exchange(density_matrix1, H1_sc)
          ! If required, calculate the first-order response of the
          ! Hartree potential.
          if (calc_H1_Ha) then
             call start_wall(walltime_Ha1_update)
             do i_dir = 1, n_dirs
                call update_Ha1(rho1(:,:,i_dir), hartree_pot1(:,i_dir))
             end do
             call stop_wall(walltime_Ha1_update)
             INT_FIRST_ORDER_HARTREE%n_dimensions = n_dirs
             call integrate(INT_FIRST_ORDER_HARTREE, H1_Ha, &
                  & timer=timestamps_Ha1, hartree_pot1=hartree_pot1)
             H1_sc = H1_sc + spread(H1_Ha,3,2)
          end if
       end if
       ! STEP 3 - Find the self-consistent parts of the U_ij
       !          coefficients.
       do i_dir = 1, n_dirs
          do i_spin = 1, n_spin
             Psi0_H1_Psi0 = Psi0_H1_Psi0_nsc(:,:,i_spin,i_dir)
             ! U0 += C0*H1_sc*C0 (for occupied-virtual only)
             if (iter > 1) then
                call start_wall(walltime_mul)
                call mul(H1_sc(:,:,i_spin,i_dir), eigenvec(:,:,i_spin), &
                     & H1_Psi0, N=n_occ_states(i_spin))
                call mul(eigenvec(:,:,i_spin), H1_psi0, Psi0_H1_Psi0, &
                     & transa='t', M=n_states_k(1)-n_occ_states(i_spin), &
                     & N=n_occ_states(i_spin), jA=n_occ_states(i_spin)+1, &
                     & iC=n_occ_states(i_spin)+1, beta=1d0)
                call stop_wall(walltime_mul)
             end if
             ! U_ij = U0_ij/(epsilon_j-epsilon_i)
             i_min = my_max_occ_row(i_spin)+1
             do i_col = 1, my_max_occ_col(i_spin)
                ! First index is over unoccupied states
                Psi0_H1_Psi0(i_min:,i_col) = Psi0_H1_Psi0(i_min:,i_col)/ &
                     & (KS_eigenvalue(my_col(i_col),i_spin,1) - &
                     & KS_eigenvalue(my_row(i_min:),i_spin,1))
             end do
             ! C1 = C0*U
             call start_wall(walltime_mul)
             call mul(eigenvec(:,:,i_spin), Psi0_H1_Psi0, &
                  & Psi1(:,:,i_spin,i_dir), N=n_states_k(1))
             call stop_wall(walltime_mul)
          end do
       end do
       ! We have the first-order wavefunctions. If this is a
       ! single-shot calculation, exit now.
       if (dfpt_iter_limit == 1 .or. &
            & imaginary_pert .and. .not. use_hartree_fock) exit main
       ! STEP 4 - Construct a new first-order density matrix from the
       !          first-order wavefunctions.
       if (imaginary_pert) then
          ! The expression for first-order density matrix is
          !   n1_ij = Sum_m^occ C_im C1_jm^* + C1_im C_jm.
          ! However, all arrays are real and thus complex conjugation
          ! has no effect. Instead, if Psi1 are imaginary, we compute
          ! n1 from
          !   n1_ij = Sum_m^occ (-C_im C1_jm + C1_im C_jm).
          beta = -1d0
       else
          beta = 1d0
       end if
       do i_dir = 1, n_dirs
          do i_spin = 1, n_spin
             call start_wall(walltime_mul)
             call mul(eigenvec(:,:,i_spin), Psi1(:,:,i_spin,i_dir), &
                  & density_matrix1(:,:,i_spin,i_dir), transb='t', &
                  & K=n_occ_states(i_spin))
             ! Here beta=-1 corresponds to the minus sign in the above
             ! formula for n1.
             call mul(Psi1(:,:,i_spin,i_dir), eigenvec(:,:,i_spin), &
                  & density_matrix1(:,:,i_spin,i_dir), transb='t', &
                  & beta=beta, K=n_occ_states(i_spin))
             call stop_wall(walltime_mul)
          end do
       end do
       ! STEP 5 - Check convergence and exit if converged
       ! Change in n1^2 over all directions and spin channels
       if (iter == 1) density_matrix1_prev = 0d0
       density_matrix1_change = (3-n_spin)* &
            & sum(abs(density_matrix1 - density_matrix1_prev)**2)
       call sync_real_number(density_matrix1_change)
       density_matrix1_change = sqrt(density_matrix1_change)
       if (iter == 1) call localorb_multi('iter   change in first-order &
            &density matrix', format='(6x, a)')
       if (density_matrix1_change > dfpt_accuracy_n1) then
          write(info_str, '(i9, 3x, 1(es14.6e2))') &
               & iter, real(density_matrix1_change,4)
       else
          write(info_str, '(i9, 3x, 1(es14.6e2), a, es9.2e2, a)') &
               & iter, real(density_matrix1_change,4), ' <', dfpt_accuracy_n1, &
               & ' (dfpt_accuracy_n1)'
          call localorb_multi(info_str)
          exit main
       end if
       call localorb_multi(info_str)
       if (iter == dfpt_iter_limit) then
          call localorb_multi('', &
               & 'WARNING: DFPT cycle did not converge in '// &
               & str(dfpt_iter_limit)//' steps (dfpt_iter_limit).', &
               & '', format='(2x, a)')
          exit main
       end if
       ! STEP 6 - Perform density matrix mixing
       call pulay_mix(density_matrix1, size(density_matrix1), iter, &
            & dfpt_pulay_steps, dfpt_linear_mix_param)
       density_matrix1_prev = density_matrix1
       ! STEP 7 - Construct a new density from the density matrix.
       if (.not. imaginary_pert .and. (flag_xc /= 0 .or. calc_H1_Ha)) then
          INT_ELECTRON_DENSITY%n_dimensions = n_dirs
          call integrate(INT_ELECTRON_DENSITY, timer=timestamps_rho1, &
               & rho1=rho1, rho1_gradient=rho1_gradient, &
               & density_matrix=density_matrix1)
       end if
       if (n_spin == 1 .and. spin_operator) then
          rho1(2,:,:) = -rho1(1,:,:)
          if (use_gga) rho1_gradient(:,2,:,:) = -rho1_gradient(:,1,:,:)
       end if
    end do main
    ! Call this after each DFPT cycle because the dimensions of the
    ! pulay work arrays could be different for the next perturbation.
    call cleanup_pulay_mixing()
  end subroutine do_DFPT

  !!  FUNCTION
  !!
  !!  Given a first-order density, calculates the corresponding
  !!  first-order Hartree potential.
  !!
  subroutine update_Ha1(rho1, hartree_pot1)
    real(dp), intent(in) :: rho1(:,:)
    real(dp), intent(out) :: hartree_pot1(:)
    character(*), parameter :: THIS_SUB = 'DFPT::update_Ha1::'
    real(dp), allocatable :: rho1_p(:,:), hartree_pot1_tmp(:)
    real(dp) :: en_elec_delta, hartree_delta_energy, &
         & hartree_multipole_correction
    integer :: output_priority_save
    character*20 :: output_level_save
    call aims_allocate(rho1_p, 2, n_full_points, THIS_SUB//'rho1_p')
    call aims_allocate(hartree_pot1_tmp, n_full_points, &
         & THIS_SUB//'hartree_pot1_tmp')
    if (use_load_balancing) then
       call permute_point_array_back(n_bp_integ, 2, rho1, rho1_p)
    else
       rho1_p = rho1
    end if
    ! Suppress output from these subroutines
    output_priority_save = output_priority
    output_priority = 2
    output_level_save = output_level
    output_level = ''
    call update_hartree_potential_p1(hartree_partition_tab, free_rho_superpos, &
         & rho1_p, 0d0, delta_v_hartree_part_at_zero, &
         & delta_v_hartree_deriv_l0_at_zero, multipole_moments, &
         & multipole_radius_sq, l_hartree_max_far_distance, 0d0, &
         & outer_potential_radius, .false., .false.)
    call sum_up_whole_potential_p1(delta_v_hartree_part_at_zero, &
         & delta_v_hartree_deriv_l0_at_zero, multipole_moments, partition_tab, &
         & rho1_p, free_hartree_superpos, free_rho_superpos, &
         & hartree_pot1_tmp, hartree_delta_energy, en_elec_delta, &
         & hartree_multipole_correction, 0d0, en_density_embed, 0d0, .false., &
         & multipole_radius_sq, l_hartree_max_far_distance, 0d0, 0d0, 0d0, &
         & outer_potential_radius, .false., .false.)
    output_priority = output_priority_save
    output_level = output_level_save
    if (use_load_balancing) then
       call permute_point_array(n_bp_integ, 1, hartree_pot1_tmp, hartree_pot1)
    else
       hartree_pot1 = hartree_pot1_tmp
    end if
    call safe_deallocate(rho1_p)
    call safe_deallocate(hartree_pot1_tmp)
  end subroutine update_Ha1

  !!  FUNCTION
  !!
  !!  Construct the exact exchange matrix from an input density
  !!  matrix. The exchange matrix is added to exchange_matrix.
  !!
  subroutine add_exact_exchange(density_matrix, exchange_matrix)
    real(dp), intent(in) :: density_matrix(:,:,:,:)
    real(dp), intent(in out) :: exchange_matrix(:,:,:,:)
    ! Work arrays for RI V
    real(dp), allocatable :: density_matrix_full(:,:,:)
    real(dp), allocatable :: exchange_matrix_full(:,:,:)
    real(dp), allocatable :: ovlp_loc(:,:)
    real(dp), allocatable :: tmp(:,:)
    character(*), parameter :: THIS_SUB = 'DFPT::add_exact_exchange::'
    real(dp)  :: dummy
    integer :: output_priority_save, i_index, mb, nb, i_spin, i_dir
    if (RI_type == 2) then ! RI V
       call aims_allocate(density_matrix_full, n_basis, n_basis, n_spin, &
            & THIS_SUB//'density_matrix_full')
       call aims_allocate(exchange_matrix_full, n_basis, n_basis, n_spin, &
            & THIS_SUB//'exchange_matrix_full')
       call aims_allocate(ovlp_loc, n_basis, n_basis, THIS_SUB//'ovlp_loc')
       call aims_allocate(tmp, n_basis, n_basis, THIS_SUB//'tmp')
       do i_dir = 1, size(density_matrix,4)
          density_matrix_full = 0d0
          density_matrix_full(my_row,my_col,:) = density_matrix(:,:,:,i_dir)
          call sync_vector(density_matrix_full, size(density_matrix_full))
          exchange_matrix_full = 0d0
          do i_index = 1, n_loc_prodbas
             ovlp_loc = 0d0
             do mb = 1, n_basis
                do nb = 1, n_basis
                   ovlp_loc(nb,mb) = ovlp_3fn(basis_nghr(nb,mb),i_index)
                end do
             end do
             do i_spin = 1, n_spin
                call mul(ovlp_loc, density_matrix_full(:,:,i_spin), tmp, &
                     & do_serial=.true.)
                call mul(tmp, ovlp_loc, exchange_matrix_full(:,:,i_spin), &
                     & alpha=-1d0, beta=1d0, do_serial=.true.)
             end do
          end do
          call sync_vector(exchange_matrix_full,size(exchange_matrix_full))
          exchange_matrix(:,:,:,i_dir) = exchange_matrix(:,:,:,i_dir) + &
               & hybrid_coeff*exchange_matrix_full(my_row,my_col,:)
       end do
       call safe_deallocate(density_matrix_full)
       call safe_deallocate(exchange_matrix_full)
       call safe_deallocate(ovlp_loc)
       call safe_deallocate(tmp)
    else ! LVL
       output_priority_save = output_priority
       output_priority = 2
       do i_dir = 1, size(density_matrix,4)
          call evaluate_exchange_matr_realspace_p0([0d0], [(0d0,0d0)], [0d0], &
               & dummy, [dummy], .false., .false., density_matrix(:,:,:,1))
          exchange_matrix(:,:,:,i_dir) = exchange_matrix(:,:,:,i_dir) - &
               & hybrid_coeff*2*hf_exchange_matr_real(:,:,1,:)
       end do
       output_priority = output_priority_save
    end if
  end subroutine add_exact_exchange

  subroutine cleanup_DFPT()
    call safe_deallocate(Psi0_H1_Psi0)
    call safe_deallocate(Psi0_H1_Psi0_nsc)
    call safe_deallocate(H1_sc)
    call safe_deallocate(H1_Psi0)
    call safe_deallocate(H1_Ha)
    call safe_deallocate(rho1)
    call safe_deallocate(rho1_gradient)
    call safe_deallocate(hartree_pot1)
    call safe_deallocate(density_matrix1)
    call safe_deallocate(density_matrix1_prev)
  end subroutine cleanup_DFPT
end module DFPT
