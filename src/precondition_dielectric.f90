!****s* FHI-aims/precondition_dielectric
!  NAME
!    precondition_dielectric
!  SYNOPSIS
subroutine precondition_dielectric(R_in, partition_tab, hartree_partition_tab)

  !  PURPOSE
  !    Charge density preconditioner using the dielectric function
  !    This routine takes an input residual R_in and preconditions it to yield some output
  !    R_prec by solving the equation
  !
  !       R_in(r) = \int dr' \epsilon(r, r'; 0) R_prec(r'),
  !
  !    where \epsilon is the dielectric function. The solution if found via gradient descent method.
  !  

  ! USES
  use dimensions
  use grids
  use geometry
  use basis
  use runtime_choices
  use mpi_utilities
  use synchronize_mpi
  use constants
  use localorb_io
  use debugmanager
  use timing

  implicit none

  ! ARGUMENTS

  real*8, dimension(n_full_points), intent(inout) :: R_in
  real*8, dimension(n_full_points), intent(in)    :: partition_tab
  real*8, dimension(n_full_points), intent(in)    :: hartree_partition_tab

  !  INPUTS
  !   o  R_in (input) -- The charge density residual over the entire grid
  !   o  partition_tab -- values of partition function
  !   o  hartree_partition_tab -- to be used by hartree_solver
  !
  !  OUTPUT
  !   o  R_in (output) -- Preconditioned density residual
  !
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Created in November 2017
  !  SOURCE

  ! Local variables
  real*8, dimension(n_full_points) :: V_hartree
  real*8, dimension(n_full_points) :: R_prec
  real*8, dimension(n_full_points) :: chi
  
  real*8, dimension(n_full_points) :: eps

  real*8, dimension(n_full_points) :: residual, p, Ap
  real*8 :: residual_norm, rsold, alpha, rsnew

  integer :: i 

  logical :: cg_con
  real*8 :: eps2(n_full_points)
  real*8 :: eps_approx(n_full_points)
  real*8 :: sync_rsold

  real*8 :: c, d, a, b, gr
  integer :: h

  character*100 :: info_str
  logical :: out_times = .false.
  
  ! Timing infrastructure
  real*8 :: hartree_time, clock_hartree_time, chi_time, clock_chi_time
  character(*), parameter :: deffmt = '2X'
  
  ! Begin work
  write(info_str,'(2X,A)') "Using dielectric preconditioning."
  call localorb_info( info_str, use_unit,'(A)', OL_norm )

  V_hartree = 0.d0

  ! Initial guess for the gradient descent method
  R_prec = R_in
  
  ! Get hartree potential from solver
  call get_timestamps (hartree_time, clock_hartree_time)
  call hartree_solver(hartree_partition_tab, R_in, V_hartree)
  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to operate with Hartree potential",&
       rtime-hartree_time, clock_rtime-clock_hartree_time, OL_norm)

  call get_timestamps (chi_time, clock_chi_time) 
  call chi_times_v_hartree(V_hartree, chi, partition_tab)
  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to operate with response function",&
       rtime-chi_time, clock_rtime-clock_chi_time, OL_norm)

  eps = R_prec - chi 
  
  ! Approximate eps to speed up calculations
  eps_approx = eps
  
  ! Gradient descent, find the preconditioned delta rho
  do i = 1, 21, 1
    cg_con = .false.
    residual = R_in - eps 
    p = residual
    rsold = sum(residual ** 2)
    
    sync_rsold = rsold
    call sync_real_number(sync_rsold)
    if (sqrt(sync_rsold) < 1d-5) cg_con = .true.
    ! We can't exit before all threads are converged.
    if (cg_con) exit

    call get_timestamps (hartree_time, clock_hartree_time)
    call hartree_solver(hartree_partition_tab, p, V_hartree)
    call get_timestamps (rtime, clock_rtime)
    if (out_times) call output_times(deffmt, "Time to operate with Hartree potential",&
         rtime-hartree_time, clock_rtime-clock_hartree_time, OL_norm)

    call get_timestamps (chi_time, clock_chi_time) 
    call chi_times_v_hartree(V_hartree, chi, partition_tab)
    call get_timestamps (rtime, clock_rtime)
    if (out_times) call output_times(deffmt, "Time to operate with response function",&
         rtime-chi_time, clock_rtime-clock_chi_time, OL_norm)

    eps2 = p - chi 
    Ap = eps2
    alpha = rsold / sum(p * Ap)
    !!! Golden section search
    ! This may slightly improve convergence in some specific cases, 
    ! like periodic systems 
    gr = (sqrt(5.0d0) + 1.0d0) / 2.0d0
    a = -1.5d0 * alpha
    b = 1.5d0 * alpha
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    do h = 1, 11, 1
      if (abs(c - d) < 1e-9) exit
      if (sum((residual - c * Ap)**2) < sum((residual - d * Ap)**2)) then
        b = d
      else
        a = c
      end if
      c = b - (b - a) / gr
      d = a + (b - a) / gr
    end do
    
    alpha = (b + a) / 2.0d0
    !!! /Gss
    eps_approx = eps_approx + eps2 * alpha

    R_prec = R_prec + alpha * p
    residual = residual - alpha * Ap
    rsnew = sum(residual ** 2)

    ! Update epsilon and start again
    ! Epsilon is a linear operator, no need to calculate it in full every step
    eps = eps_approx
  end do

  write(info_str,'(2X,A,I15)')   " | GD iterations: ", i
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  ! Output debug info
  residual = R_in - eps 
  residual_norm = sum(residual**2)
  call sync_real_number(residual_norm)
  write(info_str,'(2X,A,E16.9)') " | Norm of residual: ", sqrt(residual_norm)
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  ! Replace input delta rho with preconditioned rho
  R_in = R_prec
end subroutine precondition_dielectric
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/chi_times_v_hartree
!  NAME
!    chi_times_v_hartree
!  SYNOPSIS

subroutine chi_times_v_hartree(V_hartree, chi, partition_tab)

  !  PURPOSE
  !    Calculates int dr'' chi_0(r, r''; 0) int dr' delta_rho(r') / |r - r'|,
  !    where chi_0 is the static response function
  !  USES

  use dimensions
  use basis
  use species_data, only: l_shell_max
  use grids
  use physics, only: KS_eigenvalue, KS_eigenvector, occ_numbers, KS_eigenvector_complex
  use localorb_io
  use debugmanager
  use runtime_choices
  use mpi_utilities
  use synchronize_mpi
  use timing

  implicit none

  ! ARGUMENTS
  
  real*8, intent(in) :: V_hartree(n_full_points)
  real*8, intent(in) :: partition_tab(n_full_points)
  real*8, intent(out) :: chi(n_full_points)

  !  INPUTS
  !    o V_hartree -- Hartree potential from hartree_solver
  !    o partition_tab -- values of the partition function for integration
  !  OUTPUTS
  !    o chi -- output contaning the product of the static response function and the Hartree potential
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Created in November 2017
  !  SOURCE

  ! Local variables
  integer :: l_ylm_max
  integer, dimension(:, :), allocatable :: index_lm
  real*8, dimension(:, :, :), allocatable :: ylm_tab
  integer :: i_l, i_m

  integer :: i_full_points
  integer :: n_compute
  integer :: i_basis(n_basis)
  integer :: n_points

  real*8 :: coord_current(3)
  real*8 :: dist_tab(n_atoms, n_max_batch_size)
  real*8 :: dist_tab_sq(n_atoms, n_max_batch_size)
  real*8 :: dir_tab(3, n_atoms, n_max_batch_size)
  real*8 :: i_r(n_atoms, n_max_batch_size)
  real*8 :: trigonom_tab(4, n_atoms, n_max_batch_size)
  real*8 :: radial_wave(n_basis, n_max_batch_size)
  real*8 :: wave(n_basis, n_max_batch_size)

  real*8 :: all_waves(n_basis, n_full_points)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns * n_atoms)
  integer :: i_basis_fns_inv(n_basis_fns, n_atoms)
  integer :: i_atom_fns(n_basis_fns * n_atoms)
  
  integer :: n_compute_atoms
  integer :: atom_index(n_atoms)
  integer :: atom_index_inv(n_atoms)

  integer :: spline_array_start(n_atoms)
  integer :: spline_array_end(n_atoms)

  integer :: i_my_batch
  integer :: i_index

  integer :: j_basis, k_basis
  integer :: i_point

  real*8, allocatable :: integral(:, :)
  integer :: i_states
  integer :: n_states_occupied
  integer :: occ, un_occ

  real*8 :: phi(n_basis, n_max_batch_size)
  real*8 :: phi_tilde(n_basis, n_max_batch_size)
  real*8 :: aux_matrix(n_basis, n_basis)
  real*8 :: phi_t_x_phi(n_basis, n_basis)
  integer :: i_full_points_2
  integer :: i, j
  real*8 :: aux_vector1(n_basis)
  real*8 :: aux_vector2(n_basis)
  real*8 :: aux_vector3(n_basis)

  real*8 :: aux_matrix2(n_basis, n_basis)

  real*8, external :: ddot
  real*8 :: temp_result
  real*8 :: sub_time, clock_sub_time
  character(*), parameter :: deffmt = '2X'
  logical :: out_times = .false.
  
  ! begin work

  ! Initialize and allocate
  chi = 0.d0
  all_waves = 0d0
  phi = 0.d0
  phi_tilde = 0.d0
  phi_t_x_phi = 0.d0

  n_states_occupied = 0

  l_ylm_max = l_wave_max

  allocate(ylm_tab( (l_ylm_max + 1) ** 2, n_atoms, n_max_batch_size))
  allocate(index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max))

  ! Initialize index_lm
  i_index = 0
  do i_l = 0, l_ylm_max, 1
    do i_m = -i_l, i_l, 1
      i_index = i_index + 1
      index_lm(i_m, i_l) = i_index
    end do
  end do

  i_full_points = 0
  i_full_points_2 = 0

  call get_timestamps (sub_time, clock_sub_time)

  ! First, we need (n_basis, n_full_points) array containing the values of all basis functions at all 
  ! integration points.
  do i_my_batch = 1, n_my_batches, 1
    n_compute = 0
    i_basis = 0
    i_point = 0
    do i_index = 1, batches(i_my_batch)%size, 1
      i_full_points_2 = i_full_points_2 + 1
      if (partition_tab(i_full_points_2) .gt. 0d0 .or. .true.) then
        i_point = i_point + 1
        coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
        
        if (n_periodic > 0) then
          call map_to_center_cell(coord_current(1:3))
        end if

        ! Compute atom-centered coordinates of current integration point as viewd from all atoms
        call tab_atom_centered_coords_v2( &
          coord_current, dist_tab_sq(1, i_point), &
          dir_tab(1, 1, i_point))

        ! Determine which basis functions u(r)/r * Y_lm(theta, phi) are actually needed
        call prune_basis_v2(dist_tab_sq(1, i_point), n_compute, i_basis)
      end if
    end do ! i_index

    n_points = i_point
    if (n_compute .gt. 0) then
      i_point = 0
      do i_index = 1, batches(i_my_batch)%size, 1
        i_full_points = i_full_points + 1
        
        ! Can't ignore points with zero partition!  

          i_point = i_point + 1
          n_compute_atoms = 0
          n_compute_fns = 0
          i_basis_fns_inv = 0

          call prune_radial_basis_v2( &
            dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
            atom_index_inv, n_compute_fns, i_basis_fns, i_basis_fns_inv, &
            i_atom_fns, spline_array_start, spline_array_end )

          call tab_local_geometry( &
            dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
            dir_tab(1, 1, i_point), dist_tab(1, i_point), i_r(1, i_point))

          call evaluate_radial_functions_v2( &
            spline_array_start, spline_array_end, n_compute_atoms, &
            n_compute_fns, dist_tab(1, i_point), i_r(1, i_point), &
            atom_index, i_basis_fns_inv, basis_wave_ordered, &
            radial_wave(1, i_point), .false. )

          call tab_trigonom_v2( &
            n_compute_atoms, dir_tab(1, 1, i_point), trigonom_tab(1, 1, i_point) )

          call tab_wave_ylm_v2( &
            n_compute_atoms, atom_index, trigonom_tab(1, 1, i_point), &
            l_shell_max, l_ylm_max, ylm_tab(1, 1, i_point))

          call evaluate_waves_v2( &
            l_ylm_max, ylm_tab(1, 1, i_point), dist_tab(1, i_point), &
            index_lm, n_compute, i_basis, radial_wave(1, i_point), &
            wave(1, i_point), n_compute_atoms, atom_index_inv, i_basis_fns_inv )

          ! reset i_basis_fns_inv
          i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

          phi(:, i_point) = wave(:, i_point)

          ! phi_tilde contains V_hartree and partition_tab
          phi_tilde(:, i_point) = wave(:, i_point) * V_hartree(i_full_points) * &
            partition_tab(i_full_points)

          do j_basis = 1, n_compute, 1
            all_waves(i_basis(j_basis), i_full_points) = wave(j_basis, i_point)
          end do
          

      end do ! i_index
      
      ! Compute the phi_tilde x phi matrix
      call dgemm('N', 'T', n_compute, n_compute, n_points, 1.0d0, &
        phi_tilde, n_basis, phi, n_basis, 0.0d0, aux_matrix, n_basis)
      do i = 1, n_compute, 1
        do j = 1, n_compute, 1
          phi_t_x_phi(i_basis(i), i_basis(j)) = phi_t_x_phi(i_basis(i), i_basis(j)) + aux_matrix(i, j)
        end do
      end do

    else
      i_full_points = i_full_points + batches(i_my_batch)%size
    end if ! n_compute > 0
  end do ! i_my_batch

  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to loop over batches",&
       rtime-sub_time, clock_rtime-clock_sub_time, OL_norm)

  ! Second, we want to calculate the integral Psi_occ * Psi_unocc * V_hartree
  do i_states = n_states, 1, -1
    n_states_occupied = i_states
    if (occ_numbers(i_states, 1, 1) > 1.d-5) then
      exit
    end if
  end do
  
  allocate(integral(1:n_states_occupied,  n_states_occupied+1:n_states))
!!$  integral = 0.0d0

  call get_timestamps (sub_time, clock_sub_time)
  
!!$  do occ = 1, n_states_occupied, 1
!!$    do un_occ = n_states_occupied + 1, n_states, 1
!!$      
!!$      ! Calculate the integral as a matrix product
!!$      if (real_eigenvectors) then
!!$        aux_vector1 = KS_eigenvector(:, un_occ, 1, 1)
!!$        aux_vector3 = KS_eigenvector(:, occ, 1, 1)
!!$      else
!!$        aux_vector1 = dble(KS_eigenvector_complex(:, un_occ, 1, 1))
!!$        aux_vector3 = dble(KS_eigenvector_complex(:, occ, 1, 1))
!!$      end if
!!$      
!!$      call dgemv('N', n_basis, n_basis, 1.0d0, phi_t_x_phi, n_basis, aux_vector1, 1, 0.d0, aux_vector2, 1)
!!$      integral(occ, un_occ) = dot_product(aux_vector3, aux_vector2)
!!$
!!$    end do ! un_occ
!!$  end do ! occ

  call dgemm('N', 'N', n_basis, n_states - n_states_occupied, n_basis, 1.0d0, &
       phi_t_x_phi, n_basis, KS_eigenvector(:,n_states_occupied+1:n_states,1,1), n_basis, &
       0.0d0, aux_matrix, n_basis)
  
  call dgemm('T', 'N', n_states_occupied, n_states - n_states_occupied, n_basis, 1.0d0, &
       KS_eigenvector(:,1:n_states_occupied,1,1), n_basis, aux_matrix, n_basis, &
       0.0d0, integral(1,n_states_occupied+1), n_states_occupied )
  
  
  ! Integral needs to be synced here if we are running in parallel!
  if (use_mpi) then
    call sync_vector(integral, n_states_occupied * (n_states - n_states_occupied))
  end if

  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to calculate 'integral' array",&
       rtime-sub_time, clock_rtime-clock_sub_time, OL_norm)

  ! Third, we want the full result
  aux_matrix2 = 0.0d0

  call get_timestamps (sub_time, clock_sub_time)
  
  ! 3.1 Calculate matrix R, where
  ! R_ij = Sum_occ Sum_unocc integral(occ, unocc) * c_{occ,i} * c_{unocc,j} / (E_occ - E_unocc)
!!$  do j_basis = 1, n_basis, 1
!!$    do k_basis = 1, n_basis, 1
!!$      temp_result = 0.d0
!!$      do occ = 1, n_states_occupied, 1
!!$        do un_occ = n_states_occupied + 1, n_states, 1
!!$          if (real_eigenvectors) then
!!$            aux_vector1 = KS_eigenvector(:, un_occ, 1, 1)
!!$            aux_vector3 = KS_eigenvector(:, occ, 1, 1)
!!$          else
!!$            aux_vector1 = dble(KS_eigenvector_complex(:, un_occ, 1, 1))
!!$            aux_vector3 = dble(KS_eigenvector_complex(:, occ, 1, 1))
!!$          end if
!!$          temp_result = temp_result + integral(occ, un_occ) * aux_vector1(j_basis) * &
!!$            aux_vector3(k_basis) / (KS_eigenvalue(occ, 1, 1) - KS_eigenvalue(un_occ, 1, 1))
!!$        end do ! un_occ
!!$      end do ! occ
!!$      aux_matrix2(j_basis, k_basis) = temp_result
!!$    end do ! k_basis
!!$  end do ! j_basis

  do occ = 1, n_states_occupied, 1
     do un_occ = n_states_occupied + 1, n_states, 1
        integral(occ, un_occ) = integral(occ, un_occ) / (KS_eigenvalue(occ, 1, 1) - KS_eigenvalue(un_occ, 1, 1))
     end do
  end do

  call dgemm('N', 'T', n_states_occupied, n_basis, n_states - n_states_occupied, 1.0d0, &
       integral, n_states_occupied, KS_eigenvector(:,n_states_occupied+1:n_states,1,1), n_basis, &
       0.0d0, aux_matrix, n_basis)

  call dgemm('N', 'N', n_basis, n_basis, n_states_occupied, 1.0d0, &
       KS_eigenvector(:,1:n_states_occupied,1,1), n_basis, aux_matrix, n_basis, &
       0.0d0, aux_matrix2, n_basis)

  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to calculate matrix R",&
       rtime-sub_time, clock_rtime-clock_sub_time, OL_norm)

  call get_timestamps (sub_time, clock_sub_time)
  ! 3.2 For every integration point calculate
  ! chi(r) = phi(r)^T * R * phi(r)
  i_full_points = 0
  do i_my_batch = 1, n_my_batches, 1
    do i_index = 1, batches(i_my_batch)%size, 1
      i_full_points = i_full_points + 1
      aux_vector1 = all_waves(:, i_full_points)
      call dgemv('N', n_basis, n_basis, 1.0d0, aux_matrix2, n_basis, aux_vector1, 1, 0.d0, aux_vector2, 1)
      chi(i_full_points) = ddot(n_basis, aux_vector1, 1, aux_vector2, 1)
    end do ! i_index
  end do ! i_my_batch

  call get_timestamps (rtime, clock_rtime)
  if (out_times) call output_times(deffmt, "Time to calculate final chi",&
       rtime-sub_time, clock_rtime-clock_sub_time, OL_norm)

  
  deallocate(ylm_tab)
  deallocate(index_lm)
  deallocate(integral)

end subroutine chi_times_v_hartree
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/dielectric_func
!  NAME
!    dielectric_func
!  SYNOPSIS

subroutine dielectric_func(x, eps, partition_tab, hartree_partition_tab)

  !  PURPOSE
  !    Returns \int dr' \epsilon(r,r') * x(r')
  !    Where \epsilon is the dielectric function, and x is some function defined over
  !    the integration grid.
  !
  !    The dielectric function is defined as
  !
  !      \epsilon(r, r'; 0) = \delta(r - r') - \int dr'' \chi_0(r, r''; 0) V_hartree(r'', r')
  !  USES

  use dimensions
  implicit none

  ! ARGUMENTS
  real*8, intent(in) :: x(n_full_points)
  real*8, intent(in) :: partition_tab(n_full_points)
  real*8, intent(in) :: hartree_partition_tab(n_full_points)
  real*8, intent(out) :: eps(n_full_points)

  !  INPUTS
  !    o x -- input function, e.g. density
  !    o partition_tab -- values of the partition function for integration
  !    o hartree_partition_tab -- hartree partition tab for hartree solver
  !  OUTPUTS
  !    o eps -- the dielectric function
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Created in November 2017
  !  SOURCE

  ! Local variables
  real*8 :: V_hartree(n_full_points)

  ! Initialize variables
  eps = 0.d0
  V_hartree = 0.d0

  ! Get hartree potential from solver
  call hartree_solver(hartree_partition_tab, x, V_hartree)

  ! Get chi * V_hartree
  call chi_times_v_hartree(V_hartree, eps, partition_tab)

  eps = x - eps 

end subroutine dielectric_func
!******
