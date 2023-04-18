!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Provides a general purpose integration routine. Most integrands
!!  are stored in integrands.f90.
!!
!!  Usage:
!!    1) Call initialize_integration()
!!    2) Call integrate with the desired integration flags.
!!       See integrands.f90.
!!    3) Call cleanup_integration.
!!
!!  COMMENTS
!!
!!  Arrays whose leading dimension could change from batch to batch
!!  are sometimes declared 1-dimensional. Example: we allocate
!!  gradient_basis_wave with the size
!!  n_max_compute_ham*3*n_max_batch_size. The dimensions are then
!!  resolved in the called subroutines.
!!
!!  It is assumed that size(rho,1)==size(rho_gradient,1)==2
!!  always. See initialize_MR().
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module integration

  use, intrinsic :: iso_c_binding, only: c_bool

  use aims_memory_tracking,  only: aims_allocate
  use basis,                 only: basis_atom, basis_deriv_ordered, &
       & basis_kinetic_ordered, basis_wave_ordered
  use constants,             only: light_speed_sq
  use dimensions,            only: l_wave_max, n_basis, n_basis_fns, &
       & n_centers, n_centers_basis_I, n_centers_integrals, n_max_batch_size, &
       & n_max_compute_atoms, n_max_compute_fns_ham, n_max_compute_ham, &
       & n_my_batches, n_periodic, n_spin, use_gga
  use geometry,              only: coords, species
  use gpuMR
  use grids,                 only: batches, n_grid
  use integrands
  use load_balancing,        only: batch_perm, n_bp_integ, permute_point_array
  use mpi_tasks,             only: aims_stop
  use tools,                 only: batch_to_matrix_packed, &
       & block_cyclic_to_packed, if_present, matrix_packed_to_batch, &
       & packed_to_block_cyclic, safe_deallocate
  use pbc_lists,             only: Cbasis_to_center, centers_basis_integrals, &
       & inv_centers_basis_integrals, ptr_to_array, species_center
  use physics,               only: hartree_potential_std => hartree_potential, &
       & partition_tab, rho_std => rho, rho_gradient_std => rho_gradient
  use runtime_choices,       only: flag_rel, flag_xc, REL_atomic_zora, &
       & use_gpu, use_load_balancing, use_local_index, use_scalapack
  use scalapack_wrapper,     only: mxld, mxcol
  use species_data,          only: l_shell_max, species_z
  use synchronize_mpi_basic, only: sync_vector
  use timing,                only: start_timer, stop_timer
  use xc_f03_lib_m
  use xc_kernel_rho
  use xc_wrapper,            only: initialize_xc_wrapper, evaluate_xc_wrapper, &
       & cleanup_xc_wrapper

  ! TESTING

  use spline, only: val_spline
  use free_atoms, only: free_potential_spl

  ! END TESTING

  implicit none

  private
  public :: initialize_integration, integrate, cleanup_integration

  ! In most integration routines in aims, one can find conditionals
  ! such as "if (partition_tab(i_full_points) > 0d0)" inside a do loop
  ! over grid points. The array nonzero_indices avoids this by
  ! processing the partition function over grid points only once and
  ! then saving the indices corresponding to nonzero values of the
  ! global arrays. The integration routine can then be sped up by
  ! knowing in advance which array elements to consider. Example:
  ! nonzero_indices(2)%arr(3) corresponds to the 3rd nonzero element
  ! in batch 2.
  type(ptr_to_array), allocatable :: nonzero_indices(:)
  ! Most xc-kernel related work arrays need to be computed only
  ! once. After the first run, this variable is set to false.
  logical :: xc_kernel_rho_first_run

contains

  !!  FUNCTION
  !!
  !!  Allocates the main work variables that are required for every
  !!  call to 'integrate'.
  !!
  subroutine initialize_integration()
    character(*), parameter :: &
         & THIS_SUB = 'integration::init_integration::'
    integer, allocatable :: nonzero_indices_tmp(:)
    integer :: i_batch, i_full_points, i_point, i_index, i_tmp(1)
    call initialize_integrands()
    ! Using the partition function, determine the indices
    ! corresponding to nonzero values of any global array.
    allocate(nonzero_indices(batch_perm(n_bp_integ)%n_my_batches))
    call aims_allocate(nonzero_indices_tmp, n_max_batch_size, &
         & THIS_SUB//'nonzero_indices_tmp')
    i_full_points = 1
    do i_batch = 1, batch_perm(n_bp_integ)%n_my_batches
       nonzero_indices_tmp = 0
       i_point = 1
       do i_index = 1, batch_perm(n_bp_integ)%batches(i_batch)%size
          if (batch_perm(n_bp_integ)%partition_tab(i_full_points) > 0d0) then
             nonzero_indices_tmp(i_point) = i_full_points
             i_point = i_point + 1
          end if
          i_full_points = i_full_points + 1
       end do
       i_tmp = minloc(nonzero_indices_tmp)
       ! If the smallest element is not zero, the are no nonzero basis
       ! functions in the batch. In this case, take the whole batch.
       if (nonzero_indices_tmp(i_tmp(1)) /= 0) &
            & i_tmp = batch_perm(n_bp_integ)%batches(i_batch)%size + 1
       i_tmp(1) = i_tmp(1) - 1
       call aims_allocate(nonzero_indices(i_batch)%arr, i_tmp(1), &
            & THIS_SUB//'nonzero_indices(i_batch)%arr')
       nonzero_indices(i_batch)%arr(:i_tmp(1)) = nonzero_indices_tmp(:i_tmp(1))
    end do
    call safe_deallocate(nonzero_indices_tmp)
  end subroutine initialize_integration

  !!  FUNCTION
  !!
  !!  Integrates a function over all batches and all grid points, in
  !!  parallel. The result of integration is put into a 2D block
  !!  cyclic array matrix_BC_out. Exceptions include, for example,
  !!  first order density, in which case the result is put into rho1.
  !!
  subroutine integrate(integrand, matrix_BC_out, timer, rho1, &
       & rho1_gradient, hartree_pot1, density_matrix, points_for_distance)
    ! This determines the function to be integrated and contains
    ! information about additional quantities that might be required
    ! for integration (Nabla wave, zora factor, ...)
    type(integrand_t), intent(in) :: integrand
    ! If present, this is the block cyclic array that contains the
    ! integration result.
    real(dp), intent(out), optional :: matrix_BC_out(mxld,mxcol,*)
    ! If present, measure timing (unsynced)
    real(dp), intent(in out), optional :: timer(4)
    ! First-order electron density
    real(dp), intent(in out), optional :: rho1(:,:,:)
    ! First-order electron density gradient
    real(dp), intent(in out), optional :: rho1_gradient(:,:,:,:)
    ! First order Hartree potential
    real(dp), intent(in), optional :: hartree_pot1(:,:)
    ! An input density matrix, depending on what needs to be
    ! calculated. For instance, if the target quantity is the
    ! first-order density, this should be the first-order density
    ! matrix.
    real(dp), intent(in), optional :: density_matrix(mxld,mxcol,*)
    ! If present, distances between the given points and all
    ! integration points are calculated and passed to the
    ! integrand. For instance, some function might require distances
    ! from a given atom, in which case this argument could be passed
    ! as coords(:,i_atom). In case of two atoms, use
    ! points_for_distance=[coords(:,atom1),coords(:,atom2)].
    real(dp), intent(in), optional :: points_for_distance(:)

    character(*), parameter :: THIS_SUB = 'integration::integrate::'

    ! General purpose packed matrix for integrations.
    real(dp), allocatable :: matrix_packed(:,:,:)
    ! n1*Nabla|phi>, where n1 is the first-order density matrix
    real(dp), allocatable :: n1_phi_gradient(:,:,:)
    ! If points_for_distance is present, contains the required
    ! distances. (3,n_points,atoms of interest)
    real(dp), allocatable :: distances(:,:,:)
    ! Part of the kinetic operator on a basis function.
    real(dp), allocatable :: radial_kinetic_wave(:), kinetic_wave(:)
    ! xc-potential
    real(dp), allocatable :: vxc(:,:), vsigma(:,:)
    ! Position relative to basis function n.
    real(dp), allocatable :: r_n_basis(:)
    ! Difference of the position vectors of basis function m and n, R_m-R_n.
    real(dp), allocatable :: r_mn_basis(:)
    ! Auxiliary work arrays for GIAOs
    real(dp), allocatable :: matrix_batch_GIAO(:)
    ! Matrix times wave (e.g., H|phi>)
    real(dp), allocatable :: matrix_wave(:)
    ! Matrix elements corresponding to the current batch
    real(dp), allocatable :: matrix_batch(:)
    ! Local variable for determined how many spin channels are
    ! considered in packing.
    integer :: n_spins_for_packing
    ! Is use_load_balancing, the following arrays are permuted
    real(dp), allocatable :: rho(:,:), rho_gradient(:,:,:)
    real(dp), allocatable :: hartree_potential(:)
    ! The following arrays contains values at the batch level
    real(dp), allocatable :: local_rho(:,:), local_rho_gradient(:,:,:)
    real(dp), allocatable :: hartree_pot_loc(:)
    ! First order Hartree potential
    real(dp), allocatable :: hartree_pot1_loc(:)
    integer :: l_ylm_max
    integer, allocatable :: index_lm(:,:), i_basis(:), i_basis_fns(:), &
         & i_basis_fns_inv(:,:), i_atom_fns(:)
    real(dp), allocatable :: ylm_tab(:,:)
    real(dp), allocatable :: dylm_dtheta_tab(:,:), scaled_dylm_dphi_tab(:,:)
    real(dp), allocatable :: grid_points(:,:)
    real(dp), allocatable :: dist_tab(:,:), dist_tab_sq(:,:)
    real(dp), allocatable :: dir_tab(:,:,:)
    real(dp), allocatable :: radial_wave(:), radial_wave_deriv(:)
    real(dp), allocatable :: wave(:)
    real(dp), allocatable :: gradient_basis_wave(:)
    real(dp), allocatable :: zora_factor(:), zora_factor_tmp(:)
    real(dp), allocatable :: dist_tab_full(:)
    real(dp), allocatable :: dir_tab_full_norm(:,:)
    real(dp), allocatable :: i_r_full(:)
    ! Optimal accounting for matrix multiplications: only use points
    ! with nonzero components
    integer :: n_points
    ! For pruning of atoms, radial functions, and basis functions, to
    ! only the relevant ones ...
    integer :: n_compute, n_compute_fns
    integer :: n_compute_atoms
    integer :: atom_index(n_centers_integrals), atom_index_inv(n_centers)
    integer :: spline_array_start(n_centers_integrals)
    integer :: spline_array_end(n_centers_integrals)
    real(dp) :: i_r(n_max_compute_atoms)
    real(dp) :: coord_current(3)
    real(dp) :: trigonom_tab(4,n_max_compute_atoms)
    real(dp) :: one_over_dist_tab(n_max_compute_atoms)
    ! ZORA
    real(dp) :: partition(n_max_batch_size)
    ! Indices for basis functions that are nonzero at current point
    integer :: rad_index(n_max_compute_atoms)
    integer :: wave_index(n_max_compute_fns_ham)
    integer :: l_index(n_max_compute_fns_ham)
    integer :: l_count(n_max_compute_fns_ham)
    integer :: fn_atom(n_max_compute_fns_ham)
    ! Indices for known zero basis functions at current point
    integer :: n_zero_compute
    integer :: zero_index_point(n_max_compute_ham)
    ! Active atoms in current batch
    integer :: n_batch_centers
    integer :: batch_center(n_centers_integrals)
    ! Atomic ZORA
    integer :: i_center_2, i_center_L
    ! Counters
    integer :: i_atom, i_index, i_index2, i_l, i_m, i_point, i_spin, i_batch, &
         & i_dir, i_point_shift, i_shift
    ! For GPU
    integer :: matrix_batch_GIAO_size, r_mn_basis_size

    ! STEP 1 - Allocations and initialization
    if (integrand%spin_sensitive) then
       n_spins_for_packing = n_spin
    else
       n_spins_for_packing = 1
    end if
    ! (n_max_compute_ham,n_max_batch_size,X,n_spin), where X is equal
    ! to either n_dimensions or, in case of some GIAO integrals,
    ! GIAO_work_size.
    if (integrand%GIAO_work_size > 0) then
       i_index = n_max_compute_ham*n_max_batch_size*integrand%GIAO_work_size* &
            & n_spins_for_packing
    else
       i_index = n_max_compute_ham*n_max_batch_size*integrand%n_dimensions* &
            & n_spins_for_packing
    end if
    call aims_allocate(matrix_wave, i_index, THIS_SUB//'matrix_wave')
    if (.not. use_local_index) then
       call aims_allocate(matrix_packed, n_basis, n_basis, &
            & n_spins_for_packing*integrand%n_dimensions, &
            & THIS_SUB//'matrix_packed')
    else if (integrand%store_full_matrix) then
       call aims_allocate(matrix_packed, batch_perm(n_bp_integ)%  &
            & n_local_matrix_size, 2, n_spins_for_packing* &
            & integrand%n_dimensions, THIS_SUB//'matrix_packed')
    else
       call aims_allocate(matrix_packed, batch_perm(n_bp_integ)%  &
            & n_local_matrix_size, 1, n_spins_for_packing* &
            & integrand%n_dimensions, THIS_SUB//'matrix_packed')
    end if
    matrix_packed = 0d0
    ! (n_max_compute_ham,n_max_compute_ham)
    call aims_allocate(matrix_batch, n_max_compute_ham**2, &
         & THIS_SUB//'matrix_batch')
    call aims_allocate(dist_tab, n_centers_integrals, n_max_batch_size, &
         & THIS_SUB//'dist_tab')
    call aims_allocate(dist_tab_sq, n_centers_integrals, n_max_batch_size, &
         & THIS_SUB//'dist_tab_sq')
    call aims_allocate(dir_tab, 3, n_centers_integrals, n_max_batch_size, &
         & THIS_SUB//'dir_tab')
    call aims_allocate(i_basis_fns, n_basis_fns*n_centers_integrals, &
         & THIS_SUB//'i_basis_fns')
    call aims_allocate(i_basis_fns_inv, n_basis_fns, n_centers, &
         & THIS_SUB//'i_basis_fns_inv')
    call aims_allocate(i_atom_fns, n_basis_fns*n_centers_integrals, &
         & THIS_SUB//'i_atom_fns')
    call aims_allocate(local_rho, 2, n_max_batch_size, THIS_SUB//'local_rho')
    if (use_gga) then
       call aims_allocate(local_rho_gradient, 3, 2, n_max_batch_size, &
            & THIS_SUB//'local_rho_gradient')
    else
       allocate(local_rho_gradient(0,0,0))
    end if
    l_ylm_max = l_wave_max
    ! (n_max_compute_ham, 3, n_max_batch_size)
    call aims_allocate(gradient_basis_wave, 3*n_max_compute_ham* &
         & n_max_batch_size, THIS_SUB//'gradient_basis_wave')
    call aims_allocate(dylm_dtheta_tab, (l_ylm_max+1)**2, n_max_compute_atoms, &
         & THIS_SUB//'dylm_dtheta_tab')
    call aims_allocate(scaled_dylm_dphi_tab, (l_ylm_max+1)**2, &
         & n_max_compute_atoms, THIS_SUB//'scaled_dylm_dphi_tab')
    call aims_allocate(ylm_tab, (l_ylm_max+1)**2, n_max_compute_atoms, &
         & THIS_SUB//'ylm_tab')
    call aims_allocate(index_lm, -l_ylm_max, l_ylm_max, 0, l_ylm_max, &
         & THIS_SUB//'index_lm')
    call aims_allocate(radial_wave, n_max_compute_fns_ham, &
         & THIS_SUB//'radial_wave')
    call aims_allocate(radial_wave_deriv, n_max_compute_fns_ham, &
         & THIS_SUB//'radial_wave_deriv')
    ! (n_max_compute_ham, n_max_batch_size)
    call aims_allocate(wave, n_max_compute_ham*n_max_batch_size, &
         & THIS_SUB//'wave')
    call aims_allocate(grid_points, 3, n_max_batch_size, &
         & THIS_SUB//'grid_points')
    call aims_allocate(i_basis, n_centers_basis_I, THIS_SUB//'i_basis')
    if (flag_rel == REL_atomic_zora) then
       call aims_allocate(dist_tab_full, n_centers_integrals, &
            & THIS_SUB//'dist_tab_full')
       call aims_allocate(dir_tab_full_norm, 3, n_centers_integrals, &
            & THIS_SUB//'dir_tab_full_norm')
       call aims_allocate(i_r_full, n_centers_integrals, THIS_SUB//'i_r_full')
       call aims_allocate(zora_factor, n_max_compute_ham*n_max_batch_size, &
            & THIS_SUB//'zora_factor')
       call aims_allocate(zora_factor_tmp, n_centers_integrals, &
            & THIS_SUB//'zora_factor_tmp')
    end if
    if (present(points_for_distance)) &
         & call aims_allocate(distances, 3, n_max_batch_size, &
         & size(points_for_distance)/3, THIS_SUB//'distances')
    if (integrand%id == XC_KERNEL) &
         & xc_kernel_rho_first_run = is_first_xc_kernel_run()
    if (integrand%id == ELECTRON_DENSITY) then
       if (.not. present(rho1)) call aims_stop('rho1 must be present for &
            &''integrate''.', THIS_SUB)
       if (use_gga) call aims_allocate(n1_phi_gradient, 3, n_max_compute_ham, &
            & n_max_batch_size, THIS_SUB//'n1_phi_gradient')
    end if
    if (integrand%hamiltonian0 .or. integrand%rho) then
       call aims_allocate(rho, size(rho_std,1), &
            & batch_perm(n_bp_integ)%n_full_points, THIS_SUB//'rho')
       if (use_gga) &
            & call aims_allocate(rho_gradient, 3, size(rho_gradient_std,2), &
            & batch_perm(n_bp_integ)%n_full_points, THIS_SUB//'rho_gradient')
       if (use_load_balancing) then
          call permute_point_array(n_bp_integ, size(rho_std,1), rho_std, rho)
          if (use_gga) call permute_point_array(n_bp_integ, &
               & 3*size(rho_gradient_std,2), rho_gradient_std, rho_gradient)
       else
          rho = rho_std
          if (use_gga) rho_gradient = rho_gradient_std
       end if
    end if
    if (integrand%hamiltonian0) then
       if (flag_xc /= 0) call initialize_xc_wrapper()
       call aims_allocate(vxc, 2, n_max_batch_size, THIS_SUB//'vxc')
       call aims_allocate(radial_kinetic_wave, n_max_compute_fns_ham, &
            & THIS_SUB//'radial_kinetic_wave')
       ! (n_max_compute_ham, n_max_batch_size)
       call aims_allocate(kinetic_wave, n_max_compute_ham*n_max_batch_size, &
            & THIS_SUB//'kinetic_wave')
       if (use_gga) then
          call aims_allocate(vsigma, 3, n_max_batch_size, THIS_SUB//'vsigma')
       else
          allocate(vsigma(0,0))
       end if
       call aims_allocate(hartree_pot_loc, n_max_batch_size, &
            & THIS_SUB//'hartree_pot_loc')
       call aims_allocate(hartree_potential, batch_perm(n_bp_integ)% &
            & n_full_points, THIS_SUB//'hartree_potential')
       if (use_load_balancing) then
          call permute_point_array(n_bp_integ, 1, hartree_potential_std, &
               & hartree_potential)
       else
          hartree_potential = hartree_potential_std
       end if
    end if
    if (integrand%r_n_basis) & ! (n_max_compute_ham, 3, n_max_batch_size)
         & call aims_allocate(r_n_basis, 3*n_max_compute_ham*n_max_batch_size, &
         & THIS_SUB//'r_n_basis')
    if (integrand%r_mn_basis) then
       ! (n_max_compute_ham, n_max_compute_ham, 3)
       call aims_allocate(r_mn_basis, 3*n_max_compute_ham**2, &
            & THIS_SUB//'r_mn_basis')
       if (integrand%id == DIAMAGNETIC_MAGNETIZABILITY .and. use_gga) then
          ! With magnetizability we have integrals of type R^T x H x R,
          ! where the cross products are taken between a tensor, H,
          ! and two vectors, R. These require a lot more
          ! workspace. Particularly problematic is
          ! f_diamagnetic_magnetizability_gga (for
          ! use_gga=.true. only), which doesn't allow any sort of
          ! packing. Thus all tensor components must be kept in
          ! memory.
          if (integrand%calc_full_tensor) then
             call aims_allocate(matrix_batch_GIAO, 12*n_max_compute_ham**2, &
                  & THIS_SUB//'matrix_batch_GIAO')
          else
             call aims_allocate(matrix_batch_GIAO, 11*n_max_compute_ham**2, &
                  & THIS_SUB//'matrix_batch_GIAO')
          end if
       else if (any(integrand%id == &
            & [DIAMAGNETIC_MAGNETIZABILITY, GIAO_OVERLAP2])) then
          if (integrand%calc_full_tensor) then
             call aims_allocate(matrix_batch_GIAO, 9*n_max_compute_ham**2, &
                  & THIS_SUB//'matrix_batch_GIAO')
          else
             call aims_allocate(matrix_batch_GIAO, 8*n_max_compute_ham**2, &
                  & THIS_SUB//'matrix_batch_GIAO')
          end if
       else
          call aims_allocate(matrix_batch_GIAO, 3*n_max_compute_ham**2, &
               & THIS_SUB//'matrix_batch_GIAO')
       end if
    end if
    ! (n_max_batch_size, integrand%n_dimensions)
    if (present(hartree_pot1)) call aims_allocate(hartree_pot1_loc, &
         & integrand%n_dimensions*n_max_batch_size, &
         & THIS_SUB//'hartree_pot1_loc')
    if (present(density_matrix)) then
       do i_dir = 1, n_spins_for_packing*integrand%n_dimensions
          call block_cyclic_to_packed(density_matrix(:,:,i_dir), &
               & matrix_packed(:,:,i_dir))
       end do
    end if
    i_basis_fns_inv = 0
    i_index = 0
    do i_l = 0, l_wave_max
       do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
       end do
    end do
    i_point_shift = 0

    if (use_gpu) then
       if (allocated(matrix_batch_GIAO)) then
          matrix_batch_GIAO_size = size(matrix_batch_GIAO)
       else
          matrix_batch_GIAO_size = 0
       end if
       if (allocated(r_mn_basis)) then
          r_mn_basis_size = size(r_mn_basis)
       else
          r_mn_basis_size = 0
       end if
       associate(A => batch_perm(n_bp_integ)%i_basis_glb_to_loc, &
            & matrix_wave_size_f => &
            & size(matrix_wave)/(n_max_compute_ham*n_max_batch_size))
         call mr_initialize_gpu(size(wave), size(matrix_wave), &
              & size(matrix_batch), matrix_batch_GIAO_size, &
              & r_mn_basis_size, size(matrix_packed), size(i_basis), size(A), &
              & A, matrix_wave_size_f)
       end associate
    end if

    ! STEP 2 - Perform partitioned integration, batch by batch. It is
    !          assumed that the developer is more or less familiar
    !          with the subroutines in the following loop. For more
    !          details see the original subroutine
    !          integrate_real_hamiltonian_matrix_p2.
    if (present(timer)) call start_timer(timer)
    over_batches: do i_batch = 1, batch_perm(n_bp_integ)%n_my_batches
       ! STEP 2a - Identify the non-zero basis functions in the given
       !           batch, if there are any.
       if (i_batch > 1) i_point_shift = i_point_shift + &
            & batch_perm(n_bp_integ)%batches(i_batch-1)%size
       do i_point = 1, size(nonzero_indices(i_batch)%arr)
          coord_current = batch_perm(n_bp_integ)%batches(i_batch)% &
               & points(nonzero_indices(i_batch)%arr(i_point)-i_point_shift)% &
               & coords
          if (n_periodic > 0) call map_to_center_cell(coord_current)
          ! Compute atom-centered coordinates of current integration
          ! point, as viewed from all atoms.
          call tab_atom_centered_coords_p0(coord_current, &
               & dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
               & n_centers_integrals, centers_basis_integrals)
       end do
       n_compute = batch_perm(n_bp_integ)%batches(i_batch)%batch_n_compute
       i_basis(:n_compute) = batch_perm(n_bp_integ)%batches(i_batch)% &
            & batch_i_basis
       ! From list of n_compute active basis functions in batch,
       ! collect all atoms that are ever needed in batch.
       call collect_batch_centers_p2(n_compute, i_basis, n_centers_basis_I, &
            & n_centers_integrals, inv_centers_basis_integrals, &
            & n_batch_centers, batch_center)
       if (n_compute == 0) cycle
       n_points = size(nonzero_indices(i_batch)%arr)
       ! STEP 2b - With at least one non-zero basis function in the
       !           given batch, proceed to find the basic quantities
       !           required for integration.
       over_batch: do i_point = 1, n_points
          coord_current = batch_perm(n_bp_integ)%batches(i_batch)% &
               & points(nonzero_indices(i_batch)%arr(i_point)-i_point_shift)% &
               & coords
          if (n_periodic > 0) call map_to_center_cell(coord_current)
          grid_points(:,i_point) = coord_current
          n_compute_atoms = 0
          n_compute_fns = 0
          call prune_radial_basis_p2 (n_max_compute_atoms, &
               & n_max_compute_fns_ham, dist_tab_sq(1,i_point), &
               & dist_tab(1,i_point), dir_tab(1,1,i_point), n_compute_atoms, &
               & atom_index, atom_index_inv, n_compute_fns, i_basis_fns, &
               & i_basis_fns_inv, i_atom_fns, spline_array_start, &
               & spline_array_end, n_centers_integrals, &
               & centers_basis_integrals, n_compute, i_basis, n_batch_centers, &
               & batch_center, one_over_dist_tab, rad_index, wave_index, &
               & l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
          ! Tabulate distances, unit vectors, and inverse logarithmic
          ! grid units for all atoms currently relevant.
          call tab_local_geometry_p2(n_compute_atoms, atom_index, &
               & dist_tab(1,i_point), i_r)
          ! Compute trigonometric functions of spherical coordinate
          ! angles of current integration point, viewed from all
          ! atoms.
          call tab_trigonom_p0(n_compute_atoms, dir_tab(1,1,i_point), &
               & trigonom_tab)
          ! Tabulate those ylms needed for gradients, i.e. ylm's for
          ! l_max+1.
          call tab_gradient_ylm_p0(trigonom_tab(1,1), l_shell_max, l_ylm_max, &
               & n_compute_atoms, atom_index, ylm_tab(1,1), &
               & dylm_dtheta_tab(1,1), scaled_dylm_dphi_tab(1,1))
          ! Construct 'wave' from the compressed spline arrays.
          call evaluate_radial_functions_p0(spline_array_start, &
               & spline_array_end, n_compute_atoms, n_compute_fns, &
               & dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, &
               & basis_wave_ordered, radial_wave, .false. , n_compute, &
               & n_max_compute_fns_ham)
          call evaluate_waves_p2(n_compute, n_compute_atoms, n_compute_fns, &
               & l_ylm_max, ylm_tab, one_over_dist_tab, radial_wave, &
               & wave(n_compute*(i_point-1)+1), rad_index, wave_index, &
               & l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
          ! Construct 'gradient_basis_wave' from the compressed spline
          ! arrays.
          if (integrand%basis_gradient) then
             call evaluate_radial_functions_p0(spline_array_start, &
                  & spline_array_end, n_compute_atoms, n_compute_fns, &
                  & dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, &
                  & basis_deriv_ordered, radial_wave_deriv(1), .true., &
                  & n_compute, n_max_compute_fns_ham)
             call evaluate_wave_gradient_p2(n_compute, n_compute_atoms, &
                  & n_compute_fns, one_over_dist_tab, dir_tab(1,1,i_point), &
                  & trigonom_tab(1,1), l_ylm_max, ylm_tab, dylm_dtheta_tab, &
                  & scaled_dylm_dphi_tab, radial_wave, radial_wave_deriv, &
                  & gradient_basis_wave(3*n_compute*(i_point-1)+1), &
                  & rad_index, wave_index, l_index, l_count, fn_atom, &
                  & n_zero_compute, zero_index_point)
          end if
          ! Construct 'kinetic_wave' from the compressed spline
          ! arrays.
          if (integrand%hamiltonian0) then
             call evaluate_radial_functions_p0(spline_array_start, &
                  & spline_array_end, n_compute_atoms, n_compute_fns, &
                  & dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, &
                  & basis_kinetic_ordered, radial_kinetic_wave, .false., &
                  & n_compute, n_max_compute_fns_ham)
             call evaluate_waves_p2(n_compute, n_compute_atoms, n_compute_fns, &
                  & l_ylm_max, ylm_tab, one_over_dist_tab, &
                  & radial_kinetic_wave, &
                  & kinetic_wave(n_compute*(i_point-1)+1), rad_index, &
                  & wave_index, l_index, l_count, fn_atom, n_zero_compute, &
                  & zero_index_point)
          end if
          if (flag_rel == REL_atomic_zora) then
             call tab_global_geometry_p0(dist_tab_sq(1,i_point), &
                  & dir_tab(1,1,i_point), dist_tab_full, i_r_full, &
                  & dir_tab_full_norm, n_centers_integrals, &
                  & centers_basis_integrals)
             do i_center_L = 1, n_centers_integrals
                i_center_2 = centers_basis_integrals(i_center_L)
                zora_factor_tmp(i_center_L) = &
                     & val_spline(i_r_full(i_center_L), &
                     & free_potential_spl(1,1,species_center(i_center_2)), &
                     & n_grid(species_center(i_center_2)))
             end do
             zora_factor_tmp = &
                  & 2d0*light_speed_sq/(2d0*light_speed_sq - zora_factor_tmp)
             do i_m = 1, n_compute
                zora_factor(n_compute*(i_point-1)+i_m) = zora_factor_tmp( &
                     & inv_centers_basis_integrals(Cbasis_to_center(i_m)))
             end do
          end if
          i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
       end do over_batch
       ! STEP 2c - We know the grid points and other the basic
       !           quantities in the given batch. Next calculate
       !           additional quantities that are required for the
       !           given function before performing the actual
       !           integration.
       partition(:n_points) = batch_perm(n_bp_integ)% &
            & partition_tab(nonzero_indices(i_batch)%arr)
       if (present(hartree_pot1)) then
          do i_dir = 0, integrand%n_dimensions-1
             hartree_pot1_loc(i_dir*n_points+1:i_dir*n_points+n_points) = &
                  & hartree_pot1(nonzero_indices(i_batch)%arr,i_dir+1)
          end do
       end if
       if (integrand%hamiltonian0) then
          local_rho(:,:n_points) = rho(:,nonzero_indices(i_batch)%arr)
          if (use_gga) &
               & local_rho_gradient(:,:,:n_points) = &
               & rho_gradient(:,:,nonzero_indices(i_batch)%arr)
          if (flag_xc /= 0) then
             call evaluate_xc_wrapper(n_points, local_rho, local_rho_gradient, &
                  & vxc, vsigma)
          else
             vxc(:,:n_points) = 0d0 ! Hartree-Fock
          end if
          hartree_pot_loc(:n_points) = &
               & hartree_potential(nonzero_indices(i_batch)%arr)
       end if
       if (integrand%id == XC_KERNEL) then
          if (.not. use_gga) then
             call xc_kernel_rho_LDA(xc_kernel_rho_first_run, &
                  & integrand%n_dimensions, n_points, &
                  & nonzero_indices(i_batch)%arr, rho1, rho, local_rho)
          else
             call xc_kernel_rho_GGA(xc_kernel_rho_first_run, &
                  & integrand%n_dimensions, n_points, &
                  & nonzero_indices(i_batch)%arr, rho1, rho1_gradient, rho, &
                  & rho_gradient, local_rho, local_rho_gradient)
          end if
       end if
       ! If some points were passed as an argument, find the
       ! distances to those points.
       if (present(points_for_distance)) then
          do i_atom = 1, size(points_for_distance)/3
             do i_point = 1, n_points
                distances(:,i_point,i_atom) = grid_points(:,i_point) - &
                     & points_for_distance(3*(i_atom-1)+1:3*(i_atom-1)+3)
             end do
          end do
       end if
       if (integrand%r_n_basis) then
          do i_point = 1, n_points
             do i_dir = 1, 3
                do i_index = 1, n_compute
                   r_n_basis(i_index + (i_dir-1)*n_compute + &
                        & (i_point-1)*3*n_compute) = &
                        & grid_points(i_dir,i_point) - &
                        & coords(i_dir,basis_atom(i_basis(i_index)))
                end do
             end do
          end do
       end if
       if (integrand%r_mn_basis) then
          do i_index2 = 1, n_compute
             do i_index = 1, n_compute
                do i_dir = 1, 3
                   r_mn_basis(i_index + (i_index2-1)*n_compute + &
                        & (i_dir-1)*n_compute**2) = &
                        & coords(i_dir,basis_atom(i_basis(i_index))) - &
                        & coords(i_dir,basis_atom(i_basis(i_index2)))
                end do
             end do
          end do
       end if
       ! STEP 2d - We now have all the required quantities. Perform
       !           integration over the given batch of points and
       !           nonzero basis functions.
       select case(integrand%id)
       case(FERMI_CONTACT_SR)
          call f_fermi_contact_sr(n_compute, n_points, partition, &
               & gradient_basis_wave, distances, zora_factor, matrix_wave)
       case(PARAMAGNETIC_SO)
          call f_paramagnetic_SO(n_compute, n_points, partition, &
               & gradient_basis_wave, distances, matrix_wave)
       case(PARAMAGNETIC_SO_SR)
          call f_paramagnetic_SO_sr(n_compute, n_points, partition, &
               & gradient_basis_wave, distances, zora_factor, matrix_wave)
       case(SPIN_DIPOLE_XX_YY_ZZ)
          call f_spin_dipole_xx_yy_zz(n_compute, n_points, partition, wave, &
               & distances, matrix_wave)
       case(SPIN_DIPOLE_XY_XZ_YZ)
          call f_spin_dipole_xy_xz_yz(n_compute, n_points, partition, wave, &
               & distances, matrix_wave)
       case(DIAMAGNETIC_SO)
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_SO_full(n_compute, n_points, partition, wave, &
                  & distances, distances(1,1,2), matrix_wave)
          else
             call f_diamagnetic_SO_diag(n_compute, n_points, partition, wave, &
                  & distances, distances(1,1,2), matrix_wave)
          end if
       case(DIAMAGNETIC_SO_SR)
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_SO_sr_full(n_compute, n_points, partition, &
                  & wave, distances, distances(1,1,2), zora_factor, matrix_wave)
          else
             call f_diamagnetic_SO_sr_diag(n_compute, n_points, partition, &
                  & wave, distances, distances(1,1,2), zora_factor, matrix_wave)
          end if
       case(ELECTRIC_FIELD_GRADIENT)
          if (integrand%calc_full_tensor) then
             call f_electric_field_gradient_full(n_compute, n_points, &
                  & partition, wave, distances, matrix_wave)
          else
             call f_electric_field_gradient_diag(n_compute, n_points, &
                  & partition, wave, distances, matrix_wave)
          end if
       case(ORBITAL_ZEEMAN)
          call f_orbital_zeeman(n_compute, n_points, partition, &
               & gradient_basis_wave, distances, matrix_wave)
       case(GIAO_OVERLAP1)
          call f_giao_overlap1(n_compute, n_points, partition, wave, &
               & grid_points, matrix_wave)
          call GIAO_mult_psi_and_pack(.true., .false., .false., .false.)
          ! We already called GIAO_mult_psi_and_pack, don't do STEP 3.
          cycle
       case(GIAO_OVERLAP2)
          call f_giao_overlap2(n_compute, n_points, partition, wave, &
               & grid_points, matrix_wave)
          call GIAO_mult_psi_and_pack(.true., .true., .false., .true.)
          cycle
       case(GIAO_PARAMAGNETIC)
          if (.not. use_gga) then
             call f_giao_paramagnetic_lda(n_compute, n_points, partition, &
                  & wave, kinetic_wave, hartree_pot_loc, vxc, grid_points, &
                  & matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .false., .false., .false., 3)
          else
             call f_giao_paramagnetic_gga(n_compute, n_points, partition, &
                  & wave, kinetic_wave, hartree_pot_loc, vxc, vsigma, &
                  & grid_points, local_rho_gradient, matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .false., .false., .false., 3)
             call f_giao_paramagnetic_gga_grad(n_compute, n_points, &
                  & partition, gradient_basis_wave, vsigma, matrix_wave, &
                  & grid_points, local_rho_gradient)
             call GIAO_mult_psi_and_pack(.true., .false., .true., .false., 3)
          end if
          call f_giao_paramagnetic_0(n_compute, n_points, partition, &
               & gradient_basis_wave, r_n_basis, matrix_wave)
          call GIAO_mult_psi_and_pack(.false., .false., .false., .false.)
          cycle
       case(DIAMAGNETIC_SHIELDING_STD)
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_shielding_std_full(n_compute, n_points, &
                  & partition, wave, distances, distances(1,1,2), matrix_wave)
          else
             call f_diamagnetic_shielding_std_diag(n_compute, n_points, &
                  & partition, wave, distances, distances(1,1,2), matrix_wave)
          end if
       case(DIAMAGNETIC_SHIELDING)
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_shielding_1_full(n_compute, n_points, &
                  & partition, gradient_basis_wave, grid_points, distances, &
                  & matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .false., .false., .false.)
             call f_diamagnetic_shielding_2_full(n_compute, n_points, &
                  & partition, wave, distances, r_n_basis, matrix_wave)
          else
             call f_diamagnetic_shielding_1_diag(n_compute, n_points, &
                  & partition, gradient_basis_wave, grid_points, distances, &
                  & matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .true., .false., .false.)
             call f_diamagnetic_shielding_2_diag(n_compute, n_points, &
                  & partition, wave, distances, r_n_basis, matrix_wave)
          end if
       case(DIAMAGNETIC_MAGNETIZABILITY_STD)
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_magnetizability_std_full(n_compute, n_points, &
                  & partition, wave, distances, matrix_wave)
          else
             call f_diamagnetic_magnetizability_std_diag(n_compute, n_points, &
                  & partition, wave, distances, matrix_wave)
          end if
       case(DIAMAGNETIC_MAGNETIZABILITY)
          ! This subroutine is always called, even with GGA. The GGA
          ! terms are added later.
          call f_diamagnetic_magnetizability_lda(n_compute, n_points, &
               & partition, wave, kinetic_wave, hartree_pot_loc, vxc, &
               & grid_points, matrix_wave)
          call GIAO_mult_psi_and_pack(.true., .true., .false., .true., 6)
          if (use_gga) then
             call f_diamagnetic_magnetizability_gga(n_compute, &
                  & n_points, partition, wave, vsigma, grid_points, &
                  & local_rho_gradient, matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .false., .false., .true., 9)
             call f_diamagnetic_magnetizability_gga_grad(n_compute, &
                  & n_points, partition, gradient_basis_wave, vsigma, &
                  & grid_points, local_rho_gradient, matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .true., .true., .true., 6)
          end if
          if (integrand%calc_full_tensor) then
             call f_diamagnetic_magnetizability_R_mn_full(n_compute, n_points, &
                  & partition, gradient_basis_wave, grid_points, r_n_basis, &
                  & matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .false., .false., .false.)
             call f_diamagnetic_magnetizability_0_full(n_compute, n_points, &
                  & partition, wave, r_n_basis, matrix_wave)
          else
             call f_diamagnetic_magnetizability_R_mn_diag(n_compute, n_points, &
                  & partition, gradient_basis_wave, grid_points, r_n_basis, &
                  & matrix_wave)
             call GIAO_mult_psi_and_pack(.true., .true., .false., .false.)
             call f_diamagnetic_magnetizability_0_diag(n_compute, n_points, &
                  & partition, wave, r_n_basis, matrix_wave)
          end if
          if (n_spin == 2) then
             ! Even though f_diamagnetic_magnetizability_0 is not spin
             ! sensitive, other parts of the diamagnetic
             ! magnetizability integral are. For consistency, we force
             ! this term also to be spin sensitive by copying the
             ! values from spin 1 to spin 2.
             i_m = integrand%n_dimensions*n_compute*n_points
             matrix_wave(i_m+1:2*i_m) = matrix_wave(:i_m)
          end if
       case(XC_KERNEL)
          if (.not. use_gga) then
             do i_dir = 1, integrand%n_dimensions
                call f_xc_kernel_LDA(n_compute, n_points, partition, wave, &
                     & xc_kernel_s(:,:,i_dir), &
                     & matrix_wave((i_dir-1)*n_spin*n_compute*n_points+1:))
             end do
          else
             do i_dir = 1, integrand%n_dimensions
                call f_xc_kernel_GGA(n_compute, n_points, partition, wave, &
                     & gradient_basis_wave, xc_kernel_s(:,:,i_dir), &
                     & xc_kernel_v(:,:,:,i_dir), &
                     & matrix_wave((i_dir-1)*n_spin*n_compute*n_points+1:))
             end do
          end if
       case(FIRST_ORDER_HARTREE)
          call f_first_order_hartree(integrand%n_dimensions, n_compute, &
               & n_points, partition, wave, hartree_pot1_loc, matrix_wave)
       case(ELECTRON_DENSITY)
          do i_dir = 1, integrand%n_dimensions
             do i_spin = 1, n_spin
                ! matrix_packed contains the density matrix in the
                ! packed format.
                call matrix_packed_to_batch(n_compute, i_basis, &
                     & matrix_packed(:,:,i_spin+n_spin*(i_dir-1)), matrix_batch)
                if (.not. use_gga) then
                   call f_electron_density(n_compute, n_points, i_spin, wave, &
                        & matrix_batch, matrix_wave, local_rho)
                else
                   call f_electron_density(n_compute, n_points, i_spin, wave, &
                        & matrix_batch, matrix_wave, local_rho, &
                        & gradient_basis_wave, local_rho_gradient)
                end if
             end do
             rho1(:,nonzero_indices(i_batch)%arr,i_dir) = local_rho(:,:n_points)
             if (use_gga) &
                  & rho1_gradient(:,:,nonzero_indices(i_batch)%arr,i_dir) = &
                  & local_rho_gradient(:,:,:n_points)
          end do
          cycle ! Don't do STEP 3 in this case
       case(POSITION)
          call f_position(n_compute, n_points, partition, wave, grid_points, &
               & matrix_wave)
       end select
       ! STEP 3 - Elements of H|phi> have been locally computed. Next
       !          multiply from the left by <phi| and transfer results
       !          to a global packed array.
       if (use_gpu) then
          do i_spin = 1, n_spins_for_packing
             do i_dir = 1, integrand%n_dimensions
                ! Multiply H|phi> from the left by <phi|
                i_shift = ((i_dir-1) + (i_spin-1)*integrand%n_dimensions)* &
                     & n_compute*n_points
                call evaluate_mr_batch(n_points, n_compute, wave, &
                     & matrix_wave(1+i_shift), integrand%symmetrization, &
                     & i_dir, integrand%n_dimensions)
                i_shift = size(matrix_packed(:,:,1))* &
                     & (i_dir+(i_spin-1)*integrand%n_dimensions-1)
                call update_mr_batch_gpu(n_compute, i_shift, i_basis)
             end do
          end do
       else
          do i_spin = 1, n_spins_for_packing
             do i_dir = 1, integrand%n_dimensions
                ! Multiply H|phi> from the left by <phi|
                i_shift = ((i_dir-1) + (i_spin-1)*integrand%n_dimensions)* &
                     & n_compute*n_points
                call dgemm('n', 't', n_compute, n_compute, n_points, 1d0, &
                     & wave, n_compute, matrix_wave(1+i_shift), n_compute, &
                     & 0d0, matrix_batch, n_compute)
                ! Symmetrize if necessary
                if (integrand%symmetrization == 1) then
                   call symmetrize(n_compute, matrix_batch)
                else if (integrand%symmetrization == 2) then
                   call antisymmetrize(n_compute, matrix_batch)
                end if
                ! Do the packing
                i_shift = i_dir + (i_spin-1)*integrand%n_dimensions
                call batch_to_matrix_packed(n_compute, i_basis, matrix_batch, &
                     & matrix_packed(:,:,i_shift))
             end do
          end do
       end if
    end do over_batches
    if (present(timer)) call stop_timer(timer, .true.)

    if (use_gpu) then
       call get_matrix_packed_gpu(matrix_packed, size(matrix_packed))
       call mr_destroy_gpu()
    end if

    ! STEP 4 - Convert the result from packed to block cyclic
    if (present(matrix_BC_out)) then
       ! If packing is not employed, additional operations must be
       ! performed on matrix_packed.
       if (.not. use_local_index) call process_if_no_packing(matrix_packed)
       do i_dir = 1, integrand%n_dimensions*n_spins_for_packing
          call packed_to_block_cyclic(matrix_packed(:,:,i_dir), &
               & matrix_BC_out(:,:,i_dir))
       end do
    end if

    ! STEP 5 - Deallocations
    call safe_deallocate(matrix_wave)
    call safe_deallocate(matrix_batch)
    call safe_deallocate(dist_tab)
    call safe_deallocate(dist_tab_sq)
    call safe_deallocate(dir_tab)
    call safe_deallocate(i_basis_fns)
    call safe_deallocate(i_basis_fns_inv)
    call safe_deallocate(i_atom_fns)
    call safe_deallocate(local_rho)
    call safe_deallocate(local_rho_gradient)
    call safe_deallocate(gradient_basis_wave)
    call safe_deallocate(dylm_dtheta_tab)
    call safe_deallocate(scaled_dylm_dphi_tab)
    call safe_deallocate(ylm_tab)
    call safe_deallocate(index_lm)
    call safe_deallocate(radial_wave)
    call safe_deallocate(radial_wave_deriv)
    call safe_deallocate(wave)
    call safe_deallocate(grid_points)
    call safe_deallocate(i_basis)
    call safe_deallocate(dist_tab_full)
    call safe_deallocate(dir_tab_full_norm)
    call safe_deallocate(i_r_full)
    call safe_deallocate(zora_factor)
    call safe_deallocate(zora_factor_tmp)
    call safe_deallocate(distances)
    call safe_deallocate(n1_phi_gradient)
    call safe_deallocate(vxc)
    call safe_deallocate(radial_kinetic_wave)
    call safe_deallocate(kinetic_wave)
    call safe_deallocate(vsigma)
    call safe_deallocate(r_n_basis)
    call safe_deallocate(r_mn_basis)
    call safe_deallocate(matrix_batch_GIAO)
    call safe_deallocate(hartree_pot_loc)
    call safe_deallocate(hartree_pot1_loc)
    call safe_deallocate(hartree_potential)
    call safe_deallocate(rho)
    call safe_deallocate(rho_gradient)
    call safe_deallocate(matrix_packed)
  contains
    !!
    !!  FUNCTION
    !!
    !!  With GIAOs, some parts of STEP 3 above are complicated by the
    !!  fact that the Hamiltonian depends on both the left (bra) and
    !!  right (ket) sides. In particular, the GIAO matrix elements
    !!  often require cross products in the form R_mn x <m|H|n>. Thus,
    !!  separate subroutines are called for those cases.
    !!
    subroutine GIAO_mult_psi_and_pack(mult_R_mn_left, compact_packing, &
         & do_transpose, mult_R_mn_both, n_spin_shift)
      ! Whether to take the cross-products from the left: R_mn x <m|H|n>
      logical, intent(in) :: mult_R_mn_left
      ! In case of a diamagnetic property, when only diagonal tensor
      ! elements are requested, a compact indexing scheme is used.
      logical, intent(in) :: compact_packing
      ! Wether to symmetrize the Hamiltonian before the cross-product
      ! is taken: R_mn x <m|H+H^T|n>. This is relevant for a GGA term
      ! of a paramagnetic property.
      logical, intent(in) :: do_transpose
      ! Whether to take the cross-products from both sides:
      ! R_mn x <m|H|n> x R_mn
      logical, intent(in) :: mult_R_mn_both
      ! If the matrix elements are spin-dependent, process the two
      ! spin channels of matrix_wave separately. Spin two starts at
      ! n_spin_shift*n_compute*n_points+1 in matrix_wave. If spin
      ! up/down matrix elements are the same, but both channels still
      ! need to be explicitly recorded in the packed array, this
      ! argument may be omitted or set to 0.
      integer, intent(in), optional :: n_spin_shift
      integer :: i_spin, i_dim
      if (use_gpu) then
         call GIAO_mult_psi_and_pack_gpu(mult_R_mn_left, compact_packing, &
              & do_transpose, mult_R_mn_both, n_spin_shift)
         return
      end if
      do i_spin = 1, n_spins_for_packing
         do i_dim = 1, integrand%n_dimensions
            if (mult_R_mn_left .or. mult_R_mn_both) then
               i_shift = &
                    & (i_spin-1)*if_present(n_spin_shift,0)*n_compute*n_points+1
               call GIAO_mult_psi(compact_packing, do_transpose, &
                    & mult_R_mn_both, i_dim, integrand%n_dimensions, &
                    & n_compute, n_points, wave, matrix_wave(i_shift), &
                    & r_mn_basis, matrix_batch_GIAO, matrix_batch)
            else
               call dgemm('n', 't', n_compute, n_compute, n_points, 1d0, wave, &
                    & n_compute, matrix_wave((i_dim-1)*n_compute*n_points+1), &
                    & n_compute, 0d0, matrix_batch, n_compute)
            end if
            if (integrand%symmetrization == 1) then
               call symmetrize(n_compute, matrix_batch)
            else if (integrand%symmetrization == 2) then
               call antisymmetrize(n_compute, matrix_batch)
            end if
            i_shift = i_dim+(i_spin-1)*integrand%n_dimensions
            call batch_to_matrix_packed(n_compute, i_basis, matrix_batch, &
                 & matrix_packed(:,:,i_shift))
         end do
      end do
    end subroutine GIAO_mult_psi_and_pack

    subroutine GIAO_mult_psi_and_pack_gpu(mult_R_mn_left, compact_packing, &
         & do_transpose, mult_R_mn_both, n_spin_shift)
      logical, intent(in) :: mult_R_mn_left
      logical, intent(in) :: compact_packing
      logical, intent(in) :: do_transpose
      logical, intent(in) :: mult_R_mn_both
      integer, intent(in), optional :: n_spin_shift
      integer :: i_spin, i_dim
      logical(c_bool) :: c_compact_indexing, c_do_transpose, c_mult_R_mn_both
      c_compact_indexing = compact_packing
      c_do_transpose = do_transpose
      c_mult_R_mn_both = mult_R_mn_both
      do i_spin = 1, n_spins_for_packing
         do i_dim = 1, integrand%n_dimensions
            if (mult_R_mn_left .or. mult_R_mn_both) then
               i_shift = &
                    & (i_spin-1)*if_present(n_spin_shift,0)*n_compute*n_points+1
               call giao_mult_psi_gpu(c_compact_indexing, c_do_transpose, &
                    & c_mult_R_mn_both, i_dim, integrand%n_dimensions, &
                    & n_compute, n_points, wave, matrix_wave(i_shift), &
                    & r_mn_basis, matrix_batch_GIAO, matrix_batch, &
                    & n_spins_for_packing)
            else
               call evaluate_mr_batch_no_symm(n_points, n_compute, wave, &
                    & matrix_wave, integrand%symmetrization, i_dim, &
                    & integrand%n_dimensions)
            end if
            call symm_antisymm_gpu(integrand%symmetrization, n_compute)
            i_shift = size(matrix_packed(:,:,1))* &
                 & (i_dim+(i_spin-1)*integrand%n_dimensions-1)
            if (size(matrix_packed,2) == 2) then
               call update_mr_batch_gpu_full(n_compute, i_shift, i_basis, &
                    & size(matrix_packed,1))
            else
               call update_mr_batch_gpu(n_compute, i_shift, i_basis)
            end if
         end do
      end do
    end subroutine GIAO_mult_psi_and_pack_gpu

    subroutine process_if_no_packing(matrix)
      real(dp), intent(in out) :: matrix(:,:,:)
      integer :: i_index_2
      ! Unless the full matrix is requested to be stored explicitly,
      ! mirror the lower triangle onto the upper triangle.
      if (.not. integrand%store_full_matrix) then
         do i_index = 1, n_basis
            do i_index_2 = i_index+1, n_basis
               matrix(i_index,i_index_2,:) = matrix(i_index_2,i_index,:)
            end do
         end do
      end if
      ! If no packing, always synchronize
      call sync_vector(matrix, size(matrix))
    end subroutine process_if_no_packing
  end subroutine integrate

  pure subroutine symmetrize(n_compute, matrix_batch)
    integer, intent(in) :: n_compute
    real(dp), intent(in out) :: matrix_batch(n_compute,n_compute)
    integer :: nb
    do nb = 1, n_compute-1
       matrix_batch(nb+1:n_compute,nb) = (matrix_batch(nb+1:n_compute,nb) + &
            & matrix_batch(nb,nb+1:n_compute))/2
    end do
  end subroutine symmetrize

  pure subroutine antisymmetrize(n_compute, matrix_batch)
    integer, intent(in) :: n_compute
    real(dp), intent(in out) :: matrix_batch(n_compute,n_compute)
    integer :: nb
    do nb = 1, n_compute-1
       matrix_batch(nb+1:n_compute,nb) = (matrix_batch(nb+1:n_compute,nb) - &
            & matrix_batch(nb,nb+1:n_compute))/2
    end do
    do nb = 1, n_compute
       matrix_batch(nb,nb) = 0d0
    end do
  end subroutine antisymmetrize

  subroutine cleanup_integration()
    integer :: i_batch
    if (allocated(nonzero_indices)) then
       do i_batch = 1, size(nonzero_indices)
          call safe_deallocate(nonzero_indices(i_batch)%arr)
       end do
       deallocate(nonzero_indices)
    end if
    call cleanup_xc_kernel_rho()
    call cleanup_xc_wrapper()
  end subroutine cleanup_integration
end module integration
