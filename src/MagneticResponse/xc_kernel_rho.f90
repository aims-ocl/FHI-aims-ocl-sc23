!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Evaluates the first-order xc-potential (xc-kernel times
!!  first-order density) on a batch of grid points. The first-order
!!  potential consists of a scalar and a vectorial part:
!!
!!  xc_kernel_s[1] = v2rho2[1]n1[1] + v2rho2[2]n1[2]
!!                   + Sum_i=1^3 v2rhosigma[i]sigma1[i],
!!  xc_kernel_s[2] = v2rho2[3]n1[2] + v2rho2[2]n1[1]
!!                   + Sum_i=1^3 v2rhosigma[i+3]sigma1[i],
!!  xc_kernel_v[1] = 2 v2rhosigma[1] Nabla n[1] n1[1]
!!                 + 2 v2rhosigma[4] Nabla n[1] n1[2]
!!                 +   v2rhosigma[2] Nabla n[2] n1[1]
!!                 +   v2rhosigma[5] Nabla n[2] n1[2]
!!                 + 2 v2sigma2[1] Nabla n[1] sigma1[1]
!!                 + 2 v2sigma2[2] Nabla n[1] sigma1[2]
!!                 + 2 v2sigma2[3] Nabla n[1] sigma1[3]
!!                 +   v2sigma2[2] Nabla n[2] sigma1[1]
!!                 +   v2sigma2[4] Nabla n[2] sigma1[2]
!!                 +   v2sigma2[5] Nabla n[2] sigma1[3]
!!                 + 2 vsigma[1] Nabla n1[1] + vsigma[2] Nabla n1[2],
!!  xc_kernel_v[2] = 2 v2rhosigma[3] Nabla n[2] n1[1]
!!                 + 2 v2rhosigma[6] Nabla n[2] n1[2]
!!                 +   v2rhosigma[2] Nabla n[1] n1[1]
!!                 +   v2rhosigma[5] Nabla n[1] n1[2]
!!                 + 2 v2sigma2[3] Nabla n[2] sigma1[1]
!!                 + 2 v2sigma2[5] Nabla n[2] sigma1[2]
!!                 + 2 v2sigma2[6] Nabla n[2] sigma1[3]
!!                 +   v2sigma2[2] Nabla n[1] sigma1[1]
!!                 +   v2sigma2[4] Nabla n[1] sigma1[2]
!!                 +   v2sigma2[5] Nabla n[1] sigma1[3]
!!                 + 2 vsigma[3] Nabla n1[2] + vsigma[2] Nabla n1[1],
!!
!!  where we have used the LibXC definitions:
!!    sigma = [|Nabla n[1]|^2, Nabla n[1]Nabla n[2], |Nabla n[2]|^2]
!!    v2rho2 = \partial^2 Exc/(\partial n \partial n)
!!    vsigma = \partial Exc/\partial sigma
!!    v2rhosigma = \partial^2 Exc/(\partial n \partial sigma)
!!    v2sigma2 = \partial^2 Exc/(\partial sigma \partial sigma)
!!
!!  The indexing scheme for the latter two is
!!
!!  v2rhosigma:     v2sigma2:
!!  [1] -> 1,11   [1] -> 11,11
!!  [2] -> 1,12   [2] -> 11,12
!!  [3] -> 1,22   [3] -> 11,22
!!  [4] -> 2,11   [4] -> 12,12
!!  [5] -> 2,12   [5] -> 12,22
!!  [6] -> 2,23   [6] -> 22,22
!!
!!  The scalar and vectorial parts enter the matrix elements of the
!!  first-order response of the xc-potential as follows:
!!
!!  <i| Vxc1[alpha] |j> = Int xc_kernel_s[alpha](r) Omega_ij(r) D r
!!                    + Int xc_kernel_v[alpha](r) Nabla Omega_ij(r) D r,
!!          Omega_ij(r) = wave_i(r)*wave_j(r),
!!    Nabla Omega_ij(r) = 2 wave_i(r)*Nabla wave_j(r).
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
!!  COMMENTS
!!
!!  The below subroutines can process more than one perturbation
!!  direction simultaneously, i.e., rho1 (first-order density) can
!!  have more than one component. Currently the maximum number of
!!  perturbation directions that are processed at the same time is
!!  limited to 3 - this is built into the allocate statements of
!!  initialize_xc_kernel_rho and could easily be changed.
!!
module xc_kernel_rho

  use aims_memory_tracking, only: aims_allocate
  use dimensions,           only: n_max_batch_size, n_spin, use_gga
  use load_balancing,       only: batch_perm, n_bp_integ
  use physics,              only: rho, rho_gradient
  use runtime_choices,      only: hybrid_coeff
  use tools,                only: libxc_yes_no, norm2, safe_deallocate
  use types,                only: dp
  use xc_f03_lib_m
  use xc_wrapper,           only: contains_lyp, initialize_libxc, is_full_xc, &
       & xc_func_x, xc_func_c

  implicit none

  private
  public :: cleanup_xc_kernel_rho, initialize_xc_kernel_rho, &
       & is_first_xc_kernel_run, xc_kernel_rho_LDA, xc_kernel_rho_GGA

  ! Scalar part of thef first-order potential
  real(dp), allocatable, public :: xc_kernel_s(:,:,:)
  ! Vectorial part of the first-order potential
  real(dp), allocatable, public :: xc_kernel_v(:,:,:,:)
  ! Second derivatives of Exc. These need to be evaluated only once,
  ! as they are the same for each perturbation and each DFPT step as
  ! long as the unperturbed density doesn't change.
  real(dp), allocatable :: sigma(:,:), v2rho2(:,:), vsigma(:,:), &
       & v2rhosigma(:,:), v2sigma2(:,:)
  ! Auxiliary work arrays for the libxc subroutines. Last index loops
  ! over exchange and correlation parts.
  real(dp), allocatable :: v2rho2_tmp(:,:,:), vsigma_tmp(:,:,:), &
       & v2rhosigma_tmp(:,:,:), v2sigma2_tmp(:,:,:)
  ! First-order sigma (in literature this is often gamma).
  real(dp), allocatable :: sigma1(:,:)

contains

  !!  FUNCTION
  !!
  !!  Lets the caller know whether the second derivatives should be
  !!  evaluated on this run.
  !!
  pure logical function is_first_xc_kernel_run() result(y)
    y = .not. allocated(v2rho2)
  end function is_first_xc_kernel_run

  subroutine initialize_xc_kernel_rho()
    character(*), parameter :: &
         & THIS_SUB = 'xc_kernel_rho::initialize_xc_kernel_rho::'
    ! Determines the size of some work arrays depending on whether
    ! xc_func_c contains both exchange and correlation or just the
    ! correlation part.
    integer :: i_size
    call initialize_libxc()
    if (is_full_xc) then
       i_size = 1
    else
       i_size = 2
    end if
    call aims_allocate(v2rho2, 3, batch_perm(n_bp_integ)%n_full_points, &
         & THIS_SUB//'v2rho2')
    call aims_allocate(v2rho2_tmp, 3, n_max_batch_size, i_size, &
         & THIS_SUB//'v2rho2_tmp')
    call aims_allocate(xc_kernel_s, n_spin, n_max_batch_size, 3, &
         & THIS_SUB//'xc_kernel_s')
    if (use_gga) then
       call aims_allocate(sigma, 3, n_max_batch_size, THIS_SUB//'sigma')
       call aims_allocate(sigma1, 3, n_max_batch_size, THIS_SUB//'sigma1')
       call aims_allocate(vsigma, 3, batch_perm(n_bp_integ)%n_full_points, &
            & THIS_SUB//'vsigma')
       call aims_allocate(v2rhosigma, 6, batch_perm(n_bp_integ)%n_full_points, &
            & THIS_SUB//'v2rhosigma')
       call aims_allocate(v2sigma2, 6, batch_perm(n_bp_integ)%n_full_points, &
            & THIS_SUB//'v2sigma2')
       call aims_allocate(vsigma_tmp, 3, n_max_batch_size, i_size, &
            & THIS_SUB//'vsigma_tmp')
       call aims_allocate(v2rhosigma_tmp, 6, n_max_batch_size, i_size, &
            & THIS_SUB//'v2rhosigma_tmp')
       call aims_allocate(v2sigma2_tmp, 6, n_max_batch_size, i_size, &
            & THIS_SUB//'v2sigma2_tmp')
       call aims_allocate(xc_kernel_v, 3, n_spin, n_max_batch_size, 3, &
            & THIS_SUB//'xc_kernel_v')
    end if
  end subroutine initialize_xc_kernel_rho

  !!  FUNCTION
  !!
  !!  Computes the xc-kernel and multiplies with the first-order density.
  !!
  !!  OUTPUT
  !!
  !!  xc_kernel_s - scalar part of the first-order xc-potential
  !!
  subroutine xc_kernel_rho_LDA(xc_kernel_rho_first_run, n_perts, n_points, &
       & nonzero_indices, rho1, rho, local_rho)
    ! Note that the kernel needs to be computed only once.
    logical, intent(in) :: xc_kernel_rho_first_run
    ! Number of perturbations (different directions of the same
    ! perturbation or just different perturbations) that are processed
    ! by this call
    integer, intent(in) :: n_perts
    ! Number of points in the current batch
    integer, intent(in) :: n_points
    ! Indices of nonzero elements in the global arrays
    integer, intent(in) :: nonzero_indices(:)
    ! First-order density. Third index runs over n_perts
    real(dp), intent(in) :: rho1(:,:,:)
    ! Electron density
    real(dp), intent(in) :: rho(:,:)
    ! Work array for the density
    real(dp), intent(out) :: local_rho(:,:)
    integer :: i_pert
    if (.not. allocated(v2rho2)) call initialize_xc_kernel_rho()
    if (.not. libxc_yes_no()) then
       call xc_kernel_rho_LDA_no_libxc(n_perts, n_points, &
            & rho(:,nonzero_indices), rho1(:,nonzero_indices,:))
       return
    end if
    first_run: if (xc_kernel_rho_first_run) then
       local_rho(:,:n_points) = rho(:,nonzero_indices)
       call xc_f03_lda_fxc(xc_func_x, n_points, local_rho, v2rho2_tmp(:,:,1))
       call xc_f03_lda_fxc(xc_func_c, n_points, local_rho, v2rho2_tmp(:,:,2))
       v2rho2_tmp(:,:n_points,1) = (1d0-hybrid_coeff)*v2rho2_tmp(:,:n_points,1)
       v2rho2(:,nonzero_indices) = sum(v2rho2_tmp(:,:n_points,:),3)
    end if first_run
    do i_pert = 1, n_perts
       xc_kernel_s(1,:n_points,i_pert) = &
            & sum(v2rho2(:2,nonzero_indices)*rho1(:,nonzero_indices,i_pert),1)
       if (n_spin == 2) &
            & xc_kernel_s(2,:n_points,i_pert) = &
            & sum(v2rho2(2:,nonzero_indices)*rho1(:,nonzero_indices,i_pert),1)
    end do
  end subroutine xc_kernel_rho_LDA

  !!  FUNCTION
  !!
  !!  Compute the two parts of the xc-kernel and multiplies with the
  !!  first-order density and the first-order density gradient.
  !!
  !!  OUTPUT
  !!
  !!  xc_kernel_s - scalar part of the first-order xc-potential
  !!  xc_kernel_v - vectorial part of the first-order xc-potential
  !!
  subroutine xc_kernel_rho_GGA(xc_kernel_rho_first_run, n_perts, n_points, &
       & nonzero_indices, rho1, rho1_gradient, rho, rho_gradient_p, &
       & local_rho, local_rho_gradient)
    ! Arguments are analogous to xc_kernel_rho_LDA
    logical, intent(in) :: xc_kernel_rho_first_run
    integer, intent(in) :: n_perts, n_points, nonzero_indices(:)
    real(dp), intent(in) :: rho1(:,:,:), rho1_gradient(:,:,:,:)
    ! Electron density and density gradient
    real(dp), intent(in) :: rho(:,:), rho_gradient_p(:,:,:)
    real(dp), intent(out) :: local_rho(:,:), local_rho_gradient(:,:,:)
    integer :: i_pert, i_dir
    if (.not. allocated(v2rho2)) call initialize_xc_kernel_rho()
    local_rho_gradient(:,:,:n_points) = rho_gradient_p(:,:,nonzero_indices)
    first_run: if (xc_kernel_rho_first_run) then
       local_rho(:,:n_points) = rho(:,nonzero_indices)
       sigma(1,:n_points) = norm2(local_rho_gradient(:,1,:n_points),1)**2
       sigma(2,:n_points) = sum(local_rho_gradient(:,1,:n_points)* &
            & local_rho_gradient(:,2,:n_points),1)
       sigma(3,:n_points) = norm2(local_rho_gradient(:,2,:n_points),1)**2
       if (contains_lyp) where (local_rho(:,:n_points) < 0.2374d-10) &
            & local_rho(:,:n_points) = 0d0
       if (is_full_xc) then
          ! xc_func_c refers to both exchange and correlation (e.g., B3LYP).
          if (contains_lyp) where (local_rho(:,:n_points) < tiny(1d0)) &
               & local_rho(:,:n_points) = 1d-15
          call xc_f03_gga_vxc(xc_func_c, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,1), vsigma_tmp(:,:,1))
          call xc_f03_gga_fxc(xc_func_c, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,1), v2rhosigma_tmp(:,:,1), v2sigma2_tmp(:,:,1))
          v2rho2(:,nonzero_indices) = v2rho2_tmp(:,:n_points,1)
          v2rhosigma(:,nonzero_indices) = v2rhosigma_tmp(:,:n_points,1)
          v2sigma2(:,nonzero_indices) = v2sigma2_tmp(:,:n_points,1)
          vsigma(:,nonzero_indices) = vsigma_tmp(:,:n_points,1)
       else
          call xc_f03_gga_vxc(xc_func_x, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,1), vsigma_tmp(:,:,1))
          call xc_f03_gga_vxc(xc_func_c, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,2), vsigma_tmp(:,:,2))
          call xc_f03_gga_fxc(xc_func_x, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,1), v2rhosigma_tmp(:,:,1), v2sigma2_tmp(:,:,1))
          call xc_f03_gga_fxc(xc_func_c, n_points, local_rho, sigma, &
               & v2rho2_tmp(:,:,2), v2rhosigma_tmp(:,:,2), v2sigma2_tmp(:,:,2))
          v2rho2_tmp(:,:n_points,1) = &
               & (1d0-hybrid_coeff)*v2rho2_tmp(:,:n_points,1)
          v2rhosigma_tmp(:,:n_points,1) = &
               & (1d0-hybrid_coeff)*v2rhosigma_tmp(:,:n_points,1)
          v2sigma2_tmp(:,:n_points,1) = &
               & (1d0-hybrid_coeff)*v2sigma2_tmp(:,:n_points,1)
          vsigma_tmp(:,:n_points,1) = &
               & (1d0-hybrid_coeff)*vsigma_tmp(:,:n_points,1)
          v2rho2(:,nonzero_indices) = sum(v2rho2_tmp(:,:n_points,:),3)
          v2rhosigma(:,nonzero_indices) = sum(v2rhosigma_tmp(:,:n_points,:),3)
          v2sigma2(:,nonzero_indices) = sum(v2sigma2_tmp(:,:n_points,:),3)
          vsigma(:,nonzero_indices) = sum(vsigma_tmp(:,:n_points,:),3)
       end if
    end if first_run
    over_perturbations: do i_pert = 1, n_perts
       ! First-order sigma
       sigma1(1,:n_points) = 2d0*sum(local_rho_gradient(:,1,:n_points)* &
            & rho1_gradient(:,1,nonzero_indices,i_pert),1)
       sigma1(2,:n_points) = sum(local_rho_gradient(:,1,:n_points)* &
            & rho1_gradient(:,2,nonzero_indices,i_pert),1) + &
            & sum(local_rho_gradient(:,2,:n_points)* &
            & rho1_gradient(:,1,nonzero_indices,i_pert),1)
       sigma1(3,:n_points) = 2d0*sum(local_rho_gradient(:,2,:n_points)* &
            & rho1_gradient(:,2,nonzero_indices,i_pert),1)
       ! Scalar part
       ! v2rho2 terms
       xc_kernel_s(1,:n_points,i_pert) = &
            & sum(v2rho2(:2,nonzero_indices)*rho1(:,nonzero_indices,i_pert),1)
       ! v2rhosigma terms
       xc_kernel_s(1,:n_points,i_pert) = xc_kernel_s(1,:n_points,i_pert) + &
            & sum(v2rhosigma(:3,nonzero_indices)*sigma1(:,:n_points),1)
       ! Vectorial part
       do i_dir = 1, 3
          ! v2rhosigma terms
          xc_kernel_v(i_dir,1,:n_points,i_pert) = &
               & 2d0*local_rho_gradient(i_dir,1,:n_points)* &
               & sum(v2rhosigma(1::3,nonzero_indices)* &
               & rho1(:,nonzero_indices,i_pert),1) + &
               & local_rho_gradient(i_dir,2,:n_points)* &
               & sum(v2rhosigma(2::3,nonzero_indices)* &
               & rho1(:,nonzero_indices,i_pert),1)
          ! v2sigma2 terms
          xc_kernel_v(i_dir,1,:n_points,i_pert) = &
               & xc_kernel_v(i_dir,1,:n_points,i_pert) + &
               & 2d0*local_rho_gradient(i_dir,1,:n_points)* &
               & sum(v2sigma2(:3,nonzero_indices)*sigma1(:,:n_points),1) + &
               & local_rho_gradient(i_dir,2,:n_points)* &
               & (v2sigma2(2,nonzero_indices)*sigma1(1,:n_points) + &
               & sum(v2sigma2(4:5,nonzero_indices)*&
               & sigma1(2:3,:n_points),1))
          ! vsigma terms
          xc_kernel_v(i_dir,1,:n_points,i_pert) = &
               & xc_kernel_v(i_dir,1,:n_points,i_pert) + &
               & 2d0*vsigma(1,nonzero_indices)* &
               & rho1_gradient(i_dir,1,nonzero_indices,i_pert) + &
               & vsigma(2,nonzero_indices)* &
               & rho1_gradient(i_dir,2,nonzero_indices,i_pert)
       end do
       if (n_spin == 2) then
          ! Scalar part
          ! v2rho2 terms
          xc_kernel_s(2,:n_points,i_pert) = sum(v2rho2(2:,nonzero_indices)* &
               & rho1(:,nonzero_indices,i_pert),1)
          ! v2rhosigma terms
          xc_kernel_s(2,:n_points,i_pert) = xc_kernel_s(2,:n_points,i_pert) + &
               & sum(v2rhosigma(4:,nonzero_indices)*sigma1(:,:n_points),1)
          ! Vectorial part
          do i_dir = 1, 3
             ! v2rhosigma terms
             xc_kernel_v(i_dir,2,:n_points,i_pert) = &
                  & 2d0*local_rho_gradient(i_dir,2,:n_points)* &
                  & sum(v2rhosigma(3::3,nonzero_indices)* &
                  & rho1(:,nonzero_indices,i_pert),1) + &
                  & local_rho_gradient(i_dir,1,:n_points)* &
                  & sum(v2rhosigma(2::3,nonzero_indices)* &
                  & rho1(:,nonzero_indices,i_pert),1)
             ! v2sigma2 terms
             xc_kernel_v(i_dir,2,:n_points,i_pert) = &
                  & xc_kernel_v(i_dir,2,:n_points,i_pert) + &
                  & 2d0*local_rho_gradient(i_dir,2,:n_points)* &
                  & (v2sigma2(3,nonzero_indices)*sigma1(1,:n_points) + &
                  & sum(v2sigma2(5:6,nonzero_indices)* &
                  & sigma1(2:3,:n_points),1)) + &
                  & local_rho_gradient(i_dir,1,:n_points)* &
                  & (v2sigma2(2,nonzero_indices)*sigma1(1,:n_points)*&
                  & sum(v2sigma2(4:5,nonzero_indices)*sigma1(2:3,:n_points),1))
             ! vsigma terms
             xc_kernel_v(i_dir,2,:n_points,i_pert) = &
                  & xc_kernel_v(i_dir,2,:n_points,i_pert) + &
                  & 2d0*vsigma(3,nonzero_indices)* &
                  & rho1_gradient(i_dir,2,nonzero_indices,i_pert) + &
                  & vsigma(2,nonzero_indices)* &
                  & rho1_gradient(i_dir,1,nonzero_indices,i_pert)
          end do
       end if
    end do over_perturbations
  end subroutine xc_kernel_rho_GGA

  !!  FUNCTION
  !!
  !!  Simple finite differences implementation of the LDA kernel in
  !!  case LibXC is not available.
  !!
  subroutine xc_kernel_rho_LDA_no_libxc(n_perts, n_points, local_rho, &
       & local_rho1)
    use dimensions,      only: n_spin
    use runtime_choices, only: spin_treatment
    integer, intent(in) :: n_perts, n_points
    real(dp), intent(in) :: local_rho(:,:), local_rho1(:,:,:)
    real(dp), parameter :: epsilon = 5d-5
    real(dp) :: arg(2), tmp(2,2)
    integer :: n_spin_save, spin_treatment_save, i_pert, ip
    n_spin_save = n_spin
    n_spin = 2
    spin_treatment_save = spin_treatment
    spin_treatment = 1
    do i_pert = 1, n_perts
       do ip = 1, n_points
          arg = local_rho(:,ip) + [epsilon*local_rho1(1,ip,i_pert), 0d0]
          call evaluate_xc_loc(arg, tmp(:,1))
          arg = local_rho(:,ip) - [epsilon*local_rho1(1,ip,i_pert), 0d0]
          call evaluate_xc_loc(arg, tmp(:,2))
          xc_kernel_s(1,ip,i_pert) = tmp(1,1) - tmp(1,2)
          arg = local_rho(:,ip) + [epsilon*local_rho1(2,ip,i_pert), 0d0]
          call evaluate_xc_loc(arg, tmp(:,1))
          arg = local_rho(:,ip) - [epsilon*local_rho1(2,ip,i_pert), 0d0]
          call evaluate_xc_loc(arg, tmp(:,2))
          xc_kernel_s(1,ip,i_pert) = &
               & xc_kernel_s(1,ip,i_pert) + tmp(2,1) - tmp(2,2)
          if (n_spin_save == 2) then
             arg = local_rho(:,ip) + [0d0, epsilon*local_rho1(2,ip,i_pert)]
             call evaluate_xc_loc(arg, tmp(:,1))
             arg = local_rho(:,ip) - [0d0, epsilon*local_rho1(2,ip,i_pert)]
             call evaluate_xc_loc(arg, tmp(:,2))
             xc_kernel_s(2,ip,i_pert) = tmp(2,1) - tmp(2,2)
             arg = local_rho(:,ip) + [0d0, epsilon*local_rho1(1,ip,i_pert)]
             call evaluate_xc_loc(arg, tmp(:,1))
             arg = local_rho(:,ip) - [0d0, epsilon*local_rho1(1,ip,i_pert)]
             call evaluate_xc_loc(arg, tmp(:,2))
             xc_kernel_s(2,ip,i_pert) = xc_kernel_s(2,ip,i_pert) + &
                  & tmp(1,1) - tmp(1,2)
          end if
       end do
    end do
    xc_kernel_s(:,:n_points,:n_perts) = &
         & xc_kernel_s(:,:n_points,:n_perts)/2/epsilon
    n_spin = n_spin_save
    spin_treatment = spin_treatment
  contains
    subroutine evaluate_xc_loc(rho, local_xc_derivs)
      real(dp), intent(in) :: rho(2), local_xc_derivs(2)
      real(dp) :: dummy(30)
      call evaluate_xc(rho, 0d0, 0d0, dummy, dummy, dummy, local_xc_derivs, &
           & dummy, dummy, .false.)
    end subroutine evaluate_xc_loc
  end subroutine xc_kernel_rho_LDA_no_libxc

  subroutine cleanup_xc_kernel_rho()
    call safe_deallocate(xc_kernel_s)
    call safe_deallocate(xc_kernel_v)
    call safe_deallocate(sigma)
    call safe_deallocate(sigma1)
    call safe_deallocate(vsigma)
    call safe_deallocate(vsigma_tmp)
    call safe_deallocate(v2rho2)
    call safe_deallocate(v2rho2_tmp)
    call safe_deallocate(v2rhosigma)
    call safe_deallocate(v2rhosigma_tmp)
    call safe_deallocate(v2sigma2)
    call safe_deallocate(v2sigma2_tmp)
  end subroutine cleanup_xc_kernel_rho
end module xc_kernel_rho
