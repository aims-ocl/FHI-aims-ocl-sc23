!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Wrappers routines for cases where the code needs to decide whether
!!  to use LibXC or call internal routines such as evaluate_xc.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module xc_wrapper

  use aims_memory_tracking, only: aims_allocate
  use dimensions,           only: n_max_batch_size, n_spin, use_gga
  use mpi_tasks,            only: aims_stop
  use runtime_choices,      only: flag_xc, spin_treatment
  use tools,                only: libxc_yes_no, norm2, safe_deallocate, str
  use types,                only: dp
  use xc_f03_lib_m

  implicit none

  private
  public :: initialize_libxc, initialize_xc_wrapper, evaluate_xc_wrapper, &
       & cleanup_xc_wrapper

  ! LibXC exchange and correlation constructs
  type(xc_f03_func_t), public :: xc_func_x, xc_func_c
  real(dp), allocatable :: vxc_wrk(:,:)
  real(dp), allocatable :: vsigma_wrk(:,:)
  real(dp), allocatable :: sigma(:,:)
  ! If a LibXC functional has been initialized (LibXC doesn't have a
  ! built-in function for checking this).
  logical, public :: libxc_func_initialized = .false.
  ! Whether the exchange and correlation parts are initialized
  ! separately. For example, in case of B3LYP xc_func_c refers to both
  ! exchange and correlation, i.e. the full xc-functional.
  logical, public :: is_full_xc

  logical, public :: contains_lyp

contains

  subroutine initialize_libxc()
    character(*), parameter :: &
         & THIS_SUB = 'xc_wrapper::initialize_libxc::'
    contains_lyp = .false.
    if (libxc_yes_no()) then
       ! See evaluate_xc.f90 for deciphering the flags
       select case (flag_xc)
       case(3) ! PZ81
          call xc_f03_func_init(xc_func_x, XC_LDA_X, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_LDA_C_PZ, XC_POLARIZED)
       case(4) ! PW91
          call xc_f03_func_init(xc_func_x, XC_GGA_X_PW91, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_GGA_C_PW91, XC_POLARIZED)
       case(1, 6) ! PBE, PBE0
          call xc_f03_func_init(xc_func_x, XC_GGA_X_PBE, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_GGA_C_PBE, XC_POLARIZED)
       case(8) ! PW92
          call xc_f03_func_init(xc_func_x, XC_LDA_X, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_LDA_C_PW, XC_POLARIZED)
       case(9) ! BLYP
          call xc_f03_func_init(xc_func_x, XC_GGA_X_B88, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_GGA_C_LYP, XC_POLARIZED)
          contains_lyp = .true.
       case(10) ! B3LYP
          call xc_f03_func_init(xc_func_c, XC_HYB_GGA_XC_B3LYP, XC_POLARIZED)
          contains_lyp = .true.
       case(11) ! RPBE
          call xc_f03_func_init(xc_func_x, XC_GGA_X_RPBE, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_GGA_C_PBE, XC_POLARIZED)
       case(15) ! VWN5
          call xc_f03_func_init(xc_func_x, XC_LDA_X, XC_POLARIZED)
          call xc_f03_func_init(xc_func_c, XC_LDA_C_VWN, XC_POLARIZED)
       case default
          call aims_stop('Requested functional not available in libxc ('// &
               & str(flag_xc)//').', THIS_SUB)
       end select
       libxc_func_initialized = .true.
       is_full_xc = &
            & xc_f03_func_info_get_kind(xc_f03_func_get_info(xc_func_c)) == &
            & XC_EXCHANGE_CORRELATION
    else
       libxc_func_initialized = .false.
    end if
  end subroutine initialize_libxc

  subroutine initialize_xc_wrapper()
    character(*), parameter :: &
         & THIS_SUB = 'xc_wrapper::initialize_xc_wrapper::'
    if (allocated(vxc_wrk)) return
    if (.not. libxc_func_initialized) call initialize_libxc()
    call aims_allocate(vxc_wrk, 2, n_max_batch_size, THIS_SUB//'vxc_wrk')
    if (use_gga) then
       call aims_allocate(sigma, 3, n_max_batch_size, THIS_SUB//'sigma')
       call aims_allocate(vsigma_wrk, 3, n_max_batch_size, &
            & THIS_SUB//'vsigma_wrk')
    end if
  end subroutine initialize_xc_wrapper

  subroutine evaluate_xc_wrapper(n_points, local_rho, local_rho_gradient, vxc, &
       & vsigma)
    integer, intent(in) :: n_points
    real(dp), intent(in) :: local_rho(2,n_points)
    real(dp), intent(in) :: local_rho_gradient(3,2,n_points)
    real(dp), intent(out) :: vxc(2,n_points), vsigma(3,n_points)
    real(dp) :: dummy(6)
    integer :: n_spin_save, spin_treatment_save, ip
    if (libxc_func_initialized) then
       if (.not. use_gga) then
          call xc_f03_lda_vxc(xc_func_x, n_points, local_rho, vxc_wrk)
          call xc_f03_lda_vxc(xc_func_c, n_points, local_rho, vxc)
          vxc = vxc + vxc_wrk(:,:n_points)
       else
          sigma(1,:n_points) = norm2(local_rho_gradient(:,1,:),1)**2
          sigma(2,:n_points) = &
               & sum(local_rho_gradient(:,1,:)*local_rho_gradient(:,2,:),1)
          sigma(3,:n_points) = norm2(local_rho_gradient(:,2,:),1)**2
          call xc_f03_gga_vxc(xc_func_x, n_points, local_rho, sigma, vxc_wrk, &
               & vsigma_wrk)
          call xc_f03_gga_vxc(xc_func_c, n_points, local_rho, sigma, vxc, &
               & vsigma)
          vxc = vxc + vxc_wrk(:,:n_points)
          vsigma = vsigma + vsigma_wrk(:,:n_points)
       end if
    else
       n_spin_save = n_spin
       n_spin = 2
       spin_treatment_save = spin_treatment
       spin_treatment = 1
       do ip = 1, n_points
          call evaluate_xc(local_rho(:,ip), local_rho_gradient(:,:,ip), 0d0, &
               & dummy, dummy, dummy, vxc(:,ip), dummy, dummy, .false.)
       end do
       n_spin = n_spin_save
       spin_treatment = spin_treatment_save
    end if
  end subroutine evaluate_xc_wrapper

  subroutine cleanup_xc_wrapper()
    call safe_deallocate(vxc_wrk)
    call safe_deallocate(vsigma_wrk)
    call safe_deallocate(sigma)
    if (libxc_func_initialized) then ! If initialize_libxc has been called
       if (.not. is_full_xc) call xc_f03_func_end(xc_func_x)
       call xc_f03_func_end(xc_func_c)
       libxc_func_initialized = .false.
    end if
  end subroutine cleanup_xc_wrapper
end module xc_wrapper
