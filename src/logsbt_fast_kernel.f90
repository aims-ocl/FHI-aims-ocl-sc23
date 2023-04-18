!****h* FHI-aims/logsbt_fast_kernel
!  NAME
!    logsbt_fast_kernel
!  SYNOPSIS

module logsbt_fast_kernel

  !  PURPOSE
  !
  !    Provide a splined logSBT kernel for fast access for arbitrary offsets.
  !
  !    The logSBT kernels are nothing else than (well thought-of)
  !    approximations to spherical Bessel functions on a logarithmic grid of
  !    k*r values.  While the number of grid points in r or k and their
  !    spacing in generally stays the same, the offset depends on specific
  !    values of k or r.  Reconstructing them again and again is not
  !    particularly expensive but may be a bottleneck in some circumstances.
  !    Therefore, this module provides a splined version of the kernels which
  !    can fastly be evaluated for any offset.
  !
  !  USAGE
  !
  !    Initialize, get, cleanup.
  !
  !  USES

  use logsbt
  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  ! --- Splined kernel parameters

  ! Use four times more spline grid points than Fourier grid points for kernel
  integer, parameter :: sbt_spl_multiplier = 4

  integer, private :: sbt_Nt_spl       ! Number of Fourier components
  integer, private :: sbt_N_spl        ! Number of spline knots
  integer, private :: sbt_L_spl        ! Number of angular momenta
  real*8, private :: sbt_lnt0_spl      ! Onset
  real*8, private :: sbt_lntrange_spl  ! Range
  real*8, private :: sbt_power_spl     ! Power bias [0.d0]
  real*8, allocatable, private :: sbt_KK_spl(:,:,:)  !(4,sbt_N_spl,0:sbt_L_spl)

  integer :: sbt_n_spl_ref = 0

contains

  !----------------------------------------------------------------------------
  !****s* logsbt_fast_kernel/initialize_fast_kernel
  !  NAME
  !    initialize_fast_kernel
  !  SYNOPSIS

  subroutine initialize_fast_kernel(N, max_L, lnr0, lnk0, lnrange, power)

    !  PURPOSE
    !
    !    Initialize the module-wide splined kernel.  In most cases,
    !    power=0.d0.
    !
    !    The spline is prepared for Nt=2*N log-spacedgrid points in the
    !    interval [exp(lnr0+lnk0), exp(lnr0+lnk0+2*lnrange)[.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N, max_L
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: power

    !  INPUTS
    !    o N -- Number of Fourier components of loggrid
    !    o max_L -- Highest angular momentum the kernel is needed for
    !    o lnr0, lnk0, lnrange -- Grid ranges
    !    o power -- Power bias of kernel
    !  OUTPUTS
    !    module data
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: Nt
    real*8 :: lnt0, lntrange
    integer :: info
    character(*), parameter :: func = 'initialize_fast_kernel'

    Nt = 2*N
    lnt0 = lnr0 + lnk0
    lntrange = 2*lnrange

    if (allocated(sbt_KK_spl)) then
       if (sbt_Nt_spl /= Nt .or. &
       &   sbt_N_spl /= sbt_spl_multiplier * Nt .or. &
       &   abs(sbt_power_spl - power) > 1d-10 .or. &
       &   abs(sbt_lnt0_spl - lnt0) > 1d-10 .or. &
       &   abs(sbt_lntrange_spl - lntrange) > 1d-10) then
          call aims_stop('incompatible double initialization', func)
       else if (sbt_L_spl >= max_L) then
          return
       else
          deallocate(sbt_KK_spl)
       end if
    else
       sbt_Nt_spl = Nt
       sbt_N_spl = sbt_spl_multiplier * Nt
       sbt_power_spl = power
       sbt_lnt0_spl = lnt0
       sbt_lntrange_spl = lntrange
    end if
    sbt_L_spl = max_L
    allocate(sbt_KK_spl(4, sbt_N_spl, 0:sbt_L_spl), stat=info)
    call check_allocation(info, 'sbt_KK_spl', func)
    call sbt_splined_kernel(sbt_N_spl, sbt_Nt_spl, &
    &                       sbt_lnt0_spl, sbt_lntrange_spl,&
    &                       sbt_L_spl, power, sbt_KK_spl)

    sbt_n_spl_ref = sbt_n_spl_ref + 1

  end subroutine initialize_fast_kernel
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt_fast_kernel/cleanup_fast_kernel
  !  NAME
  !    cleanup_fast_kernel
  !  SYNOPSIS

  subroutine cleanup_fast_kernel()

    !  PURPOSE
    !    Clean up splined kernel.  Frees resources and allows another call to
    !    initialize_fast_kernel().
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'cleanup_fast_kernel'

    sbt_n_spl_ref = sbt_n_spl_ref - 1
    if (sbt_n_spl_ref < 0) then
       call aims_stop('Tried to clean up clean splined kernel', func)
    else if (sbt_n_spl_ref == 0) then
       deallocate(sbt_KK_spl)
    end if

  end subroutine cleanup_fast_kernel
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt_fast_kernel/get_fast_kernel
  !  NAME
  !    get_fast_kernel
  !  SYNOPSIS

  subroutine get_fast_kernel(N, L_min, L_max, lnt0, lnrange, power, &
  &                          use_spline, KK)

    !  PURPOSE
    !
    !    Evaluate kernel from spline prepared by intialize_fast_kernel().
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N, L_min, L_max
    real*8, intent(IN) :: lnt0, lnrange
    real*8, intent(IN) :: power
    logical, intent(IN) :: use_spline
    real*8, intent(OUT) :: KK(N, L_min:L_max)

    !  INPUTS
    !    o N -- Number of grid points [size(KK, 1)]
    !    o L -- Angular momentum
    !    o lnt0, lnrange -- Range of kernel to retrieve
    !                       Must be /within/ the corresponding arguments
    !                       to initialize_fast_kernel().
    !    o power -- Power bias [0.d0]
    !    o use_spline -- If .false., completely reconstruct kernel
    !  OUTPUTS
    !    o KK -- Kernel
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: L, i_R_0, i, i_spl
    real*8 :: Rho
    real*8 :: d_R_0, d_R_0_sq, d_R_0_cb
    real*8, allocatable :: wsave(:)
    integer :: info
    character(*), parameter :: func = 'get_fast_kernel'

    if (use_spline) then
       if (L_max > sbt_L_spl) call aims_stop('L_max too large', func)
       Rho = lnt0 - sbt_lnt0_spl
       if (Rho < 0.d0 .or. Rho+lnrange > sbt_lntrange_spl) then
          call aims_stop('Invalid |R|', func)
       end if

       ! (lnrage / N) is more coarse (larger step length) than
       ! (sbt_lntrange_spl / sbt_N_spl).
       if (abs(sbt_spl_multiplier * (sbt_lntrange_spl / sbt_N_spl) &
       &       - (lnrange / N)) > 1d-10) then
          call aims_stop('N is not multiple of sbt_N_spl', func)
       end if
       if (abs(power - sbt_power_spl) > 1d-10) then
          call aims_stop('Initialized with different power', func)
       end if

       ! --- Spline parameters

       ! invert_log_grid(): i(r) = 1 + [ln(r / r_grid_min)]/[ln(r_grid_inc)]
       d_R_0 = 1.d0 + Rho / (sbt_lntrange_spl / sbt_N_spl)
       i_R_0 = floor(d_R_0)
       d_R_0 = d_R_0 - i_R_0
       d_R_0_sq = d_R_0 * d_R_0
       d_R_0_cb = d_R_0_sq * d_R_0
       ! Security check; should have been catched before, even if occuring.
       if (i_R_0 + N - 1 > sbt_N_spl) then
          call aims_stop('Separation too large (index)')
       end if
       
       do L = L_min, L_max
          do i = 1, N
             i_spl = i_R_0 + (i - 1)*sbt_spl_multiplier
             KK(i, L) = sbt_KK_spl(1, i_spl, L) &
             &        + sbt_KK_spl(2, i_spl, L) * d_R_0 &
             &        + sbt_KK_spl(3, i_spl, L) * d_R_0_sq &
             &        + sbt_KK_spl(4, i_spl, L) * d_R_0_cb
          end do
       end do
    else
       allocate(wsave(2*N+15), stat=info)
       call check_allocation(info, 'wsave', func)
       call dffti(N, wsave)
       do L = L_min, L_max
          ! Set up t-grid approximation of kernel j_L(t) [possibly hybrid]
          call logsbt_kernel(N, KK(:,L), lnt0, lnrange, L, power, wsave,.true.)
       end do
       deallocate(wsave)
    end if

  end subroutine get_fast_kernel
   !----------------------------------------------------------------------------
  !****s* logsbt_fast_kernel/get_fast_kernel_der
  !  NAME
  !    get_fast_kernel_der
  !  SYNOPSIS

  subroutine get_fast_kernel_der(N, L_min, L_max, lnt0, lnrange, power, &
  &                          use_spline, KK)

    !  PURPOSE
    !
    !    Evaluate derivative of kernel from spline prepared by intialize_fast_kernel().
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N, L_min, L_max
    real*8, intent(IN) :: lnt0, lnrange
    real*8, intent(IN) :: power
    logical, intent(IN) :: use_spline
    real*8, intent(OUT) :: KK(N, L_min:L_max)

    !  INPUTS
    !    o N -- Number of grid points [size(KK, 1)]
    !    o L -- Angular momentum
    !    o lnt0, lnrange -- Range of kernel to retrieve
    !                       Must be /within/ the corresponding arguments
    !                       to initialize_fast_kernel().
    !    o power -- Power bias [0.d0]
    !    o use_spline -- If .false., completely reconstruct kernel
    !  OUTPUTS
    !    o KK -- Derivative of kernel with respect to lnt0
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: L, i_R_0, i, i_spl
    real*8 :: Rho
    real*8 :: d_R_0, d_R_0_sq, d_R_0_cb
    real*8, allocatable :: wsave(:)
    integer :: info
    character(*), parameter :: func = 'get_fast_kernel_der'

    if (use_spline) then
       if (L_max > sbt_L_spl) call aims_stop('L_max too large', func)
       Rho = lnt0 - sbt_lnt0_spl
       if (Rho < 0.d0 .or. Rho+lnrange > sbt_lntrange_spl) then
          call aims_stop('Invalid |R|', func)
       end if

       ! (lnrage / N) is more coarse (larger step length) than
       ! (sbt_lntrange_spl / sbt_N_spl).
       if (abs(sbt_spl_multiplier * (sbt_lntrange_spl / sbt_N_spl) &
       &       - (lnrange / N)) > 1d-10) then
          call aims_stop('N is not multiple of sbt_N_spl', func)
       end if
       if (abs(power - sbt_power_spl) > 1d-10) then
          call aims_stop('Initialized with different power', func)
       end if

       ! --- Spline parameters

       ! invert_log_grid(): i(r) = 1 + [ln(r / r_grid_min)]/[ln(r_grid_inc)]
       d_R_0 = 1.d0 + Rho / (sbt_lntrange_spl / sbt_N_spl)
       i_R_0 = floor(d_R_0)
       d_R_0 = d_R_0 - i_R_0
       d_R_0_sq = d_R_0 * d_R_0
       d_R_0_cb = d_R_0_sq * d_R_0
       ! Security check; should have been catched before, even if occuring.
       if (i_R_0 + N - 1 > sbt_N_spl) then
          call aims_stop('Separation too large (index)')
       end if
       
       do L = L_min, L_max
          do i = 1, N
             i_spl = i_R_0 + (i - 1)*sbt_spl_multiplier
             KK(i, L) = sbt_KK_spl(2, i_spl, L) &
             &     + 2* sbt_KK_spl(3, i_spl, L) * d_R_0 &
             &     + 3* sbt_KK_spl(4, i_spl, L) * d_R_0_sq
             KK(i, L) = KK(i, L) / (sbt_lntrange_spl / sbt_N_spl)
          end do
       end do
    else
       call aims_stop('get_fast_kernel_der needs use_spline to be set')
    end if

  end subroutine get_fast_kernel_der
  
  
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt_fast_kernel/sbt_splined_kernel
  !  NAME
  !    sbt_splined_kernel
  !  SYNOPSIS

  subroutine sbt_splined_kernel(N, Nt, lnt0, lntrange, max_L, power, KK_spl)

    !  PURPOSE
    !
    !    Get a spline approximation to a set of kernels.  Please note the
    !    difference between N (the number of spline nodes) and Nt (the
    !    number of Fourier terms in the kernel.  N >= Nt should always hold!
    !
    !    In case of power==0.d0, the L==0 kernel is hybrid.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    use spline, only: cubic_spline
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    integer, intent(IN) :: Nt
    real*8, intent(IN) :: lnt0, lntrange
    integer, intent(IN) :: max_L
    real*8, intent(IN) :: power
    real*8, intent(OUT) :: KK_spl(4, N, 0:max_L)

    !  INPUTS
    !    o N -- Number of spline nodes
    !    o Nt -- Number of Fourier terms (>= N)
    !    o lnt0, lntrange -- Range of loggrid: [exp(lnt0), exp(lnt0+lntrange)[
    !    o max_L -- Kernels up to this
    !    o power -- Power bias [in: sqrt(2/pi) * j_L(t) * t^power]
    !  OUTPUTS
    !    o KK_spl -- Output splines
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: erfc_width, tau0
    real*8, allocatable :: wsave_N(:)
    real*8, allocatable :: KKbar(:), KKtmp(:)
    integer :: info, L
    logical :: do_hybrid
    character*50 :: filename
    character(*), parameter :: func = 'sbt_splined_kernel'

    if (N < Nt) call aims_stop('N < Nt', func)
    allocate(wsave_N(2*N+15), stat=info)
    call check_allocation(info, 'wsave', func)
    allocate(KKbar(N), KKtmp(N), stat=info)
    call dffti(N, wsave_N)
    do_hybrid = power <= 1d-10
    erfc_width = 15 * lntrange/Nt
    tau0 = logsbt_lntonset - lnt0

    ! --- L > 0

    do L = 0, max_L
       ! Set up Nt-term kr-grid approximation of kernel j_L(kr)
       KKbar = 0.d0
       call logsbt_vanilla_kernel(Nt, KKbar, lnt0, lntrange, L, power)
       if (L == 0 .and. do_hybrid) then
          call logsbt_erfc_logFT(Nt, KKtmp, lnt0, lntrange, erfc_width, tau0)
          KKbar(1:Nt) = KKbar(1:Nt) - KKtmp(1:Nt)
       end if
       call dfftb(N, KKbar, wsave_N)    ! Throw away N+1:Nt
       KKbar = KKbar / lntrange
       if (L == 0 .and. do_hybrid) then
          ! Real-space correction only for 1:N resolution.
          call logsbt_erfc_prodspace(N, KKtmp, lnt0, lntrange, erfc_width, tau0)
          KKbar = KKbar + KKtmp * dble(N) / lntrange
       end if

       ! Spline it
       call cubic_spline(KKbar, N, KK_spl(:,:, L))
    end do
    deallocate(wsave_N, KKbar)

  end subroutine sbt_splined_kernel
  !******
end module logsbt_fast_kernel
!******



