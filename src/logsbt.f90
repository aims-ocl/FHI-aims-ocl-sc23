!****h* FHI-aims/logsbt
!  NAME
!    logsbt
!  SYNOPSIS

module logsbt

  !  PURPOSE
  !
  !    Provide subroutines to perform numerical spherical Bessel transorms
  !    (SBT) on a logarithmic grid using the algorithms of James D. Talman
  !    and making use of some concepts of Andrew J.S. Hamilton.
  !
  !    The basic idea is to calculate the SBT
  !
  !        f_l~(k) = sqrt(2/pi) \int_0^\infty r^2 j_l(kr) f_l(r)     (*)
  !
  !    of the radial component f_l(r).  It is the radial part of the Fourier
  !    transform of f(\bm r):
  !
  !        f(\bm r) = f_l(r) Y_{lm}(\hat r)
  !        f~(\bm k) = (2pi)^(-3/2) \int d^3r e^{-i\bm k \bm r} f(\bm r)
  !                  = i^-l f_l~(k) Y_{lm}(\hat k)
  !
  !    The actual integral (*) is calculated in logspace:
  !
  !        f_l~(k) = k^-power \int_0^\infty dr/r              x ! log grid
  !                            x sqrt(2pi) j_l(kr) (kr)^power x ! kernel
  !                            x r^(3-power) f_l(r)             ! array
  !
  !  THE POWER PARAMETER
  !
  !    The 'power' paramter is very important for the numerical stability of
  !    the procedure.  If f_l(r) is given with the accuracy epsilon*f_l(r) [to
  !    full relative accuracy], scaling f_l(r) does no harm.  But
  !    f_l~(k)*k**power is (because of intermediate FFTs) accurate only to
  !    machine precision in absolute terms, even if the true value is much
  !    smaller than 1e-15.  Therefore, power should be chosen according to the
  !    usage of f_l~(k).  It is assumed to be used in an expression like
  !    f_l~(k) * k**power.
  !
  !    In order to emphasize the importance of the power parameter, the
  !    scaling (multiplication with r**power or k**power) has to be done
  !    explicitly by the user of this module by the corresponding calls to
  !    logsbt_scale().  Bear in mind that every single logsbt_scale(...,
  !    power) call is in danger of magnifying numerical noise by huge factors
  !    (like exp(-lnr0*power) ~ 3e16**power).  It should therefore only be
  !    done for quantities which are given to relative numerical accuracy
  !    (which is then kept).
  !
  !    Alternatively, use logsbt_multi_double_driver, which goes along the
  !    lines of Talman [2] in performing the transform twice using different
  !    powers.  The resulting f~(k) should be accurate from f~(k) to
  !    f~(k)*k**1.5.
  !
  !  USER VISIBLE CONVENTIONS
  !
  !    Ranges in r- and k-space are given by the logarithmic onsets, lnr0 and
  !    lnk0 (lnr1 and lnk1 would have been better names), the logarithmic
  !    range lnrange, and the number of logarithmic grid points N.  The grid
  !    spacing is therefore dln=lnrange/N.  The grid points are then given by
  !    r_i = exp(lnr0 + (i-1)*dlnr), and k_i = exp(lnk0 + (i-1)*dlnr),
  !    i=1,...,N.  Please note that for algorithmic reasons, the grid spacing
  !    and the number of points must be equal both in r- and k-space.  For
  !    convenience, rho = log(r/r0), kappa = log(k/k0), and tau = log(t/t0) =
  !    log(kr/k0r0).
  !
  !    WARNING: If you do not know exactly what you are doing, you should
  !    restrict usage to logsbt_multi_double_driver, logsbt_scale, and maybe
  !    logsbt_r2ilog and logsbt_ilog2r.
  !
  !  INTERNAL CONVENTIONS
  !
  !    There is an intermediate representation of functions called logFT
  !    space, here signified by a "bar" in the variable name.  It contains
  !    Fourier coefficients as used by the dfftpack routines.  See
  !    logsbt_bar2cmplx() and logsbt_cmplx2bar() for details.
  !
  !    Without resorting to logFT space, the Kernel would be needed from
  !    lnk0+lnr0 (which is negative) to lnk0+lnr0+2*lnrange.  The algorithm
  !    assumes log-periodicity in the scaled input (r^(3-power) f(r)), and
  !    therefore the kernel is also lnrange-periodic.  For some applications,
  !    the corresponding ringing would be too severe.  In these cases (ask
  !    logsbt_opt_trange), lntrange -- the periodicity of the kernel -- is set
  !    to 2*lnrange.  But a kernel can be constructed for any onset 'lnt0',
  !    periodicty 'lntrange', and power bias 'power'.  The case of L=power=0
  !    needs special treatment by logsbt_erfc_correction, though, because
  !    j_0(kr) does not have a continuous logFT because it tends to one for
  !    ln(kr)->-\infty.
  !
  !  USES

  implicit none

  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    [1] Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !        Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !        "Ab initio simulations with Numeric Atom-Centered Orbitals:
  !        FHI-aims", Computer Physics Communications 180, 2175 (2009).
  !
  !    [2] J.D. Talman, "NumSBT: A subroutine for calculating spherical Bessel
  !        transforms numerically," Computer Physics Communications 180, 332
  !        (2009).
  !
  !    [3] A. J. S. Hamilton, "Uncorrelated modes of the non-linear power
  !        spectrum" Monthly Notices of the Royal Astronomical Society 312, 257
  !        (2000).
  !
  !  COPYRIGHT
  !    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !    e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8, parameter :: logsbt_lntonset = 0.d0 ! center of erfc correction

  integer, parameter :: i_logft_aliasing = 1
  integer, parameter :: i_sbt_aliasing_smallk = 2
  integer, parameter :: i_sbt_aliasing_largek = 3
  integer, parameter :: i_sbt_ringing_smallr = 4
  integer, parameter :: i_sbt_ringing_larger = 5
  integer, parameter :: n_sbt_errors = 5

contains

  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_multi_double_driver
  !  NAME
  !    logsbt_multi_double_driver
  !  SYNOPSIS

  subroutine logsbt_multi_double_driver(N, lnr0, lnk0, lnrange, &
  &                                     n_fn, fn_to_L, ff)

    !  PURPOSE
    !
    !    Perform the logSBT twice [1] on all radial parts to get a numerically
    !    stable transform.  The resulting f~(k) may be scaled up to
    !    f~(k)*k**1.5 without loosing (absolute) accuracy.
    !
    !  USES

    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: n_fn
    integer, intent(IN) :: fn_to_L(n_fn)
    real*8, intent(INOUT) :: ff(N, n_fn)

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- sbtgrid parameters (see module header)
    !    o n_fn -- number of radial parts
    !    o fn_to_L -- angular mumenta of radial parts
    !    o ff -- f(r) (no scaling)
    !  OUTPUTS
    !    o ff -- f~(k) (no scaling)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8, allocatable :: ffk_high(:,:), ffk_low(:,:)
    real*8 :: mindiff, this_diff, absmax, ffk1
    integer :: leftpos, rightpos, minpos
    integer :: i, i_fn
    integer :: info
    character(*), parameter :: func = 'logsbt_multi_double_driver'

    ! --- Initialization

    minpos = 1

    allocate(ffk_high(N, n_fn), ffk_low(N, n_fn), stat=info)
    call check_allocation(info, 'ffk_high, ffk_low', func)

    ! --- Transforms
    ffk_high = ff
    call logsbt_scale(N, n_fn, ffk_high, lnr0, lnrange, 1.5d0)
    ! JW: In principle, the doubling of the domain should happen on the left
    !     hand side [2] for ffk_high.
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, 1.5d0, &
    &                        n_fn, fn_to_L, ffk_high)
    call logsbt_scale(N, n_fn, ffk_high, lnk0, lnrange, -1.5d0)
    ! ffk_high == f~(k)      [error: epsilon / k**1.5]

    ffk_low = ff
    call logsbt_scale(N, n_fn, ffk_low, lnr0, lnrange, 3d0)
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, 0d0, &
    &                        n_fn, fn_to_L, ffk_low)
    ! ffk_low == f~(k)       [error: epsilon]

    ! --- Join

    do i_fn = 1, n_fn
       ! Do not search in tails
       absmax = maxval(abs(ffk_low))
       ffk1 = ffk_low(1, i_fn)
       do leftpos = 1, N
          if (abs(ffk_low(leftpos, i_fn) - ffk1) > 1d-5 * absmax) exit
       end do
       do rightpos = N, 1, -1
          if (abs(ffk_high(rightpos, i_fn)) > 1d-5 * absmax) exit
       end do
       if (leftpos > rightpos) leftpos = 1
       ! Find minimal distance
       mindiff = huge(mindiff)
       do i = leftpos, rightpos
          this_diff = abs(ffk_high(i, i_fn) - ffk_low(i, i_fn))
          if (this_diff < mindiff) then
             mindiff = this_diff
             minpos = i
          end if

       end do
       ! Join
       ff(:minpos, i_fn) = ffk_low(:minpos, i_fn)
       ff(minpos:, i_fn) = ffk_high(minpos:, i_fn)
    end do

  end subroutine logsbt_multi_double_driver
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_multi_driver
  !  NAME
  !    logsbt_multi_driver
  !  SYNOPSIS

  subroutine logsbt_multi_driver(N, lnr0, lnk0, lnrange, power_bias, &
  &                              n_fn, fn_to_L, ff, sbt_errors)

    !  PURPOSE
    !
    !    Perform SBT on n_fn radial parts in one sweep.  Input and output ff
    !    needs to be scaled.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: power_bias
    integer, intent(IN) :: n_fn, fn_to_L(n_fn)
    real*8, intent(INOUT) :: ff(N, n_fn)
    real*8, intent(OUT), optional :: sbt_errors(n_sbt_errors, n_fn)

    !  INPUTS
    !    o N -- Number of SBT grid points
    !    o lnr0, lnk0, lnrange -- SBT grid extents
    !    o power_bias -- SBT power bias (alpha)
    !    o n_fn -- Number of distinct radial parts to transform
    !    o fn_to_L -- Angular momentum of fn
    !    o ff -- f(rr) * rr**(3-power_bias)
    !  OUTPUTS
    !    o ff -- f~(kk) * rr**power_bias
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: max_L, L, i_fn
    real*8 :: lnt0
    integer, allocatable :: Nt(:)
    real*8, allocatable :: wsave(:,:), KKbar(:,:), lntrange(:)
    real*8 :: these_errors(n_sbt_errors)
    integer :: info
    character*150 :: info_str, filename
    character(*), parameter :: func = 'logsbt_multi_driver'

    max_L = maxval(fn_to_L)

    ! --- Get SBT kernels

    allocate(Nt(0:max_L), lntrange(0:max_L), stat=info)
    call check_allocation(info, 't-range', func)
    allocate(wsave(4*N+15, 0:max_L))
    call check_allocation(info, 'wsave', func)
    allocate(KKbar(2*N, 0:max_L), stat=info)
    call check_allocation(info, 'KKbar', func)

    lnt0 = lnk0 + lnr0
    do L = 0, max_L
       call logsbt_opt_trange(N, L, power_bias, lnrange, Nt(L), lntrange(L))
       call dffti(Nt(L), wsave(:,L))
       call logsbt_kernel(Nt(L), KKbar(:,L), lnt0, lntrange(L), L, power_bias,&
       &                  wsave(:,L))
    end do

    ! --- Do SBT

    ! ff == f(r)                     ! [assume full relative accuracy]
    ! call logsbt_scale(N, n_fn, ff, lnr0, lnrange, 3.d0 - power_bias)
    ! ff == f(r)*k**(3-power_bias)   ! [assume full relative accuracy]

    do i_fn = 1, n_fn
       L = fn_to_L(i_fn)
       if (L > max_L .or. L < 0) call aims_stop('Invalid L', func)
       call logsbt_expert(N, ff(:, i_fn), &
       &                  lnr0, lnk0, lnrange, L, power_bias, &
       &                  Nt(L), KKbar(:, L), lntrange(L), wsave(:, L), &
       &                  these_errors)
       if (present(sbt_errors)) sbt_errors(:, i_fn) = these_errors
    end do

    ! ff == f(k)*k**power_bias       ! [error: epsilon]
    ! call logsbt_scale(N, n_fn, ff, lnk0, lnrange, power_add)
    ! ff == f(k)*k**(p_bias+p_add)   ! [error: epsilon * k**power_add]

    ! --- Tidy up

    deallocate(Nt, wsave, KKbar)

  end subroutine logsbt_multi_driver
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_scale
  !  NAME
  !    logsbt_scale
  !  SYNOPSIS

  subroutine logsbt_scale(N, n_fn, ff, ln0, lnrange, power)

    !  PURPOSE
    !    Scale the function ff given on a logarithmic grid by rr^power.
    !  USES

    implicit none

    !  ARGUMENTS

    !f2py real*8, intent(IN,OUT,COPY) :: ff(N)
    integer, intent(IN) :: N            ! number of grid points
    integer, intent(IN) :: n_fn         ! number of functions
    real*8, intent(INOUT) :: ff(N, n_fn)! array to transform
    real*8, intent(IN) :: ln0           ! log of left border of log grid
    real*8, intent(IN) :: lnrange       ! range of log grid (dln = lnrange/N)
    real*8, intent(IN) :: power         ! power to apply

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: inc, fac, dln, ln
    integer :: i_ln, i_bunch, i_bunch_max
    integer, parameter :: n_in_bunch = 256

    if (power == 0.d0) return   ! Nothing to do

    dln = lnrange / N
    inc = exp(power * dln) ! = dr^q_eff
    i_bunch_max = (n / n_in_bunch) * n_in_bunch

    ! WPH, 11 May 2018: I've performed peeling on the main loop to circumvent
    ! a PGI bug when vectorization is used for compiling this subroutine, which
    ! caused every EXX calculation to crash.  The original code is left here
    ! but is commented out.

    !do i_bunch = 0, n-1, n_in_bunch
    do i_bunch = 0, i_bunch_max-1, n_in_bunch
       ! recalculate exact exponential to avoid noise accumulation
       ln = ln0 + i_bunch * dln
       fac = exp(power * ln) ! = r^q_eff
       !do i_ln = i_bunch+1, min(n, i_bunch+n_in_bunch)
       do i_ln = i_bunch+1, i_bunch+n_in_bunch
          ! ff(i_ln) = ff(i_ln) * exp(power * ln) ! = ff(i) * r**power
          ff(i_ln, :) = ff(i_ln, :) * fac ! = ff(i) * r**power
          fac = fac * inc
       end do
    end do

    ! WPH: This is the peeled-off loop.
    i_bunch = i_bunch + n_in_bunch
    ln = ln0 + i_bunch * dln
    fac = exp(power * ln) ! = r^q_eff
    do i_ln = i_bunch+1, n
       ff(i_ln, :) = ff(i_ln, :) * fac ! = ff(i) * r**power
       fac = fac * inc
    end do

  end subroutine logsbt_scale
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_r2ilog
  !  NAME
  !    logsbt_r2ilog
  !  SYNOPSIS

  real*8 function logsbt_r2ilog(r, ln1, dln)

    !  PURPOSE
    !
    !    Calculate log-index coordinate from real-space radius.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: r, ln1, dln

    !  INPUTS
    !    o r -- radius (coordinate)
    !    o ln1 -- log of onset
    !    o dln -- log-distence between points
    !  OUTPUTS
    !    logsbt_r2ilog -- index coordinate (log(r) == ln1 -> 1.,
    !                                       log(r) == ln1+dln -> 2.)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'logsbt_r2ilog'

    logsbt_r2ilog = (log(r) - ln1) / dln + 1.d0

  end function logsbt_r2ilog
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_ilog2r
  !  NAME
  !    logsbt_ilog2r
  !  SYNOPSIS

  real*8 function logsbt_ilog2r(ilog, ln1, dln)

    !  PURPOSE
    !
    !    Calculate log-index coordinate from real-space radius.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: ilog, ln1, dln

    !  INPUTS
    !    o ilog -- index coordinate
    !    o ln1 -- log of onset
    !    o dln -- log-distence between points
    !  OUTPUTS
    !    logsbt_ilog2r -- radius (ilog == 1. -> exp(ln1),
    !                             ilog == 2. -> exp(ln1+dln))
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'logsbt_ilog2r'

    logsbt_ilog2r = exp(ln1 + (ilog-1.d0)*dln)

  end function logsbt_ilog2r
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_driver
  !  NAME
  !    logsbt_driver
  !  SYNOPSIS

  subroutine logsbt_driver(N, ff, lnr0, lnk0, lnrange, L, power, sbt_errors)

    !  PURPOSE
    !
    !    Wrapper around logsbt_expert() to set up the kernel and the dfftpack
    !    workspace.  So in general, easier to use but less efficient.  Input
    !    and output ff need to be scaled.
    !
    !    Expects
    !      ff == f(r) * r**(3.d0 - power)
    !    as input and outputs
    !      ff == f~(k) * k**power
    !          
    !
    !  USES

    implicit none

    !  ARGUMENTS

    !f2py real*8, intent(IN,OUT,COPY) :: ff(N)
    !f2py integer, parameter :: n_sbt_errors = 5
    integer, intent(IN) :: N
    real*8, intent(INOUT) :: ff(N)
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: L
    real*8, intent(IN) :: power
    real*8, intent(OUT), optional :: sbt_errors(n_sbt_errors)

    !  INPUTS
    !    o N -- number of grid points
    !    o lnr0 -- log of onset of r-space grid
    !    o lnk0 -- log of onset of k-space grid
    !    o lnrange -- (log) range of both r- and k-space grid
    !    o L -- angular momentum
    !    o power -- power bias of the transform
    !    o ff -- f(rr) * rr**(3-power)
    !  OUTPUTS
    !    o ff -- f~(kk) * rr**power
    !  OUTPUTS
    !    o sbt_errors -- error indicators (optional)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable :: KKbar(:), wsave(:)
    integer :: Nt
    real*8 :: lnt0, lntrange

    call logsbt_opt_trange(N, L, power, lnrange, Nt, lntrange)
    lnt0 = lnk0 + lnr0

    allocate(KKbar(Nt), wsave(2*Nt+15))
    call dffti(Nt, wsave)
    call logsbt_kernel(Nt, KKbar, lnt0, lntrange, L, power, wsave)
    call logsbt_expert(N, ff, lnr0, lnk0, lnrange, L, power, &
    &                  Nt, KKbar, lntrange, wsave, sbt_errors)

    deallocate(KKbar, wsave)

  end subroutine logsbt_driver
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_single_k
  !  NAME
  !    logsbt_single_k
  !  SYNOPSIS

  subroutine logsbt_single_k(N, ff, lnr0, lnrange, L, power, k0, II)

    !  PURPOSE
    !
    !    Calculate ffk(k) for k=k0 only.  For this, a taylored kernel is
    !    constructed, transformed to kr-space and integrated using the
    !    trapezoidal rule.  Scalings are done within this routine as an
    !    exception to the general rule.
    !
    !  USES
    use constants, only: pi
    use mpi_tasks, only: aims_stop
    implicit none
    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(INOUT) :: ff(N)
    real*8, intent(IN) :: lnr0, lnrange
    integer, intent(IN) :: L
    real*8, intent(IN) :: power
    real*8, intent(IN) :: k0
    real*8, intent(OUT) :: II

    !  INPUTS
    !    o N -- number of grid points
    !    o ff -- function to be transformed (lost/scaled on output)
    !    o lnr0 -- log of onset of r-space grid
    !    o lnrange -- (log) range of both r- and k-space grid
    !    o L -- angular momentum
    !    o power -- power bias of the transform
    !  OUTPUTS
    !    o II -- ff~(k0)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: Nt
    real*8 :: lnt0, lntrange
    real*8, allocatable :: KK(:), wsave(:)

    if (k0 < 0.d0) then
       call aims_stop('Negative k0')
    else if (k0 < 1d-12) then
       ! power = 0.d0
       call logsbt_scale(N, 1, ff, lnr0, lnrange, 3.d0)
       if (L /= 0) then
          II = 0.d0
       else
          II = sqrt(2.d0 / pi) * sum(ff) * lnrange / N
       end if
    else
       ! Prepare ff
       call logsbt_scale(N, 1, ff, lnr0, lnrange, 3.d0 - power)

       Nt = N
       lntrange = lnrange
       lnt0 = lnr0 + log(k0)

       ! Get kernel [i.e. j_L(kr) (kr)^power] in product space (t-space)
       allocate(KK(Nt), wsave(2*Nt+15))
       call dffti(Nt, wsave)
       call logsbt_kernel(Nt, KK, lnt0, lntrange, L, power, wsave, .true.)

       ! Trapozoidal rule (without borders)
       II = sum(ff * KK) * lnrange

       II = II / k0**power / N

       ! Tidy up
       deallocate(KK, wsave)
    end if


  end subroutine logsbt_single_k
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_expert
  !  NAME
  !    logsbt_expert
  !  SYNOPSIS

  subroutine logsbt_expert(N, ff, lnr0, lnk0, lnrange, L, power, &
  &                        Nt, KKbar, lntrange, wsave, sbt_errors)

    !  PURPOSE
    !    Apply spherical Bessel transform to loggrid ff.  Inputs and outputs
    !    need to be scaled.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    !f2py real*8, intent(IN,OUT,COPY) :: ff(N)
    !f2py integer, parameter :: n_sbt_errors = 5
    integer, intent(IN) :: N
    real*8, intent(INOUT) :: ff(N)
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: L
    real*8, intent(IN) :: power
    integer, intent(IN) :: Nt
    real*8, intent(IN) :: KKbar(Nt)
    real*8, intent(IN) :: lntrange
    real*8, intent(INOUT) :: wsave(2*Nt+15)
    real*8, intent(OUT), optional :: sbt_errors(n_sbt_errors)

    !  INPUTS
    !    o N -- number of grid points
    !    o lnr0 -- log of onset of r-space grid
    !    o lnk0 -- log of onset of k-space grid
    !    o lnrange -- (log) range of both r- and k-space grid
    !    o L -- angular momentum
    !    o power -- power bias of the transform (*must* match kernel)
    !    o Nt -- number of grid points of kernel
    !    o KKbar -- kernel of transform, either vanilla (N) or hybrid (2*N)
    !    o lntrange -- range of kernel in log-space
    !    o wsave -- dfftpack worksapce, either (2*Nt+15)
    !    o ff -- f(rr) * rr**(3-power_bias)
    !  OUTPUTS
    !    o ff -- f~(kk) * rr**power_bias
    !    o sbt_errors -- error indicators
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: logft_aliasing
    real*8 :: A(2), Lr(2)
    real*8, allocatable :: ff2(:)
    integer :: i
    character(*), parameter :: func = 'logsbt_expert'

    if (abs(lnrange / N - lntrange / Nt) > 1d-10) then
       call aims_stop('grid mismatch of kernel')
    end if

    if (present(sbt_errors)) then
       sbt_errors(i_sbt_ringing_smallr) = maxval(abs(ff(1:4)))
       sbt_errors(i_sbt_ringing_larger) = maxval(abs(ff(N-3:N)))
    end if

    if (Nt == N) then
       call logsbt_apply_vanilla_kernel(N, ff, KKbar, wsave, logft_aliasing)
    else if (Nt > N) then
       allocate(ff2(Nt))
       ff2(1:N) = ff
       ff2(N+1:) = 0.d0
       call logsbt_apply_vanilla_kernel(Nt, ff2, KKbar, wsave, logft_aliasing)
       ff = ff2(1:N)
       deallocate(ff2)
    else
       call aims_stop('Invalid Nt', func)
    end if

    if (present(sbt_errors)) then
       call logsbt_analyze_tails(N, ff, lnk0, lnrange, A, Lr)
       sbt_errors(i_sbt_aliasing_smallk) = abs(A(1)) * exp(lnk0 * Lr(1))
       sbt_errors(i_sbt_aliasing_largek) = abs(A(2)) * &
       &                                   exp((lnk0+lnrange) * Lr(2))
       sbt_errors(i_logft_aliasing) = logft_aliasing
       if (Nt >= 2*N) then
          ! Left edge of f~(k) is taken care of by enhanced Nt
          sbt_errors(i_sbt_aliasing_smallk) = &
          & sbt_errors(i_sbt_aliasing_smallk) * exp(lnk0 * 2)
       end if
    end if

  end subroutine logsbt_expert
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_apply_vanilla_kernel
  !  NAME
  !    logsbt_apply_vanilla_kernel
  !  SYNOPSIS

  subroutine logsbt_apply_vanilla_kernel(N, ff, KKbar, wsaveN, logft_aliasing)

    !  PURPOSE
    !  
    !    Apply kernel KKbar to ff, which has to be transformed to logFT
    !    inbetween.  This (and the wrappred routine
    !    logsbt_apply_kernel_on_ffbar) is the actual workhorse of this module.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    !f2py real*8, intent(IN,OUT,COPY) :: ff(N)
    integer, intent(IN)   :: N              ! size of loggrid
    real*8, intent(INOUT) :: ff(N)          ! function values
    real*8, intent(IN)    :: KKbar(N)       ! (logFT) kernel to apply
    real*8, intent(INOUT) :: wsaveN(2*N+15) ! dfftpack workspace
    real*8, intent(OUT)   :: logft_aliasing ! Estimate of logFT-aliasing

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    call dfftf(N, ff, wsaveN)

    logft_aliasing = maxval(ff(N-3:N))

    call logsbt_apply_kernel_on_ffbar(N, ff, KKbar)

    logft_aliasing = max(logft_aliasing, maxval(ff(N-3:N)))
    logft_aliasing = logft_aliasing / sqrt(dble(N))

    call dfftb(N, ff, wsaveN)
    ff = ff / N

  end subroutine logsbt_apply_vanilla_kernel
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_apply_kernel_on_ffbar
  !  NAME
  !    logsbt_apply_kernel_on_ffbar
  !  SYNOPSIS

  subroutine logsbt_apply_kernel_on_ffbar(N, ffbar, KKbar)

    !  PURPOSE
    !    Apply kernel KKbar to ffbar, both given in logFT space.
    !    That is, return conj(ffbar) * KKbar
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN)   :: N           ! size of loggrid
    real*8, intent(INOUT) :: ffbar(N)    ! logFT of ff / fftil
    real*8, intent(IN)    :: KKbar(N)    ! (logFT) kernel to apply

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i
    real*8 :: fr, fi, kr, ki

    fr = ffbar(1)
    kr = KKbar(1)
    ffbar(1) = fr * kr
    do i = 1, N/2-1
       fr = ffbar(2*i)
       fi = - ffbar(2*i+1)   ! conj(ffbar)
       kr = KKbar(2*i)
       ki = KKbar(2*i+1)
       ffbar(2*i) = fr * kr - fi * ki
       ffbar(2*i+1)   = fr * ki + fi * kr
    end do
    if (modulo(N, 2) == 0) then
       fr = ffbar(N)
       kr = ffbar(N)
       ffbar(N) = fr * kr
    end if

  end subroutine logsbt_apply_kernel_on_ffbar
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_opt_trange
  !  NAME
  !    logsbt_opt_trange
  !  SYNOPSIS

  subroutine logsbt_opt_trange(N, L, power, lnrange, Nt, lntrange, need_erfc)

    !  PURPOSE
    !
    !    Return the optimal kernel range parameters Nt, lntrange.  See the
    !    module header INTERNAL CONVENTIONS for details.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N    ! number of grid points
    integer, intent(IN) :: L    ! angular momentum
    real*8, intent(IN) :: power ! power bias
    real*8, intent(IN) :: lnrange ! function range
    integer, intent(OUT), optional :: Nt  ! number of kernel grid points
    real*8, intent(OUT), optional  :: lntrange ! kernel range
    logical, intent(OUT), optional :: need_erfc ! need erfc correction?

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: asymp_left, asymp_right

    asymp_left = L + power
    asymp_right = 4 - power
    if (asymp_left < 3d0 .or. asymp_right < 3d0) then
       ! Suboptimal asymptotic behavior of a~(k) = f~(k) / k^power.
       ! To avoid SBT-aliasing, use doubled integration range
       ! Depending on the defaults of lnr0, lnk0 and lnrange, looser criteria
       ! could be used if the actual SBT should impose a performance issue.
       if (present(lntrange)) lntrange = 2*lnrange
       if (present(Nt)) Nt = 2*N
    else
       ! SBT decays fast enough
       if (present(lntrange)) lntrange = lnrange
       if (present(Nt)) Nt = N
    end if
    if (present(need_erfc)) need_erfc = (abs(power) < 1d-10 .and. L == 0)

  end subroutine logsbt_opt_trange
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_kernel
  !  NAME
  !    logsbt_kernel
  !  SYNOPSIS

  subroutine logsbt_kernel(Nt, KKbar, lnt0, lntrange, L, power, &
  &                        wsave, prodspace)

    !  PURPOSE
    !
    !    Wrapper around logsbt_vanilla_kernel() and logsbt_erfc_correction()
    !    to generate the kernel KKbar.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    !f2py integer, intent(IN), required :: Nt
    integer, intent(IN) :: Nt        ! grid size
    real*8, intent(OUT) :: KKbar(Nt) ! kernel (see above)
    real*8, intent(IN) :: lnt0       ! grid onset (= lnk0 + lnr0)
    real*8, intent(IN) :: lntrange   ! grid range
    integer, intent(IN) :: L         ! angular momentum
    real*8, intent(IN) :: power      ! incorporated scaling term
    real*8, intent(INOUT), optional :: wsave(2*Nt+15) ! dfft working space
    logical, intent(IN), optional :: prodspace ! default: .false.

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    logical :: need_erfc, need_wsave, my_prodspace
    real*8, allocatable :: my_wsave(:)

    ! Prepare
    need_wsave = .not. present(wsave)
    if (need_wsave) then
       allocate(my_wsave(2*Nt+15))
       call dffti(Nt, my_wsave)
    end if
    my_prodspace = .false.; if (present(prodspace)) my_prodspace = prodspace

    ! Get kernel in logFT
    call logsbt_vanilla_kernel(Nt, KKbar, lnt0, lntrange, L, power)
    need_erfc = (abs(power) < 1d-10 .and. L == 0)
    if (need_erfc) then
       if (need_wsave) then
          call logsbt_erfc_correction(Nt, KKbar, lnt0, lntrange, my_wsave)
       else
          call logsbt_erfc_correction(Nt, KKbar, lnt0, lntrange, wsave)
       end if
    end if

    ! Transfer back
    if (my_prodspace) then
       if (need_wsave) then
          call dfftb(Nt, KKbar, my_wsave)
       else
          call dfftb(Nt, KKbar, wsave)
       end if
       KKbar = KKbar / lntrange
    end if
    if (need_wsave) deallocate(my_wsave)

  end subroutine logsbt_kernel
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_erfc_correction
  !  NAME
  !    logsbt_erfc_correction
  !  SYNOPSIS

  subroutine logsbt_erfc_correction(Nt, KKbar, lnt0, lntrange, wsave)

    !  PURPOSE
    !
    !    In the general case, the kernel is purely set up in logFT space as a
    !    Fourier sum, which leads to a lnrange-periodic sequence in kr-product
    !    space.
    !
    !    This works well for decaying kernels.  But j_0 tends to 1. for small
    !    arguments, and if there is no power bias, the Fourier transform does
    !    not exist.
    !
    !    To get along this problem, a hybrid approach is used.  A smooth
    !    function containing the asymptotic part, erfc(tau-lntonset), is
    !    subtracted from the kernel.  The part corresonpding to this function
    !    is than evaluated using the trapezoidal rule.  This can also be done
    !    in logFT space, so that both terms are put together in one kernel.
    !
    !    The major drawback is that, first, the grid has to be doubled for the
    !    trapezoidal rule to be valid in the whole range, and second, the
    !    smooth function cannot be continuously Fourier transformed and the
    !    kernel is therefore only valid for this particular discretization.
    !    That is, the result of applying the kernel cannot be Fourier
    !    interpolated!
    !
    !    Could be optimized performancewise...
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: Nt               ! grid size
    real*8, intent(INOUT) :: KKbar(Nt)      ! kernel (see above)
    real*8, intent(IN) :: lnt0              ! grid onset (= lnk0 + lnr0)
    real*8, intent(IN) :: lntrange          ! grid range
    real*8, intent(INOUT) :: wsave(2*Nt+15) ! dfftpack workspace
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE
    real*8, allocatable :: KKcft(:), KKdft(:)
    real*8 :: omega, erfc_width, tau, tau0
    integer :: i

    allocate(KKcft(Nt), KKdft(Nt))
    erfc_width = 15 * lntrange/Nt     ! JW: Could be tweaked
    tau0 = logsbt_lntonset - lnt0

    ! Get analytic logFT of correction kernel
    call logsbt_erfc_logFT(Nt, KKcft, lnt0, lntrange, erfc_width, tau0)
    ! Get smooth product space version of correction kernel
    call logsbt_erfc_prodspace(Nt, KKdft, lnt0, lntrange, erfc_width, tau0)
    call dfftf(Nt, KKdft, wsave)

    ! Add up
    KKbar = KKbar - KKcft + KKdft

    ! Tidy up
    deallocate(KKcft, KKdft)

  end subroutine logsbt_erfc_correction
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_erfc_prodspace
  !  NAME
  !    logsbt_erfc_prodspace
  !  SYNOPSIS

  subroutine logsbt_erfc_prodspace(Nt, KK, lnt0, lntrange, erfc_width, tau0)

    !  PURPOSE
    !    Calculate the product space part of the erfc correction.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: Nt
    real*8, intent(OUT) :: KK(Nt)
    real*8, intent(IN) :: lnt0
    real*8, intent(IN) :: lntrange
    real*8, intent(IN) :: erfc_width, tau0

    !  INPUTS
    !    o Nt -- Number of grid points
    !    o lnt0, lntrange -- Grid is at [exp(lnt0), exp(lnt0+lntrange)[
    !    o erfc_width, tau0 -- erfc paramters
    !  OUTPUTS
    !    o KK -- real-space part of erfc-correction
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: tau
    integer :: i
    character(*), parameter :: func = 'logsbt_erfc_prodspace'

    ! Get smooth kernel on grid (logFT)
    do i = 1, Nt
       tau = (i-1) * lntrange/Nt
       KK(i) = T_erfc(tau, tau0, erfc_width) * lntrange / Nt
    end do

    ! Tidy up and leave
    return

  contains ! ------------------------------

    real*8 function T_erfc(tau, tau0, width)
      use arch_specific, only: arch_erfc
      use constants, only: pi
      implicit none

      real*8, intent(IN) :: tau
      real*8, intent(IN) :: tau0
      real*8, intent(IN) :: width
      !
      T_erfc = arch_erfc((tau-tau0) / sqrt(2.d0) / width) / sqrt(2*pi)
    end function T_erfc

  end subroutine logsbt_erfc_prodspace
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_erfc_logFT
  !  NAME
  !    logsbt_erfc_logFT
  !  SYNOPSIS

  subroutine logsbt_erfc_logFT(Nt, KKbar, lnt0, lntrange, erfc_width, tau0)

    !  PURPOSE
    !    Calculate the logFT part of the erfc correction.
    !  USES

    use constants, only: euler_gamma, pi
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: Nt
    real*8, intent(OUT) :: KKbar(Nt)
    real*8, intent(IN) :: lnt0
    real*8, intent(IN) :: lntrange
    real*8, intent(IN) :: erfc_width, tau0

    !  INPUTS
    !    o Nt -- Number of grid points
    !    o lnt0, lntrange -- Grid is at [exp(lnt0), exp(lnt0+lntrange)[
    !    o erfc_width, tau0 -- erfc paramters
    !  OUTPUTS
    !    o KKbar -- logFT part of erfc-correction in logFT
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    complex*16, allocatable :: KKcmplx(:)
    real*8 :: omega
    integer :: info, i
    character(*), parameter :: func = 'logsbt_erfc_logFT'

    allocate(KKcmplx(0:(Nt+1)/2), stat=info)
    call check_allocation(info, 'KKcmplx', func)

    ! Get analytic logFT of kernel
    KKcmplx(0) = - sqrt(2/pi) * (- lnt0 - tau0 - euler_gamma + 1.d0)
    do i = 1, (Nt+1)/2
       omega = 2*pi * i / lntrange
       KKcmplx(i) = bar_T_erfc(omega, tau0, erfc_width)
    end do
    call logsbt_cmplx2bar(Nt, KKcmplx, KKbar)

    ! Tidy up and leave
    deallocate(KKcmplx)
    return

  contains ! ------------------------------

    complex*16 function bar_T_erfc(omega, tau0, width)
      use constants, only: img_unit
      implicit none

      real*8, intent(IN) :: omega
      real*8, intent(IN) :: tau0
      real*8, intent(IN) :: width
      !
      real*8 :: gauss_val
      complex*16 :: phase
      gauss_val = exp(- 0.5d0 * width**2 * omega**2)
      phase = cmplx(cos(omega*tau0), -sin(omega*tau0), kind(0.d0))
      bar_T_erfc = img_unit * phase * sqrt(2./pi) * gauss_val / omega
    end function bar_T_erfc

  end subroutine logsbt_erfc_logFT
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_vanilla_kernel
  !  NAME
  !    logsbt_vanilla_kernel
  !  SYNOPSIS

  subroutine logsbt_vanilla_kernel(Nt, KKbar, lnt0, lntrange, L, power)

    !  PURPOSE
    !
    !    Generate the kernel KKbar.  This is a logFT
    !
    !       KKbar(i) = Kbar(lnt0 + (i-1) * lntrange/Nt)
    !       Kbar(omega) = \int_\infty^\infty K(tau) e^(-i omega tau) dtau
    !
    !    of the function
    !
    !       K(tau) = sqrt(2/pi) j_L(k0r0 e^tau) (k0r0)**power exp(power*tau).
    !
    !    By discretization in logFT space, the resulting kernel is lntrange
    !    periodic in tau.
    !
    !    By default, lnt0 = lnk0r0 and lntrange = lnrange.  But note that if
    !    both a(rho) and a~(kappa) are is lnrange-periodic, two periods of
    !    K(tau=rho+kappa) are used of the kernel (from lnk0+lnr0 to
    !    lnk0+lnrange+lnr0+lnrange).  This introduces SBT-aliasing, which is
    !    well controlled as long as the transform decays properly to both
    !    sides.  But transformations with L=0 and power=0.0 cannot be done with
    !    a logFT kernel.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    use constants, only: pi
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: Nt        ! grid size
    real*8, intent(OUT) :: KKbar(Nt) ! kernel (see above)
    real*8, intent(IN) :: lnt0       ! grid onset (= lnk0 + lnr0)
    real*8, intent(IN) :: lntrange   ! grid range
    integer, intent(IN) :: L         ! angular momentum
    real*8, intent(IN) :: power      ! incorporated scaling term

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    complex*16, allocatable :: phased_UU(:)
    complex*16 :: phasefac, this_U
    real*8 :: mu, omega, z, offset, phase
    integer :: i, i_start
    character*150 :: info_str

    mu = L + 0.5d0
    offset = power - 1.5d0
    allocate(phased_UU(0:(Nt+1)/2))

    ! Special case bias-free L=0 transform as there is no continuous logFT.
    ! See logsbt_erfc_correction.
    if (mu + 1.d0 + offset < 1d-10) then
       if (L == 0 .and. abs(power) < 1d-10) then
          phased_UU(0) = (0.d0, 0.d0)
          i_start = 1
       else
          write(info_str, "('Invalid SBT parameters: L=',I3,'; power=',F8.4)")&
          & L, power
          call aims_stop(info_str)
       end if
    else
       i_start = 0
    end if

    ! Set up actual logFT kernel
    do i = i_start, (Nt+1)/2
       omega = 2 * pi * i / lntrange
       this_U = Hamilton_U(mu, cmplx(offset, - omega, kind(0.d0)))
       phase = omega * lnt0
       phasefac = cmplx(cos(phase), sin(phase), kind(0.d0))
       phased_UU(i) = this_U * phasefac
    end do
    call logsbt_cmplx2bar(Nt, phased_UU, KKbar)
    return

  contains ! ------------------------------

    complex*16 function Hamilton_U(mu, z)
      implicit none
      real*8, intent(IN) :: mu
      complex*16, intent(IN) :: z
      ! Implementation of U_mu(z) [See Hamilton (2000)]
      !   U_mu(z) = \int_0^\infty t^z J_mu(t) dt
      !           = 2^z * Gamma(0.5*(mu+1+z)) / Gamma(0.5*(mu+1-z))

      real*8 :: mu1_2, ln2
      complex*16 :: lnnumerator, lndenominator
      complex*16 :: cdgamma

      mu1_2 = 0.5d0 * (mu + 1.)
      ln2 = log(2.d0)
      lnnumerator   = cdgamma(mu1_2 + 0.5d0*z, 1)
      lndenominator = cdgamma(mu1_2 - 0.5d0*z, 1)
      Hamilton_U = exp(ln2 * z + lnnumerator - lndenominator)

    end function Hamilton_U
    
  end subroutine logsbt_vanilla_kernel
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_cmplx2bar
  !  NAME
  !    logsbt_cmplx2bar
  !  SYNOPSIS

  subroutine logsbt_cmplx2bar(N, zz, rrbar)

    !  PURPOSE
    !    Transform complex array zz into dfftpack Fourier coefficient format.
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN)    :: N             ! grid size
    complex*16, intent(IN) :: zz(0:(N+1)/2) ! complex Fourier coefficients
    real*8, intent(OUT)    :: rrbar(N)      ! dfftpack Fourier coefficients

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer i

    rrbar(1) = real(zz(0), kind(0.d0))
    do i = 1, N/2-1
       rrbar(2*i) = real(zz(i), kind(0.d0))
       rrbar(2*i+1) = aimag(zz(i))
    end do
    if (modulo(N, 2) == 0) then
       rrbar(N) = real(zz(N/2), kind(0.d0))
    end if

  end subroutine logsbt_cmplx2bar
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_bar2cmplx
  !  NAME
  !    logsbt_bar2cmply
  !  SYNOPSIS

  subroutine logsbt_bar2cmplx(N, zz, rrbar)

    !  PURPOSE
    !    Transform dfftpack Fourier coefficient array to complex
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN)     :: N             ! grid size
    complex*16, intent(OUT) :: zz(0:(N+1)/2) ! complex Fourier coefficients
    real*8, intent(IN)      :: rrbar(N)      ! dfftpack Fourier coefficients

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer i

    zz(0) = rrbar(1)
    do i = 1, N/2-1
       zz(i) = cmplx(rrbar(2*i), rrbar(2*i+1), kind(0.d0))
    end do
    if (modulo(N, 2) == 0) then
       zz(N/2) = rrbar(N)
    end if

  end subroutine logsbt_bar2cmplx
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/logsbt_analyze_tails
  !  NAME
  !    logsbt_analyze_tails
  !  SYNOPSIS

  subroutine logsbt_analyze_tails(N, ff, ln0, lnrange, A, Lr)

    !  PURPOSE
    !
    !    Analyze the tails of ff given on a logarithmic grid.  The function is
    !    assumed to decay polynomial both for small (~ A(1) r^Lr(1)) and large
    !    (~ A(2) r^-Lr(2)) r.
    !
    !    Right now, this procedure is not particularly efficient but more a
    !    proof-of-concept implementation.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: ff(N)
    real*8, intent(IN) :: ln0, lnrange
    real*8, intent(OUT) :: A(2), Lr(2)

    !  INPUTS
    !    o N -- size of grid
    !    o ff -- function on grid
    !    o ln0, lnrange -- onset and range of log grid
    !  OUTPUTS
    !    o A -- coefficient of asymptotic behavior
    !    o Lr -- exponent of asymptotic behavior
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i, i_center, ii(1)
    real*8 :: ref, dln, cur
    character*150 :: info_str

    if (N < 10) call aims_stop('No asymptotics for N < 10')

    dln = lnrange / N

    ! --- Left border

    ref = maxval(abs(ff(1:3)))
    ref = max(1000 * ref, 1d-8)

    ! Search for point which is likely in asymptotic region
    do i = 1, N
       if (abs(ff(i)) > ref) exit
    end do
    if (i > N) then
       ! No asymptotics found
       Lr(1) = 0.d0
       ii = maxloc(abs(ff))
       A(1) = ff(ii(1))
    else
       i_center = i
       Lr(1) = log(abs(ff(i_center) / ff(i_center-1))) / dln
       cur = exp(Lr(1) * (ln0 + (i_center-1)*dln))
       A(1) = ff(i_center) / cur
    end if

    ! --- Right border

    ref = maxval(abs(ff(N-2:N)))
    ref = max(1000 * ref, 1d-10)

    ! Search for point which is likely in asymptotic region
    do i = N, 1, -1
       if (abs(ff(i)) > ref) exit
    end do
    if (i < 1) then
       ! No asymptotics found
       Lr(2) = 0.d0
       ii = maxloc(abs(ff))
       A(2) = ff(ii(1))
    else
       i_center = i
       Lr(2) = - log(abs(ff(i_center) / ff(i_center+1))) / dln
       cur = exp(Lr(2) * (ln0 + (i_center-1)*dln))
       if (cur <= 1d-100) then
          ! With such a high decay coefficient and large integration range, we
          ! are likely to get an underflow here, for a quantity which will
          ! hardly give an error.
          A(2) = 0.d0
       else
          A(2) = ff(i_center) / cur
       end if
    end if

  end subroutine logsbt_analyze_tails
  !******
end module logsbt
!******
