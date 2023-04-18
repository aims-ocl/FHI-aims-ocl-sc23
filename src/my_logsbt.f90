!****h* FHI-aims/logsbt
!  NAME
!    logsbt
!  SYNOPSIS

module my_logsbt

  implicit none

  real*8, parameter :: logsbt_lntonset = 0.d0 ! center of erfc correction

  integer, parameter :: i_logft_aliasing = 1
  integer, parameter :: i_sbt_aliasing_smallk = 2
  integer, parameter :: i_sbt_aliasing_largek = 3
  integer, parameter :: i_sbt_ringing_smallr = 4
  integer, parameter :: i_sbt_ringing_larger = 5
  integer, parameter :: n_sbt_errors = 5
  contains

  subroutine my_logsbt_multi_driver(N, lnr0, lnk0, lnrange, power_bias, &
  &                              n_fn, fn_to_L, ff, sbt_errors)

    !  PURPOSE
    !
    !    Perform SBT on n_fn radial parts in one sweep.  Input and output ff
    !    needs to be scaled.
    !
    !  USES
    use logsbt, only : logsbt_expert, logsbt_kernel, logsbt_opt_trange
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

    ! omp does not work here
!!!!$omp parallel do   &
!!!!$omp default(none) &
!!!!$omp private(i_fn, L, these_errors) &
!!!!$omp shared(n_fn, fn_to_L, max_L, N, ff, lnr0, lnk0, lnrange, power_bias, &
!!!!$omp        NT, KKbar, lntrange, wsave, sbt_errors) &
!!!!$omp schedule(static)

    print *,"multi_driver= ",n_fn
    do i_fn = 1, n_fn

       L = fn_to_L(i_fn)
       if (L > max_L .or. L < 0) call aims_stop('Invalid L', func)
       call logsbt_expert(N, ff(:, i_fn), &
       &                  lnr0, lnk0, lnrange, L, power_bias, &
       &                  Nt(L), KKbar(:, L), lntrange(L), wsave(:, L), &
       &                  these_errors)
       if (present(sbt_errors)) sbt_errors(:, i_fn) = these_errors
    end do
!!!!!$omp end parallel do


    ! ff == f(k)*k**power_bias       ! [error: epsilon]
    ! call logsbt_scale(N, n_fn, ff, lnk0, lnrange, power_add)
    ! ff == f(k)*k**(p_bias+p_add)   ! [error: epsilon * k**power_add]

    ! --- Tidy up

    deallocate(Nt, wsave, KKbar)

  end subroutine my_logsbt_multi_driver
 end module my_logsbt
!******
