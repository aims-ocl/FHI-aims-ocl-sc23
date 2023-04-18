!----------------------------------------------------------------------------
!****s* FHI-aims/hse_logsbt_integrator_grid
!  NAME
!    hse_logsbt_integrator_grid
!  SYNOPSIS

subroutine hse_logsbt_integrator_grid(i_l, wave, &
&                                     n_grid, r_grid_min, r_grid_inc, &
&                                     omega, field)

  !  PURPOSE
  !
  !    Integrate HSE field by multiplication in Fourier space using logSBT.
  !
  !  USES

  use runtime_choices, only: sbtgrid_N, sbtgrid_lnr0, &
  &                          sbtgrid_lnk0, sbtgrid_lnrange
  use constants
  use dimensions
  use bspline
  use logsbt
  use sbt_overlap
  use mpi_tasks, only: aims_stop, check_allocation
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_l
  real*8, intent(IN) :: wave(n_max_spline, n_max_grid)
  integer, intent(IN) :: n_grid
  real*8, intent(IN) :: r_grid_min, r_grid_inc
  real*8, intent(IN) :: omega
  real*8, intent(OUT) :: field(n_max_grid)

  !  INPUTS
  !    o i_l -- current angular momentum channel
  !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
  !    o n_grid -- number of logarithmic grid points
  !    o r_grid -- logrithmic grid
  !    o omega -- decay parameter in error function
  !  OUTPUTS
  !    o field -- Screened coulomb field of wave on log grid [V(r)]
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: N
  real*8 :: lnr0, lnk0, lnrange, lnrN, r_grid_max
  real*8 :: log_min, log_inc
  real*8, allocatable :: ff(:), ff_spl(:)
  real*8 :: lnr, val(1), i_r
  integer :: info
  integer :: i
  character(*), parameter :: func = 'hse_logsbt_integrator_grid'

  ! --- init

  N = sbtgrid_N
  lnr0 = sbtgrid_lnr0
  lnk0 = sbtgrid_lnk0
  lnrange = sbtgrid_lnrange
  r_grid_max = r_grid_min*r_grid_inc**(n_grid-1)
  lnrN = lnr0 + lnrange*dble(N-1)/dble(N)

  ! call aims_stop('UNTESTED', func)

  if (n_grid <= 1) call aims_stop('Non-existent grid', func)
  if (r_grid_min < exp(lnr0) .or. r_grid_max > exp(lnrN)) then
     call aims_stop('sbtgrid is smaller than usual log-grid', func)
  end if
  allocate(ff(N), ff_spl(0:N+1), stat=info)
  call check_allocation(info, 'ff, ff_spl', func)

  ! --- import

  call sbt_import_spline(N, ff, lnr0, lnrange, i_l, &
  &                      n_grid, wave, r_grid_min, r_grid_inc)

  ! --- integrate

  call hse_logsbt_integrator(i_l, N, ff, lnr0, lnk0, lnrange, omega)

  ! --- export

  call cubic_bspline_notaknot(N, 1, ff, ff_spl)

  field = 0.d0
  log_min = log(r_grid_min)
  log_inc = log(r_grid_inc)
  do i = 1, n_grid
     lnr = log_min + (i-1) * log_inc
     call val_scaled_bspline(N, 1, lnr0, lnrN, val, lnr, ff_spl)
     field(i) = val(1)
  end do

  deallocate(ff, ff_spl)

end subroutine hse_logsbt_integrator_grid
!******
