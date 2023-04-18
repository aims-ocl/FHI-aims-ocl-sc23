!----------------------------------------------------------------------------
!****s* FHI-aims/hse_logsbt_integrator
!  NAME
!    hse_logsbt_integrator
!  SYNOPSIS

subroutine hse_logsbt_integrator(L, N, ff, lnr0, lnk0, lnrange, omega)

  !  PURPOSE
  !
  !    Integrate HSE field by multiplication in Fourier space using logSBT.
  !    This should be done on an SBT grid, which is in general wider and
  !    denser than the ordinary log-grid from grids.f90.  If your function is
  !    given on the latter kind of grid, please use the wrapper
  !    hse_logsbt_integrator_grid() instead.
  !
  !  USES

  use constants
  use logsbt
  use prodbas
  use mpi_tasks, only: aims_stop, check_allocation
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: L
  integer, intent(IN) :: N
  real*8, intent(INOUT) :: ff(N)
  real*8, intent(IN) :: lnr0, lnk0, lnrange
  real*8, intent(IN) :: omega

  !  INPUTS
  !    o L -- current angular momentum channel
  !    o ff -- wave function f(r)
  !    o N, lnr0, lnk0, lnrange -- log-grid specs
  !    o omega -- decay parameter in error function
  !  OUTPUTS
  !    o ff -- Screened coulomb field of wave on log grid [V(r)]
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8 :: dln
  integer :: i
  real*8 :: k
  real*8 :: exp_arg, hse_fac, sbt_errors(n_sbt_errors)
  real*8, external :: radial_fourier_hse
  
  ! The HSE kernel is constant (pi/omega**2) for k->0 and decays k^-2 for
  ! k->oo and can therefore be applied without any problems.
  real*8, parameter :: power_bias = 1.5d0

  character(*), parameter :: func = 'hse_logsbt_integrator'

  ! JW: Right now, the logsbt kernels are reconstructed for every single
  !     transform, which is more expensive than necessary.  It would be
  !     beneficial to use logsbt_multi_driver() and do all needed transforms
  !     in one sweep.

  dln = lnrange / N
  ! ff == f(r)                     [assume full relative accuracy]

  call logsbt_scale(N, 1, ff, lnr0, lnrange, 3.d0 - power_bias)
  ! ff == f(r)**(3 - power_bias)   [assume full relative accuracy]

  call logsbt_driver(N, ff, lnr0, lnk0, lnrange, L, power_bias, sbt_errors)
  ! ff == f~(k)**power_bias        [error: epsilon]
  ! if (myid == 0) write(0,"('r->k: ', 5ES10.2)") sbt_errors

  call logsbt_scale(N, 1, ff, lnk0, lnrange, 3.d0 - power_bias)
  ! ff == f~(k)**3                 [error: epsilon * k**1.5]
  do i = 1, N
     k = logsbt_ilog2r(dble(i), lnk0, dln)
     ff(i) = ff(i) * radial_fourier_hse(k, omega)
  end do
  ! ff == V~(k)**3                 [error k->0:  epsilon * pi/omega**2 * k**1.5
  !                                       k->oo: epsilon / k**0.5] <= epsilon.

  call logsbt_driver(N, ff, lnk0, lnr0, lnrange, L, 0.d0, sbt_errors)
  ! ff = V(r)                      [error: epsilon]
  ! if (myid == 0) write(0,"('k->r: ',5ES10.2)") sbt_errors
end subroutine hse_logsbt_integrator
!******
