!****s* FHI-aims/get_rlylm_and_derivs
!  NAME
!    get_rlylm_and_derivs
!  SYNOPSIS

subroutine get_rlylm_and_derivs(relvec, Lmax, n_lm, rlylm, drlylm_drelvec)

  !  PURPOSE
  !
  !    Calculate spherical harmonics Ylm and derivatives d(Y)/d(rvec)
  !
  !  USES

  use cartesian_ylm
  use mpi_tasks
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: relvec(3)
  integer, intent(IN) :: Lmax
  integer, intent(IN) :: n_lm
  real*8, intent(OUT) :: rlylm(n_lm)
  real*8, intent(OUT) :: drlylm_drelvec(3, n_lm)

  !  INPUTS
  !    o relvec -- vector (x, y, z) defining the direction
  !    o Lmax -- maximum angular momentum for which to calculate
  !    o n_lm -- array size, needs to be >= (Lmax+1)**2
  !    o want_drylm -- do we need derivatives?
  !  OUTPUTS
  !    o rlylm -- |relvec|^l Y_lm(relvec)
  !    o drlylm_drelvec -- deriv of rlylm
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  real*8, allocatable :: cartesians(:,:)
  integer :: info
  character(*), parameter :: func = 'get_rlylm_and_derivs'

  call initialize_cartesian_ylm(Lmax)   ! Fast exit if not needed.
  allocate(cartesians(n_max_cartesian, 0:Lmax), stat=info)
  call check_allocation(info, 'cartesians', func)
  call evaluate_onecenter_cartesians(relvec, Lmax, cartesians)
  call evaluate_onecenter_cartesian_gradient_terms(Lmax, cartesians, &
  &                                                rlylm, drlylm_drelvec)
  deallocate(cartesians)

end subroutine get_rlylm_and_derivs
!******
