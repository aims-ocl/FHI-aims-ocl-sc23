!****s* FHI-aims/get_rlylm
!  NAME
!    get_rlylm
!  SYNOPSIS

subroutine get_rlylm(relvec, Lmax, n_lm, rlylm)

  !  PURPOSE
  !
  !    Calculate spherical harmonics Ylm and derivatives d(Y)/d(rvec)
  !
  !  USES

  use cartesian_ylm
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: relvec(3)
  integer, intent(IN) :: Lmax
  integer, intent(IN) :: n_lm
  real*8, intent(OUT) :: rlylm(n_lm)

  !  INPUTS
  !    o relvec -- vector (x, y, z) defining the direction
  !    o Lmax -- maximum angular momentum for which to calculate
  !    o n_lm -- array size, needs to be >= (Lmax+1)**2
  !  OUTPUTS
  !    o rlylm -- |relvec|^l Y_lm(relvec)
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

  real*8 :: r, rl
  integer :: L, M, i_LM
  character(*), parameter :: func = 'get_rlylm'

  call ylm_real(relvec, Lmax, rlylm)
  r = sqrt(sum(relvec**2))
  rl = 1.d0
  i_LM = 0
  do L = 0, Lmax
     do M = -L, L
        i_LM = i_LM + 1
        rlylm(i_LM) = rl * rlylm(i_LM)
     end do
     rl = rl * r
  end do

end subroutine get_rlylm
!******
