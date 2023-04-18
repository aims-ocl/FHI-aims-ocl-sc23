!****s* FHI-aims/get_net_force_lv
!  NAME
!   get_net_force_lv
!  SYNOPSIS

subroutine get_net_force_lv(forces_lv, net_force_lv)

! PURPOSE
! Calculates total forces on  lattice vector
! dft-forces are noisy, so it doesn't
! noise is removed by subroutine "remove_translation_and_rotation"
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  real*8, dimension(3, n_atoms), intent(in) :: forces_lv
  real*8, dimension(3), intent(out) :: net_force_lv

!  INPUTS
!   o forces_lv -- total forces on lattice vectors
!  OUTPUT
!   o net_force_lv -- total sum of forces
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



  ! counter
  integer :: i_periodic

  net_force_lv = 0.d0
  do i_periodic = 1, n_periodic, 1
     net_force_lv(:) = net_force_lv(:) + forces_lv(:, i_periodic)
  end do

end subroutine get_net_force_lv
!******
