!****s* FHI-aims/get_net_force
!  NAME
!   get_net_force
!  SYNOPSIS

subroutine get_net_force(total_forces, net_force)

! PURPOSE
! Calculates net force on the system.
! in a cluster, net force should ideally vanish
! dft-forces are noisy, so it doesn't
! noise is removed by subroutine "remove_translation_and_rotation"
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS

  real*8, dimension(3, n_atoms), intent(in) :: total_forces
  real*8, dimension(3), intent(out) :: net_force

!  INPUTS
!   o total_forces -- total forces of atoms
!  OUTPUT
!   o net_force -- total sum of forces
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
  integer :: i_atom

  net_force = 0.d0
  do i_atom = 1, n_atoms, 1
     net_force(:) = net_force(:) + total_forces(:, i_atom)
  end do

end subroutine get_net_force
!******
