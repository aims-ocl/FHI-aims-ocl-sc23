! calculates net torque on the center of mass on a system.
! in a cluster, net torque should ideally vanish
! dft-forces are noisy, so it doesn't
! noise is removed by subroutine "remove_translation_and_rotation"
!
! R.Gehrke (2007)

subroutine get_net_torque(total_forces, relative_coords, net_torque)

  use dimensions
  
  implicit none

  ! imported variables

  ! input
  real*8, dimension(3, n_atoms), intent(in) :: relative_coords
  real*8, dimension(3, n_atoms), intent(in) :: total_forces
  
  ! output
  real*8, dimension(3), intent(out) :: net_torque

  ! counter
  integer :: i_atom

  net_torque = 0.d0
  do i_atom = 1, n_atoms, 1
     net_torque(1) = net_torque(1) + relative_coords(2, i_atom) * total_forces(3, i_atom) - &
                                     relative_coords(3, i_atom) * total_forces(2, i_atom)
     net_torque(2) = net_torque(2) + relative_coords(3, i_atom) * total_forces(1, i_atom) - &
                                     relative_coords(1, i_atom) * total_forces(3, i_atom)
     net_torque(3) = net_torque(3) + relative_coords(1, i_atom) * total_forces(2, i_atom) - &
                                     relative_coords(2, i_atom) * total_forces(1, i_atom)
  end do

end subroutine get_net_torque
