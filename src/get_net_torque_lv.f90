! calculates net torque on the center of mass on a system of lattice vectors
! net torque should ideally vanish
! dft-forces are noisy, so it doesn't
! noise is removed by subroutine "remove_translation_and_rotation"
!
! R.Gehrke (2007)

subroutine get_net_torque_lv(forces_lv, relative_lv, net_torque_lv)

  use dimensions
  
  implicit none

  ! imported variables

  ! input
  real*8, dimension(3, n_periodic), intent(in) :: relative_lv
  real*8, dimension(3, n_periodic), intent(in) :: forces_lv
  
  ! output
  real*8, dimension(3), intent(out) :: net_torque_lv

  ! counter
  integer :: i_periodic

  net_torque_lv = 0.d0
  do i_periodic = 1, n_periodic, 1
     net_torque_lv(1) = net_torque_lv(1) + relative_lv(2, i_periodic) * forces_lv(3, i_periodic) - &
                                     relative_lv(3, i_periodic) * forces_lv(2, i_periodic)
     net_torque_lv(2) = net_torque_lv(2) + relative_lv(3, i_periodic) * forces_lv(1, i_periodic) - &
                                     relative_lv(1, i_periodic) * forces_lv(3, i_periodic)
     net_torque_lv(3) = net_torque_lv(3) + relative_lv(1, i_periodic) * forces_lv(2, i_periodic) - &
                                     relative_lv(2, i_periodic) * forces_lv(1, i_periodic)
  end do

end subroutine get_net_torque_lv
