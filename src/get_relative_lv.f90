! calculates relative coordinates with respect to the center of mass
! masses are set to unity
!
! R.Gehrke (2007)

subroutine get_relative_lv(lv, center_of_mass_lv, relative_lv)

  use dimensions
  
  implicit none

  ! input
  real*8, dimension(3, n_periodic), intent(in) :: lv 

  ! output
  real*8, dimension(3), intent(out) :: center_of_mass_lv
  real*8, dimension(3, n_periodic), intent(out) :: relative_lv

  ! counter
  integer :: i_periodic

  center_of_mass_lv = 0.d0
  do i_periodic = 1, n_periodic, 1
     center_of_mass_lv(:) = center_of_mass_lv(:) + lv(:, i_periodic)
  end do
  center_of_mass_lv = center_of_mass_lv / dble(n_periodic)
  do i_periodic = 1, n_periodic, 1
     relative_lv(:, i_periodic) = lv(:, i_periodic) - center_of_mass_lv(:)
  end do

! For the time beeing - do it with respect to origin:
  center_of_mass_lv = 0.d0
  relative_lv = lv


end subroutine get_relative_lv
