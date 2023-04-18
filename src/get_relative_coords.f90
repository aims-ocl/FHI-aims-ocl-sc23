! calculates relative coordinates with respect to the center of mass
! masses are set to unity
!
! R.Gehrke (2007)

subroutine get_relative_coords(coords, center_of_mass, relative_coords)

  use dimensions
  
  implicit none

  ! input
  real*8, dimension(3, n_atoms), intent(in) :: coords

  ! output
  real*8, dimension(3), intent(out) :: center_of_mass
  real*8, dimension(3, n_atoms), intent(out) :: relative_coords

  ! counter
  integer :: i_atom

  center_of_mass = 0.d0
  do i_atom = 1, n_atoms, 1
     center_of_mass(:) = center_of_mass(:) + coords(:, i_atom)
  end do
  center_of_mass = center_of_mass / dble(n_atoms)
  do i_atom = 1, n_atoms, 1
     relative_coords(:, i_atom) = coords(:, i_atom) - center_of_mass(:)
  end do

end subroutine get_relative_coords
