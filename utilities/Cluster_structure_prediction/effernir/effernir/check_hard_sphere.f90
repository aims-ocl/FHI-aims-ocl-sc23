! subroutine check_hard_sphere to check whether two atoms come
! to close together
!
! R.Gehrke (2006)
!

subroutine check_hard_sphere(coords, n_atoms, reasonable)

  use control
  use vector_array_class

  implicit none

! imported variables

! input
  type (vector_array), intent(in) :: coords
  integer, intent(in) :: n_atoms

! output
  logical, intent(out) :: reasonable

! local variables
  real*8 :: min_distance
  real*8 :: distance

! counter
  integer :: i_atom
  integer :: i_atom_2

  min_distance = 1e5
  do i_atom = 1, n_atoms, 1
     do i_atom_2 = i_atom + 1, n_atoms, 1
        distance = norm(coords%coords(i_atom) - coords%coords(i_atom_2))
        if (distance .lt. min_distance) then 
           min_distance = distance
        end if
     end do
  end do

  if (min_distance .lt. hard_sphere_radius) then
     reasonable = .false.
  else
     reasonable = .true.
  end if

end subroutine check_hard_sphere
