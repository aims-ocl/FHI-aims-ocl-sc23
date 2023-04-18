 real*8 function get_line_force(forces, search_direction, n_atoms)
  
  use vector_array_class

  implicit none

  !  imported variables
  
  !  input
  type (vector_array), intent(in) :: forces
  type (vector_array), intent(in) :: search_direction
  integer, intent(in) :: n_atoms
  
  ! local variables
  real*8 :: lineforce

  !  counter
  integer :: i_atom
  
  lineforce = 0.d0
  do i_atom = 1, n_atoms, 1
     lineforce = lineforce - forces%coords(i_atom)%x * search_direction%coords(i_atom)%x
     lineforce = lineforce - forces%coords(i_atom)%y * search_direction%coords(i_atom)%y
     lineforce = lineforce - forces%coords(i_atom)%z * search_direction%coords(i_atom)%z
  end do
  get_line_force = lineforce

end function get_line_force
