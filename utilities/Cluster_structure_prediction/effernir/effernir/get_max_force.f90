real*8 function get_max_force(forces, n_atoms)
  
  use vector_array_class

  implicit none

  !  imported variables
  
  !  input
  type (vector_array), intent(in) :: forces
  integer, intent(in) :: n_atoms
  
  !  local variables
  real*8 :: max_force
  
  !  counter
  integer :: i_atom
  
  max_force = 0.d0
  do i_atom = 1, n_atoms, 1
     if (abs(forces%coords(i_atom)%x) .gt. max_force) then
        max_force = abs(forces%coords(i_atom)%x)
     end if
     if (abs(forces%coords(i_atom)%y) .gt. max_force) then
        max_force = abs(forces%coords(i_atom)%y)
     end if
     if (abs(forces%coords(i_atom)%z) .gt. max_force) then
        max_force = abs(forces%coords(i_atom)%z)
     end if
  end do
  
  get_max_force = max_force
  
end function get_max_force
