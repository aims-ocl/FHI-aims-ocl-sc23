subroutine pre_relax_hard_sphere(coords, d_min, n_atoms)

  use vector_array_class
  use control

  implicit none

  ! imported variables
  
  ! input
  type (vector_array), intent(inout) :: coords
  integer, intent(in) :: n_atoms
  real*8, intent(in) :: d_min

  ! local variables
  type (vector) :: distance_vector
  type (vector) :: step_vector
  real*8 :: scaling_factor
  real*8 :: distance
  logical :: found
  
  ! counter
  integer :: i_atom
  integer :: i_atom_2

  found = .true.
  do while (found) 
     found = .false.
     do i_atom = 1, n_atoms, 1
        do i_atom_2 = i_atom + 1, n_atoms, 1
           
           distance_vector = coords%coords(i_atom_2) - coords%coords(i_atom)
           distance = norm(distance_vector)
           scaling_factor = 0.1 / distance
           step_vector = distance_vector * scaling_factor
           
           if (distance .lt. d_min) then 
              found = .true.
              coords%coords(i_atom)   = coords%coords(i_atom)   - step_vector
              coords%coords(i_atom_2) = coords%coords(i_atom_2) + step_vector
              if (verbose) then
                 write (150,'(A, I4, I4)') "displacement of atoms # ", i_atom, i_atom_2
              end if
           end if
           
        end do
     end do
  end do
  
end subroutine pre_relax_hard_sphere
