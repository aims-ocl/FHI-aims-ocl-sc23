!  module to map an energy surface on one coordinate
!  according to a search direction
!
!  R. Gehrke (2005)
!

module map_to_one_dimension

use vector_array_class
use cluster
use dft_file

implicit none

type (vector_array) :: zero_coords
type (vector_array) :: search_direction
type (vector_array) :: coords

contains
  subroutine allocate(n_atoms)

! imported variables

! input
    integer, intent(in) :: n_atoms

    call allocate_array(search_direction, n_atoms)
    call allocate_array(zero_coords, n_atoms)
    call allocate_array(coords, n_atoms)
  
  end subroutine allocate
  
  subroutine deallocate()

    call deallocate_array(search_direction)
    call deallocate_array(zero_coords)
    call deallocate_array(coords)

  end subroutine deallocate

  subroutine evaluate(param, forces, lineforce, energy, cosine, max_force, status)

    ! imported variables
    
    ! input
    real*8, intent(in) :: param

    ! output
    type (vector_array), intent(out) :: forces
    real*8, intent(out) :: lineforce
    real*8, intent(out) :: energy
    real*8, intent(out) :: cosine
    real*8, intent(out) :: max_force
    integer :: status

    ! functions
    real*8 :: get_max_force

    ! counter
    integer :: i_atom
    do i_atom = 1, n_atoms, 1
       coords%coords(i_atom)%x = zero_coords%coords(i_atom)%x + param * search_direction%coords(i_atom)%x
       coords%coords(i_atom)%y = zero_coords%coords(i_atom)%y + param * search_direction%coords(i_atom)%y
       coords%coords(i_atom)%z = zero_coords%coords(i_atom)%z + param * search_direction%coords(i_atom)%z
    end do
    call get_forces_and_energy(coords, energy, forces, status)

    if (status .eq. 0) then
       ! calculate minus times line_force so that line_force is real directional derivative
       lineforce = 0.d0
       do i_atom = 1, n_atoms, 1
          lineforce = lineforce - forces%coords(i_atom)%x * search_direction%coords(i_atom)%x
          lineforce = lineforce - forces%coords(i_atom)%y * search_direction%coords(i_atom)%y
          lineforce = lineforce - forces%coords(i_atom)%z * search_direction%coords(i_atom)%z
       end do
       cosine = lineforce / (array_norm(forces) * array_norm(search_direction))
       max_force = get_max_force(forces, n_atoms)
    end if
    
  end subroutine evaluate

end module map_to_one_dimension
