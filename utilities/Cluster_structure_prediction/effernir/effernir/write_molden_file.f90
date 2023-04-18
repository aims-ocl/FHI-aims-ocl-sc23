!  Subroutine write_molden_file for visualization of structure
!  with molden
!
!  R.Gehrke (2005)
!

subroutine write_molden_file(coords)

  use vector_array_class
  use cluster

  implicit none

  include "constants.f90"

!  imported variables

!  input

  type (vector_array) :: coords

!  counters
  integer :: i_atom

  write (6,'(2X,A,A)') &
       "Writing molden file for visualization to file structure.xyz ."
  open (50, file = "structure.xyz")

  write (50, '(2X,I4)') n_atoms
  write (50, *)
  
  do i_atom = 1, n_atoms, 1
     
     write (50,'(A, F10.6, 1X, F10.6, 1X, F10.6)') &
          species_name(species(i_atom)), &
          coords%coords(i_atom)%x * bohr, coords%coords(i_atom)%y * bohr, &
          coords%coords(i_atom)%z * bohr
     
  end do
  
  close(50)
  
end subroutine write_molden_file
