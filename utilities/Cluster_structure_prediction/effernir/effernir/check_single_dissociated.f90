! subroutine to check whether a single atom with index i_atom
! is dissociated (larger than d_max away from all other atoms)
!
! R. Gehrke (2009)
!

subroutine check_single_dissociated(coords, d_max, n_atoms, index_atom, dissociated)

  use vector_array_class

  implicit none

  ! imported variables

  ! input
  type (vector_array), intent(in) :: coords 
  real*8, intent(in) :: d_max
  integer, intent(in) :: n_atoms
  integer, intent(in) :: index_atom
  
  ! output
  logical, intent(out) :: dissociated

  ! local variables
  real*8 :: distance
  
  ! counter
  integer :: i_atom

  dissociated = .true.
  do i_atom = 1, n_atoms, 1
     if (i_atom .ne. index_atom) then
        distance = norm(coords%coords(i_atom) - coords%coords(index_atom))
        if (distance .lt. d_max) then
           dissociated = .false.
           exit
        end if
     end if
  end do

end subroutine check_single_dissociated
