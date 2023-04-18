!  subroutine to check whether atoms are still
!  bound to a cluster
!
!  R. Gehrke (2007)
!

subroutine check_compound(coords, d_max, dissociated)

  use vector_array_class
  use cluster

  implicit none

  ! imported variables

  ! input
  type (vector_array), intent(in) :: coords 
  real*8, intent(in) :: d_max
  
  ! output
  logical, intent(out) :: dissociated

  ! local variables
  integer, dimension(n_atoms) :: old_atoms
  integer, dimension(n_atoms) :: new_atoms
  integer, dimension(n_atoms) :: found_atoms
  integer :: n_new_found
  integer :: n_old_found
  integer :: n_found
  real*8 :: distance
  logical :: found

  !  counter
  integer :: i_atom
  integer :: i_counter
  integer :: i_counter_2

  dissociated = .true.
  n_new_found = 1
  found_atoms(:) = 0
  found_atoms(1) = 1 ! atom # 1 !!!!
  old_atoms(:) = 0
  new_atoms(:) = 0
  do i_atom = 2, n_atoms, 1
     distance = norm(coords%coords(i_atom) - coords%coords(1))
     if (distance .lt. d_max) then
        new_atoms(n_new_found) = i_atom
        n_new_found = n_new_found + 1
        found_atoms(n_new_found) = i_atom
     end if
  end do
  n_found = n_new_found
  n_new_found = n_new_found - 1
  do while (n_new_found .gt. 0)
     old_atoms(:) = new_atoms(:)
     new_atoms(:) = 0
     n_old_found = n_new_found
     n_new_found = 0
     do i_counter = 1, n_old_found, 1
        do i_atom = 1, n_atoms, 1
           if (i_atom .ne. old_atoms(i_counter)) then
              distance = norm(coords%coords(i_atom) - coords%coords(old_atoms(i_counter))) 
           end if
           if (distance .lt. d_max) then
              found = .false.
              do i_counter_2 = 1, n_found, 1
                 if (i_atom .eq. found_atoms(i_counter_2)) then
                    found = .true.
                 end if
              end do
              if (.not.found) then
                 n_new_found = n_new_found + 1
                 new_atoms(n_new_found) = i_atom
                 n_found = n_found + 1
                 found_atoms(n_found) = i_atom
              end if
           end if
        end do
     end do
  end do
  ! all atoms found?
  if (n_found .eq. n_atoms) then
     dissociated = .false.
  end if
  
end subroutine check_compound
