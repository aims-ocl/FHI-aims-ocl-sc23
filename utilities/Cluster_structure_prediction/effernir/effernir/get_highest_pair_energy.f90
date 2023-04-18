subroutine get_highest_pair_energy(current_coords, pair_energy_max, i_atom)

  use vector_array_class
  use cluster

  implicit none
  
  ! imported variables

  ! input
  type (vector_array), intent(in) :: current_coords

  ! output
  real*8, intent(out) :: pair_energy_max
  integer, intent(out) :: i_atom

  ! local variables
  real*8 :: pair_energy

  ! counter
  integer :: i_atom_2

  pair_energy_max = -1e10
  write (6,*) "pair energies:"
  i_atom = 1
  do i_atom_2 = 1, n_atoms, 1
     call get_pair_energy(current_coords, i_atom_2, pair_energy)
     write (6,*) "# ", i_atom_2, pair_energy
     if (pair_energy .gt. pair_energy_max) then
        i_atom = i_atom_2
        pair_energy_max = pair_energy
     end if
  end do
  write (6,*) "atom with highest pair energy:", i_atom

end subroutine get_highest_pair_energy
