subroutine get_pair_energy(coords, i_atom, pair_energy)
  
  use cluster
  use pair_potential
  use lennard_jones
  use vector_array_class

  implicit none

  ! imported variables

  ! input
  type (vector_array), intent(in) :: coords
  integer, intent(in) :: i_atom

  ! output
  real*8, intent(out) :: pair_energy

  ! local variables
  real*8 :: energy_contribution
  real*8 :: r_index
  real*8 :: distance

  ! counter
  integer :: i_atom_2

  pair_energy = 0.d0
!  write (6,*) "contribution:"
  do i_atom_2 = 1, n_atoms, 1
     if (i_atom .ne. i_atom_2) then
        distance = norm(coords%coords(i_atom) - coords%coords(i_atom_2))
        select case (potential_flag)

        case ('external')
           if (distance .lt. d_max(species(i_atom), species(i_atom_2))) then
              r_index = (distance - d_min(species(i_atom), species(i_atom_2))) / delta_d(species(i_atom), species(i_atom_2)) + 1
              energy_contribution = val_spline(r_index, &
                   pot_spl(1,1,species(i_atom),species(i_atom_2)), n_data(species(i_atom),species(i_atom_2)))
           else
              energy_contribution = 0.d0
           end if
        case ('LJ')
           energy_contribution = get_single_lj_energy(distance, i_atom, i_atom_2)
           
        case default
           write (6,*) "Internal inconsistency!"
           write (6,*) "(get_pair_energy)"
           write (6,*) "* Aborting."
           stop
           
        end select
        pair_energy = pair_energy + energy_contribution
        
!        write (6,*) "# ", i_atom_2, " ", energy_contribution
     end if
  end do
  
end subroutine get_pair_energy
