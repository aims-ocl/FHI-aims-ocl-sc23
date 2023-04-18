subroutine get_forces_and_energy(coords, energy, forces, status)

  use vector_array_class
  use control
  use lennard_jones

  implicit none
  
  !  imported variables
  
  !  input
  type (vector_array), intent(in) :: coords
  
  !  output
  type (vector_array), intent(out) :: forces
  real*8, intent(out) :: energy
  integer, intent(inout) :: status

  select case (potential_flag)

  case ('external')
     write (6,'(2X,A)') "Internal inconsistency."
     write (6,'(2X,A)') "(get_forces_and_energy)"
     write (6,'(2X,A)') "Internal minimizer is currently not used in connection with external code."
     write (6,'(2X,A)') "* Aborting."
     stop

  case ('LJ')
     call get_lj_forces_and_energy(coords, forces, energy)

  case default
     write (6,'(2X,A)') "Internal inconsistency."
     write (6,'(2X,A)') "(get_forces_and_energy)"
     write (6,'(2X,A)') "Chosen interaction not implemented."
     write (6,'(2X,A)') "* Aborting."
     stop

  end select
     
end subroutine get_forces_and_energy
