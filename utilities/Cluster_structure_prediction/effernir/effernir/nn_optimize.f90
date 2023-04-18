! wrapper around the external geometry optimization based upon the
! NN-code by J.Behler
!
! R.Gehrke

subroutine nn_optimize(try_coords, optimized_coords, energy, nn_atomic_energy, status)

  use vector_array_class
  use cluster

  implicit none

  include "constants.f90"

  ! imported variables

  ! input
  type (vector_array), intent(in) :: try_coords
   
  ! output
  type (vector_array), intent(out) :: optimized_coords
  real*8, intent(out) :: energy
  real*8, intent(out) :: nn_atomic_energy(n_atoms)
  integer, intent(inout) :: status

  ! local variables
  character*60 cmdtxt
  character*20 dummy

  ! counter
  integer :: i_atom

  ! write input file for external optimization code (in units of bohr)
  open (75, FILE="structure.dat")
  
  do i_atom = 1, n_atoms, 1
     write (75,'(A10,2X,3F14.6)') species_name(species(i_atom)), try_coords%coords(i_atom)%x * inv_bohr, try_coords%coords(i_atom)%y * inv_bohr, &
          try_coords%coords(i_atom)%z * inv_bohr
  end do

  close(75)
  ! call external optimization code
  cmdtxt = './geooptNN.x'
  call system(TRIM(cmdtxt))
 
  ! read output structure
  open (75, FILE="final.dat")
  
  do i_atom = 1, n_atoms, 1
     read (75,*) dummy, optimized_coords%coords(i_atom)%x, optimized_coords%coords(i_atom)%y, optimized_coords%coords(i_atom)%z
     optimized_coords%coords(i_atom)%x = optimized_coords%coords(i_atom)%x * bohr
     optimized_coords%coords(i_atom)%y = optimized_coords%coords(i_atom)%y * bohr
     optimized_coords%coords(i_atom)%z = optimized_coords%coords(i_atom)%z * bohr
  end do

  close(75)

  ! read energy
  open (75, FILE="energy.out")
  read(75,*) energy
  close(75)

  energy = energy * hartree

  ! read atomic energies
  open (75, FILE="atomenergy.out")
  do i_atom = 1, n_atoms, 1
     read (75,*) dummy, nn_atomic_energy(i_atom)  
     nn_atomic_energy(i_atom) = nn_atomic_energy(i_atom) * hartree - atomic_energy(species(i_atom))
     write (150,*) nn_atomic_energy(i_atom)
  end do
  close(75)

  status=0

end subroutine nn_optimize
