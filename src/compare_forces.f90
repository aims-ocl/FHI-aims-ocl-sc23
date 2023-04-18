!****s* FHI-aims/compare_forces
!  NAME
!    compare_forces
!  SYNOPSIS

subroutine compare_forces(pulay_forces, hellman_feynman_forces, multipole_forces, previous_forces, diff_forces)


!  PURPOSE
!   Small subroutine to evaluate the difference between two
!   set of atomic forces.
!
!  USES
!
  use constants, only: bohr, hartree
  use dimensions
  implicit none
!  ARGUMENTS

  real*8, dimension(3, n_atoms), intent(in) :: pulay_forces
  real*8, dimension(3, n_atoms), intent(in) :: hellman_feynman_forces
  real*8, dimension(3, n_atoms), intent(in) :: multipole_forces 
  real*8, dimension(3, n_atoms), intent(inout) :: previous_forces
  real*8, intent(out) :: diff_forces

!  INPUTS
!  o pulay_forces - Pulay forces
!  o hellman_feynman_forces - Hellman Feyman forces
!  o multipole_forces - Multipole forces
!  o previous_forces - Forces from previous iteration
!
!  OUTPUT
!   o diff_forces - difference between forces now and previous iterations
!   o previous_forces - Forces from this iteration
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




  ! local variables

  ! these are the forces used for the convergence check. They exclude the
  ! (expensive) GGA forces.
  real*8, dimension(3, n_atoms) :: convergence_forces

  real*8 :: conversion

  ! counters
  integer :: i_atom
  integer :: i_coords

  conversion = hartree / bohr

  diff_forces = 0.d0
  do i_atom = 1, n_atoms, 1
     do i_coords = 1, 3, 1

        convergence_forces(i_coords,i_atom) = &
            hellman_feynman_forces(i_coords, i_atom) &
          + multipole_forces(i_coords, i_atom)

        diff_forces = diff_forces + &
          (previous_forces(i_coords, i_atom) - convergence_forces(i_coords, i_atom)) ** 2

     end do
  end do
  ! norm difference to number of atoms
  diff_forces = diff_forces / n_atoms
  diff_forces = sqrt(diff_forces)

  diff_forces = diff_forces * conversion

  ! and store present forces as previous forces for next iteration
  do i_atom=1, n_atoms, 1
     do i_coords=1, 3, 1
        previous_forces(i_coords, i_atom) = &
           convergence_forces(i_coords, i_atom)
     end do
  end do

end subroutine compare_forces
!******
