!****s* FHI-aims/evaluate_nuclear_potential
!  NAME
!   evaluate_nuclear_potential
!  SYNOPSIS

subroutine evaluate_nuclear_potential ( dist_tab, &
     n_atom_list, atom_list, potential )

!  PURPOSE
!  Subroutine evaluate_free_atom_sums tabulates the sum of nuclear potentials at a given grid point
!
!  PLEASE NOTE that this subroutine will not trivially work for periodic boundary conditions:
!              The list of atoms passed here will have to be the correct list of atoms contributing
!              to an Ewald potential, NOT the list of atoms whose basis functiona are nonzero.
!  In the cluster case, both lists are the same.
! 
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use pbc_lists
  use species_data
  use spline
  use free_atoms
  use constants
  implicit none

!  ARGUMENTS

  integer:: n_atom_list
  integer:: atom_list(n_atom_list)
  real*8, dimension(n_atom_list) :: dist_tab
  real*8 :: potential


!  INPUTS
!   o n_atom_list -- number of relevant atoms
!   o atom_list -- list of relevant atoms
!   o dist_tab -- distance to atoms
!
!  OUTPUT
!   o potential - superposition of nuclear potentials
!   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE
!



  !     counters
  integer :: i_center_2, i_center_L
  integer :: i_coord
  integer :: i_spin
  integer :: i_atom_2

  !  begin work

  !       initialize
  potential = 0.d0

  !       Contributions from all atoms

  do i_center_L = 1, n_atom_list, 1

     i_center_2 = atom_list(i_center_L)

     if (.not.empty(center_to_atom(i_center_2))) then
        !         electrostatic potential contribution

        ! sign is counted negative as the sign of the electron charge also enters here
        potential = potential - dble (species_z(species_center(i_center_2)) ) / dist_tab(i_center_L)

     end if
  enddo

end subroutine evaluate_nuclear_potential
!---------------------------------------------------------------------
!******
