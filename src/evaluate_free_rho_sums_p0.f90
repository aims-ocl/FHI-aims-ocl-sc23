!****s* FHI-aims/evaluate_free_rho_sums_p0
!  NAME
!   evaluate_free_rho_sums_p0
!  SYNOPSIS

subroutine evaluate_free_rho_sums_p0 ( dist_tab, i_r, &
     free_rho_superpos,n_atom_list, atom_list )

!  PURPOSE
!  Subroutine evaluate_free_atom_sums tabulates the superposition of free-atom
!  densities on the entire integration grid. 
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
  real*8, dimension(n_atom_list) :: i_r
  real*8 :: free_rho_superpos 


!  INPUTS
!   o n_atom_list -- number of relevant atoms
!   o atom_list -- list of relevant atoms
!   o i_r -- tabulates the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!   o dist_tab -- distance to atoms
!
!  OUTPUT
!   o free_rho_superpos -- the superposition of free-atom densities
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
!



  !     counters
  integer :: i_center_2, i_center_L
  integer :: i_coord
  integer :: i_spin
  integer :: i_atom_2

  !  begin work

  !       initialize
  free_rho_superpos = 0.d0

  !       Contributions from all atoms

  do i_center_L = 1, n_atom_list, 1

     i_center_2 = atom_list(i_center_L)

     if (.not.empty(center_to_atom(i_center_2))) then
        !         electrostatic potential contribution

        if (dist_tab(i_center_L).le. &
             multipole_radius_free(species_center(i_center_2))) then

           !            density contribution
           free_rho_superpos = free_rho_superpos +  &
                val_spline &
                ( i_r(i_center_L), free_rho_spl(1,1,species_center(i_center_2)), &
                n_grid(species_center(i_center_2)) )



        end if
     end if
  enddo


end subroutine evaluate_free_rho_sums_p0
!---------------------------------------------------------------------
!******
