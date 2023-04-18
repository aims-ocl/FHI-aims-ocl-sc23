!****s* FHI-aims/tab_single_atom_centered_coords_radial_log_p0
!  NAME
!   tab_single_atom_centered_coords_radial_log_p0
!  SYNOPSIS

subroutine tab_single_atom_centered_coords_radial_log_p0 &
     ( current_center, dist_tab_sq, dir_tab,  &
     dist_tab_in, i_r, i_r_log, dir_tab_in )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates all quantitied
!  needed to access the splined Hartree multipole potential for
!  a given group of atoms
!
!  We obtain: 
!  * Distances between current integration point and all atoms
!  * Normalized direction to all atoms
!  * radial inverse of distance (for spline evaluation)
!  * log. grid inverse of distance (for spline evaluation)
!
!
!  USES

  use dimensions
  use geometry
  use grids
  use pbc_lists
    implicit none

!  ARGUMENTS

  integer :: current_center
  real*8 :: dist_tab_sq
  real*8, dimension(3) :: dir_tab
  real*8 dist_tab_in
  real*8 i_r
  real*8 i_r_log
  real*8 dir_tab_in ( 3 )


!  INPUTS
!   o current_center -- the atom center
!   o dist_tab_sq -- ( distance to the atom center)**2 
!   o dir_tab -- direction to the atom center
!    
!  OUTPUT
!   o dist_tab_in --  distance to the atom center
!   o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!   o i_r_log -- ???????????????
!   o dir_tab_in -- direnction to the atom center (normalized)
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

	





!  local variables

!     counters

!  begin work

!       tabulate dist_tab, dir_tab first

  dist_tab_in = sqrt(dist_tab_sq)
  
  dir_tab_in(:) = dir_tab(:) / dist_tab_in
  
!       tabulate logarithmic inverse distance next

  i_r_log = &
       invert_log_grid_p2 &
       ( dist_tab_in, &
         species_center(current_center) &
       )

  i_r =  &
       invert_radial_grid &
       ( dist_tab_in, &
       n_radial(species_center(current_center)), &
       scale_radial(species_center(current_center)))

!  that's it

end subroutine tab_single_atom_centered_coords_radial_log_p0
!----------------------------------------------------------------------
!******
