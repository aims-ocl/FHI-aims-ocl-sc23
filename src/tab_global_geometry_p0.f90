!****s* FHI-aims/tab_global_geometry_p0
!  NAME
!   tab_global_geometry_p0
!  SYNOPSIS

subroutine tab_global_geometry_p0 &
     ( dist_tab_sq, dir_tab, dist_tab, i_r, &
     dir_tab_norm,&
     n_atom_list, atom_list &
     )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the current integration
!  point as it appears, in atom-centered spherical coordinates, from all
!  RELEVANT atoms in the structure. 
!
!  USES

  use dimensions
  use pbc_lists
  use grids
  implicit none

!  ARGUMENTS

  integer:: n_atom_list
  integer:: atom_list(n_atom_list)
  real*8, dimension(n_atom_list) :: dist_tab_sq
  real*8, dimension(3,n_atom_list) :: dir_tab
  real*8 dist_tab ( n_atom_list)
  real*8 i_r ( n_atom_list)
  real*8, dimension(3,n_atom_list) :: dir_tab_norm

!  INPUTS
!   o n_atom_list -- number of atoms (including periodic mirror images)
!   o atom_list -- list of atoms
!   o dist_tab_sq -- (distance to atoms)**2
!   o dir_tab -- direction to atoms (not normalized)
!
!  OUTPUT
!   o dist_tab -- distance to atoms
!   o i_r --  the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!   o dir_tab_norm -- direction to atoms (normalized)
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

      integer :: i_center
      integer :: i_center_L
      integer :: i_coord

!  begin work

      do i_center_L = 1, n_atom_list, 1
         i_center = atom_list(i_center_L)
         dist_tab (i_center_L) = sqrt(dist_tab_sq(i_center_L))
         i_r(i_center_L) = &
              invert_log_grid &
              ( dist_tab(i_center_L), &
              r_grid_min(species_center(i_center)), &
              r_grid_inc(species_center(i_center)) &
              )
         dir_tab_norm(1,i_center_L) = dir_tab(1,i_center_L)/dist_tab(i_center_L)
         dir_tab_norm(2,i_center_L) = dir_tab(2,i_center_L)/dist_tab(i_center_L)
         dir_tab_norm(3,i_center_L) = dir_tab(3,i_center_L)/dist_tab(i_center_L)
      enddo

!  that's it

    end subroutine tab_global_geometry_p0
!----------------------------------------------------------------------
!******
subroutine tab_i_r_p0 &
     ( dist_r, i_center, i_r )!, &
!     n_atom_list, atom_list &
!     )

!  PURPOSE
!  gets distance dist_r from point to atom i_center in log grid coords.
!
!  USES

  use dimensions
  use pbc_lists
  use grids
  implicit none

!  ARGUMENTS

!  integer:: n_atom_list
!  integer:: atom_list(n_atom_list)
!  real*8, dimension(n_atom_list) :: dist_tab_sq
!  real*8, dimension(3,n_atom_list) :: dir_tab
!  real*8 dist_tab ( n_atom_list)
  real*8 i_r !( n_atom_list)
!  real*8, dimension(3,n_atom_list) :: dir_tab_norm
  real*8 dist_r

!  INPUTS
!   o n_atom_list -- number of atoms (including periodic mirror images)
!   i_center -- numer of atom of interest

!   o atom_list -- list of atoms
!   o dist_r -- (distance to atom)
!
!  OUTPUT
!   o i_r --  the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!
!  local variables

!     counters

      integer :: i_center
      integer :: i_center_L
      integer :: i_coord

!  begin work

         !i_r(i_center_L) = &
         i_r = &
              invert_log_grid &
              !( dist_tab(i_center_L), &
              ( dist_r, &
              r_grid_min(species_center(i_center)), &
              r_grid_inc(species_center(i_center)) &
              )

!  that's it

    end subroutine tab_i_r_p0
!----------------------------------------------------------------------
!******
