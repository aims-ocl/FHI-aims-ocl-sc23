!****s* FHI-aims/tab_atom_centered_coords_p0
!  NAME
!   tab_atom_centered_coords_p0
!  SYNOPSIS

      subroutine tab_atom_centered_coords_p0 &
      ( coord_current, &
        dist_tab_sq, dir_tab, &
        n_atoms_list, atoms_list &
      )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the current integration
!  point as it appears, in atom-centered spherical coordinates, from all
!  atoms in the structure. 
!
!  Here we first tabulate only the square of all distances, and absolute distance vectors
!  (not normalized to 1)
!
!  USES

      use dimensions
      use pbc_lists
      implicit none

!  ARGUMENTS

      real*8 coord_current(3)
      integer:: n_atoms_list
      integer:: atoms_list(n_atoms_list)
      real*8 dist_tab_sq ( n_atoms_list)
      real*8 dir_tab ( 3,n_atoms_list)

!  INPUTS
!   o coord_current -- coordinates of current point 
!   o n_atoms_list -- number of atoms in this routine
!   o atoms_list -- list of atoms
!
!  OUTPUT
!   o dist_tab_sq -- (distance to atoms)**2
!   o dir_tab -- direction to atoms (not normalized)
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

      integer :: i_coord, i_center_L

!  begin work

        do i_center_L = 1, n_atoms_list,1
           dir_tab(1,i_center_L) = coord_current(1) - &
                coords_center(1,atoms_list(i_center_L))
           dir_tab(2,i_center_L) = coord_current(2) - &
                coords_center(2,atoms_list(i_center_L))
           dir_tab(3,i_center_L) = coord_current(3) - &
                coords_center(3,atoms_list(i_center_L))
           
           dist_tab_sq(i_center_L) = &
                  dir_tab(1,i_center_L)*dir_tab(1,i_center_L) &
                + dir_tab(2,i_center_L)*dir_tab(2,i_center_L) &
                + dir_tab(3,i_center_L)*dir_tab(3,i_center_L) 

        end do

!  that's it

      end subroutine tab_atom_centered_coords_p0
!----------------------------------------------------------------------
!******



!****s* FHI-aims/tab_atom_centered_coords_p0_less
!  NAME
!   tab_atom_centered_coords_p0_less
!  SYNOPSIS

      subroutine tab_atom_centered_coords_p0_less &
      ( coord_current, &
        dist_tab_sq, &
        n_atoms_list, atoms_list, coords_center_combined &
      )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the current integration
!  point as it appears, in atom-centered spherical coordinates, from all
!  atoms in the structure. 
!
!  Here we first tabulate only the square of all distances, and absolute distance vectors
!  (not normalized to 1)
!
!  USES

      use dimensions
      use pbc_lists
      implicit none

!  ARGUMENTS

      real*8 coord_current(3)
      integer:: n_atoms_list
      integer:: atoms_list(n_atoms_list)
      real*8 dist_tab_sq ( n_atoms_list)

      real*8 coords_center_combined(3,n_centers)

      real*8 dir_tab (3)

!  INPUTS
!   o coord_current -- coordinates of current point 
!   o n_atoms_list -- number of atoms in this routine
!   o atoms_list -- list of atoms
!
!  OUTPUT
!   o dist_tab_sq -- (distance to atoms)**2
!   o dir_tab -- direction to atoms (not normalized)
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

      integer :: i_coord, i_center_L

!  begin work

        do i_center_L = 1, n_atoms_list,1
           dir_tab(1) = coord_current(1) - &
                coords_center_combined(1, i_center_L)
           dir_tab(2) = coord_current(2) - &
                coords_center_combined(2, i_center_L)
           dir_tab(3) = coord_current(3) - &
                coords_center_combined(3, i_center_L)
           
           dist_tab_sq(i_center_L) = &
                  dir_tab(1)*dir_tab(1) &
                + dir_tab(2)*dir_tab(2) &
                + dir_tab(3)*dir_tab(3) 

        end do

!  that's it

      end subroutine tab_atom_centered_coords_p0_less
!----------------------------------------------------------------------
!******
