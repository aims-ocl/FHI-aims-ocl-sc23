!****s* FHI-aims/tab_local_geometry_p0
!  NAME
!   tab_local_geometry_p0
!  SYNOPSIS

subroutine tab_local_geometry_p0 &
     ( dist_tab_sq, n_compute_atoms, atom_index, &
     dir_tab, dist_tab, i_r &
     )

!  PURPOSE
!  The subroutine calculates distances to relevant atoms in units of the logarithmic grid, i(r)
!
!  USES

  use dimensions
  use pbc_lists
  use grids
  implicit none

!  ARGUMENTS

      integer :: n_compute_atoms
      real*8, dimension(n_compute_atoms) :: dist_tab_sq
      integer :: atom_index(n_compute_atoms)
      real*8, dimension(3, n_compute_atoms) :: dir_tab
      real*8 dist_tab ( n_compute_atoms )
      real*8 i_r( n_compute_atoms )

!  INPUTS
!   o dist_tab_sq -- not used anymore
!   o n_compute_atoms -- number of relevant atoms
!   o atom_index -- list of relevant atoms
!   o dir_tab -- not used anymore
!   o dir_tab -- not used anymore
!   o dist_tab -- distance to relevant atoms
! 
!  OUTPUT
!   o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
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

      integer :: i_compute_atom
      integer :: i_coord

!  begin work

      do i_compute_atom = 1, n_compute_atoms, 1

         i_r(i_compute_atom) = &
          invert_log_grid &
          ( dist_tab(i_compute_atom), &
            r_grid_min(species_center(atom_index(i_compute_atom))), &
            r_grid_inc(species_center(atom_index(i_compute_atom))) &
          )

      enddo

!  that's it

    end subroutine tab_local_geometry_p0
!----------------------------------------------------------------------
!******
