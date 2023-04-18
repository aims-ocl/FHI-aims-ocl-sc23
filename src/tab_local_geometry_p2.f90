!****s* FHI-aims/tab_local_geometry_p2
!  NAME
!   tab_local_geometry_p2
!  SYNOPSIS

      subroutine tab_local_geometry_p2 &
      ( n_compute_atoms, atom_index, &
        dist_tab, i_r &
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
      integer :: atom_index(n_compute_atoms)
      real*8 dist_tab ( n_compute_atoms )
      real*8 i_r( n_compute_atoms )

!  INPUTS
!   o n_compute_atoms -- number of relevant atoms
!   o atom_index -- list of relevant atoms
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

!  begin work

      do i_compute_atom = 1, n_compute_atoms, 1

        i_r(i_compute_atom) = &
          invert_log_grid_p2 &
          ( dist_tab(i_compute_atom), &
            species_center(atom_index(i_compute_atom)) &
          )

      enddo

!  that's it

    end subroutine tab_local_geometry_p2
!----------------------------------------------------------------------
!******
