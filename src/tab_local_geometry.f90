!****s* FHI-aims/tab_local_geometry
!  NAME
!   tab_local_geometry
!  SYNOPSIS

subroutine tab_local_geometry &
     ( dist_tab_sq, n_compute_atoms, atom_index, &
     dir_tab, &
     dist_tab, i_r &
     )

!  PURPOSE
!  The subroutine calculates :
!  o distance to relevant atoms from table of squared distances to all atoms
!  o distances to relevant atoms in units of the logarithmic grid, i(r)
!  o normalizes directions to relevant atoms from table of directions to all atoms.
!
!  USES

  use dimensions
  use geometry
  use grids
  implicit none

!  ARGUMENTS

      real*8, dimension(n_atoms) :: dist_tab_sq
      integer :: n_compute_atoms
      integer :: atom_index(n_atoms)
      real*8, dimension(3, n_atoms) :: dir_tab
      real*8 dist_tab ( n_atoms )
      real*8 i_r ( n_atoms )

!  INPUTS
!   o dist_tab_sq -- (distance to atoms)**2 data to all atoms
!   o n_compute_atoms -- number of relevant atoms
!   o atom_index -- list of relevant atoms
!   o dir_tab -- direction to all atoms (not normalized)
! 
!  OUTPUT
!   o dir_tab -- direction to relevant atoms (normalized)
!   o dist_tab -- distance to relevant atoms
!    o i_r -- the distance from current integration point to all atoms
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

        dist_tab (i_compute_atom) = &
          sqrt(dist_tab_sq(atom_index(i_compute_atom)))

        ! rely here on the fact that atoms are ORDERED,
        ! i.e. i_compute_atom <= atom_index ALWAYS.

        dir_tab(1:3,i_compute_atom) = &
          dir_tab(1:3,atom_index(i_compute_atom)) / &
          dist_tab (i_compute_atom)

        i_r(i_compute_atom) = &
          invert_log_grid &
          ( dist_tab(i_compute_atom), &
            r_grid_min(species(atom_index(i_compute_atom))), &
            r_grid_inc(species(atom_index(i_compute_atom))) &
          )

      enddo

!  that's it

      end subroutine tab_local_geometry
!----------------------------------------------------------------------
!******
