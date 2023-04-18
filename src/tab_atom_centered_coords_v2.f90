!****s* FHI-aims/tab_atom_centered_coords_v2
!  NAME
!   tab_atom_centered_coords_v2
!  SYNOPSIS

      subroutine tab_atom_centered_coords_v2 &
      ( coord_current, &
        dist_tab_sq, dir_tab &
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
      use geometry
      implicit none

!  ARGUMENTS

      real*8 coord_current(3)
      real*8 dist_tab_sq ( n_atoms )
      real*8 dir_tab ( 3,n_atoms )

!  INPUTS
!   o coord_current -- coordinates of current point
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

      integer :: i_coord

!  begin work

!       tabulate dist_tab, dir_tab first

        dist_tab_sq (1:n_atoms) = 0.d0

        do i_coord = 1,3,1

          dir_tab(i_coord,1:n_atoms) = coord_current(i_coord) - &
            coords(i_coord,1:n_atoms)

          dist_tab_sq(1:n_atoms) = dist_tab_sq(1:n_atoms) + &
            dir_tab(i_coord, 1:n_atoms)**2.0d0

        enddo


!  that's it

      end subroutine tab_atom_centered_coords_v2
!----------------------------------------------------------------------
!******
