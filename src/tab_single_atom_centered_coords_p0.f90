!****s* FHI-aims/tab_single_atom_centered_coords_p0
!  NAME
!   tab_single_atom_centered_coords_p0
!  SYNOPSIS

subroutine tab_single_atom_centered_coords_p0( current_center, coord_current, &
     dist_tab_sq, dir_tab )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the current integration
!  point as it appears, in atom-centered spherical coordinates, from one
!  atom center.
!
!  Here we first tabulate only the square of all distances, and absolute distance vectors
!  (not normalized to 1)
!
!  USES

  use dimensions
  use pbc_lists
  implicit none

!  ARGUMENTS

      integer :: current_center
      real*8 coord_current(3)
      real*8 dist_tab_sq
      real*8 dir_tab ( 3 )

!  INPUTS
!    o current_center -- center from where distance is calculated
!    o coord_current -- current coordination point
!  OUTPUT
!    o dist_tab_sq -- (distance from atom center)**2
!    o dir_tab -- direction of atom center
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

      dist_tab_sq = 0.d0
      
      do i_coord = 1,3,1
         dir_tab(i_coord) = coord_current(i_coord) - &
              coords_center(i_coord,current_center)
      end do

      do i_coord = 1,3,1
         dist_tab_sq = dist_tab_sq + dir_tab(i_coord)**2.0d0         
      enddo

!  that's it

    end subroutine tab_single_atom_centered_coords_p0
!----------------------------------------------------------------------
!******
