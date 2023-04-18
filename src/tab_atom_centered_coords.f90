!****s*  FHI-aims/tab_atom_centered_coords
!  NAME
!   tab_atom_centered_coords
!  SYNOPSIS

      subroutine tab_atom_centered_coords &
      ( coord_current, &
        dist_tab, i_r, dir_tab &
      )

!  PURPOSE
!  Subroutine tab_atom_centered_coords tabulates the current integration
!  point as it appears, in atom-centered spherical coordinates, from all
!  atoms in the structure.
!
!  We obtain: 
!  * Distances between current integration point and all atoms
!  * logarithmic inverse of all distances (for spline evaluation)
!  * unit vectors between each atom and current point
!
!  USES

      use dimensions
      use geometry
      use grids
      implicit none

!  ARGUMENTS

      real*8 coord_current(3)
      real*8 dist_tab ( n_atoms )
      real*8 i_r ( n_atoms )
      real*8 dir_tab ( 3, n_atoms )


!  INPUTS
!   o coord_current -- coordinates of current point
!   
!  OUTPUT
!   o dist_tab -- distance to atoms
!   o i_r -- the distance from current integration point to all atoms
!         in units of the logarithmic grid, i(r)
!   o dir_tab -- direction to atoms (normalized). 
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




!  local variables

!     counters

      integer :: i_atom
      integer :: i_coord

!  begin work

      do i_atom = 1, n_atoms, 1

!       tabulate dist_tab, dir_tab first

        dist_tab (i_atom) = 0.d0

        do i_coord = 1,3,1

          dir_tab(i_coord, i_atom) = coord_current(i_coord) - &
            coords(i_coord,i_atom)

          dist_tab(i_atom) = dist_tab(i_atom) + &
            dir_tab(i_coord, i_atom)**2.0d0

        enddo

        dist_tab(i_atom) = sqrt(dist_tab(i_atom))

        if (dist_tab(i_atom).gt.1.d-12) then

          do i_coord = 1,3,1
            dir_tab(i_coord, i_atom) = &
              dir_tab(i_coord,i_atom)/dist_tab(i_atom)
          enddo

        else
!         set default direction along z axis; angle should not matter if
!         integration point is directly on top of atom i_atom

          dir_tab(1,i_atom) = 0.d0
          dir_tab(2,i_atom) = 0.d0
          dir_tab(3,i_atom) = 1.d0

        end if

!       tabulate logarithmic inverse distance next

        if (dist_tab(i_atom).gt.1.d-12) then

          i_r(i_atom) = &
          invert_log_grid &
          ( dist_tab(i_atom), &
            r_grid_min(species(i_atom)), &
            r_grid_inc(species(i_atom)) &
          )

        else
          ! safety net using a wrong value
          i_r(i_atom) = 0.d0

        end if

      enddo

!  that's it

      end subroutine tab_atom_centered_coords
!----------------------------------------------------------------------
!******
