!****s*  FHI-aims/tab_atom_centered_coords_p1
!  NAME
!   tab_atom_centered_coords_p1
!  SYNOPSIS

      subroutine tab_atom_centered_coords_p1 &
      ( coord_current, &
        dist_tab, i_r, dir_tab, &
        n_atoms_list, atoms_list, pos_atoms_list &
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
      use pbc_lists
      use grids
      implicit none

!  ARGUMENTS

      integer :: n_atoms_list
      integer :: atoms_list(n_atoms_list)
      real*8 :: pos_atoms_list(3,n_atoms_list)
      real*8 :: coord_current(3)
      real*8 :: dist_tab ( n_atoms_list, 2 )
      real*8 :: i_r ( n_atoms_list, 2 )
      real*8 :: dir_tab ( 3, n_atoms_list,2 )


!  INPUTS
!   o coord_current -- coordinates of current point
!   o n_atoms_list -- numbers of atoms w.r.t. whom the relative
!                   coordiates of current point need to be determined
!   o atoms_list -- list of atoms
!   o pos_atoms_list -- positions (in x,y,z coordinates) of the list
!                   of atoms.
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

      real*8 :: coords_atoms_opposite(3,n_atoms_list)
!     counters

      integer :: i_atom
      integer :: i_atom_1
      integer :: i_coord
      integer :: i_cell

!  begin work

! determine the position of the atomic positions in the opposite unit
! cell. These is used as a trick for the numerical integration
      coords_atoms_opposite(:,1) = pos_atoms_list(:,1)
      do i_atom = 1, n_atoms, 1
        coords_atoms_opposite(:,i_atom+1) = 2*coords(:,i_atom) - &
            pos_atoms_list(:,i_atom+1)
      enddo 

      dist_tab (:,:) = 0.d0
      do i_atom = 1, n_atoms_list, 1

!       tabulate dist_tab, dir_tab first


        do i_coord = 1, 3, 1

          dir_tab(i_coord, i_atom, 1) = coord_current(i_coord) - &
            pos_atoms_list(i_coord,i_atom)

          dist_tab(i_atom, 1) = dist_tab(i_atom,1) + &
            dir_tab(i_coord, i_atom, 1)**2.0d0

        enddo

        do i_coord = 1, 3, 1

          dir_tab(i_coord, i_atom, 2) = coord_current(i_coord) - &
            coords_atoms_opposite(i_coord,i_atom)

          dist_tab(i_atom, 2) = dist_tab(i_atom, 2) + &
            dir_tab(i_coord, i_atom, 2)**2.0d0

        enddo

      enddo

      do i_cell = 1, 2, 1
       do i_atom = 1, n_atoms_list, 1

        dist_tab(i_atom, i_cell) = sqrt(dist_tab(i_atom,i_cell))

! SVL added this condition 
        i_atom_1=mod(atoms_list(i_atom),n_atoms)
        if(i_atom_1.eq.0) then
         i_atom_1=n_atoms
        endif
        if(dist_tab(i_atom, i_cell).lt.r_grid_min(species(i_atom_1)))then
           dist_tab(i_atom, i_cell) = r_grid_min(species(i_atom_1))
        endif


        if (dist_tab(i_atom, i_cell).gt.1.d-12) then

          do i_coord = 1,3,1
            dir_tab(i_coord, i_atom, i_cell) = &
              dir_tab(i_coord,i_atom, i_cell)/dist_tab(i_atom, i_cell)
          enddo

        else
!         set default direction along z axis; angle should not matter if
!         integration point is directly on top of atom i_atom

          dir_tab(1,i_atom,i_cell) = 0.d0
          dir_tab(2,i_atom,i_cell) = 0.d0
          dir_tab(3,i_atom,i_cell) = 1.d0

        end if


! mapping the atom numbering to the first unit cell
!        i_atom_1=mod(atoms_list(i_atom),n_atoms)
!        if(i_atom_1.eq.0) then
!         i_atom_1=n_atoms
!        endif

!       tabulate logarithmic inverse distance next
        if (dist_tab(i_atom, i_cell).gt.1.d-12) then

          i_r(i_atom, i_cell) = &
          invert_log_grid &
          ( dist_tab(i_atom, i_cell), &
            r_grid_min(species(i_atom_1)), &
            r_grid_inc(species(i_atom_1)) &
          )

        else
          ! safety net using a wrong value
          i_r(i_atom, i_cell) = 0.d0

        end if
       enddo
      enddo

!  that's it

      end subroutine tab_atom_centered_coords_p1
!----------------------------------------------------------------------
!******
