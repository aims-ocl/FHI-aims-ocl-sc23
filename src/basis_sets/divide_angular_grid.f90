!****s* FHI-aims/divide_angular_grid
!  NAME
!    divide_angular_grid
!  SYNOPSIS 
      subroutine divide_angular_grid( &
           n_angular, r_angular, w_angular, &
           n_division, division_boundaries, &
           r_radial)
!  PURPOSE
!    A routine to partition one anglar grid given on the unit spehere
!    In this function the arrays r_angular and w_angular are reorganised
!    to conform the division given by n_division = number of divisions of the grid
!    (usually 1, 2, 4, or 8), and division_boundaries = [0, l1, l2, l3,..., ln_division]
!    where li = last index of the division i.
!  USES
      use dimensions,      only : n_max_angular, n_max_angular_division
      use runtime_choices, only : angular_div_cutoff, n_min_points_in_division
      use localorb_io,     only : localorb_info
      implicit none
!  ARGUMENTS
      integer :: n_angular
      real*8 :: r_angular(3,n_max_angular)
      real*8 :: w_angular(n_max_angular)
      integer :: n_division
      integer :: division_boundaries(n_max_angular_division+1)
      real*8 :: r_radial
!  INPUTS
!    o n_angular -- number of angular points in the grid
!    o r_angular -- coordinates of the angular points in the grid
!    o w_angular -- integration weights of the angular points in the grid
!    o r_radial -- radial coordnate of the angular grid
!  OUTPUT
!    o n_division -- number of divisions in the angular grid
!    o division_boundaries -- boudary indeces of the divisions
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


!     local
      integer :: division_index(n_angular)
      integer :: points_in_division(n_max_angular_division)
      real*8 :: r_permuted(3,n_angular)
      real*8 :: w_permuted(n_angular)

      integer :: i_index, i_angular, i_division

      points_in_division = 0

      if (r_radial.le.angular_div_cutoff) then
!     if we are below the cutoff raduis, do nothing
         n_division = 1
         division_index = 1
         points_in_division(1) = n_angular

      else if ((int(n_angular/n_min_points_in_division) &
           .ge.8).and. &
           (n_max_angular_division.ge.8)) then
!     divide into 8 blocks, by x-, y-, and z-coordinates
         n_division = 8
         do i_angular = 1, n_angular, 1
            if ((r_angular(1,i_angular).le.0).and. &
                 (r_angular(2,i_angular).le.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 1
               points_in_division(1) = &
                    points_in_division(1) + 1
            else if((r_angular(1,i_angular).gt.0).and. &
                    (r_angular(2,i_angular).le.0).and. &
                    (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 2
               points_in_division(2) = &
                    points_in_division(2) + 1
            else if((r_angular(1,i_angular).le.0).and. &
                    (r_angular(2,i_angular).gt.0).and. &
                    (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 3
               points_in_division(3) = &
                    points_in_division(3) + 1
            else if((r_angular(1,i_angular).le.0).and. &
                    (r_angular(2,i_angular).le.0).and. &
                    (r_angular(3,i_angular).gt.0)) then
               division_index(i_angular) = 4
               points_in_division(4) = &
                    points_in_division(4) + 1
            else if((r_angular(1,i_angular).gt.0).and. &
                    (r_angular(2,i_angular).gt.0).and. &
                    (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 5
               points_in_division(5) = &
                    points_in_division(5) + 1
            else if((r_angular(1,i_angular).gt.0).and. &
                    (r_angular(2,i_angular).le.0).and. &
                    (r_angular(3,i_angular).gt.0)) then
               division_index(i_angular) = 6
               points_in_division(6) = &
                    points_in_division(6) + 1
            else if((r_angular(1,i_angular).le.0).and. &
                    (r_angular(2,i_angular).gt.0).and. &
                    (r_angular(3,i_angular).gt.0)) then
               division_index(i_angular) = 7
               points_in_division(7) = &
                    points_in_division(7) + 1
            else
               division_index(i_angular) = 8
               points_in_division(8) = &
                    points_in_division(8) + 1
            end if
         enddo

      else if ((int(n_angular/n_min_points_in_division) &
              .ge.4).and. &
              n_max_angular_division.ge.4) then
!     divide into 4 blocks, by x- and z-coordinates
         n_division = 4
         do i_angular = 1, n_angular, 1
            if ((r_angular(1,i_angular).le.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 1
               points_in_division(1) = &
                    points_in_division(1) + 1
            else if((r_angular(1,i_angular).gt.0).and. &
                    (r_angular(3,i_angular).le.0)) then
               division_index(i_angular) = 2
               points_in_division(2) = &
                    points_in_division(2) + 1
            else if((r_angular(1,i_angular).le.0).and. &
                    (r_angular(3,i_angular).gt.0)) then
               division_index(i_angular) = 3
               points_in_division(3) = &
                    points_in_division(3) + 1
            else
               division_index(i_angular) = 4
               points_in_division(4) = &
                    points_in_division(4) + 1
            end if
         enddo

      else if ((int(n_angular/n_min_points_in_division) &
              .ge.2).and. &
              n_max_angular_division.ge.2) then
!     divide into 2 blocks, by z-coordinate
         n_division = 2
         do i_angular = 1, n_angular, 1
            if (r_angular(3,i_angular).le.0) then
               division_index(i_angular) = 1
               points_in_division(1) = &
                    points_in_division(1) + 1
            else
               division_index(i_angular) = 2
               points_in_division(2) = &
                    points_in_division(2) + 1
            end if
         enddo

      else if (n_max_angular_division.ge.1) then
!     if everything else fails, don't divide at all
         n_division = 1
         division_index = 1
         points_in_division(1) = n_angular

      else
         call localorb_info('* Error in dividing the grids.')
         stop

      end if

      division_boundaries(1) = 0

      if (n_division.gt.1) then
!       permute the angular grids points so that the first division
!       comes first, second comes second, etc.

      i_index = 0

      do i_division = 1, n_division, 1
         do i_angular = 1, n_angular, 1
            if (division_index(i_angular).eq.i_division) then
               i_index = i_index + 1
               r_permuted(:,i_index) = r_angular(:,i_angular)
               w_permuted(i_index) = w_angular(i_angular)
            end if
         enddo
         division_boundaries(i_division+1) = &
              sum(points_in_division(1:i_division))
      enddo

!     copy the permuted arrays to the original ones
      r_angular(:,1:n_angular) = r_permuted(:,1:n_angular)
      w_angular(1:n_angular) = w_permuted(1:n_angular)

      else
         ! only one division
         division_boundaries(2) = &
              n_angular

      end if

      end subroutine divide_angular_grid
!******





