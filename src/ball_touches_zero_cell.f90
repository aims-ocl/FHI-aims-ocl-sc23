!****s* FHI-aims/ball_touches_zero_cell
!  NAME
!    ball_touches_zero_cell
!  SYNOPSIS

subroutine ball_touches_zero_cell(center, radius, margin, touches)

  !  PURPOSE
  !
  !    Figure out if the sphere ("ball") around center with radius touches the
  !    0-cell (this comment used to say "real-space Wigner-Seitz cell" but
  !    what is really tested is whether we touch any of the faces of a
  !    unit cell shifted such that the coordinate system's origin is in the center).
  !
  !  USES

  use geometry
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: center(3)
  real*8, intent(IN) :: radius
  real*8, intent(IN) :: margin
  logical, intent(OUT) :: touches

  !  INPUTS
  !    o center -- center of the ball
  !    o radius -- radius of the ball
  !    o margin -- in final test, add this margin to radius
  !  OUTPUTS
  !    o touches -- .true. if it touches the 0-cell
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  integer :: max_points_in_side
  real*8 :: rmax(3), vec_1(3), vec_2(3), vec_3(3), coord(3), point_at_surface(3)
  real*8 :: dist_sq
  real*8 :: max_dist_sq
  integer :: i_lat1, i_lat2, i_lat3, i_r1, i_r2, i_r3
  character(*), parameter :: func = 'ball_touches_zero_cell'

  ! JW: Could do this analytically.  May be will in the future.  But right
  !     now, I do not want any change in behavior whatsoever.

  ! This is the number of points in discretized surface of the lattice
  ! cell.  Number is chosen quite randomly, it worked in few test cases,
  ! but some more testing would not harm. - Paula (20.2.2007)
  max_points_in_side = int(maxval(length_of_lattice_vector))*2

  rmax(1) = sum(abs(lattice_vector(1,:)))*0.5d0
  rmax(2) = sum(abs(lattice_vector(2,:)))*0.5d0
  rmax(3) = sum(abs(lattice_vector(3,:)))*0.5d0

  ! While not addressing the issue analytically (which we should)
  ! at least avoid square roots in the innermost loop. Those
  ! cost time.
  max_dist_sq = (radius + margin)**2

  touches = .false.
  if (all(abs(center(:)) < rmax(:) + radius)) then
     go_over_surfaces:  do i_lat1 = 1,3
        do i_lat2 = 1,i_lat1-1
           do i_lat3 = 1,3
              if(i_lat3 /= i_lat1 .and. i_lat3 /= i_lat2 )then
                 ! Basis for a grid on the parallel epiped surface
                 vec_1 = 0.5d0 * lattice_vector(:, i_lat1) &
                 &       / dble(max_points_in_side)
                 vec_2 = 0.5d0 * lattice_vector(:, i_lat2) &
                 &       / dble(max_points_in_side)
                 ! Vector to the surface
                 vec_3 = 0.5d0 * lattice_vector(:, i_lat3)

                 do i_r1 = -max_points_in_side,max_points_in_side,1 
                    do i_r2 = -max_points_in_side,max_points_in_side,1 
                       do i_r3 = -1, 1, 2 

                          point_at_surface =  i_r1 * vec_1 &
                          &                +  i_r2 * vec_2 &
                          &                +  i_r3 * vec_3
                          coord = point_at_surface - center
                          dist_sq = dot_product(coord, coord)

                          if (dist_sq < max_dist_sq) then
                             touches = .true.
                             exit  go_over_surfaces
                          end if
                       end do
                    end do
                 end do
              end if
           end do
        end do
     end do go_over_surfaces
  end if

end subroutine ball_touches_zero_cell
!******
