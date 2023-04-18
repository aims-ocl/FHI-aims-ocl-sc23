!  subroutine to check whether two geometries are equivalent
!  check energy and distances according to
!  
!  R. Gehrke (2005)
!

subroutine diff_distances(geometry_a, geometry_b, diff)

  use control
  use cluster
  use geometry_class

  implicit none

  ! imported variables

  ! input
  type (geometry), intent(in) :: geometry_a
  type (geometry), intent(in) :: geometry_b

  ! output
  real*8 :: diff

  ! local variables
  real*8 :: numerator
  real*8 :: denominator

  ! counter
  integer :: i_distance
  if (geometry_a%n_distances .ne. geometry_b%n_distances) then
     write (6,*) "Internal inconsistency in compare_distances()."
     write (6,*) "* Aborting."
     stop
  end if

!  write (150,*) "geometry_a:"
!  write (150,*) "distances:"
!  do i_distance = 1, geometry_a%n_distances, 1
!     write (150,*) geometry_a%distance(i_distance)
!  end do

!  write (150,*) "geometry_b:"
!  write (150,*) "distances:"
!  do i_distance = 1, geometry_b%n_distances, 1
!     write (150,*) geometry_b%distance(i_distance)
!  end do

  diff = 0.d0

  select case (diff_norm)
  ! "average-norm"
  case ('average')
     numerator   = 0.d0
     denominator = 0.d0
     do i_distance = 1, geometry_a%n_distances, 1
        numerator   = numerator + (geometry_a%distance(i_distance) - geometry_b%distance(i_distance)) * &
             (geometry_a%distance(i_distance) - geometry_b%distance(i_distance))
        denominator = denominator + geometry_a%distance(i_distance) * geometry_a%distance(i_distance) + &
             geometry_b%distance(i_distance) * geometry_b%distance(i_distance)
     end do
     diff = numerator / denominator
  case ('maximum')
     ! "maximum-norm"
     do i_distance = 1, geometry_a%n_distances, 1
        numerator   = (geometry_a%distance(i_distance) - geometry_b%distance(i_distance)) * &
             (geometry_a%distance(i_distance) - geometry_b%distance(i_distance))
        denominator = geometry_a%distance(i_distance) * geometry_a%distance(i_distance) + &
             geometry_b%distance(i_distance) * geometry_b%distance(i_distance)
        if ((numerator / denominator) .gt. diff) then
           diff = numerator / denominator
        end if
     end do
  case default

     write (6,*) "Internal inconsistency in basin-hopping!"
     write (6,*) "* Aborting."
     stop

  end select

end subroutine diff_distances
