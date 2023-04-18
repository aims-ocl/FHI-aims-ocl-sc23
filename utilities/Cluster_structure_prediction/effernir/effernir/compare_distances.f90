!  subroutine to check whether two geometries are equivalent
!  check energy and distances according to
!  
!  R. Gehrke (2005)
!

subroutine compare_distances(geometry_a, geometry_b, tolerance, equal)

  use control
  use cluster
  use geometry_class

  implicit none

  ! imported variables

  ! input
  type (geometry), intent(in) :: geometry_a
  type (geometry), intent(in) :: geometry_b
  real*8, intent(in) :: tolerance

  ! output
  logical, intent(out) :: equal

  ! local variables
  real*8 :: diff
  real*8 :: numerator
  real*8 :: denominator

  ! counter
  integer :: i_distance
  if (geometry_a%n_distances .ne. geometry_b%n_distances) then
     write (6,*) "Internal inconsistency in compare_distances()."
     write (6,*) "* Aborting."
     stop
  end if

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

  if (diff .le. tolerance) then
     equal = .true.
  else
     equal = .false.
  end if

end subroutine compare_distances
