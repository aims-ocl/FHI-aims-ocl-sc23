!  subroutine to check whether two geometries are equivalent
!  check energy and distances according to
!  
!  R. Gehrke (2005)
!

subroutine compare_angles(geometry_a, geometry_b, tolerance, equal)

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
  integer :: i_angle
  ! write (6,*) geometry_a%n_angle, geometry_b%n_angle
  if (geometry_a%n_angle .ne. geometry_b%n_angle) then
     write (6,*) geometry_a%n_angle, geometry_b%n_angle
     write (6,*) "Internal inconsistency in compare_angles()."
     write (6,*) "* Aborting."
     stop
  end if

  diff = 0.d0
  select case (diff_norm) 
     
  case ('average')
     numerator   = 0.d0
     denominator = 0.d0
     do i_angle = 1, geometry_a%n_angle, 1
        numerator   = numerator + (geometry_a%angle(i_angle) - geometry_b%angle(i_angle)) * &
             (geometry_a%angle(i_angle) - geometry_b%angle(i_angle))
        denominator = denominator + geometry_a%angle(i_angle) * geometry_a%angle(i_angle) + &
             geometry_b%angle(i_angle) * geometry_b%angle(i_angle)
     end do
     diff = numerator / denominator
     
  case ('maximum')
     
     do i_angle = 1, geometry_a%n_angle, 1
        numerator   = (geometry_a%angle(i_angle) - geometry_b%angle(i_angle)) * &
             (geometry_a%angle(i_angle) - geometry_b%angle(i_angle))
        denominator = geometry_a%angle(i_angle) * geometry_a%angle(i_angle) + &
             geometry_b%angle(i_angle) * geometry_b%angle(i_angle)
        if ((numerator / denominator) .gt. diff) then
           diff = numerator / denominator
        end if
     end do
     
  end select
  
  if (diff .le. tolerance) then
     equal = .true.
  else
     equal = .false.
  end if

end subroutine compare_angles
