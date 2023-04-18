! small subroutine to rotate a point around an arbitrary axis
! by an arbitrary angle
!
! R.Gehrke (2007)
!

subroutine rotate(axis, angle, point)

  use vector_class

  implicit none

  ! imported variables

  ! input
  type(vector), intent(in) :: axis
  real*8, intent(in) :: angle

  ! output
  type(vector), intent(inout) :: point

  ! local variables
  type(vector) :: temp_point
  real*8 :: cos_angle
  real*8 :: sin_angle

  cos_angle = cos(angle)
  sin_angle = sin(angle)

  temp_point%x = 0.d0
  temp_point%y = 0.d0
  temp_point%z = 0.d0
  
  temp_point%x = temp_point%x + (cos_angle + (1.d0 - cos_angle) * axis%x * axis%x) * point%x
  temp_point%x = temp_point%x + ((1.d0 - cos_angle) * axis%x * axis%y - axis%z * sin_angle) * point%y
  temp_point%x = temp_point%x + ((1.d0 - cos_angle) * axis%x * axis%z + axis%y * sin_angle) * point%z

  temp_point%y = temp_point%y + ((1.d0 - cos_angle) * axis%x * axis%y + axis%z * sin_angle) * point%x
  temp_point%y = temp_point%y + (cos_angle + (1.d0 - cos_angle) * axis%y * axis%y) * point%y
  temp_point%y = temp_point%y + ((1.d0 - cos_angle) * axis%y * axis%z - axis%x * sin_angle) * point%z

  temp_point%z = temp_point%z + ((1.d0 - cos_angle) * axis%x * axis%z - axis%y * sin_angle) * point%x
  temp_point%z = temp_point%z + ((1.d0 - cos_angle) * axis%y * axis%z + axis%x * sin_angle) * point%y
  temp_point%z = temp_point%z + (cos_angle + (1.d0 - cos_angle) * axis%z * axis%z) * point%z

  point = temp_point

end subroutine rotate
