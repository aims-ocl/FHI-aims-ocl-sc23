program test_cp
  use cartesian_ylm
  use constants
  use localorb_io
  implicit none
  real*8 :: rvec(3)
  real*8 :: rotmat1(3,3), rotmat2(3, 3), rotmat(3,3)
  real*8, parameter :: phi = 0.3*pi
  real*8, parameter :: theta = 0.62*pi

  call initialize_cartesian_ylm(4)

  rvec = (/0.1d0, 1.0d0, 2.2d0/)
  rotmat1(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
  rotmat1(2,:) = (/ 0.d0, cos(theta), sin(theta) /)
  rotmat1(3,:) = (/ 0.d0, -sin(theta), cos(theta) /)
  rotmat2(1,:) = (/ cos(phi), sin(phi), 0.d0 /)
  rotmat2(2,:) = (/ -sin(phi), cos(phi), 0.d0 /)
  rotmat2(3,:) = (/ 0.d0, 0.d0, 1.d0 /)
  rotmat = matmul(rotmat1, rotmat2)

  write(0,"(9F10.4)") rotmat

  call test_cartesian_rotation(rotmat, rvec)
  call test_ylm_rotation(rotmat, rvec)
  call localorb_info('Test successfull!')

  call cleanup_cartesian_ylm

end program test_cp
