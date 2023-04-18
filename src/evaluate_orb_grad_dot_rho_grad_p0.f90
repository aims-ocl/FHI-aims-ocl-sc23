! evaluates dot product between KS_orbital_gradients and
! xc_gradient_deriv (including gradient rho respectively)
! needed for gga-forces
!
! R.Gehrke (2006)

subroutine evaluate_orb_grad_dot_rho_grad_p0(KS_orbital_gradient, xc_gradient_deriv, n_points, &
     max_occ_number, orb_grad_dot_rho_grad)

  use dimensions
  use runtime_choices

  implicit none

  ! imported variables

  ! input
  real*8, dimension(n_states, n_max_batch_size, 3), intent(in) :: KS_orbital_gradient
  real*8, dimension(3, n_max_batch_size), intent(in) :: xc_gradient_deriv
  integer, intent(in) :: n_points
  integer, intent(in) :: max_occ_number

  ! output
  real*8, dimension(n_states, n_points), intent(out) :: orb_grad_dot_rho_grad

  ! local variables

  ! counter
  integer :: i_state
  integer :: i_coord
  integer :: i_point

  orb_grad_dot_rho_grad = 0.d0
  do i_coord = 1, 3, 1
     do i_point = 1, n_points, 1
        call daxpy(max_occ_number, xc_gradient_deriv(i_coord, i_point), &
             KS_orbital_gradient(1,i_point,i_coord), 1, &
             orb_grad_dot_rho_grad(1,i_point), 1)
     end do
  end do

end subroutine evaluate_orb_grad_dot_rho_grad_p0
