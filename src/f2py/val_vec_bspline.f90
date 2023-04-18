!****s* FHI-aims/val_vec_bspline
!  NAME
!    val_vec_bspline
!  SYNOPSIS

subroutine val_vec_bspline(n, n_func, n_r, a, b, vals, r, spl_param, &
&                          deriv, is_periodic, if_outside)

  !  PURPOSE
  !
  !    Evaluate an array B-splines for an array of values.
  !
  !  USES

  use bspline
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n, n_func, n_r
  real*8, intent(IN) :: a, b
  real*8, intent(OUT) :: vals(n_func, n_r)
  real*8, intent(IN) :: r(n_r)
  real*8, intent(IN) :: spl_param(0:n+1, n_func)
  integer, intent(IN), optional :: deriv
  logical, intent(IN), optional :: is_periodic
  real*8, intent(IN), optional :: if_outside

  !  INPUTS
  !    o n -- number of interpolation points
  !    o n_func -- number of splines
  !    o n_r -- number of points to evaluate at
  !    o a, b -- definition domain
  !    o r -- positions to interpolate at
  !    o spl_param -- B-spline coefficients
  !    o deriv -- if present, return deriv-th derivative
  !    o is_periodic -- if given and .true., allow sampling within [n,n+1].
  !    o if_outside -- return this if outside of domain
  !  OUTPUTS
  !    o vals -- Values at those points
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

  integer :: i_r
  logical :: my_periodic, is_outside
  character*150 :: info_str
  character(*), parameter :: func = 'val_vec_bspline'

  my_periodic = .false.; if (present(is_periodic)) my_periodic = is_periodic
  if (my_periodic) then
     do i_r = 1, n_r
        call val_periodic_scaled_bspline(n, n_func, a, b, vals(:, i_r), &
        &                                r(i_r), spl_param, deriv)
     end do
  else
     do i_r = 1, n_r
        is_outside = (r(i_r) < a - 1d-11 .or. r(i_r) > b + 1d-11)
        if (is_outside) then
           if (present(if_outside)) then
              vals(:, i_r) = if_outside
           else
              write(info_str, "(A,F12.7,A,F12.7,A,F12.7,A)") &
              & 'Try to evaluate B-Spline at', r(i_r), &
              & ', i.e. outside scope [', a, ',', b, '].'
              call aims_stop(info_str, func)
           end if
        else
           call val_scaled_bspline(n, n_func, a, b, vals(:, i_r), r(i_r), &
           &                       spl_param, deriv)
        end if
     end do
  end if
     
end subroutine val_vec_bspline
!******
