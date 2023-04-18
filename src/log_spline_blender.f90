!****s* FHI-aims/log_spline_blender
!  NAME
!    log_spline_blender
!  SYNOPSIS

recursive subroutine log_spline_blender(dist, n_r, r_min, r_inc, lmax, &
&                                       use_bwave, want_dblender, &
&                                       extrapolate_to_zero, &
&                                       i1, blender, dblender)

  !  PURPOSE
  !    Prepare evaluation of u(r)/r^{l+1} where u(r) is given as log-spline.
  !  USES

  use grids
  implicit none

  !  ARGUMENTS

  real*8, intent(IN) :: dist
  integer, intent(IN) :: n_r
  real*8, intent(IN) :: r_min, r_inc
  integer, intent(IN) :: lmax
  logical, intent(IN) :: use_bwave, want_dblender
  logical, intent(IN) :: extrapolate_to_zero
  integer, intent(OUT) :: i1
  real*8, intent(OUT) :: blender(4, 0:lmax)
  real*8, intent(OUT) :: dblender(4, 0:lmax)

  !  INPUTS
  !    o dist -- |rvec|
  !    o r_min -- smallest r of loggrid
  !    o r_inc -- scaling factor of loggrid
  !    o lmax -- maximum angular momentum (in factor 1/r^{l+1})
  !    o want_dblender -- need derivative?
  !  OUTPUTS
  !    o blender, i1 -- dot_product(blender(:,l), spl_param(i1:i1+3)) is value
  !    o dblender -- dot_product(dblender(:,l), spl_param(i1:i1+3)) is deriv
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

  real*8 :: ilog  ! real-valued point in logspace \in [1, n_r]
  real*8 :: iprime
  real*8 :: u_blender(4), du_blender(4)
  real*8 :: frmin_blender(4, 0:lmax), dfrmin_blender(4, 0:lmax)
  integer :: i_l
  character(*), parameter :: func = 'log_spline_blender'

  ilog = invert_log_grid(dist, r_min, r_inc)
  if (ilog < 1.d0-1d-11 .or. ilog > n_r+1d-11) then
     if (extrapolate_to_zero) then
        ! Evaluate f(r_min), f'(r_min) [recursion with max-depth 1]
        call log_spline_blender(r_min, n_r, r_min, r_inc, lmax, &
        &                       use_bwave, .true., .false., &
        &                       i1, frmin_blender, dfrmin_blender)
        ! Linearly extrapolate f(r) from f(r_min) to f(dist).
        blender = frmin_blender - (r_min - dist) * dfrmin_blender
        if (want_dblender) dblender = dfrmin_blender
     else
        ! For now, set region within innermost log-shell to zero.
        i1 = 1
        blender = 0.d0
        if (want_dblender) dblender = 0.d0
     end if
  else
     call spline_blender(n_r, ilog, u_blender, i1, 0, use_bwave)
     do i_l = 0, lmax
        blender(:, i_l) = u_blender / dist**(1+i_l)
     end do
     if (want_dblender) then
        ! d(u(i(r))/r)/dr = u'(i(r)) * i'(r) / r^{l+1} - (1+l) u(i(r))/r^{2+l}
        ! BL: Interface has changed 
        ! iprime = invert_log_grid_deriv(dist, r_min, r_inc)
        iprime = invert_log_grid_deriv(dist, r_inc)
        call spline_blender(n_r, ilog, du_blender, i1, 1, use_bwave)
        do i_l = 0, lmax
           dblender(:, i_l) = du_blender * iprime / dist**(1+i_l) &
           &                - (1+i_l) * u_blender / dist**(2+i_l)
        end do
     end if
  end if
end subroutine log_spline_blender
!******
