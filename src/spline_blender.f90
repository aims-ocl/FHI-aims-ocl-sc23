!****s* FHI-aims/spline_blender
!  NAME
!    spline_blender
!  SYNOPSIS

subroutine spline_blender(n, r, blender, i1, deriv, use_bspline)

  !  PURPOSE
  !     Get coefficients needed for (B-)spline evaluation.
  !  USES

  use bspline
  use mpi_tasks, only: aims_stop
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n
  real*8, intent(IN) :: r
  real*8, intent(OUT) :: blender(4)
  integer, intent(OUT) :: i1
  integer, intent(IN), optional :: deriv
  logical, intent(IN) :: use_bspline

  !  INPUTS
  !    o n -- number of interpolation points
  !    o r -- position to interpolate at ([1, n])
  !    o deriv -- if present, return deriv-th derivative
  !    o use_bspline -- Is it a B-spline (bspline.f90) or not (spline.f90)?
  !  OUTPUTS
  !    o blender, i1 -- dot_product(blender, spl_param(i1:i1+3)) is value
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

  integer :: my_deriv, i
  real*8 t, tt(4)
  character(*), parameter :: func = 'spline_blender'

  my_deriv = 0; if (present(deriv)) my_deriv = deriv

  i = int(r)
  if (r < 1d0 - 1d-10 .or. r > n + 1d-10) then
     call aims_stop('Try to evaluate B-Spline outside scope', func)
  end if
  if (i < 1) i = 1
  if (i >= n) i = n-1
  t = r - dble(i)    ! shift operation should change t

  if (use_bspline) then
     i1 = i+1          ! (0:n+1)    -> (1:n+2)
  else
     i1 = 4*(i-1) + 1  ! (1:4, 1:n) -> (1:4*n)
  end if

  select case (my_deriv)
  case(0)
     tt(1) = t**3
     tt(2) = t**2
     tt(3) = t
     tt(4) = 1.d0
  case(1)
     tt(1) = 3.d0*t**2
     tt(2) = 2.d0*t
     tt(3) = 1.d0
     tt(4) = 0.d0
  case(2)
     tt(1) = 6.d0*t
     tt(2) = 2.d0
     tt(3) = 0.d0
     tt(4) = 0.d0
  case(3)
     tt(1) = 6.d0
     tt(2) = 0.d0
     tt(3) = 0.d0
     tt(4) = 0.d0
  case default
     call aims_stop('Invalid deriv', func)
  end select
  if (use_bspline) then
     blender = matmul(tt, bspline_blending)
  else
     blender = tt(4:1:-1)
  end if

end subroutine spline_blender
!******
