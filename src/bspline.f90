!****h* FHI-aims/bspline
!  NAME
!    bspline
!  SYNOPSIS

module bspline

  !  PURPOSE
  ! 
  !  USES

  implicit none

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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  !
  !    With i = int(x) and t = x - i it follows:
  !
  !                                    (-1  3 -3  1)   ( p_{i-1} )
  ! S(t) = S_i(t) = (t^3 t^2 t^1  1) * ( 3 -6  3  0) * ( p_{i}   )
  !                                    (-3  0  3  0)   ( p_{i+1} )
  !                                    ( 1  4  1  0)   ( p_{i+2} )
  !
  !    The spline is defined from 1 to n.  Therefore, B-Splines from 0 to n+1
  !    are needed.
  
  real*8, parameter :: bspline_blending(4, 4) = reshape( &
  & (/ -1.d0,  3.d0, -3.d0,  1.d0, &
  &     3.d0, -6.d0,  0.d0,  4.d0, &
  &    -3.d0,  3.d0,  3.d0,  1.d0, &
  &     1.d0,  0.d0,  0.d0,  0.d0  /), (/4, 4/)) / 6.d0

contains

  !----------------------------------------------------------------------------
  !****s* bslpine/cubic_bspline_periodic
  !  NAME
  !    cubic_bspline_periodic
  !  SYNOPSIS

  subroutine cubic_bspline_periodic(n, n_func, ff, spl_param)

    !  PURPOSE
    !
    !    Wrapper to cubic_bspline() to generate a periodic spline (the n-th
    !    point is identified with the 0-th).  The accuracy of these kinds of
    !    splines is O(h^4) for periodic functions.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: ff(n, n_func)
    real*8, intent(OUT) :: spl_param(0:n+1, n_func)

    !  INPUTS
    !    o n -- number of grid points
    !    o n_func -- number of splines
    !    o ff -- values at grid points
    !  OUTPUTS
    !    o spl_param -- B-spline parameters
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable :: DL(:), D(:), DU(:)
    real*8, allocatable :: Amux(:,:)
    real*8 :: u_1, u_n, uAmu, uAmx(n_func)
    integer :: i_func, info
    character*150 :: info_str
    character(*), parameter :: func = 'cubic_bspline_periodic'

    !  ( 4 1     1 )        ( 3 1       )     ( 1 )
    !  ( 1 4 1     )        ( 1 4 1     )     ( 0 )
    !  (   1 4 1   )   ==   (   1 4 1   )  +  ( 0 ) * ( 1 0 0 0 1 )
    !  (     1 4 1 )        (     1 4 1 )     ( 0 )
    !  ( 1     1 4 )        (       1 3 )     ( 1 )

    ! Sherman-Morrison formula:
    !   (A + uu^t)^{-1} x = A^-1 x +  A^-1 u u^t A^-1 x / (1 + u^t A^-1 u)
    !                     = Amx - Amu * (u, Amx) / (1 + (u, Amu))
    ! => need A^-1 u, A^-1 x, and A u

    allocate(DL(n-1), D(n), DU(n-1), stat=info)
    call check_allocation(info, 'DL, D, DU', func)
    allocate(Amux(n, 0:n_func), stat=info)
    call check_allocation(info, 'Amux', func)

    ! Get u
    u_1 = sqrt(1.d0 / 6.d0)
    u_n = sqrt(1.d0 / 6.d0)

    ! Middle equations correspond to
    !   S_{i-1}(j) p_{i-1} + S_i(j) p_i + S_{i+1}(j) p_{i+1} == ff(j).
    D(2:n-1)  = 4.d0 / 6.d0    ! B-spline value at its own knot
    D(1) = 3.d0 / 6.d0         ! Diagonal value minus rank-1 correction
    D(n) = 3.d0 / 6.d0         ! Diagonal value minus rank-1 correction
    DL = 1.d0 / 6.d0           ! B-spline value at neighboring knot
    DU = 1.d0 / 6.d0           ! B-spline value at neighboring knot
    Amux(:, 0) = 0.d0
    Amux(1, 0) = u_1
    Amux(n, 0) = u_n
    Amux(:, 1:n_func) = ff
    
    ! Get A^-1 u and A^-1 x
    call dgtsv(n, 1+n_func, DL, D, DU, Amux, n, info)
    if (info /= 0) then
       write(info_str, "('dgtsv returned info =',I6)") info
       call aims_stop(info_str, func)
    end if

    ! Do not need other entries of Au
    uAmu = u_1 * Amux(1, 0) + u_n * Amux(n, 0)
    uAmx = u_1 * Amux(1, 1:n_func) + u_n * Amux(n, 1:n_func)

    do i_func = 1, n_func
       spl_param(1:n, i_func) = Amux(:, i_func) &
       &                      - Amux(:, 0) * uAmx(i_func) / (1.d0 + uAmu)
    end do
    
    ! Periodicity
    spl_param(0,   :) = spl_param(n, :)
    spl_param(n+1, :) = spl_param(1, :)

    deallocate(DL, D, DU, Amux)

  end subroutine cubic_bspline_periodic
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/cubic_bspline_fundamental
  !  NAME
  !    cubic_bspline_fundamental
  !  SYNOPSIS

  subroutine cubic_bspline_fundamental(n, n_func, ff, spl_param, dfL, dfR)

    !  PURPOSE
    !
    !    Wrapper to cubic_bspline() to generate a fundamental spline (the
    !    first derivatives at the borders are set explicitly).  The accuracy
    !    of these kinds of splines is O(h^4) if the derivatives given are
    !    correct.  Besides the periodic splines, these are the best.  But one
    !    needs to know the derivatives.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: ff(n, n_func)
    real*8, intent(OUT) :: spl_param(0:n+1, n_func)
    real*8, intent(IN) :: dfL(n_func), dfR(n_func)

    !  INPUTS
    !    o n -- number of grid points
    !    o n_func -- number of splines
    !    o ff -- values at grid points
    !    o dfL, dfR -- first derivatives at left and right borders
    !  OUTPUTS
    !    o spl_param -- B-spline parameters
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, parameter :: nRL = 3
    real*8, parameter :: cc(nRL) = (/-0.5d0, 0.d0, 0.5d0/)
    character(*), parameter :: func = 'cubic_bspline_natural'

    call cubic_bspline(n, n_func, ff, spl_param, nRL, nRL, cc, cc, dfL, dfR)

  end subroutine cubic_bspline_fundamental
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/cubic_bspline_notaknot
  !  NAME
  !    cubic_bspline_notaknot
  !  SYNOPSIS

  subroutine cubic_bspline_notaknot(n, n_func, ff, spl_param)

    !  PURPOSE
    !
    !    Wrapper to cubic_bspline() to generate a not-a-knot spline (the
    !    second and second-last knots have a continuous third derivative).
    !    The accuracy of these kinds of splines is O((2h)^4) [in practice
    !    reaches O(h^4) off the boundaries].  This is probably the best
    !    general purpose spline.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n
    integer, intent(IN) :: n_func
    real*8, intent(IN) :: ff(n, n_func)
    real*8, intent(OUT) :: spl_param(0:n+1, n_func)

    !  INPUTS
    !    o n -- number of grid points
    !    o n_func -- number of splines
    !    o ff -- values at grid points
    !  OUTPUTS
    !    o spl_param -- B-spline parameters
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, parameter :: nRL = 5
    real*8, parameter :: cc(nRL) = (/-1.d0, 4.d0, -6.d0, 4.d0, -1.d0/)
    real*8 :: zeros(n_func)
    character(*), parameter :: func = 'cubic_bspline_natural'

    zeros = 0.d0
    call cubic_bspline(n, n_func, ff, spl_param, nRL, nRL, cc, cc, zeros, zeros)

  end subroutine cubic_bspline_notaknot
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/cubic_bspline_natural
  !  NAME
  !    cubic_bspline_natural
  !  SYNOPSIS

  subroutine cubic_bspline_natural(n, n_func, ff, spl_param)

    !  PURPOSE
    !
    !    Wrapper to cubic_bspline() to generate a natural spline (second
    !    derivative set to zero at the borders).  The accuracy of these kinds
    !    of splines is O(h^2) [in practice reaches O(h^4) off the boundaries].
    !    While this spline matches the usual definition and the original
    !    motivation (smooth fit with minimized curvature), it performs rather
    !    bad at curved boundaries.  This spline is just another (equivalent)
    !    representation of the spline used in spline.f90.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: ff(n, n_func)
    real*8, intent(OUT) :: spl_param(0:n+1, n_func)

    !  INPUTS
    !    o n -- number of grid points
    !    o n_func -- number of splines
    !    o ff -- values at grid points
    !  OUTPUTS
    !    o spl_param -- B-spline parameters
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer, parameter :: nRL = 3
    real*8, parameter :: cc(nRL) = (/1.d0, -2.d0, 1.d0/)
    real*8 :: zeros(n_func)
    character(*), parameter :: func = 'cubic_bspline_natural'

    zeros = 0.d0
    call cubic_bspline(n, n_func, ff, spl_param, nRL, nRL, cc, cc, zeros, zeros)

  end subroutine cubic_bspline_natural
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/cubic_bspline
  !  NAME
  !    cubic_bspline
  !  SYNOPSIS

  subroutine cubic_bspline(n, n_func, ff, spl_param, nL, nR, cL, cR, fL, fR)

    !  PURPOSE
    !
    !     Calculate cubic B-spline coefficients interpolating ff and
    !     fulfilling generally formulated boundary conditions.
    !
    !     As an n-knot grid has n+2 B-splines, there are two additional
    !     degrees of freedoms.  These are fixed by the equations:
    !         \sum_{i=0}^{nL-1}     cL(i) * p(i) == fL
    !         \sum_{i=n-nR+2}^{n+1} cR(i) * p(i) == fR
    !
    !     This procedure is a helper function for the other interpolating
    !     procedures and is not meant to be called from outside this module
    !     (though it should be safe [but hard] to do so).
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: ff(n, n_func)
    real*8, intent(OUT) :: spl_param(0:n+1, n_func)
    integer, intent(IN) :: nL, nR
    real*8, intent(IN) :: cL(0:nL-1), cR(n+1-nR+1:n+1)
    real*8, intent(IN) :: fL(n_func), fR(n_func)

    !  INPUTS
    !    o n -- Number of knots
    !    o ff -- Function values at knots
    !  OUTPUTS
    !    o spl_param -- B-spline coefficients
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable :: mycL(:), mycR(:)
    real*8 :: myfL(n_func), myfR(n_func), fac
    real*8, allocatable :: D(:), DL(:), DU(:)
    integer :: i, info
    real*8, parameter :: stance(3) = (/1.d0, 4.d0, 1.d0/) / 6.d0
    character*150 :: info_str
    character(*), parameter :: func = 'cubic_bspline'

    allocate(mycL(0:nL-1), mycR(n+1-nR+1:n+1), stat=info)
    call check_allocation(info, 'mycL, mycR', func)
    allocate(DL(0:n), D(0:n+1), DU(1:n+1), stat=info)
    call check_allocation(info, 'DL, D, DU', func)

    ! Middle equations correspond to
    !   S_{i-1}(j) p_{i-1} + S_i(j) p_i + S_{i+1}(j) p_{i+1} == ff(j).
    D  = 4.d0 / 6.d0    ! B-spline value at its own knot
    DL = 1.d0 / 6.d0    ! B-spline value at neighboring knot
    DU = 1.d0 / 6.d0    ! B-spline value at neighboring knot
    spl_param(1:n, 1:n_func) = ff

    mycL = cL
    myfL = fL
    do i = nL-2, 1, -1
       ! Eliminate coefficient i+1
       fac = - mycL(i+1) / stance(3)
       mycL(i-1:i+1) = mycL(i-1:i+1) + fac * stance(1:3)
       myfL = myfL + fac * ff(i,1:n_func)
    end do
    D(0) = mycL(0)
    DU(1) = mycL(1)
    spl_param(0, 1:n_func) = myfL(1:n_func)

    mycR = cR
    myfR = fR
    do i = n+1-nR+2, n
       ! Eliminate coefficient i-1
       fac = - mycR(i-1) / stance(1)
       mycR(i-1:i+1) = mycR(i-1:i+1) + fac * stance(1:3)
       myfR = myfR + fac * ff(i, 1:n_func)
    end do
    D(n+1) = mycR(n+1)
    DL(n) = mycR(n)
    spl_param(n+1, 1:n_func) = myfR(1:n_func)

    call dgtsv(n+2, n_func, DL, D, DU, spl_param, n+2, info)
    if (info /= 0) then
       write(info_str, "('dgtsv returned info =',I6)") info
       call aims_stop(info_str, func)
    end if

    deallocate(DL, D, DU, mycL, mycR)

  end subroutine cubic_bspline
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/val_bspline
  !  NAME
  !    val_bspline
  !  SYNOPSIS

  subroutine val_bspline(n, n_func, val, r, spl_param, deriv, is_periodic)

    !  PURPOSE
    !     Evaluate B-spline at r.
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: r
    real*8, intent(OUT) :: val(n_func)
    real*8, intent(IN) :: spl_param(0:n+1, n_func)
    integer, intent(IN), optional :: deriv
    logical, intent(IN), optional :: is_periodic

    !  INPUTS
    !    o n -- number of interpolation points
    !    o n_func -- number of splines
    !    o r -- position to interpolate at
    !    o spl_param -- B-spline coefficients
    !    o deriv -- if present, return deriv-th derivative
    !    o is_periodic -- if given and .true., allow sampling within [n,n+1].
    !  OUTPUTS
    !    o val -- value at r
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'val_bspline'
    integer :: my_deriv
    logical :: my_periodic
    integer i
    real*8 t, tt(4), blender(4), bb(4, n_func)

    my_deriv = 0; if (present(deriv)) my_deriv = deriv
    my_periodic = .false.; if (present(is_periodic)) my_periodic = is_periodic

    i = int(r)
    if (my_periodic) then
       t = r - dble(i)    ! modulo operation should not change t
       i = modulo(i-1, n) + 1
       if (i == n) then
          bb(1:2,:) = spl_param(n-1:n,:)
          bb(3:4,:) = spl_param(1:2,  :)  ! spl_param: (0) == (n), (1) == (n+1).
       else
          bb = spl_param(i-1:i+2, :)
       end if
    else
       ! non-periodic case
       if (r < 1d0 - 1d-10 .or. r > n + 1d-10) then
          call aims_stop('Try to evaluate B-Spline outside scope', func)
       end if
       if (i < 1) i = 1
       if (i >= n) i = n-1
       bb = spl_param(i-1:i+2, :)
       t = r - dble(i)    ! shift operation should change t
    end if

    ! JW: Could use the de-Boor algorithm for efficiency and stability, but I
    ! am not convinced that we would really benefit for this equidistant grid.

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
    blender = matmul(tt, bspline_blending)
    val = matmul(blender, bb)

  end subroutine val_bspline
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/val_scaled_bspline
  !  NAME
  !    val_scaled_bspline
  !  SYNOPSIS

  subroutine val_scaled_bspline(n, n_func, a, b, val, x, spl_param, deriv)

    !  PURPOSE
    !
    !    Evaluate a scaled ([1,n] -> [a,b]) B-spline.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: a, b
    real*8, intent(OUT) :: val(n_func)
    real*8, intent(IN) :: x
    real*8, intent(IN) :: spl_param(0:n+1, n_func)
    integer, intent(IN), optional :: deriv

    !  INPUTS
    !    o n -- number of interpolation points
    !    o n_func -- number of splines
    !    o [a,b] -- domain of definition
    !    o x -- position to interpolate at (within [a,b])
    !    o spl_param -- B-spline coefficients
    !    o deriv -- if present, return deriv-th derivative
    !  OUTPUTS
    !    o val -- Values at r.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: my_deriv
    real*8 :: dr_dx, r
    character(*), parameter :: func = 'val_scaled_bspline'

    my_deriv = 0; if (present(deriv)) my_deriv = deriv
    dr_dx = (n-1) / (b - a)
    r = 1.d0 + dr_dx * (x - a)
    call val_bspline(n, n_func, val, r, spl_param, my_deriv, .false.)
    val = val * dr_dx**my_deriv

  end subroutine val_scaled_bspline
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/val_periodic_scaled_bspline
  !  NAME
  !    val_periodic_scaled_bspline
  !  SYNOPSIS

  subroutine val_periodic_scaled_bspline(n, n_func, a, b, val, x, &
  &                                      spl_param, deriv)

    !  PURPOSE
    !
    !    Evaluate a scaled ([1,n+1] -> [a,b]) B-spline.  Please note that
    !    here, in contrast to val_scaled_bspline(), the spline is defined up
    !    to n+1, which is equivalent to 1.  Therefore, a and b are equivalent
    !    points.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    real*8, intent(IN) :: a, b
    real*8, intent(OUT) :: val(n_func)
    real*8, intent(IN) :: x
    real*8, intent(IN) :: spl_param(0:n+1, n_func)
    integer, intent(IN), optional :: deriv

    !  INPUTS
    !    o n -- number of interpolation points
    !    o n_func -- number of splines
    !    o [a,b] -- domain of definition
    !    o x -- position to interpolate at (within [a,b])
    !    o spl_param -- B-spline coefficients
    !    o deriv -- if present, return deriv-th derivative
    !  OUTPUTS
    !    o val -- Values at r.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: dr_dx, r
    character(*), parameter :: func = 'val_periodic_scaled_bspline'

    dr_dx = n / (b - a)
    r = 1.d0 + dr_dx * (x - a)
    call val_bspline(n, n_func, val, r, spl_param, deriv, .true.)
    val = val * dr_dx**deriv

  end subroutine val_periodic_scaled_bspline
  !******
  !----------------------------------------------------------------------------
  !****s* bspline/get_subspline
  !  NAME
  !    get_subspline
  !  SYNOPSIS

  subroutine get_subspline(n, n_func, i1, iL, dropfac, spl_in, n_out, &
  &                         spl_out, max_drop_err)

    !  PURPOSE
    !
    !    Return a sparser spline interpolation.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n, n_func
    integer, intent(IN) :: i1, iL, dropfac
    real*8, intent(IN) :: spl_in(0:n+1, n_func)
    integer, intent(IN) :: n_out
    real*8, intent(OUT) :: spl_out(0:n_out+1, n_func)
    real*8, intent(INOUT) :: max_drop_err
    !f2py real*8, intent(IN,OUT) :: max_drop_err

    !  INPUTS
    !    o n -- Number of knots in input spline
    !    o n_func -- Number of splines
    !    o i1, iL -- First and last (input) knot to take into account
    !    o dropfac -- Only keep each dropfac-th knot. => mod(iL-i1,dropfac)==0
    !    o spl_in -- Input splines
    !    o n_out -- Size of output spline (n_out == (iL-i1)/dropfac + 1)
    !    o max_drop_err -- Old maximum dropping error
    !  OUTPUTS
    !    o spl_out -- Coarsened spline
    !    o max_drop_err -- New maximum dropping error
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: df_i1(n_func), df_iL(n_func), val_out(n_func), val_in(n_func)
    real*8 :: r_j
    real*8, allocatable :: ff_kept(:,:)
    integer :: i, j
    integer :: info
    character(*), parameter :: func = 'get_subspline'

    ! --- Initialize

    ! n_out is inclusive, so n_out-1 is the number of coarse intervals.
    if ((n_out-1)*dropfac /= iL - i1) call aims_stop('Invalid n_out', func)
    if (i1 < 1) call aims_stop('Invalid i1', func)
    if (iL > n) call aims_stop('Invalid iL', func)
    if (i1 >= iL) call aims_stop('i1 >= iL', func)

    allocate(ff_kept(n_func, n_out), stat=info)
    call check_allocation(info, 'ff_kept', func)

    ! --- Values at coarser mesh

    do i = i1, iL, dropfac
       j = (i - i1)/dropfac + 1
       call val_bspline(n, n_func, ff_kept(:,j), dble(i), spl_in)
    end do

    ! --- Border derivatives
    
    call val_bspline(n, n_func, df_i1, dble(i1), spl_in, 1)
    call val_bspline(n, n_func, df_iL, dble(iL), spl_in, 1)

    ! --- Spline to (possibly) coarser mesh

    call cubic_bspline_fundamental(n_out, n_func, ff_kept, spl_out, &
    &                              dropfac*df_i1, dropfac*df_iL)

    ! --- Check errors

    if (dropfac > 1) then
       do i = i1, iL
          r_j = (i - i1)/dble(dropfac) + 1.d0
          call val_bspline(n_out, n_func, val_out, r_j, spl_out)
          call val_bspline(n, n_func, val_in, dble(i), spl_in)
          max_drop_err = max(max_drop_err, maxval(abs(val_out - val_in)))
       end do
    end if
    deallocate(ff_kept)

  end subroutine get_subspline
  !******
end module bspline
!******
