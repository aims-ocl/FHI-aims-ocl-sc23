!****h* FHI-aims/basbas_fn_coulomb
!  NAME
!    basbas_fn_coulomb
!  SYNOPSIS

module basbas_fn_coulomb

  !  PURPOSE
  !
  !    Provide utilities for the Coulomb potential of the radial parts
  !    of the product basis functions.
  !
  !    In principle, these routines could also be applied to ordinary basis
  !    functions; there is just no point.  On the other hand, these routines
  !    should not be used for the radial Hartree potentials because the radial
  !    charge density does not include the r factor in u(r) = rP(r).
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

contains

  !----------------------------------------------------------------------------
  !****s* basbas_coulomb/trapezoidal_wave_integrator
  !  NAME
  !    trapezoidal_wave_integrator
  !  SYNOPSIS

  subroutine trapezoidal_wave_integrator(i_l, wave, n_grid, r_grid, &
  &                                      coulomb, multipole)

    !  PURPOSE
    !
    !    Calculate the Coulomb field
    !
    !         V(r) = 4*pi / (2*i_l + 1) *
    !                [\int_0^r      dr' r'**(2+i_l)   wave(r')/r' / r**(i_l+1) +
    !                 \int_r^\infty dr' r'**(2-i_l-1) wave(r')/r' * r**(i_l)]
    !                
    !    and the multipole moment M of wave.
    !
    !         M := \int_0^\infty dr * wave(r)/r' * r**(2+i_l) / (2*i_l+1)
    !            = \sum_r wave(r) * r**(2+i_l) * log(dr/r) / (2*i_l+1)
    !
    !    In principle, a higher order integrator could/should be used, but the
    !    trapezoidal rule may not be that bad in this case, actually.  OK,
    !    strictly speaking, this is not even the trapezoidal rule, as the
    !    border terms are not correct.
    !
    !    But then, the borders are in principle open and there should be
    !    essentially no charge outside the logarithmic grid.  This is surely
    !    true at the right hand side, and the error at the left hand side
    !    scales like r_grid_min**3.
    !
    !    Test calculations do not show any difference to
    !    adam_moulton_wave_integrator().
    !
    !  USES

    use constants, only: pi
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    real*8, intent(IN) :: wave(n_grid)
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid(n_grid)
    real*8, intent(OUT) :: coulomb(n_grid)
    real*8, intent(OUT) :: multipole

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o n_grid -- number of logarithmic grid points
    !    o r_grid -- logrithmic grid
    !  OUTPUTS
    !    o coulomb -- Coulomb field of wave on log grid [V(r)]
    !    o multipole -- multipole moment of this function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_grid
    real*8 :: dlnr, prefac
    real*8 :: integral_zero_r, integral_r_infinity
    character(*), parameter :: func = 'trapezoidal_wave_integrator'

    if (n_grid <= 1) call aims_stop('Non-existent grid', func)

    dlnr = log(r_grid(2)/r_grid(1))
    prefac = dlnr * 4*pi / (2*i_l + 1)

    ! --- Outward integration

    integral_zero_r = 0.d0
    do i_grid = 1, n_grid
       integral_zero_r = integral_zero_r + &
       &                 r_grid(i_grid)**(2+i_l) * wave(i_grid)
       coulomb(i_grid) = integral_zero_r * prefac / &
       &                      r_grid(i_grid)**(i_l+1)
    end do

    ! --- Inward integration

    integral_r_infinity = 0.d0
    do i_grid = n_grid, 1, -1
       coulomb(i_grid)= coulomb(i_grid) + &
       integral_r_infinity * prefac * r_grid(i_grid)**i_l

       integral_r_infinity = integral_r_infinity + &
       r_grid(i_grid)**2 * wave(i_grid)/r_grid(i_grid)**(i_l+1)
    enddo

    ! --- Multipole moment

    multipole = dlnr * integral_zero_r / (2*i_l + 1)

  end subroutine trapezoidal_wave_integrator
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/adams_moulton_wave_integrator
  !  NAME
  !    adams_moulton_wave_integrator
  !  SYNOPSIS

  subroutine adams_moulton_wave_integrator(i_l, wave, n_grid, r_grid, &
  &                                        coulomb, multipole)

    !  PURPOSE
    !
    !    Same as trapezoidal_wave_integrator(), but uses a higher order
    !    integrator.  For typical log-meshes, this does not make much of a
    !    difference, though.
    !
    !  USES

    use constants, only: pi4
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    real*8, intent(IN) :: wave(n_grid)
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid(n_grid)
    real*8, intent(OUT) :: coulomb(n_grid)
    real*8, intent(OUT) :: multipole

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o n_grid -- number of logarithmic grid points
    !    o r_grid -- logrithmic grid
    !  OUTPUTS
    !    o coulomb -- Coulomb field of wave on log grid [V(r)]
    !    o multipole -- multipole moment of this function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: dlnr, prefac
    real*8 :: integral_zero_r, integral_r_infty
    integer :: i_grid, n_grid_inf
    character(*), parameter :: func = 'adams_moulton_wave_integrator'

    if (n_grid <= 1) call aims_stop('Non-existent grid', func)

    dlnr = log(r_grid(2)/r_grid(1))
    prefac = pi4/(2*i_l+1) * dlnr

    ! --- forward integral from zero to r

    ! JW: Using f(1)*h is pretty arbitrary.  In principle, we could
    !     extrapolate wave(r) to A*r**(i_l+1) with A=wave(1)/r(1)**(i_l+1) and
    !     then integrate the extrapolation analytically from 0 to r(1).  As
    !     this term will be of the order of r_grid_min**3, the current scheme
    !     is probably sufficient.

    ! first term: f(1)*h
    integral_zero_r = wave(1) * r_grid(1)**(2+i_l)
    coulomb(1) = prefac * integral_zero_r / r_grid(1)**(i_l+1)

    ! second term: (f(1)+f(2))*h/2
    integral_zero_r = integral_zero_r + &
    &                 0.5d0 * (  wave(1) * r_grid(1)**(i_l+2) &
    &                          + wave(2) * r_grid(2)**(i_l+2))

    coulomb(2) = prefac * integral_zero_r / r_grid(2)**(i_l+1)

    ! third term: (5*f(3)+8*f(2)-f(1))*h/12
    integral_zero_r = integral_zero_r + &
    &                 1.d0/12.d0 * (  5.d0 * wave(3) * r_grid(3)**(i_l+2) &
    &                               + 8.d0 * wave(2) * r_grid(2)**(i_l+2) &
    &                               -        wave(1) * r_grid(1)**(i_l+2))

    coulomb(3) = prefac * integral_zero_r / r_grid(3)**(i_l+1)

    ! fourth term and beyond: (9*f(i)+19*f(i-1)-5f(i-2)+f(i-3))*h/24
    do i_grid = 4, n_grid, 1
       integral_zero_r = integral_zero_r + &
       &  1.d0/24.d0 * (   9.d0 * wave(i_grid)   * r_grid(i_grid  )**(i_l+2) &
       &                + 19.d0 * wave(i_grid-1) * r_grid(i_grid-1)**(i_l+2) &
       &                -  5.d0 * wave(i_grid-2) * r_grid(i_grid-2)**(i_l+2) &
       &                +         wave(i_grid-3) * r_grid(i_grid-3)**(i_l+2))

       coulomb(i_grid) = prefac * integral_zero_r / r_grid(i_grid)**(i_l+1)
    end do

    ! --- backward integral from positive infinity to r

    n_grid_inf=n_grid

    ! JW: This is again quite arbitrary.  But if wave(N)/=0, we are probably
    !     doomed, anyway.

    ! first term: f(N)*h
    integral_r_infty = wave(n_grid_inf)/ r_grid(n_grid_inf)**(i_l-1)

    ! second term: (f(N)+f(N-1))*h/2
    coulomb(n_grid_inf) = coulomb(n_grid_inf) + &
    &                     prefac * integral_r_infty*r_grid(n_grid_inf)**i_l

    integral_r_infty = integral_r_infty + &
    &   0.5d0 * (  wave(n_grid_inf)  /r_grid(n_grid_inf  )**(i_l-1) &
    &            + wave(n_grid_inf-1)/r_grid(n_grid_inf-1)**(i_l-1))

    ! third term: (5*f(N-2)+8*f(N-1)-f(N))*h/2
    coulomb(n_grid_inf-1) = coulomb(n_grid_inf-1) + &
    &                  prefac * integral_r_infty*r_grid(n_grid_inf-1)**i_l

    integral_r_infty = integral_r_infty + &
    & 1.d0/12.d0 * (  5.d0*wave(n_grid_inf-2)/r_grid(n_grid_inf-2)**(i_l-1) &
    &               + 8.d0*wave(n_grid_inf-1)/r_grid(n_grid_inf-1)**(i_l-1) &
    &               -      wave(n_grid_inf)  /r_grid(n_grid_inf  )**(i_l-1))

    coulomb(n_grid_inf-2) = coulomb(n_grid_inf-2) + &
    &                  prefac * integral_r_infty*r_grid(n_grid_inf-2)**i_l

    ! fourth term and beyond: (9*f(i)+19*f(i+1)-5*f(i+2)+f(i+3))*h/2
    do i_grid = n_grid-3, 1, -1

       integral_r_infty = integral_r_infty + &
       & 1.d0/24.d0 * (   9.d0*wave(i_grid  )/r_grid(i_grid  )**(i_l-1) &
       &               + 19.d0*wave(i_grid+1)/r_grid(i_grid+1)**(i_l-1) &
       &               -  5.d0*wave(i_grid+2)/r_grid(i_grid+2)**(i_l-1) &
       &               +       wave(i_grid+3)/r_grid(i_grid+3)**(i_l-1))

       coulomb(i_grid) = coulomb(i_grid) + &
       &                 prefac * integral_r_infty*r_grid(i_grid)**i_l

    enddo

    multipole = dlnr * integral_zero_r / (2*i_l + 1)

  end subroutine adams_moulton_wave_integrator
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/hse_wave_integrator
  !  NAME
  !    hse_wave_integrator
  !  SYNOPSIS

  subroutine hse_wave_integrator(i_l, wave, n_grid, n_max_grid, r_grid, &
  &                              hse_matrix, field)

    !  PURPOSE
    !
    !    Similar to trapezoidal_wave_integrator(), but integrates HSE field.
    !    This procedure is not really optimized for performance.
    !
    !  USES

    use constants, only: pi
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    real*8, intent(IN) :: wave(n_max_grid)
    integer, intent(IN) :: n_grid
    integer, intent(IN) :: n_max_grid
    real*8, intent(IN) :: r_grid(n_grid)
    real*8, intent(IN) :: hse_matrix(n_max_grid, n_max_grid)
    real*8, intent(OUT) :: field(n_max_grid)

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o n_grid -- number of logarithmic grid points
    !    o n_max_grid -- array dimension
    !    o r_grid -- logrithmic grid
    !    o hse_matrix -- result from integrate_errorfunction()
    !  OUTPUTS
    !    o field -- Screened coulomb field of wave on log grid [V(r)]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: alpha, hse_correction_term, dr_coef
    integer :: i_grid, j_grid
    character(*), parameter :: func = 'hse_wave_integrator'

    if (n_grid <= 1) call aims_stop('Non-existent grid', func)

    alpha =  log(r_grid(2) / r_grid(1))
    do i_grid = 1, n_grid
       hse_correction_term = 0.d0
       do j_grid = 1, n_grid
          dr_coef = r_grid(j_grid)**2 * alpha
          hse_correction_term = hse_correction_term &
          &                     + wave(j_grid) * hse_matrix(i_grid,j_grid) &
          &                       * dr_coef
       end do

       field(i_grid) = hse_correction_term * 2.0d0 * pi
    enddo

  end subroutine hse_wave_integrator
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/radial_wave_integrator
  !  NAME
  !    radial_wave_integrator
  !  SYNOPSIS

  subroutine radial_wave_integrator(i_l, n_radial, r_radial, w_radial, &
  &                                 n_grid, r_grid_min, r_grid_inc, &
  &                                 wave_spl, field)

    !  PURPOSE
    !
    !    Integrate Coulomb potential on radial (not logarithmic) grid.
    !    This is old legacy code which has not been tested for long and
    !    probably does not work.
    !
    !    I'm pretty sure that this code is buggy (in the sense of inaccurate)
    !    anyway because wave(i_radial) is completely taken into account for
    !    field(i_radial), though the integrator seems to be a mid-point rule,
    !    so that only (roughly?) half of the value should go in.
    !
    !                  -- JW
    !
    !  USES

    use constants, only: pi4
    use grids, only: invert_log_grid
    use dimensions, only: n_max_spline
    use spline
    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    integer, intent(IN) :: n_radial
    real*8, intent(IN) :: r_radial(n_radial), w_radial(n_radial)
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid_min, r_grid_inc
    real*8, intent(IN) :: wave_spl(n_max_spline, n_grid)
    real*8, intent(OUT) :: field(n_radial)
    
    !  INPUTS
    !    o i_l -- angular momentum
    !    o n_radial, r_radial, w_radial -- radial integration grid
    !    o n_grid, n_grid_min, n_grid_inc -- logarithmic grid
    !    o wave_spl -- spline of auxiliary radial function [u(r) = r*P(r)]
    !                  on log-grid
    !  OUTPUTS
    !    o field -- Coulomb field on radial grid
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_radial
    real*8 :: i_r, r_l, r_neg_l1
    real*8 :: integral_zero_r, integral_r_infty
    real*8 :: prefactor 
    real*8, allocatable :: wave(:)
    integer :: info
    character(*), parameter :: func = 'radial_wave_integrator'

    call aims_stop('Procedure not tested; probably does not work', func)

    allocate(wave(n_radial), stat=info)
    call check_allocation(info, 'wave', func)

    do i_radial = 1, n_radial
       i_r = invert_log_grid(r_radial(i_radial), r_grid_min, r_grid_inc) 
       wave(i_radial) = val_spline(i_r, wave_spl, n_grid)
    enddo

    prefactor = pi4/(2*i_l +1)
    integral_zero_r = 0.d0
    do i_radial = 1, n_radial

       r_l = r_radial(i_radial)**i_l
       r_neg_l1 = r_radial(i_radial)**(-i_l-1)

       integral_zero_r = integral_zero_r &
       &               + wave(i_radial) &
       &                 * r_radial(i_radial) * w_radial(i_radial) * r_l

       field(i_radial) = integral_zero_r * r_neg_l1
    enddo

    integral_r_infty = 0.d0
    do i_radial =  n_radial, 1, -1
       r_l = r_radial(i_radial)**i_l
       r_neg_l1 = r_radial(i_radial)**(-i_l-1)

       field(i_radial) = field(i_radial) + integral_r_infty * r_l

       integral_r_infty = integral_r_infty &
       &                + wave(i_radial) &
       &                  * r_radial(i_radial) * w_radial(i_radial) * r_neg_l1

       field(i_radial) = field(i_radial)*prefactor
    enddo

  end subroutine radial_wave_integrator
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/inspect_wave_coulomb
  !  NAME
  !    inspect_wave_coulomb
  !  SYNOPSIS

  subroutine inspect_wave_coulomb(i_l, wave, n_grid, r_grid, wave_threshold, &
  &                               charge_radius, field_radius, multipole)

    !  PURPOSE
    !
    !    Examine a given (product basis) radial function [u(r) = r\psi(r)] and
    !    return the extent of the plain wave (charge_radius) as well as its
    !    Coulomb filed (field_radius) and its multipole moment.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    real*8, intent(IN) :: wave(n_grid)
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid(n_grid)
    real*8, intent(IN) :: wave_threshold
    real*8, intent(OUT) :: charge_radius
    real*8, intent(OUT) :: field_radius
    real*8, intent(OUT) :: multipole

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o n_grid -- number of logarithmic grid points
    !    o r_grid -- logrithmic grid
    !    o wave_threshold -- threshold for radii (typically ~1d-6)
    !  OUTPUTS
    !    o charge_radius -- radius at which the auxiliary function is
    !                       wave_threshold
    !    o field_radius -- radius at which the Coulomb potential of
    !                      the auxiliary function is wave_threshold
    !    o multipole -- multipole moment of this aux function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: coulomb(n_grid)
    character(*), parameter :: func = 'inspect_wave_coulomb'

    ! --- multipole

    call adams_moulton_wave_integrator(i_l, wave, n_grid, r_grid, &
    &                                  coulomb, multipole)

    ! --- radii

    call inspect_wave_field(i_l, wave, coulomb, multipole, n_grid, r_grid, &
    &                       wave_threshold, charge_radius, field_radius)

  end subroutine inspect_wave_coulomb
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/inspect_wave_field
  !  NAME
  !    inspect_wave_field
  !  SYNOPSIS

  subroutine inspect_wave_field(i_l, wave, field, multipole, n_grid, r_grid, &
  &                             wave_threshold, charge_radius, field_radius)

    !  PURPOSE
    !
    !    Examine a given (product basis) radial function [u(r) = r\psi(r)] and
    !    return the extent of the plain wave (charge_radius) as well as its
    !    Coulomb filed (field_radius) and its multipole moment.
    !
    !  USES

    use constants, only: pi
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    real*8, intent(IN) :: wave(n_grid), field(n_grid)
    real*8, intent(IN) :: multipole
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid(n_grid)
    real*8, intent(IN) :: wave_threshold
    real*8, intent(OUT) :: charge_radius
    real*8, intent(OUT) :: field_radius

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o field -- its (possibly screened) Coulomb field
    !    o multipole -- moment for far field (if unscreened; otherwise 0.d0)
    !    o n_grid -- number of logarithmic grid points
    !    o r_grid -- logrithmic grid
    !    o wave_threshold -- threshold for radii (typically ~1d-6)
    !  OUTPUTS
    !    o charge_radius -- radius at which the auxiliary function is
    !                       wave_threshold
    !    o field_radius -- radius at which the field of
    !                      the auxiliary function is wave_threshold
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_grid
    real*8 :: rmax
    character(*), parameter :: func = 'inspect_wave_field'

    ! --- radii

    charge_radius = 0.d0
    field_radius = 0.d0
    do i_grid = 1, n_grid
       if (abs(wave(i_grid)) > wave_threshold**2) then
          charge_radius = r_grid(i_grid)
       end if
       if (abs(field(i_grid)) > wave_threshold**2) then
          field_radius = r_grid(i_grid)
       end if
    end do
    ! JW: For now, also use wave_threshold here to be on the safe side.  We do
    ! not make any use of locality right now, anyway.  As soon as we do, we
    ! have to revisit this parameter.


    ! --- update field_radius by multipole

    ! V_far(rmax) = 4*pi * multipole / rmax**(i_l+1) == wave_threshold**2

    ! I use wave_threshold**2 in this case because the situation for multipole
    ! potential times wave function is different from the ordinary overlap
    ! case.  For overlap, a criterion of 1d-6 is because both functions decay
    ! roughly at the same rate and probably drop to 1d-12 at each others core.
    ! The multipole potential, on the other hand, does not decay noteworthy
    ! between a basis function edge and its core.

    ! Won't change anything for HSE or multipole-less fields.
    rmax = (4*pi*abs(multipole) / wave_threshold**2)**(1.d0 / dble(i_l+1))
    field_radius = max(rmax, charge_radius, field_radius)

  end subroutine inspect_wave_field
  !******
  !----------------------------------------------------------------------------
  !****s* basbas_fn_coulomb/localize_basbas_fn_coulomb
  !  NAME
  !    localize_basbas_fn_coulomb
  !  SYNOPSIS

  subroutine localize_basbas_fn_coulomb(i_l, i_first_fn, i_last_fn, &
  &                                     basbas_wave, n_grid, r_grid, &
  &                                     use_AM, have_hse)

    !  PURPOSE
    !
    !    Combine product basis functions (basbas_fn) of this l-channel
    !    linearly to localize their Coulomb potential as well as possible.
    !    Described by M. Betzinger et al. [1].  This is similar in spirit to
    !    the "Poisson trick" mentioned in [2].
    !
    !    The far field of an atomic function of given (l, m) is exactly that
    !    of a multipole of corresponding lm.  Therefore, all Coulomb tails of
    !    the auxiliary basis functions of a given lm are linear dependent.  By
    !    taking linear combinations, the tail can be attributed to one
    !    auxiliary function per lm channel, effectively localizing the Coulomb
    !    field of all other functions.
    !
    !    Localization is optimized for the case that the original product basis
    !    functions are roughly sorted by their (ordinary) extent.
    !
    !    If (have_hse) then also try to concentrate higher momenta in the
    !    later radial parts.  This can only work as long as
    !    outer_radius*hse_omega_hf is smaller or at least near to one.
    !    Fortunately, it cannot do any harm apart from a potential loss of the
    !    charge-localization of the core auxiliary functions because the
    !    transform is orthogonal.  Why this should help can be understood from
    !    the series expansions in [3], although the authors do not discuss
    !    this particular set of parameters (mu*r<<1, mu*R>>1), mu=omega.
    !
    !    In practice, this procedure does not give well-localized screened
    !    fields, but strongly reduces the magnitude of fields for outer_radius
    !    < R < N/omega, where N is something like 5.  As the criterion for the
    !    field readius is rather strict, the localization procedure does not
    !    have too much influence on the field radii.
    !
    !  USES

    use dimensions, only: n_max_grid
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_l
    integer, intent(IN) :: i_last_fn, i_first_fn
    real*8, intent(INOUT) :: basbas_wave(n_max_grid, i_last_fn)
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid (n_grid)
    logical, intent(IN) :: use_AM
    logical, intent(IN) :: have_hse

    !  INPUTS
    !    o i_l -- current angular momentum channel
    !    o i_first_fn -- First auxiliary radial function
    !    o i_last_fn -- Last auxiliary radial function
    !    o basbas_wave -- the value auxiliary radial function [u(r) = r*P(r)]
    !    o n_grid -- number of logarithmic grid points
    !    o r_grid -- logrithmic grid
    !    o use_AM -- use Adoms-Moulton integrator (or trapezoidal?)
    !    o have_hse -- if this is to be used for HSE, continue localizing
    !                  higher momenta
    !  OUTPUTS
    !    o basbas_wave -- updated radial parts
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    [1] M. Betzinger, C. Friedrich, S. Bluegel: "Hybrid functionals
    !    within the all-electron FLAPW method: Implementation and applications
    !    to PBE0", Phys. Rev. B 81, 195117 (2010).
    !    [2] M. Katouda, S. Nagase, "Application of second-order
    !    Moller-Plesset perturbation theory with resolution-of-identity
    !    approximation to periodic systems", J. Chem. Phys. 133, 184103
    !    (2010).
    !    [3] J. G. Angyan, I. Gerber, M. Marsman, "Spherical harmonic
    !    expansion of short-range screened Coulomb interactions", J. Phys. A:
    !    Math. Gen. 39, 8613 (2006).
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_p, i_p, p, n_fn
    integer :: i_basbas_fn, i_grid
    real*8 :: mpm_prev, mpm_this, inv_norm
    real*8 :: mp_free, mp_all
    character(*), parameter :: func = 'localize_basbas_fn_coulomb'

    if (n_grid > n_max_grid) call aims_stop('n_grid > n_max_grid', func)
    if (have_hse) then
       n_fn = i_last_fn - i_first_fn + 1
       ! Localize as far as possible.  This is probably too much, in practice,
       ! as we don't gain much in the field for the later momenta but loose
       ! localization in the charge for the core ABFs.
       n_p = n_fn - 1
    else
       ! This is the general case to be understood first: Only localize pure
       ! multipole moment.
       n_p = 1
    end if

    ! --- Unitary transformation

    do i_p = 1, n_p
       p = 2 * (i_p-1)
       do i_basbas_fn = i_first_fn + 1, i_last_fn - (i_p - 1)
          call get_mp_moment(i_l + p, basbas_wave(:,i_basbas_fn-1), mpm_prev)
          call get_mp_moment(i_l + p, basbas_wave(:,i_basbas_fn), mpm_this)
          inv_norm = 1.d0 / sqrt(mpm_prev**2 + mpm_this**2)

          do i_grid = 1, n_grid
             mp_all = (mpm_prev * basbas_wave(i_grid, i_basbas_fn-1) + &
             &         mpm_this * basbas_wave(i_grid, i_basbas_fn))
             mp_free = (mpm_this * basbas_wave(i_grid, i_basbas_fn-1) -&
             &          mpm_prev * basbas_wave(i_grid, i_basbas_fn))
             basbas_wave(i_grid, i_basbas_fn-1) = inv_norm * mp_free
             basbas_wave(i_grid, i_basbas_fn) = inv_norm * mp_all
          end do
       end do
    end do

    return

  contains

    !****s* basbas_fn_coulomb/localize_fn_coulomb/get_mp_moment
    !  NAME
    !    get_mp_moment
    !  SYNOPSIS
    subroutine get_mp_moment(eff_L, wave, mpm)
      !  PURPOSE
      !    Calclulate something proportional to \int dr f(r)*r^(2+eff_L).
      !    For this purpose, the wave integrators are misused in two ways:
      !    First, the actual potential is thrown away.  Second, eff_L is not
      !    truly an angular momentum.
      !  USES
      implicit none
      !  ARGUMENTS
      integer, intent(IN) :: eff_L
      real*8, intent(IN) :: wave(n_max_grid)
      real*8, intent(OUT) :: mpm
      !  INPUTS
      !    o eff_L -- "effective angular momentum", specifies power
      !    o wave -- u(r) = r*f(r)
      !  OUTPUTS
      !    o mpm -- "multipole moment": \int dr f(r)*r^(2+eff_L).
      !  SOURCE
      real*8 :: coulomb_tmp(n_grid)
      if (use_AM) then
         call adams_moulton_wave_integrator(eff_L, wave, n_grid, r_grid, &
         &                                  coulomb_tmp, mpm)
      else
         call trapezoidal_wave_integrator(eff_L, wave, n_grid, r_grid, &
         &                                coulomb_tmp, mpm)
      end if
    end subroutine get_mp_moment
    !******

  end subroutine localize_basbas_fn_coulomb
  !******
  !----------------------------------------------------------------------------
  !****s* FHI-aims/sort_species_basbas_fn
  !  NAME
  !    sort_species_basbas_fn
  !  SYNOPSIS

  subroutine sort_species_basbas_fn(n_species_fn, n_grid, r_grid, &
  &                                 wave_threshold, wave, fn_l)

    !  PURPOSE
    !    Sort the radial parts of the product basis functions of one species
    !    by increasing extent.
    !  USES

    use dimensions, only: n_max_grid
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_species_fn
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: r_grid(n_max_grid)
    real*8, intent(IN) :: wave_threshold
    real*8, intent(INOUT) :: wave(n_max_grid, n_species_fn)
    integer, intent(INOUT) :: fn_l(n_species_fn)

    !  INPUTS
    !    o n_species_fn -- Number of product basis functions for this species
    !    o n_grid, r_grid -- Logarithmic grid
    !    o wave -- Radial waves [u(r) = rP(r)] on log grid
    !    o fn_l -- Angular momenta
    !  OUTPUTS
    !    o wave -- Sorted waves
    !    o fn_l -- Correspondingly permutated angular momenta
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: outer_radii(n_species_fn)
    real*8 :: charge_radius, multipole
    integer :: old2new(n_species_fn), new2old(n_species_fn)
    integer :: i_fn, i_l
    real*8 :: wave_tmp(n_max_grid, n_species_fn)
    integer :: fn_l_tmp(n_species_fn)

    character(*), parameter :: func = 'sort_species_basbas_fn'

    if (n_grid > n_max_grid) call aims_stop('n_grid > n_max_grid', func)

    do i_fn = 1, n_species_fn
       call inspect_wave_coulomb(fn_l(i_fn), wave(:, i_fn), n_grid, r_grid, &
       &                         wave_threshold, &
       &                         charge_radius, outer_radii(i_fn), multipole)
    end do

    call insertionsort(outer_radii, n_species_fn, old2new, new2old)

    do i_fn = 1, n_species_fn
       wave_tmp(:, i_fn) = wave(:, new2old(i_fn))
       fn_l_tmp(i_fn) = fn_l(new2old(i_fn))
    end do
    wave = wave_tmp
    fn_l = fn_l_tmp

  end subroutine sort_species_basbas_fn
  !******
end module basbas_fn_coulomb
!******
