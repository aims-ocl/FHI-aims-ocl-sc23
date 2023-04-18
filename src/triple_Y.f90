!****h* FHI-aims/triple_Y
!  NAME
!    triple_Y
!  SYNOPSIS

module triple_Y

  !  PURPOSE
  !    Calculate angular integrals containing three spherical harmonics,
  !    both canonical and real-valued.
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
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  integer :: max_L_YrLM = -1

  integer :: max_L_act

  real*8, allocatable :: YrLM_s(:)
  real*8, allocatable :: d_YrLM_s(:,:)

  type CrMs
    real*8, allocatable :: Cr(:,:,:)
    integer, allocatable :: Ms(:,:,:)
  end type

  type (CrMs), allocatable :: CrMs_s(:,:,:)

contains
  !----------------------------------------------------------------------------
  !****s* triple_Y/calc_YrLM_for_triple_Y
  !  NAME
  !    calc_YrLM_for_triple_Y
  !  SYNOPSIS

  subroutine calc_YrLM_for_triple_Y(Rvec, max_L)

    !  PURPOSE
    !
    !     Precalculates the YrLM values for a given Rvec
    !     (to be used in fast_sum_triple_Y_YLM_real)
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: Rvec(3)
    integer, intent(IN) :: max_L

    !  INPUTS
    !    o Rvec -- to define the direction within the spherical harmonic
    !    o max_L -- maximum L value for all basbas functions.
    !      Please note: If max_L > 10 is needed anywhere in the code,
    !      this must be specified on the first call to this routine.
    !      This is ugly and should be fixed sometimes ...
    !  OUTPUTS
    !    None, result is stored locally in YrLM_s
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'calc_YrLM_for_triple_Y'

    if(max_L_YrLM < 0) then
      ! Initialize: allocate YrLM_s and CrMs_s big enough for max_L
      max_L_YrLM = MAX(max_L,15)
      allocate(YrLM_s((2*max_L_YrLM+1)*(2*max_L_YrLM+1)))
      allocate(d_YrLM_s(3,(2*max_L_YrLM+1)*(2*max_L_YrLM+1)))
      allocate(CrMs_s(0:max_L_YrLM,0:max_L_YrLM,0:2*max_L_YrLM))
    else
      ! Not the initial call, for safety check max_L
      if(max_L > max_L_YrLM) call aims_stop('INTERNAL ERROR: max_L > max_L_YrLM', func)
    endif

    ! Calculate YLM up to L = 2*max_L
    call ylm_real_deriv(Rvec, 2*max_L, YrLM_s, d_YrLM_s)

    max_L_act = max_L

  end subroutine calc_YrLM_for_triple_Y

  !----------------------------------------------------------------------------
  !****s* triple_Y/fast_sum_triple_Y_YLM_real
  !  NAME
  !    fast_sum_triple_Y_YLM_real
  !  SYNOPSIS

  subroutine fast_sum_triple_Y_YLM_real(l1, l2, L, max_l1, max_l2, cyl, d_cyl, calc_deriv)

    !  PURPOSE
    !
    !     Calculate the sums
    !
    !        cyl(m1, m2) := \sum_{M=-L}^L C(l1, l2, L, m1, m2, M) Y_{LM}(R)
    !
    !     which are needed within the calculation of overlap integrals using
    !     spherical Bessel transforms.
    !
    !     Optionally calculates also the derivatives with respect to Rvec.
    !
    !     For this (fast_...) version, Y_{LM}(R) must have been precalculated
    !     by calc_YrLM_for_triple_Y
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l1, l2, L
    integer, intent(IN) :: max_l1, max_l2
    real*8, intent(INOUT) :: cyl(-max_l1:max_l1, -max_l2:max_l2)
    real*8, intent(INOUT) :: d_cyl(-max_l1:max_l1, -max_l2:max_l2, 3)
    logical, intent(IN) :: calc_deriv

    !  INPUTS
    !    o l1, l2, L -- Angular momentum quantum numbers
    !    o max_l1, max_l2 -- dimensions of cyl(:,:)
    !    o calc_deriv -- flag if the derivatives should also be calculated
    !  OUTPUTS
    !    o cyl -- The resulting sums.
    !    o d_cyl -- Derivatives of cyl, not touched if calc_deriv = FALSE
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'fast_sum_triple_Y_YLM_real'

    integer :: M, m1, m2, M_type, m_off

    if(L > 2*max_L_act) call aims_stop('INTERNAL ERROR: L > 2*max_L_act', func)
    if(l1 >  max_L_act) call aims_stop('INTERNAL ERROR: l1 >  max_L_act', func)
    if(l2 >  max_L_act) call aims_stop('INTERNAL ERROR: l2 >  max_L_act', func)

    ! --- Get Cr/Ms for given l1, l2, L
    !     These are not depending on the geometry (Rvec) and thus stored for
    !     further usage once they have been calculated

    if(.not.allocated(CrMs_s(l1, l2, L)%Cr)) then
      allocate(CrMs_s(l1, l2, L)%Cr(-l1:l1, -l2:l2, 2))
      allocate(CrMs_s(l1, l2, L)%Ms(-l1:l1, -l2:l2, 2))
      call triple_Y_real(l1, l2, L, l1, l2, CrMs_s(l1, l2, L)%Cr, CrMs_s(l1, l2, L)%Ms)
    endif

    ! --- Do sums

    m_off = L*L + L+1
    cyl = 0.d0
    do m2 = -l2, l2
       do m1 = -l1, l1
          M = CrMs_s(l1, l2, L)%Ms(m1, m2, 1)
          if (abs(M) <= L) then
             cyl(m1, m2) = cyl(m1, m2) + CrMs_s(l1, l2, L)%Cr(m1, m2, 1) * YrLM_s(M+m_off)
          end if
          M = CrMs_s(l1, l2, L)%Ms(m1, m2, 2)
          if (abs(M) <= L) then
             cyl(m1, m2) = cyl(m1, m2) + CrMs_s(l1, l2, L)%Cr(m1, m2, 2) * YrLM_s(M+m_off)
          end if
       end do
    end do

    ! --- Same for derivatives

    if(calc_deriv) then
       d_cyl = 0.d0
       do m2 = -l2, l2
          do m1 = -l1, l1
             M = CrMs_s(l1, l2, L)%Ms(m1, m2, 1)
             if (abs(M) <= L) then
                d_cyl(m1, m2, 1) = d_cyl(m1, m2, 1) + CrMs_s(l1, l2, L)%Cr(m1, m2, 1) * d_YrLM_s(1, M+m_off)
                d_cyl(m1, m2, 2) = d_cyl(m1, m2, 2) + CrMs_s(l1, l2, L)%Cr(m1, m2, 1) * d_YrLM_s(2, M+m_off)
                d_cyl(m1, m2, 3) = d_cyl(m1, m2, 3) + CrMs_s(l1, l2, L)%Cr(m1, m2, 1) * d_YrLM_s(3, M+m_off)
             end if
             M = CrMs_s(l1, l2, L)%Ms(m1, m2, 2)
             if (abs(M) <= L) then
                d_cyl(m1, m2, 1) = d_cyl(m1, m2, 1) + CrMs_s(l1, l2, L)%Cr(m1, m2, 2) * d_YrLM_s(1, M+m_off)
                d_cyl(m1, m2, 2) = d_cyl(m1, m2, 2) + CrMs_s(l1, l2, L)%Cr(m1, m2, 2) * d_YrLM_s(2, M+m_off)
                d_cyl(m1, m2, 3) = d_cyl(m1, m2, 3) + CrMs_s(l1, l2, L)%Cr(m1, m2, 2) * d_YrLM_s(3, M+m_off)
             end if
          end do
       end do
    endif

  end subroutine fast_sum_triple_Y_YLM_real
  !******
  !----------------------------------------------------------------------------
  !****s* triple_Y/sum_triple_Y_YLM_real
  !  NAME
  !    sum_triple_Y_YLM_real
  !  SYNOPSIS

  subroutine sum_triple_Y_YLM_real(l1, l2, L, Rvec, max_l1, max_l2, cyl)

    !  PURPOSE
    !
    !     Calculate the sums
    !
    !        cyl(m1, m2) := \sum_{M=-L}^L C(l1, l2, L, m1, m2, M) Y_{LM}(R)
    !
    !     which are needed within the calculation of overlap integrals using
    !     spherical Bessel transforms.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l1, l2, L
    real*8, intent(IN) :: Rvec(3)
    integer, intent(IN) :: max_l1, max_l2
    real*8, intent(OUT) :: cyl(-max_l1:max_l1, -max_l2:max_l2)

    !  INPUTS
    !    o l1, l2, L -- Angular momentum quantum numbers
    !    o max_l1, max_l2 -- dimensions of cyl(:,:)
    !    o Rvec -- to define the direction within the spherical harmonic
    !  OUTPUTS
    !    o cyl -- The resulting sums.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: YrLM((L+1)*(L+1))
    real*8 :: YY(-L:L)
    real*8 :: Cr(-l1:l1, -l2:l2, 2)
    integer :: Ms(-l1:l1, -l2:l2, 2)
    integer :: M, m1, m2, M_type

    ! --- Get Y_LM for given L:

    call ylm_real(Rvec, L, YrLM)
    do M = -L, L
       YY(M) = YrLM(L*L + L+M+1)
    end do

    ! --- Get Cr for given l1, l2, L

    call triple_Y_real(l1, l2, L, l1, l2, Cr, Ms)

    ! --- Do sums

    cyl = 0.d0
    do m1 = -l1, l1
       do m2 = -l2, l2
          do M_type = 1, 2
             M = Ms(m1, m2, M_type)
             if (abs(M) <= L) then
                cyl(m1, m2) = cyl(m1, m2) + Cr(m1, m2, M_type) * YY(M)
             end if
          end do
       end do
    end do

  end subroutine sum_triple_Y_YLM_real
  !******
  !----------------------------------------------------------------------------
  !****s* triple_Y/triple_Y_real
  !  NAME
  !    triple_Y_real
  !  SYNOPSIS

  subroutine triple_Y_real(l1, l2, L, max_l1, max_l2, triple_Yr, Ms)

    !  PURPOSE
    !
    !    Calculate the integrals C^r(l1, l2, L, m1, m2, M)
    !
    !        C^r = \int_\Omega d\Omega Y^r_{l1,m1} Y^r_{l2,m2} Y^r_{L,M}
    !
    !    for given l1, l2, L and for all m1, m2, and matching M.
    !    Y^r are the real-valued spherical harmonics as defined in ylm_real.f.
    !
    !    Essentially, uses equation (4.9) of the Talman paper cited below to
    !    calculate the triple Ylm integrals for real-valued spherical harmonics.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l1, l2, L
    integer, intent(IN) :: max_l1, max_l2
    real*8, intent(OUT) :: triple_Yr(-max_l1:max_l1, -max_l2:max_l2, 2)
    integer, intent(OUT) :: Ms(-max_l1:max_l1, -max_l2:max_l2, 2)

    !  INPUTS
    !    o l1, l2, L -- Angular momentum quantum numbers
    !    o max_l1, max_l2 -- Array dimensions
    !  OUTPUTS
    !    o triple_Yr(m1, m2, Ms(i)) contain C^r(l1, l2, L, m1, m2, Ms(i))
    !          C^r(..., M) == 0 if M not in Ms(:).
    !    o Ms(m1, m2, :) -- M with non-vanishing C^r.  M=-L-1 if no value.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !
    !    James D. Talman,
    !    "Numerical methods for multicenter integrals for numerically defined
    !     basis functions applied in molecular calculations",
    !    Internat. J. Quant. Chem. 93, 72 (2003).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: m1, m2, M_type, M, Mabs
    real*8 :: triple_Y(-l1:l1, -l2:l2), this_C
    integer :: mm(3), abs_m_order(3), i, j, tmp_int, sign_m(3)
    integer :: Talman_sign, aims_sign
    real*8 :: Talman_w3, prefac

    ! I have tried to make this subroutine as legible as possible.
    ! Unfortunately, the subject is very subtle.  For example, it took me
    ! quite some time to note that these integrals vanish if an odd number
    ! of ms is negative.  Talman obviously knew that and took it for granted as
    ! its formula returns garbage in this case.

    ! Initialize
    triple_Yr = 0.d0
    Ms = -L-1  ! invalid

    ! Get triple_Y (integrals of canonical spherical harmonics)
    call triple_Y_cmplx(l1, l2, L, triple_Y)

    ! Main loop
    do m1 = -l1, l1
       do m2 = -l2, l2
          do M_type = 1, 2

             ! --- Get M

             ! The absolute values of m1, m2, M must be combinable to 0
             ! for a nonzero integral.  -> |M| \in {||m1|-|m2||, |m1|+|m2|}
             if (M_type == 1) then
                Mabs = abs(abs(m1) - abs(m2))
             else if (abs(m1) > 0 .and. abs(m2) > 0) then
                Mabs = abs(m1) + abs(m2)
             else
                cycle  ! Already done in M_type == 1.
             end if
             if (Mabs > L) cycle

             ! For an odd number of negative ms, the resulting integral would
             ! be purely imaginary.  As we know it is real, it would be 0.
             ! From |M|, -|M|, only use M with even number of negs in m1,m2,M.
             if (m1 < 0 .neqv. m2 < 0) then
                if (Mabs == 0) cycle    ! Cannot get even number of negs
                M = - Mabs
             else
                M = Mabs
             end if
             ! And save index:
             Ms(m1, m2, M_type) = M


             ! --- Bubble sort for abs(m).

             mm(1) = m1
             mm(2) = m2
             mm(3) = M
             do i = 1, 3
                abs_m_order(i) = i
             end do
             do i = 1, 3
                do j = 3, i+1, -1
                   if (abs(mm(abs_m_order(j))) > abs(mm(abs_m_order(i)))) then
                      tmp_int = abs_m_order(j)
                      abs_m_order(j) = abs_m_order(i)
                      abs_m_order(i) = tmp_int
                   end if
                end do
             end do
             ! Now: abs(mm(abs_m_order(1))) >= abs(mm(abs_m_order(2))) >= ...

             ! --- Pick right element from triple_Y

             ! As can be seen from (4.9) in the Talman paper cited above, the m
             ! with the largest magnitude (m1 in (4.9)) should have negative
             ! sign, the other two positive signs in C(...).
             sign_m(abs_m_order(1)) = -1
             sign_m(abs_m_order(2)) = 1
             sign_m(abs_m_order(3)) = 1
             this_C = triple_Y(sign_m(1) * abs(m1), sign_m(2) * abs(m2))


             ! --- Get signs and prefactors

             Talman_sign = (-1)**mm(abs_m_order(1))
             if (mm(abs_m_order(2)) < 0 .and. mm(abs_m_order(3)) < 0) then
                Talman_sign = - Talman_sign
             end if
             if (mm(abs_m_order(3)) == 0) then
                Talman_w3 = 1.d0
             else
                Talman_w3 = 1.d0 / sqrt(2.d0)
             end if
             ! The definition of real-valued spherical harmonics differ between
             ! aims and Talman for m < 0 by the sign (-1)**m.  Correct for this:
             aims_sign = 1
             do i = 1, 3
                if (mm(i) < 0) aims_sign = aims_sign * (-1)**mm(i)
             end do

             ! --- Put together.

             prefac = aims_sign * Talman_sign * Talman_w3
             triple_Yr(m1, m2, M_type) = prefac * this_C

          end do
       end do
    end do

  end subroutine triple_Y_real
  !******
  !----------------------------------------------------------------------------
  !****s* triple_Y/triple_Y_cmplx
  !  NAME
  !    triple_Y_cmplx
  !  SYNOPSIS

  subroutine triple_Y_cmplx(l1, l2, L, triple_Y)

    !  PURPOSE
    !
    !    Calculate the integrals C(l1, l2, L, m1, m2, M)
    !
    !        C = \int_\Omega d\Omega Y_{l1,m1} Y_{l2,m2} Y_{L,M}
    !          = sqrt((2*l1+1)*(2*l2+1)*(2*L+1)/(4*pi)) *
    !            (( l1 l2 L ; 0 0 0 )) * (( l1 l2 L ; m1 m2 M ))
    !
    !    from the 3jm symbols (( l1 l2 l3 ; m1 m2 m3 )) as calculated in
    !    drc3jj.f for given l1, l2, L and for all m1, m2 for the canonical
    !    (complex) valued spherical harmonics.
    !
    !  USES

    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: l1, l2, L
    real*8, intent(OUT) :: triple_Y(-l1:l1, -l2:l2)

    !  INPUTS
    !    o l1, l2, L -- Angular momentum quantum numbers
    !  OUTPUTS
    !    o triple_Y(m1, m2) contains C(l1, l2, L, m1, m2, -(m1+m2))
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !     Computer Physics Communications 180, 2175 (2009).
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE


    integer :: m1, m2, m2min, m2max, ierr
    real*8 :: rl1, rl2, rl3, rm1, rm2min, rm2max
    real*8 :: thrcof(2*l2+1)
    real*8 :: prefac, C000
    real*8 :: pi
    character(*), parameter :: func = 'triple_Y_cmplx'

    if (modulo(l1 + l2 + L, 2) /= 0) then
       call aims_stop('Called with odd sum of L', func)
    end if

    pi = 4.d0 * atan(1.d0)
    triple_Y = 0.d0

    ! Prepare (real) parameters for drc3jm.
    rl1 = l1
    rl2 = l2
    rl3 = L

    ! Get 3jm symbols (l1 l2 L; m1 : : )
    do m1 = -l1, l1
       rm1 = m1
       call drc3jm(rl1, rl2, rl3, rm1, &
       &           rm2min, rm2max, thrcof, size(thrcof), ierr)
       m2min = nint(rm2min)
       m2max = nint(rm2max)
       do m2 = m2min, m2max
          triple_Y(m1, m2) = thrcof(m2 - m2min + 1)
       end do
    end do

    ! Prefactors
    C000 = triple_Y(0, 0)
    prefac = sqrt((2*l1+1) * (2*l2+1) * (2*L+1) / (4*pi))
    triple_Y = prefac * C000 * triple_Y

  end subroutine triple_Y_cmplx
  !******
  !----------------------------------------------------------------------------
  !****s* triple_Y/ylm_real_deriv
  !  NAME
  !    triple_Y_cmplx
  !  SYNOPSIS

  subroutine ylm_real_deriv(V,LMAX,Y,DY)

    !  PURPOSE
    !
    !    This is yet another variation of the YLM code which calculates also
    !    the derivative - maybe the different pieces of code dealing with
    !    spherical harmonics should be consolitated sometimes ...
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(in) :: LMAX
    real*8, intent(in)  :: V(3)
    real*8, intent(out) :: Y((LMAX+1)*(LMAX+1)), DY(3,(LMAX+1)*(LMAX+1))
!
!
      INTEGER            L, M
      INTEGER            idx, idx_1, idx_2
      DOUBLE PRECISION   A, B, C, ABC, ABCMAX
      DOUBLE PRECISION   TEMP1, TEMP2, TEMP3, TEMPA, TEMPB, PI
      DOUBLE PRECISION   X_N, Y_N, Z_N, DX_N(3), DY_N(3), DZ_N(3)

      IF(LMAX < 0) RETURN
!
      PI = (4.0D+0)*ATAN(1.0D+0)
!
!        Y(0,0)
!
      Y(1) = 1.0D+0/SQRT(4.0D+0*PI)
      DY(:,1) = 0.

      IF (LMAX == 0) RETURN

      ABCMAX = MAX(ABS(V(1)),ABS(V(2)),ABS(V(3)))
      IF (ABCMAX .GT. 0.0D+0) THEN
         A = V(1)/ABCMAX
         B = V(2)/ABCMAX
         C = V(3)/ABCMAX
         ABC = SQRT(A*A + B*B + C*C)
         X_N = A/ABC
         Y_N = B/ABC
         Z_N = C/ABC
         TEMP1 = (A*A + B*B + C*C)*ABC*ABCMAX
         DX_N = (/ B*B+C*C, -A*B, -A*C /) / TEMP1
         DY_N = (/ -B*A, A*A+C*C, -B*C /) / TEMP1
         DZ_N = (/ -C*A, -C*B, A*A+B*B /) / TEMP1
      ELSE
         ! This is an error - we continue with the following settings:
         X_N = 0.
         Y_N = 0.
         Z_N = 1.
         DX_N = (/ 1., 0., 0. /)
         DY_N = (/ 0., 1., 0. /)
         DZ_N = (/ 0., 0., 0. /)
      ENDIF
!
!        Y(1,x)
!
      TEMP1 = SQRT(3.0D+0)*Y(1)
      Y(2) =  TEMP1*Y_N
      Y(3) =  TEMP1*Z_N
      Y(4) = -TEMP1*X_N

      DY(:,2) =  TEMP1*DY_N(:)
      DY(:,3) =  TEMP1*DZ_N(:)
      DY(:,4) = -TEMP1*DX_N(:)
!
      idx_2 = 1
      idx_1 = 3

      DO L = 2, LMAX
         idx = L*L+L+1
!
!        Y(L,L) = f(Y(L-1,L-1)) ... Formula 1
!
         TEMP1 = -SQRT(DBLE(2*L+1)/DBLE(2*L))
         Y(idx+l) = TEMP1*(Y_N*Y(idx_1-l+1) + X_N*Y(idx_1+l-1))
         Y(idx-l) = TEMP1*(Y_N*Y(idx_1+l-1) - X_N*Y(idx_1-l+1))

         DY(:,idx+l) = TEMP1*(DY_N*Y(idx_1-l+1) + Y_N*DY(:,idx_1-l+1) + DX_N*Y(idx_1+l-1) + X_N*DY(:,idx_1+l-1))
         DY(:,idx-l) = TEMP1*(DY_N*Y(idx_1+l-1) + Y_N*DY(:,idx_1+l-1) - DX_N*Y(idx_1-l+1) - X_N*DY(:,idx_1-l+1))
!
!        Y(L,L-1) = f(Y(L-1,L-1)) ... Formula 2
!        (the coefficient for Y(L-2,L-1) in Formula 2 is zero)
!
         TEMP2 = SQRT(DBLE(2*L+1))
         Y(idx+l-1) = TEMP2*Z_N*Y(idx_1+l-1)
         Y(idx-l+1) = TEMP2*Z_N*Y(idx_1-l+1)

         DY(:,idx+l-1) = TEMP2*(DZ_N*Y(idx_1+l-1) + Z_N*DY(:,idx_1+l-1))
         DY(:,idx-l+1) = TEMP2*(DZ_N*Y(idx_1-l+1) + Z_N*DY(:,idx_1-l+1))
!
         TEMPA = SQRT(DBLE(4*L*L-1))
         TEMPB = -SQRT(DBLE(2*L+1)/DBLE(2*L-3))
!
         DO M = L - 2, 0, -1
!
!        Y(L,M) = f(Y(L-2,M),Y(L-1,M)) ... Formula 2
!
            TEMP1 = 1.0D+0/SQRT(DBLE((L+M)*(L-M)))
            TEMP2 = TEMPA*TEMP1
            TEMP3 = TEMPB*SQRT(DBLE((L+M-1)*(L-M-1)))*TEMP1

            Y(idx+m) = TEMP2*Z_N*Y(idx_1+m) + TEMP3*Y(idx_2+m)
            Y(idx-m) = TEMP2*Z_N*Y(idx_1-m) + TEMP3*Y(idx_2-m)

            DY(:,idx+m) = TEMP2*(DZ_N*Y(idx_1+m) + Z_N*DY(:,idx_1+m)) + TEMP3*DY(:,idx_2+m)
            DY(:,idx-m) = TEMP2*(DZ_N*Y(idx_1-m) + Z_N*DY(:,idx_1-m)) + TEMP3*DY(:,idx_2-m)

         ENDDO

         idx_2 = idx_1
         idx_1 = idx

         Y_N = -Y_N
         DY_N(:) = -DY_N(:)

      ENDDO

  end subroutine ylm_real_deriv
  !******
  !----------------------------------------------------------------------------

end module triple_Y
!******
