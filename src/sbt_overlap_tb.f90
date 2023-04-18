!****h* FHI-aims/sbt_overlap_tb
!  NAME
!    sbt_overlap_tb
!  SYNOPSIS

module sbt_overlap_tb

  !  PURPOSE
  !
  !    Provide services for tight-binding like two-center integrations of
  !    NAOs.  For two given angular momenta L1, L2 and corresponding radial
  !    NAOs given on a logarithmic grid, there are only very few (min(L1,
  !    L2)+1) radial integrals needed for the overlap (or Coulomb) matrix
  !    elements.  As these integrals do not depend on the direction of Rvec
  !    but only on its norm, they can be splined.  After that, the two-center
  !    integrals can be calculated at tight-binding cost.
  !
  !    Please note that for convenience, power_bias_int is always assumed to
  !    be zero within this module.  That is, the ffk arrays must already
  !    contain the k**3 factor (or k**1 for bare Coulomb integrals), i.e., 1.5
  !    each.
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
  !****s* sbt_overlap_tb/get_radial_integral
  !  NAME
  !    get_radial_integral
  !  SYNOPSIS

  subroutine get_radial_integral(N, lnr0, lnk0, lnrange, &
  &                              max_LL, KKbars, wsave, &
  &                           n_fn1, L1, ffk1, n_fn2, L2, ffk2, &
  &                           mp_far, moment_1, radius_1, moment_2, radius_2, &
  &                           n_Rabs, max_L, Rabs, radI)

    !  PURPOSE
    !
    !    Evaluate the integrals
    !
    !      I_L(R) = \int_0^\infty dk k^2 \sqrt{2/\pi} j_L(kR) f_1(k) f_2(k)
    !
    !    where ffk1 * ffk2 == k^3 f_1(k) f_2(k) for givens Rs.
    !
    !  USES

    use bspline, only: val_bspline
    use mpi_tasks, only: check_allocation
    use sbt_overlap, only: logsbt_r2ilog, multipole_prefac
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: max_LL
    real*8, intent(IN) :: KKbars(2*N, 0:max_LL)
    real*8, intent(INOUT) :: wsave(4*N+15)
    integer, intent(IN) :: n_fn1, n_fn2
    integer, intent(IN) :: L1, L2
    real*8, intent(IN) :: ffk1(N, n_fn1), ffk2(N, n_fn2)
    logical, intent(IN) :: mp_far
    real*8, intent(IN) :: moment_1(n_fn1), moment_2(n_fn2)
    real*8, intent(IN) :: radius_1(n_fn1), radius_2(n_fn2)
    integer, intent(IN) :: n_Rabs, max_L
    real*8, intent(IN) :: Rabs(n_Rabs)
    real*8, intent(OUT) :: radI(0:max_L, n_fn1, n_fn2, n_Rabs)

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- Grid settings.
    !    o max_LL -- ubound(KKbars, 2), max(L1+L2)
    !    o KKbars, wsave -- Cache for logsbt/dfftpack
    !    o n_fn1, n_fn2 -- Number of different radial parts of the same Li
    !    o L1, L2 -- Angular momenta
    !    o ffk1, ffk2 -- (scaled) SBTs of radial parts
    !    o mp_far -- Need Coulomb far field
    !    o max_L -- ubound(int_spl, 2) >= min(L1, L2)
    !  OUTPUTS
    !    o radI -- Integral values
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: R0, mp_pre
    real*8, allocatable :: int_spl(:,:,:,:), int_onsite(:,:), i_R(:)
    real*8, allocatable :: this_radI(:)
    integer :: i_fn1, i_fn2, i_Rabs, L, i_L, Lmax
    integer :: info
    character(*), parameter :: func = 'get_radial_integral'

    allocate(int_spl(0:N+1, 0:max_L, n_fn1, n_fn2), stat=info)
    call check_allocation(info, 'int_spl', func)
    allocate(int_onsite(n_fn1, n_fn2), stat=info)
    call check_allocation(info, 'int_onsite', func)

    allocate(i_R(n_Rabs), this_radI(0:max_L), stat=info)

    Lmax = max(L1, L2)

    R0 = exp(lnr0)
    if (any(Rabs < R0)) then
       call get_radial_integral_onsite(N, lnr0, lnk0, lnrange, &
       &                         n_fn1, L1, ffk1, n_fn2, L2, ffk2, int_onsite)
    end if
    if (any(Rabs > R0 .and. Rabs < maxval(radius_1) + maxval(radius_2))) then
       call get_radial_integral_spline(N, lnr0, lnk0, lnrange, &
       &                               max_LL, KKbars, wsave, &
       &                               n_fn1, L1, ffk1, n_fn2, L2, ffk2, &
       &                               max_L, int_spl)
       do i_Rabs = 1, n_Rabs
          i_R = logsbt_r2ilog(Rabs(i_Rabs), lnr0, lnrange/N)
       end do
    end if
    mp_pre = multipole_prefac(L1, L2)

    do i_fn1 = 1, n_fn1
       do i_fn2 = 1, n_fn2
          do i_Rabs = 1, n_Rabs
             this_radI = 0.d0
             if (Rabs(i_Rabs) < R0) then
                if (L1 == L2) then
                   this_radI(0) = int_onsite(i_fn1, i_fn2)
                end if
             else if (Rabs(i_Rabs) < radius_1(i_fn1) + radius_2(i_fn2)) then
                call val_bspline(N, 1+Lmax, this_radI(0:Lmax), i_R(i_Rabs), &
                &                int_spl(:, 0:Lmax, i_fn1, i_fn2))
             else if (mp_far) then
                this_radI(min(L1,L2)) = mp_pre / Rabs(i_Rabs)**(L1+L2+1) &
                &                       * moment_1(i_fn1) * moment_2(i_fn2)
             end if
             radI(:, i_fn1, i_fn2, i_Rabs) = this_radI
          end do
       end do
    end do
    deallocate(int_spl, int_onsite, i_R)

  end subroutine get_radial_integral
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_tb/mult_angular_integral
  !  NAME
  !    mult_angular_integral
  !  SYNOPSIS

  subroutine mult_angular_integral(Rvec, prefac, &
  &                                n_fn1, L1, max_n_fn1, max_L1, &
  &                                n_fn2, L2, max_n_fn2, max_L2, &
  &                                max_L, need_LL, radI, fullI)


    !  PURPOSE
    !
    !    Perform the angular part of 2-center overlap (or Coulomb)
    !    calculations.  The radial integrations are expected in radI(...).
    !
    !  USES

    use mpi_tasks, only: check_allocation
    use sbt_overlap, only: sum_triple_y_ylm_real
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: Rvec(3)
    real*8, intent(IN) :: prefac
    integer, intent(IN) :: n_fn1, n_fn2
    integer, intent(IN) :: L1, L2
    integer, intent(IN) :: max_n_fn1, max_n_fn2
    integer, intent(IN) :: max_L1, max_L2
    integer, intent(IN) :: max_L
    logical, intent(IN) :: need_LL(0:max_L)
    real*8, intent(IN) :: radI(0:max_L, max_n_fn1, max_n_fn2)
    real*8, intent(OUT) :: fullI(-max_L1:max_L1, -max_L2:max_L2, &
    &                            max_n_fn1, max_n_fn2)

    !  INPUTS
    !    o Rvec -- distance vector
    !    o prefac -- common prefactor
    !    o n_fni -- number of radial parts to consider in this L-channel
    !    o Li -- angular momentum channel
    !    o max_n_fni -- maximum value of n_fni (ubound(...))
    !    o max_Li -- maximum value of Li (ubound(...))
    !    o max_L -- min(max_L1, max_L2) == ubound(radI, 1)
    !    o need_LL -- which i_LL is actually needed.
    !                 need_LL(i_LL) specifies, whether LL==min(L1,L2)+i_LL*2
    !                 is needed.  Onsite elements only need LL=0, Coulomb
    !                 far-field elements only need LL=L1+L2.
    !    o radI -- radial integral input
    !  OUTPUTS
    !    o fullI -- overlap <ffk1(0) | ffk2(Rvec) >
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: thisfac, this_radI
    integer :: max_LL, LL, i_LL
    integer :: i_fn1, i_fn2, m1, m2
    real*8, allocatable :: cyl(:,:)
    integer :: info
    character(*), parameter :: func = 'sbt_atomic_ovlp'

    ! --- Prepare arrays

    allocate(cyl(-max_L1:max_L1, -max_L2:max_L2), stat=info)
    call check_allocation(info, 'cyl', func)

    ! --- Main work

    fullI = 0.d0
    do i_LL = 0, min(L1, L2)
       if (.not. need_LL(i_LL)) cycle
       LL = abs(L1 - L2) + 2*i_LL

       ! Prefactor
       thisfac = prefac * (-1)**((L1-L2-LL)/2)

       ! Angular integration [unfortunatley entails temp array for cyl]
       ! Angular integration
       call sum_triple_Y_YLM_real(L1, L2, LL, Rvec, max_L1, max_L2, cyl)

       ! Put into array:  "ovlp = ovlp + thisfac * radI * cyl"
       do i_fn2 = 1, n_fn2
          do i_fn1 = 1, n_fn1
             this_radI = radI(i_LL, i_fn1, i_fn2)
             if (abs(this_radI) /= 0.d0) then
                do m2 = -L2, L2
                   do m1 = -L1, L1
                      fullI(m1, m2, i_fn1, i_fn2) &
                      & = fullI(m1, m2, i_fn1, i_fn2) &
                      &   + thisfac * this_radI * cyl(m1, m2)
                   end do
                end do
             end if
          end do
       end do
    end do
    deallocate(cyl)

  end subroutine mult_angular_integral
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_tb/fill_overlaps
  !  NAME
  !    fill_overlaps
  !  SYNOPSIS

  subroutine fill_overlaps(n_fn1, L1, max_n_fn1, max_L1, &
  &                        n_fn2, L2, max_n_fn2, max_L2, &
  &                        fnL1_to_row, fnL2_to_col, loc_ovlp, glob_ovlp)

    !  PURPOSE
    !
    !    Fill the overlap of two radial parts into the global array indexed
    !    by fnLi_to_rowcol.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_fn1, n_fn2
    integer, intent(IN) :: L1, L2
    integer, intent(IN) :: max_n_fn1, max_n_fn2
    integer, intent(IN) :: max_L1, max_L2
    integer, intent(IN) :: fnL1_to_row(-max_L1:max_L1, max_n_fn1)
    integer, intent(IN) :: fnL2_to_col(-max_L2:max_L2, max_n_fn2)
    real*8, intent(IN) :: loc_ovlp(-max_L1:max_L1, -max_L2:max_L2, &
    &                              max_n_fn1, max_n_fn2)
    real*8, intent(INOUT) :: glob_ovlp(:,:)

    !  INPUTS
    !    o Li -- Angular momenta
    !    o max_Li -- Array dimension, >= Li
    !    o fnLi_to_rowcol -- Indexing array
    !    o loc_ovlp -- Overlap matrix for given radial parts
    !    o glob_ovlp -- Most entries are untouched.
    !  OUTPUTS
    !    o glob_ovlp -- Elements indexed by fnLi_to_rowcol are set to loc_ovlp
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_fn1, i_fn2, m1, m2, i_col, i_row
    character(*), parameter :: func = 'fill_overlaps'

    do i_fn2 = 1, n_fn2
       do m2 = -L2, L2
          i_col = fnL2_to_col(m2, i_fn2)
          if (i_col <= 0) cycle
          do i_fn1 = 1, n_fn1
             do m1 = -L1, L1
                i_row = fnL1_to_row(m1, i_fn1)
                if (i_row <= 0) cycle
                glob_ovlp(i_row, i_col) = loc_ovlp(m1, m2, i_fn1, i_fn2)
             end do
          end do
       end do
    end do
    
  end subroutine fill_overlaps
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_tb/get_radial_integral_spline
  !  NAME
  !    get_radial_integral_spline
  !  SYNOPSIS

  subroutine get_radial_integral_spline(N, lnr0, lnk0, lnrange, &
  &                                     max_LL, KKbars, wsave, &
  &                                     n_fn1, L1, ffk1, n_fn2, L2, ffk2, &
  &                                     max_L, int_spl)

    !  PURPOSE
    !
    !    Evaluate the integrals
    !
    !      I_L(R) = \int_0^\infty dk k^2 \sqrt{2/\pi} j_L(kR) f_1(k) f_2(k)
    !
    !    where ffk1 * ffk2 == k^3 f_1(k) f_2(k) for all R on a log-grid and
    !    spline the result for L = abs(L1-L2), L1+L2, 2.  Please use
    !    get_radial_integral_onsite() for the onsite integral (L==0 only).
    !
    !  USES

    use bspline, only: cubic_bspline_notaknot
    use sbt_overlap, only: logsbt_apply_vanilla_kernel
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: max_LL
    real*8, intent(IN) :: KKbars(2*N, 0:max_LL)
    real*8, intent(INOUT) :: wsave(4*N+15)
    integer, intent(IN) :: n_fn1, n_fn2
    integer, intent(IN) :: L1, L2
    real*8, intent(IN) :: ffk1(N, n_fn1), ffk2(N, n_fn2)
    integer, intent(IN) :: max_L
    real*8, intent(OUT) :: int_spl(0:N+1, 0:max_L, n_fn1, n_fn2)

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- Grid settings.
    !    o max_LL -- ubound(KKbars, 2), max(L1+L2)
    !    o KKbars, wsave -- Cache for logsbt/dfftpack
    !    o n_fn1, n_fn2 -- Number of different radial parts of the same Li
    !    o L1, L2 -- Angular momenta
    !    o ffk1, ffk2 -- (scaled) SBTs of radial parts
    !    o max_L -- ubound(int_spl, 2) >= max(L1, L2)
    !  OUTPUTS
    !    o int_spl -- Splines of integrals
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: Nt
    real*8 :: lnt0, lntrange
    real*8, allocatable :: ffk_prod(:)
    real*8 :: logft_aliasing
    integer :: i_L, L, i_fn1, i_fn2
    integer :: info
    character(*), parameter :: func = 'get_radial_integral_spline'

    allocate(ffk_prod(2*N), stat=info)
    call check_allocation(info, 'ffk_prod', func)
    do i_L = 0, min(L1, L2)
       L = abs(L1-L2) + 2*i_L
       do i_fn1 = 1, n_fn1
          do i_fn2 = 1, n_fn2
             ffk_prod(1:N) = ffk1(:, i_fn1) * ffk2(:, i_fn2)
             ffk_prod(N+1:) = 0.d0
             call logsbt_apply_vanilla_kernel(2*N, ffk_prod, KKbars(:,L), &
             &                                wsave, logft_aliasing)

             call cubic_bspline_notaknot(N, 1, ffk_prod(1:N), &
             &                           int_spl(:, i_L, i_fn1, i_fn2))
          end do
       end do
    end do
    deallocate(ffk_prod)

  end subroutine get_radial_integral_spline
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_tb/get_radial_integral_onsite
  !  NAME
  !    get_radial_integral_onsite
  !  SYNOPSIS

  subroutine get_radial_integral_onsite(N, lnr0, lnk0, lnrange, &
  &                               n_fn1, L1, ffk1, n_fn2, L2, ffk2, int_onsite)

    !  PURPOSE
    !
    !    Evaluate the integrals
    !
    !      I_L(R=0) = \int_0^\infty dk k^2 \sqrt{2/\pi} f_1(k) f_2(k)
    !
    !    where ffk1 * ffk2 == k^3 f_1(k) f_2(k).  This corresponds to the L=0
    !    integrals.  Please use get_radial_integral_spline() for general R.
    !
    !  USES

    use constants, only: pi
    use bspline, only: val_bspline
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: n_fn1, n_fn2
    integer, intent(IN) :: L1, L2
    real*8, intent(IN) :: ffk1(N, n_fn1), ffk2(N, n_fn2)
    real*8, intent(OUT) :: int_onsite(n_fn1, n_fn2)

    !  INPUTS
    !    o N, lnr0, lnk0, lnrange -- Grid settings.
    !    o n_fn1, n_fn2 -- Number of different radial parts of the same Li
    !    o L1, L2 -- Angular momenta
    !    o ffk1, ffk2 -- (scaled) SBTs of radial parts
    !  OUTPUTS
    !    o int_spl -- Splines of integrals
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: fac
    integer :: info
    character(*), parameter :: func = 'get_radial_integral_onsite'

    if (L1 /= L2) then
       int_onsite = 0.d0
    else
       fac = sqrt(2.d0/pi) * lnrange / N
       call dgemm('T', 'N', n_fn1, n_fn2, N, &
       &          fac,  ffk1, N, &
       &                ffk2, N, &
       &          0.d0, int_onsite, n_fn1)
    end if

  end subroutine get_radial_integral_onsite
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_atomic_ovlp_tb
  !  NAME
  !    sbt_atomic_ovlp_tb
  !  SYNOPSIS

  subroutine sbt_atomic_ovlp_tb(N, lnr0, lnk0, lnrange, Rvec, &
  & prefac, has_mp_far, &
  & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
  & n_basfn1, ffk1, radius_1, moment_1, n_fnL1, fnL1_to_basfn1, fnL1_to_row, &
  & n_basfn2, ffk2, radius_2, moment_2, n_fnL2, fnL2_to_basfn2, fnL2_to_col, &
  & ovlp, times)

    !  PURPOSE
    !
    !    Calculate overlap or Coulomb interaction for a given pair of sets of
    !    atomic functions.  Make sure that the sum of the power biases of ffk1
    !    and ffk2 equals 3. in the case of sbt_type==OVLP_TYPE_OVERLAP and
    !    1. in the case of sbt_type==OVLP_TYPE_COULOMB.
    !
    !    If you are shocked by the amount and complexety of the arguments, you
    !    are in principle right; sorry for that.  But keep in mind that e.g.
    !    "n_basfn1", "radius_1", and "moment_1" can be taken directly from
    !    prodbas.f90 (n_basbas_fns, charge_radius_basbas_fn,
    !    multipole_moment_basbas_fn), and "ffk1" contains just the SBTs of
    !    basbas_wave_spl() (in the same order).  Both "n_fnL1" and
    !    "fnL1_to_basfn1" can be obtained by get_fnL_to_basfn() from
    !    basbas_atom and its friends.  The indexing array "fnL1_to_row" is a
    !    little bit tricky, but if the indexing can be chosen at will, calling
    !    get_default_fnL_to_rowcol() will do.
    !
    !    This procedure does not look into prodbas.f90 because it is also used
    !    for other things than product basis functions.  It takes some
    !    arguments sorted by L-channels because pairs of L-channels should be
    !    done in one sweep for efficiency.  And the indexing arrays
    !    fnL_to_rowcol allow for correct placement of the results irrespective
    !    of the storage scheme.
    !
    !    This subroutine is a drop-in replacement for sbt_atomic_ovlp()
    !    showing the the infrastructure of this module actually works.
    !
    !  USES

    use constants, only: pi
    use mpi_tasks, only: aims_stop, check_allocation
    use sbt_overlap, only: check_need_n_radius, SBT_N_TIMES, &
        SBT_TIME_WHOLE_PAIRS, SBT_TIME_ANGUL, SBT_TIME_DGEMM, logsbt_kernel
    use timing_core, only: start_timer, stop_timer
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    real*8, intent(IN) :: Rvec(3)
    real*8, intent(IN) :: prefac
    logical, intent(IN) :: has_mp_far
    integer, intent(IN) :: max_L1, max_L2
    integer, intent(IN) :: max_n_fnL1, max_n_fnL2
    integer, intent(IN) :: n_basfn1, n_basfn2
    real*8, intent(IN) :: ffk1(N, n_basfn1), ffk2(N, n_basfn2)
    real*8, intent(IN) :: radius_1(n_basfn1), radius_2(n_basfn2)
    real*8, intent(IN) :: moment_1(n_basfn1), moment_2(n_basfn2)
    integer, intent(IN) :: n_fnL1(0:max_L1), n_fnL2(0:max_L2)
    integer, intent(IN) :: fnL1_to_basfn1(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL2_to_basfn2(max_n_fnL2, 0:max_L2)
    integer, intent(IN) :: fnL1_to_row(-max_L1:max_L1, max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL2_to_col(-max_L2:max_L2, max_n_fnL2, 0:max_L2)
    real*8, intent(INOUT) :: ovlp(:,:)
    real*8, intent(INOUT) :: times(4, SBT_N_TIMES)

    !  INPUTS
    !    o N -- Number of 1D grid points of radial parts
    !    o lnr0, lnk0, lnrange -- Loggrid spans [exp(lnk0), exp(lnk0+lnrange)[
    !    o Rvec -- distance vector
    !    o prefac -- prefactor
    !    o has_mp_far -- Defines far field behavior (R > radius_1+radias_2):
    !            .true.  -> use moment_1*moment_2/r**(LL+1)-like far field
    !            .false. -> use 0., no far field
    !    o max_Li -- maximum value of L
    !    o max_n_fnLi -- maximum number of radial parts per L-channel
    !    o n_basfni -- total number of radial parts passed
    !    o n_fnLi -- number of radial parts to consider per L-channel
    !    o fnLi_to_basfni -- (i_fnLi, Li) -> i_basfn
    !    o fnLi_to_rowcol -- (Mi, i_fnLi, Li) -> i_rowcol [in ovlp]
    !    o ovlp -- do not touch irrelevant entries
    !  OUTPUTS
    !    o ovlp -- overlap <ffk1(0) | ffk2(Rvec) >
    !              ! Scale ffki accordingly to get the Coulomb interaction.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Rabs, lnR, acc, mp_pre, use_prefac, lnt0, lntrange
    logical :: onsite
    integer :: max_L, max_LL, i_LL, LL, m1, m2, Nt
    integer :: i_fnL1, i_fnL2, L1, L2
    integer :: i_basfn1, i_basfn2, i
    real*8, allocatable :: wsave(:), KKbars(:,:)
    real*8, allocatable :: radI(:,:,:,:,:)
    logical, allocatable :: need1(:,:), need2(:,:)
    real*8 :: this_radI, this_int
    real*8 :: max_rad_1, max_rad_2
    integer :: info
    character(*), parameter :: func = 'sbt_atomic_ovlp_tb'

    call start_timer(times(:, SBT_TIME_WHOLE_PAIRS))
    if (any(n_fnL1 > max_n_fnL1)) call aims_stop('n_fnL1 > max_n_fnL1', func)
    if (any(n_fnL2 > max_n_fnL2)) call aims_stop('n_fnL2 > max_n_fnL2', func)

    ! --- Check needs

    allocate(need1(max_n_fnL1, 0:max_L1), need2(max_n_fnL2, 0:max_L2),stat=info)
    call check_allocation(info, 'needi', func)
    call check_need_n_radius(max_L1, max_n_fnL1, n_fnL1, fnL1_to_basfn1, &
    &                        n_basfn1, radius_1, fnL1_to_row, need1, max_rad_1)
    call check_need_n_radius(max_L2, max_n_fnL2, n_fnL2, fnL2_to_basfn2, &
    &                        n_basfn2, radius_2, fnL2_to_col, need2, max_rad_2)
    if (.not. (any(need1) .and. any(need2))) then
       deallocate(need1, need2)
       return
    end if
    ! Just make sure that general checks do not say "far field" where
    ! there is actually a near field at the border.
    max_rad_1 = max_rad_1 + 1d-10
    max_rad_2 = max_rad_2 + 1d-10

    use_prefac = prefac * (2.d0*pi)**1.5d0

    ! --- General quantities

    max_LL = max_L1 + max_L2
    max_L = min(max_L1, max_L2)
    Rabs = sqrt(dot_product(Rvec, Rvec))
    onsite = (Rabs < exp(lnr0))
    ! lnr0 = sbt_lnt0_spl - lnk0 is the smallest possible value for which we
    ! have the kernel splined.  It will be of the order of 1-d10.  Everything
    ! even shorter can savely be assumed to be onsite.

    ! --- Get kernels

    Nt = 2*N
    lnt0 = lnk0 + lnr0
    lntrange = 2*lnrange
    allocate(KKbars(Nt, 0:max_LL), wsave(2*Nt+15), stat=info)
    call check_allocation(info, 'KKbars, wsave', func)
    call dffti(Nt, wsave)
    do LL = 0, max_LL
       call logsbt_kernel(Nt, KKbars(:, LL), lnt0, lntrange, LL, 0.d0, wsave)
    end do

    ! --- Prepare other arrays

    allocate(radI(0:max_L, max_n_fnL1, 0:max_L1, max_n_fnL2, 0:max_L2), &
    &        stat=info)
    call check_allocation(info, 'radI', func)
    radI = 0.d0

    ! --- Main loop

    L2_LOOP: do L2 = 0, max_L2

       L1_LOOP: do L1 = 0, max_L1

          ! Do not exploit (onsite .and. L1 /= L2) because ovlp must be reset.

          LL_LOOP: do i_LL = 0, min(L1, L2)
             LL = abs(L1 - L2) + 2*i_LL

             if (onsite .and. LL /= 0) cycle

             call start_timer(times(:, SBT_TIME_DGEMM))
             ! Radial integration
             do i_fnL2 = 1, n_fnL2(L2)
                if (.not. need2(i_fnL2, L2)) cycle
                i_basfn2 = fnL2_to_basfn2(i_fnL2, L2)
                do i_fnL1 = 1, n_fnL1(L1)
                   if (.not. need1(i_fnL1, L1)) cycle
                   i_basfn1 = fnL1_to_basfn1(i_fnL1, L1)
                   call get_radial_integral(N, lnr0, lnk0, lnrange, &
                   & max_LL, KKbars, wsave, &
                   & 1, L1, ffk1(:, i_basfn1), 1, L2, ffk2(:, i_basfn2), &
                   & has_mp_far, moment_1(i_basfn1), radius_1(i_basfn1), &
                   &             moment_2(i_basfn2), radius_2(i_basfn2), &
                   & 1, max_L, (/Rabs/), radI(:, i_fnL1, L1, i_fnL2, L2))
                end do
             end do
             call stop_timer(times(:, SBT_TIME_DGEMM), .true.)

          end do LL_LOOP

       end do L1_LOOP
    end do L2_LOOP
    call stop_timer(times(:, SBT_TIME_WHOLE_PAIRS), .true.)

    call start_timer(times(:, SBT_TIME_ANGUL))
    call sbt_atomic_ovlp_tb_angular(Rvec, use_prefac, &
    & max_L1, max_L2, max_L, max_n_fnL1, max_n_fnL2, &
    & n_fnL1, fnL1_to_row, &
    & n_fnL2, fnL2_to_col, &
    & radI, ovlp)
    call stop_timer(times(:, SBT_TIME_ANGUL), .true.)

    deallocate(radI, need1, need2, KKbars, wsave)

  end subroutine sbt_atomic_ovlp_tb
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_atomic_ovlp_tb_angular
  !  NAME
  !    sbt_atomic_ovlp_tb_angular
  !  SYNOPSIS

  subroutine sbt_atomic_ovlp_tb_angular(Rvec, prefac, &
  & max_L1, max_L2, max_L, max_n_fnL1, max_n_fnL2, &
  & n_fnL1, fnL1_to_row, n_fnL2, fnL2_to_col, radI, ovlp)

    !  PURPOSE
    !
    !    Perform the angular part of 2-center overlap (or Coulomb)
    !    calculations.  The radial integrations are expected in radI(...).
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: Rvec(3)
    real*8, intent(IN) :: prefac
    integer, intent(IN) :: max_L1, max_L2, max_L
    integer, intent(IN) :: max_n_fnL1, max_n_fnL2
    integer, intent(IN) :: n_fnL1(0:max_L1), n_fnL2(0:max_L2)
    integer, intent(IN) :: fnL1_to_row(-max_L1:max_L1, max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL2_to_col(-max_L2:max_L2, max_n_fnL2, 0:max_L2)
    real*8, intent(IN) :: radI(0:max_L, max_n_fnL1, 0:max_L1, &
    &                                       max_n_fnL2, 0:max_L2)
    real*8, intent(INOUT) :: ovlp(:,:)

    !  INPUTS
    !    o Rvec -- distance vector
    !    o prefac -- common prefactor
    !    o max_Li -- maximum value of L
    !    o max_L -- min(max_L1, max_L2) == ubound(radI, 1)
    !    o max_n_fnLi -- maximum number of radial parts per L-channel
    !    o n_fnLi -- number of radial parts to consider per L-channel
    !    o fnLi_to_rowcol -- (Mi, i_fnLi, Li) -> i_rowcol [in ovlp]
    !    o radI -- radial integral input
    !    o ovlp -- do not touch irrelevant entries
    !  OUTPUTS
    !    o ovlp -- overlap <ffk1(0) | ffk2(Rvec) >
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_fn1, n_fn2, L1, L2
    real*8, allocatable :: this_ovlp(:,:,:,:)
    integer :: info
    logical :: need_LL(0:max_L)
    character(*), parameter :: func = 'sbt_atomic_ovlp_tb'

    if (any(n_fnL1 > max_n_fnL1)) call aims_stop('n_fnL1 > max_n_fnL1', func)
    if (any(n_fnL2 > max_n_fnL2)) call aims_stop('n_fnL2 > max_n_fnL2', func)

    ! --- Prepare arrays

    allocate(this_ovlp(-max_L1:max_L1, -max_L2:max_L2, &
    &                  max_n_fnL1, max_n_fnL2), stat=info)
    call check_allocation(info, 'this_ovlp', func)
    need_LL = .true.  ! Not as smart as it should be...

    ! --- Main loop

    do L2 = 0, max_L2
       n_fn2 = n_fnL2(L2)
       do L1 = 0, max_L1
          n_fn1 = n_fnL1(L1)

          call mult_angular_integral(Rvec, prefac, &
          &                          n_fn1, L1, max_n_fnL1, max_L1, &
          &                          n_fn2, L2, max_n_fnL2, max_L2, &
          &                          max_L, need_LL, radI(:,:, L1,:, L2), &
          &                          this_ovlp)

          call fill_overlaps(n_fn1, L1, max_n_fnL1, max_L1, &
          &                  n_fn2, L2, max_n_fnL2, max_L2, &
          &                  fnL1_to_row(:,:, L1), &
          &                  fnL2_to_col(:,:, L2), &
          &                  this_ovlp, ovlp)

       end do
    end do

    deallocate(this_ovlp)

  end subroutine sbt_atomic_ovlp_tb_angular
  !******
end module sbt_overlap_tb
!******
