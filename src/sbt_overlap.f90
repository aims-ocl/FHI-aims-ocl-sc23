!****h* FHI-aims/sbt_overlap
!  NAME
!    sbt_overlap
!  SYNOPSIS

module sbt_overlap

  !  PURPOSE
  !
  !    Provide infrastructure to calculate overlap integrals using
  !    spherical Bessel transforms.
  !
  !  USES

  ! WPH:  I've left the global use statements for some modules intact here,
  !       as this module may function as a wrapper module around them.
  use logsbt
  use triple_Y
  use logsbt_fast_kernel
  implicit none

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

  integer, parameter :: SBT_TIME_P2 = 1
  integer, parameter :: SBT_TIME_DGEMM = 2
  integer, parameter :: SBT_TIME_ANGUL = 3
  integer, parameter :: SBT_TIME_KERNEL = 4
  integer, parameter :: SBT_TIME_WHOLE_PAIRS = 5
  integer, parameter :: SBT_TIME_TRIPLE_MULT = 6
  integer, parameter :: SBT_TIME_TRIPLE_SBT = 7
  integer, parameter :: SBT_TIME_TRIPLE_COMB = 8
  integer, parameter :: SBT_TIME_WHOLE_TRIPLES = 9
  integer, parameter :: SBT_N_TIMES = 9

  real*8, allocatable :: ffk1(:,:)
  integer:: index_old=-1

contains

  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_atomic_ovlp
  !  NAME
  !    sbt_atomic_ovlp
  !  SYNOPSIS

  subroutine sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, Rvec, prefac, has_mp_far,&
  & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
  & n_basfn1, ffk1, radius_1, moment_1, n_fnL1, fnL1_to_basfn1, fnL1_to_row, &
  & n_basfn2, ffk2, radius_2, moment_2, n_fnL2, fnL2_to_basfn2, fnL2_to_col, &
  & ovlp, times, d_ovlp)

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
    !  USES

    use constants, only: pi
    use mpi_tasks, only: aims_stop, check_allocation
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
    real*8, intent(INOUT), optional :: d_ovlp(:,:,:)

    !  FIXME
    !
    !    o Could also take ovlpT, fn1_to_col, fn2_to_row; use spglib_symmetry.
    !    o It would be entertaining to fill ffkL2j(:,:) with exact duals:
    !      + Only the lowest N logFT components of ffkL2j are needed (rest is
    !        orthogonal to ffk1).
    !      + To those, only the lowest 2*N logFT components of the kernel
    !        contribute (convolution with N ffk2 logFT components).
    !      + The product of ffk2 with a 2*N kernel can be calculated in
    !        product space as a simple product on a 3*N grid.
    !
    !      Therefore, do the following:
    !      + LogFT transform ffk2, add zeros, transform back to 3*N grid.
    !        [Can be done outside this procedure.]
    !      + Construct 2*N logFT component kernel, add zeros,
    !        transform to 3*N grid.
    !        [Can be splined(?)]
    !      + Multiply kernel to logFT.
    !      + Transform to logFT space and throw away high components.
    !      + Backtransform to N grid.
    !      + Do dgemm part.
    !
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
    !    o ffki -- Should contain f~(k)*k**1.5.
    !              In case of Coulomb matrix elements, this effectively means
    !              to use f~(k)*k**0.5 and prefac=4*pi.
    !    o radius_i -- Radius of function ffr (or charge if (has_mp_far))
    !    o moment_i -- Multipole moment to be used if (has_mp_far)
    !    o n_fnLi -- number of radial parts to consider per L-channel
    !    o fnLi_to_basfni -- (i_fnLi, Li) -> i_basfn
    !    o fnLi_to_rowcol -- (Mi, i_fnLi, Li) -> i_rowcol [in ovlp]
    !    o ovlp -- irrelevant entries are not touched within this procedure
    !  OUTPUTS
    !    o ovlp -- overlap <ffk1(0) | ffk2(Rvec) >
    !              ! Scale ffki accordingly to get the Coulomb interaction.
    !    o d_ovlp -- derivative of ovlp with respect to Rvec
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Rabs, lnR, thisfac, acc, mp_pre
    logical :: onsite
    integer :: max_LL, LL
    integer :: i_fnL1, i_fnL2, L1, L2, m1, m2
    integer :: i_row, i_col, i_basfn1, i_basfn2, i
    real*8, allocatable :: KK(:,:), ffkL2j(:,:,:), d_KK(:,:), d_ffkL2j(:,:,:)
    real*8, allocatable :: radI(:,:), cyl(:,:), this_ovlp(:,:,:,:)
    real*8, allocatable :: d_radI(:,:), d_cyl(:,:,:), d_this_ovlp(:,:,:,:,:)
    logical, allocatable :: need1(:,:), need2(:,:)
    real*8 :: max_rad_1, max_rad_2
    integer :: info
    logical :: calc_deriv
    character(*), parameter :: func = 'sbt_atomic_ovlp'

    if (sbt_n_spl_ref <= 0) then
       call aims_stop('logsbt_fast_kernel not initialized', func)
    end if

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

    ! --- General quantities

    max_LL = max_L1 + max_L2
    Rabs = sqrt(dot_product(Rvec, Rvec))
    onsite = (Rabs < exp(lnr0))
    ! lnr0 is the smallest possible value for which we
    ! have the kernel splined.  It will be of the order of 1-d10.  Everything
    ! even shorter can savely be assumed to be onsite.

    calc_deriv = present(d_ovlp) .and. .not.onsite

    ! --- Get kernels for given |R|

    call start_timer(times(:, SBT_TIME_KERNEL))
    if (Rabs < max_rad_1 + max_rad_2) then
       allocate(ffkL2j(N, max_n_fnL2, 0:max_LL), stat=info)
       call check_allocation(info, 'ffkL2j', func)
       allocate(KK(N, 0:max_LL), stat=info)
       call check_allocation(info, 'KK', func)
       if(calc_deriv) then
          allocate(d_ffkL2j(N, max_n_fnL2, 0:max_LL), stat=info)
          call check_allocation(info, 'd_ffkL2j', func)
          allocate(d_KK(N, 0:max_LL), stat=info)
          call check_allocation(info, 'd_KK', func)
       endif
       if (onsite) then
          KK(:, 0) = sqrt(2.d0/pi)       ! sqrt(2/pi) j_0(0) = sqrt(2/pi)
          KK(:, 1:) = 0.d0
       else
          lnR = log(Rabs)
          call get_fast_kernel(N, 0, max_LL, lnR+lnk0, lnrange, &
          &                    0.d0, .true., KK)
          if(calc_deriv) then
             call get_fast_kernel_der(N, 0, max_LL, lnR+lnk0, lnrange, &
             &                        0.d0, .true., d_KK)
          endif
       end if
    end if
    call stop_timer(times(:, SBT_TIME_KERNEL), .true.)

    ! --- Prepare other arrays

    allocate(radI(max_n_fnL1, max_n_fnL2), stat=info)
    call check_allocation(info, 'radI', func)
    allocate(cyl(-max_L1:max_L1, -max_L2:max_L2), stat=info)
    call check_allocation(info, 'cyl', func)
    allocate(this_ovlp(-max_L1:max_L1, -max_L2:max_L2, max_n_fnL1, max_n_fnL2),stat=info)
    call check_allocation(info, 'this_ovlp', func)
    if(calc_deriv) then
       allocate(d_radI(max_n_fnL1, max_n_fnL2), stat=info)
       call check_allocation(info, 'd_radI', func)
       allocate(d_this_ovlp(3, -max_L1:max_L1, -max_L2:max_L2, max_n_fnL1, max_n_fnL2),stat=info)
       call check_allocation(info, 'd_this_ovlp', func)
    endif
    ! d_cyl must always be allocated since it is used as argument to fast_sum_triple_Y_YLM_real
    allocate(d_cyl(-max_L1:max_L1, -max_L2:max_L2, 3), stat=info)
    call check_allocation(info, 'd_cyl', func)

    ! Precalculate  the spherical harmonics for fast_sum_triple_Y_YLM_real
    call calc_YrLM_for_triple_Y(Rvec, max(max_L1,max_L2))

    ! --- Main loop

    L2_LOOP: do L2 = 0, max_L2

       ! prepare for given L2
       call start_timer(times(:, SBT_TIME_P2))
       do LL = 0, max_LL
          if (onsite .and. LL /= 0) cycle
          do i_fnL2 = 1, n_fnL2(L2)
             if (.not. need2(i_fnL2, L2)) cycle
             i_basfn2 = fnL2_to_basfn2(i_fnL2, L2)
             if (Rabs < max_rad_1 + radius_2(i_basfn2)) then
                ffkL2j(:, i_fnL2, LL) = ffk2(:, i_basfn2) * KK(:, LL)
                if(calc_deriv) &
                   d_ffkL2j(:, i_fnL2, LL) = ffk2(:, i_basfn2) * d_KK(:, LL)
             end if
          end do
       end do
       call stop_timer(times(:, SBT_TIME_P2), .true.)

       L1_LOOP: do L1 = 0, max_L1

          this_ovlp(-L1:L1, -L2:L2, 1:n_fnL1(L1), 1:n_fnL2(L2)) = 0.d0
          if(calc_deriv) &
             d_this_ovlp(:, -L1:L1, -L2:L2, 1:n_fnL1(L1), 1:n_fnL2(L2)) = 0.d0

          ! Do not exploit (onsite .and. L1 /= L2) because ovlp must be reset.

          LL_LOOP: do LL = abs(L1 - L2), L1 + L2, 2

             if (onsite .and. LL /= 0) cycle

             ! Prefactor
             thisfac = prefac * (2*pi)**1.5d0 * (-1)**((L1-L2-LL)/2)

             call start_timer(times(:, SBT_TIME_DGEMM))
             ! Radial integration
             do i_fnL2 = 1, n_fnL2(L2)
                if (.not. need2(i_fnL2, L2)) cycle
                i_basfn2 = fnL2_to_basfn2(i_fnL2, L2)
                do i_fnL1 = 1, n_fnL1(L1)
                   if (.not. need1(i_fnL1, L1)) cycle
                   i_basfn1 = fnL1_to_basfn1(i_fnL1, L1)
                   if (Rabs < radius_1(i_basfn1) + radius_2(i_basfn2)) then
                      acc = 0.d0
                      do i = 1, N
                         acc = acc + ffk1(i, i_basfn1) * ffkL2j(i, i_fnL2, LL)
                      end do
                      radI(i_fnL1, i_fnL2) = acc * lnrange/N
                      if(calc_deriv) then
                         acc = 0.d0
                         do i = 1, N
                            acc = acc + ffk1(i, i_basfn1) * d_ffkL2j(i, i_fnL2, LL)
                         end do
                         ! d_radI := d radI / d Rabs = (d radI / d lnR) * (d lnR / d Rabs)
                         d_radI(i_fnL1, i_fnL2) = acc * lnrange/N / Rabs
                      endif
                   else if (has_mp_far .and. LL == L1 + L2) then
                      mp_pre = multipole_prefac(L1, L2)
                      radI(i_fnL1, i_fnL2) = 1.d0 / Rabs**(LL+1) * mp_pre &
                      &              * moment_1(i_basfn1) * moment_2(i_basfn2)
                      if(calc_deriv) then
                         d_radI(i_fnL1, i_fnL2) = -(LL+1) * radI(i_fnL1, i_fnL2) / Rabs
                      endif
                   else
                      radI(i_fnL1, i_fnL2) = 0.d0
                      if(calc_deriv) d_radI(i_fnL1, i_fnL2) = 0.d0
                   end if
                end do
             end do
!!$             call dgemm('T', 'N', n_fnL1, n_fnL2, N, &
!!$             &          1.d0, ffkL1, N, &
!!$             &                ffkL2j(:,:, LL), N, &
!!$             &          0.d0, radI, max_fnL1)
             call stop_timer(times(:, SBT_TIME_DGEMM), .true.)

             call start_timer(times(:, SBT_TIME_ANGUL))
             ! Angular integration
             call fast_sum_triple_Y_YLM_real(L1, L2, LL, max_L1, max_L2, cyl, d_cyl, calc_deriv)
             call stop_timer(times(:, SBT_TIME_ANGUL), .true.)

             ! Put into array:  "ovlp = ovlp + thisfac * radI * cyl"
             do i_fnL2 = 1, n_fnL2(L2)
                if (.not. need2(i_fnL2, L2)) cycle
                do i_fnL1 = 1, n_fnL1(L1)
                   if (.not. need1(i_fnL1, L1)) cycle
                   do m2 = -L2, L2
                      do m1 = -L1, L1
                         this_ovlp(m1, m2, i_fnL1, i_fnL2) = &
                         & this_ovlp(m1, m2, i_fnL1, i_fnL2) &
                         & + thisfac * radI(i_fnL1, i_fnL2) * cyl(m1, m2)
                         if(calc_deriv) then
                            d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 1) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(1)/Rabs * cyl(m1, m2)
                            d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 2) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(2)/Rabs * cyl(m1, m2)
                            d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 3) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(3)/Rabs * cyl(m1, m2)
                         endif
                      end do
                   end do
                end do
             end do
          end do LL_LOOP

          ! Put into output array
          do i_fnL2 = 1, n_fnL2(L2)
             do m2 = -L2, L2
                i_col = fnL2_to_col(m2, i_fnL2, L2)
                if (i_col <= 0) cycle
                do i_fnL1 = 1, n_fnL1(L1)
                   do m1 = -L1, L1
                      i_row = fnL1_to_row(m1, i_fnL1, L1)
                      if (i_row <= 0) cycle
                      ovlp(i_row, i_col) = this_ovlp(m1, m2, i_fnL1, i_fnL2)
                      if(present(d_ovlp)) then
                         if(calc_deriv) then
                            d_ovlp(1, i_row, i_col) = d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2)
                            d_ovlp(2, i_row, i_col) = d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2)
                            d_ovlp(3, i_row, i_col) = d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2)
                         else
                            d_ovlp(1, i_row, i_col) = 0.
                            d_ovlp(2, i_row, i_col) = 0.
                            d_ovlp(3, i_row, i_col) = 0.
                         endif
                      endif
                   end do
                end do
             end do
          end do

       end do L1_LOOP
    end do L2_LOOP
    call stop_timer(times(:, SBT_TIME_WHOLE_PAIRS), .true.)

    if (allocated(KK)) deallocate(KK)
    if (allocated(ffkL2j)) deallocate(ffkL2j)
    if (allocated(d_KK)) deallocate(d_KK)
    if (allocated(d_ffkL2j)) deallocate(d_ffkL2j)
    deallocate(radI, cyl, this_ovlp, need1, need2)
    if(calc_deriv) deallocate(d_radI, d_this_ovlp)
    deallocate(d_cyl)

  end subroutine sbt_atomic_ovlp
  subroutine my_sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, Rvec, prefac, has_mp_far,&
  & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
  & n_basfn1, ffk1, radius_1, moment_1, n_fnL1, fnL1_to_basfn1, fnL1_to_row, &
  & n_basfn2, ffk2, radius_2, moment_2, n_fnL2, fnL2_to_basfn2, fnL2_to_col, &
  & ovlp, times, d_ovlp)

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
    !  USES
    use constants, only: pi
    use my_triple_y
    use mpi_tasks, only: aims_stop, check_allocation
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
    real*8, intent(INOUT), optional :: d_ovlp(:,:,:)

    !  FIXME
    !
    !    o Could also take ovlpT, fn1_to_col, fn2_to_row; use spglib_symmetry.
    !    o It would be entertaining to fill ffkL2j(:,:) with exact duals:
    !      + Only the lowest N logFT components of ffkL2j are needed (rest is
    !        orthogonal to ffk1).
    !      + To those, only the lowest 2*N logFT components of the kernel
    !        contribute (convolution with N ffk2 logFT components).
    !      + The product of ffk2 with a 2*N kernel can be calculated in
    !        product space as a simple product on a 3*N grid.
    !
    !      Therefore, do the following:
    !      + LogFT transform ffk2, add zeros, transform back to 3*N grid.
    !        [Can be done outside this procedure.]
    !      + Construct 2*N logFT component kernel, add zeros,
    !        transform to 3*N grid.
    !        [Can be splined(?)]
    !      + Multiply kernel to logFT.
    !      + Transform to logFT space and throw away high components.
    !      + Backtransform to N grid.
    !      + Do dgemm part.
    !
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
    !    o ffki -- Should contain f~(k)*k**1.5.
    !              In case of Coulomb matrix elements, this effectively means
    !              to use f~(k)*k**0.5 and prefac=4*pi.
    !    o radius_i -- Radius of function ffr (or charge if (has_mp_far))
    !    o moment_i -- Multipole moment to be used if (has_mp_far)
    !    o n_fnLi -- number of radial parts to consider per L-channel
    !    o fnLi_to_basfni -- (i_fnLi, Li) -> i_basfn
    !    o fnLi_to_rowcol -- (Mi, i_fnLi, Li) -> i_rowcol [in ovlp]
    !    o ovlp -- irrelevant entries are not touched within this procedure
    !  OUTPUTS
    !    o ovlp -- overlap <ffk1(0) | ffk2(Rvec) >
    !              ! Scale ffki accordingly to get the Coulomb interaction.
    !    o d_ovlp -- derivative of ovlp with respect to Rvec
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Rabs, lnR, thisfac, acc, mp_pre
    logical :: onsite
    integer :: max_LL, LL
    integer :: i_fnL1, i_fnL2, L1, L2, m1, m2
    integer :: i_row, i_col, i_basfn1, i_basfn2, i
    real*8, allocatable :: KK(:,:), ffkL2j(:,:,:), d_KK(:,:), d_ffkL2j(:,:,:)
    real*8, allocatable :: radI(:,:), cyl(:,:), this_ovlp(:,:,:,:)
    real*8, allocatable :: d_radI(:,:), d_cyl(:,:,:), d_this_ovlp(:,:,:,:,:)
    logical, allocatable :: need1(:,:), need2(:,:)
    real*8 :: max_rad_1, max_rad_2
    integer :: info
    logical :: calc_deriv
    type(YrLM_t) :: YrLM
!    type(CrMs_t), allocatable :: CrMs(:,:,:)
!    type(CrMs_t):: CrMs

    integer :: max_L_YrLM_tmp
    character(*), parameter :: func = 'my_sbt_atomic_ovlp'

    if (sbt_n_spl_ref <= 0) then
       call aims_stop('logsbt_fast_kernel not initialized', func)
    end if

    max_L_YrLM_tmp = max(max(max_L1,max_L2),15)
!    call YrLM%construct(max_L_YrLM_tmp)
    call allocate_YrLM_type_arrays(YrLM, max_L_YrLM_tmp) 

!    allocate(CrMs(0:max_L_YrLM_tmp,0:max_L_YrLM_tmp,0:2*max_L_YrLM_tmp))
!    call CrMs%construct(max_L1, max_L2)

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

    ! --- General quantities

    max_LL = max_L1 + max_L2
    Rabs = sqrt(dot_product(Rvec, Rvec))
    onsite = (Rabs < exp(lnr0))
    ! lnr0 is the smallest possible value for which we
    ! have the kernel splined.  It will be of the order of 1-d10.  Everything
    ! even shorter can savely be assumed to be onsite.

    calc_deriv = present(d_ovlp) .and. .not.onsite
    ! --- Get kernels for given |R|
    ! check this
    call start_timer(times(:, SBT_TIME_KERNEL))
    if (Rabs < max_rad_1 + max_rad_2) then
       allocate(ffkL2j(N, max_n_fnL2, 0:max_LL), stat=info)
       call check_allocation(info, 'ffkL2j', func)
       allocate(KK(N, 0:max_LL), stat=info)
       call check_allocation(info, 'KK', func)
       if(calc_deriv) then
          allocate(d_ffkL2j(N, max_n_fnL2, 0:max_LL), stat=info)
          call check_allocation(info, 'd_ffkL2j', func)
          allocate(d_KK(N, 0:max_LL), stat=info)
          call check_allocation(info, 'd_KK', func)
       endif
       if (onsite) then
          KK(:, 0) = sqrt(2.d0/pi)       ! sqrt(2/pi) j_0(0) = sqrt(2/pi)
          KK(:, 1:) = 0.d0
       else
          lnR = log(Rabs)
          call get_fast_kernel(N, 0, max_LL, lnR+lnk0, lnrange, &
          &                    0.d0, .true., KK)
          if(calc_deriv) then
             call get_fast_kernel_der(N, 0, max_LL, lnR+lnk0, lnrange, &
             &                        0.d0, .true., d_KK)
          endif
       end if
    end if
    call stop_timer(times(:, SBT_TIME_KERNEL), .true.)
    ! --- Prepare other arrays

    allocate(radI(max_n_fnL1, max_n_fnL2), stat=info)
    call check_allocation(info, 'radI', func)
    allocate(cyl(-max_L1:max_L1, -max_L2:max_L2), stat=info)
    call check_allocation(info, 'cyl', func)
    allocate(this_ovlp(-max_L1:max_L1, -max_L2:max_L2, max_n_fnL1, max_n_fnL2),stat=info)
    call check_allocation(info, 'this_ovlp', func)
    if(calc_deriv) then
       allocate(d_radI(max_n_fnL1, max_n_fnL2), stat=info)
       call check_allocation(info, 'd_radI', func)
       allocate(d_this_ovlp(3, -max_L1:max_L1, -max_L2:max_L2, max_n_fnL1, max_n_fnL2),stat=info)
       call check_allocation(info, 'd_this_ovlp', func)
    endif
    ! d_cyl must always be allocated since it is used as argument to fast_sum_triple_Y_YLM_real
    allocate(d_cyl(-max_L1:max_L1, -max_L2:max_L2, 3), stat=info)
    call check_allocation(info, 'd_cyl', func)



    ! Precalculate  the spherical harmonics for fast_sum_triple_Y_YLM_real
!    call copy_calc_YrLM_for_triple_Y(Rvec, max(max_L1,max_L2))

    call my_calc_YrLM_for_triple_Y(Rvec, max(max_L1,max_L2), YrLM)
    call allocate_CrMs_s(max_L1, max_L2, max_L1+max_L2)

    ! --- Main loop
    ! check this
    do L2 = 0, max_L2
       ! prepare for given L2
       call start_timer(times(:, SBT_TIME_P2))
       do LL = 0, max_LL
          if (onsite .and. LL /= 0) cycle
          do i_fnL2 = 1, n_fnL2(L2)
             if (.not. need2(i_fnL2, L2)) cycle
             i_basfn2 = fnL2_to_basfn2(i_fnL2, L2)
             if (Rabs < max_rad_1 + radius_2(i_basfn2)) then
                ffkL2j(:, i_fnL2, LL) = ffk2(:, i_basfn2) * KK(:, LL)
                if(calc_deriv) &
                   d_ffkL2j(:, i_fnL2, LL) = ffk2(:, i_basfn2) * d_KK(:, LL)
             end if
          end do
       end do
!     enddo  
       call stop_timer(times(:, SBT_TIME_P2), .true.)

 !!   do L2 = 0, max_L2

       do L1 = 0, max_L1

          this_ovlp(-L1:L1, -L2:L2, 1:n_fnL1(L1), 1:n_fnL2(L2)) = 0.d0
          if(calc_deriv) &
             d_this_ovlp(:, -L1:L1, -L2:L2, 1:n_fnL1(L1), 1:n_fnL2(L2)) = 0.d0

          ! Do not exploit (onsite .and. L1 /= L2) because ovlp must be reset.

          do LL = abs(L1 - L2), L1 + L2, 2

             if (onsite .and. LL /= 0) cycle

             ! Prefactor
             thisfac = prefac * (2*pi)**1.5d0 * (-1)**((L1-L2-LL)/2)

             call start_timer(times(:, SBT_TIME_DGEMM))
             ! Radial integration
             do i_fnL2 = 1, n_fnL2(L2)
                if (.not. need2(i_fnL2, L2)) cycle
                i_basfn2 = fnL2_to_basfn2(i_fnL2, L2)
                do i_fnL1 = 1, n_fnL1(L1)
                   if (.not. need1(i_fnL1, L1)) cycle
                   i_basfn1 = fnL1_to_basfn1(i_fnL1, L1)
                   if (Rabs < radius_1(i_basfn1) + radius_2(i_basfn2)) then
                      acc = 0.d0
                      do i = 1, N
                         acc = acc + ffk1(i, i_basfn1) * ffkL2j(i, i_fnL2, LL)
                      end do
                      radI(i_fnL1, i_fnL2) = acc * lnrange/N
                      if(calc_deriv) then
                         acc = 0.d0
                         do i = 1, N
                            acc = acc + ffk1(i, i_basfn1) * d_ffkL2j(i, i_fnL2, LL)
                         end do
                         ! d_radI := d radI / d Rabs = (d radI / d lnR) * (d lnR / d Rabs)
                         d_radI(i_fnL1, i_fnL2) = acc * lnrange/N / Rabs
                      endif
                   else if (has_mp_far .and. LL == L1 + L2) then
                      mp_pre = multipole_prefac(L1, L2)
                      radI(i_fnL1, i_fnL2) = 1.d0 / Rabs**(LL+1) * mp_pre &
                      &              * moment_1(i_basfn1) * moment_2(i_basfn2)
                      if(calc_deriv) then
                         d_radI(i_fnL1, i_fnL2) = -(LL+1) * radI(i_fnL1, i_fnL2) / Rabs
                      endif
                   else
                      radI(i_fnL1, i_fnL2) = 0.d0
                      if(calc_deriv) d_radI(i_fnL1, i_fnL2) = 0.d0
                   end if
                end do
             end do
!!$             call dgemm('T', 'N', n_fnL1, n_fnL2, N, &
!!$             &          1.d0, ffkL1, N, &
!!$             &                ffkL2j(:,:, LL), N, &
!!$             &          0.d0, radI, max_fnL1)
             call stop_timer(times(:, SBT_TIME_DGEMM), .true.)

             call start_timer(times(:, SBT_TIME_ANGUL))
             ! Angular integration

!             call copy_fast_sum_triple_Y_YLM_real(L1, L2, LL, max_L1, max_L2, cyl, d_cyl, calc_deriv)

             call my_fast_sum_triple_Y_YLM_real(L1, L2, LL, max_L1, max_L2, cyl, d_cyl, calc_deriv, YrLM)
             call stop_timer(times(:, SBT_TIME_ANGUL), .true.)

             ! Put into array:  "ovlp = ovlp + thisfac * radI * cyl"
             do i_fnL2 = 1, n_fnL2(L2)
                if (.not. need2(i_fnL2, L2)) cycle
                do i_fnL1 = 1, n_fnL1(L1)
                   if (.not. need1(i_fnL1, L1)) cycle
                   do m2 = -L2, L2
                      do m1 = -L1, L1
                         this_ovlp(m1, m2, i_fnL1, i_fnL2) = &
                         & this_ovlp(m1, m2, i_fnL1, i_fnL2) &
                         & + thisfac * radI(i_fnL1, i_fnL2) * cyl(m1, m2)
                         if(calc_deriv) then
                            d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 1) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(1)/Rabs * cyl(m1, m2)
                            d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 2) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(2)/Rabs * cyl(m1, m2)
                            d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2) = &
                            & d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2) &
                            & + thisfac * radI(i_fnL1, i_fnL2) * d_cyl(m1, m2, 3) &
                            & + thisfac * d_radI(i_fnL1, i_fnL2) * Rvec(3)/Rabs * cyl(m1, m2)
                         endif
                      end do
                   end do
                end do
             end do
          end do ! LL_LOOP

          ! Put into output array
          do i_fnL2 = 1, n_fnL2(L2)
             do m2 = -L2, L2
                i_col = fnL2_to_col(m2, i_fnL2, L2)
                if (i_col <= 0) cycle
                do i_fnL1 = 1, n_fnL1(L1)
                   do m1 = -L1, L1
                      i_row = fnL1_to_row(m1, i_fnL1, L1)
                      if (i_row <= 0) cycle
                      ovlp(i_row, i_col) = this_ovlp(m1, m2, i_fnL1, i_fnL2)
                      if(present(d_ovlp)) then
                         if(calc_deriv) then
                            d_ovlp(1, i_row, i_col) = d_this_ovlp(1, m1, m2, i_fnL1, i_fnL2)
                            d_ovlp(2, i_row, i_col) = d_this_ovlp(2, m1, m2, i_fnL1, i_fnL2)
                            d_ovlp(3, i_row, i_col) = d_this_ovlp(3, m1, m2, i_fnL1, i_fnL2)
                         else
                            d_ovlp(1, i_row, i_col) = 0.
                            d_ovlp(2, i_row, i_col) = 0.
                            d_ovlp(3, i_row, i_col) = 0.
                         endif
                      endif
                   end do
                end do
             end do
          end do

       end do ! L1_LOOP
    end do ! L2_LOOP

    call stop_timer(times(:, SBT_TIME_WHOLE_PAIRS), .true.)

!    call YrLM%deconstruct(max_L_YrLM_tmp)
    call deallocate_YrLM_type_arrays(YrLM,max_L_YrLM_tmp)
    call deallocate_CrMs_s(max_L1, max_L2, max_L1+max_L2)
    if (allocated(KK)) deallocate(KK)
    if (allocated(ffkL2j)) deallocate(ffkL2j)
    if (allocated(d_KK)) deallocate(d_KK)
    if (allocated(d_ffkL2j)) deallocate(d_ffkL2j)
    deallocate(radI, cyl, this_ovlp, need1, need2)
    if(calc_deriv) deallocate(d_radI, d_this_ovlp)
    deallocate(d_cyl)
    return
  end subroutine my_sbt_atomic_ovlp

  
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/output_sbt_timing
  !  NAME
  !    output_sbt_timing
  !  SYNOPSIS

  subroutine output_sbt_timing(times, fmt, priority)

    !  PURPOSE
    !    Output accumulated timings of sbt_atomic_ovlp().
    !  USES

    use synchronize_mpi_basic
    use timing_core, only: output_timer
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: times(4, SBT_N_TIMES)
    character(*), intent(IN), optional :: fmt
    integer, intent(IN), optional :: priority

    !  INPUTS
    !    o times -- Accumulated timings from sbt_atomic_ovlp()
    !    o fmt -- Format to be used for initial text (default to '2X')
    !    o priority -- Output priority (default to infinitely high)
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i
    character(*), parameter :: func = 'output_sbt_timing'

    do i = 1, SBT_N_TIMES
       call sync_timing(times(3, i))
    end do
    call output_timer('Multiplication with kernel', &
    &                 times(3:4, SBT_TIME_P2), fmt, priority)
    call output_timer('Main matrix multiplication', &
    &                 times(3:4, SBT_TIME_DGEMM), fmt, priority)
    call output_timer('Analytic angular integration', &
    &                 times(3:4, SBT_TIME_ANGUL), fmt, priority)
    call output_timer('Kernel construction', &
    &                 times(3:4, SBT_TIME_KERNEL), fmt, priority)
    call output_timer('Overall 2-center overlap time', &
    &                 times(3:4, SBT_TIME_WHOLE_PAIRS), fmt, priority)
    if (any(times(3:4, SBT_TIME_WHOLE_TRIPLES) > 0.d0)) then
       call output_timer('3-center preparation multiplications', &
       &                 times(3:4, SBT_TIME_TRIPLE_MULT), fmt, priority)
       call output_timer('3-center preparation SBTs', &
       &                 times(3:4, SBT_TIME_TRIPLE_SBT), fmt, priority)
       call output_timer('3-center combinations', &
       &                 times(3:4, SBT_TIME_TRIPLE_COMB), fmt, priority)
       call output_timer('Overal 3-center time', &
       &                 times(3:4, SBT_TIME_WHOLE_TRIPLES), fmt, priority)
    end if

  end subroutine output_sbt_timing
  !******
  !----------------------------------------------------------------------------
  !****s* logsbt/sbt_import_spline
  !  NAME
  !    sbt_import_spline
  !  SYNOPSIS

  subroutine sbt_import_spline(N, ff, lnr0, lnrange, L, &
  &                            n_grid, wave_spl, r_grid_min, r_grid_inc)

    !  PURPOSE
    !
    !    Convert a log-splined radial part u(r) as used in
    !    psi(\bm r) = u(r)/r Y_{lm}(\hat r) to a radial part on a loggrid
    !    as used in this subroutine (without '/r').
    !
    !  USES

    use spline, only: val_spline
    use localorb_io, only: localorb_info
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(OUT) :: ff(N)
    real*8, intent(IN) :: lnr0, lnrange
    integer, intent(IN) :: L
    integer, intent(IN) :: n_grid
    real*8, intent(IN) :: wave_spl(4, n_grid)
    real*8, intent(IN) :: r_grid_min, r_grid_inc

    !  INPUTS
    !    o N -- number of grid points
    !    o lnr0 -- log of onset of r-space grid
    !    o lnrange -- (log) range of both r- and k-space grid
    !    o L -- angular momentum (needed only for small-r extrapolation ~r^L)
    !    o n_grid -- number of spline grid points
    !    o wave_spl -- spline coefficients
    !    o r_grid_min -- left border of spline definition
    !    o r_grid_inc -- increment factor between grid points
    !  OUTPUTS
    !    o ff -- wave function on logsbt grid
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i, i_done
    real*8 :: lng0, dlng, dlnr
    real*8 :: A, B, err, ln, i_r, wave_g0, dwave_g0
    character*150 :: info_str
    character(*), parameter :: func = 'sbt_import_spline'

    lng0 = log(r_grid_min)
    dlng = log(r_grid_inc)
    dlnr = lnrange / N

    i_done = 0
    if (lnr0 < lng0) then
       ! Need extrapolation to the left (try r^{L_r}):

       wave_g0 = wave_spl(1, 1)                        ! value at r_grid_min
       dwave_g0 = wave_spl(2, 1) / r_grid_min / dlng   ! deriv. at r_grid_min
!!$       dwave_g0 = (wave_spl(1, 2) - wave_spl(1, 1)) / &
!!$       &          (r_grid_min*r_grid_inc - r_grid_min)

       ! Fit   A * r^(l+1) + B * r^(l+2)   to value and derivative
       A = (  (L+2) * wave_g0 - r_grid_min*dwave_g0) / r_grid_min**(L+1)
       B = (- (L+1) * wave_g0 + r_grid_min*dwave_g0) / r_grid_min**(L+2)

       ! Check fit at second grid point
       err =  A * (r_grid_min*r_grid_inc)**(L+1) &
       &   +  B * (r_grid_min*r_grid_inc)**(L+2) &
       &   - wave_spl(1, 2)
       if (err > 1d-4) then
          write(info_str,"(2X,'** ',A,'; L=',I5,'; err=',ES10.2)") &
          &  'Extrapolation trouble', L, err
          call localorb_info(info_str)
          ! myid==0 should be OK; will probably be done on each node anyway.
       end if
       do i = 1, N
          ln = lnr0 + (i-1) * dlnr
          if (ln >= lng0) exit    ! left border is done
          ff(i) = A * exp((L+1) * ln) + B * exp((L+2) * ln)
       end do
       i_done = i - 1
    end if

    do i = i_done+1, N
       ln = lnr0 + (i-1) * dlnr
       if (ln > lng0 + (n_grid-1)*dlng) exit  ! hit right border
       ! i(r) = 1 + ln(r/r_grid_min) / ln(r_grid_inc) ! from invert_log_grid()
       i_r = 1.d0 + (ln - lng0) / dlng
       ff(i) = val_spline(i_r, wave_spl, n_grid)
    end do
    i_done = i - 1

    do i = i_done+1, N
       ff(i) = 0.d0    ! exponentially localized on log grid -> e^(-e^rho)
    end do

    ! u(r) -> f(r) = u(r)/r
    call logsbt_scale(N, 1, ff, lnr0, lnrange, -1.d0)

  end subroutine sbt_import_spline
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_aims/get_fnL_to_rowcol
  !  NAME
  !    get_fnL_to_rowcol
  !  SYNOPSIS

  subroutine get_fnL_to_rowcol(i_atom, &
  &                            n_bas, bas_atom, bas_fn, bas_l, bas_m, &
  &                            max_L, max_n_fnL, n_fnL, fnL_to_basfn, &
  &                            bas2rowcol, fnL_to_rowcol)

    !  PURPOSE
    !
    !    For a given atom, get indexing array fnL_to_rowcol.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom
    integer, intent(IN) :: n_bas
    integer, intent(IN) :: bas_atom(n_bas), bas_fn(n_bas)
    integer, intent(IN) :: bas_l(n_bas), bas_m(n_bas)
    integer, intent(IN) :: max_L, max_n_fnL
    integer, intent(IN) :: n_fnL(0:max_L)
    integer, intent(IN) :: fnL_to_basfn(max_n_fnL, 0:max_L)
    integer, intent(IN) :: bas2rowcol(n_bas)
    integer, intent(OUT) :: fnL_to_rowcol(-max_L:max_L, max_n_fnL, 0:max_L)

    !  INPUTS
    !    o i_atom -- current atom
    !    o n_bas, bas_atom, bas_fn, bas_l, bas_m -- Basis function definition
    !    o max_L, max_n_fnL -- Array dimensions
    !    o n_fnL -- Number of radial parts in angular momentum channel
    !    o fnL_to_basfn -- Indexing: i_fnL, L -> i_basfn
    !    o bas2rowcol -- Indexing: i_bas -> i_rowcol
    !  OUTPUTS
    !    o fnL_to_rowcol -- Indexing: M, i_fnL, L -> i_rowcol
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: i_bas, i_bas_fn, M, L, i_fnL
    character(*), parameter :: func = 'get_fnL_to_rowcol'


    fnL_to_rowcol = 0
    do i_bas = 1, n_bas
       i_bas_fn = bas_fn(i_bas)
       if (i_atom /= bas_atom(i_bas)) cycle
       M = bas_m(i_bas)
       L = bas_l(i_bas)

       do i_fnL = 1, n_fnL(L)
          if (fnL_to_basfn(i_fnL, L) == i_bas_fn) then
             fnL_to_rowcol(M, i_fnL, L) = bas2rowcol(i_bas)
          end if
       end do
    end do

  end subroutine get_fnL_to_rowcol
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_twocenter_triples
  !  NAME
  !    sbt_twocenter_triples
  !  SYNOPSIS

  subroutine sbt_twocenter_triples(N, lnr0, lnk0, lnrange, n_Rvec, Rvec, &
  & n_basfn1,       radius_1, max_L1, max_n_fnL1, n_fnL1, fnL1_to_basfn1, &
  & n_basfn2, ffk2, radius_2, max_L2, max_n_fnL2, n_fnL2, fnL2_to_basfn2, &
  & max_L1a, max_n_fnL1a, n_fnL1a, max_L1b, max_n_fnL1b, n_fnL1b, &
  & fnL1_to_L1a, fnL1_to_fnL1a, fnL1_to_L1b, fnL1_to_fnL1b, &
  & n_basfn1a, ffr1a, fnL1a_to_basfn1a, &
  & n_basfn1b, ffr1b, fnL1b_to_basfn1b, &
  & fnL1a_to_i1, fnL1b_to_i2, fnL2_to_i3, ovlp, d_ovlp, calc_deriv, times)

    !  PURPOSE
    !
    !    Calculate triple-overlaps(!) where two centers (1a and 1b) coincide.
    !
    !    This procedure is to be used for the calculation of triple-Coulomb
    !    integrals where the auxiliary center coincides with one of the two
    !    other centers.  This functions takes the resulting product as fn1a
    !    and fn1b and calculates the overlap with fn2 at Rvec.
    !
    !  USES

    use mpi_tasks, only: aims_stop, check_allocation
    use timing_core, only: start_timer, stop_timer
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: n_Rvec
    real*8, intent(IN) :: Rvec(3, n_Rvec)

    integer, intent(IN) :: max_L1, max_n_fnL1, n_basfn1
    real*8, intent(IN) :: radius_1(n_basfn1)
    integer, intent(IN) :: n_fnL1(0:max_L1)
    integer, intent(IN) :: fnL1_to_basfn1(max_n_fnL1, 0:max_L1)

    integer, intent(IN) :: max_L2, max_n_fnL2, n_basfn2
    real*8, intent(IN) :: ffk2(N, n_basfn2), radius_2(n_basfn2)
    integer, intent(IN) :: n_fnL2(0:max_L2)
    integer, intent(IN) :: fnL2_to_basfn2(max_n_fnL2, 0:max_L2)

    integer, intent(IN) :: max_L1a, max_n_fnL1a, n_fnL1a(0:max_L1a)
    integer, intent(IN) :: max_L1b, max_n_fnL1b, n_fnL1b(0:max_L1b)

    integer, intent(IN) :: fnL1_to_L1a(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_L1b(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_fnL1a(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_fnL1b(max_n_fnL1, 0:max_L1)

    integer, intent(IN) :: n_basfn1a, fnL1a_to_basfn1a(max_n_fnL1a, 0:max_L1a)
    integer, intent(IN) :: n_basfn1b, fnL1b_to_basfn1b(max_n_fnL1b, 0:max_L1b)
    real*8, intent(IN) :: ffr1a(N, n_basfn1a), ffr1b(N, n_basfn1b)

    integer, intent(IN) :: fnL1a_to_i1(max_n_fnL1a, 0:max_L1a)
    integer, intent(IN) :: fnL1b_to_i2(max_n_fnL1b, 0:max_L1b)
    integer, intent(IN) :: fnL2_to_i3(max_n_fnL2, 0:max_L2)

    real*8, intent(INOUT) :: ovlp(:,:,:,:)
    real*8, intent(INOUT) :: d_ovlp(:,:,:,:,:)
    logical, intent(IN) :: calc_deriv
    real*8, intent(INOUT) :: times(4, SBT_N_TIMES)

    !  INPUTS
    !    Functions:
    !    ...1: Product of basis with basbas at center_1
    !    ...1a: basis at center_1
    !    ...1b: basbas at center_1
    !    ...2: basis at center_2
    !    General:
    !    o N, lnr0, lnk0, lnrange -- Grid in r- and k-space
    !    o n_Rvec -- number of different Rvecs
    !    o Rvec -- distance vector between centers for fn1 and fn2
    !    Array sizes:
    !    o max_L{1,1a,1b,2} -- Angular momenta
    !    o n_basfn{1,1a,1b,2} -- Radial functions (ffki, radius_i)
    !    o max_n_fnL{1,1a,1b,2} -- max number of radial parts in one L-channel
    !    Radial function information
    !    o ffk2 -- Scaled (by k^1.5) SBTs of radial parts
    !    o radius_{1,2} -- Extent of radial parts
    !    o ffr1a, ffr1b -- Radial parts to be multiplied to each other
    !    Information on functions
    !    o n_fnL{1,1a,1b,2} -- Number of radial parts (fns) in given L-channel
    !    o fnL{1,1a,1b,2}_to_basfn{...} -- Corresponding i_basfn [ffk, radius]
    !    o fnL1_to_L1{a,b} -- L1{a,b} of generating radial parts
    !    o fnL1_to_fnL1{a,b} -- Index of generating radial parts within L-chan.
    !    o fnL*_to_i? -- indexing of ovlp(i1, i2, i3, n_Rvec)
    !    o times -- Timing information
    !    o calc_deriv - flag if d_ovlp should be calculated
    !  OUTPUTS
    !    o ovlp(i1, i2, i3, n_Rvec) -- triple overlap integrals
    !    o d_ovlp(i1, i2, i3, 3, n_Rvec) -- derivative of ovlp with respect to Rvec
    !    o times -- Timing information
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8, allocatable :: ffk1(:,:)
    integer, allocatable :: basfn1_to_L1(:)
    real*8, allocatable :: dummy_1(:), dummy_2(:)
    integer, allocatable :: fnL1_to_row(:,:,:), fnL2_to_col(:,:,:)
    integer :: n_row, i_row, n_col, i_col
    real*8, allocatable :: ovlp2(:,:), d_ovlp2(:,:,:)

    integer :: i_Rvec
    integer :: i_basfn1, i_basfn1a, i_basfn1b
    integer :: L1a, L1b, M1a, M1b, i_fnL1a, i_fnL1b
    integer :: i1, i2, i3
    integer :: L1, i_M1, M1, i_fnL1
    integer :: L2, i_fnL2, M2
    real*8, allocatable :: triple_Yr(:,:,:)
    integer, allocatable :: M1s(:,:,:)
    real*8 :: Cfn1ab_in_fn1
    integer :: info
    character(*), parameter :: func = 'sbt_twocenter_triples'

    call start_timer(times(:, SBT_TIME_WHOLE_TRIPLES))

    allocate(dummy_1(n_basfn1), dummy_2(n_basfn2), stat=info)
    call check_allocation(info, 'dummy_i', func)
    dummy_1 = 0.d0
    dummy_2 = 0.d0

    ! --- Prepare ffk1

    allocate(ffk1(N, n_basfn1), stat=info)  ! This can get quite expensive...
    call check_allocation(info, 'ffk1', func,N,n_basfn1)
    allocate(basfn1_to_L1(n_basfn1), stat=info)
    call check_allocation(info, 'basfn1_to_L1', func)

    call start_timer(times(:, SBT_TIME_TRIPLE_MULT))
    do L1 = 0, max_L1
       do i_fnL1 = 1, n_fnL1(L1)
          L1a = fnL1_to_L1a(i_fnL1, L1)
          L1b = fnL1_to_L1b(i_fnL1, L1)
          i_fnL1a = fnL1_to_fnL1a(i_fnL1, L1)
          i_fnL1b = fnL1_to_fnL1b(i_fnL1, L1)
          i_basfn1 = fnL1_to_basfn1(i_fnL1, L1)
          i_basfn1a = fnL1a_to_basfn1a(i_fnL1a, L1a)
          i_basfn1b = fnL1b_to_basfn1b(i_fnL1b, L1b)

          ffk1(:, i_basfn1) = ffr1a(:, i_basfn1a) * ffr1b(:, i_basfn1b)
          basfn1_to_L1(i_basfn1) = L1
       end do
    end do
    call stop_timer(times(:, SBT_TIME_TRIPLE_MULT), .true.)
    call start_timer(times(:, SBT_TIME_TRIPLE_SBT))

    ! f(r)             ! [assume full relative accuracy]
    call logsbt_scale(N, n_basfn1, ffk1, lnr0, lnrange, 1.5d0)
    ! f(r)*r**1.5      ! [assume full relative accuracy]
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, 1.5d0, &
    &                        n_basfn1, basfn1_to_L1, ffk1)
    ! f(k)*k**1.5      ! [error: epsilon]

    call stop_timer(times(:, SBT_TIME_TRIPLE_SBT), .true.)

    ! --- Prepare ovlp2(:,:) and its indexing

    allocate(fnL1_to_row(-max_L1:max_L1, max_n_fnL1, 0:max_L1), stat=info)
    call check_allocation(info, 'fnL1_to_row', func)
    allocate(fnL2_to_col(-max_L2:max_L2, max_n_fnL2, 0:max_L2), stat=info)
    call check_allocation(info, 'fnL2_to_col', func)

    call get_default_fnL_to_rowcol(max_L1, max_n_fnL1, n_fnL1, &
    &                              fnL1_to_row, n_row)
    call get_default_fnL_to_rowcol(max_L2, max_n_fnL2, n_fnL2, &
    &                              fnL2_to_col, n_col)

    allocate(ovlp2(n_row, n_col), stat=info)
    call check_allocation(info, 'ovlp2', func)

    if(calc_deriv) then
       allocate(d_ovlp2(3, n_row, n_col), stat=info)
       call check_allocation(info, 'd_ovlp2', func)
    endif

    ! --- ovlp(...) = 0.d0

    call start_timer(times(:, SBT_TIME_TRIPLE_COMB))
    ZERO_1A: do L1a = 0, max_L1a
       do i_fnL1a = 1, n_fnL1a(L1a)
          do M1a = -L1a, L1a
             i1 = fnL1a_to_i1(i_fnL1a, L1a) + M1A
             if (i1 <= 0) cycle
             ZERO_1B: do L1b = 0, max_L1b
                do i_fnL1b = 1, n_fnL1b(L1b)
                   do M1b = -L1b, L1b
                      i2 = fnL1b_to_i2(i_fnL1b, L1b) + M1B
                      if (i2 <= 0) cycle
                      ZERO_2: do L2 = 0, max_L2
                         do i_fnL2 = 1, n_fnL2(L2)
                            do M2 = -L2, L2
                               i3 = fnL2_to_i3(i_fnL2, L2) + m2
                               if (i3 <= 0) cycle
                               ovlp(i1, i2, i3, :) = 0.d0
                               if(calc_deriv) d_ovlp(i1, i2, i3, :, :) = 0.d0
                            end do
                         end do
                      end do ZERO_2
                   end do
                end do
             end do ZERO_1B
          end do
       end do
    end do ZERO_1A

    allocate(triple_Yr(-max_L1a:max_L1a, -max_L1b:max_L1b, 2), stat=info)
    call check_allocation(info, 'triple_Yr', func)
    allocate(M1s(-max_L1a:max_L1a, -max_L1b:max_L1b, 2), stat=info)
    call check_allocation(info, 'M1s', func)

    RVEC_L: do i_Rvec = 1, n_Rvec

       ! --- Get overlap between fn1 and fn2
       
     if(calc_deriv) then
       call sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, &
       &                    Rvec(:, i_Rvec), 1.d0, .false., &
       & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
       & n_basfn1, ffk1, radius_1, dummy_1, n_fnL1, fnL1_to_basfn1,fnL1_to_row,&
       & n_basfn2, ffk2, radius_2, dummy_2, n_fnL2, fnL2_to_basfn2,fnL2_to_col,&
       & ovlp2, times, d_ovlp2)
     else
       call sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, &
       &                    Rvec(:, i_Rvec), 1.d0, .false., &
       & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
       & n_basfn1, ffk1, radius_1, dummy_1, n_fnL1, fnL1_to_basfn1,fnL1_to_row,&
       & n_basfn2, ffk2, radius_2, dummy_2, n_fnL2, fnL2_to_basfn2,fnL2_to_col,&
       & ovlp2, times)
     endif

       ! --- Transform to triples (fn1a * fn1b * fn2)

       L1_LOOP: do L1 = 0, max_L1
          do i_fnL1 = 1, n_fnL1(L1)
             L1a = fnL1_to_L1a(i_fnL1, L1)
             L1b = fnL1_to_L1b(i_fnL1, L1)
             i_fnL1a = fnL1_to_fnL1a(i_fnL1, L1)
             i_fnL1b = fnL1_to_fnL1b(i_fnL1, L1)

             ! Get expansion coefficient to fn1 of the product of fn1a and fn1b.
             call triple_Y_real(L1a, L1b, L1, max_L1a, max_L1b, triple_Yr, M1s)

             do M1a = -L1a, L1a
                i1 = fnL1a_to_i1(i_fnL1a, L1a) + M1a
                if (i1 <= 0) cycle
                do M1b = -L1b, L1b
                   i2 = fnL1b_to_i2(i_fnL1b, L1b) + M1b
                   if (i2 <= 0) cycle

                   do i_M1 = 1, 2
                      M1 = M1s(M1a, M1b, i_M1)
                      if (abs(M1) > L1) cycle
                      Cfn1ab_in_fn1 = triple_Yr(M1a, M1b, i_M1)

                      i_row = fnL1_to_row(M1, i_fnL1, L1)
                      if (i_row <= 0) call aims_stop('fn1 not found', func)

                      L2_LOOP: do L2 = 0, max_L2
                         do i_fnL2 = 1, n_fnL2(L2)
                            do M2 = -L2, L2
                               i3 = fnL2_to_i3(i_fnL2, L2) + M2
                               if (i3 <= 0) cycle
                               i_col = fnL2_to_col(M2, i_fnL2, L2)
                               if (i_col <= 0) then
                                  call aims_stop('fn2 not found', func)
                               end if
                               ovlp(i1, i2, i3, i_Rvec) &
                               & = ovlp(i1, i2, i3, i_Rvec) &
                               &   + Cfn1ab_in_fn1 * ovlp2(i_row, i_col)
                               if(calc_deriv) then
                                  d_ovlp(i1, i2, i3, 1:3, i_Rvec) &
                                  & = d_ovlp(i1, i2, i3, 1:3, i_Rvec) &
                                  &   + Cfn1ab_in_fn1 * d_ovlp2(:, i_row, i_col)
                               endif
                            end do
                         end do
                      end do L2_LOOP
                   end do
                end do
             end do
          end do
       end do L1_LOOP

    end do RVEC_L

    call stop_timer(times(:, SBT_TIME_TRIPLE_COMB), .true.)

    deallocate(dummy_1, dummy_2)
    deallocate(triple_Yr, M1s)
    deallocate(ovlp2, fnL1_to_row, fnL2_to_col)
    deallocate(basfn1_to_L1, ffk1)
    if(allocated(d_ovlp2)) deallocate(d_ovlp2)

    call stop_timer(times(:, SBT_TIME_WHOLE_TRIPLES), .true.)

  end subroutine sbt_twocenter_triples

  subroutine my_sbt_twocenter_triples(N, lnr0, lnk0, lnrange, n_Rvec, Rvec, &
  & n_basfn1,       radius_1, max_L1, max_n_fnL1, n_fnL1, fnL1_to_basfn1, &
  & n_basfn2, ffk2, radius_2, max_L2, max_n_fnL2, n_fnL2, fnL2_to_basfn2, &
  & max_L1a, max_n_fnL1a, n_fnL1a, max_L1b, max_n_fnL1b, n_fnL1b, &
  & fnL1_to_L1a, fnL1_to_fnL1a, fnL1_to_L1b, fnL1_to_fnL1b, &
  & n_basfn1a, ffr1a, fnL1a_to_basfn1a, &
  & n_basfn1b, ffr1b, fnL1b_to_basfn1b, &
  & fnL1a_to_i1, fnL1b_to_i2, fnL2_to_i3, ovlp, d_ovlp, calc_deriv, times, index)

    !  PURPOSE
    !
    !    Calculate triple-overlaps(!) where two centers (1a and 1b) coincide.
    !
    !    This procedure is to be used for the calculation of triple-Coulomb
    !    integrals where the auxiliary center coincides with one of the two
    !    other centers.  This functions takes the resulting product as fn1a
    !    and fn1b and calculates the overlap with fn2 at Rvec.
    !
    !  USES
    use mpi_tasks, only: aims_stop, check_allocation
    use my_logsbt, only: my_logsbt_multi_driver
    use timing_core, only: start_timer, stop_timer
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: n_Rvec
    real*8, intent(IN) :: Rvec(3, n_Rvec)

    integer, intent(IN) :: max_L1, max_n_fnL1, n_basfn1
    real*8, intent(IN) :: radius_1(n_basfn1)
    integer, intent(IN) :: n_fnL1(0:max_L1)
    integer, intent(IN) :: fnL1_to_basfn1(max_n_fnL1, 0:max_L1)

    integer, intent(IN) :: max_L2, max_n_fnL2, n_basfn2
    real*8, intent(IN) :: ffk2(N, n_basfn2), radius_2(n_basfn2)
    integer, intent(IN) :: n_fnL2(0:max_L2)
    integer, intent(IN) :: fnL2_to_basfn2(max_n_fnL2, 0:max_L2)

    integer, intent(IN) :: max_L1a, max_n_fnL1a, n_fnL1a(0:max_L1a)
    integer, intent(IN) :: max_L1b, max_n_fnL1b, n_fnL1b(0:max_L1b)

    integer, intent(IN) :: fnL1_to_L1a(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_L1b(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_fnL1a(max_n_fnL1, 0:max_L1)
    integer, intent(IN) :: fnL1_to_fnL1b(max_n_fnL1, 0:max_L1)

    integer, intent(IN) :: n_basfn1a, fnL1a_to_basfn1a(max_n_fnL1a, 0:max_L1a)
    integer, intent(IN) :: n_basfn1b, fnL1b_to_basfn1b(max_n_fnL1b, 0:max_L1b)
    real*8, intent(IN) :: ffr1a(N, n_basfn1a), ffr1b(N, n_basfn1b)

    integer, intent(IN) :: fnL1a_to_i1(max_n_fnL1a, 0:max_L1a)
    integer, intent(IN) :: fnL1b_to_i2(max_n_fnL1b, 0:max_L1b)
    integer, intent(IN) :: fnL2_to_i3(max_n_fnL2, 0:max_L2)

    real*8, intent(INOUT) :: ovlp(:,:,:,:)
    real*8, intent(INOUT) :: d_ovlp(:,:,:,:,:)
    logical, intent(IN) :: calc_deriv
    real*8, intent(INOUT) :: times(4, SBT_N_TIMES)
    integer, intent(IN) :: index

    !  INPUTS
    !    Functions:
    !    ...1: Product of basis with basbas at center_1
    !    ...1a: basis at center_1
    !    ...1b: basbas at center_1
    !    ...2: basis at center_2
    !    General:
    !    o N, lnr0, lnk0, lnrange -- Grid in r- and k-space
    !    o n_Rvec -- number of different Rvecs
    !    o Rvec -- distance vector between centers for fn1 and fn2
    !    Array sizes:
    !    o max_L{1,1a,1b,2} -- Angular momenta
    !    o n_basfn{1,1a,1b,2} -- Radial functions (ffki, radius_i)
    !    o max_n_fnL{1,1a,1b,2} -- max number of radial parts in one L-channel
    !    Radial function information
    !    o ffk2 -- Scaled (by k^1.5) SBTs of radial parts
    !    o radius_{1,2} -- Extent of radial parts
    !    o ffr1a, ffr1b -- Radial parts to be multiplied to each other
    !    Information on functions
    !    o n_fnL{1,1a,1b,2} -- Number of radial parts (fns) in given L-channel
    !    o fnL{1,1a,1b,2}_to_basfn{...} -- Corresponding i_basfn [ffk, radius]
    !    o fnL1_to_L1{a,b} -- L1{a,b} of generating radial parts
    !    o fnL1_to_fnL1{a,b} -- Index of generating radial parts within L-chan.
    !    o fnL*_to_i? -- indexing of ovlp(i1, i2, i3, n_Rvec)
    !    o times -- Timing information
    !    o calc_deriv - flag if d_ovlp should be calculated
    !  OUTPUTS
    !    o ovlp(i1, i2, i3, n_Rvec) -- triple overlap integrals
    !    o d_ovlp(i1, i2, i3, 3, n_Rvec) -- derivative of ovlp with respect to Rvec
    !    o times -- Timing information
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE


    integer, allocatable :: basfn1_to_L1(:), basfn1_to_L1_Copy(:)
    real*8, allocatable :: dummy_1(:), dummy_2(:)
    integer, allocatable :: fnL1_to_row(:,:,:), fnL2_to_col(:,:,:)
    integer :: n_row, i_row, n_col, i_col
    real*8, allocatable :: ovlp2(:,:), d_ovlp2(:,:,:)

    integer :: i_Rvec
    integer :: i_basfn1, i_basfn1a, i_basfn1b
    integer :: L1a, L1b, M1a, M1b, i_fnL1a, i_fnL1b
    integer :: i1, i2, i3
    integer :: L1, i_M1, M1, i_fnL1
    integer :: L2, i_fnL2, M2
    real*8, allocatable :: triple_Yr(:,:,:)
    integer, allocatable :: M1s(:,:,:)
    real*8 :: Cfn1ab_in_fn1
    integer :: info
    character(*), parameter :: func = 'sbt_twocenter_triples'

    integer i,j

    call start_timer(times(:, SBT_TIME_WHOLE_TRIPLES))

    allocate(dummy_1(n_basfn1), dummy_2(n_basfn2), stat=info)
    call check_allocation(info, 'dummy_i', func)
    dummy_1 = 0.d0
    dummy_2 = 0.d0

    ! --- Prepare ffk1

    if (index.ne.index_old) then
       call perfon('sbt_ffk1')
       if (allocated(ffk1)) deallocate(ffk1)
       allocate(ffk1(N, n_basfn1), stat=info)  ! This can get quite expensive...
       call check_allocation(info, 'ffk1', func,N,n_basfn1)
       allocate(basfn1_to_L1(n_basfn1), stat=info)
       call check_allocation(info, 'basfn1_to_L1', func)
       
!       allocate(ffk1_copy(N, n_basfn1), stat=info)  ! This can get quite expensive...
!       call check_allocation(info, 'ffk1_copy', func,N,n_basfn1)
!       allocate(basfn1_to_L1_copy(n_basfn1), stat=info)
!       call check_allocation(info, 'basfn1_to_L1_copy', func)
       
       call start_timer(times(:, SBT_TIME_TRIPLE_MULT))
       do L1 = 0, max_L1
          do i_fnL1 = 1, n_fnL1(L1)
             L1a = fnL1_to_L1a(i_fnL1, L1)
             L1b = fnL1_to_L1b(i_fnL1, L1)
             i_fnL1a = fnL1_to_fnL1a(i_fnL1, L1)
             i_fnL1b = fnL1_to_fnL1b(i_fnL1, L1)
             i_basfn1 = fnL1_to_basfn1(i_fnL1, L1)
             i_basfn1a = fnL1a_to_basfn1a(i_fnL1a, L1a)
             i_basfn1b = fnL1b_to_basfn1b(i_fnL1b, L1b)
             
             ffk1(:, i_basfn1) = ffr1a(:, i_basfn1a) * ffr1b(:, i_basfn1b)
             basfn1_to_L1(i_basfn1) = L1
          end do
       end do
       call stop_timer(times(:, SBT_TIME_TRIPLE_MULT), .true.)
       call start_timer(times(:, SBT_TIME_TRIPLE_SBT))
       
       ! f(r)             ! [assume full relative accuracy]
       call logsbt_scale(N, n_basfn1, ffk1, lnr0, lnrange, 1.5d0)
       ! f(r)*r**1.5      ! [assume full relative accuracy]
       
       
       call logsbt_multi_driver(N, lnr0, lnk0, lnrange, 1.5d0, &
            &                        n_basfn1, basfn1_to_L1, ffk1)

       deallocate(basfn1_to_L1)
       index_old = index
       call perfoff
    end if
!
!    call my_logsbt_multi_driver(N, lnr0, lnk0, lnrange, 1.5d0, &
!    &                        n_basfn1, basfn1_to_L1_copy, ffk1_copy)
!
!
!! check
!   do i = 1, N
!     do j=1,n_basfn1
!       if (ffk1(i,j) .ne. ffk1_copy(i,j)) then
!         print *,"ffk1 : ",i,j,ffk1(i,j), ffk1_copy(i,j)
!         stop
!       endif
!    enddo
!  enddo
!
!  do j=1,n_basfn1
!    if (basfn1_to_L1_copy(j) .ne. basfn1_to_L1(j)) then
!      print *,"basfn1_to_L1: ",j,basfn1_to_L1(j),basfn1_to_L1_copy(j)
!    endif
!  enddo

    ! f(k)*k**1.5      ! [error: epsilon]

    call stop_timer(times(:, SBT_TIME_TRIPLE_SBT), .true.)

    ! --- Prepare ovlp2(:,:) and its indexing

    allocate(fnL1_to_row(-max_L1:max_L1, max_n_fnL1, 0:max_L1), stat=info)
    call check_allocation(info, 'fnL1_to_row', func)
    allocate(fnL2_to_col(-max_L2:max_L2, max_n_fnL2, 0:max_L2), stat=info)
    call check_allocation(info, 'fnL2_to_col', func)
    call get_default_fnL_to_rowcol(max_L1, max_n_fnL1, n_fnL1, &
    &                              fnL1_to_row, n_row)
    call get_default_fnL_to_rowcol(max_L2, max_n_fnL2, n_fnL2, &
    &                              fnL2_to_col, n_col)

    allocate(ovlp2(n_row, n_col), stat=info)
    call check_allocation(info, 'ovlp2', func)

    if(calc_deriv) then
       allocate(d_ovlp2(3, n_row, n_col), stat=info)
       call check_allocation(info, 'd_ovlp2', func)
    endif


    ! --- ovlp(...) = 0.d0

    call start_timer(times(:, SBT_TIME_TRIPLE_COMB))
    ZERO_1A: do L1a = 0, max_L1a
       do i_fnL1a = 1, n_fnL1a(L1a)
          do M1a = -L1a, L1a
             i1 = fnL1a_to_i1(i_fnL1a, L1a) + M1A
             if (i1 <= 0) cycle
             ZERO_1B: do L1b = 0, max_L1b
                do i_fnL1b = 1, n_fnL1b(L1b)
                   do M1b = -L1b, L1b
                      i2 = fnL1b_to_i2(i_fnL1b, L1b) + M1B
                      if (i2 <= 0) cycle
                      ZERO_2: do L2 = 0, max_L2
                         do i_fnL2 = 1, n_fnL2(L2)
                            do M2 = -L2, L2
                               i3 = fnL2_to_i3(i_fnL2, L2) + m2
                               if (i3 <= 0) cycle
                               ovlp(i1, i2, i3, :) = 0.d0
                               if(calc_deriv) d_ovlp(i1, i2, i3, :, :) = 0.d0
                            end do
                         end do
                      end do ZERO_2
                   end do
                end do
             end do ZERO_1B
          end do
       end do
    end do ZERO_1A

    allocate(triple_Yr(-max_L1a:max_L1a, -max_L1b:max_L1b, 2), stat=info)
    call check_allocation(info, 'triple_Yr', func)
    allocate(M1s(-max_L1a:max_L1a, -max_L1b:max_L1b, 2), stat=info)
    call check_allocation(info, 'M1s', func)

    RVEC_L: do i_Rvec = 1, n_Rvec

       ! --- Get overlap between fn1 and fn2
     if(calc_deriv) then
       call my_sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, &
       &                    Rvec(:, i_Rvec), 1.d0, .false., &
       & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
       & n_basfn1, ffk1, radius_1, dummy_1, n_fnL1, fnL1_to_basfn1,fnL1_to_row,&
       & n_basfn2, ffk2, radius_2, dummy_2, n_fnL2, fnL2_to_basfn2,fnL2_to_col,&
       & ovlp2, times, d_ovlp2)
     else

       call my_sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, &
       &                    Rvec(:, i_Rvec), 1.d0, .false., &
       & max_L1, max_L2, max_n_fnL1, max_n_fnL2, &
       & n_basfn1, ffk1, radius_1, dummy_1, n_fnL1, fnL1_to_basfn1,fnL1_to_row,&
       & n_basfn2, ffk2, radius_2, dummy_2, n_fnL2, fnL2_to_basfn2,fnL2_to_col,&
       & ovlp2, times)
     endif
!!!!$omp end critical

       ! --- Transform to triples (fn1a * fn1b * fn2)

       L1_LOOP: do L1 = 0, max_L1
          do i_fnL1 = 1, n_fnL1(L1)
             L1a = fnL1_to_L1a(i_fnL1, L1)
             L1b = fnL1_to_L1b(i_fnL1, L1)
             i_fnL1a = fnL1_to_fnL1a(i_fnL1, L1)
             i_fnL1b = fnL1_to_fnL1b(i_fnL1, L1)

             ! Get expansion coefficient to fn1 of the product of fn1a and fn1b.
             call triple_Y_real(L1a, L1b, L1, max_L1a, max_L1b, triple_Yr, M1s)

             do M1a = -L1a, L1a
                i1 = fnL1a_to_i1(i_fnL1a, L1a) + M1a
                if (i1 <= 0) cycle
                do M1b = -L1b, L1b
                   i2 = fnL1b_to_i2(i_fnL1b, L1b) + M1b
                   if (i2 <= 0) cycle

                   do i_M1 = 1, 2
                      M1 = M1s(M1a, M1b, i_M1)
                      if (abs(M1) > L1) cycle
                      Cfn1ab_in_fn1 = triple_Yr(M1a, M1b, i_M1)

                      i_row = fnL1_to_row(M1, i_fnL1, L1)
                      if (i_row <= 0) call aims_stop('fn1 not found', func)

                      L2_LOOP: do L2 = 0, max_L2
                         do i_fnL2 = 1, n_fnL2(L2)
                            do M2 = -L2, L2
                               i3 = fnL2_to_i3(i_fnL2, L2) + M2
                               if (i3 <= 0) cycle
                               i_col = fnL2_to_col(M2, i_fnL2, L2)
                               if (i_col <= 0) then
                                  call aims_stop('fn2 not found', func)
                               end if
                               ovlp(i1, i2, i3, i_Rvec) &
                               & = ovlp(i1, i2, i3, i_Rvec) &
                               &   + Cfn1ab_in_fn1 * ovlp2(i_row, i_col)
                               if(calc_deriv) then
                                  d_ovlp(i1, i2, i3, 1:3, i_Rvec) &
                                  & = d_ovlp(i1, i2, i3, 1:3, i_Rvec) &
                                  &   + Cfn1ab_in_fn1 * d_ovlp2(:, i_row, i_col)
                               endif
                            end do
                         end do
                      end do L2_LOOP
                   end do
                end do
             end do
          end do
       end do L1_LOOP
    end do RVEC_L

    call stop_timer(times(:, SBT_TIME_TRIPLE_COMB), .true.)

    deallocate(dummy_1, dummy_2)
    deallocate(triple_Yr, M1s)
    deallocate(ovlp2, fnL1_to_row, fnL2_to_col)

    if(allocated(d_ovlp2)) deallocate(d_ovlp2)

    call stop_timer(times(:, SBT_TIME_WHOLE_TRIPLES), .true.)
  end subroutine my_sbt_twocenter_triples

  subroutine my_ffk1_cleanup
    if (allocated(ffk1)) deallocate(ffk1)
    index_old=-1
  end subroutine my_ffk1_cleanup

  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/get_default_fnL_to_rowcol
  !  NAME
  !    get_default_fnL_to_rowcol
  !  SYNOPSIS

  subroutine get_default_fnL_to_rowcol(max_L, max_n_fnL, n_fnL, &
  &                                    fnL_to_rowcol, n_rowcol)

    !  PURPOSE
    !
    !     Generate a default indexing for basis functions grouped
    !     by (i_fnL, L).
    !
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: max_L, max_n_fnL
    integer, intent(IN) :: n_fnL(0:max_L)
    integer, intent(OUT) :: fnL_to_rowcol(-max_L:max_L, max_n_fnL, 0:max_L)
    integer, intent(OUT) :: n_rowcol

    !  INPUTS
    !    o max_L -- Maximum angular momentum
    !    o max_n_fnL -- Maximum number of radial parts within L-channel
    !    o n_fnL -- Number of radial parts within given L-channel
    !  OUTPUTS
    !    o fnL_to_rowcol -- (M, i_fnL, L) -> i_rowcol (running index)
    !    o n_rowcol -- maxval(fnL_to_rowcol)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: n_rowcol_uptonow
    integer :: L, i_fnL, M
    character(*), parameter :: func = 'get_default_fnL_to_rowcol'

    n_rowcol_uptonow = 0
    fnL_to_rowcol = 0
    do L = 0, max_L
       do i_fnL = 1, n_fnL(L)
          do M = -L, L
             n_rowcol_uptonow = n_rowcol_uptonow + 1
             fnL_to_rowcol(M, i_fnL, L) = n_rowcol_uptonow
          end do
       end do
    end do
    n_rowcol = n_rowcol_uptonow

  end subroutine get_default_fnL_to_rowcol
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/check_need_n_radius
  !  NAME
  !    check_need_n_radius
  !  SYNOPSIS

  subroutine check_need_n_radius(max_L, max_n_fnL, n_fnL, fnL_to_basfn, &
  &                              n_basfn, radius, fnL_to_rowcol,  &
  &                              need, max_radius)

    !  PURPOSE
    !     Check for given radial parts which are needed and what their
    !     maximum radius is
    !  USES

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: max_L, max_n_fnL, n_basfn
    integer, intent(IN) :: n_fnL(0:max_L), fnL_to_basfn(max_n_fnL, 0:max_L)
    real*8, intent(IN) :: radius(n_basfn)
    integer, intent(IN) :: fnL_to_rowcol(-max_L:max_L, max_n_fnL, 0:max_L)
    logical, intent(OUT) :: need(max_n_fnL, 0:max_L)
    real*8, intent(OUT) :: max_radius

    !  INPUTS
    !    o max_L, max_n_fnL, n_basfn -- array dimensions
    !    o n_fnL -- number of radial parts in given L channel
    !    o fnL_to_basfn -- (i_fnL, L) -> i_basfn
    !    o radius -- i_basfn -> radius
    !    o fnL_to_rowcol -- (i_m, i_fnL, L) -> position in array
    !  OUTPUTS
    !    o need -- (i_fnL, L) is needed (at least one i_m with rowcol > 0)
    !    o max_radius -- maximum radius of needed radial parts
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: L, i_fnL, i_basfn
    character(*), parameter :: func = 'check_need_n_radius'

    need = .false.
    max_radius = 0.d0
    do L = 0, max_L
       do i_fnL = 1, n_fnL(L)
          i_basfn = fnL_to_basfn(i_fnL, L)
          if (any(fnL_to_rowcol(-L:L, i_fnL, L) > 0)) then
             need(i_fnL, L) = .true.
             max_radius = max(max_radius, radius(i_basfn))
          end if
       end do
    end do

  end subroutine check_need_n_radius
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/multipole_prefac
  !  NAME
  !    multipole_prefac
  !  SYNOPSIS

  real*8 function multipole_prefac(L1, L2)

    !  PURPOSE
    !
    !    Calculate the prefactor for multipole Coulomb interaction:
    !    
    !    radI = multipole_prefac(L1, L2) * moment_1 * moment_2 &
    !    &      / Rabs**(L1+L2+1).
    !
    !    V(m1, m2) = 4*pi * (2*pi)**1.5 * cyl(m1, m2) * radI
    !
    !  USES

    use constants, only: pi
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: L1, L2

    !  INPUTS
    !    o L1, L2 -- Angular momenta
    !  OUTPUTS
    !    multipole_prefac
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8, allocatable :: Besfac(:)
    integer :: L, info
    character(*), parameter :: func = 'multipole_prefac'

    allocate(Besfac(0:L1+L2), stat=info)
    call check_allocation(info, 'Besfac', func)

    Besfac(0) = sqrt(2.d0 / pi)
    do L = 1, L1+L2
       Besfac(L) = Besfac(L-1) / (2*L+1)    ! == 1 / (2*L+1)!!
    end do

    multipole_prefac = 4*pi * (2*L1+1)*(2*L2+1) * Besfac(L1) * Besfac(L2) &
    &                  / dble((2*L1+2*L2+1) * Besfac(L1+L2))

    deallocate(Besfac)

  end function multipole_prefac
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_multipole_clmb
  !  NAME
  !    sbt_multipole_clmb
  !  SYNOPSIS

  subroutine sbt_multipole_clmb(Rvec, max_L1, L1, max_L2, L2, moment_prod, V)

    !  PURPOSE
    !
    !    Calculate multipole far-field interaction matrix elements.
    !
    !  USES

    use constants, only: pi
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: Rvec(3)
    integer, intent(IN) :: max_L1, max_L2
    integer, intent(IN) :: L1, L2
    real*8, intent(IN) :: moment_prod
    real*8, intent(OUT) :: V(-max_L1:max_L1, -max_L2:max_L2)

    !  INPUTS
    !    o Rvec -- distance vector between the two centers
    !    o max_Li -- dimensionalities
    !    o Li -- angular momenta
    !    o moment_prod -- the product of the two multipole moments
    !                     (as given in multipole_basbas_fn(i_basbasfn))
    !  OUTPUTS
    !    o V -- Coulomb interaction matrix elements
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: LL
    real*8 :: Rabs, fac
    integer :: info
    character(*), parameter :: func = 'sbt_multipole_ovlp'

    LL = L1 + L2

    ! Angular integration
    call sum_triple_Y_YLM_real(L1, L2, LL, Rvec, max_L1, max_L2, V)

    ! Radial dependence
    Rabs = sqrt(dot_product(Rvec, Rvec))

    fac = 4*pi * (2*pi)**1.5d0 * (-1)**((L1-L2-LL)/2) &
    &     * moment_prod * multipole_prefac(L1, L2) / Rabs**(L1+L2+1)

    V = fac * V

  end subroutine sbt_multipole_clmb
  !******
end module sbt_overlap
!******
