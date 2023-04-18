!****h* FHI-aims/tight_binding_auxmat
!  NAME
!    tight_binding_auxmat
!  SYNOPSIS

module tight_binding_auxmat

  !  PURPOSE
  !
  !    This is an FHI-aims aware wrapper around sbt_overlap_tb.f90 with an
  !    initializer which prepares splines for all combinations of radial parts
  !    (basbasfns) and a workhorse which can calculate an atom-atom block of
  !    the Coulomb (or overlap, or HSE) matrix for a given distance vector.
  !
  !    Originally, I intended this module to be general and to depend on
  !    prodbas only upon initialization.  I myself have broken this design
  !    often enough, so far.  -- JW
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
  !    Release version, FHI-aims (2011).
  !  SOURCE

  character*20, private :: ovlp_type_name  ! For final timing

  ! Have (multipolar) far-field interaction
  logical :: have_mp_far, have_Ewald
  real*8, allocatable, private :: multipoles(:) ! n_fn
  ! Range of (nontrivial) interaction (ia).
  !   This might look overly large, but this
  !   1) simplifies the treatment of HSE.
  !   2) is general enough for monkey patching things for V_q.
  real*8, allocatable, private :: iarange(:,:)  ! n_fn, n_fn
  ! Same for whole species.
  real*8, allocatable :: iarange_species(:,:) ! n_species, n_species

  ! --- Splines and their parameters

  ! The radial integrals wil be stored as a logarithmic spline for ~0.5 AA <=
  !  |R| <= ~50 AA.  The following numbers define the actual spline properties:

  ! Be aware that, as this is a finite spline and not a periodic function, the
  ! right-hand-side /is/ one of the n_spl points.
  integer, private :: n_spl      ! Number of points
  real*8, private :: spl_ln1     ! Log of first point
  real*8, private :: spl_lnrange ! Range in log-space
  real*8, private :: spl_lnN     ! Log of last point
  real*8, private :: spl_dln
  integer :: i_minR, i_maxR      ! First and last point within sbtgrid
  integer :: use_drop_fac        ! Coarsing of splines

  integer, private :: tb_max_L   ! Maximum number of radial integrals
  !                              ! for a given pair.
  integer, private :: n_radI     ! Number of radial integrals
  integer, allocatable, private :: Lfnfn2radI(:,:,:) ! 0:tb_max_L, n_fn, n_fn

  real*8, allocatable, private :: radI_onsites(:,:,:,:) ! i_fnL,j_fnL,L,i_spec
  real*8, allocatable, private :: radI_splines(:,:)     ! 0:n_spl+1, n_radI


  ! Set by atomic_field_transform; possibly used by get_Gaussian_...().
  real*8 :: power_bias_int
  real*8 :: singularity_lifting_chi

  ! --- Buffers

  real*8, allocatable, private :: tmp_radI(:,:,:), tmp_fullI(:,:,:,:)

  real*8 :: Gaussian_gamma

  ! --- Multipole moments

  ! Each L-channel for each species should contain exactly one radial part
  ! with a non-trivial multipole moment.  For a given (L, i_species) point to:

  ! * The index within this channel:
  integer, allocatable, private :: Lsp2mp_bb_fnLsp(:,:)
  ! i_fnLsp = Lsp2mp_bb_fnLsp(L, i_species)

  ! * The global radial part number:
  integer, allocatable, private :: Lsp2mp_bb_fn(:,:)   ! (0:max_L, n_species)
  ! i_basbas_fn = Lsp2mp_bb_fn(L, i_species) == Lsp2basbas_fn(i_fnLsp, ...)

  ! * The per-atom index (i_basbas_sp) of the M=0 function:
  integer, allocatable, private :: Lsp2mp_bb_sp(:,:)   ! (0:max_L, n_species)
  ! i_basbas_sp = Lsp2mp_bb_sp(L, i_species)
  ! i_basbas = atom2basbas_off(i_atom) + i_basbas_sp + M

  ! * The actual multipole moment
  real*8, allocatable, private :: Lsp2mp_bb_moment(:,:)  ! (0:max_L,n_species)
  ! Lsp2mp_bb_moment(L, i_species) == multipole_basbas_fn(i_basbas_fn)

  ! --- Timers

  real*8, private :: time_init(4)
  real*8, private :: time_rest(4)


contains

  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/initialize_tb_auxmat
  !  NAME
  !    initialize_tb_auxmat
  !  SYNOPSIS

  subroutine initialize_tb_auxmat(drop_fac, ovlp_type)

    !  PURPOSE
    !
    !    Initialize the arrays needed for a tight-binding like auxiliary
    !    Coulomb (or overlap, or HSE) matrix construction.
    !
    !  USES

    use dimensions
    use runtime_choices
    use grids
    use prodbas
    use synchronize_mpi_basic
    use geometry
    use bravais
    use localorb_io
    use sbt_overlap
    use sbt_overlap_aims, only: sbt_atomic_field_transforms
    use timing_core, only: start_timer, stop_timer
    use bspline, only: get_subspline
    use mpi_tasks, only: aims_stop, check_allocation, myid, n_tasks
    use sbt_overlap_tb, only: get_radial_integral_onsite, &
        get_radial_integral_spline
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: drop_fac
    integer, intent(IN) :: ovlp_type

    !  INPUTS
    !    o drop_fac -- Factor by which to reduce the density of points
    !                Should probably be something in the order of 1 or 4.
    !    o ovlp_type -- Either OVLP_TYPE_OVERLAP, OVLP_TYPE_COULOMB,
    !                   OVLP_TYPE_HSE (use runtime_choices:hse_omega_hf),
    !                   or OVLP_TYPE_CUT (use ...:cutCb_*).
    !  OUTPUTS
    !    none [arrays of this module are initialized]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: lnk0, lnr0, lnt0, lnrange, lntrange, dln
    real*8 :: minR, maxR
    integer :: N, Nt, max_LL
    integer :: i
    real*8, allocatable :: wsave(:), KKbars(:,:)
    real*8, allocatable :: ffk(:,:)
    real*8, allocatable :: int_spl(:,:)
    real*8 :: drop_err, max_drop_err
    real*8 :: chrad_1, chrad_2, firad_1, firad_2, this_range
    logical :: do_bare_coulomb
    integer :: i_L, i_LL, L, LL, i_basbas_fn, i_basbas, i_spec_bb
    integer :: i_bbfn_1, i_bbfn_2, i_fnL, L1, L2, i_m, M, i_radI, i_fnL1, i_fnL2
    integer :: i_atom, i_species, i_species_1, i_species_2
    integer :: bb_off, n_spbb_uptonow, i_iarange, i_this
    real*8 :: val(1)
    integer :: n_radI_uptonow, n_jobs_uptonow, n_basbas_uptonow
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'initialize_tb_auxmat'

    ovlp_type_name =  OVLP_TYPE_NAMES(ovlp_type)
    have_mp_far = (ovlp_type == OVLP_TYPE_COULOMB)
    have_Ewald = .false.

    if(output_priority .le. OL_norm) then
      write(info_str, "(2X,A,A,A)") &
      & 'Preparing TB like auxiliary 2-center ', trim(ovlp_type_name), ' matrix'
      call localorb_info(info_str)
    endif

    time_init = 0.d0
    time_rest = 0.d0
    call start_timer(time_init)

    ! --- Prepare local convenience variables

    lnk0 = sbtgrid_lnk0
    lnr0 = sbtgrid_lnr0
    lnt0 = lnk0 + lnr0
    N = sbtgrid_N
    Nt = 2*N
    lnrange = sbtgrid_lnrange
    lntrange = 2*lnrange
    dln = lnrange/N
    tb_max_L = maxval(basbasfn_l)
    max_LL = 2*tb_max_L

    ! --- Prepare interaction ranges

    if (have_mp_far) then
       allocate(multipoles(n_basbas_fns), stat=info)
       call check_allocation(info, 'multipoles', func,n_basbas_fns)
       multipoles = multipole_basbas_fn
    end if
    allocate(iarange(n_basbas_fns, n_basbas_fns), stat=info)
    call check_allocation(info, 'iarange', func,n_basbas_fns,n_basbas_fns)
    do i_bbfn_1 = 1, n_basbas_fns
       chrad_1 = charge_radius_basbas_fn(i_bbfn_1)
       firad_1 = field_radius_basbas_fn(i_bbfn_1)
       i_species_1 = basbasfn_species(i_bbfn_1)
       do i_bbfn_2 = 1, n_basbas_fns
          chrad_2 = charge_radius_basbas_fn(i_bbfn_2)
          firad_2 = field_radius_basbas_fn(i_bbfn_2)
          i_species_2 = basbasfn_species(i_bbfn_2)

          if (ovlp_type == OVLP_TYPE_HSE) then
             this_range = min(chrad_1+firad_2, firad_1+chrad_2)
          else if (ovlp_type == OVLP_TYPE_LR) then
             this_range = min(chrad_1+firad_2, firad_1+chrad_2)
          else if (ovlp_type == OVLP_TYPE_CUT) then
             this_range = chrad_1 + chrad_2 + cutCb_rcut * cutCb_width**6
          else
             this_range = chrad_1 + chrad_2 
!             this_range = chrad_1 + chrad_2 + cutCb_rcut * cutCb_width**6
          end if

          iarange(i_bbfn_1, i_bbfn_2) = this_range
       end do
    end do

    ! --- Get "sparse log-grid" parameters.

    ! Get minimum & maximum atom-atom separation for which to prepare splines.
    call min_atomic_dist(lattice_vector, coords, minR)
    ! test write(use_unit,*) "minR: ", minR
    maxR = maxval(iarange)

    if(minR .gt. maxR) minR = 0.9*maxR

    i_minR = floor(logsbt_r2ilog(minR, lnr0, dln)) - 1
    i_maxR = ceiling(logsbt_r2ilog(maxR, lnr0, dln)) + 1
    ! Make sure that the number of steps is dividable by drop_fac:
    i_maxR = i_maxR + modulo(i_maxR - i_minR, drop_fac)
    use_drop_fac = drop_fac
    if (i_minR < 1 .or. i_maxR > N) then
       ! This should not happen with default settings because it means that
       ! atoms are either much closer to each other than 1e-4 or that the
       ! charge radii + cutting parameter is much larger than 1e4.
       call aims_stop('Insufficient sbtgrid range', func)
    end if


    n_spl = (i_maxR - i_minR) / drop_fac + 1
    spl_ln1 = log(logsbt_ilog2r(dble(i_minR), lnr0, dln))
    spl_lnN = log(logsbt_ilog2r(dble(i_maxR), lnr0, dln))
    spl_dln = (spl_lnN - spl_ln1) / (n_spl - 1)
    spl_lnrange = n_spl * spl_dln   ! Different from spl_lnN

    ! --- Prepare logSBT biases

    allocate(ffk(N, n_basbas_fns), stat=info)
    call check_allocation(info, 'ffk', func,N,n_basbas_fns)

    call sbt_atomic_field_transforms(N, lnr0, lnk0, lnrange, ovlp_type, &
    &          n_basbas_fns, basbasfn_l, basbasfn_species, basbas_wave_spl, &
    &                                ffk, power_bias_int)

!    call sbt_atomic_field_transforms(N, lnr0, lnk0, lnrange, OVLP_TYPE_CUT, &
!   &          n_basbas_fns, basbasfn_l, basbasfn_species, basbas_wave_spl, &
!   &                                ffk, power_bias_int)
    ! --- Get SBT kernels

    allocate(KKbars(Nt, 0:max_LL), wsave(2*Nt+15), stat=info)
    call check_allocation(info, 'KKbars, wsave', func,Nt,max_LL+1,2*Nt+15)
    call dffti(Nt, wsave)
    do LL = 0, max_LL
       call logsbt_kernel(Nt, KKbars(:, LL), lnt0, lntrange, LL, &
       &                  power_bias_int, wsave)
    end do


    ! --- Prepare radI index arrays

    allocate(Lfnfn2radI(0:tb_max_L, n_basbas_fns, n_basbas_fns), stat=info)
    call check_allocation(info, 'Lfnfn2radI', func,tb_max_L+1,n_basbas_fns,n_basbas_fns)
    Lfnfn2radI = 0
    n_radI_uptonow = 0
    do i_bbfn_1 = 1, n_basbas_fns
       L1 = basbasfn_l(i_bbfn_1)
       do i_bbfn_2 = 1, i_bbfn_1
          L2 = basbasfn_l(i_bbfn_2)
          do i_L = 0, min(L1, L2)
             n_radI_uptonow = n_radI_uptonow + 1
             Lfnfn2radI(i_L, i_bbfn_1, i_bbfn_2) = n_radI_uptonow
             Lfnfn2radI(i_L, i_bbfn_2, i_bbfn_1) = n_radI_uptonow
          end do
       end do
    end do
    n_radI = n_radI_uptonow

    ! --- Prepare onsites

    if (power_bias_int /= 0.d0) call aims_stop('Invalid power_bias_int', func)
    allocate(radI_onsites(max_n_basbas_fnLsp, max_n_basbas_fnLsp, &
    &                     0:tb_max_L, n_species), stat=info)
    call check_allocation(info, 'radI_onsites', func,max_n_basbas_fnLsp,max_n_basbas_fnLsp, &
                          tb_max_L+1,n_species)
    radI_onsites = 0.d0
    do i_species = 1, n_species
       do L = 0, tb_max_L
          do i_fnL1 = 1, Lsp2n_basbas_fnLsp(L, i_species)
             i_bbfn_1 = Lsp2basbas_fn(i_fnL1, L, i_species)
             do i_fnL2 = 1, i_fnL1
                i_bbfn_2 = Lsp2basbas_fn(i_fnL2, L, i_species)
                call get_radial_integral_onsite(N, lnr0, lnk0, lnrange, &
                &                               1, L, ffk(:, i_bbfn_1), &
                &                               1, L, ffk(:, i_bbfn_2), val)
                radI_onsites(i_fnL1, i_fnL2, L, i_species) = val(1)
                radI_onsites(i_fnL2, i_fnL1, L, i_species) = val(1)
             end do
          end do
       end do
    end do

    ! --- Prepare splines

    if (power_bias_int /= 0.d0) call aims_stop('Invalid power_bias_int', func)
    if(output_priority .le. OL_norm) then
       write(info_str, "(2X,'| ',A,I6,A,I6,A,I10,A)") &
       & "Need to store", n_radI, " splines of", n_spl, " points taking ", &
       & ( ((n_spl+2)*n_radI) / 2**10 + 1)*8 , " KiB."
       call localorb_info(info_str)
    endif
    allocate(radI_splines(0:n_spl+1, n_radI), stat=info)
    call check_allocation(info, 'radI_splines', func,n_spl+2,n_radI)
    radI_splines = 0.d0

    allocate(int_spl(0:N+1, 0:tb_max_L), stat=info)
    call check_allocation(info, 'int_spl', func,N+2,tb_max_L+1)
    drop_err = 0.d0

    n_jobs_uptonow = 0
    do i_bbfn_1 = 1, n_basbas_fns
       L1 = basbasfn_l(i_bbfn_1)
       do i_bbfn_2 = 1, i_bbfn_1
          L2 = basbasfn_l(i_bbfn_2)
          n_jobs_uptonow = n_jobs_uptonow + 1
          if (myid == mod(n_jobs_uptonow, n_tasks)) then
             call get_radial_integral_spline(N, lnr0, lnk0, lnrange, &
             &                               max_LL, KKbars, wsave, &
             &                               1, L1, ffk(:, i_bbfn_1), &
             &                               1, L2, ffk(:, i_bbfn_2), &
             &                               tb_max_L, int_spl)
             do i_LL = 0, min(L1, L2)
                LL = abs(L1 - L2) + 2*i_LL
                i_radI = Lfnfn2radI(i_LL, i_bbfn_1, i_bbfn_2)
                call get_subspline(N, 1, i_minR, i_maxR, drop_fac, &
                &                  int_spl(:, i_LL), n_spl, &
                &                  radI_splines(:, i_radI), drop_err)
             end do
          end if
       end do
    end do
    call sync_vector(radI_splines, (n_spl+2)*n_radI)
    call get_max_double(max_drop_err, drop_err)
    if(output_priority .le. OL_norm) then
      write(info_str, "(2X,'| ',A,I6,A,I6,A,I6,A,ES10.2)") &
      & "Spline reduction from", N, " to", n_spl, &
      & " points by (cutting and)", drop_fac, " x coarsening gives error:", &
      & max_drop_err
      call localorb_info(info_str)
    endif
    deallocate(int_spl)

    ! --- Check ranges
    
    allocate(iarange_species(n_species, n_species), stat=info)
    call check_allocation(info, 'iarange_species', func,n_species,n_species)
    iarange_species = 0.d0

    if (.not. have_mp_far) call tighten_iarange()

    ! --- Buffer storage

    allocate(tmp_radI(0:tb_max_L, max_n_basbas_fnLsp, max_n_basbas_fnLsp), &
    &        stat=info)
    call check_allocation(info, 'radI', func,tb_max_L+1,max_n_basbas_fnLsp,max_n_basbas_fnLsp)
    allocate(tmp_fullI(-tb_max_L:tb_max_L, -tb_max_L:tb_max_L, &
    &                  max_n_basbas_fnLsp, max_n_basbas_fnLsp), stat=info)
    call check_allocation(info, 'fullI', func,2*tb_max_L+1,2*tb_max_L+1,max_n_basbas_fnLsp,max_n_basbas_fnLsp)

    call stop_timer(time_init, .true.)

  end subroutine initialize_tb_auxmat
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/initialize_periodic_tb_auxmat
  !  NAME
  !    initialize_periodic_tb_auxmat
  !  SYNOPSIS

  subroutine initialize_periodic_tb_auxmat(drop_fac, Ewald_gamma)

    !  PURPOSE
    !
    !    Initialize the arrays needed for a tight-binding like auxiliary
    !    Coulomb matrix construction for periodic systems and the bare Coulomb
    !    kernel (for now).
    !
    !  USES

    use grids
    use prodbas
    use synchronize_mpi_basic
    use runtime_choices
    use geometry, only: species, recip_lattice_vector
    use dimensions
    use mpi_tasks, only: aims_stop, check_allocation, myid, n_tasks
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: drop_fac
    real*8, intent(IN) :: Ewald_gamma

    !  INPUTS
    !    o drop_fac -- Factor by which to reduce the density of points
    !                Should probably be something in the order of 1 or 4.
    !    o Ewald_gamma -- Decay coefficient (in bohr**(-2)) of Gaussian
    !                     functions in Ewald construction.
    !              This parameter should not change the final result and
    !              is important only for efficiency.  The best choice crucially
    !              depends if one is interested in V(R) or V(q).
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: N, Nt, max_LL
    real*8 :: lnk0, lnr0, lnrange, lnt0, lntrange
    integer :: i_species, L, i_fnLsp, i_basbas_fn, i_basbas_sp
    integer :: L1, L2
    real*8 :: mp_moment
    character*150 :: info_str
    integer :: i
    real*8, allocatable :: GG_onsite(:), GG_splines(:,:,:,:)
    logical, allocatable :: radI_spl_done(:)
    integer :: i_species_1, i_species_2, i_fnLsp_1, i_fnLsp_2
    integer :: i_bb_fn_1, i_bb_fn_2
    real*8 :: mp_moment_1, mp_moment_2
    integer :: i_LL, LL, i_radI
    integer :: info
    real*8, parameter :: sing_thres = 1d-8
    character(*), parameter :: func = 'initialize_periodic_tb_auxmat'

    ! --- Ordinary initialization

    call initialize_tb_auxmat(drop_fac, OVLP_TYPE_COULOMB)
    N = sbtgrid_N
    lnk0 = sbtgrid_lnk0
    lnr0 = sbtgrid_lnr0
    lnrange = sbtgrid_lnrange
    lntrange = 2*lnrange
    Nt = 2*N
    lnt0 = lnk0 + lnr0
    max_LL = 2*tb_max_L

    Gaussian_gamma = Ewald_gamma

    ! --- Prepare indexing arrays

    allocate(Lsp2mp_bb_fnLsp(0:tb_max_L, n_species), stat=info)
    call check_allocation(info, 'Lsp2mp_bb_fnLsp', func)
    allocate(Lsp2mp_bb_fn(0:tb_max_L, n_species), stat=info)
    call check_allocation(info, 'Lsp2mp_bb_fn', func)
    allocate(Lsp2mp_bb_sp(0:tb_max_L, n_species), stat=info)
    call check_allocation(info, 'Lsp2mp_bb_sp', func)
    allocate(Lsp2mp_bb_moment(0:tb_max_L, n_species), stat=info)
    call check_allocation(info, 'Lsp2mp_bb_moment', func)

    Lsp2mp_bb_fnLsp = 0
    Lsp2mp_bb_fn = 0
    Lsp2mp_bb_sp = 0
    Lsp2mp_bb_moment = 0.d0

    do i_species = 1, n_species
       do L = 0, tb_max_L
          do i_fnLsp = 1, Lsp2n_basbas_fnLsp(L, i_species)
             i_basbas_fn = Lsp2basbas_fn(i_fnLsp, L, i_species)
             mp_moment = multipole_basbas_fn(i_basbas_fn)
             if (abs(mp_moment) > 1d-8) then
                if (Lsp2mp_bb_fnLsp(L, i_species) /= 0) then
                   call aims_stop('Multiple nonzero multipoles in L-channel', &
                   &              func)
                end if
                i_basbas_sp = Lsp2basbas_sp(i_fnLsp, L, i_species)
                Lsp2mp_bb_fnLsp(L, i_species) = i_fnLsp
                Lsp2mp_bb_fn(L, i_species) = i_basbas_fn
                Lsp2mp_bb_sp(L, i_species) = i_basbas_sp
                Lsp2mp_bb_moment(L, i_species) = mp_moment
             end if
          end do
       end do
    end do

    ! --- Monkey patch radI

    allocate(GG_onsite(0:tb_max_L), stat=info)
    call check_allocation(info, 'GG_onsite', func)
    allocate(GG_splines(0:n_spl+1, 0:tb_max_L, 0:tb_max_L, 0:tb_max_L), &
    &        stat=info)
    call check_allocation(info, 'GG_splines', func)

    call get_gaussian_interaction_splines(Gaussian_gamma, GG_onsite, GG_splines)

    ! All radial splines have to be fudged by V(R) -> V(R) - p1*p2*V_Gauss(R).
    allocate(radI_spl_done(n_radI), stat=info)
    call check_allocation(info, 'radI_spl_done', func)
    radI_spl_done = .false.
    do i_species_1 = 1, n_species
       do L1 = 0, tb_max_L
          i_fnLsp_1 = Lsp2mp_bb_fnLsp(L1, i_species_1)
          if (i_fnLsp_1 <= 0) cycle
          mp_moment_1 = Lsp2mp_bb_moment(L1, i_species_1)
          i_bb_fn_1 = Lsp2mp_bb_fn(L1, i_species_1)

          radI_onsites(i_fnLsp_1, i_fnLsp_1, L1, i_species_1) = &
          & radI_onsites(i_fnLsp_1, i_fnLsp_1, L1, i_species_1) &
          & - mp_moment_1**2 * GG_onsite(L1)

          do i_species_2 = 1, n_species
             do L2 = 0, tb_max_L
                i_fnLsp_2 = Lsp2mp_bb_fnLsp(L2, i_species_2)
                if (i_fnLsp_2 <= 0) cycle
                mp_moment_2 = Lsp2mp_bb_moment(L2, i_species_2)
                i_bb_fn_2 = Lsp2mp_bb_fn(L2, i_species_2)

                do i_LL = 0, min(L1, L2)
                   LL = abs(L1 - L2) + 2*i_LL
                   i_radI = Lfnfn2radI(i_LL, i_bb_fn_1, i_bb_fn_2)
                   if (.not. radI_spl_done(i_radI)) then
                      ! Avoid double counting
                      radI_splines(:, i_radI) &
                      & = radI_splines(:, i_radI) &
                      &   - mp_moment_1*mp_moment_2*GG_splines(:, i_LL, L1, L2)
                      radI_spl_done(i_radI) = .true.
                   end if
                end do
             end do
          end do
       end do
    end do

    ! Far field should be removed, now.
    have_mp_far = .false.
    have_Ewald = .true.

!    have_mp_far = .true.
    deallocate(radI_spl_done, GG_onsite, GG_splines)

    ! --- Check iaranges (they should now match spline extents)

    call tighten_iarange()

    ! --- Get singularity_lifting_chi

    call get_singularity_lifting_chi(recip_lattice_vector, n_k_points_xyz, &
    &                                1.d0, sing_thres, singularity_lifting_chi)


  end subroutine initialize_periodic_tb_auxmat
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/deallocate_tb_auxmat
  !  NAME
  !    deallocate_tb_auxmat
  !  SYNOPSIS

  subroutine deallocate_tb_auxmat()

    !  PURPOSE
    !
    !    Free all arrays allocated by initialize_tb_auxmat().
    !
    !  USES

    use localorb_io
    use timing_core, only: output_timeheader, output_timer
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character*150 :: info_str
    character(*), parameter :: func = 'deallocate_tb_auxmat'


    if(output_priority .le. OL_norm) then
      write(info_str, "('Atomic logSBT for 2-center ',A,' matrix')") &
      & trim(ovlp_type_name)
      call output_timeheader('2X', info_str)
      call output_timer('TB initialization', time_init(3:4))
      call output_timer('TB usage', time_rest(3:4))
      call localorb_info('')
    endif

    n_spl = 0
    n_radI = 0
    if (allocated(iarange)) deallocate(iarange)
    if (allocated(iarange_species)) deallocate(iarange_species)
    if (allocated(multipoles)) deallocate(multipoles)
    if (allocated(Lfnfn2radI)) deallocate(Lfnfn2radI)
    if (allocated(radI_onsites)) deallocate(radI_onsites)
    if (allocated(radI_splines)) deallocate(radI_splines)
    if (allocated(tmp_radI)) deallocate(tmp_radI)
    if (allocated(tmp_fullI)) deallocate(tmp_fullI)
    if (allocated(Lsp2mp_bb_fnLsp)) deallocate(Lsp2mp_bb_fnLsp)
    if (allocated(Lsp2mp_bb_fn)) deallocate(Lsp2mp_bb_fn)
    if (allocated(Lsp2mp_bb_sp)) deallocate(Lsp2mp_bb_sp)
    if (allocated(Lsp2mp_bb_moment)) deallocate(Lsp2mp_bb_moment)

  end subroutine deallocate_tb_auxmat
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/calculate_tb_auxmat
  !  NAME
  !    calculate_tb_auxmat
  !  SYNOPSIS

  subroutine calculate_tb_auxmat(i_species_1, i_species_2, Rvec, auxmat)

    !  PURPOSE
    !
    !    Convenience wrapper around calculate_tb_auxmat_and_add().
    !
    !    Calculate, from the splined integrals prepared by
    !    initialize_tb_auxmat(), the auxiliary basis integrals between two
    !    atoms separated by a vector Rvec.
    !
    !  USES

    use bspline, only: val_bspline
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2
    real*8, intent(IN) :: Rvec(3)
    real*8, intent(OUT) :: auxmat(:,:)

    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Rvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !  OUTPUTS
    !    o auxmat -- Atomic part of the full 2-center auxiliary matrix.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'calculate_tb_auxmat'

    auxmat = 0.d0
    call calculate_tb_auxmat_and_add(i_species_1, i_species_2, Rvec, 1, size(auxmat,1), 1, size(auxmat,2), auxmat)

  end subroutine calculate_tb_auxmat
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/calculate_tb_auxmat_and_add
  !  NAME
  !    calculate_tb_auxmat_and_add
  !  SYNOPSIS

  subroutine calculate_tb_auxmat_and_add(i_species_1, i_species_2, Rvec, lb_1, ub_1, lb_2, ub_2, &
  &                                      auxmat, c_phases, c_auxmats)

    !  PURPOSE
    !
    !    Do the same as calculate_tb_auxmat() but add to auxmat instead of
    !    resetting it (this should be more efficient in the far field case).
    !
    !  USES

    use bspline, only: val_bspline
    use constants, only: pi
    use prodbas
    use sbt_overlap
    use timing_core, only: start_timer, stop_timer
    use mpi_tasks, only: aims_stop
    use sbt_overlap_tb, only: mult_angular_integral
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2
    real*8, intent(IN) :: Rvec(3)
    integer:: lb_1, ub_1, lb_2, ub_2
    real*8, intent(INOUT) :: auxmat(:,:)
    complex*16, intent(IN), optional :: c_phases(:)
    complex*16, intent(OUT), optional :: c_auxmats(:,:,:) ! (size(c_phases),:,:)

    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Rvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !    o c_phases (optional) -- Phase factors
    !  OUTPUTS
    !    o auxmat -- Atomic part of the full 2-center auxiliary matrix.
    !    o c_auxmats (optional) -- Replaces (real-valued) auxmat,
    !                              accumulate with different phases.
    !                              Please note that the phase index is first!
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8, parameter :: prefac = (2*pi)**1.5d0
    logical :: is_onsite
    real*8 :: Rabs, lnR, r_i_logR, val(1)
    integer :: L, L1, L2, n_fnL1, n_fnL2, i_bbfn_1, i_bbfn_2, i_radI
    integer :: i_fnL1, i_fnL2, i_LL, i_m1, i_m2
    integer :: i_spbbc_1, i_spbbc_2, i_spbb_1, i_spbb_2
    integer :: i_species
    real*8 :: r_i_lnR, mp_pre, mp_prod, mp_prod_RLL
    logical :: need_LL(0:tb_max_L)
    logical :: is_nonzero(max_n_basbas_fnLsp, max_n_basbas_fnLsp)
    character(*), parameter :: func = 'calculate_tb_auxmat_and_add'
    character*150 :: info_str
    logical :: needed
!call perfon('ctbam')
    if (.not. allocated(iarange)) call aims_stop('Module not initialized', func)
    call start_timer(time_rest)

    Rabs = sqrt(sum(Rvec**2))
    is_onsite = (Rabs < 1d-8)
    if (is_onsite) then
       if (i_species_1 /= i_species_2) then
          call aims_stop('Onsite element with differing species', func)
       end if
       i_species = i_species_1
    else
       r_i_lnR = logsbt_r2ilog(Rabs, spl_ln1, spl_dln)
       if (r_i_lnR < 1) then 
          ! This can only really happen in periodic systems, when two
          ! atoms are suddenly mapped onto one another. This _should_
          ! never happen, of course. Need good debug info to find out
          ! what went wrong if this place is reached.
          write (info_str,'(A,3(E16.8,1X),A)') &
            'This routine was called with a too small interatomic separation, ', Rvec, ' a.u.'
          call aims_stop(info_str, func)
       end if
       i_species = 0   ! Satisfy compiler
    end if

    L1_LOOP: do L1 = 0, tb_max_L
       L2_LOOP: do L2 = 0, tb_max_L
          tmp_radI = 0.d0
          n_fnL1 = Lsp2n_basbas_fnLsp(L1, i_species_1)
          n_fnL2 = Lsp2n_basbas_fnLsp(L2, i_species_2)

          if (is_onsite) then
             if (L1 /= L2) cycle  ! Short circuit
             L = L1
          else
             L = -1   ! Satisfy compiler
          end if

          !check if this combination is needed for the index range computed by this task
          needed=.false.
          do i_fnL1 = 1, n_fnL1
             i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
             do i_m1 = -L1, L1
                i_spbb_1 = i_spbbc_1 + i_m1
                if ((i_spbb_1.lt.lb_1).or.(i_spbb_1.gt.ub_1)) cycle
                do i_fnL2 = 1, n_fnL2
                   i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)
                   do i_m2 = -L2, L2
                      i_spbb_2 = i_spbbc_2 + i_m2
                      if ((i_spbb_2.lt.lb_2).or.(i_spbb_2.gt.ub_2)) cycle
                      needed=.true.
                   end do
                end do
             end do
          end do
          if (.not.needed) cycle

          ! --- Retrieve radial integrals

          is_nonzero = .false.
          need_LL = .false.
          FN1_LOOP: do i_fnL1 = 1, n_fnL1
             i_bbfn_1 = Lsp2basbas_fn(i_fnL1, L1, i_species_1)
             FN2_LOOP: do i_fnL2 = 1, n_fnL2
                i_bbfn_2 = Lsp2basbas_fn(i_fnL2, L2, i_species_2)

                IARANGE_COND: if (is_onsite) then
                   is_nonzero(i_fnL1, i_fnL2) = .true.
                   i_LL = 0
                   need_LL(i_LL) = .true.
                   tmp_radI(i_LL, i_fnL1, i_fnL2) &
                   & = radI_onsites(i_fnL1, i_fnL2, L, i_species)
                else if (Rabs < iarange(i_bbfn_1, i_bbfn_2)) then
                   is_nonzero(i_fnL1, i_fnL2) = .true.
                   if (r_i_lnR > n_spl) then
                      call aims_stop('Large separation', func)
                   end if
                   do i_LL = 0, min(L1, L2)
                      need_LL(i_LL) = .true. ! abs(val(1)) > 1d-15
                      i_radI = Lfnfn2radI(i_LL, i_bbfn_1, i_bbfn_2)
                      if (i_radI < 1) call aims_stop('i_radI<1', func)
                      call val_bspline(n_spl, 1, val, &
                      &                r_i_lnR, radI_splines(:, i_radI))
                      tmp_radI(i_LL, i_fnL1, i_fnL2) = val(1)
                   end do
!                   if(i_fnL1 .eq. 1 .and. i_fnL2 .eq. 1 .and. L1==1 .and. L2==1) then
!                     write(use_unit,'(f16.8,2I4,f16.8)') Rabs, i_fnL1, i_fnL2, tmp_radI(0, i_fnL1, i_fnL2)
!                   endif
                else if (have_mp_far) then
                   mp_prod = multipoles(i_bbfn_1) * multipoles(i_bbfn_2)
                   mp_prod_RLL = mp_prod / Rabs**(L1+L2+1)
                   if (abs(mp_prod_RLL) > 1d-18) then
                      is_nonzero(i_fnL1, i_fnL2) = .true.
                      i_LL = min(L1, L2)
                      need_LL(i_LL) = .true.
                      mp_pre = multipole_prefac(L1, L2)
                      tmp_radI(i_LL, i_fnL1, i_fnL2) = & 
                      & mp_pre * mp_prod_RLL
                   end if
                end if IARANGE_COND
             end do FN2_LOOP
          end do FN1_LOOP

          !check if tmp_fullI has to be computed
          needed=.false.
          do i_fnL1 = 1, n_fnL1
             i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
             do i_fnL2 = 1, n_fnL2
                i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)
                if (is_nonzero(i_fnL1, i_fnL2)) then
                   needed=.true.
                end if
             end do
          end do
          if (.not.needed) cycle

          ! --- Combine with angular integrations
!call perfon('mai')
          call mult_angular_integral(Rvec, prefac, &
          &                          n_fnL1, L1, max_n_basbas_fnLsp, tb_max_L, &
          &                          n_fnL2, L2, max_n_basbas_fnLsp, tb_max_L, &
          &                          tb_max_L, need_LL, tmp_radI, tmp_fullI)
!call perfoff
          ! --- Put results

          if (present(c_phases)) then
             if (.not. present(c_auxmats)) then
                call aims_stop('c_auxmats not present', func)
             else if (size(c_phases) /= size(c_auxmats, 1)) then
                call aims_stop('c_auxmats shape mismatch', func)
             end if
             do i_fnL1 = 1, n_fnL1
                i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
                do i_fnL2 = 1, n_fnL2
                   i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)
                   if (is_nonzero(i_fnL1, i_fnL2)) then
                      do i_m1 = -L1, L1
                         i_spbb_1 = i_spbbc_1 + i_m1
                         if ((i_spbb_1.lt.lb_1).or.(i_spbb_1.gt.ub_1)) cycle                           
                         do i_m2 = -L2, L2
                            i_spbb_2 = i_spbbc_2 + i_m2
                            if ((i_spbb_2.lt.lb_2).or.(i_spbb_2.gt.ub_2)) cycle
                            c_auxmats(:, i_spbb_1, i_spbb_2) &
                            & = c_auxmats(:, i_spbb_1, i_spbb_2) &
                            & + c_phases * tmp_fullI(i_m1, i_m2, i_fnL1, i_fnL2)
                         end do
                      end do
                   end if
                end do
             end do
          else
             do i_fnL1 = 1, n_fnL1
                i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
                do i_fnL2 = 1, n_fnL2
                   i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)
                   if (is_nonzero(i_fnL1, i_fnL2)) then
                      do i_m1 = -L1, L1
                         i_spbb_1 = i_spbbc_1 + i_m1
                         do i_m2 = -L2, L2
                            i_spbb_2 = i_spbbc_2 + i_m2
                            auxmat(i_spbb_1, i_spbb_2) &
                            & = auxmat(i_spbb_1, i_spbb_2) &
                            & + tmp_fullI(i_m1, i_m2, i_fnL1, i_fnL2)
                         end do
                      end do
                   end if
                end do
             end do
          end if

       end do L2_LOOP
    end do L1_LOOP

    call stop_timer(time_rest, .true.)
!call perfoff
  end subroutine calculate_tb_auxmat_and_add
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/fast_calculate_tb_auxmat
  !  NAME
  !    fast_calculate_tb_auxmat
  !  SYNOPSIS

  subroutine fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec, auxmat, d_auxmat, calc_deriv)

    !  PURPOSE
    !
    !    Do the same as calculate_tb_auxmat() but tries to do it in a faster way by
    !    - inlining the code of mult_angular_integral
    !    - using fast_sum_triple_Y_YLM_real (instead of sum_triple_Y_YLM_real)
    !
    !    Optionally calculates also the derivatives with respect to Rvec.
    !
    !  USES

    use bspline, only: val_bspline
    use constants, only: pi
    use logsbt, only: logsbt_r2ilog
    use mpi_tasks, only: aims_stop
    use prodbas, only: max_n_basbas_fnlsp, lsp2n_basbas_fnlsp, lsp2basbas_fn, &
        lsp2basbas_sp
    use sbt_overlap, only:  multipole_prefac
    use timing_core, only: start_timer, stop_timer
    use triple_y, only: calc_yrlm_for_triple_y, fast_sum_triple_y_ylm_real
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2
    real*8, intent(IN) :: Rvec(3)
    real*8, intent(INOUT) :: auxmat(:,:)
    real*8, intent(INOUT) :: d_auxmat(:,:,:)
    logical, intent(IN) :: calc_deriv

    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Rvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !    o calc_deriv -- flag if the derivatives should also be calculated
    !  OUTPUTS
    !    o auxmat -- Atomic part of the full 2-center auxiliary matrix.
    !    o d_auxmat -- Derivatives of auxmat, not touched if calc_deriv = FALSE
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8, parameter :: prefac = (2*pi)**1.5d0
    logical :: is_onsite
    real*8 :: Rabs, lnR, val(1)
    integer :: L, L1, L2, n_fnL1, n_fnL2, i_bbfn_1, i_bbfn_2, i_radI
    integer :: i_fnL1, i_fnL2, i_LL, i_m1, i_m2
    integer :: i_spbbc_1, i_spbbc_2, i_spbb_1, i_spbb_2
    integer :: i_species
    real*8 :: r_i_lnR, d_r_i_lnR(3), mp_pre, mp_prod, mp_prod_RLL
    logical :: need_LL(0:tb_max_L)
    logical :: is_nonzero(max_n_basbas_fnLsp, max_n_basbas_fnLsp)
    character(*), parameter :: func = 'fast_calculate_tb_auxmat'
    real*8, allocatable :: cyl(:,:,:), d_cyl(:,:,:,:)
    real*8, allocatable :: cyl_sum(:,:), d_cyl_sum(:,:,:)
    integer LL, ii
    real*8 :: thisfac, this_radI, sig_prefac, f
    character*150 :: info_str
    logical is_nonz

    if (.not. allocated(iarange)) call aims_stop('Module not initialized', func)
    call start_timer(time_rest)

    auxmat = 0.d0
    if(calc_deriv) d_auxmat = 0.d0

    Rabs = sqrt(sum(Rvec**2))
    is_onsite = (Rabs < 1d-8)
    if (is_onsite) then
       if (i_species_1 /= i_species_2) then
          call aims_stop('Onsite element with differing species', func)
       end if
       i_species = i_species_1
    else
       r_i_lnR = logsbt_r2ilog(Rabs, spl_ln1, spl_dln)
       if (r_i_lnR < 1) then
          ! This can only really happen in periodic systems, when two
          ! atoms are suddenly mapped onto one another. This _should_
          ! never happen, of course. Need good debug info to find out
          ! what went wrong if this place is reached.
          !       write(use_unit,*) "Rabs: ", Rabs
          !       write(use_unit,*) "spl_ln1, spl_dln: ", spl_ln1, spl_dln
          !       write(use_unit,*) "r_i_lnR ", r_i_lnR
          write (info_str,'(A,3(E16.8,1X),A)') &
            'This routine was called with a too small interatomic separation, ', Rvec, 'a.u.'
          call aims_stop(info_str, func)
       end if

       ! Get the derivative of r_i_lnR with respect to Rvec:
       ! (d r_i_lnR)/(d Rvec) = (d r_i_lnR) / (d Rabs) * (d Rabs) / (d Rvec)
       ! (d Rabs) / (d Rvec) = Rvec(i)/Rabs
       ! logsbt_r2ilog simply deliveres (log(r) - ln1) / dln + 1.d0
       ! so its derivative is 1/(r*dln)
       d_r_i_lnR(1) = Rvec(1)/(Rabs*Rabs*spl_dln)
       d_r_i_lnR(2) = Rvec(2)/(Rabs*Rabs*spl_dln)
       d_r_i_lnR(3) = Rvec(3)/(Rabs*Rabs*spl_dln)
       i_species = 0   ! Satisfy compiler
    end if

    call calc_YrLM_for_triple_Y(Rvec, tb_max_L)

    L1_LOOP: do L1 = 0, tb_max_L
       L2_LOOP: do L2 = 0, tb_max_L
          tmp_radI = 0.d0
          n_fnL1 = Lsp2n_basbas_fnLsp(L1, i_species_1)
          n_fnL2 = Lsp2n_basbas_fnLsp(L2, i_species_2)

          if (is_onsite) then
             if (L1 /= L2) cycle  ! Short circuit
             L = L1
          else
             L = -1   ! Satisfy compiler
          end if

          sig_prefac = prefac * (-1)**((L1-L2-abs(L1-L2))/2)

          allocate(cyl(-L1:L1,-L2:L2,0:tb_max_L))
          allocate(d_cyl(-L1:L1,-L2:L2,3,0:tb_max_L))
          allocate(cyl_sum(-L1:L1,-L2:L2))
          allocate(d_cyl_sum(-L1:L1,-L2:L2,3))

          ! --- Retrieve radial integrals

          need_LL = .true.
          FN1_LOOP: do i_fnL1 = 1, n_fnL1
             i_bbfn_1  = Lsp2basbas_fn(i_fnL1, L1, i_species_1)
             i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
             FN2_LOOP: do i_fnL2 = 1, n_fnL2
                i_bbfn_2  = Lsp2basbas_fn(i_fnL2, L2, i_species_2)
                i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)

                is_nonz = .false.

                IARANGE_COND: if (is_onsite) then
                   i_LL = 0
                   if(need_LL(i_LL)) then
                      LL = abs(L1 - L2) + 2*i_LL
                      ! Angular integration, derivative is not needed
                      call fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL),d_cyl(:,:,:,i_LL),.false.)
                      need_LL(i_LL) = .false.
                   endif
                   cyl_sum(:,:) = sig_prefac*radI_onsites(i_fnL1, i_fnL2, L, i_species)*cyl(:,:,i_LL)
                   if(calc_deriv) d_cyl_sum(:,:,:) = 0
                   is_nonz = .true.
                else if (Rabs < iarange(i_bbfn_1, i_bbfn_2)) then
                   if (r_i_lnR > n_spl) then
                      call aims_stop('Large separation', func)
                   end if
                   thisfac = sig_prefac
                   cyl_sum(:,:) = 0
                   if(calc_deriv) d_cyl_sum(:,:,:) = 0
                   do i_LL = 0, min(L1, L2)
                      if(need_LL(i_LL)) then
                         LL = abs(L1 - L2) + 2*i_LL
                         ! Angular integration
                         call fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL),d_cyl(:,:,:,i_LL),calc_deriv)
                         need_LL(i_LL) = .false.
                      endif
                      i_radI = Lfnfn2radI(i_LL, i_bbfn_1, i_bbfn_2)
                      if (i_radI < 1) call aims_stop('i_radI<1', func)
                      call val_bspline(n_spl, 1, val, r_i_lnR, radI_splines(:, i_radI))
                      cyl_sum(:,:) = cyl_sum(:,:) + thisfac*val(1)*cyl(:,:,i_LL)
                      if(calc_deriv) then
                         d_cyl_sum(:,:,:) = d_cyl_sum(:,:,:) + thisfac*val(1)*d_cyl(:,:,:,i_LL)
                         call val_bspline(n_spl, 1, val, r_i_lnR, radI_splines(:, i_radI), 1)
                         d_cyl_sum(:,:,1) = d_cyl_sum(:,:,1) + thisfac*val(1)*d_r_i_lnR(1)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,2) = d_cyl_sum(:,:,2) + thisfac*val(1)*d_r_i_lnR(2)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,3) = d_cyl_sum(:,:,3) + thisfac*val(1)*d_r_i_lnR(3)*cyl(:,:,i_LL)
                      endif
                      thisfac = -thisfac
                      is_nonz = .true.
                   end do
                else if (have_mp_far) then
                   mp_prod = multipoles(i_bbfn_1) * multipoles(i_bbfn_2)
                   mp_prod_RLL = mp_prod / Rabs**(L1+L2+1)
                   if (abs(mp_prod_RLL) > 1d-18) then
                      i_LL = min(L1, L2)
                      if(need_LL(i_LL)) then
                         LL = abs(L1 - L2) + 2*i_LL
                         ! Angular integration
                         call fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL),d_cyl(:,:,:,i_LL),calc_deriv)
                         need_LL(i_LL) = .false.
                      endif
                      mp_pre = multipole_prefac(L1, L2)
                      thisfac = sig_prefac
                      if(mod(i_LL,2)==1) thisfac = -thisfac
                      cyl_sum(:,:) = thisfac*mp_pre*mp_prod_RLL*cyl(:,:,i_LL)
                      if(calc_deriv) then
                         d_cyl_sum(:,:,:) = thisfac*mp_pre*mp_prod_RLL*d_cyl(:,:,:,i_LL)
                         f = -(L1+L2+1)*thisfac*mp_pre*mp_prod_RLL/Rabs**2
                         d_cyl_sum(:,:,1) = d_cyl_sum(:,:,1) + f*Rvec(1)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,2) = d_cyl_sum(:,:,2) + f*Rvec(2)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,3) = d_cyl_sum(:,:,3) + f*Rvec(3)*cyl(:,:,i_LL)
                      endif
                      is_nonz = .true.
                   end if
                end if IARANGE_COND

                if(is_nonz) then
                   auxmat(i_spbbc_1-L1:i_spbbc_1+L1,i_spbbc_2-L2:i_spbbc_2+L2) = &
                     auxmat(i_spbbc_1-L1:i_spbbc_1+L1,i_spbbc_2-L2:i_spbbc_2+L2) + cyl_sum(:,:)
                   if(calc_deriv) then
                     d_auxmat(i_spbbc_1-L1:i_spbbc_1+L1,i_spbbc_2-L2:i_spbbc_2+L2,1:3) = &
                       d_auxmat(i_spbbc_1-L1:i_spbbc_1+L1,i_spbbc_2-L2:i_spbbc_2+L2,1:3) + d_cyl_sum(:,:,:)
                   endif
                endif

             end do FN2_LOOP
          end do FN1_LOOP

          deallocate(cyl)
          deallocate(cyl_sum)
          deallocate(d_cyl)
          deallocate(d_cyl_sum)

       end do L2_LOOP
    end do L1_LOOP

    call stop_timer(time_rest, .true.)

  end subroutine fast_calculate_tb_auxmat


  subroutine my_fast_calculate_tb_auxmat(i_species_1, i_species_2, Rvec, lb1, ub1, lb2, ub2, auxmat, d_auxmat, calc_deriv, add)

    !  PURPOSE
    !
    !    Do the same as calculate_tb_auxmat() but tries to do it in a faster way by
    !    - inlining the code of mult_angular_integral
    !    - using fast_sum_triple_Y_YLM_real (instead of sum_triple_Y_YLM_real)
    !
    !    Optionally calculates also the derivatives with respect to Rvec.
    !
    !  USES

    use constants, only: pi
    use my_triple_Y
    use prodbas
    use sbt_overlap, only: logsbt_r2ilog, multipole_prefac
    use bspline, only: val_bspline
    use timing_core, only: start_timer, stop_timer
    use mpi_tasks, only: aims_stop

    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2, lb1, ub1, lb2, ub2
    real*8, intent(IN) :: Rvec(3)
    real*8, intent(INOUT) :: auxmat(lb1:ub1,lb2:ub2)
    real*8, intent(INOUT) :: d_auxmat(lb1:ub1,lb2:ub2,1:3)
    logical, intent(IN) :: calc_deriv, add

    !  INPUTS
    !    o i_species_1, i_species_2 -- Species of the two atoms involved.
    !    o Rvec -- Essentially coords(:, i_atom_2) - coords(:, i_atom_1).
    !    o calc_deriv -- flag if the derivatives should also be calculated
    !  OUTPUTS
    !    o auxmat -- Atomic part of the full 2-center auxiliary matrix.
    !    o d_auxmat -- Derivatives of auxmat, not touched if calc_deriv = FALSE
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8, parameter :: prefac = (2*pi)**1.5d0
    logical :: is_onsite, preallocated_CrMs
    real*8 :: Rabs, lnR, val(1)
    integer :: L, L1, L2, n_fnL1, n_fnL2, i_bbfn_1, i_bbfn_2, i_radI
    integer :: i_fnL1, i_fnL2, i_LL, i_m1, i_m2
    integer :: i_spbbc_1, i_spbbc_2, i_spbb_1, i_spbb_2
    integer :: i_species
    real*8 :: r_i_lnR, d_r_i_lnR(3), mp_pre, mp_prod, mp_prod_RLL
    logical :: need_LL(0:tb_max_L)
    logical :: is_nonzero(max_n_basbas_fnLsp, max_n_basbas_fnLsp)
    character(*), parameter :: func = 'fast_calculate_tb_auxmat'
    real*8, allocatable :: cyl(:,:,:), d_cyl(:,:,:,:)
    real*8, allocatable :: cyl_sum(:,:), d_cyl_sum(:,:,:)
    integer LL, ii
    real*8 :: thisfac, this_radI, sig_prefac, f
    character*150 :: info_str
    logical is_nonz
    type(YrLM_t) :: YrLM
    integer :: tb_max_L_tmp, lba1, uba1, lba2, uba2

    call perfon('tbaux')
    tb_max_L_tmp = max(tb_max_L,15)

    ! my_fast_sum_triple_Y_YLM_real needs YrLM CrMs_s to be allocated. CrMs_s is constant and therefore allocated/deallocated a few loops further out
    preallocated_CrMs=allocated(CrMs_s)
    if (.not.preallocated_CrMs) call allocate_CrMs_s(tb_max_L, tb_max_L, 2*tb_max_L)
    call allocate_YrLM_type_arrays(YrLM,tb_max_L_tmp)

    if (.not. allocated(iarange)) call aims_stop('Module not initialized', func)
    call start_timer(time_rest)

    if(.not.add) then
       auxmat = 0.d0
       if(calc_deriv) d_auxmat = 0.d0
    end if

    Rabs = sqrt(sum(Rvec**2))
    is_onsite = (Rabs < 1d-8)
    if (is_onsite) then
       if (i_species_1 /= i_species_2) then
          call aims_stop('Onsite element with differing species', func)
       end if
       i_species = i_species_1
    else
       r_i_lnR = logsbt_r2ilog(Rabs, spl_ln1, spl_dln)
       if (r_i_lnR < 1) then
          ! This can only really happen in periodic systems, when two
          ! atoms are suddenly mapped onto one another. This _should_
          ! never happen, of course. Need good debug info to find out
          ! what went wrong if this place is reached.
          !       write(use_unit,*) "Rabs: ", Rabs
          !       write(use_unit,*) "spl_ln1, spl_dln: ", spl_ln1, spl_dln
          !       write(use_unit,*) "r_i_lnR ", r_i_lnR
          write (info_str,'(A,3(E16.8,1X),A)') &
            'This routine was called with a too small interatomic separation, ', Rvec, 'a.u.'
          call aims_stop(info_str, func)
       end if

       ! Get the derivative of r_i_lnR with respect to Rvec:
       ! (d r_i_lnR)/(d Rvec) = (d r_i_lnR) / (d Rabs) * (d Rabs) / (d Rvec)
       ! (d Rabs) / (d Rvec) = Rvec(i)/Rabs
       ! logsbt_r2ilog simply deliveres (log(r) - ln1) / dln + 1.d0
       ! so its derivative is 1/(r*dln)
       d_r_i_lnR(1) = Rvec(1)/(Rabs*Rabs*spl_dln)
       d_r_i_lnR(2) = Rvec(2)/(Rabs*Rabs*spl_dln)
       d_r_i_lnR(3) = Rvec(3)/(Rabs*Rabs*spl_dln)
       i_species = 0   ! Satisfy compiler
    end if

    call my_calc_YrLM_for_triple_Y(Rvec, tb_max_L, YrLM)

    L1_LOOP: do L1 = 0, tb_max_L
       L2_LOOP: do L2 = 0, tb_max_L
          tmp_radI = 0.d0
          n_fnL1 = Lsp2n_basbas_fnLsp(L1, i_species_1)
          n_fnL2 = Lsp2n_basbas_fnLsp(L2, i_species_2)

          if (is_onsite) then
             if (L1 /= L2) cycle  ! Short circuit
             L = L1
          else
             L = -1   ! Satisfy compiler
          end if

          sig_prefac = prefac * (-1)**((L1-L2-abs(L1-L2))/2)

          allocate(cyl(-L1:L1,-L2:L2,0:tb_max_L))
          allocate(d_cyl(-L1:L1,-L2:L2,3,0:tb_max_L))
          allocate(cyl_sum(-L1:L1,-L2:L2))
          allocate(d_cyl_sum(-L1:L1,-L2:L2,3))

          ! --- Retrieve radial integrals
          need_LL = .true.
          FN1_LOOP: do i_fnL1 = 1, n_fnL1
             i_bbfn_1  = Lsp2basbas_fn(i_fnL1, L1, i_species_1)
             i_spbbc_1 = Lsp2basbas_sp(i_fnL1, L1, i_species_1)
             lba1=max(i_spbbc_1-L1,lb1)
             uba1=min(i_spbbc_1+L1,ub1)
             if(lba1.gt.uba1) cycle               !check if cycle contributes to the local part of auxmat
             FN2_LOOP: do i_fnL2 = 1, n_fnL2
                i_bbfn_2  = Lsp2basbas_fn(i_fnL2, L2, i_species_2)
                i_spbbc_2 = Lsp2basbas_sp(i_fnL2, L2, i_species_2)
                lba2=max(i_spbbc_2-L2,lb2)
                uba2=min(i_spbbc_2+L2,ub2)
                if(lba2.gt.uba2) cycle               !check if cycle contributes to the local part of auxmat
                is_nonz = .false.

                IARANGE_COND: if (is_onsite) then
                   i_LL = 0
                   if(need_LL(i_LL)) then
                      LL = abs(L1 - L2) + 2*i_LL
                      ! Angular integration, derivative is not needed
                      call my_fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL),d_cyl(:,:,:,i_LL),.false., YrLM)
                      need_LL(i_LL) = .false.
                   endif
                   cyl_sum(:,:) = sig_prefac*radI_onsites(i_fnL1, i_fnL2, L, i_species)*cyl(:,:,i_LL)
                   if(calc_deriv) d_cyl_sum(:,:,:) = 0
                   is_nonz = .true.
                else if (Rabs < iarange(i_bbfn_1, i_bbfn_2)) then
                   if (r_i_lnR > n_spl) then
                      call aims_stop('Large separation', func)
                   end if
                   thisfac = sig_prefac
                   cyl_sum(:,:) = 0
                   if(calc_deriv) d_cyl_sum(:,:,:) = 0
                   do i_LL = 0, min(L1, L2)
                      if(need_LL(i_LL)) then
                         LL = abs(L1 - L2) + 2*i_LL
                         ! Angular integration
                         call my_fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL), d_cyl(:,:,:,i_LL), calc_deriv, YrLM)

                         need_LL(i_LL) = .false.
                      endif
                      i_radI = Lfnfn2radI(i_LL, i_bbfn_1, i_bbfn_2)
                      if (i_radI < 1) call aims_stop('i_radI<1', func)

                      call val_bspline(n_spl, 1, val, r_i_lnR, radI_splines(:, i_radI))

                      cyl_sum(:,:) = cyl_sum(:,:) + thisfac*val(1)*cyl(:,:,i_LL)
                      if(calc_deriv) then
                         d_cyl_sum(:,:,:) = d_cyl_sum(:,:,:) + thisfac*val(1)*d_cyl(:,:,:,i_LL)
                         call val_bspline(n_spl, 1, val, r_i_lnR, radI_splines(:, i_radI), 1)
                         d_cyl_sum(:,:,1) = d_cyl_sum(:,:,1) + thisfac*val(1)*d_r_i_lnR(1)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,2) = d_cyl_sum(:,:,2) + thisfac*val(1)*d_r_i_lnR(2)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,3) = d_cyl_sum(:,:,3) + thisfac*val(1)*d_r_i_lnR(3)*cyl(:,:,i_LL)
                      endif
                      thisfac = -thisfac
                      is_nonz = .true.
                   end do
                else if (have_mp_far) then
                   mp_prod = multipoles(i_bbfn_1) * multipoles(i_bbfn_2)
                   mp_prod_RLL = mp_prod / Rabs**(L1+L2+1)
                   if (abs(mp_prod_RLL) > 1d-18) then
                      i_LL = min(L1, L2)
                      if(need_LL(i_LL)) then
                         LL = abs(L1 - L2) + 2*i_LL
                         ! Angular integration
                         call my_fast_sum_triple_Y_YLM_real(L1, L2, LL, L1, L2, cyl(:,:,i_LL), d_cyl(:,:,:,i_LL), &
                                                            calc_deriv, YrLM)
                         need_LL(i_LL) = .false.
                      endif
                      mp_pre = multipole_prefac(L1, L2)
                      thisfac = sig_prefac
                      if(mod(i_LL,2)==1) thisfac = -thisfac
                      cyl_sum(:,:) = thisfac*mp_pre*mp_prod_RLL*cyl(:,:,i_LL)
                      if(calc_deriv) then
                         d_cyl_sum(:,:,:) = thisfac*mp_pre*mp_prod_RLL*d_cyl(:,:,:,i_LL)
                         f = -(L1+L2+1)*thisfac*mp_pre*mp_prod_RLL/Rabs**2
                         d_cyl_sum(:,:,1) = d_cyl_sum(:,:,1) + f*Rvec(1)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,2) = d_cyl_sum(:,:,2) + f*Rvec(2)*cyl(:,:,i_LL)
                         d_cyl_sum(:,:,3) = d_cyl_sum(:,:,3) + f*Rvec(3)*cyl(:,:,i_LL)
                      endif
                      is_nonz = .true.
                   end if
                end if IARANGE_COND

                if(is_nonz) then
                   auxmat(lba1:uba1,lba2:uba2) = &
                        auxmat(lba1:uba1,lba2:uba2) + cyl_sum(lba1-i_spbbc_1:uba1-i_spbbc_1,lba2-i_spbbc_2:uba2-i_spbbc_2)                       
                   if(calc_deriv) then
                      d_auxmat(lba1:uba1,lba2:uba2,1:3) = &
                           d_auxmat(lba1:uba1,lba2:uba2,1:3) + d_cyl_sum(lba1-i_spbbc_1:uba1-i_spbbc_1,lba2-i_spbbc_2:uba2-i_spbbc_2,1:3)
                   end if
                endif

             end do FN2_LOOP
          end do FN1_LOOP
          deallocate(cyl)
          deallocate(cyl_sum)
          deallocate(d_cyl)
          deallocate(d_cyl_sum)

       end do L2_LOOP
    end do L1_LOOP
    call deallocate_YrLM_type_arrays(YrLM, tb_max_L_tmp) 
    if (.not.preallocated_CrMs) call deallocate_CrMs_s(tb_max_L, tb_max_L, 2*tb_max_L)
    call perfoff

    call stop_timer(time_rest, .true.)
  end subroutine my_fast_calculate_tb_auxmat

  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/sum_up_auxmat_qvecs_from_realspace
  !  NAME
  !    sum_up_auxmat_qvecs_from_realspace
  !  SYNOPSIS

  subroutine sum_up_auxmat_qvecs_from_realspace(n_qvec, qvec, auxmat, n_row_off, n_row_len, n_col_off, n_col_len)

    !  PURPOSE
    !
    !    Calculate auxiliary matrix for a given set of q-points by a
    !    real-space sum.  This will only work if the radial interaction
    !    splines of this modules are "short-ranged", i.e. if either
    !    (.not. have_mp_far) [bare overlap, HSE, cutCb] or if prepared for an
    !    Ewald-like treatment by initialize_periodic_tb_auxmat().
    !
    !  USES

    use bravais
    use dimensions
    use geometry
    use prodbas
    use my_triple_Y
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_qvec
    real*8, intent(IN) :: qvec(3, n_qvec)
    complex*16, intent(OUT) :: auxmat(:, :, :)
    integer, intent(IN) :: n_row_off, n_row_len, n_col_off, n_col_len

    !  INPUTS
    !    o n_qvec -- Number of qvecs.
    !    o qvec -- Different reciprocal space vectors
    !    The following parameters specify the part of auxmat which has to be set:
    !    o n_row_off -- row offset
    !    o n_row_len -- number of rows
    !    o n_col_off -- column offset
    !    o n_col_len -- number of columns
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_atom_1, i_atom_2, i_species_1, i_species_2
    integer :: bboff_1, bboff_2, n_spbb_1, n_spbb_2
    integer :: n_supercells(3), max_n_sc(3), a1, a2, a3, i_qvec, i_bb_1, i_bb_2
    real*8 :: Dvec(3), Rvec(3), Cvec(3), dummy2d(1,1), dummy3d(1,1,1), phi, Rmax, Cmax
    complex*16 :: c_phases(n_qvec)
    complex*16, allocatable :: c_auxmats(:,:,:)
    real*8, allocatable :: r_auxmat(:,:)
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'sum_up_auxmat_qvecs_from_realspace'
    integer:: lb_1, ub_1, lb_2, ub_2, i1, i2

    if (have_mp_far) then
       ! For the bare Coulomb potential, the radial interaction splines
       ! have to be fudged according to the Ewald procedure and have_mp_far
       ! has to be set to .false., afterwards.
       call aims_stop('Cannot use realspace sum with multipoles', func)
    end if

    !initialize CrMs_s of my_triple_Y so that it can be reused in the inner loop
    call allocate_CrMs_s(tb_max_L, tb_max_L, 2*tb_max_L)
    max_n_sc = 0
    do i_atom_2 = 1, n_atoms
       i_species_2 = species(i_atom_2)
       bboff_2 = atom2basbas_off(i_atom_2)
       n_spbb_2 = sp2n_basbas_sp(i_species_2)

       !check if this atom contributes to the local part of auxmat
       if(bboff_2+n_spbb_2 <= n_col_off .or. bboff_2 >= n_col_off+n_col_len) cycle

       lb_2= max(1,n_col_off-bboff_2+1) 
       ub_2= min(n_col_off+n_col_len-bboff_2,n_spbb_2)

       do i_atom_1 = 1, n_atoms    ! For now, ignore symmetry of Auxmat matrix
          i_species_1 = species(i_atom_1)
          bboff_1 = atom2basbas_off(i_atom_1)
          n_spbb_1 = sp2n_basbas_sp(i_species_1)

          !check if this atom contributes to the local part of auxmat
          if(bboff_1+n_spbb_1 <= n_row_off .or. bboff_1 >= n_row_off+n_row_len) cycle

          lb_1= max(1,n_row_off-bboff_1+1) 
          ub_1= min(n_row_off+n_row_len-bboff_1,n_spbb_1)

          allocate(r_auxmat(lb_1:ub_1, lb_2:ub_2), stat=info)
          allocate(c_auxmats(lb_1:ub_1, lb_2:ub_2, n_qvec), stat=info)

          c_auxmats = 0.d0

          Dvec = coords(:, i_atom_2) - coords(:, i_atom_1)
          Rmax = iarange_species(i_species_1, i_species_2) 
          Cmax = Rmax + sqrt(sum(Dvec**2))
          n_supercells = 0   ! Important for n_periodic < 3.
!          Cmax=200.d0

          call get_n_supercells(n_periodic, lattice_vector, Cmax, n_supercells)
          max_n_sc = max(max_n_sc, n_supercells)

          call perfon('ctbaloop')

          do a1 = -n_supercells(1), n_supercells(1)
             do a2 = -n_supercells(2), n_supercells(2)
                do a3 = -n_supercells(3), n_supercells(3)
                   Cvec = matmul(lattice_vector, (/a1, a2, a3/))
                   Rvec = Dvec + Cvec
                   if (sum(Rvec**2) > Rmax**2) cycle
                   do i_qvec = 1, n_qvec
                      phi = dot_product(qvec(:, i_qvec), Cvec)
                      c_phases(i_qvec) = cmplx(cos(phi), sin(phi), kind(0.d0))
                   end do

                   call my_fast_calculate_tb_auxmat( &
                   &      i_species_1, i_species_2, Rvec, lb_1, ub_1, lb_2, ub_2, &
                   &      r_auxmat, dummy3d, .false., .false.)                   

                   do i_qvec = 1, n_qvec
                      c_auxmats(:,:,i_qvec) = c_auxmats(:,:,i_qvec) + r_auxmat(:,:) * c_phases(i_qvec)
                   end do
                end do
!	      write(use_unit,'(3I4,5f18.8)') a1, a2, a3, sqrt(sum(Rvec**2)), c_auxmats(1,max_n_basbas_sp-1,max_n_basbas_sp-1), c_auxmats(1,max_n_basbas_sp,max_n_basbas_sp)
             end do
!	write(use_unit,'(3I4,2f18.8)') a1, a2, a3, c_auxmats(1,max_n_basbas_sp-1,max_n_basbas_sp-1)
          end do
          call perfoff

          do i_qvec = 1, n_qvec
             do i_bb_2 = lb_2, ub_2
                do i_bb_1 = lb_1, ub_1
                   auxmat(bboff_1+i_bb_1-n_row_off, bboff_2+i_bb_2-n_col_off, i_qvec) &
                   & = c_auxmats(i_bb_1, i_bb_2,i_qvec)
                enddo
             enddo
          end do
          deallocate(c_auxmats)
          deallocate(r_auxmat)
       end do
    end do
    call deallocate_CrMs_s(tb_max_L, tb_max_L, 2*tb_max_L)
  !  write(info_str, "(2X,A,': n_supercells:',3I5)") trim(func), max_n_sc
  !  call localorb_info(info_str)

  end subroutine sum_up_auxmat_qvecs_from_realspace
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/add_long_range_from_qspace
  !  NAME
  !    add_long_range_from_qspace
  !  SYNOPSIS

  subroutine add_long_range_from_qspace(n_qvec, qvec, auxmat, n_row_off, n_row_len, n_col_off, n_col_len)

    !  PURPOSE
    !
    !    Add the periodic long-range correction to auxmat.
    !
    !  USES

    use bravais
    use dimensions
    use geometry
    use prodbas
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_qvec
    real*8, intent(IN) :: qvec(3, n_qvec)
    complex*16, intent(INOUT) :: auxmat(:,:,:)
    integer, intent(IN) :: n_row_off, n_row_len, n_col_off, n_col_len


    !  INPUTS
    !    o n_qvec -- Number of qvecs.
    !    o qvec -- Different reciprocal space vectors
    !    The following parameters specify the part of auxmat which has to be set:
    !    o n_row_off -- row offset
    !    o n_row_len -- number of rows
    !    o n_col_off -- column offset
    !    o n_col_len -- number of columns
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_atom_1, i_species_1, i_fnLsp_1, i_bb_sp_1, bboff_1, i_basbas_1
    integer :: i_atom_2, i_species_2, i_fnLsp_2, i_bb_sp_2, bboff_2, i_basbas_2
    integer :: L1, M1, i_LM1
    integer :: L2, M2, i_LM2
    integer :: i_qvec, n_Dvec_all, i_Dvec_all, n_Dvec, i_Dvec
    real*8 :: mp_moment_1, mp_moment_2
    integer :: max_n_LM
    real*8, allocatable :: Dvec(:,:)
    complex*16, allocatable :: Vq(:,:,:)
    complex*16, allocatable :: tmpmat(:,:,:)
    integer, allocatable :: need_Dvec(:)
    integer :: info
    character(*), parameter :: func = 'add_long_range_from_qspace'

    ! --- Init

    max_n_LM = (tb_max_L+1)**2
    n_Dvec_all = n_atoms**2
    allocate(need_Dvec(n_Dvec_all), stat=info)
    call check_allocation(info, 'need_Dvec', func)

    allocate(tmpmat(n_row_len, n_col_len, n_qvec), stat=info)
    call check_allocation(info, 'tmpmat', func)
    

    ! --- Dvec

    !determine which Dvecs are actually needed for this task
    need_Dvec=0
    do i_atom_1 = 1, n_atoms
       bboff_1 = atom2basbas_off(i_atom_1)
       i_species_1 = species(i_atom_1)
       do i_atom_2 = 1, n_atoms
          bboff_2 = atom2basbas_off(i_atom_2)
          i_species_2 = species(i_atom_2)
          i_Dvec_all = (i_atom_2 - 1)*n_atoms + (i_atom_1 - 1)*1 + 1
          do L1 = 0, tb_max_L
             i_fnLsp_1 = Lsp2mp_bb_fnLsp(L1, i_species_1)
             if (i_fnLsp_1 <= 0) cycle
             i_bb_sp_1 = Lsp2mp_bb_sp(L1, i_species_1)          
   
             do L2 = 0, tb_max_L
                i_fnLsp_2 = Lsp2mp_bb_fnLsp(L2, i_species_2)
                if (i_fnLsp_2 <= 0) cycle
                i_bb_sp_2 = Lsp2mp_bb_sp(L2, i_species_2)
                do M1 = -L1, L1
                   i_basbas_1 = bboff_1 + i_bb_sp_1 + M1
                   if(i_basbas_1 <= n_row_off .or. i_basbas_1 > n_row_off+n_row_len) cycle
                   do M2 = -L2, L2
                      i_basbas_2 = bboff_2 + i_bb_sp_2 + M2
                      if(i_basbas_2 <= n_col_off .or. i_basbas_2 > n_col_off+n_col_len) cycle
                      !if no cycle has occured, the i_Dvec is actually needed
                      need_Dvec(i_Dvec_all)=need_Dvec(i_Dvec_all)+1
                   end do
                end do
             end do
          end do
       end do
    end do

    n_Dvec=0
    do i_Dvec_all=1, n_Dvec_all
       if (need_Dvec(i_Dvec_all).gt.0) n_Dvec=n_Dvec+1
    end do

    allocate(Dvec(3, n_Dvec), stat=info)
    call check_allocation(info, 'Dvec', func)
    allocate(Vq(max_n_LM, max_n_LM, n_Dvec), stat=info)
    call check_allocation(info, 'Vq', func)

    i_Dvec=0
    do i_atom_1 = 1, n_atoms
       do i_atom_2 = 1, n_atoms
          i_Dvec_all = (i_atom_2 - 1)*n_atoms + (i_atom_1 - 1)*1 + 1
          if (need_Dvec(i_Dvec_all).gt.0) then
             i_Dvec=i_Dvec+1
             Dvec(:, i_Dvec) = coords(:, i_atom_2) - coords(:, i_atom_1)
          end if
       end do
    end do

    ! --- Vq
    tmpmat(:,:,:)=(0.d0,0.d0)
    do i_qvec = 1, n_qvec
       call get_gaussian_Vq(qvec(:, i_qvec), recip_lattice_vector, &
       &                    singularity_lifting_chi, n_Dvec, Dvec, &
       &                    Gaussian_gamma, tb_max_L, tb_max_L, max_n_LM, Vq)

       i_Dvec=0
       ATOM_1_LOOP: do i_atom_1 = 1, n_atoms
          bboff_1 = atom2basbas_off(i_atom_1)
          i_species_1 = species(i_atom_1)

          ATOM_2_LOOP: do i_atom_2 = 1, n_atoms
             bboff_2 = atom2basbas_off(i_atom_2)
             i_species_2 = species(i_atom_2)

             i_Dvec_all = (i_atom_2 - 1)*n_atoms + (i_atom_1 - 1)*1 + 1
             if (need_Dvec(i_Dvec_all).gt.0) then
                i_Dvec = i_Dvec+1
             end if

             do L1 = 0, tb_max_L
                i_fnLsp_1 = Lsp2mp_bb_fnLsp(L1, i_species_1)
                if (i_fnLsp_1 <= 0) cycle
                mp_moment_1 = Lsp2mp_bb_moment(L1, i_species_1)
                i_bb_sp_1 = Lsp2mp_bb_sp(L1, i_species_1)
  
                do L2 = 0, tb_max_L
                   i_fnLsp_2 = Lsp2mp_bb_fnLsp(L2, i_species_2)
                   if (i_fnLsp_2 <= 0) cycle
                   mp_moment_2 = Lsp2mp_bb_moment(L2, i_species_2)
                   i_bb_sp_2 = Lsp2mp_bb_sp(L2, i_species_2)

                   M_LOOP: do M1 = -L1, L1
                      i_basbas_1 = bboff_1 + i_bb_sp_1 + M1
                      if(i_basbas_1 <= n_row_off .or. i_basbas_1 > n_row_off+n_row_len) cycle
                      i_LM1 = L1**2 + L1 + M1 + 1
                      do M2 = -L2, L2
                         i_basbas_2 = bboff_2 + i_bb_sp_2 + M2
                         if(i_basbas_2 <= n_col_off .or. i_basbas_2 > n_col_off+n_col_len) cycle
                         i_LM2 = L2**2 + L2 + M2 + 1
!                         auxmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
!                         & = auxmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         tmpmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         & = tmpmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         &   + mp_moment_1 * mp_moment_2 &
                         &     * Vq(i_LM1, i_LM2, i_Dvec)
                      end do
                   end do M_LOOP
                end do
             end do
          end do ATOM_2_LOOP
       end do ATOM_1_LOOP
    end do

    auxmat(:,:,:)=auxmat(:,:,:)+tmpmat(:,:,:)

    deallocate(need_Dvec)
    deallocate(Dvec)
    deallocate(Vq)
    deallocate(tmpmat)

!    do i_basbas_1=1, n_basbas, 1
!       write(use_unit,'(I4,4f18.8)') i_basbas_1, auxmat(i_basbas_1,i_basbas_1,1), tmpmat(i_basbas_1, i_basbas_1,1)
!    enddo

  end subroutine add_long_range_from_qspace
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/get_coulomb_long_range_coeff
  !  NAME
  !    get_coulomb_long_range_coeff
  !  SYNOPSIS

  subroutine get_coulomb_long_range_coeff(n_qvec, qvec, auxmat, n_row_off, n_row_len, n_col_off, n_col_len)

    !  PURPOSE
    !
    !    Add the periodic long-range correction to auxmat.
    !
    !  USES

    use bravais
    use dimensions
    use geometry
    use prodbas
    use mpi_tasks, only: check_allocation
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_qvec
    real*8, intent(IN) :: qvec(3, n_qvec)
    complex*16, intent(OUT) :: auxmat(:, :, :)
    integer, intent(IN) :: n_row_off, n_row_len, n_col_off, n_col_len


    !  INPUTS
    !    o n_qvec -- Number of qvecs.
    !    o qvec -- Different reciprocal space vectors
    !    The following parameters specify the part of auxmat which has to be set:
    !    o n_row_off -- row offset
    !    o n_row_len -- number of rows
    !    o n_col_off -- column offset
    !    o n_col_len -- number of columns
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_atom_1, i_species_1, i_fnLsp_1, i_bb_sp_1, bboff_1, i_basbas_1
    integer :: i_atom_2, i_species_2, i_fnLsp_2, i_bb_sp_2, bboff_2, i_basbas_2
    integer :: L1, M1, i_LM1
    integer :: L2, M2, i_LM2
    integer :: i_qvec, n_Dvec, i_Dvec
    real*8 :: mp_moment_1, mp_moment_2
    integer :: max_n_LM
    real*8, allocatable :: Dvec(:,:)
    complex*16, allocatable :: Vq_coeff(:,:,:)
    complex*16, allocatable :: tmpmat(:,:,:)
    integer :: info
    character(*), parameter :: func = 'get_coulomb_long_range_coeff'

    ! --- Init

    max_n_LM = (tb_max_L+1)**2
    n_Dvec = n_atoms**2
    allocate(Dvec(3, n_Dvec), stat=info)
    call check_allocation(info, 'Dvec', func)
    allocate(Vq_coeff(max_n_LM, max_n_LM, n_Dvec), stat=info)
    call check_allocation(info, 'Vq_coeff', func)
    allocate(tmpmat(n_basbas, n_basbas, n_qvec), stat=info)
    call check_allocation(info, 'tmpmat', func)

    ! --- Dvec

    do i_atom_1 = 1, n_atoms
       do i_atom_2 = 1, n_atoms
          i_Dvec = (i_atom_2 - 1)*n_atoms + (i_atom_1 - 1)*1 + 1
          Dvec(:, i_Dvec) = coords(:, i_atom_2) - coords(:, i_atom_1)
       end do
    end do

    ! --- Vq

    tmpmat(:,:,:)=(0.d0,0.d0)
    do i_qvec = 1, n_qvec
       call get_gaussian_Vq_coeff(qvec(:, i_qvec), recip_lattice_vector, &
       &                    singularity_lifting_chi, n_Dvec, Dvec, &
       &                    Gaussian_gamma, tb_max_L, tb_max_L, max_n_LM, Vq_coeff)
       ATOM_1_LOOP: do i_atom_1 = 1, n_atoms
          bboff_1 = atom2basbas_off(i_atom_1)
          i_species_1 = species(i_atom_1)
          do L1 = 0, tb_max_L
             i_fnLsp_1 = Lsp2mp_bb_fnLsp(L1, i_species_1)
             if (i_fnLsp_1 <= 0) cycle
             mp_moment_1 = Lsp2mp_bb_moment(L1, i_species_1)
             i_bb_sp_1 = Lsp2mp_bb_sp(L1, i_species_1)

             ATOM_2_LOOP: do i_atom_2 = 1, n_atoms
                bboff_2 = atom2basbas_off(i_atom_2)
                i_species_2 = species(i_atom_2)

                i_Dvec = (i_atom_2 - 1)*n_atoms + (i_atom_1 - 1)*1 + 1

                do L2 = 0, tb_max_L
                   i_fnLsp_2 = Lsp2mp_bb_fnLsp(L2, i_species_2)
                   if (i_fnLsp_2 <= 0) cycle
                   mp_moment_2 = Lsp2mp_bb_moment(L2, i_species_2)
                   i_bb_sp_2 = Lsp2mp_bb_sp(L2, i_species_2)

                   M_LOOP: do M1 = -L1, L1
                      i_basbas_1 = bboff_1 + i_bb_sp_1 + M1
                      if(i_basbas_1 <= n_row_off .or. i_basbas_1 > n_row_off+n_row_len) cycle
                      i_LM1 = L1**2 + L1 + M1 + 1
                      do M2 = -L2, L2
                         i_basbas_2 = bboff_2 + i_bb_sp_2 + M2
                         if(i_basbas_2 <= n_col_off .or. i_basbas_2 > n_col_off+n_col_len) cycle
                         i_LM2 = L2**2 + L2 + M2 + 1
!                         auxmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
!                         & = auxmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         tmpmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         & = tmpmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec) &
                         &   + mp_moment_1 * mp_moment_2 &
                         &     * Vq_coeff(i_LM1, i_LM2, i_Dvec)
!                           if(i_basbas_1-n_row_off .eq. i_basbas_2-n_col_off) then
!                             write(use_unit,'(I4,6f18.8)')i_basbas_1-n_row_off , mp_moment_1, mp_moment_2, Vq_coeff(i_LM1, i_LM2, i_Dvec,2), &
!                                    tmpmat(i_basbas_1-n_row_off, i_basbas_2-n_col_off, i_qvec)
!                           endif
                      end do
                   end do M_LOOP
                end do
             end do ATOM_2_LOOP
          end do
       end do ATOM_1_LOOP
    end do
    auxmat(:,:,:)=tmpmat(:,:,:)
!    do i_basbas_1=1, n_basbas, 1
!       write(use_unit,'(I4,4f18.8)') i_basbas_1, auxmat(i_basbas_1,i_basbas_1,1), tmpmat(i_basbas_1, i_basbas_1,1)
!    enddo

  end subroutine get_coulomb_long_range_coeff
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/get_qspace_auxmat
  !  NAME
  !    get_qspace_auxmat
  !  SYNOPSIS

  subroutine get_qspace_auxmat(n_qvec, qvec, auxmat, n_row_off, n_row_len, n_col_off, n_col_len)

    !  PURPOSE
    !
    !    Calculate an auxiliary matrix in q-space.  The overlap type has been
    !    defined at initialization time.  For ovlp_type=OVLP_TYPE_COULOMB,
    !    make sure that initialize_periodic_tb_auxmat() has been used.
    !
    !  USES

    use bravais
    use geometry
    use prodbas
    use mpi_tasks
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_qvec
    real*8, intent(IN) :: qvec(3, n_qvec)
    complex*16, intent(INOUT) :: auxmat(:, :, :)
    integer, intent(IN), optional :: n_row_off, n_row_len, n_col_off, n_col_len

    !  INPUTS
    !    o n_qvec -- Number of qvecs.
    !    o qvec -- Different reciprocal space vectors
    !    The following parameters need only be present if auxmat should only partially set:
    !    o n_row_off -- row offset
    !    o n_row_len -- number of rows
    !    o n_col_off -- column offset
    !    o n_col_len -- number of columns
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'get_qspace_auxmat'
    integer my_row_off, my_row_len, my_col_off, my_col_len
    integer i_row
    real*8:: checkloc, checkglob
    integer:: mpierr
call perfon('tibi')
    if(present(n_row_off)) then
      my_row_off = n_row_off
    else
      my_row_off = 0
    endif
    if(present(n_row_len)) then
      my_row_len = n_row_len
    else
      my_row_len = ubound(auxmat,1)
    endif

    if(present(n_col_off)) then
      my_col_off = n_col_off
    else
      my_col_off = 0
    endif
    if(present(n_col_len)) then
      my_col_len = n_col_len
    else
      my_col_len = ubound(auxmat,2)
    endif

    auxmat=0.
call perfon('tb_suaq')
    call sum_up_auxmat_qvecs_from_realspace(n_qvec, qvec, auxmat, my_row_off, my_row_len, my_col_off, my_col_len)
call perfoff
!    write(use_unit,*)"short-range"
!    do i_row = 1, my_row_len, 1
!      write(use_unit,'(I4,2f20.8)')i_row,  auxmat(i_row,i_row,1)
!    enddo
call perfon('tb_alr')
    if (have_Ewald) call add_long_range_from_qspace(n_qvec, qvec, auxmat, my_row_off, my_row_len, my_col_off, my_col_len)
call perfoff
!    write(use_unit,*)"add long-range"
!    do i_row =1, my_row_len, 1
!      write(use_unit,'(I4,2f20.8)')i_row,  auxmat(i_row,i_row,1)
!    enddo

    checkloc=sum(abs(auxmat))
    call mpi_allreduce(checkloc,checkglob,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpierr)
!    if (myid.eq.0) print*,'TBOUT: ', checkglob
call perfoff
  end subroutine get_qspace_auxmat
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/get_qspace_coulomb_coeff
  !  NAME
  !    get_qspace_coulomb_coeff
  !  SYNOPSIS

  subroutine get_qspace_coulomb_coeff(n_qvec, qvec, auxmat, n_row_off, n_row_len, n_col_off, n_col_len)

    !  PURPOSE
    !
    !    Calculate an auxiliary matrix in q-space.  The overlap type has been
    !    defined at initialization time.  For ovlp_type=OVLP_TYPE_COULOMB,
    !    make sure that initialize_periodic_tb_auxmat() has been used.
    !
    !  USES

    use bravais
    use geometry
    use prodbas
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: n_qvec
    real*8, intent(IN) :: qvec(3, n_qvec)
    complex*16, intent(INOUT) :: auxmat(:, :, :)
    integer, intent(IN), optional :: n_row_off, n_row_len, n_col_off, n_col_len

    !  INPUTS
    !    o n_qvec -- Number of qvecs.
    !    o qvec -- Different reciprocal space vectors
    !    The following parameters need only be present if auxmat should only partially set:
    !    o n_row_off -- row offset
    !    o n_row_len -- number of rows
    !    o n_col_off -- column offset
    !    o n_col_len -- number of columns
    !  OUTPUTS
    !    o auxmat -- overlap or (screened/cut/bare) Coulomb matrix
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'get_qspace_coulomb_coeff'
    integer my_row_off, my_row_len, my_col_off, my_col_len
    integer i_row

    if(present(n_row_off)) then
      my_row_off = n_row_off
    else
      my_row_off = 0
    endif
    if(present(n_row_len)) then
      my_row_len = n_row_len
    else
      my_row_len = ubound(auxmat,1)
    endif

    if(present(n_col_off)) then
      my_col_off = n_col_off
    else
      my_col_off = 0
    endif
    if(present(n_col_len)) then
      my_col_len = n_col_len
    else
      my_col_len = ubound(auxmat,2)
    endif

!    write(use_unit,*)"short-range"
!    do i_row = 1, my_row_len, 1
!      write(use_unit,'(I4,2f20.8)')i_row,  auxmat(i_row,i_row,1)
!    enddo
!    if (have_Ewald) call get_coulomb_long_range_coeff(n_qvec, qvec, auxmat, my_row_off, my_row_len, my_col_off, my_col_len)
     call get_coulomb_long_range_coeff(n_qvec, qvec, auxmat, my_row_off, my_row_len, my_col_off, my_col_len)
!    write(use_unit,*)"add long-range"
!    do i_row =1, my_row_len, 1
!      write(use_unit,'(I4,2f20.8)')i_row,  auxmat(i_row,i_row,1)
!    enddo
  end subroutine get_qspace_coulomb_coeff
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/integrate_auxmat_by_tb
  !  NAME
  !   integrate_auxmat_by_tb
  !  SYNOPSIS

  subroutine integrate_auxmat_by_tb(auxmat, ovlp_type, current_cell_index)

    !  PURPOSE
    !
    !    This procedure is intended to calculate the auxmat interaction
    !    matrix elements between every two auxiliary basis functions.
    !
    !  USES

    use dimensions
    use prodbas
    use geometry, only: lattice_vector, coords, species
    use localorb_io
    implicit none

    !  ARGUMENTS

    real*8, intent(OUT) :: auxmat(n_basbas, n_basbas)
    integer, intent(IN) :: ovlp_type
    integer, intent(IN), optional :: current_cell_index(3)

    !  INPUTS
    !     o ovlp_type -- Either OVLP_TYPE_OVERLAP, OVLP_TYPE_COULOMB,
    !                    or OVLP_TYPE_HSE (use runtime_choices:hse_omega_hf).
    !     o current_cell_index -- Unit cell of the column product basis funcs
    !     o ... global data in prodbas.f90 and friends
    !  OUTPUTS
    !     o auxmat -- real array, the Auxmat interaction/overlap matrix within
    !                 the auxiliary basis
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  SEE ALSO
    !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
    !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
    !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
    !    Computer Physics Communications (2008), submitted.
    !  COPYRIGHT
    !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
    !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
    !   the terms and conditions of the respective license agreement."
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE

    real*8 :: Rvec(3)
    integer :: i_atom_1, i_atom_2, i_species_1, i_species_2
    integer :: bboff_1, bboff_2, n_spbb_1, n_spbb_2
    character*150 :: info_str
    integer :: info
    character(*), parameter :: func = 'integrate_auxmat_by_tb'

    auxmat = 0.d0

    if(output_priority .le. OL_norm) then
      write(info_str, "(2X,A,A,A)") &
      & "Perparing the auxiliary ", trim(OVLP_TYPE_NAMES(ovlp_type)), &
      & " 2-center integration (by TB)."
    endif
    call localorb_info(info_str)
    call localorb_info('')

    call initialize_tb_auxmat(1, ovlp_type)

    ! --- Main loop over atom pairs

    do i_atom_2 = 1, n_atoms
       i_species_2 = species(i_atom_2)
       bboff_2 = atom2basbas_off(i_atom_2)
       n_spbb_2 = sp2n_basbas_sp(i_species_2)

       do i_atom_1 = 1, n_atoms    ! For now, ignore symmetry of Auxmat matrix
          i_species_1 = species(i_atom_1)
          bboff_1 = atom2basbas_off(i_atom_1)
          n_spbb_1 = sp2n_basbas_sp(i_species_1)

          Rvec = coords(:, i_atom_2) - coords(:, i_atom_1)
          if (present(current_cell_index)) then
             Rvec = Rvec + matmul(lattice_vector, current_cell_index)
          end if

          ! Calculate
          call calculate_tb_auxmat(i_species_1, i_species_2, Rvec, &
          &                        auxmat(bboff_1+1:bboff_1+n_spbb_1, &
          &                               bboff_2+1:bboff_2+n_spbb_2))

       end do
    end do

    ! --- Finalize

    call deallocate_tb_auxmat()

  end subroutine integrate_auxmat_by_tb
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/get_gaussian_interaction_splines
  !  NAME
  !    get_gaussian_interaction_splines
  !  SYNOPSIS

  subroutine get_gaussian_interaction_splines(gamma, GG_onsite, GG_splines)

    !  PURPOSE
    !
    !    Get the Coulomb interaction radial parts V_L(R) among Gaussians
    !    normalized to unity multipole moments
    !
    !      G_lm(r) := N_l r^l e^(-gamma r^2/2) Y_lm(r)
    !      N_l = gamma^l+1.5 / [sqrt(pi/2) (2l-1)!!]
    !
    !  USES

    use constants
    use localorb_io
    use runtime_choices
    use bspline, only: get_subspline
    use mpi_tasks, only: check_allocation
    use sbt_overlap
    use sbt_overlap_tb, only: get_radial_integral_onsite, &
        get_radial_integral_spline
    implicit none

    !  ARGUMENTS

    real*8, intent(IN) :: gamma
    real*8, intent(OUT) :: GG_onsite(0:tb_max_L)
    real*8, intent(OUT) :: GG_splines(0:n_spl+1, 0:tb_max_L, &
    &                                 0:tb_max_L, 0:tb_max_L)

    !  INPUTS
    !    o gamma -- Decay coefficient in exp(-gamma*r**2/2).
    !  OUTPUTS
    !    o GG_onsites -- Onsite Coulomb interaction (f_LM(r) | f_LM(r)).
    !    o GG_splines(:, i_LL, L1, L2) -- Spline of the i_LL-th radial part
    !                                     V_LL(R) for L1-L2 interaction.
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: N, Nt, max_LL
    real*8 :: lnk0, lnr0, lnrange, lnt0, lntrange
    real*8, allocatable :: KKbars(:,:), wsave(:)
    integer :: LL, L, i, i_LL, L1, L2
    real*8, allocatable :: Gk(:,:), normalization(:)
    real*8, allocatable :: int_spl(:,:)
    real*8 :: k_cur, exp_val, drop_err
    integer :: info
    real*8, external :: double_factorial
    character*150 :: info_str
    character(*), parameter :: func = 'get_gaussian_interaction_splines'

    ! --- Initialize conveniences & arrays

    N = sbtgrid_N
    lnk0 = sbtgrid_lnk0
    lnr0 = sbtgrid_lnr0
    lnrange = sbtgrid_lnrange
    lntrange = 2*lnrange
    Nt = 2*N
    lnt0 = lnk0 + lnr0
    max_LL = 2*tb_max_L
    allocate(KKbars(Nt, 0:max_LL), wsave(2*Nt+15), stat=info)
    call check_allocation(info, 'KKbars, wsave', func)
    allocate(Gk(N, 0:tb_max_L), stat=info)
    call check_allocation(info, 'Gk', func)
    allocate(normalization(0:tb_max_L), stat=info)
    call check_allocation(info, 'normalization', func)

    ! --- Prepare kernels

    call dffti(Nt, wsave)
    do LL = 0, max_LL
       call logsbt_kernel(Nt, KKbars(:, LL), lnt0, lntrange, LL, &
       &                  power_bias_int, wsave)
    end do

    ! --- Get Gaussian Coulomb radial splines

    ! For simplicity, simple NAO the Gaussians in Fourier space and
    ! use the same algorithms as for the actual basbas NAOs.

    do L = 0, tb_max_L
       normalization(L) = 1.d0 / sqrt(pi/2) / double_factorial(2*L-1)
    end do
    do i = 1, N
       k_cur = logsbt_ilog2r(dble(i), lnk0, lnrange/N)
       exp_val = exp(-k_cur**2 / (2*gamma)) * sqrt(4*pi) * k_cur**0.5d0
       do L = 0, tb_max_L
          Gk(i, L) = normalization(L) * k_cur**L * exp_val
       end do
    end do

    ! --- Get Gaussian interactions

    allocate(int_spl(0:N+1, 0:tb_max_L), stat=info)
    call check_allocation(info, 'int_spl', func)

    drop_err = 0.d0
    do L1 = 0, tb_max_L
       call get_radial_integral_onsite(N, lnr0, lnk0, lnrange, &
       &                               1, L1, Gk(:, L1), &
       &                               1, L1, Gk(:, L1), GG_onsite(L1))
       do L2 = 0, L1
          call get_radial_integral_spline(N, lnr0, lnk0, lnrange, &
          &                               max_LL, KKbars, wsave, &
          &                               1, L1, Gk(:, L1), &
          &                               1, L2, Gk(:, L2), &
          &                               tb_max_L, int_spl)
          do i_LL = 0, min(L1, L2)
             LL = abs(L1 - L2) + 2*i_LL
             call get_subspline(N, 1, i_minR, i_maxR, use_drop_fac, &
             &                  int_spl(:, i_LL), n_spl, &
             &                  GG_splines(:, i_LL, L1, L2), drop_err)
          end do
          if (L1 /= L2) then
             GG_splines(:,:, L2, L1) = GG_splines(:,:, L1, L2)
          end if
       end do
    end do
    if(output_priority .le. OL_norm) then
      write(info_str, "(2X,'| ',A,I6,A,I6,A,I6,A,ES10.2)") &
      & "Gaussian spline reduction from", N, " to", n_spl, &
      & " points by (cutting and)", use_drop_fac, " x coarsening gives error:", &
      & drop_err
    endif

    call localorb_info(info_str)
    deallocate(Gk, int_spl, KKbars, wsave, normalization)

  end subroutine get_gaussian_interaction_splines
  !******
  !----------------------------------------------------------------------------
  !****s* tight_binding_auxmat/tighten_iarange
  !  NAME
  !    tighten_iarange
  !  SYNOPSIS

  subroutine tighten_iarange()

    !  PURPOSE
    !
    !     Tighten iarange and set iarange_species correspondingly.
    !
    !  USES

    use dimensions
    use runtime_choices
    use prodbas
    use debug_output
    use sbt_overlap
    use localorb_io, only: localorb_info
    implicit none

    !  ARGUMENTS
    !    none
    !  INPUTS
    !    none  [radI_splines]
    !  OUTPUTS
    !    none  [iarange, iarange_species]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    integer :: i_bbfn_1, L1, i_species_1
    integer :: i_bbfn_2, L2, i_species_2
    integer :: i_this, i_LL, LL, i_radI, i_iarange, i
    real*8 :: maxerr, thismaxerr, err, this_range
    logical, parameter :: be_verbose = .false.
    logical, parameter :: out_vals = .false.
    real*8 :: vals(n_spl)
    character*150 :: info_str, filename
    character(*), parameter :: func = 'tighten_iarange'

    maxerr = 0.d0
    do i_bbfn_1 = 1, n_basbas_fns
       L1 = basbasfn_l(i_bbfn_1)
       i_species_1 = basbasfn_species(i_bbfn_1)
       do i_bbfn_2 = 1, i_bbfn_1
          L2 = basbasfn_l(i_bbfn_2)
          i_species_2 = basbasfn_species(i_bbfn_2)

          i_this = 0
          thismaxerr = 0.d0
          do i_LL = 0, min(L1, L2)
             LL = abs(L1 - L2) + 2*i_LL
             i_radI = Lfnfn2radI(i_LL, i_bbfn_1, i_bbfn_2)

             i_iarange = logsbt_r2ilog(iarange(i_bbfn_1, i_bbfn_2), &
             &                         spl_ln1, spl_dln)
             do i = 1, n_spl
                vals(i) = dot_product(radI_splines(i-1:i+1, i_radI), &
                &                     (/1.d0, 4.d0, 1.d0/) / 6.d0)
                err = abs(vals(i))
                if (i > i_iarange) thismaxerr = max(thismaxerr, err)
                if (err > wave_threshold**2) i_this = max(i_this, i)
             end do
             if (out_vals) then
                write(filename, "('radI-',I3.3,'-',I3.3,'-',I2.2,'.dat')") &
                & i_bbfn_1, i_bbfn_2, LL
                call debug_plot_log_data(filename, spl_ln1, spl_lnrange, vals)
             end if
          end do
          maxerr = max(maxerr, thismaxerr)
          this_range = logsbt_ilog2r(dble(i_this), spl_ln1, spl_dln)
          if (this_range < iarange(i_bbfn_1, i_bbfn_2)) then
             if (be_verbose) then
                write(info_str, "(2X,'| ',A,2I4,A,F10.4,A,F10.4,A)") &
                & 'Reducing iarange', i_bbfn_1, i_bbfn_2, &
                & ' from', iarange(i_bbfn_1, i_bbfn_2)*bohr, &
                & ' to', this_range*bohr, ' A.'
                call localorb_info(info_str)
             end if
             iarange(i_bbfn_1, i_bbfn_2) = this_range
             iarange(i_bbfn_2, i_bbfn_1) = this_range
          else if (be_verbose) then
             write(info_str, "(2X,'| ',A,2I4,A,F10.4,A,F10.4,A,ES10.2,A)") &
             & 'Not enlarging iarange', i_bbfn_1, i_bbfn_2, &
             & ' from', iarange(i_bbfn_1, i_bbfn_2)*bohr, &
             & ' to', this_range*bohr, ' A (error:', thismaxerr, ').'
             call localorb_info(info_str)
          end if
          iarange_species(i_species_1, i_species_2) &
          & = max(iarange_species(i_species_1, i_species_2), &
          &       iarange(i_bbfn_1, i_bbfn_2))
          iarange_species(i_species_2, i_species_1) &
          & = iarange_species(i_species_1, i_species_2)
       end do
    end do
    write(info_str, "(2X,'| Maximum value outside interaction range:',ES10.2)")&
    & maxerr
    call localorb_info(info_str)

  end subroutine tighten_iarange
  !******
end module tight_binding_auxmat
!******
