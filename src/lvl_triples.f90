!****h* FHI-aims/lvl_triples
!  NAME
!    lvl_triples
!  SYNOPSIS

module lvl_triples

  !  PURPOSE
  !
  !  USES

  use sbt_overlap, only: SBT_N_TIMES
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
  
  integer, private :: max_prods_L, max_LL
  integer, private :: max_n_prods_fnLsp
  integer, allocatable, private :: Lsp2n_prods_fnLsp(:,:)  ! (L, i_species)
  integer, allocatable, private :: Lsp2prods_fn(:,:,:) ! (i_fnL, L, i_species)
  integer, private :: max_n_prodsfn
  integer, allocatable, private :: n_prodsfn(:)   ! (i_species)
  real*8, allocatable, private :: radius_pd(:,:) ! (i_prodsfn, i_species)
  integer, allocatable, private :: fnLpd_to_Lbs(:,:,:) ! (i_fnL, L, i_species)
  integer, allocatable, private :: fnLpd_to_Lbb(:,:,:) ! (i_fnL, L, i_species)
  integer, allocatable, private :: fnLpd_to_fnLbs(:,:,:) ! (i_fnL, L, i_species)
  integer, allocatable, private :: fnLpd_to_fnLbb(:,:,:) ! (i_fnL, L, i_species)
  integer, allocatable, private :: bbfn_to_fnLbb(:) ! (i_bbfn)
  integer, allocatable, private :: bsfn_to_fnLbs(:) ! (i_bsfn)

  real*8, allocatable, private :: ffk_bs(:,:)  ! (1:N, i_bsfn)
  real*8, allocatable, private :: ffr_bs(:,:)  ! (1:N, i_bsfn)
  real*8, allocatable, private :: vvr_bb(:,:)  ! (1:N, i_bbfn)

  real*8, private :: sbt_times(4, SBT_N_TIMES)
  real*8, private :: time_init(4)
  real*8, private :: time_all(4)
  
contains

  !----------------------------------------------------------------------------
  !****s* lvl_triples/initialize_lvl_triples
  !  NAME
  !    initialize_lvl_triples
  !  SYNOPSIS

  subroutine initialize_lvl_triples(ovlp_type)

    !  PURPOSE
    !
    !    Initialize module data, i.e. calculate the logSBTs of the basis
    !    functions and the basis-basbas product functions for each species.
    !
    !  USES

    use basbas_fn_coulomb, only: adams_moulton_wave_integrator
    use basis
    use dimensions
    use runtime_choices
    use cut_coulomb_operator
    use grids, only: n_grid, r_grid_min, r_grid_inc
    use localorb_io, only: localorb_info, OL_norm, use_unit
    use mpi_tasks, only: aims_stop, check_allocation
    use prodbas
    use sbt_overlap
    use species_data, only: species_name
    use timing_core, only: start_timer, stop_timer
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: ovlp_type

    !  INPUTS
    !    o ovlp_type -- Either OVLP_TYPE_COULOMB or OVLP_TYPE_HSE
    !                   (use runtime_choices:hse_omega_hf).
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: N
    real*8 :: lnr0, lnk0, lnrange
    integer :: i_species, i_bsfn, i_bbfn
    integer :: Lbs, Lbb, Lpd, Mbs, Mbb, Mpd
    integer :: i_fnLbs, i_fnLbb, i_fnLpd
    real*8 :: dummy_mp
    real*8, parameter :: power_bias_sbt = 1.5d0
    real*8, parameter :: power_bias_int = 0.d0
    real*8, parameter :: power_bias_add = 0.5 * (3.d0 - power_bias_int) &
    &                                     - power_bias_sbt
    integer, allocatable :: Lsp2n_prods_fnLsp_uptonow(:)
    integer :: n_pdfn_uptonow
    real*8, allocatable :: ffr_bb(:), rr(:)
    integer :: L, i_fnL, M
    integer :: info
    character*150 :: info_str
    character(*), parameter :: func = 'initialize_lvl_triples'

    sbt_times = 0.d0
    time_init = 0.d0

    call start_timer(time_init)

    N = sbtgrid_N
    lnr0 = sbtgrid_lnr0
    lnk0 = sbtgrid_lnk0
    lnrange = sbtgrid_lnrange

    max_prods_L = max_basis_L + max_basbas_L
    max_LL = max_prods_L + max_basis_L

    
    ! --- basbas related stuff (bb)

    allocate(bbfn_to_fnLbb(n_basbas_fns), stat=info)
    call check_allocation(info, 'bbfn_to_fnLbb', func)

    do i_species = 1, n_species
       do L = 0, max_basbas_L
          do i_fnL = 1, Lsp2n_basbas_fnLsp(L, i_species)
             i_bbfn = Lsp2basbas_fn(i_fnL, L, i_species)
             bbfn_to_fnLbb(i_bbfn) = i_fnL
          end do
       end do
    end do

    allocate(rr(N), ffr_bb(N), vvr_bb(N, n_basbas_fns), stat=info)
    call check_allocation(info, 'vvr_bb', func)
    ! Get rr, the values of r on the log-grid
    rr = 1.d0
    call logsbt_scale(N, 1, rr, lnr0, lnrange, 1.d0)

    if (ovlp_type == OVLP_TYPE_CUT) then
       call initialize_cut_coulomb(N, lnr0, lnk0, lnrange, &
       &                           cutCb_width, cutCb_rcut)
    end if

    do i_bbfn = 1, n_basbas_fns
       L = basbasfn_L(i_bbfn)
       i_species = basbasfn_species(i_bbfn)
       call sbt_import_spline(N, ffr_bb, lnr0, lnrange, L, &
       &                      n_grid(i_species), basbas_wave_spl(:,:, i_bbfn),&
       &                      r_grid_min(i_species), r_grid_inc(i_species))
       select case (ovlp_type)
       case(OVLP_TYPE_COULOMB)
          call logsbt_scale(N, 1, ffr_bb, lnr0, lnrange, 1.d0)    ! u(r) = rP(r)
          call adams_moulton_wave_integrator(L, ffr_bb, N, rr, &
          &                                  vvr_bb(:, i_bbfn), dummy_mp)
       case(OVLP_TYPE_HSE)
          call hse_logsbt_integrator(L, N, ffr_bb, lnr0, lnk0, lnrange, &
          &                          hse_omega_hf)
          vvr_bb(:, i_bbfn) = ffr_bb
       case(OVLP_TYPE_LR)
          call hse_logsbt_integrator(L, N, ffr_bb, lnr0, lnk0, lnrange, &
          &                          hse_omega_hf)
          vvr_bb(:, i_bbfn) = ffr_bb
       case(OVLP_TYPE_CUT)
          call cut_coulomb_logsbt_integrator(N, 1, ffr_bb, (/L/), &
          &                                  lnr0, lnk0, lnrange)
          vvr_bb(:, i_bbfn) = ffr_bb
       case default
          call aims_stop('Overlap type not supported', func)
       end select
    end do
    deallocate(ffr_bb, rr)
    if (ovlp_type == OVLP_TYPE_CUT) call deallocate_cut_coulomb()


    ! --- basis related stuff (bs)

    allocate(bsfn_to_fnLbs(n_basis_fns), stat=info)
    call check_allocation(info, 'bsfn_to_fnLbs', func)
    do i_species = 1, n_species
       do L = 0, max_basis_L
          do i_fnL = 1, Lsp2n_basis_fnLsp(L, i_species)
             i_bsfn = Lsp2basis_fn(i_fnL, L, i_species)
             bsfn_to_fnLbs(i_bsfn) = i_fnL
          end do
       end do
    end do

    allocate(ffr_bs(sbtgrid_N, n_basis_fns), stat=info)
    call check_allocation(info, 'ffr_bs', func)
    allocate(ffk_bs(sbtgrid_N, n_basis_fns), stat=info)
    call check_allocation(info, 'ffk_bs', func)

    do i_bsfn = 1, n_basis_fns
       i_species = basisfn_species(i_bsfn)
       L = basisfn_l(i_bsfn)
       call sbt_import_spline(N, ffr_bs(:, i_bsfn), lnr0, lnrange, L, &
       &                      n_grid(i_species), basis_wave_spl(:,:, i_bsfn), &
       &                      r_grid_min(i_species), r_grid_inc(i_species))
    end do
    ffk_bs = ffr_bs

    ! --- Actual bs transform

    ! f(r)                             ! [assume full relative accuracy]
    call logsbt_scale(N, n_basis_fns, ffk_bs, lnr0, lnrange, 3-power_bias_sbt)
    ! f(r)*k**(3-power_bias)           ! [assume full relative accuracy]
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, power_bias_sbt, &
    &                        n_basis_fns, basisfn_l, ffk_bs)
    ! f(k)*k**power_bias               ! [error: epsilon]
    call logsbt_scale(N, n_basis_fns, ffk_bs, lnk0, lnrange, power_bias_add)
    ! f(k)*k**(power_bias+power_add)   ! [error: epsilon * k**power_add]


    ! --- product related stuff (pd)

    allocate(n_prodsfn(n_species), Lsp2n_prods_fnLsp(0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'n_prodsfn', func)

    ! count pds
    Lsp2n_prods_fnLsp = 0
    do i_species = 1, n_species
       do i_bsfn = 1, n_basis_fns
          if (basisfn_species(i_bsfn) /= i_species) cycle
          Lbs = basisfn_l(i_bsfn)
          do i_bbfn = 1, n_basbas_fns
             if (basbasfn_species(i_bbfn) /= i_species) cycle
             Lbb = basbasfn_l(i_bbfn)
             do Lpd = abs(Lbs - Lbb), Lbs + Lbb, 2
                Lsp2n_prods_fnLsp(Lpd, i_species) = Lsp2n_prods_fnLsp(Lpd, i_species) + 1
             end do
          end do
       end do
       n_prodsfn(i_species) = sum(Lsp2n_prods_fnLsp(:, i_species))
    end do
    max_n_prods_fnLsp = maxval(Lsp2n_prods_fnLsp)
    max_n_prodsfn = maxval(n_prodsfn)

    ! enumerate pds
    allocate(Lsp2prods_fn(max_n_prods_fnLsp, 0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'Lsp2prods_fn', func)
    allocate(radius_pd(max_n_prodsfn, n_species), stat=info)
    call check_allocation(info, 'radius_pd', func)
    allocate(fnLpd_to_Lbs(max_n_prods_fnLsp, 0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'fnLpd_to_Lbs', func)
    allocate(fnLpd_to_fnLbs(max_n_prods_fnLsp, 0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'fnLpd_to_fnLbs', func)
    allocate(fnLpd_to_Lbb(max_n_prods_fnLsp, 0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'fnLpd_to_Lbb', func)
    allocate(fnLpd_to_fnLbb(max_n_prods_fnLsp, 0:max_prods_L, n_species), stat=info)
    call check_allocation(info, 'fnLpd_to_fnLbb', func)

    allocate(Lsp2n_prods_fnLsp_uptonow(0:max_prods_L), stat=info)
    call check_allocation(info, 'Lsp2n_prods_fnLsp_uptonow', func)

    do i_species = 1, n_species
       Lsp2n_prods_fnLsp_uptonow = 0
       n_pdfn_uptonow = 0
       do Lbs = 0, max_basis_L
          do i_fnLbs = 1, Lsp2n_basis_fnLsp(Lbs, i_species)
             i_bsfn = Lsp2basis_fn(i_fnLbs, Lbs, i_species)

             do Lbb = 0, max_basbas_L
                do i_fnLbb = 1, Lsp2n_basbas_fnLsp(Lbb, i_species)
                   i_bbfn = Lsp2basbas_fn(i_fnLbb, Lbb, i_species)

                   do Lpd = abs(Lbs - Lbb), Lbs + Lbb, 2
                      n_pdfn_uptonow = n_pdfn_uptonow + 1
                      Lsp2n_prods_fnLsp_uptonow(Lpd) = Lsp2n_prods_fnLsp_uptonow(Lpd) + 1
                      i_fnLpd = Lsp2n_prods_fnLsp_uptonow(Lpd)

                      Lsp2prods_fn(i_fnLpd, Lpd, i_species) = n_pdfn_uptonow
                      radius_pd(n_pdfn_uptonow,i_species) = outer_radius(i_bsfn)
                      fnLpd_to_Lbs(i_fnLpd, Lpd, i_species) = Lbs
                      fnLpd_to_Lbb(i_fnLpd, Lpd, i_species) = Lbb
                      fnLpd_to_fnLbs(i_fnLpd, Lpd, i_species) = i_fnLbs
                      fnLpd_to_fnLbb(i_fnLpd, Lpd, i_species) = i_fnLbb
                   end do
                end do
             end do
          end do
       end do
    end do
    deallocate(Lsp2n_prods_fnLsp_uptonow)

    ! --- Output memory consumption

    ! In principle, the calculations could be rearranged (globally) in a way
    ! to have the outer loop over species pairs and the next loop over
    ! fn-bunches.  This would allow to reuse the product ffk1 without storing
    ! even all of them for a simple species.  But for now, ...
    write(info_str, "(A)") &
    & '| Storage requirements for basis-basbas products (processor-local):'
    call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
    do i_species = 1, n_species
       write(info_str, &
       & "('|   Species ',A,': n_prodsfn:',I6,'; storage:',F12.2,' MiB.')") &
       & trim(species_name(i_species)), n_prodsfn(i_species), &
       & 8*N*n_prodsfn(i_species) / 2.d0**20
       call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
    end do
    call localorb_info('')

    ! --- Splined kernel construction

    call initialize_fast_kernel(N, max_LL, lnr0, lnk0, lnrange, power_bias_int)

    call stop_timer(time_init)
    time_all = time_init

  end subroutine initialize_lvl_triples
  !******
  !----------------------------------------------------------------------------
  !****s* lvl_triples/cleanup_lvl_triples
  !  NAME
  !    cleanup_lvl_triples
  !  SYNOPSIS

  subroutine cleanup_lvl_triples()

    !  PURPOSE
    !    Clean up module data
    !  USES
    use localorb_io, only: localorb_info, OL_norm
    use sbt_overlap, only: output_sbt_timing
    use synchronize_mpi_basic
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
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'cleanup_lvl_triples'

    deallocate(Lsp2n_prods_fnLsp, Lsp2prods_fn)
    deallocate(n_prodsfn, radius_pd)
    deallocate(fnLpd_to_Lbs, fnLpd_to_Lbb, fnLpd_to_fnLbs, fnLpd_to_fnLbb)
    deallocate(bbfn_to_fnLbb, bsfn_to_fnLbs)
    deallocate(ffk_bs, ffr_bs, vvr_bb)

    call sync_timing(time_init(3))
    call sync_timing(time_all(3))
    call output_timeheader('2X', 'Atomic logSBT for LVL-triples', OL_norm)
    call output_sbt_timing(sbt_times, '2X', OL_norm)
    call output_timer('LVL initialization', time_init(3:4), '2X', OL_norm)
    call output_timer('Whole LVL-triples', time_all(3:4), '2X', OL_norm)
    call localorb_info('')

  end subroutine cleanup_lvl_triples
  !******
  !----------------------------------------------------------------------------
  !****s* lvl_triples/calculate_lvl_triples
  !  NAME
  !    calculate_lvl_triples
  !  SYNOPSIS

  subroutine calculate_lvl_triples(i_species_1, i_species_2, &
  &                                n_Rvec, Rvec, ovlp_3fn, d_ovlp_3fn, calc_deriv)

    !  PURPOSE
    !
    !    Calcluate the triple-overlap integrals of the basis functions of
    !       1: The basis functions at i_atom_1
    !       2: The basis functions at i_atom_2
    !       3: The auxiliary basis functions at both i_atom_1 and i_atom_2
    !
    !  USES

    use basis
    use dimensions, only: n_basbas_fns, n_basis_fns
    use mpi_tasks, only: check_allocation
    use prodbas
    use runtime_choices
    use sbt_overlap, only: sbt_twocenter_triples
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2
    integer, intent(IN) :: n_Rvec
    real*8, intent(IN) :: Rvec(3, n_Rvec)
    real*8, intent(OUT) :: ovlp_3fn(max_n_basis_sp, max_n_basis_sp, &
    &                               max_n_basbas_sp, 2, n_Rvec)
    real*8, intent(OUT) :: d_ovlp_3fn(max_n_basis_sp, max_n_basis_sp, &
    &                               max_n_basbas_sp, 2, 3, n_Rvec)
    logical, intent(IN) :: calc_deriv


    !  INPUTS
    !    o i_species_1, i_species_2 -- Atomic species numbers
    !    o n_Rvec -- Number of difference vectors
    !  OUTPUTS
    !    o ovlp_3fn -- Two-center triple-fn integrals
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    logical :: is_onsite
    real*8 :: use_Rvec(3, n_Rvec)
    real*8, allocatable :: tr_ovlp(:,:,:,:), d_tr_ovlp(:,:,:,:,:)
    integer :: i_side
    integer :: i_species_pd, i_species_ot, i_basis_pd, i_basis_ot
    integer :: i_basis_1, i_basis_2, i_basbas_pd
    integer :: info
    character(*), parameter :: func = 'calculate_lvl_triples'

    allocate(tr_ovlp(max_n_basis_sp, max_n_basbas_sp, max_n_basis_sp, n_Rvec), &
    &        stat=info)
    call check_allocation(info, 'tr_ovlp')
    if(calc_deriv) then
       allocate(d_tr_ovlp(max_n_basis_sp, max_n_basbas_sp, max_n_basis_sp, 3, n_Rvec), stat=info)
    else
       ! Dummy allocation
       allocate(d_tr_ovlp(1, 1, 1, 1, 1), stat=info)
    endif
    ovlp_3fn = 0.d0
    tr_ovlp = 0.d0
    if(calc_deriv) d_tr_ovlp = 0.d0

    do i_side = 1, 2
       if (i_side == 1) then
          i_species_pd = i_species_1
          i_species_ot = i_species_2
          use_Rvec = Rvec
       else
          i_species_pd = i_species_2
          i_species_ot = i_species_1
          use_Rvec = - Rvec
       end if

       call sbt_twocenter_triples( &
       & sbtgrid_N, sbtgrid_lnr0, sbtgrid_lnk0, sbtgrid_lnrange, &
       & n_Rvec, use_Rvec, &
       & n_prodsfn(i_species_pd), radius_pd(:, i_species_pd), &
       &   max_prods_L, max_n_prods_fnLsp, Lsp2n_prods_fnLsp(:, i_species_pd), &
       &                         Lsp2prods_fn(:,:, i_species_pd), &
       & n_basis_fns, ffk_bs, outer_radius, &
       &   max_basis_L, max_n_basis_fnLsp, Lsp2n_basis_fnLsp(:, i_species_ot), &
       &                         Lsp2basis_fn(:,:, i_species_ot), &
       & max_basis_L, max_n_basis_fnLsp, Lsp2n_basis_fnLsp(:, i_species_pd), &
       & max_basbas_L, max_n_basbas_fnLsp, Lsp2n_basbas_fnLsp(:, i_species_pd),&
       & fnLpd_to_Lbs(:,:, i_species_pd), fnLpd_to_fnLbs(:,:, i_species_pd), &
       & fnLpd_to_Lbb(:,:, i_species_pd), fnLpd_to_fnLbb(:,:, i_species_pd), &
       & n_basis_fns, ffr_bs, Lsp2basis_fn(:,:, i_species_pd), &
       & n_basbas_fns, vvr_bb, Lsp2basbas_fn(:,:, i_species_pd), &
       & Lsp2basis_sp(:,:, i_species_pd), &
       & Lsp2basbas_sp(:,:, i_species_pd), &
       & Lsp2basis_sp(:,:, i_species_ot), &
       & tr_ovlp, d_tr_ovlp, calc_deriv, sbt_times)

       do i_basis_pd = 1, sp2n_basis_sp(i_species_pd)
          do i_basis_ot = 1, sp2n_basis_sp(i_species_ot)
             if (i_side == 1) then
                i_basis_1 = i_basis_pd
                i_basis_2 = i_basis_ot
             else
                i_basis_1 = i_basis_ot
                i_basis_2 = i_basis_pd
             end if
             do i_basbas_pd = 1, sp2n_basbas_sp(i_species_pd)
                ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :) &
                & = tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :)
                if(calc_deriv) then
                   if (i_side == 1) then
                      d_ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :, :) &
                      & =  d_tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :, :)
                   else
                      d_ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :, :) &
                      & = -d_tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :, :)
                   endif
                endif
             end do
          end do
       end do
    end do

    deallocate(tr_ovlp, d_tr_ovlp)

  end subroutine calculate_lvl_triples
  !******
  !----------------------------------------------------------------------------
  !****s* lvl_triples/calculate_lvl_triples_vb
  !  NAME
  !    calculate_lvl_triples_vb
  !  SYNOPSIS

  subroutine calculate_lvl_triples_vb(i_atom_1, i_atom_2, i_species_1, i_species_2, atom2basis_len2, &
                                atom2vb_basis_off, n_atoms2, n_Rvec, Rvec, ovlp_3fn, d_ovlp_3fn, calc_deriv)

    !  PURPOSE
    !
    !    Calcluate the triple-overlap integrals of the basis functions of
    !       1: The basis functions at i_atom_1
    !       2: The basis functions at i_atom_2
    !       3: The auxiliary basis functions at both i_atom_1 and i_atom_2
    !
    !  USES

    use basis
    use dimensions, only: n_basbas_fns, n_basis_fns,n_atoms
    use mpi_tasks, only: check_allocation
    use prodbas
    use runtime_choices
    use sbt_overlap, only: sbt_twocenter_triples
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_species_1, i_species_2,i_atom_1, i_atom_2,n_atoms2
    integer, intent(IN) :: n_Rvec, atom2basis_len2(n_atoms2),atom2vb_basis_off(n_atoms2)
    real*8, intent(IN) :: Rvec(3, n_Rvec)
    real*8, intent(OUT) :: ovlp_3fn(max_n_basis_sp2, max_n_basis_sp2, &
    &                               max_n_basbas_sp, 2, n_Rvec)
    real*8, intent(OUT) :: d_ovlp_3fn(max_n_basis_sp2, max_n_basis_sp2, &
    &                               max_n_basbas_sp, 2, 3, n_Rvec)
    logical, intent(IN) :: calc_deriv


    !  INPUTS
    !    o i_species_1, i_species_2 -- Atomic species numbers
    !    o n_Rvec -- Number of difference vectors
    !  OUTPUTS
    !    o ovlp_3fn -- Two-center triple-fn integrals
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    logical :: is_onsite
    real*8 :: use_Rvec(3, n_Rvec)
    real*8, allocatable :: tr_ovlp(:,:,:,:), d_tr_ovlp(:,:,:,:,:)
    integer :: i_side
    integer :: i_species_pd, i_species_ot, i_basis_pd, i_basis_ot,i_atom_pd, i_atom_ot
    integer :: i_basis_1, i_basis_2, i_basbas_pd
    integer :: info
    character(*), parameter :: func = 'calculate_lvl_triples'

    allocate(tr_ovlp(max_n_basis_sp, max_n_basbas_sp, max_n_basis_sp, n_Rvec), &
    &        stat=info)
    call check_allocation(info, 'tr_ovlp')
    if(calc_deriv) then
       allocate(d_tr_ovlp(max_n_basis_sp, max_n_basbas_sp, max_n_basis_sp, 3, n_Rvec), stat=info)
    else
       ! Dummy allocation
       allocate(d_tr_ovlp(1, 1, 1, 1, 1), stat=info)
    endif
    ovlp_3fn = 0.d0
    tr_ovlp = 0.d0
    if(calc_deriv) d_tr_ovlp = 0.d0

    do i_side = 1, 2
       if (i_side == 1) then
          i_species_pd = i_species_1
          i_species_ot = i_species_2
          i_atom_pd = i_atom_1
          i_atom_ot = i_atom_2
          use_Rvec = Rvec
       else
          i_species_pd = i_species_2
          i_species_ot = i_species_1
          i_atom_pd = i_atom_2
          i_atom_ot = i_atom_1
          use_Rvec = - Rvec
       end if

       call sbt_twocenter_triples( &
       & sbtgrid_N, sbtgrid_lnr0, sbtgrid_lnk0, sbtgrid_lnrange, &
       & n_Rvec, use_Rvec, &
       & n_prodsfn(i_species_pd), radius_pd(:, i_species_pd), &
       &   max_prods_L, max_n_prods_fnLsp, Lsp2n_prods_fnLsp(:, i_species_pd), &
       &                         Lsp2prods_fn(:,:, i_species_pd), &
       & n_basis_fns, ffk_bs, outer_radius, &
       &   max_basis_L, max_n_basis_fnLsp, Lsp2n_basis_fnLsp(:, i_species_ot), &
       &                         Lsp2basis_fn(:,:, i_species_ot), &
       & max_basis_L, max_n_basis_fnLsp, Lsp2n_basis_fnLsp(:, i_species_pd), &
       & max_basbas_L, max_n_basbas_fnLsp, Lsp2n_basbas_fnLsp(:, i_species_pd),&
       & fnLpd_to_Lbs(:,:, i_species_pd), fnLpd_to_fnLbs(:,:, i_species_pd), &
       & fnLpd_to_Lbb(:,:, i_species_pd), fnLpd_to_fnLbb(:,:, i_species_pd), &
       & n_basis_fns, ffr_bs, Lsp2basis_fn(:,:, i_species_pd), &
       & n_basbas_fns, vvr_bb, Lsp2basbas_fn(:,:, i_species_pd), &
       & Lsp2basis_sp(:,:, i_species_pd), &
       & Lsp2basbas_sp(:,:, i_species_pd), &
       & Lsp2basis_sp(:,:, i_species_ot), &
       & tr_ovlp, d_tr_ovlp, calc_deriv, sbt_times)

       do i_basis_pd = atom2vb_basis_off(i_atom_pd)+1, atom2vb_basis_off(i_atom_pd)+atom2basis_len2(i_atom_pd)
          do i_basis_ot = atom2vb_basis_off(i_atom_ot)+1,atom2vb_basis_off(i_atom_ot)+atom2basis_len2(i_atom_ot)
             if (i_side == 1) then
                i_basis_1 = i_basis_pd-atom2vb_basis_off(i_atom_pd)
                i_basis_2 = i_basis_ot-atom2vb_basis_off(i_atom_ot)
             else
                i_basis_1 = i_basis_ot-atom2vb_basis_off(i_atom_ot)
                i_basis_2 = i_basis_pd-atom2vb_basis_off(i_atom_pd)
             end if
             do i_basbas_pd = 1, sp2n_basbas_sp(i_species_pd)
                ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :) &
                & = tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :)
                if(calc_deriv) then
                   if (i_side == 1) then
                      d_ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :, :) &
                      & =  d_tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :, :)
                   else
                      d_ovlp_3fn(i_basis_1, i_basis_2, i_basbas_pd, i_side, :, :) &
                      & = -d_tr_ovlp(i_basis_pd, i_basbas_pd, i_basis_ot, :, :)
                   endif
                endif
             end do
          end do
       end do
    end do

    deallocate(tr_ovlp, d_tr_ovlp)

  end subroutine calculate_lvl_triples_vb
  !******
  !----------------------------------------------------------------------------
  !****s* lvl_triples/calculate_lvl_triples_matrix
  !  NAME
  !    calculate_lvl_triples_matrix
  !  SYNOPSIS

  subroutine calculate_lvl_triples_matrix(i_atom_1, i_atom_2, n_row, n_col, &
  &                                       row2basbas, col2basis_pair, ovlp_3fn)

    !  PURPOSE
    !
    !    Calcluate the triple-overlap integrals of the basis functions of
    !       1: The basis functions at i_atom_1
    !       2: The basis functions at i_atom_2
    !       3: The auxiliary basis functions at both i_atom_1 and i_atom_2
    !
    !  USES

    use basis
    use geometry, only: coords, species
    use mpi_tasks, only: aims_stop, check_allocation
    use prodbas
    use timing_core, only: start_timer, stop_timer
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: i_atom_1, i_atom_2
    integer, intent(IN) :: n_row, n_col
    integer, intent(IN) :: row2basbas(n_row)
    integer, intent(IN) :: col2basis_pair(2, n_col)
    real*8, intent(INOUT) :: ovlp_3fn(:,:)

    !  INPUTS
    !    o i_atom_1, i_atom_2 -- Atom numbers
    !    o n_row -- Number of auxiliary basis functions in question
    !    o n_col -- Number of basis pairs in question
    !    o row2basbas -- List of auxiliary basis function indices
    !    o col2basis_pair -- List of basis function indices
    !    o ovlp_3fn -- Input ovlp_3fn
    !                  Things indexed by {row,col}2* get overwritten.
    !  OUTPUTS
    !    o ovlp_3fn -- Two-center triple-fn integrals
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    real*8 :: Rvec(3), dummy(1, 1, 1, 1, 1, 1)
    real*8, allocatable :: ovlp(:,:,:,:,:)
    integer :: i_species_1, i_species_2
    integer :: i_side
    integer :: i_atom_pd, i_species_pd
    integer :: i_row, i_col
    integer :: i_basbas, i_bbfn, Lbb, Mbb, i_fnLbb
    integer :: i_basis_1, i_bsfn_1, Lbs_1, Mbs_1, i_fnLbs_1
    integer :: i_basis_2, i_bsfn_2, Lbs_2, Mbs_2, i_fnLbs_2
    integer :: i1, i2, i3
    integer :: info
    character(*), parameter :: func = 'calculate_lvl_triples_matrix'

    call start_timer(time_all)

    i_species_1 = species(i_atom_1)
    i_species_2 = species(i_atom_2)
    Rvec = coords(:, i_atom_2) - coords(:, i_atom_1)

    allocate(ovlp(max_n_basis_sp, max_n_basis_sp, max_n_basbas_sp, 2, 1), &
    &        stat=info)
    call check_allocation(info, 'ovlp', func)

    call calculate_lvl_triples(i_species_1, i_species_2, 1, Rvec, ovlp, dummy, .false.)

    do i_row = 1, n_row
       i_basbas = row2basbas(i_row)
       if (i_basbas <= 0) cycle
       i_atom_pd = basbas_atom(i_basbas)
       i_species_pd = species(i_atom_pd)
       if (i_atom_pd == i_atom_1) then
          i_side = 1
       else if (i_atom_pd == i_atom_2) then
          i_side = 2
       else
          call aims_stop('invalid i_atom_pd',func)
          !stop   ! Satisfy compiler
       end if

       i_bbfn = basbas_fn(i_basbas)
       Lbb = basbas_l(i_basbas)
       Mbb = basbas_m(i_basbas)
       i_fnLbb = bbfn_to_fnLbb(i_bbfn)

       i3 = Lsp2basbas_sp(i_fnLbb, Lbb, i_species_pd) + Mbb

       do i_col = 1, n_col
          i_basis_1 = col2basis_pair(1, i_col)
          i_bsfn_1 = basis_fn(i_basis_1)
          Lbs_1 = basis_l(i_basis_1)
          Mbs_1 = basis_m(i_basis_1)
          i_fnLbs_1 = bsfn_to_fnLbs(i_bsfn_1)
          i1 = Lsp2basis_sp(i_fnLbs_1, Lbs_1, i_species_1) + Mbs_1

          i_basis_2 = col2basis_pair(2, i_col)
          i_bsfn_2 = basis_fn(i_basis_2)
          Lbs_2 = basis_l(i_basis_2)
          Mbs_2 = basis_m(i_basis_2)
          i_fnLbs_2 = bsfn_to_fnLbs(i_bsfn_2)
          i2 = Lsp2basis_sp(i_fnLbs_2, Lbs_2, i_species_2) + Mbs_2

          ovlp_3fn(i_row, i_col) = ovlp(i1, i2, i3, i_side, 1)
       end do
    end do

    call stop_timer(time_all, .true.)

  end subroutine calculate_lvl_triples_matrix
  !******
  !----------------------------------------------------------------------------
end module lvl_triples
!******
