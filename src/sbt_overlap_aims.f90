!****h* FHI-aims/sbt_overlap_aims
!  NAME
!    sbt_overlap_aims
!  SYNOPSIS

module sbt_overlap_aims

  !  PURPOSE
  !    Provide infrastructure to calculate overlap integrals using
  !    spherical Bessel transforms.  In contrast to sbt_overlap.f90,
  !    this module is aware of the global data structures outside the
  !    logSBT realm.
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
  !****s* sbt_overlap/integrate_auxmat_by_atomic_sbt
  !  NAME
  !   integrate_auxmat_by_atomic_sbt
  !  SYNOPSIS

  subroutine integrate_auxmat_by_atomic_sbt(auxmat, ovlp_type, distr_2d, &
  &                                         current_cell_index, use_tb,  &
  &                                         opt_basbas2row, opt_basbas2col)

    !  PURPOSE
    !
    !    This procedure is intended to calculate the coulomb interaction
    !    matrix elements between every two auxiliary basis functions.
    !
    !  USES

    use dimensions
    use runtime_choices
    use prodbas
    use geometry, only: coords, lattice_vector, species
    use localorb_io, only: localorb_info, OL_norm
    use sbt_overlap_tb, only: sbt_atomic_ovlp_tb
    use timing_core, only: output_timeheader
    use mpi_tasks, only: aims_stop, check_allocation
    use sbt_overlap
    implicit none

    !  ARGUMENTS

    real*8, intent(OUT) :: auxmat(:,:)
    integer, intent(IN) :: ovlp_type
    logical, intent(IN) :: distr_2d
    integer, intent(IN), optional :: current_cell_index(3)
    logical, intent(IN), optional :: use_tb
    integer, intent(IN), optional :: opt_basbas2row(n_basbas)
    integer, intent(IN), optional :: opt_basbas2col(n_basbas)

    !  INPUTS
    !     o ovlp_type -- Either OVLP_TYPE_OVERLAP, OVLP_TYPE_COULOMB,
    !                    or OVLP_TYPE_HSE (use runtime_choices:hse_omega_hf).
    !     o distr_2d -- Use 2d scalapack distribution?
    !     o current_cell_index -- Unit cell of the column product basis funcs
    !     o use_tb -- Use proof-of-principle tight-binding code (slow!)
    !     o opt_basbas2row -- user defined row distribution of output matrix
    !     o opt_basbas2col -- user defined col distribution of output matrix
    !     o ... global data in prodbas.f90 and friends
    !  OUTPUTS
    !     o auxmat -- real array, the Coulomb interaction/overlap matrix within
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

    real*8 :: times(4, SBT_N_TIMES)
    logical :: has_mp_far
    integer :: max_L, max_LL, N
    real*8 :: lnr0, lnk0, lnrange
    real*8, allocatable :: ffk(:,:), radii(:)
    real*8 :: power_bias_int
    logical :: my_use_tb

    integer :: i_basbas, i_basbas_fn, M, L, i_fnL, i
    integer :: i_atom_1, i_atom_2, i_species_1, i_species_2
    real*8 :: Rvec(3), k, hse_fac, exp_arg
    integer, allocatable :: basbas2row(:), basbas2col(:)
    integer, allocatable :: fnL1_to_row(:,:,:), fnL2_to_col(:,:,:)
    real*8, allocatable :: moments(:)
    character*150 :: info_str
    integer :: info
    character(*), parameter :: func = 'integrate_auxmat_by_atomic_sbt'

    my_use_tb = .false.; if (present(use_tb)) my_use_tb = use_tb

    auxmat = 0.d0
    times = 0.d0

    write(info_str, "(2X,A,A,A)") &
    & "Integrating the ", trim(OVLP_TYPE_NAMES(ovlp_type)), &
    & " interaction matrix for auxiliary basis functions (by atoms)"
    call localorb_info(info_str)
    call localorb_info('')

    ! --- Set procedure parameters


    max_L = max_basbas_L
    max_LL = 2*max_L   ! Worst case: max_L x max_L matrix element
    N = sbtgrid_N
    lnr0 = sbtgrid_lnr0; lnk0 = sbtgrid_lnk0; lnrange = sbtgrid_lnrange

    allocate(moments(n_basbas_fns), radii(n_basbas_fns), stat=info)
    call check_allocation(info, 'moments, radii', func)

    allocate(ffk(N, n_basbas_fns), stat=info)
    call check_allocation(info, 'ffk', func)

    call sbt_atomic_field_transforms(N, lnr0, lnk0, lnrange, ovlp_type, &
    &          n_basbas_fns, basbasfn_l, basbasfn_species, basbas_wave_spl, &
    &          ffk, power_bias_int)


    select case (ovlp_type)
    case(OVLP_TYPE_OVERLAP)
       radii = charge_radius_basbas_fn
       moments = 0.d0
       has_mp_far = .false.
    case(OVLP_TYPE_COULOMB)
       radii = charge_radius_basbas_fn
       moments = multipole_basbas_fn
       has_mp_far = .true.
    case(OVLP_TYPE_HSE)
       radii = field_radius_basbas_fn   ! This is too pessimistic.
       moments = 0.d0
       has_mp_far = .false.
    case(OVLP_TYPE_LR)
       radii = field_radius_basbas_fn   
       moments = 0.d0
       has_mp_far = .false.
    case(OVLP_TYPE_CUT)
       ! field_radius is where even the outermoste charge is cut
       ! (widthfac**5 means an argument of 5 in erfc ~ exp(-6**2) ~ 2e-16).
       ! But is too pessimistic anyway as it only drops if the fields stop
       ! overlapping.
       radii = charge_radius_basbas_fn + cutCb_rcut * cutCb_width**6
       moments = 0.d0
       has_mp_far = .false.
    case(OVLP_TYPE_CUT_ANALYTIC)
       ! IGOR for analytic truncated Coulomb
       radii = charge_radius_basbas_fn + cutCb_rcut * cutCb_width**6
       moments = 0.d0
       has_mp_far = .false.
    case(OVLP_TYPE_ERS) ! explite range separation, IGOR
       radii = field_radius_basbas_fn   ! This is too pessimistic.
       moments = 0.d0
       has_mp_far = .false.
    case default
       call aims_stop('Invalid ovlp_type', func)
    end select


    ! --- Prepare distribution index arrays

    allocate(basbas2row(n_basbas), basbas2col(n_basbas), stat=info)
    call check_allocation(info, 'basbas2xxx', func)
    if(present(opt_basbas2row) .and. present(opt_basbas2col)) then
      ! distribution is set by caller
      basbas2row(:) = opt_basbas2row(:)
      basbas2col(:) = opt_basbas2col(:)
    elseif(present(opt_basbas2row) .or. present(opt_basbas2col)) then
      ! Safety only
      call aims_stop('Either both opt_basbas2row and opt_basbas2col or none must be present')
    else
      ! set distribution according to distr_2d
      call get_basbas_to_rowcol(distr_2d, basbas2row, basbas2col)
    endif

    allocate(fnL1_to_row(-max_L:max_L, max_n_basbas_fnLsp, 0:max_L),stat=info)
    call check_allocation(info, 'basbasfn2row', func)
    allocate(fnL2_to_col(-max_L:max_L, max_n_basbas_fnLsp, 0:max_L),stat=info)
    call check_allocation(info, 'basbasfn2col', func)

    ! --- Prepare splined kernel

    call initialize_fast_kernel(N, max_LL, lnr0, lnk0, lnrange, power_bias_int)

    ! --- Main loop over atom pairs

    do i_atom_2 = 1, n_atoms
! DEBUGG: 7/Dez/12 .. wieder nur bis n_real_atoms
!    do i_atom_2 = 1, n_real_atoms
       i_species_2 = species(i_atom_2)

       call get_fnL_to_rowcol(i_atom_2, &
       & n_basbas, basbas_atom, basbas_fn, basbas_l, basbas_m, &
       & max_L, max_n_basbas_fnLsp, Lsp2n_basbas_fnLsp(:, i_species_2), &
       & Lsp2basbas_fn(:,:, i_species_2), basbas2col, fnL2_to_col)

       do i_atom_1 = 1, n_atoms    ! For now, ignore symmetry of Coulomb matrix
! DEBUGG: 7/Dez/12 .. wieder nur bis n_real_atoms
!       do i_atom_1 = 1, n_real_atoms    ! For now, ignore symmetry of Coulomb matrix
          i_species_1 = species(i_atom_1)

          call get_fnL_to_rowcol(i_atom_1, &
          & n_basbas, basbas_atom, basbas_fn, basbas_l, basbas_m, &
          & max_L, max_n_basbas_fnLsp, Lsp2n_basbas_fnLsp(:, i_species_1), &
          & Lsp2basbas_fn(:,:, i_species_1), basbas2row, fnL1_to_row)

          Rvec = coords(:, i_atom_2) - coords(:, i_atom_1)
          if (present(current_cell_index)) then
             Rvec = Rvec + matmul(lattice_vector, current_cell_index)
          end if

          if (my_use_tb) then
             call sbt_atomic_ovlp_tb(N, lnr0, lnk0, lnrange, Rvec, &
             & 1.d0, has_mp_far, &
             & max_L, max_L, max_n_basbas_fnLsp, max_n_basbas_fnLsp, &
             & n_basbas_fns, ffk, radii, moments, &
             &   Lsp2n_basbas_fnLsp(:, i_species_1), &
             &   Lsp2basbas_fn(:,:, i_species_1), fnL1_to_row, &
             & n_basbas_fns, ffk, radii, moments, &
             &   Lsp2n_basbas_fnLsp(:, i_species_2), &
             &   Lsp2basbas_fn(:,:, i_species_2), fnL2_to_col, &
             & auxmat, times)
          else
             call sbt_atomic_ovlp(N, lnr0, lnk0, lnrange, Rvec, &
             & 1.d0, has_mp_far, &
             & max_L, max_L, max_n_basbas_fnLsp, max_n_basbas_fnLsp, &
             & n_basbas_fns, ffk, radii, moments, &
             &   Lsp2n_basbas_fnLsp(:, i_species_1), &
             &   Lsp2basbas_fn(:,:, i_species_1), fnL1_to_row, &
             & n_basbas_fns, ffk, radii, moments, &
             &   Lsp2n_basbas_fnLsp(:, i_species_2), &
             &   Lsp2basbas_fn(:,:, i_species_2), fnL2_to_col, &
             & auxmat, times)
          end if
       end do
    end do

    ! --- Deallocate

    call cleanup_fast_kernel()
    deallocate(basbas2row, basbas2col, fnL1_to_row, fnL2_to_col)
    deallocate(ffk, moments)

    if(.not.present(current_cell_index))then
       write(info_str, "('Atomic logSBT for 2-center ',A,' matrix')") &
            & OVLP_TYPE_NAMES(ovlp_type)
       call output_timeheader('2X', info_str, OL_norm)
       call output_sbt_timing(times, '2X', OL_norm)
       call localorb_info('')
    endif


  end subroutine integrate_auxmat_by_atomic_sbt
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap/sbt_atomic_transforms
  !  NAME
  !    sbt_atomic_transforms
  !  SYNOPSIS

  subroutine sbt_atomic_transforms(N, lnr0, lnk0, lnrange, &
  &                 power_bias, power_add, &
  &                 n_bas_fns, basfn_l, basfn_species, basfn_wave_spl, &
  &                 fn_sbt)

    !  PURPOSE
    !    Evaluate the basfn_wave splines at the SBT grid points
    !    If needed, extrapolation to the left is done by assuming r^L_eff
    !    scaling, while to the right it is assumed to vanish.
    !  USES

    use dimensions
    use runtime_choices
    use debug_output, only: debug_plot_log_data
    use grids, only: n_grid, r_grid_min, r_grid_inc
    use localorb_io, only: localorb_info, OL_norm, use_unit
    use mpi_tasks, only: myid
    use species_data, only: species_name
    use mpi_tasks, only: check_allocation
    use sbt_overlap
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange, power_bias, power_add
    integer, intent(IN) :: n_bas_fns
    integer, intent(IN) :: basfn_l(n_bas_fns), basfn_species(n_bas_fns)
    real*8, intent(IN) :: basfn_wave_spl(n_max_spline, n_max_grid, n_bas_fns)
    real*8, intent(OUT) :: fn_sbt(N, n_bas_fns)

    !  INPUTS
    !    o N -- Number of SBT grid points
    !    o lnr0, lnk0, lnrange -- SBT grid extents
    !    o power_bias -- SBT power bias (alpha)
    !    o power_add -- fn_sbt will contain SBT * k**(power_bias+power_add)
    !                   power_add should not be too large (<= 0.5 or so).
    !    o n_bas_fns -- Number of radial parts
    !    o basfn_l -- Angular momenta
    !    o basfn_species -- Species of radial spline (for grid)
    !    o basfn_wave_spl -- Actual radial part (splined)
    !    ! o n_grid(i_species) -- Number of (spline) log-grid points
    !    ! o r_grid_min(i_species) -- First (spline) log-grid point
    !    ! o r_grid_inc(i_species) -- Increase factor of (spline) log-grid pts
    !  OUTPUTS
    !    o fn_sbt -- Interpolated quantity at SBT grid points
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer :: max_L, L, i_species, i_fn
    integer, allocatable :: n_in_channel(:,:)
    real*8, allocatable :: fn_tmp(:)
    real*8 :: sbt_errors(n_sbt_errors, n_bas_fns)
    integer :: info
    character*150 :: info_str, filename
    character(*), parameter :: func = 'sbt_atomic_transforms'

    max_L = maxval(basfn_l)

    ! --- Do SBT

    write(info_str,"(2X,'SBT integration errors (all should be ''small''):')")
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str,"(2X,'large logFT aliasing -> increase N/lnrange')")
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str,"(2X,A,' -> ',A)") &
    & 'large SBT aliasing', 'decrease lnk0 & increase lnk0+lnrange'
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str,"(2X,A,' -> ',A)") &
    & 'large SBT ringing', 'decrease lnr0 & increase lnr0+lnrange'
    call localorb_info(info_str, use_unit, '(A)', OL_norm)
    write(info_str, "(2X,'|           ',A5,A3,':',5A10)") &
    & 'El''t', 'L', 'logFT-al.', 'small-k', 'large-k', 'small-r', 'large-r'
    call localorb_info(info_str, use_unit, '(A)', OL_norm)

    if (out_basis) then
       allocate(fn_tmp(N), stat=info)
       call check_allocation(info, 'fn_tmp', func)
       allocate(n_in_channel(0:max_L, n_species), stat=info)
       call check_allocation(info, 'n_in_channel', func)
       n_in_channel = 0
    end if

    do i_fn = 1, n_bas_fns

       L = basfn_l(i_fn)
       i_species = basfn_species(i_fn)

       ! spline -> sbtgrid
       call sbt_import_spline(N, fn_sbt(:, i_fn), lnr0, lnrange, L, &
       &                      n_grid(i_species), basfn_wave_spl(:,:, i_fn), &
       &                      r_grid_min(i_species), r_grid_inc(i_species))

       ! Output radial function on sbtgrid
       if (out_basis) then
          if (myid.eq.0) then

          n_in_channel(L, i_species) = n_in_channel(L, i_species) + 1
          write(filename, "('prodbas_sbtgrid-',A,'-',I1,'-',I2.2,'.dat')") &
          & trim(species_name(i_species)), L, n_in_channel(L, i_species)
          if (use_hse .and. (hse_omega_hf /= 0.d0).and. .not. use_gw_and_hse &
              .or. lrc_pt2_started) then
             fn_tmp = fn_sbt(:, i_fn)
             call hse_logsbt_integrator(L, N, fn_tmp, &
             &                          lnr0, lnk0, lnrange, hse_omega_hf)
             call debug_plot_log_data(filename, lnr0, lnrange, &
             &                        fn_sbt(:, i_fn), fn_tmp)
          else
             call debug_plot_log_data(filename, lnr0, lnrange, fn_sbt(:, i_fn))
          end if

          end if
       end if

    end do


    ! ff == f(r)                     ! [assume full relative accuracy]
    call logsbt_scale(N, n_bas_fns, fn_sbt, lnr0, lnrange, 3.d0-power_bias)
    ! ff == f(r)*r**(3-power_bias)   ! [assume full relative accuracy]
    call logsbt_multi_driver(N, lnr0, lnk0, lnrange, power_bias, &
    &                        n_bas_fns, basfn_l, fn_sbt, sbt_errors)
    ! ff == f(k)*k**power_bias       ! [error: epsilon]
    call logsbt_scale(N, n_bas_fns, fn_sbt, lnk0, lnrange, power_add)
    ! ff == f(k)*k**(p_bias+p_add)   ! [error: epsilon * k**power_add]

    if (out_basis) n_in_channel = 0

    do i_fn = 1, n_bas_fns

       ! Output radial SBT on sbtgrid
       if (out_basis) then
          if (myid.eq.0) then

          fn_tmp = fn_sbt(:, i_fn)
          i_species = basfn_species(i_fn)
          L = basfn_l(i_fn)
          n_in_channel(L, i_species) = n_in_channel(L, i_species) + 1
          call logsbt_scale(N, 1, fn_tmp, lnk0, lnrange,-(power_bias+power_add))
          write(filename, "('prodbas_logsbt-',A,'-',I1,'-',I2.2,'.dat')") &
          & trim(species_name(i_species)), L, n_in_channel(L, i_species)
          call debug_plot_log_data(filename, lnk0, lnrange, fn_sbt(:, i_fn))

          end if
       end if

       ! check errors
!       write(info_str, "(2X,'| Errors for',A5,I3,':',5ES10.2)") &
!       & trim(species_name(basfn_species(i_fn))), basfn_l(i_fn), &
!       & sbt_errors(i_logft_aliasing, i_fn), &
!       & sbt_errors(i_sbt_aliasing_smallk, i_fn), &
!       & sbt_errors(i_sbt_aliasing_largek, i_fn), &
!       & sbt_errors(i_sbt_ringing_smallr, i_fn), &
!       & sbt_errors(i_sbt_ringing_larger, i_fn)
!       call localorb_info(info_str, use_unit, '(A)', OL_norm)

    end do
!    call localorb_info('', use_unit, '(A)', OL_norm)

    if (out_basis) then
       deallocate(n_in_channel)
       deallocate(fn_tmp)
    end if


  end subroutine sbt_atomic_transforms
  !******
  !----------------------------------------------------------------------------
  !****s* sbt_overlap_aims/sbt_atomic_field_transforms
  !  NAME
  !    sbt_atomic_field_transforms
  !  SYNOPSIS

  subroutine sbt_atomic_field_transforms(N, lnr0, lnk0, lnrange, ovlp_type, &
  &                 n_bas_fns, basfn_l, basfn_species, basfn_wave_spl, &
  &                 ffk, power_bias_int)

    !  PURPOSE
    !
    !    Calculate the SBT of the product basis functions (or their fields).
    !
    !  USES

    use dimensions
    use runtime_choices
    use prodbas
    use constants, only: pi
    use cut_coulomb_operator
    use mpi_tasks, only: aims_stop
    implicit none

    !  ARGUMENTS

    integer, intent(IN) :: N
    real*8, intent(IN) :: lnr0, lnk0, lnrange
    integer, intent(IN) :: ovlp_type
    integer, intent(IN) :: n_bas_fns
    integer, intent(IN) :: basfn_l(n_bas_fns), basfn_species(n_bas_fns)
    real*8, intent(IN) :: basfn_wave_spl(n_max_spline, n_max_grid, n_bas_fns)
    real*8, intent(OUT) :: ffk(N, n_bas_fns)
    real*8, intent(OUT) :: power_bias_int

    !  INPUTS
    !    o N -- Number of SBT grid points
    !    o lnr0, lnk0, lnrange -- SBT grid extents
    !    o ovlp_type -- One of OVLP_TYPE_{OVERLAP,COULOMB,HSE}
    !    o n_bas_fns -- Number of radial parts
    !    o basfn_l -- Angular momenta
    !    o basfn_species -- Species of radial spline (for grid)
    !    o basfn_wave_spl -- Actual radial part (splined)
    !  OUTPUTS
    !    o ffk -- SBT of aux field
    !    o power_bias_int -- 0.d0 (k**power_bias_int is still missing in
    !                              ffk(:,i)*ffk(:,j) to k**3)
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    real*8 :: add_k_power, power_bias_sbt, power_bias_add
    integer :: i
    real*8 :: k, hse_fac
    real*8, external :: radial_fourier_hse
    character(*), parameter :: func = 'sbt_atomic_field_transforms'

    select case (ovlp_type)
    case(OVLP_TYPE_OVERLAP)
       add_k_power = 0.d0
       ! As the radial parts get scaled by k^1.5 each,
       ! ffk(k) * k^1.5 should be most accurate.
       power_bias_sbt = 1.5d0
    case(OVLP_TYPE_COULOMB)
       add_k_power = -2.d0
       ! In principle, power_bias_sbt=0.5, but j_0 would probably make trouble.
       ! 0. is near enough and works because of the hybrid kernel approach.
       power_bias_sbt = 0.d0
    case(OVLP_TYPE_HSE)
       add_k_power = 0.d0
       ! HSE has no real (multipole-like) far field.  None of the "tricks"
       ! for the Coulomb potential (different scaling, far-field) are needed.
       power_bias_sbt = 1.5d0
       if (use_erfc) hse_omega_hf = erfc_omega
    case(OVLP_TYPE_LR)
       add_k_power = 0.d0
       power_bias_sbt = 1.5d0
    case(OVLP_TYPE_CUT)
       add_k_power = 0.d0
       ! Field will be scaled by 1/k^2 cut to the left to a finite value.
       ! Therefore, large k are less of a problem and a slightly reduced
       ! power_bias should be beneficial.
       power_bias_sbt = 1.d0
    case(OVLP_TYPE_CUT_ANALYTIC)
       add_k_power = 0.d0
       ! IGOR for anaytic cut coulomb
       ! Field will be scaled by 1/k^2 cut to the left to a finite value.
       ! Therefore, large k are less of a problem and a slightly reduced
       ! power_bias should be beneficial.
       power_bias_sbt = 1.d0
    case(OVLP_TYPE_ERS)
       add_k_power = 0.d0
       ! IGOR for explite range separation: Erfc[r12,rcut,widthfac]
       ! Field will be scaled by 1/k^2 cut to the left to a finite value.
       ! Therefore, large k are less of a problem and a slightly reduced
       ! power_bias should be beneficial.
       power_bias_sbt = 1.d0
    case default; call aims_stop('Invalid ovlp_type')
    end select
    power_bias_int = 0.d0
    power_bias_add = 0.5d0 * (3.d0 + add_k_power - power_bias_int) &
    &                - power_bias_sbt

    ! --- Get SBTs of product basis functions

    call sbt_atomic_transforms(N, lnr0, lnk0, lnrange, &
    &                   power_bias_sbt, power_bias_add, &
    &                   n_bas_fns, basfn_l, basfn_species, basfn_wave_spl, &
    &                   ffk)

    ! --- Prepare HSE potential if needed

    select case (ovlp_type)
    case(OVLP_TYPE_OVERLAP)
       continue
    case(OVLP_TYPE_COULOMB)
       ! The 1/k^2 part has already been done by add_k_power.
       ffk = sqrt(4.d0*pi) * ffk
    case(OVLP_TYPE_HSE)
       do i = 1, N
          k = exp(lnk0 + (i-1) * (lnrange/N))
          hse_fac = radial_fourier_hse(k, hse_omega_hf)
          ffk(i, :) = sqrt(hse_fac) * ffk(i, :)
       end do
    case(OVLP_TYPE_LR)
       do i = 1, N
          k = exp(lnk0 + (i-1) * (lnrange/N))
          hse_fac = radial_fourier_hse(k, hse_omega_hf)
          ffk(i, :) = sqrt(hse_fac) * ffk(i, :)
       end do
    case(OVLP_TYPE_CUT)
       if (.not. use_cutCb) then
          call aims_stop('Cut Coulomb operator not initialized', func)
       end if
       call initialize_cut_coulomb(N, lnr0, lnk0, lnrange, &
       &                           cutCb_width, cutCb_rcut)
       do i = 1, N
          ffk(i, :) = sqrt((2*pi)**1.5d0 * cutCb_vvk(i)) * ffk(i, :)
       end do
       call deallocate_cut_coulomb()
    case(OVLP_TYPE_CUT_ANALYTIC)
        ! IGOR: analytic cut coulomb for test
       if (.not. use_cutCb) then
          call aims_stop('Cut Coulomb operator not initialized', func)
       end if
       call initialize_cut_coulomb_analytic(N, lnr0, lnk0, lnrange, &
       &                           cutCb_width, cutCb_rcut)
       do i = 1, N
          ffk(i, :) = sqrt((2*pi)**1.5d0 * cutCb_vvk(i)) * ffk(i, :)
       end do
       call deallocate_cut_coulomb()
    case(OVLP_TYPE_ERS)
       call initialize_ers(N, lnr0, lnk0, lnrange, &
       &                           ers_width,ers_rcut, ers_omega)
       do i = 1, N
          ffk(i, :) = sqrt((2*pi)**1.5d0 * cutCb_vvk(i)) * ffk(i, :)
       end do
       call deallocate_cut_coulomb()
    case default
       call aims_stop('Invalid ovlp_type')
    end select

  end subroutine sbt_atomic_field_transforms
  !******
end module sbt_overlap_aims
!******
