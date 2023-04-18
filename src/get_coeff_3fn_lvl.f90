!****s* FHI-aims/get_coeff_3fn_lvl
!  NAME
!    get_coeff_3fn_lvl
!  SYNOPSIS

subroutine get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr)

  !  PURPOSE
  !
  !    Calculate the NON-orthonormalized product expansion coefficients in
  !    RI-LVL:
  !
  !       C_{ij}^\mu := (\phi_i \phi_j | P_\mu) L_{\mu\nu}^-1,
  !
  !    where L is the atom-pair local block of the Coulomb matrix.  The
  !    auxiliary basis is not orthonormalized because we would loose all of
  !    the sparsity, then.
  !
  !  USES

  use timing
  use prodbas
  use runtime_choices
  use sbt_overlap_aims
  use sparse_tensor
  use species_data
  use localized_basbas
  use localorb_io
  implicit none

  !  ARGUMENTS

  type(sp_ten), intent(OUT) :: coeff_3fn_ten
  real*8, intent(OUT) :: coulomb_matr(n_basbas, n_loc_prodbas)

  !  INPUTS
  !    none
  !  OUTPUTS
  !    o coeff_3fn_ten -- Expansion coefficients
  !    o coulomb_matr -- Coulomb matrix (because basbas is non-orthogonal)
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

  type(sp_ten) :: coulomb_ten
  type(sp_ten) :: want, ovlp3fn_tmp, coulomb_tmp
  real*8 :: n_bytes_o3fn, n_bytes_lcl, n_bytes_fcl, n_bytes_tot, n_bytes_full
  real*8 :: time_int2(4), time_dist2(4), time_int3(4), time_dist3(4)
  integer :: info
  real*8, allocatable :: ovlp_3fn(:,:)
  character*150 :: info_str
  character(*), parameter :: func = 'get_coeff_3fn_lvl'

  time_int3 = 0.d0; time_dist3 = 0.d0
  time_int2 = 0.d0; time_dist2 = 0.d0

  call initialize_localized_basbas()

  ! --- ovlp3fn_ten

  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn)
  if (.not. use_logsbt_lvltriples) then
     allocate(ovlp_3fn(n_basis_pairs, n_loc_prodbas), stat=info)
     call check_allocation(info, 'ovlp_3fn', func)
     ! Integrate
     call start_timer(time_int3)
     call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
     call stop_timer(time_int3)

     ! Distribute
     call localorb_info('  Redistribute 3-center integrals')
     call start_timer(time_dist3)
     call not_so_sparse_ovlp3fn(ovlp_3fn, ovlp3fn_tmp)
     deallocate(ovlp_3fn)
     call construct_LVL_3fn_want(want)
     call redistribute_sp_ten(ovlp3fn_tmp, want, '3fn-LVL')
     call dealloc_sp_ten(want)
     call simplified_tensor(ovlp3fn_tmp, ovlp3fn_tmp%val, &
     &                      coeff_3fn_ten, 'coeff3fn')
     call dealloc_sp_ten(ovlp3fn_tmp)
     call stop_timer(time_dist3)
     call localorb_info('')
  else
     ! Get local needs
     call construct_LVL_3fn_want(coeff_3fn_ten)
     ! Allocate corresponding arrays
     call upgrade_blocking(coeff_3fn_ten, SPTEN_FULL, 'coeff3fn')
  end if
  call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn)

  call debug_global_size(coeff_3fn_ten, n_bytes_o3fn)
  n_bytes_fcl = 8.d0 * n_basbas * n_max_loc_prodbas
  n_bytes_tot = n_bytes_o3fn + n_bytes_fcl
  n_bytes_full = 8.d0 * n_basis_pairs * n_max_loc_prodbas
  write(info_str, "(A,I7,' MiB x',I5,' procs =',F12.3,' GiB')") &
  & 'The resident memory for O3fn + V2fn is ', nint(n_bytes_tot/2.d0**20), &
  & n_tasks, n_tasks * n_bytes_tot/2.d0**30
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "('| ',A,I7,' MiB x',I5,' procs =',F12.3,' GiB')") &
  & 'Part O3fn (temporarily alloced  2x) is', nint(n_bytes_o3fn/2.d0**20), &
  & n_tasks, n_tasks * n_bytes_o3fn/2.d0**30
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "('| ',A,I7,' MiB x',I5,' procs =',F12.3,' GiB')") &
  & 'Part V2fn (temporarily alloced ~3x) is', nint(n_bytes_fcl/2.d0**20), &
  & n_tasks, n_tasks * n_bytes_fcl/2.d0**30
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  write(info_str, "('| ',A,I7,' MiB x',I5,' procs =',F12.3,' GiB')") &
  & 'Original full ovlp_3fn would have been', nint(n_bytes_full/2.d0**20), &
  & n_tasks, n_tasks * n_bytes_full/2.d0**30
  call localorb_info(info_str, use_unit, '(2X,A)', OL_norm)
  call localorb_info('')

  ! --- coulomb_matr

  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)

  ! Integrate
  call start_timer(time_int2)
  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
  else
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
     ! JW: We should at least symmetrize Coulomb matrix.
     call localorb_info("  * Grid Coulomb matrix for LVL not recommended.")
     call localorb_info("  * Please ask for this at least to be symmetrized.")
  end if
  call localorb_info('')
  call stop_timer(time_int2)

  ! --- coulomb_ten

  if (.not. use_logsbt_lvltriples) then
     ! Distribute for coeff3fn_ten (inversion)
     call start_timer(time_dist2)
     call localorb_info('  Distribute Coulomb matrix')
     call localorb_info('')
     call get_auxmat_descriptor(.false., coulomb_tmp, '2fn-1d_prodbas')
     call construct_LVL_coulomb_want(want)
     call redistribute_sp_ten(coulomb_tmp, want, '2fn-LVL-tmp', coulomb_matr)
     call dealloc_sp_ten(want)
     call simplified_tensor(coulomb_tmp, coulomb_tmp%val, coulomb_ten, &
     &                      '2fn-LVL')
     call dealloc_sp_ten(coulomb_tmp)
     call stop_timer(time_dist2)
     call debug_global_size(coulomb_ten, n_bytes_lcl)

     call localorb_info('')
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)

  ! --- coeff3fn_ten

  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)

  call start_timer(time_int3)
  call localorb_info('  Get 3-center expansion coefficients')
  call localorb_info('')
  if (use_logsbt_lvltriples) then
     call get_coeff3fn_ten(coeff_3fn_ten)
  else
     ! coeff_3fn := ovlp_3fn * loc_coulomb_matr^(-1)
     call ovlp3fn_to_coeff3fn_ten(coeff_3fn_ten, coulomb_ten)
     call dealloc_sp_ten(coulomb_ten)
  end if
  call stop_timer(time_int3)

  ! --- redistribute to basbas blocking

  ! JW: Could use 2D blocking, here!
  ! Do this by not using the symmetry in the first two indices
  ! and 2D distributing over one basis index and the basbas index.

  ! But do not forget to adjust evaluate_exchange_matr_LVL_eigen()!
  ! This should speed up fnKSbb -> fnKScl considerably.

  call start_timer(time_dist3)
  call localorb_info('  Redistribute 3-center expansion coefficients')
  call localorb_info('')
  call get_ovlp3fn_basbas_blocking(want)
  call redistribute_sp_ten(coeff_3fn_ten, want, 'final_coeff3fn')
  call dealloc_sp_ten(want)
  call stop_timer(time_dist3)

  call get_times(time_ovlp_multi, clock_time_ovlp_multi,tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

  call localorb_info('')
  call output_timeheader('2X', 'LVL product basis initialization')
  call output_timer('3-center integration', time_int3(3:4))
  call output_timer('3-center distribution', time_dist3(3:4))
  call output_timer('2-center integration', time_int2(3:4))
  call output_timer('2-center distribution', time_dist2(3:4))
  call localorb_info('')

end subroutine get_coeff_3fn_lvl
!******
