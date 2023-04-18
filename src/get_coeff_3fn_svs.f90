!****s* FHI-aims/get_coeff_3fn_svs
!  NAME
!    get_coeff_3fn_svs
!  SYNOPSIS

subroutine get_coeff_3fn_svs(ovlp_3fn)

  !  PURPOSE
  !
  !    Calculate the orthonormalized product expansion coefficients in RI-SVS:
  !
  !       O_{ij}^\mu := <\phi_i \phi_j . P_\mu> S_{\mu\nu}^-1 V_{\nu\lamba}^.5
  !
  !  USES

  use timing
  use dimensions
  use mpi_tasks, only: check_allocation
  use prodbas
  use runtime_choices, only: use_scalapack
  use species_data
  use sbt_overlap_aims
  implicit none

  !  ARGUMENTS

  real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)

  !  INPUTS
  !    none
  !  OUTPUTS
  !    ovlp_3fn -- product expansion coefficients.
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

  real*8, allocatable :: coulomb_matr(:,:), ovlp_prodbas(:,:)
  integer :: info
  character(*), parameter :: func = 'get_coeff_3fn_svs'
! Timing

  ! Integrate ovlp_3fn (T)
  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn)
  call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, OVLP_TYPE_OVERLAP)
  call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )


  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr )
  ! Integrate V
  allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
  else
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
  end if
  ! Integrate S
  allocate(ovlp_prodbas(n_basbas,n_loc_prodbas), stat=info)
  call check_allocation(info, 'ovlp_prodbas')
  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(ovlp_prodbas, OVLP_TYPE_OVERLAP, .false.)
  else
     call integrate_ovlp_prodbas(l_shell_max, ovlp_prodbas)
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)

  ! S^-1 V^0.5
  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if(use_scalapack) then
     call evaluate_invs_times_sqrtv_scalapack(ovlp_prodbas, coulomb_matr)
  else
     call evaluate_invs_times_sqrtv(ovlp_prodbas, coulomb_matr)
  endif
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr, tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

  ! C := T (S^-1 V^0.5)
  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_ovlp3fn_multi_sv(coulomb_matr, ovlp_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

  deallocate(coulomb_matr, ovlp_prodbas)

end subroutine get_coeff_3fn_svs
!******
