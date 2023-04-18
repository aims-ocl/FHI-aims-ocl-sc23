!****s* FHI-aims/get_coeff_3fn_lvl_full
!  NAME
!    get_coeff_3fn_lvl_full
!  SYNOPSIS

subroutine get_coeff_3fn_lvl_full(ovlp_3fn)

  !  PURPOSE
  !
  !    Calculate the orthonormalized product expansion coefficients in RI-LVL:
  !
  !       M := T L V^0.5
  !
  !    The resulting ovlp_3fn is not sparse, anymore.  This subroutine is
  !    therefore FOR TESTING PURPOSES ONLY.  Usage in production would not be
  !    sensible.
  !
  !  USES

  use timing
  use prodbas
  use sbt_overlap_aims, only: integrate_auxmat_by_atomic_sbt
  use sparse_tensor
  use localized_basbas
  implicit none

  !  ARGUMENTS

  real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)

  !  INPUTS
  !    none [atomic positions and basis/basbas settings from modules]
  !  OUTPUTS
  !    o ovlp_3fn -- product expansion coefficients.
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

  real*8, allocatable :: coulomb_matr(:,:)
  type(sp_ten) :: coeff_3fn_ten
  integer :: info
  character(*), parameter :: func = 'get_coeff_3fn_lvl_full'

  ! --- RI-LVL calculation

  allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr)

  ! --- To dense matrix

  call back_to_ovlp3fn(ovlp_3fn, coeff_3fn_ten)

  call dealloc_sp_ten(coeff_3fn_ten)

  ! --- Multiply with V^0.5

  if (use_ers .or. use_erfc) then
      call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_ers_coul, .false.)
  endif

  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if (use_scalapack) then
     call power_auxmat_scalapack(coulomb_matr, 0.5d0, 'Coulomb')
  else
     call power_auxmat_lapack(coulomb_matr, 0.5d0, 'Coulomb')
  end if
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr, tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_v_multi_ovlp3fn(coulomb_matr, ovlp_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

  ! --- Tidy up

  deallocate(coulomb_matr)

end subroutine get_coeff_3fn_lvl_full
!******
