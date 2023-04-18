!****s* FHI-aims/get_coeff_3fn_v_2d
!  NAME
!    get_coeff_3fn_v_2d
!  SYNOPSIS

subroutine get_coeff_3fn_v_2d(ovlp_3fn)

  !  PURPOSE
  !
  !    Calculate the orthonormalized product expansion coefficients in RI-V:
  !
  !       C_{ij}^\mu := (\phi_i \phi_j | P_\mu) V_{\mu\nu}^-0.5
  !
  !    Uses 2d distributed matrices internally.  The result should be
  !    identical to get_coeff_3fn_v_1d(), but this routine should be faster.
  !    It depends on scalapack, though.
  !
  !  USES

  use timing
  use prodbas
  use runtime_choices
  use dimensions
  use sbt_overlap_aims
  use species_data
  use localorb_io, only: localorb_info
  use mpi_tasks, only: aims_stop, check_allocation
  implicit none

  !  ARGUMENTS

  real*8, intent(OUT) :: ovlp_3fn(n_basis_pairs, n_loc_prodbas)

  !  INPUTS
  !    none
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

  real*8 :: t0
  real*8, allocatable :: coulomb_matr(:,:), coulomb_matr_1d(:,:)
  integer :: info
  character(*), parameter :: func = 'get_coeff_3fn_v_2d'

  if(.not.use_scalapack) then
     call aims_stop('Cannot use 2d distribution without scalapack', func)
  endif

  ! Integrate ovlp3fn (T)
  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
  call integrate_ovlp3fn(l_shell_max, ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
  call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )

  ! Integrate V
  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
  allocate(coulomb_matr(max_row_2d,max_col_2d), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .true.)
  else
     allocate(coulomb_matr_1d(n_basbas, n_loc_prodbas), stat=info)
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr_1d)
     ! t0 = mpi_wtime()
     call dist_1d_2d(n_basbas, coulomb_matr_1d, ubound(coulomb_matr_1d, 1), &
     &                         coulomb_matr, ubound(coulomb_matr,1))
     ! if(myid==0) print *,'dist_1d_2d coulomb_matr:',mpi_wtime()-t0
     deallocate(coulomb_matr_1d)
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)

  ! V^-0.5
  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if (use_asym_inv_sqrt) then
     ! transposed==.false. because get_v_multi_ovlp3fn_2() does
     !   C := TV^-0.5.
     call asym_inv_sqrt_of_auxmat_scalapack_2d(coulomb_matr, "Coulomb",.false.)
  else
     call power_auxmat_scalapack_2d(coulomb_matr, -0.5d0, "Coulomb")
  end if
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr, tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

  ! C := T V^-0.5
  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_v_multi_ovlp3fn_2(coulomb_matr, ovlp_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

  deallocate(coulomb_matr)

  call localorb_info('')

end subroutine get_coeff_3fn_v_2d
!******
