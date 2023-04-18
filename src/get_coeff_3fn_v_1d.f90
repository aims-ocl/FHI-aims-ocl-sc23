!****s* FHI-aims/get_coeff_3fn_v_1d
!  NAME
!    get_coeff_3fn_v_1d
!  SYNOPSIS

subroutine get_coeff_3fn_v_1d(ovlp_3fn)

  !  PURPOSE
  !
  !    Calculate the orthonormalized product expansion coefficients in RI-V:
  !
  !       C_{ij}^\mu := (\phi_i \phi_j | P_\mu) V_{\mu\nu}^-0.5
  !
  !    Uses 1d distributed matrices internally.
  !
  !  USES

  use timing
  use localorb_io, only: localorb_info
  use prodbas
  use dimensions
  use runtime_choices
  use sbt_overlap_aims
  use species_data
  use mpi_tasks, only: check_allocation
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
  integer :: info
  character(*), parameter :: func = 'get_coeff_3fn_v_1d'
!!
!! for testing
  integer :: i, j
!!


  ! Integrate ovlp3fn (T)
  call get_timestamps(time_ovlp3fn, clock_time_ovlp3fn )
  call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, ovlp_type_bare_or_hse_coul)
  call get_times(time_ovlp3fn, clock_time_ovlp3fn, tot_time_ovlp3fn, tot_clock_time_ovlp3fn )

!!
!  write(use_unit,*) "ovlp3fn_matr:"
!  do i = 1,n_basis_pairs
!    do j = 1, n_loc_prodbas
!        write(use_unit,*) i, j, ovlp_3fn(i,j)
!    end do
!  end do

  ! Integrate V
  call get_timestamps(time_coulomb_matr, clock_time_coulomb_matr)
  allocate(coulomb_matr(n_basbas, n_loc_prodbas), stat=info)
  call check_allocation(info, 'coulomb_matr', func)
  if (use_logsbt) then
     call integrate_auxmat_by_atomic_sbt(coulomb_matr, ovlp_type_bare_or_hse_coul, .false.)
     !do i=1, n_basbas/2
     ! write(use_unit,'(I4,2f18.6)')i, coulomb_matr(i,i), coulomb_matr(i+n_basbas/2,i)
     !enddo
  else
     call integrate_coulomb_matr_v0(l_shell_max, coulomb_matr)
  end if
  call get_times(time_coulomb_matr, clock_time_coulomb_matr, tot_time_coulomb_matr, tot_clock_time_coulomb_matr)

!!
!  write(use_unit,*) "Coulomb_matr:"
!  do i = 1,n_basbas
!    do j = 1, n_loc_prodbas
!        write(use_unit,*) i, j, coulomb_matr(i,j)
!    end do
!  end do

!!
  ! V^-0.5
  call get_timestamps(time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
  if (use_asym_inv_sqrt) then
     ! transposed==.true. because get_v_multi_ovlp3fn() does
     !   C := V^-0.5T.
     call asym_inv_sqrt_of_auxmat_scalapack(coulomb_matr, "Coulomb", .true.)
  else if (use_scalapack) then
     call power_auxmat_scalapack(coulomb_matr, -0.5d0, "Coulomb")
  else
     call power_auxmat_lapack(coulomb_matr, -0.5d0, "Coulomb")
  endif
  call get_times(time_inv_coulomb_matr, clock_time_inv_coulomb_matr, tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)

  ! C := T V^-0.5
  call get_timestamps(time_ovlp_multi, clock_time_ovlp_multi)
  call get_v_multi_ovlp3fn(coulomb_matr, ovlp_3fn)
  call get_times(time_ovlp_multi, clock_time_ovlp_multi, tot_time_ovlp_multi, tot_clock_time_ovlp_multi)

  deallocate(coulomb_matr)

  call localorb_info('')

end subroutine get_coeff_3fn_v_1d
!******
