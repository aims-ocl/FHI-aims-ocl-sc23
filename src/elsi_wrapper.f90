!****h* FHI-aims/elsi_wrapper
!  NAME
!    elsi_wrapper
!  SYNOPSIS
module elsi_wrapper
!  PURPOSE
!    This module provides interfaces to the ELSI library.
!  AUTHOR
!    Victor Yu, Duke University, FHI-aims team
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is
!    subject to the terms and conditions of the respective license
!    agreement.
!  SOURCE

  use dimensions, only: n_states,n_core_states,n_valence_basis,n_basis,&
      n_k_points,use_molecular_dynamics
  use elsi
  use free_atoms, only: free_wave_eigenval
  use localorb_io, only: use_unit
  use mpi_tasks, only: myid,mpi_comm_global
  use pbc_lists, only: k_weights
  use physics, only: n_electrons
  use rel_x2c_mod, only: dim_matrix_rel
  use runtime_choices, only: elsi_out_level,elsi_out_json,elsi_solver,&
      elsi_elpa_solver,elsi_elpa_n_single,elsi_elpa_gpu,elsi_omm_n_elpa,&
      elsi_omm_flavor,elsi_omm_tol,elsi_pexsi_np_symbo,elsi_pexsi_np_per_pole,&
      elsi_eigenexa_method,elsi_sips_n_elpa,elsi_sips_n_slice,&
      elsi_ntpoly_method,elsi_ntpoly_tol,elsi_ntpoly_filter,elsi_magma_solver,&
      elsi_filter,elsi_rw_filter,first_elsi_call,out_h_elsi,out_s_elsi,&
      use_scalapack,frozen_core_scf,relax_mode,RELAX_OFF,basis_threshold,&
      occupation_type,occupation_width,occupation_acc,n_methfessel_paxton,&
      fixed_spin_moment,flag_rel,REL_x2c,REL_4c_dks,dry_run

  implicit none

  private

  public :: aims_elsi_init_scf
  public :: aims_elsi_reinit_scf
  public :: aims_elsi_finalize_scf
  public :: aims_elsi_set_output
  public :: aims_elsi_set_illcond_check
  public :: aims_elsi_set_pexsi_mu_min
  public :: aims_elsi_set_pexsi_mu_max
  public :: aims_elsi_set_sips_ev_min
  public :: aims_elsi_set_sips_ev_max
  public :: aims_elsi_set_mu_broaden_width
  public :: aims_elsi_set_mu_spin_degen
  public :: aims_elsi_get_n_illcond
  public :: aims_elsi_get_edm
  public :: aims_elsi_get_ovlp_ev_min
  public :: aims_elsi_get_ovlp_ev_max
  public :: aims_elsi_get_pexsi_mu_min
  public :: aims_elsi_get_pexsi_mu_max
  public :: aims_elsi_ev
  public :: aims_elsi_dm
  public :: aims_elsi_orthonormalize_ev
  public :: aims_elsi_extrapolate_dm
  public :: aims_elsi_compute_mu_and_occ
  public :: aims_elsi_occ
  public :: aims_elsi_read_mat_dim
  public :: aims_elsi_read_mat
  public :: aims_elsi_write_mat
  public :: aims_elsi_stdevp

  type(elsi_handle), public :: eh_scf ! For self-consistent field
  type(elsi_rw_handle), public :: rwh_r ! For reading matrices
  type(elsi_rw_handle), public :: rwh_w ! For writing matrices

  real(8), public :: pexsi_mu_min
  real(8), public :: pexsi_mu_max

  integer :: lrow ! Local matrix size
  integer :: lcol ! Local matrix size

  interface aims_elsi_get_edm
    module procedure aims_elsi_get_edm_real
    module procedure aims_elsi_get_edm_cmplx
  end interface

  interface aims_elsi_ev
    module procedure aims_elsi_ev_real
    module procedure aims_elsi_ev_cmplx
  end interface

  interface aims_elsi_dm
    module procedure aims_elsi_dm_real
    module procedure aims_elsi_dm_cmplx
  end interface

  interface aims_elsi_orthonormalize_ev
    module procedure aims_elsi_orthonormalize_ev_real
    module procedure aims_elsi_orthonormalize_ev_cmplx
  end interface

  interface aims_elsi_extrapolate_dm
    module procedure aims_elsi_extrapolate_dm_real
    module procedure aims_elsi_extrapolate_dm_cmplx
  end interface

  interface aims_elsi_read_mat
    module procedure aims_elsi_read_mat_real
    module procedure aims_elsi_read_mat_cmplx
  end interface

  interface aims_elsi_write_mat
    module procedure aims_elsi_write_mat_real
    module procedure aims_elsi_write_mat_cmplx
  end interface

  interface aims_elsi_stdevp
    module procedure aims_elsi_stdevp_real
    module procedure aims_elsi_stdevp_cmplx
  end interface

contains

!****s* FHI-aims/aims_elsi_init_scf
!  NAME
!   aims_elsi_init_scf
!  SYNOPSIS
subroutine aims_elsi_init_scf(elsi_comm,elsi_ctxt,blk,n_lrow,n_lcol,i_kpt)

  implicit none

  integer, intent(in) :: elsi_comm
  integer, intent(in) :: elsi_ctxt
  integer, intent(in) :: blk
  integer, intent(in) :: n_lrow
  integer, intent(in) :: n_lcol
  integer, intent(in) :: i_kpt

  real(8) :: est_width
  integer :: elsi_n_basis
  integer :: elsi_n_states
  integer :: elsi_mode

  lrow = n_lrow
  lcol = n_lcol
  first_elsi_call = .true.
  out_h_elsi = .false.
  out_s_elsi = .false.

  if(frozen_core_scf) then
     elsi_n_basis = n_valence_basis
     elsi_n_states = n_states-n_core_states
  else
     if(flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) then
        elsi_n_basis = 2*dim_matrix_rel
        elsi_n_states = 2*dim_matrix_rel
     else
        elsi_n_basis = n_basis
        elsi_n_states = n_states
     end if
  end if

  if(use_scalapack) then ! MULTI_PROC mode
     elsi_mode = 1
  else ! SINGLE_PROC mode
     elsi_mode = 0
  end if

  call elsi_init(eh_scf,elsi_solver,elsi_mode,0,elsi_n_basis,n_electrons,&
       elsi_n_states)
  call elsi_init_rw(rwh_r,0,elsi_mode,n_basis,n_electrons)
  call elsi_init_rw(rwh_w,1,elsi_mode,n_basis,n_electrons)

  if(use_scalapack) then ! MULTI_PROC mode
     call elsi_set_kpoint(eh_scf,n_k_points,i_kpt,k_weights(i_kpt))
     call elsi_set_mpi(eh_scf,elsi_comm)
     call elsi_set_blacs(eh_scf,elsi_ctxt,blk)
     call elsi_set_mpi_global(eh_scf,mpi_comm_global)
     call elsi_set_rw_mpi(rwh_r,elsi_comm)
     call elsi_set_rw_blacs(rwh_r,elsi_ctxt,blk)
     call elsi_set_rw_mpi(rwh_w,elsi_comm)
     call elsi_set_rw_blacs(rwh_w,elsi_ctxt,blk)
  end if

  ! Output options
  if(myid == 0) then
     call elsi_set_output_log(eh_scf,elsi_out_json)
     call elsi_set_output(eh_scf,elsi_out_level)
  else
     call elsi_set_output_log(eh_scf,0)
     call elsi_set_output(eh_scf,0)
  end if

  call elsi_set_output_unit(eh_scf,use_unit)

  ! Check ill-conditioning
  call elsi_set_illcond_check(eh_scf,1)

  if(relax_mode /= RELAX_OFF .or. use_molecular_dynamics) then
     call elsi_set_save_ovlp(eh_scf,1)
  end if

  ! Singularity threshold
  call elsi_set_illcond_tol(eh_scf,basis_threshold)

  ! Numerical zero threshold
  call elsi_set_zero_def(eh_scf,elsi_filter)
  call elsi_set_rw_zero_def(rwh_r,elsi_rw_filter)
  call elsi_set_rw_zero_def(rwh_w,elsi_rw_filter)

  ! Broadening scheme and width
  call elsi_set_mu_broaden_scheme(eh_scf,occupation_type)
  if(occupation_type == 4) then ! Cubic polynomial
     call elsi_set_mu_broaden_scheme(eh_scf,3)
  else if(occupation_type == 5) then ! Cold
     call elsi_set_mu_broaden_scheme(eh_scf,4)
  end if
  call elsi_set_mu_broaden_width(eh_scf,occupation_width)
  call elsi_set_mu_tol(eh_scf,occupation_acc)
  if(occupation_type == 2) then ! Methfessel-Paxton
     call elsi_set_mu_mp_order(eh_scf,n_methfessel_paxton)
  end if
  if(flag_rel == REL_x2c .or. flag_rel == REL_4c_dks) then
    call elsi_set_mu_spin_degen(eh_scf,1.d0)
  end if

  ! Solver settings
  select case(elsi_solver)
  case(1)
     call elsi_set_elpa_solver(eh_scf,elsi_elpa_solver)
     call elsi_set_elpa_gpu(eh_scf,elsi_elpa_gpu)
     call elsi_set_elpa_n_single(eh_scf,elsi_elpa_n_single)
  case(2)
     call elsi_set_omm_n_elpa(eh_scf,elsi_omm_n_elpa)
     call elsi_set_omm_flavor(eh_scf,elsi_omm_flavor)
     call elsi_set_omm_tol(eh_scf,elsi_omm_tol)
  case(3)
     call elsi_set_pexsi_temp(eh_scf,occupation_width)
     call elsi_set_pexsi_np_symbo(eh_scf,elsi_pexsi_np_symbo)
     if(elsi_pexsi_np_per_pole > 0) then
        call elsi_set_pexsi_np_per_pole(eh_scf,elsi_pexsi_np_per_pole)
     end if
     est_width = abs(minval(free_wave_eigenval))
     call elsi_set_pexsi_delta_e(eh_scf,est_width)
  case(4)
     call elsi_set_eigenexa_method(eh_scf,elsi_eigenexa_method)
  case(5)
     call elsi_set_sips_n_slice(eh_scf,elsi_sips_n_slice)
     call elsi_set_sips_n_elpa(eh_scf,elsi_sips_n_elpa)
  case(6)
     call elsi_set_ntpoly_method(eh_scf,elsi_ntpoly_method)
     call elsi_set_ntpoly_tol(eh_scf,elsi_ntpoly_tol)
     call elsi_set_ntpoly_filter(eh_scf,elsi_ntpoly_filter)
  case(7)
     call elsi_set_magma_solver(eh_scf,elsi_magma_solver)
  end select

end subroutine
!******
!****s* FHI-aims/aims_elsi_reinit_scf
!  NAME
!   aims_elsi_reinit_scf
!  SYNOPSIS
subroutine aims_elsi_reinit_scf()

  implicit none

  first_elsi_call = .true.
  out_h_elsi = .false.
  out_s_elsi = .false.

  call elsi_reinit(eh_scf)

end subroutine
!******
!****s* FHI-aims/aims_elsi_finalize_scf
!  NAME
!   aims_elsi_finalize_scf
!  SYNOPSIS
subroutine aims_elsi_finalize_scf()

  implicit none

  if(.not. dry_run) then
     call elsi_finalize(eh_scf)
     call elsi_finalize_rw(rwh_r)
     call elsi_finalize_rw(rwh_w)
  end if

end subroutine
!******
!****s* FHI-aims/aims_elsi_ev_real
!  NAME
!   aims_elsi_ev_real
!  SYNOPSIS
subroutine aims_elsi_ev_real(eh,ham,ovlp,eval,evec)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: eval(eh%ph%n_basis)
  real(8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_ev_real(eh,ham,ovlp,eval,evec)

end subroutine
!******
!****s* FHI-aims/aims_elsi_ev_cmplx
!  NAME
!   aims_elsi_ev_cmplx
!  SYNOPSIS
subroutine aims_elsi_ev_cmplx(eh,ham,ovlp,eval,evec)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  complex(8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)
  complex(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: eval(eh%ph%n_basis)
  complex(8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_ev_complex(eh,ham,ovlp,eval,evec)

end subroutine
!******
!****s* FHI-aims/aims_elsi_dm_real
!  NAME
!   aims_elsi_dm_real
!  SYNOPSIS
subroutine aims_elsi_dm_real(eh,ham,ovlp,dm,energy)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: energy

  call elsi_dm_real(eh,ham,ovlp,dm,energy)

end subroutine
!******
!****s* FHI-aims/aims_elsi_dm_cmplx
!  NAME
!   aims_elsi_dm_cmplx
!  SYNOPSIS
subroutine aims_elsi_dm_cmplx(eh,ham,ovlp,dm,energy)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  complex(8), intent(inout) :: ham(eh%bh%n_lrow,eh%bh%n_lcol)
  complex(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  complex(8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: energy

  call elsi_dm_complex(eh,ham,ovlp,dm,energy)

end subroutine
!******
!****s* FHI-aims/aims_elsi_orthonormalize_ev_real
!  NAME
!   aims_elsi_orthonormalize_ev_real
!  SYNOPSIS
subroutine aims_elsi_orthonormalize_ev_real(eh,ovlp,evec)

  implicit none

  type(elsi_handle), intent(in) :: eh
  real(8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_orthonormalize_ev_real(eh,ovlp,evec)

end subroutine
!******
!****s* FHI-aims/aims_elsi_orthonormalize_ev_cmplx
!  NAME
!   aims_elsi_orthonormalize_ev_cmplx
!  SYNOPSIS
subroutine aims_elsi_orthonormalize_ev_cmplx(eh,ovlp,evec)

  implicit none

  type(elsi_handle), intent(in) :: eh
  complex(8), intent(in) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  complex(8), intent(inout) :: evec(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_orthonormalize_ev_complex(eh,ovlp,evec)

end subroutine
!******
!****s* FHI-aims/aims_elsi_extrapolate_dm_real
!  NAME
!   aims_elsi_extrapolate_dm_real
!  SYNOPSIS
subroutine aims_elsi_extrapolate_dm_real(eh,ovlp,dm)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  real(8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_extrapolate_dm_real(eh,ovlp,dm)

end subroutine
!******
!****s* FHI-aims/aims_elsi_extrapolate_dm_cmplx
!  NAME
!   aims_elsi_extrapolate_dm_cmplx
!  SYNOPSIS
subroutine aims_elsi_extrapolate_dm_cmplx(eh,ovlp,dm)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  complex(8), intent(inout) :: ovlp(eh%bh%n_lrow,eh%bh%n_lcol)
  complex(8), intent(inout) :: dm(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_extrapolate_dm_complex(eh,ovlp,dm)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_edm_real
!  NAME
!   aims_elsi_get_edm_real
!  SYNOPSIS
subroutine aims_elsi_get_edm_real(eh,edm)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_get_edm_real(eh,edm)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_edm_cmplx
!  NAME
!   aims_elsi_get_edm_cmplx
!  SYNOPSIS
subroutine aims_elsi_get_edm_cmplx(eh,edm)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  complex(8), intent(out) :: edm(eh%bh%n_lrow,eh%bh%n_lcol)

  call elsi_get_edm_complex(eh,edm)

end subroutine
!******
!******
!****s* FHI-aims/aims_elsi_set_pexsi_mu_min
!  NAME
!   aims_elsi_set_pexsi_mu_min
!  SYNOPSIS
subroutine aims_elsi_set_pexsi_mu_min(eh,mu_min)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: mu_min

  call elsi_set_pexsi_mu_min(eh,mu_min)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_pexsi_mu_max
!  NAME
!   aims_elsi_set_pexsi_mu_max
!  SYNOPSIS
subroutine aims_elsi_set_pexsi_mu_max(eh,mu_max)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: mu_max

  call elsi_set_pexsi_mu_max(eh,mu_max)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_sips_ev_min
!  NAME
!   aims_elsi_set_sips_ev_min
!  SYNOPSIS
subroutine aims_elsi_set_sips_ev_min(eh,ev_min)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: ev_min

  call elsi_set_sips_ev_min(eh,ev_min)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_sips_ev_max
!  NAME
!   aims_elsi_set_sips_ev_max
!  SYNOPSIS
subroutine aims_elsi_set_sips_ev_max(eh,ev_max)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: ev_max

  call elsi_set_sips_ev_max(eh,ev_max)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_output
!  NAME
!   aims_elsi_set_output
!  SYNOPSIS
subroutine aims_elsi_set_output(eh,out_level)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  integer, intent(in) :: out_level

  call elsi_set_output(eh,out_level)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_illcond_check
!  NAME
!   aims_elsi_set_illcond_check
!  SYNOPSIS
subroutine aims_elsi_set_illcond_check(eh,illcond_check)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  integer, intent(in) :: illcond_check

  call elsi_set_illcond_check(eh,illcond_check)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_mu_broaden_width
!  NAME
!   aims_elsi_set_mu_broaden_width
!  SYNOPSIS
subroutine aims_elsi_set_mu_broaden_width(eh,broaden_width)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: broaden_width

  call elsi_set_mu_broaden_width(eh,broaden_width)

end subroutine
!******
!****s* FHI-aims/aims_elsi_set_mu_spin_degen
!  NAME
!   aims_elsi_set_mu_spin_degen
!  SYNOPSIS
subroutine aims_elsi_set_mu_spin_degen(eh,spin_degen)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: spin_degen

  call elsi_set_mu_spin_degen(eh,spin_degen)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_n_illcond
!  NAME
!   aims_elsi_get_n_illcond
!  SYNOPSIS
subroutine aims_elsi_get_n_illcond(eh,n_illcond)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  integer, intent(out) :: n_illcond

  call elsi_get_n_illcond(eh,n_illcond)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_olvp_ev_min
!  NAME
!   aims_elsi_get_ovlp_ev_min
!  SYNOPSIS
subroutine aims_elsi_get_ovlp_ev_min(eh,ev_min)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(out) :: ev_min

  call elsi_get_ovlp_ev_min(eh,ev_min)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_olvp_ev_max
!  NAME
!   aims_elsi_get_ovlp_ev_max
!  SYNOPSIS
subroutine aims_elsi_get_ovlp_ev_max(eh,ev_max)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(out) :: ev_max

  call elsi_get_ovlp_ev_max(eh,ev_max)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_pexsi_mu_min
!  NAME
!   aims_elsi_get_pexsi_mu_min
!  SYNOPSIS
subroutine aims_elsi_get_pexsi_mu_min(eh,mu_min)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(out) :: mu_min

  call elsi_get_pexsi_mu_min(eh,mu_min)

end subroutine
!******
!****s* FHI-aims/aims_elsi_get_pexsi_mu_max
!  NAME
!   aims_elsi_get_pexsi_mu_max
!  SYNOPSIS
subroutine aims_elsi_get_pexsi_mu_max(eh,mu_max)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(out) :: mu_max

  call elsi_get_pexsi_mu_max(eh,mu_max)

end subroutine
!******
!****s* FHI-aims/aims_elsi_stdevp_real
!  NAME
!   aims_elsi_stdevp_real
!  SYNOPSIS
subroutine aims_elsi_stdevp_real(ng,nlrow,nlcol,nev,mat,eval,evec,ctxt,blk,comm)

  implicit none

  integer, intent(in) :: ng
  integer, intent(in) :: nlrow
  integer, intent(in) :: nlcol
  integer, intent(in) :: nev
  real(8), intent(inout) :: mat(nlrow,nlcol)
  real(8), intent(out) :: eval(ng)
  real(8), intent(out) :: evec(nlrow,nlcol)
  integer, intent(in) :: ctxt
  integer, intent(in) :: blk
  integer, intent(in) :: comm

  real(8) :: dummy(1,1)
  type(elsi_handle) :: eh_tmp

  call elsi_init(eh_tmp,1,1,0,ng,0.0d0,nev)
  call elsi_set_mpi(eh_tmp,comm)
  call elsi_set_blacs(eh_tmp,ctxt,blk)
  call elsi_set_unit_ovlp(eh_tmp,1)
  call elsi_ev_real(eh_tmp,mat,dummy,eval,evec)
  call elsi_finalize(eh_tmp)

end subroutine
!******
!****s* FHI-aims/aims_elsi_stdevp_cmplx
!  NAME
!   aims_elsi_stdevp_cmplx
!  SYNOPSIS
subroutine aims_elsi_stdevp_cmplx(ng,nlrow,nlcol,nev,mat,eval,evec,ctxt,blk,comm)

  implicit none

  integer, intent(in) :: ng
  integer, intent(in) :: nlrow
  integer, intent(in) :: nlcol
  integer, intent(in) :: nev
  complex(8), intent(inout) :: mat(nlrow,nlcol)
  real(8), intent(out) :: eval(ng)
  complex(8), intent(out) :: evec(nlrow,nlcol)
  integer, intent(in) :: ctxt
  integer, intent(in) :: blk
  integer, intent(in) :: comm

  complex(8) :: dummy(1,1)
  type(elsi_handle) :: eh_tmp

  call elsi_init(eh_tmp,1,1,0,ng,0.0d0,nev)
  call elsi_set_mpi(eh_tmp,comm)
  call elsi_set_blacs(eh_tmp,ctxt,blk)
  call elsi_set_unit_ovlp(eh_tmp,1)
  call elsi_ev_complex(eh_tmp,mat,dummy,eval,evec)
  call elsi_finalize(eh_tmp)

end subroutine
!******
!****s* FHI-aims/aims_elsi_compute_mu_and_occ
!  NAME
!   aims_elsi_compute_mu_and_occ
!  SYNOPSIS
subroutine aims_elsi_compute_mu_and_occ(eh,n_electrons,n_states,n_spins,n_kpts,&
  k_wt,eval,occ,mu)

  implicit none

  type(elsi_handle), intent(inout) :: eh
  real(8), intent(in) :: n_electrons
  integer, intent(in) :: n_states
  integer, intent(in) :: n_spins
  integer, intent(in) :: n_kpts
  real(8), intent(in) :: k_wt(n_kpts)
  real(8), intent(in) :: eval(n_states,n_spins,n_kpts)
  real(8), intent(out) :: occ(n_states,n_spins,n_kpts)
  real(8), intent(out) :: mu

  call elsi_compute_mu_and_occ(eh,n_electrons,n_states,n_spins,n_kpts,k_wt,&
       eval,occ,mu)

end subroutine
!******
!****s* FHI-aims/aims_elsi_occ
!  NAME
!   aims_elsi_occ
!  SYNOPSIS
subroutine aims_elsi_occ(n_electrons,n_states,n_spins,n_kpts,k_wt,eval,occ,mu)

  implicit none

  real(8), intent(in) :: n_electrons
  integer, intent(in) :: n_states
  integer, intent(in) :: n_spins
  integer, intent(in) :: n_kpts
  real(8), intent(in) :: k_wt(n_kpts)
  real(8), intent(in) :: eval(n_states,n_spins,n_kpts)
  real(8), intent(out) :: occ(n_states,n_spins,n_kpts)
  real(8), intent(out) :: mu

  type(elsi_handle) :: aux_h

  call elsi_init(aux_h,elsi_solver,1,0,n_basis,n_electrons,n_states)

  call elsi_set_mu_broaden_scheme(aux_h,occupation_type)
  if(occupation_type == 2) then ! Methfessel-Paxton
     call elsi_set_mu_mp_order(aux_h,n_methfessel_paxton)
  end if
  call elsi_set_mu_broaden_width(aux_h,occupation_width)
  call elsi_set_mu_tol(aux_h,occupation_acc)
  if(fixed_spin_moment) then
     call elsi_set_mu_spin_degen(aux_h,1.0d0)
  end if

  call elsi_compute_mu_and_occ(aux_h,n_electrons,n_states,n_spins,n_kpts,k_wt,&
       eval,occ,mu)

  call elsi_finalize(aux_h)

end subroutine
!******
!****s* FHI-aims/aims_elsi_read_mat_dim
!  NAME
!   aims_elsi_read_mat_dim
!  SYNOPSIS
subroutine aims_elsi_read_mat_dim(rwh,filename,n_electron,n_basis,n_lrow,n_lcol)

  implicit none

  type(elsi_rw_handle), intent(inout) :: rwh
  character(*), intent(in) :: filename
  real(8), intent(out) :: n_electron
  integer, intent(out) :: n_basis
  integer, intent(out) :: n_lrow
  integer, intent(out) :: n_lcol

  call elsi_read_mat_dim(rwh,filename,n_electron,n_basis,n_lrow,n_lcol)

end subroutine
!******
!****s* FHI-aims/aims_elsi_read_mat_real
!  NAME
!   aims_elsi_read_mat_real
!  SYNOPSIS
subroutine aims_elsi_read_mat_real(rwh,filename,mat)

  implicit none

  type(elsi_rw_handle), intent(inout) :: rwh
  character(*), intent(in)  :: filename
  real(8), intent(out) :: mat(rwh%bh%n_lrow,rwh%bh%n_lcol)

  call elsi_read_mat_real(rwh,filename,mat)

end subroutine
!******
!****s* FHI-aims/aims_elsi_read_mat_cmplx
!  NAME
!   aims_elsi_read_mat_cmplx
!  SYNOPSIS
subroutine aims_elsi_read_mat_cmplx(rwh,filename,mat)

  implicit none

  type(elsi_rw_handle), intent(inout) :: rwh
  character(*), intent(in)  :: filename
  complex(8), intent(out) :: mat(rwh%bh%n_lrow,rwh%bh%n_lcol)

  call elsi_read_mat_complex(rwh,filename,mat)

end subroutine
!******
!****s* FHI-aims/aims_elsi_write_mat_real
!  NAME
!   aims_elsi_write_mat_real
!  SYNOPSIS
subroutine aims_elsi_write_mat_real(rwh,filename,mat)

  implicit none

  type(elsi_rw_handle), intent(in) :: rwh
  character(*), intent(in) :: filename
  real(8), intent(in) :: mat(rwh%bh%n_lrow,rwh%bh%n_lcol)

  call elsi_write_mat_real(rwh,filename,mat)

end subroutine
!******
!****s* FHI-aims/aims_elsi_write_mat_cmplx
!  NAME
!   aims_elsi_write_mat_cmplx
!  SYNOPSIS
subroutine aims_elsi_write_mat_cmplx(rwh,filename,mat)

  implicit none

  type(elsi_rw_handle), intent(in) :: rwh
  character(*), intent(in) :: filename
  complex(8), intent(in) :: mat(rwh%bh%n_lrow,rwh%bh%n_lcol)

  call elsi_write_mat_complex(rwh,filename,mat)

end subroutine
!******
end module
!******
