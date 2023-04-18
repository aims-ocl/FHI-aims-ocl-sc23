!****s* FHI-aims/check_auxmat
!  NAME
!    check_auxmat
!  SYNOPSIS

program check_auxmat

  !  PURPOSE
  !
  !  USES

  use prodbas
  use sbt_overlap_aims
  use tight_binding_auxmat
  use debug_output
  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
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

  real*8, allocatable :: C1(:,:), C2(:,:), C3(:,:)
  integer :: ovlp_type
  integer :: info
  character*150 :: info_str
  integer :: i
  character(*), parameter :: func = 'check_auxmat'

  call init_aims()
  if (.not. use_prodbas) call aims_stop('No product basis from control.in')

  allocate(C1(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C1', func)
  allocate(C2(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C2', func)
  allocate(C3(n_basbas, n_basbas), stat=info)
  call check_allocation(info, 'C3', func)

  do i = 1, 4
     select case (i)
     case(1)
        ovlp_type = OVLP_TYPE_OVERLAP
     case(2)
        ovlp_type = OVLP_TYPE_COULOMB
     case(3)
        use_cutCb = .true.
        cutCb_width = exp(cutCb_width_factor * sbtgrid_lnrange/sbtgrid_N)
        cutCb_rcut = cutCb_rcut_factor * 8.d0/bohr
        ovlp_type = OVLP_TYPE_CUT
     case(4)
        use_hse = .true.
        hse_omega_hf = 0.11d0
        ovlp_type = OVLP_TYPE_HSE
        call cleanup_basbas()
        call initialize_prodbas()
     case default
        cycle
     end select

     write(0,*) OVLP_TYPE_NAMES(ovlp_type)
     write(0,*) 'standard'
     call integrate_auxmat_by_atomic_sbt(C1, ovlp_type, .false., use_tb=.false.)
     write(0,*) 'slow TB'
     call integrate_auxmat_by_atomic_sbt(C2, ovlp_type, .false., use_tb=.true.)
     write(0,*) 'fast TB'
     call integrate_auxmat_by_tb(C3, ovlp_type)

     write(0,"('max|standard - slow_tb|',ES10.2)") maxval(abs(C1-C2))
     write(0,"('max|standard - fast_tb|',ES10.2)") maxval(abs(C1-C3))
     write(0,"('max|slow_tb  - fast_tb|',ES10.2)") maxval(abs(C2-C3))
     write(0,*)
!!$     if (ovlp_type == OVLP_TYPE_CUT) then
!!$        call debug_array(22, C1, tofile='C1.dat')
!!$        call debug_array(22, C2, tofile='C2.dat')
!!$        call debug_array(22, C3, tofile='C3.dat')
!!$     end if
  end do


end program check_auxmat
!******
