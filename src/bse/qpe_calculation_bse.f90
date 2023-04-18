!****s*  FHI-aims/qpe_calculation
!  NAME
!    qpe_calculation
!  SYNOPSIS

subroutine qpe_calculation_bse &
(qpe)

  !  PURPOSE
  !
  !  This subroutine performs self-energy calculation, base on GW approximation
  !  and/or MP2 approximation.  It also does the anlyatic continuation and eventully
  !  evaluates the quasipariticle energy using perturbation correction of ground
  !  state calculation (e.g. DFT, HF, or hybrid functional)
  !
  !  USES

  use dimensions
  use runtime_choices
  use species_data
  use physics
  use prodbas
  use hartree_fock
  use localized_basbas
  use gw_para
  use scgw_grid
  use poles_fit 
  use constants
  use mpi_tasks
  use timing
  use localorb_io
  use synchronize_mpi, only: sync_timing
  use evaluate_self_energy_freq,  only: evaluate_self_energy_freq_1,&
                                        evaluate_self_energy_freq_2
  use quasi_particle_energies,    only: quasi_particle_energy
  use print_self_energies,        only: print_self_energy_analytic

  !  INPUTS
  !    none
  !  OUTPUT
  !    none
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



  !  begin variables

  implicit none
  !  constants

  !  Input varibales
  !    partition_tab  :  partition function table
  !    n_electrons    :  number of electrons
  !    KS_eigenvalue  :
  !    KS_eigenvector :
  !    occ_numbers    :  occupation numbers of KS states

  ! Defined now in module "physics"
  !      real*8, dimension(n_max_angular, n_max_radial, n_atoms) ::
  !     +          partition_tab
  !      real*8  :: n_electrons
  !      real*8  :: chemical_potential
  !      real*8  :: KS_eigenvector(n_basis, n_states)
  !      real*8  :: KS_eigenvalue(n_states)
  !      real*8  :: occ_numbers(n_states)
  !      real*8, dimension(n_max_angular, n_max_radial, n_atoms) ::
  !     +         rho
  !      real*8, dimension(n_max_angular, n_max_radial, n_atoms) ::
  !     +         rho_gradient
  !  Output

  !  local arrays used in GW calculation :
  !    ovlp_3KS : the two single basis functions are transformed into KS states,
  !             this name is proper here anymore, to be changed in future
  !    xc_matr :  the LDA xc potential matrix between two single basis fns
  !    xc_KS_matr : the single basis of xc_matr transformed into KS states
  !    exchange_self_energy : (diagonal) matrix elements of the Fock operator
  !    self_energy_omega : the self energy on imaginary frequency axis
  real*8, dimension(:,:,:,:), allocatable :: ovlp_3KS
  real*8, dimension(:,:,:), allocatable :: xc_matr
  real*8, dimension(:,:,:), allocatable :: x_matr
  real*8, dimension(:,:,:), allocatable :: c_matr
  real*8, dimension(:,:,:), allocatable :: xc_KS_matr
  real*8, dimension(:,:), allocatable :: x_KS_array
  real*8, dimension(:,:), allocatable :: c_KS_array
  real*8, dimension(:,:), allocatable :: exchange_self_energy
  real*8, dimension(:,:), allocatable :: exchange_self_energy_short_range !needed for G0W0@HSE
  complex*16, dimension(:,:,:), allocatable :: self_energy_omega
  complex*16, dimension(:,:,:), allocatable :: mp2_selfenergy
  complex*16, dimension(:,:,:), allocatable :: sox_selfenergy
  complex*16, dimension(:,:,:), allocatable :: soxw_selfenergy
  complex*16, dimension(:,:,:), allocatable :: sosex_selfenergy
  real*8, dimension(:,:), allocatable :: qp_energy
  real*8, dimension(:,:), allocatable :: qp_energy_old

  real*8  eex_energy
  real*8  eex_energy_real
  real*8  xc_energy
  real*8, dimension(n_spin) ::  x_energy
  real*8  c_energy
  real*8, dimension(n_spin) :: x_energy_lda
  real*8  c_energy_lda
  real*8  c_energy_ldarpa
  real*8  rpa_energy
  real*8  rpa_tot_en
  real*8  rpa_plus_energy
  real*8  rpa_plus_tot_en
  real*8  e_2oex
  real*8  e_sosex
  real*8  en_single_excit
  real*8  en_ladder_se
  real*8  en_rse_full

  logical :: keep_o3fn

  integer ::  n_lumo(n_spin)

  !      real*8 time_3KS_intg
  real*8 time_xc_intg
  real*8 time_transform
  real*8 tot_time_sigma
  real*8 time_analy_continue
  real*8 time_qp_energy

  real*8 tot_time_polar
  real*8 n_bytes_o3KS

  !  counters
  integer i_state
  integer i_basis
  integer i_spin
  integer i_freq
  integer i_iter

  ! parameters for eigenvalue self-consistent GW
  logical scgw_converged 
!  logical use_ev_scgw 
  integer max_reiteration 
  real*8  threshold_sc


  !  error variables
  integer mpierr
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'qpe_calculation'

  !this are used in the computation of the GW Total energy
  complex*16, dimension(:,:,:), allocatable :: green_fn 
  complex*16, dimension(:,:), allocatable :: green_time_0 
  real*8 exchange_energy
  real*8 GW_corr
  real*8 tot_en_GW
  
  !   bsse time keeping
  real*8  temp_time_bsse
  real*8  temp_clock_time_bsse

    real*8 gap, gap0
!add cl
  real*8, dimension(n_states), intent(out) :: qpe
!  integer :: i, j
  integer :: i_task
  character(len=1024) :: filename
  integer :: i, j, k, l
!end cl 
! Igor, RPA potentials along AC path
  real*8  :: interval_ac_path      ! interval along AC path
  if (rpa_along_ac_path_grid .gt. 0) &
      interval_ac_path = real(rpa_along_ac_path_grid)**(-1)

  ! begin work
!add cl
  if(neutral_excitation_bse) then
    n_high_state = n_states
    print*, "n_states:", n_states
    print*, "n_high_state:", n_high_state
    print*, "n_low_state:", n_low_state
    print*, "n_basis_pairs", n_basis_pairs
    print*, "n_loc_prodbas", n_loc_prodbas
    print*, "n_basbas", n_basbas
  end if
!end cl

  if (myid.eq.0) then
     write(use_unit,'(A)')
     write(use_unit,'(A)')"--------------------------------------------"
     if(use_gw) then
        write(use_unit,'(10X,A)') "GW quasiparticle calculation starts ..."
     else if(use_mp2sf) then
        write(use_unit,'(10X,A)') "MP2 quasiparticle calculation starts ..."
     end if
     if(use_rpa_ene) then
        write(use_unit,'(10X,A)') "RPA total energy calculation starts ..."
     endif
     write(use_unit,"()")
  endif

  keep_o3fn = (use_scgw .or. use_scgw0 .or. use_mp2 .or. calculate_atom_bsse.or. use_dmft_gw &
               .or. use_dmft_pbe0)

  if(use_mpi) then
     call MPI_Bcast(KS_eigenvector, &
     n_basis*n_states*n_spin*n_k_points, &
     MPI_DOUBLE_PRECISION, &
     0, mpi_comm_global, mpierr)
  endif


  if(.not.allocated(n_homo)) then
     allocate(n_homo(n_spin))
     n_homo(:) = 0
     do i_spin = 1, n_spin, 1
       do i_state = 1, n_states, 1
         if(occ_numbers(i_state,i_spin,1).gt.1.e-6) then
            n_homo(i_spin) = i_state
         endif
       enddo
      enddo
  endif
 
  if(n_low_state .gt. int((n_electrons+1)/2)) then
     n_low_state = 1
  endif
  if (flag_frozen_core_postSCF) then ! count the frozen core states
      call count_frozen_core_states(n_low_state)
  endif
  if (n_high_state .le. int(n_electrons/2)+1) then
     n_high_state = int(n_electrons/2) + 10
  endif
  if (n_high_state.gt.n_states) then
     n_high_state=n_states
  elseif (n_high_state.lt.max(n_homo(1),n_homo(n_spin))) then
     n_high_state = max(n_homo(1),n_homo(n_spin))
  endif
 
  if(use_ev_scgw)then
    n_low_state = 1
    n_high_state=n_states 
  endif

  call get_timestamps(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
  if(use_2d_corr) then
     ! Attention: ndim1_o3KS/ndim2_o3KS must be the same everywhere
     ! and ovlp_3KS must not contains undefined (NaN) entries!
     ndim1_o3KS = (n_states-1)/np1_o3KS + 1
     ndim2_o3KS = (n_high_state-1)/np2_o3KS + 1
     ! do not allocate if it is there already as in the case of BSSE
     if (.not.allocated(ovlp_3KS)) then
         allocate(ovlp_3KS(n_basbas,ndim1_o3KS,ndim2_o3KS,n_spin), stat=info)
         call check_allocation(info, 'ovlp_3KS (2D)', func)
     endif
     ovlp_3KS = 0.
     n_bytes_o3KS = 8 * n_basbas * ndim1_o3KS * ndim2_o3KS * n_spin
  else
    if(use_gws2ox) then
      allocate(ovlp_3KS(n_loc_prodbas, n_states, n_states, n_spin),stat=info)
      call check_allocation(info, 'ovlp_3KS (1D)', func)
      n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_states * n_spin
    else
      allocate(ovlp_3KS(n_loc_prodbas, n_states, n_high_state, n_spin),stat=info)
      call check_allocation(info, 'ovlp_3KS (1D)', func)
      n_bytes_o3KS = 8 * n_loc_prodbas * n_states * n_high_state * n_spin
    endif
  end if

  write(info_str, "(2X,A,I7,A,I7,A,F12.3,A)") &
  & 'The ovlp_3KS matrix takes another', nint(n_bytes_o3KS / 2.d0**20), &
  & ' MiB x', n_tasks, ' procs =', n_tasks * n_bytes_o3KS / 2.d0**30, ' GiB.'
  call localorb_info(info_str, use_unit, '(A)', OL_norm)
  call localorb_info('', use_unit, '(A)', OL_norm)

! In the case of an LC-wPBEh calculation the ovlp_3fn array
! is the LR part of the ovlp_3fn complete array.
! This information isn't needed anymore.
! If the hybrid_coff .ne. 0 during that calculation, the SR part
! is available via ovlp_3fn_SR. If not we need to calculate that first.
! After that nearly everything is the same except the calculation of 
! the exchange_self_energy (see below)
	if (use_lc_wpbeh) then
		if (allocated(ovlp_3fn_SR)) then
			ovlp_3fn = 0.0d0
			ovlp_3fn = ovlp_3fn_SR
			deallocate(ovlp_3fn_SR)
		else
			ovlp_type_bare_or_hse_coul=OVLP_TYPE_HSE
			select case (RI_type)
              case(RI_SVS)
                 call get_coeff_3fn_svs(ovlp_3fn)
              case(RI_V)
                 if (use_2d_corr) then
                    call get_coeff_3fn_v_2d(ovlp_3fn)
                 else
                    call get_coeff_3fn_v_1d(ovlp_3fn)
                 end if
              case (RI_LVL)
                 call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
              case (RI_LVL_full)
                 call get_coeff_3fn_lvl_full(ovlp_3fn)
              case (RI_LVL_2nd)
                 call get_coeff_3fn_lvl(coeff_3fn_ten, coulomb_matr_lvl)
               ! Bare Coulomb integrals to ovlp_3fn
                 call integrate_ovlp3fn(l_shell_max,ext_l_shell_max, ovlp_3fn, OVLP_TYPE_COULOMB)
              case default
           call aims_stop("Invalid version of RI (resolution of identity)", func)
          end select        
		end if
	endif

!  if(.not.use_hartree_fock)then
!add cl
     if(myid == 0) print*, "use_2d_corr, sparse_o3fn, use_gws2ox", use_2d_corr, sparse_o3fn, use_gws2ox
!end cll
     if (use_2d_corr) then
        if (sparse_o3fn) then
           ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
           call ovlp3KS_lvl_2d(n_high_state, KS_eigenvector, &
           &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
        else
!add cl
     if(myid == 0) print*, 'transform_ovlp3fn_2 in qpe_calculation_bse'
!end cl
           call transform_ovlp3fn_2(n_high_state,KS_eigenvector,ovlp_3fn,ovlp_3KS)
        end if
     else
        if (sparse_o3fn) then
           ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
           if(use_gws2ox) then
             call ovlp3KS_lvl_1d(n_states, KS_eigenvector, &
                                coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
           else
             call ovlp3KS_lvl_1d(n_high_state, KS_eigenvector, &
                                coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
           endif
        else
           ! coeff_3fn -> ovlp_3KS
          if(use_gws2ox) then
            call transform_ovlp3fn(n_states, KS_eigenvector, ovlp_3fn, ovlp_3KS)
          else
            call transform_ovlp3fn(n_high_state, KS_eigenvector, ovlp_3fn, ovlp_3KS)
          endif
        end if
     end if    
!  endif
! add cl
  if (write_ovlp_3fn) then
      open(unit = 99, file = 'ovlp_3fn')
      do i = 1, n_basis_pairs
        do j = 1, n_loc_prodbas
          write(99, *) i, j, ovlp_3fn(i, j) 
        end do
      end do
      close(99)
  end if
! end cl   
!n_basis_pairs,n_loc_prodbas
    if (myid == 0) then
      open(unit = 99, file = 'ovlp_3fn_serial')
      do i = 1, n_basis_pairs
        do j = 1, n_loc_prodbas
          write(99, *) i, j, ovlp_3fn(i, j)
        enddo
      enddo
      close(99)
!n_loc_prodbas, n_states, n_states, n_spin
      open(unit = 99, file = 'ovlp_3ks_serial')
      do i = 1, n_loc_prodbas
        do j = 1, n_states
          do k = 1, n_states
            do l = 1, n_spin
              write(99, *) i, j, k, l, ovlp_3KS(i, j, k, l)
            enddo
          enddo
        enddo
      enddo
      close(99)
    endif
!    if (myid == 1) then
!      open(unit = 99, file = 'ovlp_3fn_proc1')
!      do i = 1, n_basis_pairs
!        do j = 1, n_loc_prodbas
!          write(99, *) i, j, ovlp_3fn(i, j)
!        enddo
!      enddo
!      close(99)
!
!      open(unit = 99, file = 'ovlp_3ks_proc1')
!      do i = 1, n_basbas
!        do j = 1, ndim1_o3KS
!          do k = 1, ndim2_o3KS
!            do l = 1, n_spin
!              write(99, *) i, j, k, l, ovlp_3KS(i, j, k, l)
!            enddo
!          enddo
!        enddo
!      enddo
!      close(99)  
!    endif
! end cl
   
  if (.not. keep_o3fn) then
     ! Array ovlp_3fn is not needed any more.
     if (sparse_o3fn) then
              deallocate(coulomb_matr_lvl)
              call dealloc_sp_ten(coeff_3fn_ten)
     else
              deallocate(ovlp_3fn)
     end if
  end if
  call get_times(time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)

  !   initialize the (imaginary) time and frequency grid
  call deallocate_gw() ! NEW: for atom_bsse
  call allocate_gw()
  !print * ,'freq_grid_type:' ,freq_grid_type
  if(freq_grid_type.eq.0) then
     call tf_ini(n_freq,n_full_freq, omegamax,omegamax, &
     omega_grid,omega_full_grid, &
     womega,womega_full,.true.)
  elseif(freq_grid_type.eq.1) then
     n_freq=n_full_freq
     call tf_ini_trans(n_freq,n_full_freq, omegamax,omegamax, &
     omega_grid,omega_full_grid, &
     womega,womega_full,.true.)
  elseif (freq_grid_type.eq.2) then
     if(myid.eq.0)write(use_unit,*)"LOGARITHMIC GRID"
!     n_freq=n_full_freq
     call loggrid_freq(0d0,omegamax, omega_grid(1:n_freq),womega(1:n_freq),n_freq)
     call loggrid_freq(0d0,omegamax, omega_full_grid(1:n_full_freq),womega_full(1:n_full_freq),n_full_freq)
  elseif (freq_grid_type.eq.3) then
     if(myid.eq.0)write(use_unit,*)"HOMOGENEOUS GRID"
!     n_freq=n_full_freq 
     call tf_ini_homo2(n_freq, n_full_freq, omegamax, &
     omega_grid,omega_full_grid, &
     womega,womega_full, .true.)     
  endif

  !    calculate the Fock energy Sigma_x
  
  if(use_gw .and. use_hse ) then
     use_gw_and_hse = .true.
!     if(myid.eq.0)write(use_unit,*) 'use_gw_and_hse   ' , use_gw_and_hse
  endif

  if(use_gw_and_hse)then

     if (.not.allocated(exchange_self_energy_short_range)) then
        allocate( exchange_self_energy_short_range (n_high_state,n_spin) )
     end if
   
     if(use_2d_corr) then
        call evaluate_exchange_energy_2 &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS, &
        exchange_self_energy_short_range, &
        eex_energy &
        )
     else
        call evaluate_exchange_energy &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS, &
        exchange_self_energy_short_range, &
        eex_energy &
        )
     endif
!     if(myid.eq.0)then
!       write(use_unit,*) "  Fock operator matrix elements: "
!       do i_state=1, n_states, 1 
!         write(use_unit,*)i_state, exchange_self_energy_short_range(i_state,1)*hartree 
!       enddo
!     endif

  else

     if (.not.allocated(exchange_self_energy)) then
        allocate( exchange_self_energy (n_high_state,n_spin) )
     end if

     if(use_2d_corr) then
        call evaluate_exchange_energy_2 &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS, &
        exchange_self_energy, &
        eex_energy &
        )
     else
        call evaluate_exchange_energy &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS(:,:,1:n_high_state,:), &
        exchange_self_energy, &
        eex_energy &
        )
     endif

  endif
  if(use_gw_and_hse) then
     if(myid.eq.0)then
       write(use_unit,*) "   Reinitializing product basis for G0W0@HSE.."
     endif
!     hse_omega_pbe = 0.0d0
!     hse_omega_hf  = 0.0d0
!     hse_omega  = 0.0d0
!     use_hse = .false.

     call prepare_corr_energy_calc()

     ovlp_3KS = 0.d0 
     if (use_2d_corr) then
        if (sparse_o3fn) then
           ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
           call ovlp3KS_lvl_2d(n_high_state, KS_eigenvector, &
           &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
        else
           call transform_ovlp3fn_2(n_high_state,KS_eigenvector,ovlp_3fn,ovlp_3KS)
        end if
     else
        if (sparse_o3fn) then
           ! coeff_3fn_ten, coulomb_matr_lvl -> ovlp_3KS
           call ovlp3KS_lvl_1d(n_high_state, KS_eigenvector, &
           &                   coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
        else
           ! coeff_3fn -> ovlp_3KS
           call transform_ovlp3fn(n_high_state, KS_eigenvector, ovlp_3fn, ovlp_3KS)
        end if
     end if 
!     if(myid.eq.0)then
!       write(use_unit,*) " Value of some variable called through the subroutine: " 
!       write(use_unit,*) " use_2d_corr  " , use_2d_corr 
!       write(use_unit,*) " sparse_o3fn  " , sparse_o3fn
!       write(use_unit,*) " keep_o3fn    " , keep_o3fn
!       write(use_unit,*) " use_rpa_ene  " , use_rpa_ene 
!       write(use_unit,*) " use_hse      " , use_hse
!       write(use_unit,*) " hybrid_coeff " , hybrid_coeff
!       write(use_unit,*) " use_corr     " , use_corr
!       write(use_unit,*) " use_hartree_fock " , use_hartree_fock
!       write(use_unit,*) " ovlp_type_bare_or_hse_coul ", OVLP_TYPE_NAMES(OVLP_TYPE_COULOMB)
!     endif

!     if(allocated (ovlp_3fn)) then
!       deallocate(ovlp_3fn)
!     endif 
!     call get_coeff_3fn_v_1d(ovlp_3fn)

     if (.not.allocated(exchange_self_energy)) then
        allocate( exchange_self_energy (n_high_state,n_spin) )
     end if
   
     if(use_2d_corr) then
        call evaluate_exchange_energy_2 &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS, &
        exchange_self_energy, &
        eex_energy &
        )
     else
        call evaluate_exchange_energy &
        (  n_high_state, n_homo, &
        occ_numbers, ovlp_3KS, &
        exchange_self_energy, &
        eex_energy &
        )
     endif
! for LC-wPBEh as starting DFT only the short-range exact energy is needed
	  if(use_lc_wpbeh) then
	  		exchange_self_energy (:,:) = exchange_self_energy_short_range (:,:) &
          * (1.d0 - hybrid_coeff)
	  else
     		exchange_self_energy (:,:) = exchange_self_energy (:,:)- & 
	       exchange_self_energy_short_range (:,:) * hybrid_coeff
	  end if

     if(allocated(exchange_self_energy_short_range))then
       deallocate(exchange_self_energy_short_range)
     endif

!     if(myid.eq.0)then
!       write(use_unit,*) " hybrid_coeff", hybrid_coeff  
!       write(use_unit,*) "  Fock operator matrix elements: "
!       do i_state=1, n_states, 1 
!         write(use_unit,*)i_state, exchange_self_energy(i_state,1)*hartree 
!       enddo
!     endif

  endif

  n_lumo(:) = n_homo(:) +1

  eex_energy_real = eex_energy
  if(use_hartree_fock) then
      eex_energy =  (1.d0-2.d0*hybrid_coeff)* eex_energy
  endif
  if(use_hybrid_meta) then
     en_hf = hybrid_coeff*eex_energy
  endif
  if(use_cohsex) then
     eex_energy = 2.d0*en_xc -eex_energy
  endif

  ! Evaluating C6 coefficients at MP2/RPA level
  if(use_C6_coef) then
     call get_timestamps(time_C6_coef, clock_time_C6_coef )
     if(use_2d_corr) then
        call evaluate_C6_coeff_2 &
        ( n_electrons, &
        i_start_mp2, n_high_state, &
        n_homo, n_lumo, &
        occ_numbers, n_full_freq, &
        omega_full_grid, womega_full, &
        chemical_potential, &
        KS_eigenvalue, KS_eigenvector, &
        ovlp_3KS)
     else
        call evaluate_C6_coeff &
        ( n_electrons, &
        i_start_mp2, n_high_state, &
        n_homo, n_lumo, &
        occ_numbers, n_full_freq, &
        omega_full_grid, womega_full, &
        chemical_potential, &
        KS_eigenvalue, KS_eigenvector, &
        ovlp_3KS)
     endif
     call get_times(time_C6_coef, clock_time_C6_coef)
  endif

  !      here we calculate the RPA correlation energy if required
  if(use_rpa_ene) then

     if (.not.allocated(xc_matr)) then
        allocate( xc_matr(n_basis,n_basis,n_spin), stat=info)
        call check_allocation(info, 'xc_matr', func)
     end if

     call integrate_xc_matrix ( &
          partition_tab, &
          rho, rho_gradient, &
          kinetic_density, &
          l_shell_max, &
          xc_matr)

     call get_timestamps(time_rse_corr, clock_time_rse_corr)
     if(use_2d_corr) then
        call evaluate_single_excitation_correction_2 &
        ( n_high_state,n_homo, &
        occ_numbers, ovlp_3KS, &
        KS_eigenvalue,KS_eigenvector, &
        xc_matr, en_single_excit, en_ladder_se &
        )
        call evaluate_renormalized_single_excitation_correction_2 &
        ( n_high_state,n_homo, &
        occ_numbers, ovlp_3KS, &
        KS_eigenvalue,KS_eigenvector, &
        xc_matr, en_rse_full &
        )
     else
        call evaluate_single_excitation_correction &
        ( n_high_state,n_homo, &
        occ_numbers, ovlp_3KS, &
        KS_eigenvalue,KS_eigenvector, &
        xc_matr, en_single_excit, en_ladder_se &
        )
        call evaluate_renormalized_single_excitation_correction &
        ( n_high_state,n_homo, &
        occ_numbers, ovlp_3KS, &
        KS_eigenvalue,KS_eigenvector, &
        xc_matr, en_rse_full &
        )
     endif
     call get_times(time_rse_corr, clock_time_rse_corr)


     rpa_energy = 0.d0
     if(use_2d_corr) then
        call evaluate_rpa_correlation_energy_2 &
        ( n_electrons, &
        n_low_state, n_high_state, &
        n_homo, n_lumo, &
        occ_numbers, n_full_freq, &
        omega_full_grid, womega_full, &
        chemical_potential, &
        KS_eigenvalue, KS_eigenvector, &
        ovlp_3KS, &
        rpa_energy &
        )
     else
        call evaluate_rpa_correlation_energy &
        ( n_electrons, &
        n_low_state, n_high_state, &
        n_homo, n_lumo, &
        occ_numbers, n_full_freq, &
        omega_full_grid, womega_full, &
        chemical_potential, &
        KS_eigenvalue, KS_eigenvector, &
        ovlp_3KS, &
        rpa_energy &
        )
     endif


     if(use_rpa_plus_2ox) then
        call get_timestamps(time_2oex_corr, clock_time_2oex_corr )
        if(use_2d_corr) then
           call evaluate_2oex_energy_2 &
           ( n_low_state,n_high_state,n_homo, &
           KS_eigenvalue,ovlp_3KS,e_2oex)
        else
           call evaluate_2oex_energy &
           ( n_low_state,n_high_state,n_homo, &
           KS_eigenvalue,ovlp_3KS,e_2oex)

        endif
        call get_times(time_2oex_corr, clock_time_2oex_corr )
     endif

     if(use_rpa_plus_sosex) then
        call get_timestamps(time_sosex_corr, clock_time_sosex_corr )
        if(use_2d_corr) then
          call evaluate_sosex_energy_2 &
           ( n_low_state, n_high_state, n_homo,&
             n_full_freq, omega_full_grid, womega_full,&
             occ_numbers, KS_eigenvalue, ovlp_3KS, e_sosex )
        else
          call evaluate_sosex_energy &
           ( n_low_state, n_high_state, n_homo,&
             n_full_freq, omega_full_grid, womega_full,&
             occ_numbers, KS_eigenvalue, ovlp_3KS, e_sosex )
        endif
        call get_times(time_sosex_corr, clock_time_sosex_corr )
     endif

     call integrate_xc_energy ( &
     partition_tab, &
     rho, rho_gradient, &
     kinetic_density, &
     xc_energy, &
     x_energy, &
     c_energy, &
     x_energy_lda, &
     c_energy_lda &
     )
     call integrate_crpa (&
     partition_tab, &
     rho, rho_gradient, &
     c_energy_ldarpa &
     )

     rpa_tot_en = total_energy-en_xc + &
     eex_energy + rpa_energy

     rpa_plus_energy = rpa_energy + &
     c_energy_lda -c_energy_ldarpa

     rpa_plus_tot_en = total_energy -en_xc + &
     eex_energy + rpa_plus_energy

     post_scf_total_energy = rpa_tot_en

     if(myid.eq.0) then
        write(use_unit,*) "--------------------------------------------------------------------"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exact exchange energy             : ", eex_energy_real, &
        " Ha,", eex_energy_real*hartree, " eV"
        if (n_spin.eq.1) then
           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
           "  DFT/HF exchange energy            : ",  x_energy_lda(1), &
           " Ha,", x_energy_lda(1)*hartree, " eV"
        else
           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
           "  DFT/HF exchange energy            : ",  x_energy_lda(1)+&
           x_energy_lda(2),&
           " Ha,", (x_energy_lda(1)+x_energy_lda(2))*hartree, " eV"
        endif
!        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!       "  DFT/HF correlation en.       : ",  c_energy_lda, &
!       " Ha,", c_energy_lda*hartree, " eV"
!       write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!       "  Local-RPA correlation energy : ",  c_energy_ldarpa, &
!       " Ha,", c_energy_ldarpa*hartree, " eV"
        write(use_unit,*)
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  DFT/HF total energy               : ",  total_energy, &
        " Ha,", total_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exchange-only total energy        : ",  total_energy-en_xc+eex_energy, &
        " Ha,", (total_energy-en_xc+eex_energy)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA total energy                  : ",  rpa_tot_en, &
        " Ha,", rpa_tot_en*hartree, " eV"
        write(use_unit,*)
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA+SE total energy               : ",  rpa_tot_en+en_single_excit, &
        " Ha,", (rpa_tot_en+en_single_excit)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA+rSE (diag) total energy       : ",  rpa_tot_en+en_ladder_se, &
        " Ha,", (rpa_tot_en+en_ladder_se)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA+rSE total energy              : ",  rpa_tot_en+en_rse_full, &
        " Ha,", (rpa_tot_en+en_rse_full)*hartree, " eV"
        write(use_unit,*)
!        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!       "  RPA+ corr. energy            : ",  rpa_plus_energy, &
!       " Ha,", rpa_plus_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA+ total energy                 : ",  rpa_plus_tot_en, &
        " Ha,", rpa_plus_tot_en*hartree, " eV"
        write(use_unit,*)
        if(use_rpa_plus_2ox) then
           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
           "  RPA+2OX total energy              : ",  rpa_tot_en + e_2oex, &
           " Ha,", (rpa_tot_en + e_2oex)*hartree, " eV"
!           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!           "  RPA+2OX+SE total energy     : ", &
!          rpa_tot_en + e_2oex + en_single_excit , &
!          " Ha,", (rpa_tot_en + e_2oex + en_single_excit)*hartree, " eV"
        endif
        if(use_rpa_plus_sosex) then
           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
           "  RPA+SOSEX total energy            : ",  rpa_tot_en + e_sosex, &
           " Ha,", (rpa_tot_en + e_sosex)*hartree, " eV"
           write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
           "  rPT2(=RPA+SOSEX+rSE) total energy : ",  rpa_tot_en + e_sosex+ en_rse_full, &
           " Ha,", (rpa_tot_en + e_sosex + en_rse_full)*hartree, " eV"

           post_scf_total_energy =  rpa_tot_en + e_sosex + en_rse_full

        endif
        if (rpa_along_ac_path_grid .gt. 0) then
           write(use_unit,'(2X,A)') &
               "RPA potentials along adiabatic connection path:"
           do i_state = 1, rpa_along_ac_path_grid
               write(use_unit,'(4X,A,f21.8,A,f21.8)') &
               "  Lambda, RPA potential : ", i_state*interval_ac_path, &
               ", ", rpa_potentials_along_ac_path(i_state)
           enddo
         endif

        !                if(use_screx) then
        !                 write(use_unit,*)
        !                 write(use_unit,'(2X,A,f21.8,A,f21.8,A)')
        !     +            "  COHSEX energy correction    : ", en_screx,
        !     +            " Ha,", en_screx*hartree, " eV"
        !                endif

        write(use_unit,*)"------------------------------------------", &
        "----------------------------------"
        write(use_unit,*)
     endif
     !   record the RPA total energy, RPA+SE total energy for atom_bsse
     if (calculate_atom_bsse) then
!  if first structure
         if (current_atom_num_bsse==0) then
             BSSE_full_energy_RPA_SE= (rpa_tot_en+en_single_excit)*hartree 
             BSSE_full_energy_RPA= (rpa_tot_en)*hartree 
             call get_timestamps(temp_time_bsse, temp_clock_time_bsse )
             time_bsse =temp_time_bsse-time_bsse
             clock_time_bsse =  temp_clock_time_bsse-clock_time_bsse 
             if(myid.eq.0) then
                write(use_unit,'(2X,A,F12.3,A)') &
                "Total time taken for full structure              : ", &
                time_bsse, " s"
             endif
         else       
!  if second structure and so on
             BSSE_atom_RPA_SE(current_atom_num_bsse)=(rpa_tot_en+en_single_excit)*hartree
             BSSE_atom_RPA(current_atom_num_bsse)= (rpa_tot_en)*hartree 
             call get_timestamps(temp_time_bsse, temp_clock_time_bsse )
             time_bsse =temp_time_bsse-time_bsse
             clock_time_bsse =  temp_clock_time_bsse-clock_time_bsse 
             if(myid.eq.0) then
                write(use_unit,'(2X,A,F12.3,A)') &
                "Total time taken for the last BSSE step              : ", &
                time_bsse, " s"
             endif
         endif
     endif

     ! end of if use_rpa_ene
  endif

  if(use_gw.or.use_mp2sf) then

     !     Calculate the exchange-correlation energy matrix within local basis
     call cpu_time(time_xc_intg)
     if (.not.allocated(xc_matr)) then
        allocate( xc_matr(n_basis,n_basis,n_spin) )
     end if

     !  if(myid.eq.0) write(use_unit,*) ' HERE1 '
     if(use_screx) then
        xc_matr(:,:,:) = -fock_matr(:,:,:) + screx_matr(:,:,:)
     else
       if(use_split_xc_gw) then
         if (.not.allocated(x_matr)) then
           allocate( x_matr(n_basis,n_basis,n_spin), stat=info)
           call check_allocation(info, 'x_matr', func)
         end if

         if (.not.allocated(c_matr)) then
           allocate( c_matr(n_basis,n_basis,n_spin), stat=info)
           call check_allocation(info, 'c_matr', func)
         end if

         call integrate_split_xc_matrix ( &
                partition_tab, &
                rho, rho_gradient, &
                kinetic_density, &
                l_shell_max, &
                x_matr,c_matr)
        else

         call integrate_xc_matrix ( &
                partition_tab, &
                rho, rho_gradient, &
                kinetic_density, &
                l_shell_max, &
                xc_matr)
! end of if use_split_xc_gw
      endif
! end of if use_screw
     endif


     if(use_split_xc_gw) then
       if (.not.allocated(x_KS_array)) then
          allocate( x_KS_array(n_states,n_spin) )
       end if
       if (.not.allocated(c_KS_array)) then
          allocate( c_KS_array(n_states,n_spin) )
       end if
       call evaluate_KS_split_xc_matrix &
       (KS_eigenvector, &
        x_matr,c_matr,x_KS_array,c_KS_array)
       ! dummy allocation - in this case, this array should never be used
       if (.not.allocated(xc_KS_matr)) then
          allocate(xc_KS_matr(1,1,1))
       end if
     else 
       if (.not.allocated(xc_KS_matr)) then
          allocate(xc_KS_matr(n_states,n_states,n_spin))
       end if
      call evaluate_KS_xc_matrix &
       (KS_eigenvector, &
       xc_matr, xc_KS_matr)
       ! dummy allocations - in this case, these arrays should never be used
       if (.not.allocated(x_KS_array)) then
          allocate( x_KS_array(1,1) )
       end if
       if (.not.allocated(c_KS_array)) then
          allocate( c_KS_array(1,1) )
       end if       
     endif

!     do i_state = 1, n_states, 1
!        if(myid.eq.0)then
!          write(use_unit,*), i_state, xc_KS_matr(i_state, i_state,1)*hartree
!        endif
!     enddo

     call cpu_time(rtime)
     time_xc_intg = rtime - time_xc_intg


     !  calculate the correlation part of the self energy Sigma_c
     
     if (.not.allocated(self_energy_omega)) then

        allocate (self_energy_omega &
        ( n_low_state:n_high_state,n_freq,n_spin))

        self_energy_omega(:,:,:) = (0.d0,0.d0)
     end if


!     do i_iter = 1, n_iteration, 1 
     if(use_gw .and. .not. use_ev_scgw) then

        call get_timestamps(time_self_energy,clock_time_self_energy )
        if(use_2d_corr) then
           call evaluate_self_energy_freq_2 &
           (n_low_state, n_high_state, &
           n_homo, occ_numbers, &
           n_freq,n_full_freq, &
           omega_grid, omega_full_grid, &
           womega_full, &
           chemical_potential_spin, &
           KS_eigenvalue,KS_eigenvector, &
           partition_tab, l_shell_max, &
           ovlp_3KS, &
           self_energy_omega, out_self_energy, &
           flag_coh, &
           tot_time_sigma, tot_time_polar &
           )

           !                call evaluate_single_excitation_to_self_energy_2 &
           !                 ( n_low_state,n_high_state, &
           !                   n_homo,n_freq, omega_grid, &
           !                   occ_numbers, ovlp_3KS, &
           !                   KS_eigenvalue,KS_eigenvector, &
           !                   xc_matr, self_energy_omega &
           !                  )

        else
           call evaluate_self_energy_freq_1 &
           (n_low_state, n_high_state, &
           n_homo, &
           occ_numbers, &
           n_freq,n_full_freq, &
           omega_grid,womega, omega_full_grid, &
           womega_full, &
           chemical_potential_spin, &
           KS_eigenvalue,KS_eigenvector, &
           partition_tab, l_shell_max, &
           ovlp_3KS(:,:,1:n_high_state,:), &
           self_energy_omega, &
           flag_coh, &
           tot_time_sigma, tot_time_polar &
           )


         endif

           call get_times(time_self_energy, clock_time_self_energy)
! adding the 2nd-order exchange self-energy to G0W0 self-energy
          if(use_gw2ox) then
            call get_timestamps(time_2ox_selfenergy,clock_time_2ox_selfenergy )
            allocate (sox_selfenergy(n_freq,n_high_state,n_spin))
            sox_selfenergy = (0.d0,0.d0)

            if(use_2d_corr) then
              call evaluate_sox_selfenergy_2 &
              ( n_freq, n_low_state, n_high_state, &
                omega_grid, ovlp_3KS, &
                chemical_potential_spin, &
                n_electrons,occ_numbers, &
                KS_eigenvalue, &
                sox_selfenergy )
            else
              call evaluate_sox_selfenergy &
              ( n_freq, n_low_state, n_high_state, &
                omega_grid, ovlp_3KS, &
                chemical_potential_spin, &
                n_electrons,occ_numbers, &
                KS_eigenvalue, & 
                sox_selfenergy )
            endif

            do i_spin = 1, n_spin
              do i_freq = 1, n_freq
                do i_state = n_low_state, n_high_state
                  self_energy_omega(i_state, i_freq, i_spin) =  &
                  self_energy_omega(i_state, i_freq, i_spin) + &
                  sox_selfenergy(i_freq,i_state,i_spin)!/2/pi
                enddo
              enddo
            enddo
            call get_times(time_2ox_selfenergy, clock_time_2ox_selfenergy)

! adding the 2nd-order exchange self-energy in W to G0W0 self-energy
          elseif(use_gwsoxw) then
            call get_timestamps(time_soxw_selfenergy,clock_time_soxw_selfenergy )
            allocate (soxw_selfenergy(n_freq,n_high_state,n_spin))
            soxw_selfenergy = (0.d0,0.d0)

            if(use_2d_corr) then
                write(use_unit,*) "  The 2D parallel implementation of GW+SOXW is not working yet ... "
                stop
            else
              call evaluate_soxw_selfenergy &
              ( n_freq, n_low_state, n_high_state, &
                omega_grid, ovlp_3kS, &
                chemical_potential_spin, &
                n_electrons,occ_numbers, &
                KS_eigenvalue, &
                soxw_selfenergy &
              )

            endif

            do i_spin = 1, n_spin
              do i_freq = 1, n_freq
                do i_state = n_low_state, n_high_state
                  self_energy_omega(i_state, i_freq, i_spin) =  &
                  self_energy_omega(i_state, i_freq, i_spin) + &
                  soxw_selfenergy(i_freq,i_state,i_spin)
                enddo
              enddo
            enddo
            call get_times(time_soxw_selfenergy, clock_time_soxw_selfenergy)

        elseif(use_gws2ox) then
! adding the 2nd-order screened exchange self-energy to G0W0 self-energy
            call get_timestamps(time_s2ox_selfenergy,clock_time_s2ox_selfenergy )
           allocate (sosex_selfenergy(n_freq,n_high_state,n_spin))

           if(use_2d_corr) then
              if(myid.eq.0) then
                 write(use_unit,*) "2D GW+2SOX implementation is not available, please use the 1D version !"
              endif
              stop
           else
             if(use_sosex_selfe) then
                call evaluate_sosex_selfenergy &
                ( n_low_state, n_high_state, &
                  ovlp_3KS, &
                  n_freq,n_full_freq, &
                  omega_grid, &
                  omega_full_grid, womega_full, &
                  chemical_potential_spin, &
                  n_electrons,occ_numbers, &
                  KS_eigenvalue, &
                  sosex_selfenergy &
                 )
              elseif(use_sosex2w_selfe) then
                 call evaluate_sosex_2w_selfenergy &
                 ( n_low_state, n_high_state, &
                   ovlp_3KS, &
                   n_freq,n_full_freq, &
                   omega_grid, &
                   omega_full_grid, womega_full, &
                   chemical_potential_spin, &
                   n_electrons,occ_numbers, &
                   KS_eigenvalue, &
                  sosex_selfenergy &
                 )
               endif
           endif

           do i_spin = 1, n_spin
              do i_freq = 1, n_freq
                do i_state = n_low_state, n_high_state
                  self_energy_omega(i_state, i_freq, i_spin) =  &
                  self_energy_omega(i_state, i_freq, i_spin) + &
                  sosex_selfenergy(i_freq,i_state,i_spin)
                enddo
              enddo
            enddo
           call get_times(time_s2ox_selfenergy, clock_time_s2ox_selfenergy)
!  end of if use gw and beyond
        endif
        if (out_self_energy) then
            call output_self_energy(n_freq,n_low_state,n_high_state,omega_grid,self_energy_omega)
        endif

     elseif (use_mp2sf) then

        allocate (mp2_selfenergy &
        (n_freq,n_high_state,n_spin))
        mp2_selfenergy = (0.d0,0.d0)

        if(use_2d_corr) then
           call evaluate_mp2_selfenergy_2 &
           ( n_freq, n_low_state, n_high_state, &
           omega_grid, ovlp_3KS, &
           chemical_potential_spin, &
           n_electrons,occ_numbers, &
           KS_eigenvalue, &
           mp2_selfenergy &
           )
        else
           call evaluate_mp2_selfenergy &
           ( n_freq, n_low_state, n_high_state, &
           omega_grid, ovlp_3KS, &
           chemical_potential_spin, &
           n_electrons,occ_numbers, &
           KS_eigenvalue, &
           mp2_selfenergy &
           )
        endif

        do i_spin = 1, n_spin
           do i_freq = 1, n_freq
              do i_state = n_low_state, n_high_state
                 self_energy_omega(i_state, i_freq, i_spin) = &
                 mp2_selfenergy(i_freq,i_state,i_spin)
              enddo
           enddo
        enddo
     endif ! if(use_gw) / elseif (use_mp2sf) 


!-------------------------------------------------------------------
!-------------------------------------------------------------------
!-------------------------------------------------------------------
    if (use_ev_scgw)then !this loop is kept separated from G0W0 for safety reasons, for now

       !Loop initialization for eigenvalue self-consistent GW, 
       !at each iteration the self-energy is recalculated using the previous quasiparticle energies
       !until convergence of the HOMO level

       if(.not. allocated(qp_energy_old)) then
          allocate(qp_energy_old(1:n_states,n_spin))
          qp_energy_old(:,:) = 0.d0
       endif

       qp_energy_old(1:n_states,:) = KS_eigenvalue (1:n_states,:,1)
       i_iter = 0
       scgw_converged = .false.
       threshold_sc = 1.d-4
       n_iteration = 30
       max_reiteration = n_iteration  

       if (myid.eq.0 .and. max_reiteration .gt. 1) then
          write(use_unit,*) "  Self-Consistent loop initialization ... "
          write(use_unit,*) " "
       endif

       do while(.not. scgw_converged)

         i_iter = i_iter +1

         if (myid.eq.0 .and. max_reiteration .gt. 1) then
            write(use_unit,*) " ----------------------------------------------"
            write(use_unit,*) "             Iteration ", i_iter 
            write(use_unit,*) " ----------------------------------------------"
            write(use_unit,*) " "
         endif

         call get_timestamps(time_self_energy,clock_time_self_energy )
         if(use_2d_corr) then
           call evaluate_self_energy_freq_2 &
           ( n_low_state, n_high_state, &
           n_homo, occ_numbers, &
           n_freq,n_full_freq, &
           omega_grid, omega_full_grid, &
           womega_full, &
           chemical_potential_spin, &
           qp_energy_old,KS_eigenvector, &
           partition_tab, l_shell_max, &
           ovlp_3KS, &
           self_energy_omega, out_self_energy, &
           flag_coh, &
           tot_time_sigma, tot_time_polar &
           )

           !                call evaluate_single_excitation_to_self_energy_2 &
           !                 ( n_low_state,n_high_state, &
           !                   n_homo,n_freq, omega_grid, &
           !                   occ_numbers, ovlp_3KS, &
           !                   KS_eigenvalue,KS_eigenvector, &
           !                   xc_matr, self_energy_omega &
           !                  )

         else
           call evaluate_self_energy_freq_1 &
           (n_low_state, n_high_state, &
           n_homo, &
           occ_numbers, &
           n_freq,n_full_freq, &
           omega_grid,womega, omega_full_grid, &
           womega_full, &
           chemical_potential_spin, &
           qp_energy_old,KS_eigenvector, &
           partition_tab, l_shell_max, &
           ovlp_3KS(:,:,1:n_high_state,:), &
           self_energy_omega, &
           flag_coh, &
           tot_time_sigma, tot_time_polar &
           )


         endif

         call get_times(time_self_energy, clock_time_self_energy)

         if(.not. allocated(qp_energy)) then
            allocate(qp_energy(n_low_state:n_high_state,n_spin))
            qp_energy(:,:) = 0.d0
         endif
         if (.not.allocated(sigma_par)) then
            allocate(sigma_par(n_max_par,n_states,n_spin))
            sigma_par(:,:,:) = 0.d0
         end if
       
         call analy_continue_self_energy &
            (flag_qpe, anacon_type, &
            n_low_state,n_high_state,n_freq, &
            n_max_par, &
            qp_energy_old,chemical_potential_spin, &
            sigma_par,omega_grid, &
            self_energy_omega)
       
         if(qpe_multi_solu) then
            call quasi_particle_energy_multisolu ( anacon_type, &
              n_max_par, n_low_state, n_high_state, &
              n_homo, n_lumo, n_freq, omega_grid, &
              sigma_par, occ_numbers, KS_eigenvalue,chemical_potential_spin, &
              exchange_self_energy, &
              xc_KS_matr,x_KS_array,c_KS_array, &
              qp_energy  )
         else
            call quasi_particle_energy ( anacon_type, &
              n_max_par, n_low_state, n_high_state, &
              n_homo, n_freq, omega_grid, &
              sigma_par, occ_numbers, KS_eigenvalue,chemical_potential_spin, &
              exchange_self_energy, &
              xc_KS_matr,x_KS_array,c_KS_array, &
              qp_energy  )
          endif


         !Convergence check for eigenvalues self-consistent GW
         if (n_spin ==2)then
           if ((abs(qp_energy(n_homo(1),1)-qp_energy_old(n_homo(1),1) ).lt. threshold_sc).and.&
               (abs(qp_energy(n_homo(2),2)-qp_energy_old(n_homo(2),2) ).lt. threshold_sc))then
             scgw_converged = .true.
           endif
         elseif(n_spin==1)then
           if (abs(qp_energy(n_homo(1),1)-qp_energy_old(n_homo(1),1) ) .lt. threshold_sc)then
             scgw_converged = .true.
           endif
         endif
        
         qp_energy_old (n_low_state:n_high_state,:) = qp_energy (n_low_state:n_high_state,:)
        
         if( allocated(qp_energy)) then
            deallocate(qp_energy)
         endif
         if (allocated(sigma_par)) then
            deallocate(sigma_par)
         end if

         if ( i_iter.ge.max_reiteration ) then 
            scgw_converged = .true.
            if (myid.eq.0 .and. max_reiteration .gt.1 )then
              write(use_unit,*) "     Eigenvalue self-consistent GW did not converged to the required accuracy! Tschuess!"
            endif
         endif
       enddo !while (.not. scgw_converged )

    if( allocated(qp_energy_old)) then
       deallocate(qp_energy_old)
    endif

    endif 
 
! GW Total Energy from Galitskij-Migdal formula
    if(use_gw_energy)then
        
       call get_gw_tot_en (eex_energy, self_energy_omega)
 !      if(.not. allocated(green_fn)) then
 !        allocate(green_fn (n_states,n_freq, n_spin))
 !      endif
 !      !update the Green's function from dyson's equation
 !      exchange_energy = 0.d0
 !      GW_corr = 0.d0
 !      green_time_0 (:,:) = 0.d0
 !      do i_spin = 1, n_spin, 1
 !        do i_state =1 , n_states
 !          do i_freq = 1, n_freq, 1
 !              green_fn (i_state, i_freq, i_spin) = &
 !                1.d0/((0.d0,1.d0)*omega_grid(i_freq)-&
 !              (KS_eigenvalue(i_state, i_spin,1)-chemical_potential))
 ! 
 !          GW_corr = GW_corr + &
 !           green_fn (i_state, i_freq, i_spin)*&
 !           self_energy_omega(i_state, i_freq, i_spin)*&
 !           womega(i_freq)/pi/2.d0/n_spin
 !          enddo
 !        enddo
 !       enddo
 ! 
 !       tot_en_GW = total_energy - en_xc + &
 !        eex_energy + GW_corr
 ! 
 !      if(myid.eq.0) print *, "exchange_energy1 :", eex_energy, 'Ha', eex_energy*hartree, 'eV' 
 !      if(myid.eq.0) print *, "GW_corr          :", GW_corr, 'Ha', GW_corr*hartree, 'eV'
 !      if(myid.eq.0) print *, "GW_Tot           :", tot_en_GW, 'Ha', tot_en_GW*hartree, 'eV'
 ! 
 !      if( allocated(green_fn)) then
 !        deallocate(green_fn)
 !      endif

    endif !use_gw_energy
 

     !   quasi particle calculation
     if(myid.eq.0) then

        write(use_unit,'(A)') " --------------------------------------------"
        if (.not.allocated(sigma_par)) then
           allocate(sigma_par(n_max_par,n_states,n_spin))
        end if
        call cpu_time (time_analy_continue)
        call analy_continue_self_energy &
        (flag_qpe, anacon_type, &
        n_low_state,n_high_state,n_freq, &
        n_max_par, &
        KS_eigenvalue,chemical_potential_spin, &
        sigma_par,omega_grid, &
        self_energy_omega)
        call cpu_time (rtime)
        time_analy_continue = rtime - time_analy_continue

        call cpu_time (time_qp_energy)
        if(.not. allocated(qp_energy)) then
           allocate(qp_energy(n_low_state:n_high_state,n_spin))
           qp_energy(:,:) = 0.d0
        endif
        if(qpe_multi_solu) then
           call quasi_particle_energy_multisolu ( anacon_type, &
            n_max_par, n_low_state, n_high_state, &
            n_homo, n_lumo, n_freq, omega_grid, &
            sigma_par, occ_numbers, KS_eigenvalue,chemical_potential_spin, &
            exchange_self_energy, &
            xc_KS_matr,x_KS_array,c_KS_array, &
            qp_energy  )
        else
           call quasi_particle_energy ( anacon_type, &
            n_max_par, n_low_state, n_high_state, &
            n_homo, n_freq, omega_grid, &
            sigma_par, occ_numbers, KS_eigenvalue,chemical_potential_spin, &
            exchange_self_energy, &
            xc_KS_matr,x_KS_array,c_KS_array, &
            qp_energy  )
        endif
        call cpu_time (rtime)
        time_qp_energy = rtime - time_qp_energy
 
!        do i_spin =1 , n_spin
!           do i_state = 1, n_states, 1 
!             if(myid.eq.0) print *, i_state, i_spin, qp_energy (i_state,i_spin )
!           enddo 
!        enddo 


        if(.false.)then
           call get_qp_spectrum            &
           ( anacon_type, &
           n_max_par, n_low_state, n_high_state, &
           n_homo, n_lumo, n_freq, omega_grid, &
           sigma_par, KS_eigenvalue, &
           exchange_self_energy, &
           xc_KS_matr,"QP_spectrum.dat")
        endif
        if(flag_print_self_energy)then
           call print_self_energy_analytic &
           ( anacon_type, &
           n_max_par, n_freq, omega_grid, &
           sigma_par, state_to_print, print_range)
        endif

     endif!  if(myid.eq.0) then

     if(myid.eq.0.and.out_dos) then
        call get_species_projected_dos('GW',n_low_state, n_high_state-n_low_state+1, &
        qp_energy(n_low_state:n_high_state,1:n_spin))
        call get_species_projected_dos('KS',1,n_states,KS_eigenvalue)
     endif

     !     end of if use_gw or use_mp2sf
  endif

  if(.not.(use_mp2 .or. calculate_atom_bsse) .and. .not. use_scgw0 .and. .not. use_scgw) then
     call cleanup_hartree_fock()
     if (use_lvl) call cleanup_localized_basbas()
  endif

  !  Final tasks

  !       deallocate everything
 if (.not. calculate_atom_bsse .and. .not. use_scgw0 .and. .not. use_scgw) then
  if (allocated(ovlp_3KS)) then
     deallocate( ovlp_3KS     )
  end if
  if (allocated(self_energy_omega)) then
     deallocate( self_energy_omega )
  endif
  if (allocated(sigma_par)) then
     deallocate( sigma_par )
  endif
  if(allocated(xc_matr)) then
     deallocate(xc_matr)
  endif
  if(allocated(x_matr)) then
     deallocate(x_matr)
  endif
  if(allocated(c_matr)) then
     deallocate(c_matr)
  endif
  if(allocated(xc_KS_matr)) then
     deallocate(xc_KS_matr)
  endif
  if(allocated(x_KS_array)) then
     deallocate(x_KS_array)
  endif
  if(allocated(c_KS_array)) then
     deallocate(c_KS_array)
  endif
  if(allocated(n_homo)) then
     deallocate(n_homo)
  endif
  if (allocated(mp2_selfenergy)) then
     deallocate(mp2_selfenergy)
  endif
  if (allocated(sox_selfenergy)) then
     deallocate(sox_selfenergy)
  endif
  if (allocated(soxw_selfenergy)) then
     deallocate(soxw_selfenergy)
  endif
  if (allocated(sosex_selfenergy)) then
     deallocate(sosex_selfenergy)
  endif
 endif
!add g
  if (neutral_excitation_bse) then
    qpe = qp_energy(:, n_spin)
  end if
!end g
  !     final time accounting and exit

  if(myid.eq.0 .and. use_gw) then
     write(use_unit,'(10X,A,F12.3,A)') &
     "| Total time for transforming the 3-center integrals : ", &
     time_ovlp3fn_to_ovlp3KS, " s"
     write(use_unit,'(10X,A,/,10X,A,F12.3,A)') &
     "| Total time for calculating the exchange-correlation", &
     "                        energy matrix elements   : ", &
     time_xc_intg, " s"
     write(use_unit,'(10X,A,F12.3,A)') &
     "| Total time for calculating polarisability      : ", &
     tot_time_polar, " s"
     write(use_unit,'(10X,A,/,10X,A,F12.3,A)') &
     "| Total time for calculating self energy ", &
     "              on imaginary frequency axis        : ", &
     tot_time_sigma, " s"
     write(use_unit,'(10X,A,/,10X,A,F12.3,A)') &
     "| Total time for calculating the quasiparticle ", &
     "              energies                           : ", &
     time_qp_energy, " s"
  endif

end subroutine qpe_calculation_bse
!******
