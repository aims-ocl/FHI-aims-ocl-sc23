!****s*  FHI-aims/qpe_ospt2
!  NAME
!    qpe_ospt2
!  SYNOPSIS

subroutine qpe_ospt2 &
( )

  !  PURPOSE
  !
  !  This subroutine performs os-PT2 calculations in  the RPA-type ACDF framework.
  !
  !  USES

  use dimensions
  use runtime_choices
  use species_data
  use physics
  use prodbas
  use crpa_blacs
  use hartree_fock
  use localized_basbas
  use gw_para
  use scgw_grid
  use constants
  use mpi_tasks
  use timing
  use localorb_io
  use synchronize_mpi, only: sync_timing

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
  real*8, dimension(:,:), allocatable :: correl_energy

  real*8  eex_energy
  real*8  eex_energy_real
  real*8  xc_energy
  real*8, dimension(n_spin) ::  x_energy
  real*8  c_energy
  real*8, dimension(n_spin) :: x_energy_lda
  real*8  c_energy_lda
  real*8  c_energy_ldarpa
  real*8  rpa_energy
  real*8  rpa_energy_os
  real*8  rpa_energy_ss
  real*8  rpa_energy_tot
  real*8  rpa_energy_scs
  real*8  rpa_tot_en
  real*8  rpa_tot_en_tot
  real*8  rpa_tot_en_scs
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
  character(*), parameter :: func = 'qpe_ospt2'

  !this are used in the computation of the GW Total energy
  complex*16, dimension(:,:,:), allocatable :: green_fn 
  complex*16, dimension(:,:), allocatable :: green_time_0 
  real*8 exchange_energy
  real*8 GW_corr
  real*8 tot_en_GW

  !required for contour deformation 
  real*8, dimension(:,:,:,:), allocatable :: Wmn_freq_cd
  

    real*8 gap, gap0
  
! Igor, RPA potentials along AC path
  real*8  :: interval_ac_path      ! interval along AC path
  if (rpa_along_ac_path_grid .gt. 0) &
      interval_ac_path = real(rpa_along_ac_path_grid)**(-1)

  ! begin work


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

  if (flag_frozen_core_postSCF) then ! count the frozen core states
      call count_frozen_core_states(n_low_state)
  endif
  if(n_low_state .gt. int((n_electrons+1)/2)) then
     n_low_state = 1
  endif
  if (n_high_state .le. int(n_electrons/2)+1) then
     n_high_state = int(n_electrons/2) + 10
  endif
  if (n_high_state.gt.n_states) then
     n_high_state=n_states
  elseif (n_high_state.lt.max(n_homo(1),n_homo(n_spin))) then
     n_high_state = max(n_homo(1),n_homo(n_spin))
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

!  if(.not.use_hartree_fock)then
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
                               coeff_3fn_ten, coulomb_matr_lvl, ovlp_3KS)
        else
           ! coeff_3fn -> ovlp_3KS
           call transform_ovlp3fn(n_high_state, KS_eigenvector, ovlp_3fn, ovlp_3KS)
        end if
     end if    
!  endif
   
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

  !if(use_2d_corr) call init_crpa_blacs_distribution(n_full_freq)

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

  !      here we calculate the RPA correlation energy if required

  if(use_os_mp2_qpe) then

     rpa_energy = 0.d0
     !if(use_2d_corr) then
     !   call evaluate_rpa_correlation_energy_2 &
     !   ( n_electrons, &
     !   n_low_state, n_high_state, &
     !   n_homo, n_lumo, &
     !   occ_numbers, n_full_freq, &
     !   omega_full_grid, womega_full, &
     !   chemical_potential, &
     !   KS_eigenvalue, KS_eigenvector, &
     !   ovlp_3KS, &
     !   rpa_energy &
     !   )
     !else
        call evaluate_qpe_ospt2_energy &
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
     !endif


     rpa_tot_en = total_energy-en_xc + &
     eex_energy + rpa_energy

     post_scf_total_energy = rpa_tot_en

     if(myid.eq.0) then
        write(use_unit,*) "--------------------------------------------------------------------"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exact exchange energy             : ",  eex_energy_real, &
        " Ha,", eex_energy_real*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osPT2 correlation energy          : ",  rpa_energy, &
        " Ha,", rpa_energy*hartree, " eV"
        write(use_unit,*)
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  DFT/HF total energy               : ",  total_energy, &
        " Ha,", total_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exchange-only total energy        : ",  total_energy-en_xc+eex_energy, &
        " Ha,", (total_energy-en_xc+eex_energy)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osPT2 total energy                : ",  rpa_tot_en, &
        " Ha,", rpa_tot_en*hartree, " eV"
        write(use_unit,*)"------------------------------------------", &
        "----------------------------------"
        write(use_unit,*)
     endif
  else if (use_os_rpa_qpe) then


     rpa_energy_os = 0.d0
     rpa_energy_ss = 0.d0
     rpa_energy_tot = 0.d0
     rpa_energy_scs = 0.d0
     if(use_2d_corr) then
         call evaluate_qpe_osrpa_energy_2 ( n_electrons, &
                                          n_low_state, n_high_state, &
                                          n_homo, n_lumo, &
                                          occ_numbers, n_full_freq, &
                                          omega_full_grid, womega_full, &
                                          chemical_potential, &
                                          KS_eigenvalue, KS_eigenvector, &
                                          ovlp_3KS, &
                                          rpa_energy_os, rpa_energy_ss, &
                                          rpa_energy_tot, rpa_energy_scs &
                                          )
     else
         call evaluate_qpe_osrpa_energy ( n_electrons, &
                                          n_low_state, n_high_state, &
                                          n_homo, n_lumo, &
                                          occ_numbers, n_full_freq, &
                                          omega_full_grid, womega_full, &
                                          chemical_potential, &
                                          KS_eigenvalue, KS_eigenvector, &
                                          ovlp_3KS, &
                                          rpa_energy_os, rpa_energy_ss, &
                                          rpa_energy_tot, rpa_energy_scs &
                                          )
     endif


     rpa_tot_en = total_energy-en_xc + &
     eex_energy + rpa_energy_os

     rpa_tot_en_tot = rpa_tot_en - rpa_energy_os + rpa_energy_tot
     rpa_tot_en_scs = rpa_tot_en - rpa_energy_os + rpa_energy_scs

     post_scf_total_energy = rpa_tot_en

     if(myid.eq.0) then
        write(use_unit,*) "--------------------------------------------------------------------"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exact exchange energy             : ",  eex_energy_real, &
        " Ha,", eex_energy_real*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osRPA correlation energy          : ",  rpa_energy_os, &
        " Ha,", rpa_energy_os*hartree, " eV"
        write(use_unit,*)
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  DFT/HF total energy               : ",  total_energy, &
        " Ha,", total_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exchange-only total energy        : ",  total_energy-en_xc+eex_energy, &
        " Ha,", (total_energy-en_xc+eex_energy)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osRPA total energy                : ",  rpa_tot_en, &
        " Ha,", rpa_tot_en*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA total energy                  : ",  rpa_tot_en_tot, &
        " Ha,", rpa_tot_en_tot*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  scsRPA total energy               : ",  rpa_tot_en_scs, &
        " Ha,", rpa_tot_en_scs*hartree, " eV"
        write(use_unit,*)"------------------------------------------", &
        "----------------------------------"
        write(use_unit,*)
     endif
  else if (use_sic_rpa_qpe) then
     rpa_energy = 0.d0
     call evaluate_qpe_sicrpa_energy &
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

     rpa_tot_en = total_energy-en_xc + &
     eex_energy + rpa_energy

     post_scf_total_energy = rpa_tot_en

     if(myid.eq.0) then
        write(use_unit,*) "--------------------------------------------------------------------"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exact exchange energy             : ",  eex_energy_real, &
        " Ha,", eex_energy_real*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osRPA correlation energy          : ",  rpa_energy, &
        " Ha,", rpa_energy*hartree, " eV"
        write(use_unit,*)
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  DFT/HF total energy               : ",  total_energy, &
        " Ha,", total_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exchange-only total energy        : ",  total_energy-en_xc+eex_energy, &
        " Ha,", (total_energy-en_xc+eex_energy)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  osRPA total energy                : ",  rpa_tot_en, &
        " Ha,", rpa_tot_en*hartree, " eV"
        write(use_unit,*)"------------------------------------------", &
        "----------------------------------"
        write(use_unit,*)
     endif
  endif

end subroutine qpe_ospt2
!******
