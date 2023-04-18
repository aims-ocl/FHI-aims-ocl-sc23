!****s*  FHI-aims/rpa_calculation_p0
!  NAME
!    rpa_calculation_p0
!  SYNOPSIS

subroutine rpa_calculation_p0 (rpa_tot_en)

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
  use gw_para
  use constants
  use mpi_tasks
  use timing
  use synchronize_mpi, only: sync_timing
  use localorb_io, only: use_unit
  use pbc_lists, only: k_point_list, k_weights

!<ADDED CODE>
!  use cRPA_flow
!</ADDED CODE>

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
     real*8  rpa_tot_en

  !  local arrays used in GW calculation :
  !    ovlp_3KS : the two single basis functions are transformed into KS states,
  !             this name is proper here anymore, to be changed in future
  !    xc_matr :  the LDA xc potential matrix between two single basis fns
  !    xc_KS_matr : the single basis of xc_matr transformed into KS states
  !    exchange_self_energy : (diagonal) matrix elements of the Fock operator
  !    self_energy_omega : the self energy on imaginary frequency axis
!  real*8, dimension(:,:,:), allocatable :: xc_matr
!  real*8, dimension(:,:,:), allocatable :: xc_KS_matr
!  real*8, dimension(:,:), allocatable :: exchange_self_energy
!  complex*16, dimension(:,:,:), allocatable :: self_energy_omega
!  complex*16, dimension(:,:,:), allocatable :: mp2_selfenergy
!  real*8, dimension(:,:), allocatable :: qp_energy

  real*8, dimension(:,:), allocatable :: k_point_list_old

  real*8  eex_energy
  real*8  eex_energy_real
  real*8  energy_xc
  real*8  c_energy
  real*8  c_energy_lda
  real*8  c_energy_ldarpa
  real*8  rpa_c_energy
  real*8  rpa_plus_energy
  real*8  rpa_plus_tot_en
  real*8  e_2oex
  real*8  e_sosex
  real*8  en_single_excit
  real*8  en_ladder_se

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
  integer i_index
  integer i_cell_n, i_cell_1, i_cell_2, i_cell_3, i_k_point

  !  error variables
  integer mpierr
  integer :: info
  character*150 :: info_str
  character(*), parameter :: func = 'rpa_calculation_p0'


  if(use_rpa_ene) then
       if(myid.eq.0) then
        write(use_unit,'(A)')"--------------------------------------------"
        write(use_unit,'(2X,A)') "Periodic RPA total energy calculation starts ..."
       endif

!       if(.not.allocated(k_phase_exx))then
!          ! Allocate and calculate k_phase_exx
!          write(info_str,'(2X,A,I10,A)') '| Consuming ', (n_cells_bvk*n_k_points*16)/2**10, ' KiB for k_phase_exx.'
!          call localorb_info(info_str,use_unit,'(A)')
!          allocate(k_phase_exx( n_cells_bvk, n_k_points),stat=info)
!          call check_allocation(info, 'k_phase_exx', func)
!          k_phase_exx = (1.d0,0.d0)
!          do i_k_point = 1, n_k_points
!             i_cell_n = 1
!             do i_cell_1 = -(n_k_points_xyz(1)-1)/2, n_k_points_xyz(1)/2
!                do i_cell_2 = -(n_k_points_xyz(2)-1)/2, n_k_points_xyz(2)/2
!                   do i_cell_3 = -(n_k_points_xyz(3)-1)/2, n_k_points_xyz(3)/2
!                      
!                      if (i_cell_1 /= 0 .or. i_cell_2 /= 0 .or. i_cell_3 /= 0) then
!                         i_cell_n = i_cell_n + 1
!                         k_phase_exx( i_cell_n, i_k_point) &
!                              & = k_phase_base(1, i_k_point)**i_cell_1 &
!                              & * k_phase_base(2, i_k_point)**i_cell_2 &
!                              & * k_phase_base(3, i_k_point)**i_cell_3
!                      endif
!                   enddo
!                enddo
!             enddo
!          enddo
!       endif

        energy_xc = 0.d0
        call get_exchange_energy_p0 &
           (hybrid_coeff, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue,&
            k_weights, occ_numbers,chemical_potential,eex_energy, energy_xc)

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

      flag_k_fine_grid = .false.
      if (flag_k_fine_grid) then
! prepare the LVL triple coefficients, KS eigenvectors, and kq_point_list for a fine k mesh
        n_k_points_original = n_k_points
        allocate(k_point_list_old(n_k_points, 3),stat=i_index)
        call check_allocation(i_index, 'k_point_list_old              ')
        k_point_list_old = k_point_list

        k_enhancement_factor(1) = 1
        k_enhancement_factor(2) = 1
        k_enhancement_factor(3) = 1
        n_k_points = n_k_points * k_enhancement_factor(1)*k_enhancement_factor(2)* &
                     k_enhancement_factor(3)

        call prepare_fine_k_quantities()

        if(allocated(kq_point_list)) then
         deallocate(kq_point_list)
        endif
        allocate(kq_point_list(n_k_points,n_k_points),stat=i_index)
        call check_allocation(i_index, 'kq_point_list')
        call determine_k_minus_q_list_finemesh(n_k_points_original,k_point_list_old,kq_point_list)
      endif

      if (restart_rpa_write) then ! For restarting RPA calculations
        call evaluate_crpa_energy_kspace_restart &
           (n_low_state,occ_numbers, n_full_freq, &
            omega_full_grid, womega_full, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            rpa_c_energy  &
           )
      else
        call evaluate_crpa_energy_kspace &
           (n_low_state,occ_numbers, n_full_freq, &
            omega_full_grid, womega_full, &
            KS_eigenvalue, KS_eigenvector, &
            KS_eigenvector_complex, &
            rpa_c_energy  &
           )
      endif
!<ADDED CODE>
!        CALL calc_crpa_energy(ubound(omega_full_grid,1), &
!                             omega_full_grid, &
!                             womega_full, &
!                             rpa_c_energy)
!</ADDED CODE>

        rpa_tot_en = total_energy - en_xc +  eex_energy + rpa_c_energy
        if(myid.eq.0) then
          write(use_unit,'(A,f18.8)') "total_energy=", total_energy*hartree
          write(use_unit,'(A,f18.8)') "en_xc=", en_xc*hartree
          write(use_unit,'(A,f18.8)') "eex_energy=", eex_energy*hartree
          write(use_unit,'(A,f18.8)') "rpa_tot_en=", rpa_tot_en*hartree
        endif

  endif

  if (allocated (k_point_list_old)) then
      deallocate (k_point_list_old)
  endif

  if(myid.eq.0) then
        write(use_unit,*) "--------------------------------------------------------------------"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exact exchange energy        : ", eex_energy_real, &
        " Ha,", eex_energy_real*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  DFT/HF total energy          : ",  total_energy, &
        " Ha,", total_energy*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  Exchange-only total energy   : ",  total_energy-en_xc+eex_energy, &
        " Ha,", (total_energy-en_xc+eex_energy)*hartree, " eV"
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "  RPA total energy             : ",  rpa_tot_en, &
        " Ha,", rpa_tot_en*hartree, " eV"
        write(use_unit,*)
!        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!        "  RPA+SE total energy          : ",  rpa_tot_en+en_single_excit, &
!        " Ha,", (rpa_tot_en+en_single_excit)*hartree, " eV"
!        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
!        "  RPA+RSE total energy         : ",  rpa_tot_en+en_ladder_se, &
!        " Ha,", (rpa_tot_en+en_ladder_se)*hartree, " eV"
!        write(use_unit,*)
!
!        write(use_unit,*)"------------------------------------------", &
!        "----------------------------------"
!        write(use_unit,*)
     endif

end subroutine rpa_calculation_p0
!******
