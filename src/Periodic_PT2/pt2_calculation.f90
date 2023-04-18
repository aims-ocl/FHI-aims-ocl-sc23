!****s*  FHI-aims/pt2_calculation
!  NAME
!    pt2_calculation
!  SYNOPSIS

subroutine pt2_calculation (pt2_tot_en)

  !  PURPOSE
  !
  !  This subroutine performs periodic second-order perturbation calculation,
  !  based on SCF converged orbitals of DFT, HF or hybrid functionals.
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
  use cpt2_kspace_mod
  use cpt2_blacs, only: eex_energy_cpt2
  use synchronize_mpi, only: sync_timing
  use localorb_io, only: use_unit
  use pbc_lists, only: k_weights

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
  !
  !  Output
     real(kind=8)  pt2_tot_en

  real(kind=8)  eex_energy
  real(kind=8)  eex_energy_real
  real(kind=8)  energy_xc
  real(kind=8)  pt2_c_energy
  real(kind=8)  pt2_c_energy_os
  real(kind=8)  pt2_c_energy_ss

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
  character(*), parameter :: func = 'pt2_calculation_p0'


  if(use_mp2 .or. use_dftpt2 .or. use_os_mp2) then
      if(myid.eq.0) then
       write(use_unit,'(A)')"--------------------------------------------"
       write(use_unit,'(2X,A)') "Periodic PT2 total energy calculation starts ..."
      endif

      if (.not. use_mp2_blacs) then
        energy_xc = 0.d0
        call get_exchange_energy_p0 &
           (hybrid_coeff, KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue,&
            k_weights, occ_numbers,chemical_potential,eex_energy, energy_xc)
        eex_energy_real = eex_energy
      else
        eex_energy_real = eex_energy_cpt2
        eex_energy = eex_energy_cpt2
      end if

      if(use_hartree_fock) then
        eex_energy =  (1.d0-2.d0*hybrid_coeff)* eex_energy
      endif
      if(use_hybrid_meta) then
        en_hf = hybrid_coeff*eex_energy
      endif
      if(use_cohsex) then
        eex_energy = 2.d0*en_xc -eex_energy
      endif


      pt2_c_energy    = 0.0d0
      pt2_c_energy_os = 0.0d0
      pt2_c_energy_ss = 0.0d0
      if (use_os_mp2) then
          call evaluate_cpt2_os_energy_kspace &
             (occ_numbers, &
              pt2_c_energy_os   &
             )
      else if (use_mp2_blacs) then
          call evaluate_cpt2_energy_kspace_blacs &
             (occ_numbers, &
              pt2_c_energy,  &
              pt2_c_energy_os,  &
              pt2_c_energy_ss   &
             )
      else
          call evaluate_cpt2_energy_kspace &
             (occ_numbers, &
              pt2_c_energy,  &
              pt2_c_energy_os,  &
              pt2_c_energy_ss   &
             )
         endif

  endif

  if(myid.eq.0) then
     if(use_dftpt2) then
      write(use_unit,'(A)')
      write(use_unit,*) "--------------------------------------------------------------------"
      write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
      "| Exact exchange energy        :", eex_energy_real*dftpt2_Ex_hf, &
      " Ha,", eex_energy_real*dftpt2_Ex_hf*hartree, " eV"
      write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
      "| PT2 contribution             :", &
      (pt2_c_energy_os*dftpt2_Ec_osPT2+pt2_c_energy_ss*dftpt2_Ec_ssPT2), " Ha,",&
      (pt2_c_energy_os*dftpt2_Ec_osPT2+pt2_c_energy_ss*dftpt2_Ec_ssPT2)* hartree ," eV"
      write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
      "| DHDF/DFT Energy              :", total_energy, " Ha,", &
      total_energy* hartree," eV"
      write(use_unit,'(A)')"---------------------------------------------"
      pt2_tot_en = total_energy + eex_energy_real*dftpt2_Ex_hf &
                 + pt2_c_energy_os*dftpt2_Ec_osPT2 &
                 + pt2_c_energy_ss*dftpt2_Ec_ssPT2
      post_scf_total_energy = pt2_tot_en
      if (flag_dftpt2_dft_part .eq. 30) then
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "| Total XYG3 energy            :", pt2_tot_en, " Ha,", &
        pt2_tot_en* hartree ," eV"
      elseif (flag_dftpt2_dft_part .eq. 31) then
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "| Total xDH-PBE0 energy        :", pt2_tot_en, " Ha,", &
        pt2_tot_en* hartree ," eV"
      elseif (flag_dftpt2_dft_part .eq. 32) then
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "| Total XYGJOS energy          :", pt2_tot_en, " Ha,", &
        pt2_tot_en* hartree ," eV"
      elseif (flag_dftpt2_dft_part .eq. 33) then
        write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
        "| Total ZRS1 energy            :", pt2_tot_en, " Ha,", &
        pt2_tot_en* hartree ," eV"
      endif
      write(use_unit,'(A)')"---------------------------------------------"
      !post_scf_total_energy = total_energy + E_mp2*dftpt2_Ec_PT2
    else
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
      pt2_tot_en = total_energy - en_xc +  eex_energy + pt2_c_energy
      write(use_unit,'(2X,A,f21.8,A,f21.8,A)') &
      "  PT2 total energy             : ",  pt2_tot_en, &
      " Ha,", pt2_tot_en*hartree, " eV"
      write(use_unit,*)
    endif
  endif

end subroutine pt2_calculation
!******
