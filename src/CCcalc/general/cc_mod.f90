Module CC_corr

!  PURPOSE
!
!  This file may at some point contain all the needed variable declarations and
!  subroutines for the Coupled-Cluster (CC) FCI calculations
!
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!    This file was written by Tonghao Shen
!  SOURCE


  use dimensions
  use runtime_choices
  use physics
  use species_data
  use constants
  use gw_para
  use grids
  use geometry
  use mpi_utilities
  use synchronize_mpi
  use sbt_overlap_aims
  use localorb_io
  use basis
  use prodbas
  use hartree_fock
  use mpi_tasks
! ELPA commented out (not used anyway).
! To use ELPA, do one of the following:
! - use elpa*_2013   (2013 version)
! - use elsi         (2016/2017 version)
! Simply putting "use elpa*" won't work.
!  use elpa1
!  use elpa2
  use timing

  Implicit None

  Integer , Parameter :: CC_ip = Selected_Int_kind(8)

  Integer :: CC_valence
! The starting index of valence state

  Integer (kind = CC_ip) :: CC_n_config
! Number of configurations in CC calculations

  Integer (kind = CC_ip) :: CC_n_1a
! Number of singlet alpha excitations in CC calculations

  Integer (kind = CC_ip) :: CC_n_1b
! Number of singlet beta excitations in CC calculations

  Integer (kind = CC_ip) :: CC_n_ab
! Number of alpha-beta double-excitations in CC calculations

  Integer (kind = CC_ip) :: CC_n_2a
! Number of alpha-alpha double-excitations in CC calculations

  Integer (kind = CC_ip) :: CC_n_2b
! Number of beta-beta double-excitations in CC calculations

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_mem_w
! Memory distribution table
! CC_mem_w(:,1) --- single excitation in alpha space
! CC_mem_w(:,2) --- single excitation in beta space
! CC_mem_w(:,3) --- double excitation in alpha-beta space
! CC_mem_w(:,4) --- double excitation in alpha-alpha space
! CC_mem_w(:,5) --- double excitation in beta-beta space

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_index_w
! Global index of saving vectors
! CC_index_w(:,1) --- single excitation in alpha space
! CC_index_w(:,2) --- single excitation in beta space
! CC_index_w(:,3) --- double excitation in alpha-beta space
! CC_index_w(:,4) --- double excitation in alpha-alpha space
! CC_index_w(:,5) --- double excitation in beta-beta space


  Integer , dimension(2) :: CC_n_elec
! Numbers of alpha and beta electrons

  Integer , dimension(2) :: CC_n_occ
! Numbers of occupied(valence) orbitals

  Integer , dimension(2) :: CC_n_vir
! Numbers of virtual orbitals

  Integer , dimension(:,:) , allocatable :: CC_config_ref
! Reference (HF) configuration

  Double precision :: E_HF, CC_E_corr
! E_HF : HF energy
! CC_E_corr : Correlation energy

  Double precision , dimension(:,:) , allocatable :: CC_t_1a, CC_t_1b
  Double precision , dimension(:) , allocatable :: CC_t_2a, CC_t_ab, CC_t_2b
! Coefficient vector for CC excitations

  Double precision , dimension(:) , allocatable :: CC_w_1a, CC_w_1b, CC_w_2a, &
                                                   CC_w_ab, CC_w_2b
! Coefficient vector for CC excitations
  
  Double precision , dimension(:,:) , allocatable :: CC_t_1a_sv, CC_t_1b_sv, &
                                                     CC_t_2a_sv, CC_t_ab_sv, &
                                                     CC_t_2b_sv
! Coefficient vector saved for RLE calculation

  Double precision , dimension(:,:,:,:) , allocatable :: ovlp_3ks
! Three-center integrals (P|ij)

  Double precision , dimension(:,:) , allocatable :: CC_h_ij_a, CC_h_ij_b, &
                                                     CC_h_ba_a, CC_h_ba_b, &
                                                     CC_h_bj_a, CC_h_bj_b
! Auxiliary vectors

  Double precision , dimension(:,:) , allocatable :: CC_g_ca_a, CC_g_ca_b, &
                                                     CC_g_ik_a, CC_g_ik_b
! Auxiliary vectors

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_mem_h_ij, CC_mem_h_ba, &
                                           CC_mem_h_bj, CC_mem_g_ca, CC_mem_g_ik
! Memory distribution of auxiliary vectors
! CC_mem_h_??(:,1) --- alpha spin
! CC_mem_h_??(:,2) --- beta spin
! CC_mem_g_??(:,1) --- alpha spin
! CC_mem_g_??(:,2) --- beta spin

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_index_h_ij, CC_index_h_ba, &
                                         CC_index_h_bj, CC_index_g_ca, CC_index_g_ik
! Index of auxiliary vectors
! CC_index_h_??(:,1) --- alpha spin
! CC_index_h_??(:,2) --- beta spin
! CC_index_g_??(:,1) --- alpha spin
! CC_index_g_??(:,2) --- beta spin

  Double precision , dimension(:) , allocatable :: CC_aux_a_aa, CC_aux_a_ab, &
                                                   CC_aux_a_bb, CC_aux_b_aa, &
                                                   CC_aux_b_ab, CC_aux_b_bb, &
                                                   CC_aux_h_aaaa, CC_aux_h_abab, &
                                                   CC_aux_h_abba, CC_aux_h_baab, &
                                                   CC_aux_h_baba, CC_aux_h_bbbb
! Auxiliary vectors

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_mem_a, CC_mem_b, CC_mem_h
! Memory distribution of auxiliary vectors
! CC_mem_a(:,1) , CC_mem_b(:,1) ---- alpha-beta space
! CC_mem_a(:,2) , CC_mem_b(:,2) ---- alpha-alpha space
! CC_mem_a(:,3) , CC_mem_b(:,3) ---- beta-beta space
! CC_mem_h(:,1) ---- aaaa
! CC_mem_h(:,2) ---- abab
! CC_mem_h(:,3) ---- abba
! CC_mem_h(:,4) ---- baab
! CC_mem_h(:,5) ---- baba
! CC_mem_h(:,6) ---- bbbb

  Integer (kind = CC_ip) , dimension(:,:) , allocatable :: CC_index_a, CC_index_b, CC_index_h
! Index of auxiliary vectors
! CC_index_a(:,1) , CC_index_b(:,1) ---- alpha-beta space
! CC_index_a(:,2) , CC_index_b(:,2) ---- alpha-alpha space
! CC_index_a(:,3) , CC_index_b(:,3) ---- beta-beta space
! CC_index_h(:,1) ---- aaaa
! CC_index_h(:,2) ---- abab
! CC_index_h(:,3) ---- abba
! CC_index_h(:,4) ---- baab
! CC_index_h(:,5) ---- baba
! CC_index_h(:,6) ---- bbbb

  Integer :: CC_RLE_ndc, CC_RLE_m_bgn, CC_RLE_m_end
! CC_RLE_ndc ---- the dimension of current R matrix in RLE algorithm
! CC_RLE_m_bgn ---- the start point of saved basis functions
! CC_RLE_m_end ---- the end point of saved basis functions

  Double precision , dimension(:,:) , allocatable :: CC_RLE_t
! CC_RLE_t ---- saving the best evaluation of trial function of each iteration

  Logical :: CC_RLE_flag_singular
! the flag of sigularity of R matrix

  Logical :: CC_RLE_flag_sv_overflow
! the flag of sigularity of R matrix

End Module CC_corr
