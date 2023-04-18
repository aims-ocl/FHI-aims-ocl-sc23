!****h* FHI-aims/timing
!  NAME
!    timing
!  SYNOPSIS
      module timing
!  Module timing handles all timing information collected across the code,
!  as well as warnings summarized at the end of a run.
!
!  WPH (2018 Jan 18):  I've split the timing module into two modules:
!  - a timing_core module which handles the time collection and output
!    formatting in an aims-independent fashion
!  - the (original) timing module which functions as a wrapper around
!    timing_core and keeps track of aims-specific timings.
!  This change was made to allow other aims-independent functionality to access
!  aims' core timing subroutines without introducing dependencies on core
!  aims-specific modules such as runtime_choices and dimensions.
! USES
        use timing_core
        implicit none
!  PURPOSE
!    This module takes care of everything related to MD and requires only energies and forces
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
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

!     global variable declarations - exported to other program parts
!     these are the individual time stamps as collected throughout the code

      ! CPU time counters
      real*8       :: time_zero
      real*8       :: time_total
      real*8       :: time_sc_loop
      real*8       :: time_prep
      real*8       :: time_free_offset
      real*8       :: time_on_evs
      real*8       :: time_intg
      real*8       :: time_diag
      real*8       :: time_density
      real*8       :: time_mixing
      real*8       :: time_hartree_multi
      real*8       :: time_hartree_sum
      real*8       :: time_etot
      real*8       :: time_prodbas_total
      real*8       :: time_ovlp3fn
      real*8       :: time_ovlp3fn_to_ovlp3KS
      real*8       :: time_lvl_triple_recip
      real*8       :: time_coulomb_ten
      real*8       :: time_coulomb_matr
      real*8       :: time_inv_coulomb_matr
      real*8       :: time_ovlp_multi
      real*8       :: time_ovlp_swap
      real*8       :: time_fock
      real*8       :: time_total_gw
      real*8       :: time_self_energy
      real*8       :: time_Wmn_imagfreq
      real*8       :: time_polar_imagfreq
      real*8       :: time_band_self_energy
      real*8       :: time_2ox_selfenergy
      real*8       :: time_s2ox_selfenergy
      real*8       :: time_soxw_selfenergy
      real*8       :: time_polar
      real*8       :: time_qp_energy_cd
      real*8       :: time_self_energy_cd
      real*8       :: time_Wmn_realfreq_cd
      real*8       :: time_WPQ_realfreq_cd
      real*8       :: time_polar_realfreq_cd
      real*8       :: time_print_sigma_cd
      real*8       :: time_spec_func_cd
      real*8       :: time_rpa_corr
      real*8       :: time_rse_corr
      real*8       :: time_2oex_corr
      real*8       :: time_sosex_corr
      real*8       :: time_C6_coef
      real*8       :: time_scaled_zora
      real*8       :: time_partition
      real*8       :: time_superpos
      real*8       :: time_max_scf_cycle
      real*8       :: time_max_postprocessing
      real*8       :: time_omegabuild
      real*8       :: time_omegaresdist
      real*8       :: time_total_neutral_excitation
      real*8       :: time_total_bse
      real*8       :: time_bse_rot
      real*8       :: time_bse_solve_mat
      real*8       :: time_bse_build_V
      real*8       :: time_bse_build_W
      real*8       :: time_bse_ovlp_3fn
      real*8       :: time_bse_ovlp_3ks
      real*8       :: time_bse_init
      real*8       :: time_bse_w0
      real*8       :: time_bse_w_ovlp_3ks
      real*8       :: time_bse_prod_bas
      real*8       :: time_bse_construct_mat
      real*8       :: time_bse_osc
      real*8       :: time_omega_osc
      real*8       :: time_neutral_excitation_orbitals
      real*8       :: time_neutral_excitation_casida
      real*8       :: time_neutral_excitation_dipole
      real*8       :: time_fxc_integr
      real*8       :: time_3ks_comm
      real*8       :: time_coeffs_fxc3fn
      real*8       :: time_casida_matrix_calculation
      real*8       :: time_matvec_cont
      real*8       :: time_omega_solve
      logical      :: time_initial_max_scf_cycle
      logical      :: time_initial_max_postprocessing
      real*8       :: time_sc_tot !time for sc convergence in geometry optimization or MD
      real*8       :: time_magnetic_response
      real*8       :: time_bsse ! BSSE timings
      real*8       :: time_embedding
      real*8       :: time_vdw
      real*8       :: time_cpscf_loop
      real*8       :: time_first_order_density
      real*8       :: time_first_order_potential
      real*8       :: time_first_order_H
      real*8       :: time_first_order_DM
      real*8       :: time_Sternheimer
      real*8       :: time_first_order_S
      real*8       :: time_Hessian
      real*8       :: time_bc_init
      real*8       :: time_gpu_init
      real*8       :: time_fodft_combine
      real*8       :: time_fodft_select
      real*8       :: time_fodft_lowdin
      real*8       :: time_fodft_output

      real*8       :: tot_time_intg
      real*8       :: tot_time_diag
      real*8       :: tot_time_density
      real*8       :: tot_time_mixing
      real*8       :: tot_time_hartree_multi
      real*8       :: tot_time_hartree_sum
      real*8       :: tot_time_etot
      real*8       :: tot_time_fock
      real*8       :: tot_time_scaled_zora
      real*8       :: tot_time_wf_extra_in
      real*8       :: tot_time_wf_extra_out
      real*8       :: tot_time_partition
      real*8       :: tot_time_free_offset
      real*8       :: tot_time_on_evs
      real*8       :: tot_time_superpos
      real*8       :: time_forces_orbital_update ! stores the time of the first force update
      real*8       :: tot_time_embedding
      real*8       :: tot_time_vdw
      real*8       :: tot_time_soc
      real*8       :: tot_time_band_dos
      real*8       :: tot_time_mulliken
      real*8       :: tot_time_cube_output
      real*8       :: tot_time_dielectric
      real*8       :: tot_time_prodbas_total
      real*8       :: tot_time_ovlp3fn
      real*8       :: tot_time_coulomb_matr
      real*8       :: tot_time_inv_coulomb_matr
      real*8       :: tot_time_ovlp_multi
      real*8       :: tot_time_lvl_triple_recip
      real*8       :: tot_time_first_order_density
      real*8       :: tot_time_first_order_potential
      real*8       :: tot_time_first_order_H
      real*8       :: tot_time_first_order_DM
      real*8       :: tot_time_Sternheimer
      real*8       :: tot_time_first_order_S
      real*8       :: tot_time_Hessian
      real*8       :: tot_time_fodft
      real*8       :: tot_time_bc_init
      real*8       :: tot_time_gpu_init
      real*8       :: tot_time_matrix_io

      ! Wall clock time counters
      real*8       :: clock_time_zero
      real*8       :: clock_time_total
      real*8       :: clock_time_sc_loop
      real*8       :: clock_time_prep
      real*8       :: clock_time_intg
      real*8       :: clock_time_diag
      real*8       :: clock_time_partition
      real*8       :: clock_time_free_offset
      real*8       :: clock_time_on_evs
      real*8       :: clock_time_superpos
      real*8       :: clock_time_etot
      real*8       :: clock_time_prodbas_total
      real*8       :: clock_time_ovlp3fn
      real*8       :: clock_time_ovlp3fn_to_ovlp3KS
      real*8       :: clock_time_lvl_triple_recip
      real*8       :: clock_time_coulomb_ten
      real*8       :: clock_time_coulomb_matr
      real*8       :: clock_time_inv_coulomb_matr
      real*8       :: clock_time_ovlp_multi
      real*8       :: clock_time_ovlp_swap
      real*8       :: clock_time_total_gw
      real*8       :: clock_time_self_energy
      real*8       :: clock_time_Wmn_imagfreq
      real*8       :: clock_time_polar_imagfreq
      real*8       :: clock_time_band_self_energy
      real*8       :: clock_time_2ox_selfenergy
      real*8       :: clock_time_s2ox_selfenergy
      real*8       :: clock_time_soxw_selfenergy
      real*8       :: clock_time_polar
      real*8       :: clock_time_qp_energy_cd
      real*8       :: clock_time_self_energy_cd
      real*8       :: clock_time_Wmn_realfreq_cd
      real*8       :: clock_time_WPQ_realfreq_cd
      real*8       :: clock_time_polar_realfreq_cd
      real*8       :: clock_time_print_sigma_cd
      real*8       :: clock_time_spec_func_cd
      real*8       :: clock_time_rpa_corr
      real*8       :: clock_time_rse_corr
      real*8       :: clock_time_2oex_corr
      real*8       :: clock_time_sosex_corr
      real*8       :: clock_time_C6_coef
      real*8       :: clock_time_density
      real*8       :: clock_time_mixing
      real*8       :: clock_time_fock
      real*8       :: clock_time_hartree_multi
      real*8       :: clock_time_hartree_sum
      real*8       :: clock_time_scaled_zora
      real*8       :: clock_time_sc_tot
      real*8       :: clock_time_bsse ! BSSE timings
      real*8       :: clock_time_embedding
      real*8       :: clock_time_omegabuild
      real*8       :: clock_time_omegaresdist
      real*8       :: clock_time_total_neutral_excitation
      real*8       :: clock_time_total_bse
      real*8       :: clock_time_bse_rot
      real*8       :: clock_time_bse_solve_mat
      real*8       :: clock_time_bse_build_V
      real*8       :: clock_time_bse_build_W
      real*8       :: clock_time_bse_ovlp_3fn
      real*8       :: clock_time_bse_ovlp_3ks
      real*8       :: clock_time_bse_init
      real*8       :: clock_time_bse_w0
      real*8       :: clock_time_bse_w_ovlp_3ks
      real*8       :: clock_time_bse_prod_bas
      real*8       :: clock_time_bse_construct_mat
      real*8       :: clock_time_bse_osc
      real*8       :: clock_time_omega_osc
      real*8       :: clock_time_neutral_excitation_orbitals
      real*8       :: clock_time_neutral_excitation_casida
      real*8       :: clock_time_neutral_excitation_dipole
      real*8       :: clock_time_fxc_integr
      real*8       :: clock_time_3ks_comm
      real*8       :: clock_time_coeffs_fxc3fn
      real*8       :: clock_time_casida_matrix_calculation
      real*8       :: clock_time_matvec_cont
      real*8       :: clock_time_omega_solve
      real*8       :: clock_time_vdw
      real*8       :: clock_time_cpscf_loop
      real*8       :: clock_time_first_order_density
      real*8       :: clock_time_first_order_potential
      real*8       :: clock_time_first_order_H
      real*8       :: clock_time_first_order_DM
      real*8       :: clock_time_Sternheimer
      real*8       :: clock_time_first_order_S
      real*8       :: clock_time_Hessian
      real*8       :: clock_time_bc_init
      real*8       :: clock_time_gpu_init
      real*8       :: clock_time_magnetic_response
      real*8       :: clock_time_fodft_combine
      real*8       :: clock_time_fodft_select
      real*8       :: clock_time_fodft_lowdin
      real*8       :: clock_time_fodft_output

      real*8       :: tot_clock_time_intg
      real*8       :: tot_clock_time_diag
      real*8       :: tot_clock_time_partition
      real*8       :: tot_clock_time_free_offset
      real*8       :: tot_clock_time_on_evs
      real*8       :: tot_clock_time_superpos
      real*8       :: tot_clock_time_etot
      real*8       :: tot_clock_time_2oex_corr
      real*8       :: tot_clock_time_sosex_corr
      real*8       :: tot_clock_time_density
      real*8       :: tot_clock_time_mixing
      real*8       :: tot_clock_time_fock
      real*8       :: tot_clock_time_hartree_multi
      real*8       :: tot_clock_time_hartree_sum
      real*8       :: tot_clock_time_scaled_zora
      real*8       :: tot_clock_time_wf_extra_in
      real*8       :: tot_clock_time_wf_extra_out
      real*8       :: tot_clock_time_bsse  ! BSSE timings
      real*8       :: tot_clock_time_embedding
      real*8       :: tot_clock_time_vdw
      real*8       :: tot_clock_time_soc
      real*8       :: tot_clock_time_band_dos
      real*8       :: tot_clock_time_mulliken
      real*8       :: tot_clock_time_cube_output
      real*8       :: tot_clock_time_dielectric
      real*8       :: tot_clock_time_prodbas_total
      real*8       :: tot_clock_time_ovlp3fn
      real*8       :: tot_clock_time_coulomb_matr
      real*8       :: tot_clock_time_inv_coulomb_matr
      real*8       :: tot_clock_time_ovlp_multi
      real*8       :: tot_clock_time_lvl_triple_recip
      real*8       :: tot_clock_time_first_order_density
      real*8       :: tot_clock_time_first_order_potential
      real*8       :: tot_clock_time_first_order_H
      real*8       :: tot_clock_time_first_order_DM
      real*8       :: tot_clock_time_Sternheimer
      real*8       :: tot_clock_time_first_order_S
      real*8       :: tot_clock_time_Hessian
      real*8       :: tot_clock_time_bc_init
      real*8       :: tot_clock_time_gpu_init
      real*8       :: tot_clock_time_fodft
      real*8       :: tot_clock_time_matrix_io

      ! auxiliary variables to call cpu_time() in different parts of the code
      real*8       :: rtime
      real*8       :: time_ini
      real*8       :: time_fin
      real*8       :: clock_rtime
      real*8       :: clock_time_ini
      real*8       :: clock_time_fin

      ! Global loop counters
      integer      :: number_of_loops
      integer      :: total_number_of_loops
      integer      :: relaxation_step_number
      integer      :: MD_stepcount
      integer      :: MD_force_evaluations
      integer      :: current_atom_num_bsse
      integer      :: number_of_scf_inits

      ! XZL: added for PIMD
      integer      :: PIMD_stepcount
      integer      :: PIMD_force_evaluations
      integer      :: PIMD_beadcount

      ! memory checking
      integer      :: initial_memory
      integer      :: available_memory

      ! variables for final_warnings
      ! (only to be used in this module, no other warnings please)
      ! ATTENTION: Must initialize all warnings explicitly below in
      !            initialize_timings. Else, library-type uses of
      !            FHI-aims may not work.
      !
      logical      :: warn_integrals
      logical      :: warn_idle_cpus
      logical      :: warn_slow_scf_convergence
      logical      :: warn_mixing_parameter_range
      logical      :: warn_cube_file_format
      logical      :: warn_RI_V_grids
      logical      :: warn_cpu_consistency
      logical      :: warn_surface_dipole
      logical      :: warn_mpi_in_place
      logical      :: warn_qpe_nonconvergence
      ! Have GPUs *always* been used for the specified operation?
      ! (We set these variables to .true. initially, and change them to .false.
      ! whenever GPUs are not used)
      logical      :: gpu_density_always_used
      logical      :: gpu_hamiltonian_always_used
      logical      :: gpu_forces_always_used

      ! private variables
      character*120 , private :: info_str

!     here follows the actual work

        contains
!******

!****s* timing/initial_timings
!  NAME
!    initial_timings
!  SYNOPSIS
     subroutine initial_timings
!  PURPOSE
!    print starting time statement
!  USES
      use localorb_io, only: localorb_info, use_unit
!  INPUT
!  none
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
      ! output_datetime contains code that was originally performed here, but
      ! was moved into timing_core as it could have general utility
      call output_datetime(time_zero, clock_time_zero)

      ! initialize global timers for entire run
      time_total = 0.d0
      time_sc_loop = 0.d0
      time_prep = 0.d0
      time_free_offset = 0.d0
      time_on_evs = 0.d0
      time_intg = 0.d0
      time_diag = 0.d0
      time_density = 0.d0
      time_mixing = 0.d0
      time_hartree_multi = 0.d0
      time_hartree_sum = 0.d0
      time_etot = 0.d0
      time_prodbas_total = 0.d0
      time_ovlp3fn = 0.d0
      time_ovlp3fn_to_ovlp3KS = 0.d0
      time_lvl_triple_recip = 0.d0
      time_coulomb_ten = 0.d0
      time_coulomb_matr = 0.d0
      time_inv_coulomb_matr = 0.d0
      time_ovlp_multi = 0.d0
      time_ovlp_swap = 0.d0
      time_fock = 0.d0
      time_total_gw = 0.d0
      time_self_energy = 0.d0
      time_Wmn_imagfreq = 0.d0
      time_polar_imagfreq = 0.d0
      time_band_self_energy = 0.d0
      time_2ox_selfenergy = 0.d0
      time_s2ox_selfenergy = 0.d0
      time_soxw_selfenergy = 0.d0
      time_polar = 0.d0
      time_qp_energy_cd = 0.0d0
      time_self_energy_cd = 0.0d0
      time_Wmn_realfreq_cd = 0.0d0
      time_WPQ_realfreq_cd = 0.0d0
      time_polar_realfreq_cd = 0.0d0
      time_print_sigma_cd = 0.0d0
      time_spec_func_cd = 0.0d0
      time_rpa_corr = 0.d0
      time_rse_corr = 0.d0
      time_2oex_corr = 0.d0
      time_sosex_corr = 0.d0
      time_C6_coef = 0.d0
      time_scaled_zora = 0.d0
      time_partition = 0.d0
      time_superpos = 0.d0
      time_max_scf_cycle = 0.d0
      time_max_postprocessing = 0.d0
      time_omegabuild = 0d0
      time_omegaresdist = 0d0
      time_total_neutral_excitation = 0d0
      time_total_bse = 0d0
      time_bse_rot = 0d0
      time_bse_solve_mat = 0d0
      time_bse_build_V = 0d0
      time_bse_build_W = 0d0
      time_bse_ovlp_3fn = 0d0
      time_bse_ovlp_3ks = 0d0
      time_bse_init = 0d0
      time_bse_w0 = 0d0
      time_bse_w_ovlp_3ks = 0d0
      time_bse_prod_bas = 0d0
      time_bse_construct_mat = 0d0
      time_bse_osc = 0d0
      time_omega_osc = 0d0
      time_neutral_excitation_orbitals = 0d0
      time_neutral_excitation_casida = 0d0
      time_neutral_excitation_dipole = 0d0
      time_fxc_integr = 0d0
      time_3ks_comm = 0d0
      time_coeffs_fxc3fn = 0d0
      time_casida_matrix_calculation = 0d0
      time_matvec_cont = 0d0
      time_omega_solve = 0d0
      time_sc_tot = 0.d0
      time_magnetic_response = 0d0
      time_bsse = 0.0d0
      time_embedding = 0.d0
      time_vdw = 0.d0
      time_cpscf_loop = 0.d0
      time_first_order_density = 0.d0
      time_first_order_potential = 0.d0
      time_first_order_H = 0.d0
      time_first_order_DM = 0.d0
      time_Sternheimer = 0.d0
      time_first_order_S = 0.d0
      time_Hessian = 0.d0
      time_bc_init = 0.d0
      time_gpu_init = 0.0d0
      time_fodft_combine = 0.d0
      time_fodft_select = 0.d0
      time_fodft_lowdin = 0.d0
      time_fodft_output = 0.d0

      tot_time_intg = 0.d0
      tot_time_diag = 0.d0
      tot_time_density = 0.d0
      tot_time_mixing = 0.d0
      tot_time_hartree_multi = 0.d0
      tot_time_hartree_sum = 0.d0
      tot_time_etot = 0.d0
      tot_time_fock = 0.d0
      tot_time_scaled_zora = 0.d0
      tot_time_wf_extra_in = 0.d0
      tot_time_wf_extra_out = 0.d0
      tot_time_partition = 0.d0
      tot_time_free_offset = 0.d0
      tot_time_on_evs = 0.d0
      tot_time_superpos = 0.d0
      time_forces_orbital_update = 0.d0
      tot_time_embedding = 0.d0
      tot_time_vdw = 0.d0
      tot_time_soc = 0.d0
      tot_time_band_dos = 0.d0
      tot_time_mulliken = 0.d0
      tot_time_cube_output = 0.d0
      tot_time_dielectric = 0.d0
      tot_time_prodbas_total = 0.d0
      tot_time_ovlp3fn = 0.d0
      tot_time_coulomb_matr = 0.d0
      tot_time_inv_coulomb_matr = 0.d0
      tot_time_ovlp_multi = 0.d0
      tot_time_lvl_triple_recip = 0.d0
      tot_time_first_order_density = 0.d0
      tot_time_first_order_potential = 0.d0
      tot_time_first_order_H = 0.d0
      tot_time_first_order_DM = 0.d0
      tot_time_Sternheimer = 0.d0
      tot_time_first_order_S = 0.d0
      tot_time_Hessian = 0.d0
      tot_time_fodft = 0.d0
      tot_time_bc_init = 0.d0
      tot_time_gpu_init = 0.0d0
      tot_time_matrix_io = 0.d0

      clock_time_total = 0.d0
      clock_time_sc_loop = 0.d0
      clock_time_prep = 0.d0
      clock_time_intg = 0.d0
      clock_time_diag = 0.d0
      clock_time_partition = 0.d0
      clock_time_free_offset = 0.d0
      clock_time_on_evs = 0.d0
      clock_time_superpos = 0.d0
      clock_time_etot = 0.d0
      clock_time_prodbas_total = 0.d0
      clock_time_ovlp3fn = 0.d0
      clock_time_ovlp3fn_to_ovlp3KS = 0.d0
      clock_time_lvl_triple_recip = 0.d0
      clock_time_coulomb_ten = 0.d0
      clock_time_coulomb_matr = 0.d0
      clock_time_inv_coulomb_matr = 0.d0
      clock_time_ovlp_multi = 0.d0
      clock_time_ovlp_swap = 0.d0
      clock_time_total_gw = 0.d0
      clock_time_self_energy = 0.d0
      clock_time_Wmn_imagfreq = 0.d0
      clock_time_polar_imagfreq = 0.d0
      clock_time_band_self_energy = 0.d0
      clock_time_2ox_selfenergy = 0.d0
      clock_time_s2ox_selfenergy = 0.d0
      clock_time_soxw_selfenergy = 0.d0
      clock_time_polar = 0.d0
      clock_time_qp_energy_cd = 0.d0
      clock_time_self_energy_cd = 0.d0
      clock_time_Wmn_realfreq_cd = 0.d0
      clock_time_WPQ_realfreq_cd = 0.d0
      clock_time_polar_realfreq_cd = 0.d0
      clock_time_print_sigma_cd = 0.d0
      clock_time_spec_func_cd = 0.d0
      clock_time_rpa_corr = 0.d0
      clock_time_rse_corr = 0.d0
      clock_time_2oex_corr = 0.d0
      clock_time_sosex_corr = 0.d0
      clock_time_C6_coef = 0.d0
      clock_time_density = 0.d0
      clock_time_mixing = 0.d0
      clock_time_fock = 0.d0
      clock_time_hartree_multi = 0.d0
      clock_time_hartree_sum = 0.d0
      clock_time_scaled_zora = 0.d0
      clock_time_sc_tot = 0.d0
      clock_time_bsse = 0.d0
      clock_time_embedding = 0.d0
      clock_time_omegabuild = 0d0
      clock_time_omegaresdist = 0d0
      clock_time_total_neutral_excitation = 0d0
      clock_time_total_bse = 0d0
      clock_time_bse_rot = 0d0
      clock_time_bse_solve_mat = 0d0
      clock_time_bse_build_V = 0d0
      clock_time_bse_build_W = 0d0
      clock_time_bse_ovlp_3fn = 0d0
      clock_time_bse_ovlp_3ks = 0d0
      clock_time_bse_init = 0d0
      clock_time_bse_w0 = 0d0
      clock_time_bse_w_ovlp_3ks = 0d0
      clock_time_bse_prod_bas = 0d0
      clock_time_bse_construct_mat = 0d0
      clock_time_bse_osc = 0d0
      clock_time_omega_osc = 0d0
      clock_time_neutral_excitation_orbitals = 0d0
      clock_time_neutral_excitation_casida = 0d0
      clock_time_neutral_excitation_dipole = 0d0
      clock_time_fxc_integr = 0d0
      clock_time_3ks_comm = 0d0
      clock_time_coeffs_fxc3fn = 0d0
      clock_time_casida_matrix_calculation = 0d0
      clock_time_matvec_cont = 0d0
      clock_time_omega_solve = 0d0
      clock_time_vdw = 0.d0
      clock_time_cpscf_loop = 0.d0
      clock_time_first_order_density = 0.d0
      clock_time_first_order_potential = 0.d0
      clock_time_first_order_H = 0.d0
      clock_time_first_order_DM = 0.d0
      clock_time_Sternheimer = 0.d0
      clock_time_first_order_S = 0.d0
      clock_time_Hessian = 0.d0
      clock_time_bc_init = 0.d0
      clock_time_gpu_init = 0.0d0
      clock_time_magnetic_response = 0d0
      clock_time_fodft_combine = 0.d0
      clock_time_fodft_select = 0.d0
      clock_time_fodft_lowdin = 0.d0
      clock_time_fodft_output = 0.d0

      tot_clock_time_intg = 0.d0
      tot_clock_time_diag = 0.d0
      tot_clock_time_partition = 0.d0
      tot_clock_time_free_offset = 0.d0
      tot_clock_time_on_evs = 0.d0
      tot_clock_time_superpos = 0.d0
      tot_clock_time_etot = 0.d0
      tot_clock_time_2oex_corr = 0.d0
      tot_clock_time_sosex_corr = 0.d0
      tot_clock_time_density = 0.d0
      tot_clock_time_mixing = 0.d0
      tot_clock_time_fock = 0.d0
      tot_clock_time_hartree_multi = 0.d0
      tot_clock_time_hartree_sum = 0.d0
      tot_clock_time_scaled_zora = 0.d0
      tot_clock_time_wf_extra_in = 0.d0
      tot_clock_time_wf_extra_out = 0.d0
      tot_clock_time_bsse = 0.d0
      tot_clock_time_embedding = 0.d0
      tot_clock_time_vdw = 0.d0
      tot_clock_time_soc = 0.0d0
      tot_clock_time_band_dos = 0.d0
      tot_clock_time_mulliken = 0.0d0
      tot_clock_time_cube_output = 0.0d0
      tot_clock_time_dielectric = 0.0d0
      tot_clock_time_prodbas_total = 0.d0
      tot_clock_time_ovlp3fn = 0.d0
      tot_clock_time_coulomb_matr = 0.d0
      tot_clock_time_inv_coulomb_matr = 0.d0
      tot_clock_time_ovlp_multi = 0.d0
      tot_clock_time_lvl_triple_recip = 0.d0
      tot_clock_time_first_order_density = 0.d0
      tot_clock_time_first_order_potential = 0.d0
      tot_clock_time_first_order_H = 0.d0
      tot_clock_time_first_order_DM = 0.d0
      tot_clock_time_Sternheimer = 0.d0
      tot_clock_time_first_order_S = 0.d0
      tot_clock_time_Hessian = 0.d0
      tot_clock_time_bc_init = 0.d0
      tot_clock_time_gpu_init = 0.0d0
      tot_clock_time_fodft = 0.d0
      tot_clock_time_matrix_io = 0.d0

      total_number_of_loops           = 0
      relaxation_step_number          = 0
      MD_stepcount                    = 0
      MD_force_evaluations            = 0
      PIMD_stepcount                  = 0
      PIMD_force_evaluations          = 0
      PIMD_beadcount                  = 1
      time_initial_max_scf_cycle      = .false.
      time_initial_max_postprocessing = .false.
      current_atom_num_bsse           = 0

      number_of_scf_inits = 0

      warn_integrals = .false.
      warn_idle_cpus = .false.
      warn_slow_scf_convergence = .false.
      warn_mixing_parameter_range = .false.
      warn_cube_file_format = .false.
      warn_RI_V_grids = .false.
      warn_cpu_consistency = .false.
      warn_surface_dipole = .false.
      warn_mpi_in_place = .false.
      warn_qpe_nonconvergence = .false.

      gpu_density_always_used = .true.
      gpu_hamiltonian_always_used = .true.
      gpu_forces_always_used = .true.

   end subroutine initial_timings
!******

!****s* timing/get_total_time_elapsed
!  NAME
!    get_total_time_elapsed
!  SYNOPSIS
        subroutine get_total_time_elapsed(time_elapsed, clock_time_elapsed, &
                                          cdate, ctime)
!  PURPOSE
!    Calculate the total time that has elapsed since the calculation started.
!    Optionally, returns the current date and time in the format of the
!    date_and_time() built-in.
!  USES
          use synchronize_mpi_basic, only: sync_timing
          implicit none
!  ARGUMENTS
          real*8, intent(out) :: time_elapsed
          real*8, intent(out) :: clock_time_elapsed
          character(len=8),  optional, intent(out) :: cdate
          character(len=10), optional, intent(out) :: ctime
!  INPUTS
!    None
!  OUTPUTS
!    o time_elapsed - CPU time elapsed since beginning of calculation
!    o clock_time_elapsed - Clock time elapsed since beginning of calculation
!    o cdate - (optional) The current date, in date_and_time() format
!    o ctime - (optional) The current time, in date_and_time() format
!  AUTHORS
!    William Huhn (Duke University)
!  HISTORY
!    May 2018 - Added.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject
!    the terms and conditions of the respective license agreement."
!  SOURCE
          character(len=8)  :: cdate_temp
          character(len=10) :: ctime_temp

          character(*), parameter :: func = 'get_total_time_elapsed'

          call cpu_time(rtime)
          time_elapsed = rtime - time_zero

          call date_and_time(cdate_temp, ctime_temp)
          call time_convert (cdate_temp, ctime_temp, clock_time_elapsed)
          clock_time_elapsed = clock_time_elapsed - clock_time_zero

          call sync_timing(time_elapsed)

          if (present(cdate)) write(cdate, "(A8)")  cdate_temp
          if (present(ctime)) write(ctime, "(A10)") ctime_temp

        end subroutine get_total_time_elapsed
!******

!****s* timing/final_timings_and_memory
!  NAME
!    final_timings_and_memory
!  SYNOPSIS
        subroutine final_timings_and_memory(is_a_nice_day)
!  USES
          use dimensions
          use runtime_choices
          use aims_memory_tracking, only : aims_mem_final_output
          use localorb_io, only: localorb_info, SEPARATORLINE, use_unit
          use synchronize_mpi_basic
!  PURPOSE
!    timing calculation at the end of the run
!  ARGUMENTS
          logical, intent(IN) :: is_a_nice_day
!  INPUT
!    is_a_nice_day -- Output "Have a nice day!" to signal success?
!  OUTPUT
!  none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          character(len=8)  :: cdate
          character(len=10) :: ctime
          character(*), parameter :: deffmt = "10X"

          ! get current timestamps
          call get_total_time_elapsed(time_total, clock_time_total, cdate, ctime)

          ! If we have any final warnings, issue them here.
          call final_warnings

          ! synchronize all cumulative CPU-time timestamps
          call sync_timing(tot_time_gpu_init)
          call sync_timing(tot_time_partition)
          call sync_timing(tot_time_free_offset)
          call sync_timing(tot_time_superpos)
          call sync_timing(tot_time_on_evs)
          call sync_timing(tot_time_intg)
          call sync_timing(tot_time_fock)
          call sync_timing(tot_time_diag)
          call sync_timing(tot_time_density)
          call sync_timing(tot_time_mixing)
          call sync_timing(tot_time_hartree_multi)
          call sync_timing(tot_time_hartree_sum)
          call sync_timing(tot_time_etot)
          call sync_timing(tot_time_scaled_zora)
          call sync_timing(tot_time_wf_extra_in)
          call sync_timing(tot_time_wf_extra_out)
          call sync_timing(tot_time_embedding)
          call sync_timing(tot_time_soc)
          call sync_timing(tot_time_band_dos)
          call sync_timing(tot_time_mulliken)
          call sync_timing(tot_time_cube_output)
          call sync_timing(tot_time_dielectric)
          call sync_timing(tot_time_fodft)
          call sync_timing(time_magnetic_response)
          call sync_timing(tot_time_matrix_io)

          ! use_DFPT_polar_reduce_memory
          call sync_timing(tot_time_first_order_density)
          call sync_timing(tot_time_first_order_potential)
          call sync_timing(tot_time_first_order_H)
          call sync_timing(tot_time_first_order_DM)
          call sync_timing(tot_time_Sternheimer)

          write(info_str,'(A)') ''
          call localorb_info( info_str, use_unit )
          write(info_str,'(A)') SEPARATORLINE
          call localorb_info( info_str, use_unit )
          write(info_str,'(10X,A)') &
          "Leaving FHI-aims."
          call localorb_info( info_str, use_unit )

          write(info_str,'(10X,A,A,A,A)') "Date     :  ", cdate, &
             ", Time     :  ", ctime
          call localorb_info( info_str, use_unit )
          write(info_str,'(A)') ''
          call localorb_info( info_str, use_unit )

          write(info_str,'(10X,A)') &
          "Computational steps: "
          call localorb_info( info_str, use_unit )
          write(info_str,'(10X,A,I12)') &
          "| Number of self-consistency cycles          : ",&
          total_number_of_loops
          call localorb_info( info_str, use_unit )
          write(info_str,'(10X,A,I12)') &
          "| Number of SCF (re)initializations          : ",&
          number_of_scf_inits
          call localorb_info( info_str, use_unit )
          if (use_geo_relaxation) then
            write(info_str,'(10X,A,I12)') &
            "| Number of relaxation steps                 : ",&
            relaxation_step_number
            call localorb_info( info_str, use_unit )
          end if
          if (calculate_atom_bsse) then
            write(info_str,'(10X,A,I12)') &
            "| Number of per-atom bsse calculation steps  : ",&
            current_atom_num_bsse
            call localorb_info( info_str, use_unit )
          end if
          if (use_molecular_dynamics) then
            write(info_str,'(10X,A,I12)') &
            "| Number of molecular dynamics steps         : ", MD_stepcount
            call localorb_info( info_str, use_unit )
            write(info_str,'(10X,A,I12)') &
            "| Number of force evaluations                : ",&
            MD_force_evaluations
            call localorb_info( info_str, use_unit )
          end if
          write(info_str,'(A)') ''
          call localorb_info( info_str, use_unit )

          call output_timeheader(deffmt, 'Detailed time accounting')
          call output_times(deffmt, 'Total time', time_total, clock_time_total)
          if (use_gpu) then
            call output_times(deffmt, 'GPU initialization', tot_time_gpu_init,&
                              tot_clock_time_gpu_init)
          end if
          call output_times(deffmt, 'Preparation time', time_prep,&
                            clock_time_prep)
          call output_times(deffmt, 'Boundary condition initalization', &
          &                 tot_time_bc_init, tot_clock_time_bc_init)
          call output_times(deffmt, 'Grid partitioning', &
          &                 tot_time_partition, tot_clock_time_partition)
          call output_times(deffmt, 'Preloading free-atom quantities on grid', &
          &                 tot_time_free_offset, tot_clock_time_free_offset)
          call output_times(deffmt, 'Free-atom superposition energy', &
          &                 tot_time_superpos,  tot_clock_time_superpos)
          call output_times(deffmt, 'Total time for integrations', &
          &                 tot_time_intg, tot_clock_time_intg)
          if(use_prodbas) then
             call total_ovlp3fn_timings(deffmt)
             if(use_hartree_fock) then
                call output_times(deffmt,&
                                  'Total time evaluating exchange matrix', &
                &                 tot_time_fock, tot_clock_time_fock)
             endif
             if(use_ovlp_swap) then
                call output_times(deffmt,&
                                  'Total time swapping 3-center integrals', &
                &                 time_ovlp_swap, clock_time_ovlp_swap)
             endif
          endif
          if(use_full_spectrum) then
             call output_times(deffmt, 'Transforming ovlp_3fn to ovlp_3KS', &
             &                 time_ovlp3fn_to_ovlp3KS, &
             &                 clock_time_ovlp3fn_to_ovlp3KS)
          endif
          if(use_gw) then
             call output_times(deffmt, 'Total time for GW calculation', &
             &                 time_total_gw,clock_time_total_gw)
             call output_times(deffmt, 'GW self-energy (analytic cont.)', &
             &                 time_self_energy, clock_time_self_energy)
             call output_times(deffmt, 'Screened Coulomb matrix W_mn(iw)', &
             &                 time_Wmn_imagfreq, clock_time_Wmn_imagfreq)
             call output_times(deffmt, 'Polarisability matrix Pi_PQ(iw)', &
             &                 time_polar_imagfreq, time_polar_imagfreq)
             if(out_band) then
               call output_times(deffmt, &
                                 'Total time for GW band structure calc.', &
             &                   time_band_self_energy, &
             &                   clock_time_band_self_energy)
             endif
             if(use_gw2ox) then
               call output_times(deffmt, 'Total time for 2OX self-energy', &
               &                 time_2ox_selfenergy, clock_time_2ox_selfenergy)
             elseif(use_gws2ox) then
               call output_times(deffmt, 'Total time for 2nd-order screened exchange self-energy', &
               &                 time_s2ox_selfenergy, clock_time_s2ox_selfenergy)
             elseif(use_gwsoxw) then
               call output_times(deffmt, 'Total time for 2nd-order exchange in W self-energy', &
               &                 time_soxw_selfenergy, clock_time_soxw_selfenergy)
             endif
             if(use_contour_def_gw) then
               call output_times(deffmt, 'Contour def: | Solving qp equations', &
               &                 time_qp_energy_cd, clock_time_qp_energy_cd)
               call output_times(deffmt, 'Contour def: | GW self-energy', &
               &                 time_self_energy_cd, clock_time_self_energy_cd)
               call output_times(deffmt, 'Contour def: | Scr Coulomb matrix W_mn(w)', &
               &                 time_WPQ_realfreq_cd, clock_time_WPQ_realfreq_cd)
               call output_times(deffmt, 'Contour def: | Scr aux Coulomb mat W_PQ(w)', &
               &                 time_WPQ_realfreq_cd, clock_time_WPQ_realfreq_cd)
               call output_times(deffmt, 'Contour def: | Polarisability Pi_PQ(w)', &
               &                 time_WPQ_realfreq_cd, clock_time_WPQ_realfreq_cd)
               if(flag_print_self_energy) then
                 call output_times(deffmt, 'Contour def: | Printing self-energy CD', &
                 &                 time_print_sigma_cd, clock_time_print_sigma_cd)
               endif
               if(flag_calc_spectral_func) then
                 call output_times(deffmt, 'Contour def: | Spectral function CD', &
                 &                 time_spec_func_cd, clock_time_spec_func_cd)
               endif
             endif
          endif
          if(use_rpa_ene .or. use_os_rpa_qpe) then
             call output_times(deffmt, 'Total time for evaluating polarisability', &
             &                 time_polar, clock_time_polar)
             call output_times(deffmt, 'Total time for RPA corr. energy', &
             &                 time_rpa_corr, clock_time_rpa_corr)
             if(n_periodic .gt. 0) then
               call output_times(deffmt, 'Total time for rSE corr. energy', &
               &                 time_rse_corr, clock_time_rse_corr)
             endif
          endif
          if(use_rpa_plus_2ox) then
             call output_times(deffmt, 'Total time for calculating 2OX corr. energy', &
             &                 time_2oex_corr,  clock_time_2oex_corr)
          endif
          if(use_rpa_plus_sosex) then
             call output_times(deffmt, 'Total time for calculating SOSEX corr. energy', &
             &                 time_sosex_corr,  clock_time_sosex_corr)
          endif
          if(use_C6_coef) then
             call output_times(deffmt, 'Total time for calculating C6 coefficients', &
             &                 time_C6_coef, clock_time_C6_coef)
          endif
          call output_times(deffmt, 'Total time for solution of K.-S. equations', &
          &                 tot_time_diag,  tot_clock_time_diag)
          if (orthonormalize_evs .and. &
          &   (use_geo_relaxation.or.use_molecular_dynamics)) then
             call output_times(deffmt, 'Total time for EV reorthonormalization', &
             &                 tot_time_on_evs,  tot_clock_time_on_evs)
          end if
          if (force_potential.eq.0) then
             if (.not.use_forces) then
                call output_times(deffmt, 'Total time for density update', &
                &                 tot_time_density,  tot_clock_time_density)
             else
                call output_times(deffmt, 'Total time for density & force components',&
                &                 tot_time_density, tot_clock_time_density)
             end if
             if (use_kerker_preconditioner) then
                call output_times(deffmt, 'Total time for mixing & preconditioning',&
                &                 tot_time_mixing, tot_clock_time_mixing)
             else
                call output_times(deffmt, 'Total time for mixing', &
                &                 tot_time_mixing, tot_clock_time_mixing)
             end if
             call output_times(deffmt, 'Total time for Hartree multipole update', &
             &            tot_time_hartree_multi, tot_clock_time_hartree_multi)
             call output_times(deffmt, 'Total time for Hartree multipole sum', &
             &            tot_time_hartree_sum, tot_clock_time_hartree_sum)
             call output_times(deffmt, 'Total time for total energy evaluation', &
             &                 tot_time_etot, tot_clock_time_etot)
          end if
          if (flag_rel.gt.0) then
             call output_times(deffmt, 'Total time for scaled ZORA corrections', &
             &            tot_time_scaled_zora, tot_clock_time_scaled_zora)
          endif
          if (use_vdw_correction_hirshfeld .or. use_mbd_old &
            .or. use_mbd_dev .or. use_mbd_std .or. use_libmbd) &
          then
             call output_times(deffmt, 'Total time for vdW correction', &
             &            tot_time_vdw, tot_clock_time_vdw)
          endif
          if (use_wf_extrapolation) then
             call output_times(deffmt, 'Total time for extrapolation saves', &
             &            tot_time_wf_extra_in, tot_clock_time_wf_extra_in)
             call output_times(deffmt, 'Total time for extrapolations',&
             &            tot_time_wf_extra_out, tot_clock_time_wf_extra_out)
          endif

          if (use_embedding_pp.or.use_embedding_potential) then
             call output_times(deffmt, 'Total time for QM/MM embedding', &
             &            tot_time_embedding, tot_clock_time_embedding)
          endif

          if (magnetic_response) then
             call output_times(deffmt, 'Total time for magnetic response', &
                  & time_magnetic_response, clock_time_magnetic_response)
          end if

          if(flag_neutral_excitation) then
            call output_times(deffmt, 'Total time for neutral excitations', &
                  time_total_neutral_excitation, clock_time_total_neutral_excitation)
            call output_times(deffmt, '   * dipole moment integration', &
                  time_neutral_excitation_dipole, clock_time_neutral_excitation_dipole)
            call output_times(deffmt, '   * orbital calculations', &
                  time_neutral_excitation_orbitals, clock_time_neutral_excitation_orbitals)
            call output_times(deffmt, '      + <mu|f_xc|nu>', &
                  time_fxc_integr, clock_time_fxc_integr)
            call output_times(deffmt, '      + compute D_ia^mu', &
                  time_coeffs_fxc3fn, clock_time_coeffs_fxc3fn)
            call output_times(deffmt, '      + KS basis transformation', &
                  time_ovlp3fn_to_ovlp3KS, clock_time_ovlp3fn_to_ovlp3KS)
            call output_times(deffmt, '      + <mu|f_xc|ia>', &
                  time_matvec_cont, clock_time_matvec_cont)
            call output_times(deffmt, '   * casida formalism', &
                  time_neutral_excitation_casida, clock_time_neutral_excitation_casida)
            call output_times(deffmt, '      + TDDFT Matrix building', &
                  time_omegabuild, clock_time_omegabuild)
            call output_times(deffmt, '         - <mu|f_xc|ia> communication', &
                  time_3ks_comm, clock_time_3ks_comm)
            call output_times(deffmt, '         - casida matrix calculation', &
                  time_casida_matrix_calculation, clock_time_casida_matrix_calculation)
            call output_times(deffmt, '         - casida elements communication', &
                  time_omegaresdist, clock_time_omegaresdist)
            call output_times(deffmt, '      + solving the TDDFT Matrix', &
                  time_omega_solve, clock_time_omega_solve)
          endif

          if (neutral_excitation_bse) then
            call output_times(deffmt, 'Total time for BSE', &
                  time_total_bse + time_ovlp3fn, &
                  clock_time_total_bse + clock_time_ovlp3fn)
            call output_times(deffmt, '   * orbital calculations/basis setup', &
                  time_ovlp3fn + time_bse_ovlp_3ks + &
                  time_bse_w0 + time_bse_w_ovlp_3ks, &
                  clock_time_ovlp3fn + clock_time_bse_ovlp_3ks + &
                  clock_time_bse_w0 + clock_time_bse_w_ovlp_3ks)
            call output_times(deffmt, '      + ovlp_3fn', &
                 time_ovlp3fn, clock_time_ovlp3fn)
            call output_times(deffmt, '      + ovlp_3ks', &
                 time_bse_ovlp_3ks, clock_time_bse_ovlp_3ks)
            call output_times(deffmt, '      + W0 bas', &
                 time_bse_w0, clock_time_bse_w0)
            call output_times(deffmt, '      + W * ovlp_3ks', &
                  time_bse_w_ovlp_3ks, clock_time_bse_w_ovlp_3ks)
            call output_times(deffmt, '   * BSE Matrix building', &
                  time_bse_construct_mat - time_bse_w0 - time_bse_w_ovlp_3ks, &
                  clock_time_bse_construct_mat - clock_time_bse_w0 - &
                  clock_time_bse_w_ovlp_3ks)
            call output_times(deffmt, '   * solving the BSE Matrix', &
                  time_bse_solve_mat, clock_time_bse_solve_mat)
            call output_times(deffmt, '   * dipole moment integration', &
                  time_bse_osc, clock_time_bse_osc)
          endif

          if(output_matrix_io_timing) then
             call output_times(deffmt, "Total time for matrix I/O", &
                  tot_time_matrix_io, tot_clock_time_matrix_io)
          end if

          if(calculate_perturbative_soc) then
             call output_times(deffmt, "Total time for perturbative SOC", &
               tot_time_soc, tot_clock_time_soc)
          end if

          if (out_band.or.out_dos) then
                  call output_times(deffmt, 'Total time for band structures, DOS', &
              tot_time_band_dos, tot_clock_time_band_dos)
          end if

          if (flag_run_mulliken) then
                  call output_times(deffmt, 'Total time for Mulliken analysis', &
              tot_time_mulliken, tot_clock_time_mulliken)
          end if

          if (out_cube) then
                  call output_times(deffmt, 'Total time for cube output', &
              tot_time_cube_output, tot_clock_time_cube_output)
          end if

          if(flag_out_dielectric .or. flag_out_absorption) then
                  call output_times(deffmt, 'Total time for dielectric calculations', &
              tot_time_dielectric, tot_clock_time_dielectric)
          end if


          if(use_DFPT_reduce_memory) then
             call output_times(deffmt,&
                               'Total time firt_order_denstiy', &
             &                 tot_time_first_order_density, tot_clock_time_first_order_density)
             call output_times(deffmt,&
                               'Total time firt_order_potential', &
             &                 tot_time_first_order_potential, tot_clock_time_first_order_potential)
             call output_times(deffmt,&
                               'Total time firt_order_H', &
             &                 tot_time_first_order_H, tot_clock_time_first_order_H)
             call output_times(deffmt,&
                               'Total time solving of Sternheimer equations in DFPT', &
             &                 tot_time_Sternheimer, tot_clock_time_Sternheimer)
             call output_times(deffmt,&
                               'Total time firt_order_S', &
             &                 tot_time_first_order_S, tot_clock_time_first_order_S)
             call output_times(deffmt,&
                               'Total time for Hessian', &
             &                 tot_time_Hessian, tot_clock_time_Hessian)
          endif

          if(fo_finalise) then
             call output_times(deffmt, 'Total time for FO-DFT calculation', &
                               tot_time_fodft, tot_clock_time_fodft)
          endif

          if(use_DFPT_polar_reduce_memory) then
             call output_times(deffmt,&
                               'Total time firt_order_denstiy', &
             &                 tot_time_first_order_density, tot_clock_time_first_order_density)
             call output_times(deffmt,&
                               'Total time firt_order_DM', &
             &                 tot_time_first_order_DM, tot_clock_time_first_order_DM)
             call output_times(deffmt,&
                               'Total time firt_order_potential', &
             &                 tot_time_first_order_potential, tot_clock_time_first_order_potential)
             call output_times(deffmt,&
                               'Total time firt_order_H', &
             &                 tot_time_first_order_H, tot_clock_time_first_order_H)
             call output_times(deffmt,&
                               'Total time solving of Sternheimer equations in DFPT', &
             &                 tot_time_Sternheimer, tot_clock_time_Sternheimer)
          endif


          if (use_gpu) then
            write(info_str,'(10X,A)') ""
            call localorb_info( info_str, use_unit )
            if (use_gpu_density) then
              if (gpu_density_always_used) then
                write(info_str,'(10X,A)') "GPU acceleration has been used for charge density update."
                call localorb_info( info_str, use_unit )
              else
                write(info_str,'(10X,A)') "GPU acceleration was requested for charge density update, but &
                                          &was not always used."
                call localorb_info( info_str, use_unit )
              end if
            end if
            if (use_gpu_hamiltonian) then
              if (gpu_hamiltonian_always_used) then
                write(info_str,'(10X,A)') "GPU acceleration has been used for real-space Hamiltonian integration."
                call localorb_info( info_str, use_unit )
              else
                write(info_str,'(10X,A)') "GPU acceleration was requested for real-space Hamiltonian integration, but &
                                          &was not always used."
                call localorb_info( info_str, use_unit )
              end if
            end if
            if (use_forces.and.use_gpu_forces) then
              if (gpu_forces_always_used) then
                write(info_str,'(10X,A)') "GPU acceleration has been used for calculating contributions to forces."
                call localorb_info( info_str, use_unit )
              else
                write(info_str,'(10X,A)') "GPU acceleration was requested for calculating contributions to forces, but &
                                          &was not always used."
                call localorb_info( info_str, use_unit )
              end if
            end if
            write(info_str,'(10X,A)') "Due to GPU usage, wall_clock(cpu1) will give a more reliable estimate for timings."
            call localorb_info( info_str, use_unit )
          end if

          write(info_str,'(A)') ''
          call localorb_info( info_str, use_unit )
          call aims_mem_final_output ()

          write(info_str,'(A)') ''
          call localorb_info( info_str, use_unit )
          if (is_a_nice_day) then
             write(info_str,'(10X,A)') "Have a nice day."
             call localorb_info( info_str, use_unit )
          end if
          write(info_str,'(A)') &
          "------------------------------------------------------------"
          call localorb_info( info_str, use_unit )

        end subroutine final_timings_and_memory
!******
  !----------------------------------------------------------------------------
  !****s* timing/ovlp3fn_timings
  !  NAME
  !    ovlp3fn_timings
  !  SYNOPSIS

  subroutine ovlp3fn_timings(fmt)

    !  PURPOSE
    !
    !    Output timings related to initialize HF or correlated methods, i.e.,
    !    to calculate ovlp_3fn/(coeff_3fn_ten, coulomb_matr_lvl).
    !
    !  USES

    use dimensions, only: use_hf_realspace, use_lvl
    use runtime_choices, only: use_logsbt_lvltriples, RI_type, RI_LVL
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: fmt

    !  INPUTS
    !    o fmt -- Format string to add (e.g. '2X')
    !  OUTPUTS
    !    none [-> STDOUT]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'ovlp3fn_timings'

    call output_times(fmt, "Product basis setup: Total time", &
    &                 time_prodbas_total, clock_time_prodbas_total)

    if (.not.use_hf_realspace) then
       ! more detailed subtimings are available in this case

       if (.not. use_lvl .or. .not. use_logsbt_lvltriples) then
          call output_times(fmt, "Product basis: | 3-center integrations", &
          &                 time_ovlp3fn, clock_time_ovlp3fn)
       end if
       call output_times(fmt, "Product basis: | 2-center integrations", &
       &                 time_coulomb_matr, clock_time_coulomb_matr)
       if (RI_type /= RI_LVL) then
          call output_times(fmt, "Product basis: | 2-center linear algebra", &
          &                 time_inv_coulomb_matr, clock_time_inv_coulomb_matr)
       end if
       if (use_lvl .and. use_logsbt_lvltriples) then
          call output_times(fmt, "Product basis: | LVL 3-center & coefficients", &
          &                 time_ovlp_multi, clock_time_ovlp_multi)
       else
          call output_times(fmt, "Product basis: | 3-center x 2-center", &
         &                 time_ovlp_multi, clock_time_ovlp_multi)
       end if

    end if


  end subroutine ovlp3fn_timings
  !******
  !----------------------------------------------------------------------------
  !****s* timing/ovlp3fn_timings
  !  NAME
  !    ovlp3fn_timings
  !  SYNOPSIS

  subroutine total_ovlp3fn_timings(fmt)

    !  PURPOSE
    !
    !    Output FINAL timings related to initialize HF or correlated methods, i.e.,
    !    to calculate ovlp_3fn/(coeff_3fn_ten, coulomb_matr_lvl).
    !
    !  USES

    use dimensions, only: n_periodic, use_hf_kspace, use_rpa_ene, use_gw, &
        use_lvl, use_hf_realspace
    use runtime_choices, only: use_logsbt_lvltriples, RI_type, RI_LVL
    implicit none

    !  ARGUMENTS

    character(*), intent(IN) :: fmt

    !  INPUTS
    !    o fmt -- Format string to add (e.g. '2X')
    !  OUTPUTS
    !    none [-> STDOUT]
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE

    character(*), parameter :: func = 'ovlp3fn_timings'

    call output_times(fmt, "Product basis setup: Total time", &
    &                 tot_time_prodbas_total, tot_clock_time_prodbas_total)

    if (.not.use_hf_realspace) then
       ! more detailed subtimings are available in this case

       if (.not. use_lvl .or. .not. use_logsbt_lvltriples) then
          call output_redundant_times(fmt, "Product basis: | 3-center integrations", &
          &                 tot_time_ovlp3fn, tot_clock_time_ovlp3fn)
       end if
       call output_redundant_times(fmt, "Product basis: | 2-center integrations", &
       &                 tot_time_coulomb_matr, tot_clock_time_coulomb_matr)
       if (RI_type /= RI_LVL) then
          call output_redundant_times(fmt, "Product basis: | 2-center linear algebra", &
          &                 tot_time_inv_coulomb_matr, tot_clock_time_inv_coulomb_matr)
       end if
       if (use_lvl .and. use_logsbt_lvltriples) then
          call output_redundant_times(fmt, "Product basis: | LVL 3-center & coefficients", &
          &                 tot_time_ovlp_multi, tot_clock_time_ovlp_multi)
       else
          call output_redundant_times(fmt, "Product basis: | 3-center x 2-center", &
         &                 tot_time_ovlp_multi, tot_clock_time_ovlp_multi)
       end if

    end if

    if ((n_periodic .gt. 0) .and. (use_hf_kspace .or. use_rpa_ene .or. use_gw)) then
         call output_times(fmt, "Transforming LVL triples to recip. space", &
         &                 tot_time_lvl_triple_recip, tot_clock_time_lvl_triple_recip)
    endif

  end subroutine total_ovlp3fn_timings
  !******
!****s* timing/update_time_max_scf_cycle
!  NAME
!    update_time_max_scf_cycle
!  SYNOPSIS
        subroutine update_time_max_scf_cycle(t_this_cycle)
!  PURPOSE
!    figure out the time for the longest scf cycle so far
!  ARGUMENTS

          real*8 ::  t_this_cycle
!  INPUTS
!   o t_this_cycle -- time for the current scf cycle
!  OUTPUT
!   none
!
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          if (.not.time_initial_max_scf_cycle) then
             time_max_scf_cycle = t_this_cycle
             time_initial_max_scf_cycle = .true.
          else
             if (t_this_cycle.gt.time_max_scf_cycle) time_max_scf_cycle = t_this_cycle
          end if
        end subroutine update_time_max_scf_cycle
!******

!****s* timing/update_time_max_postprocessing
!  NAME
!    update_time_max_postprocessing
!  SYNOPSIS
        subroutine update_time_max_postprocessing(t_this_cycle)
!  PURPOSE
!    figure out the time for the longest postprocessing time so far
!
!  ARGUMENTS

          real*8 :: t_this_cycle

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    o t_this_cycle -- time for the current postprocessing cycle
!
!  OUTPUT
!  none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          if (.not.time_initial_max_postprocessing) then
             time_max_postprocessing = t_this_cycle
             time_initial_max_postprocessing = .true.
          else
             if (t_this_cycle.gt.time_max_postprocessing) time_max_postprocessing = t_this_cycle
          end if
        end subroutine update_time_max_postprocessing
!******

!****f* timing/time_check_scf_cycle
!  NAME
!    time_check_scf_cycle
!  SYNOPSIS
        logical function time_check_scf_cycle()
! PURPOSE
!    to check if there should be enough time for another scf cycle before
!    the walltime runs out.
! USES
          use runtime_choices, only: walltime_requested, walltime_limit
!  INPUTS
!   none
!  OUTPUTS
!   o time_check_scf_cycle -- enough time?
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          logical check
          ! was there a previous reading on scf cycle speed or do we have a walltime limit at all ????
          check = (.not.time_initial_max_scf_cycle).or.(.not.walltime_requested)
          ! test whether or not the time limit has been reached!
          if (.not.check) then
             call cpu_time(rtime)
             check = ((rtime+time_max_scf_cycle-time_zero).lt.walltime_limit)
          end if
          time_check_scf_cycle = check
          return
        end function time_check_scf_cycle
!******

!****f* timing/time_check_postprocessing
!  NAME
!    time_check_postprocessing
!  SYNOPSIS
        logical function time_check_postprocessing()
!  PURPOSE
!    returns true if there should be enough time for the scf postprocessing in scf_solver.f90
!  USES
     use runtime_choices, only: walltime_requested, walltime_limit
!  INPUTS
!    none
!  OUTPUTS
!    o time_check_postprocessing -- enough time for postprocessing?
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
          logical check
          check = (.not.time_initial_max_postprocessing).or.(.not.walltime_requested)
          if (.not.check) then
             call cpu_time(rtime)
             check = ((rtime+time_max_postprocessing-time_zero).lt.walltime_limit)
          end if
          time_check_postprocessing = check
          return
        end function time_check_postprocessing
!******
  !----------------------------------------------------------------------------
  !****s* timing/final_warnings
  !  NAME
  !    final_warnings
  !  SYNOPSIS

  subroutine final_warnings

    !  PURPOSE
    !
    !  Sometimes something unpleasant happens during the run.
    !  A warning is written out in the log file (but the code does not stop).
    !  Repeat some of these warnings here. (NO new warnings here, please!)
    !
    !  USES

    use applicable_citations, only: warn_citation_unknown
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks, only: n_tasks
    implicit none

    !  ARGUMENTS

    !  INPUTS
    !    none
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    if (warn_integrals) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* A warning was issued further up regarding the accuracy of onsite integrals.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  An integral with a negative eigenvalue was among the flagged integrals.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  In most standard situations, this warning is probably harmless.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  If needed, this can be verified by adjusting the radial_multiplier flag.'
       call localorb_info(info_str,use_unit,'(A)')
       end if

    if (warn_idle_cpus) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         '* A warning was issued further up regarding the CPU load.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  Some CPUs do not participate in the eigenvalue solution.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         '  In corner cases, this may be the correct behavior. Otherwise,'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         '  consider using a scalapack/elpa version of the code for a parallel treatment.'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (warn_slow_scf_convergence) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         '* A warning was issued further up regarding the speed of the s.c.f. convergence.'
       write (info_str,'(1X,A)') &
         '  One or more s.c.f. cycles took more than 50 s.c.f. iterations.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  It might be worth adjusting your convergence settings. For example, '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         '  metals need a larger broadening and a smaller "charge_mix_param" value than semiconductors.'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (warn_mixing_parameter_range) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* Warning: The chosen 'charge_mix_param' value was not between zero or one."
       write (info_str,'(1X,A)') &
         '* The value found in control.in should have caused s.c.f. convergence problems.'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (warn_cube_file_format) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* A warning regarding the cube file format (cube_content_unit legacy)"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* was issued above. The output is Angstrom-based but the cube voxel'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* definition is in atomic units. Please check carefully.'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (warn_citation_unknown) then

       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(1X,A)') &
        "* Warning - a code part requested a citation which was not known"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* to subroutine cite_reference in module applicable_citations.f90 ."
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* Likely this is a programmer's oversight. Please let us know."
       call localorb_info(info_str)

    end if

    if (warn_RI_V_grids) then

       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(1X,A)') &
        "* Warning - angular integration grids had to be adjusted for use with RI_method 'V'."
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* Background: You are using an electronic structure method that requires RI_method 'V'"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* for the resolution of identity of the Coulomb operator. (See Ren et al.,"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* New J. Phys 2012 for details.) This method requires three-center integrations"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* to be carried out on FHI-aims' normal angular and radial integration grids."
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* (See FHI-aims CPC Blum et al. 2009, p. 2183 for grid details and references.)"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* One of the angular grids specified in your control.in file was not large"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* enough for the maximal angular momenta of orbital and product basis functions"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* for its species. "
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* This can lead to atomic total energy errors which cancel in energy differences."
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* Effective Sep 22, 2014, we have changed the default behavior of the code to"
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* simply 'upgrade' the angular grids internally to a safe setting in this case."
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
        "* Please be sure to NOT mix total energies from before and after the change."
       call localorb_info(info_str)

    end if

    if (warn_cpu_consistency) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* During the run, the code found that basic geometry arrays were"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* not always the same on different MPI tasks. Please see the'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* information printed further up in this output file.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* We have observed that flips of individual bits can happen for'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* arrays that should formally be the same, but this happens only for'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* certain compilers and / or hardware platforms. Usually, this is'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* not problematic but it would be good to check that this is indeed'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '* all that is going on.'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if ( n_tasks .gt. 1000 ) then
      write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '*********************************************************************'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  Your are running FHI-aims on many cores. To achieve optimal'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  performance for your particular architecture, it might be necessary'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  to tune the following compilation and calculation settings:'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '   * Choosing the right Scalapack and MPI Libraries for your'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '      architecture is essential. If in doubt, refer to the'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '      documentation or the sysadmin of your supercomputer.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '   * Please review our documentation of Makefile choices, particularly'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '      for ELPA, in the manual.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '   * Some of the following FHI-aims keywords or combination thereof'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '      may be helpful:'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '       - FHI-aims setting: collect_eigenvectors .false. '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '       - FHI-aims setting: use_local_index .true.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '       - FHI-aims setting: load_balancing .true.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '       - FHI-aims setting: use_alltoall .false.'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '       - FHI-aims setting: elpa_settings two_step_solver'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  Be mindful that the usefulness of these keywords can depend on your'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  actual run and on the actual computer you are using. Please check'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  their description in the manual and perform specific benchmarks to'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  determine what is optimal for your computer. More detailed'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '  information and examples can be found at:'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         'http://aimsclub.fhi-berlin.mpg.de/wiki/index.php/Large_scale_settings'
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,A)') &
         '*********************************************************************'
       call localorb_info(info_str,use_unit,'(A)')
    end if

    if (warn_surface_dipole) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* A warning was triggered above, indicating that a surface dipole was"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* determined in a vacuum region that was not large enough."
       call localorb_info(info_str,use_unit,'(A)')
       write(info_str,'(1X,A)') &
         '* A surface dipole determined close to the atoms of a slab will not be accurate.'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* As a result, any numerical values that are affected by the dipole (including'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* total energies) may be subject to noticeable inaccuracies.'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* The calculation was allowed to continue, but we strongly recommend'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* to repeat it with a larger vacuum layer between slabs, large enough to make'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* this warning disappear.'
       call localorb_info(info_str)

    end if

    if (warn_mpi_in_place) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* A warning was triggered above, indicating that the MPI library used in your run"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* does not support MPI_IN_PLACE. We conducted the run without using this feature,"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* leading to a somewhat less memory efficient run."
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* However, check the results carefully, since the observed problems with MPI_IN_PLACE"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* could mean that your MPI library is outdated or otherwise problematic."
       call localorb_info(info_str,use_unit,'(A)')

    end if

    if (warn_qpe_nonconvergence) then
       write (info_str,'(1X,A)') &
         ' '
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* Warning! One or more quasiparticle eigenvalues listed above could not be"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* determined unambiguously in the list given above. This is related to the"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* default numerical fit procedure used to determine parameters for analytical"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* continuation of the self-energy from the imaginary axis to the real axis."
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* More information can be found in Sec. 3.21 in the manual."
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* In particular, the keyword 'anacon_type 0' can be used to change the default"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* analytical continuation type to the slightly less accurate 'two-pole' self-"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* energy model, which is more robust than the 16-parameter Pade approximation"
       call localorb_info(info_str,use_unit,'(A)')
       write (info_str,'(1X,1A)') &
         "* used by default. Again, please see the manual for the mathematical definitions."
       call localorb_info(info_str,use_unit,'(A)')

    end if

  end subroutine final_warnings
  !******
      end module timing
