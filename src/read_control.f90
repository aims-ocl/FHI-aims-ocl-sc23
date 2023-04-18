!****s* FHI-aims/read_control
!  NAME
!    read_control
!  SYNOPSIS

subroutine read_control ( )

!  PURPOSE
!    This is the central subroutine that reads and checks all information found
!    in input file "control.in" . The centerpiece is simply a line-by-line reading
!    of control.in, with a large "if" statement deciding what to do if a given
!    keyword is found at the beginning of a line.
!
!    Note that each keyword found in control.in must also be listed in module
!    dimensions.f90, subroutine parse_control . That subroutine is necessary
!    in order to determine certain array dimensions for correct allocation, before
!    the actual information in those arrays is read by read_control ( ).
!
!    Most of the variables _read_ in the present subroutine are globally accessible
!    variables / options, declared in module runtime_choices.f90 .
!
!    The species subsections of control.in are read and verified separately by subroutine
!    read_species_data.f90 .
!
!    The full verification of _all_ input data should happen _right after_ control.in has
!    been fully read and closed, i.e., as a part of the present subroutine. While this
!    makes for a long subroutine here, keeping all verifications close to the actual
!    reading process is immensely helpful in reducing the complexity of these checks,
!    compared to a separate, case-by-case verification in disconnected locations afterwards.
!  USES

   use types, only: dp
   use dimensions
   use runtime_choices
   use grids
   use species_data
   use plot
   use mpi_tasks
   use localorb_io
   use timing, only: warn_mixing_parameter_range, warn_RI_V_grids, warn_mpi_in_place
   use control_file
   use constraint
   use constants
   use plus_u
   use elpa1_2013, only: elpa_print_times
   use vdw_correction
   use gw_para
   use ll_vdwdf
   use synchronize_mpi
   use force_occupation, only: force_occ_pr_state, force_occ_spin, &
                               forced_occ_number, force_occ_min_KS_state, &
                               force_occ_max_KS_state, &
                               force_occupation_smearing_width, &
                               force_occ_atom, force_occ_basis_type, &
                               force_occ_basis_n, force_occ_basis_l, &
                               force_occ_basis_m, force_occ_basis_occupation, &
                               force_occ_maxred_KS_state, &
                               force_occ_step_KS_state, force_occ_autorep, &
                               force_occ_autored, fop_detect_auto
   ! force_occupation_smearing (in dimensions.f90)
   use transport
   use wf_extrapolation
   use numerical_stress
   use octree_routines
   use pseudodata
   use poles_fit
   use calculate_fock_matrix_p0, only : crit_val, coul_mat_threshold
   use pi_molecular_dynamics   ! XZL: added for PIMD
   use energy_density
   use molecular_dynamics, only: read_gle
   use dmft_para
   use xc_f03_lib_m
   use xc_library
   use esp_charges, only: read_esp_parameters
   use crpa_blacs, only : blacsdim
   use cpt2_blacs, only : blacsdim_cpt2, pair_block_cpt2
   use friction
   use debugmanager, only: activate_debugging
   use xml_write, only: xml_open_file
   use mpe_constants, only: ISC_CONST, MPE_CONST
   use fodft, only: fodft_combine, fodft_set_states
   use mbd_dev_wrapper, only: mbd_dev_parse, mbd_dev_verify, mbd_dev_flags
   use mbd_std_wrapper, only: mbd_std_flags, mbd_std_parse
   use hirshfeld, only: grid_out
   use python_interface, only: register_python_hook
   use applicable_citations
   use dimensions_soc, only : n_core_states_omit_from_soc, &
                              n_high_states_omit_from_soc
   use psi_at_nucleus_mod, only: initialize_psi_at_nucleus
   use boys, only: boys_sub_min_KS_state, boys_sub_max_KS_state, boys_sub_flags, boys_sub_flag
   use aims_memory_tracking, only: aims_mem_debug, aims_mem_sync
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team.
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!  e.V. Please note that any use of the "FHI-aims-Software" is subject to
!  the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
!@@edu>
!flag_xc =
! -1  not set
!  0  hf
!  1  pbe0
!  2
!  3  pz-lda
!  4
!  5  pbe-vdw
!  6  pbe
!  7  hse
!  8  pw-lda
!  9  blyp
! 10  b3-lyp
! 11  rpbe
! 12  revpbe
! 13  pbesol0
! 14  revpbe_vdw
! 15  vwn - this is "vwn5" (the full LDA parametrization by vwn)
! 16  vwn_gauss (this is the RPA paramerization by vwn)
! 17  pbesol
! 18
! 19
! 20  am05
! 21  pbeint
! 22  xpbe
! 23  lc_wPBEh
! 24
! 25  m06-l
! 26  m06-hf
! 27  m06
! 28  m06-2x
! 29  b1lyp
! ...
! reserve [30,50] for advanced orbital-dependent DFAs
! 30  xyg3
! 31  xdh-pbe0
! 32  xygjos
! 33  zrps
! 34  lrc-xyg6
! 35
! ...
! As requested. AJL.
! 51  tpss
! 52  revtpss
! 53  tpssloc
! 54
! 55 m08-hx
! 56 m08-so
! 57 m11
! 58 m11-l
! 59 scan
! ...
! -10 > libxc
! x+xc_dfauto_offset dfauto  (uses the normal values from above increased by
! xc_dfauto_offset)
!@@edu<

   implicit none

   logical :: modified

   logical :: flag_mixer
   logical :: flag_mixer_constraint
   logical :: flag_mix_param
   logical :: flag_prec_mix_param
   logical :: flag_spin_mix_param
   logical :: flag_pulay_iter
   logical :: flag_pulay_iter_constraint
   logical :: flag_broyden_iter
   logical :: flag_relative_fp_charge_mix
   logical :: flag_ini_linear_mixing
   logical :: flag_linear_mix_param
   logical :: flag_linear_spin_mix_param
   logical :: flag_ini_linear_mixing_constraint
   logical :: flag_sc_max
   logical :: flag_sc_init_iter
   logical :: flag_sc_init_factor
   logical :: flag_sc_abandon_etot
   logical :: flag_constraint_iter
   logical :: flag_constraint_mix
   logical :: flag_constraint_precision
   logical :: flag_partition_acc
   logical :: flag_partition_type
   logical :: flag_basis_threshold
   logical :: flag_prodbas_threshold
   logical :: flag_wave_threshold
   logical :: flag_force_lebedev
   logical :: flag_occupation_type
   logical :: flag_DFPT_occupation_type
   logical :: flag_occupation_thr
   logical :: flag_occupation_acc
   logical :: flag_mu_method
   logical :: flag_fermi_acc
   logical :: flag_max_zeroin
   logical :: flag_spin
   logical :: flag_angular_min_size
   logical :: flag_angular_div_cutoff
   logical :: flag_force
   logical :: flag_delta_numerical_stress
   logical :: flag_relax_mode
   logical :: flag_molecular_dynamics
   logical :: flag_MD_schedule
   logical :: flag_wf_named
   logical :: flag_wf_func
   logical :: flag_distributed_spline_storage
   logical :: flag_remove_unitary_force_component
   logical :: flag_energy_tolerance
   logical :: flag_aggregated_energy_tolerance
   logical :: flag_harmonic_length_scale
   logical :: flag_max_atomic_move
   logical :: flag_min_line_step
   logical :: flag_line_step_reduce
   logical :: flag_max_relaxation_steps
   logical :: flag_init_hess
   logical :: flag_distributed_hessian
   logical :: flag_points_in_batch
   logical :: flag_batch_size_limit
   logical :: flag_grid_partitioning_method
   logical :: flag_grouping_factor
   logical :: flag_preconditioner
   logical :: flag_preconditioner_turnoff
   logical :: flag_preconditioner_max_l
   logical :: flag_legacy_monopole_extrapolation
   logical :: flag_packed_matrix_format
   logical :: flag_hartree_worksize
   logical :: flag_k_points_defined
   logical :: flag_symmetry_reduced_k_grid
   logical :: flag_mixer_threshold_charge
   logical :: flag_mixer_threshold_spin
   logical :: flag_min_batch_size
   logical :: flag_use_metis_batching
   logical :: flag_initial_ev_solutions
   logical :: flag_lopcg_preconditioner
   logical :: flag_lopcg_tolerance
   logical :: flag_lopcg_block_size
   logical :: flag_max_cg_iterations
   logical :: flag_atom_ref
   logical :: flag_scsmp2
   logical :: flag_multipole_threshold
   logical :: flag_l_hartree_far_distance
   logical :: flag_density_update_method
   logical :: flag_communication_type
   logical :: flag_found_phonon
!  logical :: flag_set_vacuum_level (temporarily?) moved to runtime_choices
   logical :: flag_transport
   logical :: flag_KS_method
   logical :: flag_mbd_eigensolver
   logical :: flag_hybrid_coeff
   logical :: flag_use_local_index
   logical :: flag_use_local_index_if_scalapack
   logical :: flag_load_balancing
   logical :: flag_load_balancing_if_scalapack
   logical :: flag_use_alltoall
   logical :: flag_MB_clean_rotations_read
   logical :: flag_freq_grid_type
   logical :: flag_anacon_type
   logical :: flag_scgw_no_core
   logical :: flag_frequency_points
   logical :: flag_n_anacon_par
   logical :: flag_maximum_frequency
   logical :: flag_maximum_time
   logical :: flag_time_points
   logical :: flag_MD_QH_file_exists
   logical :: flag_use_logsbt
   logical :: flag_sbtgrid_lnr0
   logical :: flag_sbtgrid_lnk0
   logical :: flag_sbtgrid_lnrange
   logical :: flag_sbtgrid_N
   logical :: flag_use_logsbt_for_radial_hse_integration
   logical :: flag_use_2d_corr
   logical :: flag_calculate_atom_bsse
   logical :: flag_default_moment_defined
   logical :: flag_default_max_n_prodbas
   logical :: flag_default_max_l_prodbas
   logical :: flag_default_prodbas_acc
   logical :: flag_MD_tstep
   logical :: flag_hse_unit
   logical :: flag_hse03
   logical :: flag_crit_val
   logical :: flag_coul_mat_threshold
   logical :: flag_RI
   logical :: flag_restrict_kpoint_to_smp_node
   logical :: flag_max_tasks_per_smp_node
   logical :: flag_multip_moments_threshold
   logical :: flag_write_restart_geometry
   logical :: flag_hessian_to_restart_geometry
   logical :: flag_compensate_multipole_errors
   logical :: flag_normalize_initial_density
   logical :: flag_vdw_convergence_threshold
   logical :: flag_calculate_friction
   logical :: flag_cavity_restart_file_exists

   logical :: flag_elsi_method

   logical :: flag_out_mommat  ! This is only an internal flag for use in read_control.f90 - nowhere else.
                               ! It records whether a line called 'compute_momentummatrix' was found in control.in.

   logical :: flag_collect_eigenvectors ! This is only a flag to record if the keyword 'collect_eigenvectors' was
                                        ! found in control.in. See previous flags - same philosophy everywhere.

   logical, dimension(n_species) :: flag_radial
   logical, dimension(n_species) :: flag_angular
   logical, dimension(n_species) :: flag_angular_min
   logical, dimension(n_species) :: flag_angular_acc
   logical, dimension(n_species) :: flag_cut_free_atom
   logical, dimension(n_species) :: flag_max_n_prodbas
   logical, dimension(n_species) :: flag_max_l_prodbas
   logical, dimension(n_species) :: flag_prodbas_acc
   logical, dimension(n_atoms) :: flag_region

   logical :: flag_verify
   integer :: angular_new
   integer :: l_verify

   integer :: stress_for_relaxation_user
   logical :: flag_stress_for_relaxation_user
   logical :: numerical_stress_needed_for_relaxation
   logical :: compute_forces
   real :: time_read_buffer

   integer :: info
   character*800 :: info_str
   character*200 :: buffer
   character*200 :: buffer2
   integer, parameter :: buffer_arr_siz = 10
   character(len=50) :: buffer_arr(buffer_arr_siz), tag, key, val

   character*10 :: unfold_format

   integer :: i_space

   integer :: n_kpts

   integer :: temp_cube_index
!  counters

   integer :: i_cube_soc
   integer :: i_counter
   integer :: i_species
   integer :: i_cube
   integer :: cube_id
   integer :: i_write_soc_eigenvectors
   integer :: i_cube_stm
   integer :: i_flag, i_tag
   integer :: i_region
   integer :: j_region
   integer :: i_division
   integer :: i_vdw_ignore
   integer :: i_MD_QH_init_segments
   integer :: order, i_ord, n_indep
   integer :: default_max_n_prodbas, default_max_l_prodbas
   real*8 :: default_prodbas_acc
   logical :: l_dummy
   integer :: i_esp
   integer :: split_val

   integer :: i_spin
   integer :: index1, index2
   integer :: i_force_occ
   integer :: i_sub_boys
   integer :: MD_schedule_step
   integer :: i_pp_species
   character(*), parameter :: func = 'read_control'

   real*8, dimension(:), allocatable :: range_a
   real*8 :: radius_a, radius_b
   integer :: range_fragment_x
   integer :: range_fragment_y
   integer :: range_fragment_z

   ! fo-dft
   character(len=8) :: myid_string

   ! LibXC/AJL Dec16
   integer :: flag_libxc_x_id
   integer :: flag_libxc_c_id

   ! ELSI restart density matrix
   character*100 :: dm_file

!  begin work

   call localorb_info('',use_unit)
   call localorb_info( &
   "------------------------------------------------------------", &
         use_unit,'(A)')
   call localorb_info("Reading file control.in.", &
         use_unit,'(10X,A)')
   call localorb_info( &
   "------------------------------------------------------------", &
         use_unit,'(A)')

   call initialize
   call readfile
   call checkinput

   call localorb_info( &
   "------------------------------------------------------------", &
         use_unit,'(A)')
   call localorb_info('',use_unit)

   return

 contains

   subroutine initialize
!    initialize
      flag_relative_fp_charge_mix = .false.
      flag_delta_numerical_stress = .false.

      flag_KS_method = .false.
      flag_mbd_eigensolver = .false.
      flag_hybrid_coeff = .false.
      flag_use_local_index = .false.
      flag_use_local_index_if_scalapack = .false.
      flag_load_balancing = .false.
      flag_load_balancing_if_scalapack = .false.
      flag_use_alltoall = .false.
      flag_MB_clean_rotations_read = .false.
      flag_freq_grid_type = .false.
      flag_scgw_no_core = .false.
      flag_maximum_frequency = .false.
      flag_maximum_time = .false.
      flag_MD_QH_file_exists = .false.
      flag_use_logsbt = .false.
      flag_sbtgrid_lnr0 = .false.
      flag_sbtgrid_lnk0 = .false.
      flag_sbtgrid_lnrange = .false.
      flag_sbtgrid_N = .false.
      flag_use_logsbt_for_radial_hse_integration = .false.
      flag_use_2d_corr = .false.
      flag_calculate_atom_bsse = .false.
      flag_default_moment_defined = .false.
      flag_default_max_n_prodbas = .false.
      flag_default_max_l_prodbas = .false.
      flag_default_prodbas_acc = .false.
      flag_MD_tstep = .false.
      flag_hse_unit = .false.
      flag_hse03 = .false.
      flag_crit_val = .false.
      flag_coul_mat_threshold = .false.
      flag_RI = .false.
      flag_restrict_kpoint_to_smp_node = .false.
      flag_max_tasks_per_smp_node = .false.
      flag_multip_moments_threshold = .false.
      flag_write_restart_geometry = .false.
      flag_hessian_to_restart_geometry = .false.
      flag_compensate_multipole_errors = .false.
      flag_normalize_initial_density = .false.
      flag_vdw_convergence_threshold = .false.
      flag_calculate_friction = .false.
      flag_cavity_restart_file_exists = .false.
      flag_elsi_method = .false.

      flag_out_mommat = .false.

      flag_collect_eigenvectors = .false.

      flag_stress_for_relaxation_user = .false.
      numerical_stress_needed_for_relaxation = .false.

      n_kpts = 1 ! Solve one eigenproblem also in the non-periodic case.

      i_species = 0
      i_start_mp2 = 1
      i_vdw_ignore = 0
      n_wf_funcs = 0
      i_pp_species = 0

      i_cube = 0
      i_write_soc_eigenvectors = 0
      i_cube_stm = 0
      i_MD_QH_init_segments = 0
      i_force_occ = 0
      i_sub_boys = 0

      range_fragment_x = 0
      range_fragment_y = 0
      range_fragment_z = 0
!    provide reasonable defaults

      use_mpi_in_place = .true.

      if (use_constraint) then
         do i_region=1,n_region
            constraint_electrons(i_region,1)=0.d0
            constraint_electrons(i_region,2)=0.d0
         enddo
         flag_region = .false.
         constraint_mix(1)=0.5d0
         constraint_mix(2)=0.5d0
         constraint_precision=1.d-6
      end if

      j_region = 0

      charge = 0.d0
      hydro_cut = .true.
      flag_xc = -1
      flag_xc_pre = -1
      n_steps_xc_pre = 0
! Flags for libXC functionals. AJL/Aug 2016
      flag_libxc_x_id = 0
      flag_libxc_c_id = 0
!
      flag_atomic_xc = 8
      flag_rel = REL_not_set
      force_potential = 0
      !flag_xc2 = 8
      flag_fxc = -1
      flag_eg = .false.
      if (use_cube_output) then
         cube_stm(3,:) = -1.0
      end if
      ! VB: By default, always use the neutralized (accurate) version of the energy functional.
      force_new_functional = .true.

      ! Hidden parameter that enables strict enforcement of smooth
      ! basis function behavior near the cutoff potential.
      force_smooth_cutoff = .false.

      flag_force = .false.
      flag_acc_forces = .false.
      flag_acc_stress = .false.
      flag_relax_mode = .false.
      flag_molecular_dynamics = .false.
      flag_MD_schedule = .false.
      flag_prodbas_threshold = .false.

      flag_wf_named = .false.
      flag_wf_func = .false.

      ipi_dip = .false.
      ipi_hirshfeld = .false.
      ipi_ftensor = .false.
      ipi_work = .false.
      ipi_atomic_stress = .false.

      out_basis = .false.
      out_ovlp_spectrum = .false.
      out_eigenvec = .false.
      out_overlap  = .false.
      out_hamiltonian  = .false.
      out_matrices = .false.
      out_matrices_format_2005 = .false.
      out_dipole_moment = .false.
      out_quadrupole_moment = .false.
      out_grids = .false.
      out_v_eff = .false.
      out_v_eff_new = .false.
      out_dielec_func = .false.
      out_ion_dens = .false.
      out_delta_v = .false.
      out_ascii = .false.
      out_nconsistent=-1
      out_ninitialize=-1
      out_niteration=-1
      out_v_hartree = .false.
      out_rho_multipole = .false.
      out_density = .false.
      out_cube = .false.
      out_cube_soc = .false.
      out_band = .false.
      out_band_mulliken = .false.
      out_band_during_scf = .false.
      out_z2_invariant = .false.
      out_mulliken = .false.
      out_hirshfeld = .false.
      out_hirshfeld_iterative = .false.
      out_hirshfeld_always = .false.
      out_zero_multipoles = .false.
      out_vdwdf = .false.
      out_polarisability = .false.
      out_dipole_moment = .false.
      out_esp = .false.
      out_esp_full = .false.
      i_esp = 0


      full_embedding = .true.

      ! Different meaning than all the other flags!
      flag_hartree_partition_type = 0

      flag_mixer = .false.
      flag_mixer_constraint = .false.
      flag_mix_param = .false.
      flag_prec_mix_param = .false.
      flag_spin_mix_param = .false.
      flag_pulay_iter = .false.
      flag_pulay_iter_constraint = .false.
      flag_ini_linear_mixing = .false.
      flag_linear_mix_param = .false.
      flag_linear_spin_mix_param = .false.
      flag_ini_linear_mixing_constraint = .false.
      flag_acc_rho = .false.
      flag_acc_eev = .false.
      flag_acc_etot = .false.
      flag_acc_potjump = .false.
      flag_sc_max = .false.
      flag_sc_init_iter = .false.
      flag_sc_init_factor = .false.
      flag_constraint_iter = .false.
      flag_constraint_mix = .false.
      flag_constraint_precision = .false.
      flag_partition_acc = .false.
      flag_partition_type = .false.
      flag_basis_threshold = .false.
      flag_wave_threshold = .false.
      flag_force_lebedev = .false.
      flag_occupation_type = .false.
      flag_DFPT_occupation_type = .false.
      flag_occupation_thr = .false.
      flag_occupation_acc = .false.
      flag_mu_method = .false.
      flag_fermi_acc = .false.
      flag_max_zeroin = .false.
      flag_spin = .false.
      flag_angular_min_size = .false.
      flag_angular_div_cutoff = .false.
      flag_hartree_worksize = .false.
      flag_distributed_spline_storage = .false.
      flag_remove_unitary_force_component = .false.

      flag_scsmp2= .false.
      ! This is the flag for MP2 frozen core and frozen virtual treatment
      ! and declared and exported by runtime_choices
      flag_frozen_core= .false.


      flag_energy_tolerance = .false.
      flag_aggregated_energy_tolerance = .false.
      flag_harmonic_length_scale = .false.
      flag_max_atomic_move = .false.
      flag_min_line_step = .false.
      flag_line_step_reduce = .false.
      flag_max_relaxation_steps = .false.
      flag_points_in_batch = .false.
      flag_batch_size_limit = .false.
      flag_grid_partitioning_method = .false.
      flag_grouping_factor = .false.
      flag_init_hess = .false.
      flag_distributed_hessian = .false.

      flag_preconditioner = .false.
      flag_preconditioner_max_l = .false.
      flag_preconditioner_turnoff = .false.

      flag_legacy_monopole_extrapolation = .false.
      flag_packed_matrix_format = .false.
      flag_k_points_defined = .false.
      flag_symmetry_reduced_k_grid = .false.

      flag_mixer_threshold_charge = .false.
      flag_mixer_threshold_spin = .false.
      flag_prodbas_threshold = .false.

      flag_min_batch_size = .false.
      flag_use_metis_batching = .false.

      flag_initial_ev_solutions = .false.
      flag_lopcg_preconditioner = .false.
      flag_lopcg_tolerance = .false.
      flag_max_cg_iterations = .false.
      flag_lopcg_block_size = .false.

      flag_atom_ref=.false.

      flag_multipole_threshold = .false.
      flag_l_hartree_far_distance = .false.
      flag_density_update_method = .false.
      flag_communication_type = .false.

      flag_found_phonon = .false.

!      flag_set_vacuum_level = .false. !OTH: Temporarily(?) moved to runtime-choices
      flag_transport = .false.

      flag_frequency_points = .false.
      flag_anacon_type = .false.
      flag_n_anacon_par = .false.

      anacon_type_defined = .false.

      cube_ready = .false.  !SAG
      use_vdw_method = .false.
      treeflag = 3
      epsfinal = 0.000001d0
      flag_nrad_nlcorr = .false.
      flag_i_leb_nlcorr = .false.

      flag_time_points = .false.

      plus_u_petukhov_mixing_defined = .false.
      plus_u_ramping_defined = .false.
      plus_u_hydros_defined = .false.
      plus_u_occupation_matrix_control_read = .false.
      plus_u_occupation_matrix_control_write = .false.
      plus_u_matrix_error_defined = .false.
      plus_u_matrix_release_defined = .false.
      plus_u_eigenvalues_defined = .false.


      ! Handling of compute_forces and use_forces is a bit of a mess.
      ! use_forces is the general flag that  indicates the use of forces;
      ! it is already set in dimensions.f, and thus not normally reset here.
      ! However, if compute_forces .false. is explicitly requested,
      ! forces should still not be calculated even if, say, a
      ! force convergence accuracy criterion is specified. Thus, we
      ! must check here separately the need to compute forces, and
      ! possibly also reset the use_forces flag to false.

      if (use_forces) then
         ! default is compute_forces = .true. unless explicitly specified below.
         compute_forces = .true.
      end if

      if (use_initial_moment) then
         ! If any initial moments are set explicitly in geometry.in, then
         ! * by default, use a specified (non-Hund-rule) initial moment
         ! * The default moment is zero on all atoms except those that
         !   have an explicitly specified moment.
         ! The default moment can be specified explicitly in control.in, below,
         ! and can also be set to use Hund's rules (through flag_moment=false)
         flag_moment = .true.
         default_initial_moment = 0.
         flag_default_moment_defined = .true.
      else
         ! If no moments specified in geometry.in at all,
         ! set Hund's rules to initialize
         ! single-atom densities.
         ! However, this flag must be reset by default_initial_moment later -
         ! if no moment is set whatsoever, the code now stops except for a single isolated atom.
         flag_moment = .false.
         flag_default_moment_defined = .false.
      end if

      !Default: no postprocessing if SCF did not converge
      postprocess_anyway = PP_ANYWAY_NOTHING

!<<<VVG>>>
          mbd_cfdm_dip_cutoff  = 1000.0
          mbd_scs_dip_cutoff   = 600.d0
          mbd_supercell_cutoff = 10.d0
          mbd_scs_vacuum_axis(1)  = .FALSE.
          mbd_scs_vacuum_axis(2)  = .FALSE.
          mbd_scs_vacuum_axis(3)  = .FALSE.
!<<<VVG>>>

! igor, set the default number of DFTs for printing out
   n_dft_methods = 0

      if(use_contour_def_gw) then
        call set_defaults_cd(gw_cd)
      endif
   end subroutine initialize



   subroutine readfile
      implicit none

!  local variables

!  i_code : Flag to determine i/o status
!  desc_str: line descriptor

      integer :: i_code, linecount, info, idx
      integer :: i_pp_species
      logical :: flag
      logical :: explicit
!      logical :: cube_type_needs_densmat

      character*250 inputline
      character*40 desc_str
      ! Variables needed for LibXC
      ! Maybe these can be consolidated with other temp variables?
      character*40 :: libxc_x_str, libxc_c_str
      integer :: libxc_index_for_split, libxc_x_type, libxc_c_type = 0
      integer :: dummy_int, tmp
      real*8 :: dummy_real
      !
      integer :: ex_code
      !

      open (7, FILE="control.in")
      mbd_std_flags(:) = ''
      linecount = 0
      do

         read(7,'(A)',iostat=i_code) inputline
         if(i_code<0) exit               ! end of file reached
         if(i_code>0) then
            call aims_stop_coll("Unknown error reading file 'control.in'...", func)
         endif

         linecount = linecount + 1

         read(inputline,*,iostat=i_code) desc_str
         if(i_code/=0) cycle              ! skip empty line

         if (desc_str(1:1).eq.'#') cycle  ! skip comment

!     decode contents of current line

         select case (desc_str)

         case('sym_precision')
            read(inputline,*,end=88,err=99) desc_str, sym_precision
            if ( .not. use_spglib ) then
              write(info_str,'(2X,A)') &
                '*ERROR: Precision for symmetry operations specified without spglib support.'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: Please remove keyword from control.in file or compile FHI-aims"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: with spglib support, e.g. by setting 'USE_SPGLIB = yes' and 'USE_C_Files = yes'"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: in your make.sys file and by specifying a C-Compiler "
              call localorb_info(info_str,use_unit,'(A)')

              call aims_stop_coll()
            end if

         case('use_symmetric_forces')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               use_symmetric_forces = .false.
            else
               use_symmetric_forces = .true.
            endif
            if ( .not. use_spglib .and. use_symmetric_forces) then
              write(info_str,'(2X,A)') &
                '*ERROR: Asked to perform force symmetrization without spglib support.'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: Please set 'use_symmetric_forces' to false in control.in file or compile FHI-aims"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: with spglib support, e.g. by setting 'USE_SPGLIB = yes' and 'USE_C_Files = yes'"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: in your make.sys file and by specifying a C-Compiler "
              call localorb_info(info_str,use_unit,'(A)')

              call aims_stop_coll()
            end if
         case('potential_well') !SYNTAX: potential_well z_start (in A) z_end (in A) depth (in eV)
           read(inputline,*,end=88,err=99) desc_str,  potential_well_start, potential_well_end, potential_well_depth
           write(info_str,'(2X,A)') '*Potential well requested. Warning, feature is experimental.*'
           potential_well_start = potential_well_start / bohr
           potential_well_end = potential_well_end / bohr
           potential_well_depth = potential_well_depth/hartree
           potential_well_requested = .true.


         case('symmetry_reduced_k_grid_spg')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               use_symmetry_reduced_spg = .false.
            else
               use_symmetry_reduced_spg = .true.
               use_symmetry_reduced_k_grid = .true.
               flag_symmetry_reduced_k_grid = .true.
            endif
            if ( use_hf_kspace .and. use_symmetry_reduced_spg) then
              write(info_str,'(2X,A)') &
                '*ERROR: k-space periodic HF implementation and spglib k-point reduction not supported.'
              call localorb_info(info_str,use_unit,'(A)')
              call aims_stop_coll()
            end if
            if ( .not. use_spglib .and. use_symmetry_reduced_spg) then
              write(info_str,'(2X,A)') &
                '*ERROR: Asked to perform k-point reduction without spglib support.'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: Please set 'symmetry_reduced_k_grid_spg' to false in control.in file or compile FHI-aims"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: with spglib support, e.g. by setting 'USE_SPGLIB = yes' and 'USE_C_Files = yes'"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: in your make.sys file and by specifying a C-Compiler "
              call localorb_info(info_str,use_unit,'(A)')

              call aims_stop_coll()
            end if
         case('use_spg_full_Delta')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               use_spg_full_Delta = .false.
               write(info_str,'(2X,A)') &
                'Not using full Delta matrix for eigenvector rotation.'
              call localorb_info(info_str,use_unit,'(A)')
            else
               use_spg_full_Delta = .true.
               write(info_str,'(2X,A)') &
                'Using full Delta matrix for eigenvector rotation.'
              call localorb_info(info_str,use_unit,'(A)')
            endif
         case('reconstruct_proper_only')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               reconstruct_proper_only = .false.
               write(info_str,'(2X,A)') &
                'Reconstructing all eigenvectors.'
              call localorb_info(info_str,use_unit,'(A)')
            else
               reconstruct_proper_only = .true.
               write(info_str,'(2X,A)') &
                'Excluding inversion from eigenvector reconstruction .'
              call localorb_info(info_str,use_unit,'(A)')
            endif
          case('use_spg_mv_mm')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               use_spg_mv_mm = .false.
               write(info_str,'(2X,A)') &
                'Using 2 x Matrix-Matrix operation for density reconstruction.'
              call localorb_info(info_str,use_unit,'(A)')
            else
               use_spg_mv_mm = .true.
               write(info_str,'(2X,A)') &
                'Using Matrix-Vector, Matrix-Matrix operation for density reconstruction.'
              call localorb_info(info_str,use_unit,'(A)')
            endif
          case('get_full_density_first')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               get_full_density_first = .false.
               write(info_str,'(2X,A)') &
                'Full KS-density constructed on the fly for Fock-Matrix.'
              call localorb_info(info_str,use_unit,'(A)')
            else
               get_full_density_first = .true.
               write(info_str,'(2X,A)') &
                'Full KS-density constructed before Fock-Matrix evaluation.'
              call localorb_info(info_str,use_unit,'(A)')
            endif
          case('use_k_phase_exx')
	    read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               use_k_phase_exx = .false.
               write(info_str,'(2X,A)') &
                'Not storing phase factor for EXX.'
              call localorb_info(info_str,use_unit,'(A)')
            else
               use_k_phase_exx = .true.
               write(info_str,'(2X,A)') &
                'Storing phase factor for EXX.'
              call localorb_info(info_str,use_unit,'(A)')
            endif

         case('use_symmetry_analysis')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
              use_symmetry_analysis = .false.
              set_symmetry_analysis = .true.
            else
              use_symmetry_analysis = .true.
              set_symmetry_analysis = .true.
            endif
            if ( .not. use_spglib .and. use_symmetry_analysis) then
              write(info_str,'(2X,A)') &
                '*ERROR: Ask to perform symmetry analysis without spglib support.'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: Please set 'use_symmetry_analysis' to false or compile FHI-aims"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: with spglib support, e.g. by setting 'USE_SPGLIB = yes' and 'USE_C_Files = yes'"
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') &
                "*ERROR: in your make.sys file and by specifying a C-Compiler "
              call localorb_info(info_str,use_unit,'(A)')

              call aims_stop_coll()
            end if


         case('force_complex_eigenvectors') ! CC
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               write(info_str,'(2X,A)') 'Forcing the eigenvectors in the calculation to be COMPLEX.'
               call localorb_info(info_str,use_unit,'(A)')
               flag_force_complex_eigenvectors = .true.
            else
               write(info_str,'(2X,A)') 'Not enforcing the eigenvectors in the calculation to be COMPLEX.'
               call localorb_info(info_str,use_unit,'(A)')
               flag_force_complex_eigenvectors = .false.
            end if

         case('gamma_cut_coulomb')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               write(info_str,'(2X,A,A)') 'Using the cut coulomb potential to lift the gamma point singlarity', &
                                   ' in post-DFT correlation calculations.'
                gamma_cut_coulomb=.true.
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,A,A)') 'Using the mixed way to lift the gamma point singlarity', &
                                   ' in post-DFT correlation calculations.'
               call localorb_info(info_str,use_unit,'(A)')
            endif

         case('use_gw_gamma_corr')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               write(info_str,'(2X,A,A)') 'Using explicit compensating term to treat the gamma point singularity', &
                                   ' in periodic GW calculations. Only experimental at the moment!'
                use_gw_gamma_corr = .true.
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,A,A)') 'Explicit gamma point correction is not invoked in periodic', &
                                   ' GW calculations.'
               call localorb_info(info_str,use_unit,'(A)')
            endif
         case('dry_run')

               write(info_str, '(2X,A)') &
               & "Found the 'dry_run' option. This run will only read control.in and geometry.in."
               call localorb_info(info_str)
               write(info_str, '(2X,A)') &
               & "The run will stop before any work-intensive part (i.e., before initialize_scf is first invoked)."
               call localorb_info(info_str)

               dry_run = .true.

         case('el_ph_coupling')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               call aims_stop_coll("electron phonon coupling can not be used yet ...", func)
            else
               call aims_stop_coll("electron phonon coupling keyword found - experimental, stopping ...", func)
            end if

         case('charge')
            read(inputline,*,end=88,err=99) desc_str, charge
            if (charge.ne.0.d0) then
               write(info_str, '(2X,A,E14.6,A)') &
               & "Charged system requested: Charge = ", charge, "."
               call localorb_info(info_str)
            else
               write (info_str, '(2X,A,E14.6,A)') &
               & "Charge = ", charge, &
               & ": Neutral system requested explicitly."
               call localorb_info(info_str)
            end if

         case('spin')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            ! Corner case: Two spin different definitions in same control.in
            !              will not work - check first.
            if (flag_spin) then
               ! This must be a duplicate definition - spin has been set before
               ! Check if we conflict:

               if ( (desc_str.eq.'none').and.(spin_treatment.eq.1) .or. &
                    (desc_str.eq.'collinear').and.(spin_treatment.eq.0) ) then
                  ! The previous and present spin definition conflict.
                  ! We cannot silently ignore the problem.

                  write(info_str,'(1X,A)') &
                    "* Error: conflicting 'spin' keywords in input file 'control.in'. "
                  call localorb_info(info_str)

                  write(info_str,'(1X,A)') &
                    "* The offending spin definition reads: "
                  call localorb_info(info_str)

                  call localorb_info(inputline)

                  write(info_str,'(1X,A)') &
                    "* Please remove the duplicate spin definition from 'control.in' and rerun."
                  call localorb_info(info_str)

                  call aims_stop_coll('Duplicate spin definition', func)

               end if

            end if

            ! Normal case follows: No duplicate spin definition, in this case define spin as planned.

            if (desc_str.eq.'none') then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                     "Spin treatment: No spin polarisation."
               end if

               spin_treatment = 0

            else if (desc_str.eq.'collinear') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                     "Spin treatment: Spin density functional ", &
                     "theory - collinear spins."
               end if

               if (calculate_perturbative_soc) then
                 write(info_str,'(2X,A)') &
                 & "You have request a SOC run with spin-polarization."
                 call localorb_info(info_str)
                 write(info_str,'(2X,A)') &
                 & "symmetry_reduced_k_grid will be turned off."
                 call localorb_info(info_str)
                 use_symmetry_reduced_k_grid = .false.
               end if

               spin_treatment = 1

            else

               write(info_str,'(1X,A,A,A)') "* Spin treatment: ", desc_str, "?"
               call localorb_info(info_str)
               write(info_str,'(1X,A,A)') &
               & "* This type of spin treatment is not (yet?) defined - ", &
               & "please correct and rerun."
               call localorb_info(info_str)

               call aims_stop_coll('Invalid spin treatment', func)

            end if

            flag_spin = .true.

         case('default_initial_moment')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if ( (desc_str.eq.'hund') .or. (desc_str.eq.'Hund') ) then
               if (n_atoms .gt. 1 .or. n_periodic .gt. 0) then
                 write(info_str,'(1X,A)') "** You have specified default initial spin moments set by Hund's rules for a system"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** with more than one atom.  Hund's rules may give very bad initial spin moments (high"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** spin) for systems other than isolated atoms, leading to problematic initial"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** densities for SCF.  While it *could* be the case that your system consists of"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** isolated atoms, 'default_initial_moment hund' for polyatomic systems has done far"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** more harm than good for others in the past.  If you do wish to use Hund's rules to"
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** initialize specific atoms, please do so by hand.  Apologizes for the inconvenience."
                 call localorb_info(info_str)
                 write(info_str,'(1X,A)') "** Exiting."
                 call localorb_info(info_str)
                 call aims_stop_coll('', func)
               end if

               write (info_str,'(2X,A)') &
               "Initial spin moment per atom set according to Hund's rules."
               call localorb_info(info_str,use_unit,'(A)')

               flag_moment = .false.

            else

               read(inputline,*,end=88,err=99) desc_str, default_initial_moment

               write (info_str,'(2X,A,F8.4)') &
               "Initial spin moment per atom: ", &
               default_initial_moment
               call localorb_info(info_str,use_unit,'(A)')

               flag_moment = .true.

            end if

            ! record that a default initial moment was defined.
            flag_default_moment_defined = .true.

         case('use_aufbau_principle')
             read(inputline,*, iostat=i_code) desc_str, use_aufbau_principle

         case('verbatim_writeout')
              continue !Effect has been used in parse_control already


         case('compute_forces')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            ! Notice that use_forces is already checked in
            ! parse_forces(), module dimensions.f

            if (desc_str.eq.'.false.') then

               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)') &
                        "Forces will not be computed."
               end if

               compute_forces = .false.

            else if (desc_str.eq.'.true.') then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)') &
                        "Forces will be computed."
               end if

               compute_forces = .true.

            else
               if (myid .eq. 0) then
                  write (use_unit,'(1X,A,A,A)') &
                        "* compute_forces: ", desc_str, "?"
                  write (use_unit,'(1X,A,A)') &
                        "* keyword not defined - ", &
                        "please correct and rerun."
               end if
               call aims_stop_coll('', func)

            end if

            flag_force = .true.

         case('compute_numerical_stress')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'.true.') then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Numerical stress will be computed using finite differences."
               end if
            end if
         case('delta_numerical_stress')
            flag_delta_numerical_stress = .true.
            read(inputline,*,end=88,err=99) desc_str, delta_numerical_stress

            if (use_numerical_stress .and. (delta_numerical_stress<=0)) then
               delta_numerical_stress = default_delta_numerical_stress
               if (myid .eq. 0) then
               write (use_unit,'(2X,A)')  &
               "Scaling factor delta for numerical stress was set <= 0. This does not make sense (in our convention)!"
               write (use_unit,'(2X,A,E14.6,A)')  "Delta for numerical stress set to default of ", delta_numerical_stress, "."
               end if
            else if (use_numerical_stress .and. (delta_numerical_stress>0)) then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A,E14.6,A,2X)')  &
                     "Scaling factor delta for numerical stress set to ", &
                     delta_numerical_stress, "."
               end if
            end if

         case('numerical_stress_save_scf')
            read(inputline,*,end=88,err=99) desc_str, numerical_stress_save_scf
            if (numerical_stress_save_scf) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                  'Explicitly using symmetries implied by constrained relaxations to save time on numerical stress calculation.'
               end if
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') &
                  'Not using symmetries implied by constrained relaxations ', &
                  'for numerical stress calculation.'
               end if
            end if
         case('compute_analytical_stress')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               use_analytical_stress = .true.
               use_basis_gradients = .true.
            end if
         ! CC: Compute Heat Flux along MD
         case('compute_heat_flux')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               compute_heat_flux = .true.
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') "Heat flux will be computed."
               end if
            end if
         ! If symmetrization is explicitly not requested, set AS_components to 9.
         case('calc_analytical_stress_symmetrized')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               AS_components = 9
            end if

         ! FK: This keyword is currently not used. Both unsymmetrized and symmetrized analytical stress is printed.
         ! If all 9 components of the analytical stress are calculated, one can choose
         ! whether the final output is symmetrized or not.
         case('output_analytical_stress_symmetrized')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               AS_symmetric_output = .true.
            end if

         ! Read which stress was requested for unit cell relaxation.
         case('stress_for_relaxation')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if ( (desc_str.eq.'numerical') .or. (desc_str.eq.'NUMERICAL') ) then
               stress_for_relaxation_user = RELAX_NUMERICAL
            else if ( (desc_str.eq.'analytical') .or. (desc_str.eq.'ANALYTICAL') ) then
               stress_for_relaxation_user = RELAX_ANALYTICAL
            else
               write (info_str,'(2X,A,A,A)') &
                  "*** ", desc_str, " is not a defined flag for the tag stress_for_relaxation."
               call localorb_info(info_str,use_unit,'(A)')
               write (info_str,'(2X,A,A,A)') &
                  "*** Please use either analytical or numerical as flag."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if
            flag_stress_for_relaxation_user = .true.

         ! Read input for external pressure
         case('external_pressure')
            read(inputline,*,end=88,err=99) desc_str, external_pressure
            if (external_pressure .ne. 0.0d0) then
               if (n_periodic .eq. 0) then
                  write (info_str,'(2X,A,A,A)') &
                     "*** Applying external pressure does not work for non-periodic systems."
                  call localorb_info(info_str,use_unit,'(A)')
                  call aims_stop_coll('', func)
               end if

               write (info_str,'(2X,A,E10.3,A,E10.3,A)') &
                  "External pressure of ", external_pressure, &
                  " eV/A**3 (", external_pressure * giga_pascal , &
                  " GPa) is applied to system."
               call localorb_info(info_str,use_unit,'(A)')
               ! Convert unit from eV/A**3 to Ha/bohr**3
               external_pressure = external_pressure * (bohr**3)/hartree
               flag_external_pressure = .true.
            end if

         case('analytic_potential_average')
            read(inputline,*,end=88,err=99) desc_str, analytic_potential_average
            if (n_periodic .gt. 0) then
               if (analytic_potential_average) then
                  write (info_str,'(2X,A,A,A)') &
                     "Using analytically calculated rather than 3d integrated potential averages."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  write (info_str,'(2X,A,A,A)') &
                     "Using numerically integrated potential averages."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            else
                  write (info_str,'(2X,A,A,A)') &
                     "Non-periodic calculation - ignoring analytic_potential_average setting."
                  call localorb_info(info_str,use_unit,'(A)')
                  analytic_potential_average = .false.
            end if

         case('relax_geometry')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'none') then

               write (info_str,'(2X,A,A)') &
               "Geometry relaxation: ", &
               "No relaxation will be performed."
               call localorb_info(info_str,use_unit,'(A)')

               relax_mode = RELAX_OFF
               relax_geo = RELAX_OFF

            else if ((desc_str.eq.'bfgs_textbook').or.(desc_str.eq.'BFGS_TEXTBOOK')) then

               write (info_str,'(2X,A,A)') &
               "Geometry relaxation: ", &
               "Textbook BFGS (simple quadratic extrapolation)."
               call localorb_info(info_str,use_unit,'(A)')

               relax_mode = RELAX_BFGS_TB
               relax_geo = RELAX_BFGS_TB

            else if (str_lower(desc_str) .eq. 'trm_2012') then
               ! VB: Since TRM is a modified variant of BFGS and much more efficient in all tests in all tests
               !     over the course of a year, it is now the default also if BFGS is chosen.

               write (info_str,'(2X,A,A)') &
               "Geometry relaxation: ", &
               "Modified BFGS - TRM (trust radius method, i.e., trusted quadratic extrapolation)."
               call localorb_info(info_str,use_unit,'(A)')

               relax_mode = RELAX_TRM
               relax_geo = RELAX_TRM

            ! FlK: This is the adapted trm for periodic system. It uses effective strain
            ! transformations to change the lattice degrees of freedom, which fixes
            ! several issues when relaxing non-cubic systems.
            ! 29.11.2018: In accordance with Volker, this becomes the default
            ! relaxation method with the keyword bfgs or trm.
            else if ((str_lower(desc_str) .eq. 'trm') &
                .or. (str_lower(desc_str) .eq. 'bfgs')) then

               write (info_str,'(2X,2A)') &
               "Geometry relaxation: ", &
               "Modified BFGS - TRM (trust radius method) for lattice optimization."
               call localorb_info(info_str,use_unit,'(A)')

               relax_mode = RELAX_LTRM
               relax_geo = RELAX_LTRM

            ! 04.03.19: Give option to activate additional developer features by
            !           explicitly requesting `ltrm` or `lattice_trm`
            else if ((str_lower(desc_str) .eq. 'ltrm') &
                .or. (str_lower(desc_str) .eq. 'lattice_trm')) then

               write (info_str,'(2X,2A)') &
               "Geometry relaxation: ", &
               "Modified BFGS - TRM (trust radius method) for lattice optimization."
               call localorb_info(info_str,use_unit,'(A)')

               write (info_str,'(2X,2A)') &
               "** Possibly untested developer features enabled. If this is ", &
               "unwanted, please use the keyword `trm` instead of `lattice_trm`"
               call localorb_info(info_str,use_unit,'(A)')

               relax_mode = RELAX_LTRM
               relax_geo = RELAX_LTRM
               RELAX_LTRM_dev = .true.

            else

               write (info_str,'(1X,A,A,A)') &
               "* Geometry relaxation: ", desc_str, &
               " is not defined as a relaxation algorithm."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)

            end if

            if (relax_mode.ne.RELAX_OFF) then

               read(inputline,*,end=88,err=99) desc_str, desc_str, relax_accuracy_forces

               if (relax_accuracy_forces.ge.0.d0) then

                  write (info_str,'(2X,A,A,E14.6,A)') &
                  "Convergence accuracy for geometry relaxation: ", &
                  "Maximum force < ", relax_accuracy_forces," eV/A."
                  call localorb_info(info_str,use_unit,'(A)')

                  ! FlK: inform about using stress as optimization criterion
                  if ((relax_geo .eq. RELAX_LTRM) .and. (n_periodic > 0) .and. &
                    use_stress_to_check_relaxation) then

                    if (relax_accuracy_forces < 0.001d0) then
                      relax_accuracy_stress_factor = 1.d0
                    end if

                    write (info_str,'(2X,A,A,E14.6,A)') &
                    "Convergence accuracy for geometry relaxation: ", &
                    "Maximum stress component < ", relax_accuracy_stress_factor * relax_accuracy_forces," eV/A^3."
                    call localorb_info(info_str,use_unit,'(A)')
                  end if

               else

                  write (info_str,'(2X,A,E14.6,A)') &
                  "* Error: Desired relaxation accuracy ", &
                  relax_accuracy_forces, &
                  " should be greater than 0."
                  call localorb_info(info_str,use_unit,'(A)')
                  call aims_stop_coll('', func)

               end if

            end if

            flag_relax_mode = .true.

         case('relax_unit_cell')
            ! FK: If unit cell relaxation is set to none, don't stop aims for non-periodic systems.
            !     Therefore, check for periodic system will be performed after check for 'none' flag.
            read (inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'none') then
               write (info_str,'(2X,A,A)') &
                  "Unit cell relaxation: No relaxation of the unit cell will be performed."
               call localorb_info(info_str,use_unit,'(A)')
            else if (n_periodic .eq. 0) then
               write (info_str,'(1X,A)') &
                  "* Unit cell relaxation does not make sense in the cluster case!"
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll ('', func)
            ! FK: The variable for the type of unit cell relaxation was already set in dimensions.f90.
            !     Therefore, just repeat the possible flags without doing anything.
            else if (desc_str.eq.'full') then
               continue
            else if ( (desc_str.eq.'fixed_angles') .or. (desc_str.eq.'shape') ) then
               continue
            else
               write (info_str,'(1X,A,A,A)') &
                  "* Unit cell relaxation: ", desc_str, " is not a defined keyword."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll ('', func)
            end if

         case('output_vacuum_length')
            read(inputline,*,end=88,err=99) desc_str, vacuum_length
         !case('output_realspace_esp')
         !   read(inputline,*,end=88,err=99) desc_str, step_x, step_y, step_z
         !
         !   flag_out_elec_real = .true.
         !   i_cube = i_cube+1
         !   cube_type(i_cube)='hartree_potential'
         !   cube_type_needs_densmat(i_cube)=.false.
         !   cube_type_needs_eigenvectors(i_cube)=.false.
         !   !cube_filename(i_cube) = 'electrostatic_potential_realspace'
         !   call read_cube_parameters ( i_cube )
         !   out_cube = .true.
            !define the default origins for cube output
            !cube_origin(:,i_cube) = 0.0d0
            !flag_origin(i_cube) = .true.
            !define the steps and edges for cube output
            !cube_edge_steps(1,i_cube) = SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)/0.1
            !cube_edge_unit(:,1,i_cube) = lattice_vector(:,1)/cube_edge_steps(1,i_cube)
            !cube_edge_steps(2,i_cube) = SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)/0.1
            !cube_edge_unit(:,2,i_cube) = lattice_vector(:,2)/cube_edge_steps(2,i_cube)
            !cube_edge_steps(3,i_cube) = SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)/0.1
            !cube_edge_unit(:,3,i_cube) = lattice_vector(:,3)/cube_edge_steps(3,i_cube)
            !flag_edge(i_cube) = .true.
        ! case('output_locpot_atom')
        !    flag_out_locpot_atom = .true.
         case('compute_dielectric')
            read(inputline,*,end=88,err=99) desc_str, &
                                   omega_max, n_omega
            ! Check for load_balancing first, since it will automatically set use_local_index in some cases
            !if (use_load_balancing) then
            !   call aims_stop_coll ('* Dielectric calculations currently do not support load_balancing or &
            !                        &use_local_index keywords.  Exiting', func)
            !end if
            !if (use_local_index) then
            !   call aims_stop_coll ('* Dielectric calculations currently do not support use_local_index keyword.  Exiting', func)
            !end if

            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Momentum Matrix."
            end if

            flag_out_dielectric = .true.
            use_basis_gradients = .true.
            !use_scalapack = .false.
            if (flag_out_dielectric .and. .not. flag_broadening) then
               call read_plot_dielectric(inputline)
            endif

         case('dipole_trans_mat')
            dipole_trans_mat_out = .true.
            read(inputline,*,end=88,err=99) desc_str, dipole_trans_mat_k, dipole_trans_mat_orbit_num
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Output dipole transition matrix elements when calcualting dielectric function"
            end if

         case('dielectric_broadening')
            read(inputline,*,end=88,err=99) desc_str, die_method, &
                                           die_broaden_para
            ! Check for load_balancing first, since it will automatically set use_local_index in some cases
            !if (use_load_balancing) then
            !   call aims_stop_coll ('* Dielectric calculations currently do not support load_balancing or &
            !                        &use_local_index keywords.  Exiting', func)
            !end if
            !if (use_local_index) then
            !   call aims_stop_coll ('* Dielectric calculations currently do not support use_local_index keyword.  Exiting', func)
            !end if
            flag_broadening = .true.
            if (die_method == 'gaussian') then
               die_broaden_type = 1
            elseif (die_method == 'lorentzian') then
               die_broaden_type = 2
            endif

            if (flag_out_dielectric) then
               call read_plot_dielectric(inputline)
            endif



         case('compute_kubo_greenwood')
            read(inputline,*,end=88,err=99) desc_str, widthone, widthtwo, Emin,&
                                  Emax, omega_min, omega_max, n_omega, ep1, ep2
            if (desc_str.eq.'.true.') then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Greenwood-Kubo optical conductivity and Seebeck coeff."
               end if
            flag_out_greenwood = .true.
            end if

         case('greenwood_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'full') then
               greenwood_method = 1
            else if (desc_str.eq.'sparse') then
               greenwood_method = 2

            else if (desc_str.eq.'both') then
               greenwood_method = 3
               endif
          flag_greenwood_method = .true.
          ! option "both" remains for testing and comparing purposes

         case('compute_momentummatrix')
!            read(inputline,*,end=88,err=99) desc_str, n_state_min, n_state_max, k_point_mom
            read(inputline,*,end=88,err=99) desc_str, Emin, Emax, k_point_mom
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Momentum Matrix."
            end if

            out_mommat = .true.       ! This variable lives in runtim_choices.f90 and tells the rest
                                      ! of the code that the momentum matrix will be needed.

            flag_out_mommat = .true.  ! This flag is only intended for use inside read_control.f90
                                      ! and tells the code that the compute_momentummatrix keyword was
                                      ! set by the user.
         case('compute_dipolematrix')
            read(inputline,*,end=88,err=99) desc_str, Emin, Emax, k_point_mom
            if (myid .eq. 0) then
                write (use_unit,'(2X,A)')  "Compute Dipole Matrix."
            end if
            flag_out_dipmat = .true.
            use_symmetry_reduced_k_grid = .false.
            if (use_scalapack) then
              collect_eigenvectors = .true.
            endif
         case('compute_dipolematrix_k_k')
            read(inputline,*,end=88,err=99) desc_str, Emin, Emax, k_k_method
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Dipole Matrix (k,k')."
            end if
            flag_out_dipmat_k_k = .true.
            use_symmetry_reduced_k_grid = .false.

         case('compute_absorption')
            read(inputline,*,end=88,err=99) desc_str, widthone, Emin, Emax, &
                                   omega_min, omega_max, n_omega, ep1, use_gauss
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Absorption."
            end if
            flag_out_absorption = .true.
            use_basis_gradients = .true.

         case('compute_coulomb_matrix_ovl')
            read(inputline,*,end=88,err=99) desc_str, Emin, Emax, read_q_points
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute Coulomb Matrix from Overlap."
            end if
            flag_out_coulmat_ovl = .true.
            use_symmetry_reduced_k_grid = .false.

         case('output_KS_eigenvector_scalapack_hdf5')
            read(inputline,*,end=88,err=99) desc_str
            if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Paralell out put of KS eigenvectors."
            end if
            flag_out_ev_scalapack_hdf5 = .true.
            use_symmetry_reduced_k_grid = .false.
            use_full_spectrum = .true. ! Change to calculate_all_eigenstates?
         case('compute_coulomb_matrix_lvl')
            read(inputline,*,end=88,err=99) desc_str, Emin, Emax
            if (myid .eq. 0) then
                write (use_unit,'(2X,A)')  "Compute Coulomb Matrix with prod. basis"
            end if
            flag_out_coulmat_lvl = .true.
            use_lvl = .true.
            use_prodbas=.true.
            !use_hartree_fock=.true.
            use_periodic_hf=.true.
            use_hf_kspace=.true.
            use_cutCb=.true.
            use_logsbt=.true.
            RI_type=RI_LVL
            use_symmetry_reduced_k_grid = .false.
            !use_logsbt_for_radial_hse_integration = .true.
         case('compute_esp_charges')
            read(inputline,*,end=88,err=99) desc_str, esp_Emin, esp_Emax,&
                                            esp_grid_type, &
                                            esp_min_radius, esp_max_radius ,&
                                            esp_n_max_radius, esp_k_point
            if (myid .eq. 0) then
                write (use_unit,'(2X,A)')  "Computing selected ESP cahrges"
            end if
            out_esp_full = .true.
            if (use_scalapack) then
              collect_eigenvectors = .true.
            endif
         case('compute_greenwood_dc_transport')
            read(inputline,*,end=88,err=99) desc_str, n_greenenergy
            if (desc_str.eq.'.true.') then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Compute d.c. transport coefficients."
               end if
              flag_out_dclimit = .true.
            end if

         case('explicit_fermi')
            read(inputline,*,end=88,err=99) desc_str, dist_from_vbm, desc_str
            if ((desc_str.eq.'eV')  .or. (desc_str.eq.'ev') ) then
               percent_or_ev=1
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Setting Fermi level ... eV above the VBM"
               end if
             else if ((desc_str.eq.'Percent')  .or. (desc_str.eq.'percent') ) then
               percent_or_ev=2
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Setting Fermi level ... percent above the VBM into the gap."
             end if
            end if
               flag_explicit_fermi = .true.

         case('flexible_fermi_level')
            read(inputline,*,end=88,err=99) desc_str, fermistatesabove, fermistatesbelow,fermispacing
            if (desc_str.eq.'.true.') then
               if (myid .eq. 0) then
                  write (use_unit,'(2X,A)')  "Using given (artificial) Fermi-levels for Kubo-Greenwood post-processing."
               end if
               flag_flex_fermi = .true.
            end if

         case('TDDFT_run')
           read(inputline,*,end=88,err=99) desc_str, TDDFT_time, TDDFT_dt, TDDFT_write_info, TDDFT_output_every
           flag_run_tddft_real_time_propagation = .true.

         case('TDDFT_propagator')
           read(inputline,*,end=88,err=99) desc_str, TDDFT_propagator

           if(TDDFT_propagator.eq.'exact_matrix_exponential') then
             write(info_str,'(2X,2A)') 'Using exact matrix exponential for TDDFT propagation'
             call localorb_info(info_str,use_unit,'(A)')
           else if (TDDFT_propagator.eq.'crank_nicholson') then
             write(info_str,'(2X,2A)') 'Using Crank-Nicholson/Cayley propagator for TDDFT propagation'
             call localorb_info(info_str,use_unit,'(A)')
           else
             write(info_str,'(2X,2A)')'* WARNING: TDDFT propagator not known: ', TDDFT_propagator
             call localorb_info(info_str,use_unit,'(A)')
             call aims_stop_coll('', func)
           end if

         case('MD_run')
            ! Error trap! Must break here if both an MD_schedule and an MD_run are defined.
            ! The two can not be specified together.
            if (MD_use_schedule) then
               write(info_str,'(1X,A)')'* Error! Found BOTH the MD_run and the MD_schedule keywords'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)')'* in the same control.in file.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)')'* You can either specify an MD_run, or an MD_schedule,'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)')'* but not both together - this would be contradictory information.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)')'* Stopping run - please correct your control.in file and run again.'
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if

            read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble
            ! ATTENTION: Any flags that you add here, you should add to the appropriate section of
            ! the MD_schedule parser as well!!!
            flag_molecular_dynamics      = .true.
            final_forces_cleaned         = .true.
            output_in_original_unit_cell = .true.
            if (n_periodic.eq.0) then
               ! Dipole moment for non-periodic calculations only, for
               ! periodic systems we'd have to calculate it differently!
               out_dipole_moment = .true.
            end if
            if (MD_ensemble.eq.'NVE') then
               write(info_str,'(2X,2A)') 'Running Born-Oppenheimer ', &
                     'molecular dynamics in NVE ensemble.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
            else if (MD_ensemble.eq.'NVE_4th_order') then
               write(info_str,'(2X,2A)') 'Running Born-Oppenheimer ', &
                     'molecular dynamics in NVE ensemble.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | Using a fourth-order symplectic integrator'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
            else if (MD_ensemble.eq.'NVE_damped') then
               write(info_str,'(2X,2A)') 'Running Born-Oppenheimer ', &
                     'molecular dynamics in damped NVE ensemble.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               read(inputline,*,end=88,err=99) desc_str,MD_time,MD_ensemble, NVE_damping_factor
               write(info_str,'(2X,A,F14.6)') ' | damping factor = ', &
                     NVE_damping_factor
               call localorb_info(info_str,use_unit,'(A)')
            else if (MD_ensemble.eq.'NVT_berendsen') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, &
                     MD_tau_berendsen
               write(info_str,'(2X,3A)') 'Running Born-Oppenheimer ', &
                     ' molecular dynamics in NVT ensemble with ', &
                     ' Berendsen thermostat'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
                     MD_temperature,' K'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,2A,F14.6,A)') ' | thermostat ', &
                     'relaxation time = ',MD_tau_berendsen,' ps'
               call localorb_info(info_str,use_unit,'(A)')
            else if (MD_ensemble.eq.'NVT_andersen') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, &
                     MD_nu_andersen
               write(info_str,'(2X,3A)') 'Running Born-Oppenheimer ', &
                     ' molecular dynamics in NVT ensemble with ', &
                     ' Andersen thermostat'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
                     MD_temperature,' K'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,2A,F14.6,A)') ' | thermostat ', &
                     'stochastic coupling nu = ',MD_nu_andersen,' ps^-1'
               call localorb_info(info_str,use_unit,'(A)')

            else if (MD_ensemble.eq.'NVT_parrinello') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, &
                     MD_tau_BDP
               write(info_str,'(2X,3A)') 'Running Born-Oppenheimer ', &
                     ' molecular dynamics in NVT ensemble with ', &
                     ' Bussi-Donadio-Parrinello thermostat. '
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
                     MD_temperature,' K'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,2A,F14.6,A)') ' | thermostat ', &
                     'relaxation time tau = ',MD_tau_BDP,' ps'
               call localorb_info(info_str,use_unit,'(A)')

            else if (MD_ensemble.eq.'GLE_thermostat') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, MD_gle_ns!, &
!                     MD_glew
               write(info_str,'(2X,3A)') 'Running Born-Oppenheimer ', &
                     ' molecular dynamics with the Generalized Langevin Eq. ', &
                     ' colored-noise thermostat (Ceriotti et al.).'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(X,A)') 'Character of the thermostat defined below '
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(X,A)') 'WARNING: experimental version. Use at your own risk.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
                     MD_temperature,' K'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,I5)') ' | with extra number of dof  = ', &
                     MD_gle_ns
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)') &
                     "** Please do not use it if you do not understand this implementation."
               call localorb_info(info_str,use_unit,'(A)')
!               call aims_stop_coll('', func)
               call read_gle(MD_temperature, MD_gle_ns)
            else if (MD_ensemble.eq.'NVT_nose-poincare') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, &
                     MD_Q_NP
               write(info_str,'(1X,3A)') '* Error. Requested option: Born-Oppenheimer ', &
                     ' molecular dynamics in NVT ensemble with ', &
                     ' Nose-Poincare thermostat.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A,A)') &
                     "* We do not intend to develop this option any further,",&
                     " as we several other working thermostats are available."
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)') &
                     "* Please choose a different thermostat, for instance the BDP thermostat."
  !             call localorb_info(info_str,use_unit,'(A)')
  !             write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
  !                   MD_time, ' ps'
  !             call localorb_info(info_str,use_unit,'(A)')
  !             write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
  !                   MD_temperature,' K'
  !             call localorb_info(info_str,use_unit,'(A)')
  !             write(info_str,'(2X,2A,F14.6)') ' | thermostat ', &
  !                   'effective mass = ',MD_Q_NP
  !             call localorb_info(info_str,use_unit,'(A)')
   !warning to be removed once MD is tested
               call aims_stop_coll('Nose-Poincare thermostat not supported.', func)
   !end warning to be removed once MD is tested
            else if (MD_ensemble.eq.'NVT_nose-hoover') then
               read(inputline,*,end=88,err=99) desc_str, MD_time, MD_ensemble, MD_temperature, &
                     MD_Q_NP
               write(info_str,'(2X,3A)') 'Running Born-Oppenheimer ', &
                     ' molecular dynamics in NVT ensemble with ', &
                     ' Nose-Hoover thermostat'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | simulation time = ', &
                     MD_time, ' ps'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F14.6,A)') ' | at temperature  = ', &
                     MD_temperature,' K'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,2A,F14.6)') ' | thermostat ', &
                     'effective mass = ',MD_Q_NP
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,2A)')'* WARNING: molecular dynamics ', &
                     'ensemble not known: ', MD_ensemble
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if

         case('PIMD_run')         ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, PIMD_time, PIMD_ensemble, PIMD_temperature, PIMD_Q_NP

         case('normal_mode')      ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, desc_str

         case('PIMD_time_step')   ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, PIMD_tstep

         case('MD_thermostat_units')
            read(inputline,*,end=88,err=99) desc_str, MD_thermostat_units
            if ((MD_thermostat_units.eq.'amu*bohr^2').or.(MD_thermostat_units.eq.'cm^-1')) then
               write(info_str,'(2X,2A)') "Thermostat units for molecular dynamics: ", MD_thermostat_units
               call localorb_info(info_str)
            else
               write(info_str,'(2X,2A)') "Thermostat units unknown: ",MD_thermostat_units
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "Aborting."
               call localorb_info(info_str)
               call aims_stop_coll('', func)
            end if
         case('MD_gle_A')
               write(info_str,'(2X,2A)') "A matrix for GLE thermostat specified"
               call localorb_info(info_str)
         case('MD_gle_C')
               write(info_str,'(2X,2A)') "C matrix for GLE thermostat specified"
               call localorb_info(info_str)
         case('MD_time_step')
            read(inputline,*,end=88,err=99) desc_str, MD_tstep
            write(info_str,'(2X,2A,F14.6,A)') 'Molecular dynamics time ', &
                  'step = ', MD_tstep, ' ps'
            call localorb_info(info_str,use_unit,'(A)')

            flag_MD_tstep = .true.

         case('MD_RNG_seed')
            read(inputline,*,end=88,err=99) desc_str, seed
            MD_RNG_seeded = .true.
            write(info_str,'(2X,A,I8)') "Found seed for MD random number generator : ", seed
            call localorb_info(info_str,use_unit,'(A)')

         case('MD_MB_init')
            read(inputline,*,end=88,err=99) desc_str, MD_init_temperature
            MB_velocity_initialization = .true.
            write(info_str,'(2X,2A,F14.6,A)')'Initializing MD run with ', &
               'Maxwell-Boltzmann momentum distribution at T = ', &
                  MD_init_temperature, ' K'
            call localorb_info(info_str,use_unit,'(A))')

         case('MD_QH_init')
            i_MD_QH_init_segments = i_MD_QH_init_segments + 1
            read(inputline,*,end=88,err=99) desc_str, MD_QH_first_atom (i_MD_QH_init_segments), &
                              & MD_QH_last_atom  (i_MD_QH_init_segments), &
                              & MD_QH_temperature(i_MD_QH_init_segments), &
                              & MD_QH_file       (i_MD_QH_init_segments)
            if (i_MD_QH_init_segments .eq. 1) then
               write(info_str,'(2X,A)')&
               'Parameters for the initialization of the MD run on the basis of the quasi-harmonic approximation:'
               call localorb_info(info_str,use_unit,'(A))')
            end if
            write(info_str,'(2X,A,I5,A,I5,A,F9.3,A,A,A)') &
                                    '| Initializing atom', MD_QH_first_atom (i_MD_QH_init_segments), &
                                    ' to atom',            MD_QH_last_atom  (i_MD_QH_init_segments), &
                                    ' at T =',             MD_QH_temperature(i_MD_QH_init_segments), &
                                    'K on the basis of ',trim(MD_QH_file(i_MD_QH_init_segments)),"."
            call localorb_info(info_str,use_unit,'(A))')

         case('MD_clean_rotations')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.false.') then
               MB_clean_rotations = .false.
               flag_MB_clean_rotations_read = .true.
               write(info_str,'(2X,A)') 'Explicitly not using rotational cleaning for initial velocities.'
               call localorb_info(info_str,use_unit,'(A)')
            else if (desc_str.eq.'.true.') then
               if(MD_ensemble.eq.'GLE_thermostat') then
                   write(info_str,'(2X,A)') '** WARNING: This thermostat does not work with velocity cleaning.'
                   call localorb_info(info_str,use_unit,'(A)')
                   write(info_str,'(2X,A)') 'Please set MD_clean_rotations .false. in your control.in in order to proceed.'
                   call localorb_info(info_str,use_unit,'(A)')
                   call aims_stop_coll('Stopping execution due to incompatibiliy of MD_clean_rotations and thermostat.', func)
               endif
               MB_clean_rotations = .true.
               flag_MB_clean_rotations_read = .true.
               write(info_str,'(2X,A)') 'Cleaning initialized velocities from rotations.'
               call localorb_info(info_str,use_unit,'(A)')
               if ( n_periodic .gt. 0) then
                  write(info_str,'(2X,A)') '*** WARNING: Rotational cleaning should not be used for periodic simulations.'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** For periodic simulations, removing perceived rotations may actually '
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** introduce wrong physics. Since this has led to some confusion in the past,'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** the code now stops for periodic systems if MD_clean_rotations is set.'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** '
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** If you were running with this set of options on purpose (for example,'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** for an isolated molecule in a periodic box), you should re-enable the'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** option by hand, simply by commenting out the stop command in '
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** read_control.f90 and recompiling.'
                  call localorb_info(info_str,use_unit,'(A)')
                  call aims_stop_coll('Stopping execution due to questionable MD_clean_rotations setting.', func)
               end if
            else
               write(info_str,'(1X,A,A,A)') '* Unknown setting ', desc_str, ' for MD_clean_rotations in control.in.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') 'Please use .true. or .false.'
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('Stopping execution due to unknown MD_clean_rotations setting.', func)
            end if

         case('MD_maxsteps')
            read(inputline,*,end=88,err=99) desc_str, MD_maxsteps
            write(info_str,'(2X,A,I8)')'Maximal number of MD steps is ', &
                  MD_maxsteps
            call localorb_info(info_str,use_unit,'(A)')

         case('check_MD_stop')
            read(inputline,*,end=88,err=99) desc_str, check_MD_stop

         case('MD_schedule')
            ! ATTENTION: Anything you add here, you should probably also add in the MD_run section
            ! CC: TDI only shows app in the scheduler section, though
            if (MD_segments.gt.0) then
               flag_molecular_dynamics      = .true.
               flag_MD_schedule             = .true.
               final_forces_cleaned         = .true.
               output_in_original_unit_cell = .true.
               if (n_periodic.eq.0) then
               ! Dipole moment for non-periodic calculations only, for
               ! periodic systems we'd have to calculate it differently!
                  out_dipole_moment = .true.
               end if
               ! CC Begin: Add parsing of thermodynamic_integration
               ! LOOP over all MD_segments, and  read one TDI line for each MD_segment
               do MD_schedule_step = 1, MD_segments
                  read(7,'(A)',iostat=i_code) inputline
                  if(i_code<0)then
                     if (myid.eq.0) then
                        write(use_unit,*) "Error reading file 'control.in':"
                        write(use_unit,*) "EOF reached within 'MD_schedule' block"
                     endif
                     call aims_stop_coll('', func)
                  endif
                  if(i_code>0)then
                     call aims_stop_coll("Unknown error reading file 'control.in'. (expecting MD_segment)", func)
                  endif

                  read(inputline,*,end=88,err=99) &
                  desc_str, MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step)

                  !MD block parsed in here
                  select case (MD_schedule_ensemble(MD_schedule_step))
                  case('NVE')
                     continue
                  case('NVE_4th_order')
                     continue
                  case('NVE_damped')
                     read(inputline,*,end=88,err=99) desc_str, &
                         MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_damping_factor(MD_schedule_step)
                  case('NVT_berendsen')
                     read(inputline,*,end=88,err=99) desc_str, &
                         MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_temperature(MD_schedule_step), MD_schedule_tau_berendsen(MD_schedule_step)
                  case('NVT_nose-poincare')
                     read(inputline,*,end=88,err=99) desc_str, &
                         MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_temperature(MD_schedule_step), MD_schedule_Q(MD_schedule_step)
                  case('NVT_nose-hoover')
                     read(inputline,*,end=88,err=99) desc_str, &
                         MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_temperature(MD_schedule_step), MD_schedule_Q(MD_schedule_step)
                  case('NVT_andersen')
                     read(inputline,*,end=88,err=99) desc_str, &
                         MD_schedule_time(MD_schedule_step), MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_temperature(MD_schedule_step), MD_schedule_nu_andersen(MD_schedule_step)
                  case('NVT_parrinello')
                     read(inputline,*,end=88,err=99) desc_str, MD_schedule_time(MD_schedule_step), &
                         MD_schedule_ensemble(MD_schedule_step), &
                        &MD_schedule_temperature(MD_schedule_step), MD_schedule_tau_BDP(MD_schedule_step) ! , MD_schedule_random_BDP(MD_schedule_step)
                  case default
                          if (myid.eq.0) then
                              write(use_unit,*) '# *** ERROR: '
                              write(use_unit,*) '# *** MD_schedule: unknown or unimplemented ensemble: ', &
                                  trim(MD_schedule_ensemble(MD_schedule_step)),'.'
                              write(use_unit,*) '# *** Aborting execution.'
                          end if
                          call aims_stop_coll('', func)
                      end select
                  ! ^^ MD block parsed here

                  ! TDI block parsed here
                  if (use_thermodynamic_integration .or. use_reffree_AS .or. use_harmonic_pot_only) then
                     read(7,'(A)',iostat=i_code) inputline
                     if(i_code<0) then
                        if (myid.eq.0) then
                           write(use_unit,*) "# *** ERROR reading file 'control.in':"
                           write(use_unit,*) "# *** EOF reached within 'MD_schedule' block"
                        endif
                        call aims_stop_coll('', func)
                     endif
                     if(i_code>0) then
                        call aims_stop_coll&
                        ("Unknown error reading file 'control.in'. (expecting 'thermodynamic_integration')", func)
                     endif

                     read(inputline,*,end=88,err=99) desc_str
                     if((desc_str /= 'thermodynamic_integration') &
                         .and. &
                        (desc_str /= 'adiabatic_switching') &
                         .and. &
                        (desc_str /=  'harmonic_potential_only') )&
                     then
                        if (myid.eq.0) then
                           write(use_unit,*) '# *** ERROR: '
                           write(use_unit,*) '# *** MD_schedule: No thermodynamic integration specified for MD_segement: ', &
                                      MD_schedule_step,'.'
                        endif
                        call aims_stop_coll('', func)
                     endif
                  endif

                  if(use_thermodynamic_integration .and. (.not. use_harmonic_pot_only) ) then
                    read(inputline,*,end=88,err=99) desc_str,TDI_lambda_start(MD_schedule_step),TDI_lambda_end(MD_schedule_step),&
                              &TDI_QHA_file(MD_schedule_step),TDI_QHA_E0(MD_schedule_step)
                    TDI_QHA_E0(MD_schedule_step) = TDI_QHA_E0(MD_schedule_step)/hartree

                  else if (use_reffree_AS) then
                    read(inputline,*,end=88,err=99) desc_str,TDI_lambda_start(MD_schedule_step),TDI_lambda_end(MD_schedule_step)
                     !TDI_QHA_file(MD_schedule_step) = ''
                     !TDI_QHA_E0(MD_schedule_step) = 0.d0
                  else if (use_harmonic_pot_only) then
                    read(inputline,*,end=88,err=99) desc_str, TDI_QHA_file(MD_schedule_step)
                    TDI_lambda_start(MD_schedule_step) = 0.0d0
                    TDI_lambda_end(MD_schedule_step)   = 0.0d0
                    TDI_QHA_E0(MD_schedule_step)       = 0.0d0
                    skip_SCF                           = .true.
                  endif
                  ! ^^ TDI block parsed here
               end do
   ! CC End: Add parsing of thermodynamic_integration
            else ! this is the case if there are no MD segments at all
               write(info_str,'(1X,A)')'* Warning! Found MD_schedule keyword, but no MD_segments defined!'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)')'* Ignoring MD_schedule keyword!'
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('skip_SCF')
            read(inputline,*,end=88,err=99) desc_str, SKIP_SCF
            if (skip_SCF) then
               write(info_str,'(2X,A)') ' * Warning: Skipping the SCF cycles for debug purposes.'
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('MD_restart_binary')
            read(inputline,*,end=88,err=99) desc_str, MD_restart_binary
            if (MD_restart_binary) then
               write(info_str,'(2X,A)') ' * Warning: Using obsolete binary format for the MD restart file.'
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('MD_restart')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'time_restart') then
               write(info_str,'(2X,2A)') 'Reading MD restart information from file', &
                     trim(MD_restart_file)
               call localorb_info(info_str,use_unit,'(A)')
               MD_initialconf_from_restart =  .true.
            else if (desc_str.eq.'.false.') then
               write(info_str,'(2X,A)') "Explicitly requesting no molecular dynamics restart."
               call localorb_info(info_str,use_unit,'(A)')
               MD_initialconf_from_restart =  .false.
            else if (desc_str.eq.'.true.') then
!               if (MD_ensemble .eq. 'GLE_thermostat') then
!                 call aims_stop_coll('Restart for this thermostat is under implementation. Set MD_restart .false. and run again.', func)
!               endif
               MD_time_restart = .false.
               MD_initialconf_from_restart = .true.
               write(info_str,'(2X,A)') 'Continuing molecular dynamics run from previous calculation'
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,A)') 'This keyword for the restart information does not exist. Please refer to the manual.'
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if

         case('PIMD_restart')       ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               call localorb_info(info_str,use_unit,'(A)')
               PIMD_initialconf_from_restart = .true.
            end if

         case('plumed')
!           read(inputline,*,end=88,err=99) desc_str, plumed_plugin, plumed_file
            read(inputline,*,end=88,err=99) desc_str, desc_str
!           read(inputline,*,end=88,err=99) desc_str, plumed_plugin
            if (desc_str.eq.'.true.') then
               plumed_plugin = .true.
               write(info_str,'(2X,2A)') 'Using PLUMED, writing on file: ',plumed_file
!           else
!              plumed_plugin = .false.
!              write(info_str,'(2X,A)') ''

            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('plumed_new')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'.true.') then
               plumed_new_plugin = .true.
               write(info_str,'(2X,2A)') 'Using PLUMED(NEW IMPLEMENTATION), writing on file: ',plumed_file
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('use_pimd_wrapper')
            write(info_str,'(2X,A)') 'Using external wrapper (i-PI) for performing (path integral) molecular dynamics '
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') '**Attention: initial geometry.in file will be ignored!'
            call localorb_info(info_str,use_unit,'(A)')
            read(inputline,*,end=88,err=99) desc_str, pimdwrapper_server, pimdwrapper_port

         case('communicate_pimd_wrapper')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str .eq.'dipole') then
                write(info_str,'(2X,A)') 'Will pass dipole to PIMD wrapper i-PI'
                call localorb_info(info_str,use_unit,'(A)')
                out_dipole_moment = .true.
                ipi_dip = .true.
            !else if (desc_str .eq.'z-dipole') then ! MR: need to think how to warn user about the type of dipole correction, etc.
            !   write(info_str,'(2X,A)') 'Will pass z-component of dipole (periodic slab calculations) to PIMD wrapper i-PI'
            !    call localorb_info(info_str,use_unit,'(A)')
            !    use_dipole_correction = .true.
            !    dipole_correction_method = 'dipole'
            !    ipi_dip = .true.
            else if (desc_str .eq.'hirshfeld') then
                write(info_str,'(2X,A)') 'Will pass hirshfeld charge to PIMD wrapper i-PI'
                call localorb_info(info_str,use_unit,'(A)')
                out_hirshfeld_always = .true.
                out_hirshfeld = .true.
                ipi_hirshfeld = .true.
            else if (desc_str .eq.'workfunction') then
                write(info_str,'(2X,A)') 'Will pass converged work functions to PIMD wrapper i-PI'
                call localorb_info(info_str,use_unit,'(A)')
                calculate_work_function = .true.
                ipi_work = .true.
            else if (desc_str .eq.'friction') then
                write(info_str,'(2X,A)') 'Will pass non-adiabatic friction to PIMD wrapper i-PI'
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') '** WARNING: Here using numerical friction by default'
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,A)') '** WARNING: All parameters for friction calculation will have their default values unless specifically set in control.in'
                call localorb_info(info_str,use_unit,'(A)')
                use_friction = .true.
                numerical_friction = .true.
                ipi_ftensor = .true.

            ! FlK: communicate stress per atom to ipi
            else if (desc_str .eq. 'atomic_stress') then
                write(info_str,'(2X,A)') 'Will pass atomic stress to PIMD wrapper i-PI'
                call localorb_info(info_str,use_unit,'(A)')
                use_analytical_stress = .true.
                compute_heat_flux     = .true.
                ipi_atomic_stress     = .true.

            else
               if (myid.eq.0) then
                  write(use_unit,*) '"', desc_str, '"', " is not a valid communicating", &
                        " identifier for i-PI. "
               end if

               call aims_stop_coll('', func)
            end if

         case('calc_dens_superpos')
            read(inputline,*,end=88,err=99) desc_str, recalc_dens

         case('use_potential_constrain')
            read(inputline,*,end=88,err=99) desc_str, potconstype, potconspar
            use_potential_constrain=.true.
            ! MR TODO: write a catch here if someone forgets to write radius of sphere

         case('output_boys_centers')
            if (myid .eq. 0) then
                write (use_unit,'(2X,A)') &
                  "Compute Dipole Matrix and output maximally localized Boys centers for a given configuration"
            end if
            k_point_mom=1
            if (n_periodic.gt.0) then
               if (myid .eq. 0) write (use_unit,'(2X,A)') &
                  "** STOP: For periodic systems Wannier centers are needed. Under implementation... "
               call aims_stop_coll('No periodic Boys', func)
            endif
            if (n_spin.gt.1) call aims_stop_coll('Sorry, still need to add spin dimension to Boys', func)
            flag_out_boys = .true.
            use_symmetry_reduced_k_grid = .false.
            if (use_scalapack) then
              collect_eigenvectors = .true.
            endif

         case('wf_extrapolation')
            if (n_wf_funcs /= 0) then
               call aims_stop_coll('Double wf_extrapolation spec', func)
            end if
            flag_wf_named = .true.
            if (flag_wf_func) then
               call aims_stop_coll('Cannot use both wf_extrapolation <name> and wf_func ...', func)
            end if
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
            case('linear')
               wf_extra_type = WF_EXTRA_FUNC
               n_wf_funcs = 2
               call parse_general_function('constant', wf_funcs(1))
               call parse_general_function('linear', wf_funcs(2))
            case('quadratic')
               wf_extra_type = WF_EXTRA_FUNC
               n_wf_funcs = 3
               call parse_general_function('constant', wf_funcs(1))
               call parse_general_function('linear', wf_funcs(2))
               call parse_general_function('quadratic', wf_funcs(3))
            case('polynomial')
               wf_extra_type = WF_EXTRA_FUNC
               read(inputline,*,end=88,err=99) desc_str, desc_str, n_wf_funcs, order
               if (order+1 > n_wf_funcs) then
                  call aims_stop_coll('wf order too large', func)
               end if
               do i_ord = 0, order
                  write(info_str, "('polynomial ',I20)") i_ord
                  call parse_general_function(info_str, wf_funcs(i_ord+1))
               end do
               do i_ord = order+1, n_wf_funcs-1
                  ! add next odd order
                  write(info_str, "('polynomial ',I20)") &
                  & order + 1 + mod(order, 2) + 2*(i_ord - order - 1)
                  call parse_general_function(info_str, wf_funcs(i_ord+1))
               end do
            case('niklasson06')
               wf_extra_type = WF_EXTRA_NIKL
               read(inputline,*,end=88,err=99) desc_str, desc_str, n_wf_funcs
               n_indep = (n_wf_funcs+1) / 2
               order = 2*(n_indep-1) + 1
            case('none')
               wf_extra_type = WF_EXTRA_NONE
               write(info_str,'(2X,A)') &
               & 'No wf extrapolation (i.e. just use the last one)'
               call localorb_info(info_str)
            case default
               call aims_stop_coll('Type of wf_extrapolation not supported', func)
            end select

            select case (wf_extra_type)
            case(WF_EXTRA_NIKL)
               write(info_str,"(2X,3A,I3,A,I4,A)") &
               & 'Using niklasson06 wf extrapolation of order', &
               & order, ' and with up to', n_wf_funcs, ' points.'
               call localorb_info(info_str)
            case(WF_EXTRA_FUNC)
               do i_ord = 1, n_wf_funcs
                  call format_general_function(wf_funcs(i_ord), info_str)
                  write(info_str, "(2X,'Added implicit wf_func ',A)") &
                  & trim(info_str)
                  call localorb_info(info_str)
               end do
            case(WF_EXTRA_NONE)
               call localorb_info('  Wave function extrapolation disabled.')
            case default
            call aims_stop('Internal wf_extra_type error', func)
            end select

         case('wf_func')
            flag_wf_func = .true.
            if (flag_wf_named) then
               call aims_stop_coll('Cannot use both wf_extrapolation <name> and wf_func ...', func)
            end if
            wf_extra_type = WF_EXTRA_FUNC
            n_wf_funcs = n_wf_funcs + 1
            if (n_wf_funcs > n_max_wf_funcs) then
               call aims_stop_coll('Too many wf_func specs', func)
            end if
            read(inputline,"(A)") buffer
            buffer = adjustl(buffer)
            i_space = scan(buffer, ', ')
            info_str = buffer(i_space+1:)
            call parse_general_function(info_str, wf_funcs(n_wf_funcs))
            call format_general_function(wf_funcs(n_wf_funcs), info_str)
            write(info_str, "(2X,'Added explicit wf_func ',A)") trim(info_str)
            call localorb_info(info_str)

         case('wf_extra_use_densmat')
            wf_extra_use_densmat = .true.
               write(info_str,'(2X,A)') "Using extrapolated density matrix."
               call localorb_info(info_str,use_unit,'(A)')

         case('orthonormalize_eigenvectors')
            read(inputline,*,end=88,err=99) desc_str, orthonormalize_evs
            if (orthonormalize_evs) then
               write(info_str,'(2X,2A)') "Orthonormalizing eigenvectors ", &
                  "after a relaxation step."
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,2A)') "Explicitly not orthonormalizing eigenvectors ", &
                  "after a relaxation step."
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('energy_tolerance')
            read(inputline,*,end=88,err=99) desc_str, energy_tolerance
            energy_tolerance = energy_tolerance/hartree
            write(info_str,'(2A,E14.6,A)') '  Max energy increase ', &
            & 'before discarding BFGS move: ',energy_tolerance*hartree,' eV'
            call localorb_info(info_str,use_unit,'(A)')
            flag_energy_tolerance = .true.

         case('aggregated_energy_tolerance')
            read(inputline,*,end=88,err=99) desc_str, aggregated_energy_tolerance
            aggregated_energy_tolerance = aggregated_energy_tolerance/hartree
            write(info_str,'(2A,E14.6,A)') '  Max energy increase ', &
            & 'across multiple BFGS steps: ',aggregated_energy_tolerance*hartree,' eV'
            call localorb_info(info_str,use_unit,'(A)')
            flag_aggregated_energy_tolerance = .true.

         case('harmonic_length_scale')
            read(inputline,*,end=88,err=99) desc_str, harmonic_length_scale
            harmonic_length_scale = harmonic_length_scale / bohr
            write(info_str,'(A,E14.6,A)') &
            & '  Assuming qualitatively harmonic for steps shorter than ', &
            & harmonic_length_scale*bohr, ' A.'
            call localorb_info(info_str,use_unit,'(A)')
            flag_harmonic_length_scale = .true.

         case('max_atomic_move')
            read(inputline,*,end=88,err=99) desc_str, max_atomic_move
            max_atomic_move = max_atomic_move/bohr
            write(info_str,'(2A,E14.6, A)')'  Maximally allowed move', &
            ' during BFGS: ',max_atomic_move*bohr,' A'
            call localorb_info(info_str,use_unit,'(A)')
            flag_max_atomic_move = .true.

         case('min_line_step')
            read(inputline,*,end=88,err=99) desc_str, min_line_step
            write(info_str,'(2A,E14.6,A)')'  Minimum line step ', &
               ' after reduction before reinitializing Hessian in BFGS: ', &
               min_line_step,' A'
            call localorb_info(info_str,use_unit,'(A)')
            flag_min_line_step = .true.

         case('line_step_reduce')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'automatic') then
               flag_line_step_reduce = .true.
               write(info_str,'(2X,2A)') 'Line minimization in BFGS ', &
               'quadratic interpolation after rejected steps.'
            else
               read(inputline,*,end=88,err=99) desc_str, line_step_reduce
               write(info_str,'(2A,E14.6)')'  Factor by which the trial ', &
                  'line step will be reduced after rejected BFGS move: ', &
                     line_step_reduce
               call localorb_info(info_str,use_unit,'(A)')
               flag_line_step_reduce = .true.
               line_step_reduce_automatic = .false.
            end if

         case('max_relaxation_steps')
            read(inputline,*,end=88,err=99) desc_str, max_relaxation_steps
            write(info_str,'(2A,I8)') '  Maximum steps for geometry ', &
               'relaxation: ', max_relaxation_steps
            call localorb_info(info_str,use_unit,'(A)')
            flag_max_relaxation_steps = .true.

         case('init_hess')
            read(inputline,"(A)") buffer
            read(buffer, *) desc_str, desc_str
            select case (desc_str)
            case('lindh', 'Lindh')
               init_hess_type = HESS_LINDH
               ! Convert default to eV/A^2 to match control.in convention; if
               ! no parameter is given, read(..., iostat=info) doesn't change
               ! it.
               init_hess_lindh_diag = init_hess_lindh_diag * (hartree/bohr**2)
               read(buffer,*, iostat=info) desc_str, desc_str, &
               &         init_hess_lindh_diag, init_hess_lindh_thres
               init_hess_lindh_diag = init_hess_lindh_diag / (hartree/bohr**2)
               ! In case of trouble, leave parameters on default.
               write(info_str,'(2X,A,F6.2,A,E14.6,A)') &
               & 'Initial Hessian is Lindh matrix (thres =', &
               init_hess_lindh_thres, ') plus', &
               & init_hess_lindh_diag * (hartree/bohr**2), &
               & ' eV/A^2 times unity.'
               call localorb_info(info_str)
            case('reciprocal_lattice')
               write(info_str,'(2X, 3A)') &
                   & 'WARNING: Initializing the Hessian with', &
                   & ' the reciprocal lattice is requested. This is', &
                   & ' an experimental feature, please be careful!!'
               call localorb_info(info_str)

               init_hess_type = HESS_LATTICE
               ! Use reciprocal lattice to initialize Hessian for the
               ! lattice vectors, use Lindh model to initialize Hessian
               ! for the atomic positions as before
               ! Convert default to eV/A^2 to match control.in convention; if
               ! no parameter is given, read(..., iostat=info) doesn't change
               ! it.
               init_hess_lindh_diag = init_hess_lindh_diag * (hartree/bohr**2)
               ! In case of trouble, leave parameters on default.
               write(info_str,'(2X,A,E14.6,3A)') &
               & 'Initial Hessian is diagonal matrix with', &
               & init_hess_diag * (hartree/bohr**2), &
               & ' eV/A^2 times unity for atomic positions.', &
               & ' For the lattice, the reciprocal lattice matrix',&
               & ' is used.'
               call localorb_info(info_str)
            case default   ! esp. diag

               init_hess_type = HESS_DIAG

               if (desc_str == 'diag') then
                  read(buffer,*, iostat=info) desc_str, desc_str, init_hess_diag
                  ! In case of trouble, init_hess_diag is left on default
               else
                  ! Compatibility: 'init_hess 30.'
                  read(buffer,*, iostat=info) desc_str, init_hess_diag
                  if (info /= 0) then
                     call aims_stop_coll('Trouble with init_hess line', func)
                  else
                     write(info_str, "('* ',A,A,F8.4,A,A)") &
                     & 'Deprecated syntax for init_hess; ', &
                     & 'use "init_hess diag ', init_hess_diag, '" ', &
                     & 'or even try "init_hess Lindh".'
                     call localorb_info(info_str, use_unit, "(2X,A)")
                  end if
               end if
               init_hess_diag = init_hess_diag / (hartree/bohr**2)
               write(info_str,'(2X,A,E14.6,A)') &
               & 'Initial Hessian is ', init_hess_diag * (hartree/bohr**2), &
               & ' eV/A^2 times unity.'
               call localorb_info(info_str)
            end select
            if ((n_hess_blocks > 0)        .or. &
               (n_hess_blocks_lv > 0)      .or. &
               (n_hess_blocks_lv_atom > 0) .or. &
               hess_in_file) then
               init_hess_type = HESS_EXPL
               write(info_str, "(2X,'*** ',A,' ',A)") &
               & 'Specified both init_hess in', &
               & 'control.in and explicit Hessian in geometry.in.'
               call localorb_info(info_str)
               write(info_str, "(2X,'*** ',A)") &
               & 'Using the latter for Hessian initialization.'
               call localorb_info(info_str)
            end if
            flag_init_hess = .true.

         case('init_hess_lv_diag')
            read(inputline,*,end=88,err=99) desc_str, init_hess_lv_diag
            init_hess_lv_diag = init_hess_lv_diag / (hartree/bohr**2)

            if (init_hess_lv_diag.gt.0.d0) then
              write(info_str,'(2X,A,F15.8,A)') &
                 'Initial Hessian matrix for relaxation: Will use ', &
                 init_hess_lv_diag*(hartree/bohr**2), ' eV/A^2 as initial diagonal element for lattice vectors.'
              call localorb_info(info_str,use_unit,'(A)')
            else
              write(info_str,'(1X,A,F15.8,A)') &
                 '* Illegal value for init_hess_lv_diag: ', &
                 init_hess_lv_diag*(hartree/bohr**2), ' eV/A^2. Must be greater than zero (25 is a reasonable magnitude).'
              call localorb_info(info_str,use_unit,'(A)')
              call aims_stop_coll('', func)
            end if

         case('line_search_tol')
            read(inputline,*,end=88,err=99) desc_str, line_search_tol
            write(info_str,'(2X,2A,E14.6)') 'Automatic line search ', &
               'tolerance for BFGS relaxation: ', line_search_tol
            call localorb_info(info_str)

         case('bfgs_step')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'quadratic') then
               bfgs_extrapolation_on = .false.
               write(info_str,'(2X,2A)') 'Explicitly requesting ', &
               'use of purely parabolic energy surfaces for BFGS'
               call localorb_info(info_str)
            else
               read(inputline,*,end=88,err=99) desc_str, desc_str, bfgs_mixing_factor, &
                     bfgs_mixing_cap
               bfgs_extrapolation_on = .true.
               write(info_str,'(2X,2A,E14.6)') 'Using BFGS line step ', &
                     ' mixing with mixing factor ', bfgs_mixing_factor
               call localorb_info(info_str)
               write(info_str,'(2X,2A,E14.6)') '| and maximal mixing ', &
                     ' step ', bfgs_mixing_cap
            end if

         case('bfgs_extrapolation_cap')
            read(inputline,*,end=88,err=99) desc_str, bfgs_extrapolation_cap
            write(info_str,'(2X,2A,E14.6)') 'Capping all BFGS line ', &
                  'steps at ', bfgs_extrapolation_cap
               call localorb_info(info_str)

         case('walltime')
            read(inputline,*,end=88,err=99) desc_str, time_read_buffer
            walltime_requested = .true.
            walltime_limit = nint(time_read_buffer)
            write(info_str,'(2X,A,I8,A)') 'Requesting a walltime of ', &
                  walltime_limit, ' seconds.'
            call localorb_info(info_str,use_unit,'(A)')

         case('store_EV_to_disk_in_relaxation')
            read(inputline,*,end=88,err=99) desc_str, store_eigenvectors_to_disk_in_relaxation
            if (store_eigenvectors_to_disk_in_relaxation) then
               write(info_str,'(A,A)') 'Storing eigenvectors to disk ', &
               'in relaxation.'
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(A,A)') 'Storing eigenvectors in memory ', &
               'in relaxation.'
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('preconditioner')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'kerker') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               if (desc_str.eq.'off') then
                  write(info_str,'(2X,2A)') 'Explicitly turning', &
                     ' off kerker preconditioner'
                  call localorb_info(info_str)
                  use_kerker_preconditioner = .false.
               else
                  read(desc_str,*) precondition_kerker_q0
                  write(info_str,'(2X,A,E14.4)') &
                     'Using Kerker preconditioning with damping at q_0 = ', &
                     precondition_kerker_q0
                  call localorb_info(info_str)
                  use_kerker_preconditioner = .true.
               end if
               flag_preconditioner = .true.
            ! dielectric preconditioner
            else if (desc_str.eq.'dielectric') then
               use_dielectric_preconditioner = .true.
               write(info_str,'(2X,A)') &
                     'Using dielectric preconditioning'
               call localorb_info(info_str)
            else if (desc_str.eq.'turnoff') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               if (desc_str.eq.'charge') then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                        preconditioner_turnoff_charge
                  use_preconditioner_turnoff_charge = .true.
                  write(info_str,'(2X,2A,E14.6)')'Turning off ', &
                     'preconditioner after charge convergence of ', &
                     preconditioner_turnoff_charge
                  call localorb_info(info_str)
                  flag_preconditioner_turnoff = .true.
               else if (desc_str.eq.'energy') then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                        preconditioner_turnoff_energy
                  write(info_str,'(2X,2A,E14.6,A)')'Turning off ', &
                     'preconditioner after energy convergence of ', &
                     preconditioner_turnoff_energy, ' eV.'
                  call localorb_info(info_str)
                  flag_preconditioner_turnoff = .true.
               else if (desc_str.eq.'sum_ev') then
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                        preconditioner_turnoff_sum_ev
                  write(info_str,'(2X,2A,E14.6,A)')'Turning off ', &
                     'preconditioner after eigenvalue convergence of ', &
                     preconditioner_turnoff_sum_ev, ' eV.'
                  call localorb_info(info_str)
                  flag_preconditioner_turnoff = .true.
               else
                  write(info_str,'(2X,3A)') &
                        '* WARNING! Precondition turn off mechanism ', &
                        desc_str, ' not known.'
                  call localorb_info(info_str)
                  write(info_str,'(2X,A)') '* Aborting execution.'
                  call localorb_info(info_str)
                  call aims_stop_coll('', func)
               end if
            else if ( (desc_str.eq.'none') .or. (desc_str.eq.'off') ) then
               write(info_str,'(2X,2A)') 'Explicitly turning', &
                   ' off the kerker preconditioner.'
               call localorb_info(info_str)
               use_kerker_preconditioner = .false.
               flag_preconditioner = .true.
            else
               write(info_str,'(2X,3A)') &
                     '* WARNING! Preconditioner ',desc_str,' not known!'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') &
                  '*          fatal error, exiting'
               call localorb_info(info_str)
               call aims_stop_coll('', func)
            end if

         case('kerker_factor')
            read(inputline,*,end=88,err=99) desc_str, kerker_coeff
            write(info_str,'(2X,2A,E14.4,A)') 'Additional empirical factor ', &
               'for kerker preconditioner= ', &
               kerker_coeff, '.'
            call localorb_info (info_str)
            write(info_str,'(2X,2A)') 'Should be the same as damping at ', &
               'q_0 for kerker preconditioner.'
            call localorb_info (info_str)
            use_kerker_factor = .true.

         case('precondition_before_mixer')
            read(inputline,*,end=88,err=99) desc_str, precondition_before_mixer
            if (precondition_before_mixer) then
               call localorb_info('  Precondition before mixing')
               ! JW: If no one opts for this feature for a year of so from now
               ! (2011-06-23), I'd vote to remove this again.
               call localorb_info('  * Note that this is experimental and not well tested')
               call localorb_info('  * Please let us know about any experience you make.')
               call localorb_info('  * Be itgood, bad, or even indifferent.')
            else
               call localorb_info('  Preconditioning after mixing explicitly requested')
            end if

         case('precondition_max_l')
            read(inputline,*,end=88,err=99) desc_str, precondition_max_l
            write(info_str,'(2X,2A,I3)') 'Preconditioning uses ', &
               'multipole expansion up to l = ', &
               precondition_max_l
            call localorb_info (info_str)
            flag_preconditioner_max_l = .true.

         case('clean_forces')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            select case (desc_str)

            case ('none')
               remove_unitary_force_components = 0
               write(info_str,'(2X,A,A)') &
                  "Simple unitary transformations", &
                  " will not be removed from forces."
               call localorb_info ( info_str )

            case ('sayvetz')
               remove_unitary_force_components = 2
               write(info_str,'(2X,A,A)') &
                  "Cleaning forces by fulfilling the ", &
                  "sayvetz conditions."
               call localorb_info ( info_str )

            case default
               write (info_str,'(1X,A,A,A)') &
                  "* clean forces: ", desc_str, &
                  " is not defined as a cleaning method."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)

            end select

            flag_remove_unitary_force_component = .true.

         case('hydro_cut')
            ! This is a flag for testing only, which allows to use purely analytical
            ! hydrogen-like functions.
            ! If this is seriously used by anyone, it should become species-dependent,
            ! and move into read_species_data!


            read(inputline,*,end=88,err=99) desc_str, hydro_cut

            if (myid.eq.0) then
               write (use_unit,'(2X,A,L6)') "Cutoff for hydrogenic fns? ", &
                  hydro_cut
            end if

         case ('xc_pre')
             read (inputline, *, end=88, err=99) desc_str, desc_str, n_steps_xc_pre
             select case (desc_str)
             case ('m06-l')
                 flag_xc_pre = 25
             case ('pbe')
                 flag_xc_pre = 6
             end select

         case('xc')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."pz-lda") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "XC: Using Perdew-Zunger parametrisation ", &
                     "of Ceperley-Alder LDA."
               end if

               flag_xc = 3
            else if (desc_str.eq."pw-lda") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "XC: Using Perdew-Wang parametrisation ", &
                     "of Ceperley-Alder LDA."
               end if
               flag_xc = 8
            else if (desc_str.eq."vwn") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "XC: Using VWN-LDA parametrisation ", &
                     "of VWN5 form."
               end if
               flag_xc = 15
            else if (desc_str.eq."vwn-gauss") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "XC: Using VWN-LDA parametrisation ", &
                     "of VWN-RPA form."
               end if
               flag_xc = 16

            else if (desc_str.eq."hse06".or.desc_str.eq."HSE06") then

               read(inputline,*,end=188,err=99) desc_str,desc_str, hse_omega

               flag_xc = 7
               flag_hse03 = .false.

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.5d-1
               end if

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: Using HSE-functional with OMEGA =",hse_omega," <units>."
                  write (use_unit,'(2X,A)') &
                     "XC: For the actual units (bohr^-1 or Angstrom^-1), see below."
               end if

            else if (desc_str.eq."hse03".or.desc_str.eq."HSE03") then

               flag_xc = 7
               flag_hse03 = .true.
               hse_omega_pbe =0.15d0*2.0D0**(1.0d0/3.0d0)
               hse_omega_hf = 0.15d0/sqrt(2.0d0)

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.5d-1
               end if

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC:  HSE with OMEGA_PBE =" &
                     ,hse_omega_pbe, ' bohr^-1 '

                  write (use_unit,'(2X,A,E14.6,A)') &
                     "and OMEGA_HF = ",hse_omega_hf," bohr^-1."
                  write (use_unit,'(2X,A)') &
                     "This choice follows the 2006 erratum to the original 2003 Heyd, Scuseria and Ernzerhoff reference."
               end if

               ! The parameters above are fixed. No need to specify.
               flag_hse_unit = .true.

             else if (desc_str.eq."lc_wpbeh".or.desc_str.eq."LC_wPBEh") then

               read(inputline,*,end=189,err=99) desc_str,desc_str, hse_omega, hybrid_coeff

               flag_xc = 23

               if (hybrid_coeff.ne.0.0d0) then
                  flag_hybrid_coeff = .true.
               end if

               if (myid.eq.0.and.(.not.flag_hybrid_coeff)) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: Using LC-wPBE with OMEGA =",hse_omega," <units>."
                  write (use_unit,'(2X,A)') &
                     "XC: For the actual units (bohr^-1 or Angstrom^-1), see below."
               else if (myid.eq.0.and.flag_hybrid_coeff) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: Using LC-wPBEh with OMEGA =",hse_omega," <units>"
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: and alpha =",hybrid_coeff,"."
                  write (use_unit,'(2X,A)') &
                     "XC: For the actual units (bohr^-1 or Angstrom^-1), see below."
               end if

            else if (desc_str.eq."pbe") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using PBE gradient-corrected functionals."
               end if

               flag_xc = 6
           else if (desc_str.eq."r48pbe") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using mixed PBE (0.52) and RPBE(0.48) gradient-corrected functionals."
               end if

               flag_xc = 18
           else if (desc_str.eq."pw91_gga") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using PW91 gradient-corrected functionals."
               end if

               flag_xc = 4

            else if (desc_str.eq."pbe_vdw") then  !SAG

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                       "XC: Using PBE with VDW."
               end if

               flag_xc = 5
               use_vdw_method = .true.

            else if (desc_str.eq."pbesol".or.desc_str.eq."PBESOL") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using PBEsol gradient-corrected functionals."
               end if

               flag_xc = 17

            else if (desc_str.eq."pbeint") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using PBEint gradient-corrected functional."
                  write(use_unit,'(2X,A)') &
                     "Reference: E. Fabiano, L.A. Constantin, F. Della Sala, Phys. Rev. B 82, 113104 (2010)"
               end if

               flag_xc = 21

            else if (desc_str.eq."rpbe") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using RPBE gradient-corrected functionals."
               end if

               flag_xc = 11
            else if (desc_str.eq."revpbe") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using revPBE gradient-corrected functionals."
               end if

               flag_xc = 12

            else if (desc_str.eq."revpbe_vdw") then  !SAG

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                       "XC: Using revPBE with vdw."
               end if

               flag_xc = 14
               use_vdw_method = .true.

            else if (desc_str.eq."am05") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using AM05 gradient-corrected functionals."
               end if

               flag_xc = 20

            else if (desc_str.eq."xpbe") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using xPBE gradient-corrected functionals."
                  write(use_unit,'(2X,A)') &
                     "Reference: X. Xu, W.A. Goddard III, J. Chem. Phys. 9, 4068 (2004)"
               end if

               flag_xc = 22

            else if (desc_str.eq."pbe0" .or. desc_str.eq.'PBE0') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using hybrid-PBE0 functionals."
               end if

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.5d-1
               end if

               flag_xc = 1

            else if (desc_str.eq.'hartree' .or. desc_str.eq.'Hartree')then
            !FC: it sets the fock (exchange) operator to zero, and run an Hartree calculation
               if (myid.eq.0) then
                  write(use_unit,*)
                  write(use_unit,'(2X,A)') &
                     "Hartree calculation starts ...... "
                  write(use_unit,'(2X,2A)') &
                     "Attention: this version of Hartree", &
                     " employs an auxiliary basis. "
                  write(use_unit,*)
               end if
               flag_xc = 0
               hybrid_coeff = 0.d0

            else if (desc_str.eq."HF".or.desc_str.eq."hf" ) then

               if (myid.eq.0) then
                  write(use_unit,*)
                  write(use_unit,'(2X,A)') &
                     "Hartree-Fock calculation starts ...... "
                  write(use_unit,'(2X,2A)') &
                     "Attention: this version of Hartree-Fock", &
                     " employs an auxiliary basis. "
                  write(use_unit,*)
               end if

               flag_xc = 0

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 1.0d0
               end if

            else if (desc_str.eq."mp2".or.desc_str.eq."MP2") then

               if (myid.eq.0) then
                  write(use_unit,*)
                  write(use_unit,'(2X,A)') &
                     "Hartree-Fock calculation starts ...... "
                  write(use_unit,'(2X,A)') &
                     "MP2 will start after the HF calculation. "
                  write(use_unit,'(2X,A)') &
                     "Attention: this is just the experimental version."
                  write(use_unit,'(2X,A)') &
                     "Be careful with the results you obtain! "
                  write(use_unit,*)
               endif

               flag_xc = 0
            else if (desc_str.eq."blyp") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                        "XC: Using BLYP functional."
               end if
               flag_xc = 9
            else if (desc_str.eq."b3lyp") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                        "XC: Using hybrid B3LYP functional."
                  write(use_unit,'(2X,A)') &
                        "| This version of B3LYP uses the RPA parametrization of LDA given by"
                  write(use_unit,'(2X,A)') &
                        "| Vosko, Wilk, Nusair. This choice differs from Becke's original paper,"
                  write(use_unit,'(2X,A)') &
                        "| and goes back to a now de facto standard (originally, erroneous)"
                  write(use_unit,'(2X,A)') &
                        "| implementation in the Gaussian suite of codes."
               end if
               flag_xc = 10
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.2d0
               end if
            else if (desc_str.eq."b1lyp") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                        "XC: Using hybrid B1LYP functional."
               end if
               flag_xc = 29
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.5d-1
               end if
            else if (desc_str.eq."pbesol0") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                        "XC: Using hybrid-PBEsol0 functionals."
               end if

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.5d-1
               end if

               flag_xc = 13
            else if (desc_str.eq."screx") then
               if ( myid.eq.0 ) then
                  write(use_unit,'(2X,A)') &
                     "XC: Running screened exchange calculation ..."
               endif
               flag_xc = 0
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 1.d0
               end if
           else if (desc_str.eq."m06-l" .or. desc_str.eq."M06-L") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using M06-L gradient-corrected functionals."
               end if

               flag_xc = 25

           else if (desc_str.eq."m06-hf" .or. desc_str.eq."M06-HF") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using Hybrid M06-HF gradient-corrected functionals."
               end if

               flag_xc = 26

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 1.d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq."m06" .or. desc_str.eq."M06") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using Hybrid M06 gradient-corrected functionals."
               end if

               flag_xc = 27

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 2.7d-1
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq."m06-2x" .or. desc_str.eq."M06-2X") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using Hybrid M06-2X gradient-corrected functionals."
               end if

               flag_xc = 28

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 5.4d-1
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq."tpss" .or. desc_str.eq."TPSS") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using TPSS gradient-corrected functionals."
               end if

               flag_xc = 51

           else if (desc_str.eq."revtpss" .or. desc_str.eq."REVTPSS") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using revTPSS gradient-corrected functionals."
               end if

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "You should be aware that this functional has not been comparatively tested ", &
                  "against equivalent implementations in other programs, as was performed for TPSS."
                  write(use_unit,'(2X,A,A)') &
                  "HOWEVER, it is implemented through the TPSS subroutines and in an identical ", &
                  "manner. If you are uncertain about using it then please test your results!"
               end if

               flag_xc = 52

           else if (desc_str.eq."tpssloc" .or. desc_str.eq."TPSSLOC") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using TPSSloc gradient-corrected functionals."
               end if

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "You should be aware that this functional has not been comparatively tested ", &
                  "against equivalent implementations in other programs, as was performed for TPSS."
                  write(use_unit,'(2X,A,A)') &
                  "HOWEVER, it is implemented through the TPSS subroutines and in an identical ", &
                  "manner. If you are uncertain about using it then please test your results!"
               end if

               flag_xc = 53

           else if (desc_str.eq.'m08-hx' .or. desc_str.eq.'M08-HX') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using Hybrid M08-HX gradient-corrected functionals."
               end if

               flag_xc = 55

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.5223d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M08 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq.'m08-so' .or. desc_str.eq.'M08-SO') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using Hybrid M08-SO gradient-corrected functionals."
               end if

               flag_xc = 56

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.5679d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M08 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq.'m11' .or. desc_str.eq.'M11') then

               ! Hardcoding as this variable is fitted.
               hse_omega = 2.5d-1
               ! hse_omega_hf = hse_omega
               ! hse_omega_pbe = hse_omega
               hse_unit = 'b'
               flag_hse_unit = .true.

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,F15.8,A)') &
                     "XC: Using Hybrid M11 gradient-corrected functionals with OMEGA =",hse_omega," bohr^-1."
               end if

               flag_xc = 57

               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.428d0 !SR. LR = 1.0d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M11 hybrid functional has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if

           else if (desc_str.eq.'m11-l' .or. desc_str.eq.'M11-L') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC: Using M11-L gradient-corrected functionals."
               end if

               flag_xc = 58
! AJL, Jan 2018: Fixed by Bi Sheng so canonical SCAN now works.
!       else if (desc_str.eq.'scan' .or. desc_str.eq.'SCAN') then
!               call localorb_multi( &
!                   ' * At present, there are substantial concerns about the', &
!                   ' * accuracy of the potential of the canonical SCAN implementation.', &
!                   ' * Please, use the SCAN implementation provided by the dfauto', &
!                   ' * infrastructure instead. See the documentation for the xc tag', &
!                   ' * for details.' &
!                )
!               call aims_stop_coll()
            else if ((desc_str == 'scan') .or. (desc_str .eq. 'SCAN')) then

               ! Hardcoding as this variable is fitted.
               hse_omega = 0.0d-0
               ! hse_omega_hf = hse_omega
               ! hse_omega_pbe = hse_omega
               ! hse_unit = 'b'
               ! flag_hse_unit = .true.

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,F15.8,A)') &
                     "XC: Using the parameter-free meta-GGA functional SCAN."
               end if

               call cite_reference("SCAN")
               flag_xc = 59

            else if ((desc_str.eq.'lrc-xyg3').or.(desc_str.eq.'lrc-XYG3')) then
               read(inputline,*,end=189,err=99) desc_str,desc_str
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 30
               dftpt2_Ex_hf         = 8.033d-1
               dftpt2_Ec_osPT2      = 3.211d-1
               dftpt2_Ec_ssPT2      = 3.211d-1
               dftpt2_Ec_oslrcPT2   = dftpt2_Ex_hf - dftpt2_Ec_osPT2
               dftpt2_Ec_sslrcPT2   = dftpt2_Ex_hf - dftpt2_Ec_ssPT2
               use_mp2              = .true.
               use_dftpt2_and_lrc   = .true.
               lrc_pt2_unit         = 'B'
               lrc_pt2_omega        = 0.2
               !hse_omega            = 0.2
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional lrc-XYG3.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: Using lrc-PT2 with OMEGA =",lrc_pt2_omega," bohr^(-1)."
               endif
            else if ((desc_str.eq.'xyg3').or.(desc_str.eq.'XYG3')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional XYG3.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
               endif
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 30
               dftpt2_Ex_hf         = 8.033d-1
               dftpt2_Ec_osPT2      = 3.211d-1
               dftpt2_Ec_ssPT2      = 3.211d-1
               use_mp2              = .true.
            else if ((desc_str.eq.'xygjos').or.(desc_str.eq.'XYGJOS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional XYGJOS.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
               endif
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 32
               dftpt2_Ex_hf         = 7.731d-1
               dftpt2_Ec_osPT2      = 4.364d-1
               dftpt2_Ec_ssPT2      = 0.000d0
               use_os_mp2           = .true.
            else if ((desc_str.eq.'qpe-xygjos').or.(desc_str.eq.'QPE-XYGJOS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional XYGJOS.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
               endif
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 32
               dftpt2_Ex_hf         = 7.731d-1
               dftpt2_Ec_osPT2      = 4.364d-1
               dftpt2_Ec_ssPT2      = 0.000d0
               use_os_mp2_qpe       = .true.
            else if ((desc_str.eq.'xygjros').or.(desc_str.eq.'XYGJROS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional XYGJROS.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
               endif
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 32
               dftpt2_Ex_hf         = 7.731d-1
               dftpt2_Ec_osPT2      = 4.364d-1
               dftpt2_Ec_ssPT2      = 0.000d0
               use_os_rpa_qpe       = .true.
            else if ((desc_str.eq.'lrc-xyg6').or.(desc_str.eq.'lrc-XYG6')) then
               read(inputline,*,end=189,err=99) desc_str,desc_str
               flag_xc              = 10
               hybrid_coeff         = 0.2d0
               flag_dftpt2_dft_part = 34
               dftpt2_Ex_hf         = 8.780d-1
               dftpt2_Ec_osPT2      = 3.820d-1
               dftpt2_Ec_ssPT2      = 3.820d-1
               dftpt2_Ec_oslrcPT2   = 1.500d-1
               dftpt2_Ec_sslrcPT2   = 1.500d-1
               use_mp2              = .true.
               use_dftpt2_and_lrc   = .true.
               lrc_pt2_unit         = 'B'
               lrc_pt2_omega        = 0.2
               !hse_omega            = 0.2
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional lrc-XYG6.", &
                     " NOTE :: B3LYP REFERENCE IS REQUIRED"
                  write (use_unit,'(2X,A,E14.6,A)') &
                     "XC: Using lrc-PT2 with OMEGA =",lrc_pt2_omega," bohr^(-1)."
               endif
            else if ((desc_str.eq.'xdh-pbe0').or.(desc_str.eq.'xDH-PBE0')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using double-hybrid density functional xDH-PBE0", &
                     " NOTE :: PBE0 REFERENCE IS REQUIRED"
               endif
               flag_xc              = 1
               hybrid_coeff         = 2.5d-1
               flag_dftpt2_dft_part = 31
               dftpt2_Ex_hf         = 8.335d-1
               dftpt2_Ec_osPT2      = 5.428d-1
               dftpt2_Ec_ssPT2      = 0.0d0
               use_os_mp2           = .true.
            else if ((desc_str.eq.'zrps').or.(desc_str.eq.'ZRPS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                      "XC: Using one-parameter Level-5 ZRPS functional", &
                     " NOTE :: PBE0 REFERENCE IS REQUIRED"
               endif
               flag_xc              = 1
               hybrid_coeff         = 2.5d-1
               flag_dftpt2_dft_part = 33
               dftpt2_Ex_hf         = 5.000d-1
               dftpt2_Ec_osPT2      = 2.500d-1
               dftpt2_Ec_ssPT2      = 0.0d0
               flag_en_shift        = .true.
               en_shift_type        = 3
               coupling_pt2_factor  = 1.000d0
               coupling_pt2_screen  = 1.000d0
               coupling_pt2_shift   = 0.000d0
               write(info_str,'(2X,A,I6)') &
               'Energy shifting on the denominators: ', &
               en_shift_type
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Coupling factor for PT2 on the denominators: ', &
               coupling_pt2_factor
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Screening factor for coupling PT2: ', &
               coupling_pt2_screen
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Shift for screening effect in coupling PT2: ', &
               coupling_pt2_shift
               call localorb_info(info_str,use_unit,'(A)')
            else if (desc_str.eq."cohsex") then
               if ( myid.eq.0 ) then
                  write(use_unit,'(2X,2A)') &
                     "XC: Running screened exchange plus coulomb hole ", &
                     "(COHSEX) calculation ..."
               endif
               flag_xc = 0
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 1.d0
               end if
            else if (desc_str.eq."libxc") then
               ! flag_xc = 100

               ! First lets check if libxc is compiled in...
               call initialise_libxc()

               ! We need to work the ID of the functionals for initialisation
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str

               libxc_index_for_split = SCAN(desc_str,'+')
               if (libxc_index_for_split.gt.0) then
                 libxc_x_str = desc_str(1:libxc_index_for_split-1)
                 libxc_c_str = desc_str(libxc_index_for_split+1:)
               else
                 libxc_x_str = desc_str
                 libxc_c_str = ""
               endif

               flag_libxc_x_id = xc_f03_functional_get_number(libxc_x_str)
               flag_libxc_c_id = xc_f03_functional_get_number(libxc_c_str)

               ! Print out a general LibXC warning in case people get complacent
               call localorb_info('')
               write(info_str,'(2X,A)') "** PLEASE ** be careful using LibXC. Whilst most of the LDA, GGA and MGGA functionals, "
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "and hybrid variations therein, should be stable (though we always recommend testing!), "
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "there are several fancy functionals with e.g. range-separation, like HSE03/HSE06/LC-PBEh, "
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "and also MGGA functionals that depend on the laplacian of the density, like BR89, "
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "that are not fully implemented due to the infrastrucutre for these functionals being incomplete or non-generic. "
               call localorb_info(info_str)
               write(info_str,'(2X,A)') "If you want to use these functionals please contact the development team and we will work on them with you!"
               call localorb_info(info_str)
               call localorb_info('')
               if (.not.override_warning_libxc) &
                 call aims_stop_coll('Please disable this warning by adding "override_warning_libxc" in the line above your xc definition in control.in',func)

               ! Check if a hybrid has already been defined and throw error if so
               if (flag_hybrid_coeff) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid functionals from libxc have parameters optmized such that only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               endif

               ! Check if choice for libxc is valid
               if ( flag_libxc_x_id.gt.0 ) then
                 ! If so print out details
                 write(info_str,'(2x,a,a,a,i3.1,a)') 'You have selected libxc functional ' &
                               , trim(libxc_x_str),' (', flag_libxc_x_id ,')'
                 call localorb_info(info_str)
                 call check_libxc_func(flag_libxc_x_id,.false.,libxc_x_type,hybrid_coeff)
               else
                  ! Throw error message
                  write(use_unit,'(2X,4A)') &
                  "Your choice of libxc functional, ",  trim(libxc_x_str), ", doesn't seem to be recognised. ", &
                  "Please check the LIBXC manual and ensure you have made your definition correctly"
                  call aims_stop_coll('', func)
               endif

               ! Check if second part defined and if so use it as well
               if ( libxc_index_for_split.gt.0 ) then
                 ! Now check if this option is valid AND that it doesn't clash with first option
                 if ( flag_libxc_c_id.gt.0 ) then

                   write(info_str,'(2x,a,a,a,i3.1,a)') 'You have selected libxc functional ' &
                                 , trim(libxc_c_str),' (', flag_libxc_c_id ,')'
                   call localorb_info(info_str)
                   call check_libxc_func(flag_libxc_c_id,.false.,libxc_c_type,hybrid_coeff)

                   ! Check if libxc choices are both the same type - this is possible given our input construct
                   ! as e.g. one could input the correlation part first and exchange second. This
                   ! shouldn't cause an error, but if someone put e.g. 2 x exchange in then we should
                   ! catch it and tell them they have messed up.
                   !
                   ! The possible options are:
                   ! XC_EXCHANGE = 0,
                   ! XC_CORRELATION = 1,
                   ! XC_EXCHANGE_CORRELATION = 2,

                   if (libxc_x_type.eq.libxc_c_type .or. &
                      (libxc_x_type.eq.XC_EXCHANGE_CORRELATION) .or. &
                      (libxc_c_type.eq.XC_EXCHANGE_CORRELATION)) then
                     write(info_str,'(2x,A)') 'You seem to have selected incompatible LibXC options. You can either define an exchange '
                     call localorb_info(info_str)
                     write(info_str,'(2x,A)') 'functional, a correlation functional, both an exchange and a correlation functional OR an '
                     call localorb_info(info_str)
                     write(info_str,'(2x,A)') 'explicit exchange-correlation functional, but not other combinations! Please check your inputs, '
                     call localorb_info(info_str)
                     write(info_str,'(2x,4A)') trim(libxc_x_str), ' and ', trim(libxc_c_str), ' are compatible with these requirements '
                     call localorb_info(info_str)
                     call aims_stop_coll('',func)
                   endif
                else
                  ! Throw error message
                  write(use_unit,'(2X,4A)') &
                  "Your choice of libxc functional, ",  trim(libxc_c_str), ", doesn't seem to be recognised. ", &
                  "Please check the LIBXC manual and ensure you have made your definition correctly"
                  call aims_stop_coll('', func)
                endif
              endif

              ! Define flag_xc depending on the LibXC options

              ! To quote from atom_sphere_wrapper:

              ! "The first three numbers correspond to the exchange contribution to the
              ! functional.  The second three numbers correspond to the correlation
              ! contribution.  When the first three numbers are zero, that signifies
              ! that both exchange and correlation are combined, and the functional
              ! only needs to be initalized once.
              ! The negative sign is just there, because YOLO."

              ! See further down that everything is offset by another order of magnitude
              ! compared to atom_sphere_wrapper, but the general principle is the same. AJL/Dec 2016

              if (libxc_x_type.eq.XC_EXCHANGE_CORRELATION) then
                flag_xc = flag_libxc_x_id
              else if (libxc_x_type.eq.XC_EXCHANGE) then
                flag_xc = (1000*flag_libxc_x_id)+flag_libxc_c_id
              else ! Only remaining option is libxc_c actually contains the exchange option
                flag_xc = (1000*flag_libxc_c_id)+flag_libxc_x_id
              endif

              ! We need to shift this to make it identifiable. We are going to make the flag_xc negative for
              ! LibXC, and multiply by ten so it doesn't clash with flag_xc unset (= -1)
              flag_xc = flag_xc*(-10)

            else if (desc_str == "dfauto") then
                read (inputline, *, end=88, err=99) desc_str, desc_str, desc_str
                select case (str_lower(desc_str))
                case ('pw-lda'); flag_xc = xc__pw_lda
                case ('pbe'); flag_xc = xc__pbe
                case ('pbe0')
                    flag_xc = xc__pbe0
                    if (.not. flag_hybrid_coeff) hybrid_coeff = 0.25d0
                case ('tpss'); flag_xc = xc__tpss
                case ('scan'); flag_xc = xc__scan
                case ('scan0')
                    flag_xc = xc__scan0
                    if (.not. flag_hybrid_coeff) hybrid_coeff = 0.25d0
                case default
                    write (info_str, *) "Uknown dfauto functional:", desc_str
                    call aims_stop_coll(info_str)
                end select
                flag_xc = flag_xc + xc_dfauto_offset
                call cite_reference('dfauto')
            else
               if (myid.eq.0) then
                  write(use_unit,*) desc_str, " exchange-correlation is not ", &
                        "accessible at this time, feel free to change."
               end if
               call aims_stop_coll('', func)
            end if

         case('atomic_solver_xc')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq."kli" .or. desc_str.eq."KLI") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "Atomic Solver XC: Using Krieger-Li-Iafrate (KLI) ", &
                     "for generating the minimal NAO basis functions."
               endif
               flag_atomic_xc=12
            else
               if (myid.eq.0) then
                  write(use_unit,*) desc_str, &
                  " Exchange-correlation functional for the atomic solver", &
                  " is not allowed at this time, please change your", &
                  " 'atomic_solver_xc' setting."
               end if
               call aims_stop_coll('', func)
            endif

         case('lc_dielectric_constant')
            read(inputline,*,end=88,err=99) desc_str, lc_dielectric_constant
            if (lc_dielectric_constant.lt.1.0d0) then
                if (myid.eq.0) then
                    write(use_unit,'(1X,A)') &
                        "* Error: lc_dielectric_constant < 1"
                    write(use_unit,'(1X,A)') &
                        "* This does not make sense, as the values should always be larger or equal to 1 (vacuum)."
                end if
                write(info_str, "('dielectric_constant < 1 ')")
                call aims_stop_coll(info_str, func)
            end if
            if (myid.eq.0) then
               write(use_unit,'(2X,A,F18.8,A)') &
                  "lc_dielectric_constant: The dielectric constant for the LC-functionals is set to ", lc_dielectric_constant
               write(use_unit,'(2X,A,F18.8,A)') &
                  "   (See below for more information regarding the asymptotic decay of the exchange.)."
            end if

         case('hybrid_xc_coeff')
            read(inputline,*,end=88,err=99) desc_str, hybrid_coeff
            if (myid.eq.0) then
               write(use_unit,'(2X,A,F18.8,A)') &
                  "hybrid_xc_coeff: Mixing coefficient for hybrid-functional exact exchange modified to ", hybrid_coeff, " ."
               write(use_unit,'(2X,A)') &
                  "* WARNING: This change modifies the chosen exchange-correlation functional (keyword 'xc')"
               write(use_unit,'(2X,A)') &
                  "* from its literature definition. Be sure to NEVER include such a line in any default control.in files!"
            end if

            flag_hybrid_coeff = .true.
         case('hse_unit')
            read(inputline,*,end=88,err=99) desc_str, hse_unit

            if      ( (hse_unit.eq.'a') .or. (hse_unit.eq.'A') ) then
              if (myid.eq.0) then
                write(use_unit,'(2X,A)') &
                  "hse_unit: Unit for the HSE06 hybrid functional screening parameter set to A^(-1)."
              end if
              hse_unit = 'A'
            else if ( (hse_unit.eq.'b') .or. (hse_unit.eq.'B') ) then
              if (myid.eq.0) then
                write(use_unit,'(2X,A)') &
                  "hse_unit: Unit for the HSE06 hybrid functional screening parameter set to bohr^(-1)."
              end if
              hse_unit = 'B'
            else
              if (myid.eq.0) then
                write(use_unit,'(1X,A,A)') &
                  "* hse_unit: Illegal setting ", hse_unit
                write(use_unit,'(1X,A)') &
                  "* Valid settings are either 'A' (for Angstrom^-1) or 'B' (for bohr^-1)."
                write(use_unit,'(1X,A)') &
                  "* Please specify a legitimate value. Apologies for the inconvenience."
              end if
               write(info_str, "('Incorrect value for hse_unit: ',A)") trim(hse_unit)
               call aims_stop_coll(info_str, func)
            end if

            flag_hse_unit = .true.
         case('split_atoms')
             write(info_str,'(A)') 'Atoms might be split for the evaluation of the exact exchange part'
             call localorb_info(info_str)
         case('split_min_val')
            read(inputline,*,end=88,err=99) desc_str, split_min_val
         case('split_max_val')
            read(inputline,*,end=88,err=99) desc_str, split_max_val
         case('split_batch')
            read(inputline,*,end=88,err=99) desc_str, split_batch
         case('use_hf_kspace')
            read(inputline,*,end=88,err=99) desc_str, use_hf_kspace
            if(use_hf_kspace) then
               write(info_str, "(A)") "The k-space periodic HF implementation will be used."
               call localorb_info(info_str)
            else
               write(info_str, "(A)") "The real-space periodic HF implementation will be used."
               call localorb_info(info_str)
            endif
         case('KS_eigenfunc_conjg')
               write(info_str, "(A)") "Kohn-Sham eigenfunctions are complex conjugated in periodic RPA calculations."
               call localorb_info(info_str)

         case('use_hf_kspace_with_rpa')
            read(inputline,*,end=88,err=99) desc_str, use_hf_kspace_with_rpa
            if(use_hf_kspace_with_rpa) then
               write(info_str, "(A)") "The k-space periodic HF implementation will be used in EX+RPA."
               call localorb_info(info_str)
            else
               write(info_str, "(A)") "The real-space periodic HF implementation will be used in EX+RPA."
               call localorb_info(info_str)
            endif

         case('periodic_hf')
            ! keywords for periodic Hartree-Fock, HSE et al.

            read(inputline,*,end=88,err=99) desc_str, desc_str

            if ( (desc_str.eq.'screening_threshold') .or. (desc_str.eq.'crit_val') )then
              ! Approx error in exchange matrix entries:

              read(inputline,*,end=88,err=99) desc_str, desc_str, crit_val

              if ( (crit_val.gt.1d-6) .or. (crit_val.lt.0.d0) ) then
                 write(info_str, "(1X,A)") &
                   "* Periodic Hartree-Fock: Requested screening threshold for four-center integrals (crit_val)"
                 call localorb_info(info_str)
                 write(info_str, "(1X,A,E14.6)") &
                   "* is either too large, or below zero: ", crit_val
                 call localorb_info(info_str)
                 write(info_str, "(1X,A)") &
                   "* Please correct this setting before continuing."
                 call localorb_info(info_str)
                 write(info_str, "('Illegal value for crit_val: ',E14.6)") crit_val
                 call aims_stop_coll(info_str, func)
              else
                 write (info_str, '(2X,A,E14.6)') &
                   "Periodic Hartree-Fock: Screening threshold for four-center integrals: ", crit_val
                 call localorb_info(info_str)
                 flag_crit_val = .true.
              endif

            else if ( (desc_str.eq.'coulomb_threshold') .or. (desc_str.eq.'coul_mat_threshold') ) then
               ! Rows/columns in the coulomb matrices which are completely below coul_mat_threshold are left away!

              read(inputline,*,end=88,err=99) desc_str, desc_str, coul_mat_threshold

              if ( (coul_mat_threshold.gt.1d-6) .or. (coul_mat_threshold.lt.0.d0) ) then
                 write(info_str, "(1X,A)") &
                   "* Periodic Hartree-Fock: Requested screening threshold for Coulomb matrix (coul_mat_threshold)"
                 call localorb_info(info_str)
                 write(info_str, "(1X,A,E14.6)") &
                   "* is either too large, or below zero: ", coul_mat_threshold
                 call localorb_info(info_str)
                 write(info_str, "(1X,A)") &
                   "* Please correct this setting before continuing."
                 call localorb_info(info_str)
                 write(info_str, "('Illegal value for coul_mat_threshold: ',E14.6)") coul_mat_threshold
                 call aims_stop_coll(info_str, func)
              else
                 write (info_str, '(2X,A,E14.6)') &
                   "Periodic Hartree-Fock: Screening threshold for Coulomb matrix elements: ", coul_mat_threshold
                 call localorb_info(info_str)
                 flag_coul_mat_threshold = .true.
              endif

            else if ( (desc_str.eq.'n_cells_bunch') ) then
               ! number of cells per bunch: no longer active, just kept in order to not stop anyone's old control.in file needlessly
               write (info_str, '(1X,A)') &
                 "* Periodic Hartree-Fock: Keyword n_cells_bunch is no longer active and will be ignored."
               call localorb_info(info_str)

            else

               write(info_str, "('Periodic Hartree-Fock - unknown keyword: ',A)") trim(desc_str)
               call aims_stop_coll(info_str, func)

            end if

         case('RI_method')
            !
            ! If you add a new RI_method, please do not forget to
            ! modify the corresponding block in dimensions.f90 .
            !
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'V'.or.desc_str.eq.'v') then
               RI_type = RI_V
               write(info_str, "(2A)") " The 'V' version of RI (resolution of identity)", &
                  & " technique is used."
               call localorb_info(info_str)
            else if (desc_str.eq.'LVL'.or.desc_str.eq.'lvl'.or.desc_str.eq.'LVL_fast'.or.desc_str.eq.'lvl_fast') then
               RI_type = RI_LVL
               write(info_str, "(2X,A,' ',A)") "The 'LVL' version of RI", &
                  & "(resolution of identity) technique is used."
               call localorb_info(info_str)
               call localorb_info("  ** This is a localized expansion of all basis function products.")
               call localorb_info("     Its intended use are Hartree-Fock and hybrid functionals.")
               call localorb_info("     For RPA, MP2, GW, etc., LVL is not yet as accurate. In those cases, check carefully.")

            else if (desc_str.eq.'SVS'.or.desc_str.eq.'svs') then
               RI_type = RI_SVS
               write(info_str, "(2A)") " The 'SVS' version of RI (resolution of identity)", &
                  & " technique is used."
               call localorb_info(info_str)
            else if (desc_str.eq.'LVL_full'.or.desc_str.eq.'lvl_full') then
               RI_type = RI_LVL_full
               write(info_str, "(2X,A)") &
               & "Using RI-LVL with dense matrices (LVL_full)."
               call localorb_info(info_str)
               write(info_str, "(2A)") &
               & "  ** Please note that this gathers the drawbacks", &
               & " of both RI-V and RI-LVL"
               call localorb_info(info_str)
            else if (desc_str.eq.'LVL_2nd'.or.desc_str.eq.'lvl_2nd') then
               RI_type = RI_LVL_2nd
               write(info_str, "(2X,A)") &
               & "Using RI-LVL with first order correction."
               call localorb_info(info_str)
               write(info_str, "(2X,2A)") &
               & "** Please note that its implementation is ", &
               & "for accuracy testing, only."
               call localorb_info(info_str)
            else
               write(info_str, "('Unknown RI version: ',A)") trim(desc_str)
               call aims_stop_coll(info_str, func)
            endif
            sparse_o3fn = (RI_type == RI_LVL)
            ! use_lvl is now read and set in dimensions.f90 and should not be changed later.
            ! use_lvl = (RI_type == RI_LVL .or. RI_type == RI_LVL_full .or. &
            ! &                                 RI_type == RI_LVL_2nd)

            flag_RI = .true. ! RI_method was explicitly requested in control.in




         case('qpe_calc')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if ((desc_str.eq."gw").or.(desc_str.eq.'GW')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "GW quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if (desc_str.eq."ev_scgw") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Eigenvalue self-consistent GW quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if (desc_str.eq."ev_scgw0") then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Eigenvalue self-consistent GW0 quasiparticle calculation (without MO update) ", &
                     "of excited state will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."gw+2ox").or.(desc_str.eq.'GW+2OX')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "GW+2OX quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."gw+soxw").or.(desc_str.eq.'GW+SOXW')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "GW plus second-order exchange in W quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."gw+s2ox").or.(desc_str.eq.'GW+S2OX')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "GW+S2OX (GW plus screened 2nd-order exchange) quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."gw+sosex2w").or.(desc_str.eq.'GW+SOSEX2W')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,3A)') &
                  "GW+SOSEX2W (GW plus screened 2nd-order exchange with two screend Coulomb lines) ", &
                    " quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."mp2").or.(desc_str.eq.'MP2')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "MP2 quasiparticle calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq."gw_expt").or.(desc_str.eq.'GW_EXPT')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Periodic GW calculations will be started if finite k point sampling is also set. ", &
                     "Please be warned that periodic GW is only at the experimental stage and is still under testing. Be cautious about your results."
               endif
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "no alternative method other than GW, GW+2OX, GW+2SOX, GW+SOSEX2W and MP2 for excited ", &
                     "states calculation available, change back to GW/MP2."
               endif
               call aims_stop_coll('', func)
            end if

         case('path_integral')        ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if ((desc_str.eq."pimd").or.(desc_str.eq.'PIMD')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "PIMD simulation for the statistical quantum nuclear ", &
                     "properties will be started."
               endif
            else if ((desc_str.eq."cmd").or.(desc_str.eq.'CMD')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "CMD simulation for the statistical and dynamical ", &
                     "quantum nuclear properties will be started."
               endif
            else if ((desc_str.eq."rpmd").or.(desc_str.eq.'RPMD')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "RPMD simulation for the statistical and dynamical ", &
                     "quantum nuclear properties will be started."
               endif
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "no alternative method other than PIMD/CMD/RPMD for quantum ", &
                     "nuclear effect calculation available."
               endif
               call aims_stop_coll('', func)
            end if
            if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Features connected to PIMD are highly experimental. Please do not use them unless you really", &
                     "know what you are doing. For using them, comment the stop flag in read_control.f90."
            endif
            call aims_stop_coll('', func)

         case('use_split_xc_gw')
            read(inputline,*,end=88,err=99) desc_str, use_split_xc_gw

         case('sc_self_energy')
            read(inputline,*,end=88,err=99) desc_str, desc_str
           if ((desc_str.eq."scgw").or.(desc_str.eq.'SCGW')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Self Consistent GW calculation of excited states ", &
                     "will be started after the DFT/HF calculation."
               endif
          elseif ((desc_str.eq."scgw0").or.(desc_str.eq.'SCGW0'))then
             if (myid.eq.0) then
                write(use_unit,'(2X,A,A)') &
                 "Self Consistent GW calculation of excited states ", &
                  "will be started after the DFT/HF calculation."
             endif

            elseif ((desc_str.eq."dmft_pbe0").or.(desc_str.eq.'DMFT_PBE0')&
                    .or.(desc_str.eq.'dmft_PBE0').or.(desc_str.eq.'DMFT_pbe0')) then
               read(inputline,*,end=88,err=99)desc_str,desc_str,dmft_alpha, hybrid_alpha

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DMFT embeding of PBE0 self-energy calculation ", &
                     "will be started after the DFT calculation."
               endif
            elseif ((desc_str.eq."dmft_gw").or.(desc_str.eq.'DMFT_GW')&
                    .or.(desc_str.eq.'dmft_GW').or.(desc_str.eq.'DMFT_gw')) then
               read(inputline,*,end=88,err=99)desc_str,desc_str,dmft_alpha

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DMFT embeding of GW self-energy calculation ", &
                     "will be started after the DFT calculation."
               endif

           endif
         case('total_energy_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if ((desc_str.eq.'HF').or.(desc_str.eq.'hf')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "HF exchange energy calculation ", &
                     "will be started after the DFT calculation."
               endif
               ! I don't really like the fact we have different behaviour for HF calculations
               ! and hybrid-DFT here. For the HF calculation, post-processing is just turned off.
               ! For hybrid-DFT, we throw an error. It'd be better if both had the same behaviour.
               if (use_hartree_fock .and. flag_xc.ne.0) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'RPA').or.(desc_str.eq.'rpa').or.(desc_str.eq.'klein')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "RPA correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'osRPA').or.(desc_str.eq.'osrpa')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "os-RPA correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'osPT2').or.(desc_str.eq.'ospt2')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Quasiparticle os-PT2 correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'SIC-RPA').or.(desc_str.eq.'sic-rpa')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Self-interaction-corrected RPA correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'GW').or.(desc_str.eq.'gw')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "GW Total energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'RPA+2OX').or.(desc_str.eq.'rpa+2ox')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "RPA+2OX correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'RPA+SOSEX').or.(desc_str.eq.'rpa+sosex') &
                     .or.(desc_str.eq.'SOSEX').or.(desc_str.eq.'sosex')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "RPA+SOSEX correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'rPT2').or.(desc_str.eq.'rpt2') ) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "rPT2 (=RPA+SOSEX+rSE) correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'mp2').or.(desc_str.eq.'MP2')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "MP2 correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'mp2_blacs').or.(desc_str.eq.'MP2_BLACS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "MP2 correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'os-pt2').or.(desc_str.eq.'OS-PT2')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Only the opposite-spin component of PT2 correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'ci').or.(desc_str.eq.'CI')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Configuration interaction correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'cc').or.(desc_str.eq.'CC')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Coupled cluster correlation energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'lrc-xyg3').or.(desc_str.eq.'lrc-XYG3')) then
               if (myid .eq. 0) then
                  write(use_unit,'(2X,A)') &
                      "Running post-KS double-hybrid lrc-XYG3 calculation"
               endif
               if (flag_xc .ne. 10) then
                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,A)') &
                         "Error: Double-hybrid density functional lrc-XYG3 ", &
                         "must be used together with B3LYP reference."
                      write(use_unit,'(2X,A)') &
                         "Usage: xc                   B3LYP"
                      write(use_unit,'(2X,A)') &
                         "       total_energy_method  lrc-XYG3"
                      write(use_unit,'(2X,A)') &
                         "   or: xc                   lrc-XYG3"
                   endif
                   call aims_stop_coll('', func)
               endif
               flag_dftpt2_dft_part = 30
               dftpt2_Ex_hf         = 8.033d-1
               dftpt2_Ec_osPT2      = 3.211d-1
               dftpt2_Ec_ssPT2      = 3.211d-1
               dftpt2_Ec_oslrcPT2   = dftpt2_Ex_hf - dftpt2_Ec_osPT2
               dftpt2_Ec_sslrcPT2   = dftpt2_Ex_hf - dftpt2_Ec_ssPT2
               use_mp2              = .true.
               use_dftpt2_and_lrc   = .true.
               lrc_pt2_unit         = 'B'
               lrc_pt2_omega        = 0.2
               flag_hse_unit        = .true.
            else if ((desc_str.eq.'xyg3').or.(desc_str.eq.'XYG3')) then
               if (myid .eq. 0) then
                  write(use_unit,'(2X,A)') &
                      "Running post-KS double-hybrid XYG3 calculation"
               endif
               if (flag_xc .ne. 10) then
                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,A)') &
                         "Error: Double-hybrid density functional XYG3 ", &
                         "must be used together with B3LYP reference."
                      write(use_unit,'(2X,A)') &
                         "Usage: xc                   B3LYP"
                      write(use_unit,'(2X,A)') &
                         "       total_energy_method  XYG3"
                      write(use_unit,'(2X,A)') &
                         "   or: xc                   XYG3"
                   endif
                   call aims_stop_coll('', func)
               endif
               flag_dftpt2_dft_part = 30
               dftpt2_Ex_hf         = 8.033d-1
               dftpt2_Ec_osPT2      = 3.211d-1
               dftpt2_Ec_ssPT2      = 3.211d-1
               use_mp2              = .true.
            else if ((desc_str.eq.'xygjos').or.(desc_str.eq.'XYGJOS')) then
               if (myid .eq. 0) then
                  write(use_unit,'(2X,A)') &
                      "Running post-KS double-hybrid XYGJOS calculation"
               endif
               if (flag_xc .ne. 10) then
                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,A)') &
                         "Error: Double-hybrid density functional XYGJOS ", &
                         "must be used together with B3LYP reference."
                      write(use_unit,'(2X,A)') &
                         "Usage: xc                   B3LYP"
                      write(use_unit,'(2X,A)') &
                         "       total_energy_method  XYGJOS"
                      write(use_unit,'(2X,A)') &
                         "   or: xc                   XYGJOS"
                   endif
                   call aims_stop_coll('', func)
               endif
               flag_dftpt2_dft_part = 32
               dftpt2_Ex_hf         = 7.731d-1
               dftpt2_Ec_osPT2      = 4.364d-1
               dftpt2_Ec_ssPT2      = 0.000d0
               use_os_mp2           = .true.
            else if ((desc_str.eq.'xdh-pbe0').or.(desc_str.eq.'xDH-PBE0')) then
               if (myid .eq. 0) then
                  write(use_unit,'(2X,A)') &
                      "Running post-KS double-hybrid xDH-PBE0 calculation"
               endif
               if (flag_xc .ne. 1 .and. flag_xc .ne. 7) then
                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,A)') &
                         "Error: Double-hybrid density functional xDH-PBE0 ", &
                         "must be used together with PBE0 or HSE06 reference."
                      write(use_unit,'(2X,A)') &
                         "Usage: xc                   PBE0 .or. HSE06"
                      write(use_unit,'(2X,A)') &
                         "       total_energy_method  xDH-PBE0"
                      write(use_unit,'(2X,A)') &
                         "   or: xc                   xDH-PBE0"
                   endif
                   call aims_stop_coll('', func)
               endif
               flag_dftpt2_dft_part = 31
               dftpt2_Ex_hf         = 8.335d-1
               dftpt2_Ec_osPT2      = 5.428d-1
               dftpt2_Ec_ssPT2      = 0.0d0
               use_os_mp2           = .true.
            else if ((desc_str.eq.'zrps').or.(desc_str.eq.'ZRPS')) then
               if (myid .eq. 0) then
                  write(use_unit,'(2X,A)') &
                      "Running post-KS double-hybrid ZRPS calculation"
               endif
               if (flag_xc .ne. 1 .and. flag_xc .ne. 7) then
                   if (myid.eq.0) then
                      write(use_unit,'(2X,A,A)') &
                         "Error: One-parameter level-5 functional ZRPS ", &
                         "must be used together with PBE0 or HSE06 reference."
                      write(use_unit,'(2X,A)') &
                         "Usage: xc                   PBE0 .or. HSE06"
                      write(use_unit,'(2X,A)') &
                         "       total_energy_method  ZRPS"
                      write(use_unit,'(2X,A)') &
                         "   or: xc                   ZRPS"
                   endif
                   call aims_stop_coll('', func)
               endif
               use_os_mp2           = .true.
               flag_dftpt2_dft_part = 33
               dftpt2_Ex_hf         = 5.000d-1
               dftpt2_Ec_osPT2      = 2.500d-1
               dftpt2_Ec_ssPT2      = 0.0d0
               flag_en_shift        = .true.
               en_shift_type        = 3
               coupling_pt2_factor  = 1.000d0
               coupling_pt2_screen  = 1.000d0
               coupling_pt2_shift   = 0.000d0
               write(info_str,'(2X,A,I6)') &
               'Energy shifting on the denominators: ', &
               en_shift_type
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Coupling factor for PT2 on the denominators: ', &
               coupling_pt2_factor
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Screening factor for coupling PT2: ', &
               coupling_pt2_screen
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,F16.8)') &
               'Shift for screening effect in coupling PT2: ', &
               coupling_pt2_shift
               call localorb_info(info_str,use_unit,'(A)')
            else if (desc_str.eq.'C6_coef') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "C6 dispersion coefficients at MP2/RPA level ", &
                     "will be calculated after the DFT/HF calculation."
               endif
            else if ((desc_str.eq.'ll_vdwdf').or.(desc_str.eq.'LL_vdwdf')) then
               if ((flag_xc.eq.6).or.(flag_xc.eq.12)) then
                  if (.not.out_vdwdf) then
                     if (myid.eq.0) then
                        write(use_unit,'(2X,A,A)') &
                        "Langreth-Lundqvist van der Waals density functional ",&
                        "calculation will be started after the DFT/HF calculation."
                     endif
                     call read_ll_vdw_parameters
                     out_vdwdf=.true.
                  else
                     if (myid.eq.0) then
                        write(info_str,'(2X,A)') &
                              'll_vdwdf parameters have been already provided: '
                     endif
                     call aims_stop_coll(info_str, func)
                  endif
               else

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "* LL_vdwdf not possible with current xc-functional ", &
                           "-- this is desinged for revpbe."
                     write (use_unit,'(2X,A,A)') &
                           "* Please correct it and rerun."
                  endif
                  call aims_stop_coll('', func)
               endif
            else if ((desc_str.eq.'m06-l').or.(desc_str.eq.'M06-L')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga M06-L energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 1
            else if ((desc_str.eq.'m06-hf').or.(desc_str.eq.'M06-HF')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid meta-gga M06-HF energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 9
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 1.0d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if
               if (use_hartree_fock) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'m06').or.(desc_str.eq.'M06')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid meta-gga M06 energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 2
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.27d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if
               if (use_hartree_fock) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'m06-2x').or.(desc_str.eq.'M06-2X')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid meta-gga M06-2X energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 5
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.54d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M06 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if
               if (use_hartree_fock) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'m08-hx').or.(desc_str.eq.'M08-HX')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid meta-gga M08-HX energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 10
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.5223d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M08 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if
               if (use_hartree_fock) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'m08-so').or.(desc_str.eq.'M08-SO')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Hybrid meta-gga M08-SO energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 11
               if (.not.flag_hybrid_coeff) then
                  hybrid_coeff = 0.5679d0
               else
                  write(use_unit,'(2X,A,A)') &
                  "The M08 family of hybrid functionals has parameters optmized such that it only ", &
                  "the specified hybrid mixing works with them. Don't use any other hybrid coeff."
                  call aims_stop_coll('', func)
               end if
               if (use_hartree_fock) then
                  write(use_unit,'(2X,A,A)') &
                  "The usage of a hybrid xc functional with ", &
                  "hybrid post-processing is not yet supported, sorry."
                  call aims_stop_coll('', func)
               end if
            else if ((desc_str.eq.'m11').or.(desc_str.eq.'M11')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga M11 energy calculation ", &
                     "will be started after the DFT/HF calculation."
                  !! I haven't managed to get this working, to have to throw an error for now. To-do. AJL
                  write(use_unit,'(2X,A,A)') &
                  "The long-range corrected M11 functional has not been implemented yet for post-processing.", &
                  " Please use an alternative for now."
                  call aims_stop_coll('', func)
               endif
               flag_post_xc = 12
            else if ((desc_str.eq.'m11-l').or.(desc_str.eq.'M11-L')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga M11-L energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 13
            elseif ( (desc_str.eq.'NLCORR').or.(desc_str.eq.'nlcorr')  )then !SAG
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                       "Nonlocal correlation energy calculation ", &
                       "will be started after the DFT/HF calculation."
               endif
            elseif ( (desc_str.eq.'pbe_vdw').or.(desc_str.eq.'PBE_VDW')  )then !SAG
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                       "PBE with VDW calculation ", &
                       "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 3
            elseif ( (desc_str.eq.'revpbe_vdw').or.(desc_str.eq.'REVPBE_VDW')  )then !SAG
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                       "REVPBE with VDW calculation ", &
                       "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 4
!@@edu>> tpss
            else if ((desc_str.eq.'tpss').or.(desc_str.eq.'TPSS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga TPSS energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 6

            else if ((desc_str.eq.'revtpss').or.(desc_str.eq.'revTPSS')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga revTPSS energy calculation ", &
                     "will be started after the DFT/HF calculation."
               endif
               flag_post_xc = 7

            else if ((desc_str.eq.'tpssloc').or.(desc_str.eq.'TPSSloc')) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Meta-gga TPSSloc energy calculation ", &
                     "will be started after the DFT/HF calculation."
                  write(use_unit,'(2X,A,A)')&
         &           "Reference: L.A. Constantin, E. Fabiano, F.Della Sala,",&
         &           " Phys. Rev. B. 86, 035130 (2012)"
               endif
               flag_post_xc = 8
!@@edu<< tpss


            else if ((desc_str == 'scan') .or. (desc_str .eq. 'SCAN')) then
                if (myid == 0) then
                    write (use_unit, '(2X,A,A)') &
                        "Meta-gga SCAN energy calculation ", &
                        "will be started after the DFT/HF calculation."
                endif
                call cite_reference("SCAN")
                flag_post_xc = 14

! AJL, 5/3/2019
! Not currently enabled in evaluate_post_xc.f90 and so diabling here to prevent
! spurious use. It would be nice long-term to remove this file completely and just use the
! evaluate_xc infrastructure, which does not differ significantly.
!
!            else if ((desc_str.eq.'b3lyp').or.(desc_str.eq.'B3LYP')) then
!               if (myid.eq.0) then
!                    write (use_unit, '(2X,A,A)') &
!                        "Hybrid B3LYP energy calculation ", &
!                        "will be started after the DFT/HF calculation."
!                  write(use_unit,'(2X,A)') &
!                        "| This version of B3LYP uses the RPA parametrization of LDA given by"
!                  write(use_unit,'(2X,A)') &
!                        "| Vosko, Wilk, Nusair. This choice differs from Becke's original paper,"
!                  write(use_unit,'(2X,A)') &
!                        "| and goes back to a now de facto standard (originally, erroneous)"
!                  write(use_unit,'(2X,A)') &
!                        "| implementation in the Gaussian suite of codes."
!               endif
!               flag_post_xc = 15
!               if (.not.flag_hybrid_coeff) then
!                  hybrid_coeff = 0.2d0
!               else
!                  write(use_unit,'(2X,A,A)') &
!                  "The B3LYP hybrid functional is parameterized such that only ", &
!                  "the specified hybrid mixing works. Don't use any other hybrid coeff."
!                  call aims_stop_coll('', func)
!               end if
!               if (use_hartree_fock) then
!                  write(use_unit,'(2X,A,A)') &
!                  "The usage of a hybrid xc functional with ", &
!                  "hybrid post-processing is not yet supported, sorry."
!                  call aims_stop_coll('', func)
!               end if

            else if (desc_str == "dfauto") then
                read (inputline, *, end=88, err=99) desc_str, desc_str, desc_str
                select case (str_lower(desc_str))
                case ('pw-lda'); flag_post_xc = xc__pw_lda
                case ('pbe'); flag_post_xc = xc__pbe
                case ('pbe0')
                    flag_post_xc = xc__pbe0
                    if (.not. flag_hybrid_coeff) hybrid_coeff = 0.25d0
                case ('tpss'); flag_post_xc = xc__tpss
                case ('scan'); flag_post_xc = xc__scan
                case ('scan0')
                    flag_post_xc = xc__scan0
                    if (.not. flag_hybrid_coeff) hybrid_coeff = 0.25d0
                case default
                    write (info_str, *) "Unknown dfauto functional:", desc_str
                    call aims_stop_coll(info_str)
                end select
                flag_post_xc = flag_post_xc + xc_dfauto_offset
                call cite_reference('dfauto')
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "This specified method is not available for accurate ", &
                     "total energy calculation, stop!"
               endif
               call aims_stop_coll('', func)
            endif

         case('k_grid')
            read(inputline,*,end=88,err=99) desc_str, n_k_points_xyz(1), n_k_points_xyz(2), &
                  n_k_points_xyz(3)

            n_kpts = n_k_points_xyz(1) *  n_k_points_xyz(2) * n_k_points_xyz(3)

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I10,I10,I10)') &
                  'Found k-point grid:', &
                  n_k_points_xyz(1), n_k_points_xyz(2), &
                  n_k_points_xyz(3)
            end if

            if (n_kpts.le.0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') "* Error: Less than 1 k point requested."
                  write(use_unit,'(1X,A)') "* This is unphysical, please correct."
               end if
               call aims_stop_coll('', func)

            end if

            flag_k_points_defined = .true.


            ! SVL/AS: Perturbative DOS
         case('dos_kgrid_factors')
            read(inputline,*,end=88,err=99) desc_str, n_k_points_dos_xyz(1), n_k_points_dos_xyz(2), &
                  n_k_points_dos_xyz(3)

            if (myid.eq.0) then
               write(use_unit,'(2X,A,I10,I10,I10)') &
                  'Found k-point grid multipliers for perturbative DOS :', &
                  n_k_points_dos_xyz(1), n_k_points_dos_xyz(2), &
                  n_k_points_dos_xyz(3)
            end if
            if (n_k_points_dos_xyz(1)*n_k_points_dos_xyz(2)*n_k_points_dos_xyz(3) .ge. 1) then
               pert_dos_on = .true.
            else
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                     '* Erroneous dos_kgrid_factors? Their product is less than one.'
                  write(use_unit,'(1X,A)') &
                     '* Stopping due to suspected typo or misunderstanding.'
                  call aims_stop_coll('Invalid dos_kgrid_factors', func)
               end if
               pert_dos_on = .false.
            end if
         case('k_offset')
            read(inputline,*,end=88,err=99) desc_str, k_points_offset(1), k_points_offset(2), &
                  k_points_offset(3)

            if (myid.eq.0) then
               write(use_unit,'(2X,A,F10.5,F10.5,F10.5)') &
                  'K-point grid offset (additive shift of all k-points):', &
                  k_points_offset(1), k_points_offset(2), &
                  k_points_offset(3)
            end if

            if (use_hartree_fock) then
               ! An offset for the k-point grid is not safe in this case.
               ! Note that I believe that an offset could well be made to work,
               ! but anyone trying this should please look very carefully at the
               ! way in which k-points and k-point grid symmetry assumptions enter the
               ! exact exchange code.
               ! Also, if anyone has experiences to the contrary, please let me know. - VB
               if ( (k_points_offset(1).ne.0.d0) .or. (k_points_offset(2).ne.0.d0) &
                    .or. (k_points_offset(3).ne.0.d0) ) &
                  then
                  if (myid.eq.0) then
                      write(use_unit,'(1X,A,A)') &
                        '* Error. Non-zero offsets to the k-point grid are not supported ',&
                        'for functionals including exact exchange at present.'
                      write(use_unit,'(1X,A)') &
                        '* Stopping to avoid potentially wrong results.'
                  end if
                  call aims_stop_coll('Not safe to use nonzero k_offset with exact exchange.', func)
               end if
            end if

         case('symmetry_reduced_k_grid')
            read(inputline,*,end=88,err=99) &
            & desc_str,  use_symmetry_reduced_k_grid
            if (use_symmetry_reduced_k_grid) then
               write(info_str,'(2X,A)') &
               & "Explicitly requesting symmetry reduced k-points."
               call localorb_info(info_str)
               write(info_str,'(2X,A)') &
               & "Currently, only time-reversal symmetry is supported."
               call localorb_info(info_str)
            else
               write(info_str,'(2X,A)') &
                  "Request to not use symmetry reduced k-point grids."
               call localorb_info(info_str)
            end if
            if (n_periodic == 0) then
               write(info_str,'(2X,A)') &
               & "* Not a periodic system - this option should have no effect."
               call localorb_info(info_str)
            end if
            flag_symmetry_reduced_k_grid = .true.

         case('symmetry_thresh')
            read(inputline,*,end=88,err=99) desc_str, symmetry_thresh
            write(info_str,'(2X,A,E14.6)') "Found a symmetry threshold of ", symmetry_thresh
            call localorb_info(info_str)

         case('use_full_symmetry')
            read(inputline,*,end=88,err=99) desc_str, use_full_symmetry
            write(info_str,'(2X,A)') "Using full point group symmetry in periodic RPA/GW calculations."
            call localorb_info(info_str)

         case('k_points_external')

            if (use_hartree_fock) then
               ! An external k-point list is presently not supported for any functionals
               ! using exact exchange. Even for a Gamma-centered grid, the code currently
               ! does not fill some necessary variables with meaningful values.
               ! Presently, only Gamma-centered grids using the k_grid tag are supported for periodic exact exchange.
               if (myid.eq.0) then
                   write(use_unit,'(1X,A,A)') &
                     "* Error - 'k_points_external' is not supported ",&
                     "for functionals including exact exchange at present."
               end if
               call aims_stop_coll("Cannot use 'k_points_external' for functionals that include nonlocal exchange.", func)
            end if

            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
               "k_points_external: Expecting k-points in file 'k_list.in'."
            end if

            ! check if k_list.in exists
            open(88,file='k_list.in')

            read(88,*,iostat=i_code) desc_str
            if (i_code.ne.0) then
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)') "* Empty input file 'k_list.in'. "
                  write (use_unit,'(1X,A)') "* Needed for keyword 'k_points_external'."
               end if
               call aims_stop_coll('', func)
            end if

            read(88,*) n_kpts
            if (myid.eq.0) then
               write (use_unit,'(2X,A,I8,A)') &
                  "Header of 'k_list.in' anticipates ", n_kpts, " k-points."
            end if

            close(88)

            if (n_kpts.le.0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') "* Error: Less than 1 k-point requested."
                  write(use_unit,'(1X,A)') "* This is unphysical, please correct."
               end if
               call aims_stop_coll('', func)

            end if

            flag_k_points_defined = .true.
            read_k_points_from_file = .true.
            out_k_point_list = .true.

         case('legacy_monopole_extrapolation')
            read(inputline,*,end=88,err=99) desc_str, &
            & legacy_monopole_extrapolation
            if (legacy_monopole_extrapolation) then
               call localorb_info('  Legacy monopole extrapolation explicitly required.')
            else
               call localorb_info('  Using new monopole extrapolation.')
            end if
            flag_legacy_monopole_extrapolation = .true.

         case('packed_matrix_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
                  packed_matrix_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of packed matrixes:', &
                     packed_matrix_threshold
            end if

         case('multip_radius_free_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
                  multipole_radius_free_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of hartree multipole radius', &
            'of free atoms:', &
                     multipole_radius_free_threshold
            end if

         case('multip_radius_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
                  far_distance_hartree_multipole_radius_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of hartree multipole radius:', &
                  far_distance_hartree_multipole_radius_threshold
            end if

         case('multip_moments_rad_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
                  far_distance_hartree_multipole_moment_radius_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of hartree multipole moments radius:', &
                  far_distance_hartree_multipole_moment_radius_threshold
            end if

         case('multip_moments_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
               far_distance_hartree_multipole_moment_threshold
            if (myid.eq.0) then
              write(use_unit,*) &
                 '| Threshold below which residual multipole moments in ', &
                 'the Hartree potential will be set to zero: ', &
                 far_distance_hartree_multipole_moment_threshold
              write(use_unit,*) '| (Affects periodic systems only. See subroutine integrate_hartree_log_grid for comments.)'
            end if

            flag_multip_moments_threshold = .true.

         case('hartree_radius_threshold')
            read(inputline,*,end=88,err=99) desc_str, far_distance_hartree_radius_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of hartree real space radius:', &
               far_distance_hartree_radius_threshold
            end if

         case('adaptive_hartree_radius_th')
            read(inputline,*,end=88,err=99) desc_str, far_distance_adaptive_hartree_radius_threshold
            if (myid.eq.0) then
            write(use_unit,*) 'Threshold of adaptive hartree real scpace radius:', &
                  far_distance_adaptive_hartree_radius_threshold
            end if

         case('hartree_fourier_part_th')
            read(inputline,*,end=88,err=99) desc_str, &
               hartree_potential_fourier_part_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,A,E14.6)') &
               'Threshold of the Fourier part of the Hartree potential: ', &
               hartree_potential_fourier_part_threshold
            end if

         case('hartree_convergence_parameter','ewald_radius','Ewald_radius')
            read(inputline,*,end=88,err=99) desc_str, &
                  desc_str
            if (desc_str.eq.'automatic') then

              Ewald_radius_automatic = .true.

              if (myid.eq.0) then
                 write(use_unit,'(2X,A,A)') &
                       'Range separation parameter for Ewald summation (hartree_convergence_parameter) : ', &
                       "Determined automatically."
              end if

            else
              ! expect actual numerical value

              read(inputline,*,end=88,err=99) desc_str, &
                  Ewald_radius

              Ewald_radius_automatic = .false.

              ! VB: Safety net
              if (Ewald_radius.lt.1.d0) then
                 if (myid.eq.0) then
                    write(use_unit,'(1X,A,E14.6,A)') &
                       '* Attention: hartree_convergence_parameter choice ', &
                       Ewald_radius, ' bohr'
                    write(use_unit,'(1X,A)') &
                       '* is smaller than what we consider safely tested. If you wish to use this anyway,'
                    write(use_unit,'(1X,A)') &
                       '* go ahead with care by commenting out this warning and stop in the source code.'
                 end if
                 call aims_stop_coll('', func)
              else if (Ewald_radius.gt.5.d0) then
                 if (myid.eq.0) then
                    write(use_unit,'(1X,A,E14.6,A)') &
                       '* Attention: hartree_convergence_parameter choice ', &
                       Ewald_radius, ' bohr'
                    write(use_unit,'(1X,A)') &
                       '* is larger than what we consider safely tested. In essence, you must make sure that'
                    write(use_unit,'(1X,A)') &
                       '* the outermost radial integration shells _and_ integration partition tables'
                    write(use_unit,'(1X,A)') &
                       '* (defined implicitly through cut_free_atom) reach far enough out to capture the Ewald'
                    write(use_unit,'(1X,A)') &
                       '* compensating Gaussian charge density, whose extent is defined by hartree_convergence_parameter.'
                    write(use_unit,'(1X,A)') &
                       '* If you wish to use this anyway,'
                    write(use_unit,'(1X,A)') &
                       '* go ahead with care by commenting out this warning and stop in the source code.'
                 end if
                 call aims_stop_coll('', func)
              end if
              if (myid.eq.0) then
                 write(use_unit,'(2X,A,E14.6,A)') &
                       'Range separation parameter for Ewald summation (hartree_convergence_parameter) : ', &
                       Ewald_radius, "bohr."
              end if

            end if ! hartree_convergence parameter not 'automatic'

         case('fast_Ewald')
            read(inputline,*,end=88,err=99) desc_str, &
                  fast_Ewald
            if (fast_Ewald) then
               write(info_str,'(2X,A)') &
               "New faster Ewald potential for periodic systems requested."
            else
               write(info_str,'(2X,A)') &
               "fast_Ewald: Explicit request to use older Ewald potential evaluation."
            end if
            call localorb_info(info_str)

         case('hartree_fp_function_splines')
            read(inputline,*,end=88,err=99) desc_str, hartree_fp_function_splines
            if (hartree_fp_function_splines) then
               write(info_str,'(2X,A)') "Splining the Hartree Fp-functions during calculation."
            else
               write(info_str,'(2X,A)') "Explicitly NOT splining the Hartree Fp-functions during calculation."
            end if
            call localorb_info(info_str)

         case('multipole_threshold')
            read(inputline,*,end=88,err=99) desc_str, &
                  multipole_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,A,E14.6)') &
                     'Cutoff for q(lm)/r^(l+1) in long-range Hartree potential : ', &
                     multipole_threshold
            end if

            flag_multipole_threshold = .true.

         case('transport')
            flag_transport = .true.
            call transport_read_control_in(inputline)

         case('l_hartree_far_distance')
            read(inputline,*,end=88,err=99) desc_str, &
                  l_hartree_far_distance
            if (myid.eq.0) then
               write(use_unit,'(2X,A,E14.6)') &
                     'Multipole components q(lm)/r^(l+1) in long-range Hartree potential NOT computed above l= ', &
                     l_hartree_far_distance
            end if

            flag_l_hartree_far_distance = .true.


         case('use_hartree_non_periodic_ewald')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            select case(desc_str)

              case('.true.')
                use_hartree_non_periodic_ewald = .true.

              case('gridspacing')
                use_hartree_non_periodic_ewald = .true.
                read(inputline,*,end=88,err=99) desc_str, desc_str, gridwidth_in_pm ! actually in A
                gridwidth_in_pm  =  100 * gridwidth_in_pm

              case('.false.')
               call localorb_info( "  Switch 'use_hartree_non_periodic_ewald' is set to 'false' (default)." )
               call localorb_info( "  | Do not use Ewald's decomposition of the Hartree potential &
                                                                        &in the non-periodic case." )
              case default
                call localorb_info( "  | STOP: Wrong argument(s) for switch 'use_hartree_non_periodic_ewald'!" )
                call aims_stop_coll('', func)

            end select  ! case(desc_str)


            if ( use_hartree_non_periodic_ewald ) then
               call localorb_info( "  Switch 'use_hartree_non_periodic_ewald' is enabled." )
               call localorb_info( "  | Use Ewald's decomposition of the Hartree potential in the &
                                                                        &non-periodic case." )
               call localorb_info( "  | Warning: This feature is experimental!" )
               if ( n_periodic > 0 ) then
                  call localorb_info( "  | STOP: This feature works only in the non-periodic case!" )
                  call aims_stop_coll('', func)
               end if
            end if

         case('hartree_d_matrix_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case(desc_str)
              case('pseudo_inverse')
               hartree_d_matrix_method = 2
               call localorb_info( "  | Use pseudo inverse method for d matrix in Hartree potential ")
              case('submatrix')
               hartree_d_matrix_method = 1
               call localorb_info( "  | Use submatrix method for d matrix in Hartree potential ")
            end select  ! case(desc_str)


         case('set_vacuum_level')
! OTH: On suggestion of JW and VB:
           write(info_str,*) '************************************************************************************'
           call localorb_info(info_str)
           write(info_str,*) 'ATTENTION: Since the keyword set_vacuum_level is actually geometry-related'
           call localorb_info(info_str)
           write(info_str,*) 'it is moved now to the geometry.in file.'
           call localorb_info(info_str)
           write(info_str,*) 'We apologize for any inconvenience.'
           call localorb_info(info_str)
           write(info_str,*) '************************************************************************************'
           call localorb_info(info_str)
           write(info_str,*) 'Update: It is now no longer strictly necessary to set this keyword.'
           call localorb_info(info_str)
           write(info_str,*) 'In proper circumstances (i.e., when dipole correction or work function is requested,'
           call localorb_info(info_str)
           write(info_str,*) 'the vacuum level can be determined automatically.'
           call localorb_info(info_str)
           write(info_str,*) 'Alas, the functionality is still experimental, so be sure to check it every now and then'
           call localorb_info(info_str)
           write(info_str,*) '*********************************************************'
           call localorb_info(info_str)
           call aims_stop_coll ('Keyword moved',func)

            !read(inputline,*,end=88,err=99) desc_str, vacuum_z_level

            !vacuum_z_level = vacuum_z_level/bohr
            !flag_set_vacuum_level = .true.
!
!            if (myid.eq.0) then
!               write(use_unit,'(2X,A,E14.6)')'Vacuum level in z-axis : ',  vacuum_z_level*bohr
!            end if

         case('evaluate_work_function')
               read(inputline,*, iostat=i_code) desc_str, desc_str
               calculate_work_function = .true.
               if (i_code == 0 .and. desc_str /= '') then
                  read(inputline,*,end=88,err=99) &
                  & desc_str, calculate_work_function
               end if
            if (myid.eq.0) then
               if (calculate_work_function) write(use_unit,'(A)')'Evaluating the work function.'
            end if

         case('use_dipole_correction')
               read(inputline,*, iostat=i_code) desc_str, desc_str
               use_dipole_correction = .true.
               if (i_code == 0 .and. desc_str /= '') then
                  read(inputline,*,end=88,err=99) &
                  & desc_str, use_dipole_correction
               end if
            if (myid.eq.0) then
               if (use_dipole_correction) then
                  write(use_unit,'(2X,A)')'Using dipole correction.'
               else
                  write(use_unit,'(2X,A)')'Dipole correction explicitly switched off.'
               end if
            end if

         case ('evaluate_dipole_correction_method')
               read(inputline,*, iostat=i_code) desc_str, desc_str
               if(trim(desc_str).eq.'potential') then
                    dipole_correction_method='potential'
               elseif(trim(desc_str).eq.'dipole') then
                    dipole_correction_method='dipole'
               else
                    call aims_stop('Unknown dipole correctiod method')
               endif

         case ('nomad_comment')
               continue
         case ('overwrite_existing_cube_files')
               read(inputline,*,end=88,err=99) desc_str, overwrite_existing_cube_files
         case ('cube_default_size_safeguard')
               read(inputline,*,end=88,err=99) desc_str, cube_default_size_safeguard

         case ('output_cube_nth_iteration')
                 if (calculate_perturbative_soc) then
                   call localorb_info("")
                   write(info_str, "(2X,A)") '* The keyword output_cube_nth_iteration cannot be used with &
                                             & spin-orbit coupling, as SOC is applied after the SCF cycle.  Exiting.'
                   call aims_stop_coll(info_str, func)
                 endif

                 read(inputline,*,end=88,err=99) desc_str, cube_output_nth_iteration
                 out_cube_nth_iteration=.true.
                 write(use_unit,'(2X,A)') 'Cube files will be written during SCF'
                 write(use_unit,'(2X,A)') 'This feature is highly experimental and could'
                 write(use_unit,'(2X,A)') 'break anything in the code, including giving wrong'
                 write(use_unit,'(2X,A)') 'number here and there.'

         case ('use_old_cube_routine')
                 if (calculate_perturbative_soc) then
                   call localorb_info("")
                   write(info_str, "(2X,A)") '* The keyword use_old_cube_routine is not supported with &
                                             & spin-orbit coupling, as the old cube routines are considered &
                                             & legacy functionality.  Exiting.'
                   call aims_stop_coll(info_str, func)
                 endif

                 read(inputline,*,end=88,err=99) desc_str,use_old_cube_routine

         case ('cube_content_unit')
                 read(inputline,*,end=88,err=99) desc_str,cube_content_unit
                 if(trim(cube_content_unit).ne.'legacy'.and.trim(cube_content_unit).ne.'bohr') then
                    call aims_stop('Unknown cube_content_unit. Please correct.')
                 endif

         case ('distribute_leftover_charge')
                 read(inputline,*,end=88,err=99) desc_str,distribute_leftover_charge

            if (myid.eq.0) then
               if (distribute_leftover_charge) then

                  write(use_unit,'(1X,A)') &
                    '* Error. The "distribute_leftover_charge" option is no longer recommended.'
                  write(use_unit,'(1X,A)') &
                    '* A much better version is the "compensate_multipole_errors" keyword,'
                  write(use_unit,'(1X,A)') &
                    '* which is now recommended instead. '
                  write(use_unit,'(1X,A)') &
                    '* If you really insist on using "distribute_leftover_charge", the option '
                  write(use_unit,'(1X,A)') &
                    '* can be brought back by simply commenting out the "aims_stop" line below,'
                  write(use_unit,'(1X,A)') &
                    '* BUT if you do so, you should have a really good reason. Especially for'
                  write(use_unit,'(1X,A)') &
                    '* "light" integration grid settings, "distribute_leftover_charge" can cause'
                  write(use_unit,'(1X,A)') &
                    '* noticeable total energy inaccuracies and discontinuities.'

!                  write(use_unit,'(2X,A)') 'Small charge integration errors on the 3D grid will be omitted'
!                  write(use_unit,'(2X,A)') 'from the far-field multipole components by added small compensating'
!                  write(use_unit,'(2X,A)') 'charges.'
!                  write(use_unit,'(1X,A)') "*** Warning. As implemented right now, the 'distribute_leftover_charge'"
!                  write(use_unit,'(1X,A)') "*** option WILL introduce small discontinuitues in the energy landscape"
!                  write(use_unit,'(1X,A)') "*** as a function of the nuclear positions. Thus, especially for light"
!                  write(use_unit,'(1X,A)') "*** integration grids, any structure relaxation may show difficulties."
!                  write(use_unit,'(1X,A)') "*** Only use this option to obtain a well defined work function and"
!                  write(use_unit,'(1X,A)') "*** dipole correction in the vacuum of a slab - be careful otherwise."
!                  write(use_unit,'(1X,A)') "*** A better version may be forthcoming soon."

                  call aims_stop_coll ('Option no longer recommended.', func)
              end if
            end if

         case ('compensate_multipole_errors')
                 read(inputline,*,end=88,err=99) desc_str,compensate_multipole_errors

            if (myid.eq.0) then
               if (compensate_multipole_errors) then
                  write(use_unit,'(2X,A)') 'Charge integration errors on the 3D integration grid will be compensated'
                  write(use_unit,'(2X,A)') 'by explicit normalization and distribution of residual charges.'
                  write(use_unit,'(2X,A)') 'This choice is critical if you do DFPT calculations, since it is not yet'
                  write(use_unit,'(2X,A)') 'implemented in the DFPT routines (including friction and magnetic_response).'
               else
                  write(use_unit,'(2X,A)') 'Switching off the compensation of charge integration errors on the'
                  write(use_unit,'(2X,A)') '3D integration grid.'
               end if
            end if

            flag_compensate_multipole_errors = .true.

         case ('normalize_initial_density')
                 read(inputline,*,end=88,err=99) desc_str,normalize_initial_density

            if (myid.eq.0) then
               if (normalize_initial_density) then
                  write(use_unit,'(2X,A)') &
                  'Initial density will be explicitly normalized to have the right charge on the 3D integration grid.'
               else
                  write(use_unit,'(2X,A)') &
                  'Initial density will not be normalized to have the right charge on the 3D integration grid.'
                  write(use_unit,'(2X,A)') &
                  'The normalization used will pertain to free atoms normalized on a dense logarithmic grid.'
               end if
            end if

            flag_normalize_initial_density = .true.

         case('force_new_functional')
            read(inputline,*,end=88,err=99) desc_str, force_new_functional

            if (myid.eq.0)then
               if (force_new_functional)then
               write(use_unit,'(2X,A)') 'Using neutralized energy functional.'
               else
		if (n_periodic.eq.0) then
		    write(use_unit,'(2X,A)') 'Using old cluster energy functional.'
		else
		    write(use_unit,'(2X,A,A)') &
		    'force_new_functional: Cannot use ', &
		    'old cluster energy functional in periodic case.'
		    write(use_unit,'(2X,A)')' Ignoring this request.'
		end if
               end if
            end if

            if ( (n_periodic.gt.0) .and. (.not.(force_new_functional)) ) &
               then

               force_new_functional = .true.

            end if

         case('use_density_matrix_hf')
            use_density_matrix_hf = .true.
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
               'Using density matrix for exchange operator in HF.'
            end if

         case('compute_kinetic')
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') 'Kinetic Energy matrix will be computed.'
            end if
            flag_compute_kinetic = .true.

         case('Adams_Moulton_integrator')
            read(inputline,*,end=88,err=99) desc_str, &
                  Adams_Moulton_integrator
            if (myid.eq.0) then
               if (Adams_Moulton_integrator) then
                  write(use_unit,'(2X,A,A)') &
                  'Using Adams-Moulton linear multistep integrator in', &
                  'Hartree potential integrals.'
               else
                  write(use_unit,'(1X,A,A)') &
                  '* Warning - not using Adams Moulton linear multistep integrator', &
                  'in Hartree potential integrals.'
                  write(use_unit,'(1X,A,A)') &
                  '* Do this only for test purposes, ', &
                  'not for production.'
               end if
            end if

         case('density_update_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'orbital') then
               density_update_method = UPDATE_ORBITAL
            else if (desc_str.eq.'density_matrix') then
               density_update_method = UPDATE_DENS_MATRIX
            else if (desc_str .eq. 'split_update_methods') then
               density_update_method = UPDATE_SPLIT
            else if (desc_str.eq.'automatic') then
               density_update_method = UPDATE_AUTO
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Requested density update method ", &
                        "is not available."
               end if
               call aims_stop_coll('', func)
            end if
            flag_density_update_method = .true.

         case('collect_eigenvectors')
            read(inputline,*,end=88,err=99) desc_str, collect_eigenvectors

            ! record that the keyword 'collect_eigenvectors' was set specifically within control.in by the user.
            flag_collect_eigenvectors = .true.

         case('packed_matrix_format')

            flag_packed_matrix_format = .true.

            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'none') then
               packed_matrix_format = PM_none
            else if (desc_str.eq.'skyline') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                  & "Warning: Skyline code has been removed, use index."
               end if
               packed_matrix_format = PM_index
            else if (desc_str.eq.'index') then
               packed_matrix_format = PM_index
            end if
            if (myid.eq.0) then
               write(use_unit,*) 'Using packed matrix format', &
                     packed_matrix_format
            end if

         case('hartree_worksize')
            read(inputline,*,end=88,err=99) desc_str, hartree_worksize

            if (myid.eq.0) then
               write(use_unit,'(2X,A,E14.6,A)') &
               "Allowed work space size for Hartree potential: ", &
               hartree_worksize, " MB."
            end if

            flag_hartree_worksize = .true.

         case('hartree_partition_type')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."rho_r2") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in Hartree potential ", &
                        "calculations: rho / r^2"
               end if

               flag_hartree_partition_type = 1

            else if (desc_str.eq."rho_r") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in Hartree potential ", &
                        "calculations: rho / r"
               end if

               flag_hartree_partition_type = 2

            else if (desc_str.eq."rho") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in Hartree potential ", &
                        " calculations: rho"
               end if

               flag_hartree_partition_type = 3

            else if (desc_str.eq."fermi") then
   !            Fermi function
   !            1 / ( 1 + exp ( (r-p1)/p2 ) )
   !            where p1 and p2 are given in Angstroms

               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                     hartree_partition_parameters(1), &
                     hartree_partition_parameters(2)

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A,F10.5,A,F10.5,A)') &
                        "Partition function in Hartree potential ", &
                        "calculations: Fermi, radius = ", &
                        hartree_partition_parameters(1), " A, width = ", &
                        hartree_partition_parameters(2), " A."
               end if

               hartree_partition_parameters(1) = &
                  hartree_partition_parameters(1) / bohr

               hartree_partition_parameters(2) = &
                  hartree_partition_parameters(2) / bohr

               flag_hartree_partition_type = 4

            else if (desc_str.eq.'stratmann') then
               write(info_str,'(2X,2A)') 'Partition function in Hartree ', &
                     'potential is Stratmann.'
               call localorb_info(info_str)
               flag_hartree_partition_type = 5
            else if (desc_str.eq.'rho_r2_lda') then
               write(info_str,'(2X,2A)') 'Partition function in Hartree ', &
                     'potential: rho(LDA) / r^2 .'
               call localorb_info(info_str)
               flag_hartree_partition_type = 6
            else if (desc_str.eq.'stratmann_smooth') then
               write(info_str,'(2X,2A)') 'Partition function in Hartree ', &
                     'potential is Stratmann, but with a smooth decay to zero'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| at the outer radius of each atom as given by basis functions, free-atom potential etc.'
               call localorb_info(info_str)
               flag_hartree_partition_type = 7
            else if (desc_str.eq.'stratmann_smoother') then
               write(info_str,'(2X,2A)') 'Partition function in Hartree ', &
                     'potential is Stratmann, with an additionally restricted atom list.'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| At each point, only atoms whose outer radius of each atom as given by basis functions, '
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| free-atom potential etc. may contribute.'
               call localorb_info(info_str)
               flag_hartree_partition_type = 8
!LN: added for consistency
            else if (desc_str.eq.'stratmann_sparse') then
               write(info_str,'(2X,2A)') 'Partition function in Hartree ', &
                     'potential is Stratmann, but only stores the O(N) needed interatomic distances.'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| At each point, only atoms whose outer radius of each atom as given by basis functions, '
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| free-atom potential etc. may contribute.'
               call localorb_info(info_str)
               flag_hartree_partition_type = 9
            else

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Requested partition type for Hartree potential ", &
                        "does not exist!"
               end if
               call aims_stop_coll('', func)

            end if

         case('partition_type')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."rho_r2") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in integrals ", &
                        "calculations: rho / r^2"
               end if

               partition_type = 1

            else if (desc_str.eq."rho_r") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in integrals ", &
                        "calculations: rho / r"
               end if

               partition_type = 2

            else if (desc_str.eq."rho") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Partition function in integrals ", &
                        " calculations: rho"
               end if

               partition_type = 3

            else if (desc_str.eq."fermi") then


   !            Fermi function
   !            1 / ( 1 + exp ( (r-p1)/p2 ) )
   !            where p1 and p2 are given in Angstroms

               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                     integral_partition_parameters(1), &
                     integral_partition_parameters(2)

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A,F10.5,A,F10.5,A)') &
                        "Partition function in integrals: ", &
                        "Fermi, radius = ", &
                        integral_partition_parameters(1), " A, width = ", &
                        integral_partition_parameters(2), " A."
               end if

               integral_partition_parameters(1) = &
                  integral_partition_parameters(1) / bohr

               integral_partition_parameters(2) = &
                  integral_partition_parameters(2) / bohr

               partition_type = 4


             else if (desc_str.eq.'stratmann') then
               write(info_str,'(2X,2A)') 'Partition function in ', &
                     'integrals is Stratmann.'
               call localorb_info(info_str)
               partition_type = 5
            else if (desc_str.eq.'rho_r2_lda') then
               write(info_str,'(2X,2A)') 'Partition function in ', &
                     'integrals: rho(LDA) / r^2 .'
               call localorb_info(info_str)
               partition_type = 6
             else if (desc_str.eq.'stratmann_smooth') then
               write(info_str,'(2X,2A)') 'Partition function in ', &
                     'integrals is Stratmann, but with a smooth decay to zero'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') &
                     '| at the outer radius of each atom as given by basis functions, free-atom potential etc.'
               call localorb_info(info_str)
               partition_type = 7
             else if (desc_str.eq.'stratmann_smoother') then
               write(info_str,'(2X,2A)') 'Partition function in ', &
                     'integrals is Stratmann, with an additionally restricted atom list.'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| At each point, only atoms whose outer radius of each atom as given by basis functions, '
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| free-atom potential etc. may contribute.'
               call localorb_info(info_str)
               partition_type = 8
   !       LN: add new keyword to be able to compare results with new and old version of the atom_atom_tab
             else if (desc_str.eq.'stratmann_sparse') then
               write(info_str,'(2X,2A)') 'Partition function in ', &
                     'integrals is Stratmann, with only the O(N) relevant interatomic distances stored.'
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| At each point, only atoms whose outer radius of each atom as given by basis functions, '
               call localorb_info(info_str)
               write(info_str,'(2X,A)') '| free-atom potential etc. may contribute.'
               call localorb_info(info_str)
               partition_type = 9
            else

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Requested integral partition function ", &
                        "does not exist (yet)!"
               end if
               call aims_stop_coll('', func)

            end if

            flag_partition_type = .true.

         case('stratmann_a')
            read(inputline,*,end=88,err=99)  desc_str, stratmann_a
            write(info_str,'(2X,2A,F10.6)') 'Using Stratmann ', &
                  'partition tab parameter a = ', stratmann_a
            call localorb_info(info_str)

         case('relativistic')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."none") then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                        "Non-relativistic treatment of kinetic energy."
               end if

               flag_rel = REL_none
            else if (desc_str.eq."zora") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, zora_threshold
               if (desc_str.eq.'scalar') then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                           "Scalar relativistic treatment of kinetic ", &
                           "energy: ZORA."
                     write(use_unit,'(2X,A,E14.6)') &
                           "Threshold value for ZORA: ", &
                           zora_threshold
                  end if

                  flag_rel = REL_zora
               else if (desc_str.eq.'spinor') then
                  if (myid.eq.0) then
                        write(use_unit,'(2X,A,A)') &
                           "Relativistic treatment of kinetic ", &
                           "energy: ZORA."
                        write(use_unit,'(2X,A,E14.6)') &
                           "Threshold value for ZORA: ", &
                           zora_threshold
                  end if

                  flag_rel = REL_zora_spinor
               else

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A,A)') &
                     "ZORA relativistic treatment of kinetic energy: ", &
                           desc_str, ", unknown variant."
                  end if

                  call aims_stop_coll('', func)
               end if
            else if (desc_str.eq."atomic_zora") then
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               if (desc_str.eq.'scalar') then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Scalar relativistic treatment of kinetic ", &
                     'energy: on-site free-atom approximation to ZORA.'
                  end if

                  flag_rel = REL_atomic_zora
               else if (desc_str.eq.'spinor') then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Relativistic treatment of kinetic ", &
                     'energy: on-site free-atom approximation to ZORA.'
                  end if

                  flag_rel = REL_at_zora_spinor
               else

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A,A)') &
                     "Relativistic treatment of kinetic energy: ", &
                     desc_str, ", unknown type of treatment."
                  end if

                  call aims_stop_coll('', func)
               end if
            else if (desc_str.eq."own") then

               read(inputline,*,end=88,err=99) desc_str, desc_str

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Scalar relativistic treatment of kinetic ", &
                     'energy: Valence ZORA plus scalar core.'
                  end if

                  flag_rel = REL_own

            else if (desc_str.eq."kolning_harmon") then

               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, zora_threshold

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Kolning Harmon type of relativistic treatment"
                  end if

                  flag_rel = REL_KOLNING_HARMON
            else if (desc_str.eq."x2c") then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Relativistic treatment of kinetic ", &
                     'energy: X2C.'
                  end if

                  flag_rel = REL_x2c
            else if (desc_str.eq."4c_dks") then

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                     "Relativistic treatment of kinetic ", &
                     'energy: four-component DKS.'
                  end if

                  flag_rel = REL_4c_dks
            else

               if (myid.eq.0) then
                  write(use_unit,*) desc_str, " relativistic treatment is not ", &
                        "implemented at this time, feel free to change."
               end if

               call aims_stop_coll('', func)
            end if

         case('override_relativity')
            read(inputline,*,end=88,err=99) desc_str, override_relativity

            if (override_relativity) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "override_relativity: Explicitly overriding any built-in relativity checks."
                  write(use_unit,'(2X,A,A)') &
                  "If you use this flag, you should really know what you are doing."
               end if
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "override_relativity: Built-in relativity checks will be honored."
               end if
            end if

          case('override_warning_nobasis')
               override_warning_nobasis = .true.

          case('override_warning_libxc')
               override_warning_libxc = .true.

          case('override_warning_negative_nucleus')
              override_warning_negative_nucleus = .true.

         case('postprocess_anyway')
            read(inputline,*,iostat=info) desc_str, flag
            if (info == 0) then
               if (flag) then
                  postprocess_anyway = PP_ANYWAY_EVERYTHING
               else
                  postprocess_anyway = PP_ANYWAY_NOTHING
               end if
            else
               ! Retry to read integer
               read(inputline,*,end=88,err=99) desc_str, i_code
               if (i_code > 0) then
                  postprocess_anyway = PP_ANYWAY_EVERYTHING
               else
                  postprocess_anyway = PP_ANYWAY_NOTHING
               end if
            end if
            select case (postprocess_anyway)
            case(PP_ANYWAY_EVERYTHING)
               write(info_str,'(2X,A,A)') &
               "Postprocessing will be performed even if SCF does not converge."
            case(PP_ANYWAY_NOTHING)
               write(info_str,'(2X,A,A)') &
               "No postprocessing will be performed if SCF does not converge"
            case default
               call aims_stop('Invalid internal postprocess_anything flag', func)
            end select
            call localorb_info(info_str)

         case('use_small_component')
            read(inputline,*,end=88,err=99) desc_str, use_small_component

            if (myid.eq.0) then
               if (use_small_component)then
                  write(use_unit,*) 'Using the small component in KH core states'
               else
                  write(use_unit,*) 'NOT using the small component in KH core states'
               end if
            end if

         case('empty_states')
            read(inputline,*,end=88,err=99) desc_str, n_empty

            if (n_empty.le.0) then
               if (myid.eq.0) then
                  write (use_unit,'(1X,A,I8)') &
                     "Warning: Number of requested empty states per atom: ", n_empty
                  write (use_unit,'(1X,A)') &
                     "To properly determine the Fermi level, a number greater than zero is required!"
               end if

               call aims_stop_coll('', func)
            end if

            if (myid.eq.0) then
               write (use_unit,'(2X,A,I8)') &
                     "Number of empty states per atom: ", n_empty
            end if

         case('use_full_spectrum')
            read(inputline,*,end=88,err=99) desc_str
            if (myid.eq.0) then
              write (use_unit,'(2X,A,I8)') &
                   "Calculating the full spectrum for the basis set (i.e. includes all empty states.)"
            end if
            use_full_spectrum = .true.

         case('calculate_all_eigenstates')
            ! The difference between this keyword and use_full_spectrum is that this
            ! keyword does not assume that the user is using a correlated method
            read(inputline,*,end=88,err=99) desc_str
            if (myid.eq.0) then
              write (use_unit,'(2X,A,I8)') &
                   "Calculating all eigenstates for the basis set (i.e. includes all empty states.)"
            end if
            calculate_all_eigenstates = .true.

         case('full_embedding')
            read(inputline,*,end=88,err=99) desc_str, full_embedding

            if (myid.eq.0) then
               if (full_embedding) then
                  write (use_unit,'(2X,A,A)') &
                  "Embedding potential will change the self-consistent ", &
                  "charge density."
               else
                  write (use_unit,'(2X,A,A)') &
                  "An external embedding potential will not act on the ", &
                  "charge density. Calculating embedding energy only."
               end if
            end if

         case('KH_post_correction')
            read(inputline,*,end=88,err=99) desc_str, flag_KH_post_correction

         case('KS_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'lapack_old') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Kohn-Sham eigenvalues and eigenfunctions ", &
                        "calculated by LAPACK diagonalization."
               end if

               flag_KS = -1
            else if (desc_str.eq.'lapack_2010') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Kohn-Sham eigenvalues and eigenfunctions ", &
                        "calculated by LAPACK diagonalization (D&C)."
               end if

               use_lapack_fast = .true.
               flag_KS = -1
            else if ((desc_str.eq.'lapack_fast').or.(desc_str.eq.'serial')) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Kohn-Sham eigenvalues and eigenfunctions ", &
                        "calculated by LAPACK via ELSI."
               end if

               use_lapack_fast = .true.
               use_elsi = .true.
               flag_KS = -1
            else if (desc_str.eq.'svd') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Kohn-Sham eigenvalues and eigenfunctions ", &
                        "calculated by singular value decomposition."
               end if

               flag_KS = -2
            else if ((desc_str.eq.'lapack').or.(desc_str.eq.'scalapack')) then
               call localorb_multi(&
                    "*",&
                    "* The former keyword option 'KS_method lapack/scalapack' in control.in has been",&
                    "* changed to read 'KS_method serial/parallel'.  The reason is that these options",&
                    "* have indicated other linear algebra than custom (Sca)LAPACK in FHI-aims for a",&
                    "* long time.  For example, 'KS_method scalapack' actually employed the ELPA",&
                    "* eigenvalue solver, currently called through the ELSI infrastructure.  The new",&
                    "* naming scheme 'KS_method serial/parallel' indicates linear algebra default",&
                    "* options that may develop further over time, but these default options will not",&
                    "* necessarily be tied to one and the same library going forward.  Please update",&
                    "* your setting in control.in.",&
                    "*",&
                    "* Stopping the code. - Sorry about any inconvenience this may cause (hopefully a",&
                    "* one-time change only).",&
                    unit=use_unit)

               call aims_stop_coll('', func)
            else if (desc_str.eq.'scalapack_old') then
               write (info_str,'(2X,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by original ScaLAPACK (legacy!)."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               flag_KS = -1
            else if ((desc_str.eq.'scalapack_fast').or.(desc_str.eq.'elpa').or.&
                     (desc_str.eq.'elsi').or.(desc_str.eq.'parallel')) then
               write (info_str,'(2X,A)') &
                     "Kohn-Sham electronic structure solved by the ELSI library."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               use_elpa = .true.
               use_elsi = .true.
               flag_KS = -1
            else if (desc_str.eq.'elpa_2013') then
               write (info_str,'(2X,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by the ELPA eigensolver library (2013 version)."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               use_elpa = .true.
               flag_KS = -1
            else if (desc_str.eq.'scalapack_dc') then
               write (info_str,'(2X,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by the ELPA eigensolver library (DC version)."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               use_elpa = .true.
               flag_KS = -1
            else if (desc_str.eq.'lopcg') then
               write (info_str,'(2X,A,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by locally optimal preconditioned ", &
                     "conjugate gradient method."
               call localorb_info(info_str,use_unit,'(A)')

               use_cg = .true.
               flag_KS = -1

            else if (desc_str.eq.'scalapack+lopcg') then
               write (info_str,'(2X,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by ScaLAPACK and LOPCG method."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               use_cg = .true.
               flag_KS = -1

            else if (desc_str.eq.'elpa+lopcg') then
               write (info_str,'(2X,A,A)') &
                     "Kohn-Sham eigenvalues and eigenfunctions ", &
                     "calculated by ScaLAPACK (fast) and LOPCG method."
               call localorb_info(info_str,use_unit,'(A)')

               use_scalapack = .true.
               use_elpa = .true.
               use_cg = .true.
               flag_KS = -1

            else

               if (myid.eq.0) then
                  write(use_unit,*) 'KS_method: "', desc_str, '"', &
                  " is not a valid treatment", &
                  " type for the Kohn-Sham eigenvalue problem. "
               end if

               call aims_stop_coll('', func)
            end if

            flag_KS_method = .true.

         case('elpa_settings')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'auto') then

              write(info_str,'(2X,A)') &
                  "Settings for the ELPA eigensolver: automatic"
              call localorb_info(info_str,use_unit,'(A)')
              solver_method = 0
              elsi_elpa_solver = 2

            else if (desc_str.eq.'one_step_solver') then

              write(info_str,'(2X,A)') &
                  "Settings for the ELPA eigensolver: Always use the one-step solver."
              call localorb_info(info_str,use_unit,'(A)')
              solver_method = 1
              elsi_elpa_solver = 1

            else if (desc_str.eq.'two_step_solver') then

              write(info_str,'(2X,A)') &
                  "Settings for the ELPA eigensolver: Always use the two-step solver."
              call localorb_info(info_str,use_unit,'(A)')
              solver_method = 2
              elsi_elpa_solver = 2

            else

              write(info_str,'(1X,A)') &
                  "* Settings for the ELPA eigensolver: Unknown selection"
              call localorb_info(info_str,use_unit,'(A)')
              if (myid.eq.0) then
                write(use_unit,'(A)') inputline
              end if
              write(info_str,'(1X,A)') &
                  "* Please correct before continuing."
              call localorb_info(info_str,use_unit,'(A)')
              call aims_stop_coll('Unknown choice for "elpa_settings" keyword', func)

            end if

         case('override_illconditioning')
            read(inputline,*,end=88,err=99) desc_str, override_illconditioning

            if (override_illconditioning) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                  "override_illconditioning: Explicitly overriding any built-in checks for an ill-conditioned overlap matrix."
                  write(use_unit,'(2X,A)') &
                  "*** WARNING: If you use this flag, you should really know what you are doing."
                  write(use_unit,'(2X,A)') &
                  "*** DO NOT keep this flag set by default in all your control.in files."
               end if
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                  "override_illconditioning: Built-in checks for ill-conditioning will be honored."
               end if
            end if

         case('onsite_accuracy_threshold')
            read (inputline,*) desc_str, onsite_accuracy_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,A,F20.8)') 'onsite_accuracy_threshold set to ', onsite_accuracy_threshold, ' eV.'
            end if
            onsite_accuracy_threshold = onsite_accuracy_threshold/hartree

         case('override_integration_accuracy')
            read (inputline,*) desc_str, override_integration_accuracy
            if (myid.eq.0) then
               if (override_integration_accuracy) then
                  write(use_unit,'(2X,A)') &
                    'override_integration_accuracy: Will NOT stop regardless of integration accuracy concerns.'
               else
                  write(use_unit,'(2X,A)') &
                    'override_integration_accuracy: Code will stop if there are integration accuracy concerns.'
               end if
            end if

         case('condition_penalty')
            read (inputline,*) desc_str, condition_penalty
            if (myid == 0) then
               write(use_unit,'(2X,A,ES10.2,A)') &
               & "Set penalty for ill-conditioned basis functions to", &
               & condition_penalty, " Ha."
            end if

         case('initial_ev_solutions')
            read(inputline,*,end=88,err=99) desc_str, initial_ev_solutions

            write(info_str,'(2X,A,I4)') &
                  "Initial EV solutions before lopcg: ", &
                  initial_ev_solutions
            call localorb_info(info_str,use_unit,'(A)')

            flag_initial_ev_solutions = .true.

         case('max_lopcg_iterations')
            read(inputline,*,end=88,err=99) desc_str, max_cg_iterations

            write(info_str,'(2X,A,I4)') &
                  "Maximal number of iterations for lopcg: ", &
                  max_cg_iterations
            call localorb_info(info_str,use_unit,'(A)')

            flag_max_cg_iterations = .true.

         case('lopcg_tolerance')
            read(inputline,*,end=88,err=99) desc_str, lopcg_tol

            write(info_str,'(2X,A,E12.5)') &
                  "Tolerance for residual in lopcg: ", &
                  lopcg_tol
            call localorb_info(info_str,use_unit,'(A)')

            flag_lopcg_tolerance = .true.

         case('lopcg_start_tolerance')
            read(inputline,*,end=88,err=99) desc_str, lopcg_start_tol

            write(info_str,'(2X,A,A,E12.5)') &
                  "Starting lopcg when change in ", &
                  "sum of eigenvalues falls below: ", &
                  lopcg_start_tol
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_preconditioner')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."diagonal") then
               cg_preconditioner = prec_diagonal

               write(info_str,'(2X,A)') &
                     "Using diagonal preconditioner with lopcg."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq."ovlp_inverse") then
               cg_preconditioner = prec_inv_ovlp

               write(info_str,'(2X,A,A)') &
                     "Using inverse of the overlap matrix as a ", &
                     "preconditioner for lopcg."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq."ham_ovlp_inverse") then
               cg_preconditioner = prec_inv_ham_ovlp

               write(info_str,'(2X,A,A)') &
                     "Using inverse of the hamiltonian + overlap matrix as a ", &
                     "preconditioner for lopcg."
               call localorb_info(info_str,use_unit,'(A)')

            else

               write(info_str,'(1X,A)') &
                     "* Unknown preconditioner for lopcg."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)

            end if

            flag_lopcg_preconditioner = .true.

         case('lopcg_omega')
            read(inputline,*,end=88,err=99) desc_str, lopcg_omega
            write(info_str,'(2X,A,E11.4)') &
                  "Lopcg omega set to value: ", &
                  lopcg_omega
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_block_size')
            read(inputline,*,end=88,err=99) desc_str, lopcg_block_size
            write(info_str,'(2X,A,I6)') &
                  "Block size for lopcg set to: ", &
                  lopcg_block_size
            call localorb_info(info_str,use_unit,'(A)')

            flag_lopcg_block_size = .true.

         case('lopcg_use_all_evs')
            read(inputline,*,end=88,err=99) desc_str, lopcg_use_all_evs
            write(info_str,'(2X,A)') &
                  "Using all eigenvectors in lopcg regardless of block size."
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_skip_rate')
            read(inputline,*,end=88,err=99) desc_str, lopcg_skip_rate
            write(info_str,'(2X,A,I4,A)') &
                  "Skipping LOPCG every ", lopcg_skip_rate, " cycle."
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_slide_ev_block')
            read(inputline,*,end=88,err=99) desc_str, lopcg_slide_ev_block
            write(info_str,'(2X,A)') &
                  "Enabling sliding of blocks in LOPCG."
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_ritz_scalapack')
            read(inputline,*,end=88,err=99) desc_str, use_ritz_scalapack
            write(info_str,'(2X,A)') &
                  "Using ScaLAPACK for solving the Ritz problem."
            call localorb_info(info_str,use_unit,'(A)')

         case('lopcg_adaptive_tolerance')
            read(inputline,*,end=88,err=99) desc_str, lopcg_adaptive_tolerance
            write(info_str,'(2X,A)') &
                  "Using adaptice tolerance in lopcg."
            call localorb_info(info_str,use_unit,'(A)')

            flag_lopcg_block_size = .true.

         case('lopcg_auto_blocksize')
            read(inputline,*,end=88,err=99) desc_str, lopcg_auto_blocksize
            write(info_str,'(2X,A)') &
                  "Using automatic blocksize detection in lopcg."
            call localorb_info(info_str,use_unit,'(A)')

         case('basis_threshold')
            read(inputline,*,end=88,err=99) desc_str, basis_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,A,E11.4)') &
                     "Threshold for basis singularities: ", basis_threshold
            end if
            flag_basis_threshold = .true.

         case('tau_threshold')
         ! This is an attempt to try and assist meta-GGA convergence.
         ! When tau is very small we get singularities in the potential for
         ! the meta-GGA SCF, as most meta-GGAs contain a 1/tau term.
         ! Currently I'm trying this hard cutoff, which is manually applied
         ! when there is a convergence problem,  but what would be ideal is to
         ! have a damping function like in J. Chem. Phys. 127 214103 (2007)
         ! (doi:10.1063/1.2800011), however they way this paper is written makes
         ! me think that the damping is included in the functional derivation,
         ! rather than externally, because the damping is applied to f_xc[rho]
         ! instead of to the derivatives (Eqn. 6). It'd be great if someone who
         ! understands this could read and tell me if I'm interpreting things
         ! correctly! AJL/Dec2016
            read(inputline,*,end=88,err=99) desc_str, tau_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,A,E11.4)') &
                     "Threshold for tau cutoff (meta-GGA): ", tau_threshold
            end if

         case('prodbas_nb')
            if (use_prodbas) then
               read(inputline,*,end=88,err=99) desc_str, prodbas_nb
               if (myid.eq.0) then
                  if (prodbas_nb < 0) then
                     call aims_stop_coll('prodbas_nb < 0', func)
                  else if (prodbas_nb == 0) then
                     write(use_unit,'(2X,A)') &
                     & "Distribute product basis in chunks of automatic size."
                  else if (prodbas_nb <= 64) then
                     write(use_unit,'(2X,A,I6,A)') &
                     & "Distribute product basis in chunks of ", &
                     &  prodbas_nb, "."
                  else
                     write(use_unit,'(2X,A,I6,A)') &
                     & "** Distribute product basis in chunks of ", &
                     &  prodbas_nb, "."
                     write(use_unit,'(2X,A,I6,A)') &
                     & "** Be aware that values above 64 are not sensible."
                  end if
               end if
            end if

         case('prodbas_threshold')
            if (use_prodbas) then
               read(inputline,*,end=88,err=99) desc_str, prodbas_threshold

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,E11.4)') &
                     "Threshold for auxiliary basis singularities: ", &
                        prodbas_threshold
               end if

               flag_prodbas_threshold = .true.
            else
               call localorb_info('* Ignoring prodbas_threshold', use_unit, '(2X,A)')
            endif

         case('default_prodbas_acc')
            read(inputline,*,end=88,err=99) desc_str, default_prodbas_acc
            flag_default_prodbas_acc = .true.
         case('default_max_n_prodbas')
            read(inputline,*,end=88,err=99) desc_str, default_max_n_prodbas
            flag_default_max_n_prodbas = .true.
            call localorb_info('** Now ignoring max_n_prodbas.', use_unit, '(2X,A)')
         case('default_max_l_prodbas')
            read(inputline,*,end=88,err=99) desc_str, default_max_l_prodbas
            flag_default_max_l_prodbas = .true.

         case('use_smooth_prodbas_threshold')
            read(inputline,*,end=88,err=99) desc_str, use_smooth_prodbas_threshold
            if (myid.eq.0) then
               write(use_unit,'(2X,2A)') "Using smooth prodbas_threshold ", &
               & "for auxiliary basis singularities"
            end if

         case('use_asym_inv_sqrt')
            read(inputline,*,end=88,err=99) desc_str, use_asym_inv_sqrt
            if (myid.eq.0) then
               write(use_unit,'(2X,2A)') "Using asymmetric inversion ", &
               & "for the auxiliary basis Coulomb matrix"
            end if

         case('nbeads')         ! XZL: added for PIMD
            read(inputline,*,end=88,err=99) desc_str, n_beads
            if (myid.eq.0) then
               write(use_unit,'(2X,2A)') "Using a finite number of beads ", &
               & "for the PIMD simulation"
            end if

         case('mixer')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq."linear") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Using linear charge density mixing."
               end if

               mixer = MIX_LINEAR
            else if (desc_str.eq."pulay") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Using pulay charge density mixing."
               end if

               mixer = MIX_PULAY

            else if (desc_str.eq."broyden") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Using broyden charge density mixing."
               end if

               mixer = MIX_BROYDEN

               if (myid.eq.0) then
                 write (use_unit,'(1X,A)') &
                        "* Warning! Broyden mixer is under development!"
               end if

            else

               if (myid.eq.0) then
                  write (use_unit,'(1X,A,A,A)') &
                        "* Bogus mixer requested? Mixer ", &
                        desc_str, " is not defined - sorry!"
               end if

               call aims_stop_coll('', func)
            end if

            flag_mixer = .true.

         case('mixer_constraint')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq."linear") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Using linear charge density mixing ", &
                        "in spin constrained DFT."
               end if

               mixer_constraint = 0

            else if (desc_str.eq."pulay") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Using pulay charge density mixing ", &
                        "in spin constrained DFT."
               end if

               mixer_constraint = 1

            else if (desc_str.eq."bfgs") then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Applying BFGS-minimization routine ", &
                        "in spin constrained DFT."
               end if

               mixer_constraint = 2

            else

               if (myid.eq.0) then
                  write (use_unit,'(1X,A,A,A,A)') &
                        "* Bogus mixer requested in spin constrained DFT?", &
                        " Mixer ", desc_str, " is not defined - sorry !"
               end if

               call aims_stop_coll('', func)
            end if

            flag_mixer_constraint = .true.

         case('ini_linear_mixing')
            read(inputline,*,end=88,err=99) desc_str, ini_linear_mixing

            if (myid.eq.0) then
               write (use_unit,'(2X,A,A,I5)') &
                     "Number of initial iterations with linear ", &
                     "charge density mixing: ", &
                     ini_linear_mixing
            end if

            flag_ini_linear_mixing = .true.

         case('ini_linear_mix_param')
            read(inputline,*,end=88,err=99) desc_str, linear_mix_param(1)

            if (myid.eq.0) then
               write (use_unit,'(2X,A,A,E12.5)') &
                     "Charge density mixing parameter ", &
                     "for initial linear mixing: ", &
                     linear_mix_param(1)
            end if

            flag_linear_mix_param = .true.

         case('ini_spin_mix_param')
            if (n_spin .eq. 2) then


               read(inputline,*,end=88,err=99) desc_str, linear_mix_param(2)

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A,E12.5)') &
                        "Spin density mixing parameter ", &
                        "for initial linear mixing: ", &
                        linear_mix_param(2)
               end if

               flag_linear_spin_mix_param = .true.

            else

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A,A)') &
                        "Spin density mixing parameter ", &
                        "for initial linear mixing ", &
                        "ignored since no spin-polarized ", &
                        "calculation selected."

               end if

            end if

         case('ini_linear_mixing_constraint')
            read(inputline,*,end=88,err=99) desc_str, ini_linear_mixing_constraint

            if (myid.eq.0) then
               write (use_unit,'(2X,A,A,I5)') &
                     "Number of initial iterations with linear ", &
                     "charge density mixing in spin constrained DFT: ", &
                     ini_linear_mixing_constraint
            end if

            flag_ini_linear_mixing_constraint = .true.

         case('adjust_scf')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'once') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, adjust_scf_iteration

               if (adjust_scf_iteration.lt.0) then
                 if (myid.eq.0) then
                    write (use_unit,'(1X,A,I5,A)') &
                        "* Error: 'ajust_scf once' cannot be called with ", &
                      adjust_scf_iteration, " iterations."
                 end if
                 call aims_stop_coll('Please correct keyword adjust_scf.', func)
               else
                 if (myid.eq.0) then
                    write (use_unit,'(2X,A,I5,A)') &
                        "S.C.F. convergence parameters will be adjusted in iteration number ", &
                    adjust_scf_iteration, " of the first full s.c.f. cycle only."
                 end if
               end if
               adjust_scf = .true.
               adjust_scf_always = .false.

            else if (desc_str.eq.'always') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, adjust_scf_iteration

               if (adjust_scf_iteration.lt.0) then
                 if (myid.eq.0) then
                    write (use_unit,'(1X,A,I5,A)') &
                        "* Error: 'ajust_scf always' cannot be called with ", &
                      adjust_scf_iteration, " iterations."
                 end if
                 call aims_stop_coll('Please correct keyword adjust_scf.', func)
               else
                 if (myid.eq.0) then
                    write (use_unit,'(2X,A,I5,A)') &
                        "S.C.F. convergence parameters will be adjusted in iteration number ", &
                    adjust_scf_iteration, " of each new s.c.f. cycle."
                 end if
               end if
               adjust_scf = .true.
               adjust_scf_always = .true.

            else if (desc_str.eq.'never') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                      "S.C.F. convergence parameters will never be adjusted atomatically."
               end if
               adjust_scf = .false.
               adjust_scf_always = .false.

            else

               if (myid.eq.0) then
                 write (use_unit,'(1X,A,A,A)') &
                     "* Error: 'ajust_scf always' cannot be called with option ", &
                   desc_str, "."
               end if
               call aims_stop_coll('Please correct keyword adjust_scf.', func)

            end if

         case('charge_mix_param')
            read(inputline,*,end=88,err=99) desc_str, charge_mix_param(1)

            if (myid.eq.0) then
               write (use_unit,'(2X,A,F10.4)') &
                     "Charge density mixing - mixing parameter: ", &
                     charge_mix_param(1)
            end if

            if ( (charge_mix_param(1).le.0.d0) .or. (charge_mix_param(1).gt.1.d0) ) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                    "* Warning: The chosen 'charge_mix_param' value is not between zero or one."
                  write(use_unit,'(1X,A)') &
                    '* The value found in control.in will likely lead to no s.c.f. convergence: '
                  write(use_unit,'(1X,A)') &
                    inputline
                  write(use_unit,'(1X,A)') &
                    '* Please check. '
               end if

               warn_mixing_parameter_range = .true.

            end if

            flag_mix_param = .true.

            charge_mix_param_set = .true.

         case('prec_mix_param')
            read(inputline,*,end=88,err=99) desc_str, prec_mix_param(1)

            if (myid.eq.0) then
               write (use_unit,'(2X,A,F10.4)') &
                     "Charge density mixing - mixing parameter: ", &
                     prec_mix_param(1)
            end if

            flag_prec_mix_param = .true.

         case('spin_mix_param')
            if (n_spin .eq. 2) then

               read(inputline,*,end=88,err=99) desc_str, charge_mix_param(2)

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,F10.4)') &
                        "Spin density mixing - mixing parameter: ", &
                        charge_mix_param(2)
               end if

               flag_spin_mix_param = .true.

               charge_mix_param_set = .true.

            else
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A,A)') &
                        "Spin density mixing parameter ", &
                        "ignored since no spin-polarized ", &
                        "calculation selected."

               end if
            end if

         case('switch_external_pert')
            read(inputline,*,end=88,err=99) desc_str, ext_cycle, ext_safe

            if (myid.eq.0) then
               write (use_unit,'(2X, A, I10)') &
                     "External perturbation switched on up to cycle : ", &
                     ext_cycle
            end if

            flag_atom_ref = .true.

         case('scs_mp2_parameters')
            read(inputline,*,end=88,err=99) desc_str, pt_mp2, ps_mp2

            use_scsmp2 = .true.

            if (myid.eq.0) then
               write (use_unit,'(2X, A, F10.6, F10.6)') &
            "SCS MP2 approximation with parallel/opposite components ", &
                     pt_mp2, ps_mp2
            end if

            flag_scsmp2 =.true.

         case('frozen_core_scf')
            read(inputline,*,end=88,err=99) desc_str, frozen_core_scf

         case('frozen_core_scf_cutoff')
            read(inputline,*,end=88,err=99) desc_str, frozen_core_scf_cutoff
            frozen_core_scf_cutoff = frozen_core_scf_cutoff/hartree

         case('frozen_core_scf_factor')
            read(inputline,*,end=88,err=99) desc_str, frozen_core_scf_factor

         case('frozen_core')
            allocate (n_fc_shell(n_species),stat=i_spin)
            call check_allocation(i_spin, 'n_fc_shell                    ')

            n_fc_shell=0

            read(inputline,*,end=88,err=99) desc_str, i_start_mp2
            if (i_start_mp2.gt.1) then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Frozen core treatment."
               endif

               flag_frozen_core= .true.

               n_low_state = i_start_mp2
            endif

         case('n_max_pulay')
            read(inputline,*,end=88,err=99) desc_str, n_max_pulay

            if (myid.eq.0) then
               write (use_unit,'(2X,A,I4)') &
                     "Pulay mixing - number of memorized iterations: ", &
                     n_max_pulay
            end if

            if (n_max_pulay.lt.1) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                    "* Error: Found zero or fewer density iterations in Pulay mixing."
                  write(use_unit,'(1X,A)') &
                    '* Requested n_max_pulay value in control.in is invalid: '
                  write(use_unit,'(1X,A)') &
                    inputline
                  write(use_unit,'(1X,A)') &
                    '* Please correct. '
               end if

               call aims_stop_coll('', func)

            end if

            flag_pulay_iter = .true.

         case('n_max_pulay_constraint')
            read(inputline,*,end=88,err=99) desc_str, n_max_pulay_constraint

            if (myid.eq.0) then
               write (use_unit,'(2X,A,A,I4)') &
                     "Pulay mixing for constrained potential ", &
                     "- number of memorized iterations: ", &
                     n_max_pulay_constraint
            end if

            flag_pulay_iter_constraint = .true.

         case('n_max_broyden')
            read(inputline,*,end=88,err=99) desc_str, n_max_broyden

            if (myid.eq.0) then
               write (use_unit,'(2X,A,I4)') &
                     "Broyden mixing - number of memorized iterations: ", &
                     n_max_broyden
            end if

            if (n_max_broyden.lt.1) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                    "* Error: Found zero or fewer density iterations in Pulay mixing."
                  write(use_unit,'(1X,A)') &
                    '* Requested n_max_broyden value in control.in is invalid: '
                  write(use_unit,'(1X,A)') &
                    inputline
                  write(use_unit,'(1X,A)') &
                    '* Please correct. '
               end if

               call aims_stop_coll('', func)

            end if

            flag_broyden_iter = .true.

         case('relative_fp_charge_mix')
            read(inputline,*,end=88,err=99) desc_str, relative_fp_charge_mix

            if (myid.eq.0) then
               write (use_unit,'(2X,A,E14.6)') &
                     "Broyden mixing - Relative fixed point under-relaxation: ", &
                     relative_fp_charge_mix
            end if

            flag_relative_fp_charge_mix = .true.

         case('mixer_swap_boundary')
            read(inputline,*,end=88,err=99) desc_str, i_code

            call localorb_info('*** Ignoring mixer_swap_boundary as the features is removed.')
            call localorb_info('*** Please contact us if you really need it.')

         case('mixer_threshold')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."charge") then

               read(inputline,*,end=88,err=99) desc_str, desc_str, max_rho_change(1)

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                  "Found threshold value for charge density mixing step:", &
                  max_rho_change(1), " ."
                  write (use_unit,'(2X,A)') &
               "Charge density steps with a larger norm will be truncated."
               end if

               flag_mixer_threshold_charge = .true.

            else if (desc_str.eq."spin") then

               if (n_spin.eq.2) then

               read(inputline,*,end=88,err=99) desc_str, desc_str, max_rho_change(2)

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,E14.6,A)') &
                     "Found threshold value for spin density mixing step:", &
                     max_rho_change(2), " ."
                     write (use_unit,'(2X,A)') &
               "Spin density steps with a larger norm will be truncated."
                  end if

                  flag_mixer_threshold_spin = .true.
               else

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                     "Threshold value for spin density mixing step ", &
                     "will be ignored - unpolarized calculation."
                  end if

               end if

            else

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A)') &
                  "* mixer_threshold: Unknown mixer threshold tag ", &
                  desc_str, "."
                  write(use_unit,'(1X,A)') &
                  '* Supported choices are "charge" or "spin".'
               end if

               call aims_stop_coll('', func)

            end if

         case('sc_iter_limit')
            read(inputline,*,end=88,err=99) desc_str, sc_iter_limit

            if (myid.eq.0) then
               write (use_unit,'(2X, A, I5)') &
                     "Maximum number of s.-c. iterations  : ", &
                     sc_iter_limit
            end if

            flag_sc_max = .true.

         case('sc_init_iter')

            read(inputline,*,end=88,err=99) desc_str, sc_init_iter

            if (sc_init_iter.le.1) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A)') &
                  "* sc_init_iter: The specified value ", &
                  desc_str, " is not feasible."
                  write(use_unit,'(1X,A)') &
                  '* At least a single scf iteration must be completed for any total energy calculation.'
               end if

               call aims_stop_coll('', func)
            end if

            if (myid.eq.0) then
               write (use_unit,'(2X, A, I5)') &
                     "Number of s.-c. iterations prior to a full restart of mixing : ", &
                     sc_init_iter
            end if

            scf_restarted = .false.

            flag_sc_init_iter = .true.

         case('sc_init_factor')

            read(inputline,*,end=88,err=99) desc_str, sc_init_factor

            if (sc_init_factor.le.1.d0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A)') &
                  "* sc_init_factor: The specified value ", &
                  desc_str, " is likely not safe (at least untested)."
                  write(use_unit,'(1X,A,A)') &
                  "* Please remove the following stop in read_control.f90 by hand", &
                  " and recompile if you wish to proceed."
               end if

               call aims_stop_coll('', func)
            end if

            if (myid.eq.0) then
               ! This relates to sc_init_iter keyword
               write (use_unit,'(2X, A, I5)') &
                     "SCF density residual factor below which no scf restart will be triggered : ", &
                     sc_init_factor
            end if

            flag_sc_init_factor = .true.

         case('sc_abandon_etot')

            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."never") then
              check_sc_abandon_etot = .false.

              if (myid.eq.0) then
                write (use_unit,'(2X, A)') &
                     "The s.c.f. cycle will not be abandoned even for very large jumps between iterations."
              end if

            else
              read(inputline,*,end=88,err=99) desc_str, sc_abandon_etot_iter, sc_abandon_etot_threshold

              ! Syntax: sc_abandon_etot <iterations> <threshold>
              !         The threshold is in eV/atom !

              if (myid.eq.0) then
                write (use_unit,'(2X, A, I5, A)') &
                     "The s.c.f. cycle will be abandoned if the total energy change in ", &
                     sc_abandon_etot_iter, " successive iterations"
                write (use_unit,'(2X, A, F20.8, A)') &
                     "exceeded ", &
                     sc_abandon_etot_threshold, " eV/atom ."
              end if

              check_sc_abandon_etot = .true.

            end if

            flag_sc_abandon_etot = .true.

         case('constraint_it_lim')
            read(inputline,*,end=88,err=99) desc_str, constraint_it_lim

            write (info_str,'(2X, A, I5)') &
               "Maximum number of constraint iterations  : ", &
               constraint_it_lim
            call localorb_info(info_str,use_unit,'(A)')

            flag_constraint_iter= .true.

         case('constraint_mix')
            read(inputline,*,end=88,err=99) desc_str,constraint_mix(1),constraint_mix(2)

            write (info_str,'(2X,A,F10.5,1X,F10.5)') &
               "Mixing parameters for constraint potential : ", &
               constraint_mix(1),constraint_mix(2)

            call localorb_info(info_str,use_unit,'(A)')

            flag_constraint_mix= .true.

         case('constraint_precision')
            read(inputline,*,end=88,err=99) desc_str,constraint_precision

            write (info_str,'(2X,A,1X,E14.6)') &
               "Convergence criterion for constraint_electrons : ", &
               constraint_precision
            call localorb_info(info_str,use_unit,'(A)')

            flag_constraint_precision= .true.

         case('constraint_debug')
            read(inputline,*,end=88,err=99) desc_str,constraint_debug

            if (constraint_debug) then
               write (info_str,'(2X,A)') &
               "Extra output requested for locally constrained DFT."
            else
               write (info_str,'(2X,A)') &
               "No extra output requested for locally constrained DFT."
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('constraint_electrons')
            j_region=j_region+1
            read(inputline,*,end=88,err=99) desc_str,i_region

            write(info_str,'(2X,A,I5,A)') &
               "Found constraint for subsystem ", i_region, ":"
            call localorb_info(info_str,use_unit,'(A)')

            if (i_region.gt.n_atoms) then
               write(info_str,'(1X,A,A)') &
                  "* Label of subsystem must not be", &
                  " larger than total number of atoms."
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,'(1X,A)') "* Please correct."
               call localorb_info(info_str,use_unit,'(A)')

               call aims_stop_coll('', func)
            endif

            if (flag_region(i_region))then

               write(info_str,'(1X,A,A,I5,A)') &
                  "* Warning: constraint_electrons doubly specified", &
                  " for subsystem ",i_region,"."
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,'(1X,A,A,I5,A)') &
                  "* Overriding previous settings ", &
                  "for subsystem ",i_region,"."
               call localorb_info(info_str,use_unit,'(A)')

            else

               constraint_region_number(i_region)=j_region
               constraint_region_label(j_region)=i_region

               ! initialize stored constraint potential to zero here
               constraint_potential(j_region,1:n_spin) = 0.d0

               flag_region(i_region)=.true.

            end if


            read(inputline,*,end=88,err=99) desc_str,i_region, &
            ( constraint_electrons &
                  (constraint_region_number(i_region),i_spin), &
                  i_spin = 1, n_spin, 1 )

            do i_spin = 1, n_spin, 1

               if ( constraint_electrons &
                     (constraint_region_number(i_region),i_spin) &
                     .lt.0.d0 ) then

                  write(info_str,'(1X,A,A,I3,A)') &
                     "* Negative number of electrons specified for spin ", &
                     "channel ", i_spin, "."
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(1X,A)') &
                     "* Impossible - please correct."
                  call localorb_info(info_str,use_unit,'(A)')
                  call aims_stop_coll('', func)

               else

                  write(info_str,'(2X,A,I5,A,F10.5)') &
                  "| Electrons in spin channel ", i_spin, ": ", &
                     constraint_electrons &
                     ( constraint_region_number(i_region),i_spin)
                  call localorb_info(info_str,use_unit,'(A)')

               end if

            enddo
   ! wrapper to transform the multiplicity flag into the spin_constraint flag
         case('multiplicity')
            if (n_spin.eq.1) then

               write(info_str,'(1X,A)') &
                  '* Found multiplicity request, but this is a non-spinpolarized'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(1X,A)') &
                  '* calculation. Ignoring request.'
               call localorb_info(info_str,use_unit,'(A)')

            else if (n_periodic.ne.0) then
               fixed_spin_moment = .true.
               read(inputline,*,end=88,err=99)desc_str, hf_multiplicity
               write(info_str,'(2X,A,I4,A)') &
                     'Found multiplicity = ', hf_multiplicity, &
                     ' for Hartree-Fock and DFT calculations in periodic calculation.'
               call localorb_info(info_str,use_unit,'(A)')
!              write(info_str,'(1X,A)') &
!                '* Found multiplicity request. Unfortunately, this is not yet'
!              call localorb_info(info_str,use_unit,'(A)')
!              write(info_str,'(1X,A)') &
!                '* supported in periodic boundary conditions.'
!              call localorb_info(info_str,use_unit,'(A)')
!              write(info_str,'(1X,A)') &
!                "* It's high on the to-do list, please contact us if needed."
!              call localorb_info(info_str,use_unit,'(A)')

!              call aims_stop_coll('', func)

            else
               ! Spin is requested, thus multiplicity makes sense, cluster case

               read(inputline,*,end=88,err=99)desc_str, hf_multiplicity
               write(info_str,'(2X,A,I4,A)') &
                     'Found multiplicity = ', hf_multiplicity, &
                     ' for Hartree-Fock and DFT calculations.'
               call localorb_info(info_str,use_unit,'(A)')

               ! this is handled through the constraint infrastructure. FSM would work as well.
               ! In this case we need a Mulliken analysis later for trivial reasons (unshift eigenvalues)
               flag_run_mulliken = .true.

               j_region=j_region+1
               i_region=1

               write(info_str,'(2X,A,I5,A)') &
                  "Found constraint for subsystem ", i_region, ":"
               call localorb_info(info_str,use_unit,'(A)')

               if (i_region.gt.n_atoms) then
                  write(info_str,'(1X,A,A)') &
                     "* Label of subsystem must not be", &
                     " larger than total number of atoms."
                  call localorb_info(info_str,use_unit,'(A)')

                  write(info_str,'(1X,A)') "* Please correct."
                  call localorb_info(info_str,use_unit,'(A)')

                  call aims_stop_coll('', func)
               endif

               if (flag_region(i_region))then

                  write(info_str,'(1X,A,A,I5,A)') &
                     "* Warning: constraint_electrons doubly specified", &
                     " for subsystem ",i_region,"."
                  call localorb_info(info_str,use_unit,'(A)')

                  write(info_str,'(1X,A,A,I5,A)') &
                     "* Overriding previous settings ", &
                     "for subsystem ",i_region,"."
                  call localorb_info(info_str,use_unit,'(A)')

               else

                  constraint_region_number(i_region)=j_region
                  constraint_region_label(j_region)=i_region

                  ! initialize stored constraint potential to zero here
                  constraint_potential(j_region,1:n_spin) = 0.d0

                  flag_region(i_region)=.true.

               end if

            end if

         case('distributed_spline_storage')
            read(inputline,*,end=88,err=99) desc_str, use_distributed_spline_storage
            flag_distributed_spline_storage = .true.

            if (use_distributed_spline_storage .and. .not.use_mpi) then
               use_distributed_spline_storage = .false.
               call localorb_info('  * No distributed spline storage on serial run.')
            else
               if (myid.eq.0) then
                  if (use_distributed_spline_storage) then
                     write (use_unit,'(2X,A)') &
                     "Requested to store only required components of Hartree potential."
                  else
                     write (use_unit,'(2X,A,A)') &
                     "Full copy of Hartree potential on each thread ", &
                     "requested explicitly."
                  end if
               end if
            end if

         case('communication_type')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."standard") then
               communication_type = standard_mpi
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Standard (broadcast) MPI communication ", &
                  "requested explicitly."
               end if
            else if (desc_str.eq."bluegene") then
               communication_type = BlueGene_mpi
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "BlueGene-type (send/receive) MPI communication ", &
                  "requested explicitly."
               end if
            else if (desc_str.eq."cray") then
               communication_type = Cray_mpi
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Cray-type (receive/send) MPI communication ", &
                  "requested explicitly."
               end if
            else if (desc_str.eq."mpi_1coeff") then
               communication_type = mpi_1coeff
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Broadcast MPI communication for only 1st coeff ", &
                  "requested explicitly."
               end if
            else if (desc_str.eq."calc_hartree") then
               communication_type = calc_hartree
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "New calculation of Hartree spline coefficients ", &
                  "requested explicitly."
               end if
            else if (desc_str.eq."shmem") then
               call aims_shm_implemented(info)
               if (info == 0) then
                  if (myid.eq.0) then
                     write(use_unit,"(2X,A)") "* Shared memory not compiled in."
                     write(use_unit,"(2X,A)") "* Ignoring 'communication_type shmem'."
                     write(use_unit,"(2X,A)") "* Setting to calc_hartree."
                  end if
                  communication_type = calc_hartree
               else
                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                     "Shared memory storage of spline coefficients ", &
                     "requested explicitly."
                  end if
                  communication_type = shmem_comm
               end if
            else
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A,A)') &
                  "Error: Unknown MPI communication type ", &
                  communication_type, &
                  "requested explicitly."
               end if
               call aims_stop_coll('', func)
            end if
            if(communication_type /= calc_hartree .and. communication_type /= shmem_comm) then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') '* Warning: The requested communication type has been removed '// &
                                     'since it is inefficient'
                  write (use_unit,'(2X,A)') '* Using calc_hartree instead'
               endif
               communication_type = calc_hartree
            endif

            flag_communication_type = .true.

         case('occupation_acc')
            read(inputline,*,end=88,err=99) desc_str, occupation_acc

            if (myid.eq.0) then
               write (use_unit,'(2X, A, E11.4)') &
                     "Accuracy of occupation numbers  : ", &
                     occupation_acc
            end if

            flag_occupation_acc = .true.

         case('occupation_thr')
            read(inputline,*,end=88,err=99) desc_str, occupation_thr

            if (myid.eq.0) then
               write (use_unit,'(2X, A, E11.4)') &
               "Threshold below which occupation numbers are zeroed  : ", &
               occupation_thr
            end if

            flag_occupation_thr = .true.

         case('mu_determination_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq."zeroin") then
               mu_method = 0

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Chemical potential determination method  :  zeroin."
               end if
            else if (desc_str.eq."bisection") then
               mu_method = 1

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Chemical potential determination method  :  bisection."
               end if
            else
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)') &
                  "Error: Unknown chemical potential determination method."
               end if

               call aims_stop_coll('', func)
            end if

            flag_mu_method = .true.

         case('fermi_acc')
            read(inputline,*,end=88,err=99) desc_str, fermi_acc

            if (myid.eq.0) then
               write (use_unit,'(2X, A, E11.4)') &
                     "Accuracy of fermi level  : ", &
                     fermi_acc
            end if

            flag_fermi_acc = .true.

         case('max_zeroin')
            read(inputline,*,end=88,err=99) desc_str, max_zeroin

            if (myid.eq.0) then
               write (use_unit,'(2X, A, I4)') &
                     "Maximum # of zeroin steps  : ", &
                     max_zeroin
            end if

            flag_max_zeroin = .true.

         case('set_blacsdim')
            read(inputline,*,end=88,err=99) desc_str, desc_str, blacsdim
            if (desc_str.eq.'true') then
               set_blacsdim = .true.
            else
               set_blacsdim = .false.
               blacsdim = 1
            endif
            if (myid.eq.0) then
               print *,"Blacsdim in evalulate_crpa_energy_kspace will be set to ",blacsdim
            endif

         case('set_blacsdim_cpt2')
            read(inputline,*,end=88,err=99) desc_str, desc_str, blacsdim_cpt2
            if (desc_str.eq.'true') then
               set_blacsdim_cpt2 = .true.
            else
               set_blacsdim_cpt2 = .false.
               blacsdim_cpt2 = 1
            endif
            if (myid.eq.0) then
               print *,"Blacsdim in evalulate_cpt2_energy_kspace will be set to ",blacsdim_cpt2
            endif
         case('set_pair_block_cpt2')
            read(inputline,*,end=88,err=99) desc_str, desc_str, pair_block_cpt2
            if (desc_str.eq.'true') then
               set_pair_block_cpt2 = .true.
            else
               set_pair_block_cpt2 = .false.
               pair_block_cpt2 = 1
            endif
            if (myid.eq.0) then
               print *,"pair_block in evalulate_cpt2_energy_kspace will be set to ",pair_block_cpt2
            endif

         case('use_threadsafe_gwinit')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'true') then
               use_threadsafe_gwinit = .true.
            else
               use_threadsafe_gwinit = .false.
            endif
            if (myid.eq.0) then
               print *,"Threadsafe variant of gwinit will be used"
            endif


         case('occupation_type')
            read(inputline,*,end=88,err=99) desc_str, desc_str, occupation_width

            if (occupation_width.le.0.d0) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A)') &
                  "* Warning. Zero or negative width found for occupation_type keyword ", &
                  desc_str, "."
                  write(use_unit,'(1X,A)') &
                  "* This choice will lead to numerical problems, including possible divisions by zero."
                  write(use_unit,'(1X,A)') &
                  "* Please choose a finite (if narrow) occupation_width parameter."
                  write(use_unit,'(1X,A)') &
                  "* "
                  write(use_unit,'(1X,A)') &
                  "* If you are trying to avoid fractional occupation numbers for degenerate ground states:"
                  write(use_unit,'(1X,A)') &
                  "* If there _is_ a symmetry-breaking, non-degenerate ground state, a very small width"
                  write(use_unit,'(1X,A)') &
                  "* (0.00001 or similar) will often help. If this does not solve the problem:"
                  write(use_unit,'(1X,A)') &
                  "* Fractional occupation numbers for degenerate ground states can sometimes be avoided."
                  write(use_unit,'(1X,A)') &
                  "* by placing explicit constraints, using the force_occupation_* infrastructure."
                  write(use_unit,'(1X,A)') &
                  "* We realize that this is tricky and requires a good physical understanding of what the"
                  write(use_unit,'(1X,A)') &
                  "* target state should be. But remember that fractional occupation numbers are, fundamentally,"
                  write(use_unit,'(1X,A)') &
                  "* an allowed feature of generalized Kohn-Sham schemes. It is the density and total energy that"
                  write(use_unit,'(1X,A)') &
                  "* matters in the (hypothetical) exact theory, not the fictitious effective single particle orbitals."
               end if

               call aims_stop_coll('', func)

            end if

            if (desc_str.eq.'gaussian') then
               occupation_type = 0

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                        "Occupation type: Gaussian broadening, width = ", &
                        occupation_width, " eV."
               end if

               occupation_width = occupation_width / hartree
            else if (desc_str.eq.'fermi') then
               occupation_type = 1

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                        "Occupation type: Fermi broadening, width = ", &
                        occupation_width, " eV."
               end if

               occupation_width = occupation_width / hartree
            else if (desc_str.eq.'methfessel-paxton') then
               occupation_type = 2
               read(inputline,*,end=88,err=99) desc_str, desc_str, occupation_width, &
                     n_methfessel_paxton

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A,I5)') &
                  "Occupation type: Methfessel-Paxton broadening, width = ", &
                  occupation_width, " eV, N =", &
                  n_methfessel_paxton
               end if

               occupation_width = occupation_width / hartree
           else if (desc_str.eq.'integer') then
               occupation_type = 3
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Occupation type: Integer occupation."
               end if

               occupation_width = occupation_width / hartree
           else if (desc_str.eq.'cubic') then
               occupation_type = 4
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                  "Occupation type: Cubic polynomial broadening, width = ", &
                  occupation_width, " eV."
               end if

               occupation_width = occupation_width / hartree
           else if (desc_str.eq.'cold') then
               occupation_type = 5
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,E14.6,A)') &
                  "Occupation type: Marzari-Vanderbilt (cold) broadening, width = ", &
                  occupation_width, " eV."
               end if

               occupation_width = occupation_width / hartree
           else

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A)') &
                  "* Unknown type of orbital occupation (smearing) - ", &
                  desc_str, "."
                  write(use_unit,'(1X,A)') &
                  "* Please fix in control.in."
               end if

               call aims_stop_coll('', func)
            end if

            flag_occupation_type = .true.

            occupation_width_set = .true.

         case('partition_acc')
            read(inputline,*,end=88,err=99) desc_str, partition_acc

            if (myid.eq.0) then
               write (use_unit,'(2X,A,E11.4)') &
                     "Lower limit of norm for partition fn: ", &
                     partition_acc
            end if

            if (partition_acc.gt.1.d-12) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A)') &
                        "* Warning - partition fn ", &
                        "accuracy exceeds 1.d-12! "
                  write(use_unit,'(1X,A,A,A,A)') &
                        "* This is an untested accuracy ", &
                        "modification, so do not ", &
                        "disable this warning unless ", &
                        "you know what you're doing."
               end if

               call aims_stop_coll('', func)
            end if

            flag_partition_acc = .true.

         case('wave_threshold')
            read(inputline,*,end=88,err=99) desc_str, wave_threshold

            if (myid.eq.0) then
               write (use_unit,'(2X,A,E11.4)') &
                     "Threshold value for wave function integration: ", &
                     wave_threshold
               end if

            if (wave_threshold.gt.1.d-5) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A)') &
                  "* Warning - threshold for evaluation of u(r) in integrals", &
                  " is rather high."
                  write(use_unit,'(1X,A,A)') &
                  "* This is an untested accuracy modification, so do not ", &
                  "disable this warning unless you know what you're doing."
               end if

!             call aims_stop_coll('', func)
            end if

            flag_wave_threshold = .true.

         case('force_potential')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'sc') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Potential: Self-consistent Kohn-Sham potential."
               end if

               force_potential = 0
            else if (desc_str.eq.'superpos_pot') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Potential: Non-selfconsistent superposition", &
                        " of free-atom potentials."
               end if

               force_potential = 1
            else if (desc_str.eq.'superpos_rho') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Potential: Non-selfconsistent superposition", &
                        " of free-atom densities."
               end if

               force_potential = 2
            else if (desc_str.eq.'non-selfconsistent') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Potential: Non-selfconsistent evaluation requested.", &
                        " Only eigenvalues will be calculated."
               end if

               force_potential = 3
            end if

         case('force_smooth_cutoff')
            ! Parameter that enforces smooth convergence of all basis functions to zero
            ! at the cutoff potential. If the basis function, derivative or kinetic energy
            ! value on the last two non-zero points exceeds the specfied value of
            ! smooth_cutoff_limit , then the calculation will simply be stopped ...
            !
            ! See subroutine shrink_fixed_basis.

            read(inputline,*,end=88,err=99) desc_str, smooth_cutoff_limit

            if (smooth_cutoff_limit.gt.0.d0) then
               force_smooth_cutoff = .true.
            end if

         case('use_angular_division')
            read(inputline,*,end=88,err=99) desc_str, use_angular_division

         case('points_in_batch')
            read(inputline,*,end=88,err=99) desc_str, n_points_in_batch
            flag_points_in_batch = .true.

         case('batch_size_limit')
            read(inputline,*,end=88,err=99) desc_str, batch_size_hard_limit
            flag_batch_size_limit = .true.

            write (info_str,'(2X,A,I6)') &
                  "Setting batch size limit to ", &
                  batch_size_hard_limit
            call localorb_info(info_str,use_unit,'(A)')

         case('grid_partitioning_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            flag_grid_partitioning_method = .true.

            if (desc_str.eq.'octree') then

               grid_partitioning_method = 1
               write (info_str,'(2X,A)') &
                     "Using octree method for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'parallel_octree') then

               grid_partitioning_method = 1
               if (.not.use_mpi) then
                  write (info_str,'(2X,A,A)') &
                        "Parallel octree method is available only in ", &
                        "parallel runs. Using serial octree."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  parallel_grid_partitioning = .true.
                  write (info_str,'(2X,A,A)') &
                        "Using parallel octree method ", &
                        "for grid partitioning."
                  call localorb_info(info_str,use_unit,'(A)')
               end if

            else if (desc_str.eq.'parallel_maxmin') then

               grid_partitioning_method = 5
               if (.not.use_mpi) then
                  write (info_str,'(2X,A,A)') &
                        "Parallel maxmin method is available only in ", &
                        "parallel runs. Using serial maxmin."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  parallel_grid_partitioning = .true.
                  write (info_str,'(2X,A,A)') &
                        "Using modified parallel maxmin method ", &
                        "for grid partitioning."
                  call localorb_info(info_str,use_unit,'(A)')
               end if

            else if (desc_str.eq.'parallel_hash') then

               grid_partitioning_method = 6
               if (.not.use_mpi) then
                  write (info_str,'(2X,A,A)') &
                        "Parallel hash method is available only in ", &
                        "parallel runs. Using serial maxmin."
                  call localorb_info(info_str,use_unit,'(A)')
                  grid_partitioning_method = 5
               else
                  parallel_grid_partitioning = .true.
                  write (info_str,'(2X,A,A)') &
                        "Using parallel hash method ", &
                        "for grid partitioning."
                  call localorb_info(info_str,use_unit,'(A)')
               end if

            else if (desc_str.eq.'parallel_hash+maxmin') then

               grid_partitioning_method = 6
               if (.not.use_mpi) then
                  write (info_str,'(2X,A,A)') &
                        "Parallel hash method is available only in ", &
                        "parallel runs. Using serial maxmin."
                  call localorb_info(info_str,use_unit,'(A)')
                  grid_partitioning_method = 5
               else
                  parallel_grid_partitioning = .true.
                  use_hashed_maxmin = .true.
                  write (info_str,'(2X,A,A)') &
                        "Using parallel hash+maxmin method ", &
                        "for grid partitioning."
                  call localorb_info(info_str,use_unit,'(A)')
               end if

            else if (desc_str.eq.'tetgen+metis') then

               grid_partitioning_method = 2
               use_tetgen = .true.
               write (info_str,'(2X,A)') &
                     "Using TetGen + METIS for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'qhull+metis') then

               grid_partitioning_method = 2
               use_qhull = .true.
               write (info_str,'(2X,A)') &
                     "Using qhull + METIS for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'nearest+metis') then

               grid_partitioning_method = 2
               use_nearest_neighbours = .true.
               write (info_str,'(2X,A,A)') &
                     "Using nearest neighbours ", &
                     "+ METIS for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'octant') then

               grid_partitioning_method = 3
               write (info_str,'(2X,A)') &
                     "Using octant based grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'group') then

               grid_partitioning_method = 4
               write (info_str,'(2X,A)') &
                     "Using grouping method for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else if (desc_str.eq.'maxmin') then

               grid_partitioning_method = 5
               write (info_str,'(2X,A)') &
                     "Using maxmin method for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')

            else

               write (info_str,'(2X,A)') &
                     "Unknown method for grid partitioning."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)

            end if

         case('use_integer_hash')
            read(inputline,*,end=88,err=99) desc_str, use_integer_hash
            if (use_integer_hash) then
               write (info_str,'(2X,A)') &
                    "Using integers in grid point hashing."
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('grouping_factor')
            read(inputline,*,end=88,err=99) desc_str, grouping_factor
            flag_grouping_factor = .true.

         case('min_batch_size')
            read(inputline,*,end=88,err=99) desc_str, min_batch_size
            flag_min_batch_size = .true.

            write (info_str,'(2X,A,I6)') &
                  "Minimal batch size set to ", min_batch_size
            call localorb_info(info_str,use_unit,'(A)')

         case('prune_basis_once')
            read(inputline,*,end=88,err=99) desc_str, prune_basis_once

         case('force_mpi_virtual_topo')
            read(inputline,*,end=88,err=99) desc_str, force_mpi_virtual_topo

            if (force_mpi_virtual_topo) then
               if (.not.use_mpi) then
                  write (info_str,'(2X,A)') &
                        "Virtual topologies apply only to MPI runs."
                  call localorb_info(info_str,use_unit,'(A)')
                  write (info_str,'(2X,A)') &
                        "Disabling use of virtual topologies."
                  call localorb_info(info_str,use_unit,'(A)')
                  force_mpi_virtual_topo = .false.
               else
                  write (info_str,'(2X,A)') &
                        "Forcing virtual MPI SMP-cluster topology."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            end if

         case('batch_distribution_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if (desc_str.eq.'qhull+metis') then

               use_metis_batch_distribution = .true.
               flag_use_metis_batching = .true.

               if (.not.use_mpi) then
                  write (info_str,'(2X,A,A)') &
                        "qhull+METIS grid batch distribution ", &
                        "applies only to MPI runs."
                  call localorb_info(info_str,use_unit,'(A)')
                  write (info_str,'(2X,A)') &
                        "Disabling use of qhull+METIS."
                  call localorb_info(info_str,use_unit,'(A)')
                  use_metis_batch_distribution = .false.
                  flag_use_metis_batching = .false.
               else
                  write (info_str,'(2X,A)') &
                        "Using qhull+METIS to distribute grid batches."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            else
               write (info_str,'(2X,A)') &
                     "Unknown batch distribution method."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if

         case('recompute_batches_in_relaxation')
            read(inputline,*,end=88,err=99) desc_str, recompute_batches_in_relaxation

            if (.not.recompute_batches_in_relaxation) then
               write (info_str,'(2X,A,A)') &
                     "Turning off recomputation of grid batches ", &
                     "when relaxing."
               call localorb_info(info_str,use_unit,'(A)')
            else
               write (info_str,'(2X,A)') &
                     "Recomputing grid batches during relaxation."
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('force_lebedev')
   !         this is also an option to guarantee legacy compatibility;
   !         we have two variants available for Lebedev integration grids,
   !         one as originally published by Lebedev/Laikov, and the other
   !         kindly provided by Bernard Delley. Delley s is the more
   !         accurate version, should always be used for production.


            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'estd') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Using ESTD angular ", &
                        "integration grids. "
               end if

               force_lebedev = 2
            else if (desc_str.eq.'original') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Using original Lebedev-Laikov angular ", &
                        "integration grids. "
               end if

               force_lebedev = 1
            else if ((desc_str.eq.'Delley').or.(desc_str.eq.'delley')) &
                     then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A,A)') &
                        "Using Bernard Delley's corrected Lebedev-Laikov ", &
                        "angular integration grids - ", &
                        "please cite his 1996 paper."
               end if

               force_lebedev = 0
            else

               if (myid.eq.0) then
                  write (use_unit,'(1X,A,A,A)') &
                        "* force_lebedev ", desc_str," ?"
                  write (use_unit,'(1X,A,A,A)') &
                        "* Unknown type of Lebedev-Laikov ", &
                        "angular integration ", &
                        "grids requested:"
                  write (use_unit,'(1X,A,A,A)') &
                        "* force_lebedev ", desc_str," ?"
                  write (use_unit,'(1X,A)') &
                        "* Did you intend this?."
               end if

               call aims_stop_coll('', func)
            end if

            flag_force_lebedev = .true.

         case('density_restart')
            read(inputline,*,end=88,err=99) desc_str, restart_from_density
         case('keep_restart_density_frozen')
            read(inputline,*,end=88,err=99) desc_str, freeze_restart_density

         ! FO-DFT keywords
         case('fo_folders')
            read(inputline,*,end=88,err=99) desc_str, buffer, buffer2
            fo_folder1 = trim(buffer)
            fo_folder2 = trim(buffer2)

         case('fo_dft')
            read(inputline,*,end=88,err=99) desc_str, fo_dft_option
            fo_dft = .true.
            if(fo_dft_option.eq.'fragment') then
                restart_read = .true.
                restart_write = .true.
                keep_restart_info  = .true.
                fo_fragment = .true.
                restart_write_file = 'restart.frag'
                fo_info_file = 'info.frag'
                restart_read_file = restart_write_file
                inquire(FILE=restart_read_file,EXIST=restart_file_exists)

            elseif(fo_dft_option.eq.'final') then
                restart_read_file = 'restart.combined'
                inquire(FILE=restart_read_file,EXIST=restart_file_exists)
                    if(.not.restart_file_exists) then
                        write(info_str, '(2X,A)') "No previous restart file was found for FO_DFT run:"
                        call localorb_info(info_str)
                        write(info_str, '(2X,A)') "... Generating new file from restart.f1 and restart.f2!"
                        call localorb_info(info_str)
                        !call fodft_combine()
                    else
                        write(info_str, '(2X,A)') "Found previous restart file (restart.combined),"
                        call localorb_info(info_str)
                        write(info_str, '(2X,A)') "won't create a new file from fragments!"
                        call localorb_info(info_str)
                    end if
                restart_read = .true.
                keep_restart_info  = .true.
                restart_file_exists = .true.
                fo_finalise = .true.

            else
              write(info_str,'(2X,A)') "ERROR: No fo_dft option has been set! Please choose an appropriate option!"
              call localorb_info(info_str)
              call aims_stop()
            end if

        case('fo_verbosity')
            read(inputline,*,end=88,err=99) desc_str, fo_verbosity

        case('fo_flavour')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'reset') then
                occ_fo = .true.
                write(info_str,'(2X,A)') "Selected H(2n-1)@DA / H(2n+1)@D-A- FO-DFT."
                call localorb_info(info_str)
            end if

         case('fo_deltaplus')
            read(inputline,*,end=88,err=99) desc_str, fo_deltaplus
                if (fo_deltaplus) then
                    use_fo_potential = .true.
                    write(unit=myid_string, fmt='(i8.8)') myid
                    inquire(FILE="output_grid_"//myid_string, EXIST=fo_file_exists)
                    if (fo_file_exists) then
                        fo_potential_file = .true.
                    end if
                endif

         case('fo_orbitals')
            read(inputline,*,end=88,err=99) desc_str, fo_orb1, fo_orb2, fo_range1, fo_range2, fo_type

         case('restart')

            read(inputline,*,end=88,err=99) desc_str, restart_read_file
            keep_restart_info  = .true.
            restart_read       = .true.
            restart_write      = .true.

            if ((n_periodic.gt.0.or.use_hf_realspace).and..not.(flag_force_single_restartfile &
               .or.flag_restart_eigenvectors_periodic))then !SVL: restart
               if (myid < 10)then
                  write(restart_read_file,'(A,A,I1)') &
                        trim(restart_read_file), '00',myid
               else if (myid < 100)then
                  write(restart_read_file,'(A,A,I2)') &
                        trim(restart_read_file), '0',myid
               else if (myid < 1000)then
                  write(restart_read_file,'(A,I3)') &
                        trim(restart_read_file), myid
               else if (myid < 10000) then
                  write(restart_read_file,'(A,A,I4)') &
                       trim(restart_read_file), '00',myid
               else if (myid < 100000) then
                  write(restart_read_file,'(A,A,I5)') &
                       trim(restart_read_file), '0',myid
               else if (myid < 1000000) then
                  write(restart_read_file,'(A,I6)') &
                       trim(restart_read_file), myid
               end if
            end if


             restart_write_file = restart_read_file

!            MC: old version -> fixed below
!             if(n_spin == 1 .or. .not. use_scalapack)then
!                inquire(FILE=restart_read_file,EXIST=restart_file_exists)
!                write(info_str,'(2X,2A)') 'Restart file specified as : ', &
!                     restart_read_file
!             else
!
!
!                write(info_str,'(2X,2A)') 'Restart file specified as : ', &
!                     restart_read_file
!
!                spin_restart: do i_spin = 1,n_spin
!                   write(desc_str,'(A,A,I1)') trim(restart_read_file),'.',i_spin
!                   inquire(FILE=desc_str,EXIST=restart_file_exists)
!                   if(.not.restart_file_exists) exit spin_restart
!                end do spin_restart
!
!            end if

            write(info_str,'(2X,2A)') 'Restart file specified as : ', &
                     restart_read_file
            call localorb_info(info_str,use_unit,'(A)')

            inquire(FILE=restart_read_file,EXIST=restart_file_exists)
            if(.not.restart_file_exists)then
                spin_restart: do i_spin = 1,n_spin
                   write(desc_str,'(A,A,I1)') trim(restart_read_file),'.',i_spin
                   inquire(FILE=desc_str,EXIST=restart_file_exists)
                   if(.not.restart_file_exists) exit spin_restart
                end do spin_restart
            end if

            if (restart_file_exists) then
               write(info_str,'(2X,2A)') '| Reading information from ', &
                     'existing file. '
            else
               write(info_str,'(2X,2A)') '| This file does not exist, ', &
                     'will create it after the first complete scf cycle.'
               restart_read = .false.
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_read_only')
            read(inputline,*,end=88,err=99) desc_str, restart_read_file
            keep_restart_info  = .true.
            restart_read       = .true.

            if ((n_periodic.gt.0.or.use_hf_realspace).and..not.(flag_force_single_restartfile &
               .or.flag_restart_eigenvectors_periodic)) then !SVL: restart
               if (myid < 10)then
                  write( restart_read_file,'(A,A,I1)') &
                        trim( restart_read_file), '00',myid
               else if (myid < 100)then
                  write( restart_read_file,'(A,A,I2)') &
                        trim( restart_read_file), '0',myid
               else if (myid < 1000)then
                  write(restart_read_file,'(A,I3)') &
                        trim( restart_read_file), myid
               else if (myid < 10000) then
                  write(restart_read_file,'(A,A,I4)') &
                       trim(restart_read_file), '00',myid
               else if (myid < 100000) then
                  write(restart_read_file,'(A,A,I5)') &
                       trim(restart_read_file), '0',myid
               else if (myid < 1000000) then
                  write(restart_read_file,'(A,I6)') &
                       trim(restart_read_file), myid
               end if
            end if

!            MC: old version -> fixed below
!             if(n_spin == 1 .or. .not. use_scalapack)then
!                inquire(FILE=restart_read_file,EXIST=restart_file_exists)
!             else
!                spin_restart_read: do i_spin = 1,n_spin
!                   write(desc_str,'(A,A,I1)') trim(restart_read_file),'.',i_spin
!                   inquire(FILE=desc_str,EXIST=restart_file_exists)
!                   if(.not.restart_file_exists) exit spin_restart_read
!                end do spin_restart_read
!            end if

            inquire(FILE=restart_read_file,EXIST=restart_file_exists)
            if(.not.restart_file_exists)then
                spin_restart_read: do i_spin = 1,n_spin
                   write(desc_str,'(A,A,I1)') trim(restart_read_file),'.',i_spin
                   inquire(FILE=desc_str,EXIST=restart_file_exists)
                   if(.not.restart_file_exists) exit spin_restart_read
                end do spin_restart_read
            end if


            if (restart_file_exists) then
               write(info_str,'(2X,3A)')'Reading restart information ', &
                     'from file ',restart_read_file
            else
               write(info_str,'(2X,3A)') '* WARNING: Requested restart ', &
                  'file does not exist: ', restart_read_file
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,3A)') '* Ignoring input'
               restart_read       = .false.
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_write_only')
            read(inputline,*,end=88,err=99) desc_str, restart_write_file
            keep_restart_info  = .true.
            restart_write      = .true.

            if ((n_periodic.gt.0.or.use_hf_realspace).and..not.(flag_force_single_restartfile &
               .or.flag_restart_eigenvectors_periodic)) then !SVL: restart
               if (myid < 10)then
                  write( restart_write_file,'(A,A,I1)') &
                        trim( restart_write_file), '00',myid
               else if (myid < 100)then
                  write( restart_write_file,'(A,A,I2)') &
                        trim( restart_write_file), '0',myid
               else if (myid < 1000)then
                  write(restart_write_file,'(A,I3)') &
                        trim( restart_write_file), myid
               else if (myid < 10000) then
                  write(restart_write_file,'(A,A,I4)') &
                       trim(restart_write_file), '00',myid
               else if (myid < 100000) then
                  write(restart_write_file,'(A,A,I5)') &
                       trim(restart_write_file), '0',myid
               else if (myid < 1000000) then
                  write(restart_write_file,'(A,I6)') &
                       trim(restart_write_file), myid
               end if
            end if

            write(info_str,'(2X,3A)') 'Writing restart information to ', &
               'file ', restart_write_file
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_save_iterations')
            read(inputline,*,end=88,err=99) desc_str, restart_save_iterations
            use_restart_save_iterations = .true.
            keep_restart_info           = .true.
            restart_write               = .true.
            write(info_str,'(2X,2A,I4,A)') 'Saving restart information', &
                  ' every ', restart_save_iterations, ' iterations.'
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_exx_write')
            read(inputline,*,end=88,err=99) desc_str, restart_exx_write
            write(info_str,'(2X,2A,I4,A)') 'Saving restart information', &
                  ' for exact exchange calculations.'
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_exx_read')
            read(inputline,*,end=88,err=99) desc_str, restart_exx_read
            write(info_str,'(2X,2A,I4,A)') 'Reading restart information', &
                  ' for exact exchange calculations.'
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_exx_3c_int')
            read(inputline,*,end=88,err=99) desc_str, restart_exx_3c_int
            write(info_str,'(2X,2A,I4,A)') 'Reading 3-center integrals', &
                  ' for exact exchange calculations.'
            call localorb_info(info_str,use_unit,'(A)')

         case('restart_relaxations')
            read(inputline,*,end=88,err=99) desc_str, restart_relaxations
            if (restart_relaxations) then
               write(info_str,'(2X,A)') &
                     'Restarting and saving the relaxation trajectory.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                     '*** Attention: The resulting "relaxation_restart_file.FHIaims" is only useful'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                     '*** for the "bfgs_textbook" relaxation algorithm (not our default).'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                     '*** Otherwise, just use the "geometry.in.next_step" file as the new "geometry.in" file.'
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,A)') &
                     'Explicitly not restarting and saving the relaxation trajectory.'
               call localorb_info(info_str,use_unit,'(A)')
            end if

         case('hessian_to_restart_geometry')
            read(inputline,*,end=88,err=99) desc_str, hessian_to_restart_geometry
            if (hessian_to_restart_geometry) then
               write(info_str,'(2X,A)') &
                     'Hessian information is written to geometry.in.next_step for relaxation restart.'
            else
               write(info_str,'(2X,2A)') &
                     'Hessian information is explicitly not written to ', &
                     'geometry.in.next_step for relaxation restart.'
            end if
            call localorb_info(info_str,use_unit,'(A)')

            flag_hessian_to_restart_geometry = .true.

         case('distributed_hessian')
            read(inputline,*,end=88,err=99) desc_str, use_distributed_hessian

            if (myid == 0) then
               if (use_distributed_hessian) then
                  write (info_str,'(2X,A)') &
                     "Requested to distribute Hessian matrix across MPI tasks."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  write (info_str,'(2X,A)') &
                     "Full copy of Hessian matrix on each MPI task requested explicitly."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            end if

            flag_distributed_hessian = .true.

         case('write_restart_geometry')
            read(inputline,*,end=88,err=99) desc_str, write_restart_geometry
            if (write_restart_geometry) then
               write(info_str,'(2X,A)') &
                     'A file "geometry.in.next_step" will be written after each relaxation step.'
            else
               write(info_str,'(2X,2A)') &
                     'A file "geometry.in.next_step" will NOT be written after any relaxation step.'
            end if
            call localorb_info(info_str,use_unit,'(A)')

            flag_write_restart_geometry = .true.

         case('squeeze_memory')
            read(inputline,*,end=88,err=99) desc_str, l_dummy
            write(info_str,'(A)') '  *** Warning: squeeze_memory is not needed any more, '// &
                                  'you should remove it from your input file!'
            call localorb_info(info_str,use_unit,'(A)')

         case('use_local_index')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str == 'if_scalapack') then
               flag_use_local_index_if_scalapack = .true.
            else
               read(inputline,*,end=88,err=99) desc_str, use_local_index
               if (use_local_index) then
                  write (info_str,'(2X, A, E11.4)') &
                  "Storing only matrices local to each thread."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  write (info_str,'(2X, A, E11.4)') &
                  "use_local_index requested false explicitly."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               flag_use_local_index = .true.
            end if

         case('load_balancing')
            ! WPH, 2018-01-22:  load_balancing now automatically sets use_local_index if the
            ! user has not already set it.  This was done because users were frequently setting
            ! load_balancing only, in which case nothing will happen, and they'll conclude that
            ! the keyword does nothing.

            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str == 'if_scalapack') then
               flag_load_balancing_if_scalapack = .true.
               ! If use_local_index has not yet been specified, automatically turn on use_local_index if_scalapack
               if (.not.flag_use_local_index.and..not.flag_use_local_index_if_scalapack) then
                  write (info_str,'(2X, A, E11.4)') "Since 'load_balancing if_scalapack' is present in control.in, &
                                                    &automatically turning on 'use_local_index if_scalapack'."
                  call localorb_info(info_str,use_unit,'(A)')
                  flag_use_local_index_if_scalapack = .true.
               end if
            else
               read(inputline,*,end=88,err=99) desc_str, use_load_balancing
               if (use_load_balancing) then
                  ! If use_local_index has not yet been specified, automatically turn on use_local_index .true.
                  if (.not.flag_use_local_index.and..not.flag_use_local_index_if_scalapack) then
                     write (info_str,'(2X, A, E11.4)') "Since 'load_balancing .true.' is present in control.in, &
                                                       &automatically turning on 'use_local_index .true.'."
                     call localorb_info(info_str,use_unit,'(A)')
                     use_local_index = .true.
                  end if
                  write (info_str,'(2X, A)') &
                     "Using load balancing for integrations."
                  call localorb_info(info_str,use_unit,'(A)')
               else
                  write (info_str,'(2X, A)') &
                     "load_balancing requested false explicitly."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               flag_load_balancing = .true.
            end if

         case('heterostructure')
            call aims_stop_coll("The keyword 'heterostructure' is no longer supported", func )

         case('use_alltoall')
            read(inputline,*,end=88,err=99) desc_str, use_alltoall

            if (use_alltoall) then
               write (info_str,'(2X, A)') &
                  "Usage of alltoall communication explicitly requested."
               call localorb_info(info_str,use_unit,'(A)')
            else
               write (info_str,'(2X, A)') &
                  "Usage of alltoall communication explicitly disabled."
               call localorb_info(info_str,use_unit,'(A)')
            end if

            flag_use_alltoall = .true.


            case('use_mpi_in_place')
               read(inputline,*,end=88,err=99) desc_str, use_mpi_in_place

            if (use_mpi_in_place) then
               write (info_str,'(2X, A)') &
                  "Usage of MPI_IN_PLACE communication explicitly requested."
               call localorb_info(info_str,use_unit,'(A)')
            else
               write (info_str,'(2X, A)') &
                  "Usage of MPI_IN_PLACE communication explicitly disabled."
               call localorb_info(info_str,use_unit,'(A)')
            end if


         case('use_embedding_pp')
         case('pp_charge')
         case('local_component')
         case('esp_constraint')
!

         case('output')

            read(inputline,*,end=88,err=99) desc_str, desc_str



            if (desc_str.eq.'basis') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Numerical basis data will be written out."
               end if

               out_basis = .true.

            else if (desc_str .eq. 'realspace_esp') then

               if (myid .eq.0) then
                  write (use_unit, '(2X,A)') &
                     "Realspace electrostatic potential will be written out."
               end if
               read(inputline,*,end=88,err=99) desc_str, desc_str, step_x, step_y, step_z

               flag_out_elec_real = .true.
               i_cube = i_cube+1
               cube_type(i_cube)='hartree_potential'
               cube_type_needs_densmat(i_cube)=.false.
               cube_type_needs_eigenvectors(i_cube)=.false.
               !cube_filename(i_cube) = 'electrostatic_potential_realspace'
               call read_cube_parameters ( i_cube )
               out_cube = .true.
            else if (desc_str .eq. 'on-site_esp') then
               if (myid .eq. 0) then
                  write (use_unit, '(2X,A)') &
                     "atom avearge electrostatic potential will be written out."
               endif
               flag_out_locpot_atom = .true.
            else if (desc_str.eq.'onsite_integrands') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Onsite integrands f(r) for <phi_i|H_atom|phi_i> = int dr f(r) will be written out."
               end if

               out_onsite = .true.


            else if (desc_str.eq.'ovlp_spectrum') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "EV spectrum of the overlap matrix ", &
                        "will be written out."
               end if

               out_ovlp_spectrum = .true.

            else if (desc_str.eq.'json_log' .or. desc_str.eq.'scf_timings_file') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "JSON-formatted output will be written to aims.log ."
                  write (use_unit,'(2X,A)') &
                        "- Only task 0 will output JSON text in this fashion."
               end if
               out_aims_json_log = .true.

            else if (desc_str.eq.'eigenvectors') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Kohn-Sham eigenvectors will be written out."
                  if ( (n_periodic.ne.0) .and. (n_plot_band.le.0) ) then
                    write (use_unit,'(1X,A)') &
                        "* Warning. This option will have no effect since you did not request any band structure output."
                  end if
               end if

               out_eigenvec = .true.

            else if (desc_str.eq.'soc_eigenvectors') then
               if (.not.allocated(write_soc_eigenvectors)) then
                 allocate( write_soc_eigenvectors(n_write_soc_eigenvectors) )
                 i_write_soc_eigenvectors = 0
               end if

               i_write_soc_eigenvectors = i_write_soc_eigenvectors + 1

               read(inputline,*,end=88,err=99) desc_str, desc_str, write_soc_eigenvectors(i_write_soc_eigenvectors)

            else if (desc_str.eq.'overlap_matrix') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "The overlap matrix will be written out for the k-points requested by 'output band'."
                  if (n_periodic.eq.0) then
                    write (use_unit,'(1X,A)') &
                        "* Warning. This option will have no effect in a non-periodic calculation."
                  else if (n_plot_band.le.0) then
                    write (use_unit,'(1X,A)') &
                        "* Warning. This option will have no effect since you did not request any band structure output."
                  end if

               end if

               out_overlap = .true.

            else if (desc_str.eq.'hamiltonian_matrix') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "The hamiltonian matrix will be written out for the k-points requested by 'output band'."
                  if (n_periodic.eq.0) then
                    write (use_unit,'(1X,A)') &
                        "* Warning. This option will have no effect in a non-periodic calculation."
                  else if (n_plot_band.le.0) then
                    write (use_unit,'(1X,A)') &
                        "* Warning. This option will have no effect since you did not request any band structure output."
                  end if

               end if

               out_hamiltonian = .true.

            else if ( desc_str == 'elpa_timings' .or. desc_str == 'elpa_2013_timings') then

               call localorb_info( &
               "ELPA (2013 version) timings will be written out.",&
               6,'(2X,A)' )

               elpa_print_times = .true.


            else if ( desc_str == 'k_point_list' ) then

               call localorb_info( "List of k-points will be written out.",use_unit,'(2X,A)' )

               out_k_point_list = .true.

            else if (desc_str.eq.'friction_matrices') then
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)') &
                        "* Friction matrices,First order H1 and S1, will be written to file"
               endif
               output_first_order_matrices = .true.

            else if (desc_str.eq.'friction_eigenvectors') then
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)') &
                        "* Friction eigenvectors will be printed in JMOL format"
               endif
               output_friction_eigenvectors = .true.

            else if (desc_str.eq.'matrices') then

               if (myid.eq.0) then
                  write (use_unit,'(1X,A)') &
                        "* The output for Hamilton and ovelap matrices has changed."
                  write (use_unit,'(1X,A)') &
                        "* Modern output compatible to many third-party applications"
                  write (use_unit,'(1X,A)') &
                        "* can be obtained with 'output h_s_matrices'. This results"
                  write (use_unit,'(1X,A)') &
                        "* in a format where each line of ouput has data for one entry:"
                  write (use_unit,'(1X,A)') &
                        "* <row> <column> <value>."
                  write (use_unit,'(1X,A)') &
                        "* Legacy output can be accessed with 'output matrices_2005'."
               end if

               call aims_stop_coll("'output matrices' is no longer supported", func )

            else if (desc_str.eq.'h_s_matrices') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Hamilton and overlap matrices ", &
                        "will be written out in a_ij-format."
               end if

               out_matrices = .true.

            else if (desc_str.eq.'ks_coulomb_integral') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, coulomb_integral_format
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "KS Coulomb integrals ", &
                        "will be written out in ijkl-format."
               end if
               call localorb_info( "    The chosen format is '"  &
                                        // trim(coulomb_integral_format) // "'."  )
               if (coulomb_integral_format /= 'full' .and. coulomb_integral_format /= 'qtmp' ) then
                 call aims_stop_coll(&
                    "Syntax error: Format has to be either 'full' or 'qtmp' (for quantum package use).", func )
               end if
               output_ks_coulomb_integral = .true.
            else if (desc_str.eq.'matrices_2005') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Hamilton and overlap matrices ", &
                        "will be written out in legacy format."
               end if

               out_matrices = .true.
               out_matrices_format_2005 = .true.

!           >>> AB: feb 2012
!           additional keyword is added to output overlap matrix and KS-eigenvectors
!           in a format compatible with "aitranss" module (Karlsruhe transport code)
            else if (desc_str.eq.'aitranss') then

               if (myid.eq.0) then
                  write (use_unit,'(2x,a,/,2x,a)') &
                        'Overlap matrix and KS-eigenvectors will be written out', &
                        ' in a format compatible with "aitranss" module '
               end if
               out_aitranss = .true.

!           <<< done with insert: AB, jan 2012

            else if ( desc_str == 'matrices_parallel' ) then

               read(inputline,*,end=88,err=99) desc_str, desc_str, out_mat_par

               call localorb_info(                                        &
                       "  Output option 'matrices_parallel' is set to '"  &
                          // trim(out_mat_par) // "'."                     )

               read(inputline,*,end=50,err=99) desc_str, desc_str, &
                                                 out_mat_par, out_mat_par_format

               call localorb_info( "    The chosen format is '"  &
                                        // trim(out_mat_par_format) // "'."  )

               if ( out_mat_par_format /= 'asc' .and. out_mat_par_format /= 'bin' ) then
                 call aims_stop_coll(                                            &
                    "Syntax error: Format has to be either 'asc' or 'bin'.", func )
               end if

               50 continue


            else if (desc_str.eq.'nuclear_potential_matrix') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Nuclear potential matrix ", &
                        "will be written out."
               end if

               out_nuc_matrix = .true.

            else if (desc_str.eq.'t_plus_v_matrix') then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "The sum of kinetic energy and nuclear potential matrix ", &
                        "will be written out. This sum is also called H(core) sometimes."
               end if

               out_t_plus_v_matrix = .true.
            else if (desc_str.eq.'nuc_nuc_repulsion') then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "The nuclear nuclear repulsion energy will be written &
                         &out in aims output. (needed by Quantum Package)"
               end if
               out_nuc_nuc_repulsion = .true.

            else if (desc_str.eq.'bare_coulomb_integrals') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Bare moleculare Coulomb Integrals ", &
                        "will be written out."
               end if

               use_bare_ci = .true.


            else if (desc_str.eq.'dipole') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Dipole moment will be written out."
               end if

               out_dipole_moment = .true.

               ! Physical error trap - must let user know
               if (n_periodic.ne.0) then
                  if (myid.eq.0) then
                     write (use_unit,'(1X,A)') &
                           "* Attention. You requested to output a dipole moment, but this is a periodic system!"
                     write (use_unit,'(1X,A)') &
                           "* The dipole moment of a single unit cell of a periodic system is strictly undefined "
                     write (use_unit,'(1X,A)') &
                           "* (you can change it for the same structure, just by shifting the unit cell around), "
                     write (use_unit,'(1X,A)') &
                           "* and so this output makes no sense - at least not in the way that we would calculate "
                     write (use_unit,'(1X,A)') &
                           "* the dipole moment right now. We know that there are good ways around this (Berry phase), "
                     write (use_unit,'(1X,A)') &
                           "* but this is presently not implemented - sorry! "
                     write (use_unit,'(1X,A)') &
                           "* We stop the code here to alert you to this fact. If you require the dipole moment"
                     write (use_unit,'(1X,A)') &
                           "* for periodic systems, please let us know!"
                  end if
                  call aims_stop_coll('', func)
               end if
            else if (desc_str.eq.'quadrupole') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "quadrupole moment will be written out."
               end if

               out_quadrupole_moment = .true.
            else if (desc_str.eq.'grids') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Radial integration grids will be written out."
               end if

               out_grids = .true.
            else if (desc_str.eq.'v_eff') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Effective potential will be written out."
               end if

               out_v_eff = .true.

            else if (desc_str.eq.'v_eff_new') then
!
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Effective potential will be written out in new format."
               end if

               out_v_eff_new = .true.

            else if (desc_str.eq.'dielec_func') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Dielectric function used by (L/M)PBE solver will be written out."
               end if

               out_dielec_func = .true.

            else if (desc_str.eq.'ion_dens') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "(L/M)PBE ionic charge density will be written out.&
                        &WARNING! This cube output is experimental. No guarantee for results, always recheck&
			&with delta v on FHI-aims grid. Output can take long time, try first smaller cube."
               end if

               out_ion_dens = .true.

            else if (desc_str.eq.'delta_v') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Delta v potential will be written out.&
			&WARNING! This cube output is experimental. No guarantee for results, always recheck&
			&with delta v on FHI-aims grid. Output can take long time, try first smaller cube."
               end if

               out_delta_v = .true.

            else if (desc_str.eq.'ascii') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Some output will be written in ascii."
               end if

               out_ascii = .true.

            else if (desc_str.eq.'out_nconsistent') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, out_nconsistent

               if (myid.eq.0) then
                  if (out_nconsistent .eq. -1) then
                     write (use_unit,'(2X,A)') &
                       "No output will be written after consistent step."
                  else
                     write (use_unit,'(2X,A)') &
                          "Some output will be written every ",out_nconsistent,"th consistent step"
                  endif
               end if

            else if (desc_str.eq.'out_niteration') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, out_niteration

               if (myid.eq.0) then
                  if (out_niteration .eq. -1) then
                     write (use_unit,'(2X,A)') &
                       "No output will be written after iteration step."
                  else
                     write (use_unit,'(2X,A)') &
                          "Some output will be written every ",out_niteration,"th iteration step"
                  endif
               end if

            else if (desc_str.eq.'out_ninitialize') then
               read(inputline,*,end=88,err=99) desc_str, desc_str, out_ninitialize

               if (myid.eq.0) then
                  if (out_niteration .eq. -1) then
                     write (use_unit,'(2X,A)') &
                       "No output will be written after initialization step."
                  else
                     write (use_unit,'(2X,A)') &
                          "Some output will be written every ",out_ninitialize,"th initialization step"
                  endif
               end if

            else if (desc_str.eq.'v_hartree') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Partitioned Hartree potentials will be ", &
                        "written out."
               end if

               out_v_hartree = .true.

            else if (desc_str.eq.'effective_potential') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "The effective potential will be written out."
                  write (use_unit,'(2X,A,A)') &
                        "Attention: Output happens only on FHI-aims internal grid,"
                  write (use_unit,'(2X,A,A)') &
                        "not on an even-spaced grid. Mostly intended for debug purposes."

               end if

               flag_output_effective_potential = .true.


            else if (desc_str.eq.'rho_multipole') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                        "Multipole components of partitioned density ", &
                        "will be written out."
               end if

               out_rho_multipole = .true.

            else if (desc_str.eq.'k_eigenvalue')then

               read(inputline,*,end=88,err=99) desc_str, desc_str, out_k_points_eigenvalues
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,I8)') &
                        "Number of k-sets of eigenvalues printed out:", &
                        out_k_points_eigenvalues

               end if

            else if (desc_str.eq.'postscf_eigenvalues') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                    'Found postscf_eigenvalues request to write the Kohn-Sham eigenvalues'
                  write(use_unit,'(2X,A)') &
                     '| across the Brillouin zone on a potentially dense k-space grid.'
                  write(use_unit,'(2X,A)') &
                      '| Currently, "output postscf_eigenvalues" is only supported when there are'
                  write(use_unit,'(2X,A)') &
                     '| more k-points in the s.c.f. "k_grid" than there are CPUs (serial eigenvalue solution).'
               end if

               postscf_eigenvalues = .true.

            else if (desc_str.eq.'band')then

               out_band = .true.
               call read_plot_band(inputline)
            else if (desc_str.eq.'band_mulliken') then
               out_band_mulliken = .true.
!               out_band = .true.
               write(info_str,'(1X,A)') "Mulliken analysis of all k points in band will be performed"
               call localorb_info(info_str)
               call read_plot_band(inputline)
            !TZ
            else if (desc_str .eq. 'dielectric') then
               flag_out_dielectric = .true.
               flag_dielectric_test = .true.
               call read_plot_dielectric(inputline)

         !   else if (desc_str .eq. 'absorption') then
         !      flag_out_absorption = .true.
         !      call read_plot_absorption(inputline)
            !TZ
            else if (desc_str.eq.'DFPT_phonon_band')then

               call read_plot_DFPT_phonon_band(inputline)

            else if (desc_str.eq.'band_during_scf')then

               out_band_during_scf = .true.
               call read_plot_band_during_scf(inputline)

            else if (desc_str.eq.'rrs_pbc_band')then

               !rrs_pbc_out_band = .true.
               call read_rrs_pbc_plot_band(inputline)

            ! CC: Z2 by CC / Adapted from Cmera
            else if (desc_str.eq.'Z2_invariant')then
               out_z2_invariant = .true.
               read(inputline,*) desc_str, desc_str, Z2_n_planes, Z2_n_plot, Z2_n_k_points
               write(info_str,'(2X,A)') '* Calculation of Z2 invariant requested:'
               call localorb_info(info_str, use_unit, '(A)')
               write(info_str,'(2X,A,I2)') '| Computing Wilson loop defined via index: ', Z2_n_planes
               call localorb_info(info_str, use_unit, '(A)')
               write(info_str,'(2X,A,I5,A)') '| Calculating Wannier centers at ', Z2_n_plot, ' k-points'
               call localorb_info(info_str, use_unit, '(A)')
               if ( Z2_n_k_points .lt. 0 ) then
                 write(info_str,'(2X,A)') "|   using the default number of k-points for each Wannier center."
                 call localorb_info(info_str, use_unit, '(A)')
                 write(info_str,'(2X,A)') "*** WARNING: This default might not be quantitatively safe."
                 call localorb_info(info_str, use_unit, '(A)')
               else
                 write(info_str,'(2X,A,I5,A)') "|   using ", Z2_n_k_points, " k-points for each Wannier center."
                 call localorb_info(info_str, use_unit, '(A)')
               end if

            else if (desc_str.eq.'density') then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Full charge density will be written out."
               end if

               out_density = .true.

            else if (desc_str.eq.'sorted_eigenvalue_list') then
                   out_sorted_eigenvalue_list = .true.

            else if (desc_str.eq.'cube') then
   !           Cube style output required; read a number of lines
   !           with the relevant data next ...

               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               if (desc_str.eq.'total_density') then
                  i_cube = i_cube + 1

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " Total charge density. "
                  end if

                  cube_type(i_cube)='total_density'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0

               else if (desc_str.eq.'spin_density') then
                  i_cube = i_cube + 1

                  if (spin_treatment .eq. 0) then
                     write (use_unit,'(2X,A)') "Spin density not available for unpolarized calculations."
                     call aims_stop_coll('', func)
                  endif

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " spin density. "
                  end if
                  cube_type(i_cube)='spin_density'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0

               else if (desc_str.eq.'delta_density') then
                  i_cube = i_cube + 1

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " Delta charge density. "
                  end if

                  cube_type(i_cube)='delta_density'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0

               else if (desc_str.eq.'ion_dens') then
                  i_cube = i_cube + 1

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " Ionic charge density from (L/M)PBE calculation. "
                  end if

                  cube_type(i_cube)='ion_dens'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0

               else if (desc_str.eq.'delta_v') then
                  i_cube = i_cube + 1

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " Delta Hartree potential. "
                  end if

                  cube_type(i_cube)='delta_v'
                  cube_type_needs_densmat(i_cube) = .false.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0

               else if (desc_str.eq.'dielec_func') then
                  i_cube = i_cube + 1

                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                           "'cube'-style plottable output required:", &
                           " Empirical Dielectric Function. "
                  end if

                  cube_type(i_cube)='dielec_func'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  cube_index(i_cube) = 0
! eigenstate_soc, 4 channels in 1 command: spin up/down, real/imag part
!               else if (desc_str.eq.'eigenstate_soc') then
!! not support 'homo/lumo' yet, only numbers
!                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, temp_cube_index
!                  do i_cube_soc = 1, 4 ! 4 channels
!                    i_cube = i_cube + 1
!                    cube_type(i_cube) = 'eigenstate'
!                    cube_type_needs_densmat(i_cube) = .false.
!                    cube_type_needs_eigenvectors(i_cube) = .true.
!                    cube_index(i_cube) = temp_cube_index
!!                    call read_cube_parameters ( i_cube )
!                  end do
!

! eigenstate_non_soc
               else if (desc_str.eq.'eigenstate_non_soc') then
                  i_cube = i_cube + 1
                  cube_type(i_cube)= desc_str
                  cube_type_needs_densmat(i_cube) = .false.
                  cube_type_needs_eigenvectors(i_cube) = .true.
                  cube_type_needs_soc_eigenvectors(i_cube) = .false.
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                          cube_index(i_cube)
!
               else if (desc_str.eq.'eigenstate' .or. &
                        desc_str.eq.'eigenstate_imag' .or. & !Deactived but kept for debug purposes
                        desc_str.eq.'eigenstate_density') then
                  i_cube = i_cube + 1
                  cube_type(i_cube)= desc_str
                  cube_type_needs_densmat(i_cube) = .false.
                  cube_type_needs_eigenvectors(i_cube) = .true.
                  cube_type_needs_soc_eigenvectors(i_cube) = .true. ! default to be true, needs to be used together with out_cub_soc in output_cube_p2.f90
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, desc_str
                  if ((desc_str(1:4).eq.'lumo').or.(desc_str(1:4).eq.'LUMO')) then
                     if (len_trim(desc_str).eq.4) then
                        cube_index(i_cube) = 1
                     else
                        desc_str = desc_str(5:len_trim(desc_str))
                        read(desc_str,*,end=88,err=99) cube_index(i_cube)
                        cube_index(i_cube) = cube_index(i_cube)+1
                     end if
                     cube_fermi_ref(i_cube) = .true.
                  elseif ((desc_str(1:4).eq.'homo').or.(desc_str(1:4).eq.'HOMO')) then
                     if (len_trim(desc_str).eq.4) then
                        cube_index(i_cube) = 0
                     else
                        desc_str = desc_str(5:len_trim(desc_str))
                        read(desc_str,*,end=88,err=99) cube_index(i_cube)
                        cube_index(i_cube) = cube_index(i_cube)
                     end if
                     cube_fermi_ref(i_cube) = .true.
                  else
                     read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                          cube_index(i_cube)
                  endif
                  if ((cube_index(i_cube).le.0).and.(.not.cube_fermi_ref(i_cube))) then

                     if (myid.eq.0) then
                        write(use_unit,'(1X,A,A,I4,A)') &
                           "* Illegal eigenstate number specified for ", &
                              "cube output # ", i_cube,"."
                        write(use_unit,'(1X,A,A,A)') &
                              "* There are no eigenstates with ", &
                              "numbers less ", &
                              "or equal zero - please correct control.in."
                     end if

                     call aims_stop_coll('', func)
                  end if
      !             Second sanity test for max. number of states has
      !             to be done in read_geo(), as n_states is only
      !             known after reading geometry.in.

                  if (.not.cube_fermi_ref(i_cube)) then
                     write(desc_str,'(I8)') cube_index(i_cube)
                  else
                     if (cube_index(i_cube) .eq. 0) then
                        write(desc_str,'(A4)') 'HOMO'
                     else if (cube_index(i_cube) .eq. 1) then
                        write(desc_str,'(A4)') 'LUMO'
                     else if (cube_index(i_cube) .gt. 1) then
                        write(desc_str,'(A5,I3)') 'LUMO+',cube_index(i_cube)
                     else
                        write(desc_str,'(A5,I3)') 'HOMO-',-cube_index(i_cube)
                     end if
                  end if
                  write (info_str,'(2X,3A)') &
                       "'cube'-style plottable output required:", &
                       " Wavefunction for eigenstate ", &
                       desc_str
                  call localorb_info(info_str)

               else if (desc_str.eq.'stm') then
                  i_cube = i_cube + 1
                  i_cube_stm = i_cube_stm + 1
                  cube_index(i_cube) = i_cube_stm
                  cube_type(i_cube)='stm'
                  cube_type_needs_densmat(i_cube) = .true.
                  cube_type_needs_eigenvectors(i_cube) = .false.
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                        cube_stm(2,i_cube)
                  cube_stm(1,i_cube)=0.0
                  write(info_str,'(2X,A,2(1X,F11.4))') &
                        "| STM width in eV: ", &
                     cube_stm(2,i_cube)
                  call localorb_info(info_str)
                  cube_stm(2,i_cube)=cube_stm(2,i_cube)/hartree
                  cube_stm(3,i_cube)=1.
          !Sanity check: system must be periodic
                  if (n_periodic.eq.0) then
                      call aims_stop_coll('STM not supported for non-periodic systems',func)
                  endif
               elseif (desc_str.eq.'long_range_potential') then
                  !write(info_str,*) ' *** WARNING: Cube files for potentials only implemented for'
                  !call localorb_info(info_str)
                  !write(info_str,*) ' test purposes. Do not use for production runs unless you '
                  !call localorb_info(info_str)
                  !write(info_str,*) ' exactly know what you are doing. (Esp.: do not trust the '
                  !call localorb_info(info_str)
                  !write(info_str,*) ' results in regions where there is still electron density.) '
                  !call localorb_info(info_str)
                  !i_cube = i_cube+1
                  !cube_type(i_cube)='long_range_potential'
                  !cube_type_needs_densmat(i_cube)=.false.
                  !Sanity check
                  !if (n_periodic.eq.0) then
                  !    call aims_stop_coll('Cube output for potential not supported for non-periodic systems',func)
                  !endif
                  write(info_str,*) 'Error: The long range potential keyword is now deprecated and no longer recommended'
                  call localorb_info(info_str)
                  write(info_str,*) 'Please use output cube hartree_potential instead'
                  call aims_stop_coll(info_str)
               !elseif (desc_str.eq.'potential') then
               !     i_cube = i_cube+1
               !     cube_type(i_cube)='potential'
               !     cube_type_needs_densmat(i_cube)=.false.
               !     cube_type_needs_eigenvectors(i_cube)=.false.
                  write(info_str,*) 'Error: The potential keyword is now deprecated and no longer recommended'
                  call localorb_info(info_str)
                  write(info_str,*) 'Please use output cube hartree_potential or xc_potential instead'
                  call aims_stop_coll(info_str)
               elseif (desc_str.eq.'xc_potential') then
                    i_cube = i_cube+1
                    cube_type(i_cube)='xc_potential'
                    cube_type_needs_densmat(i_cube)=.false.
                    cube_type_needs_eigenvectors(i_cube)=.false.
                    if (flag_xc.ne.6)then
                      call aims_stop_coll('Cube output for xc potential only supported for PBE',func)
                    endif
               elseif (desc_str.eq.'hartree_potential') then
                    i_cube = i_cube+1
                    cube_type(i_cube)='hartree_potential'
                    cube_type_needs_densmat(i_cube)=.false.
                    cube_type_needs_eigenvectors(i_cube)=.false.
               elseif (desc_str.eq.'elf') then
                i_cube = i_cube+1
                cube_type(i_cube)='elf'
                cube_type_needs_densmat(i_cube)=.false.
                cube_type_needs_eigenvectors(i_cube)=.false.
                cube_type_needs_soc_eigenvectors(i_cube) = .false.


               else
                  write(info_str,'(1X,A,A,A)') &
                  & "* Unknown cube style output required: ", &
                  & trim(desc_str), ". Please correct file control.in."
                  call aims_stop_coll(info_str, func)
               end if
!           At this point, we expect the necessary cube specification
!           (4 lines).
               call read_cube_parameters ( i_cube )
               out_cube = .true.

!Sanity check for cube files - TODO: Move to plot.f90
              ! Check for load_balancing first, since it will automatically set use_local_index in some cases
             ! if (use_load_balancing .and. cube_type_needs_densmat(i_cube)) then
             !    call localorb_info("")
             !    write(info_str, "(2X,A,' ',A,' ',A)") &
             !    & 'Cube type', trim(cube_type(i_cube)), &
             !    & 'not implemented for load_balancing or use_local_index keywords.  Exiting.'
             !    call aims_stop_coll(info_str, func)
             ! endif
             ! if (use_local_index .and. cube_type_needs_densmat(i_cube)) then
             !    call localorb_info("")
             !    write(info_str, "(2X,A,' ',A,' ',A)") &
             !    & 'Cube type', trim(cube_type(i_cube)), &
             !    & 'not implemented for the use_local_index keyword.  Exiting.'
             !    call aims_stop_coll(info_str, func)
             ! endif
             ! if (calculate_perturbative_soc .and. cube_type_needs_densmat(i_cube)) then
             !    call localorb_info("")
             !    write(info_str, "(2X,A,' ',A,' ',A)") &
             !    & 'Cube type', trim(cube_type(i_cube)), &
             !    & 'not implemented for spin-orbit coupling.  Exiting.'
             !    call aims_stop_coll(info_str, func)
             ! endif

!End sanity check for cube files

            else if (desc_str.eq.'vacuum_potential')then
               out_vacuum_potential =.true.

               read(inputline,*,end=88,err=99) desc_str, desc_str, out_vacuum_potential_z1,  out_vacuum_potential_z2, &
                     out_vacuum_potential_z_grid, out_vacuum_potential_x_grid, out_vacuum_potential_y_grid

               out_vacuum_potential_z1 = out_vacuum_potential_z1/bohr
               out_vacuum_potential_z2 = out_vacuum_potential_z2/bohr


            else if (desc_str.eq.'embedding_potential_z')then
               out_embedding_potential =.true.

               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                     out_embedding_potential_z1, &
                     out_embedding_potential_z2, &
                     out_embedding_potential_z_grid

               out_embedding_potential_z1 = out_embedding_potential_z1/bohr
               out_embedding_potential_z2 = out_embedding_potential_z2/bohr




            else if ( (desc_str.eq.'mulliken') .or. &
                     (desc_str.eq.'Mulliken')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Full Mulliken charge analysis will be written ", &
                  "to file 'Mulliken.out'."
               end if

               out_mulliken = .true.
               flag_run_mulliken = .true.

            else if ( (desc_str.eq.'eigenvec_ovlp') .or. &
                     (desc_str.eq.'Eigenvec_ovlp')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Output eigenvector and overlap matrix ", &
                  "to file 'eigenvec.out' and 'ovlp_mat.out'."
               end if

               use_out_eigenvec_ovlp = .true.
               use_full_spectrum = .true. ! Change to calculate_all_eigenstates?
               use_symmetry_reduced_k_grid = .false.
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "WARNING: use_out_eigenvec_ovlp will force: ", &
                  "use_full_spectrum = .true. and use_symmetry_reduced_k_grid = .false."
               end if

               if (use_scalapack) then
                  collect_eigenvectors = .true.
                  if (myid.eq.0) then
                     write (use_unit,'(2X,A,A)') &
                     "use_out_eigenvec_ovlp and use_scalapack will force: ", &
                     "collect_eigenvectors = .true."
                  end if
               end if

            else if ( (desc_str.eq.'lowdin') .or. &
                     (desc_str.eq.'Lowdin')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A,A)') &
                  "Full Loewdin charge analysis will be written ", &
                  "to file 'Loewdin.out'."
               end if

               out_lowdin = .true.
               flag_run_lowdin = .true.



            else if ( (desc_str.eq.'hirshfeld') .or. &
                     (desc_str.eq.'Hirshfeld')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Hirshfeld charge analysis will be performed."
               end if

               out_hirshfeld = .true.

            else if ( (desc_str.eq.'hirshfeld_always') .or. &
                     (desc_str.eq.'Hirshfeld_always')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Hirshfeld charge analysis will be written for every geometry step,"
                  write (use_unit,'(2X,A)') &
                  "| even for output_level 'MD_light'."
               end if

               out_hirshfeld_always = .true.
               out_hirshfeld = .true.

            else if ( (desc_str.eq.'hirshfeld-I') .or. &
                     (desc_str.eq.'Hirshfeld-I')      ) then

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Hirshfeld-I charge analysis will be performed."
               end if

               out_hirshfeld_iterative = .true.


            else if (desc_str.eq.'zero_multipoles') then
               read(inputline,*,end=88,err=99) desc_str,desc_str,out_zero_multipoles

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Atomic charges (l=0 multipole moments) of the Hartree potential will be written out."
               end if

            else if (desc_str.eq.'hartree_multipoles') then

               out_hartree_multipoles = .true.

               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "All multipole moments of the Hartree potential will be written out."
               end if

            else if (desc_str.eq.'dos') then
               out_dos =.true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,dos_low_energy, dos_high_energy, &
                     dos_n_en_points, dos_alpha
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)')'Calculating density of states with'
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        dos_low_energy, dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & broadening :', &
                        dos_n_en_points, dos_alpha
               end if
            else if (desc_str.eq.'dos_tetrahedron') then
               out_dos_tetrahedron =.true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,dos_low_energy, dos_high_energy, &
                     dos_n_en_points
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)')'Calculating density of states with'
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        dos_low_energy, dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & tetrahedron method :', &
                        dos_n_en_points
               end if

            else if (desc_str.eq.'atom_proj_dos') then
               flag_run_mulliken = .true.
   !             flag_run_lowdin = .true.
               out_atom_dos =.true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,atom_dos_low_energy, atom_dos_high_energy, &
                     atom_dos_n_en_points, atom_dos_alpha
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)')'Calculating atom-projected ', &
                                    'density of states with'
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        atom_dos_low_energy, atom_dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & broadening :', &
                        atom_dos_n_en_points, atom_dos_alpha
               end if

            else if (desc_str.eq.'atom_proj_dos_tetrahedron') then
               flag_run_mulliken = .true.
   !             flag_run_lowdin = .true.
               out_atom_dos_tetrahedron =.true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,atom_dos_low_energy, atom_dos_high_energy, &
                     atom_dos_n_en_points
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)')'Calculating atom-projected ', &
                                    'density of states with'
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        atom_dos_low_energy, atom_dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & tetrahedron :', &
                        atom_dos_n_en_points
               end if

            else if (desc_str.eq.'species_proj_dos') then
               flag_run_mulliken = .true.
   !             flag_run_lowdin = .true.
               out_l_proj_dos = .true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,l_proj_dos_low_energy, l_proj_dos_high_energy, &
                     l_proj_dos_n_en_points, l_proj_dos_alpha
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') 'Calculating angular-momentum', &
                        ' projected DOS with '
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        l_proj_dos_low_energy, l_proj_dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & broadening :', &
                        l_proj_dos_n_en_points, l_proj_dos_alpha
               end if
            else if (desc_str.eq.'species_proj_dos_tetrahedron') then
               flag_run_mulliken = .true.
   !             flag_run_lowdin = .true.
               out_l_proj_dos_tetrahedron = .true.
               read(inputline,*,end=88,err=99)desc_str,desc_str,l_proj_dos_low_energy, l_proj_dos_high_energy, &
                     l_proj_dos_n_en_points
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') 'Calculating angular-momentum', &
                        ' projected DOS with '
                  write(use_unit,'(2X,A,2F10.6)')'| energy interval :', &
                        l_proj_dos_low_energy, l_proj_dos_high_energy
                  write(use_unit,'(2X,A,I6,F10.6)')'| points & tetrahedron :', &
                        l_proj_dos_n_en_points
               end if
            else if (desc_str.eq.'dipole_polarisability') then
               out_polarisability = .true.
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') 'MP2/RPA dipole polarisability will', &
                        ' will be written out. '
               endif
            else if (desc_str.eq.'self_energy') then
               out_self_energy = .true.
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') 'Quasiparticle self-energy will be', &
                        'be written out. '
                  call system('mkdir self_energy')
               endif
            else if (desc_str.eq.'gw_regular_kgrid') then
               out_gw_regular_kgrid = .true.
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') 'Calculate the periodic GW self-energy on', &
                        ' regular k grid. '
               endif
            else if (desc_str.eq.'hessian' .or. desc_str.eq.'Hessian') then
               out_hessian = .true.
               call localorb_info('  Hessian will be written out.')
            else if (desc_str.eq.'memory' ) then
               call aims_stop('The "output memory" keyword is no longer supported.', func)
            else if (desc_str.eq.'memory_tracking' ) then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                        "Verbose memory tracking output has been activated."
               end if
               aims_mem_debug = .true.
            else if (desc_str.eq.'fermisurface') then
                write(info_str,*) 'Output of fermisurface requested.'
                call localorb_info(info_str)
                write(info_str,*) 'output will be given as .bxsf file'
                call localorb_info(info_str)
                write(info_str,*) 'See XCrysden Manual for details'
                call localorb_info(info_str)
                write(info_str,*) 'Note: You might want to set symmetry_reduced_k_grid to false.'
                call localorb_info(info_str)
                write(info_str,*) 'Note: This is an untested feature'
                call localorb_info(info_str)
                out_fermisurface = .true.

            else if ( (desc_str.eq.'ESP') .or. &
                     (desc_str.eq.'esp')      ) then
               ! Check for load_balancing first, since it will automatically set use_local_index in some cases
             !  if (use_load_balancing) then
             !    write(info_str, "(2X,A)") &
             !    & 'Esp charges not implemented for load_balancing or use_local_index keywords'
             !    call aims_stop_coll(info_str, func)
             !  endif
             !  if (use_local_index) then
             !    write(info_str, "(2X,A)") &
             !    & 'Esp charges not implemented for the use_local_index keyword'
             !    call aims_stop_coll(info_str, func)
             !  endif
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                  "Fitting ESP-charges to Potential. "
               end if
               i_esp = i_esp + 1
               out_esp = .true.
               !if (use_scalapack) then
               !  collect_eigenvectors = .true.
               !endif
               call read_esp_parameters(i_esp)
            else if ( desc_str.eq.'soc_eigenvalues' .or. &
                      desc_str.eq.'SOC_eigenvalues' ) then
               out_soc_eigenvalues = .true.
               call localorb_info('  Printing of SOC-perturbed eigenvalues at each SCF k-point requested.')
               call localorb_info('    If SOC is enabled, results will be output to an SOC_eigenvalues.dat file.')
               call localorb_info('    If SOC is not enabled, this request will be ignored.')
            else if ( desc_str.eq.'batch_statistics' ) then
               out_batch_statistics = .true.
               call localorb_info( &
                    '  Statistics about the batch distribution used for real-space operations will be output to disk.')
            else if ( desc_str.eq.'matrices_elsi') then
               out_matrices_elsi = .true.
            ! CC: DGrid ouput
            else if ((desc_str.eq.'DGrid').or.(desc_str.eq.'Dgrid').or.&
                     (desc_str.eq.'dgrid').or.(desc_str.eq.'DGRID')) then
               out_dgrid = .TRUE.
               call localorb_info( &
                    '  Separate outputfile dgrid_aims.dat will be generated for later DGrid analysis.')
            else
               if (myid.eq.0) then
                  write(use_unit,*) '"', desc_str, '"', " is not a valid output", &
                        " identifier. "
               end if

               call aims_stop_coll('', func)
            end if

         case('exx_band_structure_version')
            read(inputline,*,end=88,err=99) desc_str, exx_band_structure_version

         case('xml_file')
            out_xml = .true.
            read(inputline, *, end=88, err=99) desc_str, xml_file
            if (myid == 0) then
                call xml_open_file (xml_file, 'results')
            end if

         case('output_level')
            read(inputline,*,end=88,err=99) desc_str, output_level
            if ((output_level .eq. 'full'  ) .or. &
               (output_level .eq. 'normal') .or. &
               (output_level .eq. 'MD_light')) then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') "Requested output level: ", output_level
               end if
            else
               if (myid.eq.0) then
                  write(use_unit,'(2X,2A)') "*** WARNING: output_level not known - ", output_level
                  write(use_unit,'(2X,A)' ) "*** Please check! "
               end if
               call aims_stop_coll('', func)
            end if

         case('species')
            i_species = i_species + 1

            read(inputline,*,end=88,err=99) desc_str, species_name(i_species)

            call read_species_data &
               ( i_species, flag_radial(i_species), &
               flag_angular(i_species), flag_angular_min(i_species), &
               flag_angular_acc(i_species), flag_cut_free_atom(i_species), &
               flag_prodbas_acc(i_species), flag_max_n_prodbas(i_species), &
               flag_max_l_prodbas(i_species) &
               )


         case('auxil_basis')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'opt'.or.desc_str.eq.'optimal') then
               flag_auxil_basis = PRODBAS_OPT
               call localorb_info("  Optimized auxiliary basis is used.")
            else if (desc_str.eq.'full') then
               flag_auxil_basis = PRODBAS_FULL
               call localorb_info("  Full auxiliary basis is used.")
            else if (desc_str.eq.'svd') then
               flag_auxil_basis = PRODBAS_SVD
               call localorb_info("  SVD auxiliary basis is used.")
               call localorb_info("  *** Warning: Onsite SVD (auxil_basis svd)&
                                  & is EXPERIMENTAL and inefficient.")
            else
               write(info_str,*) "The present choice of the auxiliary basis", &
               & trim(desc_str), "is not known."
               call aims_stop_coll(info_str, func)
            endif

         case('final_forces_cleaned')
            read(inputline,*,end=88,err=99) desc_str, final_forces_cleaned
            if (final_forces_cleaned) then
               write(info_str,'(2X,A)') &
                  'Cleaning forces before final output. '
            else
               write(info_str,'(2X,A)') &
               'Explicitly NOT cleaning forces before final output.'
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('output_in_original_unit_cell')
            read(inputline,*,end=88,err=99) desc_str, output_in_original_unit_cell
            if (output_in_original_unit_cell) then
               write(info_str,'(2X,A)') &
                  "Printing relaxed geometry in original unit cell."
            else
               write(info_str,'(2X,A)') &
               "Explicit request to output atoms in internal orientation."
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_petukhov_mixing')
            plus_u_petukhov_mixing_defined = .true.
            read(inputline,*,end=88,err=99) desc_str, plus_u_petukhov_mixing
            write(info_str,'(2X,A,F16.12)') &
                     'Fixing petukhov mixing factor to ',plus_u_petukhov_mixing
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_ramping_accuracy')
            plus_u_ramping_defined = .true.
            read(inputline,*,end=88,err=99) desc_str, plus_u_ramping_accuracy
            write(info_str,'(2X,A,F16.12)') &
                     'DFT+U ramping is ON: U value will be increased until target value is reached'
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_occ_mat_error')
            plus_u_matrix_error_defined = .true.
            write(info_str,'(2X,A,F16.12)') &
                     'DFT+U also calculates the idempotence error of the occupation matrix (Tr(n-n*n))'
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_use_hydros')
            plus_u_hydros_defined = .true.
            write(info_str,'(2X,A,F16.12)') &
                     'hydrogenic basis funktions also used for DFT+U reference states'
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_matrix_control')
            plus_u_occupation_matrix_control_read = .true.
            plus_u_occupation_matrix_control_write = .true.
            inquire(file='occupation_matrix_control.txt',exist=occ_mat_file_exists)
            if(occ_mat_file_exists)then
               write(info_str,'(2X,A,F16.12)') &
                     'Reading DFT+U occupation matrix, from file'
               call localorb_info(info_str,use_unit,'(A)')
            else
               plus_u_occupation_matrix_control_read = .false.
            endif

         case('plus_u_matrix_release')
            plus_u_occupation_matrix_control_read = .true.
            plus_u_occupation_matrix_control_write = .true.
            inquire(file='occupation_matrix_control.txt',exist=occ_mat_file_exists)
            if(occ_mat_file_exists)then
               write(info_str,'(2X,A,F16.12)') &
                     'Reading DFT+U occupation matrix, from file'
               call localorb_info(info_str,use_unit,'(A)')
            else ! lets use the inital density, however write out occupation matrix in the end
               plus_u_occupation_matrix_control_read = .false.
            endif
            plus_u_matrix_release_defined = .true.
            read(inputline,*,end=88,err=99) desc_str, plus_u_matrix_release
            write(info_str,'(2X,A,F16.12)') &
                  'DFT+U occupation matrix, is again calculated if total energy is converged within ', plus_u_matrix_release
            call localorb_info(info_str,use_unit,'(A)')

         case('plus_u_out_eigenvalues')
            plus_u_eigenvalues_defined = .true.
            write(info_str,'(2X,A,F16.12)') &
                     'All eigenvalues of the self-consistent DFT+U occupation matrix are printed out after each run'
            call localorb_info(info_str,use_unit,'(A)')


         case('plus_u_use_mulliken')
            if(use_forces) then
               call aims_stop_coll('No forces available for DFT+U with mulliken charges.', func)
            else
               plus_u_mulliken_charges = .true.
               write(info_str,'(2X,A,F16.12)') &
                     'DFT+U occupation matrix is now calculated using the dual-representation'
               call localorb_info(info_str,use_unit,'(A)')
            endif

         case('plus_u_use_full')
            if(use_forces) then
              call aims_stop_coll('No forces available for DFT+U with full representation.', func)
            else
              !call aims_stop_coll('WARNING: DFT+U with full representation, experimental only!', func)
              plus_u_full = .true.
              write(info_str,'(2X,A,F16.12)') &
                    'DFT+U occupation matrix is now calculated using the full-representation'
              call localorb_info(info_str,use_unit,'(A)')
            endif

         case('vdw_correction')
            write(info_str,'(2X,A,I4)') &
                     'Using empirical van der Waals correction: ', &
                     vdw_pairs
            call localorb_info(info_str,use_unit,'(A)')
            call read_vdw(vdw_pairs)

        case ('many_body_dispersion', 'many_body_dispersion_nl')
            call localorb_info( &
                '  Using libMBD for calculating MBD and TS energies', &
                use_unit &
            )
            select case (desc_str)
            case ('many_body_dispersion')
            case ('many_body_dispersion_nl')
                mbd_input%method = 'mbd-nl'
            end select
            buffer_arr(:) = ''
            read (inputline, *, iostat=i_code) desc_str, buffer_arr
            if (i_code > 0) go to 99
            do i_tag = 1, size(buffer_arr)
                tag = buffer_arr(i_tag)
                if (len(trim(tag)) == 0) exit
                idx = index(tag, '=')
                if (idx == 0) then
                    select case (str_lower(tag))
                        case default
                            go to 99
                    end select
                else
                    key = tag(1:idx-1)
                    val = str_replace(tag(idx+1:), ':', ' ')
                    select case (str_lower(key))
                        case ('ts_ene_acc')
                            read (val, *, iostat=i_code) mbd_input%ts_ene_acc
                        case ('ts_f_acc')
                            read (val, *, iostat=i_code) mbd_input%ts_f_acc
                        case ('k_grid')
                            read (val, *, iostat=i_code) mbd_input%k_grid
                        case ('k_shift')
                            read (val, *, iostat=i_code) mbd_input%k_grid_shift
                        case ('beta')
                            read (val, *, iostat=i_code) mbd_beta
                        case ('freq_grid')
                            read (val, *, iostat=i_code) mbd_input%n_omega_grid
                        case ('zero_negative')
                            read (val, *, iostat=i_code) mbd_input%zero_negative_eigvals
                        case ('grid_out')
                            read (val, *, iostat=i_code) grid_out
                        case ('debug')
                            read (val, *, iostat=i_code) mbd_input%debug
                        case ('do_rpa')
                            read (val, *, iostat=i_code) mbd_input%do_rpa
                        case ('rpa_orders')
                            read (val, *, iostat=i_code) mbd_input%rpa_orders
                        case ('rpa_rescale_eigs')
                            read (val, *, iostat=i_code) mbd_input%rpa_rescale_eigs
                        case ('parallel_mode')
                            read (val, *, iostat=i_code) mbd_input%parallel_mode
                        case default
                            go to 99
                    end select
                    if (i_code /= 0) go to 99
                end if
            end do
        case ('many_body_dispersion_alt')
             if(use_mbd_old) then
             write(info_str,'(2X,2A)') &
            'Self-consistent screening and van der Waals energy (TS, TS+SCS,MBD) module'
             call localorb_info(info_str,use_unit,'(A)')
             endif
        case ('many_body_dispersion_dev')
            write (info_str, '(2X,A)') &
                'Using the experimental MBD implementation'
            call localorb_info(info_str, use_unit, '(A)')
            mbd_dev_flags(:) = ''
            read (inputline, *, iostat=i_code) desc_str, mbd_dev_flags
            if (i_code > 0) go to 99
            i_code = mbd_dev_parse()
            if (i_code /= 0) go to 99
        case ('many_body_dispersion_pre2019')
            write (info_str, '(2X,A)') &
                'Using the default MBD implementation with analytical forces'
            call localorb_info(info_str, use_unit, '(A)')
            do i_flag = 1, 10
                if (mbd_std_flags(i_flag) == '') exit
            end do
            read (inputline, *, iostat=i_code) desc_str, mbd_std_flags(i_flag:)
            if (i_code > 0) go to 99
            ! the mbd flags are actualy parsed after this master do loop
            ! because we need to know the k-grid and the xc functional
         case ('mbd_cfdm_dip_cutoff') !=100
           if (n_periodic.eq.0) then
               if (myid.eq.0) then
               write (use_unit,'(1X,A)')&
                     "* WARNING - found 'mbd_cfdm_dip_cutoff' keyword."
               write (use_unit,'(1X,A)')&
                     "* this keyword is only useful in case of periodic system"
               write (use_unit,'(1X,A)')&
                     "* Ignoring 'mbd_cfdm_dip_cutoff' in control.in ."
               end if
            else
              read(inputline,*,end=88,err=99) desc_str,mbd_cfdm_dip_cutoff
              if(mbd_cfdm_dip_cutoff.lt.10.0) then
                 write(info_str,'(2X,2A,f6.2,A)') &
                 '***Warning: Radial distance used to sum dipole field in ', &
                 'MBD calculation is too small',mbd_cfdm_dip_cutoff, " Angstrom"
                 call localorb_info(info_str,use_unit,'(A)')
                 write(info_str,'(2X,2A)') &
                 'MBD energy may not converge please check convergece of MBD ', &
                 'energy by increasing value of parameter mbd_cfdm_dip_cutoff '
                  call localorb_info(info_str,use_unit,'(A)')
              else
                  write(info_str,'(2X,A,f6.2,A)') &
                  '| Using radial distance of ',mbd_cfdm_dip_cutoff, &
                  ' Angstrom to sum dipole field in MBD calculation '
                  call localorb_info(info_str,use_unit,'(A)')
             endif
           endif
         case ('mbd_scs_dip_cutoff')  !=120
           if (n_periodic.eq.0) then
               if (myid.eq.0) then
               write (use_unit,'(1X,A)')&
                     "* WARNING - found 'mbd_scs_dip_cutoff' keyword."
               write (use_unit,'(1X,A)')&
                     "* this keyword is only useful in case of periodic system"
               write (use_unit,'(1X,A)')&
                     "* Ignoring 'mbd_scs_dip_cutoff' in control.in ."
               end if
            else
              read(inputline,*,end=88,err=99) desc_str,mbd_scs_dip_cutoff
              if(mbd_scs_dip_cutoff.lt.10.0) then
                write(info_str,'(2X,2A,f6.2,A)') &
                   '***Warning: Radial distance used to sum dipole field in ', &
                   'SCS calculation is too small',mbd_scs_dip_cutoff," Angstrom"
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,3A)') &
                   'C6 coefficients may not converge please check ', &
                   'convergence of C6 coefficients by increasing value of ', &
                   'parameter mbd_scs_dip_cutoff'
                call localorb_info(info_str,use_unit,'(A)')
              else
                write(info_str,'(2X,A,f6.2,A)') &
                   '| Using radial distance of ', mbd_scs_dip_cutoff, &
                   ' Angstrom to sum dipole field in SCS calculation '
                call localorb_info(info_str,use_unit,'(A)')
              endif
           endif
         case ('mbd_supercell_cutoff')!=25
           if (n_periodic.eq.0) then
               if (myid.eq.0) then
               write (use_unit,'(1X,A)')&
                     "* WARNING - found 'mbd_supercell_cutoff' keyword."
               write (use_unit,'(1X,A)')&
                     "* this keyword is only useful in case of periodic system"
               write (use_unit,'(1X,A)')&
                     "* Ignoring 'mbd_supercell_cutoff' in control.in ."
               end if
            else
              read(inputline,*,end=88,err=99) desc_str,mbd_supercell_cutoff
              if(mbd_supercell_cutoff.lt.5.0) then
                write(info_str,'(2X,2A,f6.2,A)') &
                   '***Warning: Radial cutoff distance used to creat super ', &
                   'cell in MBD calculation is too small', &
                   mbd_supercell_cutoff," Angstrom"
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(2X,3A)') &
                   'please check convergece of MBD energy using larger ', &
                   'super cell by increasing value of parameter ', &
                   'mbd_supercell_cutoff'
                call localorb_info(info_str,use_unit,'(A)')
              else
                write(info_str,'(2X,A,f6.2,A)') &
               '| Using cutoff distance of ',mbd_supercell_cutoff,' Angstrom to creat super cell in MBD calculation'
                call localorb_info(info_str,use_unit,'(A)')
              endif
            endif
         case('mbd_scs_vacuum_axis')
           if (n_periodic.eq.0) then
               if (myid.eq.0) then
               write (use_unit,'(1X,A)')&
                     "* WARNING - found 'mbd_cfdm_dip_cutoff' keyword."
               write (use_unit,'(1X,A)')&
                     "* this keyword is only useful in case of periodic system"
               write (use_unit,'(1X,A)')&
                     "* Ignoring 'mbd_cfdm_dip_cutoff' in control.in ."
               end if
            else
               read(inputline,*,end=88,err=99) desc_str,mbd_scs_vacuum_axis(1),mbd_scs_vacuum_axis(2),mbd_scs_vacuum_axis(3)

               if(mbd_scs_vacuum_axis(1)) then
                  write(info_str,'(2X,A)') '| Using vaccum in X-direction in MBD/SCS calculation'
                  call localorb_info(info_str,use_unit,'(A)')

               if(mbd_scs_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in X-direction is less than mbd_scs_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in X-direction larger than mbd_scs_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif
               if(mbd_cfdm_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in X-direction is less than mbd_cfdm_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in X-direction larger than mbd_cfdm_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               endif

               if (mbd_scs_vacuum_axis(2))then
                  write(info_str,'(2X,A)') '| Using vaccum in Y-direction in MBD/SCS calculation'
                  call localorb_info(info_str,use_unit,'(A)')

               if(mbd_scs_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in Y-direction is less than mbd_scs_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in Y-direction larger than mbd_scs_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif
               if(mbd_cfdm_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in Y-direction is less than mbd_cfdm_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in Y-direction larger than mbd_cfdm_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif

               endif

               if(mbd_scs_vacuum_axis(3))then
                   write(info_str,'(2X,A)') '| Using vaccum in Z-direction in MBD/SCS calculation'
                   call localorb_info(info_str,use_unit,'(A)')

               if(mbd_scs_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in Z-direction is less than mbd_scs_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in Z-direction larger than mbd_scs_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif
               if(mbd_cfdm_dip_cutoff/bohr.gt.SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) then
                  write(info_str,'(2X,A)') '*** Warning detected the length vaccum level in Z-direction is less than mbd_cfdm_dip_cutoff'
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A)') '*** Please set vaccum in Z-direction larger than mbd_cfdm_dip_cutoff to avoid artificial interactions with periodic images'
                  call localorb_info(info_str,use_unit,'(A)')
               endif
               endif

               if(mbd_scs_vacuum_axis(1).AND.mbd_scs_vacuum_axis(2).AND.mbd_scs_vacuum_axis(3)) then
                   write(info_str,'(A)') '***Warning vaccum is specified in all three direction in MBD/SCS calculation.'
                   call localorb_info(info_str,use_unit,'(A)')
                   write(info_str,'(A)') 'Please check the keyword mbd_scs_vacuum_axis for consistency. Aborting current run.....'
                   call localorb_info(info_str,use_unit,'(A)')
                   call aims_stop_coll ()
               endif
           endif
         case('mbd_eigensolver')
            flag_mbd_eigensolver = .true.
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq.'lapack') then
              write(info_str,'(2X,2A)') '| Using LAPACK for solving ', &
                'MBD eigenvalue problem.'
              call localorb_info(info_str,use_unit,'(A)')
              use_scalapack_mbd = .false.
              use_elpa_mbd      = .false.
            else if (desc_str.eq.'scalapack') then
              write(info_str,'(2X,2A)') '| Using scaLAPACK for solving ', &
                'MBD eigenvalue problem.'
              call localorb_info(info_str,use_unit,'(A)')
              use_scalapack_mbd = .true.
              use_elpa_mbd      = .false.
            else if (desc_str.eq.'elpa') then
              write(info_str,'(2X,2A)') '| Using ELPA for solving ', &
                'MBD eigenvalue problem.'
              call localorb_info(info_str,use_unit,'(A)')
              use_scalapack_mbd = .true.
              use_elpa_mbd      = .true.
            else
               write(info_str,'(2X,3A)') '| Using ', desc_str, ' for solving ', &
                'MBD eigenvalue problem is not a valid treatment.'
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop()
            end if

         case('vdw_correction_hirshfeld')
            if (use_vdw_correction_hirshfeld) then
               write(info_str,'(2X,2A)') &
                       'Using non-empirical van der Waals correction', &
                       ' (Hirshfeld partitioning): '
            else
               write(info_str,'(2X,2A)') &
                       'Non-empirical van der Waals correction', &
                       ' (Hirshfeld partitioning) explicitly disabled.'
            end if
         case('vdw_correction_hirshfeld_sc')
            if (use_vdw_correction_hirshfeld_sc) then
               write(info_str,'(2X,2A)') &
                       'Using non-empirical van der Waals correction', &
                       ' (Hirshfeld partitioning): '
            else
               write(info_str,'(2X,2A)') &
                       'Non-empirical van der Waals correction', &
                       ' (Hirshfeld partitioning) explicitly disabled.'
            end if
            call localorb_info(info_str,use_unit,'(A)')

         case('vdw_convergence_threshold')
             ! Convergence threshold (in eV) for supercell sum
             ! in case the Tkatchenko-Scheffler vdW correction is used
             ! in periodic systems
             read(inputline,*,end=88,err=99) desc_str, vdw_convergence_threshold

             if (vdw_convergence_threshold.le.0.d0) then
                write(info_str,'(1X,A,E13.6,A)') &
                       '*** Error. The specified "vdw_convergence_threshold" value ', &
                       vdw_convergence_threshold, ' eV could lead to an infinite loop.'
                call localorb_info(info_str,use_unit,'(A)')
                write(info_str,'(1X,A)') &
                       '*** Please choose a positive value.'
                call localorb_info(info_str,use_unit,'(A)')
                call aims_stop('*** Invalid vdw_convergence_threshold', func)
             else

                write(info_str,'(2X,A,E13.6,A)') &
                       'Found vdw_convergence_threshold value: ', &
                       vdw_convergence_threshold, ' eV.'
                call localorb_info(info_str,use_unit,'(A)')

             end if

             ! Conversion to Hartree units
             vdw_convergence_threshold = vdw_convergence_threshold/hartree

             flag_vdw_convergence_threshold = .true.

         case('vdw_cells')
            read(inputline,*,end=88,err=99) desc_str, vdw_n1, vdw_n2, vdw_n3

            write (info_str,'(2X,A,I3,I3,I3,A)') &
                  "VdW energy will be summed over: (", vdw_n1, vdw_n2, &
                  vdw_n3, ") cells"
            call localorb_info(info_str,use_unit,'(A)')

         case('vdw_pair_ignore')
            i_vdw_ignore = i_vdw_ignore + 1
            read(inputline,*,end=88,err=99) desc_str, vdw_pairs_ignore(1,i_vdw_ignore), vdw_pairs_ignore(2,i_vdw_ignore)
            write(info_str,'(2X,4A)') "Ignoring vdw interaction between the species pair ",&
                  trim(vdw_pairs_ignore(1,i_vdw_ignore)), &
                  ' and ',trim(vdw_pairs_ignore(2,i_vdw_ignore))

         case('vdw_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            if(desc_str.eq.'octree')then
               read (inputline,*,end=88,err=99) desc_str, desc_str, epsfinal
               treeflag=int(3)

               write(info_str,'(2X,A,F16.12)')"Using octrees for vdw. Epsilon:", epsfinal
               call localorb_info(info_str,use_unit,'(A)')
               if(epsfinal.eq.0d0)then
                  write(use_unit,*)"****Must specify nonzero epsilon for octrees!"
                  stop
               endif

            elseif(desc_str.eq.'mixed')then
               read (inputline,*,end=88,err=99) desc_str, desc_str, epsfinal
               treeflag=4
               write(info_str,'(2X,A,F16.12)')"Using octrees and multipoles for vdw. Epsilon:", epsfinal
               call localorb_info(info_str,use_unit,'(A)')

            elseif(desc_str.eq.'multipoles')then
               treeflag=int(1)
               write(info_str,'(2X,A)')"Using multipoles for vdw."
               call localorb_info(info_str,use_unit,'(A)')
            endif

         case('nlcorr_nrad')
            read(inputline,*,end=88,err=99) desc_str, nrad_nlcorr
            write (info_str,'(2X, A, I5)') &
                 "radial grid points for nonlocal correlation integrals:", &
                 nrad_nlcorr
            call localorb_info(info_str,use_unit,'(A)')
            flag_nrad_nlcorr = .true.

         case('nlcorr_i_leb')
            read(inputline,*,end=88,err=99) desc_str, i_leb_nlcorr
            write (info_str,'(2X, A, I5)') &
                 "lebedev grid index for nonlocal correlation integrals:", &
                 i_leb_nlcorr
            call localorb_info(info_str,use_unit,'(A)')
            flag_i_leb_nlcorr = .true.

         case('qmmm')
            write(info_str,'(2X,A,I4)') &
                     'QM/MM embedding on.'
            call localorb_info(info_str,use_unit,'(A)')


         case('add_embedding_grids')
            read(inputline,*,end=88,err=99) desc_str, add_embedding_grids

            if(.not.(add_embedding_grids) ) then

               write(info_str,'(2X,A)') &
                         'The global integration grid is NOT extended around the pseudocores!'
               call localorb_info(info_str,use_unit,'(A)')

               write(info_str,'(2X,A)') &
                         'This is dangerous, only use this flag when you are absolutely sure, what you are doing!'
               call localorb_info(info_str,use_unit,'(A)')

            else
               write(info_str,'(2X,A)') &
                     'The global integration grid is extended around the pseudocores.'
               call localorb_info(info_str,use_unit,'(A)')
            endif


         case('hf_version')
            read(inputline,*,end=88,err=99)desc_str, buffer
            select case (buffer)
            case('0', 'density_matrix')
               hf_version = HF_DM
               write(info_str,'(2X,A,A)') &
               & 'Using density matrix for evaluating', &
               & ' exchange matrix (new version).'
               call localorb_info(info_str)
            case('1', 'eigencoefficients', 'overlap')
               hf_version = HF_EIGEN
               write(info_str,'(2X,A,A)') &
               & 'Using transformed overlap matrix for evaluating', &
               & ' exchange matrix (old version).'
               call localorb_info(info_str)
            case default
               call aims_stop('*** Invalid hf_version', func)
            end select

         case('use_logsbt')
            read(inputline,*,end=88,err=99) desc_str, use_logsbt
            if (use_logsbt) then
               write(info_str,'(2X,A)') &
                     'Using 1D integrations for auxiliary 2-center integrals.'
               call localorb_info(info_str)
            else
               write(info_str,'(2X,A)') &
                     'Using 3D grid integrations for auxiliary 2-center integrals.'
               call localorb_info(info_str)
            end if
            flag_use_logsbt = .true.

         case('sbtgrid_lnr0')
            read(inputline,*,end=88,err=99) desc_str, sbtgrid_lnr0
            write(info_str,'(2X,A,F16.12)') &
            & 'Onset of logarithmic r-grid for SBT set to', sbtgrid_lnr0
            call localorb_info(info_str)
            flag_sbtgrid_lnr0 = .true.

         case('sbtgrid_lnk0')
            read(inputline,*,end=88,err=99) desc_str, sbtgrid_lnk0
            write(info_str,'(2X,A,F16.12)') &
            & 'Onset of logarithmic k-grid for SBT set to', sbtgrid_lnk0
            call localorb_info(info_str)
            flag_sbtgrid_lnk0 = .true.

         case('sbtgrid_lnrange')
            read(inputline,*,end=88,err=99) desc_str, sbtgrid_lnrange
            write(info_str,'(2X,A,F16.12)') &
            & 'Range of logarithmic r- and k-grids for SBT set to',&
            & sbtgrid_lnrange
            call localorb_info(info_str)
            flag_sbtgrid_lnrange = .true.

         case('sbtgrid_N')
            read(inputline,*,end=88,err=99) desc_str, sbtgrid_N
            write(info_str,'(2X,A,I7)') &
            & 'Number of logarithmic r- and k-grid points for SBT set to',&
            & sbtgrid_N
            call localorb_info(info_str)
            flag_sbtgrid_N = .true.

         case('cutCb_width')
            read(inputline,*,end=88,err=99) desc_str, cutCb_width
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cut-width to', cutCb_width, ' * rcut.'
            call localorb_info(info_str)
            if (flag_cutCb_width /= CUTCB_FROM_DEFAULT) then
               call aims_stop_coll('Inconsistent cutCb_width', func)
            end if
            flag_cutCb_width = CUTCB_FROM_EXPLICIT

         case('cutCb_rcut')
            read(inputline,*,end=88,err=99) desc_str, cutCb_rcut
            cutCb_rcut = cutCb_rcut / bohr
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cutting radius rcut to', &
            & cutCb_rcut*bohr, ' A.'
            call localorb_info(info_str)
            if (flag_cutCb_rcut /= CUTCB_FROM_DEFAULT) then
               call aims_stop_coll('Inconsistent cutCb_rcut', func)
            end if
            flag_cutCb_rcut = CUTCB_FROM_EXPLICIT

         case('cutCb_width_factor')
            read(inputline,*,end=88,err=99) desc_str, cutCb_width_factor
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cut-width factor to', cutCb_width_factor, '.'
            call localorb_info(info_str)
            if (flag_cutCb_width /= CUTCB_FROM_DEFAULT) then
               call aims_stop_coll('Inconsistent cutCb_width', func)
            end if
            flag_cutCb_width = CUTCB_FROM_FACTOR

         case('cutCb_rcut_factor')
            read(inputline,*,end=88,err=99) desc_str, cutCb_rcut_factor
            write(info_str,'(2X,A,F16.12,A)') &
            & 'Setting Coulomb cutting radius rcut factor to', &
            & cutCb_width_factor, '.'
            call localorb_info(info_str)
            if (flag_cutCb_rcut /= CUTCB_FROM_DEFAULT) then
               call aims_stop_coll('Inconsistent cutCb_rcut', func)
            end if
            flag_cutCb_rcut = CUTCB_FROM_FACTOR

         case('use_logsbt_for_radial_hse_integration')
            read(inputline,*,end=88,err=99) desc_str, use_logsbt_for_radial_hse_integration
            if (use_logsbt_for_radial_hse_integration) then
               call localorb_info('Using logSBT for radial HSE integration', &
               &                  use_unit, '(2X,A)')
            else
               call localorb_info('Using hse_matrix for radial HSE integration', &
               &                  use_unit, '(2X,A)')
            end if
            flag_use_logsbt_for_radial_hse_integration = .true.

         case('use_logsbt_lvltriples')
            read(inputline,*,end=88,err=99) desc_str, use_logsbt_lvltriples
            if (use_logsbt_lvltriples) then
               call localorb_info('Using logSBT for LVL triples', &
               &                  use_unit, '(2X,A)')
            else
               call localorb_info('Using grid integration for LVL triples', &
               &                  use_unit, '(2X,A)')
            end if
            ! Use .false. here only for debugging purposes.

         case('use_2d_corr')
            read(inputline,*,end=88,err=99) desc_str, use_2d_corr
            if (use_2d_corr) then
               write(info_str,'(2X,A)') 'Efficient 2D distribution for correlation code.'
               call localorb_info(info_str)
            else
               write(info_str,'(2X,A)') 'Old 1D distribution for correlation code.'
               call localorb_info(info_str)
            end if
            flag_use_2d_corr = .true.
         case('use_ovlp_swap')
               write(info_str,'(2X,A,A)') &
                  '3-center overlap matrix (ovlp_3fn) is written to', &
                  ' disk and deallocated before MP2 calculation starts. '
               call localorb_info(info_str,use_unit,'(A)')
               use_ovlp_swap = .true.

         case('phonon')
            if (.not.flag_found_phonon) then
               write(info_str,'(2X,2A)') "Found phonon prescription, ", &
                  "this only works if called from within phonon script."
               call localorb_info(info_str,use_unit,'(A)')
               flag_found_phonon = .true.
            end if

         case('unfold')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (.not. flag_found_unfold) then
               write(info_str,'(2X 2A)') "Found prescription, ", &
                  "for band structure unfolding, still under development."
               call localorb_info(info_str,use_unit,'(A)')
               flag_found_unfold = .true.
               out_overlap  = .true.
               out_hamiltonian  = .false.
               out_eigenvec = .true.
            end if

            if(desc_str.eq.'format')then

               read (inputline,*,end=88,err=99) desc_str, desc_str, unfold_format
               if (unfold_format.eq.'hdf5') then
                  write(info_str,'(2X 3A)')"| Band structure unfolding: ", &
                       "The output format is set to ", unfold_format
                  call localorb_info(info_str,use_unit,'(A)')
                  flag_hdf5_unfold  = .true.
                  flag_ASCII_unfold = .false.
               else if (unfold_format.eq.'ASCII') then
                  write(info_str,'(2X 3A)')"| Band structure unfolding: ", &
                       "The output format is set to ", unfold_format
                  call localorb_info(info_str,use_unit,'(A)')
                  flag_ASCII_unfold = .true.
                  flag_hdf5_unfold  = .false.
               end if
            end if

         case('vibrations')
            write(info_str,'(2X,A)') "Found vibrations prescription"
            call localorb_info(info_str,use_unit,'(A)')

         case('scgw_no_core')
               read(inputline,*,end=88,err=99)desc_str, first_state_included
                flag_scgw_no_core = .true.
               write(info_str,'(2X,A)') &
               'Excluding part of the core-states for self-consistent GW.'
         case('plot_self_energy')
              read(inputline,*,end=88,err=99)desc_str, selfe_matr_el
              plot_self_energy = .true.
         case('freq_grid_type')
               read(inputline,*,end=88,err=99)desc_str, desc_str
            if (desc_str.eq."GL" .or. desc_str.eq."Gauss-Legendre") then
               write(info_str,'(2X,A)') &
               'Using Gauss-Legendre frequency grid. '
               freq_grid_type=0
               call localorb_info(info_str,use_unit,'(A)')
               flag_freq_grid_type = .true.
            else if (desc_str.eq."mod-GL".or.desc_str.eq."mod-Gauss-Legendre") then
               write(info_str,'(2X,A)') &
               'Using modified Gauss-Legendre frequency grid. '
               freq_grid_type=1
               call localorb_info(info_str,use_unit,'(A)')
               flag_freq_grid_type = .true.
            else if (desc_str.eq."log".or.desc_str.eq."logarithmic") then
               write(info_str,'(2X,A)') &
               'Using logarithmic frequency grid. '
               freq_grid_type=2
               call localorb_info(info_str,use_unit,'(A)')
               flag_freq_grid_type = .true.
            else if (desc_str.eq."homo".or.desc_str.eq."homogeneous") then
               write(info_str,'(2X,A)') &
               'Using homogeneous frequency grid. '
               freq_grid_type=3
               call localorb_info(info_str,use_unit,'(A)')
               flag_freq_grid_type = .true.

            else
               write(info_str,'(2X,3A)') &
               "This type of frequency grid is not implemented. ", &
               "Please change to the Gausse-Legendre grid (seting 'GL') ", &
               "or its modified version ('mod-GL')"
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            endif
          case("n_poles")
             read(inputline,*,end=88,err=99)desc_str, n_poles
             write(info_str,'(2X,A,I6)') &
             'Number of poles used for fitting the Green function: ', &
             n_poles
             call localorb_info(info_str,use_unit,'(A)')
          case("scgw_chem_pot")
             read(inputline,*,end=88,err=99)desc_str, scgw_chem_pot
             specify_mu = .true.
          case("scgw_print_all_spectrum")
             scgw_print_all_spectrum = .true.
          case("scgw_scf_accuracy")
             read(inputline,*,end=88,err=99)desc_str, scgw_scf_accuracy
          case("scgw_it_limit")
             read(inputline,*,end=88,err=99)desc_str, scgw_it_limit
          case("scgw_mix_param")
             read(inputline,*,end=88,err=99)desc_str, scgw_mix_param
          case("dist_pole_type")
             read(inputline,*,end=88,err=99)desc_str, pole_dist_type
          case("pole_max")
             read(inputline,*,end=88,err=99)desc_str, pole_max
          case("pole_min")
             read(inputline,*,end=88,err=99)desc_str, pole_min
          case("print_self_energy")
             read(inputline,*,end=88,err=99)desc_str, state_to_print
             flag_print_self_energy = .true.
             call get_optional_parameter(inputline,param_position=3,real_parameter=print_range(1),&
                                         explicit=explicit)
             call get_optional_parameter(inputline,param_position=4,real_parameter=print_range(2),&
                                         explicit=explicit)
             if(.not.explicit) then
               print_range(1)=-10.0d0
               print_range(2)=10.0d0
             endif
             if(print_range(1) > print_range(2)) call aims_stop('check range for Sigma-Printing',func)
         case('n_max_coeff_3fn')
            read(inputline,*,end=88,err=99)desc_str, n_max_coeff_3fn
            write(info_str,'(2X,A,I6)') &
            'Number of unit cells of the array "coeff_3fn" in periodic exact-exchange matrix evaluation: ', &
            n_max_coeff_3fn
            call localorb_info(info_str,use_unit,'(A)')
            flag_max_coeff_3fn=.true.

         case('frequency_points')
            read(inputline,*,end=88,err=99)desc_str, n_full_freq
            write(info_str,'(2X,A,I6)') &
            'Number of frequency points used for the self-energy calculation: ', &
            n_full_freq
            call localorb_info(info_str,use_unit,'(A)')
            flag_frequency_points=.true.

         case("time_points")
            read(inputline,*,end=88,err=99)desc_str, n_full_time
            write(info_str,'(2X,A,I6)') &
            'Number of time points used for the self-energy calculation: ', &
            n_full_time
            call localorb_info(info_str,use_unit,'(A)')
            flag_time_points=.true.

         case('maximum_frequency')
            read(inputline,*,end=88,err=99)desc_str, omegamax
            write(info_str,'(2X,A,f14.2,A)') &
            'Maximum frequency for the self-energy: ', &
            omegamax, '  Ha.'
            call localorb_info(info_str,use_unit,'(A)')
            flag_maximum_frequency = .true.
         case('maximum_time')
               read(inputline,*,end=88,err=99)desc_str, taumax
            write(info_str,'(2X,A,f14.2,A)') &
            'Maximum time for the self-energy: ', &
            taumax
            call localorb_info(info_str,use_unit,'(A)')
            flag_maximum_time = .true.
         case('state_lower_limit')
               read(inputline,*,end=88,err=99)desc_str, n_low_state
            write(info_str,'(2X,A,I6)') &
            'Lower limit of the eigenstates for the self-energy correction : ', &
            n_low_state
            call localorb_info(info_str,use_unit,'(A)')

         case('state_upper_limit')
               read(inputline,*,end=88,err=99)desc_str, n_high_state
            write(info_str,'(2X,A,I6)') &
            'Upper limit of the eigenstates for the self-energy correction : ', &
            n_high_state
            call localorb_info(info_str,use_unit,'(A)')

         case('n_anacon_par')
               read(inputline,*,end=88,err=99)desc_str, n_max_par
            write(info_str,'(2X,A,I6)') &
            'Number of fitting parameters for analytical continuation : ', &
            n_max_par
            call localorb_info(info_str,use_unit,'(A)')
            flag_n_anacon_par = .true.

         case('anacon_type')
            read(inputline,*,end=88,err=99)desc_str, desc_str
            if ((desc_str .eq. '0').or.(desc_str.eq.'two-pole')) then
               anacon_type = 0
               write(info_str,'(2X,A)') &
                  'Using two-pole fitting for analytical continuation.  '
               call localorb_info(info_str,use_unit,'(A)')
               flag_anacon_type = .true.
            else if ((desc_str .eq. '1').or.(desc_str.eq.'pade')) then
               anacon_type = 1
               write(info_str,'(2X,A)') &
                  'Using Pade approximation for analytical continuation. '
               call localorb_info(info_str,use_unit,'(A)')
               flag_anacon_type = .true.
            else
               write(info_str,'(2X,A,A)') &
                  "The requested type of analytical continuation (anacon_type) does not exist.", &
                  "Please set anacon_type to 'two-pole' or 'pade'."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            endif
            anacon_type_defined = .true.
         case('contour_def_gw')
              do i_spin = 1, n_spin
                 if(i_spin ==  1) then
                   read(inputline,*,end=88,err=99)desc_str,&
                   gw_cd%contour_def_start(i_spin), gw_cd%contour_def_end(i_spin)
                 else
                   call get_optional_parameter(inputline,param_position=4,&
                                               int_parameter=gw_cd%contour_def_start(i_spin),&
                                               explicit=explicit)
                   if(.not.explicit) gw_cd%contour_def_start(i_spin)=gw_cd%contour_def_start(1)
                   call get_optional_parameter(inputline,param_position=5,&
                                               int_parameter=gw_cd%contour_def_end(i_spin),&
                                               explicit=explicit)
                   if(.not.explicit) gw_cd%contour_def_end(i_spin)=gw_cd%contour_def_end(1)
                 endif
                 write(info_str,'(2X,A42,I6,T56,A2,T60,I6)') &
                  'Contour deformation in GW for the levels: ', gw_cd%contour_def_start(i_spin),&
                  'to', gw_cd%contour_def_end(i_spin)
                 call localorb_info(info_str,use_unit,'(A)')
              enddo
         case('contour_spin_channel')
              read(inputline,*,end=88,err=99)desc_str, tmp
              gw_cd%spin_channel(1:n_spin) = tmp
              write(info_str,'(2X,A46,T60,I6)') &
               'Restrict contour deformation to spin channel: ', tmp
              call localorb_info(info_str,use_unit,'(A)')
         case('contour_eta')
              read(inputline,*,end=88,err=99)desc_str,gw_cd%eta
              write(info_str,'(2X,A26,F14.6)') &
               'Broadening parameter eta: ', gw_cd%eta
              call localorb_info(info_str,use_unit,'(A)')
         case('contour_zshot_offset')
              read(inputline,*,end=88,err=99)desc_str,&
              gw_cd%contour_def_offset
              write(info_str,'(2X,A16,F14.6)') &
               'Zshot parameter: ', gw_cd%contour_def_offset
              call localorb_info(info_str,use_unit,'(A)')
         case('contour_restart')
              read(inputline,*,end=88,err=99)desc_str,&
              gw_cd%restart
              write(info_str,'(2X,A)') &
               'Use restart option for GW calculation with contour deformation'
              call localorb_info(info_str,use_unit,'(A)')
         case('pre_adjust_qp_relativistic')
              flag_adjust_qp_pre = .true.
              write(info_str,'(2X,A)') &
                   'Adjusting KS eigenvalues before quasi-particle calculation.'
              call localorb_info(info_str,use_unit,'(A)')
              call read_relativistic_qp_corrections(inputline)
         case('post_adjust_qp_relativistic')
              flag_adjust_qp_post = .true.
              write(info_str,'(2X,A)') 'Adjusting quasi-particles'
              call localorb_info(info_str,use_unit,'(A)')
              call read_relativistic_qp_corrections(inputline)
         case('calc_spectral_func')
              flag_calc_spectral_func = .true.
              read(inputline,*,end=88,err=99)desc_str, &
              gw_cd%min_freq_spec, gw_cd%max_freq_spec
              call get_optional_parameter(inputline,param_position=4,&
                                          real_parameter=gw_cd%spec_resolution,&
                                          explicit=explicit)
              if(.not.explicit) gw_cd%spec_resolution = 0.001d0
              write(info_str,'(2X,A42,F10.3,3X,A2,3X,F10.3,2X A2)') &
               'Spectral function will be calculated from: ', gw_cd%min_freq_spec,&
                'to',  gw_cd%max_freq_spec, 'eV'
              call localorb_info(info_str,use_unit,'(A)')
         case('spectral_func_state')
              gw_cd%calc_single_state_spec = .true.
              read(inputline,*,end=88,err=99)desc_str, gw_cd%state_spec
              write(info_str,'(2X,A32,I6)') &
               'Calculate spectral function only for state:', gw_cd%state_spec
         case('full_cmplx_sigma')
              read(inputline,*,end=88,err=99)desc_str,gw_cd%full_cmplx_sigma
               write(info_str,'(2X,A)') &
                'Evaluate full complex self energy for CD (including complex integral term). '
              call localorb_info(info_str,use_unit,'(A)')
         case('iterations_sc_cd')
              read(inputline,*,end=88,err=99)desc_str,gw_cd%n_iter_sc
              write(info_str,'(2X,A53,I6)') &
              'Number of iterations for ev_scGW or ev_scGW0 with CD: ', gw_cd%n_iter_sc
              call localorb_info(info_str,use_unit,'(A)')
         case('nocc_sc_cd')
              do i_spin =1, n_spin
                 read(inputline,*,end=88,err=99)desc_str,gw_cd%sc_env(i_spin)%n_occ
                 write(info_str,'(2X,A78,I6,A9,I4)') &
                 'Number of occupied states explicity calculated in ev_scGW or ev_scGW0 with CD: ',&
                  gw_cd%sc_env(i_spin)%n_occ, 'for spin', i_spin
                 call localorb_info(info_str,use_unit,'(A)')
              enddo
         case('nvirt_sc_cd')
              do i_spin =1, n_spin
                 read(inputline,*,end=88,err=99)desc_str,gw_cd%sc_env(i_spin)%n_virt
                 write(info_str,'(2X,A78,I6,A9,I4)') &
                 'Number of virtual states explicity calculated in ev_scGW or ev_scGW0 with CD: ',&
                  gw_cd%sc_env(i_spin)%n_virt, 'for spin', i_spin
                 call localorb_info(info_str,use_unit,'(A)')
              enddo
         case('sc_reiterate')
              do i_spin =1, n_spin
                 read(inputline,*,end=88,err=99)desc_str,gw_cd%sc_env(i_spin)%reiterate
                 write(info_str,'(2X,A)') &
                 'Reiterate the state of interest after convergence of ev-scGW or ev-scGW0'
                 call localorb_info(info_str,use_unit,'(A)')
              enddo
         case('try_zshot')
              read(inputline,*,end=88,err=99)desc_str,gw_cd%try_zshot
              write(info_str,'(2X,A)') &
               'Try zshot if QP energies do not converge within ev_scGW or ev_scGW0 with CD '
              call localorb_info(info_str,use_unit,'(A)')
         case('gw_hedin_shift')
              read(inputline,*,end=88,err=99)desc_str
              call get_optional_parameter(inputline,param_position=2,int_parameter=state_to_shift,&
                                          explicit=explicit)
              use_hedin_shift = .true.
              if(explicit) then
                 write(info_str,'(2X,A50,I6)') &
                 "Use Hedin's poor-man's self-consistency for state:", state_to_shift
              else
                 write(info_str,'(2X,A)') &
                 "Use Hedin's poor-man's self-consistency and termine shift for each state"
              endif
              call localorb_info(info_str,use_unit,'(A)')
         case('qpe_multi_solu')
            read(inputline,*,end=88,err=99)desc_str, qpe_multi_solu
            if (qpe_multi_solu) then
                write(info_str,'(2X,A)') &
                 'Multi solutions of the quasiparticle energy will be determined, if they exist.  '
                call localorb_info(info_str,use_unit,'(A)')
            endif
         case('gw_zshot')
            read(inputline,*,end=88,err=99)desc_str, gw_zshot
            if (gw_zshot) then
                write(info_str,'(2X,A)') &
                 'Using tailor expansion of the quasiparticle equations instead of iterative scheme.'
                call localorb_info(info_str,use_unit,'(A)')
            endif
         case('force_occupation_projector')
            i_force_occ = i_force_occ + 1
            read(inputline,*,end=88,err=99) desc_str, force_occ_pr_state(i_force_occ), &
               force_occ_spin(i_force_occ), &
               forced_occ_number(i_force_occ),&
               force_occ_min_KS_state(i_force_occ),force_occ_max_KS_state(i_force_occ)
            if (myid.eq.0) then
               write (use_unit,'(2X,A,I3,1X,I8,1X,F12.8)') &
               "Occupation forced using projector ",force_occ_spin(i_force_occ),force_occ_pr_state(i_force_occ),&
               forced_occ_number(i_force_occ)
            endif

         case('force_occupation_projector_redsub')
            read(inputline,*,end=88,err=99) desc_str, force_occ_maxred_KS_state(i_force_occ),force_occ_step_KS_state(i_force_occ)
            if (myid.eq.0) then
               write (use_unit,'(2X,A,I3,1X,I3)') &
               "Mom subspace reduction parameters (maxred, step): ", &
               force_occ_maxred_KS_state(i_force_occ),force_occ_step_KS_state(i_force_occ)
            endif

         case('force_occupation_projector_auto')
            read(inputline,*,end=88,err=99) desc_str, force_occ_autored(i_force_occ),force_occ_autorep(i_force_occ)
            fop_detect_auto = 0
            if (myid.eq.0) then
               write (use_unit,'(2X,A,I3,1X,I3)') &
               "Automatic subspace reduction (start-step, required repetitions): ", &
               force_occ_autored(i_force_occ), force_occ_autorep(i_force_occ)
            endif

         case('force_occupation_smearing')
            force_occupation_smearing=.true.
            read(inputline,*,end=88,err=99) desc_str, force_occupation_smearing_width
            if (myid.eq.0) then
               write (use_unit,'(2X,A,1X,F12.8,1X,A)') &
               "Occupation forced with smearing width of ", force_occupation_smearing_width, "eV"
            end if
            force_occupation_smearing_width = force_occupation_smearing_width / hartree

         case('force_occupation_basis')
            i_force_occ = i_force_occ + 1
            read(inputline,*,end=88,err=99) desc_str,force_occ_atom(i_force_occ),&
               force_occ_spin(i_force_occ),&
               force_occ_basis_type(i_force_occ),&
               force_occ_basis_n(i_force_occ),&
               force_occ_basis_l(i_force_occ),&
               force_occ_basis_m(i_force_occ),&
               force_occ_basis_occupation(i_force_occ),&
               force_occ_max_KS_state(i_force_occ)
!              ,force_occ_basis_coeff_rank(i_force_occ)
            if (myid.eq.0) then
               write (use_unit,'(2X,A)')&
                  "Occupation forced using basis function selection"
               write (use_unit,'(4X,A,I5)')&
               "selected atom: ",force_occ_atom(i_force_occ)
               write (use_unit,'(4X,A,I3)')&
               "selected spin channel: ",force_occ_spin(i_force_occ)
               write (use_unit,'(4X,A,I5,I5,I5)')&
               "selected basis function - n l m: ",force_occ_basis_n(i_force_occ),&
                  force_occ_basis_l(i_force_occ),force_occ_basis_m(i_force_occ)
               write (use_unit,'(4X,A,I8)')&
               "maximum number KS state for forcing occupation: ",force_occ_max_KS_state(i_force_occ)
            endif

         case('apply_boys')
            i_sub_boys = i_sub_boys + 1
            if (n_spin.eq.2) then
              read(inputline,*,end=88,err=99) desc_str, boys_sub_min_KS_state(i_sub_boys,1), &
                 boys_sub_max_KS_state(i_sub_boys,1), boys_sub_min_KS_state(i_sub_boys,2), &
                 boys_sub_max_KS_state(i_sub_boys,2), boys_sub_flags(i_sub_boys)
            else
              read(inputline,*,end=88,err=99) desc_str, boys_sub_min_KS_state(i_sub_boys,1), &
                 boys_sub_max_KS_state(i_sub_boys,1), boys_sub_flags(i_sub_boys)
            end if

            boys_sub_flag = maxval(boys_sub_flags)
            ! Quickfix: Also set Emin and Emax (runtime_choices) such as that we compute the entire dipmat
            Emin = -100000000
            Emax = +100000000
            if (myid.eq.0) then
               do i_spin = 1, n_spin
                 write (use_unit,'(2X,A,I1,1X,I2,1X,I5,1X,I5)') &
                 "Boys localization applied to spin/subspace ", i_spin, i_sub_boys, boys_sub_min_KS_state(i_sub_boys,i_spin), &
                 boys_sub_max_KS_state(i_sub_boys,i_spin)
               end do
            endif

            ! Not entirely sure we need this but probably for boys_centers
            if (use_scalapack) then
              collect_eigenvectors = .true.
            endif

         case('calculate_perturbative_soc_old_cluster')
            write(info_str,'(2X,A)') "*** The keyword calculate_perturbative_soc_old_cluster is no longer "
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** supported. Please use the include_spin_orbit instead.  Exiting."
            call localorb_info(info_str)
            call aims_stop_coll('',func)

         case('include_spin_orbit_sc')
            write(info_str,'(2X,A)') "*** The keyword include_spin_orbit_sc is no longer "
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** supported. Please use the include_spin_orbit instead.  Exiting."
            call localorb_info(info_str)
            call aims_stop_coll('',func)

         case('calculate_perturbative_soc')
            write(info_str,'(2X,A)') "*** The keyword calculate_perturbative_soc is no longer "
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** supported. Please use the include_spin_orbit instead.  Exiting."
            call localorb_info(info_str)
            call aims_stop_coll('',func)

         case('include_spin_orbit')
            if (calculate_perturbative_soc) then
              write(info_str,'(2X,A)') "*** You have specified a handling of spin-orbit coupling more than once."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** Exiting."
              call localorb_info(info_str)
              call aims_stop_coll('',func)
            end if

            read(inputline,*,end=89,err=99) desc_str, desc_str
            calculate_perturbative_soc = .true.
   ! This hack is necessary to trick read_control into supporting optional flags.
   ! If read_control hasn't set calculate_perturbative_soc, then it must have hit EOF
   ! and jumped directly to 89, and thus no flag was specified.
   89 continue
            if (.not.calculate_perturbative_soc) then
              desc_str = ""
              calculate_perturbative_soc = .true.
            end if

            ! Determine which handling of SOC the user wants, and ensure
            ! special conditions are met for each
            if (desc_str .eq. "old_cluster") then
              ! Old cluster SOC code is no longer supported
              write(info_str,'(2X,A)') "*** You have specified the legacy 'include_spin_orbit old_cluster' keyword."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** The spin-orbit coupling code no longer uses quasi-degenerate perturbation"
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** theory.  Please use the 'include_spin_orbit' keyword instead.  Exiting."
              call localorb_info(info_str)
              call aims_stop_coll('',func)
            else if (desc_str .eq. "self_consistent") then
              ! Self-consistent SOC
              write(info_str,'(2X,A)') ""
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** Self-consistent perturbative spin-orbit coupling is no longer available,"
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** as it was experimental functionality that was unmaintained."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** Exiting."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') ""
              call aims_stop_coll('Self-consistent SOC no longer supported.',func)
            else if ((desc_str .eq. "") .or. (desc_str .eq.  "non_self_consistent")) then
              ! By assumption, if no special flag is set, then the non-self-consistent variant is being requested
              if (myid.eq.0) then
                write(use_unit,'(2X,A)') "Calculating non-self-consistent second-variational spin-orbit coupling after scf-cycle"
              end if
            else
              write(info_str,'(2X,A)') "** You have specified an incorrect handling of SOC.  Exiting."
              call localorb_info(info_str)
              call aims_stop_coll('',func)
            end if

            if (spin_treatment.eq.1 .and. use_symmetry_reduced_k_grid) then
              write(info_str,'(2X,A)') &
                   & "You have request a SOC run with spin-polarization."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') &
                   & "symmetry_reduced_k_grid will be turned off."
              call localorb_info(info_str)
              use_symmetry_reduced_k_grid = .false.
            end if

            ! All consistency checks (as well as the setting of save_soc_perturbed_eigenvectors) will be done in check_consistency

         case('soc_subspace')
            write(info_str,'(2X,A)') "*** You have specified the legacy 'soc_subspace' keyword."
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** The spin-orbit coupling code no longer uses quasi-degenerate perturbation"
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** theory.  Please use the 'include_spin_orbit' keyword instead.  Exiting."
            call localorb_info(info_str)
            call aims_stop_coll('',func)

         case('soc_perturb_order')
            write(info_str,'(2X,A)') "*** You have specified the legacy 'soc_perturb_order' keyword."
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** The spin-orbit coupling code no longer uses quasi-degenerate perturbation"
            call localorb_info(info_str)
            write(info_str,'(2X,A)') "*** theory.  Please use the 'include_spin_orbit' keyword instead.  Exiting."
            call localorb_info(info_str)
            call aims_stop_coll('',func)

         case('include_pw_lda_in_v_soc')
            if (n_spin.ne.1) then
              write(info_str,'(2X,A)') "*** include_pw_lda_in_v_soc does not support spin-polarized calculations."
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** For (collinearly) spin-polarized potentials, it is unclear how to handle"
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** the non-spin-diagonal terms in the SOC operators (though, to be fair, we"
              call localorb_info(info_str)
              write(info_str,'(2X,A)') "*** haven't put that much thought into it.)  Exiting."
              call localorb_info(info_str)
              call aims_stop_coll('',func)
            end if
            read(inputline,*,end=88,err=99) desc_str, include_pw_lda_in_v_soc

         case('n_core_states_omit_from_soc')
            read(inputline,*,end=88,err=99) desc_str, n_core_states_omit_from_soc

         case('min_energy_include_in_soc')
            min_energy_include_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, min_energy_include_in_soc

         case('min_energy_save_in_soc')
            min_energy_save_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, min_energy_save_in_soc

         case('n_high_states_omit_from_soc')
            read(inputline,*,end=88,err=99) desc_str, n_high_states_omit_from_soc

         case('max_energy_include_in_soc')
            max_energy_include_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, max_energy_include_in_soc

         case('max_energy_save_in_soc')
            max_energy_save_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, max_energy_save_in_soc

         case('gap_for_min_energy_in_soc')
            gap_for_min_energy_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, gap_for_min_energy_in_soc

         case('gap_for_saved_min_energy_in_soc')
            gap_for_saved_min_energy_in_soc_set = .true.
            read(inputline,*,end=88,err=99) desc_str, gap_for_saved_min_energy_in_soc

         case('atomic_solver')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if ( desc_str .eq. "sratom" ) then
              atomic_solver = ATOMIC_SOLVER_SRATOM
              use_atomic_solver_fns_as_minimal_basis = .false.
              write(info_str,'(2X,A)') "Using sratom as the radial atomic solver to generate basis functions"
              call localorb_info(info_str)
            else if (desc_str .eq. "atom_sphere" ) then
              atomic_solver = ATOMIC_SOLVER_ATOM_SPHERE
              use_atomic_solver_fns_as_minimal_basis = .true.
              call cite_reference("atom_sphere")
              write(info_str,'(2X,A)') "Using atom_sphere as the radial atomic solver to generate basis functions"
              call localorb_info(info_str)
            else
              call aims_stop_coll('The argument to atomic_solver in control.in was not recognized, exiting.',func)
            end if

         case('finite_nuclear_radius')
            read(inputline,*,end=88,err=99) desc_str, finite_nuclear_radius
            write (info_str,'(2X,A,F15.8,A)')&
                 'For functionality that supports it, using a finite nuclear radius of ', finite_nuclear_radius, ' bohr'
            call localorb_info(info_str,use_unit, '(A)')

         case('magnetic_response')
            magnetic_response = .true.
            use_full_spectrum = .true.
            if (index(inputline, '#') /= 0) &
                 & inputline = adjustl(inputline(:index(inputline,'#')-1))
            magnetic_response_pars = inputline

         case('mr_experimental')
            mr_experimental = .true.

         case('dfpt_accuracy_n1')
            read(inputline,*,end=88,err=99) desc_str, dfpt_accuracy_n1
            if (dfpt_accuracy_n1 < 0d0) then
               call localorb_info('Provided dfpt_accuracy_n1 is negative. Exiting.', format='(2x, a)')
               call aims_stop_coll('',func)
            end if

         case('dfpt_iter_limit')
            read(inputline,*,end=88,err=99) desc_str, dfpt_iter_limit
            if (dfpt_iter_limit < 1) then
               call localorb_info('dfpt_iter_limit must be at least 1. Exiting.', format='(2x, a)')
               call aims_stop_coll('',func)
            end if

         case('dfpt_linear_mix_param')
            read(inputline,*,end=88,err=99) desc_str, dfpt_linear_mix_param
            if (dfpt_linear_mix_param < 0) then
               write(info_str,'(2X,A)') ""
               call localorb_info('Provided dfpt_linear_mix_param is negative. Exiting.', format='(2x, a)')
               call aims_stop_coll('',func)
            end if

         case('dfpt_pulay_steps')
            read(inputline,*,end=88,err=99) desc_str, dfpt_pulay_steps
            use_dfpt_pulay = .true.
            if (dfpt_pulay_steps < 1) then
               write(info_str,'(2X,A)') ""
               call localorb_info('dfpt_pulay_steps must be at least 1. Exiting.', format='(2x, a)')
               call aims_stop_coll('',func)
            end if

         case('output_sxml')
            output_sxml = .true.
            inputline = inputline
            inputline = inputline(index(inputline,'m'):) ! In case of tabs
            if (index(inputline, '#') /= 0) inputline = adjustl(inputline(:index(inputline,'#')-1))
            inputline = inputline(11:)
            if (len(trim(inputline)) /= 0) &
                 & sxml_filename = trim(adjustl(inputline))

         case('mr_gauge_origin')
            read(inputline,*,end=88,err=99) desc_str, mr_gauge_origin
            mr_gauge_origin = mr_gauge_origin/bohr

         case('force_n_electrons')
            force_n_electrons = .true.
            read(inputline,*,end=88,err=99) desc_str,forced_n_electrons
            if (myid.eq.0) then
               write (use_unit,'(2X,A,F15.8)')&
               "Forcing number of electrons to ",forced_n_electrons
               write (use_unit,'(2X,A)')&
               "YOU MUST KNOW WHAT YOU ARE DOING !!!"
            endif

         case('fixed_spin_moment')
            if (n_spin.eq.2) then
               fixed_spin_moment = .true.
               read(inputline,*,end=88,err=99) desc_str,spin_moment
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,F15.8)')&
                  "Fixed spin moment calculation, N_up - N_down = ",spin_moment
               endif
            else ! only one spin channel
               fixed_spin_moment = .false.
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)')&
                  "* WARNING - found 'fixed_spin_moment' keyword, but no spin polarization requested."
                  write (use_unit,'(1X,A)')&
                  "* Ignoring 'fixed_spin_moment' in control.in ."
               endif



            end if

          !FRICTION STUFF
          case('calculate_friction')
             read(inputline,*,end=88,err=99) desc_str, buffer
             use_friction = .true.
             if (buffer.eq.'numerical_friction') then
                 numerical_friction = .true.
                 if (myid.eq.0) write (use_unit,'(2X,A)')&
                 "Calculating first order S and H matrices by finite difference displacement"
             else if (buffer .eq. 'DFPT') then
                 numerical_friction = .false.
                 if (myid.eq.0) write (use_unit,'(2X,A)')&
                 "Calculating first order S and H matrices with DFPT"
             else
               if (myid.eq.0) then
                  write (use_unit,'(1X,A)')&
                      "TYPE OF FRICTION UNKNOWN - setting use_friction to .false."
               endif
               use_friction = .false.
             endif
             if (use_scalapack .and. (.not.collect_eigenvectors) ) then
               collect_eigenvectors = .true.
               write(info_str,'(2X,A)') '* '
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') &
                 '*   Warning: Calculation of friction tensor was requested that needs'
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') &
                 '*   "collect_eigenvectors .true.", but ".false." was explicitly requested.'
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') &
                 '*   Overriding this choice and collecting eigenvectors after all.'
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') &
                 '*   However, for very large systems, this change may lead to memory problems.'
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') &
                 '*   Please inform us if this becomes an issue for you.'
               call localorb_info(info_str,use_unit, '(A)')
               write(info_str,'(2X,A)') '* '
               call localorb_info(info_str,use_unit, '(A)')
             end if
             if (myid.eq.0) write (use_unit,'(2X,A,I8)')&
             "Nonadabatic Friction tensor will be calculated!"
         case('friction_numeric_disp')
             read(inputline,*,end=88,err=99) desc_str, friction_numeric_disp
             friction_numeric_disp = friction_numeric_disp / bohr
             if (myid.eq.0) write (use_unit,'(2X,A,F10.4)')&
             "Numerical finite-difference displacment for friction set to: ",friction_numeric_disp
         case('friction_iter_limit')
             read(inputline,*,end=88,err=99) desc_str, friction_iter_limit
             if (myid.eq.0) write (use_unit,'(2X,A,I8)')&
             "Iteration Limit of DFPT calculation for Friction Matrices set to: ",friction_iter_limit
         case('friction_accuracy_etot')
             read(inputline,*,end=88,err=99) desc_str, friction_accuracy_etot
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Convergence accuracy of total energy: ", &
                 friction_accuracy_etot
         case('friction_accuracy_rho')
             read(inputline,*,end=88,err=99) desc_str, friction_accuracy_rho
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Convergence accuracy of density: ", &
                 friction_accuracy_rho
         case('friction_accuracy_eev')
             read(inputline,*,end=88,err=99) desc_str, friction_accuracy_eev
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Convergence accuracy of sum over eigenvalues: ", &
                 friction_accuracy_eev
         case('friction_accuracy_potjump')
             read(inputline,*,end=88,err=99) desc_str, friction_accuracy_potjump
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Convergence accuracy of potential: ", &
                 friction_accuracy_potjump
         case('friction_delta_type')
             read(inputline,*,end=88,err=99) desc_str, friction_delta_type
             if (myid.eq.0) write (use_unit,'(2X,A,A,A)') &
                "Calculating friction tensor with a ", &
                 friction_delta_type," window function"
         case('friction_window_size')
             read(inputline,*,end=88,err=99) desc_str, friction_window_size
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Friction window size is ", &
                 friction_window_size
             friction_window_size = friction_window_size / hartree
         case('friction_output_spectrum')
             read(inputline,*,end=88,err=99) desc_str, friction_output_spectrum
         case('friction_broadening_width')
             read(inputline,*,end=88,err=99) desc_str, friction_broadening_width
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4)') &
                "Friction broadening width is ", &
                 friction_broadening_width
             friction_broadening_width = friction_broadening_width / hartree
         case('friction_discretization_length')
             read(inputline,*,end=88,err=99) desc_str, friction_discretization_length
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4,A)') &
                "Friction spectrum is discretized with ", &
                friction_discretization_length, " eV"
             friction_discretization_length = friction_discretization_length / hartree
         case('friction_max_energy')
             read(inputline,*,end=88,err=99) desc_str, friction_max_energy
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4,A)') &
                "Friction is included from excitations up to  ", &
                friction_max_energy, " eV"
             friction_max_energy = friction_max_energy / hartree
         case('friction_temperature')
             read(inputline,*,end=88,err=99) desc_str, friction_temperature
             if (myid.eq.0) write (use_unit,'(2X,A,E11.4,A)') &
                "Friction is calculated at temperature of  ", &
                friction_temperature, " K"
         case('friction_q_integration')
             read(inputline,*,end=88,err=99) desc_str, friction_knotk
             if (friction_knotk) then
                 if (myid.eq.0) write (use_unit,'(2X,A)') &
                    "Friction is calculated without momentum conservation"
             else
                 if (myid.eq.0) write (use_unit,'(2X,A)') &
                    "Friction is calculated with momentum conservation"
             endif
         case('friction_read_matrices')
             read(inputline,*,end=88,err=99) desc_str, friction_read_matrices
             if (friction_read_matrices) then
                 if (myid.eq.0) write (use_unit,'(2X,A)') &
                    "First order matrices are read from file"
             endif
         case('friction_coupling_matrix_mode')
             read(inputline,*,end=88,err=99) desc_str, friction_coupling_matrix_mode
             if (myid.eq.0) then
                 write (use_unit,'(2X,A,I6)') &
                "Friction coupling matrix mode is ", &
                friction_coupling_matrix_mode
             endif
         !END FRICTION STUFF

          case('DFPT') !shanghui
             read(inputline,*,end=88,err=99) desc_str, desc_str
             if(desc_str.eq.'vibration') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif


               if( flag_xc.ne.3) then
               call aims_stop_coll("Now DFPT only support pz-lda, please install libxc if you want to use other functionals",func)
               endif


             else if(desc_str.eq.'vibration_reduce_memory') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT_reduce_memory calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( flag_xc.ne.3.and. flag_xc.ne.6) then
               call aims_stop_coll("Now DFPT only support pz-lda and pbe",func)
               endif
               !-----------a warning for DFPT-PBE--------------
               if( flag_xc.eq.6) then
                 if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                  "Warning: Now 'DFPT vibration_reduce_memory'(PBE) can NOT get right Hessian, only first_order_H is right. "
                 endif
               endif

             else if(desc_str.eq.'vibration_with_moving_grid_effect') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( partition_type.ne.1)  then
               write(info_str,'(2X,3A)') &
                  "Now DFPT vibration_with_moving_grid_effect only support partition_type = rho_r2, please ", &
                  "set 'partition_type  rho_r2' in control.in  before set ", &
                  "'DFPT vibration_with_moving_grid_effect' "
               call aims_stop_coll(info_str,func)
               endif

               if( flag_xc.ne.3) then
               call aims_stop_coll("Now DFPT only support pz-lda, please install libxc if you want to use other functionals",func)
               endif


             else if(desc_str.eq.'vibration_without_moving_grid_effect') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT vibration_without_moving_grid_effect caculation will begin,", &
                  "only served as comparison with vibration_with_moving_grid_effect", &
                  "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( partition_type.ne.1)  then
               write(info_str,'(2X,3A)') &
                  "Now DFPT vibration_without_moving_grid_effect only support partition_type = rho_r2, please ", &
                  "set 'partition_type  rho_r2' in control.in  before set ", &
                  "'DFPT vibration_without_moving_grid_effect' "
               call aims_stop_coll(info_str,func)
               endif

               if( flag_xc.ne.3) then
               call aims_stop_coll("Now DFPT only support pz-lda, please install libxc if you want to use other functionals",func)
               endif




             else if(desc_str.eq.'polarizability') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif
               if(use_hartree_fock) then
                 if (.not.flag_RI) then
                    call aims_stop_coll("DFPT polarizability currently only supports RI_V. Please add RI_method V to your control.in file.",func)
                 endif
               endif

             else if(desc_str.eq.'polar_reduce_memory') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT ( use_scalapack and PM_index)  calculation."
               endif

               if ( .not. flag_packed_matrix_format)then
                  packed_matrix_format = PM_index
                  if (myid.eq.0) then
                     write(use_unit,'(2X,A)') 'DFPT_polar_reduce_memory:'
                     write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
                  end if
               else if (packed_matrix_format .ne. PM_index) then
                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A,A)') &
                     'Error: Must use packed matrix format when ', &
                     'DFPT_polar_reduce_memory ', &
                     'are both requested.'
                  end if
                  call aims_stop_coll('', func)
               end if

               if(.not. use_scalapack) then 
                    call aims_stop_coll("DFPT polar_reduce_memory currently only supports scalpack",func)
               endif 

               if(use_hartree_fock) then
                    call aims_stop_coll("DFPT polar_reduce_memory currently NOT  supports Hartree-Fock",func)
               endif


             else if(desc_str.eq.'dielectric') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               use_basis_gradients = .true. !needed for computing the momentum matrix


             else if(desc_str.eq.'phonon_gamma') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( flag_xc.ne.3) then
               call aims_stop_coll("Now DFPT only support pz-lda, please install libxc if you want to use other functionals",func)
               endif

               if(packed_matrix_format .ne. PM_none) then
                call aims_stop_coll("the PM_index version is coming soon",func)
               endif


             else if(desc_str.eq.'phonon') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT phonon calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( flag_xc.ne.3) then
               call aims_stop_coll("Now DFPT only support pz-lda, please install libxc if you want to use other functionals",func)
               endif



              ! Check if ScaLapack exists in the binary
                call BLACS_Pinfo(index1,index2)
                if (index2.gt.0) then
              ! This means for both DFT(scalpack) and DFT(lapack), we will always use
              !  DFPT(scalapack) in case scalapck is compiled in.
                   use_scalapack_DFPT_phonon = .true.
                   write(info_str,'(2A)') 'Since scalapack support is enabled,', &
                                        ' use scalpack version for DFPT_phonon'
                   call localorb_info(info_str,use_unit,'(A)')
                else
                   use_scalapack_DFPT_phonon = .false.
                   write(info_str,'(3A)') '* WARNING: Requested Scalapack run', &
                   ' without having it compiled into the binary,', &
                   ' use lapack version for DFPT_phonon.'
                   call localorb_info(info_str,use_unit,'(A)')
                endif

               if (.not.use_mpi) then
               ! Check if this even an MPI run, in  serial version, we swich to lapack version.
                  write (info_str,'(1X,A,A)') &
                  '* DFPT_phonon_ScaLAPACK  makes sense only with mpi runs. '
                  call localorb_info(info_str,use_unit,'(A)')
                  use_scalapack_DFPT_phonon = .false.
               end if



             else if(desc_str.eq.'phonon_reduce_memory') then

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "DFPT phonon (reciprocal space method) calculation will begin  ", &
                     "after the DFT/HF calculation."
               endif

               if( spin_treatment.ne.0) then
               call aims_stop_coll("Now DFPT only support n_spin=1",func)
               endif

               if( flag_xc.ne.3.and. flag_xc.ne.6) then
               call aims_stop_coll("Now DFPT only support pz-lda and pbe",func)
               endif
             endif


          case('DFPT_sc_accuracy_dm')
             read(inputline,*,end=88,err=99) desc_str, DFPT_sc_accuracy_dm
             write (info_str,'(2X, A, E14.6)') &
                 "DFPT_sc_accuracy_dm =", DFPT_sc_accuracy_dm
             call localorb_info(info_str,use_unit,'(A)')

          case('DFPT_mixing')
             read(inputline,*,end=88,err=99) desc_str, DFPT_mixing
             write (info_str,'(2X, A, E14.6)') &
                 "DFPT_mixing =", DFPT_mixing
             call localorb_info(info_str,use_unit,'(A)')

          case('DFPT_restart_read_DM1')
             read(inputline,*,end=88,err=99) desc_str, desc_str
             if (desc_str.eq.".true.") then
                restart_read_DM1 = .true.
               write (info_str,'(2X, A)') &
                 "DFPT_restart_read_DM1 is set to true"
               call localorb_info(info_str,use_unit,'(A)')
             endif

          case('DFPT_restart_write_DM1')
             read(inputline,*,end=88,err=99) desc_str,desc_str
             if (desc_str.eq.".true.") then
                restart_write_DM1 = .true.
              write (info_str,'(2X, A)') &
                 "DFPT_restart_write_DM1 is set to true"
              call localorb_info(info_str,use_unit,'(A)')
             endif

          case('DFPT_width')
             read(inputline,*,end=88,err=99) desc_str, DFPT_width
             write (info_str,'(2X, A, E14.6, A)') &
                 "DFPT_width =", DFPT_width, "eV"
             call localorb_info(info_str,use_unit,'(A)')
             DFPT_width = DFPT_width/hartree
             flag_DFPT_occupation_type = .true.
         ! atomization_bsse subroutine
         case('calculate_atom_bsse')
            read(inputline,*,end=88,err=99) desc_str, calculate_atom_bsse
            if (calculate_atom_bsse) then
               write(info_str,'(2X,A)') &
                  'Atomization BSSE calculation requested.'
               call localorb_info(info_str)
            end if
            flag_calculate_atom_bsse = .true.
         ! end

         ! implicit solvent effects
         case('solvent')
            ! implemented electrostatic implicit solvent models:
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if (desc_str.eq."mpe") then
            ! === Multipole expansion method ===
               write(info_str,'(2X,A)') &
                  "Solvent effects will be calculated using the MPE method."
               call localorb_info(info_str)
               solvent_method = SOLVENT_MPE
	    else if (desc_str.eq."mpb") then
	       write(info_str,'(2X,A)') " "
               call localorb_info(info_str)

               write(info_str,'(2X,A,A,E14.6,A)') &
                  "Solvent effects will be calculated by using MPB method."
               call localorb_info(info_str)

	       write(info_str,'(2X,A)') "################################################################################################"
               call localorb_info(info_str)
           write(info_str,'(3X,A)') "MPB solvent effects have been requested. Please read the"
               call localorb_info(info_str)
           write(info_str,'(3X,A)') "following publications (and cite them) for a correct use and also check the Manual (!):"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "---------------------------------------------------------------------------------------------------"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "Ringe, S., Oberhofer, H., Hille, C., Matera, S., Reuter, K., Function-Space oriented Solution Scheme for the"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT, JCTC 2016, 12, 4052-4066."
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "DOI: 10.1021/acs.jctc.6b00435"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "Ringe, S., Oberhofer, H., Reuter, K., Transferable ionic parameters for first-principles"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "Poisson-Boltzmann solvation calculations: Neutral solutes in"
               call localorb_info(info_str)
           write(info_str,'(3X,A)') "aqueous monovalent salt solutions, JCP 2017, 13, 134103."
               call localorb_info(info_str)
           write(info_str,'(3X,A)') "DOI: 10.1063/1.4978850"
               call localorb_info(info_str)
	       write(info_str,'(3X,A)') "---------------------------------------------------------------------------------------------------"
               call localorb_info(info_str)

               solvent_method = SOLVENT_MPB
	       call read_mpb_data() !read in all necessary parameters

	       write(info_str,'(2X,A)') "################################################################################################"
            else
            ! === unknown method requested ===
               write(info_str,'(A,2(A,1X))') &
                  "Valid choices for implicit solvent effects are: ", &
                  "mpe", "mpb"
               call aims_stop_coll(info_str, func)
            endif

         case('mpe_nonelectrostatic_model')
            ! implemented nonelectrostatic models
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
            case("linear_OV")
               mpe_nonel_model = MPE_CONST % NONEL_linOV
               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                  mpe_nonel_alpha, mpe_nonel_beta
               write(info_str,'(2X,A)') &
                  "MPE: Non-electrostatic contribution will be calculated as:"
               call localorb_info(info_str)
               write(info_str,'(4X,E12.5,1X,A,A,E12.5,1X,A,A)') &
                  mpe_nonel_alpha, "eV/AA^2"," x (cavity surface) + ", &
                  mpe_nonel_beta, "eV/AA^3"," x (cavity volume)."
               call localorb_info(info_str)
               mpe_nonel_alpha = mpe_nonel_alpha * bohr**2 / hartree
               mpe_nonel_beta = mpe_nonel_beta * bohr**3 / hartree
            case default
               ! unknown model requested
               write(info_str,'(2X,3A)') &
                  "The requested nonelectrostatic model '", desc_str, &
                  "' is unknown."
               call aims_stop_coll(info_str, func)
            end select

         case('mpe_nonelectrostatic_interface')
            read(inputline,*,end=88,err=99) desc_str, mpe_nonel_alpha_if
            write(info_str,'(2X,A)') &
               "MPE: Non-electrostatic contribution (interface) will be calculated as:"
            call localorb_info(info_str)
            write(info_str,'(4X,E12.5,1X,A,A,E12.5,1X,A,A)') &
               -mpe_nonel_alpha_if, "eV/AA^2"," x (cutout interfacial surface)"
            call localorb_info(info_str)
            mpe_nonel_alpha_if = mpe_nonel_alpha_if * bohr**2 / hartree

         case('mpe_nonelectrostatic_model_bulk')
            ! implemented nonelectrostatic models
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
            case("linear_OV")
               mpe_nonel_model_bulk = MPE_CONST % NONEL_linOV
               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                  mpe_nonel_alpha_bulk, mpe_nonel_beta_bulk
               write(info_str,'(2X,A)') &
                  "MPE: Non-electrostatic contribution (bulk) will be calculated as:"
               call localorb_info(info_str)
               write(info_str,'(4X,E12.5,1X,A,A,E12.5,1X,A,A)') &
                  mpe_nonel_alpha_bulk, "eV/AA^2"," x (cavity surface) + ", &
                  mpe_nonel_beta_bulk, "eV/AA^3"," x (cavity volume)."
               call localorb_info(info_str)
               mpe_nonel_alpha_bulk = mpe_nonel_alpha_bulk * bohr**2 / hartree
               mpe_nonel_beta_bulk = mpe_nonel_beta_bulk * bohr**3 / hartree
            case default
               ! unknown model requested
               write(info_str,'(2X,3A)') &
                  "The requested nonelectrostatic model '", desc_str, &
                  "' is unknown."
               call aims_stop_coll(info_str, func)
            end select


         case('isc_cavity_type')
            ! implemented cavity generation methods
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
            case("overlapping_spheres")
               isc_cavity_type = ISC_CONST % CAVITY_OvlpSph
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str
               select case (desc_str)
               case ("radius")
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                     isc_initial_radius
                  write(info_str,'(2X,2A)') &
                     "The implicit solvent cavity is defined ", &
                     "by overlapping spheres"
                  call localorb_info(info_str)
                  write(info_str,'(4X,A,1X,E12.5,1X,A)') &
                     "of radius", isc_initial_radius, "AA."
                  call localorb_info(info_str)
                  isc_initial_radius = isc_initial_radius / bohr
               case ("rho")
                  read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                     isc_isodensity_value
                  write(info_str,'(2X,2A)') &
                     "The implicit solvent cavity is defined ", &
                     "by overlapping spheres"
                  call localorb_info(info_str)
                  write(info_str,'(4X,2A,1X,E12.5,1X,A)') &
                     "defined by atomic electron densities ", &
                     "at value of", isc_isodensity_value, "e/AA^3."
                  call localorb_info(info_str)
                  isc_isodensity_value = isc_isodensity_value * (bohr**3)
               case default
                  ! unknown method requested
                  write(info_str,'(2X,3A)') &
                     "Subkey '", desc_str, &
                     "' for the overlapping_spheres cavity is unknown."
                  call aims_stop_coll(info_str, func)
               end select

            case("rho_free")
               isc_cavity_type = ISC_CONST % CAVITY_RhoFree
               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                  isc_isodensity_value
               if (isc_isodensity_value.le.0) then
                  write(info_str,'(A)') &
                     "The iso-density value has to be positive!"
                  call aims_stop_coll(info_str, func)
               endif
               write(info_str,'(2X,2A)') &
                  "The implicit solvent cavity is defined ", &
                  "as an iso-density surface of the "
               call localorb_info(info_str)
               write(info_str,'(4X,2A,1X,E12.5,1X,A)') &
                  "superposed density of free atoms ", &
                  "with a value of", isc_isodensity_value, "e/AA^3."
               call localorb_info(info_str)
               isc_isodensity_value = isc_isodensity_value * (bohr**3)
            case("rho_multipole_static")
               isc_cavity_type = ISC_CONST % CAVITY_RhoMPStat
               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                  isc_isodensity_value
               if (isc_isodensity_value.le.0) then
                  write(info_str,'(A)') &
                     "The iso-density value has to be positive!"
                  call aims_stop_coll(info_str, func)
               endif
               write(info_str,'(2X,2A)') &
                  "The implicit solvent cavity is defined ", &
                  "as an iso-density surface of the "
               call localorb_info(info_str)
               write(info_str,'(4X,2A,1X,E12.5,1X,A)') &
                  "multipole-expanded electron density at initialization ", &
                  "with a value of", isc_isodensity_value, "e/AA^3."
               call localorb_info(info_str)
               isc_isodensity_value = isc_isodensity_value * (bohr**3)
            case("rho_multipole_dynamic")
               isc_cavity_type = ISC_CONST % CAVITY_RhoMPDyn
               read(inputline,*,end=88,err=99) desc_str, desc_str, &
                  isc_isodensity_value
               if (isc_isodensity_value.le.0) then
                  write(info_str,'(A)') &
                     "The iso-density value has to be positive!"
                  call aims_stop_coll(info_str, func)
               endif
               write(info_str,'(2X,2A)') &
                  "The implicit solvent cavity is defined ", &
                  "as an iso-density surface of the "
               call localorb_info(info_str)
               write(info_str,'(4X,2A,1X,E12.5,1X,A)') &
                  "self-consistent multipole-expanded electron density ", &
                  "with a value of", isc_isodensity_value, "e/AA^3."
               call localorb_info(info_str)
               isc_isodensity_value = isc_isodensity_value * (bohr**3)
            case default
               ! unknown method requested
               write(info_str,'(2X,3A)') &
                  "The requested cavity generation method '", desc_str, &
                  "' is unknown."
               call aims_stop_coll(info_str, func)
            end select

         case('mpe_xml_logging')
            ! get filename
            read(inputline,*,end=88,err=99) desc_str, mpe_xml_logfile
            ! get logging level
            read(inputline,*,end=88102,err=99) desc_str, desc_str, desc_str
            select case (desc_str)
            case('off')
               mpe_xml_loglevel = MPE_CONST % XML_NOLOG
            case('basic')
               mpe_xml_loglevel = MPE_CONST % XML_BASIC
            case('medium')
               mpe_xml_loglevel = MPE_CONST % XML_MEDIUM
            case('detailed')
               mpe_xml_loglevel = MPE_CONST % XML_DETAILED
            case('custom')
               read(inputline,*,end=88,err=99) desc_str, desc_str, desc_str, &
                  mpe_xml_loglevel
               mpe_xml_loglevel = max(mpe_xml_loglevel, MPE_CONST % XML_NOLOG)
            case default
               write(info_str,'(2X,3A)') &
                  "Encountered unknown xml logging level '", desc_str, &
                  "'. Stopping."
               call aims_stop_coll(info_str, func)
            end select
88102       continue
            ! report
            select case (mpe_xml_loglevel)
            case (MPE_CONST % XML_NOLOG)
               write(info_str,'(2X,A,I3,A)') &
                  "[MPE] xml logging is disabled (", mpe_xml_loglevel, ")."
               call localorb_info(info_str)
            case (MPE_CONST % XML_BASIC)
               write(info_str,'(2X,3A,I3,A)') &
                  "[MPE] Using xml logfile '", mpe_xml_logfile, &
                  "' with loglevel 'basic' (", mpe_xml_loglevel, ")."
               call localorb_info(info_str)
            case (MPE_CONST % XML_MEDIUM)
               write(info_str,'(2X,3A,I3,A)') &
                  "[MPE] Using xml logfile '", mpe_xml_logfile, &
                  "' with loglevel 'medium' (", mpe_xml_loglevel, ")."
               call localorb_info(info_str)
            case (MPE_CONST % XML_DETAILED)
               write(info_str,'(2X,3A,I3,A)') &
                  "[MPE] Using xml logfile '", mpe_xml_logfile, &
                  "' with loglevel 'detailed' (", mpe_xml_loglevel, ")."
               call localorb_info(info_str)
            case default
               write(info_str,'(2X,3A,I3,A)') &
                  "[MPE] Using xml logfile '", mpe_xml_logfile, &
                  "' with loglevel 'custom' (", mpe_xml_loglevel, ")."
               call localorb_info(info_str)
            end select

         case('isc_cavity_restart')
            read(inputline,*,end=88,err=99) desc_str, isc_cavity_restart_read_file
            isc_cavity_restart_write_file = isc_cavity_restart_read_file
            isc_cavity_restart_write = .true.
            ! check if file exists
            inquire(file=isc_cavity_restart_read_file, &
               exist=flag_cavity_restart_file_exists)
            if ( flag_cavity_restart_file_exists ) then
               isc_cavity_restart_read = .true.
               write(info_str,'(2X,2A)') &
                  "Reading cavity from restart file ", isc_cavity_restart_read_file
               call localorb_info(info_str)
               write(info_str,'(2X,A)') &
                  "This file will be overwritten with the new cavity."
               call localorb_info(info_str)
            else
               isc_cavity_restart_read = .false.
               write(info_str,'(2X,2A)') &
                  "Writing cavity to restart file ", isc_cavity_restart_read_file
               call localorb_info(info_str)
            endif

         case('isc_cavity_restart_read')
            read(inputline,*,end=88,err=99) desc_str, isc_cavity_restart_read_file
            ! check if file exists
            inquire(file=isc_cavity_restart_read_file, &
               exist=flag_cavity_restart_file_exists)
            if ( flag_cavity_restart_file_exists ) then
               isc_cavity_restart_read = .true.
               write(info_str,'(2X,2A)') &
                  "Reading cavity from restart file ", isc_cavity_restart_read_file
               call localorb_info(info_str)
            else
               isc_cavity_restart_read = .false.
               write(info_str,'(2X,3A)') &
                  "Specified cavity restart file ", isc_cavity_restart_read_file, &
                        " does not exists."
               call localorb_info(info_str)
               call aims_stop_coll(info_str, func)
            endif

         case('isc_cavity_restart_write')
            read(inputline,*,end=88,err=99) desc_str, isc_cavity_restart_write_file
            isc_cavity_restart_write = .true.
            ! check if file exists
            inquire(file=isc_cavity_restart_read_file, &
               exist=flag_cavity_restart_file_exists)
            if ( flag_cavity_restart_file_exists ) then
               write(info_str,'(2X,3A)') "The cavity restart file ", &
                     isc_cavity_restart_write_file, &
                     " will be overwritten with the new cavity."
               call localorb_info(info_str)
            else
               write(info_str,'(2X,2A)') &
                  "Writing cavity to restart file ", isc_cavity_restart_read_file
               call localorb_info(info_str)
            endif

         case('mpe_solvent_permittivity')
            read(inputline,*,end=88,err=99) desc_str, mpe_solvent_permittivity
            if (mpe_solvent_permittivity.lt.1) then
               write(info_str,"(A)") &
                  "mpe_solvent_permittivity must be at least 1"
               call aims_stop_coll(info_str, func)
            endif
            write(info_str,'(2X,A,F14.4)') &
                  "The permittivity of the solvent is set to ", &
                  mpe_solvent_permittivity
            call localorb_info(info_str)

         case('mpe_bulk_permittivity')
            read(inputline,*,end=88,err=99) desc_str, mpe_bulk_permittivity
            if (mpe_bulk_permittivity.lt.1) then
               write(info_str,"(A)") &
                  "mpe_bulk_permittivity must be at least 1"
               call aims_stop_coll(info_str, func)
            endif

         case('mpe_factorization_type')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
            case("qr")
               mpe_factorization_type = MPE_CONST % FACTZN_QR
               write(info_str,'(2X,2A,1X,A)') &
                  "Factorization method of the MPE equations' coefficient ", &
                  "matrix is:", "QR"
               call localorb_info(info_str)
               write(info_str,'(1X,A)') &
                  "* Warning (MPE): Statistical analysis is only reliable with "// &
                  "mpe_factorization_type qr+svd  or  svd"
               call localorb_info(info_str)
               call aims_stop_coll("Error (MPE): mpe_factorization_type qr "//&
                       "temporarily disabled until statistical analysis is available")
            case("qr+svd")
               mpe_factorization_type = MPE_CONST % FACTZN_QRpSVD
               write(info_str,'(2X,2A,1X,A)') &
                  "Factorization method of the MPE equations' coefficient ", &
                  "matrix is:", "QR, followed by SVD of R"
               call localorb_info(info_str)
            case("svd")
               mpe_factorization_type = MPE_CONST % FACTZN_SVD
               write(info_str,'(2X,2A,1X,A)') &
                  "Factorization method of the MPE equations' coefficient ", &
                  "matrix is:", "SVD"
               call localorb_info(info_str)
            case default
               mpe_factorization_type = MPE_CONST % FACTZN_INVAL
               write(info_str,'(2X,2A,1X,A)') &
                  "Factorization method of the MPE equations' coefficient ", &
                  "matrix is:", "UNKNOWN"
               call aims_stop_coll(info_str, func)
            end select

         case('mpe_skip_first_n_scf_steps')
            read(inputline,*,end=88,err=99) desc_str, mpe_skip_first_n_scf_steps
            if (mpe_skip_first_n_scf_steps.lt.0) then
               call aims_stop_coll("mpe_skip_first_n_scf_steps may not be negative!", func)
            endif

         case('mpe_degree_of_determination')
            read(inputline,*,end=88,err=99) desc_str, mpe_degree_of_determination
            if (mpe_degree_of_determination.lt.1.e0_dp) then
               call aims_stop_coll("mpe_degree_of_determination has to be at least 1!", func)
            endif

         case('mpe_f_sparsity_threshold')
            call aims_stop_coll("Error (MPE): mpe_f_sparsity_threshold "//&
                       "temporarily disabled until statistical analysis is available")
!            read(inputline,*,end=88,err=99) desc_str, mpe_f_sparsity_threshold
!            if (mpe_f_sparsity_threshold.lt.0.e0_dp) then
!               call aims_stop_coll("mpe_f_sparsity_threshold may not be negative!", func)
!            endif

         case('mpe_tol_adjR2')
            read(inputline,*,end=88,err=99) desc_str, mpe_tol_adjR2
            if (mpe_tol_adjR2.lt.0.e0_dp .or. mpe_tol_adjR2.gt.1.e0_dp ) then
               call aims_stop_coll("mpe_tol_adjR2 has to lie in the range [0,1]", func)
            endif

         case('mpe_tol_adjR2_wait_scf')
            read(inputline,*,end=88,err=99) desc_str, mpe_tol_adjR2_wait_scf

         case('mpe_lmax_rf')
            read(inputline,*,end=88,err=99) desc_str, mpe_lmax_rf
            if (mpe_lmax_rf.lt.0) then
               call aims_stop_coll("mpe_lmax_rf may not be negative!", func)
            endif

         case('mpe_lmax_ep')
            read(inputline,*,end=88,err=99) desc_str, mpe_lmax_ep
            if (mpe_lmax_ep.lt.0) then
               call aims_stop_coll("mpe_lmax_ep may not be negative!", func)
            endif

         case('mpe_lmax_ep_from_basis')
            read(inputline,*,end=88,err=99) desc_str, mpe_lmax_ep_from_basis

         case('mpe_lmax_ep_additional_order')
            read(inputline,*,end=88,err=99) desc_str, mpe_lmax_ep_additional_order

         case('mpe_lmax_ep_increase_cluster')
            read(inputline,*,end=88,err=99) desc_str, mpe_lmax_ep_increase_cluster

         case('mpe_n_centers_ep')
            read(inputline,*,end=88,err=99) desc_str, mpe_n_centers_ep
            if (mpe_n_centers_ep.le.0) then
               call aims_stop_coll("mpe_n_centers_ep has to be positive!", func)
            elseif (mpe_n_centers_ep.gt.n_atoms) then
               call aims_stop_coll("mpe_n_centers_ep may not exceed the "//&
                                   "total number of atoms n_atoms!", func)
            endif

         case('mpe_n_boundary_conditions')
            read(inputline,*,end=88,err=99) desc_str, mpe_n_boundary_conditions
            if ((mpe_n_boundary_conditions.ne.2).and.&
                (mpe_n_boundary_conditions.ne.4)) then
               call aims_stop_coll("valid choices for mpe_n_boundary_conditions are 2 or 4!", func)
            endif

         case('isc_record_cavity_creation')
            read(inputline,*,end=88,err=99) desc_str, &
               isc_record_cavity_creation_file, isc_record_every_n_steps

         case('isc_rho_rel_deviation_threshold')
            read(inputline,*,end=88,err=99) desc_str, isc_rho_rel_deviation_threshold
            if (isc_rho_rel_deviation_threshold.le.0) then
               call aims_stop_coll("isc_rho_rel_deviation_threshold has to be positive!", func)
            endif

         case('isc_dt')
            read(inputline,*,end=88,err=99) desc_str, isc_dt
            if (isc_dt.le.0.e0_dp) then
               call aims_stop_coll("isc_dt has to be positive!", func)
            endif

         case('isc_rho_k')
            read(inputline,*,end=88,err=99) desc_str, isc_rho_k
            if (isc_rho_k.le.0.e0_dp) then
               call aims_stop_coll("isc_rho_k has to be positive!", func)
            endif

         case('isc_g_k')
            read(inputline,*,end=88,err=99) desc_str, isc_g_k
            if (isc_g_k.le.0.e0_dp) then
               call aims_stop_coll("isc_g_k has to be positive!", func)
            endif

         case('isc_rep_k')
            read(inputline,*,end=88,err=99) desc_str, isc_rep_k
            if (isc_rep_k.le.0.e0_dp) then
               call aims_stop_coll("isc_rep_k has to be positive!", func)
            endif

         case('isc_gradient_threshold')
            read(inputline,*,end=88,err=99) desc_str, isc_gradient_threshold
            if (isc_gradient_threshold.le.0.e0_dp) then
               call aims_stop_coll("isc_gradient_threshold has to be positive!", func)
            endif

         case('isc_dynamics_friction')
            read(inputline,*,end=88,err=99) desc_str, isc_dynamics_friction
            if (isc_dynamics_friction.lt.0.e0_dp .or. isc_dynamics_friction.gt.1.e0_dp) then
               call aims_stop_coll("isc_dynamics_friction has to be in the range of [0,1]!", func)
            endif

         case('isc_update_nlist_interval')
            read(inputline,*,end=88,err=99) desc_str, isc_update_nlist_interval
            if (isc_update_nlist_interval.le.0) then
               call aims_stop_coll("isc_update_nlist_interval has to be positive!", func)
            endif

         case('isc_max_dyn_steps')
            read(inputline,*,end=88,err=99) desc_str, isc_max_dyn_steps
            if (isc_max_dyn_steps.le.0) then
               call aims_stop_coll("isc_max_dyn_steps has to be positive!", func)
            endif

         case('isc_calculate_surface_and_volume')
            read(inputline,*,end=88,err=99) desc_str, &
                           isc_calculate_surface_and_volume

         case('isc_surface_curvature_correction')
            read(inputline,*,end=88,err=99) desc_str, &
                           isc_surface_curvature_correction

         case('isc_kill_ratio')
            read(inputline,*,end=88,err=99) desc_str, isc_kill_ratio
            if (isc_kill_ratio.lt.0.e0_dp) then
               call aims_stop_coll("isc_kill_ratio may not be negative!", func)
            endif

         case('isc_try_restore_convergence')
            read(inputline,*,end=88,err=99) desc_str, isc_try_restore_convergence

         case('ifp_plane')
                 ifp = .true.
            read(inputline,*,end=88,err=99) desc_str, desc_str
            select case (desc_str)
               case ('x-')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/-1.e0_dp, 0.e0_dp, 0.e0_dp/)
               case ('x+')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/ 1.e0_dp, 0.e0_dp, 0.e0_dp/)
               case ('y-')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/ 0.e0_dp,-1.e0_dp, 0.e0_dp/)
               case ('y+')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/ 0.e0_dp, 1.e0_dp, 0.e0_dp/)
               case ('z-')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/ 0.e0_dp, 0.e0_dp,-1.e0_dp/)
               case ('z+')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist
                  ifp_normal = (/ 0.e0_dp, 0.e0_dp, 1.e0_dp/)
               case ('hesse')
                  read(inputline,*,end=88,err=99) desc_str, desc_str, ifp_dist, ifp_normal(1:3)
                  ! check if normal is not zero vector
                  if ( all(ifp_normal.eq.0.e0_dp) ) &
                     call aims_stop_coll("Specify non-zero normal vector"//&
                        " for the permittivity dividing plane!", func)
               case default
                  ! invalid choice, print short help
                  write(info_str,'(2X,A,A)') 'Supported choices for the ', &
                     'permittivity dividing plane are:'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'x+', &
                     'plane parallel to yz plane at the chosen value of x. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with smaller x',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'x-', &
                     'plane parallel to yz plane at the chosen value of x. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with larger x',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'y+', &
                     'plane parallel to zx plane at the chosen value of y. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with smaller y',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'y-', &
                     'plane parallel to zx plane at the chosen value of y. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with larger y',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'z+', &
                     'plane parallel to xy plane at the chosen value of z. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with smaller z',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'z-', &
                     'plane parallel to xy plane at the chosen value of z. '
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space with larger z',&
                     ' coordinates is considered to be vacuum region', &
                     ' (permittivity = 1)'
                  call localorb_info(info_str)

                  write(info_str,'(4X,A,2X,A)') 'hesse', &
                     'plane in Hesse normal form (dist_to_origin nx ny nz)'
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'the half-space where the',&
                     ' normal vector points to is considered to be vacuum',&
                     ' region (permittivity = 1)'
                  call localorb_info(info_str)
                  write(info_str,'(8X,A,A,A)') 'annotation: the normal vector',&
                     ' will be normalized internally without influence on the',&
                     ' distance.'
                  call localorb_info(info_str)
            endselect
            ifp_dist = ifp_dist / bohr

         case('ifp_rmin')
            read(inputline,*,end=88,err=99) desc_str, ifp_rmin
            if (ifp_rmin.le.0) then
               call aims_stop_coll("ifp_rmin has to be positive!", func)
            endif
            ifp_rmin = ifp_rmin / bohr

         case('ifp_rmax')
            read(inputline,*,end=88,err=99) desc_str, ifp_rmax
            if (ifp_rmax.le.0) then
               call aims_stop_coll("ifp_rmax has to be positive!", func)
            endif
            ifp_rmax = ifp_rmax / bohr

         case('ifp_manual')
            read(inputline,*,end=88,err=99) desc_str, ifp_manual

         case('ifp_rlog')
            read(inputline,*,end=88,err=99) desc_str, ifp_rlog
            if (ifp_rlog.le.0) then
               call aims_stop_coll("ifp_rlog has to be positive!", func)
            endif
            ifp_rlog = ifp_rlog / bohr

         case('ifp_n_shells_log')
            read(inputline,*,end=88,err=99) desc_str, ifp_n_shells_log
            if (ifp_n_shells_log.le.0) then
               call aims_stop_coll("ifp_n_shells_log has to be positive!", func)
            endif

         case('ifp_area_per_point')
            read(inputline,*,end=88,err=99) desc_str, ifp_area_per_point
            if (ifp_area_per_point.le.0) then
               call aims_stop_coll("ifp_area_per_point has to be positive!", func)
            endif
            ifp_area_per_point = ifp_area_per_point / (bohr**2)

         case('ifp_n_angular_max')
            read(inputline,*,end=88,err=99) desc_str, ifp_n_angular_max
            if (ifp_n_angular_max.le.0) then
               call aims_stop_coll("ifp_n_angular_max has to be positive!", func)
            endif

         case('ifp_min_dr')
            read(inputline,*,end=88,err=99) desc_str, ifp_min_dr
            if (ifp_min_dr.le.0) then
               call aims_stop_coll("ifp has to be positive!", func)
            endif
            ifp_min_dr = ifp_min_dr / bohr

         case('ifp_n_shells_lin')
            read(inputline,*,end=88,err=99) desc_str, ifp_n_shells_lin
            if (ifp_n_shells_lin.le.0) then
               call aims_stop_coll("ifp_n_shells_lin has to be positive!", func)
            endif

         case('energy_density')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if ( desc_str .eq. 'harris_foulkes' )  then
               write(info_str,'(2X,A)') 'Calculating Harris-Foulkes-type energy density.'
               call localorb_info(info_str)
               flag_harris_foulkes_energy_density = .true.
               flag_energy_density                = .true.
            elseif ( desc_str .eq. 'chetty_martin' )  then
               write(info_str,'(2X,A)') 'Calculating Chetty-Martin-type energy density.'
               call localorb_info(info_str)
               flag_chetty_martin_energy_density = .true.
               flag_energy_density               = .true.
            else
               write(info_str,'(2X,A,A,A)') " *** Option ",desc_str," unknown for keyword 'energy_density'."
               call localorb_info(info_str)
               call aims_stop_coll(' *** Aborting execution.', func)
            end if

         case('neutral_excitation')
            read(inputline,*,end=88,err=99) desc_str, buffer
            if (buffer.eq."rpa") then
              write(info_str,'(2X,A)') "Calculation of optical excitation energies with RPA will be done."
              call localorb_info(info_str)
              flag_neutral_excitation = .true.
	      neutral_excitation_rpa = .true.
            elseif (buffer.eq."tdhf") then
	      write(info_str,'(2x,a)') "Calculation of optical excitation energies with TDHF will be done."
              call localorb_info(info_str)
              flag_neutral_excitation = .true.
	      neutral_excitation_tdhf = .true.
            elseif (buffer.eq."tddft") then
	      write(info_str,'(2x,a)') "Calculation of optical excitation energies with TDDFT will be done."
              call localorb_info(info_str)
              flag_neutral_excitation = .true.
	      neutral_excitation_tddft = .true.
            elseif (buffer.eq."bse") then
              write(info_str,'(2x,a)') "Calculation of optical excitation energies with BSE will be done."
              call localorb_info(info_str)
              neutral_excitation_bse = .true.
	    else
              write(info_str,'(2X,A)') "ERROR: No calculation mode of optical properties has been set!"
              call localorb_info(info_str)
              call aims_stop()
            end if
            if ( use_geo_relaxation ) then
              write(info_str,'(2X,A)') "ERROR: TDDFT calculations do not work with geometry relaxation."
              call localorb_info(info_str)
              call aims_stop("Please disable keywords that trigger geometry relaxation (e.g. relax_geometry).")
            endif
            ! the correct variable to use here is use_forces, not compute_forces
            if ( use_forces ) then
              write(info_str,'(2X,A)') "ERROR: TDDFT calculations do not work with force calculations."
              call localorb_info(info_str)
              call aims_stop("Please disable keywords that trigger force calculations (e.g. sc_accuracy_forces).")
            endif
! bse related
         case('read_write_qpe')
            read(inputline,*,end=88,err=99) desc_str, buffer
            if (buffer.eq."w") then
              write(info_str,'(2X,A)') "write qpe in file energy_qp"
              call localorb_info(info_str)
              write_qpe = .true.
            elseif (buffer.eq."r") then
              write(info_str,'(2X,A)') "read qpe from file energy_qp"
              call localorb_info(info_str)
              read_qpe = .true.
            elseif ((buffer.eq."wr") .or. (buffer.eq."rw")) then
              write(info_str,'(2X,A)') "write and read qpe from file energy_qp"
              call localorb_info(info_str)
              write_qpe = .true.
              read_qpe = .true.
            endif
         case('bse_s_t')
            read(inputline,*,end=88,err=99) desc_str, buffer
            if (buffer.eq."singlet") then
              write(info_str,'(2X,A)') "calculating BSE singlet states"
              call localorb_info(info_str)
              bse_singlet = .true.
            elseif (buffer.eq."triplet") then
              write(info_str,'(2X,A)') "calculating BSE triplet states"
              call localorb_info(info_str)
              bse_triplet = .true.
            endif
         case('bse_lower_limit')
            read(inputline,*,end=88,err=99) desc_str, bse_lower_limit
!            write(info_str,'(2x, A, I6)') &
!               'BSE: lower limit of states in BSE is set to ', bse_lower_limit
!            call localorb_info(info_str)
         case('bse_upper_limit')
            read(inputline,*,end=88,err=99) desc_str, bse_upper_limit
!            write(info_str,'(2x, A, I6)') &
!               'BSE: upper limit of states in BSE is set to ', bse_upper_limit
!            call localorb_info(info_str)
         case('bse_reduce_matrix')
            bse_reduce_matrix = .true.
!            write(info_str,'(2X, A)') 'BSE reduce matrix is set'
!            call localorb_info(info_str)
         case('bse_reduce_unocc')
            read(inputline,*,end=88,err=99) desc_str, bse_reduce_unocc
!            write(info_str,'(2x,2a,f7.2,a)') &
!               'BSE: Unoccupied states for BSE Matrix will be reduced ', &
!               'to under ', bse_reduce_unocc, 'Ha.'
!            call localorb_info(info_str)
         case('bse_reduce_occ')
            read(inputline,*,end=88,err=99) desc_str, bse_reduce_occ
!            write(info_str,'(2x,2a,f7.2,a)') &
!               'BSE: Occupied states for BSE Matrix will be reduced ', &
!               'to over ', bse_reduce_occ, 'Ha.'
!            call localorb_info(info_str)
         case('band_mulliken_orbit_num')
            read(inputline,*,end=88,err=99) desc_str, band_mulliken_orbit_num
!
         case('excited_states')
	    read(inputline,*,end=88,err=99) desc_str, num_excited_states
	    if(num_excited_states.lt.-1) then
	      if(myid==0) then
	      write(info_str,'(2x,a)') 'The number of excited states that were specified if lower then 0.'
	      call localorb_info(info_str)
	      write(info_str,'(2x,a)') 'That makes no sense, stopping here.'
	      call localorb_info(info_str)
	      call aims_stop('Check the excited_states keyword!')
	      endif
	    endif
	    write(info_str,'(2x,a,i4.1)') 'The number of excited states printed from TDDFT is ', num_excited_states
	    call localorb_info(info_str)
	 case('excited_mode')
	    read(inputline,*,end=88,err=99) desc_str, desc_str
	    if(desc_str.eq.'singlet') then
	      excited_mode_triplet=.false.
	    elseif(desc_str.eq.'triplet') then
	      excited_mode_singlet=.false.
	    elseif(desc_str.eq.'both') then
	      excited_mode_singlet=.true.
	      excited_mode_triplet=.true.
	    else
	      write(info_str,'(a,a)') '**** Error in excited state mode specification: ', desc_str
	      call localorb_info(info_str)
	      call aims_stop('Error in control.in')
	    endif
	 case('tddft_x')
	     if(.not.use_libxc_tddft) then
	      write(info_str,'(2x,a)') 'You cannot specify other f_xc kernels without libxc!'
	      call localorb_info(info_str)
	      call aims_stop('Please install libxc or use the kernel from aims!')
	     else
	      read(inputline,*,end=88,err=99) desc_str, desc_str
              ! Outdated as of LibXC v3.0.0. New function saves our bookkeeping. AJL.
	      ! call determine_func_number(desc_str,libxc_tddft_x)
              libxc_tddft_x = xc_f03_functional_get_number(desc_str)
	      write(info_str,'(2x,a,a,a,i3.1,a)') &
            'You have selected libxc functional ', trim(desc_str), &
            ' (',libxc_tddft_x ,') for TDDFT exchange.'
	      call localorb_info(info_str)
	      if(libxc_tddft_x.eq.0) call aims_stop('Error in tddft_x specification. Please Check!')
	      call check_libxc_func(libxc_tddft_x,.true.,dummy_int,dummy_real)
	     endif
	 case('tddft_c')
	     if(.not.use_libxc_tddft) then
	      write(info_str,'(2x,a)') 'You cannot specify other f_xc kernels without libxc!'
	      call localorb_info(info_str)
	      call aims_stop('Please install libxc or use the kernel from aims!')
	     else
	      read(inputline,*,end=88,err=99) desc_str, desc_str
              ! Outdated as of LibXC v3.0.0. New function saves our bookkeeping. AJL.
	      ! call determine_func_number(desc_str,libxc_tddft_c)
              libxc_tddft_c = xc_f03_functional_get_number(desc_str)
	      write(info_str,'(2x,a,a,a,i3.1,a)') &
            'You have selected libxc functional ', trim(desc_str),' (', &
            libxc_tddft_c ,') for TDDFT correlation.'
	      call localorb_info(info_str)
	      if(libxc_tddft_c.eq.0) call aims_stop('Error in tddft_c specification. Please Check!')
	      call check_libxc_func(libxc_tddft_c,.true.,dummy_int,dummy_real)
	     endif
	 case('tddft_xc')
	     if(.not.use_libxc_tddft) then
	      write(info_str,'(2x,a)') 'You cannot specify other f_xc kernels without libxc!'
	      call localorb_info(info_str)
	      call aims_stop('Please install libxc or use the kernel from aims!')
	     else
	      if(libxc_tddft_chosen_x.or.libxc_tddft_chosen_c) then
	        write(info_str,'(2x,2a)') &
              '***** Error. You have to choose either one X and one C or ', &
              'one XC functional. Both have been set now.'
	        call aims_stop('Please corrent the input file for ( tddft_x .AND. tddft_c ) .XOR. tddft_xc')
	      endif

              ! Doesn't it make more sense to set this as the text definition, and then transfer that in to
              ! a code using the libxc functionals like we do for tddft_x and tddft_c, above?
              ! Updating so this is the case as the setup is otherwise inconsistent, and this option is
              ! definitely more user-friendly. AJL/August 2016

	      !read(inputline,*,end=88,err=99) desc_str, info
	      !libxc_tddft_xc = info
	      !write(info_str,'(2x,a,i3.1,a)') &
            !'You have selected libxc functional number ', libxc_tddft_xc , &
            !' for TDDFT exchange and correlation.'

              read(inputline,*,end=88,err=99) desc_str, desc_str
              libxc_tddft_xc = xc_f03_functional_get_number(desc_str)
              write(info_str,'(2x,a,a,a,i3.1,a)') &
            'You have selected libxc functional ', trim(desc_str),' (', &
            libxc_tddft_xc ,') for TDDFT exchange and correlation.'

	      call localorb_info(info_str)
	      call check_libxc_func(libxc_tddft_xc,.true.,dummy_int,dummy_real)
	     endif
	 case('tddft_kernel')
             if (use_libxc_tddft) then
               ! Update: I've moved this into a subroutine in optical_reponse/libxc_function.f90
               !         I still need to test if it is needed but now the pointers are in this part of the code
               !         What this does do is catch anyone trying to use LibXC without linking to it

               call initialise_libxc()

               ! call xc_f03_func_init(xc_func, xc_info, 1, XC_UNPOLARIZED)
               ! call xc_f03_func_end(xc_func)
             else
               read(inputline,*,end=88,err=99) desc_str, desc_str
               if ( desc_str .eq. 'pw-lda' )  then
                 write(info_str,'(2x,a,a)') &
                   'Kernel of Perdew-Wang parametrisation of ', &
                   'Ceperley-Alder LDA will be used for TDDFT calculations.'
	         call localorb_info(info_str)
                 flag_fxc = 8 ! pw-lda
               elseif ( desc_str .eq. 'pz-lda' )  then
                 write(info_str,'(2x,a,a)') &
                   'Kernel of Perdew-Zunger parametrisation of ', &
                   'Ceperley-Alder LDA will be used for TDDFT calculations.'
	         call localorb_info(info_str)
                 flag_fxc = 3 ! pz-lda
               else
                 write(info_str,'(2x,a,a)') &
                   'At the moment, the only built-in kernels for ', &
                   'TDDFT calculations are pw-lda or pz-lda.'
	         call localorb_info(info_str)
                 call aims_stop('Please select either pw-lda or pz-lda kernel.')
               endif
             endif
         case('casida_reduce_matrix')
            read(inputline,*,end=88,err=99) desc_str, casida_reduce_matrix
!            if(.not.casida_reduce_matrix.or.casida_reduce_matrix) then
!              write(info_str,'(2x,a)') 'Error in control.in. Please set casida_reduce_matrix to .true. or .false. only!'
!              call localorb_info(info_str)
!              call aims_stop('*** ABORT - Syntex error in control.in')
!            endif
         case('casida_reduce_occ')
            read(inputline,*,end=88,err=99) desc_str, casida_occ_limit
            if(casida_occ_limit.lt.0) then
              write(info_str,'(2x,2a,f7.2,a)') &
                 'TDDFT: Occupied states for Omega Matrix will be reduced ', &
                 'to over ', casida_occ_limit, 'Ha.'
              call localorb_info(info_str)
            endif
         case('casida_reduce_unocc')
            read(inputline,*,end=88,err=99) desc_str, casida_unocc_limit
            if(casida_unocc_limit.gt.0) then
              write(info_str,'(2x,a,f7.2,a)') 'TDDFT: Virtual states for Omega Matrix will be reduced to under ', casida_unocc_limit, 'Ha.'
              call localorb_info(info_str)
            endif
         case('casida_setup_scalapack')
            casida_setup_scalapack = .true.
         ! invoke RRS-PBC subroutine, igor
         case('rrs_pbc')
            write(info_str,'(2X,A)') &
               'RRS-PBC :: invoke RRS-PBC calculation.'
            call localorb_info(info_str)

         case('rrs_pbc_unit_cell')
            read(inputline,*,end=88,err=99) desc_str, &
                rrs_pbc_n_center_atom, rrs_pbc_n_equal
            call read_rrs_pbc_unit_cell(7,0)
            write(info_str,'(2X,A)') &
               'RRS-PBC :: loading unit cell index'
            call localorb_info(info_str)

         case('rrs_pbc_k_grid')
            read(inputline,*,end=88,err=99) desc_str, &
                rrs_pbc_n_k_points_xyz(1), &
                rrs_pbc_n_k_points_xyz(2), &
                rrs_pbc_n_k_points_xyz(3)
            write(info_str,'(2X,A)') &
               'RRS-PBC :: loading k grids'
            call localorb_info(info_str)
         ! Frozen core implementation for post-SCF calculations, including RPA,
         ! MP2 and so on.
         case('frozen_core_postscf','frozen_core_postSCF')
            read(inputline,*,end=88,err=99)desc_str, n_frozen_shell
            write(info_str,'(2X,A,I6)') &
            'Core electrons are frozen in the post-SCF calculation :', &
            n_frozen_shell
            flag_frozen_core_postSCF = .true.
            call localorb_info(info_str,use_unit,'(A)')
         ! Frozen virtual implementation for post-SCF calculations, including RPA,
         ! MP2 and so on.
         case('frozen_virtual_postSCF')
            read(inputline,*,end=88,err=99)desc_str, n_frozen_virtual_states
            write(info_str,'(2X,A,I6)') &
            'Several highest virtual states are frozen in the post-SCF calculation :', &
            n_frozen_virtual_states
            flag_frozen_virtual_postSCF = .true.
            call localorb_info(info_str,use_unit,'(A)')
         ! EN-shift on MP2, RPA and GW and so on. Igor
         ! MP2 and so on.
         ! en_shift_type = 1: for hole-hole EN2
         !                 2: for constant energy shift
         !                 3: for (screening) coupling PT2
         !                 4: for (screening) coupling PT2 + the first-order term
         case('en_shift')
            read(inputline,*,end=88,err=99)desc_str, en_shift_type
            flag_en_shift       = .true.
            write(info_str,'(2X,A,I6)') &
            'Energy shifting on the denominators: ', &
            en_shift_type
            call localorb_info(info_str,use_unit,'(A)')
         ! en_shift_constant  for en_shift = 2
         case('en_shift_constant')
            read(inputline,*,end=88,err=99)desc_str, en_shift_constant
            write(info_str,'(2X,A,F16.8)') &
            'Energy shifting amount on the denominators: ', &
            en_shift_constant
            call localorb_info(info_str,use_unit,'(A)')
         ! coupling_pt2_factor,
         ! coupling_pt2_screen and coupling_pt2_shift for en_shift = 3
         case('coupling_pt2_factor')
            read(inputline,*,end=88,err=99)desc_str, coupling_pt2_factor
            write(info_str,'(2X,A,F16.8)') &
            'Coupling factor for PT2 on the denominators: ', &
            coupling_pt2_factor
            call localorb_info(info_str,use_unit,'(A)')
         case('coupling_pt2_screen')
            read(inputline,*,end=88,err=99)desc_str, coupling_pt2_screen
            write(info_str,'(2X,A,F16.8)') &
            'Screening factor for coupling PT2: ', &
            coupling_pt2_screen
            call localorb_info(info_str,use_unit,'(A)')
         case('coupling_pt2_shift')
            read(inputline,*,end=88,err=99)desc_str, coupling_pt2_shift
            write(info_str,'(2X,A,F16.8)') &
            'Shift for screening effect in coupling PT2: ', &
            coupling_pt2_shift
            call localorb_info(info_str,use_unit,'(A)')
         ! rescaled coefficients for scsRPA
         case('scsrpa_coeff')
            read(inputline,*,end=88,err=99)desc_str, scsrpa_coeff(1),scsrpa_coeff(2)
            write(info_str,'(2X,A,2F12.3)') &
            'Coefficients for spin-component-scaled RPA (scs-RPA): ', &
            scsrpa_coeff(1),scsrpa_coeff(2)
         ! approximate renormalization of energy gaps for ospt2 and osrpa
         case('renormalized_eg')
            read(inputline,*,end=88,err=99)desc_str, renormalized_eg(1),renormalized_eg(2),renormalized_eg(3),renormalized_eg(4)
            write(info_str,'(2X,A,4F12.3)') &
            'Approximate renormalization on the energy gaps for osRPA: ', &
            renormalized_eg(1),renormalized_eg(2),renormalized_eg(3), &
            renormalized_eg(4)
            call localorb_info(info_str,use_unit,'(A)')
         case('osrpa_threshold')
            read(inputline,*,end=88,err=99)desc_str, c_osrpa_threshold
         ! evaluate osrpa in higher orders
         case('osrpa_orders')
            read(inputline,*,end=88,err=99)desc_str, c_osrpa_order
            write(info_str,'(2X,A,I3,A)') &
            'Evaluate the osRPA correlation up to the ',c_osrpa_order,'th order.'
            call localorb_info(info_str,use_unit,'(A)')
            !if (c_osrpa_order.eq.30) then
            !    write(info_str,'(2X,A,A)') &
            !    "At present, we do not want to evaluate the osRPA correlation more than 30th order. ",&
            !    "In consequence, the requirement of higher-order osRPA correlation is rejected."
            !    c_osrpa_order = 30
            !    call localorb_info(info_str,use_unit,'(A)')
            !end if
         ! Number of Gauss-Legendre grid for the coupling-constant integration in osRPA
         case('n_lambda_osrpa')
            read(inputline,*,end=88,err=99)desc_str, n_lambda_osrpa
            write(info_str,'(2X,A,I3)') &
            'Number of Gauss-Legendre grid for the coupling-constant integration in osRPA: ', &
            n_lambda_osrpa
            call localorb_info(info_str,use_unit,'(A)')
         ! the highest excitation for the truncated configuration
         ! interaction (CI) calculation
         case('ci_truncation')
            read(inputline,*,end=88,err=99)desc_str, ci_truncation
            write(info_str,'(2X,A,I1,A)') &
            'Configuration interaction calculations consider the configurations with at most ', ci_truncation, &
            '-electron excitations'
            call localorb_info(info_str,use_unit,'(A)')
         ! determine whether to consider the single excitation configurations
         case('ci_single_excitations')
            read(inputline,*,end=88,err=99)desc_str, flag_single_excitation
            if (flag_single_excitation) then
                write(info_str,'(2X,A,A)') &
                'The contributions of single-excitation configurations are ', &
                'taken into account in CI calculations'
                call localorb_info(info_str,use_unit,'(A)')
            else
                write(info_str,'(2X,A,A)') &
                'The contributions of single-excitation configurations are not ', &
                'taken into account'
                call localorb_info(info_str,use_unit,'(A)')
            endif
         ! determine whether to store the two-electron integrals for alpha spin
         ! channel only
         case('ci_single_spin_channel')
            read(inputline,*,end=88,err=99)desc_str, flag_single_spin_channel
            if (flag_single_spin_channel) then
                write(info_str,'(2X,A,A)') &
                'Only two-electron integrals for alpha spin channel are stored', &
                ' in FCIDUMP (similar to ROHF) for the next FCI-QMC simulation'
                call localorb_info(info_str,use_unit,'(A)')
            else
                write(info_str,'(2X,A,A)') &
                'Two-electron integrals for both alpha and beta spin channel are stored', &
                ' in FCIDUMP (UHF) for the next FCI-QMC simulation'
                call localorb_info(info_str,use_unit,'(A)')
            endif
         ! the maximum iteraction number for the CI calculation
         case('n_scf_ci')
            read(inputline,*,end=88,err=99)desc_str, n_scf_ci
            write(info_str,'(2X,A,I6)') &
            'The maximum iteraction number for the CI calculation is ', &
            n_scf_ci
            call localorb_info(info_str,use_unit,'(A)')
         ! the maximum iteraction number for the CI calculation
         case('ci_accuracy_etot')
            read(inputline,*,end=88,err=99)desc_str, ci_accuracy_etot
            write(info_str,'(2X,A,E13.6)') &
            'The convergence criterion for the CI total energy is ', &
            ci_accuracy_etot
            call localorb_info(info_str,use_unit,'(A)')
         ! the criterion for the printout of configuration coefficients
         case('ci_coeff_threshold')
            read(inputline,*,end=88,err=99)desc_str, ci_coeff_threshold
            write(info_str,'(2X,A,E13.6)') &
            'The configuration coefficients are printout, if the absolute value is larger than ', &
            ci_coeff_threshold
            call localorb_info(info_str,use_unit,'(A)')
         ! the flag to turn on the calculations of Coupled-cluster methods
         case('cc_method')
             read(inputline,*,end=88,err=99) desc_str,CC_method
             if ((CC_method.eq.'ccsd').or.(CC_method.eq.'CCSD')) then
               CC_calc_method = 1
               write(info_str,'(2X,A)') 'Coupled cluster calculation: CCSD'
             else if ((CC_method.eq.'CCSD(T)').or.(CC_method.eq.'ccsd(t)')) then
               CC_calc_method = 2
               write(info_str,'(2X,A)') 'Coupled cluster calculation: CCSD(T)'
             end if
             call localorb_info(info_str,use_unit,'(A)')
         case('cc_general')
             flag_cc_general = .true.
             write(info_str,'(2X,A)') 'CC_gerenal'
             call localorb_info(info_str,use_unit,'(A)')
         case('cc_collective')
             CC_collective = .true.
         case('cc_solver')
             read(inputline,*,end=88,err=99) desc_str,CC_solver
         case('cc_use_disk')
             CC_use_disk = .true.
         case('cc_scf_max_cycle')
             read(inputline,*,end=88,err=99) desc_str,CC_max_cycle
             write(info_str,'(2X,A,I4)') 'CC calculation max scf cycles:', &
                                         CC_max_cycle
             call localorb_info(info_str,use_unit,'(A)')
         case('cc_acc_rle')
             CC_acc_method = 1
             write(info_str,'(2X,A,I4)') 'CC acceleration method: RLE'
         case('cc_acc_diis')
             CC_acc_method = 2
             write(info_str,'(2X,A,I4)') 'CC acceleration method: DIIS'
         case('cc_rle_vector_save')
             read(inputline,*,end=88,err=99) desc_str,CC_n_RLE_sv
             write(info_str,'(2X,A,I4)') 'CC calculation RLE basis saved', &
                                         CC_n_RLE_sv
         case('cc_rle_thresh')
             read(inputline,*,end=88,err=99) desc_str,CC_RLE_bas_thresh
             write(info_str,'(2X,A,E15.8)') 'CC calculation RLE basis threshold', &
                                         CC_RLE_bas_thresh
         case('cc_converg_thresh')
             read(inputline,*,end=88,err=99) desc_str,CC_converg_thresh
             write(info_str,'(2X,A,E15.8)') 'CC calculation converg threshold', &
                                         CC_converg_thresh
         case('cc_n_domain')
             read(inputline,*,end=88,err=99) desc_str,CC_n_domain
             write(info_str,'(2X,A,I6,A)') 'MPI tasks are divided into ', &
                                            CC_n_domain, 'groups'
         case('cc_mpi_max_len')
             read(inputline,*,end=88,err=99) desc_str,CC_MPI_max_len
             write(info_str,'(2X,A,I10)') 'Max length of MPI array is set to', &
                                            CC_MPI_max_len
         case('cc_work_tmp_size')
             read(inputline,*,end=88,err=99) desc_str,CC_work_tmp_size
             write(info_str,'(2X,A,I10)') 'Temporary array size is set to', &
                                            CC_work_tmp_size
          case('cc_sv_control')
             read(inputline,*,end=88,err=99) desc_str,CC_sv_control
             write(info_str,'(2X,A,A,A)') 'Using ', CC_sv_control, 'saving strategy'
          case('cc_restart')
             read(inputline,*,end=88,err=99) desc_str,CC_restart_step
             CC_restart_flag = .true.
             write(info_str,'(2X,A)') 'CC restart'
          case('cc_pt_n_out')
              read(inputline,*,end=88,err=99) desc_str,CC_PT_n_out
              write(info_str,'(2X,A,I3)') 'Number of outsteps in PT calculation: ', &
                                          CC_PT_n_out
          case('cc_diis_step')
              read(inputline,*,end=88,err=99) desc_str,CC_DIIS_step
              write(info_str,'(2X,A,I4)') 'CC DIIS_step',CC_DIIS_step
          case('cc_diis_n_sv')
              read(inputline,*,end=88,err=99) desc_str,CC_DIIS_n_sv
              write(info_str,'(2X,A,I4)') 'CC DIIS_n_sv',CC_DIIS_n_sv
          case('cc_init_level_shift')
              read(inputline,*,end=88,err=99) desc_str,CC_lv_shift
              write(info_str,'(2X,A,F14.8,A)') 'CC lv_shift',CC_lv_shift,' Ha.'
          case('cc_check')
              CC_check_flag = .true.
              write(info_str,'(2X,A)') 'CC check'
          case('cc_abcd_not_save')
              CC_abcd_sv_flag = .false.
              write(info_str,'(2X,A)') 'CC integral abcd tensor is not saved'
          case('cc_gamma')
              CC_gamma_flag = .true.
              write(info_str,'(2X,A)') 'Gamma only calculations are taken by cluster model'
          case('cc_mem_check')
              CC_mem_check = .true.
              write(info_str,'(2X,A)') 'Just a memory check for CC calculation'
          case('cc_super_system')
              CC_super_system_flag = .true.
              write(info_str,'(2X,A)') 'CC calculation for extremely large system'
         ! the flag to turn on the calculations of RPA potentials along the adiabatic connection path,
         case('rpa_along_ac_path')
            read(inputline,*,end=88,err=99)desc_str, rpa_along_ac_path_grid
            write(info_str,'(2X,A,I3,A)') &
            'Sampling RPA potentials along the adiabatic connection path with ', &
            rpa_along_ac_path_grid, ' grid points.'
            call localorb_info(info_str,use_unit,'(A)')
        ! the flag to specify the dft components for printing out
        case('printout_dft_components')
            read(inputline,*,end=88,err=99)desc_str, desc_str
            n_dft_methods = n_dft_methods + 1
            if (n_dft_methods .gt. 10) then
                write(info_str,'(2X,A,A,A)') &
                    'WARNNING: The number of DFT methods for postprocessing XC evaluation is in excess of 10!!! Therefore, ', &
                    desc_str, ' will be ignored.'
            end if
            if (desc_str .eq. 'pbe' .or. desc_str .eq. 'PBE') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: PBE"
               end if
               dft_print_out(n_dft_methods) = 6
           else if (desc_str .eq. 'blyp') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: BLYP"
               end if
               dft_print_out(n_dft_methods) = 9
           else if (desc_str .eq. 'vwn-gauss') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: vwn-gauss"
               end if
               dft_print_out(n_dft_methods) = 16
           else if (desc_str .eq. 'xpbe') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: xPBE"
               end if
               dft_print_out(n_dft_methods) = 22
           else if (desc_str .eq. 'scan'.or.desc_str.eq.'SCAN') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: SCAN"
               end if
               dft_print_out(n_dft_methods) = 59
           else if (desc_str .eq. 'tpss'.or.desc_str.eq.'TPSS') then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') &
                     "XC for printing out: TPSS"
               end if
               dft_print_out(n_dft_methods) = 51
           end if
        ! turn on the restarting utility for RPA, Igor
        case('restart_rpa')
           read(inputline,*,end=88,err=99) desc_str, restart_rpa_read_file
           !keep_restart_info_rpa  = .true.
           restart_rpa_read       = .true.
           restart_rpa_write      = .true.
           restart_rpa_write_file = restart_rpa_read_file
           write(info_str,'(2X,2A)') 'Restart file specified as : ', &
                    restart_rpa_read_file
           call localorb_info(info_str,use_unit,'(A)')
           inquire(FILE=restart_rpa_read_file,EXIST=restart_rpa_file_exists)
           if (restart_rpa_file_exists) then
              write(info_str,'(2X,2A)') '| Reading information from ', &
                    'existing file. '
           else
              write(info_str,'(2X,2A)') '| This file does not exist, ', &
                    'will create it after the first complete frequence calculation.'
              restart_rpa_read = .false.
           end if
           call localorb_info(info_str,use_unit,'(A)')
        ! turn on the restarting utility for PT2, Igor
        case('restart_pt2')
           read(inputline,*,end=88,err=99) desc_str, restart_pt2_read_file
           !keep_restart_info_rpa  = .true.
           restart_pt2_read       = .true.
           restart_pt2_write      = .true.
           write(info_str,'(2X,2A)') 'Restart file specified as : ', &
                    restart_pt2_read_file
           restart_pt2_write_file = restart_pt2_read_file
           call localorb_info(info_str,use_unit,'(A)')
           inquire(FILE=restart_pt2_read_file,EXIST=restart_pt2_file_exists)
           if (restart_pt2_file_exists) then
              write(info_str,'(2X,2A)') '| Reading information from ', &
                    'existing file. '
           else
              write(info_str,'(2X,2A)') '| This file does not exist, ', &
                    'will create it after a batch of calculations.'
              restart_pt2_read = .false.
           end if
           call localorb_info(info_str,use_unit,'(A)')
        ! parameters for explict range separation scheme (ERS), Igor
        case('ers_width')
           read(inputline,*,end=88,err=99) desc_str, ers_width
           write(info_str,'(2X,A,F16.12,A)') &
           & 'Setting explict range separation width to', ers_width, ' * rcut.'
           call localorb_info(info_str)
           use_ers = .true.
        case('ers_rcut')
           read(inputline,*,end=88,err=99) desc_str, ers_rcut
           ers_rcut = ers_rcut / bohr
           write(info_str,'(2X,A,F16.12,A)') &
           & 'Setting explict range separation shift to', &
           & ers_rcut*bohr, ' A.'
           call localorb_info(info_str)
           use_ers = .true.
        case('ers_omega')
           read(inputline,*,end=88,err=99) desc_str, ers_omega
           write(info_str,'(2X,A,F16.12,A)') &
           & 'Setting HSE width for ERS to', &
           & ers_omega, ' omega.'
           call localorb_info(info_str)
           use_ers = .true.
        case('erfc_omega_postSCF')
           read(inputline,*,end=88,err=99) desc_str, erfc_omega
           write(info_str,'(2X,A,F16.12,A)') &
           & 'Setting Erfc width ', &
           & erfc_omega, ' omega.'
           call localorb_info(info_str)
           use_erfc = .true.
         ! to skip the SCF procedure for post-SCF calculations only
         case('skip_scf_for_postscf')
           read(inputline,*,end=88,err=99) desc_str, skip_scf_for_postscf
           if (skip_scf_for_postscf) then
              write(info_str,'(2X,A)') &
                  ' * Warning: Skipping the SCF cycles for postSCF calculations.'
              call localorb_info(info_str,use_unit,'(A)')
           end if
         case('restart_for_postscf')
           read(inputline,*,end=88,err=99) desc_str, restart_for_postscf
           if (restart_for_postscf) then
              write(info_str,'(2X,A)') &
                  ' * Note: generate the restart file for postscf calculations.'
              call localorb_info(info_str,use_unit,'(A)')
           end if
         ! to froze virtual orbitals for post-SCF calculations
         case('frozen_virtual_file')
            read(inputline,*,end=88,err=99) desc_str, fv_orbs_file
            fv_filter  = .true.
            write(info_str,'(2X,2A)') 'Read in the frozen virtual orbitals from ', &
                     fv_orbs_file
            call localorb_info(info_str,use_unit,'(A)')
            inquire(FILE=fv_orbs_file,EXIST=restart_pt2_file_exists)
            if (restart_pt2_file_exists) then
               write(info_str,'(2X,2A)') '| Reading information from ', &
                     'existing file. '
               call localorb_info(info_str,use_unit,'(A)')
            else
               write(info_str,'(2X,3A)') '| The file of frozen virtual orbitals ', &
                   'does not exist:', fv_orbs_file
               call aims_stop(info_str)
            end if
         ! k-point parallelism logic that avoids splitting k-points over SMP nodes
         case('restrict_kpoint_to_smp_node')
            read(inputline,*,end=88,err=99) desc_str, restrict_kpoint_to_smp_node
            if (restrict_kpoint_to_smp_node) then
               write(info_str,'(2X,A)') &
                  'Restriction of each k-point to a single SMP node requested.'
               call localorb_info(info_str)
            end if
            flag_restrict_kpoint_to_smp_node = .true.

         ! maximum tasks per SMP node
         case('max_tasks_per_smp_node')
            read(inputline,*,end=88,err=99) desc_str, max_tasks_per_smp_node
            if (restrict_kpoint_to_smp_node) then
               write(info_str,'(2X,A,I10)') &
                  'Maximum number of tasks per node set to ', max_tasks_per_smp_node
               call localorb_info(info_str)
            end if
            flag_max_tasks_per_smp_node = .true.

         case ('python_hook')
             buffer_arr(:) = ''
             read (inputline, *, iostat=i_code) desc_str, buffer_arr
             if (i_code > 0) go to 99
             i_code = register_python_hook(buffer_arr(1), buffer_arr(2), buffer_arr(3))
             if (i_code >= 0) go to 99

         case('memtrace_use')

         case('memtrace_summary')

         case('memtrace_targetdir')

         ! Use GPU-accelerated methods that are considered stable
         case('use_gpu')
            use_gpu = .true.
            use_gpu_density = .true.
            use_gpu_hamiltonian = .true.
            use_gpu_forces = .true.

         ! Use GPU acceleration for density-matrix-based charge density update
         case('gpu_density')
            use_gpu = .true.
            use_gpu_density = .true.

         ! Use GPU acceleration for real-space Hamiltonian matrix integration
         case('gpu_hamiltonian')
            use_gpu = .true.
            use_gpu_hamiltonian = .true.

         ! Use GPU acceleration for calculation of various contributions to forces
         ! and analytical stress tensor
         case('gpu_forces')
            use_gpu = .true.
            use_gpu_forces = .true.

         ! Tag, set by user, to identify this run as belonging to a larger group
         ! of calculations
         case('calculation_set')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if(allocated(calculation_set)) deallocate(calculation_set)
            calculation_set = trim(adjustl(desc_str))

         ! Tag, set by user, to identify this run as belonging to a particular
         ! subset of a larger group of calculations
         case('calculation_subset')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if(allocated(calculation_subset)) deallocate(calculation_subset)
            calculation_subset = trim(adjustl(desc_str))

         case('check_cpu_consistency')
            read(inputline,*,end=88,err=99) desc_str,desc_str
            if (desc_str.eq.'.true.') then
               flag_cpu_consistency = .true.
               write(info_str, '(2X,A)') &
         &     "Consistency of geometry between different MPI tasks will be checked."
               call localorb_info(info_str)
            else if (desc_str.eq.'.false.') then
               flag_cpu_consistency = .false.
               write(info_str, '(2X,A)') &
         &     "Consistency of geometry between different MPI tasks will NOT be checked."
               call localorb_info(info_str)
            else
               write(info_str, '(1X,A)') &
         &     "*** Error in file control.in: Unknown value for keyword 'check_cpu_consistency'."
               call localorb_info(info_str)
               write(info_str, '(1X,A)') &
         &     "*** Keyword check_cpu_consistency can only be '.true.' or '.false'."
               call localorb_info(info_str)
               call aims_stop('Unknown value for keyword check_cpu_consistency in control.in.')
            end if

         case('cpu_consistency_threshold')
            read(inputline,*,end=88,err=99) desc_str,cpu_consistency_threshold
               write(info_str, '(2X,A, E16.8,A)') &
         &     "cpu_consistency_threshold set to ", cpu_consistency_threshold, '.'
               call localorb_info(info_str)

         case('scalapack_block_size')
            read(inputline,*,end=88,err=99) desc_str, scalapack_block_size
            if(scalapack_block_size .le. 0) then
               write(info_str, '(1X,A)') &
                  "*** Error in file control.in: 'scalapack_block_size' must be positive."
               call localorb_info(info_str)
               call aims_stop('Invalid value for keyword scalapack_block_size.')
            else
               write(info_str,'(2X,A,I6)') &
                  "ScaLAPACK block size set to ", scalapack_block_size
               call localorb_info(info_str)
            endif

         case('ELPA_AEO')
            read(inputline,*,end=88,err=99) desc_str, elpa_aeo
            if(elpa_aeo.lt.4) then
               if(myid==0) then
                  write(info_str,'(2x,a)') 'Run ELPA in single precision mode.'
                  call localorb_info(info_str)
               end if
            else
               if(myid==0) then
                  write(info_str,'(2x,a)') 'Run ELPA in double precision mode.'
                  call localorb_info(info_str)
               end if
            end if

         case('elsi_method')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if(desc_str.eq.'ev') then
               use_elsi_ev = .true.
               use_elsi_dm = .false.
            elseif(desc_str.eq.'dm') then
               use_elsi_ev = .false.
               use_elsi_dm = .true.
            else
               call aims_stop('Invalid elsi_method.')
            endif

            flag_elsi_method = .true.

         case('elsi_solver')
            read(inputline,*,end=88,err=99) desc_str, desc_str

            select case(desc_str)
            case('elpa','ELPA')
               elsi_solver = 1
            case('omm','OMM','libOMM')
               elsi_solver = 2
            case('pexsi','PEXSI')
               elsi_solver = 3
            case('eigenexa','EigenExa')
               elsi_solver = 4
            case('sips','slepc','SLEPc')
               elsi_solver = 5
            case('ntpoly','NTPoly')
               elsi_solver = 6
            case('magma','MAGMA')
               elsi_solver = 7
            case('auto','AUTO')
               elsi_solver = 0
            case default
               elsi_solver = 0
               call aims_stop('Invalid elsi_solver.')
            end select

         case('elsi_restart')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if(desc_str.eq.'write') then
               elsi_write_dm = .true.
               read(inputline,*,end=88,err=99) desc_str, desc_str, elsi_write_dm_freq
            elseif(desc_str.eq.'read') then
               elsi_read_dm = .true.
            elseif(desc_str.eq.'read_and_write') then
               elsi_write_dm = .true.
               read(inputline,*,end=88,err=99) desc_str, desc_str, elsi_write_dm_freq

               ! Read only if restart file exists
               write(dm_file,'(A,I2.2,A,I6.6,A)') "D_spin_",1,"_kpt_",1,".csc"
               inquire(file=dm_file,exist=elsi_read_dm)
            else
               call aims_stop('Invalid elsi_restart.')
            endif

            output_matrix_io_timing = .true.

         case('elsi_restart_use_overlap')
            read(inputline,*,end=88,err=99) desc_str, elsi_extrap_dm

         case('elsi_filter')
            read(inputline,*,end=88,err=99) desc_str, elsi_filter

         case('elsi_rw_filter')
            read(inputline,*,end=88,err=99) desc_str, elsi_rw_filter

         case('elsi_output')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if(desc_str.eq.'none') then
               elsi_out_level = 0
            elseif(desc_str.eq.'light') then
               elsi_out_level = 1
            elseif(desc_str.eq.'detail') then
               elsi_out_level = 2
            elseif(desc_str.eq.'debug') then
               elsi_out_level = 3
            elseif(desc_str.eq.'json') then
               elsi_out_json = 1
            else
               elsi_out_level = 1
            endif

         case('elsi_elpa_solver')
            read(inputline,*,end=88,err=99) desc_str, elsi_elpa_solver
            if((elsi_elpa_solver < 1) .or. (elsi_elpa_solver > 2)) then
               call aims_stop('Invalid elsi_elpa_solver.')
            endif

         case('elsi_elpa_n_single')
            read(inputline,*,end=88,err=99) desc_str, elsi_elpa_n_single

         case('elsi_elpa_gpu')
            read(inputline,*,end=88,err=99) desc_str, elsi_elpa_gpu

         case('elsi_omm_n_elpa')
            read(inputline,*,end=88,err=99) desc_str, elsi_omm_n_elpa
            if(elsi_omm_n_elpa < 0) then
               call aims_stop('Invalid elsi_omm_n_elpa.')
            endif

         case('elsi_omm_flavor')
            read(inputline,*,end=88,err=99) desc_str, elsi_omm_flavor
            if((elsi_omm_flavor /= 0) .and. (elsi_omm_flavor /= 2)) then
               call aims_stop('Invalid elsi_omm_flavor.')
            endif

         case('elsi_omm_tol')
            read(inputline,*,end=88,err=99) desc_str, elsi_omm_tol
            if(elsi_omm_tol < 0.d0) then
               call aims_stop('Invalid elsi_omm_tol.')
            endif

         case('elsi_pexsi_np_per_pole')
            read(inputline,*,end=88,err=99) desc_str, elsi_pexsi_np_per_pole
            if(elsi_pexsi_np_per_pole < 0) then
               call aims_stop('Invalid elsi_pexsi_np_per_pole.')
            endif

         case('elsi_pexsi_np_symbo')
            read(inputline,*,end=88,err=99) desc_str, elsi_pexsi_np_symbo
            if(elsi_pexsi_np_symbo < 0) then
               call aims_stop('Invalid elsi_pexsi_np_symbo.')
            endif

         case('elsi_eigenexa_method')
            read(inputline,*,end=88,err=99) desc_str, elsi_eigenexa_method
            if((elsi_eigenexa_method < 1) .or. (elsi_eigenexa_method > 2)) then
               call aims_stop('Invalid elsi_eigenexa_method.')
            endif

         case('elsi_sips_n_slice')
            read(inputline,*,end=88,err=99) desc_str, elsi_sips_n_slice
            if(elsi_sips_n_slice < 0) then
               call aims_stop('Invalid elsi_sips_n_slice.')
            endif

         case('elsi_sips_slice_type')
            read(inputline,*,end=88,err=99) desc_str, elsi_sips_slice_type
            if(elsi_sips_slice_type < 0) then
               call aims_stop('Invalid elsi_sips_slice_type.')
            endif

         case('elsi_sips_n_elpa')
            read(inputline,*,end=88,err=99) desc_str, elsi_sips_n_elpa
            if(elsi_sips_n_elpa < 0) then
               call aims_stop('Invalid elsi_sips_n_elpa.')
            endif

         case('elsi_ntpoly_method')
            read(inputline,*,end=88,err=99) desc_str, elsi_ntpoly_method
            if((elsi_ntpoly_method < 0) .or. (elsi_ntpoly_method > 3)) then
               call aims_stop('Invalid elsi_ntpoly_method.')
            endif

         case('elsi_ntpoly_tol')
            read(inputline,*,end=88,err=99) desc_str, elsi_ntpoly_tol
            if(elsi_ntpoly_tol < 0.d0) then
               call aims_stop('Invalid elsi_ntpoly_tol.')
            endif

         case('elsi_ntpoly_filter')
            read(inputline,*,end=88,err=99) desc_str, elsi_ntpoly_filter

         case('elsi_magma_solver')
            read(inputline,*,end=88,err=99) desc_str, elsi_magma_solver
            if((elsi_magma_solver < 1) .or. (elsi_magma_solver > 2)) then
               call aims_stop('Invalid elsi_magma_solver.')
            endif

         case("debug_module")
            read(inputline,*,end=88,err=99) desc_str, desc_str
            call activate_debugging(desc_str)

         case("global_memory_tracking")
            read(inputline,*,end=88,err=99) desc_str, aims_mem_sync

            if (.not. aims_mem_sync) then
               if (myid == 0) then
                  write (use_unit,'(2X,A)') &
                    "FHI-aims will only report memory usage on task 0."
               end if
            end if

         case('switch_rho_multipole')
            read(inputline,*,end=88,err=99) desc_str,switch_rho_multipole

         if (switch_rho_multipole) then
            write(info_str, '(2X,A)') &
         &  "This run will use rho_multipole and rho_multipole_gradient for rho and rho_gradient."
            call localorb_info(info_str)
         end if

         case('force_single_restartfile')
            read(inputline,*,end=88,err=99) desc_str, force_single_restartfile
            write(info_str, '(2X,A)') "Requested restart from rotated KS_eigenvector (restart_rotated .true.)"
            call localorb_info(info_str)

         case('restart_eigenvectors_periodic')
            read(inputline,*,end=88,err=99) desc_str, restart_eigenvectors_periodic
            write(info_str, '(2X,A)') "Requested restart from KS_eigenvectors (periodic)"
            call localorb_info(info_str)

         case('rlsy_symmetry')
            read(inputline,*,end=88,err=99) desc_str, desc_str
            if ( index(desc_str,'all') > 0  ) then
                use_rlsy_symmetry=.true.
            endif


         case default

             ! A number of other keywords are parsed in module control_file.f90,
             ! most notably the sc_accuracy_* keywords.
             ! This enables sharing the parser when reading a control file
             ! with a different name (control.update.in) to update
             ! certain criteria manually along the way (e.g., after an s.c.f. step).

             ex_code = run_all_parsers(desc_str, inputline)
             select case (ex_code) ! case (0) is normal exit, nothing is done
             case (88)
                 go to 88
             case (99)
                 go to 99
             case (1) ! no descriptor was found
                 if (myid.eq.0) then
                     write(use_unit,*) "Unknown descriptor ", desc_str, &
                         " in file control.in."
                 end if

                 call aims_stop_coll('', func)
             end select

         end select

      end do

!  close control.in
      ! read the default data for dielectric function
! reading pseudopot data
! use_embedding_pp will be false if there are no actual
!   pseudopotentials in geometry.in. In this case, we do
!   not read the data even if they are defined in control.in


      if (use_embedding_pp) then
        i_pp_species = 0
        do i_species = 1, n_species, 1
           if(species_pseudoized(i_species)) then
              i_pp_species = i_pp_species + 1
              if(species_nonlinear_core(i_species) ) use_nonlinear_core = .true.
              species2pp_species(i_species) = i_pp_species
              pp_species2species(i_pp_species) = i_species

              call parse_pseudopot_dimensions(pp_path_string(i_species), &
                              n_pseudofn, n_points_pseudofn)
              n_pp_fns(i_pp_species) = n_pseudofn
              n_points_pp_fn(i_pp_species) = n_points_pseudofn

           endif
        enddo


        do i_species = 1, n_species
          if(species_pseudoized(i_species)) then
             call read_pseudopotdata(i_species)
          end if
        end do
      end if

      if (linecount == 0) then
         call aims_stop_coll('Attention - empty input file control.in.', func)
      end if

      close(7)

      if (use_mbd_std) then
          if (mbd_std_parse() > 0) call aims_stop_coll()
      end if

      if (myid.eq.0) then
         write(use_unit,*)
         write (use_unit,'(2X,A)') &
         "Finished reading input file 'control.in'. Consistency checks are next."
         write(use_unit,*)
      end if

      return

   88 continue
      if (myid == 0) then
         write(use_unit,*) "Syntax error reading 'control.in'."
         write(use_unit,*) "The keyword is valid but it is missing one or more mandatory arguments."
         write(use_unit,*) "Check the input file, the manual or the source code in read_control.f90 if needed."
         write(use_unit,*) "line: '"//trim(inputline)//"'"
      end if
      call aims_stop_coll('Missing keyword argument in control.in', func)

  188 continue
      ! For the HSE06 functional, we add a special explanatory error message
      ! if the "omega" frequency parameter was missing
      if (myid.eq.0) then
        write(use_unit,'(1X,A)') "* Input error in input line (control.in):"
        write(use_unit,'(1X,A,A)') "* ",inputline
        write(use_unit,'(1X,A)') "* By definition, the HSE06 functional requires that a"
        write(use_unit,'(1X,A)') "* frequency parameter omega be specified in the following form:"
        write(use_unit,'(1X,A)') "*   xc hse06 <omega> "
        write(use_unit,'(1X,A)') "*   hse_unit <unit> "
        write(use_unit,'(1X,A)') "* where typical values of omega vary in the literature."
        write(use_unit,'(1X,A)') "* In practice, omega = 0.11 bohr^-1, approximately equal to "
        write(use_unit,'(1X,A)') "* omega = 0.208 Angstrom^-1, is sometimes found in the literature."
        write(use_unit,'(1X,A)') "* In present versions of FHI-aims, you MUST set the unit "
        write(use_unit,'(1X,A)') "* explicitly, using the 'hse_unit' tag, to avoid any confusion with "
        write(use_unit,'(1X,A)') "* the varying conventions followed by other codes and authors."
        write(use_unit,'(1X,A)') "* Allowed choices for <unit> are 'A' (Angstrom^-1) or 'B' (bohr^-1)."
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* Some comments on the choice of omega in the (early) literature:"
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* The original value of 0.15 bohr^-1 by Heyd, Scuseria and Ernzerhoff 2003"
        write(use_unit,'(1X,A)') "* was never true - see their 2006 erratum. In FHI-aims, the 'hse03' functional"
        write(use_unit,'(1X,A)') "* implements their actual choice."
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* Krukau, Vydrov, Izmaylov and Scuseria 2006 clarify the distinction between"
        write(use_unit,'(1X,A)') "* 'hse03' and 'hse06' (in addition to the Erratum mentioned above). Their"
        write(use_unit,'(1X,A)') "* conclusion is that omega=0.11 bohr^-1 is a reasonable choice."
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* Vydrov, Heyd, Krukau and Scuseria in 2006 appear to favor omega=0.25 bohr^-1,"
        write(use_unit,'(1X,A)') "* but with a mixing parameter ('hybrid_coeff') of 0.5 for the short-range exchange."
        write(use_unit,'(1X,A)') "* (The default in FHI-aims is hybrid_coeff=0.25, i.e., only a quarter of HF-like exchange.)"
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* You get the idea. As much as we would like to, we can not specify a single omega"
        write(use_unit,'(1X,A)') "* parameter for HSE06 by default - the choice is up to you. Apologies for the inconvenience."
      end if
      call aims_stop_coll('Attention - missing frequency parameter for HSE06 functional.', func)

  189 continue
      ! For the LC_wPBEh functional, we add a special explanatory error message
      ! if the "omega" frequency parameter or the hybrid_coeff was missing
      if (myid.eq.0) then
        write(use_unit,'(1X,A)') "* Input error in input line (control.in):"
        write(use_unit,'(1X,A,A)') "* ",inputline
        write(use_unit,'(1X,A)') "* By definition, the LC_wPBE functional requires that a"
        write(use_unit,'(1X,A)') "* frequency parameter omega as well as a hybrid coefficient alpha be specified in the following form:"
        write(use_unit,'(1X,A)') "*   xc LC_wPBEh <omega> <alpha>"
        write(use_unit,'(1X,A)') "*   hse_unit <unit> "
        write(use_unit,'(1X,A)') "* In present versions of FHI-aims, you MUST set the unit "
        write(use_unit,'(1X,A)') "* explicitly, using the 'hse_unit' tag, to avoid any confusion with "
        write(use_unit,'(1X,A)') "* the varying conventions followed by other codes and authors."
        write(use_unit,'(1X,A)') "* Allowed choices for <unit> are 'A' (Angstrom^-1) or 'B' (bohr^-1)."
        write(use_unit,'(1X,A)') "* "
        write(use_unit,'(1X,A)') "* Special cases: alpha = 0; omega > 0 --> LC-wPBE"
        write(use_unit,'(1X,A)') "*                alpha > 0; omega = 0 --> PBE0 (use 'xc pbe0' instead)"
        write(use_unit,'(1X,A)') "*                alpha = 0; omega = 0 --> PBE (use 'xc pbe' instead)"
      end if
      call aims_stop_coll('Attention - missing hybrid_coefficient parameter for LC_wPBEh functional.', func)

      99 continue
         if (myid == 0) then
            write(use_unit,*) "Syntax error reading 'control.in'"
            write(use_unit,*) "The keyword is valid but one or more of its mandatory arguments could not be read."
            write(use_unit,*) "Check the input file, the manual or the source code in read_control.f90 if needed."
            write(use_unit,*) "line: '"//trim(inputline)//"'"
         end if
         call aims_stop_coll('Error reading keyword argument in control.in', func)

   end subroutine readfile


   subroutine checkinput


     implicit none

!  check whether relevant quantities have been set - if not, use defaults:

!  has the xc been set?
    if (flag_xc==-1) call aims_stop(' No xc specified. What shall I do?')

    if (flag_xc==0 .and. use_hf_post) then
    ! This error says that post-processing HF with HF has no effect
    ! What is hides is that it also doesn't work currently.
    ! It'd be good to get this working as it'd be a nice sanity check for post-processing
       write(info_str,'(2X,A)') &
            'Hartree-Fock calculation requested both using the xc-functional'
       call localorb_info(info_str)
       write(info_str,'(2X,A)') &
            'and as a post-processing step. In this case post-processing has no effect.'
       call localorb_info(info_str)
       write(info_str,'(2X,A)') &
            'Turning post-processing off.'
       call localorb_info(info_str)
       use_hf_post = .false.

    end if

    if (flag_xc.eq.7) then
      ! special case, new: Check if the units for the HSE omega parameter were specified.
      if (.not.flag_hse_unit) then
         if (myid.eq.0) then
            write(use_unit,'(1X,A)') &
               "* Error. HSE06 functional requested, but 'hse_unit' not specified in control.in ."
            write(use_unit,'(1X,A)') &
               "* Explanation: The HSE06 functional can be used with an explicit screening parameter."
            write(use_unit,'(1X,A)') &
               "* The literature default value is 0.2 Angstrom^-1 or approx. 0.11 bohr^-1."
            write(use_unit,'(1X,A)') &
               "* In FHI-aims, this value had to be set by hand in the past. However, the pertinent unit"
            write(use_unit,'(1X,A)') &
               "* was bohr^-1 - in line with the reference cited in the manual, but contrary to the"
            write(use_unit,'(1X,A)') &
               "* conventional FHI-aims policy, which requires all input units to be in eV, Angstrom etc."
            write(use_unit,'(1X,A)') &
               "* To retain backwards compatibility, we have introduced a new keyword 'hse_unit', which"
            write(use_unit,'(1X,A)') &
               "* needs to be set when using the 'hse06' XC functional. Legitimate values for 'hse_unit'"
            write(use_unit,'(1X,A)') &
               "* are either 'A' (for Angstrom^-1) or 'B' (for inverse bohr radii)."
            write(use_unit,'(1X,A)') &
               "* We apologize for stopping, but rather than present our users with a silent change of units"
            write(use_unit,'(1X,A)') &
               "* and inconsistent publications, we'd much rather alert you to the question once too often."
         end if
         write(info_str, "('Missing value for hse_unit: ',A)") trim(hse_unit)
         call aims_stop_coll(info_str, func)
      else if (.not.flag_hse03) then
         ! this is the HSE06 functional - must be careful NOT to reparameterize HSE03 !
         if (hse_unit.eq.'A') then
           ! convert omega for HSE06 to bohr^-1
           hse_omega = hse_omega * bohr
           if (myid.eq.0) then
             write(use_unit,'(2X,A,F15.8,A)') &
               "HSE06 functional: Converting omega to   ", hse_omega, " bohr^-1."
           end if
         end if
         hse_omega_pbe = hse_omega
         hse_omega_hf = hse_omega
         if (myid.eq.0) then
           write(use_unit,'(2X,A,F15.8,A)') &
             "HSE06 screening parameters: omega_HF  = ", hse_omega_hf, " bohr^-1, and"
           write(use_unit,'(2X,A,F15.8,A)') &
             "                            omega_PBE = ", hse_omega_pbe, " bohr^-1."
         end if
      end if
    end if

     ! Pre-convergence steps in case of difficult functionals for SCF convergence
     ! AL: This can be meta-GGA *or* hybrid-DFT
     if (flag_xc_pre > 0) then
         flag_old_xc = flag_xc
         !hybrid_coeff_old = hybrid_coeff
         flag_xc = flag_xc_pre
         !hybrid_coeff = 0.d0
         use_hartree_fock = .false.
         use_periodic_hf = .false.
     end if

    if (use_lc_wpbeh) then
        if (.not.flag_hse_unit) then
            if (myid.eq.0) then
                write(use_unit,'(1X,A)') &
                    "* Error. LC-wPBEh functional requested, but 'hse_unit' not specified in control.in ."
                write(use_unit,'(1X,A)') &
                    "* Explanation: The LC-wPBEh functional is used with an explicit screening parameter."
                write(use_unit,'(1X,A)') &
                    "* The conventional FHI-aims policy requires all input units to be in eV, Angstrom etc."
                write(use_unit,'(1X,A)') &
                    "* To retain backwards compatibility, the keyword 'hse_unit' needs to be set when using"
                write(use_unit,'(1X,A)') &
                    "* the 'lc_wpbeh' XC functional. Legitimate values for 'hse_unit'"
                write(use_unit,'(1X,A)') &
                    "* are either 'A' (for Angstrom^-1) or 'B' (for inverse bohr radii)."
                write(use_unit,'(1X,A)') &
                    "* We apologize for stopping, but rather than present our users with a silent change of units"
                write(use_unit,'(1X,A)') &
                    "* and inconsistent publications, we'd much rather alert you to the question once too often."
            end if
            write(info_str, "('Missing value for hse_unit: ',A)") trim(hse_unit)
            call aims_stop_coll(info_str, func)
        else
            if (.not.flag_use_logsbt) then
                use_logsbt = .true.
                if (myid.eq.0) then
                    write (use_unit,'(2X,A)') &
                      "Coulomb operator: Defaulting to logarithmic spherical Bessel transform integrals (use_logsbt .true.)."
                end if
            else if (.not.use_logsbt) then
                 if (myid.eq.0) then
                     write(use_unit,'(1X,A)') &
                         "* Error: use_logsbt has to be .true. for any LC-wPBEh calculations for now. Sry"
                 end if
                 write(info_str, '(A)') "use_logsbt = .false."
                 call aims_stop_coll(info_str, func)
            end if
            if (hse_omega.eq.0.0d0) then
                 if (myid.eq.0) then
                     write(use_unit,'(1X,A)') &
                         "* Error: screening parameter = 0."
                     write(use_unit,'(1X,A)') &
                         "* Please use the already implemented functionals PBE or PBE0 for this case."
                 end if
                 write(info_str, "('omega = 0 ')")
                 call aims_stop_coll(info_str, func)
             end if
             if (hse_unit.eq.'A') then
             ! convert omega for lc-wpbeh to bohr^-1
                 hse_omega = hse_omega * bohr
                 if (myid.eq.0) then
                     write(use_unit,'(2X,A,F15.8,A)') &
                         "LC-wPBEh functional: Converting omega to   ", hse_omega, " bohr^-1."
                 end if
             end if
             hse_omega_pbe = hse_omega
             hse_omega_hf = hse_omega
             if (myid.eq.0) then
                 write(use_unit,'(2X,A,F15.8,A)') &
                     "LC-wPBEh screening parameters: omega  = ", hse_omega_hf, " bohr^-1, and"
                 write(use_unit,'(2X,A,F15.8,A)') &
                     "                               alpha = ", hybrid_coeff, "."
             end if
             if (lc_dielectric_constant.ne.1.0d0) then
                 if (myid.eq.0) then
                     write(use_unit,'(1X,A,F15.8,A)') &
                         "The dielectric constant is set to ", lc_dielectric_constant, ", thus the asymptotic decay of"
                     write(use_unit,'(1X,A,F15.8,A)') &
                         "  the exchange is now 1/(", lc_dielectric_constant, "*r)."
                 end if
             end if
        end if
        if (n_periodic /= 0 .and. use_hf_kspace .and. out_band) then
        		if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                   "* Error: You are trying to output a band structure with use_hf_kspace for a long-range corrected functional."
               write(use_unit,'(1X,A)') &
                   "* This is not implemented yet. Sorry for any inconvenience!"
           end if
           write(info_str, "('Band structure for use_hf_kspace not yet implemented for long-range corrected functionals')")
           call aims_stop_coll(info_str, func)
        end if
    end if

! Check consistency of resolution of identity method with
! correlated methods. This definition needs to be revisited later.

    if (use_lvl_fast .and. (use_corr .or. use_qpe) .and. (n_periodic.eq.0) ) then

       write(info_str,'(1X,A)') &
         '* Error. The "LVL_fast" (synonymous with "LVL") version of resolution of identity was requested'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* together with a correlated method such as RPA, GW, or MP2.'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* Currently, this combination is not supported. Please choose another version of'
       call localorb_info(info_str)
       write(info_str,'(1X,A)') &
         '* the RI_method keyword.'
       call localorb_info(info_str)
       call aims_stop('RI_method is not consistent with correlated treatment.')

    end if
    if (flag_en_shift .and. use_2d_corr) then
        write(use_unit,'(1X,A)')&
               "* Error: Currently only lapack implementation is"
        write(use_unit,'(1X,A)')&
               "* avaliable for (s)BGE2. Aborting sBGE2 part"
       call aims_stop('KS_method is not consistent with correlated treatment.')
    endif


! check consistency of use_mpi_in_place flag
    if (use_mpi_in_place) then
       l_dummy = .true.
       call check_mpi_in_place ( l_dummy )
       if (.not.l_dummy) then
          write(info_str,'(1X,A)') &
            '* MPI_IN_PLACE does not seem to work with this MPI implementation.'
          call localorb_info(info_str)
          write(info_str,'(1X,A)') &
            '* The MPI library you linked against may be outdated or otherwise problematic. Please check carefully.'
          call localorb_info(info_str)
          write(info_str,'(1X,A)') &
            '* Setting use_mpi_in_place to .false. (see manual).'
          call localorb_info(info_str)
          use_mpi_in_place = .false.
          warn_mpi_in_place = .true.
       else
          write(info_str,'(2X,A)') &
            'MPI_IN_PLACE appears to work with this MPI implementation.'
          call localorb_info(info_str)
          write(info_str,'(2X,A)') &
            "| Keeping use_mpi_in_place .true. (see manual)."
          call localorb_info(info_str)
       end if
    end if

! check for duplicate species definition

      do i_species = 1,n_species,1
        do i_counter = i_species+1, n_species,1
          if (species_name(i_species)==species_name(i_counter)) then
              write(info_str,*) '*** Error: Duplicate species definition found for ', species_name(i_species), '.'
              call localorb_info(info_str)
              call aims_stop( &
              'Please remove or rename one of the two definitions of the same species in your control.in file.')
          endif
        enddo
      enddo
!     verify integration grids
!       must check consistency of species integration grids here
      do i_species = 1, n_species, 1
         if (.not.flag_radial(i_species)) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,A,A)') &
                  "* Species ", trim(species_name(i_species)), ": ", &
                  "Missing radial integration grid."
            end if

            call aims_stop_coll('', func)
         end if
         if (.not.flag_angular(i_species)) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,A,A)') &
                  "Species ", trim(species_name(i_species)), ": ", &
                  "Missing angular integration grid."
            end if

            call aims_stop_coll('', func)
         end if
         if (.not.flag_angular_acc(i_species)) then

            if (myid.eq.0) then
               write(use_unit,'(1X,A,A,A,A)') &
                  "Species ", trim(species_name(i_species)), ": ", &
                  "Missing angular grid accuracy."
            end if

            call aims_stop_coll('', func)
         end if
         if(use_prodbas .and. n_periodic .gt. 0) then
           use_lvl = .true.
         endif
         if(use_prodbas .and. (.not. flag_RI)) then
             if(myid.eq.0) then
                write (use_unit,'(2X,A)') &
                 "Explicit treatment of the Coulomb operator is needed, but the RI method is not specified."
             endif
             if(flag_default_lvl .or. n_periodic .gt. 0 ) then
                 RI_type = RI_LVL
                 use_lvl=.true.
                 use_lvl_fast=.true.
                 if (myid.eq.0) then
                   write (use_unit,'(2X,A)') &
                  "Coulomb operator: Defaulting to 'LVL fast' resolution of identity - RI_type = LVL_fast "
                 endif
                 max_l_prodbas = 20
             else
                 RI_type = RI_V
                 if (myid.eq.0) then
                    write (use_unit,'(2X,A)') &
                    "Coulomb operator: Defaulting to 'RI-V' resolution of identity - RI_type = V "
                 end if
             endif
             sparse_o3fn = (RI_type == RI_LVL)
         endif
         if (use_prodbas) then
            if (.not.flag_prodbas_acc(i_species)) then
               if (flag_default_prodbas_acc) then
                  prodbas_acc(i_species) = default_prodbas_acc
               else if (use_lvl) then
                  prodbas_acc(i_species) = 1.0d-4
               else
                  if (species_z(i_species).le.10.d0) then
                    prodbas_acc(i_species) = 1.0d-2
                  else if (species_z(i_species).le.18.d0) then
                    prodbas_acc(i_species) = 1.0d-3
                  else
                    prodbas_acc(i_species) = 1.0d-4
                  end if
               end if
               if (myid.eq.0) then
                  write(use_unit,'(2X,4A,ES14.6,A)') &
                  & "Species ", trim(species_name(i_species)), ": ", &
                  & "Using default value for prodbas_acc = ", &
                  & prodbas_acc(i_species), "."
               end if
            end if
            if (.not.flag_max_n_prodbas(i_species)) then
               if (flag_default_max_n_prodbas) then
                  max_n_prodbas(i_species) = default_max_n_prodbas
               else if (use_lvl) then
                  max_n_prodbas(i_species) = 20
               else
                  max_n_prodbas(i_species) = 20
               end if
               ! JW: Do not output as long as not used:
               ! if (myid.eq.0) then
               !    write(use_unit,'(2X,4A,I5,A)') &
               !    & "Species ", trim(species_name(i_species)), ": ", &
               !    & "Using default value max_n_prodbas = ", &
               !    & max_n_prodbas(i_species), "."
               ! end if
            end if

            if (.not.flag_max_l_prodbas(i_species)) then
               if (flag_default_max_l_prodbas) then
                  max_l_prodbas(i_species) = default_max_l_prodbas
               else if (use_lvl) then
                  max_l_prodbas(i_species) = 20
               else
                  ! this is the case of RI-V
                  ! can use much lower max_l_prodbas for good accuracy.
                  ! 5 is good enough for most elements,
                  ! but must at least catch on-atom products of all
                  ! valence electrons ... i.e. l=3+3 for f electron systems
                  if (species_z(i_species).le.54) then
                     max_l_prodbas(i_species) = 5
                  else
                     max_l_prodbas(i_species) = 6
                  end if
               end if
               if (myid.eq.0) then
                  write(use_unit,'(2X,4A,I5,A)') &
                  & "Species ", trim(species_name(i_species)), ": ", &
                  & "Using default value max_l_prodbas = ", &
                  & max_l_prodbas(i_species), "."
               end if
            end if
         end if  ! use_prodbas
      enddo

!  make sure that any angular integration grid settings are consistent with the
!  specified Hartree potential:

      if (force_potential.eq.0) then
      ! in this case, the self-consistency loop will actually be executed
      ! only then do we need to compute the Hartree potential

         do i_species = 1, n_species, 1

   !         verify angular_limit - notice that this MUST simultaneously be done
   !         in dimensions.f
            call verify_angular_grid &
            ( angular_limit(i_species), l_hartree(i_species), &
            angular_new, flag_verify )

            if (flag_verify) then

               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A,A,A,I5,A)') &
                     "* Species ", trim(species_name(i_species)), &
                     ": Specified max. number of angular ", &
                     "integration points", &
                     " is ", angular_limit(i_species), ","
                  write(use_unit,'(1X,A,I3,A,I5,A,A,I5,A)') &
                     "* but required value for l_hartree = ", &
                     l_hartree(i_species), &
                     " is ", angular_new, ". Increasing ", &
                     "angular_limit to ", &
                     angular_new, "."
               end if

               angular_limit(i_species) = angular_new
            end if

            if (angular_acc(i_species).ne.0.d0) then
   !           also verify angular_min
               call verify_angular_grid &
               ( angular_min(i_species), l_hartree(i_species), &
                  angular_new, flag_verify )

               if (flag_verify) then

                  if (myid.eq.0) then
                     write(use_unit,'(1X,A,A,A,A,A,I5,A)') &
                           "* Species ", trim(species_name(i_species)), &
                           ": Specified min. number of angular ", &
                           "integration points", &
                           " is ", angular_min(i_species), ","
                     write(use_unit,'(1X,A,I3,A,I5,A,A,I5,A)') &
                           "* but required value for l_hartree = ", &
                           l_hartree(i_species), &
                           " is ", angular_new, ". ", &
                           "Increasing angular_min to ", &
                           angular_new, "."
                  end if

                  angular_min(i_species) = angular_new
               end if

            end if

         enddo

      end if

      ! Now, adjust specified grids, if any
      if (use_specified_grid) then
         do i_species = 1, n_species, 1
            if (specified_grid(i_species)) then
               do i_division = 1, n_ang_shells(i_species)

                  if ( n_ang_points(i_division,i_species) .lt. &
                        angular_min(i_species) ) then

                     if (myid.eq.0) then
                     write(use_unit,'(1X,A,A,A,I5,A,I5,A)') &
                        "* Species ", trim(species_name(i_species)), &
                        ": Grid shell ",i_division," contains only ", &
                        n_ang_points(i_division, i_species), &
                        " angular integration points."
                     write(use_unit,'(1X,A,I5,A)') &
                        "* Increasing to angular_min = ", &
                        angular_min(i_species), " integration points."
                     end if

                     n_ang_points(i_division, i_species) = &
                     angular_min(i_species)

                  else if ( n_ang_points(i_division,i_species) .gt. &
                        angular_limit(i_species) ) then

                     if (myid.eq.0) then
                     write(use_unit,'(1X,A,A,A,A)') &
                        "* Species ", trim(species_name(i_species)), &
                        ": Specified grid shell contains more ", &
                        "angular integration points than angular_limit."
                     write(use_unit,'(1X,A,I5,A)') &
                        "* Increasing angular_limit to = ", &
                        n_ang_points(i_division, i_species), &
                        " integration points."
                     end if

                     angular_limit(i_species) = &
                     n_ang_points(i_division, i_species)

                  end if

               enddo
            end if
         enddo
      end if

!     In the case of RI-V, we perform three-center integrals that run over the specified
!     angular grids. In the case of Au and in practice (empirically), this requires
!     Lebedev grid shells of 110 points or above for absolute total energies.
!     See FHI-aims CPC paper Blum et al 2009, p. 2183, for an explanation and references
!     to the meaning of the Lebedev grids.
!     Grid shells of 50 points or below produce large total energy errors for Au when
!     using RI-V. These errors cancel in energy differences. Nonetheless, we now increase
!     the grid shell density for RI_method 'V' to be absolutely sure that no errors happen.
!     We also include a 'stop' for heavy elements - where large differences are expected -
!     to alert users that total energies before and after the change should not be compared.
!     This 'stop' should be removed after February 1, 2015.

      if (use_prodbas.and.(RI_type.eq.RI_V)) then
         ! verify angular grids for three-center integrals
         do i_species = 1, n_species, 1

            ! maximum sum of angular momenta that should ever enter
            ! three-center integral (probably overestimated)
            l_verify = max_l_prodbas(i_species) + 2*l_shell_max(i_species)

            ! verify angular_limit - this is the largest number of grid points
            ! found on any shell for adaptive grids. If this is too small, we
            ! will not proceed. Since this is an almost impossible corner case today,
            ! we will stop the calculation instead.
            !
            ! Note that angular_limit CANNOT be increased without increasing it
            ! in dimensions.f90 first.
            !
            call verify_angular_grid &
            ( angular_limit(i_species), l_verify, &
            angular_new, flag_verify )

            if (flag_verify) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A,A,A,I5,A)') &
                     "* Species ", trim(species_name(i_species)), &
                     ": Specified max. number of angular ", &
                     "integration points", &
                     " is ", angular_limit(i_species), "."
                  write(use_unit,'(1X,A,I5,A,A,I5,A)') &
                     "* The angular momenta for this species require ", &
                      angular_new, "for RI_type 'V'. Please increase ", &
                     "your maximum angular grid to ", &
                     angular_new, " in control.in ."
               end if
               call aims_stop_coll('Please increase max. number of angular grid points.', func)
            end if

            ! Next, verify the minimum number of grid shells, angular_min.
            ! This is done both for the adaptive and the specified grid case.
            ! Since adaptive grids are a rare corner case since about 2007,
            ! there will be no special measures for adaptive grids.
            call verify_angular_grid &
            (  angular_min(i_species), l_verify, &
               angular_new, flag_verify )

            if (flag_verify) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A,A,A,A,I5,A)') &
                        "* Species ", trim(species_name(i_species)), &
                        ": Specified min. number of angular ", &
                        "integration points", &
                        " is ", angular_min(i_species), ""
                  write(use_unit,'(1X,A,I5,A,A,I5,A)') &
                     "* The angular momenta for this species require ", &
                      angular_new, " for RI_type 'V'. ", &
                     "Increasing angular_min to ", &
                     angular_new, "."
               end if
               angular_min(i_species) = angular_new
            end if

            ! For the normal case of 'specified' grids, increase
            ! the Lebedev grids on any shells that are too small.
            if (use_specified_grid) then
               if (specified_grid(i_species)) then
                  do i_division = 1, n_ang_shells(i_species)

                     if ( n_ang_points(i_division,i_species) .lt. angular_min(i_species) ) then

                        if (myid.eq.0) then
                           write(use_unit,'(1X,A,A,A,I5,A,I5,A)') &
                           "* Species ", trim(species_name(i_species)), &
                           ": Grid shell ",i_division," contains only ", &
                           n_ang_points(i_division, i_species), &
                           " angular integration points."
                           write(use_unit,'(1X,A,I5,A)') &
                           "* Increasing to angular_min = ", &
                           angular_min(i_species), " integration points."
                           write(use_unit,'(1X,A)') &
                           "* This change is related to the three-center integrals used in RI_method 'V'."
                        end if

                        n_ang_points(i_division, i_species) = &
                        angular_min(i_species)

                        warn_RI_V_grids = .true.

                        ! For the case of heavy elements, we introduce the
                        ! above-mentioned 'stop' to alert users of the total energy
                        ! change. This stop should be removed after Feb 15, 2015.
                        ! The default behavior for the code should be to adjust the
                        ! grids as needed.
                        if (species_z(i_species).gt.54) then
                           if (myid.eq.0) then
                              write(use_unit,'(1X,A)') &
                              "* You are using an electronic structure method that requires RI_method 'V'"
                              write(use_unit,'(1X,A)') &
                              "* for the resolution of identity of the Coulomb operator. (See Ren et al.,"
                              write(use_unit,'(1X,A)') &
                              "* New J. Phys 2012 for details.) This method requires three-center integrations"
                              write(use_unit,'(1X,A)') &
                              "* to be carried out on FHI-aims' normal angular and radial integration grids."
                              write(use_unit,'(1X,A)') &
                              "* (See FHI-aims CPC Blum et al. 2009, p. 2183 for grid details and references.)"
                              write(use_unit,'(1X,A)') &
                              "* One of the angular grids specified in your control.in file was not large"
                              write(use_unit,'(1X,A)') &
                              "* enough for the maximal angular momenta of orbital and product basis functions"
                              write(use_unit,'(1X,A,A,A)') &
                              "* for species ", trim(species_name(i_species)), ". "
                              write(use_unit,'(1X,A)') &
                              "* This can lead to atomic total energy errors which cancel in energy differences."
                              write(use_unit,'(1X,A)') &
                              "* Effective Sep 22, 2014, we have changed the default behavior of the code to"
                              write(use_unit,'(1X,A)') &
                              "* simply 'upgrade' the angular grids to a safe setting."
                              write(use_unit,'(1X,A)') &
                              "* However, you are also using a heavy element, where the total energy change can"
                              write(use_unit,'(1X,A)') &
                              "* lead to changes in energy differences between calculations done before and after"
                              write(use_unit,'(1X,A)') &
                              "* the change."
                              write(use_unit,'(1X,A)') &
                              "* Until February 15, 2015, the code will therefore STOP here to alert you to the change."
                              write(use_unit,'(1X,A,I5)') &
                              "* Please increase the angular grid shells in your calculation to a minimum value", &
                              angular_min(i_species)
                              write(use_unit,'(1X,A)') &
                              "* in control.in by hand before proceeding. Apologies for this inconvenience."
                           end if
                           call aims_stop_coll('Please increase max. number of angular grid points.', func)
                        end if
                     end if

                     if ( n_ang_points(i_division,i_species) .gt. angular_limit(i_species) ) then
                        ! Safety check, paranoia only.
                        ! PLEASE adjust this check if it ever triggers.
                        ! This can be done much better if the problem is more than just a corner case.
                        if (myid.eq.0) then
                        write(use_unit,'(1X,A,A,A,A)') &
                           "* Species ", trim(species_name(i_species)), &
                           ": Specified grid shell contains more ", &
                           "angular integration points than angular_limit."
                        write(use_unit,'(1X,A)') &
                           "* Please increase the maximum number of angular grid points"
                        write(use_unit,'(1X,A)') &
                           "* in your calculation."
                        end if
                        call aims_stop_coll('Please increase max. number of angular grid points.', func)
                     end if

                  enddo
               end if
            end if

         enddo
      end if

!     Verify the angular division of the grids and set a default,
!     if necessary. Not needed for the default method of creating grid batches later.

      if (.not.use_angular_division) then

         if (myid.eq.0) then
            write(use_unit,'(2X,A)') &
               "Switching off subdivision of angular shells."
         end if

      end if

   !     check whether angular grids have been requested explicitly
   !     and warn if an error will arise.

      if (flag_force_lebedev) then
         if (force_lebedev.eq.0) then
            do i_species = 1, n_species, 1
               if (angular_limit(i_species).ge.1454) then
   !             Delley grid was explicitly requested, but the largest
   !             one has only 1202 points.
               angular_limit(i_species) = 1202
               end if
            enddo
         end if
      end if

      if (.not.flag_points_in_batch) then

        if (use_gpu) then
          n_points_in_batch = 200
          write (info_str,'(2X,A,A,A,I4)') &
                "Target number of points in a grid batch ", &
                "is not set.  Since GPU acceleration is enabled, ", &
                "defaulting to ", n_points_in_batch
          call localorb_info(info_str,use_unit,'(A)')
        else
          n_points_in_batch = 100
          write (info_str,'(2X,A,A,I4)') &
                "Target number of points in a grid batch ", &
                "is not set. Defaulting to ", n_points_in_batch
          call localorb_info(info_str,use_unit,'(A)')
        end if
      else

         if (n_points_in_batch.le.0) then
            write (info_str,'(1X,A,A,I4)') &
                  "* Target number of points in a grid batch ", &
                  "set to ", n_points_in_batch
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(1X,A)') &
                  "* This value is unreasonable."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)

         else
            write (info_str,'(2X,A,A,I4)') &
                  "Target number of points in a grid batch ", &
                  "set to ", n_points_in_batch
            call localorb_info(info_str,use_unit,'(A)')
         end if
      end if

      if (.not.flag_grid_partitioning_method) then

         if (.not.use_mpi) then
            grid_partitioning_method = 5
            write (info_str,'(2X,A,A)') &
               "Method for grid partitioning is not set. ", &
               "Defaulting to 'maxmin' partitioning."
            call localorb_info(info_str,use_unit,'(A)')
         else
            ! parallel default is now parallel_hash+maxmin
            grid_partitioning_method = 6
            parallel_grid_partitioning = .true.
            use_hashed_maxmin = .true.
            write (info_str,'(2X,A,A)') &
               "Method for grid partitioning is not set. ", &
               "Defaulting to parallel hash+maxmin partitioning."
            call localorb_info(info_str,use_unit,'(A)')
         end if
      end if

      if ((grid_partitioning_method == 4).and. &
            (.not.flag_grouping_factor)) then
         write (info_str,'(2X,A,A,I6)') &
               "Factor for grid batch grouping is not set. ", &
               "Defaulting to ", grouping_factor
         call localorb_info(info_str,use_unit,'(A)')

      else if ((grid_partitioning_method == 4).and. &
               (flag_grouping_factor)) then

         if (grouping_factor > 1) then
            write (info_str,'(2X,A,I6)') &
                  "Factor for grid batch grouping set to ", &
                  grouping_factor
            call localorb_info(info_str,use_unit,'(A)')
         else
            write (info_str,'(1X,A,I6)') &
                  "* Illegal grid batch grouping factor: ", &
                  grouping_factor
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

      else if ((grid_partitioning_method /= 4).and. &
               (flag_grouping_factor)) then
         write (info_str,'(2X,A,A)') &
               "Factor for grid batch grouping doesn't apply. ", &
               "Ignoring setting."
         call localorb_info(info_str,use_unit,'(A)')
      end if

      if (flag_min_batch_size) then
         if ((grid_partitioning_method /= 1).or. &
               (parallel_grid_partitioning)) then
            write (info_str,'(2X,A,A,A)') &
                  "Minimal batch size applies ", &
                  "only to octree method. ", &
                  "Omitting setting for minimal batch size."
            call localorb_info(info_str,use_unit,'(A)')
         end if
      else
         if ((grid_partitioning_method == 1).and. &
               (.not.(parallel_grid_partitioning))) then
            write (info_str,'(2X,A,A,I6)') &
                  "Minimal batch size for octree method is not set. ", &
                  "Defaulting to ", min_batch_size
            call localorb_info(info_str,use_unit,'(A)')
         end if
      end if

      if (.not. flag_batch_size_limit) then
         write (info_str,'(2X,A,A,I6)') &
               "Batch size limit is not set. ", &
               "Defaulting to ", batch_size_hard_limit
            call localorb_info(info_str,use_unit,'(A)')
      end if



      if (prune_basis_once) then
         write (info_str,'(2X,A,A)') &
            "By default, will store active basis functions ", &
            "for each batch."
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(2X,A,A)') &
            "If in need of memory, prune_basis_once .false. ", &
            "can be used to disable this option."
         call localorb_info(info_str,use_unit,'(A)')
      end if

      ! VB: The communication type for the Hartree potential and the infrastructure for
      !     distributed_spline_storage have changed. We now set the new defaults that are
      !     fastest under any circumstances.

      if (.not.flag_communication_type) then
         ! Communication type not set ... in essence, we can do what we want
         communication_type = calc_hartree
         if (myid.eq.0) then
            write(use_unit,'(2X,A)') "communication_type for Hartree potential was not specified."
            write(use_unit,'(2X,A)') "Defaulting to calc_hartree ."
            ! RJ: I think the following hint is obsolete ... reenable if shmem is found to be much better in some cases!
            ! write(use_unit,'(2X,A)') "Note that you may want to use 'communication_type shmem' instead,"
            ! write(use_unit,'(2X,A)') "but only if you have built a code version with shared memory support."
         end if
      end if

      if (communication_type.eq.shmem_comm) then
         if (use_distributed_spline_storage) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "* Attention: "
               write(use_unit,'(1X,A)') "* Using communication_type 'shmem' for the Hartree potential."
               write(use_unit,'(1X,A)') "* This communication type works ONLY if distributed_spline_storage is .false."
               ! write(use_unit,'(1X,A)') "* Your setting seems to be .true. and must be a legacy setting."
               ! write(use_unit,'(1X,A)') "* We stop the code execution here to alert you to the fact that your 'control.in' file"
               ! write(use_unit,'(1X,A)') "* should be updated. If you REALLY wish to run with distributed_spline_storage .true. anyway,"
               ! write(use_unit,'(1X,A)') "* please simply go to subroutine read_control.f90, comment out the 'stop' command, and recompile the code."
               write(use_unit,'(1X,A)') "* Setting distributed_spline_storage to .false."
            end if
            ! call aims_stop_coll('', func)
            use_distributed_spline_storage = .false.
         end if
         flag_distributed_spline_storage = .true.
      else if (communication_type.eq.calc_hartree) then

         ! If the default of use_distributed_spline_storage is .false. and the user doesn't request it ...
         if (.not.flag_distributed_spline_storage .and. .not.use_distributed_spline_storage .and. .not.use_mpi) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "* Attention: "
               write(use_unit,'(1X,A)') "* You may want to set distributed_spline_storage to .true. in order to save memory."
               write(use_unit,'(1X,A)') "* This doesn't harm performance and is perfectly ok unless:"
               write(use_unit,'(1X,2A)') &
                  "* - you get error messages in postprocessing ", &
                  "(about missing atoms in get_rho_multipole_spl)"
               write(use_unit,'(1X,A)') "* - you need full spline storage in hirshfeld_analysis"
            endif
         end if

      end if

!     Verify settings for self-consistent part

      if (force_potential.eq.0) then

         if (.not.flag_mixer) then
   !         use default.

            if (myid.eq.0) then
               write (use_unit,'(2X,A)') &
                     "Defaulting to Pulay charge density mixer."
            end if

            mixer = MIX_PULAY
         end if

         if (mixer.eq.MIX_PULAY) then
            if (.not.flag_pulay_iter) then
               n_max_pulay = 8

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Pulay mixer: Number of relevant ", &
                        "iterations not set."
                  write(use_unit,'(2X,A,I4,A)') &
                        "Defaulting to ", n_max_pulay, &
                        " iterations."
               end if

            end if
            if (.not.flag_ini_linear_mixing) then
               ini_linear_mixing = 0

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "Pulay mixer: Number of initial linear mixing ", &
                     "iterations not set."
                  write(use_unit,'(2X,A,I4,A)') &
                     "Defaulting to ", ini_linear_mixing, &
                     " iterations."
               end if

            end if
         end if

         if (mixer.eq.MIX_BROYDEN) then
            if (.not.flag_broyden_iter) then
               n_max_broyden = 8

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Broyden mixer: Number of relevant ", &
                        "iterations not set."
                  write(use_unit,'(2X,A,I4,A)') &
                        "Defaulting to ", n_max_broyden, &
                        " iterations."
               end if

            end if
            if (.not.flag_relative_fp_charge_mix) then
               relative_fp_charge_mix = 0.05

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Broyden mixer: Relative fixed point ", &
                        "under-relaxation not set."
                  write(use_unit,'(2X,A,F10.4,A)') &
                        "Defaulting to ", relative_fp_charge_mix, &
                        " iterations."
               end if

            end if
            if (.not.flag_pulay_iter) then
               n_max_pulay = 8

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Broyden mixer: Number of relevant ", &
                        "iterations not set."
                  write(use_unit,'(2X,A,I4,A)') &
                        "Defaulting to ", n_max_pulay, &
                        " iterations."
               end if

            end if
            if (.not.flag_ini_linear_mixing) then
               ini_linear_mixing = 0

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                     "Broyden mixer: Number of initial linear mixing ", &
                     "iterations not set."
                  write(use_unit,'(2X,A,I4,A)') &
                     "Defaulting to ", ini_linear_mixing, &
                     " iterations."
               end if

            end if
         end if

         ! If we use Fock matrix mixing, then the density matrix must be mixed.
         use_density_matrix_hf = .true.
         ! JW (2010-08-24):
         ! So, use_density_matrix_hf is only .false. if (force_potential /= 0)??

         if (flag_mixer_threshold_charge) then
            if ( (spin_treatment.eq.1) .and. &
                  (.not.flag_mixer_threshold_spin) ) then

               max_rho_change(2) = max_rho_change(1)

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Mixer threshold: Charge mixer threshold defined, ", &
                  "but spin density mixer threshold not defined."
                     write(use_unit,'(2X,A,E14.6,A)') &
                     "Defaulting to ", max_rho_change(2), &
                     " threshold for spin density mixing."

               end if
            end if
         end if

         if (flag_mixer_threshold_spin) then
            if ( (.not.flag_mixer_threshold_charge) ) then

               max_rho_change(1) = 0.d0

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                  "Mixer threshold: Spin mixer threshold defined, ", &
                  "but charge density mixer threshold not defined."
                     write(use_unit,'(2X,A)') &
                     "Thresholding spin density mixing only."
               end if
            end if
         end if

         if (.not.flag_hartree_worksize) then

            hartree_worksize = 200
            if (myid.eq.0) then
               write(use_unit,'(2X,A,A)') &
               "Work space size for distributed Hartree potential ", &
               "not set."
               write(use_unit,'(2X,A,E14.6,A)') &
               "Defaulting to ", hartree_worksize, &
               " MB."
            end if

         else
            ! check the chosen value to make sure it will not cause problems

            if (hartree_worksize.lt.50) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A,A)') &
                  "* Selected work space size for distributed Hartree ", &
                  "potential is unusually small."
                  write(use_unit,'(1X,A,A)') &
                  "* This could impact performance. Please remove this ", &
                  "warning from the code if you need this value."
               end if
               call aims_stop_coll('', func)
            end if

         end if

      ! Set reasonable defaults for charge density mixing.
      ! This is a pain.

      ! Set mixing parameter for linear mixing
         if (mixer.eq.MIX_LINEAR) then

            if (flag_mix_param.and.flag_linear_mix_param) then

               ! use default
               linear_mix_param(1) = charge_mix_param(1)

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A,A)') &
                        "Mixing parameters for charge density mixing ", &
                        "and for initial linear mixing", &
                        " have both been set."
                  write(use_unit,'(2X,A,A)') &
                        "You are only doing linear mixing - ", &
                        "can only use one value."
                  write(use_unit,'(2X,A,F10.4,A)') &
                        "Defaulting to charge_mix_param = ", &
                        linear_mix_param(1), "."
               end if

            else if (flag_mix_param) then

               linear_mix_param(1) = charge_mix_param(1)

            else if (flag_linear_mix_param) then

               if (myid.eq.0) then
               write(use_unit,'(2X,A,A)') &
                     "Linear mixing: charge_mix_param not set, ", &
                     "but ini_linear_mix_param set."
               write(use_unit,'(2X,A,F10.4,A)') &
                     "Defaulting to use ini_linear_mix_param = ", &
                     linear_mix_param(1), " throughout."
               end if

            else
               ! use default

               linear_mix_param(1) = 0.1

               if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                     "Linear mixing: Mixing parameter not set."
               write(use_unit,'(2X,A,F10.4,A)') &
                     "Defaulting to ", &
                     linear_mix_param(1), "."
               end if

            end if

            ! if spin polarization required, check the same for spin density mixing
            if (n_spin.gt.1) then
               if (flag_linear_spin_mix_param.and.flag_spin_mix_param) then

                  linear_mix_param(2) = charge_mix_param(2)
                  if (myid.eq.0) then
                     write(use_unit,'(2X,A)') &
                     "Linear spin density mixing: Two parameters are set."
                     write(use_unit,'(2X,A,F10.5,A)') &
                     "Defaulting to spin_mix_param: ", charge_mix_param(2), &
                     "."
                  end if

               else if (flag_spin_mix_param) then

                  linear_mix_param(2) = charge_mix_param(2)

               else if (flag_linear_spin_mix_param) then

                  charge_mix_param(2) = linear_mix_param(2)

               else
               ! use default

                  linear_mix_param(2) = linear_mix_param(1)

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A)') &
                        "Linear mixing: Spin mixing parameter not set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Defaulting to ", &
                        linear_mix_param(2), "."
                  end if

               end if
            end if

         end if

         ! Set parameters for Pulay mixing and possible ini. linear mixing
         if (mixer.eq.MIX_PULAY) then

            if (.not.flag_mix_param) then
               ! use default

               if (.not.adjust_scf) then
                 charge_mix_param(1) = 0.2
               else
                 !more cautious first step
                 charge_mix_param(1) = 0.05
               end if

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for charge density mixing ", &
                        "has not been set."
                  write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: charge_mix_param = ", &
                        charge_mix_param(1), "."
               end if

               charge_mix_param_set = .false.

               if (adjust_scf) then
                  if (adjust_scf_always) then
                     if (myid.eq.0) then
                        write (use_unit,'(2X,A,I5,A)') &
                        "The mixing parameter will be adjusted in iteration number ", &
                        adjust_scf_iteration, " of each s.c.f. cycle (i.e., each geometry step)."
                     end if
                  else
                     if (myid.eq.0) then
                        write (use_unit,'(2X,A,I5,A)') &
                        "The mixing parameter will be adjusted in iteration number ", &
                        adjust_scf_iteration, " of the first full s.c.f. cycle only."
                     end if
                  end if
               end if

            end if

            if (ini_linear_mixing.gt.0) then

               if (.not.flag_linear_mix_param) then
               ! use default
                  linear_mix_param(1) = 0.1

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for initial linear ", &
                        "charge density mixing has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: ", &
                        linear_mix_param(1), "."
               end if

               end if

            end if

            if (n_spin.gt.1) then

               if (.not.flag_spin_mix_param) then

                  charge_mix_param(2) = charge_mix_param(1)

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for spin density mixing ", &
                        "has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using charge_mix_param as default: ", &
                        charge_mix_param(2), "."
                  end if

               end if

               if ( (ini_linear_mixing.gt.0) .and. &
                  (.not.flag_linear_spin_mix_param) ) then

                  linear_mix_param(2) = linear_mix_param(1)

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for initial linear ", &
                        "spin density mixing has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: ", &
                        linear_mix_param(2), "."
                  end if

               end if

            end if

         end if

         ! Set parameters for Broyden mixing and possible ini. linear mixing
         if (mixer.eq.MIX_BROYDEN) then

            if (.not.flag_mix_param) then
               ! use default
               charge_mix_param(1) = 0.2

               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for charge density mixing ", &
                        "has not been set."
                  write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: charge_mix_param = ", &
                        charge_mix_param(1), "."
               end if

            end if

            if (ini_linear_mixing.gt.0) then

               if (.not.flag_linear_mix_param) then
               ! use default
                  linear_mix_param(1) = 0.1

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for initial linear ", &
                        "charge density mixing has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: ", &
                        linear_mix_param(1), "."
               end if

               end if

            end if

            if (n_spin.gt.1) then

               if (.not.flag_spin_mix_param) then

                  charge_mix_param(2) = charge_mix_param(1)

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for spin density mixing ", &
                        "has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using charge_mix_param as default: ", &
                        charge_mix_param(2), "."
                  end if

               end if

               if ( (ini_linear_mixing.gt.0) .and. &
                  (.not.flag_linear_spin_mix_param) ) then

                  linear_mix_param(2) = linear_mix_param(1)

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Mixing parameter for initial linear ", &
                        "spin density mixing has not been set."
                     write(use_unit,'(2X,A,F10.4,A)') &
                        "Using default: ", &
                        linear_mix_param(2), "."
                  end if

               end if

            end if

         end if

         if (mixer.eq.MIX_PULAY) then

            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                  "Algorithm-dependent basis array size parameters:"
               write(use_unit,'(2X,A,I8)') &
                  "| n_max_pulay                         : ", &
                  n_max_pulay
            end if

         end if

         if (mixer.eq.MIX_BROYDEN) then

            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                  "Algorithm-dependent basis array size parameters:"
               write(use_unit,'(2X,A,I8)') &
                  "| n_max_broyden                         : ", &
                  n_max_broyden
            end if
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                  "Algorithm-dependent relative under-relaxation:"
               write(use_unit,'(2X,A,E14.6)') &
                  "| relative_fp_charge_mix                         : ", &
                  relative_fp_charge_mix
            end if

         end if

         if (use_constraint) then
            if (.not.flag_mixer_constraint) then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A)') &
                     "No mixer set for ", &
                     "spin constrained DFT."
                  write (use_unit,'(2X,A)') &
                     "Defaulting to BFGS-minimization routine ", &
                     "in spin constrained DFT."
               end if
               mixer_constraint = 2
            end if
            if (mixer_constraint.eq.1) then
               if (.not.flag_ini_linear_mixing_constraint) then
                  ini_linear_mixing_constraint = 0

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                        "Spin constrained DFT: Pulay mixer: Number of ", &
                        " initial linear mixing iterations not set."
                     write(use_unit,'(2X,A,I4,A)') &
                        "Defaulting to ", ini_linear_mixing_constraint, &
                        " iterations."
                  end if

               end if
               if (.not.flag_pulay_iter_constraint) then
                  n_max_pulay_constraint = 8

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A)') &
                           "Spin constrained DFT: Pulay mixer: Number of ", &
                           "relevant iterations not set."
                     write(use_unit,'(2X,A,I4,A)') &
                           "Defaulting to ", n_max_pulay_constraint, &
                           " iterations for no good reason."
                  end if

               end if

            end if

         end if

!SK: This part has moved down the file, since defaults of sc_accuracy_rho depend
!    now on if forces or stress is calucalated. Thus, it has to be moved to the 
!    last change of use_forces.

!         if (.not.flag_acc_rho) then
!            ! Decide on a default value for sc_accuracy_rho that reflects
!            ! empirical experience for good convergence.
!            if (n_atoms.le.6) then
!               sc_accuracy_rho = 1.0d-6
!            else if (use_forces .or. use_analytical_stress) then
!               if (n_atoms.le.6000) then
!                  sc_accuracy_rho = 1.0d-6 * sqrt(dble(n_atoms)/6.d0)
!               else ! safety check because sc_accuracy_rho really was not tested
!                    ! for situations above 6000 atoms
!                  sc_accuracy_rho = 1.0d-6 * sqrt(1000.d0)
!               end if
!            else
!               if (n_atoms.le.6000) then
!                  sc_accuracy_rho = 1.0d-5 * dble(n_atoms)/60.d0
!               else ! safety check because sc_accuracy_rho really was not tested
!                    ! for situations above 6000 atoms
!                  sc_accuracy_rho = 1.0d-3
!               end if
!            end if
!
!            write(info_str,'(1X,A,A)') &
!                  "* Notice: The s.c.f. convergence criterion sc_accuracy_rho was", &
!                  " not provided."
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A)') &
!                  "* We used to stop in this case, and ask the user to provide reasonable"
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A)') &
!                  "* scf convergence criteria. However, this led some users to employ criteria"
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A)') &
!                  "* that led to extremely long run times, e.g., for simple relaxation."
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A)') &
!                  "* We now preset a default value for sc_accuracy_rho if it is not set."
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A)') &
!                  "* You may still wish to check if this setting is too tight or too loose for your needs."
!            call localorb_info(info_str,use_unit,'(A)')
!            write(info_str,'(1X,A,E13.6,A)') &
!                  "* Based on n_atoms, FHI-aims chose sc_accuracy_rho = ", sc_accuracy_rho, " ."
!            call localorb_info(info_str,use_unit,'(A)')
!
!            ! Check the accuracy of the forces and shout or stop if a waste of time can happen
!            if ( flag_acc_forces ) then
!               if (sc_accuracy_forces .gt. 0 .and. &
!                   sc_accuracy_forces .lt. 1.0d-6) then
!                  write(info_str,'(2X,A,1X,A,1X,E10.3,A)') &
!                      "*** You did choose a very tight ", &
!                      "force convergence criterion (sc_accuracy_forces) of ", &
!                      sc_accuracy_forces, " ,"
!                  call localorb_info(info_str,use_unit,'(A)')
!                  write(info_str,'(2X,A,1X,A,1X,A)') &
!                      "*** but did not specify sc_accuracy_rho. This", &
!                      "would most certainly lead to unnecessary many force"
!                  call localorb_info(info_str,use_unit,'(A)')
!                  write(info_str,'(2X,A)') &
!                      "*** evaluations, slowing the computation down (see manual)."
!                  call localorb_info(info_str,use_unit,'(A)')
!                  write (info_str,'(2X,A,A)') &
!                      "*** Please carefully choose an sc_accuracy_rho criterion ", &
!                      "and modify your input file control.in accordingly."
!                  call localorb_info(info_str,use_unit,'(A)')
!                  call localorb_info('',use_unit,'(A)')
!                  call aims_stop_coll('', func)
!               else if (sc_accuracy_forces .gt. 0 .and. &
!                   sc_accuracy_forces .lt. 1.0d-4) then
!                  write(info_str,'(2X,A,1X,A,1X,E10.3,A)') &
!                      "*** WARNING: You did choose a rather tight", &
!                      "force convergence criterion (sc_accuracy_forces) of", &
!                      sc_accuracy_forces, " ,"
!                  call localorb_info(info_str,use_unit,'(A)')
!                  write(info_str,'(2X,A,1X,A)') &
!                      "*** but did not specify sc_accuracy_rho. Please", &
!                      "be aware that this can lead to unnecessary many force"
!                  call localorb_info(info_str,use_unit,'(A)')
!                  write(info_str,'(2X,A)') &
!                      "*** evaluations, slowing the computation down (see manual)."
!                  call localorb_info(info_str,use_unit,'(A)')
!               end if
!            end if
!
!
!         end if
         if (.not.flag_acc_eev) then
            sc_accuracy_eev = 0.d0
         end if
         if (.not.flag_acc_etot) then
            sc_accuracy_etot = 0.d0
         end if

         if (.not.flag_sc_max) then
   !         use default
            sc_iter_limit = 1000

            write(info_str,'(2X,A,A)') &
                  "Maximum number of self-consistency ", &
                  "iterations not provided."
            call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,1X,I5,1X,A)') &
                  "Presetting", sc_iter_limit, "iterations."
            call localorb_info(info_str,use_unit,'(A)')
         end if

         if (.not.flag_sc_init_iter) then
            ! default value is set here
            sc_init_iter = 1001 ! increased from 40 to 1001 on Feb 10, 2018, after observing a
                                ! failure due to sc_init_iter that affects large slab systems.
                                ! This should be revisited with more precise logic to catch this behavior.
            scf_restarted = .false.
            write(info_str,'(2X,A,1X,I8,1X,A)') &
                  "Presetting ", sc_init_iter, "iterations before the initial mixing cycle"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
                  "is restarted anyway using the sc_init_iter criterion / keyword."
            call localorb_info(info_str,use_unit,'(A)')
         end if

         if (.not.flag_sc_init_factor) then
            ! default value is set here
            sc_init_factor = 1.d0
            write(info_str,'(2X,A,1X,F10.3,1X,A)') &
                  "Presetting a factor", sc_init_factor, "between actual scf density residual"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
                  "and density convergence criterion sc_accuracy_rho below which sc_init_iter"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
                  "takes no effect."
            call localorb_info(info_str,use_unit,'(A)')

         end if

      else
      ! force_potential .ne. 0 - see way above for corresponding if!

         if (.not.flag_sc_max) then
!          need only scf initialization anyway
            sc_iter_limit = 0
         end if

      end if

      if ( (force_potential.eq.0) .and. &
         (.not.flag_acc_rho)    .and. &
         (.not.flag_acc_eev)    .and. &
         (.not.flag_acc_etot)   .and. &
         (.not.flag_acc_potjump) .and. &
         (.not.(use_forces .and. sc_accuracy_forces > 0.d0))) then

            write(info_str,'(1X,A)') &
                  "* No s.c.f. convergence criteria (sc_accuracy_*) were provided in control.in."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
                  "* The run will proceed with a reasonable default guess, but please check whether."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
                  "* the s.c.f. cycles seem to take too much or too little time."
            call localorb_info(info_str,use_unit,'(A)')

      end if

      if ( .not.( flag_force .or. flag_acc_forces .or. flag_relax_mode &
         .or. flag_molecular_dynamics ) ) then
         if (myid .eq. 0) then
            write(use_unit,'(2X,A,A)') &
            "Calculation of forces was not defined in control.in. ", &
            "No forces will be calculated."
         end if

         compute_forces = .false.

      end if

         ! Relaxation and molecular dynamics make no sense in the presence of
         ! "empty" sites (unless a full BSSE correction also of forces is
         ! performed). Disabling this for now.
      if ((flag_relax_mode .and. (relax_mode.ne.RELAX_OFF)) .or. &
      &   flag_molecular_dynamics) then
         if (n_empty_atoms.gt.0) then
            write (info_str,'(1X,A,A)') &
            "* 'empty' sites in the structure may lead to physically ", &
            "meaningless relaxation or molecular dynamics."
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(1X,A,A)') &
            "* If you know what you are doing, disable this warning ", &
            "and recompile."
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(1X,A,A)') &
            "* Else, please perform this run either without empty sites, ", &
            "or without relaxation / MD."
            call localorb_info(info_str,use_unit,'(A)')
 !DB 092612 .. this is just a little hack .. remove that afterwards
  !          call aims_stop_coll('', func)
         end if
      end if

         ! Check all other items associated with relaxation
      if (.not.flag_relax_mode) then

         relax_mode = RELAX_OFF
         relax_geo = RELAX_OFF
         write (info_str,'(2X,A,A)') &
         "Geometry relaxation not requested:", &
         " no relaxation will be performed."
         call localorb_info(info_str,use_unit,'(A)')

      else
         if ( (relax_mode.ne.RELAX_OFF) .and. (.not.compute_forces) ) then

            call localorb_info('',use_unit,'(A)')
            write (info_str,'(1X,A,A)') &
            "* Note: You explicitly requested no force computation, ", &
            "but explicitly asked for geometry relaxation."
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(1X,A,A)') &
            "* This combination of requests is outrightly ", &
            "contradictory - please modify your input file control.in."
            call localorb_info(info_str,use_unit,'(A)')
            call localorb_info('',use_unit,'(A)')
            call aims_stop_coll('', func)

         end if
         if (relax_mode.ne.RELAX_OFF) then
            ! inform about choices that were not made wrt BFGS relaxation,
            ! where defaults are being used
            if (.not.flag_max_relaxation_steps) then
               write(info_str,'(2A,I8)')'  No maximum number', &
               ' of relaxation steps, defaulting to ',max_relaxation_steps
               call localorb_info(info_str,use_unit,'(A)')
            end if
            if (.not.flag_init_hess) then
               if ((n_hess_blocks > 0)        .or. &
                  (n_hess_blocks_lv > 0)      .or. &
                  (n_hess_blocks_lv_atom > 0) .or. &
                  hess_in_file) then
                  init_hess_type = HESS_EXPL
               else
                  ! Default Hessian corresponds to "init_hess Lindh 2.".  This
                  ! is different from "init_hess Lindh" in order to avoid
                  ! problems with systems containing many loosely bound
                  ! solvation molecules.
                  !
                  ! For lattice vectors, the additional columns/rows of the
                  ! Hessian are diagonal with a default of their own.
                  init_hess_type = HESS_LINDH
                  init_hess_lindh_diag = 2.d0 / (hartree/bohr**2)
                  init_hess_lv_diag = 25.d0 / (hartree/bohr**2)
               end if
               select case (init_hess_type)
               case(HESS_EXPL)
                  write(info_str,'(2X,2A,E14.6)') &
                  & 'Using explicit Hessian from geometry.in.'
                  call localorb_info(info_str)
               case(HESS_DIAG)
                  write(info_str,'(2X,A,E14.6,A)') &
                  & 'Default initial Hessian is ', &
                  & init_hess_diag*(hartree/bohr**2), &
                  & ' eV/A^2 times unity.'
                  call localorb_info(info_str)
               case(HESS_LINDH)
                  write(info_str,'(2X,A,F6.2,A,E14.6,A)') &
                  & 'Default initial Hessian is Lindh matrix (thres =', &
                  init_hess_lindh_thres, ') plus', &
                  & init_hess_lindh_diag*(hartree/bohr**2), &
                  & ' eV/A^2 times unity.'
                  call localorb_info(info_str)
               case(HESS_LATTICE)
                 ! In case of trouble, leave parameters on default.
                 write(info_str,'(2X,A,F6.2,A,E14.6,3A)') &
                 & 'Default initial Hessian is Lindh matrix (thres =', &
                 init_hess_lindh_thres, ') plus', &
                 & init_hess_lindh_diag*(hartree/bohr**2), &
                 & ' eV/A^2 times unity for atomic positions.', &
                 & ' For the lattice, the reciprocal lattice matrix',&
                 & ' is used.'
                 call localorb_info(info_str)
               case default
                  call aims_stop('Assertion: Invalid Hessian type', func)
               end select
            end if
            if (.not.flag_energy_tolerance) then
               write(info_str,'(2A,E14.6)')'  No maximum energy', &
               & ' tolerance for TRM/BFGS moves, defaulting to ', &
               & energy_tolerance*hartree
               call localorb_info(info_str,use_unit,'(A)')
            end if
            if (.not.flag_aggregated_energy_tolerance) then
               write(info_str,'(2A,E14.6)')'  Maximum energy', &
               & ' tolerance by which TRM/BFGS trajectory may increase over multiple steps: ', &
               & aggregated_energy_tolerance*hartree
               call localorb_info(info_str,use_unit,'(A)')
            end if
            if ((relax_mode == RELAX_TRM) .or. (relax_mode == RELAX_LTRM)) then
               if (.not.flag_harmonic_length_scale) then
                  write(info_str,'(2A,E14.6,A)')'  No harmonic length scale.', &
                  & ' Defaulting to ', harmonic_length_scale*bohr, ' A.'
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               if (.not.flag_max_atomic_move) then
                  write(info_str,'(2A,E14.6,A)') &
                  & '  No trust radius initializer.', &
                  & ' Defaulting to ', max_atomic_move*bohr, ' A.'
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            else   ! BFGS
               if (.not.flag_max_atomic_move) then
                  write(info_str,'(2A,E14.6)')'  No maximum atomic', &
                  ' movement during BFGS, defaulting to ', max_atomic_move*bohr
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               if (.not.flag_min_line_step) then
                  write(info_str,'(2A,E14.6)')'  No explicit', &
                  ' minimum BFGS line step, defaulting to ', min_line_step
                  call localorb_info(info_str,use_unit,'(A)')
               end if
               if (.not.flag_line_step_reduce) then
                  write(info_str,'(2X,2A)')'Defaulting to automatic ', &
                  'line searching during relaxation.'
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            end if
         end if
      end if

      if (relax_mode .eq. RELAX_OFF) then
         ! someone might have set relevant flags for BFGS relaxation, but then
         ! forgotten to actually request the relaxation itself, inform that the
         ! BFGS parameters are ignored if that is the case.
         if (flag_max_relaxation_steps) then
         write(info_str,'(2A)') 'Requested max number of relaxation ', &
               'steps, but no relevant relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_init_hess) then
            write(info_str,'(2A)') 'Requested initial Hessian,', &
            ' but no relevant relaxation.  Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_max_atomic_move) then
         write(info_str,'(2A)') 'Requested max value for single atom ', &
            'move, but no relevant relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_energy_tolerance) then
         write(info_str,'(A,A)')'Requested energy tolerance for BFGS,', &
            ' but no relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_aggregated_energy_tolerance) then
         write(info_str,'(A,A)')'Requested aggregated energy tolerance for BFGS,', &
            ' but no relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_min_line_step) then
         write(info_str,'(A,A)')'Requested minimum line step for BFGS,', &
            ' but no relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (flag_line_step_reduce) then
         write(info_str,'(A,A)')'Requested line step reduction for ', &
            'BFGS, but no relaxation. Ignoring input. '
            call localorb_info(info_str,use_unit,'(A)')
         end if
      end if

      ! FK: Here, we check whether the user settings are compatible with the unit cell relaxation.
      !     We determine which stress (numerical or analytical) is used, too.
      if (relax_unit_cell .ne. 0) then

         ! Check for structure optimizer
         if (relax_mode == RELAX_OFF) then
            write (info_str,'(2X,A)') &
               "*** Unit cell relaxation: Need to specify structure optimizer for unit cell relaxation!"
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A)') &
               "***                       Please specify in control.in the following: relax_geometry <type> <tolerance>"
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! If unit cell relaxation was requested, determine stress_for_relaxation,
         ! either automatically or use user choice (if it is reasonable)

         ! The following line sets, for the first time, the stress_for_relaxation default
         ! if it was not user specified in control.in.
         ! As written, it alters the global initialization of this variable in
         ! runtime_choices.f90.
         ! This snippet used to be different for the non-load_balancing vs. load_balancing cases.
         ! However, at the time of writing it is my belief that the load_balancing case has
         ! been fixed to give the correct stress tensor also then. In my tests, this is the case.
         ! Please keep an eye on this, however, and please let me know if there are cases in
         ! which the analytical stress tensor fails due to the use of load_balancing.
         ! - VB
         if (.not. numerical_stress_needed_for_relaxation) then
            stress_for_relaxation = RELAX_ANALYTICAL
         end if

         ! Check for vdW
         if (.not.(use_vdw_method.or.use_vdw_correction.or.use_mbd_old) ) then
            if (.not. numerical_stress_needed_for_relaxation) then
               stress_for_relaxation = RELAX_ANALYTICAL
            end if
         else
            if ( flag_stress_for_relaxation_user .and. (stress_for_relaxation_user.eq.RELAX_ANALYTICAL) ) then
               write(info_str,'(2X,3A)') &
                  "*** Unit cell relaxation: Analytical stress was ", &
                  "specified for relaxation but this is not possible with ", &
                  "vdW corrections besides Hirschfeld."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if
            stress_for_relaxation = RELAX_NUMERICAL
            numerical_stress_needed_for_relaxation = .true.
         end if

         ! Check for used functional
         if (flag_xc.eq. 0 .or. & ! hf
             flag_xc.eq. 1 .or. & ! pbe0
             flag_xc.eq. 3 .or. & ! pz-lda
             flag_xc.eq. 6 .or. & ! pbe
             flag_xc.eq. 7 .or. & ! hse
             flag_xc.eq. 8 .or. & ! pw-lda
             flag_xc.eq. 9 .or. & ! blyp
             flag_xc.eq.10 .or. & ! b3lyp
             flag_xc.eq.11 .or. & ! rpbe
             flag_xc.eq.12 .or. & ! revpbe
             flag_xc.eq.13 .or. & ! pbesol0
             flag_xc.eq.15 .or. & ! vwn
             flag_xc.eq.16 .or. & ! vwn-gauss
             flag_xc.eq.17 .or. & ! pbesol
             flag_xc.eq.18 .or. & ! r48pbe
             flag_xc.eq.20 .or. & ! am05
             flag_xc.eq.21 .or. & ! pbeint
             flag_xc.eq.22 .or. & ! xpbe
             flag_xc.eq.25 .or. & ! m06-l
             flag_xc.eq.26 .or. & ! m06
             flag_xc.eq.27 .or. & ! m06-2x
             flag_xc.eq.28 .or. & ! m06-hf
             flag_xc.eq.29 .or. & ! b1lyp
             flag_xc.eq.51 .or. & ! tpss
             flag_xc.eq.52 .or. & ! revtpss
             flag_xc.eq.53 .or. & ! tpssloc
             flag_xc.eq.55 .or. & ! m08-hx
             flag_xc.eq.56 .or. & ! m08-so
             flag_xc.eq.58 .or. & ! m11-l
             flag_xc.eq.59 .or. & ! scan
             flag_xc.le.-10)    & ! libxc
         then
            if (.not. numerical_stress_needed_for_relaxation) then
               stress_for_relaxation = RELAX_ANALYTICAL
            end if
         else
            if ( flag_stress_for_relaxation_user .and. (stress_for_relaxation_user.eq.RELAX_ANALYTICAL) ) then
               write(info_str,'(2X,3A)') &
                  "*** Unit cell relaxation: Analytical stress was ", &
                  "specified for relaxation but this is only possible with ", &
                  "LDA, GGA, Meta-GGA and hybrid functionals."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            end if
            stress_for_relaxation = RELAX_NUMERICAL
            numerical_stress_needed_for_relaxation = .true.
         end if

         ! Check for relativistic settings
         if (flag_rel .eq. REL_zora) then
            write(info_str,'(1X,A)') '* Error: No forces available for "relativistic zora".'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') '* For relaxations, please use "relativistic atomic_zora scalar".'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') '* Please correct your control.in file before continuing.'
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! We checked above whether the user choice was reasonable or not.
         ! Therefore, set the global variable stress_for_relaxation according to the user choice.
         if (flag_stress_for_relaxation_user) then
            stress_for_relaxation = stress_for_relaxation_user
         end if

         ! Print which unit cell relaxation is used
         if (relax_unit_cell .eq. 1) then
             write (info_str,'(2X,A,A)') &
                "Unit cell relaxation: Unit cell will be relaxed fully."
             call localorb_info(info_str,use_unit,'(A)')
         else if (relax_unit_cell .eq. 2) then
             write (info_str,'(2X,A,A)') &
                "Unit cell relaxation: Restricted relaxation with fixed angles between lattice vectors."
             call localorb_info(info_str,use_unit,'(A)')
         end if

         ! Print which stress is used for unit cell relaxation and set variables accorrdingly
         if (stress_for_relaxation .eq. RELAX_NUMERICAL) then
            write (info_str,'(2X,A,A)') &
               "Unit cell relaxation: Numerical stress will be used."
            call localorb_info(info_str,use_unit,'(A)')
            if (.not.use_numerical_stress) then
               use_numerical_stress = .true.
            end if
         else if (stress_for_relaxation .eq. RELAX_ANALYTICAL) then
            write (info_str,'(2X,A,A)') &
               "Unit cell relaxation: Analytical stress will be used."
            call localorb_info(info_str,use_unit,'(A)')
            if (.not.use_analytical_stress) then
               use_basis_gradients = .true.
               use_analytical_stress = .true.
            end if
         end if
      end if !relax_unit_cell .ne. 0

      ! Check use_distributed_hessian
      if (relax_mode /= RELAX_OFF) then
         if (.not. flag_distributed_hessian) then
             use_distributed_hessian = .true.
         end if

         ! Must have MPI
         if (.not. use_mpi) then
            use_distributed_hessian = .false.
         end if

         ! Must have BLACS
         call BLACS_Pinfo(index1,index2)

         if (index2 <= 0) then
            use_distributed_hessian = .false.
         end if

         ! Textbook BFGS not implemented
         if (relax_mode == RELAX_BFGS_TB) then
            use_distributed_hessian = .false.
         end if

         ! HESS_LATTICE not implemented
         if (init_hess_type == HESS_LATTICE) then
            use_distributed_hessian = .false.
         end if
      else ! No relaxation
         use_distributed_hessian = .false.
      end if

      if (flag_acc_forces .and. (.not.compute_forces)) then
         ! if force computation was explicitly prohibited and no
         ! relaxation was requested, prevent force computation
         ! even if a convergence accuracy for the forces was
         ! explicitly specified.

         use_forces = .false.

      end if

      if ((.not. flag_delta_numerical_stress) .and. (use_numerical_stress)) then
         delta_numerical_stress = default_delta_numerical_stress
         if (myid.eq.0) then
         write (use_unit,'(2X,A,A,E14.6,A)') &
               "No scaling factor delta for numerical stress given. ", &
               "Defaulting to ", delta_numerical_stress, "."
         end if
      end if

      ! CC: Switch on analytical stress for HF calc.
      if (compute_heat_flux) then
        use_analytical_stress = .true.
      end if

      if (use_analytical_stress) then
         ! Print usage.
         write(info_str,'(2X,A)') &
            "Analytical stress will be computed."
         call localorb_info(info_str,use_unit,'(A)')

         ! Print if symmetrization was chosen or not.
         if (AS_components .eq. 9) then
            write(info_str,'(2X,A)') &
               "Analytical stress calculation: All 9 components are calculated."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               "                               Final output is given unsymmetrized and symmetrized."
            call localorb_info(info_str,use_unit,'(A)')
         else
            write(info_str,'(2X,A)') &
               "Analytical stress calculation: Only the upper triangle is calculated."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               "                               Final output is symmetrized."
            call localorb_info(info_str,use_unit,'(A)')
         end if


         ! Enforce force calculation if computation of analytical stress is requested.
         if (.not. use_forces) then
            use_forces = .true.
            write(info_str,'(2X,A,A)') &
               "Analytical stress calculation: Force calculation was not requested. ", &
               "However, forces are needed and will be computed therefore."
            call localorb_info(info_str,use_unit,'(A)')
         end if

         ! Set default sc_accuracy_stress if it is not given but computation of analytical stress is requested.
         if (.not. flag_acc_stress) then
            sc_accuracy_stress = -5d0
            write(info_str,'(2X,A)') &
               "Analytical stress calculation: scf convergence accuracy of stress not set. "
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               "                               Analytical stress self-consistency will not be checked explicitly."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               "                               Be sure to set other criteria like sc_accuracy_rho tight enough."
            call localorb_info(info_str,use_unit,'(A)')
         end if

         ! Check for non periodic system
         if (n_periodic .eq. 0) then
            call localorb_info('',use_unit,'(A)')
            write(info_str,'(1X,A)') &
               "*** Analytical stress calculation: Stress calculations are not possible for non-periodic systems,"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
               "***                                since non-periodic systems do not have a unit cell and thus the"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
               "***                                concept of a stress tensor has no meaning.  Please modify your "
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
               "***                                control.in to remove the request to calculate the stress tensor."
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(1X,A)') &
               "***                                Exiting."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! Check for vdw
         if (use_vdw_method .or. use_vdw_correction) then
            write(info_str,'(2X,A)') &
               "*** Analytical stress calculation: Stress calculation only possible with the Hirschfeld vdW correction."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! Current problem with meta-GGAs and spin-polarised systems. AJL, Jan 2018
         if (n_spin.gt.1.and. &
             (flag_xc.eq.25 .or. & ! m06-l
              flag_xc.eq.26 .or. & ! m06
              flag_xc.eq.27 .or. & ! m06-2x
              flag_xc.eq.28 .or. & ! m06-hf
              flag_xc.eq.51 .or. & ! tpss
              flag_xc.eq.52 .or. & ! revtpss
              flag_xc.eq.53 .or. & ! tpssloc
              flag_xc.eq.55 .or. & ! m08-hx
              flag_xc.eq.56 .or. & ! m08-so
              flag_xc.eq.58 .or. & ! m11-l
              flag_xc.eq.59)) then ! scan
               write(info_str,'(2X,A)') &
               "*** Analytical stress calculation: Stress calculation only possible with spin-paired meta-GGAs."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! Check for used functional
         if (flag_xc.ne. 0 .and. & ! hf
             flag_xc.ne. 1 .and. & ! pbe0
             flag_xc.ne. 3 .and. & ! pz-lda
             flag_xc.ne. 6 .and. & ! pbe
             flag_xc.ne. 7 .and. & ! hse
             flag_xc.ne. 8 .and. & ! pw-lda
             flag_xc.ne. 9 .and. & ! blyp
             flag_xc.ne.10 .and. & ! b3lyp
             flag_xc.ne.11 .and. & ! rpbe
             flag_xc.ne.12 .and. & ! revpbe
             flag_xc.ne.13 .and. & ! pbesol0
             flag_xc.ne.15 .and. & ! vwn
             flag_xc.ne.16 .and. & ! vwn-gauss
             flag_xc.ne.17 .and. & ! pbesol
             flag_xc.ne.18 .and. & ! r48pbe
             flag_xc.ne.20 .and. & ! am05
             flag_xc.ne.21 .and. & ! pbeint
             flag_xc.ne.22 .and. & ! xpbe
             flag_xc.ne.25 .and. & ! m06-l
             flag_xc.ne.26 .and. & ! m06
             flag_xc.ne.27 .and. & ! m06-2x
             flag_xc.ne.28 .and. & ! m06-hf
             flag_xc.ne.29 .and. & ! b1lyp
             flag_xc.ne.51 .and. & ! tpss
             flag_xc.ne.52 .and. & ! revtpss
             flag_xc.ne.53 .and. & ! tpssloc
             flag_xc.ne.55 .and. & ! m08-hx
             flag_xc.ne.56 .and. & ! m08-so
             flag_xc.ne.58 .and. & ! m11-l
             flag_xc.ne.59 .and. & ! scan
             flag_xc.gt.-10 &      ! libxc
             .and. .not. is_dfauto(flag_xc) &
         ) then
            write(info_str,'(2X,A)') &
               "*** Analytical stress calculation: Stress calculation only possible for LDA, GGA, Meta-GGA and hybrid functionals."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if

         ! Check for relativistic settings and set use_AS_Jac_in_pulay
         if (flag_rel .eq. REL_zora) then
            write(info_str,'(2X,3A)') &
               "*** Analytical stress calculation: Stress calculation only", &
               " possible in the non-relativistic case or with atomic ZORA", &
               " (scalar)."
            call localorb_info(info_str,use_unit,'(A)')
            call aims_stop_coll('', func)
         end if
         ! FIXME: use_AS_Jac_in_pulay should become a user definable
         !        variable, since critical to accuracy for large systems
         if (flag_rel .eq. REL_none) then
            use_AS_Jac_in_pulay = .true.
         else if ( flag_rel .eq. REL_atomic_zora) then
            use_AS_Jac_in_pulay = .true.
         end if

      end if

      !SK: Since we decided to couple sc_accuracy to use_forces, this is the
      !    place where use_forces is not altered anymore. So, I put it here.

      if (.not.flag_acc_rho) then
         ! Decide on a default value for sc_accuracy_rho that reflects
         ! empirical experience for good convergence.
         if (n_atoms.le.6) then
            sc_accuracy_rho = 1.0d-6
         else if (use_forces .or. use_analytical_stress) then
            if (n_atoms.le.6000) then
               sc_accuracy_rho = 1.0d-6 * sqrt(dble(n_atoms)/6.d0)
            else ! safety check because sc_accuracy_rho really was not tested
                 ! for situations above 6000 atoms
               sc_accuracy_rho = 1.0d-6 * sqrt(1000.d0)
            end if
         else
            if (n_atoms.le.6000) then
               sc_accuracy_rho = 1.0d-5 * dble(n_atoms)/60.d0
            else ! safety check because sc_accuracy_rho really was not tested
                 ! for situations above 6000 atoms
               sc_accuracy_rho = 1.0d-3
            end if
         end if

         write(info_str,'(1X,A,A)') &
               "* Notice: The s.c.f. convergence criterion sc_accuracy_rho was", &
               " not provided."
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A)') &
               "* We used to stop in this case, and ask the user to provide reasonable"
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A)') &
               "* scf convergence criteria. However, this led some users to employ criteria"
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A)') &
               "* that led to extremely long run times, e.g., for simple relaxation."
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A)') &
               "* We now preset a default value for sc_accuracy_rho if it is not set."
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A)') &
               "* You may still wish to check if this setting is too tight or too loose for your needs."
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(1X,A,E13.6,A)') &
               "* Based on n_atoms, FHI-aims chose sc_accuracy_rho = ", sc_accuracy_rho, " ."
         call localorb_info(info_str,use_unit,'(A)')

         ! Check the accuracy of the forces and shout or stop if a waste of time can happen
         if ( flag_acc_forces ) then
            if (sc_accuracy_forces .gt. 0 .and. &
                sc_accuracy_forces .lt. 1.0d-6) then
               write(info_str,'(2X,A,1X,A,1X,E10.3,A)') &
                   "*** You did choose a very tight ", &
                   "force convergence criterion (sc_accuracy_forces) of ", &
                   sc_accuracy_forces, " ,"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,1X,A,1X,A)') &
                   "*** but did not specify sc_accuracy_rho. This", &
                   "would most certainly lead to unnecessary many force"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                   "*** evaluations, slowing the computation down (see manual)."
               call localorb_info(info_str,use_unit,'(A)')
               write (info_str,'(2X,A,A)') &
                   "*** Please carefully choose an sc_accuracy_rho criterion ", &
                   "and modify your input file control.in accordingly."
               call localorb_info(info_str,use_unit,'(A)')
               call localorb_info('',use_unit,'(A)')
               call aims_stop_coll('', func)
            else if (sc_accuracy_forces .gt. 0 .and. &
                sc_accuracy_forces .lt. 1.0d-4) then
               write(info_str,'(2X,A,1X,A,1X,E10.3,A)') &
                   "*** WARNING: You did choose a rather tight", &
                   "force convergence criterion (sc_accuracy_forces) of", &
                   sc_accuracy_forces, " ,"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A,1X,A)') &
                   "*** but did not specify sc_accuracy_rho. Please", &
                   "be aware that this can lead to unnecessary many force"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
                   "*** evaluations, slowing the computation down (see manual)."
               call localorb_info(info_str,use_unit,'(A)')
            end if
         end if


      end if !(.not.flag_acc_rho)


      if (use_forces) then

         modified = .false.
         if (.not.flag_acc_forces) then
            ! The default behavior is now to no longer check the force
            ! convergence behavior explicitly and rely on the sc_accuracy_rho default choice instead.
            ! Past user-based choices could lead to immensely long run times, but
            ! some users remained unaware that their own input choice had this effect.
            sc_accuracy_forces = -5.d0 ! this is equivalent to 'not_checked' in control.in
            modified = .true.
            write(info_str,'(2X,A,A)') &
               "Force calculation: scf convergence accuracy ", &
               "of forces not set."
            call localorb_info ( info_str )
         end if

         ! This trap is now no longer enforced. Convergence will no longer be checked explicitly.
         !
         !if (relax_mode.ne.RELAX_OFF) then
         !
         !   if (relax_accuracy_forces.lt.(sc_accuracy_forces*5.d0)) then
         !      sc_accuracy_forces = relax_accuracy_forces/5.d0
         !      modified = .true.
         !      write(info_str,'(2X,A,A)') &
         !         "Required force accuracy for relaxation overrides ", &
         !         "required scf force convergence accuracy."
         !      call localorb_info ( info_str )
         !   end if
         !
         !end if

         if (modified) then
            write(info_str,'(2X,A)') &
               "Defaulting to 'sc_accuracy_forces not_checked'."
            call localorb_info ( info_str )
         end if

         if ((relax_mode.ne.RELAX_OFF) .or. &
         &   use_molecular_dynamics .or. &
         &   final_forces_cleaned) then

            if ((.not.use_embedding_potential) .and. .not. &
               (flag_remove_unitary_force_component)) then
               remove_unitary_force_components = 2
               write(info_str,'(2X,A,A)') &
                  "Handling of forces: Unphysical translation and rotation", &
                  " will be removed from forces."
               call localorb_info ( info_str )
            else if (use_embedding_potential) then
               remove_unitary_force_components = 0
               write(info_str,'(2X,A,A)') &
                  "Handling of forces in the presence of an external ", &
                  "potential:"
               call localorb_info ( info_str )
               write(info_str,'(2X,A,A)') &
                  "Simple unitary transformations", &
                  " will not be removed from forces."
               call localorb_info ( info_str )
            end if

         end if

         if(use_embedding_pp.and.(.not.add_embedding_grids).and.use_forces) then
            ! in this case, we get indexing problems in the force computation (e.g. with
            ! 'species2pp_species(species(center_to_atom(current_center)))' in sum_up_whole_potential.
          !  call aims_stop_coll('Attention - always set add_embedding_grids .true. when computing QM/MM forces.', func)
         endif

      end if

      ! Checks concerning MD
      if (flag_molecular_dynamics) then

         ! ask for time step
         if (.not.flag_MD_tstep) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') '* Error: Please specify a time step for your molecular dynamics run.'
               write(use_unit,'(1X,A)') '* We used to set a default time step of 1 fs, equivalent to the request'
               write(use_unit,'(1X,A)') '* "MD_time_step 0.001" in control.in. 1 fs will work for most purposes,'
               write(use_unit,'(1X,A)') '* but we do prefer that you specify the desired step (=accuracy) for'
               write(use_unit,'(1X,A)') '* your run, instead of us. For example, very accurate hydrogen stretch'
               write(use_unit,'(1X,A)') '* motions in Hydrogen bonded systems will require less than 1 fs, '
               write(use_unit,'(1X,A)') '* while MD involving only very heavy elements (low vibrational frequencies)'
               write(use_unit,'(1X,A)') '* may be perfectly happy with 3 or 4 fs for the time step.'
            end if
            call aims_stop_coll('', func)
         end if

         ! More than 0.1 ps for the MD time step must be a misprint.
         if (MD_tstep.gt.0.1) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A,F14.5,A)') '* Attention! You have requested a time step of ', MD_tstep, &
                 "ps (!) for molecular dynamics."
               write(use_unit,'(1X,A,F14.5,A)') '* That would be ', MD_tstep*1000, 'fs.'
               write(use_unit,'(1X,A)') '* While you are free to as you wish, this sounds like a misprint,'
               write(use_unit,'(1X,A)') '* perhaps a mix-up of femtoseconds and picoseconds (the unit expected'
               write(use_unit,'(1X,A)') '* in control.in is picoseconds). We stop here to make sure that such'
               write(use_unit,'(1X,A)') '* typos have no immediate consequences. Apologies for any inconvenience.'
            end if
            call aims_stop_coll('', func)
         end if

         ! Ensure initialization of velocities

         if ((.not.MB_velocity_initialization) .and. &
         &   (.not.MD_initialconf_from_restart) .and. &
         &   (.not.use_MD_QH_init) .and. &
         &   (.not.MD_velocity_in_input)) then

            if (.not.flag_MD_schedule) then
               select case (MD_ensemble)
               case('NVT_berendsen', 'NVT_nose-hoover', 'NVT_nose-poincare', &
               &    'NVT_andersen', 'NVT_parrinello')
                  write(info_str,'(2X,2A)') &
                  & '* WARNING: missing velocity initialization for MD run!', &
                  & ' Using default temperature.'
                  call localorb_info(info_str)
                  MB_velocity_initialization = .true.
                  MD_init_temperature = MD_temperature
               case('NVE', 'NVE_4th_order', 'NVE_damped')
                  ! absolutely no velocity specs and nothing to get it from, need to break here:
                  write(info_str,'(2X,2A)') &
                  & "* WARNING: missing velocity specification in an NVE molecular dynamics.",&
                  & " Please fix. Aborting."
                  call localorb_info(info_str)
                  call aims_stop_coll('', func)
               case default
                  write(info_str, "(2X,2A)") &
                  & '*** Internal error; unknown ensemble ', trim(MD_ensemble)
                  call localorb_info(info_str)
               end select
            else
               select case (MD_schedule_ensemble(1))
               case('NVT_berendsen', 'NVT_nose-hoover', 'NVT_nose-poincare', &
               &    'NVT_andersen', 'NVT_parrinello')
                  write(info_str,'(2X,2A)') &
                  & '* WARNING: missing velocity initialization for MD run!', &
                  & ' Using default temperature from first segment.'
                  call localorb_info(info_str)
                  MB_velocity_initialization = .true.
                  MD_init_temperature = MD_schedule_temperature(1)
               case('NVE', 'NVE_4th_order', 'NVE_damped')
                  ! absolutely no velocity specs and nothing to get it from, need to break here:
                  write(info_str,'(2X,2A)') &
                  & "* WARNING: missing velocity specification in an NVE molecular dynamics.",&
                  & " Please fix. Aborting."
                  call localorb_info(info_str)
                  call aims_stop_coll('', func)
               case default
                  write(info_str, "(2X,2A)") &
                  & '*** Internal error; unknown ensemble ', trim(MD_ensemble)
                  call localorb_info(info_str)
               end select
            end if
         end if

         ! change thermostat mass to internal MD units

         if (MD_thermostat_units.eq.'cm^-1') then
            if (.not. flag_MD_schedule) then
               if ((MD_ensemble.eq.'NVT_nose-hoover') .or. &
               &   (MD_ensemble.eq.'NVT_nose-poincare')) then
                  MD_Q_NP = MD_Q_conversion*dble(3*n_atoms)*MD_temperature/(MD_Q_NP*MD_Q_NP)
                  write(info_str,'(2X,A,E14.4E4,A)') 'Changed MD thermostat mass units, Q = ',&
                  MD_Q_NP,' amu*bohr^2'
                  call localorb_info(info_str)
               end if
            else
               do MD_schedule_step = 1, MD_segments  ! check each segmend for Nose-thermostats
                  if ((MD_schedule_ensemble(MD_schedule_step).eq."NVT_nose-hoover").or.&
                  &   (MD_schedule_ensemble(MD_schedule_step).eq."NVT_nose-poincare"))  then
                     MD_schedule_Q(MD_schedule_step) = MD_Q_conversion*dble(3*n_atoms) &
                     &                                 * MD_schedule_temperature(MD_schedule_step) &
                     &                                 / (MD_schedule_Q(MD_schedule_step)**2)
                     write(info_str,'(2X,A,I3,A,E14.4E4,A)') &
                     & 'Changed MD thermostat mass units for schedule segment ', &
                     & MD_schedule_step,', Q = ',MD_schedule_Q(MD_schedule_step),' amu*bohr^2'
                     call localorb_info(info_str)
                  end if
               end do
            end if
         end if

         ! rescale NVT_parrinello tau parameter

         if (.not. flag_MD_schedule) then
            if (MD_ensemble.eq.'NVT_parrinello') then
               MD_tau_BDP = MD_tau_BDP / MD_tstep ! rescaling tau
               write(info_str,'(2X,A,E14.4E4)') &
               & 'Bussi-Donadio-Parrinello thermostat relax. time rescaled with MD time step; tau = ',&
               & MD_tau_BDP
               call localorb_info(info_str)
            end if
         else
            do MD_schedule_step = 1, MD_segments  ! check each segmend for BDP thermostats and change units of tau
               if (MD_schedule_ensemble(MD_schedule_step).eq."NVT_parrinello") then
                  MD_schedule_tau_BDP(MD_schedule_step) = MD_schedule_tau_BDP(MD_schedule_step) / MD_tstep
                  write(info_str,'(2X,A,I3,A,E14.4E4)') 'BDP thermostat tau rescaled with MD time step for sched. segm.', &
                  MD_schedule_step, '; tau = ',MD_schedule_tau_BDP(MD_schedule_step)
                  call localorb_info(info_str)
               end if
            end do
         end if
      end if ! (flag_molecular_dynamics)

   !     check thermodynamic integration setting
      if (use_thermodynamic_integration) then
         !Safety check, user FHI-AIMS should already have died in this case
         if (MD_segments .ne. TDI_segments) then
            if (myid.eq.0) then
               write(use_unit,*) '  *** ERROR: '
               write(use_unit,*) &
               '  *** MD_schedule: Either all or none of the MD_segments must have a thermodynamic_integration specification.'
               write(use_unit,*) '  *** Aborting execution.'
            end if
            call aims_stop
         end if

         !Iterate over all MD_segments
         do MD_schedule_step = 1, MD_segments
            !Check for lambda values
            if ( (TDI_lambda_start(MD_schedule_step).lt.0d0).or.(TDI_lambda_start(MD_schedule_step).gt.1d0) ) then
               if (myid.eq.0) then
                  write(use_unit,*) '  *** ERROR: '
                  write(use_unit,*) '  *** thermodynamic_integration: Starting value for interpolation parameter'
                  write(use_unit,*) '  *** lambda is not between 0.0 and 1.0 in segment ',MD_schedule_step, '.'
                  write(use_unit,*) '  *** Aborting execution.'
               end if
               call aims_stop
            end if
            if ( (TDI_lambda_end(MD_schedule_step).lt.0d0).or.(TDI_lambda_end(MD_schedule_step).gt.1d0) ) then
               if (myid.eq.0) then
                  write(use_unit,*) '  *** ERROR: '
                  write(use_unit,*) '  *** thermodynamic_integration: Stopping value for interpolation parameter'
                  write(use_unit,*) '  *** lambda is not between 0.0 and 1.0 in segment ',MD_schedule_step, '.'
                  write(use_unit,*) '  *** Aborting execution.'
               end if
               call aims_stop
            end if
            !Check if force constants file exists, unless 'null'  is specified for reference free AS (only for clusters)
            inquire(FILE=TDI_QHA_file(MD_schedule_step),EXIST=TDI_QHA_file_exists(MD_schedule_step))
            if (.not.TDI_QHA_file_exists(MD_schedule_step)) then
               if (myid.eq.0) then
                  write(use_unit,*) '  *** ERROR: '
                  write(use_unit,*) '  *** thermodynamic_integration: File ',trim(TDI_QHA_file(MD_schedule_step))
                  write(use_unit,*) '  *** required in MD_segment ',MD_schedule_step, 'does not exist.'
                  write(use_unit,*) '  *** Aborting execution.'
               end if
               call aims_stop
            end if
         end do
      end if

   !     check AS setting

      ! Check for the availability of classical fields when
      ! enabling skip_SCF
      if (SKIP_SCF) then

         if ( (.not.use_thermodynamic_integration) .and. (.not.use_harmonic_pot_only) ) then
            if (myid.eq.0) then
               write(use_unit,'(A)') '  *** ERROR: '
               write(use_unit,'(A)') '  *** Option skip_SCF can only be used if an alternative'
               write(use_unit,'(A)') '  *** (classical) force field is specified.'
               write(use_unit,'(A)') '  *** Aborting execution.'
            end if
            call aims_stop
         end if

      end if


   !     check flag for partition function

      if (.not.flag_partition_acc) then
         partition_acc = 1.d-15

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,E11.4,A)') &
                  "No accuracy limit for integral partition fn. given. ", &
                  "Defaulting to ", partition_acc, "."
         end if

      end if

   !     check flag for threshold value of u(r)

      if (.not.flag_wave_threshold) then
         wave_threshold = 1.d-6

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,E11.4,A)') &
                  "No threshold value for u(r) in integrations given. ", &
                  "Defaulting to ", wave_threshold, "."
         end if

      end if

      if (.not.flag_scsmp2) then
         pt_mp2=1.d0
         ps_mp2=1.d0
      endif


      if (flag_atom_ref) then
      ! This flag forces specific occupation parameters, and a temporary
      ! field for symmetry breaking. The following must therefore appear instead
      ! of the default settings for smearing.

         if ((ext_safe.eq.'safe').or.(.not.flag_occupation_type)) then
            ! Overwrite default occupation settings
            if (occupation_type .ne. 3) then  ! if integer occupation is specified, do not overwrite IGOR
                occupation_type = 0
                occupation_width=1.d-5
                occupation_width_set = .true.
            end if
         end if

         if (myid.eq.0) then
            write (use_unit,*)
            write (use_unit,'(2X,A,A)') &
            "Attention! Symmetry-breaking settings ", &
            "requested above."

            write (use_unit,'(2X,A)') &
            "Enforcing the following defaults: "

            write (use_unit,'(4X,A,A)') &
            "Occupation (smearing) type: ", &
            "Gaussian broadening"

            write (use_unit,'(4X,A,E14.6,A)') &
            "Gaussian broadening width : " &
            ,occupation_width, " eV."

            write (use_unit,*)

         end if

         occupation_width = occupation_width / hartree

      else if (.not.flag_occupation_type) then
   !       use default: Gaussian smearing
         occupation_type = 0
         occupation_width = 0.01

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A,E11.4,A)') &
                  "No occupation type (smearing scheme) given. ", &
                  "Defaulting to Gaussian broadening, width = ", &
                  occupation_width, " eV."
         end if

         occupation_width = occupation_width/hartree

         if (adjust_scf) then
            if (adjust_scf_always) then
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,I5,A)') &
                   "The width will be adjusted in iteration number ", &
                   adjust_scf_iteration, " of each s.c.f. cycle (i.e., each geometry step)."
               end if
           else
               if (myid.eq.0) then
                  write (use_unit,'(2X,A,I5,A)') &
                  "The width will be adjusted in iteration number ", &
                  adjust_scf_iteration, " of the first full s.c.f. cycle only."
               end if
            end if
         end if

         if (myid.eq.0) then
             write (use_unit,'(2X,A,I5,A)') &
               "S.C.F. convergence parameters will be adjusted in iteration number ", &
             adjust_scf_iteration, " of the first full s.c.f. cycle only."
         end if

      end if

      !if (.not.flag_DFPT_occupation_type) then
      !   DFPT_width = 0.01
      !   write (info_str,'(2X, A, E14.6, A)') &
      !       "No width specified for DFPT occupation numbers. Using default: DFPT_width =", DFPT_width, "eV"
      !   call localorb_info(info_str,use_unit,'(A)')
      !   DFPT_width = DFPT_width/hartree
      !end if

      ! ELSI occupation not working with use_constraint, force_occupation_*, and
      ! 'integer' occupation type.
      ! FHI-aims occupation not working with 'cubic' and 'cold' occupation types.
      if (.not. flag_mu_method) then
         if(use_constraint .or. force_occupation_basis &
            .or. force_occupation_projector .or. occupation_type == 3) then
            mu_method = 0 ! zeroin
         else
            mu_method = 1 ! bisection
         endif
      else
         if(mu_method == 0) then
            if(occupation_type == 4) then
               if(myid == 0) then
                  write(use_unit,'(1X,A)') &
                     "* 'cubic' occupation type only available with ELSI occupation."
               endif

               call aims_stop_coll('',func)
            elseif(occupation_type == 5) then
               if(myid == 0) then
                  write(use_unit,'(1X,A)') &
                     "* 'cold' (Marzari-Vanderbilt) occupation type only available with ELSI occupation."
               endif

               call aims_stop_coll('',func)
            endif
         elseif(mu_method == 1) then
            if(use_constraint) then
               if(myid == 0) then
                  write(use_unit,'(1X,A)') &
                     "* ELSI occupation not supported with 'use_constraint'."
               endif

               call aims_stop_coll('',func)
            endif

            if(force_occupation_basis .or. force_occupation_projector) then
               if(myid == 0) then
                  write(use_unit,'(1X,A)') &
                     "* ELSI occupation not supported with 'force_occupation_*'."
               endif

               call aims_stop_coll('',func)
            endif

            if(occupation_type == 3) then
               if(myid == 0) then
                  write(use_unit,'(1X,A)') &
                     "* ELSI occupation not supported with 'integer' occupation type."
               endif

               call aims_stop_coll('',func)
            endif
         endif
      endif

   !     check flag for occupation accuracy settings
      if (.not.flag_occupation_acc) then
         if(mu_method == 1) then
            occupation_acc = 1e-13
         else
            occupation_acc = 1e-8
         endif

         if (myid.eq.0) then
            write (use_unit,'(2X,A,A,E11.4,A)') &
                  "No accuracy for occupation numbers given. ", &
                  "Defaulting to ", occupation_acc, "."
         end if

      end if
      if (.not.flag_occupation_thr) then
         occupation_thr = 0.d0

         if (myid.eq.0) then
            write (use_unit,'(2X,A,A,E11.4,A)') &
                  "No threshold value for occupation numbers given. ", &
                  "Defaulting to ", occupation_thr, "."
         end if

      end if
      if (.not.flag_fermi_acc) then
         fermi_acc = 1e-20

         if (myid.eq.0) then
            write (use_unit,'(2X,A,A,E11.4,A)') &
                  "No accuracy for fermi level given. ", &
                  "Defaulting to ", fermi_acc, "."
         end if

      end if
      if (.not.flag_max_zeroin) then
         max_zeroin = 200

         if (myid.eq.0) then
            write (use_unit,'(2X,A,A,I4,A)') &
                  "Maximum # of iterations to find E_F not set. ", &
                  "Defaulting to ", max_zeroin, "."
         end if

      end if

   !     Set eigenvalue solver defaults in a sane way:
      if (.not.flag_KS_method) then
         write (info_str,'(2X,A)') &
               "Preferred method for the eigenvalue solver ('KS_method') not specified in 'control.in'."
         call localorb_info(info_str,use_unit,'(A)')
         flag_KS = -1
         use_elsi = .true.
         if ( (use_mpi) .and. (n_tasks.gt.n_kpts) ) then
            ! switch to elpa by default when mpi is compiled in
            ! and more than one processor is used

            call localorb_info('  Calling BLACS routine to test compilation state')
            call MPI_Barrier(mpi_comm_global, info)

            ! check if scalapack_support was compiled in
            call BLACS_Pinfo(index1,index2)
            if (index2.gt.0) then
               ! scalapack support is there - use elpa
               use_scalapack = .true.
               use_elpa = .true.
               write (info_str,'(2X,A)') &
                     "Since ScaLAPACK support is enabled, defaulting to ELPA (via ELSI)."
               call localorb_info(info_str,use_unit,'(A)')
            else
               ! use lapack
               use_lapack_fast = .true.
               write (info_str,'(2X,A)') &
                     "The best choice would be ELPA (via ELSI), but your version of FHI-aims is compiled without ScaLAPACK support."
               call localorb_info(info_str,use_unit,'(A)')
               write (info_str,'(2X,A)') &
                     "Defaulting to LAPACK (via ELSI) instead (serial solution, time will not improve with more CPUs)."
               call localorb_info(info_str,use_unit,'(A)')
            end if
         else
            ! use lapack by default
            if (n_tasks.gt.1) then
               write (info_str,'(2X,A)') &
               "Defaulting to serial version, LAPACK (via ELSI), since more k-points than CPUs available."
               call localorb_info(info_str,use_unit,'(A)')
               use_lapack_fast = .true.
            else
               write (info_str,'(2X,A)') &
               "Defaulting to serial version LAPACK (via ELSI)."
               call localorb_info(info_str,use_unit,'(A)')
               use_lapack_fast = .true.
            end if
         end if
      else
         ! This is the case where KS_method was set explicitly. Must now check if the reqested setting makes sense.

         !     If ScaLapack was requested, perform basic checks to determine if
         !     this even makes sense. Else, switch to Lapack.

         if ((use_scalapack).and.(.not.use_mpi)) then
            ! Check if this even an MPI run
            if (.not.use_elpa) then
               write (info_str,'(1X,A,A)') &
               '* ScaLAPACK solver makes sense only with MPI runs. ', &
               'Switching to LAPACK (via ELSI) instead.'
               call localorb_info(info_str,use_unit,'(A)')
               use_scalapack = .false.
               use_lapack_fast = .true.
            end if
         end if

         if (use_scalapack) then
            ! Check if ScaLapack exists in the binary
            call BLACS_Pinfo(index1,index2)
            if (index2.le.0) then
               use_scalapack = .false.
               write(info_str,'(2A)') '* WARNING: Requested ScaLAPACK run', &
               ' without having it compiled into the binary.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2A)') '* Using LAPACK (via ELSI) instead.'
               call localorb_info(info_str,use_unit,'(A)')
               use_lapack_fast = .true.
            end if
         end if

         if (use_scalapack) then
            ! Check if we have "few enough" k-points to even run ScaLapack.
            if (flag_k_points_defined) then
               if ( n_tasks .lt. n_kpts ) then
                  write(info_str, '(1X,A,I5,A,I5,A)') '* You have ', n_kpts, &
                     ' k-points in your system, but only ', n_tasks, &
                     ' MPI-tasks available.'
                  call localorb_info(info_str,use_unit, '(A)')

                  write(info_str,'(1X,A,A)') &
                     '* In this case using ScaLAPACK ', &
                    'is not necessary / supported.'
                  call localorb_info(info_str,use_unit, '(A)')

                  write(info_str,'(1X,A)') &
                     '* Switching to LAPACK (via ELSI) instead.'
                  call localorb_info(info_str,use_unit, '(A)')

                  use_scalapack = .false.
                  use_lapack_fast = .true.

                  ! Check for load_balancing first, since it will automatically set use_local_index in some cases
                  if (use_load_balancing) then
                     write(info_str,'(1X,A)') &
                        '* Also switching off load_balancing and use_local_index.'
                        call localorb_info(info_str,use_unit, '(A)')
                        use_load_balancing = .false.
                        use_local_index = .false.
                  end if
                  if (use_local_index) then
                     write(info_str,'(1X,A)') &
                        '* Also switching off use_local_index.'
                        call localorb_info(info_str,use_unit, '(A)')
                        use_local_index = .false.
                  end if
               else if (.not.flag_max_tasks_per_smp_node) then
                  max_tasks_per_smp_node = 512
                  write(info_str,'(2X,A,A)') &
                        "Maximum number of tasks per SMP node", &
                        "not specified."
                  call localorb_info(info_str,use_unit,'(A)')
                  write(info_str,'(2X,A,1X,I5,1X,A)') &
                        "Presetting", max_tasks_per_smp_node, "tasks."
                  call localorb_info(info_str,use_unit,'(A)')
               end if
            end if
         end if
      end if ! KS_method set in control.in?

      if(use_elsi) then
         if(use_constraint) then
            use_elsi = .false.

            write (info_str,'(2X,A)') &
               "Switching off ELSI because of 'use_constraint'."
            call localorb_info(info_str,use_unit,'(A)')
         else
            if(.not. flag_elsi_method) then
               use_elsi_ev = .true.
               use_elsi_dm = .false.
            endif
         endif

         if(output_level == "MD_light") then
            elsi_out_level = 0
         endif
      else
         elsi_out_level = 0
      endif

      if (use_mbd_old .and. (.not. flag_mbd_eigensolver) ) then
       write (info_str,'(2X,A)') &
               "Preferred method for the MBD eigenvalue solver ('mbd_eigensolver') not specified in 'control.in'."
         call localorb_info(info_str,use_unit,'(A)')
         if ( (use_mpi) .and. (n_tasks.gt.1) ) then
            ! switch to scalapack by default when mpi is compiled in
            ! and more than one processor is used

            call localorb_info('  Calling blacs routine to test compilation state')
            call MPI_Barrier(mpi_comm_global, info)

            ! check if scalapack_support was compiled in
            call BLACS_Pinfo(index1,index2)
            if (index2.gt.0) then
               ! scalapack support is there - use scalapack
               use_scalapack_mbd = .true.
               use_elpa_mbd      = .false.
               write (info_str,'(2X,A)') &
                     "Since scalapack support is enabled, defaulting to 'scalapack'."
               call localorb_info(info_str,use_unit,'(A)')
            else
               ! use lapack
               use_scalapack_mbd = .false.
               use_elpa_mbd      = .false.
               write (info_str,'(2X,2A)') &
                 "The best choice would be 'elpa', but your version of FHI-aims ", &
                 "is compiled without scalapack support."
               call localorb_info(info_str,use_unit,'(A)')
            end if
        else
          ! use lapack
               use_scalapack_mbd = .false.
               use_elpa_mbd      = .false.
               write (info_str,'(2X,A)') &
                 "Serial version of FHI-aims: dafault is LAPACK."
               call localorb_info(info_str,use_unit,'(A)')

        end if
      end if

      if (use_mbd_dev) then
          if (mbd_dev_verify() /= 0) call aims_stop_coll('Incosistent MBD settings', func)
      end if

      ! Make sure load_balancing, when enabled, matches use_local_index
      if (flag_load_balancing_if_scalapack.and..not.flag_use_local_index_if_scalapack) then
         call aims_stop_coll("The keyword 'load_balancing if_scalapack' in control.in requires &
                             &the keyword 'use_local_index if_scalapack' also be present.  Exiting.", func)
      end if
      if (use_load_balancing.and..not.use_local_index) then
         call aims_stop_coll("The keyword 'load_balancing .true.' in control.in requires &
                             &the keyword 'use_local_index .true.' also be present.  Exiting.", func)
      end if

      ! Determine the status of domain decomposition and load balancing when using the if_scalapack keyword
      if (.not. flag_use_local_index) then
         if (flag_use_local_index_if_scalapack) then
            use_local_index = use_scalapack
            if (use_local_index) then
               call localorb_info("  Have scalapack: Storing only matrices local to each thread.")
            else
               call localorb_info("  No scalapack: No use_local_index.")
            end if
         end if
      end if
      if (.not. flag_load_balancing) then
         if (flag_load_balancing_if_scalapack) then
            use_load_balancing = use_scalapack
            if (use_load_balancing) then
               call localorb_info("  Have scalapack: Using load balancing for integrations.")
            else
               call localorb_info("  No scalapack: No load_balancing.")
            end if
         end if
      end if

      ! When running on a system with many CPUs, it can be much faster to use mpi_alltoallv
      ! than doing n_tasks time a sendrecv call; however this costs a significant amount of memory.
      ! Thus we switch to using mpi_alltoallv when using 1024 CPUs or more (the number 1024 may be adapted of course!)

      if (.not. flag_use_alltoall) then
         if (n_tasks>=1024) then
            call localorb_info("  Will use alltoall communication since running on >= 1024 CPUs.")
            call localorb_info("  *** If you encounter memory problems, try setting this to .false.")
            use_alltoall = .true.
         else
            call localorb_info("  Will not use alltoall communication since running on < 1024 CPUs.")
            use_alltoall = .false.
         end if
      end if


   !     After eigensolver is set, must check whether we are requiring eigenvectors to be written out
   !     for any reason.
      if (out_eigenvec) then

         write(info_str,'(2X,A)') &
         'Output of Kohn-Sham eigenvector components into separate files was explicitly requested.'
         call localorb_info(info_str,use_unit, '(A)')

         if (n_periodic.gt.0) then
         if (.not.out_band) then
            write(info_str,'(1X,A)') &
               "* Warning: Eigenvector output for periodic systems happens only for specified bands, "
            call localorb_info(info_str,use_unit, '(A)')
            write(info_str,'(1X,A)') &
               "* but no band plotting was requested. Not writing any eigenvectors as a result."
            call localorb_info(info_str,use_unit, '(A)')
            out_eigenvec = .false.
         else
            write(info_str,'(2X,A)') &
               "| Periodic system: Eigenvector output will happen at the very end of the run, for "
            call localorb_info(info_str,use_unit, '(A)')
            write(info_str,'(2X,A)') &
               "| those bands whose output was requested. "
            call localorb_info(info_str,use_unit, '(A)')
         end if
         end if

         if (use_scalapack .and. (.not.collect_eigenvectors) ) then
         collect_eigenvectors = .true.
         write(info_str,'(1X,A)') &
            "* Warning: Eigenvector output requires 'collect_eigenvectors .true.', but '.false' "
         call localorb_info(info_str,use_unit, '(A)')
         write(info_str,'(1X,A)') &
            "* was explicitly requested. Overriding this choice and collecting eigenvectors after all,"
         call localorb_info(info_str,use_unit, '(A)')
         write(info_str,'(1X,A)') &
            "* but for very large systems, this change may lead to memory problems."
         call localorb_info(info_str,use_unit, '(A)')
         write(info_str,'(1X,A)') &
            "* Please inform us if this becomes a problem."
         call localorb_info(info_str,use_unit, '(A)')
         end if

      end if


!     check basis threshold flag

      if (flag_KS.eq.-1) then
         if (.not.flag_basis_threshold) then
   !         use default
            basis_threshold =1.d-05

            if (myid.eq.0) then
               write(use_unit,'(2X,A,E11.4)') &
                     "Threshold for basis singularities not set."
               write(use_unit,'(2X,A,E11.4)') &
                     "Default threshold for basis singularities: ", &
                     basis_threshold
            end if

         end if
      end if

      if (.not.flag_partition_type) then
         ! Set default partition_type (stratmann_sparse) here
         partition_type = 9

         if (myid.eq.0) then
            write (use_unit,'(2X,A)') &
               "partition_type (choice of integration weights) for integrals was not specified."
            write (use_unit,'(2X,A)') &
               "| Using a version of the partition function of Stratmann and coworkers ('stratmann_sparse')."
            write (use_unit,'(2X,A)') &
               "| At each grid point, the set of atoms used to build the partition table is smoothly restricted to"
            write (use_unit,'(2X,A)') &
               "| only those atoms whose free-atom density would be non-zero at that grid point."
         end if

      end if

      if (flag_hartree_partition_type.eq.0) then
         flag_hartree_partition_type = partition_type

         if (myid.eq.0) then
         write (use_unit,'(2X,A,A)') &
            "Partitioning for Hartree potential was not ", &
            "defined. Using partition_type for integrals."
         end if

      else if (flag_hartree_partition_type.ne.partition_type) then
         ! Temporary patch: Although theoretically possible, different
         ! partition functions for the Hartree potential and integrations
         ! do not work. Stop user from falling into this trap.

         if (myid.eq.0) then
         write (use_unit,'(1X,A)') &
            "* Warning - different partition functions for the Hartree potential and the integrations "
         write (use_unit,'(1X,A)') &
            "* were explicitly requested. This is theoretically possible, but not implemented"
         write (use_unit,'(1X,A)') &
            "* at this time. While not hard to do, please contact us before trying."
         write (use_unit,'(1X,A)') &
            "* Aborting the calculation for now - sorry."
         end if
         call aims_stop_coll('', func)

      end if

      ! error trap - check and complain about incorrect partition type explicitly for
      !              self-adapting angular grids
      if ( (partition_type .eq. 5) .or. (partition_type .eq. 7) .or. (partition_type .eq. 8) .or. (partition_type .eq. 9) ) then
         ! Stratmann partition type is not yet enabled for self-adapting angular grids.
         ! A fix for this is simple but requires some programming, as the relevant evaluate_partition_* subroutine
         ! in initialize_integrals simply does not support Stratmann.

         ! For those who really need this fixed - there is also a cop-out - just use an auxiliary partition_type
         ! value that is set to 1 (rho_r2) if the actual partition_type is 5 ... but that's a hack, not a fix.

         if (.not. (ALL(specified_grid).and.ALL(angular_acc == 0.0d0)) ) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A,A)') &
                  "* Adaptive integration grids required, but these are", &
                  " not yet compatible with the Stratmann partition table."
               write(use_unit,'(1X,A,A)') &
                  "* Please use 'partition_type rho_r2' instead, and CONTACT us if you need this feature."
               write(use_unit,'(1X,A,A)') &
                  "* See also the comment in read_control.f90. Stopping program execution for now."
            end if
            call aims_stop_coll('', func)
         end if

      end if

      ! After the partition_type is definitely set, we revisit the threshold at which multipole moments are cut off
      ! in periodic systems. For the stratmann table, at least, the default must be adjusted.
      if ( ( .not.flag_multip_moments_threshold) .and. &
           ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) .or. (partition_type.eq.9) ) )  then
         ! adjust:
         far_distance_hartree_multipole_moment_threshold         = 1e-12
         if (myid.eq.0) then
            write(use_unit,'(2X,A,E20.8)') &
               "| Adjusted default value of keyword multip_moments_threshold to: ", far_distance_hartree_multipole_moment_threshold
            write(use_unit,'(2X,A)') &
               "| This value may affect high angular momentum components of the Hartree potential in periodic systems."
         end if
      else
         if (myid.eq.0) then
            write(use_unit,'(2X,A,E20.8)') &
               "| Reporting present value of keyword multip_moments_threshold  : ", far_distance_hartree_multipole_moment_threshold
            write(use_unit,'(2X,A)') &
               "| This value may affect high angular momentum components of the Hartree potential in periodic systems."
         end if
      end if

   !     verify that spin-polarisation is handled correctly
      if (.not.flag_spin) then

         if (myid.eq.0) then
            write(use_unit,'(2X,A,A)') &
               "Spin handling was not defined in control.in.", &
               " Defaulting to unpolarized case."
         end if

         spin_treatment = 0
      end if

      if (spin_treatment.ne.0) then

         ! GSM 2016/03/24 This breaks the density-matrix based spin collinear restart case due to a different size of the
         ! Hamiltonian when we use Hund's rule rather than a specified spin moment.
         !if(restart_read .and. restart_file_exists)then
            !flag_moment = .false.
            !flag_default_moment_defined = .true.
         !end if

         if (.not.flag_default_moment_defined .and. .not.(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)) then
            write(info_str,'(2X,A)') &
               "*** WARNING! This is a spin-polarized calculation, but the initial spin density"
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               "*** was not specified."
            call localorb_info(info_str,use_unit,'(A)')
            if ((n_periodic.gt.0).or.(n_atoms.gt.1)) then
               ! For a spin-polarized treatment, we now require that the user
               ! defines the initial spin polarization. Else, we stop. The only exception is a single
               ! isolated atom, for which we simply set Hund's rules.
               write(info_str,'(2X,A)') &
               "*** The spin density for the initialization of the s.c.f. cycle can either be"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** defined by setting the default_initial_moment value in control.in, or by"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** setting atom-specific initial_moment keywords in geometry.in (see the manual"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** for more information). However, what is not a good idea is to let the code"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** guess the spin initialization on its own. For example, a high-spin initialization"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** can lead to extremely slow s.c.f. convergence if the actual system is mostly"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** unpolarized. Or, a system might have multiple magnetic states - think Cr, which"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** is antiferromagnetic in the ground state. FHI-aims should not guess for you"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** which particular spin state you are aiming to describe. We therefore stop"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** the calculation and ask you to define any kind of spin initialization at all -"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** either by initial_moment, or by default_initial_moment, or both, whichever you"
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** think useful. Apologies for any inconvenience."
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') &
               "*** Exiting."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            else
               ! This is the case of a single isolated atom, the only exception, where
               ! we specify a Hund's rule initialization unless otherwise requested.
               write(info_str,'(2X,A)') &
               "*** In this case (a single, isolated atom), we default to Hund's rules (high spin)."
               call localorb_info(info_str,use_unit,'(A)')
            end if
         end if
      end if

      if (flag_moment) then
         if (spin_treatment.eq.0) then
            if (default_initial_moment.ne.0) then
               write(info_str,'(2X,A,A)') &
            "No spin polarisation - ignoring explicitly requested", &
            " default initial moment."
               call localorb_info(info_str,use_unit,'(A)')
            end if
         else
            do i_species = 1, n_species, 1
               if (abs(default_initial_moment) &
                  .gt.(species_z(i_species))) then

                  write(info_str,'(1X,A,A,A,A)') &
            "* Default initial moment per atom is larger than ", &
            "number of electrons on species ", &
                     trim(species_name(i_species)),"."
                  call localorb_info(info_str,use_unit,'(A)')

                  write(info_str,'(1X,A)') &
            "* Hazardous setting - please correct."
                  call localorb_info(info_str,use_unit,'(A)')

                  call aims_stop_coll('', func)

               end if
            enddo
         end if
      end if

!     verify settings for possible locally constrained DFT

      if (use_constraint) then

   !        if (spin_treatment.eq.0) then
   !          write(info_str,'(1X,A,A)')
   !     +      "* Locally constrained DFT for unpolarized case",
   !     +      " not yet enabled."
   !          call localorb_info(info_str,use_unit,'(A)')
   !
   !          write(info_str,'(1X,A)')
   !     +      "* Use collinear spin for constraint DFT."
   !          call localorb_info(info_str,use_unit,'(A)')
   !
   !          stop
   !        end if

         if (.not.flag_constraint_iter) then

            constraint_it_lim = 200

            write(info_str,'(2X,A)') &
            "Maximum number of constraint iterations not provided."
            call localorb_info(info_str,use_unit,'(A)')

            write(info_str,'(2X,A,1X,I5,1X,A)') &
         "  Presetting", constraint_it_lim, "iterations."
            call localorb_info(info_str,use_unit,'(A)')

         end if

      end if

      if (flag_use_metis_batching) then
         if (parallel_grid_partitioning) then
            write(info_str,'(2X,A,A,A)') &
               "Parallel grid partitioning routines don't ", &
               "support qhull+METIS batch distribution.", &
               "Disabling feature."
            call localorb_info(info_str,use_unit,'(A)')
            use_metis_batch_distribution = .false.
         end if
      end if

      if (use_cg) then
         if (.not.flag_initial_ev_solutions) then
            write(info_str,'(2X,2A,I4)') &
               'Number of inital EV solutions ', &
               'for lopcg iterations is not set. Defaulting to ', &
               initial_ev_solutions
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (.not.flag_max_cg_iterations) then
            write(info_str,'(2X,2A,I4)') 'Number of max lopcg ', &
               'iterations is not set. Defaulting to ', &
               max_cg_iterations
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (.not.flag_lopcg_tolerance) then
            write(info_str,'(2X,2A,E11.4)') 'Tolerance for lopcg ', &
               'iterations is not set. Defaulting to ', &
               lopcg_tol
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (.not.flag_lopcg_preconditioner) then
            write(info_str,'(2X,3A)') 'Preconditioner for lopcg ', &
               'iterations is not set. Defaulting to inverse ', &
               'of the overlap matrix.'
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (.not.flag_lopcg_block_size) then
            write(info_str,'(2X,2A,I4)') 'Block size for lopcg ', &
               'iterations is not set. Defaulting to ', &
               lopcg_block_size
            call localorb_info(info_str,use_unit,'(A)')
         end if
      end if

      ! consistency checks for preconditioning input
      if (.not.flag_preconditioner) then
         ! Kerker preconditioner will be used by default in periodic systems.
         ! The only exception are periodic Hartree-Fock based methods, as
         ! here, the mixed and preconditioned real-space density would be
         ! inconsistent with the mixed (but not preconditioned) density matrix.
         ! (most likely fixable by simple maths, but not yet relevant.)
         if ( (n_periodic.gt.0) .and. .not.(use_hartree_fock)) then
            use_kerker_preconditioner = .true.
         end if
      end if

      if (use_kerker_preconditioner) then
         if (precondition_kerker_q0.eq.0d0) then
            write(info_str,'(2X,2A)') '* WARNING: The precondition wave', &
                  ' vector is equal to zero, this is numerical nonsense'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') '* in the current implementation.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') '* Turning off preconditioner.'
            call localorb_info(info_str,use_unit,'(A)')
            use_kerker_preconditioner = .false.
         end if
         if (.not.flag_preconditioner_max_l) then
            precondition_max_l = 0
            ! precondition_max_l = 111
            ! do i_species = 1, n_species
            !   precondition_max_l = min(l_hartree(i_species), &
            !        precondition_max_l)
            ! end do
            write(info_str,'(2X,2A)') 'Angular momentum expansion for', &
               ' Kerker preconditioner not set explicitly.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A,I3)') '| Using default value of ', &
               precondition_max_l
            call localorb_info(info_str,use_unit,'(A)')
         end if
         if (.not.flag_preconditioner_turnoff) then
            preconditioner_turnoff_charge = sc_accuracy_rho
            use_preconditioner_turnoff_charge = .true.
            if (flag_acc_eev) then
              preconditioner_turnoff_sum_ev = sc_accuracy_eev
            end if
            if (flag_acc_etot) then
              preconditioner_turnoff_energy = sc_accuracy_etot
            end if
            write(info_str,'(2X,2A)') 'No explicit requirement for ', &
               'turning off preconditioner.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,2A)') '| By default, it will be ', &
         'turned off when the charge convergence reaches'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A,E14.6)')'| sc_accuracy_rho  = ', &
                  preconditioner_turnoff_charge
            call localorb_info(info_str,use_unit,'(A)')
            if (flag_acc_eev) then
              write(info_str,'(2X,A,E14.6)')'| sc_accuracy_eev  = ', &
                  preconditioner_turnoff_sum_ev
              call localorb_info(info_str,use_unit,'(A)')
            end if
            if (flag_acc_etot) then
              write(info_str,'(2X,A,E14.6)')'| sc_accuracy_etot = ', &
                  preconditioner_turnoff_energy
              call localorb_info(info_str,use_unit,'(A)')
            end if
         end if
         if (.not.flag_prec_mix_param) then
            if (myid.eq.0) then
               write(use_unit,'(2X,A,A)') &
                  "No special mixing parameter while Kerker", &
                  " preconditioner is on."
               write(use_unit,'(2X,A,F10.4,A)') &
                  "Using default: charge_mix_param = ", &
                  charge_mix_param(1), "."
            end if
            prec_mix_param(1) = charge_mix_param(1)
         end if
      end if
      if (flag_preconditioner_max_l.and.(.not.use_kerker_preconditioner)) then
         write(info_str,'(2X,2A)')'* Kerker preconditioning not ', &
         'requested, only a value for its angular momentum expansion.'
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(2X,2A)') '* Please check!'
         call localorb_info(info_str,use_unit,'(A)')
      end if
      if (flag_preconditioner_turnoff.and.(.not.use_kerker_preconditioner)) then
         write(info_str,'(2X,2A)')'* Kerker preconditioning not ', &
         'requested, only a stopping criterion.'
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,'(2X,2A)') '* Please check!'
         call localorb_info(info_str,use_unit,'(A)')
      end if

!  special treatment of the hybrid_coeff
      if (.not.flag_hybrid_coeff) then
         if (use_hartree_fock .and. .not.use_pbe0 &
            .and..not.use_hse .and. .not.use_b3lyp &
            .and..not.use_b1lyp &
            .and. .not.use_pbesol0 .and. .not.use_minnesota &
            .and. .not.use_libxc) then
            hybrid_coeff = 1.d0
         endif
      end if

      if (.not. use_hartree_fock .and. use_hf_kspace) then
         use_hf_kspace = .false.
      endif

      !SVL commented out, to enable realspace HF to be used in EX+RPA
      !if (use_hartree_fock .and. n_periodic .gt. 0 .and. use_rpa_ene) then
      !   use_hf_kspace = .true.
      !   use_hf_realspace = .false.
      !endif
!      write(use_unit,*) "use_hf_kspace 1", use_hf_kspace
!      if (use_corr .and. (use_hse .and. hse_omega_hf /= 0.d0)) then
!         call aims_stop_coll('Correlated methods after HSE not implemented.', func)
!      end if

      if (.not.flag_prodbas_threshold) then
         prodbas_threshold = 1.d-05
         if (use_prodbas) then
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
                  "Threshold for auxiliary basis singularities not set."
               write(use_unit,'(2X,A,E11.4)') &
                  "Default threshold for auxiliary basis singularities: ", &
                  prodbas_threshold
            end if
         endif
      endif

! Begin settings relevant for GW et al.

      if(.not. flag_freq_grid_type) then
          freq_grid_type = 1
      endif

      if (use_gw_expt.or.use_gw.or.use_mp2sf) then
         if (.not. flag_anacon_type) then
             write(info_str,'(2X,A)') &
                 "* "
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* Input verification for quasiparticle calculation:"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* Missing keyword in control.in."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* This appears to be a quasiparticle calculation, but the type of analytical continuation"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* of the self-energy from the imaginary to the real axis was not set."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* Please use the 'anacon_type' keyword to set the analytical continuation type"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* and restart your calculation. The possible options are 'two-pole' or 'pade'."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* "
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* FHI-aims used to set this default silently, but the analytical continuation type"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* can have a noticeable effect on the quasiparticle eigenvalues from, e.g., a GW calculation."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* The Pade approximation is known to be more accurate for well-behaved systems, but can show"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* numerical instabilities in some cases. This is, therefore, a choice of which anyone"
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* using GW or the MP2 self-energy should be aware."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* For more information, please see the manual (keyword anacon_type)."
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)') &
                 "* Thank you & sorry for the inconvenience."
             call localorb_info(info_str,use_unit,'(A)')
             call aims_stop_coll('Please set the anacon_type keyword.', func)
         endif
      end if

      if(.not. flag_n_anacon_par) then
          if(anacon_type .eq. 0) then
             n_max_par=5
          elseif(anacon_type .eq. 1) then
             n_max_par=16
          else
             write(info_str,'(2X,A,A)') &
                 "The type of analytical continuation is not correctly set; ", &
                 "please set anacon_type to 'two-pole' or 'pade'."
             call localorb_info(info_str,use_unit,'(A)')
             call aims_stop_coll('', func)
          endif
      endif

      if (.not.flag_frequency_points) then
         if (freq_grid_type.eq.0) then
            if(anacon_type .eq. 0) then
               n_full_freq=80
               n_freq=40
            elseif(anacon_type .eq. 1) then
               n_full_freq=200
               n_freq=100
            else
             write(info_str,'(2X,A,A)') &
                 "The type of analytical continuation is not correctly set; ", &
                 "please set anacon_type to 'two-pole' or 'pade'."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            endif
         else if (freq_grid_type.eq.1) then
            if(anacon_type .eq. 0) then
               n_full_freq=40
               n_freq=40
            elseif(anacon_type .eq. 1) then
               n_full_freq=100
               n_freq=100
            else
               write(info_str,'(2X,A,A)') &
                 "The type of analytical continuation is not correctly set; ", &
                 "please set anacon_type to 'two-pole' or 'pade'."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            endif
         else if (freq_grid_type.eq.2) then
            if(anacon_type .eq. 0) then
               n_full_freq=100
               n_freq=100
            elseif(anacon_type .eq. 1) then
               n_full_freq=200
               n_freq=200
            else
               write(info_str,'(2X,A,A)') &
                 "The type of analytical continuation is not correctly set; ", &
                 "please set anacon_type to 'two-pole' or 'pade'."
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop_coll('', func)
            endif
         endif
         if (use_scgw .or. use_scgw0)then
              n_full_freq=60
         elseif (use_gw_energy)then
              n_full_freq=100
              n_freq=100
         endif
      endif

      if(.not.flag_time_points) then
           n_full_time=80

         if (use_scgw .or. use_scgw0)then
           n_full_time=60
         endif
      endif
      if(use_scgw .or. use_scgw0)then
        if (.not.flag_maximum_frequency) then
          omegamax = 7000
          taumax = 1000
        endif
      else
        if(freq_grid_type == 2)then
           omegamax = 5000
           taumax   = 5000
        endif
      endif

      if(use_scgw .or. use_scgw0) then
         use_scalapack = .false.
         use_2d_corr = .false.
         use_lapack_fast = .true.
         flag_KS = -1
         if(flag_KS_method)then
           write(info_str,'(2X,2A)')'* Self-consistent GW and self-consistent GW0 ', &
           'is currently available only with Lapack.'
           call localorb_info(info_str,use_unit,'(A)')
         endif
      endif

      if(use_gw .and. (n_periodic .gt. 0) .and. (.not. use_gw_expt)) then
             write(info_str,'(2X,A)')'* Periodic GW implementation is still at experimental stage. If you want to do periodic GW calculations, '
             call localorb_info(info_str,use_unit,'(A)')
             write(info_str,'(2X,A)')'* please set explicitly "qpe_calc gw_expt" in the "control.in" file.'
             call localorb_info(info_str,use_unit,'(A)')
             call aims_stop_coll('', func)
      endif

      if(use_gw_expt .and. n_periodic .gt. 0) then
             write(info_str,'(2X,A)')'* Periodic GW calculations will start ...... '
             call localorb_info(info_str,use_unit,'(A)')
             use_periodic_gw = .true.
             use_basis_gradients = .true.
! to ensure the entire spectrum is included when computing the dielectric function along the imaginary axis.
             Emin = -1.e9
             Emax =  1.e9
             n_plot_dielectric = 3
      endif

      !if(use_full_spectrum) then
      !   write(info_str,'(2X,2A)')'* Doing correlated calculations, ', &
      !   'so all empty single-particle states will be included.'
      !   call localorb_info(info_str,use_unit,'(A)')
      !   n_empty = 2000
      !endif
      ! Igor :: For illconditioning problems, it would be good if we can specify
      !         the number of empty states by hand even for correlated
      !         calculations
      !         Note that, in "read_geo.in", 'use_full_spectrum' will be used
      !         again to defined the variable 'n_state=n_max_basis', we don't need to
      !         specify a very large 'n_empty' here.
      if (n_empty .ne. -1) then
          if(use_full_spectrum) then
             write(info_str,'(2X,2A)')'* Doing correlated calculations, ', &
             'so all empty single-particle states will be included.'
             call localorb_info(info_str,use_unit,'(A)')
             n_empty = -1
          endif
      !else
      !    use_full_spectrum = .false.
      endif

! End settings relevant for GW et al.

      if (.not.flag_multipole_threshold) then
         ! Use default
         multipole_threshold = 1.d-10

         if (myid.eq.0) then
            write(use_unit,'(2X,A)') &
            'No q(lm)/r^(l+1) cutoff set for long-range Hartree potential.'
            write(use_unit,'(2X,A,E14.6,A)') &
            '| Using default value of', multipole_threshold, " ."
            write(use_unit,'(2X,A)') &
            '| Verify using the multipole_threshold keyword.'
         end if

      end if

      if (.not. flag_legacy_monopole_extrapolation) then
         if (legacy_monopole_extrapolation) then
            call localorb_info('  Defaulting to legacy monopole extrapolation.')
         else
            call localorb_info('  Defaulting to new monopole extrapolation.')
         end if
      end if

      ! Check l_hartree_far_distance and warn if needed
      if (maxval(l_hartree(:)).gt.l_hartree_far_distance) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A,I5,A)') &
            '* Warning: Maximum requested Hartree potential components are l_hartree = ', maxval(l_hartree(:)), ", "
            write(use_unit,'(1X,A, I5)') &
            '* BUT any Hartree potential components above l = ', l_hartree_far_distance
            write(use_unit,'(1X,A)') &
            '* will be cut off in the analytical (far-field) part of the Hartree potential.'
            write(use_unit,'(1X,A)') &
            '* Use parameter l_hartree_far_distance in control.in to change this behavior - but check for numerical noise.'
         end if

      end if

      ! check relativity related settings
      if (flag_rel==REL_zora .and. compute_forces)then

         if (myid.eq.0) then

            write(use_unit,'(1X,A)') '* Error: No forces available for "relativistic zora".'
            write(use_unit,'(1X,A)') '* For relaxations, please use "relativistic atomic_zora scalar".'

            if ((flag_acc_forces).and.(.not.flag_relax_mode)) then
               write(use_unit,'(1X,A)') '* Perhaps you did not intend to relax anything, '
               write(use_unit,'(1X,A)') '* but set the s.c.f. convergence criterion sc_acc_forces ? '
               write(use_unit,'(1X,A)') '* In that case, FHI-aims must attempt to compute forces anyway,'
               write(use_unit,'(1X,A)') '* do not use with "relativistic zora". '
            end if

            write(use_unit,'(1X,A)') '* Please correct your control.in file before continuing. '

         end if

         call aims_stop_coll('', func)

      end if
!   check  setting related to new relativity threatment
      if (flag_rel==REL_zora_spinor .or.flag_rel==REL_at_zora_spinor) then
      !       check that we have two component in calculations
         if (spin_treatment.eq.0) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') '* Error: Only one spin channel in two- or four-component relativity treatment'
               write(use_unit,'(1X,A)') '* "spin collinear" is required for two- or four-component treatment'
               write(use_unit,'(1X,A)') '* Please correct your control.in file before continuing. '
            end if
            call aims_stop_coll('', func)
         end if

      end if
      if (use_molecular_dynamics) then
         if (.not.flag_wf_func .and. .not.flag_wf_named) then
            if (.not. MD_use_schedule .and. MD_ensemble == 'NVE') then
               ! polynomial 3 1
               wf_extra_type = WF_EXTRA_FUNC
               n_wf_funcs = 3
               call parse_general_function('constant', wf_funcs(1))
               call format_general_function(wf_funcs(1), buffer)
               write(info_str, "(2X,'Added default wf_func ',A)") trim(buffer)
               call localorb_info(info_str)
               call parse_general_function('linear', wf_funcs(2))
               call format_general_function(wf_funcs(2), buffer)
               write(info_str, "(2X,'Added default wf_func ',A)") trim(buffer)
               call localorb_info(info_str)
               call parse_general_function('cubic', wf_funcs(3))
               call format_general_function(wf_funcs(3), buffer)
               write(info_str, "(2X,'Added default wf_func ',A)") trim(buffer)
               call localorb_info(info_str)
            else
               if (n_wf_funcs /= 0) call aims_stop('n_wf_funcs mismatch', func)
               if (wf_extra_type /= WF_EXTRA_NONE) then
                  call aims_stop('why wf_type?', func)
               end if
            end if
         else
            if (wf_extra_type /= WF_EXTRA_NONE .and. &
            &   MD_ensemble == 'NVE_4th_order') then
               call aims_stop('No wf extrapolation with 4th order NVE', func)
            end if
         end if
      end if
      use_wf_extrapolation = (wf_extra_type /= WF_EXTRA_NONE)

!  choose density update method
      select case (density_update_method)
      case(UPDATE_ORBITAL)
         write(info_str,'(2X,A)') &
               'Density update method: Loop over occupied states selected for charge density update.'
         call localorb_info(info_str)
         use_density_matrix = .false.
         if (n_periodic > 0) then
            write(info_str,'(1X,A,A)') &
                  '* Error: periodic system requires density matrix based ', &
                  'density update '
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (wf_extra_use_densmat) then
            write(info_str,'(1X,A,A)') &
                  '* Error: use of keyword wf_extra_use_densmat requires density ', &
                  'matrix based density update'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (use_plus_u) then
            use_density_matrix = .true.
            write(info_str,'(1X,A)') '* Error: DFT+U requires using density matrix based charge density update.'
            call localorb_info(info_str,use_unit, '(A)')
            call aims_stop_coll('', func)
         end if
         if (.not. collect_eigenvectors) then
            write(info_str,'(1X,A,A)') &
                  '* Error: with the orbital-based density update method, ', &
                  'collect_eigenvectors .false. is not possible '
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if

      case(UPDATE_DENS_MATRIX)
         write(info_str,'(2X,A)') &
               'Density update method: density matrix based density update selected.'
         call localorb_info(info_str)
         use_density_matrix = .true.
         if (use_embedding_pp) then
            write(info_str,'(1X,A)') &
                  '* Error: When pseudopotentials are used, only the Kohn-Sham orbital based'
            call localorb_info(info_str)
            write(info_str,'(1X,A,A,A)') &
                  '* density update method is currently allowed.'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (use_hartree_fock .and. use_scalapack .and. (n_periodic.eq.0) .and. (RI_type .ne. RI_LVL )) then
            write(info_str,'(1X,A,A)') '* Only KS-based density update method is allowed ', &
                                       'as Hartree-Fock with RI-V needs packed_matrix_format "none"'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
      case(UPDATE_SPLIT)
         write(info_str,'(2X,A,A)') &
            'Density update method: split updates selected - force will be updated via density matrix ', &
            'and electron density via Kohn Sham orbitals.'
         call localorb_info(info_str)
         use_density_matrix = .false.
         split_updates = .true.
         if (n_periodic > 0) then
            write(info_str,'(1X,A,A)') &
                  '* Error: periodic system requires density matrix based ', &
                  'density update '
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (wf_extra_use_densmat) then
            write(info_str,'(1X,A,A)') &
                  '* Error: use of keyword wf_extra_use_densmat requires density ', &
                  'matrix based density update'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (use_plus_u) then
            use_density_matrix = .true.
            write(info_str,'(1X,A)') '* Error: DFT+U requires using density matrix based charge density update.'
            call localorb_info(info_str,use_unit, '(A)')
            call aims_stop_coll('', func)
         end if
         if (use_embedding_pp) then
            write(info_str,'(1X,A)') &
                  '* Error: When pseudopotentials are used, only the Kohn-Sham orbital based'
            call localorb_info(info_str)
            write(info_str,'(1X,A,A,A)') &
                  '* density update method is currently allowed.'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (use_hartree_fock .and. use_scalapack .and. (n_periodic.eq.0).and. (RI_type .ne. RI_LVL )) then
            write(info_str,'(1X,A,A)') '* Only KS-based density update method is allowed ',&
                                       'as Hartree-Fock with RI-V needs packed_matrix_format "none"'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if (split_updates .and. use_scalapack) then
            packed_matrix_format = PM_index
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') 'Scalapack and split updates requested:'
               write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
            end if
         end if
      case(UPDATE_AUTO)
         write(info_str,'(2X,A)') &
            'Density update method: automatic selection selected.'
         call localorb_info(info_str)
         autoselect_density_method = .true.
         if (.not. (ALL(specified_grid).and.ALL(angular_acc == 0.0d0))) then
            write(info_str,'(1X,A,A)') &
                  '* Error: adapted grids: automatic selection of density update ', &
                  'method not possible'
            call localorb_info(info_str)
            call aims_stop_coll('', func)
         end if
         if ((n_periodic > 0) .or. wf_extra_use_densmat) then
            use_density_matrix = .true.
            autoselect_density_method = .false.
            write(info_str,'(2X,A)') 'Using density matrix based charge density update.'
            call localorb_info(info_str,use_unit, '(A)')
         end if
         if (use_plus_u) then
            use_density_matrix = .true.
            write(info_str,'(1X,A)') '* Warning: DFT+U requires using density matrix based charge density update.'
            call localorb_info(info_str,use_unit, '(A)')
            autoselect_density_method = .false.
         end if
         if (use_hartree_fock .and. use_scalapack .and. (n_periodic.eq.0) .and. (RI_type .ne. RI_LVL )) then
            write(info_str,'(2X,A,A)') 'KS-based density update method is selected as ',&
                                       'as Hartree-Fock with RI-V needs packed_matrix_format "none"'
            call localorb_info(info_str)
            use_density_matrix = .false.
            autoselect_density_method = .false.
         end if
         if (.not. collect_eigenvectors) then
            write(info_str,'(1X,A,A)') &
                  '* Warning: with the KS-based density update method, ', &
                  'collect_eigenvectors .false. is not possible. '
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                  '* The density matrix based density update will be used.'
            call localorb_info(info_str)
            use_density_matrix = .true.
            autoselect_density_method = .false.
         end if
         if (use_embedding_pp) then
            write(info_str,'(1X,A)') &
                  '* Warning: When pseudopotentials are used, only the Kohn-Sham orbital based'
            call localorb_info(info_str)
            write(info_str,'(1X,A,A,A)') &
                  '* density update method is currently allowed. For very large QM systems, this'
            call localorb_info(info_str)
            write(info_str,'(1X,A,A,A)') &
                  '* will lead to O(N^2) scaling, please contact us if this is a practical problem.'
            call localorb_info(info_str)
            use_density_matrix = .false.
            autoselect_density_method = .false.
         end if
         if (flag_packed_matrix_format .and. (packed_matrix_format .ne. PM_index) .and. use_scalapack) then
            write(info_str,'(1X,A,A)') &
                  '* Warning: If the automatic selection of the density update method ', &
                  'selects the density matrix based update'
            call localorb_info(info_str)
            write(info_str,'(1X,A)') &
                  '           or split updates, a packed matrix format ("index") will be used.'
            call localorb_info(info_str)
         end if
      end select

      if (use_elsi_dm) then
         if (.not. flag_density_update_method) then
            density_update_method = UPDATE_DENS_MATRIX
            use_density_matrix = .true.
         end if

         if (n_periodic == 0 .and. density_update_method /= UPDATE_DENS_MATRIX) then
            write(info_str,"(1X,2A)") "* ERROR: 'elsi_method dm' requires",&
               " 'density_update_method density_matrix'."
            call localorb_info(info_str)
            call aims_stop_coll()
         end if

         if (collect_eigenvectors) then
            write(info_str,"(1X,2A)") "* ERROR: 'elsi_method dm' requires",&
               " 'collect_eigenvectors .false.'."
            call localorb_info(info_str)
            call aims_stop_coll()
         end if
      end if

      if (elsi_extrap_dm) then
         if (elsi_solver /= 1) then
            elsi_extrap_dm = .false.
         else if (.not. elsi_write_dm .and. .not. elsi_read_dm) then
            elsi_extrap_dm = .false.
         end if
      end if

      if (wf_extra_use_densmat) then
         use_density_matrix = .true.
         write(info_str,'(2X,A)') 'Using density matrix based charge density update.'
         call localorb_info(info_str,use_unit, '(A)')
      end if

      if (use_plus_u) then
         if(.not.use_density_matrix)then
            ! This will overwrite the embedding statement, this is totally ok!
            use_density_matrix = .true.
            write(info_str,'(1X,A)') '* DFT+U requires using density matrix based charge density update.'
            call localorb_info(info_str,use_unit, '(A)')
            write(info_str,'(1X,A)') '* Switched the update method to the density matrix based update.'
            call localorb_info(info_str,use_unit, '(A)')
         endif
         ! Check for load_balancing first, since it will automatically set use_local_index in some cases
         if(use_load_balancing)then
            write(info_str,'(1X,A)') '* ERROR: DFT+U not yet implemented for load_balancing or use_local_index.'
            call localorb_info(info_str,use_unit, '(A)')
            call aims_stop_coll('',func)
         endif
         if(use_local_index)then
            write(info_str,'(1X,A)') '* ERROR: DFT+U not yet implemented for use_local_index.'
            call localorb_info(info_str,use_unit, '(A)')
            call aims_stop_coll('',func)
         endif
      end if

! Periodic boundaries check

      if (n_periodic > 0)then

         use_density_matrix = .true.
         if (myid.eq.0) then
            write(use_unit,'(2X,A)') &
            'Using density matrix based charge density update.'
         end if

         if ( .not. flag_packed_matrix_format)then
            packed_matrix_format = PM_index
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') &
               'Using packed matrix style: index .'
            end if
         end if

         if (.not. flag_k_points_defined)then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
               '* Error: k-grid not specified!'
            end if
            call aims_stop_coll('', func)
         end if

         if (use_prodbas) then
            if (use_corr .and. (.not. use_rpa_ene) .and. &
                (.not. use_gw) .and. (.not. use_ci) .and. (.not. use_mp2) &
                .and. (.not. use_dftpt2) .and. (.not. use_os_mp2)) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                     '* Error. You are attempting to run a periodic structure with a correlation treatment'
                  write(use_unit,'(1X,A)') &
                     '* which requires a two-electron Coulomb operator (e.g., GW, MP2 et al.), '
                  write(use_unit,'(1X,A)') &
                     '*  We are actively working on all these for periodic boundary conditions, but the first '
                  write(use_unit,'(1X,A)') &
                     '* implementation is not yet finished. Apologies for the inconvenience.'
               end if
               call aims_stop_coll('', func)
            else if(use_mp2 .or. use_dftpt2 .or. use_os_mp2) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                     '* WARNING: This is an experimental version of the periodic PT2 implemenation.'
                  write(use_unit,'(1X,A)') &
                     '* We are still actively working on it to improve the efficiency and reliability. '
                  write(use_unit,'(1X,A)') &
                     '* Do not rely on any results you obtained from periodic PT2 calculation.'
               endif
            else if(use_rpa_ene) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                     '* WARNING: This is an experimental version of the periodic RPA implemenation.'
                  write(use_unit,'(1X,A)') &
                     '* We are still actively working on it to improve the efficiency and reliability. '
                  write(use_unit,'(1X,A)') &
                     '* Do not rely on any results you obtained from periodic RPA calculation.'
               endif
            else if(use_gw) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                     '* WARNING: This is an experimental version of the periodic GW implemenation.'
                  write(use_unit,'(1X,A)') &
                     '* We are still actively working on it to improve the efficiency and reliability. '
                  write(use_unit,'(1X,A)') &
                     '* Do not rely on any results you obtained from periodic GW calculation.'
               endif
            else if (use_periodic_hf) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') &
                  '* NOTE: You are attempting to run a periodic structure with an exchange treatment'
                  write(use_unit,'(1X,A)') &
                  '* which requires a two-electron Coulomb operator (Hartree-Fock or a hybrid functional).'
                  write(use_unit,'(1X,A)') &
                  '* The present implementation is reliable to our knowledge, but some things'
                  write(use_unit,'(1X,A)') &
                  '* (like band structure output, or perhaps the auxilary basis set for small'
                  write(use_unit,'(1X,A)') &
                  '* orbital basis sets - see the "for_aux" keyword) require extra care.'
                  write(use_unit,'(1X,A)') '* Please let us know if you encounter any trouble.'
               endif
            endif
         end if

      else  ! Cluster case

         if ( flag_k_points_defined) then
            if (myid.eq.0) then
            write(use_unit,'(1X,A)') &
               '* Error: Cluster calculation cannot have k-grid!'
            end if
            call aims_stop_coll('', func)
         end if

         if (out_band) then
            out_band = .false.
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') '* Warning: Requesting band structure output in', &
                  'non-periodic systems has no effect.'
            end if
         end if

         if (out_band_during_scf) then
            out_band_during_scf = .false.
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') '* Warning: Requesting band structure output in', &
                  'non-periodic systems has no effect.'
            end if
         end if

         ! CC: Z2 by CC / Adapted from Cmera
         if (out_z2_invariant) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') '* Warning: Requesting Z2 invariant output in', &
                  'non-periodic systems has no effect.'
            end if
         end if


         if ( (use_density_matrix) .and. (use_scalapack) ) then

            if ( .not. flag_packed_matrix_format)then
               packed_matrix_format = PM_index
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') 'Scalapack and density-matrix based density update requested:'
                  write(use_unit,'(2X,A)') 'Switching to packed matrix format "index" by default.'
               end if
            else if (packed_matrix_format .ne. PM_index) then
               if (myid.eq.0) then
                  write(use_unit,'(2X,A,A,A)') &
                  'Error: Must use packed matrix format when ', &
                  'density-matrix based density update and ', &
                  'Scalapack are both requested.'
               end if
               call aims_stop_coll('', func)
            end if

            if (use_plus_u .and. use_embedding_pp ) then
                if (myid.eq.0) then
                   write(use_unit,'(2X,A)') 'packed_matrix_format was set to PM_index'
                   write(use_unit,'(2X,A)') 'This does not work with the QM/MM embedding infrastructure'
                end if
                call aims_stop_coll('', func)
            endif

         end if

         if (use_density_matrix .and. compute_forces &
            .and. .not. force_new_functional)then

            if (myid.eq.0) then
               write(use_unit,*) 'Error: cluster forces work with', &
                  'density matrix only with new energy functional'
            end if
            call aims_stop_coll('', func)

         end if

      end if


      if (use_periodic_hf) then
        ! print out any unmodified default parameters and set necessary defaults
        if (.not.flag_crit_val) then
          if (myid.eq.0) then
            write (use_unit,'(2X,A,E14.6)') &
              "Coulomb operator: Defaulting to screening threshold for four-center integrals: ", crit_val
          end if
        end if
        if (.not.flag_coul_mat_threshold) then
          if (myid.eq.0) then
            write (use_unit,'(2X,A,E14.6)') &
              "Coulomb operator: Defaulting to screening threshold for Coulomb matrix: ", coul_mat_threshold
          end if
        end if

        if (.not.flag_RI) then
           ! Periodic Hartree-Fock only really works with LVL RI, for now
           RI_type = RI_LVL
          if (myid.eq.0) then
            write (use_unit,'(2X,A)') &
              "Coulomb operator: Defaulting to 'LVL_fast' resolution of identity - RI_type = LVL_fast "
            write (use_unit,'(2X,A)') &
              "No other resolution of identity version is supported for periodic boundary conditions at present."
          end if
        else if (RI_type .ne. RI_LVL) then
          if (myid.eq.0) then
            write (use_unit,'(1X,A,I5)') &
              "* Error: Coulomb operator. Resolution of identity is not LVL - RI_type = ", RI_type
            write (use_unit,'(1X,A)') &
              "* At present, only RI_method LVL is supported for production calculations."
            write (use_unit,'(1X,A)') &
              "* If you are actively working on another method, please disable the following 'stop' in the source code"
            write (use_unit,'(1X,A)') &
              "* for your own version. We stop here to prevent unpredictable results for anyone else."
            write (use_unit,'(1X,A)') &
              "* Apologies for the inconvenience."
          end if
          call aims_stop_coll('Unsupported RI_method choice for periodic Hartree-Fock exchange', func)
        end if

        if (.not.flag_use_logsbt) then
           ! Periodic Hartree-Fock should use logsbt integrations
           use_logsbt = .true.
           if (myid.eq.0) then
             write (use_unit,'(2X,A)') &
               "Coulomb operator: Defaulting to logarithmic spherical Bessel transform integrals (use_logsbt .true.)."
           end if
        else if (.not.use_logsbt) then
          if (myid.eq.0) then
            write (use_unit,'(1X,A,I5)') &
              "* Error: use_logsbt .false. is not officially supported for RI_method LVL_fast."
          end if
          call aims_stop_coll('RI_method LVL_fast needs logsbt integrals (use_logsbt .true.).', func)
        end if

        if (use_hse .and. .not. use_lc_wpbeh) then
          ! periodic HSE wants logsbt as well
          if (.not.flag_use_logsbt_for_radial_hse_integration) then
             ! Periodic Hartree-Fock should use logsbt integrations
             use_logsbt_for_radial_hse_integration = .true.
             if (myid.eq.0) then
               write (use_unit,'(2X,A)') &
                 "HSE with LVL_fast: Also defaulting to logarithmic spherical Bessel transform integrals."
             end if
          else if (.not.use_logsbt_for_radial_hse_integration) then
            if (myid.eq.0) then
              write (use_unit,'(1X,A,I5)') &
                "* Error: use_logsbt_for_radial_hse_integration .true. required for HSE / LVL_fast."
            end if
            call aims_stop_coll('LVL_fast: HSE needs logsbt integrals.', func)
          end if
        end if

        if (use_lc_wpbeh) then
          ! periodic lc_wpbeh wants logsbt as well
          if (.not.flag_use_logsbt_for_radial_hse_integration) then
             ! Periodic Hartree-Fock should use logsbt integrations
             use_logsbt_for_radial_hse_integration = .true.
             if (myid.eq.0) then
               write (use_unit,'(2X,A)') &
                 "LC-wPBEh with LVL_fast: Also defaulting to logarithmic spherical Bessel transform integrals."
             end if
          else if (.not.use_logsbt_for_radial_hse_integration) then
            if (myid.eq.0) then
              write (use_unit,'(1X,A,I5)') &
                "* Error: use_logsbt_for_radial_hse_integration .true. required for LC-wPBEh / LVL_fast."
            end if
            call aims_stop_coll('LVL_fast: LC-wPBEh needs logsbt integrals.', func)
          end if
        end if

      end if  ! end if use_periodic_hf

      if (store_eigenvectors_to_disk_in_relaxation .and. use_scalapack) then
         call aims_stop_coll('Option store_EV_to_disk_in_relaxation not allowed with scalapack.', func)
      end if
      if ( (.not. collect_eigenvectors) .and. (.not. use_scalapack) )then
         if (myid.eq.0) then
            write(use_unit,'(1X,A)') '* Warning: Since this is not a scalapack run, collect_eigenvectors .false. is not possible.'
            write(use_unit,'(1X,A)') '* Eigenvectors will be collected.'
         end if
         collect_eigenvectors = .true.
      end if

      if ((.not. collect_eigenvectors).and.(.not.use_density_matrix))then
         if (myid.eq.0) then
            write(use_unit,'(1X,2A)') &
               '* Warning: with the KS-based density update method, ', &
               'collect_eigenvectors .false. is not possible.'
            write(use_unit,'(1X,A)') '* Eigenvectors should be collected. Aborting the present run.'
            write(use_unit,'(1X,A)') '* Please EITHER set collect_eigenvectors .true. in control.in (if you have enough memory!),'
            write(use_unit,'(1X,A)') '*        OR     use the density matrix based density update (especially if you lack memory!!)'
         end if
         call aims_stop_coll('', func)
      end if


      if (use_local_index .or. use_load_balancing) then
         ! WPH, 2018-01-22: I've added "and/or load_balancing" to all output messages
         ! here, since the load_balancing keyword now automatically sets use_local_index
         ! unless the user has set it otherwise.  Since load_balancing relies on
         ! use_local_index, any error message for the latter is automatically an error
         ! message for the former.


         ! TZ, 2018-02-01: adding check-inconsistency here

         if (flag_out_dielectric) then
            call aims_stop_coll ('* Dielectric calculations currently do not support load_balancing or &
                                                &use_local_index keywords.Exiting', func)
         endif

         if (out_cube) then
            do cube_id = 1, i_cube, 1
               if (cube_type_needs_densmat(cube_id)) then
                  call localorb_info("")
                  write(info_str, "(2X,A,' ',A,' ',A)") &
                     & 'Cube type', trim(cube_type(cube_id)), &
                     & 'not implemented for load_balancing or use_local_index keywords.  Exiting.'
                  call aims_stop_coll(info_str, func)
               endif
            enddo
         endif

         if (out_esp) then
             call aims_stop_coll ('* ESP charges calculations currently do not support load_balancing or &
                                                &use_local_index keywords.Exiting', func)
         endif

         if (out_matrices) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing ', &
               'does not work with out_matrices.'
            endif
            call aims_stop_coll('', func)
         endif

         if (use_constraint) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing ', &
               'does not work with use_constraint'
            endif
            call aims_stop_coll('', func)
         endif

         if (use_hartree_fock .and. (.not. use_periodic_hf)) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing ', &
                    'works only with the LVL resolution of identity method.'
            endif
            call aims_stop_coll('', func)
         endif

         if (.not. use_scalapack) then
            if (myid.eq.0) then
               write(use_unit,*) ''
               write(use_unit,*) 'Error: use_local_index and/or load_balancing ', &
               'requires a parallel eigensolver to be used.  You are seeing this'
               write(use_unit,*) 'message because you are ', &
               'running with too few MPI tasks to enable a parallel eigensolver.'
               write(use_unit,*) 'KS_method parallel will be enabled automatically', &
               'if you are using more MPI tasks than S.C.F. k-points.'
               write(use_unit,*) 'Either use more MPI tasks ', &
               'or omit the use_local_index and load_balancing keywords.'
            endif
            call aims_stop_coll('', func)
         endif

         if (packed_matrix_format /= PM_index) then
            if (.not.flag_packed_matrix_format) then
            packed_matrix_format = PM_index
            if (myid.eq.0) then
               write(use_unit,'(1X,A,A)') '* Warning: use_local_index and/or load_balancing', &
               ' needs packed_matrix_format = index'
               write(use_unit,'(1X,A,A)') '* Changing the default setting to ', &
               'packed_matrix_format = index'
            endif
            else
               if (myid.eq.0) then
                  write(use_unit,*) 'Error: use_local_index and/or load_balancing ', &
                  'needs packed_matrix_format = index'
               endif
               call aims_stop_coll('', func)
            end if
         endif

         if (sum(angular_acc).ne.0.d0) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing works only with fixed grids'
            endif
            call aims_stop_coll('', func)
         end if
      endif

!     Check compatibility of 'output aitranss' with other flags
      if(out_aitranss) then
         if(n_periodic > 0) then
            out_matrices_elsi = .true.

            write(info_str,'(2X,A)') &
               '*  WARNING: Output of the overlap matrix and KS-eigenvectors was requested'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  in a format compatible with the "AITRANSS" code. This option conflicts'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  with the present periodic geometry. The overlap matrix will be written'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  to an ELSI CSC file.'
            call localorb_info(info_str,use_unit,'(A)')
         endif

         if(use_scalapack .and. (.not.collect_eigenvectors)) then
            collect_eigenvectors = .true.

            write(info_str,'(2X,A)') '* '
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  WARNING: Output of the overlap matrix and KS-eigenvectors was requested'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  in a format compatible with the "AITRANSS" code. This requires'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  "collect_eigenvectors .true.", but ".false." was explicitly requested.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  Overriding this choice and collecting eigenvectors after all.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  However, for very large systems, this change may lead to memory problems.'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  Please inform us if this becomes an issue for you.'
            call localorb_info(info_str,use_unit,'(A)')
         endif

         if(use_scalapack .and. (packed_matrix_format /= PM_none)) then
            out_matrices_elsi = .true.

            write(info_str,'(2X,A)') &
               '*  WARNING: Output of the overlap matrix and KS-eigenvectors was requested'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  in a format compatible with the "AITRANSS" code. This option conflicts'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  with a usage of the packed matrix format, which has been explicitly'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  requested. The overlap matrix will be written to an ELSI CSC file.'
            call localorb_info(info_str,use_unit,'(A)')
         endif

         if(use_local_index.or.use_load_balancing) then
            out_matrices_elsi = .true.

            write(info_str,'(2X,A)') &
               '*  WARNING: Output of the overlap matrix and KS-eigenvectors was requested'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  in a format compatible with the "AITRANSS" code. This option conflicts'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  with "use_local_index" and/or "load_balancing", which have been explicitly'
            call localorb_info(info_str,use_unit,'(A)')
            write(info_str,'(2X,A)') &
               '*  requested.  The overlap matrix will be written to an ELSI CSC file.'
            call localorb_info(info_str,use_unit,'(A)')
         endif
      endif ! out_aitranss

      if (restart_relaxations) then
         if (flag_write_restart_geometry) then
            if (.not.write_restart_geometry) then
               write(info_str,'(1X,2A)') &
               '* Warning: "restart_relaxations" requested, but "write_restart_geometry" is explicitly false. '
               call localorb_info(info_str)
               write(info_str,'(1X,2A)') &
               '* This sounds like conflicting settings - please check. '
               call localorb_info(info_str)
               call aims_stop_coll('', func)
            end if
         else
            write_restart_geometry = .true.
         end if
      end if

      if (((relax_mode.eq.RELAX_TRM) .and. (relax_geo.eq.RELAX_TRM)) .or. &
        ((relax_mode.eq.RELAX_LTRM) .and. (relax_geo.eq.RELAX_LTRM))) then
         if (.not.flag_write_restart_geometry) then
            write_restart_geometry = .true.
            write(info_str,'(2X,2A)') &
              'Geometry relaxation: A file "geometry.in.next_step" is written out by default after each step.'
            call localorb_info(info_str)
            write(info_str,'(2X,2A)') &
              '| This file contains the geometry of the current relaxation step as well as'
            call localorb_info(info_str)
            write(info_str,'(2X,2A)') &
              '| the current Hessian matrix needed to restart the relaxation algorithm.'
            call localorb_info(info_str)
            write(info_str,'(2X,2A)') &
              '| If you do not want part or all of this information, use the keywords'
            call localorb_info(info_str)
            write(info_str,'(2X,2A)') &
              '| "write_restart_geometry" or "hessian_to_restart_geometry" to switch the output off.'
            call localorb_info(info_str)
         end if
      end if

      if (write_restart_geometry) then
         if (.not.flag_hessian_to_restart_geometry) then
            hessian_to_restart_geometry = .true.
         end if
      end if

      if (use_local_index) then

         if (.not. prune_basis_once) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing needs prune_basis_once'
            endif
            call aims_stop_coll('', func)
         endif

         if (.not. parallel_grid_partitioning) then
            if (myid.eq.0) then
               write(use_unit,*) 'Error: use_local_index and/or load_balancing needs parallel_grid_partitioning'
            endif
            call aims_stop_coll('', func)
         endif

         if (restart_read .or. restart_write) then
            if (myid.eq.0) then
               call aims_stop_coll('Restart is not supported with use_local_index and/or load_balancing.', func)
            end if
            restart_read = .false.
            restart_write = .false.
         end if

      endif

      if(use_load_balancing .and. .not.use_local_index) then
         ! WPH, 2018-01-22:  Code now stops if load balancing is enabled without domain decomposition.
         call aims_stop_coll('Load balancing requires that use_local_index be enabled.  Exiting.', func)
      endif
      if(use_load_balancing .and. .not.use_mpi) then
         call localorb_info('** Turning off load_balancing and use_local_index because we are lacking MPI.')
         use_local_index = .false.
         use_load_balancing = .false.
      endif
      if(use_local_index .and. .not.use_mpi) then
         call localorb_info('** Turning off use_local_index because we are lacking MPI.')
         use_local_index = .false.
      endif

      if (condition_penalty > 0.d0) then
         if (use_load_balancing) then
            call aims_stop_coll('No condition_penalty with load_balancing or use_local_index', func)
         end if
         if (use_local_index) then
            call aims_stop_coll('No condition_penalty with use_local_index', func)
         end if
         if (.not. collect_eigenvectors) then
            call aims_stop_coll('No condition_penalty without local eigenvecs', func)
         end if
      endif

      if (far_distance_hartree_radius_threshold < 1e-15)then
         if (myid.eq.0) then
            write(use_unit,*) 'Error: hartree_radius_threshold can not be smaller than 1e-15'
         end if
         call aims_stop_coll('', func)
      end if

      if(use_split_xc_gw .and. use_gw &
           .and. (flag_xc .ne. 1 .and. flag_xc .ne. 6 .and. flag_xc .ne.8 .and. flag_xc .ne. 23) ) then
         if (myid.eq.0) then
            write(use_unit,*) '* Error: Only for LDA, PBE, PBE0 or LC-wPBEh reference, one can separate '
            write(use_unit,*) '*    the DFT exchange and correlation contributions in G0W0 calculations.'
         end if
         call aims_stop_coll('', func)
      endif

      if (n_periodic > 0) then
         if (.not. flag_symmetry_reduced_k_grid) then
            if (use_prodbas) then
               ! Does not work?
               use_symmetry_reduced_k_grid = .false.
            end if
            ! Otherwise, stick to the default.
            if (use_symmetry_reduced_k_grid) then
               call localorb_info('  Defaulting to use time-reversal symmetry for k-point grid.')
            else
               call localorb_info('  Defaulting to use full k-point grid (no symmetry reduction).')
            end if
         else
            if (use_symmetry_reduced_k_grid) then
               if (use_prodbas) then
                  call localorb_info('Testing: Use Inversion Symmetry for ProdBas. Please report errors.')
               end if
            end if
         end if
      end if

      if ((out_l_proj_dos.or.out_atom_dos.or.out_l_proj_dos_tetrahedron.or.out_atom_dos_tetrahedron)) then
         if (out_lowdin) then
            out_mulliken = .false.
            flag_run_mulliken = .false.
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') 'The projected DOS plots are based on Loewdin analysis.'
            endif
         else if (out_mulliken) then
            flag_run_lowdin = .false.
            if (myid.eq.0) then
               write(use_unit,'(2X,A)') 'The projected DOS plots are based on Mulliken analysis.'
            endif
         else
           if (n_periodic==0) then !Disabled for periodic systems, since this file can become disturbingly large.
               if (myid.eq.0) then
                  write(use_unit,'(2X,A)') 'Did not ask for Mulliken analysis, but requiring it anyway due to projected DOS plots.'
                  !write(use_unit,'(2X,A)') '| Enabling the output of the Mulliken analysis for completeness anyway.'
               end if
               flag_run_mulliken = .true.
              endif
         endif
      end if

      if (flag_run_mulliken .and. use_scalapack .and. (.not. collect_eigenvectors) .and. packed_matrix_format /= PM_index)then
         if (myid.eq.0) then
            write(use_unit,*) 'Error: the combination: Mulliken analysis + scalapack + collect_eigenvectors .false.'
            write(use_unit,*) '       works only with packed_matrix_format index'
         end if
         call aims_stop_coll('', func)
      end if

      if (out_lowdin .and. n_periodic > 0) then
         if (myid.eq.0) then
            write(use_unit,*) 'Error: Loewdin analysis is not yet working for periodic system.'
         endif
         call aims_stop_coll('', func)
      endif


      if ( (n_periodic.gt.0) .and. use_dipole_correction .and. (charge.ne.0.d0) ) then
         if (myid.eq.0) then
            write(use_unit,'(1X,A)') "* Attention. You are trying to use the dipole correction for surface"
            write(use_unit,'(1X,A)') "* (slab) systems. While there is no formal reason that this should not"
            write(use_unit,'(1X,A)') "* work, there is a practical reason why some extra care is needed."
            write(use_unit,'(1X,A)') "* The Coulomb potential from a charged surface reaches far into the"
            write(use_unit,'(1X,A)') "* vacuum. If you are not careful, your dipole correction will not"
            write(use_unit,'(1X,A)') "* correct the surface dipole, but rather depend randomly on the "
            write(use_unit,'(1X,A)') "* chosen location of the vacuum level (set_vacuum_level) in the slab."
            write(use_unit,'(1X,A)') "* If you know what you are doing, you can simply disable the following"
            write(use_unit,'(1X,A)') "* 'stop' command in read_control.f90 and recompile the code."
            write(use_unit,'(1X,A)') "* In any case, we advise you to first plot the long-range part"
            write(use_unit,'(1X,A)') "* of the Coulomb potential in the vacuum layer before continuing."
            write(use_unit,'(1X,A)') "* For now, code execution will be stopped in order to alert you to."
            write(use_unit,'(1X,A)') "* the problem. Apologies for the inconvenience."
         end if
         !call aims_stop
      end if

      flag_acc_potjump = flag_acc_potjump .and. use_dipole_correction  !Simple ignore potential jump convergence citeria if it's not even calculated.

      if (.not.flag_compensate_multipole_errors) then
         ! If there was no specific setting in control.in, multipole errors should be always compensated.
         ! Except a DFPT calculation is reuqested.
         if (use_DFPT.or.use_DFPT_reduce_memory & 
             .or.use_DFPT_polarizability.or.use_DFPT_polar_reduce_memory &
             .or.use_DFPT_dielectric.or.use_DFPT_phonon_gamma.or.use_DFPT_phonon &
             .or.use_DFPT_phonon_reduce_memory.or.use_friction.or.magnetic_response) then

            if (myid.eq.0) then
                  write(use_unit,'(2X,A)') 'WARNING: Switching off the compensation of charge integration errors'
                  write(use_unit,'(2X,A)') 'on the 3D integration grid, since you have requested a DFPT calculation.'
                  write(use_unit,'(2X,A)') 'This is unfortunately necessary since "compensate_multipole_errors" is'
                  write(use_unit,'(2X,A)') 'not yet implemented for DFPT.'
                  write(use_unit,'(2X,A)') 'Use the "compensate_multipole_errors" flag to overwrite this behaviour.'
            end if

           compensate_multipole_errors = .false.

         else
            if (myid.eq.0) then
                  write(use_unit,'(2X,A)') 'Charge integration errors on the 3D integration grid will be compensated'
                  write(use_unit,'(2X,A)') 'by explicit normalization and distribution of residual charges.'
                  write(use_unit,'(2X,A)') 'Use the "compensate_multipole_errors" flag to change this behaviour.'
            end if

           compensate_multipole_errors = .true.

          end if
      end if

      if ((use_hartree_fock .or. use_hf_post ).and. sparse_o3fn) then
         if (hf_version == HF_DM) then
            call localorb_info('  LVL not implemented with density matrix yet')
            call localorb_info('  -> switching to transformed overlap')
            hf_version = HF_EIGEN
         end if
         if (.not. flag_use_logsbt) then
            use_logsbt = .true.
            call localorb_info('  LVL probably works better with logsbt')
         elseif(.not. use_logsbt) then
            call localorb_info('  * Be advised that LVL probably works better with use_logsbt .true.')
         end if
         if (RI_type == RI_LVL_2nd .and. use_corr) then
            call aims_stop('RI_LVL_2nd not implemented for correlated methods', func)
         end if
      end if


      ! One last error trap, independent of all others:
      ! hydro_cut .false. (testing only!) requires a larger free-atom cutoff radius ...
      ! If a product basis is used, the cutoff radius is infinite anyway ... else,
      ! use the same default as for Gaussian basis functions.
      if ((.not.hydro_cut).and.(.not.use_prodbas)) then
         do i_species = 1, n_species,1
            ! if there are any hydrogenic functions defined for this species ...
            if (n_hydro(i_species).gt.0) then
               if (.not.flag_cut_free_atom(i_species)) then

                  cut_free_atom(i_species) = .true.
                  free_r_cut(i_species) = 10.d0 / bohr  ! default: 10 A

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A,A,A)') &
                     "Species ", trim(species_name(i_species)), ": ", &
                     "Analytical hydrogen-like functions requested (hydro_cut = .false.)!"
                     write(use_unit,'(2X,A)') &
                     "Setting particularly large cutoff onset for free atom density."
                     write(use_unit,'(2X,A)') &
                     "Nonetheless, verify explicitly that this is enough if you suspect any "
                     write(use_unit,'(2X,A)') &
                     "integration inaccuracies."
                     write(use_unit,'(2X,A,E15.8,A)') &
                     "Default cutoff onset for free atom density etc. : ", &
                     free_r_cut(i_species)*bohr, &
                     " AA."
                  end if

               else

                  if (myid.eq.0) then
                     write(use_unit,'(2X,A,A,A,A)') &
                     "Species ", trim(species_name(i_species)), ": ", &
                     "Analytical hydrogen-like functions requested (hydro_cut = .false.)!"
                     write(use_unit,'(2X,A)') &
                     "Free-atom cutoff radius was set explicitly for this species, but "
                     write(use_unit,'(2X,A)') &
                     "should be (near) infinite when hydrogen-like functions have no "
                     write(use_unit,'(2X,A)') &
                     "explicit cutoff."
                     write(use_unit,'(1X,A)') &
                     "* Since the current settings were explicitly requested, the calculation will run, "
                     write(use_unit,'(1X,A)') &
                     "* but be sure to check any results VERY carefully and test your cut_free_atom settings."
                  end if

               end if
            end if
         enddo
      end if

      if (use_prodbas) then
         if (use_gw .and. use_mp2sf) then
            call aims_stop_coll('Cannot do both GW and MP2 spectra', func)
         end if
      end if

      if (use_scalapack .and. use_density_matrix)then

         call bcast_logical( keep_restart_info, 0)
         call bcast_logical( restart_read, 0)
         call bcast_logical( restart_write, 0)
         call bcast_logical(restart_file_exists,0)


      end if


      if (.not. flag_use_logsbt .and. use_prodbas) then
         if (use_logsbt) then
            write(info_str,'(2X,A)') &
            'Default to 1D ("use_logsbt") integrations for auxiliary 2-center integrals.'
            call localorb_info(info_str)
         else
            write(info_str,'(2X,A)') &
            'Default to 3D grid integrations for auxiliary 2-center integrals.'
            call localorb_info(info_str)
         end if
      end if
      if (use_logsbt .and. use_prodbas) then
         if (.not. flag_sbtgrid_lnr0) then
            write(info_str,'(2X,A,F16.12)') &
            & 'Default onset of logarithmic r-grid for SBT is', sbtgrid_lnr0
            call localorb_info(info_str)
         end if
         if (.not. flag_sbtgrid_lnk0) then
            write(info_str,'(2X,A,F16.12)') &
            & 'Default onset of logarithmic k-grid for SBT is', sbtgrid_lnk0
            call localorb_info(info_str)
         end if
         if (.not. flag_sbtgrid_lnrange) then
            write(info_str,'(2X,A,F16.12)') &
            & 'Default range of logarithmic r- and k-grid for SBT is', &
            & sbtgrid_lnrange
            call localorb_info(info_str)
         end if
         if (.not. flag_sbtgrid_N) then
            write(info_str,'(2X,A,I6)') &
            & 'Default number of logarithmic r- and k-grid for SBT is', &
            & sbtgrid_N
            call localorb_info(info_str)
         end if
      end if
      if (.not. flag_use_logsbt_for_radial_hse_integration .and. &
      &   (use_hse .and. hse_omega_hf /= 0.d0)) then
         ! Make a consistent choice.
         if (use_logsbt) then
            use_logsbt_for_radial_hse_integration = .true.
            call localorb_info('  * Use logSBT for radial HSE integration, too')
         else
            use_logsbt_for_radial_hse_integration = .false.
            call localorb_info('  Use old radial HSE integrator.')
         end if
      end if
      if (.not. flag_use_2d_corr) then
         use_2d_corr = use_scalapack
         if (use_scalapack .and. use_prodbas) then
            if (use_2d_corr) then
               write(info_str,'(2X,2A)') &
               & 'Default of use_2d_corr is to use the efficient ', &
               & '2D distribution where possible'
            else
               write(info_str,'(2X,A)') &
               & 'Default of use_2d_corr is to use the old 1D distribution'
            end if
            call localorb_info(info_str)
         end if
      else
         if (use_2d_corr .and. (use_scgw.or.use_scgw0)) then
            write(info_str,'(2X,A)') &
            & 'Please note that SCGW still uses the old 1D distribution.'
            call localorb_info(info_str)
         end if
         if (use_2d_corr .and. .not. use_scalapack) then
            write(info_str,'(2X,A)') &
            & '*** The 2D distribution only works with scalapack; disabling.'
            call localorb_info(info_str)
            use_2d_corr = .false.
         end if
      end if
      if((use_contour_def_gw) .and. use_mpi) then
         !** check if scalapack is compiled in
         call BLACS_Pinfo(index1,index2)
         if (index2.gt.0) then
            use_scalapack = .true.
            if(.not.use_2d_corr) then
              use_2d_corr = .true.
              write(info_str,'(2X,A)') &
              & '*** Contour Deformation: Only support 2D parallization. Switching.'
              call localorb_info(info_str,use_unit,'(A)')
            endif
         else
            call aims_stop('Only scalapack support', func)
         endif
      endif
      if ((use_scgw .or. use_scgw) .and. sparse_o3fn) then
         call aims_stop('SCGW with RI-LVL not implemented', func)
      end if
      if(use_2d_corr .and. use_ovlp_swap) then
         call aims_stop('use_ovlp_swap incompatible with use_2d', func)
      endif

      if (force_n_electrons .and. n_periodic.eq.0) then
         write (info_str,'(2X,A)') &
               "force_n_electrons not supported for cluster calculations."
         call localorb_info(info_str,use_unit,'(A)')
         call aims_stop_coll('', func)
      end if

      if (force_occupation_projector .and. packed_matrix_format /= PM_index &
         .and. n_periodic.ne.0) then
               write (info_str,'(2X,A)') &
            "force_occupation_projector needs packed_matrix_format index in periodic case."
         call localorb_info(info_str,use_unit,'(A)')
         call aims_stop_coll('', func)
      end if

      if (force_occupation_projector .and. use_scalapack) then
         write (info_str,'(2X,A)') &
            "*** WARNING: force_occupation_projector currently only works with lapack."
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(2X,A)')"*** eigenvalue solver. Please fix. Aborting."
         call localorb_info(info_str,use_unit,'(A)')
!         call aims_stop_coll('', func)
      end if

!      if (force_occupation_projector .and. .not. restart_file_exists) then
      if ((force_occupation_projector .and. &
          .not. restart_file_exists) .and. &
          myid.eq.0) then
         write (info_str,'(2X,A)') &
            "*** WARNING: force_occupation_projector requires existing restart files."
         call localorb_info(info_str,use_unit,'(A)')
         write (info_str,'(2X,A)')"*** Please provide some. Aborting."
         call localorb_info(info_str,use_unit,'(A)')
         call aims_stop('', func)
      end if

      if (flag_transport)then
         call  transport_check_control_in
      end if


      if (output_level .eq. 'full') output_priority = 0
!     else if (flag_molecular_dynamics .and. (.not.flag_output_level)) then
!         output_level = 'MD_light'
!         output_priority = 2
!      end if
      if (out_t_plus_v_matrix .and. n_periodic.gt.0) then
         write(use_unit,'(1X,2A)') &
           "* Error! You requested the output of t_plus_v matrix (Hcore)."
         write(use_unit,'(1X,A)') &
           "* This feature is not implemented for periodic boundary conditions."
        call aims_stop_coll('', func)
      end if

      if (out_nuc_matrix) then
         if (n_periodic.gt.0) then
            if (myid.eq.0) then
               write(use_unit,'(1X,2A)') &
                  "* Error! You requested the output of the bare nuclear ", &
                  "potential matrix."
               write(use_unit,'(1X,A)') "* This feature is not implemented for periodic boundary conditions."
               write(use_unit,'(1X,2A)') &
                  "* It could be done using an Ewald compensating charge ", &
                  "background, but you will have to notify us."
               write(use_unit,'(1X,A)') "* Stopping for now to alert you - sorry! VB"
            end if
            call aims_stop_coll('', func)
         elseif(packed_matrix_format .ne. PM_none) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "* Error! You requested the output of the bare nuclear potential matrix."
               write(use_unit,'(1X,A)') "* This feature is not implemented for packed matrix storage."
               write(use_unit,'(1X,A)') "* It could be done with relative ease if really needed, but you will have to notify us."
               write(use_unit,'(1X,A)') "* Stopping for now to alert you - sorry! VB"
            end if
            call aims_stop_coll('', func)
         end if
      end if

      ! Set decent defaults if MD_clean_rotations nit specified
      if (.not.flag_MB_clean_rotations_read) then
         if (n_periodic .gt. 0) then
            MB_clean_rotations = .false.
         else
            MB_clean_rotations = .true.
         end if
         if (use_relaxation_constraints) then
            MB_clean_rotations = .false.
         end if
      end if

      ! set zero C6 coefficients for C6 corrections for the relevant pairs if necessary
      if (n_vdw_pairs_ignore .gt. 0) then
         call build_C6_prefactor_matrix()
      end if

      if (use_smooth_prodbas_threshold) then
         if (.not. use_asym_inv_sqrt) then
            call aims_stop('use_smooth_prodbas_threshold only with use_asym_inv_sqrt')
         end if
      end if
      if (use_asym_inv_sqrt) then
         if (.not. use_prodbas) then
            call aims_stop('No prodbas for use_asym_inv_sqrt')
         elseif(RI_type /= RI_V) then
            call aims_stop('use_asym_inv_sqrt only for RI-V')
         elseif(.not. use_scalapack) then
            call aims_stop('use_asym_inv_sqrt only for use_scalpack')
         elseif(use_logsbt) then
            call localorb_info('** use_asym_inv_sqrt for use_logsbt is questionable.', use_unit, "(2X,A)")
         end if
      end if
      ! relaxation and beyond LDA/GGA
      if ( ( use_geo_relaxation .or. compute_forces ) .and. use_prodbas) then
        if (use_full_spectrum) then
          if(myid.eq.0) then
            write(use_unit,*)
            write(use_unit, '(1X,A)') "* Error. Geometry relaxation does not work for correlated methods "
            write(use_unit, '(1X,A)') "* such as MP2, RPA, etc. We apologize for any possible inconvience."
          endif
          call aims_stop_coll ()

        elseif (use_lc_wpbeh) then
          if(myid.eq.0) then
            write(use_unit,*)
            write(use_unit, '(1X,A)') "* Error. Geometry relaxation does not work for range-separated "
            write(use_unit, '(1X,A)') "* functionals such as LC-wpbeh and M11, where the long-range "
            write(use_unit, '(1X,A)') "* contribution is 100% hartree-fock. We apologize for any possible inconvience."
          endif
          call aims_stop_coll ()

        elseif ( (RI_type .ne. RI_LVL) .and. (n_periodic.eq.0) ) then
          ! Stop with a warning in the cluster case.
          ! In periodic systems, we only have RI-LVL available anyway.
          if(myid.eq.0) then
            write(use_unit,*)
            write(use_unit, '(1X,A)') "* Error. Geometry relaxation for hybrid density functionals "
            write(use_unit, '(1X,A)') "* only works with the 'LVL_fast' version of RI."
            write(use_unit, '(1X,A)') "* Please set the keyword 'RI_method' to 'lvl_fast' or 'LVL_fast' . "
            write(use_unit, '(1X,A)') "* "
            write(use_unit, '(1X,A)') "* As you probably guessed, we could simply set 'RI_method lvl_fast' "
            write(use_unit, '(1X,A)') "* at this point even if you did not set the keyword explicitly in"
            write(use_unit, '(1X,A)') "* control.in and then run your calculation without disturbing you. "
            write(use_unit, '(1X,A)') "* Why do we not just do this? "
            write(use_unit, '(1X,A)') "* "
            write(use_unit, '(1X,A)') "* We will, in a future version. However, at present, we have two"
            write(use_unit, '(1X,A)') "* different variants of resolution of identity. RI-V is more"
            write(use_unit, '(1X,A)') "* accurate in general, and supports perturbative approaches,"
            write(use_unit, '(1X,A)') "* such as RPA, MP2, rPT2, GW, etc. RI-LVL does not support the"
            write(use_unit, '(1X,A)') "* latter methods yet, but is much faster for large systems"
            write(use_unit, '(1X,A)') "* with Hartree-Fock and hybrid functionals. In periodic systems,"
            write(use_unit, '(1X,A)') "* RI-LVL is even the only option."
            write(use_unit, '(1X,A)') "* "
            write(use_unit, '(1X,A)') "* So, we have two versions. At present, we do not wish to just"
            write(use_unit, '(1X,A)') "* silently choose for you. The risk is that someone might inadvertently"
            write(use_unit, '(1X,A)') "* calculate energy differences with different settings."
            write(use_unit, '(1X,A)') "* This is why we would like to make you aware of the two different"
            write(use_unit, '(1X,A)') "* approaches. However, sorry about the inconvenience of having."
            write(use_unit, '(1X,A)') "* stopped your run."
          endif
          call aims_stop_coll ()

        endif
      endif

      ! relaxation and molecular dynamics

      if (use_geo_relaxation .and. use_molecular_dynamics) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A)') "* Error. Both relaxation and molecular dynamics were requested within the same control.in."
            write(use_unit,'(1X,A)') "* At present, these two options are mutually exclusive - the behaviour of the code with"
            write(use_unit,'(1X,A)') "* this combination is not defined."
            write(use_unit,'(1X,A)') "* If you intend to first prerelax a structure, and then run molecular dynamics based on"
            write(use_unit,'(1X,A)') "* the resulting structure, please do this in two separate steps."
            write(use_unit,'(1X,A)') "* The present run will be aborted to alert you to the problem."
         end if
         call aims_stop_coll ()
      end if

      if (out_band .and. use_periodic_hf) then
         if(exx_band_structure_version .eq. 0)then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') ""
               write(use_unit,'(1X,A)') "* Error. Band structure and exact exchange are requested at the same time,"
               write(use_unit,'(1X,A)') "* but the keyword 'exx_band_structure_version' is not set."
               write(use_unit,'(1X,A)') "* We require that this keyword be set explicitly to make sure that our"
               write(use_unit,'(1X,A)') "* users understand a specific aspect of the current band structure"
               write(use_unit,'(1X,A)') "* calculation framework for hybrid density functionals with FHI-aims."
               write(use_unit,'(1X,A)') "* This has been used successfully by many groups for published work,"
               write(use_unit,'(1X,A)') "* but please be sure to understand the current usage. "
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "* There are two variants of the band structure calculation for hybrid"
               write(use_unit,'(1X,A)') "* density functionals: "
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "* --- exx_band_structure_version 1 :"
               write(use_unit,'(1X,A)') "*     This version is evaluated in a real-space formalism. It is fast"
               write(use_unit,'(1X,A)') "*     and the preferred version. HOWEVER, it does require a sufficiently"
               write(use_unit,'(1X,A)') "*     dense 'k_grid' to be used for the preceding s.c.f. cycle."
               write(use_unit,'(1X,A)') "*     If the 'k_grid used for s.c.f. is NOT dense enough, the real-"
               write(use_unit,'(1X,A)') "*     space approach will produce bands that are OBVIOUSLY incorrect."
               write(use_unit,'(1X,A)') "*     You cannot miss this - they won't look 'slightly' correct. In"
               write(use_unit,'(1X,A)') "*     that case, please try to increase the 'k_grid' parameters used"
               write(use_unit,'(1X,A)') "*     for the s.c.f. cycle. In particular, please avoid k_grid dimensions "
               write(use_unit,'(1X,A)') "*     of '1' (one) in any dimension. We apologize for this inconvenience."
               write(use_unit,'(1X,A)') "*     On the bright side, if you use exx_band_structure_version 1 "
               write(use_unit,'(1X,A)') "*     correctly, it will give reliable results without much overhead."
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "* --- exx_band_structure_version 2 : "
               write(use_unit,'(1X,A)') "*     This version is evaluated completely in reciprocal space and does"
               write(use_unit,'(1X,A)') "*     not have any problems with sparse s.c.f. 'k_grid' settings. "
               write(use_unit,'(1X,A)') "*     HOWEVER: exx_band_structure_version 2 is computationally VERY expensive."
               write(use_unit,'(1X,A)') "*     Please only use it for verification and/or if, for some reason, you"
               write(use_unit,'(1X,A)') "*     cannot make exx_band_structure_version 1 work at all. The time and memory"
               write(use_unit,'(1X,A)') "*     demand of exx_band_structure_version 2 is very significant."
               write(use_unit,'(1X,A)') "*     Ideally, at least work with a low number of k-points within band segments,"
               write(use_unit,'(1X,A)') "*     for example, 11 or less."
            endif
            call aims_stop_coll('', func)
         elseif(exx_band_structure_version .eq. 1)then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "* NOTE: Band structure and exact exchange are requested at the same time,"
               write(use_unit,'(1X,A)') "* and exx_band_structure_version is set to 1. "
               write(use_unit,'(1X,A)') "* exx_band_structure_version 1 is evaluated in a real-space formalism. It is fast"
               write(use_unit,'(1X,A)') "* and the preferred version. HOWEVER, it does require a sufficiently"
               write(use_unit,'(1X,A)') "* dense 'k_grid' to be used for the preceding s.c.f. cycle."
               write(use_unit,'(1X,A)') "* If the 'k_grid used for s.c.f. is NOT dense enough, the real-"
               write(use_unit,'(1X,A)') "* space approach will produce bands that are OBVIOUSLY incorrect."
               write(use_unit,'(1X,A)') "* You cannot miss this - they won't look 'slightly' correct. In"
               write(use_unit,'(1X,A)') "* that case, please try to increase the 'k_grid' parameters used"
               write(use_unit,'(1X,A)') "* for the s.c.f. cycle. In particular, please avoid k_grid dimensions "
               write(use_unit,'(1X,A)') "* of '1' (one) in any dimension. We apologize for this inconvenience."
               write(use_unit,'(1X,A)') "* On the bright side, if you use exx_band_structure_version 1 "
               write(use_unit,'(1X,A)') "* correctly, it will give reliable results without much overhead."
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "* If your band structure does not look reasonable, increase the s.c.f. k-point"
               write(use_unit,'(1X,A)') "* grid density. Again, it is especially recommended that you avoid dimensions of 1 when"
               write(use_unit,'(1X,A)') "* specifying k_grid and using this option, even for large unit cells where a"
               write(use_unit,'(1X,A)') "* gamma-point-only k-grid would be sufficient to converge the accuracy of the"
               write(use_unit,'(1X,A)') "* calculation."
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "* Please test (e.g., by comparing results from different k_grid settings) when in any doubt."
            endif
         elseif(exx_band_structure_version .eq. 2)then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "*** WARNING: Band structure and exact exchange are requested at the same time,"
               write(use_unit,'(1X,A)') "*** and exx_band_structure_version is set to 2. This implies that you request"
               write(use_unit,'(1X,A)') "*** a band structure calculation based on explicit calculation of the Fock matrix"
               write(use_unit,'(1X,A)') "*** for each k-point."
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "*   !! In general, exx_band_structure_version 1 should be preferred (but read the"
               write(use_unit,'(1X,A)') "*   comments associated with its use!)."
               write(use_unit,'(1X,A)') "* "
               write(use_unit,'(1X,A)') "*** exx_band_structure_version 2 is intended as a fallback only and is computationally expensive."
               write(use_unit,'(1X,A)') "*** The number of k-points within band segments should not exceed around 11."
               write(use_unit,'(1X,A)') "*** The program will continue."
            endif
            if(use_scalapack)then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') ""
                  write(use_unit,'(1X,A)') "* Error. Scalapack version of exx_band_structure_version = 2 is not yet implemented."
                  write(use_unit,'(1X,A)') "* Use exx_band_structure_version = 1. The program will stop."
               endif
               call aims_stop_coll('', func)
            endif
            if(use_local_index.or.use_load_balancing)then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)') ""
                  write(use_unit,'(1X,A)') "* Error. exx_band_structure_version = 2 does not support use_local_index or load "
                  write(use_unit,'(1X,A)') "* balancing, since it was not written with ScaLAPACK support, and thus neither of "
                  write(use_unit,'(1X,A)') "* these methods are possible."
                  write(use_unit,'(1X,A)') "* That being said, if you are seeing this keyword, it is possible that someone has"
                  write(use_unit,'(1X,A)') "* added ScaLAPACK support.  Support for use_local_index and load balancing is then"
                  write(use_unit,'(1X,A)') "* possible, but someone needs to implement it."
                  write(use_unit,'(1X,A)') "* Use exx_band_structure_version = 1. The program will stop."
               endif
               call aims_stop_coll('', func)
            endif

         else
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') "* Error. Unknown exx_band_structure_version. The program will stop."
            endif
            call aims_stop_coll('', func)
         endif
      endif

      if (out_band .or. out_band_mulliken) then ! add mlk
         real_eigenvectors = .false.   ! XR: this appears to me a proper setting, or did I miss anything?
      endif

      if (postscf_eigenvalues) then
         if (n_periodic.eq.0) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* Warning. Found "output postscf_eigenvalues" request for a non-periodic geometry.'
               write(use_unit,'(1X,A)') &
                  '* In this case and unless the output_level is set to "light", the eigenvalues'
               write(use_unit,'(1X,A)') &
                  '* are always written to the standard output file. No separate file will be attempted.'
            end if
            postscf_eigenvalues = .false.
         else if (.not.pert_dos_on) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* Warning. Found "output postscf_eigenvalues" request but no dos_kgrid_factors.'
               write(use_unit,'(1X,A)') &
                  '* Switching "output_dos_kgrid_factors on with minimal values (1). '
               write(use_unit,'(1X,A)') &
                  '* However, additional restrictions currently apply along this code path and may'
               write(use_unit,'(1X,A)') &
                  '* lead to errors below. In that case, please work with us to identify them.'
            end if
            pert_dos_on = .true.
            n_k_points_dos_xyz(1) = 1
            n_k_points_dos_xyz(2) = 1
            n_k_points_dos_xyz(3) = 1
         end if
      end if

      if (pert_dos_on) then
         ! There are, unfortunately, a number of restrictions for which the post-processed density of
         ! states calculation has not yet been enabled or proven.
         ! We catch some of these here. In all these cases, help enabling the feature would be
         ! appreciated.
         if ( (.not.out_dos) .and. (.not.postscf_eigenvalues)) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* Warning. Found "dos_kgrid_factors" request but neither a DOS nor'
               write(use_unit,'(1X,A)') &
                  '* postscf_eigenvalues are requested. dos_kgrid_factors will be ignored.'
            end if
            pert_dos_on = .false.
         end if
         if (n_periodic.eq.0) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* Warning. Found "dos_kgrid_factors" request for a non-periodic geometry.'
               write(use_unit,'(1X,A)') &
                  '* In this case and unless the output_level is set to "light", the eigenvalues'
               write(use_unit,'(1X,A)') &
                  '* are always written to the standard output file. No separate output will be written.'
            end if
            pert_dos_on = .false.
         end if
         if (use_load_balancing) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)')
               write(use_unit,'(1X,A)') &
                  '* Error. Neither "use_local_index" nor "load_balancing" support perturbed DOS output.'
               write(use_unit,'(1X,A)') &
                  '* Either omit the dos_kgrid_factors keyword to output a non-perturbed DOS or omit the'
               write(use_unit,'(1X,A)') &
                  '* use_local_index and load_balancing keywords.'
               write(use_unit,'(1X,A)') &
                  '* It would not be rocket science to change this in the code.'
               write(use_unit,'(1X,A)') &
                  '* Please support us if you really need this combination.'
            end if
            call aims_stop_coll('Keyword mismatch', func)
            pert_dos_on = .false.
         end if
         if (use_local_index) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)')
               write(use_unit,'(1X,A)') &
                  '* Error. "use_local_index" does not support perturbed DOS output.'
               write(use_unit,'(1X,A)') &
                  '* Either omit the dos_kgrid_factors keyword to output a non-perturbed'
               write(use_unit,'(1X,A)') &
                  '* DOS or omit the use_local_index keyword.'
               write(use_unit,'(1X,A)') &
                  '* It would not be rocket science to change this in the code.'
               write(use_unit,'(1X,A)') &
                  '* Please support us if you really need this combination.'
            end if
            call aims_stop_coll('Keyword mismatch', func)
            pert_dos_on = .false.
         end if
         if (use_symmetry_reduced_spg) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* Error. Reduction of k-points with spglib does not support DOS output.'
               write(use_unit,'(1X,A)') &
                  '* Only output based on irreducible set of k-point.'
               write(use_unit,'(1X,A)') &
                  '* It would not be rocket science to change this in the code.'
               write(use_unit,'(1X,A)') &
                  '* Please support us if you really need this combination.'
            end if
            call aims_stop_coll('Keyword mismatch', func)
         end if
      endif

      ! MD_QH_init
      if (use_MD_QH_init) then

         ! Check for MD_restart:
         if (MD_initialconf_from_restart) then
            if (myid.eq.0) then
               write(use_unit,'(1X,A)')      "*** ERROR:"
               write(use_unit,'(1X,A)')      "*** Keywords MD_restart and MD_QH_init are not compatible."
            end if
            call aims_stop_coll ()
         end if

         do i_MD_QH_init_segments=1,MD_QH_init_segments

            ! Check for atom numbering
            if ((MD_QH_first_atom(i_MD_QH_init_segments) .lt. 1) .OR. &
               (MD_QH_first_atom(i_MD_QH_init_segments) .gt. n_atoms) .OR. &
               (MD_QH_last_atom(i_MD_QH_init_segments)  .lt. 1) .OR. &
               (MD_QH_last_atom(i_MD_QH_init_segments) .gt. n_atoms)) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)')      "*** ERROR:"
                  write(use_unit,'(1X,A)')      "*** Files control.in and geometry.in contain inconsistent specifications:"
                  write(use_unit,'(1X,A,I5,A)') "*** MD_QH_init statement",i_MD_QH_init_segments," refers to non-existing atoms."
               end if
               call aims_stop_coll ()
            end if

            ! Check for atom sequence
            if (MD_QH_first_atom(i_MD_QH_init_segments) .gt. MD_QH_last_atom(i_MD_QH_init_segments)) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)')      "*** ERROR:"
                  write(use_unit,'(1X,A)')      "*** Index of first atom is larger than index of last atom in"
                  write(use_unit,'(1X,A,I5,A)') "*** MD_QH_init statement",i_MD_QH_init_segments," ."
               end if
               call aims_stop_coll ()
            end if

            ! Check for negative temperatures
            if (MD_QH_temperature(i_MD_QH_init_segments) .lt. 0.0d0 ) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)')      "*** ERROR:"
                  write(use_unit,'(1X,A,I5,A)') "*** Temperature specified in MD_QH_init statement",i_MD_QH_init_segments," is"
                  write(use_unit,'(1X,A)')      "*** negative. Please specify the temperature in Kelvin."
               end if
               call aims_stop_coll ()
            end if

            ! Check for negative temperatures
            flag_MD_QH_file_exists = .false.
            inquire(FILE=MD_QH_file(i_MD_QH_init_segments),EXIST=flag_MD_QH_file_exists)
            if (.not. flag_MD_QH_file_exists ) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)')      "*** ERROR:"
                  write(use_unit,'(1X,A,A,A)')  "*** File ",trim(MD_QH_file(i_MD_QH_init_segments))," specified"
                  write(use_unit,'(1X,A,I5,A)') "*** in MD_QH_init statement",i_MD_QH_init_segments," does not exist."
               end if
               call aims_stop_coll ()
            end if

            ! Check for too high temperatures (WARNING ONLY)
            if (MD_QH_temperature(i_MD_QH_init_segments) .gt. 1000.0d0 ) then
               if (myid.eq.0) then
                  write(use_unit,'(1X,A)')      "*** WARNING:"
                  write(use_unit,'(1X,A,I5,A)') "*** Temperature specified in MD_QH_init statement",i_MD_QH_init_segments," is"
                  write(use_unit,'(1X,A)')      "*** higher than 1000K. Most probably, the quasi-harmonic approximation breaks "
                  write(use_unit,'(1X,A)')      "*** down at such high temperatures."
               end if
            end if

         end do
      end if

      if (out_cube .and. .not. flag_wave_threshold .and. .not.restart_read) then
          write(info_str,*) ' * Volumetric cube output requested. Tightening default wave_threshold to 1d-8'
          call localorb_info(info_str)
          wave_threshold=1d-8
      elseif (out_cube .and. wave_threshold.gt.1d-7) then
         if (myid.eq.0) then
          write(info_str,*) '* WARNING: Volumetric cube output requested, but wave_threshold is loose. *'
          call localorb_info(info_str)
          write(info_str,*) '* In some cases, this can lead to kinked cube files. For production calc- *'
          call localorb_info(info_str)
          write(info_str,*) '* ulations, it is recommended to run the calculation with wave_threshold  *'
          call localorb_info(info_str)
          write(use_unit,*) '* 1d-7 or tighter. Please note that changing wave_threshold breaks RESTART-files  *'
          call localorb_info(info_str)
          write(use_unit,*) '* so these can not be used. We apologize for any inconvenience.                  *'
          call localorb_info(info_str)
          endif
      endif

      if (n_periodic.eq.0.and.abs(charge).gt.1e-8.and.distribute_leftover_charge) then
             write(use_unit,*) 'distribute_leftover_charge incompatible with'
             write(use_unit,*) 'charged systems. Please correct'
             call aims_stop
      endif

      if (restart_from_density.and.use_hartree_fock) then
           call aims_stop('density restart cannot be used with hybrid functionals')
      endif

      if (mixer.eq.MIX_BROYDEN.and.flag_ini_linear_mixing) then
          call aims_stop('Broyden mixer does not support initial linear mixing')
      endif

      if (use_mbd_old.and.use_vdw_correction_hirshfeld) then
          write(info_str,'(1X,A)') "*** Keyword conflict between 'many_body_dispersion' and 'vdw_correction_hirshfeld'."
          call localorb_info(info_str)
          call aims_stop("Please correct control.in")
      endif

      if ( flag_neutral_excitation .and. neutral_excitation_tddft ) then
      if ( flag_fxc==-1 .and. .not.use_libxc_tddft ) then
          write(info_str,'(2X,2A)') "*** Calculation of optical excitation energies with TDDFT ", &
                                    "('neutral_excitation tddft') explicitly requires to set the TDDFT kernel."
          call localorb_info(info_str)
          call aims_stop("Please set keyword 'tddft_kernel' in control.in")
      endif
      endif

      if (.not.flag_vdw_convergence_threshold) then
         ! This earlier choice was a trap, since it could lead
         ! to unreasonably long vdW summation times by a non-obvious
         ! connection to sc_accuracy_etot. Left comment here, but a new
         ! default, which also accounts for the number of atoms in the
         ! structure, is proposed below.
         !if (flag_acc_etot) then
         !   ! use old default
         !   if (sc_accuracy_etot.gt.0.d0) then
         !      vdw_convergence_threshold = sc_accuracy_etot/hartree
         !   end if
         !end if

         ! The default vdw_convergence_threshold in runtime_choices.f90
         ! is now taken to be per atom. In addition, a lower floor of 10^-6 eV is set
         ! to avoid unreasonably lengthy lattice sums for small unit cells.
         if (n_atoms.gt.0) then
           vdw_convergence_threshold = vdw_convergence_threshold * dble(n_atoms)
           if (vdw_convergence_threshold*hartree .lt. 1.d-6) then
              vdw_convergence_threshold = 1.d-6/hartree
           end if
         end if
         if (use_vdw_correction_hirshfeld.and.(n_periodic.gt.0)) then
            write(info_str,'(2X,A)') "Using the following vdw_convergence_threshold for supercell sums"
            call localorb_info(info_str)
            write(info_str,'(2X,A,E13.6,A)') "over the Tkatchenko-Scheffler correction: ", &
                 vdw_convergence_threshold*hartree, ' eV.'
             call localorb_info(info_str)
         end if
      end if

      ! CC: Force eigenvectors to be complex makes no sense in cluster case
      if ( ( flag_force_complex_eigenvectors ) .and. (n_periodic .eq. 0 ) ) then
        write(info_str,'(2X,A)') "*** Forcing eigenvectors to be complex is ONLY supported under periodic boundary conditions"
        call localorb_info(info_str,use_unit,'(A)')
        call aims_stop("*** Please correct your control.in file")
      end if

     !Check for the self-consistent TS-vdW implementation
      if (use_vdw_correction_hirshfeld_sc) then
     !The vdW forces are not implemented self-consistently, the code stops
         if (use_forces) then
            write(info_str,'(2X,4A)') '*** WARNING: Forces are not implemented for self-consistent vdW effects. Use the post-process vdW correction to relax the structure.'
            call localorb_info(info_str)
            write(info_str,'(2X,A)') '*** Please use the correct flag in control.in. Aborting execution. '
            call localorb_info(info_str)
            call aims_stop()
         end if
      endif

      if (flag_found_unfold) then
         write(info_str,'()') ! end line
         call localorb_info(info_str, use_unit,'(A)')
         write (info_str,'(2X,A)') "Band structure unfolding:"
         call localorb_info(info_str,use_unit,'(A)')
         if (n_periodic.eq.0) then
            write (info_str,'(2X,A)') "*** Warning. The band structure unfolding request will have no effect in a non-periodic calculation."
            call localorb_info(info_str,use_unit,'(A)')
         else if (n_plot_band.le.0) then
            write (info_str,'(1X,A)') "*** Warning. The band structure unfolding request will have no effect since you did not specify a band structure output."
            call localorb_info(info_str,use_unit,'(A)')
         else if (use_scalapack) then
            flag_hdf5_unfold  = .true.
            flag_ASCII_unfold = .false.
            write(info_str,'(2X 2A)')"| Band structure unfolding + ScaLapack: ", &
                 "The output format is set to HDF5"
            call localorb_info(info_str,use_unit,'(A)')
         else
            write(info_str,'(2X A)') "| The overlap matrix will be written out for the k-points requested by 'output band'."
            call localorb_info(info_str,use_unit,'(A)')
            write (info_str,'(2X,A)')"| Kohn-Sham eigenvectors will be written for the k-points requested by 'output band'."
            call localorb_info(info_str,use_unit,'(A)')
         end if
         write(info_str,'()') ! end line
         call localorb_info(info_str, use_unit,'(A)')
      end if

      ! Pre-calc options that are important for FO-DFT! These have to be set before read_control is
      ! finished!
      if (fo_dft.and.fo_dft_option.eq.'final') then
          call fodft_set_states(n_empty)
          write(info_str,'(2X,A, I4)') "Number of empty states calculated from fragments: ", n_empty
          call localorb_info(info_str)
          call fodft_combine()
      endif
      ! Checks for FO-DFT
      if (fo_dft.and.(n_periodic.ne.0)) then
          call aims_stop('FO-DFT not implemented for periodic calculations!')
      endif

      if (fo_dft.and.occ_fo) then
          write(info_str,*) '* WARNING: FO-DFT based on occupation reset is experimental and NOT tested *'
          call localorb_info(info_str)
          write(info_str,*) '* very well. Use with great care! (or not at all)*'
          call localorb_info(info_str)
      endif

      if (fo_finalise) then

      if (use_fo_potential.and.(n_empty_atoms.eq.0)) then
          call aims_stop("No empty sites for other fragment found in geometry.in. This is a requirement for polarized FODFT. See the manual!")
      end if

      if (fo_finalise.and.use_load_balancing) then
          call aims_stop("Load balancing (load_balancing) and/or local storage of matrices (use_local_index) is not supported for FODFT.")
      endif
      if (fo_finalise.and.use_local_index) then
          call aims_stop("Local storage of matrices (use_local_index) is not supported for FODFT.")
      endif
      endif !fo_finalise

      if (force_single_restartfile) then
              if (n_kpts>1) then
                 call aims_stop("force_single_restartfile option not possible for periodic calculations with more than one k-point!")
              endif
              if (.not.collect_eigenvectors) then
                 call aims_stop("force_single_restartfile option not possible with 'collect_eigenvectors .false.'!")
              endif
      endif

      if (restart_eigenvectors_periodic) then
          if (.not.collect_eigenvectors) then
              call aims_stop("restart_eigenvectors_periodic option is not possible with 'collect_eigenvectors .false.'!")
          endif
          if (force_single_restartfile) then
              call aims_stop("'restart_eigenvectors_periodic' and 'force_single_restartfile' can not be used at the same time!")
          endif
      endif

      ! Consistency checking for spin-orbit coupling
      if (calculate_perturbative_soc) then
         !TZ add the inconsistency here
        if (out_cube) then
         do cube_id = 1, i_cube, 1
            if (cube_type_needs_densmat(cube_id)) then
                  call localorb_info("")
                  write(info_str, "(2X,A,' ',A,' ',A)") &
                     & 'Cube type', trim(cube_type(cube_id)), &
                     & 'not implemented for load_balancing or use_local_index keywords.  Exiting.'
                  call aims_stop_coll(info_str, func)
            endif
         enddo
        endif


        ! Check for conditions common to all handlings of SOC
        if(out_cube) then
          write(info_str,'(2X,A)') "*** At present, periodic spin-orbit coupling does not output cube"
          call localorb_info(info_str)
          write(info_str,'(2X,A)') "*** files.  This is unlikely to change soon.  Exiting."
          call localorb_info(info_str)
!          call aims_stop_coll('',func)
        end if

        out_cube_soc = out_cube .and. calculate_perturbative_soc
        if (flag_run_mulliken .or. flag_out_dielectric .or. flag_out_absorption .or. &
             flag_out_dipmat .or. &
             out_cube_soc .or. n_write_soc_eigenvectors .gt. 0) then
          save_soc_perturbed_eigenvectors = .true.
        end if
      end if

      if (n_write_soc_eigenvectors.gt.0 .and. .not.calculate_perturbative_soc) then
        write(info_str,'(2X,A)') "*** To write out spin-orbit-coupled eigenvectors, one must set the include_spin_orbit keyword"
        call localorb_info(info_str)
        write(info_str,'(2X,A)') "*** in control.in."
        call localorb_info(info_str)
        write(info_str,'(2X,A)') "*** Exiting."
        call localorb_info(info_str)
        call aims_stop_coll('', func)
      end if

      if(gamma_cut_coulomb .and. use_gw_gamma_corr .and. use_gw) then
        if(myid.eq.0) then
               write(use_unit,'(1X,A)') &
                  '* The explicit gamma point correction treatment and the cut Coulomb operator for '
               write(use_unit,'(1X,A)') &
                  '* cannot be both true for periodic GW calculations; setting "gamma_cut_coulomb" to be false. '
        endif
        gamma_cut_coulomb = .false.
      endif

      ! MPE implicit solvation
      select case (solvent_method)
      case (SOLVENT_MPE)
         if (n_periodic > 0) then
            write(info_str,'(2X,A)') &
               "The current implementation of MPE does not support "//&
               "periodic calculations."
            call localorb_info(info_str)
            call aims_stop_coll(info_str, func)
         endif
         if (  flag_xc .eq. 0 .or. &
               is_xc_lda(normalize_flag_xc(flag_xc)) ) then
            write(info_str,'(2X,A)') &
               "The current implementation of MPE needs an analytical density "//&
               "gradient and thus cannot be used with the specified "//&
               "XC functional."
            call localorb_info(info_str)
            call aims_stop_coll(info_str, func)
         endif
         if ( (isc_cavity_type.eq.ISC_CONST % CAVITY_UNDEF) &
               .or.(isc_cavity_type.eq.ISC_CONST % CAVITY_INVAL) ) then
            write(info_str,'(2X,A)') &
               "isc_cavity_type needs to be specified for MPE calculations."
            call localorb_info(info_str)
            call aims_stop_coll(info_str, func)
         endif
         if ( (.not. isc_calculate_surface_and_volume) .and. &
               (mpe_nonel_model .eq. MPE_CONST % NONEL_linOV) ) then
            write(info_str,'(2X,A)') &
               "The chosen non-electrostatic model needs the surface area "//&
               "and/or the volume of the cavity. Switching on the "//&
               "corresponding flag."
            call localorb_info(info_str)
            isc_calculate_surface_and_volume = .true.
         endif
      case (SOLVENT_INVALID)
         write(info_str,'(2X,A)') &
            "Illegal solvent method has been requested."
         call localorb_info(info_str)
         call aims_stop_coll(info_str, func)
      end select

      if (any(.not. flag_cut_free_atom .and. .not. (use_prodbas .and. &
           & .not. (use_lvl .and. use_logsbt)) .and. n_gaussian > 0)) &
           call localorb_multi('Distances at which the Gaussian orbitals are &
           &below wave_threshold are the following:', format='(2x, a)')
      do i_species = 1, n_species
         if (.not. flag_cut_free_atom(i_species) .and. &
              & .not. (use_prodbas .and. .not. (use_lvl .and. use_logsbt)) &
              & .and. n_gaussian(i_species) > 0) &
              & call set_gaussian_free_cut(i_species, wave_threshold, &
              & free_r_cut(i_species))
      end do

      if ( use_symmetry_reduced_spg .and. use_analytical_stress ) then
        call localorb_info("")
        write(info_str, "(2X,A,A)") &
           & '*** Symmetry-reduced k-grid and analytical stress not yet', &
           & ' supported. Exiting.'
        call aims_stop_coll(info_str, func)
      end if

      if ( use_rlsy_symmetry ) then
        if ( n_periodic .ne. 3 ) then
          call localorb_info("")
          write(info_str, "(2X,A,A)") &
             & '*** Symmetry-reduction only available for periodic systems.', &
             & ' Exiting.'
          call aims_stop_coll(info_str, func)
        endif
      endif


      call initialize_psi_at_nucleus()

      if ( solvent_method .eq. SOLVENT_MPE .and. &
              use_forces ) then
         write(info_str, '(A)') '*** Forces not implemented for MPE implicit solvation'
         call aims_stop_coll(info_str, func)
      endif

      if (out_mommat) then
         ! check restrictions to the implementation of output_momentummatrix and make sure
         ! they are fulfilled

         if (use_scalapack) then
            if (.not.flag_collect_eigenvectors) then
            ! In this case, the user did not set a particular value for collect_eigenvectors.
            ! We need to set it to the value needed for compute_momentummatrix .
               if (.not.collect_eigenvectors) then
                  write(info_str,'(2X,A)') &
                    "compute_momentummatrix requested: Setting 'collect_eigenvectors' to be '.true.'."
                  call localorb_info(info_str)
                  write(info_str,'(2X,A)') &
                    "This will increase memory usage in large parallel runs."
                  call localorb_info(info_str)
                  collect_eigenvectors = .true.
               endif
            else
               ! In this case, the user set a particular value for collect_eigenvectors. We
               ! need to make sure that the user's request does not clash with
               ! compute_momentummatrix . If it does, we provide an informative warning and stop.
               if (.not.collect_eigenvectors) then
                  write(info_str,'(2X,A)') &
                    "** compute_momentummatrix requested but 'collect_eigenvectors' was set to be '.false.' in control.in."
                  call localorb_info(info_str)
                  write(info_str,'(2X,A)') &
                    "** collect_eigenvectors .false. is not supported by the compute_momentummatrix keyword."
                  call localorb_info(info_str)
                  write(info_str,'(2X,A)') &
                    "** Please correct your input file or contact the developers."
                  call localorb_info(info_str)
                  call aims_stop_coll('collect_eigenvectors .false. not supported by compute_momentummatrix.', func)
               end if
            end if
         endif ! end if (use_scalapack)

         if (.not.flag_symmetry_reduced_k_grid) then
            if (use_symmetry_reduced_k_grid) then
                 write(info_str,'(2X,A)') &
                 & "You have requested the compute_momentummatrix keyword."
                 call localorb_info(info_str)
                 write(info_str,'(2X,A)') &
                 & "symmetry_reduced_k_grid will be turned off."
                 call localorb_info(info_str)
            end if
            use_symmetry_reduced_k_grid = .false.
         else
            if (use_symmetry_reduced_k_grid) then
                 write(info_str,'(2X,A)') &
                 & "** compute_momentummatrix keyword and symmetry_reduced_k_grid both requested in control.in."
                 call localorb_info(info_str)
                 write(info_str,'(2X,A)') &
                 & "** symmetry_reduced_k_grid is not supported for compute_momentummatrix. Please correct."
                 call localorb_info(info_str)
                 call aims_stop_coll('symmetry_reduced_k_grid not supported by compute_momentummatrix.', func)

            end if
         end if ! end checking flag_symmetry_reduced_k_grid

      end if  ! end checking of the 'out_mommat flag for consistency

   end subroutine checkinput

! **************************************************************************************************
!> brief reads in relativistic corrections for QP energies
! o n_qps -- number of correct QP energies
! o state -- corrected states
! o values -- relativistic corrections for the respective states
! **************************************************************************************************
   subroutine read_relativistic_qp_corrections(inputline)
  
      character(len=250), intent(in)             :: inputline

      logical                                    :: explicit
      integer                                    :: i_counter
      real(kind=8)                               :: temp_float
 
      i_counter=-1 ! so value is correct on exit
      explicit=.true.
      do while (explicit)
         i_counter=i_counter+1 
         
         call get_optional_parameter(inputline, param_position=3 + 2 * i_counter, &
              real_parameter=temp_float, explicit=explicit)
      end do
      gw_cd%qp_adjust%n_qps = i_counter

      write(info_str,'(2X,A)') &
           'Relativistic correction currently indicated only for 1s core states of 2nd period elements.'
      call localorb_info(info_str,use_unit,'(A)')
      call allocate_qp_adjust(gw_cd)
      do i_counter = 1, gw_cd%qp_adjust%n_qps
         call get_optional_parameter(inputline, param_position=2 * i_counter, &
              int_parameter=gw_cd%qp_adjust%state(i_counter),&
              explicit=explicit)
         
         call get_optional_parameter(inputline, param_position=2 * i_counter + 1, &
              real_parameter=gw_cd%qp_adjust%values(i_counter),&
              explicit=explicit)

         write(info_str,'(2X,A,I3,A,F10.6,A)') &
           'Adjusting state ',gw_cd%qp_adjust%state(i_counter),' by relativistic correction ', &
            gw_cd%qp_adjust%values(i_counter),' eV'
         call localorb_info(info_str,use_unit,'(A)')
      end do
   end subroutine read_relativistic_qp_corrections

! **************************************************************************************************
!> brief reads additional input of the type "keyword parameter1 parameter2", where parameter2 is
!        optional
! o inputline --  whole inputline
! o param_position -- position of the optional parameter; i.e. position of parameter2 would be "3"
! o real_parameter -- optional real parameter
! o int_parameter -- optional integer parameter
! o explicit -- if parameter is given
! **************************************************************************************************
   subroutine get_optional_parameter(inputline,param_position,real_parameter,int_parameter,explicit)

      implicit none

      character(len=250), intent(in)             :: inputline
      integer, intent(in)                        :: param_position
      real(kind=8), intent(out),optional         :: real_parameter
      integer, intent(out),optional              :: int_parameter
      logical, intent(out)                       :: explicit

      integer                                    :: lastpos, icount, i
      character(len=1)                           :: ch
      real(kind=8)                               :: temp_parameter

      lastpos=-1
      icount=0
      temp_parameter=0.d0

      do i=1,len(inputline)
         ch=inputline(i:i)
         if((ch==' '.or.i==len(inputline)).and.lastpos/= -1) then
           icount = icount + 1
           if(icount == param_position) then
             read(inputline(lastpos:i-1),*) temp_parameter
             exit
           endif
           lastpos = -1
         endif
         if(ch/=' '.and.lastpos==-1) then
           lastpos=i
         endif
      enddo

      if(present(real_parameter)) then
        real_parameter = temp_parameter
      endif
      if(present(int_parameter)) then
        int_parameter = NINT(temp_parameter)
      endif

      explicit = .false.
      if(icount == param_position) then
        explicit = .true.
      endif

   end subroutine get_optional_parameter

end subroutine read_control
!******
