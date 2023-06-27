# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

# CREATE TARGET
# First we create a target with no source files. The source will
# follow later.
if (BUILD_LIBRARY)
  add_library(aims main.f90)
else()
  add_executable(aims main.f90)
endif()
set_target_properties(aims
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

# The linker language should be Fortran. Otherwise, if USE_CXX_FILES
# is ON, it tries to create a C++ executable which can lead to some
# very obscure errors during linking.
set_target_properties(aims PROPERTIES LINKER_LANGUAGE Fortran)

target_include_directories(aims PRIVATE ${INC_PATHS})
set_target_properties(aims PROPERTIES OUTPUT_NAME ${TARGET_NAME})

# Location of the module files is 'modules' in the build directory
set_target_properties(aims PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
target_include_directories(aims PRIVATE ${PROJECT_BINARY_DIR}/modules)

# INSTALL
# This is pretty basic for now and can be expanded when the structure
# of aims becomes more complex.
if (NOT BUILD_LIBRARY)
  install(TARGETS aims
    EXPORT aimsConfig
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  install(EXPORT aimsConfig NAMESPACE aims:: DESTINATION lib/cmake/aims)
  include(CMakePackageConfigHelpers)
  write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/aimsConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/aimsConfigVersion.cmake
    DESTINATION lib/cmake/aims)
  target_include_directories(aims INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
  export(EXPORT aimsConfig NAMESPACE aims::
    FILE ${PROJECT_BINARY_DIR}/aimsConfig.cmake)
endif()
# When installed, aims can be included in a different project by
# calling find_package(aims x.x.x), where x.x.x (optional) is the
# minimum required version. If found, linking is as simple as
# target_link_libraries(Foo aims). If not found, see the CMake manual
# for an additional PATHS argument to find_package.

configure_file(cmake_info.f90.in cmake_info.f90)
target_include_directories(aims PRIVATE ${PROJECT_BINARY_DIR})

# SOURCE
# Finally, include the source files. The PRIVATE keyword of
# target_sources means that the sources are added to the aims target
# only and are not directly included in other targets when linking
# against the aims target. Simply put, always use PRIVATE unless there
# is a good reason to use PUBLIC or INTERFACE (see the CMake manual
# for more information).
if (NOT BUILD_LIBRARY)
  # Whether to split the calculation to multiple aims instances
  if (MULTIAIMS)
    _target_sources(aims PRIVATE multiaims.f90)
  else()
    _target_sources(aims PRIVATE aims.f90)
  endif()
endif()
# Create the auxiliary targets
add_library_f77(aims6)
add_library_normal(aims1 aims6)
add_library_f77(aims8 aims1)
add_library_normal(aims7 aims8)
add_library_reduced(aims2 aims7)
add_library_normal(aims3 aims2)
add_library_minimal(aims4 aims3)
add_library_normal(aims5 aims4)
target_link_libraries(aims PRIVATE aims5)
static_linking_fix(aims aims5 LIBS)
target_include_directories(aims5 PRIVATE ${INC_PATHS})

# MPI stubs
if (NOT USE_MPI)
  target_include_directories(aims1 PUBLIC external) # For mpif.h
  target_sources(aims5 PRIVATE mpi_stub.f90)
endif()

# Scalapack stubs
if (NOT USE_SCALAPACK)
  target_sources(aims5 PRIVATE scalapack_stub.f90)
endif()

# Architecture
if (ARCHITECTURE STREQUAL "PGI")
  target_sources(aims1 PRIVATE arch_pgi.f90)
else()
  target_sources(aims1 PRIVATE arch_generic.f90)
endif()

# C files
if (USE_C_FILES)
  target_sources(aims5 PRIVATE
    get_stacksize.c
    check_stacksize.f90
    change_directory.f90
    change_directory_c.c
    unlimit_stack.c)
else()
  target_sources(aims5 PRIVATE
    check_stacksize_stub.f90
    change_directory_stub.f90)
endif()

# ELSI
if (EXTERNAL_ELSI_PATH)
  find_package(elsi 2.4 REQUIRED NO_DEFAULT_PATH PATHS ${EXTERNAL_ELSI_PATH})
  target_link_libraries(aims PRIVATE elsi::elsi)
  static_linking_fix(aims elsi::elsi LIBS)
  mark_as_advanced(elsi_DIR)
  set(USE_EXT_ELPA ON)
else()
  # External ELPA using internal ELSI
  if (USE_EXTERNAL_ELPA)
    set(USE_EXT_ELPA ON)
  endif()
endif()

# SPGlib (see also external/spglib/src)
if (USE_SPGLIB)
  target_sources(aims1 PRIVATE
    spglib_symm.f90
    spglib_f2c.f90)
else()
  target_sources(aims1 PRIVATE spglib_symm_stub.f90)
endif()

# Python 2 or 3 interface to FHI-aims
if (USE_CFFI)
  include(cffi/CMakeLists.txt)
else()
  target_sources(aims3 PRIVATE python_interface_stub.f90)
endif()

# i-Pi Python interface
if (NOT USE_iPI)
  target_sources(aims5 PRIVATE ipi_stubs.f90)
endif()

# HDF5
if (NOT USE_HDF5)
  target_sources(aims5 PRIVATE hdf5_stubs.f90)
endif()

if (USE_LIBXC)
  target_sources(aims5 PRIVATE libxc_yes.f90)
else()
  target_sources(aims5 PRIVATE libxc_no.f90)
endif()
if (USE_RLSY)
   add_subdirectory(rlsy)
else()
   add_subdirectory(rlsy_stubs)
endif()

# Subdirectories
add_subdirectory(basis_sets)
add_subdirectory(ipc)
add_subdirectory(elpa)
add_subdirectory(optical_response)
add_subdirectory(MagneticResponse)
add_subdirectory(soc)
add_subdirectory(SCGW)
add_subdirectory(Periodic_PT2)
add_subdirectory(LRC_PT2)
add_subdirectory(mbd-std)
add_subdirectory(mbd-dev)
add_subdirectory(DMFT_embed)
add_subdirectory(DFPT)
add_subdirectory(DFPT_reduce_memory)
add_subdirectory(DFPT_polarizability)
add_subdirectory(DFPT_polar_reduce_memory)
add_subdirectory(DFPT_phonon)
add_subdirectory(DFPT_phonon_reduce_memory)
add_subdirectory(DFPT_dielectric)
add_subdirectory(DFPT_phonon_gamma)
add_subdirectory(sBGE2)
add_subdirectory(RRS-PBC)
add_subdirectory(xc_dfauto)
add_subdirectory(FCIQMC)
add_subdirectory(bse)
add_subdirectory(CCcalc)
add_subdirectory(qpe_osPT2)
add_subdirectory(friction)
add_subdirectory(external)
add_subdirectory(relativity)
if (NOT EXTERNAL_ELSI_PATH)
  target_link_libraries(aims PRIVATE elsi)
  static_linking_fix(aims elsi LIBS)
endif()
# External projects to link against
if (USE_LIBXC AND LIBXC_VERSION)
  target_link_libraries(aims PRIVATE libxc-custom)
endif()
# Subdirectories for which we cannot use add_subdirectory because they
# contain commands such as enable_language(), which can only be used
# in the topmost CMakeLists.txt (technically, where the aims target is
# defined).
include(DFPT_OpenCL/CMakeLists.txt)
include(gpu/CMakeLists.txt)
target_link_libraries(aims PRIVATE ${LIBS})
# Source files in the root source directory
target_sources(aims1 PRIVATE
  mpi_tasks.f90
  synchronize_mpi_basic.f90
  aims_memory_tracking.f90
  mpe_constants.f90
  types.f90
  runtime_choices.f90
  constants.f90
  localorb_io.f90
  dimensions.f90
  lapack_wrapper.f90
  full_local_mat_lapack.f90
  evaluate_moment_p2.f90
  output_periodic_dipole_moment.f90
  tab_trigonom2.f90
  get_relative_coords.f90
  hartree_non_periodic_ewald.f90
  analytic_multipole_coefficients.f90
  hartree_potential_real_p0.f90
  hartree_potential_storage.f90
  general_function.f90
  octree_routines.f90
  debugmanager.f90
  wf_extrapolation.f90
  symmetry_reduce_k_points.f90
  gaussian_lattice_gamma.f90
  gaussian_lattice_sum.f90
  get_gaussian_Vq_coeff.f90
  get_gaussian_Vq.f90
  get_singularity_lifting_chi.f90
  get_relative_lv.f90
  get_net_torque_lv.f90
  get_aux_gaussian_basis_fns.f90
  limit_rho_change.f90
  get_inner_max.f90
  pulay_grad_store.f90
  pulay_store.f90
  heapsort_index.f90
  shrink_opt_auxil_basis.f90
  shrink_full_auxil_basis_aux.f90
  shrink_svd_auxil_basis.f90
  wf_save_overlap.f90
  ball_touches_zero_cell.f90
  k_point_symmetry_check.f90
  stop_illconditioning.f90
  numerical_stress.f90
  bspline.f90
  precondition.f90
  trust_radius_method.f90
  debug_output.f90
  lindh.f90
  scalapack_matrix_type.f90
  matrix_inversion.f90
  symmetry_constrained_relaxation.f90
  distributed_hessian.f90
  energy_density.f90
  sparse_tensor.f90
  hartree_fock.f90
  generate_aims_uuid.f90
  elsi_wrapper.f90
  density_matrix_evaluation.f90
  scalapack_utils.f90
  synchronize_mpi.f90
  prodbas.f90
  cg_scalapack.f90
  basbas_fn_coulomb.f90
  scalapack_wrapper.f90
  restart_elsi.f90
  mpi_utilities.f90
  triple_Y.f90
  load_balancing.f90
  analyze_arrays.f90
  applicable_citations.f90
  timing_core.f90
  timing.f90
  xc_library.f90
  physics.f90
  check_environment.f90
  Hartree_F_p_functions.f90
  numerical_utilities.f90
  bravais.f90
  geometry.f90
  pbc_lists.f90
  )
target_sources(aims8 PRIVATE
  zeroin.f
  )
target_sources(aims7 PRIVATE
  my_triple_Y.f90
  lpb_solver_utilities.f90
  tight_binding_auxmat.f90
  relaxation.f90
  sbt_overlap_aims.f90
  lvl_triples.f90
  pseudodata.f90
  sbt_overlap_tb.f90
  sbt_overlap.f90
  logsbt_fast_kernel.f90
  my_logsbt.f90
  cut_coulomb_operator.f90
  logsbt.f90
  initialize_pulay_matrix_p1.f90
  xc.f90
  pulay_mix_p1.f90
  broyden_mix_p1.f90
  mixing.f90
  hartree_fock_p0.f90
  sym_base.f90
  analytical_stress.f90
  hartree_potential_recip.f90
  hartree_solver.f90
  precondition_dielectric.f90
  )
target_sources(aims2 PRIVATE
  calculate_fock_matrix_p0.f90
  )
target_sources(aims3 PRIVATE
  read_fixed_grid.f90
  mpb_solver_utilities.f90
  DFPT_pulay_mixing.f90
  thermodynamic_integration.f90
  check_mpi_in_place.f90
  set_gaussian_free_cut.f90
  read_plot_band.f90
  read_plot_dielectric.f90
  read_plot_band_during_scf.f90
  read_specified_grid.f90
  read_species_data.f90
  read_mpb_data.f90
  TDI_change_schedule_step.f90
  evaluate_KS_density_v2.f90
  evaluate_KS_orbital_density.f90
  cuba_stub.f90
  get_net_torque.f90
  binomial_coefficient.f90
  prune_basis_v2.f90
  prune_radial_basis_v2.f90
  tab_local_geometry.f90
  evaluate_radial_functions_v2.f90
  tab_trigonom_v2.f90
  tab_wave_ylm_v2.f90
  evaluate_waves_v2.f90
  get_net_force_lv.f90
  get_net_force.f90
  TDI_map_to_center_cell.f90
  verify_angular_grid.f90
  plot.f90
  hirshfeld.f90
  mbd_std_interface.f90
  mbd_std_wrapper.f90
  mbd_dev_interface.f90
  mbd_dev_wrapper.f90
  linalg.f90
  fodft.f90
  xml_write.f90
  friction.f90
  mod_crpa_blacs.f90
  esp_grids.f90
  esp_charges.f90
  molecular_dynamics.f90
  pi_molecular_dynamics.f90
  transport.f90
  separate_core_states.f90
  mixing_constraint.f90
  bfgs.f90
  force_occupation.f90
  cartesian_ylm.f90
  ll_vdwdf.f90
  contour_def_gw_types.f90
  gw_para.f90
  vdw_correction.f90
  heat_flux.f90
  plus_u.f90
  c_helper.f90
  aims_gpu.f90
  statistics.f90
  batch_statistics.f90
  json_output.f90
  tetra.f90
  output_mulliken.f90
  constraint.f90
  control_file.f90
  )
target_sources(aims4 PRIVATE
  read_control.f90
  )
target_sources(aims5 PRIVATE
  out_plot_band_soc_mulliken_scalapack.f90
  scs_cfdm.f90
  scf_solver.f90
  mbd_rsscs_numerical_forces.f90
  add_const_to_hamiltonian.f90
  add_gradient_part_to_H_p0.f90
  add_nonlocal_pot.f90
  add_zora_gradient_part_p0.f90
  add_zora_matrix_p1.f90
  adjust_scf_settings.f90
  advance_KS_solution.f90
  analy_continue_self_energy.f90
  analy_continue_self_energy_p0.f90
  asym_inv_sqrt_of_auxmat_scalapack.f90
  atom_bsse_results.f90
  aufbau_principle.f90
  auxiliary_matrix_multi.f90
  baswave_pp_overlap2.f90
  baswave_pp_overlap3.f90
  baswave_pp_overlap.f90
  binomial_coefficient_neg.f90
  broyden_store.f90
  calculate_coulmat_ovl.f90
  calculate_dipolemat.f90
  calculate_dipolemat_k_k.f90
  calculate_kinetic_density.f90
  calculate_kinetic_density_on_arbitrary_points.f90
  calculate_mommat_base.f90
  calculate_whole_potential_on_arbitrary_points.f90
  cg.f90
  check_consistency_of_keywords.f90
  check_cpu_consistency.f90
  check_geometry_convergence.f90
  check_norm_constraint.f90
  check_norm.f90
  check_norm_p0.f90
  check_norm_periodic_v2.f90
  check_norm_v2.f90
  check_n_states.f90
  check_occs.f90
  check_shell_convergence_p0.f90
  classical_field.f90
  collect_batch_centers_p2.f90
  compare_forces.f90
  compute_g_times_w_p0.f90
  compute_g_times_w_single_kpoint.f90
  condense_basis_pairs.f90
  construct_hamiltonian_and_ovl.f90
  construct_hamiltonian.f90
  construct_overlap.f90
  construct_xc_matr_kspace.f90
  contour_def_gw.f90
  contour_def_gw_environment.f90
  convert_l_str.f90
  count_frozen_core_states.f90
  create_ini_type_lookup.f90
  crpapw.f90
  cube_output.f90
  density_matrix_mixing.f90
  density_response.f90
  determine_center_of_molecule.f90
  determine_irreducible_k_grid.f90
  determine_k_minus_q_list.f90
  determine_k_minus_q_list_finemesh.f90
  determine_vacuum_level.f90
  dftpt2_dft_part_energy.f90
  diagonalize_auxmat_lapack.f90
  diagonalize_auxmat_scalapack_complex.f90
  diagonalize_auxmat_scalapack_supercell.f90
  dielectric_analysis.f90
  directories.f90
  distribute_irreducible_eigenvectors.f90
  dont_stop_no_parser.f90
  double_factorial.f90
  embedding_potential.f90
  esp_grid_storage.f90
  esp_partition_grid.f90
  estimate_Ewald_radius.f90
  estimate_n_compute_maxes.f90
  evaluate_2oex_energy.f90
  evaluate_2ox_selfenergy.f90
  evaluate_all_basis_waves.f90
  evaluate_bare_ci.f90
  evaluate_C6_coeff.f90
  evaluate_coulombhole_shell.f90
  evaluate_coulomb_matr_shell_p0.f90
  evaluate_coulomb_matr_shell_v0.f90
  evaluate_crpa_energy_kspace.f90
  evaluate_crpa_energy_kspace_restart.f90
  evaluate_density_gradient_complex_p1.f90
  evaluate_density_gradient_denmat.f90
  evaluate_density_gradient_p1.f90
  evaluate_densmat_part.f90
  evaluate_dipole_polarisability.f90
  evaluate_eex_energy.f90
  evaluate_ex_and_xc_matr_kspace.f90
  evaluate_exchange_energy.f90
  evaluate_exchange_matr_kspace_p0.f90
  evaluate_exchange_matr_kspace_single_kpoint_p0.f90
  evaluate_exchange_matr_LVL_eigen.f90
  evaluate_exchange_matr_v0.f90
  evaluate_exchange_matr_v1.f90
  evaluate_ext_charge_hellman_feynman_forces_p0.f90
  evaluate_exx_matr_kspace.f90
  evaluate_fitting_function.f90
  evaluate_free_atom_sums_p0.f90
  evaluate_free_atom_sums_p2.f90
  evaluate_free_rho_sums_p0.f90
  evaluate_gga_forces_p0.f90
  evaluate_gw_selfenergy_band_kpoint.f90
  evaluate_hamiltonian_shell_p1.f90
  evaluate_hellman_feynman_forces_p0.f90
  evaluate_h_minus_e_times_psi_v2.f90
  evaluate_H_psi_p0.f90
  evaluate_H_psi_p2.f90
  evaluate_invs_times_sqrtv.f90
  evaluate_invs_times_sqrtv_scalapack.f90
  evaluate_kinetic_density_densmat.f90
  evaluate_KS_density_complex_p0.f90
  evaluate_KS_density_densmat.f90
  evaluate_KS_density_p0.f90
  evaluate_KS_orbital_density_complex_p0.f90
  evaluate_KS_orbital_density_p1.f90
  evaluate_density_kinetic_p1.f90
  evaluate_density_kinetic_complex_p1.f90
  evaluate_KS_orbital_gradients_complex_p1.f90
  evaluate_KS_orbital_gradients_p1.f90
  evaluate_KS_split_xc_matrix.f90
  evaluate_KS_xc_matrix.f90
  evaluate_localpot_mult_pseudoatoms_v2.f90
  evaluate_loc_prod_waves_p0.f90
  evaluate_mgga_contribution_and_add_to_hamiltonian_shell.f90
  evaluate_moment_p1.f90
  evaluate_mp2_correlation_energy.f90
  evaluate_mp2_selfenergy.f90
  evaluate_nlcc_correction_energy.f90
  evaluate_nlcc_forces.f90
  evaluate_nonlocal_energy.f90
  evaluate_nuclear_gradients_p0.f90
  evaluate_nuclear_potential.f90
  evaluate_O_2bs1HF.f90
  evaluate_O_bbbsKS_part.f90
  evaluate_orb_grad_dot_rho_grad_p0.f90
  evaluate_orb_hess_dot_rho_grad_p0.f90
  evaluate_osmp2_correlation_energy.f90
  evaluate_ovlp3fn_shell_v0.f90
  evaluate_ovlp_3MO_EN.f90
  evaluate_ovlp_3MO.f90
  evaluate_ovlp_shell_p0.f90
  evaluate_partition_p2.f90
  evaluate_partition_tab_2atoms.f90
  evaluate_partition_tab_2atoms_gen_p0.f90
  evaluate_partition_tab_2atoms_p0.f90
  evaluate_partition_tab_3atoms.f90
  evaluate_partition_tab_3atoms_p0.f90
  evaluate_partition_tab_p0.f90
  evaluate_partition_tab_p2.f90
  evaluate_periodic_gw_selfenergy.f90
  evaluate_polarisability_freq.f90
  evaluate_polarisability_kspace.f90
  evaluate_post_xc.f90
  evaluate_potential_cube.f90
  evaluate_pot_superpos_p0.f90
  evaluate_pp_hellman_feynman_forces_p0.f90
  evaluate_pp_nlcc_forces.f90
  evaluate_pp_nonlocal_forces.f90
  evaluate_prod_waves.f90
  evaluate_prod_waves_p0.f90
  evaluate_pulay_forces_p0.f90
  evaluate_quadrupole_moment.f90
  evaluate_radial_functions_deriv_p0.f90
  evaluate_radial_functions_p0.f90
  evaluate_renormalized_single_excitation_correction.f90
  evaluate_rho_from_densmat.f90
  evaluate_rho_gradient.f90
  evaluate_rho_gradient_from_densmat.f90
  evaluate_rho_multipole_gradient.f90
  evaluate_rpa_correlation_energy.f90
  evaluate_rpa_integrand_along_ac_path.f90
  evaluate_rpa_integrand.f90
  evaluate_scaled_KS_orbitals_v1.f90
  evaluate_screx_matrix.f90
  evaluate_self_energy_freq.f90
  evaluate_single_excitation_correction.f90
  evaluate_single_excitation_correction_p0.f90
  evaluate_single_excitation_to_self_energy.f90
  evaluate_sosex_2w_selfenergy.f90
  evaluate_sosex_energy.f90
  evaluate_sosex_selfenergy.f90
  evaluate_soxw_selfenergy.f90
  evaluate_T_psi.f90
  evaluate_v_hartree_and_rho_multipole_gradient.f90
  evaluate_v_hartree_gradient.f90
  evaluate_V_psi_p2.f90
  evaluate_T_plus_V_psi_p2.f90
  evaluate_v_times_waves.f90
  evaluate_v_times_waves_p0.f90
  evaluate_wavefunc_cmplx.f90
  evaluate_wavefunc.f90
  evaluate_wavefunc_real.f90
  evaluate_wave_gradient_cartesian_p1.f90
  evaluate_wave_gradient_cartesian_p2.f90
  evaluate_wave_gradient_p0.f90
  evaluate_wave_gradient_p2.f90
  evaluate_wave_gradient_prodbas_v1.f90
  evaluate_wave_gradient_v1.f90
  evaluate_wave_hessian_cartesian_p0.f90
  evaluate_wave_hessian_cartesian_p2.f90
  evaluate_wave_psi.f90
  evaluate_waves_mult_point_center_fn.f90
  evaluate_waves_p0.f90
  evaluate_waves_p2.f90
  evaluate_waves_v0.f90
  evaluate_xc_energy_shell.f90
  evaluate_xc.f90
  evaluate_xc_matr_kspace.f90
  evaluate_xc_matr_shell.f90
  evaluate_xc_psi_p2.f90
  evaluate_xc_split.f90
  evaluate_zora_vector_p1.f90
  exchange_ev.f90
  exchange_trico.f90
  external_pressure_on_system.f90
  factorial.f90
  f_external_constraint.f90
  final_deallocations.f90
  find_homo_lumo_gap.f90
  final_energy_output.f90
  forces_densmat.f90
  gather_auxmat.f90
  generate_basbasfn.f90
  gen_inv_auxmat_scalapack_2d.f90
  get_absorption.f90
  get_atomic_charges.f90
  get_coeff_3fn_lvl.f90
  get_coeff_3fn_lvl_full.f90
  get_coeff_3fn_svs.f90
  get_coeff_3fn_v_1d.f90
  get_coeff_3fn_v_2d.f90
  get_constraint_fermi.f90
  get_coulelement_lvl_v0.f90
  get_coulelement_ovl.f90
  get_coulomb_coeff_blacs.f90
  get_coulomb_matr_recip.f90
  get_dielectric_function.f90
  get_dipolemat_k_k.f90
  get_dipolematrix.f90
  get_distance.f90
  get_entropy_correction_p1.f90
  get_esp_charges.f90
  get_ev_scalapack_hdf5.f90
  get_ev_sum_p0.f90
  get_exchange_energy_p0.f90
  get_fnKSbb.f90
  get_fock_energy.f90
  get_free_superpos_energy_p1.f90
  get_gw_band_struct_info.f90
  get_gw_tot_en.f90
  get_hf_hamiltonian_complex_p0.f90
  get_hf_hamiltonian.f90
  get_hf_hamiltonian_real_p0.f90
  get_kinetic_energy.f90
  get_KS_orbitals_bandplot.f90
  get_lda_free_atom_densities.f90
  get_lvl_tricoeff_bravais.f90
  get_lvl_tricoeff_recip.f90
  get_lvl_tricoeff_recip_general.f90
  get_machine_precision.f90
  get_max_force.f90
  get_memory.f90
  get_momentummatrix.f90
  get_n_compute_maxes_p1.f90
  get_n_max_for_densities.f90
  get_occupation_numbers.f90
  get_occupation_numbers_fsm.f90
  get_occupation_numbers_p0.f90
  get_occupation_numbers_single_channel.f90
  get_occupation_numbers_v2.f90
  get_ovlp3fn_multi_sv.f90
  get_pade_approx_para.f90
  get_pairwise_coeff_3fn.f90
  get_penalty_energy.f90
  get_projected_dos.f90
  get_qp_spectrum.f90
  get_real_selfenergy.f90
  get_rho_mp.f90
  get_rho_mp_single_center.f90
  get_rlylm_and_derivs.f90
  get_rlylm.f90
  get_screx_hamiltonian.f90
  get_sym_overlap.f90
  get_system_symmetries.f90
  get_total_energy.f90
  get_total_forces.f90
  get_v_multi_ovlp3fn.f90
  get_xc.f90
  hdf5_output.f90
  hdf5_tools.f90
  heapsort.f90
  heapsort_general.f90
  hf_postproc.f90
  hirshfeld_analysis.f90
  hirshfeld_analysis_iterative.f90
  hse_logsbt_integrator.f90
  hse_logsbt_integrator_grid.f90
  ifcore_stub.f90
  ifp_dielectric_interface_plane.f90
  improve_complex_eigenfunctions.f90
  improve_eigenfunctions.f90
  increment_ylm_deriv.f90
  increment_ylm_deriv_forvdw.f90
  init_citations.f90
  init_debug.f90
  initialize_cube_grid_storage.f90
  initialize_fitting_parameter.f90
  initialize_fixed_grids.f90
  initialize_grid_storage.f90
  initialize_hartree_fock.f90
  initialize_hartree_fock_p0.f90
  initialize_integrals_p0.f90
  initialize_noscf.f90
  initialize_scf.f90
  inner_product.f90
  integrate_coulombhole.f90
  integrate_coulomb_matr_p1.f90
  integrate_coulomb_matr_v0.f90
  integrate_dipmom_pairstates.f90
  integrate_errorfunction.f90
  integrate_hamiltonian_matrix_p2.f90
  integrate_hartree_log_grid.f90
  integrate_kinetic_matrix.f90
  integrate_log_coulomb_metric.f90
  integrate_log_mesh.f90
  integrate_nlcorr.f90
  integrate_ovlp3fn.f90
  integrate_ovlp_matrix_p1.f90
  integrate_ovlp_matrix_p2.f90
  integrate_ovlp_prodbas.f90
  integrate_post_xc_energy.f90
  integrate_split_xc_matrix.f90
  integrate_v_external_matrix.f90
  integrate_t_plus_v_matrix.f90
  integrate_v_psi.f90
  integrate_v_times_radialwaves.f90
  integrate_v_times_radialwaves_v2.f90
  integrate_xc_energy.f90
  integrate_xc_matrix.f90
  integrate_xc_realspace_p2.f90
  ion_pseudopot_potential_v1.f90
  ion_pseudopot_potential_v2.f90
  isc_implicit_solvent_cavity.f90
  isc_constraint_dynamics.f90
  isc_common.f90
  isc_projected_voronoi.f90
  KH_core_states.f90
  KS_optical_properties.f90
  KS_optical_properties_tetrahedron.f90
  ks_wrapper.f90
  lapack_solver.f90
  lc_wpbeh.f90
  localized_basbas.f90
  log_spline_blender.f90
  look_at_rho.f90
  lowdin_analysis.f90
  lpb_solver.f90
  map_to_center_cell.f90
  map_to_first_octant.f90
  MD_QH_init.f90
  meta_gga_postp.f90
  metis_stub.f90
  mp2_calculation.f90
  mpe_interface.f90
  mpe_reaction_field.f90
  mpe_types.f90
  mpe_dielectric_interfaces.f90
  mpe_dielectric_interfaces_common.f90
  mpi_routines_stub.f90
  mr_min_lsq.f90
  mulliken_analysis.f90
  multi_ovlp3KS_sqrtw.f90
  my_lvl_triples.f90
  nlcorr.f90
  noscf_solver.f90
  obtain_initial_dimensions.f90
  omp_stub.f90
  orthonormalize_eigenvectors.f90
  orthonormalize_prodbas_fn.f90
  WannierCenters.f90
  out_plot_band.f90
  out_plot_band_hf_k_space.f90
  out_plot_band_mulliken.f90
  out_plot_band_mulliken_scalapack.f90
  out_plot_band_scalapack.f90
  out_plot_band_soc_mulliken.f90
  output_bands_during_scf.f90
  output_bxsf.f90
  output_cube_files_p1.f90
  output_cube_files_p2.f90
  output_delta_v.f90
  output_delta_v_step.f90
  output_density_p1.f90
  output_dielec_func.f90
  output_dipole_moment.f90
  output_effective_potential.f90
  output_eigenfunctions.f90
  output_eigenvec_ovlp_p1.f90
  output_energy_and_forces.f90
  output_ev_sum.f90
  output_gradient_basis.f90
  output_handling.f90
  output_hartree_pot.f90
  output_ka_overlap.f90
  output_KS_dos.f90
  output_KS_dos_pert.f90
  output_KS_dos_tetrahedron.f90
  output_lowdin.f90
  output_matrices_aij.f90
  output_matrices.f90
  output_molecular_dynamics.f90
  output_pi_molecular_dynamics.f90
  output_potential_p1.f90
  output_potential_p1_new.f90
  output_quadrupole_moment.f90
  output_rho.f90
  output_rho_multipole.f90
  output_self_energy.f90
  output_sorted_eigenvalues.f90
  output_v_external_matrix.f90
  output_whole_potential_on_arbitrary_points.f90
  ovlp3KS_lvl.f90
  packed_matrix_transform_complex.f90
  packed_matrix_transform.f90
  partial_core_dens2grid_v2.f90
  partition_grid.f90
  perform_symmetry_analysis.f90
  perf_stubs.f90
  physics_ext.f90
  plumed_new_interface.f90
  plumed_stub.f90
  post_scf_correlation_treatment_p0.f90
  pot_constrain.f90
  power_auxmat_complex.f90
  power_auxmat_lapack.f90
  power_auxmat_lapack_supercell.f90
  power_auxmat_scalapack_complex.f90
  power_auxmat_scalapack.f90
  power_auxmat_scalapack_supercell.f90
  power_genmat_lapack.f90
  predict_new_geometry.f90
  prepare_corr_energy_calc.f90
  prepare_embedding_grids.f90
  prepare_fine_k_quantities.f90
  prepare_partition_tabs.f90
  prepare_scf.f90
  prepare_system.f90
  print_illconditioning.f90
  print_self_energies.f90
  prune_basis_p0.f90
  prune_basis_p2.f90
  prune_basis_v1.f90
  prune_density_matrix.f90
  prune_density_matrix_sparse.f90
  prune_force_atoms_v1.f90
  prune_general_basis.f90
  prune_prodbas_v1.f90
  prune_radial_basis_p0.f90
  prune_radial_basis_p2.f90
  qpe_calculation.f90
  quasi_particle_energies.f90
  quasi_particle_energy_multisolu.f90
  quasi_particle_energy_p0.f90
  quicksort_index.f90
  radial_fourier_hse.f90
  read_control_update.f90
  read_geo.f90
  read_input_data.f90
  read_MD_restart.f90
  read_multiaims.f90
  read_multipole_data.f90
  read_PIMD_restart.f90
  reinitialize_noscf.f90
  reinitialize_scf.f90
  renormalize_density.f90
  report_out_of_walltime_error.f90
  reset_grid_partitions.f90
  reset_hamiltonian.f90
  restart.f90
  restart_rpa.f90
  rho_multipole_evaluation.f90
  rpa_calculation_p0.f90
  run_mpb_solver.f90
  run_newton_mpb.f90
  scalapack_generic_wrapper.f90
  scaled_zora_transform.f90
  scatter_auxmat.f90
  screened_coulomb_interaction.f90
  setup_cube_defaults.f90
  shm_stub.f90
  shrink_full_auxil_basis.f90
  shrink_full_auxil_basis_v2.f90
  smooth_partition_edge.f90
  solve_KS_eigen.f90
  spectral_func_cd.f90
  SPE_solver.f90
  spline_blender.f90
  stratmann_partition_weight.f90
  stratmann_weight_restricted.f90
  sum_up_whole_potential_p1.f90
  tab_atom_centered_coords.f90
  tab_atom_centered_coords_p0.f90
  tab_atom_centered_coords_p1.f90
  tab_atom_centered_coords_v2.f90
  tab_global_geometry_p0.f90
  tab_gradient_ylm.f90
  tab_gradient_ylm_p0.f90
  tab_interatomic_distances.f90
  tab_local_geometry_p0.f90
  tab_local_geometry_p2.f90
  tab_single_atom_centered_coords_p0.f90
  tab_single_atom_centered_coords_radial_log_p0.f90
  tab_single_gradient_ylm_p0.f90
  tab_single_gradient_ylm_p2.f90
  tab_single_gradient_ylm_p2_forvdw.f90
  tab_single_radial_weights_v2.f90
  tab_single_trigonom_p0.f90
  tab_single_wave_ylm_p0.f90
  tab_single_wave_ylm_p2.f90
  tab_trigonom.f90
  tab_trigonom_p0.f90
  tab_wave_ylm.f90
  tab_wave_ylm_p0.f90
  tab_wave_ylm_p1.f90
  tab_ylm.f90
  tddft_real_time_propagation.f90
  TDI_update_lambda.f90
  test_pdtran.f90
  tf_ini.f90
  threshold_occ_numbers_constraint.f90
  threshold_occ_numbers.f90
  transform_ovlp3fn.f90
  unpack_matrices.f90
  unshift_eigenvalues.f90
  update_density_and_forces_orbital.f90
  update_density_densmat.f90
  update_full_matrix_p0.f90
  update_grid_coords.f90
  update_hartree_potential_p1.f90
  update_missing_density_densmat.f90
  update_missing_density_orbital.f90
  v_hartree_multipole_evaluation.f90
  VTune_stubs.f90
  write_cube_header.f90
  write_lattice_parameters.f90
  zeroin_constraint.f90
  zeroin_p0.f90
  zeroin_periodic_v2.f90
  zeroin_v2.f90
  write_version_stamp_stub.f90
  evaluate_wavefunc_and_grad_real.f90
  evaluate_wavefunc_and_grad_cmplx.f90
  evaluate_wavefunc_and_grad.f90
  output_dgrid.f90
  ${CMAKE_CURRENT_BINARY_DIR}/cmake_info.f90
  set_aims_defaults.f90
  set_mpe_defaults.f90
  write_preamble.f90
  DFPT_wrapper.f90
  )
