# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/poles_fit.f90
  ${CMAKE_CURRENT_LIST_DIR}/scgw_grid.f90
  ${CMAKE_CURRENT_LIST_DIR}/poles.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/analy_continue_green_fn.f90
  ${CMAKE_CURRENT_LIST_DIR}/check_the_error.f90
  ${CMAKE_CURRENT_LIST_DIR}/compute_scgw_dipole.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_green_fn.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_self_en.f90
  ${CMAKE_CURRENT_LIST_DIR}/E_kin.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_ovlp_NAO_KS.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_exchange_self_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_green_function_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_green_function_time.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_hartree_pot.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_polar.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_real_space_density.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_self_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_spectrum.f90
  ${CMAKE_CURRENT_LIST_DIR}/gt.f90
  ${CMAKE_CURRENT_LIST_DIR}/invert_green.f90
  ${CMAKE_CURRENT_LIST_DIR}/invert_overlap_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/read_and_write_G.f90
  ${CMAKE_CURRENT_LIST_DIR}/scgw_allocations.f90
  ${CMAKE_CURRENT_LIST_DIR}/scgw_densmat_analysis.f90
  ${CMAKE_CURRENT_LIST_DIR}/self_consistent_gw.f90
  ${CMAKE_CURRENT_LIST_DIR}/single_particle_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_dyson_equation.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_dyson_equation_re.f90
  ${CMAKE_CURRENT_LIST_DIR}/spectrum.f90
  ${CMAKE_CURRENT_LIST_DIR}/test_exchange_self_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/test_self_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/total_energy_calculation.f90
  ${CMAKE_CURRENT_LIST_DIR}/transform_and_interpolate_green_fn_to_time.f90
  ${CMAKE_CURRENT_LIST_DIR}/transform_exchange_energy.f90)
