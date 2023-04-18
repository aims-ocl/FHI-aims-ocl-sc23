# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/construct_first_order_matrix_k_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_matrix_complex.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_matrix_complex_no_symmetry.f90
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/debug_compare.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Born_effective_charges_Pulay.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_dielectric_constant.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_EDM_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_gradient_rho_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_Q_k_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_momentum_matrix_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Omega_MO_diag.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Omega_MO_v1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Omega_MO_v2.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Omega_part_1_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_r_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_DM_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_EDM_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_momentummatrix_B.f90
  ${CMAKE_CURRENT_LIST_DIR}/get_momentum_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_Born_effective_charges_Hellman_Feynman.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_Born_effective_charges_MP.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_Born_effective_charges_MP_force.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_Born_effective_charges_Pulay.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_momentum_matrix_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_Omega_part_1_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_polarizability_dielectric_constant.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_r_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/prune_density_matrix_sparse_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_shanghui_dielectric.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_shanghui_dielectric.f90)
