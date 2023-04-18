# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/calculate_first_order_S_and_H_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_first_order_matrix_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_first_order_S_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_Born_effective_charges_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_EDM_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_gradient_rho_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_GGA_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_hellman_feynman_delta_part_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_pulay_dynamical_matrix_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_DM_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_EDM_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_dynamical_matrix_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_dynamical_matrix_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_sum_gradient_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_hellam_dynamical_matrix_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/prune_density_matrix_sparse_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_shanghui_phonon_reduce_memory.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_shanghui_phonon_reduce_memory.f90)
