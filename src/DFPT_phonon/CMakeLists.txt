# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims3 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/add_matrix_to_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_phonon_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/DFPT_phonon_supercell.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_supercell_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_EDM_supercell_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_pulay_hessian_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_DM_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_zero_order_EDM_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_delta_hellman_hessian_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_hessian_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_sum_gradient_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_pulay_hessian_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/prune_first_order_density_matrix_sparse.f90
  ${CMAKE_CURRENT_LIST_DIR}/read_plot_DFPT_phonon_band.f90
  ${CMAKE_CURRENT_LIST_DIR}/solve_eigen_supercell.f90
  ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_shanghui_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/tab_two_center_coords_PBC.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_circle_to_cubic.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_cubic_to_circle.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_first_order_sparse_to_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_hessian_to_dynamical_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_second_order_sparse_to_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/trans_sparse_to_matrix.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_full_matrix_p1.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_full_matrix_p1_for_fd_benchmark.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_at_zero.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_shanghui_p1.f90)
