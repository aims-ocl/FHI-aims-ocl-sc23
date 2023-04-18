# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/construct_first_order_H_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/construct_first_order_S_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_phonon_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_pulay_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_S_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_zero_order_DM_EDM_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_H_pulay_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_second_order_S_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_S_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_hessian_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_free_atom_sum_gradient_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_pulay_hessian_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/sum_up_whole_potential_shanghui_p0.f90
  ${CMAKE_CURRENT_LIST_DIR}/update_hartree_potential_shanghui_p0.f90)
