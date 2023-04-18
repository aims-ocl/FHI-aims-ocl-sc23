# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/add_HF_to_first_order_H_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/cpscf_solver_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_DM_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_gradient_rho_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_H_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_rho_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_first_order_U_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_xc_DFPT.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_H_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_first_order_rho_polarizability.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_polarizability.f90)
