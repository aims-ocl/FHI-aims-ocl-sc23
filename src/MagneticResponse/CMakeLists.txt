# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/tools.f90
  ${CMAKE_CURRENT_LIST_DIR}/psi_at_nucleus_mod.f90
  )
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/compute_dipolar_couplings.f90
  ${CMAKE_CURRENT_LIST_DIR}/DFPT.f90
  ${CMAKE_CURRENT_LIST_DIR}/GIAO_mult_psi.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrands.f90
  ${CMAKE_CURRENT_LIST_DIR}/integrate_FC.f90
  ${CMAKE_CURRENT_LIST_DIR}/integration.f90
  ${CMAKE_CURRENT_LIST_DIR}/MR_core.f90
  ${CMAKE_CURRENT_LIST_DIR}/MR_global.f90
  ${CMAKE_CURRENT_LIST_DIR}/MR_main.f90
  ${CMAKE_CURRENT_LIST_DIR}/MR_output.f90
  ${CMAKE_CURRENT_LIST_DIR}/nuclear_data.f90
  ${CMAKE_CURRENT_LIST_DIR}/xc_kernel_rho.f90
  ${CMAKE_CURRENT_LIST_DIR}/xc_wrapper.f90)
