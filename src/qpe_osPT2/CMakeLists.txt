# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_ospt2_integrand.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_ospt2_polarizability_freq.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_osrpa_integrand.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_osrpa_integrand_2.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_qpe_ospt2_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_qpe_osrpa_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_qpe_sicrpa_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_special_radius_x0.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_special_radius_x2.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_sicrpa_integrand.f90
  ${CMAKE_CURRENT_LIST_DIR}/logarithmic_matrix_lapack.f90
  ${CMAKE_CURRENT_LIST_DIR}/diagonalize_auxmat_scalapack_real.f90
  ${CMAKE_CURRENT_LIST_DIR}/qpe_osPT2.f90)
