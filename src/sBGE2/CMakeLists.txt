# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_en_mp2_correlation_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/nsolver_iepa_mp2.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_iepa_mp2_correlation_energy.f90
  ${CMAKE_CURRENT_LIST_DIR}/nsolver_iepa_mp2_v02.f90
  ${CMAKE_CURRENT_LIST_DIR}/evaluate_iepa_mp2_correlation_energy_v02.f90)
