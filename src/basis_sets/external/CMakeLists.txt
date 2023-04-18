# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims6 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/anginc.f
  ${CMAKE_CURRENT_LIST_DIR}/atomini.f
  ${CMAKE_CURRENT_LIST_DIR}/Lebedev-Laikov.f
  ${CMAKE_CURRENT_LIST_DIR}/HeoXu.f)

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/estd.f90)
