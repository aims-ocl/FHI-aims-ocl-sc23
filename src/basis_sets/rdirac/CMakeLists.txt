# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/fderiv.f90
  ${CMAKE_CURRENT_LIST_DIR}/polynom.f90
  ${CMAKE_CURRENT_LIST_DIR}/rdiracdme.f90
  ${CMAKE_CURRENT_LIST_DIR}/rdirac.f90
  ${CMAKE_CURRENT_LIST_DIR}/rdiracint.f90
  ${CMAKE_CURRENT_LIST_DIR}/spline.f90
  ${CMAKE_CURRENT_LIST_DIR}/dftatom.f90)
