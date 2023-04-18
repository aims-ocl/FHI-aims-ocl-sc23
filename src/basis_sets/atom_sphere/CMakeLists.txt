# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_sources(aims1 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/allcommon.f90
  ${CMAKE_CURRENT_LIST_DIR}/atomsolver.f90
  ${CMAKE_CURRENT_LIST_DIR}/zorasolver.f90
  ${CMAKE_CURRENT_LIST_DIR}/driveXC.f90
  ${CMAKE_CURRENT_LIST_DIR}/xcfunction.f90)
