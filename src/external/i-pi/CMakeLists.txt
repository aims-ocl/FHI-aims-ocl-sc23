# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

if (USE_iPI)
  target_sources(aims5 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/run_driver.f90
    ${CMAKE_CURRENT_LIST_DIR}/sockets.c)
endif()
