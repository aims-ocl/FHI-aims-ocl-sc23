# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

if (NOT USE_LIBXC)
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/libxc_stubs.f90)
endif()
