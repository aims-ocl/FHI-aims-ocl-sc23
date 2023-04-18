# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

target_include_directories(aims1 PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
if (USE_IPC)
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ipc.f90
    ${CMAKE_CURRENT_LIST_DIR}/ipc_c.c
    ${CMAKE_CURRENT_LIST_DIR}/ipc_packet.c
    ${CMAKE_CURRENT_LIST_DIR}/ipc_packet_test.c
    ${CMAKE_CURRENT_LIST_DIR}/ipc_stdio.c
    ${CMAKE_CURRENT_LIST_DIR}/ipc_stdio_test.c
    ${CMAKE_CURRENT_LIST_DIR}/ipc_test.c)
else()
  target_sources(aims1 PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/ipc_stub.f90)
endif()
target_sources(aims5 PRIVATE
  ${CMAKE_CURRENT_LIST_DIR}/ipc_hamiltonian_and_ovl_transfer.f90)
