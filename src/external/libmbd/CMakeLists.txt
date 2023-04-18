# Copyright: Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
# e.V. Please note that any use of the "FHI-aims-Software" is subject
# to the terms and conditions of the respective license agreement.

set(MBDDIR "${CMAKE_CURRENT_LIST_DIR}/src")

target_sources(aims1 PRIVATE
  ${MBDDIR}/mbd_constants.f90
  ${MBDDIR}/mbd_mpi.f90
  ${MBDDIR}/mbd.f90
  ${MBDDIR}/mbd_blacs.f90
  ${MBDDIR}/mbd_damping.f90
  ${MBDDIR}/mbd_dipole.f90
  ${MBDDIR}/mbd_elsi.f90
  ${MBDDIR}/mbd_formulas.f90
  ${MBDDIR}/mbd_geom.f90
  ${MBDDIR}/mbd_gradients.f90
  ${MBDDIR}/mbd_hamiltonian.f90
  ${MBDDIR}/mbd_lapack.f90
  ${MBDDIR}/mbd_linalg.f90
  ${MBDDIR}/mbd_matrix.f90
  ${MBDDIR}/mbd_methods.f90
  ${MBDDIR}/mbd_rpa.f90
  ${MBDDIR}/mbd_scalapack.f90
  ${MBDDIR}/mbd_scs.f90
  ${MBDDIR}/mbd_ts.f90
  ${MBDDIR}/mbd_utils.f90
  ${MBDDIR}/mbd_vdw_param.f90
)
