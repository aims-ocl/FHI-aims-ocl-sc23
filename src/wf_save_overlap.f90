!----------------------------------------------------------------------------
!****s* FHI-aims/wf_save_overlap
!  NAME
!    wf_save_overlap
!  SYNOPSIS

subroutine wf_save_overlap(n_dim1, n_dim2, ovlp_scalapack)

  !  PURPOSE
  !
  !     Saves the overlap matrix to module data.
  !
  !     In this module, the plain overlap matrix is used.  In
  !     scalapack_wrapper.f90, the overlap matrix is overwritten by its
  !     factorization.  Therefore, this routine is called at the right point
  !     of time from scalapack_wrapper.f90.
  !
  !     Not in wf_extrapolation.f90 because scalapack_wrapper cannot use
  !     wf_extrapolation.f90.
  !
  !  USES

  use runtime_choices
  use mpi_tasks
  use wf_extrapolation, only: wf_scalapack_overlap
  use scalapack_wrapper, only: set_full_matrix_real
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_dim1, n_dim2
  real*8, intent(IN) :: ovlp_scalapack(n_dim1, n_dim2)

  !  INPUTS
  !    o n_dim1, n_dim2 -- sizes of local storage
  !    o ovlp_scalapack -- distributed overlap matrix in 'U' format
  !  OUTPUTS
  !    none
  !    wf_scalapack_overlap is set to the full overlap matrix
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'wf_save_overlap'

  if (.not. real_eigenvectors) call aims_stop('.not. real_eigenvectors', func)

  wf_scalapack_overlap = ovlp_scalapack

  call set_full_matrix_real(wf_scalapack_overlap)

end subroutine wf_save_overlap
!******
!----------------------------------------------------------------------------
!****s* FHI-aims/wf_save_overlap_cmplx
!  NAME
!    wf_save_overlap_cmplx
!  SYNOPSIS

subroutine wf_save_overlap_cmplx(n_dim1, n_dim2, ovlp_scalapack_cmplx)

  !  PURPOSE
  !
  !     See wf_save_overlap() for the complex case.
  !
  !  USES

  use runtime_choices
  use mpi_tasks
  use wf_extrapolation, only: wf_scalapack_overlap_cmplx
  use scalapack_wrapper, only: set_full_matrix_complex
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_dim1, n_dim2
  complex*16, intent(IN) :: ovlp_scalapack_cmplx(n_dim1, n_dim2)

  !  INPUTS
  !    o n_dim1, n_dim2 -- sizes of local storage
  !    o ovlp_scalapack -- distributed overlap matrix in 'U' format
  !  OUTPUTS
  !    none
  !    wf_scalapack_overlap is set to the full overlap matrix
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'wf_save_overlap_cmplx'

  if (real_eigenvectors) call aims_stop('real_eigenvectors', func)

  wf_scalapack_overlap_cmplx = ovlp_scalapack_cmplx

  call set_full_matrix_complex(wf_scalapack_overlap_cmplx)

end subroutine wf_save_overlap_cmplx
!******
