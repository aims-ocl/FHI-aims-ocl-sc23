! --------------------------------------------------------------------------------------------------
!
! This file contains a stub routine for the elpa2 kernels if ELPA with AVX optimizations are
! used
!
! Reason: ELPA AVX needs to compile c and c++ files, however, the FHI-aims makefile nevertheless
! expects to compile a FORTRAN file for ELPA kernels
!
! --------------------------------------------------------------------------------------------------



subroutine elpa2_kernels_real_stubs

  implicit none

  return

end subroutine elpa2_kernels_real_stubs
