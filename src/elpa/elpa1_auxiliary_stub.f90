!****h* FHI-aims/elpa/ELPA1_auxiliary_stub
!  NAME
!    ELPA1_auxiliary_stub
!  SYNOPSIS

module ELPA1_auxiliary_2013

!  PURPOSE
!  Contains dummy routines to enable linking to the ELPA Library
!
!  AUTHOR
!    Lydia Nemec, FHI-aims team, Chair for Theoretical Chemistry, Technical University Munich
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180, 2175-2196 (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!
!  Uses

   implicit none
!  SOURCE

! By default all variables below are private
  PRIVATE ! set default to private

! public routines
  public :: mult_at_b_real_2013             ! Multiply real matrices A**T * B
  public :: mult_ah_b_complex_2013          ! Multiply complex matrices A**H * B
  public :: cholesky_complex_2013           ! Cholesky factorization of a complex matrix
  public :: cholesky_real_2013              ! Cholesky factorization of a real matrix
  public :: invert_trm_complex_2013         ! Invert complex triangular matrix
  public :: invert_trm_real_2013            ! Invert real triangular matrix

  contains

  subroutine mult_ah_b_complex_2013(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

      implicit none

      character*1                   :: uplo_a, uplo_c
      integer                       :: na, ncb, lda, ldb, nblk, mpi_comm_rows, mpi_comm_cols, ldc
      complex*16                    :: a(lda,*), b(ldb,*), c(ldc,*) ! remove assumed size!

  end subroutine mult_ah_b_complex_2013


  subroutine mult_at_b_real_2013(uplo_a, uplo_c, na, ncb, a, lda, b, ldb, nblk, mpi_comm_rows, mpi_comm_cols, c, ldc)

      implicit none

      character*1                   :: uplo_a, uplo_c
      integer                       :: na, ncb, lda, ldb, nblk, mpi_comm_rows, mpi_comm_cols, ldc
      real*8                        :: a(lda,*), b(ldb,*), c(ldc,*) ! remove assumed size!

  end subroutine  mult_at_b_real_2013


  subroutine cholesky_complex_2013(na, a, lda, nblk, mxcol, mpi_comm_rows, mpi_comm_cols, wantDebug, success)
  
     implicit none

     integer,intent(IN)   :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols
     complex*16           :: a(:,:)
     logical              :: success
     logical              :: wantDebug

     integer mxcol 
     mxcol = SIZE(a,2)

  end subroutine cholesky_complex_2013

  subroutine cholesky_real_2013(na, a, lda, nblk, mxcol, mpi_comm_rows, mpi_comm_cols, wantDebug, success)
  
     implicit none

     integer,intent(IN)   :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols
     real*8               :: a(:,:)
     logical              :: success
     logical              :: wantDebug
     integer mxcol 
     mxcol = SIZE(a,2)

  end subroutine cholesky_real_2013


  subroutine invert_trm_complex_2013(na, a, lda, nblk, mxcol, mpi_comm_rows, mpi_comm_cols, wantDebug, success)

     implicit none

     integer,intent(IN)  :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols
     complex*16          :: a(:,:)
     logical             :: success
     logical             :: wantDebug
     integer mxcol 
     mxcol = SIZE(a,2)

  end subroutine invert_trm_complex_2013


  subroutine  invert_trm_real_2013(na, a, lda, nblk, mxcol, mpi_comm_rows, mpi_comm_cols, wantDebug, success)

     implicit none

     integer,intent(IN)  :: na, lda, nblk, mpi_comm_rows, mpi_comm_cols
     real*8              :: a(:,:)
     logical             :: success
     logical             :: wantDebug
     integer mxcol 
     mxcol = SIZE(a,2)

  end subroutine


end module ELPA1_auxiliary_2013
