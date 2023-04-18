!****h* FHI-aims/matrix_inversion
!  NAME
!    Routines for performing matrix inversion of general
!    (not necessarily positive definite) matrices
!  SYNOPSIS

module matrix_inversion

   use scalapack_matrix_type,      only: matrix_structure_type
   use mpi_tasks,                  only: aims_stop
 
   implicit none

   private
   public :: invmat_real_lapack, invmat_cmplx_lapack, invmat_real_scalapack,&
             invmat_cmplx_scalapack
   public :: invert_real_matrix, invert_cmplx_matrix, set_matrix_elements


contains

! **************************************************************************************************
!> \brief returns the inverse of a symmetric matrix using LU factorization using lapack or scalapack 
!>  o a real symmetric matrix
!>  o matrix_struct structure of a
!>  o do_scalapack flag whether lapack or scalapack is called
! **************************************************************************************************
   subroutine invert_real_matrix(matrix_a, matrix_struct, do_scalapack)
      real(kind=8), dimension(:,:), intent(inout)        :: matrix_a
      type(matrix_structure_type), intent(in)            :: matrix_struct
      logical, intent(in)                                :: do_scalapack

      if(do_scalapack) then
        call invmat_real_scalapack(matrix_a, matrix_struct)
      else
        call invmat_real_lapack(matrix_a)
      endif
   
   end subroutine invert_real_matrix

! **************************************************************************************************
!> \brief returns the inverse of a real matrix using LU factorization (scalapack)
!>  o a complex symmetric matrix
!>  o matrix_struct structure of a
!>  o do_scalapack flag whether lapack or scalapack is called
! **************************************************************************************************
   subroutine invert_cmplx_matrix(matrix_a, matrix_struct, do_scalapack)
      complex(kind=8), dimension(:,:), &
        intent(inout)                             :: matrix_a
      type(matrix_structure_type), intent(in)     :: matrix_struct
      logical, intent(in)                         :: do_scalapack
 
      if(do_scalapack) then
        call invmat_cmplx_scalapack(matrix_a, matrix_struct) 
      else
        call invmat_cmplx_lapack(matrix_a)
      endif

   end subroutine invert_cmplx_matrix

! **************************************************************************************************
!> \brief returns the inverse of a symmetric matrix using LU factorization (lapack)
!> o a real symmetric matrix 
! **************************************************************************************************
   subroutine invmat_real_lapack(a)
      real(kind=8), dimension(:, :), intent(inout)       :: a


      integer                                            :: lwork, n, info
      integer, dimension(:),  allocatable                :: ipiv
      real(kind=8), dimension(:), allocatable            :: work

      n = size(a, 1)
      lwork = 20*n
      allocate (ipiv(n))
      allocate (work(lwork))
      ipiv = 0
      work = 0.0d0
      info = 0
      call dgetrf(n, n, a, n, ipiv, info)
      if (info /=0 ) then
        call aims_stop('LU decomposition has failed')
      endif
      if (info == 0) THEN
         call dgetri(n, a, n, ipiv, work, lwork, info)
      end if
      if (info /=0 ) then
        call aims_stop('LU inversion has failed')
      endif
      deallocate (ipiv, work)
   end subroutine invmat_real_lapack

! **************************************************************************************************
!> \brief returns the inverse of a complex matrix using LU factorization (lapack)
!>  o a complex symmetric matrix
! **************************************************************************************************
   subroutine invmat_cmplx_lapack(a)
      complex(kind=8), dimension(:,:), intent(inout)     :: a


      integer                                            :: lwork, n, info
      integer, dimension(:), allocatable                 :: ipiv
      complex(kind=8), dimension(:), allocatable         :: work

      n = size(a, 1)
      lwork = 20*n
      allocate (ipiv(n))
      allocate (work(lwork))
      ipiv = 0
      work = (0.0d0, 0.0d0)
      info = 0
      call zgetrf(n, n, a, n, ipiv, info)
      if (info /=0 ) then
        call aims_stop('LU decomposition has failed')
      endif
      if (info == 0) THEN
         call zgetri(n, a, n, ipiv, work, lwork, info)
      end if
      if (info /=0 ) then
        call aims_stop('LU inversion has failed')
      endif
      deallocate (ipiv, work)
   end subroutine invmat_cmplx_lapack

! **************************************************************************************************
!> \brief returns the inverse of a real matrix using LU factorization (scalapack)
!>  o a real symmetric matrix
!>  o matrix_struct structure of a
! **************************************************************************************************
   subroutine invmat_real_scalapack(matrix_a, matrix_struct)
      real(kind=8), dimension(:,:), intent(inout)        :: matrix_a
      type(matrix_structure_type), intent(in)            :: matrix_struct
 
      integer                                            :: liwork, lwork
      integer                                            :: n, info
      integer, dimension(9)                              :: desca
      integer, dimension(:), allocatable                 :: ipivot 
      integer, dimension(:), allocatable                 :: iwork
      real(kind=8), dimension(:), allocatable            :: work
 
      n = matrix_struct%nrow_global 
      desca(:) = matrix_struct%descriptor
 
      allocate(ipivot(n+matrix_struct%nrow_block))
      ipivot(:) = 0

      call pdgetrf(n, n, matrix_a(1, 1), 1, 1, desca, ipivot, info)
      
      if (info /=0 ) then
        call aims_stop('LU decomposition has failed')
      endif
      !*** do work size query
      allocate(work(1),iwork(1))
      call pdgetri(n, matrix_a, 1, 1, desca, ipivot, work, -1, iwork, -1, info)
      lwork = int(work(1))
      liwork = int(iwork(1))
      deallocate(work, iwork)

      !*** do actual inversion
      allocate(work(lwork),iwork(liwork))
      call pdgetri(n, matrix_a, 1, 1, desca, ipivot, work, lwork, iwork, liwork, info)

      if (info /=0 ) then
        call aims_stop('LU inversion has failed')
      endif
     
     deallocate(work,iwork,ipivot)
 
   end subroutine invmat_real_scalapack

! **************************************************************************************************
!> \brief returns the inverse of a real matrix using LU factorization (scalapack)
!>  o a complex symmetric matrix
!>  o matrix_struct structure of a
! **************************************************************************************************
   subroutine invmat_cmplx_scalapack(matrix_a, matrix_struct)
      complex(kind=8), dimension(:,:), &
        intent(inout)                             :: matrix_a
      type(matrix_structure_type), intent(in)     :: matrix_struct
 
      integer                                     :: liwork, lwork
      integer                                     :: n, info
      integer, dimension(9)                       :: desca
      integer, dimension(:), allocatable          :: ipivot 
      integer, dimension(:), allocatable          :: iwork
      complex(kind=8), dimension(:), allocatable  :: work
 
      n = matrix_struct%nrow_global 
      desca(:) = matrix_struct%descriptor
 
      allocate(ipivot(n+matrix_struct%nrow_block))
      ipivot(:) = 0

      call pzgetrf(n, n, matrix_a(1, 1), 1, 1, desca, ipivot, info)
      
      if (info /=0 ) then
        call aims_stop('LU decomposition has failed')
      endif

      !*** do work size query
      allocate(work(1),iwork(1))
      call pzgetri(n, matrix_a, 1, 1, desca, ipivot, work, -1, iwork, -1, info)
      lwork = int(work(1))
      liwork = int(iwork(1))
      deallocate(work, iwork)

      !*** do actual inversion
      allocate(work(lwork),iwork(liwork))
      call pzgetri(n, matrix_a, 1, 1, desca, ipivot, work, lwork, iwork, liwork, info)

      if (info /=0 ) then
        call aims_stop('LU inversion has failed')
      endif
     
     deallocate(work,iwork,ipivot)
 
   end subroutine invmat_cmplx_scalapack

! **************************************************************************************************
!> \brief set diagonal and off-diagonal elements of a matrix
!>  o matrix_a -- real matrix
!>  o matrix_struct -- structure of matrix
!>  o off_val -- value of off-diagonal elements
!>  o diag_val --- value of diagonal elements
!>  o do_scalapack -- if we use scalapack or not 
! **************************************************************************************************
   subroutine set_matrix_elements(matrix_a, matrix_struct, off_val, diag_val, do_scalapack)
      real(kind=8), dimension(:,:), intent(inout)        :: matrix_a
      type(matrix_structure_type), intent(in)            :: matrix_struct
      real(kind=8), intent(in)                           :: off_val, diag_val
      logical, intent(in)                                :: do_scalapack
      
      integer                                            :: nrow, ncol
      integer, dimension(9)                              :: desca

      nrow = matrix_struct%nrow_global 
      ncol = matrix_struct%ncol_global 
      desca(:) = matrix_struct%descriptor

      if(do_scalapack) then
        call pdlaset('Full', nrow, ncol, off_val, diag_val, matrix_a, 1, 1, desca)
      else
        call dlaset('Full', nrow, ncol, off_val, diag_val, matrix_a, nrow)
      endif

   end subroutine set_matrix_elements
   
end module matrix_inversion
