!-------------------------------------------------------------------------------
!****s* FHI-aims/complex_lapack_solver_fast
!  NAME
!    complex_lapack_solver_fast
!  SYNOPSIS

subroutine rrs_pbc_complex_lapack_solver_fast(n_basis, n_states, overlap_matrix, ham,  &
                                      KS_eigenvalue, KS_eigenvector,info)

!  PURPOSE
!  Subroutine real_lapack_solver_fast solves the generalized eigenvalue problem
!  (Cholesky decomposition of overlap matrix, transformed EVP by normal diagonalisation,
!  backtransform) directly using standard lapack.
!
!  This version uses faster LAPACK routines than the standard one
!
!  USES
   use mpi_tasks
   use elpa1_2013
   use dimensions, ONLY: n_k_points
   use lapack_wrapper, ONLY: ovlp_complex

   implicit none

!  ARGUMENTS

   integer n_basis
   integer n_states
   complex*16, target :: overlap_matrix(n_basis,n_basis)
   complex*16 ham(n_basis,n_basis)

  ! JW: This subroutine is called without explicit interface (it is not a
  ! module subroutine) and there is no check that overlap_matrix actually has
  ! the target attribute.  Indeed, it does not, at least not explicitly.
  ! While I do not really expect trouble, this might be a possible
  ! explanation if there are problems with some compilers.

   real*8 KS_eigenvalue (n_states)
   complex*16 KS_eigenvector (n_basis, n_states)

!  INPUTS
!  o n_basis -- number of basis functions
!  o n_states -- number of states we want to calculate
!  o overlap_matrix -- overlap matrix, destroyed on exit !!!
!  o ham -- Hamiltonian matrix, destroyed on exit !!!
! 
!  OUTPUT
!  o KS_eigenvalue -- eigenvalues
!  o KS_eigenvector -- eigenvectors
! 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

   integer, parameter :: nblk = 128 ! Block size for matrix multiplies

   complex*16, parameter :: C_ZERO = (0.d0, 0.d0), C_ONE = (1.d0, 0.d0)

   complex*16, pointer :: tmp(:, :) ! this points to the working space

   complex*16, allocatable :: tau(:)
   real*8, allocatable :: d(:), e(:), evtmp(:, :)

   integer i, j, n, info, nwork

   character(*), parameter :: func = 'complex_lapack_solver_fast'

   ! If ovlp_complex is not allocated, then allocate it and calculate the
   ! inverse upper cholesky factor.
   ! If it is already allocated, nothing is to do since it contains
   ! the correct data

   if(.not. allocated(ovlp_complex)) then

      allocate(ovlp_complex(n_basis,n_basis))
      ovlp_complex(:,:) = 0

      ! Set upper part of ovlp_complex from overlap_matrix

      do i=1,n_basis
         do j=1,i
            ovlp_complex(j,i) = overlap_matrix(j,i)
         enddo
      enddo

      ! Cholesky factorization ovlp_complex = U**H * U

      call ZPOTRF('U', n_basis, ovlp_complex, n_basis, info)
      call check_info(info, 'ZPOTRF')

      ! Calculate the inverse of U and save it in ovlp_complex

      call ZTRTRI('U', 'N', n_basis, ovlp_complex, n_basis, info)
      call check_info(info, 'ZTRTRI')

   endif

   ! Since overlap_matrix is not needed any more in this routine
   ! (and it is also not needed after the call), we can use
   ! it as temprorary working space

   tmp => overlap_matrix

   allocate(d(n_basis))
   allocate(e(n_basis))
   allocate(tau(n_basis))

   ! Only the upper part of ham is set when this routine is called
   ! but the full matrix is needed for the multiplications below.

   do i=1,n_basis-1
      do j=i+1,n_basis
         ham(j,i) = conjg(ham(i,j))
      enddo
   enddo

   ! Transform problem to standard eigenvalue problem
   ! Compute: U**-H * Ham * U**-1

   ! Step 1: tmp = Ham * U**-1, only the upper triangle of blocks is needed.

   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call ZGEMM('N','N',n+nwork-1,nwork,n+nwork-1,C_ONE,ham(1,1),UBOUND(ham,1), &
                 ovlp_complex(1,n),UBOUND(ovlp_complex,1),C_ZERO,tmp(1,n),UBOUND(tmp,1))
   enddo


   ! Step 2: ham = U**-H * tmp, only the upper triangle of blocks is needed.

   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call ZGEMM('C','N',nwork,n_basis-n+1,n+nwork-1,C_ONE, &
                 ovlp_complex(1,n),UBOUND(ovlp_complex,1), &
                 tmp(1,n),UBOUND(tmp,1),C_ZERO,ham(n,n),UBOUND(ham,1))
   enddo

   ! Transform ham to a tridiagonal matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!

   call ZHETRD('U', n_basis, ham, UBOUND(ham,1), d, e, tau, tmp, size(tmp), info)
   call check_info(info, 'ZHETRD')

   ! Calculate eigenvalues of tridiagonal matrix.
   ! We use solve_tridi (from elpa1) instead of DSTEDC (Lapack)
   ! since solve_tridi calculates only the eigenvectors needed.
   ! Please note that the full space for all eigenvectors must be provided
   ! and thus we don't use KS_eigenvector for the eigenvectors

   allocate(evtmp(n_basis, n_basis))
   call solve_tridi_2013(n_basis, n_states, d, e, evtmp, UBOUND(evtmp,1), 64,&
        mpi_comm_self, mpi_comm_self)

   ! Store eigenvalues/eigenvectors

   KS_eigenvalue(1:n_states) = d(1:n_states)
   KS_eigenvector(1:n_basis,1:n_states) = evtmp(1:n_basis,1:n_states)
   deallocate(evtmp)

   ! Backtransform eigenvectors to eigenvectors of full matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!

   call ZUNMTR('L', 'U', 'N', n_basis, n_states, ham, UBOUND(ham,1), tau, &
               KS_eigenvector, UBOUND(KS_eigenvector,1), tmp, size(tmp), info)
   call check_info(info, 'ZUNMTR')

   ! Backtransform eigenvectors to the original (generalized) problem

   call ZTRMM('L','U','N','N', n_basis, n_states, C_ONE, ovlp_complex, &
              UBOUND(ovlp_complex,1), KS_eigenvector, UBOUND(KS_eigenvector,1))

   deallocate(d, e, tau)

   ! If n_k_points > n_tasks we cannot save the factored overlap matrix
   ! since every processor has more than one k point
   if(n_k_points > n_tasks) deallocate(ovlp_complex)

contains

   subroutine check_info(info, name)
      integer info
      character*(*) name
      character*150 :: info_str

      if(info/=0) then
         write(info_str,*) name,' failed, info= ', info
         !call aims_stop(info_str, func)
      endif

   end subroutine check_info

 end subroutine rrs_pbc_complex_lapack_solver_fast
!******	
!-------------------------------------------------------------------------------
