!****s* FHI-aims/diagonalize_auxmat_lapack
!  NAME
!   diagonalize_auxmat_lapack
!  SYNOPSIS

      subroutine diagonalize_auxmat_lapack_complex &
      ( n_basis, auxmat, safe_minimum, threshold, &
        n_nonsingular, eigenvalues, auxmat_transform, &
        name &
      )

!  PURPOSE
!  diagonalize the overlap or Coulomb interaction matrix within the
!  auxiliary basis using standard lapack routines.
!
!  USES
      use localorb_io, only: use_unit
      use mpi_tasks
!      use physics 
      implicit none


!  imported variables

!     input
!      n_basbas=n_basis
      integer, intent(in) :: n_basis
      complex*16, intent(IN) :: auxmat(n_basis,n_basis)
      real*8, intent(IN) :: safe_minimum
      real*8, intent(IN) :: threshold
      character*(*), intent(IN) :: name

!     output

!     n_nonsingular: Number of non-singular eigenvectors of the overlap matrix
!     eigenvalues: Eigenvalues of the overlap matrix
!     auxmat_transform: Non-singular eigenvectors of the ovlp matrix

      integer, intent(OUT) :: n_nonsingular
      real*8, dimension(n_basis), intent(out) :: eigenvalues
      complex*16, intent(OUT) :: auxmat_transform(n_basis,n_basis)

!  INPUTS
!  o  n_basbas -- dimension of the (global) auxiliary basis
!  o  auxmat -- the overlap or Coulomb matrix within the auxiliary basis to
!         be diagonalized
!  o  safe_minimum -- the absolute tolerance for the eigenvalue precision needed by
!         the lapack eigenvalue solver.
!  o  threshold -- the cutoff threshold for the eigenvalues of the concerned
!         overlap or Coulomb matrix, only those above this threshold will be included.
!  o  name -- name of matrix (for output only); Choose '' for no output.
!
!  OUTPUTS
!  o  n_nonsingular -- real number,  number of nonsingular eigenvales larger than
!          a certain (positive) cutoff threshold.
!  o  eigenvalues -- Eigenvalues of the Coulomb matrix
!  o  auxmat_transform -- Non-singular eigenvectors of the Coulomb matrix
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

!     dspevx input

!     jobz     : calculate only eigenvalues, or eigenvalues and eigenvectors?
!     range    : specifies how we enter the range of eigenvalues / eigenvectors to be determined.
!     i_lower, i_upper : All eigenvalues between (in ascending order) number i_lower and i_upper
!             will be calculated.
!     uplo     : Determines the shape of the input (hermitian) matrices: Is the lower or the
!             upper triangular half given?
!     v_lower, v_upper : Dummy arguments which are unused for range = 'I'
!     abs_tol :  Absolute tolerance for eigenvalue precision in dsygvx

      character jobz
      character range
      integer i_lower, i_upper
      character uplo
      integer lwork
      real*8 v_lower, v_upper
      real*8 abs_tol

      character*100 :: info_str

!     dspevx output

!     n_found: total number of (auxiliary matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer, dimension(:), allocatable :: ifail
      integer info

!     dspevx work arrays

      complex*16, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: iwork
      real*8, dimension(:), allocatable :: rwork

!     counters

      integer i_eval

!  begin work

      if(myid.eq.0 .and. name /= '') then
         write(use_unit,'(2X,4A)') &
         & "Diagonalizing the ", trim(name), " matrix ", &
         & " and checking for singularities."
         write(use_unit,*)
      endif

!     initialize input for dspevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate all eigenvalues greater than a threshold.
!     i_lower, i_upper, v_upper are dummies for this case.
!  FIXME: The value for v_upper is not clear. The maximum eigenvalues of
!  FIXME: S_ij are in the range of 2, but I am not sure whether or not eigenvectors
!  FIXME: with eigenvalues cose to or above 2 are valid.
      range = 'A'
      i_lower = 1 
      i_upper = n_basis
      v_lower = threshold
      v_upper = 1.d5

!     in general, provide upper triangular input matrices
!     FIXME: Currently, we provide every single matrix element (a factor of 2 too expensive)
      uplo = 'U'

!     Set absolute tolerance for eigenvalue precision
!     2*safe_minimum is the lapack-recommended most precise choise.
!     However, safe_minimum may not work for obscure reasons.
!     (e.g. dmach() gives zero if overoptimized by the compiler.)
!     Try a lowermost value of 1.d-12 instead.
!      abs_tol = 2.0d0 * safe_minimum
      abs_tol = 2.0d0 *  1.d-12 

!     allocate necessary space

!      if (.not.allocated(work)) then
        allocate(work(10*n_basis))
!      end if

!      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basis))
!      end if
        allocate(rwork(7*n_basis))

!      if (.not.allocated(ifail)) then
        allocate(ifail(n_basis))
!      end if

      lwork = 8*n_basis

!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

       call zheevx &
      ( jobz, range, uplo, n_basis, auxmat, n_basis, v_lower, &
        v_upper,i_lower, i_upper, abs_tol, n_found, eigenvalues, &
        auxmat_transform, n_basis, work, lwork,rwork, iwork, ifail, info)

!      call zhpevx &
!      ( jobz, range, uplo, n_basbas, auxmat, v_lower, &
!        v_upper,i_lower, i_upper, abs_tol, n_found, eigenvalues, &
!        auxmat_transform, n_basbas, work,rwork, iwork, ifail, info)

!          call zheevx &
!          ( jobz, range, uplo, n_basbas, auxmat, n_basbas, v_lower, &
!            '',i_lower, '', abs_tol, n_found, eigenvalues, &
!            auxmat_transform, n_basbas, work, lwork,rwork, iwork, ifail, info)

!  consistency checks

      if (info.lt.0) then
        write(use_unit,*) "Eigenvalue solver zheevx():"
        write(use_unit,*) "The ", -info, "th argument in zheevx() had an ", &
          "illegal value. Check."
      else if (info.gt.0) then
        write(use_unit,*) "Eigenvalue solver  zheevx():"
        write(use_unit,*) info, " eigenvectors failed to converge."
        write(use_unit,*) "ifail: ", ifail
        stop
      end if

      if (n_found.lt.n_basis .and. myid.eq.0 .and. name/='') then
        write(use_unit,'(2X,2A,I8,A)') trim(name), " matrix is singular:"
        write(use_unit,'(2X,A,I8,A,I8,A)') &
          "| Using ", n_found, &
          " out of a possible ", n_basis, &
          " specified basis functions."
      end if

      n_nonsingular = n_found

!test write out eigenvalues and eigenvectors
!      write(use_unit,*)
!      write(use_unit,*) " Auxiliary eigenvalues:"
!      do i_eval = 1, n_basbas, 1
!        write(use_unit,*) i_eval, eigenvalues(i_eval)
!      enddo
!      write(use_unit,*)
!      write(use_unit,*) " Auxiliary eigenvectors:"
!      do i_eval = 1, n_basbas, 1
!        write(use_unit,*) (KS_eigenvector(i_basis, i_eval),
!     +              i_basis = 1, n_basbas, 1)
!      enddo
!      write(use_unit,*)
!test end
       if (allocated(work)) then
         deallocate(work)
       endif
       if (allocated(iwork)) then
         deallocate(iwork)
       endif
       if (allocated(rwork)) then
         deallocate(rwork)
       endif
       if (allocated(ifail)) then
         deallocate(ifail)
       endif


!  that's all folks
      end subroutine diagonalize_auxmat_lapack_complex
!-----------------------------------------------------------------------
!******
