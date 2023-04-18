!****s* FHI-aims/logarithmic_matrix_lapack
!  NAME
!   logarithmic_matrix_lapack
!  SYNOPSIS

      subroutine logarithmic_matrix_lapack &
      ( n_basbas, auxmat, safe_minimum, threshold, name  &
      )
      !  n_nonsingular, eigenvalues, auxmat_transform, &
      !  name &
      !)

!  PURPOSE
!  diagonalize the overlap or Coulomb interaction matrix within the
!  auxiliary basis using standard lapack routines.
!
!  USES
      use mpi_tasks
      use localorb_io, only : use_unit
      implicit none


!  imported variables

!     input

      integer n_basbas
      real*8 auxmat(n_basbas,n_basbas)
      real*8 safe_minimum
      real*8 threshold
      character*(*), intent(IN) :: name

!     output

!     n_nonsingular: Number of non-singular eigenvectors of the overlap matrix
!     eigenvalues: Eigenvalues of the overlap matrix
!     auxmat_transform: Non-singular eigenvectors of the ovlp matrix

      integer n_nonsingular
      real*8, dimension(n_basbas) :: eigenvalues
!      real*8 auxmat_transform_2(n_basbas,n_basbas)
!      real*8 temp_auxmat(n_basbas,n_basbas)

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
!  o  coulomb_transform -- Non-singular eigenvectors of the Coulomb matrix
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
      integer i_basbas

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

      real*8, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: iwork

      real*8, dimension(:,:), allocatable :: auxmat_transform

!     counters

      integer i_eval

      real*8, parameter :: my_thres = -1d10  ! ~ - huge(my_thres)

!  begin work

!  DB: excluding the case that we are dealing with a pair of ghost atoms
!      without any basis function
      if(n_basbas.eq.0) return


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
      range = 'V'
      i_lower = 1
      i_upper = n_basbas
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
      abs_tol = 2.0d0 * safe_minimum

!     allocate necessary space

!      if (.not.allocated(work)) then
        allocate(work(8*n_basbas))
!      end if

!      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basbas))
!      end if

!      if (.not.allocated(ifail)) then
        allocate(ifail(n_basbas))
!      end if

      lwork = 8*n_basbas

      allocate(auxmat_transform(n_basbas,n_basbas))


!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

!      temp_auxmat = 0.0d0

      !call dsyevx &
      !( jobz, range, uplo, n_basbas, auxmat, n_basbas, v_lower, &
      !  v_upper,i_lower, i_upper, abs_tol, n_found, eigenvalues, &
      !  auxmat_transform, n_basbas, work, lwork, iwork, ifail, info &
      !)

      call diagonalize_auxmat_lapack(n_basbas, auxmat, safe_minimum, &
      & my_thres, n_nonsingular, eigenvalues, auxmat_transform, '')

      write(use_unit,*) "Igor debug", n_nonsingular, eigenvalues

!  consistency checks

      !if (info.lt.0) then
      !  write(use_unit,*) "Eigenvalue solver dsyevx():"
      !  write(use_unit,*) "The ", -info, "th argument in dspevx() had an ", &
      !    "illegal value. Check."
      !else if (info.gt.0) then
      !  write(use_unit,*) "Eigenvalue solver dsyevx():"
      !  write(use_unit,*) info, " eigenvectors failed to converge."
      !  write(use_unit,*) "ifail: ", ifail
      !  stop
      !end if

      !if (n_found.lt.n_basbas .and. myid.eq.0 .and. name/='') then
      !  write(use_unit,'(2X,2A,I8,A)') trim(name), " matrix is singular:"
      !  write(use_unit,'(2X,A,I8,A,I8,A)') &
      !    "| Using ", n_found, &
      !    " out of a possible ", n_basbas, &
      !    " specified basis functions."
      !end if

      !! should be used later
      !n_nonsingular = n_found

      auxmat(:,:) = 0.0d0
      do i_basbas = 1, n_basbas
        auxmat(i_basbas,i_basbas) = log(eigenvalues(i_basbas))
      enddo

      call dgemm &
      ( 'N','N',n_basbas,n_basbas,n_basbas,1.0d0,auxmat_transform,&
        n_basbas, auxmat, n_basbas, 0.0d0, auxmat, n_basbas&
      )

      call power_genmat_lapack(n_basbas, auxmat_transform, -1, &
      &                        safe_minimum, threshold, '')

      call dgemm &
      ( 'N','N',n_basbas,n_basbas,n_basbas,1.0d0, auxmat,&
        n_basbas, auxmat_transform, n_basbas, 0.0d0, auxmat, n_basbas &
      )

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
       if (allocated(ifail)) then
         deallocate(ifail)
       endif
       if (allocated(auxmat_transform)) then
         deallocate(auxmat_transform)
       endif


!  that's all folks

      end subroutine logarithmic_matrix_lapack
!-----------------------------------------------------------------------
!******
