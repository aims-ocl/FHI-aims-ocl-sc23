!****s* FHI-aims/power_auxmat_complex
!  NAME
!   power_auxmat_complex
!  SYNOPSIS

      subroutine power_auxmat_complex(n_basbas,power,auxmat,name)

!  PURPOSE
!  power the (screened) Coulomb interaction matrix within the
!  auxiliary basis using standard lapack routines.
!
!  USES
      use mpi_tasks
      use runtime_choices
      use localorb_io, only : use_unit
      implicit none


!  imported variables

!     input

      integer n_basbas
      real*8  power
      character*(*), intent(IN) :: name


!     input/output
      complex*16 auxmat(n_basbas,n_basbas)



!  INPUTS
!  o  n_basbas -- dimension of the (global) auxiliary basis
!  o  name --  name of matrix (for output only); Choose '' for no output.
!  o  power --  tell what power (.e.g., square root) of the input matrix should be computed
!  o  auxmat -- matrix of dimension (n_basbas,n_basbas)
!
!  OUTPUTS
!  o  auxmat -- the inverse of the input matrix
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

!   n_nonsingular    : number of eigenvalues larger than a chosen threshold
!   eigenvalues      : eigenvalues of the input matrix
!   auxmat_transform : eigenvectors of the input matrix

      integer n_nonsingular
      real*8  threshold
      real*8, dimension(:), allocatable :: eigenvalues
      complex*16, dimension(:,:), allocatable :: auxmat_transform
      complex*16, dimension(:), allocatable :: auxmat_input

!      real*8, dimension(n_basbas) :: eigenvalues
!      complex*16 auxmat_transform(n_basbas,n_basbas)

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
      real*8, dimension(:), allocatable :: rwork
      integer, dimension(:), allocatable :: iwork

!     counters

      integer i_eval
      integer i_prodbas_1, i_prodbas_2
      integer i_index

!  begin work

!      if(myid.eq.0 .and. name /= '') then
!         write(use_unit,'(2X,4A)') &
!         & "Diagonalizing the ", trim(name), " matrix ", &
!         & " and checking for singularities."
!         write(use_unit,*)
!      endif

      allocate(auxmat_input(n_basbas*(n_basbas+1)/2))
      allocate(eigenvalues(n_basbas))
      allocate(auxmat_transform(n_basbas,n_basbas))

      if(power .lt. 0) then
        threshold = prodbas_threshold
      else
        threshold = 0.d0
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

      i_index = 0
      do i_prodbas_2=1, n_basbas, 1
         do i_prodbas_1 =1, i_prodbas_2, 1
            i_index=i_index+1
            auxmat_input(i_index) = auxmat(i_prodbas_1,i_prodbas_2)
         enddo
      enddo

!     Set absolute tolerance for eigenvalue precision
!     2*safe_minimum is the lapack-recommended most precise choise.
!     However, safe_minimum may not work for obscure reasons.
!     (e.g. dmach() gives zero if overoptimized by the compiler.)
!     Try a lowermost value of 1.d-12 instead.
      abs_tol = 2.0d0 * safe_minimum

!     allocate necessary space

      allocate(work(2*n_basbas))
      allocate(rwork(7*n_basbas))
      allocate(iwork(5*n_basbas))
      allocate(ifail(n_basbas))


!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

     call zhpevx &
      ( jobz, range, uplo, n_basbas, auxmat_input, v_lower, v_upper, &
        i_lower, i_upper, abs_tol, n_found, eigenvalues, &
        auxmat_transform, n_basbas, work, rwork, iwork, ifail, info &
      )

!  consistency checks

      if (info.lt.0) then
        write(use_unit,*) "Eigenvalue solver zhpevx():"
        write(use_unit,*) "The ", -info, "th argument in zhpevx() had an ", &
          "illegal value. Check."
      else if (info.gt.0) then
        write(use_unit,*) "Eigenvalue solver zhpevx():"
        write(use_unit,*) info, " eigenvectors failed to converge."
        write(use_unit,*) "ifail: ", ifail
        stop
      end if

      if (n_found.lt.n_basbas .and. myid.eq.0 .and. name/='') then
        write(use_unit,'(2X,2A,I8,A)') trim(name), " matrix is singular:"
        write(use_unit,'(2X,A,I8,A,I8,A)') &
         "| Using ", n_found, &
         " out of a possible ", n_basbas, &
         " specified basis functions."
      end if

      n_nonsingular = n_found

!      write(use_unit,*) n_basbas, n_nonsingular
!      write(use_unit,*) eigenvalues(:)
      do i_eval = 1, n_nonsingular, 1
         auxmat_transform(:,i_eval) =  &
         auxmat_transform(:,i_eval)*(sqrt(eigenvalues(i_eval)))**power
      enddo

      call zgemm('N', 'C', n_basbas, n_basbas, n_nonsingular, &
         (1.0d0,0.d0), auxmat_transform(:,1:n_nonsingular), n_basbas, &
         auxmat_transform(:,1:n_nonsingular), n_basbas, &
         (0.d0,0.d0), auxmat, n_basbas)

        
      if (allocated(auxmat_input)) then
         deallocate(auxmat_input)
      endif
      if (allocated(eigenvalues)) then
         deallocate(eigenvalues)
      endif
      if (allocated(auxmat_transform)) then
         deallocate(auxmat_transform)
      endif
      if (allocated(work)) then
         deallocate(work)
      endif
      if (allocated(rwork)) then
         deallocate(rwork)
      endif
      if (allocated(iwork)) then
         deallocate(iwork)
      endif
      if (allocated(ifail)) then
         deallocate(ifail)
      endif
      return

!  that's all folks

      end subroutine power_auxmat_complex
!-----------------------------------------------------------------------
!******
