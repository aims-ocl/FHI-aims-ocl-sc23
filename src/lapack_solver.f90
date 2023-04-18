!****s* FHI-aims/complex_lapack_solver
!  NAME
!    complex_lapack_solver
!  SYNOPSIS

subroutine complex_lapack_solver &
     (n_basis, n_states, overlap_matrix, hamiltonian, safe_min, &
     dim_work, KS_eigenvalue, KS_eigenvector)

!  PURPOSE
!  Subroutine complex_lapack_solver solves the generalized eigenvalue problem
!
!  This is a driver routine which calls either complex_lapack_solver_old
!  or complex_lapack_solver_fast
!
!  USES

  use mpi_tasks
  use runtime_choices
  use localorb_io,only:use_unit
  use aims_memory_tracking, only : aims_allocate, aims_deallocate

  implicit none

!  ARGUMENTS

  integer n_basis
  complex*16 overlap_matrix(n_basis, n_basis)
  complex*16 hamiltonian(n_basis, n_basis)
  integer n_states
  real*8 safe_min
  integer dim_work
  real*8 KS_eigenvalue (n_states)
  complex*16 KS_eigenvector (n_basis, n_states)

! INPUTS
! o n_basis -- number of basis functions
! o n_states -- number of eigenstates
! o overlap_matrix -- overlap matrix
! o hamiltonian -- hamiltonian matrix
! o safe_minimum --  tolerance for eigenvalue precision
! o dim_work -- dimension of work memory table
!  
!  OUTPUT
! o KS_eigenvalue -- eigenvalues
! o KS_eigenvector -- eigenvectors
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


    if(.not. use_lapack_fast) then
    ! VB: This output would be written for every k-point ... no good.
    !   if(myid==0) write(use_unit,'(2X,A)') &
    !      'Solving complex generalised eigenvalue problem by standard LAPACK (legacy)'
       call complex_lapack_solver_old &
          (n_basis, n_states, overlap_matrix, hamiltonian, safe_min, &
           dim_work, KS_eigenvalue, KS_eigenvector)
    else
    ! VB: This output would be written for every k-point ... no good.
    !   if(myid==0) write(use_unit,'(2X,A)') &
    !      'Solving complex generalised eigenvalue problem by standard LAPACK (fast)'
       call complex_lapack_solver_fast &
          (n_basis, n_states, overlap_matrix, hamiltonian, &
           KS_eigenvalue, KS_eigenvector)
    endif

    end subroutine complex_lapack_solver

!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/complex_lapack_solver_old
!  NAME
!    complex_lapack_solver_old
!  SYNOPSIS

subroutine complex_lapack_solver_old &
     (n_basis, n_states, overlap_matrix, hamiltonian, safe_minimum, &
     dim_work, KS_eigenvalue, KS_eigenvector)

!  PURPOSE
!  Subroutine complex_lapack_solver_old solves the generalized eigenvalue problem
!  (Cholesky decomposition of overlap matrix, transformed EVP by normal diagonalisation,
!  backtransform) directly using standard lapack.
!
!  This subroutine applies strictly only if the structure in question has no inversion
!  symmetry, i.e. if the Hamiltonian is actually complex. Otherwise, the "real"
!  lapack solver would suffice.
!
!  USES

   use localorb_io,only:use_unit
   use aims_memory_tracking, only : aims_allocate, aims_deallocate
   implicit none

!  ARGUMENTS

  integer n_basis
  complex*16 overlap_matrix(n_basis, n_basis)
  complex*16 hamiltonian(n_basis, n_basis)
  integer n_states
  real*8 safe_minimum
  integer dim_work
  real*8 KS_eigenvalue (n_states)
  complex*16 KS_eigenvector (n_basis, n_states)

! INPUTS
! o n_basis -- number of basis functions
! o n_states -- number of eigenstates
! o overlap_matrix -- overlap matrix
! o hamiltonian -- hamiltonian matrix
! o safe_minimum --  tolerance for eigenvalue precision
! o dim_work -- dimension of work memory table
!  
!  OUTPUT
! o KS_eigenvalue -- eigenvalues
! o KS_eigenvector -- eigenvectors
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

!  i_type   : Type of generalized EVP to be handled by zhegvx (lapack subroutine)
!             We want i_type=1, i.e. A*x = (lambda)*B*x .
!  jobz     : calculate only eigenvalues, or eigenvalues and eigenvectors?
!  range    : specifies how we enter the range of eigenvalues / eigenvectors to be determined.
!  i_lower, i_upper : All eigenvalues between (in ascending order) number i_lower and i_upper
!             will be calculated.
!  uplo     : Determines the shape of the input (hermitian) matrices: Is the lower or the
!             upper triangular half given?
!  matrix_dim : Passes actual matrix dimensions of overlap_matrix, hamiltonian to zhegvx().
!             In this way, we can use hamiltonian and overlap_matrix directly, instead of
!             copying them onto smaller working arrays first.
!  v_lower, v_upper : Dummy arguments which are unused for range = 'I'
!  abs_tol :  Absolute tolerance for eigenvalue precision in zhegvx
!  KS_dim :   passes actual dimension of KS_eigenvector to zhegvx
!  l_work :    dimension for workspace work

!     zhegvx input

      integer i_type
      character jobz
      character range
      integer i_lower, i_upper
      character uplo
      integer matrix_dim
      real*8  v_lower, v_upper
      real*8 abs_tol
      integer KS_dim
      integer l_work

!     zhegvx output

!  n_found :  Number of eigenvalues found (should be i_upper - i_lower +1)
!  output_eigenval : Array of output eigenvalues
!  i_fail : lists the eigenvectors which failed (or not) to converge.
!  info : Indicator of success (or not) of zhegvx

      integer n_found
      double precision output_eigenval (n_basis)
      integer i_fail(n_basis)
      integer info

!     zhegvx workspaces

!  work  :     complex workspace for zhegvx
!  rwork :     double precision workspace for zhegvx; dimension is 7*n_basis
!              (the actual dimension of the matrices overlap_matrix, hamiltonian)
!  iwork :     integer workspace for zhegvx; dimension is 5*N

!      complex*16 work(dim_work)
      complex*16,dimension(:),allocatable:: work
      double precision rwork(7*n_basis)
      integer iwork(5*n_basis)

!  counters

      integer i_basis
      integer i_state

!  begin work

!      write(use_unit,'(2X,A,A)')
!     +  "Solving complex generalised eigenvalue problem ",
!     +  "by standard LAPACK."

!  initialize input needed by zhegvx

!     Solve generalized EVP of type A*x = (lambda)*B*x .
      i_type = 1

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     zhgevx will calculate eigenvalues between i_lower and i_upper (in ascending order)
!     v_lower, v_upper are dummies for this case.
      range = 'I'
!test
!      range = 'A'
!test end
      i_lower = 1
      i_upper = n_states
      v_lower = 0.d0
      v_upper = 0.d0
      n_found = n_states

!     in general, provide upper triangular input matrices
!     FIXME: Currently, we provide every single matrix element (a factor of 2 too expensive)
      uplo = 'U'

!     store actual dimension of overlap_matrix, hamiltonian
!     FIXME : In f77, this allows to use overlap_matrix, hamiltonian directly, instead of
!     copying them to working arrays first. HOWEVER, this dim must be set equal to n_basis
!     once we actually use allocatable arrays.
      matrix_dim = n_basis

!     Set absolute tolerance for eigenvalue precision
!     2*safe_minimum is the lapack-recommended most precise choise.
!     However, safe_minimum may not work for obscure reasons. (dmach() gives zero.)
!     Try a lowermost value of 1.d-12 instead.
      abs_tol = 2.0d0 * safe_minimum

!     Store actual dimension of KS_eigenvector array
      KS_dim = n_basis

!  zhegvx is the lapack driver (expert mode) which solves the generalized
!  eigenvalue problem by standard diagonalisation. See the extensive comments in the
!  header of that subroutine for detailed explanations.

!  the first call only determines the optimum dimension l_work.
!     perform workspace size query first
      l_work = -1

!test
!      write(use_unit,*) "* Before Lapack 1."
!test end
      call aims_allocate( work, dim_work, "work" )

      call zhegvx &
       ( i_type, jobz, range, uplo, n_basis, hamiltonian, matrix_dim, &
         overlap_matrix, matrix_dim, v_lower, v_upper, i_lower, &
         i_upper, abs_tol, n_found, output_eigenval, KS_eigenvector, &
         KS_dim, work, l_work, rwork, iwork, i_fail, info &
       )

!test
!      write(use_unit,*)
!      write(use_unit,'(2X,A,A)') "Determination of workspace dimension for ",
!     +  "zhegvx() completed."
!      write(use_unit,'(2X,A,I4)') "Error indicator: info = ", info
!      write(use_unit,*)
!      write(use_unit,'(2X,A,A,2G12.4)') "Optimum work space dimension lwork, ",
!     +  "according to zhegvx:", work(1)
!      write(use_unit,'(2X,A,I6)') "Actual work space dimension (2*N-1) :",
!     +  (2*max_basis-1)
!test end

!     set work space dimension lwork correctly
      l_work = NINT(DBLE(work(1)))

      if(l_work > dim_work)then
         call aims_deallocate( work, "work" )
         call aims_allocate( work, l_work, "work" )
      end if

!test
!$$$      if (l_work.gt.dim_work) then
!$$$        write(use_unit,*) "Optimum dimension l_work =", l_work,
!$$$     +    ", determined by zhegvx(),"
!$$$        write(use_unit,*) "is larger than array dimension dim_work = ",
!$$$     +    dim_work, " determined from ilaenv()."
!$$$        write(use_unit,*) "Please check before proceeding."
!$$$       ! stop
!$$$      end if
!test end

!  zhegvx is the lapack driver (expert mode) which solves the generalized
!  eigenvalue problem by standard diagonalisation. See the extensive comments in the
!  header of that subroutine for detailed explanations.

!test
!      write(use_unit,*) "* Before Lapack 2."
!test end

!  Now really solve the eigenvalue problem.



      call zhegvx &
       ( i_type, jobz, range, uplo, n_basis, hamiltonian, matrix_dim, &
         overlap_matrix, matrix_dim, v_lower, v_upper, i_lower, &
         i_upper, abs_tol, n_found, output_eigenval, KS_eigenvector, &
         KS_dim, work, l_work, rwork, iwork, i_fail, info &
       )

      call aims_deallocate( work, "work" )

!test
!      write(use_unit,*) "* After Lapack 2."
!test end

!  consistency checks

      if (info.lt.0) then
        write(use_unit,*) "Generalized eigenvalue problem solver zhegvx():"
        write(use_unit,*) "The ", -info, "th argument in zhegvx() had an ", &
          "illegal value. Check."
        stop
      else if (info.gt.n_basis) then
        write(use_unit,*) "Generalized eigenvalue problem solver zhegvx():"
        write(use_unit,*) "The leading minor of order", (info-n_basis), &
          " of overlap_matrix() is not positive definite."
        stop
      else if (info.gt.0) then
        write(use_unit,*) "Generalized eigenvalue problem solver zhegvx():"
        write(use_unit,*) info, " eigenvectors failed to converge."
        write(use_unit,*) "i_fail: ", i_fail
        stop
      end if

      if (n_found.ne.n_states) then
        write(use_unit,*) "Inconsistent number of eigenvalues from lapack ", &
          "subroutine zhegvx:"
        write(use_unit,*) "n_found = ", n_found, ", should have been", &
        n_states, "."
        stop
      end if

!test write out eigenvalues and eigenvectors
!      write(use_unit,*)
!      write(use_unit,*) " Eigenvalues:"
!      do i_state = 1, n_states, 1
!        write(use_unit,*) i_state, output_eigenval(i_state)
!      enddo
!      write(use_unit,*)
!      write(use_unit,*) " Eigenvectors:"
!      do i_state = 1, n_states, 1
!        write(use_unit,*) (KS_eigenvector(i_basis, i_state),
!     +              i_basis = 1, n_basis, 1)
!      enddo
!      write(use_unit,*)
!test end

!  Save results in appropriate arrays

      do i_state = 1, n_states, 1
        KS_eigenvalue(i_state) = output_eigenval(i_state)
      enddo

      return
    end subroutine complex_lapack_solver_old
!******
!-------------------------------------------------------------------------------
!****s* FHI-aims/complex_lapack_solver_fast
!  NAME
!    complex_lapack_solver_fast
!  SYNOPSIS

subroutine complex_lapack_solver_fast(n_basis, n_states, ovlp, ham,  &
                                      KS_eigenvalue, KS_eigenvector)

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
   use aims_memory_tracking, only : aims_allocate, aims_deallocate

   implicit none

!  ARGUMENTS

   integer n_basis
   integer n_states
   complex*16 ovlp(n_basis,n_basis)
   complex*16 ham(n_basis,n_basis)

   real*8 KS_eigenvalue (n_states)
   complex*16 KS_eigenvector (n_basis, n_states)

!  INPUTS
!  o n_basis -- number of basis functions
!  o n_states -- number of states we want to calculate
!  o ovlp -- overlap matrix, destroyed on exit !!!
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

   complex*16, allocatable :: tau(:)
   real*8, allocatable :: d(:), e(:), evtmp(:, :)

   integer i, j, n, info, nwork

   character(*), parameter :: func = 'complex_lapack_solver_fast'

   ! If ovlp_complex is not allocated, then allocate it and calculate the
   ! inverse upper cholesky factor.
   ! If it is already allocated, nothing is to do since it contains
   ! the correct data

   if(.not. allocated(ovlp_complex)) then

      call aims_allocate ( ovlp_complex, n_basis,n_basis, "ovlp_complex" )
      ovlp_complex(:,:) = 0

      ! Set upper part of ovlp_complex from ovlp

      do i=1,n_basis
         do j=1,i
            ovlp_complex(j,i) = ovlp(j,i)
         enddo
      enddo

      ! Cholesky factorization ovlp_complex = U**H * U

      call ZPOTRF('U', n_basis, ovlp_complex, n_basis, info)
      call check_info(info, 'ZPOTRF')

      ! Calculate the inverse of U and save it in ovlp_complex

      call ZTRTRI('U', 'N', n_basis, ovlp_complex, n_basis, info)
      call check_info(info, 'ZTRTRI')

   endif

   ! From here on, *ovlp* is temporary storage!!!

   call aims_allocate( d, n_basis,     "d" )
   call aims_allocate( e, n_basis,     "e" )
   call aims_allocate( tau, n_basis, "tau" )

   ! Only the upper part of ham is set when this routine is called
   ! but the full matrix is needed for the multiplications below.

   do i=1,n_basis-1
      do j=i+1,n_basis
         ham(j,i) = conjg(ham(i,j))
      enddo
   enddo

   ! Transform problem to standard eigenvalue problem
   ! Compute: U**-H * Ham * U**-1

   ! Step 1: ovlp = Ham * U**-1, only the upper triangle of blocks is needed.

   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call ZGEMM('N','N',n+nwork-1,nwork,n+nwork-1,C_ONE,ham(1,1),UBOUND(ham,1), &
                 ovlp_complex(1,n),UBOUND(ovlp_complex,1),C_ZERO,ovlp(1,n),UBOUND(ovlp,1))
   enddo


   ! Step 2: ham = U**-H * ovlp, only the upper triangle of blocks is needed.

   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call ZGEMM('C','N',nwork,n_basis-n+1,n+nwork-1,C_ONE, &
                 ovlp_complex(1,n),UBOUND(ovlp_complex,1), &
                 ovlp(1,n),UBOUND(ovlp,1),C_ZERO,ham(n,n),UBOUND(ham,1))
   enddo

   ! Transform ham to a tridiagonal matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!

   call ZHETRD('U', n_basis, ham, UBOUND(ham,1), d, e, tau, ovlp, size(ovlp), info)
   call check_info(info, 'ZHETRD')

   ! Calculate eigenvalues of tridiagonal matrix.
   ! We use solve_tridi (from elpa) instead of DSTEDC (Lapack)
   ! since solve_tridi calculates only the eigenvectors needed.
   ! Please note that the full space for all eigenvectors must be provided
   ! and thus we don't use KS_eigenvector for the eigenvectors

   call aims_allocate( evtmp, n_basis, n_basis, "evtmp" )
   call solve_tridi_2013(n_basis, n_states, d, e, evtmp, UBOUND(evtmp,1), 64,&
        mpi_comm_self, mpi_comm_self)

   ! Store eigenvalues/eigenvectors

   KS_eigenvalue(1:n_states) = d(1:n_states)
   KS_eigenvector(1:n_basis,1:n_states) = evtmp(1:n_basis,1:n_states)
   call aims_deallocate( evtmp,                 "evtmp" )

   ! Backtransform eigenvectors to eigenvectors of full matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!

   call ZUNMTR('L', 'U', 'N', n_basis, n_states, ham, UBOUND(ham,1), tau, &
               KS_eigenvector, UBOUND(KS_eigenvector,1), ovlp, size(ovlp), info)
   call check_info(info, 'ZUNMTR')

   ! Backtransform eigenvectors to the original (generalized) problem

   call ZTRMM('L','U','N','N', n_basis, n_states, C_ONE, ovlp_complex, &
              UBOUND(ovlp_complex,1), KS_eigenvector, UBOUND(KS_eigenvector,1))

   call aims_deallocate( d,     "d" )
   call aims_deallocate( e,     "e" )
   call aims_deallocate( tau, "tau" )

   ! If n_k_points > n_tasks we cannot save the factored overlap matrix
   ! since every processor has more than one k point
   if(n_k_points > n_tasks) call aims_deallocate( ovlp_complex, "ovlp_complex" )

contains

   subroutine check_info(info, name)
      integer info
      character*(*) name
      character*150 :: info_str

      if(info/=0) then
         write(info_str,*) name,' failed, info= ', info
         call aims_stop(info_str, func)
      endif

   end subroutine check_info

 end subroutine complex_lapack_solver_fast
!******	
!-------------------------------------------------------------------------------
