!****h* FHI-aims/lapack_wrapper
!  NAME
!   lapack_wrapper
!  SYNOPSIS

module lapack_wrapper

!  PURPOSE
!  Module lapack_wrapper provides wrappers around the high-level Lapack calls
!  which the program needs. The purpose is to allocate / deallocate workspaces
!  efficiently for use by different lapack routines.
!
!  USES

  use mpi_tasks
  use localorb_io
  use dimensions
  use runtime_choices

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

  !  dspvex, dspgvx, dgesvd, (and zhpevx) work arrays

  integer, dimension(:), allocatable, private :: ifail ! n_basis
  real*8, dimension(:), allocatable, private :: work   ! 8*n_basis or 7*n_basis
  integer, dimension(:), allocatable, private :: iwork ! 5*n_basis
  ! additional zhpevx work array
  complex*16, dimension(:), allocatable, private :: zwork   ! 2*n_basis

  !  dsysvx work arrays
  real*8, dimension(:,:), allocatable, private :: tmp_matrix
  integer, dimension(:),  allocatable, private :: i_piv
  real*8, dimension(:),   allocatable, private :: pulay_work

  !  Storage for the factored and inverted overlap matrix used in real_lapack_solver_fast
  !  This is saved from call to call if n_k_points <= n_tasks
  real*8, allocatable, private :: ovlp(:, :)
  complex*16, allocatable :: ovlp_complex(:, :) ! this is used in lapack_solver.f90

  integer,parameter:: BASIS_NON_SINGULAR = 0
  integer,parameter:: BASIS_SINGULAR = 1
  integer,parameter:: BASIS_SINGULAR_NOT_TESTED = -1

  logical :: suppress_further_ill_cond_output = .false. ! Added to prevent the situation where very dense k-grids
                                              ! would cause gigabytes of text to be output to stdout,
                                              ! making the calculation functionally useless

contains
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/allocate_work_space
!  NAME
!   allocate_work_space
!  SYNOPSIS

  subroutine allocate_work_space(n_basis)

! PURPOSE
! The subroutine allocates the module variables.
!
!  USES
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    integer:: n_basis

    if (.not.allocated(work)) then
       allocate(work(8*n_basis))
    end if

    if (.not.allocated(iwork)) then
       allocate(iwork(5*n_basis))
    end if

    if (.not.allocated(ifail)) then
       allocate(ifail(n_basis))
    end if

  end subroutine allocate_work_space

!******
!----------------------------------------------------------------------------
!****s* lapack_wrapper/am_i_output_thread
!  NAME
!    am_i_output_thread
!  SYNOPSIS

logical function am_i_output_thread(i_k_point) result(output_thread)

  !  PURPOSE
  !
  !    This is an awkward construction for the output. We want the output only
  !    once, but from a task that is actually used (if fewer k-points that MPI
  !    tasks, not all tasks actually call improve_eigenfunctions).  So we
  !    either write on task 0 or on task 1, depending on the large distinction
  !    of cases in get_KS_orbitals_p0.f90. Please see in that subroutine for
  !    details. This is not the cleanest construct in the code, and also
  !    error-prone.
  !
  !  USES

  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: i_k_point

  !  INPUTS
  !    o i_k_point -- k-point id
  !  OUTPUTS
  !    o output_thread -- Am I the output thread?
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  character(*), parameter :: func = 'am_i_output_thread'

  output_thread = .false.

  if (((n_periodic.eq.0) .and. (packed_matrix_format.ne.PM_index)) .or. &
  &   (n_tasks.eq.1)) then
     if ((myid.eq.0) .and. (i_k_point.eq.1)) then
        output_thread = .true.
     end if
  else
     ! in periodic systems, task 1 is the first used for a k-point,
     ! if available! except in serial runs!
     if ((myid.eq.1).and.(i_k_point.eq.1)) then
        output_thread = .true.
     end if
  end if

end function am_i_output_thread
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_overlap
!  NAME
!   diagonalize_overlap
!  SYNOPSIS

      subroutine diagonalize_overlap &
      ( n_basis, overlap_matrix, safe_minimum, basis_threshold, &
        n_nonsingular, eigenvalues, overlap_transform, i_k_point &
      )

!  PURPOSE
!  Subroutine diagonalize_overlap is a wrapper around the expert Lapack
!  driver to determine the eigenvalues (and hence, singular values) of
!  the overlap matrix. The output is a matrix overlap_tramsform, which
!  contains only the non-singular eigenvectors of the overlap matrix.
!
!  USES
    implicit none
!  ARGUMENTS

    integer n_basis
    real*8 overlap_matrix(n_basis*(n_basis+1)/2)
    real*8 safe_minimum
    real*8 basis_threshold
    integer n_nonsingular
    real*8, dimension(n_basis) :: eigenvalues
    real*8 overlap_transform(n_basis,n_basis)

    integer :: i_k_point

! INPUTS
! o n_basis -- number of basis functions
! o overlap_matrix -- overlap matrix
! o safe_minimum -- absolute tolerance for eigenvalue precision
! o basis_threshold -- minimum accepted eigenvalue
! o i_k_point -- provided here ONLY for output purposes - is not used and SHOULD NOT be
!                used anywhere else in this routine
!
! OUTPUT
! o n_nonsingular -- Number of non-singular eigenvectors of the overlap matrix
! o eigenvalues -- Eigenvalues of the overlap matrix
! o overlap_transform -- Non-singular eigenvectors of the ovlp matrix
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
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
      real*8 v_lower, v_upper
      real*8 abs_tol

!     dspevx output

!     n_found: total number of (overlap matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer info

      logical :: output_thread
      character*150 :: info_str
      character(*), parameter :: func = 'diagonalize_overlap'

!  begin work

      output_thread = am_i_output_thread(i_k_point)

      if (output_thread .and. (output_priority .le. OL_norm)) then
          write(use_unit,'(2X,A)') "Transforming overlap matrix and checking for singularities."
      end if

!     initialize input for dspevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate all eigenvalues greater than a threshold.
!     i_lower, i_upper, v_upper are dummies for this case.
!  FIXME: The value for v_upper is not clear.
      range = 'V'
      i_lower = 1
      i_upper = n_basis
      v_lower = basis_threshold
      v_upper = 10.d6

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

      if (.not.allocated(work)) then
        allocate(work(8*n_basis))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basis))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_basis))
      end if

!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call dspevx &
      ( jobz, range, uplo, n_basis, overlap_matrix, v_lower, v_upper, &
        i_lower, i_upper, abs_tol, n_found, eigenvalues, &
        overlap_transform, n_basis, work, iwork, ifail, info &
      )

!  consistency checks

      if (info /= 0) then
         write(info_str, "('Error from dspevx, info =',I8,'; i_kpoint =',I8)")&
         & info, i_k_point
         call aims_stop(info_str, func)
      end if

      if (n_found.lt.n_basis) then

         if (n_periodic.gt.0) then
           write(use_unit,'(2X,A,I8,A)') "Warning! At k-point ", i_k_point, ":"
         end if

         if ( (n_periodic.gt.0) .or. (output_thread .and. output_priority.le.OL_norm) ) then
            write(use_unit,'(2X,A,I8,A)') "Overlap matrix is singular:"
            write(use_unit,'(2X,A,I8,A,I8,A)') &
                 "| Using ", n_found, &
                 " out of a possible ", n_basis, &
                 " specified basis functions."
            write(use_unit,'(2X,A,E13.6)') &
                 "| Lowest remaining eigenvalue: ", eigenvalues(1)
         end if

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if

         if(n_states > n_found) then
            n_states = n_found
         endif
      else if (n_found.eq.n_basis) then

         if (output_thread.and.(output_priority.le.OL_norm)) then
            write(use_unit,'(2X,A,I8,A)') "Overlap matrix is nonsingular:"
            if (n_k_points.gt.1) then
              write(use_unit,'(2X,A,I8,A,E13.6)') &
                 "| Lowest eigenvalue at k-point ", i_k_point, ": ", eigenvalues(1)
            else
              write(use_unit,'(2X,A,E13.6)') &
                 "| Lowest eigenvalue: ", eigenvalues(1)
            end if
         end if

         if (eigenvalues(1).le.1.d-5) then
           if (n_periodic.gt.0) then
             write(use_unit,'(1X,A,I8,A)') "* At k-point ", i_k_point, ":"
           end if
           if (.not.suppress_further_ill_cond_output.and. ((n_periodic.gt.0) .or. output_thread) ) then
             write(use_unit,'(1X,A)') "* Warning! Overlap matrix is near-singular!"
             write(use_unit,'(1X,A,E13.6)') "* Lowest eigenvalue: ", eigenvalues(1)
             write(use_unit,'(1X,A)') "* Consider using a larger value of basis_threshold to ensure numerical stability!"
             write(use_unit,'(1X,A,I8,A)') "* Further output regarding ill-conditioning on task ", myid, " will be suppressed."
             suppress_further_ill_cond_output = .true.
           end if

           if (.not.override_illconditioning) then
             call stop_illconditioning()
           end if

         end if

      end if

      n_nonsingular = n_found

      end subroutine diagonalize_overlap
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_overlap_complex
! NAME
!  diagonalize_overlap_complex
! SYNOPSIS

      subroutine diagonalize_overlap_complex &
      ( n_basis, overlap_matrix, safe_minimum, basis_threshold, &
        n_nonsingular, eigenvalues, overlap_transform, i_k_point &
      )

! PURPOSE
!  Subroutine diagonalize_overlap_complex is a wrapper around the expert Lapack
!  driver to determine the eigenvalues (and hence, singular values) of
!  the overlap matrix. The output is a matrix overlap_tramsform, which
!  contains only the non-singular eigenvectors of the overlap matrix.
!
! USES
      implicit none
! ARGUMENTS

      integer n_basis
      complex*16 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 safe_minimum
      real*8 basis_threshold
      integer i_k_point

      integer n_nonsingular
      real*8, dimension(n_basis) :: eigenvalues
      complex*16 overlap_transform(n_basis,n_basis)

! INPUTS
! o n_basis -- number of basis functions
! o overlap_matrix -- overlap matrix
! o safe_minimum -- absolute tolerance for eigenvalue precision
! o basis_threshold -- minimum accepted eigenvalue
! o i_k_point -- k-point index
!
! OUTPUT
! o n_nonsingular -- Number of non-singular eigenvectors of the overlap matrix
! o eigenvalues -- Eigenvalues of the overlap matrix
! o overlap_transform -- Non-singular eigenvectors of the ovlp matrix
! AUTHOR
!   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
! HISTORY
!   Release version, FHI-aims (2008).
! SOURCE

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
      real*8 v_lower, v_upper
      real*8 abs_tol

!     dspevx output

!     n_found: total number of (overlap matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer info

      logical :: output_thread
      character(*), parameter :: func = 'diagonalize_overlap_complex'

      integer i_state

!  begin work

      output_thread = am_i_output_thread(i_k_point)

      if (output_thread .and. (output_priority .le. OL_norm)) then
        write(use_unit,'(2X,A)') &
          "Transforming overlap matrix and checking for singularities."
      end if

!     initialize input for dspevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate all eigenvalues greater than a threshold.
!     i_lower, i_upper, v_upper are dummies for this case.
!  FIXME: The value for v_upper is not clear.
      range = 'V'
!      range = 'I'
      i_lower = 1
      i_upper = n_basis
      v_lower = basis_threshold
      v_upper = 10.d6

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
      if (.not.allocated(zwork)) then
        allocate(zwork(2*n_basis))
      end if

      if (.not.allocated(work)) then
        allocate(work(7*n_basis))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basis))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_basis))
      end if

!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call zhpevx &
      ( jobz, range, uplo, n_basis, overlap_matrix, v_lower, v_upper, &
        i_lower, i_upper, abs_tol, n_found, eigenvalues, &
        overlap_transform, n_basis, zwork, work, iwork, ifail, info &
      )

!     write(use_unit,*) "basis_threshold:", basis_threshold
!     do i_state = 1, n_basis, 1
!      write(use_unit,'(I4,2f16.8)') i_state, eigenvalues(i_state)
!     enddo
!  consistency checks

      if (info.lt.0) then

         write(use_unit,*) "Eigenvalue solver zhpevx():"
         write(use_unit,*) "The ", -info, "th argument in zhpgvx() had an ", &
              "illegal value for k-point ", i_k_point," . Check."
         call aims_stop('', func)

      else if (info.gt.0) then

         write(use_unit,*) "Eigenvalue solver zhpevx():"
         write(use_unit,*) info, " eigenvectors failed to converge"
         write(use_unit,*) " for k-point ", i_k_point
         write(use_unit,*) "ifail: ", ifail
         call aims_stop('', func)

      end if

      if (n_found.lt.n_basis) then

         if (.not.suppress_further_ill_cond_output) then
           write(use_unit,'(2X,A,A,I8)') "Overlap matrix is singular ", &
                "at k-point: ", i_k_point
           write(use_unit,'(2X,A,I8,A,I8,A)') &
                "| Using ", n_found, &
                " out of a possible ", n_basis, &
                " specified basis functions."
           write(use_unit,'(2X,A,E13.6)') &
                   "| Lowest eigenvalue: ", eigenvalues(1)
           write(use_unit,'(2X,A,I8,A)') "Further output regarding ill-conditioning on task ", myid, " will be suppressed."
           suppress_further_ill_cond_output = .true.
         end if

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if
      else

         if (output_thread .and. (output_priority .le. OL_norm)) then
            write(use_unit,'(2X,A,I8,A)') "Overlap matrix is nonsingular at first k-point:"
            write(use_unit,'(2X,A,E13.6)') &
                 "| Lowest eigenvalue: ", eigenvalues(1)
        end if

        if (eigenvalues(1).le.1d-5) then
          if (.not.suppress_further_ill_cond_output) then
            write(use_unit,'(1X,A,I8,A)') "* At k-point ", i_k_point, ":"
            write(use_unit,'(1X,A)') "* Warning! Overlap matrix is near-singular!"
            write(use_unit,'(1X,A,E13.6)') "* Lowest eigenvalue: ", eigenvalues(1)
            write(use_unit,'(1X,A)') "* Consider using a larger value of basis_threshold to ensure numerical stability!"
            write(use_unit,'(1X,A,I8,A)') "* Further output regarding ill-conditioning on task ", myid, " will be suppressed."
            suppress_further_ill_cond_output = .true.
          end if

          if (.not.override_illconditioning) then
            call stop_illconditioning()
          end if

        end if

      end if

      n_nonsingular = n_found

      end subroutine diagonalize_overlap_complex
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_overlap_complex_check
! NAME
!  diagonalize_overlap_complex_check
! SYNOPSIS

      subroutine diagonalize_overlap_complex_check &
      ( n_basis, overlap_matrix, safe_minimum, basis_threshold, &
        n_nonsingular, eigenvalues, overlap_transform, i_k_point &
      )

! PURPOSE
!  Subroutine diagonalize_overlap_complex_check is a wrapper around the expert Lapack
!  driver to determine the eigenvalues (and hence, singular values) of
!  the overlap matrix. The output is a matrix overlap_tramsform, which
!  contains only the non-singular eigenvectors of the overlap matrix.
!
! USES
      implicit none
! ARGUMENTS

      integer n_basis
      complex*16 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 safe_minimum
      real*8 basis_threshold
      integer i_k_point

      integer n_nonsingular
      real*8, dimension(n_basis) :: eigenvalues
      complex*16 overlap_transform(n_basis,n_basis)

! INPUTS
! o n_basis -- number of basis functions
! o overlap_matrix -- overlap matrix
! o safe_minimum -- absolute tolerance for eigenvalue precision
! o basis_threshold -- minimum accepted eigenvalue
! o i_k_point -- k-point index
!
! OUTPUT
! o n_nonsingular -- Number of non-singular eigenvectors of the overlap matrix
! o eigenvalues -- Eigenvalues of the overlap matrix
! o overlap_transform -- Non-singular eigenvectors of the ovlp matrix
! AUTHOR
!   FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
! HISTORY
!   Release version, FHI-aims (2008).
! SOURCE

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
      real*8 v_lower, v_upper
      real*8 abs_tol

!     dspevx output

!     n_found: total number of (overlap matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer info

      logical :: output_thread
      character(*), parameter :: func = 'diagonalize_overlap_complex_check'

!     counters

!      integer i_eval
!      integer i_basis
!      integer i_offset
      integer i_state

!  begin work

      output_thread = am_i_output_thread(i_k_point)

      if (output_thread .and. (output_priority .le. OL_norm)) then
        write(use_unit,'(2X,A)') &
          "Transforming overlap matrix and checking for singularities."
      end if

!     initialize input for dspevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate all eigenvalues greater than a threshold.
!     i_lower, i_upper, v_upper are dummies for this case.
!  FIXME: The value for v_upper is not clear.
!      range = 'V'
      range = 'I'
      i_lower = 1
      i_upper = n_basis
      v_lower = basis_threshold
      v_upper = 10.d6

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
      if (.not.allocated(zwork)) then
        allocate(zwork(2*n_basis))
      end if

      if (.not.allocated(work)) then
        allocate(work(7*n_basis))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basis))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_basis))
      end if

!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call zhpevx &
      ( jobz, range, uplo, n_basis, overlap_matrix, v_lower, v_upper, &
        i_lower, i_upper, abs_tol, n_found, eigenvalues, &
        overlap_transform, n_basis, zwork, work, iwork, ifail, info &
      )

!  consistency checks

      if (info.lt.0) then

         write(use_unit,*) "Eigenvalue solver zhpevx():"
         write(use_unit,*) "The ", -info, "th argument in zhpgvx() had an ", &
              "illegal value for k-point ", i_k_point," . Check."
         call aims_stop('', func)

      else if (info.gt.0) then

         write(use_unit,*) "Eigenvalue solver zhpevx():"
         write(use_unit,*) info, " eigenvectors failed to converge"
         write(use_unit,*) " for k-point ", i_k_point
         write(use_unit,*) "ifail: ", ifail
         call aims_stop('', func)

      end if

      if (n_found.lt.n_basis) then

         write(use_unit,'(2X,A,A,I8)') "For checking purpose only: overlap matrix is singular ", &
              "at k-point: ", i_k_point
         write(use_unit,'(2X,A,I8,A,I8,A)') &
              "| Using ", n_found, &
              " out of a possible ", n_basis, &
              " specified basis functions."
         write(use_unit,'(2X,A,E13.6)') &
                 "| Lowest eigenvalue: ", eigenvalues(1)

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if

        if(n_states > n_found) then
           n_states= n_found
        endif

      else

         if (output_thread .and. (output_priority .le. OL_norm)) then
            write(use_unit,'(2X,A,I8,A)') "Overlap matrix is nonsingular at first k-point:"
            write(use_unit,'(2X,A,E13.6)') &
                 "| Lowest eigenvalue: ", eigenvalues(1)
        end if

        if (eigenvalues(1).le.1d-5) then
          write(use_unit,'(1X,A,I8,A)') "* At k-point ", i_k_point, ":"
          write(use_unit,'(1X,A)') "* Warning! Overlap matrix is near-singular!"
          write(use_unit,'(1X,A,E13.6)') "* Lowest eigenvalue: ", eigenvalues(1)
          write(use_unit,'(1X,A)') "* Consider using a larger value of basis_threshold to ensure numerical stability!"

          if (.not.override_illconditioning) then
            call stop_illconditioning()
          end if

        end if

      end if

      n_nonsingular = n_found

      end subroutine diagonalize_overlap_complex_check
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/real_lapack_solver
!  NAME
!   real_lapack_solver
!  SYNOPSIS

      subroutine real_lapack_solver &
        (n_basis, n_states, overlap_matrix, hamiltonian, safe_minimum, &
         t_out,  KS_eigenvalue, KS_eigenvector)

!  PURPOSE
!  Subroutine real_lapack_solver solves the generalized eigenvalue problem
!
!  This is a driver routine which calls either real_lapack_solver_old
!  or real_lapack_solver_fast
!
!  USES

        implicit none
!  ARGUMENTS

      integer n_basis
      real*8 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 hamiltonian(n_basis*(n_basis+1)/2)
      integer n_states
      real*8 safe_minimum

      logical :: t_out

      real*8 KS_eigenvalue (n_states)
      real*8 KS_eigenvector (n_basis, n_states)

! INPUTS
! o n_basis -- number of basis functions
! o n_states -- number of states we want to calculate
! o overlap_matrix -- overlap matrix
! o hamiltonian -- Hamiltonian matrix
! o safe_minimum -- accuracy of the eigenvalues
! o t_out -- is the information printed out or not?
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

      character*100 :: info_str

    if(.not. use_lapack_fast) then
       if (t_out) then
          write (info_str,'(2X,2A)') &
          & 'Solving real generalised eigenvalue problem ', &
          & 'by standard LAPACK (legacy)'
          call localorb_info(info_str,use_unit,'(A)',OL_norm)
       end if
       call real_lapack_solver_old &
          (n_basis, n_states, overlap_matrix, hamiltonian, safe_minimum, &
           KS_eigenvalue, KS_eigenvector)
    else
       if (t_out) then
          write (info_str,'(2X,2A)') &
          & 'Solving real generalised eigenvalue problem ', &
          & 'by standard LAPACK (fast)'
          call localorb_info(info_str,use_unit,'(A)',OL_norm)
       end if
       call real_lapack_solver_fast &
          (n_basis, n_states, overlap_matrix, hamiltonian, &
           KS_eigenvalue, KS_eigenvector)
    endif

    end subroutine real_lapack_solver
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/real_lapack_solver_old
!  NAME
!   real_lapack_solver_old
!  SYNOPSIS

      subroutine real_lapack_solver_old &
        (n_basis, n_states, overlap_matrix, hamiltonian, safe_minimum, &
         KS_eigenvalue, KS_eigenvector)

!  PURPOSE
!  Subroutine real_lapack_solver_old solves the generalized eigenvalue problem
!  (Cholesky decomposition of overlap matrix, transformed EVP by normal diagonalisation,
!  backtransform) directly using standard lapack.
!
!  Note: The real version uses packed overlap/hamiltonian matrices, the complex version does not.
!
!  USES
        implicit none
!  ARGUMENTS

      integer n_basis
      real*8 overlap_matrix(n_basis*(n_basis+1)/2)
      real*8 hamiltonian(n_basis*(n_basis+1)/2)
      integer n_states
      real*8 safe_minimum

      real*8 KS_eigenvalue (n_states)
      real*8 KS_eigenvector (n_basis, n_states)

! INPUTS
! o n_basis -- number of basis functions
! o n_states -- number of states we want to calculate
! o overlap_matrix -- overlap matrix
! o hamiltonian -- Hamiltonian matrix
! o safe_minimum -- accuracy of the eigenvalues
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

!     dspgvx input

!  i_type   : Type of generalized EVP to be handled by dsygvx (lapack subroutine)
!             We want i_type=1, i.e. A*x = (lambda)*B*x .
!  jobz     : calculate only eigenvalues, or eigenvalues and eigenvectors?
!  range    : specifies how we enter the range of eigenvalues / eigenvectors to be determined.
!  i_lower, i_upper : All eigenvalues between (in ascending order) number i_lower and i_upper
!             will be calculated.
!  uplo     : Determines the shape of the input (hermitian) matrices: Is the lower or the
!             upper triangular half given?
!  v_lower, v_upper : Dummy arguments which are unused for range = 'I'
!  abs_tol :  Absolute tolerance for eigenvalue precision in dsygvx

      integer i_type
      character jobz
      character range
      integer i_lower, i_upper
      character uplo
      real*8  v_lower, v_upper
      real*8 abs_tol

!     dspgvx output

!  n_found :  Number of eigenvalues found (should be i_upper - i_lower +1)
!  output_eigenval : Array of output eigenvalues
!  info : Indicator of success (or not) of dsygvx

      integer n_found
      double precision output_eigenval (n_basis)
      integer info

      integer i_state

      character*100 :: info_str
!  begin work

!     allocate necessary space

      if (.not.allocated(work)) then
        allocate(work(8*n_basis))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_basis))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_basis))
      end if

!  initialize input needed by dsygvx

!     Solve generalized EVP of type A*x = (lambda)*B*x .
      i_type = 1

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate eigenvalues between i_lower and i_upper (in ascending order)
!     v_lower, v_upper are dummies for this case.
      range = 'I'

      i_lower = 1
      i_upper = n_states
      v_lower = 0.d0
      v_upper = 0.d0
      n_found = n_states

!     in general, provide upper triangular input matrices
!     FIXME: Currently, we provide every single matrix element (a factor of 2 too expensive)
      uplo = 'U'

!     Set absolute tolerance for eigenvalue precision
!     2*safe_minimum is the lapack-recommended most precise choise.
!     However, safe_minimum may not work for obscure reasons. (dmach() gives zero.)
!     Try a lowermost value of 1.d-12 instead.
      abs_tol = 2.0d0 * safe_minimum

!  dspgvx is the lapack driver (expert mode) which solves the generalized
!  eigenvalue problem by standard diagonalisation for packed overlap and
!  Hamiltonian matrices. See the extensive comments in the header of that
!  subroutine for detailed explanations.

      call dspgvx &
       ( i_type, jobz, range, uplo, n_basis, hamiltonian, &
         overlap_matrix, v_lower, v_upper, i_lower, &
         i_upper, abs_tol, n_found, output_eigenval, KS_eigenvector, &
         n_basis, work, iwork, ifail, info &
       )

!  consistency checks

      if (info.lt.0) then

         if (myid.eq.0) then
            write(use_unit,*) "* Generalized eigenvalue problem solver dspgvx():"
            write(use_unit,*) "* The ", -info, "th argument in dspgvx() had an ", &
                 "illegal value. Check."
         end if

         stop
      else if (info.gt.n_basis) then

         if (myid.eq.0) then
            write(use_unit,*) "* Generalized eigenvalue problem solver dspgvx():"
            write(use_unit,*) "* The leading minor of order", (info-n_basis), &
                 " of overlap_matrix() is not positive definite."
         end if

         stop
      else if (info.gt.0) then

         if (myid.eq.0) then
            write(use_unit,*) "* Generalized eigenvalue problem solver dspgvx():"
            write(use_unit,*) "*", info, " eigenvectors failed to converge."
            write(use_unit,*) "* ifail: ", ifail
         end if

         stop
      end if

      if (n_found.ne.n_states) then

         if (myid.eq.0) then
            write(use_unit,*) "* Inconsistent number of eigenvalues from lapack ", &
                 "subroutine dspgvx:"
            write(use_unit,*) "* n_found = ", n_found, ", should have been", &
                 n_states, "."
         end if

         stop
      end if

!  Save results in appropriate arrays

      do i_state = 1, n_states, 1
        KS_eigenvalue(i_state) = output_eigenval(i_state)
      enddo

      return
      end subroutine real_lapack_solver_old
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/real_lapack_solver_fast
!  NAME
!   real_lapack_solver_fast
!  SYNOPSIS

subroutine real_lapack_solver_fast(n_basis, n_states, overlap_matrix, hamiltonian, &
                                   KS_eigenvalue, KS_eigenvector)

!  PURPOSE
!  Subroutine real_lapack_solver_fast solves the generalized eigenvalue problem
!  (Cholesky decomposition of overlap matrix, transformed EVP by normal diagonalisation,
!  backtransform) directly using standard lapack.
!
!  This version unpacks the packed matrices and uses the unpacked matrices
!  for LAPACK calls since this is faster (although needing more memory)
!
!  USES
   use elpa1_2013
   use aims_memory_tracking, only : aims_allocate, aims_deallocate

   implicit none
   include "mpif.h"

!  ARGUMENTS

   integer n_basis
   integer n_states
   real*8 overlap_matrix(n_basis*(n_basis+1)/2)
   real*8 hamiltonian(n_basis*(n_basis+1)/2)

   real*8 KS_eigenvalue (n_states)
   real*8 KS_eigenvector (n_basis, n_states)

!  INPUTS
!  o n_basis -- number of basis functions
!  o n_states -- number of states we want to calculate
!  o overlap_matrix -- overlap matrix
!  o hamiltonian -- Hamiltonian matrix
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

   real*8, parameter :: R_ZERO = 0.d0, R_ONE = 1.d0

   real*8, allocatable :: ham(:, :), tmp(:, :)
   real*8, allocatable :: d(:), e(:), tau(:)
   real*8              :: b_time,e_time,t_time
   integer i, j, n, info, nwork

   ! If ovlp is not allocated, then allocate it and calculate the
   ! inverse upper cholesky factor.
   ! If it is already allocated, nothing is to do since it contains
   ! the correct data

   if(.not.allocated(ovlp)) then

      call aims_allocate( ovlp, n_basis,n_basis, "ovlp" )
      ovlp = 0

      ! Fill upper triangle with packed matrix
      n = 0
      do i=1,n_basis
         do j=1,i
            n = n+1
            ovlp(j,i) = overlap_matrix(n)
         enddo
      enddo

      ! Cholesky factorization ovlp = U**T * U
      call DPOTRF('U', n_basis, ovlp, n_basis, info)
      call check_info(info, 'DPOTRF')

      ! Calculate the inverse of U and save it in ovlp

      call DTRTRI('U', 'N', n_basis, ovlp, n_basis, info)
      call check_info(info, 'DTRTRI')
   endif

   ! Set ham (full matrix is needed for the multiplications below)

   call aims_allocate( ham, n_basis,n_basis, "ham" )
   call aims_allocate( tmp, n_basis,n_basis, "tmp" )
   call aims_allocate( d, n_basis,             "d" )
   call aims_allocate( e, n_basis,             "e" )
   call aims_allocate( tau, n_basis,         "tau" )

   n = 0
   do i=1,n_basis
      do j=1,i
         n = n+1
         ham(j,i)  = hamiltonian(n)
         ham(i,j)  = hamiltonian(n)
      enddo
   enddo

   ! Transform problem to standard eigenvalue problem
   ! Compute: U**-T * Ham * U**-1

   ! Step 1: tmp = Ham * U**-1, only the upper triangle of blocks is needed.
   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call DGEMM('N','N',n+nwork-1,nwork,n+nwork-1,R_ONE,ham(1,1),n_basis, &
                 ovlp(1,n),n_basis,R_ZERO,tmp(1,n),n_basis)
   enddo

   ! Step 2: ham = U**-T * tmp, only the upper triangle of blocks is needed.

   do n=1,n_basis,nblk
      nwork = nblk
      if(n+nwork-1>n_basis) nwork = n_basis-n+1
      call DGEMM('T','N',nwork,n_basis-n+1,n+nwork-1,R_ONE,ovlp(1,n),n_basis, &
                 tmp(1,n),n_basis,R_ZERO,ham(n,n),n_basis)
   enddo

   ! Transform ham to a tridiagonal matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!
   call DSYTRD('U', n_basis, ham, UBOUND(ham,1), d, e, tau, tmp, size(tmp), info)
   call check_info(info, 'DSYTRD')

   ! Calculate eigenvalues of tridiagonal matrix.
   ! We use solve_tridi (from elpa) instead of DSTEDC (Lapack)
   ! since solve_tridi calculates only the eigenvectors needed.
   ! Please note that the full space for all eigenvectors must be provided
   ! and thus we don't use KS_eigenvector for the eigenvectors

   call solve_tridi_2013(n_basis, n_states, d, e, tmp, UBOUND(tmp,1), 64,&
        mpi_comm_self, mpi_comm_self)

   ! Store eigenvalues/eigenvectors

   KS_eigenvalue(1:n_states) = d(1:n_states)
   KS_eigenvector(1:n_basis,1:n_states) = tmp(1:n_basis,1:n_states)

   ! Backtransform eigenvectors to eigenvectors of full matrix
   ! The provided workspace (n_basis*n_basis) should be enough for optimum work!
   call DORMTR('L', 'U', 'N', n_basis, n_states, ham, UBOUND(ham,1), tau, &
               KS_eigenvector, UBOUND(KS_eigenvector,1), tmp, size(tmp), info)
   call check_info(info, 'DORMTR')

   ! Backtransform eigenvectors to the original (generalized) problem

   call DTRMM('L','U','N','N', n_basis, n_states, R_ONE, ovlp, n_basis, KS_eigenvector, n_basis)

   ! Allocatable array that are tracked
   call aims_deallocate( ham, "ham" )
   call aims_deallocate( tmp, "tmp" )
   call aims_deallocate( d,     "d" )
   call aims_deallocate( e,     "e" )
   call aims_deallocate( tau, "tau" )

   ! If n_k_points > n_tasks we cannot save the factored overlap matrix
   ! since every processor has more than one k point
   if(n_k_points > n_tasks) call aims_deallocate( ovlp, "ovlp" )

contains

   subroutine check_info(info, name)
      integer info
      character*(*) name

      if(info/=0) then
         write(use_unit,*) name,' failed in real_lapack_solver_fast, info= ',info
         call aims_stop
      endif

   end subroutine check_info

end subroutine real_lapack_solver_fast
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/clear_lapack_overlap_matrix
!  NAME
!   clear_lapack_overlap_matrix
!  SYNOPSIS

subroutine clear_lapack_overlap_matrix

!  PURPOSE
!  Deallocates a saved overlap matrix after a relaxation step
   use aims_memory_tracking, only : aims_deallocate
   implicit none

   if(allocated(ovlp))         call aims_deallocate( ovlp,                 "ovlp" )
   if(allocated(ovlp_complex)) call aims_deallocate( ovlp_complex, "ovlp_complex" )

end subroutine clear_lapack_overlap_matrix
!
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_hamiltonian
!  NAME
!   diagonalize_hamiltonian
!  SYNOPSIS

     subroutine diagonalize_hamiltonian &
      ( n_nonsingular, n_states, trafo_hamiltonian, &
        safe_minimum, KS_eigenvalue, aux_eigenvector &
      )

!  PURPOSE
!  Subroutine diagonalize_hamiltonian is a wrapper around the expert Lapack
!  driver to determine the eigenvalues of the transformed (non-singular) Hamiltonian.
!
!  USES
      implicit none
!  ARGUMENTS


!  imported variables

!     input

      integer n_nonsingular
      integer n_states
      real*8 trafo_hamiltonian(n_nonsingular*(n_nonsingular+1)/2)
      real*8 safe_minimum

!     output

      real*8 KS_eigenvalue(n_states)
      real*8 aux_eigenvector(n_nonsingular,n_states)

! INPUTS
! o n_nonsingular -- number of non-signular basis functions
! o n_states -- number of eigestates we want to calculate
! o trafo_hamiltonian -- Hamiltonian matrix in basis where singular basis functions are projected out.
! o safe_minimum -- accuracy of the eigenvalues
!
!  OUTPUT
! o KS_eigenvalue -- eigenvalues
! o aux_eigenvector -- eigenvectors in basis where singular basis functions are projected out.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
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
      real*8 v_lower, v_upper
      real*8 abs_tol

!     dspevx output

!     n_found: total number of (overlap matrix) eigenvalues found
!     eigenvalues : The desired eigenvalues, including those to be discarded
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer info

!  begin work

      if (myid.eq.0) then
         write(use_unit,'(2X,A,A)') "Solving standard eigenvalue problem for transformed Hamiltonian."
      end if

      if(n_states > n_nonsingular)then
         write(use_unit,*) '* Error: The overlap matrix is near-singular, and a reduction of the Hamiltonian size'
         write(use_unit,*) '* to ', n_nonsingular, ' is required after singular value decomposition.'
         write(use_unit,*) '* Unfortunately, this conflicts with the presently set number of states to be found '
         write(use_unit,*) '* (occupied plus empty states; for example, in GW, MP2 et al., you should ask for '
         write(use_unit,*) '* exactly all possible states back, so that the sum of occupied and empty states '
         write(use_unit,*) '* equals exactly the Hamiltonian size after singular value dcomposition, ', n_nonsingular, '.'
         write(use_unit,*) '* The number of states asked, ', n_states,', is now larger than the number of nonsingular basis', &
              ' functions', n_nonsingular
         write(use_unit,*) '* Please reduce empty_states and restart calculation.'
         write(use_unit,*) '* Aborting.'
         stop
      end if

!     initialize input for dspevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     dsygvx will calculate the lowest n_states eigenvalues.
!     v_lower, v_upper are dummies for this case.
      range = 'I'
      i_lower = 1
      i_upper = n_states
      v_lower = 0.d0
      v_upper = 0.d0

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

!     notice that these allocations can create a problem if the basis size changes during the run

      if (.not.allocated(work)) then
        allocate(work(8*n_nonsingular))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_nonsingular))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_nonsingular))
      end if

!     dspevx is the lapack driver (expert mode) which solves a standard
!     real symmetric eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call dspevx &
      ( jobz, range, uplo, n_nonsingular, trafo_hamiltonian, &
        v_lower, v_upper,  i_lower, i_upper, abs_tol, &
        n_found, KS_eigenvalue, aux_eigenvector, n_nonsingular, &
        work, iwork, ifail, info &
      )

!  consistency checks

      if (info.lt.0) then
         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver dspevx():"
            write(use_unit,*) "The ", -info, "th argument in dspevx() had an ", &
                 "illegal value. Check."
         end if
         stop

      else if (info.gt.0) then

         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver dspevx():"
            write(use_unit,*) info, " eigenvectors failed to converge."
            write(use_unit,*) "ifail: ", ifail
         end if

         stop
      end if

      if (n_found.ne.n_states) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A)') "* Not enough eigenvalues found."
         end if

         stop
      end if

!  that's all folks

      end subroutine diagonalize_hamiltonian
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_hamiltonian_complex
!  NAME
!   diagonalize_hamiltonian_complex
!  SYNOPSIS

      subroutine diagonalize_hamiltonian_complex &
      ( n_nonsingular,  n_states, trafo_hamiltonian, &
        safe_minimum, KS_eigenvalue, aux_eigenvector &
      )

!  PURPOSE
!  Subroutine diagonalize_hamiltonian_complex is a wrapper around the expert Lapack
!  driver to determine the eigenvalues of the transformed (non-singular) Hamiltonian.
!
!  USES
      implicit none
!  ARGUMENTS

!     input

      integer n_nonsingular
      integer n_states
      complex*16 trafo_hamiltonian(n_nonsingular*(n_nonsingular+1)/2)
      real*8 safe_minimum

!     output

      real*8 KS_eigenvalue(n_states)
      complex*16 aux_eigenvector(n_nonsingular,n_states)

! INPUTS
! o n_nonsingular -- number of non-singular basis functions
! o n_states -- number of eigenstates
! o trafo_hamiltonian -- Hamiltonian matrix in basis where singular basis functions are projected out.
! o safe_minimum -- accuracy of the eigenvalues
!
!  OUTPUT
! o KS_eigenvalue -- eigenvalues
! o aux_eigenvector -- eigenvectors in basis where singular basis functions are projected out.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

!     zhpevx input

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

!     zhpevx output

!     n_found: total number of eigenvalues found
!     eigenvalues : The desired eigenvalues
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer n_states_at_work
      integer info
      integer i_state

!  begin work

      if (.not.suppress_further_ill_cond_output.and.myid.eq.0) then
         write(use_unit,'(2X,A,A)') "Solving standard eigenvalue problem for transformed Hamiltonian."
      end if

      n_states_at_work = n_states
      if(n_states > n_nonsingular)then
! XR: unfortunately we have to deal with this strange situation. My current fix is: when this happens
!     set the eigenvectors between n_nonsingular and n_states to be zero, so that the post-DFT calculations
!     are not affected by these states

        n_states_at_work = n_nonsingular

        if (.not.suppress_further_ill_cond_output) then
          write(use_unit,*) '* Warning: states from ', n_nonsingular+1, ' to ', n_states, ' will be effectively '
          write(use_unit,*) '* excluded in post-DFT calculations because of ill-conditioning.'
        end if

!         write(use_unit,*) '* Error: The overlap matrix is near-singular, and a reduction of the Hamiltonian size'
!         write(use_unit,*) '* to ', n_nonsingular, ' is required after singular value decomposition.'
!         write(use_unit,*) '* Unfortunately, this conflicts with the presently set number of states to be found '
!         write(use_unit,*) '* (occupied plus empty states; for example, in GW, MP2 et al., you should ask for '
!         write(use_unit,*) '* exactly all possible states back, so that the sum of occupied and empty states '
!         write(use_unit,*) '* equals exactly the Hamiltonian size after singular value decomposition, ', n_nonsingular, '.'
!         write(use_unit,*) '* The number of states asked, ', n_states,', is now larger than the number of nonsingular basis', &
!              ' functions', n_nonsingular
!         write(use_unit,*) '* Please reduce empty_states and restart calculation.'
!         write(use_unit,*) '* Aborting.'
!         stop
      end if

      if (.not.suppress_further_ill_cond_output) then
        write(use_unit,'(A,I8,A)') "* Further output regarding ill-conditioning on task ", myid, " will be suppressed."
        suppress_further_ill_cond_output = .true.
      end if

!     initialize input for zhpevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     zhpevx will calculate the lowest n_states eigenvalues.
!     v_lower, v_upper are dummies for this case.
      range = 'I'
      i_lower = 1
      i_upper = n_states_at_work
      v_lower = 0.d0
      v_upper = 0.d0

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

!     notice that these allocations can create a problem if the basis size changes during the run
      if (.not.allocated(zwork)) then
        allocate(zwork(2*n_nonsingular))
      end if

      if (.not.allocated(work)) then
        allocate(work(7*n_nonsingular))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_nonsingular))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_nonsingular))
      end if

!     zhpevx is the lapack driver (expert mode) which solves a standard
!     hermitian eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call zhpevx &
      ( jobz, range, uplo, n_nonsingular, trafo_hamiltonian, &
        v_lower, v_upper,  i_lower, i_upper, abs_tol, &
        n_found, KS_eigenvalue, aux_eigenvector, n_nonsingular, &
        zwork, work, iwork, ifail, info &
      )

!  consistency checks

      if (info.lt.0) then
         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver zhpevx():"
            write(use_unit,*) "The ", -info, "th argument in zhpevx() had an ", &
                 "illegal value. Check."
         end if
         stop

      else if (info.gt.0) then

         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver zhpevx():"
            write(use_unit,*) info, " eigenvectors failed to converge."
            write(use_unit,*) "ifail: ", ifail
         end if

         stop
      end if

      if (n_found.ne.n_states_at_work) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A)') "* Not enough eigenvalues found."
         end if

         stop
      end if

      if (n_found .lt. n_states) then
        do i_state = n_found+1, n_states
           KS_eigenvalue(i_state) = 1.e9
           aux_eigenvector(:,i_state)=(0.d0,0.d0)
        enddo
      endif

      end subroutine diagonalize_hamiltonian_complex
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/diagonalize_rse_hamiltonian_complex
!  NAME
!   diagonalize_rse_hamiltonian_complex
!  SYNOPSIS

      subroutine diagonalize_rse_hamiltonian_complex &
      ( n_nonsingular,  n_states, trafo_hamiltonian, &
        safe_minimum, KS_eigenvalue, aux_eigenvector &
      )

!  PURPOSE
!  Subroutine diagonalize_hamiltonian_complex is a wrapper around the expert Lapack
!  driver to determine the eigenvalues of the transformed (non-singular) Hamiltonian.
!
!  USES
      implicit none
!  ARGUMENTS

!     input

      integer n_nonsingular
      integer n_states
      complex*16 trafo_hamiltonian(n_nonsingular*(n_nonsingular+1)/2)
      real*8 safe_minimum

!     output

      real*8 KS_eigenvalue(n_states)
      complex*16 aux_eigenvector(n_nonsingular,n_states)

! INPUTS
! o n_nonsingular -- number of non-singular basis functions
! o n_states -- number of eigenstates
! o trafo_hamiltonian -- Hamiltonian matrix in basis where singular basis functions are projected out.
! o safe_minimum -- accuracy of the eigenvalues
!
!  OUTPUT
! o KS_eigenvalue -- eigenvalues
! o aux_eigenvector -- eigenvectors in basis where singular basis functions are projected out.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

!     zhpevx input

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

!     zhpevx output

!     n_found: total number of eigenvalues found
!     eigenvalues : The desired eigenvalues
!     ifail : If there are non-converged eigenvectors, these are their indices.
!     info : Status value to indicate potential errors

      integer n_found
      integer n_states_at_work
      integer info
      integer i_state

!  begin work

      if (.not.suppress_further_ill_cond_output.and.myid.eq.0) then
         write(use_unit,'(2X,A,A)') "Solving standard eigenvalue problem for transformed Hamiltonian."
      end if

      n_states_at_work = n_states
      if(n_states > n_nonsingular)then
! XR: unfortunately we have to deal with this strange situation. My current fix is: when this happens
!     set the eigenvectors between n_nonsingular and n_states to be zero, so that the post-DFT calculations
!     are not affected by these states

        n_states_at_work = n_nonsingular

        if (.not.suppress_further_ill_cond_output) then
          write(use_unit,*) '* Warning: states from ', n_nonsingular+1, ' to ', n_states, ' will be effectively '
          write(use_unit,*) '* excluded in post-DFT calculations because of ill-conditioning.'
        end if

!         write(use_unit,*) '* Error: The overlap matrix is near-singular, and a reduction of the Hamiltonian size'
!         write(use_unit,*) '* to ', n_nonsingular, ' is required after singular value decomposition.'
!         write(use_unit,*) '* Unfortunately, this conflicts with the presently set number of states to be found '
!         write(use_unit,*) '* (occupied plus empty states; for example, in GW, MP2 et al., you should ask for '
!         write(use_unit,*) '* exactly all possible states back, so that the sum of occupied and empty states '
!         write(use_unit,*) '* equals exactly the Hamiltonian size after singular value decomposition, ', n_nonsingular, '.'
!         write(use_unit,*) '* The number of states asked, ', n_states,', is now larger than the number of nonsingular basis', &
!              ' functions', n_nonsingular
!         write(use_unit,*) '* Please reduce empty_states and restart calculation.'
!         write(use_unit,*) '* Aborting.'
!         stop
      end if

      if (.not.suppress_further_ill_cond_output) then
        write(use_unit,'(A,I8,A)') "* Further output regarding ill-conditioning on task ", myid, " will be suppressed."
        suppress_further_ill_cond_output = .true.
      end if

!     initialize input for zhpevx

!     Compute both eigenvalues and eigenvectors.
      jobz = 'V'

!     zhpevx will calculate the lowest n_states eigenvalues.
!     v_lower, v_upper are dummies for this case.
      range = 'I'
      i_lower = 1
      i_upper = n_states_at_work
      v_lower = 0.d0
      v_upper = 0.d0

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

!     notice that these allocations can create a problem if the basis size changes during the run
      if (.not.allocated(zwork)) then
        allocate(zwork(2*n_nonsingular))
      end if

      if (.not.allocated(work)) then
        allocate(work(7*n_nonsingular))
      end if

      if (.not.allocated(iwork)) then
        allocate(iwork(5*n_nonsingular))
      end if

      if (.not.allocated(ifail)) then
        allocate(ifail(n_nonsingular))
      end if

!     zhpevx is the lapack driver (expert mode) which solves a standard
!     hermitian eigenvalue problem using packed input matrices.
!     See the extensive comments in the header of that
!     subroutine for detailed explanations.

      call zhpevx &
      ( jobz, range, uplo, n_nonsingular, trafo_hamiltonian, &
        v_lower, v_upper,  i_lower, i_upper, abs_tol, &
        n_found, KS_eigenvalue, aux_eigenvector, n_nonsingular, &
        zwork, work, iwork, ifail, info &
      )

!  consistency checks

      if (info.lt.0) then
         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver zhpevx():"
            write(use_unit,*) "The ", -info, "th argument in zhpevx() had an ", &
                 "illegal value. Check."
         end if
         stop

      else if (info.gt.0) then

         if (myid.eq.0) then
            write(use_unit,*) "Eigenvalue solver zhpevx():"
            write(use_unit,*) info, " eigenvectors failed to converge."
            write(use_unit,*) "ifail: ", ifail
         end if

         stop
      end if

      if (n_found.ne.n_states_at_work) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A)') "* Not enough eigenvalues found."
         end if

         stop
      end if

      if (n_found .lt. n_states) then
        do i_state = n_found+1, n_states
           KS_eigenvalue(i_state) = 1.e9
           aux_eigenvector(:,i_state)=(0.d0,0.d0)
        enddo
      endif

      if (allocated(zwork)) then
        deallocate(zwork)
      end if

      if (allocated(work)) then
        deallocate(work)
      end if

      if (allocated(iwork)) then
        deallocate(iwork)
      end if

      if (allocated(ifail)) then
        deallocate(ifail)
      end if

      end subroutine diagonalize_rse_hamiltonian_complex
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/overlap_svd
!  NAME
!   overlap_svd
!  SYNOPSIS

      subroutine overlap_svd &
      ( n_basis,  basis_threshold, &
        n_nonsingular, eigenvalues, overlap_transform &
      )

!  PURPOSE
!  Subroutine overlap_svd is a wrapper around the expert Lapack
!  driver to determine the eigenvalues (and hence, singular values) of
!  the overlap matrix. The output is a matrix overlap_transform, which
!  contains only the non-singular eigenvectors of the overlap matrix.
!
!  USES
      implicit none
!  ARGUMENTS

      integer n_basis
      real*8 basis_threshold
      integer n_nonsingular
      real*8, dimension(n_basis) :: eigenvalues
      real*8 overlap_transform(n_basis,n_basis)

!  INPUTS
!  o n_basis -- number of basis functions
!  o basis_threshold -- minimum accepted eigenvalue
!  o overlap_transform -- On input, the fully expanded overlap matrix
!
!  OUTPUT
!  o n_nonsingular -- Number of non-singular eigenvectors of the overlap matrix
!  o eigenvalues -- Eigenvalues of the overlap matrix
!  o overlap_transform -- On output, overwritten by non-singular eigenvectors of the ovlp matrix
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

!  local variables

!     dgesvd input

!     matrices range_space and null_space are never referenced,
!     hence not allocated

      character :: jobu
      character :: jobvt

      real*8, dimension(:,:), allocatable :: range_space
      real*8, dimension(:,:), allocatable :: null_space

      integer :: l_work

!     dgesvd output

      integer :: info

      integer :: i_eval

!  begin work

      call localorb_info( &
           "Overlap matrix -- singular value decomposition.", &
           6,'(2X,A)' )

!     initialize input for dgesvd

!     return basis vectors which span the range of overlap_matrix in
!     working array overlap_transform
      jobu = 'O'

!     do not return the nullspace of overlap_matrix
      jobvt = 'N'

!     dummy allocations for matrices U, V**T, which are never referenced
      allocate (range_space(1,1))
      allocate (null_space(1,1))

      if (.not.allocated(work)) then
        allocate(work(1))
      end if

!     determine optimum working dimensions
      l_work = -1

      call dgesvd &
      ( jobu, jobvt, n_basis, n_basis, overlap_transform, &
        n_basis, eigenvalues, range_space, 1, null_space, 1, &
        work, l_work, info &
      )

      if (info.lt.0) then

         if (myid.eq.0) then
            write(use_unit,'(1X,A,I3,2A)') &
                 "* SVD work space query: The ", -info, &
                 "th argument to dgesvt had an illegal value. ",&
                 "Please check."
         end if

         stop
      else if (info.gt.0) then

         if (myid.eq.0) then
           write(use_unit,'(1X,A,I5,A,A)') &
           "* SVD work space query: ", info, &
           " superdiagonals of intermediate bilinear form failed to ", &
           "converge."
           write(use_unit,'(1X,A,A)') &
           "* This error message is almost certainly incorrect, as ", &
           "no SVD was yet performed - please check the code."
         end if

         stop
      end if

!     actual singular value decomposition
      l_work = work(1)
      if (l_work.gt.size(work)) then
        deallocate(work)
        allocate (work(l_work))
      end if

      call dgesvd &
      ( jobu, jobvt, n_basis, n_basis, overlap_transform, &
        n_basis, eigenvalues, range_space, 1, null_space, 1, &
        work, l_work, info &
      )

!     check proper execution of SVD
      if (info.lt.0) then

         if (myid.eq.0) then
           write(use_unit,'(1X,A,I3,A)') &
           "* SVD work space query: The ", -info, &
           "th argument to dgesvt had an illegal value. Please check."
         end if

         stop
      else if (info.gt.0) then

         if (myid.eq.0) then
           write(use_unit,'(1X,A,I5,A,A)') &
          "* SVD work space query: ", info, &
          " superdiagonals of intermediate bilinear form failed to ", &
          "converge."
         end if

      end if

!     determine the number of trustworthy (non-singular) eigenvectors of the overlap matrix
!     we could do this via the condition number but for some reason the largest singular
!     values are always 2 anyway so we can use basis_threshold directly.

      n_nonsingular = 0
      do i_eval = 1, n_basis, 1
        if (eigenvalues(i_eval).gt.basis_threshold) then
          n_nonsingular = n_nonsingular+1
        end if
      enddo

      if (n_nonsingular.lt.n_basis) then

         if (myid.eq.0) then
            write(use_unit,'(2X,A,I8,A)') "Overlap matrix is singular:"
            write(use_unit,'(2X,A,I8,A,I8,A)') &
                 "| Using ", n_nonsingular, &
                 " out of a possible ", n_basis, &
                 " specified basis functions."
         end if

         if (.not.override_illconditioning) then
            call stop_illconditioning()
         end if

         if(n_states.gt.n_nonsingular) then
           n_states =  n_nonsingular
         endif

      end if

      deallocate(range_space)
      deallocate(null_space)

      end subroutine overlap_svd
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/solve_pulay
!  NAME
!    solve_pulay
!  SYNOPSIS

      subroutine solve_pulay(n_max_pulay, pulay_saved_iter, &
      &                      pulay_matrix, pulay_vector, mixing_factor)

        !  PURPOSE
        !
        !    Subroutine pulay_solve is a wrapper around the lapack expert
        !    solver for a symmetric indefinite problem. We solve Eq. (90) of
        !    Kresse and Furthmuellers Comp. Mat. Sci. 6 (1996) 15-50
        !
        !    There is no real reason for this routine to be present in the
        !    lapack_wrapper module as it does not need much workspace, unlike
        !    the other solvers which concern the basis. It is kept here only
        !    for consistency.
        !
        !    For numerical stability in border cases, scale the matrix to have
        !    unity diagonals and use least squares.
        !
        !  USES
        use aims_memory_tracking, only : aims_allocate, aims_deallocate
        implicit none

        !  ARGUMENTS

        integer, intent(IN) :: n_max_pulay
        integer, intent(IN) :: pulay_saved_iter
        real*8, intent(IN) :: pulay_matrix(n_max_pulay, n_max_pulay)
        real*8, intent(IN) :: pulay_vector(pulay_saved_iter)
        real*8, intent(OUT) :: mixing_factor(pulay_saved_iter)

        ! INPUTS
        !   o n_max_pulay -- Actual size of arrays
        !   o pulay_saved_iter -- Logical size of arrays
        !   o pulay_matrix -- Matrix of LSE
        !   o pulay_vector -- Right hand side of LSE
        !
        ! OUTPUT
        !   o mixing_factor -- Solution of LSE
        !  AUTHOR
        !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
        !  HISTORY
        !    Release version, FHI-aims (2008).
        !  SOURCE

        !  local variables

        real*8, parameter :: r_cond = 1d-14
        integer :: l_work, i, j, rank
        real*8 :: work_tmp(1)
        integer info
        character*150 :: info_str
        character(*), parameter :: func = 'solve_pulay'

        !  begin work

        if (allocated(tmp_matrix)) then
           if (n_max_pulay**2 > size(tmp_matrix)) then
              call aims_deallocate( tmp_matrix, "tmp_matrix" )
           end if
        end if
        if (.not.allocated(tmp_matrix)) then
           call aims_allocate( tmp_matrix, n_max_pulay, n_max_pulay, "tmp_matrix" )
        end if

        if (allocated(i_piv)) then
           if (n_max_pulay > size(i_piv)) call aims_deallocate( i_piv, "i_piv" )
        end if
        if (.not.allocated(i_piv)) then
           call aims_allocate( i_piv, n_max_pulay, "i_piv" )
           i_piv = 0    ! dgelsy would permute according to this, otherwise.
        end if

        ! tmp_matrix  <-  scaled Pulay matrix
        ! The scaling is equivalent to using normalized charge differences.
        do i = 1, pulay_saved_iter
           if (pulay_matrix(i, i) <= 0.d0) then   ! Taken care of by caller
              call aims_stop('Ivalid diagonal in Pulay matrix', func)
           end if
           do j = i, pulay_saved_iter
              tmp_matrix(i, j) = pulay_matrix(i, j) / &
              &                  sqrt(pulay_matrix(i, i) * pulay_matrix(j, j))
              tmp_matrix(j, i) = tmp_matrix(i, j)
           end do
           mixing_factor(i) = pulay_vector(i) / sqrt(pulay_matrix(i, i))
        end do

        ! determine optimal work space size for lapack solver

        call DGELSY(pulay_saved_iter, pulay_saved_iter, 1, &
        &           tmp_matrix, n_max_pulay, &
        &           mixing_factor, n_max_pulay, &
        &           i_piv, r_cond, rank, &
        &           work_tmp, -1, info)
        if (info /= 0) then
           write(info_str, "('dgelsy workspace query: info =',I6)") info
           call aims_stop(info_str, func)
        end if

        l_work = work_tmp(1)
        if (allocated(pulay_work)) then
           if (l_work > size(pulay_work)) call aims_deallocate( pulay_work, "pulay_work" )
        end if
        if (.not.allocated(pulay_work)) then
           call aims_allocate( pulay_work, l_work, "pulay_work" )
        end if

        call DGELSY(pulay_saved_iter, pulay_saved_iter, 1, &
        &           tmp_matrix, n_max_pulay, &
        &           mixing_factor, n_max_pulay, &
        &           i_piv, r_cond, rank, &
        &           pulay_work, l_work, info)

        if (info /= 0) then
           write(info_str, "('dgelsy: info =',I6)") info
           call aims_stop(info_str, func)
        end if

        if (rank < pulay_saved_iter) then
           write(info_str, "(2X,'* ',A,'; rank =',I6,'; saved_iter =',I6)") &
           & 'Rank deficiency in solve_pulay', rank, pulay_saved_iter
           call localorb_info(info_str)
        end if

        ! Scale back to unnormalized charge differences
        do i = 1, pulay_saved_iter
           mixing_factor(i) = mixing_factor(i) / sqrt(pulay_matrix(i, i))
        end do

        ! Output Pulay mixing factors as diagnostic information for output_level full
        ! This is information that is rarely needed but will be helpful to understand
        ! cases where s.c.f. convergence is subtly problematic.

        if (output_priority.lt.1) then
           write(info_str, "(2X,A)") &
              & 'Diagnostic output of the density mixing factors determined by the Pulay scheme:'
              call localorb_info(info_str,use_unit,'(A)',OL_low)

           write(info_str, "(2X,A,A)") &
              & '| Note that these factors are applied to density residuals, ', &
              & 'i.e., to density differences that are normalized to zero.'
              call localorb_info(info_str,use_unit,'(A)',OL_low)

           write(info_str, "(2X,A)") &
              & '| Iteration         Mixing factor  '
              call localorb_info(info_str,use_unit,'(A)',OL_low)

           do i = 1, pulay_saved_iter
              write(info_str, "(2X,'| N minus ',I3,4X,E15.6)") &
                 & i, mixing_factor(i)
                 call localorb_info(info_str,use_unit,'(A)',OL_low)
              
           end do
           
        end if
        
      end subroutine solve_pulay
!******
!-------------------------------------------------------------------------------
!****s* lapack_wrapper/cleanup_lapack_wrapper
!  NAME
!   cleanup_lapack_wrapper
!  SYNOPSIS

      subroutine cleanup_lapack_wrapper &
       ( )
!  PURPOSE
!  This is the real purpose of this module, to allow for an organized deallocation
!  of auxiliary arrays
!
!  USES
     use aims_memory_tracking, only : aims_deallocate
     implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      if (allocated(zwork)) then
        deallocate(zwork)
      end if
      if (allocated(work)) then
        deallocate(work)
      end if
      if (allocated(iwork)) then
        deallocate(iwork)
      end if
      if (allocated(ifail)) then
        deallocate(ifail)
      end if

      if (allocated(tmp_matrix)) then
        call aims_deallocate( tmp_matrix, "tmp_matrix" )
      end if
      if (allocated(i_piv)) then
        call aims_deallocate( i_piv,           "i_piv" )
      end if
      if (allocated(pulay_work)) then
        call aims_deallocate( pulay_work, "pulay_work" )
      end if

      end subroutine cleanup_lapack_wrapper
!******
!-------------------------------------------------------------------------------
      end module lapack_wrapper
