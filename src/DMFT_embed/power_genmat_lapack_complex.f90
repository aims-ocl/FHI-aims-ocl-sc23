!****s* FHI-aims/power_genmat_lapack
!  power_genmat_lapack
!    power_genmat_lapack
!  SYNOPSIS

subroutine power_genmat_lapack_complex(n_size, matrix, power, safe_minimum, threshold, name)

  !  PURPOSE
  !     Take the power of matrix by diagonalizing.  Throw away all eigenmodes
  !     with eigenvalues smaller than threshold.
  !  USES

  use localorb_io, only: use_unit
  use mpi_tasks
  implicit none

  !  ARGUMENTS

  integer, intent(IN) :: n_size
  complex*16, intent(INOUT) :: matrix(n_size, n_size)
  real*8, intent(IN) :: power
  real*8, intent(IN) :: safe_minimum, threshold
  character*(*), intent(IN) :: name

  !  INPUTS
  !  o  n_size -- dimension of the matrix
  !  o  matrix -- the actual (symmetric) matrix
  !  o  power -- the power to take
  !  o  safe_minimum -- the absolute tolerance for the eigenvalue precision needed by
  !         the lapack eigenvalue solver.
  !  o  threshold -- the cutoff threshold for the eigenvalues of the concerned 
  !         overlap or Coulomb matrix, only those above this threshold will be included.
  !  o  name -- name of matrix (for output only); Choose '' for no output.
  !
  !  OUTPUTS
  !  o  matrix -- number,  number of nonsingular eigenvales larger than
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
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2010).
  !  SOURCE

  real*8, allocatable :: eigenvalues(:)
  complex*16, allocatable :: transform(:,:)
  real*8 :: ev_sqrt
  real*8, parameter :: my_thres = -1d10  ! ~ - huge(my_thres)
  integer :: n_nonsingular, i_basis_1, i_basis_2
  character(*), parameter :: func = 'power_genmat_lapack'

  allocate(eigenvalues(n_size))
  allocate(transform(n_size,n_size))
  eigenvalues=0.d0
  transform=0.d0

  matrix = - matrix  ! lapack sorts from small to large

  call diagonalize_auxmat_lapack_complex(n_size, matrix, safe_minimum, &
  & my_thres, n_nonsingular, eigenvalues, transform, '')

  eigenvalues = - eigenvalues  ! now we get from large to small

  if (n_nonsingular /= n_size) then
     write(use_unit,*) n_size, n_nonsingular
     write(use_unit,*) eigenvalues
     call aims_stop('Unphysical eigenvalues', func)
  end if

  do i_basis_1 = 1, n_size
     if (eigenvalues(i_basis_1) <= threshold) exit
  end do
  n_nonsingular = min(i_basis_1-1, n_size)

  if(name /= '') then
     write(use_unit,'(2X,A,I6,3A,ES14.4,A,ES14.4,A)') &
     & "Task", myid, &
     & ": Eigenvalues of the ", trim(name), " matrix range from", &
     & eigenvalues(1), " to", eigenvalues(n_size), "."
     write(use_unit,'(2X,A,I6,A,I8,A,I8,3A)') &
     & "Task", myid, ": Using ", n_nonsingular, "   eigenvalues out of rank ", &
     & n_size, "   ", trim(name), " matrix (auxiliary basis)."
     if (n_nonsingular < n_size) then
        write(use_unit,'(2X,A,I6,A,ES14.4,A,ES14.4,3A)') &
        & "Task", myid, &
        & ": Still using eigenvalue ", eigenvalues(n_nonsingular), &
        & " while cutting ", eigenvalues(n_nonsingular+1), &
        & " in ", trim(name), " matrix."
     end if
     write(use_unit,*)
  end if
  if (eigenvalues(n_nonsingular).lt. 1.d-5 .and. power < 0.d0) then
     write(use_unit,'(2X,A,I6,3A)') &
     & "Task", myid, &
     & ": Be careful! The ", trim(name), " matrix may be ill-conditioned."
     write(use_unit,*)
  end if

  do i_basis_1 = 1, n_nonsingular
     ev_sqrt = sqrt(eigenvalues(i_basis_1))
     do i_basis_2 = 1, n_size, 1
        transform(i_basis_2, i_basis_1) = &
        transform(i_basis_2,i_basis_1) * ev_sqrt**power
     enddo
  enddo

  matrix = (0.d0,0.d0)
  call zgemm('N', 'C', n_size, n_size, n_nonsingular, &
  &          (1.0d0,0.d0), transform(:,1:n_nonsingular), n_size, &
  &                 transform(:,1:n_nonsingular), n_size, &
  &           (0.d0,0.d0), matrix, n_size)

  deallocate(eigenvalues, transform)

end subroutine power_genmat_lapack_complex
!******
