subroutine solve_eigen_supercell(S,H, &
           KS_eigenvalue_supercell,KS_eigenvector_supercell)
!DSYGV computes all the eigenvalues, and optionally, the eigenvectors
! of a real generalized symmetric-definite eigenproblem,

use dimensions

implicit none

real*8 S(n_basis_sc_DFPT,n_basis_sc_DFPT),H(n_basis_sc_DFPT,n_basis_sc_DFPT)
integer itype,i
character*1 jobz
real*8 KS_eigenvalue_supercell(n_basis_sc_DFPT)
real*8 KS_eigenvector_supercell(n_basis_sc_DFPT,n_basis_sc_DFPT)
complex*16, dimension(:),allocatable :: work
real*8, dimension(:), allocatable :: rwork
integer lwork
integer info

itype=1 !   A*x = (lambda)*B*x
jobz='V' !  Compute eigenvalues and eigenvectors
lwork=-1
if(.not.allocated(work))  allocate(work(1))
call dsygv(itype,jobz,'L',n_basis_sc_DFPT,H,n_basis_sc_DFPT,S,n_basis_sc_DFPT,&
           KS_eigenvalue_supercell,work,lwork,info)

lwork = NINT(DBLE(work(1)))
deallocate(work)
allocate(work(lwork))
call dsygv(itype,jobz,'L',n_basis_sc_DFPT,H,n_basis_sc_DFPT,S,n_basis_sc_DFPT,&
           KS_eigenvalue_supercell,work,lwork,info)

do i = 1,n_basis_sc_DFPT
KS_eigenvector_supercell(1:n_basis_sc_DFPT,i)=H(1:n_basis_sc_DFPT,i)
enddo

end subroutine

