subroutine construct_solve_bse_mat(val, cond, qpe, coulomb_4ks_vcvc, coulomb_4ks_vccv, screened_coulomb_4ks_vvcc, screened_coulomb_4ks_vccv)
  use dimensions
  use physics
  use prodbas
  use constants
  use runtime_choices
  use synchronize_mpi
  use quicksort_index, only: dquicksort_indexlist
  implicit none

  integer :: val, cond
  real*8, dimension(n_states) :: qpe
  real*8, dimension(val, cond, val, cond) :: coulomb_4ks_vcvc
  real*8, dimension(val, cond, cond, val) :: coulomb_4ks_vccv
  real*8, dimension(val, val, cond, cond) :: screened_coulomb_4ks_vvcc
  real*8, dimension(val, cond, cond, val) :: screened_coulomb_4ks_vccv

  integer :: dms, diag, v, c, dms2, v1, c1, v2, c2, row, col, i
  real*8, allocatable, dimension(:, :) :: A, B, AB
  real*8 :: temp_qpe_v
! parameters for dsyev in lapack
  integer :: errorflag
  integer, parameter :: LWORK = 10000000
  real*8, allocatable, dimension(:) :: eigvals, WR, WI
  real*8, allocatable, dimension(:, :) :: VL, VR
  real*8, dimension(LWORK) :: WORK
  real*8 :: bse_alpha
  integer, allocatable, dimension(:) :: sort_ind
! reduce
  integer :: reduced_occ_states, reduced_unocc_states
  integer :: occ_states, unocc_states, i_state, i_spin
! new diagonalizatino
 integer                   :: t_kl, t_ij, nmat
 real*8,allocatable :: amb_eigval(:),bigomega(:)
 real*8,allocatable :: bigx(:, :),bigy(:, :)
 real*8,allocatable :: amb_matrix(:, :)
 real*8,allocatable :: temp_amb_matrix(:, :)
 real*8,allocatable :: apb_matrix(:, :)  ! apb_matrix constains (A+B) in the input, however it used a matrix buffer after
 real*8,allocatable       :: eigenvalue(:)

! singlet or triplet
  if(bse_singlet) then
    bse_alpha = 2
  elseif(bse_triplet) then
    bse_alpha = 0
  else
    call aims_stop('bse_s_t keyword needs to be set either singlet or triplet')
  end if
  if(myid == 0) print*, 'bse_alpha', bse_alpha

! bse_lower_limit, bse_upper_limit
  bse_lower_limit = max(1, bse_lower_limit)
  bse_upper_limit = min(bse_upper_limit, n_states)

! reduced_occ_states, reduced_unocc_states
  i_spin = 1
  occ_states = 0
  unocc_states = 0
  do i_state = 1, n_states
    if(occ_numbers(i_state, i_spin, n_k_points) < 1e-6) then
      unocc_states = unocc_states + 1
    else
      occ_states = occ_states + 1
    end if
  end do
  reduced_occ_states = occ_states - bse_lower_limit + 1
  reduced_unocc_states = bse_upper_limit - occ_states


!  dms = val * cond
  dms = reduced_occ_states * reduced_unocc_states
  allocate(A(dms, dms))
  allocate(B(dms, dms))
  A = 0
  B = 0
  diag = 0
!  do v = 1, val
  do v = bse_lower_limit, occ_states
!    do c = 1, cond
    do c = 1, reduced_unocc_states
      diag = diag + 1
!      A(diag, diag) = qpe(val + c) - qpe(v)
      A(diag, diag) = qpe(occ_states + c) - qpe(v)
    end do
  end do

  row = 0
!  open(unit = 99, file = 'bse_mat_qpe_v')
!  do v1 = 1, val
  do v1 = bse_lower_limit, occ_states
!    do c1 = 1, cond
    do c1 = 1, reduced_unocc_states
      row = row + 1
      col = 0
!      do v2 = 1, val
      do v2 = bse_lower_limit, occ_states
!        do c2 = 1, cond
        do c2 = 1, reduced_unocc_states
          col = col + 1
!          temp_qpe_v = A(row, col) + 2 * coulomb_4ks_vcvc(v1, c1, v2, c2)
          A(row, col) = A(row, col) - screened_coulomb_4ks_vvcc(v1, v2, c1, c2)+&
                        bse_alpha * coulomb_4ks_vcvc(v1, c1, v2, c2)
!           write(99, *) temp_qpe_v
!          if(temp_qpe_v > 1e-6) write(99, *) row, col,temp_qpe_v
          B(row, col) = B(row, col) - screened_coulomb_4ks_vccv(v1, c2, c1, v2)+&
                        bse_alpha * coulomb_4ks_vccv(v1, c1, c2, v2)
        end do
      end do
    end do
  end do
!  close(99)
! check A
!  open(unit = 99, file = 'bse_mat')
!  do v = 1, dms
!    do c = 1, dms
!      if(abs(A(v, c)) > 1e-6) write(99, *) v, c, A(v, c)
!    end do
!  end do
!  close(99)

  dms2 = 2 * dms
  allocate(AB(dms2, dms2))
  do v = 1, dms
    do c = 1, dms
      AB(v, c) = A(v, c)
      AB(v + dms, c) = B(v, c)
      AB(v, c + dms) = (-1) * B(v, c)
      AB(v + dms, c + dms) = (-1) * A(v, c)
    end do
  end do

  nmat = dms
  allocate(amb_matrix(nmat,nmat))
  allocate(temp_amb_matrix(nmat,nmat))
  allocate(apb_matrix(nmat,nmat))
  apb_matrix = A + B
  amb_matrix = A - B

  print*, 'solve BSE matrix'
  allocate(eigvals(dms))
  call dsyev('V','U',dms,A,dms,eigvals,WORK,LWORK,errorflag)
  allocate(WR(dms2))
  allocate(WI(dms2))
  allocate(VL(dms2, dms2))
  allocate(VR(dms2, dms2))
  call dgeev('N','V', dms2, AB, dms2, WR, WI, VL, dms2, VR, dms2, WORK, LWORK, errorflag)
  print*, "************************************************************"
  print*, "bse result TDA, lowest 20"
  do i = 1, min(20, dms)
    print*, eigvals(i) * hartree
  end do
  print*, "bse result without TDA, lowest 20"
! **** sort
  allocate(sort_ind(dms2))
  do i = 1, dms2
    sort_ind(i) = i
  end do

  call dquicksort_indexlist(dms2, abs(WR), sort_ind)
  print*, "sorted:"
  do i = 1, min(40, dms2)
    if (WR(sort_ind(i)) > 0) print*, WR(sort_ind(i)) * hartree
  end do
  deallocate(sort_ind)
!  call sort(WR, dms2)
!  print*, "sorted:"
!  do i = 1, 20
!    print*, WR(i) * hartree
!  end do


! ********************* another diagonalization *********************
! nmat = dms
!! allocate(amb_matrix(nmat,nmat))
!! allocate(temp_amb_matrix(nmat,nmat))
!! allocate(apb_matrix(nmat,nmat))
! allocate(eigenvalue(nmat))
! allocate(bigx(nmat,nmat))
! allocate(bigy(nmat,nmat))
!!  apb_matrix = A + B
!!  amb_matrix = A - B
! ! First symmetrize the matrices since only the lower triangle was calculated
! do t_kl=1,nmat
!   do t_ij=t_kl+1,nmat
!     amb_matrix(t_kl,t_ij) = amb_matrix(t_ij,t_kl)
!     apb_matrix(t_kl,t_ij) = apb_matrix(t_ij,t_kl)
!   enddo
! enddo
!
!  allocate(amb_eigval(dms))
!!  print*, "amb_matrix", amb_matrix
!  print*, "calculating A - B eigenvalue"
!!  temp_amb_matrix = amb_matrix
!!  call ssyev('V','U',dms,amb_matrix,dms,amb_eigval,WORK,LWORK,errorflag)
!  call dsyev('V','U',dms,amb_matrix,dms,amb_eigval,WORK,LWORK,errorflag)
!!  print*, "amb_eigval", amb_eigval
! ! bigx contains the (A-B)**1/2
! ! bigy contains the (A-B)**-1/2
!!  amb_matrix = temp_amb_matrix
!!  print*, "amb_matrix", amb_matrix
!
! forall(t_kl=1:nmat)
!   bigx(:,t_kl) = amb_matrix(:,t_kl)*SQRT(amb_eigval(t_kl))
!   bigy(:,t_kl) = amb_matrix(:,t_kl)/SQRT(amb_eigval(t_kl))
! end forall
! deallocate(amb_eigval)
!
! amb_matrix = TRANSPOSE( amb_matrix )
! bigx(:,:) = MATMUL( bigx(:,:) , amb_matrix(:,:) )
! bigy(:,:) = MATMUL( bigy(:,:) , amb_matrix(:,:) )
!
! ! Use amb_matrix as a temporary matrix here:
! amb_matrix(:,:) = MATMUL( apb_matrix , bigx )
! apb_matrix(:,:)  = MATMUL( bigx, amb_matrix )
!
! allocate(bigomega(nmat))
! print*, "calculating big Omega"
!! call ssyev('V','U',dms,apb_matrix,dms,bigomega,WORK,LWORK,errorflag)
! call dsyev('V','U',dms,apb_matrix,dms,bigomega,WORK,LWORK,errorflag)
!! print*, "bigomega square"
!! print*, bigomega
! bigomega(:) = SQRT(bigomega(:))
!! print*, "bigomega"
!! print*, bigomega
!
! forall(t_kl=1:nmat)
!   apb_matrix(:,t_kl) = apb_matrix(:,t_kl) / SQRT(bigomega(t_kl))
!   eigenvalue(t_kl) = bigomega(t_kl)
! end forall
!
!  print*, "another bse result no TDA, lowest 20"
!  allocate(sort_ind(nmat))
!  do i = 1, nmat
!    sort_ind(i) = i
!  end do
!
!  call dquicksort_indexlist(nmat, eigenvalue, sort_ind)
!  do i = 1, 20
!    print*, eigenvalue(sort_ind(i)) * hartree
!  end do
!  deallocate(sort_ind)
!
! deallocate(amb_matrix)
!! deallocate(temp_amb_matrix)
! deallocate(apb_matrix)
! deallocate(eigenvalue)
! deallocate(bigx)
! deallocate(bigy)
!*************************
  print*, "* end of subroutine bse in bse/bse.f90"
  print*, "************************************************************"
  deallocate(WR)
  deallocate(WI)
  deallocate(VL)
  deallocate(VR)
  deallocate(eigvals)
  deallocate(A)
  deallocate(B)
  deallocate(AB)
end subroutine

subroutine sort(x, dms)
  implicit none
  integer, intent(in) :: dms
  real*8, dimension(dms), intent(inout) :: x
  integer :: location, i, j, newdms
  real*8 :: minimum, temp

  j = 1
  do i = 1, dms
    if (x(i) > 0.1) then
      x(j) = x(i)
      j = j + 1
    end if
  end do
  print*, "new dms, dms:", j - 1, dms

  newdms = j - 1
  do i = 1, newdms - 1
    minimum = x(i)
    location = i
    do j = i + 1, newdms
      if(x(j) < minimum) then
        minimum = x(j)
        location = j
      end if
    end do
    if(location /= i) then
      temp = x(i)
      x(i) = x(location)
      x(location) = temp
    end if
  end do
end subroutine sort
