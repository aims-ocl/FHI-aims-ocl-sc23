!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Subroutines for performing Pulay mixing for some quantity X. X can
!!  be a density (with an arbitrary number of spin channels), density
!!  matrix, or in general any array. This module was originally
!!  designed for the DFPT* and MagneticResponse parts, but should be
!!  generic enough to be used elsewhere. ATTENTION: for the original
!!  Pulay mixer used in the initial SCF cycle, see mixing.f90. There
!!  might be other mixers lurking in the code. If you are unsure of
!!  which one to use, contact any of the authors of these subroutines
!!  (see git log or git blame).
!!
!!  Here, the implementation is based on Kresse and Furthmüller,
!!  Comput. Mater. Sci. 6, 15 (1996) [Eqs. (88)-(92)]. We recast those
!!  equations into a slightly different form, which makes array
!!  indexing easier. The optimal input value for X is expressed as
!!
!!       X_in^opt = X_in^1 + G R^1
!!             + Sum_{i=1}^{Nmax-1} alpha^i(Delta_X^i + G Delta_R^i),
!!      Delta_X^i = X_in^i - X_in^{i+1},
!!            R^i = X_out^i[X_in^i] - X_in^i,
!!      Delta_R^i = R^i - R^{i+1},
!!        A alpha + <Delta_R | R^1> = 0,
!!           A_ij = <Delta_R^j | Delta_R^i>,
!!
!!  where the superscripts count away from the current SCF loop. For
!!  instance, X_in^1 is the input value of the current loop, X_in^2
!!  that of the previous loop, and so on. 'A alpha' refers to a
!!  matrix-vector product. Nmax is the number of Pulay steps
!!  (n_pulay_steps below). G is the linear mixing parameter
!!  (linear_mix_param argument to pulay_mix). Other quantities in
!!  these equations are explained below.
!!
!!  Usage:
!!
!!   i) Call pulay_mix when mixing is required. The first call
!!      initalizes the global work variables.
!!  ii) Call cleanup_pulay_mixing.
!!
!!  COMMENTS
!!
!!  For linear mixing, call pulay_mix with n_pulay_steps=1.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module DFPT_pulay_mixing

  use aims_memory_tracking,  only: aims_allocate
  use mpi_tasks,             only: check_allocation
  use synchronize_mpi_basic, only: sync_vector
  use tools,                 only: safe_deallocate

  implicit none

  private
  public :: pulay_mix, cleanup_pulay_mixing

  ! X_in^i - input values of the current and previous loops. Note that
  ! the last dimension is 2 instead of n_pulay_steps here, i.e., only
  ! this and the previous input value are required for the
  ! algorithm. See the X_size argument to pulay_mix. (X_size,2)
  real*8, allocatable :: X_in(:,:)
  ! R^i - residuals corresponding to the current and previous
  ! loop. Note again that the last dimension is 2 here.
  ! (X_size,2)
  real*8, allocatable :: residuals(:,:)
  ! Delta_X^i - differences of input values. See the n_pulay_steps
  ! argument to pulay_mix. (X_size,n_pulay_steps-1)
  real*8, allocatable :: X_in_diffs(:,:)
  ! Delta_R^i - differences of residuals. (X_size,n_pulay_steps-1)
  real*8, allocatable :: residual_diffs(:,:)
  ! A_ij. (n_pulay_steps-1,n_pulay_steps-1)
  real*8, allocatable :: A_matrix(:,:)
  ! alpha_i - The mixing coefficients. (n_pulay_steps-1)
  real*8, allocatable :: alphas(:)
  ! This array is used for indexing some of the work arrays. Instead
  ! of explicitly permuting entries in the work arrays, we permute the
  ! indexing array itself at the end of each call to
  ! pulay_mix. (n_pulay_steps-1)
  integer, allocatable :: indices(:)

contains

  !!  FUNCTION
  !!
  !!  Allocates the work variables. This subroutine is automatically
  !!  called with the first call to pulay_mix.
  !!
  subroutine initialize_pulay_mixing(X_size, n_pulay_steps)
    integer, intent(in) :: X_size, n_pulay_steps
    character(*), parameter :: &
         & THIS_SUB = 'DFPT_pulay_mixing::initialize_pulay_mixing::'
    integer :: i_tmp
    call aims_allocate(X_in, X_size, 2, THIS_SUB//'X_in')
    call aims_allocate(residuals, X_size, 2, THIS_SUB//'residuals')
    call aims_allocate(X_in_diffs, X_size, n_pulay_steps-1, &
         & THIS_SUB//'X_in_diffs')
    call aims_allocate(residual_diffs, X_size, n_pulay_steps-1, &
         & THIS_SUB//'residual_diffs')
    call aims_allocate(indices, n_pulay_steps-1, THIS_SUB//'indices')
    X_in = 0d0
    indices = [(i_tmp, i_tmp = 1, size(indices))]
  end subroutine initialize_pulay_mixing

  !!  FUNCTION
  !!
  !!  Applies the Pulay algorithm for finding the new optimal input
  !!  value for quantity X. Also permutes the indices array for the
  !!  next iteration.
  !!
  !!  COMMENTS
  !!
  !!  Below there are two calls to sync_vector, which implicitly
  !!  assumes that X_out is a block-cyclic matrix and X_size \propto
  !!  size(eigenvec). One might suppose then that when Scalapack is
  !!  not used, this leads to wrong results because syncing across
  !!  cpus is not necessary. Actually, since the A matrix and alphas
  !!  appear as a product in the above equations, syncing both arrays
  !!  does not hurt and we still get the correct results. This
  !!  subroutine is thus completely general and works with arrays of
  !!  any memory layout (simple and 2D block-cyclic).
  !!
  subroutine pulay_mix(X_out, X_size, iter, n_pulay_steps, linear_mix_param)
    integer, intent(in) :: X_size, iter, n_pulay_steps
    ! On entry, X_out is the output of solving equations of interest
    ! (e.g., the KS equations for a ground state calculation). On
    ! exit, this is the mixed quantity.
    real*8, intent(in out) :: X_out(X_size)
    real*8, intent(in) :: linear_mix_param
    character(*), parameter :: THIS_SUB = 'DFPT_pulay_mixing::pulay_mix::'
    integer :: i_row, i_col
    ! dgelsy work variables
    integer :: i_tmp, i_tmp_1, i_tmp_2
    integer, allocatable :: i_work(:)
    real*8, allocatable :: r_work(:),S_matr(:)
    ! Initialize work variables if not yet initialized
    if (.not. allocated(X_in)) &
         & call initialize_pulay_mixing(X_size, n_pulay_steps)
    ! Update R
    residuals(:,2) = residuals(:,1)
    residuals(:,1) = X_out - X_in(:,1)
    not_first_iter: if (iter > 1 .and. n_pulay_steps > 1) then
       ! Update Delta_X and Delta_R
       X_in_diffs(:,indices(1)) = X_in(:,1) - X_in(:,2)
       residual_diffs(:,indices(1)) = residuals(:,1) - residuals(:,2)
       ! Set up the Pulay matrix
       if (iter <= n_pulay_steps) then
          call safe_deallocate(A_matrix)
          call safe_deallocate(alphas)
          call aims_allocate(A_matrix, iter-1, iter-1, THIS_SUB//'A_matrix')
          call aims_allocate(alphas, iter-1, THIS_SUB//'alphas')
       end if
       do i_col = 1, size(A_matrix,1)
          do i_row = 1, size(A_matrix,2)
             A_matrix(i_row,i_col) = sum(residual_diffs(:,indices(i_row))* &
                  & residual_diffs(:,indices(i_col)))
          end do
       end do
       call sync_vector(A_matrix, size(A_matrix))
       ! Solve a linear least squares problem to find the new
       ! coefficients alpha_i.
       !allocate(i_work(size(alphas)), r_work(1))
       allocate(i_work(1), r_work(1),S_matr(size(alphas)))
       do i_row = 1, size(alphas,1)
          ! At this point, 'alphas' refers to <Delta_R | R^1>.
          alphas(i_row) = -sum(residual_diffs(:,indices(i_row))*residuals(:,1))
       end do
       call dgelsd(size(alphas), size(alphas), 1, A_matrix, size(alphas), &
            & alphas, size(alphas), S_matr, 1d-14, i_tmp_1, r_work, -1, i_work, &
            & i_tmp_2)
       i_tmp = int(r_work(1),4)
       deallocate(r_work)
       allocate(r_work(i_tmp))
       !call dgelsy(size(alphas), size(alphas), 1, A_matrix, size(alphas), &
       !     & alphas, size(alphas), i_work, 1d-14, i_tmp, r_work, &
       !     & size(r_work), i_tmp)
       deallocate(i_work)
       allocate(i_work(3*size(alphas)*10 + 11*size(alphas)))
       call dgelsd(size(alphas), size(alphas), 1, A_matrix, size(alphas), &
            & alphas, size(alphas), S_matr, 1d-14, i_tmp_1, r_work, size(r_work), i_work, i_tmp_2)
       call sync_vector(alphas, size(alphas))
    end if not_first_iter
    ! Perform mixing
    X_out = X_in(:,1)+linear_mix_param*residuals(:,1)
    do i_row = 1, min(iter,n_pulay_steps)-1
       X_out = X_out + alphas(i_row)*(X_in_diffs(:,indices(i_row)) + &
            & linear_mix_param*residual_diffs(:,indices(i_row)))
    end do
    ! The mixed value becomes the next input value
    X_in(:,2) = X_in(:,1)
    X_in(:,1) = X_out
    ! Permute the indices for the next iteration
    if (size(indices) > 1) &
         & indices = [indices(size(indices)), indices(:size(indices)-1)]
  end subroutine pulay_mix

  subroutine cleanup_pulay_mixing()
    call safe_deallocate(X_in)
    call safe_deallocate(residuals)
    call safe_deallocate(X_in_diffs)
    call safe_deallocate(residual_diffs)
    call safe_deallocate(A_matrix)
    call safe_deallocate(alphas)
    call safe_deallocate(indices)
  end subroutine cleanup_pulay_mixing
end module DFPT_pulay_mixing
