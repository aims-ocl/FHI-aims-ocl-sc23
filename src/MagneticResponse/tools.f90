!!  COPYRIGHT
!!
!!  Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!  e.V. Please note that any use of the "FHI-aims-Software" is
!!  subject to the terms and conditions of the respective license
!!  agreement.
!!
!!  FUNCTION
!!
!!  Procedures for converting between packed and block cyclic matrices
!!  (both ways), wrappers for some of the commonly used Scalapack
!!  routines, and otherwise useful helper functions.
!!
!!  AUTHORS
!!
!!  FHI-aims team
!!
module tools

  use aims_memory_tracking,  only: aims_deallocate
  use dimensions,            only: n_basis, n_spin
  use load_balancing,        only: batch_perm, n_bp_integ, &
      & get_full_local_matrix, set_full_local_ham
  use mpi_tasks,             only: mpi_comm_global
  use pbc_lists,             only: column_index_hamiltonian, index_hamiltonian
  use scalapack_wrapper,     only: get_sparse_local_matrix_scalapack, &
       & ham, l_col, l_row, my_col, my_row, n_my_rows, n_my_cols, sc_desc, &
       & set_sparse_local_ham_scalapack, set_full_matrix_real
  use synchronize_mpi_basic, only: sync_vector
  use runtime_choices,       only: output_level, use_load_balancing, &
       & use_local_index, use_scalapack
  use types,                 only: dp

  implicit none

  private
  public :: antisymmetrize, batch_to_matrix_packed, block_cyclic_to_packed, &
       & if_present, libxc_yes_no, matrix_packed_to_batch, mul, newunit, &
       & norm2, packed_to_block_cyclic, safe_deallocate, start_wall, &
       & stop_wall, str, walltime_mat_conv

  ! Converts a number to string
  interface str
     module procedure int_to_str, real64_to_str
  end interface str

  ! Performs C = alpha*A*B + beta*C with Scalapack
  interface mul
     module procedure               &
          & pdgemm_wrapper,         &
          & pzgemm_wrapper,         &
          & pzgemm_im_re_wrapper,   &
          & pzgemm_re_im_wrapper,   &
          & pzgemm_re_re_im_wrapper
  end interface mul

  ! Calculates the norm2 of an array of dimension n
  interface norm2
     module procedure norm2_vector, norm2_matrix
  end interface norm2

  interface
     pure logical function libxc_yes_no() result(y)
     end function libxc_yes_no
  end interface

  ! Convenience functions for initializing optional arguments.
  interface if_present
     module procedure          &
          & if_present_int,    &
          & if_present_int_1d, &
          & if_present_bool,   &
          & if_present_char,   &
          & if_present_real,   &
          & if_present_complex
  end interface if_present

  ! Wrapper around the aims_deallocate routines
  interface safe_deallocate
     module procedure &
          & safe_deallocate_r1, safe_deallocate_r2, safe_deallocate_r3, &
          & safe_deallocate_r4, safe_deallocate_r5, safe_deallocate_i1, &
          & safe_deallocate_i2, safe_deallocate_i3
  end interface safe_deallocate

  ! Measures the wall time spent on subroutines that do conversion
  ! from one matrix type to another, i.e. packed_to_block_cyclic and
  ! block_cyclic_to_packed. Note: this should be reset manually when
  ! desirable and the value needs to be divided by the count rate
  ! obtained from system_clock afterwards.
  integer :: walltime_mat_conv

contains

  pure function int_to_str(x) result(y)
    integer, intent(in) :: x
    character(:), allocatable :: y
    character(128) :: tmp
    write(tmp, '(i0)') x
    y = trim(tmp)
  end function int_to_str

  pure function real64_to_str(x) result(y)
    real(dp), intent(in) :: x
    character(:), allocatable :: y
    character(128) :: tmp
    write(tmp, '(f20.4)') x
    y = trim(adjustl(tmp))
  end function real64_to_str

  !!  FUNCTION
  !!
  !!  Finds an unused file handle. (F2008 intrinsic)
  !!
  integer function newunit(unit) result(io_unit)
    integer, intent(out) :: unit
    logical :: unit_open
    do io_unit = 10, 100
       inquire (unit=io_unit, opened=unit_open)
       if (.not. unit_open) then
          unit = io_unit
          return
       end if
    end do
    io_unit = -1
  end function newunit

  !!  FUNCTION
  !!
  !!  Antisymmetrizes a block cyclic matrix.
  !!
  pure subroutine antisymmetrize(matrix, n_spins)
    real(dp), intent(in out) :: matrix(:,:,:)
    integer, intent(in) :: n_spins
    integer :: i_row, i_col
    do i_col = 1, size(matrix,2)
       do i_row = 1, size(matrix,1)
          if (my_row(i_row) < my_col(i_col)) &
               & matrix(i_row,i_col,:n_spins) = -matrix(i_row,i_col,:n_spins)
       end do
    end do
  end subroutine antisymmetrize

  !!  FUNCTION
  !!
  !!  Wrappers for matrix multiplications. If
  !!  use_scalapack==.true. and do_serial=.false., pdgemm or pzgemm
  !!  are called. Otherwise, dgemm or zgemm are called.
  !!
  pure subroutine pdgemm_wrapper(A, B, C, transa, transb, M, N, K, alpha, &
       & iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
    interface
       pure subroutine pdgemm(transa, transb, M, N, K, alpha, A, ia, ja, &
            & desca, B, ib, jb, descb, beta, C, ic, jc, descc)
         import dp
         real(dp), intent(in) :: A(*), B(*), alpha, beta
         real(dp), intent(in out) :: C(*)
         character, intent(in) :: transa, transb
         integer, dimension(9), intent(in) :: desca, descb, descc
         integer, intent(in) :: M, N, K, ia, ja, ib, jb, ic, jc
       end subroutine pdgemm
       pure subroutine dgemm(transa,transb,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
         import dp
         real(dp), intent(in) :: A, B, alpha, beta
         real(dp), intent(in out) :: C
         character, intent(in) :: transa, transb
         integer, intent(in) :: M, N, K, LDA, LDB, LDC
       end subroutine dgemm
    end interface
    real(dp), intent(in) :: A(:,:), B(:,:)
    real(dp), intent(in out) :: C(:,:)
    character, intent(in), optional :: transa, transb
    real(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: M, N, K, iA, jA, iB, jB, iC, jC
    integer, intent(in), optional :: descA(9), descB(9), descC(9)
    logical, intent(in), optional :: do_serial
    character :: transa_loc, transb_loc
    real(dp) :: alpha_loc, beta_loc
    integer :: M_loc, N_loc, K_loc, iA_loc, jA_loc, descA_loc(9), iB_loc, &
         & jB_loc, descB_loc(9), iC_loc, jC_loc, descC_loc(9)
    logical :: do_parallel
    do_parallel = use_scalapack .and. .not. if_present(do_serial, .false.)
    alpha_loc = if_present(alpha, 1d0)
    beta_loc = if_present(beta, 0d0)
    transa_loc = if_present(transa, 'n')
    transb_loc = if_present(transb, 'n')
    descA_loc = if_present(descA, sc_desc)
    descB_loc = if_present(descB, sc_desc)
    descC_loc = if_present(descC, sc_desc)
    if (do_parallel) then
       if (transa_loc == 'n' .or. transa_loc == 'N') then
          M_loc = if_present(M, descA_loc(3))
          N_loc = if_present(N, descB_loc(4))
          K_loc = if_present(K, descA_loc(4))
       else
          M_loc = if_present(M, descA_loc(4))
          N_loc = if_present(N, descB_loc(3))
          K_loc = if_present(K, descA_loc(3))
       end if
    else
       M_loc = if_present(M, n_basis)
       N_loc = if_present(N, n_basis)
       K_loc = if_present(K, n_basis)
    end if
    iA_loc = if_present(iA, 1)
    jA_loc = if_present(jA, 1)
    iB_loc = if_present(iB, 1)
    jB_loc = if_present(jB, 1)
    iC_loc = if_present(iC, 1)
    jC_loc = if_present(jC, 1)
    if (do_parallel) then
       call pdgemm(transa_loc, transb_loc, M_loc, N_loc, K_loc, alpha_loc, A, &
            & iA_loc, jA_loc, descA_loc, B, iB_loc, jB_loc, descB_loc, &
            & beta_loc, C, iC_loc, jC_loc, descC_loc)
    else
       call dgemm(transa_loc, transb_loc, M_loc, N_loc, K_loc, alpha_loc, &
            & A(iA_loc,jA_loc), size(A,1), B(iB_loc,jB_loc), size(B,1), &
            & beta_loc, C(iC_loc,jC_loc), size(C,1))
    end if
  end subroutine pdgemm_wrapper

  pure subroutine pzgemm_wrapper(A, B, C, transa, transb, M, N, K, alpha, &
       & iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
    interface
       pure subroutine pzgemm(transa, transb, M, N, K, alpha, A, ia, ja, &
            & desca, B, ib, jb, descb, beta, C, ic, jc, descc)
         import dp
         complex(dp), intent(in) :: A(*), B(*), alpha, beta
         complex(dp), intent(in out) :: C(*)
         character, intent(in) :: transa, transb
         integer, intent(in) :: descA(9), descB(9), descC(9)
         integer, intent(in) :: M, N, K, ia, ja, ib, jb, ic, jc
       end subroutine pzgemm
       pure subroutine zgemm(transa,transb,M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
         import dp
         complex(dp), intent(in) :: A, B, alpha, beta
         complex(dp), intent(in out) :: C
         character, intent(in) :: transa, transb
         integer, intent(in) :: M, N, K, LDA, LDB, LDC
       end subroutine zgemm
    end interface
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(in out) :: C(:,:)
    character, intent(in), optional :: transa, transb
    complex(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: M, N, K, iA, jA, iB, jB, iC, jC
    integer, intent(in), optional :: descA(9), descB(9), descC(9)
    logical, intent(in), optional :: do_serial
    character :: transa_loc, transb_loc
    complex(dp) :: alpha_loc, beta_loc
    integer :: M_loc, N_loc, K_loc, iA_loc, jA_loc, descA_loc(9), iB_loc, &
         & jB_loc, descB_loc(9), iC_loc, jC_loc, descC_loc(9)
    logical :: do_parallel
    do_parallel = use_scalapack .and. .not. if_present(do_serial, .false.)
    alpha_loc = if_present(alpha, (1d0,0d0))
    beta_loc = if_present(beta, (0d0,0d0))
    transa_loc = if_present(transa, 'n')
    transb_loc = if_present(transb, 'n')
    descA_loc = if_present(descA, sc_desc)
    descB_loc = if_present(descB, sc_desc)
    descC_loc = if_present(descC, sc_desc)
    if (do_parallel) then
       if (transa_loc == 'n' .or. transa_loc == 'N') then
          M_loc = if_present(M, descA_loc(3))
          N_loc = if_present(N, descB_loc(4))
          K_loc = if_present(K, descA_loc(4))
       else
          M_loc = if_present(M, descA_loc(4))
          N_loc = if_present(N, descB_loc(3))
          K_loc = if_present(K, descA_loc(3))
       end if
    else
       M_loc = if_present(M, n_basis)
       N_loc = if_present(N, n_basis)
       K_loc = if_present(K, n_basis)
    end if
    iA_loc = if_present(iA, 1)
    jA_loc = if_present(jA, 1)
    iB_loc = if_present(iB, 1)
    jB_loc = if_present(jB, 1)
    iC_loc = if_present(iC, 1)
    jC_loc = if_present(jC, 1)
    if (do_parallel) then
       call pzgemm(transa_loc, transb_loc, M_loc, N_loc, K_loc, alpha_loc, A, &
            & iA_loc, jA_loc, descA_loc, B, iB_loc, jB_loc, descB_loc, &
            & beta_loc, C, iC_loc, jC_loc, descC_loc)
    else
       call zgemm(transa_loc, transb_loc, M_loc, N_loc, K_loc, alpha_loc, &
            & A(iA_loc,jA_loc), size(A,1), B(iB_loc,jB_loc), size(B,1), &
            & beta_loc, C(iC_loc,jC_loc), size(C,1))
    end if
  end subroutine pzgemm_wrapper

  pure subroutine pzgemm_im_re_wrapper(A, B, C, transa, transb, M, N, K, &
       & alpha, iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
    complex(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: B(:,:)
    complex(dp), intent(in out) :: C(:,:)
    character, intent(in), optional :: transa, transb
    complex(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: M, N, K, iA, jA, iB, jB, iC, jC
    integer, intent(in), optional :: descA(9), descB(9), descC(9)
    logical, intent(in), optional :: do_serial
    call pzgemm_wrapper(A, cmplx(B,0d0,8), C, transa, transb, M, N, K, alpha, &
       & iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
  end subroutine pzgemm_im_re_wrapper

  pure subroutine pzgemm_re_im_wrapper(A, B, C, transa, transb, M, N, K, &
       & alpha, iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
    real(dp), intent(in) :: A(:,:)
    complex(dp), intent(in) :: B(:,:)
    complex(dp), intent(in out) :: C(:,:)
    character, intent(in), optional :: transa, transb
    complex(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: M, N, K, iA, jA, iB, jB, iC, jC
    integer, intent(in), optional :: descA(9), descB(9), descC(9)
    logical, intent(in), optional :: do_serial
    call pzgemm_wrapper(cmplx(A,0d0,8), B, C, transa, transb, M, N, K, alpha, &
       & iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
  end subroutine pzgemm_re_im_wrapper

  pure subroutine pzgemm_re_re_im_wrapper(A, B, C, transa, transb, M, N, K, &
       & alpha, iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, do_serial)
    real(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(in out) :: C(:,:)
    character, intent(in), optional :: transa, transb
    complex(dp), intent(in), optional :: alpha, beta
    integer, intent(in), optional :: M, N, K, iA, jA, iB, jB, iC, jC
    integer, intent(in), optional :: descA(9), descB(9), descC(9)
    logical, intent(in), optional :: do_serial
    call pzgemm_wrapper(cmplx(A,0d0,8), cmplx(B,0d0,8), C, transa, transb, M, &
         & N, K, alpha, iA, jA, descA, iB, jB, descB, beta, iC, jC, descC, &
         & do_serial)
  end subroutine pzgemm_re_re_im_wrapper

  pure integer function if_present_int(A, B) result(Y)
    integer, intent(in), optional :: A
    integer, intent(in) :: B
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_int

  pure function if_present_int_1d(A, B) result(Y)
    integer, intent(in), optional :: A(:)
    integer, intent(in) :: B(:)
    integer :: Y(size(B))
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_int_1d

  pure logical function if_present_bool(A, B) result(Y)
    logical, intent(in), optional :: A
    logical, intent(in) :: B
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_bool

  pure character function if_present_char(A, B) result(Y)
    character, intent(in), optional :: A
    character, intent(in) :: B
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_char

  pure real(dp) function if_present_real(A, B) result(Y)
    real(dp), intent(in), optional :: A
    real(dp), intent(in) :: B
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_real

  pure complex(dp) function if_present_complex(A, B) result(Y)
    complex(dp), intent(in), optional :: A
    complex(dp), intent(in) :: B
    if (present(A)) then
       Y = A
    else
       Y = B
    end if
  end function if_present_complex

  !!  FUNCTION
  !!
  !!  Reimplementation of the norm2 intrinsic (F2008).
  !!
  pure real(dp) function norm2_vector(x) result(y)
    real(dp), intent(in) :: x(:)
    y = sqrt(sum(x**2))
  end function norm2_vector

  pure function norm2_matrix(x, dim) result(y)
    real(dp), intent(in) :: x(:,:)
    ! Note that here dim can only be 1, so this is not exactly the
    ! same as the intrinsic function.
    integer, intent(in) :: dim
    real(dp) :: y(size(x,2))
    integer :: i
    if (dim == 1) y = [(norm2_vector(x(:,i)), i=1,size(x,2))]
  end function norm2_matrix

  !!  FUNCTION
  !!
  !!  Transfer elements from a local matrix corresponding to a given
  !!  batch into a global packed array.
  !!
  pure subroutine batch_to_matrix_packed(n_compute, i_basis, matrix_batch, &
       & matrix_packed)
    ! n_compute - number of nonzero basis functions in a given batch
    integer, intent(in) :: n_compute, i_basis(n_compute)
    real(dp), intent(in) :: matrix_batch(n_compute,n_compute)
    real(dp), intent(in out) :: matrix_packed(:,:)
    integer :: i_range(2), i_new_start ! For iterating over index_hamiltonian
    integer :: mb, nb ! Reduced indices for basis functions
    integer :: i_index, i_off
    if (.not. use_local_index) then
       ! No packing here
       matrix_packed(i_basis,i_basis) = matrix_packed(i_basis,i_basis) + &
            & matrix_batch
       return
    end if
    if (use_load_balancing) then
       if (size(matrix_packed,2) == 1) then
          do mb = 1, n_compute
             ! Local/global indices are mapped according to
             ! batch_perm(n_bp_integ)%i_basis_glb_to_loc. This scheme
             ! was first used in integrate_real_hamiltonian_matrix_p2.
             i_off = batch_perm(n_bp_integ)%i_basis_glb_to_loc(i_basis(mb))
             i_off = (i_off*(i_off-1))/2
             do nb = 1, mb
                i_index = batch_perm(n_bp_integ)% &
                     & i_basis_glb_to_loc(i_basis(nb)) + i_off
                matrix_packed(i_index,1) = &
                     & matrix_packed(i_index,1) + matrix_batch(mb,nb)
             end do
          end do
       else
          ! If the full matrix is to be stored, matrix_packed(:,1)
          ! contains the lower triangle and matrix_packed(:,2) is the
          ! upper triangle. Diagonal elements are counted twice here
          ! and will later be divided by two.
          do mb = 1, n_compute
             i_off = batch_perm(n_bp_integ)%i_basis_glb_to_loc(i_basis(mb))
             i_off = (i_off*(i_off-1))/2
             do nb = 1, mb
                i_index = batch_perm(n_bp_integ)% &
                     & i_basis_glb_to_loc(i_basis(nb)) + i_off
                matrix_packed(i_index,1) = &
                     & matrix_packed(i_index,1) + matrix_batch(mb,nb)
                matrix_packed(i_index,2) = &
                     & matrix_packed(i_index,2) + matrix_batch(nb,mb)
             end do
          end do
       end if
    else
       if (size(matrix_packed,2) == 1) then
          do mb = 1, n_compute
             ! The range of column indices corresponding to
             ! i_basis(mb) as the row element
             i_range = index_hamiltonian(:, 1, i_basis(mb))
             i_new_start = 1
             do nb = 1, mb ! We work on the lower triangle
                do i_index = i_range(1), i_range(2)
                   ! See if any element in the indexing array
                   ! corresponds to the current i_basis(mb) and
                   ! i_basis(nb).
                   if (column_index_hamiltonian(i_index) == i_basis(nb)) then
                      matrix_packed(i_index,1) = matrix_packed(i_index,1) + &
                           & matrix_batch(mb,nb)
                      i_new_start = i_index + 1
                      exit
                   else if (column_index_hamiltonian(i_index) > i_basis(nb)) &
                        & then
                      ! The remaining indices in the indexing array
                      ! are necessarily larger than the given colum
                      ! index. There is no reason to continue this
                      ! loop.
                      i_new_start = i_index
                      exit
                   end if
                end do
                ! Given the ordering of indices in the indexing array,
                ! the next suitable column index cannot be found at a
                ! location that is less than that of the previously
                ! found element. We may thus start the search at
                ! i_new_start without losing any information.
                i_range(1) = i_new_start
             end do
          end do
       else
          do mb = 1, n_compute
             i_range = index_hamiltonian(:, 1, i_basis(mb))
             i_new_start = 1
             do nb = 1, mb
                do i_index = i_range(1), i_range(2)
                   if (column_index_hamiltonian(i_index) == i_basis(nb)) then
                      matrix_packed(i_index,1) = matrix_packed(i_index,1) + &
                           & matrix_batch(mb,nb)
                      matrix_packed(i_index,2) = matrix_packed(i_index,2) + &
                           & matrix_batch(nb,mb)
                      i_new_start = i_index + 1
                      exit
                   else if (column_index_hamiltonian(i_index) > i_basis(nb)) &
                        & then
                      i_new_start = i_index
                      exit
                   end if
                end do
                i_range(1) = i_new_start
             end do
          end do
       end if
    end if
  end subroutine batch_to_matrix_packed

  !!  FUNCTION
  !!
  !!  Unpacks the contents of matrix_packed into
  !!  matrix_batch(:n_compute,:n_compute). It is basically
  !!  batch_to_matrix_packed in reverse, except here it is assumed
  !!  that only the lower triangle is ever saved in the packed array.
  !!
  pure subroutine matrix_packed_to_batch(n_compute, i_basis, matrix_packed, &
       & matrix_batch)
    integer, intent(in) :: n_compute, i_basis(n_compute)
    real(dp), intent(in) :: matrix_packed(*)
    real(dp), intent(out) :: matrix_batch(n_compute,n_compute)
    integer :: mb, nb, i_range(2), i_new_start, i_index, i_off
    if (.not. use_local_index) then
       ! No packing here
       do nb = 1, n_compute
          matrix_batch(:,nb) = matrix_packed(i_basis+(i_basis(nb)-1)*n_basis)
       end do
       return
    end if
    if (use_load_balancing) then
       do mb = 1, n_compute
          ! See the comments in batch_to_matrix_packed
          i_off = batch_perm(n_bp_integ)%i_basis_glb_to_loc(i_basis(mb))
          i_off = (i_off*(i_off-1))/2
          do nb = 1, mb
             i_index = batch_perm(n_bp_integ)% &
                  & i_basis_glb_to_loc(i_basis(nb)) + i_off
             matrix_batch(mb,nb) = matrix_packed(i_index)
             matrix_batch(nb,mb) = matrix_packed(i_index)
          end do
       end do
    else
       do mb = 1, n_compute
          i_range = index_hamiltonian(:, 1, i_basis(mb))
          i_new_start = 1
          do nb = 1, mb
             do i_index = i_range(1), i_range(2)
                if (column_index_hamiltonian(i_index) == i_basis(nb)) then
                   matrix_batch(mb,nb) = matrix_packed(i_index)
                   matrix_batch(nb,mb) = matrix_packed(i_index)
                   i_new_start = i_index + 1
                   exit
                else if (column_index_hamiltonian(i_index) > i_basis(nb)) then
                   i_new_start = i_index
                   exit
                end if
             end do
             i_range(1) = i_new_start
          end do
       end do
    end if
  end subroutine matrix_packed_to_batch

  !!  FUNCTION
  !!
  !!  Procedures for converting an array that is stored in the packed
  !!  format into a 2D block cyclic array or vice
  !!  versa. matrix_packed(:,1) contains the upper triangle while
  !!  matrix_packed(:,2) (only if size(matrix_packed,2)==2) contains
  !!  the lower triangle.
  !!
  subroutine packed_to_block_cyclic(matrix_packed, matrix_BC)
    real(dp), intent(in) :: matrix_packed(:,:)
    real(dp), intent(out) :: matrix_BC(:,:)
    integer :: mb, n_spin_save
    character*20 :: output_level_save
    if (.not. use_local_index) then
       ! No packing here
       matrix_BC = matrix_packed(my_row,my_col)
       return
    end if
    call start_wall(walltime_mat_conv)
    output_level_save = output_level
    output_level = 'MD_light'
    n_spin_save = n_spin
    n_spin = 1
    if (use_load_balancing) then
       call set_full_local_ham(matrix_packed)
    else
       call set_sparse_local_ham_scalapack(matrix_packed)
    end if
    if (size(matrix_packed,2) == 1) then
       call set_full_matrix_real(ham(:,:,1))
       matrix_BC = ham(:,:,1)
    else
       matrix_BC = ham(:,:,1)
       call pdtran(n_basis, n_basis, 1d0, matrix_BC, 1, 1, sc_desc, 0d0, ham, &
            & 1, 1, sc_desc)
       matrix_BC = ham(:,:,1)
       if (use_load_balancing) then
          call set_full_local_ham(matrix_packed(:,2))
       else
          call set_sparse_local_ham_scalapack(matrix_packed(:,2))
       end if
       do mb = 1, n_basis
          if (l_row(mb) > 0 .and. l_col(mb) > 0) &
               & ham(l_row(mb),l_col(mb),1) = 0d0
       end do
       matrix_BC = matrix_BC + ham(:,:,1)
    end if
    n_spin = n_spin_save
    output_level = output_level_save
    call stop_wall(walltime_mat_conv)
  end subroutine packed_to_block_cyclic

  subroutine block_cyclic_to_packed(matrix_BC, matrix_packed)
    real(dp), intent(in) :: matrix_BC(:,:)
    real(dp), intent(out) :: matrix_packed(:,:)
    character*20 :: output_level_save
    output_level_save = output_level
    output_level = 'MD_light'
    if (.not. use_local_index) then
       ! No packing here
       matrix_packed = 0d0
       matrix_packed(my_row,my_col) = matrix_BC
       if (use_scalapack) call sync_vector(matrix_packed, size(matrix_packed))
       return
    end if
    call start_wall(walltime_mat_conv)
    ham(:,:,1) = matrix_BC
    if (use_load_balancing) then
       call get_full_local_matrix(matrix_packed,1)
    else
       call get_sparse_local_matrix_scalapack(matrix_packed,1)
    end if
    output_level = output_level_save
    call stop_wall(walltime_mat_conv)
  end subroutine block_cyclic_to_packed

  subroutine start_wall(walltime)
    integer, intent(in out) :: walltime
    integer :: time_tmp, i_err
    call mpi_barrier(mpi_comm_global, i_err)
    call system_clock(time_tmp)
    walltime = walltime - time_tmp
  end subroutine start_wall

  subroutine stop_wall(walltime)
    integer, intent(in out) :: walltime
    integer :: time_tmp
    call system_clock(time_tmp)
    walltime = walltime + time_tmp
  end subroutine stop_wall

  subroutine safe_deallocate_r1(arr)
    real(dp), allocatable, intent(in out) :: arr(:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_r1

  subroutine safe_deallocate_r2(arr)
    real(dp), allocatable, intent(in out) :: arr(:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_r2

  subroutine safe_deallocate_r3(arr)
    real(dp), allocatable, intent(in out) :: arr(:,:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_r3

  subroutine safe_deallocate_r4(arr)
    real(dp), allocatable, intent(in out) :: arr(:,:,:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_r4

  subroutine safe_deallocate_r5(arr)
    real(dp), allocatable, intent(in out) :: arr(:,:,:,:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_r5

  subroutine safe_deallocate_i1(arr)
    integer, allocatable, intent(in out) :: arr(:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_i1

  subroutine safe_deallocate_i2(arr)
    integer, allocatable, intent(in out) :: arr(:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_i2

  subroutine safe_deallocate_i3(arr)
    integer, allocatable, intent(in out) :: arr(:,:,:)
    if (allocated(arr)) call aims_deallocate(arr, ' ')
  end subroutine safe_deallocate_i3
end module tools
