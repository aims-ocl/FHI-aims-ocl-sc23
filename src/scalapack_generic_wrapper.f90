! Over time, scalapack_wrapper has become 12000+ lines of subroutines, many of
! which are copy-pastes of one another, referring to one special case:  matrices
! with an identical structure to hamiltonian including a possible
! real-space-to-Bloch-representation conversion.
!
! That is all well and good, but it is inapplicable to the simple case of
! "I have a generic matrix at this k-point, the leading dimension is not
! n_basis, and I want to use it with ScaLAPACK/ELPA/ELSI."
!
! This collection of subroutines (not a module) is intended to simplify the set-up
! of a ScaLAPACK distribution of a general matrix that is defined on top of aims'
! pre-existing k-point distribution/MPI communicators.  This is a special case, but
! Bloch's theorem guarantees that this is the most common special case.
!
! Assumptions:
! o The matrix, in the periodic case, has a diagonal k-point dependence; that
!   is, no k-point coupling. (In the non-periodic case, this assumption is trivial.)
! o We are using the same k-point distribution across ranks as the Hamiltonian,
!   with the associated MPI communicator for the current rank being my_scalapack_comm_all.
! o We are using the same BLACS grid/context (my_blacs_ctxt) as the Hamiltonian.
! o The indexing functions support rectangular matrices, but the ScaLAPACK descriptor
!   initalization requires square matrices.  Adding support for rectangular matrices
!   to the later should be straightforward, and I (WPH) have constructed the interface
!   to support this in the future.  I just don't have the time to debug that functionality
!   right now.
!
!-----------------------------------------------------------------------------------
!****s* initialize_scalapack_descriptor_matrix
! NAME
!    initialize_scalapack_descriptor_matrix
  subroutine initialize_scalapack_descriptor_matrix( n_global_rows, n_global_cols, block_size, sc_desc, &
                  mb, nb, mxld, mxcol, n_my_rows, n_my_cols, l_row, l_col, my_row, my_col )
!  PURPOSE
!    Create a 2D block cyclic distribution for a matrix, as well as the indexing
!    arrays between the local matrix and the global matrix.
!
!    This function is a fork of the middle section of initialize_scalapack,
!    where the ScaLAPACK descriptor and various indexing variables are created
!    for the Hamiltonian and similarly-packed matrices.
! USES
  use runtime_choices, only : use_elpa, scalapack_block_size
  use localorb_io, only : localorb_info, use_unit
  use synchronize_mpi_basic, only : sync_find_max
  use mpi_tasks, only : aims_stop, myid, mpi_comm_global
  use scalapack_wrapper, only : my_blacs_ctxt, dlen_, my_scalapack_id, nprow, npcol, myprow, mypcol
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: n_global_rows
  integer, intent(in)  :: n_global_cols
! OUTPUTS
  integer, intent(out) :: block_size
  integer, intent(out), dimension(dlen_) :: sc_desc
  integer, intent(out) :: mb
  integer, intent(out) :: nb
  integer, intent(out) :: mxld
  integer, intent(out) :: mxcol
  integer, intent(out) :: n_my_rows
  integer, intent(out) :: n_my_cols
  integer, intent(out), dimension(n_global_rows) :: l_row
  integer, intent(out), dimension(n_global_cols) :: l_col
  integer, intent(out), dimension(n_global_rows) :: my_row ! n_my_rows elements will be valid
  integer, intent(out), dimension(n_global_cols) :: my_col ! n_my_rows elements will be valid
! AUTHOR
!   William Huhn
! HISTORY
!   November 2016, fork off of initialize_scalapack
! SOURCE
  integer, parameter :: csrc = 0
  integer, parameter :: rsrc = 0

  integer :: max_npcol, max_nprow, matrix_dim
  integer :: info
  integer :: mpierr
  integer :: i, lc, lr

  integer, external :: numroc

  character*200 :: info_str

  if ( n_global_rows .ne. n_global_cols ) then
    call aims_stop("initialize_scalapack_descriptor_matrix does not currently support non-rectangular matrices, exiting.")
  end if
  matrix_dim = n_global_rows

  mb = 64
  nb = 64

  ! My guess is that, should one want to generalize this code to rectangular
  ! matrices, that you use the smaller of n_global_rows and n_global_cols to
  ! determine block_size
  ! Do not use different block sizes for the rows and columns!  That is, mb
  ! should not be different from nb.  Some libraries (notably ELPA) do not
  ! support this.

  ! Calculate blocksize (if the user does not specify) based on matrix_dim
  ! and nprow/npcol
  ! We want that every processor owns a part of the matrix,
  ! i.e. mxld>0 and mxcol>0 everywhere.
  ! This is needed to catch a "bug" in ScaLAPACK:
  ! If a processor doesn't own a part of a matrix, the results are not
  ! distributed to this one
  ! Theoretically, nprow/npcol can be different for different k-points,
  ! so we have to get the global maximum.

  call sync_find_max(nprow,max_nprow)
  call sync_find_max(npcol,max_npcol)

  block_size = 1 ! Minimum permitted size
  write(info_str, *) ' Calculating blocksize based on matrix dimension = ',matrix_dim, &
                       ' max_nprow = ',max_nprow,' max_npcol = ',max_npcol
  call localorb_info(info_str, use_unit, '(A)')

  if(block_size*MAX(max_nprow,max_npcol) > matrix_dim) then
    write(info_str, *) 'ERROR: matrix dimension = ',matrix_dim,' too small for this processor grid'
    call localorb_info(info_str, use_unit, '(A)')
    call aims_stop
  endif

  ! Increase blocksize to maximum possible size or 64
  do while (2*block_size*MAX(max_nprow,max_npcol) <= matrix_dim &
       .and. block_size<64)
    block_size = 2*block_size
  end do

  ! use_elpa works best with SMALL blocksizes!
  if(use_elpa .and. block_size>16) block_size = 16

  ! If block_size is defined in control.in file, check it here.
  if(scalapack_block_size > 0) then ! Set by user.
    if(scalapack_block_size*(MAX(max_nprow,max_npcol)-1) .ge. matrix_dim) then
      ! See comments above!
      write(info_str, *) &
           'ERROR: User defined block size ',scalapack_block_size,&
           ' too large for this process grid and the matrix dimension.'
      call localorb_info(info_str, use_unit, '(A)')
      call aims_stop
    else
      ! (Hopefully) reasonable value is accepted here.
      block_size = scalapack_block_size
      write(info_str, *) &
           ' Use block size defined in control.in file.'
      call localorb_info(info_str, use_unit, '(A)')
    endif
  endif

  nb = block_size
  mb = block_size
  write(info_str, *) ' Scalapack blocksize set to: ',block_size
  call localorb_info(info_str, use_unit, '(A)')

  ! initialize the Scalapack descriptor

  if(my_scalapack_id<npcol*nprow) then
     mxld = numroc( matrix_dim, mb, myprow, rsrc, nprow )
     mxcol = numroc( matrix_dim, nb, mypcol, csrc, npcol )

! RJ: If mxld/mxcol are too small, they *might* trigger an error in the
! Intel/Scalapack implementation, so set them to a at least 64:

! BL: This is potentially dangerous in case the local
!     dimensions are smaller than matrix_dim: Horror for parallel writing based on
!     patterns

!       if(mxld  < 64) mxld  = 64
!       if(mxcol < 64) mxcol = 64

    call descinit( sc_desc, matrix_dim, matrix_dim, mb, nb, rsrc, csrc, &
         my_blacs_ctxt, MAX(1,mxld), info )

    ! Safety check only, the following should never happen
    if(mxld<=0 .or. mxcol<=0) then
      print *,'ERROR Task #',myid,' mxld= ',mxld,' mxcol= ',mxcol
      call mpi_abort(mpi_comm_global,1,mpierr)
    endif

  else
    mxld = 1
    mxcol = 1
  endif

  ! Mapping of global rows/cols to local

  l_row(:) = 0
  l_col(:) = 0

  ! ATTENTION: The following code assumes rsrc==0 and csrc==0 !!!!
  ! For processors outside the working set, l_row/l_col will stay
  ! completely at 0

  lr = 0 ! local row counter
  lc = 0 ! local column counter

  do i = 1, matrix_dim

    if( MOD((i-1)/mb,nprow) == myprow) then
      ! row i is on local processor
      lr = lr+1
      l_row(i) = lr
    endif

    if( MOD((i-1)/nb,npcol) == mypcol) then
      ! column i is on local processor
      lc = lc+1
      l_col(i) = lc
    endif

  enddo

  ! Mapping of local rows/cols to global

  n_my_rows = lr
  n_my_cols = lc
  lr = 0
  lc = 0

  do i = 1, matrix_dim
    if(l_row(i)>0) then; lr = lr+1; my_row(lr) = i; endif
    if(l_col(i)>0) then; lc = lc+1; my_col(lc) = i; endif
  enddo
end subroutine initialize_scalapack_descriptor_matrix

!-----------------------------------------------------------------------------------
!****s* convert_global_to_local_matrix
! NAME
!    convert_global_to_local_matrix
! SYNOPSIS
subroutine convert_global_to_local_matrix( n_global_rows, n_global_cols, l_row, l_col, global_matrix, &
                mxld, mxcol, local_matrix )
!  PURPOSE
!    This subroutine fills in the elements of a real local (2D block-cyclic) matrix from
!    the global matrix.  Ideally, this subroutine shouldn't be needed; one should work with
!    the local matrix only.  But when working with legacy code, this may be necessary to
!    avoid a substantial code rewrite.
!
!    The reverse operation is done by convert_local_to_global_matrix.
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: l_row(n_global_rows)
  integer, intent(in)  :: l_col(n_global_cols)
  integer, intent(in)  :: n_global_rows
  integer, intent(in)  :: n_global_cols
  real*8,  intent(in)  :: global_matrix( n_global_rows,n_global_cols )
  integer, intent(in)  :: mxld
  integer, intent(in)  :: mxcol
! OUTPUTS
  real*8,  intent(out) :: local_matrix( mxld, mxcol )
! AUTHOR
!   William Huhn
! HISTORY
!   November 2016, modified version of setup_scalapack_rmatrix
! SOURCE
  integer :: i_col, i_row, lr, lc

  local_matrix = 0.0d0

  do i_col = 1, n_global_cols, 1
    lc = l_col(i_col) ! local column number
    if(lc>0) then
      do i_row = 1, n_global_rows, 1
        lr = l_row(i_row) ! local row number
        if(lr>0) then
          local_matrix(lr,lc) = global_matrix(i_row, i_col)
        endif
      enddo
    endif
  end do

end subroutine convert_global_to_local_matrix

!-----------------------------------------------------------------------------------
!****s* convert_global_to_local_matrix_complex
! NAME
!    convert_global_to_local_matrix_complex
! SYNOPSIS
subroutine convert_global_to_local_matrix_complex( n_global_rows, n_global_cols, l_row, l_col, global_matrix, &
                mxld, mxcol, local_matrix )
!  PURPOSE
!    This subroutine fills in the elements of a complex local (2D block-cyclic) matrix from
!    the global matrix.  Ideally, this subroutine shouldn't be needed; one should work with
!    the local matrix only.  But when working with legacy code, this may be necessary to
!    avoid a substantial code rewrite.
!
!    The reverse operation is done by convert_local_to_global_matrix_complex.
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: n_global_rows
  integer, intent(in)  :: n_global_cols
  integer, intent(in)  :: l_row(n_global_rows)
  integer, intent(in)  :: l_col(n_global_cols)
  complex*16,  intent(in)  :: global_matrix( n_global_rows,n_global_cols )
  integer, intent(in)  :: mxld
  integer, intent(in)  :: mxcol
! OUTPUTS
  complex*16,  intent(out) :: local_matrix( mxld, mxcol )
! AUTHOR
!   William Huhn
! HISTORY
!   November 2016, modified version of setup_scalapack_rmatrix
! SOURCE
  integer :: i_col, i_row, lr, lc

  local_matrix = (0.0d0, 0.0d0)

  do i_col = 1, n_global_cols, 1
    lc = l_col(i_col) ! local column number
    if(lc>0) then
      do i_row = 1, n_global_rows, 1
        lr = l_row(i_row) ! local row number
        if(lr>0) then
          local_matrix(lr,lc) = global_matrix(i_row, i_col)
        endif
      enddo
    endif
  end do

end subroutine convert_global_to_local_matrix_complex

!-----------------------------------------------------------------------------------
!****s* convert_local_to_global_matrix
! NAME
!    convert_local_to_global_matrix
! SYNOPSIS
subroutine convert_local_to_global_matrix( n_my_rows, n_my_cols, my_row, my_col, &
                mxld, mxcol, local_matrix, n_global_rows, n_global_cols, global_matrix )
!  PURPOSE
!    This subroutine fills in the elements for a real global matrix from a local
!    (2D block-cyclic) matrix.  Ideally, this subroutine shouldn't be needed; one should
!    work with the local matrix only.  But when working with legacy code, this
!    may be necessary to avoid a substantial code rewrite.
!
!    The reverse operation is done by convert_global_to_local_matrix_complex.
!
!    Note:  unlike the original subroutine, we do *not* sync the matrix afterwards.  I can
!    envision debug cases where one would not want to do this.  Since MPI calls are a
!    significant bottleneck that should be avoided as much as possible, it's good practice
!    to force the programmer to do it themselves.
!
!    "Only the upper half of the matrices is set in the code below [WPH: unless you set the
!     lower half in the local matrix yourself], so the parameter UPLO to the Scalapack routines
!     must always be 'U'!"
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: n_my_rows
  integer, intent(in)  :: n_my_cols
  integer, intent(in)  :: my_row(n_my_rows)
  integer, intent(in)  :: my_col(n_my_cols)
  integer, intent(in)  :: mxld
  integer, intent(in)  :: mxcol
  real*8,  intent(in)  :: local_matrix( mxld, mxcol )
  integer, intent(in)  :: n_global_rows
  integer, intent(in)  :: n_global_cols
! OUTPUTS
  real*8,  intent(out)  :: global_matrix( n_global_rows, n_global_cols )
! AUTHOR
!   William Huhn
! HISTORY
!   November 2016
! SOURCE
  integer :: i_col, i_row

  global_matrix = 0.0d0

  do i_row = 1, n_my_rows
    do i_col = 1, n_my_cols
      ! The last conditional is needed, as often we will use a ScaLAPACK
      ! descriptor with a matrix with a certain leading dimension but fewer
      ! columns than the original matrix specified by the ScaLAPACK descriptor
      ! (for example, using sc_desc for a given spin channel for both hamiltonian
      ! and KS_eigenvector, even though the former is an n_basis*n_basis matrix
      ! and the later is an n_basis*n_states matrix)
      if( my_row(i_row) > 0 .and. my_col(i_col) > 0 .and. my_col(i_col) .le. n_global_cols ) then
        global_matrix( my_row(i_row), my_col(i_col) ) = local_matrix(i_row,i_col)
      endif
    end do
  end do

!  call sync_complex_vector(global_matrix, matrix_dim*matrix_dim, my_scalapack_comm_all)
end subroutine convert_local_to_global_matrix

!-----------------------------------------------------------------------------------
!****s* convert_local_to_global_matrix_complex
! NAME
!    convert_local_to_global_matrix_complex
! SYNOPSIS
subroutine convert_local_to_global_matrix_complex( n_my_rows, n_my_cols, my_row, my_col, &
                mxld, mxcol, local_matrix, n_global_rows, n_global_cols, global_matrix )
!  PURPOSE
!    This subroutine fills in the elements for a complex global matrix from a local
!    (2D block-cyclic) matrix.  Ideally, this subroutine shouldn't be needed; one should
!    work with the local matrix only.  But when working with legacy code, this
!    may be necessary to avoid a substantial code rewrite.
!
!    The reverse operation is done by convert_global_to_local_matrix_complex.
!
!    Note:  unlike the original subroutine, we do *not* sync the matrix afterwards.  I can
!    envision debug cases where one would not want to do this.  Since MPI calls are a
!    significant bottleneck that should be avoided as much as possible, it's good practice
!    to force the programmer to do it themselves.
!
!    "Only the upper half of the matrices is set in the code below [WPH: unless you set the
!     lower half in the local matrix yourself], so the parameter UPLO to the Scalapack routines
!     must always be 'U'!"
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: n_my_rows
  integer, intent(in)  :: n_my_cols
  integer, intent(in)  :: my_row(n_my_rows)
  integer, intent(in)  :: my_col(n_my_cols)
  integer, intent(in)  :: mxld
  integer, intent(in)  :: mxcol
  complex*16,  intent(in)  :: local_matrix( mxld, mxcol )
  integer, intent(in)  :: n_global_rows
  integer, intent(in)  :: n_global_cols
! OUTPUTS
  complex*16,  intent(out)  :: global_matrix( n_global_rows, n_global_cols )
! AUTHOR
!   William Huhn
! HISTORY
!   November 2016
! SOURCE
  integer :: i_col, i_row

  global_matrix = (0.0d0, 0.0d0)

  do i_row = 1, n_my_rows
    do i_col = 1, n_my_cols
      ! The last conditional is needed, as often we will use a ScaLAPACK
      ! descriptor with a matrix with a certain leading dimension but fewer
      ! columns than the original matrix specified by the ScaLAPACK descriptor
      ! (for example, using sc_desc for a given spin channel for both hamiltonian
      ! and KS_eigenvector, even though the former is an n_basis*n_basis matrix
      ! and the later is an n_basis*n_states matrix)
      if( my_row(i_row) > 0 .and. my_col(i_col) > 0 .and. my_col(i_col) .le. n_global_cols ) then
        global_matrix( my_row(i_row), my_col(i_col) ) = local_matrix(i_row,i_col)
      endif
    end do
  end do

!  call sync_complex_vector(global_matrix, matrix_dim*matrix_dim, my_scalapack_comm_all)
end subroutine convert_local_to_global_matrix_complex

!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_square_matrix
!  NAME
!    set_full_square_matrix
!  SYNOPSIS
  subroutine set_full_square_matrix(matrix_dim, nb, sc_desc, l_row, l_col, mxld, mxcol, mat )
!  PURPOSE
!    For square matrix, sets the lower half of a distributed matrix from the upper half.
!  USES
    use scalapack_wrapper, only : dlen_
    implicit none
!  ARGUMENTS
!  INPUTS
    integer, intent(in)                               :: matrix_dim
    integer, intent(in)                               :: nb
    integer, intent(in), dimension(dlen_)             :: sc_desc
    integer, intent(in), dimension(matrix_dim)        :: l_row
    integer, intent(in), dimension(matrix_dim)        :: l_col
    integer, intent(in)                               :: mxld
    integer, intent(in)                               :: mxcol
!  OUTPUTS
    real*8,  intent(inout), dimension(mxld, mxcol)    :: mat
!  AUTHOR
!    William Huhn
!  HISTORY
!    November 2016, forked from set_full_matrix_real
!  SOURCE

    integer i_col, i_row
    real*8, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb)) ! no idea whats really needed

    ! This routine is called only from the working set, so no need to check here

    call pdtran(matrix_dim,matrix_dim,1.d0,mat,1,1,sc_desc,0.d0,tmp2,1,1,sc_desc)

    do i_col=1,matrix_dim-1
       if(l_col(i_col)==0) cycle
       do i_row=i_col+1,matrix_dim
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    deallocate(tmp2)

  end subroutine set_full_square_matrix


!******
!-----------------------------------------------------------------------------------
!****s* scalapack_wrapper/set_full_square_matrix_complex
!  NAME
!    set_full_square_matrix_complex
!  SYNOPSIS
  subroutine set_full_square_matrix_complex(matrix_dim, nb, sc_desc, l_row, l_col, mxld, mxcol, mat )
!  PURPOSE
!    For square matrix, sets the lower half of a distributed matrix from the upper half.
!  USES
    use scalapack_wrapper, only : dlen_
    implicit none
!  ARGUMENTS
!  INPUTS
    integer, intent(in)                               :: matrix_dim
    integer, intent(in)                               :: nb
    integer, intent(in), dimension(dlen_)             :: sc_desc
    integer, intent(in), dimension(matrix_dim)        :: l_row
    integer, intent(in), dimension(matrix_dim)        :: l_col
    integer, intent(in)                               :: mxld
    integer, intent(in)                               :: mxcol
!  OUTPUTS
    complex*16, intent(inout), dimension(mxld, mxcol) :: mat
!  AUTHOR
!    William Huhn
!  HISTORY
!    November 2016, forked from set_full_matrix_complex
!  SOURCE

    integer i_col, i_row
    complex*16, allocatable :: tmp2(:,:)

    ! Allocate tmp2 bigger than necessary to catch overwrites in pdtran

    allocate(tmp2(mxld,mxcol+2*nb)) ! no idea whats really needed

    ! This routine is called only from the working set, so no need to check here
    call pztranc(matrix_dim,matrix_dim,(1.d0,0.d0),mat,1,1,sc_desc,(0.d0,0.d0),tmp2,1,1,sc_desc)

    do i_col=1,matrix_dim-1
       if(l_col(i_col)==0) cycle
       do i_row=i_col+1,matrix_dim
          if(l_row(i_row)>0) mat(l_row(i_row),l_col(i_col)) = tmp2(l_row(i_row),l_col(i_col))
       enddo
    enddo

    ! For safety: Make diagonal real

    do i_col=1,matrix_dim
       if(l_col(i_col)==0 .or. l_row(i_col)==0) cycle
       mat(l_row(i_col),l_col(i_col)) = dble(mat(l_row(i_col),l_col(i_col)))
    enddo

    deallocate(tmp2)

  end subroutine set_full_square_matrix_complex

!-----------------------------------------------------------------------------------
!****s* convert_local_to_global_up_tri_matrix
! NAME
!    convert_local_to_global_up_tri_matrix
! SYNOPSIS
subroutine convert_local_to_global_up_tri_matrix( n_my_rows, n_my_cols, my_row, my_col, &
                mxld, mxcol, local_matrix, global_dim, global_matrix )
!  PURPOSE
!    This subroutine fills in the elements for a real (square) global matrix stored in a 1D
!    upper triangular format from a local (2D block-cyclic) matrix.  Ideally, this
!    subroutine shouldn't be needed; one should work with the local matrix only.  But when
!    working with legacy code, this may be necessary to avoid a substantial code rewrite.
!
!    Note:  unlike the original subroutine, we do *not* sync the matrix afterwards.  I can
!    envision debug cases where one would not want to do this.  Since MPI calls are a
!    significant bottleneck that should be avoided as much as possible, it's good practice
!    to force the programmer to do it themselves.
!
!    Assumes that the upper triangular half of the local matrix is set.
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer, intent(in)  :: n_my_rows
  integer, intent(in)  :: n_my_cols
  integer, intent(in)  :: my_row(n_my_rows)
  integer, intent(in)  :: my_col(n_my_cols)
  integer, intent(in)  :: mxld
  integer, intent(in)  :: mxcol
  real*8,  intent(in)  :: local_matrix( mxld, mxcol )
  integer, intent(in)  :: global_dim
! OUTPUTS
  real*8,  intent(out) :: global_matrix( global_dim * (global_dim+1)/2 )
! AUTHOR
!   William Huhn
! HISTORY
!   September 2017 - Created.
! SOURCE
  integer :: i_col, i_col_global, i_row, i_row_global, i_index

  global_matrix = 0.0d0

  do i_row = 1, n_my_rows
    do i_col = 1, n_my_cols
      ! The last conditional is needed, as often we will use a ScaLAPACK
      ! descriptor with a matrix with a certain leading dimension but fewer
      ! columns than the original matrix specified by the ScaLAPACK descriptor
      ! (for example, using sc_desc for a given spin channel for both hamiltonian
      ! and KS_eigenvector, even though the former is an n_basis*n_basis matrix
      ! and the later is an n_basis*n_states matrix)
      i_row_global = my_row(i_row)
      i_col_global = my_col(i_col)
      if( i_row_global > 0 .and. i_col_global > 0 .and. (i_row_global .le. i_col_global) .and. &
           (i_col_global .le. global_dim) ) then
        ! Offset for the current column
        i_index = ((i_col_global-1)*(i_col_global))/2
        ! Offset for the current row
        i_index = i_index + i_row_global

        global_matrix( i_index ) = local_matrix(i_row,i_col)
      endif
    end do
  end do

end subroutine convert_local_to_global_up_tri_matrix

!-----------------------------------------------------------------------------------
!****s* convert_local_to_global_up_tri_matrix_complex
! NAME
!    convert_local_to_global_up_tri_matrix_complex
! SYNOPSIS
subroutine convert_local_to_global_up_tri_matrix_complex( n_my_rows, n_my_cols, my_row, my_col, &
                mxld, mxcol, local_matrix, global_dim, global_matrix )
!  PURPOSE
!    This subroutine fills in the elements for a complex (square) global matrix stored in a 1D
!    upper triangular format from a local (2D block-cyclic) matrix.  Ideally, this
!    subroutine shouldn't be needed; one should work with the local matrix only.  But when
!    working with legacy code, this may be necessary to avoid a substantial code rewrite.
!
!    Note:  unlike the original subroutine, we do *not* sync the matrix afterwards.  I can
!    envision debug cases where one would not want to do this.  Since MPI calls are a
!    significant bottleneck that should be avoided as much as possible, it's good practice
!    to force the programmer to do it themselves.
!
!    Assumes that the upper triangular half of the local matrix is set.
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer,    intent(in)  :: n_my_rows
  integer,    intent(in)  :: n_my_cols
  integer,    intent(in)  :: my_row(n_my_rows)
  integer,    intent(in)  :: my_col(n_my_cols)
  integer,    intent(in)  :: mxld
  integer,    intent(in)  :: mxcol
  complex*16, intent(in)  :: local_matrix( mxld, mxcol )
  integer,    intent(in)  :: global_dim
! OUTPUTS
  complex*16, intent(out) :: global_matrix( global_dim * (global_dim+1)/2 )
! AUTHOR
!   William Huhn
! HISTORY
!   September 2017 - Created.
! SOURCE
  integer :: i_col, i_col_global, i_row, i_row_global, i_index

  global_matrix = (0.0d0, 0.0d0)

  do i_row = 1, n_my_rows
    do i_col = 1, n_my_cols
      ! The last conditional is needed, as often we will use a ScaLAPACK
      ! descriptor with a matrix with a certain leading dimension but fewer
      ! columns than the original matrix specified by the ScaLAPACK descriptor
      ! (for example, using sc_desc for a given spin channel for both hamiltonian
      ! and KS_eigenvector, even though the former is an n_basis*n_basis matrix
      ! and the later is an n_basis*n_states matrix)
      i_row_global = my_row(i_row)
      i_col_global = my_col(i_col)
      if( i_row_global > 0 .and. i_col_global > 0 .and. (i_row_global .le. i_col_global) .and. &
           (i_col_global .le. global_dim) ) then
        ! Offset for the current column
        i_index = ((i_col_global-1)*(i_col_global))/2
        ! Offset for the current row
        i_index = i_index + i_row_global

        global_matrix( i_index ) = local_matrix(i_row,i_col)
      endif
    end do
  end do
end subroutine convert_local_to_global_up_tri_matrix_complex

!-----------------------------------------------------------------------------------
!****s* convert_local_to_global_low_tri_matrix_complex
! NAME
!    convert_local_to_global_up_tri_matrix_complex
! SYNOPSIS
subroutine convert_local_to_global_low_tri_matrix_complex( n_my_rows, n_my_cols, my_row, my_col, &
                mxld, mxcol, local_matrix, global_dim, global_matrix )
!  PURPOSE
!    This subroutine fills in the elements for a complex (square) global matrix stored in a 1D
!    lower triangular format from a local (2D block-cyclic) matrix.  Ideally, this
!    subroutine shouldn't be needed; one should work with the local matrix only.  But when
!    working with legacy code, this may be necessary to avoid a substantial code rewrite.
!
!    Note:  unlike the original subroutine, we do *not* sync the matrix afterwards.  I can
!    envision debug cases where one would not want to do this.  Since MPI calls are a
!    significant bottleneck that should be avoided as much as possible, it's good practice
!    to force the programmer to do it themselves.
!
!    Assumes that the UPPER triangular half of the local matrix is set.
! USES
  implicit none
! ARGUMENTS
! INPUTS
  integer,    intent(in)  :: n_my_rows
  integer,    intent(in)  :: n_my_cols
  integer,    intent(in)  :: my_row(n_my_rows)
  integer,    intent(in)  :: my_col(n_my_cols)
  integer,    intent(in)  :: mxld
  integer,    intent(in)  :: mxcol
  complex*16, intent(in)  :: local_matrix( mxld, mxcol )
  integer,    intent(in)  :: global_dim
! OUTPUTS
  complex*16, intent(out) :: global_matrix( global_dim * (global_dim+1)/2 )
! AUTHOR
!   William Huhn
! HISTORY
!   September 2017 - Created.
! SOURCE
  integer :: i_col, i_col_global, i_row, i_row_global, i_index

  global_matrix = (0.0d0, 0.0d0)

  do i_row = 1, n_my_rows
    do i_col = 1, n_my_cols
      ! The last conditional is needed, as often we will use a ScaLAPACK
      ! descriptor with a matrix with a certain leading dimension but fewer
      ! columns than the original matrix specified by the ScaLAPACK descriptor
      ! (for example, using sc_desc for a given spin channel for both hamiltonian
      ! and KS_eigenvector, even though the former is an n_basis*n_basis matrix
      ! and the later is an n_basis*n_states matrix)
      i_row_global = my_row(i_row)
      i_col_global = my_col(i_col)
      if( i_row_global > 0 .and. i_col_global > 0 .and. (i_col_global .le. i_row_global) .and. &
           (i_col_global .le. global_dim) ) then
        ! Offset for the current column
        i_index = (global_dim * (global_dim+1))/2 - ((global_dim-i_col_global+1)*((global_dim-i_col_global+1)+1))/2
        ! Offset for the current row
        i_index = i_index + (i_row_global - i_col_global + 1)

        global_matrix( i_index ) = local_matrix(i_row,i_col)
      endif
    end do
  end do
end subroutine convert_local_to_global_low_tri_matrix_complex
