module rlsy_scalapack_helper
!!
!! Wrapper to help with basic Scalapack things
!!
use rlsy_constants, only: r8,rl_hugeint,rl_exitcode_mpi
use rlsy_helpers, only: tochar
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully

implicit none
private
public :: rl_blacs_helper

!> BLACS parameters, essentially the same as the MPI parameters.
type rl_blacs_helper
    !> context, suppose it's like an MPI communicator
    integer :: icontxt=-rl_hugeint
    !> number of processors in this context
    integer :: n=-rl_hugeint
    !> dimensions of processor grid
    integer :: n_proc_row=-rl_hugeint
    integer :: n_proc_col=-rl_hugeint
    !> where in the grid am I?
    integer :: proc_row=-rl_hugeint
    integer :: proc_col=-rl_hugeint
    !> error code
    integer :: error=-rl_hugeint
    contains
        !> create the blacs handle
        procedure :: init=>init_blacs
        !> destroy blacs handle
        procedure :: destroy=>destroy_blacs
        !> decide on a sensible block size
        procedure :: get_blocksize
        !> return number of local rows for distributed matrix
        procedure :: count_row_local
        !> return number of local columns for distributed matrix
        procedure :: count_col_local
        !> return rank for a pair of global indices
        procedure :: rank_from_global_indices
        !> convert indices from local to global and vice versa
        procedure :: local_row_from_global
        procedure :: local_col_from_global
        procedure :: global_row_from_local
        procedure :: global_col_from_local
end type

! I always start indexing my processors at zero. See no point in ever doing anything else.
integer, parameter :: proc_row_zero=0
integer, parameter :: proc_col_zero=0

contains

!> initialize BLACS
subroutine init_blacs(bw,mw)
    !> blacs helper
    class(rl_blacs_helper), intent(out) :: bw
    !> MPI communicator to use for this BLACS context
    type(rl_mpi_helper), intent(inout) :: mw

    integer, dimension(2) :: dims

    ! get neat dimensions of the processor grid.
    bw%error=0
    dims=0
    call MPI_Dims_create(mw%n,2,dims,mw%error)
    bw%n_proc_row=dims(1)
    bw%n_proc_col=dims(2)

    ! set the context to a color to create sensible groups?
    bw%icontxt=mw%comm
    call blacs_gridinit(bw%icontxt,'R',bw%n_proc_row,bw%n_proc_col,bw%error)
    if ( bw%error .ne. 0 ) then
        call rl_stop_gracefully(['blacs gridinit exit code: '//tochar(bw%error)],rl_exitcode_mpi,mw%comm)
    endif

    call blacs_gridinfo(bw%icontxt,bw%n_proc_row,bw%n_proc_col,bw%proc_row,bw%proc_col,bw%error)
    if ( bw%error .ne. 0 ) then
        call rl_stop_gracefully(['blacs gridinfo exit code: '//tochar(bw%error)],rl_exitcode_mpi,mw%comm)
    endif

    ! Sanity check
    bw%n=bw%n_proc_row*bw%n_proc_col
    if ( bw%n .ne. mw%n ) then
        call rl_stop_gracefully(['Failed creating BLACS context in a sensible way.'],rl_exitcode_mpi,mw%comm)
    endif
end subroutine

!> Local row index from global
function local_row_from_global(bw,blocksize,i) result(li)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> block size, same in both directions
    integer, intent(in) :: blocksize
    !> row index (global)
    integer, intent(in) :: i
    !> row index (local)
    integer :: li

    li  = ( (i-1)/(bw%n_proc_row*blocksize) )*blocksize + mod(i-1,blocksize) + 1
end function

!> Local column index from global
function local_col_from_global(bw,blocksize,j) result(li)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> block size, same in both directions
    integer, intent(in) :: blocksize
    !> col index (global)
    integer, intent(in) :: j
    !> col index (local)
    integer :: li

    li  = ( (j-1)/(bw%n_proc_col*blocksize) )*blocksize + mod(j-1,blocksize) + 1
end function

!> Global row index from local
function global_row_from_local(bw,blocksize,li) result(i)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> block size, same in both directions
    integer, intent(in) :: blocksize
    !> row index (local)
    integer, intent(in) :: li
    !> row index (global)
    integer :: i

    ! I pray this is correct. Pretty sure about it.
    i = (( ((li-1)/blocksize) * bw%n_proc_row) + bw%proc_row) *blocksize + mod(li-1,blocksize) + 1
end function

!> Global column index from local
function global_col_from_local(bw,blocksize,lj) result(j)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> block size, same in both directions
    integer, intent(in) :: blocksize
    !> col index (local)
    integer, intent(in) :: lj
    !> col index (global)
    integer :: j

    ! I pray this is correct. Pretty sure about it.
    j = (( ((lj-1)/blocksize) * bw%n_proc_col) + bw%proc_col) *blocksize + mod(lj-1,blocksize) + 1
end function

!> Returns the rank an element should be stored on from a pair of global indices. Same blocksize in all directions.
function rank_from_global_indices(bw,nrow_global,ncol_global,blocksize,i,j) result(rnk)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> number of rows in distributed matrix.
    integer, intent(in) :: nrow_global
    !> number of columns in distributed matrix.
    integer, intent(in) :: ncol_global
    !> block size, same in both directions
    integer, intent(in) :: blocksize
    !> row index (global)
    integer, intent(in) :: i
    !> column index (global)
    integer, intent(in) :: j
    !> resulting rank
    integer :: rnk

    integer :: pi,pj

    pi=mod((i-1)/blocksize,bw%n_proc_row)
    pj=mod((j-1)/blocksize,bw%n_proc_col)
    rnk=pi*bw%n_proc_col+pj
end function

!> Counts number of rows that end up on a given rank.
function count_row_local(bw, nrow_global, blocksize) result(nrow)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> The number of rows/columns in distributed matrix.
    integer, intent(in) :: nrow_global
    !> Block size, size of the blocks the distributed matrix is split into.
    integer, intent(in) :: blocksize
    !> Number of local rows
    integer :: nrow

    integer :: extrablks, nblocks

    ! Figure the total number of whole blocks nrow_global is split up into
    nblocks = nrow_global / blocksize
    ! Figure the minimum number of rows/cols a process can have
    nrow = (nblocks/bw%n_proc_row) * blocksize
    ! See if there are any extra blocks
    extrablks = mod( nblocks, bw%n_proc_row )
    ! If I have an extra block
    if( bw%proc_row .lt. extrablks ) then
        nrow = nrow + blocksize
    elseif( bw%proc_row .eq. extrablks ) then
        ! If I have last block, it may be a partial block
        nrow = nrow + mod( nrow_global, blocksize )
    end if
end function

!> Counts number of columns that end up on a given rank.
function count_col_local(bw, ncol_global, blocksize) result(ncol)
    !> BLACS helper
    class(rl_blacs_helper), intent(in) :: bw
    !> The number of rows/columns in distributed matrix.
    integer, intent(in) :: ncol_global
    !> Block size, size of the blocks the distributed matrix is split into.
    integer, intent(in) :: blocksize
    !> Number of local rows
    integer :: ncol

    integer :: extrablks, nblocks

    ! Figure the total number of whole blocks nrow_global is split up into
    nblocks = ncol_global / blocksize
    ! Figure the minimum number of rows/cols a process can have
    ncol = (nblocks/bw%n_proc_col) * blocksize
    ! See if there are any extra blocks
    extrablks = mod( nblocks, bw%n_proc_col )
    ! If I have an extra block
    if( bw%proc_col .lt. extrablks ) then
        ncol = ncol + blocksize
    elseif( bw%proc_col .eq. extrablks ) then
        ! If I have last block, it may be a partial block
        ncol = ncol + mod( ncol_global, blocksize )
    end if
end function

!> decide on a sensible block size for BLACS matrices. Same in both direction.
function get_blocksize(bw,nr,nc) result(bz)
    !> blacs helper
    class(rl_blacs_helper), intent(in) :: bw
    !> rows in matrix
    integer, intent(in) :: nr
    !> columns in matrix
    integer, intent(in) :: nc
    !> resulting blocksize
    integer :: bz

    integer :: nelem
    integer :: ir,ic,ni,nj,l,i

    ! Total number of elements in the matrix
    nelem=nr*nc
    ! If there are fewer elements than ranks, it is problematic.
    ! Probably better to tell the user to change settings than
    ! to add tons of code to handle it. It is a weird case.
    if ( nelem .lt. bw%n ) then
        bz=-1
        return
    endif


    bz=16   ! Seems to be a good choice for ELPA.
    sizeloop: do i=1,5
        ! This counts the minimum number of elements per rank
        l=rl_hugeint
        do ir=0,bw%n_proc_row-1
        do ic=0,bw%n_proc_col-1
            ni=dummycount(ir,nr,bw%n_proc_row,bz)
            nj=dummycount(ic,nc,bw%n_proc_col,bz)
            l=min(l,ni*nj)
        enddo
        enddo
        ! Check if it's all nonzero?
        if ( l .eq. 0 ) then
            ! That means some rank becomes empty. That is likely not a good idea.
            bz=bz/2
        else
            ! All ranks have something. Good enough!
            return
        endif
        if ( i .eq. 5 ) then
            ! This means everything failed miserably. Should never happen.
            bz=-1
            return
        endif
    enddo sizeloop

    contains

    function dummycount(i,n,np,bz) result(nrc)
        integer, intent(in) :: i,n,np,bz
        integer :: nrc

        integer :: nblocks,extrablks
        ! Figure the total number of whole blocks nrow_global is split up into
        nblocks = n / bz
        ! Figure the minimum number of rows/cols a process can have
        nrc = (nblocks/np) * bz
        ! See if there are any extra blocks
        extrablks = mod( nblocks, np )
        ! If I have an extra block
        if( i .lt. extrablks ) then
            nrc = nrc + bz
        elseif( i .eq. extrablks ) then
            ! If I have last block, it may be a partial block
            nrc = nrc + mod( n, bz )
        end if
    end function

end function

!> destroy BLACS.
subroutine destroy_blacs(bw)
    !> blacs helper
    class(rl_blacs_helper), intent(inout) :: bw

    ! the "1" means that I should be able to do normal MPI after this.
    ! not sure I should use this as the destructor, maybe it is too angry.
    ! call blacs_exit(1)

    ! Free buffers. Not sure if the call below already does this. The 0 means don't wait.
    call blacs_freebuff( bw%icontxt, 0 )
    ! Destroys this specific context
    call blacs_gridexit( bw%icontxt )
    ! Make sure everything is nonsense.
    bw%icontxt=-rl_hugeint
    bw%n=-rl_hugeint
    bw%n_proc_row=-rl_hugeint
    bw%n_proc_col=-rl_hugeint
    bw%proc_row=-rl_hugeint
    bw%proc_col=-rl_hugeint
    bw%error=-rl_hugeint
end subroutine

end module
