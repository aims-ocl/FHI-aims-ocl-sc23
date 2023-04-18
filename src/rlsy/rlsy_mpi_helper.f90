module rlsy_mpi_helper
!!
!! Wrapper for MPI things. Makes code extremely verbose without it.
!! Also I can never remember the syntax and always mess stuff up.
!!
use rlsy_constants, only: r8,rl_iou,rl_hugeint,rl_tiny,rl_exitcode_mpi
use mpi_tasks

implicit none
private

public :: rl_mpi_helper
public :: rl_stop_gracefully
public :: mpi_wtime
! Expose some MPI constants, in case someone needs them
public :: MPI_CHARACTER,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION
public :: MPI_INTEGER,MPI_LOGICAL
public :: MPI_MAX,MPI_MIN,MPI_SUM,MPI_IN_PLACE

!> Helper that keeps track of all things MPI.
type rl_mpi_helper
    !> which communicator
    integer :: comm=-rl_hugeint
    !> current rank
    integer :: r=-rl_hugeint
    !> total number of ranks
    integer :: n=-rl_hugeint
    !> exit status
    integer :: error=-rl_hugeint
    !> is this rank allowed to write to stdout?
    logical :: talk=.false.
    !> current status
    integer, dimension(:), allocatable :: stat
    contains
        !> initialize
        procedure :: init=>init_mpi
        !> close
        procedure :: destroy=>kill_mpi
        !> split communicator
        procedure :: split=>split_communicator
        !> barrier
        procedure :: barrier
        !> mpi allreduce
        generic :: allreduce=>&
            allreduce_int,&
            allreduce_1d_int,&
            allreduce_2d_int,&
            allreduce_3d_int,&
            allreduce_real,&
            allreduce_1d_real,&
            allreduce_2d_real,&
            allreduce_3d_real,&
            allreduce_4d_real,&
            allreduce_complex,&
            allreduce_2d_complex,&
            allreduce_3d_complex

        procedure, private :: allreduce_int
        procedure, private :: allreduce_1d_int
        procedure, private :: allreduce_2d_int
        procedure, private :: allreduce_3d_int
        procedure, private :: allreduce_real
        procedure, private :: allreduce_1d_real
        procedure, private :: allreduce_2d_real
        procedure, private :: allreduce_3d_real
        procedure, private :: allreduce_4d_real
        procedure, private :: allreduce_complex
        procedure, private :: allreduce_2d_complex
        procedure, private :: allreduce_3d_complex
        !> mpi broadcast
        generic :: bcast=>&
            broadcast_int,&
            broadcast_1d_int,&
            broadcast_2d_int

        procedure, private :: broadcast_int
        procedure, private :: broadcast_1d_int
        procedure, private :: broadcast_2d_int

        !> check and synchronize variables @TODO add tolerance for check and sync
        generic :: check_and_sync=>&
            check_and_sync_int,&
            check_and_sync_1d_int,&
            check_and_sync_1d_real,&
            check_and_sync_2d_real

        procedure, private :: check_and_sync_int
        procedure, private :: check_and_sync_1d_int
        procedure, private :: check_and_sync_1d_real
        procedure, private :: check_and_sync_2d_real
end type

contains

!> initialize MPI
subroutine init_mpi(mw,communicator,notalk)
   !> communicator helper
   class(rl_mpi_helper), intent(out) :: mw
   !> specific communicator, if not world
   integer, intent(in), optional :: communicator
   !> are all ranks forbidden to talk?
   logical, intent(in), optional :: notalk

   integer :: mpierr
   logical :: init

   ! start MPI, if not already done so
   call mpi_initialized(init,mpierr)
   if ( init .eqv. .false. ) call mpi_init(mpierr)

   ! communicator
   if ( present(communicator) ) then
       mw%comm=communicator
   else
       mw%comm=MPI_COMM_WORLD
   endif

   ! current rank
   call mpi_comm_rank(mw%comm,mw%r,mpierr)
   ! number of ranks
   call mpi_comm_size(mw%comm,mw%n,mpierr)
   ! which rank can talk? Always rank 0.
   mw%talk=.false.
   if ( mw%r .eq. 0 ) mw%talk=.true.
   if ( present(notalk) ) then
       if ( notalk ) then
           mw%talk=.false.
       endif
   endif

   ! status
   allocate(mw%stat(MPI_STATUS_SIZE))
   mw%stat=0
end subroutine

!> kill MPI
subroutine kill_mpi(mw)
   !> communicator helper
   class(rl_mpi_helper), intent(inout) :: mw

   integer :: mpierr
   ! just make sure everything is destroyed
   call mpi_finalize(mpierr)
   mw%comm=-1
   mw%r=-1
   mw%n=-1
   mpierr=-1
   mw%stat=-1
end subroutine

!> MPI barrier
subroutine barrier(mw)
    !> communicator
    class(rl_mpi_helper), intent(in) :: mw

    integer :: mpierr
    call MPI_Barrier(mw%comm,mpierr)
    if ( mpierr .ne. 0 ) then
        call rl_stop_gracefully(['MPI_Barrier exit code: '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
    endif
end subroutine

!> split an MPI communicator into chunks
subroutine split_communicator(mw,mn,color)
   !> communicator to split
   class(rl_mpi_helper), intent(in) :: mw
   !> split communicator
   type(rl_mpi_helper), intent(out) :: mn
   !> per-rank color, ranks with equal color end up in the same communicator
   integer, intent(in) :: color

   integer :: mpierr
   logical :: dlog

   ! Check that MPI is running
   call mpi_initialized(dlog,mpierr)
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['MPI error when checking initialization'],rl_exitcode_mpi,mw%comm)
   endif
   if ( dlog .eqv. .false. ) then
       call rl_stop_gracefully(['No MPI initialized when splitting.'],rl_exitcode_mpi,mw%comm)
   endif

   ! No error checks are done on the color. Not sure if good.
   call MPI_comm_split(mw%comm,color,mw%r,mn%comm,mpierr)
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['MPI error when splitting.'],rl_exitcode_mpi,mw%comm)
   endif

   ! current rank
   call mpi_comm_rank(mn%comm,mn%r,mn%error)
   ! number of ranks
   call mpi_comm_size(mn%comm,mn%n,mn%error)
   ! which rank can talk? None of them.
   mn%talk=.false.

   ! status
   allocate(mn%stat(MPI_STATUS_SIZE))
   mn%stat=0
   mn%error=0
end subroutine

!> create processor grid dimensions?

!> Will check that an int is synced across ranks. If it is not, maybe throw a warning/error.
subroutine check_and_sync_int(mw,x,ref,bigdeal,vname)
   !> mpi helper
   class(rl_mpi_helper), intent(in) :: mw
   !> int
   integer, intent(inout) :: x
   !> which rank to use as reference
   integer, intent(in) :: ref
   !> how big of an issue?
   integer, intent(in), optional :: bigdeal
   !> name of variable for message
   character(len=*), intent(in), optional :: vname

   integer :: j,loudness,mpierr

   if ( present(bigdeal) ) then
       loudness=bigdeal
   else
       loudness=2 ! Default to angry error
   endif

   j=0
   call mpi_allreduce(x,j,1,MPI_INTEGER,MPI_SUM,mw%comm,mpierr)
   j=j/mw%n
   if ( abs(x-j) .gt. 0 ) then
       select case(loudness)
       case(0)
           ! Just sync and move on with life
           call MPI_Bcast(x,1,MPI_INTEGER,ref,mw%comm,mpierr)
       case(1)
           ! Sync and throw warning
           call MPI_Bcast(x,1,MPI_INTEGER,ref,mw%comm,mpierr)
           if ( mw%talk ) then
               if ( present(vname) ) then
                   write(rl_iou,*) 'WARNING: "'//trim(vname)//'" out of sync:'
               else
                   write(rl_iou,*) 'WARNING: out of sync:'
               endif
               write(rl_iou,*) '    average:',j
               write(rl_iou,*) '     rank 0:',x
           endif
       case(2)
           ! KILL
           if ( present(vname) ) then
               call rl_stop_gracefully(['"'//trim(vname)//'" out of sync'],rl_exitcode_mpi,mw%comm)
           else
               call rl_stop_gracefully(['out of sync'],rl_exitcode_mpi,mw%comm)
           endif
       end select
   endif
end subroutine

subroutine check_and_sync_1d_int(mw,x,ref,bigdeal,vname)
   !> mpi helper
   class(rl_mpi_helper), intent(in) :: mw
   !> integer to check
   integer, dimension(:), intent(inout) :: x
   !> which rank to use as reference
   integer, intent(in) :: ref
   !> how big of an issue?
   integer, intent(in), optional :: bigdeal
   !> name of variable for message
   character(len=*), intent(in), optional :: vname

   integer, dimension(:), allocatable :: dr
   integer :: loudness,mpierr

   if ( present(bigdeal) ) then
       loudness=bigdeal
   else
       loudness=2
   endif
   allocate(dr(size(x)))
   dr=0
   call mpi_allreduce(x,dr,size(x),MPI_INTEGER,MPI_SUM,mw%comm,mpierr)
   dr=dr/mw%n
   if ( sum(abs(x-dr)) .ne. 0 ) then
       select case(loudness)
       case(0)
           ! Just sync and move on with life
           call MPI_Bcast(x,size(x),MPI_INTEGER,ref,mw%comm,mpierr)
       case(1)
           ! Sync and throw warning
           call MPI_Bcast(x,size(x),MPI_INTEGER,ref,mw%comm,mpierr)
           if ( mw%talk ) then
               if ( present(vname) ) then
                   write(rl_iou,*) 'WARNING: "'//trim(vname)//'" out of sync:'
               else
                   write(rl_iou,*) 'WARNING: out of sync:'
               endif
           endif
       case(2)
           ! KILL
           if ( present(vname) ) then
               call rl_stop_gracefully(['"'//trim(vname)//'" out of sync'],rl_exitcode_mpi,mw%comm)
           else
               call rl_stop_gracefully(['out of sync'],rl_exitcode_mpi,mw%comm)
           endif
       end select
   endif
end subroutine

! subroutine check_and_sync_real(mw,x,ref,bigdeal,vname,filename,linenumber)
!    !> mpi helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> real to check
!    real(r8), intent(inout) :: x
!    !> which rank to use as reference
!    integer, intent(in) :: ref
!    !> how big of an issue?
!    integer, intent(in), optional :: bigdeal
!    !> name of variable for message
!    character(len=*), intent(in), optional :: vname
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!
!    real(r8) :: dr
!    integer :: loudness
!
!    if ( present(bigdeal) ) then
!        loudness=bigdeal
!    else
!        loudness=2 ! Default to angry error
!    endif
!
!    dr=0.0_r8
!    call mpi_allreduce(x,dr,1,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mpierr)
!    dr=dr/mw%n
!    if ( abs(x-dr) .gt. 1E-10_r8 ) then
!        select case(loudness)
!        case(0)
!            ! Just sync and move on with life
!            call MPI_Bcast(x,1,MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
!        case(1)
!            ! Sync and throw warning
!            call MPI_Bcast(x,1,MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
!            if ( mw%talk ) then
!                if ( present(vname) ) then
!                    write(rl_iou,*) 'WARNING: "'//trim(vname)//'" out of sync:'
!                else
!                    write(rl_iou,*) 'WARNING: out of sync:'
!                endif
!                write(rl_iou,*) '    average:',dr
!                write(rl_iou,*) '     rank 0:',x
!            endif
!        case(2)
!            ! KILL
!            if ( present(vname) .and. present(filename) .and. present(linenumber) ) then
!                call rl_stop_gracefully(['"'//trim(vname)//'" out of sync'],rl_exitcode_mpi,filename,linenumber,mw%comm)
!            else
!                call rl_stop_gracefully(['out of sync'],rl_exitcode_mpi,mw%comm)
!            endif
!        end select
!    endif
! end subroutine

subroutine check_and_sync_1d_real(mw,x,ref,bigdeal,vname)
   !> mpi helper
   class(rl_mpi_helper), intent(in) :: mw
   !> real to check
   real(r8), dimension(:), intent(inout) :: x
   !> which rank to use as reference
   integer, intent(in) :: ref
   !> how big of an issue?
   integer, intent(in), optional :: bigdeal
   !> name of variable for message
   character(len=*), intent(in), optional :: vname

   real(r8), dimension(:), allocatable :: dr
   integer :: loudness,mpierr

   if ( present(bigdeal) ) then
       loudness=bigdeal
   else
       loudness=2
   endif
   allocate(dr(size(x)))
   dr=0.0_r8
   call mpi_allreduce(x,dr,size(x),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mpierr)
   dr=dr/mw%n
   if ( sum(abs(x-dr))/size(x) .gt. 1E-10_r8 ) then
       select case(loudness)
       case(0)
           ! Just sync and move on with life
           call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
       case(1)
           ! Sync and throw warning
           call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
           if ( mw%talk ) then
               if ( present(vname) ) then
                   write(rl_iou,*) 'WARNING: "'//trim(vname)//'" out of sync:'
               else
                   write(rl_iou,*) 'WARNING: out of sync:'
               endif
           endif
       case(2)
           ! KILL
           if ( present(vname) ) then
               call rl_stop_gracefully(['"'//trim(vname)//'" out of sync'],rl_exitcode_mpi,mw%comm)
           else
               call rl_stop_gracefully(['out of sync'],rl_exitcode_mpi,mw%comm)
           endif
       end select
   endif
end subroutine

subroutine check_and_sync_2d_real(mw,x,ref,bigdeal,vname)
   !> mpi helper
   class(rl_mpi_helper), intent(in) :: mw
   !> real to check
   real(r8), dimension(:,:), intent(inout) :: x
   !> which rank to use as reference
   integer, intent(in) :: ref
   !> how big of an issue?
   integer, intent(in), optional :: bigdeal
   !> name of variable for message
   character(len=*), intent(in), optional :: vname

   real(r8), dimension(:,:), allocatable :: dr
   integer :: loudness,mpierr

   if ( present(bigdeal) ) then
       loudness=bigdeal
   else
       loudness=2
   endif

   allocate(dr(size(x,1),size(x,2)))
   dr=x
   call MPI_Bcast(dr,size(dr),MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
   if ( sum(abs(x-dr)) .gt. rl_tiny ) then
       select case(loudness)
       case(0)
           ! Just sync and move on with life
           call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
       case(1)
           ! Sync and throw warning
           call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mpierr)
           if ( mw%talk ) then
               if ( present(vname) ) then
                   write(rl_iou,*) 'WARNING: "'//trim(vname)//'" out of sync:'
               else
                   write(rl_iou,*) 'WARNING: out of sync.'
               endif
           endif
       case(2)
           ! KILL
           if ( present(vname) ) then
               call rl_stop_gracefully(['"'//trim(vname)//'" out of sync'],rl_exitcode_mpi,mw%comm)
           else
               call rl_stop_gracefully(['out of sync'],rl_exitcode_mpi,mw%comm)
           endif
       end select
   endif
end subroutine

!> MPI Bcast for one integer
subroutine broadcast_int(mw,i,from)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> Integer to broadcast
   integer, intent(inout) :: i
   !> Where to broadcast from
   integer, intent(in) :: from

   integer :: mpierr
   call MPI_Bcast(i,1,MPI_INTEGER,from,mw%comm,mpierr)
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine

!> MPI Bcast for 1D ints
subroutine broadcast_1d_int(mw,d,from)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> Array to broadcast
   integer, dimension(:), intent(inout) :: d
   !> Where to broadcast from
   integer, intent(in) :: from

   integer :: mpierr
   call MPI_Bcast(d,size(d),MPI_INTEGER,from,mw%comm,mpierr)
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine

!> MPI Bcast for 2D ints
subroutine broadcast_2d_int(mw,d,from)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> Array to broadcast
   integer, dimension(:,:), intent(inout) :: d
   !> Where to broadcast from
   integer, intent(in) :: from

   integer :: mpierr
   call MPI_Bcast(d,size(d),MPI_INTEGER,from,mw%comm,mpierr)
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine

! !> MPI Bcast for one real
! subroutine broadcast_real(mw,d,from,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> Integer to broadcast
!    real(r8), intent(inout) :: d
!    !> Where to broadcast from
!    integer, intent(in) :: from
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!
!    call MPI_Bcast(d,1,MPI_DOUBLE_PRECISION,from,mw%comm,mpierr)
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine
! !> MPI Bcast for 1D reals
! subroutine broadcast_1d_real(mw,d,from,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> Array to broadcast
!    real(r8), dimension(:), intent(inout) :: d
!    !> Where to broadcast from
!    integer, intent(in) :: from
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!    ! Broadcast it
!    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mpierr)
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine
! !> MPI Bcast for 2D reals
! subroutine broadcast_2d_real(mw,d,from,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> Array to broadcast
!    real(r8), dimension(:,:), intent(inout) :: d
!    !> Where to broadcast from
!    integer, intent(in) :: from
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!    ! Broadcast it
!    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mpierr)
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine
! !> MPI Bcast for 3D reals
! subroutine broadcast_3d_real(mw,d,from,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> Array to broadcast
!    real(r8), dimension(:,:,:), intent(inout) :: d
!    !> Where to broadcast from
!    integer, intent(in) :: from
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!    ! Broadcast it
!    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mpierr)
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine
! !> MPI Bcast for 4D reals
! subroutine broadcast_4d_real(mw,d,from,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> Array to broadcast
!    real(r8), dimension(:,:,:,:), intent(inout) :: d
!    !> Where to broadcast from
!    integer, intent(in) :: from
!    !> filename we call from for debugging
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from for debugging
!    integer, intent(in), optional :: linenumber
!    ! Broadcast it
!    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mpierr)
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_bcast exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine

!> Allreduce for a single integer
subroutine allreduce_int(mw,operation,i,j)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> integer to allreduce
   integer, intent(inout) :: i
   !> destination to allgather to, if omitted default to in-place
   integer, intent(inout), optional :: j

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(j) ) then
       j=0
       call mpi_allreduce(i,j,1,MPI_INTEGER,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,i,1,MPI_INTEGER,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 1D integer array
subroutine allreduce_1d_int(mw,operation,i,j)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> integer to allreduce
   integer, intent(inout), dimension(:) :: i
   !> destination to allgather to, if omitted default to in-place
   integer, intent(inout), dimension(:), optional :: j

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(j) ) then
       j=0
       call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 2D integer array
subroutine allreduce_2d_int(mw,operation,i,j)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> integer to allreduce
   integer, intent(inout), dimension(:,:) :: i
   !> destination to allgather to, if omitted default to in-place
   integer, intent(inout), dimension(:,:), optional :: j

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(j) ) then
       j=0
       call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 3D integer array
subroutine allreduce_3d_int(mw,operation,i,j)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> integer to allreduce
   integer, intent(inout), dimension(:,:,:) :: i
   !> destination to allgather to, if omitted default to in-place
   integer, intent(inout), dimension(:,:,:), optional :: j

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(j) ) then
       j=0
       call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine

!> Allreduce for a single double-precision number
subroutine allreduce_real(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> number to allreduce
   real(r8), intent(inout) :: x
   !> destination to allreduce to, if omitted default to in-place
   real(r8), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       y=0.0_r8
       call mpi_allreduce(x,y,1,MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,1,MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 1D double-precision array
subroutine allreduce_1d_real(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(inout) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   real(r8), dimension(:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   real(r8), dimension(:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 2D double-precision array
subroutine allreduce_2d_real(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   real(r8), dimension(:,:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   real(r8), dimension(:,:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 3D double-precision array
subroutine allreduce_3d_real(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   real(r8), dimension(:,:,:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   real(r8), dimension(:,:,:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 4D double-precision array
subroutine allreduce_4d_real(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   real(r8), dimension(:,:,:,:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   real(r8), dimension(:,:,:,:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine

! !> Allreduce for a 4D double-precision array
! subroutine allreduce_4d_real(mw,operation,x,y,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> what operation to do
!    character(len=*), intent(in) :: operation
!    !> array to allreduce
!    real(r8), dimension(:,:,:,:), intent(inout) :: x
!    !> array to allreduce to, if omitted default to in-place
!    real(r8), dimension(:,:,:,:), intent(out), optional :: y
!    !> filename we call from
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from
!    integer, intent(in), optional :: linenumber
!
!    integer :: mpiop
!    ! Fetch proper operation code
!    mpiop=mpi_operation_code(operation)
!    ! Do the actual communication
!    if ( present(y) ) then
!        if ( size(x) .ne. size(y) ) then
!            call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        endif
!        y=0.0_r8
!        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
!    else
!        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
!    endif
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine
! subroutine allreduce_5d_real(mw,operation,x,y,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> what operation to do
!    character(len=*), intent(in) :: operation
!    !> array to allreduce
!    real(r8), dimension(:,:,:,:,:), intent(inout) :: x
!    !> array to allreduce to, if omitted default to in-place
!    real(r8), dimension(:,:,:,:,:), intent(out), optional :: y
!    !> filename we call from
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from
!    integer, intent(in), optional :: linenumber
!
!    integer :: mpiop
!    ! Fetch proper operation code
!    mpiop=mpi_operation_code(operation)
!    ! Do the actual communication
!    if ( present(y) ) then
!        if ( size(x) .ne. size(y) ) then
!            call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        endif
!        y=0.0_r8
!        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
!    else
!        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mpierr)
!    endif
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine

!> Allreduce for a single complex number
subroutine allreduce_complex(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   complex(r8), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   complex(r8), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       y=0.0_r8
       call mpi_allreduce(x,y,1,MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,1,MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 2D double-complex array
subroutine allreduce_2d_complex(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   complex(r8), dimension(:,:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   complex(r8), dimension(:,:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
!> Allreduce for a 3D double-complex array
subroutine allreduce_3d_complex(mw,operation,x,y)
   !> MPI helper
   class(rl_mpi_helper), intent(in) :: mw
   !> what operation to do
   character(len=*), intent(in) :: operation
   !> array to allreduce
   complex(r8), dimension(:,:,:), intent(inout) :: x
   !> array to allreduce to, if omitted default to in-place
   complex(r8), dimension(:,:,:), intent(out), optional :: y

   integer :: mpiop,mpierr
   mpiop=mpi_operation_code(operation)
   if ( present(y) ) then
       if ( size(x) .ne. size(y) ) then
           call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,mw%comm)
       endif
       y=0.0_r8
       call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   else
       call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
   endif
   if ( mpierr .ne. 0 ) then
       call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
   endif
end subroutine
! !> Allreduce for a 3D double-complex array
! subroutine allreduce_5d_complex(mw,operation,x,y,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> what operation to do
!    character(len=*), intent(in) :: operation
!    !> array to allreduce
!    complex(r8), dimension(:,:,:,:,:), intent(inout) :: x
!    !> array to allreduce to, if omitted default to in-place
!    complex(r8), dimension(:,:,:,:,:), intent(out), optional :: y
!    !> filename we call from
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from
!    integer, intent(in), optional :: linenumber
!
!    integer :: mpiop
!    ! Fetch proper operation code
!    mpiop=mpi_operation_code(operation)
!    ! Do the actual communication
!    if ( present(y) ) then
!        if ( size(x) .ne. size(y) ) then
!            call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        endif
!        y=0.0_r8
!        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
!    else
!        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mpierr)
!    endif
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_allreduce exit code '//tochar(mpierr)],rl_exitcode_mpi,__FILE__,__LINE__,mw%comm)
!        endif
!    endif
! end subroutine
!
! !> Reduce 2D real array
! subroutine reduce_2d_real(mw,operation,x,recvrank,y,filename,linenumber)
!    !> MPI helper
!    class(rl_mpi_helper), intent(inout) :: mw
!    !> what operation to do
!    character(len=*), intent(in) :: operation
!    !> array to reduce
!    real(r8), dimension(:,:), intent(inout) :: x
!    !> rank to reduce to
!    integer, intent(in) :: recvrank
!    !> array to reduce to, if omitted default to in-place
!    real(r8), dimension(:,:), intent(out), optional :: y
!    !> filename we call from
!    character(len=*), intent(in), optional :: filename
!    !> line number we call from
!    integer, intent(in), optional :: linenumber
!
!    integer :: mpiop
!    ! Fetch proper operation code
!    mpiop=mpi_operation_code(operation)
!    ! Do the actual communication
!    if ( present(y) ) then
!        if ( size(x) .ne. size(y) ) then
!            call rl_stop_gracefully(['mpi_allreduce inconsistent array sizes'],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        endif
!        y=0.0_r8
!        call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mpierr)
!    else
!        if ( mw%r .eq. recvrank ) then
!            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mpierr)
!        else
!            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mpierr)
!        endif
!    endif
!    ! Check that things went ok
!    if ( mpierr .ne. 0 ) then
!        if ( present(filename) .and. present(linenumber) ) then
!            call rl_stop_gracefully(['mpi_reduce exit code '//tochar(mpierr)],rl_exitcode_mpi,filename,linenumber,mw%comm)
!        else
!            call rl_stop_gracefully(['mpi_reduce exit code '//tochar(mpierr)],rl_exitcode_mpi,mw%comm)
!        endif
!    endif
! end subroutine

!> return the proper MPI operation from my string thingy.
function mpi_operation_code(operation) result(code)
   !> string describing the operation
   character(len=*), intent(in) :: operation
   !> proper MPI code
   integer :: code

   ! Decide what operation to use
   select case(trim(adjustl(operation)))
   case('sum')
       code=MPI_SUM
   case('max')
       code=MPI_MAX
   case('min')
       code=MPI_MIN
   case default
       call rl_stop_gracefully(['Unknown MPI operation: '//operation],rl_exitcode_mpi)
   end select
end function

!> Stop gracefully, but this version handles and kills MPI as well. Had to repeat it, only way to not get circular dependencies I think.
subroutine rl_stop_gracefully(msg,exitcode,communicator)
   !> message
   character(len=*), dimension(:), intent(in) :: msg
   !> what exit code to give?
   integer, intent(in) :: exitcode
   !> if it's an MPI application, send in the communicator
   integer, intent(in), optional :: communicator

   integer :: rank,nrank,i,errcode

   write(rl_iou,*) ''
   write(rl_iou,*) 'ERROR'
   select case(exitcode)
   case(1)
       write(rl_iou,*) 'exit code 1: unexpected dimensions'
   case(2)
       write(rl_iou,*) 'exit code 1: blas/lapack returned nonzero exitcode'
   case(3)
       write(rl_iou,*) 'exit code 3: unphysical value detected'
   case(4)
       write(rl_iou,*) 'exit code 4: symmetry error'
   case(5)
       write(rl_iou,*) 'exit code 5: bad parameters sent to routine'
   case(6)
       write(rl_iou,*) 'exit code 6: I/O error'
   case(7)
       write(rl_iou,*) 'exit code 7: MPI error'
   case(8)
       write(rl_iou,*) 'exit code 8: Memory error'
   end select
   write(rl_iou,*) ''
   do i=1,size(msg)
       write(rl_iou,*) trim(msg(i))
   enddo
   write(rl_iou,*) ''

   if ( present(communicator) ) then
       call mpi_comm_rank(communicator,rank,errcode)
       call mpi_comm_size(communicator,nrank,errcode)
       call mpi_finalize(errcode)
   endif

   ! Seems error codes are f2008, don't use those
   stop

   ! Seems I can not use a variable as the errorcode with Ifort.
   ! No, worries, got a brute force solution:
   ! select case(exitcode)
   ! case(1)
   !     error stop 1
   ! case(2)
   !     error stop 2
   ! case(3)
   !     error stop 3
   ! case(4)
   !     error stop 4
   ! case(5)
   !     error stop 5
   ! case(6)
   !     error stop 6
   ! case(7)
   !     error stop 7
   ! case default
   !     error stop
   ! end select
end subroutine

! Need a copy of tochar here (from rlsy_helpers) to
! avoid circular dependencies. Ugly -- I know.
pure function tochar(i,padding) result(s)
    !> integer to convert
    integer, intent(in) :: i
    !> pad the integer? Positive number means zero-padding, negative means pad with whitespace
    integer, intent(in), optional :: padding
    !> resulting string
    character(len=:), allocatable :: s

    character(len=range(i)) :: tmp
    character(len=range(i)+5) :: ttmp
    integer :: j,k

    if ( present(padding) ) then
        write(tmp,'(i0)') i ! get i to a string
        ! build a string that contains the padding
        do j=1,range(i)+5
            if ( padding > 0 ) then
                ttmp(j:j)='0'
            else
                ttmp(j:j)=' '
            endif
        enddo
        ! wrap it together to a nice string
        k=len(trim(adjustl(tmp))) ! how many digits where relevant
        s=ttmp(1:abs(padding)-k)//trim(adjustl(tmp))
    else
        ! much easier
        write(tmp,'(i0)') i
        s=trim(adjustl(tmp))
    endif
end function

end module
