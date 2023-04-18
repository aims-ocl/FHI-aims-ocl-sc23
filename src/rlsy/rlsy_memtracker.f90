module rlsy_memtracker
!! Barebones memory tracker
use rlsy_constants, only: i4,i8,r8,rl_iou,rl_hugeint
implicit none
private
public :: rl_memtracker

!> keep track of the number of allocations/deallocations and memory used.
type rl_memtracker
    !> number of calls to allocate
    integer(i8) :: n_allocate=-rl_hugeint
    !> number of calls to deallocate
    integer(i8) :: n_deallocate=-rl_hugeint
    !> current amount of memory allocated
    integer(i8) :: persistent_scalable=-rl_hugeint
    integer(i8) :: persistent_nonscalable=-rl_hugeint
    integer(i8) :: temporary_scalable=-rl_hugeint
    integer(i8) :: temporary_nonscalable=-rl_hugeint
    !> keep track of peak memory of each kind?
    integer(i8) :: peak_persistent_scalable=-rl_hugeint
    integer(i8) :: peak_persistent_nonscalable=-rl_hugeint
    integer(i8) :: peak_temporary_scalable=-rl_hugeint
    integer(i8) :: peak_temporary_nonscalable=-rl_hugeint
    !> space for catching error codes
    integer, private :: ierr=0
    contains
        !> initialize the tracker
        procedure :: init

        !> allocate memory
        generic :: allocate=>&
            allocate_1d_i4,&
            allocate_2d_i4,&
            allocate_1d_r8,&
            allocate_2d_r8,&
            allocate_3d_r8,&
            allocate_4d_r8,&
            allocate_2d_c8,&
            allocate_3d_c8

        !> free memory
        generic :: deallocate=>&
            deallocate_1d_i4,&
            deallocate_2d_i4,&
            deallocate_1d_r8,&
            deallocate_2d_r8,&
            deallocate_3d_r8,&
            deallocate_4d_r8,&
            deallocate_2d_c8,&
            deallocate_3d_c8

        procedure, private :: allocate_1d_i4
        procedure, private :: allocate_2d_i4

        procedure, private :: allocate_1d_r8
        procedure, private :: allocate_2d_r8
        procedure, private :: allocate_3d_r8
        procedure, private :: allocate_4d_r8

        procedure, private :: allocate_2d_c8
        procedure, private :: allocate_3d_c8



        procedure, private :: deallocate_1d_i4
        procedure, private :: deallocate_2d_i4

        procedure, private :: deallocate_1d_r8
        procedure, private :: deallocate_2d_r8
        procedure, private :: deallocate_3d_r8
        procedure, private :: deallocate_4d_r8

        procedure, private :: deallocate_2d_c8
        procedure, private :: deallocate_3d_c8

        !> dump info to stdout
        procedure :: dump
end type

! Should probably throw an error if I try to allocate something very large.
! This can be supressed, but is good for catching errors when you allocate
! something you think is small but it turns out to be large!
integer, parameter :: size_for_throwing_error=1000000
contains

!> create the memory tracker
subroutine init(mem)
    !> memory tracker
    class(rl_memtracker), intent(out) :: mem

    ! Just set everything to zero
    mem%n_allocate=0
    mem%n_deallocate=0
    mem%persistent_scalable=0
    mem%persistent_nonscalable=0
    mem%temporary_scalable=0
    mem%temporary_nonscalable=0
    mem%peak_persistent_scalable=0
    mem%peak_persistent_nonscalable=0
    mem%peak_temporary_scalable=0
    mem%peak_temporary_nonscalable=0
    mem%ierr=0
end subroutine

!> dump snapshot @TODO add file and line
subroutine dump(mem)
    !> memory tracker
    class(rl_memtracker), intent(in) :: mem

    write(rl_iou,*) '      persistent scalable:',mem%persistent_scalable
    write(rl_iou,*) '   persistent nonscalable:',mem%persistent_nonscalable
    write(rl_iou,*) '       temporary scalable:',mem%temporary_scalable
    write(rl_iou,*) '    temporary nonscalable:',mem%temporary_nonscalable
end subroutine

subroutine allocate_1d_i4(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    integer(i4), dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( n .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_i4(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    integer(i4), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_1d_r8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( n .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_r8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_3d_r8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(3), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_4d_r8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(4), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3),n(4)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_c8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_3d_c8(mem,x,n,persistent,scalable,supress_error)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(3), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(rl_iou,*) '       that could not have been intended.'
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(rl_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(rl_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(rl_iou,*) 'ERROR: failed allocation.'
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine

subroutine deallocate_1d_i4(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    integer(i4), dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_i4(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    integer(i4), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_1d_r8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_r8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_3d_r8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_4d_r8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_c8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_3d_c8(mem,x,persistent,scalable)
    !> memory tracker
    class(rl_memtracker), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable

    integer :: n

    if ( .not.allocated(x) ) then
        write(rl_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(rl_iou,*) 'ERROR: failed deallocation.'
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine


end module rlsy_memtracker
