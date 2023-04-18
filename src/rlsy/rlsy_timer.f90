module rlsy_timer
!! Barebones memory tracker
use rlsy_constants, only: i4,i8,r8,rl_huge,rl_hugeint,rl_iou
use rlsy_helpers, only: rl_mean,rl_stddev
use rlsy_mpi_helper, only: rl_mpi_helper,mpi_wtime

implicit none
private
public :: rl_timer
public :: rl_timer_eigenvalue
public :: rl_timer_density
public :: rl_timer_overlap
public :: rl_timer_multipole

!> timer for density update
type rl_timer_density
    real(r8) :: init        =0.0_r8
    real(r8) :: idle        =0.0_r8
    real(r8) :: total       =0.0_r8
    real(r8) :: fetchdm     =0.0_r8
    real(r8) :: reset       =0.0_r8
    real(r8) :: tab_coord   =0.0_r8
    real(r8) :: prune_basis =0.0_r8
    real(r8) :: tab_basis   =0.0_r8
    real(r8) :: prune_dm    =0.0_r8
    real(r8) :: matrixop    =0.0_r8
end type

!> timings for eigenvalue solution
type rl_timer_eigenvalue
    real(r8) :: init=0.0_r8
    real(r8) :: idle=0.0_r8
    real(r8) :: total=0.0_r8
    real(r8) :: forward_ft=0.0_r8
    real(r8) :: solve=0.0_r8
    real(r8) :: occnumbers=0.0_r8
    real(r8) :: densitymatrix=0.0_r8
end type

!> timings for overlap integration
type rl_timer_overlap
    real(r8) :: init=0.0_r8
    real(r8) :: idle=0.0_r8
    real(r8) :: total=0.0_r8
    real(r8) :: pairlocator=0.0_r8
    real(r8) :: reset=0.0_r8
    real(r8) :: tab_coord=0.0_r8
    real(r8) :: tab_basis=0.0_r8
    real(r8) :: matrixop=0.0_r8
    real(r8) :: communication=0.0_r8
end type

!> timings for multipole expansion
type rl_timer_multipole
    real(r8) :: init=0.0_r8
    real(r8) :: idle=0.0_r8
    real(r8) :: total=0.0_r8
    real(r8) :: pairlocator=0.0_r8
    real(r8) :: reset=0.0_r8
    real(r8) :: tab_coord=0.0_r8
    real(r8) :: tab_basis=0.0_r8
    real(r8) :: matrixop=0.0_r8
    real(r8) :: communication=0.0_r8
end type

!> timing information, collected in one place
type rl_timer
    type(rl_timer_eigenvalue) :: eigenvalue
    type(rl_timer_density) :: density
    type(rl_timer_overlap) :: overlap
    type(rl_timer_multipole) :: multipole
    real(r8) :: tstart=-rl_huge
    real(r8) :: tstop=-rl_huge
    contains
        !> initialize the timer
        procedure :: init
        !> dump timing information
        procedure :: dump
end type

! Should probably throw an error if I try to allocate something very large.
! This can be supressed, but is good for catching errors when you allocate
! something you think is small but it turns out to be large!
contains

!> create the memory tracker
subroutine init(mem)
    !> memory tracker
    class(rl_timer), intent(out) :: mem
end subroutine

!> dump timings
subroutine dump(tmr,mw,iou)
    !> timers
    class(rl_timer), intent(in) :: tmr
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> unit to dump to
    integer, intent(in) :: iou

    integer, parameter :: n_timer_things=1000 ! Or something
    integer, parameter :: n_major_sections=3
    real(r8), dimension(n_major_sections) :: secbuf,secpercent
    real(r8), dimension(100) :: brkbuf,brkpercent
    real(r8), dimension(:,:), allocatable :: buf
    real(r8), dimension(:), allocatable :: mean,dev,maxv,minv
    integer :: i

    ! Accumulate timings for all the things?
    allocate(buf(mw%n,n_timer_things))
    allocate(mean(n_timer_things))
    allocate(dev(n_timer_things))
    allocate(maxv(n_timer_things))
    allocate(minv(n_timer_things))
    buf=0.0_r8
    mean=0.0_r8
    dev=0.0_r8

    ! Eigenvalue
    buf(mw%r+1,1)=tmr%eigenvalue%total
    buf(mw%r+1,2)=tmr%eigenvalue%init
    buf(mw%r+1,3)=tmr%eigenvalue%idle
    buf(mw%r+1,4)=tmr%eigenvalue%forward_ft
    buf(mw%r+1,5)=tmr%eigenvalue%solve
    buf(mw%r+1,6)=tmr%eigenvalue%occnumbers
    buf(mw%r+1,7)=tmr%eigenvalue%densitymatrix

    ! Overlap
    buf(mw%r+1,11)=tmr%overlap%total
    buf(mw%r+1,12)=tmr%overlap%init
    buf(mw%r+1,13)=tmr%overlap%idle
    buf(mw%r+1,14)=tmr%overlap%pairlocator
    buf(mw%r+1,15)=tmr%overlap%reset
    buf(mw%r+1,16)=tmr%overlap%tab_coord
    buf(mw%r+1,17)=tmr%overlap%tab_basis
    buf(mw%r+1,18)=tmr%overlap%matrixop

    ! Density update
    buf(mw%r+1,21)=tmr%density%total
    buf(mw%r+1,22)=tmr%density%init
    buf(mw%r+1,23)=tmr%density%idle
    buf(mw%r+1,24)=tmr%density%fetchdm
    buf(mw%r+1,25)=tmr%density%reset
    buf(mw%r+1,26)=tmr%density%tab_coord
    buf(mw%r+1,27)=tmr%density%prune_basis
    buf(mw%r+1,28)=tmr%density%tab_basis
    buf(mw%r+1,29)=tmr%density%prune_dm
    buf(mw%r+1,30)=tmr%density%matrixop

    ! Get the mean and standard deviation for each thing
    call mw%allreduce('sum',buf)
    do i=1,n_timer_things
        mean(i)=rl_mean(buf(:,i))
        dev(i)=rl_stddev(buf(:,i))
        maxv(i)=maxval(buf(:,i))
        minv(i)=minval(buf(:,i))
    enddo

    ! ... dump timings
    if ( mw%talk ) then

        write(iou,*) ''
        write(iou,*) 'TIMING INFORMATION:'

        ! First dump major breakdown:
        secbuf(1)=maxv(1)    ! Eigenvalue solution
        secbuf(2)=maxv(11)   ! Overlap integration
        secbuf(3)=maxv(21)   ! Density update
        secpercent=100.0_r8*secbuf/sum(secbuf)

        write(iou,"(1X,A20,1X,F14.5,'s',4X,'(',F6.2,' %)')") 'Eigenproblem:',secbuf(1),secpercent(1)
        write(iou,"(1X,A20,1X,F14.5,'s',4X,'(',F6.2,' %)')") 'Overlap:',secbuf(2),secpercent(2)
        write(iou,"(1X,A20,1X,F14.5,'s',4X,'(',F6.2,' %)')") 'Density update:',secbuf(3),secpercent(3)

        ! Breakdown per part
        brkbuf=0.0_r8
        brkpercent=0.0_r8
        brkbuf(2:7)=maxv(2:7)
        brkpercent=100.0_r8*brkbuf/sum(brkbuf)
        write(iou,*) ''
        write(iou,*) 'Breakdown of eigenproblem timings:'
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")    'initialization:',mean(2),maxv(2)-minv(2),brkpercent(2)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")              'idle:',mean(3),maxv(3)-minv(3),brkpercent(3)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')") 'Fourier transform:',mean(4),maxv(4)-minv(4),brkpercent(4)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")       'ELSI solver:',mean(5),maxv(5)-minv(5),brkpercent(5)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')") 'Occupation number:',mean(6),maxv(6)-minv(6),brkpercent(6)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")     'Densitymatrix:',mean(7),maxv(7)-minv(7),brkpercent(7)

        brkbuf=0.0_r8
        brkpercent=0.0_r8
        brkbuf(2:8)=maxv(12:18)
        brkpercent=100.0_r8*brkbuf/sum(brkbuf)
        write(iou,*) ''
        write(iou,*) 'Breakdown of overlap integration timings:'
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")        'initialization:',mean(12),maxv(12)-minv(12),brkpercent(2)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")                  'idle:',mean(13),maxv(13)-minv(13),brkpercent(3)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")     'Build pairlocator:',mean(14),maxv(14)-minv(14),brkpercent(4)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")                 'Reset:',mean(15),maxv(15)-minv(15),brkpercent(5)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')") 'Calculate coordinates:',mean(16),maxv(16)-minv(16),brkpercent(6)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")     'Evaluate basis fn:',mean(17),maxv(17)-minv(17),brkpercent(7)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,' %)')")     'Matrix operations:',mean(18),maxv(18)-minv(18),brkpercent(8)

        brkbuf=0.0_r8
        brkpercent=0.0_r8
        brkbuf(2:10)=maxv(22:30)
        brkpercent=100.0_r8*brkbuf/sum(brkbuf)
        write(iou,*) ''
        write(iou,*) 'Breakdown of density update timings:'
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")        'initialization:',mean(22),maxv(22)-minv(22),brkpercent(2)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")                  'idle:',mean(23),maxv(23)-minv(23),brkpercent(3)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")           'Get full DM:',mean(24),maxv(24)-minv(24),brkpercent(4)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")                 'Reset:',mean(25),maxv(25)-minv(25),brkpercent(5)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')") 'Calculate coordinates:',mean(26),maxv(26)-minv(26),brkpercent(6)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")       'Prune basis set:',mean(27),maxv(27)-minv(27),brkpercent(7)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")     'Evaluate basis fn:',mean(28),maxv(28)-minv(28),brkpercent(8)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")              'Prune DM:',mean(29),maxv(29)-minv(29),brkpercent(9)
        write(iou,"(1X,A20,1X,F14.5,'s',1X,F14.5,4X,'(',F6.3,'%)')")     'Matrix operations:',mean(30),maxv(30)-minv(30),brkpercent(10)

    endif
end subroutine

end module rlsy_timer
