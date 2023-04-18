module rlsy_verletlist
!!
!! Verlet lists for fast distance lookups.
!! @TODO update for periodic boundary conditions when needed
!!
use rlsy_constants, only: r8,rl_hugeint,rl_huge,rl_sqtol,rl_exitcode_symmetry
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_stop_gracefully

implicit none
private
public :: rl_verletbox

!> list of points in a verlet box
type rl_verletbox_box
    !> how many points in this box
    integer :: n=-rl_hugeint
    !> indices to points in this box
    integer, dimension(:), allocatable :: ind
end type

!> minimal Verlet-box to generate distancelists and things like that
type rl_verletbox
    !> box divisions
    integer :: nx=-rl_hugeint,ny=-rl_hugeint,nz=-rl_hugeint
    !> lower bounds
    real(r8), dimension(3) :: rmin=rl_huge
    !> upper bounds
    real(r8), dimension(3) :: rmax=rl_huge
    !> scalefactor per dimension
    real(r8), dimension(3) :: ird=rl_huge
    !> boxes with points in them
    type(rl_verletbox_box), dimension(:,:,:), allocatable :: box
    contains
        !> stuff the particles into boxes
        procedure :: generate=>add_particles_in_boxes
        !> box-indices from a point
        procedure :: boxind=>boxind_from_coordinate
        !> locate index of a point
        procedure :: locate=>locate_index_of_point
        !> size in memory
        procedure :: size_in_mem
        !> destroy
        procedure :: destroy
end type

contains

!> measure size in memory, in bytes
pure function size_in_mem(vb) result(mem)
    !> verlet boxes
    class(rl_verletbox), intent(in) :: vb
    !> size in memory, in bytes
    integer :: mem

    integer :: i,j,k
    mem=0
    mem=mem+storage_size(vb)
    if ( allocated(vb%box) ) then
        do k=1,size(vb%box,3)
        do j=1,size(vb%box,2)
        do i=1,size(vb%box,1)
            mem=mem+storage_size(vb%box(i,j,k))
            if ( allocated(vb%box(i,j,k)%ind) ) mem=mem+storage_size(vb%box(i,j,k)%ind)*size(vb%box(i,j,k)%ind)
        enddo
        enddo
        enddo
    endif
    mem=mem/8
end function

!> destroy the object
subroutine destroy(vb,mem)
    !> verlet boxes
    class(rl_verletbox), intent(inout) :: vb
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    integer :: i,j,k

    if ( allocated(vb%box) ) then
        do k=1,size(vb%box,3)
        do j=1,size(vb%box,2)
        do i=1,size(vb%box,1)
            if ( allocated(vb%box(i,j,k)%ind) ) then
                call mem%deallocate(vb%box(i,j,k)%ind,persistent=.false.,scalable=.false.)
            endif
        enddo
        enddo
        enddo
        deallocate(vb%box)
    endif
end subroutine

!> get the indices of a box from the coordinate of a point
subroutine boxind_from_coordinate(vb,r,bi,bj,bk)
    !> verlet boxes
    class(rl_verletbox), intent(in) :: vb
    !> coordinates of point
    real(r8), dimension(3), intent(in) :: r
    !> box-indices
    integer, intent(out) :: bi,bj,bk

    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
end subroutine

!> Locate the index of a point using Verlet boxes. Index -1 means the point does not exist.
function locate_index_of_point(vb,points,r,singlebox) result(ind)
    !> verlet boxes
    class(rl_verletbox), intent(in) :: vb
    !> all possible points
    real(r8), dimension(:,:), intent(in) :: points
    !> coordinates of point
    real(r8), dimension(3), intent(in) :: r
    !> check only the first box that comes to mind
    logical, intent(in), optional :: singlebox
    !> index of points
    integer :: ind

    real(r8), dimension(3) :: w
    real(r8) :: f0
    integer :: bi,bj,bk,i,j,k,l,ii

    ! Assume I don't find anything
    ind=-1
    ! First some sanity tests
    j=0
    do i=1,3
        if ( r(i) .lt. vb%rmin(i) ) j=j+1
        if ( r(i) .gt. vb%rmax(i) ) j=j+1
    enddo
    if ( j .ne. 0 ) then
        ! This means the point is not in any of the boxes.
        return
    endif

    ! Which box to start looking?
    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
    ! First look in this box, decent odds the point is there
    do l=1,vb%box(bi,bj,bk)%n
        ii=vb%box(bi,bj,bk)%ind(l)
        w=points(:,ii)-r
        f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
        if ( f0 .lt. rl_sqtol ) then
            ! Yup, found it
            ind=ii
            return
        endif
    enddo

    if ( present(singlebox) ) then
        ! Maybe I only want to check one box at a time.
        if ( singlebox ) return
    endif

    ! Then look for it in the adjacent boxes
    do i=max(1,bi-1),min(vb%nx,bi+1)
    do j=max(1,bj-1),min(vb%ny,bj+1)
    do k=max(1,bk-1),min(vb%nz,bk+1)
        do l=1,vb%box(i,j,k)%n
            ii=vb%box(i,j,k)%ind(l)
            w=points(:,ii)-r
            f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
            if ( f0 .lt. rl_sqtol ) then
                ind=ii
                return
            endif
        enddo
    enddo
    enddo
    enddo
end function

!> Put particles in Verlet boxes
subroutine add_particles_in_boxes(vb,r,ndim,mem)
    !> verlet boxes
    class(rl_verletbox), intent(out) :: vb
    !> particles
    real(r8), dimension(:,:), intent(in) :: r
    !> number of boxes
    integer, dimension(3), intent(in) :: ndim
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    real(r8), dimension(3) :: v0
    integer :: n,i,j,k,l

    ! number of particles
    n=size(r,2)

    ! First get the bounds
    vb%rmin=rl_huge
    vb%rmax=-rl_huge
    do i=1,n
    do j=1,3
        vb%rmin(j)=min(vb%rmin(j),r(j,i))
        vb%rmax(j)=max(vb%rmax(j),r(j,i))
    enddo
    enddo
    ! Slight tolerance, should be fine.
    v0=(vb%rmax-vb%rmin)*rl_sqtol+rl_sqtol
    vb%rmin=vb%rmin-v0
    vb%rmax=vb%rmax+v0
    ! Scaling factor
    vb%ird=(ndim*1.0_r8)/(vb%rmax-vb%rmin)
    ! Dimensions
    vb%nx=ndim(1)
    vb%ny=ndim(2)
    vb%nz=ndim(3)
    ! Make space
    allocate(vb%box( vb%nx,vb%ny,vb%nz ))

    ! Reset the counter
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        vb%box(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! Count members per box
    do l=1,n
        v0=r(:,l)
        i=floor( (v0(1)-vb%rmin(1))*vb%ird(1) )+1
        j=floor( (v0(2)-vb%rmin(2))*vb%ird(2) )+1
        k=floor( (v0(3)-vb%rmin(3))*vb%ird(3) )+1
        vb%box(i,j,k)%n=vb%box(i,j,k)%n+1
    enddo
    ! Make some space and reset counter
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        if ( vb%box(i,j,k)%n .gt. 0 ) then
            call mem%allocate(vb%box(i,j,k)%ind,vb%box(i,j,k)%n,persistent=.false.,scalable=.false.)
        endif
        vb%box(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! Store members
    do l=1,n
        v0=r(:,l)
        i=floor( (v0(1)-vb%rmin(1))*vb%ird(1) )+1
        j=floor( (v0(2)-vb%rmin(2))*vb%ird(2) )+1
        k=floor( (v0(3)-vb%rmin(3))*vb%ird(3) )+1
        vb%box(i,j,k)%n=vb%box(i,j,k)%n+1
        vb%box(i,j,k)%ind( vb%box(i,j,k)%n )=l
    enddo
    ! Check that everything got stored
    l=0
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        l=l+vb%box(i,j,k)%n
    enddo
    enddo
    enddo
    if ( l .ne. n ) then
        call rl_stop_gracefully(['Failed putting all particles in boxes.'],rl_exitcode_symmetry)
    endif
end subroutine

end module rlsy_verletlist
