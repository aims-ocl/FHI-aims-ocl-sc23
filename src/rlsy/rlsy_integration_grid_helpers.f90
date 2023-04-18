module rlsy_integration_grid_helpers
!!
!! The integration grid module got long, thought I should move some functions
!! here. All routines and types here are useless outside the context of
!! rlsy_integration_grid and should only be used there.
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_pi,rl_twopi,rl_sqtol,rl_exitcode_mpi,&
    rl_exitcode_symmetry,rl_exitcode_memory
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: tochar,rl_sqnorm,norm2
use rlsy_sorting, only: rl_qsort
use rlsy_linalg, only: rl_symmetric_eigensystem_3x3matrix
use rlsy_mpi_helper, only: rl_mpi_helper,MPI_INTEGER,MPI_DOUBLE_PRECISION,rl_stop_gracefully
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list,rl_alltoall_distribute_2d_real,rl_alltoall_distribute_2d_int

implicit none
private
public :: rl_chunking_helper
public :: rl_set_of_points
public :: rl_distribute_sets
public :: rl_split_set
public :: rl_get_irr_points_per_rank

!> helper type to keep track of how the division into smaller pieces is going when I do it in parallel.
type rl_chunking_helper
    !> local number of members of this chunk
    integer :: n_local=-rl_hugeint
    !> global number of members of this chunk
    integer :: n_global=-rl_hugeint
    !> center of mass
    real(r8), dimension(3) :: com=-rl_huge
    !> bounding sphere from the center of mass
    real(r8) :: boundingsphere=-rl_huge
    !> normal of dividing plane
    real(r8), dimension(3) :: normal
    !> position of dividing plane
    real(r8) :: pos=-rl_huge
    !> safe lower bound of dividing plane position
    real(r8) :: pos_lb=-rl_huge
    !> safe upper bound of dividing plane position
    real(r8) :: pos_ub=-rl_huge
    !> verdict what is supposed to be done with this chunk
    integer :: howslice=-rl_hugeint
end type

!> helper type to assist with serial partitioning
type rl_set_of_points
    !> how many points
    integer :: n
    !> coordinates of the points
    real(r8), dimension(:,:), allocatable :: r
    !> some indices that can be useful to keep track of
    integer, dimension(:,:), allocatable :: ind
end type

!> helper type to aid with symmetry reduction
type rl_radshell
    !> which irreducible atom is it attached to?
    integer :: index_irr_atom=-rl_hugeint
    !> which atom is it attached to
    integer :: index_atom=-rl_hugeint
    !> which species does the grid come from
    integer :: index_species=-rl_hugeint
    !> which radial shell on that species did it come from
    integer :: index_radial=-rl_hugeint
    !> how many points in this shell
    integer :: n_point=-rl_hugeint
    !> how many irreducible points in this shell?
    integer :: n_irr_point=-rl_hugeint
    !> actual points, all of them
    real(r8), dimension(:,:), allocatable :: r
    !> irreducible points
    real(r8), dimension(:,:), allocatable :: ir
    !> which irreducible point is each point
    integer, dimension(:), allocatable :: irr_ind
    !> what is the original angular index of each irreducible point?
    integer, dimension(:), allocatable :: irr_angular_ind
    !> which operation takes the irreducible to the full array
    integer, dimension(:), allocatable :: irr_operation
    !> for each irreducible point, how many of the full points does it fold out to?
    integer, dimension(:), allocatable :: irr_ctr
    !> for each irreducible point, which of full points does it fold out to?
    integer, dimension(:,:), allocatable :: irr_unfold_ind
    !> for each irreducible point, which operation is used to unfold it?
    integer, dimension(:,:), allocatable :: irr_unfold_op
end type

! Helper types to speed up symmetry detection of radial grids
type rl_sphere_boxes_box
    integer :: n
    integer, dimension(:), allocatable :: ind
end type
type rl_sphere_boxes
    !> number of boxes
    integer :: n
    !> actual boxes
    type(rl_sphere_boxes_box), dimension(:), allocatable :: b
    contains
        !> sort points into boxes
        procedure :: sort_points_into_boxes
        !> get the box index from the coordinate
        procedure :: boxind
end type

contains

!> distribute sets from one rank to all the ranks in the subcommunicator
subroutine rl_distribute_sets(set,n_set,ml,mw)
    !> set of points
    type(rl_set_of_points), dimension(:), allocatable, intent(inout) :: set
    !> how many sets
    integer, intent(inout) :: n_set
    !> MPI helper. The ml instead of mw, like normal, indicate that this is a semilocal communicator.
    type(rl_mpi_helper), intent(inout) :: ml
    !> MPI helper, the usual world communicator
    type(rl_mpi_helper), intent(inout) :: mw


    real(r8), dimension(:,:), allocatable :: rbuf0,rbuf1
    integer, dimension(:,:), allocatable :: ibuf0,ibuf1,ctr_for_batch
    integer, dimension(:), allocatable :: sendcount,sendoffset
    integer :: nset_semilocal,nset_per_rank,nrow
    integer :: i,j,l,ii,jj

    ! First count the number of sets in this chunk.
    allocate(sendcount(ml%n))
    sendcount=0
    sendcount(ml%r+1)=n_set
    call ml%allreduce('sum',sendcount)
    nset_semilocal=sum(sendcount)

    ! Sanity test: the design should be that there is a constant number
    ! of batches per rank. Check that it is really the case.
    if ( mod(nset_semilocal,ml%n) .ne. 0 ) then
        call rl_stop_gracefully(['The points did not get divided as one would expect. Stopping now.'],rl_exitcode_symmetry,mw%comm)
    endif
    ! There should also only be one rank per chunk that holds any points
    if ( count(sendcount/=0) .ne. 1 ) then
        call rl_stop_gracefully(['The points did not get divided as one would expect. Stopping now.'],rl_exitcode_symmetry,mw%comm)
    endif
    deallocate(sendcount)

    ! Get a count for number of points per set and rank.
    nset_per_rank=nset_semilocal/ml%n
    allocate(ctr_for_batch(nset_per_rank,ml%n))
    if ( n_set .gt. 0 ) then
        l=0
        do i=1,ml%n
        do j=1,nset_per_rank
            l=l+1
            ctr_for_batch(j,i)=set(l)%n
        enddo
        enddo
    else
        ctr_for_batch=0
    endif
    call ml%bcast(ctr_for_batch,0)

    ! I also need to know the number or rows in the index array on all ranks.
    if ( n_set .gt. 0 ) then
        nrow=size(set(1)%ind,1)
    else
        nrow=0
    endif
    call ml%bcast(nrow,0)

    ! That made things easier, I think. I will just divide the sets in the easiest possible way.
    ! I will likely permute them afterwards in any case. First we have to move the sets into a
    ! large array for communication.
    if ( n_set .gt. 0 ) then
        allocate(sendcount(ml%n))
        allocate(sendoffset(ml%n))
        sendcount=0
        sendoffset=0

        ! Count total number of points
        j=sum(ctr_for_batch)
        ! Space for send buffers
        allocate(rbuf0(3,j))
        allocate(ibuf0(nrow,j))
        ! Store things in buffer
        l=0
        do i=1,n_set
            do j=1,set(i)%n
                l=l+1
                rbuf0(:,l)=set(i)%r(:,j)
                ibuf0(:,l)=set(i)%ind(:,j)
            enddo
        enddo
        ! Figure out how many things to send to each rank
        ii=0
        do i=1,ml%n
            jj=sum(ctr_for_batch(:,i))
            sendcount(i)=jj
            sendoffset(i)=ii
            ii=ii+jj
        enddo
    else
        allocate(sendcount(1))
        allocate(sendoffset(1))
        allocate(rbuf0(1,1))
        allocate(ibuf0(1,1))
        sendcount=-rl_hugeint
        sendoffset=-rl_hugeint
        rbuf0=-rl_huge
        ibuf0=-rl_hugeint
    endif

    ! Space for recieve buffers
    j=sum(ctr_for_batch(:,ml%r+1))
    allocate(rbuf1(3,j))
    allocate(ibuf1(nrow,j))
    rbuf1=0.0_r8
    ibuf1=0

    ! Distribute the data!
    call MPI_Scatterv(rbuf0, sendcount*3, sendoffset*3,&
                      MPI_DOUBLE_PRECISION, rbuf1, size(rbuf1),&
                      MPI_DOUBLE_PRECISION, 0, ml%comm, ml%error)
    if ( ml%error .ne. 0 ) then
        call rl_stop_gracefully(['mpi_scatterv exit code '//tochar(ml%error)],rl_exitcode_mpi,mw%comm)
    endif

    call MPI_Scatterv(ibuf0, sendcount*nrow, sendoffset*nrow,&
                      MPI_INTEGER, ibuf1, size(ibuf1),&
                      MPI_INTEGER, 0, ml%comm, ml%error)
    if ( ml%error .ne. 0 ) then
        call rl_stop_gracefully(['mpi_scatterv exit code '//tochar(ml%error)],rl_exitcode_mpi,mw%comm)
    endif

    ! Partial cleanup
    deallocate(rbuf0)
    deallocate(ibuf0)
    deallocate(sendcount)
    deallocate(sendoffset)
    ! Destroy existing sets
    do i=1,n_set
        deallocate(set(i)%r)
        deallocate(set(i)%ind)
    enddo
    if ( allocated(set) ) deallocate(set)
    n_set=0

    ! Allocate sets again, and populate
    n_set=nset_per_rank
    allocate(set(nset_per_rank))
    l=0
    do i=1,n_set
        set(i)%n=ctr_for_batch(i,ml%r+1)
        allocate(set(i)%r(3,set(i)%n))
        allocate(set(i)%ind(nrow,set(i)%n))
        set(i)%r=0.0_r8
        set(i)%ind=0
        do j=1,set(i)%n
            l=l+1
            set(i)%r(:,j)=rbuf1(:,l)
            set(i)%ind(:,j)=ibuf1(:,l)
        enddo
    enddo

    ! And final cleanup
    deallocate(rbuf1)
    deallocate(ibuf1)
    deallocate(ctr_for_batch)
end subroutine

!> Split a set of points into several sets
subroutine rl_split_set(set,is,nset,splitcount)
    !> all sets of points
    type(rl_set_of_points), dimension(:), intent(inout) :: set
    !> which is the current set we are working on
    integer, intent(in) :: is
    !> How many sets do I currently have
    integer, intent(in) :: nset
    !> How should it be split
    integer, dimension(:), intent(in) :: splitcount

    real(r8), dimension(:,:), allocatable :: dr
    real(r8), dimension(:), allocatable :: dv0,dv1
    real(r8), dimension(3,3) :: m0,m1
    real(r8), dimension(3) :: com,v0,normal
    real(r8) :: f0
    integer, dimension(:,:), allocatable :: di
    integer, dimension(:), allocatable :: dj
    integer :: i,j,l,js,im

    ! Some dummy space needed
    allocate(di(size(set(is)%ind,1),set(is)%n))
    allocate(dj(set(is)%n))
    allocate(dr(3,set(is)%n))
    allocate(dv0(set(is)%n))
    allocate(dv1(set(is)%n))
    di=0
    dj=0
    dr=0.0_r8
    dv0=0.0_r8
    dv1=0.0_r8

    ! Calculate the center of mass
    com=0.0_r8
    do i=1,set(is)%n
        com=com+set(is)%r(:,i)
    enddo
    com=com/real(set(is)%n,r8)

    ! Get the normal of the plane that passes through the center of mass.
    ! That will split things as evenly as possible. First get the vectors
    ! that point from the center of mass to the grid.
    do i=1,set(is)%n
        dr(:,i)=set(is)%r(:,i)-com
    enddo
    ! calculate m0 = dr0*dr0^T
    m0=0.0_r8
    call dsyrk('U','N',3,set(is)%n,1.0_r8,dr,3,0.0_r8,m0,3)
    ! Fill out the missing pieces of the matrix, don't think
    ! it's necessary but does not hurt.
    m0(2,1)=m0(1,2)
    m0(3,1)=m0(1,3)
    m0(3,2)=m0(2,3)
    ! get eigenvalues
    call rl_symmetric_eigensystem_3x3matrix(m0,v0,m1)
    ! set the normal to the eigenvector corresponding to largest eigenvalue
    f0=-rl_huge
    do i=1,3
        if ( v0(i) .gt. f0 ) then
            f0=v0(i)
            normal=m1(:,i)/norm2(m1(:,i))
        endif
    enddo

    ! calculate distances to this plane
    f0=-dot_product(normal,com)
    do i=1,set(is)%n
        dv0(i)=dot_product(normal,set(is)%r(:,i))+f0
    enddo

    ! sort the distances
    dv1=dv0
    call rl_qsort(dv1,dj)

    ! Keep all the points in a temporary place
    dr=set(is)%r
    di=set(is)%ind
    ! Destroy current points
    deallocate(set(is)%r)
    deallocate(set(is)%ind)

    ! Store the new sets.
    l=0
    do im=1,size(splitcount)
        ! What is the index of the new set?
        js=is+(im-1)*nset
        set(js)%n=splitcount(im)
        allocate(set(js)%r(3,splitcount(im)))
        allocate(set(js)%ind(size(di,1),splitcount(im)))
        do i=1,splitcount(im)
            j=dj(l+i)
            set(js)%r(:,i)=dr(:,j)
            set(js)%ind(:,i)=di(:,j)
        enddo
        l=l+splitcount(im)
    enddo

    ! And a little cleanup
    deallocate(dr)
    deallocate(dv0)
    deallocate(dv1)
    deallocate(di)
    deallocate(dj)
end subroutine

!> sort radial shells into boxes
subroutine sort_points_into_boxes(box,r,n)
    !> box helper
    class(rl_sphere_boxes), intent(out) :: box
    !> points on the unit sphere
    real(r8), dimension(:,:), intent(in) :: r
    !> how many boxes to construct
    integer, intent(in) :: n

    integer :: i,j

    box%n=n
    allocate(box%b(box%n*box%n))
    do i=1,box%n**2
        box%b(i)%n=0
    enddo

    do i=1,size(r,2)
        j=box%boxind(r(:,i))
        box%b(j)%n=box%b(j)%n+1
    enddo
    do i=1,box%n**2
        if ( box%b(i)%n .gt. 0 ) then
            allocate(box%b(i)%ind( box%b(i)%n ) )
            box%b(i)%ind=0
        endif
        box%b(i)%n=0
    enddo
    do i=1,size(r,2)
        j=box%boxind(r(:,i))
        box%b(j)%n=box%b(j)%n+1
        box%b(j)%ind( box%b(j)%n )=i
    enddo
end subroutine

!> get the box index from coordinate
function boxind(box,r) result(ind)
    !> box helper
    class(rl_sphere_boxes), intent(in) :: box
    !> point
    real(r8), dimension(3), intent(in) :: r
    !> box index
    integer :: ind

    real(r8), parameter :: sh=1E-7_r8
    real(r8), parameter :: iptheta=1.0_r8/(rl_pi+2*sh)
    real(r8), parameter :: ipphi=1.0_r8/(rl_twopi+2*sh)
    real(r8) :: theta,phi
    integer :: i,j

    theta=acos(r(3))
    phi=atan2(r(2),r(1))+rl_pi
    i=floor((theta+sh)*box%n*iptheta)+1
    j=floor((phi+sh)*box%n*ipphi)+1
    ! These should never trigger, I think
    i=min(i,box%n)
    i=max(i,1)
    j=min(j,box%n)
    j=max(j,1)
    ind=(i-1)*box%n+j
end function

!> Reduce the integration grids with symmetry, and return a list of points per rank to keep working on
subroutine rl_get_irr_points_per_rank(p,sym,mw,mem,n_radial,n_angular,r_angular,&
    n_irr_point_local,irr_ind)
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> symmetry handle
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> grid information from AIMS
    integer, dimension(:), intent(in) :: n_radial
    integer, dimension(:,:), intent(in) :: n_angular
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    !> how many irreducible points do I have on this rank
    integer, intent(out) :: n_irr_point_local
    !> index arrays that help map out symmetry
    integer, dimension(:,:), allocatable, intent(out) :: irr_ind


    ! Helper container for the irreducible shells that don't take up too much space.
    ! Have to be really careful with memory this early in the division.
    integer :: nshell_local
    type(rl_radshell), dimension(:), allocatable :: sh
    !init: block
        integer :: a1,i1,s1,ish,ctr
        integer :: i,j,k
    !reduce: block
        integer, parameter :: prime1=19,prime2=17
        type(rl_sphere_boxes) :: bx1,bx2
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(3) :: v0,v1,w
        real(r8) :: f0
        integer, dimension(:), allocatable :: dj !,safeoperations
        integer, dimension(:,:), allocatable :: di
        integer, dimension(sym%n_operation) :: ops,ok_op
        !integer :: ish,nop,a1,ctr_ok_op
        integer :: nop,ctr_ok_op
        !integer :: i,j,k,l,o,oo,ii
        integer :: l,o,oo,ii
    !packanddistribute: block
        !integer :: ish,maxmul,nrowind
        integer :: maxmul,nrowind
        integer :: iirr,iatm,irad,iang,ispc
        !integer :: i,j,l

    ! count some basic things
    !init: block

        ! Count how many shells should be reduced with symmetry, on all ranks and on this rank
        ctr=0
        nshell_local=0
        do i1=1,sym%n_irreducible_atom
            a1=sym%irr_to_all(i1)
            s1=p%species(a1)
            do k=1,n_radial(s1)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .eq. mw%r ) nshell_local=nshell_local+1
            enddo
        enddo

        ! Start populating the shells
        allocate(sh(nshell_local))
        ish=0
        ctr=0
        do i1=1,sym%n_irreducible_atom
            a1=sym%irr_to_all(i1)
            s1=p%species(a1)
            do i=1,n_radial(s1)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                ish=ish+1
                ! Store away this shell. First some helpers to keep
                ! track of where it comes from so that things can be
                ! reconstructed later, if needed.
                sh(ish)%index_irr_atom = i1
                sh(ish)%index_atom = a1
                sh(ish)%index_species = s1
                sh(ish)%index_radial = i
                sh(ish)%n_point = n_angular(i,s1)
                ! Radius of this shell
                !sh(ish)%shell_radius = r_radial(i,s1)
                ! Then the actual points. It is important that the coordinates
                ! are not downfolded to the irreducible, and that they have
                ! norm 1 for now, to make trigonometrics faster in the symmetry stuff.
                call mem%allocate(sh(ish)%r,[3,sh(ish)%n_point],persistent=.true.,scalable=.true.)
                do j=1,n_angular(i,s1)
                    sh(ish)%r(:,j)=r_angular(:,j,i,s1)
                enddo
            enddo
        enddo
    !end block init

    ! So, at this point I have created all the shells of all the unique
    ! atoms, with the shells distributed evenly across ranks. I have to
    ! distribute them because it can be quite many points. Now I can
    ! proceed serially, to reduce the local shells of points with symmetry.
    !reduce: block

        shloop: do ish=1,nshell_local
            ! Symmetry reduction essentially boil down to taking a point in the shell,
            ! rotating it with a 3x3-matrix, and see if the result is some other point
            ! in the shell. To avoid N^2 scaling I do some hashing of the points into
            ! boxes to make lookups fast.

            ! But first I have to decide on the pool of possible operations to test with,
            ! the mesh should only be reduced with the symmetry operations that leave
            ! this atom invariant.
            ops=-1
            nop=0
            a1=sh(ish)%index_atom
            do o=1,sym%n_operation
                if ( sym%op(o)%fmap(a1) .eq. a1 ) then
                    nop=nop+1
                    ops(nop)=o
                endif
            enddo

            !@TODO insert fast exit here in case there is just the identity.

            ! Now put the points into boxes to facilitate fast lookup.
            ! I use two boxes, that way I don't have to look in neighbouring
            ! boxes and I don't have to deal with annoying tolerances.
            call bx1%sort_points_into_boxes(sh(ish)%r,prime1)
            call bx2%sort_points_into_boxes(sh(ish)%r,prime2)

            ! Some temporary space
            call mem%allocate(dr,[3,sh(ish)%n_point],persistent=.false.,scalable=.true.)
            call mem%allocate(di,[sh(ish)%n_point,nop],persistent=.false.,scalable=.true.)
            call mem%allocate(dj,sh(ish)%n_point,persistent=.false.,scalable=.true.)
            dr=0.0_r8
            di=0
            dj=0
            ! Now start checking with respect to all points.
            ok_op=-1
            ctr_ok_op=0
            oploop1: do oo=1,nop
                dj=0
                o=ops(oo)
                ! Rotate the points
                call dgemm('N','N',3,sh(ish)%n_point,3,1.0_r8,sym%op(o)%m,3,sh(ish)%r,3,0.0_r8,dr,3)
                ! See where they end up
                l=0
                il1: do i=1,sh(ish)%n_point
                    ii=0
                    ! check first box
                    k=bx1%boxind( dr(:,i) )
                    do j=1,bx1%b(k)%n
                        w=sh(ish)%r(:,bx1%b(k)%ind(j) )-dr(:,i)
                        f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
                        if ( f0 .lt. rl_sqtol ) then
                            l=l+1
                            ii=ii+1
                            dj(i)=bx1%b(k)%ind(j)
                            cycle il1
                        endif
                    enddo
                    ! maybe check second box, should be very unlikely.
                    if ( ii .eq. 0 ) then
                        k=bx2%boxind( dr(:,i) )
                        do j=1,bx2%b(k)%n
                            w=sh(ish)%r(:,bx2%b(k)%ind(j) )-dr(:,i)
                            f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
                            if ( f0 .lt. rl_sqtol ) then
                                l=l+1
                                ii=ii+1
                                dj(i)=bx2%b(k)%ind(j)
                                cycle il1
                            endif
                        enddo
                        ! now if we still did not find it, this is not a good operation. Therefore
                        ! we can kill the search early. With good grids it will never happen.
                        if ( ii .eq. 0 ) cycle oploop1
                    endif
                enddo il1

                ! In case all went well, make note of that.
                if ( l .eq. sh(ish)%n_point ) then
                    ! Yep, this was a good operation
                    ctr_ok_op=ctr_ok_op+1
                    ok_op(ctr_ok_op)=o
                    di(:,ctr_ok_op)=dj
                else
                    ! This should be impossible.
                    call rl_stop_gracefully(['Clearly I do not understand symmetry'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo oploop1
            call mem%deallocate(dj,persistent=.false.,scalable=.true.)
            call mem%deallocate(dr,persistent=.false.,scalable=.true.)

            ! Destroy the boxes
            do i=1,bx1%n
                if ( allocated(bx1%b(i)%ind) ) deallocate(bx1%b(i)%ind)
            enddo
            deallocate(bx1%b)
            do i=1,bx2%n
                if ( allocated(bx2%b(i)%ind) ) deallocate(bx2%b(i)%ind)
            enddo
            deallocate(bx2%b)

            ! Note to self: once I have proper symmetric grids, I should never need
            ! the ctr_ok_op, they should all be ok. Will fix eventually.

            ! Now I have a permutation array index-wise, in the sense
            ! that I know what action some of the symmetry operations
            ! have on the shell of points, i.e. how they become permuted.
            ! And cooking-show style, I know how to sort the symmetries
            ! from an array like that, since I already prepared it!
            call rl_coordination_shells_from_permutation_list(di(:,1:ctr_ok_op),sh(ish)%irr_ctr,sh(ish)%irr_unfold_ind,sh(ish)%irr_ind,mem,dj)
            call mem%deallocate(di,persistent=.false.,scalable=.true.)

            ! Number of irreducible points
            sh(ish)%n_irr_point = size(sh(ish)%irr_ctr)
            ! Store the coordinates of the irreducible points to avoid having to ninja-index everything.
            allocate(sh(ish)%ir(3,sh(ish)%n_irr_point))
            sh(ish)%ir=0.0_r8
            do i=1,sh(ish)%n_irr_point
                j=dj(i)
                sh(ish)%ir(:,i)=sh(ish)%r(:,j)
            enddo

            ! I did not, however, make it supergeneral so we still have to
            ! sort out which operation does what.
            call mem%allocate(sh(ish)%irr_operation  ,sh(ish)%n_point,persistent=.true.,scalable=.true.)
            call mem%allocate(sh(ish)%irr_unfold_op  ,[size(sh(ish)%irr_unfold_ind,1), sh(ish)%n_irr_point],persistent=.true.,scalable=.true.)
            call mem%allocate(sh(ish)%irr_angular_ind,sh(ish)%n_irr_point,persistent=.true.,scalable=.true.)
            sh(ish)%irr_operation=0
            sh(ish)%irr_unfold_op=0
            sh(ish)%irr_angular_ind=dj

            ! I don't think this takes too long, so I can do it kinda stupid I think.
            ! Also this serves as a sensible sanity test, so that's always something.
            ! Anyway, this gets you which operation takes the irreducible to the full.
            ! This is used when you stand at a full point and want to fetch data from
            ! the irreducible.
            do i=1,sh(ish)%n_point
                j=sh(ish)%irr_ind(i)
                v0=sh(ish)%r(:,i)
                v1=sh(ish)%ir(:,j)
                l=0
                do oo=1,ctr_ok_op
                    o=ok_op(oo)
                    w=v0-matmul(sym%op(o)%m,v1)
                    if ( rl_sqnorm(w) .lt. rl_sqtol ) then
                        l=l+1
                        sh(ish)%irr_operation(i)=o
                        exit
                    endif
                enddo
                ! Check that I was not an idiot.
                if ( l .ne. 1 ) then
                    call rl_stop_gracefully(['Clearly I do not understand symmetry'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo

            ! And now the opposite, the data you need if you stand at an irreducible point and want to
            ! know which full points are equivalent.
            do i=1,sh(ish)%n_irr_point
                v0=sh(ish)%ir(:,i)
                do j=1,sh(ish)%irr_ctr(i)
                    k=sh(ish)%irr_unfold_ind(j,i)
                    v1=sh(ish)%r(:,k)
                    o=sh(ish)%irr_operation(k)
                    if ( rl_sqnorm(matmul(sym%op(o)%m,v0)-v1) .lt. rl_sqtol ) then
                        sh(ish)%irr_unfold_op(j,i)=o
                    else
                        call rl_stop_gracefully(['Clearly I do not understand symmetry'],rl_exitcode_symmetry,mw%comm)
                    endif
                enddo
            enddo

            ! some cleanup
            call mem%deallocate(dj,persistent=.true.,scalable=.false.)
        enddo shloop
    !end block reduce

    ! Now pack the points together, and distribute them evenly across ranks.
    !packanddistribute: block

        ! First count number of points on this rank, and the largest
        ! multiplicity of any points
        n_irr_point_local=0
        maxmul=0
        do ish=1,nshell_local
            n_irr_point_local=n_irr_point_local+sh(ish)%n_irr_point
            maxmul=max(maxmul,maxval(sh(ish)%irr_ctr))
        enddo
        call mw%allreduce('max',maxmul)
        nrowind=maxmul*2+5

        ! Make some space to store points
        if ( n_irr_point_local .gt. 0 ) then
            allocate(irr_ind(nrowind,n_irr_point_local))
            irr_ind=0
        else
            allocate(irr_ind(1,1))
            irr_ind=-rl_hugeint
        endif

        ! Go through all shells and store the symmetry information.
        l=0
        do ish=1,nshell_local
            ! many indices are constant per shell
            iirr=sh(ish)%index_irr_atom
            iatm=sym%irr_to_all(iirr)
            ispc=p%species(iatm)
            irad=sh(ish)%index_radial
            do i=1,sh(ish)%n_irr_point
                l=l+1
                iang=sh(ish)%irr_angular_ind(i)
                irr_ind(1,l)=iirr
                irr_ind(2,l)=iatm
                irr_ind(3,l)=irad
                irr_ind(4,l)=iang
                irr_ind(5,l)=sh(ish)%irr_ctr(i)
                do j=1,sh(ish)%irr_ctr(i)
                    irr_ind(5+j,l)=sh(ish)%irr_unfold_op(j,i)
                    irr_ind(5+sh(ish)%irr_ctr(i)+j,l)=sh(ish)%irr_unfold_ind(j,i)
                enddo
            enddo
        enddo

        ! Distribute the points more evenly across ranks.
        call rl_alltoall_distribute_2d_int(irr_ind,nrowind,n_irr_point_local,mw)

        ! And cleanup
        do ish=1,nshell_local
            deallocate(sh(ish)%ir)
            deallocate(sh(ish)%irr_ind)
            call mem%deallocate(sh(ish)%r,persistent=.true.,scalable=.true.)
            call mem%deallocate(sh(ish)%irr_operation,persistent=.true.,scalable=.true.)
            call mem%deallocate(sh(ish)%irr_unfold_op,persistent=.true.,scalable=.true.)
            call mem%deallocate(sh(ish)%irr_angular_ind,persistent=.true.,scalable=.true.)
            call mem%deallocate(sh(ish)%irr_ctr,persistent=.true.,scalable=.false.)
            call mem%deallocate(sh(ish)%irr_unfold_ind,persistent=.true.,scalable=.false.)
        enddo
        if ( allocated(sh) ) deallocate(sh)
    !end block packanddistribute

    ! Check that everything is cleared properly. Might remove this check and move it outside later.
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
end subroutine

end module
