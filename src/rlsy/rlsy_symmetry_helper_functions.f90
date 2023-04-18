module rlsy_symmetry_helper_functions
!!
!! OH: Some quite general helper functions that are used every now and then.
!! I put the functions here once I use them in more than one place.
!!
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_tol,rl_sqtol,&
    rl_exitcode_symmetry,rl_status,rl_exitcode_mpi,rl_exitcode_param
use rlsy_memtracker, only: rl_memtracker
use rlsy_sorting, only: rl_qsort,rl_return_unique,rl_return_unique_indices
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,MPI_DOUBLE_PRECISION,MPI_INTEGER

implicit none
private
public :: rl_coordination_shells_from_permutation_list
public :: rl_alltoall_distribute_2d_real
public :: rl_alltoall_distribute_2d_int
public :: rl_find_prime_factors
public :: rl_find_divisors
public :: rl_distribute_weighted

contains

!> distribute a set of weighted things across ranks
subroutine rl_distribute_weighted(mw,weight,n_elem,elem,mem)
    !> MPI communicator to distribute across
    type(rl_mpi_helper), intent(in) :: mw
    !> weight of things I am to distribute
    integer, dimension(:), intent(in) :: weight
    !> how many elements does this rank get?
    integer, intent(out) :: n_elem
    !> which elements does this rank get
    integer, dimension(:), allocatable, intent(out) :: elem
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    integer, dimension(:), allocatable :: buf_wt,buf_ind
    integer :: n,i,j,ihi,ilo,iter,niter

    ! The purpose of this routine is to distribute things over MPI ranks. Each of the things has a certain
    ! weight, or cost, associated with it. We want to ensure that we get a distribution of sums of weights
    ! over the ranks that is about as constant as possible. Pretty sure it is an NP-hard problem in general,
    ! e.g. box-packing, so I'm just going for something that is not horrible. If load-balancing becomes an
    ! issue eventually, I would have to revise.

    ! Number of things we have to distribute.
    n=size(weight)
    ! We need a little temporary space
    call mem%allocate(buf_wt,n,persistent=.false.,scalable=.false.)
    call mem%allocate(buf_ind,n,persistent=.false.,scalable=.false.)
    buf_wt=0
    buf_ind=0

    ! Start by sorting the things by weight
    buf_wt=weight
    call rl_qsort(buf_wt,buf_ind)

    ! Now count things per rank.
    n_elem=0
    do i=1,n
        if ( mod(i,mw%n) .eq. mw%r ) n_elem=n_elem+1
    enddo
    if ( n_elem .gt. 0 ) then
        call mem%allocate(elem,n_elem,persistent=.true.,scalable=.true.)
        elem=0
    else
        ! Dummy allocate, and return early in case there are no elements on this rank.
        call mem%allocate(elem,1,persistent=.true.,scalable=.true.)
        elem=-rl_hugeint
        return
    endif

    ! A little convoluted way of looping over values, but should help with load balancing.
    niter=0
    do
        if ( niter*mw%n .gt. n ) then
            exit
        else
            niter=niter+1
        endif
    enddo

    ! So I go through the sorted values, and alternate direction every nrank steps.
    ! that should neatly fill up the ranks in a somewhat evenly distributed manner.
    j=0
    n_elem=0
    do iter=1,niter
        ilo=(iter-1)*mw%n+1
        ihi=min(iter*mw%n,n)
        if ( mod(iter,2) .eq. 0 ) then
            do i=ihi,ilo,-1
                j=j+1
                if ( mod(j,mw%n) .eq. mw%r ) then
                    n_elem=n_elem+1
                    elem(n_elem)=buf_ind(i)
                endif
            enddo
        else
            do i=ilo,ihi
                j=j+1
                if ( mod(j,mw%n) .eq. mw%r ) then
                    n_elem=n_elem+1
                    elem(n_elem)=buf_ind(i)
                endif
            enddo
        endif
    enddo

    ! Cleanup
    call mem%deallocate(buf_wt,persistent=.false.,scalable=.false.)
    call mem%deallocate(buf_ind,persistent=.false.,scalable=.false.)
end subroutine

!> Find a neat set of divisors to distribute stuff over two levels.
subroutine rl_find_divisors(nrank,npts,nchunk,n_pts_per_batch,n_pts_per_chunk,n_divisors_per_chunk,divisors_per_chunk)
    !> total number of ranks
    integer, intent(in) :: nrank
    !> total number of points
    integer, intent(in) :: npts
    !> total number of chunks
    integer, intent(in) :: nchunk
    !> desired number of points per batch
    integer, intent(in) :: n_pts_per_batch
    !> number of points per chunk
    integer, dimension(:), allocatable, intent(out) :: n_pts_per_chunk
    !> number of divisors per chunk
    integer, dimension(:), allocatable :: n_divisors_per_chunk
    !> divisors per chunk
    integer, dimension(:,:), allocatable :: divisors_per_chunk

    integer, parameter :: max_n_divisors=256 ! Enough for 10^77 points.
    real(r8) :: f0,f1,f2,f3
    integer, dimension(:), allocatable :: n_ranks_per_chunk
    integer, dimension(:), allocatable :: di,dj,dk
    integer, dimension(4) :: listofprimes
    integer :: i,j,k,l
    integer :: ii !,jj,kk

    ! Some primes I will try to use for partitioning.
    listofprimes=[1,2,3,5]

    ! Get the number of ranks per chunk
    allocate(n_ranks_per_chunk(nchunk))
    n_ranks_per_chunk=0
    do i=1,nrank
        j=mod(i,nchunk)+1
        n_ranks_per_chunk(j)=n_ranks_per_chunk(j)+1
    enddo

    ! I can get the number of points/chunk. This should be
    ! somewhat reasonable, I think, and give a decent standard
    ! deviation of the number or points/batch.
    allocate(n_pts_per_chunk(nchunk))
    n_pts_per_chunk=0
    j=0
    do i=1,nchunk
        if ( i .lt. nchunk ) then
            f0=real(npts,r8)/real(nrank,r8)     ! points/rank
            f0=f0*n_ranks_per_chunk(i)          ! ideal points/chunk
            if ( mod(i,2) .eq. 0 ) then
                n_pts_per_chunk(i)=ceiling(f0)
            else
                n_pts_per_chunk(i)=floor(f0)
            endif
            j=j+n_pts_per_chunk(i)
        else
            n_pts_per_chunk(i)=npts-j           ! the rest of the points
        endif
    enddo

    ! Ok, decent start. Now I need a series of divisors such
    ! that I get a constant number of batches per rank, with
    ! somewhat similar size. Start with the requirement that
    ! we get at least one batch per rank.
    allocate(n_divisors_per_chunk(nchunk))
    allocate(divisors_per_chunk(max_n_divisors,nchunk))
    divisors_per_chunk=0
    n_divisors_per_chunk=0

    do i=1,nchunk
        call rl_find_prime_factors(n_ranks_per_chunk(i),di)
        n_divisors_per_chunk(i)=size(di)
        divisors_per_chunk(1:size(di),i)=di
        deallocate(di)
    enddo

    ! Now try a few different ways to divide it, and see how close
    ! I can get to the target number of points per batch.
    allocate(di(max_n_divisors))
    allocate(dj(max_n_divisors))
    allocate(dk(max_n_divisors))
    di=0
    dj=0
    ! First we divide out the baseline of what is there.
    f0=0.0_r8
    do i=1,nchunk
        f0=f0+real(n_pts_per_chunk(i),r8)/product(divisors_per_chunk(1:n_divisors_per_chunk(i),i))
    enddo
    f0=f0/nchunk
    ! f0 now holds the current number of points/rank

    ! Find sensible series of divisions that give something close to the desired number of points.
    ! Can be made more flexible, at the cost of dividing by larger primes. Not sure what is best.
    f2=rl_huge
    f3=rl_huge
    j=0
    k=0
    do ii=1,3
    !do jj=1,3
    !do kk=1,4
        di=2
        di(1)=listofprimes(ii)
        !di(2)=listofprimes(jj)
        !di(3)=listofprimes(kk)
        do i=1,max_n_divisors
            f1=f0/product(di(1:i))
            if ( abs(f1-n_pts_per_batch) .lt. f2 ) then
                f2=abs(f1-n_pts_per_batch)
                dj=di
                j=i
            endif
            if ( abs(f1-n_pts_per_batch) .lt. f3 .and. f1 .gt. real(n_pts_per_batch,r8) ) then
                f3=abs(f1-n_pts_per_batch)
                dk=di
                k=i
            endif
            if ( f1 .lt. n_pts_per_batch*0.75_r8 ) exit
        enddo
    enddo
    !enddo
    !enddo

    ! Pick the one above the minimum if available, if not pick the closest.
    di=0
    if ( k .gt. 0 ) then
        ii=k
        di(1:ii)=dk(1:ii)
    else
        ii=j
        di(1:ii)=dj(1:ii)
    endif

    ! Now add this to the divisors.
    do i=1,nchunk
        l=0
        dj=0
        do j=1,n_divisors_per_chunk(i)
            l=l+1
            dj(l)=divisors_per_chunk(j,i)
        enddo
        do j=1,ii
            if ( di(j) .eq. 1 ) cycle
            l=l+1
            dj(l)=di(j)
        enddo
        ! Sort so that the large divisions are first.
        dj(1:l)=-dj(1:l)
        call rl_qsort(dj(1:l))
        dj(1:l)=-dj(1:l)
        ! Store divisions
        n_divisors_per_chunk(i)=l
        divisors_per_chunk(1:l,i)=dj(1:l)
    enddo

    ! And a little cleanup.
    deallocate(di)
    deallocate(dj)
    deallocate(dk)
    deallocate(n_ranks_per_chunk)
end subroutine

!> redistribute columns of a large array across mpi ranks. The starting point is that it is already distributed, but badly.
subroutine rl_alltoall_distribute_2d_real(A,nr,n,mw)
    !> local part of the large array
    real(r8), dimension(:,:), allocatable, intent(inout) :: A
    !> number of rows in the matrix
    integer, intent(in) :: nr
    !> how many columns on this rank? could be zero
    integer, intent(inout) :: n
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw

    real(r8), dimension(:,:), allocatable :: buf
    integer, dimension(:), allocatable :: di,dj,dk
    integer, dimension(:), allocatable :: orig_count,new_count
    integer, dimension(:), allocatable :: sendcounts,sendoffset
    integer, dimension(:), allocatable :: recvcounts,recvoffset

    integer :: i,j,k,l,ct1,ct2,m,nc

    ! Everybody needs to know the global number of columns, both in
    ! the current distribution and in the resulting distribution.
    allocate(orig_count(mw%n))
    allocate(new_count(mw%n))
    orig_count=0
    new_count=0
    orig_count(mw%r+1)=n
    call mw%allreduce('sum',orig_count)
    ! And the total number of points
    nc=sum(orig_count)

    ! Sanity test to ensure the dimensions make sense
    if ( n .gt. 0 ) then
    if ( size(A,2) .ne. orig_count(mw%r+1) ) then
        call rl_stop_gracefully(['Inconsistent dimensions in alltoall'],rl_exitcode_param,mw%comm)
    endif
    endif

    ! Figure out how many points should end up on each rank
    new_count=0
    do i=1,nc
    do j=1,mw%n
        if ( mod(i,mw%n) .eq. j-1 ) new_count(j)=new_count(j)+1
    enddo
    enddo

    ! Number of columns on this rank in the new distribution
    ! Beware, it could be zero in some edge cases.
    m=new_count(mw%r+1)

    ! Annoying index algebra.
    if ( n .gt. 0 ) then
        allocate(di(n))
    else
        allocate(di(1))
    endif
    if ( m .gt. 0 ) then
        allocate(dj(m))
    else
        allocate(dj(1))
    endif
    di=-1
    dj=-1

    l=0 ! semi-global counter
    k=0 ! counter which rank this global point should end up on
    do i=1,mw%n
    do j=1,orig_count(i)
        l=l+1
        if ( i-1 .eq. mw%r ) di(j)=k ! Where should this point be sent?
        if ( k .eq. mw%r )   dj(l)=i-1 ! Where should this point be recieved from
        if ( l .eq. new_count(k+1) ) then
            k=k+1
            l=0
        endif
    enddo
    enddo

    ! If I did my algebra correctly, the di and dj arrays
    ! should be sorted, kinda. At least elements should be
    ! grouped by rank. Maybe not a necessary test, but nice
    ! to catch things here, much harder later.
    if ( n .gt. 0 ) then
        ! Count the number of unique, two ways.
        allocate(dk(mw%n))
        dk=-1
        ct1=0
        ct2=0
        k=-1
        il1: do i=1,n
            ! First way of counting
            if ( di(i) .ne. k ) then
                ct1=ct1+1
                k=di(i)
            endif

            ! Second way of counting
            do j=1,ct2
                if ( dk(j) .eq. di(i) ) cycle il1
            enddo

            ct2=ct2+1
            if ( ct2 .gt. mw%n ) then
                call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
            else
                dk(ct2)=di(i)
            endif
        enddo il1
        if ( ct1 .ne. ct2 ) then
            call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
        endif
        deallocate(dk)
    endif
    if ( m .gt. 0 ) then
        ! Count the number of unique, two ways.
        allocate(dk(mw%n))
        dk=-1
        ct1=0
        ct2=0
        k=-1
        il2: do i=1,m
            ! First way of counting
            if ( dj(i) .ne. k ) then
                ct1=ct1+1
                k=dj(i)
            endif

            ! Second way of counting
            do j=1,ct2
                if ( dk(j) .eq. dj(i) ) cycle il2
            enddo

            ct2=ct2+1
            if ( ct2 .gt. mw%n ) then
                call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
            else
                dk(ct2)=dj(i)
            endif
        enddo il2
        if ( ct1 .ne. ct2 ) then
            call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
        endif
        deallocate(dk)
    endif

    ! Start preparing for communication
    allocate(sendcounts(mw%n))
    allocate(recvcounts(mw%n))
    allocate(sendoffset(mw%n))
    allocate(recvoffset(mw%n))
    sendcounts=0
    recvcounts=0
    sendoffset=0
    recvoffset=0

    ! I should rewrite these to be O(N) at some point.
    ! First get the send counts?
    l=0
    do i=1,mw%n
        k=0
        do j=1,n
            if ( di(j) .eq. i-1 ) k=k+1
        enddo
        sendcounts(i)=k
        sendoffset(i)=l
        l=l+k
    enddo
    ! Then the recieve counts.
    l=0
    do i=1,mw%n
        k=0
        do j=1,m
            if ( dj(j) .eq. i-1 ) k=k+1
        enddo
        recvcounts(i)=k
        recvoffset(i)=l
        l=l+k
    enddo

    if ( sum(recvcounts) .ne. m ) then
        call rl_stop_gracefully(['Could not count properly'],rl_exitcode_mpi,mw%comm)
    endif
    if ( sum(sendcounts) .ne. n ) then
        call rl_stop_gracefully(['Could not count properly'],rl_exitcode_mpi,mw%comm)
    endif

    ! Space for the recive buffer?
    if ( m .gt. 0 ) then
        allocate(buf(nr,m))
        buf=0.0_r8
    else
        ! dummy that should not get referenced
        allocate(buf(1,1))
        buf=-rl_huge
    endif

    ! Begin operation actual operation.
    call MPI_ALLTOALLV(A, sendcounts*nr, sendoffset*nr, MPI_DOUBLE_PRECISION, buf, recvcounts*nr, recvoffset*nr, MPI_DOUBLE_PRECISION, mw%comm, rl_status)
    if ( rl_status .ne. 0 ) then
        call rl_stop_gracefully(['Failed mpi_alltoallv'],rl_exitcode_mpi,mw%comm)
    endif

    ! return the results in the same array
    if ( allocated(A) ) deallocate(A)
    call move_alloc(buf,A)
    n=m

    ! And some cleanup
    deallocate(sendcounts)
    deallocate(recvcounts)
    deallocate(sendoffset)
    deallocate(recvoffset)
    deallocate(orig_count)
    deallocate(new_count)
    deallocate(di)
    deallocate(dj)
end subroutine

!> redistribute columns of a large array across mpi ranks. The starting point is that it is already distributed, but badly.
subroutine rl_alltoall_distribute_2d_int(A,nr,n,mw)
    !> local part of the large array
    integer, dimension(:,:), allocatable, intent(inout) :: A
    !> number of rows in the matrix
    integer, intent(in) :: nr
    !> how many columns on this rank? could be zero
    integer, intent(inout) :: n
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw

    integer, dimension(:,:), allocatable :: buf
    integer, dimension(:), allocatable :: di,dj,dk
    integer, dimension(:), allocatable :: orig_count,new_count
    integer, dimension(:), allocatable :: sendcounts,sendoffset
    integer, dimension(:), allocatable :: recvcounts,recvoffset

    integer :: i,j,k,l,ct1,ct2,m,nc

    ! Everybody needs to know the global number of columns, both in
    ! the current distribution and in the resulting distribution.
    allocate(orig_count(mw%n))
    allocate(new_count(mw%n))
    orig_count=0
    new_count=0
    orig_count(mw%r+1)=n
    call mw%allreduce('sum',orig_count)
    ! And the total number of points
    nc=sum(orig_count)

    ! Sanity test to ensure the dimensions make sense
    if ( n .gt. 0 ) then
        if ( size(A,2) .ne. orig_count(mw%r+1) ) then
            call rl_stop_gracefully(['Inconsistent dimensions in alltoall'],rl_exitcode_param,mw%comm)
        endif
        if ( size(A,1) .ne. nr ) then
            call rl_stop_gracefully(['Inconsistent dimensions in alltoall'],rl_exitcode_param,mw%comm)
        endif
    endif

    ! Figure out how many points should end up on each rank
    new_count=0
    do i=1,nc
    do j=1,mw%n
        if ( mod(i,mw%n) .eq. j-1 ) new_count(j)=new_count(j)+1
    enddo
    enddo

    ! Number of columns on this rank in the new distribution
    ! Beware, it could be zero in some edge cases.
    m=new_count(mw%r+1)

    ! Annoying index algebra.
    if ( n .gt. 0 ) then
        allocate(di(n))
    else
        allocate(di(1))
    endif
    if ( m .gt. 0 ) then
        allocate(dj(m))
    else
        allocate(dj(1))
    endif
    di=-1
    dj=-1

    l=0 ! semi-global counter
    k=0 ! counter which rank this global point should end up on
    do i=1,mw%n
    do j=1,orig_count(i)
        l=l+1
        if ( i-1 .eq. mw%r ) di(j)=k ! Where should this point be sent?
        if ( k .eq. mw%r )   dj(l)=i-1 ! Where should this point be recieved from
        if ( l .eq. new_count(k+1) ) then
            k=k+1
            l=0
        endif
    enddo
    enddo

    ! If I did my algebra correctly, the di and dj arrays
    ! should be sorted, kinda. At least elements should be
    ! grouped by rank. Maybe not a necessary test, but nice
    ! to catch things here, much harder later.
    if ( n .gt. 0 ) then
        ! Count the number of unique, two ways.
        allocate(dk(mw%n))
        dk=-1
        ct1=0
        ct2=0
        k=-1
        il1: do i=1,n
            ! First way of counting
            if ( di(i) .ne. k ) then
                ct1=ct1+1
                k=di(i)
            endif

            ! Second way of counting
            do j=1,ct2
                if ( dk(j) .eq. di(i) ) cycle il1
            enddo

            ct2=ct2+1
            if ( ct2 .gt. mw%n ) then
                call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
            else
                dk(ct2)=di(i)
            endif
        enddo il1
        if ( ct1 .ne. ct2 ) then
            call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
        endif
        deallocate(dk)
    endif
    if ( m .gt. 0 ) then
        ! Count the number of unique, two ways.
        allocate(dk(mw%n))
        dk=-1
        ct1=0
        ct2=0
        k=-1
        il2: do i=1,m
            ! First way of counting
            if ( dj(i) .ne. k ) then
                ct1=ct1+1
                k=dj(i)
            endif

            ! Second way of counting
            do j=1,ct2
                if ( dk(j) .eq. dj(i) ) cycle il2
            enddo

            ct2=ct2+1
            if ( ct2 .gt. mw%n ) then
                call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
            else
                dk(ct2)=dj(i)
            endif
        enddo il2
        if ( ct1 .ne. ct2 ) then
            call rl_stop_gracefully(['Could not sort properly'],rl_exitcode_mpi,mw%comm)
        endif
        deallocate(dk)
    endif

    ! Start preparing for communication
    allocate(sendcounts(mw%n))
    allocate(recvcounts(mw%n))
    allocate(sendoffset(mw%n))
    allocate(recvoffset(mw%n))
    sendcounts=0
    recvcounts=0
    sendoffset=0
    recvoffset=0

    ! I should rewrite these to be O(N) at some point.
    ! First get the send counts?
    l=0
    do i=1,mw%n
        k=0
        do j=1,n
            if ( di(j) .eq. i-1 ) k=k+1
        enddo
        sendcounts(i)=k
        sendoffset(i)=l
        l=l+k
    enddo
    ! Then the recieve counts.
    l=0
    do i=1,mw%n
        k=0
        do j=1,m
            if ( dj(j) .eq. i-1 ) k=k+1
        enddo
        recvcounts(i)=k
        recvoffset(i)=l
        l=l+k
    enddo

    if ( sum(recvcounts) .ne. m ) then
        call rl_stop_gracefully(['Could not count properly'],rl_exitcode_mpi,mw%comm)
    endif
    if ( sum(sendcounts) .ne. n ) then
        call rl_stop_gracefully(['Could not count properly'],rl_exitcode_mpi,mw%comm)
    endif

    ! Space for the recive buffer?
    if ( m .gt. 0 ) then
        allocate(buf(nr,m))
        buf=0
    else
        ! dummy that should not get referenced
        allocate(buf(1,1))
        buf=-rl_hugeint
    endif

    ! Begin operation actual operation.
    call MPI_ALLTOALLV(A, sendcounts*nr, sendoffset*nr, MPI_INTEGER, buf, recvcounts*nr, recvoffset*nr, MPI_INTEGER, mw%comm, rl_status)
    if ( rl_status .ne. 0 ) then
        call rl_stop_gracefully(['Failed mpi_alltoallv'],rl_exitcode_mpi,mw%comm)
    endif

    ! return the results in the same array
    if ( allocated(A) ) deallocate(A)
    call move_alloc(buf,A)
    n=m

    ! And some cleanup
    deallocate(sendcounts)
    deallocate(recvcounts)
    deallocate(sendoffset)
    deallocate(recvoffset)
    deallocate(orig_count)
    deallocate(new_count)
    deallocate(di)
    deallocate(dj)
end subroutine

!> build coordination shells from a list of transformation rules
subroutine rl_coordination_shells_from_permutation_list(perm,shell_ctr,shell_member,shell_index,mem,prototype_shell,mw)
    !> list of transformation rules
    integer, dimension(:,:), intent(in) :: perm
    !> counter for the number of elements in each shell
    integer, dimension(:), allocatable, intent(out) :: shell_ctr
    !> the members of each shell
    integer, dimension(:,:), allocatable, intent(out) :: shell_member
    !> which shell does each point belong to
    integer, dimension(:), allocatable, intent(out) :: shell_index
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> which are the unique points
    integer, dimension(:), allocatable, intent(out), optional :: prototype_shell
    !> mpi helper
    type(rl_mpi_helper), intent(in), optional :: mw

    integer, dimension(:,:), allocatable :: pt,pu
    integer, dimension(:), allocatable :: di,dk
    integer :: i,j,np,no,n_shell

    ! This is a fairly general routine to aid in symmetry reduction. The things
    ! that are to be reduced are given as a list of integers. The action of a
    ! symmetry operation is to permute these indices. From the list of all
    ! permutaitons, in array perm(n_things,n_operations), we can deduce what the
    ! irreducible things are, and how they should be divided into irreducible groups.

    ! Some temporary space
    np=size(perm,1) ! number of elements
    no=size(perm,2) ! number of operations
    call mem%allocate(pt,[no,np],persistent=.false.,scalable=.false.,supress_error=.true.)
    call mem%allocate(pu,[no,np],persistent=.false.,scalable=.false.,supress_error=.true.)
    pt=0
    pu=transpose(perm)

    ! some temporary space
    ! It can either be done in parallel or serially, so two version are given here.
    if ( present(mw) ) then
        ! So, sort each column.
        do i=1,size(pt,2)
            ! parallel over particles
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            pt(:,i)=pu(:,i)
            call rl_qsort(pt(:,i))
        enddo
        call mw%allreduce('sum',pt)

        ! The unique columns will be the coordination shells, I think.
        ! Should really make a parallel version of this somehow.
        call mem%allocate(shell_index,np,persistent=.true.,scalable=.false.)
        shell_index=0
        call rl_return_unique_indices(pt,di,redind=shell_index)

        ! return which the unique are
        if ( present(prototype_shell) ) then
            call mem%allocate(prototype_shell,size(di),persistent=.true.,scalable=.false.)
            prototype_shell=di
        endif

        ! Store things as shells
        n_shell=size(di)
        call mem%allocate(shell_member,[no,n_shell],persistent=.true.,scalable=.false.)
        call mem%allocate(shell_ctr,n_shell,persistent=.true.,scalable=.false.)
        shell_member=0
        shell_ctr=0
        do i=1,n_shell
            ! parallel over shells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            j=di(i)
            call rl_return_unique(pt(:,j),dk,mem)
            shell_ctr(i)=size(dk)
            shell_member(1:shell_ctr(i),i)=dk
            call mem%deallocate(dk,persistent=.true.,scalable=.false.)
        enddo
        ! Reduce over MPI.
        call mw%allreduce('sum',shell_ctr)
        call mw%allreduce('sum',shell_member)

        ! Small sanity test that I think always must hold
        if ( sum(shell_ctr) .ne. np ) then
            call rl_stop_gracefully(['Lost particle when dividing into shells.'],rl_exitcode_symmetry,mw%comm)
        endif
        ! Then I could think of one more test
        deallocate(di)
        !call mem%deallocate(di,persistent=.true.,scalable=.false.)

        call mem%allocate(di,np,persistent=.false.,scalable=.false.)
        di=-1
        do i=1,n_shell
        do j=1,shell_ctr(i)
            di( shell_member(j,i) )=di( shell_member(j,i) )+1
        enddo
        enddo
        if ( sum(abs(di)) .ne. 0 ) then
            call rl_stop_gracefully(['Cluster does not divide cleanly into shells.'],rl_exitcode_symmetry,mw%comm)
        endif
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    else
        ! Same thing as above, but not parallel.
        do i=1,size(pt,2)
            pt(:,i)=pu(:,i)
            call rl_qsort(pt(:,i))
        enddo

        ! The unique columns will be the coordination shells, I think.
        allocate(shell_index(np))
        shell_index=0
        call rl_return_unique_indices(pt,di,redind=shell_index)

        ! return which the unique are
        if ( present(prototype_shell) ) then
            call mem%allocate(prototype_shell,size(di),persistent=.true.,scalable=.false.)
            prototype_shell=di
        endif

        ! Store things as shells
        n_shell=size(di)
        call mem%allocate(shell_member,[no,n_shell],persistent=.true.,scalable=.false.)
        call mem%allocate(shell_ctr,n_shell,persistent=.true.,scalable=.false.)
        shell_member=0
        shell_ctr=0
        do i=1,n_shell
            j=di(i)
            call rl_return_unique(pt(:,j),dk,mem)
            shell_ctr(i)=size(dk)
            shell_member(1:shell_ctr(i),i)=dk
            call mem%deallocate(dk,persistent=.true.,scalable=.false.)
        enddo

        ! Small sanity test that I think always must hold
        if ( sum(shell_ctr) .ne. np ) then
            call rl_stop_gracefully(['Cluster does not divide cleanly into shells.'],rl_exitcode_symmetry,mw%comm)
        endif
        ! Then I could think of one more test
        !call mem%deallocate(di,persistent=.true.,scalable=.false.)
        deallocate(di)

        call mem%allocate(di,np,persistent=.false.,scalable=.false.)
        di=-1
        do i=1,n_shell
        do j=1,shell_ctr(i)
            di( shell_member(j,i) )=di( shell_member(j,i) )+1
        enddo
        enddo
        if ( sum(abs(di)) .ne. 0 ) then
            call rl_stop_gracefully(['Cluster does not divide cleanly into shells.'],rl_exitcode_symmetry,mw%comm)
        endif
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    endif

!    if ( mw%talk ) then
!    endif

    ! And cleanup
    call mem%deallocate(pt,persistent=.false.,scalable=.false.)
    call mem%deallocate(pu,persistent=.false.,scalable=.false.)

! write(*,*) 'exit symshells:'
! call mem%dump()
end subroutine

! Keep a heap-sort routine here, in case I want to switch to that at some point.
! I have not run into the crossover between heapsort and quicksort, not yet
! at least.
! subroutine heapsort(a)
!   real(r8), intent(inout) :: a(0:)
!   integer :: start, n, bottom
!   real(r8) :: temp
!
!   n = size(a)
!   do start = (n - 2) / 2, 0, -1
!     call siftdown(a, start, n);
!   end do
!
!   do bottom = n - 1, 1, -1
!     temp = a(0)
!     a(0) = a(bottom)
!     a(bottom) = temp;
!     call siftdown(a, 0, bottom)
!   end do
! end subroutine heapsort
!
! subroutine siftdown(a, start, bottom)
!  real(r8), intent(in out) :: a(0:)
!  integer, intent(in) :: start, bottom
!  integer :: child, root
!  real(r8) :: temp
!
!  root = start
!  do while(root*2 + 1 < bottom)
!    child = root * 2 + 1
!
!    if (child + 1 < bottom) then
!      if (a(child) < a(child+1)) child = child + 1
!    end if
!
!    if (a(root) < a(child)) then
!      temp = a(child)
!      a(child) = a (root)
!      a(root) = temp
!      root = child
!    else
!      return
!    end if
!  end do
! end subroutine siftdown

!> Prime factorization of numbers that are not particularly large.
subroutine rl_find_prime_factors(n, d)
    !> number to factorize
    integer, intent(in) :: n
    !> array with factors such that product(d)=n, sorted largest to smallest.
    integer, dimension(:), allocatable, intent(out) :: d

    integer :: div, next, rest
    integer :: i,ctr

    ctr=0
    i = 1
    div=2
    next=3
    rest=n
    do while ( rest .ne. 1 )
        do while ( mod(rest, div) == 0 )
            ctr=ctr+1
            i=i+1
            rest=rest/div
        enddo
        div=next
        next=next+2
    end do

    allocate(d(ctr))
    d=0
    i = 1
    div=2
    next=3
    rest=n
    do while ( rest .ne. 1 )
        do while ( mod(rest, div) == 0 )
            d(i) = div
            i=i+1
            rest=rest/div
        enddo
        div=next
        next=next+2
    end do

    ! And sort it largest to smallest.
    d=-d
    call rl_qsort(d)
    d=-d
end subroutine

end module
