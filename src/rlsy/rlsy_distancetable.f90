module rlsy_distancetable
!!
!! Type to handle calculation of distance tables, with or without
!! periodic boundary conditions for arbitrary sets of points. It is
!! somewhat fast, but not classical-md-fast.
!!
use rlsy_constants, only: r8,rl_iou,rl_huge,rl_hugeint,rl_tol,rl_sqtol,&
                          rl_status,rl_exitcode_param,rl_exitcode_symmetry
use rlsy_helpers, only: tochar,rl_clean_fractional_coordinates,rl_chop,&
                        rl_sqnorm,rl_cross,rl_invert3x3matrix,norm2
use rlsy_sorting, only: rl_qsort
use rlsy_geometry, only: rl_inscribed_sphere_in_box,rl_bounding_sphere_of_box
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime

implicit none
private
public :: rl_distancetable

!> particle in a distance table
type rl_distancetable_particle
    !> how many neighbours are there?
    integer :: n=-rl_hugeint
    !> What is the vector pointing towards these particles?
    real(r8), dimension(:,:), allocatable :: v
    !> Vector that wraps the periodic boundary things
    real(r8), dimension(:,:), allocatable :: lv
    !> What is the index in the original cell of the neighbouring particle?
    integer, dimension(:), allocatable :: ind
    !> What is the distance
    real(r8), dimension(:), allocatable :: d
end type

! A distance table for a set of points. Handles periodic and non-periodic
! boundary conditions, and should be reasonably fast. Not classical-MD-fast,
! but fast enough for most uses. With sufficiently many points and short
! cutoff, it automatically switches to O(N) Verlet lists. Or at least it
! will, once I need it. Right now it just throws an error instead.
type rl_distancetable
    !> Number of particles
    integer :: np=-rl_hugeint
    !> The cutoff
    real(r8) :: cutoff=-rl_huge
    !> The list originating from each particle
    type(rl_distancetable_particle), dimension(:), allocatable :: particle
    contains
        !> create the distance table
        procedure :: generate
        !> largest number of neighbours
        procedure :: max_number_of_neighbours
        !> smallest number of neighbours
        procedure :: min_number_of_neighbours
        !> count neighbours within cutoff
        procedure :: count_neighbours_within_cutoff
        !> remove neighbours outside some cutoff
        procedure :: prune
        !> size in memory
        procedure :: size_in_mem
        !> destroy
        procedure :: destroy
end type

contains

!> generate distance table, periodic boundary conditions.
subroutine generate(dt,particles,basis,cutoff,verbosity,mw,tolerance)
    !> The distance table
    class(rl_distancetable), intent(out) :: dt
    !> The list of particles in fractional coordinates
    real(r8), dimension(:,:), intent(in) :: particles
    !> The basis vectors for the particles
    real(r8), dimension(3,3), intent(in) :: basis
    !> Maximum distance to consider
    real(r8), intent(in) :: cutoff
    !> How much to talk
    integer, intent(in) :: verbosity
    !> MPI helper
    type(rl_mpi_helper), intent(inout), optional :: mw
    !> Tolerance
    real(r8), intent(in), optional :: tolerance

    ! with fever particles than this it's not worth using a verlet list.
    integer, parameter :: verletlist_crossover=800
    ! largest number of Verlet boxes. Should be plenty for anything reasonable.
    integer, parameter :: maxnbox=40
    ! it is bad if the cutoff is exactly on a shell, so buffer it a little.
    real(r8) :: cutoffbuf

    integer :: np,algo
    integer, dimension(3) :: nrep,nbox
    real(r8), dimension(:,:), allocatable :: cleanpositions,rcart
    real(r8) :: timer,rc,sqrc,tol,sqtol
    logical :: verletlist,mpi
    !init: block
    real(r8), dimension(3,3) :: m0
    integer :: i,j
    !repalgo: block
    real(r8), dimension(:,:,:), allocatable :: dumr,dumlv
    real(r8), dimension(:,:), allocatable :: dumd,latvec,llv
    real(r8), dimension(:), allocatable :: dr
    real(r8), dimension(3) :: v0,v1,v2
    real(r8) :: buf
    integer, dimension(:,:), allocatable :: dumind,dumjnd
    integer, dimension(:), allocatable :: dumctr,atind,lai,rnkctr,offset
    !integer :: maxnn,ncells,ctr,ii,jj,kk,a1,a2,l,i,j
    integer :: maxnn,ncells,ctr,ii,jj,kk,a1,a2,l
    !repalgoser: block
    !real(r8), dimension(:,:,:), allocatable :: dumr
    !real(r8), dimension(:,:,:), allocatable :: dumlv
    !real(r8), dimension(:,:), allocatable :: dumd,latvec
    !real(r8), dimension(:), allocatable :: dr
    !real(r8), dimension(3) :: v0,v1,v2
    !real(r8) :: buf
    !integer, dimension(:,:), allocatable :: dumind,dumjnd
    !integer, dimension(:), allocatable :: dumctr,atind
    !integer :: maxnn,ncells,a1,a2,i,j,ii,jj,kk,l
    !dirnn: block
    !real(r8), dimension(:,:,:), allocatable :: dumr
    !real(r8), dimension(:,:,:), allocatable :: dumlv
    !real(r8), dimension(:,:), allocatable :: dumd
    !real(r8), dimension(3) :: v0,v1,v2
    !real(r8) :: buf
    !integer, dimension(:,:), allocatable :: dumind
    !integer, dimension(:), allocatable :: sortind,dumctr,offset
    integer, dimension(:), allocatable :: sortind
    !integer :: maxnn,a1,a2,i,j,l
    !dirnnser: block
    integer, parameter :: countcrossover=100
    !real(r8), dimension(:,:,:), allocatable :: dumr,dumlv
    !real(r8), dimension(:,:), allocatable :: dumd
    !real(r8), dimension(3) :: v0,v1,v2
    !real(r8) :: buf
    !integer, dimension(:,:), allocatable :: dumind
    !integer, dimension(:), allocatable :: sortind,dumctr
    !integer :: a1,a2,l,i,j,maxnn

    ! set some basic stuff
    !init: block

        ! Start timer if talking
        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
        endif

        ! Tolerance?
        if ( present(tolerance) ) then
            tol=tolerance
            sqtol=tol**2
        else
            tol=rl_tol
            sqtol=rl_sqtol
        endif

        ! Should it be done in parallel?
        if ( present(mw) ) then
            mpi=.true.
        else
            mpi=.false.
        endif

        ! number of particles
        np=size(particles,2)
        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... building distance table for '//tochar(np)//' points'
        endif
        ! Add buffer to cutoff, will explain later why I use a buffer.
        cutoffbuf=100*tol
        rc=cutoff+cutoffbuf
        sqrc=(cutoff+cutoffbuf)**2

        ! Clean the positions to make sure they are [0-1[
        ! This does not affect the resulting distance table in any
        ! way, only makes evaluation faster.
        allocate(cleanpositions(3,np))
        cleanpositions=rl_clean_fractional_coordinates(particles)
        if ( verbosity .gt. 0 ) write(rl_iou,*) '... cleaned the positions'

        ! The cutoff is likely larger than the periodic box the particles
        ! live in, so we have to know how many times I have to repeat the
        ! box to catch all neighbours. I try to do this in a sensibly clever
        ! way to ensure that I repeat the box the minimum number of times
        ! required, that could be different number of repetitions in each
        ! direction.
        nrep=0
        if ( rl_inscribed_sphere_in_box(basis) .lt. rc ) then
            do j=1,100
                do i=1,3
                    m0(:,i)=basis(:,i)*(2*nrep(i)+1)
                enddo
                if ( rl_inscribed_sphere_in_box(m0)-rl_bounding_sphere_of_box(basis) .gt. cutoff+cutoffbuf ) exit
                nrep=increment_dimensions(nrep,basis)
            enddo
        endif

        ! Now, if no repetitions are needed, it might be worth trying a Verlet-list type calculation of distances.
        ! For that to be any kind of speedup, I must be able to divide the box into at least 4 sub-boxes, where
        ! each sub-box must be large enough to contain the cutoff. Also, the number of particles needs to be large
        ! enough to validate the extra overhead of switching to an O(NlogN) algorithm.
        !@TODO Actually implement this. Has never been needed so far, but when the first error message comes I know what to do.
        if ( np .gt. verletlist_crossover ) then
            nbox=0
            do i=1,3
                do j=maxnbox,1,-1
                    m0=basis
                    m0(:,i)=m0(:,i)/j
                    if ( rl_inscribed_sphere_in_box(m0) .gt. rc ) then
                        nbox(i)=j
                        exit
                    endif
                enddo
            enddo
            ! Now I can check if the divisions are large enough to use.
            if ( minval(nbox) .gt. 3 ) then
                verletlist=.true.
            else
                verletlist=.false.
            endif
        else
            ! no point in even trying with this few particles.
            verletlist=.false.
        endif

        ! Ok, now I have enough information to determine wich algorithm to use.
        algo=0
        if ( sum(nrep) .eq. 0 ) then
            if ( mpi ) then
                if ( verletlist ) then
                    algo=0
                else
                    algo=3
                endif
            else
                if ( verletlist ) then
                    algo=0
                else
                    algo=4
                endif
            endif
        else
            if ( mpi ) then
                algo=1
            else
                algo=2
            endif
        endif

        if ( algo .eq. 0 ) then
            call rl_stop_gracefully(['Nag on Olle to implement Verlet lists.'],&
                                     rl_exitcode_param)
        endif

        ! For some algorithms, I pre-convert to Cartesian coordinates
        if ( algo .le. 2 ) then
            allocate(rcart(3,np))
            call dgemm('N','N',3,np,3,1.0_r8,basis,3,cleanpositions,3,0.0_r8,rcart,3)
        endif

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... using',tochar(nrep),' repetitions to build the distance table'
        endif
    !end block init

    ! Now for the individual algorithms
    select case(algo)
    case(1)
    !repalgo: block
        ! Use filtered repetitions, parallel over MPI

        ! Parallel things
        allocate(rnkctr(mw%n))
        allocate(offset(mw%n))
        rnkctr=0
        offset=0

        ! First step is to generate a list of latticevectors and atoms, local
        ! to each rank
        ncells=0
        ctr=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                v1=matmul(basis,[ii,jj,kk]*1.0_r8)
                ! First the fast test
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! Then angrier test. Disabled for now, was not stable.
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
            enddo
            enddo
            enddo
        enddo

        allocate(llv(3,ncells))
        allocate(lai(ncells))
        llv=0.0_r8
        lai=0
        ncells=0
        ctr=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                v1=matmul(basis,[ii,jj,kk]*1.0_r8)
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
                lai(ncells)=a1
                llv(:,ncells)=v1
            enddo
            enddo
            enddo
        enddo

        ! Figure out how many cells I have on each rank, and the offset
        rnkctr=0
        offset=0
        rnkctr(mw%r+1)=ncells
        call mw%allreduce('sum',rnkctr)
        do i=1,mw%n-1
            offset(i+1)=sum(rnkctr(1:i))
        enddo

        ! Make large arrays to store all cells on all ranks?
        allocate(atind(sum(rnkctr)))
        allocate(latvec(3,sum(rnkctr)))
        atind=0
        latvec=0.0_r8
        ! Store the latticevectors and atom indices in the right place
        do i=1,ncells
            j=offset(mw%r+1)+i
            atind(j)=lai(i)
            latvec(:,j)=llv(:,i)
        enddo
        deallocate(lai)
        deallocate(llv)
        call mw%allreduce('sum',atind)
        call mw%allreduce('sum',latvec)
        ! And make sure total number of cells is known
        ncells=sum(rnkctr)

        ! That was the first level of parallelism. Now count the number of neighbours? Not sure what
        ! the best way is, perhaps a big flat array thing?
        allocate(dumctr(np))
        dumctr=0
        do i=1,ncells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( rl_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
            enddo
        enddo

        ! This would be the size of my buffer thingies.
        ! Pretty sure it's ok, a little waste of memory
        ! but should be fine.
        maxnn=maxval(dumctr)
        call mw%allreduce('max',maxnn)
        do i=1,mw%n
            offset(i)=(i-1)*maxnn
        enddo
        maxnn=maxnn*mw%n
        allocate(dumr(3,maxnn,np))
        allocate(dumlv(3,maxnn,np))
        allocate(dumd(maxnn,np))
        allocate(dumind(maxnn,np))
        allocate(dumjnd(maxnn,np))
        allocate(dr(maxnn))
        dumr=0.0_r8
        dumlv=0.0_r8
        dumd=1E10_r8*rc
        dumind=0
        dumjnd=0
        dumctr=0

        ! Now store vectors and stuff
        do i=1,ncells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( rl_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
                j=dumctr(a1)+offset(mw%r+1)
                dumr(:,j,a1)=v2
                dumlv(:,j,a1)=v1
                dumd(j,a1)=norm2( v2 )
                dumind(j,a1)=a2
            enddo
        enddo
        deallocate(atind)
        deallocate(latvec)
        ! Communicate everywhere
        call mw%allreduce('sum',dumr)
        call mw%allreduce('sum',dumlv)
        call mw%allreduce('min',dumd)
        call mw%allreduce('sum',dumind)
        call mw%allreduce('sum',dumctr)

        ! Sort things by distance.
        dumjnd=0
        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            dr=dumd(:,a1)
            call rl_qsort(dr,dumjnd(:,a1))
        enddo
        call mw%allreduce('sum',dumjnd)

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_r8
        do j=1,1000
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(dumjnd(i,a1),a1) ) .lt. 2*tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+tol
            endif
        enddo
        if ( verbosity .gt. 0 .and. buf .gt. sqtol ) then
            write(rl_iou,*) '... increased cutoff a little to accomodate full shells'
        endif

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_r8
        allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            allocate(dt%particle(a1)%d(l))
            allocate(dt%particle(a1)%v(3,l))
            allocate(dt%particle(a1)%lv(3,l))
            allocate(dt%particle(a1)%ind(l))
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) then
                    l=l+1
                    dt%particle(a1)%d(l)=rl_chop( dumd(j,a1), sqtol )
                    dt%particle(a1)%v(:,l)=rl_chop( dumr(:,j,a1), sqtol )
                    dt%particle(a1)%lv(:,l)=rl_chop( dumlv(:,j,a1), sqtol )
                    dt%particle(a1)%ind(l)=dumind(j,a1)
                endif
            enddo
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        deallocate(dumd)
        deallocate(dumr)
        deallocate(dumlv)
        deallocate(dumind)
        deallocate(dumjnd)
        deallocate(dumctr)
        deallocate(rnkctr)
        deallocate(offset)
        deallocate(dr)
    !end block repalgo
    case(2)
    !repalgoser: block
        ! Serial version of filtered repetitions

        ! First step is to generate a list of latticevectors and atoms, local
        ! to each rank
        ncells=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                v1=matmul(basis,[ii,jj,kk]*1.0_r8)
                ! First the fast test
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! Then angrier test
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
            enddo
            enddo
            enddo
        enddo

        allocate(latvec(3,ncells))
        allocate(atind(ncells))
        latvec=0.0_r8
        atind=0
        ncells=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                v1=matmul(basis,[ii,jj,kk]*1.0_r8)
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
                atind(ncells)=a1
                latvec(:,ncells)=v1
            enddo
            enddo
            enddo
        enddo

        ! That was the first level of parallelism. Now count the number of neighbours? Not sure what
        ! the best way is, perhaps a big flat array thing?
        allocate(dumctr(np))
        dumctr=0
        do i=1,ncells
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( rl_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
            enddo
        enddo

        ! This would be the size of my buffer thingies.
        ! Pretty sure it's ok, a little waste of memory
        ! but should be fine.
        maxnn=maxval(dumctr)
        allocate(dumr(3,maxnn,np))
        allocate(dumlv(3,maxnn,np))
        allocate(dumd(maxnn,np))
        allocate(dumind(maxnn,np))
        allocate(dumjnd(maxnn,np))
        allocate(dr(maxnn))
        dumr=0.0_r8
        dumlv=0.0_r8
        dumd=1E10_r8*rc
        dumind=0
        dumjnd=0
        dumctr=0

        ! Now store vectors and stuff
        do i=1,ncells
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( rl_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
                j=dumctr(a1)
                dumr(:,j,a1)=v2
                dumlv(:,j,a1)=v1
                dumd(j,a1)=norm2( v2 )
                dumind(j,a1)=a2
            enddo
        enddo
        ! Sort things by distance
        dumjnd=0
        do a1=1,np
            dr=dumd(:,a1)
            call rl_qsort(dr,dumjnd(:,a1))
        enddo

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_r8
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(dumjnd(i,a1),a1) ) .lt. 2*tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+rl_tol
            endif
        enddo
        if ( verbosity .gt. 0 .and. buf .gt. sqtol ) then
            write(rl_iou,*) '... increased cutoff a little to accomodate full shells'
        endif

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_r8
        allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            allocate(dt%particle(a1)%d(l))
            allocate(dt%particle(a1)%v(3,l))
            allocate(dt%particle(a1)%lv(3,l))
            allocate(dt%particle(a1)%ind(l))
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) then
                    l=l+1
                    dt%particle(a1)%d(l)=rl_chop( dumd(j,a1), sqtol )
                    dt%particle(a1)%v(:,l)=rl_chop( dumr(:,j,a1), sqtol )
                    dt%particle(a1)%lv(:,l)=rl_chop( dumlv(:,j,a1), sqtol )
                    dt%particle(a1)%ind(l)=dumind(j,a1)
                endif
            enddo
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        deallocate(dumd)
        deallocate(dumr)
        deallocate(dumlv)
        deallocate(dumind)
        deallocate(dumjnd)
        deallocate(dumctr)
        deallocate(dr)
        deallocate(latvec)
        deallocate(atind)
    !end block repalgoser
    case(3)
    !dirnn: block
        ! Direct calculation, in parallel

        ! Count to get an upper bound on the number of neighbours
        allocate(dumctr(np))
        allocate(offset(mw%n))
        dumctr=0
        offset=0
        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( rl_sqnorm(v1) .gt. sqrc ) cycle
                l=l+1
            enddo
            dumctr(a1)=l
        enddo

        ! Make some space for buffers
        maxnn=maxval(dumctr)
        call mw%allreduce('max',maxnn)
        do i=1,mw%n
            offset(i)=(i-1)*maxnn
        enddo
        maxnn=maxnn*mw%n
        allocate(dumr(3,maxnn,np))
        allocate(dumlv(3,maxnn,np))
        allocate(dumd(maxnn,np))
        allocate(dumind(maxnn,np))
        allocate(sortind(maxnn))
        dumr=0.0_r8
        dumlv=0.0_r8
        dumd=1E10_r8*rc
        dumind=0
        dumctr=0
        sortind=0

        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( rl_sqnorm(v1) .gt. sqrc ) cycle
                ! lattice vector
                v2=v1-v0
                v2=matmul(basis,anint(v2))
                l=l+1
                dumr(:,l,a1)=v1
                dumlv(:,l,a1)=v2
                dumd(l,a1)=norm2( v1 )
                dumind(l,a1)=a2
            enddo
            dumctr(a1)=l
            ! Sort by distance
            call rl_qsort(dumd(1:l,a1),sortind(1:l))
            dumr(:,1:l,a1)=dumr(:,sortind(1:l),a1)
            dumlv(:,1:l,a1)=dumlv(:,sortind(1:l),a1)
            dumind(1:l,a1)=dumind(sortind(1:l),a1)
        enddo
        call mw%allreduce('sum',dumctr)
        call mw%allreduce('sum',dumr)
        call mw%allreduce('min',dumd)
        call mw%allreduce('sum',dumlv)
        call mw%allreduce('sum',dumind)

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_r8
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(i,a1) ) .lt. 2*rl_tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+rl_tol
            endif
        enddo

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_r8
        allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                if ( dumd(i,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            allocate(dt%particle(a1)%d(l))
            allocate(dt%particle(a1)%v(3,l))
            allocate(dt%particle(a1)%lv(3,l))
            allocate(dt%particle(a1)%ind(l))
            dt%particle(a1)%d=rl_chop( dumd(1:l,a1),sqtol )
            dt%particle(a1)%v=rl_chop( dumr(:,1:l,a1),sqtol )
            dt%particle(a1)%lv=rl_chop( dumlv(:,1:l,a1),sqtol )
            dt%particle(a1)%ind=dumind(1:l,a1)
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        deallocate(dumd)
        deallocate(dumr)
        deallocate(dumlv)
        deallocate(dumind)
        deallocate(dumctr)
        deallocate(sortind)
    !end block dirnn
    case(4)
    !dirnnser: block
        ! Direct calculation, no repetition, serial

        ! Upper bound for number of neighbours
        if ( np .gt. countcrossover ) then
            maxnn=0
            do a1=1,np
                l=0
                do a2=1,np
                    v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                    v1=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
                    v1=matmul(basis,v1)
                    if ( rl_sqnorm(v1) .gt. sqrc ) cycle
                    l=l+1
                enddo
                maxnn=max(maxnn,l)
            enddo
        else
            maxnn=np
        endif

        allocate(dumctr(np))
        allocate(dumr(3,maxnn,np))
        allocate(dumlv(3,maxnn,np))
        allocate(dumd(maxnn,np))
        allocate(dumind(maxnn,np))
        allocate(sortind(maxnn))
        dumr=0.0_r8
        dumlv=0.0_r8
        dumd=1E10_r8*rc
        dumind=0
        dumctr=0
        sortind=0
        do a1=1,np
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( rl_sqnorm(v1) .gt. sqrc ) cycle
                ! lattice vector
                v2=v1-v0
                v2=matmul(basis,anint(v2))
                l=l+1
                dumr(:,l,a1)=v1
                dumlv(:,l,a1)=v2
                dumd(l,a1)=norm2(v1)
                dumind(l,a1)=a2
            enddo
            dumctr(a1)=l
            ! Sort by distance
            call rl_qsort(dumd(1:l,a1),sortind(1:l))
            dumr(:,1:l,a1)=dumr(:,sortind(1:l),a1)
            dumlv(:,1:l,a1)=dumlv(:,sortind(1:l),a1)
            dumind(1:l,a1)=dumind(sortind(1:l),a1)
        enddo

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_r8
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(i,a1) ) .lt. 2*rl_tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+rl_tol
            endif
        enddo
        ! now store everything
        dt%np=np
        dt%cutoff=0.0_r8
        allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                if ( dumd(i,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            allocate(dt%particle(a1)%d(l))
            allocate(dt%particle(a1)%v(3,l))
            allocate(dt%particle(a1)%lv(3,l))
            allocate(dt%particle(a1)%ind(l))
            dt%particle(a1)%d=rl_chop( dumd(1:l,a1), rl_sqtol )
            dt%particle(a1)%v=rl_chop( dumr(:,1:l,a1), rl_sqtol )
            dt%particle(a1)%lv=rl_chop( dumlv(:,1:l,a1), rl_sqtol )
            dt%particle(a1)%ind=dumind(1:l,a1)
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        deallocate(dumd)
        deallocate(dumr)
        deallocate(dumlv)
        deallocate(dumind)
        deallocate(dumctr)
        deallocate(sortind)
    !end block dirnnser
    end select

    ! And some final cleanup
    if ( allocated(cleanpositions) ) deallocate(cleanpositions)
    if ( allocated(rcart) ) deallocate(rcart)
    if ( verbosity .gt. 0 ) write(rl_iou,*) '... done generating distance table (',tochar(mpi_wtime()-timer),'s)'
end subroutine

pure function shortest_distance_cell_to_point(latticevector,point,basis) result(distance)
    real(r8), dimension(3), intent(in) :: latticevector
    real(r8), dimension(3), intent(in) :: point
    real(r8), dimension(3,3), intent(in) :: basis
    real(r8) :: distance

    real(r8), dimension(3) :: v0
    real(r8) :: r0

    ! Bounding sphere of box
    r0=rl_bounding_sphere_of_box(basis)
    ! Get the center of the box
    v0=latticevector+matmul(basis,[0.5_r8,0.5_r8,0.5_r8])-point
    distance=norm2(v0)-r0
end function

!> Get the largest number of neighours of any particle
pure function max_number_of_neighbours(dt) result(n)
    !> distance table
    class(rl_distancetable), intent(in) :: dt
    !> number
    integer :: n

    integer :: i

    n=-1
    do i=1,dt%np
        n=max(dt%particle(i)%n,n)
    enddo
end function

!> Get the smallest number of neighours of any particle
pure function min_number_of_neighbours(dt) result(n)
    !> distance table
    class(rl_distancetable), intent(in) :: dt
    !> number
    integer :: n

    integer :: i

    n=huge(n)
    do i=1,dt%np
        n=min(dt%particle(i)%n,n)
    enddo
end function

!> remove unwanted neighbours from the distance table
subroutine prune(dt,cutoff_per_particle)
    !> distance table
    class(rl_distancetable), intent(inout) :: dt
    !> max cutoff, per particle
    real(r8), dimension(dt%np), intent(in) :: cutoff_per_particle

    integer :: i,j,n
    real(r8), dimension(:,:), allocatable :: v
    real(r8), dimension(:,:), allocatable :: lv
    real(r8), dimension(:), allocatable :: d
    real(r8), dimension(:), allocatable :: weight
    integer, dimension(:), allocatable :: ind
    !
    i=dt%max_number_of_neighbours()
    allocate(v(3,i))
    allocate(lv(3,i))
    allocate(d(i))
    allocate(weight(i))
    allocate(ind(i))

    do i=1,dt%np
        ! see how many I will keep. I assume that the distance table is sorted properly.
        n=-1
        do j=1,dt%particle(i)%n
            if ( dt%particle(i)%d(j)-rl_tol .gt. cutoff_per_particle(i) ) then
                n=j
                exit
            endif
        enddo
        ! maybe no pruning needed?
        if ( n .eq. dt%particle(i)%n ) cycle
        ! make a copy and reallocate
        if( allocated(dt%particle(i)%v) ) then
            v(:,1:n)=dt%particle(i)%v(:,1:n)
            deallocate(dt%particle(i)%v)
            !allocate(dt%particle(i)%v,source=v(:,1:n))
            allocate(dt%particle(i)%v(3,1:n))
            dt%particle(i)%v=v(:,1:n)
        endif
        if( allocated(dt%particle(i)%lv) ) then
            lv(:,1:n)=dt%particle(i)%lv(:,1:n)
            deallocate(dt%particle(i)%lv)
            !allocate(dt%particle(i)%lv,source=lv(:,1:n))
            allocate(dt%particle(i)%lv(3,1:n))
            dt%particle(i)%lv=lv(:,1:n)
        endif
        if( allocated(dt%particle(i)%ind) ) then
            ind(1:n)=dt%particle(i)%ind(1:n)
            deallocate(dt%particle(i)%ind)
            !allocate(dt%particle(i)%ind,source=ind(1:n))
            allocate(dt%particle(i)%ind(1:n))
            dt%particle(i)%ind=ind(1:n)
        endif
        ! if( allocated(dt%particle(i)%weight) ) then
        !     weight(1:n)=dt%particle(i)%weight(1:n)
        !     deallocate(dt%particle(i)%weight)
        !     allocate(dt%particle(i)%weight,source=weight(1:n))
        ! endif
        ! and store the new number of neighbours
        dt%particle(i)%n=n
    enddo
end subroutine

!> increase supercell dimensions in reasonably clever way
function increment_dimensions(dimin,box) result(dimut)
    integer, dimension(3), intent(in) :: dimin
    real(r8), dimension(3,3), intent(in) :: box
    integer, dimension(3) :: dimut
    !
    real(r8), dimension(3,3) :: m0
    integer, dimension(3) :: di
    integer :: i,j,k
    real(r8) :: f0,f1,ff0

    ! try with the sphere thing. First get a baseline
    do j=1,3
        m0(:,j)=box(:,j)*(2*dimin(j)+1)
    enddo
    ff0=rl_inscribed_sphere_in_box(m0)
    ! Increment the dimension that gives the biggest increase for
    ! the radius of the inscribed sphere.
    f0=0.0_r8
    dimut=0
    do i=1,3
        di=dimin
        di(i)=di(i)+1
        do j=1,3
            m0(:,j)=box(:,j)*(2*di(j)+1)
        enddo
        f1=rl_inscribed_sphere_in_box(m0)
        if ( f1 .gt. f0 .and. abs(f1-ff0) .gt. rl_tol ) then
            dimut=di
            f0=f1
        endif
    enddo

    ! if nothing helped, increment the lowest number.
    ! If all dimensions equal, just pick one to increase.
    if ( dimut(1) .eq. 0 ) then
        j=rl_hugeint
        k=0
        do i=1,3
            if ( di(i) .lt. j ) then
                j=di(i)
                k=i
            endif
        enddo
        dimut=dimin
        dimut(k)=dimut(k)+1
    endif
end function

!> Count the number of neighbours per atom within a certain cutoff. Not including itself.
subroutine count_neighbours_within_cutoff(dt,cutoff,ctr)
    !> distance table
    class(rl_distancetable), intent(in) :: dt
    !> cutoff
    real(r8), intent(in) :: cutoff
    !> how many neighbours for each atom?
    integer, dimension(:), intent(out) :: ctr

    integer :: i,j,k

    ctr=0
    ! it might be really simple
    if ( cutoff .gt. dt%cutoff-rl_tol ) then
        do i=1,dt%np
            ctr(i)=dt%particle(i)%n-1
        enddo
        return
    endif
    do i=1,dt%np
        k=0
        do j=2,dt%particle(i)%n
            if ( dt%particle(i)%d(j) .lt. cutoff+rl_tol ) then
                k=k+1
            else
                exit
            endif
        enddo
        ctr(i)=k
    enddo
end subroutine

!> measure size in memory, in bytes
function size_in_mem(dt) result(mem)
    !> distance table
    class(rl_distancetable), intent(in) :: dt
    !> memory
    integer :: mem

    integer :: i

    mem=0
    mem=mem+storage_size(dt)
    if ( allocated(dt%particle) ) then
        do i=1,size(dt%particle)
            mem=mem+storage_size(dt%particle(i))
            if ( allocated(dt%particle(i)%v  ) ) mem=mem+storage_size(dt%particle(i)%v  )*size(dt%particle(i)%v  )
            if ( allocated(dt%particle(i)%lv ) ) mem=mem+storage_size(dt%particle(i)%lv )*size(dt%particle(i)%lv )
            if ( allocated(dt%particle(i)%ind) ) mem=mem+storage_size(dt%particle(i)%ind)*size(dt%particle(i)%ind)
            if ( allocated(dt%particle(i)%d  ) ) mem=mem+storage_size(dt%particle(i)%d  )*size(dt%particle(i)%d  )
        enddo
    endif
    mem=mem/8
end function

!> deallocate everything
subroutine destroy(dt)
    !> distance table
    class(rl_distancetable), intent(inout) :: dt

    integer :: i

    if ( allocated(dt%particle) ) then
        do i=1,size(dt%particle)
            if ( allocated(dt%particle(i)%v  ) ) deallocate(dt%particle(i)%v  )
            if ( allocated(dt%particle(i)%lv ) ) deallocate(dt%particle(i)%lv )
            if ( allocated(dt%particle(i)%ind) ) deallocate(dt%particle(i)%ind)
            if ( allocated(dt%particle(i)%d  ) ) deallocate(dt%particle(i)%d  )
        enddo
        deallocate(dt%particle)
    endif
end subroutine

end module
