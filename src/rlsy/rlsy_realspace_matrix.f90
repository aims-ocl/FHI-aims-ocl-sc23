module rlsy_realspace_matrix
!!
!! This module handles the basic storage of a realspace matrix that is kinda large, like
!! the overlap, hamiltonian or density matrix. Makes it distributed, organized via pairs.
!! I think. Have to make up my mind how to do this in a neat way. Will probably iterate
!! many times before I come to any sensible conclusions.
!!
!! I think three schemes of distribution make sense:
!!
!! 1) Everyone has a copy. Fastest.
!!
!! 2) Each group of ranks that handle one k-point has a copy. Avoids all-to-all, communication quite
!!    simple, I think. Maybe useless intermediate thing, not sure.
!!
!! 3) Completely distributed. Slightly annoying communication when Fourier-transforming, but that's
!!    not too bad, I think.
!!
!! They are not implemented all of them, but preliminary tests seem quite promising. After the
!! initial setup, which you don't really have to understand, it is rather easy to work with.
!! There are substantially fewer if-statements than there normally is in AIMS, with pretty much
!! the same functionality.
!!
use rlsy_constants, only: r8,i8,rl_huge,rl_hugeint,rl_tol,rl_sqtol,rl_twopi,&
    rl_exitcode_symmetry,rl_exitcode_param,rl_exitcode_memory,rl_iou
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: tochar,rl_sqnorm,rl_chop,tochar,rl_mom_real,norm2
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_distancetable, only: rl_distancetable
use rlsy_extended_cluster, only: rl_extended_cluster
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_verletlist, only: rl_verletbox
use rlsy_spacegroup, only: rl_spacegroup,rl_spacegroup_operation
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list,rl_distribute_weighted
use rlsy_kpointmesh, only: rl_kpoint_mesh
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc

implicit none
private
public :: rl_realspace_matrix
public :: rl_realspace_matrix_notdistributed
public :: rl_realspace_matrix_mediumdistributed
public :: rl_realspace_matrix_fulldistributed
public :: rl_setup_realspace_storage

!> minimal information about a pair.
type rl_realspace_matrix_pair
    !> in the unit cell, what is the index of the first atom
    integer :: a1=-rl_hugeint
    !> in the unit cell, what is the index of the second atom
    integer :: a2=-rl_hugeint
    !> vector associated with the pair
    real(r8), dimension(3) :: v=-rl_huge
    !> what is the lattice vector associated with the pair?
    real(r8), dimension(3) :: lv=-rl_huge
    !> number of basis functions on atom a1
    integer :: nb1=-rl_hugeint
    !> number of basis functions on atom a2
    integer :: nb2=-rl_hugeint
    contains
        !> rotate a block associated with this pair
        procedure :: rotate=>rl_rotate_realspace_pair_block
end type

!> general pair that can be constructed from an irreducible pair.
type, extends(rl_realspace_matrix_pair) :: rl_realspace_matrix_full_pair
    !> what is the index of the transpose of this pair?
    integer :: index_transpose=-rl_hugeint
    !> which unique pair does this pair correspond to
    integer :: index_unique=-rl_hugeint
    !> which operation takes the unique pair here. Negative index means it should be transposed as well.
    integer :: index_operation=-rl_hugeint
    !> index in the AIMS array, so that I can fetch and set values. Should be retired once not needed anymore.
    integer, dimension(:,:), allocatable :: idx

    ! below pertains to the indexing in the Fourier transform thingy

    !> how many elements are locally relevant for this pair
    integer :: n_ft_element=-rl_hugeint
    !> the relevant local elements
    integer, dimension(:,:), allocatable :: ft_element
end type

!> irreducible pair, all the pairs can be built by operating on this guy
type, extends(rl_realspace_matrix_pair) :: rl_realspace_matrix_irr_pair
    !> how many full pairs does this pair unfold to
    integer :: n_unfold_pair=-rl_hugeint
    !> which pairs (indices) does it unfold to
    integer, dimension(:), allocatable :: unfold_index
    !> which operation does it unfold with, negative index means an additional transpose
    integer, dimension(:), allocatable :: unfold_operation
    !> index in the AIMS array, so that I can fetch and set values. Should be retired once not needed anymore.
    integer, dimension(:,:), allocatable :: idx
    !> Hamiltonian
    real(r8), dimension(:,:,:), allocatable :: hamiltonian
    !> Density matrix
    real(r8), dimension(:,:,:), allocatable :: densitymatrix
    !> Overlap matrix
    real(r8), dimension(:,:), allocatable :: overlap
end type

!> contains a realspace matrix. It may or may not be distributed. Later problem.
type :: rl_realspace_matrix
    !> total number of pairs on this rank, may or may not be the same as the global number.
    integer :: n_full_pair=-rl_hugeint
    !> number of irreducible pairs
    integer :: n_irr_pair=-rl_hugeint
    !> all pairs
    type(rl_realspace_matrix_full_pair), dimension(:), allocatable :: full_pair
    !> irreducible pairs
    type(rl_realspace_matrix_irr_pair), dimension(:), allocatable :: irr_pair
    contains
        !> build the index map from the AIMS hamiltonian to mine
        procedure :: build_index
        !> measure size in memory, approximately
        procedure :: size_in_mem=>mtx_size_in_mem
        !> destroy
        procedure :: destroy=>destroy_mtx

        !> collect data from AIMS
        procedure :: collect_overlap
        procedure :: collect_hamiltonian
        procedure :: collect_densitymatrix
        !> inject data into AIMS
        procedure :: inject_hamiltonian
        procedure :: inject_densitymatrix
end type

!> Distribution scheme when everyone has a copy of everything. Will likely hold something eventually.
type, extends(rl_realspace_matrix) :: rl_realspace_matrix_notdistributed
end type

!> Distribution scheme when each group has a copy of everything.
type, extends(rl_realspace_matrix) :: rl_realspace_matrix_mediumdistributed
    !> local number of irreducible pairs
    integer :: n_irr_pair_global=-rl_hugeint
    !> local number of full pairs
    integer :: n_full_pair_global=-rl_hugeint
    !> offset for irreducible pair indices
    integer :: irr_offset=-rl_hugeint
    !> offset for full pair indices
    integer :: full_offset=-rl_hugeint
end type

!> Distribution scheme that is completely distributed.
type, extends(rl_realspace_matrix) :: rl_realspace_matrix_fulldistributed
    !> local number of irreducible pairs
    integer :: n_irr_pair_global=-rl_hugeint
    !> local number of full pairs
    integer :: n_full_pair_global=-rl_hugeint
    !> offset for irreducible pair indices
    integer :: irr_offset=-rl_hugeint
    !> offset for full pair indices
    integer :: full_offset=-rl_hugeint
end type

! Helper type to sort out the symmetries
type pair_sym_helper
    integer :: n
    real(r8), dimension(:,:), allocatable :: v,w
    integer, dimension(:), allocatable :: ind,jnd,a2
end type

! Magical parameter, decides on the crossover between distributed and not. Should be made smarter.
real(r8), parameter :: storeall_safely=400.0_r8 ! Or some other number, how would I know.

contains

!> create space to hold the distributed matrix. Also decide on distribution pattern.
subroutine rl_setup_realspace_storage(rmtx,p,basis,ec,sym,KS,mw,verbosity,mem)
    !> distributed matrix
    class(rl_realspace_matrix), allocatable, intent(out) :: rmtx
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> eigenproblem handle
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    real(r8) :: timer,t0,t1
    real(r8), dimension(:,:), allocatable :: dr
    integer, dimension(:,:), allocatable :: sizeofpair,di
    integer, dimension(:), allocatable :: tri

    !init: block
        type(rl_distancetable) :: dt
        real(r8), dimension(:), allocatable :: atomcutoff
        real(r8) :: cutoff
        integer, dimension(:), allocatable :: offset,rnkctr
        integer :: i,j,l,ctr,a1,s1,s2

    !symreduceident: block
        type(pair_sym_helper), dimension(:), allocatable :: sh
        type(rl_verletbox), dimension(:), allocatable :: vb
        real(r8), dimension(3) :: v0
        real(r8) :: f0,f1
        integer, dimension(:,:), allocatable :: dk
        !integer :: npair,i,j,k,s1,s2
        integer :: npair,k

    !symreducenormal: block
        !type(pair_sym_helper), dimension(:), allocatable :: sh
        !type(rl_verletbox), dimension(:), allocatable :: vb
        !real(r8), dimension(3) :: v0,v1
        real(r8), dimension(3) :: v1
        !real(r8) :: f0,f1
        !integer, dimension(:,:), allocatable :: permutation,shell_member,dk
        integer, dimension(:,:), allocatable :: permutation,shell_member
        integer, dimension(:), allocatable :: shell_ctr,shell_index,prot_index,dj,ind
        !integer :: i,j,k,l,npair,iop,ia,ja,s1,s2,ii,p1,p2
        integer :: iop,ia,ja,ii,p1,p2

        !integer :: ctr0,ctr1,ctr2

    !setstorage: block
        !integer :: ipair,a1,a2,nb1,nb2,s1,s2
        integer :: ipair,a2,nb1,nb2

    !grablatticevecs: block
        !real(r8), dimension(3) :: v0
        !integer :: i

    !sanitychecksym: block
        type(rl_verletbox) :: vbb
        !real(r8), dimension(3) :: v0,v1
        !real(r8) :: f0
        !integer :: a1,a2,ie,s1,s2,p1,p2
        integer :: ie
        !integer :: i,j,k,l

    !indexingsingleproc: block
        !integer :: ipair

    !indexingmultiproc: block
        !integer :: ipair,a1,a2,s1,s2,nb1,nb2,b1,b2,ob1,ob2
        integer :: b1,b2,ob1,ob2
        !integer :: gi,gj,li,lj,ctr !,bigctr
        integer :: gi,gj,li,lj !,bigctr

    ! Set up some general things.
    !@TODO Add explanation what I do.
    !init: block

        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
            t0=timer
            t1=t0
            write(rl_iou,*) ''
            write(rl_iou,*) 'SETTING UP REALSPACE STORAGE'
        endif

        ! Build the distance table, with a cutoff that should be a safe upper bound for sure.
        cutoff=basis%longest_cutoff*2+10*rl_tol
        call dt%generate(p%fractional_coordinate,p%latticevectors,cutoff,verbosity=-1,mw=mw)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got distance table (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Then I need the cutoff per atom
        call mem%allocate(atomcutoff,p%n_atom,scalable=.false.,persistent=.false.)
        atomcutoff=0
        do a1=1,p%n_atom
            s1=p%species(a1)
            atomcutoff(a1)=maxval(basis%species(s1)%basis_cutoff)
        enddo

        ! As a helper thing, I want to know how much memory each kind of pair consumes.
        call mem%allocate(sizeofpair,[p%n_species,p%n_species],scalable=.false.,persistent=.false.)
        sizeofpair=0
        do s1=1,p%n_species
        do s2=1,p%n_species
            sizeofpair(s1,s2)=basis%species(s1)%n_basis*basis%species(s2)%n_basis
        enddo
        enddo

        ! Count the number of pairs.
        call mem%allocate(offset,mw%n,persistent=.false.,scalable=.false.)
        call mem%allocate(rnkctr,mw%n,persistent=.false.,scalable=.false.)
        offset=0
        rnkctr=0
        ctr=0
        l=0
        do i=1,dt%np
        do j=1,dt%particle(i)%n
            ! Make it MPI parallel mostly to make sure it is synced
            l=l+1
            if ( mod(l,mw%n) .ne. mw%r ) cycle
            cutoff=atomcutoff(i)+atomcutoff(dt%particle(i)%ind(j))
            if ( dt%particle(i)%d(j) .lt. cutoff ) then
                ! Counter for the number of pairs
                ctr=ctr+1
            endif
        enddo
        enddo
        ! Get the offset per rank
        rnkctr(mw%r+1)=ctr
        call mw%allreduce('sum',rnkctr)
        call mw%allreduce('sum',ctr)
        l=0
        do i=1,mw%n
            offset(i)=l
            l=l+rnkctr(i)
        enddo

        ! Make a little space.
        if ( ctr .gt. 10000000 ) then
            ! That is more than 10 million pairs. That is something like
            ! 400MB. Should probably throw a warning/error and deal with
            ! that in some other way. I suppose there are bigger issues
            ! before I hit that wall though.
            call rl_stop_gracefully(['There are more than 10 million pairs. Everything will break. Tell Olle to think about that.'],&
                                    rl_exitcode_symmetry,mw%comm)
        endif
        call mem%allocate(dr,[3,ctr],persistent=.false.,scalable=.false.)
        call mem%allocate(di,[2,ctr],persistent=.false.,scalable=.false.)

        dr=0.0_r8
        di=0
        ! Then store the relevant pairs.
        ctr=0
        l=0
        do i=1,dt%np
        do j=1,dt%particle(i)%n
            ! Make it MPI parallel mostly to make sure it is synced
            l=l+1
            if ( mod(l,mw%n) .ne. mw%r ) cycle
            cutoff=atomcutoff(i)+atomcutoff(dt%particle(i)%ind(j))
            if ( dt%particle(i)%d(j) .lt. cutoff ) then
                ! Store vector of pairs and indices
                ctr=ctr+1
                dr(:,ctr+offset(mw%r+1))=dt%particle(i)%v(:,j)
                di(:,ctr+offset(mw%r+1))=[i,dt%particle(i)%ind(j)]
            endif
        enddo
        enddo
        ! And sum it up over ranks.
        call mw%allreduce('sum',dr)
        call mw%allreduce('sum',di)

        ! A little cleanup
        call mem%deallocate(offset,persistent=.false.,scalable=.false.)
        call mem%deallocate(rnkctr,persistent=.false.,scalable=.false.)
        call mem%deallocate(atomcutoff,scalable=.false.,persistent=.false.)
        call dt%destroy()

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... fetched list of relevant pairs (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block init

    ! First I need a list of pairs, reduced by symmetry, to be able to make
    ! any kind of sensible decision on storage and strategy.
    select case(sym%n_operation)
    case(1)
    ! This means only the identity symmetry operation is valid, so we don't have to bother
    ! with a lot of fancy stuff. Treating this as a special case also removes a lot of corner
    ! cases from the normal symmetry reduction, so it's a win-win. There is still the transposition
    ! symmetry to consider: H_ij = H_ji, so that has to be worked out properly.
    !
    ! There is a lot of copy-paste between here and the procedure below, but that is to be a bit
    ! future-safe, in the case of a bajillion atoms I might have to get a bit more clever about
    ! how I do this and not just brute-force it. I'm just going to assume that in the case of
    ! fadrillions of atoms there will not be any particular symmetry.
    !symreduceident: block

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... using transposition symmetry to reduce pairs'
        endif

        ! I will make use of the symmetry helper again, but only partially. The
        ! only thing it will be used for is to keep track of transposes. I use it
        ! to avoid the O(N^2) step in finding the transpose of each pair.
        npair=size(dr,2)
        allocate(sh(p%n_atom))
        do i=1,p%n_atom
            sh(i)%n=0
        enddo
        do i=1,npair
            j=di(1,i)
            sh(j)%n=sh(j)%n+1
        enddo
        do i=1,p%n_atom
            call mem%allocate(sh(i)%v,[3,sh(i)%n],persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%w,[3,sh(i)%n],persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%ind,  sh(i)%n,persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%jnd,  sh(i)%n,persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%a2,   sh(i)%n,persistent=.false.,scalable=.false.)
            sh(i)%v=0.0_r8
            sh(i)%w=0.0_r8
            sh(i)%ind=0
            sh(i)%jnd=0
            sh(i)%a2=0
            sh(i)%n=0
        enddo
        do i=1,npair
            j=di(1,i)
            sh(j)%n=sh(j)%n+1
            sh(j)%v(:,sh(j)%n)=dr(:,i)
            sh(j)%ind(sh(j)%n)=i
            sh(j)%a2(sh(j)%n)=di(2,i)
        enddo
        allocate(vb(p%n_atom))
        do i=1,p%n_atom
            call vb(i)%generate(sh(i)%v,[7,7,7],mem)
        enddo

        ! First thing I will do is to locate the transpose of each vector, that will
        ! always come handy, and is always relevant.
        call mem%allocate(tri,npair,persistent=.false.,scalable=.true.)
        tri=0
        do i=1,npair
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            v0=-dr(:,i) ! backwards vector
            j=di(2,i)   ! the atom this pair points to
            ! Locate backwards vector
            k=vb(j)%locate(sh(j)%v,v0)
            if ( k .gt. 0 ) then
                tri(i)=sh(j)%ind(k)
            else
                call rl_stop_gracefully(['Could not locate transpose of vector, should never happen.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo
        call mw%allreduce('sum',tri)
        ! Good. Now tri(i) holds the index to the transpose of each vector. Cleanup the
        ! temporary things I used to get the transpose index thing.
        do i=1,p%n_atom
            call vb(i)%destroy(mem)
            call mem%deallocate(sh(i)%v,  persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%w,  persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%ind,persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%jnd,persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%a2, persistent=.false.,scalable=.false.)
        enddo
        deallocate(vb)
        deallocate(sh)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... mapped transpose (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Start counting unique pairs, and size of storage:
        j=0
        f0=0.0_r8
        do i=1,npair
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            if ( tri(i) .ge. i ) then
                j=j+1
                s1=p%species( di(1,i) )
                s2=p%species( di(2,i) )
                f0=f0+sizeofpair(s1,s2)
            endif
        enddo
        call mw%allreduce('sum',j)
        call mw%allreduce('sum',f0)
        f0=f0*(4*8*KS%n_spin)/1024_r8**2  ! This should get me some approximate number in MiB
        select type(KS)
        type is(rl_kspace_eigenproblem_multiproc)
            ! In multi-proc mode, I check the max memory per rank?
            f1=f0/KS%ml%n
            call mw%allreduce('max',f1)
        class default
            ! Sanity check. It is weird to have single-proc mode for the KS eigenproblem
            ! and a distributed Hamiltonian. That seems like a pointless middle-ground
            ! something that makes little sense in the end. It's probably better to stop
            ! and tell the user to do something better.
            f1=rl_huge
            f0=0.0_r8
        end select

        ! Use the approximate size to decide on storage:
        if ( f0 .lt. storeall_safely ) then
            ! This means that I can safely store a copy of the realspace data on each rank.
            allocate(rl_realspace_matrix_notdistributed::rmtx)
        elseif ( f1 .lt. storeall_safely) then
            ! This means every group of ranks that deal with one eigenproblem holds a full copy.
            ! It is half-distributed if you will, and avoids annoying all-to-all communications.
            allocate(rl_realspace_matrix_mediumdistributed::rmtx)
        else
            ! Guess we have no choice except to distribute it completely over all ranks.
            allocate(rl_realspace_matrix_fulldistributed::rmtx)
        endif

        ! Start storing the pairs. Will depend on distribution scheme.
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! This is not distributed at all, everyone has a copy of the realspace
            ! Hamiltonian. Start by counting pairs.
            m%n_full_pair = npair
            m%n_irr_pair = 0
            do i=1,npair
                if ( tri(i) .ge. i ) m%n_irr_pair=m%n_irr_pair+1
            enddo
            allocate(m%full_pair(m%n_full_pair))
            allocate(m%irr_pair(m%n_irr_pair))

            ! Start populating. This should populate
            ! the irreducible pairs.
            j=0
            do i=1,m%n_full_pair
                if ( tri(i) .lt. i ) cycle
                j=j+1
                m%irr_pair(j)%a1=di(1,i)
                m%irr_pair(j)%a2=di(2,i)
                m%irr_pair(j)%v=dr(:,i)
                if ( tri(i) .eq. i ) then
                    ! self-term, does not unfold to transpose
                    m%irr_pair(j)%n_unfold_pair = 1
                    allocate(m%irr_pair(j)%unfold_index( 1 ))
                    allocate(m%irr_pair(j)%unfold_operation( 1 ))
                    m%irr_pair(j)%unfold_index=i
                    m%irr_pair(j)%unfold_operation=1
                    m%full_pair(i)%index_unique=j
                    m%full_pair(i)%index_operation=1
                else
                    ! normal term, unfolds to it's transpose
                    m%irr_pair(j)%n_unfold_pair = 2
                    allocate(m%irr_pair(j)%unfold_index( 2 ))
                    allocate(m%irr_pair(j)%unfold_operation( 2 ))
                    m%irr_pair(j)%unfold_index=[i,tri(i)]
                    m%irr_pair(j)%unfold_operation=[1,-1]
                    m%full_pair(i)%index_unique=j
                    m%full_pair(i)%index_operation=1
                    m%full_pair(tri(i))%index_unique=j
                    m%full_pair(tri(i))%index_operation=-1
                endif
            enddo

            ! And the rest of the missing information for the full
            ! pairs that was not set above.
            do i=1,m%n_full_pair
                m%full_pair(i)%a1=di(1,i)
                m%full_pair(i)%a2=di(2,i)
                m%full_pair(i)%v=dr(:,i)
                m%full_pair(i)%index_transpose=tri(i)
            enddo

            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... built '//tochar(m%n_irr_pair)//' non-distributed unique pairs out of '//tochar(m%n_full_pair)//' (',tochar(t1-t0),'s)'
                t0=t1
            endif

        type is(rl_realspace_matrix_mediumdistributed)
            ! Distribute the full Hamiltonian across the groups that deal with the same eigenproblem
            select type(KS)
            type is(rl_kspace_eigenproblem_multiproc)

                ! Count irreducible and full pairs per rank.
                j=0
                k=0
                do i=1,npair
                    if ( mod(i,KS%ml%n) .ne. KS%ml%r ) cycle
                    if ( tri(i) .lt. i ) cycle
                    j=j+1
                    if ( i .eq. tri(i) ) then
                        k=k+1
                    else
                        k=k+2
                    endif
                enddo
                m%n_irr_pair=j
                m%n_full_pair=k

                ! Make space for pairs if there are any. There is a very high risk that there are
                ! pairs on each rank, but in extreme corners cases, maybe not.
                if ( m%n_irr_pair .gt. 0 ) then
                    allocate(m%irr_pair(m%n_irr_pair))
                    allocate(m%full_pair(m%n_full_pair))
                endif

                ! Now go through this guy again and store pairs:
                j=0
                k=0
                do i=1,npair
                    if ( mod(i,KS%ml%n) .ne. KS%ml%r ) cycle
                    if ( tri(i) .lt. i ) cycle
                    j=j+1 ! Counter for irreducible
                    m%irr_pair(j)%a1=di(1,i)
                    m%irr_pair(j)%a2=di(2,i)
                    m%irr_pair(j)%v=dr(:,i)
                    if ( tri(i) .eq. i ) then
                        k=k+1 ! Counter for full pairs
                        ! Store the full pair
                        m%full_pair(k)%a1=di(1,i)
                        m%full_pair(k)%a2=di(2,i)
                        m%full_pair(k)%v=dr(:,i)
                        m%full_pair(k)%index_unique=j
                        m%full_pair(k)%index_operation=1
                        m%full_pair(k)%index_transpose=k
                        ! Store how it unfolds
                        m%irr_pair(j)%n_unfold_pair = 1
                        allocate(m%irr_pair(j)%unfold_index( 1 ))
                        allocate(m%irr_pair(j)%unfold_operation( 1 ))
                        m%irr_pair(j)%unfold_index=k
                        m%irr_pair(j)%unfold_operation=1
                    else
                        k=k+1 ! Counter for full pairs
                        ! First store non-transposed
                        m%full_pair(k)%a1=di(1,i)
                        m%full_pair(k)%a2=di(2,i)
                        m%full_pair(k)%v=dr(:,i)
                        m%full_pair(k)%index_unique=j
                        m%full_pair(k)%index_operation=1
                        m%full_pair(k)%index_transpose=k+1
                        k=k+1 ! Counter for full pairs
                        ! Store transposed
                        m%full_pair(k)%a1=di(1,tri(i))
                        m%full_pair(k)%a2=di(2,tri(i))
                        m%full_pair(k)%v=dr(:,tri(i))
                        m%full_pair(k)%index_unique=j
                        m%full_pair(k)%index_operation=-1
                        m%full_pair(k)%index_transpose=k-1
                        ! Store unfold information
                        m%irr_pair(j)%n_unfold_pair = 2
                        allocate(m%irr_pair(j)%unfold_index( 2 ))
                        allocate(m%irr_pair(j)%unfold_operation( 2 ))
                        m%irr_pair(j)%unfold_index=[k-1,k]
                        m%irr_pair(j)%unfold_operation=[1,-1]
                    endif
                enddo

                ! And that's it. Now I just need the offsets to convert between global and local
                ! indices, as well as the global counters.
                call KS%ml%allreduce('sum',m%n_irr_pair,m%n_irr_pair_global)
                call KS%ml%allreduce('sum',m%n_full_pair,m%n_full_pair_global)

                call mem%allocate(dk,[2,KS%ml%n],persistent=.false.,scalable=.false.)
                dk=0
                dk(:,KS%ml%r+1)=[m%n_irr_pair,m%n_full_pair]
                call KS%ml%allreduce('sum',dk)

                ! Some small sanity tests never hurt anyone.
                if ( sum(dk(2,:)) .ne. npair ) then
                    call rl_stop_gracefully(['I have not distributed all full pairs properly'],rl_exitcode_symmetry,mw%comm)
                endif

                ! Determine offsets:
                j=0
                k=0
                do i=1,KS%ml%n
                    if ( KS%ml%r .eq. i-1 ) then
                        m%irr_offset=j
                        m%full_offset=k
                    endif
                    j=j+dk(1,i)
                    k=k+dk(2,i)
                enddo

                ! Some cleanup
                call mem%deallocate(dk,persistent=.false.,scalable=.false.)

                if ( verbosity .gt. 0 ) then
                    t1=mpi_wtime()
                    write(rl_iou,*) '... built '//tochar(m%n_irr_pair_global)//' medium-distributed unique pairs out of '//tochar(m%n_full_pair_global)//' (',tochar(t1-t0),'s)'
                    t0=t1
                endif

            class default
                ! This means single-proc + distributed Hamiltonian, strange middle-ground that should never happen.
                call rl_stop_gracefully(['Distributed Hamitlonian in single-proc mode, makes no sense'],rl_exitcode_param,mw%comm)
            end select
        type is(rl_realspace_matrix_fulldistributed)
            ! Full Hamiltonian distributed everywhere
            ! Should be very similar do the one above, just a different communicator.
            write(*,*) 'fixme fulldistribution'
            call mw%destroy()
            stop
        end select

        ! Cleanup
        call mem%deallocate(sizeofpair,scalable=.false.,persistent=.false.)
        call mem%deallocate(tri,persistent=.false.,scalable=.true.)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.)
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    !end block symreduceident
    case default
    ! This means there is more than one symmetry operation, and we have to be a bit more careful
    ! when reducing it with symmetry
    !symreducenormal: block

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '... using '//tochar(sym%n_operation)//' symmetry operations to reduce pairs'
        endif

        ! Build some helper array thingies to make this fast enough
        ! to be negligible. If one implements this completely naively
        ! it can take very long. The extra helpers is to make it O(N).
        npair=size(dr,2)
        allocate(sh(p%n_atom))
        do i=1,p%n_atom
            sh(i)%n=0
        enddo
        do i=1,npair
            j=di(1,i)
            sh(j)%n=sh(j)%n+1
        enddo
        do i=1,p%n_atom
            call mem%allocate(sh(i)%v,[3,sh(i)%n],persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%w,[3,sh(i)%n],persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%ind,  sh(i)%n,persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%jnd,  sh(i)%n,persistent=.false.,scalable=.false.)
            call mem%allocate(sh(i)%a2,   sh(i)%n,persistent=.false.,scalable=.false.)
            sh(i)%v=0.0_r8
            sh(i)%w=0.0_r8
            sh(i)%ind=0
            sh(i)%jnd=0
            sh(i)%a2=0
            sh(i)%n=0
        enddo
        do i=1,npair
            j=di(1,i)
            sh(j)%n=sh(j)%n+1
            sh(j)%v(:,sh(j)%n)=dr(:,i)
            sh(j)%ind(sh(j)%n)=i
            sh(j)%a2(sh(j)%n)=di(2,i)
        enddo

        ! So, now I have a list of vectors for each atom, stored in sh(atom), together
        ! with the index in the long list of pairs. A symmetry operation will take an
        ! atom to another atom, and all pairs to some other pairs. The helper type makes
        ! it O(N) to figure out which guy goes to which guy with what operation with the
        ! help of Verlet boxes that I generate below:
        allocate(vb(p%n_atom))
        do i=1,p%n_atom
            call vb(i)%generate(sh(i)%v,[7,7,7],mem)
        enddo

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... prepared symmetryreduction (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! First thing I will do is to locate the transpose of each vector, that will
        ! always come handy, and is always relevant.
        call mem%allocate(tri,npair,persistent=.false.,scalable=.true.)
        tri=0
        do i=1,npair
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            v0=-dr(:,i) ! backwards vector
            j=di(2,i)   ! the atom this pair points to
            ! Locate backwards vector
            k=vb(j)%locate(sh(j)%v,v0)
            if ( k .gt. 0 ) then
                tri(i)=sh(j)%ind(k)
            else
                call rl_stop_gracefully(['Could not locate transpose of vector, should never happen.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo
        call mw%allreduce('sum',tri)
        ! Good. Now tri(i) holds the index to the transpose of each vector.

        ! Go through all operations and see how the pairs get permuted.
        call mem%allocate(permutation,[npair,sym%n_operation*2],persistent=.false.,scalable=.false.)
        permutation=0
        l=0
        do iop=1,sym%n_operation
        do ia=1,p%n_atom
           ! Make it parallel
           l=l+1
           if ( mod(l,mw%n) .ne. mw%r ) cycle
           ! This means that atom ia gets transformed to atom ja with this operation
           ja=sym%op(iop)%fmap(ia)
           ! Rotate the vectors, but with a small sanity test attached.
           if ( sh(ia)%n .eq. sh(ja)%n ) then
               call dgemm('N','N',3,sh(ia)%n,3,1.0_r8,sym%op(iop)%m,3,sh(ia)%v,3,0.0_r8,sh(ja)%w,3)
           else
               call rl_stop_gracefully(['Invalid symmetry operation, equivalent atoms have different number of neighbours.'],rl_exitcode_symmetry)
           endif
           ! Now see which pair the rotated pairs become
           do i=1,sh(ja)%n
               j=vb(ja)%locate(sh(ja)%v,sh(ja)%w(:,i))
               if ( j .gt. 0 ) then
                   ! Found the transformed pair! Make a note of it in the permutation array.
                   permutation( sh(ia)%ind(i),iop )=sh(ja)%ind(j)
                   ! Also make a note of where the transpose is!
                   k=tri( sh(ja)%ind(j) )
                   permutation( sh(ia)%ind(i),iop+sym%n_operation )=k
               else
                   ! This means the pair could not be located. Throw error, I think.
                   call rl_stop_gracefully(['Invalid symmetry operation, could not locate transformed pair.'],rl_exitcode_symmetry)
               endif
           enddo
        enddo
        enddo
        ! Add together over ranks
        call mw%allreduce('sum',permutation)

        ! ! SLOW sanity test. Only enably for debugging.
        ! dumdum: block
        !     real(r8), dimension(3) :: v0,v1,v2
        !     integer :: p1,p2
        !
        !     do iop=1,sym%n_operation
        !     do i=1,npair
        !         ! Normal
        !         j=permutation(i,iop)
        !         v0=dr(:,i)
        !         v1=dr(:,j)
        !         v2=matmul(sym%op(iop)%m,v0)-v1
        !         s1=p%species( di(1,i) )
        !         s2=p%species( di(2,i) )
        !         p1=p%species( di(1,j) )
        !         p2=p%species( di(2,j) )
        !         if ( rl_sqnorm(v2) .gt. rl_sqtol ) then
        !             call rl_stop_gracefully(['Invalid symmetry operation.'],rl_exitcode_symmetry)
        !         endif
        !         if ( abs(s1-p1)+abs(s2-p2) .gt. 0 ) then
        !             call rl_stop_gracefully(['Invalid symmetry operation.'],rl_exitcode_symmetry)
        !         endif
        !         ! Not normal
        !         j=permutation(i,iop+sym%n_operation)
        !         v0=dr(:,i)
        !         v1=dr(:,j)
        !         v2=matmul(sym%op(iop)%m,v0)+v1
        !         s1=p%species( di(1,i) )
        !         s2=p%species( di(2,i) )
        !         p1=p%species( di(1,j) )
        !         p2=p%species( di(2,j) )
        !         if ( rl_sqnorm(v2) .gt. rl_sqtol ) then
        !             call rl_stop_gracefully(['Invalid perm+symmetry operation.'],rl_exitcode_symmetry)
        !         endif
        !         if ( abs(s2-p1)+abs(s1-p2) .gt. 0 ) then
        !             call rl_stop_gracefully(['Invalid perm+symmetry operation.'],rl_exitcode_symmetry)
        !         endif
        !     enddo
        !     enddo
        ! end block dumdum

        ! A permutation array can be turned into coordination shells! Do that.
        call rl_coordination_shells_from_permutation_list(permutation,shell_ctr,shell_member,shell_index,mem,prot_index,mw=mw)
        call mem%deallocate(permutation,persistent=.false.,scalable=.false.)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... constructed shells (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! So, now I can try to decide how I am to distribute things. I will do a crude measure
        ! how much memory it cost to store the irreducible pairs, and based on that choose the
        ! distribution scheme. This should probably be revised into something smarter eventually.
        f0=0.0_r8
        do i=1,size(prot_index) ! loop over irreducible
            s1=p%species( di(1,prot_index(i)) )
            s2=p%species( di(2,prot_index(i)) )
            f0=f0+sizeofpair(s1,s2)
        enddo
        f0=f0*(4*8*KS%n_spin)/1024_r8**2  ! This should get me some approximate number in MiB
        select type(KS)
        type is(rl_kspace_eigenproblem_multiproc)
            ! In multi-proc mode, I check the max memory per rank?
            f1=f0/KS%ml%n
            call mw%allreduce('max',f1)
        class default
            ! Sanity check. It is weird to have single-proc mode for the KS eigenproblem
            ! and a distributed Hamiltonian. That seems like a pointless middle-ground
            ! something that makes little sense in the end. It's probably better to stop
            ! and tell the user to do something better.
            f1=rl_huge
            f0=0.0_r8
        end select

        ! This warrants some explanation. f0 holds (approximately) the amount of space needed
        ! to store the realspace Hamiltonian and some other stuff. Probably not very exact, but
        ! at least it's a number directly proportional to the real number, of the correct order
        ! of magnitude. I have set a threshold of how much I can afford to allocate per rank.
        ! This should be enough information to decide on a distribution.
        if ( f0 .lt. storeall_safely ) then
            ! This means that I can safely store a copy of the realspace data on each rank.
            allocate(rl_realspace_matrix_notdistributed::rmtx)
        elseif ( f1 .lt. storeall_safely) then
            ! This means every group of ranks that deal with one eigenproblem holds a full copy.
            ! It is half-distributed if you will, and avoids annoying all-to-all communications.
            allocate(rl_realspace_matrix_mediumdistributed::rmtx)
        else
            ! Guess we have no choice except to distribute it completely over all ranks.
            allocate(rl_realspace_matrix_fulldistributed::rmtx)
        endif

        ! Start storing the pairs. Will depend on distribution scheme.
        select type(m=>rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! Here we have the full Hamiltonian on every rank.
            m%n_full_pair = npair
            m%n_irr_pair = size(prot_index)
            allocate(m%full_pair(m%n_full_pair))
            allocate(m%irr_pair(m%n_irr_pair))

            ! Start populating
            do i=1,m%n_full_pair
                m%full_pair(i)%a1=di(1,i)
                m%full_pair(i)%a2=di(2,i)
                m%full_pair(i)%v=dr(:,i)
                m%full_pair(i)%index_unique=shell_index(i)
                m%full_pair(i)%index_transpose=tri(i)
                ! Will figure out operation later.
                m%full_pair(i)%index_operation=0
            enddo
            do j=1,m%n_irr_pair
                i=prot_index(j)
                m%irr_pair(j)%a1=di(1,i)
                m%irr_pair(j)%a2=di(2,i)
                m%irr_pair(j)%v=dr(:,i)
                m%irr_pair(j)%n_unfold_pair = shell_ctr(j)
                allocate(m%irr_pair(j)%unfold_index( shell_ctr(j) ))
                allocate(m%irr_pair(j)%unfold_operation( shell_ctr(j) ))
                m%irr_pair(j)%unfold_index=shell_member(1:shell_ctr(j),j)
                ! Will figure out the operations below
                m%irr_pair(j)%unfold_operation=0
            enddo

            ! Then the backwards list, how to unfold the irreducible to the full. In theory, this
            ! could have been done when I worked out the shells above, but it serves as a really
            ! good unit test to do it once again in a slightly different way.
            l=maxval(shell_ctr)
            call mem%allocate(dk,[l,m%n_irr_pair],persistent=.false.,scalable=.false.)
            dk=0
            do i=1,m%n_irr_pair
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                ! The pair vector for the irreducible pair
                v0=m%irr_pair(i)%v
                ! Species of the irreducible pair
                s1=p%species( m%irr_pair(i)%a1 )
                s2=p%species( m%irr_pair(i)%a2 )
                ! Here is a sensible place to really make sure that
                ! transposition symmetry holds. Hmm Hmm Double Hmm.
                do j=1,m%irr_pair(i)%n_unfold_pair
                    ! Skip if we have already sorted out this connection
                    if ( dk(j,i) .ne. 0 ) cycle
                    ! Find the operation
                    k=m%irr_pair(i)%unfold_index(j)
                    ! Pair vector of unfolded pair
                    v1=m%full_pair(k)%v
                    p1=p%species( m%full_pair(k)%a1 )
                    p2=p%species( m%full_pair(k)%a2 )
                
                    l=0
                    opl1: do iop=1,sym%n_operation                       
                        if ( rl_sqnorm(matmul(sym%op(iop)%m,v0)-v1) .lt. rl_sqtol ) then    
                            ! Don't transpose unless I say so!
                            if ( abs(s1-p1)+abs(s2-p2) .gt. 0 ) cycle
                            l=iop
!if ( abs(s1-p1)+abs(s2-p2) .gt. 0 ) then
!    write(rl_iou,*) 'OOOPS OPERATION'
!    write(rl_iou,*) s1,s2,p1,p2
!    write(rl_iou,*) 'v0',v0 
!    write(rl_iou,*) 'v1',v1
!    
!endif
                            exit opl1
                        endif
                    enddo opl1
                    ! Then test with transpose
                    if ( l .eq. 0 ) then
                        opl2: do iop=1,sym%n_operation
                            if ( rl_sqnorm(matmul(sym%op(iop)%m,v0)+v1) .lt. rl_sqtol ) then
                                ! Here the species must be transposed!
                                if ( abs(s1-p2)+abs(s2-p1) .gt. 0 ) cycle
                                l=-iop
                                exit opl2
                            endif
                        enddo opl2
                    endif
                    if ( l .ne. 0 ) then
                        dk(j,i)=l
                    else
                        call rl_stop_gracefully(['I do not understand vectors.'],rl_exitcode_symmetry,mw%comm)
                    endif
                    ! Now, it makes perfect sense to find the transpose right away:
!                    ufl2: do ii=1,m%irr_pair(i)%n_unfold_pair
!                        if ( ii .eq. j ) cycle
!                        if ( m%irr_pair(i)%unfold_index(ii) .eq. m%full_pair(k)%index_transpose ) then
!                            dk(ii,i)=-l
!                            exit ufl2
!                        endif
!                    enddo ufl2
                enddo
            enddo
            call mw%allreduce('sum',dk)

            do i=1,m%n_irr_pair
            do j=1,m%irr_pair(i)%n_unfold_pair
                l=dk(j,i)
                if ( l .ne. 0 ) then
                    m%irr_pair(i)%unfold_operation(j)=l
                else
                    call rl_stop_gracefully(['I do not understand vectors.'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo
            enddo
            call mem%deallocate(dk,persistent=.false.,scalable=.false.)

            ! It is fast enough to get the inverse of the list of operations
            ! above, how to unfold if we are at some full pair.
            do i=1,m%n_irr_pair
            do j=1,m%irr_pair(i)%n_unfold_pair
                k=m%irr_pair(i)%unfold_index(j)
                l=m%irr_pair(i)%unfold_operation(j)
                m%full_pair(k)%index_unique=i
                m%full_pair(k)%index_operation=l
            enddo
            enddo

            ! Report that we are done with this step
            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... built '//tochar(m%n_irr_pair)//' non-distributed unique pairs out of '//tochar(m%n_full_pair)//' (',tochar(t1-t0),'s)'
                t0=t1
            endif
        type is(rl_realspace_matrix_mediumdistributed)
            ! This is the case where I want to distribute the realspace matrices
            ! over the sub-communicators associated with each eigenvalue problem.
            !
            ! Please not that in this case I assume we are in multi-proc mode,
            ! and will throw an error if not.
            select type(KS)
            type is(rl_kspace_eigenproblem_multiproc)

            ! Full Hamiltonian per group that deals with the same eigenproblem. First
            ! I have to decide on a sensible distribution. I will try to distribute
            ! the irreducible such that we have an approximately constant number
            ! of unfolded pairs per rank. The array ind hold the indices to the
            ! irreducible pairs that should end up on this rank.
            call rl_distribute_weighted(KS%ml,shell_ctr,m%n_irr_pair,ind,mem)
            ! Now I should know the number of irreducible on this rank,
            ! get the global number of irreducible:
            call KS%ml%allreduce('sum',m%n_irr_pair,m%n_irr_pair_global)

            if ( m%n_irr_pair .gt. 0 ) then
                ! Sorting out the information once distributed needs a few passes over
                ! all the pairs to get everything right. bear with me, it's not too
                ! complicated, I think. First we store the relevant irreducible pairs,
                ! and in the same pass count the number of full pairs that can be
                ! constructed from those:
                allocate(m%irr_pair(m%n_irr_pair))

                l=0
                do j=1,m%n_irr_pair
                    i=ind(j)        ! global irreducible index
                    i=prot_index(i) ! global index to pair
                    m%irr_pair(j)%a1=di(1,i)
                    m%irr_pair(j)%a2=di(2,i)
                    m%irr_pair(j)%v=dr(:,i)
                    m%irr_pair(j)%n_unfold_pair = shell_ctr( ind(j) )
                    ! prepare space for unfolding information
                    allocate(m%irr_pair(j)%unfold_index( m%irr_pair(j)%n_unfold_pair ))
                    allocate(m%irr_pair(j)%unfold_operation( m%irr_pair(j)%n_unfold_pair ))
                    m%irr_pair(j)%unfold_index=0
                    m%irr_pair(j)%unfold_operation=0
                    l=l+m%irr_pair(j)%n_unfold_pair ! This counts the number of full pairs on this rank.
                enddo

                ! Make space for the full pairs
                m%n_full_pair=l
                allocate(m%full_pair(m%n_full_pair))
                ! Now pass through again, and fill out the list of full pairs:
                l=0
                do j=1,m%n_irr_pair
                    do i=1,m%irr_pair(j)%n_unfold_pair
                        l=l+1 ! Counter for full pairs
                        k=shell_member(i,ind(j)) ! global index to full pair that this irreducible pair can unfold to
                        m%full_pair(l)%a1=di(1,k)
                        m%full_pair(l)%a2=di(2,k)
                        m%full_pair(l)%v=dr(:,k)
                        m%full_pair(l)%index_unique=j    ! Unique index points to rank-local pair
                        m%full_pair(l)%index_operation=0 ! Sort out operations later.
                        ! Fill out the unfold index
                        m%irr_pair(j)%unfold_index(i)=l
                    enddo
                enddo

                write(rl_iou,*) 'FIXME SPECIES MEDIUMDISTRIBUTED:'
                call mw%destroy()
                stop

                ! Get the operations that unfold the irreducible to the full points.
                do i=1,m%n_irr_pair
                    v0=m%irr_pair(i)%v
                    do j=1,m%irr_pair(i)%n_unfold_pair
                        ! Skip if already determined
                        if ( m%irr_pair(i)%unfold_operation(j) .ne. 0 ) cycle
                        k=m%irr_pair(i)%unfold_index(j)
                        v1=m%full_pair(k)%v
                        l=0
                        ! Test with normal operation
                        do iop=1,sym%n_operation
                            if ( rl_sqnorm(matmul(sym%op(iop)%m,v0)-v1) .lt. rl_sqtol ) then
                                l=iop
                                exit
                            endif
                        enddo
                        ! Then test with transpose
                        if ( l .eq. 0 ) then
                            do iop=1,sym%n_operation
                                if ( rl_sqnorm(matmul(sym%op(iop)%m,v0)+v1) .lt. rl_sqtol ) then
                                    l=-iop
                                    exit
                                endif
                            enddo
                        endif
                        if ( l .ne. 0 ) then
                            ! Found the connecting operation
                            m%irr_pair(i)%unfold_operation(j)=l
                        else
                            ! This means there is no connecting operation. Should never happen.
                            call rl_stop_gracefully(['I do not understand vectors.'],rl_exitcode_symmetry,mw%comm)
                        endif
                        ! This ensures transposition symmetry, I think.
                        do ii=1,m%irr_pair(i)%n_unfold_pair
                            if ( ii .eq. j ) cycle
                            if ( m%irr_pair(i)%unfold_index(ii) .eq. m%full_pair(k)%index_transpose ) then
                                m%irr_pair(i)%unfold_operation(ii)=-l
                                exit
                            endif
                        enddo
                    enddo
                enddo
                ! Then build the inverse of this, simple enough:
                do i=1,m%n_irr_pair
                do j=1,m%irr_pair(i)%n_unfold_pair
                    k=m%irr_pair(i)%unfold_index(j)
                    l=m%irr_pair(i)%unfold_operation(j)
                    m%full_pair(k)%index_unique=i
                    m%full_pair(k)%index_operation=l
                enddo
                enddo
                ! That's it, all the pairs are distributed as they should! I hope.
            else
                ! This means n_irr_pair=0, that is there are no pairs on this rank. Could happen in some odd
                ! cases with small systems with very many ranks, I suppose. Not particularly efficient, but
                ! not sure what to do about it. Should maybe throw a warning saying that the user is wasting
                ! their time?
                m%n_irr_pair=0
                m%n_full_pair=0
            endif

            ! At this point, all pairs are distributed as they should be. All that remains is to
            ! keep track of the offsets, such that I can construct a global index from the local
            ! pair indices. Note that this is with respect to the sub-communicator in KS%ml.
            call mem%allocate(dk,[2,KS%ml%n],persistent=.false.,scalable=.false.)
            dk=0
            dk(:,KS%ml%r+1)=[m%n_irr_pair,m%n_full_pair]
            call KS%ml%allreduce('sum',dk)

            ! Some small sanity tests never hurt anyone.
            if ( sum(dk(2,:)) .ne. npair ) then
                call rl_stop_gracefully(['I have not distributed all full pairs properly'],rl_exitcode_symmetry,mw%comm)
            endif
            if ( sum(dk(1,:)) .ne. size(prot_index) ) then
                call rl_stop_gracefully(['I have not distributed all irreducible pairs properly'],rl_exitcode_symmetry,mw%comm)
            endif

            ! Determine offsets:
            j=0
            k=0
            do i=1,KS%ml%n
                if ( KS%ml%r .eq. i-1 ) then
                    m%irr_offset=j
                    m%full_offset=k
                endif
                j=j+dk(1,i)
                k=k+dk(2,i)
            enddo

            ! Some cleanup
            call mem%deallocate(ind,persistent=.true.,scalable=.true.)
            call mem%deallocate(dk,persistent=.false.,scalable=.false.)

            ! Report that we are done with this step
            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... built '//tochar(m%n_irr_pair_global)//' medium-distributed unique pairs out of '//tochar(m%n_full_pair_global)//' (',tochar(t1-t0),'s)'
                t0=t1
            endif
        class default
            ! This means single-proc + distributed Hamiltonian, strange middle-ground that should never happen.
            call rl_stop_gracefully(['Distributed Hamiltonian in single-proc mode, makes no sense'],rl_exitcode_param,mw%comm)
        end select
        type is(rl_realspace_matrix_fulldistributed)
            ! Full Hamiltonian distributed everywhere
            ! Should be very similar do the one above, just a different communicator.
            write(*,*) 'fixme fulldistribution'
            call mw%destroy()
            stop
        end select

        ! Some intermediate cleanup?
        do i=1,p%n_atom
            call vb(i)%destroy(mem)
            call mem%deallocate(sh(i)%v,  persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%w,  persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%ind,persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%jnd,persistent=.false.,scalable=.false.)
            call mem%deallocate(sh(i)%a2, persistent=.false.,scalable=.false.)
        enddo
        deallocate(vb)
        deallocate(sh)
        call mem%deallocate(shell_member,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_ctr,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_index,persistent=.true.,scalable=.false.)
        call mem%deallocate(prot_index,persistent=.true.,scalable=.false.)
        call mem%deallocate(sizeofpair,scalable=.false.,persistent=.false.)
        call mem%deallocate(tri,persistent=.false.,scalable=.true.)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.)
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
    !end block symreducenormal
    end select

    ! No matter how I symmetry-reduced or distributed, I need to sort out
    ! the lattice vectors
    !grablatticevecs: block
        ! No matter how I did it, I have to make sure we have the lattice vectors.
        do i=1,rmtx%n_full_pair
            v0=rmtx%full_pair(i)%v+p%cartesian_coordinate(:,rmtx%full_pair(i)%a1)-p%cartesian_coordinate(:,rmtx%full_pair(i)%a2)
            v0=matmul(p%inv_latticevectors,v0)
            if ( rl_sqnorm(v0-anint(v0)) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['I do not understand vectors.'],rl_exitcode_symmetry,mw%comm)
            else
                rmtx%full_pair(i)%lv=matmul(p%latticevectors,anint(v0))
            endif
        enddo
        do i=1,rmtx%n_irr_pair
            v0=rmtx%irr_pair(i)%v+p%cartesian_coordinate(:,rmtx%irr_pair(i)%a1)-p%cartesian_coordinate(:,rmtx%irr_pair(i)%a2)
            v0=matmul(p%inv_latticevectors,v0)
            if ( rl_sqnorm(v0-anint(v0)) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['I do not understand vectors.'],rl_exitcode_symmetry,mw%comm)
            else
                rmtx%irr_pair(i)%lv=matmul(p%latticevectors,anint(v0))
            endif
        enddo
    !end block grablatticevecs

    ! We should also make note of the basis set used, and make space for Hamiltonians
    ! overlap matrices and density matrices.
    !setstorage: block

        do ipair=1,rmtx%n_full_pair
            ! which two atoms are involved
            a1=rmtx%full_pair(ipair)%a1
            a2=rmtx%full_pair(ipair)%a2
            ! which two species are involved
            s1=p%species(a1)
            s2=p%species(a2)
            ! how many basis functions do these two atoms have
            nb1=basis%species(s1)%n_basis
            nb2=basis%species(s2)%n_basis
            ! store the dimensions of this block
            rmtx%full_pair(ipair)%nb1=nb1
            rmtx%full_pair(ipair)%nb2=nb2
            ! Make space for the indexing, the dummy one where I grab stuff from AIMS
            allocate(rmtx%full_pair(ipair)%idx(nb1,nb2))
            rmtx%full_pair(ipair)%idx=0
        enddo

        ! And for the irreducible part, here we do the same thing, but also
        ! make space for the actual Hamiltonian/Overlap
        do ipair=1,rmtx%n_irr_pair
            ! which two atoms are involved
            a1=rmtx%irr_pair(ipair)%a1
            a2=rmtx%irr_pair(ipair)%a2
            ! which two species are involved
            s1=p%species(a1)
            s2=p%species(a2)
            ! how many basis functions do these two atoms have
            nb1=basis%species(s1)%n_basis
            nb2=basis%species(s2)%n_basis
            ! store the dimensions of this block
            rmtx%irr_pair(ipair)%nb1=nb1
            rmtx%irr_pair(ipair)%nb2=nb2
            ! Make space for the indexing, the dummy one where I grab stuff from AIMS?
            allocate(rmtx%irr_pair(ipair)%idx(nb1,nb2))
            rmtx%irr_pair(ipair)%idx=0
            allocate(rmtx%irr_pair(ipair)%hamiltonian( nb1 , nb2 , KS%n_spin))
            allocate(rmtx%irr_pair(ipair)%densitymatrix( nb1 , nb2 , KS%n_spin))
            allocate(rmtx%irr_pair(ipair)%overlap( nb1 , nb2 ))
            rmtx%irr_pair(ipair)%hamiltonian   = 0.0_r8
            rmtx%irr_pair(ipair)%densitymatrix = 0.0_r8
            rmtx%irr_pair(ipair)%overlap       = 0.0_r8
        enddo
    !end block setstorage

    ! It does not cost too much to sanity-test the symmetries of the graph-based thing
    ! so we might as well do it here? That should not be too expensive. Can skip this
    ! once things stabilize.
    !sanitychecksym: block

        ! Create a verlet-box for the extended cluster
        call vbb%generate(ec%cartesian_coordinate,[15,15,15],mem)

        do i=1,rmtx%n_irr_pair
            ! this checks that folding out from the irreducible pairs works as it should.
            v0=rmtx%irr_pair(i)%v
            s1=p%species( rmtx%irr_pair(i)%a1 )
            s2=p%species( rmtx%irr_pair(i)%a2 )

            do l=1,rmtx%irr_pair(i)%n_unfold_pair
                j=rmtx%irr_pair(i)%unfold_index(l)
                k=rmtx%irr_pair(i)%unfold_operation(l)
                v1=rmtx%full_pair(j)%v
                p1=p%species( rmtx%full_pair(j)%a1 )
                p2=p%species( rmtx%full_pair(j)%a2 )
                if ( k .gt. 0 ) then                
                    f0=rl_sqnorm(v1-matmul(sym%op(k)%m,v0))                     ! Check that vectors match                    
                    f0=f0+abs(s1-p1)+abs(s2-p2)                                 ! Check species
                    f0=f0+abs( rmtx%irr_pair(i)%nb1 - rmtx%full_pair(j)%nb1 )   ! Check dimensions of buffers
                    f0=f0+abs( rmtx%irr_pair(i)%nb2 - rmtx%full_pair(j)%nb2 )
                else
                    f0=rl_sqnorm(v1+matmul(sym%op(-k)%m,v0))
                    f0=f0+abs(s1-p2)+abs(s2-p1)
                    f0=f0+abs( rmtx%irr_pair(i)%nb1 - rmtx%full_pair(j)%nb2 )
                    f0=f0+abs( rmtx%irr_pair(i)%nb2 - rmtx%full_pair(j)%nb1 )
                endif
                if ( f0 .gt. rl_sqtol ) then
                    call rl_stop_gracefully(['Bad symmetry operation irreducible pair list.'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo
            ! Check that we can find this pair in the extended cluster.
            a1=rmtx%irr_pair(i)%a1
            a2=rmtx%irr_pair(i)%a2
            v0=ec%cartesian_coordinate(:,ec%locate_unit_cell_atom_in_cluster(a1))
            v0=v0+rmtx%irr_pair(i)%v
            ie=vbb%locate(ec%cartesian_coordinate,v0)
            if ( ie .lt. 1 ) then
                call rl_stop_gracefully(['Could not find pair in extended cluster.'],rl_exitcode_symmetry,mw%comm)
            else
                if ( ec%index_unit_cell(ie) .ne. a2 ) then
                    call rl_stop_gracefully(['Inconsistent labelling of pair.'],rl_exitcode_symmetry,mw%comm)
                endif
            endif
        enddo

        do i=1,rmtx%n_full_pair
            ! this checks that I can rotate the unique pair to the full pair
            j=rmtx%full_pair(i)%index_unique
            k=rmtx%full_pair(i)%index_operation
            v1=rmtx%full_pair(i)%v
            v0=rmtx%irr_pair(j)%v
            s1=p%species( rmtx%irr_pair(j)%a1 )
            s2=p%species( rmtx%irr_pair(j)%a2 )
            p1=p%species( rmtx%full_pair(i)%a1 )
            p2=p%species( rmtx%full_pair(i)%a2 )
            if ( k .gt. 0 ) then
                f0=rl_sqnorm(v1-matmul(sym%op(k)%m,v0))
                f0=f0+abs(s1-p1)+abs(s2-p2)
                f0=f0+abs( rmtx%irr_pair(j)%nb1 - rmtx%full_pair(i)%nb1 )
                f0=f0+abs( rmtx%irr_pair(j)%nb2 - rmtx%full_pair(i)%nb2 )
            else
                f0=rl_sqnorm(v1+matmul(sym%op(-k)%m,v0))
                f0=f0+abs(s1-p2)+abs(s2-p1)
                f0=f0+abs( rmtx%irr_pair(j)%nb1 - rmtx%full_pair(i)%nb2 )
                f0=f0+abs( rmtx%irr_pair(j)%nb2 - rmtx%full_pair(i)%nb1 )
            endif
            if ( f0 .gt. rl_sqtol ) then
                call rl_stop_gracefully(['Bad symmetry operation in full pair list.'],rl_exitcode_symmetry,mw%comm)
            endif
            ! check that the transpose is really the transpose
            j=rmtx%full_pair(i)%index_transpose
            v0=rmtx%full_pair(j)%v
            f0=rl_sqnorm(v1+v0)
            if ( f0 .gt. rl_sqtol ) then
                call rl_stop_gracefully(['Bad transpose in full pair list.'],rl_exitcode_symmetry,mw%comm)
            endif
            ! Check that we can find this pair in the extended cluster.
            a1=rmtx%full_pair(i)%a1
            a2=rmtx%full_pair(i)%a2
            v0=ec%cartesian_coordinate(:,ec%locate_unit_cell_atom_in_cluster(a1))
            v0=v0+rmtx%full_pair(i)%v
            ie=vbb%locate(ec%cartesian_coordinate,v0)
            if ( ie .lt. 1 ) then
                call rl_stop_gracefully(['Could not find pair in extended cluster.'],rl_exitcode_symmetry,mw%comm)
            else
                if ( ec%index_unit_cell(ie) .ne. a2 ) then
                    call rl_stop_gracefully(['Inconsistent labelling of pair.'],rl_exitcode_symmetry,mw%comm)
                endif
            endif
        enddo

        call vbb%destroy(mem)
    !end block sanitychecksym

    ! So, now I have a certain distribution of the realspace storage. This is as good
    ! a time as any to sort out the indexing. The indexing I talk about here is for the
    ! Fourier transform. The realspace graph-based storage I have as pairs is not useful
    ! for throwing into any kind of linalg solver. I'm not completely clear on how to
    ! store this information without branching everything brutally and/or making things
    ! too slow. I'm reasonably sure that the Fourier transform step is not a particular
    ! bottleneck in most calculations, and can afford to be a little sloppy here. Might
    ! have to revise that later though. Anyway, the indexing will surely depend on how
    ! the KS problem is distributed.
    select type(KS)
    type is(rl_kspace_eigenproblem_singleproc)
    !indexingsingleproc: block
        ! For single-proc mode the indexing is pretty much non-existents. Just keep
        ! track of which atoms and the number of basis functions and all that.
        do ipair=1,rmtx%n_full_pair
            ! Make sure to not that we don't work with the fancy indexing:
            rmtx%full_pair(ipair)%n_ft_element=-rl_hugeint
            allocate(rmtx%full_pair(ipair)%ft_element(1,1))
            rmtx%full_pair(ipair)%ft_element=-rl_hugeint
        enddo
    !end block indexingsingleproc
    type is(rl_kspace_eigenproblem_multiproc)
    !indexingmultiproc: block

        ! Now, I'm pretty sure the indexing thingy does not depend on how the realspace
        ! things are distributed. Maybe. Also not sure if it's worth doing this in parallel
        ! in the case of full distribution, or if that only will add to the confusion.
        do ipair=1,rmtx%n_full_pair
            ! which two atoms are involved
            a1=rmtx%full_pair(ipair)%a1
            a2=rmtx%full_pair(ipair)%a2
            ! which two species are involved
            s1=p%species(a1)
            s2=p%species(a2)
            ! how many basis functions do these two atoms have
            nb1=basis%species(s1)%n_basis
            nb2=basis%species(s2)%n_basis
            ! what is the offset for the basis functions, to get a global index thing
            ob1=basis%offset(a1)
            ob2=basis%offset(a2)
            ! Now for the real indexing:
            ctr=0
            do b1=1,nb1
            do b2=1,nb2
                gi=b1+ob1 ! global row index
                gj=b2+ob2 ! global column index
                ! See if this element is relevant on this pair?
                if ( KS%bl%rank_from_global_indices(KS%n_basis,KS%n_basis,KS%buf%blocksize,gi,gj) .ne. KS%ml%r ) cycle
                ! Count number of elements for this rank.
                ctr=ctr+1
            enddo
            enddo
            ! Make space for element indexing
            rmtx%full_pair(ipair)%n_ft_element=ctr
            if ( ctr .gt. 0 ) then
                allocate(rmtx%full_pair(ipair)%ft_element(4,ctr))
                rmtx%full_pair(ipair)%ft_element=0
            else
                allocate(rmtx%full_pair(ipair)%ft_element(1,1))
                rmtx%full_pair(ipair)%ft_element=-rl_hugeint
            endif
            ! Store the local indices
            ctr=0
            do b1=1,nb1
            do b2=1,nb2
                gi=b1+ob1 ! global row index
                gj=b2+ob2 ! global column index
                if ( KS%bl%rank_from_global_indices(KS%n_basis,KS%n_basis,KS%buf%blocksize,gi,gj) .ne. KS%ml%r ) cycle
                ctr=ctr+1
                li=KS%bl%local_row_from_global(KS%buf%blocksize,gi)
                lj=KS%bl%local_col_from_global(KS%buf%blocksize,gj)
                rmtx%full_pair(ipair)%ft_element(:,ctr)=[b1,b2,li,lj]
            enddo
            enddo
        enddo
    !end block indexingmultiproc
    end select

    ! Check if I know how to type 'deallocate' for every time I type 'allocate'. Harder than you would think.
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

    ! Report that we are done with this step
    if ( verbosity .gt. 0 ) then
        t1=mpi_wtime()
        write(rl_iou,*) '... worked out indexing for Fourier transform (',tochar(t1-t0),'s)'
        t0=t1
    endif
end subroutine

!> Construct the map from the hamiltonian variable thing in aims to my format. Will have to update this to agree with all possible ways of indexing and distributing the Hamiltonian
subroutine build_index(rmtx,p,basis,mw,mem,&
        coords_center, center_to_cell, center_to_atom, &
        index_hamiltonian, column_index_hamiltonian, cbasis_to_basis, &
        cbasis_to_center, centers_basis_integrals, frac_coords, lattice_vector,&
        position_in_hamiltonian,&
        verbosity)
    !> realspace hamiltonian
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis
    type(rl_lcao_basis_set), intent(in) :: basis
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> AIMS coords_centers, Where is center i
    real(r8), intent(in), dimension(:,:) :: coords_center
    !> AIMS center_to_cell. What vector moves center i back to the unit cell
    integer, intent(in), dimension(:) :: center_to_cell
    !> AIMS center_to_atom, what index in the unit cell is center i?
    integer, intent(in), dimension(:) :: center_to_atom
    !> AIMS index Hamiltonian
    integer, intent(in), dimension(:,:,:) :: index_hamiltonian
    !> AIMS column index Hamiltonian
    integer, intent(in), dimension(:) :: column_index_hamiltonian
    !> AIMS indices that map the all the basis functions to something I can understand.
    integer, intent(in), dimension(:) :: cbasis_to_basis
    integer, intent(in), dimension(:) :: cbasis_to_center
    integer, intent(in), dimension(:) :: centers_basis_integrals
    !> AIMS coordinates of atoms
    real(r8), intent(in), dimension(:,:) :: frac_coords
    real(r8), intent(in), dimension(3,3) :: lattice_vector
    !> AIMS position_in_hamiltonian
    integer, dimension(:,:), intent(in) :: position_in_hamiltonian
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:,:), allocatable :: crd
    real(r8) :: timer,t0,t1
    !init: block
        integer :: i
    !matchvectors: block
        type(rl_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(3) :: v0
        !integer :: i,j,k,a1,a2,s1,s2,b1,b2,bb1,bb2,ipair,celli,hi,ctr
        integer :: j,k,a1,a2,s1,s2,b1,b2,bb1,bb2,ipair,celli,hi,ctr
    !finalize: block
        !integer :: i,j,k,a1,a2,s1,s2,b1,b2,bb1,bb2

    ! Set up some initial things
    !init: block

        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
            t0=timer
            write(rl_iou,*) ''
            write(rl_iou,*) 'BUILDING REALSPACE MATRIX INDEX'
        endif

        ! Long list of sanity checks to make sure things are what I think they are
        if ( p%n_atom .ne. size(frac_coords,2) ) then
            call rl_stop_gracefully(['Inconsistent geometry when building index'],rl_exitcode_param,mw%comm)
        endif
        if ( basis%n_basis .ne. size(index_hamiltonian,3) ) then
            call rl_stop_gracefully(['Inconsistent basis when building index'],rl_exitcode_param,mw%comm)
        endif

        ! In case there already is some indexing, make sure we destroy it
        ! And in case there is nothing, make space for the indices
        do i=1,rmtx%n_full_pair
            rmtx%full_pair(i)%idx=0
        enddo
        do i=1,rmtx%n_irr_pair
            rmtx%irr_pair(i)%idx=0
        enddo

        ! The Cartesian coordinates seem to change during the scf cycle.
        ! Or maybe it is during initialization. Who knows. This routine
        ! should be retired anyway. Ugly workaround to fix it:
        call mem%allocate(crd,[3,p%n_atom],persistent=.false.,scalable=.false.)
        crd=0.0_r8
        do i=1,p%n_atom
            crd(:,i)=matmul(lattice_vector,frac_coords(:,i))
        enddo
    !end block init

    ! Match the AIMS vectors with the ones I have generated
    !matchvectors: block

        ! store relevant center coordinates in a Verlet list, for fast lookup.
        call mem%allocate(dr,[3,size(centers_basis_integrals)],persistent=.false.,scalable=.false.)
        dr=0.0_r8
        do i=1,size(centers_basis_integrals)
            j=centers_basis_integrals(i)
            dr(:,i)=coords_center(:,j)
        enddo
        call vb%generate(dr,[15,15,15],mem)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... sorted points ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! This depends on how the realspace matrices are distributed, and also on how the
        ! AIMS matrix is distributed. Very annoying.
        select type(rmtx)
        type is(rl_realspace_matrix_notdistributed)
            ! All ranks have everything, means I will index in parallel and communicate in the end.
            ctr=0
            do ipair=1,rmtx%n_full_pair
                if ( mod(ipair,mw%n) .ne. mw%r ) cycle
                ! I need the coordinate of atom a2 to be able to locate which cell I should look in.
                a1=rmtx%full_pair(ipair)%a1
                a2=rmtx%full_pair(ipair)%a2
                v0=rmtx%full_pair(ipair)%v+crd(:,a1)
                i=vb%locate(dr,v0)
                if ( i .gt. 0 ) then
                    j=centers_basis_integrals(i)
                    if ( center_to_atom(j) .eq. a2 ) then
                        celli=center_to_cell(j)
                    else
                        call rl_stop_gracefully(['Matched pair wrong'],rl_exitcode_symmetry,mw%comm)
                    endif
                else
                    call rl_stop_gracefully(['Matched pair wrong'],rl_exitcode_symmetry,mw%comm)
                endif

                ! If all went well, now I know that this pair is in cell 'celli'
                do bb2=1,basis%n_basis ! Number of basis functions/unit cell
                    hi=index_hamiltonian(1,celli,bb2)
                    if ( hi .le. 0 ) cycle
                    ! j is the linear index in the Hamiltonian, and "1" above and "2" below
                    ! indicate what I should loop over.
                    j=hi-1
                    do i=hi,index_hamiltonian(2,celli,bb2)
                        j=j+1 ! linear index in Hamiltonian
                        bb1=column_index_hamiltonian(j) ! Basis function, in the unit cell.
                        ! Unflattened indices in the Hamiltonian
                        if ( a1 .ne. cbasis_to_center(bb1) ) cycle
                        if ( a2 .ne. cbasis_to_center(bb2) ) cycle
                        b1=cbasis_to_basis(bb1)-basis%offset(a1)
                        b2=cbasis_to_basis(bb2)-basis%offset(a2)
                        rmtx%full_pair(ipair)%idx(b1,b2)=j
                    enddo
                enddo
                ! Accumulate counter
                ctr=ctr+1
            enddo
            ! Make sure we counted everything properly.
            call mw%allreduce('sum',ctr)
            if ( ctr .ne. rmtx%n_full_pair ) then
                call rl_stop_gracefully(['Could not index all pairs'],rl_exitcode_symmetry,mw%comm)
            endif

            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... fetched some indices ('//tochar(t1-t0)//'s)'
                t0=t1
            endif


            ! Now communicate things to everywhere
            call rl_communicate_hamiltonian_index(rmtx,basis,mw)

            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... distributed indices ('//tochar(t1-t0)//'s)'
                t0=t1
            endif
        class default
            ! This is either medium or full distributed. In any case, I can safely index serially
            ! and not worry about communicating the collected indices. This would have to be revised
            ! for USE_LOCAL_INDEX, but that is too annoying to bother with right now.
            do ipair=1,rmtx%n_full_pair
                ! I need the coordinate of atom a2 to be able to locate which cell I should look in.
                a1=rmtx%full_pair(ipair)%a1
                a2=rmtx%full_pair(ipair)%a2
                v0=rmtx%full_pair(ipair)%v+crd(:,a1)
                i=vb%locate(dr,v0)
                if ( i .gt. 0 ) then
                    j=centers_basis_integrals(i)
                    if ( center_to_atom(j) .eq. a2 ) then
                        celli=center_to_cell(j)
                    else
                        call rl_stop_gracefully(['Matched pair wrong'],rl_exitcode_symmetry,mw%comm)
                    endif
                else
                    call rl_stop_gracefully(['Matched pair wrong'],rl_exitcode_symmetry,mw%comm)
                endif

                ! If all went well, now I know that this pair is in cell 'celli'
                do bb2=1,basis%n_basis ! Number of basis functions/unit cell
                    hi=index_hamiltonian(1,celli,bb2)
                    if ( hi .le. 0 ) cycle
                    ! j is the linear index in the Hamiltonian, and "1" above and "2" below
                    ! indicate what I should loop over.
                    j=hi-1
                    do i=hi,index_hamiltonian(2,celli,bb2)
                        j=j+1 ! linear index in Hamiltonian
                        bb1=column_index_hamiltonian(j) ! Basis function, in the unit cell.
                        ! Unflattened indices in the Hamiltonian
                        if ( a1 .ne. cbasis_to_center(bb1) ) cycle
                        if ( a2 .ne. cbasis_to_center(bb2) ) cycle
                        b1=cbasis_to_basis(bb1)-basis%offset(a1)
                        b2=cbasis_to_basis(bb2)-basis%offset(a2)
                        rmtx%full_pair(ipair)%idx(b1,b2)=j
                    enddo
                enddo
            enddo

            if ( verbosity .gt. 0 ) then
                t1=mpi_wtime()
                write(rl_iou,*) '... collected indices ('//tochar(t1-t0)//'s)'
                t0=t1
            endif
        end select

        ! Some intermediate cleanup
        call vb%destroy(mem)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.)
        call mem%deallocate(crd,persistent=.false.,scalable=.false.)
    !end block matchvectors

    !finalize: block

        ! Fill out the transpose
        do i=1,rmtx%n_full_pair
            a1=rmtx%full_pair(i)%a1
            a2=rmtx%full_pair(i)%a2
            s1=p%species(a1)
            s2=p%species(a2)
            j=rmtx%full_pair(i)%index_transpose
            do b2=1,basis%species(s2)%n_basis
            do b1=1,basis%species(s1)%n_basis
                bb1=rmtx%full_pair(i)%idx(b1,b2)
                bb2=rmtx%full_pair(j)%idx(b2,b1)
                k=max(bb1,bb2) ! should be safe, since the guys that are not matched are set to zero.
                if ( rmtx%full_pair(i)%idx(b1,b2) .eq. 0 ) rmtx%full_pair(i)%idx(b1,b2)=k
                if ( rmtx%full_pair(j)%idx(b2,b1) .eq. 0 ) rmtx%full_pair(j)%idx(b2,b1)=k
            enddo
            enddo
        enddo

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... matched transpositions ('//tochar(t1-t0)//'s)'
            t0=t1
        endif

        ! Also a good idea to index the irreducible perhaps?
        do i=1,rmtx%n_irr_pair
            ! Find the pair that unfolds with the identity operation
            k=0
            do j=1,rmtx%irr_pair(i)%n_unfold_pair
                if ( rmtx%irr_pair(i)%unfold_operation(j) .eq. 1 ) then
                    k=rmtx%irr_pair(i)%unfold_index(j)
                    exit
                endif
            enddo
            if ( k .ne. 0 ) then
                rmtx%irr_pair(i)%idx=rmtx%full_pair(k)%idx
            else
                call rl_stop_gracefully(['Matched pair wrong'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo
    !end block finalize

    ! Check for stupidity
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],   rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],    rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'], rl_exitcode_memory,mw%comm)

    if ( verbosity .gt. 0 ) then
        write(rl_iou,*) 'Done building index (',tochar(mpi_wtime()-timer),'s)'
    endif

end subroutine

!> rotate a pair block of the hamiltonian/overlap/densitymatrix
subroutine rl_rotate_realspace_pair_block(pair,basis,structure,op,forward,transposition,original_block,rotated_block)
    !> realspace hamiltonian pair to be rotated
    class(rl_realspace_matrix_pair), intent(in) :: pair
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: structure
    !> which operation should we rotate with
    type(rl_spacegroup_operation), intent(in) :: op
    !> forward or backward
    logical, intent(in) :: forward
    !> should the block be transposed?
    logical, intent(in) :: transposition
    !> original block
    real(r8), dimension(:,:), intent(in) :: original_block
    !> rotated block
    real(r8), dimension(:,:), intent(out) :: rotated_block

    integer :: a1,a2,s1,s2,nb1,nb2
    !rot1: block
        real(r8), dimension(pair%nb1,pair%nb1) :: rm1
        real(r8), dimension(pair%nb2,pair%nb2) :: rm2
        real(r8), dimension(pair%nb1,pair%nb2) :: dm0
        real(r8), dimension(pair%nb2,pair%nb1) :: dt0

    ! some heuristics
    !init: block
        a1=pair%a1
        a2=pair%a2
        s1=structure%species(a1)
        s2=structure%species(a2)
        nb1=pair%nb1
        nb2=pair%nb2
    !end block init

    ! Rotate. This can and should be made faster with the block-diagonal stuff.
    ! Other options is to rotate a bunch of things at once. Later problem.
    !rot1: block

        ! Grab correct rotation matrices
        if ( forward ) then
            if ( transposition ) then
                call basis%realspace_rotation_matrix_for_one_atom(s1,op%im,rm1)
                call basis%realspace_rotation_matrix_for_one_atom(s2,op%im,rm2)

                ! call basis%realspace_rotation_matrix_for_one_atom(s1,-op%im,rm1)
                ! call basis%realspace_rotation_matrix_for_one_atom(s2,-op%im,rm2)
                ! call realspace_rotation_matrix_for_one_atom(es,a1,-op%im,rm1)
                ! call realspace_rotation_matrix_for_one_atom(es,a2,-op%im,rm2)
            else
                call basis%realspace_rotation_matrix_for_one_atom(s1,op%im,rm1)
                call basis%realspace_rotation_matrix_for_one_atom(s2,op%im,rm2)

                ! call realspace_rotation_matrix_for_one_atom(es,a1,op%im,rm1)
                ! call realspace_rotation_matrix_for_one_atom(es,a2,op%im,rm2)
            endif
        else
            if ( transposition ) then
                call basis%realspace_rotation_matrix_for_one_atom(s1,op%m,rm1)
                call basis%realspace_rotation_matrix_for_one_atom(s2,op%m,rm2)
                ! call basis%realspace_rotation_matrix_for_one_atom(s1,-op%m,rm1)
                ! call basis%realspace_rotation_matrix_for_one_atom(s2,-op%m,rm2)
                ! call realspace_rotation_matrix_for_one_atom(es,a1,-op%m,rm1)
                ! call realspace_rotation_matrix_for_one_atom(es,a2,-op%m,rm2)
            else
                call basis%realspace_rotation_matrix_for_one_atom(s1,op%m,rm1)
                call basis%realspace_rotation_matrix_for_one_atom(s2,op%m,rm2)
                ! call realspace_rotation_matrix_for_one_atom(es,a1,op%m,rm1)
                ! call realspace_rotation_matrix_for_one_atom(es,a2,op%m,rm2)
            endif
        endif

        ! Rotate? @TODO GEMM THIS
        if ( transposition ) then
            ! dt0=transpose(original_block)
            ! dt0=matmul(dt0,transpose(rm1))
            ! dt0=matmul(rm2,dt0)
            ! rotated_block=dt0

            ! dt0=matmul(transpose(original_block),transpose(rm1))
            ! dt0=matmul(rm2,dt0)
            ! rotated_block=dt0
            dm0=matmul(original_block,transpose(rm2))
            dm0=matmul(rm1,dm0)
            rotated_block=transpose(dm0)
        else
            dm0=matmul(original_block,transpose(rm2))
            dm0=matmul(rm1,dm0)
            rotated_block=dm0
        endif
    !end block rot1
end subroutine

!> Spread out information from pairs.
subroutine rl_communicate_hamiltonian_index(rmtx,basis,mw)
    !> extended cluster
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw

    integer, parameter :: sensible_size=500000 ! Should decide this on-the-fly
    integer, dimension(:,:,:), allocatable :: ibf
    integer :: b1,b2,iter,niter,i,l,ihi,ilo,ipair,nb,np

    ! What is a decent size of buffer to use?
    nb=0
    do i=1,size(basis%species)
        nb=max(nb,basis%species(i)%n_basis)
    enddo
    np=ceiling(real(sensible_size,r8)/real(nb*nb,r8))
    np=min(np,rmtx%n_full_pair)
    niter=ceiling(real(rmtx%n_full_pair,r8)/np)

    ! Idiot check? Might remove at some point.
    call mw%check_and_sync(nb,0,bigdeal=2,vname='nbasis')
    call mw%check_and_sync(np,0,bigdeal=2,vname='npair')
    call mw%check_and_sync(niter,0,bigdeal=2,vname='niter')

    allocate(ibf(nb,nb,np))
    ibf=0
    do iter=1,niter
        ilo=(iter-1)*np+1
        ihi=min(iter*np,rmtx%n_full_pair)

        ibf=0
        l=0
        do ipair=ilo,ihi
            if ( ipair .lt. (iter-1)*np ) cycle
            if ( ipair .gt. iter*np ) cycle
            l=l+1
            do b2=1,rmtx%full_pair(ipair)%nb2
            do b1=1,rmtx%full_pair(ipair)%nb1
                ibf(b1,b2,l)=rmtx%full_pair(ipair)%idx(b1,b2)
            enddo
            enddo
        enddo
        call mw%allreduce('sum',ibf)
        l=0
        do ipair=ilo,ihi
            l=l+1
            do b2=1,rmtx%full_pair(ipair)%nb2
            do b1=1,rmtx%full_pair(ipair)%nb1
                rmtx%full_pair(ipair)%idx(b1,b2)=ibf(b1,b2,l)
            enddo
            enddo
        enddo
    enddo
end subroutine

! !> collect the hamiltonian from AIMS. Should be retired.
! subroutine collect_overlap(rmtx,basis,overlap_matrix)
!     !> realspace hamiltonian/overlap/densitymatrix
!     class(rl_realspace_matrix), intent(inout) :: rmtx
!     !> basis set
!     type(rl_lcao_basis_set), intent(in) :: basis
!     !> overlap matrix as it comes from AIMS
!     real(r8), dimension(:), intent(in) :: overlap_matrix
!
!     integer :: i,j,b1,b2
!
!     do i=1,rmtx%n_irr_pair
!         do b1=1,rmtx%irr_pair(i)%nb1
!         do b2=1,rmtx%irr_pair(i)%nb2
!             j=rmtx%irr_pair(i)%idx(b1,b2)
!             if ( j .le. 0) cycle
!             rmtx%irr_pair(i)%overlap(b1,b2)=overlap_matrix(j)
!         enddo
!         enddo
!     enddo
! end subroutine

!> collect the hamiltonian from AIMS. Should be retired.
subroutine collect_overlap(rmtx,basis,p,KS,sym,mw,mem,overlap_matrix)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> kohnsham problem helper
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> overlap matrix as it comes from AIMS
    real(r8), dimension(:), intent(in) :: overlap_matrix

    type(rl_mom_real), dimension(:,:), allocatable :: obuf,rbuf
    type(rl_mom_real), dimension(:,:), allocatable :: abuf
    real(r8) :: pref
    integer :: ipair,jpair,iop,ctr
    integer :: i,j,b1,b2,a1,a2,s1,s2,p1,p2

    ! Maybe make some space for rotation buffers
    allocate(obuf(p%n_species,p%n_species))
    allocate(rbuf(p%n_species,p%n_species))
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%allocate(obuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        call mem%allocate(rbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        obuf(s1,s2)%m=0.0_r8
        rbuf(s1,s2)%m=0.0_r8
    enddo
    enddo

    ! Set irreducible to nothing
    do i=1,rmtx%n_irr_pair
        rmtx%irr_pair(i)%overlap=0.0_r8
    enddo

    select case(sym%n_operation)
    case(1)
        ! Simplest way of collecting, no symmetry to worry about at all.
        do i=1,rmtx%n_irr_pair
            do b1=1,rmtx%irr_pair(i)%nb1
            do b2=1,rmtx%irr_pair(i)%nb2
                j=rmtx%irr_pair(i)%idx(b1,b2)
                if ( j .le. 0) cycle
                rmtx%irr_pair(i)%overlap(b1,b2)=overlap_matrix(j)
            enddo
            enddo
        enddo
    case default
        ! Start collecting where we have to transform the elements
        do ipair=1,rmtx%n_irr_pair
            pref=1.0_r8/real(rmtx%irr_pair(ipair)%n_unfold_pair,r8)
            s1=p%species( rmtx%irr_pair(ipair)%a1 )
            s2=p%species( rmtx%irr_pair(ipair)%a2 )

            do i=1,rmtx%irr_pair(ipair)%n_unfold_pair
                jpair=rmtx%irr_pair(ipair)%unfold_index(i)
                iop=rmtx%irr_pair(ipair)%unfold_operation(i)
                p1=p%species( rmtx%full_pair(jpair)%a1 )
                p2=p%species( rmtx%full_pair(jpair)%a2 )

                ! Collect the full element
                obuf(p1,p2)%m=0.0_r8
                do b2=1,rmtx%full_pair(jpair)%nb2
                do b1=1,rmtx%full_pair(jpair)%nb1
                    j=rmtx%full_pair(jpair)%idx(b1,b2)
                    if ( j .le. 0) cycle
                    obuf(p1,p2)%m(b1,b2)=overlap_matrix(j)
                enddo
                enddo
                ! Rotate it to the right place
                if ( iop .lt. 0 ) then
                    call rmtx%full_pair(jpair)%rotate( basis,p,sym%op(-iop),&
                        forward=.false.,transposition=.true.,&
                        original_block=obuf(p1,p2)%m,&
                        rotated_block=rbuf(s1,s2)%m )
                else
                    call rmtx%full_pair(jpair)%rotate( basis,p,sym%op(iop),&
                        forward=.false.,transposition=.false.,&
                        original_block=obuf(p1,p2)%m,&
                        rotated_block=rbuf(s1,s2)%m )
                endif
                rmtx%irr_pair(ipair)%overlap=rmtx%irr_pair(ipair)%overlap+rbuf(s1,s2)%m*pref
            enddo
        enddo
    end select

!     ! Here is a dummy thingy to try and figure out how symmetry works, and if I got it right.
!     ! I'm a little skeptical for some reason.
!     do ipair=1,rmtx%n_full_pair
!         jpair=rmtx%full_pair(ipair)%index_unique
!         v0=rmtx%full_pair(ipair)%v
!         v1=rmtx%irr_pair(jpair)%v
!
!         s1=p%species( rmtx%full_pair(ipair)%a1 )
!         s2=p%species( rmtx%full_pair(ipair)%a2 )
!
!         p1=p%species( rmtx%irr_pair(jpair)%a1 )
!         p2=p%species( rmtx%irr_pair(jpair)%a2 )
!
!         obuf(p1,p2)%m=0.0_r8
!         rbuf(s1,s2)%m=0.0_r8
!         abuf(s1,s2)%m=0.0_r8
!
!         do b2=1,rmtx%irr_pair(jpair)%nb2
!         do b1=1,rmtx%irr_pair(jpair)%nb1
!             j=rmtx%irr_pair(jpair)%idx(b1,b2)
!             if ( j .le. 0) cycle
!             obuf(p1,p2)%m(b1,b2)=overlap_matrix(j)
!         enddo
!         enddo
!
! if ( norm2(obuf(p1,p2)%m) .lt. 1E-7_r8 ) cycle
!
!         ! Average?
!         ctr=0
!         do iop=1,sym%n_operation
!             if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)-v0 ) .lt. rl_sqtol ) then
!                 ctr=ctr+1
!                 call rmtx%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
!                     forward=.true.,transposition=.false.,&
!                     original_block=obuf(p1,p2)%m,&
!                     rotated_block=rbuf(s1,s2)%m )
!                 abuf(s1,s2)%m=abuf(s1,s2)%m+rbuf(s1,s2)%m
!             endif
!
!         enddo
!         ! do iop=1,sym%n_operation
!         !     if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)+v0 ) .lt. rl_sqtol ) then
!         !         ctr=ctr+1
!         !         call rmtx%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
!         !             forward=.true.,transposition=.true.,&
!         !             original_block=obuf(p1,p2)%m,&
!         !             rotated_block=rbuf(s1,s2)%m )
!         !         abuf(s1,s2)%m=abuf(s1,s2)%m+rbuf(s1,s2)%m
!         !     endif
!         ! enddo
!         abuf(s1,s2)%m=abuf(s1,s2)%m/real(ctr,r8)
!
!         ! Do it again and compare
!         ! Average?
! write(*,*) 'pair',ipair,jpair
!         ctr=0
!         do iop=1,sym%n_operation
!             if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)-v0 ) .lt. rl_sqtol ) then
!                 ctr=ctr+1
!                 call rmtx%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
!                     forward=.true.,transposition=.false.,&
!                     original_block=obuf(p1,p2)%m,&
!                     rotated_block=rbuf(s1,s2)%m )
!                 write(*,*) iop,norm2(abuf(s1,s2)%m-rbuf(s1,s2)%m),norm2(abuf(s1,s2)%m),norm2(rbuf(s1,s2)%m)
!             endif
!         enddo
!         do iop=1,sym%n_operation
!             if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)+v0 ) .lt. rl_sqtol ) then
!                 ctr=ctr+1
!                 call rmtx%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
!                     forward=.true.,transposition=.true.,&
!                     original_block=obuf(p1,p2)%m,&
!                     rotated_block=rbuf(s1,s2)%m )
!                 write(*,*) -iop,norm2(abuf(s1,s2)%m-rbuf(s1,s2)%m),norm2(abuf(s1,s2)%m),norm2(rbuf(s1,s2)%m)
!             endif
!         enddo
!
!         ! And the opposite operation maybe?
!         do iop=1,sym%n_operation
!             if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)-v0 ) .lt. rl_sqtol ) then
!                 ctr=ctr+1
!                 call rmtx%full_pair(ipair)%rotate( basis,p,sym%op(iop),&
!                     forward=.false.,transposition=.false.,&
!                     original_block=abuf(s1,s2)%m,&
!                     rotated_block=rbuf(s1,s2)%m )
!                 write(*,*) iop,norm2(obuf(s1,s2)%m-rbuf(s1,s2)%m),norm2(obuf(s1,s2)%m),norm2(rbuf(s1,s2)%m)
!             endif
!         enddo
!         do iop=1,sym%n_operation
!             if ( rl_sqnorm( matmul(sym%op(iop)%m,v1)+v0 ) .lt. rl_sqtol ) then
!                 ctr=ctr+1
!                 call rmtx%full_pair(ipair)%rotate( basis,p,sym%op(iop),&
!                     forward=.false.,transposition=.true.,&
!                     original_block=abuf(p1,p2)%m,&
!                     rotated_block=rbuf(s1,s2)%m )
!                 write(*,*) -iop,norm2(obuf(s1,s2)%m-rbuf(s1,s2)%m),norm2(obuf(s1,s2)%m),norm2(rbuf(s1,s2)%m)
!             endif
!         enddo
!
!
!     enddo
!
! write(*,*) 'stopping here'
! !call mw%destroy()
! stop

    ! Cleanup
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%deallocate(obuf(s1,s2)%m,persistent=.false.,scalable=.false.)
        call mem%deallocate(rbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
    enddo
    enddo
    deallocate(obuf)
    deallocate(rbuf)
end subroutine

!> collect the hamiltonian from AIMS. Should be retired.
subroutine collect_densitymatrix(rmtx,basis,densitymatrix)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> overlap matrix as it comes from AIMS
    real(r8), dimension(:), intent(in) :: densitymatrix

    integer :: i,j,b1,b2

    do i=1,rmtx%n_irr_pair
        do b1=1,rmtx%irr_pair(i)%nb1
        do b2=1,rmtx%irr_pair(i)%nb2
            j=rmtx%irr_pair(i)%idx(b1,b2)
            if ( j .le. 0) cycle
            rmtx%irr_pair(i)%densitymatrix(b1,b2,1)=densitymatrix(j)
        enddo
        enddo
    enddo
end subroutine

!> collect the hamiltonian from AIMS. Should be retired.
subroutine collect_hamiltonian(rmtx,basis,p,KS,sym,mw,mem,hamiltonian)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> kohnsham problem helper
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> hamiltonian as it comes from AIMS
    real(r8), dimension(:,:), intent(in) :: hamiltonian

    type(rl_mom_real), dimension(:,:), allocatable :: hbuf,rbuf
    real(r8) :: pref
    integer :: ipair,jpair,iop,ispin
    integer :: i,j,b1,b2,a1,a2,s1,s2,p1,p2

    ! Maybe make some space for rotation buffers
    allocate(hbuf(p%n_species,p%n_species))
    allocate(rbuf(p%n_species,p%n_species))
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%allocate(hbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        call mem%allocate(rbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        hbuf(s1,s2)%m=0.0_r8
        rbuf(s1,s2)%m=0.0_r8
    enddo
    enddo

    ! Set irreducible to nothing
    do i=1,rmtx%n_irr_pair
        rmtx%irr_pair(i)%hamiltonian=0.0_r8
    enddo

    ! Collect everything
    select case(sym%n_operation)
    case(1)
        do i=1,rmtx%n_irr_pair
            do b1=1,rmtx%irr_pair(i)%nb1
            do b2=1,rmtx%irr_pair(i)%nb2
                j=rmtx%irr_pair(i)%idx(b1,b2)
                if ( j .le. 0) cycle
                rmtx%irr_pair(i)%hamiltonian(b1,b2,:)=hamiltonian(j,:)
            enddo
            enddo
        enddo
    case default
        ! Start collecting
        do ipair=1,rmtx%n_irr_pair
            pref=1.0_r8/real(rmtx%irr_pair(ipair)%n_unfold_pair,r8)
            s1=p%species( rmtx%irr_pair(ipair)%a1 )
            s2=p%species( rmtx%irr_pair(ipair)%a2 )

            do i=1,rmtx%irr_pair(ipair)%n_unfold_pair
                jpair=rmtx%irr_pair(ipair)%unfold_index(i)
                iop=rmtx%irr_pair(ipair)%unfold_operation(i)
                p1=p%species( rmtx%full_pair(jpair)%a1 )
                p2=p%species( rmtx%full_pair(jpair)%a2 )

                do ispin=1,KS%n_spin
                    ! Collect the full element
                    hbuf(p1,p2)%m=0.0_r8
                    do b2=1,rmtx%full_pair(jpair)%nb2
                    do b1=1,rmtx%full_pair(jpair)%nb1
                        j=rmtx%full_pair(jpair)%idx(b1,b2)
                        if ( j .le. 0) cycle
                        hbuf(p1,p2)%m(b1,b2)=hamiltonian(j,ispin)
                    enddo
                    enddo
                    ! Rotate it to the right place
                    if ( iop .lt. 0 ) then
                        call rmtx%full_pair(jpair)%rotate( basis,p,sym%op(-iop),&
                            forward=.false.,transposition=.true.,&
                            original_block=hbuf(p1,p2)%m,&
                            rotated_block=rbuf(s1,s2)%m )
                    else
                        call rmtx%full_pair(jpair)%rotate( basis,p,sym%op(iop),&
                            forward=.false.,transposition=.false.,&
                            original_block=hbuf(p1,p2)%m,&
                            rotated_block=rbuf(s1,s2)%m )
                    endif
                    rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin)=rmtx%irr_pair(ipair)%hamiltonian(:,:,ispin)+rbuf(s1,s2)%m*pref
                enddo
            enddo
        enddo
    end select

    ! Cleanup
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%deallocate(hbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
        call mem%deallocate(rbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
    enddo
    enddo
    deallocate(hbuf)
    deallocate(rbuf)
end subroutine

!> inject Hamiltonian into AIMS. Should be retired.
subroutine inject_hamiltonian(rmtx,basis,p,KS,sym,mw,mem,hamiltonian)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> kohnsham problem helper
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> hamiltonian as it comes from AIMS
    real(r8), dimension(:,:), intent(inout) :: hamiltonian

    type(rl_mom_real), dimension(:,:), allocatable :: hbuf
    integer :: j,b1,b2,s1,s2,a1,a2
    integer :: ipair,jpair,iop,ispin

    ! Maybe make some space for rotation buffers
    allocate(hbuf(p%n_species,p%n_species))
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%allocate(hbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        hbuf(s1,s2)%m=0.0_r8
    enddo
    enddo

    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        ! This should be easy enough, maybe.
        do ipair=1,m%n_full_pair
            jpair=m%full_pair(ipair)%index_unique
            iop=m%full_pair(ipair)%index_operation
            a1=m%full_pair(ipair)%a1
            a2=m%full_pair(ipair)%a2
            s1=p%species(a1)
            s2=p%species(a2)
            do ispin=1,KS%n_spin
                select case(sym%n_operation)
                case(1)
                    ! Trivial operations
                    select case(iop)
                    case(-1)
                        hbuf(s1,s2)%m=transpose(m%irr_pair(jpair)%hamiltonian(:,:,ispin))
                    case(1)
                        hbuf(s1,s2)%m=m%irr_pair(jpair)%hamiltonian(:,:,ispin)
                    end select
                case default
                    ! Have to rotate perhaps
                    if ( iop .lt. 0 ) then
                        call m%irr_pair(jpair)%rotate( basis,p,sym%op(-iop),&
                            forward=.true.,transposition=.true.,&
                            original_block=m%irr_pair(jpair)%hamiltonian(:,:,ispin),&
                            rotated_block=hbuf(s1,s2)%m )
                    else
                        call m%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
                            forward=.true.,transposition=.false.,&
                            original_block=m%irr_pair(jpair)%hamiltonian(:,:,ispin),&
                            rotated_block=hbuf(s1,s2)%m )
                    endif
                end select
                ! Now we can place it
                do b2=1,m%full_pair(ipair)%nb2
                do b1=1,m%full_pair(ipair)%nb1
                    j=rmtx%full_pair(ipair)%idx(b1,b2)
                    if ( j .le. 0) cycle
                    hamiltonian(j,ispin)=hbuf(s1,s2)%m(b1,b2)
                enddo
                enddo
            enddo
        enddo
    type is(rl_realspace_matrix_mediumdistributed)
        write(*,*) 'FIXME INJECT'
        stop
    type is(rl_realspace_matrix_fulldistributed)
    end select

    ! Cleanup
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%deallocate(hbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
    enddo
    enddo
    deallocate(hbuf)

end subroutine

!> inject densitymatrix into AIMS. Should be retired.
subroutine compare_densitymatrix(rmtx,basis,p,KS,sym,mw,mem,density_matrix_sparse,ispin)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> kohnsham problem helper
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> density matrix from AIMS
    real(r8), dimension(:), intent(inout) :: density_matrix_sparse
    !> spin index from AIMS
    integer, intent(in) :: ispin

    type(rl_mom_real), dimension(:,:), allocatable :: dbuf,abuf
    integer :: j,b1,b2,s1,s2,a1,a2
    integer :: ipair,jpair,iop

    ! Maybe make some space for rotation buffers
    allocate(dbuf(p%n_species,p%n_species))
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%allocate(dbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        dbuf(s1,s2)%m=0.0_r8
    enddo
    enddo

    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        ! This should be easy enough, maybe.
        do ipair=1,m%n_full_pair
            jpair=m%full_pair(ipair)%index_unique
            iop=m%full_pair(ipair)%index_operation
            a1=m%full_pair(ipair)%a1
            a2=m%full_pair(ipair)%a2
            s1=p%species(a1)
            s2=p%species(a2)
            select case(sym%n_operation)
            case(1)
                ! Trivial operations
                select case(iop)
                case(-1)
                    dbuf(s1,s2)%m=transpose(m%irr_pair(jpair)%densitymatrix(:,:,ispin))
                case(1)
                    dbuf(s1,s2)%m=m%irr_pair(jpair)%densitymatrix(:,:,ispin)
                end select
            case default
                ! Have to rotate perhaps
                if ( iop .lt. 0 ) then
                    call m%irr_pair(jpair)%rotate( basis,p,sym%op(-iop),&
                        forward=.true.,transposition=.true.,&
                        original_block=m%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                        rotated_block=dbuf(s1,s2)%m )
                else
                    call m%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
                        forward=.true.,transposition=.false.,&
                        original_block=m%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                        rotated_block=dbuf(s1,s2)%m )
                endif
            end select
            ! Now we can place it
            do b2=1,m%full_pair(ipair)%nb2
            do b1=1,m%full_pair(ipair)%nb1
                j=rmtx%full_pair(ipair)%idx(b1,b2)
                if ( j .le. 0) cycle
                density_matrix_sparse(j)=dbuf(s1,s2)%m(b1,b2)
            enddo
            enddo
        enddo
    type is(rl_realspace_matrix_mediumdistributed)
        write(*,*) 'FIXME INJECT'
        stop
    type is(rl_realspace_matrix_fulldistributed)
    end select

    ! Cleanup
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%deallocate(dbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
    enddo
    enddo
    deallocate(dbuf)
end subroutine


!> inject densitymatrix into AIMS. Should be retired.
subroutine inject_densitymatrix(rmtx,basis,p,KS,sym,mw,mem,density_matrix_sparse,ispin)
    !> realspace hamiltonian/overlap/densitymatrix
    class(rl_realspace_matrix), intent(inout) :: rmtx
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> kohnsham problem helper
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> density matrix from AIMS
    real(r8), dimension(:), intent(inout) :: density_matrix_sparse
    !> spin index from AIMS
    integer, intent(in) :: ispin

    type(rl_mom_real), dimension(:,:), allocatable :: dbuf
    integer :: j,b1,b2,s1,s2,a1,a2
    integer :: ipair,jpair,iop

    ! Maybe make some space for rotation buffers
    allocate(dbuf(p%n_species,p%n_species))
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%allocate(dbuf(s1,s2)%m,[basis%species(s1)%n_basis,basis%species(s2)%n_basis],persistent=.false.,scalable=.false.)
        dbuf(s1,s2)%m=0.0_r8
    enddo
    enddo

    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        ! This should be easy enough, maybe.
        do ipair=1,m%n_full_pair
            jpair=m%full_pair(ipair)%index_unique
            iop=m%full_pair(ipair)%index_operation
            a1=m%full_pair(ipair)%a1
            a2=m%full_pair(ipair)%a2
            s1=p%species(a1)
            s2=p%species(a2)
            select case(sym%n_operation)
            case(1)
                ! Trivial operations
                select case(iop)
                case(-1)
                    dbuf(s1,s2)%m=transpose(m%irr_pair(jpair)%densitymatrix(:,:,ispin))
                case(1)
                    dbuf(s1,s2)%m=m%irr_pair(jpair)%densitymatrix(:,:,ispin)
                end select
            case default
                ! Have to rotate perhaps
                if ( iop .lt. 0 ) then
                    call m%irr_pair(jpair)%rotate( basis,p,sym%op(-iop),&
                        forward=.true.,transposition=.true.,&
                        original_block=m%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                        rotated_block=dbuf(s1,s2)%m )
                else
                    call m%irr_pair(jpair)%rotate( basis,p,sym%op(iop),&
                        forward=.true.,transposition=.false.,&
                        original_block=m%irr_pair(jpair)%densitymatrix(:,:,ispin),&
                        rotated_block=dbuf(s1,s2)%m )
                endif
            end select
            ! Now we can place it
            do b2=1,m%full_pair(ipair)%nb2
            do b1=1,m%full_pair(ipair)%nb1
                j=rmtx%full_pair(ipair)%idx(b1,b2)
                if ( j .le. 0) cycle
                density_matrix_sparse(j)=dbuf(s1,s2)%m(b1,b2)
            enddo
            enddo
        enddo
    type is(rl_realspace_matrix_mediumdistributed)
        write(*,*) 'FIXME INJECT'
        stop
    type is(rl_realspace_matrix_fulldistributed)
    end select

    ! Cleanup
    do s1=1,p%n_species
    do s2=1,p%n_species
        call mem%deallocate(dbuf(s1,s2)%m,persistent=.false.,scalable=.false.)
    enddo
    enddo
    deallocate(dbuf)
end subroutine

!> measure size in memory in bytes, approximately. Should be larger than the real usage, maybe.
function mtx_size_in_mem(rmtx) result(mem)
    !> extended cluster
    class(rl_realspace_matrix), intent(in) :: rmtx
    !> memory in bytes
    integer(i8) :: mem

    integer :: i

    mem=0
    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        mem=mem+storage_size(m)
    type is(rl_realspace_matrix_mediumdistributed)
        mem=mem+storage_size(m)
    type is(rl_realspace_matrix_fulldistributed)
        mem=mem+storage_size(m)
    end select

    ! These things always exit
    if ( allocated(rmtx%full_pair) ) then
        do i=1,size(rmtx%full_pair)
            mem=mem+storage_size(rmtx%full_pair(i))
            if ( allocated(rmtx%full_pair(i)%idx) ) mem=mem+storage_size(rmtx%full_pair(i)%idx)*size(rmtx%full_pair(i)%idx)
            if ( allocated(rmtx%full_pair(i)%ft_element) ) mem=mem+storage_size(rmtx%full_pair(i)%ft_element)*size(rmtx%full_pair(i)%ft_element)
        enddo
    endif
    if ( allocated(rmtx%irr_pair) ) then
        do i=1,size(rmtx%irr_pair)
            mem=mem+storage_size(rmtx%irr_pair(i))
            if ( allocated(rmtx%irr_pair(i)%unfold_index) ) mem=mem+storage_size(rmtx%irr_pair(i)%unfold_index)*size(rmtx%irr_pair(i)%unfold_index)
            if ( allocated(rmtx%irr_pair(i)%unfold_operation) ) mem=mem+storage_size(rmtx%irr_pair(i)%unfold_operation)*size(rmtx%irr_pair(i)%unfold_operation)
            if ( allocated(rmtx%irr_pair(i)%idx) ) mem=mem+storage_size(rmtx%irr_pair(i)%idx)*size(rmtx%irr_pair(i)%idx)
            if ( allocated(rmtx%irr_pair(i)%hamiltonian) ) mem=mem+storage_size(rmtx%irr_pair(i)%hamiltonian)*size(rmtx%irr_pair(i)%hamiltonian)
            if ( allocated(rmtx%irr_pair(i)%densitymatrix) ) mem=mem+storage_size(rmtx%irr_pair(i)%densitymatrix)*size(rmtx%irr_pair(i)%densitymatrix)
            if ( allocated(rmtx%irr_pair(i)%overlap) ) mem=mem+storage_size(rmtx%irr_pair(i)%overlap)*size(rmtx%irr_pair(i)%overlap)
            !if ( allocated(XXX) ) mem=mem+storage_size(XXX)*size(XXX)
        enddo
    endif
    mem=mem/8
end function

subroutine destroy_mtx(rmtx)
    !> extended cluster
    class(rl_realspace_matrix), intent(inout) :: rmtx

    integer :: i

    ! These things always exit
    if ( allocated(rmtx%full_pair) ) then
        do i=1,size(rmtx%full_pair)
            if ( allocated(rmtx%full_pair(i)%idx) )             deallocate(rmtx%full_pair(i)%idx)
            if ( allocated(rmtx%full_pair(i)%ft_element) )      deallocate(rmtx%full_pair(i)%ft_element)
        enddo
    endif
    if ( allocated(rmtx%irr_pair) ) then
        do i=1,size(rmtx%irr_pair)
            if ( allocated(rmtx%irr_pair(i)%unfold_index) )     deallocate(rmtx%irr_pair(i)%unfold_index)
            if ( allocated(rmtx%irr_pair(i)%unfold_operation) ) deallocate(rmtx%irr_pair(i)%unfold_operation)
            if ( allocated(rmtx%irr_pair(i)%idx) )              deallocate(rmtx%irr_pair(i)%idx)
            if ( allocated(rmtx%irr_pair(i)%hamiltonian) )      deallocate(rmtx%irr_pair(i)%hamiltonian)
            if ( allocated(rmtx%irr_pair(i)%densitymatrix) )    deallocate(rmtx%irr_pair(i)%densitymatrix)
            if ( allocated(rmtx%irr_pair(i)%overlap) )          deallocate(rmtx%irr_pair(i)%overlap)
        enddo
    endif
end subroutine

end module
