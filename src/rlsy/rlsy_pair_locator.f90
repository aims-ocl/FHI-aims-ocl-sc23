module rlsy_pair_locator
!!
!! Somewhat useful helper object: it is a quite slow way of taking two indices
!! in the extended cluster and returning an index to a pair, full or irreducible.
!! Why this is so annoyingly complex is that the number of pairs in the extended
!! cluster can quickly become very large, so one can not just go for the naive
!! implementation of storing a huge array.
!!
!! What we have is instead a list, per extended atom, of which pairs it has, sort
!! of. By selecting the first index in the extended cluster, ie, I do a binary search
!! matching it to those je that are associated with a pair. Fast enough to not be a
!! bottleneck, I think.
!!
!! Sould maybe move this to the extended_cluster module, but that one is already
!! quite long.
!!
use rlsy_constants, only: r8,rl_hugeint,rl_exitcode_param
use rlsy_memtracker, only: rl_memtracker
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully
use rlsy_sorting, only: rl_qsort
use rlsy_extended_cluster, only: rl_extended_cluster
use rlsy_realspace_matrix, only: rl_realspace_matrix,rl_realspace_matrix_notdistributed,&
    rl_realspace_matrix_mediumdistributed,rl_realspace_matrix_fulldistributed
use rlsy_kspace_eigenproblem, only: rl_kspace_eigenproblem,&
    rl_kspace_eigenproblem_singleproc,rl_kspace_eigenproblem_multiproc
use rlsy_verletlist, only: rl_verletbox

implicit none
private

public :: rl_pair_locator

!> list of pairs for an atom
type rl_pair_locator_atom
    !> how many pairs does this atom have?
    integer :: npair=-rl_hugeint
    !> indices to extended atoms
    integer, dimension(:), allocatable :: index_ext
    !> indices to extended pair
    integer, dimension(:), allocatable :: index_pair
end type
!> helper that locates pairs from the extended cluster
type rl_pair_locator
    !> list of pairs for each atom
    type(rl_pair_locator_atom), dimension(:), allocatable :: atom
    !> number of pairs that can be located
    integer :: n_pair=-rl_hugeint
    contains
        !> build pair locator
        procedure :: generate=>create_pair_locator
        !> locate the index of a pair
        procedure :: locate=>find_pair_index
        !> destroy and deallocate
        !procedure :: destroy
end type

contains

!> find the pair index of a generic pair in the supercluster
function find_pair_index(pl,ie,je) result(ipair)
    !> pair locator
    class(rl_pair_locator), intent(in) :: pl
    !> first atom in extended cluster
    integer, intent(in) :: ie
    !> second atom in extended cluster
    integer, intent(in) :: je
    !> index to pair
    integer :: ipair

    integer :: i
    integer :: lo,hi,mid

    ! dumbest possible way
    ipair=-1
    ! Could be a really quick exit
    if ( pl%atom(ie)%npair .eq. 0 ) then
        return
    endif

    ! See if we are below the range, or touching
    i=pl%atom(ie)%index_ext(1)
    if ( je .lt. i ) then
        return
    elseif ( je .eq. i ) then
        ipair=pl%atom(ie)%index_pair(1)
        return
    endif

    ! See if we are above the range, or touching
    i=pl%atom(ie)%index_ext( pl%atom(ie)%npair )
    if ( je .gt. i ) then
        return
    elseif ( je .eq. i ) then
        ipair=pl%atom(ie)%index_pair( pl%atom(ie)%npair )
        return
    endif

    ! Stupid linear search that for some reason is quite fast
    ! do i=1,pl%atom(ie)%npair
    !     if ( pl%atom(ie)%index_ext(i) .eq. je ) then
    !         ipair=pl%atom(ie)%index_pair(i)
    !         return
    !     endif
    ! enddo

    ! Slightly smarter binary search that is slow and stupid. Hmmm.
    lo=1
    hi=size(pl%atom(ie)%index_ext)
    do while(lo<=hi)
        mid=(hi+lo)/2
        if ( pl%atom(ie)%index_ext(mid) .eq. je ) then
            ipair=pl%atom(ie)%index_pair(mid)
            return
        elseif ( pl%atom(ie)%index_ext(mid) .gt. je ) then
            hi=mid-1
        else
            lo=mid+1
        endif
    enddo
end function

!> create a pair locator function thingy.
subroutine create_pair_locator(pl,rmtx,KS,ec,irreducible,mw,mem)
    !> pair locator
    class(rl_pair_locator), intent(out) :: pl
    !> distributed realspace matrices
    class(rl_realspace_matrix), intent(in) :: rmtx
    !> Kohn-Sham solution thingy handle.
    class(rl_kspace_eigenproblem), intent(in) :: KS
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> should we locate full or irreducible pairs?
    logical, intent(in) :: irreducible
    !> MPI helper
    type(rl_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    type(rl_verletbox) :: vb
    integer, parameter :: boxdim=9
    real(r8), dimension(:,:), allocatable :: dv
    real(r8), dimension(3) :: v0,v1
    integer, dimension(:,:), allocatable :: di,dj
    integer, dimension(:), allocatable :: ctr,ind,atom1,atom2
    integer :: ie,je,ipair,i,a1,a2
    integer :: npair

    ! The point here is that I want to be able to take an arbitrary pair of atoms
    ! from the extended cluster, and say which of the pairs in rmtx that correspond
    ! to. What makes it a little complicated is that the number of pairs in the
    ! extended cluster can be very many, so care has to be taken such that this
    ! does not eat up all memory available.

    ! Put extended cluster into boxes
    !@TODO Be a little smarter when selecting box sizes here.
    !@TODO Make the whole thing a little memory-distributed for good measure.
    call vb%generate(ec%cartesian_coordinate,[boxdim,boxdim,boxdim],mem)

    ! Make space for counter
    call mem%allocate(ctr,ec%n_extended_atom,persistent=.false.,scalable=.false.)
    ctr=0

    ! First determine how many pairs we are going to work with:
    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        if ( irreducible ) then
            npair=m%n_irr_pair
        else
            npair=m%n_full_pair
        endif
    type is(rl_realspace_matrix_mediumdistributed)
        if ( irreducible ) then
            npair=m%n_irr_pair_global
        else
            npair=m%n_full_pair_global
        endif
    type is(rl_realspace_matrix_fulldistributed)
        if ( irreducible ) then
            npair=m%n_irr_pair_global
        else
            npair=m%n_full_pair_global
        endif
    end select

    ! Get some dummy space
    call mem%allocate(dv,[3,npair],persistent=.false.,scalable=.false.)
    call mem%allocate(atom1,npair,persistent=.false.,scalable=.false.)
    call mem%allocate(atom2,npair,persistent=.false.,scalable=.false.)
    dv=0.0_r8
    atom1=0
    atom2=0

    ! Fetch all the pairs
    select type(m=>rmtx)
    type is(rl_realspace_matrix_notdistributed)
        if ( irreducible ) then
            do ipair=1,npair
                dv(:,ipair)=m%irr_pair(ipair)%v
                atom1(ipair)=m%irr_pair(ipair)%a1
                atom2(ipair)=m%irr_pair(ipair)%a2
            enddo
        else
            do ipair=1,npair
                dv(:,ipair)=m%full_pair(ipair)%v
                atom1(ipair)=m%full_pair(ipair)%a1
                atom2(ipair)=m%full_pair(ipair)%a2
            enddo
        endif
    type is(rl_realspace_matrix_mediumdistributed)
        if ( irreducible ) then
            do i=1,m%n_irr_pair
                ipair=m%irr_offset+i
                dv(:,ipair)=m%irr_pair(i)%v
                atom1(ipair)=m%irr_pair(i)%a1
                atom2(ipair)=m%irr_pair(i)%a2
            enddo
        else
            do i=1,m%n_full_pair
                ipair=m%full_offset+i
                dv(:,ipair)=m%full_pair(i)%v
                atom1(ipair)=m%full_pair(i)%a1
                atom2(ipair)=m%full_pair(i)%a2
            enddo
        endif
        select type(KS)
        type is(rl_kspace_eigenproblem_multiproc)
            call KS%ml%allreduce('sum',dv)
            call KS%ml%allreduce('sum',atom1)
            call KS%ml%allreduce('sum',atom2)
        class default
            call rl_stop_gracefully(['Single-proc + distributed should not be able to happen.'],rl_exitcode_param,mw%comm)
        end select
    type is(rl_realspace_matrix_fulldistributed)
        call rl_stop_gracefully(['Not done with all distributions in pair locator.'],rl_exitcode_param,mw%comm)
    end select

    ! Count the number of pairs per extended atom
    do ie=1,ec%n_extended_atom
        if ( mod(ie,mw%n) .ne. mw%r ) cycle
        v0=ec%cartesian_coordinate(:,ie)
        a1=ec%index_unit_cell(ie)
        do ipair=1,npair
            if ( a1 .ne. atom1(ipair) ) cycle
            v1=v0+dv(:,ipair)
            je=vb%locate(ec%cartesian_coordinate,v1)
            if ( je .gt. 0 ) then
                a2=ec%index_unit_cell(ie)
                ! Think it should be impossible that a2 atom atom2(ipair) are not
                ! equal, might as well test for that.
                if ( a2 .ne. atom2(ipair) ) then
                    call rl_stop_gracefully(['Error matching pair, inconsistent species in pair locator.'],rl_exitcode_param,mw%comm)
                else
                    ctr(ie)=ctr(ie)+1
                endif
            endif
        enddo
    enddo
    call mw%allreduce('sum',ctr)

    ! Make some temporary space
    i=maxval(ctr)
    call mem%allocate(di,[i,ec%n_extended_atom],persistent=.false.,scalable=.false.)
    call mem%allocate(dj,[i,ec%n_extended_atom],persistent=.false.,scalable=.false.)
    call mem%allocate(ind,i,persistent=.false.,scalable=.false.)
    di=0
    dj=0
    ind=0
    ctr=0
    ! Store pair index things
    do ie=1,ec%n_extended_atom
        if ( mod(ie,mw%n) .ne. mw%r ) cycle
        a1=ec%index_unit_cell(ie)
        v0=ec%cartesian_coordinate(:,ie)
        do ipair=1,npair
            ! need the same starting unit cell atom
            if ( a1 .ne. atom1(ipair) ) cycle
            v1=v0+dv(:,ipair)
            je=vb%locate(ec%cartesian_coordinate,v1)
            if ( je .gt. 0 ) then
                ctr(ie)=ctr(ie)+1
                di(ctr(ie),ie)=ipair
                dj(ctr(ie),ie)=je
            endif
        enddo

        if ( ctr(ie) .gt. 0 ) then
            ! Sort the list according to je to make lookup faster later.
            call rl_qsort(dj(1:ctr(ie),ie),ind(1:ctr(ie)))
            di(1:ctr(ie),ie)=di(ind(1:ctr(ie)),ie)
        endif
    enddo
    call mw%allreduce('sum',ctr)
    call mw%allreduce('sum',di)
    call mw%allreduce('sum',dj)

    ! Start storing things
    allocate(pl%atom(ec%n_extended_atom))
    do ie=1,ec%n_extended_atom
        pl%atom(ie)%npair=ctr(ie)
        if ( pl%atom(ie)%npair .gt. 0 ) then
            ! Store indices to je and ipair
            allocate(pl%atom(ie)%index_ext( pl%atom(ie)%npair ))
            allocate(pl%atom(ie)%index_pair( pl%atom(ie)%npair ))
            pl%atom(ie)%index_ext=dj(1:ctr(ie),ie)
            pl%atom(ie)%index_pair=di(1:ctr(ie),ie)
        else
            ! No pairs originating from this atom.
            allocate(pl%atom(ie)%index_ext( 1 ))
            allocate(pl%atom(ie)%index_pair( 1 ))
            pl%atom(ie)%index_ext=-rl_hugeint
            pl%atom(ie)%index_pair=-rl_hugeint
        endif
    enddo

    ! Also store the number of pairs, globally
    pl%n_pair=npair

    ! And a little cleanup
    call vb%destroy(mem)
    call mem%deallocate(di,persistent=.false.,scalable=.false.)
    call mem%deallocate(dj,persistent=.false.,scalable=.false.)
    call mem%deallocate(ind,persistent=.false.,scalable=.false.)
    call mem%deallocate(ctr,persistent=.false.,scalable=.false.)
    call mem%deallocate(dv,persistent=.false.,scalable=.false.)
    call mem%deallocate(atom1,persistent=.false.,scalable=.false.)
    call mem%deallocate(atom2,persistent=.false.,scalable=.false.)
end subroutine

end module
