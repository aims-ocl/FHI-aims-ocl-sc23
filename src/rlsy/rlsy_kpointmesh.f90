module rlsy_kpointmesh
!!
!! This module handles meshes in reciprocal space
!!
use rlsy_memtracker, only: rl_memtracker
use rlsy_constants, only: r8,rl_huge,rl_hugeint,rl_exitcode_param,rl_iou,rl_exitcode_memory,rl_exitcode_symmetry,rl_sqtol
use rlsy_helpers, only: rl_chop,tochar,rl_clean_fractional_coordinates,rl_sqnorm
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list
use rlsy_sorting, only: rl_return_unique
use rlsy_verletlist, only: rl_verletbox
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup

implicit none
private
public :: rl_kpoint_mesh

!> general k-point
type rl_kpoint
    !> coordinate of k-point in Cartesian reciprocal coordinates
    real(r8), dimension(3) :: r=-rl_huge
    !> integration weight of this k-point. Should sum to 1.
    real(r8) :: integration_weight=-rl_huge

    ! Operations invariant at each k-point. Later problem.
    ! !> how many symmetry operations keep leave point invariant
    ! integer :: n_invariant_operation=-rl_hugeint
    ! !> which are the invariant operations
    ! integer, dimension(:), allocatable :: invariant_operation
!end type
!
!!> a k-point in the full mesh
!type, extends(rl_kpoint) :: rl_full_kpoint
    !> which irreducible k-point is it equivalent to
    integer :: irreducible_index=-rl_hugeint
    !> which operation transforms the irreducible to this point
    integer :: operation_from_irreducible=-rl_hugeint
!end type
!
!!> an irreducible k-point
!type, extends(rl_kpoint) :: rl_irr_kpoint
    !> which index in the full grid does it correspond to?
    integer :: full_index=-rl_hugeint
    !> how many points in the full mesh can this k-point transform to
    integer :: n_full_point=-rl_hugeint
    !> which full points can it transform to
    integer, dimension(:), allocatable :: index_full_point
    !> which symmetry operation does the transformation to the index specified above
    integer, dimension(:), allocatable :: operation_full_point
end type

!> everything about the k-point mesh that might be useful
type rl_kpoint_mesh
    !> how many k-points in total
    integer :: n_full_kpoint=-rl_hugeint
    !> how many irreducible k-points in total
    integer :: n_irr_kpoint=-rl_hugeint
    !> all k-points
    type(rl_kpoint), dimension(:), allocatable :: ap
    !> irreducible k-points
    type(rl_kpoint), dimension(:), allocatable :: ip
    !> index array that can map from the AIMS k-points to my full k-points, and vice versa
    integer, dimension(:), allocatable :: ind_aims_to_rl
    integer, dimension(:), allocatable :: ind_rl_to_aims
    contains
        !> create the kspace handle
        procedure :: generate
        !> measure size in memory, approximately
        !procedure :: size_in_mem
end type

contains

!> creates the k-point mesh. Uses the data from AIMS for now.
subroutine generate(kmesh,p,sym,mw,mem,verbosity,&
    n_k_points_xyz,k_point_list,k_points_offset)
    !> k-space mesh
    class(rl_kpoint_mesh), intent(out) :: kmesh
    !> structure
    type(rl_crystalstructure), intent(inout) :: p
    !> space group
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> data from AIMS
    integer, dimension(3), intent(in) :: n_k_points_xyz
    real(r8), dimension(:,:), intent(in) :: k_point_list
    real(r8), dimension(3), intent(in) :: k_points_offset

    type(rl_verletbox) :: vb
    real(r8) :: timer,t0,t1
    real(r8), dimension(:,:), allocatable :: kv

    !init: block
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(3) :: v0
        integer :: i,l
    !reduce: block
        !real(r8), dimension(:,:), allocatable :: dr,dn
        real(r8), dimension(:,:), allocatable :: dn
        integer, dimension(:,:), allocatable :: perm,shell_member
        integer, dimension(:), allocatable :: shell_ctr,shell_index,prototype_shell,di,dj
        !integer :: iop,ct1,ct2,ctr_op,i,k
        integer :: iop,ct1,ct2,ctr_op,k
    !fillinformation: block
        !real(r8), dimension(3) :: v0,v1,v2
        real(r8), dimension(3) :: v1,v2
        !integer :: i,j,k,l
        integer :: j
    !mapaimstomine: block
        real(r8), dimension(:,:), allocatable :: dr0,dr1
        !integer :: i,j,k,l

    !init: block

        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
            t0=timer
            t1=timer
            write(rl_iou,*) ''
            write(rl_iou,*) 'REDUCING KPOINTS WITH SYMMETRY'
        endif

        ! So, I won't bother with distributing the k-mesh or anything like that,
        ! the size of the meshes in DFT are negligible. So I will start by just
        ! grabbing the list of k-vectors. But I want all of them, not just a few
        ! at random, which is what you get from AIMS after their pseudo-reduction.
        l=size(k_point_list,1)
        call mem%allocate(dr,[3,2*l],persistent=.false.,scalable=.false.)
        dr=0.0_r8
        do i=1,l
            dr(:,i)=k_point_list(i,:)
            dr(:,i+l)=-dr(:,i)
        enddo
        do i=1,size(dr,2)
            dr(:,i)=rl_chop(dr(:,i),1E-10_r8)
            dr(:,i)=rl_clean_fractional_coordinates(dr(:,i))
            dr(:,i)=rl_chop(dr(:,i),1E-10_r8)
            dr(:,i)=rl_clean_fractional_coordinates(dr(:,i))
        enddo
        call rl_return_unique(dr,kv)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.)

        ! Check that this makes sense
        if ( size(kv,2) .ne. product(n_k_points_xyz) ) then
            call rl_stop_gracefully(['Clearly I do not understand symmetry.'],rl_exitcode_symmetry,mw%comm)
        endif

        ! Good, now I know the total number of k-points.
        kmesh%n_full_kpoint=size(kv,2)

        ! Think I want Gamma as the first point, in case it is included in the mesh.
        l=0
        do i=1,size(kv,2)
            if ( rl_sqnorm(kv(:,i)) .lt. rl_sqtol ) then
                l=i
                exit
            endif
        enddo
        if ( l .gt. 0 ) then
            kv(:,l)=kv(:,1)
            kv(:,1)=0.0_r8
        endif

        ! Make a Verlet box for locating points
        call vb%generate(kv,[15,15,15],mem)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... expanded number of k-points to ',tochar(size(kv,2)),' from ',tochar(size(k_point_list,1))
            t0=timer
        endif
    !end block init

    !reduce: block

        ! Create the permutation list
        call mem%allocate(perm,[kmesh%n_full_kpoint,2*sym%n_operation],persistent=.false.,scalable=.false.,supress_error=.true.)
        call mem%allocate(di,kmesh%n_full_kpoint,persistent=.false.,scalable=.false.)
        call mem%allocate(dj,kmesh%n_full_kpoint,persistent=.false.,scalable=.false.)
        call mem%allocate(dr,[3,kmesh%n_full_kpoint],persistent=.false.,scalable=.false.)
        call mem%allocate(dn,[3,kmesh%n_full_kpoint],persistent=.false.,scalable=.false.)
        perm=0
        di=0
        dj=0
        dr=0.0_r8
        dn=0.0_r8
        ctr_op=0
        do iop=1,sym%n_operation
            ! Reset counters
            ct1=0
            ct2=0
            di=0
            dj=0
            ! Rotate the points
            call dgemm('N','N',3,kmesh%n_full_kpoint,3,1.0_r8,sym%op(iop)%rfm,3,kv,3,0.0_r8,dr,3)
            !dr=rl_clean_fractional_coordinates(rl_chop(dr,1E-10_r8))
            dr=rl_chop(dr,1E-10_r8)
            dr=rl_clean_fractional_coordinates(dr(:,:))
            !dr=rl_clean_fractional_coordinates(rl_chop(dr,1E-10_r8))
            dr=rl_chop(dr,1E-10_r8)
            dr=rl_clean_fractional_coordinates(dr(:,:))
            ! Then we have the same operation but with time-reversal added on.
            dn=-dr
            !dn=rl_clean_fractional_coordinates(rl_chop(dn,1E-10_r8))
            dn=rl_chop(dr,1E-10_r8)
            dn=rl_clean_fractional_coordinates(dn(:,:))
            !dn=rl_clean_fractional_coordinates(rl_chop(dn,1E-10_r8))
            dn=rl_chop(dr,1E-10_r8)
            dn=rl_clean_fractional_coordinates(dn(:,:))
            ! Use the Verlet boxes to check where each point ended up.
            do i=1,kmesh%n_full_kpoint
                ! locate index of the transformed point in the original array
                k=vb%locate(kv,dr(:,i))
                if ( k .gt. 0 ) then
                    ct1=ct1+1
                    di(i)=k
                else
                    ! Make note that it's not a good operation
                    ct1=0
                endif

                k=vb%locate(kv,dn(:,i))
                if ( k .gt. 0 ) then
                    ct2=ct2+1
                    dj(i)=k
                else
                    ! Make note that it's not a good operation
                    ct1=0
                endif
            enddo
            ! Check if the operation was ok, and in that case store it.
            if ( ct1 .eq. kmesh%n_full_kpoint ) then
                ctr_op=ctr_op+1
                perm(:,ctr_op)=di
            endif
            if ( ct2 .eq. kmesh%n_full_kpoint ) then
                ctr_op=ctr_op+1
                perm(:,ctr_op)=dj
            endif
        enddo

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... built permutation list (',tochar(t1-t0),')'
            t0=timer
        endif

        ! With the permutation list I can slice the points into shells and stuff.
        call rl_coordination_shells_from_permutation_list(perm(:,1:ctr_op),shell_ctr,shell_member,shell_index,mem,prototype_shell,mw=mw)

        ! Now start storing information.
        kmesh%n_irr_kpoint=size(prototype_shell)
        allocate(kmesh%ap(kmesh%n_full_kpoint))
        allocate(kmesh%ip(kmesh%n_irr_kpoint))
        do i=1,kmesh%n_full_kpoint
            kmesh%ap(i)%r=kv(:,i)
            kmesh%ap(i)%integration_weight=0.0_r8
            kmesh%ap(i)%irreducible_index=shell_index(i)
            kmesh%ap(i)%operation_from_irreducible=0
        enddo
        do i=1,kmesh%n_irr_kpoint
            kmesh%ip(i)%full_index=prototype_shell(i)
            kmesh%ip(i)%r=kv(:,kmesh%ip(i)%full_index)
            kmesh%ip(i)%integration_weight=0.0_r8
            kmesh%ip(i)%n_full_point=shell_ctr(i)
            allocate(kmesh%ip(i)%index_full_point( kmesh%ip(i)%n_full_point ))
            allocate(kmesh%ip(i)%operation_full_point( kmesh%ip(i)%n_full_point ))
            kmesh%ip(i)%index_full_point=shell_member(1:shell_ctr(i),i)
            kmesh%ip(i)%operation_full_point=0
        enddo

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... found ',tochar(kmesh%n_irr_kpoint),' irreducible points out of ',tochar(kmesh%n_full_kpoint),' (',tochar(t1-t0),')'
            t0=timer
        endif

        ! Cleanup
        call mem%deallocate(perm,persistent=.false.,scalable=.false.)
        call mem%deallocate(di,persistent=.false.,scalable=.false.)
        call mem%deallocate(dj,persistent=.false.,scalable=.false.)
        call mem%deallocate(dr,persistent=.false.,scalable=.false.)
        call mem%deallocate(dn,persistent=.false.,scalable=.false.)
        call mem%deallocate(shell_ctr,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_index,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_member,persistent=.true.,scalable=.false.)
        call mem%deallocate(prototype_shell,persistent=.true.,scalable=.false.)
    !end block reduce

    ! Now figure out the rest of the stuff, how things are related and all that.
    !fillinformation: block

        ! Get all the operations sorted out
        ! First how to build full from the irreducible
        do i=1,kmesh%n_full_kpoint
            v0=kmesh%ap(i)%r
            v1=kmesh%ip( kmesh%ap(i)%irreducible_index )%r
            ! First try normal operations
            l=0
            do j=1,sym%n_operation
                v2=matmul(sym%op(j)%rfm,v1)
                v2=rl_clean_fractional_coordinates(rl_chop(rl_clean_fractional_coordinates(v2),1E-10_r8))
                if ( rl_sqnorm(v2-v0) .lt. rl_sqtol ) then
                    l=j
                    exit
                endif
            enddo
            ! Then with time-reversal
            if ( l .eq. 0 ) then
                do j=1,sym%n_operation
                    v2=-matmul(sym%op(j)%rfm,v1)
                    v2=rl_clean_fractional_coordinates(rl_chop(rl_clean_fractional_coordinates(v2),1E-10_r8))
                    if ( rl_sqnorm(v2-v0) .lt. rl_sqtol ) then
                        l=-j
                        exit
                    endif
                enddo
            endif
            if ( l .ne. 0 ) then
                kmesh%ap(i)%operation_from_irreducible=l
            else
                call rl_stop_gracefully(['Clearly I do not understand symmetry.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo

        ! Then almost the same thing again. Could probably reuse the thing above
        ! but I don't have the energy. This takes no time anyway.
        do i=1,kmesh%n_irr_kpoint
            v0=kmesh%ip(i)%r
            do k=1,kmesh%ip(i)%n_full_point
                v1=kmesh%ap( kmesh%ip(i)%index_full_point(k) )%r
                l=0
                do j=1,sym%n_operation
                    v2=matmul(sym%op(j)%rfm,v0)
                    v2=rl_clean_fractional_coordinates(rl_chop(rl_clean_fractional_coordinates(v2),1E-10_r8))
                    if ( rl_sqnorm(v2-v1) .lt. rl_sqtol ) then
                        l=j
                        exit
                    endif
                enddo
                ! Then with time-reversal
                if ( l .eq. 0 ) then
                    do j=1,sym%n_operation
                        v2=-matmul(sym%op(j)%rfm,v0)
                        v2=rl_clean_fractional_coordinates(rl_chop(rl_clean_fractional_coordinates(v2),1E-10_r8))
                        if ( rl_sqnorm(v2-v1) .lt. rl_sqtol ) then
                            l=-j
                            exit
                        endif
                    enddo
                endif
                if ( l .ne. 0 ) then
                    kmesh%ip(i)%operation_full_point(k)=l
                else
                    call rl_stop_gracefully(['Clearly I do not understand symmetry.'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo
        enddo

        ! Set integration weights. Not they sum to 1, so they include the volume factor thingy.
        do i=1,kmesh%n_full_kpoint
            kmesh%ap(i)%integration_weight=1.0_r8/real(kmesh%n_full_kpoint,r8)
            j=kmesh%ap(i)%irreducible_index
            kmesh%ip(j)%integration_weight=kmesh%ip(j)%integration_weight+1.0_r8
        enddo
        do i=1,kmesh%n_irr_kpoint
            kmesh%ip(i)%integration_weight=kmesh%ip(i)%integration_weight/real(kmesh%n_full_kpoint,r8)
        enddo

        ! And convert coordinates to Cartesian
        do i=1,kmesh%n_full_kpoint
            kmesh%ap(i)%r=rl_chop(matmul(p%reciprocal_latticevectors,kmesh%ap(i)%r),1E-13_r8)
        enddo
        do i=1,kmesh%n_irr_kpoint
            kmesh%ip(i)%r=rl_chop(matmul(p%reciprocal_latticevectors,kmesh%ip(i)%r),1E-13_r8)
        enddo

        !@TODO Insert some agressive sanity tests here.
    !end block fillinformation

    ! A little cleanup
    deallocate(kv)
    call vb%destroy(mem)

    ! Get an index array that can convert from an AIMS kpoint index to one of my k-points
    !mapaimstomine: block

        ! Get two lists of k-points, dr0 holds the list from AIMS, dr1 holds mine.
        call mem%allocate(dr0,[size(k_point_list,2),size(k_point_list,1)],persistent=.false.,scalable=.false.)
        call mem%allocate(dr1,[3,kmesh%n_full_kpoint],persistent=.false.,scalable=.false.)
        dr0=transpose(k_point_list)
        dr1=0.0_r8
        do i=1,kmesh%n_full_kpoint
            dr1(:,i)=matmul(p%inv_reciprocal_latticevectors,kmesh%ap(i)%r)
        enddo
        ! Make sure we really have clean coordinates. Should not be necessary, but does not really hurt.
        dr0=rl_clean_fractional_coordinates(dr0(:,:))
        dr0=rl_chop(dr0,1E-10_r8)
        dr0=rl_clean_fractional_coordinates(dr0(:,:))
        dr0=rl_chop(dr0,1E-10_r8)
        dr1=rl_clean_fractional_coordinates(dr1(:,:))
        dr1=rl_chop(dr1,1E-10_r8)
        dr1=rl_clean_fractional_coordinates(dr1(:,:))
        dr1=rl_chop(dr1,1E-10_r8)
        ! Create a Verlet box for fast lookup
        call vb%generate(dr1,[13,13,13],mem)

        ! Space for the index map
        allocate(kmesh%ind_aims_to_rl(size(dr0,2)))
        kmesh%ind_aims_to_rl=0
        do i=1,size(dr0,2)
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            j=vb%locate(dr1,dr0(:,i))
            if ( j .gt. 0 ) then
                kmesh%ind_aims_to_rl(i)=j
            else
                call rl_stop_gracefully(['Could not map from AIMS k-points to mine. Should never happen.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo
        call mw%allreduce('sum',kmesh%ind_aims_to_rl)

        ! Then we might as well get the reverse map as well.
        call vb%destroy(mem)
        call vb%generate(dr0,[13,13,13],mem)
        allocate(kmesh%ind_rl_to_aims(kmesh%n_full_kpoint))
        kmesh%ind_rl_to_aims=0
        do i=1,kmesh%n_full_kpoint
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            j=vb%locate(dr0,dr1(:,i))
            if ( j .gt. 0 ) then
                kmesh%ind_rl_to_aims(i)=j
            endif
        enddo
        call mw%allreduce('sum',kmesh%ind_rl_to_aims)
        call vb%destroy(mem)

        ! And a little cleanup
        call mem%deallocate(dr0,persistent=.false.,scalable=.false.)
        call mem%deallocate(dr1,persistent=.false.,scalable=.false.)

    !end block mapaimstomine
    ! Check that I did not do anything stupid
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)

    if ( verbosity .gt. 0 ) then
        t1=mpi_wtime()
        write(rl_iou,*) '... generated k-point mesh (',tochar(mpi_wtime()-timer),')'
        t0=timer
    endif

end subroutine

end module
