module rlsy_extended_cluster
!!
!! This module handles the extended cluster built up by all the atoms that have any
!! overlap with the central cell, and ensures that the cluster obeys symmetry.
!!
!!               o   o   o   o
!!
!!           o   o   o   o   o   o
!!
!!       o   o   o   o   o   o   o   o
!!                 ┌───────┐
!!       o   o   o │ o   o │ o   o   o
!!                 │   x   │
!!       o   o   o │ o   o │ o   o   o
!!                 └───────┘
!!       o   o   o   o   o   o   o   o
!!
!!           o   o   o   o   o   o
!!
!!               o   o   o   o
!!
!! So, to make it clear. The box above is the unit cell. This routine determines
!! all the little rings, that is the atoms whose basis function can touch any part
!! of the box. The 'x' marks the origin in the coordinate system.
!!
!! It is done in two steps, first a super-Wigner-Seitz cell is built, and used to work out
!! the symmetry of the supercluster. After the symmetry is figured out, the cluster gets
!! pruned down to just the absolute necessary. Should end up with exactly the same supercluster
!! that AIMS always uses, but with some metadata attached.
!!
!! For non-periodic systems this routine would need minor modifications, but is also not
!! really necessary, since the cluster would just be a copy of the normal list of atoms.
!!
!! In addition, there is a helper type defined here that can return a pruned list of
!! how many atoms in the supercluster are within some distance of an arbitrary point in
!! the unit cell, useful for speeding up lookups and distance calculations.
!!
use rlsy_constants, only: r8,rl_iou,rl_huge,rl_hugeint,rl_tol,rl_sqtol,rl_exitcode_symmetry,&
    rl_exitcode_param,rl_exitcode_memory
use rlsy_memtracker, only: rl_memtracker
use rlsy_helpers, only: tochar,rl_sqnorm,rl_clean_fractional_coordinates, norm2
use rlsy_sorting, only: rl_qsort,rl_return_unique
use rlsy_geometry, only: rl_bounding_sphere_of_box,rl_inscribed_sphere_in_box
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime
use rlsy_voronoi, only: rl_voronoi_diagram
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_distancetable, only: rl_distancetable
use rlsy_verletlist, only: rl_verletbox
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list

implicit none
private
public :: rl_extended_cluster
public :: rl_extended_cluster_hashedspace

!> coordination shell of extended atoms
type rl_extended_cluster_shell
    !> how many atoms in this shell
    integer :: n_extended_atom_in_shell=-rl_hugeint
    !> which atoms belong to this shell
    integer, dimension(:), allocatable :: index_ext_atom_in_shell
end type

!> an extended cluster embedded in a periodic system
type rl_extended_cluster
    !> how many extended atoms are there
    integer :: n_extended_atom=-rl_hugeint
    !> how many shells of extended atoms are there
    integer :: n_extended_shell=-rl_hugeint
    !> cartesian coordinates of the extended atoms
    real(r8), dimension(:,:), allocatable :: cartesian_coordinate
    !> unit cell index of the extended atoms
    integer, dimension(:), allocatable :: index_unit_cell
    !> coordination shells of extended atoms
    type(rl_extended_cluster_shell), dimension(:), allocatable :: shell
    !> cluster is only invariant under the symmorphic space group, store that here as well.
    type(rl_spacegroup) :: pointgroup
    !> indices to good prototype unit cell atoms in the cluster
    integer, dimension(:), allocatable :: locate_unit_cell_atom_in_cluster

    ! The variables below are only available after pruning.
    ! Possible a good idea to move them somewhere else. Or not.
    ! Maybe add a status variable telling us wether we are pruned
    ! or not. Time will tell what the best thing is.

    ! !> longest cutoff for confinement potential per atom
    ! real(r8), dimension(:), allocatable :: confinement_cutoff_sq
    ! !> basis index offset per atom.
    ! integer, dimension(:), allocatable :: basis_index_offset

    contains
        !> create the cluster
        procedure :: generate
        !> size in memory
        procedure :: size_in_mem => ec_size_in_mem
end type

!> one box in the hashes space.
type rl_extended_cluster_hashedspace_box
    !> number of atoms that can touch this box
    integer :: n
    !> indices to atoms that can touch this box
    integer, dimension(:), allocatable :: ind
end type
!> helper that divides space into sensible regions counts atoms.
type rl_extended_cluster_hashedspace
    !> how many boxes are there in each direction
    integer :: nx=-rl_hugeint,ny=-rl_hugeint,nz=-rl_hugeint
    !> bounding sphere of a small box
    real(r8) :: bounding_sphere_radius=-rl_huge
    !> boxes with lists
    type(rl_extended_cluster_hashedspace_box), dimension(:), allocatable :: box
    !> largest number of points in list
    integer :: max_n_pts=-rl_hugeint
    !> which cutoff was used to generate this?
    real(r8) :: cutoff=-rl_huge
    !> temporary buffer that holds positions
    real(r8), dimension(:,:), allocatable :: temp_cartesian_coordinate
    !> temporary buffer that holds unitcell index
    integer, dimension(:), allocatable :: temp_unitcell_index
    !> temporary buffer that holds cluster index
    integer, dimension(:), allocatable :: temp_cluster_index
    !> temporary buffer that holds squared distance to the point
    real(r8), dimension(:), allocatable :: temp_dist_sq
    !> temporary buffer that holds distances to the point
    real(r8), dimension(:), allocatable :: temp_dist
    contains
        !> create
        procedure :: generate=>create_hashedspace
        !> get the index of the box some point is in
        procedure :: box_index
        !> return relevant atoms; coordinates, unit cell indices and original indices
        procedure :: fetch_relevant_atoms
end type

contains

!> create everything related to an an extended supercluster embedded in a periodic system.
!> well not everything everything, but a decent start.
subroutine generate(ec,p,cutoff,mw,verbosity,nosym,mem)
    !> extended cluster
    class(rl_extended_cluster), intent(out) :: ec
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> radial cutoff
    real(r8), intent(in) :: cutoff
    !> mpi helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> skip symmetry considerations
    logical, intent(in) :: nosym
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem

    real(r8), dimension(:,:), allocatable :: ext_atom
    integer, dimension(:,:), allocatable :: shell_member,perm
    integer, dimension(:), allocatable :: shell_ctr,shell_index
    type(rl_voronoi_diagram) :: voro
    real(r8) :: t0,t1
    integer :: np
    !init: block
    type(rl_spacegroup) :: esym
    !integer, parameter :: maxiter=2000 ! Or some other number. Should never be an issue.
    real(r8), dimension(3,3) :: m0
    real(r8), dimension(3,1) :: r0
    real(r8), dimension(3) :: v0
    real(r8) :: buffer
    integer, dimension(3) :: ssdim
    integer :: i,j,l,op,iter
    logical :: okdimension
    !buildcluster: block
    integer, parameter :: maxiter=4000 ! Or some other number. Should never be an issue.
    !real(r8), dimension(3,3) :: m0
    !real(r8), dimension(3) :: v0
    real(r8) :: f0,rminsq,rmaxsq
    integer, dimension(3) :: nrep
    !integer :: i,iter,ctr,ii,jj,kk
    integer :: ctr,ii,jj,kk
    !checksym: block
    type(rl_distancetable) :: dt
    type(rl_verletbox) :: bx
    real(r8), dimension(:,:), allocatable :: dr,dr0,dr1
    !real(r8), dimension(3) :: v0
    integer, dimension(:), allocatable :: di
    !integer :: op,ct1,ct2,sh,i,j,l,k
    integer :: ct1,ct2,sh,k
    !sortatoms: block
    !real(r8), dimension(:,:), allocatable :: dr
    !integer, dimension(:), allocatable :: di,dj
    integer, dimension(:), allocatable :: dj
    !integer :: i,j,k,l
    !store: block
    !real(r8), dimension(3) :: v0
    !real(r8) :: f0
    !integer :: i,j,k,sh,nv
    integer :: nv

    ! Set some basic things
    !init: block

        if ( verbosity .gt. 0 ) then
            ! start the local timer
            t0=mpi_wtime()
        endif

        ! First thing I need to have is a Wigner-Seitz supercell, large enough to
        ! fit the naive cluster. I will make it simple for now, I think. Or a little
        ! worse, depending on how fast things are. This gets a little convoluted, but
        ! first I need the set of symmetry operations for the empty lattice.
        r0=0.0_r8
        call esym%generate(&
            latticevectors=p%latticevectors,&
            timereversal=.true.,&
            fractional_coordinate=r0,&
            species=[1],symmorphic=.true.,tol=1E-5_r8,&
            verbosity=-1,mw=mw,nosym=nosym,mem=mem)

        ! And to confuse things further, I want the symmorphic part of the space
        ! group as well. Think this is the only thing that is valid. Maybe.
        call ec%pointgroup%generate(&
            latticevectors=p%latticevectors,&
            timereversal=.true.,&
            fractional_coordinate=p%fractional_coordinate,&
            species=p%species,symmorphic=.true.,tol=1E-5_r8,&
            verbosity=-1,mw=mw,nosym=nosym,mem=mem)
        !call ec%pointgroup%get_character_table(verbosity=-1,mw=mw)

        ! Buffer for the radius to be on the safe side.
        buffer=rl_bounding_sphere_of_box(p%latticevectors)

        ! Now, get a series of increasingly larger supercells. For each of these cells,
        ! check if the inscribed sphere in the Voronoi cell is large enough to safely
        ! fit the cutoff. As an additional safety, check that the Voronoi cell satisifies
        ! the symmetries of the parent lattice.

        !@TODO Make this parallel. Later problem. Fast anyway and has constant cost.
        ssdim=1
        r0=0.0_r8
        iterl: do iter=1,maxiter
            ! Sanity check
            if ( iter .eq. maxiter ) then
                call rl_stop_gracefully(['Could not find Wigner-Seitz supercell. Odd.'],rl_exitcode_symmetry,mw%comm)
            endif
            ! get dimensions for supercell
            do i=1,3
                m0(:,i)=p%latticevectors(:,i)*ssdim(i)
            enddo
            ! get the Voronoi diagram
            call voro%generate(r0,m0,verbosity=-1)

            ! If too small, get a larger cell.
            if ( voro%cell(1)%rmin .lt. cutoff+buffer ) then
                ! yup, too small
                ssdim=increment_dimensions_single(ssdim,p%latticevectors)
                cycle iterl
            endif

            ! test if it satisfies all the point symmetries of the lattice
            okdimension=.true.
            opl: do op=1,esym%n_operation
                do i=1,voro%cell(1)%n_node
                    l=0
                    v0=matmul(esym%op(op)%m,voro%cell(1)%node(:,i))
                    do j=1,voro%cell(1)%n_node
                        if ( rl_sqnorm(v0-voro%cell(1)%node(:,j)) .lt. rl_sqtol ) l=l+1
                    enddo
                    if ( l .ne. 1 ) then
                        okdimension=.false.
                        exit opl
                    endif
                enddo
            enddo opl

            ! Check if this was good enough.
            if ( okdimension ) then
                ! yup, large enough
                exit iterl
            else
                ! Not good enough, increase size of supercell
                ssdim=increment_dimensions_single(ssdim,p%latticevectors)
                cycle iterl
            endif
        enddo iterl

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got a voronoi cell! dim: '//tochar(ssdim)//' (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block init

    ! Now generate the initial, huge, cluster.
    !buildcluster: block

        ! First find the number of repetitions needed to generate a massive cluster
        ! that I will later slice with the Voronoi cell.
        nrep=0
        do iter=1,maxiter
            do i=1,3
                m0(:,i)=p%latticevectors(:,i)*(2*nrep(i)+1)
            enddo
            f0=rl_inscribed_sphere_in_box(m0)
            if ( f0 .gt. voro%cell(1)%rmax ) then
                ! we are done!
                exit
            else
                ! pick more repetitions
                nrep=increment_dimensions_double(nrep,p%latticevectors)
            endif
        enddo

        ! First count how many things we get
        rminsq=voro%cell(1)%rmin**2
        rmaxsq=voro%cell(1)%rmax**2
        ctr=0
        do ii=-nrep(1),nrep(1)
        do jj=-nrep(2),nrep(2)
        do kk=-nrep(3),nrep(3)
            do i=1,p%n_atom
                v0=[ii,jj,kk]
                v0=matmul(p%latticevectors,v0)+p%cartesian_coordinate(:,i)
                f0=rl_sqnorm(v0)
                ! Now start checking. First the cheap tests
                if ( f0 .gt. rmaxsq ) then
                    ! way outside
                    cycle
                endif
                if ( f0 .lt. rminsq ) then
                    ! safely inside
                    ctr=ctr+1
                    cycle
                endif
                ! Ok not sure. have to check seriously.
                if ( voro%cell(1)%is_point_inside(v0,rl_tol) ) then
                    ! yup, was inside
                    ctr=ctr+1
                    cycle
                else
                    ! not inside
                    cycle
                endif
            enddo
        enddo
        enddo
        enddo

        ! Make space for the extended atoms
        np=ctr
        allocate(ext_atom(3,np))
        ext_atom=0.0_r8

        ! And populate
        ctr=0
        do ii=-nrep(1),nrep(1)
        do jj=-nrep(2),nrep(2)
        do kk=-nrep(3),nrep(3)
            do i=1,p%n_atom
                v0=[ii,jj,kk]
                v0=matmul(p%latticevectors,v0)+p%cartesian_coordinate(:,i)
                f0=rl_sqnorm(v0)
                ! Now start checking. First the cheap tests
                if ( f0 .gt. rmaxsq ) then
                    ! way outside
                    cycle
                endif
                if ( f0 .lt. rminsq ) then
                    ! safely inside
                    ctr=ctr+1
                    ext_atom(:,ctr)=v0
                    cycle
                endif
                ! Ok not sure. have to check seriously.
                if ( voro%cell(1)%is_point_inside(v0,rl_tol) ) then
                    ! yup, was inside
                    ctr=ctr+1
                    ext_atom(:,ctr)=v0
                    cycle
                else
                    ! not inside
                    cycle
                endif
            enddo
        enddo
        enddo
        enddo

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... built supercluster, n=',tochar(np),' (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block buildcluster

    ! Now check the symmetry of this humongous cluster of atoms
    !checksym: block

        ! Sort the points into boxes, useful later for fast lookup.
        ! The number of boxes chosen via black magic.
        i=ceiling( real(np,r8)**(1.0_r8/3.0_r8) )+5
        i=max(i,10)
        i=min(i,50)
        call bx%generate(ext_atom,[i,i,i],mem)

        ! Check with all the possible operations. Why not?
        ! Can see later if it is enough to check with the symmorphic.
        allocate(dr(3,np))
        allocate(perm(np,ec%pointgroup%n_operation))
        perm=0
        dr=0.0_r8
        ct1=0
        do op=1,ec%pointgroup%n_operation
            ! rotate all the atoms
            call dgemm('N','N',3,np,3,1.0_r8,ec%pointgroup%op(op)%m,3,ext_atom,3,0.0_r8,dr,3)
            ! Use the Verlet boxes to check where each point ended up.
            do i=1,np
                ! make parallel over particles. Not sure if actual speedup.
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                ! locate index of the transformed point in the original array
                k=bx%locate(ext_atom,dr(:,i))
                if ( k .gt. 0 ) then
                    ct1=ct1+1
                    perm(i,op)=k
                else
                    ! Report if I find bad operation
                    ct1=0
                endif
            enddo
        enddo
        deallocate(dr)

        ! Add up from all ranks.
        call mw%allreduce('sum',ct1)
        call mw%allreduce('sum',perm)

        ! Test that it went fine
        if ( ct1 .ne. np*ec%pointgroup%n_operation ) then
            call rl_stop_gracefully(['Huge cluster is not invariant under the symmorphic group. It should be.'],rl_exitcode_symmetry,mw%comm)
        endif

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... built permutations (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Good. Now generate coordination shells from this.
        call rl_coordination_shells_from_permutation_list(perm,shell_ctr,shell_member,shell_index,mem,mw=mw)
        deallocate(perm)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... built initial shells (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Get the naive cluster. Only shells that are touched by the naive
        ! cluster gets stored, the rest are thrown away.
        call dt%generate(p%fractional_coordinate,p%latticevectors,cutoff,-1,mw)
        l=sum(dt%particle(:)%n)
        allocate(dr0(3,l))
        dr0=0.0_r8
        l=0
        do i=1,dt%np
            v0=rl_clean_fractional_coordinates(p%fractional_coordinate(:,i)+0.5_r8)-0.5_r8
            v0=matmul(p%latticevectors,v0)
            do j=1,dt%particle(i)%n
                l=l+1
                dr0(:,l)=dt%particle(i)%v(:,j)+v0
            enddo
        enddo
        call rl_return_unique(dr0,dr1)

        ! see which shells are populated
        l=size(shell_ctr)
        allocate(di(l))
        di=0
        do i=1,size(dr1,2)
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            k=bx%locate(ext_atom,dr1(:,i))
            if ( k .gt. 0 ) then
                j=shell_index(k)
                di(j)=di(j)+1
            endif
        enddo
        call mw%allreduce('sum',di)

        ! Sanity test, should never trigger.
        if ( sum(di) .ne. size(dr1,2) ) then
            call rl_stop_gracefully(['Cluster does not contain the distance table.'],rl_exitcode_symmetry,mw%comm)
        endif
        deallocate(dr1)

        ! Now store the pruned points, and do all the symmetry and shell division once again
        ! but now on the pruned set that is significantly smaller. There is a good reason I
        ! do this twice, this makes sure that I only prune complete coordination shells,
        ! as soon as you have a half-full shell things will get very wonky later.

        ! First count the pruned number of points
        ct2=np
        ct1=0
        do sh=1,size(shell_ctr)
            if ( di(sh) .eq. 0 ) cycle
            ct1=ct1+shell_ctr(sh)
        enddo
        ! Then store the pruned points
        allocate(dr(3,ct1))
        l=0
        do sh=1,size(shell_ctr)
            if ( di(sh) .eq. 0 ) cycle
            do i=1,shell_ctr(sh)
                l=l+1
                dr(:,l)=ext_atom(:,shell_member(i,sh))
            enddo
        enddo
        deallocate(ext_atom)
        np=l
        allocate(ext_atom(3,np))
        ext_atom=dr

        ! Divide into boxes again
        call bx%destroy(mem)
        i=ceiling( real(np,r8)**(1.0_r8/3.0_r8) )+5
        i=max(i,10)
        i=min(i,50)
        call bx%generate(ext_atom,[i,i,i],mem)
        ! Reduce by symmetry, again
        allocate(perm(np,ec%pointgroup%n_operation))
        perm=0
        dr=0.0_r8
        ct1=0
        do op=1,ec%pointgroup%n_operation
            call dgemm('N','N',3,np,3,1.0_r8,ec%pointgroup%op(op)%m,3,ext_atom,3,0.0_r8,dr,3)
            do i=1,np
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                k=bx%locate(ext_atom,dr(:,i))
                if ( k .gt. 0 ) then
                    ct1=ct1+1
                    perm(i,op)=k
                else
                    ct1=0
                endif
            enddo
        enddo
        call bx%destroy(mem)
        deallocate(dr)
        call mw%allreduce('sum',ct1)
        call mw%allreduce('sum',perm)

        ! Test that it went fine
        if ( ct1 .ne. np*ec%pointgroup%n_operation ) then
            call rl_stop_gracefully(['Pruned cluster is not invariant under the symmorphic group. It should be.'],rl_exitcode_symmetry,mw%comm)
        endif

        call mem%deallocate(shell_ctr,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_member,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_index,persistent=.true.,scalable=.false.)

        ! And divide into shells again
        call rl_coordination_shells_from_permutation_list(perm,shell_ctr,shell_member,shell_index,mem,mw=mw)
        deallocate(perm)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... pruned shells ',tochar(ct2),' > ',tochar(np),' (',tochar(t1-t0),'s)'
            t0=t1
        endif
        deallocate(di)
    !end block checksym

    ! It makes life a lot easier later on if all the atoms are
    ! in a consistent order, that way I don't have to index(index)
    ! all the time.
    !sortatoms: block

        allocate(di(np))
        allocate(dj(np))
        allocate(dr(3,np))
        di=0
        dj=0
        l=0
        do i=1,size(shell_ctr)
        do j=1,shell_ctr(i)
            l=l+1
            k=shell_member(j,i)
            ! This is now a permutation such that old index -> new index
            di(l)=k
            dj(k)=l
        enddo
        enddo

        ! Apply this permutation everywhere
        do i=1,np
            dr(:,i)=ext_atom(:,di(i))
        enddo
        ext_atom=dr

        do i=1,size(shell_ctr)
        do j=1,shell_ctr(i)
            k=shell_member(j,i)
            shell_member(j,i)=dj(k)
        enddo
        enddo

        ! cleanup
        deallocate(di,dj,dr)
    !end block sortatoms

    ! Now store all this information in the handy object. Also perform plenty of sanity
    ! checks to make life easier to debug in the future.
    !store: block

        ! Number of extended atoms
        ec%n_extended_atom = np
        ! Coordinates
        allocate(ec%cartesian_coordinate(3,np))
        ec%cartesian_coordinate=ext_atom

        ! Index in the unit cell
        allocate(ec%index_unit_cell(np))
        ec%index_unit_cell=0
        do i=1,np
            v0=matmul(p%inv_latticevectors,ec%cartesian_coordinate(:,i))
            v0=v0-anint(v0)
            v0=rl_clean_fractional_coordinates(v0)
            k=0
            do j=1,p%n_atom
                if ( rl_sqnorm(v0-p%fractional_coordinate(:,j)) .lt. rl_sqtol ) then
                    k=j
                    exit
                endif
            enddo
            ! Make sure I found this atom in the unit cell
            if ( k .eq. 0 ) then
                call rl_stop_gracefully(['Could not find extended atoms in unit cell.'],rl_exitcode_symmetry,mw%comm)
            else
                ec%index_unit_cell(i)=j
            endif

            ! Sanity check to be really sure things make sense.
            v0=ec%cartesian_coordinate(:,i)-p%cartesian_coordinate(:,ec%index_unit_cell(i))
            v0=matmul(p%inv_latticevectors,v0)
            if ( rl_sqnorm(v0-anint(v0)) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['Clearly I do no understand vectors.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo

        ! Number of coordination shells
        ec%n_extended_shell=size(shell_ctr)

        ! Store information about the shells
        allocate(ec%shell(ec%n_extended_shell))
        do sh=1,ec%n_extended_shell
            ! Store some indices of the shell
            nv=shell_ctr(sh)
            ec%shell(sh)%n_extended_atom_in_shell = nv
            allocate(ec%shell(sh)%index_ext_atom_in_shell( nv ))
            ec%shell(sh)%index_ext_atom_in_shell = shell_member(1:shell_ctr(sh),sh)

            ! Check that shell only contains one kind of atom?
            j=ec%shell(sh)%index_ext_atom_in_shell(1)
            j=ec%index_unit_cell(j)
            j=p%species(j)
            do i=1,ec%shell(sh)%n_extended_atom_in_shell
                k=ec%shell(sh)%index_ext_atom_in_shell(i)
                k=ec%index_unit_cell(k)
                k=p%species(k)
                if ( j .ne. k ) then
                    call rl_stop_gracefully(['Shell contains member of different species. Should never happen.'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo
        enddo

        ! Finally, locate the unit cell atoms in the cluster for easy lookup later in life.
        allocate(ec%locate_unit_cell_atom_in_cluster(p%n_atom))
        ec%locate_unit_cell_atom_in_cluster=0
        do i=1,p%n_atom
            v0=rl_clean_fractional_coordinates(p%fractional_coordinate(:,i)+0.5_r8)-0.5_r8
            v0=matmul(p%latticevectors,v0)
            k=-1
            do j=1,ec%n_extended_atom
                if ( rl_sqnorm(v0-ec%cartesian_coordinate(:,j)) .lt. rl_sqtol ) then
                    k=j
                    exit
                endif
            enddo
            if ( k .le. 0 ) then
                call rl_stop_gracefully(['Clearly I do no understand vectors.'],rl_exitcode_symmetry,mw%comm)
            else
                ec%locate_unit_cell_atom_in_cluster(i)=k
            endif
        enddo

        ! Some cleanup
        call mem%deallocate(shell_ctr,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_member,persistent=.true.,scalable=.false.)
        call mem%deallocate(shell_index,persistent=.true.,scalable=.false.)

        ! And we are done.
        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... done constructing initial cluster (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block store

    ! Check that memory is cleared properly.
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
end subroutine

!> increase supercell dimensions in reasonably clever way
function increment_dimensions_single(dimin,box) result(dimut)
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
        m0(:,j)=box(:,j)*dimin(j)
    enddo
    ff0=rl_inscribed_sphere_in_box(m0)
    ! Increment the dimension that gives the biggest increase in radii
    f0=0.0_r8
    dimut=0
    do i=1,3
        di=dimin
        di(i)=di(i)+1
        do j=1,3
            m0(:,j)=box(:,j)*di(j)
        enddo
        f1=rl_inscribed_sphere_in_box(m0)
        if ( f1 .gt. f0 .and. abs(f1-ff0) .gt. rl_tol ) then
            dimut=di
            f0=f1
        endif
    enddo

    ! if nothing helped, increment the lowest number
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

!> increase supercell dimensions in reasonably clever way
function increment_dimensions_double(dimin,box) result(dimut)
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
    ! Increment the dimension that gives the biggest increase in radii
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

    ! if nothing helped, increment the lowest number
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

!> returns the hashing of space
subroutine create_hashedspace(hs,ec,p,cutoff,mw)
    !> space divided into boxes
    class(rl_extended_cluster_hashedspace), intent(out) :: hs
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> cutoff
    real(r8), intent(in) :: cutoff
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw

    integer, dimension(3) :: n_division

    !init: block
    integer, parameter :: approximate_ncells=200
    real(r8), dimension(3,3) :: m0
    real(r8) :: fa,fb,fc,f0,minvol
    integer, dimension(3) :: div
    integer :: na,nb,nc,ii,jj,kk,i

    !chopup: block
    real(r8), dimension(:,:), allocatable :: centers
    real(r8), dimension(3) :: v0,v1
    !real(r8) :: rcsq,f0
    real(r8) :: rcsq
    integer, dimension(:,:), allocatable :: di
    integer, dimension(:), allocatable :: ctr
    !integer :: i,j,k,l,ll
    integer :: j,k,l,ll
    integer :: nbox

    ! Decide on how to divide
    !init: block
        ! A somewhat arbitrary number, but sensibly small.
        ! 200 cells with 5000 atoms yields roughly 1 million integers to store,
        ! if the cutoff is really long. If the cutoff is sensible, which it most
        ! likely is, the storage will be much less.

        ! So, I want to divide space into equal-ish boxes. Ideally, the boxes should be
        ! square, but that is not going to happen. Instead I will guesstimate what a good
        ! number of divisions is, to yield boxes that are not too strangely shaped.
        fa=norm2(p%latticevectors(:,1))
        fb=norm2(p%latticevectors(:,2))
        fc=norm2(p%latticevectors(:,3))
        f0=(approximate_ncells/(fa*fb*fc))**(1.0_r8/3.0_r8)
        na=ceiling(fa*f0)
        nb=ceiling(fb*f0)
        nc=ceiling(fc*f0)
        na=max(na,3)
        nb=max(nb,3)
        nc=max(nc,3)
        ! This is a decent guesstimate. Now search around these parameters
        ! to pick the one that gives the smallest boxes.
        minvol=rl_huge
        do ii=-2,2
        do jj=-2,2
        do kk=-2,2
            div=[na+ii,nb+jj,nc+kk]
            if ( product(div) .gt. 2*na*nb*nc ) cycle   ! Don't want too many divisions.
            do i=1,3
                m0(:,i)=p%latticevectors(:,i)/real(div(i),r8)
            enddo
            f0=rl_bounding_sphere_of_box(m0)

            if ( f0 .lt. minvol ) then
                minvol=f0
                n_division=div
            endif
        enddo
        enddo
        enddo

        ! Also make note how many boxes there are in each direction
        hs%nx=n_division(1)
        hs%ny=n_division(2)
        hs%nz=n_division(3)
        ! Store som aux information.
        do i=1,3
            m0(:,i)=p%latticevectors(:,i)/real(n_division(i),r8)
        enddo
        hs%bounding_sphere_radius=rl_bounding_sphere_of_box(m0)
        hs%cutoff=cutoff
    !end block init

    ! Start dividing
    !chopup: block

        nbox=product(n_division)
        allocate(centers(3,nbox))
        allocate(ctr(nbox))
        centers=0.0_r8
        ctr=0

        ! Collect the coordinates of the centers of the boxes
        l=0
        do i=1,hs%nx
        do j=1,hs%ny
        do k=1,hs%nz
            l=l+1
            ! this should yield the center of the box
            v0=[i,j,k]-0.5_r8
            v0(1)=v0(1)/real(hs%nx,r8)
            v0(2)=v0(2)/real(hs%ny,r8)
            v0(3)=v0(3)/real(hs%nz,r8)
            v0=v0-0.5_r8
            ! idiot check that I can index properly, unit-testing the locator function.
            ll=hs%box_index(v0)
            if ( ll .ne. l ) then
                call rl_stop_gracefully(['I do not know how to index'],rl_exitcode_symmetry,mw%comm)
            endif
            centers(:,l)=matmul(p%latticevectors,v0)
        enddo
        enddo
        enddo

        ! Start counting atoms within the cutoff for each center.
        rcsq=(cutoff+hs%bounding_sphere_radius+1E-6_r8)**2
        ctr=0
        do i=1,nbox
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            v0=centers(:,i)
            do j=1,ec%n_extended_atom
                v1=v0-ec%cartesian_coordinate(:,j)
                f0=v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
                if ( f0 .lt. rcsq ) then
                    ctr(i)=ctr(i)+1
                endif
            enddo
        enddo
        call mw%allreduce('sum',ctr)

        ! Now do it again, but keep the indices
        j=maxval(ctr)
        allocate(di(j,nbox))
        di=0
        ctr=0
        do i=1,nbox
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            v0=centers(:,i)
            do j=1,ec%n_extended_atom
                v1=v0-ec%cartesian_coordinate(:,j)
                f0=v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)
                if ( f0 .lt. rcsq ) then
                    ctr(i)=ctr(i)+1
                    di( ctr(i),i )=j
                endif
            enddo
        enddo
        call mw%allreduce('sum',di)

        ! Now store this in the appropriate place
        allocate(hs%box(nbox))
        hs%max_n_pts=0
        do i=1,nbox
            hs%box(i)%n=count(di(:,i)>0)
            allocate(hs%box(i)%ind(hs%box(i)%n))
            hs%box(i)%ind=di(1:hs%box(i)%n,i)
            hs%max_n_pts=max(hs%max_n_pts,hs%box(i)%n)
        enddo

        ! And make space for buffers to collect things to
        allocate(hs%temp_cartesian_coordinate(3,hs%max_n_pts))
        allocate(hs%temp_unitcell_index(hs%max_n_pts))
        allocate(hs%temp_cluster_index(hs%max_n_pts))
        allocate(hs%temp_dist_sq(hs%max_n_pts))
        allocate(hs%temp_dist(hs%max_n_pts))
        hs%temp_cartesian_coordinate=-rl_huge
        hs%temp_unitcell_index=-rl_hugeint
        hs%temp_cluster_index=-rl_hugeint
        hs%temp_dist_sq=-rl_huge
        hs%temp_dist=-rl_huge

        ! And some cleanup!
        deallocate(di)
        deallocate(centers)
        deallocate(ctr)
    !end block chopup
end subroutine

!> return the index of a box from the coordinates
function box_index(hs,r) result(ijk)
    !> divided space
    class(rl_extended_cluster_hashedspace), intent(in) :: hs
    !> point, in fractional coordinates -0.5 < r < 0.5
    real(r8), dimension(3), intent(in) :: r
    !> index to box
    integer :: ijk

    real(r8), dimension(3) :: v0
    integer :: i,j,k

    v0=r+0.5_r8
    v0=max(v0,1E-15_r8)
    v0=min(v0,0.999999999999_r8)
    i=ceiling(v0(1)*hs%nx)
    j=ceiling(v0(2)*hs%ny)
    k=ceiling(v0(3)*hs%nz)
    ijk=k+(j-1)*hs%nz+(i-1)*hs%nz*hs%ny
end function

!> return a list of atoms and indices that are relevant for the current point
subroutine fetch_relevant_atoms(hs,ec,p,r,cutoff,n_point)
    !> divided space
    class(rl_extended_cluster_hashedspace), intent(inout) :: hs
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> point in fractional -0.5<r<0.5 coordinates
    real(r8), dimension(3), intent(in) :: r
    !> cutoff, fetch points within this distance from r, in Cartesian coordinates
    real(r8), intent(in) :: cutoff
    !> number of points I fetched
    integer, intent(out) :: n_point

    real(r8), dimension(3) :: rcart,v,w
    real(r8) :: f0,rcsq
    integer :: ibox,i,j,l

    ! Figure out which box we are in?
    ibox=hs%box_index(r)
    if ( cutoff .gt. hs%cutoff ) then
        call rl_stop_gracefully(['Cutoff larger than generating cutoff, thinking was wrong.'],rl_exitcode_param)
    endif

    ! Initialize to nonsense
    hs%temp_cartesian_coordinate=-rl_huge
    hs%temp_unitcell_index=-rl_hugeint
    hs%temp_cluster_index=-rl_hugeint
    hs%temp_dist_sq=-rl_huge
    hs%temp_dist=-rl_huge

    ! I need some Cartesian coordinates as well
    rcart=matmul(p%latticevectors,r)

    ! Start fetching
    rcsq=cutoff**2
    l=0
    do i=1,hs%box(ibox)%n
        j=hs%box(ibox)%ind(i)
        v=ec%cartesian_coordinate(:,j)
        w=v-rcart
        f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
        if ( f0 .le. rcsq ) then
            l=l+1
            hs%temp_cartesian_coordinate(:,l)=v
            hs%temp_unitcell_index(l)=ec%index_unit_cell(j)
            hs%temp_cluster_index(l)=j
            hs%temp_dist_sq(l)=f0
            hs%temp_dist(l)=sqrt(f0)
        endif
    enddo
    ! Store the number of points fetched
    n_point=l
end subroutine

!> measure size in memory in bytes, approximately. Should be larger than the real usage, maybe.
function ec_size_in_mem(ec) result(mem)
    !> extended cluster
    class(rl_extended_cluster), intent(in) :: ec
    !> memory in bytes
    integer :: mem,i

    mem=0
    mem=mem+storage_size(ec)
    if ( allocated(ec%cartesian_coordinate) )  mem=mem+storage_size(ec%cartesian_coordinate)*size(ec%cartesian_coordinate)
    if ( allocated(ec%index_unit_cell) )       mem=mem+storage_size(ec%index_unit_cell)*size(ec%index_unit_cell)
    !if ( allocated(ec%confinement_cutoff_sq) ) mem=mem+storage_size(ec%confinement_cutoff_sq)*size(ec%confinement_cutoff_sq)
    !if ( allocated(ec%basis_index_offset) )    mem=mem+storage_size(ec%basis_index_offset)*size(ec%basis_index_offset)
    if ( allocated(ec%locate_unit_cell_atom_in_cluster) )    mem=mem+storage_size(ec%locate_unit_cell_atom_in_cluster)*size(ec%locate_unit_cell_atom_in_cluster)

    if ( allocated(ec%shell) ) then
        do i=1,size(ec%shell)
            mem=mem+storage_size(ec%shell(i))
            if ( allocated(ec%shell(i)%index_ext_atom_in_shell) ) mem=mem+storage_size(ec%shell(i)%index_ext_atom_in_shell)*size(ec%shell(i)%index_ext_atom_in_shell)
        enddo
    endif
    mem=mem/8
    mem=mem+ec%pointgroup%size_in_mem()
end function

end module
