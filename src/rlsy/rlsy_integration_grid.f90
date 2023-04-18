module rlsy_integration_grid
!!
!! This module chops a realspace grid into manageable batches. It is the same as
!!
!! Havu, V., Blum, V., Havu, P. & Scheffler, M.
!! Efficient O (N) integration for all-electron electronic structure calculation using numeric basis functions.
!! J. Comput. Phys. 228, 8367–8379 (2009).
!!
!! With some modifications. Not too much modifications. Ok quite a lot of changes, but the idea is the same.
!!
use rlsy_constants, only: r8,i8,rl_iou,rl_huge,rl_hugeint,rl_tol,rl_twopi,rl_sqtol,rl_tiny,&
    rl_exitcode_symmetry,rl_exitcode_mpi,rl_exitcode_param,rl_exitcode_memory
use rlsy_memtracker, only: rl_memtracker
use rlsy_sorting, only: rl_qsort
use rlsy_helpers, only: rl_chop,tochar,rl_sqnorm,rl_clean_fractional_coordinates,norm2
use rlsy_linalg, only: rl_symmetric_eigensystem_3x3matrix
use rlsy_mpi_helper, only: rl_mpi_helper,rl_stop_gracefully,mpi_wtime,MPI_DOUBLE_PRECISION,MPI_INTEGER,mpi_wtime
use rlsy_crystalstructure, only: rl_crystalstructure
use rlsy_spacegroup, only: rl_spacegroup
use rlsy_symmetry_helper_functions, only: rl_coordination_shells_from_permutation_list,&
    rl_alltoall_distribute_2d_real,rl_alltoall_distribute_2d_int,rl_find_prime_factors,rl_find_divisors
use rlsy_extended_cluster, only: rl_extended_cluster,rl_extended_cluster_hashedspace
use rlsy_basis_set, only: rl_lcao_basis_set
use rlsy_free_atom, only: rl_free_atom
use rlsy_distancetable, only: rl_distancetable
use rlsy_partition_function, only: rl_stratmann_smoother_partition_function
use rlsy_integration_grid_helpers, only: rl_chunking_helper,rl_set_of_points,rl_distribute_sets,&
    rl_get_irr_points_per_rank,rl_split_set

implicit none
private
public :: rl_integration_grid
public :: rl_integration_batch
!public :: rl_integration_batch_irr
!public :: rl_integration_batch_full

!> things that are defined for all kinds of batches
type rl_integration_batch
    !> how many integration points in this batch
    integer :: n_point=-rl_hugeint
    !> radial coordinate from originating atom
    real(r8), dimension(:,:), allocatable :: radial_coordinate
    !> coordinate of integration point in Cartesian coordinates, folded to the first cell
    real(r8), dimension(:,:), allocatable :: folded_coordinate
    !> integration weight, sum(integration_weight*partition_function) = volume, approximately
    real(r8), dimension(:), allocatable :: integration_weight
    !> partition function
    real(r8), dimension(:), allocatable :: partition_function

    !> index to atom that owns the point
    integer, dimension(:), allocatable :: index_atom
    !> radial index used to generate the point
    integer, dimension(:), allocatable :: index_radial
    !> angular index used to generate the point
    integer, dimension(:), allocatable :: index_angular

    !> for pruning, it's great to know the center of mass
    real(r8), dimension(3) :: center_of_mass=-rl_huge
    !> for pruning, useful to have the radius of the bounding sphere centered at the center of mass
    real(r8) :: bounding_sphere_radius=-rl_huge
    !> number of atoms whose basis functions can touch this batch
    integer :: n_relevant_ext_atom=-rl_hugeint
    !> list of relevant atoms for this batch
    integer, dimension(:), allocatable :: relevant_ext_atom
!end type
!
!!> batch of integration points on the irreducible grid.
!type, extends(rl_integration_batch) :: rl_integration_batch_irr
    !> offset for semi-local irreducible index, irr_offset+local point index = semilocal index
    integer :: semilocal_irr_offset=-rl_hugeint
    !> index to irreducible atom
    integer, dimension(:), allocatable :: index_irr_atom

    !> information how to unfold each grid point to the full grid
    integer, dimension(:),   allocatable :: unfold_ctr
    integer, dimension(:,:), allocatable :: unfold_operation
    integer, dimension(:,:), allocatable :: unfold_index_radial
    integer, dimension(:,:), allocatable :: unfold_index_angular
    integer, dimension(:,:), allocatable :: unfold_atom
!end type
!
!!> batch of integration points of the full grid
!type, extends(rl_integration_batch) :: rl_integration_batch_full
    !> offset for semi-local irreducible index, irr_offset+local point index = semilocal index
    integer :: semilocal_full_offset=-rl_hugeint
end type

!> Small helper types to index basis functions globally without wasting a lot of memory.
type rl_integration_grid_index_atom
    integer, dimension(:), allocatable :: radial
end type
type rl_integration_grid_index
    type(rl_integration_grid_index_atom), dimension(:), allocatable, private :: atom
    contains
        procedure :: global_index
end type

!> Small helper to fetch radial or angular integration weights for a points. Can add the
!> actual coordinates here if we need to at some point.
type rl_integration_grid_weights_species_radial
    real(r8) :: radial_weight=-rl_huge
    real(r8), dimension(:), allocatable :: angular_weight
end type
type rl_integration_grid_weights_species
    real(r8) :: radial_scale=-rl_huge
    type(rl_integration_grid_weights_species_radial), dimension(:), allocatable :: rad
end type
type rl_integration_grid_weights
    type(rl_integration_grid_weights_species), dimension(:), allocatable, private :: species
    contains
        procedure :: n_radial=>fetch_n_radial
        procedure :: radial_scale=>fetch_radial_scale
        procedure :: radial_weight=>fetch_radial_weight
        procedure :: angular_weight=>fetch_angular_weight
end type

!> Everything related to integration grids. Well not everything, but most of it.
type rl_integration_grid
    !> global number of points
    integer :: n_point_global=-rl_hugeint
    !> global number of full points
    integer :: n_point_full_global=-rl_hugeint
    !> global number of irreducible points
    integer :: n_point_irr_global=-rl_hugeint
    !> global number of empty irreducible points, empty meaning zero partition function
    integer :: n_point_irr_zero_global=-rl_hugeint
    !> global number of chunks of points
    integer :: n_chunk_global=-rl_hugeint

    !> semi-local classification, each chunk has it's own communicator
    integer :: chunk_index=-rl_hugeint
    !> semi-local MPI helper that handles all ranks with the same chunk index.
    type(rl_mpi_helper) :: ml

    !> number of batches of irreducible points with nonzero weight.
    integer :: n_irr_batch=-rl_hugeint
    !> batch of irreducible points with nonzero weight. This is the normal set of points to loop over.
    type(rl_integration_batch), dimension(:), allocatable :: irr_batch

    !> number of batches of irreducible points with zero weight
    integer :: n_irr_zero_batch=-rl_hugeint
    !> batch of irreducible points with zero weight. You will likely not need to use these.
    type(rl_integration_batch), dimension(:), allocatable :: irr_zero_batch

    !> number of full batches
    integer :: n_full_batch=-rl_hugeint
    !> full batches (with nonzero weight)
    type(rl_integration_batch), dimension(:), allocatable :: full_batch

    !> largest number of basis functions per batch (locally)
    integer :: max_n_basis_per_batch=-rl_hugeint
    !> largest number of extended atoms per batch (locally)
    integer :: max_n_atom_per_batch=-rl_hugeint
    !> largest number of integration points per batch (locally)
    integer :: max_n_point_per_batch=-rl_hugeint

    !> global indexing helper
    type(rl_integration_grid_index) :: idx
    !> integration weight helper
    type(rl_integration_grid_weights) :: wts
    contains
        !> create the irreducible grid from the grid in AIMS
        procedure :: generate=>rl_partition_integration_grid
end type

contains

!> Partition the integration grid into batches
subroutine rl_partition_integration_grid(grid,p,sym,ec,basis,free,mw,mem,verbosity,&
    n_radial,n_angular,r_radial,r_angular,w_radial,w_angular,outer_partition_radius,scale_radial)
    !> actual grid
    class(rl_integration_grid), intent(out) :: grid
    !> crystal structure
    type(rl_crystalstructure), intent(in) :: p
    !> symmetry handle
    type(rl_spacegroup), intent(in) :: sym
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> free atom quantities
    type(rl_free_atom), intent(in) :: free
    !> mpi helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(rl_memtracker), intent(inout) :: mem
    !> grid information from AIMS
    integer, dimension(:), intent(in) :: n_radial
    integer, dimension(:,:), intent(in) :: n_angular
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:,:), intent(in) :: w_radial
    real(r8), dimension(:,:,:), intent(in) :: w_angular
    real(r8), dimension(:), intent(in) :: outer_partition_radius
    real(r8), dimension(:), intent(in) :: scale_radial
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! This should probably be made not parameters, but later problem.
    integer, parameter :: n_pts_per_batch=200
    integer, parameter :: max_n_irr_on_rank=1500
    ! timers
    real(r8) :: timer,t0,t1

    ! how many times to divide into chunks. Explain later.
    integer :: n_first_division

    ! Helper arrays for the division of the irreducible points
    integer :: n_irr_point_local
    integer :: n_irr_zero_point_local
    integer, dimension(:,:), allocatable :: irr_ind
    integer, dimension(:,:), allocatable :: irr_zero_ind

    !init: block
        integer :: i,ispc,iatm,irad
    !globind: block
        integer :: i_atom,i_species,i_radial,offset
    !sortweights: block
        integer :: is,ir,ia
    !irr_distribution: block
        type(rl_set_of_points), dimension(:), allocatable :: irr_set,irr_zero_set
        integer :: n_irr_set,n_irr_zero_set
        real(r8), dimension(:,:), allocatable :: irr_r,irr_zero_r
        integer, dimension(:,:), allocatable :: divisors_per_chunk
        integer, dimension(:), allocatable :: n_pts_per_chunk
        integer, dimension(:), allocatable :: n_divisors_per_chunk
        !integer :: local_my_chunk,i,j
        integer :: local_my_chunk,j
    !full_distribution: block
        type(rl_set_of_points), dimension(:), allocatable :: set
        integer :: n_set
        real(r8), dimension(:,:), allocatable :: local_r
        !integer, dimension(:,:), allocatable :: divisors_per_chunk
        !integer, dimension(:), allocatable :: n_pts_per_chunk
        !integer, dimension(:), allocatable :: n_divisors_per_chunk
        integer, dimension(:,:), allocatable :: local_ind
        !integer :: ibatch,ip,ctr,i,l
        integer :: ibatch,ip,ctr,l
        !integer :: local_npts,local_my_chunk
        integer :: local_npts
    !batchstat: block
        integer :: nbmax,namax,npmax
        !integer :: ib,i,l,ie,iu,is
        integer :: ib,ie,iu
    ! Set some basic things, start timers and so on.
    !init: block

        ! start timers
        if ( verbosity .gt. 0 ) then
            timer=mpi_wtime()
            t0=timer
            write(rl_iou,*) ''
            write(rl_iou,*) 'SYMMETRIZING AND PARTITIONING INTEGRATION GRID'
        endif

        ! count the total number of points? Might come handy eventually.
        i=0
        do iatm=1,p%n_atom
            ispc=p%species(iatm)
            do irad=1,n_radial(ispc)
                i=i+n_angular(irad,ispc)
            enddo
        enddo
        grid%n_point_global=i

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... initialized grid division (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Return the irreducible points, disitributed evenly across ranks in no
        ! particular order. They are stored in the index array irr_ind, which contains
        ! the bare minimum information to reconstruct the whole grid.
        call rl_get_irr_points_per_rank(p,sym,mw,mem,n_radial,n_angular,r_angular,&
            n_irr_point_local,irr_ind)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... reduced mesh with symmetry (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! At this stage, we have some number of irreducible points per rank and an
        ! index array that allows me to reconstruct the full grid. Now calculate
        ! the partition function for the irreducible points and divide points into
        ! those with or without any contribution to the partition function. The
        ! partition function will be calculated again, later, and used as a sanity
        ! test that I did not mess anything up while distributing points all over
        ! the place.
        call sort_by_partition_function(p,ec,basis,free,sym,mw,r_radial,r_angular,outer_partition_radius,&
            n_irr_point_local,n_irr_zero_point_local,irr_ind,irr_zero_ind)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... divided mesh with partition function (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now it gets a little tricky. I want to divide the irreducible points
        ! into pieces. Each of those pieces gets their own MPI communicator,
        ! such that every full point in that subcommunicator shares the irreducible
        ! with all the other points. Then memory is scalable, and I trade some
        ! complexity now for simplicity later.

        ! First I need to know the total number of irreducible:
        grid%n_point_irr_global=n_irr_point_local
        grid%n_point_irr_zero_global=n_irr_zero_point_local
        call mw%allreduce('sum',grid%n_point_irr_global)
        call mw%allreduce('sum',grid%n_point_irr_zero_global)

        ! Then decide on the number of chunks.
        i=0
        do
            ! I can not divide onto more ranks than I have.
            if ( 2**(i+1) .gt. mw%n ) then
                n_first_division=i
                exit
            endif
            ! This might be enough division?
            if ( grid%n_point_irr_global/(2**i) .lt. max_n_irr_on_rank ) then
                n_first_division=i
                exit
            endif
            ! if not, try more division!
            i=i+1
        enddo
        ! Now I know the number of chunks:
        grid%n_chunk_global=2**(n_first_division)

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '    total n point:',grid%n_point_global
            write(rl_iou,*) '      irr n point:',grid%n_point_irr_global,' ',tochar(100.0_r8*grid%n_point_irr_global/grid%n_point_global),'%'
            write(rl_iou,*) ' zero irr n point:',grid%n_point_irr_zero_global
            write(rl_iou,*) '         n chunks:',grid%n_chunk_global
        endif
    !end block init

    ! Now for a quick intermission, it is nice to figure out a way to index the points. Sort
    ! that out before we continue with dividing the points into batches and all that.
    !globind: block
        ! I need a consistent way to turn
        !   index_atom
        !   index_radial
        !   index_angular
        ! into a linear index that runs from 1 to the total number
        ! of points in the mesh. It should also be reasonably fast.
        ! Start by setting the offsets to nothing.
        allocate(grid%idx%atom(p%n_atom))
        do i_atom=1,p%n_atom
            i_species=p%species(i_atom)
            allocate(grid%idx%atom(i_atom)%radial(n_radial(i_species)))
            grid%idx%atom(i_atom)%radial=0
        enddo

        ! Start counting basis functions and store the offset.
        offset=0
        do i_atom=1,p%n_atom
            i_species=p%species(i_atom)
            do i_radial=1,n_radial(i_species)
                grid%idx%atom(i_atom)%radial(i_radial)=offset
                offset=offset+n_angular(i_radial,i_species)
            enddo
        enddo
        ! The way this works is that
        ! grid%idx%atom(i_atom)%radial(i_radial)+i_angular = global index!
        ! It should not need too much memory. n_atom*n_radial, that should
        ! be somewhat manageable.
    !end block globind

    ! Another intermission, it is nice to be able to fetch the partial integration weights
    ! for an arbitrary point, especially when debugging.
    !sortweights: block

        ! Make a little space
        allocate(grid%wts%species(p%n_species))
        do is=1,p%n_species
            grid%wts%species(is)%radial_scale=scale_radial(is)
            allocate(grid%wts%species(is)%rad( n_radial(is) ))
            do ir=1,n_radial(is)
                grid%wts%species(is)%rad(ir)%radial_weight=w_radial(ir,is)
                allocate(grid%wts%species(is)%rad(ir)%angular_weight( n_angular(ir,is) ))
                do ia=1,n_angular(ir,is)
                    grid%wts%species(is)%rad(ir)%angular_weight(ia)=w_angular(ia,ir,is)
                enddo
            enddo
        enddo
    !end block sortweights

    ! Now I start by distributing the irreducible points, both the zero and nonzero.
    ! I don't really know when the nonzero are needed, but neat to stash them away
    ! somewhere for safekeeping.
    !irr_distribution: block
        ! sets of points divided into pieces

        ! So, I know how many chunks I am going to divide the set of points into.
        ! Each chunk then gets it's own communicator. I need to know how many ranks
        ! each communicator will get. Then there is the whole matter of getting an
        ! even number batches per rank, and then an even number of points per batch.
        ! I think this will make load-balancing easier in the end.
        call rl_find_divisors(mw%n,grid%n_point_irr_global,grid%n_chunk_global,n_pts_per_batch,&
                              n_pts_per_chunk,n_divisors_per_chunk,divisors_per_chunk)

        if ( verbosity .gt. 0 ) then
            write(*,*) 'FIXME OUTPUT GRID THINGS LARGE SYSTEM'
            write(rl_iou,*) '        pts/chunk:',tochar(n_pts_per_chunk)
            do i=1,min(grid%n_chunk_global,10)
            write(rl_iou,*) '         divisors: ',tochar(i),' div ',tochar(divisors_per_chunk(1:n_divisors_per_chunk(i),i))
            enddo
        endif

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... decided on division (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! So, now I have a decent idea of how the data should be sliced, go
        ! on with the actual slicing. The first few steps are done in parallel,
        ! to generate the chunks with fixed size.
        call parallel_minmax_predivision(p,mw,r_radial,r_angular,irr_r,irr_ind,n_irr_point_local,&
             ndiv=n_first_division,verbosity=verbosity,n_pts_per_chunk=n_pts_per_chunk,&
             local_my_chunk=local_my_chunk)
        ! Quick sanity test to make sure this went as expected:
        if ( mw%r+1 .le. grid%n_chunk_global ) then
        if ( n_pts_per_chunk(mw%r+1) .ne. n_irr_point_local ) then
            call rl_stop_gracefully(['The points did not get divided as one would expect. Stopping now.'],rl_exitcode_symmetry,mw%comm)
        endif
        endif

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... split into chunks (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! If we reach this point safely, it means I have chunks the huge mesh down into smaller
        ! pieces. The smaller pieces can be handled serially, pretty much, since mesh division
        ! is NlogN, kinda. There is also the magic parameter to tune, i.e. how large the chunks
        ! are maximally allowed to be. If that is a smaller number, it gets more parallel.
        ! Will have to test on a large computer to see what makes the most sense, and how it
        ! scales in practice. So, this is a decent place to figure out which ranks deal with
        ! which chunk and then split the communicator.
        do i=1,mw%n
            j=mod(i-1,grid%n_chunk_global)+1
            if ( mw%r+1 .eq. i ) grid%chunk_index=j
        enddo
        call mw%split(grid%ml,grid%chunk_index)

        ! Test that things are the way I think they are.
        if ( grid%ml%r .eq. 0 ) then
        if ( n_pts_per_chunk(grid%chunk_index) .ne. n_irr_point_local ) then
            call rl_stop_gracefully(['The points did not get divided as one would expect. Stopping now.'],rl_exitcode_symmetry,mw%comm)
        endif
        endif

        ! Anyhow, in all the frustrating index algebra above, I figured out how many times each
        ! chunk should be divided, and by which factor it should be divided each time. So now do
        ! just that.
        if ( grid%ml%r .eq. 0 ) then
            call serial_minmax_division(irr_r,irr_ind,n_irr_point_local,&
                divisors_per_chunk(1:n_divisors_per_chunk(grid%chunk_index),grid%chunk_index),&
                irr_set,n_irr_set,verbosity)
        else
            ! allocate a token empty set maybe?
            allocate(irr_set(1))
            n_irr_set=0
        endif
        ! Distribute the sets semilocally
        call rl_distribute_sets(irr_set,n_irr_set,grid%ml,mw)
        ! Create batches from the sets
        call convert_set_into_irreducible_batches(p,sym,irr_set,n_irr_set,&
        r_radial,r_angular,w_radial,w_angular,n_angular,&
        grid%irr_batch,grid%n_irr_batch,grid%ml,mw)

        ! That was the division of the nonzero irreducible points. Now do the same thing to the
        ! irreducible points with zero partition function, in case I want them later at some
        ! point in the future. They will likely be really weirdly distributed, so I pick the
        ! batch size to be 1/4 of the usual batch size for these guys.
        deallocate(n_pts_per_chunk)
        deallocate(n_divisors_per_chunk)
        deallocate(divisors_per_chunk)
        call rl_find_divisors(mw%n,grid%n_point_irr_zero_global,grid%n_chunk_global,n_pts_per_batch/2,&
                              n_pts_per_chunk,n_divisors_per_chunk,divisors_per_chunk)

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '        pts/chunk:',tochar(n_pts_per_chunk)
            do i=1,grid%n_chunk_global
            write(rl_iou,*) '         divisors: ',tochar(i),' div ',tochar(divisors_per_chunk(1:n_divisors_per_chunk(i),i))
            enddo
        endif
        call parallel_minmax_predivision(p,mw,r_radial,r_angular,irr_zero_r,irr_zero_ind,n_irr_zero_point_local,&
             ndiv=n_first_division,verbosity=verbosity,n_pts_per_chunk=n_pts_per_chunk,&
             local_my_chunk=local_my_chunk)

        ! Quick sanity test to make sure this went as expected:
        if ( mw%r+1 .le. grid%n_chunk_global ) then
        if ( n_pts_per_chunk(mw%r+1) .ne. n_irr_zero_point_local ) then
            call rl_stop_gracefully(['The zero-weight points did not get divided as one would expect. Stopping now.'],rl_exitcode_symmetry,mw%comm)
        endif
        endif
        ! Divide serially
        if ( grid%ml%r .eq. 0 ) then
            call serial_minmax_division(irr_zero_r,irr_zero_ind,n_irr_zero_point_local,&
                divisors_per_chunk(1:n_divisors_per_chunk(grid%chunk_index),grid%chunk_index),&
                irr_zero_set,n_irr_zero_set,verbosity)
        else
            ! allocate a token empty set maybe?
            allocate(irr_zero_set(1))
            n_irr_zero_set=0
        endif
        ! Distribute the sets semilocally
        call rl_distribute_sets(irr_zero_set,n_irr_zero_set,grid%ml,mw)
        ! Create batches from the sets
        call convert_set_into_irreducible_batches(p,sym,irr_zero_set,n_irr_zero_set,&
        r_radial,r_angular,w_radial,w_angular,n_angular,&
        grid%irr_zero_batch,grid%n_irr_zero_batch,grid%ml,mw)

        !@TODO Here it would be great to move the batches around between ranks.
        ! What I would want, ideally, is to have the batches with the same nonzero
        ! atoms on the same rank to as large degree as possible. That will help
        ! a lot later in the density update, so that I only have to fetch the
        ! densitymatrix for a small subset of atom pairs for each rank. It will
        ! speed things up I believe.

        ! A little cleanup
        if ( allocated(irr_r)                ) deallocate(irr_r)
        if ( allocated(irr_ind)              ) deallocate(irr_ind)
        if ( allocated(irr_set)              ) deallocate(irr_set)
        if ( allocated(irr_zero_r)           ) deallocate(irr_zero_r)
        if ( allocated(irr_zero_ind)         ) deallocate(irr_zero_ind)
        if ( allocated(irr_zero_set)         ) deallocate(irr_zero_set)
        if ( allocated(divisors_per_chunk)   ) deallocate(divisors_per_chunk)
        if ( allocated(n_pts_per_chunk)      ) deallocate(n_pts_per_chunk)
        if ( allocated(n_divisors_per_chunk) ) deallocate(n_divisors_per_chunk)

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... generated irreducible batches from points (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block irr_distribution

    ! That was the irreducible batches. Thing is, I will likely need the set of
    ! batches corresponding to the full grid as well. I'm not sure what the best
    ! thing to do is here. In case there is no symmetry, the set of full and
    ! irreducible are the same, that needs to be handled smoothly somehow.
    !full_distribution: block

        ! First count the total number of full nonzero points
        ctr=0
        do ibatch=1,grid%n_irr_batch
        do ip=1,grid%irr_batch(ibatch)%n_point
            ctr=ctr+grid%irr_batch(ibatch)%unfold_ctr(ip)
        enddo
        enddo
        call mw%allreduce('sum',ctr)
        ! ctr now holds the total number of points.
        call rl_find_divisors(mw%n,ctr,grid%n_chunk_global,n_pts_per_batch,&
                              n_pts_per_chunk,n_divisors_per_chunk,divisors_per_chunk)

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '        pts/chunk:',tochar(n_pts_per_chunk)
            do i=1,min(grid%n_chunk_global,10)
                write(rl_iou,*) '         divisors: ',tochar(i),' div ',tochar(divisors_per_chunk(1:n_divisors_per_chunk(i),i))
            enddo
        endif

        ! Now I know how to divide it, just have to expand all the points properly.
        ! start by counting the number of unfolded points. This number is local per rank.
        local_npts=0
        do ibatch=1,grid%n_irr_batch
        do ip=1,grid%irr_batch(ibatch)%n_point
            local_npts=local_npts+grid%irr_batch(ibatch)%unfold_ctr(ip)
        enddo
        enddo

        ! Then make a little space for the local index thing that specify the point.
        ! And store information about the full set of points.
        allocate(local_ind(4,local_npts))
        local_ind=0
        l=0
        do ibatch=1,grid%n_irr_batch
        do ip=1,grid%irr_batch(ibatch)%n_point
            do i=1,grid%irr_batch(ibatch)%unfold_ctr(ip)
                l=l+1
                local_ind(1,l)=0
                local_ind(2,l)=grid%irr_batch(ibatch)%unfold_atom(i,ip)
                local_ind(3,l)=grid%irr_batch(ibatch)%unfold_index_radial(i,ip)
                local_ind(4,l)=grid%irr_batch(ibatch)%unfold_index_angular(i,ip)
            enddo
        enddo
        enddo

        ! Do the first division step, in parallel
        call parallel_minmax_predivision(p,mw,r_radial,r_angular,local_r,local_ind,local_npts,n_first_division,verbosity,n_pts_per_chunk,local_my_chunk)
        ! Divide serially
        if ( grid%ml%r .eq. 0 ) then
            call serial_minmax_division(local_r,local_ind,local_npts,&
                divisors_per_chunk(1:n_divisors_per_chunk(grid%chunk_index),grid%chunk_index),&
                set,n_set,verbosity)
        else
            ! allocate a token empty set maybe?
            allocate(set(1))
            n_set=0
        endif
        ! Distribute the sets semilocally
        call rl_distribute_sets(set,n_set,grid%ml,mw)
        ! And finally convert the sets to proper batches.
        call convert_set_into_full_batches(p,set,n_set,r_radial,r_angular,w_radial,w_angular,n_angular,grid%full_batch,grid%n_full_batch,grid%ml,mw)
    !end block full_distribution

    ! Get some statistics from the batches
    !batchstat: block
        ! get statistics for the irreducible batches
        call get_aux_information_for_batches(grid%irr_batch,p,ec,basis)
        call get_aux_information_for_batches(grid%irr_zero_batch,p,ec,basis)
        call get_aux_information_for_batches(grid%full_batch,p,ec,basis)

        ! Get some upper bounds on things to make life easier later.
        nbmax=0
        namax=0
        npmax=0
        do ib=1,grid%n_irr_batch
            ! Count maximal number of basis functions per batch?
            l=0
            do i=1,grid%irr_batch(ib)%n_relevant_ext_atom
                ie=grid%irr_batch(ib)%relevant_ext_atom(i)
                iu=ec%index_unit_cell(ie)
                is=p%species(iu)
                l=l+basis%species(is)%n_basis
            enddo
            nbmax=max(l,nbmax)
            ! and count max number of atoms per batch
            namax=max(grid%irr_batch(ib)%n_relevant_ext_atom,namax)
            ! and max number of points per batch
            npmax=max(grid%irr_batch(ib)%n_point,npmax)
        enddo
        do ib=1,grid%n_full_batch
            l=0
            do i=1,grid%full_batch(ib)%n_relevant_ext_atom
                ie=grid%full_batch(ib)%relevant_ext_atom(i)
                iu=ec%index_unit_cell(ie)
                is=p%species(iu)
                l=l+basis%species(is)%n_basis
            enddo
            nbmax=max(l,nbmax)
            namax=max(grid%full_batch(ib)%n_relevant_ext_atom,namax)
            npmax=max(grid%full_batch(ib)%n_point,npmax)
        enddo
        ! Store, locally, the max extent of several arrays.
        grid%max_n_basis_per_batch=nbmax
        grid%max_n_atom_per_batch =namax
        grid%max_n_point_per_batch=npmax

        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got batch statistics (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block batchstat

    ! Get the partition functions
    !partitionfn: block
        ! This got quite long, moved it to it's own routine.
        call evaluate_partition_function(grid%irr_batch,p,basis,free,ec,mw,&
            r_radial,r_angular,outer_partition_radius)
        call evaluate_partition_function(grid%irr_zero_batch,p,basis,free,ec,mw,&
            r_radial,r_angular,outer_partition_radius)
        call evaluate_partition_function(grid%full_batch,p,basis,free,ec,mw,&
            r_radial,r_angular,outer_partition_radius)
        if ( verbosity .gt. 0 ) then
            t1=mpi_wtime()
            write(rl_iou,*) '... got partition function (',tochar(t1-t0),'s)'
            t0=t1
        endif
    !end block partitionfn

    ! Check that I did not do anything stupid
    if ( mem%persistent_scalable .ne. 0 )    call rl_stop_gracefully(['Persistent scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%persistent_nonscalable .ne. 0 ) call rl_stop_gracefully(['Persistent nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_scalable .ne. 0 )     call rl_stop_gracefully(['Temporary scalable memory not cleared.'],rl_exitcode_memory,mw%comm)
    if ( mem%temporary_nonscalable .ne. 0 )  call rl_stop_gracefully(['Temporary nonscalable memory not cleared.'],rl_exitcode_memory,mw%comm)
end subroutine

!> In case I can not fit all integration points into memory, we have to do the grid partitioning in parallel. That is mighty annoying.
subroutine parallel_minmax_predivision(p,mw,r_radial,r_angular,local_r,local_ind,local_npts,ndiv,verbosity,n_pts_per_chunk,local_my_chunk)
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    ! Coordinates from AIMS
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    !> resulting point set (local on this rank)
    real(r8), dimension(:,:), allocatable, intent(out) :: local_r
    !> resulting indices (local on this rank)
    integer, dimension(:,:), allocatable, intent(inout) :: local_ind
    !> how many points did this rank get
    integer, intent(inout) :: local_npts
    !> how many times should I divide it?
    integer, intent(in) :: ndiv
    !> how much to talk
    integer, intent(in) :: verbosity
    !> pre-determined number of points per chunk
    integer, dimension(:), intent(in) :: n_pts_per_chunk
    !> which chunk ended up on this rank?
    integer, intent(out) :: local_my_chunk

    ! Timers
    real(r8) :: timer,t0,t1,tmr(5)
    ! Global number of points
    integer :: global_npts
    ! Number of rows in the index array
    integer :: nrowind
    ! Actual point, (3,local_npts)
    real(r8), dimension(:,:), allocatable :: r
    ! Arrays that keep track of how things are divided (local_npts)
    integer, dimension(:), allocatable :: current_div,next_div
    ! Index for looping over divisions.
    integer :: idiv
    ! Helper to keep track of things.
    type(rl_chunking_helper), dimension(:), allocatable :: ch
    integer :: nchunks
    integer, dimension(:,:), allocatable :: chunk_ctr


    !init: block
        real(r8), dimension(:,:), allocatable :: shifted_rcart
        real(r8), dimension(3) :: v0
        integer :: i,j,k,iatm,irad,iang,ispc
    !initchop: block
        real(r8), dimension(:,:,:), allocatable :: ata
        real(r8), dimension(:,:), allocatable :: dr
        real(r8), dimension(:), allocatable :: ds
        real(r8), dimension(3,3) :: eigenvectors
        real(r8), dimension(3) :: eigenvalues
        real(r8) :: f0
        integer, dimension(:,:), allocatable :: dii
        integer, dimension(:), allocatable :: ctr1
        !integer :: ip,ic,l,i
        integer :: ip,ic,l
    !findplanepos: block
        integer, parameter :: nplanepos=13,nplaneposhalf=6,maxiniter=100

        real(r8), dimension(:,:), allocatable :: planepos,rbuf
        real(r8), dimension(:), allocatable :: ct_pos,ub_pos,lb_pos
        !real(r8), dimension(3) :: v0
        !real(r8) :: f0,x1,x2,y1,y2
        real(r8) :: x1,x2,y1,y2
        integer, dimension(:,:), allocatable :: planectr
        integer, dimension(:), allocatable :: chunk_done,n_below,n_above
        !integer :: i,j,k,l,ic,ip,buf,ii1,ii2,ii3
        integer :: buf,ii1,ii2,ii3
        integer :: il,iu,jl,ju,kl,ku
        integer :: iter
    !slicechunks: block
        !real(r8) :: f0
        integer, dimension(:,:,:), allocatable :: octr1,octr2
        integer, dimension(:,:), allocatable :: mctr1,mctr2
        !integer, dimension(:), allocatable :: n_below,n_above,di
        integer, dimension(:), allocatable :: di
        !integer :: ic,ip,ir,ii,jj,kk,ll,i
        integer :: ir,ii,jj,kk,ll
    !distributepoints: block
        real(r8), dimension(:,:), allocatable :: rbuf0,rbuf1
        integer, dimension(:,:), allocatable :: ibuf0,ibuf1

        integer, dimension(:,:), allocatable :: ctr2
        integer, dimension(:), allocatable :: sendcounts,sendoffset
        integer, dimension(:), allocatable :: recvcounts,recvoffset
        !integer, dimension(:), allocatable :: di,dj
        integer, dimension(:), allocatable :: dj
        !integer :: ic,ip,i,j
        integer :: my_chunk_size,my_chunk


    ! Figure out some basic things and distribute the points across ranks.
    !init: block

        timer=mpi_wtime()
        t0=timer
        t1=0.0_r8
        tmr=0.0_r8

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) ''
            write(rl_iou,*) 'DISTRIBUTED MESH DIVISION'
        endif

        ! I will likely need shifted Cartesian coordinates
        allocate(shifted_rcart(3,p%n_atom))
        do i=1,p%n_atom
            v0=rl_clean_fractional_coordinates(p%fractional_coordinate(:,i)+0.5_r8)-0.5_r8
            shifted_rcart(:,i)=matmul(p%latticevectors,v0)
        enddo

        ! Get the global number of points
        call mw%allreduce('sum',local_npts,global_npts)
        ! Make sure we agree on at least some basic things.
        if ( sum(n_pts_per_chunk) .ne. global_npts ) then
            call rl_stop_gracefully(['Inconsistent number of points'],rl_exitcode_param,mw%comm)
        endif
        if ( 2**ndiv .ne. size(n_pts_per_chunk) ) then
            call rl_stop_gracefully(['Inconsistent number of chunks'],rl_exitcode_param,mw%comm)
        endif

        ! Figure out the number of rows in the index array
        if ( local_npts .gt. 0 ) then
            nrowind=size(local_ind,1)
        else
            nrowind=0
        endif
        call mw%allreduce('max',nrowind)

        ! Make space for the actual points and some helper arrays
        if ( local_npts .gt. 0 ) then
            allocate(r(3,local_npts))
            do i=1,local_npts
                iatm=local_ind(2,i)
                ispc=p%species(iatm)
                irad=local_ind(3,i)
                iang=local_ind(4,i)
                v0=shifted_rcart(:,iatm)+r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)
                v0=matmul(p%inv_latticevectors,v0)
                v0=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
                v0=matmul(p%latticevectors,v0)
                r(:,i)=v0
            enddo
        else
            allocate(r(1,1))
            r=-rl_huge
        endif

        ! Initialize the chunk divisions
        if ( local_npts .gt. 0 ) then
            allocate(current_div(local_npts))
            allocate(next_div(local_npts))
            current_div=1
            next_div=0
        else
            allocate(current_div(1))
            allocate(next_div(1))
            current_div=-rl_hugeint
            next_div=-rl_hugeint
        endif

        ! It is a neat idea to determine in advance how many points the set should be divided
        ! into at each iteration.
        if ( ndiv .gt. 0 ) then
            allocate(chunk_ctr(2**ndiv,ndiv+1))
            chunk_ctr=-1
            ! I start from the bottom:
            chunk_ctr(:,ndiv+1)=n_pts_per_chunk
            do i=ndiv,1,-1
                k=2**(i-1)
                do j=1,k
                    chunk_ctr(j,i)=chunk_ctr(j,i+1)+chunk_ctr(j+k,i+1)
                enddo
            enddo
        else
            ! No division.
            allocate(chunk_ctr(1,1))
            chunk_ctr=global_npts
        endif

        ! Accumulate timer
        t1=mpi_wtime()
        tmr(1)=tmr(1)+t1-t0
        t0=t1

        ! Start the counter
        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) '    divided mesh, iteration 00  ( pts/chunk: ',tochar(global_npts),' )'
        endif
    !end block init

    ! Start building the dividing planes. This gets a little messy, but no worries.
    ! I am extra careful because I really want the tree to be balanced in the end.
    divloop: do idiv=1,ndiv

    ! At this point, I have a division, or a way to divide the points into even
    ! sized pieces. It is stored in current_div(local_npts). For each of these
    ! divisions, how many they might be, I want to slice each of them into two
    ! new chunks. First I need a sensible guess for where to put the slicing plane,
    ! as well as some statistics about each chunk that I can later use.
    !initchop: block

        ! How many chunks do we have.
        nchunks=2**(idiv-1)

        ! Make space for the helper
        if ( allocated(ch) ) deallocate(ch)
        allocate(ch(nchunks))

        allocate(ctr1(nchunks))
        allocate(dr(3,nchunks))
        ctr1=0
        dr=0.0_r8
        ! Calculate the center of mass per chunk.
        do ip=1,local_npts
            ic=current_div(ip)
            dr(:,ic)=dr(:,ic)+r(:,ip)
            ctr1(ic)=ctr1(ic)+1
        enddo
        ! Store the local count
        do ic=1,nchunks
            ch(ic)%n_local=ctr1(ic)
        enddo
        ! Accumulate for the average
        call mw%allreduce('sum',dr)
        call mw%allreduce('sum',ctr1)
        ! Store the center of mass.
        do ic=1,nchunks
            ch(ic)%n_global=ctr1(ic)
            ch(ic)%com = dr(:,ic)/real(ch(ic)%n_global,r8)
        enddo
        deallocate(dr)

        ! Now we have the center of mass per chunk. Next step is to calculate the
        ! covariance matrix so that I can decide on where to put the normal for the
        ! plane that should divide each chunk of points.

        ! Get local indices per chunk
        l=maxval(ch(:)%n_local)
        allocate(dii(l,nchunks))
        allocate(dr(3,l))
        dr=0.0_r8
        ctr1=0
        dii=-1
        do ip=1,local_npts
            ic=current_div(ip)
            ctr1(ic)=ctr1(ic)+1
            dii(ctr1(ic),ic)=ip
        enddo

        ! get covariance matrix to get the dividing planes
        allocate(ds(nchunks))
        allocate(ata(3,3,nchunks))
        ds=0.0_r8
        ata=0.0_r8
        ! With the center of mass, get the normals of the dividing planes
        do ic=1,nchunks
            dr=0.0_r8
            do i=1,ctr1(ic)
                ip=dii(i,ic)
                dr(:,i)=r(:,ip)-ch(ic)%com
                ! Also keep track of a bounding sphere per chunk, will come handy later
                ds(ic)=max(ds(ic),rl_sqnorm(dr(:,i)))
            enddo
            if ( ctr1(ic) .gt. 0 ) then
                ! This accumulates the covariance matrix
                call dsyrk('U','N',3,ctr1(ic),1.0_r8,dr(:,1:ctr1(ic)),3,0.0_r8,ata(:,:,ic),3)
                ata(2,1,ic)=ata(1,2,ic)
                ata(3,1,ic)=ata(1,3,ic)
                ata(3,2,ic)=ata(2,3,ic)
            else
                ata(:,:,ic)=0.0_r8
            endif
        enddo
        call mw%allreduce('sum',ata)
        call mw%allreduce('max',ds)

        ! Calculate normals
        deallocate(dr)
        allocate(dr(3,nchunks))
        dr=0.0_r8
        do ic=1,nchunks
            ! Pick the normal to be the one with the largest eigenvalue of the covariance matrix
            call rl_symmetric_eigensystem_3x3matrix(ata(:,:,ic),eigenvalues,eigenvectors)
            f0=-rl_huge
            do i=1,3
                if ( eigenvalues(i) .gt. f0 ) then
                    f0=eigenvalues(i)
                    dr(:,ic)=eigenvectors(:,i)/norm2(eigenvectors(:,i))
                endif
            enddo
        enddo

        ! Make sure everyone has the same normals, will get really weird
        ! if there are som bit shifts or something.
        call mw%check_and_sync(dr,0,bigdeal=1,vname='normals')
        do ic=1,nchunks
            ! Store the normals
            ch(ic)%normal=dr(:,ic)
            ! Also store the bounding sphere
            ch(ic)%boundingsphere=sqrt(ds(ic))
        enddo

        ! Final cleanup.
        deallocate(ds)
        deallocate(dr)
        deallocate(ctr1)
        deallocate(dii)
        deallocate(ata)

        ! Accumulate the timer
        t1=mpi_wtime()
        tmr(2)=tmr(2)+t1-t0
        t0=t1
    !end block initchop

    ! So, where are we now. We are in the middle of trying to divide nchunks set of
    ! points into 2*nchunks set of points. We have the centers of mass of the old
    ! chunks, we have the normals of the dividing planes, but we don't quite know
    ! exactly where to put the planes so that we divide the point set into two
    ! equal halves. Normally this is done via sorting by the distance from the plane,
    ! but parallel sorting is quite annoying and requires a lot of communication, so
    ! I solve it with a bisection algorithm instead. I think it works because I can
    ! get a decent guess for where the dividing plane should be, and work from there.
    !findplanepos: block
        ! Parameters, to be checked what makes sense later.

        ! I know, from the recursion thingy in the init block, how many
        ! points I want above and below the plane for each chunk. I store
        ! that here to make things easier.
        allocate(n_below(nchunks))
        allocate(n_above(nchunks))
        do ic=1,nchunks
            n_below(ic)=chunk_ctr(ic,idiv+1)
            n_above(ic)=chunk_ctr(ic+nchunks,idiv+1)
        enddo

        ! First we have to guess where to put the plane:
        allocate(ct_pos(nchunks))
        allocate(ub_pos(nchunks))
        allocate(lb_pos(nchunks))
        do ic=1,nchunks
            f0=-dot_product(ch(ic)%com,ch(ic)%normal)
            ! Guess for where the exact division might be
            ct_pos(ic)=f0
            ! Safe upper and lower bounds
            ub_pos(ic)=f0+(ch(ic)%boundingsphere+rl_tol)
            lb_pos(ic)=f0-(ch(ic)%boundingsphere+rl_tol)
            ! Stow away safe bounds
            ch(ic)%pos_lb=ub_pos(ic)
            ch(ic)%pos_ub=lb_pos(ic)
        enddo

        ! Space for bisection algo
        allocate(planepos(nplanepos,nchunks))
        allocate(planectr(nplanepos,nchunks))
        allocate(chunk_done(nchunks))
        planepos=0.0_r8
        planectr=0
        chunk_done=0
        ! Start the bisection algorithm
        itrl1: do iter=1,maxiniter

            ! First thing is to choose good positions for the planes.
            do ic=1,nchunks
                ! Don't bother with the chunks that are already done
                if ( abs(chunk_done(ic)) .eq. 1 ) cycle
                ! Space the points from the guess to the lower and upper bounds.
                do k=-nplaneposhalf,nplaneposhalf
                    i=nplaneposhalf+k+1
                    f0=(abs(k))/real(nplaneposhalf,r8)
                    ! I space them logarithmically the first iteration to speed things up.
                    ! Actually I don't. Nevermind. Did not help particularly.
                    !if ( iter .eq. 1 ) then
                    !    f0=f0**planeposexponent
                    !endif
                    if ( k .lt. 0 ) then
                        planepos(i,ic)=ct_pos(ic)+f0*(ub_pos(ic)-ct_pos(ic))
                    elseif (k .gt. 0 ) then
                        planepos(i,ic)=ct_pos(ic)+f0*(lb_pos(ic)-ct_pos(ic))
                    else
                        planepos(i,ic)=ct_pos(ic)
                    endif
                enddo
            enddo

            ! For each of these positions, count the number of points below the plane.
            planectr=0
            do ip=1,local_npts
                ic=current_div(ip)
                if ( abs(chunk_done(ic)) .eq. 1 ) cycle
                v0=r(:,ip)
                do i=1,nplanepos
                    f0=dot_product(v0,ch(ic)%normal)+planepos(i,ic)
                    if ( f0 .lt. 0.0_r8 ) then
                        planectr(i,ic)=planectr(i,ic)+1
                    endif
                enddo
            enddo
            call mw%allreduce('sum',planectr)

            ! Good place to do a sanity check, my bounds have to work properly.
            if ( iter .eq. 1 ) then
                do ic=1,nchunks
                    if ( minval(planectr(:,ic)) .ne. 0 ) then
                        call rl_stop_gracefully(['Did not bound the set of points, should be impossible'],rl_exitcode_symmetry,mw%comm)
                    endif
                    if ( maxval(planectr(:,ic)) .ne. ch(ic)%n_global ) then
                        call rl_stop_gracefully(['Did not bound the set of points, should be impossible'],rl_exitcode_symmetry,mw%comm)
                    endif
                enddo
            endif

            ! Now evaluate how this went
            chl1: do ic=1,nchunks
                if ( abs(chunk_done(ic)) .eq. 1 ) cycle
                ! Go through the different plane positions and check how we are doing.
                ii1=rl_hugeint
                ii2=-rl_hugeint
                ii3=rl_hugeint
                l=0
                il=0
                iu=0
                do i=1,nplanepos
                    j=planectr(i,ic)-n_below(ic)
                    ! This masks out the best one, keep it for now!
                    if ( abs(j) .lt. ii1 ) then
                        ii1=abs(j)
                        ch(ic)%pos=planepos(i,ic)
                        l=i
                    endif
                    ! This should mask out the lower bound
                    if ( j .le. 0 .and. j .gt. ii2 ) then
                        ii2=j
                        il=i
                    endif
                    ! This should mask out the upper bound
                    if ( j .gt. 0 .and. j .lt. ii3 ) then
                        ii3=j
                        iu=i
                    endif
                enddo

                ! Small sanity test to catch. Should be impossible to trigger.
                if ( il .eq. 0 ) then
                    call rl_stop_gracefully(['Could not determine lower bound'],rl_exitcode_symmetry,mw%comm)
                endif
                if ( iu .eq. 0 ) then
                    call rl_stop_gracefully(['Could not determine upper bound'],rl_exitcode_symmetry,mw%comm)
                endif

                ! I want safe upper and lower bounds such that I eventually
                ! can get a balanced tree, which makes life easier.
                jl=abs(planectr(il,ic)-n_below(ic))
                ju=abs(planectr(iu,ic)-n_below(ic))
                buf=max(jl,ju) ! I need at least this margin in both directions.

                kl=0
                ku=0
                ii2=-rl_hugeint
                ii3=rl_hugeint
                do i=1,nplanepos
                    j=planectr(i,ic)-n_below(ic)
                    ! This should mask out the lower bound
                    if ( j .lt. 0 .and. j .gt. ii2 .and. abs(j) .ge. buf ) then
                        ii2=j
                        kl=i
                    endif
                    ! This should mask out the upper bound
                    if ( j .gt. 0 .and. j .lt. ii3 .and. abs(j) .ge. buf ) then
                        ii3=j
                        ku=i
                    endif
                enddo

                ! This should trigger quite often, but not always. Hopefully at least once.
                ! I will have to rethink my life in case this fails at some point. Pretty
                ! sure this should be safe. Almost always.
                if ( kl .gt. 0 .and. ku .gt. 0 ) then
                    ch(ic)%pos_lb=planepos(kl,ic)
                    ch(ic)%pos_ub=planepos(ku,ic)
                endif

                ! ! Check if the division is really good!
                if ( planectr(l,ic)-n_below(ic) .eq. 0 ) then
                    ! Yep, this is a good dividing plane.
                    chunk_done(ic)=1
                else
                    ! It is not good enough.
                    lb_pos(ic)=planepos(il,ic)
                    ub_pos(ic)=planepos(iu,ic)
                    ! Linearly interpolate between these to points to get a guess for the new center.
                    x1=lb_pos(ic)
                    x2=ub_pos(ic)
                    y1=real(planectr(il,ic)-n_below(ic),r8)
                    y2=real(planectr(iu,ic)-n_below(ic),r8)
                    if ( abs(y1-y2) .gt. 1E-11_r8 ) then
                        f0=y1/(y1-y2)
                        f0=max(f0,1E-8_r8)
                        f0=min(f0,1.0_r8-1E-8_r8)
                    else
                        f0=0.5_r8
                    endif
                    ct_pos(ic)=x1+f0*(x2-x1)
                endif

                if ( chunk_done(ic) .eq. 0 ) then
                    ! Also, this chunk could have gotten stuck, with no good place
                    ! to put the dividing plane. This will be solved, eventually.
                    ! So I count the number of different divisions I can get.
                    l=0
                    do i=1,nplanepos-1
                        if ( abs(planectr(i,ic)-planectr(i+1,ic)) .ne. 0 ) l=l+1
                    enddo
                    if ( l .eq. 1 ) then
                        ! Yep, this chunk is stuck. Deal with it later.
                        chunk_done(ic)=-1
                    endif
                endif
            enddo chl1

            ! Check for convergence?
            if ( sum(abs(chunk_done)) .eq. nchunks ) then
                exit itrl1
            endif
        enddo itrl1

        ! So, once the iterative solution is done, make sure everyone
        ! agrees on everything. Just want to be really sure.
        allocate(rbuf(3,nchunks))
        do ic=1,nchunks
            rbuf(:,ic)=[ch(ic)%pos,ch(ic)%pos_lb,ch(ic)%pos_ub]
        enddo
        call mw%check_and_sync(chunk_done,0,bigdeal=1,vname='chunk_done')
        call mw%check_and_sync(rbuf,0,bigdeal=1,vname='rbuf')
        ! Small sanity check, all chunks need a verdit.
        if ( count(chunk_done==0) .gt. 0 ) then
            call rl_stop_gracefully(['No verdict on this chunk, should be impossible'],rl_exitcode_symmetry,mw%comm)
        endif

        ! Seems it went fine.
        do ic=1,nchunks
            ch(ic)%howslice = chunk_done(ic)
            ch(ic)%pos = rbuf(1,ic)
            ch(ic)%pos_lb = rbuf(2,ic)
            ch(ic)%pos_ub = rbuf(3,ic)
        enddo

        ! And some cleanup
        deallocate(rbuf)
        deallocate(planepos)
        deallocate(ct_pos)
        deallocate(ub_pos)
        deallocate(lb_pos)
        deallocate(planectr)
        deallocate(chunk_done)
        deallocate(n_below)
        deallocate(n_above)


        ! Accumulate the timer
        t1=mpi_wtime()
        tmr(3)=tmr(3)+t1-t0
        t0=t1
    !end block findplanepos

    ! Now that I have the planes, I just have to divide the chunks into two
    ! new chunks. The only reason this is so awfully confusing is that positioning
    ! the dividing plane fails sometimes, which means I have a set of points stuck
    ! in limbo, that I have to assign. It happens when the dividing plane goes exactly
    ! through many points. Then it does not really matter where you put the points that
    ! are undecided, so I just make sure to divide them evenly.
    ! I think I now have a sensible algorithm to do so.
    !slicechunks: block

        ! I know nothing about the next division.
        next_div=0
        ! Except I know something, I know how many should end up above and below the division.
        allocate(n_below(nchunks))
        allocate(n_above(nchunks))
        do ic=1,nchunks
            n_below(ic)=chunk_ctr(ic,idiv+1)
            n_above(ic)=chunk_ctr(ic+nchunks,idiv+1)
        enddo

        ! Initialize some counters.
        allocate(octr1(3,nchunks,mw%n))
        octr1=0     ! Big counter. Holds ([nup,ndown,nmiss],chunk,rank)

        ! Now we are hopefully somewhat converged. Do the actual division.
        do ip=1,local_npts
            ic=current_div(ip)
            if ( ch(ic)%howslice .eq. 1 ) then
                ! If the chunk finished cleanly, easy to divide
                f0=dot_product(r(:,ip),ch(ic)%normal)+ch(ic)%pos
                if ( f0 .lt. 0.0_r8 ) then
                    next_div(ip)=ic
                    octr1(1,ic,mw%r+1)=octr1(1,ic,mw%r+1)+1
                else
                    next_div(ip)=ic+nchunks
                    octr1(2,ic,mw%r+1)=octr1(2,ic,mw%r+1)+1
                endif
            elseif ( ch(ic)%howslice .eq. -1 ) then
                ! This could not be divided easily, have to think.
                f0=dot_product(r(:,ip),ch(ic)%normal)+ch(ic)%pos_lb
                if ( f0 .lt. 0.0_r8 ) then
                    next_div(ip)=ic
                    octr1(1,ic,mw%r+1)=octr1(1,ic,mw%r+1)+1
                endif
                f0=dot_product(r(:,ip),ch(ic)%normal)+ch(ic)%pos_ub
                if ( f0 .gt. 0.0_r8 ) then
                    next_div(ip)=ic+nchunks
                    octr1(2,ic,mw%r+1)=octr1(2,ic,mw%r+1)+1
                endif
            endif
            ! Count how many points did not get assigned, rank-resolved.
            if ( next_div(ip) .eq. 0 ) then
                octr1(3,ic,mw%r+1)=octr1(3,ic,mw%r+1)+1
            endif
        enddo

        ! Sync up counters to be able to decide if anything annoying will have to be dealt with.
        call mw%allreduce('sum',octr1)

        ! Number of points that got assigned
        ii=sum(octr1(1:2,:,:))
        ! Number of points not assigned
        jj=sum(octr1(3,:,:))
        ! Always good with a sanity check
        if ( ii+jj .ne. global_npts ) then
            call rl_stop_gracefully(['Clearly I can not count'],rl_exitcode_symmetry,mw%comm)
        endif

        ! Now se if everything got assigned. If not, deal with it.
        ! As I mentioned above, sometimes the division with planes does
        ! leave some points in limbo, the code block below is what deals
        ! with the points in limbo.
        if ( jj .ne. 0 ) then
            ! First, count how many points should be assigned up, and
            ! how many should be assigned down.
            allocate(mctr1(2,nchunks))
            allocate(mctr2(2,nchunks))
            ! Counter for the unassigned points
            mctr1=0
            mctr2=0
            ! Global counter for how points was assigned.
            do ic=1,nchunks
            do i=1,2
                mctr2(i,ic)=sum(octr1(i,ic,:))
            enddo
            enddo
            ! Go through all points and count how they should be assigned.
            do ic=1,nchunks
                ! Global number of unassigned for this chunk.
                ii=sum(octr1(3,ic,:))
                ! Go through and count how it should be assigned.
                do i=1,ii
                    if ( mctr2(1,ic) .lt. n_below(ic) ) then
                        mctr1(1,ic)=mctr1(1,ic)+1
                        mctr2(1,ic)=mctr2(1,ic)+1
                    else
                        mctr1(2,ic)=mctr1(2,ic)+1
                        mctr2(2,ic)=mctr2(2,ic)+1
                    endif
                enddo
            enddo

            ! So I know, globally, how many above and below I should assign.
            ! This has to be spread out over ranks, to keep the global tree
            ! completely balanced, or at least pretty well balanced.
            allocate(octr2(2,nchunks,mw%n))
            octr2=0
            do ic=1,nchunks
                ! For this chunk, how many points should get assigned up or down?
                ! This is held in mctr1(1,j) for down, mctr1(2,j) for up.
                ! So I will increment these counters as I assign, and compare
                ! at the end.
                ii=0
                jj=0
                do ir=1,mw%n
                    ! On this rank, how many unassigned are there?
                    do i=1,octr1(3,ic,ir)
                        if ( ii .lt. mctr1(1,ic) ) then
                            ! Assign this point below
                            ii=ii+1
                            octr2(1,ic,ir)=octr2(1,ic,ir)+1
                        else
                            ! Assign it above
                            jj=jj+1
                            octr2(2,ic,ir)=octr2(2,ic,ir)+1
                        endif
                    enddo
                enddo

                ! Sanity check this, that all points got assigned.
                if ( ii .ne. mctr1(1,ic) ) then
                    call rl_stop_gracefully(['Clearly I can not count'],rl_exitcode_symmetry,mw%comm)
                endif
                if ( jj .ne. mctr1(2,ic) ) then
                    call rl_stop_gracefully(['Clearly I can not count'],rl_exitcode_symmetry,mw%comm)
                endif

                ! Stronger sanity check, I think
                ii=sum(octr1(1,ic,:))
                jj=sum(octr1(2,ic,:))
                kk=sum(octr2(1,ic,:))
                ll=sum(octr2(2,ic,:))
                if ( ii+jj+kk+ll .ne. ch(ic)%n_global ) then
                    call rl_stop_gracefully(['Clearly I can not count'],rl_exitcode_symmetry,mw%comm)
                endif
            enddo

            ! And finally, I can start assigning.
            mctr2=0
            do ip=1,local_npts
                if ( next_div(ip) .ne. 0 ) cycle
                ic=current_div(ip)
                if ( mctr2(1,ic) .lt. octr2(1,ic,mw%r+1) ) then
                    next_div(ip)=ic
                    mctr2(1,ic)=mctr2(1,ic)+1
                else
                    next_div(ip)=ic+nchunks
                    mctr2(2,ic)=mctr2(2,ic)+1
                endif
            enddo

            ! Some cleanup
            deallocate(mctr1)
            deallocate(mctr2)
            deallocate(octr2)
        endif

        ! And update the divisions for the next iteration
        current_div=next_div
        next_div=0
        ! Cleanup
        deallocate(octr1)
        ! Accumulate timers
        t1=mpi_wtime()
        tmr(4)=tmr(4)+t1-t0
        t0=t1

        ! Maybe report how it's going?
        if ( verbosity .gt. 0 ) then
            allocate(di(size(n_below)*2))
            ll=0
            do ic=1,size(n_below)
                ll=ll+1
                di(ll)=n_below(ic)
            enddo
            do ic=1,size(n_above)
                ll=ll+1
                di(ll)=n_above(ic)
            enddo
            if ( ll .le. 8 ) then
            write(rl_iou,*) '    divided mesh, iteration ',tochar(idiv,2),'  ( pts/chunk: ',tochar(di),' )'
            else
            write(rl_iou,*) '    divided mesh, iteration ',tochar(idiv,2),'  ( pts/chunk: ',tochar(di(1:8)),' ... )'
            endif
            if ( idiv .eq. ndiv ) then
            write(rl_iou,*) '    =============================================================='
            endif
            deallocate(di)
        endif
        deallocate(n_below)
        deallocate(n_above)
    !end block slicechunks

    enddo divloop

    ! Now the grids are divided into manageable chunks. Now I have to
    ! distribute them across ranks, in some manner. Not sure how to do
    ! that, but here is a good start, I suppose?
    !distributepoints: block

        ! So, how many divisions are there?
        nchunks=maxval(current_div)
        call mw%allreduce('max',nchunks)

        ! Sanity check my algorithm.
        if ( nchunks .gt. mw%n ) then
            call rl_stop_gracefully(['More chunks than MPI ranks, that should never happen.'],rl_exitcode_symmetry,mw%comm)
        endif

        ! This is a little annoying, but I have to sort the points
        ! according to chunk, otherwise all-to-all will be way too complicated.
        ! At the end, rbuf0 will hold the points, ibuf0 the indices, sorted by chunk.
        if ( local_npts .gt. 0 ) then
            allocate(di(local_npts))
            allocate(dj(local_npts))
            allocate(rbuf0(3,local_npts))
            allocate(ibuf0(nrowind+1,local_npts))
        else
            allocate(di(1))
            allocate(dj(1))
            allocate(rbuf0(1,1))
            allocate(ibuf0(1,1))
        endif
        di=0
        dj=0
        rbuf0=0.0_r8
        ibuf0=0
        di=current_div
        call rl_qsort(di,dj)
        do i=1,local_npts
            j=dj(i)
            rbuf0(:,i)=r(:,j)
            ibuf0(1:nrowind,i)=local_ind(:,j)
            ibuf0(nrowind+1,i)=current_div(j)
        enddo

        allocate(sendcounts(mw%n))
        allocate(recvcounts(mw%n))
        allocate(sendoffset(mw%n))
        allocate(recvoffset(mw%n))
        sendcounts=0
        recvcounts=0
        sendoffset=0
        recvoffset=0

        ! Do some counting to help with distribution.
        allocate(ctr2(mw%n,nchunks))
        ctr2=0
        do ip=1,local_npts
            ic=di(ip)
            ctr2(mw%r+1,ic)=ctr2(mw%r+1,ic)+1
        enddo
        call mw%allreduce('sum',ctr2)

        ! How many points should this rank recieve?
        if ( mw%r+1 .gt. nchunks ) then
            my_chunk_size=0
            my_chunk=-1
        else
            my_chunk_size=sum(ctr2(:,mw%r+1))
            my_chunk=mw%r+1
        endif

        ! count numbers to send
        j=0
        do i=1,nchunks
            sendcounts(i)=ctr2(mw%r+1,i)
            sendoffset(i)=j
            j=j+ctr2(mw%r+1,i)
        enddo
        ! count numbers to recieve
        if ( my_chunk_size .gt. 0 ) then
            j=0
            do i=1,mw%n
                recvcounts(i)=ctr2(i,mw%r+1)
                recvoffset(i)=j
                j=j+ctr2(i,mw%r+1)
            enddo
        else
            recvcounts=0
            recvoffset=0
        endif

        ! Prepare output!
        if ( my_chunk_size .gt. 0 ) then
            allocate(rbuf1(3,my_chunk_size))
            allocate(ibuf1(nrowind+1,my_chunk_size))
            rbuf1=0.0_r8
            ibuf1=0
        else
            ! allocate token things? Maybe
            allocate(rbuf1(1,1))
            allocate(ibuf1(1,1))
            rbuf1=rl_huge
            ibuf1=rl_hugeint
        endif

        ! Send everything everywhere. Copied the signature here, always forget how it works.
        !
        ! MPI_ALLTOALLV(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, ierr)
        ! IN sendbuf	starting address of send buffer (choice)
        ! IN sendcounts	integer array equal to the group size specifying the number of elements to send to each processor
        ! IN sdispls	integer array (of length group size). Entry j specifies the displacement (relative to sendbuf from which to take the outgoing data destined for process j
        ! IN sendtype	data type of send buffer elements (handle)
        ! OUT recvbuf	address of receive buffer (choice)
        ! IN recvcounts	integer array equal to the group size specifying the number of elements that can be received from each processor
        ! IN rdispls	integer array (of length group size). Entry i specifies the displacement (relative to recvbuf at which to place the incoming data from process i
        ! IN recvtype	data type of receive buffer elements (handle)
        ! IN comm	communicator (handle)
        !

        ! Check that the sizes make sense, so that it dies with an easier error message than the cryptic MPI ones.
        if ( sum(sendcounts) .ne. size(ibuf0,2) ) then
            call rl_stop_gracefully(['Inconsistent sizes of arrays when communicating'],rl_exitcode_symmetry,mw%comm)
        endif
        if ( sum(sendcounts) .ne. size(rbuf0,2) ) then
            call rl_stop_gracefully(['Inconsistent sizes of arrays when communicating'],rl_exitcode_symmetry,mw%comm)
        endif
        if ( sum(recvcounts) .ne. size(ibuf1,2) .and. my_chunk_size .gt. 0 ) then
            call rl_stop_gracefully(['Inconsistent sizes of arrays when communicating'],rl_exitcode_symmetry,mw%comm)
        endif
        if ( sum(recvcounts) .ne. size(rbuf1,2) .and. my_chunk_size .gt. 0 ) then
            call rl_stop_gracefully(['Inconsistent sizes of arrays when communicating'],rl_exitcode_symmetry,mw%comm)
        endif

        ! Begin operation actual operation.
        call MPI_ALLTOALLV(rbuf0, sendcounts*3, sendoffset*3, MPI_DOUBLE_PRECISION, rbuf1, recvcounts*3, recvoffset*3, MPI_DOUBLE_PRECISION, mw%comm, mw%error)
        if ( mw%error .ne. 0 ) then
            call rl_stop_gracefully(['Failed distributing points.'],rl_exitcode_symmetry,mw%comm)
        endif
        call MPI_ALLTOALLV(ibuf0, sendcounts*(nrowind+1), sendoffset*(nrowind+1), MPI_INTEGER, ibuf1, recvcounts*(nrowind+1), recvoffset*(nrowind+1), MPI_INTEGER, mw%comm, mw%error)
        if ( mw%error .ne. 0 ) then
            call rl_stop_gracefully(['Failed distributing points.'],rl_exitcode_symmetry,mw%comm)
        endif

        ! Another sanity test, just to make sure I divided things right! Naturally, this will never trigger.
        do i=1,my_chunk_size
            if ( ibuf1(nrowind+1,i) .ne. my_chunk ) then
                call rl_stop_gracefully(['Failed distributing points.'],rl_exitcode_symmetry,mw%comm)
            endif
        enddo

        ! Some pre-cleanup?
        deallocate(r)
        deallocate(current_div)
        deallocate(next_div)
        if ( allocated(ch) ) deallocate(ch)
        deallocate(rbuf0)
        deallocate(ibuf0)
        deallocate(sendcounts)
        deallocate(recvcounts)
        deallocate(sendoffset)
        deallocate(recvoffset)

        ! ! Now I can actually output clean arrays with points and indices!
        local_npts=my_chunk_size
        if ( local_npts .gt. 0 ) then
            ! This means this ranks has some points
            deallocate(local_ind)
            allocate(local_r(3,local_npts))
            allocate(local_ind(nrowind,local_npts))
            local_r = rbuf1
            do i=1,local_npts
                local_ind(:,i)=ibuf1(1:nrowind,i)
            enddo
            local_my_chunk=my_chunk
        else
            ! This rank has no points, but allocate minimal anyway.
            deallocate(local_ind)
            allocate(local_r(1,1))
            allocate(local_ind(1,1))
            local_r=-rl_huge
            local_ind=-rl_hugeint
            local_my_chunk=-1
        endif

        ! Accumulate timers
        t1=mpi_wtime()
        tmr(5)=tmr(5)+t1-t0
        t0=t1
    !end block distributepoints

    ! Short timing report maybe?
    ! if ( verbosity .gt. 0 ) then
    !     write(rl_iou,*) 'Timers for parallel grid division:'
    !     write(rl_iou,"(1X,A,F12.8,F10.5)") '          init:',tmr(1),100.0_r8*tmr(1)/sum(tmr)
    !     write(rl_iou,"(1X,A,F12.8,F10.5)") '    com/normal:',tmr(2),100.0_r8*tmr(2)/sum(tmr)
    !     write(rl_iou,"(1X,A,F12.8,F10.5)") '      position:',tmr(3),100.0_r8*tmr(3)/sum(tmr)
    !     write(rl_iou,"(1X,A,F12.8,F10.5)") '         slice:',tmr(4),100.0_r8*tmr(4)/sum(tmr)
    !     write(rl_iou,"(1X,A,F12.8,F10.5)") '    final comm:',tmr(5),100.0_r8*tmr(5)/sum(tmr)
    !     write(rl_iou,"(1X,A,F12.8)")       '         total:',sum(tmr)
    ! endif
end subroutine

!> take an arbitrary list of points and slice it into batches according to a set of divisors.
subroutine serial_minmax_division(r,ind,npts,divisors,set,nset,verbosity)
    !> list of points
    real(r8), dimension(:,:), allocatable, intent(inout) :: r
    !> indices to keep track where the points came from
    integer, dimension(:,:), intent(in) :: ind
    !> how many points
    integer, intent(in) :: npts
    !> which numbers should I divide with
    integer, dimension(:), intent(in) :: divisors
    !> resulting division of points
    type(rl_set_of_points), dimension(:), allocatable, intent(out) :: set
    !> how many sets did I get?
    integer, intent(out) :: nset
    !> talk a lot?
    integer, intent(in) :: verbosity

    integer, dimension(:,:), allocatable :: setctr

    !init: block
        integer, dimension(:), allocatable :: di,dj,dk
        integer :: i,j,k,ll,nval,nmult
    !divide: block
        integer, dimension(:), allocatable :: splitcount
        integer :: mult
        !integer :: i,j,is,iter
        integer :: is,iter

    ! Set some quick things
    !init: block

        ! Now I know how many sets I should get. Then I can figure out the number of points
        ! per set, and then recursively figure out how the points should be divided to get
        ! an equal number per batch.
        nset=product(divisors)
        ll=size(divisors)+1
        allocate(setctr(nset,ll))
        allocate(di(nset))
        allocate(dj(ll))
        allocate(dk(ll))

        ! Number of points per set
        di=0
        do i=1,npts
            j=mod(i,nset)+1
            di(j)=di(j)+1
        enddo
        ! Initialize the counter
        setctr=0
        setctr(:,ll)=di
        ! dummy thingy for the divisors
        dj(1)=1
        dj(2:ll)=divisors
        ! dummy thing with the number of sets per iteration
        dk=0
        j=1
        do i=1,ll
            j=j*dj(i)
            dk(i)=j
        enddo
        ! figure out the number of points per set backwards, from the final distribution.
        do i=ll-1,1,-1
            nmult=dj(i+1)
            nval=dk(i)
            do j=1,nval
            do k=1,nmult
                setctr(j,i)=setctr(j,i)+setctr(j+(k-1)*nval,i+1)
            enddo
            enddo
        enddo
        ! And cleanup the dummys
        deallocate(di)
        deallocate(dj)
        deallocate(dk)

        ! Space for the sets of points, initialized to empty.
        allocate(set(nset))
        do i=1,nset
            set(i)%n=0
        enddo

        ! In the beginning, we only have one set of points, so fill it up with the input
        allocate(set(1)%r(3,npts))
        allocate(set(1)%ind(size(ind,1),npts))
        set(1)%n=npts
        set(1)%r=r
        set(1)%ind=ind
        nset=1

        if ( verbosity .gt. 0 ) then
            write(rl_iou,*) ''
            write(rl_iou,*) 'SERIAL DIVISION OF MESH'
        endif
    !end block init

    ! Now start dividing the set of points
    !divide: block

        iterloop: do iter=1,size(divisors)
            ! How many times should each set be split?
            mult=divisors(iter)
            allocate(splitcount(mult))

            do is=1,nset
                ! Fetch how it should be split
                do i=1,mult
                    j=is+(i-1)*nset
                    splitcount(i)=setctr(j,iter+1)
                enddo
                ! Split it!
                call rl_split_set(set,is,nset,splitcount)
            enddo
            deallocate(splitcount)

            ! Update the number of sets!
            nset=nset*mult

            ! Maybe report how it's going?
            if ( verbosity .gt. 0 ) then
                allocate(splitcount(nset))
                splitcount=set(1:nset)%n
                if ( nset .le. 8 ) then
                    write(rl_iou,*) '    divided mesh, iteration ',tochar(iter,2),'  ( pts/batch: ',tochar(splitcount),' )'
                else
                    write(rl_iou,*) '    divided mesh, iteration ',tochar(iter,2),'  ( pts/batch: ',tochar(splitcount(1:8)),' ... )'
                endif
                if ( iter .eq. size(divisors) ) then
                write(rl_iou,*) '    =============================================================='
                endif
                deallocate(splitcount)
            endif

        enddo iterloop

        ! A small sanity test perhaps. Never hurts.
        do is=1,nset
            if ( set(is)%n .ne. setctr(is,size(divisors)+1) ) then
                call rl_stop_gracefully(['Could not divide the points the way I thought I could'],rl_exitcode_symmetry)
            endif
        enddo

    !end block divide

    !@TODO insert better division here.
    ! My idea was to use the centers of mass of the current distribution
    ! and calculate the Voronoi diagram. This should then get reweighted
    ! until it has a constant number of points per cell. But that is too
    ! annoying to code up now. Should not be complicated or slow, only
    ! takes time to write up.
end subroutine



!> evaluate the partition function
subroutine evaluate_partition_function(batch,p,basis,free,ec,mw,&
    r_radial,r_angular,outer_partition_radius)
    !> batch of points
    class(rl_integration_batch), dimension(:), intent(inout) :: batch
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> free atom quantities
    type(rl_free_atom), intent(in) :: free
    !> extended cluster of atoms
    type(rl_extended_cluster), intent(in) :: ec
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:), intent(in) :: outer_partition_radius

    ! local
    real(r8), parameter :: tinybuf=1E-9_r8
    type(rl_extended_cluster_hashedspace) :: hs
    type(rl_distancetable) :: dt
    real(r8), dimension(:), allocatable :: inner_cutoff_sq,outer_cutoff
    real(r8), dimension(3) :: va
    real(r8) :: f0
    integer :: ib,ip,ie,iatm,ispc,irad,iang

    ! Get the outer cutoff radius per unitcell atom
    allocate(inner_cutoff_sq(p%n_atom))
    allocate(outer_cutoff(p%n_atom))
    inner_cutoff_sq=0.0_r8
    outer_cutoff=0.0_r8
    do iatm=1,p%n_atom
        ispc=p%species(iatm)
        ! Pick the largest number out of several to choose from:
        f0=0.0_r8
        ! First the largest basis function cutoff?
        f0=max(f0,maxval(basis%species(ispc)%basis_cutoff))
        ! Then the magic number from AIMS?
        f0=max(f0,outer_partition_radius(ispc))
        outer_cutoff(iatm)=f0
    enddo
    ! Calculate the unitcell distance table. Might retire.
    call dt%generate(p%fractional_coordinate,p%latticevectors,maxval(outer_cutoff)+tinybuf,verbosity=-1,mw=mw)
    ! Get the inner cutoff
    do iatm=1,p%n_atom
        ispc=p%species(iatm)
        ! Smallest cutoff: first check smallest grid radius
        f0=free%species(ispc)%rho%r_min
        ! Then safe Stratmann sphere thingy. 0.64 is the Stratmann a parameter.
        if ( dt%particle(iatm)%n .gt. 1 ) then
            f0=max(f0,dt%particle(iatm)%d(2)*0.5_r8*(1-0.64_r8))
            inner_cutoff_sq(iatm)=f0**2
        else
            inner_cutoff_sq(iatm)=f0**2
        endif
    enddo
    ! Create hashed space for faster lookup. Maybe squared cutoff? Not sure.
    call hs%generate(ec,p,maxval(outer_cutoff)+tinybuf,mw)

    do ib=1,size(batch)
    do ip=1,batch(ib)%n_point
        ! fetch the point in sensible coordinates
        iatm=batch(ib)%index_atom(ip)
        ispc=p%species(iatm)
        irad=batch(ib)%index_radial(ip)
        iang=batch(ib)%index_angular(ip)
        ie=ec%locate_unit_cell_atom_in_cluster(iatm)
        va=r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)+ec%cartesian_coordinate(:,ie)
        call rl_stratmann_smoother_partition_function(ec,hs,p,iatm,va,outer_cutoff(iatm),inner_cutoff_sq,batch(ib)%partition_function(ip))
    enddo
    enddo

    !@TODO insert destructors
    deallocate(inner_cutoff_sq)
    deallocate(outer_cutoff)
end subroutine

!> take a set of irreducible points divide them into zero and nonzero partition function
subroutine sort_by_partition_function(p,ec,basis,free,sym,mw,r_radial,r_angular,outer_partition_radius,&
    n_irr_point_local,n_irr_zero_point_local,irr_ind,irr_zero_ind)
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis
    !> free atom quantities
    type(rl_free_atom), intent(in) :: free
    !> symmetry handle
    type(rl_spacegroup), intent(in) :: sym
    !> MPI helper
    type(rl_mpi_helper), intent(inout) :: mw
    !> grid information from AIMS
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:), intent(in) :: outer_partition_radius
    !> how many irreducible points do I have on this rank
    integer, intent(inout) :: n_irr_point_local
    !> how many points with zero partition function do I have on this rank
    integer, intent(out) :: n_irr_zero_point_local
    !> index arrays that help map out symmetry
    integer, dimension(:,:), allocatable, intent(inout) :: irr_ind
    !> index arrays for the points with zero partition function
    integer, dimension(:,:), allocatable, intent(out) :: irr_zero_ind

    ! temp space for the partition function
    real(r8), dimension(:), allocatable :: irr_partitionfunction

    !getpf: block
        real(r8), parameter :: tinybuf=1E-9_r8
        type(rl_extended_cluster_hashedspace) :: hs
        type(rl_distancetable) :: dt
        real(r8), dimension(:), allocatable :: inner_cutoff_sq,outer_cutoff
        real(r8), dimension(3) :: va
        real(r8) :: f0
        integer :: i,j,iirr,iatm,ispc,irad,iang

    !divideanddistr: block
        integer, dimension(:,:), allocatable :: zerobuf,nonzerobuf
        !integer :: n_nonzero,n_zero,nrow,i
        integer :: n_nonzero,n_zero,nrow
        integer :: n_nonzero_global,n_zero_global

    ! First calculate the partition function for all points
    !getpf: block

        ! Get the outer cutoff radius per unitcell atom
        allocate(inner_cutoff_sq(p%n_atom))
        allocate(outer_cutoff(p%n_atom))
        inner_cutoff_sq=0.0_r8
        outer_cutoff=0.0_r8
        do iatm=1,p%n_atom
            ispc=p%species(iatm)
            ! Pick the largest number out of several to choose from:
            f0=0.0_r8
            ! First the largest basis function cutoff?
            f0=max(f0,maxval(basis%species(ispc)%basis_cutoff))
            ! Then the magic number from AIMS?
            f0=max(f0,outer_partition_radius(ispc))
            outer_cutoff(iatm)=f0
        enddo
        ! Calculate the unitcell distance table. Might retire.
        call dt%generate(p%fractional_coordinate,p%latticevectors,maxval(outer_cutoff)+tinybuf,verbosity=-1,mw=mw)
        ! Get the inner cutoff
        do iatm=1,p%n_atom
            ispc=p%species(iatm)
            ! Smallest cutoff: first check smallest grid radius
            f0=free%species(ispc)%rho%r_min
            ! Then safe Stratmann sphere thingy. 0.64 is the Stratmann a parameter.
            if ( dt%particle(iatm)%n .gt. 1 ) then
                f0=max(f0,dt%particle(iatm)%d(2)*0.5_r8*(1-0.64_r8))
                inner_cutoff_sq(iatm)=f0**2
            else
                inner_cutoff_sq(iatm)=f0**2
            endif
        enddo
        ! Create hashed space for faster lookup. Maybe squared cutoff? Not sure.
        call hs%generate(ec,p,maxval(outer_cutoff)+tinybuf,mw)
        ! Some space for the partition function
        if ( n_irr_point_local .gt. 0 ) then
            allocate(irr_partitionfunction(n_irr_point_local))
            irr_partitionfunction=0.0_r8
        else
            allocate(irr_partitionfunction(1))
            irr_partitionfunction=0.0_r8
        endif

        ! Now calculate the actual partition function
        do i=1,n_irr_point_local
            ! Coordinate of the point in sensible coordinates
            iirr=irr_ind(1,i)
            iatm=sym%irr_to_all(iirr)
            ispc=p%species(iatm)
            irad=irr_ind(3,i)
            iang=irr_ind(4,i)
            j=ec%locate_unit_cell_atom_in_cluster(iatm)
            va=r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)+ec%cartesian_coordinate(:,j)
            ! and the partition function
            call rl_stratmann_smoother_partition_function(ec,hs,p,iatm,va,outer_cutoff(iatm),inner_cutoff_sq,irr_partitionfunction(i))
        enddo

        !@TODO insert destructors
        deallocate(inner_cutoff_sq)
        deallocate(outer_cutoff)
    !end block getpf

    ! Now divide the points into those with zero or nonzero partition function.
    !divideanddistr: block

        ! Count number of zero or nonzero
        n_nonzero=0
        n_zero=0
        do i=1,n_irr_point_local
            if ( irr_partitionfunction(i) .gt. rl_tiny ) then
                n_nonzero=n_nonzero+1
            else
                n_zero=n_zero+1
            endif
        enddo

        ! Figure out number of rows in the index array
        if ( n_irr_point_local .gt. 0 ) then
            nrow=size(irr_ind,1)
        else
            nrow=0
        endif
        call mw%allreduce('max',nrow)

        ! space for temporary buffers
        if ( n_nonzero .gt. 0 ) then
            allocate(nonzerobuf(nrow,n_nonzero))
            nonzerobuf=0
        else
            allocate(nonzerobuf(1,1))
            nonzerobuf=-rl_hugeint
        endif
        if ( n_zero .gt. 0 ) then
            allocate(zerobuf(nrow,n_zero))
            zerobuf=0
        else
            allocate(zerobuf(1,1))
            zerobuf=-rl_hugeint
        endif

        ! store indices in buffers
        n_nonzero=0
        n_zero=0
        do i=1,n_irr_point_local
            if ( irr_partitionfunction(i) .gt. rl_tiny ) then
                n_nonzero=n_nonzero+1
                nonzerobuf(:,n_nonzero)=irr_ind(:,i)
            else
                n_zero=n_zero+1
                zerobuf(:,n_zero)=irr_ind(:,i)
            endif
        enddo

        call mw%allreduce('sum',n_nonzero,n_nonzero_global)
        call mw%allreduce('sum',n_zero,n_zero_global)

        ! Distribute these points evenly. Might be overkill, but you never know.
        if ( n_nonzero_global .gt. 0 ) then
            call rl_alltoall_distribute_2d_int(nonzerobuf,nrow,n_nonzero,mw)
            ! Update the output!
            deallocate(irr_ind)
            if ( n_nonzero .gt. 0 ) then
                allocate(irr_ind(nrow,n_nonzero))
                irr_ind=nonzerobuf
                n_irr_point_local=n_nonzero
            else
                allocate(irr_ind(1,1))
                irr_ind=-rl_hugeint
                n_irr_point_local=0
            endif
        else
            call rl_stop_gracefully(['No nonzero points, can not be'],rl_exitcode_symmetry,mw%comm)
        endif

        if ( n_zero_global .gt. 0 ) then
            call rl_alltoall_distribute_2d_int(zerobuf,nrow,n_zero,mw)
            if ( n_nonzero .gt. 0 ) then
                allocate(irr_zero_ind(nrow,n_zero))
                irr_zero_ind=zerobuf
                n_irr_zero_point_local=n_zero
            else
                allocate(irr_zero_ind(1,1))
                irr_zero_ind=-rl_hugeint
                n_irr_zero_point_local=0
            endif
        else
            allocate(irr_zero_ind(1,1))
            irr_zero_ind=-rl_hugeint
            n_irr_zero_point_local=0
        endif

        ! And a little cleanup
        deallocate(zerobuf)
        deallocate(nonzerobuf)
    !end block divideanddistr
end subroutine

!> convert a set of points to a proper integration batch, this version is for irreducible batches.
subroutine convert_set_into_irreducible_batches(p,sym,set,n_set,r_radial,r_angular,w_radial,w_angular,n_angular,batch,n_batch,ml,mw)
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> spacegroup
    type(rl_spacegroup), intent(in) :: sym
    !> set of points
    type(rl_set_of_points), dimension(:), allocatable, intent(in) :: set
    !> how many sets
    integer, intent(inout) :: n_set
    !> grid generating information
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:,:), intent(in) :: w_radial
    real(r8), dimension(:,:,:), intent(in) :: w_angular
    integer, dimension(:,:), intent(in) :: n_angular
    !> batches of points
    type(rl_integration_batch), dimension(:), allocatable, intent(out) :: batch
    !> number of batches
    integer, intent(out) :: n_batch
    !> MPI helper, semilocal
    type(rl_mpi_helper), intent(inout) :: ml
    !> MPI helper, global
    type(rl_mpi_helper), intent(inout) :: mw

    integer, dimension(:), allocatable :: ni_per_rank,oi_per_rank
    real(r8), dimension(3) :: v0,iradvec,jradvec
    integer :: ib,ni,i,j,k,ii,jj !,l
    integer :: irad,iang,jang,ispc,iirr,iatm,jatm,imul,jmul,iopr,jopr,kopr,n1
    integer :: irr_off_ctr

    ! At least we know the number of batches and points now.
    n_batch=n_set
    allocate(batch(n_batch))

    ! Do a first pass and unpack the data
    do ib=1,n_batch
        ! shorthand for the number of irreducible points in this batch
        ni=set(ib)%n
        ! figure out how many unfolded points there are in the batch
        n1=0
        do i=1,ni
            !l=l+1
            iirr=set(ib)%ind(1,i)      ! irreducible atom
            imul=set(ib)%ind(5,i)      ! number of outfolded points for this atom
            jmul=sym%irr_unfold_ctr(iirr)  ! adjust with number of atoms it folds out to
            n1=max(n1,imul*jmul)           ! accumulate largest number of unfolded atoms per irreducible atom
        enddo
        ! Make space in the batch
        batch(ib)%n_point = ni
        allocate(batch(ib)%folded_coordinate(3,ni))
        allocate(batch(ib)%radial_coordinate(3,ni))
        allocate(batch(ib)%index_irr_atom(ni))
        allocate(batch(ib)%index_atom(ni))
        allocate(batch(ib)%index_radial(ni))
        allocate(batch(ib)%index_angular(ni))
        allocate(batch(ib)%integration_weight(ni))
        allocate(batch(ib)%partition_function(ni))
        allocate(batch(ib)%unfold_ctr(ni))
        allocate(batch(ib)%unfold_operation(n1,ni))
        allocate(batch(ib)%unfold_index_radial(n1,ni))
        allocate(batch(ib)%unfold_index_angular(n1,ni))
        allocate(batch(ib)%unfold_atom(n1,ni))
        batch(ib)%folded_coordinate=0.0_r8
        batch(ib)%radial_coordinate=0.0_r8
        batch(ib)%index_irr_atom=0
        batch(ib)%index_atom=0
        batch(ib)%index_radial=0
        batch(ib)%index_angular=0
        batch(ib)%partition_function=0.0_r8
        batch(ib)%integration_weight=0.0_r8
        batch(ib)%unfold_ctr=0
        batch(ib)%unfold_operation=0
        batch(ib)%unfold_index_radial=0
        batch(ib)%unfold_index_angular=0
        batch(ib)%unfold_atom=0
    enddo

    ! To get the offsets properly, I need the number of points per rank, and their offset.
    allocate(ni_per_rank(ml%n))
    allocate(oi_per_rank(ml%n))
    ni_per_rank=0
    oi_per_rank=0
    ni_per_rank(ml%r+1)=sum(batch(:)%n_point)
    call ml%allreduce('sum',ni_per_rank)
    n1=0
    do i=1,ml%n
        oi_per_rank(i)=n1
        n1=n1+ni_per_rank(i)
    enddo

    ! Now do a second pass and actually store data everywhere.
    irr_off_ctr=oi_per_rank(ml%r+1)  ! offsets in semi-global arrays for irreducible points
    do ib=1,n_batch
        ! shorthand for the number of points in this batch
        ni=batch(ib)%n_point
        ! Store the offsets
        batch(ib)%semilocal_irr_offset=irr_off_ctr
        irr_off_ctr=irr_off_ctr+ni
        do i=1,ni
            !l=l+1
            iirr=set(ib)%ind(1,i)       ! irreducible atom index
            iatm=sym%irr_to_all(iirr)   ! atom index
            ispc=p%species(iatm)        ! species index
            irad=set(ib)%ind(3,i)       ! radial index
            iang=set(ib)%ind(4,i)       ! angular index

            ! Vector that points radially from the irreducible atom, deduced
            ! from the indices and the arrays from AIMS.
            iradvec=r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)

            ! Actually, seems I'm not quite done with tests just yet.
            ! Check that the unfolding to the full grid makes sense.
            imul=set(ib)%ind(5,i)          ! number of outfolded points for this atom
            jmul=sym%irr_unfold_ctr(iirr)  ! adjust with number of atoms it folds out to
            ii=0
            do j=1,imul
            do k=1,jmul
                iopr=set(ib)%ind(5+j,i)                  ! operation to unfold the gridpoint
                jopr=sym%irr_unfold_operation(k,iirr)    ! operation to unfold atom
                kopr=sym%multiplication_table(jopr,iopr) ! resulting operation
                ! Locate angular index. If it eventually works, will replace with
                ! box-thingy to make it significantly faster.
                jang=-1
                v0=r_angular(:,iang,irad,ispc)
                v0=matmul(sym%op(kopr)%m,v0)
                do jj=1,n_angular(irad,ispc)
                    if ( rl_sqnorm(v0-r_angular(:,jj,irad,ispc)) .lt. rl_sqtol ) then
                        jang=jj
                        exit
                    endif
                enddo
                if ( jang .le. 0 ) then
                    call rl_stop_gracefully(['Mis-counted symmetry information.'],rl_exitcode_symmetry,mw%comm)
                endif
                jatm=sym%irr_unfold_index(k,iirr)        ! index of rotated atom

                ! Fetch transformed radial vector, I hope.
                jradvec=r_angular(:,jang,irad,ispc)*r_radial(irad,ispc)
                ! Check that transformation is actually the transformation I think it is.
                if ( rl_sqnorm(matmul(sym%op(kopr)%m,iradvec)-jradvec) .gt. rl_sqtol ) then
                    call rl_stop_gracefully(['Mis-indexed symmetry information.'],rl_exitcode_symmetry,mw%comm)
                endif
                if ( rl_sqnorm(matmul(sym%op(kopr)%im,jradvec)-iradvec) .gt. rl_sqtol ) then
                    call rl_stop_gracefully(['Mis-indexed symmetry information.'],rl_exitcode_symmetry,mw%comm)
                endif
                ! Ok, seems I was not an idiot if we make it here.
                ii=ii+1
                batch(ib)%unfold_operation(ii,i)=kopr     ! operation to unfold with
                batch(ib)%unfold_index_angular(ii,i)=jang ! angular index to unfold to
                batch(ib)%unfold_index_radial(ii,i)=irad  ! radial index to unfold to
                batch(ib)%unfold_atom(ii,i)=jatm          ! atom to unfold to
            enddo
            enddo
            batch(ib)%unfold_ctr(i)=imul*jmul

            ! Now for the moment of truth, does this agree with the point that I have?
            ! Once I am satisfied I can probably remove this check. But it does not cost
            ! anything in particular anyway.
            v0=iradvec+p%cartesian_coordinate(:,iatm)
            v0=matmul(p%inv_latticevectors,v0)
            v0=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
            v0=matmul(p%latticevectors,v0)
            if ( rl_sqnorm(v0-set(ib)%r(:,i)) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['I did not keep track of my indices properly. Shame on me.'],rl_exitcode_symmetry,mw%comm)
            endif

            ! Store the vectors and some indices
            batch(ib)%index_irr_atom(i)=iirr
            batch(ib)%index_atom(i)=iatm
            batch(ib)%index_radial(i)=irad
            batch(ib)%index_angular(i)=iang
            batch(ib)%folded_coordinate(:,i)=v0
            batch(ib)%radial_coordinate(:,i)=iradvec
            ! This seems to be the correct way of getting the total weight per point.
            ! Factor 2*pi confuses me severely. But defined this way, integration_weight*partition_function integrates to volume.
            batch(ib)%integration_weight(i)=w_radial(irad,ispc)*w_angular(iang,irad,ispc)*(r_radial(irad,ispc)**2)*2*rl_twopi
            batch(ib)%integration_weight(i)=batch(ib)%integration_weight(i)*batch(ib)%unfold_ctr(i)
        enddo
    enddo
end subroutine

!> convert a set of points to a proper integration batch, this version is for irreducible batches.
subroutine convert_set_into_full_batches(p,set,n_set,r_radial,r_angular,w_radial,w_angular,n_angular,batch,n_batch,ml,mw)
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> set of points
    type(rl_set_of_points), dimension(:), allocatable, intent(in) :: set
    !> how many sets
    integer, intent(inout) :: n_set
    !> grid generating information
    real(r8), dimension(:,:), intent(in) :: r_radial
    real(r8), dimension(:,:,:,:), intent(in) :: r_angular
    real(r8), dimension(:,:), intent(in) :: w_radial
    real(r8), dimension(:,:,:), intent(in) :: w_angular
    integer, dimension(:,:), intent(in) :: n_angular
    !> batches of points
    type(rl_integration_batch), dimension(:), allocatable, intent(out) :: batch
    !> number of batches
    integer, intent(out) :: n_batch
    !> MPI helper, semilocal
    type(rl_mpi_helper), intent(inout) :: ml
    !> MPI helper, global
    type(rl_mpi_helper), intent(inout) :: mw

    real(r8), dimension(3) :: v0,iradvec
    integer :: ib,ni,i
    integer :: irad,iang,ispc,iirr,iatm

    ! At least we know the number of batches and points now.
    n_batch=n_set
    allocate(batch(n_batch))

    ! Do a first pass and unpack the data
    do ib=1,n_batch
        ! shorthand for the number of irreducible points in this batch
        ni=set(ib)%n
        ! Make space in the batch
        batch(ib)%n_point = ni
        allocate(batch(ib)%folded_coordinate(3,ni))
        allocate(batch(ib)%radial_coordinate(3,ni))
        allocate(batch(ib)%index_atom(ni))
        allocate(batch(ib)%index_radial(ni))
        allocate(batch(ib)%index_angular(ni))
        allocate(batch(ib)%integration_weight(ni))
        allocate(batch(ib)%partition_function(ni))
        batch(ib)%folded_coordinate=0.0_r8
        batch(ib)%radial_coordinate=0.0_r8
        batch(ib)%index_atom=0
        batch(ib)%index_radial=0
        batch(ib)%index_angular=0
        batch(ib)%partition_function=0.0_r8
        batch(ib)%integration_weight=0.0_r8
    enddo

    ! Now do a second pass and actually store data everywhere.
    do ib=1,n_batch
        ! shorthand for the number of points in this batch
        ni=batch(ib)%n_point
        ! Store the offsets
        do i=1,ni
            iatm=set(ib)%ind(2,i)       ! atom index
            irad=set(ib)%ind(3,i)       ! radial index
            iang=set(ib)%ind(4,i)       ! angular index
            ispc=p%species(iatm)        ! species index

            ! Vector that points radially from the irreducible atom, deduced
            ! from the indices and the arrays from AIMS.
            iradvec=r_angular(:,iang,irad,ispc)*r_radial(irad,ispc)

            ! Now for the moment of truth, does this agree with the point that I have?
            ! Once I am satisfied I can probably remove this check. But it does not cost
            ! anything in particular anyway.
            v0=iradvec+p%cartesian_coordinate(:,iatm)
            v0=matmul(p%inv_latticevectors,v0)
            v0=rl_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
            v0=matmul(p%latticevectors,v0)
            if ( rl_sqnorm(v0-set(ib)%r(:,i)) .gt. rl_sqtol ) then
                call rl_stop_gracefully(['I did not keep track of my indices properly. Shame on me.'],rl_exitcode_symmetry,mw%comm)
            endif

            ! Store the vectors and some indices
            batch(ib)%index_atom(i)=iatm
            batch(ib)%index_radial(i)=irad
            batch(ib)%index_angular(i)=iang
            batch(ib)%folded_coordinate(:,i)=v0
            batch(ib)%radial_coordinate(:,i)=iradvec
            ! This seems to be the correct way of getting the total weight per point.
            ! Factor 4*pi confuses me severely. But defined this way, integration_weight*partition_function integrates to volume.
            batch(ib)%integration_weight(i)=w_radial(irad,ispc)*w_angular(iang,irad,ispc)*(r_radial(irad,ispc)**2)*2*rl_twopi
        enddo
    enddo
end subroutine

!> count number of atoms within cutoff and calculate other auxiliary things for a set of batches
subroutine get_aux_information_for_batches(batch,p,ec,basis)
    !> batches of points
    class(rl_integration_batch), dimension(:), intent(inout) :: batch
    !> structure
    type(rl_crystalstructure), intent(in) :: p
    !> extended cluster
    type(rl_extended_cluster), intent(in) :: ec
    !> basis set
    type(rl_lcao_basis_set), intent(in) :: basis

    real(r8), parameter :: tinybuf=1E-7_r8
    real(r8), dimension(:), allocatable :: cutoff_per_atom
    real(r8), dimension(3) :: v0
    real(r8) :: f0,f1,sqrc,rc
    integer, dimension(:), allocatable :: included,di
    integer :: iatm,ispc,ib,ctr,i,j,k

    allocate(included(ec%n_extended_atom))
    allocate(di(ec%n_extended_atom))
    allocate(cutoff_per_atom(p%n_atom))
    included=0
    di=0
    cutoff_per_atom=0.0_r8

    do iatm=1,p%n_atom
        ispc=p%species(iatm)
        ! Pick the largest number out of several to choose from:
        f0=0.0_r8
        ! First the largest basis function cutoff?
        f0=max(f0,maxval(basis%species(ispc)%basis_cutoff))
        ! Something else I will think of later ...
        cutoff_per_atom(iatm)=f0 ! Add a small buffer to be on the safe side.
    enddo

    ! cutoff, with a slight buffer
    rc=maxval(cutoff_per_atom)+tinybuf
    sqrc=rc**2

    do ib=1,size(batch)
        ! First get the bounding sphere real stupid, just get the center of mass
        ! and then draw the bounding sphere from that. It might be good enough?
        v0=0.0_r8
        do i=1,batch(ib)%n_point
            v0=v0+batch(ib)%folded_coordinate(:,i)
        enddo
        v0=v0/real(batch(ib)%n_point,r8)
        f1=0.0_r8
        do i=1,batch(ib)%n_point
            f0=rl_sqnorm(v0-batch(ib)%folded_coordinate(:,i))
            f1=max(f0,f1)
        enddo
        ! Now I have the center of mass and bounding sphere.
        batch(ib)%center_of_mass=v0
        batch(ib)%bounding_sphere_radius=sqrt(f1)

        ! Pre-fetch a list atoms that can touch the bounding sphere?
        di=-1
        ctr=0
        do i=1,ec%n_extended_atom
            ! Kinda dumb I know but if you read this far you deserve to judge me.
            f0=norm2( ec%cartesian_coordinate(:,i)-batch(ib)%center_of_mass )-2*batch(ib)%bounding_sphere_radius
            f0=f0**2
            !if ( f0 .lt. ec%confinement_cutoff_sq(i) ) then
            if ( f0 .lt. sqrc ) then
                ctr=ctr+1
                di(ctr)=i
            endif
        enddo

        ! store those that can be within the cutoff. Not a bad thing to have.
        included=0
        do i=1,batch(ib)%n_point
        do k=1,ctr
            j=di(k)
            f0=rl_sqnorm( batch(ib)%folded_coordinate(:,i) - ec%cartesian_coordinate(:,j) )
            !if ( f0 .lt. ec%confinement_cutoff_sq(j) ) then
            if ( f0 .lt. sqrc ) then
                included(j)=1
            endif
        enddo
        enddo

        ! can not understand why this was not done before.
        ! Does not cost too much memory. Or maybe it does? I dunno.
        batch(ib)%n_relevant_ext_atom=sum(included)
        allocate(batch(ib)%relevant_ext_atom( batch(ib)%n_relevant_ext_atom ))
        batch(ib)%relevant_ext_atom=-1
        j=0
        do i=1,ec%n_extended_atom
            if ( included(i) .ne. 0 ) then
                j=j+1
                batch(ib)%relevant_ext_atom(j)=i
            endif
        enddo
        batch(ib)%n_relevant_ext_atom=j
    enddo

    deallocate(included)
    deallocate(di)
    deallocate(cutoff_per_atom)
end subroutine

!> get the global index of an integration point
elemental function global_index(idx,i_atom,i_radial,i_angular) result(i_global)
    !> indexing helper
    class(rl_integration_grid_index), intent(in) :: idx
    !> index to atom in the unit cell
    integer, intent(in) :: i_atom
    !> radial index
    integer, intent(in) :: i_radial
    !> angular index
    integer, intent(in) :: i_angular
    !> global index
    integer :: i_global

    i_global=idx%atom(i_atom)%radial(i_radial)+i_angular
end function

!> get the radial integration weight of an integration points
elemental function fetch_radial_weight(wts,i_species,i_radial) result(w)
    !> weight helper
    class(rl_integration_grid_weights), intent(in) :: wts
    !> index to species
    integer, intent(in) :: i_species
    !> radial index
    integer, intent(in) :: i_radial
    !> radial integration weight
    real(r8) :: w
    ! just fetch
    w=wts%species(i_species)%rad(i_radial)%radial_weight
end function

!> get the radial integration weight of an integration points
elemental function fetch_radial_scale(wts,i_species) result(s)
    !> weight helper
    class(rl_integration_grid_weights), intent(in) :: wts
    !> index to species
    integer, intent(in) :: i_species
    !> radial scale
    real(r8) :: s
    ! just fetch
    s=wts%species(i_species)%radial_scale
end function

!> get the radial integration weight of an integration points
elemental function fetch_n_radial(wts,i_species) result(n)
    !> weight helper
    class(rl_integration_grid_weights), intent(in) :: wts
    !> index to species
    integer, intent(in) :: i_species
    !> number of radial shells
    integer :: n
    n=size(wts%species(i_species)%rad)
end function

!> get the angular integration weight of an integration points
elemental function fetch_angular_weight(wts,i_species,i_radial,i_angular) result(w)
    !> weight helper
    class(rl_integration_grid_weights), intent(in) :: wts
    !> index to species
    integer, intent(in) :: i_species
    !> radial index
    integer, intent(in) :: i_radial
    !> angular index
    integer, intent(in) :: i_angular
    !> angular integration weight
    real(r8) :: w
    ! Just fetch it.
    w=wts%species(i_species)%rad(i_radial)%angular_weight(i_angular)
end function

end module
