!****s* FHI-aims/partition_grid
!  NAME
!    partition_grid
!  SYNOPSIS

subroutine partition_grid( )

!  PURPOSE
!    Partitions the integration grid into manageable batches used
!    in integration of the Hamilton matrix and in the update of the
!    electron density.
!
!  USES

  use batch_statistics, only: batch_sizes_stats
  use dimensions
  use runtime_choices
  use applicable_citations
  use grids
  use geometry
  use pbc_lists
  use mpi_utilities
  use synchronize_mpi
  use localorb_io, only: output_priority, localorb_info, use_unit, OL_high, &
      OL_norm, OL_low
  use statistics, only: stats_int, calc_stats_array
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  
  implicit none

!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none -- batch information is stored in the variable batches on
!            each thread
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


  integer :: i_atom, i_radial, i_angular
  integer :: i_coord
  integer :: i_point, i_batch
  integer :: i_index
  integer :: n_points_in_grid, points_out
  integer :: i_my_batch
  integer :: current_batch
  integer :: i_my_point, n_my_points
  
  integer, dimension(:), allocatable :: batch_sizes
  real*8, dimension(:,:), allocatable :: batch_coords

  real*8, dimension(:,:), allocatable :: all_coords
  real*8, dimension(:), allocatable :: point_weight
  real*8, dimension(3) :: coord_current
  real*8, dimension(3) :: coord_current_temp
  integer, dimension(3) :: coord_current_int

  real*8 :: min_point_weight, max_point_weight
  real*8 :: total_min_point_weight, total_max_point_weight

  real*8 dist_tab_sq(n_centers_integrals)
  real*8 dir_tab(3,n_centers_integrals)
  integer :: i_basis(n_centers_basis_I)
  integer :: n_compute_a, n_compute_c
  integer :: max_n_max_compute
  integer :: n_max_compute_part

  character*100 :: info_str
  integer:: info

  integer, dimension(:), allocatable :: grid_partition

  integer, dimension(:), allocatable :: adj_rows
  integer, dimension(:), allocatable :: adjacency
  integer :: n_total_edges

  type tree_node
     integer :: size
     integer, dimension(:), allocatable :: point_indexes
     integer, dimension(:), allocatable :: global_point_indexes
     integer :: n_children
     type(tree_node), pointer, dimension(:) :: children
  end type tree_node

  type point_group
     integer :: size
     integer, pointer, dimension(:) :: point_indexes
     integer :: level
  end type point_group

  type(tree_node) :: root

  type(point_group), pointer, dimension(:) :: base_group

  real*8, dimension(3) :: coord_max, coord_min, bb_dim, h
  integer, dimension(3) :: n_bins, idx, h_int, coord_min_int, remainder

  real*8 time0, time_work, time2, time_work2

  integer :: hash_value, hash_method
  real*8 :: hash_multiplier = 1.0e6

  logical :: debug = .false.
  integer :: max_n_my_points, min_n_my_points

  integer i_center_ii
  real*8, dimension(:,:), allocatable:: coords_center_combined    ! coords_center_combined(3,n_centers)

  !----------------------------------------------------------------------------

  ! this variable now lives in module grids.f
  if (grid_partitioned) return

  ! References using the FHI-aims internal grids should cite:
  call cite_reference("FHI-aims-grids")

  ! if we perform a supercell calculation:
  if (n_periodic > 0 .and. flag_cpu_consistency) then
     ! check for consistent lattice vectors
     call check_cpu_consistency_matrix(lattice_vector, &
        "lattice for partition grid", 3, n_periodic)

     ! check for consistent map_to_center_cell_matrix
     call check_cpu_consistency_matrix(map_to_center_cell_matrix, &
        "map_to_center_cell_matrix for partition grid", 3, n_periodic)

     ! check for consistent coords
     call check_cpu_consistency_matrix (coords,"coords",3,n_atoms)

  endif

  ! then, let's count the total number of points in the grid
  n_points_in_grid = 0
  do i_atom = 1, n_atoms, 1
     do i_radial = 1, n_radial(species(i_atom)), 1
        n_points_in_grid = n_points_in_grid + &
             n_angular( i_radial,species(i_atom) )
     end do
  end do

  ! the case of parallel grid partitioning
  if (parallel_grid_partitioning) then

     if (grid_partitioning_method == 6) then
        call find_bb_dimensions(coord_max, coord_min, bb_dim)
        if (use_integer_hash) coord_min_int = INT(coord_min*hash_multiplier)
        call find_bins(bb_dim, n_grid_batches, n_bins, h)
        if (use_integer_hash) h_int = INT(h*hash_multiplier)
     end if

     if (use_hashed_maxmin) then
        hash_method = 2
     else
        hash_method = 1
     end if

     n_my_points = 0

     select case(hash_method)
     case(1)

        do i_point = 1, n_points_in_grid
           if (myid == MOD(i_point, n_tasks)) n_my_points = n_my_points + 1
        end do

     case(2)

        do i_atom = 1, n_atoms, 1
           do i_radial = 1, n_radial(species(i_atom)), 1
              do i_angular = 1, n_angular( i_radial,species(i_atom) )

                 coord_current(:) = coords_center( :,i_atom ) + &
                      r_angular(:, i_angular, i_radial, species(i_atom)) * &
                      r_radial(i_radial, species(i_atom))
        
                 coord_current_temp = coord_current 

                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current_temp(1:3) )
                 end if

                 if (use_integer_hash) then

                    coord_current_int = INT(coord_current_temp*hash_multiplier)
                    
                    do i_index = 1, 3, 1
                       remainder(i_index) = MODULO(coord_current_int(i_index) - coord_min_int(i_index), h_int(i_index))
                       idx(i_index) = &
                            (coord_current_int(i_index) - coord_min_int(i_index) - remainder(i_index)) / h_int(i_index) + 1
                    end do

                 else

                    idx = CEILING((coord_current_temp - coord_min)/h)

                 end if

                 idx = MIN(idx, n_bins)
                 idx = MAX(idx, 1)
        
                 hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*n_bins(2)*n_bins(1)
        
                 if (myid+1 == hash_value) n_my_points = n_my_points + 1
              end do
           end do
        end do

     end select

     ! allocate grid storage for the current thread only
     call aims_allocate(all_coords, 3,n_my_points, "all_coords")

     root%size = n_my_points
     allocate( root%point_indexes(root%size),stat=info)
     call check_allocation(info, ' root%point_indexes           ')

     allocate( root%global_point_indexes(root%size),stat=info )
     call check_allocation(info, 'root%global_point_indexes     ')

     ! allocate array for the weights for the points in the current thread
     call aims_allocate(point_weight, n_my_points, "point_weight")


     if(myid .eq. 0) print*, "n_centers_integrals=", n_centers_integrals, " n_centers=", n_centers
     if(n_centers .ne. n_centers_integrals) then
         print*, "n_centers_integrals .ne. n_centers! Not applicable to current optimization !", &
                 "n_centers_integrals=", n_centers_integrals, " n_centers=", n_centers
         call aims_stop
     endif
     allocate(coords_center_combined(3,n_centers_integrals),stat=info)
     call check_allocation(info, ' coords_center_combined  ')

     do i_center_ii = 1, n_centers_integrals, 1
         coords_center_combined(1,i_center_ii) = coords_center(1,centers_basis_integrals(i_center_ii))
         coords_center_combined(2,i_center_ii) = coords_center(2,centers_basis_integrals(i_center_ii))
         coords_center_combined(3,i_center_ii) = coords_center(3,centers_basis_integrals(i_center_ii))
     enddo

     time2 = mpi_wtime()

     n_max_compute_part = 0
     i_point = 0
     i_my_point = 0
     do i_atom = 1, n_atoms, 1
        do i_radial = 1, n_radial(species(i_atom)), 1
           do i_angular = 1, n_angular( i_radial,species(i_atom) )
                 
              i_point = i_point + 1
              coord_current(:) = coords_center( :,i_atom ) + &
                   r_angular(:, i_angular, i_radial, species(i_atom)) * &
                   r_radial(i_radial, species(i_atom))
              
              coord_current_temp = coord_current 

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current_temp(1:3) )
              end if

              ! for the parallel routines it is best if the points are distributed
              ! as uniformly as possible over the threads, hence this distribution and not
              ! the one based on radial shells
              select case(hash_method)
              case(1)
                 hash_value = MOD(i_point, n_tasks) + 1
              case(2)

                 if (use_integer_hash) then

                    coord_current_int = INT(coord_current_temp*hash_multiplier)
                    
                    do i_index = 1, 3, 1
                       remainder(i_index) = MODULO(coord_current_int(i_index) - coord_min_int(i_index), h_int(i_index))
                       idx(i_index) = &
                            (coord_current_int(i_index) - coord_min_int(i_index) - remainder(i_index)) / h_int(i_index) + 1
                    end do

                 else

                    idx = CEILING((coord_current_temp - coord_min)/h)

                 end if

                 idx = MIN(idx, n_bins)
                 idx = MAX(idx, 1)
                 
                 hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*n_bins(2)*n_bins(1)
              end select

              if (myid+1 == hash_value) then
                 i_my_point = i_my_point + 1
                    
                 all_coords(:,i_my_point) = coord_current(:)
                 root%point_indexes(i_my_point) = i_my_point
                 root%global_point_indexes(i_my_point) = i_point
                 
                 n_compute_c = 0
                 n_compute_a = 0
                 i_basis = 0
                 
                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current(1:3) )
                 end if
              
                 ! compute atom-centered coordinates of current integration point,
                 ! as viewed from all atoms
               !   call tab_atom_centered_coords_p0 &
               !        ( coord_current,  dist_tab_sq, dir_tab, &
               !        n_centers_integrals, centers_basis_integrals )
                 call tab_atom_centered_coords_p0_less &
                      ( coord_current,  dist_tab_sq, &
                      n_centers_integrals, centers_basis_integrals, coords_center_combined)
              
                 ! determine which basis functions are relevant at current integration point,
                 ! and tabulate their indices
                 call prune_basis_p0X ( dist_tab_sq, &
                      n_compute_a, n_compute_c, i_basis,  &
                      n_centers_basis_I, n_centers_integrals, &
                      inv_centers_basis_integrals, i_my_point )
!NEC_CB Using version which saves massively indirect addressed arrays

                 
                 ! store the relaitve weight of current integration point
                 n_max_compute_part = MAX(n_max_compute_part, n_compute_a)
                 point_weight(i_my_point) = dble(n_compute_a)
                    
              end if
                 
           end do
        end do
     end do

     if(myid .eq. 0) print*, "myid=", myid, " partition_grid:part1=", mpi_wtime()-time2

     if (allocated(coords_center_combined)) then
        deallocate(coords_center_combined)
     end if

     if (debug) then
        call sync_find_max(n_my_points, max_n_my_points)
        n_my_points = -n_my_points
        call sync_find_max(n_my_points, min_n_my_points)
        n_my_points = -n_my_points
        min_n_my_points = -min_n_my_points
        write(use_unit,'(A,I5,A,I7,A,I7,A,I7)') 'Myid', myid, ' n_my_points', n_my_points, &
             ' max_points', max_n_my_points, ' min_points', min_n_my_points
     end if

     ! normalize the weights
     call sync_find_max(n_max_compute_part, max_n_max_compute)
     point_weight = point_weight / dble(max_n_max_compute)

     max_point_weight = MAXVAL(point_weight)
     min_point_weight = MINVAL(point_weight)

     call get_max_double( total_max_point_weight, max_point_weight )
     call get_min_double( total_min_point_weight, min_point_weight )

     ! allocate array for the grid partition
     call aims_allocate(grid_partition, n_my_points, "grid_partition")

     grid_partition = 0

     select case(grid_partitioning_method)
        
        ! octree
     case(1)
        write (info_str,'(2X,A,A)') &
             "Partitioning the integration grid into batches ", &
             "with parallel octree method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        write (info_str,'(2X,A,F9.3)') &
             "| Maximal weight for a single point: ", total_max_point_weight
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        write (info_str,'(2X,A,F9.3)') &
             "| Minimal weight for a single point: ", total_min_point_weight
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        n_grid_batches = 0
        call build_partition_tree( root )
!        call remove_small_batches()
        
     case(5)
        write (info_str,'(2X,A,A)') &
             "Partitioning the integration grid into batches ", &
             "with parallel modified maxmin method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        write (info_str,'(2X,A,F9.3)') &
             "| Maximal weight for a single point: ", total_max_point_weight
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        write (info_str,'(2X,A,F9.3)') &
             "| Minimal weight for a single point: ", total_min_point_weight
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        n_grid_batches = 0
        call build_maxmin_partition_tree( root )

     case(6)

        if (.not.use_hashed_maxmin) then
           write (info_str,'(2X,A,A)') &
                "Partitioning the integration grid into batches ", &
                "with parallel hashing method."
        else
           write (info_str,'(2X,A,A)') &
                "Partitioning the integration grid into batches ", &
                "with parallel hashing+maxmin method."           
        end if
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        call hash_points( root )
     end select

     if (allocated(root%point_indexes)) then
        deallocate(root%point_indexes)
     end if
     if (allocated(root%global_point_indexes)) then
        deallocate(root%global_point_indexes)
     end if
     if (allocated(point_weight)) then
        call aims_deallocate(point_weight, "point_weight")
     end if

     ! allocate the batch data size array
     call aims_allocate(batch_sizes, n_grid_batches, "batch_sizes")
     call aims_allocate(batch_coords, 3,n_grid_batches, "batch_coords")

     ! compute the sizes of the batches according to the distribution scheme
     batch_sizes = 0
     batch_coords = 0
     do i_point = 1, n_my_points, 1
        batch_sizes(grid_partition(i_point)) = &
             batch_sizes(grid_partition(i_point)) + 1
        coord_current(:) = all_coords(:,i_point)
        if(n_periodic > 0)then
           call map_to_center_cell( coord_current(1:3) )
        end if
        batch_coords(:,grid_partition(i_point)) = &
           batch_coords(:,grid_partition(i_point)) + coord_current
     end do
     call sync_integer_vector(batch_sizes, n_grid_batches)
     call sync_vector(batch_coords, 3*n_grid_batches)
     batch_coords(1,:) = batch_coords(1,:) / MAX(batch_sizes(:),1)
     batch_coords(2,:) = batch_coords(2,:) / MAX(batch_sizes(:),1)
     batch_coords(3,:) = batch_coords(3,:) / MAX(batch_sizes(:),1)

     ! print out some information on the batches
     call print_out_batch_report( )

     ! While we still have the batch sizes in memory, calculate their statistics
     ! for later reuse
     if (out_aims_json_log) then
       call calc_stats_array(batch_sizes, batch_sizes_stats, n_fractiles = 10, &
                             n_bins = 10)
     end if

     ! create the batch tasking list
     if(use_local_index) then
        call distribute_batch_tasks_by_location( batch_sizes, batch_coords )
     else
        call distribute_batch_tasks( batch_sizes )
     endif
     ! create the batches
     time0 = mpi_wtime()
     if(myid .eq. 0) print*, "=== call create_grid_batches, time_now=", time0
     call create_grid_batches( )
     if(myid .eq. 0) print*, "=== finished create_grid_batches, time_used=", mpi_wtime() - time0
     grid_partitioned = .true.

     if (allocated(all_coords)) then
        call aims_deallocate(all_coords, "all_coords")
     end if
     if (allocated(batch_sizes)) then
        call aims_deallocate(batch_sizes, "batch_sizes")
     end if
     if (allocated(batch_coords)) then
        call aims_deallocate(batch_coords, "batch_coords")
     end if
     if (allocated(grid_partition)) then
        call aims_deallocate(grid_partition, "grid_partition")
     end if

     return
  end if
  
  ! partition grid only on thread 0, the arrays are allocated only on thread 0
  if (myid == 0) then
     allocate(all_coords(3,n_points_in_grid),stat=info)
     call check_allocation(info, 'all_coords                    ')

     root%size = n_points_in_grid
     allocate( root%point_indexes(root%size),stat=info ) 
     call check_allocation(info, 'root%point_indexes            ')
     

     i_point = 0
   
     do i_atom = 1, n_atoms, 1
        do i_radial = 1, n_radial(species(i_atom)), 1
           do i_angular = 1, n_angular( i_radial,species(i_atom) )
              
              i_point = i_point + 1
              
              coord_current(:) = coords_center( :,i_atom ) + &
                   r_angular(:, i_angular, i_radial, species(i_atom)) * &
                   r_radial(i_radial, species(i_atom))
              
              all_coords(:,i_point) = coord_current(:)
              root%point_indexes(i_point) = i_point
              
           end do
        end do
     end do
  end if

  ! allocate array for the grid partition
  allocate(grid_partition(n_points_in_grid),stat=info)
  call check_allocation(info, 'grid_partition                ')


  ! allocate array for the weights for the points
  allocate(point_weight(n_points_in_grid),stat=info)
  call check_allocation(info, 'point_weight                  ')


  time2 = mpi_wtime()

  ! compute the weights for the partitioning
  n_max_compute_part = 0
  point_weight = 0.0d0
  i_index = 0
  do i_atom = 1, n_atoms, 1
     do i_radial = 1, n_radial(species(i_atom)), 1
        if (myid == radial_task_list( i_radial,i_atom )) then

           do i_angular = 1, n_angular( i_radial,species(i_atom) )

              i_index = i_index + 1
              coord_current(:) = coords_center( :,i_atom ) + &
                   r_angular(:, i_angular, i_radial, species(i_atom)) * &
                   r_radial(i_radial, species(i_atom))
              
              n_compute_c = 0
              n_compute_a = 0
              i_basis = 0

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if
              
              ! compute atom-centered coordinates of current integration point,
              ! as viewed from all atoms
              call tab_atom_centered_coords_p0 &
                   ( coord_current,  &
                   dist_tab_sq,  &
                   dir_tab, &
                   n_centers_integrals, centers_basis_integrals )
              
              ! determine which basis functions are relevant at current integration point,
              ! and tabulate their indices
              call prune_basis_p0 &
                   ( dist_tab_sq, &
                   n_compute_a, n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, &
                   inv_centers_basis_integrals )

              ! store the relaitve weight of current integration point
              n_max_compute_part = MAX(n_max_compute_part, n_compute_a)
!              point_weight(my_indexes(i_my_point)) = dble(n_compute_a)
              point_weight(i_index) = dble(n_compute_a)
             
           end do
        else
           i_index = i_index + n_angular( i_radial, species(i_atom) )
        end if
     end do
  end do

  if(myid .eq. 0) print*, "myid=", myid, " partition_grid:part2=", mpi_wtime()-time2

  ! normalize the weights
  call sync_find_max(n_max_compute_part, max_n_max_compute)
  call sync_vector( point_weight, n_points_in_grid )
  point_weight = point_weight / dble(max_n_max_compute)

  !  point_weight = 1.0d0

  ! partition the grid only on the master thread
  if (myid == 0) then
     select case(grid_partitioning_method)

     ! octree
     case(1)
        write (info_str,'(2X,A)') &
             "Partitioning the integration grid into batches with octree method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        write (info_str,'(2X,A,F9.3)') &
             "| Maximal weight for a single point: ", MAXVAL(point_weight)
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        write (info_str,'(2X,A,F9.3)') &
             "| Minimal weight for a single point: ", MINVAL(point_weight)
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        n_grid_batches = 0
        call build_partition_tree( root )
        call remove_small_batches()

     ! versions of METIS
     case(2)

        n_grid_batches = CEILING(dble(n_points_in_grid)/dble(n_points_in_batch))
        points_out = n_points_in_grid
        grid_partition = 0

        if (use_tetgen) then
           write (info_str,'(2X,A)') &
                "Partitioning the integration grid into batches with TetGen and METIS."
           call localorb_info(info_str,use_unit,'(A)',OL_high)
           call metis_tetgen_wrapper(all_coords(1,:), all_coords(2,:), all_coords(3,:), &
                grid_partition, points_out, n_grid_batches)

        else if (use_qhull) then
           write (info_str,'(2X,A)') &
                "Partitioning the integration grid into batches with qhull and METIS."
           call localorb_info(info_str,use_unit,'(A)',OL_high)
           call metis_qhull_wrapper(all_coords(1,:), all_coords(2,:), all_coords(3,:), &
                grid_partition, points_out, n_grid_batches)

        else if (use_nearest_neighbours) then
           write (info_str,'(2X,A)') &
                "Partitioning the integration grid into batches using nearest neighbours and METIS."
           call localorb_info(info_str,use_unit,'(A)',OL_high)

           call build_nn_adjacency_graph()
           call metis_nearest_wrapper(adj_rows, adjacency, grid_partition, points_out, &
                n_grid_batches)
           
        else
           write (info_str,'(1X,A)') &
                "* Unknown graph generation method specified."
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           stop
        end if

        ! check if the meshing succeeded or if there were points that got lost
        if (points_out /= n_points_in_grid) then
           write (info_str,'(1X,A)') &
                "* METIS grid partitioning failed. Patching."
           call localorb_info(info_str,use_unit,'(A)',OL_high)
           ! if so, try patching the partition
           call patch_partition()
        end if

        ! finally, check if the partition has legal values, if not, switch to octree
        if ((MINVAL(grid_partition,1) < 1) &
             .or.(MAXVAL(grid_partition,1) > n_grid_batches)) then

           write (info_str,'(1X,A)') &
                "* Patching METIS grid partitioning failed. Switching to octree."
           call localorb_info(info_str,use_unit,'(A)',OL_high)

           ! if so, switch back to octree as a grid partitioning method
           write (info_str,'(2X,A)') &
                "Partitioning the integration grid into batches with octree method."
           call localorb_info(info_str,use_unit,'(A)',OL_high)
           
           write (info_str,'(2X,A,F9.3)') &
                "| Maximal weight for a single point: ", MAXVAL(point_weight)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           
           write (info_str,'(2X,A,F9.3)') &
                "| Minimal weight for a single point: ", MINVAL(point_weight)
           call localorb_info(info_str,use_unit,'(A)',OL_norm)
           
           n_grid_batches = 0
           grid_partition = 0
           call build_partition_tree( root )
        end if

     ! octants
     case(3)

        write (info_str,'(2X,A)') &
             "Partitioning the integration grid into batches using octants."
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        n_grid_batches = 0
        grid_partition = 0
        i_point = 0
        do i_atom = 1, n_atoms, 1
           do i_radial = 1, n_radial(species(i_atom)), 1

              call partition_angular_shell( r_angular(:, :, i_radial, species(i_atom)), &
                   n_angular( i_radial,species(i_atom) ), i_point )

           end do
        end do

     ! grouping
     case(4)

        write (info_str,'(2X,A)') &
             "Partitioning the integration grid into batches using grouping method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)
        write (info_str,'(2X,A,I6)') &
             "Using grouping factor ", grouping_factor
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        allocate(base_group(n_points_in_grid),stat=info)
        call check_allocation(info, 'base_group                    ')

        i_point = 0
        do i_atom = 1, n_atoms, 1
           do i_radial = 1, n_radial(species(i_atom)), 1
              do i_angular = 1, n_angular( i_radial,species(i_atom) )
                 i_point = i_point + 1
                 base_group(i_point)%level = 1
                 base_group(i_point)%size = 1
                 allocate(base_group(i_point)%point_indexes(base_group(i_point)%size),stat=info)
                 call check_allocation(info, 'base_groupent_basis_wave      ')

                 base_group(i_point)%point_indexes(1) = i_point
              end do
           end do
        end do

        call merge_groups(base_group, n_points_in_grid)

        do i_point = 1, n_points_in_grid, 1
           deallocate(base_group(i_point)%point_indexes)
        end do
        deallocate(base_group)

     ! maxmin
     case(5)
        write (info_str,'(2X,A)') &
             "Partitioning the integration grid into batches with maxmin method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)

        write (info_str,'(2X,A,F9.3)') &
             "| Maximal weight for a single point: ", MAXVAL(point_weight)
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        write (info_str,'(2X,A,F9.3)') &
             "| Minimal weight for a single point: ", MINVAL(point_weight)
        call localorb_info(info_str,use_unit,'(A)',OL_norm)

        n_grid_batches = 0
        call build_maxmin_partition_tree( root )
!        call remove_small_batches()

     case default
        write (info_str,'(1X,A)') &
             "* Unknown grid partitioning method."
        call localorb_info(info_str,use_unit,'(A)',OL_high)
        stop

     end select
  end if

  if (allocated(all_coords)) then
     deallocate(all_coords)
  end if
  if (allocated(root%point_indexes)) then
     deallocate(root%point_indexes)
  end if

  ! broadcast the partition to other threads
  call sync_grid_partition( grid_partition, n_grid_batches, n_points_in_grid )

  ! allocate the batch data size array
  allocate(batch_sizes(n_grid_batches),stat=info)
  call check_allocation(info, 'batch_sizes                   ')


  ! compute the sizes of the batches according to the distribution scheme
  batch_sizes = 0
  do i_point = 1, n_points_in_grid, 1
     batch_sizes(grid_partition(i_point)) = batch_sizes(grid_partition(i_point)) + 1
  end do

  ! print some information on the resulting batches
  call print_out_batch_report( )

  ! While we still have the batch sizes in memory, calculate their statistics
  ! for later reuse
  if (out_aims_json_log) then
     call calc_stats_array(batch_sizes, batch_sizes_stats, n_fractiles = 10, &
                           n_bins = 10)
  end if

  ! create the batch tasking list
  if (use_metis_batch_distribution) then
     call distribute_batch_tasks_metis( batch_sizes, grid_partition, n_points_in_grid )
  else
     call distribute_batch_tasks( batch_sizes )
  end if

  ! create the batches
  call create_grid_batches( )

  grid_partitioned = .true.

  if (allocated(batch_sizes)) then
     deallocate(batch_sizes)
  end if
  if (allocated(grid_partition)) then
     deallocate(grid_partition)
  end if
  if (allocated(point_weight)) then
     deallocate(point_weight)
  end if


!!$  do i_batch = 1, n_my_batches, 1
!!$     if (i_batch == 8) then
!!$     print *, 'step1: i_batch: ', i_batch, 'size: ', batches(i_batch) % size
!!$     do i_point = 1, batches(i_batch) % size, 1
!!$        print *, batches(i_batch)%points(i_point)%coords(1), &
!!$             batches(i_batch)%points(i_point)%coords(2),batches(i_batch)%points(i_point)%coords(3)
!!$        print *, 'index_atom: ', batches(i_batch)%points(i_point)%index_atom
!!$     end do
!!$     end if
!!$  end do
!!$  stop


  contains
    !******
!----------------------------------------------------------------------------------
!****s* partition_grid/create_grid_batches
!  NAME
!    create_grid_batches
!  SYNOPSIS

    subroutine create_grid_batches()

!  PURPOSE
!    Creates the batches for each thread. At this point the following variables are set:
!    o n_grid_batches - total number of batches in the grid
!    o n_my_batches - the number of batches the thread has
!    o batch_sizes - size of each batch
!    o batch_task_list - distribution of the batches over the threads
!    o grid_partition - partition of the points of the grid into batches
!    In the end the array batches is created that stores all the information threadwise
!  USES

      implicit none
!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE


      integer :: i_my_batch, i_batch, i_point
      integer :: i_atom, i_radial, i_angular, info, i_atom_o
      integer :: current_batch, my_current_batch
      integer :: n_atom_points, max_atom_points, i_my_point, i_offset
      integer :: n_atoms_tile_size = 128
      integer, allocatable :: global_grid_partition(:)
      integer, allocatable :: global_grid_partition_tile(:)
      integer, allocatable :: i_offset_tile(:)
      integer, dimension(:), allocatable :: batch_point
      integer, dimension(:), allocatable :: batch_index

      allocate(batches(n_my_batches),stat=info)
      call check_allocation(info, 'batches                       ')

      allocate(batch_index(n_grid_batches),stat=info)
      call check_allocation(info, 'batch_index                   ')


      batch_index = 0

      ! distribute the the batches to threads according to the task list
      ! and allocate the storage space for each thread
      i_my_batch = 0
      do i_batch = 1, n_grid_batches, 1
         if (myid.eq.batch_task_list(i_batch)) then
            i_my_batch = i_my_batch + 1
            batches(i_my_batch) % size = batch_sizes(i_batch)
            batch_index(i_batch) = i_my_batch
            if (batch_sizes(i_batch).gt.0) then
               allocate( batches(i_my_batch) % points(batch_sizes(i_batch)),stat=info )
               call check_allocation(info, 'batches                       ')

            end if
         end if
      end do

      allocate(batch_point(n_my_batches),stat=info)
      call check_allocation(info, 'batch_point                   ')

      batch_point = 0

      ! Get the maximum number of points per atom

      max_atom_points = 0
      do i_atom = 1, n_atoms, 1
         n_atom_points = 0
         do i_radial = 1, n_radial(species(i_atom)), 1
            n_atom_points = n_atom_points + n_angular( i_radial,species(i_atom) )
         enddo
         max_atom_points = max(max_atom_points,n_atom_points)
      enddo

      allocate(global_grid_partition(max_atom_points),stat=info)
      call check_allocation(info, 'global_grid_partition         ')
      allocate(global_grid_partition_tile(max_atom_points * n_atoms_tile_size),stat=info)
      call check_allocation(info, 'global_grid_partition_tile    ')
      allocate(i_offset_tile(n_atoms_tile_size),stat=info)
      call check_allocation(info, 'i_offset_tile         ')

      if(myid .eq. 0) print*, "max_atom_points=", max_atom_points

      ! store all required grid information to the thread that is treating the batch
      i_point = 0
      i_my_point = 0
      do i_atom_o = 1, n_atoms, n_atoms_tile_size
         if (parallel_grid_partitioning) then
            global_grid_partition_tile = 0
         end if

         i_offset = i_point
         do i_atom = i_atom_o, min(n_atoms, i_atom_o+n_atoms_tile_size-1), 1
            ! i_offset_tile(i_atom - i_atom_o + 1) = i_offset

            if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then

               ! global_grid_partition = 0
               do i_radial = 1, n_radial(species(i_atom)), 1
                  do i_angular = 1, n_angular( i_radial,species(i_atom) )
                     i_point = i_point + 1
                     if (myid == MOD(i_point, n_tasks)) then
                        i_my_point = i_my_point + 1
                        ! global_grid_partition(i_point-i_offset) = grid_partition(i_my_point)
                        global_grid_partition_tile(i_point-i_offset) = grid_partition(i_my_point)
                     endif
                  end do
               end do
               ! call sync_integer_vector(global_grid_partition,i_point-i_offset)

            elseif (parallel_grid_partitioning.and.use_hashed_maxmin) then
               ! global_grid_partition = 0

               do i_radial = 1, n_radial(species(i_atom)), 1
                  do i_angular = 1, n_angular( i_radial,species(i_atom) )
                     i_point = i_point + 1

                     coord_current(:) = coords_center( :,i_atom ) + &
                          r_angular(:, i_angular, i_radial, species(i_atom)) * &
                          r_radial(i_radial, species(i_atom))

                    coord_current_temp = coord_current 

                    if(n_periodic > 0)then
                       call map_to_center_cell(coord_current_temp(1:3) )
                    end if

                    if (use_integer_hash) then

                       coord_current_int = INT(coord_current_temp*hash_multiplier)

                       do i_index = 1, 3, 1
                          remainder(i_index) = MODULO(coord_current_int(i_index) - coord_min_int(i_index), h_int(i_index))
                          idx(i_index) = &
                               (coord_current_int(i_index) - coord_min_int(i_index) - remainder(i_index)) / h_int(i_index) + 1
                       end do

                    else

                       idx = CEILING((coord_current_temp - coord_min)/h)

                    end if
                 
                    idx = MIN(idx, n_bins)
                    idx = MAX(idx, 1)

                    hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*n_bins(2)*n_bins(1)

                    if (hash_value == myid+1) then
                       i_my_point = i_my_point + 1
                       ! global_grid_partition(i_point-i_offset) = grid_partition(i_my_point)
                       global_grid_partition_tile(i_point-i_offset) = grid_partition(i_my_point)
                    endif
                 end do
              end do
              ! call sync_integer_vector(global_grid_partition,i_point-i_offset)

            else

               do i_radial = 1, n_radial(species(i_atom)), 1
                  do i_angular = 1, n_angular( i_radial,species(i_atom) )
                     i_point = i_point + 1
                     global_grid_partition_tile(i_point-i_offset) = grid_partition(i_point)
                  end do
               end do

            endif
         
         end do ! end do i_atom = i_atom_o, min(n_atoms, i_atom_o+n_atoms_tile_size-1), 1

         if (parallel_grid_partitioning) then
            call sync_integer_vector(global_grid_partition_tile, i_point-i_offset)
         end if

         i_point = i_offset ! reset it for following loop
         do i_atom = i_atom_o, min(n_atoms, i_atom_o+n_atoms_tile_size-1), 1
            
            ! i_point = i_offset_tile(i_atom - i_atom_o + 1)

            do i_radial = 1, n_radial(species(i_atom)), 1
               do i_angular = 1, n_angular( i_radial,species(i_atom) )
                  i_point = i_point + 1
                  ! current_batch = global_grid_partition(i_point-i_offset)
                  current_batch = global_grid_partition_tile(i_point-i_offset)
                  !
                  !  VB: There appears to be a possibility that, during parallel grid partitioning,
                  !      a grid point does not get assigned to any batch of grid points.
                  !      A way for this to happen is a tiny deviation between the values of
                  !      the system lattice vectors (usually seen for unit cell optimizations)
                  !      or coordinates as stored by different MPI tasks.
                  !      The most frequent occurrence is due to an off-by-one-bit discrepancy
                  !      which is tiny and evidently related to specific compiler versions of
                  !      ifort. In practice, this particular reason is harmless.
                  !      We still stop a run in case this happens but it is possible to 
                  !      override the stop, as outlined below.
                  !
                  if (current_batch.eq.0) then
                     ! must use explicit write statements, not localorb_info.
                     ! localorb_info will only write output for MPI task myid=0, not for any others.
                     write(use_unit,*) '*** Error creating the grid point partitions on task myid ', myid, ': Point ', i_point, 'not found.'
                     write(use_unit,*) "*** "
                     write(use_unit,*) "*** In our experience, this may indicate a difference between the values"
                     write(use_unit,*) "*** for the system geometry (coords and lattice_vectors) stored on different"
                     write(use_unit,*) "*** MPI tasks."
                     write(use_unit,*) "*** Obviously, this should not happen. However, the most common occurrence"
                     write(use_unit,*) "*** of this phenomenon appears to he harmless."
                     write(use_unit,*) "*** It seems that certain compilers can produce off-by-one-bit"
                     write(use_unit,*) "*** geometry values between different MPI tasks (processors)."
                     write(use_unit,*) "*** "
                     write(use_unit,*) "*** A possible solution is to set the flag 'check_cpu_consistency'"
                     write(use_unit,*) "*** in control.in. In this case, the code will simply update the"
                     write(use_unit,*) "*** geometry to be the same geometry on all CPU tasks."
                     write(use_unit,*) "*** "
                     write(use_unit,*) "*** This flag could be set by default, at least for a very small discrepancy."
                     write(use_unit,*) "*** We do not do so in order to avoid masking any potential real discrepancies "
                     write(use_unit,*) "*** that could arise in the future."
                     write(use_unit,*) "*** It you set the 'check_cpu_consistency' flag, please check the output file"
                     write(use_unit,*) "*** to ensure that the magnitude of the inconsistency is indeed practically zero."
                     call aims_stop("Grid point not in any batch","create_grid_batches in partition_grid")
                  end if
                  ! test end
                  my_current_batch = batch_index(current_batch)
                  coord_current(:) = coords_center( :,i_atom ) + &
                       r_angular(:, i_angular, i_radial, species(i_atom)) * &
                       r_radial(i_radial, species(i_atom))
                  if (my_current_batch /= 0) then
                     batch_point(my_current_batch) = batch_point(my_current_batch) + 1
                     batches(my_current_batch)%points(batch_point(my_current_batch))%coords(:) = &
                          coord_current(:)
                     batches(my_current_batch)%points(batch_point(my_current_batch))%index_atom = &
                          i_atom
                     batches(my_current_batch)%points(batch_point(my_current_batch))%index_radial = &
                          i_radial
                     batches(my_current_batch)%points(batch_point(my_current_batch))%index_angular = &
                          i_angular
                  end if
               end do ! end do i_angular = 1, n_angular( i_radial,species(i_atom) )
            end do ! end do i_radial = 1, n_radial(species(i_atom)), 1
         end do ! end do i_atom = i_atom_o, min(n_atoms, i_atom_o+n_atoms_tile_size-1), 1
      end do ! end do i_atom = 1, n_atoms, 1

      deallocate(batch_index)
      deallocate(batch_point)
      deallocate(global_grid_partition)
      deallocate(global_grid_partition_tile)
      deallocate(i_offset_tile)
      
    end subroutine create_grid_batches
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/print_out_batch_report
!  NAME
!    print_out_batch_report
!  SYNOPSIS
    subroutine print_out_batch_report()

!  PURPOSE
!    Prints out a report on the batches created by the grid partitiong routine
!  USES

      implicit none
!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
! SOURCE


      real*8 :: avg_batch_size, batch_size_variance

      ! print some information on the resulting batches
      write (info_str,'(2X,A,I7)') &
           "| Number of batches:  ", n_grid_batches
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      n_max_batch_size = MAXVAL(batch_sizes) 
      write (info_str,'(2X,A,I7)') &
           "| Maximal batch size: ", n_max_batch_size
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      write (info_str,'(2X,A,I7)') &
           "| Minimal batch size: ", MINVAL(batch_sizes)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      avg_batch_size = dble(n_points_in_grid)/dble(n_grid_batches)
      write (info_str,'(2X,A,F11.3)') &
           "| Average batch size: ", avg_batch_size
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      batch_size_variance = 0.0d0
      do i_batch = 1, n_grid_batches, 1
         batch_size_variance = batch_size_variance + &
              (batch_sizes(i_batch) - avg_batch_size)**2
      end do
      batch_size_variance = batch_size_variance / dble(n_grid_batches)
      write (info_str,'(2X,A,F11.3)') &
           "| Standard deviation of batch sizes: ", sqrt(batch_size_variance)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
    end subroutine print_out_batch_report
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/build_partition_tree
!  NAME
!    build_partition_tree
!  SYNOPSIS

    recursive subroutine build_partition_tree( node )

!  PURPOSE
!    Computes a partition of the grid points using a recursive octree method.
!    As a result the array grid_partition is filled with the batch indexes.
!  USES

      implicit none

!  ARGUMENTS

      type(tree_node) :: node

!  INPUTS
!    o node -- a node in the partition octree of the grid points
!  OUTPUT
!    none
!  SOURCE


      integer :: i_index
      real*8, dimension(3) :: center_of_mass
      real*8, dimension(3) :: coord_current
      real*8 :: max_point_weight
      integer, dimension(:), allocatable :: size_children
      integer, dimension(:), allocatable :: i_point
      logical, dimension(3) :: side_coord
      integer :: total_node_size
      real*8 :: total_max_point_weight

      max_point_weight = 0.0d0
      do i_index = 1, node%size, 1
         max_point_weight = MAX(max_point_weight, &
              point_weight(node%point_indexes(i_index)))
      end do

      if (parallel_grid_partitioning) then
         call get_max_double( total_max_point_weight, max_point_weight)
      else
         total_max_point_weight = max_point_weight
      end if

      total_node_size = node%size
      if (parallel_grid_partitioning) then
         call sync_one_integer( total_node_size )
      end if

      ! if the current node is too large
      if ((total_node_size*total_max_point_weight > n_points_in_batch) .or. &
           (total_node_size > batch_size_hard_limit)) then

         ! compute the center of mass of the current node
         center_of_mass = 0.0d0
         do i_index = 1, node%size, 1
            coord_current(:) = all_coords(:,node%point_indexes(i_index))
            if(n_periodic > 0)then
               call map_to_center_cell( coord_current(1:3) )
            end if
            center_of_mass = center_of_mass + coord_current
         end do

         if (parallel_grid_partitioning) then
            call sync_vector( center_of_mass, 3 )
         end if

         center_of_mass = center_of_mass / DBLE(total_node_size)

         node % n_children = 8
         allocate( node%children(node % n_children),stat=info )
         call check_allocation(info, ' node%children                ')


         allocate( size_children(node % n_children),stat=info )
         call check_allocation(info, 'size_children                 ')
      

         allocate( i_point(node % n_children) ,stat=info)
         call check_allocation(info, 'i_point                       ')

         ! compute the sizes of the child nodes
         size_children = 0
         do i_index = 1, node%size, 1

            coord_current(:) = all_coords(:,node%point_indexes(i_index))
            if (n_periodic > 0) then
               call map_to_center_cell( coord_current(1:3) )
            end if
            side_coord(:) = ( coord_current(:).gt.center_of_mass(:) )

            if ( side_coord(1).and.side_coord(2).and.side_coord(3) ) then
               size_children(1) = size_children(1) + 1
            elseif ( side_coord(1).and.side_coord(2).and.(.not.side_coord(3)) ) then
               size_children(2) = size_children(2) + 1
            elseif ( side_coord(1).and.(.not.side_coord(2)).and.side_coord(3) ) then
               size_children(3) = size_children(3) + 1
            elseif ( (.not.side_coord(1)).and.side_coord(2).and.side_coord(3) ) then
               size_children(4) = size_children(4) + 1
            elseif ( side_coord(1).and.(.not.side_coord(2)).and.(.not.side_coord(3)) ) then
               size_children(5) = size_children(5) + 1
            elseif ( (.not.side_coord(1)).and.side_coord(2).and.(.not.side_coord(3)) ) then
               size_children(6) = size_children(6) + 1
            elseif ( (.not.side_coord(1)).and.(.not.side_coord(2)).and.side_coord(3) ) then
               size_children(7) = size_children(7) + 1
            else
               size_children(8) = size_children(8) + 1
            end if

         end do

         do i_index = 1, node % n_children, 1
            node % children(i_index) % size = size_children(i_index)
            if (size_children(i_index) > 0) then
               allocate( node % children(i_index) % point_indexes(size_children(i_index)),stat=info)
               call check_allocation(info, 'node % children               ')

               if (parallel_grid_partitioning) then
                  allocate( node % children(i_index) % global_point_indexes(size_children(i_index)),stat=info)
                  call check_allocation(info, 'node % children               ')

               end if
            end if
         end do

         ! create the child nodes by geometric split using the center of mass as origin
         i_point = 0
         do i_index = 1, node%size, 1

            coord_current(:) = all_coords(:,node%point_indexes(i_index))
            if (n_periodic > 0) then
               call map_to_center_cell( coord_current(1:3) )
            end if
            side_coord(:) = ( coord_current(:).gt.center_of_mass(:) )
            
            if ( side_coord(1).and.side_coord(2).and.side_coord(3) ) then
               i_point(1) = i_point(1) + 1
               node % children(1) % point_indexes(i_point(1)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(1) % global_point_indexes(i_point(1)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( side_coord(1).and.side_coord(2).and.(.not.side_coord(3)) ) then
               i_point(2) = i_point(2) + 1
               node % children(2) % point_indexes(i_point(2)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(2) % global_point_indexes(i_point(2)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( side_coord(1).and.(.not.side_coord(2)).and.side_coord(3) ) then
               i_point(3) = i_point(3) + 1
               node % children(3) % point_indexes(i_point(3)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(3) % global_point_indexes(i_point(3)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( (.not.side_coord(1)).and.side_coord(2).and.side_coord(3) ) then
               i_point(4) = i_point(4) + 1
               node % children(4) % point_indexes(i_point(4)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(4) % global_point_indexes(i_point(4)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( side_coord(1).and.(.not.side_coord(2)).and.(.not.side_coord(3)) ) then
               i_point(5) = i_point(5) + 1
               node % children(5) % point_indexes(i_point(5)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(5) % global_point_indexes(i_point(5)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( (.not.side_coord(1)).and.side_coord(2).and.(.not.side_coord(3)) ) then
               i_point(6) = i_point(6) + 1
               node % children(6) % point_indexes(i_point(6)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(6) % global_point_indexes(i_point(6)) = &
                       node % global_point_indexes(i_index)
               end if
            elseif ( (.not.side_coord(1)).and.(.not.side_coord(2)).and.side_coord(3) ) then
               i_point(7) = i_point(7) + 1
               node % children(7) % point_indexes(i_point(7)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(7) % global_point_indexes(i_point(7)) = &
                       node % global_point_indexes(i_index)
               end if
            else
               i_point(8) = i_point(8) + 1
               node % children(8) % point_indexes(i_point(8)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(8) % global_point_indexes(i_point(8)) = &
                       node % global_point_indexes(i_index)
               end if
            end if

         end do

         ! for each child, do a recursive split
         do i_index = 1, node % n_children, 1
            call build_partition_tree( node % children(i_index) )
         end do

         do i_index = 1, node % n_children, 1
            if (allocated( node%children(i_index)%point_indexes )) then
               deallocate( node%children(i_index)%point_indexes )
            end if
            if (allocated(node%children(i_index)%global_point_indexes)) then
               deallocate( node%children(i_index)%global_point_indexes )
            end if
         end do
         deallocate( node%children )
         deallocate( size_children )
         deallocate( i_point )

      else

         ! if the node is small enough (and > 0) then add it as a batch to the system
         if (total_node_size > 0) then
            n_grid_batches = n_grid_batches + 1
            do i_index = 1, node%size, 1
               grid_partition(node%point_indexes(i_index)) = n_grid_batches
            end do
         end if

      end if
    end subroutine build_partition_tree
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/build_maxmin_partition_tree
!  NAME
!    build_maxmin_partition_tree
!  SYNOPSIS
    recursive subroutine build_maxmin_partition_tree( node )
!  PURPOSE
!    Computes a partition of the grid points using a recursive adapted cut-plane method.
!    As a result the array grid_partition is filled with the batch indexes.
!  USES
      implicit none
!  ARGUMENTS
      type(tree_node) :: node
!  INPUTS
!    o node -- a node in the partition binary tree of the grid points
!  OUTPUT
!    none
! SOURCE


      integer :: i_index, i_index_2
      real*8, dimension(3) :: center_of_mass
      real*8, dimension(3) :: coord_current
      real*8 :: max_point_weight

      integer, dimension(:), allocatable :: size_children
      integer, dimension(:), allocatable :: i_point
      logical, dimension(3) :: side_coord

      real*8, dimension(3) :: normal
      real*8, dimension(3,3) :: matrix
      real*8, allocatable, dimension(:,:) :: coords_vector
      real*8, allocatable, dimension(:) :: vector
      real*8, allocatable, dimension(:) :: temp_vector
      real*8 :: median

      real*8, dimension(3) :: eigenvalues
      integer, parameter :: lwork = 8
      real*8 , dimension(lwork) :: work
      integer :: info
      integer :: total_node_size
      real*8 :: total_max_point_weight

      real*8, external :: ddot

      max_point_weight = 0.0d0
      do i_index = 1, node%size, 1
         max_point_weight = MAX(max_point_weight, &
              point_weight(node%point_indexes(i_index)))
      end do

      if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
         call get_max_double( total_max_point_weight, max_point_weight)
      else
         total_max_point_weight = max_point_weight
      end if

      total_node_size = node%size
      if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
         call sync_one_integer( total_node_size )
      end if

      ! if the current node is too large
      if ((total_node_size*total_max_point_weight > n_points_in_batch) .or. &
           (total_node_size > batch_size_hard_limit)) then

         ! compute the center of mass of the current node
         center_of_mass = 0.0d0
         do i_index = 1, node%size, 1
            coord_current(:) = all_coords(:,node%point_indexes(i_index))
            if(n_periodic > 0)then
               call map_to_center_cell( coord_current(1:3) )
            end if
            center_of_mass = center_of_mass + coord_current
         end do

         if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
            call sync_vector( center_of_mass, 3 )
         end if

         center_of_mass = center_of_mass / DBLE(total_node_size)

         node % n_children = 2
         allocate( node%children(node % n_children) ,stat=info)
         call check_allocation(info, 'node%children                 ')

         allocate( size_children(node % n_children),stat=info )
         call check_allocation(info, 'size_children                 ')

         allocate( i_point(node % n_children),stat=info )
         call check_allocation(info, ' i_point                      ')

         ! compute the normal of the spltting plane
         matrix = 0.0d0

         if (node%size > 0) then
            allocate(coords_vector(node%size,3),stat=info)
            call check_allocation(info, 'coords_vector                 ')

            do i_index = 1, node%size, 1
               coord_current(:) = all_coords(:,node%point_indexes(i_index))
               if(n_periodic > 0)then
                  call map_to_center_cell( coord_current(1:3) )
               end if
               coords_vector(i_index, :) = coord_current(:) - center_of_mass(:)
            end do
            
            do i_index = 1, 3, 1
               do i_index_2 = 1, i_index, 1
                  matrix(i_index_2, i_index) = &
                       ddot( node%size, coords_vector(:,i_index), 1, coords_vector(:,i_index_2), 1 )
               end do
            end do
         end if

         if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
            do i_index = 1, 3, 1
               call sync_vector( matrix(:,i_index), 3 )
            end do
         end if

         call dsyev( 'V', 'U', 3, matrix, 3, eigenvalues, work, lwork, info )
         normal = matrix(:,3)

         if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
!            median = ddot( 3, normal, 1, center_of_mass, 1 )
            median = 0.0d0
            if (node%size > 0) then
               allocate(vector( node%size ),stat=info)
               call check_allocation(info, 'vector                        ')

               allocate(temp_vector( node%size ),stat=info)
               call check_allocation(info, 'temp_vector                   ')

               call DGEMV('N', node%size, 3, 1.0d0, coords_vector, node%size, normal, 1, 0.0d0, vector, 1)
               temp_vector = vector
               call heapsort(temp_vector, node%size)
               if (node%size > 1) then
                  median = temp_vector(INT(dble(node%size)/2.0))
               else
                  median = temp_vector(1)
               end if
            end if
            call sync_real_number(median)
            median = median / dble(n_tasks)
         else
            ! compute n^T x, sort and find the median
            allocate(vector( node%size ),stat=info)
            call check_allocation(info, 'vector                        ')

            allocate(temp_vector( node%size ),stat=info)
            call check_allocation(info, 'temp_vector                   ')

            call DGEMV('N', node%size, 3, 1.0d0, coords_vector, node%size, normal, 1, 0.0d0, vector, 1)

            temp_vector = vector
            call heapsort(temp_vector, node%size)
            median = temp_vector(INT(dble(node%size)/2.0))
         end if

         ! compute the sizes of the child nodes
         size_children = 0
         do i_index = 1, node%size, 1
            if (vector(i_index) < median) then
               size_children(1) = size_children(1) + 1
            else
               size_children(2) = size_children(2) + 1
            end if
         end do

         ! allocate the child nodes
         do i_index = 1, node % n_children, 1
            node % children(i_index) % size = size_children(i_index)
            if (size_children(i_index) > 0) then
               allocate( node % children(i_index) % point_indexes(size_children(i_index)),stat=info)
               call check_allocation(info, 'node % children               ')

               if (parallel_grid_partitioning) then
                  allocate( node % children(i_index) % global_point_indexes(size_children(i_index)),stat=info)
                  call check_allocation(info, 'node % children               ')

               end if
            end if
         end do

         ! create the child nodes by geometric split with the median of the plane coordinates
         i_point = 0
         do i_index = 1, node%size, 1
            if (vector(i_index) < median) then
               i_point(1) = i_point(1) + 1
               node % children(1) % point_indexes(i_point(1)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(1) % global_point_indexes(i_point(1)) = &
                       node % global_point_indexes(i_index)
               end if
            else
               i_point(2) = i_point(2) + 1
               node % children(2) % point_indexes(i_point(2)) = node % point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(2) % global_point_indexes(i_point(2)) = &
                       node % global_point_indexes(i_index)
               end if
            end if
         end do

         if (allocated(vector)) deallocate(vector)
         if (allocated(temp_vector)) deallocate(temp_vector)
         if (allocated(coords_vector)) deallocate(coords_vector)

         ! for each child, do a recursive split
         do i_index = 1, node % n_children, 1
            call build_maxmin_partition_tree( node % children(i_index) )
         end do

         do i_index = 1, node % n_children, 1
            if (allocated( node%children(i_index)%point_indexes )) then
               deallocate( node%children(i_index)%point_indexes )
            end if
            if (allocated(node%children(i_index)%global_point_indexes)) then
               deallocate( node%children(i_index)%global_point_indexes )
            end if
         end do

         deallocate( node%children )
         deallocate( size_children )
         deallocate( i_point )

      else

         ! if the node is small enough (and > 0) then add it as a batch to the system
         if (total_node_size > 0) then
            n_grid_batches = n_grid_batches + 1
            do i_index = 1, node%size, 1
               grid_partition(node%point_indexes(i_index)) = n_grid_batches
            end do
         end if

      end if
    end subroutine build_maxmin_partition_tree
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/remove_small_batches
!  NAME
!    remove_small_batches
!  SYNOPSIS
    recursive subroutine remove_small_batches()
!  PURPOSE 
!    Removes batches that contain too few points from the partition by merging
!    them with a closest batch.
!  USES
!  ARGUMENTS
!  INPUTS
!    none -- array grid_partition must contain a valid partition of the grid
!  OUTPUT
!    none -- as a result the array grid_partition is modified
!  SOURCE


      integer, allocatable, dimension(:) :: batch_sizes
      real*8, allocatable, dimension(:,:) :: batch_com
      real*8, allocatable, dimension(:) :: dist_sq

      integer :: active_batch, new_batch

      integer :: i_point, i_batch

      ! This doesn't work for parallel grid partitioning
      if (parallel_grid_partitioning) return

      allocate(batch_sizes(n_grid_batches),stat=info)
      call check_allocation(info, 'batch_sizes                   ')

      batch_sizes = 0

      ! check if any of the batches is too small
      do i_point = 1, n_points_in_grid, 1
         batch_sizes(grid_partition(i_point)) = batch_sizes(grid_partition(i_point)) + 1
      end do
      if (MINVAL(batch_sizes,1) >= min_batch_size) then
         deallocate(batch_sizes)
         return
      end if

      ! compute the center of mass for each of the batches
      allocate(batch_com(3,n_grid_batches),stat=info)
      call check_allocation(info, 'batch_com                     ')

      batch_com = 0.0d0
      do i_point = 1, n_points_in_grid, 1
         batch_com(:,grid_partition(i_point)) = batch_com(:,grid_partition(i_point)) + &
              all_coords(:,i_point)
      end do
      do i_batch = 1, n_grid_batches, 1
         batch_com(:,i_batch) = batch_com(:,i_batch) / dble(batch_sizes(i_batch))
      end do

      ! find a batch that is too small
      do i_batch = 1, n_grid_batches, 1
         if (batch_sizes(i_batch) < min_batch_size) then
            active_batch = i_batch
            exit
         end if
      end do

      ! find the batch that is closest to the selected batch
      allocate(dist_sq(n_grid_batches),stat=info)
      call check_allocation(info, 'dist_sq                       ')

      dist_sq(:) = (batch_com(1,:) - batch_com(1,active_batch))**2 + &
           (batch_com(2,:) - batch_com(1,active_batch))**2 + &
           (batch_com(3,:) - batch_com(1,active_batch))**2
      dist_sq(active_batch) = 1.0d6
      new_batch = MINLOC(dist_sq,1)

      ! merge the small batch to the close one
      do i_point = 1, n_points_in_grid, 1
         if (grid_partition(i_point) == active_batch) grid_partition(i_point) = new_batch
      end do
      do i_point = 1, n_points_in_grid, 1
         if (grid_partition(i_point) >= active_batch) then
            grid_partition(i_point) = grid_partition(i_point) - 1
         end if
      end do
      n_grid_batches = n_grid_batches - 1

      deallocate(batch_sizes)
      deallocate(batch_com)
      deallocate(dist_sq)

      call remove_small_batches()

    end subroutine remove_small_batches
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/partition_angular_shell
!  NAME
!    partition_angular_shell
!  SYNOPSIS
    subroutine partition_angular_shell(r_angular, n_angular, i_offset)
!  PURPOSE
!    Partitions a single origin-centered angular Lebedev-shell into batches.
!  USES

!  ARGUMENTS
      integer :: n_angular
      real*8 :: r_angular(3,n_angular)
      integer :: i_offset
!  INPUTS
!    o n_angular -- number of points in the angular shell
!    o r_angular -- coordinates of the points in the angular shell
!    o i_offset -- current offset of the array grid_partition
!  OUTPUT
!    none -- as a result the array grid_partition is updated
!  SOURCE


      integer :: n_division
      integer :: local_batch_index(n_angular)

      integer :: i_end
      integer :: i_angular

      if (n_angular.le.110) then
!     smallest grids: only one division
         n_division = 1
         local_batch_index = n_grid_batches + 1
      else if (n_angular.le.302) then
!     divide into 2 blocks, by z-coordinate
         n_division = 2
         do i_angular = 1, n_angular, 1
            if (r_angular(3,i_angular).le.0) then
               local_batch_index(i_angular) = n_grid_batches + 1
             else
                local_batch_index(i_angular) = n_grid_batches + 2
             end if
          enddo
       else if (n_angular.le.770) then
!      divide into 4 blocks, by x- and z-coordinates
          n_division = 4
          do i_angular = 1, n_angular, 1
             if ((r_angular(1,i_angular).le.0).and. &
                  (r_angular(3,i_angular).le.0)) then
                local_batch_index(i_angular) = n_grid_batches + 1
             else if((r_angular(1,i_angular).gt.0).and. &
                  (r_angular(3,i_angular).le.0)) then
               local_batch_index(i_angular) = n_grid_batches + 2
            else if((r_angular(1,i_angular).le.0).and. &
                 (r_angular(3,i_angular).gt.0)) then
               local_batch_index(i_angular) = n_grid_batches + 3
            else
               local_batch_index(i_angular) = n_grid_batches + 4
            end if
         enddo

      else
!     divide into 8 blocks, by x-, y-, and z-coordinates
         n_division = 8
         do i_angular = 1, n_angular, 1
            if ((r_angular(1,i_angular).le.0).and. &
                 (r_angular(2,i_angular).le.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               local_batch_index(i_angular) = n_grid_batches + 1
            else if((r_angular(1,i_angular).gt.0).and. &
                 (r_angular(2,i_angular).le.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               local_batch_index(i_angular) = n_grid_batches + 2
            else if((r_angular(1,i_angular).le.0).and. &
                 (r_angular(2,i_angular).gt.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               local_batch_index(i_angular) = n_grid_batches + 3
            else if((r_angular(1,i_angular).le.0).and. &
                 (r_angular(2,i_angular).le.0).and. &
                 (r_angular(3,i_angular).gt.0)) then
               local_batch_index(i_angular) = n_grid_batches + 4
            else if((r_angular(1,i_angular).gt.0).and. &
                 (r_angular(2,i_angular).gt.0).and. &
                 (r_angular(3,i_angular).le.0)) then
               local_batch_index(i_angular) = n_grid_batches + 5
            else if((r_angular(1,i_angular).gt.0).and. &
                 (r_angular(2,i_angular).le.0).and. &
                 (r_angular(3,i_angular).gt.0)) then
               local_batch_index(i_angular) = n_grid_batches + 6
            else if((r_angular(1,i_angular).le.0).and. &
                   (r_angular(2,i_angular).gt.0).and. &
                   (r_angular(3,i_angular).gt.0)) then
               local_batch_index(i_angular) = n_grid_batches + 7
            else
               local_batch_index(i_angular) = n_grid_batches + 8
            end if
         enddo

      end if

      n_grid_batches = n_grid_batches + n_division
      do i_angular = 1, n_angular, 1
         i_offset = i_offset + 1
         grid_partition(i_offset) = local_batch_index(i_angular)
      end do

    end subroutine partition_angular_shell
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/patch_partition
!  NAME
!    patch_partition
!  SYNOPSIS
    subroutine patch_partition()
!  PURPOSE
!    In a given partition, finds points that don't have a legal index
!    and merges them with nearest legal batch,
!  USES
!  ARGUMENTS
!  INPUTS
!    none -- the array grid_partition must be set
!  OUTPUT
!    none -- the array grid_partition may be modified
!  SOURCE


      real*8, dimension(:,:), allocatable :: batch_com
      real*8, dimension(:), allocatable :: dist_sq
      integer, dimension(:), allocatable :: batch_size
      integer, dimension(:), allocatable :: extra_points

      integer :: i_point, i_extra_point, i_batch
      integer :: n_extra_points

      allocate(batch_com(3,n_grid_batches),stat=info)
      call check_allocation(info, 'batch_com                     ')

      batch_com = 0.0d0
      allocate(batch_size(n_grid_batches),stat=info)
      call check_allocation(info, 'batch_size                    ')
      
      batch_size = 0

      ! count the points having an illegal index
      n_extra_points = 0
      do i_point = 1, n_points_in_grid, 1
         if ((grid_partition(i_point) >= 1).and. &
              (grid_partition(i_point) <= n_grid_batches)) then
            batch_size(grid_partition(i_point)) = &
                 batch_size(grid_partition(i_point))+1
            batch_com(:,grid_partition(i_point)) = &
                 batch_com(:,grid_partition(i_point)) + &
                 all_coords(:,i_point)
         else
            n_extra_points = n_extra_points + 1
         end if
      end do

      ! find the centers of mass for the batches that are ok
      batch_com(1,:) = batch_com(1,:) / dble(batch_size(:))
      batch_com(2,:) = batch_com(2,:) / dble(batch_size(:))
      batch_com(3,:) = batch_com(3,:) / dble(batch_size(:))

      ! if there are any illegal points
      if (n_extra_points > 0) then
         allocate(extra_points(n_extra_points),stat=info)
         call check_allocation(info, 'extra_points                  ')


         i_extra_point = 0
         do i_point = 1, n_points_in_grid, 1
            if ((grid_partition(i_point) < 1).or. &
                 (grid_partition(i_point) > n_grid_batches)) then
               i_extra_point = i_extra_point + 1
               extra_points(i_extra_point) = i_point
            end if
         end do

         ! put the illegal points to the nearest batch
         allocate(dist_sq(n_grid_batches),stat=info)
         call check_allocation(info, 'dist_sq                       ')
       
         
         do i_extra_point = 1, n_extra_points, 1
            dist_sq(:) = (batch_com(1,:) - all_coords(1,extra_points(i_extra_point)))**2 + &
                 (batch_com(2,:) - all_coords(2,extra_points(i_extra_point)))**2 + &
                 (batch_com(3,:) - all_coords(3,extra_points(i_extra_point)))**2
            i_batch = MINLOC(dist_sq,1)
            grid_partition(extra_points(i_extra_point)) = i_batch
         end do

         deallocate(extra_points)
      end if

      deallocate(batch_com)
      deallocate(batch_size)

    end subroutine patch_partition
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/merge_groups
!  NAME
!    merge_groups
!  SYNOPSIS
    recursive subroutine merge_groups( groups, n_groups )
!  PURPOSE
!    In the grouping algorithm, recursively merges batches until they
!    become large enough.
!  USES
!  ARGUMENTS
      integer :: n_groups
      type(point_group), dimension(n_groups) :: groups
!  INPUTS
!    o n_groups -- current number of point groups
!    o groups -- the point groups
!  OUTPUT
!    none -- in the end the array grid partition is set and the groups
!            have become the new grid batches
!  SOURCE


      integer :: max_group_size
      integer :: i_group, i_point, i_offset

      logical, dimension(:), allocatable :: assigned

      type(point_group), pointer, dimension(:) :: new_groups

      real*8, dimension(:,:), allocatable :: group_com
      real*8, dimension(:), allocatable :: group_dist_sq

      integer :: n_new_groups, group_one, group_two
      integer :: i_new_group, new_group_size
      integer :: n_nonassigned, n_ungrouped

      integer, dimension(:), allocatable :: selected_groups

      max_group_size = 0
      do i_group = 1, n_groups, 1
         max_group_size = MAX(max_group_size, groups(i_group)%size)
      end do

      ! if the maximal group doesn't have enough points yet
      if (max_group_size < n_points_in_batch) then

         allocate(assigned(n_groups),stat=info)
         call check_allocation(info, 'assigned                      ')

         assigned = .false.

         allocate(group_com(3,n_groups),stat=info)
         call check_allocation(info, 'group_com                     ')

         group_com = 0.0d0

         ! compute centers of masses for all groups
         do i_group = 1, n_groups, 1
            do i_point = 1, groups(i_group)%size, 1
               group_com(:,i_group) = group_com(:,i_group) + &
                    all_coords(:,groups(i_group)%point_indexes(i_point))
            end do
            group_com(:,i_group) = group_com(:,i_group) / dble(groups(i_group)%size)
         end do
         
         ! find the number of new groups and the number of leftover groups
         n_ungrouped =  MOD(n_groups,grouping_factor)
         n_new_groups = (n_groups - n_ungrouped) / grouping_factor + n_ungrouped
             
         allocate(new_groups(n_new_groups),stat=info)
         call check_allocation(info, 'new_groups                    ')

         allocate(group_dist_sq(n_groups),stat=info)
         call check_allocation(info, 'group_dist_sq                 ')


         i_new_group = 0

         allocate(selected_groups(grouping_factor),stat=info)
         call check_allocation(info, 'selected_groups               ')


         ! loop as long as there are groups to merge
         do while (COUNT(.not.assigned) >= grouping_factor)

            ! find the first group to be assigned to act as a "seed"
            selected_groups = 0
            loop: do i_group = 1, n_groups, 1
               if (.not.assigned(i_group)) then
                  group_one = i_group
                  selected_groups(1) = i_group
                  exit loop
               end if
            end do loop
      
            group_dist_sq(:) = (group_com(1,:) - group_com(1,selected_groups(1)))**2 + &
                 (group_com(2,:) - group_com(2,selected_groups(1)))**2 + &
                 (group_com(3,:) - group_com(3,selected_groups(1)))**2
            assigned(selected_groups(1)) = .true.

            ! find the closest groups to the "seed"
            new_group_size = groups(selected_groups(1))%size
            do i_group = 2, grouping_factor,1
               selected_groups(i_group) = MINLOC(group_dist_sq,1,(.not.assigned))
               assigned(selected_groups(i_group)) = .true.
               new_group_size = new_group_size + groups(selected_groups(i_group))%size
            end do

            i_new_group = i_new_group + 1

            new_groups(i_new_group)%size = new_group_size
            allocate(new_groups(i_new_group)%point_indexes(new_group_size),stat=info)
            call check_allocation(info, 'new_groups                    ')

            ! merge the groups
            i_offset = 1
            do i_group = 1, grouping_factor,1
               new_groups(i_new_group)% &
                    point_indexes(i_offset:i_offset+groups(selected_groups(i_group))%size-1) = &
                    groups(selected_groups(i_group))%point_indexes(:)
               i_offset = i_offset + groups(selected_groups(i_group))%size
            end do

            new_groups(i_new_group)%level = groups(selected_groups(1))%level + 1
           
         end do
         
         selected_groups = 0
         ! just copy the leftover groups to the next level
         do while ( COUNT(.not.assigned) > 0)

            loop2: do i_group = 1, n_groups, 1
               if (.not.assigned(i_group)) then
                  selected_groups(1) = i_group
                  exit loop2
               end if
            end do loop2
            
            new_group_size = groups(selected_groups(1))%size

            i_new_group = i_new_group + 1
            new_groups(i_new_group)%size = new_group_size
            allocate(new_groups(i_new_group)%point_indexes(new_group_size),stat=info)
            call check_allocation(info, 'new_groups                    ')


            new_groups(i_new_group)%point_indexes(1:groups(selected_groups(1))%size) = &
                 groups(selected_groups(1))%point_indexes(:)

            new_groups(i_new_group)%level = groups(selected_groups(1))%level + 1
            assigned(selected_groups(1)) = .true.

         end do
         
         deallocate(group_com)
         deallocate(group_dist_sq)
         deallocate(assigned)
         deallocate(selected_groups)

         ! recursively merge the new groups again
         call merge_groups( new_groups, n_new_groups )

         do i_group = 1, n_new_groups, 1
            deallocate(new_groups(i_group)%point_indexes)
         end do
         deallocate(new_groups)

      else

         ! if the grouping is done, unload the data to grid_partition
         do i_group = 1, n_groups, 1
            do i_point = 1, groups(i_group)%size, 1
               grid_partition(groups(i_group)%point_indexes(i_point)) = i_group
            end do
         end do
         n_grid_batches = n_groups

      end if

    end subroutine merge_groups
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/build_nn_adjacency_graph
!  NAME
!    build_nn_adjacency_graph
!  SYNOPSIS
    subroutine build_nn_adjacency_graph()
!  PURPOSE
!    Builds an adjacency graph of the points in the grid where each node has
!    the index n_nearest_neighbours
!  USES
      implicit none
!  ARGUMENTS
!  INPUTS
!    none -- the array all_coords must contain the coordinates of the grid points
!  OUTPUT
!    none -- in the end the edges of the graph are stored in the array adjacency and
!            the indeces of the rows are in the array adj_rows
!  SOURCE


      real*8, allocatable, dimension(:) :: dist_sq

      integer, allocatable, dimension(:) :: n_edges
      integer, allocatable, dimension(:,:) :: edges
      integer, allocatable, dimension(:) :: edge_index

      integer :: i_point, i_point_2, i_point_3, i_nearest, current_edge
      logical :: already_set

      allocate(dist_sq(n_points_in_grid),stat=info)
      call check_allocation(info, 'dist_sq                       ')
          
      allocate(edges(n_nearest_neighbours, n_points_in_grid),stat=info)
      call check_allocation(info, 'bedges                        ')

      edges = 0
      allocate(n_edges(n_points_in_grid),stat=info)
      call check_allocation(info, 'n_edges                       ')

      n_edges = 0

      do i_point = 1, n_points_in_grid, 1
         dist_sq(:) = (all_coords(1,:) - all_coords(1,i_point))**2 + &
              (all_coords(2,:) - all_coords(2,i_point))**2 + &
              (all_coords(3,:) - all_coords(3,i_point))**2
         dist_sq(i_point) = 1.0d6

         do i_point_2 = 1, n_nearest_neighbours, 1
            i_nearest = MINLOC(dist_sq,1)
            if (i_nearest > i_point) then
               edges(i_point_2, i_point) = i_nearest
               n_edges(i_point) = n_edges(i_point) + 1
               n_edges(i_nearest) = n_edges(i_nearest) + 1
            else
               already_set = .false.
               do i_point_3 = 1, n_nearest_neighbours, 1
                  if (edges(i_point_3,i_nearest) == i_point) already_set = .true.
               end do
               if (.not.already_set) then
                  edges(i_point_2, i_point) = i_nearest
                  n_edges(i_point) = n_edges(i_point) + 1
                  n_edges(i_nearest) = n_edges(i_nearest) + 1
               end if
            end if
            dist_sq(i_nearest) = 1.0d6
         end do

      end do

      deallocate(dist_sq)

      n_total_edges = SUM(n_edges)

      allocate(adj_rows(n_points_in_grid+1),stat=info)
      call check_allocation(info, 'adj_rows p                    ')

      allocate(adjacency(n_total_edges),stat=info)
      call check_allocation(info, 'adjacency                     ')

      allocate(edge_index(n_points_in_grid),stat=info)
      call check_allocation(info, 'edge_index                    ')

      edge_index = 0

      adj_rows(1) = 1
      do i_point = 2, n_points_in_grid+1, 1
         adj_rows(i_point) = adj_rows(i_point-1) + n_edges(i_point-1)
      end do

      deallocate(n_edges)

      do i_point = 1, n_points_in_grid, 1
         do i_point_2 = 1, n_nearest_neighbours, 1
            current_edge = edges(i_point_2, i_point)
            if ( current_edge /= 0) then
               adjacency(adj_rows(i_point) + edge_index(i_point)) = current_edge
               adjacency(adj_rows(current_edge) + edge_index(current_edge)) = i_point
               edge_index(i_point) = edge_index(i_point) + 1
               edge_index(current_edge) = edge_index(current_edge) + 1
            end if
         end do
      end do

      deallocate(edge_index)
      deallocate(edges)

    end subroutine build_nn_adjacency_graph
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/hash_points
!  NAME
!    hash_points
!  SYNOPSIS
    subroutine hash_points( root )
!  PURPOSE
!  USES
      implicit none
!  ARGUMENTS
      type(tree_node) :: root
!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE
      integer, dimension(n_tasks) :: grid_batches_per_task

      integer :: i_atom, i_radial, i_angular, i_point, i_offset
      integer :: hash_value

      if (use_hashed_maxmin) then

         n_grid_batches = 0
         call build_maxmin_partition_tree( root )

         grid_batches_per_task = 0
         grid_batches_per_task(myid+1) = n_grid_batches
         call sync_integer_vector(grid_batches_per_task, n_tasks)

         i_offset = SUM(grid_batches_per_task(1:myid))
         grid_partition = grid_partition + i_offset

         n_grid_batches = SUM(grid_batches_per_task)
         
      else
         
         do i_my_point = 1, n_my_points, 1
            coord_current(:) = all_coords(:,i_my_point)

            if(n_periodic > 0)then
               call map_to_center_cell(coord_current(1:3) )
            end if

            if (use_integer_hash) then

               coord_current_int = INT(coord_current*hash_multiplier)
                    
               do i_index = 1, 3, 1
                  remainder(i_index) = MODULO(coord_current_int(i_index) - coord_min_int(i_index), h_int(i_index))
                  idx(i_index) = &
                       (coord_current_int(i_index) - coord_min_int(i_index) - remainder(i_index)) / h_int(i_index) + 1
               end do

            else

               idx = CEILING((coord_current - coord_min)/h)

            end if

            idx = MIN(idx, n_bins)
            idx = MAX(idx, 1)
            
            grid_partition(i_my_point) = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*n_bins(2)*n_bins(1)
            
         end do

      end if

    end subroutine hash_points
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/find_bb_dimensions
!  NAME
!    find_bb_dimensions
!  SYNOPSIS
    subroutine find_bb_dimensions(coord_max, coord_min, bb_dim)
!  PURPOSE
!    Find the extent of the bounding box for the grid points 
!  USES
      implicit none
!  ARGUMENTS
      real*8, dimension(3) :: coord_max, coord_min, bb_dim
!  INPUTS
!    none
!  OUTPUT
!    o coord_max -- maximal bounding box coordinates in x, y, and z
!    o coord_min -- minimal bounding box coordinates in x, y, and z
!    o bb_dim -- bounding box dimension in x, y, and z
!  SOURCE
      real*8 :: x_max_loc, x_min_loc, y_max_loc, y_min_loc, z_max_loc, z_min_loc
      integer :: i_atom, i_radial, i_angular, i_point
      real*8, dimension(3) :: coord_current

      x_max_loc = -1.0d6
      y_max_loc = -1.0d6
      z_max_loc = -1.0d6
      x_min_loc = 1.0d6
      y_min_loc = 1.0d6
      z_min_loc = 1.0d6

      i_point = 0
      do i_atom = 1, n_atoms, 1
         do i_radial = 1, n_radial(species(i_atom)), 1
            do i_angular = 1, n_angular( i_radial,species(i_atom) )
               i_point = i_point + 1
               if (myid == MOD(i_point, n_tasks)) then
                  coord_current(:) = coords_center( :,i_atom ) + &
                       r_angular(:, i_angular, i_radial, species(i_atom)) * &
                       r_radial(i_radial, species(i_atom))
                  
                  if(n_periodic > 0)then
                     call map_to_center_cell( coord_current(1:3) )
                  end if
                  
                  x_max_loc = MAX(x_max_loc, coord_current(1))
                  x_min_loc = MIN(x_min_loc, coord_current(1))
                  y_max_loc = MAX(y_max_loc, coord_current(2))
                  y_min_loc = MIN(y_min_loc, coord_current(2))
                  z_max_loc = MAX(z_max_loc, coord_current(3))
                  z_min_loc = MIN(z_min_loc, coord_current(3))
               end if
            end do
         end do
      end do
      call get_max_double(coord_max(1), x_max_loc)
      call get_max_double(coord_max(2), y_max_loc)
      call get_max_double(coord_max(3), z_max_loc)
      call get_min_double(coord_min(1), x_min_loc)
      call get_min_double(coord_min(2), y_min_loc)
      call get_min_double(coord_min(3), z_min_loc)
      
      bb_dim = coord_max - coord_min

    end subroutine find_bb_dimensions
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/find_bins
!  NAME
!    find_bins
!  SYNOPSIS
    subroutine find_bins(bb_dim, n_grid_batches, n_bins, h)
!  PURPOSE
!    Computes the number of bins in x, y, and and the corresponding bin dimensions
!  USES
      implicit none
!  ARGUMENTS      
      integer, dimension(3) :: n_bins
      integer :: n_grid_batches
      real*8, dimension(3) :: h, bb_dim
!  INPUTS
!    o bb_dim -- bounding box dimension in x, y, and z
!  OUTPUT
!    o n_grid_batches -- number of grid batches if simple hashing is used
!    o n_bins -- number of bin divisions in x, y, and z
!    o h -- bin dimension in x, y, and z
!  SOURCE
      integer :: n_batches, n_x, n_y, n_z, div1, div2, n_factors
      integer, allocatable, dimension(:) :: prime_factors
      real*8 :: optimal_value, trial_value

      if (use_hashed_maxmin) then
         call find_prime_factors(n_tasks, n_factors)
         allocate(prime_factors(n_factors))
         call find_prime_factors(n_tasks, n_factors, prime_factors)
         n_batches = n_tasks
         n_bins(1) = &
              NINT( EXP(1.0/3.0*LOG((bb_dim(1)/bb_dim(2))*(bb_dim(1)/bb_dim(3))*dble(n_batches))))
         n_bins(2) = &
              NINT( EXP(1.0/3.0*LOG((bb_dim(2)/bb_dim(3))*(bb_dim(2)/bb_dim(1))*dble(n_batches))))
         n_bins(3) = &
              NINT( EXP(1.0/3.0*LOG((bb_dim(3)/bb_dim(2))*(bb_dim(3)/bb_dim(1))*dble(n_batches))))

         if (PRODUCT(n_bins(:)) /= n_tasks) then

            optimal_value = 1.0d9
            do div1 = 0, n_factors
               do div2 = div1, n_factors
                  n_x = PRODUCT(prime_factors(1:div1))
                  n_y = PRODUCT(prime_factors(div1+1:div2))
                  n_z = n_tasks/(n_x*n_y)
                  trial_value = (bb_dim(1)/dble(n_x) - bb_dim(2)/dble(n_y))**2 + &
                       (bb_dim(1)/dble(n_x) - bb_dim(3)/dble(n_z))**2 + &
                       (bb_dim(2)/dble(n_y) - bb_dim(3)/dble(n_z))**2
                  if (trial_value < optimal_value) then
                     optimal_value = trial_value
                     n_bins(1) = n_x
                     n_bins(2) = n_y
                     n_bins(3) = n_z
                  end if
               end do
            end do
         end if

      else

         n_batches = NINT(dble(n_points_in_grid)/dble(n_points_in_batch))
         n_bins(1) = &
              NINT( ((bb_dim(1))/(bb_dim(2))*(bb_dim(1))/(bb_dim(3))*dble(n_batches))**(1.0d0/3.0d0) )
         n_bins(2) = &
              NINT( ((bb_dim(2))/(bb_dim(3))*(bb_dim(2))/(bb_dim(1))*dble(n_batches))**(1.0d0/3.0d0) )
         n_bins(3) = &
              NINT( ((bb_dim(3))/(bb_dim(2))*(bb_dim(3))/(bb_dim(1))*dble(n_batches))**(1.0d0/3.0d0) )
      end if

      n_grid_batches = n_bins(1)*n_bins(2)*n_bins(3)
      
      h = bb_dim / dble(n_bins)

    end subroutine find_bins
!******
!----------------------------------------------------------------------------------
!****s* partition_grid/find_prime_factors
!  NAME
!    find_prime_factors
!  SYNOPSIS
    subroutine find_prime_factors(n, n_factors, prime_factors)
!  PURPOSE
!    Computes the prime factorization of a positive integer
!  USES
      implicit none
!  ARGUMENTS      
      integer :: n, n_factors
      integer, optional :: prime_factors(n_factors)
!  INPUTS
!    o n -- integer to be factorized
!  OUTPUT
!    o n_factors -- number of factors found
!    o prime_factors -- optional, prime factors of the integer n, must be allocated before call
!  SOURCE
      integer :: i, top, n_loc, n_found
      logical :: found

      n_loc = n
      
      if (.not.(PRESENT(prime_factors))) then
         n_factors = 1
         found = .true.
         do while(found)
            found = .false.
            top = FLOOR(SQRT(dble(n_loc)))
            do i = 2, top
               if (MOD(n_loc,i) == 0) then
                  found = .true.
                  n_factors = n_factors + 1
                  n_loc = n_loc/i
                  exit
               end if
            end do
         end do
      else
         n_found = 0
         do while(n_found < n_factors)
            found = .false.
            top = FLOOR(SQRT(dble(n_loc)))
            do i = 2, top
               if (MOD(n_loc,i) == 0) then
                  found = .true.
                  n_found = n_found + 1
                  n_loc = n_loc/i
                  prime_factors(n_found) = i
                  exit
               end if
            end do
            if (.not.found) then
               n_found = n_found + 1
               prime_factors(n_found) = n_loc
            end if
         end do
      end if

    end subroutine find_prime_factors
!******
end subroutine partition_grid
