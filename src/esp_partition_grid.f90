!****s* FHI-aims/esp_partition_grid
!  NAME
!    esp_partition_grid
!  SYNOPSIS

subroutine esp_partition_grid( radius_esp_min_new, radius_esp_max_new, &
                               grid_output,rm,equal_grid,cube_rad_grid,&
           i_esp_equal_grid_x_in, i_esp_equal_grid_y_in, &
           i_esp_equal_grid_z_in,real_cube,  cube_edge_unit, off_in)

!  PURPOSE
!    Partitions the integration grid into manageable batches used
!    in integration sumation of Hartree potential for ESP charege fitting.
!
!  USES

  use dimensions, only: n_atoms, n_periodic, n_centers_integrals, &
                        n_centers_basis_I, n_species
  use esp_grids, only: grid_partitioned_esp,n_radial_esp,n_angular_esp,r_angular_esp,&
                       r_radial_esp,n_grid_batches_esp,batches_esp,n_my_batches_esp,&
                       n_max_batch_size_esp,n_full_points_esp
  use geometry, only: species, lattice_vector
  use species_data, only:  species_name
  use runtime_choices, only: use_hashed_maxmin, use_local_index, &
      batch_size_hard_limit
  use pbc_lists, only: coords_center, centers_basis_integrals, &
      inv_centers_basis_integrals
  use synchronize_mpi, only:sync_grid_partition
  use synchronize_mpi_basic, only: sync_find_max, sync_integer_vector, &
                                   sync_one_integer, sync_vector, &
                                   sync_real_number, get_max_double, &
                                   get_min_double
  use localorb_io, only: localorb_info, localorb_allinfo, &
                         OL_high, OL_low, OL_norm, output_priority, use_unit
  use constants, only: bohr
  use mpi_tasks, only: check_allocation, myid, n_tasks, use_mpi

  implicit none

  real*8,  dimension(n_species), INTENT(IN) :: radius_esp_min_new
  real*8,  dimension(n_species), INTENT(IN) :: radius_esp_max_new
  logical :: grid_output
  integer :: rm 
  logical :: equal_grid
  logical :: cube_rad_grid
  integer :: i_esp_equal_grid_x_in, i_esp_equal_grid_y_in, i_esp_equal_grid_z_in
  logical :: parallel_grid_partitioning = .true.
  logical,INTENT(IN)  :: real_cube
  real*8, dimension(3,3),INTENT(IN) ::        cube_edge_unit
  real*8, dimension(3), INTENT(IN) ::        off_in
!  ARGUMENTS
!    none
!  INPUTS
!   o radius_esp_min -- esp_min*vdw_radius min radius around atom where points
!                       will be created
!   o radius_esp_max -- esp_max*vdw_radius max radius around atom where points
!                       will be created
!   o grid_output -- output grids to file
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
  integer :: i_point, i_batch
  integer :: n_points_in_grid
  integer :: i_my_point, n_my_points
  
  integer, dimension(:), allocatable :: batch_sizes
  real*8, dimension(:,:), allocatable :: batch_coords

  real*8, dimension(:,:), allocatable :: all_coords
  real*8, dimension(:), allocatable :: point_weight
  real*8, dimension(3) :: coord_current
  real*8, dimension(3) :: coord_current_temp

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
  integer, dimension(:), allocatable :: batch_point
  integer, dimension(:), allocatable :: batch_index

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

  real*8, dimension(3) :: coord_max, coord_min, bb_dim, h
  integer, dimension(3) :: n_bins, idx

  integer :: hash_value, hash_method

  logical :: debug = .false.
  integer :: max_n_my_points, min_n_my_points

  integer, dimension(:), allocatable :: batch_task_list_esp

  real*8  :: dir_r(3)
  real*8  :: r_lattvec(3)
  real*8  :: lattice_vector_new(3,3)
  real*8  :: off(3)
  integer :: r_nvec(3)
  integer :: rn_x,rn_y,rn_z
  integer :: i_x,i_y,i_z,n_x,n_y,n_z
  real*8 :: dist
  integer ::  counter, counter_2, j_atom, i
  real*8, dimension(:,:), allocatable :: grid_out
  real*8, dimension(3) :: coord_new
  real*8, dimension(3) :: d_r
  real*8, dimension(3) :: cube_coord_min
  real*8, dimension(3) :: cube_coord_max
  integer :: cells_atoms
  integer :: n_points_in_batch_esp
  integer :: k, j  
  !----------------------------------------------------------------------------
  n_x = i_esp_equal_grid_x_in
  n_y = i_esp_equal_grid_y_in
  n_z = i_esp_equal_grid_z_in
  cube_coord_min = 0d0
  cube_coord_max = 0d0
  n_points_in_batch_esp = 100
  off(1) = sum(lattice_vector(1,:))*0.5
  off(2) = sum(lattice_vector(2,:))*0.5
  off(3) = sum(lattice_vector(3,:))*0.5
  ! this variable now lives in module grids.f
  !OTH & AJ: The next line seems to have been causing problems when plotting multiple
  !hartree_potential files. My understanding is that all that this line does is
  !preventing a re-calculation of the grid. The only reason for that possible is
  !performance. Hence, removing that line might make the whole thing slower, but
  !if that means the results are correct, I think that's worth it. 

  !if (grid_partitioned_esp)  return
  ! then, let's count the total number of points in the grid
  if (n_periodic > 0) then
    cells_atoms = ((rm*2+1)**3)*n_atoms-1
  else
    cells_atoms = (n_atoms-1)
  endif

  do k = 1, 3, 1
    do j = 1, 3, 1
       lattice_vector_new(k,j) = lattice_vector(k,j)
    enddo
  enddo

  if ((equal_grid.and.n_periodic.eq.0))then
    lattice_vector_new(:,:) = 0d0
    cube_coord_min = 0d0
    cube_coord_max = 0d0
    do i_atom = 1, n_atoms, 1
       if(i_atom.eq.1)then
         cube_coord_min(1) = coords_center( 1,i_atom )-&
                                  radius_esp_max_new(species(i_atom))
         cube_coord_min(2) = coords_center( 2,i_atom )-&
                                  radius_esp_max_new(species(i_atom))
         cube_coord_min(3) = coords_center( 3,i_atom )-&
                                  radius_esp_max_new(species(i_atom))
         cube_coord_max(1) = coords_center( 1,i_atom )+&
                                  radius_esp_max_new(species(i_atom))
         cube_coord_max(2) = coords_center( 2,i_atom )+&
                                  radius_esp_max_new(species(i_atom))
         cube_coord_max(3) = coords_center( 3,i_atom )+&
                                  radius_esp_max_new(species(i_atom))
       else
         cube_coord_min(1) = min(cube_coord_min(1),coords_center( 1,i_atom )-&
                                  radius_esp_max_new(species(i_atom)))
         cube_coord_min(2) = min(cube_coord_min(2),coords_center( 2,i_atom )-&
                                  radius_esp_max_new(species(i_atom)))
         cube_coord_min(3) = min(cube_coord_min(3),coords_center( 3,i_atom )-&
                                  radius_esp_max_new(species(i_atom)))
         cube_coord_max(1) = max(cube_coord_max(1),coords_center( 1,i_atom )+&
                                  radius_esp_max_new(species(i_atom)))
         cube_coord_max(2) = max(cube_coord_max(2),coords_center( 2,i_atom )+&
                                  radius_esp_max_new(species(i_atom)))
         cube_coord_max(3) = max(cube_coord_max(3),coords_center( 3,i_atom )+&
                                  radius_esp_max_new(species(i_atom)))
       endif
    enddo
    lattice_vector_new(1,1) = (cube_coord_max(1)-&
                               cube_coord_min(1))
    lattice_vector_new(2,2) = (cube_coord_max(2)-&
                               cube_coord_min(2))
    lattice_vector_new(3,3) = (cube_coord_max(3)-&
                               cube_coord_min(3))
    off(1) = cube_coord_min(1)
    off(2) = cube_coord_min(2)
    off(3) = cube_coord_min(3)
  endif
  
  if ((equal_grid.and.n_periodic.ne.0.and.cube_rad_grid)) then
    lattice_vector_new(3,3) = maxval(coords_center( 3,1:n_atoms )) &
                            + maxval(radius_esp_max_new)
    off(3)                  = minval(coords_center( 3,1:n_atoms )) &
                            - maxval(radius_esp_max_new)
    lattice_vector_new(3,3) = abs(lattice_vector_new(3,3) - off(3))

    if (abs(lattice_vector_new(3,3)) > abs(lattice_vector(3,3))) then
      off(3)                  = -0.5 * lattice_vector(3,3)
      lattice_vector_new(3,3) = lattice_vector(3,3)
    endif
  endif
  if (real_cube) then
    do k = 1, 3, 1
        lattice_vector_new(k,1) = cube_edge_unit(k,1)*n_x
        lattice_vector_new(k,2) = cube_edge_unit(k,2)*n_y
        lattice_vector_new(k,3) = cube_edge_unit(k,3)*n_z
        off(k) = -off_in(k)
    enddo
  endif
  n_points_in_grid = 0
  counter = 0
    ! first count number of points, generated by superposition of radial grids
  if (equal_grid)then
    if(sum(radius_esp_min_new).gt.0)then
      do i_x = 0, n_x-1, 1 
      do i_y = 0, n_y-1, 1
      do i_z = 0, n_z-1, 1
	d_r(1) = dble(i_x)/n_x
	d_r(2) = dble(i_y)/n_y
	d_r(3) = dble(i_z)/n_z
	coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
	call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
	if (counter_2.eq.(cells_atoms+1))then
	      counter=counter+1
	endif
      enddo
      enddo
      enddo
    else
      if (n_periodic > 0) then
        call map_to_center_cell(coord_new(1:3))
      endif
      counter=n_x*n_y*n_z
    endif
  else
    do i_atom = 1, n_atoms, 1
      do i_radial = 1, n_radial_esp(species(i_atom)), 1
	do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
	  coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	  r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	  r_radial_esp(i_radial, species(i_atom))
          call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	  if (counter_2.eq.cells_atoms)then
	      counter=counter+1
	  endif
	enddo
      enddo
    enddo
  endif
  n_points_in_grid=counter
  if(grid_output)then
      counter = 0
      if (equal_grid)then
        if(sum(radius_esp_min_new).gt.0)then
	do i_x = 0, n_x-1, 1 
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
          d_r(1) = dble(i_x)/n_x
          d_r(2) = dble(i_y)/n_y
          d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
	  call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
	  if (counter_2.eq.(cells_atoms+1))then
		counter=counter+1
	  endif
	enddo
	enddo
	enddo
        else
          if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
          endif
          counter=n_x*n_y*n_z
        endif
      else
        do i_atom = 1, n_atoms, 1
	  do i_radial = 1, n_radial_esp(species(i_atom)), 1
	    do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
              coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	      r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	      r_radial_esp(i_radial, species(i_atom))
              call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	      if (counter_2.eq.cells_atoms)then
	        counter=counter+1
	      endif
	    enddo
	  enddo
        enddo
      endif
      allocate(grid_out(3, counter),stat=info)
      call check_allocation(info, 'grid_out               ')
      counter = 0
      if (equal_grid)then
	do i_x = 0, n_x-1, 1 
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
          d_r(1) = dble(i_x)/n_x
          d_r(2) = dble(i_y)/n_y
          d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
            call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
          else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2 = cells_atoms+1
          endif
	  if (counter_2.eq.(cells_atoms+1))then
		counter=counter+1
                grid_out(:,counter) = coord_new(:)
	  endif
	enddo
	enddo
	enddo
      else
        do i_atom = 1, n_atoms, 1
	  do i_radial = 1, n_radial_esp(species(i_atom)), 1
	    do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
              coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	      r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	      r_radial_esp(i_radial, species(i_atom))
              call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	      if (counter_2.eq.cells_atoms)then
	        counter=counter+1
	        grid_out(:,counter) = coord_new(:)
	      endif
	    enddo
	  enddo
        enddo
      endif
      open (UNIT=50, FILE='grid_out2.dat', STATUS='REPLACE')
      write(50,"(A12, 3A18)")'# x', 'y', 'z' 
      do i=1, counter
	write(50,"(3f18.6)")  grid_out(1,i)*bohr, grid_out(2,i)*bohr, &
                              grid_out(3,i)*bohr
      end do
      close(50)
      deallocate(grid_out)
  endif
  if (parallel_grid_partitioning) then

     call find_bb_dimensions(coord_max, coord_min, bb_dim)
     call find_bins(bb_dim, n_grid_batches_esp, n_bins, h)
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
      if (equal_grid)then
	do i_x = 0, n_x-1, 1 
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
         d_r(1) = dble(i_x)/n_x
         d_r(2) = dble(i_y)/n_y
         d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
	    call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
          else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2=(cells_atoms+1)
          endif
	  if (counter_2.eq.(cells_atoms+1))then
                 coord_current(:) = coord_new(:)        
                 coord_current_temp = coord_current 
                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current_temp(1:3) )
                 end if

                 idx = CEILING((coord_current_temp - coord_min)/h)
                 idx = MIN(idx, n_bins)
                 idx = MAX(idx, 1)
        
                 hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
                              n_bins(2)*n_bins(1)
        
                 if (myid+1 == hash_value) n_my_points = n_my_points + 1
	  endif
	enddo
	enddo
	enddo
      else
	do i_atom = 1, n_atoms, 1
	  do i_radial = 1, n_radial_esp(species(i_atom)), 1
	    do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
	      coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	      r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	      r_radial_esp(i_radial, species(i_atom))
              call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	      if (counter_2.eq.cells_atoms)then
                 coord_current(:) = coord_new(:)        
                 coord_current_temp = coord_current 
                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current_temp(1:3) )
                 end if

                 idx = CEILING((coord_current_temp - coord_min)/h)
                 idx = MIN(idx, n_bins)
                 idx = MAX(idx, 1)
        
                 hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
                              n_bins(2)*n_bins(1)
        
                 if (myid+1 == hash_value) n_my_points = n_my_points + 1
	      endif
	    enddo
	  enddo
	enddo
       endif
     end select
     ! allocate grid storage for the current thread only
     allocate(all_coords(3,n_my_points),stat=info)
     call check_allocation(info, 'all_coords                    ')

     root%size = n_my_points
     allocate( root%point_indexes(root%size),stat=info)
     call check_allocation(info, ' root%point_indexes           ')

     allocate( root%global_point_indexes(root%size),stat=info )
     call check_allocation(info, 'root%global_point_indexes     ')

     ! allocate array for the weights for the points in the current thread
     allocate(point_weight(n_my_points),stat=info)
     call check_allocation(info, 'point_weight                  ')

     n_max_compute_part = 0
     i_point = 0
     i_my_point = 0
     if (equal_grid)then
	do i_x = 0, n_x-1, 1 
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
         d_r(1) = dble(i_x)/n_x
         d_r(2) = dble(i_y)/n_y
         d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
	    call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
          else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2=(cells_atoms+1)
          endif
	 if (counter_2.eq.(cells_atoms+1))then
		  i_point = i_point + 1
		  coord_current(:) = coord_new(:)
		  
		  coord_current_temp = coord_current 

		  if(n_periodic > 0)then
		    call map_to_center_cell(coord_current_temp(1:3) )
		  end if

		  ! for the parallel routines it is best if the points are 
		  ! distributed as uniformly as possible over the threads, hence
		  !  this distribution and not the one based on radial shells
		  select case(hash_method)
		  case(1)
		    hash_value = MOD(i_point, n_tasks) + 1
		  case(2)
		    idx = CEILING((coord_current_temp - coord_min)/h)
		    idx = MIN(idx, n_bins)
		    idx = MAX(idx, 1)
		    
		    hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
                                 n_bins(2)*n_bins(1)
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
		  
		    ! compute atom-centered coordinates of current integration 
		    ! point, as viewed from all atoms
		    call tab_atom_centered_coords_p0 &
			  ( coord_current,  dist_tab_sq, dir_tab, &
			  n_centers_integrals, centers_basis_integrals )
		  
		    ! determine which basis functions are relevant at current 
		    ! integration point, and tabulate their indices
		    call prune_basis_p0X ( dist_tab_sq, &
			  n_compute_a, n_compute_c, i_basis,  &
			  n_centers_basis_I, n_centers_integrals, &
			  inv_centers_basis_integrals, i_my_point )
    !NEC_CB Using version which saves massively indirect addressed arrays

		    
		    ! store the relaitve weight of current integration point
		    n_max_compute_part = MAX(n_max_compute_part, n_compute_a)
		    point_weight(i_my_point) = dble(n_compute_a)
			
		  end if
	  endif
	enddo
	enddo
	enddo
     else
      do i_atom = 1, n_atoms, 1
        do i_radial = 1, n_radial_esp(species(i_atom)), 1
           do i_angular = 1, n_angular_esp( i_radial,species(i_atom) )
	      coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	      r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	      r_radial_esp(i_radial, species(i_atom))
              call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	      if (counter_2.eq.cells_atoms)then
		  i_point = i_point + 1
		  coord_current(:) = coord_new(:)
		  
		  coord_current_temp = coord_current 

		  if(n_periodic > 0)then
		    call map_to_center_cell(coord_current_temp(1:3) )
		  end if

		  ! for the parallel routines it is best if the points are 
		  ! distributed as uniformly as possible over the threads, hence
		  !  this distribution and not the one based on radial shells
		  select case(hash_method)
		  case(1)
		    hash_value = MOD(i_point, n_tasks) + 1
		  case(2)
		    idx = CEILING((coord_current_temp - coord_min)/h)
		    idx = MIN(idx, n_bins)
		    idx = MAX(idx, 1)
		    
		    hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
                                 n_bins(2)*n_bins(1)
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
		  
		    ! compute atom-centered coordinates of current integration 
		    ! point, as viewed from all atoms
		    call tab_atom_centered_coords_p0 &
			  ( coord_current,  dist_tab_sq, dir_tab, &
			  n_centers_integrals, centers_basis_integrals )
		  
		    ! determine which basis functions are relevant at current 
		    ! integration point, and tabulate their indices
		    call prune_basis_p0X ( dist_tab_sq, &
			  n_compute_a, n_compute_c, i_basis,  &
			  n_centers_basis_I, n_centers_integrals, &
			  inv_centers_basis_integrals, i_my_point )
    !NEC_CB Using version which saves massively indirect addressed arrays

		    
		    ! store the relaitve weight of current integration point
		    n_max_compute_part = MAX(n_max_compute_part, n_compute_a)
		    point_weight(i_my_point) = dble(n_compute_a)
			
		  end if
  	      endif              
           end do
        end do
      end do
     endif
     if (debug) then
        call sync_find_max(n_my_points, max_n_my_points)
        n_my_points = -n_my_points
        call sync_find_max(n_my_points, min_n_my_points)
        n_my_points = -n_my_points
        min_n_my_points = -min_n_my_points
        write(use_unit,'(A,I5,A,I7,A,I7,A,I7)') 'Myid', myid, ' n_my_points', &
             n_my_points, &
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
     allocate(grid_partition(n_my_points),stat=info)
     call check_allocation(info, 'grid_partition                ')

     grid_partition = 0


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
     if (allocated(root%point_indexes)) then
        deallocate(root%point_indexes)
     end if
     if (allocated(root%global_point_indexes)) then
        deallocate(root%global_point_indexes)
     end if
     if (allocated(point_weight)) then
        deallocate(point_weight)
     end if


     ! allocate the batch data size array
     allocate(batch_sizes(n_grid_batches_esp),stat=info)
     call check_allocation(info, 'batch_sizes                   ')

     allocate(batch_coords(3,n_grid_batches_esp),stat=info)
     call check_allocation(info, 'batch_coords                  ')


     ! compute the sizes of the batches according to the distribution scheme
     batch_sizes = 0
     batch_coords = 0d0
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
     call sync_integer_vector(batch_sizes, n_grid_batches_esp)
     call sync_vector(batch_coords, 3*n_grid_batches_esp)
     batch_coords(1,:) = batch_coords(1,:) / MAX(batch_sizes(:),1)
     batch_coords(2,:) = batch_coords(2,:) / MAX(batch_sizes(:),1)
     batch_coords(3,:) = batch_coords(3,:) / MAX(batch_sizes(:),1)
     ! print out some information on the batches
     call print_out_batch_report( )
     
     ! create the batch tasking list
     if(use_local_index) then
        call distribute_batch_tasks_by_location_esp( batch_sizes, batch_coords)
     else
        call distribute_batch_tasks_esp( batch_sizes )
     endif
     ! create the batches
     call create_grid_batches( )
     grid_partitioned_esp = .true.

     if (allocated(all_coords)) then
        deallocate(all_coords)
     end if
     if (allocated(batch_sizes)) then
        deallocate(batch_sizes)
     end if
     if (allocated(batch_coords)) then
        deallocate(batch_coords)
     end if
     if (allocated(grid_partition)) then
        deallocate(grid_partition)
     end if

     return
  end if
  

      if (allocated(batch_task_list_esp)) then
	deallocate(batch_task_list_esp)
      end if

  contains
!******
!----------------------------------------------------------------------------------
!****s* esp_partition_grid/check_point
!  NAME
!    check_point
!  SYNOPSIS

    subroutine check_point(rm,lattice_vector,radius_esp_min_new,coord_new,counter_2)


!  PURPOSE
!    Check if point will be included
!    o n_grid_batches_esp - total number of batches in the grid
!    o n_my_batches_esp - the number of batches the thread has
!    o batch_sizes - size of each batch
!    o batch_task_list_esp - distribution of the batches over the threads
!    o grid_partition - partition of the points of the grid into batches
!    In the end the array batches is created that stores all the information 
!    threadwise
!  USES
  use dimensions, only: n_atoms, n_periodic, n_species
  use geometry, only: species
  use pbc_lists, only: coords_center

      implicit none
!  ARGUMENTS
      integer,INTENT(IN)	:: rm
      real*8,INTENT(IN)		:: lattice_vector(3,3)
      real*8,INTENT(IN)		:: radius_esp_min_new(n_species)
      real*8,INTENT(IN)		:: coord_new(3)
      integer,INTENT(OUT)	:: counter_2

!  INPUTS
!    
!  OUTPUT
!    
!  SOURCE


      integer :: rn_x,rn_y,rn_z, j_atom
      real*8 :: r_nvec(3),r_lattvec(3)

      counter_2 = 0
      if (n_periodic > 0) then
         call map_to_center_cell(coord_new(1:3))
      endif
      do j_atom = 1, n_atoms, 1
         if (n_periodic > 0) then
            do rn_x=-rm,rm,1
            do rn_y=-rm,rm,1
            do rn_z=-rm,rm,1
	        r_nvec(1) = rn_x
	        r_nvec(2) = rn_y
	        r_nvec(3) = rn_z
	        r_lattvec(1:3) = matmul(lattice_vector,&
					r_nvec(1:3))
	        dir_r(:)=coord_new(:)-(coords_center( :,j_atom )+&
				r_lattvec(:))
		dist=sqrt(sum(dir_r**2))
		if (dist>=radius_esp_min_new(species(j_atom)))then
			  counter_2=counter_2+1
		endif
            enddo
            enddo
            enddo
         else
	    dist =  sqrt(sum((coord_new(:)-coords_center( :,j_atom ))**2))
	     if (dist>=radius_esp_min_new(species(j_atom)))then
		  counter_2=counter_2+1
	     endif
         endif
      enddo
    end subroutine check_point
!******
!----------------------------------------------------------------------------------
!****s* esp_partition_grid/check_point
!  NAME
!    check_point
!  SYNOPSIS

    subroutine check_point_rad(rm,lattice_vector,radius_esp_min_new,coord_new,counter_2)


!  PURPOSE
!    Check if point will be included
!    o n_grid_batches_esp - total number of batches in the grid
!    o n_my_batches_esp - the number of batches the thread has
!    o batch_sizes - size of each batch
!    o batch_task_list_esp - distribution of the batches over the threads
!    o grid_partition - partition of the points of the grid into batches
!    In the end the array batches is created that stores all the information 
!    threadwise
!  USES
  use dimensions, only: n_atoms, n_periodic, n_species
  use geometry, only: species
  use pbc_lists, only: coords_center

      implicit none
!  ARGUMENTS
      integer,INTENT(IN)	:: rm
      real*8,INTENT(IN)		:: lattice_vector(3,3)
      real*8,INTENT(IN)		:: radius_esp_min_new(n_species)
      real*8,INTENT(IN)		:: coord_new(3)
      integer,INTENT(OUT)	:: counter_2

!  INPUTS
!    
!  OUTPUT
!    
!  SOURCE


      integer :: rn_x,rn_y,rn_z, j_atom
      real*8 :: r_nvec(3),r_lattvec(3)

	  counter_2=0
	  do j_atom = 1, n_atoms, 1
		if (n_periodic > 0) then
		      do rn_x=-rm,rm,1
		      do rn_y=-rm,rm,1
		      do rn_z=-rm,rm,1
			r_nvec(1) = rn_x
			r_nvec(2) = rn_y
			r_nvec(3) = rn_z
			if (i_atom.ne.j_atom.or.(sum(abs(r_nvec)).ne.0)) then
			  r_lattvec(1:3) = matmul(lattice_vector,&
					  r_nvec(1:3))
			  dir_r(1:3)=coord_new(1:3)-(coords_center( 1:3,j_atom )+&
				  r_lattvec(1:3))
			  dist=sqrt(sum(dir_r**2))
			  if (dist>=radius_esp_min_new(species(j_atom)))then
			    counter_2=counter_2+1
			  endif
			endif
		      enddo
		      enddo
		      enddo
                     if (n_periodic > 0) then
		      call map_to_center_cell(coord_new(1:3))
                     endif
		else
		    if (i_atom.ne.j_atom) then
		      dist =  sqrt(sum((coord_new(1:3) &
			    -coords_center( 1:3,j_atom ))**2))
		      if (dist>=radius_esp_min_new(species(j_atom)))then
			counter_2=counter_2+1
		      endif
		    endif    
		endif
	  enddo
    end subroutine check_point_rad
!******
!----------------------------------------------------------------------------------
!****s* esp_partition_grid/create_grid_batches
!  NAME
!    create_grid_batches
!  SYNOPSIS

    subroutine create_grid_batches()

!  PURPOSE
!    Creates the batches for each thread. At this point the following variables 
!    are set:
!    o n_grid_batches_esp - total number of batches in the grid
!    o n_my_batches_esp - the number of batches the thread has
!    o batch_sizes - size of each batch
!    o batch_task_list_esp - distribution of the batches over the threads
!    o grid_partition - partition of the points of the grid into batches
!    In the end the array batches is created that stores all the information 
!    threadwise
!  USES

      use mpi_tasks, only: aims_stop
      implicit none
!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  SOURCE


      integer :: i_my_batch, i_batch, i_point
      integer :: i_atom, i_radial, i_angular, info
      integer :: current_batch, my_current_batch
      integer :: n_atom_points, max_atom_points, i_my_point, i_offset
      integer, allocatable :: global_grid_partition(:)

      allocate(batches_esp(n_my_batches_esp),stat=info)
      call check_allocation(info, 'batches_esp                       ')

      allocate(batch_index(n_grid_batches_esp),stat=info)
      call check_allocation(info, 'batch_index                   ')

      batch_index = 0

      ! distribute the the batches to threads according to the task list
      ! and allocate the storage space for each thread
      i_my_batch = 0
      do i_batch = 1, n_grid_batches_esp, 1
         if (myid.eq.batch_task_list_esp(i_batch)) then
            i_my_batch = i_my_batch + 1
            batches_esp(i_my_batch) % size_esp = batch_sizes(i_batch)
            batch_index(i_batch) = i_my_batch
            if (batch_sizes(i_batch).gt.0) then
               allocate( batches_esp(i_my_batch) % points_esp(&
                         batch_sizes(i_batch)),stat=info )
               call check_allocation(info, 'batches_esp                       ')

            end if
         end if
      end do

      allocate(batch_point(n_my_batches_esp),stat=info)
      call check_allocation(info, 'batch_point                   ')

      batch_point = 0

      ! Get the maximum number of points per atom

      max_atom_points = 0
      if (equal_grid)then
	do i_x = 0, n_x-1, 1 
        n_atom_points = 0
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
         d_r(1) = dble(i_x)/n_x
         d_r(2) = dble(i_y)/n_y
         d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
	    call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
           else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2=(cells_atoms+1)
           endif
	   if (counter_2.eq.(cells_atoms+1))then
                n_atom_points = n_atom_points + 1
	   endif
	enddo
	enddo
   !     max_atom_points = ceiling(dble(n_atom_points)/n_atoms)
        max_atom_points = max(max_atom_points,n_atom_points)
	enddo
      else
      do i_atom = 1, n_atoms, 1
        n_atom_points = 0
	do i_radial = 1, n_radial_esp(species(i_atom)), 1
	  do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
	    coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	    r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	    r_radial_esp(i_radial, species(i_atom))
            call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	    if (counter_2.eq.cells_atoms)then
               n_atom_points = n_atom_points + 1
	    endif
	  enddo
	enddo
        max_atom_points = max(max_atom_points,n_atom_points)
      enddo
      endif
      allocate(global_grid_partition(max_atom_points),stat=info)
      call check_allocation(info, 'global_grid_partition         ')

      ! store all required grid information to the thread that is treating the 
      ! batch
      i_point = 0
      i_my_point = 0
      i_offset = 0
      if (equal_grid)then
	      do i_x = 0, n_x-1, 1 
            i_offset = i_point
            if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then
              ! do i_x = 0, n_x -1, 1 
              ! i_offset = i_point   
               global_grid_partition = 0
	            do i_y = 0, n_y-1, 1
	            do i_z = 0, n_z-1, 1
	               d_r(1) = dble(i_x)/n_x
	               d_r(2) = dble(i_y)/n_y
	               d_r(3) = dble(i_z)/n_z
	               coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
	               if(sum(radius_esp_min_new).gt.0)then
		               call check_point(rm,lattice_vector,radius_esp_min_new,&
			               coord_new,counter_2)
	               else
		               if (n_periodic > 0) then
		                  call map_to_center_cell(coord_new(1:3))
		               endif
		               counter_2=(cells_atoms+1)
	               endif
	               if (counter_2.eq.(cells_atoms+1))then
			            i_point = i_point + 1
			            if (myid == MOD(i_point, n_tasks)) then
			               i_my_point = i_my_point + 1
			               global_grid_partition(i_point-i_offset) = &
			               grid_partition(i_my_point)
			            endif
	               endif
	            enddo
	            enddo
               !enddo
               call sync_integer_vector(global_grid_partition,i_point-i_offset)
            elseif (parallel_grid_partitioning.and.use_hashed_maxmin) then
               !do i_x = 0, n_x -1, 1 
               !i_offset = i_point    
               global_grid_partition = 0
	            do i_y = 0, n_y-1, 1
	            do i_z = 0, n_z-1, 1
	               d_r(1) = dble(i_x)/n_x
	               d_r(2) = dble(i_y)/n_y
	               d_r(3) = dble(i_z)/n_z
	               coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
	               if(sum(radius_esp_min_new).gt.0)then
		               call check_point(rm,lattice_vector,radius_esp_min_new,&
			                  coord_new,counter_2)
	               else
		               if (n_periodic > 0) then
		                  call map_to_center_cell(coord_new(1:3))
		               endif
		               counter_2=(cells_atoms+1)
	               endif
	               if (counter_2.eq.(cells_atoms+1))then
			            i_point = i_point + 1
			            coord_current(:) = coord_new(:)

			            coord_current_temp = coord_current 
                     if (n_periodic > 0) then 
                        call map_to_center_cell(coord_current_temp(1:3))
                     endif

			            idx = CEILING((coord_current_temp - coord_min)/h)
			            idx = MIN(idx, n_bins)
			            idx = MAX(idx, 1)
			            hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
			            n_bins(2)*n_bins(1)

			            if (myid+1 == hash_value) then
			               i_my_point = i_my_point + 1
			               global_grid_partition(i_point-i_offset) = &
			                  grid_partition(i_my_point)
			            endif
	               endif
	            enddo
	            enddo
              ! enddo
               !if (debug) then
               !   call sync_find_max(n_my_points, max_n_my_points)
               !   n_my_points = -n_my_points
               !   call sync_find_max(n_my_points, min_n_my_points)
               !   n_my_points = -n_my_points
               !    min_n_my_points = -min_n_my_points
               !    write(use_unit,'(A,I5,A,I7)') 'Myid', myid, ' i_my_points', &
               !       i_my_point
               !end if


               call sync_integer_vector(global_grid_partition,i_point-i_offset)
            else
             !  do i_x = 0, n_x-1, 1 
             !  i_offset = i_point    
	            do i_y = 0, n_y-1, 1
	            do i_z = 0, n_z-1, 1
                  d_r(1) = dble(i_x)/n_x
                  d_r(2) = dble(i_y)/n_y
                  d_r(3) = dble(i_z)/n_z
	               coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
                  if(sum(radius_esp_min_new).gt.0)then
	                  call check_point(rm,lattice_vector,radius_esp_min_new,&
			               coord_new,counter_2)
                  else
                     if (n_periodic > 0) then
                        call map_to_center_cell(coord_new(1:3))
                     endif
                     counter_2=(cells_atoms+1)
                  endif
	               if (counter_2.eq.(cells_atoms+1))then
                    i_point = i_point + 1
                    global_grid_partition(i_point-i_offset) = &
                    grid_partition(i_point)
	               endif
	            enddo
	            enddo
            !   enddo 
            endif
           ! do i_x = 0, n_x-1, 1 
           ! i_offset = i_point    
            i_point = i_offset ! reset it for following loop
	do i_y = 0, n_y-1, 1
	do i_z = 0, n_z-1, 1
         d_r(1) = dble(i_x)/n_x
         d_r(2) = dble(i_y)/n_y
         d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
	    call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
          else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2=(cells_atoms+1)
          endif
	  if (counter_2.eq.(cells_atoms+1))then
		 i_point = i_point + 1
                     current_batch = global_grid_partition(i_point-i_offset)
                  
!test
!
! VB: This test statement kept for now, on August 11, 2012. 
!
!     There is a specific machine on which I can produce a segfault here during
!     unit cell optimization (fcc Al test case, 4 mpi tasks), but I am not able 
!     to reproduce the problem on a formally identical machine with a new  
!     compile one day later. Nothing related to this place has changed. My
!     suspicion is that I am either facing a particularly nasty case of an 
!     uninitialized variable which only happensfor the unit cell optimization 
!     case - or that the first test machine produced a different compilation 
!     (cpu bug?). What I do see is that there are tiny numerical differences
!     between both compiled versions, same settings, on the same machine.
!
!     As I am unable to reproduce the problem in a consistent manner, 
!     I am placing this documentation here for now.
!
               if (current_batch.eq.0) then
                  write(use_unit,*) '*** Error on myid ', myid, ': Point ', i_point, '.'
                  call aims_stop("Batch not found","create_grid_batches in esp_partition_grid")
               end if
!test end
               my_current_batch = batch_index(current_batch)

               coord_current(:) = coord_new(:)
               if (my_current_batch /= 0) then
                  batch_point(my_current_batch)=batch_point(my_current_batch)+ 1
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%coords_esp(:) = &
                       coord_current(:)
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_atom_esp = &
                       i_x
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_radial_esp = &
                       i_y
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_angular_esp = &
                       i_z
               end if
	  endif
	enddo
	enddo
	enddo


      else

      do i_atom = 1, n_atoms, 1
         i_offset = i_point

         if (parallel_grid_partitioning.and.(.not.use_hashed_maxmin)) then

            global_grid_partition = 0
            do i_radial = 1, n_radial_esp(species(i_atom)), 1
               do i_angular = 1, n_angular_esp( i_radial,species(i_atom) )
		  coord_new(1:3)=coords_center( 1:3,i_atom )+ &
		  r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
		  r_radial_esp(i_radial, species(i_atom))
                 call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
		  if (counter_2.eq.cells_atoms)then
		    i_point = i_point + 1
		    if (myid == MOD(i_point, n_tasks)) then
		      i_my_point = i_my_point + 1
		      global_grid_partition(i_point-i_offset) = &
                      grid_partition(i_my_point)
		    endif
                  endif
               end do
            end do
            call sync_integer_vector(global_grid_partition,i_point-i_offset)

         elseif (parallel_grid_partitioning.and.use_hashed_maxmin) then
            global_grid_partition = 0

            do i_radial = 1, n_radial_esp(species(i_atom)), 1
               do i_angular = 1, n_angular_esp( i_radial,species(i_atom) )
		  coord_new(1:3)=coords_center( 1:3,i_atom )+ &
		  r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
		  r_radial_esp(i_radial, species(i_atom))
                 call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
		  if (counter_2.eq.cells_atoms)then
		      i_point = i_point + 1

		      coord_current(:) = coord_new(:)

		    coord_current_temp = coord_current 

		    if(n_periodic > 0)then
			call map_to_center_cell(coord_current_temp(1:3) )
		    end if

		      idx = CEILING((coord_current_temp - coord_min)/h)
		      idx = MIN(idx, n_bins)
		      idx = MAX(idx, 1)
		      hash_value = idx(1) + (idx(2)-1)*n_bins(1) + (idx(3)-1)*&
                      n_bins(2)*n_bins(1)

		      if (hash_value == myid+1) then
			i_my_point = i_my_point + 1
			global_grid_partition(i_point-i_offset) = &
                        grid_partition(i_my_point)
		      endif
                   endif
               end do
            end do
            call sync_integer_vector(global_grid_partition,i_point-i_offset)
         else
            do i_radial = 1, n_radial_esp(species(i_atom)), 1
               do i_angular = 1, n_angular_esp( i_radial,species(i_atom) )
		  coord_new(1:3)=coords_center( 1:3,i_atom )+ &
		  r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
		  r_radial_esp(i_radial, species(i_atom))
                 call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
		  if (counter_2.eq.cells_atoms)then
                    i_point = i_point + 1
                    global_grid_partition(i_point-i_offset) = &
                    grid_partition(i_point)
                  endif
               end do
            end do

         endif
         i_point = i_offset ! reset it for following loop
         do i_radial = 1, n_radial_esp(species(i_atom)), 1
            do i_angular = 1, n_angular_esp( i_radial,species(i_atom) )
		  coord_new(1:3)=coords_center( 1:3,i_atom )+ &
		  r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
		  r_radial_esp(i_radial, species(i_atom))
                 call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
		  if (counter_2.eq.cells_atoms)then
                     i_point = i_point + 1
                     current_batch = global_grid_partition(i_point-i_offset)
                  
!test
!
! VB: This test statement kept for now, on August 11, 2012. 
!
!     There is a specific machine on which I can produce a segfault here during
!     unit cell optimization (fcc Al test case, 4 mpi tasks), but I am not able 
!     to reproduce the problem on a formally identical machine with a new  
!     compile one day later. Nothing related to this place has changed. My
!     suspicion is that I am either facing a particularly nasty case of an 
!     uninitialized variable which only happensfor the unit cell optimization 
!     case - or that the first test machine produced a different compilation 
!     (cpu bug?). What I do see is that there are tiny numerical differences
!     between both compiled versions, same settings, on the same machine.
!
!     As I am unable to reproduce the problem in a consistent manner, 
!     I am placing this documentation here for now.
!
               if (current_batch.eq.0) then
                  write(use_unit,*) '*** Error on myid ', myid, ': Point ', i_point, '.'
                  stop
               end if
!test end
               my_current_batch = batch_index(current_batch)
               coord_current(:) = coord_new(:)
               if (my_current_batch /= 0) then
                  batch_point(my_current_batch)=batch_point(my_current_batch)+ 1
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%coords_esp(:) = &
                       coord_current(:)
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_atom_esp = &
                       i_atom
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_radial_esp = &
                       i_radial
                  batches_esp(my_current_batch)%points_esp(&
                       batch_point(my_current_batch))%index_angular_esp = &
                       i_angular
               end if
             endif
            end do
         end do
      end do
endif
      deallocate(batch_index)
      deallocate(batch_point)
      deallocate(global_grid_partition)
      
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
           "| Number of batches:  ", n_grid_batches_esp
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      n_max_batch_size_esp = MAXVAL(batch_sizes)
      write (info_str,'(2X,A,I7)') &
           "| Maximal batch size: ", n_max_batch_size_esp
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      write (info_str,'(2X,A,I7)') &
           "| Minimal batch size: ", MINVAL(batch_sizes)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      avg_batch_size = dble(n_points_in_grid)/dble(n_grid_batches_esp)
      write (info_str,'(2X,A,F11.3)') &
           "| Average batch size: ", avg_batch_size
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
      batch_size_variance = 0.0d0
      do i_batch = 1, n_grid_batches_esp, 1
         batch_size_variance = batch_size_variance + &
              (batch_sizes(i_batch) - avg_batch_size)**2
      end do
      batch_size_variance = batch_size_variance / dble(n_grid_batches_esp)
      write (info_str,'(2X,A,F11.3)') &
           "| Standard deviation of batch sizes: ", sqrt(batch_size_variance)
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      
    end subroutine print_out_batch_report
!******
!-------------------------------------------------------------------------------
!****s* partition_grid/build_maxmin_partition_tree
!  NAME
!    build_maxmin_partition_tree
!  SYNOPSIS
    recursive subroutine build_maxmin_partition_tree( node )
!  PURPOSE
!    Computes a partition of the grid points using a recursive adapted cut-plane
!     method. As a result the array grid_partition is filled with the batch 
!     indexes.
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
      if ((total_node_size*total_max_point_weight > n_points_in_batch_esp) .or.&
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
                       ddot( node%size, coords_vector(:,i_index), 1, &
                       coords_vector(:,i_index_2), 1 )
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

               call DGEMV('N', node%size, 3, 1.0d0, coords_vector, node%size, &
                    normal, 1, 0.0d0, vector, 1)
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

            call DGEMV('N', node%size, 3, 1.0d0, coords_vector, node%size, &
                 normal, 1, 0.0d0, vector, 1)

            temp_vector = vector
            call heapsort(temp_vector, node%size)
            median = temp_vector(INT(dble(node%size)/2.0))
         end if

         ! compute the sizes of the child nodes
         size_children = 0
         do i_index = 1, node%size, 1
            if (vector(i_index) <= median) then
               size_children(1) = size_children(1) + 1
            else
               size_children(2) = size_children(2) + 1
            end if
         end do
         ! allocate the child nodes
         do i_index = 1, node % n_children, 1
            node % children(i_index) % size = size_children(i_index)
            if (size_children(i_index) > 0) then
               allocate( node % children(i_index) % point_indexes(&
               size_children(i_index)),stat=info)
               call check_allocation(info, 'node % children               ')

               if (parallel_grid_partitioning) then
                  allocate( node % children(i_index) % global_point_indexes&
                  (size_children(i_index)),stat=info)
                  call check_allocation(info, 'node % children               ')

               end if
            end if
         end do

         ! create the child nodes by geometric split with the median of the 
         ! plane coordinates
         i_point = 0
         do i_index = 1, node%size, 1
            if (vector(i_index) <= median) then
               i_point(1) = i_point(1) + 1
               node % children(1) % point_indexes(i_point(1)) = node % &
               point_indexes(i_index)
               if (parallel_grid_partitioning) then
                  node % children(1) % global_point_indexes(i_point(1)) = &
                       node % global_point_indexes(i_index)
               end if
            else
               i_point(2) = i_point(2) + 1
               node % children(2) % point_indexes(i_point(2)) = node % &
               point_indexes(i_index)
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

         ! if the node is small enough (and > 0) then add it as a batch to the 
         ! system
         if (total_node_size > 0) then
            n_grid_batches_esp = n_grid_batches_esp + 1
            do i_index = 1, node%size, 1
               grid_partition(node%point_indexes(i_index)) = n_grid_batches_esp
            end do
         end if

      end if
    end subroutine build_maxmin_partition_tree
!****
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/hash_points
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

      integer :: i_offset

      if (use_hashed_maxmin) then

         n_grid_batches_esp = 0
         call build_maxmin_partition_tree( root )
!      print *, 'MYID: ', myid, 'n_grid_batches_esp: ', n_grid_batches_esp

         grid_batches_per_task = 0
         grid_batches_per_task(myid+1) = n_grid_batches_esp
         call sync_integer_vector(grid_batches_per_task, n_tasks)

         i_offset = SUM(grid_batches_per_task(1:myid))
         grid_partition = grid_partition + i_offset

         n_grid_batches_esp = SUM(grid_batches_per_task)
!     print *, 'MYID: ', myid, 'n_grid_batches_esp, NOW: ',n_grid_batches_esp
         
      else
         
         do i_my_point = 1, n_my_points, 1
            coord_current(:) = all_coords(:,i_my_point)

            if(n_periodic > 0)then
               call map_to_center_cell(coord_current(1:3) )
            end if
            
            idx = CEILING((coord_current - coord_min)/h)
            idx = MIN(idx, n_bins)
            idx = MAX(idx, 1)
            
            grid_partition(i_my_point) = idx(1) + (idx(2)-1)*n_bins(1) + &
            (idx(3)-1)*n_bins(2)*n_bins(1)
            
         end do

      end if
    end subroutine hash_points
!******
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/find_bb_dimensions
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
  if (equal_grid)then
    do i_x = 0, n_x-1, 1 
    do i_y = 0, n_y-1, 1
    do i_z = 0, n_z-1, 1
      d_r(1) = dble(i_x)/n_x
      d_r(2) = dble(i_y)/n_y
      d_r(3) = dble(i_z)/n_z
	  coord_new(1:3) = matmul(lattice_vector_new,d_r(1:3)) + off(1:3)
          if(sum(radius_esp_min_new).gt.0)then
	    call check_point(rm,lattice_vector,radius_esp_min_new,&
			coord_new,counter_2)
          else
            if (n_periodic > 0) then
             call map_to_center_cell(coord_new(1:3))
            endif
            counter_2=(cells_atoms+1)
          endif
      if (counter_2.eq.(cells_atoms+1))then
               i_point = i_point + 1
               !if (myid == MOD(i_point, n_tasks)) then
                  coord_current(:) = coord_new(:)
                  
                  if(n_periodic > 0)then
                     call map_to_center_cell( coord_current(1:3) )
                  end if
                  
                  x_max_loc = MAX(x_max_loc, coord_current(1))
                  x_min_loc = MIN(x_min_loc, coord_current(1))
                  y_max_loc = MAX(y_max_loc, coord_current(2))
                  y_min_loc = MIN(y_min_loc, coord_current(2))
                  z_max_loc = MAX(z_max_loc, coord_current(3))
                  z_min_loc = MIN(z_min_loc, coord_current(3))
               !end if
      endif
    enddo
    enddo
    enddo
      coord_max(1) = x_max_loc
      coord_max(2) = y_max_loc
      coord_max(3) = z_max_loc
      coord_min(1) = x_min_loc
      coord_min(2) = y_min_loc
      coord_min(3) = z_min_loc
    else
      do i_atom = 1, n_atoms, 1
	do i_radial = 1, n_radial_esp(species(i_atom)), 1
	  do i_angular = 1, n_angular_esp( i_radial,species(i_atom) ), 1
	    coord_new(1:3)=coords_center( 1:3,i_atom )+ &
	    r_angular_esp(1:3, i_angular, i_radial, species(i_atom)) * &
	    r_radial_esp(i_radial, species(i_atom))
            call check_point_rad(rm,lattice_vector,radius_esp_min_new,&
                       coord_new,counter_2)
	    if (counter_2.eq.cells_atoms)then
               i_point = i_point + 1
               if (myid == MOD(i_point, n_tasks)) then
                  coord_current(:) = coord_new(:)
                  
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
	    endif
	  enddo
	enddo
      enddo
      call get_max_double(coord_max(1), x_max_loc)
      call get_max_double(coord_max(2), y_max_loc)
      call get_max_double(coord_max(3), z_max_loc)
      call get_min_double(coord_min(1), x_min_loc)
      call get_min_double(coord_min(2), y_min_loc)
      call get_min_double(coord_min(3), z_min_loc)
      endif
      
      bb_dim = coord_max - coord_min

    end subroutine find_bb_dimensions
!******
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/find_bins
!  NAME
!    find_bins
!  SYNOPSIS
    subroutine find_bins(bb_dim, n_grid_batches_esp, n_bins, h)
!  PURPOSE
!    Computes the number of bins in x, y, and and the corresponding bin 
!    dimensions
!  USES
      implicit none
!  ARGUMENTS      
      integer, dimension(3) :: n_bins
      integer :: n_grid_batches_esp
      real*8, dimension(3) :: h, bb_dim
!  INPUTS
!    o bb_dim -- bounding box dimension in x, y, and z
!  OUTPUT
!    o n_grid_batches_esp -- number of grid batches if simple hashing is used
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
              NINT( EXP(1.0/3.0*LOG((bb_dim(1)/bb_dim(2))*(bb_dim(1)/bb_dim(3))&
              *dble(n_batches))))
         n_bins(2) = &
              NINT( EXP(1.0/3.0*LOG((bb_dim(2)/bb_dim(3))*(bb_dim(2)/bb_dim(1))&
              *dble(n_batches))))
         n_bins(3) = &
              NINT( EXP(1.0/3.0*LOG((bb_dim(3)/bb_dim(2))*(bb_dim(3)/bb_dim(1))&
              *dble(n_batches))))

!!$         if (myid == 0) then
!!$            print *, 'PRE'
!!$            print *, n_tasks
!!$            print *, n_bins
!!$            print *, bb_dim
!!$            print *, coord_min
!!$            print *, coord_max
!!$         end if

         if (PRODUCT(n_bins(:)) /= n_tasks) then

            optimal_value = 1.0d9
            do div1 = 0, n_factors
               do div2 = div1, n_factors
                  n_x = PRODUCT(prime_factors(1:div1))
                  n_y = PRODUCT(prime_factors(div1+1:div2))
                  n_z = n_tasks/(n_x*n_y)
                  trial_value = (bb_dim(1)/dble(n_x) - bb_dim(2)/dble(n_y))**2+&
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

         n_batches = NINT(dble(n_points_in_grid)/dble(n_points_in_batch_esp))
         n_bins(1) = &
              NINT( ((bb_dim(1))/(bb_dim(2))*(bb_dim(1))/(bb_dim(3))*&
              dble(n_batches))**(1.0d0/3.0d0) )
         n_bins(2) = &
              NINT( ((bb_dim(2))/(bb_dim(3))*(bb_dim(2))/(bb_dim(1))*&
              dble(n_batches))**(1.0d0/3.0d0) )
         n_bins(3) = &
              NINT( ((bb_dim(3))/(bb_dim(2))*(bb_dim(3))/(bb_dim(1))*&
              dble(n_batches))**(1.0d0/3.0d0) )
      end if

      n_grid_batches_esp = n_bins(1)*n_bins(2)*n_bins(3)
      
      h = bb_dim / dble(n_bins)

    end subroutine find_bins
!******
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/find_prime_factors
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
!    o prime_factors -- optional, prime factors of the integer n, must be 
!                       allocated before call
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
     subroutine distribute_batch_tasks_esp( batch_sizes)
!  PURPOSE
!    Distributes the grid batches to different tasks by balancing the workload.
!  USES
!
!  ARGUMENTS
      integer :: batch_sizes( n_grid_batches_esp )
!  INPUTS
!    o batch_sizes -- number of grid points in each of the grid batches
!  OUTPUT
!    none -- at exit the array batch_task_list is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: i_batch, index
      integer :: batch_task_work_table(n_tasks)
      character*100 :: info_str

      if (.not.allocated(batch_task_list_esp)) then
         allocate(batch_task_list_esp(n_grid_batches_esp))
      end if

      batch_task_list_esp = 0
      n_full_points_esp = 0

      batch_task_work_table = 0

      n_my_batches_esp = 0

      do i_batch = 1, n_grid_batches_esp, 1

         index = MINLOC(batch_task_work_table,1)
         batch_task_list_esp(i_batch) = index - 1
         batch_task_work_table(index) = &
              batch_task_work_table(index) + batch_sizes(i_batch)

         if (myid.eq.batch_task_list_esp(i_batch)) then
            n_my_batches_esp = n_my_batches_esp + 1
            n_full_points_esp = n_full_points_esp + batch_sizes(i_batch)
         end if

      end do

      if (use_mpi) then
!         if (myid.eq.0) then
            write(info_str,'(A)') ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A,A,I5,A)') &
            "Integration load balanced ", &
            "across ", n_tasks, " MPI tasks. "
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') &
                 "Work distribution over tasks is as follows:"
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            if (output_priority .le. OL_low) then 
               do index = 1, n_tasks, 1
                  write(info_str,'(2X,A,I5,A,I10,A)') &
                       "Task ", index-1, &
                       " has ", &
                       batch_task_work_table(index), &
                       " integration points."
                  call localorb_info(info_str,use_unit,'(A)',OL_low)
               end do
            end if
!         end if
      end if

      end subroutine distribute_batch_tasks_esp
!****** 
!------------------------------------------------------------------------------
!****s* esp_partition_grid/distribute_batch_tasks_by_location_esp
!  NAME
!    distribute_batch_tasks_by_location_esp
!  SYNOPSIS
      subroutine distribute_batch_tasks_by_location_esp(batch_sizes,coords)
!  PURPOSE
!    Distributes the grid batches to different tasks based on their location.
!  USES
!
!  ARGUMENTS 
      integer :: batch_sizes( n_grid_batches_esp )
      real*8  :: coords(3,n_grid_batches_esp)
!  INPUTS
!    o batch_sizes -- number of grid points in each of the grid batches
!    o coords -- coordinates of the grid batches
!  OUTPUT
!    none -- at exit the array batch_task_list_esp is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


      real*8, allocatable :: batch_desc(:,:)
      integer :: i_batch, my_off, my_len
      character*100 :: info_str

      if (.not.allocated(batch_task_list_esp)) then
         allocate(batch_task_list_esp(n_grid_batches_esp))
      end if

      batch_task_list_esp = 0

      ! Set up batch description

      allocate(batch_desc(5,n_grid_batches_esp))

      do i_batch=1,n_grid_batches_esp
         batch_desc(1:3,i_batch) = coords(1:3,i_batch)
         batch_desc(4,i_batch)   = batch_sizes(i_batch)
         batch_desc(5,i_batch)   = i_batch
      enddo

      ! Distribute batches
      my_off = 0
      my_len = n_grid_batches_esp
      call distribute_batches_by_location_esp(batch_desc, 0, n_tasks-1, &
           my_off, my_len)

      ! Set my batches in batch_task_list_esp

      n_full_points_esp = 0
      n_my_batches_esp = my_len

      batch_task_list_esp(:) = 0

      do i_batch = my_off+1, my_off+my_len

         batch_task_list_esp(nint(batch_desc(5,i_batch))) = myid
         n_full_points_esp = n_full_points_esp + batch_sizes(&
         nint(batch_desc(5,i_batch)))

      end do

      call sync_integer_vector(batch_task_list_esp, n_grid_batches_esp)

      if (use_mpi) then
         write(info_str,'(A)') ''
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,A,I6,A)') &
              "Integration load balanced across ", n_tasks, " MPI tasks."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,A)') "Work distribution over tasks is as follows:"
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,A,I6,A,I10,A)') &
              "Task ", myid, " has ", n_full_points_esp, " integration points."
         call localorb_allinfo(info_str,use_unit,'(A)',OL_low)
      end if

      deallocate(batch_desc)

      end subroutine distribute_batch_tasks_by_location_esp
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/distribute_batches_by_location_esp
!  NAME
!    distribute_batches_by_location_esp
!  SYNOPSIS
      recursive subroutine distribute_batches_by_location_esp &
                           (batch_desc,cpu_a,cpu_b,my_off,my_len)
!  PURPOSE
!    Recursively distributes the batches to the tasks by their location
!  USES

      implicit none
!  ARGUMENTS 
      real*8, intent(inout) :: batch_desc(:,:)
      integer, intent(in) :: cpu_a, cpu_b
      integer, intent(inout) :: my_off, my_len
!  INPUTS
!    o batch_desc(1,.) -- x-coord of batch (typically coord of some point in 
!                         batch)
!    o batch_desc(2,.) -- y-coord of batch
!    o batch_desc(3,.) -- z-coord of batch
!    o batch_desc(4,.) -- weight of batch
!    o batch_desc(5,.) -- number of batch, must be 1,2,3... on top level call
!    o cpu_a -- number of first CPU for distribution
!    o cpu_b -- number of last CPU for distribution
!    o my_off -- offset of batches on current CPU set (cpu_a to cpu_b) within 
!                batch_desc
!                must be 0 on top level call 
!    o my_len -- number of batches on current CPU set
!                must be the number of batches on top level call
!  OUTPUT
!    o my_off -- offset of my batches on my CPU within batch_desc
!    o my_len -- number of my batches
!    o batch_desc -- is sorted so that batch_desc(:,my_off+1:my_off+my_len) 
!                    contains the corresponding values for my batches.
!                    Normally only batch_desc(5,:) is relevant on output.
!                    Values outside my_off+1:my_off+my_len are undefined and 
!                    should not be used!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      real*8 total_weight, target_weight
      integer median, cpu_m, n, icoord
      real*8 xmax(3), xmin(3), d

      ! If there is only 1 CPU for distribution, we are done

      if(cpu_a==cpu_b) return

      ! Get geometric extensions and total number of points of all batches on 
      ! current CPU set

      total_weight = 0
      xmax = -1.d99
      xmin =  1.d99
      do n=my_off+1,my_off+my_len
         xmax(:) = max(xmax(:),batch_desc(1:3,n))
         xmin(:) = min(xmin(:),batch_desc(1:3,n))
         total_weight = total_weight+batch_desc(4,n)
      enddo

      ! Get dimension with biggest distance from min to max

      d = 0
      icoord = 1
      do n=1,3
         if(xmax(n)-xmin(n) > d) then
            d = xmax(n)-xmin(n)
            icoord = n
         endif
      enddo

      ! CPU number where to split CPU set

      cpu_m = (cpu_a+cpu_b-1)/2

      ! If the number of CPUs is not even, we have to split the set
      ! of batches accordingly into to differently sized halfes
      ! in order to end with an equal weight on every CPU.

      target_weight = total_weight*(cpu_m-cpu_a+1) / (cpu_b-cpu_a+1)

      ! divide in two equal halfes with the weight of the left half approx.
      ! equal to target_weight

      call divide_values_esp(batch_desc(:,my_off+1:my_off+my_len), 5, my_len, &
                         icoord, &
                         target_weight, median)

      ! Set my_off and my_len for next distribution step

      if(myid<=cpu_m) then
         ! my_off stays the same
         my_len = median
      else
         my_off = my_off+median
         my_len = my_len-median
      endif

      ! If there are only two CPUs, we are done

      if(cpu_b == cpu_a+1) return

      ! Further divide both halves recursively

      if(myid<=cpu_m) then
         call distribute_batches_by_location_esp(batch_desc,cpu_a,  &
              cpu_m,my_off,my_len)
      else
         call distribute_batches_by_location_esp(batch_desc,cpu_m+1,cpu_b,&
              my_off,my_len)
      endif


      end subroutine distribute_batches_by_location_esp
!******
!-------------------------------------------------------------------------------
!****s* esp_partition_grid/divide_values
!  NAME
!    divide_values
!  SYNOPSIS
   subroutine divide_values_esp(x, ldx, num, icoord, target, ndiv)
!  PURPOSE
!    Divides an array x(:,num) so that
!    - all x(icoord,1:ndiv) are smaller than all x(icoord,ndiv+1:num)
!    - SUM(x(4,1:ndiv)) == target (approximatly)
!  USES
    use mpi_tasks, only: aims_stop
    implicit none
!  ARGUMENTS 
    real*8, intent(inout) :: x(ldx,num)
    integer, intent(in)   :: ldx, num, icoord
    real*8, intent(in)    :: target
    integer, intent(out)  :: ndiv
!  INPUTS
!    o x -- on entry contains unsorted coords, weights, additional values
!    o ldx -- leading dimension of x
!    o num -- number of data points in x
!    o icoord -- number of coord (row in x) to be used for sorting
!    o target -- target weight for left part of sorted values
!  OUTPUT
!    o x -- on exit contains sorted values as described under PURPOSE
!    o ndiv -- number of values in left (lower) part
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2011).
!  SOURCE

    integer, parameter :: iwgt = 4 ! where to get the weight from

    integer ns, ne, nl, nr, i
    real*8 xmin, xmax, xmid, sl_tot, sr_tot, sl, sr, xmax_l, xmin_r
    real*8 xtmp(ldx)

    ! Safety check only
    if(num<=0) then
      print *,'Error divide_values: num = ',num
      call aims_stop ('No data points','divide_values_esp')
    endif

    xmin = minval(x(icoord,:))
    xmax = maxval(x(icoord,:))

    ns = 1
    ne = num

    sl_tot = 0
    sr_tot = 0

    do

      xmid = 0.5*(xmin+xmax)
      ! Make sure that at least 1 value is below xmid and 1 value is above xmid.
      ! Otherways xmin and xmax are equal or only 1 bit apart and there is no 
      ! more separation possible
      if(xmid<=xmin .or. xmid>=xmax) exit

      ! divide values between ns..ne into values <= xmid and values > xmid
      ! with a quicksort like algorithm

      nl = ns
      nr = ne
      sl = 0
      sr = 0

      xmax_l = xmin
      xmin_r = xmax

      do while(nl<=nr)

        ! skip values in front <= xmid

        do while(nl <= nr)
          if(x(icoord,nl) > xmid) exit
          sl = sl + x(iwgt,nl)
          xmax_l = max(xmax_l,x(icoord,nl))
          nl = nl+1
        enddo

        ! skip values in back > xmid

        do while(nl <= nr)
          if(x(icoord,nr) <= xmid) exit
          sr = sr + x(iwgt,nr)
          xmin_r = min(xmin_r,x(icoord,nr))
          nr = nr-1
        enddo

        if(nl>nr) exit

        ! Exchange elements at nr/nl

        xtmp(:) = x(:,nl)

        x(:,nl) = x(:,nr)
        sl = sl + x(iwgt,nl)
        xmax_l = max(xmax_l,x(icoord,nl))
        nl = nl+1

        x(:,nr) = xtmp(:)
        sr = sr + x(iwgt,nr)
        xmin_r = min(xmin_r,x(icoord,nr))
        nr = nr-1

      enddo

      ! Safety check:
      ! Check that at least 1 value is in both halves, this must be always the 
      ! case! Otherways there is something screwed up with the program logic and
      ! we better exit
      if(nl==ns .or. nr==ne) then
        print *,'INTERNAL error in program logic of divide_values'
        call aims_stop ('INTERNAL ERROR','divide_values_esp')
      endif

      ! we can keep one half of the sorted values as is whereas the other half 
      ! has to be sorted again

      if(sl_tot+sl < target) then
        ! Left is ok, right must be sorted again
        ns = nl
        sl_tot = sl_tot+sl
        xmin = xmin_r
      else
        ! Right is ok, left must be sorted again
        ne = nr
        sr_tot = sr_tot+sr
        xmax = xmax_l
      endif

    enddo

    ! Safety check
    if(ns>ne) then
      print *,'INTERNAL error in divide_values: ns=',ns,' ne=',ne
      call aims_stop ('INTERNAL ERROR','divide_values_esp')
    endif

    ! Now the value searched must be somewhere between ns an ne
    ! Please note that the coords of x(icoord,ns:ne) are (nearly) the same

    do i = ns, ne
      sl_tot = sl_tot + x(iwgt,i)
      if(sl_tot >= target) then
        ndiv = i
        return
      endif
    enddo

    ! Well, we should never come here unless there is something wrong with the 
    ! weights ...

    ndiv = ne

  end subroutine divide_values_esp
!******
end subroutine esp_partition_grid
