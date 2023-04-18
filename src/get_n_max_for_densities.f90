!****s* FHI-aims/get_n_max_for_densities
!  NAME
!    get_n_max_for_densities
!  SYNOPSIS 
subroutine get_n_max_for_densities( partition_tab, hartree_partition_tab)
!  PURPOSE
!    Computes the maximum (density related) dimensions for the hirshfeld evaluations
!  USES
  use dimensions
  use grids
  use runtime_choices
  use localorb_io
  use mpi_utilities
  use synchronize_mpi
  use species_data
  use pbc_lists

  implicit none
!  ARGUMENTS
  real*8, dimension(n_full_points) :: partition_tab
  real*8, dimension(n_full_points) :: hartree_partition_tab
!  INPUTS
!    o partition_tab -- the partition tab
!  OUTPUT
!    n_max_compute_missing_dens, n_max_compute_fns_missing_dens, n_max_compute_missing_atoms and
!    n_avg_compute_missing_dens
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


 ! real*8, dimension( n_full_points_hamiltonian_integrals) :: hamiltonian_partition_tab
  
  ! locals
  integer :: i_center, i_center_L
  integer :: i_my_batch, i_index
  integer :: i_point
  logical :: point_on_atom

  integer :: n_compute_c, n_compute_a
  integer :: i_full_points_A, i_full_points_C, i_full_points_2C
  integer :: i_full_points, i_full_points_2, i_full_points_3

  real*8 :: coord_current(3)
  real*8 :: dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 :: dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 :: dir_tab(3,n_centers_integrals, n_max_batch_size)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

  integer :: n_centers_max, n_points

  integer,dimension(:),allocatable :: i_basis

  character*100 :: info_str
  logical, dimension(:), allocatable :: my_functions

  real*8, dimension(n_species) :: r_grid_min_sq



  write (info_str,'(2X,A)') "Obtaining max. number of non-zero basis functions in each batch ."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  n_centers_max = MAX(n_centers_basis_I, n_centers_basis_T)

  allocate(i_basis(n_centers_max))

  n_max_compute_missing_dens = 0
  n_max_compute_fns_missing_dens = 0

  n_avg_compute_missing_dens = 0
  n_max_compute_missing_atoms = 0


  i_full_points_C = 0
  i_full_points_2C = 0
  i_full_points_A = 0

  i_basis_fns_inv = 0

  r_grid_min_sq(:) = r_grid_min(:) * r_grid_min(:)

  if (use_metis_batch_distribution) then
     allocate(my_functions(n_centers_basis_I))
     my_functions = .false.
  end if

  do i_my_batch = 1, n_my_batches, 1


        n_compute_c = 0
        n_compute_a = 0
        i_basis = 0
     
        i_point = 0

        ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points_2C = i_full_points_2C + 1
           
!           Criteria must be like it is in update_missing_density_*
!           if (partition_tab(i_full_points_2C).gt.0.d0) then
            if (max(partition_tab(i_full_points_2C),&
                   hartree_partition_tab(i_full_points_2C)).le.0.d0) then 

              i_point = i_point+1

              ! get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
              
              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if
              
              ! compute atom-centered coordinates of current integration point,
              ! as viewed from all atoms
              call tab_atom_centered_coords_p0 &
                   ( coord_current,  &
                   dist_tab_sq(1,i_point),  &
                   dir_tab(1,1,i_point), &
                   n_centers_integrals, centers_basis_integrals )
              
              ! determine which basis functions are relevant at current integration point,
              ! and tabulate their indices
              
              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed

              call prune_basis_p0 &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_a, n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
              
           end if
        enddo  ! end loop over one part of the angular division              

        n_points = i_point
        
!        No prune_basis_once here
!        if (prune_basis_once) then
!           batches(i_my_batch)%batch_n_compute = n_compute_a
!           allocate(batches(i_my_batch)%batch_i_basis(n_compute_a))
!           batches(i_my_batch)%batch_i_basis = i_basis(1:n_compute_a)
!        end if

        if (use_metis_batch_distribution) then
           my_functions(i_basis(1:n_compute_a)) = .true.
        end if
        
        n_max_compute_missing_dens = MAX(n_compute_c, n_max_compute_missing_dens)
        n_avg_compute_missing_dens = n_avg_compute_missing_dens + dble(n_compute_c*n_points)

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_c.gt.0) then

           i_point = 0

           ! loop over one division of the angular grid
           do i_index = 1, batches(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points_C = i_full_points_C + 1
              
!              Criteria must be like it is in update_missing_density_*
!              if (partition_tab(i_full_points_C).gt.0.d0) then
               if (max(partition_tab(i_full_points_C),&
                   hartree_partition_tab(i_full_points_C)).le.0.d0) then   
              
                 i_point = i_point+1

                 n_compute_atoms = 0
                 n_compute_fns = 0
!                 i_basis_fns_inv = 0

!                Check if point on atom: that is inside the innermost log grid shell of a given atom
                 point_on_atom = .false.
                  do i_center = 1, n_centers_integrals, 1                    
                    if ( dist_tab_sq(i_center,i_point).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
                      point_on_atom = .true.
                      exit ! exit do loop
                    end if
                  end do

!                Do pruning of radial basis only if point not on atom,
!                to avoid divisions by zero
                 if (.not.point_on_atom) then
                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                 ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                 ! without any copying and without doing any unnecessary operations. 
                 ! The price is that the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                   call prune_radial_basis_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dist_tab(1,i_point), &
                      dir_tab(1,1,i_point), &
                      n_compute_atoms, atom_index, atom_index_inv, &
                      n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                      i_atom_fns, spline_array_start, spline_array_end, &
                      n_centers_integrals, centers_basis_integrals)
                 
                   n_max_compute_fns_missing_dens = MAX(n_compute_fns, n_max_compute_fns_missing_dens)
                   n_max_compute_missing_atoms = MAX(n_compute_atoms, n_max_compute_missing_atoms)
                 endif !(.not.point_on_atom)

               end if! partition tab if
 
            end do ! end loop over a batch
           
        else
               ! must increment grid counter in this case too, else we get an inconsistent
               ! partition_tab and all sorts of trouble

               i_full_points_C = i_full_points_C + batches(i_my_batch)%size

        end if ! end if (n_compute.gt.0)
        
     ! end if ! end distribution of tasks over threads
     
  end do ! end loop over batches

!  if (use_metis_batch_distribution) then
!     write(use_unit,'(2X,A,I8,A,A,I8,A,I3)') "| Using ", COUNT(my_functions), &
!          " distinct basis functions ", "out of ", n_centers_basis_I, " in task ", myid
!  end if

  call sync_n_avg_compute( n_avg_compute_missing_dens )
  n_avg_compute_missing_dens = n_avg_compute_missing_dens / dble(n_int_points_total)
  
  if( allocated(i_basis))then
     deallocate(i_basis)
  end if
  if (allocated(my_functions)) then
     deallocate(my_functions)
  end if
  
  
end subroutine get_n_max_for_densities
!******
