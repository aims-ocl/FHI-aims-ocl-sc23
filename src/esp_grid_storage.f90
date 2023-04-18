!****s* FHI-aims/esp_grid_storage
!  NAME
!   esp_grid_storage
!  SYNOPSIS

    subroutine esp_grid_storage( partition_tab, free_hartree_superpos_trans, &
                                 equal_grid )

!  PURPOSE
!   calculates relevant quantities for integration on the grid, as well as 
!   initializing a number of free atom references for the calculation of esp 
!   charges.
!
!  USES
!
      use dimensions,only:n_species,n_multipoles,n_centers_basis_integrals,&
                          n_spin,n_periodic
      use esp_grids,only:n_my_batches_esp,batches_esp,invert_log_grid_esp,&
                         w_radial_esp,w_angular_esp,&
                         r_grid_min_esp,r_grid_inc_esp,n_full_points_esp
      use grids,only:n_grid
      use runtime_choices,only:stratmann_a
      use geometry,only:multipole_coords,empty,species
      use spline,only:val_spline
      use localorb_io,only:OL_high,OL_low,localorb_info,use_unit
      use free_atoms,only:hartree_partition_rho_spl
      use species_data,only:cut_free_atom,multipole_radius_free_sq,&
                            outer_partition_radius,free_r_cut,w_cutoff,&
                            multipole_radius_free
      use pbc_lists,only:coords_center,centers_basis_integrals,species_center,&
                         center_to_atom
      use synchronize_mpi,only:MPI_INTEGER,MPI_SUM,&
                               MPI_COMM_GLOBAL,sync_workload
      use mpi_tasks,only:check_allocation,myid
      use constants,only:bohr
      use physics, only: weight_tab


      implicit none

!  ARGUMENTS
      real*8, dimension(n_full_points_esp),INTENT(OUT) :: partition_tab
      real*8, dimension(n_full_points_esp),INTENT(OUT) :: free_hartree_superpos_trans
      logical,INTENT(IN)  :: equal_grid

! INPUTS
!   none
! OUTPUTS
! o partition_tab -- grid integration weight
! o free_hartree_superpos_trans -- Hartree potential of free atoms
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
!  SOURCE
!


!     local variables
      real*8, dimension(3) :: coord_current
      real*8 :: dist_tab_sq(n_centers_basis_integrals)
      real*8 :: dist_tab(n_centers_basis_integrals)
      real*8 :: dir_tab(3,n_centers_basis_integrals)
      real*8 :: dir_tab_norm(3,n_centers_basis_integrals)
      real*8 :: i_r(n_centers_basis_integrals)
      integer :: n_compute_atoms
      integer, dimension(n_centers_basis_integrals) :: center_index
      real*8, dimension(n_centers_basis_integrals) :: temp_free_rho
      integer :: n_compute_occ_atoms
      integer, dimension(n_centers_basis_integrals) :: i_occ2i_compute
      real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho
      real*8, dimension(:, :), allocatable :: atom_atom_tab
      real*8, dimension(:),    allocatable :: min_atom_atom_tab
      integer, dimension(:),   allocatable :: atom_atom_index
      real*8 :: aux_spl
      integer :: n_atom_atom_tab
      real*8, dimension(:), allocatable :: atom_atom_dist_list
      integer, dimension(:), allocatable :: atom_idx_A
      integer, dimension(:), allocatable :: atom_idx_B

      logical :: point_on_atom, points_on_multipole

      character*140 :: info_str

!     counters
      integer :: i_atom, i_species, i_point,  &
           i_atom_2, i_index, i_multipole
      integer :: current_atom, current_radial, current_angular
      integer :: i_my_batch
      integer :: i_center_L, i_center, i_center_L2

      real*8, dimension(3) :: dir_current
      real*8 :: dist_current, dist_current_sq, dens_current, r_temp
      real*8 ::  dist_to_multipole_sq

      real*8 :: free_rho_superpos_tmp
      real*8 :: free_rho_gradient_superpos_tmp(3)
      real*8 :: rho_tmp(n_spin)
      real*8 :: local_rho_gradient(3,n_spin)
      !SR: laplacian of nfree
      !real*8 :: dummy

      integer :: n_int_points,n_full_points_total,n_int_points_total
  n_full_points_total = 0
  n_int_points_total = 0
  rho_tmp = 0d0
  weight_tab = 0d0
      write (info_str,'(2X,A,A)') & 
      "Initializing partition tables, free-atom densities, potentials, etc. ", &
      "across the integration grid (initialize_grid_storage)."
      call localorb_info(info_str,use_unit,'(A)',OL_high)

      n_int_points = 0
      i_point = 0

      ! We begin by preparing the partition table of Stratmann and coworkers. 
      ! This part became lengthy over time, please simply take it as one block 
      ! (the end of this part is also marked by a comment below)

        ! All these are variants of the partition table by Stratmann and 
        ! coworkers. Incremental improvements during development, keeping them 
        ! does not hurt. 

        ! Determine what will be the outermost distance at which the stratmann 
        ! partition table for a given atom should be non-zero

        ! outer_partition_radius is the outermost radius at which any integrand
        ! or density component will be partitioned to belong to a given atom - 
        ! splease ee Delley (1990) for a very natural and intuitive explanation 
        ! of the partitioning procedure that is used for the integrals and 
        ! densities in codes like DMol or FHI-aims.

        if (cut_free_atom(1)) then  ! there has to be at least one species. 
                                    ! Referring to cut_free_atom(1) is an ugly 
                                    !hack, but the output is explanatory only 
                                    ! anyway.
            write (info_str,'(2X,A)') &
              "| initialize_grid_storage: Actual outermost partition radius vs. multipole_radius_free"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            write (info_str,'(2X,A)') &
              "| (-- VB: in principle, multipole_radius_free should be larger, hence this output)"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
        end if

        do i_species = 1, n_species, 1
          if (cut_free_atom(i_species)) then
            ! This is the case where the free-atom radius is bounded somehow -
            ! should be the normal case. If the free atom radius is not bounded,
            ! we have a problem for periodic systems - very extended 
            ! integration regions must be accounted for.
            outer_partition_radius(i_species) = free_r_cut (i_species) + &
            w_cutoff (i_species)
            ! ... and add some diagnostic output here.
            write (info_str,'(2X,A,I8,A,F30.15,A,F30.15,A)') &
              "| Species ", i_species, ": outer_partition_radius = ", &
              outer_partition_radius(i_species)*bohr, &
              " AA, multipole_radius_free = ", multipole_radius_free(i_species)&
              *bohr, " AA."
            call localorb_info(info_str,use_unit,'(A)',OL_low)
          else
            ! The free-atom radius is not bounded. Well, we must then use the 
            ! outermost radius of the free atom density (very large), for 
            ! better or for worse.
            outer_partition_radius(i_species) = multipole_radius_free(i_species)
          end if
        end do


        ! for stratmann partitioning scheme, make available only if necessary
        ! tabulate the interatomic distances for all relevant atoms
        allocate(atom_atom_tab(n_centers_basis_integrals,&
                 n_centers_basis_integrals),stat=i_point)
        call check_allocation(i_point, 'atom_atom_tab                 ')

      
        ! Not needed, allocate dummy
        n_atom_atom_tab = 1
        allocate(atom_atom_dist_list(1))
        allocate(atom_idx_A(1))
        allocate(atom_idx_B(1))


      allocate(atom_atom_index(n_centers_basis_integrals),stat=i_point)
      call check_allocation(i_point, 'atom_atom_index               ')

      allocate(min_atom_atom_tab(n_centers_basis_integrals),stat=i_point)
      call check_allocation(i_point, 'min_atom_atom_tab             ')


        ! FIXME:
        ! (1) count distance^2 and figure out which atoms could potentially be 
        !     important in each other's stratmann partition tab
        ! (2) calculate min_distance^2
        ! (3) only allocate as much data as is required, otherwise the arrays 
        !     may get too big ...
        ! (4) tabulate interatomic distances ONLY for those atoms needed.
        ! (5) remember how many relevant atoms there are for each possible 
        !     i_atom
        ! the resulting list should be enough to build the stratmann list from 
        ! it...
        call tab_interatomic_distances(n_centers_basis_integrals, &
             centers_basis_integrals, atom_atom_tab)
        ! for each atom, determine distance to next neighbour:
        do i_atom = 1, n_centers_basis_integrals
          min_atom_atom_tab(i_atom) = 1d100
          do i_atom_2 = 1, n_centers_basis_integrals
            if (i_atom.ne.i_atom_2) then
              min_atom_atom_tab(i_atom) = min(min_atom_atom_tab(i_atom),&
              atom_atom_tab(i_atom,i_atom_2))
            end if
          end do
        end do
        ! need this for comparison ...
        min_atom_atom_tab(:) = (1d0-stratmann_a)*min_atom_atom_tab(:)/2d0

      ! End all preparations for the Stratmann and coworkers partition table.


      do i_my_batch = 1, n_my_batches_esp, 1

        do i_index = 1, batches_esp(i_my_batch)%size_esp, 1

          i_point = i_point + 1


          coord_current(:) = batches_esp(i_my_batch) % points_esp(i_index) % &
                             coords_esp(:)
                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current(1:3) )
                 end if
if(.not.equal_grid)then
          current_atom    = batches_esp(i_my_batch) % points_esp(i_index) %  &
                            index_atom_esp
          current_radial  = batches_esp(i_my_batch) % points_esp(i_index) %  &
                            index_radial_esp
          current_angular = batches_esp(i_my_batch) % points_esp(i_index) %  &
                            index_angular_esp
          ! remember where this point is with respect to its own atom as well
          ! as its free atom density
          dir_current(:)  = coord_current(:)- &
                            coords_center(:,&
                            centers_basis_integrals(current_atom))
          dist_current_sq = dir_current(1)*dir_current(1) &
                          + dir_current(2)*dir_current(2) &
                          + dir_current(3)*dir_current(3)
          dist_current    = sqrt(dist_current_sq)
          dir_current(:)  = dir_current(:)/dist_current
          r_temp          = invert_log_grid_esp &
                            ( dist_current, &
                            r_grid_min_esp(species_center(&
                            centers_basis_integrals(current_atom))), &
                            r_grid_inc_esp(species_center(&
                            centers_basis_integrals(current_atom))))

          dens_current    =  val_spline &
                ( r_temp, hartree_partition_rho_spl(1,1,species_center(&
                  centers_basis_integrals(current_atom))), &
                  n_grid(species_center(centers_basis_integrals(current_atom))))
endif
!              tabulate current integration point as it appears from spherical
!              coordinates centered at each atom

          call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, abs(n_centers_basis_integrals), centers_basis_integrals )

          ! To avoid a floating-point exception, check here whether this 
          ! integration point happens to sit on another atom.
          ! VB: We could also check here whether we are inside the innermost 
          !     logarithmic grid shell of an atom, as is done further down in 
          !     evaluate_partition_tab. This would be even better, since we
          !     could then remove the same check fromevaluate_partition 
          !     altogether ...
          point_on_atom = .false.
          do i_center_L = 1, n_centers_basis_integrals, 1
            if ( dist_tab_sq(i_center_L).eq.0.d0) then
              point_on_atom = .true.
              exit ! exit do loop
            end if
          end do

          ! also check for multipole singularities
          points_on_multipole = .false.
          do i_multipole = 1, n_multipoles
            dist_to_multipole_sq = &
               (coord_current(1) - multipole_coords(1, i_multipole))**2 &
              +(coord_current(2) - multipole_coords(2, i_multipole))**2 &
              +(coord_current(3) - multipole_coords(3, i_multipole))**2 
            if(dist_to_multipole_sq.eq.0.d0) then
              points_on_multipole = .true.
              exit
            endif
          enddo
          if (.not. point_on_atom .and. .not.points_on_multipole) then
            ! This is the normal case. If our grid point does not sit on an 
            ! atom, we keep it for later use.

            ! For the following operations, we only need those atoms that have 
            ! a non-zero free-atom charge density at the current point. We 
            ! assemble that list here explicitly ...
            n_compute_atoms = 0
            n_compute_occ_atoms = 0

            do i_center_L = 1, n_centers_basis_integrals, 1
              i_center = centers_basis_integrals(i_center_L)
if(.not.equal_grid)then
              if ((i_center.eq.current_atom).or. &
                  (dist_tab_sq(i_center_L).lt.multipole_radius_free_sq(&
                                               species_center(i_center)) )) then
                  ! this center has a non-zero free-atom density, or it belongs
                  ! to current integration point
                  n_compute_atoms                  = n_compute_atoms + 1
                  center_index(n_compute_atoms)    = i_center
                  dist_tab_sq(n_compute_atoms)     = dist_tab_sq(i_center_L)
                  dir_tab(:,n_compute_atoms)       = dir_tab(:,i_center_L)   
                  ! This ensures consistent handling later, 
                  ! note n_compute_atoms <= i_center_L
                  atom_atom_index(n_compute_atoms) = i_center_L             
                  ! indexing for later use of atom_atom_tab, which is NOT
                  ! recomputed here for speed reasons
              end if
else
                  n_compute_atoms                  = n_compute_atoms + 1
                  center_index(n_compute_atoms)    = i_center
                  dist_tab_sq(n_compute_atoms)     = dist_tab_sq(i_center_L)
                  dir_tab(:,n_compute_atoms)       = dir_tab(:,i_center_L)   
                  atom_atom_index(n_compute_atoms) = i_center_L   
endif
            end do

            call tab_global_geometry_p0 &
                ( dist_tab_sq,         &
                  dir_tab,             &
                  dist_tab,            &
                  i_r,                 &
                  dir_tab_norm,        &
                  n_compute_atoms,     &
                  center_index )
            ! calculate the free-atom density only for the (now) known atoms ...
            do i_center_L2 = 1, n_compute_atoms
              i_center = center_index(i_center_L2)
if(.not.equal_grid)then
              if (i_center.eq.current_atom) i_center_L = i_center_L2 
endif
                  ! remember the center we are currently at!
              aux_spl = val_spline &
                  ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,&
                    species_center(i_center)), &
                    n_grid(species_center(i_center)) )
              temp_free_rho(i_center_L2) = aux_spl

              if (.not.empty(center_to_atom(i_center))) then
                  n_compute_occ_atoms = n_compute_occ_atoms + 1
                  i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
                  temp_occ_rho(n_compute_occ_atoms) = aux_spl
              end if

            end do
        if(.not.equal_grid)then
              call evaluate_partition_tab_p2  &
                ( current_atom,              &
                  i_center_L,                &
                  dist_current,              &
                  dist_current_sq,           &
                  dens_current,              &
                  dist_tab,                  &
                  i_r,                       &
                  w_radial_esp(current_radial, species(current_atom)),  &
                  w_angular_esp(current_angular, current_radial, &
                                             species(current_atom)),  &
                  partition_tab(i_point),    &
                  weight_tab(i_point),       &
                  n_centers_basis_integrals, &
                  n_compute_atoms,           &
                  center_index,              &
                  temp_free_rho,             &
                  dist_tab_sq,               &
                  atom_atom_tab,             &
                  atom_atom_index,           &
                  min_atom_atom_tab,         &
                  n_atom_atom_tab,           &
                  atom_atom_dist_list,          &
                  atom_idx_A,                &
                  atom_idx_B)
          endif


!            write (81,'(2X,F15.8,2X,F15.8,2X,F15.8,2X,F15.8,2X, I3,A)')  &
!                   coord_current(1), coord_current(2), coord_current(3), &
!              partition_tab(i_point),current_atom, " coord, partition_tab"
            call evaluate_free_atom_sums_p2  &
                ( i_r, dir_tab_norm,  &
                free_hartree_superpos_trans(i_point),  &
                free_rho_superpos_tmp,  &
                free_rho_gradient_superpos_tmp, &
                rho_tmp,  &
                local_rho_gradient, n_compute_atoms, center_index, &
                n_compute_occ_atoms, i_occ2i_compute, temp_occ_rho) !SR: for laplacian of free el. dens: ,dist_tab_sq,&
                !dummy)
            if(.not.equal_grid)then
              if (partition_tab(i_point).gt.0) then
                n_int_points = n_int_points + 1
              end if
            else
              n_int_points = n_int_points + 1
            endif
          else
          ! this is the case where the current integration point sits exactly 
          ! on top of another atom must never be considered in later dealings!
            if(.not.equal_grid)then
              partition_tab(i_point) = 0.d0
            endif
          end if

        end do ! loop over batch
      end do ! loop over batches
      call sync_workload( n_full_points_esp, n_int_points,   &
           n_full_points_total, n_int_points_total )
      if(equal_grid)then
         n_int_points_total = n_full_points_total
      endif
      write (info_str,'(2X,A,I8)')  &
              "| Net number of integration points: ",   &
              n_full_points_total
      call localorb_info(info_str,use_unit,'(A)', OL_high)
      write (info_str,'(2X,A,I8)')  &
             "| of which are non-zero points    : ", n_int_points_total
      call localorb_info(info_str,use_unit,'(A)', OL_high)


      ! clean up the odd allocated variable
      deallocate(min_atom_atom_tab)
      deallocate(atom_atom_index)
      deallocate(atom_atom_tab)
      deallocate(atom_atom_dist_list)
      deallocate(atom_idx_A)
      deallocate(atom_idx_B)
! test
!      close(80)
!      close(81)
! test end

    end subroutine esp_grid_storage
!******
