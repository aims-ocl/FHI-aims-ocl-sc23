!****h* FHI-aims/initialize_grid_storage
!  NAME
!   initialize_grid_storage
!  SYNOPSIS
module initialize_grid_storage

      use aims_memory_tracking, only: aims_allocate, aims_deallocate
      use dimensions
      use grids
      use runtime_choices
      use geometry
      use spline
      use free_atoms
      use pbc_lists
      use species_data
      use mpi_utilities
      use synchronize_mpi
      use localorb_io
      use pseudodata
      use physics,         only: rho_pce, rho_pce_gradient
      use rel_x2c_mod

     implicit none

     private

     public :: initialize_grid_storage_p1

contains

! **************************************************************************************************
    subroutine initialize_grid_storage_p1( &
          partition_tab, weight_tab, hartree_partition_tab, hartree_weight_tab, hartree_partition_deriv, &
          partition_deriv_delley,   &  !shanghui add 
          hartree_potential,   &
          free_hartree_superpos, free_rho_superpos,  &
          free_rho_gradient_superpos, &
          rho, rho_gradient, &
          pot_ion_embed &
    )

!  PURPOSE
!   calculates relevant quantities for integration on the grid, as well as initializing
!   a number of free atom references.
!
!   A specifically important item that is tabulated here for all grid points is the
!   so-called "partition table" or "partition_tab". There are two separate
!   partition tables, one for integrations and one for the multi-center multipole
!   decomposition of the density from which the Hartree potential is determined.
!
!   For integrations, "partition_tab" is simply the integration weight at each point.
!   To understand the definition and idea behind partition_tab, please see Eqs. (15-17)
!   of the FHI-aims CPC reference, 
!
!     Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!     Xinguo Ren, Karsten Reuter, and Matthias Scheffler
!     'Ab initio molecular simulations with numeric atom-centered orbitals'
!     Computer Physics Communications 180, 2175-2196 (2009)
!     http://dx.doi.org/10.1016/j.cpc.2009.06.022
!
!   The default in FHI-aims is a MODIFIED version of the partition table suggested
!   by Stratmann and coworkers, Ref. 25 of the FHI-aims CPC. The modification restricts
!   the maximum radius of the partition table for a given atom to a finite value,
!   otherwise periodic systems would become exceedingly difficult to handle. This
!   modification is described in words in the manual, as well as in equations in
!   a forthcoming paper on the stress tensor in FHI-aims by
!
!     Franz Knuth, Christian Carbogno, Viktor Atalla, Volker Blum, Matthias Scheffler,
!     to be submitted (2014).
!
!   For the Hartree potential multicenter multipole decomposition, see the 
!   FHI-aims CPC, Eq. (31)-(36). 
!
!   The work on the partition table in FHI-aims is based on Bernard Delley's 
!   1990 publication in DMol (Ref. 10 of the FHI-aims CPC). It is highly recommended
!   to read the description in Delley's paper - the version in the FHI-aims CPC 
!   is more condensed since Delley's writeup is already excellent. Delley's suggested
!   partition table is also accessible in FHI-aims as 'rho_r2'.
!
      use heat_flux
      implicit none
!  ARGUMENTS

      real*8, dimension(n_full_points) :: partition_tab
      real*8, dimension(n_full_points) :: weight_tab
      real*8, dimension(n_full_points) :: hartree_weight_tab
      real*8, dimension(n_full_points) :: hartree_partition_tab
      real*8, dimension(n_full_points) :: hartree_partition_deriv
      real*8, dimension(3,n_atoms,n_full_points) :: partition_deriv_delley
      real*8, dimension(n_full_points) :: hartree_potential
      real*8, dimension(n_full_points) :: free_hartree_superpos 
      real*8, dimension(n_full_points) :: free_rho_superpos
      real*8, dimension(n_full_points) :: pot_ion_embed
      real*8, dimension(n_spin,n_full_points) :: rho
      real*8, dimension(3, n_full_points) :: free_rho_gradient_superpos
      real*8, dimension(3, n_spin, n_full_points) :: rho_gradient

! INPUTS
!   none
! OUTPUTS
! o partition_tab -- grid integration weight
! o hartree_partition_tab -- grid integration weight for hartree potential
! o hartree_partition_deriv -- partition tab derivative
! o hartree_potential -- overall hartree potential (including any external ions)
! o free_hartree_superpos -- superposition of hartree potential of free atoms; reference
! o free_rho_superpos -- superposition of free atom densities
! o pot_ion_embed -- ion embedding potential
! o rho -- density
! o free_rho_gradient_superpos -- free atom density gradient
! o rho_gradient -- density gradient
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
      real*8, dimension(3) :: embedding_force
      real*8, dimension(3) :: coord_current
      real*8 :: dist_tab_sq(n_centers_basis_integrals)
      real*8 :: dist_tab(n_centers_basis_integrals)
      real*8 :: dir_tab(3,n_centers_basis_integrals)
      real*8 :: dir_tab_norm(3,n_centers_basis_integrals)
      real*8 :: i_r(n_centers_basis_integrals)
      real*8 :: local_rho_gradient(3,n_spin)
!      real*8 :: local_xc_gradient_deriv(3,n_spin)
!      real*8 :: local_xc_density_deriv
!      real*8 :: en_density_xc
      integer :: n_compute_atoms
      integer, dimension(n_centers_basis_integrals) :: center_index
      real*8, dimension(n_centers_basis_integrals) :: temp_free_rho
      integer :: n_compute_occ_atoms
      integer, dimension(n_centers_basis_integrals) :: i_occ2i_compute
      real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho
      real*8, dimension(n_centers_basis_integrals) :: temp_occ_rho_rel
      real*8, dimension(:, :), allocatable :: atom_atom_tab
      real*8, dimension(:),    allocatable :: min_atom_atom_tab
      integer, dimension(:),   allocatable :: atom_atom_index
      real*8 :: aux_spl, rel_rho_spl
      integer :: n_atom_atom_tab
      real*8, dimension(:), allocatable :: atom_atom_dist_list
      integer, dimension(:), allocatable :: atom_idx_A
      integer, dimension(:), allocatable :: atom_idx_B

      logical :: point_on_atom, points_on_multipole, point_on_pseudocore

      character*140 :: info_str

!     counters
      integer :: i_atom, i_species, i_coord, i_point, i_spin,  &
           i_atom_2, i_index, i_multipole
      integer :: current_atom, current_radial, current_angular
      integer :: i_my_batch
      integer :: i_center_L, i_center, i_center_L2

      real*8, dimension(3) :: dir_current
      real*8 :: dist_current, dist_current_sq, dens_current, r_temp
      real*8 ::  dist_to_multipole_sq, dist_to_pseudocore_sq
      integer :: i_pp_atom

      ! kept only for convenience, will not hurt production, may simplify
      ! future debugging efforts:
      integer :: write_out = 0 

!       real*8 :: part_type_at1(n_full_points),part_type_at2(n_full_points)
! test
!      character*20 :: file_name, file_name2
! test end
      weight_tab = 0d0
      hartree_weight_tab = 0d0
      if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) rho_pce = 0.d0
!     initialize partition_tab, free-atom potential and density superposition
      if (grid_storage_initialized) return

      ! CC: Store per atom rho_free on grid
      if (compute_heat_flux) then
        if (.not.allocated(HF_rho_free_per_atom)) allocate(HF_rho_free_per_atom(n_full_points,n_atoms))
        HF_rho_free_per_atom(:,:) = 0.0d0
      end if

      write (info_str,'(2X,A,A)') & 
        "Initializing partition tables, free-atom densities, potentials, etc. ", &
        "across the integration grid (initialize_grid_storage)."
      call localorb_info(info_str,use_unit,'(A)',OL_high)

      n_int_points = 0
      i_point = 0

      ! We begin by preparing the partition table of Stratmann and coworkers. This part
      ! became lengthy over time, please simply take it as one block 
      ! (the end of this part is also marked by a comment below)
      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) &
            .or. (partition_type.eq.9) ) then
        ! All these are variants of the partition table by Stratmann and coworkers.
        ! Incremental improvements during development, keeping them does not hurt. 

        ! Determine what will be the outermost distance at which the stratmann partition
        ! table for a given atom should be non-zero

        ! outer_partition_radius is the outermost radius at which any integrand
        ! or density component will be partitioned to belong to a given atom.
        ! See the extended comment at the top of this subroutines for Equation
        ! numbers in the FHI-aims CPC and other applicable references.

        if (cut_free_atom(1)) then  
        ! there has to be at least one species. Referring to cut_free_atom(1) is an ugly hack, 
        ! but the output is explanatory only anyway.
            write (info_str,'(2X,A)') &
              "| initialize_grid_storage: Actual outermost partition radius vs. multipole_radius_free"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            write (info_str,'(2X,A)') &
              "| (-- VB: in principle, multipole_radius_free should be larger, hence this output)"
            call localorb_info(info_str,use_unit,'(A)',OL_low)
        end if

        do i_species = 1, n_species, 1
          ! If the Stratmann et al. partition table is used to determine the integration weights, then
          ! make sure that the outermost radius where the partition table elements for a given atom
          ! can be non-zero is bounded by the maximum extent of basis functions / densities associated
          ! with that atom.
          !
          ! See the comment at the top of this subroutine for references and equation numbers.
          !
          if (cut_free_atom(i_species)) then
            ! This is the case where the free-atom radius is bounded somehow - should be the normal
            ! case. If the free atom radius is not bounded, we have a problem for periodic systems -
            ! very extended integration regions must be accounted for.
            outer_partition_radius(i_species) = free_r_cut (i_species) + w_cutoff (i_species)
            ! ... and add some diagnostic output here.
            write (info_str,'(2X,A,I8,A,F30.15,A,F30.15,A)') &
              "| Species ", i_species, ": Confinement radius = ", outer_partition_radius(i_species)*bohr, &
              " AA, multipole_radius_free = ", multipole_radius_free(i_species)*bohr, " AA."
            call localorb_info(info_str,use_unit,'(A)',OL_low)
            !
            ! In the case of analytically tabulated Gaussian functions or hydrogenic orbitals,
            ! the extent of basis functions is not guarded by the confinement potential.
            ! In these cases (especially diffuse Gaussian functions) we may need to set a 
            ! larger outer boundary.
            !
            ! multipole_radius_free(i_species) is made sure to be larger than any radial function centered
            ! at this atom in subroutine shrink_fixed_basis_phi_thresh.f90 .
            !
            ! Hence, we employ multipole_radius_free here if it is larger. In the normal case,
            ! where both the free atom and the basis functions are subject to the same cutoff
            ! potential, the difference should be very small.
            !
            ! And yes, I purposely do not employ the 'max' function here. 
            if (multipole_radius_free(i_species).gt.outer_partition_radius(i_species)) then
               outer_partition_radius(i_species) = multipole_radius_free(i_species)
            end if
          else
            ! The free-atom radius is not bounded. Well, we must then use the outermost
            ! radius of the free atom density (very large), for better or for worse.
            outer_partition_radius(i_species) = multipole_radius_free(i_species)
          end if
          write (info_str,'(2X,A,I8,A,F30.15,A)') &
            "| Species ", i_species, ": outer_partition_radius set to ", outer_partition_radius(i_species)*bohr, &
            " AA ."
          call localorb_info(info_str,use_unit,'(A)')

        end do

      end if

      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) ) then
        ! for stratmann partitioning scheme, make available only if necessary
        ! tabulate the interatomic distances for all relevant atoms
        call aims_allocate(atom_atom_tab, n_centers_basis_integrals, n_centers_basis_integrals, "atom_atom_tab")
      else
        ! Not needed, allocate dummy
        call aims_allocate(atom_atom_tab, 1, 1, "atom_atom_tab")
      endif
      
      if (partition_type.eq.9) then
        ! for stratmann partitioning scheme, make available only if necessary
        ! tabulate the interatomic distances for all relevant atoms
        call tab_interatomic_distances_count_entries(n_centers_basis_integrals, centers_basis_integrals, n_atom_atom_tab)

        call aims_allocate(atom_atom_dist_list, n_atom_atom_tab, "atom_atom_dist_list")
        call aims_allocate(atom_idx_A, n_atom_atom_tab, "atom_idx_A")
        call aims_allocate(atom_idx_B, n_centers_basis_integrals+1, "atom_idx_B")
      else
        ! Not needed, allocate dummy
        n_atom_atom_tab = 1

        call aims_allocate(atom_atom_dist_list, 1, "atom_atom_dist_list")
        call aims_allocate(atom_idx_A, 1, "atom_idx_A")
        call aims_allocate(atom_idx_B, 1, "atom_idx_B")
      endif


      call aims_allocate(atom_atom_index, n_centers_basis_integrals, "atom_atom_index")
      call aims_allocate(min_atom_atom_tab, n_centers_basis_integrals, "min_atom_atom_tab")

      if ( (partition_type.eq.5) .or. (partition_type.eq.7) .or. (partition_type.eq.8) ) then
        ! FIXME:
        ! (1) count distance^2 and figure out which atoms could potentially be important in each other's stratmann partition tab
        ! (2) calculate min_distance^2
        ! (3) only allocate as much data as is required, otherwise the arrays may get too big ...
        ! (4) tabulate interatomic distances ONLY for those atoms needed.
        ! (5) remember how many relevant atoms there are for each possible i_atom
        ! the resulting list should be enough to build the stratmann list from it...
        call tab_interatomic_distances(n_centers_basis_integrals, centers_basis_integrals, atom_atom_tab)
        ! for each atom, determine distance to next neighbour:
        do i_atom = 1, n_centers_basis_integrals
          min_atom_atom_tab(i_atom) = 1d100
          do i_atom_2 = 1, n_centers_basis_integrals
            if (i_atom.ne.i_atom_2) then
              min_atom_atom_tab(i_atom) = min(min_atom_atom_tab(i_atom),atom_atom_tab(i_atom,i_atom_2))
            end if
          end do
        end do
        ! need this for comparison ...
        min_atom_atom_tab(:) = (1d0-stratmann_a)*min_atom_atom_tab(:)/2d0

      else if (partition_type.eq.9) then
        call tab_interatomic_distances_local(n_centers_basis_integrals, centers_basis_integrals, &
                                             n_atom_atom_tab, atom_atom_dist_list, atom_idx_A, atom_idx_B)
        ! for each atom, determine distance to next neighbour:
        do i_atom = 1, n_centers_basis_integrals
          min_atom_atom_tab(i_atom) = 1d100
          do i_atom_2 = atom_idx_B(i_atom),atom_idx_B(i_atom+1)-1
            min_atom_atom_tab(i_atom) = min(min_atom_atom_tab(i_atom),atom_atom_dist_list(i_atom_2))
          end do
        end do
        ! need this for comparison ...
        min_atom_atom_tab(:) = (1d0-stratmann_a)*min_atom_atom_tab(:)/2d0
      end if

      ! End all preparations for the Stratmann and coworkers partition table.


! test
!      write(file_name, '(A,I3)') "partition_data.", myid
!      file_name=trim(file_name)
!      open(80,FILE=file_name)
!      write(file_name2, '(A,I3)') "partition_data2.", myid
!      file_name2=trim(file_name2)
!      open(81,FILE=file_name2)
! test end 

      ! for periodic systems, from module free_atoms.f90: average free-atom electrostatic potential in unit cell
      average_free_es_pot = 0.d0

      if (use_partition_deriv) then
      !----------shanghui add for gradient_partition-------------- 
      partition_deriv_delley(1:3,1:n_atoms,1:n_full_points)=0.0d0
      !----------shanghui end add for gradient_partition-------------- 
      endif

      do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

          i_point = i_point + 1


          coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
          current_atom    = batches(i_my_batch) % points(i_index) % index_atom
          current_radial  = batches(i_my_batch) % points(i_index) % index_radial
          current_angular = batches(i_my_batch) % points(i_index) % index_angular

          ! remember where this point is with respect to its own atom as well as its free atom density
          dir_current(:)  = coord_current(:)- &
                            coords_center(:,centers_basis_integrals(current_atom))
          dist_current_sq = dir_current(1)*dir_current(1) &
                          + dir_current(2)*dir_current(2) &
                          + dir_current(3)*dir_current(3)
          dist_current    = sqrt(dist_current_sq)
          dir_current(:)  = dir_current(:)/dist_current
          r_temp          = invert_log_grid &
                            ( dist_current, &
                            r_grid_min(species_center(centers_basis_integrals(current_atom))), &
                            r_grid_inc(species_center(centers_basis_integrals(current_atom))))
          dens_current    =  val_spline &
                ( r_temp, hartree_partition_rho_spl(1,1,species_center(centers_basis_integrals(current_atom))), &
                  n_grid(species_center(centers_basis_integrals(current_atom))) )

!              tabulate current integration point as it appears from spherical
!              coordinates centered at each atom

          call tab_atom_centered_coords_p0 &
              ( coord_current, &
              dist_tab_sq, &
              dir_tab, &
              n_centers_basis_integrals, centers_basis_integrals )

          ! To avoid a floating-point exception, check here whether this integration point happens to sit on
          ! another atom.
          ! VB: We could also check here whether we are inside the innermost logarithmic
          !     grid shell of an atom, as is done further down in evaluate_partition_tab.
          !     This would be even better, since we could then remove the same check from
          !     evaluate_partition altogether ...
          point_on_atom = .false.
          do i_center_L = 1, n_centers_basis_integrals, 1
            if ( dist_tab_sq(i_center_L).eq.0.d0) then
              point_on_atom = .true.
              exit ! exit do loop
            end if
          end do

          point_on_pseudocore = .false.
          do i_pp_atom = 1, n_pp_atoms, 1
            dist_to_pseudocore_sq = &
               (coord_current(1) - pp_coords(1, i_pp_atom))**2 &
              +(coord_current(2) - pp_coords(2, i_pp_atom))**2 &
              +(coord_current(3) - pp_coords(3, i_pp_atom))**2 
            if(dist_to_pseudocore_sq.eq.0.d0) then
              point_on_pseudocore = .true.
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

          if (.not. point_on_atom .and. .not.points_on_multipole .and. .not.point_on_pseudocore) then
            ! This is the normal case. If our grid point does not sit on an atom, we keep it for later use.

            ! For the following operations, we only need those atoms that have a non-zero
            ! free-atom charge density at the current point. We assemble that list here explicitly ...
            n_compute_atoms = 0
            n_compute_occ_atoms = 0

            do i_center_L = 1, n_centers_basis_integrals, 1
              i_center = centers_basis_integrals(i_center_L)

              if ((i_center.eq.current_atom).or. &
                  (dist_tab_sq(i_center_L).lt.multipole_radius_free_sq(species_center(i_center)) )) then
                  ! this center has a non-zero free-atom density, or it belongs to current integration point
                  n_compute_atoms                  = n_compute_atoms + 1
                  center_index(n_compute_atoms)    = i_center
                  dist_tab_sq(n_compute_atoms)     = dist_tab_sq(i_center_L)
                  dir_tab(:,n_compute_atoms)       = dir_tab(:,i_center_L)   ! This ensures consistent handling later, note n_compute_atoms <= i_center_L
                  atom_atom_index(n_compute_atoms) = i_center_L              ! indexing for later use of atom_atom_tab, which is NOT recomputed here for speed reasons
              end if

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
              if (i_center.eq.current_atom) i_center_L = i_center_L2 ! remember the center we are currently at!
              aux_spl = val_spline &
                  ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )
              if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                  rel_rho_spl = val_spline &
                  ( i_r(i_center_L2), free_den_diff_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )
              endif

              temp_free_rho(i_center_L2) = aux_spl
!DB101813
               if ((.not.empty(center_to_atom(i_center))).and.&
                  (.not.(species_pseudoized(species(center_to_atom(i_center)))))) then
!              if (.not.empty(center_to_atom(i_center))) then
                  n_compute_occ_atoms = n_compute_occ_atoms + 1
                  i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
                  temp_occ_rho(n_compute_occ_atoms) = aux_spl
                  if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)&
                    temp_occ_rho_rel(n_compute_occ_atoms) = rel_rho_spl
              end if

            end do


!               evaluate Hartree partition tab first
            call evaluate_partition_p2 &
            ( flag_hartree_partition_type,      &
              current_atom,                     &
              i_center_L,                       &
              dist_current,                     &
              dir_current,                      &
              dist_current_sq,                  &
              dens_current,                     &
              dist_tab,                         &
              i_r,                              &
              dir_tab_norm,                     &
              w_angular(current_angular, current_radial, species(current_atom)),  &
              hartree_partition_tab(i_point),   &
              hartree_weight_tab(i_point),      &
              hartree_partition_deriv(i_point), &
              n_centers_basis_integrals,        &
              n_compute_atoms,                  &
              center_index,                     &
              temp_free_rho,                    &
              dist_tab_sq,                      &
              atom_atom_tab,                    &
              atom_atom_index,                  &
              min_atom_atom_tab,                &
              n_atom_atom_tab,                  &
              atom_atom_dist_list,                 &
              atom_idx_A,                       &
              atom_idx_B,                       &
              write_out, & 
              partition_deriv_delley(1:3,1:n_atoms,i_point) )  
!           evaluate partition_tab
            if (partition_type.ne.flag_hartree_partition_type) then
              ! Partition functions for integration and Hartree potential are different

              ! re-evaluate pieces for the partition tab to be on the safe side
              if ( (partition_type .ne. 5) .and. (partition_type .ne. 7) &
                   .and. (partition_type .ne. 8) .and. (partition_type .ne. 9) ) then ! no density needed for stratmann
                dens_current    =  val_spline &
                  ( r_temp, partition_rho_spl(1,1,species_center(centers_basis_integrals(current_atom))), &
                    n_grid(species_center(centers_basis_integrals(current_atom))) )
                do i_center_L2 = 1, n_compute_atoms
                  i_center = center_index(i_center_L2)
                  if (i_center.eq.current_atom) i_center_L = i_center_L2 ! remember the center we are currently at!
                  aux_spl = val_spline &
                  ( i_r(i_center_L2), hartree_partition_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )

                  temp_free_rho(i_center_L2) = aux_spl

                end do

              end if


  
              call evaluate_partition_tab_p2  &
                ( current_atom,              &
                  i_center_L,                &
                  dist_current,              &
                  dist_current_sq,           &
                  dens_current,              &
                  dist_tab,                  &
                  i_r,                       &
                  w_radial(current_radial, species(current_atom)),  &
                  w_angular(current_angular, current_radial, species(current_atom)),  &
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


            else
              ! Partition function for integration is the same as for Hartree potential
              ! and just needs some integration weights ...


              weight_tab(i_point) = hartree_weight_tab(i_point) * &
                w_radial(current_radial, species(current_atom)) * &
                dist_current_sq

              partition_tab(i_point) = hartree_partition_tab(i_point) * &
                w_radial(current_radial, species(current_atom)) * &
                dist_current_sq


!            write (81,'(2X,F15.8,2X,F15.8,2X,F15.8,2X,F15.8,2X, I3,A)')  coord_current(1), coord_current(2), coord_current(3), &
!              partition_tab(i_point),current_atom, " coord, partition_tab"

              if (use_partition_deriv) then
              !----------shanghui add for gradient_partition-------------- 
              partition_deriv_delley(1:3,1:n_atoms,i_point) = &  
              partition_deriv_delley(1:3,1:n_atoms,i_point) * & 
                 w_radial(current_radial, species(current_atom)) * &
                 dist_current_sq
              !----------shanghui end add for gradient_partition--------------
              endif
            end if

            ! Fallout from the partition_tab fix rho_r2_lda: More complexity
            if ( flag_hartree_partition_type.eq.6 ) then
              ! must retabulate densities of non-empty atoms for the initial electron density here;
              ! the previous tabulation was LDA, possibly not the correct XC.
              n_compute_occ_atoms = 0
              do i_center_L2 = 1, n_compute_atoms
                i_center = center_index(i_center_L2)
! DB 092812:  we should also exclude pseudocores at this place.
!             won't change the result, but should give a minimal speed-up.
               if ((.not.empty(center_to_atom(i_center))).and.&
                  (.not.(species_pseudoized(species(center_to_atom(i_center)))))) then

!                if (.not.empty(center_to_atom(i_center))) then

                  aux_spl = val_spline &
                  ( i_r(i_center_L2), free_rho_spl(1,1,species_center(i_center)), &
                    n_grid(species_center(i_center)) )
                  if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                    rel_rho_spl = val_spline &
                    ( i_r(i_center_L2), free_den_diff_spl(1,1,species_center(i_center)), &
                      n_grid(species_center(i_center)) )
                  endif

                  n_compute_occ_atoms = n_compute_occ_atoms + 1
                  i_occ2i_compute(n_compute_occ_atoms) = i_center_L2
                  temp_occ_rho(n_compute_occ_atoms) = aux_spl
                  if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)&
                    temp_occ_rho_rel(n_compute_occ_atoms) = rel_rho_spl
                end if

              end do
            end if

            !
            ! FIXME: The sum over free&relevant atoms might not necessarily be the same as the sum
            !        over free and relevant ions, especially if they are negatively charged!
            !        This might be inconsistent in the current code version.
            !


            call evaluate_free_atom_sums_p2  &
                ( i_r, dir_tab_norm,  &
                free_hartree_superpos(i_point),  &
                free_rho_superpos(i_point),  &
                free_rho_gradient_superpos(1, i_point), &
                rho(1,i_point),  &
                local_rho_gradient, n_compute_atoms, center_index, &
                n_compute_occ_atoms, i_occ2i_compute, temp_occ_rho, i_point )
                !SR: for evaluation of laplacian of nfree:
                !,dist_tab_sq,&
                !free_rho_laplace_superpos(i_point))
            if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
              call evaluate_free_atom_sums_rel (i_r, dir_tab_norm, rho_pce(1,i_point), rho_pce_gradient(1,1,i_point), &
                n_compute_atoms, center_index, n_compute_occ_atoms, i_occ2i_compute, temp_occ_rho_rel)
            end if


            if (n_periodic.gt.0) then
                average_free_es_pot = average_free_es_pot + &
                free_hartree_superpos(i_point) * partition_tab(i_point)
            end if

            if (use_density_gradient) then
              do i_spin = 1, n_spin, 1
                do i_coord = 1,3,1
                  rho_gradient(i_coord,i_spin,i_point) =   &
                        local_rho_gradient(i_coord,i_spin)
                enddo
              enddo
            end if

            if (use_embedding_potential) then
              call embedding_potential  &
                    ( coord_current,  &
                    pot_ion_embed(i_point), embedding_force )

            end if

            hartree_potential(i_point) =   &
              free_hartree_superpos(i_point)

            if (use_embedding_potential) then
              if (full_embedding) then
                  hartree_potential(i_point) =   &
                  hartree_potential(i_point) +  &
                  pot_ion_embed(i_point)
              end if
            end if

!DB - when pseudopot are used, we have to add the local potential initially 

            if (use_embedding_pp) then
                  hartree_potential(i_point) =   &
                  hartree_potential(i_point) +  &
                  whole_local_pseudpot_on_intgrid(i_point)
            end if


            if (partition_tab(i_point).gt.0) then
              n_int_points = n_int_points + 1
            end if


          else
          ! this is the case where the current integration point sits exactly on top of another atom
          ! must never be considered in later dealings!

            hartree_partition_tab(i_point) = 0.d0
            if (multipole_interpolation_style.eq.1) then
              ! legacy only, no longer used
              hartree_partition_deriv(i_point) = 0.d0
            end if
            partition_tab(i_point) = 0.d0

            free_hartree_superpos(i_point) = 0.d0
            free_rho_superpos(i_point) = 0.d0
            if ((multipole_interpolation_style.eq.1) .or. (use_density_gradient)) then
              free_rho_gradient_superpos(:,i_point) = 0.d0
            end if
            rho(:,i_point) = 0.d0
            if (use_density_gradient) then
              rho_gradient(:,:,i_point) = 0.d0
            end if
            if (use_embedding_potential) then
              pot_ion_embed(i_point) = 0.d0
            end if
            if (use_embedding_pp) then
              whole_local_pseudpot_on_intgrid(i_point) = 0.d0
              if(use_nonlinear_core) then
                 partial_core_rho(i_point) = 0.d0
                 partial_core_rho_grad(:,i_point) = 0.d0
              endif
            end if
            hartree_potential(i_point) = 0.d0
          end if
! 	  if (current_atom ==1) then
! 	    part_type_at1(i_point) = hartree_partition_tab(i_point)
! 	  else
! 	    part_type_at1(i_point) = 100
! 	  end if
! 	  if (current_atom ==2) then
! 	    part_type_at2(i_point) = hartree_partition_tab(i_point)
! 	  else
! 	    part_type_at2(i_point) = 100
! 	  end if
        end do ! loop over batch
      end do ! loop over batches

      call sync_workload( n_full_points, n_int_points,   &
           n_full_points_total, n_int_points_total )

      write (info_str,'(2X,A,I8)')  &
              "| Net number of integration points: ",   &
              n_full_points_total
      call localorb_info(info_str,use_unit,'(A)', OL_high)
      write (info_str,'(2X,A,I8)')  &
             "| of which are non-zero points    : ", n_int_points_total
      call localorb_info(info_str,use_unit,'(A)', OL_high)

      if (n_periodic.gt.0) then
         call sync_average_potential ( average_free_es_pot )
         average_free_es_pot = average_free_es_pot/cell_volume
         write (info_str,'(2X,A,1X,F15.8,A)')  &
            "| Numerical average free-atom electrostatic potential    :",   &
            average_free_es_pot*hartree, " eV"
         call localorb_info(info_str,use_unit,'(A)', OL_norm)

         if (analytic_potential_average) then
            ! In this case, we switch to the analytical average which we can
            ! determine from the free atoms
            average_free_es_pot = 0.d0
            do i_atom = 1, n_atoms, 1
               average_free_es_pot = average_free_es_pot + free_atom_average_es_pot(species(i_atom))
            enddo
            average_free_es_pot = average_free_es_pot/cell_volume
            write (info_str,'(2X,A,1X,F15.8,A)')  &
               "| Analytical average free-atom electrostatic potential:   ",   &
               average_free_es_pot*hartree, " eV"
            call localorb_info(info_str,use_unit,'(A)', OL_norm)
         end if
      end if

      grid_storage_initialized = .true.

      ! clean up the odd allocated variable
      call aims_deallocate(min_atom_atom_tab, "min_atom_atom_tab")
      call aims_deallocate(atom_atom_index, "atom_atom_index")
      call aims_deallocate(atom_atom_tab, "atom_atom_tab")
      call aims_deallocate(atom_atom_dist_list, "atom_atom_dist_list")
      call aims_deallocate(atom_idx_A, "atom_idx_A")
      call aims_deallocate(atom_idx_B, "atom_idx_B")


if (use_embedding_pp) then
!! from here on we need to set some free atom properties to zero if (species_pseudoized)
  do i_species =1, n_species

    if (species_pseudoized(i_species).or.no_basis(i_species)) then
        free_rho(:,i_species) = 0.d0
        free_rho_spl(:,:,i_species) = 0.d0
        if (allocated(free_drho_dr_spl)) then
           free_drho_dr_spl(:,:,i_species) = 0.d0
        end if
    end if
  end do
end if


! test
!      close(80)
!      close(81)
! test end

    end subroutine initialize_grid_storage_p1
!******
end module initialize_grid_storage
