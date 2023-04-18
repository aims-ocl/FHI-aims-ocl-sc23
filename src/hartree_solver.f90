!****s* FHI-aims/hartree_solver
!  NAME
!    hartree_solver
!  SYNOPSIS

subroutine hartree_solver &
  ( partition_tab, rho, delta_v_hartree )

  !  PURPOSE
  !  Calculates hartree potential given input (delta) rho at integration points using 
  !  multipole expansion.
  !  
  !  Most of the code here is repurposed from update_hartree_potential_p1 and
  !  sum_up_whole_potential_p1
  !
  !  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use analytical_stress
  use hartree_potential_storage, only : rho_multipole_index
  use hartree_potential_real_p0
  use hartree_potential_recip
  use hartree_non_periodic_ewald
  implicit none

  !  ARGUMENTS

  real*8, dimension(n_full_points), intent(in)                 :: partition_tab  
  real*8, dimension(n_full_points), intent(in)                 :: rho
  real*8, dimension(n_full_points), intent(out)                :: delta_v_hartree

  !  INPUTS
  !   o  partition_tab -- values of partition function
  !   o  rho -- electron density
  !
  !  OUTPUT
  !   o  delta_v_hartree --  Hartree potential at integration points
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
  !    Created in November 2017
  !  SOURCE


  !  rho_multipole is the angular integral over
  !  Y_lm(Omega) * partition_fn(at,r,Omega) * delta_rho(r,Omega),
  !  Delley 1990, Eq. (11) 
  !  i.e. rho_multipole = rho_multipole(r) . For convenience

  real*8, dimension(n_atoms) :: multipole_radius_sq
  real*8, dimension((l_pot_max + 1) ** 2, n_atoms) :: multipole_moments
  integer, dimension(n_atoms) :: l_hartree_max_far_distance
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius
  real*8, dimension(n_atoms) :: delta_v_hartree_part_at_zero
  real*8, dimension(3, n_atoms) :: delta_v_hartree_deriv_l0_at_zero


  real*8 dir_tab(3)
  real*8 temp_rho_new
  real*8 multipole_radius
  real*8 ylm_tab((l_pot_max+1)**2)
  real*8, allocatable :: current_rho_multipole_spl(:,:,:), current_rho_multipole(:,:)

  integer l_h_dim(n_atoms)
  real*8 :: drho

  !  counters

  integer i_atom
  integer i_index
  integer current_atom, current_radial, current_angular
  integer i_full_points
  integer i_my_batch
  integer i_atom_index

  integer i_l, i_m

  ! Parameters
  integer, parameter :: n_coeff_hartree = 2 

  ! Allocatable
  real*8, allocatable :: rho_multipole(:, :, :)
  real*8, allocatable :: delta_v_multipole_spl(:, :, :, :)

  ! Reconstruction
  real*8 :: dist_tab_sq, ddot
  real*8 :: coord_current(3), trigonom_tab(4)
  integer :: i_center, current_spl_atom, current_center

  real*8, allocatable :: aux_delta_v_hartree_result(:)
  integer, allocatable :: index_lm(:, :)
  real*8, allocatable :: max_lm_spl_rad_sq(:)
  real*8, allocatable :: current_delta_v_hart_part_spl(:, :, :)

  integer :: atom_of_splines
  real*8 :: delta_v_hartree_aux
  integer :: l_h_dim_2
  integer :: i_atom_2
  real*8 :: dist_tab_in
  real*8 :: i_r
  real*8 :: i_r_log
  real*8 :: dir_tab_in(3)
  real*8, dimension((l_pot_max + 1) ** 2) :: rho_multipole_component
  integer :: l_atom_max
  integer :: i_lm
  integer :: i_batch
  logical :: forces_on = .false. ! Work-around
  real*8 :: dist_tab_out
  real*8 :: dip_gradient, dip_origin 
  real*8 :: dip_length, dip_coord_current
  real*8, dimension((l_pot_max + 1) ** 2) :: delta_v_hartree_multipole_component
  real*8, dimension(:), allocatable :: adap_outer_radius_sq

  !  begin work

  !  Allocations & initializations:
  allocate(current_rho_multipole((l_pot_max+1)**2, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole')

  allocate(current_rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=i_index)
  call check_allocation(i_index, 'current_rho_multipole_spl')

  allocate(rho_multipole((l_pot_max+1)**2, n_max_radial+2, n_atoms))
  rho_multipole = 0

  allocate(delta_v_multipole_spl((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, n_atoms))
  delta_v_multipole_spl = 0

  allocate(current_delta_v_hart_part_spl((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid))

  allocate(aux_delta_v_hartree_result((l_pot_max+1)**2))
  allocate(max_lm_spl_rad_sq((l_pot_max+1)**2))

  allocate(adap_outer_radius_sq(n_atoms))

  allocate(index_lm(-l_pot_max:l_pot_max, 0:l_pot_max))
  index_lm = 0
  i_index = 0
  do i_l = 0, l_pot_max, 1
    do i_m = -i_l, i_l
      i_index = i_index + 1
      index_lm(i_m, i_l) = i_index
    end do
  end do


  do i_atom = 1, n_atoms, 1
    l_h_dim(i_atom) = (l_hartree(species(i_atom))+1)**2
  end do

  rho_multipole = 0.0d0
  i_full_points = 0

  do i_my_batch = 1, n_my_batches, 1
    do i_index = 1, batches(i_my_batch)%size, 1
           
      current_atom = batches(i_my_batch) % points(i_index) % index_atom
      current_radial = batches(i_my_batch) % points(i_index) % index_radial
      current_angular = batches(i_my_batch) % points(i_index) % index_angular

      ! run for the first time through the integration grid to tabulate 
      ! integral_1_zero_r and thus calculate integral_zero_infinity
      !
      ! meanwhile tabulate first greens-function solution of partitioned potential
      ! and difference charge density
      
      i_full_points = i_full_points + 1

      ! execute only if partition_tab.gt.0 here, i.e. if the integration point
      ! makes sense
      !!!!! This is NOT the usual integration partition tab
      !!!!! This is the HARTREE partition tab. 
      !!!!! Thus this criterion is not the same as in all other integrations.
      !!!!! We must fix that one day.
      if (partition_tab(i_full_points).gt.0.d0) then

        ! Check if the atom has been correctly detected as mine
        i_atom_index = rho_multipole_index(current_atom)
        if(i_atom_index<=0) then
          print '(2(a,i5))','ID ',myid,' INTERNAL ERROR hartree_solver - need atom! Atom: ',current_atom
          call aims_stop
        endif

        ! compute atom-centered coordinates of current integration point,
        ! BUT HERE ONLY FOR PRESENT ATOM!
        dir_tab(:) = r_angular( : , current_angular, current_radial, species(current_atom))

        ylm_tab (1:l_h_dim(current_atom)) = &
             local_ylm_tab(1:l_h_dim(current_atom),current_angular, &
             lebedev_grid_index(current_radial,species(current_atom)))

        ! calculate contribution to the angular parts of the two integrals which
        ! are equal 

        ! notice that the radial and angular integration weights are already part of 
        ! partition_tab; we need not multiply them in separately.

        ! consider only difference density for hartree potential

        !!!temp_rho_new = 0.d0
        temp_rho_new = rho(i_full_points)

        ! implied loop over all (l,m) from (0,0) to (l_hartree,l_hartree)
        rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) = &
             rho_multipole(1:l_h_dim(current_atom), current_radial+1, i_atom_index) + &
             ylm_tab(1:l_h_dim(current_atom)) * partition_tab(i_full_points) * temp_rho_new 

      end if
    end do !    end loop over a batch
  end do !     end loop over batches

  ! Add all rho_multipole contributions

  do i_atom = 1, n_atoms

    i_atom_index = rho_multipole_index(i_atom)
    if(i_atom_index > 0) then
      current_rho_multipole(:,:) = rho_multipole(:,:,i_atom_index)
    else
      current_rho_multipole(:,:) = 0.
    endif

    call sync_vector(current_rho_multipole, ((l_pot_max+1)**2)*(n_max_radial+2))

    if(i_atom_index > 0) then
      rho_multipole(:,:,i_atom_index) = current_rho_multipole(:,:)
    endif

  enddo


  ! Set boundaries on rho_multipole

  do i_atom_index = 1, n_atoms

    i_atom = i_atom_index

    ! At Infinity (n_radial+2)
    ! enforce zero charge at infinity explicitly by a trick:
    ! (used in splines later on)
    rho_multipole(1:l_h_dim(i_atom), n_radial(species(i_atom))+2, i_atom_index) = 0.d0

    ! At zero
    drho = (rho_multipole(1,2,i_atom_index) - rho_multipole(1,3,i_atom_index)) &
    &     /    (r_radial(1,species(i_atom)) - r_radial(2,species(i_atom)))

    if (legacy_monopole_extrapolation) then
       ! Backwards compatible but inaccurate choice of sign.
       ! Please note that the influence of this boundary rapidly decreases
       ! with increasing grid density.
       rho_multipole(1,1,i_atom_index) = rho_multipole(1,2,i_atom_index) &
       &                                 + drho * r_radial(1,species(i_atom))
    else
       rho_multipole(1,1,i_atom_index) = rho_multipole(1,2,i_atom_index) &
       &                                 - drho * r_radial(1,species(i_atom))
    end if
    rho_multipole(2:l_h_dim(i_atom),1,i_atom_index) = 0.d0

  enddo

  ! Calculate output variables
  ! the variables must be set to 0 since they are sync'd at the end
  delta_v_hartree_part_at_zero(:) = 0
  delta_v_hartree_deriv_l0_at_zero(:, :) = 0
  multipole_moments(:, :) = 0
  multipole_radius_sq(:) = 0
  l_hartree_max_far_distance(:) = 0
  outer_potential_radius(:, :) = 0

  do i_atom = 1, n_atoms

    if (mod(i_atom-1,n_tasks) == myid) then

      call get_rho_multipole_spl_dev(current_rho_multipole_spl, i_atom, rho_multipole)
    
    if (force_hartree_log_grid) then
      call integrate_hartree_log_grid( &
        i_atom, current_rho_multipole_spl, &
        delta_v_hartree_part_at_zero(i_atom), &
        delta_v_hartree_deriv_l0_at_zero(1:3, i_atom), &
        multipole_moments(1, i_atom), &
        multipole_radius, &
        l_hartree_max_far_distance(i_atom), &
        outer_potential_radius(0:l_pot_max, i_atom) )

      multipole_radius_sq(i_atom) = multipole_radius ** 2
    else
      write(use_unit,*) 'Only log integrals are supported'
      stop
    end if ! force_hartree_log_grid

      ! Integrate splined rho multipoles to get delta v splines
      !call integrate_delta_v_hartree(current_rho_multipole_spl, delta_v_multipole_spl(:, :, :, i_atom), n_coeff_hartree, i_atom)
      !     end distribution over threads
    end if

    !     end loop over atoms 
  end do

  ! synchronize results

  call sync_vector(delta_v_hartree_part_at_zero, n_atoms)
  call sync_vector(delta_v_hartree_deriv_l0_at_zero, 3 * n_atoms)
  call sync_vector(multipole_moments, (l_pot_max + 1) ** 2 * n_atoms)
  call sync_vector(multipole_radius_sq, n_atoms)
  call sync_integer_vector(l_hartree_max_far_distance, n_atoms)
  call sync_vector(outer_potential_radius, (l_pot_max + 1) * n_atoms)

  ! The rest of the code is repurposed from sum_up_whole_potential_p1

  call hartree_potential_real_coeff( &
    index_lm, multipole_moments, &
    l_hartree_max_far_distance, n_centers_hartree_potential )

  if (n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
    call update_outer_radius_l( &
      outer_potential_radius, multipole_moments, &
      multipole_radius_sq, l_hartree_max_far_distance, &
      index_lm )
    if ( .not. use_hartree_non_periodic_ewald ) then
      do i_atom_2 = 1, n_atoms
        adap_outer_radius_sq(i_atom_2) = maxval(outer_potential_radius(:, i_atom_2))
      end do
    else
      adap_outer_radius_sq(:) = multipole_radius_sq(:)
    end if ! .not. use_hartree_non_periodic_ewald
  else
    adap_outer_radius_sq = 1e8
  end if !n_periodic > 0 .or. use_hartree_non_periodic_ewald

  if (n_periodic > 0) then
    call evaluate_hartree_recip_coef(index_lm, multipole_moments, l_hartree_max_far_distance, .false.)
  else if (use_hartree_non_periodic_ewald) then
    call calculate_extended_hartree_non_periodic_ewald( &
      l_hartree_max_far_distance, multipole_radius_sq, &
      adap_outer_radius_sq )
  end if

  ! intialize variables: 
  current_spl_atom = 0
  delta_v_hartree  = 0d0

  dip_gradient = 0
  dip_origin = 0
  dip_length = 0

  ! for any given atom, sum over all centers that relate to this atom 
  atom_of_splines = 0
  do i_center = 1, n_centers_hartree_potential, 1

    current_center = centers_hartree_potential(i_center)
    current_spl_atom = center_to_atom(current_center)
    
    if ((n_periodic > 0) .or. force_new_functional) then
      
      ! In this case, use multipole_spline to compute Hartree potential components on each individual
      ! atomic nucleus, before doing anything else.

      if (mod(current_spl_atom-1, n_tasks) == myid) then
        ! For this case rho_multipole for current_spl atom is on myid

        do i_atom_2 = 1, n_atoms, 1
          if (empty(i_atom_2)) cycle

          ! Reinitialize the physical quantities of interest
          delta_v_hartree_aux = 0.d0
          coord_current(:) = coords(:, i_atom_2)
          
          ! Tabulate distances and directions to all atoms -
          ! including relative positions on logarithmic and
          ! radial integration grids.

          call tab_single_atom_centered_coords_p0( &
            current_center, coord_current, &
            dist_tab_sq, dir_tab)

          if (i_atom_2 .eq. current_center) then
            dist_tab_sq = r_grid_min(species(i_atom_2)) ** 2 + 1e-15
          end if

          ! At each integration point, the Hartree potential components coming from
          ! different atoms are split into two groups:
          ! Any components coming from atoms close by are evaluated by explicit
          ! numerical splines; far away atoms are simply represented by
          ! an analytical long-distance multipole potential.
          if (dist_tab_sq .lt. multipole_radius_sq(current_spl_atom)) then
            ! Begin with everything related to n_atoms_in
            
            if (current_spl_atom /= atom_of_splines) then
              call get_rho_multipole_spl_dev(current_rho_multipole_spl, current_spl_atom, rho_multipole)
              
              call integrate_delta_v_hartree( &
                current_rho_multipole_spl, current_delta_v_hart_part_spl, &
                n_coeff_hartree, current_spl_atom)
              
              atom_of_splines = current_spl_atom

            end if ! current_spl_atom /= atom_of_splines

            call tab_single_atom_centered_coords_radial_log_p0( &
              current_center, dist_tab_sq, &
              dir_tab, dist_tab_in, &
              i_r, i_r_log, &
              dir_tab_in )

            ! For all inner atoms, we need ylm functions and their gradients explicitly
            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)
            
            call tab_single_wave_ylm_p0( &
              current_center, trigonom_tab, &
              l_hartree, l_pot_max, &
              ylm_tab)

            ! Now loop over all inner atoms (those which are close enough
            ! to the current integration point so we need explicit numerical splines)
            !   partitioned
            !   hartree potentials need to be summed up
            !   according to Delley (eq. 12c)

            l_h_dim_2 = (l_hartree(species_center(current_center)) + 1) ** 2
            ! obtain spline-interpolated values of the multipole components
            ! of the partitioned Hartree potential, splined on the logarithmic
            ! integration grid

            call spline_vector_v2( &
              i_r_log, current_delta_v_hart_part_spl, &
              (l_pot_max + 1) ** 2, n_coeff_hartree, &
              n_hartree_grid, n_grid(species_center(current_center)), &
              l_h_dim_2, delta_v_hartree_multipole_component ) 

            if (i_atom_2 .eq. current_center) then
              delta_v_hartree_multipole_component(1) = delta_v_hartree_part_at_zero(current_spl_atom)

              do i_lm = 2, l_h_dim_2
                delta_v_hartree_multipole_component(i_lm) = &
                  2 * current_delta_v_hart_part_spl(i_lm, 1, 1) - &
                  current_delta_v_hart_part_spl(i_lm, 1, 2)
              end do ! i_lm
            end if ! i_atom_2 = current_center

            ! Sum up the Hartree potential contribution from the present atom i_atom_in
            delta_v_hartree_aux = delta_v_hartree_aux + &
              ddot( l_h_dim_2, &
              delta_v_hartree_multipole_component, 1, ylm_tab, 1 )
            
            if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
              call far_distance_hartree_Fp_periodic_single_atom &
                (current_spl_atom, i_center, &
                dist_tab_in, l_hartree_max_far_distance, .true., forces_on, &
                multipole_radius_sq(current_spl_atom),                      &
                sqrt( adap_outer_radius_sq(current_spl_atom) )                 )
            end if

          else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom) ) then

            ! Tabulate distances only for outer part
            dist_tab_out = sqrt(dist_tab_sq)

            if (n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald) then

              ! Recursively tabulate the radial behavior of all multipole components
              ! These radial functions are a private variable in module
              ! hartree_potential_real.f90, and are reused there later in
              ! subroutine far_distance_hartree_potential_real_single_atom
              call far_distance_hartree_Fp_cluster_single_atom_p2 &
                   ( dist_tab_out, &
                   l_hartree_max_far_distance( current_spl_atom), forces_on )

              call far_distance_real_hartree_potential_single_atom_p2 &
                   ( current_spl_atom, delta_v_hartree_aux, &
                   l_hartree_max_far_distance(current_spl_atom), coord_current )

            
            else
              ! Now sum up the potential contributions from all far-field atoms
              ! (analytical multipole potentials only ..)

              call far_distance_hartree_Fp_periodic_single_atom( &
                current_spl_atom, i_center, &
                dist_tab_out, l_hartree_max_far_distance, &
                .false., forces_on, &
                multipole_radius_sq(current_spl_atom), sqrt(adap_outer_radius_sq(current_spl_atom)) )
            end if

          end if ! Multipole radius

          if (dist_tab_sq .lt. max(adap_outer_radius_sq(current_spl_atom), &
            multipole_radius_sq(current_spl_atom)) ) then

            if (n_periodic > 0 .or. use_hartree_non_periodic_ewald) then
              
              ! Now sum up the far distance parts of the Hartree potential
              call far_distance_real_hartree_potential_single_atom( &
                current_center, i_center, &
                delta_v_hartree_aux, l_hartree_max_far_distance, &
                coord_current )

            end if ! n_periodic
          end if ! dist_tab_sq ...

        end do ! i_atom_2
      end if ! myid
    end if ! n_periodic
    !
    ! Now follows the real work: Summing the multipole potential
    ! on the integration grid
    !

    ! Reset grid counter for current Hartree potential center
    i_full_points = 0

    do i_batch = 1, n_my_batches
      
      ! loop over one batch
      do i_index = 1, batches(i_batch)%size, 1

        ! i_full_points is the index that indicates where we are in the entire grid 
        ! (for external quantities like rho, potential, ...)
        i_full_points = i_full_points + 1

        ! Only execute if partition_tab is > 0, else we can run into 1/0 !!!
        if (partition_tab(i_full_points) .gt. 0.d0) then
          
          ! Get current integration point coordinate
          coord_current(:) = batches(i_batch) % points(i_index) % coords(:)

          ! Tabulate distances and directions to a single atom -
          ! including relative positions on logarithmic and radial
          ! integration grids.

          call tab_single_atom_centered_coords_p0( &
            current_center, coord_current, &
            dist_tab_sq, dir_tab )

          l_atom_max = l_hartree(species(current_spl_atom))
          do while ( (outer_potential_radius(l_atom_max, current_spl_atom) .lt. dist_tab_sq ) &
            .and. (l_atom_max .gt. 0 ) )

            l_atom_max = l_atom_max - 1
          end do ! while

          ! At each integration point, the Hartree potential components coming from
          ! different atoms are split into two groups:
          ! Any components coming from atoms close by are evaluated by explicit
          ! numerical splines; far away atoms are simply represented by
          ! an analytical long-distance multipole potential.

          ! We treat first the atoms close by
          if (dist_tab_sq .lt. multipole_radius_sq(current_spl_atom) ) then

            if (current_spl_atom /= atom_of_splines) then
            
              call get_rho_multipole_spl_dev(current_rho_multipole_spl, current_spl_atom, rho_multipole)
              call integrate_delta_v_hartree( &
                current_rho_multipole_spl, current_delta_v_hart_part_spl, &
                n_coeff_hartree, current_spl_atom )
              atom_of_splines = current_spl_atom
            end if ! current_spl_atom /= atom_of_splines

            call tab_single_atom_centered_coords_radial_log_p0( &
              current_center, dist_tab_sq, &
              dir_tab, dist_tab_in, &
              i_r, i_r_log, &
              dir_tab_in )

            ! For an inner atom we need ylm functions and their gradients explicitly
            
            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

            call tab_single_wave_ylm_p2( &
              trigonom_tab, l_atom_max, &
              l_pot_max, ylm_tab)

            
            ! For an inner atoms (those which are close enough
            ! to the current integration point so we need explicit numerical splines)
            !     partitioned
            !     hartree potentials need to be summed up
            !     according to Delley (eq. 12c)

            l_h_dim_2 = (l_atom_max + 1) ** 2

            ! Obtain spline-interpolated values of the multipole components
            ! of the partitioned Hartree potential, splined on the logarithmic integration grid
            call spline_vector_v2( &
              i_r_log, current_delta_v_hart_part_spl, &
              (l_pot_max + 1) ** 2, n_coeff_hartree, &
              n_hartree_grid, n_grid(species_center(current_center)), &
              l_h_dim_2, delta_v_hartree_multipole_component )

            ! Sum up the Hartree potential contribution from present inner atom
            delta_v_hartree_aux = ddot( l_h_dim_2, delta_v_hartree_multipole_component, 1, ylm_tab, 1)

            delta_v_hartree(i_full_points) = delta_v_hartree(i_full_points) + delta_v_hartree_aux

            ! Obtain spline-interpolated values of the multipole density,
            ! this time splined on the radial integration grid

            call spline_vector_v2( &
              i_r + 1, current_rho_multipole_spl, &
              (l_pot_max + 1) ** 2, n_max_spline, &
              n_max_radial + 2, n_radial(species_center(current_center)) + 2, &
              l_h_dim_2, rho_multipole_component)

            ! Far-distance treatment for the center i_center for the periodic system
            if (n_periodic > 0 .or. use_hartree_non_periodic_ewald) then
              call far_distance_hartree_Fp_periodic_single_atom( &
                current_spl_atom, i_center, &
                dist_tab_in, l_hartree_max_far_distance, .true., &
                forces_on, multipole_radius_sq(current_spl_atom), &
                sqrt(adap_outer_radius_sq(current_spl_atom)) )
            end if

          else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then
            
            ! The current center is in the far-distance part
            ! Tabulate distance for the outer part
            dist_tab_out = sqrt(dist_tab_sq)

            ! Now sum up the potential contributions from a far-field atom
            ! (analytical multipole potentials only ...)
            if (n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald) then

              ! Recursively tabulate the radial behaviour of all multipole components
              ! These radial functions are a private variable in module
              ! hartree_potential_real.f90, and are reused there later in 
              ! subroutine far_distance_hartree_potential_real_single_atom
              call far_distance_hartree_Fp_cluster_single_atom_p2( &
                dist_tab_out, l_atom_max, &
                forces_on )

              call far_distance_real_hartree_potential_single_atom_p2( &
                i_center, delta_v_hartree(i_full_points), &
                l_atom_max, coord_current )

            else ! Periodic system or non-periodic Ewald method

              call far_distance_hartree_Fp_periodic_single_atom( &
                current_spl_atom, i_center, &
                dist_tab_out, l_hartree_max_far_distance, &
                .false., forces_on, &
                multipole_radius_sq(current_spl_atom), sqrt(adap_outer_radius_sq(current_spl_atom)) )

                ! Note: 'far_distance_real_hartree_potential_single_atom' is called in the periodic
                ! systems later, because then we have Fp in inner and outer multipole radius

            end if ! (n_periodic == 0 .and. use_hartree_non_periodic_ewald)

          end if ! dist_tab_sq

          ! Now sum up the far distance part of the Hartree potential if Ewald's
          ! decomposition is performed. In the opposite case, this was already done

          if (n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
          
            if (dist_tab_sq .lt. &
              max(adap_outer_radius_sq(current_spl_atom), multipole_radius_sq(current_spl_atom))) then
              
              ! Far distance analytic real-space part of the potential

              call far_distance_real_hartree_potential_single_atom( &
                current_center, i_center, &
                delta_v_hartree(i_full_points), l_hartree_max_far_distance, &
                coord_current )

            end if ! dist_tab_sq >

          end if ! n_periodic > 0 .or. use_hartree_non_periodic_ewald 

        end if ! partition_tab(i_full_points) .gt. 0.d0

      end do ! i_index

    end do ! i_batch

  end do          ! loop over relevant Hartree centers

  ! Next, we integrate the actual energy 

  i_full_points = 0
  do i_batch = 1, n_my_batches
    
    ! Loop over one batch
    do i_index = 1, batches(i_batch)%size, 1

      ! i_full_points is the index that indicates where we are in the entire grid (for external 
      ! quantities like rho, potential, ...)
      i_full_points = i_full_points + 1

      if (partition_tab(i_full_points) .gt. 0.d0) then
        
        ! The reciprocal space contribution needs to be treated separately
        ! from the centers above

        if (n_periodic .gt. 0) then
          ! Get current integration point coordinate
          coord_current(:) = batches(i_batch) % points(i_index) % coords(:)

          if (fast_Ewald) then
            call update_hartree_potential_recip_v2( &
              coord_current, delta_v_hartree(i_full_points) )
          else
            call update_hartree_potential_recip_v2( &
              coord_current, delta_v_hartree(i_full_points) )
          end if ! fast_Ewald

          if ( use_dipole_correction) then

            i_m = int(floor((coord_current(3) - vacuum_z_level) / dip_length))
            dip_coord_current = coord_current(3) - i_m * dip_length

            if (dip_coord_current < vacuum_z_level) then
              dip_coord_current = dip_coord_current + dip_length
            end if

            delta_v_hartree(i_full_points) = delta_v_hartree(i_full_points) - &
              (dip_coord_current - dip_origin) * dip_gradient

          end if ! use_dipole_correction

        else if (use_hartree_non_periodic_ewald) then
          call interpolate_extended_hartree_non_periodic_ewald( &
            batches(i_batch) % points(i_index) % coords, &
            delta_v_hartree(i_full_points) )

        end if ! n_periodic .gt. 0

      end if ! partition_tab(i_full_points) .gt. 0.d0

    end do ! i_index
  end do ! i_batch

  deallocate(current_rho_multipole)
  deallocate(current_rho_multipole_spl)
  deallocate(rho_multipole)
  deallocate(aux_delta_v_hartree_result)
  deallocate(index_lm)
  deallocate(max_lm_spl_rad_sq)
  deallocate(adap_outer_radius_sq)

end subroutine hartree_solver
!******
!------------------------------------------------------------------------------
!****s* FHI-aims/get_rho_multipole_spl_dev
!  NAME
!    get_rho_multipole_spl_dev
!  SYNOPSIS

subroutine get_rho_multipole_spl_dev(rho_multipole_spl, spl_atom, rho_multipole)

  !  PURPOSE
  !    Delivers the spline coefficients for rho_multipole
  !    Performs the same function as get_rho_multipole_spl, but does not entirely relay on
  !    hartree_potential_storage
  !  USES

  use mpi_utilities
  use dimensions
  use grids
  use species_data
  use spline
  use hartree_potential_storage, only : rho_multipole_index
  use geometry, only: species

  implicit none

  !  ARGUMENTS

  real*8, intent(out) :: rho_multipole_spl((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
  integer, intent(in) :: spl_atom
  real*8, intent(in)  :: rho_multipole((l_pot_max+1)**2, n_max_radial+2, n_atoms)

  !  INPUTS
  !    o spl_atom -- the atom for which the spline coefficients are needed
  !    o rho_multipole -- multipole components for rho
  !  OUTPUTS
  !    o rho_multipole_spl -- the spline coefficients
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  HISTORY
  !    Created in November 2017
  !  SOURCE

  integer i_radial, l_h_dim, n_rad, j, i_atom_index

  real*8 :: i_r_outer
  real*8 :: delta, delta_2, delta_3

  ! Get index of spl_atom in rho_multipole and check if rho_multipole for this atom is stored

  i_atom_index = rho_multipole_index(spl_atom)
  if(i_atom_index<=0) then
    print '(2(a,i5))','ID ',myid,' INTERNAL ERROR get_rho_multipole_spl_dev - need atom! Atom: ',spl_atom
    call aims_stop
  endif

  l_h_dim = (l_hartree(species(spl_atom))+1)**2
  n_rad   = n_radial(species(spl_atom))

  ! First, simply produce spline coefficients from the charge density as
  ! given on the radial integration grid.

  rho_multipole_spl(1:l_h_dim,1,1:n_rad+2) = rho_multipole(1:l_h_dim,1:n_rad+2,i_atom_index)

  ! Spline interpolation
  call cubic_spline_v2(rho_multipole_spl, (l_pot_max+1)**2, n_max_spline, n_max_radial+2, &
                       n_rad+2, l_h_dim)


  ! Splines are now tabulated up to r = infinity in principle.
  ! NOW, "doctor" all splines in the far field:
  ! We know that no density must occur outside the radius of the free atom,
  ! because the partition table is zero there.
  ! Therefore, extrapolate from last radial shell inside free atom radius to 
  ! become zero at multipole_radius_free

  ! find outermost radial grid point that is possibly non-zero
  i_radial = n_rad
  do while ( ( r_radial(i_radial,species(spl_atom)) .ge. &
       multipole_radius_free(species(spl_atom)) ) &
       .and.(i_radial.gt.1) )

     rho_multipole_spl(1:l_h_dim,:,i_radial+1) = 0.d0

     i_radial = i_radial - 1
  enddo

  ! Outermost atom radius in units of the radial integration grid
  i_r_outer = invert_radial_grid &
       ( multipole_radius_free(species(spl_atom)), &
       n_rad, &
       scale_radial(species(spl_atom)) )

  delta = dble(i_r_outer - i_radial)
  delta_2 = delta*delta
  delta_3 = delta_2*delta

  ! This is an ugly hack because the element i_radial+1 in rho_multipole_spl
  ! now corresponds to radial shell r_radial(i_radial). What a mess.
  i_radial = i_radial + 1

  ! doctor the spline coefficients at the outermost finite value
  ! i_radial to go smoothly to zero at multipole_radius_free

  do j = 1, l_h_dim

    rho_multipole_spl( j, 3, i_radial) = &
             - 3.d0 / delta_2 * rho_multipole_spl( j, 1, i_radial) &
             - 2.d0 / delta   * rho_multipole_spl( j, 2, i_radial)

    rho_multipole_spl( j, 4, i_radial) = &
               2.d0 / delta_3 * rho_multipole_spl( j, 1, i_radial) &
             + 1.d0 / delta_2 * rho_multipole_spl( j, 2, i_radial)

  enddo
end subroutine get_rho_multipole_spl_dev
!******
