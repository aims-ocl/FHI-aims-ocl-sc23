!****s* FHI-aims/evaluate_potential_cube
!  NAME
!   evaluate_potential_cube
!  AUTHOR
!   Stefan Ringe
!  REFERENCE
!   This routine is part of the MPBE solver for the modeling of 
!   electrolytes with FHI-aims. Large parts of this code are still 
!   experimental, that is why we highly recommend you to contact
!   the main authors of the corresponding paper (below) before starting
!   to calculate anything. We are highly motivated to help and 
!   cooperate with interested FHI-aims users.
!  SEE ALSO
!    Ringe, S., Oberhofer, H., Hille, C., Matera, S., Reuter, K., "Function-Space oriented Solution Scheme for the 
!    Size-Modified Poisson-Boltzmann Equation in Full-Potential DFT", JCTC, 2016, accepted.
!    DOI: 10.1021/acs.jctc.6b00435
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2016).
!  SYNOPSIS

subroutine evaluate_potential_cube(n_cube_points,points,gradient_output,potential_cube,gradient_cube_abs)
!  PURPOSE
!    Evaluate the delta v potential on the cubic grid from the multipole components that were calculated
!    before in the SPE solver. This is maybe not implemented in the best way and can cost a lot of time 
!    so try first with smaller cube files.
!  SOURCE

!purpose: evaluate delta v by the use of the splined multipole moments
!obtained by kerker_mpb.f90
!on the grid required for cube output
!for more information look into kerker_mpb sum_up_mpe

!currently only supported for solve_in_mpe = false
!TODO: put this into the kerker_mpb

  use grids
  use pbc_lists
  use constants
  use dimensions
  use synchronize_mpi
  use physics, only:  rho  
  use SPE_solver, only: R_multipole_spl_mpb, &
    limit_l_max_SPE_species,max_lm_spl_rad_sq_mpb, index_lm_mpb, partition_tab_mpb,&
    aux_R_prec_result_mpb,limit_l_dim_SPE_species,limit_l_dim_SPE_species_current,n_max_spline_mpb,&
    n_spline_grid_dim
  use runtime_choices, only: extra_adding_to_hartree_potential_distance, communication_type,shmem_comm
 use lpb_solver_utilities, only: l_dim_SPE_public
  use hartree_potential_storage, only : get_rho_multipole_spl
  use hartree_potential_real_p0
  use hartree_potential_recip
  use hartree_non_periodic_ewald
  use physics, only : multipole_radius_sq_public
  use spline, only: spline_vector, spline_vector_v2, spline_deriv_vector_v2
  use localorb_io, only: use_unit

  implicit none

  integer :: n_cube_points
  integer :: i_cube_points, prev_atom, current_spl_atom, current_center, i_l,&
    i_m, i_grid, i_spline, i_center
  real*8 :: max_spl_rad_sq, delta, multipole_thresh
  real*8, dimension(n_cube_points) :: potential_cube
  real*8, dimension(3,n_cube_points) :: gradient_cube
  real*8, dimension(n_cube_points) :: gradient_cube_abs
  real*8, dimension(3,n_cube_points) :: points
  logical :: set_mult_rad_to_addtomrad
  integer :: precondition_max_l_mpb
  integer :: l_atom_max
  
  real*8 :: dist_tab_sq
  real*8 :: dist_tab
  real*8, dimension(3) :: dir_tab
  real*8, dimension(3) :: coord_current
  real*8 :: radius_rad, log_weight, radial_weight, prec_times_rgrid,&
    trigonom_tab(4), ddot, radius
  integer :: lmax_cut = 0
  real*8, dimension(:), allocatable :: ylm_tab 
  real*8, dimension((l_pot_max+1)**2) :: ylm_tab2 

  real*8, dimension(n_species) :: addtomprad!100.0
  logical :: use_analytic_far_field = .true.
  integer n_bytes
  integer, parameter :: n_coeff_hartree = 2
  integer :: info
  integer l_h_dim
  real*8, dimension(3) :: dir_tab_in,dir_tab_out
  real*8 :: radius_cut = 15d0
 ! real*8, dimension(n_full_points) :: test_pot

  integer :: i_my_batch, i_index
  real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl
  real*8, dimension(:,:,:), allocatable :: current_delta_v_hart_part_spl
  real*8,  dimension(:,:,:)   , allocatable :: current_R_multipole_spl_mpb ! to store the multipole splines of one given atom for calculation of final answer
  real*8 i_r_log,i_r
  real*8, dimension((l_pot_max+1)**2) :: delta_v_hartree_multipole_component
  real*8 dist_tab_in,dist_tab_out
  real*8, dimension(:), allocatable :: delta_v_hartree_multipole_deriv
  real*8 :: d_v_hartree_free_d_r
  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
  logical , intent(in):: gradient_output
  real*8, dimension(3) :: v_hartree_gradient_temp
  real*8, dimension(n_atoms) :: adap_outer_radius_sq
  logical :: forces_on 
  integer :: i,j,k
  
  forces_on = gradient_output


  if (solvent_method.eq.SOLVENT_MPB) then
    if (.not.allocated(current_R_multipole_spl_mpb)) allocate(current_R_multipole_spl_mpb(l_dim_SPE_public, n_max_spline_mpb, n_spline_grid_dim))
    if (.not.allocated(ylm_tab)) allocate(ylm_tab(l_dim_SPE_public))
    current_R_multipole_spl_mpb = 0d0
  else
!     write(use_unit,*) 'l_pot_max',l_pot_max
!     if (.not.allocated(ylm_tab)) then
!       allocate(ylm_tab((l_pot_max+1)**2)) !(l_pot_max+1)**2))
!       call check_allocation(info, 'ylm_tab ')
!       ylm_tab = 0d0
!     end if
    if (.not.allocated(current_delta_v_hart_part_spl)) then
      allocate(current_delta_v_hart_part_spl &
	    ((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid),stat=info)
      call check_allocation(info, 'current_delta_v_hart_part_spl ')
      current_delta_v_hart_part_spl = 0.d0
    end if
    if (.not.allocated(current_rho_multipole_spl)) then
      allocate(current_rho_multipole_spl &
	    ((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=info)
      call check_allocation(info, 'current_rho_multipole_spl     ')
    end if
    current_delta_v_hart_part_spl = 0d0

    if (gradient_output) then
     if (.not.allocated(delta_v_hartree_multipole_deriv)) then
        allocate(delta_v_hartree_multipole_deriv((l_pot_max+1)**2),stat=info)
        call check_allocation(info, 'delta_v_hartree_multipole_deri')
     end if
      if (.not.allocated(dylm_dtheta_tab)) then
	  allocate(dylm_dtheta_tab((l_pot_max+1)**2, n_centers_hartree_potential),stat=info)
	  call check_allocation(info, 'dylm_dtheta_tab               ')
      end if
      if (.not.allocated(scaled_dylm_dphi_tab)) then
	  allocate(scaled_dylm_dphi_tab((l_pot_max+1)**2,n_centers_hartree_potential),stat=info)
	  call check_allocation(info, 'scaled_dylm_dphi_tab          ')
      end if
     end if
  end if
  adap_outer_radius_sq = 1e8

!   if (use_analytic_far_field) then
!     set_mult_rad_to_addtomrad = .true.
!     do i_species = 1, n_species, 1
!       addtomprad(i_species) = (multipole_radius_free_SPE(i_species) & 
! 	+extra_adding_to_hartree_potential_distance)**2 + 20d0
!     end do
!   end if  
!     
!   set_mult_rad_to_addtomrad = .false.
!   multipole_thresh = 1d-15
  
  !if (myid==0) write(use_unit,*) 'starting to evaluate potential cube'
  !the outer loop goes over all multipole centers, the inner over all cube points

  potential_cube = 0d0
  gradient_cube = 0d0

  do i_center = 1, n_centers_hartree_potential, 1

  
      prev_atom        = current_spl_atom
      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)                     ! this is the actual atom number, .ne. current_center for PBC's
      
     if (current_spl_atom.ne.prev_atom) then
	  ! distribute the information about that particular atom to all the threads
	  ! simultaneously, each one is going to use it in order to compute its part 
	  ! of the new residual ... 

	  if (solvent_method.eq.SOLVENT_MPB) then
	    call get_R_prec_spline_mpi(current_spl_atom,l_dim_SPE_public,&
	      current_R_multipole_spl_mpb,R_multipole_spl_mpb, n_spline_grid_dim)	
	  else
	    call get_rho_multipole_spl(current_rho_multipole_spl, current_spl_atom)
	    call integrate_delta_v_hartree( current_rho_multipole_spl, current_delta_v_hart_part_spl, &
					    n_coeff_hartree, current_spl_atom )
	  end if
     end if !current_spl_atom.ne.prev_atom

!       i_cube_points = 0
!       do i_my_batch = 1, n_my_batches, 1
! 	    do i_index = 1, batches(i_my_batch)%size, 1
! 		i_cube_points   = i_cube_points + 1
! 		  coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

     do i_cube_points = 1,n_cube_points,1
	coord_current(:) = points(1:3,i_cube_points)
	call  tab_single_atom_centered_coords_p0( current_center, coord_current, &
	  dist_tab_sq, dir_tab )


	  
	if (solvent_method.eq.SOLVENT_MPB) then

	    call tab_single_trigonom_p0( dir_tab, trigonom_tab ) 

	    call tab_single_wave_ylm_p0(current_center,trigonom_tab, &
	      limit_l_max_SPE_species,precondition_max_l_mpb, ylm_tab )

	    radius = SQRT(dist_tab_sq)    ! this is the actual distance from point to atom: spline this ... 
	    if (radius.gt.radius_cut.and.lmax_cut.gt.0) then
	      limit_l_dim_SPE_species_current = (lmax_cut + 1)**2
	    else
	      limit_l_dim_SPE_species_current = limit_l_dim_SPE_species
	    end if! 			end if
	    radius = invert_log_grid(radius, &
	      r_grid_min(species(current_spl_atom)), &
	      r_grid_inc(species(current_spl_atom)))

	    call spline_vector                             &
	      ( radius,                                 &
		current_R_multipole_spl_mpb,                &  
		n_spline_grid_dim,                      &
		l_dim_SPE_public,                             &   ! FIXME: this is only necessary for some maximum l, should be faster
		n_grid(species_center(current_center)), &
		limit_l_dim_SPE_species(species_center(current_center)), &
		aux_R_prec_result_mpb)
	    potential_cube(i_cube_points) = potential_cube(i_cube_points) + &
	      ddot ( limit_l_dim_SPE_species_current(species_center(current_center)), &
	      aux_R_prec_result_mpb(:), 1, ylm_tab(:), 1 ) 
	else
	    l_atom_max = l_hartree(species(current_spl_atom))
!             do while ( (outer_potential_radius(l_atom_max, current_spl_atom) .lt. dist_tab_sq ) & 
!                  .and. (l_atom_max.gt.0) ) 
!               l_atom_max = l_atom_max - 1
!             enddo
 	  if (dist_tab_sq.lt.multipole_radius_sq_public(current_spl_atom)) then
! 	    write(use_unit,*) 'coord_current', coord_current , dist_tab_sq
	    
	    call tab_single_atom_centered_coords_radial_log_p0 &
	      ( current_center, dist_tab_sq, dir_tab,  &
	      dist_tab_in, i_r, i_r_log, dir_tab_in )
	      
! 	    write(use_unit,*) 'irs', dist_tab_in, i_r, i_r_log,dir_tab_in
	    
	    if (gradient_output) then
                call tab_single_radial_weights_v2 &
                     ( current_spl_atom, dist_tab_in, i_r, &
                     log_weight, radial_weight )
	    end if
	    
! 	    write(use_unit,*) 'weights', log_weight, radial_weight 	    

            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

!             write(use_unit,*) 'trigonom_tab', trigonom_tab
            
              if (gradient_output) then 
                call tab_single_gradient_ylm_p2 &
                     ( trigonom_tab, l_atom_max, l_pot_max,  &
                     ylm_tab2, dylm_dtheta_tab, scaled_dylm_dphi_tab)
              else
                call tab_single_wave_ylm_p2 &
                     ( trigonom_tab, l_atom_max,  &
                     l_pot_max, ylm_tab2)
              end if

            l_h_dim = (l_atom_max + 1)**2

            write(use_unit,*) 'current_delta_v_hart_part_spl', current_delta_v_hart_part_spl(1,1,1:5)
            
            call spline_vector_v2 &
	      ( i_r_log, &
	      current_delta_v_hart_part_spl, &
	      (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, &
	      n_grid(species_center(current_center)), &
	      l_h_dim, &
	      delta_v_hartree_multipole_component)

	    if (gradient_output) then
	      call spline_deriv_vector_v2 &
		    ( i_r_log, &
		    current_delta_v_hart_part_spl, &
		    (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid,  &
		    n_grid(species_center(current_center)), &
		    l_h_dim, &
		    delta_v_hartree_multipole_deriv)

	      delta_v_hartree_multipole_deriv(1:l_h_dim) = &
		    delta_v_hartree_multipole_deriv(1:l_h_dim) * log_weight

	      d_v_hartree_free_d_r =  0d0

	      v_hartree_gradient_temp = 0.0d0
	      call evaluate_v_hartree_gradient &
		    ( dist_tab_in, dir_tab_in,  &
		    trigonom_tab,  &
		    ylm_tab2,  &
		    dylm_dtheta_tab,  &
		    scaled_dylm_dphi_tab,  &
		    delta_v_hartree_multipole_component,  &
		    delta_v_hartree_multipole_deriv, &
		    d_v_hartree_free_d_r, &
		    l_h_dim, &
		    v_hartree_gradient_temp)
	      gradient_cube(:,i_cube_points) = &
		gradient_cube(:,i_cube_points) + v_hartree_gradient_temp
	    else
	      potential_cube(i_cube_points) = &
		potential_cube(i_cube_points) + &
		ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab2, 1 )
	    end if !gradient_output
	    stop
	  else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then
	    if (gradient_output) then
	      dist_tab_out = sqrt(dist_tab_sq)
	      dir_tab_out = dir_tab / dist_tab_out

	      call far_distance_hartree_Fp_cluster_single_atom_p2 &
		    ( dist_tab_out, &
		    l_atom_max, forces_on )

	      v_hartree_gradient_temp =  0.d0

	      call far_distance_real_gradient_hartree_potential_single_atom_p2 &
		    ( current_spl_atom, dir_tab, &
		    v_hartree_gradient_temp, &
		    l_atom_max )    

	      gradient_cube(:,i_cube_points) = &
		gradient_cube(:,i_cube_points) + v_hartree_gradient_temp
	    end if
	  end if !far field distiguishing
	end if !solvent or vacuum calculation?
    enddo      !cube points
   end do !loop over all centers

   gradient_cube_abs = gradient_cube(1,:) !sqrt(sum(gradient_cube**2,DIM=1))

!   if (solvent_method.eq.SOLVENT_MPB) then
!     if (allocated(current_R_multipole_spl_mpb)) deallocate(current_R_multipole_spl_mpb)
!   else
!     if (allocated(current_delta_v_hart_part_spl)) deallocate(current_delta_v_hart_part_spl)
!     if (allocated(current_rho_multipole_spl)) deallocate(current_rho_multipole_spl)
!   end if   
!    
!   if (allocated(ylm_tab)) deallocate(ylm_tab)

end subroutine evaluate_potential_cube
