!****h* FHI-aims/v_hartree_multipole_evaluation
!  NAME
!    v_hartree_multipole_evaluation
!  SYNOPSIS

module v_hartree_multipole_evaluation

!  PURPOSE 
!    This module contains several routines that evaluate
!
!     o the multipole expanded hartree potential (as sum of 
!              free + delta contributions) and its gradient
!
!    at arbitrary points in real space.
!
!    However, there is NO check implemented whether requested evaluation
!    points lie on top of atoms / multipole centers. Thus it is the user's
!    responsibility to exclude such points when calling the following routines:
!
!     o get_v_hartree_multipole_and_gradient
!
!  USES

   use types, only: dp
   use mpi_tasks, only: aims_stop
   use runtime_choices, only: &
         communication_type, &
         shmem_comm, &
         force_new_functional, &
         use_hartree_non_periodic_ewald
   use physics, only: &
         multipole_radius_sq, &
         outer_potential_radius, &
         l_hartree_max_far_distance, &
         multipole_moments
   use analytic_multipole_coefficients, only: &
         initialize_analytic_multipole_coefficients, &
         cleanup_analytic_multipole_coefficients
   use hartree_potential_storage, only: get_rho_multipole_spl
   use hartree_potential_real_p0, only: &
         hartree_potential_real_coeff, &
         update_outer_radius_l, &
         hartree_force_l_add, &
         far_distance_hartree_Fp_periodic_single_atom, &
         far_distance_hartree_Fp_cluster_single_atom_p2, &
         far_distance_real_hartree_potential_single_atom, &
         far_distance_real_hartree_potential_single_atom_p2, &
         far_distance_real_gradient_hartree_potential_single_atom, &
         far_distance_real_gradient_hartree_potential_single_atom_p2
   use dimensions, only: &
         l_pot_max, &
         n_atoms, &
         n_centers_hartree_potential, &
         n_max_spline, &
         n_max_radial, &
         n_hartree_grid, &
         n_periodic, &
         use_forces
   use grids, only: n_grid, n_radial
   use geometry, only: species, empty
   use free_atoms, only: free_pot_es_spl
   use spline, only: &
         spline_vector_v2, &
         spline_deriv_vector_v2, &
         val_spline, &
         val_spline_deriv
   use pbc_lists, only: &
         centers_hartree_potential, &
         species_center, &
         center_to_atom
   use species_data, only: &
         species_z, &
         l_hartree

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   implicit none

   ! make everything private by default
   private

   ! specifically declare public routines
   public :: get_v_hartree_multipole_and_gradient
   public :: get_v_nuc_and_gradient

contains


subroutine get_v_nuc_and_gradient( &
      n_points, points, &
      v_nuc, v_nuc_gradient )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded hartree potential 
!    and its gradient at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real(dp), intent(in) :: points(3,n_points)

   real(dp), intent(out) :: v_nuc(n_points)

   real(dp), intent(out) :: v_nuc_gradient(3,n_points)

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where density and density
!               gradient are evaluated
!  OUTPUT
!   o v_nuc -- electrostatic potential of all nuclei at points
!   o v_nuc_gradient -- gradients of the nulear electrostatic potential
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   integer :: i_center
   integer :: current_center
   integer :: current_spl_atom
   integer :: i_point

   real(dp) :: dist_sq, inv_dist_sq, inv_dist
   real(dp) :: relvec(3)

   v_nuc = 0.e0_dp

   do i_center = 1, n_centers_hartree_potential

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)

      if (.not.empty(current_spl_atom)) then
         do i_point = 1, n_points

            call tab_single_atom_centered_coords_p0( &
                 current_center, &
                 points(:,i_point),  &
                 dist_sq,  &
                 relvec )
            inv_dist_sq = 1/dist_sq
            inv_dist = sqrt(inv_dist_sq)

            v_nuc(i_point) = v_nuc(i_point) + &
                  species_z(species_center(current_center)) * inv_dist

            v_nuc_gradient(:,i_point) = v_nuc_gradient(:,i_point) - &
                  species_z(species_center(current_center)) * &
                  inv_dist_sq * inv_dist * relvec(:)
         enddo ! i_point
      endif
   enddo ! i_center

end subroutine get_v_nuc_and_gradient


subroutine get_v_hartree_multipole_and_gradient( &
      n_points, points, &
      v_hartree_mp, v_hartree_mp_gradient )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded hartree potential 
!    and its gradient at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real(dp), intent(in) :: points(3,n_points)

   real(dp), intent(out) :: v_hartree_mp(n_points)

   real(dp), intent(out), optional :: v_hartree_mp_gradient(3,n_points)

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where density and density
!               gradient are evaluated
!  OUTPUT
!   o v_hartree_mp -- values of the multipole expanded Hartree potential
!   o v_hartree_mp_gradient -- gradients of the multipole expanded Hartree potential
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   integer, parameter :: n_coeff_hartree = 2 ! # of spline coeffs for current_delta_v_hart_part_spl

   real(dp), external :: ddot

   logical :: forces_on

   integer :: atom_of_splines
   integer :: i_center
   integer :: current_center
   integer :: current_spl_atom
   integer :: i_point
   integer :: l_atom_max
   integer :: n_bytes
   integer :: l_h_dim
   integer :: i_index
   integer :: i_l, i_m
   integer :: index_lm(-l_pot_max:l_pot_max, 0:l_pot_max)
   integer :: i_atom

   real(dp) :: v_hartree_gradient_temp(3)
   real(dp) :: coord_current(3)
   real(dp) :: dist_tab_sq
   real(dp) :: dist_tab_in
   real(dp) :: dist_tab_out
   real(dp) :: dir_tab(3)
   real(dp) :: dir_tab_in(3)
   real(dp) :: dir_tab_out(3)
   real(dp) :: i_r, i_r_log
   real(dp) :: log_weight
   real(dp) :: radial_weight
   real(dp) :: trigonom_tab(4)
   real(dp) :: ylm_tab((l_pot_max+1)**2)
   real(dp) :: delta_v_hartree_multipole_component((l_pot_max+1)**2)
   real(dp) :: v_hartree_free_aux
   real(dp) :: delta_v_hartree_aux
   real(dp) :: d_v_hartree_free_d_r

   real(dp), allocatable :: adap_outer_radius_sq(:)
   real(dp), allocatable :: current_rho_multipole_spl(:,:,:)
   real(dp), allocatable :: current_delta_v_hart_part_spl(:,:,:)
   real(dp), allocatable :: dylm_dtheta_tab(:,:)
   real(dp), allocatable :: scaled_dylm_dphi_tab(:,:)
   real(dp), allocatable :: delta_v_hartree_multipole_deriv(:)


   if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
      call aims_stop("Not implemented", "get_v_hartree_multipole_and_gradient")
   endif

   forces_on = present(v_hartree_mp_gradient)

   allocate(adap_outer_radius_sq(n_atoms))
   allocate(current_rho_multipole_spl &
          ((l_pot_max+1)**2, n_max_spline, n_max_radial+2))
   allocate(current_delta_v_hart_part_spl &
          ((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid))
   current_delta_v_hart_part_spl = 0.e0_dp
   if (forces_on) then
      allocate(dylm_dtheta_tab((l_pot_max+1)**2, n_centers_hartree_potential))
      allocate(scaled_dylm_dphi_tab((l_pot_max+1)**2, n_centers_hartree_potential))
      allocate(delta_v_hartree_multipole_deriv((l_pot_max+1)**2))
   endif


   if (forces_on) then
      ! MS: It really breaks my heart to write this here. Seriously, was there
      ! no better way to design this in the first place???
      if (.not. use_forces) then
         call cleanup_analytic_multipole_coefficients()
         use_forces = .true.
         call initialize_analytic_multipole_coefficients()
         use_forces = .false.
      endif
      hartree_force_l_add = 1
   else
      hartree_force_l_add = 0
   endif

   ! MS: To be on the safe side, we repeat the following initialization
   i_index = 0
   do i_l = 0, l_pot_max, 1
      do i_m = -i_l, i_l
         i_index = i_index + 1
         index_lm(i_m, i_l) = i_index
      enddo
   enddo

   call hartree_potential_real_coeff(index_lm, multipole_moments, &
         l_hartree_max_far_distance, n_centers_hartree_potential )

   if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
      call update_outer_radius_l( outer_potential_radius, multipole_moments, multipole_radius_sq, & 
                                 l_hartree_max_far_distance, index_lm)

      if ( .not. use_hartree_non_periodic_ewald ) then
         do i_atom = 1, n_atoms
            adap_outer_radius_sq(i_atom) = maxval(outer_potential_radius(:,i_atom))
         end do
      else
         adap_outer_radius_sq(:) = multipole_radius_sq(:)
      end if
   else
      adap_outer_radius_sq = 1e8_dp
   end if


   v_hartree_mp = 0.e0_dp
   if (forces_on) v_hartree_mp_gradient = 0.e0_dp

   atom_of_splines = 0
   do i_center = 1, n_centers_hartree_potential, 1

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)

      do i_point = 1, n_points

         v_hartree_gradient_temp = 0.e0_dp

         ! get current integration point coordinate
         coord_current(:) = points(:,i_point)

         ! Tabulate distances and directions to a single atom -
         ! including relative positions on logarithmic and
         ! radial integration grids.

         call tab_single_atom_centered_coords_p0 &
              ( current_center, &
              coord_current,  &
              dist_tab_sq,  &
              dir_tab )

         ! VB: For uniformity, determine the maximum angular momentum required for the present atom 
         ! right here

         l_atom_max = l_hartree(species(current_spl_atom))
         do while ( (outer_potential_radius(l_atom_max, current_spl_atom) .lt. dist_tab_sq ) & 
              .and. (l_atom_max.gt.0) ) 
            l_atom_max = l_atom_max - 1
         enddo
         ! At each integration point, the Hartree potential components coming from
         ! different atoms are split into two groups:
         ! Any components coming from atoms close by are evaluated by explicit
         ! numerical splines; far away atoms are simply represented by
         ! an analytical long-distance multipole potential.

         ! We treat first the atoms close by
         if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom)) then

            if(current_spl_atom /= atom_of_splines) then
               call get_rho_multipole_spl(current_rho_multipole_spl, current_spl_atom)

               if(communication_type.eq.shmem_comm) then
                  n_bytes = (l_pot_max+1)**2 * n_coeff_hartree * n_hartree_grid * 8
                  call aims_shm_get(current_delta_v_hart_part_spl, (current_spl_atom-1)*n_bytes, n_bytes)
               else
                  call integrate_delta_v_hartree( current_rho_multipole_spl, current_delta_v_hart_part_spl, &
                                             n_coeff_hartree, current_spl_atom )
               endif
               atom_of_splines = current_spl_atom
            endif

            call tab_single_atom_centered_coords_radial_log_p0 &
                ( current_center, dist_tab_sq, dir_tab,  &
                dist_tab_in, i_r, i_r_log, dir_tab_in )

            if (forces_on) then 
               ! need radial derivatives di/dr for the benefit of
               ! spline derivative evaluation later
               call tab_single_radial_weights_v2 &
                  ( current_spl_atom, dist_tab_in, i_r, &
                  log_weight, radial_weight )
            end if

            ! for an inner atom we need ylm functions and their gradients explicitly

            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

            if (forces_on) then 
               call tab_single_gradient_ylm_p2 &
                  ( trigonom_tab, l_atom_max, l_pot_max,  &
                  ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
            else
               call tab_single_wave_ylm_p2 &
                  ( trigonom_tab, l_atom_max,  &
                  l_pot_max, ylm_tab)
            end if

            ! For an inner atoms (those which are close enough
            ! to the current integration point so we need explicit numerical splines)
            !     partitioned
            !     hartree potentials need to be summed up
            !     according to Delley (eq. 12c)
            l_h_dim = (l_atom_max + 1)**2

            ! obtain spline-interpolated values of the multipole components
            ! of the partitioned Hartree potential, splined on the logarithmic
            ! integration grid
            call spline_vector_v2 &
                 ( i_r_log, &
                 current_delta_v_hart_part_spl, &
                 (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, &
                 n_grid(species_center(current_center)), &
                 l_h_dim, &
                 delta_v_hartree_multipole_component)

            ! MS: calculate superposition of free atoms
            v_hartree_free_aux = val_spline( &
               i_r_log, free_pot_es_spl(1,1,species_center(current_center)), &
               n_grid(species_center(current_center)) )

            delta_v_hartree_aux = &
                  ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab, 1 )

            ! MS: add up contributions to total hartree potential
            v_hartree_mp(i_point) = v_hartree_mp(i_point) + &
                  v_hartree_free_aux + delta_v_hartree_aux

            ! If needed, obtain multipole force correction for present inner atom
            if (forces_on) then 
               ! get all needed derivatives for forces from splines

               call spline_deriv_vector_v2 &
                  ( i_r_log, &
                  current_delta_v_hart_part_spl, &
                  (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid,  &
                  n_grid(species_center(current_center)), &
                  l_h_dim, &
                  delta_v_hartree_multipole_deriv)

               ! splines are now derivatives df/di where i is the grid point index
               ! must convert to df/dr = df/di * di/dr
               delta_v_hartree_multipole_deriv(1:l_h_dim) = &
                     delta_v_hartree_multipole_deriv(1:l_h_dim) * log_weight

               ! obtain gradients of the free atom density and potential explicitly

               ! radial derivative of hartree potential of free atom # i_atom (without core potential!!!)
               ! subtract core potential => - ( -Z/r ) = + Z/r , d/dr => - Z/r^2

               if(n_periodic == 0 .and.(.not. force_new_functional))then
                  d_v_hartree_free_d_r =  &
                     val_spline_deriv( i_r_log, &
                     free_pot_es_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom))) * log_weight &
                     - species_z(species(current_spl_atom)) / (dist_tab_sq)
               else
                  d_v_hartree_free_d_r = &
                     val_spline_deriv( i_r_log, &
                     free_pot_es_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom))) * log_weight
               end if

               ! obtain gradients for present atom ...

               call evaluate_v_hartree_gradient &
                   ( dist_tab_in, dir_tab_in,  &
                   trigonom_tab,  &
                   ylm_tab,  &
                   dylm_dtheta_tab,  &
                   scaled_dylm_dphi_tab,  &
                   delta_v_hartree_multipole_component,  &
                   delta_v_hartree_multipole_deriv, &
                   d_v_hartree_free_d_r, &
                   l_h_dim, &
                   v_hartree_gradient_temp )

            end if ! (forces_on)

            ! far-distance treatment for then center i_center for the periodic system
            if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then

               !  VB: FIXME!!!!! 
               !      I kept the original call here because I was not sure I was doing the completely
               !      right thing. We should cut the periodic part down to l_atom_max as well below,
               !      but that requires a change to the periodic infrastructure first!
               !
               call far_distance_hartree_Fp_periodic_single_atom &
                   (current_spl_atom, i_center, &
                   dist_tab_in, l_hartree_max_far_distance, .true., forces_on,  &
                   multipole_radius_sq(current_spl_atom),                       &
                   sqrt( adap_outer_radius_sq(current_spl_atom) )                )

               ! VB: Paula - this is what the function call should look like but the function
               !             needs to be changed!
               !                    call far_distance_hartree_Fp_periodic_single_atom &
               !                        (i_center, &
               !                        dist_tab_in, l_atom_max, .true., forces_on )
            end if


         else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then

            ! the current center is in the far-distance part-----------------
            ! Tabulate distance for the outer part
            dist_tab_out = sqrt(dist_tab_sq)

            dir_tab_out = dir_tab / dist_tab_out


            ! Now sum up the potential contributions from a far-field atom
            ! (analytical multipole potentials only ...)
            if ( n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald ) then

               ! if not periodic, need electronic v_hartree_free part for total energy

               v_hartree_free_aux = 0.e0_dp
               delta_v_hartree_aux = 0.e0_dp

               if (.not.(empty(current_spl_atom))) then
                  if(.not. force_new_functional)then
                     v_hartree_free_aux = &
                           species_z(species(current_spl_atom)) / dist_tab_out
                  end if
               end if

               ! Recursively tabulate the radial behavior of all multipole components
               ! These radial functions are a private variable in module
               ! hartree_potential_real.f90, and are reused there later in
               ! subroutine far_distance_hartree_potential_real_single_atom
               call far_distance_hartree_Fp_cluster_single_atom_p2 &
                    ( dist_tab_out, &
                    l_atom_max, forces_on )

               call far_distance_real_hartree_potential_single_atom_p2 &
                         ( i_center, delta_v_hartree_aux, &
                         l_atom_max, coord_current )


               v_hartree_mp(i_point) = v_hartree_mp(i_point) + &
                     v_hartree_free_aux + delta_v_hartree_aux

               ! ... and compute force contribution of a far-distance atom
               if (forces_on) then 

                  ! All that's left of the free Hartree potential is the
                  ! nuclear z/r part

                  if(.not. force_new_functional)then
                     d_v_hartree_free_d_r =  &
                           - species_z(species(current_spl_atom))/dist_tab_sq
                  else
                     d_v_hartree_free_d_r = 0.e0_dp
                  end if

                  ! VB: FIXME: This needs to go for the new functional part!
                  v_hartree_gradient_temp(:) =  &
                           dir_tab_out(:) * d_v_hartree_free_d_r

                  ! add multipole potential gradients to free-atom gradient
                  ! The periodic case is calculated later because it uses both outer and inner atoms.
                  call far_distance_real_gradient_hartree_potential_single_atom_p2 &
                     ( current_spl_atom, dir_tab, &
                     v_hartree_gradient_temp, &
                     l_atom_max )                           

               end if ! (forces_on)

            else ! Periodic system or non-periodic Ewald method

               ! nuclear z/r part
               !  d_v_hartree_free_d_r =  &
               !       - 0*species_z(species(current_spl_atom))/dist_tab_sq

               ! d_v_hartree_free_d_r = 0.d0

               ! v_hartree_gradient(:,i_center,i_point) =  &
               !      dir_tab_out(:) * d_v_hartree_free_d_r

               call far_distance_hartree_Fp_periodic_single_atom &
                  (current_spl_atom, i_center, &
                  dist_tab_out, l_hartree_max_far_distance, .false.,forces_on,  &
                  multipole_radius_sq(current_spl_atom),                        &
                  sqrt( adap_outer_radius_sq(current_spl_atom) )                 )

               ! Note: 'far_distance_real_hartree_potential_single_atom' is called in the periodic
               ! systems later, because then we have Fp in inner and outer multipole radius.

            end if  ! n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald

         end if  ! end if for separating current_center either far-distance or near part


         ! Now sum up the far distance part of the Hartree potential if Ewald's
         ! decomposition is performed. In the opposite case, this was already done.
         if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then


            if (dist_tab_sq.lt. &
                max(adap_outer_radius_sq(current_spl_atom),multipole_radius_sq(current_spl_atom))) then

               ! Far distance analytic real-space part of the potential

               ! VB: FIXME: l_atom_max not implemented here although it should be!!

               delta_v_hartree_aux = 0.e0_dp
               call far_distance_real_hartree_potential_single_atom &
                       ( current_center, i_center, delta_v_hartree_aux, &
                       l_hartree_max_far_distance, coord_current )

               v_hartree_mp(i_point) = v_hartree_mp(i_point) + &
                     delta_v_hartree_aux

               if (forces_on) then

                  ! Far distance analytic real-space part of the potential gradient

                  call far_distance_real_gradient_hartree_potential_single_atom &
                      ( current_spl_atom, i_center, dir_tab, &
                      v_hartree_gradient_temp, &
                      l_hartree_max_far_distance, &
                      dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )

               end if ! (forces_on)

            end if !dist_tab .lt....

         end if ! end periodic system

         if (forces_on) then
            v_hartree_mp_gradient(:,i_point) = &
               v_hartree_mp_gradient(:,i_point) + v_hartree_gradient_temp(:)
         endif

      end do ! i_point

   end do ! i_center

end subroutine

end module v_hartree_multipole_evaluation

