!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/sum_up_whole_potential
!  NAME
!   Here I remove the variables that I do not want, e.g, free_part. Only loop over atoms.
!  SYNOPSIS

subroutine integrate_Born_effective_charges_Hellman_Feynman &
     ( delta_v_hartree_part_at_zero,  delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, &
     multipole_radius_sq, &
     l_hartree_max_far_distance, &
     outer_potential_radius, &
     hellman_feynman_forces, &
     test_Born_effective_charges_by_print_force)

!  PURPOSE
!    The subroutine get Helman-Feynman part used in Born-effect Charge.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use free_atoms
  use spline
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use hartree_potential_real_p0
  use hartree_potential_recip
  use hartree_non_periodic_ewald
  use pbc_lists
  ! rho_multipole from hartree_potential_storage conflicts with the local variable here,
  ! so just only use get_rho_multipole_spl() from this module.
  ! Maybe the local variable name should be changed ...
  use hartree_potential_storage, only : get_rho_multipole_spl

  implicit none

!  ARGUMENTS


  real*8, dimension(n_atoms)                     :: delta_v_hartree_part_at_zero
  real*8, dimension(3, n_atoms)                  :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments

  real*8, dimension(n_atoms)              :: multipole_radius_sq
  integer, dimension( n_atoms)            :: l_hartree_max_far_distance
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius
  real*8, dimension(3, n_atoms)           :: hellman_feynman_forces
  logical :: test_Born_effective_charges_by_print_force

!  INPUTS
! o delta_v_hartree_part_at_zero -- Hartree potential at origin of the atoms from different multipole components
! o delta_v_hartree_deriv_l0_at_zero -- Derivative of Hartree potential at origin of the atoms
! o multipole_moments -- mutlpole moments of the Hartree potential
! o partition_tab_std -- values of partition function
! o rho_std -- electron density
! o multipole_radius_sq -- outer radius of multipole components
! o l_hartree_max_far_distance -- maximum l-components of for the far distance Hartree potential (periodic systems)
! o rho_multipole_old_std -- multipole components for delta charge
!
!  OUTPUT
! o potential_std -- Hartree potential
! o hellman_feynman_forces -- Helman-Feyman force components
! o outer_potential_radius -- outer radius of the real part of the hartree potential (periodic systems)
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

  real*8 :: hartree_multipole_error


  !  local variables

  integer index_lm(-l_pot_max:l_pot_max, 0:l_pot_max )

  real*8 coord_current(3)
  real*8 dist_tab_sq
  real*8 dist_tab_in
  real*8 dist_tab_out
  real*8 dir_tab(3)
  real*8 dir_tab_in(3)
  real*8 dir_tab_out(3)
  real*8 log_weight
  real*8 radial_weight
  real*8 trigonom_tab(4)
  real*8 i_r
  real*8 i_r_log
  real*8 ylm_tab((l_pot_max+1)**2)



  !     for spline_vector
  real*8, dimension((l_pot_max+1)**2) :: delta_v_hartree_multipole_component
  real*8, dimension((l_pot_max+1)**2) :: rho_multipole_component
  integer l_h_dim

  integer, parameter :: n_coeff_hartree = 2 ! Number of spline coeffs for current_delta_v_hart_part_spl

  real*8, dimension(:), allocatable :: delta_v_hartree_multipole_deriv
  real*8, dimension(:), allocatable :: rho_multipole_deriv
  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
  real*8 :: v_hartree_gradient_temp(3)
  real*8 :: rho_multipole_gradient_temp(3)

  real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl
  real*8, dimension(:,:,:), allocatable :: current_delta_v_hart_part_spl

  real*8, dimension(:),allocatable :: adap_outer_radius_sq
   
  real*8 d_v_hartree_free_d_r
  real*8 delta_v_hartree_aux

  integer :: current_spl_atom
  integer :: current_center
  integer :: atom_of_splines

  integer :: l_atom_max

  !  counters

  integer i_atom_2
  integer i_center
  integer i_batch
  integer i_l
  integer i_m
  integer i_coord
  integer i_index
  integer i_lm
  integer :: i_spin

  integer :: i_lat
  integer :: i_full_points

  !  external functions
  real*8, external :: ddot
  character*100 :: info_str
  real*8:: HF_forces(3)

  real*8 :: HF_temp(3,n_atoms)

  integer :: mpierr
  integer :: info



   integer n_bytes

  ! Load balancing stuff

  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array

  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:)
  real*8, pointer :: potential(:)

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


   call localorb_info("Integrating Born_effective_charges_HF matrix.", use_unit,'(2X,A)', OL_norm )

    hartree_force_l_add = 1


  allocate(adap_outer_radius_sq(n_atoms),stat=info)
  call check_allocation(info, 'adap_outer_radius_sq          ')


  if (.not.allocated(current_rho_multipole_spl)) then
     allocate(current_rho_multipole_spl &
          ((l_pot_max+1)**2, n_max_spline, n_max_radial+2),stat=info)
     call check_allocation(info, 'current_rho_multipole_spl     ')

  end if

  if (.not.allocated(current_delta_v_hart_part_spl)) then
     allocate(current_delta_v_hart_part_spl &
          ((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid),stat=info)
     call check_allocation(info, 'current_delta_v_hart_part_spl ')

  end if

  if (.not.allocated(delta_v_hartree_multipole_deriv)) then
     allocate(delta_v_hartree_multipole_deriv((l_pot_max+1)**2),stat=info)
     call check_allocation(info, 'delta_v_hartree_multipole_deri')

  end if
  if (.not.allocated(rho_multipole_deriv)) then
     allocate(rho_multipole_deriv((l_pot_max+1)**2),stat=info)
     call check_allocation(info, 'rho_multipole_deriv           ')

  end if
  if (.not.allocated(dylm_dtheta_tab)) then
     allocate(dylm_dtheta_tab((l_pot_max+1)**2, n_centers_hartree_potential),stat=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

  end if
  if (.not.allocated(scaled_dylm_dphi_tab)) then
     allocate(scaled_dylm_dphi_tab((l_pot_max+1)**2,n_centers_hartree_potential),stat=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')

  end if



  !  initialize index_lm
  i_index = 0
  do i_l = 0, l_pot_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo




  call hartree_potential_real_coeff(index_lm, multipole_moments, &
       l_hartree_max_far_distance, n_centers_hartree_potential )

  if ( n_periodic > 0  ) then
     call update_outer_radius_l( outer_potential_radius, multipole_moments, & 
                                 multipole_radius_sq, l_hartree_max_far_distance, index_lm)

     do i_atom_2 = 1, n_atoms
        adap_outer_radius_sq(i_atom_2) = maxval(outer_potential_radius(:,i_atom_2))
     end do    
  else
     adap_outer_radius_sq = 1e8    
  end if



  if(n_periodic > 0)then
     call evaluate_hartree_recip_coef( index_lm, multipole_moments, &
                                       l_hartree_max_far_distance )
  end if


  ! initialize force components
  if(n_periodic > 0 .or. force_new_functional  )then
        hellman_feynman_forces = 0.0d0
  end if

  
  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

    ! First loop over the atoms: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the integration grid

    atom_of_splines = 0
    do i_center = 1, n_centers_hartree_potential, 1

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)

      if (( (n_periodic > 0) .or. force_new_functional) )then
        ! in this case, use multipole_spline to compute
        ! Hartree potential components ON each individual atomic nucleus, 
        ! before doing anything else.

        if (mod(current_spl_atom-1,n_tasks) == myid) then
          ! for this case rho_multipole for current_spl_atom is on myid

          do i_atom_2 = 1, n_occ_atoms, 1

            ! reinitialize the physical quantities of interest
            delta_v_hartree_aux = 0.d0

            v_hartree_gradient_temp = 0.d0


            ! get current integration point coordinate
            do i_coord = 1, 3, 1
              coord_current(i_coord) = coords(i_coord, i_atom_2)
            end do

            ! Tabulate distances and directions to all atoms -
            ! including relative positions on logarithmic and
            ! radial integration grids.

            call tab_single_atom_centered_coords_p0 &
                 ( current_center, coord_current,  &
                 dist_tab_sq,  &
                 dir_tab )

            if (i_atom_2.eq. current_center ) then
              dist_tab_sq = r_grid_min(species(i_atom_2))**2 + 1e-15
            end if



            ! VB: FIXME - At this point, must determine the maximum angular momentum
            !             used for atom i_atom_2 in current location!!!



            ! At each integration point, the Hartree potential components coming from
            ! different atoms are split into two groups:
            ! Any components coming from atoms close by are evaluated by explicit
            ! numerical splines; far away atoms are simply represented by
            ! an analytical long-distance multipole potential.
            if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom)) then
               ! begin with everything related to n_atoms_in ...

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

               ! for all inner atoms, we need ylm functions and their gradients explicitly
               call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)


                  call tab_single_gradient_ylm_p0 &
                       ( trigonom_tab, l_hartree, l_pot_max,  &
                       current_center,  &
                       ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
              

               ! Now loop over all inner atoms (those which are close enough
               ! to the current integration point so we need explicit numerical splines)
               !     partitioned
               !     hartree potentials need to be summed up
               !     according to Delley (eq. 12c)

               l_h_dim = (l_hartree(species_center(current_center))+1)**2

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

               if (i_atom_2.eq. current_center ) then
                  delta_v_hartree_multipole_component(1) =  delta_v_hartree_part_at_zero(current_spl_atom)

                  do i_lm = 2, l_h_dim


                     delta_v_hartree_multipole_component(i_lm) = 2* current_delta_v_hart_part_spl(i_lm,1,1) &
                          - current_delta_v_hart_part_spl(i_lm,1,2)

                  end do
               end if

               ! sum up the Hartree potential contribution from the present i_atom_in
               delta_v_hartree_aux = delta_v_hartree_aux + &
                    ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, &
                    ylm_tab, 1 )

               if ( n_periodic > 0  ) then
                  call far_distance_hartree_Fp_periodic_single_atom &
                       (current_spl_atom, i_center, &
                       dist_tab_in, l_hartree_max_far_distance, .true., .true., &
                       multipole_radius_sq(current_spl_atom),                      &
                       sqrt( adap_outer_radius_sq(current_spl_atom) )                 )
               end if

              
                  ! call spline vector derivative
                  ! dot priduct.
                  ! abs(V_radial_deriv) * dir(:,i_center_L)

                  call tab_single_radial_weights_v2 &
                       ( current_spl_atom, dist_tab_in, i_r, &
                       log_weight, radial_weight )

                  ! splines are now derivatives df/di where i is the grid point index
                  ! must convert to df/dr = df/di * di/dr

                  call spline_deriv_vector_v2 &
                       ( i_r_log, current_delta_v_hart_part_spl, &
                       (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid,  &
                       n_grid(species_center(current_center)), &
                       l_h_dim, &
                       delta_v_hartree_multipole_deriv)

                  delta_v_hartree_multipole_deriv(:) = &
                       delta_v_hartree_multipole_deriv(:) * log_weight


                  if(i_atom_2 == i_center) then

                     d_v_hartree_free_d_r =  0.d0

                     ! l=1 components

                     delta_v_hartree_multipole_deriv(2) = delta_v_hartree_deriv_l0_at_zero(1, i_atom_2)
                     delta_v_hartree_multipole_deriv(3) = delta_v_hartree_deriv_l0_at_zero(2, i_atom_2)
                     delta_v_hartree_multipole_deriv(4) = delta_v_hartree_deriv_l0_at_zero(3, i_atom_2)

                  else


                   if(test_Born_effective_charges_by_print_force) then 
!----------------d_v_hartree_free_d_r is non-zero for HF force----                      
                     d_v_hartree_free_d_r =  &
                          val_spline_deriv( i_r_log, &
                          free_pot_es_spl(1,1,species(current_spl_atom)),  &
                          n_grid(species(current_spl_atom))) * log_weight
                   else
!----------------d_v_hartree_free_d_r is zero for Born_effect_charge_HF----                      
                     d_v_hartree_free_d_r =  0.d0
                   endif 
 

                  end if

                  ! Remove the self interaction == nucleus interaction to itself.
                  v_hartree_gradient_temp = 0.0d0

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
                       v_hartree_gradient_temp)


            else if (dist_tab_sq.lt. adap_outer_radius_sq(current_spl_atom) )then

! VB: FIXME - is this correct for the cluster case? any dist_tab_sq should be the right one for the cluster case!

               ! Tabulate distances only for outer part
               dist_tab_out = sqrt(dist_tab_sq)

                  ! Now sum up the potential contributions from all far-field atoms
                  ! (analytical multipole potentials only ...)
                  call far_distance_hartree_Fp_periodic_single_atom &
                       (current_spl_atom, i_center, &
                       dist_tab_out, l_hartree_max_far_distance, .false. , .true., &
                       multipole_radius_sq(current_spl_atom),                         &
                       sqrt( adap_outer_radius_sq(current_spl_atom) )                  )

            end if ! multipole radius


            if (dist_tab_sq.lt. max(adap_outer_radius_sq(current_spl_atom), &
                 multipole_radius_sq(current_spl_atom)) )then


               if ( n_periodic > 0  ) then
                  ! Now sum up the far distance parts of the Hartree potential
                  call far_distance_real_hartree_potential_single_atom &
                       ( current_center, i_center, delta_v_hartree_aux, &
                       l_hartree_max_far_distance, coord_current )


                  
                     call far_distance_real_gradient_hartree_potential_single_atom &
                          ( current_spl_atom, i_center, dir_tab, &
                          v_hartree_gradient_temp, &
                          l_hartree_max_far_distance, &
                          dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )
                     ! hF forces = HF + Z (dV - Z/r^2 * dir)



               end if
            end if



               if(i_atom_2 == i_center )then
                  HF_forces =  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               else
                  HF_forces=  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               end if

               hellman_feynman_forces(1:3, i_atom_2 ) =  hellman_feynman_forces(1:3, i_atom_2) &
                    + HF_forces(1:3)


               v_hartree_gradient_temp = 0.d0


          end do ! loop over atoms
        end if
      end if ! (( (n_periodic > 0) .or. force_new_functional) .and. i_iter == 1 )

    end do  ! end loop over source atoms



  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for potential: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  if(n_periodic > 0  )then
     ! The Fourier part of the potential - in this loop, we no longer compute anything on the grid,
     ! but rather we compute the reciprocal-space part of the potential on each nucleus, and 
     ! add this to the energy of the nuclei in the Hartree potential of the electrons.

    
        ! The most expensive part of the calculations is evaluate_hartree_recip_gradient_coef_atom
        ! Therefore we take care that it is called only once for each atom and
        ! work is distributed over the tasks

        HF_temp = 0.d0

        do i_center = 1, n_atoms

           if (myid.eq.task_list(i_center)) then

              call  evaluate_hartree_recip_gradient_coef_atom( index_lm, multipole_moments, &
                   l_hartree_max_far_distance, i_center)

              do i_atom_2 = 1, n_occ_atoms, 1

                 ! get current integration point coordinate
                 coord_current(:) = coords(:, i_atom_2)

                 v_hartree_gradient_temp = 0.d0
                 call update_hartree_gradient_recip( coord_current,  v_hartree_gradient_temp, i_center)


                 HF_temp(1:3, i_atom_2 ) =  HF_temp(1:3, i_atom_2) &
                      + species_z(species(i_atom_2)) *  v_hartree_gradient_temp(1:3)
              end do

           end if

        enddo

        call sync_vector(HF_temp,3*n_atoms)



     do i_atom_2 = 1, n_occ_atoms, 1
        if (myid.eq.task_list(i_atom_2)) then

           ! reinitialize the physical quantities of interest
           delta_v_hartree_aux = 0.d0

           ! get current integration point coordinate
           coord_current(:) = coords(:, i_atom_2)

           delta_v_hartree_aux = 0.d0
           call update_hartree_potential_recip( coord_current, delta_v_hartree_aux )

           
              hellman_feynman_forces(1:3, i_atom_2 ) =  hellman_feynman_forces(1:3, i_atom_2) &
                   + HF_temp(1:3, i_atom_2 )

              !          write(use_unit,*) i_atom_2, HF_temp(1,i_atom_2)

        end if
     end do
     
    
  end if  ! n_periodic > 0

   !-------shanghui begin parallel------  
  call sync_vector(hellman_feynman_forces,3*n_atoms) 
   !-------shanghui end parallel------  





  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)


  if (allocated(current_rho_multipole_spl)) then
     deallocate(current_rho_multipole_spl)
  end if
  if (allocated(current_delta_v_hart_part_spl)) then
     deallocate(current_delta_v_hart_part_spl)
  end if

  if (allocated(delta_v_hartree_multipole_deriv)) then
     deallocate(delta_v_hartree_multipole_deriv)
  end if
  if (allocated(rho_multipole_deriv)) then
     deallocate(rho_multipole_deriv)
  end if
  if (allocated(dylm_dtheta_tab)) then
     deallocate(dylm_dtheta_tab)
  end if
  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate(scaled_dylm_dphi_tab)
  end if


end subroutine integrate_Born_effective_charges_Hellman_Feynman
!******
!------------------------------------------------------------------------------------------------------------
