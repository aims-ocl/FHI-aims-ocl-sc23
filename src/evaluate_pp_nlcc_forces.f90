! Hellman-Feynman forces on the pseudocores
! AT

subroutine evaluate_pp_nlcc_forces(temp_rho, partition, temp_rho_gradient, &
                                   temp_kinetic_density)

  use dimensions
  use runtime_choices
  use grids
  use spline
  use geometry
  use pseudodata
  use species_data
  use mpi_utilities
  use localorb_io
  use constants
  use synchronize_mpi
  use physics

      use xc

  implicit none
  ! imported variables
  
  ! input
  real*8, dimension(n_spin, n_full_points) :: temp_rho
  real*8, dimension(n_full_points) :: partition
  real*8, dimension(3,n_spin, n_full_points) :: temp_rho_gradient
  real*8, dimension(n_spin, n_full_points) :: temp_kinetic_density

! local variable

!  real*8, dimension(3, n_pp_atoms) :: pp_nlcc_forces 

  real*8 :: point_term
  real*8 :: atomic_term

  integer :: i_multipole
  integer :: i_coord
  integer :: i_spin

  integer :: i_index
  integer :: i_my_batch

  integer :: i_full_points

  real*8, dimension(3) :: direction
  real*8, dimension(3) :: coord_current
  real*8 :: distance_squared, dist

  real*8 :: en_density_c 

  real*8 :: rho_inc_partialcore, local_xc_derivs(2), en_density_xc, en_density_x(2)
  real*8 :: rho_gradient_inc_partialcore(3)
  real*8 :: d_partial_core_rho_d_R

  real*8 :: conversion

  integer :: i_atom, i_pp_atom, i_pp_species

  real*8 :: Z, q, i_r_log

  character*120 :: info_str
  character*40 :: func = 'evaluate_pp_hellman_feynman_forces_p0'

  ! output
!  real*8, dimension(3, n_atoms) :: pseudocore_forces 

  real*8 :: local_xc_gradient_deriv(3,n_spin)
  real*8 :: local_xc_tau_deriv(n_spin)

  conversion = hartree / bohr  

  ! first we add up the interaction with all nulcei

!  allocate(xc_gradient_deriv(3,n_spin, n_full_points))


  pp_nlcc_forces = 0.d0
  nlcc_forces = 0.d0
!  dummy = 0 

  i_full_points = 0


  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1
           if (partition_tab(i_full_points).gt.0.d0) then
              !     get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

              i_atom = batches(i_my_batch)%points(i_index)%index_atom


              do i_spin = 1, n_spin, 1
                rho_inc_partialcore = 0.d0
                rho_gradient_inc_partialcore = 0.d0
                do i_pp_atom = 1,n_pp_atoms


                 i_pp_species = pp_species(i_pp_atom)

                       
                 rho_inc_partialcore = &
                     temp_rho(i_spin,i_full_points) + partial_core_rho(i_full_points) !0.5

                 if(use_density_gradient) then 
                       rho_gradient_inc_partialcore(1:3) = &
                          temp_rho_gradient(1:3,i_spin,i_full_points) + & 
                          partial_core_rho_grad(1:3,i_full_points)
                 endif

                 call evaluate_xc  &
                      ( rho_inc_partialcore,   &
                        rho_gradient_inc_partialcore(1:3),  &
                        temp_kinetic_density, &
                        en_density_xc,   &
                        en_density_x, en_density_c,  &
                        local_xc_derivs(i_spin),  &
                        local_xc_gradient_deriv(1:3,i_spin), &
                        local_xc_tau_deriv(i_spin), &
                        .false. &
                      )


!                point_term = partition(i_full_points) * en_density_xc
                point_term = partition(i_full_points) * local_xc_derivs(i_spin)

              ! now find distance of pseudocore to integration point

                 distance_squared = 0.d0
                 do i_coord = 1,3,1

                   direction(i_coord) = coord_current(i_coord) - pp_coords(i_coord, i_pp_atom) 

                   distance_squared = distance_squared + &
                      (direction(i_coord))**2

                 enddo

                   dist = sqrt(distance_squared)


                 if(dist.gt.(localpot_outer_radius(i_pp_species))) then
                    d_partial_core_rho_d_R = 0.d0
                 elseif (dist.lt.(pp_r_grid_min(i_pp_species))) then
                    d_partial_core_rho_d_R = 0.d0
                 else

                    d_partial_core_rho_d_R = 0.d0
 
                    i_r_log = invert_log_grid(dist, &
                    pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))
 

                    d_partial_core_rho_d_R = val_spline( i_r_log, &
                         partial_core_dens_deriv_spl(1,1,i_pp_species), &
                         n_points_pp_fn(i_pp_species))


                  
                 endif
       
 
                 do i_coord = 1,3,1
                 pp_nlcc_forces(i_coord, i_pp_atom) = pp_nlcc_forces(i_coord, i_pp_atom) +  &
                         point_term*partial_core_rho_grad(i_coord,i_full_points)



!                 nlcc_forces(i_coord, i_atom) = nlcc_forces(i_coord, i_atom) + &
!                         point_term*partial_core_rho_grad(i_coord,i_full_points)
! + &
!                         partition(i_full_points) * local_xc_derivs_2(i_spin)*partial_core_rho_grad(i_coord,i_full_points)
!                         point_term*partial_core_rho_grad(i_coord,i_full_points)




                 enddo

!                 if(use_density_gradient) then
!                 do i_coord = 1,3,1
!                   pp_nlcc_forces(i_coord, i_pp_atom) = pp_nlcc_forces(i_coord, i_pp_atom) +  &
!                      partition(i_full_points) 

!                 enddo
!                 endif

              enddo !i_pp_atom
             enddo !i_spin

          endif
        enddo

  enddo

  call sync_pp_charge_forces(pp_nlcc_forces)
  call sync_nlcc_forces(nlcc_forces)


end subroutine evaluate_pp_nlcc_forces
