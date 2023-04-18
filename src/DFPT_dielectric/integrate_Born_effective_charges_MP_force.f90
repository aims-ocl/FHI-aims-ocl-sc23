!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/sum_up_whole_potential
!  NAME
!   Here I remove the variables that I do not want. Only loop over grid to get MP-part.
!  SYNOPSIS

subroutine integrate_Born_effective_charges_MP_force &
     (multipole_moments, &
     partition_tab_std, rho_std, rho_multipole, free_rho_superpos_std, &
     v_hartree_gradient_multipole, &      
     multipole_radius_sq, &
     l_hartree_max_far_distance, &
     outer_potential_radius, &
     multipole_forces )

!  PURPOSE
!    The subroutine get MP part in Born effective charge.
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
  use pbc_lists
  use load_balancing
  ! rho_multipole from hartree_potential_storage conflicts with the local variable here,
  ! so just only use get_rho_multipole_spl() from this module.
  ! Maybe the local variable name should be changed ...
  use hartree_potential_storage, only : get_rho_multipole_spl

  implicit none

!  ARGUMENTS


  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments

  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(n_full_points) :: rho_std
  real*8, dimension(n_full_points), intent(inout)  :: rho_multipole 
  real*8, target, dimension(n_full_points)        :: free_rho_superpos_std
  real*8, dimension(3, n_atoms, n_full_points), intent(out)  :: v_hartree_gradient_multipole 


  real*8, dimension(n_atoms)              :: multipole_radius_sq
  integer, dimension( n_atoms)            :: l_hartree_max_far_distance
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius
  real*8, dimension(3, n_atoms)           :: multipole_forces

!  INPUTS
! o multipole_moments -- mutlpole moments of the Hartree potential
! o partition_tab_std -- values of partition function
! o rho_std -- electron density
! o multipole_radius_sq -- outer radius of multipole components
! o l_hartree_max_far_distance -- maximum l-components of for the far distance Hartree potential (periodic systems)
!
!  OUTPUT
! o multipole_forces -- mutipole forces = the force components of multipole energy correction
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

  ! work arrays
  real*8, allocatable :: delta_v_hartree(:)
  real*8, allocatable :: v_hartree_free(:)

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


  real*8 :: total_rho

  real*8 :: force_component(n_centers_hartree_potential)

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

  real*8 :: d_v_hartree_free_d_r
  real*8 :: d_rho_free_d_r

  real*8 rho_aux
  real*8 pot_aux
  real*8 rho_multipole_aux
  real*8 delta_v_hartree_aux
  real*8 v_hartree_free_aux

  ! In periodic systems (only), this will be the correct potential zero for the electrostatic potential
  ! The average electrostatic potential from the previous iteration will be saved, as the 
  !     current "chemical potential" (the Fermi level) was shifted by the average electrostatic
  !     potential from the previous iteration.
  real*8 :: average_delta_v_hartree_real
  real*8, save :: previous_average_delta_v_hartree_real = 0.d0

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
  integer :: i_iter, n_iter

  integer :: i_lat
  integer :: i_full_points

  !  external functions
  real*8, external :: ddot
  character*100 :: info_str


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
  real*8, pointer :: free_rho_superpos(:)

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  if(use_batch_permutation > 0) then
    call localorb_info("Integrating Born_effective_charges_MP matrix with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Integrating Born_effective_charges_MP matrix.", use_unit,'(2X,A)', OL_norm )
  endif

    hartree_force_l_add = 1

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    n_full_points_work = batch_perm(n_bp)%n_full_points

    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(rho(n_full_points_work))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    allocate(free_rho_superpos(n_full_points_work))
    call permute_point_array(n_bp,1,free_rho_superpos_std,free_rho_superpos)

  else

    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    free_rho_superpos => free_rho_superpos_std

  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------


  allocate(delta_v_hartree(n_full_points_work),stat=info)
  call check_allocation(info, 'delta_v_hartree')

  allocate(v_hartree_free(n_full_points_work),stat=info)
  call check_allocation(info, 'v_hartree_free')

  

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

  if ( n_periodic > 0 ) then
     call update_outer_radius_l( outer_potential_radius, multipole_moments, & 
                                 multipole_radius_sq, l_hartree_max_far_distance, index_lm)

     do i_atom_2 = 1, n_atoms
        adap_outer_radius_sq(i_atom_2) = maxval(outer_potential_radius(:,i_atom_2))
     end do    
  else
     adap_outer_radius_sq = 1e8    
  end if



  ! initialize energies and force components
  multipole_forces(:,:) = 0.d0

 

  ! There are two separate loops over integration points:
  ! 1) A loop over the integration grid to tabulate the multipole components and densities at each
  !    grid point one by one
  ! 2) A second loop over the grid to integrate over products rho_multipole(r)*v_multipole(r)

  ! initialize the physical quantities that are tabulated over the entire integration grid
  rho_multipole    = 0.d0
  v_hartree_gradient_multipole = 0.0d0
  delta_v_hartree  = 0.d0
  v_hartree_free   = 0.d0

  ! The following averages are only needed in periodic systems
  average_delta_v_hartree_real = 0.d0
  
  n_iter = 2
  
  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_iter = 1, n_iter

    ! loop over the grid: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the integration grid

    atom_of_splines = 0
    do i_center = 1, n_centers_hartree_potential, 1

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)


      !
      ! Now follows the real work: Summing the multipole potential, density, and their
      ! derivatives on the integration grid
      !

      ! Reset grid counter for current Hartree potential center
      i_full_points = 0

      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          ! Only execute if partition_tab is .gt. zero, else
          ! we can run into 1/0 !!!
          if (partition_tab(i_full_points).gt.0.d0) then

            if(i_iter == 2) then
              v_hartree_gradient_temp = 0
              rho_multipole_gradient_temp = 0
            endif

            ! get current integration point coordinate
            coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)

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
            if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom) ) then

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

              if (i_iter == 2) then
                ! need radial derivatives di/dr for the benefit of
                ! spline derivative evaluation later

                call tab_single_radial_weights_v2 &
                     ( current_spl_atom, dist_tab_in, i_r, &
                     log_weight, radial_weight )
              end if

              ! for an inner atom we need ylm functions and their gradients explicitly

              call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

              if (i_iter == 2) then
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

              if(i_iter == 1) then
                ! sum up the Hartree potential contribution from the present inner atom
                delta_v_hartree_aux = &
                       ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab, 1 )

                delta_v_hartree(i_full_points) = &
                       delta_v_hartree(i_full_points) + delta_v_hartree_aux
              endif

              ! Obtain spline-interpolated values of the multipole density,
              ! this time splined on the radial integration grid
              call spline_vector_v2 &
                   ( i_r+1, current_rho_multipole_spl, &
                   (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                   n_radial(species_center(current_center))+2,  &
                   l_h_dim, &
                   rho_multipole_component)

              if(i_iter == 1) then
                ! sum up the multipole density contribution from the present inner atom
                rho_multipole_aux = &
                       ddot ( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)

                rho_multipole(i_full_points) = rho_multipole(i_full_points) + rho_multipole_aux


              endif

              ! If needed, obtain multipole force correction for present inner atom
              if (i_iter == 2) then
                ! get all needed derivatives for forces from splines

                call spline_deriv_vector_v2 &
                     ( i_r_log, &
                     current_delta_v_hart_part_spl, &
                     (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid,  &
                     n_grid(species_center(current_center)), &
                     l_h_dim, &
                     delta_v_hartree_multipole_deriv)

                call spline_deriv_vector_v2 &
                     ( i_r+1,  &
                     current_rho_multipole_spl, &
                     (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                     n_radial(species_center(current_center))+2,  &
                     l_h_dim, &
                     rho_multipole_deriv)

                ! splines are now derivatives df/di where i is the grid point index
                ! must convert to df/dr = df/di * di/dr
                rho_multipole_deriv(1:l_h_dim) = rho_multipole_deriv(1:l_h_dim) * radial_weight
                                             
                delta_v_hartree_multipole_deriv(1:l_h_dim) = &
                     delta_v_hartree_multipole_deriv(1:l_h_dim) * log_weight

                ! gradients of the free atom density and potential is set to 0.0d0.
!                d_v_hartree_free_d_r = 0.0d0 
!                d_rho_free_d_r = 0.0d0
                  d_v_hartree_free_d_r =  &
                        val_spline_deriv( i_r_log, &
                        free_pot_es_spl(1,1,species(current_spl_atom)),  &
                        n_grid(species(current_spl_atom))) * log_weight


                d_rho_free_d_r = pi4_inv *  &
                     val_spline_deriv(i_r_log,  &
                     free_rho_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom)))  * log_weight

                ! obtain gradients for present atom ...
                call evaluate_v_hartree_and_rho_multipole_gradient &
                     ( dist_tab_in, dir_tab_in,  &
                     trigonom_tab,  &
                     ylm_tab,  &
                     dylm_dtheta_tab,  &
                     scaled_dylm_dphi_tab,  &
                     delta_v_hartree_multipole_component,  &
                     delta_v_hartree_multipole_deriv, &
                     rho_multipole_component, &
                     rho_multipole_deriv,  &
                     d_v_hartree_free_d_r, d_rho_free_d_r, &
                     l_h_dim, &
                     v_hartree_gradient_temp, &
                     rho_multipole_gradient_temp)

              end if ! iter=2


              ! far-distance treatment for then center i_center for the periodic system
              if ( n_periodic > 0 ) then

                !  VB: FIXME!!!!! 
                !      I kept the original call here because I was not sure I was doing the completely
                !      right thing. We should cut the periodic part down to l_atom_max as well below,
                !      but that requires a change to the periodic infrastructure first!
                !
                call far_distance_hartree_Fp_periodic_single_atom &
                      (current_spl_atom, i_center, &
                      dist_tab_in, l_hartree_max_far_distance, .true., .true.,  &
                      multipole_radius_sq(current_spl_atom),                       &
                      sqrt( adap_outer_radius_sq(current_spl_atom) )                )

              end if


            else if (dist_tab_sq .lt. adap_outer_radius_sq(current_spl_atom)) then
              ! the current center is in the far-distance part-----------------

              ! Tabulate distance for the outer part
              dist_tab_out = sqrt(dist_tab_sq)

              if (i_iter == 2) then
                dir_tab_out = dir_tab / dist_tab_out
              end if



                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                ! nuclear z/r part
                !  d_v_hartree_free_d_r =  &
                !       - 0*species_z(species(current_spl_atom))/dist_tab_sq

                ! d_v_hartree_free_d_r = 0.d0

                ! v_hartree_gradient(:,i_center,i_full_points) =  &
                !      dir_tab_out(:) * d_v_hartree_free_d_r

                call far_distance_hartree_Fp_periodic_single_atom &
                     (current_spl_atom, i_center, &
                     dist_tab_out, l_hartree_max_far_distance, .false.,.true.,  &
                     multipole_radius_sq(current_spl_atom),                        &
                     sqrt( adap_outer_radius_sq(current_spl_atom) )                 )

                ! Note: 'far_distance_real_hartree_potential_single_atom' is called in the periodic
                ! systems later, because then we have Fp in inner and outer multipole radius.

            end if  ! end if for separating current_center either far-distance or near part


            ! Now sum up the far distance part of the Hartree potential if Ewald's
            ! decomposition is performed. In the opposite case, this was already done.
            if ( n_periodic > 0  ) then


              if (dist_tab_sq.lt. &
                   max(adap_outer_radius_sq(current_spl_atom),multipole_radius_sq(current_spl_atom))) then

                ! Far distance analytic real-space part of the potential

                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                if(i_iter == 1) then
                  call far_distance_real_hartree_potential_single_atom &
                          ( current_center, i_center, delta_v_hartree(i_full_points), &
                          l_hartree_max_far_distance, coord_current )

                endif

                if(i_iter == 2)then

                  ! Far distace analytic real-space part of the potential gradient

                  call far_distance_real_gradient_hartree_potential_single_atom &
                         ( current_spl_atom, i_center, dir_tab, &
                         v_hartree_gradient_temp, &
                         l_hartree_max_far_distance, &
                         dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )

                    

                end if
              end if

            end if ! end periodic system



            ! Force calculations in second iteration

            if (i_iter == 2) then

              total_rho = 0.d0
              do i_spin = 1, n_spin, 1
                total_rho = total_rho + rho(i_full_points)
              enddo
              

             ! if(i_center.eq.1.and.i_full_points.le.2) then 
             ! write(use_unit,*) 'shanghui test:', i_iter,  rho(i_full_points),rho_multipole(i_full_points), &
             !                              pi4_inv * free_rho_superpos(i_full_points), &
             !                              v_hartree_gradient_temp(1)
             ! endif 


                if (.not.empty(center_to_atom(current_center))) then

                  do i_coord = 1, 3, 1


                    multipole_forces(i_coord,center_to_atom(current_center)) = &
                            multipole_forces(i_coord,center_to_atom(current_center)) &
                            + (total_rho -  rho_multipole(i_full_points)) &
                            * v_hartree_gradient_temp(i_coord) &
                            * partition_tab(i_full_points)

                    v_hartree_gradient_multipole(i_coord,center_to_atom(current_center), & 
                                                 i_full_points) = & 
                    v_hartree_gradient_multipole(i_coord,center_to_atom(current_center), & 
                                                 i_full_points) + & 
                    v_hartree_gradient_temp(i_coord)

                  end do
                endif ! end cycle over ghost atoms


            end if ! if (i_iter == 2)

          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    end do  ! end loop over source atoms

    if(i_iter == 1)then
      ! next, we integrate the actual energy and force quantities
      i_full_points = 0
      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          if (partition_tab(i_full_points).gt.0.d0) then

            rho_multipole(i_full_points) = &
                 rho_multipole(i_full_points) + &
                 pi4_inv * free_rho_superpos(i_full_points)

          end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    endif ! if(i_iter == 1)

  enddo ! Loop over n_iter




  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for potential: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  if(get_batch_weights) then
    call set_batch_weights(n_bp,batch_times)
    deallocate(batch_times)
  endif


   !-------shanghui begin parallel------  
  call sync_vector(multipole_forces,3*n_atoms)
   !-------shanghui end parallel------  





  if(use_batch_permutation > 0) then


    deallocate(rho)

  endif

  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  if(allocated(delta_v_hartree)) deallocate(delta_v_hartree)
  if(allocated(v_hartree_free))  deallocate(v_hartree_free)

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


end subroutine integrate_Born_effective_charges_MP_force
!******
!------------------------------------------------------------------------------------------------------------
