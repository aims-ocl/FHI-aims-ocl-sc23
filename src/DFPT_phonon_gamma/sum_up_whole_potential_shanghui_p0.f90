!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/sum_up_whole_potential_shanghui_p0
!  NAME
!   sum_up_whole_potential_shanghui_p0
!  SYNOPSIS

subroutine sum_up_whole_potential_shanghui_p0 &
     ( delta_v_hartree_part_at_zero,  delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, &
     partition_tab_std, delta_rho_std, &
     delta_potential_std,  &
     forces_on,multipole_radius_sq, &
     l_hartree_max_far_distance,  &
     outer_potential_radius, &
     hellman_feynman_term)

!  PURPOSE
!    The subroutine sums up the hartree potential from the multipole components
!    which have been calculated before hand. The actual summation is done here.
!
!    Using  rho_multiple(*****delta_rho******)  calcualted in previous subroutine
!    update_hartree_potential_p2_shanghui.f90  by get_rho_multipole_spl
!    and then output (******delta_potential*****)
!
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
  use load_balancing
  ! rho_multipole from hartree_potential_storage conflicts with the local variable here,
  ! so just only use get_rho_multipole_spl() from this module.
  ! Maybe the local variable name should be changed ...
  use hartree_potential_storage, only : get_rho_multipole_spl

  implicit none

!  ARGUMENTS


  real*8, dimension(n_atoms)                     :: delta_v_hartree_part_at_zero
  real*8, dimension(3, n_atoms)                  :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments

  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(n_full_points) :: delta_rho_std
  real*8, target, dimension(n_full_points)        :: delta_potential_std


  logical :: forces_on
  real*8, dimension(n_atoms)              :: multipole_radius_sq
  integer, dimension( n_atoms)            :: l_hartree_max_far_distance
  real*8, dimension(3, n_atoms)           :: hellman_feynman_term
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius

!  INPUTS
! o delta_v_hartree_part_at_zero -- Hartree potential at origin of the atoms from different multipole components
! o delta_v_hartree_deriv_l0_at_zero -- Derivative of Hartree potential at origin of the atoms
! o multipole_moments -- mutlpole moments of the Hartree potential
! o partition_tab_std -- values of partition function
! o delta_rho_std -- delta_electron density
! o forces_on -- are the forces calculated in this round or not?
! o multipole_radius_sq -- outer radius of multipole components
! o l_hartree_max_far_distance -- maximum l-components of for the far distance Hartree potential (periodic systems)
!
!  OUTPUT
! o delta_potential_std -- delta_ Hartree potential
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
  real*8, allocatable :: rho_multipole(:)

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

  integer :: i_lat
  integer :: i_full_points

  !  external functions
  real*8, external :: ddot
  character*100 :: info_str
  real*8:: HF_forces(3)

  real*8 :: HF_temp(3,n_atoms)

  integer :: mpierr
  integer :: info


  real*8:: out_vacuum_potential_z, dip_gradient, dip_origin, dip_lenght, dip_coord_current
  integer:: i_z

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

  if(use_batch_permutation > 0) then
    call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  endif

  if (forces_on) then
    hartree_force_l_add = 1
  else
    hartree_force_l_add = 0
  endif

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
!  if(use_batch_permutation > 0) then
!
!    n_my_batches_work = batch_perm(n_bp)%n_my_batches
!    n_full_points_work = batch_perm(n_bp)%n_full_points
!
!    batches_work => batch_perm(n_bp)%batches
!    partition_tab => batch_perm(n_bp)%partition_tab
!
!    allocate(rho(n_full_points_work))
!    call permute_point_array(n_bp,delta_rho_std,rho)
!!
!
!
!    allocate(potential(n_full_points_work))

!  else

    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std
    rho => delta_rho_std
    potential => delta_potential_std

!  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !-----------------------------------------------------------------------------


  allocate(delta_v_hartree(n_full_points_work),stat=info)
  call check_allocation(info, 'delta_v_hartree')

  allocate(v_hartree_free(n_full_points_work),stat=info)
  call check_allocation(info, 'v_hartree_free')

  allocate(rho_multipole(n_full_points_work),stat=info)
  call check_allocation(info, 'rho_multipole')
  

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

  if (forces_on) then
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
  end if

  if(n_periodic > 0 .or. force_new_functional  )then
     if(forces_on) then
        hellman_feynman_term = 0.0d0
     end if
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

  if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
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
  else if (use_hartree_non_periodic_ewald) then
     call calculate_extended_hartree_non_periodic_ewald( l_hartree_max_far_distance, &
                                           multipole_radius_sq, adap_outer_radius_sq  )
  end if

  ! Initialize potential (which is an output variable and runs over all grid points
  potential   = 0.d0

  ! initialize energies and force components
  hartree_multipole_error = 0.d0



  
  dip_gradient = 0
  dip_origin = 0
  dip_lenght = 0

  if(n_periodic > 0  )then
     
     if(use_dipole_correction .or. calculate_work_function)then
        if(dipole_correction_method=='potential') then
             write(info_str,*) 'Calculating dipole correction via potential gradient'
             call evaluate_dipole_correction & 
            ( previous_average_delta_v_hartree_real, dip_gradient, dip_origin, dip_lenght,.true. )
        elseif(dipole_correction_method=='dipole') then
             write(info_str,*) 'Calculating dipole correction via slab dipole moment'
           call evaluate_dipole_correction_from_dipole &
           (previous_average_delta_v_hartree_real, dip_gradient,dip_origin,dip_lenght,.true.)
        endif
     end if
  end if
 

  ! There are two separate loops over integration points:
  ! 1) A loop over the integration grid to tabulate the multipole components and densities at each
  !    grid point one by one
  ! 2) A second loop over the grid to integrate over products rho_multipole(r)*v_multipole(r)

  ! initialize the physical quantities that are tabulated over the entire integration grid
  rho_multipole    = 0.d0
  delta_v_hartree  = 0.d0
  v_hartree_free   = 0.d0

  ! The following averages are only needed in periodic systems
  average_delta_v_hartree_real = 0.d0


  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()



    ! First loop over the grid: We run over the Hartree potential center by center, adding the Hartree
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

            if(forces_on)then
              v_hartree_gradient_temp = 0.d0
            end if



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

              !write(use_unit,*) 'stop in sum_up' 
              !stop !shanghui

               call tab_single_atom_centered_coords_radial_log_p0 &
                    ( current_center, dist_tab_sq, dir_tab,  &
                    dist_tab_in, i_r, i_r_log, dir_tab_in )

               ! for all inner atoms, we need ylm functions and their gradients explicitly
               call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

               if (forces_on) then
                  call tab_single_gradient_ylm_p0 &
                       ( trigonom_tab, l_hartree, l_pot_max,  &
                       current_center,  &
                       ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
               else
                  call tab_single_wave_ylm_p0 &
                       ( current_center,  &
                       trigonom_tab, l_hartree,  &
                       l_pot_max, ylm_tab)
               end if

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

               if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
                  call far_distance_hartree_Fp_periodic_single_atom &
                       (current_spl_atom, i_center, &
                       dist_tab_in, l_hartree_max_far_distance, .true., forces_on, &
                       multipole_radius_sq(current_spl_atom),                      &
                       sqrt( adap_outer_radius_sq(current_spl_atom) )                 )
               end if

               if(forces_on)then
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

                  ! obtain gradients of the free atom density and potential explicitly

                  ! radial derivative of hartree potential of free atom # i_atom (without core potential!!!)
                  ! subtract core potential => - ( - Z/r ) = + Z/r , d/dr => - Z/r^2

                  if(i_atom_2 == i_center) then

                     d_v_hartree_free_d_r =  0.d0

                     ! l=1 components

                     delta_v_hartree_multipole_deriv(2) = delta_v_hartree_deriv_l0_at_zero(1, i_atom_2)
                     delta_v_hartree_multipole_deriv(3) = delta_v_hartree_deriv_l0_at_zero(2, i_atom_2)
                     delta_v_hartree_multipole_deriv(4) = delta_v_hartree_deriv_l0_at_zero(3, i_atom_2)

                  else
                  d_v_hartree_free_d_r =  &
                          val_spline_deriv( i_r_log, &
                          free_pot_es_spl(1,1,species(current_spl_atom)),  &
                          n_grid(species(current_spl_atom))) * log_weight


                  end if

                  ! Remove the self interaction == nucleus interaction to itself.
                  v_hartree_gradient_temp = 0.0d0

                  ! obtain gradients for present atom ...

                  !--------shanghui add for new_hellman_feynman_hessian---------------
                   d_v_hartree_free_d_r = 0.0d0
                  !--------shanghui end add for new_hellman_feynman_hessian---------

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

               end if !(forces_on)

            else if (dist_tab_sq.lt. adap_outer_radius_sq(current_spl_atom) )then

! VB: FIXME - is this correct for the cluster case? any dist_tab_sq should be the right one for the cluster case!

               ! Tabulate distances only for outer part
               dist_tab_out = sqrt(dist_tab_sq)

               if ( n_periodic ==0 .and. .not. use_hartree_non_periodic_ewald ) then


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

                  ! ... and compute force contribution of a far-distance atom
                  if (forces_on) then

                     ! With the neutralized energy functional. there is no separate nuclear component to the 
                     ! gradient of the potential here.
                     d_v_hartree_free_d_r = 0.d0
!                     v_hartree_gradient_temp =  dir_tab_out(:) * d_v_hartree_free_d_r
                     v_hartree_gradient_temp = 0.d0

                     ! add multipole potential gradients to free-atom gradient
                     ! The periodic case is calculated later because it uses both outer and inner atoms.
                     call far_distance_real_gradient_hartree_potential_single_atom_p2 &
                          ( current_spl_atom, dir_tab, &
                          v_hartree_gradient_temp, &
                          l_hartree_max_far_distance(current_spl_atom) )

                  end if

               else
                  ! Now sum up the potential contributions from all far-field atoms
                  ! (analytical multipole potentials only ...)
                  call far_distance_hartree_Fp_periodic_single_atom &
                       (current_spl_atom, i_center, &
                       dist_tab_out, l_hartree_max_far_distance, .false. , forces_on, &
                       multipole_radius_sq(current_spl_atom),                         &
                       sqrt( adap_outer_radius_sq(current_spl_atom) )                  )

               end if
            end if ! multipole radius


            if (dist_tab_sq.lt. max(adap_outer_radius_sq(current_spl_atom), &
                 multipole_radius_sq(current_spl_atom)) )then


               if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
                  ! Now sum up the far distance parts of the Hartree potential
                  call far_distance_real_hartree_potential_single_atom &
                       ( current_center, i_center, delta_v_hartree_aux, &
                       l_hartree_max_far_distance, coord_current )


                  if(forces_on)then
                     call far_distance_real_gradient_hartree_potential_single_atom &
                          ( current_spl_atom, i_center, dir_tab, &
                          v_hartree_gradient_temp, &
                          l_hartree_max_far_distance, &
                          dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )
                     ! hF forces = HF + Z (dV - Z/r^2 * dir)

                     if(use_dipole_correction .and. i_atom_2 == i_center )then
                       v_hartree_gradient_temp(3) = v_hartree_gradient_temp(3) - dip_gradient
                     end if

                  end if
               end if
            end if


            if(forces_on)then
               if(i_atom_2 == i_center )then
                  HF_forces =  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               else
                  HF_forces=  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               end if

               hellman_feynman_term(1:3, i_atom_2 ) =  hellman_feynman_term(1:3, i_atom_2) &
                    - HF_forces(1:3) ! shanghui change to '-', we need hessian. 

               v_hartree_gradient_temp = 0.d0

            end if
          end do ! loop over atoms
        end if  ! parallel over i_center_spl_atom
      end if ! (( (n_periodic > 0) .or. force_new_functional) .and. i_iter == 1 )

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

              ! for an inner atom we need ylm functions and their gradients explicitly

              call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)

              call tab_single_wave_ylm_p2 &
                     ( trigonom_tab, l_atom_max,  &
                     l_pot_max, ylm_tab)


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

                ! sum up the Hartree potential contribution from the present inner atom
                delta_v_hartree_aux = &
                       ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab, 1 )

                delta_v_hartree(i_full_points) = &
                       delta_v_hartree(i_full_points) + delta_v_hartree_aux

              ! Obtain spline-interpolated values of the multipole density,
              ! this time splined on the radial integration grid
              call spline_vector_v2 &
                   ( i_r+1, current_rho_multipole_spl, &
                   (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                   n_radial(species_center(current_center))+2,  &
                   l_h_dim, &
                   rho_multipole_component)

                ! sum up the multipole density contribution from the present inner atom
                rho_multipole_aux = &
                       ddot ( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)

                rho_multipole(i_full_points) = rho_multipole(i_full_points) + rho_multipole_aux

                if (n_periodic.eq.0 .and. (.not. force_new_functional)) then
                  ! non-periodic calculation; need electron-only Hartree potential later on.
                  if (.not.empty(  current_spl_atom)) then
                    v_hartree_free(i_full_points) = v_hartree_free(i_full_points) + &
                         species_z(species_center(current_center)) / &
                         dist_tab_in
                  endif
                end if


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



              ! Now sum up the potential contributions from a far-field atom
              ! (analytical multipole potentials only ...)
              if ( n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald ) then

                ! if not periodic, need electronic v_hartree_free part for total energy
                ! FIXME - THIS PART SHOULD SIMPLY GO ENTIRELY ONE DAY, WHEN THE FORMALISM
                ! IS UNIFIED WITH THE POERIODIC FORMALISM

                  if (.not.(empty(current_spl_atom))) then

                    if(.not. force_new_functional)then
                      v_hartree_free_aux = &
                               species_z(species(current_spl_atom)) / dist_tab_out
                    else
                      v_hartree_free_aux = 0.d0
                    end if

                    v_hartree_free(i_full_points) = &
                            v_hartree_free(i_full_points) + v_hartree_free_aux
                  end if

                ! Recursively tabulate the radial behavior of all multipole components
                ! These radial functions are a private variable in module
                ! hartree_potential_real.f90, and are reused there later in
                ! subroutine far_distance_hartree_potential_real_single_atom
                call far_distance_hartree_Fp_cluster_single_atom_p2 &
                     ( dist_tab_out, &
                       l_atom_max, forces_on )

                call far_distance_real_hartree_potential_single_atom_p2 &
                     ( i_center, delta_v_hartree(i_full_points), &
                       l_atom_max, coord_current )

                ! ... and compute force contribution of a far-distance atom

              else ! Periodic system


                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                ! nuclear z/r part
                !  d_v_hartree_free_d_r =  &
                !       - 0*species_z(species(current_spl_atom))/dist_tab_sq

                ! d_v_hartree_free_d_r = 0.d0

                ! v_hartree_gradient(:,i_center,i_full_points) =  &
                !      dir_tab_out(:) * d_v_hartree_free_d_r

                call far_distance_hartree_Fp_periodic_single_atom &
                     (current_spl_atom, i_center, &
                     dist_tab_out, l_hartree_max_far_distance, .false.,forces_on,  &
                     multipole_radius_sq(current_spl_atom),                        &
                     sqrt( adap_outer_radius_sq(current_spl_atom) )                 )

                ! Note: 'far_distance_real_hartree_potential_single_atom' is called in the periodic
                ! systems later, because then we have Fp in inner and outer multipole radius.

              end if ! if peridic system.
            end if  ! end if for separating current_center either far-distance or near part


            ! Now sum up the far distance part of the Hartree potential if Ewald's
            ! decomposition is performed. In the opposite case, this was already done.
            if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then


              if (dist_tab_sq.lt. &
                   max(adap_outer_radius_sq(current_spl_atom),multipole_radius_sq(current_spl_atom))) then

                ! Far distance analytic real-space part of the potential

                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                  call far_distance_real_hartree_potential_single_atom &
                          ( current_center, i_center, delta_v_hartree(i_full_points), &
                          l_hartree_max_far_distance, coord_current )


              end if

            end if ! end periodic system




          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    end do  ! end loop over source atoms

      ! next, we integrate the actual energy and force quantities
      i_full_points = 0
      do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

          if (partition_tab(i_full_points).gt.0.d0) then

            ! the reciprocal space contribution needs to be treated separately
            ! from the centers above
            if (n_periodic.gt.0) then
              ! get current integration point coordinate
              coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)


              ! before the reciprocal-space component is added, average over the real-space part of 
              ! delta_v_hartree
              average_delta_v_hartree_real = average_delta_v_hartree_real + &
              delta_v_hartree(i_full_points) * partition_tab(i_full_points)

              call update_hartree_potential_recip &
                    ( coord_current, delta_v_hartree(i_full_points) )

              if( use_dipole_correction)then

                i_m = int(floor((coord_current(3) - vacuum_z_level)/dip_lenght))
                dip_coord_current = coord_current(3) - i_m * dip_lenght
                  
                if( dip_coord_current <  vacuum_z_level)then
                  dip_coord_current = dip_coord_current + dip_lenght
                end if
      
                delta_v_hartree(i_full_points)  =  delta_v_hartree(i_full_points) &
                       -  (dip_coord_current-dip_origin) * dip_gradient

              end if

            else if ( use_hartree_non_periodic_ewald ) then
              call interpolate_extended_hartree_non_periodic_ewald(              &
                              batches_work(i_batch) % points(i_index) % coords,  &
                              delta_v_hartree(i_full_points)                      )
            end if ! n_periodic > 0


            rho_multipole(i_full_points) = &
                 rho_multipole(i_full_points) 

            potential(i_full_points) = &
                 delta_v_hartree(i_full_points) 
                 ! shanghui changed here for p2_shanghui
                 !+ free_hartree_superpos(i_full_points)


              total_rho = rho( i_full_points)

            !            Add contribution to difference hartree-energy
            !            Use correction to the Hartree potential which reduces the error
            !            from linear to quadratic order in the multipole expansion:
            !
            !            E_H = int [V_H * (rho_full - 1/2 rho_multipole)]
            !
            !            Dunlap et al JCP 71, 3396 (1979) or Bastug et al, CPL 211, 119 (1993)
            !
            !            But in essence this is trivial as it only means calc. E_H = int(rho_mult*v_mult)

            


            hartree_multipole_error = &
                 hartree_multipole_error  + &
                 ( total_rho - rho_multipole(i_full_points) )**2 &
                 * partition_tab(i_full_points)


            ! finally, after all is said and done, we add a possible external
            ! embedding potential to the total electrostatic potential
            !------shanghui says: we do not need embedding potential here
            !!if (use_embedding_potential) then
            !!  if (full_embedding) then
            !!    potential(i_full_points) = potential(i_full_points) +  &
            !!           pot_ion_embed(i_full_points)
            !!  end if
            !!  en_density_embed = en_density_embed + partition_tab(i_full_points) * &
            !!        pot_ion_embed(i_full_points) * total_rho
            !!end if

          end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches


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


  if(n_periodic > 0  )then
     ! The Fourier part of the potential - in this loop, we no longer compute anything on the grid,
     ! but rather we compute the reciprocal-space part of the potential on each nucleus, and 
     ! add this to the energy of the nuclei in the Hartree potential of the electrons.

     if(forces_on)then

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

     endif

     do i_atom_2 = 1, n_occ_atoms, 1
        if (myid.eq.task_list(i_atom_2)) then

           ! reinitialize the physical quantities of interest
           delta_v_hartree_aux = 0.d0

           ! get current integration point coordinate
           coord_current(:) = coords(:, i_atom_2)

           delta_v_hartree_aux = 0.d0
           call update_hartree_potential_recip( coord_current, delta_v_hartree_aux )

           if( use_dipole_correction)then

              i_m = int(floor((coord_current(3) - vacuum_z_level)/dip_lenght))
              coord_current(3) = coord_current(3) - i_m * dip_lenght
              
              if( coord_current(3) <  vacuum_z_level)then
                 coord_current(3) = coord_current(3) + dip_lenght
              end if
              
              delta_v_hartree_aux = delta_v_hartree_aux -  (coord_current(3)-dip_origin) * dip_gradient

           end if



           if(forces_on)then
              hellman_feynman_term(1:3, i_atom_2 ) =  hellman_feynman_term(1:3, i_atom_2) &
                   - HF_temp(1:3, i_atom_2 ) ! shanghui change to '-', we need hessian
          
              !          write(use_unit,*) i_atom_2, HF_temp(1,i_atom_2)
           end if

        end if
     end do
     
  else if (use_hartree_non_periodic_ewald) then

     do i_atom_2 = 1, n_occ_atoms
        if ( myid == task_list(i_atom_2) ) then

           ! reinitialize the physical quantities of interest
           delta_v_hartree_aux = 0

           ! get current integration point coordinate
           coord_current(:) = coords(:, i_atom_2)


           if (forces_on) then

             call interpolate_extended_hartree_non_periodic_ewald( coord_current, delta_v_hartree_aux, &
                                                                   v_hartree_gradient_temp             )


           else

             call interpolate_extended_hartree_non_periodic_ewald( coord_current, delta_v_hartree_aux )

           end if


        end if
     end do
    
  end if  ! n_periodic > 0

   !-------shanghui begin parallel------  
  call sync_sum_up_whole_potential_shanghui(  &
       hartree_multipole_error)
  call sync_vector(hellman_feynman_term,3*n_atoms) 
   !-------shanghui end parallel------  

  

  ! write root-mean square error of Hartree potential
  hartree_multipole_error =   &
       sqrt(hartree_multipole_error)

     write(info_str,'(2X,A,1X,E14.6)')  &
          "| RMS delta_rho error from multipole expansion :",   &
          hartree_multipole_error
     call localorb_info(info_str, use_unit,'(A)', OL_norm )
     

  ! And finally, for periodic systems, we must remember to shift the potential zero 
  ! to the average real-space potential in periodic systems, also for energy-related quantities
  if (n_periodic.gt.0) then

      call sync_average_potential ( average_delta_v_hartree_real )

      average_delta_v_hartree_real = average_delta_v_hartree_real/cell_volume
      previous_average_delta_v_hartree_real = average_delta_v_hartree_real

      write (info_str,'(2X,A,1X,F15.8,A)')  &
           "| Average real-space part of the electrostatic potential :",   &
           average_delta_v_hartree_real*hartree, " eV"
      call localorb_info(info_str, use_unit,'(A)', OL_norm )

      !test
      !  Ewald_zero_level = 0.d0
      !  do i_atom_2 = 1, n_atoms
      !    Ewald_zero_level = Ewald_zero_level & 
      !      - sqrt(pi4) * multipole_moments(1, i_atom_2)*pi*Ewald_radius**2/cell_volume
      !  end do
      !  if (myid.eq.0) then
      !    write(use_unit,'(2X,A,1X,F15.8,A)')  &
      !      "| Ewald zero level :",   &
      !      Ewald_zero_level*hartree, " eV"
      !  end if
      !test end

      potential(:) = potential(:) - average_delta_v_hartree_real



 

       ! If requested, we can now write the average electrostatic potential 
       ! in the vacuum region of a slab calculation
       ! Notice that only the reciprocal-space part and the average real-space potential
       ! contribute here.
       if( out_vacuum_potential)then

         write(info_str,*) ' | Potential at the vacuum (eV):'
         call localorb_info(info_str, use_unit,'(A)')    

         do i_z = 1, out_vacuum_potential_z_grid, 1

            if(out_vacuum_potential_z_grid==1) then
               out_vacuum_potential_z = out_vacuum_potential_z1
            else
              out_vacuum_potential_z = out_vacuum_potential_z1 + (i_z-1)*(out_vacuum_potential_z2 - out_vacuum_potential_z1) &
                /dble(out_vacuum_potential_z_grid-1)
            endif

            call evaluate_potential_at_vacuum( out_vacuum_potential_z,  out_vacuum_potential_x_grid, &
              out_vacuum_potential_y_grid, dip_gradient,  dip_origin, dip_lenght, average_delta_v_hartree_real)

         end do

       end if

  end if

!  if(use_batch_permutation > 0) then
!
!    call permute_point_array_back(n_bp,1,potential,potential_std)
!    deallocate(potential)
!
!    deallocate(rho)
!
!  endif

  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

  if(allocated(delta_v_hartree)) deallocate(delta_v_hartree)
  if(allocated(v_hartree_free))  deallocate(v_hartree_free)
  if(allocated(rho_multipole))   deallocate(rho_multipole)

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


end subroutine sum_up_whole_potential_shanghui_p0
!******
!------------------------------------------------------------------------------------------------------------
