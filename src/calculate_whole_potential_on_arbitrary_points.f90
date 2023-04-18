!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/ calculate_whole_potential_on_arbitrary_points
!  NAME
!    calculate_whole_potential_on_arbitrary_points
!  SYNOPSIS

subroutine  calculate_whole_potential_on_arbitrary_points &
     ( n_points, cube_points, cube_potential, &
     delta_v_hartree_part_at_zero,  delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, &
     my_partition_tab, &
      my_rho, &
     my_free_hartree_superpos, my_free_rho_superpos, & 
     hartree_delta_energy, en_elec_delta, hartree_multipole_correction, &
!     pot_ion_embed_std,
      en_density_embed, &
     multipole_radius_sq, &
     l_hartree_max_far_distance, &
!rho_multipole_old_std,  &
     outer_potential_radius &
     )

!  PURPOSE
!    The subroutine sums up the hartree potential from the multipole components
!    which have been calculated before hand. The actual summation is done here.
!   WARNING - THIS ROUTINE HAS BEEN ADJUSTED FROM sum_up_whole_potential.f90
!   WARNING - IT HAS BEEN TEST TO GIVE REASONABLE CUBE_POTENTIAL
!   WARNING - ALL (!!!) OTHER VALUES ARE KEPT FOR CONVENIENCE, TO ALLOW
!   WARNING - THE ROUTINE LATER ON. 
!   WARNING
!   WARNING - DO
!   WARNING - NOT 
!   WARNING - TRUST
!   WARNING - THEM
! one particular problem here is that the average potential is calculated on the provided grid, which is
! neither sufficently dense nor necessarily encompasses all of the unit cell

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
  use pseudodata

  implicit none

!  ARGUMENTS


   integer :: n_points !Number of points 
   real*8,dimension(3,n_points) :: cube_points !Coordinates of points
   real*8, dimension(n_points) :: cube_potential !Potential 

  real*8, dimension(n_atoms)                     :: delta_v_hartree_part_at_zero
  real*8, dimension(3, n_atoms)                  :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments
  real*8, dimension(n_points) :: my_partition_tab

  real*8, dimension(n_spin,n_points) :: my_rho
  real*8, dimension(n_points)        :: my_free_hartree_superpos
  real*8, dimension(n_points)        :: my_free_rho_superpos
  real*8, dimension(n_points)        :: my_pot_ion_embed !not used at all yet
  real*8, dimension(n_points)        :: my_rho_multipole_old



  real*8 :: hartree_delta_energy
  real*8 :: en_elec_delta
  real*8 :: hartree_multipole_correction

  real*8 :: en_density_embed
  real*8, dimension(n_atoms)              :: multipole_radius_sq
  integer, dimension( n_atoms)            :: l_hartree_max_far_distance
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius

  real*8, dimension(3) :: dummy !dummy variable

!  INPUTS
! o delta_v_hartree_part_at_zero -- Hartree potential at origin of the atoms from different multipole components
! o delta_v_hartree_deriv_l0_at_zero -- Derivative of Hartree potential at origin of the atoms
! o multipole_moments -- mutlpole moments of the Hartree potential
! o partition_tab -- values of partition function
! o rho_std -- electron density
! o free_hartree_superpos_std -- superposition of the free atoms Hartree potential
! o free_rho_superpos_std -- superposition of the free atoms charge density 
! o pot_ion_embed_std -- embedded potential of ions
! o en_density_embed -- embedded electron density
! o multipole_radius_sq -- outer radius of multipole components
! o l_hartree_max_far_distance -- maximum l-components of for the far distance Hartree potential (periodic systems)
! o rho_multipole_old_std -- multipole components for delta charge
!
!  OUTPUT
! o potential_std -- Hartree potential
! o hartree_delta_energy -- Hartree energy minus free atoms Hartree energy
! o en_elec_delta -- Hartree energy from electron density only
! o hartree_multipole_correction -- multipole correction of Hartree energy
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
!  real*8, allocatable :: delta_v_hartree(:)
!  real*8, allocatable :: v_hartree_free(:)
!  real*8, allocatable :: rho_multipole(:)

  real*8, allocatable :: my_delta_v_hartree(:)
  real*8, allocatable :: my_v_hartree_free(:)
  real*8, allocatable :: my_rho_multipole(:)
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

  real*8 :: dist_to_center !Distance betwen current integratio coordinate and current center

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
  integer :: i_point

  integer :: i_lat
 ! integer :: i_full_points

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

!  integer n_my_batches_work  ! Number of batches actually used
!  integer n_full_points_work ! Number of integration points actually used
!  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array

  !real*8, pointer :: partition_tab(:)
  !real*8, pointer :: rho(:,:)
  !real*8, pointer :: rho_multipole_old(:)
  !real*8, pointer :: free_rho_superpos(:)
  !real*8, pointer :: free_hartree_superpos(:)
  !real*8, pointer :: pot_ion_embed(:)
  !real*8, pointer :: potential(:)

   character*100, parameter :: func='calaculate_whole_potential_on_arbitrary_points'

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  !if(use_batch_permutation > 0) then
  !  call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  !else
  !  call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  !endif
  
    hartree_force_l_add = 0


  !OTH: Block not supported
  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  !n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then
     call aims_stop('Load balancing not supported in this function. Sorry for the inconvenience.',func)

!     n_my_batches_work = batch_perm(n_bp)%n_my_batches
!    n_full_points_work = batch_perm(n_bp)%n_full_points
!
!    batches_work => batch_perm(n_bp)%batches
!    partition_tab => batch_perm(n_bp)%partition_tab
!
!    allocate(rho(n_spin,n_full_points_work))
!    call permute_point_array(n_bp,n_spin,rho_std,rho)
!
!    if(flag_delta_rho_in_multipole) then
!       allocate(rho_multipole_old(n_full_points_work))
!       call permute_point_array(n_bp,1,rho_multipole_old_std,rho_multipole_old)
!    else
!       nullify(rho_multipole_old)
!    endif
!
!    allocate(free_rho_superpos(n_full_points_work))
!    call permute_point_array(n_bp,1,free_rho_superpos_std,free_rho_superpos)
!
!    allocate(free_hartree_superpos(n_full_points_work))
!    call permute_point_array(n_bp,1,free_hartree_superpos_std,free_hartree_superpos)
!
!    if (use_embedding_potential) then
!       allocate(pot_ion_embed(n_full_points_work))
!       call permute_point_array(n_bp,1,pot_ion_embed_std,pot_ion_embed)
!    else
!       nullify(pot_ion_embed)
!    endif
!
!    allocate(potential(n_full_points_work))
!
!  else
!
!    n_my_batches_work = n_my_batches
!    n_full_points_work = n_full_points
!    batches_work => batches
!    partition_tab => partition_tab_std
!    rho => rho_std
!    rho_multipole_old => rho_multipole_old_std
!    free_rho_superpos => free_rho_superpos_std
!    free_hartree_superpos => free_hartree_superpos_std
!    pot_ion_embed => pot_ion_embed_std
!    potential => potential_std
!
  endif

!  if(get_batch_weights) then
!    allocate(batch_times(n_my_batches_work))
!    batch_times(:) = 0
!  endif

  !-----------------------------------------------------------------------------
  ! OTH --------------------END NOT SUPPORTED BLOCK


  !allocate(delta_v_hartree(n_full_points_work),stat=info)
  !call check_allocation(info, 'delta_v_hartree')

  !allocate(v_hartree_free(n_full_points_work),stat=info)
  !call check_allocation(info, 'v_hartree_free')

  !allocate(rho_multipole(n_full_points_work),stat=info)
  !call check_allocation(info, 'rho_multipole')
  
  allocate(my_delta_v_hartree(n_points),stat=info)
  call check_allocation(info, 'my_delta_v_hartree')

  allocate(my_v_hartree_free(n_points),stat=info)
  call check_allocation(info, 'my_v_hartree_free')

  allocate(my_rho_multipole(n_points),stat=info)
  call check_allocation(info, 'my_rho_multipole')
  

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

  !  initialize index_lm
  i_index = 0
  do i_l = 0, l_pot_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

!  call hartree_potential_real_coeff(index_lm, multipole_moments, &
!       l_hartree_max_far_distance, n_centers_hartree_potential )

!  if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
!     call update_outer_radius_l( outer_potential_radius, multipole_moments, multipole_radius_sq, & 
!                                 l_hartree_max_far_distance, index_lm)
!
!     do i_atom_2 = 1, n_atoms
!        adap_outer_radius_sq(i_atom_2) = maxval(outer_potential_radius(:,i_atom_2))
!     end do    
!  else
!     adap_outer_radius_sq = 1e8    
!  end if



!  if(n_periodic > 0)then
!     call evaluate_hartree_recip_coef( index_lm, multipole_moments, &
!                                       l_hartree_max_far_distance )
!  else if (use_hartree_non_periodic_ewald) then
!     call calculate_extended_hartree_non_periodic_ewald( l_hartree_max_far_distance, &
!                                           multipole_radius_sq, adap_outer_radius_sq  )
!  end if

  ! Initialize cube_potential (which is an output variable and runs over cube_points(3,n_points))
  cube_potential = 0.d0
  my_pot_ion_embed=0d0


  ! initialize energies and force components
!  hartree_delta_energy = 0.d0
!  en_elec_delta = 0.d0
!  hartree_multipole_correction = 0.d0
!  hartree_multipole_error = 0.d0
!  en_density_embed = 0.d0


  dip_gradient = 0
  dip_origin = 0
  dip_lenght = 0

  if(n_periodic > 0  )then
     
     if(use_dipole_correction .or. calculate_work_function)then
        if (dipole_correction_method=='potential') then
            call evaluate_dipole_correction & 
           ( previous_average_delta_v_hartree_real, dip_gradient, dip_origin, dip_lenght, .false. )
        elseif(dipole_correction_method=='dipole') then
            call evaluate_dipole_correction_from_dipole &
             (previous_average_delta_v_hartree_real, dip_gradient, dip_origin, dip_lenght,.false.)
        endif
     end if
  end if
 

  ! There are two separate loops over integration points:
  ! 1) A loop over the integration grid to tabulate the multipole components and densities at each
  !    grid point one by one
  ! 2) A second loop over the grid to integrate over products rho_multipole(r)*v_multipole(r)

  ! initialize the physical quantities that are tabulated over the entire integration grid
!  rho_multipole    = 0.d0
!  delta_v_hartree  = 0.d0
!  v_hartree_free   = 0.d0

  my_rho_multipole    = 0.d0
  my_delta_v_hartree  = 0.d0
  my_v_hartree_free   = 0.d0

  ! The following averages are only needed in periodic systems
  average_delta_v_hartree_real = 0.d0


!  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
!  time0 = mpi_wtime()




    ! First loop over the grid: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the supplied grid
 

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

          !do i_atom_2 = 1, n_occ_atoms, 1
          !rather than running over the integration grid, we run over points
           do i_point = 1,n_points,1

            ! reinitialize the physical quantities of interest
            delta_v_hartree_aux = 0.d0



            ! get current integration point coordinate
!  OTH Change: running over points

            coord_current(:)=cube_points(:,i_point)
            !do i_coord = 1, 3, 1
            !  coord_current(i_coord) = coords(i_coord, i_atom_2)
            !end do


            ! Tabulate distances and directions to all atoms -
            ! including relative positions on logarithmic and
            ! radial integration grids.

            call tab_single_atom_centered_coords_p0 &
                 ( current_center, coord_current,  &
                 dist_tab_sq,  &
                 dir_tab )

            ! if (i_atom_2.eq. current_center ) then
            ! if (dist_tab_sq.lt.1e-6) then
            !  dist_tab_sq = r_grid_min(species(i_atom_2))**2 + 1e-15
            !end if



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


             !  call tab_single_atom_centered_coords_radial_log_p0 &
             !       ( current_center, dist_tab_sq, dir_tab,  &
             !       dist_tab_in, i_r, i_r_log, dir_tab_in )

               ! for all inner atoms, we need ylm functions and their gradients explicitly
             !  call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)


               ! Now loop over all inner atoms (those which are close enough
               ! to the current integration point so we need explicit numerical splines)
               !     partitioned
               !     hartree potentials need to be summed up
               !     according to Delley (eq. 12c)

             !  l_h_dim = (l_hartree(species_center(current_center))+1)**2

               ! obtain spline-interpolated values of the multipole components
               ! of the partitioned Hartree potential, splined on the logarithmic
               ! integration grid

             !  call spline_vector_v2 &
             !       ( i_r_log, &
             !       current_delta_v_hart_part_spl, &
             !       (l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid, &
             !       n_grid(species_center(current_center)), &
             !       l_h_dim, &
             !       delta_v_hartree_multipole_component)

               dist_to_center=sqrt((coord_current(1)-coords_center(1,current_center))**2 + &
                              (coord_current(2)-coords_center(2,current_center))**2 + &
                              (coord_current(3)-coords_center(3,current_center))**2)
            !  if (i_atom_2.eq. current_center ) then
                if (dist_to_center<1e-6) then
                  delta_v_hartree_multipole_component(1) =  delta_v_hartree_part_at_zero(current_spl_atom)

                  do i_lm = 2, l_h_dim


                     delta_v_hartree_multipole_component(i_lm) = 2* current_delta_v_hart_part_spl(i_lm,1,1) &
                          - current_delta_v_hart_part_spl(i_lm,1,2)

                  end do
               end if

               ! sum up the Hartree potential contribution from the present i_atom_in
            !   delta_v_hartree_aux = delta_v_hartree_aux + &
            !        ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, &
            !        ylm_tab, 1 )
!
!               if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then
!                  call far_distance_hartree_Fp_periodic_single_atom &
!                       (current_spl_atom, i_center, &
!                       dist_tab_in, l_hartree_max_far_distance, .true., .false., &
!                       multipole_radius_sq(current_spl_atom),                      &
!                       sqrt( adap_outer_radius_sq(current_spl_atom) )                 )
!               end if

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
                       l_hartree_max_far_distance( current_spl_atom), .false. )

                  call far_distance_real_hartree_potential_single_atom_p2 &
                       ( current_spl_atom, delta_v_hartree_aux, &
                       l_hartree_max_far_distance(current_spl_atom), coord_current )

               else
                  ! Now sum up the potential contributions from all far-field atoms
                  ! (analytical multipole potentials only ...)
                  call far_distance_hartree_Fp_periodic_single_atom &
                       (current_spl_atom, i_center, &
                       dist_tab_out, l_hartree_max_far_distance, .false. , .false., &
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


               end if
            end if

            hartree_delta_energy = hartree_delta_energy +  &
                 species_z(species(i_atom_2))*delta_v_hartree_aux

          end do ! loop over atoms
        end if
      end if ! (( (n_periodic > 0) .or. force_new_functional) )

      !
      ! Now follows the real work: Summing the multipole potential, density, and their
      ! derivatives on the integration grid
      ! FIXME: OTH: Change comment to points
      !

      ! Reset grid counter for current Hartree potential center
      !i_full_points = 0

      !do i_batch = 1, n_my_batches_work

        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
      !  do i_index = 1, batches_work(i_batch)%size, 1
         do i_point=1,n_points,1


          ! Only execute if partition_tab is .gt. zero, else
          ! we can run into 1/0 !!!
          if (my_partition_tab(i_point).gt.0.d0) then

            ! get current integration point coordinate
            !coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)
             coord_current(:) = cube_points(:,i_point)


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

                my_delta_v_hartree(i_point) = &
                       my_delta_v_hartree(i_point) + delta_v_hartree_aux

              ! Obtain spline-interpolated values of the multipole density,
              ! this time splined on the radial integration grid
!              call spline_vector_v2 &
!                   ( i_r+1, current_rho_multipole_spl, &
!                   (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
!                   n_radial(species_center(current_center))+2,  &
!                   l_h_dim, &
!                   rho_multipole_component)

                ! sum up the multipole density contribution from the present inner atom
!                rho_multipole_aux = &
!                       ddot ( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)

!                my_rho_multipole(i_point) = my_rho_multipole(i_point) + rho_multipole_aux


                if (n_periodic.eq.0 .and. (.not. force_new_functional)) then
                  ! non-periodic calculation; need electron-only Hartree potential later on.
                  if (.not.empty(  current_spl_atom)) then
                    my_v_hartree_free(i_point) = my_v_hartree_free(i_point) + &
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
                      dist_tab_in, l_hartree_max_far_distance, .true., .false.,  &
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

                    my_v_hartree_free(i_point) = &
                            my_v_hartree_free(i_point) + v_hartree_free_aux
                  end if

                ! Recursively tabulate the radial behavior of all multipole components
                ! These radial functions are a private variable in module
                ! hartree_potential_real.f90, and are reused there later in
                ! subroutine far_distance_hartree_potential_real_single_atom
                call far_distance_hartree_Fp_cluster_single_atom_p2 &
                     ( dist_tab_out, &
                     l_atom_max, .false. )

                  call far_distance_real_hartree_potential_single_atom_p2 &
                         ( i_center, my_delta_v_hartree(i_point), &
                         l_atom_max, coord_current )

              else ! Periodic system


                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                ! nuclear z/r part
                !  d_v_hartree_free_d_r =  &
                !       - 0*species_z(species(current_spl_atom))/dist_tab_sq

                ! d_v_hartree_free_d_r = 0.d0

                ! v_hartree_gradient(:,i_center,i_full_points) =  &
                !      dir_tab_out(:) * d_v_hartree_free_d_r

!                call far_distance_hartree_Fp_periodic_single_atom &
!                     (current_spl_atom, i_center, &
!                     dist_tab_out, l_hartree_max_far_distance, .false.,.false.,  &
!                     multipole_radius_sq(current_spl_atom),                        &
!                     sqrt( adap_outer_radius_sq(current_spl_atom) )                 )

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
                          ( current_center, i_center, my_delta_v_hartree(i_point), &
                          l_hartree_max_far_distance, coord_current )


              end if

            end if ! end periodic system



          end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch --> i_full_points
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
     ! end do ! end loop over batches
    end do  ! end loop over source atoms


!     i_full_points = 0
!      do i_batch = 1, n_my_batches_work

!        if(get_batch_weights) time_start = mpi_wtime()

        ! loop over one batch
!        do i_index = 1, batches_work(i_batch)%size, 1
         do i_point=1,n_points,1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
!          i_full_points = i_full_points + 1

           if (my_partition_tab(i_point).gt.0.d0) then

            ! the reciprocal space contribution needs to be treated separately
            ! from the centers above
            if (n_periodic.gt.0) then
              ! get current integration point coordinate
              !coord_current(:) = batches_work(i_batch) % points(i_index) % coords(:)
               coord_current(:) = cube_points(:,i_point)


              ! before the reciprocal-space component is added, average over the real-space part of 
              ! delta_v_hartree
              average_delta_v_hartree_real = average_delta_v_hartree_real + &
              my_delta_v_hartree(i_point) * my_partition_tab(i_point)

              call update_hartree_potential_recip &
                    ( coord_current, my_delta_v_hartree(i_point) )

              if( use_dipole_correction)then

                i_m = int(floor((coord_current(3) - vacuum_z_level)/dip_lenght))
                dip_coord_current = coord_current(3) - i_m * dip_lenght
                  
                if( dip_coord_current <  vacuum_z_level)then
                  dip_coord_current = dip_coord_current + dip_lenght
                end if
      
                my_delta_v_hartree(i_point)  =  my_delta_v_hartree(i_point) &
                       -  (dip_coord_current-dip_origin) * dip_gradient

              end if

            else if ( use_hartree_non_periodic_ewald ) then
              call aims_stop ('Sorry, non_periodic ewald not allowed for output of cube potential yet. Feel free to change.') 
              ! in principle, the change should work, but I did not test it, so the code stops for now
              call interpolate_extended_hartree_non_periodic_ewald(                       &
                             ! batches_work(i_batch) % points(i_index) % coords,             &
                              coord_current, &
                              !delta_v_hartree(i_full_points)  )
                              my_delta_v_hartree(i_point)  )
            end if ! n_periodic > 0

 !           if(flag_delta_rho_in_multipole)then
 !             my_rho_multipole_old(i_point) =  my_rho_multipole_old(i_point)*(1.d0-multipole_feed_back_parameter) &
 !                   + multipole_feed_back_parameter * my_rho_multipole(i_point)
 !           end if

 !           my_rho_multipole(i_point) = &
 !                my_rho_multipole(i_point) + &
 !                pi4_inv * my_free_rho_superpos(i_point)

            cube_potential(i_point) = &
                 my_delta_v_hartree(i_point) + &
                 my_free_hartree_superpos(i_point)

            !v_hartree_free(i_full_points) = &
            !     v_hartree_free(i_full_points) + &
            !     free_hartree_superpos(i_full_points)
 !           my_v_hartree_free(i_point) = &
 !                my_v_hartree_free(i_point) + &
 !                my_free_hartree_superpos(i_point)

 !           total_rho = 0.d0
 !           do i_spin = 1, n_spin, 1
 !             total_rho = total_rho + my_rho(i_spin, i_point)
 !           enddo

            !            Add contribution to difference hartree-energy
            !            Use correction to the Hartree potential which reduces the error
            !            from linear to quadratic order in the multipole expansion:
            !
            !            E_H = int [V_H * (rho_full - 1/2 rho_multipole)]
            !
            !            Dunlap et al JCP 71, 3396 (1979) or Bastug et al, CPL 211, 119 (1993)
            !
            !            But in essence this is trivial as it only means calc. E_H = int(rho_mult*v_mult)

 !           hartree_delta_energy = hartree_delta_energy + &
 !                ( (my_v_hartree_free(i_point) + &
 !                my_delta_v_hartree(i_point)) * &
 !                ( my_rho_multipole(i_point) ) &
 !                - pi4_inv * my_free_rho_superpos(i_point) * &
 !                my_v_hartree_free(i_point) ) * my_partition_tab(i_point)
            
 !           en_elec_delta = en_elec_delta + &
 !                ( (my_v_hartree_free(i_point) + &
 !                my_delta_v_hartree(i_point)) * &
 !                ( my_rho_multipole(i_point) ) &
 !               - pi4_inv * my_free_rho_superpos(i_point) * &
 !                my_v_hartree_free(i_point) ) * my_partition_tab(i_point)

            !hartree_multipole_correction = &
            !     hartree_multipole_correction + &
            !     ( (my_v_hartree_free(i_point) + &
            !     my_delta_v_hartree(i_point)) * &
            !     ( total_rho - my_rho_multipole(i_point) ) ) &
            !     *my_partition_tab(i_point)

            !hartree_multipole_error = &
            !     hartree_multipole_error  + &
            !     ( total_rho - my_rho_multipole(i_point) )**2 &
            !     * my_partition_tab(i_point)



            ! finally, after all is said and done, we add a possible external
            ! embedding potential to the total electrostatic potential
            if (use_embedding_potential) then
               call embedding_potential(cube_points(:,i_point),my_pot_ion_embed(i_point),dummy)
               cube_potential(i_point)=cube_potential(i_point)+my_pot_ion_embed(i_point)
            !  call aims_stop('Use of embedding potential on regular grid not implemented yet. Feel free to change')
            !  if (full_embedding) then
            !    write(798,*) pot_ion_embed(i_full_points)
            !    potential(i_full_points) = potential(i_full_points) +  &
            !           pot_ion_embed(i_full_points)
!
!              end if
              !en_density_embed = en_density_embed + my_partition_tab(i_point) * &
              !      pot_ion_embed(i_point) * total_rho
            end if

            if (use_embedding_pp) then 
!              call aims_stop('Use of Pseudopotentials on regular grid not implemented yet. Feel free to change')

                !potential(i_full_points) = potential(i_full_points) +  &
                !       whole_local_pseudpot_on_intgrid(i_full_points)

            end if

           end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch --> loop over points
   !     if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      !end do ! end loop over batches






  ! Get work time and total time after barrier
!  time_work = mpi_wtime()-time0
!  call mpi_barrier(mpi_comm_global,info)
!  time_all = mpi_wtime()-time0
!  call sync_real_number(time_work)
!  call sync_real_number(time_all)
!  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for potential: real work ', &
!     time_work,' s, elapsed ',time_all,' s'
!  call localorb_info(info_str, use_unit,, "(A)", OL_norm)


!  if(get_batch_weights) then
!    call set_batch_weights(n_bp,batch_times)
!    deallocate(batch_times)
!  endif

  if(n_periodic > 0  )then
     ! The Fourier part of the potential - in this loop, we no longer compute anything on the grid,
     ! but rather we compute the reciprocal-space part of the potential on each nucleus, and 
     ! add this to the energy of the nuclei in the Hartree potential of the electrons.

!FIXME: OTH: Change from i_atom_2 i_point necessary? 
!      do i_point=1,n_points,1
!     !do i_atom_2 = 1, n_occ_atoms, 1
!        if (myid.eq.task_list(i_atom_2)) then !paralellization
!
!           ! reinitialize the physical quantities of interest
!           delta_v_hartree_aux = 0.d0
!
!           ! get current integration point coordinate
!           !coord_current(:) = coords(:, i_atom_2)
!            coord_current(:) = cube_points(:,i_point)
!
!           delta_v_hartree_aux = 0.d0
!           call update_hartree_potential_recip( coord_current, delta_v_hartree_aux )
!           if( use_dipole_correction)then
!
!              i_m = int(floor((coord_current(3) - vacuum_z_level)/dip_lenght))
!              coord_current(3) = coord_current(3) - i_m * dip_lenght
!              
!              if( coord_current(3) <  vacuum_z_level)then
!                 coord_current(3) = coord_current(3) + dip_lenght
!              end if
!              
!              delta_v_hartree_aux = delta_v_hartree_aux -  (coord_current(3)-dip_origin) * dip_gradient
!
!           end if !paralellization
!
!
!           hartree_delta_energy = hartree_delta_energy +  &
!                species_z(species(i_atom_2))*delta_v_hartree_aux
!
!        end if
!!     end do
     
  else if (use_hartree_non_periodic_ewald) then
     call aims_stop('This part of the code has not been adjusted yet')

     !Fixme: OTH: Adjust to points??
     do i_atom_2 = 1, n_occ_atoms
        if ( myid == task_list(i_atom_2) ) then

           ! reinitialize the physical quantities of interest
           delta_v_hartree_aux = 0

           ! get current integration point coordinate
           coord_current(:) = coords(:, i_atom_2)

           call interpolate_extended_hartree_non_periodic_ewald(  &
                                      coord_current, delta_v_hartree_aux )

           hartree_delta_energy = hartree_delta_energy +  &
                species_z(species(i_atom_2)) * delta_v_hartree_aux


        end if
     end do
    
  end if  ! n_periodic > 0

  ! Now we synchronize all total energy related quantities
  call sync_sum_up_whole_potential(  &
       hartree_delta_energy,  &
       hartree_multipole_correction,  &
       hartree_multipole_error, &
       en_density_embed, &
       en_elec_delta )

  ! write root-mean square error of Hartree potential
  !hartree_multipole_error =   &
  !     sqrt(hartree_multipole_error)
!
  !   write(info_str,'(2X,A,1X,E14.6)')  &
  !        "| RMS charge density error from multipole expansion :",   &
  !        hartree_multipole_error
  !   call localorb_info ( info_str, use_unit,'(A)', OL_norm )

  ! And finally, for periodic systems, we must remember to shift the potential zero 
  ! to the average real-space potential in periodic systems, also for energy-related quantities
  if (n_periodic.gt.0) then

!      call sync_average_potential ( average_delta_v_hartree_real )
!
!      average_delta_v_hartree_real = average_delta_v_hartree_real/cell_volume
!      previous_average_delta_v_hartree_real = average_delta_v_hartree_real
!
!      write (info_str,'(2X,A,1X,F15.8,A)')  &
!           "| Average real-space part of the electrostatic potential :",   &
!           average_delta_v_hartree_real*hartree, " eV"
!      call localorb_info ( info_str, use_unit,'(A)', OL_norm )

      !test
      !  Ewald_zero_level = 0.d0
      !  do i_atom_2 = 1, n_atoms
      !    Ewald_zero_level = Ewald_zero_level & 
      !      - sqrt(pi4) * multipole_moments(1, i_atom_2)*pi*hartree_conv**2/cell_volume
      !  end do
      !  if (myid.eq.0) then
      !    write(use_unit,'(2X,A,1X,F15.8,A)')  &
      !      "| Ewald zero level :",   &
      !      Ewald_zero_level*hartree, " eV"
      !  end if
      !test end

      !potential(:) = potential(:) - average_delta_v_hartree_real
      ! cube_potential(:) = cube_potential(:)  - average_delta_v_hartree_real

       ! The following addition are needed to treat charged systems correctly. In effect, we here add the interaction
       ! energy between a constant, neutralizing charge density q/volume and the multipole electrostatic potential 
       ! in the unit cell.
       !
       ! Because the neutralizing charge background is constant, no integral is needed, other than the electrostatic
       ! potential averages in the unit cell. Note, however, that there are other prescriptions for non-constant
       ! neutralizing charge densities. In those cases, explicit integrals between the neutralizing density
       ! and the electrostatic potential would be required. 

       ! Notice that we here use the actual charge from the multipole decomposition of the charge density, rather than the formal charge.
       ! This is a matter of choice; however, we thus guarantee
       ! that any non-neutral terms that enter the calculation of neutral systems by way of numerical noise are properly cancelled.

!       hartree_delta_energy = hartree_delta_energy & 
!         -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_delta_v_hartree_real

!       en_elec_delta = en_elec_delta &
!         -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_delta_v_hartree_real
 
       !!! hartree_delta_energy = hartree_delta_energy & 
       !!!  +  charge * average_delta_v_hartree_real
 
       ! interaction of neutralizing charge background with average_free_es_pot as 
       ! computed once and for all in initialize_grid_storage()

       ! This term should affect only the electrostatic energy of the neutralizing background, and not
       ! en_elec_delta, which is concerned only with actual electrons

!        hartree_delta_energy = hartree_delta_energy & 
!          -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_free_es_pot

       !!! hartree_delta_energy = hartree_delta_energy & 
       !!!   + charge * average_free_es_pot

!       do i_atom_2 = 1, n_occ_atoms, 1
          ! The shift below is needed to generate a consistent total energy. 
          !
          ! For a consistent potential zero, the nuclei also need to 
          ! be shifted correctly. This is done in Hartree delta energy although
          ! the term really does not belong there.
          ! The factor 2 arises since hartree_delta_energy gets a 1/2 due to 
          ! double counting further down the road.
!          hartree_delta_energy = hartree_delta_energy -  &
!                2.d0 * species_z(species(i_atom_2)) * average_delta_v_hartree_real

          ! This accounts for the electrostatic shift of the sum of eigenvalues
          ! The factor is 1.0 as this quantity enters the expression for the kinetic energy
          ! without any prefactors (unlike hartree_delta_energy, which gets a 0.5 in the total energy)
!          en_elec_delta = en_elec_delta -  &
!                1.d0 * species_z(species(i_atom_2)) * average_delta_v_hartree_real

!       enddo


  end if

  if(use_batch_permutation > 0) then
    call aims_stop('Load balancing not supported in this function. Sorry for the inconvenience.', func)

!    call permute_point_array_back(n_bp,1,potential,potential_std)
!    deallocate(potential)
!
!    if(flag_delta_rho_in_multipole) then
!      call permute_point_array_back(n_bp,1,rho_multipole_old,rho_multipole_old_std)
!      deallocate(rho_multipole_old)
!    endif
!
!    deallocate(rho)
!    deallocate(free_rho_superpos)
!    deallocate(free_hartree_superpos)
!    if (use_embedding_potential) then
!       deallocate(pot_ion_embed)
!    endif
!
  endif

  ! Write one last line to bound output
  !write(info_str,*) ' '
  !call localorb_info(info_str,use_unit,'(A)',OL_norm)

  if(allocated(my_delta_v_hartree)) deallocate(my_delta_v_hartree)
  if(allocated(my_v_hartree_free))  deallocate(my_v_hartree_free)
  if(allocated(my_rho_multipole))   deallocate(my_rho_multipole)

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


end subroutine calculate_whole_potential_on_arbitrary_points
!******
!------------------------------------------------------------------------------------------------------------
