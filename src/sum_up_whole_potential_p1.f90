!------------------------------------------------------------------------------------------------------------

!****s* FHI-aims/sum_up_whole_potential
!  NAME
!   sum_up_whole_potential
!  SYNOPSIS

subroutine sum_up_whole_potential_p1 &
     ( delta_v_hartree_part_at_zero,  delta_v_hartree_deriv_l0_at_zero, &
     multipole_moments, &
     partition_tab_std, rho_std, &
     free_hartree_superpos_std, free_rho_superpos_std, potential_std,  &
     hartree_delta_energy, en_elec_delta, hartree_multipole_correction, &
     pot_ion_embed_std, en_density_embed, &
     multipole_forces, forces_on,multipole_radius_sq, &
     l_hartree_max_far_distance, hellman_feynman_forces,energy_deriv_tress_tensor, rho_multipole_old_std,  &
     outer_potential_radius, AS_stress_on, add_free_hartree, local_fo_potential)
!
!  PURPOSE
!    The subroutine sums up the hartree potential from the multipole components
!    which have been calculated before hand. The actual summation is done here.
!
!  USES

  use types, only: dp
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
  use load_balancing, only: use_batch_permutation, batch_perm, get_batch_weights, &
      permute_point_array, permute_point_array_back, set_batch_weights
  ! rho_multipole from hartree_potential_storage conflicts with the local variable here,
  ! so just only use get_rho_multipole_spl() from this module.
  ! Maybe the local variable name should be changed ...
  use hartree_potential_storage, only : get_rho_multipole_spl
  use pseudodata
  use energy_density
  use analytical_stress
  use heat_flux
  use mpe_interface, only: &
         mpe_dc_indvalues => IND_DC, &
         mpe_get_continuum_name, &
         mpe_dc_indices, &
         mpe_check_initialization, &
         mpe_calculate_interaction_energy_with_nuclei, &
         mpe_calculate_nonelectrostatic_energy, &
         mpe_store_energy_contributions, &
         mpe_calculate_reaction_field, &
         mpe_charge_per_dc
  use timing
  use physics, only: multipole_radius_sq_public, n_electrons,elec_delta_atom
  use lpb_solver_utilities, only:  mpb_solver_started,&
    debug_lpb_solver, delta_v_MPB, delta_v_MPB_gradient, &
    v_free_gradient_public, &
    l_dim_SPE_public,&
    l_dim_species_SPE_public,qSPE_plus_delta_rho,epsinf_mpb,dielec_func_mpb,&
    c_mpb, z_mpb, kBT_mpb,phi_zero_mpb,a_mpb,&
    kappainf_mpb,get_cut_f,deps_drho, dalpha_drho,solve_lpbe_only,&
    use_mpbe_free_energy, en_eps_rho,en_alpha_rho, mpbe_no_ion,&!mpb_forces_on, forces_hf_mpb,&
    freeen_LPBE_nonelstat_solv_energy, alpha_func_mpb,solve_debye_only,& 
    use_separate_alpha, mpb_forces_on,&
    delta_v_MPB_gradient_atom,&
    v_free_gradient_public_atom, forces_mpb_force_off,&
    evaluate_v_free_gradient, evaluate_mpb_energies
    !reg_method=='vacandconst'
    !delta_v_MPB_vac_gradient,delta_v_MPB_vac
  use mpb_solver_utilities, only: freeen_MPBE_nonelstat_solv_energy, f_function,&
    ln_cosh_function, cosh_cosh_function, solve_pbe_only, KS_mpb_corr1, KS_mpb_corr2,&
    h_function
  use SPE_solver, only: evaluate_pot_lpb_at_zero, get_rho_multipole_spl_mpb,&
    delta_v_MPB_multipole_component_at_zero, delta_v_MPB_multipole_deriv_at_zero,&
    delta_v_MPB_gradient_at_zero
  use fodft, only: fodft_in_potential
  use rlsy_interface, only: rlsy_h
  use batch_rlsy_interface, only: batch_irr, &
                                  compute_batch_rlsy_full_to_irr_points, &
                                  compute_batch_rlsy_full_to_irr_points_gradient, &
                                  compute_batch_rlsy_irr_to_full_points, &
                                  compute_batch_rlsy_irr_to_full_points_gradient
  implicit none

!  ARGUMENTS


  real*8, dimension(n_atoms)                     :: delta_v_hartree_part_at_zero
  real*8, dimension(3, n_atoms)                  :: delta_v_hartree_deriv_l0_at_zero
  real*8, dimension( ( l_pot_max+1)**2, n_atoms) :: multipole_moments

  real*8, target, dimension(n_full_points)        :: partition_tab_std
  real*8, target, dimension(n_spin,n_full_points) :: rho_std
  real*8, target, dimension(n_full_points)        :: free_hartree_superpos_std
  real*8, target, dimension(n_full_points)        :: free_rho_superpos_std
  real*8, target, dimension(n_full_points)        :: potential_std

  real*8 :: hartree_delta_energy
  real*8 :: hartree_delta_energy_atom
  real*8 :: en_elec_delta
  real*8 :: hartree_multipole_correction

  real*8, target, dimension(n_full_points):: pot_ion_embed_std
  real*8 :: en_density_embed
  real*8, dimension(3, n_atoms)           :: multipole_forces
  logical :: forces_on
  real*8, dimension(n_atoms)              :: multipole_radius_sq
  integer, dimension( n_atoms)            :: l_hartree_max_far_distance
  real*8, dimension(3, n_atoms)           :: hellman_feynman_forces
  real*8, dimension(3, 3)                 :: energy_deriv_tress_tensor
  real*8, target, dimension(n_full_points):: rho_multipole_old_std
  real*8, dimension(0:l_pot_max, n_atoms) :: outer_potential_radius
  logical                                 :: AS_stress_on
  logical, intent(in)                     :: add_free_hartree
!  INPUTS
! o delta_v_hartree_part_at_zero -- Hartree potential at origin of the atoms from different multipole components
! o delta_v_hartree_deriv_l0_at_zero -- Derivative of Hartree potential at origin of the atoms
! o multipole_moments -- multipole moments of the Hartree potential
! o partition_tab_std -- values of partition function
! o rho_std -- electron density
! o free_hartree_superpos_std -- superposition of the free atoms Hartree potential
! o free_rho_superpos_std -- superposition of the free atoms charge density 
! o pot_ion_embed_std -- embedded potential of ions
! o en_density_embed -- embedded electron density
! o forces_on -- are the forces calculated in this round or not?
! o multipole_radius_sq -- outer radius of multipole components
! o l_hartree_max_far_distance -- maximum l-components of for the far distance Hartree potential (periodic systems)
! o rho_multipole_old_std -- multipole components for delta charge
! o add_free_hartree -- whether to add the superposition of free atoms Hartree potential at the end of
!                       calculaton. If false, only delta_V_Ha is returned.
!
!  OUTPUT
! o potential_std -- Hartree potential
! o hartree_delta_energy -- Hartree energy minus free atoms Hartree energy
! o en_elec_delta -- Hartree energy from electron density only
! o hartree_multipole_correction -- multipole correction of Hartree energy
! o multipole_forces -- mutipole forces = the force components of multipole energy correction
! o hellman_feynman_forces -- Helman-Feyman force components
! o energy_deriv_stress_tensor -- Hellmann-Feynman and multipole energy contribution to energy derivatives
!                                 respent of the lattice vectors (used if stress tensor is calculated).
!                                 (This variable does absolutely nothing.)
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

! DB 231112
  real*8, dimension(3) :: direction
  real*8 :: distance_squared, dist
  real*8 :: d_V_loc_dr, i_r_log_2
  integer :: i_pp_species, i_coord_2
!

  real*8 rho_aux
  real*8 pot_aux
  real*8 rho_multipole_aux
  real*8 delta_v_hartree_aux
  real*8 v_hartree_free_aux

  ! MPE solvation 
  real(dp) :: timestamps_mpe_pot(4)
  real(dp) :: mpe_reaction_field_at_point
  real(dp) :: mpe_interaction_energy_nuc
  real(dp) :: mpe_interaction_energy_total_rho

  integer :: mpe_i_dc
  ! start from NA to get invalid values as well, see mpe_dc_indvalues
  real(dp) :: mpe_total_charge

  !SR: MPB solvation
  !we will regularize with free_hartree_superpos, which is the FULL (el+nuc) free atoms potential!!!
  real*8 :: lpb_zero_energy !double counting and nuc-nuc interaction correctionn
  real*8 :: total_rho_mpb
  real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl_mpb
  real*8, dimension((l_pot_max+1)**2) :: rho_multipole_component_mpb
  real*8 :: rho_multipole_aux_mpb
  real*8 :: rho_multipole_gradient_temp_save(n_full_points,3)
  real*8 :: rho_multipole_gradient_temp_mpb(3)

  ! In periodic systems (only), this will be the correct potential zero for the electrostatic potential
  ! The average electrostatic potential from the previous iteration will be saved, as the 
  !     current "chemical potential" (the Fermi level) was shifted by the average electrostatic
  !     potential from the previous iteration.
  real*8 :: average_delta_v_hartree_real
  real*8, save :: previous_average_delta_v_hartree_real = 0.d0

  ! Optionally, we can compute the average electrostatic potential of the real-space part
  ! (i.e., the actual electrostatic potential plus the Ewald compensating potential)
  ! analytically based on the logarithmic grids of the free atoms. The following variables
  ! are used for this purpose.
  logical, dimension(:), allocatable :: have_atom_average_es_pot
  real*8,  dimension(:), allocatable :: atom_average_es_pot
  

  integer :: current_spl_atom
  integer :: current_center
  integer :: atom_of_splines

  integer :: l_atom_max

  !  counters

  integer i_atom_2,i_atom_3
  logical :: doneyet = .false.
  integer i_center
  integer i_batch
  integer i_l
  integer i_m
  integer i_coord
  integer i_coord2
  integer i_index
  integer i_lm
  integer :: i_spin

  integer :: i_lat
  integer :: i_full_points

  !  external functions
  real*8, external :: ddot
  character*120 :: info_str
  real*8:: HF_forces(3)

  real*8 :: HF_temp(3,n_atoms)

  integer :: mpierr
  integer :: info

  integer :: i_iter, n_iter

  real*8:: out_vacuum_potential_z, dip_gradient, dip_origin, dip_lenght, dip_coord_current, dip_gradient2
  integer:: i_z

   integer n_bytes

  ! Load balancing stuff

  integer n_my_batches_work  ! Number of batches actually used
  integer n_full_points_work ! Number of integration points actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array

  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_multipole_old(:)
  real*8, pointer :: free_rho_superpos(:)
  real*8, pointer :: free_hartree_superpos(:)
  real*8, pointer :: pot_ion_embed(:)
  real*8, pointer :: potential(:)
  integer, pointer :: dc_indices_use(:)

  integer n_bp

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  real*8 time_coef, clock_time_coef
  real*8 time_coef_fin, clock_time_coef_fin

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all
  real*8 time1, time_work1, time_all1
  
  real*8 :: cut_pot_limit
  real*8 :: ion_freeen_1, ion_freeen_2
  real*8 :: dummy, pot_current
 
  real*8 :: freen_solv_plus_ion,freeen_sol_plus_solv_plus_ion, freen_est

  ! FODFT
  real*8, dimension(n_full_points) :: local_fo_potential
  integer :: index_mpb

  ! rlsy
  integer :: i_op, i_atom, j_atom
  real*8, dimension(3, n_atoms)           :: multipole_forces_temp
  real*8 :: AS_dde_v_rho_mp_a_sum(3,3)
  real*8 :: AS_dde_v_rho_f_a_sum(3,3)
  real*8 :: NEW_AS_rho_stress_sum(3,3)
  real*8 :: AS_rho_f_f_self_correction_sum(3,3)

  if(use_batch_permutation > 0) then
    call localorb_info("Summing up the Hartree potential with load balancing.", use_unit,'(2X,A)', OL_norm )
  else
    call localorb_info("Summing up the Hartree potential.", use_unit,'(2X,A)', OL_norm )
  endif

   timestamps_mpe_pot = 0
   if (solvent_method .eq. SOLVENT_MPE) then
      mpe_interaction_energy_total_rho = 0
      mpe_charge_per_dc = 0
   endif

  ! if mpb_forces_on ???
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
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    n_full_points_work = batch_perm(n_bp)%n_full_points

    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(rho(n_spin,n_full_points_work))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    if(flag_delta_rho_in_multipole) then
       allocate(rho_multipole_old(n_full_points_work))
       call permute_point_array(n_bp,1,rho_multipole_old_std,rho_multipole_old)
    else
       nullify(rho_multipole_old)
    endif

    allocate(free_rho_superpos(n_full_points_work))
    call permute_point_array(n_bp,1,free_rho_superpos_std,free_rho_superpos)

    allocate(free_hartree_superpos(n_full_points_work))
    call permute_point_array(n_bp,1,free_hartree_superpos_std,free_hartree_superpos)

    if (use_embedding_potential) then
       allocate(pot_ion_embed(n_full_points_work))
       call permute_point_array(n_bp,1,pot_ion_embed_std,pot_ion_embed)
    else
       nullify(pot_ion_embed)
    endif

    allocate(potential(n_full_points_work))

    if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then
       allocate(dc_indices_use(n_full_points_work))
       call permute_point_array(n_bp,mpe_dc_indices,dc_indices_use)
    else
       nullify(dc_indices_use)
    end if
  
  elseif (use_rlsy_symmetry) then

    n_my_batches_work = batch_irr%n_my_batches
    n_full_points_work = batch_irr%n_full_points
    batches_work => batch_irr%batches
    partition_tab => batch_irr%partition_tab

    allocate(rho(n_spin,batch_irr%n_full_points))
    do i_spin = 1, n_spin
      call compute_batch_rlsy_full_to_irr_points(rho_std(i_spin,:),rho(i_spin,:))
    end do

    if(flag_delta_rho_in_multipole) then
      allocate(rho_multipole_old(batch_irr%n_full_points))
      call compute_batch_rlsy_full_to_irr_points(rho_multipole_old_std(:), &
                                                 rho_multipole_old(:))
    else
      nullify(rho_multipole_old)
    end if

    allocate(free_rho_superpos(batch_irr%n_full_points))
    call compute_batch_rlsy_full_to_irr_points(free_rho_superpos_std(:), &
                                               free_rho_superpos(:))
    allocate(free_hartree_superpos(batch_irr%n_full_points))
    call compute_batch_rlsy_full_to_irr_points(free_hartree_superpos_std(:), &
                                               free_hartree_superpos(:))

    if (use_embedding_potential) then
      allocate(pot_ion_embed(batch_irr%n_full_points))
      call compute_batch_rlsy_full_to_irr_points(pot_ion_embed_std(:), &
                                                 pot_ion_embed(:))
    else
      nullify(pot_ion_embed)
    end if

    allocate(potential(batch_irr%n_full_points))
    potential_std = 0.0d0
    !call compute_batch_rlsy_full_to_irr_points(potential_std(:), &
    !                                           potential(:))

    if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then
       ! impossible to reach this point
    else
       nullify(dc_indices_use)
    end if

       
  else

    n_my_batches_work = n_my_batches
    n_full_points_work = n_full_points
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    rho_multipole_old => rho_multipole_old_std
    free_rho_superpos => free_rho_superpos_std
    free_hartree_superpos => free_hartree_superpos_std
    pot_ion_embed => pot_ion_embed_std
    potential => potential_std

    if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then 
       dc_indices_use => mpe_dc_indices
    else
       nullify(dc_indices_use)
    end if
       
  endif

  if(get_batch_weights) then
    allocate(batch_times(n_my_batches_work))
    batch_times(:) = 0
  endif

  !------------------------------------------------------------------------------


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
  
  

  if (solvent_method.eq.SOLVENT_MPB) then
    freeen_MPBE_nonelstat_solv_energy = 0d0
    freeen_LPBE_nonelstat_solv_energy = 0d0
    freen_solv_plus_ion = 0d0
    freeen_sol_plus_solv_plus_ion = 0d0
    freen_est = 0d0
!     call get_cut_f(cut_pot_limit) ! routine to calculate actual value of z*v/kBT at which limit is reached
    if (.not.allocated(current_rho_multipole_spl_mpb)) then
      allocate(current_rho_multipole_spl_mpb &
	    (l_dim_SPE_public, n_max_spline, n_max_radial+2),stat=info)
      call check_allocation(info, 'current_rho_multipole_spl_mpb     ')

    end if
    en_eps_rho = 0d0
    en_alpha_rho = 0d0
  end if


  if (.not.allocated(current_delta_v_hart_part_spl)) then
     allocate(current_delta_v_hart_part_spl &
          ((l_pot_max+1)**2, n_coeff_hartree, n_hartree_grid),stat=info)
     call check_allocation(info, 'current_delta_v_hart_part_spl ')
     current_delta_v_hart_part_spl = 0.d0
  end if

  if (analytic_potential_average) then
     allocate(have_atom_average_es_pot(n_atoms))
     allocate(atom_average_es_pot(n_atoms))
     have_atom_average_es_pot = .false.
     atom_average_es_pot = 0.d0
  end if

  if (forces_on.or.&
      ((solvent_method.eq.SOLVENT_MPB).and.mpb_solver_started)) then
!do we really need this if forces are off? check!
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
  
  if (flag_out_locpot_atom) then
     elec_delta_atom  = 0.0d0
    ! elec_delta_atom_real = 0.0d0
    ! elec_delta_atom_recip = 0.0d0 
  endif

  ! CC: Allocate arrays for energy density
  if (flag_energy_density) then
    if (.not. allocated(ed_local_potential))              allocate(ed_local_potential              (n_full_points) )
    if (.not. allocated(ed_hartree_delta_energy))         allocate(ed_hartree_delta_energy         (n_full_points) )
    if (.not. allocated(ed_hartree_delta_energy_nuclei))  allocate(ed_hartree_delta_energy_nuclei  (n_atoms)       ) 
    if (.not. allocated(ed_hartree_multipole_correction)) allocate(ed_hartree_multipole_correction (n_full_points) )    
     
    ed_local_potential(:)                  = 0.0d0
    ed_hartree_delta_energy(:)             = 0.0d0
    ed_hartree_delta_energy_nuclei(:)      = 0.0d0
    ed_hartree_multipole_correction(:)     = 0.0d0
  end if

  if(n_periodic > 0 .or. force_new_functional  )then
     if(forces_on) then
        hellman_feynman_forces = 0.0d0
     end if
  end if

  ! CC: Allocate Arrays for analytical stress / virial 
  !if ( (forces_on) .and. ( use_analytical_stress ) ) then 
  if ( ( AS_stress_on ) ) then 
    call allocate_and_or_zero_all_AS_arrays(n_full_points_work)
    if (compute_heat_flux) then
      if (.not.allocated(HF_stress_per_atom_MP_CO))  allocate(HF_stress_per_atom_MP_CO  (1:3,1:3,n_atoms))
      if (.not.allocated(HF_stress_per_atom_MP_AT))  allocate(HF_stress_per_atom_MP_AT  (1:3,1:3,n_atoms))
      if (.not.allocated(HF_stress_per_atom_MP_EL))  allocate(HF_stress_per_atom_MP_EL  (1:3,1:3,n_atoms))
      if (.not.allocated(HF_rho_multipole_per_atom)) allocate(HF_rho_multipole_per_atom (n_full_points,n_atoms))
      if (.not.allocated(HF_dde_rho               )) allocate(HF_dde_rho                (1:3,1:3,n_full_points,n_atoms))
      HF_stress_per_atom_MP_CO(:,:,:) = 0.0d0
      HF_stress_per_atom_MP_AT(:,:,:) = 0.0d0
      HF_stress_per_atom_MP_EL(:,:,:) = 0.0d0
      HF_rho_multipole_per_atom(:,:)  = 0.0d0
      HF_dde_rho(:,:,:,:)             = 0.0d0
    end if
  end if


  !! ! CC: Debug: Write/read dipole moments to file
  !! if ( AS_stress_on .and. AS_flag_write_mm_moments) then
  !!   call AS_write_mm_moments( multipole_moments, n_atoms, (l_pot_max+1)**2 )
  !! end if
  !! if ( AS_stress_on .and. AS_flag_read_mm_moments) then
  !!   if ( .not.allocated(AS_mm_save) ) allocate(AS_mm_save(1:n_atoms,1:(l_pot_max+1)**2))
  !!   AS_mm_save(:,:) = 0.0d0
  !!   AS_mm_save      = multipole_moments
  !!   call AS_read_mm_moments( multipole_moments, n_atoms, (l_pot_max+1)**2 )
  !! end if

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
     call update_outer_radius_l( outer_potential_radius, multipole_moments, multipole_radius_sq, & 
                                 l_hartree_max_far_distance, index_lm)

     if ( .not. use_hartree_non_periodic_ewald ) then
        do i_atom_2 = 1, n_atoms
           adap_outer_radius_sq(i_atom_2) = maxval(outer_potential_radius(:,i_atom_2))
        end do
     else
        adap_outer_radius_sq(:) = multipole_radius_sq(:)
     end if

  else
     adap_outer_radius_sq = 1e8    
  end if


  if(n_periodic > 0)then
     ! CC: The so called B and Lambda Matrices for the reciprocal space contr. 
     !     to the stress are determined within this routine

     if (output_priority<=OL_low) then
       call get_timestamps(time_coef,clock_time_coef)
     end if

     if ( AS_stress_on ) then
       call evaluate_hartree_recip_coef( index_lm, multipole_moments, l_hartree_max_far_distance , .true.  )
     else
       call evaluate_hartree_recip_coef( index_lm, multipole_moments, l_hartree_max_far_distance , .false. )
     end if 

     if (output_priority<=OL_low) then
       call get_timestamps(time_coef_fin,clock_time_coef_fin)

       write(info_str,'(2X,A,I5,A,F10.3,A,F10.3,A)') &
         '| Ewald reciprocal-space coefficients, myid = ', myid," : cpu ",time_coef_fin-time_coef, &
         ' s, wall ', clock_time_coef_fin-clock_time_coef,' s.'
       call localorb_allinfo(info_str,use_unit,'(A)',OL_low)
     end if

  else if (use_hartree_non_periodic_ewald) then
     call calculate_extended_hartree_non_periodic_ewald( l_hartree_max_far_distance, &
                                           multipole_radius_sq, adap_outer_radius_sq  )
  end if

  ! Initialize potential (which is an output variable and runs over all grid points
  potential   = 0.d0

  ! initialize energies and force components
  hartree_delta_energy = 0.d0
  hartree_delta_energy_atom = 0.d0

  en_elec_delta = 0.d0
  hartree_multipole_correction = 0.d0
  hartree_multipole_error = 0.d0
  if (forces_on) then
     multipole_forces(:,:) = 0.d0
  end if
  en_density_embed = 0.d0

  
  dip_gradient = 0
  dip_origin = 0
  dip_lenght = 0

  if(n_periodic > 0  )then
     
     if(use_dipole_correction .or. calculate_work_function)then
        if(dipole_correction_method=='potential') then
             write(info_str,"(2X,A)") 'Calculating dipole correction via potential gradient'
             call localorb_info(info_str, use_unit, "(A)", OL_norm)
             call evaluate_dipole_correction & 
            ( previous_average_delta_v_hartree_real, dip_gradient, dip_origin, dip_lenght,.true. )
        elseif(dipole_correction_method=='dipole') then
             write(info_str,"(2X,A)") 'Calculating dipole correction via slab dipole moment'
             call localorb_info(info_str, use_unit, "(A)", OL_norm)
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

  ! If forces_on is set, we have to make 2 iterations, else 1

  if(forces_on) then
    n_iter = 2
  else
    n_iter = 1
  endif

  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_iter = 1, n_iter

    ! First loop over the grid: We run over the Hartree potential center by center, adding the Hartree
    ! contribution of that atom at each point of the integration grid

    
    atom_of_splines = 0
    do i_center = 1, n_centers_hartree_potential, 1

      current_center   = centers_hartree_potential(i_center)
      current_spl_atom = center_to_atom(current_center)

      if (( (n_periodic > 0) .or. force_new_functional) .and. i_iter == 1 ) then

        ! in this case, use multipole_spline to compute
        ! Hartree potential components ON each individual atomic nucleus, 
        ! before doing anything else.

        if (mod(current_spl_atom-1,n_tasks) == myid) then
          ! for this case rho_multipole for current_spl_atom is on myid

          do i_atom_2 = 1, n_atoms,1
            if (empty(i_atom_2)) cycle
            ! if (species_pseudoized(species(i_atom_2))) cycle

            ! do i_atom_2 = 1, n_occ_atoms, 1

            ! reinitialize the physical quantities of interest
              delta_v_hartree_aux = 0.d0

            if (forces_on) then
              v_hartree_gradient_temp = 0.d0
            end if

            ! get current integration point coordinate
            coord_current(:) = coords(:,i_atom_2)
            ! do i_coord = 1, 3, 1
            ! coord_current(i_coord) = coords(i_coord, i_atom_2)
            ! end do

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

            ! At each integration point, the Hartree potential components coming from
            ! different atoms are split into two groups:
            ! Any components coming from atoms close by are evaluated by explicit
            ! numerical splines; far away atoms are simply represented by
            ! an analytical long-distance multipole potential.
            if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom)) then
               ! begin with everything related to n_atoms_in ...

               if(current_spl_atom /= atom_of_splines) then
                 call get_rho_multipole_spl(current_rho_multipole_spl, current_spl_atom)

		             if (solvent_method.eq.SOLVENT_MPB.and.&
		               mpb_solver_started) then
		               call get_rho_multipole_spl_mpb(current_rho_multipole_spl_mpb,&
		                 current_spl_atom)
		             end if
		  


                 !! ! CC DEBUG: Write/read mm components and splines
                 !! if ( AS_stress_on .and. AS_flag_write_mm_moments ) then
                 !!   call AS_write_mm_splines(current_rho_multipole_spl, (l_pot_max+1)**2 ,  n_max_spline, n_max_radial+2, current_spl_atom )
                 !! end if
                 !! if ( AS_stress_on .and. AS_flag_read_mm_moments ) then
                 !!   if ( .not.allocated(AS_mm_splines_save) ) then 
                 !!     allocate(AS_mm_splines_save(1:(l_pot_max+1)**2 ,  1:n_max_spline, 1:n_max_radial+2 ))
                 !!     AS_mm_splines_save = current_rho_multipole_spl 
                 !!   end if
                 !!   call AS_read_mm_splines(current_rho_multipole_spl, (l_pot_max+1)**2 ,  n_max_spline, n_max_radial+2 , current_spl_atom )
                 !! end if

                 if(communication_type.eq.shmem_comm) then
                   n_bytes = (l_pot_max+1)**2 * n_coeff_hartree * n_hartree_grid * 8
                   call aims_shm_get(current_delta_v_hart_part_spl, (current_spl_atom-1)*n_bytes, n_bytes)
                 else
                   call integrate_delta_v_hartree( current_rho_multipole_spl, current_delta_v_hart_part_spl, &
                                                   n_coeff_hartree, current_spl_atom )
                 endif
                 atom_of_splines = current_spl_atom
                 ! This is the first time that we have the actual splined Hartree potential components
                 ! that will be used below. If we are to compute the average electrostatic potential
                 ! analytically (based on the actual splines on the logarithmic grid) we must do it here.
                 if (analytic_potential_average) then
                    if (.not.have_atom_average_es_pot(current_spl_atom)) then
                       call integrate_average_atom_potential &
                          (current_delta_v_hart_part_spl,n_coeff_hartree, &
                           multipole_radius_sq(current_spl_atom),adap_outer_radius_sq(current_spl_atom), &
                           multipole_moments(1,current_spl_atom), &
                           atom_average_es_pot(current_spl_atom),current_spl_atom )
                       have_atom_average_es_pot(current_spl_atom) = .true.
                    end if
                 end if
               endif

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

               ! CC Save potential@atom
               if (AS_stress_on) then
                 ! Num. pot. 
                 AS_v_at_n_temp            = ddot ( l_h_dim, delta_v_hartree_multipole_component, 1, ylm_tab, 1 )
                 AS_v_at_n(i_atom_2)       = AS_v_at_n(i_atom_2)         + AS_v_at_n_temp

                 ! Free pot.
                 if (i_atom_2 == current_center) then
                   AS_v_at_f_temp = free_pot_es_at_zero(species(i_atom_2))
                 else
                   AS_v_at_f_temp = &
                     val_spline( i_r_log, &
                     free_pot_es_spl(1,1,species(current_spl_atom)),  &
                     n_grid(species(current_spl_atom))) 
                 end if
                 AS_v_at_f(i_atom_2) = AS_v_at_f(i_atom_2) + AS_v_at_f_temp
               end if

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

                  ! AJL/Feb2014
                  ! Moved this variable assignment outside the if loop otherwise the variable
                  ! is unassigned for empty sites and we get numerical noise
                  d_v_hartree_free_d_r =  0.0d0

                  if(i_atom_2 == i_center) then

                     ! d_v_hartree_free_d_r =  0.d0

                     ! l=1 components

                     delta_v_hartree_multipole_deriv(2) = delta_v_hartree_deriv_l0_at_zero(1, i_atom_2)
                     delta_v_hartree_multipole_deriv(3) = delta_v_hartree_deriv_l0_at_zero(2, i_atom_2)
                     delta_v_hartree_multipole_deriv(4) = delta_v_hartree_deriv_l0_at_zero(3, i_atom_2)

                  else
                    !  DB 100112: exclude ghost atoms at this point, since those atoms should not have any free_pot_es.
                    !             Otherwise this would lead to a bug when doing force evaluations.

                     if(.not.(empty(current_spl_atom))) then
		                   d_v_hartree_free_d_r =  &
                          val_spline_deriv( i_r_log, &
                          free_pot_es_spl(1,1,species(current_spl_atom)),  &
                          n_grid(species(current_spl_atom))) * log_weight
                     endif
                  end if

                  ! Remove the self interaction == nucleus interaction to itself.
                  v_hartree_gradient_temp = 0.0d0
                  ! obtain gradients for present atom ...
                if ( mpb_forces_on ) then

                  call evaluate_v_hartree_gradient &
                       ( dist_tab_in, dir_tab_in,  &
                       trigonom_tab,  &
                       ylm_tab,  &
                       dylm_dtheta_tab,  &
                       scaled_dylm_dphi_tab,  &
                       delta_v_MPB_multipole_component_at_zero(:,i_center,i_atom_2),  &
                       delta_v_MPB_multipole_deriv_at_zero(:,i_center,i_atom_2), &
                       d_v_hartree_free_d_r, &
                       l_dim_SPE_public, &
                       v_hartree_gradient_temp)

	               else

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
	             
		            endif
	
                 
                ! CC: Save grad_v @atom
                if (AS_stress_on) then
                  do i_coord=1,3
                    do i_coord2=1,3
                      if(i_atom_2 == i_center) then
                        AS_dde_v_at_f_temp(i_coord,i_coord2)     = 0.0d0
                        AS_dde_v_at_n_temp(i_coord,i_coord2)     = 0.0d0
                      else
                        AS_dde_v_at_f_temp(i_coord,i_coord2) = &
                             d_v_hartree_free_d_r &
                           * dir_tab_in(i_coord) * dir_tab(i_coord2)

                        AS_dde_v_at_n_temp(i_coord,i_coord2) = &
                             v_hartree_gradient_temp(i_coord) &
                           * dir_tab(i_coord2) &
                           - AS_dde_v_at_f_temp(i_coord,i_coord2) 
                      end if
                      AS_dde_v_at_f(i_coord,i_coord2,i_atom_2) = &
                           AS_dde_v_at_f(i_coord,i_coord2,i_atom_2) &
                         + AS_dde_v_at_f_temp(i_coord,i_coord2)
                      AS_dde_v_at_n(i_coord,i_coord2,i_atom_2) = &
                           AS_dde_v_at_n(i_coord,i_coord2,i_atom_2) &
                         + AS_dde_v_at_n_temp(i_coord,i_coord2) 
                    end do
                  end do
                end if

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
                     ! v_hartree_gradient_temp =  dir_tab_out(:) * d_v_hartree_free_d_r
                     v_hartree_gradient_temp = 0.d0

		                 if (solvent_method.eq.SOLVENT_MPB.and.mpb_forces_on) then
			                  !we have already calculated the gradient of delta_v_MPB at the position of the
			                  !nuclei in SPE_solver, so just place it here
			                  v_hartree_gradient_temp = delta_v_MPB_gradient_at_zero(:,current_spl_atom,i_atom_2)
		                 else
		                     ! add multipole potential gradients to free-atom gradient
		                     ! The periodic case is calculated later because it uses both outer and inner atoms.
		                     call far_distance_real_gradient_hartree_potential_single_atom_p2 &
			                   ( current_spl_atom, dir_tab, &
			                   v_hartree_gradient_temp, &
			                   l_hartree_max_far_distance(current_spl_atom) )
		                    end if
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

                  ! CC: Since far_distance_real_gradient_hartree_potential_single_atom directly
                  !     add the real space grad v to the existing v_hartree_gradient_temp, we
                  !     save the value here
                  if ( AS_stress_on ) then
                    AS_delta_v_hartree_aux_save     = delta_v_hartree_aux
                    AS_v_hartree_gradient_temp_save = v_hartree_gradient_temp
                  end if 

                  ! Now sum up the far distance parts of the Hartree potential
                  call far_distance_real_hartree_potential_single_atom &
                       ( current_center, i_center, delta_v_hartree_aux, &
                       l_hartree_max_far_distance, coord_current )

                  ! CC : Save potential@atom
                  if ( AS_stress_on ) then
                    AS_v_at_r_temp      = delta_v_hartree_aux - AS_delta_v_hartree_aux_save
                    AS_v_at_r(i_atom_2) = AS_v_at_r(i_atom_2) + AS_v_at_r_temp
                  end if 


                  if(forces_on)then
                     call far_distance_real_gradient_hartree_potential_single_atom &
                          ( current_spl_atom, i_center, dir_tab, &
                          v_hartree_gradient_temp, &
                          l_hartree_max_far_distance, &
                          dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )
                     ! hF forces = HF + Z (dV - Z/r^2 * dir)

                     ! CC: Save grad_v @atom
                     if ( AS_stress_on ) then
                       do i_coord=1,3,1
                         do i_coord2=1,3,1
                           AS_dde_v_at_r_temp(i_coord,i_coord2) = &
                                ( v_hartree_gradient_temp(i_coord) &
                                - AS_v_hartree_gradient_temp_save(i_coord) ) &
                              * dir_tab(i_coord2) 
                           AS_dde_v_at_r(i_coord,i_coord2,i_atom_2) = &
                                AS_dde_v_at_r(i_coord,i_coord2,i_atom_2) &
                              + AS_dde_v_at_r_temp(i_coord,i_coord2)
                         end do
                       end do
                     end if 

                     if(use_dipole_correction .and. i_atom_2 == i_center )then
                       v_hartree_gradient_temp(3) = v_hartree_gradient_temp(3) - dip_gradient
                     end if

                  end if
               end if
            end if

	          !eq. (66) in Comp.Phys.Comm. E_{nuc-nuc} - E_{el-nuc}
	          if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
	           !we already calculated this energy term in SPE_solver.f90
	            continue
	          else
	              hartree_delta_energy = hartree_delta_energy +  &
		          species_z(species(i_atom_2))*delta_v_hartree_aux
	              hartree_delta_energy_atom = hartree_delta_energy_atom +  &
		          species_z(species(i_atom_2))*delta_v_hartree_aux
	          end if
            
            if (flag_out_locpot_atom) then
               ! if (i_center .eq. 1) then 
                  elec_delta_atom(i_atom_2) = elec_delta_atom(i_atom_2)+delta_v_hartree_aux
               !else 
               !   elec_delta_atom(i_atom_2) = elec_delta_atom(i_atom_2)+delta_v_hartree_aux
               !endif
            endif

            ! CC: Save nuclear contributions to hartree_delta_energy
            ! This is the electrostatice energy of i_atom_2 in the field of all other atoms (index i_center)
            if (flag_energy_density) then
              if (i_center .eq. 1) then 
                ed_hartree_delta_energy_nuclei(i_atom_2) = species_z(species(i_atom_2))*delta_v_hartree_aux
              else
                ed_hartree_delta_energy_nuclei(i_atom_2) = &
                     ed_hartree_delta_energy_nuclei(i_atom_2) &
                   + species_z(species(i_atom_2)) * delta_v_hartree_aux
              end if
            end if

            if(forces_on)then

              ! AJL/Feb2014: What is the purpose of this if statement??

               if(i_atom_2 == i_center )then
                  HF_forces =  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               else
                  HF_forces=  species_z(species(i_atom_2)) * v_hartree_gradient_temp(1:3)
               end if

               hellman_feynman_forces(1:3, i_atom_2 ) = &
                    hellman_feynman_forces(1:3, i_atom_2) &
                  + HF_forces(1:3)               

               ! CC: Sum everything up:
               if (AS_stress_on) then

                 ! Here comes the HF stress for the nuclei: NB: no Jacobian, since sum and not integration!
                 do i_coord = 1, 3, 1
                   do i_coord2 = 1, 3, 1
                     NEW_AS_at_stress(i_coord,i_coord2) = &
                          NEW_AS_at_stress(i_coord,i_coord2) &
                        + ( 0.5d0 * species_z(species(i_atom_2)) &
                        * ( AS_dde_v_at_f_temp(i_coord,i_coord2) &
                          + AS_dde_v_at_n_temp(i_coord,i_coord2) &
                          + AS_dde_v_at_r_temp(i_coord,i_coord2) ) )
                     if (compute_heat_flux) then
                       HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) = & 
                         HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) - &
                         ( 0.5d0 * species_z(species(i_atom_2)) * & 
                          ( AS_dde_v_at_f_temp(i_coord,i_coord2) + & 
                            AS_dde_v_at_n_temp(i_coord,i_coord2) + & 
                            AS_dde_v_at_r_temp(i_coord,i_coord2) ) )
                     end if
                    end do
                 end do

                 ! CC: reset temporary variables 
                 AS_v_at_f_temp          = 0.0d0
                 AS_v_at_n_temp          = 0.0d0
                 AS_v_at_r_temp          = 0.0d0
                 AS_dde_v_at_f_temp      = 0.0d0
                 AS_dde_v_at_n_temp      = 0.0d0
                 AS_dde_v_at_r_temp      = 0.0d0
               end if

               v_hartree_gradient_temp = 0.d0
               
            end if

          end do ! loop over atoms (i_atom_2)
        end if ! parallelization over current_spl_atom
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
	  
            if(i_iter == 2) then 
              v_hartree_gradient_temp = 0
              rho_multipole_gradient_temp = 0

              ! CC: Should be zero anyway
              if ( AS_stress_on ) then
                 AS_v_rho_n_temp     = 0.0d0
                 AS_v_rho_r_temp     = 0.0d0
                 AS_dde_v_rho_f_temp = 0.0d0
                 AS_dde_v_rho_n_temp = 0.0d0
                 AS_dde_v_rho_r_temp = 0.0d0
              end if
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
            if (dist_tab_sq.lt.multipole_radius_sq(current_spl_atom)) then


              if(current_spl_atom /= atom_of_splines) then
                

		            call get_rho_multipole_spl(current_rho_multipole_spl, current_spl_atom)

		            if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
		              call get_rho_multipole_spl_mpb(current_rho_multipole_spl_mpb,&
		                current_spl_atom)
		            end if
                !! ! CC DEBUG: Read/write mm components and splines
                !! if ( AS_stress_on .and. AS_flag_read_mm_moments ) then
                !!   call AS_read_mm_splines(current_rho_multipole_spl, (l_pot_max+1)**2 ,  n_max_spline, n_max_radial+2 , current_spl_atom )
                !! end if
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

                ! CC : Save  v @point
                ! VB : It is possible that this evaluation can go entirely (at least, the spline evaluation)
                if ( AS_stress_on ) then
                  AS_v_rho_n_temp           = delta_v_hartree_aux
                  AS_v_rho_n(i_full_points) = AS_v_rho_n(i_full_points) + AS_v_rho_n_temp

                  ! Free Pot & Free rho
                  AS_v_rho_f_temp = &
                    val_spline( i_r_log, &
                    free_pot_es_spl(1,1,species(current_spl_atom)),  &
                    n_grid(species(current_spl_atom))) 
                  AS_v_rho_f(i_full_points) = AS_v_rho_f(i_full_points) + AS_v_rho_f_temp

                end if

              endif


              ! Obtain spline-interpolated values of the multipole density,
              ! this time splined on the radial integration grid
	            if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
		            ! we need only the multipole rho from the lpb solver
		            call spline_vector_v2 &
		                ( i_r+1, current_rho_multipole_spl_mpb, &
		                l_dim_SPE_public, n_max_spline, n_max_radial+2,  &
		                n_radial(species_center(current_center))+2,  &
		                l_dim_species_SPE_public(species(current_spl_atom)), &
		                rho_multipole_component)
                ! 		else if (regularize_lpb_with_hartree.and.reg_method.eq.'rxnpot') then
                ! 		  ! we need both multipole rho's
                ! ! 		  call spline_vector_v2 &
                ! ! 		      ( i_r+1, current_rho_multipole_spl_mpb, &
                ! ! 		      l_dim_SPE_public, n_max_spline, n_max_radial+2,  &
                ! ! 		      n_radial(species_center(current_center))+2,  &
                ! ! 		      l_dim_species_SPE_public(species(current_spl_atom)), &
                ! ! 		      rho_multipole_component)
                ! 		  call spline_vector_v2 &
                ! 		      ( i_r+1, current_rho_multipole_spl, &
                ! 		      (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
                ! 		      n_radial(species_center(current_center))+2,  &
                ! 		      l_h_dim, &
                ! 		      rho_multipole_component)
                ! 		end if
	            else
		            call spline_vector_v2 &
		                ( i_r+1, current_rho_multipole_spl, &
		                (l_pot_max+1)**2, n_max_spline, n_max_radial+2,  &
		                n_radial(species_center(current_center))+2,  &
		                l_h_dim, &
		                rho_multipole_component)
	            end if

	      

              if(i_iter == 1) then
  
		            if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
		                rho_multipole_aux = &
		            	  ddot ( l_dim_species_SPE_public(species(current_spl_atom)),&
		            	  rho_multipole_component, 1, ylm_tab, 1)
		            else
		              rho_multipole_aux = &
		            	ddot ( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)
		            end if

		            rho_multipole(i_full_points) = rho_multipole(i_full_points) + rho_multipole_aux
                if (AS_stress_on .and. compute_heat_flux) then
                  HF_rho_multipole_per_atom(i_full_points,current_spl_atom) = HF_rho_multipole_per_atom(i_full_points,current_spl_atom) + rho_multipole_aux
                end if

                if (n_periodic.eq.0 .and. (.not. force_new_functional)) then
                  ! non-periodic calculation; need electron-only Hartree potential later on.
                  if (.not.empty(  current_spl_atom)) then
                    v_hartree_free(i_full_points) = v_hartree_free(i_full_points) + &
                         species_z(species_center(current_center)) / &
                         dist_tab_in
                  endif
                end if
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


                d_rho_free_d_r = pi4_inv *  &
                     val_spline_deriv(i_r_log,  &
                     renormalized_free_rho_spl(1,1,species(current_spl_atom)),  &
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
                     
                     
		            if (solvent_method.eq.SOLVENT_MPB) then
		  
                  ! delta_v_MPB_vac_gradient(:,i_full_points) = &
                  !   delta_v_MPB_vac_gradient(:,i_full_points) + v_hartree_gradient_temp - dir_tab_in(:) * d_v_hartree_free_d_r

                  ! v_hartree_gradient_temp(:) = -delta_v_MPB_gradient(i_full_points,:)  +&
                  !   d_v_hartree_free_d_r*dir_tab_in(:)


		              if (mpb_solver_started) then
		                  !if we regularize_lpb_with vacuum hartree we will do the following later
		                  call spline_deriv_vector_v2 &
		              	( i_r+1,  &
		              	current_rho_multipole_spl_mpb, &
		              	l_dim_SPE_public, n_max_spline, n_max_radial+2,  &
		              	n_radial(species_center(current_center))+2,  &
		              	l_dim_species_SPE_public(species(current_spl_atom)), &
		              	rho_multipole_deriv)

		                  ! splines are now derivatives df/di where i is the grid point index
		                  ! must convert to df/dr = df/di * di/dr
		                  rho_multipole_deriv(1:l_dim_species_SPE_public(species(current_spl_atom))) =&
		                    rho_multipole_deriv(1:l_dim_species_SPE_public(species(current_spl_atom))) * &
		                    radial_weight

		                    call evaluate_rho_gradient &
		              	  ( dist_tab_in, dir_tab_in,  &
		              	  trigonom_tab,  &
		              	  ylm_tab,  &
		              	  dylm_dtheta_tab,  &
		              	  scaled_dylm_dphi_tab,  &
		              	  rho_multipole_component, &
		              	  rho_multipole_deriv,  &
		              	  d_rho_free_d_r, &
		              	  l_dim_species_SPE_public(species(current_spl_atom)), &
		              	  rho_multipole_gradient_temp)
		              end if

		            end if !solvent_method.eq.SOLVENT_MPB)
              
                ! CC: Save grad v / grad rho @point
                if ( AS_stress_on ) then

                  ! Free Pot & Free rho at current point for current center
                  AS_rho_f_temp = pi4_inv * &
                    val_spline( i_r_log, &
                    renormalized_free_rho_spl(1,1,species(current_spl_atom)),  &
                    n_grid(species(current_spl_atom))) 
                  AS_v_rho_f_temp = &
                    val_spline( i_r_log, &
                    free_pot_es_spl(1,1,species(current_spl_atom)),  &
                    n_grid(species(current_spl_atom))) 

                  ! We also determine the gradient of rho_free
                  do i_coord=1,3,1
                    do i_coord2=1,3,1
                      AS_dde_v_rho_f_temp(i_coord,i_coord2) = &
                           d_v_hartree_free_d_r &
                         * dir_tab_in(i_coord) &
                         * dir_tab(i_coord2)
                      AS_dde_rho_free_temp(i_coord,i_coord2) = &
                           d_rho_free_d_r &
                         * dir_tab_in(i_coord) &
                         * dir_tab(i_coord2)
                      AS_dde_v_rho_n_temp(i_coord,i_coord2) = &
                        (v_hartree_gradient_temp(i_coord) * dir_tab(i_coord2)) &
                        - AS_dde_v_rho_f_temp(i_coord,i_coord2)
                      AS_dde_rho_mp_temp(i_coord,i_coord2) = &
                           (rho_multipole_gradient_temp(i_coord) &
                          * dir_tab(i_coord2)) &
                         - AS_dde_rho_free_temp(i_coord,i_coord2)
                      AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) = &
                           AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) &
                         + AS_dde_v_rho_f_temp(i_coord,i_coord2)
                      AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) = &
                           AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) &
                         + AS_dde_v_rho_n_temp(i_coord,i_coord2)
                      AS_dde_rho_free(i_coord,i_coord2,i_full_points) = &
                           AS_dde_rho_free(i_coord,i_coord2,i_full_points) &
                         + AS_dde_rho_free_temp(i_coord,i_coord2) 
                      AS_dde_rho_mp(i_coord,i_coord2,i_full_points) = &
                           AS_dde_rho_mp(i_coord,i_coord2,i_full_points) &
                         + AS_dde_rho_mp_temp(i_coord,i_coord2) 

                      AS_rho_f_f_self_correction(i_coord,i_coord2) = &
                           AS_rho_f_f_self_correction(i_coord,i_coord2) &
                         - ( 0.5d0 * (( AS_rho_f_temp &
                                   * AS_dde_v_rho_f_temp(i_coord,i_coord2)) &
                         + ( AS_dde_rho_free_temp(i_coord,i_coord2) &
                           * AS_v_rho_f_temp ) ) &
                         * partition_tab(i_full_points) )
                      if (compute_heat_flux) then
                        HF_dde_rho(i_coord,i_coord2,i_full_points,current_spl_atom) = & 
                          HF_dde_rho(i_coord,i_coord2,i_full_points,current_spl_atom) + &
                          AS_dde_rho_free_temp(i_coord,i_coord2) + AS_dde_rho_mp_temp(i_coord,i_coord2)
                        HF_stress_per_atom_MP_CO(i_coord,i_coord2,current_spl_atom) = & 
                          HF_stress_per_atom_MP_CO(i_coord,i_coord2,current_spl_atom) + ( 0.5d0 * ( &      
                           ( AS_rho_f_temp * & 
                             AS_dde_v_rho_f_temp(i_coord,i_coord2)               )   &
                          +( AS_dde_rho_free_temp(i_coord,i_coord2)            * & 
                             AS_v_rho_f_temp ) ) &
                          * partition_tab(i_full_points) )
                      end if
                         
                    end do
                    AS_rho_f_f_self_correction(i_coord,i_coord) = &
                         AS_rho_f_f_self_correction(i_coord,i_coord) &
                       - ( 0.5d0 * ((AS_rho_f_temp * AS_v_rho_f_temp)) &
                       * partition_tab(i_full_points))
                    if (compute_heat_flux) then
                      HF_stress_per_atom_MP_CO(i_coord,i_coord,current_spl_atom) = & 
                        HF_stress_per_atom_MP_CO(i_coord,i_coord,current_spl_atom) + ( 0.5d0 * (    &
                          ( AS_rho_f_temp * &
                            AS_v_rho_f_temp ) ) &
                         * partition_tab(i_full_points) )
                    end if
                  end do
                end if

              end if ! (i_iter == 2) - happens only when forces are on

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

              if (i_iter == 2) then
                dir_tab_out = dir_tab / dist_tab_out
              end if


              ! Now sum up the potential contributions from a far-field atom
              ! (analytical multipole potentials only ...)
              if ( n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald ) then

                ! if not periodic, need electronic v_hartree_free part for total energy

                if(i_iter == 1) then
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
                endif

                ! Recursively tabulate the radial behavior of all multipole components
                ! These radial functions are a private variable in module
                ! hartree_potential_real.f90, and are reused there later in
                ! subroutine far_distance_hartree_potential_real_single_atom
                call far_distance_hartree_Fp_cluster_single_atom_p2 &
                     ( dist_tab_out, &
                     l_atom_max, forces_on )

                if(i_iter == 1) then

                  call far_distance_real_hartree_potential_single_atom_p2 &
                         ( i_center, delta_v_hartree(i_full_points), &
                         l_atom_max, coord_current )
                endif

                ! ... and compute force contribution of a far-distance atom
                if (i_iter == 2) then 

                  ! All that's left of the free Hartree potential is the
                  ! nuclear z/r part

                  if(.not. force_new_functional)then
                    d_v_hartree_free_d_r =  &
                           - species_z(species(current_spl_atom))/dist_tab_sq
                  else
                    d_v_hartree_free_d_r =  0.d0
                  end if

                  ! VB: FIXME: This needs to go for the new functional part!
                  v_hartree_gradient_temp(:) =  &
                           dir_tab_out(:) * d_v_hartree_free_d_r

	  
                  ! add multipole potential gradients to free-atom gradient
                  ! The periodic case is calculated later because it uses both outer and inner atoms.
                  call far_distance_real_gradient_hartree_potential_single_atom_p2_c1 &
                        ( current_spl_atom, dir_tab, &
                        v_hartree_gradient_temp, &
                        l_atom_max )                           
                  !SR: for MPB solvation, reg_method='vacandconst'                           
                  ! 		  delta_v_MPB_vac_gradient(:,i_full_points) = &
                  ! 		    delta_v_MPB_vac_gradient(:,i_full_points) + v_hartree_gradient_temp - dir_tab_out(:) *d_v_hartree_free_d_r
                        

                end if !i_iter = 2 or solvent_method.eq.SOLVENT_MPB

              else ! Periodic system or non-periodic Ewald method


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

              end if  ! n_periodic == 0 .and. .not. use_hartree_non_periodic_ewald

            end if  ! end if for separating current_center either far-distance or near part


            ! Now sum up the far distance part of the Hartree potential if Ewald's
            ! decomposition is performed. In the opposite case, this was already done.
            if ( n_periodic > 0 .or. use_hartree_non_periodic_ewald ) then


              if (dist_tab_sq.lt. &
                   max(adap_outer_radius_sq(current_spl_atom),multipole_radius_sq(current_spl_atom))) then

                ! Far distance analytic real-space part of the potential

                ! VB: FIXME: l_atom_max not implemented here although it should be!!

                ! CC: Since far_distance_real_gradient_hartree_potential_single_atom directly
                !     adds the real space grad v to the existing v_hartree_gradient_temp, we
                !     save the value here (cf. delta_v_hartree)
                if ( AS_stress_on ) then 
                  AS_delta_v_hartree_aux_save     = delta_v_hartree(i_full_points)
                  AS_v_hartree_gradient_temp_save = v_hartree_gradient_temp
                end if

                if(i_iter == 1) then


                  call far_distance_real_hartree_potential_single_atom &
                          ( current_center, i_center, delta_v_hartree(i_full_points), &
                          l_hartree_max_far_distance, coord_current )

                  ! CC : Save v @ point   
                  if ( AS_stress_on ) then 
                    AS_v_rho_r_temp           = delta_v_hartree(i_full_points) - AS_delta_v_hartree_aux_save
                    AS_v_rho_r(i_full_points) = AS_v_rho_r(i_full_points) + AS_v_rho_r_temp
                  end if

                endif

                if (i_iter == 2) then

                  ! Far distance analytic real-space part of the potential gradient

                  call far_distance_real_gradient_hartree_potential_single_atom &
                         ( current_spl_atom, i_center, dir_tab, &
                         v_hartree_gradient_temp, &
                         l_hartree_max_far_distance, &
                         dist_tab_sq, adap_outer_radius_sq(current_spl_atom)  )

                  ! CC: Save grad_v 
                  if ( AS_stress_on ) then
                    do i_coord=1,3,1
                      do i_coord2=1,3,1
                        AS_dde_v_rho_r_temp(i_coord,i_coord2) = &
                             (v_hartree_gradient_temp(i_coord) &
                           - AS_v_hartree_gradient_temp_save(i_coord)) &
                           * dir_tab(i_coord2)
                        AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) = &
                             AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) &
                           + AS_dde_v_rho_r_temp(i_coord,i_coord2)
                      end do
                    end do
                  end if 

                    
                  if (use_dipole_correction) then
                   v_hartree_gradient_temp(3) = &
                      v_hartree_gradient_temp(3) - dip_gradient
                  end if

                end if !i_iter = 2 or MPB solvation

              end if !dist_tab .lt....

            end if ! end periodic system


            ! Force calculations in second iteration

            if (i_iter == 2) then

              total_rho = 0.d0
              do i_spin = 1, n_spin, 1
                total_rho = total_rho + rho(i_spin, i_full_points)
              enddo

              if (n_periodic .eq. 0) then

                rho_aux = (total_rho - 0.5d0 * rho_multipole(i_full_points))
                pot_aux = 0.5d0*(v_hartree_free(i_full_points) + delta_v_hartree(i_full_points) )


                do i_coord = 1, 3, 1

                  if(force_new_functional)then

                    ! ! either we have real atoms
                    ! if (mpb_forces_on) then
                    ! force_to_check_hf2(i_coord,center_to_atom(current_center))=&
                    ! force_to_check_hf2(i_coord,center_to_atom(current_center))+&
                    ! total_rho * (v_free_gradient_public_atom(current_center,i_coord,i_full_points)&
                    ! +delta_v_MPB_gradient_atom(current_center,i_coord,i_full_points))*&
                    ! partition_tab(i_full_points)
                    ! end if

                    if(.not.species_pseudoized(species(center_to_atom(current_center)))) then
		                  if (solvent_method.eq.SOLVENT_MPB.and.mpb_forces_on) then
			                    total_rho_mpb = qSPE_plus_delta_rho(i_full_points)
			                    multipole_forces(i_coord,center_to_atom(current_center)) =&
			                      multipole_forces(i_coord,center_to_atom(current_center))&
			                      + (total_rho - rho_multipole(i_full_points) - total_rho_mpb )&
 			                      * (v_free_gradient_public_atom(current_center,i_coord,i_full_points)&
 			                      +delta_v_MPB_gradient_atom(current_center,i_coord,i_full_points))&
			                      * partition_tab(i_full_points)
		                  else    
 			                    multipole_forces(i_coord,center_to_atom(current_center)) = &
			                      multipole_forces(i_coord,center_to_atom(current_center)) &
 			                    	+ (total_rho -  rho_multipole(i_full_points)) &
 			                    	* v_hartree_gradient_temp(i_coord) &
 			                    	* partition_tab(i_full_points)
            
		                  end if

                    else

                      ! or we have to deal with pseudopotentials
                      ! then we have take v_hartree_gradient_temp + the gradient of the local
                      ! pseudopotential

                      i_pp_species = species2pp_species(species(center_to_atom(current_center)))

                      distance_squared = 0.d0
                      do i_coord_2 = 1,3,1

                         direction(i_coord_2) = &
                           coord_current(i_coord_2) - coords(i_coord_2, center_to_atom(current_center)) 

                         distance_squared = distance_squared + &
                           (direction(i_coord_2))**2

                      enddo


                      if(distance_squared.gt.(localpot_outer_radius(i_pp_species))**2 ) then
                        ! extrapolate the local potential according to Coulomb character
	        	            d_V_loc_dr = -pp_charge(i_pp_species)/distance_squared**1.5d0


                      elseif (distance_squared.lt.(pp_r_grid_min(i_pp_species))**2) then
                            d_V_loc_dr = 0.d0
                      else
  
                            dist = sqrt(distance_squared)
                            i_r_log = invert_log_grid(dist, &
                            pp_r_grid_min(i_pp_species), pp_r_grid_inc(i_pp_species))
 
                        ! chain rule: df/dr = df/di *di/dr
 
                        !df/di:
         	              d_V_loc_dr = &
                            val_spline_deriv( i_r_log, &
                            local_pseudopot_spl(1,1,i_pp_species),  &
                            n_points_pp_fn(i_pp_species))

                        !df/di *di/dr:
                            d_V_loc_dr = d_V_loc_dr/log(pp_r_grid_inc(i_pp_species))/dist

                        !here we have to scale with 1/r, since vec(direction) is not normalized
                            d_V_loc_dr = d_V_loc_dr/dist

                      endif

                       multipole_forces(i_coord,center_to_atom(current_center)) = &
                                multipole_forces(i_coord,center_to_atom(current_center)) &
                                - (total_rho -  rho_multipole(i_full_points)) &
                                * (d_V_loc_dr*direction(i_coord)- v_hartree_gradient_temp(i_coord))&
                                * partition_tab(i_full_points)

                    end if
                  else
                  
		      

                    ! Atomic Potential gradient is nonzero for all atoms
                    force_component(i_center) = &
                            rho_aux * &
                            v_hartree_gradient_temp(i_coord)
                       
                    ! Atomic density gradient is zero for all atoms outside
                    ! multipole radius -  consider only inner atoms
                    if(dist_tab_sq.lt.multipole_radius_sq(center_to_atom(current_center))) then
                          
                      force_component(i_center) = &
                               force_component(i_center) - pot_aux * &
                               rho_multipole_gradient_temp(i_coord)
                       
                    end if 

                    if(i_center <= n_occ_atoms) then
                       multipole_forces(i_coord,i_center) = &
                                multipole_forces(i_coord,i_center) &
                                +  force_component(i_center)* partition_tab(i_full_points)
                    endif
                  end if


                end do


              else ! n_periodic > 0

                if (.not.empty(center_to_atom(current_center))) then

                  do i_coord = 1, 3, 1


                    multipole_forces(i_coord,center_to_atom(current_center)) = &
                            multipole_forces(i_coord,center_to_atom(current_center)) &
                            + (total_rho -  rho_multipole(i_full_points)) &
                            * v_hartree_gradient_temp(i_coord) &
                            * partition_tab(i_full_points)


                    ! CC FIXME To avoid saving everything, we could 
                    ! compute part of the stress here. However, this  
                    ! should not be memory or CPU critical wrt. Pulay terms
                    if ( AS_stress_on ) then
                      AS_v_rho_n_temp          = 0.0d0
                      AS_v_rho_r_temp          = 0.0d0
                      AS_dde_v_rho_f_temp(:,:) = 0.0d0
                      AS_dde_v_rho_n_temp(:,:) = 0.0d0
                      AS_dde_v_rho_r_temp(:,:) = 0.0d0
                    end if

                  end do
                endif ! end cycle over ghost atoms
              end if ! n_periodic

            end if ! if (i_iter == 2)

           end if ! end if (partition_tab.gt.0.d0)
        end do  ! end loop over points in a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    end do  ! end loop over source atoms
    if(i_iter == 1)then
      ! next, we integrate the actual energy and force quantities

      if (output_priority<=OL_low) then
        call get_timestamps(time_coef,clock_time_coef)
      end if

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

              if ( AS_stress_on ) then

                AS_delta_v_hartree_aux_save = delta_v_hartree(i_full_points)
                AS_dde_v_rho_k_temp(:,:) = 0.0d0

                ! VB: Flag that is intended to be temporary unless
                ! problems are found with the faster '_v2' version
                if (fast_Ewald) then
                  call update_hartree_potential_recip_v2 &
                      ( coord_current, delta_v_hartree(i_full_points), AS_dde_v_rho_k_temp )
                else
                  call update_hartree_potential_recip &
                      ( coord_current, delta_v_hartree(i_full_points), AS_dde_v_rho_k_temp )
                end if

                ! CC save v @ point
                AS_v_rho_k_temp = delta_v_hartree(i_full_points) - AS_delta_v_hartree_aux_save
                AS_v_rho_k(i_full_points)     = AS_v_rho_k_temp
                do i_coord=1,3,1
                  do i_coord2=1,3,1
                    AS_dde_v_rho_k(i_coord,i_coord2,i_full_points) = AS_dde_v_rho_k_temp(i_coord,i_coord2)
                  end do
                end do

              else
                ! VB: Flag that is intended to be temporary unless
                ! problems are found with the faster '_v2' version
                if (fast_Ewald) then
                  call update_hartree_potential_recip_v2 &
                      ( coord_current, delta_v_hartree(i_full_points) )
                else
                  call update_hartree_potential_recip &
                      ( coord_current, delta_v_hartree(i_full_points) )
                end if
              end if

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

	          if(flag_delta_rho_in_multipole)then
	            rho_multipole_old(i_full_points) =  rho_multipole_old(i_full_points)*(1.d0-multipole_feed_back_parameter) &
	                  + multipole_feed_back_parameter * rho_multipole(i_full_points)
	          end if

	          if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
	            rho_multipole(i_full_points) = &
	                rho_multipole(i_full_points) + & !*1.d0/dielec_func_mpb(i_full_points) + &
	                pi4_inv * free_rho_superpos(i_full_points)
	          else 
	            rho_multipole(i_full_points) = &
	                rho_multipole(i_full_points) + &
	                pi4_inv * free_rho_superpos(i_full_points)
	          end if

            if (AS_stress_on .and. compute_heat_flux) then
              HF_rho_multipole_per_atom(i_full_points,1:n_atoms) = & 
                  HF_rho_multipole_per_atom(i_full_points,1:n_atoms) &
                  + pi4_inv*HF_rho_free_per_atom(i_full_points,1:n_atoms)
            end if

            if (add_free_hartree) then
               potential(i_full_points) = &
                    delta_v_hartree(i_full_points) + &
                    free_hartree_superpos(i_full_points)
            else
               potential(i_full_points) = delta_v_hartree(i_full_points)
            end if

	          !normally v_hartree_free is zero until here
            v_hartree_free(i_full_points) = &
                 v_hartree_free(i_full_points) + &
                 free_hartree_superpos(i_full_points)

            total_rho = 0.d0
            do i_spin = 1, n_spin, 1
              total_rho = total_rho + rho(i_spin, i_full_points)
            enddo

	          if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
	            potential(i_full_points) = &
                      (free_hartree_superpos(i_full_points)+delta_v_MPB(i_full_points))
                    if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
                      index_mpb=1
                    else
                      index_mpb=i_full_points
                    end if
	            call evaluate_mpb_energies(hartree_delta_energy, en_elec_delta,&
		              hartree_multipole_correction,hartree_multipole_error,&
		              potential(i_full_points),&
		              total_rho,rho_multipole(i_full_points),free_rho_superpos(i_full_points),&
		              v_hartree_free(i_full_points),partition_tab(i_full_points),& 
		              !changed at output: energies:
		              freeen_MPBE_nonelstat_solv_energy,&
		              freeen_LPBE_nonelstat_solv_energy,&
		              freen_solv_plus_ion,&
		              freeen_sol_plus_solv_plus_ion,&
		              freen_est,f_function(index_mpb),&
                      h_function(index_mpb),&
                      ln_cosh_function(index_mpb),&
                      cosh_cosh_function(index_mpb),&
                      KS_mpb_corr1,KS_mpb_corr2,i_full_points)
	          else
              ! Add contribution to difference hartree-energy
              ! Use correction to the Hartree potential which reduces the error
              ! from linear to quadratic order in the multipole expansion:
              !
              ! E_H = int [V_H * (rho_full - 1/2 rho_multipole)]
              !
              ! Dunlap et al JCP 71, 3396 (1979) or Bastug et al, CPL 211, 119 (1993)
              !
              ! But in essence this is trivial as it only means calc. E_H = int(rho_mult*v_mult)
                   !non-solvent case
                    hartree_delta_energy = hartree_delta_energy + &
                         ( (v_hartree_free(i_full_points) + &
                         delta_v_hartree(i_full_points)) * &
                         ( rho_multipole(i_full_points) ) &
                         - pi4_inv * free_rho_superpos(i_full_points) * &
                         v_hartree_free(i_full_points) ) * partition_tab(i_full_points)
                    
                    en_elec_delta = en_elec_delta + &
                         ( (v_hartree_free(i_full_points) + &
                         delta_v_hartree(i_full_points)) * &
                         ( rho_multipole(i_full_points) ) &
                         - pi4_inv * free_rho_superpos(i_full_points) * &
                         v_hartree_free(i_full_points) ) * partition_tab(i_full_points)
        
                    hartree_multipole_correction = &
                         hartree_multipole_correction + &
                         ( (v_hartree_free(i_full_points) + &
                         delta_v_hartree(i_full_points)) * &
                         ( total_rho - rho_multipole(i_full_points) ) ) &
                         *partition_tab(i_full_points)
        
                    hartree_multipole_error = &
                         hartree_multipole_error  + &
                         ( total_rho - rho_multipole(i_full_points) )**2 &
                         * partition_tab(i_full_points)
	          end if

            ! CC: Saving the different contributions hartree_delta_energy / hartree_multipole_correction
            ! as well as the local potential (without v_xc, as usual)
            if (flag_energy_density) then

              ed_local_potential(i_full_points) = &
              &  ( delta_v_hartree(i_full_points) + &
              &  free_hartree_superpos(i_full_points) ) * total_rho

              ed_hartree_delta_energy(i_full_points) = &
                   ((v_hartree_free(i_full_points) &
                  +  delta_v_hartree(i_full_points)) &
                 * (rho_multipole(i_full_points)) &
                 - (v_hartree_free(i_full_points)) &
                 * ( free_rho_superpos(i_full_points) * pi4_inv ))

              ed_hartree_multipole_correction(i_full_points) = &
                   ((v_hartree_free(i_full_points) &
                  +  delta_v_hartree(i_full_points)) &
                 * (total_rho - rho_multipole(i_full_points)))

            end if



            ! finally, after all is said and done, we add a possible external
            ! embedding potential to the total electrostatic potential
            if (use_embedding_potential) then
              if (full_embedding) then
                potential(i_full_points) = potential(i_full_points) +  &
                       pot_ion_embed(i_full_points)

              end if
              en_density_embed = en_density_embed + partition_tab(i_full_points) * &
                    pot_ion_embed(i_full_points) * total_rho
            end if

            if (use_embedding_pp) then 
                potential(i_full_points) = potential(i_full_points) +  &
                       whole_local_pseudpot_on_intgrid(i_full_points)

                en_density_embed = en_density_embed + partition_tab(i_full_points) * &
                    whole_local_pseudpot_on_intgrid(i_full_points) * total_rho

            end if

            if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then
               call start_timer(timestamps_mpe_pot(1:2))
               call mpe_calculate_reaction_field( &
                     coord=batches_work(i_batch) % points(i_index) % coords(:), &
                     rho_free=free_rho_superpos(i_full_points)*pi4_inv, &
                     rho_mp=rho_multipole(i_full_points), &
                     rho_total=total_rho, &
                     hartree_potential= v_hartree_free(i_full_points) + &
                        delta_v_hartree(i_full_points), &
                     dc_index=dc_indices_use(i_full_points), &
                     potential_correction=mpe_reaction_field_at_point )
               potential(i_full_points) = potential(i_full_points) +  &
                     mpe_reaction_field_at_point
               ! add incremental energy contributions
               hartree_delta_energy = hartree_delta_energy + &
                     mpe_reaction_field_at_point &
                     * rho_multipole(i_full_points) &
                     * partition_tab(i_full_points)
               mpe_interaction_energy_total_rho = &
                  mpe_interaction_energy_total_rho + &
                     partition_tab(i_full_points) * &
                     mpe_reaction_field_at_point * &
                     total_rho
               call stop_timer(timestamps_mpe_pot, unsynced=.true.)

               mpe_i_dc = dc_indices_use(i_full_points)
               mpe_charge_per_dc(mpe_i_dc) = mpe_charge_per_dc(mpe_i_dc) + &
                  partition_tab(i_full_points) * total_rho
            endif

            ! add the hartree potential for the embedded fragment here
            if (use_fo_potential.and.fo_file_exists) then
                potential(i_full_points) = potential(i_full_points) + &
                    local_fo_potential(i_full_points)
            end if


          end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
        if(get_batch_weights) batch_times(i_batch) = batch_times(i_batch) + mpi_wtime() - time_start
      end do ! end loop over batches
    end if !if i_iter = 1
    

    ! CC: Also compute average of stress derivative of potential (real part)
    if ( (i_iter == 2) .and. (AS_stress_on)) then

      i_full_points = 0
      do i_batch = 1, n_my_batches_work

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

            ! Average over the real-space part of  dv_de
            do i_coord=1,3,1
              do i_coord2=1,3,1
                AS_dde_v_rho_f_a(i_coord,i_coord2) = AS_dde_v_rho_f_a(i_coord,i_coord2) & 
                  + ( AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) * partition_tab(i_full_points) )
                AS_dde_v_rho_mp_a(i_coord,i_coord2) = AS_dde_v_rho_mp_a(i_coord,i_coord2) & 
                  + ( ( AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) & 
                      + AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) ) * partition_tab(i_full_points) )
              end do
            end do 
            !
            AS_average_delta_v_hartree_real = &
                 AS_average_delta_v_hartree_real &
               + ((AS_v_rho_n(i_full_points) + AS_v_rho_r(i_full_points)) &
                 * partition_tab(i_full_points))

            AS_average_v_hartree_free = AS_average_v_hartree_free &
               + (AS_v_rho_f(i_full_points )* partition_tab(i_full_points))

           end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
      end do ! end loop over batches

      ! Average
      if (use_mpi) then
        call sync_matrix(AS_dde_v_rho_mp_a, 3, 3 )
        call sync_matrix(AS_dde_v_rho_f_a,  3, 3 )
        call sync_real_number(AS_average_delta_v_hartree_real)
        call sync_real_number(AS_average_v_hartree_free)
      end if
      if (AS_stress_on .and. use_rlsy_symmetry) then
        AS_dde_v_rho_mp_a_sum = 0.0d0
        do i_op = 1, rlsy_h%spacegroup%n_operation
          AS_dde_v_rho_mp_a_sum(:,:) = AS_dde_v_rho_mp_a_sum(:,:) + &
              matmul( &
              matmul( rlsy_h%spacegroup%op(i_op)%m(:,:), &
                      AS_dde_v_rho_mp_a(:,:)), &
                      transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
        end do
        AS_dde_v_rho_mp_a = AS_dde_v_rho_mp_a_sum

        AS_dde_v_rho_f_a_sum = 0.0d0
        do i_op = 1, rlsy_h%spacegroup%n_operation
          AS_dde_v_rho_f_a_sum(:,:) = AS_dde_v_rho_f_a_sum(:,:) + &
              matmul( &
              matmul( rlsy_h%spacegroup%op(i_op)%m(:,:), &
                      AS_dde_v_rho_f_a(:,:)), &
                      transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
        end do
        AS_dde_v_rho_f_a = AS_dde_v_rho_f_a_sum
      end if
      do i_coord=1,3,1
        do i_coord2=1,3,1
          AS_dde_v_rho_mp_a(i_coord,i_coord2) = AS_dde_v_rho_mp_a(i_coord,i_coord2) / cell_volume
          AS_dde_v_rho_f_a(i_coord,i_coord2)  = AS_dde_v_rho_f_a(i_coord,i_coord2)  / cell_volume
        end do
      end do 
      AS_average_delta_v_hartree_real  = AS_average_delta_v_hartree_real / cell_volume
      AS_average_v_hartree_free        = AS_average_v_hartree_free       / cell_volume

      ! Now that everything is in place, we start computing
      i_full_points = 0
      do i_batch = 1, n_my_batches_work

        ! loop over one batch
        do i_index = 1, batches_work(i_batch)%size, 1

          ! i_full_points is the index that indicates where we are in the entire grid (for external quanities like rho, potential, ...)
          i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

            total_rho = 0.d0
            do i_spin = 1, n_spin, 1
              total_rho = total_rho + rho(i_spin, i_full_points)
            enddo

            ! Here we sum up all the relevant terms
            do i_coord=1,3,1
              do i_coord2=1,3,1
                ! Individual derivatives
                NEW_AS_rho_stress(i_coord,i_coord2) = &
                     NEW_AS_rho_stress(i_coord,i_coord2) &
                   - ( 0.5d0 * (((rho_multipole(i_full_points)) &
                          * ( AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) &
                            + AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) & 
                            + AS_dde_v_rho_k(i_coord,i_coord2,i_full_points) &
                            + AS_dde_v_rho_f(i_coord,i_coord2,i_full_points))) &
                       + ((AS_dde_rho_mp(i_coord,i_coord2,i_full_points) &
                         + AS_dde_rho_free(i_coord,i_coord2,i_full_points)) &
                         * ( AS_v_rho_n(i_full_points) &
                           + AS_v_rho_r(i_full_points) &
                           + AS_v_rho_k(i_full_points) &
                           + v_hartree_free(i_full_points)))) &
                         * partition_tab(i_full_points) ) 
                if (compute_heat_flux) then
                  HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) = HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) - ( 0.5d0 * ( &
                      ( ( HF_rho_multipole_per_atom(i_full_points,1:n_atoms)                                                 ) * & 
                           ( AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) + AS_dde_v_rho_r(i_coord,i_coord2,i_full_points)     & 
                           + AS_dde_v_rho_k(i_coord,i_coord2,i_full_points) + AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) ) ) &
                    + ( ( HF_dde_rho(i_coord,i_coord2,i_full_points,1:n_atoms)                                               ) * & 
                           ( AS_v_rho_n(i_full_points)                      + AS_v_rho_r(i_full_points)                          &
                           + AS_v_rho_k(i_full_points)                      + v_hartree_free(i_full_points)                  ) ) ) * partition_tab(i_full_points) ) 
                end if
                   
                if (use_AS_Jac_in_pulay) then
                  AS_dde_potential(i_coord,i_coord2,i_full_points) = &
                       AS_dde_potential(i_coord,i_coord2,i_full_points) &
                     + AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) &
                     + AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) &
                     + AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) & 
                     + AS_dde_v_rho_k(i_coord,i_coord2,i_full_points) &
                     - AS_dde_v_rho_mp_a(i_coord,i_coord2) 
                else
                  AS_pu_rho_dde_fnrka(i_coord,i_coord2) = &
                       AS_pu_rho_dde_fnrka(i_coord,i_coord2) &
                     + (((total_rho &
                     * ( AS_dde_v_rho_f(i_coord,i_coord2,i_full_points) &
                       + AS_dde_v_rho_n(i_coord,i_coord2,i_full_points) &
                       + AS_dde_v_rho_r(i_coord,i_coord2,i_full_points) & 
                       + AS_dde_v_rho_k(i_coord,i_coord2,i_full_points) &
                       - AS_dde_v_rho_mp_a(i_coord,i_coord2)))) &
                     * partition_tab(i_full_points))
                end if
              end do
              !Jacobi term
              NEW_AS_rho_stress(i_coord,i_coord) = &
                   NEW_AS_rho_stress(i_coord,i_coord) &
                 - ( 0.5d0 * (((rho_multipole(i_full_points)) &
                      * ( AS_v_rho_n(i_full_points) &
                        + AS_v_rho_r(i_full_points) & 
                        + AS_v_rho_k(i_full_points) &
                        + v_hartree_free(i_full_points)))) &
                      * partition_tab(i_full_points)) 
              if (compute_heat_flux) then
                HF_stress_per_atom_MP_EL(i_coord,i_coord,1:n_atoms) = HF_stress_per_atom_MP_EL(i_coord,i_coord,1:n_atoms)  - ( 0.5d0 * ( &
                     ( ( HF_rho_multipole_per_atom(i_full_points,1:n_atoms) )  * & 
                          ( AS_v_rho_n(i_full_points)                      + AS_v_rho_r(i_full_points)                      & 
                          + AS_v_rho_k(i_full_points)                      + v_hartree_free(i_full_points)                  ) ) ) * partition_tab(i_full_points) ) 
              end if
            end do
         
           
           end if ! if (partition_tab(i_full_points).gt.0)

        end do ! end loop over a batch
      end do ! end loop over batches
      
    endif ! if(i_iter == 2).and.(use_analytical_stress)

  enddo ! Loop over n_iter

  if(use_rlsy_symmetry) then
    hartree_delta_energy         = (hartree_delta_energy - hartree_delta_energy_atom ) &
                                   * rlsy_h%spacegroup%n_operation + hartree_delta_energy_atom
    en_elec_delta                = en_elec_delta                * rlsy_h%spacegroup%n_operation
    hartree_multipole_correction = hartree_multipole_correction * rlsy_h%spacegroup%n_operation
    hartree_multipole_error      = hartree_multipole_error      * rlsy_h%spacegroup%n_operation
    average_delta_v_hartree_real = average_delta_v_hartree_real * rlsy_h%spacegroup%n_operation
  end if

  !write(info_str,'(a)') ' hartree potential'
  !call localorb_info(info_str, use_unit, "(A)", OL_norm)
  !do i_full_points = 1, n_full_points
  !  write(info_str,'((f30.15))') , potential_std(i_full_points)
  !  call localorb_info(info_str, use_unit, "(A)", OL_norm)
  !end do

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
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

                 ! get current nucleus coordinate
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

           ! get current nucleus coordinate
           coord_current(:) = coords(:, i_atom_2)

           ! CC: The reciprocal part of the stress does not require the gradients wrt to the positions,
           !     but with respect to k. We have computed the respective derivatives of the hartree_coefficients in
           !        evaluate_hartree_recip_coef at the beginning of this subroutine and now use them here
           delta_v_hartree_aux     = 0.0d0
           AS_dde_v_at_k_temp(:,:) = 0.0d0
           if ( AS_stress_on ) then

             if (fast_Ewald) then
               call update_hartree_potential_recip_v2(coord_current, delta_v_hartree_aux, AS_dde_v_at_k_temp)
             else
               call update_hartree_potential_recip(coord_current, delta_v_hartree_aux, AS_dde_v_at_k_temp)
             end if

             AS_v_at_k(i_atom_2) = delta_v_hartree_aux
             do i_coord=1,3,1
               do i_coord2=1,3,1
                 AS_dde_v_at_k(i_coord,i_coord2,i_atom_2) = AS_dde_v_at_k_temp(i_coord,i_coord2)
               end do
             end do

           else
             if (fast_Ewald) then
               call update_hartree_potential_recip_v2(coord_current, delta_v_hartree_aux)
             else
               call update_hartree_potential_recip(coord_current, delta_v_hartree_aux)
             end if
           end if

           if( use_dipole_correction)then

              i_m = int(floor((coord_current(3) - vacuum_z_level)/dip_lenght))
              coord_current(3) = coord_current(3) - i_m * dip_lenght
              
              if( coord_current(3) <  vacuum_z_level)then
                 coord_current(3) = coord_current(3) + dip_lenght
              end if
              
              delta_v_hartree_aux = delta_v_hartree_aux -  (coord_current(3)-dip_origin) * dip_gradient

           end if


           hartree_delta_energy = hartree_delta_energy +  &
                species_z(species(i_atom_2))*delta_v_hartree_aux
                
           if (flag_out_locpot_atom) then 
               elec_delta_atom(i_atom_2) = &
                  elec_delta_atom(i_atom_2) + & 
                   delta_v_hartree_aux
           endif
           ! CC: Save nuclear contributions to hartree_delta_energy
           ! This is the electrostatice energy of i_atom_2 in the field of all other PBC atoms 
           if (flag_energy_density) then
               ed_hartree_delta_energy_nuclei(i_atom_2) = &
                    ed_hartree_delta_energy_nuclei(i_atom_2) &
                  + species_z(species(i_atom_2)) * delta_v_hartree_aux
           end if

           ! CC: See comment above
           if ( AS_stress_on ) then
             do i_coord=1,3
               do i_coord2=1,3
                 NEW_AS_at_stress(i_coord,i_coord2) = NEW_AS_at_stress(i_coord,i_coord2) + &
                    ( 0.5d0 * species_z(species(i_atom_2)) * AS_dde_v_at_k(i_coord,i_coord2,i_atom_2) )
                 if (compute_heat_flux) then
                   HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) = HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) - &
                    ( 0.5d0 * species_z(species(i_atom_2)) * AS_dde_v_at_k(i_coord,i_coord2,i_atom_2) )
                 end if
               end do
             end do
           end if
 
           if(forces_on)then

              hellman_feynman_forces(1:3, i_atom_2 ) =  hellman_feynman_forces(1:3, i_atom_2) &
                   + HF_temp(1:3, i_atom_2 )

           end if

        end if
     end do

  else if (use_hartree_non_periodic_ewald) then
     do i_atom_2 = 1, n_occ_atoms
        if ( myid == task_list(i_atom_2) ) then

           ! Set the following variable to zero because the potential is _added_ to it.
           delta_v_hartree_aux = 0

           ! get current nucleus coordinate
           coord_current(:) = coords(:, i_atom_2)


           if (forces_on) then

             call interpolate_extended_hartree_non_periodic_ewald( coord_current, delta_v_hartree_aux, &
                                                                   v_hartree_gradient_temp             )

             hellman_feynman_forces(:,i_atom_2)  =  hellman_feynman_forces(:,i_atom_2)  +  &
                                             species_z(species(i_atom_2)) * v_hartree_gradient_temp(:)

           else

             call interpolate_extended_hartree_non_periodic_ewald( coord_current, delta_v_hartree_aux )

           end if

           hartree_delta_energy = hartree_delta_energy +  &
                species_z(species(i_atom_2)) * delta_v_hartree_aux

           ! CC: Save nuclear contributions to hartree_delta_energy
           ! This is the electrostatice energy of i_atom_2 in the field of all other PBC atoms 
           if (flag_energy_density) then
               ed_hartree_delta_energy_nuclei(i_atom_2) = &
                    ed_hartree_delta_energy_nuclei(i_atom_2) &
                  + species_z(species(i_atom_2)) * delta_v_hartree_aux
           end if
           if (flag_out_locpot_atom) then 
               elec_delta_atom(i_atom_2) = &
                  elec_delta_atom(i_atom_2)+ &
                  delta_v_hartree_aux
           endif
           
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

  ! MS: Add the interaction energy of the solvent potential with all nuclei
  !     after synchronization. Otherwise, it is added for every task.
   if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then
      call start_timer(timestamps_mpe_pot(1:2))
      mpe_interaction_energy_nuc = mpe_calculate_interaction_energy_with_nuclei()
      call stop_timer(timestamps_mpe_pot, unsynced=.true.)
      hartree_delta_energy = hartree_delta_energy + &
               mpe_interaction_energy_nuc
      call output_timer('MPE Reaction field evaluation', &
         timestamps_mpe_pot(3:4))

      call sync_real_number(mpe_interaction_energy_total_rho)

      call mpe_store_energy_contributions( &
         rho_V=mpe_interaction_energy_total_rho, &
         nuc_V=mpe_interaction_energy_nuc, &
         nonel=mpe_calculate_nonelectrostatic_energy() )

      call sync_vector(mpe_charge_per_dc, size(mpe_charge_per_dc))
      mpe_total_charge = sum(mpe_charge_per_dc)
      write(info_str,"(4X,A)") &
            "integrated electron density in dielectric regions"
      call localorb_info(info_str, use_unit)
      write(info_str,"(6X,A2,2X,A24,2X,A15,2X,A6)") &
            "i", "name of region", "electrons", "%"
      call localorb_info(info_str, use_unit)
      write(info_str,"(6X,A)") repeat('-', 53)
      call localorb_info(info_str, use_unit)
      do mpe_i_dc = lbound(mpe_charge_per_dc,1), ubound(mpe_charge_per_dc,1)
         if (mpe_charge_per_dc(mpe_i_dc) .eq. 0.e0_dp) cycle
         write(info_str,"(6X,I2,2X,A24,2X,F15.3,2X,F6.2)") &
            mpe_i_dc, trim(mpe_get_continuum_name(mpe_i_dc)), &
            mpe_charge_per_dc(mpe_i_dc), &
            1.e2_dp*mpe_charge_per_dc(mpe_i_dc)/mpe_total_charge
         call localorb_info(info_str, use_unit)
      enddo
      write(info_str,"(6X,A)") repeat('-', 53)
      call localorb_info(info_str, use_unit)
      write(info_str,"(6X,A28,2X,F15.3,2X,F6.2)") &
            "total integrated charge:", mpe_total_charge, 1.e2_dp
      call localorb_info(info_str, use_unit)
      write(info_str,"(6X,A28,2X,F15.3)") &
            "expected number of electrons:", n_electrons
      call localorb_info(info_str, use_unit)
      write(info_str,"(6X,A28,2X,E15.8)") &
            "discrepancy: ", n_electrons - mpe_total_charge
      call localorb_info(info_str, use_unit)
   endif
   
  if (solvent_method.eq.SOLVENT_MPB.and.mpb_solver_started) then
      call sync_sum_up_whole_potential(en_eps_rho,dummy,dummy,dummy,dummy )
      call sync_sum_up_whole_potential(en_alpha_rho,dummy,dummy,dummy,dummy )
      if (.not.(solve_lpbe_only.and..not.use_mpbe_free_energy).and..not.mpbe_no_ion) then
	call sync_sum_up_whole_potential(freeen_MPBE_nonelstat_solv_energy,dummy,dummy,dummy,dummy )
	call sync_sum_up_whole_potential(freen_solv_plus_ion,dummy,dummy,dummy,dummy )
	call sync_sum_up_whole_potential(freeen_sol_plus_solv_plus_ion,dummy,dummy,dummy,dummy )
	call sync_sum_up_whole_potential(freen_est,dummy,dummy,dummy,dummy )
	write(info_str,'(2X,A,E14.6)')  &
	      "| Nonelstat. MPBE Solvation Energy", freeen_MPBE_nonelstat_solv_energy
	call localorb_info ( info_str, use_unit,'(A)', OL_norm )
	write(info_str,'(2X,A,E14.6)')  &
	      "| Solute+Solvent+Ions Free Energy", freeen_sol_plus_solv_plus_ion
	call localorb_info ( info_str, use_unit,'(A)', OL_norm )
	write(info_str,'(2X,A,E14.6)')  &
	      "| Solvent+Ions Free Energy", freen_solv_plus_ion
	call localorb_info ( info_str, use_unit,'(A)', OL_norm )
	write(info_str,'(2X,A,E14.6)')  &
	      "| Estimate for Ionic Effect", freen_est
	call localorb_info ( info_str, use_unit,'(A)', OL_norm )
      else if (solve_lpbe_only.and..not.use_mpbe_free_energy) then
	call sync_sum_up_whole_potential(freeen_LPBE_nonelstat_solv_energy,dummy,dummy,dummy,dummy )
	write(info_str,'(2X,A,E14.6)')  &
	      "| Nonelstat. LPBE Solvation Energy", freeen_LPBE_nonelstat_solv_energy
	call localorb_info ( info_str, use_unit,'(A)', OL_norm )
      end if
      call evaluate_pot_lpb_at_zero(lpb_zero_energy) 
      if (output_level.ne.'MD_light') then
        write(use_unit,*) 'hartree_delta_energy', hartree_delta_energy
      endif
      hartree_delta_energy= hartree_delta_energy + lpb_zero_energy
      if (output_level.ne.'MD_light') then
        write(use_unit,*) 'lpb_zero_energy',lpb_zero_energy
      endif
  end if    

  if (debug_lpb_solver.and.solvent_method.eq.SOLVENT_MPB.and.&
    mpb_solver_started.and.(output_level.ne.'MD_light')) then
    write(use_unit,'(2X,A,F15.7)') '| LPB zero energy', lpb_zero_energy
  end if

  if (debug_lpb_solver.and.solvent_method.eq.SOLVENT_MPB) then
    write(use_unit,'(2X,A,F6.2)') '| Total Hartree delta energy correction', hartree_delta_energy
  end if

  ! CC: Sync nuclear energy density used to ensure the 0 shift below
  if (flag_energy_density) then
    call sync_vector(ed_hartree_delta_energy_nuclei,n_atoms)
  end if
  
  if (flag_out_locpot_atom) then 
     call sync_vector(elec_delta_atom, n_atoms)
   !  call sync_vector(elec_delta_atom_recip, n_atoms)
  endif

  ! CC: Sync stress before adding average
  if (AS_stress_on) then
    call sync_matrix( NEW_AS_at_stress,  3, 3 )
    call sync_matrix( NEW_AS_rho_stress, 3, 3 )
    call sync_matrix( AS_rho_f_f_self_correction, 3, 3 )
    if (compute_heat_flux) then
      call sync_matrix( HF_stress_per_atom_MP_AT(1,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_AT(2,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_AT(3,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_EL(1,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_EL(2,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_EL(3,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_CO(1,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_CO(2,:,:) , 3, n_atoms )
      call sync_matrix( HF_stress_per_atom_MP_CO(3,:,:) , 3, n_atoms )
    end if
    if (use_rlsy_symmetry) then
      NEW_AS_rho_stress_sum = 0.0d0
      do i_op = 1, rlsy_h%spacegroup%n_operation
        NEW_AS_rho_stress_sum(:,:) = NEW_AS_rho_stress_sum(:,:) + &
            matmul( &
            matmul( rlsy_h%spacegroup%op(i_op)%m(:,:), &
                    NEW_AS_rho_stress(:,:)), &
                    transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
      end do
      NEW_AS_rho_stress = NEW_AS_rho_stress_sum

      AS_rho_f_f_self_correction_sum = 0.0d0
      do i_op = 1, rlsy_h%spacegroup%n_operation
        AS_rho_f_f_self_correction_sum(:,:) = AS_rho_f_f_self_correction_sum(:,:) + &
            matmul( &
            matmul( rlsy_h%spacegroup%op(i_op)%m(:,:), &
                    AS_rho_f_f_self_correction(:,:)), &
                    transpose(rlsy_h%spacegroup%op(i_op)%m(:,:)))
      end do
      AS_rho_f_f_self_correction = AS_rho_f_f_self_correction_sum
    end if
  end if

  ! write root-mean square error of Hartree potential
  hartree_multipole_error =   &
       sqrt(hartree_multipole_error)

     write(info_str,'(2X,A,1X,E14.6)')  &
          "| RMS charge density error from multipole expansion :",   &
          hartree_multipole_error
     call localorb_info ( info_str, use_unit,'(A)', OL_norm )

  ! And finally, for periodic systems, we must remember to shift the potential zero 
  ! to the average real-space potential in periodic systems, also for energy-related quantities
  if (n_periodic.gt.0) then

      call sync_average_potential ( average_delta_v_hartree_real )

      average_delta_v_hartree_real = average_delta_v_hartree_real/cell_volume
      previous_average_delta_v_hartree_real = average_delta_v_hartree_real

      write (info_str,'(2X,A,1X,F15.8,A)')  &
           "| Average real-space part of the electrostatic potential :",   &
           average_delta_v_hartree_real*hartree, " eV"
      call localorb_info ( info_str, use_unit,'(A)', OL_norm )

      ! If the average electrostatic potential is really computed analytically,
      ! we replace the above value as follows:
      if (analytic_potential_average) then
         average_delta_v_hartree_real = 0.d0
         do i_atom_2 = 1,n_atoms,1
            average_delta_v_hartree_real = &
            average_delta_v_hartree_real + atom_average_es_pot(i_atom_2)
         enddo
         average_delta_v_hartree_real = average_delta_v_hartree_real / cell_volume
         previous_average_delta_v_hartree_real = average_delta_v_hartree_real
         call sync_average_potential ( average_delta_v_hartree_real )
         write (info_str,'(2X,A,1X,F15.8,A)')  &
           "| Analytical average real-space  electrostatic potential :",   &
           average_delta_v_hartree_real*hartree, " eV"
         call localorb_info ( info_str, use_unit,'(A)', OL_norm )
      end if

      !CC: Adjust local potential for energy density
      if (flag_chetty_martin_energy_density) then
        do i_index=lbound(ed_local_potential,1),ubound(ed_local_potential,1)
          ed_local_potential(i_index) = ed_local_potential(i_index) - average_delta_v_hartree_real
        end do
      end if

      potential(:) = potential(:) - average_delta_v_hartree_real  

      !if (flag_out_locpot_atom) then 
      !   elec_delta_atom_real(:) = elec_delta_atom_real(:) - average_delta_v_hartree_real
      !endif


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

       hartree_delta_energy = hartree_delta_energy & 
         -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_delta_v_hartree_real

       ! CC: Since this term has been reduced to a sum, no Jacobian term appears
       if ( AS_stress_on ) then
         do i_coord=1,3
           do i_coord2=1,3
             NEW_AS_rho_stress(i_coord,i_coord2) = NEW_AS_rho_stress(i_coord,i_coord2) + &
               ( 0.5d0*sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * AS_dde_v_rho_mp_a(i_coord,i_coord2) )
             if (compute_heat_flux) then
               HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) = HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) + &
                 ( 0.5d0*   ( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * AS_dde_v_rho_mp_a(i_coord,i_coord2) )
             end if
           end do
         end do
       end if


       ! Save nuclear contributions to hartree_delta_energy
       ! This is the electrostatice energy of i_atom_2 in the field of all other atoms (index i_center)
       if (flag_energy_density) then
         do i_atom_2=1,n_atoms,1
                ed_hartree_delta_energy_nuclei(i_atom_2) = ed_hartree_delta_energy_nuclei(i_atom_2) & 
                 & -  multipole_moments(1,i_atom_2) * sqrt(pi4) * average_delta_v_hartree_real 
         end do
       end if
       
     !  if (flag_out_locpot_atom) then 
     !     do i_atom_2 = 1, n_atoms, 1 
     !        elec_delta_atom(i_atom_2) = elec_delta_atom(i_atom_2)/species_z(species(i_atom_2))
     !     enddo
     !  endif
       
       en_elec_delta = en_elec_delta &
         -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_delta_v_hartree_real
 
       !!! hartree_delta_energy = hartree_delta_energy & 
       !!!  +  charge * average_delta_v_hartree_real
 
       ! interaction of neutralizing charge background with average_free_es_pot as 
       ! computed once and for all in initialize_grid_storage()

       ! This term should affect only the electrostatic energy of the neutralizing background, and not
       ! en_elec_delta, which is concerned only with actual electrons

        hartree_delta_energy = hartree_delta_energy & 
          -  sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * average_free_es_pot

        ! CC: see comment given above
        if ( AS_stress_on ) then
          do i_coord=1,3,1
            do i_coord2=1,3,1
              NEW_AS_rho_stress(i_coord,i_coord2) = NEW_AS_rho_stress(i_coord,i_coord2) + &
                ( 0.5d0*sum( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * ( AS_dde_v_rho_f_a(i_coord,i_coord2)  ) )
               if (compute_heat_flux) then
                 HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) = HF_stress_per_atom_MP_EL(i_coord,i_coord2,1:n_atoms) + &
                   ( 0.5d0*   ( multipole_moments(1,1:n_atoms)) * sqrt(pi4) * AS_dde_v_rho_f_a(i_coord,i_coord2) )
               end if
            end do
          end do
        end if

        ! CC: Save nuclear contributions to hartree_delta_energy
        if (flag_energy_density) then
          do i_atom_2=1,n_atoms,1
                 ed_hartree_delta_energy_nuclei(i_atom_2) = ed_hartree_delta_energy_nuclei(i_atom_2) & 
                  & -  multipole_moments(1,i_atom_2) * sqrt(pi4) * average_free_es_pot
          end do
       end if

       !!! hartree_delta_energy = hartree_delta_energy & 
       !!!   + charge * average_free_es_pot

       do i_atom_2 = 1, n_occ_atoms, 1
          ! The shift below is needed to generate a consistent total energy. 
          !
          ! For a consistent potential zero, the nuclei also need to 
          ! be shifted correctly. This is done in Hartree delta energy although
          ! the term really does not belong there.
          ! The factor 2 arises since hartree_delta_energy gets a 1/2 due to 
          ! double counting further down the road.
          hartree_delta_energy = hartree_delta_energy -  &
                2.d0 * species_z(species(i_atom_2)) * average_delta_v_hartree_real

          ! CC: Save nuclear contributions to hartree_delta_energy
          if (flag_energy_density) then
            ed_hartree_delta_energy_nuclei(i_atom_2) = ed_hartree_delta_energy_nuclei(i_atom_2) & 
             & - 2.d0 * species_z(species(i_atom_2)) * average_delta_v_hartree_real 
          end if

          ! CC: Stress:
          ! see comments given above
          if (AS_stress_on) then
            do i_coord=1,3,1
              do i_coord2=1,3,1
                NEW_AS_at_stress(i_coord,i_coord2) = NEW_AS_at_stress(i_coord,i_coord2) - &
                   ( species_z(species(i_atom_2)) * AS_dde_v_rho_mp_a(i_coord,i_coord2) )
                if (compute_heat_flux) then
                  HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) = HF_stress_per_atom_MP_AT(i_coord,i_coord2,i_atom_2) + &
                   ( species_z(species(i_atom_2)) * AS_dde_v_rho_mp_a(i_coord,i_coord2) )
                end if
              end do
            end do
          end if


          ! This accounts for the electrostatic shift of the sum of eigenvalues
          ! The factor is 1.0 as this quantity enters the expression for the kinetic energy
          ! without any prefactors (unlike hartree_delta_energy, which gets a 0.5 in the total energy)
          en_elec_delta = en_elec_delta -  &
                1.d0 * species_z(species(i_atom_2)) * average_delta_v_hartree_real

       enddo

       ! If requested, we can now write the average electrostatic potential 
       ! in the vacuum region of a slab calculation
       ! Notice that only the reciprocal-space part and the average real-space potential
       ! contribute here.
       if ( out_embedding_potential ) then
         call output_embedding_potential_z(out_embedding_potential_z1, &
               out_embedding_potential_z2, out_embedding_potential_z_grid) 
       end if

       if ( out_vacuum_potential ) then
         write(info_str,*) ' | Potential at the vacuum:'
         call localorb_info(info_str,use_unit,'(A)')
         write(info_str,*) ' z (A)  avgPot (eV)   dipPot (eV) '
         call localorb_info(info_str,use_unit,'(A)')
    

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
  
  
  
  multipole_radius_sq_public = multipole_radius_sq
  if(use_batch_permutation > 0) then

    call permute_point_array_back(n_bp,1,potential,potential_std)
    deallocate(potential)

    if (solvent_method.eq.SOLVENT_MPE .and. mpe_check_initialization()) then
       call permute_point_array_back(n_bp,dc_indices_use,mpe_dc_indices)
       deallocate(dc_indices_use)
    end if
       
    if(flag_delta_rho_in_multipole) then
      call permute_point_array_back(n_bp,1,rho_multipole_old,rho_multipole_old_std)
      deallocate(rho_multipole_old)
    endif

    deallocate(rho)
    deallocate(free_rho_superpos)
    deallocate(free_hartree_superpos)
    if (use_embedding_potential) then
       deallocate(pot_ion_embed)
    endif

  endif

  if(use_rlsy_symmetry) then


    if (forces_on) then
      multipole_forces_temp = multipole_forces
      multipole_forces = 0.0d0
      do i_op = 1, rlsy_h%spacegroup%n_operation
        do i_atom = 1, n_atoms
          j_atom = rlsy_h%spacegroup%op(i_op)%fmap(i_atom)
          multipole_forces(:,j_atom) = multipole_forces(:,j_atom) + &
              matmul(rlsy_h%spacegroup%op(i_op)%m, multipole_forces_temp(:,i_atom))
        end do
      end do
    end if

    call compute_batch_rlsy_irr_to_full_points(potential(:), potential_std(:))


    if(flag_delta_rho_in_multipole) then
      deallocate(rho_multipole_old)
    endif
    deallocate(rho)
    deallocate(free_rho_superpos)
    deallocate(free_hartree_superpos)
    if (use_embedding_potential) then
       deallocate(pot_ion_embed)
    endif
    deallocate(potential)
  endif

  ! Write one last line to bound output
  write(info_str,*) ' '
  call localorb_info(info_str,use_unit,'(A)',OL_norm)



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

  if (allocated(have_atom_average_es_pot)) then
     deallocate(have_atom_average_es_pot)
  end if
  if (allocated(atom_average_es_pot)) then
     deallocate(atom_average_es_pot)
  end if
 

  if (AS_stress_on) then
    call deallocate_all_AS_arrays()
  end if


  average_potential=previous_average_delta_v_hartree_real

end subroutine sum_up_whole_potential_p1
!******
!------------------------------------------------------------------------------------------------------------
