!****s* FHI-aims/update_density_and_forces_orbital
!  NAME
!   update_density_and_forces_orbital
!  SYNOPSIS

subroutine update_density_and_forces_orbital &
     ( KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, partition_tab,  &
     hartree_partition_tab, rho, rho_gradient, kinetic_density, hartree_potential, &
     basis_l_max, delta_rho_KS, delta_rho_gradient, delta_kinetic_density, &
     rho_change, hellman_feynman_forces, pulay_forces, gga_forces, nlcc_forces, &
     forces_on, gga_forces_on, nlcc_forces_on, meta_gga_forces_on &
     )

  !  PURPOSE
  !  Subroutine update_density obtains the new KS density from the eigenvectors
  !  and occupation numbers; no mixing.
  !
  !  For the benefit of Pulay mixing later on, we obtain the change between
  !  the input and output density, delta_rho_KS = rho_KS - rho_in,
  !  rather than overwriting rho itself.
  !
  !  This routine corresponds to Eq. (26) in Ref.
  !  Blum et al., Computer Physics Communications 180 (2009) 2175-2196
  !  which is beneficial for non-periodic systems of relatively small size.
  !
  !  Note that this routine is therefore only ever used for non-periodic systems 
  !  FHI-aims. For periodic systems, we use update_density_densmat and
  !  update_forces_densmat instead (which relies on Eq. (28) of Blum et al.,
  !  Computer Physics Communications 180 (2009) 2175-2196).
  !
  !  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use basis
  use cartesian_ylm
  use constants
  use KH_core_states
  use pseudodata
  use energy_density
  use lpb_solver_utilities, only: evaluate_sc_Gnonmf, evaluate_Gnonmf, &
    evaluate_Gnonmf_energy_shell, en_pot_Gnonmf, mpb_solver_started, &
    en_Gnonmf, surface_and_volume_calc, surface_mpb, volume_mpb,&
    alphaplusgamma_mpb, beta_mpb, evaluate_cavity_surface_and_volume,&
    Gnonmf_forces, evaluate_Gnonmf_forces, Gnonmf_forces_on,&
    evaluate_orb_hess_dot_rho_grad_Gnonmf, atomic_MERM, mpb_solver_started
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use aims_gpu, only: gpu_forces_used, gpu_density_used
  use timing, only: gpu_forces_always_used, gpu_density_always_used
  implicit none

  !  ARGUMENTS


  !  input 

  real*8,     dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points), intent(IN) :: KS_eigenvector_complex

  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: occ_numbers
  real*8, dimension(n_states, n_spin, n_k_points), intent(IN) :: KS_eigenvalue      

  real*8, dimension(n_full_points) ::  partition_tab
  real*8, dimension(n_full_points) ::  hartree_partition_tab
  !     rho is only input; what we store is the density residual (i.e. the change in the density)
  real*8, dimension(n_spin, n_full_points), intent(IN) :: rho
  real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
  real*8, dimension(n_spin, n_full_points), intent(IN) :: kinetic_density
  real*8, dimension(n_full_points) :: hartree_potential

  integer basis_l_max (n_species)
  logical :: forces_on
  logical :: gga_forces_on
  logical :: nlcc_forces_on
  logical :: meta_gga_forces_on

  !  output

  real*8, intent(OUT) :: delta_rho_KS(n_full_points, n_spin)
  real*8, intent(OUT) :: delta_rho_gradient(3, n_full_points, n_spin)
  real*8, intent(OUT) :: delta_kinetic_density(n_full_points, n_spin)
  real*8, dimension(3, n_atoms) :: hellman_feynman_forces
  real*8, dimension(3, n_atoms) :: pulay_forces
  real*8, dimension(3, n_atoms) :: gga_forces
  real*8, dimension(3, n_atoms) :: nlcc_forces
  real*8, dimension(n_spin) :: rho_change

  !  INPUTS
  !   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
  !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
  !   o occ_numbers -- occupations of eigenstates
  !   o KS_eigenvalue -- Kohn-Sham eigenvalues
  !   o partition_tab -- values of partition function
  !   o hartree_partition_tab -- values of partition function used in multipole expansion of charges in Hartree potential.
  !   o rho -- electron density, 
  !            rho is only input; what we store is the density residual (i.e. the change in the density)
  !   o rho_gradient -- gradient of electron density
  !   o kinetic_density -- kinetic-energy density
  !   o hartree_potential -- Hartree potential
  !   o basis_l_max -- maximum l of basis functions
  !   o forces_on -- are we calculating forces or not
  !   o gga_forces_on -- are we calculating gga forces or not
  !   o nlcc_forces_on -- are we calculating nlcc forces or not 
  !   o meta_gga_forces_on -- are we calculating meta_gga forces or not
  !
  !  OUTPUT
  !   o delta_rho_KS -- calculated charge density  -  old rho (input)
  !   o delta_rho_gradient -- calculated gradient  -  old rho_gradient (input)
  !   o hellman_feynman_forces -- Hellman Feyman forces
  !                               NOTE:  This subroutine no longer calculates
  !                               Hellman-Feynman forces by default.  This is
  !                               done in sum_up_whole_potential_p1 instead.
  !   o pulay_forces -- Pulay forces
  !   o gga_forces -- gga forces
  !   o nlcc_forces -- nlcc forces
  !   o rho_change -- change of charge density
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



  !  local variables

  real*8 coord_current(3)
  real*8 dist_tab(n_centers_basis_integrals  , n_max_batch_size)
  real*8 dist_tab_sq(n_centers_basis_integrals, n_max_batch_size)
  real*8 i_r(n_centers_basis_integrals)
  real*8 dir_tab(3, n_centers_basis_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_basis_integrals)

  !  real*8 radial_wave(n_centers_basis_T, n_max_batch_size)
  !  real*8 radial_wave_deriv(n_centers_basis_T, n_max_batch_size)
  !  real*8 radial_wave_2nd_deriv(n_centers_basis_T, n_max_batch_size)
  !  real*8 wave(n_centers_basis_T, n_max_batch_size)
  !  real*8 kinetic_wave(n_centers_basis_T, n_max_batch_size)

  real*8,dimension(:),allocatable:: radial_wave
  real*8,dimension(:),allocatable:: radial_wave_deriv
  real*8,dimension(:),allocatable:: radial_wave_2nd_deriv
  real*8,dimension(:,:),allocatable:: wave
  real*8,dimension(:),allocatable:: kinetic_wave
  real*8,dimension(:)  ,allocatable:: kinetic_wave_deriv
  real*8, dimension(:,:,:), allocatable :: kinetic_gradient_basis_wave


  real*8 :: local_potential_parts(n_spin)
  real*8 :: xc_gradient_deriv(3, n_spin, n_max_batch_size)
  real*8 :: xc_tau_deriv(n_spin, n_max_batch_size)
  real*8 :: en_density_xc(n_max_batch_size)
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8 :: local_xc_derivs(n_spin, n_max_batch_size)
  real*8 :: nlcc_local_xc_derivs(n_spin, n_max_batch_size)
  real*8 :: nlcc_xc_gradient_deriv(3, n_spin, n_max_batch_size)

  ! allocate only in case of force-calculations with gga
  real*8, dimension(:,:,:), allocatable :: cartesians

  !!VB: Order of sum_gradient changed in _p2 version
  real*8, dimension(:,:,:), allocatable :: sum_gradient
  real*8, dimension(:,:,:), allocatable :: sum_hessian
  real*8, dimension(:,:,:), allocatable :: orb_grad_dot_rho_grad 
  real*8, dimension(:,:,:,:), allocatable :: &
       orb_hess_dot_rho_grad

  ! mGGA temporary storage
  real*8, dimension(:,:,:), allocatable :: orb_hess_dot_orb_grad

  !     pruning of atoms, radial functions, and basis functions

  integer :: n_compute_c
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_basis_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_basis_integrals)
  integer :: i_center
  integer :: n_compute_atoms
  integer :: atom_index(n_centers_basis_integrals)
  integer :: atom_index_inv(n_centers)

  ! pruning for force calculation
  integer :: global_atom(n_atoms)
  integer :: basis_offset(n_atoms+1)

  integer :: spline_array_start(n_centers_basis_integrals)
  integer :: spline_array_end(n_centers_basis_integrals)

  ! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_dens)
  integer :: l_index(n_max_compute_fns_dens)
  integer :: l_count(n_max_compute_fns_dens)
  integer :: fn_atom(n_max_compute_fns_dens)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_dens)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  !     other local variables
  integer, dimension(n_spin) :: max_occ_number
  real*8 :: occ_numbers_sqrt

  integer :: l_ylm_max
  integer :: n_compute_force_atoms
  integer :: n_local_compute
  integer :: n_points

  logical :: write_out

  real*8,     dimension(:,:,:),allocatable :: KS_vec_times_occ_sqrt
  complex*16, dimension(:,:,:),allocatable :: KS_vec_times_occ_sqrt_complex



  ! (n_basis, n_states, n_spin) on purpose
  ! beware: inside the subroutines (evaluate_KS_density_v1 for instance)
  ! dimension is (n_compute, max_occ_number)
  ! that's a trick to get a continuous data flow and not
  ! bad programming
  real*8, dimension(n_species) :: r_grid_min_sq

  real*8,     dimension(:,:,:),allocatable :: KS_ev_compute
  complex*16, dimension(:,:,:),allocatable :: KS_ev_compute_complex

  real*8,    dimension(:,:,:),allocatable :: KS_orbital
  complex*16,dimension(:,:,:),allocatable :: KS_orbital_complex

  real*8 temp_rho(n_max_batch_size,n_spin)
  real*8,dimension(:,:),allocatable :: temp_rho_small

  real*8 temp_rho_gradient(3,n_max_batch_size,n_spin)
  real*8 temp_kinetic_density(n_max_batch_size,n_spin)
  real*8 partition(n_max_batch_size)

  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:,:), allocatable :: dir_tab_global
  real*8, dimension(:,:), allocatable :: dist_tab_sq_global

  !     only allocated and referenced for gradient functionals
  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave
  real*8, dimension(:,:,:), allocatable :: hessian_basis_wave
  real*8, dimension(:,:,:,:), allocatable :: KS_orbital_gradient
  complex*16, dimension(:,:,:,:), allocatable :: KS_orbital_gradient_complex
  real*8, dimension(:,:,:,:), allocatable :: nuclear_gradient
  real*8, dimension(:,:,:),   allocatable :: h_times_wave
  real*8, dimension(:,:,:),   allocatable ::h_minus_e_times_psi

  real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
  real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore


!  real*8, dimension(:,:), allocatable :: partial_core_rho_grad_in_batch

  !     counters

  integer i_atom
  integer i_l
  integer i_m
  integer i_coord
  integer :: i_state
  integer :: i_point
  integer :: i_bas

  integer :: i_my_batch

  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_3

  ! i_spin for future use in a spin-polarized version
  integer :: i_spin = 1

  integer ::  i_index, i_k_point, i_k_point_group

  !     mpi

  !      integer, dimension(:,:), allocatable :: n_points_mpi
  integer :: i_atom_2
  integer :: info
  integer:: n_k_group, i_k_point_g

  real*8 :: dummy  

  !//MPB solvation
  real*8, dimension(:), allocatable :: local_Gnonmf_derivs
  real*8, dimension(:,:), allocatable :: Gnonmf_gradient_deriv  
  real*8, dimension(:,:,:), allocatable :: orb_grad_dot_rho_grad_Gnonmf 
  real*8, dimension(:,:,:,:), allocatable :: orb_hess_dot_rho_grad_Gnonmf

  ! mGGA temporary storage
  real*8, dimension(:,:,:), allocatable :: orb_hess_dot_orb_grad_Gnonmf
  
  logical :: gga_or_Gnonmf, updating_non_partition_point, point_on_atom

  character*200 :: info_str

  !     begin work

  ! GPU acceleration is only implemented for density-matrix-based density update
  if (use_gpu_density .or. use_gpu_forces) then
    write(info_str,'(2X,A)') "Since orbital-based density update is being used,&
         &GPU acceleration will not be employed for density or forces."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
  end if
  gpu_density_used = .false.
  gpu_forces_used = .false.
  gpu_density_always_used = .false.
  gpu_forces_always_used = .false.

  if (gga_forces_on.or.Gnonmf_forces_on) then
    gga_or_Gnonmf = .true.
  else
    gga_or_Gnonmf = .false.
  end if

  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  if (Gnonmf_forces_on) then
    allocate( local_Gnonmf_derivs(n_max_batch_size),stat=info)
    call check_allocation(info, 'local_Gnonmf_derivs               ')
    allocate( Gnonmf_gradient_deriv(3,n_max_batch_size),stat=info)
    call check_allocation(info, 'Gnonmf_gradient_deriv             ')
  end if  
!  if(use_embedding_pp.and.use_nonlinear_core) then
!     if(.not.(allocated(partial_core_rho_grad_in_batch))) then
!        allocate(partial_core_rho_grad_in_batch(3,n_max_batch_size))
!     endif
!  endif

!  if(use_embedding_pp.and.use_nonlinear_core) then
  if(use_embedding_pp.and.nlcc_forces_on) then
      if(.not. allocated(rho_inc_partialcore)) allocate(rho_inc_partialcore(n_spin,n_full_points))
      if(.not. allocated(rho_gradient_inc_partialcore)) allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
      do i_spin = 1,n_spin
         rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
         rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
         if(use_density_gradient) then
            rho_gradient_inc_partialcore(:,i_spin,:) = &
               rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
         endif
      enddo
      ! Not implemented for meta-GGAs as the PSPs do not contain the necessary information, currently
  endif



  if (.not.forces_on) then
     call localorb_info( &
          "Evaluating new KS density.", use_unit,'(2X,A)', OL_norm )
  else
     call localorb_info( &
          "Evaluating new KS density and force components.",  &
          use_unit,'(2X,A)', OL_norm )
  end if

  ! CC: Skip check_occs in the case of energy density calculations
  !     The computation of the EEV density requires the occs to 
  !     be unequal the number of electrons.
  !     See comments in calling scf_solver.f90
  if ( .not. flag_harris_foulkes_energy_density ) then
    call check_occs('update_density_and_forces_orbital', occ_numbers, .true.)
  end if

  if(use_small_component)then
     allocate(temp_rho_small(n_max_batch_size,n_spin),stat=info)
     call check_allocation(info, 'temp_rho_small                ') 
  end if

  if(real_eigenvectors)then

     if(.not. allocated( KS_vec_times_occ_sqrt))then
        allocate( KS_vec_times_occ_sqrt(n_centers_basis_T,(n_states*n_k_points_group), n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_vec_times_occ_sqrt'
           stop
        end if
     end if

     if(.not. allocated( KS_ev_compute))then
        allocate( KS_ev_compute(n_states*n_k_points_group, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute'
           stop
        end if
     end if

     if(.not. allocated( KS_orbital))then
        allocate( KS_orbital(n_states*n_k_points_group,n_max_batch_size,n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_orbital'
           stop
        end if
     end if

  else

     if(.not. allocated( KS_vec_times_occ_sqrt_complex))then
        allocate( KS_vec_times_occ_sqrt_complex(n_centers_basis_T, (n_states*n_k_points_group), n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation:  KS_vec_times_occ_sqrt_complex'
           stop
        end if
     end if

     if(.not. allocated( KS_ev_compute_complex))then
        allocate( KS_ev_compute_complex(n_states*n_k_points_group, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute_complex'
           stop
        end if
     end if

     if(.not. allocated( KS_orbital_complex))then
        allocate( KS_orbital_complex(n_states*n_k_points_group,n_max_batch_size,n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_orbital_complex'
           stop
        end if
     end if
  end if

  if(.not. allocated(i_basis))then
     allocate(i_basis(n_centers_basis_T),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: i_basis'
        stop
     end if
  end if

  if(.not.allocated(radial_wave))then
     allocate(radial_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
  end if

  if(.not. allocated(radial_wave_deriv))then
     allocate(radial_wave_deriv(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave_deriv'
        stop
     end if
  end if

  if(.not. allocated(radial_wave_2nd_deriv))then
     allocate(radial_wave_2nd_deriv(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave_2nd_deriv'
        stop
     end if
  end if

  if(.not. allocated(wave))then
     call aims_allocate( wave, n_max_compute_dens, n_max_batch_size, "wave")
  end if

  if(.not. allocated(kinetic_wave))then
     allocate(kinetic_wave(n_max_compute_fns_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: kinetic_wave'
        stop
     end if
  end if

  l_ylm_max = l_wave_max

  if(.not. allocated( ylm_tab))then
     allocate( ylm_tab( (l_ylm_max+1)**2,n_centers_basis_integrals),stat=info )
     if(info/=0)then
        write(use_unit,*)'Error in allocation: ylm_tab'
        stop
     end if
  end if


  if(.not. allocated( index_lm))then
     allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info) 
     if(info/=0)then
        write(use_unit,*)'Error in allocation: index_lm'
        stop
     end if
  end if


  if (use_density_gradient .or. forces_on) then
     !     allocate local arrays needed for gradients
     call aims_allocate ( gradient_basis_wave, n_max_compute_dens, 3, n_max_batch_size, "gradient_basis_wave" )

     if (.not. gga_forces_on.and..not.Gnonmf_forces_on) then 
        allocate( dylm_dtheta_tab((l_ylm_max+1)**2, n_centers_basis_integrals),stat=info )
        if(info/=0)then
           write(use_unit,*)'Error in allocation: dylm_dtheta_tab'
           stop
        end if
        allocate( scaled_dylm_dphi_tab((l_ylm_max+1)**2,n_centers_basis_integrals),stat=info )
        if(info/=0)then
           write(use_unit,*)'Error in allocation: scaled_dylm_dphi_tab'
           stop
        end if
     end if
  end if

  if (.not.allocated(KS_orbital_gradient)) then
     if ( use_density_gradient .and. (real_eigenvectors) ) then
           allocate(KS_orbital_gradient(n_states*n_k_points_group,n_max_batch_size,3,n_spin),stat=info)
           if(info/=0)then
              write(use_unit,*)'Error in allocation: KS_orbital_gradient'
              stop
           end if
           KS_orbital_gradient = 0.d0

     else ! dummy allocation - needed for either spin channel
           allocate(KS_orbital_gradient(1,1,1,n_spin),stat=info)
     end if
  end if
  
  ! Since this routine is never used for periodic boundary conditions in practice, the complex infrastructure
  ! is never presently used. Note that this could change if we were to introduce magnetic fields in the future.
  if (.not.allocated(KS_orbital_gradient_complex)) then
     if ( use_density_gradient .and. (.not.real_eigenvectors) ) then
        allocate(KS_orbital_gradient_complex(n_states*n_k_points_group,n_max_batch_size,3,n_spin),stat=info)
        if(info/=0)then
              write(use_unit,*)'Error in allocation: KS_orbital_gradient_complex'
              stop
        end if
           KS_orbital_gradient_complex = (0.d0,0.d0)
     else ! dummy allocation
        allocate(KS_orbital_gradient_complex(1,1,1,n_spin),stat=info)
     end if
  end if

  if (forces_on) then
     if (.not.allocated(h_times_wave)) then
        call aims_allocate( h_times_wave, n_basis, n_max_batch_size, n_spin, "h_times_wave" )
     end if
     if (.not.allocated(dir_tab_global)) then
        allocate(dir_tab_global(3, n_atoms, n_max_batch_size))
     end if
     if (.not.allocated(dist_tab_sq_global)) then
        allocate(dist_tab_sq_global(n_atoms, n_max_batch_size))
     end if
     if (.not.allocated(nuclear_gradient)) then
        allocate(nuclear_gradient &
             (n_states,n_max_batch_size,3,n_spin))
     end if
     if (.not.allocated(h_minus_e_times_psi)) then
        allocate(h_minus_e_times_psi(n_states,n_max_batch_size,n_spin))
     end if
     if (gga_forces_on.or.Gnonmf_forces_on) then
        if (.not.allocated(cartesians)) then
           allocate(cartesians(n_max_cartesian, 0:l_wave_max, n_atoms))
        end if
        if (.not.allocated(sum_gradient)) then
           allocate(sum_gradient((l_wave_max+1) ** 2, 3, n_atoms))
        end if
        if (.not.allocated(sum_hessian)) then
           allocate(sum_hessian(6, (l_wave_max+1) ** 2, n_atoms))
        end if
        if (.not.allocated(hessian_basis_wave)) then
           call aims_allocate( hessian_basis_wave, n_max_compute_dens, 6, n_max_batch_size, "hessian_basis_wave" )
        end if
        if (gga_forces_on) then
	  if (.not.allocated(orb_grad_dot_rho_grad)) then
	    allocate(orb_grad_dot_rho_grad &
		  (n_states, n_max_batch_size, n_spin))
	  end if
	  if (.not.allocated(orb_hess_dot_rho_grad)) then
	    allocate(orb_hess_dot_rho_grad &
		  (n_states, n_max_batch_size, 3, n_spin))
	  end if
	  ! As I pass this into the subroutines this needs to be
	  ! explicitly defined, and not dependent on meta_gga_forces_on. AJL
	  if (.not.allocated(orb_hess_dot_orb_grad)) then
	    allocate(orb_hess_dot_orb_grad &
		  (3, n_max_batch_size, n_spin))
	  end if
	end if
	if (Gnonmf_forces_on) then
	  if (.not.allocated(orb_grad_dot_rho_grad_Gnonmf)) then
	    allocate(orb_grad_dot_rho_grad_Gnonmf &
		  (n_states, n_max_batch_size, n_spin))
	  end if
	  if (.not.allocated(orb_hess_dot_rho_grad_Gnonmf)) then
	    allocate(orb_hess_dot_rho_grad_Gnonmf &
		  (n_states, n_max_batch_size, 3, n_spin))
	  end if
	  if (.not.allocated(orb_hess_dot_orb_grad_Gnonmf)) then
	    allocate(orb_hess_dot_orb_grad_Gnonmf &
		  (3, n_max_batch_size, n_spin))
	  end if	  
	end if
     end if !gga or Gnonmf
  end if !force calculation

  if(flag_rel.eq.REL_atomic_zora)then
     allocate (kinetic_gradient_basis_wave(n_max_compute_dens,3,n_max_batch_size ))
     allocate(kinetic_wave_deriv(n_max_compute_fns_dens))
  end if

  !     initialize index_lm
  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

  if (use_forces) then
     do i_atom = 1, n_atoms, 1
        do i_coord = 1, 3, 1

           if(.not. force_new_functional) hellman_feynman_forces(i_coord, i_atom) = 0.d0
           pulay_forces          (i_coord, i_atom) = 0.d0
           if (use_gga) then
              gga_forces          (i_coord, i_atom) = 0.d0
           end if
           if (Gnonmf_forces_on) then
	      Gnonmf_forces 	(i_coord,i_atom) = 0d0
           end if
           if (use_nonlinear_core) then
              nlcc_forces          (i_coord, i_atom) = 0.d0
           end if
        end do
     end do
  end if

  !     initialize charge density convergence criterion
  rho_change = 0.d0
  delta_rho_KS = 0.d0

  if (use_density_gradient) then
     delta_rho_gradient = 0.0d0
     if (use_meta_gga) then
       delta_kinetic_density = 0.0d0
     end if
  end if

  !     find the maximal occupation number
  n_k_group = ceiling(real(n_k_points)/real(n_k_points_group))

  !  write(use_unit,*) n_k_group,n_k_points_group

  do i_k_point_group = 1,n_k_group

     if(real_eigenvectors)then

        do i_spin = 1, n_spin, 1

           max_occ_number(i_spin) = 0
           i_index = 0


           do i_k_point_g = 1, n_k_points_group

              i_k_point = (i_k_point_group-1) * n_k_points_group + i_k_point_g

              if(i_k_point <= n_k_points)then
                 do i_state = 1, n_states, 1
                    if (occ_numbers(i_state,i_spin, i_k_point).gt.0.d0) then

                       i_index = i_index + 1

                       occ_numbers_sqrt =  sqrt(occ_numbers(i_state,i_spin,i_k_point))


                       do i_bas = 1,  n_centers_basis_T, 1



                          KS_vec_times_occ_sqrt(i_bas,i_index,i_spin) =  &
                               KS_eigenvector(Cbasis_to_basis(i_bas),i_state,i_spin, i_k_point) * &
                               occ_numbers_sqrt * &
                               dble(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),i_k_point))


                       end do
                    end if
                 end do
              end if
           enddo
           max_occ_number(i_spin) = i_index        
        enddo

     else

        do i_spin = 1, n_spin, 1


           max_occ_number(i_spin) = 0
           i_index = 0

           do i_k_point_g = 1, n_k_points_group

              i_k_point = (i_k_point_group-1) * n_k_points_group + i_k_point_g

              if(i_k_point <= n_k_points)then

                 do i_state = 1, n_states, 1
                    if (occ_numbers(i_state,i_spin, i_k_point).gt.0.d0) then

                       i_index = i_index + 1

                       occ_numbers_sqrt =  sqrt(occ_numbers(i_state,i_spin,i_k_point))


                       do i_bas = 1,  n_centers_basis_T, 1



                          KS_vec_times_occ_sqrt_complex(i_bas,i_index,i_spin) =  &
                               KS_eigenvector_complex(Cbasis_to_basis(i_bas),i_state,i_spin, i_k_point) * &
                               occ_numbers_sqrt * &
                               dconjg(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),i_k_point))


                       end do
                    end if
                 end do
              end if
           end do
           max_occ_number(i_spin) = i_index

        enddo
     end if

     !     loop over integration grid

     i_basis_fns_inv = 0


     i_full_points = 0
     i_full_points_2 = 0
     i_full_points_3 = 0

     write_out = .false.
     do i_my_batch = 1, n_my_batches, 1

           n_compute_c = 0
           i_basis = 0

           i_point = 0


           ! loop over one batch
           do i_index = 1, batches(i_my_batch)%size, 1

              i_full_points = i_full_points + 1

              updating_non_partition_point = .False.
              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).le.0.d0) then
                   updating_non_partition_point = .True.
              end if
              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then

                 i_point = i_point+1


                 !     get current integration point coordinate
                 coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

                 if(n_periodic > 0)then
                    call map_to_center_cell(coord_current(1:3) )
                 end if

                 ! compute atom-centered coordinates of current integration point,
                 ! as viewed from all atoms
                 call tab_atom_centered_coords_p0 &
                      ( coord_current,  &
                      dist_tab_sq(1,i_point),  &
                      dir_tab(1,1,i_point), &
                      n_centers_basis_integrals, centers_basis_integrals )

                 ! in case of force-calculations global geometry information is required 
                 ! (for the hellman-feynman part)

                 if (forces_on) then !.and..not.updating_non_partition_point) then
                    dist_tab_sq_global(:,i_point) = dist_tab_sq(:,i_point)
                    dir_tab_global(:,:,i_point) = dir_tab(:,:,i_point)
                    partition(i_point) = partition_tab(i_full_points)
                 end if

                 !     determine which basis functions are relevant at current integration point,
                 !     and tabulate their indices

                 ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                 if (.not.prune_basis_once .or. (atomic_MERM .and. mpb_solver_started) ) then
                    call prune_basis_p2 &
                         ( dist_tab_sq(1,i_point), &
                         n_compute_c, i_basis,  &
                         n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals  )
                 end if
              end if
           enddo ! end loop over the angular shell
           if (prune_basis_once .and. .not. (atomic_MERM .and. mpb_solver_started) ) then
              n_compute_c = batches(i_my_batch)%batch_n_compute
              i_basis(1:n_compute_c) = batches(i_my_batch)%batch_i_basis
           end if

           ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
           call collect_batch_centers_p2 &
                ( n_compute_c, i_basis, n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals, &
                n_batch_centers, batch_center &
                )

           n_points = i_point

           if (n_compute_c.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0
              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1

                 updating_non_partition_point = .False.
                 if (max(partition_tab(i_full_points_2),&
                      hartree_partition_tab(i_full_points_2)).le.0.d0) then
                      updating_non_partition_point = .True.
                 end if
                 
                 if (max(partition_tab(i_full_points_2),&
                      hartree_partition_tab(i_full_points_2)).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then

                    i_point = i_point+1
                    n_compute_atoms = 0
                    n_compute_fns = 0
                    
                    point_on_atom = .false.

                    if ((atomic_MERM .and. mpb_solver_started).and.updating_non_partition_point) then 
                      !in this case we calculate points which are usually not calculated. we need
                      !to avoid division by zero (cf. update_missing_density comments) check if any grid point lies on atom
                    
                        do i_center = 1, n_centers_integrals, 1
                            if ( dist_tab_sq(i_center,i_point).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
                                point_on_atom = .true.
                                exit ! exit do loop
                            end if
                        end do
                     end if
                     if (point_on_atom) then
                    ! set waves (the only quantity to be computed) to zero, density not needed!

                        wave(1:n_compute_c,i_point) = 0.d0

                     else
                     
                    
                        ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                        ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                        ! without any copying and without doing any unnecessary operations. 
                        ! The price is that the interface is no longer explicit in terms of physical 
                        ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                        !                     dir_tab(:,:,i_point,i_division) = dir_tab_full(:,:,i_angular)
                        call prune_radial_basis_p2 &
                            ( n_max_compute_atoms, n_max_compute_fns_dens, &
                            dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                            n_compute_atoms, atom_index, atom_index_inv, &
                            n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                            i_atom_fns, spline_array_start, spline_array_end, &
                            n_centers_basis_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                            n_batch_centers, batch_center, &
                            one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                            fn_atom, n_zero_compute, zero_index_point &
                            )


                        ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                        ! for all atoms which are actually relevant
                        call tab_local_geometry_p2 &
                            ( n_compute_atoms, atom_index, &
                            dist_tab(1,i_point),  &
                            i_r &
                            )


                        ! Determine all needed radial functions from efficient splines

                        ! Now evaluate radial functions u(r) from the previously stored compressed 
                        ! spline arrays  
                        call evaluate_radial_functions_p0 &
                            (   spline_array_start, spline_array_end, &
                            n_compute_atoms, n_compute_fns,  &
                            dist_tab(1,i_point), i_r(1), &
                            atom_index, i_basis_fns_inv, &
                            basis_wave_ordered, radial_wave(1), &
                            .false., n_compute_c, n_max_compute_fns_dens   &
                            )
                      end if



                    ! for forces or density gradient, radial derivatives are required 
                    if (use_density_gradient .or. forces_on) then
                       call evaluate_radial_functions_p0 &
                            (spline_array_start, spline_array_end, &
                            n_compute_atoms, n_compute_fns,  &
                            dist_tab(1,i_point), i_r(1), &
                            atom_index, i_basis_fns_inv, &
                            basis_deriv_ordered,  &
                            radial_wave_deriv(1),  &
                            .true., n_compute_c, n_max_compute_fns_dens  &
                            )
                    end if

                    ! for the gga-force-correction term, hessians of basis-functions are required
                    ! => 2nd derivative of radial basis function is calculated by spline-derivative
                    ! of first analytical derivative which is splined on the logarithmic grid
                    ! (not really consistent, 2nd derivative of spline polynomial of radial function 
                    ! should be evaluated, but since the logarithmic grid is dense enough, there is 
                    ! basically no difference)
                    if (forces_on .and. (gga_forces_on.or.Gnonmf_forces_on) .and..not.updating_non_partition_point) then

                       call evaluate_radial_functions_deriv_p0 &
                            (spline_array_start, spline_array_end, &
                            n_compute_atoms, n_compute_fns,  &
                            dist_tab(1,i_point), i_r(1), &
                            atom_index, i_basis_fns_inv, &
                            basis_deriv_ordered,  &
                            radial_wave_2nd_deriv(1),  &
                            .true., n_max_compute_fns_dens &
                            )

                    end if

                    !     Determine all needed ylm functions, and from these
                    !     * wave functions
                    !     * wave function derivatives
                    !     * wave function hessians

                    ! What exactly is needed depends heavily on what we are doing - forces or not, 
                    ! gga or not, etc.
                    if (forces_on .and. (gga_forces_on.or.Gnonmf_forces_on)) then


                       ! Basis functions, their gradients and hessians are needed.
                       ! They are evaluated on the basis of a cartesian expansion of the ylm-functions -
                       ! hence, if we do need all this, we get a lot of stuff directly that we'd have to
                       ! compute separately otherwise (see below)

                       ! first, corresponding cartesian terms x^l_x*y^l_y*z^l_z are evaluated
                       ! for the current integration point 
                       call evaluate_cartesians &
                            (dir_tab(1,1,i_point), l_shell_max, l_wave_max,  &
                            atom_index, n_compute_atoms, cartesians(1,0,1) )

                       ! further "ingredients" of the hessian are evaluated based on the cartesians
                       ! ylm's are evaluated on the fly
                       call evaluate_cartesian_hessian_and_gradient_terms_p2 &
                            (l_shell_max, l_wave_max, cartesians(1,0,1), &
                            atom_index, n_compute_atoms,  &
                            ylm_tab,  &
                            sum_gradient,  &
                            sum_hessian(1,1,1))

                       ! now, all hessians of the basis functions are evaluated based on the
                       ! "ingredients" (ylm, sum_two_gradient, sum_hessian)
                       ! (gradient and ylm-functions themselves would be for free, we have
                       ! to check the performance)
                       call evaluate_wave_gradient_cartesian_p2 &
                            ( n_compute_c, n_compute_atoms, n_compute_fns, &
                            one_over_dist_tab, dir_tab(1,1,i_point), l_wave_max,  &
                            ylm_tab, radial_wave, radial_wave_deriv, &
                            sum_gradient, gradient_basis_wave(1,1,i_point), &
                            rad_index, wave_index, l_index, l_count, fn_atom, &
                            n_zero_compute, zero_index_point,n_max_compute_dens &
                            )

                       if (.not.updating_non_partition_point) then
                           ! ugh ... the hessian itself.
                           call evaluate_wave_hessian_cartesian_p2 &
                                (dist_tab(1,i_point), i_r(1),  &
                                dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                                n_compute_c, i_basis, atom_index_inv,  &
                                i_basis_fns_inv, radial_wave(1),  &
                                radial_wave_deriv(1),  &
                                radial_wave_2nd_deriv(1),  &
                                ylm_tab,  &
                                sum_gradient,  &
                                sum_hessian(1,1,1), &
                                hessian_basis_wave(1,1,i_point), n_compute_atoms)
                       end if

                    else if (.not. point_on_atom .and. .not. &
                          (forces_on.and.(gga_forces_on.or.Gnonmf_forces_on).and.updating_non_partition_point)) then
                       !if forces are off 
                       ! if (.not.(gga_forces_on .and. forces_on)), then we compute the needed
                       ! ylm pieces separately from the cartesian evaluation that is done
                       ! for the hessians ...

                       ! compute trigonometric functions of spherical coordinate angles
                       ! of current integration point, viewed from all atoms

                       call tab_trigonom_p0 &
                            ( n_compute_atoms, dir_tab(1,1,i_point),  &
                            trigonom_tab(1,1) &
                            )

                       if (use_density_gradient .or. forces_on) then
                          ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                          call tab_gradient_ylm_p0 &
                               (trigonom_tab(1,1), basis_l_max,  &
                               l_ylm_max, n_compute_atoms, atom_index, &
                               ylm_tab(1,1),  &
                               dylm_dtheta_tab(1,1),  &
                               scaled_dylm_dphi_tab(1,1) &
                               )
                            
                          ! evaluate the wave function gradient directly here                  
                          call evaluate_wave_gradient_p2  &
                               ( n_compute_c, n_compute_atoms, n_compute_fns, &
                               one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                               l_ylm_max, ylm_tab,  &
                               dylm_dtheta_tab,  &
                               scaled_dylm_dphi_tab,  &
                               radial_wave,  &
                               radial_wave_deriv,  &
                               gradient_basis_wave (1:n_compute_c,1:3,i_point),  &
                               rad_index, wave_index, l_index, l_count, fn_atom, &
                               n_zero_compute, zero_index_point &
                               )

                       else

                          ! tabulate distance and Ylm's w.r.t. other atoms            
                          call tab_wave_ylm_p0 &
                               ( n_compute_atoms, atom_index,  &
                               trigonom_tab(1,1), basis_l_max,  &
                               l_ylm_max, &
                               ylm_tab(1,1) )

                       end if
                    end if !forces_on
                    ! tabulate total wave function value for each basis function in all cases -
                    ! but only now are we sure that we have ylm_tab ...

                    ! tabulate total wave function value for each basis function
                    if (.not. point_on_atom) then
                        call evaluate_waves_p2  &
                            ( n_compute_c, n_compute_atoms, n_compute_fns, &
                            l_ylm_max, ylm_tab, one_over_dist_tab,   &
                            radial_wave, wave(1,i_point), &
                            rad_index, wave_index, l_index, l_count, fn_atom, &
                            n_zero_compute, zero_index_point &
                            )
                    end if

                    !     Finally, must re-evaluate the Hamiltonian at current point if forces
                    !     are needed (for straightforward Pulay force terms ...)
                    if (forces_on) then ! .and..not.updating_non_partition_point) then

                       call evaluate_radial_functions_p0  &
                            ( spline_array_start, spline_array_end,  &
                            n_compute_atoms, n_compute_fns,   &
                            dist_tab(1,i_point), i_r(1),  &
                            atom_index, i_basis_fns_inv,  &
                            basis_kinetic_ordered, kinetic_wave(1),  &
                            .false., n_compute_c, n_max_compute_fns_dens &
                            )


                       if(flag_rel == REL_atomic_zora)then
                          kinetic_wave = kinetic_wave *0.5
                          ! write(use_unit,*) kinetic_wave(1) 

                       end if
                       

                       coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
                       if(n_periodic > 0)then
                          call map_to_center_cell(coord_current(1:3) )
                       end if
                       if(nlcc_forces_on) then
! 210113 DB: adding the partial core density here does not change the total energy but only the pulay + GGA force

!                           rho_inc_partialcore(i_spin,i_full_points_2) = &
!                              rho(i_spin,i_full_points_2) + partial_core_rho(i_full_points_2)


!                           if(use_density_gradient) then 
!                              rho_gradient_inc_partialcore(1:3,i_spin,i_full_points_2) = &
!                                rho_gradient(1:3,i_spin,i_full_points_2) + & 
!                                partial_core_rho_grad(1:3,i_full_points_2)
!                           endif

! TODO: check whether i_full_points instead of i_point would matter
!                           partial_core_rho_grad_in_batch(1:3,i_point)  = &
!                              partial_core_rho_grad(1:3,i_full_points_2)


                          call evaluate_xc  &
                               ( rho_inc_partialcore(1,i_full_points_2),   &
                               rho_gradient_inc_partialcore(1,1,i_full_points_2),  &
                               kinetic_density(1,i_full_points_2), &
                               en_density_xc(i_point),   &
                               en_density_x, en_density_c,  &
                               local_xc_derivs(1,i_point),  &
                               xc_gradient_deriv(1,1,i_point),  &
                               xc_tau_deriv(1,i_point), &
                               .false., &
                               coord_current)


                      else

                          call evaluate_xc  &
                               ( rho(1,i_full_points_2),   &
                               rho_gradient(1,1,i_full_points_2),  &
                               kinetic_density(1,i_full_points_2), &
                               en_density_xc(i_point),   &
                               en_density_x, en_density_c,  &
                               local_xc_derivs(1,i_point),  &
                               xc_gradient_deriv(1,1,i_point),  &
                               xc_tau_deriv(1,i_point), &
                               .false., &
                               coord_current)

                      endif
                      
		              if (Gnonmf_forces_on) then
			            call evaluate_Gnonmf(rho(:,i_full_points_2),   &
				            rho_gradient(:,:,i_full_points_2),  &
				            local_Gnonmf_derivs(i_point),  &
				            Gnonmf_gradient_deriv(:,i_point), &
				            i_full_points_2)
		              end if                      

                      do i_spin = 1, n_spin, 1

                          local_potential_parts(i_spin) =   &
                               hartree_potential(i_full_points_2) +   &
                               local_xc_derivs(i_spin,i_point)
                               
            			  if (Gnonmf_forces_on) then
            			      local_potential_parts(i_spin) = &
            				  local_potential_parts(i_spin) + &
            				  local_Gnonmf_derivs(i_point)
            			  end if                               

                          call evaluate_H_psi_p2  &
                               ( n_compute_c, n_compute_atoms, n_compute_fns, &
                               l_ylm_max, ylm_tab, one_over_dist_tab,  &
                               radial_wave, h_times_wave(1, i_point, i_spin),  &
                               local_potential_parts(i_spin),  &
                               kinetic_wave, 0.5d0, &
                               rad_index, wave_index, l_index, l_count, fn_atom, &
                               n_zero_compute, zero_index_point &
                               )

                       end do
                       if ((flag_rel.eq.REL_atomic_zora)) then

                          call evaluate_radial_functions_deriv_p0 &
                               (spline_array_start, spline_array_end, &
                               n_compute_atoms, n_compute_fns,  &
                               dist_tab(1,i_point), i_r(1), &
                               atom_index, i_basis_fns_inv, &
                               basis_kinetic_ordered,  &
                               kinetic_wave_deriv(1),  &
                               .true., n_max_compute_fns_dens &
                               )

                          kinetic_wave_deriv =  kinetic_wave_deriv *0.5
                          ! write(use_unit,*) kinetic_wave(1) 

                          if(.not. gga_forces_on .and..not.Gnonmf_forces_on)then ! GGA-force part is not calculated

                             ! and finally, assemble the actual gradients
                             call evaluate_wave_gradient_p2  &
                                  ( n_compute_c, n_compute_atoms, n_compute_fns, &
                                  one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                                  l_ylm_max, ylm_tab,  &
                                  dylm_dtheta_tab,  &
                                  scaled_dylm_dphi_tab,  &
                                  kinetic_wave,  &
                                  kinetic_wave_deriv,  &
                                  kinetic_gradient_basis_wave  (1:n_compute_c,1:3,i_point),  &
                                  rad_index, wave_index, l_index, l_count, fn_atom, &
                                  n_zero_compute, zero_index_point &
                                  )


                          else

                             call evaluate_wave_gradient_cartesian_p2 &
                                  ( n_compute_c, n_compute_atoms, n_compute_fns, &
                                  one_over_dist_tab, dir_tab(1,1,i_point), l_wave_max,  &
                                  ylm_tab, kinetic_wave, kinetic_wave_deriv, &
                                  sum_gradient, kinetic_gradient_basis_wave(1,1,i_point), &
                                  rad_index, wave_index, l_index, l_count, fn_atom, &
                                  n_zero_compute, zero_index_point,n_max_compute_dens &
                                  )


                          end if

                       end if !flag_rel == atomic_zora



                    end if ! forces_on
                    if (.not. point_on_atom) then
                        if(use_small_component )then

                        do i_spin = 1, n_spin, 1
                            call small_component(n_compute_c, i_basis,  n_compute_atoms, dist_tab_sq(1,i_point), &
                                atom_index_inv(1), atom_index(1), gradient_basis_wave(1,1,i_point), &
                                temp_rho_small(i_point,i_spin), n_max_compute_dens, i_spin)
                        end do
                        end if



                        ! reset i_basis_fns_inv
                        i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
                    end if

                 end if ! end if (partition_tab.gt.0)
              enddo ! end loop over one batch of the grid
              ! - all quantities that are evaluated pointwise are now known ...

              if (n_points.gt.0) then
                 ! Now perform all operations which are done across the entire integration
                 ! shell at once.

                 if (forces_on) then
                    ! Determine the atoms to whose forces the present integration
                    ! shell will actually contribute
                    call prune_force_atoms_v1(n_compute_c, i_basis,   &
                         global_atom, basis_offset, n_compute_force_atoms)
                 end if


                 do i_spin = 1, n_spin, 1

                    ! VB: Zero max_occ_number can happen if one spin channel is constrained
                    !     to be completely empty. 
                    if (max_occ_number(i_spin).gt.0) then

                       ! New density is always evaluated
                       if(real_eigenvectors)then

                          call evaluate_KS_density_p0  &
                               (  n_points, wave(1,1), n_compute_c,   &
                               i_basis, KS_vec_times_occ_sqrt(1,1,i_spin),   &
                               KS_ev_compute(1,1,i_spin),  &
                               max_occ_number(i_spin),  &
                               KS_orbital(1,1,i_spin),   &
                               temp_rho(1,i_spin), n_max_compute_dens, &
                               n_centers_basis_T )


                          if(use_small_component )then
                             temp_rho = temp_rho + temp_rho_small
                             temp_rho_small = 0.d0
                          end if




                       else

                          call evaluate_KS_density_complex_p0  &
                               (  n_points, wave(1,1), n_compute_c,   &
                               i_basis, KS_vec_times_occ_sqrt_complex(1,1,i_spin),   &
                               KS_ev_compute_complex(1,1,i_spin),  &
                               max_occ_number(i_spin),  &
                               KS_orbital_complex(1,1,i_spin),   &
                               temp_rho(1,i_spin), n_max_compute_dens, &
                               n_centers_basis_T )

                       end if

                       ! The order of evaluation depends heavily on whether or not forces are computed
                       ! We can exploit some synergy here, specifically for KS_orbital_gradient 
                       ! forces are evaluated term by term in order to save a little memory

                       if (forces_on) then


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ! Hellman-Feynman forces
                          if(.not. force_new_functional)then
                             call evaluate_hellman_feynman_forces_p0  &
                                  (hellman_feynman_forces,   &
                                  temp_rho(1,i_spin), partition,   &
                                  dist_tab_sq_global,   &
                                  dir_tab_global, n_points)
                          end if


!                         if(nlcc_forces_on) then
!                             call evaluate_nlcc_forces &
!                                (local_xc_derivs(1,1),&
!                                partition,partial_core_rho_grad_in_batch, dist_tab_sq_global, &
!                                dir_tab_global, n_points)
!                         endif



                          if(flag_rel.eq.REL_atomic_zora )then

                             if (use_density_gradient) then
                                ! must initialize KS_orbital_gradient once ...
                                ! w/o density gradient, we never need the full KS orbital gradient
                                KS_orbital_gradient = 0.d0
                             end if

                             call evaluate_wave_psi &
                                  (n_points,   &
                                  wave(1,1), n_compute_c,  &
                                  KS_ev_compute(1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  KS_eigenvalue(1,i_spin,1), KS_orbital(1,1,i_spin),i_spin,  &
                                  h_minus_e_times_psi(1,1,i_spin))
  
                             do i_atom_2 = 1, n_compute_force_atoms, 1
                                ! index the basis functions which are attached to current atom
                                n_local_compute = basis_offset(i_atom_2 + 1) -   &
                                     basis_offset(i_atom_2)

				if (Gnonmf_forces_on.or.gga_forces_on) then
				  ! Nuclear gradients are always evaluated, but KS_orbital_gradient
				  ! only stored if (use_density_gradient) is actually true.
				  call evaluate_nuclear_gradients_p0(n_points,  &
					  kinetic_gradient_basis_wave, n_compute_c,   &
					  KS_ev_compute(1,1,i_spin),   &
					  max_occ_number(i_spin),  &
					  i_atom_2, basis_offset, n_local_compute,  &
					  nuclear_gradient(1,1,1,i_spin),  &
					  .true., &
					  KS_orbital_gradient(1,1,1,i_spin))
				else
				  call evaluate_nuclear_gradients_p0(n_points,  &
					  kinetic_gradient_basis_wave, n_compute_c,   &
					  KS_ev_compute(1,1,i_spin),   &
					  max_occ_number(i_spin),  &
					  i_atom_2, basis_offset, n_local_compute,  &
					  nuclear_gradient(1,1,1,i_spin),  &
					  .false., &
					  KS_orbital_gradient(1,1,1,i_spin))
                                end if

                                call evaluate_pulay_forces_p0  &
                                  (i_atom_2, h_minus_e_times_psi(1,1,i_spin),  &
                                  nuclear_gradient(1,1,1,i_spin),  &
                                  partition, max_occ_number(i_spin),  &
                                  n_points, n_compute_force_atoms,   &
                                  global_atom,  &
                                  pulay_forces)
                             end do

                          end if ! atomic_zora
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          ! preparation for Pulay and gga forces
                          if (use_density_gradient) then
                             ! must initialize KS_orbital_gradient once ...
                             ! w/o density gradient, we never need the full KS orbital gradient
                             KS_orbital_gradient = 0.d0
                          end if

                          call evaluate_h_minus_e_times_psi_v2 &
                               (n_points,   &
                               h_times_wave(1,1,i_spin), n_compute_c,  &
                               KS_ev_compute(1,1,i_spin),   &
                               max_occ_number(i_spin),   &
                               KS_eigenvalue(1,i_spin,1), KS_orbital(1,1,i_spin),i_spin,  &
                               h_minus_e_times_psi(1,1,i_spin))  

                         if (gga_forces_on.or.Gnonmf_forces_on) then

                             call evaluate_KS_orbital_gradients_p1  &
                                  (n_points, gradient_basis_wave, n_compute_c,  &
                                  KS_ev_compute(1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  KS_orbital_gradient(1,1,1,i_spin), &
                                  n_max_compute_dens)
                                  
                               !evaluate (grad psi_i) * xc_gradient_deriv/Gnonmf_gradient_deriv
                               if (gga_forces_on) then

				call evaluate_orb_grad_dot_rho_grad_p0  &
				    (KS_orbital_gradient(1,1,1,i_spin),  &
				    xc_gradient_deriv(1:3,i_spin,1:n_points),  &
				    n_points, max_occ_number(i_spin),  &
				    orb_grad_dot_rho_grad(1,1,i_spin))
				    
			       end if
			       if (Gnonmf_forces_on) then

				call evaluate_orb_grad_dot_rho_grad_p0  &
				    (KS_orbital_gradient(1,1,1,i_spin),  &
				    Gnonmf_gradient_deriv(1:3,1:n_points),  &
				    n_points, max_occ_number(i_spin),  &
				    orb_grad_dot_rho_grad_Gnonmf(1,1,i_spin))
				    			       
			       end if
                          endif

                          do i_atom_2 = 1, n_compute_force_atoms, 1
                             ! index the basis functions which are attached to current atom
                             n_local_compute = basis_offset(i_atom_2 + 1) -   &
                                  basis_offset(i_atom_2)

                             ! Nuclear gradients are always evaluated, but KS_orbital_gradient
                             ! only stored if (use_density_gradient) is actually true.
                             call evaluate_nuclear_gradients_p0(n_points,  &
                                gradient_basis_wave, n_compute_c,   &
                                KS_ev_compute(1,1,i_spin),   &
                                max_occ_number(i_spin),  &
                                i_atom_2, basis_offset, n_local_compute,  &
                                nuclear_gradient(1,1,1,i_spin), &
                                gga_or_Gnonmf,  &
                                KS_orbital_gradient(1,1,1,i_spin))

                            ! Pulay forces

                            call evaluate_pulay_forces_p0  &
                               (i_atom_2, h_minus_e_times_psi(1,1,i_spin),  &
                               nuclear_gradient(1,1,1,i_spin),  &
                               partition, max_occ_number(i_spin),  &
                               n_points, n_compute_force_atoms,   &
                               global_atom,  &
                               pulay_forces)

                            ! GGA forces
                            if (gga_forces_on) then
                                call evaluate_orb_hess_dot_rho_grad_p0  &
                                     (hessian_basis_wave, n_compute_c,   &
                                     KS_ev_compute(1,1,i_spin),   &
                                     xc_gradient_deriv(1:3,i_spin,1:n_points),  &
                                     n_points, max_occ_number(i_spin),  &
                                     i_atom_2, basis_offset, n_local_compute,  &
                                     orb_hess_dot_rho_grad(1,1,1,i_spin), &
                                     ! Additional terms for mGGA
                                     KS_orbital_gradient(1,1,1,i_spin), &
                                     xc_tau_deriv(i_spin,1:n_points), &
                                     orb_hess_dot_orb_grad(1,1,i_spin), &
                                     meta_gga_forces_on)

                                call evaluate_gga_forces_p0  &
                                     (i_atom_2,KS_orbital(1,1,i_spin),  &
                                     nuclear_gradient(1,1,1,i_spin),  &
                                     orb_grad_dot_rho_grad(1,1,i_spin),  &
                                     orb_hess_dot_rho_grad(1,1,1,i_spin),  &
                                     partition, max_occ_number(i_spin),  &
                                     n_points, global_atom, gga_forces, &
                                     ! Additional term for mGGA
                                     orb_hess_dot_orb_grad(1,1,i_spin), &
                                     meta_gga_forces_on)

                                ! AJL: Original version that didn't include MGGA
                                !
                                ! call evaluate_orb_hess_dot_rho_grad_p0  &
                                !      (hessian_basis_wave, n_compute_c,   &
                                !      KS_ev_compute(1,1,i_spin),   &
                                !      xc_gradient_deriv(1:3,i_spin,1:n_points),  &
                                !      n_points, max_occ_number(i_spin),  &
                                !      i_atom_2, basis_offset, n_local_compute,  &
                                !      orb_hess_dot_rho_grad  &
                                !      (1,1,1,i_spin) ) 
                                !
                                ! call evaluate_gga_forces_p0  &
                                !      (i_atom_2,KS_orbital(1,1,i_spin),  &
                                !      nuclear_gradient(1,1,1,i_spin),  &
                                !      orb_grad_dot_rho_grad(1,1,i_spin),  &
                                !      orb_hess_dot_rho_grad(1,1,1,i_spin),  &
                                !      partition, max_occ_number(i_spin),  &
                                !      n_points, n_compute_force_atoms,   &
                                !      global_atom, gga_forces)

                            end if
                            
                            if (Gnonmf_forces_on) then
				!SR: 1st tabulate (nabla_at nabla psi_i) *Gnonmf_gradient_deriv
				!(orb_hess_dot_rho_grad_Gnonmf)
                                call evaluate_orb_hess_dot_rho_grad_Gnonmf  &
                                     (hessian_basis_wave, n_compute_c,   &
                                     KS_ev_compute(1,1,i_spin),   &
                                     Gnonmf_gradient_deriv(1:3,1:n_points),  &
                                     n_points, max_occ_number(i_spin),  &
                                     i_atom_2, basis_offset, n_local_compute,  &
                                     orb_hess_dot_rho_grad_Gnonmf(1,1,1,i_spin))
				!now calculate the sum over the 2 contributions
				!(nabla_at nabla psi)psi + (nabla_at psi)(nabla psi)
                                call evaluate_Gnonmf_forces  &
                                     (i_atom_2,KS_orbital(1,1,i_spin),  &
                                     nuclear_gradient(1,1,1,i_spin),  &
                                     orb_grad_dot_rho_grad_Gnonmf(1,1,i_spin),  &
                                     orb_hess_dot_rho_grad_Gnonmf(1,1,1,i_spin),  &
                                     partition, max_occ_number(i_spin),  &
                                     n_points, global_atom, Gnonmf_forces)
                            end if
                            

                        end do



                       else if (use_density_gradient) then
                          ! If we do not compute forces but need the density
                          ! gradient anyway, compute KS_orbital_gradient directly

                          if(real_eigenvectors)then

                             call evaluate_KS_orbital_gradients_p1  &
                                  (n_points, gradient_basis_wave, n_compute_c,  &
                                  KS_ev_compute(1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  KS_orbital_gradient(1,1,1,i_spin), &
                                  n_max_compute_dens)

                          else

                             call evaluate_KS_orbital_gradients_complex_p1  &
                                  (n_points, gradient_basis_wave, n_compute_c,  &
                                  KS_ev_compute_complex(1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  KS_orbital_gradient_complex(1,1,1,i_spin), &
                                  n_max_compute_dens)

                          end if

                       end if ! if force on

                       ! Get the density gradient ...
                       if (use_density_gradient) then

                          if(real_eigenvectors)then


                             call evaluate_density_gradient_p1  &
                                  (n_points, n_compute_c, KS_orbital(1,1,i_spin),    &
                                  KS_orbital_gradient(1,1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  temp_rho_gradient(1,1,i_spin) &
                                  )


                          else

                             call evaluate_density_gradient_complex_p1  &
                                  (n_points, n_compute_c, KS_orbital_complex(1,1,i_spin),    &
                                  KS_orbital_gradient_complex(1,1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  temp_rho_gradient(1,1,i_spin) &
                                  )

                          end if

                          if (use_meta_gga) then
                            if(real_eigenvectors)then

                             call evaluate_density_kinetic_p1  &
                                  (n_points, n_compute_c,    &
                                  KS_orbital_gradient(1,1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  temp_kinetic_density(1,i_spin) )

                            else

                             call evaluate_density_kinetic_complex_p1  &
                                  (n_points, n_compute_c,     &
                                  KS_orbital_gradient_complex(1,1,1,i_spin),   &
                                  max_occ_number(i_spin),   &
                                  temp_kinetic_density(1,i_spin) )

                            end if
                          end if ! use_meta_gga
                       end if ! use_density_gradient

                    else
                       ! case of max_occ_number = 0
                       ! use default value zero for everything

                       temp_rho(1:n_points,i_spin) = 0.d0
                       temp_rho_gradient(1:3,1:n_points,i_spin) = 0.d0
                       temp_kinetic_density(1:n_points,i_spin) = 0.d0

                       ! Note: Forces are addititive and initialized to zero,
                       !       hence no extra treatment of force terms needed

                    end if
                 end do
              end if ! (n_points>0)

              ! calculate change in electron density
              i_point = 0
              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_3 = i_full_points_3 + 1

                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).gt.0 .or. (atomic_MERM .and. mpb_solver_started)) then

                    i_point = i_point + 1

                    ! and local change in charge density and gradient for Pulay mixing
                    if (spin_treatment.eq.0) then

                       delta_rho_KS(i_full_points_3, 1) =   &
                            delta_rho_KS(i_full_points_3, 1) & 
                            + temp_rho(i_point,1) - rho(1, i_full_points_3)/dble(n_k_group)

                    elseif (spin_treatment.eq.1) then

                       delta_rho_KS(i_full_points_3, 1) =   &
                            delta_rho_KS(i_full_points_3, 1)   &
                            +( temp_rho(i_point, 1) +   &
                            temp_rho(i_point, 2) )-  &
                            ( rho(1,i_full_points_3) + rho(2,i_full_points_3) )/dble(n_k_group)

                       delta_rho_KS(i_full_points_3, 2) =   &
                            delta_rho_KS(i_full_points_3, 2)   &
                            +( temp_rho(i_point, 1) -   &
                            temp_rho(i_point, 2) )-  &
                            ( rho(1,i_full_points_3) - rho(2, i_full_points_3) )/dble(n_k_group)

                    end if

                    ! prepare charge density root-mean-square distance

                    if(i_k_point_group == n_k_group)then
                       do i_spin = 1, n_spin, 1

                          rho_change(i_spin) = rho_change(i_spin) +  &
                               partition_tab(i_full_points_3) *  &
                               delta_rho_KS(i_full_points_3, i_spin)**2

                       enddo
                    end if

                    if (use_density_gradient) then
                       !prepare delta_rho_gradient for later mixing - this
                       !will be mixed in exactly the same way as delta_rho_KS
                       if (spin_treatment.eq.0) then

                          do i_coord = 1,3,1
                             delta_rho_gradient(i_coord, i_full_points_3, 1) =   &
                                  delta_rho_gradient(i_coord, i_full_points_3, 1) &   
                                  + temp_rho_gradient(i_coord, i_point,1)-  &
                                  rho_gradient(i_coord,1,i_full_points_3)/n_k_group
                          enddo

                          if (use_meta_gga) then

                             delta_kinetic_density(i_full_points_3, 1) =   &
                                  delta_kinetic_density(i_full_points_3, 1) + &
                                  temp_kinetic_density(i_point,1) - &
                                  kinetic_density(1,i_full_points_3)/n_k_group                                  
                          endif

                       elseif (spin_treatment.eq.1) then

                          do i_coord = 1,3,1

                             delta_rho_gradient(i_coord, i_full_points_3, 1) = &   
                                  delta_rho_gradient(i_coord, i_full_points_3, 1) &
                                  + (temp_rho_gradient(i_coord, i_point,1)+  &
                                  temp_rho_gradient(i_coord, i_point,2))- &
                                  (rho_gradient(i_coord,1,i_full_points_3)+ &
                                  rho_gradient(i_coord,2,i_full_points_3))/n_k_group

                             delta_rho_gradient(i_coord, i_full_points_3, 2) = &
                                  delta_rho_gradient(i_coord, i_full_points_3, 2) &
                                  + (temp_rho_gradient(i_coord, i_point,1)- &
                                  temp_rho_gradient(i_coord, i_point,2))- &
                                  (rho_gradient(i_coord,1,i_full_points_3)- &
                                  rho_gradient(i_coord,2,i_full_points_3))/n_k_group

                          enddo

                          if (use_meta_gga) then

                             delta_kinetic_density(i_full_points_3, 1) = &
                                  delta_kinetic_density(i_full_points_3, 1) + &
                                  (temp_kinetic_density(i_point,1) + &
                                   temp_kinetic_density(i_point,2)) - &
                                  (kinetic_density(1,i_full_points_3) + &
                                  kinetic_density(2,i_full_points_3))/n_k_group

                             delta_kinetic_density(i_full_points_3, 2) = &
                                  delta_kinetic_density(i_full_points_3, 2) + &
                                  (temp_kinetic_density(i_point,1) - &
                                   temp_kinetic_density(i_point,2)) - &
                                  (kinetic_density(1,i_full_points_3)- &
                                  kinetic_density(2,i_full_points_3))/n_k_group

                          end if

                       end if
                    end if
                 endif
              end do

           else
              ! Even if n_compute .eq. 0 for the entire current batch of grid points, we still need to
              ! make sure that the density _change_ at this point is ( zero minus previous density ).
              ! This ensures that. even for a zero KS density, a potentially non-zero initialization density
              ! is subtracted to properly account for the density change ....

              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1
                 i_full_points_3 = i_full_points_3 + 1

                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then

                    ! local change in charge density and gradient for Pulay mixing

                    if (spin_treatment.eq.0) then

                       delta_rho_KS(i_full_points_3, 1) =   &
                            delta_rho_KS(i_full_points_3, 1) & 
                            - rho(1, i_full_points_3)/dble(n_k_group)

                    elseif (spin_treatment.eq.1) then

                       delta_rho_KS(i_full_points_3, 1) =   &
                            delta_rho_KS(i_full_points_3, 1) -  &
                            ( rho(1,i_full_points_3) + rho(2,i_full_points_3) )/dble(n_k_group)

                       delta_rho_KS(i_full_points_3, 2) =   &
                            delta_rho_KS(i_full_points_3, 2) -  &
                            ( rho(1,i_full_points_3) - rho(2, i_full_points_3) )/dble(n_k_group)

                    end if

                    ! prepare charge density root-mean-square distance
                    if(i_k_point_group == n_k_group)then
                       do i_spin = 1, n_spin, 1

                          rho_change(i_spin) = rho_change(i_spin) +  &
                               partition_tab(i_full_points_3) *  &
                               delta_rho_KS(i_full_points_3, i_spin)**2

                       enddo
                    end if

                    if (use_density_gradient) then
                       !prepare delta_rho_gradient for later mixing - this
                       !will be mixed in exactly the same way as delta_rho_KS
                       if (spin_treatment.eq.0) then

                          do i_coord = 1,3,1
                             delta_rho_gradient(i_coord, i_full_points_3, 1) =   &
                                  delta_rho_gradient(i_coord, i_full_points_3, 1) - &
                                  rho_gradient(i_coord,1,i_full_points_3)/n_k_group
                          enddo

                          if (use_meta_gga) then

                             delta_kinetic_density(i_full_points_3, 1) =   &
                             delta_kinetic_density(i_full_points_3, 1) - &
                             kinetic_density(1,i_full_points_3)/n_k_group

                          end if ! use_meta_gga

                       elseif (spin_treatment.eq.1) then

                          do i_coord = 1,3,1

                             delta_rho_gradient(i_coord, i_full_points_3, 1) = &   
                                  delta_rho_gradient(i_coord, i_full_points_3, 1) - &
                                  (rho_gradient(i_coord,1,i_full_points_3)+ &
                                  rho_gradient(i_coord,2,i_full_points_3))/n_k_group

                             delta_rho_gradient(i_coord, i_full_points_3, 2) = &
                                  delta_rho_gradient(i_coord, i_full_points_3, 2) - &
                                  (rho_gradient(i_coord,1,i_full_points_3)- &
                                  rho_gradient(i_coord,2,i_full_points_3))/n_k_group

                          enddo

                          if (use_meta_gga) then

                             delta_kinetic_density(i_full_points_3, 1) = &
                             delta_kinetic_density(i_full_points_3, 1) - &
                             (kinetic_density(1,i_full_points_3) + &
                             kinetic_density(2,i_full_points_3))/n_k_group

                             delta_kinetic_density(i_full_points_3, 2) = &
                             delta_kinetic_density(i_full_points_3, 2) - &
                             (kinetic_density(1,i_full_points_3) - &
                             kinetic_density(2,i_full_points_3))/n_k_group
                                  
                          end if ! use_meta_gga
                       end if
                    end if
                 endif

              enddo
           end if  ! end if (n_compute.gt.0)
        ! end if ! end distribution over threads
     end do ! end loop over batches
  end do ! end loop over groups of k-points

  if (forces_on) then
     ! Pulay forces still need an extra factor of 2
     pulay_forces = 2.d0*pulay_forces
     if (gga_forces_on) then
        ! gga forces still need an extra factor of 4
        gga_forces = 4.d0*gga_forces
     end if
     if (Gnonmf_forces_on) then
	    !need a factor of 2 here 
	    Gnonmf_forces = 2.d0*Gnonmf_forces
     end if
  end if

!  call sync_pp_charge_forces(pp_nlcc_forces)

  ! broadcast the result to all threads
  call sync_density(rho_change)

  do i_spin = 1, n_spin, 1
     rho_change(i_spin) = sqrt(rho_change(i_spin))
  enddo

  ! Allocatable arrays that are tracked
  if (allocated(wave))                call aims_deallocate( wave,                               "wave" )
  if (allocated(gradient_basis_wave)) call aims_deallocate( gradient_basis_wave, "gradient_basis_wave" )
  if (allocated(h_times_wave))        call aims_deallocate( h_times_wave,               "h_times_wave" )
  if (allocated(hessian_basis_wave))  call aims_deallocate( hessian_basis_wave,   "hessian_basis_wave" )

  ! finally, deallocate stuff.

  if (allocated( temp_rho_small       )) deallocate( temp_rho_small       )

  if (allocated(dylm_dtheta_tab)) then
     deallocate (dylm_dtheta_tab)
  end if
  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate (scaled_dylm_dphi_tab)
  end if
  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if (allocated(KS_orbital_gradient)) then
     deallocate(KS_orbital_gradient)
  end if
  if (allocated(nuclear_gradient)) then
     deallocate(nuclear_gradient)
  end if
  if (allocated(dir_tab_global)) then
     deallocate(dir_tab_global)
  end if
  if (allocated(dist_tab_sq_global)) then
     deallocate(dist_tab_sq_global)
  end if
  if (allocated(cartesians)) then
     deallocate(cartesians)
  end if
  if (allocated(sum_gradient)) then
     deallocate(sum_gradient)
  end if
  if (allocated(sum_hessian)) then
     deallocate(sum_hessian)
  end if
  if (allocated(h_minus_e_times_psi)) then
     deallocate(h_minus_e_times_psi)
  end if
  if (allocated(orb_grad_dot_rho_grad)) then
     deallocate(orb_grad_dot_rho_grad)
  end if
  if (allocated(orb_hess_dot_rho_grad)) then
     deallocate(orb_hess_dot_rho_grad)
  end if
  if (allocated(orb_hess_dot_orb_grad)) then
     deallocate(orb_hess_dot_orb_grad)
  end if
  if(allocated(radial_wave))then
     deallocate(radial_wave)
  end if
  if(allocated(radial_wave_deriv))then
     deallocate(radial_wave_deriv)
  end if
  if(allocated(radial_wave_2nd_deriv))then
     deallocate(radial_wave_2nd_deriv)
  end if
  if(allocated(kinetic_wave))then
     deallocate(kinetic_wave)
  end if
  if(allocated(i_basis))then
     deallocate(i_basis)
  end if
  if( allocated(KS_vec_times_occ_sqrt))then
     deallocate(KS_vec_times_occ_sqrt)
  end if
  if(allocated(KS_vec_times_occ_sqrt_complex))then
     deallocate(KS_vec_times_occ_sqrt_complex)
  end if
  if(allocated( KS_ev_compute))then
     deallocate( KS_ev_compute)
  end if
  if(allocated(KS_ev_compute_complex))then
     deallocate(KS_ev_compute_complex)
  end if
  if(allocated(KS_ev_compute_complex))then
     deallocate(KS_ev_compute_complex)
  end if
  if(allocated( KS_orbital))then
     deallocate( KS_orbital)
  end if
  if(allocated( KS_orbital_complex))then
     deallocate( KS_orbital_complex)
  end if

  if(allocated( rho_inc_partialcore  )) then 
     deallocate( rho_inc_partialcore  )
  endif

  if(allocated( rho_gradient_inc_partialcore )) then 
     deallocate( rho_gradient_inc_partialcore )
  endif

  if (Gnonmf_forces_on) then
    if (allocated(orb_hess_dot_rho_grad_Gnonmf)) then
      deallocate(orb_hess_dot_rho_grad_Gnonmf)
    end if
    if (allocated(orb_grad_dot_rho_grad_Gnonmf)) then
      deallocate(orb_grad_dot_rho_grad_Gnonmf)
    end if
    if (allocated(orb_hess_dot_orb_grad_Gnonmf)) then
      deallocate(orb_hess_dot_orb_grad_Gnonmf)
    end if
    if (allocated(local_Gnonmf_derivs)) then
      deallocate(local_Gnonmf_derivs)
    end if
    if (allocated(Gnonmf_gradient_deriv)) then
      deallocate(Gnonmf_gradient_deriv)
    end if
  end if
  if(allocated(kinetic_gradient_basis_wave)) deallocate(kinetic_gradient_basis_wave)
  if(allocated(kinetic_wave_deriv)) deallocate(kinetic_wave_deriv)


end subroutine update_density_and_forces_orbital
!******
