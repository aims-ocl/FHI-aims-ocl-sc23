!****s* FHI-aims/initialize_integrals_p0
!  NAME
!   initialize_integrals_p0
!  SYNOPSIS

subroutine initialize_integrals_p0 &
     ( basis_l_max, overlap_matrix, hamiltonian &
     )

!  PURPOSE
!  The subroutine integrates Hamiltonian and overlapmatrix at first time.
!  If asked the angular integration grid is adaptively optimized using
!  Hamiltonian and overlap matrix as a test sample.
!
!  Structure:
!  
!   Loop over atoms in structure
!   |
!   | Loop over all radial integration shells
!   | |
!   | |
!   | | if maximum angular grid required then
!   | | |
!   | | | * simply tabulate everything on fixed grids
!   | | |
!   | | else
!   | | |
!   | | | Loop over increasing Lebedev grids
!   | | | |
!   | | | | * tabulate partition_tab
!   | | | | * add up all overlap matrix contributions from present shell
!   | | | |
!   | | | until all overlap matrix contributions from present radial shell are converged 
!   | | | with respect to Lebedev grid density
!   | | |
!   | | | Tabulate all Hamiltonian matrix contributions for present Lebedev grid
!   | | |
!   | | | Loop over increasing Lebedev grids
!   | | | |
!   | | | | * tabulate partition_tab
!   | | | | * tabulate potential components / potential
!   | | | | * add up all Hamiltonian matrix contributions from present shell
!   | | | |
!   | | | until all Hamiltonian matrix contributions from present radial shell are converged 
!   | | | with respect to Lebedev grid density
!   | | | 
!   | | | Store integration grid data for present shell. 
!   | | | 
!   | | end if
!   | | 
!   | next radial shell
!   |
!   next atom in structure
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use localorb_io
  use synchronize_mpi
  use constants
  use pbc_lists
  use species_data
  use read_fixed_grid
  use analyze_arrays
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none


!  ARGUMENTS

  integer basis_l_max (n_species)
  real*8 overlap_matrix( n_hamiltonian_matrix_size)
  real*8 hamiltonian(n_hamiltonian_matrix_size, n_spin )

!  INPUTS
!   o basis_l_max -- maximum l index of basis functions
!  OUTPUT
!   o overlap_matrix -- overlap matrix
!   o hamiltonian -- Hamiltonian matrix
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




  ! local variables

  ! variables for local geometry description
  ! trigonom_tab contains, for each atom, the trigonometric functions
  !       cos(theta), sin(theta), cos(phi), sin(phi) of the 
  !       current integration point

  real*8 coord_current(3, n_max_angular)


  real*8,dimension(:,:,:),allocatable:: dist_tab

!  real*8 dist_tab( n_centers_basis_integrals, n_max_angular, n_max_angular_division)
  real*8 dist_tab_full( n_centers_integrals, n_max_angular)
  real*8 dist_tab_par_tab( n_centers_basis_integrals)
  real*8 dist_tab_sq( n_centers_integrals, n_max_angular)
  real*8 dist_tab_sq_full( n_centers_integrals, n_max_angular)
  real*8 dist_tab_sq_par_tab( n_centers_basis_integrals)
  real*8, dimension(:,:,:,:),allocatable::  dir_tab

!  real*8 dir_tab(3, n_centers_basis_integrals,n_max_angular, n_max_angular_division)

  real*8 dir_tab_full(3, n_centers_integrals,n_max_angular)

  real*8 dir_tab_full_norm(3, n_centers_integrals,n_max_angular)

  real*8 dir_tab_par_tab(3, n_centers_basis_integrals)
  real*8 dir_tab_par_tab_norm(3, n_centers_basis_integrals)

!  real*8 i_r( n_centers_basis_integrals, n_max_angular, n_max_angular_division)
  real*8,dimension( :,:,:),allocatable::  i_r
  real*8 i_r_full( n_centers_integrals, n_max_angular)
  real*8 i_r_par_tab( n_centers_basis_integrals)

!  real*8 trigonom_tab(4, n_centers_integrals,n_max_angular, &
!       n_max_angular_division)
  real*8,dimension(:,:,:,:),allocatable:: trigonom_tab

  ! variables for grid quantities over one angular shell

  real*8, dimension(n_max_angular) :: partition_tab

  ! variables which need no longer be stored over one angular shell
  real*8 :: free_hartree_superpos 
  !  Nadia
  real*8 :: pot_ion_embed_local
  real*8 :: free_rho_superpos 
  real*8, dimension(3) :: free_rho_gradient_superpos 
  real*8, dimension(n_spin) :: rho 
  real*8, dimension(n_spin,n_max_angular) :: local_potential_parts 
  real*8, dimension(n_spin,n_max_angular) :: zora_potential_parts 
  ! auxiliary variables for density gradient at a given point, to
  ! avoid referencing rho_gradient directly if it is not allocated
  real*8 local_rho_gradient(3,n_spin)
  real*8 local_kinetic_density(n_spin)
  real*8 local_xc_gradient_deriv(3,n_spin)
  real*8 local_xc_tau_deriv(n_spin)
  real*8 local_gradient(3,n_spin,n_max_angular)

  ! need ylm_tab allocatable to switch between calculations with / without 
  ! wave function gradient. The allocation dimension is l_ylm_max, globally

  integer :: l_ylm_max
  real*8, dimension(:,:,:,:), allocatable :: ylm_tab
  integer, dimension(:,:), allocatable :: index_lm

  ! dylm_dtheta_tab contains dy_(lm)/d(theta) 
  ! scaled_dylm_dphi_tab contains 1/sin(theta)*dy_(lm)/d(phi), the
  ! scaled derivative which is needed to compute cartesian
  ! gradients and avoids the apparent coordinate singularity
  !  on the z axis ( sin(theta)=0 ). 
  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

  ! variables for potential

  real*8 :: local_xc_derivs(n_spin)
  real*8 :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c

  ! variables for non-relativistic part of Hamiltonian

  real*8, dimension(:,:,:), allocatable :: radial_wave
  real*8, dimension(:,:,:), allocatable :: wave
  real*8, dimension(:), allocatable :: kinetic_wave
  real*8, dimension(:), allocatable :: radial_wave_deriv
  real*8, dimension(:,:,:), allocatable :: H_times_psi

  integer :: n_compute_atoms

  integer :: atom_index( n_centers_integrals)
  integer :: atom_index_inv( n_centers)

  ! integer :: n_compute(n_max_angular_division)
  integer :: n_compute_a(n_max_angular_division)
  integer :: n_compute_c(n_max_angular_division)
  integer :: i_basis(n_centers_basis_I,n_max_angular_division)

  integer :: n_compute_fns

  integer :: i_basis_fns    (n_basis_fns*n_centers_integrals)
  ! integer :: i_basis_fns_inv(n_basis_fns,n_centers_basis_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns     (n_basis_fns*n_centers_integrals)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

  ! optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points(n_max_angular_division)
  integer :: try_n_points(n_max_angular_division)
  integer :: n_rel_points
  integer :: try_n_rel_points

  ! and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_angular)

  ! gradient_basis_wave stores the full (3d) gradient of each basis function
  real*8, dimension(:,:), allocatable :: gradient_basis_wave

  ! variables for ZORA operator
  ! zora_operator is the first order kinetic energy correction between 
  !                   the wave function derivatives
  ! t_zora determines whether zora is needed

  real*8, dimension(n_spin) :: zora_operator
  logical :: t_zora
  !      real*8, dimension(:,:,:,:), allocatable :: zora_vector

  ! grid- or grid-optimisation related variables

  logical shell_converged
  logical first_loop
  logical:: full_shell(n_max_radial, n_atoms)


  ! n_lebedev is the index (number) of the Lebedev grid that we have currently
  ! reached for a given atom and a given radial shell

  integer, dimension(:,:), allocatable :: n_lebedev

  ! variables for trial angular grids
  integer :: try_n_angular
  real*8, dimension(:,:), allocatable :: try_r_angular
  real*8, dimension(:),   allocatable :: try_w_angular
  real*8, dimension(:),   allocatable :: try_partition_tab
  real*8, dimension(:),   allocatable :: try_partition
  integer :: try_n_division
  integer, dimension(:), allocatable :: try_division_boundaries

  ! derived variables for trial grids - needed because we converge the Hamiltonian matrix
  ! after the overlap matrix is already done
  real*8, dimension(:,:),     allocatable :: try_coord
  real*8, dimension(:,:,:),   allocatable :: try_dist_tab
  real*8, dimension(:,:),     allocatable :: try_dist_tab_full
  real*8, dimension(:,:),     allocatable :: try_dist_tab_sq
  real*8, dimension(:,:),     allocatable :: try_dist_tab_sq_full
  real*8, dimension(:,:,:),   allocatable :: try_i_r
  real*8, dimension(:,:),     allocatable :: try_i_r_full
  real*8, dimension(:,:,:,:), allocatable :: try_dir_tab
  real*8, dimension(:,:,:),   allocatable :: try_dir_tab_full
  real*8, dimension(:,:,:),   allocatable :: try_dir_tab_full_norm
  real*8, dimension(:,:,:,:), allocatable :: try_trigonom_tab
  real*8, dimension(:,:,:,:), allocatable :: try_ylm_tab
  real*8, dimension(:,:,:),   allocatable :: try_radial_wave
  real*8, dimension(:,:,:),   allocatable :: try_wave

  real*8, dimension(:,:,:),  allocatable :: try_local_potential_parts 
  real*8, dimension(:,:,:),  allocatable :: try_zora_potential_parts 
  real*8, dimension(:,:,:,:),allocatable :: try_local_gradient

  integer :: try_n_compute_c(n_max_angular_division)
  integer :: try_n_compute_a(n_max_angular_division)
  integer, dimension(:,:), allocatable :: try_i_basis

  ! matrix_shell is used for individual integration steps ...
  real*8, dimension(:,:), allocatable :: matrix_shell

  !     variables for trial matrices
  real*8, dimension(:,:), allocatable :: prev_shell
  real*8, dimension(:,:), allocatable :: new_shell

  real*8 :: dummy_force(3)

  real*8 :: pp_potential

  ! counters

  integer i_species
  integer i_atom
  integer i_radial
  integer i_lebedev
  integer i_angular
  integer i_index, i_l, i_m
  integer i_coord
  integer i_division, info

  integer :: i_grid_shell

  integer :: i_spin
  integer :: i_point
  character*100 :: info_str



  !  begin work

  !  write(use_unit,*) myid

  i_spin = 1
  t_zora = .false.
  local_xc_gradient_deriv = 0


  call localorb_info( "",use_unit,'(A)')
  write(info_str,'(2X,A,A)') "Initial 3D integrations: Overlap and Hamiltonian matrix."
  call localorb_info(info_str,use_unit,'(A)')

  write(info_str,'(2X,A,A)') "| Adapting angular integration grids if requested."
  call localorb_info(info_str,use_unit,'(A)')

  ! begin with general allocations
  if (((flag_rel.eq.0).or.(flag_rel.eq.REL_atomic_zora).or.(flag_rel.eq.REL_own)).and.(.not.use_gga)) then
     !       no gradients needed
     l_ylm_max = l_wave_max
  else if ((flag_rel.eq.1).or.use_gga) then
     !       wave function gradients needed
     l_ylm_max = l_wave_max 

     call aims_allocate( gradient_basis_wave, n_centers_basis_I, 3,             "gradient_basis_wave" )

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info )
     call check_allocation(info, 'dylm_dtheta_tab               ') 

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ') 
  end if

  !  allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_integrals, n_max_angular, &
  !       n_max_angular_division) )
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info )
  call check_allocation(info, 'index_lm                      ') 


  ! More allocations of routinely needed variables
  if (.not.allocated(radial_wave)) then
     allocate ( radial_wave(n_max_compute_fns_ham, n_max_points_per_div, n_max_angular_division),stat=info )
     call check_allocation(info, 'radial_wave                   ') 
  end if
  if (.not.allocated(wave)) then
     call aims_allocate( wave, n_max_compute_ham, n_max_points_per_div,n_max_angular_division, "wave" )
  end if
  if (.not.allocated(kinetic_wave)) then
     allocate (kinetic_wave(n_max_compute_fns_ham),stat=info)
     call check_allocation(info, 'kinetic_wave                  ') 
  end if
  if (.not.allocated(radial_wave_deriv)) then
     allocate (radial_wave_deriv(n_max_compute_fns_ham),stat=info)
     call check_allocation(info, 'radial_wave_deriv             ') 
  end if
  if (.not.allocated(H_times_psi)) then
     call aims_allocate( H_times_psi, n_max_compute_ham, n_max_points_per_div, n_spin,  "H_times_psi" )
  end if

  ! always need matrix_shell
  if (.not.allocated(matrix_shell)) then
     call aims_allocate( matrix_shell, n_max_compute_ham, n_max_compute_ham,           "matrix_shell" )
  end if



  ! if needed, allocate and initialize variables for adjustable integration grids
  do i_species = 1, n_species, 1

     if (angular_acc(i_species).gt.0.d0) then
        ! Allocations are needed ...

        if(.not. allocated(i_r))then
           allocate(i_r( n_centers_integrals, n_max_points_per_div, n_max_angular_division),stat=info)
           call check_allocation(info, 'i_r                           ') 
        end if

        if(.not. allocated(trigonom_tab))then
           allocate(trigonom_tab(4, n_centers_integrals,n_max_points_per_div, n_max_angular_division),stat=info)
           call check_allocation(info, 'trigonom_tab                  ') 
        end if

        if( .not. allocated( ylm_tab))then
           allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info )
           call check_allocation(info, 'ylm_tab                       ') 
        end if

        if(.not. allocated(dist_tab))then
           allocate(dist_tab( n_centers_integrals, n_max_points_per_div, n_max_angular_division),stat=info)
           call check_allocation(info, 'mdist_tab                     ') 
        end if

        if(.not. allocated(dir_tab))then
           allocate( dir_tab(3, n_centers_integrals,n_max_points_per_div, n_max_angular_division),stat=info)
           call check_allocation(info, 'dir_tab                       ') 
        end if




        if(.not. allocated(try_local_potential_parts))then
           allocate(try_local_potential_parts(n_spin,n_max_angular,2),stat=info)
           call check_allocation(info, 'try_local_potential_parts     ') 
        end if

        if(.not. allocated( try_zora_potential_parts))then
           allocate( try_zora_potential_parts(n_spin,n_max_angular,2),stat=info)
           call check_allocation(info, 'try_zora_potential_parts      ') 
        end if
     
        if(.not. allocated( try_local_gradient))then
           allocate( try_local_gradient(3,n_spin,n_max_angular,2),stat=info)
           call check_allocation(info, 'try_local_gradient            ') 
        end if

        if (.not.allocated(try_r_angular)) then
           allocate ( try_r_angular(3,n_max_angular),stat=info)
           call check_allocation(info, 'try_r_angular                 ') 
        end if
        if (.not.allocated(try_w_angular)) then
           allocate ( try_w_angular(n_max_angular),stat=info )
           call check_allocation(info, 'try_w_angular                 ') 
        end if
        if (.not.allocated(try_division_boundaries)) then
           allocate ( try_division_boundaries( n_max_angular_division+1),stat=info )
           call check_allocation(info, 'try_division_boundaries       ') 
        end if

        if (.not.allocated(try_partition_tab)) then
           allocate ( try_partition_tab(n_max_angular),stat=info)
           call check_allocation(info, 'try_partition_tab             ') 
        end if
        if (.not.allocated(try_partition)) then
           allocate ( try_partition(n_max_angular),stat=info )
           call check_allocation(info, 'try_partition                 ') 
        end if

        if (.not.allocated(try_coord)) then
           allocate(try_coord(3, n_max_angular),stat=info)
           call check_allocation(info, 'try_coord                     ') 
        end if

        if (.not.allocated(try_dist_tab)) then
           allocate(try_dist_tab(n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_dist_tab                  ') 
         end if
        if (.not.allocated(try_dist_tab_full)) then
           allocate(try_dist_tab_full(n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_dist_tab_full             ') 
        end if
        if (.not.allocated(try_dist_tab_sq)) then
           allocate(try_dist_tab_sq(n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_dist_tab_sq               ') 
        end if
        if (.not.allocated(try_dist_tab_sq_full)) then
           allocate(try_dist_tab_sq_full(n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_dist_tab_sq_full          ') 
        end if
        if (.not.allocated(try_i_r)) then
           allocate(try_i_r(n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_i_r                       ') 
        end if
        if (.not.allocated(try_i_r_full)) then
           allocate(try_i_r_full(n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_i_r_full                  ') 
        end if
        if (.not.allocated(try_dir_tab_full)) then
           allocate(try_dir_tab_full(3, n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_dir_tab_full              ') 
        end if
        if (.not.allocated(try_dir_tab_full_norm)) then
           allocate(try_dir_tab_full_norm(3, n_centers_integrals, n_max_angular),stat=info)
           call check_allocation(info, 'try_dir_tab_full_norm         ') 
        end if
        if (.not.allocated(try_dir_tab)) then
           allocate(try_dir_tab(3, n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_dir_tab                   ') 
        end if
        if (.not.allocated(try_trigonom_tab)) then
           allocate(try_trigonom_tab(4, n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_trigonom_tab              ') 
        end if
        if (.not.allocated(try_ylm_tab)) then
           allocate( try_ylm_tab( (l_ylm_max+1)**2, n_centers_integrals, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_ylm_tab                   ') 
        end if
        if (.not.allocated(try_i_basis)) then
           allocate(try_i_basis(n_centers_basis_I,n_max_angular_division),stat=info)
           call check_allocation(info, 'try_i_basis                   ') 
           try_i_basis = 0
        end if
        if (.not.allocated(try_radial_wave)) then
           allocate(try_radial_wave(n_max_compute_fns_ham, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'mtry_radial_wave              ') 
        end if
        if (.not.allocated(try_wave)) then
           allocate(try_wave(n_max_compute_ham, n_max_points_per_div, &
                n_max_angular_division),stat=info)
           call check_allocation(info, 'try_wave                      ') 
        end if

        if (.not.allocated(prev_shell)) then
           call aims_allocate ( prev_shell, n_hamiltonian_matrix_size, n_spin, "prev_shell" ) 
        end if

        if (.not.allocated(new_shell)) then
           call aims_allocate ( new_shell, n_hamiltonian_matrix_size, n_spin,   "new_shell" ) 
        end if
     end if
  enddo


  ! if angular_acc == 0 for all species, THEN the following
  ! allocations will be made - else all these arrays are
  ! already allocated so that the entire section is silently skipped.
  if(.not. allocated(i_r))then
     allocate(i_r( n_centers_integrals, 1, 1),stat=info)
     call check_allocation(info, 'i_r                           ') 
  end if
  if(.not. allocated(trigonom_tab))then
     allocate(trigonom_tab(4, n_centers_integrals,1,1),stat=info)
     call check_allocation(info, 'trigonom_tab                  ') 
  end if
  if( .not. allocated( ylm_tab))then
     allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_integrals, 1, 1),stat=info)
     call check_allocation(info, 'ylm_tab                       ') 
  end if
  if(.not. allocated(dist_tab))then
     allocate(dist_tab( n_centers_integrals, 1, 1),stat=info)
     call check_allocation(info, 'dist_tab                      ') 
  end if
  if(.not. allocated(dir_tab))then
     allocate(dir_tab(3, n_centers_integrals,1,1),stat=info)
     call check_allocation(info, 'dir_tab                       ') 
  end if


  ! this one we need anyway
  if (.not.allocated(n_lebedev)) then
     allocate( n_lebedev(n_max_radial, n_species),stat=info )
     call check_allocation(info, 'n_lebedev                     ') 
  end if


  ! initialize matrices

  overlap_matrix = 0.
  hamiltonian    = 0.
  full_shell = .false.
  local_gradient        = 0.d0
  local_potential_parts = 0.d0
  zora_potential_parts  = 0.d0

  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  ! initialize basis indecies

  i_basis = 0



  ! if needed, preload adjustable integration grids for each species
  do i_species = 1, n_species, 1
     ! initialize n_lebedev and angular grids

     ! VB: THIS belongs into a very different place in the code - that is, into prepare_scf,
     !     pretty much where get_external_grids is located.
     !
     !     It MUST not stay here
     !

!!$     if(species_fixed_grid(i_species))then
!!$        write(use_unit,*)
!!$        write(use_unit,*) "VB: PLEASE SEE MY REQUEST TO MOVE THIS PART INTO ANOTHER PART OF THE CODE."
!!$        write(use_unit,*) "    READING INPUT FILES SHOULD NOT HAPPEN AT THIS STAGE - IT BELONGS"
!!$        write(use_unit,*) "    INTO PREPARE_SCF. REALLY."
!!$        write(use_unit,*) "    Of course, you can remove the stop for now, but I wanted to "
!!$        write(use_unit,*) "    make sure that my statement was read. Sorry - Volker."
!!$        write(use_unit,*)
!!$        stop
!!$        call fixed_grids_open_for_angular_data(species_name(i_species), 88)
!!$     end if

     do i_radial = 1, n_radial(i_species), 1

        if (specified_grid(i_species)) then

          i_grid_shell = 1
          do while ( (r_ang_shell(i_grid_shell,i_species).lt.r_radial(i_radial, i_species)) .and. &
                     (i_grid_shell.lt.n_ang_shells(i_species)) )
            i_grid_shell = i_grid_shell + 1
          enddo

          n_lebedev(i_radial,i_species) = grid_ceil(n_ang_points(i_grid_shell,i_species))

          n_angular(i_radial,i_species) = &
                n_angular_lebedev(n_lebedev(i_radial,i_species))
          r_angular(:,:,i_radial,i_species) = &
                r_angular_lebedev(:,:,n_lebedev(i_radial,i_species))
          w_angular(:,i_radial,i_species) = &
                w_angular_lebedev(:,n_lebedev(i_radial,i_species))

          if (n_angular(i_radial,i_species).gt.n_max_angular) then
            ! smallest grid is already too large, abort calculation
            write(use_unit,'(1X,A,I5,A,I5,A)') "* Species ", i_species, " radial shell ",    i_radial, ":"
            write(use_unit,'(1X,A)')           "* No suitable angular grid available. "
            write(use_unit,'(1X,A,A)')         "* Adjust the prescribed number of grid points ", &
                 "in control.in."
            stop
          end if

        else if (angular_acc(i_species).gt.0.d0) then

           ! get minimum number of Lebedev grid requested for each species

           n_lebedev(i_radial,i_species) = grid_ceil(angular_min(i_species))

!!$           call get_angular_grid &
!!$                ( n_lebedev(i_radial,i_species), &
!!$                n_angular(i_radial,i_species), &
!!$                r_angular(1,1,i_radial,i_species), &
!!$                w_angular(1,i_radial,i_species) &
!!$                )

           n_angular(i_radial,i_species) = &
                n_angular_lebedev(n_lebedev(i_radial,i_species))
           r_angular(:,:,i_radial,i_species) = &
                r_angular_lebedev(:,:,n_lebedev(i_radial,i_species))
           w_angular(:,i_radial,i_species) = &
                w_angular_lebedev(:,n_lebedev(i_radial,i_species))

           if (n_angular(i_radial,i_species).gt.n_max_angular) then
              ! smallest grid is already too large, abort calculation
              write(use_unit,'(1X,A,I5,A,I5,A)') "* Species ", i_species, " radial shell ",    i_radial, ":"
              write(use_unit,'(1X,A)')           "* No suitable angular grid available. "
              write(use_unit,'(1X,A,A)')         "* Adjust the prescribed number of grid points ", &
                   "in control.in."
              stop
           end if

        else

!!$           if(species_fixed_grid(i_species))then
!!$              read(88,*) angular_limit(i_species)
!!$           end if
!           if(species_fixed_grid(i_species))then
!              n_lebedev(i_radial,i_species) = fixed_grid_index(i_radial,i_species)
!           else
              n_lebedev(i_radial,i_species) = grid_floor(angular_limit(i_species))
!           end if

           n_angular(i_radial,i_species) = &
                n_angular_lebedev(n_lebedev(i_radial,i_species))
           r_angular(:,:,i_radial,i_species) = &
                r_angular_lebedev(:,:,n_lebedev(i_radial,i_species))
           w_angular(:,i_radial,i_species) = &
                w_angular_lebedev(:,n_lebedev(i_radial,i_species))

        end if

!!$        call divide_angular_grid( &
!!$             n_angular(i_radial, i_species), &
!!$             r_angular(1,1,i_radial,i_species), &
!!$             w_angular(1,i_radial,i_species), &
!!$             n_division(i_radial, i_species), &
!!$             division_boundaries(1,i_radial,i_species), &
!!$             r_radial(i_radial,i_species) &
!!$             )

        n_division(i_radial,i_species) = &
             n_division_lebedev(n_lebedev(i_radial,i_species))
        division_boundaries(:,i_radial,i_species) = &
             division_boundaries_lebedev(:,n_lebedev(i_radial,i_species))
     enddo

!!$     if(species_fixed_grid(i_species))then
!!$        close(88)
!!$     end if
  enddo

  ! perform partitioned integration, atom by atom, and point by point 
  ! This will be the outermost loop, to save evaluations of the Y_lm functions

  do i_atom = 1, n_atoms,1

     do i_radial = 1, n_radial(species(i_atom)), 1

        if (myid.eq.radial_task_list(i_radial, i_atom)) then

           local_gradient        = 0.d0
           local_potential_parts = 0.d0
           zora_potential_parts  = 0.d0

           if (angular_acc(species(i_atom)).eq.0.d0 .and.  .not. full_shell(i_radial, i_atom)) then


              ! We are using the default grid.
              ! Do not loop over Lebedev grids.
              ! The angular grid has already been set outside this subroutine,
              ! either in get_grids() or as an externally provided grid.

              do i_division = 1, n_division(i_radial,species(i_atom)), 1

                 ! find all non-zero wave functions in this angular shell

                 n_compute_a(i_division) = 0
                 n_compute_c(i_division) = 0
                 i_basis(:,i_division) = 0

                 i_point = 0
                 do i_angular = division_boundaries(i_division, &
                      i_radial,species(i_atom))+1, &
                      division_boundaries(i_division+1, &
                      i_radial,species(i_atom)), 1


                    ! get current integration point coordinate
                    do i_coord = 1, 3, 1
                       coord_current(i_coord, i_angular) = &
                            coords_center(i_coord,i_atom ) + &
                            r_angular(i_coord,i_angular, i_radial,species(i_atom)) * &
                            r_radial(i_radial, species(i_atom))
                    enddo



                    if(n_periodic == 0 )then                       


                       ! tabulate current integration point as it appears from spherical
                       ! coordinates centered at each atom
                       call tab_atom_centered_coords_p0 &
                            ( coord_current(1,i_angular),  &
                            dist_tab_sq(1,i_angular), &
                            dir_tab_full(1,1,i_angular), &
                            n_centers_integrals, centers_basis_integrals )


                       call tab_global_geometry_p0 &
                            ( dist_tab_sq(1,i_angular), &
                            dir_tab_full(1,1,i_angular), &
                            dist_tab_full(1,i_angular), &
                            i_r_full(1,i_angular), &
                            dir_tab_full_norm(1,1,i_angular), &
                            n_centers_integrals,  centers_basis_integrals)

                       dist_tab_sq_full(:,i_angular) =  dist_tab_sq(:,i_angular) 


                       ! evaluate partition_tab
                       call evaluate_partition_tab_p0 &
                            ( i_atom, i_atom, dist_tab_full(1,i_angular), &
                            i_r_full(1,i_angular), &
                            w_radial(i_radial, species(i_atom)), &
                            w_angular(i_angular, i_radial, species(i_atom)), &
                            partition_tab(i_angular), &
                            n_centers_basis_integrals, &
                            centers_basis_integrals )



                    else


                       ! tabulate current integration point as it appears from spherical
                       ! coordinates centered at each atom
                       call tab_atom_centered_coords_p0 &
                            ( coord_current(1,i_angular),  &
                            dist_tab_sq_par_tab(1), &
                            dir_tab_par_tab(1,1), &
                            n_centers_basis_integrals, centers_basis_integrals )


                       call tab_global_geometry_p0 &
                            ( dist_tab_sq_par_tab(1), &
                            dir_tab_par_tab(1,1), &
                            dist_tab_par_tab(1), &
                            i_r_par_tab(1), &
                            dir_tab_par_tab_norm(1,1), &
                            n_centers_basis_integrals, centers_basis_integrals )



                       ! evaluate partition_tab
                       call evaluate_partition_tab_p0 &
                            ( i_atom, i_atom, dist_tab_par_tab(1), &
                            i_r_par_tab(1), &
                            w_radial(i_radial, species(i_atom)), &
                            w_angular(i_angular, i_radial, species(i_atom)), &
                            partition_tab(i_angular), &
                            n_centers_basis_integrals, &
                            centers_basis_integrals )



                       call map_to_center_cell(coord_current(1:3, i_angular) )

                       call tab_atom_centered_coords_p0 &
                            ( coord_current(1,i_angular),  &
                            dist_tab_sq(1,i_angular), &
                            dir_tab_full(1,1,i_angular), &
                            n_centers_integrals, centers_basis_integrals )


                       call tab_global_geometry_p0 &
                            ( dist_tab_sq(1,i_angular), &
                            dir_tab_full(1,1,i_angular), &
                            dist_tab_full(1,i_angular), &
                            i_r_full(1,i_angular), &
                            dir_tab_full_norm(1,1,i_angular), &
                            n_centers_integrals,  centers_basis_integrals)

                       dist_tab_sq_full(:,i_angular) =  dist_tab_sq(:,i_angular) 


                    end if




                    ! execute only if partition_tab.gt.0 here, i.e. if the integration point
                    ! makes sense
                    if (partition_tab(i_angular).gt.0.d0)then

                       i_point = i_point + 1
                       ! determine which basis functions are relevant at current integration point,
                       ! and tabulate their indices
                       call prune_basis_p0(dist_tab_sq(1,i_angular), &
                            n_compute_a(i_division), n_compute_c(i_division), &
                            i_basis(1,i_division), &
                            n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals )



                    end if
                 enddo ! do i_angular .....

                 n_points(i_division) = i_point

                 !write(use_unit,*) i_radial, n_compute_c(i_division)



                 if (n_compute_a(i_division).gt.0) then




                    n_rel_points = 0
                    i_point = 0

                    ! collect all wave function components for this part of the angular shell
                    do i_angular = division_boundaries(i_division, &
                         i_radial,species(i_atom))+1, &
                         division_boundaries(i_division+1, &
                         i_radial,species(i_atom)), 1

                       ! execute only if partition_tab.gt.0 here, i.e. if the integration point
                       ! makes sense
                       if (partition_tab(i_angular).gt.0.d0) then

                          ! i_point counts the number of integration points which are actually needed
                          ! in later matrix multiplications
                          i_point = i_point + 1
                          partition(i_point) = partition_tab(i_angular)

                          n_compute_atoms = 0
                          n_compute_fns = 0
                          i_basis_fns_inv = 0

                          ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                          ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                          ! without any copying and without doing any unnecessary operations. 
                          ! The price is that the interface is no longer explicit in terms of physical 
                          ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                          dir_tab(:,:,1,1) = dir_tab_full(:,:,i_angular)

                          call prune_radial_basis_p0 &
                               ( dist_tab_sq(1,i_angular), & 
                               !                                  dist_tab(1,i_point,i_division), &
                               dist_tab(1,1,1), &
                               dir_tab(1,1,1,1), &
                               n_compute_atoms, &
                               atom_index,  &
                               atom_index_inv, &
                               n_compute_fns, &
                               i_basis_fns, &
                               i_basis_fns_inv, &
                               i_atom_fns, &
                               spline_array_start, &
                               spline_array_end, &
                               n_centers_integrals, centers_basis_integrals)


                          ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                          ! for all atoms which are actually relevant
                          ! We need to unload dir_tab from its full version since
                          ! only tab_local_geometry normalises its value


                          call tab_local_geometry_p0 &
                               ( dist_tab_sq(1, i_angular), &
                               n_compute_atoms, &
                               atom_index, &
                               dir_tab(1,1,1,1), &
                               dist_tab(1,1,1), &
                               i_r(1,1,1) &
                               )


                          ! tabulate trigonometric functions sin(theta), cos(theta), sin(phi), cos(phi)
                          ! of current integration point as seen from each atom
                          call tab_trigonom_p0 &
                               ( n_compute_atoms,  &
                               dir_tab(1,1,1,1), &
                               trigonom_tab(1,1,1,1) &
                               )


                          ! tabulate distance and Ylm's w.r.t. other atoms            
                          call tab_wave_ylm_p0 &
                               ( n_compute_atoms, &
                               atom_index, &
                               trigonom_tab(1,1,1,1), &
                               basis_l_max, &
                               l_ylm_max, ylm_tab(1,1,1,1))


                          ! Now evaluate radial functions, radial kinetic energy terms, and
                          ! possibly radial derivatives from the previously stored compressed spline arrays  
                          call evaluate_radial_functions_p0 &
                               (   &
                               spline_array_start, &
                               spline_array_end, &
                               n_compute_atoms, &
                               n_compute_fns, &
                               dist_tab(1,1,1), &
                               i_r(1,1,1), &
                               atom_index, &
                               i_basis_fns_inv, &
                               basis_wave_ordered, &
                               radial_wave(1,i_point,i_division), &
                               .false.,  n_compute_c(i_division), n_max_compute_fns_ham  &
                               )


                          ! tabulate total wave function value for each basis function
                          call evaluate_waves_p0 &
                               ( l_ylm_max,  &
                               ylm_tab(1,1,1,1), &
                               dist_tab(1,1,1), &
                               index_lm, n_compute_c(i_division), &
                               i_basis(1,i_division), &
                               radial_wave(1,i_point,i_division), &
                               wave(1,i_point,i_division), &
                               n_compute_atoms, &
                               atom_index_inv, &
                               n_compute_fns, &
                               i_basis_fns_inv, n_max_compute_fns_ham &
                               )



                          ! now tackle specifics to do with the potential

                             ! get the local potentials

                             if (force_potential.eq.1) then
                                ! only superposition of free-atom effective potentials

                                ! at present, we compute only spin-averaged superposition potential
                                call evaluate_pot_superpos_p0 &
                                     (  &
                                     i_r_full(1,i_angular), &
                                     local_potential_parts(1,i_angular), &
                                     n_centers_integrals, centers_basis_integrals )


                             else
                                ! normal case: Superposition of atomic densities defines initial potential
                                ! TODO: calculate the local_kinetic_density in this subroutine

                                call evaluate_free_atom_sums_p0 &
                                     ( &
                                     dist_tab_full(1,i_angular), &
                                     i_r_full(1,i_angular), &
                                     dir_tab_full_norm(1,1,i_angular), &
                                     free_hartree_superpos, &
                                     free_rho_superpos, &
                                     free_rho_gradient_superpos, &
                                     rho, &
                                     local_rho_gradient, n_centers_integrals, centers_basis_integrals &
                                     )

                                  call evaluate_xc &
                                       ( rho, &
                                         local_rho_gradient, &
                                         local_kinetic_density, &
                                         en_density_xc, &
                                         en_density_x, en_density_c, &
                                         local_xc_derivs, &
                                         local_xc_gradient_deriv, &
                                         local_xc_tau_deriv, &
                                         use_hartree_fock .or. use_meta_gga &
                                       )

                                do i_spin = 1, n_spin, 1
                                   local_potential_parts(i_spin,i_angular) = &
                                        free_hartree_superpos + &
                                        local_xc_derivs(i_spin)

                                   local_gradient(1:3,i_spin,i_angular) = &
                                        local_xc_gradient_deriv(1:3,i_spin)*4

                                enddo

                                ! Nadia: if (use_embedding_potential) then add pot_ion_embed_local
                                if (use_embedding_potential) then

                                   call embedding_potential &
                                        ( coord_current(1,i_angular), &
                                        pot_ion_embed_local, dummy_force &
                                        )

                                   do i_spin = 1, n_spin, 1
                                      local_potential_parts(i_spin,i_angular) = &
                                           local_potential_parts(i_spin,i_angular) + &
                                           pot_ion_embed_local
                                   end do
                                end if

                                if (use_embedding_pp) then
                                   call ion_pseudopot_potential_v1 &
                                        ( coord_current(1,i_angular), &
                                        pp_potential,0)

                                   do i_spin = 1, n_spin, 1
                                      local_potential_parts(i_spin,i_angular) = &
                                           local_potential_parts(i_spin,i_angular) + &
                                           pp_potential
                                   end do
                                end if

                             end if

                          ! check whether we need scalar relativistic terms ...
                          ! FIXME: With gradient functionals, we never have the entire 
                          ! Kohn-Sham potential available; the gradient terms are always 
                          ! missing because they are treated by integration by parts instead.
                          ! Hence, our "ZORA operator" can only be based on LDA-like terms. 
                          ! Hopefully this will be good enough for scalar relativity.
                          ! Alternatively, we could remember those (hopefully few) grid points where ZORA
                          ! is required, and calculate and Pulay-mix the density Hessian only there, so as
                          ! to compute the local potential explicitly for only those few points.

                          if (flag_rel.eq.1) then



                             call evaluate_pot_superpos_p0 &
                                  ( &
                                  i_r_full(1,i_angular), &
                                  zora_potential_parts(1,i_angular), &
                                  n_centers_integrals, centers_basis_integrals )


                             do i_spin = 1, n_spin, 1

                                zora_operator(i_spin) = &
                                     2.d0 * light_speed_sq / &
                                     (2 * light_speed_sq - &
                                     zora_potential_parts(i_spin,i_angular))
                             enddo

                             ! if zora required for one spin component, do it also for the other
                             t_zora = .false.
                             do i_spin = 1, n_spin, 1
                                t_zora = t_zora .or. &
                                     (abs(zora_operator(i_spin)-0.5d0) &
                                     .gt.zora_threshold)
                             enddo

                             t_zora = .true.

                          end if ! (flag_rel.eq.1)


                          if (use_gga .or. t_zora) then 
                             ! we require the gradient of each basis function

                             ! tabulate radial derivatives of those radial functions 
                             ! which are actually non-zero at current point, using vectorized splines


                             call evaluate_radial_functions_p0 &
                                  ( &
                                  spline_array_start, &
                                  spline_array_end, &
                                  n_compute_atoms, &
                                  n_compute_fns, &
                                  dist_tab(1,1,1), &
                                  i_r(1,1,1), &
                                  atom_index, &
                                  i_basis_fns_inv, &
                                  basis_deriv_ordered, &
                                  radial_wave_deriv, .true., &
                                  n_compute_c(i_division), n_max_compute_fns_ham )


                             ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                             call tab_gradient_ylm_p0 &
                                  ( trigonom_tab(1,1,1,1), &
                                  basis_l_max, &
                                  l_ylm_max, n_compute_atoms, &
                                  atom_index, &
                                  ylm_tab(1,1,1,1), &
                                  dylm_dtheta_tab(1,1), &
                                  scaled_dylm_dphi_tab(1,1) &
                                  )

                             call evaluate_wave_gradient_p0 &
                                  ( dist_tab(1,1,1), &
                                  dir_tab(1,1,1,1),  &
                                  trigonom_tab(1,1,1,1), &
                                  l_ylm_max, ylm_tab(1,1,1,1), &
                                  dylm_dtheta_tab(1,1), &
                                  scaled_dylm_dphi_tab(1,1), &
                                  index_lm, n_compute_c(i_division), &
                                  i_basis(1,i_division), &
                                  radial_wave(1,i_point,i_division), &
                                  radial_wave_deriv, &
                                  gradient_basis_wave, &
                                  n_compute_atoms, &
                                  atom_index_inv, &
                                  n_compute_fns, &
                                  i_basis_fns_inv, n_max_compute_fns_ham   )

                          end if

                          ! Now, evaluate vector of components H*phi(i,r)
                          ! Local potential parts first; in the case of GGA, 
                          ! the real gradient parts are added further below
                          !                   if (flag_rel == 0 ) then
                          !                   if (.not. t_zora) then
                          ! Non-relativistic treatment - simply evaluate 
                          ! H*phi(i,r) all in one

                          ! First, obtain radial kinetic energy terms from vectorized splines
                          call evaluate_radial_functions_p0 &
                               ( &
                               spline_array_start, &
                               spline_array_end,   &
                               n_compute_atoms, &
                               n_compute_fns, &
                               dist_tab(1,1,1), &
                               i_r(1,1,1), &
                               atom_index, &
                               i_basis_fns_inv, &
                               basis_kinetic_ordered, kinetic_wave, &
                               .false., n_compute_c(i_division), n_max_compute_fns_ham  &
                               )


                          do i_spin = 1, n_spin, 1
                             call evaluate_H_psi_p0 &
                                  ( l_ylm_max, &
                                  ylm_tab(1,1,1,1), &
                                  dist_tab(1, 1, 1), &
                                  index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), &
                                  H_times_psi(1,i_point, i_spin), &
                                  radial_wave(1,i_point,i_division), &
                                  local_potential_parts(i_spin,i_angular), &
                                  n_compute_c(i_division), &
                                  i_basis(1,i_division), &
                                  n_compute_atoms,  &
                                  atom_index_inv, &
                                  n_compute_fns, &
                                  i_basis_fns_inv, &
                                  kinetic_wave, &
                                  zora_operator(i_spin), n_max_compute_fns_ham &
                                  )
                          enddo


                          if ( t_zora ) then

                             n_rel_points = n_rel_points + 1

                             do i_spin = 1, n_spin, 1

                                zora_operator(i_spin) = &
                                     light_speed_sq / &
                                     (2 * light_speed_sq - &
                                     zora_potential_parts(i_spin,i_angular))**2

                                call  add_zora_gradient_part_p0( &
                                     local_gradient(1,i_spin,i_angular), &
                                     i_r_full(1,i_angular), &
                                     dir_tab_full_norm(1,1,i_angular),  &
                                     dist_tab_full(1,i_angular), &
                                     zora_operator(i_spin),  &
                                     n_centers_integrals, centers_basis_integrals  )

                             enddo
                          end if ! end ZORA preparations


                          ! If using a GGA, add the true gradient terms to the Hamiltonian vector
                          if (use_gga .or.n_rel_points.gt.0 ) then

                             do i_spin = 1, n_spin, 1

                                call add_gradient_part_to_H_p0 &
                                     ( n_compute_c(i_division), &
                                     gradient_basis_wave, &
                                     local_gradient(1,i_spin,i_angular), &
                                     H_times_psi(1, i_point,i_spin))
                             enddo
                          end if
                       end if ! (partition_tab >0)
                    enddo ! end angular integration loop for given grids


                    ! now increment overlap and Hamiltonian matrix

                    if (n_points(i_division).gt.0) then

                       ! add contributions to the overlap matrix elements

                       call evaluate_ovlp_shell_p0 &
                            ( n_points(i_division),  &
                            partition, &
                            n_compute_a(i_division),n_compute_c(i_division), &
                            wave(1,1,i_division), &
                            matrix_shell, n_max_compute_ham  &
                            )

                       call update_full_matrix_p0 &
                            ( n_compute_c(i_division), n_compute_a(i_division),  &
                            i_basis(1 ,i_division), &
                            matrix_shell,  &
                            overlap_matrix &
                            )


                       ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
                       ! work separately for each spin channel
                       do i_spin = 1, n_spin, 1

                          ! add full non-relativistic contributions and (for relativistic points)
                          ! all contributions from the potential to the Hamiltonian matrix elements
                          call evaluate_hamiltonian_shell_p1 &
                               ( n_points(i_division),  &
                               partition,  &
                               n_compute_c(i_division), &
                               H_times_psi(1,1,i_spin), &
                               n_max_compute_ham, wave(1,1,i_division),  &
                               matrix_shell  )


                          call update_full_matrix_p0 &
                               ( n_compute_c(i_division), n_compute_a(i_division), &
                               i_basis(1,i_division), &
                               matrix_shell,  &
                               hamiltonian(1,i_spin) &
                               )

                       enddo !i_spin = ....
                    end if !(n_poits >0)
                    ! finished incrementing overlap and Hamiltonian matrices 
                 end if ! (n_compute_a >0)
              enddo ! end loop over partitions of angular shells




              !pppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppppp

           else



              ! Angular grids will be optimized. Try successive Lebedev grids.

              ! set i_lebedev to current value for given radial shell and species
              i_lebedev = n_lebedev(i_radial,species(i_atom))

              ! preload current grid data
              try_n_angular = n_angular(i_radial, species(i_atom)) 
              do i_angular = 1,try_n_angular
                 do i_coord = 1,3,1
                    try_r_angular(i_coord, i_angular) = &
                         r_angular( i_coord, i_angular, &
                         i_radial, species(i_atom) )

                 enddo


                 try_w_angular(i_angular) = &
                      w_angular(i_angular, i_radial, species(i_atom))
              enddo

              try_n_division = n_division(i_radial, &
                   species(i_atom))
              do i_division = 1, try_n_division + 1, 1
                 try_division_boundaries(i_division) = &
                      division_boundaries(i_division, &
                      i_radial,species(i_atom))
              enddo

              first_loop = .true.
              shell_converged = .false.


              !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
              !
              ! Optimizing the grids using overlap matrix
              !
              !


              ! begin loop over Lebedev grids   
              do while &
                   ( .not.shell_converged &
                   )


                 new_shell(:,1) = 0.d0

                 do i_division = 1, try_n_division, 1

                    ! find again all non-zero wave functions in the shell
                    try_n_compute_a(i_division) = 0
                    try_n_compute_c(i_division) = 0
                    try_i_basis(:,i_division) = 0
                    i_point = 0


                    ! integrate all overlap matrix terms on present part of the shell (i_radial,i_atom)
                    do i_angular = try_division_boundaries(i_division)+1, &
                         try_division_boundaries(i_division+1), 1 

                       ! get current integration point coordinate
                       do i_coord = 1, 3, 1
                          try_coord(i_coord,i_angular) =  &
                               coords_center(i_coord,i_atom ) + &
                               try_r_angular(i_coord,i_angular) * &
                               r_radial(i_radial, species(i_atom))
                       enddo


                       if(n_periodic == 0)then

                          ! tabulate current integration point as it appears from spherical
                          ! coordinates centered at each atom
                          call tab_atom_centered_coords_p0 &
                               ( try_coord(1,i_angular), &
                               try_dist_tab_sq(1,i_angular), &
                               try_dir_tab_full(1,1,i_angular),&
                               n_centers_integrals, centers_basis_integrals )


                          call tab_global_geometry_p0 &
                               ( try_dist_tab_sq(1,i_angular), &
                               try_dir_tab_full(1,1,i_angular), &
                               try_dist_tab_full(1,i_angular), &
                               try_i_r_full(1,i_angular), &
                               try_dir_tab_full_norm(1,1,i_angular), &
                               n_centers_integrals,  centers_basis_integrals)

                          try_dist_tab_sq_full(:,i_angular) = try_dist_tab_sq(:,i_angular)

                          ! evaluate partition_tab
                          call evaluate_partition_tab_p0 &
                               ( i_atom, i_atom, try_dist_tab_full(1, i_angular), &
                               try_i_r_full(1, i_angular), &
                               w_radial(i_radial, species(i_atom)), &
                               try_w_angular(i_angular), &
                               try_partition_tab(i_angular), &
                               n_centers_basis_integrals, centers_basis_integrals )


                       else




                          ! tabulate current integration point as it appears from spherical
                          ! coordinates centered at each atom
                          call tab_atom_centered_coords_p0 &
                               ( try_coord(1,i_angular),  &
                               dist_tab_sq_par_tab(1), &
                               dir_tab_par_tab(1,1), &
                               n_centers_basis_integrals, centers_basis_integrals )


                          call tab_global_geometry_p0 &
                               ( dist_tab_sq_par_tab(1), &
                               dir_tab_par_tab(1,1), &
                               dist_tab_par_tab(1), &
                               i_r_par_tab(1), &
                               dir_tab_par_tab_norm(1,1), &
                               n_centers_basis_integrals, centers_basis_integrals )



                          ! evaluate partition_tab
                          call evaluate_partition_tab_p0 &
                               ( i_atom, i_atom, dist_tab_par_tab(1), &
                               i_r_par_tab(1), &
                               w_radial(i_radial, species(i_atom)), &
                               try_w_angular(i_angular), &
                               try_partition_tab(i_angular), &
                               n_centers_basis_integrals, &
                               centers_basis_integrals )



                          call map_to_center_cell(try_coord(1:3,i_angular))

                          call tab_atom_centered_coords_p0 &
                               ( try_coord(1,i_angular), &
                               try_dist_tab_sq(1,i_angular), &
                               try_dir_tab_full(1,1,i_angular),&
                               n_centers_integrals, centers_basis_integrals )

                          call tab_global_geometry_p0 &
                               ( try_dist_tab_sq(1,i_angular), &
                               try_dir_tab_full(1,1,i_angular), &
                               try_dist_tab_full(1,i_angular), &
                               try_i_r_full(1,i_angular), &
                               try_dir_tab_full_norm(1,1,i_angular), &
                               n_centers_integrals,  centers_basis_integrals)

                          try_dist_tab_sq_full(:,i_angular) = try_dist_tab_sq(:,i_angular)


                       end if


                       ! execute only if try_partition_tab .gt. 0 and determine
                       ! relevant basis functions
                       if (try_partition_tab(i_angular) .gt. 0.0d0) &
                            then

                          i_point = i_point + 1

                          call prune_basis_p0( &
                               try_dist_tab_sq(1,i_angular), &
                               try_n_compute_a(i_division), try_n_compute_c(i_division),  &
                               try_i_basis(1,i_division), &
                               n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals )
                       end if
                    enddo



                    try_n_points(i_division) = i_point



                    ! collect the wave functions in this shell

                    if (try_n_compute_a(i_division) .gt. 0) then

                       i_point = 0

                       do i_angular = try_division_boundaries( &
                            i_division)+1, &
                            try_division_boundaries(i_division+1), 1 


                          ! execute only if try_partition_tab.gt.0 here, i.e. if the integration point
                          ! makes sense
                          if (try_partition_tab(i_angular).gt.0.d0) then

                             ! i_point counts the integration points in this shell
                             i_point = i_point + 1
                             try_partition(i_point) = &
                                  try_partition_tab(i_angular)

                             ! tabulate trigonometric functions sin(theta), cos(theta), sin(phi), cos(phi)
                             ! of current integration point as seen from each atom

                             n_compute_atoms = 0
                             n_compute_fns = 0
                             i_basis_fns_inv = 0

                             ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                             ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                             ! without any copying and without doing any unnecessary operations. 
                             ! The price is that the interface is no longer explicit in terms of physical 
                             ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                             try_dir_tab(:,:,i_point,i_division) = try_dir_tab_full(:,:,i_angular)

                             call prune_radial_basis_p0 &
                                  ( try_dist_tab_sq(1,i_angular), &
                                  try_dist_tab(1,i_point,i_division), &
                                  try_dir_tab(1,1,i_point,i_division), &
                                  n_compute_atoms, &
                                  atom_index, &
                                  atom_index_inv, &
                                  n_compute_fns, &
                                  i_basis_fns, &
                                  i_basis_fns_inv, &
                                  i_atom_fns, &
                                  spline_array_start, &
                                  spline_array_end, &
                                  n_centers_integrals, centers_basis_integrals)


                             ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                             ! for all atoms which are actually relevant
                             ! We need to unload dir_tab from its full version since
                             ! only tab_local_geometry normalises its value

                             call tab_local_geometry_p0 &
                                  ( try_dist_tab_sq(1,i_angular), &
                                  n_compute_atoms, &
                                  atom_index, &
                                  try_dir_tab(1,1,i_point,i_division), &
                                  try_dist_tab(1,i_point,i_division),  &
                                  try_i_r(1,i_point,i_division) &
                                  )


                             ! tabulate trigonometric functions sin(theta), cos(theta), sin(phi), cos(phi)
                             ! of current integration point as seen from each atom
                             call tab_trigonom_p0 &
                                  ( n_compute_atoms, &
                                  try_dir_tab(1,1,i_point,i_division), &
                                  try_trigonom_tab(1,1,i_point,i_division) &
                                  )


                             ! tabulate distance and Ylm's w.r.t. other atoms            
                             call tab_wave_ylm_p0 &
                                  ( n_compute_atoms, &
                                  atom_index,  &
                                  try_trigonom_tab(1,1,i_point,i_division), &
                                  basis_l_max, l_ylm_max,  &
                                  try_ylm_tab(1,1,i_point,i_division) )


                             ! Now evaluate radial functions, radial kinetic energy terms, and
                             ! possibly radial derivatives from the previously stored compressed spline arrays  
                             call evaluate_radial_functions_p0 &
                                  ( spline_array_start, &
                                  spline_array_end, &
                                  n_compute_atoms,  &
                                  n_compute_fns, &
                                  try_dist_tab(1,i_point,i_division), &
                                  try_i_r(1,i_point,i_division), &
                                  atom_index,  &
                                  i_basis_fns_inv, &
                                  basis_wave_ordered, &
                                  try_radial_wave(1,i_point,i_division), &
                                  .false., try_n_compute_c(i_division), n_max_compute_fns_ham  & 
                                  )

                             ! tabulate total wave function value for each basis function
                             call evaluate_waves_p0 &
                                  ( l_ylm_max, &
                                  try_ylm_tab(1,1,i_point,i_division), &
                                  try_dist_tab(1,i_point,i_division),  &
                                  index_lm, try_n_compute_c(i_division), &
                                  try_i_basis(1,i_division), &
                                  try_radial_wave(1,i_point,i_division), &
                                  try_wave(1,i_point,i_division), &
                                  n_compute_atoms, &
                                  atom_index_inv, &
                                  n_compute_fns, &
                                  i_basis_fns_inv, n_max_compute_fns_ham  &
                                  )

                          end if !(try_partition_tab >0)
                       enddo  ! end angular integration loop


                       ! add contributions to the overlap matrix elements

                       call evaluate_ovlp_shell_p0 &
                            ( try_n_points(i_division), &
                            try_partition, &
                            try_n_compute_a(i_division), try_n_compute_c(i_division), &
                            try_wave(1,1,i_division), &
                            matrix_shell, n_max_compute_ham  &
                            )

                    endif !(try_n_compute_a > 0)


                    ! Now we add the contributions from matrix_shell, which is the contribution
                    ! from the current division of an angular shell, to new_shell, which is
                    ! the contribution from the whole current integration shell
                    call update_full_matrix_p0 &
                         ( try_n_compute_c(i_division), try_n_compute_a(i_division), &
                         try_i_basis(1,i_division), &
                         matrix_shell,  &
                         new_shell(1,1) &
                         )

                 enddo ! end partitions of the angular shell


                 !  check convergence
                 !  To do this, we need at least two loops.
                 if (.not.first_loop) then

                    ! Initialize shell_convergence
                    shell_converged = .true.

                    !  use the maximum deviation of a single overlap matrix term 
                    !  as the convergence criterion
                    call check_shell_convergence_p0 &
                         ( prev_shell(:,1), new_shell(:,1),  &
                         angular_acc(species(i_atom)), shell_converged )

                 else
                    first_loop = .false.
                 end if

                 if (.not.shell_converged) then
                    !  store present angular grid data for potential future use
                    !  also check (below) whether shell is actually converged according to other data

                    if ( n_angular(i_radial, species(i_atom)).lt.try_n_angular ) then

                       !  only store species-dependent data if grid changed

                       n_lebedev(i_radial,species(i_atom)) = i_lebedev
                       n_angular(i_radial,species(i_atom)) = try_n_angular
                       do i_angular = 1,n_angular(i_radial,species(i_atom))
                          do i_coord = 1,3,1
                             r_angular( i_coord, i_angular,  &
                                  i_radial, species(i_atom) ) =  &
                                  try_r_angular(i_coord, i_angular)
                          enddo
                          w_angular(i_angular, i_radial, species(i_atom)) =  &
                               try_w_angular(i_angular)
                       enddo

                       n_division(i_radial, species(i_atom)) = &
                            try_n_division
                       do i_division = 1, try_n_division + 1, 1 
                          division_boundaries(i_division, &
                               i_radial,species(i_atom)) = &
                               try_division_boundaries(i_division)
                       enddo
                    end if ! ( n_angular(i_radial, species(i_atom)) 



                    ! store atom-dependent data in any case.
                    do i_angular = 1,n_angular(i_radial,species(i_atom))
                       partition_tab(i_angular) = try_partition_tab(i_angular)
                    end do

                    ! store grid geometry and wave function components to be reused for Hamiltonian, below

                    coord_current     = try_coord
                    dist_tab_sq_full  = try_dist_tab_sq_full
                    dist_tab_sq       = try_dist_tab_sq
                    dist_tab          = try_dist_tab
                    dist_tab_full     = try_dist_tab_full
                    i_r               = try_i_r
                    i_r_full          = try_i_r_full
                    dir_tab           = try_dir_tab
                    dir_tab_full      = try_dir_tab_full
                    dir_tab_full_norm = try_dir_tab_full_norm
                    trigonom_tab      = try_trigonom_tab
                    ylm_tab           = try_ylm_tab
                    wave              = try_wave
                    n_compute_a       = try_n_compute_a
                    n_compute_c       = try_n_compute_c
                    i_basis           = try_i_basis
                    radial_wave       = try_radial_wave

                    ! spline_array_start = try_spline_array_start
                    ! spline_array_end = try_spline_array_end
                    ! n_compute_fns = try_n_compute_fns
                    ! n_compute_atoms = try_n_compute_atoms
                    ! i_atom_fns = try_i_atom_fns
                    ! i_basis_fns = try_i_basis_fns
                    ! i_basis_fns_inv = try_i_basis_fns_inv
                    ! atom_index = try_atom_index
                    ! atom_index_inv = try_atom_index_inv

                    ! store overlap matrix terms from preceding shell
                    prev_shell(:,1) = new_shell(:,1)

                    ! try next Lebedev grid
                    i_lebedev = i_lebedev + 1

                    if (i_lebedev.le.n_grids) then
                       ! get next angular grid and grid weights
!!$                          call get_angular_grid &
!!$                               ( i_lebedev, try_n_angular, try_r_angular,  &
!!$                               try_w_angular &
!!$                               )

                       try_n_angular = n_angular_lebedev(i_lebedev)
                       try_r_angular = r_angular_lebedev(:,:,i_lebedev)
                       try_w_angular = w_angular_lebedev(:,i_lebedev)

                       if (try_n_angular.gt.angular_limit(species(i_atom)) .or. full_shell(i_radial,i_atom))then

                          ! maximum allowed grid size exceeded, shell is "converged"
                          shell_converged = .true.
                          full_shell(i_radial,i_atom) = .true.

                       else
                          ! partition the grid
!!$                             call divide_angular_grid( &
!!$                                  try_n_angular, &
!!$                                  try_r_angular,  &
!!$                                  try_w_angular, &
!!$                                  try_n_division, &
!!$                                  try_division_boundaries, &
!!$                                  r_radial(i_radial, &
!!$                                  species(i_atom)) &
!!$                                  )

                          try_n_division = n_division_lebedev(i_lebedev)
                          try_division_boundaries = division_boundaries_lebedev(:,i_lebedev)

                       end if

                    else

                       ! if we're out of Lebedev grids to try, shell is "converged" anyway
                       shell_converged = .true.
                       full_shell(i_radial, i_atom) = .true.

                    end if
                 end if !(.not. shell_converged)
              enddo  ! end loop over possible Lebedev grids.



              ! increment overlap_matrix for present radial shell

              overlap_matrix =  overlap_matrix +  new_shell(:,1)



              !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
              !          
              !     Next, check the convergence of the Hamiltonian matrix for the present Lebedev grids

              !     bootstrap: compute Hamiltonian elements for the present chosen grid



              prev_shell = 0.d0

              do i_division = 1, n_division(i_radial,species(i_atom)), 1

                 i_point = 0
                 n_rel_points = 0

                 if (n_compute_a(i_division) .gt. 0) then

                    do i_angular = division_boundaries( &
                         i_division,i_radial,species(i_atom))+1, &
                         division_boundaries(i_division+1,i_radial, &
                         species(i_atom)), 1 

                       if (partition_tab(i_angular).gt.0.d0) then

                          ! i_point counts the number of integration points in this batch
                          i_point = i_point + 1
                          partition(i_point) = partition_tab(i_angular)

                          ! Pruning of radial functions will be repeated here.
                          ! We already did precisely the same operation for the ovlp. matrix,
                          ! but repeating it here is cheaper than storing all the index
                          ! arrays across the current integration shell.

                          n_compute_atoms = 0
                          n_compute_fns = 0
                          i_basis_fns_inv = 0

                          ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                          ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                          ! without any copying and without doing any unnecessary operations. 
                          ! The price is that the interface is no longer explicit in terms of physical 
                          ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.


                          dir_tab(:,:,i_point,i_division) = &
                               dir_tab_full(:,:,i_angular)
                          dist_tab_sq(:,i_angular) = dist_tab_sq_full(:,i_angular)



                          call prune_radial_basis_p0 &
                               ( dist_tab_sq(1,i_angular),  &
                               dist_tab(1,i_point,i_division), &
                               dir_tab(1,1,i_point,i_division), &
                               n_compute_atoms,  &
                               atom_index,  &
                               atom_index_inv, &
                               n_compute_fns,  &
                               i_basis_fns,  &
                               i_basis_fns_inv, &
                               i_atom_fns,  &
                               spline_array_start,  &
                               spline_array_end, &
                               n_centers_integrals, centers_basis_integrals)



                          ! evaluate initial potential at current integration point


                          if (force_potential.eq.1) then

                             ! only superposition of free-atom effective potentials

                             call evaluate_pot_superpos_p0 &
                                  ( i_r_full(1,i_angular),  &
                                  local_potential_parts(1,i_angular), &
                                  n_centers_integrals, centers_basis_integrals )

                          else  ! Normal case

                             ! TODO: Evaluate local_kinetic_density in this subroutine
 
                             call evaluate_free_atom_sums_p0 &
                                  ( dist_tab_full(1,i_angular),  &
                                  i_r_full(1,i_angular),  &
                                  dir_tab_full_norm(1,1,i_angular), &
                                  free_hartree_superpos, &
                                  free_rho_superpos, &
                                  free_rho_gradient_superpos, &
                                  rho, &
                                  local_rho_gradient, &
                                  n_centers_integrals, centers_basis_integrals )
  
                                  call evaluate_xc &
                                       ( rho, &
                                         local_rho_gradient, &
                                         local_kinetic_density, &
                                         en_density_xc, &
                                         en_density_x, en_density_c, &
                                         local_xc_derivs, &
                                         local_xc_gradient_deriv, &
                                         local_xc_tau_deriv, &
                                         use_hartree_fock .or. use_meta_gga &
                                       )

                             do i_spin = 1, n_spin, 1
                                local_potential_parts(i_spin,i_angular) =  &
                                     free_hartree_superpos + &
                                     local_xc_derivs(i_spin)

                                local_gradient(1:3,i_spin,i_angular) =  &
                                     4*local_xc_gradient_deriv(1:3,i_spin)
                             enddo


                             ! Nadia: if (use_embedding_potential) then add pot_ion_embed_local
                             if (use_embedding_potential) then

                                call embedding_potential &
                                     ( coord_current(1,i_angular),  &
                                     pot_ion_embed_local, dummy_force &
                                     )

                                do i_spin = 1, n_spin, 1
                                   local_potential_parts(i_spin,i_angular) =  &
                                        local_potential_parts(i_spin,i_angular) + &
                                        pot_ion_embed_local
                                end do
                             end if


                             if (use_embedding_pp) then
                                call ion_pseudopot_potential_v1 &
                                     ( coord_current(1,i_angular), &
                                     pp_potential,0)

                                do i_spin = 1, n_spin, 1
                                    local_potential_parts(i_spin,i_angular) = &
                                      local_potential_parts(i_spin,i_angular) + &
                                      pp_potential
                                end do

                             end if



                          end if  !  end evaluation of initial potential




                          ! check whether we need scalar relativistic terms ...
                          ! FIXME: With gradient functionals, we never have the entire 
                          ! Kohn-Sham potential available; the gradient terms are always 
                          ! missing because they are treated by integration by parts instead.
                          ! Hence, our "ZORA operator" can only be based on LDA-like terms. 
                          ! Hopefully this will be good enough for scalar relativity.
                          ! Alternatively we could remember those (hopefully few) grid points where ZORA
                          ! is required, and calculate and Pulay-mix the density Hessian only there, so as
                          ! to compute the local potential explicitly for only those few points.

                          if (flag_rel.eq.1) then


                             call evaluate_pot_superpos_p0 &
                                  (  &
                                  i_r_full(1,i_angular),  &
                                  zora_potential_parts(1,i_angular), &
                                  n_centers_integrals, centers_basis_integrals )


                             do i_spin = 1, n_spin, 1

                                zora_operator(i_spin) = &
                                     2.d0 * light_speed_sq / &
                                     ( 2 * light_speed_sq - &
                                     zora_potential_parts(i_spin,i_angular) )
                             enddo

                             ! if zora required for one spin component, do it also for the other
                             t_zora = .false.
                             do i_spin = 1, n_spin, 1
                                t_zora = t_zora .or. &
                                     (abs(zora_operator(i_spin)-0.5d0) &
                                     .gt.zora_threshold)
                             enddo
                             t_zora = .true.
                          end if !(flag_rel.eq.1)



                          if ((use_gga).or.t_zora) then

                             call evaluate_radial_functions_p0 &
                                  ( spline_array_start,  &
                                  spline_array_end, &
                                  n_compute_atoms,  &
                                  n_compute_fns,  &
                                  dist_tab(1,i_point,i_division),  &
                                  i_r(1,i_point,i_division), &
                                  atom_index,  &
                                  i_basis_fns_inv, &
                                  basis_deriv_ordered,  &
                                  radial_wave_deriv, .true., &
                                  n_compute_c(i_division), n_max_compute_fns_ham &
                                  )

                             call tab_gradient_ylm_p0 &
                                  ( trigonom_tab(1,1,i_point,i_division),  &
                                  basis_l_max,  &
                                  l_ylm_max, n_compute_atoms,  &
                                  atom_index, &
                                  ylm_tab(1,1,i_point,i_division),  &
                                  dylm_dtheta_tab(1,1),  &
                                  scaled_dylm_dphi_tab(1,1) &
                                  )

                             call evaluate_wave_gradient_p0 &
                                  ( dist_tab(1,i_point,i_division), &
                                  dir_tab(1,1,i_point,i_division), &
                                  trigonom_tab(1,1,i_point,i_division), &
                                  l_ylm_max, ylm_tab(1,1,i_point,i_division), &
                                  dylm_dtheta_tab(1,1), &
                                  scaled_dylm_dphi_tab(1,1), &
                                  index_lm, n_compute_c(i_division), &
                                  i_basis(1,i_division), &
                                  radial_wave(1, i_point, i_division), &
                                  radial_wave_deriv, &
                                  gradient_basis_wave, &
                                  n_compute_atoms, &
                                  atom_index_inv, &
                                  n_compute_fns, &
                                  i_basis_fns_inv, n_max_compute_fns_ham )
                          end if

                          ! for the purpose of grid adaptation, use both spin components.
                          ! this could be done with only one spin component, but then we
                          ! create all sorts of corner cases (like a single H atom which
                          ! is forced to be fully spin-polarized spin down)

                          ! Now, evaluate vector of components H*phi(i,r)
                          ! Local potential parts first; in the case of GGA, 
                          ! the real gradient parts are added further below


                          call evaluate_radial_functions_p0 &
                               ( spline_array_start,  &
                               spline_array_end, &
                               n_compute_atoms,  &
                               n_compute_fns,  &
                               dist_tab(1,i_point,i_division),  &
                               i_r(1,i_point,i_division), &
                               atom_index,  &
                               i_basis_fns_inv, &
                               basis_kinetic_ordered, kinetic_wave, &
                               .false., n_compute_c(i_division), n_max_compute_fns_ham &
                               )

                          do i_spin = 1, n_spin, 1

                             call evaluate_H_psi_p0 &
                                  ( l_ylm_max, &
                                  ylm_tab(1,1,i_point,i_division), &
                                  dist_tab(1,i_point,i_division), &
                                  index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), &
                                  H_times_psi(1, i_point, i_spin), &
                                  radial_wave(1,i_point,i_division), &
                                  local_potential_parts(i_spin,i_angular), &
                                  n_compute_c(i_division), &
                                  i_basis(1,i_division), &
                                  n_compute_atoms,  &
                                  atom_index_inv, &
                                  n_compute_fns,  &
                                  i_basis_fns_inv,  &
                                  kinetic_wave, &
                                  zora_operator(i_spin), n_max_compute_fns_ham &
                                  )
                          enddo



                          if(t_zora)then
                             ! ZORA scalar relativistic treatment required

                             n_rel_points = n_rel_points + 1

                             do i_spin = 1, n_spin, 1

                                zora_operator(i_spin) = &
                                     light_speed_sq / &
                                     (2 * light_speed_sq - &
                                     zora_potential_parts(i_spin,i_angular))**2

                                call  add_zora_gradient_part_p0(  &
                                     local_gradient(1,i_spin,i_angular), &
                                     i_r_full(1,i_angular),  &
                                     dir_tab_full_norm(1,1,i_angular),  &
                                     dist_tab_full(1,i_angular), &
                                     zora_operator(i_spin), &
                                     n_centers_integrals, centers_basis_integrals  )
                             end do
                          end if ! end ZORA preparations

                          if (use_gga .or. n_rel_points.gt.0 ) then

                             do i_spin = 1, n_spin, 1

                                !                                  local_gradient(:,i_spin,i_angular) = 0.0

                                call add_gradient_part_to_H_p0 &
                                     ( n_compute_c(i_division),  &
                                     gradient_basis_wave, &
                                     local_gradient(1,i_spin,i_angular), &
                                     H_times_psi(1, i_point,i_spin) )
                             enddo
                          end if
                       end if ! end if partition_tab.gt.0.d0 
                    enddo  ! end "i_angular" loop over one division of the angular shell



                    n_points(i_division) = i_point

                    ! evaluate the hamiltonian shell for this batch of n_points

                    if (n_points(i_division).gt.0) then

                       do i_spin = 1, n_spin, 1         

                          call evaluate_hamiltonian_shell_p1 &
                               ( n_points(i_division),  &
                               partition,  &
                               n_compute_c(i_division), &
                               H_times_psi(1,1,i_spin), &
                               n_max_compute_ham, wave(1,1,i_division),  &
                               matrix_shell  )

                          call update_full_matrix_p0 &
                               ( n_compute_c(i_division), n_compute_a(i_division), &
                               i_basis(1,i_division),  &
                               matrix_shell,   &
                               prev_shell(1,i_spin) &
                               )

                       enddo
                    end if
                 end if  !  end if (n_compute .gt 0)
              enddo !   enddo over divisions of the angular shell



              ! p  p  p  p  p  p p  p  p  p  p  p p  p  p  p  p  p p  p  p  p  p  p p  p  p  p  p  p p  p  p  p  p  p p  p  p  p  p  p

              if(n_periodic >0)then

                 try_local_potential_parts(:,:,1) = local_potential_parts(:,:)
                 try_zora_potential_parts(:,:,1)  = zora_potential_parts(:,:)
                 try_local_gradient(:,:,:,1)      = local_gradient(:,:,:)

              end if

              !aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa

              ! Check convergence of Hamiltonian unless there are no denser grids available

              !                 full_shell(i_radial, i_atom) = .true.

              if (.not.full_shell(i_radial, i_atom)) then

                 ! i_lebedev, try_grids etc are already set from overlap_matrix evaluation
                 ! can simply continue with last test grid and see how it performs for the
                 ! Hamiltonian ...


                 first_loop = .true.
                 shell_converged = .false.


                 ! begin loop over Lebedev grids   
                 do while ( .not.shell_converged )

                    new_shell = 0.d0

                    do i_division = 1, try_n_division, 1

                       if (.not.first_loop) then
                          ! determine try_n_compute, try_i_basis using prune_basis

                          try_n_compute_a(i_division) = 0
                          try_n_compute_c(i_division) = 0
                          try_i_basis(:,i_division) = 0

                          ! prune shell division (i_radial,i_atom)
                          do i_angular = try_division_boundaries( i_division)+1, &
                               try_division_boundaries(i_division+1), 1 

                             ! get current integration point coordinate
                             do i_coord = 1, 3, 1
                                try_coord(i_coord,i_angular) =  &
                                     coords_center(i_coord,i_atom ) + &
                                     try_r_angular(i_coord,i_angular) * &
                                     r_radial(i_radial, species(i_atom))
                             enddo

                             if(n_periodic == 0)then
                                ! tabulate current integration point as it appears from spherical
                                ! coordinates centered at each atom
                                call tab_atom_centered_coords_p0 &
                                     ( try_coord(1,i_angular),  &
                                     try_dist_tab_sq(1,i_angular), &
                                     try_dir_tab_full(1,1,i_angular), &
                                     n_centers_integrals, centers_basis_integrals )


                                call tab_global_geometry_p0 &
                                     ( try_dist_tab_sq(1,i_angular), &
                                     try_dir_tab_full(1,1,i_angular), &
                                     try_dist_tab_full(1,i_angular), &
                                     try_i_r_full(1,i_angular), &
                                     try_dir_tab_full_norm(1,1,i_angular), &
                                     n_centers_integrals,  centers_basis_integrals)

                                try_dist_tab_sq_full(:,i_angular) = try_dist_tab_sq(:,i_angular)


                                call evaluate_partition_tab_p0 &
                                     ( i_atom, i_atom, try_dist_tab_full(1, i_angular),  &
                                     try_i_r_full(1, i_angular),  &
                                     w_radial(i_radial, species(i_atom)), &
                                     try_w_angular(i_angular), &
                                     try_partition_tab(i_angular), &
                                     n_centers_basis_integrals, &
                                     centers_basis_integrals )

                             else 



                                ! coordinates centered at each atom
                                call tab_atom_centered_coords_p0 &
                                     ( try_coord(1,i_angular),  &
                                     dist_tab_sq_par_tab(1), &
                                     dir_tab_par_tab(1,1), &
                                     n_centers_basis_integrals, centers_basis_integrals )


                                call tab_global_geometry_p0 &
                                     ( dist_tab_sq_par_tab(1), &
                                     dir_tab_par_tab(1,1), &
                                     dist_tab_par_tab(1), &
                                     i_r_par_tab(1), &
                                     dir_tab_par_tab_norm(1,1), &
                                     n_centers_basis_integrals, centers_basis_integrals )



                                ! evaluate partition_tab
                                call evaluate_partition_tab_p0 &
                                     ( i_atom, i_atom, dist_tab_par_tab(1), &
                                     i_r_par_tab(1), &
                                     w_radial(i_radial, species(i_atom)), &
                                     try_w_angular(i_angular), &
                                     try_partition_tab(i_angular), &
                                     n_centers_basis_integrals, &
                                     centers_basis_integrals )


                                call map_to_center_cell(try_coord(1:3,i_angular)  )

                                call tab_atom_centered_coords_p0 &
                                     ( try_coord(1,i_angular),  &
                                     try_dist_tab_sq(1,i_angular), &
                                     try_dir_tab_full(1,1,i_angular), &
                                     n_centers_integrals, centers_basis_integrals )

                                call tab_global_geometry_p0 &
                                     ( try_dist_tab_sq(1,i_angular), &
                                     try_dir_tab_full(1,1,i_angular), &
                                     try_dist_tab_full(1,i_angular), &
                                     try_i_r_full(1,i_angular), &
                                     try_dir_tab_full_norm(1,1,i_angular), &
                                     n_centers_integrals,  centers_basis_integrals)

                                try_dist_tab_sq_full(:,i_angular) = try_dist_tab_sq(:,i_angular)

                             end if






!!$
!!$                                else
!!$                                   try_partition_tab(i_angular) = 0.d0
!!$                                end if


                             ! execute only if try_partition_tab .gt. 0 and determine
                             ! relevant basis functions
                             if (try_partition_tab(i_angular) .gt. 0.0d0)then

                                i_point = i_point + 1

                                call prune_basis_p0( &
                                     try_dist_tab_sq(1,i_angular), &
                                     try_n_compute_a(i_division), try_n_compute_c(i_division), &
                                     try_i_basis(1,i_division), &
                                     n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals ) 
                             end if
                          enddo ! i_angular = try_division_boundaries.....

                          try_n_points(i_division) = i_point

                       end if ! (.not.first_loop)



                       ! reinitialize Hamiltonian matrix terms for this shell
                       try_n_points(i_division) = 0
                       try_n_rel_points = 0

                       if (try_n_compute_a(i_division) .gt. 0) then

                          i_point = 0

                          do i_angular = try_division_boundaries( &
                               i_division)+1, &
                               try_division_boundaries(i_division+1), 1 

                             ! now compute next term in Hamiltonian integral

                             ! execute only if partition_tab.gt.0 here, i.e. if the integration point
                             ! makes sense
                             if (try_partition_tab(i_angular).gt.0.d0) then

                                i_point = i_point + 1
                                try_partition(i_point) = try_partition_tab(i_angular)

                                ! At this point, simply re-evaluate the pruning of
                                ! radial grids, on a per-point basis. In the first loop,
                                ! we would have this already through the overlap matrix integral before,
                                ! but would also have to store everything across the grid ...

                                n_compute_atoms = 0
                                n_compute_fns = 0
                                i_basis_fns_inv = 0

                                ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                                ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                                ! without any copying and without doing any unnecessary operations. 
                                ! The price is that the interface is no longer explicit in terms of physical 
                                ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.


                                try_dir_tab(:,:,i_point,i_division) = try_dir_tab_full(:,:,i_angular)
                                try_dist_tab_sq(:,i_angular)        = try_dist_tab_sq_full(:,i_angular)


                                call prune_radial_basis_p0 &
                                     ( try_dist_tab_sq(1,i_angular),  &
                                     try_dist_tab(1,i_point,i_division), &
                                     try_dir_tab(1,1,i_point,i_division), &
                                     n_compute_atoms,  &
                                     atom_index,  &
                                     atom_index_inv, &
                                     n_compute_fns,  &
                                     i_basis_fns,  &
                                     i_basis_fns_inv, &
                                     i_atom_fns,  &
                                     spline_array_start,  &
                                     spline_array_end, n_centers_integrals, centers_basis_integrals )


                                if (.not.first_loop) then
                                   ! tabulate trigonometric functions sin(theta), cos(theta), sin(phi), cos(phi)
                                   ! of current integration point as seen from each atom

                                   ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                                   ! for all atoms which are actually relevant
                                   ! We need to unload dir_tab from its full version since
                                   ! only tab_local_geometry normalises its value

                                   call tab_local_geometry_p0 &
                                        ( try_dist_tab_sq(1, i_angular),  &
                                        n_compute_atoms,  &
                                        atom_index, &
                                        try_dir_tab(1,1,i_point,i_division),  &
                                        try_dist_tab(1,i_point,i_division),  &
                                        try_i_r(1,i_point,i_division) &
                                        )

                                   ! tabulate trigonometric functions sin(theta), cos(theta), sin(phi), cos(phi)
                                   ! of current integration point as seen from each atom
                                   call tab_trigonom_p0 &
                                        ( n_compute_atoms,  &
                                        try_dir_tab(1,1,i_point,i_division),  &
                                        try_trigonom_tab(1,1,i_point,i_division) &
                                        )

                                   ! tabulate distance and Ylm's w.r.t. other atoms            
                                   call tab_wave_ylm_p0 &
                                        ( n_compute_atoms,  &
                                        atom_index,  &
                                        try_trigonom_tab(1,1,i_point,i_division),  &
                                        basis_l_max,  &
                                        l_ylm_max, &
                                        try_ylm_tab(1,1,i_point,i_division) )

                                   ! Now evaluate radial functions, radial kinetic energy terms, and
                                   ! possibly radial derivatives from the previously stored compressed spline arrays  
                                   call evaluate_radial_functions_p0 &
                                        ( spline_array_start,  &
                                        spline_array_end, &
                                        n_compute_atoms,  &
                                        n_compute_fns,  &
                                        try_dist_tab(1,i_point,i_division),  &
                                        try_i_r(1,i_point,i_division), &
                                        atom_index,  &
                                        i_basis_fns_inv, &
                                        basis_wave_ordered,  &
                                        try_radial_wave(1,i_point,i_division), &
                                        .false., try_n_compute_c(i_division), n_max_compute_fns_ham   &
                                        )

                                   ! tabulate total wave function value for each basis function
                                   call evaluate_waves_p0 &
                                        ( l_ylm_max,  &
                                        try_ylm_tab(1,1,i_point,i_division),  &
                                        try_dist_tab(1,i_point,i_division),  &
                                        index_lm, try_n_compute_c(i_division),  &
                                        try_i_basis(1,i_division),  &
                                        try_radial_wave(1,i_point,i_division),  &
                                        try_wave(1,i_point,i_division),  &
                                        n_compute_atoms,  &
                                        atom_index_inv,  &
                                        n_compute_fns, &
                                        i_basis_fns_inv, n_max_compute_fns_ham  &
                                        )

                                end if !(.not.first_loop)



                                ! evaluate initial potential at current integration point
                                if (force_potential.eq.1) then
                                   ! only superposition of free-atom effective potentials

                                   call evaluate_pot_superpos_p0 &
                                        ( try_i_r_full(1,i_angular),  &
                                        local_potential_parts(1,i_angular), &
                                        n_centers_integrals, centers_basis_integrals )

                                else  !  normal case

                                   ! TODO: Evaluate local kinetic_density in this subroutine

                                   call evaluate_free_atom_sums_p0 &
                                        ( try_dist_tab_full(1,i_angular), &
                                        try_i_r_full(1,i_angular),  &
                                        try_dir_tab_full_norm(1,1,i_angular), &
                                        free_hartree_superpos, &
                                        free_rho_superpos, &
                                        free_rho_gradient_superpos, &
                                        rho, &
                                        local_rho_gradient , n_centers_integrals, centers_basis_integrals &
                                        )
  
                                     call evaluate_xc &
                                         ( rho, &
                                           local_rho_gradient, &
                                           local_kinetic_density, &
                                           en_density_xc, &
                                           en_density_x, en_density_c, &
                                           local_xc_derivs, &
                                           local_xc_gradient_deriv, &
                                           local_xc_tau_deriv, &
                                           use_hartree_fock .or. use_meta_gga &
                                          )

                                   do i_spin = 1, n_spin, 1
                                      local_potential_parts &
                                           (i_spin,i_angular) =  &
                                           free_hartree_superpos &
                                           + local_xc_derivs(i_spin)


                                      local_gradient(1:3,i_spin,i_angular)  = &
                                           4*local_xc_gradient_deriv(1:3,i_spin)

                                   enddo

                                   ! Nadia: if (use_embedding_potential) then add pot_ion_embed_local
                                   if (use_embedding_potential) then

                                      call embedding_potential &
                                           ( try_coord(1,i_angular),  &
                                           pot_ion_embed_local, dummy_force  &
                                           )

                                      do i_spin = 1, n_spin, 1
                                         local_potential_parts(i_spin,i_angular) =   &
                                              local_potential_parts(i_spin,i_angular) +  &
                                              pot_ion_embed_local
                                      end do
                                   end if

                                   if (use_embedding_pp) then
                                      call ion_pseudopot_potential_v1 &
                                           ( try_coord(1,i_angular), &
                                           pp_potential,0)

                                      do i_spin = 1, n_spin, 1
                                         local_potential_parts(i_spin,i_angular) = &
                                              local_potential_parts(i_spin,i_angular) + &
                                              pp_potential
                                      end do
                                   end if


                                end if  !  end evaluation of initial potential



                                ! check whether we need scalar relativistic terms ...
                                if (flag_rel.eq.1) then


                                   call evaluate_pot_superpos_p0  &
                                        (   &
                                        try_i_r_full(1,i_angular),   &
                                        zora_potential_parts(1,i_angular),  &
                                        n_centers_integrals, centers_basis_integrals )

                                   do i_spin=1,n_spin,1

                                      zora_operator(i_spin) =  &
                                           2.d0* light_speed_sq  &
                                           / ( 2 * light_speed_sq -   &
                                           zora_potential_parts(i_spin,i_angular) )
                                   enddo

                                   t_zora = .false.
                                   do i_spin = 1, n_spin, 1
                                      t_zora = t_zora.or.  &
                                           ( abs(zora_operator(i_spin)-0.5d0)  &
                                           .gt.zora_threshold )
                                   enddo
                                   t_zora = .true.
                                end if !(flag_rel.eq.1)



                                if (use_gga.or. t_zora  ) then

                                   call evaluate_radial_functions_p0  &
                                        (   &
                                        spline_array_start,   &
                                        spline_array_end,  &
                                        n_compute_atoms,   &
                                        n_compute_fns,   &
                                        try_dist_tab(1,i_point,i_division),   &
                                        try_i_r(1,i_point,i_division),  &
                                        atom_index,   &
                                        i_basis_fns_inv,  &
                                        basis_deriv_ordered,   &
                                        radial_wave_deriv, .true.,try_n_compute_c(i_division), &
                                        n_max_compute_fns_ham &
                                        )

                                   ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                                   call tab_gradient_ylm_p0  &
                                        ( try_trigonom_tab(1,1,i_point,i_division),   &
                                        basis_l_max,   &
                                        l_ylm_max,   &
                                        n_compute_atoms,   &
                                        atom_index,  &
                                        try_ylm_tab(1,1,i_point,i_division),   &
                                        dylm_dtheta_tab(1,1),   &
                                        scaled_dylm_dphi_tab(1,1) &
                                        )

                                   !                                      write(use_unit,*) try_dir_tab(:,1,i_point,i_division)

                                   call evaluate_wave_gradient_p0 &
                                        ( try_dist_tab(1,i_point,i_division), &
                                        try_dir_tab(1,1,i_point,i_division), &
                                        try_trigonom_tab(1,1,i_point,i_division), &
                                        l_ylm_max,  &
                                        try_ylm_tab(1,1,i_point,i_division), &
                                        dylm_dtheta_tab(1,1), &
                                        scaled_dylm_dphi_tab(1,1), &
                                        index_lm, try_n_compute_c(i_division), &
                                        try_i_basis(1,i_division), &
                                        try_radial_wave(1, i_point, i_division), &
                                        radial_wave_deriv, &
                                        gradient_basis_wave, &
                                        n_compute_atoms, &
                                        atom_index_inv, &
                                        n_compute_fns, &
                                        i_basis_fns_inv, n_max_compute_fns_ham  )

                                end if

                                ! add non-relativistic contributions to the Hamiltonian matrix elements
                                call evaluate_radial_functions_p0 &
                                     ( spline_array_start,  &
                                     spline_array_end, &
                                     n_compute_atoms,  &
                                     n_compute_fns,  &
                                     try_dist_tab(1,i_point,i_division),  &
                                     try_i_r(1,i_point,i_division), &
                                     atom_index,  &
                                     i_basis_fns_inv, &
                                     basis_kinetic_ordered,  &
                                     kinetic_wave, &
                                     .false., try_n_compute_c(i_division), n_max_compute_fns_ham  &
                                     )

                                do i_spin = 1, n_spin, 1

                                   H_times_psi(1, i_point, i_spin) = 0.0d0

                                   call evaluate_H_psi_p0 &
                                        ( l_ylm_max, &
                                        try_ylm_tab(1,1,i_point,i_division), &
                                        try_dist_tab(1,i_point,i_division), &
                                        index_lm(-l_ylm_max:l_ylm_max, &
                                        0:l_ylm_max), &
                                        H_times_psi(1, i_point, i_spin), &
                                        try_radial_wave(1, i_point, i_division), &
                                        local_potential_parts(i_spin,i_angular), &
                                        try_n_compute_c(i_division), &
                                        try_i_basis(1,i_division), &
                                        n_compute_atoms,  &
                                        atom_index_inv, &
                                        n_compute_fns,  &
                                        i_basis_fns_inv,  &
                                        kinetic_wave, &
                                        zora_operator(i_spin), n_max_compute_fns_ham   &
                                        )


                                enddo ! i_spin = 1,....





                                if (t_zora) then
                                   !  ZORA scalar relativistic treatment required

                                   try_n_rel_points = try_n_rel_points + 1

                                   do i_spin = 1, n_spin, 1

                                      zora_operator(i_spin) = &
                                           light_speed_sq / &
                                           (2 * light_speed_sq - &
                                           zora_potential_parts(i_spin,i_angular))**2

                                      call  add_zora_gradient_part_p0(  &
                                           local_gradient(1:3,i_spin,i_angular), &
                                           try_i_r_full(1,i_angular), &
                                           try_dir_tab_full_norm(1,1,i_angular),  &
                                           try_dist_tab_full(1,i_angular), &
                                           zora_operator(i_spin),  n_centers_integrals, centers_basis_integrals )

                                   enddo
                                end if ! end ZORA preparations


                                if (use_gga .or.try_n_rel_points.gt. 0 ) then

                                   do i_spin = 1, n_spin, 1

                                      call add_gradient_part_to_H_p0 &
                                           ( try_n_compute_c(i_division), &
                                           gradient_basis_wave, &
                                           local_gradient(1,i_spin,i_angular), &
                                           H_times_psi(1,i_point,i_spin) )
                                   enddo
                                end if
                             end if ! (try_partition_tab > 0)
                          enddo !  end i_angular loop over part of the present trial angular shell



                          try_n_points(i_division) = i_point

                          ! evaluate the hamiltonian shell for this batch of n_points

                          if (try_n_points(i_division) .gt. 0) then
                             do i_spin = 1, n_spin, 1

                                call evaluate_hamiltonian_shell_p1 &
                                     ( try_n_points(i_division),  &
                                     try_partition,  &
                                     try_n_compute_c(i_division), &
                                     H_times_psi(1,1,i_spin), &
                                     n_max_compute_ham, try_wave(1,1,i_division),  &
                                     matrix_shell )

                                call update_full_matrix_p0 &
                                     ( try_n_compute_c(i_division),  try_n_compute_a(i_division), &
                                     try_i_basis(1,i_division),  &
                                     matrix_shell,   &
                                     new_shell(1,i_spin) &
                                     )

                             enddo

                          end if
                       end if ! (try_n_compute .gt. 0)
                    enddo ! end loop over angular divisions


                    first_loop = .false.

                    ! now check convergence
                    ! initialize shell_converge to true ...
                    shell_converged = .true.

                    ! use the maximum deviation of a single overlap matrix term 
                    ! as the convergence criterion
                    do i_spin = 1, n_spin, 1
                       call check_shell_convergence_p0 &
                            ( prev_shell(:,i_spin), new_shell(:,i_spin),  &
                            angular_acc(species(i_atom)), shell_converged )
                    enddo

                    if (.not.shell_converged) then
                       ! store present angular grid data for possible future use
                       ! also check (below) whether shell is actually converged according to other data

                       ! store grid data

                       if(n_periodic > 0)then

                          try_local_potential_parts(:,:,1) = local_potential_parts(:,:)
                          try_zora_potential_parts(:,:,1)  = zora_potential_parts(:,:)
                          try_local_gradient(:,:,:,1)      = local_gradient(:,:,:)

                       end if

                       n_lebedev(i_radial,species(i_atom)) = i_lebedev
                       n_angular(i_radial,species(i_atom)) = try_n_angular

                       do i_angular = 1,n_angular(i_radial, &
                            species(i_atom))
                          do i_coord = 1,3,1
                             r_angular( i_coord, i_angular,  &
                                  i_radial, species(i_atom) ) =  &
                                  try_r_angular(i_coord, i_angular)
                          enddo
                          w_angular(i_angular, i_radial, & 
                               species(i_atom)) =  &
                               try_w_angular(i_angular)
                       enddo

                       n_division(i_radial, species(i_atom)) = try_n_division
                       do i_division = 1, n_division(i_radial,  &
                            species(i_atom)) + 1, 1
                          division_boundaries(i_division,  &
                               i_radial, species(i_atom)) = &
                               try_division_boundaries(i_division)
                       enddo

                       do i_angular = 1,n_angular(i_radial, &
                            species(i_atom))
                          partition_tab(i_angular) =  &
                               try_partition_tab(i_angular)
                       end do

                       ! store hamiltonian matrix terms from preceding shell
                       prev_shell = new_shell

                       n_rel_points = try_n_rel_points

                       ! try next Lebedev grid
                       i_lebedev = i_lebedev + 1

                       if (i_lebedev.le.n_grids) then
                          ! get next angular grid and grid weights
!!$                             call get_angular_grid &
!!$                                  ( i_lebedev, try_n_angular,  &
!!$                                  try_r_angular,  &
!!$                                  try_w_angular &
!!$                                  )

                          try_n_angular = n_angular_lebedev(i_lebedev)
                          try_r_angular = r_angular_lebedev(:,:,i_lebedev)
                          try_w_angular = w_angular_lebedev(:,i_lebedev)

                          if (try_n_angular.gt. &
                               angular_limit(species(i_atom)) .or. full_shell(i_radial,(i_atom)) )then

                             !  maximum allowed grid size exceeded, shell is "converged"
                             shell_converged = .true.

                          else
                             ! divide the angular grid
!!$                                call divide_angular_grid( &
!!$                                     try_n_angular, &
!!$                                     try_r_angular,  &
!!$                                     try_w_angular, &
!!$                                     try_n_division, &
!!$                                     try_division_boundaries, &
!!$                                     r_radial(i_radial,species(i_atom))      &
!!$                                     )

                             try_n_division = n_division_lebedev(i_lebedev)
                             try_division_boundaries = division_boundaries_lebedev(:,i_lebedev)

                          end if

                       else ! if we're out of Lebedev grids to try, shell is "converged" anyway

                          shell_converged = .true.

                       end if
                    end if ! (shell_converged)

                 enddo  ! end loop over Lebedev grids for Hamiltonian matrix



                 if(n_periodic > 0)then

                    try_local_potential_parts(:,:,2) = local_potential_parts(:,:)
                    try_zora_potential_parts (:,:,2) = zora_potential_parts(:,:)
                    try_local_gradient     (:,:,:,2) = local_gradient(:,:,:)

                    local_potential_parts(:,:)       = try_local_potential_parts(:,:,1)
                    zora_potential_parts(:,:)        = try_zora_potential_parts (:,:,1)
                    local_gradient(:,:,:)            = try_local_gradient     (:,:,:,1)

                 end if

                 ! after convergence, increment hamiltonian for present radial shell
                 do i_spin = 1, n_spin, 1

                    hamiltonian(:,i_spin) = hamiltonian(:,i_spin)  &
                         + new_shell(:,i_spin)
                 enddo

              else
                 ! we already had the full grid, no further grids tried
                 ! simply increment hamiltonian for present radial shell, using prev_shell

                 do i_spin = 1, n_spin, 1

                    hamiltonian(:,i_spin) = hamiltonian(:,i_spin) +  &
                         prev_shell(:,i_spin)
                 enddo

              end if ! (.not. full_shell)
           end if !(angular_acc...
        end if   ! end MPI task distribution

     end do ! end radial integration loop

     ! IF MPI is used, sync all radial data before the next atom
     ! IF MPI is not used (normal case), all the synch suboutines simply return.

     if (angular_acc(species(i_atom)).gt.0.d0) then

!!$        call sync_initialize_integrals_grids(i_atom,  &
!!$             n_radial(species(i_atom)), &
!!$             n_lebedev(1:n_radial(species(i_atom)), species(i_atom)),  &
!!$             n_angular(1:n_radial(species(i_atom)), species(i_atom)), &
!!$             r_angular(:,:,  &
!!$             1:n_radial(species(i_atom)), species(i_atom)),  &
!!$             w_angular(:, 1:n_radial(species(i_atom)), species(i_atom)), &
!!$             division_boundaries(:, 1:n_radial(species(i_atom)), &
!!$             species(i_atom)), n_division(1:n_radial(species(i_atom)), &
!!$             species(i_atom)) &
!!$             )

        call sync_initialize_integrals_grids_p0(i_atom,  &
             n_radial(species(i_atom)), &
             n_lebedev(1:n_radial(species(i_atom)), species(i_atom)) )


        do i_radial = 1, n_radial(species(i_atom)), 1

           i_lebedev = n_lebedev(i_radial, species(i_atom))
           n_angular(i_radial, species(i_atom)) = n_angular_lebedev(i_lebedev)
           r_angular(:,:,i_radial, species(i_atom)) = &
                r_angular_lebedev(:,:,i_lebedev)
           w_angular(:, i_radial, species(i_atom)) = w_angular_lebedev(:, i_lebedev)
           division_boundaries(:, i_radial, species(i_atom)) = &
                division_boundaries_lebedev(:, i_lebedev)
           n_division(i_radial, species(i_atom)) = n_division_lebedev(i_lebedev)
           
        enddo

     end if
  enddo ! end integration loop over atoms

  ! store the lebedev grid indexes for further use
  lebedev_grid_index = 0
  do i_species = 1, n_species, 1
     do i_radial = 1, n_radial(i_species), 1
        
        lebedev_grid_index(i_radial, i_species) = n_lebedev(i_radial, i_species )
        
     end do
  end do

  ! sync the hamiltonian and overlap matrices
  call sync_initialize_integrals_matrices( &
       hamiltonian, overlap_matrix)

  if (out_grids) then

     ! write the entire integration grid in the format required 
     ! for a re-readin in the file grids.dat
     if (myid.eq.0) then
        open (50,file="grids_out.dat")
        write(50,'(I5)') n_atoms
        do i_atom = 1, n_atoms, 1
           write(50,'(2X,I5,1X,I5)') i_atom, &
                n_radial(species(i_atom))
           do i_radial = 1,n_radial(species(i_atom)) 
              write(50,'(4X,I5,1X,E30.15,1X,E30.15,1X,I5)') &
                   i_radial, r_radial(i_radial,species(i_atom)),  &
                   w_radial(i_radial,species(i_atom)),  &
                   n_angular(i_radial,species(i_atom))
              do i_angular = 1, n_angular(i_radial, species(i_atom))
                 write(50,'(6X,I5,1X,3(1X,E30.15),1X,E30.15)')  &
                      i_angular,  &
                      (r_angular(i_coord, i_angular, i_radial,  &
                      species(i_atom)), i_coord = 1,3,1), &
                      w_angular(i_angular, i_radial,  &
                      species(i_atom))
              enddo
           enddo
        enddo
        close(50)
     end if

  end if


  ! write integration grid per species in the format needed for control.in
  if (myid.eq.0) then
       write(use_unit,*) 
       write(use_unit,'(2X,A)') "Output of integration grids in suitable form for copy-paste into control.in:"
       do i_species = 1, n_species, 1
         write(use_unit,*)
         write(use_unit,'(2X,A,A,A)') "Species ", species_name(i_species), ":"
         do i_radial = 2, n_radial(i_species), 1
           if (n_angular(i_radial,i_species).ne.(n_angular(i_radial-1,i_species))) then
             write(use_unit,'(6X,A,1X,F8.4,1X,I4)') &
             "division", ( r_radial(i_radial-1,i_species) + 0.001 ) * bohr, &
             n_angular(i_radial-1,i_species)
           end if
         enddo
         write(use_unit,'(6X,A,1X,I4)') &
         "outer_grid", n_angular(i_radial-1,i_species)
       enddo
       write(use_unit,*)
  end if

  ! Allocatable arrays that are tracked
  if (allocated(prev_shell))          call aims_deallocate( prev_shell,                   "prev_shell" )
  if (allocated(new_shell))           call aims_deallocate( new_shell,                     "new_shell" )
  if (allocated(gradient_basis_wave)) call aims_deallocate( gradient_basis_wave, "gradient_basis_wave" )
  if (allocated(matrix_shell))        call aims_deallocate( matrix_shell,               "matrix_shell" )
  if (allocated(wave))                call aims_deallocate( wave,                               "wave" )
  if (allocated(H_times_psi))         call aims_deallocate( H_times_psi,                 "H_times_psi" )

  ! deallocate

  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(dylm_dtheta_tab)) then
     deallocate(dylm_dtheta_tab)
  end if
  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate(scaled_dylm_dphi_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if (allocated(n_lebedev)) then
     deallocate (n_lebedev)
  end if
  if (allocated(try_r_angular)) then
     deallocate (try_r_angular)
  end if
  if (allocated(try_w_angular)) then
     deallocate (try_w_angular)
  end if
  if (allocated(try_division_boundaries)) then
     deallocate (try_division_boundaries)
  end if

  if (allocated(try_partition_tab)) then
     deallocate (try_partition_tab)
  end if
  if (allocated(try_partition)) then
     deallocate (try_partition)
  end if
  if (allocated(try_coord)) then
     deallocate(try_coord)
  end if
  if (allocated(try_dist_tab)) then
     deallocate(try_dist_tab)
  end if
  if (allocated(try_i_r)) then
     deallocate(try_i_r)
  end if
  if (allocated(try_i_r_full)) then
     deallocate(try_i_r_full)
  end if
  if (allocated(try_dir_tab)) then
     deallocate(try_dir_tab)
  end if
  if (allocated(try_dir_tab_full)) then
     deallocate(try_dir_tab_full)
  end if
  if (allocated(try_dir_tab_full_norm)) then
     deallocate(try_dir_tab_full_norm)
  end if
  if (allocated(try_trigonom_tab)) then
     deallocate(try_trigonom_tab)
  end if
  if (allocated(try_ylm_tab)) then
     deallocate(try_ylm_tab)
  end if
  if (allocated(try_i_basis)) then
     deallocate(try_i_basis)
  end if
  if (allocated(try_radial_wave)) then
     deallocate(try_radial_wave)
  end if
  if (allocated(try_wave)) then
     deallocate(try_wave)
  end if
  if (allocated(radial_wave)) then
     deallocate ( radial_wave )
  end if
  if (allocated(radial_wave_deriv)) then
     deallocate ( radial_wave_deriv )
  end if
  if (allocated(kinetic_wave)) then
     deallocate ( kinetic_wave )
  end if

  if (allocated(try_dist_tab_sq_full))then
     deallocate (try_dist_tab_sq_full)
  end if
  if (allocated(try_dist_tab_sq)) then
     deallocate(try_dist_tab_sq)
  end if
  if (allocated(try_dist_tab_full)) then
     deallocate(try_dist_tab_full)
  end if
  if (allocated(try_local_gradient)) then
     deallocate(try_local_gradient)
  end if
  if (allocated(try_local_potential_parts)) then
     deallocate(try_local_potential_parts)
  end if
  if (allocated(try_zora_potential_parts)) then
     deallocate(try_zora_potential_parts)
  end if

  if( allocated(i_r))then
     deallocate(i_r)
  end if
  if(allocated(trigonom_tab))then
     deallocate(trigonom_tab)
  end if
  if( allocated( ylm_tab))then
     deallocate( ylm_tab)
  end if
  if(allocated(dist_tab))then
     deallocate(dist_tab)
  end if
  if(allocated(dir_tab))then
     deallocate(dir_tab)
  end if



 


end subroutine initialize_integrals_p0

!----------------------------------------------------------------------
!****** 
