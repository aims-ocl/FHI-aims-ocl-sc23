!****s* FHI-aims/integrate_xc_realspace_p2
!  NAME
!   integrate_xc_realspace_p2
!  SYNOPSIS

subroutine integrate_xc_realspace_p2 &
     ( hartree_potential_std, rho_std, rho_gradient_std,  &
     kinetic_density_std, &
     partition_tab_std, basis_l_max, en_xc, en_pot_xc, &
     xc_realspace, &
     x_realspace, &
     c_realspace &
     )

!  PURPOSE
!  Integrates the exchange-correlation elements for the Hamiltonian matrix,
!  using a fixed basis set. The subroutine also calculates xc-energy.
!
!  We only import the Hartree potential across the grid, and evaluate
!  the XC potential on the fly. Hence, it is convenient to compute also
!  the XC energy and the average XC potential in this subroutine.
!
!  XR (23.08.12): Modified from "integrate_hamiltonian_matrix_p2.f90",
!        here we don't need to care about the relativistic correction.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use plus_u
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data, only: species_name
  use load_balancing
  use physics, only: deallocate_vdw_splines
  use pseudodata
  use energy_density
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: hartree_potential_std
  real*8, target, dimension(n_spin, n_full_points)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points) :: rho_gradient_std
  real*8, target, dimension(n_full_points)            :: partition_tab_std
  integer ::  basis_l_max (n_species)
  real*8  :: en_xc
  real*8  :: en_pot_xc
  ! when this routine is called, hamiltonian has either the dimension
  ! (n_hamiltonian_matrix_size, n_spin) or (n_local_matrix_size, n_spin)
  ! so we declare it here as a 1D assumed size array
  real*8  :: xc_realspace(*)
  real*8  :: x_realspace(*)
  real*8  :: c_realspace(*)
  real*8, target, dimension(n_spin, n_full_points)    :: kinetic_density_std

!  INPUTS
!  o hartree_potential_std -- Hartree potential
!  o rho_std -- electron density
!  o rho_gradient_std -- gradient of electron density.
!    These should only ever be referenced if (use_gga)
!    import dimensions from above (if not used, all dimensions=1)
!  o partition_tab_std -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!  o kinetic_density_std -- kinetic density of system
!    Should only ever be referenced if (use_meta_gga)
!
!  OUTPUT
!  o en_xc -- xc-energy energy
!  o en_pot_xc -- xc-potential energy
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











  !  local variables

!  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

  real*8 coord_current(3)

!  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
!  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)

  real*8,dimension(:,:),allocatable:: dist_tab
  real*8,dimension(:,:),allocatable:: dist_tab_sq

  real*8 i_r(n_max_compute_atoms)

!  real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
  real*8, dimension(:,:,:),allocatable:: dir_tab


  real*8 trigonom_tab(4,n_max_compute_atoms)

!  real*8,dimension(:,:,:),allocatable:: H_times_psi
  real*8,dimension(:,:,:),allocatable:: xc_times_psi
  real*8,dimension(:,:,:),allocatable:: x_times_psi
  real*8,dimension(:,:,:),allocatable:: c_times_psi
  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave


  real*8, dimension(:),    allocatable :: en_density_xc
  !real*8, dimension(:),    allocatable :: en_density_x
  !real*8, dimension(:),    allocatable :: en_density_c
  real*8, dimension(n_spin) :: en_density_x
  real*8, dimension(n_spin) :: en_density_c
  real*8, dimension(:, :), allocatable :: local_xc_derivs
  real*8, dimension(:, :), allocatable :: local_x_derivs
  real*8, dimension(:, :), allocatable :: local_c_derivs
  real*8, dimension(:,:,:),allocatable :: xc_gradient_deriv
  real*8, dimension(:,:,:),allocatable :: x_gradient_deriv
  real*8, dimension(:,:,:),allocatable :: c_gradient_deriv
  real*8, dimension(:, :), allocatable :: xc_tau_deriv
  real*8, dimension(:, :), allocatable :: x_tau_deriv
  real*8, dimension(:, :), allocatable :: c_tau_deriv

  real*8, dimension(:,:),  allocatable :: local_rho
  real*8, dimension(:,:,:),allocatable :: local_rho_gradient
  real*8, dimension(:,:), allocatable  :: local_kinetic_density


  !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:), allocatable :: xc_realspace_shell
  real*8, dimension(:), allocatable :: x_realspace_shell
  real*8, dimension(:), allocatable :: c_realspace_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)
  real*8 :: energy_partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave

  !     Following is all that is needed for the handling of ZORA scalar relativity

!  real*8, dimension(n_spin) :: zora_operator
!  logical, dimension(n_spin) :: t_zora
!  real*8, dimension(n_spin) :: zora_potential_parts

  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  ! This term contains contributions from the xc potential and the
  ! zora formalism (if applicable) which are summed up using Gauss' law:
  ! < grad(phi_i) | local_gradient_sum |grad(phi_j) >
  real*8, dimension(3,n_spin) :: sum_of_local_gradients
  real*8, dimension(3,n_spin) :: sum_of_local_gradients_x
  real*8, dimension(3,n_spin) :: sum_of_local_gradients_c

  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_c, n_compute_a
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns

  integer,dimension(:),  allocatable :: i_basis_fns
  integer,dimension(:,:),allocatable :: i_basis_fns_inv
  integer,dimension(:),  allocatable :: i_atom_fns

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute
  integer :: zero_index_point(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers
  integer :: batch_center(n_centers_integrals)

  !     for splitting of angular shells into "octants"

  integer division_low
  integer division_high

  !  counters

  integer i_basis_1
  integer i_basis_2
  integer i_atom, i_atom_2
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_spin
  character*200 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer ld_hamiltonian  ! leading dimension of hamiltonian in calling routine

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_gradient(:,:,:)
  real*8, pointer :: hartree_potential(:)
  real*8, pointer :: kinetic_density(:,:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i_off, i, j, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

  real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
  real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore

  ! begin work



  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating local/semilocal XC matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating local/semilocal matrix: batch-based integration."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  ! begin with general allocations

  allocate(dist_tab(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab                      ')

  allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dist_tab_sq                   ')

  allocate(dir_tab(3,n_centers_integrals, n_max_batch_size),stat=info)
  call check_allocation(info, 'dir_tab                       ')

  allocate(i_basis_fns(n_basis_fns*n_centers_integrals), stat=info)
  call check_allocation(info, 'i_basis_fns                   ')

  allocate(i_basis_fns_inv(n_basis_fns,n_centers), stat=info)
  call check_allocation(info, 'i_basis_fns_inv               ')

  allocate(i_atom_fns(n_basis_fns*n_centers_integrals),stat=info)
  call check_allocation(info, 'i_atom_fns                    ')

  allocate( en_density_xc(n_max_batch_size),stat=info)
  call check_allocation(info, 'en_density_xc                 ')

!  allocate( en_density_x(n_max_batch_size),stat=info)
!  call check_allocation(info, 'en_density_x                 ')

!  allocate( en_density_c(n_max_batch_size),stat=info)
!  call check_allocation(info, 'en_density_c                 ')

  allocate( local_xc_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_xc_derivs               ')

  allocate( local_x_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_x_derivs               ')

  allocate( local_c_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_c_derivs               ')

  allocate( xc_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_gradient_deriv             ')

  allocate( x_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'x_gradient_deriv             ')

  allocate( c_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'c_gradient_deriv             ')

  allocate( xc_tau_deriv(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_tau_deriv             ')

  allocate( x_tau_deriv(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'x_tau_deriv                  ')

  allocate( c_tau_deriv(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'c_tau_deriv                  ')

  allocate( local_rho(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho                     ')

  allocate( local_rho_gradient(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho_gradient            ')

  allocate( local_kinetic_density(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_kinetic_density         ')

  if ((flag_rel.eq.REL_none.or.flag_rel==REL_atomic_zora.or.flag_rel.eq.REL_own).and.(.not.(use_gga))) then
     !       no gradients needed
     l_ylm_max = l_wave_max
  else if ((flag_rel.eq.REL_zora).or.(use_gga).or.(flag_rel==REL_KOLNING_HARMON)) then
     l_ylm_max = l_wave_max
     allocate (gradient_basis_wave(n_max_compute_ham,3),STAT=info)
     call check_allocation(info, 'gradient_basis_wave           ')

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')

  end if

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')

  allocate ( xc_realspace_shell(n_max_compute_ham*n_max_compute_ham), STAT=info )
  call check_allocation(info, 'xc_realspace_shell                  ')

  allocate ( x_realspace_shell(n_max_compute_ham*n_max_compute_ham), STAT=info )
  call check_allocation(info, 'x_realspace_shell                  ')

  allocate ( c_realspace_shell(n_max_compute_ham*n_max_compute_ham), STAT=info )
  call check_allocation(info, 'c_realspace_shell                  ')

  allocate(xc_times_psi(n_max_compute_ham, n_max_batch_size, n_spin), STAT=info )
  call check_allocation(info, 'xc_times_psi                  ')

  allocate(x_times_psi(n_max_compute_ham, n_max_batch_size, n_spin), STAT=info )
  call check_allocation(info, 'x_times_psi                  ')

  allocate(c_times_psi(n_max_compute_ham, n_max_batch_size, n_spin), STAT=info )
  call check_allocation(info, 'c_times_psi                  ')


  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  allocate(wave(n_max_compute_ham, n_max_batch_size), STAT=info )
  call check_allocation(info, 'wave                          ')

  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')

  ! CC: Save local (!) v_xc(r) * rho(r) if Harris-Foulkes-like energy density shall be computed
  if (flag_harris_foulkes_energy_density) then
    if(.not.allocated(ed_xc_pot_energy_density)) allocate(ed_xc_pot_energy_density(n_full_points))
    ed_xc_pot_energy_density(:) = 0.0d0
  end if
  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)


  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab

    allocate(rho(n_spin,batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    allocate(hartree_potential(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,hartree_potential_std,hartree_potential)

    if(use_density_gradient) then
      allocate(rho_gradient(3,n_spin,batch_perm(n_bp)%n_full_points))
      call permute_point_array(n_bp,3*n_spin,rho_gradient_std,rho_gradient)
    else
      ! Even though rho_gradient_std is allocated to a dummy size in this case,
      ! the array rho_gradient is used below as a dummy argument in full size
      ! (calls to evaluate_xc).
      ! rho_gradient therefore shouldn't be a dangling or nullified pointer
      ! since this will generated errors when in bounds checking mode.
      ! So we have to allocate it here, although it isn't needed actually.
      allocate(rho_gradient(3,n_spin,batch_perm(n_bp)%n_full_points))
    endif

    if(use_meta_gga) then
      allocate(kinetic_density(n_spin,batch_perm(n_bp)%n_full_points))
      call permute_point_array(n_bp,n_spin,kinetic_density_std,kinetic_density)
    else
      ! Same as above? Need to allocate so no dangling pointers
      allocate(kinetic_density(n_spin,batch_perm(n_bp)%n_full_points))
    endif

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

    ld_hamiltonian = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    hartree_potential => hartree_potential_std
    rho_gradient => rho_gradient_std
    kinetic_density => kinetic_density_std

    ld_hamiltonian = n_hamiltonian_matrix_size

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------


  if(use_embedding_pp.and.use_nonlinear_core) then
      allocate(rho_inc_partialcore(n_spin,n_full_points))
      allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
      do i_spin = 1,n_spin
         rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
         rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
         if(use_density_gradient) then
            rho_gradient_inc_partialcore(:,i_spin,:) = &
               rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
         endif
      enddo
  endif


!  hamiltonian(1:ld_hamiltonian*n_spin) = 0.d0
  xc_realspace(1:ld_hamiltonian*n_spin) = 0.d0
  x_realspace(1:ld_hamiltonian*n_spin) = 0.d0
  c_realspace(1:ld_hamiltonian*n_spin) = 0.d0

  i_basis_fns_inv = 0

  en_xc = 0.d0
  en_pot_xc = 0.d0

  ! initialize index_lm
  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_full_points = 0
  i_full_points_2 = 0

  ! perform partitioned integration, batch by batch of integration point.
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions

  call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
  time0 = mpi_wtime()

  do i_my_batch = 1, n_my_batches_work, 1
     if(get_batch_weights) time_start = mpi_wtime()

     n_compute_c = 0
     n_compute_a = 0
     i_basis = 0

     i_point = 0

     ! loop over one batch
     do i_index = 1, batches_work(i_my_batch)%size, 1

        i_full_points_2 = i_full_points_2 + 1

        if (partition_tab(i_full_points_2).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
           coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

           if(n_periodic > 0)then
              call map_to_center_cell(coord_current(1:3) )
           end if

           ! compute atom-centered coordinates of current integration point,
           ! as viewed from all atoms
           call tab_atom_centered_coords_p0 &
                ( coord_current,  &
                dist_tab_sq(1,i_point),  &
                dir_tab(1,1,i_point), &
                n_centers_integrals, centers_basis_integrals )

           ! determine which basis functions are relevant at current integration point,
           ! and tabulate their indices

           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           if (.not.prune_basis_once) then
              call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
           endif

        end if
     enddo


     if (prune_basis_once) then
        n_compute_c = batches_work(i_my_batch)%batch_n_compute
        i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
     end if

     ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
     call collect_batch_centers_p2 &
     ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals, &
       n_batch_centers, batch_center &
     )

     n_points = i_point

     ! Perform actual integration if more than 0 basis functions
     ! are actually relevant on the present angular shell ...
     if (n_compute_c.gt.0) then

        n_rel_points = 0
        i_point = 0

        ! CC: UGLY! The call evaluate_xc_energy_shell() is not aware of i_full_point (global index)
        !           but just of the batch-wise non-zero-partition non-zero-density integration points
        !           To recover the local (!) v_xc(r) * rho(r) term entering the Harris-Foulkes-like 
        !           energy density, we construct a i_point <-> i_full_point map  
        if (flag_harris_foulkes_energy_density) then
          if(.not.allocated(ed_i_point_to_i_full_point_map)) allocate(ed_i_point_to_i_full_point_map(batches(i_my_batch)%size))
          ed_i_point_to_i_full_point_map(:) = - 1
        end if

        ! loop over one batch of integration points
        do i_index = 1, batches_work(i_my_batch)%size, 1

           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)!SAG
              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if


!              if (flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON) then
!
!                 call tab_global_geometry_p0 &
!                      ( dist_tab_sq(1,i_point), &
!                      dir_tab(1,1,i_point), &
!                      dist_tab_full, &
!                      i_r_full, &
!                      dir_tab_full_norm, &
!                      n_centers_integrals,  centers_basis_integrals)
!
!              end if

              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)
              energy_partition(i_point) = partition_tab(i_full_points)

              ! CC: i_point <-> i_full_point map : See comment above!
              if (flag_harris_foulkes_energy_density) then
                ed_i_point_to_i_full_point_map(i_point) = i_full_points
              end if
              
              ! for vectorized xc
              do i_spin = 1, n_spin, 1
                 local_rho(i_spin,i_point) = rho(i_spin,i_full_points)
              enddo
              if (use_gga) then
                 do i_spin = 1, n_spin, 1
                    do i_coord = 1,3,1
                       local_rho_gradient(i_coord,i_spin,i_point) = &
                            rho_gradient(i_coord,i_spin,i_full_points)
                    enddo
                 enddo
              end if
              if (use_meta_gga) then
                 do i_spin = 1, n_spin, 1
                    local_kinetic_density(i_spin,i_point) = &
                          kinetic_density(i_spin,i_full_points)
                 enddo
              endif

              n_compute_atoms = 0
              n_compute_fns = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
              ! without any copying and without doing any unnecessary operations.
              ! The price is that the interface is no longer explicit in terms of physical
              ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

              call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
                     n_compute_atoms, atom_index, atom_index_inv, &
                     n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                     i_atom_fns, spline_array_start, spline_array_end, &
                     n_centers_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )


              ! Tabulate distances, unit vectors, and inverse logarithmic grid units
              ! for all atoms which are actually relevant
              call tab_local_geometry_p2 &
                   ( n_compute_atoms, atom_index, &
                     dist_tab(1,i_point), i_r )

              ! compute trigonometric functions of spherical coordinate angles
              ! of current integration point, viewed from all atoms
              call tab_trigonom_p0 &
                   ( n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab )

!              if ((use_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON) ) then
              if (use_gga) then
                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )

              else
                ! tabulate distance and Ylm's w.r.t. other atoms
                call tab_wave_ylm_p0 &
                   ( n_compute_atoms, atom_index,  &
                   trigonom_tab, basis_l_max,  &
                   l_ylm_max, ylm_tab )
              end if

              ! Now evaluate radial functions
              ! from the previously stored compressed spline arrays
              call evaluate_radial_functions_p0  &
                   (   spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_wave_ordered, radial_wave,  &
                   .false. , n_compute_c, n_max_compute_fns_ham )


              ! tabulate total wave function value for each basis function
              call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, &
                     l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave, wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )


              ! in the remaining part of the subroutine, some decisions (scalar
              !  relativity) depend on the potential; must therefore evaluate the
              ! potential and derived quantities right here

              ! Local exchange-correlation parts of the potential are evaluated
              ! right here, to avoid having to store them separately elsewhere.
              ! For large systems, savings are significant
              if ((use_hartree_fock.or.use_meta_gga).and.&
              ! if ((use_hartree_fock).and.&
                  (first_integration)) then
                 call evaluate_xc &
                      ( rho(1,i_full_points), &
                      rho_gradient(1,1,i_full_points), &
                      kinetic_density(1,i_full_points), &
                      en_density_xc(i_point), &
                      en_density_x, en_density_c, &
                      local_xc_derivs(1,i_point), &
                      xc_gradient_deriv(1,1,i_point), &
                      xc_tau_deriv(1,i_point), &
                      .true., coord_current &
                               )  !SAG
              else

                 if(use_embedding_pp.and.use_nonlinear_core) then

                    call evaluate_xc  &
                         ( rho_inc_partialcore(1,i_full_points),   &
                         rho_gradient_inc_partialcore(1,1,i_full_points),  &
                         kinetic_density(1,i_full_points), &
                         en_density_xc(i_point), &
                         en_density_x, en_density_c, &
                         local_xc_derivs(1,i_point),  &
                         xc_gradient_deriv(1,1,i_point), &
                         xc_tau_deriv(1,i_point), &
                         .false., coord_current  &
                         ) 

                 else

                    call evaluate_xc  &
                         ( rho(1,i_full_points),   &
                         rho_gradient(1,1,i_full_points),  &
                         kinetic_density(1,i_full_points), &
                         en_density_xc(i_point), &
                         en_density_x, en_density_c, &
                         local_xc_derivs(1,i_point),  &
                         xc_gradient_deriv(1,1,i_point), &
                         xc_tau_deriv(1,i_point), &
                         .false., coord_current  &
                         ) 
!---------------W.Chibani---------------------------------------------------------
  
                   call evaluate_xc_split &
                        ( rho(1,i_full_points),&
                        rho_gradient(1,1,i_full_points), &
                        kinetic_density(1,i_full_points), &
                        en_density_xc(i_point), &
                        en_density_x, en_density_c, &
                        local_x_derivs(1,i_point),&
                        local_c_derivs(1,i_point), &
                        x_gradient_deriv(1,1,i_point),&
                        c_gradient_deriv(1,1,i_point), &
                        x_tau_deriv(1,i_point), &
                        c_tau_deriv(1,i_point), &
                        coord_current &
                        )

!--------------------------------------------------------------------------------

                 endif

              end if

              do i_spin = 1, n_spin, 1
!                 local_potential_parts(i_spin) =   &
!                      hartree_potential(i_full_points) +   &
!                      local_xc_derivs(i_spin,i_point)

                 if (use_gga) then
                    sum_of_local_gradients(1:3,i_spin) =   &
                         xc_gradient_deriv(1:3,i_spin,i_point)*4.d0
                    sum_of_local_gradients_x(1:3,i_spin) =   &
                         x_gradient_deriv(1:3,i_spin,i_point)*4.d0
                    sum_of_local_gradients_c(1:3,i_spin) =   &
                         c_gradient_deriv(1:3,i_spin,i_point)*4.d0

                 else
                    sum_of_local_gradients(1:3,i_spin) = 0.d0
                 end if


              enddo

              if (use_gga) then
                 ! we require the gradient of each basis function

                 ! tabulate radial derivatives of those radial functions
                 ! which are actually non-zero at current point, using vectorized splines
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end,  &
                      n_compute_atoms, n_compute_fns,   &
                      dist_tab(1,i_point), i_r,  &
                      atom_index, i_basis_fns_inv,  &
                      basis_deriv_ordered,   &
                      radial_wave_deriv(1), .true.,  &
                      n_compute_c, n_max_compute_fns_ham )

                 ! and finally, assemble the actual gradients


                 call evaluate_wave_gradient_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                   l_ylm_max, ylm_tab,  &
                   dylm_dtheta_tab,  &
                   scaled_dylm_dphi_tab,  &
                   radial_wave,  &
                   radial_wave_deriv,  &
                   gradient_basis_wave,  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )


              end if


              ! Now, evaluate vector of components H*phi(i,r)
              ! Local potential parts first; in the case of GGA,
              ! the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate
              ! H*phi(i,r) all in one

              ! First, obtain radial kinetic energy terms from vectorized splines
              call evaluate_radial_functions_p0  &
                   ( spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_kinetic_ordered, kinetic_wave(1),  &
                   .false., n_compute_c, n_max_compute_fns_ham )


              do i_spin = 1, n_spin, 1
                 
                 call evaluate_xc_psi_p2 &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab, &
                   radial_wave, xc_times_psi(1:n_compute_c,i_point,i_spin), &
                   local_xc_derivs (i_spin,i_point),  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point )

!---------------W.Chibani---------------------------------------------------------

                 call evaluate_xc_psi_p2 &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab, &
                   radial_wave, x_times_psi(1:n_compute_c,i_point,i_spin), &
                   local_x_derivs (i_spin,i_point),  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point )

                 call evaluate_xc_psi_p2 &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab, &
                   radial_wave, c_times_psi(1:n_compute_c,i_point,i_spin), &
                   local_c_derivs (i_spin,i_point),  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point )

!----------------------------------------------------------------------------------

                 if(use_gga) then
                   call add_gradient_part_to_H_p0 &
                   ( n_compute_c, gradient_basis_wave, &
                     sum_of_local_gradients(1,i_spin),  &
                     xc_times_psi(1, i_point, i_spin) )


!---------------W.Chibani---------------------------------------------------------
                   call add_gradient_part_to_H_p0 &
                   ( n_compute_c, gradient_basis_wave, &
                     sum_of_local_gradients_x(1,i_spin),  &
                     x_times_psi(1, i_point, i_spin) )

                   call add_gradient_part_to_H_p0 &
                   ( n_compute_c, gradient_basis_wave, &
                     sum_of_local_gradients_c(1,i_spin),  &
                     c_times_psi(1, i_point, i_spin) )
!----------------------------------------------------------------------------------

                 endif
                   
              enddo

              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0


           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch

        ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
        ! work separately for each spin channel
        do i_spin = 1, n_spin, 1

!XR: this evaluates the integration of the xc matrix for a given batch
           call evaluate_hamiltonian_shell_p1  &
                ( n_points, partition(1), n_compute_c, &
                xc_times_psi(1,1,i_spin),  &
                n_max_compute_ham, wave(1,1),  &
                xc_realspace_shell )

!---------------W.Chibani---------------------------------------------------------
           call evaluate_hamiltonian_shell_p1  &
                ( n_points, partition(1), n_compute_c, &
                x_times_psi(1,1,i_spin),  &
                n_max_compute_ham, wave(1,1),  &
                x_realspace_shell )

           call evaluate_hamiltonian_shell_p1  &
                ( n_points, partition(1), n_compute_c, &
                c_times_psi(1,1,i_spin),  &
                n_max_compute_ham, wave(1,1),  &
                c_realspace_shell )
!---------------------------------------------------------------------------------


           if(use_batch_permutation > 0) then
              ! If use_batch_permutation > 0 is set, the local hamiltonian is always stored
              ! in full form for the local basis functions

              ! Get position of basis functions of current batch within local hamiltonian
              do i=1,n_compute_c
                 ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
              enddo

              ! Insert hamiltonian_shell of current batch
              do i=1,n_compute_c
                 i_off = (i_spin-1)*ld_hamiltonian + (ins_idx(i)*(ins_idx(i)-1))/2
                 do j=1,i ! n_compute_c
!                    hamiltonian(ins_idx(j)+i_off) = hamiltonian(ins_idx(j)+i_off) &
!                                                  + hamiltonian_shell(j+(i-1)*n_compute_c)

                    xc_realspace(ins_idx(j)+i_off) = xc_realspace(ins_idx(j)+i_off) &
                                                  + xc_realspace_shell(j+(i-1)*n_compute_c)
!---------------W.Chibani---------------------------------------------------------
                    x_realspace(ins_idx(j)+i_off) = x_realspace(ins_idx(j)+i_off) &
                                                  + x_realspace_shell(j+(i-1)*n_compute_c)

                    c_realspace(ins_idx(j)+i_off) = c_realspace(ins_idx(j)+i_off) &
                                                  + c_realspace_shell(j+(i-1)*n_compute_c)
!---------------------------------------------------------------------------------

                 enddo
              enddo
           else

              call update_full_matrix_p0X(  &
                   n_compute_c, n_compute_c, i_basis(1), xc_realspace_shell,    &
                   xc_realspace((i_spin-1)*ld_hamiltonian+1) )

!---------------W.Chibani---------------------------------------------------------
              call update_full_matrix_p0X(  &
                   n_compute_c, n_compute_c, i_basis(1), x_realspace_shell,    &
                   x_realspace((i_spin-1)*ld_hamiltonian+1) )

              call update_full_matrix_p0X(  &
                   n_compute_c, n_compute_c, i_basis(1), c_realspace_shell,    &
                   c_realspace((i_spin-1)*ld_hamiltonian+1) )
!---------------------------------------------------------------------------------
           endif

        enddo
        call evaluate_xc_energy_shell  &
             ( n_points, energy_partition, en_density_xc, local_xc_derivs,  &
             xc_gradient_deriv, xc_tau_deriv, local_rho, local_rho_gradient,  &
             local_kinetic_density, en_xc, en_pot_xc )

        ! CC: Everything has been mapped in evaluate_xc_energy_shell above
        !     i_point <-> i_full_point map (see comment above) not required anymore
        if (flag_harris_foulkes_energy_density) then
          if(allocated(ed_i_point_to_i_full_point_map)) deallocate(ed_i_point_to_i_full_point_map)
        end if
        
     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches


  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) &
    info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)
  call sync_update_xc_potential( en_xc, en_pot_xc )

  !     synchronise the hamiltonian
!  if(.not. use_local_index) call sync_integrate_hamiltonian( hamiltonian )
  if(.not. use_local_index) call sync_integrate_hamiltonian( xc_realspace )

!---------------W.Chibani---------------------------------------------------------
  if(.not. use_local_index) call sync_integrate_hamiltonian( x_realspace )
  if(.not. use_local_index) call sync_integrate_hamiltonian( c_realspace )
!---------------------------------------------------------------------------------


  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)

  if(allocated( i_r_full             )) deallocate( i_r_full             )
  if(allocated( dir_tab_full_norm    )) deallocate( dir_tab_full_norm    )
  if(allocated( dist_tab_full        )) deallocate( dist_tab_full        )
  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( wave                 )) deallocate( wave                 )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( xc_times_psi         )) deallocate( xc_times_psi         )
  if(allocated( x_times_psi         )) deallocate( x_times_psi         )
  if(allocated( c_times_psi         )) deallocate( c_times_psi         )
  if(allocated( xc_realspace_shell     )) deallocate( xc_realspace_shell     )
  if(allocated( x_realspace_shell     )) deallocate( x_realspace_shell     )
  if(allocated( c_realspace_shell     )) deallocate( c_realspace_shell     )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if(allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if(allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
  if(allocated( local_kinetic_density )) deallocate( local_kinetic_density )
  if(allocated( local_rho_gradient   )) deallocate( local_rho_gradient   )
  if(allocated( local_rho            )) deallocate( local_rho            )
  if(allocated( xc_tau_deriv         )) deallocate( xc_tau_deriv         )
  if(allocated( x_tau_deriv          )) deallocate( x_tau_deriv          )
  if(allocated( c_tau_deriv          )) deallocate( c_tau_deriv          )
  if(allocated( xc_gradient_deriv    )) deallocate( xc_gradient_deriv    )
  if(allocated( x_gradient_deriv    )) deallocate( x_gradient_deriv    )
  if(allocated( c_gradient_deriv    )) deallocate( c_gradient_deriv    )
  if(allocated( local_xc_derivs      )) deallocate( local_xc_derivs      )
  if(allocated( local_x_derivs      )) deallocate( local_x_derivs      )
  if(allocated( local_c_derivs      )) deallocate( local_c_derivs      )
  if(allocated( en_density_xc        )) deallocate( en_density_xc        )
!  if(allocated( en_density_x        )) deallocate( en_density_x        )
!  if(allocated( en_density_c        )) deallocate( en_density_c        )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )

  if(allocated( rho_inc_partialcore  )) then 
     deallocate( rho_inc_partialcore  )
  endif

  if(allocated( rho_gradient_inc_partialcore )) then 
     deallocate( rho_gradient_inc_partialcore )
  endif

  
  if(use_batch_permutation > 0) then
    deallocate(rho)
    deallocate(hartree_potential)
    deallocate(rho_gradient) ! always allocated
    deallocate(kinetic_density) ! always allocated
    deallocate(ins_idx)
  endif

end subroutine integrate_xc_realspace_p2
!******
