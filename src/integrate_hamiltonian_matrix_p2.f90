!****** FHI-aims/integrate_real_hamiltonian_matrix_p2
!  NAME
!   integrate_real_hamiltonian_matrix_p2
!  SYNOPSIS

subroutine integrate_real_hamiltonian_matrix_p2 &
     ( hartree_potential_std, rho_std, rho_gradient_std, kinetic_density_std, &
     partition_tab_std, basis_l_max, en_xc, en_pot_xc, hamiltonian, &
     en_vdw, en_pot_vdw )

!  PURPOSE
!
!  Integrates the matrix elements for the Hamiltonian matrix,
!  using a fixed basis set. The subroutine also calculates the xc-energy.
!
!  We only import the Hartree potential across the grid, and evaluate
!  the XC potential on the fly. Hence, it is convenient to compute also
!  the XC energy and the average XC potential in this subroutine.
!
!  If you wish to understand this routine, DO NOT assume that it is optional
!  to read the papers that describe the integration philosophy and real-space grid handling
!  in FHI-aims. There are two such papers and, yes, you really do need to read them both.
!  If you do not, do not be surprised why you do not understand that there are certain
!  odd index arrays and loops that you did not expect. The following papers do not explain
!  all the details - but they are a necessary start.
!
!  [1] Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,"
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler"
!    'Ab initio molecular simulations with numeric atom-centered orbitals'"
!    Computer Physics Communications 180, 2175-2196 (2009)"
!    http://dx.doi.org/10.1016/j.cpc.2009.06.022"
!
!  [2] Ville Havu, Volker Blum, Paula Havu, and Matthias Scheffler,"
!    'Efficient O(N) integration for all-electron electronic structure calculation'"
!    'using numerically tabulated basis functions'"
!    Journal of Computational Physics 228, 8367-8379 (2009)."
!    http://dx.doi.org/10.1016/j.jcp.2009.08.008"
!
!  This routine and its subroutines perform an essential task in FHI-aims - an
!  efficient, domain-decomposed, parallel integration of the central matrix elements
!  of DFT in 3D space. The routine is complex for a reason. Some operations could be
!  written much simpler, but would cost far more time as a result. This explains
!  particularly the use of index arrays, which map the non-zero basis functions
!  and other items back onto a particular "group" of points. You may want to understand
!  particularly the reasons for the use of packed matrices and the choice of the
!  basis functions that are needed in periodic boundary conditions.
!
!  For the
!  periodic case, it is imperative to understand EXACTLY the meaning and index ranges
!  of Eqs. (22)-(25) of Ref. [1] above. This routine handles ONLY Eq. (24), the real-
!  space integral over all localized basis functions that touch the zeroth unit cell.
!
!  The number of basis functions that touch the zeroth unit cell is given by
!  n_centers_basis_I . This could be a large number, and has already been evaluated in
!  pbc_lists.f90 . The identities of these basis functions are stored in individual
!  arrays. For example, the array
!
!     Cbasis_to_center(i_basis_1)
!
!  is used directly in subroutine prune_basis_p2 . For each of the
!  n_centers_basis_I basis functions, this array remembers the number of
!  the atomic center at which that basis function is located. If you want to
!  understand how the properties of individual basis functions are stored, you
!  thus have to go to the definitions of these arrays elsewhere in the code
!  and see how they are set up.
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
  use synchronize_mpi_basic, only: sync_vector
  use localorb_io
  use constants
  use species_data, only: species_name, species_z
  use load_balancing
  use physics, only: deallocate_vdw_splines, rho_pce, rho_pce_gradient
  use pseudodata
  use energy_density
  use timing
  use rel_x2c_mod
  use vdw_correction
  use mbd_std_wrapper, only: &
      mbd_std_potential, mbd_std_cleanup, mbd_self_consistent
  use lpb_solver_utilities, only: evaluate_sc_Gnonmf, evaluate_Gnonmf, &
    evaluate_Gnonmf_energy_shell, en_pot_Gnonmf, mpb_solver_started, &
    en_Gnonmf, surface_and_volume_calc, surface_mpb, volume_mpb,&
    alphaplusgamma_mpb, beta_mpb, evaluate_cavity_surface_and_volume
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use pbc_lists
  use aims_gpu, only: gpu_hamiltonian_used
  use c_helper, only: fort_log_to_c_int

  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: hartree_potential_std
  real*8, target, dimension(n_spin, n_full_points)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points) :: rho_gradient_std
  real*8, target, dimension(n_full_points)            :: partition_tab_std
  integer ::  basis_l_max (n_species)
  real*8  :: en_xc
  real*8  :: en_pot_xc
  real*8  :: en_vdw
  real*8  :: en_pot_vdw
  real*8, target, dimension(n_spin, n_full_points)    :: kinetic_density_std

  ! when this routine is called, hamiltonian has either the dimension
  ! (n_hamiltonian_matrix_size, n_spin) or (n_local_matrix_size, n_spin)
  ! so we declare it here as a 1D assumed size array
  real*8  :: hamiltonian(*)

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
!  o hamiltonian -- Hamiltonian matrix
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

  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab, ylm_tab_s

  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab, dylm_dtheta_tab_s
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab, scaled_dylm_dphi_tab_s

  real*8 coord_current(3)

!  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
!  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)

  integer ::  basis_l_max_s (n_species)
  real*8,dimension(:,:),allocatable:: dist_tab, dist_tab_s
  real*8,dimension(:,:),allocatable:: dist_tab_sq, dist_tab_sq_s

  real*8 i_r(n_max_compute_atoms), i_r_s(n_max_compute_atoms)

!  real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
  real*8, dimension(:,:,:),allocatable:: dir_tab, dir_tab_s


  real*8 trigonom_tab(4,n_max_compute_atoms), trigonom_tab_s(4,n_max_compute_atoms)

  real*8,dimension(:,:,:),allocatable:: H_times_psi
  real*8,dimension(:,:,:),allocatable:: V_times_psi ! for X2C, it saves V_times_psi, but for Q4C, it saves (V+T)_times_psi
  real*8,dimension(:,:,:),allocatable:: T_times_psi ! for X2C, it saves T_times_psi, but for Q4C, it is empty
  real*8,dimension(:,:,:),allocatable:: W_times_psi ! for X2C, it saves V_times_psi_small, but for Q4C, it saves (V+T)_times_psi_small
  real*8,dimension(:,:,:),allocatable:: one_times_psi ! for generating V matrix (for X2C and 4C-DKS)
  real*8,dimension(:,:,:),allocatable:: sq_psi

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_s ! small component
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv_s ! small comp.
  real*8,dimension(:)  ,allocatable:: kinetic_wave   ! radial kinetic wave
  real*8,dimension(:)  ,allocatable:: kinetic_wave_s ! small comp. radial kinetic wave
  real*8,dimension(:,:),allocatable:: wave
  real*8,dimension(:,:),allocatable:: wave_s
  real*8,dimension(:,:),allocatable:: wave_ts ! small comp. kinetic wave

  real*8, dimension(:),    allocatable :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8, dimension(:, :), allocatable :: local_xc_derivs
  real*8, dimension(:,:,:),allocatable :: xc_gradient_deriv
  real*8, dimension(:, :), allocatable :: xc_tau_deriv

  real*8, dimension(:,:),  allocatable :: local_rho
  real*8, dimension(:,:,:),allocatable :: local_rho_gradient
  real*8, dimension(:,:),  allocatable :: local_kinetic_density

  real*8, dimension(:),  allocatable :: vdw_potential


  !     Auxiliary Hamiltonian matrix, to sum up contributions from only a
  !     single batch of integration points.
  real*8, dimension(:), allocatable :: hamiltonian_shell
  real*8, dimension(:), allocatable :: dirac_v_shell ! see the comments for V_times_psi
  real*8, dimension(:), allocatable :: dirac_t_shell
  real*8, dimension(:), allocatable :: dirac_w_shell
  real*8, dimension(:), allocatable :: dirac_ss_shell
  integer :: i_compute

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)
  real*8 :: energy_partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave, gradient_basis_wave_s

  !     necessary components for meta-GGA SCF. AJL
  real*8, dimension(:,:), allocatable :: gradient_basis_wave_store
  real*8, dimension(:,:,:), allocatable :: left_side_of_mgga_dot_product

  !     Following is all that is needed for the handling of ZORA scalar relativity

  real*8, dimension(n_spin) :: zora_operator
  logical, dimension(n_spin) :: t_zora
  real*8, dimension(n_spin) :: zora_potential_parts

  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  real*8, dimension(:,:,:,:), allocatable :: zora_vector1
  real*8, dimension(:,:,:,:), allocatable :: zora_vector2

  ! This term contains contributions from the xc potential and the
  ! zora formalism (if applicable) which are summed up using Gauss' law:
  ! < grad(phi_i) | local_gradient_sum |grad(phi_j) >
  real*8, dimension(3,n_spin) :: sum_of_local_gradients

  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_c, n_compute_small
!  integer :: i_basis(n_centers_basis_I)
  integer,dimension(:),allocatable :: i_basis, i_basis_small

  integer :: n_compute_fns, n_compute_fns_s

!  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
!  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
!  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer,dimension(:),  allocatable :: i_basis_fns, i_basis_fns_s
  integer,dimension(:,:),allocatable :: i_basis_fns_inv, i_basis_fns_inv_s
  integer,dimension(:),  allocatable :: i_atom_fns, i_atom_fns_s

  integer :: n_compute_atoms, n_compute_atoms_s
  integer :: atom_index(n_centers_integrals), atom_index_s(n_centers_integrals)
  integer :: atom_index_inv(n_centers), atom_index_inv_s(n_centers)

  integer :: spline_array_start(n_centers_integrals), spline_array_start_s(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals), spline_array_end_s(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_atoms), one_over_dist_tab_s(n_max_compute_atoms)

  ! indices for basis functions that are nonzero at current point

  integer :: rad_index(n_max_compute_atoms), rad_index_small(n_max_compute_atoms)
  integer :: wave_index(n_max_compute_fns_ham), wave_index_small(n_max_compute_fns_ham)
  integer :: l_index(n_max_compute_fns_ham), l_index_small(n_max_compute_fns_ham)
  integer :: l_count(n_max_compute_fns_ham), l_count_small(n_max_compute_fns_ham)
  integer :: fn_atom(n_max_compute_fns_ham), fn_atom_small(n_max_compute_fns_ham)

  ! indices for known zero basis functions at current point
  integer :: n_zero_compute, n_zero_compute_s
  integer :: zero_index_point(n_max_compute_ham), zero_index_point_s(n_max_compute_ham)

  ! active atoms in current batch
  integer :: n_batch_centers, n_batch_centers_s
  integer :: batch_center(n_centers_integrals), batch_center_s(n_centers_integrals)

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

  ! The actual size of full point arrays
  ! Use this for local variable allocation once it is set AFTER the load balancing part
  integer :: my_n_full_points

  integer :: i_spin
  character*200 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer ld_hamiltonian  ! leading dimension of hamiltonian in calling routine
  integer ld_dirac_large, ld_dirac_small

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

  ! Timers
  character(*), parameter :: deffmt = "10X"

  real*8       :: time_total_local
  real*8       :: clock_total_local

  logical :: gpu_save, index_on_gpu

  real*8, dimension(:), allocatable :: local_Gnonmf_derivs
  real*8, dimension(:,:), allocatable :: Gnonmf_gradient_deriv
  logical :: Gnonmf_calc
  real*8 :: dummy

  !BL: map for hamiltonian matrix elements for cuda
  integer, dimension(:), allocatable :: map

  real*8 :: c2

  ! begin work
  if (solvent_method.eq.SOLVENT_MPB.and.evaluate_sc_Gnonmf.and.mpb_solver_started.and.surface_and_volume_calc) then
    Gnonmf_calc = .True.
  else
    Gnonmf_calc = .False.
  end if

  gpu_save = use_gpu

  if (use_gpu_hamiltonian .and. .not. use_gpu) then
    ! A failsafe.  Only works for myid=0, though.
    write(info_str,'(2X,A)')
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "You have request GPU acceleration of Hamiltonian matrix integration, but no GPU"
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "is available.  Turning off GPU acceleration for Hamiltonian matrix integration."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    use_gpu = .false.
  end if

  use_gpu  = use_gpu .and. use_gpu_hamiltonian

  call get_timestamps (time_total_local, clock_total_local)

  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating Hamiltonian matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating Hamiltonian matrix: batch-based integration."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  if (use_gpu) then
    write(info_str,'(2X,A)') "GPU acceleration will be used when integrating the Hamiltonian matrix."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    gpu_hamiltonian_used = .true.
  else
    gpu_hamiltonian_always_used = .false.
    gpu_hamiltonian_used = .false.
  end if

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

  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! Fully-relativistic cases:
    allocate(dist_tab_s(n_centers_integrals, n_max_batch_size),stat=info)
    call check_allocation(info, 'dist_tab_s                    ')
    allocate(dist_tab_sq_s(n_centers_integrals, n_max_batch_size),stat=info)
    call check_allocation(info, 'dist_tab_sq_s                 ')
    allocate(dir_tab_s(3,n_centers_integrals, n_max_batch_size),stat=info)
    call check_allocation(info, 'dir_tab_s                     ')
    allocate(i_basis_fns_s(n_basis_fns*n_centers_integrals), stat=info)
    call check_allocation(info, 'i_basis_fns_s                 ')
    allocate(i_basis_fns_inv_s(n_basis_fns,n_centers), stat=info)
    call check_allocation(info, 'i_basis_fns_inv_s             ')
    allocate(i_atom_fns_s(n_basis_fns*n_centers_integrals),stat=info)
    call check_allocation(info, 'i_atom_fns_s                  ')
  endif

  allocate( en_density_xc(n_max_batch_size),stat=info)
  call check_allocation(info, 'en_density_xc                 ')

  allocate( local_xc_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_xc_derivs               ')

  if (Gnonmf_calc) then
    allocate( local_Gnonmf_derivs(n_max_batch_size),stat=info)
    call check_allocation(info, 'local_Gnonmf_derivs               ')
    local_Gnonmf_derivs=0d0
    allocate( Gnonmf_gradient_deriv(3,n_max_batch_size),stat=info)
    call check_allocation(info, 'Gnonmf_gradient_deriv             ')
    Gnonmf_gradient_deriv=0d0
  end if

  allocate( xc_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_gradient_deriv             ')

  allocate( xc_tau_deriv(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_tau_deriv             ')

  allocate( local_rho(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho                     ')

  allocate( local_rho_gradient(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho_gradient            ')

  allocate( local_kinetic_density(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_kinetic_density         ')

  if ((flag_rel.eq.REL_none.or.flag_rel==REL_atomic_zora.or.flag_rel.eq.REL_own).and.(.not.(use_gga)).and.(.not.Gnonmf_calc)) then
     ! no gradients needed
     l_ylm_max = l_wave_max
  else if (flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then ! Fully relativistic methods
     c2 = 1.d0/(2*light_speed)
     basis_l_max_s(:) = basis_l_max(:)+1

     l_ylm_max = l_wave_max+1
     allocate( ylm_tab_s( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
     call check_allocation(info, 'ylm_tab_s                     ')

     if(use_gga)then
        call aims_allocate( gradient_basis_wave,   n_max_compute_ham,3, "gradient_basis_wave" )
        call aims_allocate( gradient_basis_wave_s, n_max_compute_ham,3, "gradient_basis_wave_s" )
        allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
        call check_allocation(info, 'dylm_dtheta_tab               ')
        allocate( dylm_dtheta_tab_s( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
        call check_allocation(info, 'dylm_dtheta_tab_s             ')
        allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
        call check_allocation(info, 'scaled_dylm_dphi_tab          ')
        allocate( scaled_dylm_dphi_tab_s( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
        call check_allocation(info, 'scaled_dylm_dphi_tab_s        ')
     endif
  else if ((use_gga) .or. (use_meta_gga) .or. (Gnonmf_calc) .or. (flag_rel.eq.REL_zora) .or. (flag_rel==REL_KOLNING_HARMON)) then
     l_ylm_max = l_wave_max
     call aims_allocate( gradient_basis_wave, n_max_compute_ham,3,            "gradient_basis_wave" )
     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')
     if (use_meta_gga) then
        call aims_allocate ( gradient_basis_wave_store, n_max_compute_ham,3*n_max_batch_size,                  "gradient_basis_wave_store" )
        call aims_allocate ( left_side_of_mgga_dot_product, n_max_compute_ham, 3*n_max_batch_size, n_spin, "left_side_of_mgga_dot_product" )
     endif
  endif

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')

  call aims_allocate( hamiltonian_shell, n_max_compute_ham*n_max_compute_ham,   "hamiltonian_shell" )
  hamiltonian_shell = 0.d0
  call aims_allocate( H_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,       "H_times_psi" )
  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! Fully-relativistic cases:
    ld_dirac_large = (n_centers_basis_I+1)*n_centers_basis_I/2
    ld_dirac_small = (n_centers_basis_I_small+1)*n_centers_basis_I_small/2
    call aims_allocate( dirac_v_shell, n_max_compute_ham*n_max_compute_ham,           "dirac_v_shell" )
    call aims_allocate( dirac_w_shell, n_max_compute_ham*n_max_compute_ham,           "dirac_w_shell" )
    call aims_allocate( dirac_ss_shell, n_max_compute_ham*n_max_compute_ham,         "dirac_ss_shell" )
    call aims_allocate( V_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,       "V_times_psi" )
    call aims_allocate( W_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,       "W_times_psi" ) !(Rundong) n_max_compute_ham should be corrected for small comp. basis
    call aims_allocate( one_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,   "one_times_psi" )
    if(flag_rel.eq.REL_x2c)then
      call aims_allocate( dirac_t_shell, n_max_compute_ham*n_max_compute_ham,       "dirac_t_shell" )
      call aims_allocate( T_times_psi, n_max_compute_ham, n_max_batch_size, n_spin, "T_times_psi" )
    endif
    dirac_v_sum = 0.d0
    if(scf_iteration.eq.0.or.upw)then
       call aims_allocate( dirac_w_sum, ld_dirac_small, "dirac_w_sum" )
       dirac_w_sum = 0.d0
    endif
    if(scf_iteration.eq.0)then
       call aims_allocate( dirac_ss_sum, ld_dirac_small, "dirac_ss_sum" )
       dirac_ss_sum = 0.d0 !; dirac_ssc2_sum=0.d0
    endif
    if(scf_iteration.eq.0.and.flag_rel.eq.REL_x2c) dirac_t_sum = 0.d0
  endif

  call aims_allocate( sq_psi, n_max_compute_ham, n_max_batch_size, n_spin,                 "sq_psi" )

  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_s(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_s                 ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  call aims_allocate( wave,    n_max_compute_ham, n_max_batch_size,      "wave" )

  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')

  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
     call aims_allocate( wave_s,  n_max_compute_ham, n_max_batch_size,  "wave_s" )
     call aims_allocate( wave_ts,  n_max_compute_ham, n_max_batch_size,  "wave_ts" )
     allocate(i_basis_small(n_centers_basis_I_small))
     allocate(kinetic_wave_s(n_max_compute_fns_ham), STAT=info )
     call check_allocation(info, 'kinetic_wave_s                ')
     allocate(radial_wave_deriv_s(n_max_compute_fns_ham), STAT=info )
     call check_allocation(info, 'radial_wave_deriv_s            ')
  endif

  if (flag_rel.eq.REL_zora.or.flag_rel==REL_KOLNING_HARMON ) then
     ! allocate all arrays relevant for ZORA

     if (.not.allocated(dist_tab_full)) then
        allocate(dist_tab_full(n_centers_integrals),STAT=info )
        call check_allocation(info, 'dist_tab_full                 ')

     end if
     if (.not.allocated(dir_tab_full_norm)) then
        allocate(dir_tab_full_norm(3,n_centers_integrals),STAT=info )
        call check_allocation(info, 'dir_tab_full_norm             ')
     end if
     if (.not.allocated(i_r_full)) then
        allocate(i_r_full(n_centers_integrals),STAT=info )
        call check_allocation(info, 'i_r_full                      ')

     end if

     if (.not.allocated(zora_vector1)) then
        call aims_allocate( zora_vector1, n_max_compute_ham,3,n_max_batch_size,n_spin, "zora_vector1" )
     end if
     if (.not.allocated(zora_vector2)) then
        call aims_allocate( zora_vector2, n_max_compute_ham,3,n_max_batch_size,n_spin, "zora_vector2" )
     end if

  end if

  !-----------------------------------------------------------------------------


  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing) or to standard batches / arrays (no load balancing)

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

    ! CC: Save local (!) v_xc(r) * rho(r) if Harris-Foulkes-like energy density shall be computed
    if (flag_harris_foulkes_energy_density) then
      if(.not.allocated(ed_xc_pot_energy_density)) allocate(ed_xc_pot_energy_density(batch_perm(n_bp)%n_full_points))
      ed_xc_pot_energy_density(:) = 0.0d0
    end if

    my_n_full_points = batch_perm(n_bp)%n_full_points

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    hartree_potential => hartree_potential_std
    rho_gradient => rho_gradient_std
    kinetic_density => kinetic_density_std

    ld_hamiltonian = n_hamiltonian_matrix_size  ! n_hamiltonian_matrix_size = (n_centers_basis_I+1)*n_centers_basis_I/2

    my_n_full_points = n_full_points

  endif

  ! BL: For some reason n_full_points and the sum of
  ! full points for the load balance case can differ.
  ! This is potentially dangerous but seems
  ! necessary to keep the load balancing working:
  ! When batches are swapped the temporary array for switching
  ! need both full_points information.
  ! However: This means whenever a new array with the dimension n_full_points
  ! is created, one has to ensure that for the lb case n_full_points has to
  ! been taken from batch_perm()%n_full_points rather then just taking the
  ! n_full_points
  ! As a work around I introduced a local variable my_n_full_points, which
  ! give allocated arrays the right dimension.

  ! CC: Save local (!) v_xc(r) * rho(r) if Harris-Foulkes-like energy density shall be computed
  if (flag_harris_foulkes_energy_density) then
    if(.not.allocated(ed_xc_pot_energy_density)) allocate(ed_xc_pot_energy_density(my_n_full_points))
    ed_xc_pot_energy_density(:) = 0.0d0
  end if

  allocate(vdw_potential(my_n_full_points),stat=info)
  call check_allocation(info, 'vdw_potential                 ')

  allocate(batch_times(n_my_batches_work))

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

  hamiltonian(1:ld_hamiltonian*n_spin) = 0.d0
  ! (Rundong) dirac_v t w correspond to the "hamiltonian" in nonrel cases.
  ! We save them in module rel_x2c_mod.
  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
    prune_basis_once = .false.
    dirac_v = (0.d0,0.d0)
    if(scf_iteration.eq.0)then
      dirac_ss   = (0.d0,0.d0)
      dirac_ssc2 = (0.d0,0.d0)
    endif
    if(scf_iteration.eq.0 .or. upw)then
      dirac_w = (0.d0,0.d0)
    endif
    if(scf_iteration.eq.0.and.flag_rel.eq.REL_x2c)then
      dirac_t = (0.d0,0.d0)
    endif
  endif

  if (use_gpu) then
    index_on_gpu = use_batch_permutation > 0 .or. packed_matrix_format == PM_none;

    ! Write out GPU memory usage to stdout
    ! TODO:  Update this block to use output_mem_array_gpu in aims_gpu module
    write (info_str,'(A)') 'Reporting GPU memory usage by myid 0.  Note that it is likely not the only MPI task binding to its GPU!'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(n_max_compute_ham) *8/1.d6,' MB on GPU for hamiltonian_shell'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(n_max_batch_size)*8/1.d6,' MB on GPU for wave'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_batch_size) *8/1.d6,' MB on GPU for partition'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(n_max_batch_size) *8/1.d6,' MB on GPU for h_times_psi'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    if (use_meta_gga) then
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(3*n_max_batch_size) *8/1.d6,' MB on GPU for left_side_of_mgga_dot_product'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(3*n_max_batch_size) *8/1.d6,' MB on GPU for gradient_basis_wave_store'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if
    if (flag_rel.eq.REL_zora.or.flag_rel==REL_KOLNING_HARMON) then
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*3.0d0*dble(n_max_batch_size) *8/1.d6,' MB on GPU for zora_vector1'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*3.0d0*dble(n_max_batch_size) *8/1.d6,' MB on GPU for zora_vector2'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if
    if (index_on_gpu) then
      ! Indexing of the Hamiltonian shell into the Hamiltonian will be done on the GPU
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(ld_hamiltonian)*dble(n_spin)*8/1.d6,' MB on GPU for hamiltonian'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') 'Allocating ',dble(n_max_compute_ham)*dble(n_max_compute_ham) *4/1.d6,' MB on GPU for map'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if

    call hamiltonian_create_gpu (n_max_compute_ham, n_max_batch_size, ld_hamiltonian, n_spin, fort_log_to_c_int(use_meta_gga), &
         fort_log_to_c_int(flag_rel.eq.REL_zora.or.flag_rel==REL_KOLNING_HARMON), & ! Whether to use ZORA
         fort_log_to_c_int(index_on_gpu))

    if (index_on_gpu) then
      ! Indexing of the Hamiltonian shell into the Hamiltonian will be done on the GPU
      call set_hamiltonian_gpu(hamiltonian, ld_hamiltonian*n_spin)
      allocate ( map(n_centers_basis_I * n_centers_basis_I) )
    else if (packed_matrix_format == PM_index) then
      ! The Hamiltonian shell will be communicated back to CPU, and indexing will be done on the CPU
      allocate ( map(1) )
    else
      ! Shouldn't occur, but just in case.
      call aims_stop("Indicted packed matrix format not supported by CUDA code, exiting.")
    end if
  end if


  i_basis_fns_inv = 0; if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) i_basis_fns_inv_s = 0

  en_xc = 0.d0
  en_pot_xc = 0.d0

  if (solvent_method.eq.SOLVENT_MPB) then
    en_Gnonmf = 0d0
    !SR: en_Gnonmf is still calculated via n_full_points summation
    !this should be integrated in the code below
    if (surface_and_volume_calc.and.mpb_solver_started) then
      !evaluate surface cavity and volume
      call evaluate_cavity_surface_and_volume(rho,rho_gradient)
      en_Gnonmf = 0.0000006423049116300181*alphaplusgamma_mpb*surface_mpb+0.0000339882737736419*beta_mpb*volume_mpb
      !en_Gnonmf = 50.*en_Gnonmf
      if (evaluate_sc_Gnonmf) then
        !this is actually calculated during this routine
        en_pot_Gnonmf = 0d0
      end if
    end if
  end if
  ! The evaluation of SC vdW is done here, before the loop over the batches, then the potential is summed to the total
  ! hamiltonian component by component inside the loop. Another approach could be to compute the vdW potential on the
  ! fly inside the loop over the batches below. In this way we can discard the loop currently needed to compute SC vdW.
  !Initialization
  en_pot_vdw = 0.d0
  en_vdw = 0.d0
  vdw_potential = 0.d0
  if (use_vdw_correction_hirshfeld_sc) then
     !Compute the vdW energy
     call en_vdw_hirshfeld(en_vdw)
     !Compute the vdW potential (a vector), en_pot_vdw is just the integral of vdw_potential over the batches, is performed below
     call calc_vdw_potential(vdw_potential)
  endif

  if (use_mbd_std .and. mbd_self_consistent) then
      call mbd_std_potential(rho_std, en_vdw, vdw_potential)
      call mbd_std_cleanup()
  endif

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
     time_start = mpi_wtime()

     n_compute_c = 0
     n_compute_small = 0

     i_basis = 0; if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) i_basis_small = 0
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
           call tab_atom_centered_coords_p0 (coord_current, dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
                n_centers_integrals, centers_basis_integrals)
           if (flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
              dist_tab_sq_s(:,i_point)=dist_tab_sq(:,i_point); dir_tab_s(:,:,i_point)=dir_tab(:,:,i_point)
           endif

           ! determine which basis functions are relevant at current integration point, and tabulate their indices
           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           if(.not.prune_basis_once)then
              call prune_basis_p2 (dist_tab_sq(1,i_point), n_compute_c, i_basis, &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals)
              if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks) then ! For relativistic small component basis:
                 call prune_basis_p2_rel (dist_tab_sq_s(1,i_point), n_compute_small, i_basis_small, &
                      n_centers_basis_I_small, n_centers_integrals, inv_centers_basis_integrals)
              endif
           endif

        end if
     enddo

     if (prune_basis_once) then ! (Rundong) Need to be further changed for fully-relativistic basis. How?
        n_compute_c = batches_work(i_my_batch)%batch_n_compute
        i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
     end if

     ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
     !call get_timestamps(time_batch_centers, clock_batch_centers)
     call collect_batch_centers_p2 (n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, &
          inv_centers_basis_integrals, n_batch_centers, batch_center)
     if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
       call collect_batch_centers_rel (n_compute_small, i_basis_small, n_centers_basis_I_small, n_centers_integrals, &
            inv_centers_basis_integrals, n_batch_centers_s, batch_center_s)
     endif
     !call get_times(time_batch_centers, clock_batch_centers, tot_time_batch_centers, tot_clock_batch_centers,.true.)
     n_points = i_point

     ! Perform actual integration if more than 0 basis functions are actually relevant on the present angular shell ...

! DB 14.Feb.13. modified
     if (n_compute_c.gt.0 .or. n_compute_small.gt.0) then
!     if (n_compute_c.gt.0 .or. (use_nonlinear_core)) then

        n_rel_points = 0
        i_point = 0

        ! Reset meta_gga variables
        if (use_meta_gga) then
           left_side_of_mgga_dot_product = 0.0d0
           gradient_basis_wave_store = 0.d0
        end if

        ! CC: UGLY! The call evaluate_xc_energy_shell() is not aware of
        !           i_full_point (global index)
        !           but just of the batch-wise non-zero-partition
        !           non-zero-density integration points
        !           To recover the local (!) v_xc(r) * rho(r) term entering
        !           the Harris-Foulkes-like
        !           energy density, we construct a
        !           i_point <-> i_full_point map

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

              if(flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON)then
                 call tab_global_geometry_p0 (dist_tab_sq(1,i_point), dir_tab(1,1,i_point), dist_tab_full, i_r_full, &
                      dir_tab_full_norm, n_centers_integrals,  centers_basis_integrals)
              end if
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
              if (use_gga.or.Gnonmf_calc) then
                 do i_spin = 1, n_spin, 1
                   do i_coord = 1,3,1
                      local_rho_gradient(i_coord,i_spin,i_point) = rho_gradient(i_coord,i_spin,i_full_points)
                   enddo
                 enddo
              end if
              if (use_meta_gga) then
                 do i_spin = 1, n_spin, 1
                    local_kinetic_density(i_spin,i_point) = kinetic_density(i_spin,i_full_points)
                 enddo
              endif

              n_compute_atoms = 0; n_compute_atoms_s = 0
              n_compute_fns = 0; n_compute_fns_s = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed) are stored in a compact spline array that can be accessed
              ! by spline_vector_waves, without any copying and without doing any unnecessary operations.
              ! The price is that the interface is no longer explicit in terms of physical objects. See shrink_fixed_basis() for details regarding
              ! the reorganized spline arrays.
              call prune_radial_basis_p2 (n_max_compute_atoms, n_max_compute_fns_ham, dist_tab_sq(1,i_point), dist_tab(1,i_point), &
                   dir_tab(1,1,i_point), n_compute_atoms, atom_index, atom_index_inv, n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                   i_atom_fns, spline_array_start, spline_array_end, n_centers_integrals, centers_basis_integrals, n_compute_c, &
                   i_basis, n_batch_centers, batch_center, one_over_dist_tab, rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point)

              ! Tabulate distances, unit vectors, and inverse logarithmic grid units for all atoms which are actually relevant
              call tab_local_geometry_p2 (n_compute_atoms, atom_index, dist_tab(1,i_point), i_r)

              ! compute trigonometric functions of spherical coordinate angles of current integration point, viewed from all atoms
              call tab_trigonom_p0 (n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab)

              if((use_gga) .or. (flag_rel.eq.REL_zora) .or. (flag_rel==REL_KOLNING_HARMON) .or. Gnonmf_calc)then
                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0 (trigonom_tab(1,1), basis_l_max, l_ylm_max, n_compute_atoms, atom_index, &
                      ylm_tab(1,1), dylm_dtheta_tab(1,1), scaled_dylm_dphi_tab(1,1) )
              else
                ! tabulate distance and Ylm's w.r.t. other atoms
                call tab_wave_ylm_p0 (n_compute_atoms, atom_index, trigonom_tab, basis_l_max, l_ylm_max, ylm_tab)
              end if

              ! Now evaluate radial functions from the previously stored compressed spline arrays
              call evaluate_radial_functions_p0 (spline_array_start, spline_array_end, n_compute_atoms, n_compute_fns, &
                   dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, basis_wave_ordered, radial_wave, .false., &
                   n_compute_c, n_max_compute_fns_ham)

              ! tabulate total wave function value for each basis function
              call evaluate_waves_p2 (n_compute_c, n_compute_atoms, n_compute_fns, l_ylm_max, ylm_tab, one_over_dist_tab, &
                  radial_wave, wave(1,i_point), rad_index, wave_index, l_index, l_count, fn_atom, n_zero_compute, zero_index_point)

              if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks).and.(scf_iteration.eq.0.or.upw))then ! Small component counterpart.
               ! (Rundong) For fully-relativistic cases, there are more shells than nonrel cases: each p, d, f, g, etc. is split to two shells.
               ! This was already considered in the basis generating process.

               ! For X2C and 4c-DKS, subroutine prune_radial_basis_rel is used for generating small component part only; and in this case,
               ! prune_radial_basis_p2 is used to generate the large comp. part.
                call prune_radial_basis_rel (n_max_compute_atoms, n_max_compute_fns_ham, dist_tab_sq_s(1,i_point), dist_tab_s(1,i_point), &
                  dir_tab_s(1,1,i_point), n_compute_atoms_s, atom_index_s, atom_index_inv_s, n_compute_fns_s, i_basis_fns_s, i_basis_fns_inv_s, &
                  i_atom_fns_s, spline_array_start_s, spline_array_end_s, n_centers_integrals, centers_basis_integrals, n_compute_small, &
                  i_basis_small, n_batch_centers_s, batch_center_s, one_over_dist_tab_s, rad_index_small, wave_index_small, &
                  l_index_small, l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s)

               ! Tabulate distances, unit vectors, and inverse logarithmic grid units for all atoms which are actually relevant
                call tab_local_geometry_p2 (n_compute_atoms_s, atom_index_s, dist_tab_s(1,i_point), i_r_s)
                ! compute trigonometric functions of spherical coordinate angles of current integration point, viewed from all atoms
                call tab_trigonom_p0 (n_compute_atoms_s, dir_tab_s(1,1,i_point), trigonom_tab_s)

                if(use_gga)then
                  ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                   call tab_gradient_ylm_p0 (trigonom_tab_s(1,1), basis_l_max_s, l_ylm_max, n_compute_atoms_s, atom_index_s, &
                      ylm_tab_s(1,1), dylm_dtheta_tab_s(1,1), scaled_dylm_dphi_tab_s(1,1) )
                else
                  ! tabulate distance and Ylm's w.r.t. other atoms
                   call tab_wave_ylm_p0 (n_compute_atoms_s, atom_index_s, trigonom_tab_s, basis_l_max_s, l_ylm_max, ylm_tab_s)
                endif

               ! Now evaluate radial functions from the previously stored compressed spline arrays
                call evaluate_radial_functions_p0 (spline_array_start_s, spline_array_end_s, n_compute_atoms_s, n_compute_fns_s, &
                  dist_tab_s(1,i_point), i_r_s, atom_index_s, i_basis_fns_inv_s, basis_wave_s_ordered, radial_wave_s, .false., &
                  n_compute_small, n_max_compute_fns_ham)
               ! tabulate total wave function value for each basis function
                call evaluate_waves_p2 (n_compute_small, n_compute_atoms_s, n_compute_fns_s, l_ylm_max, ylm_tab_s, one_over_dist_tab_s, &
                  radial_wave_s, wave_s(1,i_point), rad_index_small, wave_index_small, l_index_small, l_count_small, &
                  fn_atom_small, n_zero_compute_s, zero_index_point_s)

                if(scf_iteration.eq.0)then
                call evaluate_radial_functions_p0 (spline_array_start_s, spline_array_end_s, n_compute_atoms_s, n_compute_fns_s, &
                   dist_tab_s(1,i_point), i_r_s, atom_index_s, i_basis_fns_inv_s, basis_kinetic_s_ordered, kinetic_wave_s, .false., &
                   n_compute_small, n_max_compute_fns_ham)
                call evaluate_waves_p2 (n_compute_small, n_compute_atoms_s, n_compute_fns_s, l_ylm_max, ylm_tab_s, one_over_dist_tab_s, &
                  kinetic_wave_s, wave_ts(1,i_point), rad_index_small, wave_index_small, l_index_small, l_count_small, &
                  fn_atom_small, n_zero_compute_s, zero_index_point_s)
                endif
              endif

              ! in the remaining part of the subroutine, some decisions (scalar relativity) depend on the potential;
              ! must therefore evaluate the potential and derived quantities right here

              ! Local exchange-correlation parts of the potential are evaluated
              ! right here, to avoid having to store them separately elsewhere.
              ! For large systems, savings are significant

              if((use_hartree_fock.or.use_meta_gga).and.(first_integration)) then
                 call evaluate_xc (rho(1,i_full_points), rho_gradient(1,1,i_full_points), kinetic_density(1,i_full_points), &
                   en_density_xc(i_point), en_density_x, en_density_c, local_xc_derivs(1,i_point), xc_gradient_deriv(1,1,i_point), &
                   xc_tau_deriv(1,i_point), .true., coord_current)  !SAG
              else
                 if(use_embedding_pp.and.use_nonlinear_core) then
                   call evaluate_xc (rho_inc_partialcore(1,i_full_points), rho_gradient_inc_partialcore(1,1,i_full_points), &
                        kinetic_density(1,i_full_points), en_density_xc(i_point), en_density_x, en_density_c, local_xc_derivs(1,i_point),  &
                        xc_gradient_deriv(1,1,i_point), xc_tau_deriv(1,i_point), .false., coord_current)
                 else
                    call evaluate_xc (rho(1,i_full_points), rho_gradient(1,1,i_full_points), kinetic_density(1,i_full_points), &
                         en_density_xc(i_point), en_density_x, en_density_c, local_xc_derivs(1,i_point), xc_gradient_deriv(1,1,i_point), &
                         xc_tau_deriv(1,i_point), .false., coord_current)
                 endif !using pp embedding
              end if !first_integration.and.orbital_dependent

              if(Gnonmf_calc)then
                call evaluate_Gnonmf (rho(:,i_full_points), rho_gradient(:,:,i_full_points), local_Gnonmf_derivs(i_point), &
                     Gnonmf_gradient_deriv(:,i_point), i_full_points)
              endif

              ! add the judgement by igor for the 2-index integrals for the ci calculations
              if(.not. flag_turn_off_hartree)then
                do i_spin = 1, n_spin, 1
                   local_potential_parts(i_spin) =  hartree_potential(i_full_points) + local_xc_derivs(i_spin,i_point) + vdw_potential(i_full_points)
                   if(Gnonmf_calc)then
                     local_potential_parts(i_spin) = local_potential_parts(i_spin) + local_Gnonmf_derivs(i_point)
                   end if
                   if(use_gga)then
                   !
                   ! You might look at this and scratch your head about where the factor of 4 (below) is from - I did, and
                   ! it's not in the CPC paper. It falls out from the method of calculating the contribution to the potential from the GGA density gradient:
                   !
                   ! The contribution from the gradient of the density, grad(rho), is calculated as 2(xc_gradient_deriv.grad(rho)).grad(\chi_i*\chi_j),
                   ! where \chi_i and \chi_j are basis functions with indices i and j.
                   !
                   ! (Wondering where grad rho is multiplied with xc_gradient_deriv? i.e. Equation 29/30 of the Comp. Phys. Comm. 2009 paper?
                   ! Have a look in here --> evaluate_xc.f90 --> xc.f90
                   ! This is dealt with individually for each xc, so xc_gradient_deriv is ACTUALLY xc_gradient_deriv*grad_rho. AJL)
                   !
                   ! grad ( \chi_i*\chi_j ) = [ \chi_i*grad(\chi_j) + grad(\chi_i)*\chi_j ]
                   !
                   ! Therefore 2(xc_gradient_deriv).grad(\chi_i*\chi_j) = 2(xc_gradient_deriv).(\chi_i*grad(\chi_j)) + 2(xc_gradient_deriv).(grad(\chi_i)*\chi_j)
                   ! and as indices i and j are summed over the same basis functions, this boils down to: 4(xc_gradient_deriv).(grad(\chi_i)*\chi_j)
                   !
                   ! Hence the factor of 4! AJL
                   !
                      sum_of_local_gradients(1:3,i_spin) = xc_gradient_deriv(1:3,i_spin,i_point)*4.d0
                   else
                      sum_of_local_gradients(1:3,i_spin) = 0.d0
                   end if

                   if(Gnonmf_calc)then
                      ! SR: we need only a factor of 2 due to the reason AJL explained above
                      ! the other factor of 2 considered in GGA was specific for th GGA formular (29),
                      ! since the derivative of f_xc was taken with respect to |nabla n|^2 instead of
                      ! nabla n as one would do it for a normal functional derivative
                      sum_of_local_gradients(1:3,i_spin) = sum_of_local_gradients(1:3,i_spin) + Gnonmf_gradient_deriv(1:3,i_point) * 2.d0
                   endif
                enddo
              else
                do i_spin = 1, n_spin, 1
                   local_potential_parts(i_spin) = 0.0d0
                   do i_atom = 1, n_atoms, 1
                     local_potential_parts(i_spin) = local_potential_parts(i_spin) - species_z(species(i_atom))/dsqrt(dist_tab_sq(i_atom,i_point))
                   enddo
                   sum_of_local_gradients(1:3,i_spin) = 0.d0
                enddo
              endif

              ! Integration of the vdW potential to obtain the energy.
              ! Maybe it's not necessary to divide the two cases.
              if(use_embedding_pp.and.use_nonlinear_core) then
                 do i_spin = 1, n_spin, 1
                   en_pot_vdw = en_pot_vdw + vdw_potential(i_full_points) * rho_inc_partialcore(i_spin,i_full_points) * partition_tab(i_full_points)
                 end do
              else
                 do i_spin = 1, n_spin, 1
                   en_pot_vdw = en_pot_vdw + vdw_potential(i_full_points) * rho(i_spin,i_full_points) * partition_tab(i_full_points)
                 end do
              endif

              ! Check whether relativistic corrections are needed at the present point.
              ! The check is based entirely on the local parts of the potential - i.e. in a GGA, the terms due to d(rho*exc)/d(|grad(rho|^2)
              ! is not evaluated. Hopefully this approximation to the full ZORA energy is small.
              if (flag_rel.eq.REL_zora.or. (flag_rel==REL_KOLNING_HARMON)) then
                 ! if we need ZORA, must get the _full_ local geometry in  order to create the superposition of atomic potentials which is used
                 ! to estimate the potential gradient for ZORA
                 call evaluate_pot_superpos_p0 (i_r_full, zora_potential_parts(1), n_centers_integrals, centers_basis_integrals)

                 do i_spin = 1, n_spin, 1
                    ! factor 2.d0 required because a factor 1/2 is already included in kinetic_wave later ...
                    zora_operator(i_spin) = 2.d0 * light_speed_sq / ( 2 * light_speed_sq - zora_potential_parts(i_spin) )
                 enddo
              end if

              if ((use_gga) .or. (use_meta_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON).or. (Gnonmf_calc)) then
                 ! We require the gradient of each basis function.
                 ! Tabulate radial derivatives of those radial functions which are actually non-zero at current point, using vectorized splines.
                 call evaluate_radial_functions_p0 (spline_array_start, spline_array_end, n_compute_atoms, n_compute_fns, &
                      dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, basis_deriv_ordered, radial_wave_deriv(1), .true., &
                      n_compute_c, n_max_compute_fns_ham)

                 ! and finally, assemble the actual gradients
                 ! Note "gradient_basis_wave" should never be used in the calling routines, because the leading dimension could change in the subroutine
                 ! where it is actually calculated.
                 call evaluate_wave_gradient_p2 (n_compute_c, n_compute_atoms, n_compute_fns, one_over_dist_tab, dir_tab(1,1,i_point), &
                   trigonom_tab(1,1), l_ylm_max, ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab, radial_wave, radial_wave_deriv, &
                   gradient_basis_wave(1,1), rad_index, wave_index, l_index, l_count, fn_atom, n_zero_compute, zero_index_point)

                 ! Small component counterpart:
                 if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks).and.(scf_iteration.eq.0.or.upw))then
                    call evaluate_radial_functions_p0 (spline_array_start_s, spline_array_end_s, &
                         n_compute_atoms_s, n_compute_fns_s, dist_tab_s(1,i_point), i_r_s, atom_index_s, &
                         i_basis_fns_inv_s, basis_deriv_s_ordered, radial_wave_deriv_s, .true., &
                         n_compute_small, n_max_compute_fns_ham)
                    call evaluate_wave_gradient_p2 (n_compute_small, n_compute_atoms_s, n_compute_fns_s, &
                         one_over_dist_tab_s, dir_tab_s(1,1,i_point), trigonom_tab_s(1,1), l_ylm_max, &
                         ylm_tab_s, dylm_dtheta_tab_s, scaled_dylm_dphi_tab_s, radial_wave_s, radial_wave_deriv_s, &
                         gradient_basis_wave_s(1,1), rad_index_small, wave_index_small, l_index_small, l_count_small, &
                         fn_atom_small, n_zero_compute_s, zero_index_point_s)
                 endif
              end if

              ! Now, evaluate vector of components H*phi(i,r)
              ! Local potential parts first; in the case of GGA, the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate H*phi(i,r) all in one

              ! First, obtain radial kinetic energy terms from  vectorized splines
              call evaluate_radial_functions_p0 (spline_array_start, spline_array_end, n_compute_atoms, n_compute_fns, &
                   dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, basis_kinetic_ordered, kinetic_wave(1), .false., &
                   n_compute_c, n_max_compute_fns_ham)

              do i_spin = 1, n_spin, 1
                 if(flag_rel.eq.REL_4c_dks)then ! Q4C algorithm:
                   ! T+V matrix integration:
                   call evaluate_H_psi_p2 (n_compute_c, n_compute_atoms, n_compute_fns,&
                   l_ylm_max, ylm_tab, one_over_dist_tab, radial_wave, V_times_psi(1,i_point,i_spin),&
                   local_potential_parts(i_spin), kinetic_wave, zora_operator(i_spin), rad_index, wave_index,&
                   l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
                   ! W+Ts matrix integration:
                   if(scf_iteration.eq.0.or.upw)then ! if update W matrix
                      call evaluate_H_psi_p2 (n_compute_small, n_compute_atoms_s, n_compute_fns_s,&
                      l_ylm_max, ylm_tab_s, one_over_dist_tab_s, radial_wave_s, W_times_psi(1,i_point,i_spin),&
                      local_potential_parts(i_spin), kinetic_wave_s, zora_operator(i_spin), rad_index_small, wave_index_small,&
                      l_index_small, l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s)
                   endif
                   if(scf_iteration.eq.0)then ! Ss matrix is only needed to be generated once.
                      ! Ss matrix integration
                      call evaluate_VTW_psi_rel (3,n_compute_small, n_compute_atoms_s, n_compute_fns_s, &
                      l_ylm_max, ylm_tab_s, one_over_dist_tab_s, radial_wave_s, one_times_psi(1,i_point,i_spin),&
                      local_potential_parts(i_spin), kinetic_wave, rad_index_small, wave_index_small, &
                      l_index_small, l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s)
                   endif
                 elseif(flag_rel.eq.REL_x2c)then ! X2C algorithm:
                   ! V matrix integration
                   call evaluate_VTW_psi_rel (1,n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab, radial_wave, V_times_psi(1,i_point,i_spin),&
                   local_potential_parts(i_spin), kinetic_wave, rad_index, wave_index, &
                   l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
                   ! W matrix integration:
                   if(scf_iteration.eq.0.or.upw)then ! if update W matrix
                      call evaluate_VTW_psi_rel (4,n_compute_small, n_compute_atoms_s, n_compute_fns_s,&
                      l_ylm_max, ylm_tab_s, one_over_dist_tab_s, radial_wave_s, W_times_psi(1,i_point,i_spin),&
                      local_potential_parts(i_spin), kinetic_wave_s, rad_index_small, wave_index_small, &
                      l_index_small, l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s)
                   endif
                   if(scf_iteration.eq.0)then ! no need to generate TS matrice in the SCF loop.
                      ! T matrix integration
                      call evaluate_VTW_psi_rel (2,n_compute_c, n_compute_atoms, n_compute_fns, &
                      l_ylm_max, ylm_tab, one_over_dist_tab,  radial_wave, T_times_psi(1,i_point,i_spin),&
                      local_potential_parts(i_spin), kinetic_wave, rad_index, wave_index, &
                      l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
                      ! Ss matrix integration
                      call evaluate_VTW_psi_rel (3,n_compute_small, n_compute_atoms_s, n_compute_fns_s, &
                      l_ylm_max, ylm_tab_s, one_over_dist_tab_s, radial_wave_s, one_times_psi(1,i_point,i_spin),&
                      local_potential_parts(i_spin), kinetic_wave, rad_index_small, wave_index_small, &
                      l_index_small, l_count_small, fn_atom_small, n_zero_compute_s, zero_index_point_s)
                   endif
                 else
                   call evaluate_H_psi_p2 (n_compute_c, n_compute_atoms, n_compute_fns,&
                   l_ylm_max, ylm_tab, one_over_dist_tab, radial_wave, H_times_psi(1,i_point,i_spin),&
                   local_potential_parts(i_spin), kinetic_wave, zora_operator(i_spin), rad_index, wave_index,&
                   l_index, l_count, fn_atom, n_zero_compute, zero_index_point)
                 endif
              enddo

              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
              if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks) i_basis_fns_inv_s(:,atom_index_s(1:n_compute_atoms_s)) = 0

              if ((flag_rel.eq.REL_zora).or. flag_rel==REL_KOLNING_HARMON) then

                 ! Scalar relativistic treatment.
                 ! count number of "truly" relativistic points for ZORA treatment of kinetic energy ...
                 do i_spin = 1, n_spin, 1
                    zora_operator(i_spin) = light_speed_sq / (2 * light_speed_sq - zora_potential_parts(i_spin))**2
                    call add_zora_gradient_part_p0 (sum_of_local_gradients(1,i_spin), i_r_full, dir_tab_full_norm, dist_tab_full, &
                         zora_operator(i_spin), n_centers_integrals, centers_basis_integrals)
                 end do

                 do i_spin = 1, n_spin, 1
                    ! Evaluate difference of scalar relativistic kinetic energy operator for the true potential and the superposition of free atom
                    ! potentials separately, and only for all relativistic points in shell. Here, use partially integrated version, leading to a vector:
                    ! zora_operator(r)*grad(phi(r,i))
                    zora_operator(i_spin) = light_speed_sq * (local_potential_parts(i_spin) - zora_potential_parts(i_spin)) /  &
                        ( 2 * light_speed_sq - local_potential_parts(i_spin)) / ( 2 * light_speed_sq - zora_potential_parts(i_spin))
                    call evaluate_zora_vector_p1 (zora_operator(i_spin), partition_tab(i_full_points), gradient_basis_wave(1,1),  &
                         n_compute_c, zora_vector1(1, 1, n_rel_points+1, i_spin), zora_vector2(1, 1, n_rel_points+1, i_spin), &
                         n_max_compute_ham, t_zora(i_spin))
                 enddo

                 if (n_spin.eq.1) then
                   if(t_zora(1)) then
                      n_rel_points = n_rel_points + 1
                   end if
                 else if (n_spin.eq.2) then
                   if(t_zora(1).or.t_zora(2)) then
                      n_rel_points = n_rel_points + 1
                   end if
                 end if

              end if  ! end ZORA preparations

              ! If using a GGA, add the true gradient terms to the Hamiltonian vector
              ! if (use_gga .or. (n_rel_points.gt.0)) then
              if (use_gga .or.(flag_rel.eq.REL_zora) .or. flag_rel==REL_KOLNING_HARMON.or. (Gnonmf_calc)) then
                 ! Calculates sum_of_local_gradients*grad(\chi)
                 ! Multiplied by \chi later on to complete the Hamiltonian. AJL
                 do i_spin = 1, n_spin, 1
                    if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then ! Fully relativistic cases
                       call add_gradient_part_to_H_p0 (n_compute_c, gradient_basis_wave(1,1), &
                          sum_of_local_gradients(1,i_spin), V_times_psi(1, i_point, i_spin))
                       ! This is valid for Q4C. I don't guarantee X2C for now:
                       if(scf_iteration.eq.0.or.upw) call add_gradient_part_to_H_p0 &
                          (n_compute_small, gradient_basis_wave_s(1,1), &
                          sum_of_local_gradients(1,i_spin), W_times_psi(1, i_point, i_spin))
                    else
                       call add_gradient_part_to_H_p0 (n_compute_c, gradient_basis_wave(1,1), &
                            sum_of_local_gradients(1,i_spin), H_times_psi(1, i_point, i_spin))
                    endif
                 enddo
              end if

              ! We store the components here to calculate grad(chi).grad(chi)
              ! Left side contains xc_tau_deriv and partition weight. AJL
              if (use_meta_gga.and..not.first_integration) then

                 ! Prepare left side of dot product with gradient_basis_wave
                 do i_spin = 1, n_spin, 1
                   call evaluate_mgga_left_side_of_dot_product (n_compute_c, n_points, i_point, partition(i_point), xc_tau_deriv(i_spin, i_point), &
                      gradient_basis_wave(1,1), left_side_of_mgga_dot_product(1,1,i_spin) )
                 enddo

                 ! Store right side of dot product as gradient_basis_wave_store
                 call store_gradient_basis_wave (n_compute_c, n_points, i_point, gradient_basis_wave(1,1), gradient_basis_wave_store(1,1))
              endif ! use_meta_gga
           end if  ! end if (partition_tab(i_full_points).gt.0)
        enddo ! end loop over a batch

        ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications work separately for each spin channel
        do i_spin = 1, n_spin, 1

           ! add full non-relativistic contributions and (for relativistic points) all contributions from the potential to the Hamiltonian matrix elements
           if (use_gpu) then
              call evaluate_hamiltonian_shell_gpu (n_points, partition(1), n_compute_c, &
                   H_times_psi(1,1,i_spin), n_max_compute_ham, wave(1,1), hamiltonian_shell)
           else
              if(flag_rel.eq.REL_4c_dks)then ! Q4C:
                 if(i_spin.eq.2)cycle ! i_spin.eq.2 is not used for now, an interface is left here
                 call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_c, &
                      V_times_psi(1,1,i_spin), n_max_compute_ham, wave(1,1), dirac_v_shell)
                 if(scf_iteration.eq.0.or.upw)then
                    call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_small, &
                         W_times_psi(1,1,i_spin), n_max_compute_ham, wave_s(1,1), dirac_w_shell)
                 endif
                 if(scf_iteration.eq.0)then
                    call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_small, &
                         one_times_psi(1,1,i_spin), n_max_compute_ham, wave_s(1,1), dirac_ss_shell)
                 endif
              elseif(flag_rel.eq.REL_x2c)then ! X2C:
                 if(i_spin.eq.2)cycle ! i_spin.eq.2 is not used for now, an interface is left here
                 call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_c, &
                      V_times_psi(1,1,i_spin), n_max_compute_ham, wave(1,1), dirac_v_shell)
                 if(scf_iteration.eq.0.or.upw)then
                    W_times_psi(1:n_max_compute_ham,1:n_points,i_spin) = &
                         W_times_psi(1:n_max_compute_ham,1:n_points,i_spin)*0.25d0/light_speed_sq
                    call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_small, &
                         W_times_psi(1,1,i_spin), n_max_compute_ham, wave_ts(1,1), dirac_w_shell)
                 endif
                 if(scf_iteration.eq.0)then
                    call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_c, &
                         T_times_psi(1,1,i_spin), n_max_compute_ham, wave(1,1), dirac_t_shell)
                    call evaluate_hamiltonian_shell_p1 (n_points, partition, n_compute_small, &
                         one_times_psi(1,1,i_spin), n_max_compute_ham, wave_s(1,1), dirac_ss_shell)
                 endif
              else ! Non-rel cases:
                call evaluate_hamiltonian_shell_p1 (n_points, partition(1), n_compute_c, &
                     H_times_psi(1,1,i_spin), n_max_compute_ham, wave(1,1), hamiltonian_shell)
              endif
           end if

           ! AJL. For meta_gga, add d/d(tau) contributions
           if (use_meta_gga.and..not.first_integration) then
             if (use_gpu) then
                call mgga_contribution_gpu (n_compute_c, n_compute_c, n_points, left_side_of_mgga_dot_product(1,1,i_spin), gradient_basis_wave_store(1,1) )
             else
                call evaluate_mgga_contribution_and_add_to_hamiltonian_shell (n_compute_c, n_compute_c, n_points, left_side_of_mgga_dot_product(1,1,i_spin), &
                    gradient_basis_wave_store(1,1), hamiltonian_shell )
             endif
           endif

           ! For all relativistic (ZORA only) points, add kinetic energy contributions
           if (n_rel_points.gt.0) then
              if (use_gpu) then
                 call add_zora_matrix_gpu (zora_vector1(1,1,1,i_spin), zora_vector2(1,1,1,i_spin), n_max_compute_ham, n_rel_points, n_compute_c)
              else
                 call add_zora_matrix_p1 (zora_vector1(1,1,1,i_spin), zora_vector2(1,1,1,i_spin), n_max_compute_ham, n_rel_points, n_compute_c, hamiltonian_shell)
              end if
           end if

           if(use_batch_permutation > 0) then
              ! If use_batch_permutation > 0 is set, the local hamiltonian is always stored in full form for the local basis functions
              ! Get position of basis functions of current batch within local hamiltonian
              do i=1,n_compute_c
                 ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
              enddo

              ! Insert hamiltonian_shell of current batch
              if(use_gpu) then
                 ! For load-balancing Hamiltonians, we update the local part of the Hamiltonian using the batch matrix directly on the GPU,
                 ! as the indexing scheme is trivial.  We communicate the local part of the Hamiltonian back to the CPU at the end of this
                 ! subroutine, after all batches have been processed
                 call update_batch_matrix_gpu (i_spin, ld_hamiltonian, ins_idx, n_compute_c)
              else
                 do i=1,n_compute_c
                    i_off = (i_spin-1)*ld_hamiltonian + (ins_idx(i)*(ins_idx(i)-1))/2
                    do j=1,i ! n_compute_c
                       if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! Fully-relativistic cases:
                         !dirac_vtw_sum(ins_idx(j)+i_off) = dirac_v_sum(ins_idx(j)+i_off) + dirac_v_shell(j+(i-1)*n_compute_c)
                       else
                          hamiltonian(ins_idx(j)+i_off) = hamiltonian(ins_idx(j)+i_off) + hamiltonian_shell(j+(i-1)*n_compute_c)
                       endif
                       if (1 == 2) then
                         print *,"i = ", i, " j = ", j, " ShellIdx = ", j+(i-1)*n_compute_c, " MatIdx = ", ins_idx(j)+i_off, &
                         " Shell = ", hamiltonian_shell(j+(i-1)*n_compute_c), " Matrix = ",  hamiltonian(ins_idx(j)+i_off)
                      end if
                    enddo
                 enddo
              end if
           else
              if (use_gpu) then
                  if (packed_matrix_format == PM_none) then
                    ! For non-packed Hamiltonians, we update the Hamiltonian using the batch matrix directly on the GPU, as the indexing
                    ! scheme is trivial.  We communicate the Hamiltonian back to the CPU at the end of this subroutine, after all batches
                    ! have been processed
                    !
                    ! This could introduce a memory bottleneck on the GPU, since the full Hamiltonian is stored on the GPU.
                    call get_update_full_matrix_map (i_spin, ld_hamiltonian, n_compute_c, i_basis, map)
                    call update_full_matrix_via_map_gpu (hamiltonian, ld_hamiltonian*n_spin, hamiltonian_shell, n_compute_c, n_compute_c, &
                         map, n_compute_c * n_compute_c)
                  else if (packed_matrix_format == PM_index) then
                    ! For packed and local-indexed-but-not-load-balanced Hamiltonians, the indexing scheme from the batch matrix to
                    ! the Hamiltonian is non-trivial.  Furthermore, for packed Hamiltonians, the Hamiltonian matrix will quickly become a
                    ! memory bottleneck.
                    !
                    ! Accordingly, in this case we transfer the batch matrix back to the CPU and let it do the indexing back to the
                    ! Hamiltonian matrix.  This introduces an additional data transfer each cycle and will impact performance.
                    call get_hamiltonian_shell_gpu(hamiltonian_shell, n_max_compute_ham)
                    call update_full_matrix_p0X (n_compute_c, n_compute_c, i_basis(1), hamiltonian_shell, hamiltonian((i_spin-1)*ld_hamiltonian+1))
                  else
                    call aims_stop("Packed matrix format not supported by CUDA code, exiting.")
                  end if
              else
                 if(flag_rel.eq.REL_4c_dks)then ! Q4C:
                    call update_full_matrix_rel (n_compute_c, i_basis, ld_dirac_large, &
                         dirac_v_shell, dirac_v_sum((i_spin-1)*ld_dirac_large+1))
                    if(scf_iteration.eq.0.or.upw)then
                       call update_full_matrix_rel (n_compute_small, i_basis_small, ld_dirac_small, &
                            dirac_w_shell, dirac_w_sum((i_spin-1)*ld_dirac_small+1))
                    endif
                    if(scf_iteration.eq.0)then
                       call update_full_matrix_rel (n_compute_small, i_basis_small, ld_dirac_small, &
                            dirac_ss_shell, dirac_ss_sum((i_spin-1)*ld_dirac_small+1))
                    endif
                 elseif(flag_rel.eq.REL_x2c)then ! X2C:
                    call update_full_matrix_rel (n_compute_c, i_basis, ld_dirac_large, &
                         dirac_v_shell, dirac_v_sum((i_spin-1)*ld_dirac_large+1))
                    if(scf_iteration.eq.0.or.upw)then
                       call update_full_matrix_rel (n_compute_small, i_basis_small, ld_dirac_small, &
                            dirac_w_shell, dirac_w_sum((i_spin-1)*ld_dirac_small+1))
                    endif
                    if(scf_iteration.eq.0)then
                       call update_full_matrix_rel (n_compute_c, i_basis, ld_dirac_large, &
                            dirac_t_shell, dirac_t_sum((i_spin-1)*ld_dirac_large+1))
                       call update_full_matrix_rel (n_compute_small, i_basis_small, ld_dirac_small, &
                            dirac_ss_shell, dirac_ss_sum((i_spin-1)*ld_dirac_small+1))
                    endif
                 else
                   call update_full_matrix_p0X (n_compute_c, n_compute_c, i_basis, hamiltonian_shell, &
                        hamiltonian((i_spin-1)*ld_hamiltonian+1))
                 endif ! end if flag_rel
              end if ! end of use_gpu
           endif
        enddo ! end of i_spin

        ! Hamiltonian is now complete.
        !
        ! Since we already have the pieces, add terms of XC energy here.
        ! Notice that these terms are not added for ANY shell where n_compute happens to be zero. This should be correct
        ! because all wave functions are zero here anyway, i.e. also the density.

        call evaluate_xc_energy_shell (n_points, energy_partition, en_density_xc, local_xc_derivs, &
             xc_gradient_deriv, xc_tau_deriv, local_rho, local_rho_gradient, local_kinetic_density, en_xc, en_pot_xc)

        if(Gnonmf_calc)then
          call evaluate_Gnonmf_energy_shell (n_points, energy_partition, local_Gnonmf_derivs, Gnonmf_gradient_deriv, local_rho, local_rho_gradient, &
               en_pot_Gnonmf)
        end if

        ! CC: Everything has been mapped in evaluate_xc_energy_shell above
        !     i_point <-> i_full_point map (see comment above)
        !     not required anymore
        if (flag_harris_foulkes_energy_density) then
          if(allocated(ed_i_point_to_i_full_point_map)) deallocate(ed_i_point_to_i_full_point_map)
        end if

     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches


  !Get the total vdW potential energy
  if ( use_vdw_correction_hirshfeld_sc .or. (use_mbd_std .and. mbd_self_consistent) ) then
     call sync_real_number(en_pot_vdw)
  endif

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  ! Synchronise the hamiltonian
  if (use_gpu) then
     ! When using GPU acceleration for either the non-packed or load-balanced Hamiltonian, the Hamiltonian is stored and updated directly on the GPU.
     ! So we transfer it back to the CPU here.
     !
     ! For the packed or local-indexed-but-non-load-balanced Hamiltonian, we get the batch matrix back to the CPU after every evaluation and do the
     ! indexing into the Hamiltonian on the CPU, so no transfer needs to be done here.
     if (index_on_gpu) then
       call get_hamiltonian_gpu( hamiltonian, ld_hamiltonian*n_spin )
     end if

     write (info_str,'(A)')'Deallocating memory used on GPU for Hamiltonian integration.'
     call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
     call hamiltonian_destroy_gpu()
  end if
  use_gpu = gpu_save

  if(.not. use_local_index)then
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
        call sync_vector(dirac_v_sum,    ld_dirac_large)
        if(scf_iteration.eq.0.or.upw)                    call sync_vector(dirac_w_sum,  ld_dirac_small)
        if(scf_iteration.eq.0)                           call sync_vector(dirac_ss_sum, ld_dirac_small)
        if(scf_iteration.eq.0 .and. flag_rel.eq.REL_x2c) call sync_vector(dirac_t_sum,  ld_dirac_large)
     else
        call sync_integrate_hamiltonian( hamiltonian )
     endif
  endif

  if (condition_penalty > 0.d0) then
     if(use_batch_permutation > 0) call aims_stop('condition_penalty not implemented yet for load balancing')
     call add_const_to_hamiltonian(hamiltonian, condition_penalty)
  end if

   if(get_batch_weights) then
      call set_batch_weights(n_bp, batch_times)
   endif

  ! synchronise the XC energy / potential contributions
  call sync_update_xc_potential( en_xc, en_pot_xc )

  if (Gnonmf_calc) then
  !en_Gnonmf does not need to be synced but we should integrate the evaluation of en_Gnonmf in the code above should be more efficient
  !than going over the grid several times
    call sync_update_xc_potential( dummy, en_pot_Gnonmf )
  end if

  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! Now, transform the obtained scalar integrations to spinor integrations
    write(info_str,'(2X,A)')'spinor V matrix'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    call transform_scalar2spinor(1, 1, n_centers_basis_I, ld_dirac_large, dirac_v_sum, dirac_v)
    if(scf_iteration.eq.0.or.upw)then
      write(info_str,'(2X,A)')'spinor W matrix'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      call transform_scalar2spinor(4, 1, n_centers_basis_I_small, ld_dirac_small, dirac_w_sum, dirac_w)
      call aims_deallocate( dirac_w_sum, "dirac_w_sum" ) ! This could be a very large array, release to save the memory
    endif
    if(scf_iteration.eq.0)then
      write(info_str,'(2X,A)')'spinor Ss matrix'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      call transform_scalar2spinor(5, 1, n_centers_basis_I_small, ld_dirac_small, dirac_ss_sum, dirac_ss)
      write(info_str,'(2X,A)')'spinor Ssc2 matrix'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
     ! Now, dirac_Ss is generated, use dirac_ss_sum to save the values of dirac_ssc2_sum
      dirac_ss_sum = -2.d0*dirac_ss_sum*light_speed_sq
      call transform_scalar2spinor(5, 1, n_centers_basis_I_small, ld_dirac_small, dirac_ss_sum, dirac_ssc2)
      call aims_deallocate( dirac_ss_sum, "dirac_ss_sum" ) ! This could be a very large array, release to save the memory
    endif
    if(scf_iteration.eq.0 .and. flag_rel.eq.REL_x2c)then
      write(info_str,'(2X,A)')'spinor T matrix'
      call localorb_info(info_str,use_unit,'(A)',OL_norm)
      call transform_scalar2spinor(2, 1, n_centers_basis_I, ld_dirac_large, dirac_t_sum, dirac_t)
    endif
  endif

  if (first_integration) then

     if (out_grids) then
        ! write the entire integration grid in the format required for a re-readin in the file grids.dat
        if (myid.eq.0) then
           open (50,file="grids_out.dat")
           write(50,'(I5)') n_atoms
           do i_atom = 1, n_atoms, 1
              write(50,'(2X,I5,1X,I5)') i_atom, n_radial(species(i_atom))
              do i_radial = 1,n_radial(species(i_atom))
                 write(50,'(4X,I5,1X,E30.15,1X,E30.15,1X,I5)') i_radial, r_radial(i_radial,species(i_atom)),  &
                      w_radial(i_radial,species(i_atom)), n_angular(i_radial,species(i_atom))
                 do i_angular = 1, n_angular(i_radial, species(i_atom))
                    write(50,'(6X,I5,1X,3(1X,E30.15),1X,E30.15)') i_angular, (r_angular(i_coord, i_angular, i_radial,  &
                         species(i_atom)), i_coord = 1,3,1), w_angular(i_angular, i_radial, species(i_atom))
                 enddo
              enddo
           enddo
           close(50)
        end if

        ! write integration grid per species in the format needed for control.in
        if (myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A,A)') "Output of integration grids in suitable form ", "for copy-paste into control.in:"
           do i_species = 1, n_species, 1
              write(use_unit,*)
              write(use_unit,'(2X,A,A,A)') "Species ", species_name(i_species), ":"
              do i_radial = 2, n_radial(i_species), 1
                 if (n_angular(i_radial,i_species).ne.(n_angular(i_radial-1,i_species))) then
                    write(use_unit,'(6X,A,1X,F8.4,1X,I4)') "division", ( r_radial(i_radial-1,i_species) + 0.001 ) * bohr, n_angular(i_radial-1,i_species)
                 end if
              enddo
              write(use_unit,'(6X,A,1X,I4)') "outer_grid", n_angular(i_radial-1,i_species)
           enddo
           write(use_unit,*)
        end if

     end if

  end if

  ! Allocatable arrays that are tracked
  if(allocated( zora_vector2         ))         call aims_deallocate( zora_vector2,                                   "zora_vector2" )
  if(allocated( zora_vector1         ))         call aims_deallocate( zora_vector1,                                   "zora_vector1" )
  if(allocated( H_times_psi          ))         call aims_deallocate( H_times_psi,                                     "H_times_psi" )
  if(allocated( V_times_psi          ))         call aims_deallocate( V_times_psi,                                     "V_times_psi" )
  if(allocated( T_times_psi          ))         call aims_deallocate( T_times_psi,                                     "T_times_psi" )
  if(allocated( W_times_psi          ))         call aims_deallocate( W_times_psi,                                     "W_times_psi" )
  if(allocated( one_times_psi        ))         call aims_deallocate( one_times_psi,                                 "one_times_psi" )
  if(allocated( hamiltonian_shell    ))         call aims_deallocate( hamiltonian_shell,                         "hamiltonian_shell" )
  if(allocated( dirac_v_shell        ))         call aims_deallocate( dirac_v_shell,                                 "dirac_v_shell" )
  if(allocated( dirac_t_shell        ))         call aims_deallocate( dirac_t_shell,                                 "dirac_t_shell" )
  if(allocated( dirac_w_shell        ))         call aims_deallocate( dirac_w_shell,                                 "dirac_w_shell" )
  if(allocated( dirac_ss_shell       ))         call aims_deallocate( dirac_ss_shell,                               "dirac_ss_shell" )
  if(allocated( gradient_basis_wave_store))     call aims_deallocate( gradient_basis_wave_store,         "gradient_basis_wave_store" )
  if(allocated( left_side_of_mgga_dot_product)) call aims_deallocate( left_side_of_mgga_dot_product, "left_side_of_mgga_dot_product" )
  if(allocated( wave                 ))         call aims_deallocate( wave,                                                   "wave" )
  if(allocated( wave_s               ))         call aims_deallocate( wave_s,                                               "wave_s" )
  if(allocated( wave_ts              ))         call aims_deallocate( wave_ts,                                             "wave_ts" )
  if(allocated( sq_psi               ))         call aims_deallocate( sq_psi,                                               "sq_psi" )
  if(allocated( gradient_basis_wave  ))         call aims_deallocate( gradient_basis_wave,                     "gradient_basis_wave" )
  if(allocated( gradient_basis_wave_s))         call aims_deallocate( gradient_basis_wave_s,                 "gradient_basis_wave_s" )

  if(allocated( i_r_full             )) deallocate( i_r_full             )
  if(allocated( dir_tab_full_norm    )) deallocate( dir_tab_full_norm    )
  if(allocated( dist_tab_full        )) deallocate( dist_tab_full        )
  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( i_basis_small        )) deallocate( i_basis_small        )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave_deriv_s  )) deallocate( radial_wave_deriv_s  )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( radial_wave_s        )) deallocate( radial_wave_s        )
  ! MGGA allocations
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( ylm_tab_s            )) deallocate( ylm_tab_s            )
  if(allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if(allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if(allocated( dylm_dtheta_tab_s    )) deallocate( dylm_dtheta_tab_s    )
  if(allocated( scaled_dylm_dphi_tab_s))deallocate(scaled_dylm_dphi_tab_s)
  if(allocated( local_kinetic_density)) deallocate( local_kinetic_density)
  if(allocated( local_rho_gradient   )) deallocate( local_rho_gradient   )
  if(allocated( local_rho            )) deallocate( local_rho            )
  if(allocated( xc_tau_deriv         )) deallocate( xc_tau_deriv         )
  if(allocated( xc_gradient_deriv    )) deallocate( xc_gradient_deriv    )
  if(allocated( local_xc_derivs      )) deallocate( local_xc_derivs      )
  if(allocated( en_density_xc        )) deallocate( en_density_xc        )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_atom_fns_s         )) deallocate( i_atom_fns_s         )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns_inv_s    )) deallocate( i_basis_fns_inv_s    )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( i_basis_fns_s        )) deallocate( i_basis_fns_s        )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )
  if(allocated( dir_tab_s            )) deallocate( dir_tab_s            )
  if(allocated( dist_tab_sq_s        )) deallocate( dist_tab_sq_s        )
  if(allocated( dist_tab_s           )) deallocate( dist_tab_s           )
  if(allocated( vdw_potential        )) deallocate( vdw_potential        )
  if(allocated( map                  )) deallocate( map                  )

  if (Gnonmf_calc) then
    if (allocated( local_Gnonmf_derivs )) deallocate( local_Gnonmf_derivs)
    if (allocated( Gnonmf_gradient_deriv)) deallocate ( Gnonmf_gradient_deriv)
  end if

  if(allocated( rho_inc_partialcore  )) then
     deallocate( rho_inc_partialcore  )
  endif

  if(allocated( rho_gradient_inc_partialcore )) then
     deallocate( rho_gradient_inc_partialcore )
  endif

  if(use_batch_permutation > 0) then
    deallocate(rho)
    deallocate(hartree_potential)
    deallocate(kinetic_density) ! always allocated
    deallocate(rho_gradient) ! always allocated
    deallocate(ins_idx)
  endif

  if (use_nlcorr_in_xc) then
     call deallocate_vdw_splines()
  end if

  if(allocated(batch_times)) deallocate(batch_times)

  if(.not.use_local_index) then
    ! Please note if add_plus_u_to_hamiltonian should be implemented for use_local_index:
    ! The size of the hamiltonian may be different from n_hamiltonian_matrix_size
    ! if loadbalancing is in effect!
    if(use_plus_u)then
      if(plus_u_mulliken_charges .eqv. .true.) then
         call add_plus_u_mulliken_to_hamiltonian(hamiltonian)
      elseif(plus_u_full .eqv. .true.) then
         call add_plus_u_full_to_hamiltonian(hamiltonian)
      else ! use the on-site representation (default)
         call add_plus_u_to_hamiltonian(hamiltonian)
      endif
    endif
  endif


   call get_times (time_total_local, clock_total_local)

   ! synchronize all cumulative CPU-time timestamps
   call sync_timing(time_total_local)

   if (output_level .eq. "full") then
      call output_timeheader(deffmt, 'Integrate hamiltonian matrix')
      call output_times(deffmt, 'Total time', time_total_local, clock_total_local)
   endif
end subroutine integrate_real_hamiltonian_matrix_p2
!******
