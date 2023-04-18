!****s* FHI-aims/update_density_densmat
!  NAME
!    update_density_densmat
!  SYNOPSIS

subroutine update_density_densmat &
     ( KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
       partition_tab_std, hartree_partition_tab_std, rho_std, &
       rho_gradient_std,  kinetic_density_std, basis_l_max, delta_rho_KS_std, &
       delta_rho_gradient_std, delta_kinetic_density_std, &
       rho_change,  density_matrix_sparse )

!  PURPOSE
!
!  The subroutine obtains the new KS density using the density
!  matrix.
!
!  For the benefit of Pulay mixing later on, we obtain the change between
!  the input and output density, delta_rho_KS = rho_KS - rho_in,
!  rather than overwriting rho itself. So the results are in variable delta_rho_KS.
!
!  This routine corresponds to Eq. (28) in Ref.
!  Blum et al., Computer Physics Communications 180 (2009) 2175-2196
!  and is used for all periodic systems as well as for large non-periodic systems
!  in FHI-aims. For small non-periodic systems, see
!  update_density_and_forces_orbital instead.
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
  use analyze_arrays
  use scalapack_wrapper
  use wf_extrapolation
  use KH_core_states
  use plus_u
  use density_matrix_evaluation
  use load_balancing
  use timing
  use sym_base,only:evaluate_densmat_sym, occ_numbers_nosym
  use physics,only:KS_eigenvalue,rho_pce,rho_pce_gradient
  use rel_x2c_mod
  use lpb_solver_utilities, only: atomic_MERM, mpb_solver_started
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use aims_gpu, only: gpu_density_used
  use c_helper, only: fort_log_to_c_int
  use rlsy_interface, only: rlsy_h
  use batch_rlsy_interface, only: batch_irr, &
                                  compute_batch_rlsy_full_to_irr_points, &
                                  compute_batch_rlsy_full_to_irr_points_gradient, &
                                  compute_batch_rlsy_irr_to_full_points, &
                                  compute_batch_rlsy_irr_to_full_points_gradient
  implicit none

!  ARGUMENTS

  real*8,     dimension(n_basis,n_states,n_spin,n_k_points_task), intent(IN):: KS_eigenvector 
  complex*16, intent(IN):: KS_eigenvector_complex(*) ! (n_basis,n_states,n_spin,n_k_points_task) or (2*dim_matrix_rel,n_states,n_spin,n_k_points_task)
  real*8,     dimension(n_states, n_spin, n_k_points)             , intent(IN):: occ_numbers
  ! NOTE:  For this routine to give the correct results, the occ_numbers variable passed in should have already been
  !        properly k-weighted (and are thus not the "true" occupation numbers.)
  !        In scf_solver, this is done by the kweight_occs function (and undone at the end by de_kweight_occs)
  real*8, target, dimension(n_full_points), intent(IN)            :: partition_tab_std
  real*8, target, dimension(n_full_points), intent(IN)            :: hartree_partition_tab_std
  real*8, target, dimension(n_spin, n_full_points), intent(IN)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points), intent(IN) :: rho_gradient_std
  real*8, target, dimension(n_spin, n_full_points), intent(IN) :: kinetic_density_std
  integer, dimension(n_species), intent(IN)                 :: basis_l_max
  real*8, target, intent(OUT) :: delta_rho_KS_std(n_full_points, n_spin)
  real*8, target, intent(OUT) :: delta_rho_gradient_std(3, n_full_points, n_spin)
  real*8, target, intent(OUT) :: delta_kinetic_density_std(n_full_points, n_spin)
  real*8,  dimension(n_spin), intent (OUT)                    :: rho_change
  ! when this routine is called, density_matrix_sparse has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8,  dimension(*)                          :: density_matrix_sparse

! INPUTS
! o KS_eigenvector -- eigenvectors if real eigenvectors are in use
! o KS_eigenvector_complex -- eigenvectors is complex eigenvectors are in use
! o occ_numbers -- occupation numbers of states, with k-weighting already applied
! o partition_tab -- partition function values
! o hartree_partition_tab -- hartree potentials partition function values
! o rho -- electron density, rho is only input; what we store is the density residual (i.e. the change in the density)
! o rho_gradient -- gradient of the electron density
! o kinetic_density -- kinetic electron density
! o basis_l_max -- maximum basis l value of the basis functions
! o density_matrix_sparse -- this is the works space for density matrix, it is used if packed matrixes are in use.
!
! OUTPUT
!   o delta_rho_KS -- calculated charge density  -  old rho (input)
!   o delta_rho_gradient -- calculated gradient  -  old rho_gradient (input)
!   o delta_kinetic_density -- calculated kinetic density - old kinetic_density (input)
!   o rho_change -- how much the electron density changed during the update.
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




  real*8, dimension(:),     allocatable:: dist_tab
  real*8, dimension(:,:),   allocatable:: dist_tab_sq
  real*8, dimension(:),     allocatable:: i_r
  real*8, dimension(:,:,:), allocatable:: dir_tab
  real*8, dimension(:,:),   allocatable:: trigonom_tab

  real*8,dimension(:),  allocatable:: radial_wave
  real*8,dimension(:),  allocatable:: radial_wave_deriv
  real*8,dimension(:),  allocatable:: radial_wave_2nd_deriv
  real*8,dimension(:,:),allocatable:: wave
  real*8,dimension(:),  allocatable:: kinetic_wave

  integer,dimension(:),allocatable :: i_basis
  integer :: n_compute_c

  real*8 coord_current(3)

  integer :: n_compute_fns, i_center
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

  real*8 :: one_over_dist_tab(n_max_compute_atoms)

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

  integer :: l_ylm_max
  integer :: n_points

  ! (n_basis, n_states, n_spin) on purpose
  ! beware: inside the subroutines (evaluate_KS_density_v1 for instance)
  ! dimension is (n_compute, max_occ_number)
  ! that's a trick to get a continuous data flow and not
  ! bad programming

  real*8 :: temp_rho(n_max_batch_size,n_spin)
  real*8 :: temp_rho_gradient(3,n_max_batch_size,n_spin)
  real*8 :: temp_kinetic_density(n_max_batch_size,n_spin)
  real*8,dimension(:,:),allocatable :: temp_rho_small

  integer, dimension(:,:), allocatable :: index_lm
  real*8,  dimension(:,:), allocatable :: ylm_tab

  !     only allocated and referenced for gradient functionals
  real*8, dimension(:,:),  allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:),  allocatable :: scaled_dylm_dphi_tab
  real*8, dimension(:,:,:),allocatable :: gradient_basis_wave

  real*8, dimension(:,:),  allocatable :: density_matrix
  real*8, dimension(:),  allocatable :: density_matrix_con
  real*8, dimension(:,:),  allocatable :: work
  real*8, dimension(:,:),  allocatable :: grad_work
  complex*16, dimension(:,:),allocatable :: work_complex
  real*8, allocatable :: density_matrix_extrapolated(:,:,:)
  real*8, allocatable :: density_matrix_sparse_extrapolated(:,:)

  !     counters
  integer :: i_l
  integer :: i_m
  integer :: i_coord
  integer :: i_state
  integer :: i_point

  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_3

  integer :: i_spin,  i_spin_2
  integer :: i_index, i_k_point
  integer :: i_my_batch
  integer :: info
  logical :: is_restarted, is_extrapolated
  character*120 :: info_str

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: hartree_partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_gradient(:,:,:)
  real*8, pointer :: kinetic_density(:,:)
  real*8, pointer :: delta_rho_KS(:,:)
  real*8, pointer :: delta_rho_gradient(:,:,:)
  real*8, pointer :: delta_kinetic_density(:,:)
  integer :: my_full_points 

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8, allocatable :: batch_times(:)
  real*8 time_start
  real*8 time0, time_work, time_all, time_include_irr_to_full, time1

  ! Timers
  character(*), parameter :: deffmt = "10X"


  real*8       :: time_total_local 
  real*8       :: clock_total_local
  real*8       :: time_post 
  real*8       :: clock_post

   
  logical :: gpu_save

  logical :: updating_non_partition_point, point_on_atom
  
  !shanghui debug 
  !( this debug can be used together with DFPT_phonon/integrate_first_order_rho_p1.f90's debug_one_point
  !  to give all matrix elements' contribution to rho(r) at one point) 
  logical, parameter :: debug_one_point = .false.
  real*8, dimension(n_species) :: r_grid_min_sq
  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  gpu_save = use_gpu

  if (use_gpu_density .and. .not. use_gpu) then
    ! A failsafe.  Only works for myid=0, though.
    write(info_str,'(2X,A)')
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "You have request GPU acceleration of the charge density update, but no GPU"
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    write(info_str,'(2X,A)') "is available.  Turning off GPU acceleration for the charge density update."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    use_gpu = .false.
  end if

  use_gpu = use_gpu .and. use_gpu_density

  time_total_local = 0d0;
  clock_total_local = 0d0;
  
  call get_timestamps(time_total_local, clock_total_local)
  
  call localorb_info( &
       "Evaluating new KS density using the density matrix", &
       use_unit,'(2X,A)', OL_norm )

  if (use_gpu) then
    write(info_str,'(2X,A)') "GPU acceleration will be used when evaluting the KS density using the density matrix."
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    gpu_density_used = .true.
  else
    gpu_density_always_used = .false.
    gpu_density_used = .false.
  end if

  if(use_small_component)then
     allocate(temp_rho_small(n_max_batch_size,n_spin),stat=info)
     call check_allocation(info, 'temp_rho_small                ')
  end if

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)


  n_bp = use_batch_permutation

  if(n_bp > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab
    my_full_points = batch_perm(n_bp)%n_full_points

    allocate(hartree_partition_tab(batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,1,hartree_partition_tab_std,hartree_partition_tab)

    allocate(rho(n_spin,batch_perm(n_bp)%n_full_points))
    call permute_point_array(n_bp,n_spin,rho_std,rho)

    if (use_density_gradient) then
      allocate(rho_gradient(3,n_spin,batch_perm(n_bp)%n_full_points))
      call permute_point_array(n_bp,3*n_spin,rho_gradient_std,rho_gradient)
    else
      allocate(rho_gradient(1,1,1))
    endif

    if (use_meta_gga) then
      allocate(kinetic_density(n_spin,batch_perm(n_bp)%n_full_points))
      call permute_point_array(n_bp,n_spin,kinetic_density_std,kinetic_density)
    else
      allocate(kinetic_density(1,1))
    endif

    allocate(delta_rho_KS(batch_perm(n_bp)%n_full_points,n_spin))

    if (use_density_gradient) then
      allocate(delta_rho_gradient(3,batch_perm(n_bp)%n_full_points,n_spin))
    else
      allocate(delta_rho_gradient(1,1,1))
    end if

    if (use_meta_gga) then
      allocate(delta_kinetic_density(batch_perm(n_bp)%n_full_points,n_spin))
    else
      allocate(delta_kinetic_density(1,1))
    endif

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

    ! evaluate_densmat uses the BLACS-distributed Hamiltonian matrix as a work
    ! array to calculate the density matrix, so we need to reinitialize
    ! index arrays for local_index communication
    call init_comm_full_local_matrix( &
         batch_perm(n_bp)%n_basis_local, &
         batch_perm(n_bp)%i_basis_local )

  elseif (use_rlsy_symmetry) then
    n_my_batches_work = batch_irr%n_my_batches
    batches_work => batch_irr%batches
    partition_tab => batch_irr%partition_tab
    hartree_partition_tab => batch_irr%partition_tab
    !call compute_batch_rlsy_full_to_irr(hartree_potential_std, rho_std, &
    !                         rho_gradient_std, kinetic_density_std)
    my_full_points = batch_irr%n_full_points

    allocate(rho(n_spin,batch_irr%n_full_points))
    do i_spin = 1, n_spin
      call compute_batch_rlsy_full_to_irr_points(rho_std(i_spin,:),rho(i_spin,:))
    end do

    if (use_density_gradient) then
      allocate(rho_gradient(3,n_spin,batch_irr%n_full_points))
      do i_spin = 1, n_spin
        call compute_batch_rlsy_full_to_irr_points_gradient(rho_gradient_std(:,i_spin,:),rho_gradient(:,i_spin,:))
      end do
    else
      allocate(rho_gradient(1,1,1))
    endif

    if (use_meta_gga) then
      allocate(kinetic_density(n_spin,batch_irr%n_full_points))
      do i_spin = 1, n_spin
        call compute_batch_rlsy_full_to_irr_points(kinetic_density_std(i_spin,:),kinetic_density(i_spin,:))
      end do
    else
      allocate(kinetic_density(1,1))
    endif

    allocate(delta_rho_KS(batch_irr%n_full_points,n_spin))
    if (use_density_gradient) then
      allocate(delta_rho_gradient(3,batch_irr%n_full_points,n_spin))
    else
      allocate(delta_rho_gradient(1,1,1))
    end if
    if (use_meta_gga) then
      allocate(delta_kinetic_density(batch_irr%n_full_points,n_spin))
    else
      allocate(delta_kinetic_density(1,1))
    endif


  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    hartree_partition_tab => hartree_partition_tab_std
    rho => rho_std
    rho_gradient => rho_gradient_std
    kinetic_density => kinetic_density_std
    delta_rho_KS => delta_rho_KS_std
    delta_rho_gradient => delta_rho_gradient_std
    delta_kinetic_density => delta_kinetic_density_std
    my_full_points = n_full_points

  endif

  ! BL: Debug
  !if (output_level .eq. "full") then 
  !   write (info_str,'(A,I5,A,I5)') "  | CPU ", myid, &
  !          ": Number of batches = " , n_my_batches_work
  !   call localorb_allinfo(info_str)
  !endif

  allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize charge density convergence criterion
  rho_change = 0.d0

  delta_rho_KS = 0.d0

  if (use_density_gradient) then
     delta_rho_gradient = 0.0d0
  end if

  if (use_meta_gga) then
     delta_kinetic_density = 0.0d0
  endif

  if(packed_matrix_format == PM_none)then
     if(.not. allocated(density_matrix))then
        call aims_allocate( density_matrix, n_centers_basis_T, n_centers_basis_T, "density_matrix" )
     end if
  else
     ! Allocate dummy since otherways compiling with checking doesn't work
    call aims_allocate(density_matrix, 1,1,                                      "density_matrix" )
  end if

  if (use_wf_extrapolation .and. wf_extra_use_densmat) then

     ! Please note:
     ! Load balancing doesn't work with wf_extrapolation because:
     ! 1. density_matrix_sparse is stored full in reality
     ! 2. the content of density_matrix_sparse is different because of batch permutation
     ! 3. the batch permutations currently change after 1st SCF cycle
     !
     ! BTW: Currently wf_extrapolate_density_matrix() doesn't work with Scalapack at all ...
     !
     if(use_load_balancing) call aims_stop('ERROR: Load balancing currently doesn''t work with wf_extrapolation for densmat')

     if (wf_extrapolation_has_been_done()) then
        is_extrapolated = .false.
     else
        if(packed_matrix_format == PM_none)then
           call aims_allocate( density_matrix_extrapolated, n_centers_basis_T, n_centers_basis_T, n_spin,    "density_matrix_extrapolated" )
        else
           call aims_allocate( density_matrix_sparse_extrapolated, n_hamiltonian_matrix_size, n_spin, "density_matrix_sparse_extrapolated" )
        end if

        ! JW: Unfortunately, I need to calculate this here because
        ! it is implemented for both spins simultaneously.  This is also
        ! slightly more efficient because the Cholesky decomposition of the
        ! overlap has to be done only once.
        ! Additionally, I hope that MD runs tend to be less memory critical
        ! than other operation modes.
        call wf_extrapolate_density_matrix(density_matrix_extrapolated, &
        &                                  density_matrix_sparse_extrapolated, &
        &                                  .false., is_extrapolated)
     end if
  else
     is_extrapolated = .false.
  end if

  ! --- start actual calculation

  time_work = 0
  time_all  = 0

  if (use_gpu) then
    ! TODO:  Update this block to use output_mem_array_gpu in aims_gpu module
    write (info_str,'(A)') &
         'Reporting GPU memory usage by myid 0.  Note that it is likely not the only MPI task binding to its GPU!'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_full_points)*dble(n_spin)*8/1.d6,' MB on GPU for rho'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_full_points)*dble(n_spin)*8/1.d6,' MB on GPU for delta_rho_KS'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_max_batch_size)*8/1.d6,' MB on GPU for temp_rho'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_full_points)*8/1.d6,' MB on GPU for partition_tab'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_full_points)*8/1.d6,' MB on GPU for hartree_partition_tab'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_max_compute_dens)*dble(n_max_batch_size)*8/1.d6,' MB on GPU for wave'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_max_compute_dens)*dble(n_max_compute_dens)*8/1.d6,' MB on GPU for density_matrix'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_max_batch_size)*4/1.d6,' MB on GPU for dev_batchPoint2FullPoint'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(n_max_compute_dens)*dble(3*n_max_batch_size)*8/1.d6,' MB on GPU for dev_work'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    write (info_str,'(A,F12.3,A,A)') &
         'Allocating ',dble(3*n_max_batch_size)*dble(n_max_batch_size)*8/1.d6,' MB on GPU for dev_resultMat'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    if (use_density_gradient) then
      write (info_str,'(A,F12.3,A,A)') &
           'Allocating ',dble(3*n_full_points)*dble(n_spin)*8/1.d6,' MB on GPU for rho_gradient'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') &
           'Allocating ',dble(3*n_full_points)*dble(n_spin)*8/1.d6,' MB on GPU for delta_rho_gradient'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') &
           'Allocating ',dble(3*n_max_batch_size)*8/1.d6,' MB on GPU for temp_rho_gradient'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
      write (info_str,'(A,F12.3,A,A)') &
           'Allocating ',dble(3*n_max_compute_dens)*dble(n_max_batch_size)*8/1.d6,' MB on GPU for gradient_basis_wave'
      call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    end if

    ! AJL for WPH: Todo: extend for meta-gga. I think you need just the kinetic_density, temp_ and delta_ arrays

    ! Note: For system with large vacuum regions, n_max_compute_den can be zero
    !       on certain tasks!  This means that the GPU arrays for wave,
    !       density_matrix, and gradient_basis_wave may be unallocated.
    call density_create_gpu( &
          n_max_compute_dens, &
          n_max_batch_size, &
          my_full_points, &
          n_spin, &
          fort_log_to_c_int(use_density_gradient))
    call set_delta_rho_gpu(delta_rho_KS, my_full_points, n_spin)
    call set_rho_gpu(rho, my_full_points, n_spin)
    call set_rho_change_gpu(rho_change, n_spin)
    call set_partition_tab_gpu(partition_tab, my_full_points)
    call set_hartree_partition_tab_gpu(hartree_partition_tab, my_full_points)
    if (use_density_gradient) then
      call set_delta_rho_gradient_gpu(delta_rho_gradient, my_full_points, n_spin)
      call set_rho_gradient_gpu(rho_gradient, my_full_points, n_spin)
    end if
    ! AJL for WPH: Obviously needs extending for kinetic density
  end if

  do i_spin = 1, n_spin

     ! --- density matrix construction/administration

     is_restarted = .false.
     if (use_scalapack .and. restart_read .and. .not.(use_local_index.or.force_single_restartfile &
        .or.restart_eigenvectors_periodic)) then
        call restart_scalapack_read(density_matrix_sparse, i_spin, is_restarted)
     elseif (fo_density) then
        call restart_scalapack_read(density_matrix_sparse, i_spin, is_restarted)
     else if (use_scalapack .and. restart_read .and. (force_single_restartfile.or. &
        restart_eigenvectors_periodic) .and..not.use_local_index) then
        call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
           &                     density_matrix, density_matrix_sparse, i_spin, .false.)
     is_restarted = .true.
     end if

     if (.not. is_restarted) then
        if (is_extrapolated) then
           write(info_str, "('  Using extrapolated density matrix')")
           call localorb_info(info_str, use_unit, '(A)', OL_norm)
           if(packed_matrix_format == PM_none) then
              density_matrix = density_matrix_extrapolated(:,:, i_spin)
           else
              density_matrix_sparse(1:n_hamiltonian_matrix_size) = density_matrix_sparse_extrapolated(:, i_spin)
           end if
        else

! BB: experimental symmetry reconstruction of density. Will be cleaned up when
!     operational and fully tested. Should do no harm unless keyword:
!     "symmetry_reduced_k_grid_spg .true." is set.
           if (use_symmetry_reduced_spg) then
              call evaluate_densmat_sym(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
                      density_matrix, density_matrix_sparse, i_spin, .false.)
           else
              do_elsi_rw = .true.
              call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
                      density_matrix, density_matrix_sparse, i_spin, .false.)
              do_elsi_rw = .false.
           end if
        end if
     end if


     if (use_scalapack .and. restart_write .and. .not.(use_local_index.or.force_single_restartfile &
        .or.restart_eigenvectors_periodic)) then
        call restart_scalapack_write(density_matrix_sparse, i_spin)
     end if

     ! ---------------- end construct density matrix -------------------------

     if(.not. allocated(dist_tab))then
        allocate(dist_tab(n_centers_integrals),stat=info)
        call check_allocation(info, 'dist_tab                      ')
     end if
     if(.not. allocated(dist_tab_sq))  then
        allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dist_tab_sq                   ')
     end if
     if(.not. allocated(i_r))then
        allocate(i_r(n_max_compute_atoms),stat=info)
        call check_allocation(info, 'i_r                           ')
     end if
     if(.not. allocated(dir_tab)) then
        allocate(dir_tab(3, n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dir_tab                       ')
     end if
     if(.not. allocated(trigonom_tab)) then
        allocate(trigonom_tab(4, n_max_compute_atoms),stat=info)
        call check_allocation(info, 'trigonom_tab                  ')
     end if
     if(.not. allocated(density_matrix_con))then
        call aims_allocate( density_matrix_con, n_max_compute_dens*n_max_compute_dens, "density_matrix_con" )
     end if
     if(.not. allocated(work))then
        call aims_allocate( work, n_max_compute_dens, n_max_batch_size,                              "work" )
     end if
     if(.not. allocated(grad_work))then
        call aims_allocate( grad_work, n_max_compute_dens, n_max_batch_size,                    "grad_work" )
     end if
     if(.not. allocated(i_basis))then
        allocate(i_basis(n_centers_basis_T),stat=info)
        call check_allocation(info, 'i_basis                       ')
     end if
     if(.not.allocated(radial_wave))then
        allocate(radial_wave(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave                   ')
     end if
     if(.not. allocated(radial_wave_deriv))then
        allocate(radial_wave_deriv(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave_deriv             ')
     end if
     if(.not. allocated(radial_wave_2nd_deriv))then
        allocate(radial_wave_2nd_deriv(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'radial_wave_2nd_deriv         ')
     end if
     if(.not. allocated(wave))then
        call aims_allocate( wave, n_max_compute_dens, n_max_batch_size,                             "wave" )
     end if
     if(.not. allocated(kinetic_wave))then
        allocate(kinetic_wave(n_max_compute_fns_dens),stat=info)
        call check_allocation(info, 'kinetic_wave                  ')
     end if

     l_ylm_max = l_wave_max

     if(.not. allocated( ylm_tab))then
        allocate( ylm_tab( (l_ylm_max+1)**2,n_max_compute_atoms),stat=info )
        call check_allocation(info, 'ylm_tab                       ')
     end if
     if(.not. allocated( index_lm))then
        allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info)
        call check_allocation(info, 'index_lm                      ')
      end if


     if (use_density_gradient) then   !     allocate local arrays needed for gradients
        if(.not. allocated (gradient_basis_wave))then
           call aims_allocate( gradient_basis_wave, n_max_compute_dens, 3, n_max_batch_size, "gradient_basis_wave" )
        end if
        if(.not. allocated( dylm_dtheta_tab))then
           allocate( dylm_dtheta_tab((l_ylm_max+1)**2, n_max_compute_atoms),stat=info )
           call check_allocation(info, ' dylm_dtheta_tab              ')
         end if
        if(.not. allocated( scaled_dylm_dphi_tab))then
           allocate( scaled_dylm_dphi_tab((l_ylm_max+1)**2,n_max_compute_atoms),stat=info )
           call check_allocation(info, 'scaled_dylm_dphi_tab          ')
          end if
     end if


     ! ---------- calculate LDA/GGA+U occupation numbers from the density matrix -----------
     !
     ! this is totally independent of the force calculation, but
     ! it does need the density matrix. Which is only available at this point
     if(.not.use_local_index) then
        ! Please note if occ_numbers_plus_u should be implemented for use_local_index:
        ! The size of density_matrix_sparse may be different from n_hamiltonian_matrix_size
        ! if loadbalancing is in effect!
        if(use_plus_u)then

           ! take the occupation matrix from a file 
           if (plus_u_occupation_matrix_control_read .eqv. .true.) then
               call plus_u_read_occ_mat()
           else ! calculate it
               if (plus_u_mulliken_charges .eqv. .true.) then
                  call occ_numbers_plus_u_mulliken(density_matrix_sparse, density_matrix, i_spin)
               elseif (plus_u_full .eqv. .true.) then
                  call occ_numbers_plus_u_full(density_matrix_sparse, density_matrix, i_spin)
               else ! use on-site representation (default)
                  call occ_numbers_plus_u(density_matrix_sparse, density_matrix, i_spin)
               endif
           endif

           ! check if everything is correct with the occupation matrix
           call plus_u_check_occupation()

        endif
     endif

     ! --------------------------------------------------------------------------------------

     ! initialize index_lm
     i_index = 0
     do i_l = 0, l_ylm_max, 1
        do i_m = -i_l, i_l
           i_index = i_index + 1
           index_lm(i_m, i_l) = i_index
        enddo
     enddo

     call mpi_barrier(mpi_comm_global,info) ! Barrier is for correct timing!!!
     time0 = mpi_wtime()

     ! -------- go over the grids and construct the density using density matrix -------------


     i_basis_fns_inv = 0

     i_full_points = 0
     i_full_points_2 = 0
     i_full_points_3 = 0

     do i_my_batch = 1, n_my_batches_work, 1

           time_start = mpi_wtime()

           n_compute_c = 0
           i_basis = 0
           i_point = 0

           ! loop over one batch
           do i_index = 1, batches_work(i_my_batch)%size, 1

              i_full_points = i_full_points + 1
              
              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).gt.0.d0 .or. (atomic_MERM .and. mpb_solver_started)) then

                  i_point = i_point+1

                 !     get current integration point coordinate
                  coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

                  if(n_periodic>0)then
                     call map_to_center_cell( coord_current)
                  end if

                 ! compute atom-centered coordinates of current integration point,
                 ! as viewed from all atoms
                 
                  call tab_atom_centered_coords_p0 &
                      ( coord_current,  &
                      dist_tab_sq(1,i_point),  &
                      dir_tab(1,1,i_point), &
                      n_centers_integrals, centers_basis_integrals )


                 !     determine which basis functions are relevant at current integration point,
                 !     and tabulate their indices
                 
                    ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                  if (.not.prune_basis_once .or. (atomic_MERM .and. mpb_solver_started)) then
                      call prune_basis_p2 &
                      ( dist_tab_sq(1,i_point), &
                          n_compute_c, i_basis,  &
                          n_centers_basis_T, n_centers_basis_integrals, &
                          inv_centers_basis_integrals  )

                  end if
              end if
           enddo ! end loop over the angular shell


           if (prune_basis_once .and. .not.(atomic_MERM .and. mpb_solver_started)) then
              n_compute_c = batches_work(i_my_batch)%batch_n_compute
              i_basis(1:n_compute_c) = batches_work(i_my_batch)%batch_i_basis
           end if

           if(packed_matrix_format /= PM_none )then
              if(n_bp > 0) then
                do i=1,n_compute_c
                  ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
                enddo
                density_matrix_con(1:n_compute_c*n_compute_c) = 0
                do i=1,n_compute_c
                  i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                  do j=1,i
                    density_matrix_con(j+(i-1)*n_compute_c) = density_matrix_sparse(ins_idx(j)+i_off)
                  enddo
                enddo
              else
                call  prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
                     n_compute_c, i_basis)
              endif
           else
              call  prune_density_matrix(density_matrix, density_matrix_con, &
                   n_compute_c, i_basis)
           end if

           ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
           call collect_batch_centers_p2 &
           ( n_compute_c, i_basis, n_centers_basis_T, &
             n_centers_basis_integrals, inv_centers_basis_integrals, &
             n_batch_centers, batch_center &
           )

           n_points = i_point


           if (n_compute_c.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0
              do i_index = 1, batches_work(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1

                 updating_non_partition_point = .False.
                 if (max(partition_tab(i_full_points_2),&
                     hartree_partition_tab(i_full_points_2)).le.0.d0) then
                     updating_non_partition_point = .True.
                 end if                   
                 
                 if (max(partition_tab(i_full_points_2),&
                      hartree_partition_tab(i_full_points_2)).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

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
                        gradient_basis_wave (1:n_compute_c,1:3,i_point)=0d0

                     else                    

                        ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                        ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
                        ! without any copying and without doing any unnecessary operations.
                        ! The price is that the interface is no longer explicit in terms of physical
                        ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                        call prune_radial_basis_p2 &
                        ( n_max_compute_atoms, n_max_compute_fns_dens, &
                        dist_tab_sq(1,i_point), dist_tab, dir_tab(1,1,i_point), &
                        n_compute_atoms, atom_index, atom_index_inv, &
                        n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                        i_atom_fns, spline_array_start, spline_array_end, &
                        n_centers_basis_integrals, centers_basis_integrals, &
                        n_compute_c, i_basis, &
                        n_batch_centers, batch_center, &
                        one_over_dist_tab, rad_index, wave_index, l_index,&
                        l_count, &
                        fn_atom, n_zero_compute, zero_index_point &
                        )


                        ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                        ! for all atoms which are actually relevant
                        call tab_local_geometry_p2 &
                        ( n_compute_atoms, atom_index, &
                            dist_tab,  &
                            i_r &
                        )

                        ! Determine all needed radial functions from efficient splines

                        ! Now evaluate radial functions u(r) from the previously stored compressed
                        ! spline arrays

                        call evaluate_radial_functions_p0 &
                            (   spline_array_start, spline_array_end, &
                            n_compute_atoms, n_compute_fns,  &
                            dist_tab, i_r(1), &
                            atom_index, i_basis_fns_inv, &
                            basis_wave_ordered, radial_wave(1), &
                            .false., n_compute_c, n_max_compute_fns_dens   &
                            )
                            
                        ! for forces or density gradient, radial derivatives are required
                        if (use_density_gradient) then
                        call evaluate_radial_functions_p0 &
                                (spline_array_start, spline_array_end, &
                                n_compute_atoms, n_compute_fns,  &
                                dist_tab, i_r(1), &
                                atom_index, i_basis_fns_inv, &
                                basis_deriv_ordered,  &
                                radial_wave_deriv(1),  &
                                .true., n_compute_c, n_max_compute_fns_dens  &
                                )
                        end if


                        ! We compute the needed
                        ! ylm pieces separately from the cartesian evaluation 
                        ! that is done
                        ! for the hessians ...

                        ! compute trigonometric functions of spherical 
                        ! coordinate angles
                        ! of current integration point, viewed from all atoms
                        call tab_trigonom_p0 &
                                ( n_compute_atoms, dir_tab(1,1,i_point),  &
                                trigonom_tab(1,1) &
                                )

                        if (use_density_gradient) then
                            ! tabulate those ylms needed for gradients, 
                            ! i.e. ylm's for l_max+1
                        call tab_gradient_ylm_p0 &
                                (trigonom_tab(1,1), basis_l_max,  &
                                l_ylm_max, n_compute_atoms, atom_index, &
                                ylm_tab(1,1),  &
                                dylm_dtheta_tab(1,1),  &
                                scaled_dylm_dphi_tab(1,1) &
                                )

                        call evaluate_wave_gradient_p2  &
                            ( n_compute_c, n_compute_atoms, n_compute_fns, &
                                one_over_dist_tab, dir_tab(1,1,i_point), &
                                trigonom_tab(1,1),  &
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

                        ! tabulate total wave function value for each basis
                        ! function in all cases -
                        ! but only now are we sure that we have ylm_tab ...

                        ! tabulate total wave function value for each basis function

                        call evaluate_waves_p2  &
                        ( n_compute_c, n_compute_atoms, n_compute_fns, &
                            l_ylm_max, ylm_tab, one_over_dist_tab,   &
                            radial_wave, wave(1,i_point), &
                            rad_index, wave_index, l_index, l_count, fn_atom, &
                            n_zero_compute, zero_index_point &
                        )

                        if(use_small_component )then

                        call small_component(n_compute_c, i_basis, &
                                n_compute_atoms, dist_tab_sq(1,i_point), &
                                atom_index_inv(1), atom_index(1), &
                                gradient_basis_wave(1,1,i_point), &
                                temp_rho_small(i_point,i_spin), &
                                n_max_compute_dens, i_spin)
                        end if

                        ! reset i_basis_fns_inv
                        i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

                     end if ! point_on_atom
                 end if ! end if (partition_tab.gt.0)
              enddo ! end loop over a batch
              ! - all quantities that are evaluated pointwise are now known ...
 
              if (n_points.gt.0) then
                 ! Now perform all operations which are done across the 
                 ! entire integration shell at once.

                 ! New density is always evaluated

                 if (use_gpu) then
                    call set_density_matrix_gpu(density_matrix_con, n_compute_c)
                    call evaluate_KS_density_densmat_gpu ( &
                      n_points, wave, n_compute_c, n_max_compute_dens)
                 else
                    call evaluate_KS_density_densmat &
                      (  n_points, wave(1,1), n_compute_c,   &
                      temp_rho(1,i_spin), n_max_compute_dens, &
                      n_centers_basis_T, density_matrix_con, work(1,1) )
                 end if

                 if(use_small_component )then
                    temp_rho = temp_rho + temp_rho_small
                    temp_rho_small = 0.d0
                    if (use_gpu) then
                       call aims_stop("*** TODO: SMALL COMPONENT HAS TO BE IMPLEMENTED FOR CUDA",&
                            "update_density_densmat.f90")
                    end if
                 end if
                 
                 ! Finally, get the density gradient ...
                 if (use_density_gradient) then
                    if (use_gpu) then
                       call eval_density_grad_densmat_gpu ( &
                         n_points, gradient_basis_wave, &
                         n_compute_c, n_max_compute_dens)
                    else
                       call evaluate_density_gradient_denmat  &
                         (n_points, gradient_basis_wave, n_compute_c,  &
                         density_matrix_con, &
                         temp_rho_gradient(1,1,i_spin),  wave(1,1),  &
                         n_max_compute_dens, work, grad_work)
                    end if

                    if (use_meta_gga) then
                      if (use_gpu) then
                        ! AJL for WPH: needs an equivalent call to evaluate kinetic density
                      else
                        call evaluate_kinetic_density_denmat  &
                             (n_points, gradient_basis_wave, n_compute_c,  &
                             density_matrix_con,  temp_kinetic_density(1,i_spin),  &
                             n_max_compute_dens)
                      end if
                    end if ! use_meta_gga

                 end if ! use_density_gradient
              end if ! (n_points>0)


              ! calculate change in electron density

              if (use_gpu) then
                 call update_delta_rho_KS_gpu( &
                       i_full_points_3, &
                       batches_work(i_my_batch)%size, &
                       my_full_points, &
                       n_points, &
                       i_spin, &
                       n_spin)

                 if (use_density_gradient) then
                    call update_grad_delta_rho_KS_gpu( &
                          my_full_points, &
                          n_points, &
                          i_spin, &
                          n_spin)

                   if (use_meta_gga) then
                     ! AJL for WPH: Needs a corresponding call for delta_kinetic_density
                   end if ! use_meta_gga
                 end if ! use_density_gradient

                 i_full_points_3 = i_full_points_3 + &
                    batches_work(i_my_batch)%size

              else
                 i_point = 0

                 do i_index = 1, batches_work(i_my_batch)%size, 1

                    i_full_points_3 = i_full_points_3 + 1

                    if (max(partition_tab(i_full_points_3),&
                       hartree_partition_tab(i_full_points_3)).gt.0.or.(atomic_MERM .and. mpb_solver_started)) then

                       i_point = i_point + 1

                       ! and local change in charge density and gradient for 
                       ! Pulay mixing

                       if (spin_treatment.eq.0) then

                          delta_rho_KS(i_full_points_3, i_spin) =   &
                            delta_rho_KS(i_full_points_3, i_spin) &
                            + temp_rho(i_point,i_spin) - rho(i_spin, i_full_points_3)
                          ! (Rundong) PCE correction:
                          if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                            delta_rho_KS(i_full_points_3,i_spin) = delta_rho_KS(i_full_points_3,i_spin) + rho_pce(i_spin,i_full_points_3)
                          endif

                       elseif (spin_treatment.eq.1) then

                          ! Note that delta_rho_KS mixes the spin components of
                          ! the density - but, in the dm
                          ! formalism, we do first an _outermost_ loop over 
                          ! spin up, and only then over spin down.

                          ! AJL: there must be a way to shrink this if statement down with smart algebra?
                          if(i_spin == 1)then

                             delta_rho_KS(i_full_points_3, 1) &
                                = delta_rho_KS(i_full_points_3, 1)   &
                                + temp_rho(i_point, 1) &
                                - rho(1,i_full_points_3) 

                             delta_rho_KS(i_full_points_3, 2) &
                                = delta_rho_KS(i_full_points_3, 2)   &
                                + temp_rho(i_point, 1) &
                                - rho(1,i_full_points_3)

                             if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                               delta_rho_KS(i_full_points_3,1) = delta_rho_KS(i_full_points_3,1) + rho_pce(1,i_full_points_3)
                               delta_rho_KS(i_full_points_3,2) = delta_rho_KS(i_full_points_3,2) - rho_pce(2,i_full_points_3)
                             endif
                       
                          else

                             delta_rho_KS(i_full_points_3, 1) &
                                = delta_rho_KS(i_full_points_3, 1)   &
                                + temp_rho(i_point, 2) &
                                - rho(2,i_full_points_3)

                             delta_rho_KS(i_full_points_3, 2) &
                                = delta_rho_KS(i_full_points_3, 2)   &
                                - temp_rho(i_point, 2) &
                                + rho(2, i_full_points_3)

                          end if
                       end if

                       ! prepare charge density root-mean-square distance

                       if (i_spin == n_spin) then
                          do i_spin_2 = 1,n_spin

                             rho_change(i_spin_2) = rho_change(i_spin_2) +  &
                               partition_tab(i_full_points_3) *  &
                               delta_rho_KS(i_full_points_3, i_spin_2)**2
                          end do
                       end if

                       if (use_density_gradient) then
                          !prepare delta_rho_gradient for later mixing - this
                          !will be mixed in exactly the same way as delta_rho_KS
                          if (spin_treatment.eq.0) then

                             do i_coord = 1,3,1
                               delta_rho_gradient(i_coord, i_full_points_3, i_spin) &
                               = delta_rho_gradient(i_coord, i_full_points_3, i_spin) &
                               + temp_rho_gradient(i_coord, i_point, i_spin) &
                               - rho_gradient(i_coord, i_spin, i_full_points_3)
                             enddo
                             if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                               delta_rho_gradient(:,i_full_points_3,i_spin) = &
                               delta_rho_gradient(:,i_full_points_3,i_spin) + &
                               rho_pce_gradient(:,i_spin,i_full_points_3)
                             endif

                             if (use_meta_gga) then
                               delta_kinetic_density(i_full_points_3, 1) &
                               = delta_kinetic_density(i_full_points_3, 1) &
                               + temp_kinetic_density(i_point,1) &
                               - kinetic_density(1,i_full_points_3)
                             end if

                          elseif (spin_treatment.eq.1) then

                            ! AJL: shrink this to remove the if statement, as above?
                            if(i_spin == 1)then
                              do i_coord = 1,3,1

                                delta_rho_gradient(i_coord, i_full_points_3, 1) &
                                = delta_rho_gradient(i_coord, i_full_points_3, 1) &
                                + temp_rho_gradient(i_coord, i_point,1)  &
                                - rho_gradient(i_coord,1,i_full_points_3)

                                delta_rho_gradient(i_coord, i_full_points_3, 2) &
                                = delta_rho_gradient(i_coord, i_full_points_3, 2) &
                                + temp_rho_gradient(i_coord, i_point,1) &
                                - rho_gradient(i_coord,1,i_full_points_3)

                              enddo

                              if (use_meta_gga) then
                                delta_kinetic_density(i_full_points_3, 1) &
                                = delta_kinetic_density(i_full_points_3, 1) &
                                + temp_kinetic_density(i_point,1)  &
                                - kinetic_density(1,i_full_points_3)

                                delta_kinetic_density(i_full_points_3, 2) &
                                = delta_kinetic_density(i_full_points_3, 2) &
                                + temp_kinetic_density(i_point,1) &
                                - kinetic_density(1,i_full_points_3)
                              end if

                            else

                              do i_coord = 1,3,1

                                delta_rho_gradient(i_coord, i_full_points_3, 1) &
                                = delta_rho_gradient(i_coord, i_full_points_3, 1) &
                                + temp_rho_gradient(i_coord, i_point,2) &
                                - rho_gradient(i_coord,2,i_full_points_3)

                                delta_rho_gradient(i_coord, i_full_points_3, 2) &
                                = delta_rho_gradient(i_coord, i_full_points_3, 2) &
                                - temp_rho_gradient(i_coord, i_point,2) &
                                + rho_gradient(i_coord,2,i_full_points_3)

                              enddo

                              if (use_meta_gga) then

                                 delta_kinetic_density(i_full_points_3, 1) &
                                 = delta_kinetic_density(i_full_points_3, 1) &
                                 + temp_kinetic_density(i_point,2) &
                                 - kinetic_density(2,i_full_points_3)

                                 delta_kinetic_density(i_full_points_3, 2) &
                                 = delta_kinetic_density(i_full_points_3, 2) &
                                 - temp_kinetic_density(i_point,2) &
                                 + kinetic_density(2,i_full_points_3)                                      

                              end if ! use_meta_gga

                            end if ! i_spin

                          end if ! spin_treatment
                      end if ! use_density_gradient
                   endif
                end do

                end if ! cublas
           else 

                ! Even if n_compute .eq. 0 for the entire current batch of grid points,
                ! we still need to
                ! make sure that the density _change_ at this point is 
                ! ( zero minus previous density ).
                ! This ensures that. even for a zero KS density, a potentially non-zero
                ! initialization density
                ! is subtracted to properly account for the density change ....

                if (use_gpu) then
                   ! In CUDA we initialize the density with zero and perform 
                   ! the neccessary operations exactly as for n_compute > 0
                   call evaluate_KS_density_densmat_gpu ( &
                      batches_work(i_my_batch)%size, wave, 0, n_max_compute_dens)

                   call update_delta_rho_KS_gpu( &
                       i_full_points_3, &
                       batches_work(i_my_batch)%size, &
                       my_full_points, &
                       batches_work(i_my_batch)%size, &
                       i_spin, &
                       n_spin)

                 if (use_density_gradient) then
                    call eval_density_grad_densmat_gpu (&
                          batches_work(i_my_batch)%size, &
                          gradient_basis_wave, 0, n_max_compute_dens)

                    call update_grad_delta_rho_KS_gpu( &
                          my_full_points, &
                          batches_work(i_my_batch)%size, &
                          i_spin, &
                          n_spin)

                    if (use_meta_gga) then
                      ! AJL for WPH: needs equivalent routines for mgga
                    end if

                 end if

                    i_full_points_2 = i_full_points_2 + &
                         batches_work(i_my_batch)%size
                    i_full_points_3 = i_full_points_3 + &
                         batches_work(i_my_batch)%size

                else ! not use_gpu

                    do i_index = 1, batches_work(i_my_batch)%size, 1

                       i_full_points_2 = i_full_points_2 + 1
                       i_full_points_3 = i_full_points_3 + 1

                       if (max(partition_tab(i_full_points_3),&
                            hartree_partition_tab(i_full_points_3)).gt.0.d0.or.(atomic_MERM .and. mpb_solver_started)) then

                          ! only evaluate once both spin components are known
                          if (i_spin.eq.n_spin) then

                             ! local change in charge density and gradient for Pulay mixing

                             if (spin_treatment.eq.0) then

                                do i_spin_2 = 1, n_spin, 1
                                   delta_rho_KS(i_full_points_3, i_spin_2) =   &
                                        delta_rho_KS(i_full_points_3, i_spin_2) &
                                        - rho(i_spin_2, i_full_points_3)

                                   ! (Rundong) PCE correction:
                                   if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                                     delta_rho_KS(i_full_points_3,i_spin_2) = &
                                     delta_rho_KS(i_full_points_3,i_spin_2) + &
                                     rho_pce(i_spin_2,i_full_points_3)
                                   endif
                                enddo

                             elseif (spin_treatment.eq.1) then

                                delta_rho_KS(i_full_points_3, 1) =   &
                                     delta_rho_KS(i_full_points_3, 1) -  &
                                     ( rho(1,i_full_points_3) + rho(2,i_full_points_3) )

                                delta_rho_KS(i_full_points_3, 2) =   &
                                     delta_rho_KS(i_full_points_3, 2) -  &
                                     ( rho(1,i_full_points_3) - rho(2, i_full_points_3) )

                                !(Rundong) I'am confused here. Check later.
                                if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                                  delta_rho_KS(i_full_points_3,1) = delta_rho_KS(i_full_points_3,1) + &
                                     rho_pce(1,i_full_points_3) + rho_pce(2,i_full_points_3)
                                endif

                             end if

                             ! prepare charge density root-mean-square distance
                             do i_spin_2 = 1, n_spin, 1

                                rho_change(i_spin_2) = rho_change(i_spin_2) +  &
                                     partition_tab(i_full_points_3) *  &
                                     delta_rho_KS(i_full_points_3, i_spin_2)**2

                             enddo

                             if (use_density_gradient) then
                                !prepare delta_rho_gradient for later mixing - this
                                !will be mixed in exactly the same way as delta_rho_KS
                                if (spin_treatment.eq.0) then

                                   do i_spin_2 = 1, n_spin, 1
                                     do i_coord = 1,3,1
                                        delta_rho_gradient(i_coord, i_full_points_3, i_spin_2) =   &
                                           delta_rho_gradient(i_coord, i_full_points_3, i_spin_2) - &
                                           rho_gradient(i_coord,i_spin_2,i_full_points_3)
                                     enddo
                                     if(flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks)then
                                       delta_rho_gradient(:,i_full_points_3,i_spin_2) = &
                                       delta_rho_gradient(:,i_full_points_3,i_spin_2) + &
                                       rho_pce_gradient(:,i_spin_2,i_full_points_3)
                                     endif
                                   enddo

                                   if (use_meta_gga) then
                                      delta_kinetic_density(i_full_points_3, 1) =   &
                                           delta_kinetic_density(i_full_points_3, 1) - &
                                           kinetic_density(1,i_full_points_3)                                           
                                   end if

                                elseif (spin_treatment.eq.1) then

                                   do i_coord = 1,3,1

                                      delta_rho_gradient(i_coord, i_full_points_3, 1) = &
                                           delta_rho_gradient(i_coord, i_full_points_3, 1) - &
                                           (rho_gradient(i_coord,1,i_full_points_3)+ &
                                           rho_gradient(i_coord,2,i_full_points_3))

                                      delta_rho_gradient(i_coord, i_full_points_3, 2) = &
                                           delta_rho_gradient(i_coord, i_full_points_3, 2) - &
                                           (rho_gradient(i_coord,1,i_full_points_3)- &
                                           rho_gradient(i_coord,2,i_full_points_3))

                                   enddo

                                   if (use_meta_gga) then

                                      delta_kinetic_density(i_full_points_3, 1) = &
                                           delta_kinetic_density(i_full_points_3, 1) - &
                                           (kinetic_density(1,i_full_points_3)+ &
                                           kinetic_density(2,i_full_points_3))

                                      delta_kinetic_density(i_full_points_3, 2) = &
                                           delta_kinetic_density(i_full_points_3, 2) - &
                                           (kinetic_density(1,i_full_points_3)- &
                                           kinetic_density(2,i_full_points_3))                                           

                                   end if ! use_meta_gga

                                end if ! spin_treatment
                             end if ! use_density_gradient

                          end if ! (i_spin.eq.n_spin)

                       endif !(partition_tab)

                    enddo ! (loop over all grid points, for empty batches)
                end if ! cublas

           end if  ! end if (n_compute.gt.0)

           ! end if ! end distribution over threads

        batch_times(i_my_batch) = mpi_wtime() - time_start

     end do ! end loop over grid batches

     ! Get work time and total time after barrier
     time_work = time_work + mpi_wtime()-time0
     call mpi_barrier(mpi_comm_global,info)
     time_all  = time_all  +mpi_wtime()-time0

     !--------------------------------- end go over the grids ------------------------------
  end do ! end i_spin

  call get_timestamps(time_post, clock_post)

  ! If we use a GPU, fetch the calculated data back in CPU memory
  if (use_gpu) then
    ! For density update, we always perform indexing on the GPU, so we always
    ! need to retrieve it here
    call calculate_rho_change_gpu(my_full_points, n_spin)
    call get_delta_rho_gpu(delta_rho_KS, my_full_points, n_spin)
    if (use_density_gradient) then
      call get_delta_rho_gradient_gpu(delta_rho_gradient, my_full_points, n_spin)
      ! AJL for WPH: presumably needs an equivalent for GPU
    end if
    call get_rho_change_gpu(rho_change, n_spin)

    write (info_str,'(A)') &
         'Deallocating memory used on GPU for charge density update via density matrices.'
    call localorb_info ( info_str,use_unit,'(2X,A)', OL_norm  )
    call density_destroy_gpu()
  end if
  use_gpu = gpu_save

  call get_times(time_post, clock_post)

  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') &
     '  Time summed over all CPUs for getting density from density matrix: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)
  

  ! broadcast the result to all threads
  call sync_density(rho_change)

  ! (Rundong) I don't distinguish spin channel currently:
  if((flag_rel.eq.REL_x2c.or.flag_rel.eq.REL_4c_dks).and.n_spin.eq.2)rho_change(2)=0.d0

  do i_spin = 1, n_spin, 1
     rho_change(i_spin) = sqrt(rho_change(i_spin))
  enddo

  if(n_bp > 0) then

    do i_spin = 1, n_spin
      call permute_point_array_back(n_bp,1,delta_rho_KS(:,i_spin),delta_rho_KS_std(:,i_spin))
      if(use_density_gradient) then
        call permute_point_array_back(n_bp,3,delta_rho_gradient(:,:,i_spin),delta_rho_gradient_std(:,:,i_spin))
        if (use_meta_gga) then
          call permute_point_array_back(n_bp,1,delta_kinetic_density(:,i_spin),delta_kinetic_density_std(:,i_spin))
        end if
      endif
    enddo

    deallocate(hartree_partition_tab)
    deallocate(rho)
    deallocate(rho_gradient)
    deallocate(kinetic_density)

    deallocate(delta_rho_KS)
    deallocate(delta_rho_gradient)
    deallocate(delta_kinetic_density)

    deallocate(ins_idx)

  endif

  !wyj: TODO debug
  !print *, 'SCF_rho=', myid, delta_rho_KS_std(1:5,1), 'sum=', sum(delta_rho_KS_std(:,1))

  if(use_rlsy_symmetry) then

    time1 = mpi_wtime()

    rho_change = 0.0d0
    do i_spin = 1, n_spin
      call compute_batch_rlsy_irr_to_full_points(delta_rho_KS(:,i_spin), delta_rho_KS_std(:,i_spin))
      if(use_density_gradient) then
        call compute_batch_rlsy_irr_to_full_points_gradient(delta_rho_gradient(:,:,i_spin),delta_rho_gradient_std(:,:,i_spin))
        if (use_meta_gga) then
          call compute_batch_rlsy_irr_to_full_points(delta_kinetic_density(:,i_spin), delta_kinetic_density_std(:,i_spin))
        end if
      endif
      rho_change(i_spin) = dot_product(partition_tab_std(:),delta_rho_KS_std(:,i_spin)**2)
      call sync_density(rho_change)
      rho_change(i_spin) = sqrt(rho_change(i_spin))
    enddo


    deallocate(rho)
    deallocate(rho_gradient)
    deallocate(kinetic_density)

    deallocate(delta_rho_KS)
    deallocate(delta_rho_gradient)
    deallocate(delta_kinetic_density)

    time_include_irr_to_full = 0.0d0
    call mpi_barrier(mpi_comm_global,info)
    time_include_irr_to_full  = time_include_irr_to_full  +mpi_wtime()-time1
    call sync_real_number(time_include_irr_to_full)
!    write(info_str,'(a,f12.4,a)') "YY time for density irr to full ", &
!                                   time_include_irr_to_full, " s"
!    call localorb_info(info_str, use_unit, "(A)", OL_norm)

  endif

    time_include_irr_to_full = 0.0d0
    call mpi_barrier(mpi_comm_global,info)
    time_include_irr_to_full  = time_include_irr_to_full  +mpi_wtime()-time0
    call sync_real_number(time_include_irr_to_full)
!    write(info_str,'(a,f12.4,a)') "YY time for density update w/o density matrix update ", &
!                                   time_include_irr_to_full, " s"
!    call localorb_info(info_str, use_unit, "(A)", OL_norm)

  if (get_batch_weights) then
    call set_batch_weights(n_bp, batch_times)
  endif

  !---------- finally, deallocate stuff -----------------
  if(allocated (batch_times)) deallocate(batch_times)

  ! Allocatable arrays that are tracked
  if (allocated( density_matrix                     )) call aims_deallocate( density_matrix,           "density_matrix" )
  if (allocated( density_matrix_extrapolated        )) &
       call aims_deallocate( density_matrix_extrapolated,                                 "density_matrix_extrapolated" )
  if (allocated( density_matrix_sparse_extrapolated )) &
       call aims_deallocate( density_matrix_sparse_extrapolated,                   "density_matrix_sparse_extrapolated" )
  if (allocated( density_matrix_con                 )) call aims_deallocate( density_matrix_con,   "density_matrix_con" )
  if (allocated( work                               )) call aims_deallocate( work,                               "work" )
  if (allocated( grad_work                          )) call aims_deallocate( grad_work,                     "grad_work" )
  if (allocated( wave                               )) call aims_deallocate( wave,                               "wave" )
  if (allocated( gradient_basis_wave                )) call aims_deallocate( gradient_basis_wave, "gradient_basis_wave" )

  if (allocated( temp_rho_small       )) deallocate( temp_rho_small       )
  if (allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if (allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if (allocated( index_lm             )) deallocate( index_lm             )
  if (allocated( ylm_tab              )) deallocate( ylm_tab              )
  if (allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if (allocated( radial_wave_2nd_deriv)) deallocate( radial_wave_2nd_deriv)
  if (allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if (allocated( radial_wave          )) deallocate( radial_wave          )
  if (allocated( i_basis              )) deallocate( i_basis              )
  if (allocated( trigonom_tab         )) deallocate( trigonom_tab         )
  if (allocated( dir_tab              )) deallocate( dir_tab              )
  if (allocated( i_r                  )) deallocate( i_r                  )
  if (allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if (allocated( dist_tab             )) deallocate( dist_tab             )
  if (allocated( work_complex         )) deallocate( work_complex         )

  call get_times(time_total_local, clock_total_local)

   ! synchronize all cumulative CPU-time timestamps
   call sync_timing(time_total_local)

  if (output_level .eq. "full") then
     call output_timeheader(deffmt, 'Density and Forces (densmat)')
     call output_times(deffmt,'Total time', &
           time_total_local,clock_total_local)
     call output_times(deffmt,'Fetch Density Variables time', &
                   time_post,clock_post)
  endif
end subroutine update_density_densmat
!******
