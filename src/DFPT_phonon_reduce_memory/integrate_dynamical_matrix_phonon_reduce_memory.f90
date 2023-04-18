!****s* FHI-aims/integrate_dynamical_matrix_phonon_reduce_memory
!  NAME
!   integrate_dynamical_matrix_phonon_reduce_memory
!  SYNOPSIS

subroutine integrate_dynamical_matrix_phonon_reduce_memory &
     ( hartree_potential_std,first_order_potential_std, &  
       rho_std,rho_gradient_std, first_order_rho_std,  &
       partition_tab_std, basis_l_max,     &
       density_matrix_sparse, first_order_density_matrix_sparse,  & 
       energy_density_matrix_sparse, first_order_energy_density_matrix_sparse,  &  
       i_q_point, j_atom, j_coord, &
       hellman_feynman_dynamical_matrix_delta_part, &
       pulay_dynamical_matrix & 
     ) 

!  PURPOSE
! for pbc_pulay_hessian phonon_gamma
! shanghui 2013.12.30 @ Berlin

! for pbc_pulay_dynamical_matrix phonon_reduce_meomery
! shanghui 2015.08 

!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use synchronize_mpi
  use localorb_io
  use constants
  use species_data 
  use load_balancing
  use cartesian_ylm  !shanghui add for Hessian
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: hartree_potential_std
  complex*16, target, dimension(n_full_points) :: first_order_potential_std
  real*8, target, dimension(n_spin, n_full_points)    :: rho_std
  real*8, target, dimension(3, n_spin, n_full_points) :: rho_gradient_std
  complex*16, target, dimension(n_full_points) ::  first_order_rho_std
  real*8, target, dimension(n_full_points)            :: partition_tab_std

  integer ::  basis_l_max (n_species)

  integer, intent(in) ::  i_q_point
  integer, intent(in) ::  j_atom
  integer, intent(in) ::  j_coord

  real*8, dimension(n_hamiltonian_matrix_size),intent(in):: density_matrix_sparse
  real*8, dimension(n_hamiltonian_matrix_size),intent(in):: energy_density_matrix_sparse
  complex*16, dimension(n_hamiltonian_matrix_size),intent(in):: first_order_density_matrix_sparse
  complex*16, dimension(n_hamiltonian_matrix_size),intent(in):: first_order_energy_density_matrix_sparse
 
  complex*16, dimension(3, n_atoms, 3,n_atoms) :: hellman_feynman_dynamical_matrix_delta_part
  complex*16, dimension(3, n_atoms, 3,n_atoms) :: pulay_dynamical_matrix


!  INPUTS
!  o hartree_potential_std -- Hartree potential
!  o rho_std -- electron density
!  o rho_gradient_std -- gradient of electron density.
!    These should only ever be referenced if (use_gga)
!    import dimensions from above (if not used, all dimensions=1)
!  o partition_tab_std -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!
!  OUTPUT
!  o pulay_dynamical_matrix
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

  real*8,dimension(:,:,:),allocatable:: H_times_psi
  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave

  real*8, dimension(:),    allocatable :: en_density_xc
  real*8, dimension(n_spin) :: en_density_x
  real*8 :: en_density_c
  real*8, dimension(:, :), allocatable :: local_xc_derivs
  real*8, dimension(:, :), allocatable :: local_dVxc_drho
  real*8, dimension(:,:,:),allocatable :: xc_gradient_deriv

  real*8, dimension(:,:),  allocatable :: local_rho
  real*8, dimension(:,:,:),allocatable :: local_rho_gradient

  complex*16, dimension(:),  allocatable :: local_first_order_potential
  complex*16, dimension(:),  allocatable :: local_first_order_rho



  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)
  real*8 :: energy_partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave
  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave_npoints

!-------For pulay_hessian--------------------------------------

    real*8, dimension(:,:), allocatable   :: coords_npoints
    real*8, dimension(:,:,:), allocatable :: hessian_basis_wave
    real*8, dimension(:,:,:), allocatable :: sum_hessian
    real*8, dimension(:,:,:), allocatable :: cartesians
    real*8, dimension(:,:,:), allocatable :: sum_gradient
    real*8,dimension(:),allocatable:: radial_wave_2nd_deriv
    real*8, dimension(:,:,:,:), allocatable :: H_times_gradient_psi
    real*8,dimension(:)  ,allocatable:: kinetic_wave_deriv
    real*8, dimension(:,:,:), allocatable :: kinetic_gradient_basis_wave

    real*8, dimension(:),  allocatable :: density_matrix_con
    complex*16, dimension(:,:),  allocatable :: first_order_density_matrix_con
         
    real*8, dimension(:),  allocatable :: energy_density_matrix_con
    complex*16, dimension(:,:),  allocatable :: first_order_energy_density_matrix_con
!-------End for pulay_hessian----------------------------------





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

  integer :: n_compute_c, n_compute_a
!  integer :: i_basis(n_centers_basis_I)
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns

!  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
!  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
!  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

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

  integer i_basis_1, j_basis
  integer i_atom, i_atom_2 
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord, i_coord_1, i_coord_2
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2
   
  integer i_compute 
 
  integer i_counter
  integer index_hessian(3,3)

  integer :: i_spin
  character*200 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used


  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)
  real*8, pointer :: rho(:,:)
  real*8, pointer :: rho_gradient(:,:,:)
  real*8, pointer :: hartree_potential(:)
  complex*16, pointer :: first_order_rho(:)
  complex*16, pointer :: first_order_potential(:)


  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i_off, i, j, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  ! begin work




  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating Hessian matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating Hessian matrix: batch-based integration."
  endif
  call localorb_info(info_str, use_unit,'(A)',OL_norm)

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

  allocate( local_xc_derivs(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_xc_derivs               ')
  allocate( local_dVxc_drho(n_spin, n_max_batch_size),stat=info)
  call check_allocation(info, 'local_dVxc_drho            ')

  allocate( xc_gradient_deriv(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'xc_gradient_deriv             ')

  allocate( local_rho(n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho                     ')


  allocate( local_rho_gradient(3,n_spin,n_max_batch_size),stat=info)
  call check_allocation(info, 'local_rho_gradient            ')



!shanghui-gga: we do not gga, but we always need gradient. 
!  if ((flag_rel.eq.REL_none.or.flag_rel==REL_atomic_zora.or.flag_rel.eq.REL_own).and.(.not.(use_gga))) then
     !       no gradients needed
!     l_ylm_max = l_wave_max
!  else if ((flag_rel.eq.REL_zora).or.(use_gga).or.(flag_rel==REL_KOLNING_HARMON)) then
     l_ylm_max = l_wave_max
     allocate (gradient_basis_wave(n_max_compute_ham,3),STAT=info)
     call check_allocation(info, 'gradient_basis_wave           ')
     allocate (gradient_basis_wave_npoints(n_max_compute_ham,3,n_max_batch_size),STAT=info)
     call check_allocation(info, 'gradient_basis_wave_npoints        ')

     allocate (coords_npoints(3,n_max_batch_size),STAT=info)
     call check_allocation(info, 'coords_npoints        ')

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')

!  end if !-----shanghui always calcualet gradients-------------------------

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')


  allocate(H_times_psi(n_max_compute_ham, n_max_batch_size, n_spin), STAT=info )
  call check_allocation(info, 'H_times_psi                   ')

  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  allocate(wave(n_max_compute_ham, n_max_batch_size), STAT=info )
  call check_allocation(info, 'wave                          ')

  allocate( local_first_order_potential(n_max_batch_size),stat=info)
  call check_allocation(info, 'local_first_order_potential                     ')

  allocate( local_first_order_rho(n_max_batch_size),stat=info)
  call check_allocation(info, 'local_first_order_rho                     ')



  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')


!-------For pulay_hessian--------------------------------------
    if (.not.allocated(cartesians)) then
       allocate(cartesians(n_max_cartesian, 0:l_wave_max,n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'cartesians                    ')
    end if
    if (.not.allocated(sum_gradient)) then
       allocate(sum_gradient(3, (l_wave_max+1) ** 2, n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'sum_gradient                  ')
    end if
    if (.not.allocated(sum_hessian)) then
       allocate(sum_hessian(6, (l_wave_max+1) ** 2, n_max_compute_atoms ),stat=info)
       call check_allocation(info, 'sum_hessian                   ')
    end if
    if (.not.allocated(hessian_basis_wave)) then
       allocate(hessian_basis_wave(n_max_compute_ham, 6, n_max_batch_size),stat=info)
       call check_allocation(info, 'hessian_basis_wave            ')
    end if
    if(.not. allocated(radial_wave_2nd_deriv))then
       allocate(radial_wave_2nd_deriv(n_max_compute_fns_ham),stat=info)
       call check_allocation(info, 'radial_wave_2nd_deriv         ')
    end if
    if(.not.allocated(H_times_gradient_psi)) then
       allocate(H_times_gradient_psi(n_max_compute_ham,3, n_max_batch_size, n_spin), STAT=info )
       call check_allocation(info, 'H_times_gradient_psi                   ')
    endif
    allocate (kinetic_gradient_basis_wave(n_max_compute_ham,3,n_max_batch_size ),stat=info)
    call check_allocation(info, 'kinetic_gradient_basis_wave   ')

    allocate(kinetic_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave_deriv            ')

  if(.not. allocated(density_matrix_con))then
  allocate(density_matrix_con(n_max_compute_dens*n_max_compute_dens),stat=info)
  call check_allocation(info, 'density_matrix_con            ')
  end if

  if(.not. allocated(first_order_density_matrix_con))then
  allocate(first_order_density_matrix_con(n_max_compute_dens,n_max_compute_dens),stat=info)
  call check_allocation(info, 'first_order_density_matrix_con            ')
  end if

  if(.not. allocated(energy_density_matrix_con))then
  allocate(energy_density_matrix_con(n_max_compute_dens*n_max_compute_dens),stat=info)
  call check_allocation(info, 'energy_density_matrix_con            ')
  end if

  if(.not. allocated(first_order_energy_density_matrix_con))then
  allocate(first_order_energy_density_matrix_con(n_max_compute_dens,n_max_compute_dens),stat=info)
  call check_allocation(info, 'first_order_energy_density_matrix_con            ')
  end if

!-------End for pulay_hessian--------------------------------------




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
        allocate(zora_vector1(n_max_compute_ham,3,n_max_batch_size,n_spin),STAT=info )
        call check_allocation(info, 'zora_vector1                  ')
     end if
     if (.not.allocated(zora_vector2)) then
        allocate(zora_vector2(n_max_compute_ham,3,n_max_batch_size,n_spin),STAT=info )
        call check_allocation(info, 'zora_vector2                  ')
     end if

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

    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))


  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std
    rho => rho_std
    hartree_potential => hartree_potential_std
    rho_gradient => rho_gradient_std

    first_order_rho => first_order_rho_std
    first_order_potential => first_order_potential_std

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize
  !hellman_feynman_dynamical_matrix_delta_part(1:3,1:n_atoms,1:3,1:n_atoms) = (0.0d0, 0.0d0)
   pulay_dynamical_matrix(1:3,1:n_atoms,j_coord,j_atom) = (0.0d0, 0.0d0)



  i_basis_fns_inv = 0


  ! initialize index_lm

  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo


  !--------shanghui add for hessian_index-----------------
    i_counter = 0
    index_hessian = 0 
    do i_coord_1 = 1,3
       i_counter = i_counter +1
       index_hessian(i_coord_1, i_coord_1) = i_counter
       do i_coord_2 = i_coord_1+1,3
          i_counter = i_counter +1
          index_hessian(i_coord_1, i_coord_2) = i_counter
       end do
    end do

    do i_coord_1 = 1,3
       do i_coord_2 = 1,3
          index_hessian(i_coord_1, i_coord_2) = max(index_hessian(i_coord_1, i_coord_2), index_hessian(i_coord_2, i_coord_1))
       end do
    end do
  !--------shanghui end add for hessian_index--------------



  


  i_full_points = 0
  i_full_points_2 = 0


  ! perform partitioned integration, batch by batch of integration point.
  ! This will be the outermost loop, to save evaluations of the potential.
  ! and the Y_lm functions

  call mpi_barrier(mpi_comm_world,info) ! Barrier is for correct timing!!!
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

    !-----------shanghui begin test prune_matrix here-------------------------
     call  prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
                     n_compute_c, i_basis)  ! here only give upper matrix
         do i=1,n_compute_c
         do j=1,n_compute_c
            !if(i>j) density_matrix_con(j,i)= density_matrix_con(i,j)
            if(i>j) density_matrix_con((j-1)*n_compute_c+i)= density_matrix_con((i-1)*n_compute_c+j)
         enddo
         enddo

     call  prune_density_matrix_sparse_phonon_reduce_memory(first_order_density_matrix_sparse, &
                     first_order_density_matrix_con, &
                     n_compute_c, i_basis)

     call  prune_density_matrix_sparse(energy_density_matrix_sparse, energy_density_matrix_con, &
                     n_compute_c, i_basis)  ! here only give upper matrix
         do i=1,n_compute_c
         do j=1,n_compute_c
            !if(i>j) energy_density_matrix_con(j,i)= energy_density_matrix_con(i,j)
            if(i>j) energy_density_matrix_con((j-1)*n_compute_c+i)=  & 
                    energy_density_matrix_con((i-1)*n_compute_c+j)
         enddo
         enddo

     call  prune_density_matrix_sparse_phonon_reduce_memory(first_order_energy_density_matrix_sparse, &
                     first_order_energy_density_matrix_con, &
                     n_compute_c, i_basis)

     !-----------shanghui end test prune_matrix here-------------------------





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

        ! loop over one batch of integration points
        do i_index = 1, batches_work(i_my_batch)%size, 1

           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)!SAG
              coords_npoints(1:3,i_point) = batches_work(i_my_batch) % points(i_index) % coords(1:3)

              if (flag_rel.eq.REL_zora.or. flag_rel==REL_KOLNING_HARMON) then

                 call tab_global_geometry_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dir_tab(1,1,i_point), &
                      dist_tab_full, &
                      i_r_full, &
                      dir_tab_full_norm, &
                      n_centers_integrals,  centers_basis_integrals)

              end if

              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)
  

              energy_partition(i_point) = partition_tab(i_full_points)

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
  
              local_first_order_potential(i_point)=      &
              first_order_potential(i_full_points)


              local_first_order_rho(i_point)=      &
              first_order_rho(i_full_points)

              n_compute_atoms = 0
              n_compute_fns = 0

              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! Are stored in a compact spline array that can be accessed by spline_vector_waves,
              ! without any copying and without doing any unnecessary operations.
              ! The price is that the interface is no longer explicit in terms of physical
              ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
              !write(use_unit,*) dir_tab(:,1,i_point)
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


!shanghui-gga: we alwayse need gradients
!              if ((use_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON) ) then
                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )

!              else
!                ! tabulate distance and Ylm's w.r.t. other atoms
!                call tab_wave_ylm_p0 &
!                   ( n_compute_atoms, atom_index,  &
!                   trigonom_tab, basis_l_max,  &
!                   l_ylm_max, ylm_tab )
!              end if

              ! Now evaluate radial functions
              ! from the previously stored compressed spline arrays
              call evaluate_radial_functions_p0  &
                   ( spline_array_start, spline_array_end,  &
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
              if ((use_hartree_fock).and.(first_integration)) then
                 call evaluate_xc &
                      ( rho(1,i_full_points), &
                      rho_gradient(1,1,i_full_points), &
                      en_density_xc(i_point), local_xc_derivs(1,i_point), &
                      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), .true., coord_current &
                               )  !SAG
              else
                !---------We need local_dVxc_drho, but evaluate_xc do not have--
                !---------so we need evaluate_xc_shanghui below----------------
                !---------just like integrate_first_order_H-------------------- 
                ! call evaluate_xc  &
                !      ( rho(1,i_full_points),   &
                !      rho_gradient(1,1,i_full_points),  &
                !      en_density_xc(i_point), &
                !      en_density_x, en_density_c, &
                !      local_xc_derivs(1,i_point),  &
                !      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), &
                !      coord_current &
                !      )           !SAG

                 call evaluate_xc_shanghui  &
                      ( rho(1,i_full_points),   &
                      rho_gradient(1,1,i_full_points),  &
                      en_density_xc(i_point), &
                      en_density_x, en_density_c, &
                      local_xc_derivs(1,i_point),  &
                      xc_gradient_deriv(1,1,i_point), local_dVxc_drho(1,i_point), &
                      coord_current &
                      )  

              end if


              do i_spin = 1, n_spin, 1
                 local_potential_parts(i_spin) =   &
                     !--------V_Ze------------------------------
                     !-1.0d0/dsqrt(dist_tab_sq(1,i_point)) &
                    ! -1.0d0/dsqrt(dist_tab_sq(2,i_point))!+& 
                    !--------end VZe---------------------------
                    !---------V_ee-----------------------------
                      hartree_potential(i_full_points)   +   &
                    !  1.0d0/dsqrt(dist_tab_sq(1,i_point)) &
                    ! +1.0d0/dsqrt(dist_tab_sq(2,i_point))!+& 
                    !---------end V_ee-------------------------
                    !---------V_xc-----------------------------
                      local_xc_derivs(i_spin,i_point)
                    !---------end V_xc-------------------------


                 if (use_gga) then
                    sum_of_local_gradients(1:3,i_spin) =   &
                         xc_gradient_deriv(1:3,i_spin,i_point)*4.d0
                 else
                    sum_of_local_gradients(1:3,i_spin) = 0.d0
                 end if


              enddo

              ! Check whether relativistic corrections are needed at the present point.
              ! The check is based entirely on the local parts of the potential - i.e.
              ! in a GGA, the terms due to d(rho*exc)/d(|grad(rho|^2) is not evaluated.
              ! Hopefully this approximation to the full ZORA energy is small.
              if (flag_rel.eq.REL_zora.or. (flag_rel==REL_KOLNING_HARMON)) then

                 ! if we need ZORA, must get the _full_ local geometry in order to
                 ! create the superposition of atomic potentials which is used to estimate
                 ! the potential gradient for ZORA

                 call evaluate_pot_superpos_p0  &
                      (   &
                      i_r_full,   &
                      zora_potential_parts(1),  &
                      n_centers_integrals, centers_basis_integrals )

                 do i_spin = 1, n_spin, 1

                    ! factor 2.d0 required because a factor 1/2 is already included in kinetic_wave later ...
                    zora_operator(i_spin) =  &
                         2.d0 * light_speed_sq /  &
                         ( 2 * light_speed_sq -  &
                         zora_potential_parts(i_spin) )

                 enddo

              end if

!              if ((use_gga) .or. (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON)) then
!----------shanghui-gga: always calculate gradient---------------------
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

               if(.false.) then  !shanghui add for hessian_pulay, we always need hessian

                 ! and finally, assemble the actual gradients
                 call evaluate_wave_gradient_p2  &
                    ( n_compute_c, n_compute_atoms, n_compute_fns, &
                      one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                      l_ylm_max, ylm_tab,  &
                      dylm_dtheta_tab,  &
                      scaled_dylm_dphi_tab,  &
                      radial_wave,  &
                      radial_wave_deriv,  &
                      gradient_basis_wave_npoints(1:n_compute_c,1:3,i_point),  &
                      rad_index, wave_index, l_index, l_count, fn_atom, &
                      n_zero_compute, zero_index_point &
                    )
 
               else
 
                  call evaluate_radial_functions_deriv_p0 &
                      (spline_array_start, spline_array_end, &
                       n_compute_atoms, n_compute_fns,  &
                       dist_tab(1,i_point), i_r(1), &
                       atom_index, i_basis_fns_inv, &
                       basis_deriv_ordered,  &
                       radial_wave_2nd_deriv(1),  &
                       .true., n_max_compute_fns_dens &
                      )
        
                  call evaluate_cartesians &
                       (dir_tab(1,1,i_point), basis_l_max, l_wave_max,  &
                        atom_index, n_compute_atoms, cartesians(1,0,1) )
        
                 call evaluate_cartesian_hessian_and_gradient_terms & 
                       (basis_l_max, l_wave_max, cartesians(1,0,1), &
                        atom_index, n_compute_atoms,  &
                        ylm_tab(1,1),  &
                        sum_gradient(1,1,1),  &
                        sum_hessian(1,1,1))


                 call evaluate_wave_gradient_cartesian_p1 &
                     (dist_tab(1,i_point), i_r(1),  &
                     dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                     n_compute_c, i_basis, atom_index_inv,  &
                     i_basis_fns_inv, radial_wave(1),  &
                     radial_wave_deriv(1), &
                     ylm_tab(1,1),  &
                     sum_gradient(1,1,1),  &
                     gradient_basis_wave_npoints(1:n_compute_c,1:3,i_point) & 
                     ,n_compute_atoms)

                 call evaluate_wave_hessian_cartesian_p0 &
                       (dist_tab(1,i_point), i_r(1),  &
                       dir_tab(1,1,i_point), index_lm, l_wave_max,  &
                       n_compute_c, i_basis, atom_index_inv,  &
                       i_basis_fns_inv, radial_wave(1),  &
                       radial_wave_deriv(1),  &
                       radial_wave_2nd_deriv(1),  &
                       ylm_tab(1,1),  &
                       sum_gradient(1,1,1),  &
                       sum_hessian(1,1,1), &
                       hessian_basis_wave(1,1,i_point), n_compute_atoms)

               endif ! calc_hessian_pulay



!---------------shanghui calclate kinetic-gradient_phi---------------------
              ! First, obtain radial kinetic energy terms from vectorized splines
                  call evaluate_radial_functions_p0  &
                       ( spline_array_start, spline_array_end,  &
                       n_compute_atoms, n_compute_fns,   &
                       dist_tab(1,i_point), i_r,  &
                       atom_index, i_basis_fns_inv,  &
                       basis_kinetic_ordered, kinetic_wave(1),  &
                       .false., n_compute_c, n_max_compute_fns_ham )
     
                      !kinetic_wave =  kinetic_wave *0.5

                  call evaluate_radial_functions_deriv_p0 &
                       (spline_array_start, spline_array_end, &
                       n_compute_atoms, n_compute_fns,  &
                       dist_tab(1,i_point), i_r(1), &
                       atom_index, i_basis_fns_inv, &
                       basis_kinetic_ordered,  &
                       kinetic_wave_deriv(1),  &
                       .true., n_max_compute_fns_dens &
                       )

                      !kinetic_wave_deriv =  kinetic_wave_deriv *0.5
                  if(flag_rel == REL_atomic_zora)then   !!!DSB
                     kinetic_wave_deriv= kinetic_wave_deriv*0.5
                     kinetic_wave = kinetic_wave *0.5   !!!DSB preparation for atomic zora.
                  endif

                  call evaluate_wave_gradient_cartesian_p1 &
                      (dist_tab(1,i_point), i_r(1), &
                      dir_tab(1,1,i_point), index_lm,  l_wave_max, &
                      n_compute_c, i_basis, atom_index_inv, &
                      i_basis_fns_inv, kinetic_wave(1), &
                      kinetic_wave_deriv(1), &
                      ylm_tab(1,1), &
                      sum_gradient(1,1,1), &
                      kinetic_gradient_basis_wave(1:n_compute_c,1:3,i_point), & 
                      n_compute_atoms)
                 
!---------------shanghui end calclate kinetic_gradient_phi---------------------


                do i_spin=1,n_spin
                   do i_compute=1,n_compute_c
                   H_times_gradient_psi(i_compute,1:3,i_point,i_spin)= & 
                   kinetic_gradient_basis_wave(i_compute,1:3,i_point)+& 
                   local_potential_parts(i_spin)*  & 
                   gradient_basis_wave_npoints(i_compute,1:3,i_point)
                   enddo
                enddo 


              ! Now, evaluate vector of components H*phi(i,r)
              ! Local potential parts first; in the case of GGA,
              ! the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate
              ! H*phi(i,r) all in one

!!!              ! First, obtain radial kinetic energy terms from vectorized splines
!!!              call evaluate_radial_functions_p0  &
!!!                   ( spline_array_start, spline_array_end,  &
!!!                   n_compute_atoms, n_compute_fns,   &
!!!                   dist_tab(1,i_point), i_r,  &
!!!                   atom_index, i_basis_fns_inv,  &
!!!                   basis_kinetic_ordered, kinetic_wave(1),  &
!!!                   .false., n_compute_c, n_max_compute_fns_ham )


              do i_spin = 1, n_spin, 1
                 call evaluate_H_psi_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab,  &
                   radial_wave, H_times_psi(1, i_point, i_spin),  &
                   local_potential_parts(i_spin),  &
                   kinetic_wave, zora_operator(i_spin), &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )


              enddo

              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0


              if ((flag_rel.eq.REL_zora).or. flag_rel==REL_KOLNING_HARMON) then

                 ! Scalar relativistic treatment.
                 ! count number of "truly" relativistic points for ZORA treatment
                 ! of kinetic energy ...

                 do i_spin = 1, n_spin, 1

                    zora_operator(i_spin) =  &
                         light_speed_sq /  &
                         (2 * light_speed_sq -  &
                         zora_potential_parts(i_spin))**2

                    call  add_zora_gradient_part_p0(   &
                         sum_of_local_gradients(1,i_spin),  &
                         i_r_full,  &
                         dir_tab_full_norm,   &
                         dist_tab_full,  &
                         zora_operator(i_spin), &
                         n_centers_integrals, centers_basis_integrals )

                 end do

                 do i_spin = 1, n_spin, 1

                    ! Evaluate difference of scalar relativistic kinetic energy operator for the
                    ! true potential and the superposition of free atom potentials separately, and
                    ! only for all relativistic points in shell. Here, use partially
                    ! integrated version, leading to a vector:
                    ! zora_operator(r)*grad(phi(r,i))

                    zora_operator(i_spin) =  &
                         light_speed_sq *  &
                         (local_potential_parts(i_spin) -  &
                         zora_potential_parts(i_spin))/  &
                         ( 2 * light_speed_sq -  &
                         local_potential_parts(i_spin))/  &
                         ( 2 * light_speed_sq -  &
                         zora_potential_parts(i_spin))

                    call evaluate_zora_vector_p1  &
                         ( zora_operator(i_spin),  &
                         partition_tab(i_full_points),  &
                         gradient_basis_wave(1,1),  &
                         n_compute_c,  &
                         zora_vector1(1, 1, n_rel_points+1, i_spin),  &
                         zora_vector2(1, 1, n_rel_points+1, i_spin), &
                         n_max_compute_ham, t_zora(i_spin)  )

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
              !if (use_gga .or. (n_rel_points.gt.0)) then
              if (use_gga .or.(flag_rel.eq.REL_zora).or. flag_rel==REL_KOLNING_HARMON) then

                 do i_spin = 1, n_spin, 1
                    call add_gradient_part_to_H_p0  &
                         ( n_compute_c,   &
                         gradient_basis_wave(1,1),  &
                         sum_of_local_gradients(1,i_spin),  &
                         H_times_psi(1, i_point, i_spin) )
                 enddo
              end if
   

           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch


        do i_spin = 1, n_spin, 1

!          call  evaluate_hellman_feynman_delta_part_phonon_reduce_memory(   &
!                n_points, partition,      &
!                local_rho(i_spin,1:n_points),local_first_order_rho(1:n_points), &
!                coords_npoints(:,1:n_points), &
!                i_q_point, j_atom,j_coord, &
!                hellman_feynman_dynamical_matrix_delta_part )


            call evaluate_pulay_dynamical_matrix_phonon_reduce_memory( & 
                 partition,  &
                 n_points,n_compute_c, i_basis,  &
                 wave(:,1:n_points), gradient_basis_wave_npoints(:,:,1:n_points), &
                 index_hessian, hessian_basis_wave, &
                 H_times_psi(:,1:n_points,i_spin), &
                 local_first_order_rho(1:n_points), & 
                 local_first_order_potential(1:n_points), & 
                 H_times_gradient_psi(:,:,1:n_points,i_spin),  &
                 local_dVxc_drho(i_spin,1:n_points), &
                 coords_npoints(:,1:n_points), &
                 i_q_point, j_atom,j_coord, &
                 density_matrix_con,first_order_density_matrix_con,  &
                 energy_density_matrix_con, first_order_energy_density_matrix_con, &
                 pulay_dynamical_matrix,kinetic_gradient_basis_wave) 
   
        enddo


     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches

   
    !-------shanghui begin parallel------
!    if(.not. use_local_index) call sync_vector_complex(hellman_feynman_dynamical_matrix_delta_part,3*n_atoms*3*n_atoms)
    if(.not. use_local_index) call sync_vector_complex(pulay_dynamical_matrix(:,:,j_coord,j_atom),3*n_atoms)
    !-------shanghui end parallel------




  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_world,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  if(time_all>time_work*1.3 .and. .not.use_load_balancing) &
    info_str = trim(info_str) // ' => Consider using load balancing!'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)


  !     synchronise the hamiltonian
!  if(.not. use_local_index) call sync_integrate_hamiltonian( hamiltonian )
  

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)



  if(allocated( zora_vector2         )) deallocate( zora_vector2         )
  if(allocated( zora_vector1         )) deallocate( zora_vector1         )
  if(allocated( i_r_full             )) deallocate( i_r_full             )
  if(allocated( dir_tab_full_norm    )) deallocate( dir_tab_full_norm    )
  if(allocated( dist_tab_full        )) deallocate( dist_tab_full        )
  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( wave                 )) deallocate( wave                 )
  if(allocated( local_first_order_rho      )) deallocate( local_first_order_rho      )
  if(allocated( local_first_order_potential      )) deallocate( local_first_order_potential     )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( H_times_psi          )) deallocate( H_times_psi          )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if(allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if(allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
  if(allocated( gradient_basis_wave_npoints)) deallocate( gradient_basis_wave_npoints)
  if(allocated( coords_npoints)) deallocate( coords_npoints)
  if(allocated( local_rho_gradient   )) deallocate( local_rho_gradient   )
  if(allocated( local_rho            )) deallocate( local_rho            )
  if(allocated( xc_gradient_deriv    )) deallocate( xc_gradient_deriv    )
  if(allocated( local_xc_derivs      )) deallocate( local_xc_derivs      )
  if(allocated( local_dVxc_drho      )) deallocate( local_dVxc_drho      )
  if(allocated( en_density_xc        )) deallocate( en_density_xc        )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )

  if (allocated( density_matrix_con   )) deallocate( density_matrix_con   )
  if (allocated( first_order_density_matrix_con   )) deallocate( first_order_density_matrix_con   )
  if (allocated( energy_density_matrix_con   )) deallocate( energy_density_matrix_con   )
  if (allocated( first_order_energy_density_matrix_con   )) deallocate( first_order_energy_density_matrix_con   )

  if(use_batch_permutation > 0) then
    deallocate(rho)
    deallocate(hartree_potential)
    deallocate(rho_gradient) ! always allocated
    deallocate(ins_idx)
  endif

  if(get_batch_weights) deallocate(batch_times)



end subroutine integrate_dynamical_matrix_phonon_reduce_memory
!******
