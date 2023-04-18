!****h* FHI-aims/meta_gga_postp
!  NAME
!    meta_gga_postp - calculates meta_gga functioanls as post-processing
!  SYNOPSIS

module meta_gga_postp 

!  PURPOSE
!  This module contains all routines related to the evaluation of the meta-gga xc energy
!  as post processing
!
!  AUTHOR
!    Mariana Rossi
!  HISTORY
!    Development version, FHI-aims (2010).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Develoment version, FHI-aims (2010).
!  SOURCE
  implicit none

!save
! the kinetic density of the system \tau_{sigma} = \sum |grad . psi|^2
!real*8, dimension(:,:), allocatable :: kinetic_density

contains
!!******	
!!------------------------------------------------------------------------------
!!****s* molecular_dynamics/allocate_kdens
!!  NAME
!!    allocate_kdens
!!  SYNOPSIS
!subroutine allocate_kdens
!!  PURPOSE
!!    allocation of kinetic density
!!  USES
!  use dimensions
!  use runtime_choices
!!  ARGUMENTS
!!  INPUTS
!!    none
!!  OUTPUT
!!    none
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!
!  if (.not.allocated(kinetic_density)) then
!     allocate(kinetic_density(n_spin, n_full_points))
!  end if
!
!end subroutine allocate_kdens
!!******

!!------------------------------------------------------------------------------
!!****s* molecular_dynamics/clean_kdens
!  NAME
!!    clean_kdens - deallocation of kinetic_density
!!  SYNOPSIS
!subroutine clean_kdens
!!  PURPOSE
!!!    deallocation of local variables
!!  USES
!   use runtime_choices
!!  AUTHOR
!!    FHI-aims team.
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SEE ALSO
!!    FHI-aims CPC publication (in copyright notice above)
!!  SOURCE
!
!  implicit none
!  if (allocated(kinetic_density))  deallocate(kinetic_density)
!
!end subroutine clean_kdens
!******
!!****s* FHI-aims/evaluate_density_gradient_denmat
!!  NAME
!!   evaluate_density_gradient_denmat
!!  SYNOPSIS
!
!subroutine evaluate_kinetic_density_denmat( & 
!     n_points, gradient_basis_wave, n_compute, &
!     density_matrix,kinetic_density_own, n_basis_list)
!
!!  PURPOSE
!!  Evaluates kinetic density for one grid batch using density matrix formalism.
!!
!!  USES

!  use dimensions
!  use runtime_choices
!  implicit none
!
!!  ARGUMENTS
!
!  integer :: n_basis_list
!  integer :: n_points ! number of points in this batch
!  integer :: n_compute ! number of non zero basis functions in this batch
!
!  real*8, dimension(n_basis_list, 3, n_points) :: gradient_basis_wave ! input: the gradient of the basis
!  real*8, dimension(n_compute, n_compute) :: density_matrix ! input: the part of the density matrix relevant for this batch
!
!  real*8, dimension(n_points) :: kinetic_density_own ! output: the kinetic density for this batch
!
!  real*8, dimension(3, n_points) :: kinetic_density_vec ! temporary variable that stores the kinetic density before summing up in the coordinates
!  real*8:: vec_work(n_compute, n_points) ! temporary variable that stores the gradient_basis_wave for each coordinate
!  real*8:: work(n_compute, n_points) ! temporary variable that stores the result of DM * grad phi

!
!! INPUTS
!! o gradient_basis_wave
!! o density_matrix - see descriptions above
!  OUTPUT
!! o temp_kinetic_density - see descriptions above
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  SEE ALSO
!!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!!   the terms and conditions of the respective license agreement."
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!
!  !     counters
!  integer :: i_coord
!  integer :: i_point
!
!  !  begin work

!  kinetic_density_own(:)=0.d0
!
!  do i_coord = 1, 3, 1
!
!     vec_work =  gradient_basis_wave(1:n_compute, i_coord, 1:n_points)
!
!     ! do DM * grad phi for all basis n
!
!     call dsymm('L','U', n_compute, n_points, 1.0d0,  &
!!          density_matrix, n_compute, vec_work, n_compute, &
!          0.0d0, work, n_compute)
!
!     do i_point = 1, n_points,1
!  
!!        ! now do (DM * grad phi) * grad phi for all basis m, for each point
!
!        kinetic_density_vec(i_coord, i_point) = & 
!          dot_product(work(1:n_compute, i_point), gradient_basis_wave(1:n_compute, i_coord, i_point))
!
!       ! now sum over the coordinates, for each point
!
!        kinetic_density_own(i_point) = kinetic_density_own(i_point)+kinetic_density_vec(i_coord, i_point)

!     end do
!
!  end do
!
!
!end subroutine evaluate_kinetic_density_denmat
!!---------------------------------------------------------------------
!!******	 
!!-----------------------------------------------------
!****s* FHI-aims/calculate_xc_energy_meta_gga
!  NAME
!   calculate_xc_energy_meta_gga
!  SYNOPSIS
! AJL, Feb2018: Disabling as this is a duplicate routine (see calculate_kinetic_density.f90)
!!subroutine calculate_kinetic_density &
!     ( rho, rho_gradient, hartree_partition_tab, &
!     partition_tab, basis_l_max, &
!     KS_eigenvector, KS_eigenvector_complex, occ_numbers, &
!     kinetic_density)
!
!!  PURPOSE
!!  Calculates the basis derivatives and density matrix required to
!!  build the kinetic density. This is a bit of code duplication, but
!!  a lot of modification would have to be done in the code overall to
!!  get these quantities otherwise. 
!!  This actually modifies the variable kinetic_density
!!  USES
!
!  use dimensions
!  use runtime_choices
!  use grids
!  use geometry
!  use mixing
!  use species_data
!  use mpi_tasks
!  use mpi_utilities
!  use synchronize_mpi
!  use localorb_io
!  use basis
!  use cartesian_ylm
!  use constants
!  use analyze_arrays
!  use scalapack_wrapper
!  use density_matrix_evaluation
!  use sym_base, only: evaluate_densmat_sym
!  use pbc_lists
!  implicit none
!
!!  ARGUMENTS
!
!  real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
!  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
!  real*8,     dimension(n_states, n_spin, n_k_points)             :: occ_numbers
!
!  real*8,  dimension(n_full_points)             ::  partition_tab
!  real*8,  dimension(n_full_points)             ::  hartree_partition_tab
!  real*8,  dimension(n_spin, n_full_points)     :: rho
!  real*8,  dimension(3, n_spin, n_full_points)  :: rho_gradient
!  integer, dimension(n_species)                 :: basis_l_max
!  real*8,  dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse
!  real*8,  dimension(n_spin, n_full_points)     :: kinetic_density
!
!! INPUTS
!! o KS_eigenvector -- eigenvectors if real eigenvectors are in use
!! o KS_eigenvector_complex -- eigenvectors is complex eigenvectors are in use
!! o occ_numbers -- occupation numbers of states
!! o partition_tab -- partition function values
!! o rho -- electron density, rho is only input 
!! o rho_gradient -- gradient of the electron density 
!! o basis_l_max -- maximum basis l value of the basis functions
!! OUTPUTS
!! o kinetic_density -- kinetic density calcualted from orbitals
!!
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!  SEE ALSO
!!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!!    Computer Physics Communications (2008), submitted.
!!  COPYRIGHT
!!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!!   the terms and conditions of the respective license agreement."
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!
!
!  real*8, dimension(:),     allocatable:: dist_tab
!  real*8, dimension(:,:),   allocatable:: dist_tab_sq
!  real*8, dimension(:),     allocatable:: i_r
!  real*8, dimension(:,:,:), allocatable:: dir_tab
!  real*8, dimension(:,:),   allocatable:: trigonom_tab
!
!  real*8,dimension(:),  allocatable:: radial_wave
!  real*8,dimension(:),  allocatable:: radial_wave_deriv
!  real*8,dimension(:),  allocatable:: radial_wave_2nd_deriv
!  real*8,dimension(:,:),allocatable:: wave
!  real*8, dimension(n_spin) :: x_energy 
!  real*8 :: c_energy
!  integer,dimension(:),allocatable :: i_basis
!  integer :: n_compute_c
!  real*8 coord_current(3)
!
!  integer :: n_compute_fns
!  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
!  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
!  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)
!
!  integer :: n_compute_atoms
!  integer :: atom_index(n_centers_integrals)
!  integer :: atom_index_inv(n_centers)
!
!  integer :: spline_array_start(n_centers_integrals)
!  integer :: spline_array_end(n_centers_integrals)
!
!! VB - renewed index infrastructure starts here
!
!      real*8 :: one_over_dist_tab(n_max_compute_atoms)
!
!      ! indices for basis functions that are nonzero at current point
!
!      integer :: rad_index(n_max_compute_atoms)
!      integer :: wave_index(n_max_compute_fns_dens)
!      integer :: l_index(n_max_compute_fns_dens)
!      integer :: l_count(n_max_compute_fns_dens)
!      integer :: fn_atom(n_max_compute_fns_dens)
!
!      ! indices for known zero basis functions at current point
!      integer :: n_zero_compute
!      integer :: zero_index_point(n_max_compute_dens)
!
!      ! active atoms in current batch
!      integer :: n_batch_centers
!      integer :: batch_center(n_centers_integrals)
!
!  !     other local variables
!
!  integer :: l_ylm_max
!  integer :: n_points
!
!  ! (n_basis, n_states, n_spin) on purpose
!  ! beware: inside the subroutines (evaluate_KS_density_v1 for instance)
!  ! dimension is (n_compute, max_occ_number)
!  ! that's a trick to get a continuous data flow and not
!  ! bad programming
!
!
!  integer, dimension(:,:), allocatable :: index_lm
!  real*8,  dimension(:,:), allocatable :: ylm_tab
!
!  !     things for gradient functionals
!  real*8, dimension(:,:),  allocatable :: dylm_dtheta_tab
!  real*8, dimension(:,:),  allocatable :: scaled_dylm_dphi_tab
!  real*8, dimension(:,:,:),allocatable :: gradient_basis_wave
!  real*8, dimension(:,:),  allocatable :: density_matrix
!  real*8, dimension(:,:),  allocatable :: density_matrix_con
!  real*8, dimension(:,:),  allocatable :: work
!  real*8, dimension(:,:),  allocatable :: grad_work
!  complex*16, dimension(:,:),allocatable :: work_complex
!  real*8, dimension(:), allocatable :: temp_kinetic_density
!
!
!  !     counters
!  integer :: i_l
!  integer :: i_m
!  integer :: i_coord
!  integer :: i_state
!  integer :: i_point
!
!  integer :: i_full_points
!  integer :: i_full_points_2
!  integer :: i_full_points_3
!
!  integer :: i_spin,  i_spin_2
!  integer :: i_index, i_k_point
!  integer :: i_my_batch
!  integer :: info
!  logical :: is_restarted
!
!! temporary debug stuff
!!    real*8 :: i_full_points_4 
!!    real*8, dimension(n_spin) :: debug_kinetic_energy
!
!
!  call localorb_info( &
!       "Calculating meta-GGA energy as post processing, using DM", use_unit,'(2X,A)', OL_norm )
!
!  ! initialize variables
!  !  en_post_xc = 0.d0
!
!! begin allocations
!
!  if(packed_matrix_format == PM_none)then
!     if(.not. allocated(density_matrix))then
!        allocate(density_matrix(n_centers_basis_T, n_centers_basis_T),stat=info)
!        call check_allocation(info, 'density_matrix                ') 
!     end if
!  end if
!
!! begin outermost loop over spin
!  do i_spin = 1, n_spin
!
!! build the density matrix for each spin
!     is_restarted = .false.
!     if (use_scalapack .and. restart_read) then
!        call restart_scalapack_read(density_matrix_sparse, i_spin, is_restarted)
!     end if
!
!     if (.not. is_restarted) then
!	if(use_symmetry_reduced_spg)then
!          call evaluate_densmat_sym(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
!          &                     density_matrix, density_matrix_sparse, i_spin, .false.)
!        else
!          call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
!          &                     density_matrix, density_matrix_sparse, i_spin, .false.)
!        endif
!     end if
!
!     if (use_scalapack .and. restart_write) then
!        call restart_scalapack_write(density_matrix_sparse, i_spin)
!     end if
!
!
!! ---------------- end build density matrix -------------------------
!! continue allocating stuff
!     if(.not. allocated(dist_tab))then
!        allocate(dist_tab(n_centers_integrals),stat=info)
!        call check_allocation(info, 'dist_tab                      ') 
!     end if
!     if(.not. allocated(dist_tab_sq))  then
!        allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
!        call check_allocation(info, 'dist_tab_sq                   ')
!     end if
!     if(.not. allocated(i_r))then
!        allocate(i_r(n_max_compute_atoms),stat=info)
!        call check_allocation(info, 'i_r                           ')
!     end if
!     if(.not. allocated(dir_tab)) then
!        allocate(dir_tab(3, n_centers_integrals, n_max_batch_size),stat=info)
!        call check_allocation(info, 'dir_tab                       ')
!     end if
!     if(.not. allocated(trigonom_tab)) then
!        allocate(trigonom_tab(4, n_max_compute_atoms),stat=info)
!        call check_allocation(info, 'trigonom_tab                  ')
!     end if
!     if(.not. allocated(density_matrix_con))then
!        allocate(density_matrix_con(n_max_compute_dens, n_max_compute_dens),stat=info)
!        call check_allocation(info, 'density_matrix_con            ')
!     end if
!     if(.not. allocated(work))then
!        allocate(work(n_max_compute_dens, n_max_batch_size),stat=info)
!        call check_allocation(info, 'work                          ')
!     end if
!     if(.not. allocated(grad_work))then
!        allocate(grad_work(n_max_compute_dens, n_max_batch_size),stat=info)
!        call check_allocation(info, 'grad_work                     ')
!     end if
!     if(.not. allocated(temp_kinetic_density))then
!        allocate(temp_kinetic_density(n_max_batch_size),stat=info)
!        call check_allocation(info, 'temp_kinetic_density          ')
!     end if
!     if(.not. allocated(i_basis))then
!        allocate(i_basis(n_centers_basis_T),stat=info)
!        call check_allocation(info, 'i_basis                       ')
!     end if
!     if(.not.allocated(radial_wave))then
!        allocate(radial_wave(n_max_compute_fns_dens),stat=info)
!        call check_allocation(info, 'radial_wave                   ')
!     end if
!     if(.not. allocated(radial_wave_deriv))then
!        allocate(radial_wave_deriv(n_max_compute_fns_dens),stat=info)
!        call check_allocation(info, 'radial_wave_deriv             ')
!     end if
!     if(.not. allocated(radial_wave_2nd_deriv))then
!        allocate(radial_wave_2nd_deriv(n_max_compute_fns_dens),stat=info)
!        call check_allocation(info, 'radial_wave_2nd_deriv         ')
!     end if
!     if(.not. allocated(wave))then
!        allocate(wave(n_max_compute_dens, n_max_batch_size),stat=info)
!        call check_allocation(info, 'wave                          ')
!     end if
!
!     l_ylm_max = l_wave_max
!
!     if(.not. allocated( ylm_tab))then
!        allocate( ylm_tab( (l_ylm_max+1)**2,n_max_compute_atoms),stat=info )
!        call check_allocation(info, 'ylm_tab                       ')
!     end if
!     if(.not. allocated( index_lm))then
!        allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info)
!        call check_allocation(info, 'index_lm                      ')
!     end if
!
!
!!   allocate local arrays needed for gradients
!        if(.not. allocated (gradient_basis_wave))then
!           allocate (gradient_basis_wave(n_max_compute_dens, 3, n_max_batch_size),stat=info)
!           call check_allocation(info, 'gradient_basis_wave           ')
!        end if
!        if(.not. allocated( dylm_dtheta_tab))then
!           allocate( dylm_dtheta_tab((l_ylm_max+1)**2, n_max_compute_atoms),stat=info )
!           call check_allocation(info, ' dylm_dtheta_tab              ')
!         end if
!        if(.not. allocated( scaled_dylm_dphi_tab))then
!           allocate( scaled_dylm_dphi_tab((l_ylm_max+1)**2,n_max_compute_atoms),stat=info )
!           call check_allocation(info, 'scaled_dylm_dphi_tab          ')
!        end if
!! finally all allocations are done
!
!
!     ! initialize index_lm
!     i_index = 0
!     do i_l = 0, l_ylm_max, 1
!        do i_m = -i_l, i_l
!           i_index = i_index + 1
!           index_lm(i_m, i_l) = i_index
!        enddo
!     enddo
!
!! -------- go over the grids and construct the kinetic density using density matrix -------------
!! first initialize variables
!
!     kinetic_density(i_spin,:) = 0.d0
!     i_basis_fns_inv = 0
!
!     i_full_points = 0
!     i_full_points_2 = 0
!     i_full_points_3 = 0
!
!! go over the batches
!     do i_my_batch = 1, n_my_batches, 1
!
!           temp_kinetic_density(:)=0.d0
!
!           n_compute_c = 0
!           i_basis = 0
!
!           i_point = 0
!
!! go over the points inside each batch
!           do i_index = 1, batches(i_my_batch)%size, 1
!
!!debug
!! write(use_unit,*)"batch sizes", batches(i_my_batch)%size
!!end debug
!              i_full_points = i_full_points + 1
!
!!               if (partition_tab(i_full_points).gt.0) then
!              if (max(partition_tab(i_full_points),&
!                   hartree_partition_tab(i_full_points)).gt.0.d0) then
!
!                 i_point = i_point+1
!
!                 !     get current integration point coordinate
!                 coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
!
!                 if(n_periodic>0)then
!                    call map_to_center_cell( coord_current)
!                 end if
!
!                 ! compute atom-centered coordinates of current integration point,
!                 ! as viewed from all atoms
!                 call tab_atom_centered_coords_p0 &
!                      ( coord_current,  &
!                      dist_tab_sq(1,i_point),  &
!                      dir_tab(1,1,i_point), &
!                      n_centers_integrals, centers_basis_integrals )
!
!
!
!                 !     determine which basis functions are relevant at current integration point,
!                 !     and tabulate their indices
!
!                 ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
!                 if (.not.prune_basis_once) then
!                   call prune_basis_p2 &
!                   ( dist_tab_sq(1,i_point), &
!                     n_compute_c, i_basis,  &
!                     n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals  )
!                 end if !prune basis
!              end if !partition_tab
!           enddo ! end loop over the angular shell
!
!
!           if (prune_basis_once) then
!              n_compute_c = batches(i_my_batch)%batch_n_compute
!              i_basis(1:n_compute_c) = batches(i_my_batch)%batch_i_basis
!           end if
!
!           if(packed_matrix_format /= PM_none )then
!              call  prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
!                   n_compute_c, i_basis)
!           else
!              call  prune_density_matrix(density_matrix, density_matrix_con, &
!                   n_compute_c, i_basis)
!           end if
!
!
!           ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
!           call collect_batch_centers_p2 &
!           ( n_compute_c, i_basis, n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals, &
!             n_batch_centers, batch_center &
!           )
!
!           n_points = i_point
!
!           if (n_compute_c.gt.0) then
!
!              ! Determine all radial functions, ylm functions and their derivatives that
!              ! are best evaluated strictly locally at each individual grid point.
!              i_point = 0
!              do i_index = 1, batches(i_my_batch)%size, 1
!
!                 i_full_points_2 = i_full_points_2 + 1
!
!!                  if (partition_tab(i_full_points_2).gt.0) then
!                 if (max(partition_tab(i_full_points_2),&
!                   hartree_partition_tab(i_full_points_2)).gt.0.d0) then
!
!                    i_point = i_point+1
!                    n_compute_atoms = 0
!                    n_compute_fns = 0
!
!
!                    ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
!                    ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
!                    ! without any copying and without doing any unnecessary operations. 
!                    ! The price is that the interface is no longer explicit in terms of physical 
!                    ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
!
!
!                    call prune_radial_basis_p2 &
!                    ( n_max_compute_atoms, n_max_compute_fns_dens, &
!                       dist_tab_sq(1,i_point), dist_tab, dir_tab(1,1,i_point), &
!                       n_compute_atoms, atom_index, atom_index_inv, &
!                       n_compute_fns, i_basis_fns, i_basis_fns_inv, &
!                       i_atom_fns, spline_array_start, spline_array_end, &
!                       n_centers_basis_integrals, centers_basis_integrals, n_compute_c, i_basis, &
!                       n_batch_centers, batch_center, &
!                       one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
!                       fn_atom, n_zero_compute, zero_index_point &
!                    )
!
!                    ! Tabulate distances, unit vectors, and inverse logarithmic grid units
!                    ! for all atoms which are actually relevant
!                    call tab_local_geometry_p2 &
!                      ( n_compute_atoms, atom_index, &
!                        dist_tab,  &
!                        i_r &
!                      )
!
!                    ! Determine all needed radial functions from efficient splines
!
!                    ! Now evaluate radial functions u(r) from the previously stored compressed 
!                    ! spline arrays  
!                    call evaluate_radial_functions_p0 &
!                         (   spline_array_start, spline_array_end, &
!                         n_compute_atoms, n_compute_fns,  &
!                         dist_tab, i_r(1), &
!                         atom_index, i_basis_fns_inv, &
!                         basis_wave_ordered, radial_wave(1), &
!                         .false., n_compute_c, n_max_compute_fns_dens   &
!                         )
!
!                    ! for forces or density gradient, radial derivatives are required. Here always required. 
!                       call evaluate_radial_functions_p0 &
!                            (spline_array_start, spline_array_end, &
!                            n_compute_atoms, n_compute_fns,  &
!                            dist_tab, i_r(1), &
!                            atom_index, i_basis_fns_inv, &
!                            basis_deriv_ordered,  &
!                            radial_wave_deriv(1),  &
!                            .true., n_compute_c, n_max_compute_fns_dens  &
!                            )
!                    ! We compute the needed
!                    ! ylm pieces separately from the cartesian evaluation that is done
!                    ! for the hessians ...
!                    
!                    ! compute trigonometric functions of spherical coordinate angles
!                    ! of current integration point, viewed from all atoms
!
!                       call tab_trigonom_p0 &
!                            ( n_compute_atoms, dir_tab(1,1,i_point),  &
!                            trigonom_tab(1,1) &
!                            )
!
!                          ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
!                       call tab_gradient_ylm_p0 &
!                               (trigonom_tab(1,1), basis_l_max,  &
!                               l_ylm_max, n_compute_atoms, atom_index, &
!                               ylm_tab(1,1),  &
!                               dylm_dtheta_tab(1,1),  &
!                               scaled_dylm_dphi_tab(1,1) &
!                               )
!
!                          ! evaluate the wave function gradient directly here                  
!                       call evaluate_wave_gradient_p2  &
!                          ( n_compute_c, n_compute_atoms, n_compute_fns, &
!                            one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
!                            l_ylm_max, ylm_tab,  &
!                            dylm_dtheta_tab,  &
!                            scaled_dylm_dphi_tab,  &
!                            radial_wave,  &
!                            radial_wave_deriv,  &
!                            gradient_basis_wave (1:n_compute_c,1:3,i_point),  &
!                            rad_index, wave_index, l_index, l_count, fn_atom, &
!                            n_zero_compute, zero_index_point &
!                          )
!
!
!                    ! reset i_basis_fns_inv
!                    i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
!
!
!                 end if ! end if (partition_tab.gt.0)
!              enddo ! end loop over a batch
!
!! - all quantities that are evaluated pointwise are now known ...
!
!              if (n_points.gt.0) then
!! Now perform all operations which are done across the entire integration
!! shell at once.
!
!! In this case, really finally evaluate the kinetic density
!                    call evaluate_kinetic_density_denmat  &
!                         (n_points, gradient_basis_wave, n_compute_c,  &
!                         density_matrix_con,  temp_kinetic_density,  &
!                         n_max_compute_dens)
!
!              end if ! (n_points>0)
!
!! assign quantities to their right place
!              i_point = 0
!              do i_index = 1, batches(i_my_batch)%size, 1
!
!                 i_full_points_3 = i_full_points_3 + 1
!
!
!!                  if (partition_tab(i_full_points_3).gt.0) then
!                 if (max(partition_tab(i_full_points_3),&
!                      hartree_partition_tab(i_full_points_3)).gt.0) then
!
!                    i_point = i_point + 1                    
!
!                    kinetic_density(i_spin, i_full_points_3) =   &
!                            temp_kinetic_density(i_point)
!
!                 end if
!               end do
!
!
!       else
!
!        i_full_points_2 = i_full_points_2 + batches(i_my_batch)%size
!        i_full_points_3 = i_full_points_3 + batches(i_my_batch)%size
!
!      end if ! end if (n_compute.gt.0) then
!
!    end do ! end loop over grid batches
!     !--------------------------------- end go over the grids ------------------------------
!
!  end do ! end i_spin
!
!! DEBUG CHECK: calculate kinetic energy by integrating kinetic energy density
!! I leave this piece of code here because it is very useful if one wants to get the kinetic energy of the system 
!! using the kinetic density -- MR. 
!!
!!do i_spin = 1, n_spin, 1
!!  debug_kinetic_energy(i_spin)=0.d0
!!  i_full_points_4 = 0
!!
!!  do  i_my_batch = 1, n_my_batches, 1
!!
!!        do i_index = 1, batches(i_my_batch)%size, 1
!!
!!           i_full_points_4 = i_full_points_4 + 1
!!           if(partition_tab(i_full_points_4) .gt. 0.d0) then
!!
!!                 
!!                 debug_kinetic_energy(i_spin) = debug_kinetic_energy(i_spin) + 0.5*kinetic_density(i_spin, i_full_points_4)*partition_tab(i_full_points_4)
!!
!!           endif
!!
!!        enddo ! end of mpi distribution
!!   
!!  enddo! end integration loop over batches
!!
!!enddo ! i_spin
!!
!!
!!  if (use_mpi) then
!!        call sync_vector(debug_kinetic_energy,n_spin)
!!  end if
!!  if (myid==0) then
!!    write(use_unit,*) "   Kinetic energy from kinetic density: ", debug_kinetic_energy
!!  endif
!! END OF DEBUG CHECK  
! !---------- finally, deallocate stuff -----------------
!  if (allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
!  if (allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
!  if (allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
!  if (allocated( index_lm             )) deallocate( index_lm             )
!  if (allocated( ylm_tab              )) deallocate( ylm_tab              )
!  if (allocated( wave                 )) deallocate( wave                 )
!  if (allocated( radial_wave_2nd_deriv)) deallocate( radial_wave_2nd_deriv)
!  if (allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
!  if (allocated( radial_wave          )) deallocate( radial_wave          )
!  if (allocated( i_basis              )) deallocate( i_basis              )
!  if (allocated( work                 )) deallocate( work                 )
!  if (allocated( grad_work            )) deallocate( grad_work            )
!  if (allocated( density_matrix_con   )) deallocate( density_matrix_con   )
!  if (allocated( trigonom_tab         )) deallocate( trigonom_tab         )
!  if (allocated( dir_tab              )) deallocate( dir_tab              )
!  if (allocated( i_r                  )) deallocate( i_r                  )
!  if (allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
!  if (allocated( dist_tab             )) deallocate( dist_tab             )
!  if (allocated( work_complex         )) deallocate( work_complex         )
!  if (allocated( density_matrix       )) deallocate( density_matrix       )
!  if (allocated( temp_kinetic_density )) deallocate( temp_kinetic_density )
!
!  return
!
!end subroutine calculate_kinetic_density
!******
! End AJL Disabled.
!---------------------------------------------------------------------------
!****s* FHI-aims/integrate_post_xc_energy
!  NAME
!   integrate_post_xc_energy
!  SYNOPSIS

  subroutine integrate_post_xc_energy &
       ( partition_tab, &
         rho,rho_gradient, &
         en_post_xc, &
         x_energy_post, &
         c_energy_post, &
         kinetic_density )

!  PURPOSE
! Subroutine designed to integrate xc energy density for meta-gga post processing
! kinetic_density is used inside but is not modified
!
!  USES

    use dimensions
    use runtime_choices
    use grids
    use xc
    use mpi_tasks
    use mpi_utilities
    use synchronize_mpi
    implicit none

!  ARGUMENTS
!     input
      real*8, dimension(n_full_points)              :: partition_tab
      real*8, dimension(n_spin,n_full_points)       :: rho
      real*8, dimension(3,n_spin,n_full_points)     :: rho_gradient
      real*8, optional, dimension(n_spin, n_full_points)      :: kinetic_density
!     output
      real*8 :: en_post_xc
      real*8 :: x_energy_post
      real*8 :: c_energy_post

!  INPUTS
!  o partition_tab -- the partition function for real space integration     
!  o rho -- the charge density
!  o rho_gradient -- the gradient of the charge density
!  o kinetic_density -- the kinetic density

!  OUTPUTS
!  o xc_energy -- the exchange-correlation energy    
!  o x_energy -- the exchange energy    
!  o c_energy -- the correlation energy
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


  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  real*8 :: en_density_xc
  real*8 :: en_density_x
  real*8 :: en_density_c
!  real*8 :: local_xc_density_deriv(n_spin)
!  real*8 :: local_xc_gradient_deriv(3,n_spin)

  real*8 :: xc_energy_mpi
  real*8, dimension(n_spin) :: pot_xc

  integer :: mpierr

  logical :: xc_undefined

  !  counters

  integer i_my_batch
  integer i_index
  integer i_spin


  integer i_full_points
  real*8 :: coord_current(3)

  !  begin work
  !     initialize

  en_post_xc = 0.d0
  x_energy_post = 0.d0
  c_energy_post = 0.d0

  i_full_points = 0
  !     perform partitioned integration, atom by atom, and point by point
  !     This will be the outermost loop, to save evaluations of the potential.
  !     and the Y_lm functions

  do  i_my_batch = 1, n_my_batches, 1



        do i_index = 1, batches(i_my_batch)%size, 1



           i_full_points = i_full_points + 1


           if(partition_tab(i_full_points) .gt. 0.d0) then


              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)


              call evaluate_post_xc &
                   ( rho(1:n_spin,i_full_points), &
                   rho_gradient(1:3,1:n_spin, i_full_points), &
                   en_density_xc,  en_density_x, en_density_c, &
!                   local_xc_density_deriv(1:n_spin), &
!                   local_xc_gradient_deriv(1:3,1:n_spin), &
                   coord_current(:), &
                   kinetic_density(1:n_spin, i_full_points) &
                   )

              do i_spin = 1, n_spin, 1
                 en_post_xc = en_post_xc + &
                      rho(i_spin,i_full_points)* &
                      en_density_xc * &
                      partition_tab(i_full_points)

                 x_energy_post = x_energy_post + &
                      rho(i_spin,i_full_points)* &
                      en_density_x * &
                      partition_tab(i_full_points)

                 c_energy_post = c_energy_post + &
                      rho(i_spin,i_full_points)* &
                      en_density_c * &
                      partition_tab(i_full_points)
!

              enddo !i_spin

           endif! partition_tab

        enddo !i_index
        !       end of mpi distribution

     !     end integration loop over batches
  enddo !i_my_batch

  if(use_mpi) then
     call sync_real_number(en_post_xc)
     call sync_real_number(x_energy_post)
     call sync_real_number(c_energy_post)
  endif

  return
end subroutine integrate_post_xc_energy

!!----------------------------------------------------------------------
!!************
!!****s* FHI-aims/evaluate_post_xc
!!  NAME
!!  evaluate_post_xc
!!  SYNOPSIS
!      subroutine evaluate_post_xc_old &
!      ( rho, rho_gradient, en_density_xc, en_density_x, en_density_c, &
!!        local_xc_derivs, xc_gradient_deriv,
!        coord_current, kinetic_density_pointwise  &
!      )
!
!!  PURPOSE
!!  Subroutine evaluate_post_xc evaluates the exchange correlation potential
!!  and energy contribution for a given density at one integration point.
!!!
!!
!!  USES
!!
!      use runtime_choices
!      use dimensions
!      use xc
!      use constants
!      implicit none
!
!!  ARGUMENTS
!
!!      real*8 rho(n_spin)
!!      real*8 rho_gradient(3,n_spin)
!!      real*8 en_density_xc, en_density_x, en_density_c
!!      real*8 en_density_x2(n_spin)
!!      real*8 local_xc_derivs(n_spin)
!!      real*8 xc_gradient_deriv(3,n_spin)
!!      real*8 :: coord_current(3)
!!!      real*8, optional :: kinetic_density_pointwise(n_spin)
!!      !real*8, optional :: coord_current(3)
!!      real*8 xc_tau_deriv(n_spin)
!
!      
!!  INPUTS
!!   o xc_gradient_deriv : Partial derivative of the exchange-correlation
!!               energy functional by the modulus square of the density gradient
!!  AUTHOR
!!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!!    Updated by Andrew Logsdail, Jan 2015, for SCF Meta-GGA calculations
!!  SEE ALSO
!!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!!    Computer Physics Communications (2008), submitted.
!!  COPYRIGHT
!!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!!   the terms and conditions of the respective license agreement."
!!  HISTORY
!!    Release version, FHI-aims (2008).
!!  SOURCE
!
!!  local variables
!
!      logical :: xc_undefined
!
!!     counters
!
!      integer :: i_spin
!
!!  begin work
!
!        xc_undefined = .true.
!
!        do i_spin = 1, n_spin, 1
!          ! This is true if both densities are below zero
!          xc_undefined = xc_undefined .and. (rho(i_spin).le.0.d0)
!        enddo
!
!        do i_spin = 1, n_spin, 1
!          ! This is true if one of the densities is less than zero
!          xc_undefined = xc_undefined .or. (rho(i_spin).lt.0.d0)
!        enddo
!
!        en_density_xc = 0.d0
!        en_density_x =0.d0
!        en_density_c=0.d0
!        local_xc_derivs = 0.d0
!        xc_gradient_deriv = 0.d0
!        xc_tau_deriv = 0.d0
!
!        if (xc_undefined) then
!
!           continue
!
!!          en_density_xc = 0.d0
!!          en_density_x =0.d0
!!          en_density_c=0.d0
!!          local_xc_derivs = 0.d0
!!          xc_gradient_deriv = 0.d0
!!          xc_tau_deriv = 0.d0
!
!        else
!
!!!         calculate xc-potential, thus initializing new potential
!          
!          select case(flag_post_xc)
!           case(0)
!!        Hartree-Fock calculation, set the XC contribution to zero here.
!!            en_density_xc = 0.d0
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!          case(1)
!!           M06L functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_partials_m06l &
!            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c, &
!              local_xc_derivs, xc_gradient_deriv, xc_tau_deriv)
!
!          case(2)
!!           M06 functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_m06 &
!            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c)
!
!         case(3)
!!          pbe_vdw XC - post processing
!            call xc_partials_pbe_vdw &
!                 ( rho, rho_gradient, en_density_xc,  en_density_x2, &
!                 en_density_c, local_xc_derivs, xc_gradient_deriv, &
!                 coord_current)
!            en_density_x = en_density_x2(1)+ en_density_x2(2) !strictly, en_density_x should be an array...
!            
!            
!         case(4)
!!          revpbe_vdw XC - post processing
!            call xc_partials_revpbe_vdw &
!                 ( rho, rho_gradient, en_density_xc, &
!                 en_density_x2, en_density_c, local_xc_derivs, &
!                 xc_gradient_deriv, coord_current )
!            en_density_x = en_density_x2(1)+en_density_x2(2)
!
!          case(5)
!!           M06-2X functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_m062x &
!            ( rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c)
!            
! !@@edu>> tpss
!          case(6)
!!           TPSS functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_tpss &
!            (1, rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c)
!
!          case(7)
!!           revTPSS functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_tpss &
!            (2, rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c)
!
!          case(8)
!!           TPSSloc functional - post processing
!!            local_xc_derivs = 0.d0
!!            xc_gradient_deriv = 0.d0
!
!            call xc_tpss &
!            (3, rho, rho_gradient, kinetic_density_pointwise, en_density_xc, en_density_x, en_density_c)
!                 
!!@@edu<<  tpss
!           
!         case default
!            write(use_unit,*) "Chosen type of XC post-processing is not yet implemented."
!            stop
!          end select
!
!        end if
!
!
!      end subroutine evaluate_post_xc_old
!!---------------------------------------------------------------------

end module meta_gga_postp
