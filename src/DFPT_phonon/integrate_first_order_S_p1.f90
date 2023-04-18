!****s* FHI-aims/integrate_first_order_S_p1
!  NAME
!   integrate_first_order_S_p1
!  SYNOPSIS

subroutine integrate_first_order_S_p1 &
     (partition_tab_std, basis_l_max, first_order_S_sparse) 

!  PURPOSE
!  
!  Integrates the matrix elements for the first_order overlap(S) matrix,
!  using a fixed basis set. 

!  called by a SCF subroutine

!  shanghui 2012.05.08 : created 
!  shanghui 2013.12.20 : change to phonon_gamma version (p0)
!  shanghui 2014.11.11 : change to phonon only with sparse matrix.
!  
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
  use pbc_lists
  use species_data, only: species_name
  use load_balancing
  use physics, only : overlap_matrix
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points)            :: partition_tab_std
  integer ::  basis_l_max (n_species)
  !shanghui------------------------------------------------------------------------
  real*8, dimension(3, n_centers_in_sc_DFPT,n_hamiltonian_matrix_size),intent(inout) :: & 
          first_order_S_sparse
  !shanghui------------------------------------------------------------------------

!  INPUTS
!  o partition_tab_std -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!
!  OUTPUT
!  o first_order_S_sparse -- first_order overlap matrix
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

  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave


  real*8, dimension(:,:,:,:), allocatable :: matrix_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave
  real*8, dimension(:,:,:), allocatable :: gradient_basis_wave_npoints



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

  integer i_basis_1
  integer i_basis_2
  integer i_atom, i_atom_2, i_center
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_spin
  character*2000 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

!  integer ld_hamiltonian  ! leading dimension of hamiltonian in calling routine

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i_off, i, j, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

!----------------------------shanghui debug DFPT_phonon: the supercell fd-benchmark------------------------
  integer i_cell_1, i_want, n_basis_uc
  real*8  hamiltonian_i_want(1000)
  logical, parameter :: supercell_fd_benchmark = .false.
!----------------------------shanghui end debug DFPT_phonon: the supercell fd-benchmark--------------------

  ! begin work

  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A)') "Integrating first_order_S_sparse matrix: batch-based integration with load balancing"
  else
    write(info_str,'(2X,A)') "Integrating first_order_S_sparse matrix: batch-based integration."
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




     l_ylm_max = l_wave_max
     allocate (gradient_basis_wave(n_max_compute_ham,3),STAT=info)
     call check_allocation(info, 'gradient_basis_wave           ')
     allocate (gradient_basis_wave_npoints(n_max_compute_ham,3, & 
     n_max_batch_size),STAT=info)
     call check_allocation(info, 'gradient_basis_wave_npoints     ')

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')




  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info )
  call check_allocation(info, 'index_lm                      ')


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

!  write(info_str,'(a,2(f12.3,a))') ' in evaluate_first_order_S_p1: ', &
!                  8*dble(3*n_max_compute_ham*n_max_compute_ham*n_max_compute_ham)/1048576,  'MN'
!  call localorb_info ( info_str,use_unit,'(a)', ol_norm )


  allocate ( matrix_shell(3, n_max_compute_ham,n_max_compute_ham,n_max_compute_ham) )


  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points (for load balancing)
  ! or to standard batches / arrays (no load balancing)

  n_bp = use_batch_permutation
  if(use_batch_permutation > 0) then

    n_my_batches_work = batch_perm(n_bp)%n_my_batches
    batches_work => batch_perm(n_bp)%batches
    partition_tab => batch_perm(n_bp)%partition_tab



    allocate(ins_idx(batch_perm(n_bp)%n_basis_local))

    !ld_hamiltonian = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    !ld_hamiltonian = n_hamiltonian_matrix_size

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize
  first_order_S_sparse= 0.0d0

  i_basis_fns_inv = 0

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


              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)


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

                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )


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
                   n_zero_compute, zero_index_point  & 
                 )



              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0


           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch


!---------------------shanghui begin first_order_S--------------------
          call evaluate_first_order_S_p1    & 
              ( n_points, partition(1:n_points), &
               n_compute_c,i_basis,   & 
               wave, gradient_basis_wave_npoints, & 
               first_order_S_sparse)    
!          call evaluate_first_order_S_p1    & 
!              (matrix_shell(1:3,1:n_compute_c,1:n_compute_c,1:n_compute_c), & 
!               n_points, partition(1:n_points), &
!               n_compute_c,i_basis,   & 
!               wave, gradient_basis_wave_npoints)    

!          call update_full_matrix_p1 &
!              (n_compute_c, n_compute_c,  &
!               i_basis, &
!               matrix_shell(1:3,1:n_compute_c,1:n_compute_c,1:n_compute_c),  &
!               first_order_S_sparse)
!---------------------shanghui end first_order_S--------------------

     else

       i_full_points = i_full_points + batches_work(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

     if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  end do ! end loop over batches

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


  !-------shanghui begin parallel-----
  do i_coord = 1,3 
     do i_center = 1, n_centers_in_sc_DFPT 
  if(.not. use_local_index) call sync_sparse_matrix( first_order_S_sparse(i_coord,i_center,:) )
     enddo 
  enddo 
  !-------shanghui end parallel------

!
!  if (condition_penalty > 0.d0) then
!     if(use_batch_permutation > 0) call aims_stop('condition_penalty not implemented yet for load balancing')
!     call add_const_to_hamiltonian(hamiltonian, condition_penalty)
!  end if


  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)

  if (first_integration) then

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

        ! write integration grid per species in the format needed for control.in
        if (myid.eq.0) then
           write(use_unit,*)
           write(use_unit,'(2X,A,A)') "Output of integration grids in suitable form ", &
                "for copy-paste into control.in:"
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

     end if

  end if



  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( wave                 )) deallocate( wave                 )
  if(allocated( kinetic_wave         )) deallocate( kinetic_wave         )
  if(allocated( radial_wave_deriv    )) deallocate( radial_wave_deriv    )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( scaled_dylm_dphi_tab )) deallocate( scaled_dylm_dphi_tab )
  if(allocated( dylm_dtheta_tab      )) deallocate( dylm_dtheta_tab      )
  if(allocated( gradient_basis_wave  )) deallocate( gradient_basis_wave  )
  if(allocated( gradient_basis_wave_npoints)) deallocate( gradient_basis_wave_npoints)
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )

  if(supercell_fd_benchmark) then
!----------------------------shanghui debug DFPT_phonon: the supercell fd-benchmark------------------------ 
          n_basis_uc = 10 ! H-tier1 
          i_want = 0
 
          do i_cell_1 = 1, 1
          do i_basis_1 = 1,n_basis
 
             do i_center = index_hamiltonian(1,i_cell_1,i_basis_1), & 
                           index_hamiltonian(2,i_cell_1,i_basis_1)
                  if(column_index_hamiltonian(i_center).le.n_basis_uc ) then
                    
                     if( mod(i_basis_1-1,n_basis_uc)+1 .ge. column_index_hamiltonian(i_center)) then 
                     i_want = i_want + 1
 
                     write(info_str,'(a,2i5)') 'i_want:',i_want,i_center
                     call localorb_info(info_str, use_unit,'(A)', OL_norm)
 
                     write(info_str,'(a,f10.6)') 's:', overlap_matrix(i_center)
                     call localorb_info(info_str, use_unit,'(A)', OL_norm)
 
                     hamiltonian_i_want(i_want)= overlap_matrix(i_center)
                     endif
 
                  endif
             enddo 
 
          end do
          end do
 
   write(info_str,'(a,i5)') 'shanghui for python:',i_want
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
 
   write(info_str,'(30f20.15)') (hamiltonian_i_want(i_basis_1),i_basis_1=1,30)
   call localorb_info(info_str, use_unit,'(A)', OL_norm)
 
!----------------------------shanghui end debug DFPT_phonon: the supercell fd-benchmark------------------------ 
   endif
  if(use_batch_permutation > 0) then
    deallocate(ins_idx)
  endif

  if(get_batch_weights) deallocate(batch_times)



end subroutine integrate_first_order_S_p1
!******
