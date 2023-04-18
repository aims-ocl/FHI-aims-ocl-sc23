!****s* FHI-aims/integrate_ovlp_matrix_p2
!  NAME
!   integrate_ovlp_matrix_p2
!  SYNOPSIS

subroutine integrate_ovlp_matrix_p2( partition_tab_std, basis_l_max, &
     overlap_matrix )

!  PURPOSE
!  The subroutine integrates the overlap matrix
!  using a fixed basis set (no adaptive modifications).
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use basis
  use mpi_utilities
  use mpi_tasks
  use synchronize_mpi
  use localorb_io
  use constants
  use load_balancing
  use aims_memory_tracking
  use rel_x2c_mod
!----------------------------shanghui debug DFPT_phonon: the supercell fd-benchmark------------------------
  use pbc_lists
!----------------------------shanghui end debug DFPT_phonon: the supercell fd-benchmark------------------------
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points) :: partition_tab_std
  integer basis_l_max (n_species)
  ! when this routine is called, overlap_matrix has either the dimension
  ! (n_hamiltonian_matrix_size) or (n_local_matrix_size)
  ! so we declare it here as a 1D assumed size array
  real*8 :: overlap_matrix(*)


!  INPUTS
!  o partition_tab_std -- values of partition function
!  o basis_l_max -- maximum l component of basis functions
!
!  OUTPUT
!  o overlap_matrix -- overlap matrix
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

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8 coord_current(3)
  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_centers_integrals)
  real*8 dir_tab(3, n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_integrals)

  real*8 radial_wave(n_max_compute_fns_ham)
  real*8 wave(n_max_compute_ham, n_max_batch_size)

  !    Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:), allocatable :: matrix_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  !     for pruning of atoms, radial functions, and basis functions, to only the relevant ones ...

  integer :: n_compute_a, n_compute_c
  integer :: i_basis(n_centers_basis_I)

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_max_compute_atoms)
  integer :: spline_array_end(n_max_compute_atoms)

! VB - renewed index infrastructure starts here

      real*8 one_over_dist_tab(n_max_compute_atoms)

      ! indices for basis functions that are nonzero at current point
!
! VB: Warning - n_max_compute_fns_ham is outrightly wrong for the density
!     update!!
!
!     This should MUCH rather be counted in prune_basis for each batch and then
!     passed here, that would be the appropriate maximum!
!

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
  integer i_grid
  integer i_index, i_l, i_m
  integer i_coord
  integer i_center, i_center_L
  integer i_division

  integer i_species

  integer i_point
  integer :: i_full_points
  integer :: i_full_points_2

  integer :: i_my_batch

  character*200 :: info_str
  integer :: info

  ! Load balancing stuff

  integer n_my_batches_work ! Number of batches actually used
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches actually used

  integer dim_overlap_matrix ! actual dimension of overlap matrix
  integer ld_dirac_large ! actual dimension of overlap matrix (for SCALAR integration, 
                         ! instead of the FINAL form!) in fully-relativistic cases

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all

!----------------------------shanghui debug DFPT_phonon: the supercell fd-benchmark------------------------
  integer i_cell_1, i_want
  real*8  overlap_matrix_i_want(40)
!----------------------------shanghui end debug DFPT_phonon: the supercell fd-benchmark------------------------

  !  begin work

  if(use_batch_permutation > 0) then
    write(info_str,'(2X,A,A)')"Integrating overlap matrix with loadbalancing"
  else
    write(info_str,'(2X,A,A)')"Integrating overlap matrix."
  endif
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  !     begin with general allocations
  l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ), stat=info )
  call check_allocation(info,'ylm_tab','integrate_ovlp_matrix_p2',  (l_ylm_max+1)**2,  n_max_compute_atoms)
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), stat=info )
  call check_allocation(info,'index_lm','integrate_ovlp_matrix_p2')

  allocate ( matrix_shell(n_max_compute_ham*n_max_compute_ham), stat=info )
  call check_allocation(info,'matrix_shell','integrate_ovlp_matrix_p2')

  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
    ld_dirac_large = (n_centers_basis_I+1)*n_centers_basis_I/2
    call aims_allocate( dirac_s_sum, ld_dirac_large, 'dirac_s_sum')
    dirac_s_sum = 0.d0
  endif

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

    dim_overlap_matrix = batch_perm(n_bp)%n_local_matrix_size

  else

    n_my_batches_work = n_my_batches
    batches_work => batches
    partition_tab => partition_tab_std

    dim_overlap_matrix = n_hamiltonian_matrix_size

  endif

  if(get_batch_weights) allocate(batch_times(n_my_batches_work))

  !-----------------------------------------------------------------------------

  ! initialize

  overlap_matrix(1:dim_overlap_matrix) = 0.d0

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

              !     get current integration point coordinate
              coord_current(:) = batches_work(i_my_batch) % points(i_index) % coords(:)

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if

              !     compute atom-centered coordinates of current integration point,
              !     as viewed from all atoms
              call tab_atom_centered_coords_p0( coord_current, dist_tab_sq(1,i_point), &
                   dir_tab(1,1,i_point), n_centers_integrals, centers_basis_integrals )

              !    determine which basis functions are relevant at current integration point,
              !     and tabulate their indices

              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
              if (.not.prune_basis_once) then
                 call prune_basis_p2( dist_tab_sq(1,i_point), n_compute_c, &
                      i_basis, n_centers_basis_I, n_centers_integrals, &
                      inv_centers_basis_integrals )
              end if
           end if

        enddo ! end loop over one batch

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

           i_point = 0

           ! loop over one batch of integration points
           do i_index = 1, batches_work(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points = i_full_points + 1

              if (partition_tab(i_full_points).gt.0.d0) then

                 i_point = i_point+1

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
                     n_centers_integrals, centers_basis_integrals, &
                     n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
                     fn_atom, n_zero_compute, zero_index_point &
                    )


                 ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                 ! for all atoms which are actually relevant
                 call tab_local_geometry_p2 &
                      ( n_compute_atoms, atom_index, dist_tab(1,i_point), i_r )

                 ! compute trigonometric functions of spherical coordinate angles
                 ! of current integration point, viewed from all atoms
                 call tab_trigonom_p0 &
                      ( n_compute_atoms, dir_tab(1,1,i_point), trigonom_tab )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_p0 &
                      ( n_compute_atoms, atom_index, trigonom_tab, basis_l_max, l_ylm_max, ylm_tab )


                 ! Now evaluate radial functions
                 ! from the previously stored compressed spline arrays
                 call evaluate_radial_functions_p0  &
                      ( spline_array_start, spline_array_end, n_compute_atoms, n_compute_fns,&
                      dist_tab(1,i_point), i_r, atom_index, i_basis_fns_inv, basis_wave_ordered, radial_wave,  &
                      .false. , n_compute_c, n_max_compute_fns_ham )


                 ! tabulate total wave function value for each basis function
                 call evaluate_waves_p2  &
                   ( n_compute_c, n_compute_atoms, n_compute_fns, l_ylm_max, ylm_tab, one_over_dist_tab,   &
                     radial_wave(1), wave(1,i_point), rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point )

                 ! Reset i_basis_fns_inv
                 i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

              end if

              ! end loop over one batch
           enddo

           ! add full non-relativistic contributions and (for relativistic points)
           ! all contributions from the potential to the Hamiltonian matrix elements
           call evaluate_ovlp_shell_p0 ( n_points, partition, n_compute_c, n_compute_c, &
                wave(1,1), matrix_shell, n_max_compute_ham )

           if(use_batch_permutation > 0) then
              ! If use_batch_permutation > 0, the local overlap is always stored
              ! in full form for the local basis functions

              ! Get position of basis functions of current batch within local overlap
              do i=1,n_compute_c
                 ins_idx(i) = batch_perm(n_bp)%i_basis_glb_to_loc(i_basis(i))
              enddo

              ! Insert hamiltonian_shell of current batch
              do i=1,n_compute_c
                 i_off = (ins_idx(i)*(ins_idx(i)-1))/2
                 do j=1,i ! n_compute_c
                    overlap_matrix(ins_idx(j)+i_off) = overlap_matrix(ins_idx(j)+i_off) &
                                                     + matrix_shell(j+(i-1)*n_compute_c)
                 enddo
              enddo
           else
              if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
                call update_full_matrix_rel (n_compute_c, i_basis, ld_dirac_large, matrix_shell, dirac_s_sum)
              else
                call update_full_matrix_p0X ( n_compute_c, n_compute_c,  &
                   i_basis, matrix_shell, overlap_matrix )
              endif
           endif
        else
           i_full_points = i_full_points + batches_work(i_my_batch)%size
        end if ! n_compute.gt.0

        if(get_batch_weights) batch_times(i_my_batch) = mpi_wtime() - time_start

  enddo ! end loop over batches
 

  ! Get work time and total time after barrier
  time_work = mpi_wtime()-time0
  call mpi_barrier(mpi_comm_global,info)
  time_all = mpi_wtime()-time0
  call sync_real_number(time_work)
  call sync_real_number(time_all)
  write(info_str,'(a,2(f12.3,a))') '  Time summed over all CPUs for integration: real work ', &
     time_work,' s, elapsed ',time_all,' s'
  call localorb_info(info_str, use_unit, "(A)", OL_norm)

  !     synchronise the hamiltonian
  if(.not. use_local_index)then
     if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then
        call sync_vector(dirac_s_sum,ld_dirac_large)
     else
        call sync_integrate_ovlp( overlap_matrix )
     endif
  endif

  if(flag_rel.eq.REL_x2c .or. flag_rel.eq.REL_4c_dks)then ! Now, transform the obtained scalar integrations to spinor integrations
    write(info_str,'(2X,A)')'spinor S matrix'
    call localorb_info(info_str,use_unit,'(A)',OL_norm)
    call transform_scalar2spinor(3, 1, n_centers_basis_I, ld_dirac_large, dirac_s_sum, dirac_s)
  endif

  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)

  if (allocated(dirac_s_sum)) call aims_deallocate( dirac_s_sum, 'dirac_s_sum')

  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if

  if (allocated(matrix_shell)) then
     deallocate( matrix_shell )
  end if

  if (allocated(ins_idx)) then
     deallocate(ins_idx)
  endif

  if (allocated(batch_times)) then
     deallocate(batch_times)
  endif

end subroutine integrate_ovlp_matrix_p2

!----------------------------------------------------------------------
!******
