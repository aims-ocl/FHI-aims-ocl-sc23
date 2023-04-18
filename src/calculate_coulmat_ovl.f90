!****h* FHI-aims/meta_gga_postp
!  NAME
!    meta_gga_postp - calculates meta_gga functioanls as post-processing
!  SYNOPSIS


module calculate_coulmat_ovl
!  PURPOSE
!  This module contains all routines related to the evaluation of the Coulomb-
!  matarix-elements as postprocessing partioned in its fourier coefficents, 
!  calculated from overlap integrals.
!
!  AUTHOR
!    Bjoern Bieniek
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

save
! the coulomb matrix as function of space C = <psi |exp(iq)|  psi>

complex*16, dimension(:,:,:), allocatable :: coulmat_full_w_k
complex*16, dimension(:,:), allocatable :: coulelement_k
contains
!******	

!!!!!!!!!!!! Allocation routines !!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine allocate_coulmat_w
!  PURPOSE
!    allocation of coulomb matrix fourier component, integrated, Fourier 
!    transformed, basis transformed
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  integer :: info
  if (.not.allocated(coulmat_full_w_k)) then
     allocate( coulmat_full_w_k ( n_k_points, n_k_points, n_basis*&
                                                    (n_basis+1)/2),stat=info)
     call check_allocation(info, 'coulmat_full_w_k             ')
  end if
end subroutine allocate_coulmat_w

subroutine allocate_coulelement(n_state_min_in, n_state_max_in)
!  PURPOSE
!    allocation of coulomb matrix fourier component, integrated, Fourier 
!    transformed, basis transformed
!  USES
  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: check_allocation
!  ARGUMENTS
  integer, Intent(in) :: n_state_min_in
  integer, Intent(in) :: n_state_max_in
!  INPUTS
! o n_state_min_in -- Minimum state concidered
! o n_state_max_in -- Maximum state concidered
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
  integer :: info

  if (.not.allocated(coulelement_k)) then
     allocate(coulelement_k(((n_state_max_in-n_state_min_in+1)+1)* &
                             (n_state_max_in-n_state_min_in+1)/2, &
                              n_k_points),stat=info)
     call check_allocation(info, 'coulelement_k                ')
  end if
end subroutine allocate_coulelement

subroutine clean_coulmat_w
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(coulmat_full_w_k))  deallocate(coulmat_full_w_k)

end subroutine clean_coulmat_w

subroutine clean_coulelement
!  PURPOSE
!    deallocation of local variables
!  USES
   use runtime_choices
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

  implicit none
  if (allocated(coulelement_k)) deallocate(coulelement_k)

end subroutine clean_coulelement

subroutine calculate_coulmat( partition_tab_std, basis_l_max,coulmat_full_w, &
                              q_vec )

!  PURPOSE
!  The subroutine integrates the fourier components of the coulomb matrix
!  using a fixed basis set (no adaptive modifications). The integrated matrix 
!  is fourier-transformed right away 
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
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, target, dimension(n_full_points) :: partition_tab_std
  integer basis_l_max (n_species)
  real*8 :: q_vec(3)
  complex*16::     coulmat_full_w(n_k_points, n_k_points, n_basis*(n_basis+1)/2)

!  INPUTS
!  o partition_tab_std -- values of partition function
!  o basis_l_max -- maximum l component of basis functions
!  o q_vec -- q-point vector
!  o k_point -- k-point
!  o k_strich_point -- k-strich-point
!
!  OUTPUT
!  o coulmat_full_w -- coulomb matrix fourier component
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
  real*8 coord_current_new(3)

  complex*16 :: coul(n_max_batch_size)
  real*8 :: q_vec_kspace(3)

  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_centers_integrals)
  real*8 dir_tab(3, n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_integrals)

  real*8 radial_wave(n_max_compute_fns_ham)
  real*8 wave(n_max_compute_ham, n_max_batch_size)

  !    Auxiliary Hamiltonian matrix, to sum up contributions from only a single 
  !    integration shell.
  complex*16, dimension(:), allocatable :: matrix_shell

  !     optimal accounting for matrix multiplications: only use points with 
  !     nonzero components
  integer :: n_points

  !     and condensed version of partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

  !     for pruning of atoms, radial functions, and basis functions, to only 
  !     the relevant ones ...

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
  type (batch_of_points), pointer :: batches_work(:) ! Pointer to batches 
                                                     ! actually used

  integer dim_overlap_matrix ! actual dimension of overlap matrix

  ! Pointers to the actually used array
  real*8, pointer :: partition_tab(:)

  ! Timing
  real*8, allocatable :: batch_times(:)
  real*8 time_start

  integer i, j, i_off, n_bp
  integer, allocatable :: ins_idx(:)

  ! Timings for analyzing work imbalance
  real*8 time0, time_work, time_all


  !  begin work


  !     begin with general allocations
  l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ), stat=info )
  call check_allocation(info,'ylm_tab','integrate_ovlp_matrix_p2', &
                        (l_ylm_max+1)**2,  n_max_compute_atoms)
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), stat=info )
  call check_allocation(info,'index_lm','integrate_ovlp_matrix_p2')

  allocate ( matrix_shell(n_max_compute_ham*n_max_compute_ham), stat=info )
  call check_allocation(info,'matrix_shell','integrate_ovlp_matrix_p2')

  !-----------------------------------------------------------------------------

  ! Initialize load balancing:
  ! Set pointers either to permuted batches / arrays over integration points 
  ! (for load balancing)
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

  coulmat_full_w = (0.d0,0.d0)
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

  coul = (1.d0, 0.d0)

  q_vec_kspace(1)= q_vec(1)*recip_lattice_vector(1,1)+q_vec(2)*&
                   recip_lattice_vector(2,1)+q_vec(3)*recip_lattice_vector(3,1)
  q_vec_kspace(2)= q_vec(1)*recip_lattice_vector(1,2)+q_vec(2)*&
                   recip_lattice_vector(2,2)+q_vec(3)*recip_lattice_vector(3,2)
  q_vec_kspace(3)= q_vec(1)*recip_lattice_vector(1,3)+q_vec(2)*&
                   recip_lattice_vector(2,3)+q_vec(3)*recip_lattice_vector(3,3)

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
              coord_current(:) = batches_work(i_my_batch) % points(i_index) & 
                                 % coords(:)

              if(n_periodic > 0)then
                 call map_to_center_cell(coord_current(1:3) )
              end if

              !     compute atom-centered coordinates of current integration 
              !     point, as viewed from all atoms
              call tab_atom_centered_coords_p0( coord_current, &
                   dist_tab_sq(1,i_point), dir_tab(1,1,i_point), &
                   n_centers_integrals, centers_basis_integrals )

              !    determine which basis functions are relevant at current 
              !    integration point, and tabulate their indices

              ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) 
              ! are actually needed
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

        ! from list of n_compute active basis functions in batch, collect all 
        ! atoms that are ever needed in batch.
        call collect_batch_centers_p2 &
        ( n_compute_c, i_basis, n_centers_basis_I, n_centers_integrals, &
          inv_centers_basis_integrals, n_batch_centers, batch_center &
        )

        n_points = i_point

        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_c.gt.0) then

           i_point = 0

           ! loop over one batch of integration points
           do i_index = 1, batches_work(i_my_batch)%size, 1

              ! Increment the (global) counter for the grid, to access storage 
              ! arrays
              i_full_points = i_full_points + 1

              if (partition_tab(i_full_points).gt.0.d0) then

                 i_point = i_point+1

                 ! for all integrations
                 partition(i_point) = partition_tab(i_full_points)

                 n_compute_atoms = 0
                 n_compute_fns = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if 
                 ! needed) Are stored in a compact spline array that can be 
                 ! accessed by spline_vector_waves, without any copying and 
                 ! without doing any unnecessary operations. The price is that 
                 ! the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the 
                 ! reorganized spline arrays.

                 call prune_radial_basis_p2 &
                   ( n_max_compute_atoms, n_max_compute_fns_ham, &
                     dist_tab_sq(1,i_point), dist_tab(1,i_point), &
                     dir_tab(1,1,i_point), n_compute_atoms, atom_index, &
                     atom_index_inv, n_compute_fns, i_basis_fns, &
                     i_basis_fns_inv, i_atom_fns, spline_array_start, &
                     spline_array_end, n_centers_integrals, &
                     centers_basis_integrals, n_compute_c, i_basis, &
                     n_batch_centers, batch_center, &
                     one_over_dist_tab, rad_index, wave_index, l_index, &
                     l_count, fn_atom, n_zero_compute, zero_index_point &
                    )


                 ! Tabulate distances, unit vectors, and inverse logarithmic 
                 ! grid units for all atoms which are actually relevant
                 call tab_local_geometry_p2 &
                      ( n_compute_atoms, atom_index, &
                        dist_tab(1,i_point),  &
                      i_r )

                 ! compute trigonometric functions of spherical coordinate 
                 ! angles of current integration point, viewed from all atoms
                 call tab_trigonom_p0 &
                      ( n_compute_atoms, dir_tab(1,1,i_point),  &
                      trigonom_tab )

                 ! tabulate distance and Ylm's w.r.t. other atoms
                 call tab_wave_ylm_p0 &
                      ( n_compute_atoms, atom_index,  &
                      trigonom_tab, basis_l_max,  &
                      l_ylm_max, &
                      ylm_tab )


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
                     radial_wave(1), wave(1,i_point), &
                     rad_index, wave_index, l_index, l_count, fn_atom, &
                     n_zero_compute, zero_index_point &
                   )

		 coord_current_new(:) = batches_work(i_my_batch) % & 
                                        points(i_index) % coords(:)
		 if(n_periodic > 0)then
		   call map_to_center_cell(coord_current_new(1:3) )
		 end if
                 coul(i_point)=exp((0.d0,1.d0)*&
                         dot_product(coord_current_new(1:3),q_vec_kspace(1:3)))
!		 call get_phases &
!		  ( n_compute_c, n_compute_atoms, n_compute_fns, &
!		    dir_tab(1,1,i_point), dist_tab(1,i_point), &
!                    coord_current_new, q_vec_kspace, coul(i_point), &
!		    wave_index, l_count, fn_atom, &
!		    n_zero_compute, zero_index_point )
                 ! Reset i_basis_fns_inv
                 i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
              end if

              ! end loop over one batch
           enddo

           ! add full non-relativistic contributions and (for relativistic 
           ! points) all contributions from the potential to the Hamiltonian 
           ! matrix elements

	   call evaluate_coulmat( & 
		  n_points, partition(1), n_compute_c, coul, &
		  n_max_compute_ham, wave(1,1),matrix_shell)

           call update_full_matrix_and_ft &
                   ( n_compute_c, n_compute_c,  &
                   i_basis, &
                   matrix_shell,  &
                   coulmat_full_w)
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

  !     synchronise the hamiltonian
  if(.not. use_local_index) call sync_vector_complex( coulmat_full_w , &
                                n_k_points*n_k_points*n_basis*(n_basis+1)/2 )
  if(get_batch_weights) call set_batch_weights(n_bp, batch_times)


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

end subroutine calculate_coulmat
subroutine check_q_point_list(n_q_point_in)

  !  PURPOSE
  !  Checks if file 'q_list.in' exists and reads number of q points
  !
  !  USES

  use dimensions
  use runtime_choices
  use mpi_tasks
  use localorb_io
  implicit none

  !  ARGUMENTS

  integer, Intent(out) :: n_q_point_in

  ! INPUTS
  ! o file 'q_list.in'
  !  OUTPUT
  ! o n_q_point_in -- number of q points in file
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


  !  begin work
  integer :: i_code
  character*40 desc_str
  character(*), parameter :: func = 'check_q_point_list'

  if (myid.eq.0) then
    write(use_unit,'(2X,A)') &
    "q_point_external: Expecting the q points in file 'q_list.in' ."
  end if

    ! check if q_list.in exists
    open(888,file='q_list.in')

    read(888,*,iostat=i_code) desc_str
    if (i_code.ne.0) then
       if (myid.eq.0) then
           write(use_unit,'(1X,A)') "* Empty input file 'q_list.in'. "
           write(use_unit,'(1X,A)') "* Needed for 'q_point_external' ."
       end if
       call aims_stop_coll('', func)
    end if

    read(888,*) n_q_point_in
    if (myid.eq.0) then
       write(use_unit,'(2X,A,I8,A)') &
        "Header of 'q_list.in' anticipates ", n_q_point_in, " q points."
    end if

    close(888)

    if (n_q_point_in.le.0) then

      if (myid.eq.0) then
         write(use_unit,'(1X,A)') "* Error: Less than 1 q point requested."
         write(use_unit,'(1X,A)') "* This is unphysical, please correct."
      end if
      call aims_stop_coll('', func)

    end if

end subroutine check_q_point_list

subroutine get_q_point_list(q_point_list_in, n_q_point_in)

!  PURPOSE
!  Read q-point list from file 'q_list.in'
!
!  USES

  use dimensions
  use runtime_choices
  use mpi_tasks
  use localorb_io
  implicit none

!  ARGUMENTS

  real*8, Intent(out) :: q_point_list_in(n_q_point_in, 3)
  integer, Intent(inout) :: n_q_point_in

  ! INPUTS
  ! o file 'q_list.in'
  ! o n_q_point_in -- number of q points in file
  !  OUTPUT
  ! o n_q_point_in -- number of q points in file
  ! o q_point_list_in - List of q-points in units or reciprocal lattice vectors
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


  !  begin work
  real*8, allocatable :: q_weight(:)
  Integer :: n_q_point_xyz(3)
  Integer :: i_q_point
  real*8 :: r_x, r_y, r_z
  integer :: info
  character*140 :: info_str
  character(*), parameter :: func = 'get_q_point_list'

  write(info_str,'(2X,A)')  'External q-points from file q_list.in:'
  call localorb_info(info_str,use_unit,'(A)') 

  open(888,file='q_list.in')

  read(888,*)  n_q_point_xyz(1), n_q_point_xyz(2), n_q_point_xyz(3)

  write(info_str,'(2X,A,3I10)')'| k-points in every direction:',  &
  n_q_point_xyz(1), n_q_point_xyz(2), n_q_point_xyz(3)
  call localorb_info(info_str,use_unit,'(A)') 

  read(888,*)  n_q_point_in

  write(info_str,'(2X,A,I10)') '| Number of q-points:', n_q_point_in
  call localorb_info(info_str,use_unit,'(A)') 

  allocate(q_weight(n_q_point_in),stat=info)
  call check_allocation(info, 'q_weight', func)

  q_point_list_in(:,:) = 0.d0
  do  i_q_point = 1,  n_q_point_in

  read(888,*)  r_x, r_y, r_z,  q_weight(i_q_point)

  q_point_list_in (i_q_point,1) = r_x
  q_point_list_in (i_q_point,2) = r_y
  q_point_list_in (i_q_point,3) = r_z

  end do  ! i_q_point


  close(888)
  real_eigenvectors = .false.

  if(abs(sum(q_weight)-1) >1e-5)then
    write(use_unit,*) 'Error: sum of q-vector weights is not one!', &
    sum(q_weight)
    stop
  end if

  ! There can be numerical inaccuracy. This takes care of them.
  q_weight = q_weight/  sum(q_weight)

end subroutine get_q_point_list

subroutine get_phases &
      ( n_compute, n_compute_atoms, n_compute_fns, &
        dir_tab, dist_tab, coord_current_new, q_vec_kspace, coul_phase, &
        wave_index, l_count, fn_atom, &
        n_zero_compute, zero_index_point )

!  PURPOSE
!     Get the fourier component of the coulomb operater at current point
!
!  USES

      implicit none

!  ARGUMENTS


      integer, intent(IN) :: n_compute_atoms
      integer, intent(IN) :: n_compute_fns
      integer, intent(IN) :: n_compute

      real*8, intent(IN) :: dir_tab(3,n_compute_atoms)
      real*8, intent(IN) :: dist_tab(n_compute_atoms)
      real*8, intent(IN) :: q_vec_kspace(3)
      real*8, intent(IN) :: coord_current_new(3)
      integer, intent(IN) :: wave_index(n_compute_fns)
      integer, intent(IN) :: l_count(n_compute_fns)
      integer, intent(IN) :: fn_atom(n_compute_fns)
      integer, intent(IN) :: n_zero_compute
      integer, intent(IN) :: zero_index_point(n_compute)
      complex*16, intent(OUT)  :: coul_phase

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o n_compute_fns -- number of non-zero basis fns
!    o n_compute -- number of non-zero basis functions 
!    o dir_tab -- r/|r|
!    o dist_tab -- |r|
!    o q_vec_kspace -- q-point vector
!    o coord_current_new -- current point on grid
!
!    indices for basis functions that are nonzero at current point
!     o l_count
!     o fn_atom
!    indices for known zero basis functions at current point
!     o n_zero_compute
!     o zero_index_point
!
!  OUTPUT
!   o coul_phase -- fourier component at curremt point
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


!     counters

      integer :: i_compute_point

      integer :: i_compute
      integer :: i_compute_fn



!     begin work
!    coul_phase = (1.0d0, 0.0d0)
    coul_phase=exp((0.d0,1.d0)*&
                         dot_product(coord_current_new(1:3),q_vec_kspace(1:3)))

!   Atom centered coordinates
!     first, the nonzero functions
!     do i_compute_fn = 1, n_compute_fns, 1
!       coul_phase ( wave_index(i_compute_fn):wave_index(i_compute_fn)+&
!                     l_count(i_compute_fn) ) = exp((0.d0,1.d0)* &
!                      dot_product(dist_tab(fn_atom(i_compute_fn)) * &
!                      dir_tab(1:3,fn_atom(i_compute_fn)),q_vec_kspace(1:3)))
!     enddo
      !coul_phase(1:n_compute_fns)=exp((0.d0,1.d0)*&
      !                  dot_product(coord_current_new(1:3),q_vec_kspace(1:3)))
      ! then, the zero functions
      !do i_compute_point = 1, n_zero_compute, 1
      !  i_compute = zero_index_point(i_compute_point)

      !  coul_phase(i_compute) = (0.0d0, 0.0d0)

      !enddo


end subroutine get_phases

subroutine evaluate_coulmat( & 
     n_points, partition, n_compute, coul_in, &
     n_basis_list, wave, coulmat_shell_out)

!  PURPOSE
!  Evaluates coulmat for one grid batch,
!
!  USES

  use dimensions
  use runtime_choices
  implicit none

!  ARGUMENTS

  integer, Intent(in) :: n_basis_list
  integer, Intent(in) :: n_points
  integer, Intent(in) :: n_compute
  real*8, Intent(in)  :: partition(n_points)
  real*8, Intent(in)  :: wave(n_basis_list, n_points)
  complex*16, Intent(in)  :: coul_in(n_points)
  complex*16, Intent(out)  :: coulmat_shell_out( n_compute, n_compute)

! INPUTS
! o n_points -- number of points in this batch
! o partition -- partition table
! o n_compute -- number of non zero basis functions in this batch
! o coul_in -- Fourier component of Coulomb operator at points in batch
! o n_basis_list -- basis function list
! o wave -- wavefunctions at points in batch
!  OUTPUT
! o coulmat_shell_out - coulmat integrated over current batch
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


  !  begin work
      complex*16 :: wave_compute_a(n_compute,n_points)
      complex*16 :: wave_compute_b(n_compute,n_points)
!     counters

      integer :: i_compute, i_point, i, j
!     begin work

      coulmat_shell_out = (0.0d0,0.0d0)
      wave_compute_a = (0.0d0,0.0d0)
      wave_compute_b = (0.0d0,0.0d0)
!     Condense basis functions to only those that are used at present 
!     integration points

      do i_point = 1, n_points, 1
         wave_compute_a(1:n_compute, i_point) =  &
               (partition(i_point))*coul_in(i_point)*&
                wave(1:n_compute, i_point)  
         wave_compute_b(1:n_compute, i_point) =  &
                wave(1:n_compute, i_point) 
      enddo
      call zgemm('N', 'T', n_compute, n_compute, n_points, (1.0d0,0.0d0), &
               wave_compute_a, n_compute, wave_compute_b, n_compute, &
               (0.0d0,0.0d0), coulmat_shell_out, &
               n_compute ) 

!       call zsyr2k('U', 'N', n_compute, n_points, (0.5d0,0.0d0),  &
!                  wave_compute_a, n_compute, &
!                  (1.0, 0.0)*wave(1:n_compute, 1:n_points), n_compute, &
!                  (0.0d0,0.0d0), coulmat_shell, n_compute)
!       call zsyrk('U','N', n_compute,n_points, (1.0d0,0.0d0), &
!                  wave_compute_a, n_compute,  (0.0d0,0.0d0), coulmat_shell, &
!                  n_compute )
!write(use_unit,*)  wave_compute_a(5,10), wave(5,10), partition(10), coul_in(10)
end subroutine evaluate_coulmat

subroutine update_full_matrix_and_ft &
     ( n_compute_c, n_compute_a, i_basis, &
     matrix_shell,coulmat_full_w &
     )

!  PURPOSE
!  Subroutine adds a part of the integrals in a matrix and fourier transforms it
!  (only for the n_compute basis functions that are nonzero at the
!  current integration shell) to the full matrix (which contains
!  all n_basis basis functions). The link between i_compute = 1 ... n_compute
!  and i_basis = 1 ... n_basis is provided by the index array 
!  i_basis(i_compute).
!
!  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  use mpi_tasks, only: aims_stop, check_allocation
  use localorb_io, only: localorb_info
  implicit none

!  ARGUMENTS

  integer n_compute_c
  integer n_compute_a
  integer i_basis(n_compute_c)
  complex*16 matrix_shell(n_compute_c,n_compute_a)
  complex*16:: coulmat_full_w(n_k_points,n_k_points,n_basis*(n_basis+1)/2)
!  INPUTS
!   o n_compute_c == n_compute_a -- number of relevant basis functions
!   o i_basis -- list of relevant basis functions
!   o matrix_shell -- hamiltonian / overlap_matrix of relevant basis functions
!
!  OUTPUT
!   o coulmat_full_w -- date from matrix_shell is added here
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

  !     auxiliary matrices for Level 3 Blas matrix multiplications

  !     counters

  integer :: i_compute_1
  integer :: i_compute_2
  integer :: i_offset
  integer :: i_index_real
  integer :: i_offset_first_part
  integer :: i_cell_index, i_cell_1
  integer :: i_max_basis, i_min_basis
  integer :: i_one_part, i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell,&
             i_index
  integer :: offset(n_cells) !_in_hamiltonian)
  integer :: offset_end(n_cells) !_in_hamiltonian)
  integer :: help
  integer :: i_k, j_k
!  integer :: direct_citing(n_compute_a,n_compute_a )
!  real*8 :: data(n_cells_in_hamiltonian, n_basis)
!  real*8,dimension(:,:), allocatable :: data
  integer:: i_cell_old

!NEC_CB
  integer::help1(n_compute_a)
  integer::help2(n_compute_a)

  integer,dimension(:,:),allocatable :: index_b
  !     begin work

  !      now add the aux. ham. matrix to the actual ham. matrix
  !      this requires translating between the actually computed matrix elements
  !      and the full ham. matrix ...

  ! When basis is smaller than n_basis

  !      write(use_unit,*) n_compute_c,n_compute_a

  select case(packed_matrix_format)

  case(PM_none)

     i_index_real = 0
     do i_compute_2 = 1, n_compute_a, 1

        i_offset = (i_basis(i_compute_2)-1)*i_basis(i_compute_2)/2


        do i_compute_1 = 1,i_compute_2,1


           i_index_real = i_offset + i_basis(i_compute_1) 

           !matrix(i_index_real) = matrix(i_index_real) &
           !     + matrix_shell(i_compute_1, i_compute_2)


        enddo
     enddo


  case(PM_index) !--------------------------------------------------------------

     if(n_periodic == 0)then

        do i_compute_1 = 1, n_compute_a, 1

           
!           write(use_unit,*) 'i_compute_2', i_compute_2
           
           i_start =  index_hamiltonian(1,1, i_basis(i_compute_1))
           i_end   =  index_hamiltonian(2,1, i_basis(i_compute_1))


!           write(use_unit,*) '-'

           do i_compute_2 = 1,i_compute_1,1
              
              i_basis_2 = i_basis(i_compute_2)
              
              place: do i_place = i_start, i_end, 1
                 
                 if( column_index_hamiltonian( i_place) == i_basis_2)then
                       
!NEC_CB                    matrix(i_place) = matrix(i_place) + & 
!                                         matrix_shell(i_compute_2, i_compute_1)
                    if (i_compute_2.le.i_compute_1) then
                       coulmat_full_w(1,1,index_b(i_basis_2, i_basis_1))= &
                       coulmat_full_w(1,1,index_b(i_basis_2, i_basis_1)) +&
                                (1.0,0.0)*matrix_shell(i_compute_2, i_compute_1) 
                    else
                       coulmat_full_w(1,1,index_b(i_basis_2, i_basis_1))= &
                       coulmat_full_w(1,1,index_b(i_basis_2, i_basis_1)) +& 
                                (1.0,0.0)*matrix_shell(i_compute_1, i_compute_2)
                    end if

                    i_index_real = i_place
                    exit place 

                 
                 else if(column_index_hamiltonian( i_place) > i_basis_2)then
                    i_index_real = i_place
                    exit place  

                 end if
              end do place
              i_start = i_index_real
           end do
        end do
        


           
        
     else ! Periodic case----------------

     allocate(index_b(n_basis,n_basis), stat=i_index)
     call check_allocation(i_index, 'index_b                       ')
     index_b = 0
     i_index = 0
     do i_basis_2 = 1,n_basis, 1
        do i_basis_1 = 1,i_basis_2,1
           i_index = i_index + 1
           index_b(i_basis_1, i_basis_2) = i_index 
        end do
     end do
        ! Unfortunately the periodic systems can not use the searching routine 
        ! used now in the clusters. This is because the peridic systems have 
        ! extra packing for supercell information.

        do i_compute_1 = 1, n_compute_a, 1
          help1(i_compute_1)=Cbasis_to_basis(i_basis(i_compute_1))
          help2(i_compute_1)=&
                      center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))
        end do

        do i_compute_1 = 1, n_compute_a, 1

           i_basis_1 = help1(i_compute_1) !Cbasis_to_basis(i_basis(i_compute_1))
           i_cell_old = help2(i_compute_1)
                         !center_to_cell(Cbasis_to_center(i_basis(i_compute_1)))


           offset_end = -1
           offset = -1

           do i_cell_1 = 1, n_cells 


              i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

              offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

           end do



           do i_compute_2 = 1, n_compute_a
              
              i_basis_2 = help1(i_compute_2)
                                       !Cbasis_to_basis(i_basis(i_compute_2))



              if(i_basis_2 <= i_basis_1)then

                 i_cell    =  help2(i_compute_2)
                       !center_to_cell(Cbasis_to_center(i_basis(i_compute_2)))


                 place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

!                    write(use_unit,*) offset(i_cell), offset_end(i_cell)
                 
                    if( column_index_hamiltonian( i_place) == i_basis_2)then
                      do i_k=1, n_k_points,1
                        do j_k=1, n_k_points,1
                          if (i_compute_2.le.i_compute_1) then
                            coulmat_full_w(i_k, j_k, index_b(i_basis_2, &
                                           i_basis_1))=&
                            coulmat_full_w(i_k, j_k, index_b(i_basis_2, &
                                           i_basis_1)) +&
                                           (k_phase(i_cell_old,i_k)*&
                                           conjg(k_phase(i_cell,j_k)))&
                               *(1.0,0.0)*matrix_shell(i_compute_2, i_compute_1)
                           else
                             coulmat_full_w(i_k, j_k, index_b(i_basis_2, &
                                            i_basis_1))=&
                             coulmat_full_w(i_k, j_k, index_b(i_basis_2, &
                                            i_basis_1)) +&
                                            (k_phase(i_cell_old,i_k)*&
                                            conjg(k_phase(i_cell,j_k)))&
                               *(1.0,0.0)*matrix_shell(i_compute_1, i_compute_2)
                           end if
                         enddo
                       enddo
                       exit place_2
                 
                    else if(column_index_hamiltonian( i_place) > i_basis_2)then
                       exit place_2  


                    end if
                 end do place_2

                 offset(i_cell) = i_place

              end if
           end do
        end do

     end if ! n_periodic == 0 - else


  case default

     call localorb_info('Invalid packing')
     call aims_stop

  end select

end subroutine update_full_matrix_and_ft

subroutine calc_coulelement(coulelement, coulmat_full_w, KS_vec_one, &
                            KS_vec_complex_one, KS_vec_two, KS_vec_complex_two,&
                            n_state_min_in, n_state_max_in)


  !  PURPOSE
  !   Calculete Coulomb-Matrix fourier component elements by summing up of 
  !   KS EV coefficents
  !
  ! USES
  use pbc_lists
  use geometry
  use dimensions
  use runtime_choices
  use basis
  implicit none
  !  ARGUMENTS

  complex*16, INTENT(IN)::     KS_vec_complex_one(n_basis, &
                                     n_state_max_in-n_state_min_in+1, n_spin)
  real*8, INTENT(IN)::     KS_vec_one(n_basis, &
                                      n_state_max_in-n_state_min_in+1, n_spin)
  complex*16, INTENT(IN)::     KS_vec_complex_two(n_basis, &
                                      n_state_max_in-n_state_min_in+1, n_spin)
  real*8, INTENT(IN)::     KS_vec_two(n_basis, &
                                      n_state_max_in-n_state_min_in+1, n_spin)
  complex*16, INTENT(IN)::     coulmat_full_w(n_basis*(n_basis+1)/2)
  complex*16, INTENT(OUT) :: coulelement(((n_state_max_in-n_state_min_in+1)+1)*&
                                          (n_state_max_in-n_state_min_in+1)/2)
  integer, INTENT(IN):: n_state_min_in
  integer, INTENT(IN):: n_state_max_in
  !  INPUTS
  !    o KS_ve_complex_one -- KS coefficents at k-point, complex
  !    o KS_vec_one -- KS coefficents at k-point, real
  !    o KS_ve_complex_two -- KS coefficents at k_strich-point, complex
  !    o KS_vec_two -- KS coefficents at k_strich-point, real
  !    o coulmat_full_w -- Coulomb-Matrix fourier component, real spaces basis
  !    o n_state_min_in -- Minimum state concidered
  !    o n_state_max_in -- Maximum state concidered
  !  OUTPUT
  !    o coulelement -- Coulomb-Matrix fourier component, KS basis
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  ! SOURCE

  integer:: i_spin,i_basis_2, i_basis_1, n_state, m_state, i_basis
  integer::  num, num_base
     i_spin = 1
     num = 0
     coulelement = (0.0,0.0)
     do n_state = 1, n_state_max_in-n_state_min_in+1 !n_state_min_in, 
                                                     !n_state_max_in
       do m_state = n_state, n_state_max_in-n_state_min_in+1 !n_state_max_in
          num = num + 1
          coulelement(num) = (0.0,0.0)
          num_base = 0
        do i_basis_2 = 1, n_basis
           do i_basis_1 = 1, i_basis_2
                 num_base = num_base + 1
                 if(real_eigenvectors)then
                    if (i_basis_1==i_basis_2)then
                       coulelement(num) = coulelement(num) &
                            + (KS_vec_one(i_basis_1,n_state,i_spin)) &
			      *coulmat_full_w(num_base)*&
                               KS_vec_two(i_basis_2,m_state,i_spin) 
            	    else
                       coulelement(num) = coulelement(num) &
                            + (KS_vec_one(i_basis_1,n_state,i_spin)) &
			      *coulmat_full_w(num_base)*&
                               KS_vec_two(i_basis_2,m_state,i_spin) &
                            + conjg((KS_vec_one(i_basis_1,n_state,i_spin)) &
			      *coulmat_full_w(num_base)*&
                               KS_vec_two(i_basis_2,m_state,i_spin))
                    endif
                 else
                    if (i_basis_1==i_basis_2)then
                      coulelement(num) = coulelement(num)&
                          +conjg(KS_vec_complex_one(i_basis_1,n_state,i_spin)) &
  		          *coulmat_full_w(num_base)*&
                          (KS_vec_complex_two(i_basis_2,m_state,i_spin))
                    else
                      coulelement(num) = coulelement(num)+&
                           conjg(KS_vec_complex_one(i_basis_1,n_state,i_spin)) &
  		           *coulmat_full_w(num_base)*&
                           (KS_vec_complex_two(i_basis_2,m_state,i_spin))+&
                          ((conjg(KS_vec_complex_two(i_basis_2,n_state,i_spin))&
			   *conjg(coulmat_full_w(num_base))*&
                           (KS_vec_complex_one(i_basis_1,m_state,i_spin))))
                    endif
                 end if
           end do
         end do
        end do
      end do
end subroutine calc_coulelement

end module calculate_coulmat_ovl

