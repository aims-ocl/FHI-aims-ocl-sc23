!****s* FHI-aims/integrate_v_external_matrix_p2
!  NAME
!   integrate_v_external_matrix_p2
!  SYNOPSIS

subroutine integrate_v_external_matrix_p2 &
     (  &
     partition_tab, basis_l_max, hamiltonian &
     )

!  PURPOSE
!  Integrates the matrix elements ONLY of the external (nuclear) potential.
!  using a fixed basis set
!
!  WARNING: Does not work for periodic boundary conditions (an extra Ewald
!           compensating charge density would have to be added to neutralize the
!           otherwise infinite sum of nuclear potentials).
!           Requires special output for packed matrices!
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
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  implicit none

!  ARGUMENTS

  real*8, dimension( n_full_points)           :: partition_tab
  integer ::  basis_l_max (n_species)
  real*8  :: hamiltonian( n_hamiltonian_matrix_size, n_spin )

!  INPUTS
!  o partition_tab -- values of partition functions
!  o basis_l_max -- maximum l of basis functions.
!
!  OUTPUT
!  o hamiltonian -- Hamiltonian matrix
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2009), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
!  SOURCE

!  local variables

  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8 coord_current(3)

  real*8,dimension(:,:),allocatable:: dist_tab
  real*8,dimension(:,:),allocatable:: dist_tab_sq

  real*8 i_r(n_max_compute_atoms)
 
  real*8, dimension(:,:,:),allocatable:: dir_tab 
  
  real*8 trigonom_tab(4,n_max_compute_atoms)

  real*8, dimension(:), allocatable :: dist_tab_full
  real*8, dimension(:,:), allocatable :: dir_tab_full_norm
  real*8, dimension(:), allocatable :: i_r_full

  real*8,dimension(:,:,:),allocatable:: V_times_psi
  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:,:)  ,allocatable:: wave


  !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:,:), allocatable :: hamiltonian_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)

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
  character*100 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info

  !  begin work

  write(info_str,'(2X,A)') "Integrating matrix elements of electron-nuclear potential."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  !     begin with general allocations

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

     !       no gradients needed
     l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info ) 
  call check_allocation(info, 'index_lm                      ')

  call aims_allocate( hamiltonian_shell, n_max_compute_ham, n_max_compute_ham, "hamiltonian_shell" ) 
  call aims_allocate( V_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,      "V_times_psi" )

  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  call aims_allocate( wave, n_max_compute_ham, n_max_batch_size,                            "wave" )

  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')

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

  !     initialize

  hamiltonian = 0.d0

  i_basis_fns_inv = 0

  !     initialize index_lm
  i_index = 0
  do i_l = 0, l_wave_max, 1
     do i_m = -i_l, i_l
        i_index = i_index+1
        index_lm(i_m,i_l) = i_index
     enddo
  enddo

  i_full_points = 0
  i_full_points_2 = 0

  !     perform partitioned integration, batch by batch of integration point.
  !     This will be the outermost loop, to save evaluations of the potential.
  !     and the Y_lm functions

  do i_my_batch = 1, n_my_batches, 1

     n_compute_c = 0
     n_compute_a = 0
     i_basis = 0

     i_point = 0

     ! loop over one batch
     do i_index = 1, batches(i_my_batch)%size, 1

        i_full_points_2 = i_full_points_2 + 1


        if (partition_tab(i_full_points_2).gt.0.d0) then

           i_point = i_point+1

           ! get current integration point coordinate
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
                n_centers_integrals, centers_basis_integrals )

           ! determine which basis functions are relevant at current integration point,
           ! and tabulate their indices
           
           ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
           if (.not.prune_basis_once) then
              call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                   n_compute_c, i_basis,  &
                   n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
           end if

        end if
     enddo  ! end loop over one batch

     if (prune_basis_once) then
        n_compute_c = batches(i_my_batch)%batch_n_compute
        i_basis(1:n_compute_c) = batches(i_my_batch)%batch_i_basis
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
        do i_index = 1, batches(i_my_batch)%size, 1

           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1

               call tab_global_geometry_p0 &
                  ( dist_tab_sq(1,i_point), &
                    dir_tab(1,1,i_point), &
                    dist_tab_full, &
                    i_r_full, &
                    dir_tab_full_norm, &
                    n_centers_integrals,  centers_basis_integrals)


              ! for all integrations
              partition(i_point) = partition_tab(i_full_points)

              n_compute_atoms = 0
              n_compute_fns = 0
              
              ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
              ! are stored in a compact spline array that can be accessed by spline_vector_waves, 
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

                ! tabulate distance and Ylm's w.r.t. other atoms            
                call tab_wave_ylm_p0 &
                   ( n_compute_atoms, atom_index,  &
                   trigonom_tab, basis_l_max,  &
                   l_ylm_max, ylm_tab )
     
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

              i_spin = 1 !VB: If spin dependence needs to be restored one day, we'll see.

              !!! The LIST of atoms passed into this subroutine will NOT work for periodic boundary conditions!
              !!! ... even if an Ewald part were added. Must then also adjust the list of contributing atoms.
              call evaluate_nuclear_potential &
                   ( &
                     dist_tab_full,  &
                     n_centers_integrals, centers_basis_integrals,  &
                     local_potential_parts(i_spin) &
                   )

!              do i_spin = 1, n_spin, 1
                 call evaluate_V_psi_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab,  &
                   radial_wave, V_times_psi(1, i_point, i_spin),  &
                   local_potential_parts(i_spin),  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )
!             enddo

              ! Reset i_basis_fns_inv
              i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch
        
        ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
        ! work separately for each spin channel
        ! do i_spin = 1, n_spin, 1

           i_spin = 1 !VB: If spin dependence needs to be restored one day, we'll see.

           ! add full non-relativistic contributions and (for relativistic points)
           ! all contributions from the potential to the Hamiltonian matrix elements
           call evaluate_hamiltonian_shell_p1  &
                ( n_points, partition(1), n_compute_c, &
                V_times_psi(1,1,i_spin),  &
                n_max_compute_ham, wave(1,1),  &
                hamiltonian_shell )

           call update_full_matrix_p0(  &
                n_compute_c, n_compute_c, i_basis(1), hamiltonian_shell,    &
                hamiltonian(1,i_spin) )

        ! enddo

        ! Hamiltonian is now complete.

     else

       i_full_points = i_full_points + batches(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then

end do ! end loop over bathces


  !     synchronize the hamiltonian
  if(.not. use_local_index) call sync_integrate_hamiltonian(  hamiltonian )

  ! Allocatable arrays that are tracked
  if(allocated( hamiltonian_shell    )) call aims_deallocate( hamiltonian_shell, "hamiltonian_shell" )
  if(allocated( V_times_psi          )) call aims_deallocate( V_times_psi,             "V_times_psi" )
  if(allocated( wave                 )) call aims_deallocate( wave,                           "wave" )

  if(allocated( i_basis              )) deallocate( i_basis              )
  if(allocated( radial_wave          )) deallocate( radial_wave          )
  if(allocated( index_lm             )) deallocate( index_lm             )
  if(allocated( ylm_tab              )) deallocate( ylm_tab              )
  if(allocated( i_atom_fns           )) deallocate( i_atom_fns           )
  if(allocated( i_basis_fns_inv      )) deallocate( i_basis_fns_inv      )
  if(allocated( i_basis_fns          )) deallocate( i_basis_fns          )
  if(allocated( dir_tab              )) deallocate( dir_tab              )
  if(allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if(allocated( dist_tab             )) deallocate( dist_tab             )

  if(allocated( dist_tab_full        )) deallocate( dist_tab_full        )
  if(allocated( i_r_full             )) deallocate( i_r_full             )
  if(allocated( dir_tab_full_norm    )) deallocate( dir_tab_full_norm    )

end subroutine integrate_v_external_matrix_p2
!******
