!****s* FHI-aims/integrate_ovlp_matrix_p1
!  NAME
!   integrate_ovlp_matrix_p1
!  SYNOPSIS

subroutine integrate_ovlp_matrix_p1( partition_tab, basis_l_max, &
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
  use synchronize_mpi
  use localorb_io
  use constants
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: partition_tab
  integer basis_l_max (n_species)
  real*8 :: overlap_matrix( n_hamiltonian_matrix_size )


!  INPUTS
!  o partition_tab -- values of partition function
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
  real*8, dimension(:,:), allocatable :: matrix_shell

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
  
  character*100 :: info_str

  !  begin work

  write(info_str,'(2X,A,A)')"Integrating overlap matrix."
  call localorb_info(info_str,use_unit,'(A)',OL_norm)

  !     begin with general allocations
  l_ylm_max = l_wave_max

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) )
  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max) ) 
  
  allocate ( matrix_shell(n_max_compute_ham,n_max_compute_ham) ) 

  !     initialize

  overlap_matrix = 0.d0

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

              !     get current integration point coordinate
              coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
           
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
                 call prune_basis_p0( dist_tab_sq(1,i_point), n_compute_a, n_compute_c, &
                      i_basis, n_centers_basis_I, n_centers_integrals, &
                      inv_centers_basis_integrals )
              end if
           end if

        enddo ! end loop over one batch

        if (prune_basis_once) then
           n_compute_a = batches(i_my_batch)%batch_n_compute
           n_compute_c = n_compute_a
           i_basis(1:n_compute_a) = batches(i_my_batch)%batch_i_basis
        end if

        n_points = i_point
     
        ! Perform actual integration if more than 0 basis functions
        ! are actually relevant on the present angular shell ...
        if (n_compute_a.gt.0) then

           i_point = 0
        
           ! loop over one batch of integration points
           do i_index = 1, batches(i_my_batch)%size, 1
           
              ! Increment the (global) counter for the grid, to access storage arrays
              i_full_points = i_full_points + 1

              if (partition_tab(i_full_points).gt.0.d0) then

                 i_point = i_point+1
               
                 ! for all integrations
                 partition(i_point) = partition_tab(i_full_points)
               
                 n_compute_atoms = 0
                 n_compute_fns = 0
                 !!! i_basis_fns_inv = 0

                 ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                 ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                 ! without any copying and without doing any unnecessary operations. 
                 ! The price is that the interface is no longer explicit in terms of physical 
                 ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.
                       
                 call prune_radial_basis_p0 &
                      ( dist_tab_sq(1,i_point), &
                      dist_tab(1,i_point), &
                      dir_tab(1,1,i_point), &
                      n_compute_atoms, atom_index, atom_index_inv, &
                      n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                      i_atom_fns, spline_array_start, spline_array_end, &
                      n_centers_integrals, centers_basis_integrals)

                 ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                 ! for all atoms which are actually relevant
                 call tab_local_geometry_p0 &
                      ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                      dir_tab(1,1,i_point), dist_tab(1,i_point),  &
                      i_r )
                       
                 !              compute trigonometric functions of spherical coordinate angles
                 !              of current integration point, viewed from all atoms
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
                 call evaluate_waves_p0  &
                      ( l_ylm_max,   &
                      ylm_tab, dist_tab(1,i_point),   &
                      index_lm, n_compute_c,   &
                      i_basis, radial_wave,   &
                      wave(1,i_point), n_compute_atoms,   &
                      atom_index_inv, n_compute_fns,  &
                      i_basis_fns_inv,  n_max_compute_fns_ham )

                 ! Reset i_basis_fns_inv
                 i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0
              
                 !            end if (partition_tab.gt.0)
              end if
              
              ! end loop over one batch
           enddo
        
           ! add full non-relativistic contributions and (for relativistic points)
           ! all contributions from the potential to the Hamiltonian matrix elements
           call evaluate_ovlp_shell_p0 &
                ( n_points,  &
                partition, &
                n_compute_a, n_compute_c, &
                wave(1,1), &
                matrix_shell, n_max_compute_ham  &
                )
           
           call update_full_matrix_p0 &
                ( n_compute_c, n_compute_a,  &
                i_basis, &
                matrix_shell,  &
                overlap_matrix &
                )
        else
           i_full_points = i_full_points + batches(i_my_batch)%size
           !      end if (n_compute.gt.0) then
        end if

        ! end distribution over batches
        !     end if
           
     !     end loop over bathces
  enddo

  !     synchronise the hamiltonian
  if(.not. use_local_index) call sync_integrate_ovlp( overlap_matrix )

  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if

  if (allocated(matrix_shell)) then
     deallocate( matrix_shell ) 
  end if





end subroutine integrate_ovlp_matrix_p1

!----------------------------------------------------------------------
!******     
