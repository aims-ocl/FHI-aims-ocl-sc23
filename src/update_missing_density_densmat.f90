!****s* FHI-aims/update_missing_density_densmat
!  NAME
!    update_missing_density_densmat
!  SYNOPSIS

subroutine update_missing_density_densmat &
     ( KS_eigenvector, KS_eigenvector_complex, occ_numbers, partition_tab,  &
     hartree_partition_tab, rho, basis_l_max, &
     density_matrix_sparse )

!  PURPOSE
!  Subroutine update_missing_density obtains the KS density from the eigenvectors
!  and occupation numbers, but ONLY on those grid points where the partition_tab is
!  normally zero
!
!  Notice that the output density is an _unmixed_ electron density, i.e., it is only 
!  consistent with the remaining density in the fully self-consistent case. 
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use species_data
  use mpi_utilities
  use localorb_io
  use basis
  use constants
  use analyze_arrays
  use scalapack_wrapper
  use KH_core_states
  use density_matrix_evaluation
  use sym_base, only: evaluate_densmat_sym
  implicit none

!  ARGUMENTS

  real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task):: KS_eigenvector_complex
  real*8,     dimension(n_states, n_spin, n_k_points)             :: occ_numbers

  real*8,  dimension(n_full_points)             ::  partition_tab
  real*8,  dimension(n_full_points)             ::  hartree_partition_tab
  real*8,  dimension(n_spin, n_full_points)     :: rho
  real*8,  dimension(3, n_spin, n_full_points)  :: rho_gradient
  integer, dimension(n_species)                 :: basis_l_max
  real*8,  dimension(n_spin)                    :: rho_change
  real*8,  dimension(n_hamiltonian_matrix_size) :: density_matrix_sparse

! INPUTS
! o KS_eigenvector -- eigenvectors if real eigenvectors are in use
! o KS_eigenvector_complex -- eigenvectors is complex eigenvectors are in use
! o occ_numbers -- occupation numbers of states
! o partition_tab -- partition function values
! o hartree_partition_tab -- hartree potentials partition function values
! o rho -- electron density, rho is only input; what we store is the density residual (i.e. the change in the density) 
! o rho_gradient -- gradient of the electron density 
! o basis_l_max -- maximum basis l value of the basis functions
! o density_matrix_sparse -- this is the works space for density matrix, it is used if packed matrixes are in use.
!
! OUTPUT
! o rho_change -- how much the electron density changed during the update.
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
  real*8,dimension(:,:),allocatable:: wave

  integer,dimension(:),allocatable :: i_basis
  integer :: n_compute_c

  real*8 coord_current(3)
  real*8, dimension(n_species) :: r_grid_min_sq

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_integrals)
  integer :: spline_array_end(n_centers_integrals)

! VB - renewed index infrastructure starts here

      real*8 :: one_over_dist_tab(n_max_compute_missing_atoms)

      ! indices for basis functions that are nonzero at current point

      integer :: rad_index(n_max_compute_missing_atoms)
      integer :: wave_index(n_max_compute_fns_missing_dens)
      integer :: l_index(n_max_compute_fns_missing_dens)
      integer :: l_count(n_max_compute_fns_missing_dens)
      integer :: fn_atom(n_max_compute_fns_missing_dens)

      ! indices for known zero basis functions at current point
      integer :: n_zero_compute
      integer :: zero_index_point(n_max_compute_missing_dens)

      ! active atoms in current batch
      integer :: n_batch_centers
      integer :: batch_center(n_centers_integrals)

  !     other local variables

  integer :: l_ylm_max
  integer :: n_points

  logical :: point_on_atom
  integer :: i_center

  real*8 :: temp_rho(n_max_batch_size,n_spin)
  real*8,dimension(:,:),allocatable :: temp_rho_small

  integer, dimension(:,:), allocatable :: index_lm
  real*8,  dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:),  allocatable :: density_matrix
  real*8, dimension(:,:),  allocatable :: density_matrix_con
  real*8, dimension(:,:),  allocatable :: work
  complex*16, dimension(:,:),allocatable :: work_complex

  !     counters
  integer :: i_l
  integer :: i_m
  integer :: i_coord
  integer :: i_state
  integer :: i_point
  integer :: i_bas, i_bas1, i_bas2

  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_3

  integer :: i_spin,  i_spin_2,  i_cell
  integer :: i_index, i_k_point, i_k
  integer :: i_my_batch
  integer :: info
  integer :: i_index_real, i_size


  call localorb_info( &
       "Evaluating KS electron density on previously omitted grid points, using the density matrix.", use_unit,'(2X,A)',OL_norm )

  if(use_small_component)then
     allocate(temp_rho_small(n_max_batch_size,n_spin),stat=info)
     call check_allocation(info, 'temp_rho_small                ') 
  end if

  if(packed_matrix_format == PM_none)then
     if(.not. allocated(density_matrix))then
        allocate(density_matrix(n_centers_basis_T, n_centers_basis_T),stat=info)
        call check_allocation(info, 'density_matrix                ') 
     end if
  else
     ! Allocate a dummy density matrix since otherwise the code will crash
     ! at the call of evaluate_densmat if checking is enabled
     if(.not. allocated(density_matrix)) allocate(density_matrix(1,1))
  end if

  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  do i_spin = 1, n_spin
  
     if(use_symmetry_reduced_spg)then
          call evaluate_densmat_sym(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
          &                     density_matrix, density_matrix_sparse, i_spin, .false.)
     else
          call evaluate_densmat(KS_eigenvector, KS_eigenvector_complex, occ_numbers,&
          &                     density_matrix, density_matrix_sparse, i_spin, .false.)
     endif

     if(.not. allocated(dist_tab))then
        allocate(dist_tab(n_centers_integrals),stat=info)
        call check_allocation(info, 'dist_tab                      ') 
     end if
     if(.not. allocated(dist_tab_sq))  then
        allocate(dist_tab_sq(n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dist_tab_sq                   ')
     end if
     if(.not. allocated(i_r))then
        allocate(i_r(n_max_compute_missing_atoms),stat=info)
        call check_allocation(info, 'i_r                           ')
     end if
     if(.not. allocated(dir_tab)) then
        allocate(dir_tab(3, n_centers_integrals, n_max_batch_size),stat=info)
        call check_allocation(info, 'dir_tab                       ')
     end if
     if(.not. allocated(trigonom_tab)) then
        allocate(trigonom_tab(4, n_max_compute_missing_atoms),stat=info)
        call check_allocation(info, 'trigonom_tab                  ')
     end if
     if(.not. allocated(density_matrix_con))then
        allocate(density_matrix_con(n_max_compute_missing_dens, n_max_compute_missing_dens),stat=info)
        call check_allocation(info, 'density_matrix_con            ')
     end if
     if(.not. allocated(work))then
        allocate(work(n_max_compute_missing_dens, n_max_batch_size),stat=info)
        call check_allocation(info, 'work                          ')
     end if
     if(.not. allocated(i_basis))then
        allocate(i_basis(n_centers_basis_T),stat=info)
        call check_allocation(info, 'i_basis                       ')
     end if
     if(.not.allocated(radial_wave))then
        allocate(radial_wave(n_max_compute_fns_missing_dens),stat=info)
        call check_allocation(info, 'radial_wave                   ')
     end if
     if(.not. allocated(wave))then
        allocate(wave(n_max_compute_missing_dens, n_max_batch_size),stat=info)
        call check_allocation(info, 'wave                          ')
     end if

     l_ylm_max = l_wave_max

     if(.not. allocated( ylm_tab))then
        allocate( ylm_tab( (l_ylm_max+1)**2,n_max_compute_missing_atoms),stat=info )
        call check_allocation(info, 'ylm_tab                       ')
     end if
     if(.not. allocated( index_lm))then
        allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info)
        call check_allocation(info, 'index_lm                      ')
      end if

     ! initialize index_lm
     i_index = 0
     do i_l = 0, l_ylm_max, 1
        do i_m = -i_l, i_l
           i_index = i_index + 1
           index_lm(i_m, i_l) = i_index
        enddo
     enddo

     ! -------- go over the grids and construct the density using density matrix -------------

     i_basis_fns_inv = 0

     i_full_points = 0
     i_full_points_2 = 0
     i_full_points_3 = 0

     do i_my_batch = 1, n_my_batches, 1

           n_compute_c = 0
           i_basis = 0

           i_point = 0

           ! loop over one batch
           do i_index = 1, batches(i_my_batch)%size, 1

              i_full_points = i_full_points + 1

              ! VB: This criterion must be exactly the inverse criterion of that used in update_density_densmat!
              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).le.0.d0) then

                 i_point = i_point+1

                 !     get current integration point coordinate
                 coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)

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
                 !
                 ! prune_basis_once must not be used for missing density - points are different!
                 ! if (.not.prune_basis_once) then
                   call prune_basis_p2 &
                   ( dist_tab_sq(1,i_point), &
                     n_compute_c, i_basis,  &
                     n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals  )
                 ! end if
              end if
           enddo ! end loop over the angular shell

           ! prune_basis_once must not be used for missing density - points are different!
           ! if (prune_basis_once) then
           !   n_compute_c = batches(i_my_batch)%batch_n_compute
           !   i_basis(1:n_compute_c) = batches(i_my_batch)%batch_i_basis
           ! end if

           if(packed_matrix_format /= PM_none )then
              call  prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
                   n_compute_c, i_basis)
           else
              call  prune_density_matrix(density_matrix, density_matrix_con, &
                   n_compute_c, i_basis)
           end if


           ! from list of n_compute active basis functions in batch, collect all atoms that are ever needed in batch.
           call collect_batch_centers_p2 &
           ( n_compute_c, i_basis, n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals, &
             n_batch_centers, batch_center &
           )

           n_points = i_point

           if (n_compute_c.gt.0) then

              ! Determine all radial functions, ylm functions and their derivatives that
              ! are best evaluated strictly locally at each individual grid point.
              i_point = 0
              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1

                 ! VB: This criterion must be exactly the inverse criterion of that used in update_density_densmat!
                 !     See above!
                 !     but now additionally make sure that we are not _on_ the nucleus of another atom,
                 !     to avoid a division by zero
                 if (max(partition_tab(i_full_points_2),&
                      hartree_partition_tab(i_full_points_2)).le.0.d0) then

                    i_point = i_point+1
                    n_compute_atoms = 0
                    n_compute_fns = 0

                   point_on_atom = .false.
                   do i_center = 1, n_centers_integrals, 1
                      ! same check as in the setting up of the partition tab: a point must not be within the 
                      ! innermost grid shell of any other atom 
                     if ( dist_tab_sq(i_center,i_point).lt.r_grid_min_sq(species_center(centers_basis_integrals(i_center)))) then
                       point_on_atom = .true.
                       exit ! exit do loop
                     end if
                   end do

                   if (point_on_atom) then
                   ! set waves (the only quantity to be computed) to zero, density not needed!

                     wave(1:n_compute_c,i_point) = 0.d0

                   else
                   ! go ahead with normal calculation

                    ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                    ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                    ! without any copying and without doing any unnecessary operations. 
                    ! The price is that the interface is no longer explicit in terms of physical 
                    ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                    call prune_radial_basis_p2 &
                    ( n_max_compute_missing_atoms, n_max_compute_fns_missing_dens, &
                       dist_tab_sq(1,i_point), dist_tab, dir_tab(1,1,i_point), &
                       n_compute_atoms, atom_index, atom_index_inv, &
                       n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                       i_atom_fns, spline_array_start, spline_array_end, &
                       n_centers_basis_integrals, centers_basis_integrals, n_compute_c, i_basis, &
                       n_batch_centers, batch_center, &
                       one_over_dist_tab, rad_index, wave_index, l_index, l_count, &
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
                         .false., n_compute_c, n_max_compute_fns_missing_dens   &
                         )

                    ! compute trigonometric functions of spherical coordinate angles
                    ! of current integration point, viewed from all atoms

                    call tab_trigonom_p0 &
                         ( n_compute_atoms, dir_tab(1,1,i_point),  &
                         trigonom_tab(1,1) &
                         )


                    ! tabulate distance and Ylm's w.r.t. other atoms            
                    call tab_wave_ylm_p0 &
                      ( n_compute_atoms, atom_index,  &
                      trigonom_tab(1,1), basis_l_max,  &
                      l_ylm_max, &
                      ylm_tab(1,1) )

                    ! tabulate total wave function value for each basis function in all cases -
                    ! but only now are we sure that we have ylm_tab ...

                    ! tabulate total wave function value for each basis function
                    call evaluate_waves_p2  &
                    ( n_compute_c, n_compute_atoms, n_compute_fns, &
                      l_ylm_max, ylm_tab, one_over_dist_tab,   &
                      radial_wave, wave(1,i_point), &
                      rad_index, wave_index, l_index, l_count, fn_atom, &
                      n_zero_compute, zero_index_point &
                    )

                    ! Not yet production-ready! Kept for completeness only!
                    !if(use_small_component )then
                    !
                    !   call small_component(n_compute_c, i_basis,  n_compute_atoms, dist_tab_sq(1,i_point), &
                    !        atom_index_inv(1), atom_index(1), gradient_basis_wave(1,1,i_point), &
                    !        temp_rho_small(i_point,i_spin), n_max_compute_dens, i_spin)
                    !end if

                    ! reset i_basis_fns_inv
                    i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

                   end if ! point_on_atom

                 end if ! end if (partition_tab.gt.0)
              enddo ! end loop over a batch
              ! - all quantities that are evaluated pointwise are now known ...


              if (n_points.gt.0) then
                 ! Now perform all operations which are done across the entire integration
                 ! shell at once.

                 call evaluate_KS_density_densmat &
                      (  n_points, wave(1,1), n_compute_c,   &
                      temp_rho(1,i_spin), n_max_compute_missing_dens, &
                      n_centers_basis_T, density_matrix_con, work(1,1) )

                 if(use_small_component )then
                    temp_rho = temp_rho + temp_rho_small
                    temp_rho_small = 0.d0
                 end if
             
              end if ! (n_points>0)

              ! calculate change in electron density
              i_point = 0
              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_3 = i_full_points_3 + 1

                 ! again, the "inverted" criterion for the density
                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).le.0) then

                    i_point = i_point + 1

                    rho(i_spin,i_full_points_3) = temp_rho(i_point, i_spin)

                 endif

              end do

           else
           ! Even if n_compute .eq. 0 for the entire current batch of grid points, we still need to
           ! make sure that the density _change_ at this point is ( zero minus previous density ).
           ! This ensures that. even for a zero KS density, a potentially non-zero initialization density
           ! is subtracted to properly account for the density change ....

             do i_index = 1, batches(i_my_batch)%size, 1

               i_full_points_2 = i_full_points_2 + 1
               i_full_points_3 = i_full_points_3 + 1

               if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).le.0.d0) then

                 rho(i_spin,i_full_points_3) = 0.d0

               endif !(partition_tab)

             enddo ! (loop over all grid points, for empty batches)

           end if  ! end if (n_compute.gt.0)
        ! end if ! end distribution over threads
     end do ! end loop over grid batches
     !--------------------------------- end go over the grids ------------------------------
  end do ! end i_spin

  !---------- finally, deallocate stuff -----------------

  if (allocated( temp_rho_small       )) deallocate( temp_rho_small       )
  if (allocated( index_lm             )) deallocate( index_lm             )
  if (allocated( ylm_tab              )) deallocate( ylm_tab              )
  if (allocated( wave                 )) deallocate( wave                 )
  if (allocated( radial_wave          )) deallocate( radial_wave          )
  if (allocated( i_basis              )) deallocate( i_basis              )
  if (allocated( work                 )) deallocate( work                 )
  if (allocated( density_matrix_con   )) deallocate( density_matrix_con   )
  if (allocated( trigonom_tab         )) deallocate( trigonom_tab         )
  if (allocated( dir_tab              )) deallocate( dir_tab              )
  if (allocated( i_r                  )) deallocate( i_r                  )
  if (allocated( dist_tab_sq          )) deallocate( dist_tab_sq          )
  if (allocated( dist_tab             )) deallocate( dist_tab             )
  if (allocated( work_complex         )) deallocate( work_complex         )
  if (allocated( density_matrix       )) deallocate( density_matrix       )

end subroutine update_missing_density_densmat
!******
