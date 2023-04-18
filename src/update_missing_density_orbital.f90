!****s* FHI-aims/update_missing_density_orbital
!  NAME
!   update_missing_density_orbital
!  SYNOPSIS

subroutine update_missing_density_orbital &
     ( KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, partition_tab,  &
     hartree_partition_tab, rho, basis_l_max &
     )

  !  PURPOSE
  !  Subroutine update_missing density obtains the KS density from the eigenvectors
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
  use KH_core_states
  implicit none

  !  ARGUMENTS


  !  input 

  real*8,     dimension(n_basis, n_states, n_spin,n_k_points):: KS_eigenvector
  complex*16, dimension(n_basis, n_states, n_spin,n_k_points):: KS_eigenvector_complex

  real*8, dimension(n_states, n_spin, n_k_points) :: occ_numbers
  real*8, dimension(n_states, n_spin, n_k_points) :: KS_eigenvalue      

  real*8, dimension(n_full_points) ::  partition_tab
  real*8, dimension(n_full_points) ::  hartree_partition_tab
  !     rho is only input; what we store is the density residual (i.e. the change in the density)
  real*8, dimension(n_spin, n_full_points) :: rho

  integer basis_l_max (n_species)

  !  INPUTS
  !   o KS_eigenvector -- Kohn-Sham eigenvectors (real format)
  !   o KS_eigenvector_complex -- Kohn-Sham eigenvectors (complex format)
  !   o occ_numbers -- occupations of eigenstates
  !   o KS_eigenvalue -- Kohn-Sham eigenvalues
  !   o partition_tab -- values of partition function
  !   o hartree_partition_tab -- values of partition function used in multipole expansion of charges in Hartree potential.
  !   o basis_l_max -- maximum l of basis functions
  !    
  !  OUTPUT
  !   o rho -- electron density - rho is here also output, we only update it where rho was not computed before.
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
  !

  !  local variables

  real*8 coord_current(3)
  real*8 dist_tab(n_centers_basis_integrals  , n_max_batch_size)
  real*8 dist_tab_sq(n_centers_basis_integrals, n_max_batch_size)
  real*8 i_r(n_centers_basis_integrals)
  real*8 dir_tab(3, n_centers_basis_integrals, n_max_batch_size)
  real*8 trigonom_tab(4, n_centers_basis_integrals)

  real*8, dimension(n_species)        :: r_grid_min_sq
  real*8, dimension(:),   allocatable :: radial_wave
  real*8, dimension(:,:), allocatable :: wave

  !     pruning of atoms, radial functions, and basis functions

  integer :: n_compute_c
  integer,dimension(:),allocatable :: i_basis

  integer :: n_compute_fns
  integer :: i_basis_fns(n_basis_fns*n_centers_basis_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_basis_integrals)

  integer :: n_compute_atoms
  integer :: atom_index(n_centers_basis_integrals)
  integer :: atom_index_inv(n_centers)

  integer :: spline_array_start(n_centers_basis_integrals)
  integer :: spline_array_end(n_centers_basis_integrals)

  ! VB - renewed index infrastructure starts here

  real*8 one_over_dist_tab(n_max_compute_missing_atoms)

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
  integer, dimension(n_spin) :: max_occ_number
  real*8 :: occ_numbers_sqrt

  integer :: l_ylm_max
  integer :: n_compute_force_atoms
  integer :: n_local_compute
  integer :: n_points

  real*8,     dimension(:,:,:),allocatable :: KS_vec_times_occ_sqrt
  complex*16, dimension(:,:,:),allocatable :: KS_vec_times_occ_sqrt_complex

  real*8,     dimension(:,:,:),allocatable :: KS_ev_compute
  complex*16, dimension(:,:,:),allocatable :: KS_ev_compute_complex

  real*8,    dimension(:,:,:),allocatable :: KS_orbital
  complex*16,dimension(:,:,:),allocatable :: KS_orbital_complex

  real*8 temp_rho(n_max_batch_size,n_spin)
  real*8,dimension(:,:),allocatable :: temp_rho_small

  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  !     counters

  integer i_atom
  integer i_l
  integer i_m
  integer i_coord
  integer :: i_state
  integer :: i_point
  integer :: i_bas

  integer :: i_my_batch

  integer :: i_full_points
  integer :: i_full_points_2
  integer :: i_full_points_3

  logical :: point_on_atom
  integer :: i_center

  ! i_spin for future use in a spin-polarized version
  integer :: i_spin = 1

  integer ::  i_index, i_k_point, i_k_point_group

  !      integer, dimension(:,:), allocatable :: n_points_mpi
  integer :: i_atom_2
  integer :: info
  integer:: n_k_group, i_k_point_g

  !     begin work


  call localorb_info( &
    "Evaluating KS electron density on previously omitted grid points.", use_unit,'(2X,A)',OL_norm )

  if(use_small_component)then
     allocate(temp_rho_small(n_max_batch_size,n_spin),stat=info)
     call check_allocation(info, 'temp_rho_small                ') 
  end if

  r_grid_min_sq(:) = r_grid_min(:)*r_grid_min(:)

  if(real_eigenvectors)then

     if(.not. allocated( KS_vec_times_occ_sqrt))then
        allocate( KS_vec_times_occ_sqrt(n_centers_basis_T,(n_states*n_k_points_group), n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_vec_times_occ_sqrt'
           stop
        end if
     end if

     if(.not. allocated( KS_ev_compute))then
        allocate( KS_ev_compute(n_states*n_k_points_group, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute'
           stop
        end if
     end if

     if(.not. allocated( KS_orbital))then
        allocate( KS_orbital(n_states*n_k_points_group,n_max_batch_size,n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_orbital'
           stop
        end if
     end if

  else

     if(.not. allocated( KS_vec_times_occ_sqrt_complex))then
        allocate( KS_vec_times_occ_sqrt_complex(n_centers_basis_T, (n_states*n_k_points_group), n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation:  KS_vec_times_occ_sqrt_complex'
           stop
        end if
     end if

     if(.not. allocated( KS_ev_compute_complex))then
        allocate( KS_ev_compute_complex(n_states*n_k_points_group, n_centers_basis_T, n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_ev_compute_complex'
           stop
        end if
     end if

     if(.not. allocated( KS_orbital_complex))then
        allocate( KS_orbital_complex(n_states*n_k_points_group,n_max_batch_size,n_spin),stat=info)
        if(info/=0)then
           write(use_unit,*)'Error in allocation: KS_orbital_complex'
           stop
        end if
     end if
  end if

  if(.not. allocated(i_basis))then
     allocate(i_basis(n_centers_basis_T),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: i_basis'
        stop
     end if
  end if

  if(.not.allocated(radial_wave))then
     allocate(radial_wave(n_max_compute_fns_missing_dens),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: radial_wave'
        stop
     end if
  end if

  if(.not. allocated(wave))then
     allocate(wave(n_max_compute_missing_dens, n_max_batch_size),stat=info)
     if(info/=0)then
        write(use_unit,*)'Error in allocation: wave'
        stop
     end if
  end if

  l_ylm_max = l_wave_max

  if(.not. allocated( ylm_tab))then
     allocate( ylm_tab( (l_ylm_max+1)**2,n_centers_basis_integrals),stat=info )
     if(info/=0)then
        write(use_unit,*)'Error in allocation: ylm_tab'
        stop
     end if
  end if


  if(.not. allocated( index_lm))then
     allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info) 
     if(info/=0)then
        write(use_unit,*)'Error in allocation: index_lm'
        stop
     end if
  end if

  call check_occs('update_missing_density_orbital', occ_numbers, .true.)

  !     initialize index_lm
  i_index = 0
  do i_l = 0, l_ylm_max, 1
     do i_m = -i_l, i_l
        i_index = i_index + 1
        index_lm(i_m, i_l) = i_index
     enddo
  enddo

  ! find the maximal occupation number
  ! VB: In reality, any periodic calculations are done using the density
  !     matrix, and do not ever get here. Therefore, the use of n_k_group
  !     is somewhat academic here.
  n_k_group = ceiling(real(n_k_points)/real(n_k_points_group))

  do i_k_point_group = 1,n_k_group

     if(real_eigenvectors)then

        do i_spin = 1, n_spin, 1

           max_occ_number(i_spin) = 0
           i_index = 0


           do i_k_point_g = 1, n_k_points_group

              i_k_point = (i_k_point_group-1) * n_k_points_group + i_k_point_g

              if(i_k_point <= n_k_points)then
                 do i_state = 1, n_states, 1
                    if (occ_numbers(i_state,i_spin, i_k_point).gt.0.d0) then

                       i_index = i_index + 1

                       occ_numbers_sqrt =  sqrt(occ_numbers(i_state,i_spin,i_k_point))


                       do i_bas = 1,  n_centers_basis_T, 1



                          KS_vec_times_occ_sqrt(i_bas,i_index,i_spin) =  &
                               KS_eigenvector(Cbasis_to_basis(i_bas),i_state,i_spin, i_k_point) * &
                               occ_numbers_sqrt * &
                               dble(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),i_k_point))


                       end do
                    end if
                 end do
              end if
           enddo
           max_occ_number(i_spin) = i_index        
        enddo

     else

        do i_spin = 1, n_spin, 1


           max_occ_number(i_spin) = 0
           i_index = 0

           do i_k_point_g = 1, n_k_points_group

              i_k_point = (i_k_point_group-1) * n_k_points_group + i_k_point_g

              if(i_k_point <= n_k_points)then

                 do i_state = 1, n_states, 1
                    if (occ_numbers(i_state,i_spin, i_k_point).gt.0.d0) then

                       i_index = i_index + 1

                       occ_numbers_sqrt =  sqrt(occ_numbers(i_state,i_spin,i_k_point))


                       do i_bas = 1,  n_centers_basis_T, 1



                          KS_vec_times_occ_sqrt_complex(i_bas,i_index,i_spin) =  &
                               KS_eigenvector_complex(Cbasis_to_basis(i_bas),i_state,i_spin, i_k_point) * &
                               occ_numbers_sqrt * &
                               dconjg(k_phase(center_to_cell(Cbasis_to_center(  i_bas  )),i_k_point))


                       end do
                    end if
                 end do
              end if
           end do
           max_occ_number(i_spin) = i_index

        enddo
     end if

     !     loop over integration grid

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

              ! VB: This criterion must be exactly the inverse criterion of that used in update_density_and_forces_orbital!
              if (max(partition_tab(i_full_points),&
                   hartree_partition_tab(i_full_points)).le.0.d0) then

                 i_point = i_point+1

                 !     get current integration point coordinate
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
                      n_centers_basis_integrals, centers_basis_integrals )

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

                 ! VB: This criterion must be exactly the inverse criterion of that used in update_density_and_forces_orbital!
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
                         dist_tab_sq(1,i_point), dist_tab(1,i_point), dir_tab(1,1,i_point), &
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
                         dist_tab(1,i_point),  &
                         i_r &
                         )

                    ! Determine all needed radial functions from efficient splines

                    ! Now evaluate radial functions u(r) from the previously stored compressed 
                    ! spline arrays  
                    call evaluate_radial_functions_p0 &
                         (   spline_array_start, spline_array_end, &
                         n_compute_atoms, n_compute_fns,  &
                         dist_tab(1,i_point), i_r(1), &
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


                    ! Experimental only, kept here for completeness. Not yet production ready.
                    ! if(use_small_component )then
                    ! 
                    !   do i_spin = 1, n_spin, 1
                    !      call small_component(n_compute_c, i_basis,  n_compute_atoms, dist_tab_sq(1,i_point), &
                    !           atom_index_inv(1), atom_index(1), gradient_basis_wave(1,1,i_point), &
                    !           temp_rho_small(i_point,i_spin), n_max_compute_dens, i_spin)
                    !   end do
                    ! end if

                    ! reset i_basis_fns_inv
                    i_basis_fns_inv(:,atom_index(1:n_compute_atoms)) = 0

                   end if ! point_on_atom

                 end if ! end if (partition_tab.gt.0)
              enddo ! end loop over one batch of the grid
              ! - all quantities that are evaluated pointwise are now known ...

              if (n_points.gt.0) then
                 ! Now perform all operations which are done across the entire integration
                 ! shell at once.

                 do i_spin = 1, n_spin, 1

                    ! VB: Zero max_occ_number can happen if one spin channel is constrained
                    !     to be completely empty. 
                    if (max_occ_number(i_spin).gt.0) then

                       ! New density is always evaluated
                       if(real_eigenvectors)then

                          call evaluate_KS_density_p0  &
                               (  n_points, wave(1,1), n_compute_c,   &
                               i_basis, KS_vec_times_occ_sqrt(1,1,i_spin),   &
                               KS_ev_compute(1,1,i_spin),  &
                               max_occ_number(i_spin),  &
                               KS_orbital(1,1,i_spin),   &
                               temp_rho(1,i_spin), n_max_compute_missing_dens, &
                               n_centers_basis_T )


                          if(use_small_component )then
                             temp_rho = temp_rho + temp_rho_small
                             temp_rho_small = 0.d0
                          end if

                       else

                          call evaluate_KS_density_complex_p0  &
                               (  n_points, wave(1,1), n_compute_c,   &
                               i_basis, KS_vec_times_occ_sqrt_complex(1,1,i_spin),   &
                               KS_ev_compute_complex(1,1,i_spin),  &
                               max_occ_number(i_spin),  &
                               KS_orbital_complex(1,1,i_spin),   &
                               temp_rho(1,i_spin), n_max_compute_missing_dens, &
                               n_centers_basis_T )

                       end if

                    else
                       ! case of max_occ_number = 0
                       ! use default value zero for everything

                       temp_rho(1:n_points,i_spin) = 0.d0
                    
                    end if
                 end do
              end if ! (n_points>0)

              ! store the electron density
              i_point = 0
              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_3 = i_full_points_3 + 1

                 ! again, the "inverted" criterion for the density
                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).le.0) then

                    i_point = i_point + 1

                    if (i_k_point_group.eq.1) then
                      rho(1:n_spin,i_full_points_3) = temp_rho(i_point,1:n_spin)
                    else
                      rho(1:n_spin,i_full_points_3) = temp_rho(i_point,1:n_spin) + rho(1:n_spin,i_full_points_3) 
                    end if

                 endif ! partition_tab
              end do

           else
              ! Even if n_compute .eq. 0 for the entire current batch of grid points, we still need to
              ! fill in a zero density.

              do i_index = 1, batches(i_my_batch)%size, 1

                 i_full_points_2 = i_full_points_2 + 1
                 i_full_points_3 = i_full_points_3 + 1

                 if (max(partition_tab(i_full_points_3),&
                      hartree_partition_tab(i_full_points_3)).le.0.d0) then

                      rho(1:n_spin,i_full_points_3) = 0.d0

                 endif ! partition_tab

              enddo
           end if  ! end if (n_compute.gt.0)
        ! end if ! end distribution over threads
     end do ! end loop over batches
  end do ! end loop over groups of k-points

  ! finally, deallocate stuff.

  if (allocated( temp_rho_small       )) deallocate( temp_rho_small       )

  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if(allocated(radial_wave))then
     deallocate(radial_wave)
  end if
  if(allocated(wave))then
     deallocate(wave)
  end if
  if(allocated(i_basis))then
     deallocate(i_basis)
  end if
  if( allocated(KS_vec_times_occ_sqrt))then
     deallocate(KS_vec_times_occ_sqrt)
  end if
  if(allocated(KS_vec_times_occ_sqrt_complex))then
     deallocate(KS_vec_times_occ_sqrt_complex)
  end if
  if(allocated( KS_ev_compute))then
     deallocate( KS_ev_compute)
  end if
  if(allocated(KS_ev_compute_complex))then
     deallocate(KS_ev_compute_complex)
  end if
  if(allocated(KS_ev_compute_complex))then
     deallocate(KS_ev_compute_complex)
  end if
  if(allocated( KS_orbital))then
     deallocate( KS_orbital)
  end if
  if(allocated( KS_orbital_complex))then
     deallocate( KS_orbital_complex)
  end if

end subroutine update_missing_density_orbital
!******
