!****s* FHI-aims/main
!  NAME
!    integrate_real_kinetic_matrix -- integrates the kinetic energy matrix
!    
!  SYNOPSIS
subroutine integrate_real_kinetic_matrix &
     (rho, &
      partition_tab, basis_l_max, &
      kinetic &
      )
!  PURPOSE
!    Calculates kinetic matrix elements for use in the kinetic energy
!    evaluation. Still in test phase, needs some comments cleaning. 
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
  use species_data, only: species_name 
  use energy_density
  use aims_memory_tracking, only : aims_allocate, aims_deallocate
  use pbc_lists
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    kinetic 
!  SEE ALSO
!    FHI-aims CPC publication (see copyright notice below)
!  COPYRIGHT
!    Copyright by FHI-aims Team,
!    Theory Department, Fritz Haber Institute, Berlin (2008).
!
!    FHI-aims comes without warranty of any kind. We hope that it will be
!    useful for your work, and if you contact us, will try to help solve
!    any problems as best as we can.
!
!    FHI-aims has been distributed to you including its source code, 
!    but the terms of the license agreement apply. In particular, do 
!    not redistribute FHI-aims. Instead, please refer all interested 
!    parties to
!
!      aims-coordinators@fhi-berlin.mpg.de
!
!    instead. If you use FHI-aims, cite
!
!      Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!      Xinguo Ren, Karsten Reuter, and Matthias Scheffler, 
!      "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!      Computer Physics Communications (2008), submitted.
!
!    (or preferably the final version of this article.)
!******
!  Any other comments (not for documentation) follow below.

!  Declaration of variables



  implicit none

!  real*8, dimension(n_full_points) :: hartree_potential

  real*8, dimension(n_spin, n_full_points) :: rho

  !     These should only ever be referenced if (use_gga)
  !     import dimensions from above (if not used, all dimensions=1)

  real*8, dimension( n_full_points) :: partition_tab

  integer basis_l_max (n_species)

  !  output


  real*8 kinetic(n_hamiltonian_matrix_size, n_spin )

  !  local variables

  real*8, dimension(n_spin) :: local_potential_parts

  integer :: l_ylm_max
  integer, dimension(:,:), allocatable :: index_lm
  real*8, dimension(:,:), allocatable :: ylm_tab

  real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
  real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

  real*8 coord_current(3)
  real*8 dist_tab(n_centers_integrals, n_max_batch_size)
  real*8 dist_tab_sq(n_centers_integrals, n_max_batch_size)
  real*8 i_r(n_max_compute_atoms)
  real*8 dir_tab(3,n_centers_integrals, n_max_batch_size)
  real*8 trigonom_tab(4,n_max_compute_atoms)


  real*8,dimension(:,:,:),allocatable:: T_times_psi
  real*8,dimension(:)  ,allocatable:: radial_wave
  real*8,dimension(:)  ,allocatable:: radial_wave_deriv
  real*8,dimension(:)  ,allocatable:: kinetic_wave
  real*8,dimension(:,:)  ,allocatable:: wave



  real*8, dimension(n_spin,n_max_batch_size) :: local_rho

  !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
  !     The hope is that such a separate treatment will allow to minimize numerical noise
  !     introduced through ZORA
  real*8, dimension(:,:), allocatable :: kinetic_shell

  !     optimal accounting for matrix multiplications: only use points with nonzero components
  integer :: n_points
  integer :: n_rel_points

  !     and condensed version of hamiltonian_partition_tabs on angular grids
  real*8 :: partition(n_max_batch_size)
  real*8 :: energy_partition(n_max_batch_size)

  real*8, dimension(:,:), allocatable :: gradient_basis_wave

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
  integer :: i_basis_fns(n_basis_fns*n_centers_integrals)
  integer :: i_basis_fns_inv(n_basis_fns,n_centers)
  integer :: i_atom_fns(n_basis_fns*n_centers_integrals)

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

  ! CC: Variables used for the pointwise iteration
  !     in the case of energy density
  integer                :: i_point_2
  integer                :: i_index_2
  integer,dimension(1:2) :: i_full_points_3

  integer :: i_spin
  character*100 :: info_str

  integer :: i_my_batch

  integer :: i_radial, i_angular, info


  !  begin work

  write(info_str,'(2X,A)') "Integrating kinetic matrix: batch-based integration."
  call localorb_info(info_str,use_unit,'(A)')

  !     begin with general allocations
  if ((flag_rel.eq.REL_none.or.flag_rel==REL_atomic_zora.or.flag_rel.eq.REL_own).and.(.not.(use_gga))) then
     !       no gradients needed
     l_ylm_max = l_wave_max
  else if ((flag_rel.eq.REL_zora).or.(use_gga).or.(flag_rel==REL_KOLNING_HARMON)) then
     l_ylm_max = l_wave_max
     call aims_allocate (gradient_basis_wave, n_max_compute_ham,3, "gradient_basis_wave" )

     allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info)
     call check_allocation(info, 'dylm_dtheta_tab               ')

     allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_max_compute_atoms ) ,STAT=info)
     call check_allocation(info, 'scaled_dylm_dphi_tab          ')

  end if

  allocate( ylm_tab( (l_ylm_max+1)**2, n_max_compute_atoms ),STAT=info )
  call check_allocation(info, 'ylm_tab                       ')

  allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max), STAT=info ) 
  call check_allocation(info, 'index_lm                      ')

  call aims_allocate( kinetic_shell, n_max_compute_ham, n_max_compute_ham,      "kinetic_shell" ) 
  call aims_allocate( T_times_psi, n_max_compute_ham, n_max_batch_size, n_spin,   "T_times_psi" )

  allocate(radial_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave                   ')

  allocate(radial_wave_deriv(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'radial_wave_deriv             ')

  allocate(kinetic_wave(n_max_compute_fns_ham), STAT=info )
  call check_allocation(info, 'kinetic_wave                  ')

  call aims_allocate( wave, n_max_compute_ham, n_max_batch_size,                         "wave" )

  allocate(i_basis(n_centers_basis_I), STAT=info)
  call check_allocation(info, 'i_basis                       ')

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
        call aims_allocate(zora_vector1, n_max_compute_ham,3,n_max_batch_size,n_spin, "zora_vector1" )
     end if
     if (.not.allocated(zora_vector2)) then
        call aims_allocate(zora_vector2, n_max_compute_ham,3,n_max_batch_size,n_spin, "zora_vector2" )
     end if

  end if
  
  ! CC: Allocate Chetty & Martin kinetic energy density
  if (flag_chetty_martin_energy_density) then

    if (.not. allocated(ed_kinetic_energy_density) ) then
      allocate(ed_kinetic_energy_density(n_full_points), STAT=info)
      call check_allocation(info, 'ed_kinetic_energy_density ')
    end if
    if (.not. allocated(ed_kinetic_batch) ) then
      allocate(ed_kinetic_batch(n_hamiltonian_matrix_size, n_spin,n_full_points), STAT=info)
      call check_allocation(info, 'ed_kinetic_batch ')
    end if

    ed_kinetic_energy_density(:)  = 0.0d0
    ed_kinetic_batch(:,:,:)       = 0.0d0
    ! CC: Here, kinetic is computed through recursive integration - hence make sure that the
    ! starting value is 0
    kinetic(:,:)                  = 0.0d0
  end if
     


  !     initialize

!flag_rel =0

  kinetic = 0.d0

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
  i_full_points_3(:) = 0 


  !     perform partitioned integration, atom by atom, and point by point 
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
                 
        n_rel_points = 0
        i_point = 0


        ! loop over one batch of integration points
        do i_index = 1, batches(i_my_batch)%size, 1



           ! Increment the (global) counter for the grid, to access storage arrays
           i_full_points = i_full_points + 1

           if (partition_tab(i_full_points).gt.0.d0) then

              i_point = i_point+1





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


              n_compute_atoms = 0
              n_compute_fns = 0
              !!! i_basis_fns_inv = 0


              
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

              if ( (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON) ) then

                 ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                 call tab_gradient_ylm_p0  &
                      ( trigonom_tab(1,1), basis_l_max,   &
                      l_ylm_max, n_compute_atoms, atom_index,  &
                      ylm_tab(1,1),   &
                      dylm_dtheta_tab(1,1),   &
                      scaled_dylm_dphi_tab(1,1)  )
                               
              else
                ! tabulate distance and Ylm's w.r.t. other atoms            
                call tab_wave_ylm_p0 &
                   ( n_compute_atoms, atom_index,  &
                   trigonom_tab, basis_l_max,  &
                   l_ylm_max, ylm_tab )
              end if                  
     
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



              ! in the remaining part of the subroutine, some decisions (scalar
              !  relativity) depend on the potential; must therefore evaluate the 
              ! potential and derived quantities right here


             


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






              if ( (flag_rel.eq.REL_zora).or.(flag_rel==REL_KOLNING_HARMON)) then
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
                 
                 ! and finally, assemble the actual gradients
                 call evaluate_wave_gradient_p2  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   one_over_dist_tab, dir_tab(1,1,i_point), trigonom_tab(1,1),  &
                   l_ylm_max, ylm_tab,  &
                   dylm_dtheta_tab,  &
                   scaled_dylm_dphi_tab,  &
                   radial_wave,  &
                   radial_wave_deriv,  &
                   gradient_basis_wave,  &
                   rad_index, wave_index, l_index, l_count, fn_atom, &
                   n_zero_compute, zero_index_point &
                 )

              end if



              
              ! Now, evaluate vector of components H*phi(i,r)
              ! Local potential parts first; in the case of GGA, 
              ! the real gradient parts are added further below
              !               if ( (flag_rel/=1)) then
              ! Non-relativistic treatment - simply evaluate 
              ! H*phi(i,r) all in one
              
              ! First, obtain radial kinetic energy terms from vectorized splines
              call evaluate_radial_functions_p0  &
                   ( spline_array_start, spline_array_end,  &
                   n_compute_atoms, n_compute_fns,   &
                   dist_tab(1,i_point), i_r,  &
                   atom_index, i_basis_fns_inv,  &
                   basis_kinetic_ordered, kinetic_wave(1),  &
                   .false., n_compute_c, n_max_compute_fns_ham )



              do i_spin = 1, n_spin, 1
                 call evaluate_T_psi  &
                 ( n_compute_c, n_compute_atoms, n_compute_fns, &
                   l_ylm_max, ylm_tab, one_over_dist_tab,  &
                   T_times_psi(1, i_point, i_spin),  &
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



           end if  ! end if (hamiltonian_partition_tab.gt.0)
        enddo ! end loop over a batch
        


        ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
        ! work separately for each spin channel
        do i_spin = 1, n_spin, 1

           
           ! CC: Do "on-the-fly" integration over batch if no energy density is requested ...
           if (.not. flag_chetty_martin_energy_density) then

             ! Compute for whole batch
           ! add full non-relativistic contributions and (for relativistic points)
           ! all contributions from the potential to the Hamiltonian matrix elements
           call evaluate_hamiltonian_shell_p1  &
                ( n_points, partition(1), n_compute_c, &
                T_times_psi(1,1,i_spin),  &
                n_max_compute_ham, wave(1,1),  &
                kinetic_shell )

           ! For all relativistic points, add kinetic energy contributions
           if (n_rel_points.gt.0) then
              

              call add_zora_matrix_p1 ( &
                   zora_vector1(1,1,1,i_spin),   &
                   zora_vector2(1,1,1,i_spin),   &
                   n_max_compute_ham, n_rel_points, &
                   n_compute_c, &
                   kinetic_shell )

             end if

             call update_full_matrix_p0(  &
                   n_compute_c, n_compute_c, i_basis(1), kinetic_shell,    &
                   kinetic(1,i_spin) )
             
           ! CC: ... work point by point if the computation of the kinetic energy density is requested!
           else
             ! Loop over batch
             i_point_2          = 0
             do i_index=1,batches(i_my_batch)%size, 1

               ! Increment the (global) counter for the grid, to access storage arrays
               i_full_points_3(i_spin) = i_full_points_3(i_spin) + 1
               
               if (partition_tab(i_full_points_3(i_spin)).gt.0.d0) then
                 i_point_2 = i_point_2+1
                 
                 ! Compute for single point in batch
                 kinetic_shell(:,:) = 0.0d0
                 call evaluate_hamiltonian_shell_p1  &
                      ( 1, 1.0d0, n_compute_c, &
                      T_times_psi(1,i_point_2,i_spin),  &
                      n_max_compute_ham, wave(1,i_point_2),  &
                      kinetic_shell )
 
                 ! FIXME For all relativistic points, add kinetic energy contributions
                 if (n_rel_points.gt.0) then
                   write(use_unit,*) " *** WARNINING *** Zora not (yet) implemented in the calculation of kinetic energy densities"
                   call aims_stop()
                 end if

                 ed_kinetic_batch(:,i_spin,i_full_points_3(i_spin)) = 0.0d0
                 call update_full_matrix_p0(  &
                     n_compute_c, n_compute_c, i_basis(1), kinetic_shell,    &
                     ed_kinetic_batch(:,i_spin,i_full_points_3(i_spin)) )

                 !Integrate to get complete kinetic energy
                 do i_index_2=1,n_hamiltonian_matrix_size
                     kinetic(i_index_2,i_spin) = kinetic(i_index_2,i_spin) &
                        + ed_kinetic_batch(i_index_2,i_spin, &
                                           i_full_points_3(i_spin)) &
                        * partition_tab(i_full_points_3(i_spin))
                 end do

               end if
             end do

           end if



        enddo

        ! Hamiltonian is now complete.
        !
        ! Since we already have the pieces, add terms of XC energy here. 
        ! Notice that these terms are not added for ANY shell
        ! where n_compute happens to be zero. This should be correct because all wave functions
        ! are zero here anyway, i.e. also the density.
 
     else

       i_full_points = i_full_points + batches(i_my_batch)%size

     end if ! end if (n_compute.gt.0) then


  ! end if  ! end mpi work distribution
end do ! end loop over bathces


  !     synchronise the kinetic
  if(.not. use_local_index) call sync_integrate_hamiltonian( kinetic )

  ! Allocatable arrays that are tracked
  if (allocated(gradient_basis_wave)) call aims_deallocate( gradient_basis_wave, "gradient_basis_wave" )
  if (allocated(kinetic_shell))       call aims_deallocate( kinetic_shell,             "kinetic_shell" ) 
  if (allocated(T_times_psi))         call aims_deallocate( T_times_psi,                 "T_times_psi" )
  if (allocated(wave))                call aims_deallocate( wave,                               "wave" )
  if (allocated(zora_vector1))        call aims_deallocate( zora_vector1,               "zora_vector1" )
  if (allocated(zora_vector2))        call aims_deallocate( zora_vector2,               "zora_vector2" )
  
  if (allocated(scaled_dylm_dphi_tab)) then
     deallocate (scaled_dylm_dphi_tab)
  end if
  if (allocated(dylm_dtheta_tab)) then
     deallocate (dylm_dtheta_tab)
  end if

  if (allocated(ylm_tab)) then
     deallocate(ylm_tab)
  end if
  if (allocated(index_lm)) then
     deallocate(index_lm)
  end if
  if (allocated(dist_tab_full)) then
     deallocate(dist_tab_full)
  end if
  if (allocated(dir_tab_full_norm)) then
     deallocate(dir_tab_full_norm)
  end if
  if (allocated(i_r_full)) then
     deallocate(i_r_full)
  end if
  if(allocated(radial_wave))then
     deallocate(radial_wave)
  end if
  if(allocated(radial_wave_deriv))then
     deallocate(radial_wave_deriv)
  end if  
  if(allocated(kinetic_wave))then
     deallocate(kinetic_wave)
  end if
 if( allocated(i_basis))then
    deallocate(i_basis)
 end if

end subroutine integrate_real_kinetic_matrix

!----------------------------------------------------------------------
