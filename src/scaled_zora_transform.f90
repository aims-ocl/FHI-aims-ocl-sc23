!****h* FHI-aims/scaled_zora_transform
!  NAME
!    scaled_zora_transform
!  SYNOPSIS

module scaled_zora_transform 

!  USES

  use dimensions
  use scalapack_wrapper
  use pseudodata
  implicit none

!  PURPOSE
!
!     The module includes routines for the scaled ZORA postprocessing. If the scalapack is in use
!     the separate subroutines in the module scalapack_wrapper are called. This is because the then
!     the eigenvectors are in the scalapack format and separate routine is needed.
!
!     In order to use the scalad zora, the self consistant calculations have to be done with ZORA
!     aproximation on.
!
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
!  SOURCE



  real*8, dimension(:,:), allocatable,private  :: SZ_hamiltonian

!******

contains
  !-----------------------------------------------------------
!****s* scaled_zora_transform/allocate_scaled_zora_transform
!  NAME
!    allocate_scaled_zora_transform
!  SYNOPSIS

  subroutine allocate_scaled_zora_transform

!  PURPOSE
!    Allocates the works space for scaled scaled zora transformation.
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

   implicit none
    integer:: info

    if (.not.allocated(SZ_hamiltonian)) then
       allocate( SZ_hamiltonian(n_hamiltonian_matrix_size, n_spin),stat=info)
       call check_allocation(info, 'SZ_hamiltonian                ')
    end if


  end subroutine allocate_scaled_zora_transform
!******
!---------------------------------------------------------------------------------
!****s* scaled_zora_transform/deallocate_scaled_zora_transform
!  NAME
!    deallocate_scaled_zora_transform
!  SYNOPSIS

  subroutine deallocate_scaled_zora_transform

!  PURPOSE
!    Deallocates the works space of the module scaled_zora_transform.
!  INPUTS
!    none 
!  OUTPUT
!    none 
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

    implicit none

    if (allocated(SZ_hamiltonian)) then
       deallocate( SZ_hamiltonian)
    end if


  end subroutine deallocate_scaled_zora_transform
!******
!---------------------------------------------------------------------------------
!****s* scaled_zora_transform/integrate_scaled_zora_transf_p2
!  NAME
!    integrate_scaled_zora_transf_p2
!  SYNOPSIS

  subroutine integrate_scaled_zora_transf_p2( rho, rho_gradient, kinetic_density, &
       hartree_potential, partition_tab, basis_l_max )

!  PURPOSE
!
!    The routine integrates the scaled zora integrals and saves the results to SZ_hamiltonian.
!    After this the actual Scaled zora transformation
!    can be do to the eigenvalues by calling the subroutine evaluate_scaled_zora_transf_p1.
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

    real*8, dimension(n_spin, n_full_points)    :: rho
    real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
    real*8, dimension(n_spin, n_full_points)    :: kinetic_density
    real*8, dimension(n_full_points)            :: hartree_potential
    real*8, dimension( n_full_points)           :: partition_tab
    integer, dimension(n_species)               :: basis_l_max

!  INPUTS
!
!   o rho  -- electron density
!   o rho_gradient  -- electron density gradient, 
!                     These should only ever be referenced if (use_gga)
!                     import dimensions from above (if not used, all dimensions=1)
!   o kinetic_density -- kinetic-energy density. Only used if use_meta_gga.
!   o hartree_potential  -- Hartree potential
!   o partition_tab  -- Partition function
!   o basis_l_max    -- The maximum value of l of the basis
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE




    real*8, dimension(n_spin) :: local_potential_parts

    integer :: l_ylm_max
    integer, dimension(:,:), allocatable :: index_lm
    real*8, dimension(:,:), allocatable :: ylm_tab

    real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
    real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

    real*8 coord_current(3)
    real*8 dist_tab(n_centers_integrals)
    real*8 dist_tab_sq(n_centers_integrals, n_max_compute_ang)
    real*8 i_r(n_centers_integrals)
    real*8 dir_tab(3,n_centers_integrals, n_max_compute_ang)
    real*8 trigonom_tab(4,n_centers_integrals)

    !paula


    real*8,dimension(:,:,:),allocatable:: H_times_psi
    real*8,dimension(:)  ,allocatable:: radial_wave
    real*8,dimension(:)  ,allocatable:: radial_wave_deriv
    real*8,dimension(:)  ,allocatable:: kinetic_wave
    real*8,dimension(:,:)  ,allocatable:: wave


    !     for XC
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv
    real*8, dimension(n_spin) :: xc_tau_deriv

    ! Unused?
    ! real*8, dimension(n_spin) :: local_rho
    ! real*8, dimension(3,n_spin) :: local_rho_gradient


    !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
    !     The hope is that such a separate treatment will allow to minimize numerical noise
    !     introduced through ZORA
    real*8, dimension(:,:), allocatable :: hamiltonian_shell

    !     optimal accounting for matrix multiplications: only use points with nonzero components
    integer :: n_points
    integer :: n_rel_points

    !     and condensed version of hamiltonian_partition_tabs on angular grids
    real*8 :: partition(n_max_compute_ang)


    real*8, dimension(:,:), allocatable :: gradient_basis_wave

    !     Following is all that is needed for the handling of ZORA scalar relativity

    real*8, dimension(n_spin) :: zora_operator
    logical :: t_zora(2)
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

    !     for splitting of angular shells into "octants"

    integer division_low
    integer division_high

    !  counters

    integer i_basis_1
    integer i_basis_2
    integer i_grid
    integer i_index, i_l, i_m
    integer i_coord
    integer i_division

    integer i_species

    integer i_point, info

    integer :: i_my_batch

    integer :: i_full_points_C
    integer :: i_full_points_2C

    integer :: i_full_points_A
    !  integer :: i_full_points_2A

    integer :: i_spin, i
    character*100 :: info_str

    real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
    real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore

    !test
    !  integer i_compute_1, i_compute_2, i_state, i_k_point, info

    !test end
    !  integer,dimension(:,:),allocatable :: index_b

    !  begin work


    t_zora = .false.

    write(info_str,'(2X,A,A)') & 
         "Computing correction matrix for scaled ZORA eigenvalues."
    call localorb_info(info_str,use_unit,'(A)')

    !     begin with general allocations
    if ((flag_rel.eq.0).and.(.not.(use_gga))) then
       !       no gradients needed
       l_ylm_max = l_wave_max
    else if ((flag_rel.eq.1).or.(use_gga)) then
       l_ylm_max = l_wave_max
       allocate (gradient_basis_wave(n_max_compute_ham,3),stat=info)
       call check_allocation(info, 'gradient_basis_wave           ') 
       allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info )
       call check_allocation(info, 'dylm_dtheta_tab               ') 
       allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_centers_integrals) ,stat=info)
       call check_allocation(info, 'scaled_dylm_dphi_tab          ')  
   end if

    allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info )
    call check_allocation(info, 'ylm_tab                       ') 
    allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info ) 
    call check_allocation(info, 'index_lm                      ') 
    allocate ( hamiltonian_shell(n_max_compute_ham, n_max_compute_ham),stat=info ) 
    call check_allocation(info, 'hamiltonian_shell             ')

    allocate(H_times_psi(n_max_compute_ham, n_max_compute_ang, n_spin),stat=info)
    call check_allocation(info, 'H_times_psi                   ') 
    allocate(radial_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave                   ') 
    allocate(radial_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave_deriv             ') 
    allocate(kinetic_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave                  ') 
    allocate(wave(n_max_compute_ham, n_max_compute_ang),stat=info)
    call check_allocation(info, 'wave                          ') 
    allocate(i_basis(n_centers_basis_I),stat=info)
    call check_allocation(info, 'i_basis                       ') 

    if (flag_rel.eq.1) then
       ! allocate all arrays relevant for ZORA

       if (.not.allocated(dist_tab_full)) then
          allocate(dist_tab_full(n_centers_integrals),stat=info)
          call check_allocation(info, 'dist_tab_full                 ')
       end if
       if (.not.allocated(dir_tab_full_norm)) then
          allocate(dir_tab_full_norm(3,n_centers_integrals),stat=info)
          call check_allocation(info, 'dir_tab_full_norm             ') 
       end if
       if (.not.allocated(i_r_full)) then
          allocate(i_r_full(n_centers_integrals),stat=info)
          call check_allocation(info, 'i_r_full                      ') 
       end if

       if (.not.allocated(zora_vector1)) then
          allocate(zora_vector1(n_max_compute_ham,3,n_max_compute_ang,n_spin),stat=info)
          call check_allocation(info, 'zora_vector1                  ') 
       end if
       if (.not.allocated(zora_vector2)) then
          allocate(zora_vector2(n_max_compute_ham,3,n_max_compute_ang,n_spin),stat=info)
          call check_allocation(info, 'zora_vector2                  ') 
       end if

    end if

    SZ_hamiltonian = 0.d0

    if(use_embedding_pp.and.use_nonlinear_core) then
        allocate(rho_inc_partialcore(n_spin,n_full_points))
        allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
        do i_spin = 1,n_spin
           rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
           rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
           if(use_density_gradient) then
              rho_gradient_inc_partialcore(:,i_spin,:) = &
                 rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
           endif
        enddo
    endif


    !     initialize index_lm
    i_index = 0
    do i_l = 0, l_wave_max, 1
       do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
       enddo
    enddo



    i_full_points_C = 0
    i_full_points_2C = 0
    i_full_points_A = 0

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


             i_full_points_2C = i_full_points_2C + 1


             if (partition_tab(i_full_points_2C).gt.0.d0) then

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
                   call prune_basis_p0 &
                        ( dist_tab_sq(1,i_point), &
                        n_compute_a, n_compute_c, i_basis,  &
                        n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
                end if
             end if
          enddo  ! end loop over a batch of points

          if (prune_basis_once) then
             n_compute_a = batches(i_my_batch)%batch_n_compute
             n_compute_c = n_compute_a
             i_basis(1:n_compute_a) = batches(i_my_batch)%batch_i_basis
          end if

          n_points = i_point


          ! Perform actual integration if more than 0 basis functions
          ! are actually relevant on the present angular shell ...
          if (n_compute_c.gt.0) then

             n_rel_points = 0
             i_point = 0

             do i_index = 1, batches(i_my_batch)%size, 1

                ! Increment the (global) counter for the grid, to access storage arrays
                i_full_points_C = i_full_points_C + 1
                i_full_points_A = i_full_points_A + 1

                if (partition_tab(i_full_points_C).gt.0.d0) then

                   i_point = i_point+1


                   call tab_global_geometry_p0 &
                        ( dist_tab_sq(1,i_point), &
                        dir_tab(1,1,i_point), &
                        dist_tab_full, &
                        i_r_full, &
                        dir_tab_full_norm, &
                        n_centers_integrals,  centers_basis_integrals)

                   ! for all integrations
                   partition(i_point) = partition_tab(i_full_points_C)
                   ! energy_partition(i_point) = partition_tab(i_full_points_A)

                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                   ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                   ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                   ! without any copying and without doing any unnecessary operations. 
                   ! The price is that the interface is no longer explicit in terms of physical 
                   ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.


                   call prune_radial_basis_p0 &
                        ( dist_tab_sq(1,i_point), &
                        dist_tab, &
                        dir_tab(1,1,i_point), &
                        n_compute_atoms, atom_index, atom_index_inv, &
                        n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                        i_atom_fns, spline_array_start, spline_array_end, &
                        n_centers_integrals, centers_basis_integrals)

                   ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                   ! for all atoms which are actually relevant
                   call tab_local_geometry_p0 &
                        ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                        dir_tab(1,1,i_point), dist_tab,  &
                        i_r &
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


                   ! Now evaluate radial functions
                   ! from the previously stored compressed spline arrays  
                   call evaluate_radial_functions_p0  &
                        (   spline_array_start, spline_array_end,  &
                        n_compute_atoms, n_compute_fns,   &
                        dist_tab, i_r,  &
                        atom_index, i_basis_fns_inv,  &
                        basis_wave_ordered, radial_wave(1),  &
                        .false. , n_compute_c, n_max_compute_fns_ham )

                   ! tabulate total wave function value for each basis function
                   call evaluate_waves_p0  &
                        ( l_ylm_max,   &
                        ylm_tab(1,1), dist_tab,   &
                        index_lm, n_compute_c,   &
                        i_basis, radial_wave(1),   &
                        wave(1,i_point), n_compute_atoms,   &
                        atom_index_inv, n_compute_fns,  &
                        i_basis_fns_inv,  n_max_compute_fns_ham )


                   ! in the remaining part of the subroutine, some decisions (scalar
                   !  relativity) depend on the potential; must therefore evaluate the 
                   ! potential and derived quantities right here

                   ! Local exchange-correlation parts of the potential are evaluated
                   ! right here, to avoid having to store them separately elsewhere.
                   ! For large systems, savings are significant
                   ! This one is required for vdW-DF
                   coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)!SAG
                   if(n_periodic > 0)then
                      call map_to_center_cell(coord_current(1:3) )
                   end if

                   if(use_embedding_pp.and.use_nonlinear_core) then
! here we might add partial core densities from pseudocores
! only for the evaluation of the xc-energy

                        call evaluate_xc  &
                            ( rho_inc_partialcore(1,i_full_points_A),   &
                            rho_gradient_inc_partialcore(1,1,i_full_points_A),  &
                            kinetic_density(1,i_full_points_A), &
                            en_density_xc, &
                            en_density_x, en_density_c,  &
                            local_xc_derivs(1),  &
                            xc_gradient_deriv(1,1),  &
                            xc_tau_deriv(1), &
                            .false., &
                            coord_current(:) &
                            )

                   else

                        call evaluate_xc  &
                            ( rho(1,i_full_points_A),   &
                            rho_gradient(1,1,i_full_points_A),  &
                            kinetic_density(1,i_full_points_A), &
                            en_density_xc, &
                            en_density_x, en_density_c,  &
                            local_xc_derivs(1),  &
                            xc_gradient_deriv(1,1),  &
                            xc_tau_deriv(1), &
                            .false., &
                            coord_current(:) &
                            )

                   endif

                   do i_spin = 1, n_spin, 1
                      local_potential_parts(i_spin) =   &
                           hartree_potential(i_full_points_A) +   &
                           local_xc_derivs(i_spin)

                      if (use_gga) then
                         sum_of_local_gradients(1:3,i_spin) =  0.d0
                      else
                         sum_of_local_gradients = 0.d0
                      end if


                   enddo

                   ! Check whether relativistic corrections are needed at the present point.
                   ! The check is based entirely on the local parts of the potential - i.e. 
                   ! in a GGA, the terms due to d(rho*exc)/d(|grad(rho|^2) is not evaluated.
                   ! Hopefully this approximation to the full ZORA energy is small.
                   if (flag_rel.eq.1) then

                      ! if we need ZORA, must get the _full_ local geometry in order to 
                      ! create the superposition of atomic potentials which is used to estimate
                      ! the potential gradient for ZORA


                      call evaluate_pot_superpos_p0  &
                           (   &
                           i_r_full,   &
                           zora_potential_parts(1),  &
                           n_centers_integrals, centers_basis_integrals ) 

                      do i_spin = 1, n_spin, 1

                         zora_operator(i_spin) =  &
                            2* light_speed_sq /  &
                              ( 2 * light_speed_sq -  &
                              zora_potential_parts(i_spin) )**2
                      enddo

                   end if

                   if ((use_gga) .or. (flag_rel.eq.1)) then
                      ! we require the gradient of each basis function

                      ! tabulate radial derivatives of those radial functions 
                      ! which are actually non-zero at current point, using vectorized splines
                      call evaluate_radial_functions_p0  &
                           ( spline_array_start, spline_array_end,  &
                           n_compute_atoms, n_compute_fns,   &
                           dist_tab, i_r,  &
                           atom_index, i_basis_fns_inv,  &
                           basis_deriv_ordered,   &
                           radial_wave_deriv(1), .true.,  &
                           n_compute_c, n_max_compute_fns_ham )

                      ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                      call tab_gradient_ylm_p0  &
                           ( trigonom_tab(1,1), basis_l_max,   &
                           l_ylm_max, n_compute_atoms, atom_index,  &
                           ylm_tab(1,1),   &
                           dylm_dtheta_tab(1,1),   &
                           scaled_dylm_dphi_tab(1,1)  )


                      ! and finally, assemble the actual gradients
                      call evaluate_wave_gradient_p0  &
                           ( dist_tab,  &
                           dir_tab(1,1,i_point),  &
                           trigonom_tab(1,1),  &
                           l_ylm_max, ylm_tab(1,1),  &
                           dylm_dtheta_tab(1,1),  &
                           scaled_dylm_dphi_tab(1,1),  &
                           index_lm, n_compute_c,  &
                           i_basis(1:n_compute_c),  &
                           radial_wave(1),  &
                           radial_wave_deriv(1),  &
                           gradient_basis_wave (1,1),  &
                           n_compute_atoms,  &
                           atom_index_inv,  &
                           n_compute_fns,  &
                           i_basis_fns_inv, n_max_compute_fns_ham   )

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
                        dist_tab, i_r,  &
                        atom_index, i_basis_fns_inv,  &
                        basis_kinetic_ordered, kinetic_wave(1),  &
                        .false., n_compute_c, n_max_compute_fns_ham )


                   do i_spin = 1, n_spin, 1
                      call evaluate_H_psi_p0  &
                           ( l_ylm_max,  &
                           ylm_tab(1:(l_ylm_max+1)**2,1:n_centers_integrals),  &
                           dist_tab,  &
                           index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),  &
                           H_times_psi(1, i_point, i_spin),  &
                           radial_wave(1),  &  
                                !                               local_potential_parts(i_spin),  &
                           0.d0,  &
                           n_compute_c,  &
                           i_basis(1:n_max_compute_ham),  &
                           n_compute_atoms, atom_index_inv,  &
                           n_compute_fns,   &  
                           i_basis_fns_inv,   &
                           kinetic_wave(1),  &
                           zora_operator(i_spin), n_max_compute_fns_ham )
                   enddo

                   if ((flag_rel.eq.1)) then

                      ! Scalar relativistic treatment. 
                      ! count number of "truly" relativistic points for ZORA treatment
                      ! of kinetic energy ...

                      do i_spin = 1, n_spin, 1

                         zora_operator(i_spin) =  &
                              2*light_speed_sq /  &
                              (2 * light_speed_sq -  &
                              zora_potential_parts(i_spin))**3

                         call  add_zora_gradient_part_p0(   &
                              sum_of_local_gradients(1:3,i_spin),  &
                              i_r_full,  &
                              dir_tab_full_norm(1,1),   &
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
                              light_speed_sq /  &
                              (2 * light_speed_sq -  &
                             local_potential_parts(i_spin))**2  &
                             - &
                             light_speed_sq /  &
                             (2 * light_speed_sq -  &
                             zora_potential_parts(i_spin))**2

                         call evaluate_zora_vector_p1  &
                              ( zora_operator(i_spin),  &
                              partition_tab(i_full_points_C),  &
                              gradient_basis_wave(1,1),  &
                              n_compute_c,  &
                              zora_vector1(1, 1, n_rel_points+1, i_spin),  &
                              zora_vector2(1, 1, n_rel_points+1, i_spin), &
                              n_max_compute_ham, t_zora(i_spin)  )

                      enddo

                      n_rel_points = n_rel_points + 1

                   end if  ! end ZORA preparations

                   ! If using a GGA, add the true gradient terms to the Hamiltonian vector
                   !if (use_gga .or. (n_rel_points.gt.0)) then
                   if (use_gga .or.(flag_rel.eq.1)) then

                      do i_spin = 1, n_spin, 1  
                         call add_gradient_part_to_H_p0  &
                              ( n_compute_c,   &
                              gradient_basis_wave(1,1),  &
                              sum_of_local_gradients(1,i_spin),  &
                              H_times_psi(1, i_point, i_spin) )
                      enddo
                   end if
                end if  ! end if (hamiltonian_partition_tab.gt.0)
             enddo ! end loop over a batch of points


             ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
             ! work separately for each spin channel
             do i_spin = 1, n_spin, 1

                ! add full non-relativistic contributions and (for relativistic points)
                ! all contributions from the potential to the Hamiltonian matrix elements
                call evaluate_hamiltonian_shell_p1  &
                     ( n_points, partition(1), n_compute_c, &
                     H_times_psi(1,1,i_spin),  &
                     n_max_compute_ham, wave(1,1),  &
                     hamiltonian_shell )

                ! For all relativistic points, add kinetic energy contributions
                if (n_rel_points.gt.0) then

                   call add_zora_matrix_p1 ( &
                        zora_vector1(1,1,1,i_spin),   &
                        zora_vector2(1,1,1,i_spin),   &
                        n_max_compute_ham, n_rel_points,  &
                        n_compute_c, &
                        hamiltonian_shell )


                end if


                call update_full_matrix_p0(  &
                     n_compute_c, n_compute_a, i_basis(1), hamiltonian_shell,    &
                     SZ_hamiltonian(1,i_spin) )


             enddo

             ! Hamiltonian is now complete.
             !
             ! Since we already have the pieces, add terms of XC energy here. 
             ! Notice that these terms are not added for ANY shell
             ! where n_compute happens to be zero. This should be correct because all wave functions
             ! are zero here anyway, i.e. also the density.

          else
             do i_index = 1, batches(i_my_batch)%size, 1
                i_full_points_C = i_full_points_C + 1
                i_full_points_A = i_full_points_A + 1
             enddo


          end if ! end if (n_compute.gt.0) then
       ! end if ! end mpi work distribution
    enddo ! end loop over batches

    !     synchronise the hamiltonian
    if(.not. use_local_index) then
       call sync_integrate_hamiltonian(  SZ_hamiltonian )
    endif

    if (allocated(gradient_basis_wave)) then
       deallocate (gradient_basis_wave)
    end if
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

    if (allocated(hamiltonian_shell)) then
       deallocate( hamiltonian_shell ) 
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
    if (allocated(zora_vector1)) then
       deallocate(zora_vector1)
    end if
    if (allocated(zora_vector2)) then
       deallocate(zora_vector2)
    end if
    if(allocated(H_times_psi))then
       deallocate(H_times_psi)
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
    if(allocated(wave))then
       deallocate(wave)
    end if
    if( allocated(i_basis))then
       deallocate(i_basis)
    end if

    if(allocated( rho_inc_partialcore  )) then 
       deallocate( rho_inc_partialcore  )
    endif

    if(allocated( rho_gradient_inc_partialcore )) then 
       deallocate( rho_gradient_inc_partialcore )
    endif

  end subroutine integrate_scaled_zora_transf_p2
!******

!---------------------------------------------------------------------------------
!****s* scaled_zora_transform/int_scaled_zora_tra_atomic_zora
!  NAME
!    int_scaled_zora_tra_atomic_zora
!  SYNOPSIS

  subroutine int_scaled_zora_tra_atomic_zora( rho, rho_gradient, hartree_potential,    &
       partition_tab, basis_l_max )

!  PURPOSE
!
!    The routine integrates the scaled zora integrals using so called atomic zora apprizimation.
!    The results are save to the results to SZ_hamiltonian.
!    After this the actual Scaled zora transformation can be do to the eigenvalues by calling
!    the subroutine evaluate_scaled_zora_transf_p1.
!
!    Note: the atomic zora approximation has not given a good results in our tests.
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

    real*8, dimension(n_spin, n_full_points)    :: rho
    real*8, dimension(3, n_spin, n_full_points) :: rho_gradient
    real*8, dimension(n_full_points)            :: hartree_potential
    real*8, dimension( n_full_points)           :: partition_tab
    integer, dimension(n_species)               :: basis_l_max

!  INPUTS
!
!   o rho  -- electron density
!   o rho_gradient  -- electron density gradient, 
!                     These should only ever be referenced if (use_gga)
!                     import dimensions from above (if not used, all dimensions=1)
!   o hartree_potential  -- Hartree potential
!   o partition_tab  -- Partition function
!   o basis_l_max    -- The maximum value of l of the basis
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCCE




    real*8, dimension(n_spin) :: local_potential_parts

    integer :: l_ylm_max
    integer, dimension(:,:), allocatable :: index_lm
    real*8, dimension(:,:), allocatable :: ylm_tab

    real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
    real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab

    real*8 coord_current(3)
    real*8 dist_tab(n_centers_integrals)
    real*8 dist_tab_sq(n_centers_integrals, n_max_compute_ang)
    real*8 i_r(n_centers_integrals)
    real*8 dir_tab(3,n_centers_integrals, n_max_compute_ang)
    real*8 trigonom_tab(4,n_centers_integrals)

    !paula


    real*8,dimension(:,:,:),allocatable:: H_times_psi
    real*8,dimension(:)  ,allocatable:: radial_wave
    real*8,dimension(:)  ,allocatable:: radial_wave_deriv
    real*8,dimension(:)  ,allocatable:: kinetic_wave
    real*8,dimension(:,:)  ,allocatable:: wave


    !     for XC
    real*8 :: en_density_xc
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8, dimension(n_spin) :: local_xc_derivs
    real*8, dimension(3,n_spin) :: xc_gradient_deriv

    ! Unused?
    ! real*8, dimension(n_spin) :: local_rho
    ! real*8, dimension(3,n_spin) :: local_rho_gradient


    !     Auxiliary Hamiltonian matrix, to sum up contributions from only a single integration shell
    !     The hope is that such a separate treatment will allow to minimize numerical noise
    !     introduced through ZORA
    real*8, dimension(:,:), allocatable :: hamiltonian_shell

    !     optimal accounting for matrix multiplications: only use points with nonzero components
    integer :: n_points
    integer :: n_rel_points

    !     and condensed version of hamiltonian_partition_tabs on angular grids
    real*8 :: partition(n_max_compute_ang)


    real*8, dimension(:,:), allocatable :: gradient_basis_wave

    !     Following is all that is needed for the handling of ZORA scalar relativity

    real*8, dimension(n_spin) :: zora_operator
    logical :: t_zora(2)
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

    !     for splitting of angular shells into "octants"

    integer division_low
    integer division_high

    !  counters

    integer i_basis_1
    integer i_basis_2
    integer i_grid
    integer i_index, i_l, i_m
    integer i_coord
    integer i_division

    integer i_species

    integer i_point, info

    integer :: i_my_batch

    integer :: i_full_points_C
    integer :: i_full_points_2C

    integer :: i_full_points_A
    !  integer :: i_full_points_2A

    integer :: i_spin
    character*100 :: info_str

    !test
    !  integer i_compute_1, i_compute_2, i_state, i_k_point, info

    !test end
    !  integer,dimension(:,:),allocatable :: index_b

    !  begin work

    t_zora = .false.

    write(info_str,'(2X,A,A)') & 
         "Computing correction matrix for scaled atomic-ZORA eigenvalues."
    call localorb_info(info_str,use_unit,'(A)')

     l_ylm_max = l_wave_max

    !     begin with general allocations
    if ((flag_rel.eq.0).and.(.not.(use_gga))) then
       !       no gradients needed
       l_ylm_max = l_wave_max
    else if ((flag_rel.eq.1).or.(use_gga)) then
       l_ylm_max = l_wave_max
       allocate (gradient_basis_wave(n_max_compute_ham,3),stat=info)
       call check_allocation(info, 'gradient_basis_wave           ') 
       allocate( dylm_dtheta_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info )
       call check_allocation(info, 'dylm_dtheta_tab               ') 
       allocate( scaled_dylm_dphi_tab( (l_ylm_max+1)**2, n_centers_integrals),stat=info )
       call check_allocation(info, 'scaled_dylm_dphi_tab          ')  

   end if

    allocate( ylm_tab( (l_ylm_max+1)**2, n_centers_integrals ),stat=info )
    call check_allocation(info, 'ylm_tab                          ') 
    allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info ) 
    call check_allocation(info, 'index_lm                         ') 
    allocate ( hamiltonian_shell(n_max_compute_ham, n_max_compute_ham),stat=info ) 
    call check_allocation(info, 'hamiltonian_shell                ') 
    allocate(H_times_psi(n_max_compute_ham, n_max_compute_ang, n_spin),stat=info)
    call check_allocation(info, 'H_times_psi                      ') 
    allocate(radial_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave                      ') 
    allocate(radial_wave_deriv(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'radial_wave_deriv                ') 
    allocate(kinetic_wave(n_max_compute_fns_ham),stat=info)
    call check_allocation(info, 'kinetic_wave                     ') 
    allocate(wave(n_max_compute_ham, n_max_compute_ang),stat=info)
    call check_allocation(info, 'wave                             ') 
    allocate(i_basis(n_centers_basis_I),stat=info)
    call check_allocation(info, 'i_basis                          ') 

!    if (flag_rel.eq.1) then
       ! allocate all arrays relevant for ZORA

       if (.not.allocated(dist_tab_full)) then
          allocate(dist_tab_full(n_centers_integrals),stat=info)
    call check_allocation(info, 'dist_tab_full                    ') 
       end if
       if (.not.allocated(dir_tab_full_norm)) then
          allocate(dir_tab_full_norm(3,n_centers_integrals),stat=info)
    call check_allocation(info, 'dir_tab_full_norm                ') 
       end if
       if (.not.allocated(i_r_full)) then
          allocate(i_r_full(n_centers_integrals),stat=info)
    call check_allocation(info, 'i_r_full                         ') 
       end if

       if (.not.allocated(zora_vector1)) then
          allocate(zora_vector1(n_max_compute_ham,3,n_max_compute_ang,n_spin),stat=info)
    call check_allocation(info, 'zora_vector1                     ') 
       end if
       if (.not.allocated(zora_vector2)) then
          allocate(zora_vector2(n_max_compute_ham,3,n_max_compute_ang,n_spin),stat=info)
    call check_allocation(info, 'zora_vector2                     ') 
       end if

!    end if

    SZ_hamiltonian = 0.d0




    !     initialize index_lm
    i_index = 0
    do i_l = 0, l_wave_max, 1
       do i_m = -i_l, i_l
          i_index = i_index+1
          index_lm(i_m,i_l) = i_index
       enddo
    enddo



    i_full_points_C = 0
    i_full_points_2C = 0
    i_full_points_A = 0

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


             i_full_points_2C = i_full_points_2C + 1


             if (partition_tab(i_full_points_2C).gt.0.d0) then

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
                   call prune_basis_p0 &
                        ( dist_tab_sq(1,i_point), &
                        n_compute_a, n_compute_c, i_basis,  &
                        n_centers_basis_I, n_centers_integrals, inv_centers_basis_integrals  )
                end if
             end if
          enddo  ! end loop over a batch of points

          if (prune_basis_once) then
             n_compute_a = batches(i_my_batch)%batch_n_compute
             n_compute_c = n_compute_a
             i_basis(1:n_compute_a) = batches(i_my_batch)%batch_i_basis
          end if
          n_points = i_point

  !        write(use_unit,*) 'eka 3'


          ! Perform actual integration if more than 0 basis functions
          ! are actually relevant on the present angular shell ...
          if (n_compute_c.gt.0) then

             n_rel_points = 0
             i_point = 0

             do i_index = 1, batches(i_my_batch)%size, 1

                ! Increment the (global) counter for the grid, to access storage arrays
                i_full_points_C = i_full_points_C + 1
                i_full_points_A = i_full_points_A + 1

                if (partition_tab(i_full_points_C).gt.0.d0) then

                   i_point = i_point+1

!                   call tab_global_geometry_p0 &
!                        ( dist_tab_sq(1,i_point), &
!                        dir_tab(1,1,i_point), &
!                        dist_tab_full, &
!                        i_r_full, &
!                        dir_tab_full_norm, &
!                        n_centers_integrals,  centers_basis_integrals)


                   ! for all integrations
                   partition(i_point) = partition_tab(i_full_points_C)
                   ! energy_partition(i_point) = partition_tab(i_full_points_A)

                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                   ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                   ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                   ! without any copying and without doing any unnecessary operations. 
                   ! The price is that the interface is no longer explicit in terms of physical 
                   ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.


                   call prune_radial_basis_p0 &
                        ( dist_tab_sq(1,i_point), &
                        dist_tab, &
                        dir_tab(1,1,i_point), &
                        n_compute_atoms, atom_index, atom_index_inv, &
                        n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                        i_atom_fns, spline_array_start, spline_array_end, &
                        n_centers_integrals, centers_basis_integrals)

                   ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                   ! for all atoms which are actually relevant
                   call tab_local_geometry_p0 &
                        ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                        dir_tab(1,1,i_point), dist_tab,  &
                        i_r &
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


                   ! Now evaluate radial functions
                   ! from the previously stored compressed spline arrays  
                   call evaluate_radial_functions_p0  &
                        (   spline_array_start, spline_array_end,  &
                        n_compute_atoms, n_compute_fns,   &
                        dist_tab, i_r,  &
                        atom_index, i_basis_fns_inv,  &
                        basis_wave_ordered, radial_wave(1),  &
                        .false. , n_compute_c, n_max_compute_fns_ham )

                   ! tabulate total wave function value for each basis function
                   call evaluate_waves_p0  &
                        ( l_ylm_max,   &
                        ylm_tab(1,1), dist_tab,   &
                        index_lm, n_compute_c,   &
                        i_basis, radial_wave(1),   &
                        wave(1,i_point), n_compute_atoms,   &
                        atom_index_inv, n_compute_fns,  &
                        i_basis_fns_inv,  n_max_compute_fns_ham )




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
                        dist_tab, i_r,  &
                        atom_index, i_basis_fns_inv,  &
                        basis_kinetic_scaled_zora_ordered, kinetic_wave(1),  &
                        .false., n_compute_c, n_max_compute_fns_ham )


                   do i_spin = 1, n_spin, 1
                      call evaluate_H_psi_p0  &
                           ( l_ylm_max,  &
                           ylm_tab(1:(l_ylm_max+1)**2,1:n_centers_integrals),  &
                           dist_tab,  &
                           index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),  &
                           H_times_psi(1, i_point, i_spin),  &
                           radial_wave(1),  &  
                                !                               local_potential_parts(i_spin),  &
                           0.d0,  &
                           n_compute_c,  &
                           i_basis(1:n_max_compute_ham),  &
                           n_compute_atoms, atom_index_inv,  &
                           n_compute_fns,   &  
                           i_basis_fns_inv,   &
                           kinetic_wave(1),  &
                           zora_operator(i_spin), n_max_compute_fns_ham )
                   enddo

                end if  ! end if (hamiltonian_partition_tab.gt.0)
             enddo ! end loop over a batch of points





             ! Now add all contributions to the full Hamiltonian, by way of matrix multiplications
             ! work separately for each spin channel
             do i_spin = 1, n_spin, 1

                ! add full non-relativistic contributions and (for relativistic points)
                ! all contributions from the potential to the Hamiltonian matrix elements
                call evaluate_hamiltonian_shell_p1  &
                     ( n_points, partition(1), n_compute_c, &
                     H_times_psi(1,1,i_spin),  &
                     n_max_compute_ham, wave(1,1),  &
                     hamiltonian_shell )


                call update_full_matrix_p0(  &
                     n_compute_c, n_compute_a, i_basis(1), hamiltonian_shell,    &
                     SZ_hamiltonian(1,i_spin) )
             enddo

             ! Hamiltonian is now complete.
             !
             ! Since we already have the pieces, add terms of XC energy here. 
             ! Notice that these terms are not added for ANY shell
             ! where n_compute happens to be zero. This should be correct because all wave functions
             ! are zero here anyway, i.e. also the density.

          else
             do i_index = 1, batches(i_my_batch)%size, 1
                i_full_points_C = i_full_points_C + 1
                i_full_points_A = i_full_points_A + 1
             enddo

          end if ! end if (n_compute.gt.0) then
       ! end if ! end mpi work distribution
    enddo ! end loop over batches


    !     synchronise the hamiltonian
    if(.not. use_local_index) then
       call sync_integrate_hamiltonian(  SZ_hamiltonian )
    endif

    if (allocated(gradient_basis_wave)) then
       deallocate (gradient_basis_wave)
    end if
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

    if (allocated(hamiltonian_shell)) then
       deallocate( hamiltonian_shell ) 
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
    if (allocated(zora_vector1)) then
       deallocate(zora_vector1)
    end if
    if (allocated(zora_vector2)) then
       deallocate(zora_vector2)
    end if
    if(allocated(H_times_psi))then
       deallocate(H_times_psi)
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
    if(allocated(wave))then
       deallocate(wave)
    end if
    if( allocated(i_basis))then
       deallocate(i_basis)
    end if
 
  end subroutine int_scaled_zora_tra_atomic_zora
!******
!-----------------------------------------------------------------------------

!****s* scaled_zora_transform/evaluate_scaled_zora_transf_p1
!  NAME
!    evaluate_scaled_zora_transf_p1
!  SYNOPSIS

  subroutine evaluate_scaled_zora_transf_p1( &
       KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue)


!  PURPOSE
!
!    The subroutine applies the "scaled ZORA" transformation to the eigenvalues.
!    This is post processing routine. Before calling this, the subroutine
!    integrate_scaled_zora_transf_p2 (or int_scaled_zora_tra_atomic_zora
!    if atomic zora is in use) must have been called.
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
    use scalapack_wrapper
    use localorb_io
    use KH_core_states
    implicit none

!  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points_task) :: KS_eigenvector_complex
    real*8,     dimension(n_states, n_spin, n_k_points)              :: KS_eigenvalue

!  INPUTS
!    o KS_eigenvector -- KS-eigenvectors if real eigenvector format is in use
!    o KS_eigenvector_complex -- KS-eigenvectors if complex eigenvector format is in use
!    o KS_eigenvalue -- KS-eigenvalues before scaled zora tranformation
!  OUTPUT
!    o KS_eigenvalue -- KS-eigenvalues after scaled zora tranformation
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



    real*8:: scaled_zora_integral

    real*8,    dimension(:,:),allocatable :: hamiltonian_w
    complex*16,dimension(:,:),allocatable :: hamiltonian_w_complex

    integer :: i_basis_1, i_basis_2, i_index, i_spin, info
    integer :: i_compute_1, i_compute_2, i_state, i_k_point, i_k

    !test end
    integer,dimension(:,:),allocatable :: index_b
    character*100 :: info_str
    real*8, dimension(:,:,:),allocatable :: work_ham

    ! n_core_states renamed with a suffix "_loc" as the original name conflicts
    ! with a module variable in the "dimensions" module. Jan 22, 2018
    integer:: n_core_states_loc


    if(output_priority .le. OL_norm) then
       write(info_str,'(2X,A,A)') & 
           "Evaluating scaled zora correction for each Kohn-Sham eigenstate."
       call localorb_info(info_str,use_unit,'(A)')
    endif

    if(flag_KH_core_states)then
       n_core_states_loc = n_KH_core_states
    else
       n_core_states_loc = 0 
    end if


    if(use_scalapack)then
      call evaluate_scaled_zora_tra_scalapack( KS_eigenvalue, SZ_hamiltonian )
      return
   end if




    allocate(index_b(n_basis,n_basis),stat=info)
    call check_allocation(info, 'index_b                       ')
    index_b = 0
    i_index = 0
    do i_basis_2 = 1,n_basis, 1
       do i_basis_1 = 1,i_basis_2,1
          i_index = i_index + 1
          index_b(i_basis_1, i_basis_2) = i_index
       end do
    end do


    if(n_periodic .eq. 0 .and. packed_matrix_format == PM_none)then

       ! VB: FIXME - not parallel!

       i_k_point = 1

       do i_spin = 1, n_spin, 1

          do i_state = 1+n_core_states_loc, n_states

             scaled_zora_integral = 0.0d0

             do i_compute_2 = 1,n_basis, 1
                do i_compute_1 = 1, n_basis,1


                   i_index =  max(index_b( i_compute_1, i_compute_2),index_b( i_compute_2, i_compute_1))

                   if(i_index > 0)then

                      scaled_zora_integral =  scaled_zora_integral &
                           + ( (KS_eigenvector(i_compute_1, i_state, i_spin,1)) &
                           * (KS_eigenvector(i_compute_2, i_state, i_spin,1)) &
                           *  SZ_hamiltonian(i_index,i_spin))
                   end if


                end do
             end do

             KS_eigenvalue(i_state,i_spin,i_k_point) =  KS_eigenvalue(i_state,i_spin,i_k_point)/ &
                  (1+ scaled_zora_integral)

          end do
       end do




    elseif(real_eigenvectors)then


       allocate(hamiltonian_w   (n_basis*(n_basis+1)/2,n_spin),stat=info)
       call check_allocation(info, 'hamiltonian_w                 ')
       ! complex not needed, dummy only for call below
       allocate(hamiltonian_w_complex   (1,1),stat=info)
       call check_allocation(info, 'hamiltonian_w_complex         ')
       ! work_ham not needed, dummy only for call below
       allocate(work_ham(1, 1, 1),stat=info)
       call check_allocation(info, 'work_ham                      ') 

       i_k = 0
       do i_k_point = 1, n_k_points, 1

          if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
             i_k = i_k + 1


             call construct_hamiltonian(SZ_hamiltonian, hamiltonian_w, &
                  hamiltonian_w_complex, i_k_point, work_ham )


             do i_spin = 1, n_spin, 1

                do i_state = 1+n_core_states_loc, n_states

                   scaled_zora_integral = 0.d0

                   do i_compute_2 = 1,n_basis, 1
                      do i_compute_1 = 1, n_basis,1


                         i_index =  max(index_b( i_compute_1, i_compute_2),index_b( i_compute_2, i_compute_1))

                         if(i_index > 0)then

                            scaled_zora_integral =  scaled_zora_integral &
                                 + ( (KS_eigenvector(i_compute_1, i_state, i_spin,i_k)) &
                                 * (KS_eigenvector(i_compute_2, i_state, i_spin,i_k)) &
                                 *  hamiltonian_w(i_index,i_spin))

                         end if


                      end do
                   end do

                   KS_eigenvalue(i_state,i_spin,i_k_point) =  KS_eigenvalue(i_state,i_spin,i_k_point)/ &
                        (1+ scaled_zora_integral)

                end do
             end do

          else

             KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0

          end if
       end do

       deallocate(hamiltonian_w)
       deallocate(hamiltonian_w_complex)
       deallocate(work_ham)

    else ! complex eigenvectors


       allocate(hamiltonian_w_complex   (n_basis*(n_basis+1)/2,n_spin),stat=info)
       call check_allocation(info, 'hamiltonian_w_complex         ') 
       ! real array not needed, dummy only for call below
       allocate(hamiltonian_w   (1,1),stat=info)
       call check_allocation(info, 'hamiltonian_w                 ') 

       if(packed_matrix_format == PM_none)then
          allocate(work_ham(1:n_centers_basis_I, 1:n_centers_basis_I, n_spin),stat=info)
          call check_allocation(info, 'work_ham                      ') 
       else 
          ! work_ham only needed for non-packed case, here: dummy only for call below
          allocate(work_ham(1, 1, 1),stat=info)
          call check_allocation(info, 'work_ham                      ') 
       end if


       i_k = 0
       do i_k_point = 1, n_k_points, 1

          if (myid.eq.  MOD(i_k_point, n_tasks) .and. myid <= n_k_points ) then
             i_k = i_k + 1


             call construct_hamiltonian(SZ_hamiltonian, hamiltonian_w, &
                  hamiltonian_w_complex, i_k_point, work_ham )



             do i_spin = 1, n_spin, 1

                do i_state = 1+n_core_states_loc, n_states

                   scaled_zora_integral = 0.d0

                   do i_compute_2 = 1,n_basis, 1
                      do i_compute_1 = 1, n_basis,1


                         i_index =  max(index_b( i_compute_1, i_compute_2),index_b( i_compute_2, i_compute_1))

                         if(i_index > 0)then

                            if(i_compute_1<=i_compute_2) then


                               scaled_zora_integral =  scaled_zora_integral &
                                    + dble( dconjg(KS_eigenvector_complex(i_compute_1, i_state, i_spin,i_k)) &
                                    * KS_eigenvector_complex(i_compute_2, i_state, i_spin,i_k) &
                                    *  hamiltonian_w_complex(i_index,i_spin))

                            else

                               scaled_zora_integral =  scaled_zora_integral &
                                    + dble( dconjg(KS_eigenvector_complex(i_compute_1, i_state, i_spin,i_k)) &
                                    * KS_eigenvector_complex(i_compute_2, i_state, i_spin,i_k)  &
                                    *  dconjg(hamiltonian_w_complex(i_index,i_spin)))

                            end if
                         end if


                      end do
                   end do


                   KS_eigenvalue(i_state,i_spin,i_k_point) =  KS_eigenvalue(i_state,i_spin,i_k_point)/ &
                        (1+ scaled_zora_integral)

                end do
             end do
          else

             KS_eigenvalue(1:n_states,1:n_spin,i_k_point) =  0.d0


          end if
       end do

       deallocate(hamiltonian_w_complex)
       deallocate(hamiltonian_w)
       deallocate(work_ham)

    end if



    ! FIXME - Patch by VB: Cluster case is not parallel, should be!
    if (n_periodic.ne.0 .or. packed_matrix_format==PM_index) then
      call sync_eigenvalues( KS_eigenvalue)
    end if  


  end subroutine evaluate_scaled_zora_transf_p1
!******
!--------------------------------------------------------------------------------
!****s* scaled_zora_transform/relativistic_correction_term_p1
!  NAME
!    relativistic_correction_term_p1
!  SYNOPSIS

  subroutine relativistic_correction_term_p1 &
       ( KS_eigenvector, KS_eigenvector_complex, KS_eigenvalue, occ_numbers, partition_tab,  &
       rho, rho_gradient, kinetic_density, hartree_potential, basis_l_max)

!  PURPOSE
!
!    The subroutine have been used for testing relativistic post process method. 
!    The method is not currently supported.
!
!  USES

    use dimensions
    use runtime_choices
    use grids
    use geometry
    use mixing
    use species_data
    use mpi_utilities
    use synchronize_mpi
    use localorb_io
    use basis
    use cartesian_ylm
    use constants
    use pbc_lists
    implicit none

!  ARGUMENTS

    real*8,     dimension(n_basis, n_states, n_spin,n_k_points):: KS_eigenvector
    complex*16, dimension(n_basis, n_states, n_spin,n_k_points):: KS_eigenvector_complex
    real*8,     dimension(n_states, n_spin, n_k_points)        :: KS_eigenvalue      
    real*8,     dimension(n_states, n_spin, n_k_points)        :: occ_numbers
    real*8,     dimension(n_full_points)                       :: partition_tab
    real*8,     dimension(n_spin, n_full_points)               :: rho
    real*8,     dimension(3, n_spin, n_full_points)            :: rho_gradient
    real*8,     dimension(n_spin, n_full_points)               :: kinetic_density
    real*8,     dimension(n_full_points)                       :: hartree_potential
    integer :: basis_l_max (n_species)

!  INPUTS
!
!   o KS_eigenvector -- KS-eigenvectors if real eigenvector format is in use
!   o KS_eigenvector_complex -- KS-eigenvectors if complex eigenvector format is in use
!   o KS_eigenvalue -- KS-eigenvalues before post process
!   o occ_numbers -- occupation weights of KS-eigenstates
!   o partition_tab  -- Partition function
!   o rho  -- electron density
!   o rho_gradient  -- electron density gradient, 
!                     These should only ever be referenced if (use_gga)
!                     import dimensions from above (if not used, all dimensions=1)
!   o kinetic_density -- Kinetic-energy density. Only used for meta_gga
!   o hartree_potential  -- Hartree potential
!   o basis_l_max    -- The maximum value of l of the basis
!
!  OUTPUT
!
!   o KS_eigenvalue -- KS-eigenvalues after post process
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE






    !  imported variables

    !  input 




    real*8 coord_current(3)
    real*8 dist_tab(n_centers_basis_integrals  , n_max_batch_size)
    real*8 dist_tab_sq(n_centers_basis_integrals, n_max_batch_size)
    real*8 i_r(n_centers_basis_integrals)
    real*8 dir_tab(3, n_centers_basis_integrals, n_max_batch_size)
    real*8 trigonom_tab(4, n_centers_basis_integrals)

    real*8,dimension(:),allocatable:: radial_wave
    real*8,dimension(:),allocatable:: radial_wave_deriv
    real*8,dimension(:,:),allocatable:: wave
    real*8 :: local_potential_parts(n_spin)
    real*8 :: xc_gradient_deriv(3, n_spin, n_max_batch_size)
    real*8 :: xc_tau_deriv(n_spin)
    real*8 :: en_density_xc(n_max_batch_size)
    real*8, dimension(n_spin) :: en_density_x
    real*8 :: en_density_c
    real*8 :: local_xc_derivs(n_spin, n_max_batch_size)

    integer :: n_compute_c, n_compute_a
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

    integer :: l_ylm_max
    integer :: n_points
    integer::  i_compute, i_compute2
    logical :: write_out
    integer, dimension(:,:), allocatable :: index_lm
    real*8, dimension(:,:), allocatable :: ylm_tab

    real*8, dimension(:,:), allocatable :: dylm_dtheta_tab
    real*8, dimension(:,:), allocatable :: scaled_dylm_dphi_tab
    real*8, dimension(:,:,:), allocatable :: gradient_basis_wave

    integer i_atom
    integer i_l, info
    integer i_m
    integer i_coord
    integer :: i_state
    integer :: i_point

    integer :: i_my_batch

    integer :: i_full_points
    integer :: i_full_points_2
    integer :: i_full_points_3


    integer :: i_spin = 1
    integer ::  i_index, i_k_point, i_ind

    real*8:: kinetic
    real*8, dimension(:,:), allocatable :: KS_eigenvalue_post

    real*8, dimension(:,:), allocatable  :: rho_inc_partialcore
    real*8, dimension(:,:,:), allocatable  :: rho_gradient_inc_partialcore



    !     begin work


    call localorb_info( &
         "Relativistic KH energy correction", 6,'(2X,A)' )


    if (.not.allocated(KS_eigenvalue_post)) then
       allocate( KS_eigenvalue_post(n_states,n_spin),stat=info )
       call check_allocation(info, 'KS_eigenvalue_post            ')
    end if
    KS_eigenvalue_post = 0.d0

    if(.not. allocated(i_basis))then
       allocate(i_basis(n_centers_basis_T),stat=info)
       call check_allocation(info, 'i_basis                       ') 
    end if

    if(.not.allocated(radial_wave))then
       allocate(radial_wave(n_max_compute_fns_dens),stat=info)
       call check_allocation(info, 'radial_wave                   ')
    end if

    if(.not. allocated(radial_wave_deriv))then
       allocate(radial_wave_deriv(n_max_compute_fns_dens),stat=info)
       call check_allocation(info, 'radial_wave_deriv             ') 
    end if
    if(.not. allocated(wave))then
       allocate(wave(n_max_compute_dens, n_max_batch_size),stat=info)
       call check_allocation(info, 'wave                          ') 
    end if

    l_ylm_max = l_wave_max

    if(.not. allocated( ylm_tab))then
       allocate( ylm_tab( (l_ylm_max+1)**2,n_centers_basis_integrals),stat=info )
       call check_allocation(info, 'ylm_tab                       ') 
    end if

    if(.not. allocated( index_lm))then
       allocate( index_lm( -l_ylm_max:l_ylm_max, 0:l_ylm_max),stat=info) 
       call check_allocation(info, 'index_lm                      ')
    end if

    allocate (gradient_basis_wave(n_max_compute_dens, 3, n_max_batch_size),stat=info)
    call check_allocation(info, 'gradient_basis_wave           ') 

    allocate( dylm_dtheta_tab((l_ylm_max+1)**2, n_centers_basis_integrals),stat=info )
    call check_allocation(info, 'dylm_dtheta_tab               ') 

    allocate( scaled_dylm_dphi_tab((l_ylm_max+1)**2,n_centers_basis_integrals),stat=info )
    call check_allocation(info, 'scaled_dylm_dphi_tab          ') 


    if(use_embedding_pp.and.use_nonlinear_core) then
        allocate(rho_inc_partialcore(n_spin,n_full_points))
        allocate(rho_gradient_inc_partialcore(3,n_spin,n_full_points))
        do i_spin = 1,n_spin
           rho_inc_partialcore(i_spin,:) = rho(i_spin,:) + partial_core_rho(:)
           rho_gradient_inc_partialcore(:,i_spin,:) = rho_gradient(:,i_spin,:)
           if(use_density_gradient) then
              rho_gradient_inc_partialcore(:,i_spin,:) = &
                 rho_gradient_inc_partialcore(:,i_spin,:) + partial_core_rho_grad(:,:)
           endif
        enddo
    endif



    !     initialize index_lm
    i_index = 0
    do i_l = 0, l_ylm_max, 1
       do i_m = -i_l, i_l
          i_index = i_index + 1
          index_lm(i_m, i_l) = i_index
       enddo
    enddo


    if (use_density_gradient) then
       delta_rho_gradient = 0.0d0
    end if

    i_full_points = 0
    i_full_points_2 = 0
    i_full_points_3 = 0
    i_k_point = 1
    write_out = .false.




    do i_my_batch = 1, n_my_batches, 1

          !       write(use_unit,*) i_batch,  energy_correction*hartree
          !       write(use_unit,*) i_batch !,  energy_correction*hartree

          n_compute_c = 0
          n_compute_a = 0
          i_basis = 0

          i_point = 0

          ! loop over one batch
          do i_index = 1, batches(i_my_batch)%size, 1


             i_full_points = i_full_points + 1

             if (partition_tab(i_full_points).gt.0.d0) then

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

                ! next, determine which basis functions u(r)/r*Y_lm(theta,phi) are actually needed
                if (.not.prune_basis_once) then
                   call prune_basis_p0 &
                        (dist_tab_sq(1,i_point), &
                        n_compute_a, n_compute_c, i_basis,  &
                        n_centers_basis_T, n_centers_basis_integrals, inv_centers_basis_integrals )
                end if
             end if
          enddo ! end loop over the angular shell

          if (prune_basis_once) then
             n_compute_a = batches(i_my_batch)%batch_n_compute
             n_compute_c = n_compute_a
             i_basis(1:n_compute_a) = batches(i_my_batch)%batch_i_basis
          end if

          n_points = i_point

          if (n_compute_c.gt.0) then

             ! Determine all radial functions, ylm functions and their derivatives that
             ! are best evaluated strictly locally at each individual grid point.
             i_point = 0
             do i_index = 1, batches(i_my_batch)%size, 1

                i_full_points_2 = i_full_points_2 + 1

                if (partition_tab(i_full_points_2).gt.0.d0) then

                   i_point = i_point+1
                   n_compute_atoms = 0
                   n_compute_fns = 0
                   i_basis_fns_inv = 0

                   ! All radial functions (i.e. u(r), u''(r)+l(l+2)/r^2, u'(r) if needed)
                   ! Are stored in a compact spline array that can be accessed by spline_vector_waves, 
                   ! without any copying and without doing any unnecessary operations. 
                   ! The price is that the interface is no longer explicit in terms of physical 
                   ! objects. See shrink_fixed_basis() for details regarding the reorganized spline arrays.

                   !                     dir_tab(:,:,i_point,i_division) = dir_tab_full(:,:,i_angular)

                   call prune_radial_basis_p0 &
                        ( dist_tab_sq(1,i_point),  &
                        dist_tab(1,i_point), &
                        dir_tab(1,1,i_point), &
                        n_compute_atoms, atom_index, atom_index_inv, &
                        n_compute_fns, i_basis_fns, i_basis_fns_inv, &
                        i_atom_fns, spline_array_start, spline_array_end, &
                        n_centers_basis_integrals, centers_basis_integrals)

                   ! Tabulate distances, unit vectors, and inverse logarithmic grid units
                   ! for all atoms which are actually relevant
                   call tab_local_geometry_p0 &
                        ( dist_tab_sq(1, i_point), n_compute_atoms, atom_index, &
                        dir_tab(1,1,i_point), dist_tab(1,i_point),  &
                        i_r(1) &
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
                        .false., n_compute_c, n_max_compute_fns_dens   &
                        )

                   ! for forces or density gradient, radial derivatives are required 
                   call evaluate_radial_functions_p0 &
                        (spline_array_start, spline_array_end, &
                        n_compute_atoms, n_compute_fns,  &
                        dist_tab(1,i_point), i_r(1), &
                        atom_index, i_basis_fns_inv, &
                        basis_deriv_ordered,  &
                        radial_wave_deriv(1),  &
                        .true., n_compute_c, n_max_compute_fns_dens  &
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

                   !                    if (use_density_gradient .or. forces_on) then
                   ! tabulate those ylms needed for gradients, i.e. ylm's for l_max+1
                   call tab_gradient_ylm_p0 &
                        (trigonom_tab(1,1), basis_l_max,  &
                        l_ylm_max, n_compute_atoms, atom_index, &
                        ylm_tab(1,1),  &
                        dylm_dtheta_tab(1,1),  &
                        scaled_dylm_dphi_tab(1,1) &
                        )

                   ! evaluate the wave function gradient directly here                  
                   call evaluate_wave_gradient_p0 &
                        (dist_tab(1,i_point), &
                        dir_tab(1,1,i_point), &
                        trigonom_tab(1,1), &
                        l_ylm_max, ylm_tab(1,1), &
                        dylm_dtheta_tab(1,1), &
                        scaled_dylm_dphi_tab(1,1), &
                        index_lm, n_compute_c, &
                        i_basis(1:n_compute_c), &
                        radial_wave(1), &
                        radial_wave_deriv(1), &
                        gradient_basis_wave (1:n_compute_c,1:3,i_point), &
                        n_compute_atoms, &
                        atom_index_inv, &
                        n_compute_fns, &
                        i_basis_fns_inv, n_max_compute_fns_dens  )


                   ! tabulate total wave function value for each basis function in all cases -
                   ! but only now are we sure that we have ylm_tab ...

                   call evaluate_waves_p0 &
                        (l_ylm_max, ylm_tab(1,1),  &
                        dist_tab(1,i_point),  &
                        index_lm, n_compute_c,  &
                        i_basis, radial_wave(1),  &
                        wave(1,i_point), n_compute_atoms,   &
                        atom_index_inv, n_compute_fns,  &
                        i_basis_fns_inv, n_max_compute_fns_dens  &
                        )

                   !     Finally, must re-evaluate the Hamiltonian at current point if forces
                   !     are needed (for straightforward Pulay force terms ...)



                end if ! end if (partition_tab.gt.0)
             enddo ! end loop over one batch of the grid
             ! - all quantities that are evaluated pointwise are now known ...

             if (n_points.gt.0) then
                ! Now perform all operations which are done across the entire integration
                ! shell at once.
                i_point = 0
                do i_index = 1, batches(i_my_batch)%size, 1

                   i_full_points_3 = i_full_points_3 + 1

                   if (partition_tab(i_full_points_3).gt.0.d0) then

                      i_point = i_point+1

                      ! This one is required for vdW-DF
                      coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)!SAG
                      if(n_periodic > 0)then
                         call map_to_center_cell(coord_current(1:3) )
                      end if

                      if(use_embedding_pp.and.use_nonlinear_core) then
! here we might add partial core densities from pseudocores
! only for the evaluation of the xc-energy
                         call evaluate_xc  &
                              ( rho_inc_partialcore(1,i_full_points_3),   &
                              rho_gradient_inc_partialcore(1,1,i_full_points_3),  &
                              kinetic_density(1,i_full_points_3), &
                              en_density_xc(i_point),   &
                              en_density_x, en_density_c,  &
                              local_xc_derivs(1,i_point),  &
                              xc_gradient_deriv(1,1,i_point),  &
                              xc_tau_deriv(1), &
                              .false., &
                              coord_current(:) &
                              )

                      else

                         call evaluate_xc  &
                              ( rho(1,i_full_points_3),   &
                              rho_gradient(1,1,i_full_points_3),  &
                              kinetic_density(1,i_full_points_3), &
                              en_density_xc(i_point),   &
                              en_density_x, en_density_c,  &
                              local_xc_derivs(1,i_point),  &
                              xc_gradient_deriv(1,1,i_point),  &
                              xc_tau_deriv(1), &
                              .false., &
                              coord_current(:) &
                              )

                      endif

                      do i_spin = 1, n_spin, 1

                         local_potential_parts(i_spin) =   &
                              hartree_potential(i_full_points_3) +   &
                              local_xc_derivs(i_spin,i_point)


                      end do


                      !                            write(use_unit,*) 'loopit'

                      i_ind = 0
                      ! do i_k_point_g = 1, n_k_points_group

                      ! i_k_point = (i_k_point_group-1) * n_k_points_group + i_k_point_g

                      ! if(i_k_point <= n_k_points)then
                      do i_state = 1, n_states, 1
                         do i_spin = 1, n_spin, 1

                            if (occ_numbers(i_state,i_spin, i_k_point).gt.0.d0) then

                               i_ind = i_ind + 1


                               kinetic = &
                                    light_speed_sq /  &
                                    ( 2 * light_speed_sq   &
                                    - local_potential_parts(i_spin)) &
                                    - &
                                    light_speed_sq /  &
                                    ( 2 * light_speed_sq +  &
                                    KS_eigenvalue(i_state,i_spin,i_k_point) - local_potential_parts(i_spin) )


                               if(real_eigenvectors)then

                                  do i_coord = 1, 3, 1

                                     do i_compute = 1, n_compute_c, 1
                                        do i_compute2 = 1, n_compute_c, 1



                                           KS_eigenvalue_post(i_state,i_spin) = &
                                                KS_eigenvalue_post(i_state,i_spin) &
                                                - gradient_basis_wave(i_compute,i_coord,i_point)  &
                                                * gradient_basis_wave(i_compute2,i_coord,i_point) &
                                                * kinetic  &
                                                * KS_eigenvector(i_basis(i_compute),i_state,i_spin, i_k_point) &
                                                * KS_eigenvector(i_basis(i_compute2),i_state,i_spin, i_k_point) & 
                                                * partition_tab(i_full_points_3)

                                        enddo
                                     end do
                                  end do

                               else



                                  do i_coord = 1, 3, 1

                                     do i_compute = 1, n_compute_c, 1
                                        do i_compute2 = 1, n_compute_c, 1

                                           KS_eigenvalue_post(i_state,i_spin) = &
                                                KS_eigenvalue_post(i_state,i_spin) &
                                                - gradient_basis_wave(i_compute,i_coord,i_point)  &
                                                * gradient_basis_wave(i_compute2,i_coord,i_point) &
                                                * kinetic  &
          * KS_eigenvector_complex(Cbasis_to_basis(i_basis(i_compute)),i_state,i_spin, i_k_point) &
          * dconjg(KS_eigenvector_complex(Cbasis_to_basis(i_basis(i_compute2)),i_state,i_spin, i_k_point)) & 
                                                * partition_tab(i_full_points_3)


                                        end do
                                     end do
                                  end do



                               end if

                            end if
                         end do
                      end do
                      !  end if
                      !  end do




                   end if
                end do
             end if ! (n_points>0)


          else

             do i_index = 1, batches(i_my_batch)%size, 1

                i_full_points_2 = i_full_points_2 + 1
                i_full_points_3 = i_full_points_3 + 1

             enddo
          end if  ! end if (n_compute.gt.0)
       ! end if ! end distribution over threads
    end do ! end loop over batches
    !end do ! end loop over groups of k-points



    do i_k_point  = 1, n_k_points     

       KS_eigenvalue(:,:,i_k_point)  =  KS_eigenvalue(:,:,i_k_point) &
            +  KS_eigenvalue_post(:,:)

    end do



    ! finally, deallocate stuff.
    if (allocated(dylm_dtheta_tab)) then
       deallocate (dylm_dtheta_tab)
    end if
    if (allocated(scaled_dylm_dphi_tab)) then
       deallocate (scaled_dylm_dphi_tab)
    end if
    if (allocated(ylm_tab)) then
       deallocate(ylm_tab)
    end if
    if (allocated(index_lm)) then
       deallocate(index_lm)
    end if
    if (allocated(gradient_basis_wave)) then
       deallocate(gradient_basis_wave)
    end if
    if(allocated(radial_wave))then
       deallocate(radial_wave)
    end if
    if(allocated(radial_wave_deriv))then
       deallocate(radial_wave_deriv)
    end if
    if(allocated(wave))then
       deallocate(wave)
    end if
    if(allocated(i_basis))then
       deallocate(i_basis)
    end if
    if (allocated(KS_eigenvalue_post)) then
       deallocate( KS_eigenvalue_post )
    end if

    if(allocated( rho_inc_partialcore  )) then 
       deallocate( rho_inc_partialcore  )
    endif

    if(allocated( rho_gradient_inc_partialcore )) then 
       deallocate( rho_gradient_inc_partialcore )
    endif


  end subroutine relativistic_correction_term_p1
  !******
  !-------------------------------------------------------------------------



end module scaled_zora_transform
!----------------------------------------------------------------------
