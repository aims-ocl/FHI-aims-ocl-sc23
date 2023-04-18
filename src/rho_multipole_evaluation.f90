!****h* FHI-aims/rho_multipole_evaluation
!  NAME
!    rho_multipole_evaluation
!  SYNOPSIS

module rho_multipole_evaluation

!  PURPOSE 
!    This module contains several routines that evaluate
!
!     o the free atoms' superposition electron density and its gradient
!     o the multipole expanded delta density and its gradient
!     o the multipole expanded density (as sum of free + delta) and its gradient
!
!    at arbitrary points in real space.
!
!    However, there is NO check implemented whether requested evaluation
!    points lie on top of atoms / multipole centers. Thus it is the user's
!    responsibility to exclude such points when calling the following routines:
!
!     o get_rho_multipole
!     o get_rho_multipole_and_gradient
!     o get_free_rho
!     o get_free_rho_and_gradient
!     o get_delta_rho_multipole
!     o get_delta_rho_multipole_and_gradient
!
!  USES

   use constants,    only: pi4_inv
   use mpi_tasks,    only: myid, use_mpi, n_tasks
   use physics,      only: multipole_radius_sq, outer_potential_radius
   use hartree_potential_storage, only : get_rho_multipole_spl
   use dimensions,   only: l_pot_max, n_spin, n_centers_basis_integrals, &
                           n_centers_hartree_multipole, &
                           n_max_spline, n_max_radial
   use grids,        only: n_grid, n_radial
   use geometry,     only: species, empty
   use free_atoms,   only: renormalized_free_rho_spl, &
                           renormalized_free_drho_dr_spl, &
                           free_rho_spl, free_drho_dr_spl
   use spline,       only: spline_vector_v2, spline_deriv_vector_v2, &
                           val_spline, val_spline_deriv
   use pbc_lists,    only: centers_basis_integrals, &
                           centers_hartree_multipole, &
                           species_center, center_to_atom
   use species_data, only: multipole_radius_free_sq, &
                           species_pseudoized, &
                           l_hartree

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   implicit none

   ! make everything private by default
   private

   ! specifically declare public routines
   public :: get_rho_multipole
   public :: get_rho_multipole_and_gradient
   public :: get_free_rho
   public :: get_free_rho_and_gradient
   public :: get_delta_rho_multipole
   public :: get_delta_rho_multipole_and_gradient
   
contains


!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_rho_multipole_and_gradient
!  NAME
!    get_rho_multipole_and_gradient
!  SYNOPSIS
subroutine get_rho_multipole_and_gradient( &
      n_points, &
      points, &
      rho_multipole, &
      rho_multipole_gradient )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded density and its gradient
!    at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: rho_multipole
   real*8, dimension(3, n_points), intent(out) :: rho_multipole_gradient

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where density and density
!               gradient are evaluated
!  OUTPUT
!   o rho_multipole -- values of the multipole expanded density
!   o rho_multipole_gradient -- gradients of the multipole expanded density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE


   ! temporary variables
   real*8, dimension(:), allocatable :: delta_rho_multipole
   real*8, dimension(:,:), allocatable :: delta_rho_multipole_gradient
   
   ! index variables
   integer :: i_point


   ! ALLOCATION
   
   allocate(delta_rho_multipole(n_points))
   allocate(delta_rho_multipole_gradient(3,n_points))


   ! evaluate free atoms superposition
   ! store directly in output array
   call get_free_rho_and_gradient(n_points, points, &
            rho_multipole, rho_multipole_gradient)


   ! evaluate delta density
   ! store in temporary array
   call get_delta_rho_multipole_and_gradient(n_points, points, &
            delta_rho_multipole, delta_rho_multipole_gradient)


   do i_point = 1, n_points, 1

      ! add up free and delta contributions
      rho_multipole(i_point) = rho_multipole(i_point) + &
                  delta_rho_multipole(i_point)
      
      rho_multipole_gradient(:,i_point) = rho_multipole_gradient(:,i_point) + &
                  delta_rho_multipole_gradient(:,i_point)

   enddo ! i_point


   ! DEALLOCATION
   
   if (allocated(delta_rho_multipole)) &
      deallocate(delta_rho_multipole)

   if (allocated(delta_rho_multipole_gradient)) &
      deallocate(delta_rho_multipole_gradient)


end subroutine get_rho_multipole_and_gradient




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_rho_multipole
!  NAME
!    get_rho_multipole
!  SYNOPSIS
subroutine get_rho_multipole( &
      n_points, &
      points, &
      rho_multipole )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded density
!    at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: rho_multipole

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where the density is evaluated
!  OUTPUT
!   o rho_multipole -- values of the multipole expanded density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   
   ! temporary variables
   real*8, dimension(:), allocatable :: delta_rho_multipole
   
   ! index variables
   integer :: i_point


   ! ALLOCATION
   
   allocate(delta_rho_multipole(n_points))


   ! evaluate free atoms superposition
   ! save directly in output array
   call get_free_rho(n_points, points, &
            rho_multipole)


   ! evaluate delta density
   ! save in temporary array
   call get_delta_rho_multipole(n_points, points, &
            delta_rho_multipole)


   do i_point = 1, n_points, 1

      ! add up free and delta contributions
      rho_multipole(i_point) = rho_multipole(i_point) + &
                  delta_rho_multipole(i_point)
      
   enddo ! i_point


   ! DEALLOCATION
   
   if (allocated(delta_rho_multipole)) &
      deallocate(delta_rho_multipole)

end subroutine get_rho_multipole




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_free_rho_and_gradient
!  NAME
!    get_free_rho_and_gradient
!  SYNOPSIS
subroutine get_free_rho_and_gradient( &
      n_points, &
      points, &
      free_rho, &
      free_rho_gradient )
!  PURPOSE 
!    This subroutine evaluates the free electron density (density of the free 
!    atoms superposition) and its gradient at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: free_rho
   real*8, dimension(3, n_points), intent(out) :: free_rho_gradient

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where density and density
!               gradient are evaluated
!  OUTPUT
!   o free_rho -- values of the free electron density
!   o free_rho_gradient -- gradients of the free electron density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE


   ! distance and direction
   real*8, dimension(:), allocatable :: dist_tab_sq_free !(n_centers_basis_integrals)
   real*8, dimension(:), allocatable :: dist_tab_free !(n_centers_basis_integrals)
   real*8, dimension(:,:), allocatable :: dir_tab_free !(3,n_centers_basis_integrals)
   real*8, dimension(:,:), allocatable :: dir_tab_norm !(3,n_centers_basis_integrals)
   
   ! logarithmic grid shell indices
   real*8, dimension(:), allocatable :: i_r_free !(n_centers_basis_integrals)
   
   ! density gradient of free atom
   real*8 :: atomic_gradient
   
   ! helpers
   integer :: current_center_free

   ! external function
   real*8, external :: ddot
   
   ! sparse representation of multipole center indices
   integer :: n_compute_atoms
   integer, dimension(:), allocatable :: center_index !n_centers_basis_integrals

   ! index variables
   integer :: i_point
   integer :: i_center_L, i_center_L2


   ! ALLOCATION
   
   allocate(dist_tab_sq_free(n_centers_basis_integrals))
   allocate(dist_tab_free(n_centers_basis_integrals))
   allocate(dir_tab_free(3,n_centers_basis_integrals))
   allocate(dir_tab_norm(3,n_centers_basis_integrals))
   allocate(i_r_free(n_centers_basis_integrals))
   allocate(center_index(n_centers_basis_integrals))


   free_rho = 0.d0
   free_rho_gradient = 0.d0

   ! evaluate free atoms superposition
   do i_point = 1, n_points, 1

      !Begin free atom part (as in initialize_grid_storage_p1)
      call tab_atom_centered_coords_p0 &
          ( points(:,i_point), &
            dist_tab_sq_free, &
            dir_tab_free, &
            n_centers_basis_integrals, centers_basis_integrals )
      
      ! relevant atoms
      n_compute_atoms = 0

      ! get relevant centers that are within the free atom's charge radius
      do i_center_L = 1, n_centers_basis_integrals, 1
      
         current_center_free = centers_basis_integrals(i_center_L)
         
         ! within radius?
         if ( (dist_tab_sq_free(i_center_L) .lt. &
               multipole_radius_free_sq(species_center(current_center_free)) ) .and.&
         ! not empty?
              (.not. empty(center_to_atom(current_center_free))) .and. &
         ! not pseudo-atom?
              (.not. species_pseudoized(species(center_to_atom(current_center_free)))) &
         ! then go
               ) then
            n_compute_atoms = n_compute_atoms + 1
            dist_tab_sq_free(n_compute_atoms) = dist_tab_sq_free(i_center_L) 
            dir_tab_free(:,n_compute_atoms) = dir_tab_free(:,i_center_L) 
            center_index(n_compute_atoms) = current_center_free
         endif
      enddo

      call tab_global_geometry_p0 &
          ( dist_tab_sq_free, &
            dir_tab_free, &
            dist_tab_free, &
            i_r_free, &
            dir_tab_norm, &
            n_compute_atoms, &
            center_index )
      
      ! calculate the free-atom density only for the (now) known atoms ...
      do i_center_L2 = 1, n_compute_atoms, 1

         current_center_free = center_index(i_center_L2)
         
         ! sum up density
         free_rho(i_point) = free_rho(i_point) + val_spline &
                ( i_r_free(i_center_L2), &
                  renormalized_free_rho_spl(1,1,species_center(current_center_free)), &
                  n_grid(species_center(current_center_free)) )
         
         ! get absolute value of gradient
         atomic_gradient = val_spline &
                ( i_r_free(i_center_L2), &
                  renormalized_free_drho_dr_spl(1,1,species_center(current_center_free)), &
                  n_grid(species_center(current_center_free)) )
         
         ! sum up gradient
         free_rho_gradient(:,i_point) = &
               free_rho_gradient(:,i_point) + &
               atomic_gradient * dir_tab_norm(:,i_center_L2)

         ! As far as I understand, a call to evaluate_free_atom_sums_p2 is
         ! unnecessary. The free atom density and its gradient have already
         ! been calculated from the interpolants.
      
      enddo

      ! divide out factor of pi4 to get free atom density and gradient
      free_rho(i_point) = free_rho(i_point) * pi4_inv
      free_rho_gradient(:,i_point) = &
                  free_rho_gradient(:,i_point) * pi4_inv
   
   enddo ! i_point


   ! DEALLOCATION
   
   if (allocated(dist_tab_sq_free)) &
      deallocate(dist_tab_sq_free)

   if (allocated(dist_tab_free)) &
      deallocate(dist_tab_free)

   if (allocated(dir_tab_free)) &
      deallocate(dir_tab_free)

   if (allocated(dir_tab_norm)) &
      deallocate(dir_tab_norm)

   if (allocated(i_r_free)) &
      deallocate(i_r_free)

   if (allocated(center_index)) &
      deallocate(center_index)


end subroutine get_free_rho_and_gradient




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_delta_rho_multipole_and_gradient
!  NAME
!    get_delta_rho_multipole_and_gradient
!  SYNOPSIS
subroutine get_delta_rho_multipole_and_gradient( &
      n_points, &
      points, &
      delta_rho_multipole, &
      delta_rho_multipole_gradient )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded delta electron density 
!    (i.e. the difference between the multipole expanded density and the free 
!    atoms superposition) and its gradient at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: delta_rho_multipole
   real*8, dimension(3, n_points), intent(out) :: delta_rho_multipole_gradient

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where density and density
!               gradient are evaluated
!  OUTPUT
!   o delta_rho_multipole -- values of the multipole expanded delta density
!   o delta_rho_multipole_gradient -- gradients of the multipole expanded
!                                     delta density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   
   ! distance and direction
   real*8 :: dist_tab_sq
   real*8 :: dist_tab_in
   real*8, dimension(3) :: dir_tab
   real*8, dimension(3) :: dir_tab_in
   
   ! trigonometric function values
   real*8, dimension(4) :: trigonom_tab
   real*8, dimension(:), allocatable :: ylm_tab !(l_pot_max+1)**2
   real*8, dimension(:), allocatable :: dylm_dtheta_tab !(l_pot_max+1)**2
   real*8, dimension(:), allocatable :: scaled_dylm_dphi_tab !(l_pot_max+1)**2
   
   ! temporary variables
   real*8, dimension(3) :: rho_multipole_gradient_temp
   
   ! logarithmic grid shell indices
   real*8 :: i_r
   real*8 :: i_r_log
   
   ! infinitesimal weights
   real*8 :: log_weight
   real*8 :: radial_weight
   
   ! interpolated and partitioned representation of delta density
   real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl !((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
   real*8, dimension(:), allocatable :: rho_multipole_component !(l_pot_max+1)**2
   real*8, dimension(:), allocatable :: rho_multipole_deriv !(l_pot_max+1)**2
   
   ! array dimensions
   integer :: l_p_dim ! (l_pot_max+1)**2
   integer :: l_atom_max
   integer :: l_h_dim ! (l_atom_max+1)**2, no constant

   ! helpers
   integer :: atom_of_splines
   integer :: current_spl_atom
   integer :: current_center

   ! external function
   real*8, external :: ddot
   
   ! index variables
   integer :: i_center
   integer :: i_point


   l_p_dim = (l_pot_max + 1)**2

   ! ALLOCATION
   
   allocate(ylm_tab(l_p_dim))
   allocate(dylm_dtheta_tab(l_p_dim))
   allocate(scaled_dylm_dphi_tab(l_p_dim))
   allocate(current_rho_multipole_spl(l_p_dim, n_max_spline, n_max_radial+2))
   allocate(rho_multipole_component(l_p_dim))
   allocate(rho_multipole_deriv(l_p_dim))


   delta_rho_multipole = 0.d0
   delta_rho_multipole_gradient = 0.d0

   ! this variable indicates for which center the density is currently
   ! splined. set to 0 in the beginning to force calculation in the first step.
   atom_of_splines = 0

   ! sum up contribution of multipole centers
   do i_center = 1, n_centers_hartree_multipole, 1

      ! translate center index
      current_center = centers_hartree_multipole(i_center)
      current_spl_atom = center_to_atom(current_center)

      ! check if the current atom has already been splined or
      ! if an update is necessary
      if (current_spl_atom .ne. atom_of_splines) then
         call get_rho_multipole_spl(current_rho_multipole_spl, &
                                    current_spl_atom)
         atom_of_splines = current_spl_atom
      endif

      ! now loop over all points
      do i_point = 1, n_points, 1

         ! get distance to current center
         call tab_single_atom_centered_coords_p0 &
             ( current_center, &
               points(:,i_point), &
               dist_tab_sq, &
               dir_tab )
         
         ! Is the point within the multipole radius of the current atom?
         if ( dist_tab_sq .lt. multipole_radius_sq(current_spl_atom) ) then
            
            ! determine the maximum angular momentum required for the present atom 
            l_atom_max = l_hartree(species(current_spl_atom))
            do while ( (outer_potential_radius(l_atom_max, current_spl_atom) &
                        .lt. dist_tab_sq ) &
                        .and. (l_atom_max .gt. 0) )
               l_atom_max = l_atom_max - 1
            enddo

            l_h_dim = (l_atom_max + 1)**2
            
            call tab_single_atom_centered_coords_radial_log_p0 &
                ( current_center, dist_tab_sq, dir_tab,  &
                  dist_tab_in, i_r, i_r_log, dir_tab_in )
            
            call tab_single_radial_weights_v2 &
                ( current_spl_atom, dist_tab_in, i_r, &
                  log_weight, radial_weight )
              
            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)
              
            call tab_single_gradient_ylm_p2 &
                ( trigonom_tab, l_atom_max, l_pot_max,  &
                  ylm_tab, dylm_dtheta_tab, scaled_dylm_dphi_tab)
            
            ! spline delta rho multipole, obtain partitioned rho multipole
            call spline_vector_v2 &
                ( i_r+1,  current_rho_multipole_spl, &
                  l_p_dim, &
                  n_max_spline, &
                  n_max_radial+2,  &
                  n_radial(species_center(current_center))+2,  &
                  l_h_dim, &
                  rho_multipole_component)
            
            ! calculate delta density
            delta_rho_multipole(i_point) = delta_rho_multipole(i_point) + &
                  ddot( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)

            ! spline the gradient of delta rho multipole
            call spline_deriv_vector_v2 &
                ( i_r+1,  &
                  current_rho_multipole_spl,&                           
                  l_p_dim, &
                  n_max_spline, &
                  n_max_radial+2,  &
                  n_radial(species_center(current_center))+2,  &
                  l_h_dim, &
                  rho_multipole_deriv)
            
            ! multiply by some radial weight (?)
            rho_multipole_deriv(1:l_h_dim) = &
                  rho_multipole_deriv(1:l_h_dim) * radial_weight
            
            
            !! This would again be the free atoms part of the gradient
            !! BUT
            !! the free atoms part is already covered in a different routine
            !! so we don't need to compute its gradient here.
            !d_rho_free_d_r = pi4_inv * log_weight * val_spline_deriv( &
            !         i_r_log,  &
            !         renormalized_free_rho_spl(1,1,species(current_spl_atom)),  &
            !         n_grid(species(current_spl_atom)) &
            !         )
            !! side note: only approximately (?) equivalent to
            !d_rho_free_d_r = pi4_inv * val_spline &
            !    ( i_r_log, &
            !      renormalized_free_drho_dr_spl(1,1,species_center(current_center)), &
            !      n_grid(species_center(current_center)) )

            
            ! now calculate the gradient
            call evaluate_rho_multipole_gradient &
                ( dist_tab_in, dir_tab_in, &
                  trigonom_tab, &
                  ylm_tab, &
                  dylm_dtheta_tab, &
                  scaled_dylm_dphi_tab, &
                  rho_multipole_component, &
                  rho_multipole_deriv, &
                  0.0d0, & !d_rho_free_d_r
                  l_h_dim, &
                  rho_multipole_gradient_temp )

            ! add up gradient
            delta_rho_multipole_gradient(:,i_point) = &
                  delta_rho_multipole_gradient(:,i_point) + &
                  rho_multipole_gradient_temp(:)

         endif ! dist_tab_sq .lt. multipole_radius_sq

      enddo ! i_point
      
   enddo ! i_center


   ! DEALLOCATION
   
   if (allocated(ylm_tab)) &
      deallocate(ylm_tab)

   if (allocated(dylm_dtheta_tab)) &
      deallocate(dylm_dtheta_tab)

   if (allocated(scaled_dylm_dphi_tab)) &
      deallocate(scaled_dylm_dphi_tab)

   if (allocated(current_rho_multipole_spl)) &
      deallocate(current_rho_multipole_spl)

   if (allocated(rho_multipole_component)) &
      deallocate(rho_multipole_component)

   if (allocated(rho_multipole_deriv)) &
      deallocate(rho_multipole_deriv)


end subroutine get_delta_rho_multipole_and_gradient




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_free_rho
!  NAME
!    get_free_rho
!  SYNOPSIS
subroutine get_free_rho( &
      n_points, &
      points, &
      free_rho )
!  PURPOSE 
!    This subroutine evaluates the free electron density (density of the free 
!    atoms superposition) at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: free_rho

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where the density is evaluated
!  OUTPUT
!   o free_rho -- values of the free electron density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE


   ! distance and direction
   real*8, dimension(:), allocatable :: dist_tab_sq_free !(n_centers_basis_integrals)
   real*8, dimension(:), allocatable :: dist_tab_free !(n_centers_basis_integrals)
   real*8, dimension(:,:), allocatable :: dir_tab_free !(3,n_centers_basis_integrals)
   real*8, dimension(:,:), allocatable :: dir_tab_norm !(3,n_centers_basis_integrals)
   
   ! logarithmic grid shell indices
   real*8, dimension(:), allocatable :: i_r_free !(n_centers_basis_integrals)
   
   ! helpers
   integer :: current_center_free

   ! external function
   real*8, external :: ddot
   
   ! sparse representation of multipole center indices
   integer :: n_compute_atoms
   integer, dimension(:), allocatable :: center_index !n_centers_basis_integrals

   ! index variables
   integer :: i_point
   integer :: i_center_L, i_center_L2


   ! ALLOCATION
   
   allocate(dist_tab_sq_free(n_centers_basis_integrals))
   allocate(dist_tab_free(n_centers_basis_integrals))
   allocate(dir_tab_free(3,n_centers_basis_integrals))
   allocate(dir_tab_norm(3,n_centers_basis_integrals))
   allocate(i_r_free(n_centers_basis_integrals))
   allocate(center_index(n_centers_basis_integrals))


   free_rho = 0.d0

   ! evaluate free atoms superposition
   do i_point = 1, n_points, 1

      !Begin free atom part (as in initialize_grid_storage_p1)
      call tab_atom_centered_coords_p0 &
          ( points(:,i_point), &
            dist_tab_sq_free, &
            dir_tab_free, &
            n_centers_basis_integrals, centers_basis_integrals )
      
      ! relevant atoms
      n_compute_atoms = 0

      ! get relevant centers that are within the free atom's charge radius
      do i_center_L = 1, n_centers_basis_integrals, 1
      
         current_center_free = centers_basis_integrals(i_center_L)
         
         ! within radius?
         if ( (dist_tab_sq_free(i_center_L) .lt. &
               multipole_radius_free_sq(species_center(current_center_free)) ) .and.&
         ! not empty?
              (.not. empty(center_to_atom(current_center_free))) .and. &
         ! not pseudo-atom?
              (.not. species_pseudoized(species(center_to_atom(current_center_free)))) &
         ! then go
               ) then
            n_compute_atoms = n_compute_atoms + 1
            dist_tab_sq_free(n_compute_atoms) = dist_tab_sq_free(i_center_L) 
            dir_tab_free(:,n_compute_atoms) = dir_tab_free(:,i_center_L) 
            center_index(n_compute_atoms) = current_center_free
         endif
      enddo

      call tab_global_geometry_p0 &
          ( dist_tab_sq_free, &
            dir_tab_free, &
            dist_tab_free, &
            i_r_free, &
            dir_tab_norm, &
            n_compute_atoms, &
            center_index )
      
      ! calculate the free-atom density only for the (now) known atoms ...
      do i_center_L2 = 1, n_compute_atoms, 1

         current_center_free = center_index(i_center_L2)
         
         ! sum up density
         free_rho(i_point) = free_rho(i_point) + val_spline &
                ( i_r_free(i_center_L2), &
                  renormalized_free_rho_spl(1,1,species_center(current_center_free)), &
                  n_grid(species_center(current_center_free)) )
         
         ! As far as I understand, a call to evaluate_free_atom_sums_p2 is
         ! unnecessary. The free atom density has already
         ! been calculated from the interpolant.
      
      enddo

      ! divide out factor of pi4 to get free atom density and gradient
      free_rho(i_point) = free_rho(i_point) * pi4_inv
   
   enddo ! i_point


   ! DEALLOCATION
   
   if (allocated(dist_tab_sq_free)) &
      deallocate(dist_tab_sq_free)

   if (allocated(dist_tab_free)) &
      deallocate(dist_tab_free)

   if (allocated(dir_tab_free)) &
      deallocate(dir_tab_free)

   if (allocated(dir_tab_norm)) &
      deallocate(dir_tab_norm)

   if (allocated(i_r_free)) &
      deallocate(i_r_free)

   if (allocated(center_index)) &
      deallocate(center_index)


end subroutine get_free_rho




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/get_delta_rho_multipole
!  NAME
!    get_delta_rho_multipole
!  SYNOPSIS
subroutine get_delta_rho_multipole( &
      n_points, &
      points, &
      delta_rho_multipole )
!  PURPOSE 
!    This subroutine evaluates the multipole expanded delta electron density 
!    (i.e. the difference between the multipole expanded density and the free 
!    atoms superposition) at arbitrary points in real space.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

   real*8, dimension(n_points), intent(out) :: delta_rho_multipole

!  INPUTS
!   o n_points -- number of points at which the evaluation takes place
!   o points -- real space coordinates of points where the density is evaluated
!  OUTPUT
!   o delta_rho_multipole -- values of the multipole expanded delta density
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   
   ! distance and direction
   real*8 :: dist_tab_sq
   real*8 :: dist_tab_in
   real*8, dimension(3) :: dir_tab
   real*8, dimension(3) :: dir_tab_in
   
   ! trigonometric function values
   real*8, dimension(4) :: trigonom_tab
   real*8, dimension(:), allocatable :: ylm_tab !(l_pot_max+1)**2
   
   ! logarithmic grid shell indices
   real*8 :: i_r
   real*8 :: i_r_log
   
   ! interpolated and partitioned representation of delta density
   real*8, dimension(:,:,:), allocatable :: current_rho_multipole_spl !((l_pot_max+1)**2, n_max_spline, n_max_radial+2)
   real*8, dimension(:), allocatable :: rho_multipole_component !(l_pot_max+1)**2
   
   ! array dimensions
   integer :: l_p_dim ! (l_pot_max+1)**2
   integer :: l_atom_max
   integer :: l_h_dim ! (l_atom_max+1)**2, no constant

   ! helpers
   integer :: atom_of_splines
   integer :: current_spl_atom
   integer :: current_center

   ! external function
   real*8, external :: ddot
   
   ! index variables
   integer :: i_center
   integer :: i_point


   l_p_dim = (l_pot_max + 1)**2

   ! ALLOCATION
   
   allocate(ylm_tab(l_p_dim))
   allocate(current_rho_multipole_spl(l_p_dim, n_max_spline, n_max_radial+2))
   allocate(rho_multipole_component(l_p_dim))


   delta_rho_multipole = 0.d0

   ! this variable indicates for which center the density is currently
   ! splined. set to 0 in the beginning to force calculation in the first step.
   atom_of_splines = 0

   ! sum up contribution of multipole centers
   do i_center = 1, n_centers_hartree_multipole, 1

      ! translate center index
      current_center = centers_hartree_multipole(i_center)
      current_spl_atom = center_to_atom(current_center)

      ! check if the current atom has already been splined or
      ! if an update is necessary
      if (current_spl_atom .ne. atom_of_splines) then
         call get_rho_multipole_spl(current_rho_multipole_spl, &
                                    current_spl_atom)
         atom_of_splines = current_spl_atom
      endif

      ! now loop over all points
      do i_point = 1, n_points, 1

         ! get distance to current center
         call tab_single_atom_centered_coords_p0 &
             ( current_center, &
               points(:,i_point), &
               dist_tab_sq, &
               dir_tab )

         ! Is the point within multipole radius of the current atom?
         if ( dist_tab_sq .lt. multipole_radius_sq(current_spl_atom) ) then
            
            l_atom_max = l_hartree(species(current_spl_atom))
            do while ( (outer_potential_radius(l_atom_max, current_spl_atom) &
                        .lt. dist_tab_sq ) &
                        .and. (l_atom_max .gt. 0) )
               l_atom_max = l_atom_max - 1
            enddo

            l_h_dim = (l_atom_max + 1)**2
            
            call tab_single_atom_centered_coords_radial_log_p0 &
                ( current_center, dist_tab_sq, dir_tab,  &
                  dist_tab_in, i_r, i_r_log, dir_tab_in )
            
            call tab_single_trigonom_p0(dir_tab_in, trigonom_tab)
            
            call tab_single_wave_ylm_p2 &
                ( trigonom_tab, l_atom_max, l_pot_max,  &
                  ylm_tab)
            
            ! spline delta rho multipole, obtain partitioned rho multipole
            call spline_vector_v2 &
                ( i_r+1,  current_rho_multipole_spl, &
                  l_p_dim, &
                  n_max_spline, &
                  n_max_radial+2,  &
                  n_radial(species_center(current_center))+2,  &
                  l_h_dim, &
                  rho_multipole_component)
            
            ! calculate delta density
            delta_rho_multipole(i_point) = delta_rho_multipole(i_point) + &
                  ddot( l_h_dim, rho_multipole_component, 1, ylm_tab, 1)

         endif ! dist_tab_sq .lt. multipole_radius_sq

      enddo ! i_point
      
   enddo ! i_center


   ! DEALLOCATION
   
   if (allocated(ylm_tab)) &
      deallocate(ylm_tab)

   if (allocated(current_rho_multipole_spl)) &
      deallocate(current_rho_multipole_spl)

   if (allocated(rho_multipole_component)) &
      deallocate(rho_multipole_component)


end subroutine get_delta_rho_multipole




!******
!----------------------------------------------------------------------------------------------
!****s* rho_multipole_evaluation/test_rho_multipole_evaluation_routines
!  NAME
!    test_rho_multipole_evaluation_routines
!  SYNOPSIS
subroutine test_rho_multipole_evaluation_routines( &
      n_points, &
      points )
!  PURPOSE 
!    This routine calls all density evaluation routines in this module and
!    writes their output to disk.
!  USES
   implicit none
!  ARGUMENTS
   integer, intent(in) :: n_points
   real*8, dimension(3, n_points), intent(in) :: points

!  INPUTS
!   o n_points -- number of points at which the evaluations take place
!   o points -- real space coordinates of points where the density and
!               its gradient are evaluated
!  OUTPUT
!   writing to disk
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   
   ! temporary variables
   real*8, dimension(:), allocatable :: free_rho
   real*8, dimension(:,:), allocatable :: free_rho_gradient
   real*8, dimension(:), allocatable :: delta_rho_multipole
   real*8, dimension(:,:), allocatable :: delta_rho_multipole_gradient
   real*8, dimension(:), allocatable :: rho_multipole
   real*8, dimension(:,:), allocatable :: rho_multipole_gradient
   
   ! index variables
   integer :: i_point


   ! ALLOCATION
   
   allocate(free_rho(n_points))
   allocate(delta_rho_multipole(n_points))
   allocate(rho_multipole(n_points))
   allocate(free_rho_gradient(3,n_points))
   allocate(delta_rho_multipole_gradient(3,n_points))
   allocate(rho_multipole_gradient(3,n_points))


   ! OUTPUT points
   if (myid.eq.0) call write_debug1("points.dat", n_points, 3, points)

   ! TEST get_free_rho
   call get_free_rho(n_points, points, free_rho)
   if (myid.eq.0) &
      call write_debug1("get_free_rho.dat", n_points, 1, free_rho)

   ! TEST get_delta_rho_multipole
   call get_delta_rho_multipole(n_points, points, delta_rho_multipole)
   if (myid.eq.0) &
      call write_debug1("get_delta_rho_multipole.dat", n_points, &
                        1, delta_rho_multipole)

   ! TEST get_rho_multipole
   call get_rho_multipole(n_points, points, rho_multipole)
   if (myid.eq.0) &
      call write_debug1("get_rho_multipole.dat", n_points, 1, rho_multipole)

   ! TEST get_free_rho_and_gradient
   call get_free_rho_and_gradient(n_points, points, free_rho, free_rho_gradient)
   if (myid.eq.0) &
      call write_debug2("get_free_rho_and_gradient.dat", n_points, &
                        1, free_rho, 3, free_rho_gradient)

   ! TEST get_delta_rho_multipole_and_gradient
   call get_delta_rho_multipole_and_gradient(n_points, points, &
                           delta_rho_multipole, delta_rho_multipole_gradient)
   if (myid.eq.0) &
      call write_debug2("get_delta_rho_multipole_and_gradient.dat", n_points, &
                        1, delta_rho_multipole, 3, delta_rho_multipole_gradient)

   ! TEST get_rho_multipole_and_gradient
   call get_rho_multipole_and_gradient(n_points, points, &
                           rho_multipole, rho_multipole_gradient)
   if (myid.eq.0) &
      call write_debug2("get_rho_multipole_and_gradient.dat", n_points, &
                        1, rho_multipole, 3, rho_multipole_gradient)


   ! DEALLOCATION

   if (allocated(free_rho)) &
      deallocate(free_rho)

   if (allocated(delta_rho_multipole)) &
      deallocate(delta_rho_multipole)

   if (allocated(rho_multipole)) &
      deallocate(rho_multipole)

   if (allocated(free_rho_gradient)) &
      deallocate(free_rho_gradient)

   if (allocated(delta_rho_multipole_gradient)) &
      deallocate(delta_rho_multipole_gradient)

   if (allocated(rho_multipole_gradient)) &
      deallocate(rho_multipole_gradient)


   contains
   
   
   subroutine write_debug1(filename, n_points, dim1, val1)
   
      implicit none
      
      character(*), intent(in) :: filename
      integer, intent(in) :: n_points
      integer, intent(in) :: dim1
      real*8, dimension(dim1, n_points), intent(in) :: val1
      
      integer, parameter :: myunit = 36
      
      integer :: i_point
      character*64 :: fmtstring
      
      ! write format string
      write(fmtstring,'(A,I1,A)') '(I4,1X,', dim1,'(1X,E18.8))'
      
      ! open file
      open(unit=myunit, file=trim(filename), status='unknown', action='write')
      
      ! write debug output
      do i_point = 1, n_points, 1
         write(myunit,fmtstring) i_point, val1(:,i_point)
      enddo ! i_point
      
      ! close file
      close(unit=myunit)
      
   end subroutine write_debug1


   subroutine write_debug2(filename, n_points, dim1, val1, dim2, val2)
   
      implicit none
      
      character(*), intent(in) :: filename
      integer, intent(in) :: n_points
      integer, intent(in) :: dim1
      real*8, dimension(dim1, n_points), intent(in) :: val1
      integer, intent(in) :: dim2
      real*8, dimension(dim2, n_points), intent(in) :: val2
      
      integer, parameter :: myunit = 36
      
      integer :: i_point
      character*64 :: fmtstring
      
      ! write format string
      write(fmtstring,'(A,I1,A)') '(I4,1X,', dim1+dim2,'(1X,E18.8))'
      
      ! open file
      open(unit=myunit, file=trim(filename), status='unknown', action='write')
      
      ! write debug output
      do i_point = 1, n_points, 1
         write(myunit,fmtstring) i_point, val1(:,i_point), val2(:,i_point)
      enddo ! i_point
      
      ! close file
      close(unit=myunit)
      
   end subroutine write_debug2
   
end subroutine test_rho_multipole_evaluation_routines


end module rho_multipole_evaluation

