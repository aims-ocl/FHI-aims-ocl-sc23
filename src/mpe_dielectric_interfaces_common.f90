!****h* FHI-aims/mpe_dielectric_interfaces_common
!  NAME
!    mpe_dielectric_interfaces_common
!  SYNOPSIS

module mpe_dielectric_interfaces_common

!  PURPOSE
!    This module contains subroutines and functions of general use for
!    dielectric interfaces in the MPE implicit solvation model. Subroutines for
!    specific types of interfaces are contained in respective modules.
!    Subroutines to be called from the outside world are contained here and in
!    mpe_dielectric_interfaces.f90
!
!    Schematically, the hierarchy of dielectric interface modules goes as
!    follows:
!
!    1.                      mpe_dielectric_interfaces
!
!       Contains subroutines to be called "from the outside world", currently
!       in mpe_interface
!                           /                        \
!                          /                          \
!    2. isc_implicit_solvent_cavity          ifp_dielectric_interface_plane
!
!       Contains subroutines for the         Contains subroutines for the
!       isodensity cavity around the         (optional) interface between
!       solute                       \__     two solvents
!          |                            \__
!          |                               \
!    3. isc_projected_voronoi                isc_constraint_dynamics
!
!       Contains relatively self-con-        Contains the MD-like cavity
!       tained subroutines used in           optimization procedure
!       isc_implicit_solvent_cavity      __/
!          |                          __/
!          |                         /
!    4. isc_common
!
!       Contains types and functions
!       used in multiple isc modules
!
!
!    5. mpe_dielectric_interfaces_common
!              (this module)
!       Contains subroutines used in
!       multiple higher-level modules
!
!    To avoid circular dependencies modules exclusively use subroutines and
!    variables from lower level modules
!
!    Short description of public subroutines and functions (please refer to
!    headers for more specific details):
!
!     o get_perpendicular_d3
!       TASKS
!       - returns a vector which is perpendicular to one or two input vectors
!
!     o cross_d3
!       TASKS
!       - returns the cross product of two input vectors
!
!     o normalized_d3
!       TASKS
!       - returns normalized vector in same direction as input vector
!
!     o dot_d3
!       TASKS
!       - returns the dot product of two input vectors
!
!     o iswap
!       TASKS
!       - swaps the values of two input integers
!
!     o calculate_distance_to_plane
!       TASKS
!       - calculates distance of a point to a plane
!
!     o determine_species_radius
!       TASKS
!       - assigns an atomic radius to a species, either by user choice or by
!         inversion of splined atomic density and choosing radius with
!         isodensity value
!
!     o merge_interfaces
!       TASKS
!       - gathers interface points of two interfaces in the first
!
!     o get_atoms_and_radii
!       TASKS
!       - allocates and fills the arrays grid_atoms and R_atom
!
!  USES

   use localorb_io, only: localorb_info
   use constants, only: pi4, bohr
   use types, only: dp
   use free_atoms, only: renormalized_free_rho_spl
   use grids, only: n_grid, r_grid
   use spline, only: cubic_spline, val_spline
   use runtime_choices, only: &
         isc_initial_radius, &
         isc_isodensity_value, &
         mpe_lmax_ep_additional_order
   use dimensions, only: n_atoms, n_species
   use geometry, only: coords, species, empty
   use species_data, only: species_name, species_pseudoized
   use mpe_types, only: &
         DielectricInterface, &
         InterfacePoint
   use species_data, only: l_shell_max

   implicit none
   
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications 180 (2009), 2175-2196.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! make everything private by default
   private

   ! PUBLIC subroutines and functions
   public get_perpendicular_d3
   public cross_d3
   public normalized_d3
   public dot_d3
   public iswap
   public calculate_distance_to_plane
   public determine_species_radius
   public merge_interfaces
   public get_atoms_and_radii

   real(dp), public :: center(3) ! center of atom coordinates
   
   real(dp), allocatable, public :: grid_atoms(:,:), R_atom(:) ! coordinates and radii of atoms
   
   real(dp), public :: charlen_for_if_gen
   
   contains

   
!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/get_perpendicular_d3
!  NAME
!    get_perpendicular_d3
!  SYNOPSIS

pure function get_perpendicular_d3(a, specified_direction) result(p)

!  PURPOSE
!    This functions returns a perpendicular vector to the 3-dimensional
!    input vector. If possible, the vector is chosen such that it is also
!    perpendicular to specified_direction.
!
!  USES
   implicit none

   real(dp), dimension(3) :: p
   real(dp), dimension(3), intent(in) :: a
   real(dp), dimension(3), intent(in), optional :: specified_direction
!  INPUTS
!   o a -- three dimensional vector
!   o b -- three dimensional vector
!  RESULT
!   o p -- cross-product of {a} and {b}
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   real(dp), parameter :: eps = 1.e-15_dp, default_pd(3) = (/1.e0_dp,0.e0_dp,0.e0_dp/)
   real(dp) :: pd(3), normsq

   pd = default_pd
   if (present(specified_direction)) then
      pd = specified_direction
   endif

   p = cross_d3(a, pd)
   normsq = dot_d3(p,p)

   if (normsq .gt. eps*eps) then
      ! p is perpendicular to a and specified_direction, normalize
      p = p / sqrt(normsq)
   else
      ! a and specified_direction are collinear, workaround:
      ! determine non-collinear vector wrt a
      if ( abs(a(1)) .lt. eps ) then
         ! vector in x-direction is perpendicular
         pd = (/ 1.e0_dp, 0.e0_dp, 0.e0_dp /)
      elseif ( abs(a(2)) .lt. eps ) then
         ! vector in y-direction is perpendicular
         pd = (/ 0.e0_dp, 1.e0_dp, 0.e0_dp /)
      else
         ! vector "p*(-1, 1, 1)" is at least not collinear
         ! under the preconditions made before
         pd(1) = -a(1)
         pd(2) =  a(2)
         pd(3) =  a(3)
      end if

      p = normalized_d3(cross_d3(a, pd))
   endif

end function get_perpendicular_d3



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/cross_d3
!  NAME
!    cross_d3
!  SYNOPSIS

pure function cross_d3(a, b) result(c)

!  PURPOSE
!    This functions returns the cross-product of the three-dimensional
!    input vectors.
!
!  USES
   implicit none

   real(dp), intent(in) :: a(3), b(3)
   real(dp) :: c(3)
!  INPUTS
!   o a -- three dimensional vector
!   o b -- three dimensional vector
!  OUTPUT
!   o c -- cross-product of {a} and {b}
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   c(1) = a(2)*b(3) - a(3)*b(2)
   c(2) = a(3)*b(1) - a(1)*b(3)
   c(3) = a(1)*b(2) - a(2)*b(1)

end function cross_d3



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/normalized_d3
!  NAME
!    normalized_d3
!  SYNOPSIS

pure function normalized_d3(a, threshold) result(n)

!  PURPOSE
!    This functions returns the normalized three-dimensional input vector.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: a(3)
   real(dp), intent(in), optional :: threshold
   real(dp) :: n(3)

!  INPUTS
!   o a -- three dimensional vector
!  RESULT
!   o n -- normalized vector
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
   real(dp), parameter :: default_threshold = 1.e-10_dp
   real(dp) :: threshsq, normsq

   threshsq = default_threshold*default_threshold
   if (present(threshold)) threshsq = threshold*threshold

   normsq = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
   if (normsq.lt.threshsq) then
      n = a
   else
      n = a / sqrt(normsq)
   endif

end function normalized_d3



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/dot_d3
!  NAME
!    dot_d3
!  SYNOPSIS

pure function dot_d3(a, b)

!  PURPOSE
!    This functions returns the dot-product of the three-dimensional
!    input vectors.
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), dimension(3), intent(in) :: a
   real(dp), dimension(3), intent(in) :: b
   real(dp) :: dot_d3

!  INPUTS
!   o a -- three dimensional vector
!   o b -- three dimensional vector
!  OUTPUT
!   o dot_d3 -- dot-product of {a} and {b}
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   dot_d3 = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

end function dot_d3


!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/iswap
!  NAME
!    iswap
!  SYNOPSIS

subroutine iswap(a, b)

!  PURPOSE
!    This functions swaps the value of two integer variables.
!
!  USES
   implicit none

   integer, intent(inout) :: a
   integer, intent(inout) :: b

!  INPUTS
!   o a -- integer with value a
!   o b -- integer with value b
!  OUTPUT
!   o a -- integer with value b
!   o b -- integer with value a
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: itemp

   itemp = a
   a = b
   b = itemp

end subroutine iswap


!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/calculate_distance_to_plane
!  NAME
!    calculate_distance_to_plane
!  SYNOPSIS

pure function calculate_distance_to_plane(point, plane_normal, plane_dist) &
   result(distance)

!  PURPOSE
!    This functions calculates the distance of a point to a plane given in
!    hesse form
!
!  USES
   implicit none
   real(dp), intent(in) :: point(3), plane_normal(3), plane_dist

!  INPUTS
!   o point -- coordinates of point
!   o plane_normal -- normal vector of plane
!   o plane_dist -- distance of plane to origin
!  OUTPUT
!   o distance -- distance of point to plane
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   real(dp) :: distance

   distance = dot_d3(point, plane_normal) - plane_dist
end function calculate_distance_to_plane



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/determine_species_radius
!  NAME
!    determine_species_radius
!  SYNOPSIS

function determine_species_radius(i_species, rho_iso) result(radius)

!  PURPOSE
!    Determine the radii of the initial superposition of spheres
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: i_species
   real(dp), intent(in) :: rho_iso

   real(dp) :: radius
!  INPUTS
!   o i_spec -- species index
!   o rho_iso -- desired iso-density value
!  RETURNS
!   o radius -- radius of sphere around this species
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   if (isc_initial_radius.gt.0) then
      radius = isc_initial_radius
   else
      ! attention: internally the free atoms density and its spline
      !            are scaled by a factor of 4*Pi, thus also multiply
      !            the desired iso-density value by that
      radius = get_free_atom_radius_from_density_spline( &
                  rho_iso*pi4, &
                  n_grid(i_species), &
                  r_grid(1,i_species), &
                  renormalized_free_rho_spl(1,1,i_species))
   endif

   contains

   function get_free_atom_radius_from_density_spline(rho_iso, &
                     n_grid, r_grid, spline) &
      result(radius)

      !use spline, only: cubic_spline, val_spline
      implicit none

      real(dp), intent(in) :: rho_iso
      integer, intent(in) :: n_grid
      real(dp), dimension(n_grid), intent(in) :: r_grid
      real(dp), dimension(4,n_grid), intent(in) :: spline

      real(dp) :: radius


      integer :: i_grid, n_inv
      real(dp), dimension(:,:), allocatable :: inv_spl_param
      real(dp) :: rho_frac

      ! find invertible region (exclude all zeros except for first one)
      do i_grid = n_grid, 1, -1
         if (spline(1,i_grid).gt.0.e0_dp) exit
      enddo
      n_inv = min(n_grid, i_grid + 1)

      ! get matching window for fractional coordinate for spline
      do i_grid = n_inv, 1, -1
         if (spline(1,i_grid).ge.rho_iso) exit
      enddo

      ! catch corner cases
      if (i_grid.eq.n_inv) then
         ! density was probably zero, so yield right border
         radius = r_grid(i_grid)
      elseif (i_grid.eq.0) then
         ! density on left border is still too small,
         ! but we just yield the left border
         radius = r_grid(1)
      else
         ! spline inverse function
         allocate(inv_spl_param(4,n_inv))
         call cubic_spline(r_grid(1), n_inv, inv_spl_param)

         rho_frac = real(i_grid,dp) + &
                     (spline(1,i_grid) - rho_iso) / &
                     (spline(1,i_grid) - spline(1,i_grid+1))

         ! get inverse
         radius = val_spline(rho_frac, inv_spl_param, n_inv)
         deallocate(inv_spl_param)
      endif

   end function get_free_atom_radius_from_density_spline

end function determine_species_radius



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/merge_interfaces
!  NAME
!    merge_interfaces
!  SYNOPSIS

subroutine merge_interfaces(interface1, interface2)

!  PURPOSE
!    Gathers interface points of two interfaces in the first of the two
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(inout) :: interface1, interface2

!  INPUTS
!   o interface1 -- first dielectric interface
!   o interface2 -- second dielectric interface
!  OUTPUT
!   o interface1 -- first dielectric interface, now containing all points
!   o interface2 -- second dielectric interface, now empty
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'merge_interfaces'

   integer :: n_p(2), n_total
   type(InterfacePoint), allocatable :: tmp(:)

   ! get sizes
   n_p(1) = size(interface1%p)
   n_p(2) = size(interface2%p)
   n_total = sum(n_p)
   ! increase size of first interface to hold all points
   allocate(tmp(n_p(1)))
   tmp = interface1%p
   deallocate(interface1%p)
   allocate(interface1%p(n_total))
   interface1%p(1:n_p(1)) = tmp
   deallocate(tmp)
   ! fill with remaining points of other interface(s)
   interface1%p(1+n_p(1):n_p(2)+n_p(1)) = interface2%p
   deallocate(interface2%p)
   allocate(interface2%p(0))

end subroutine merge_interfaces





!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces_common/get_atoms_and_radii
!  NAME
!    get_atoms_and_radii
!  SYNOPSIS

subroutine get_atoms_and_radii(lmax_center)

!  PURPOSE
!    This subroutine allocates and fills the arrays grid_atoms and R_atom
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(out), allocatable, optional :: lmax_center(:)

!  INPUT
!   none
!  OUTPUT
!   none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   real(dp), allocatable :: R_spec(:)

   integer :: i_species, j_atom, i_grid_atom, n_grid_atoms

   character(132) :: info_str

   integer, allocatable :: lmax_species(:)

   allocate(R_spec(n_species))
   if (present(lmax_center)) allocate(lmax_species(n_species))
   ! get radii of all atoms
   if (isc_initial_radius.gt.0) then
      write(info_str,'(6X,A,A)') 'Use specified radius for ', &
            'the initial spheres'
      call localorb_info(info_str)
   else
      write(info_str,'(6X,A,A)') 'Calculate initial radii from ', &
            'the inversion of the splined free atoms density'
      call localorb_info(info_str)
   endif

   if (present(lmax_center)) then
      write(info_str,'(8X,A7,2X,A4,2X,A11,2X,A4)') 'species', 'name', 'radius [AA]', 'lmax'
   else
      write(info_str,'(8X,A7,2X,A4,2X,A11)') 'species', 'name', 'radius [AA]'
   endif
   call localorb_info(info_str)
   do i_species = 1, n_species, 1
      if (present(lmax_center)) then
         R_spec(i_species) = determine_species_radius(i_species, &
                                 isc_isodensity_value)
         lmax_species(i_species) = l_shell_max(i_species)

         write(info_str,'(8X,I7,2X,A4,2X,F11.7,2X,I4)') i_species, &
               species_name(i_species), R_spec(i_species)*bohr, &
               lmax_species(i_species)
         call localorb_info(info_str)
      else
         R_spec(i_species) = determine_species_radius(i_species, &
                                 isc_isodensity_value)

         write(info_str,'(8X,I7,2X,A4,2X,F11.7)') i_species, &
               species_name(i_species), R_spec(i_species)*bohr
         call localorb_info(info_str)
      endif
   enddo ! i_species


   ! define centers for grids and create a list of radii for all atoms
   ! here, we just take "real" atomic centers and not pseudo-potentials
   n_grid_atoms = 0
   do j_atom = 1, n_atoms
      if (.not. (empty(j_atom).or.species_pseudoized(species(j_atom))) ) then
         n_grid_atoms = n_grid_atoms + 1
      endif
   enddo ! j_atom
   allocate(grid_atoms(3,n_grid_atoms))
   allocate(R_atom(n_grid_atoms))
   if (present(lmax_center)) allocate(lmax_center(n_grid_atoms))
   i_grid_atom = 0
   do j_atom = 1, n_atoms
      if (.not. (empty(j_atom).or.species_pseudoized(species(j_atom))) ) then
         i_grid_atom = i_grid_atom + 1
         grid_atoms(:,i_grid_atom) = coords(:,j_atom)
         R_atom(i_grid_atom) = R_spec(species(j_atom))
         if (present(lmax_center)) lmax_center(i_grid_atom) = &
                        lmax_species(species(j_atom)) + mpe_lmax_ep_additional_order
      endif
   enddo ! j_atom

   deallocate(R_spec)
   if (present(lmax_center)) deallocate(lmax_species)


end subroutine get_atoms_and_radii





end module mpe_dielectric_interfaces_common
