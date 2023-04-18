!****h* FHI-aims/isc_common
!  NAME
!    isc_common
!  SYNOPSIS

module isc_common

!  PURPOSE
!    This module contains subroutines, functions and types used in multiple
!    isc modules.
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
!              (this module)
!       Contains types and functions
!       used in multiple isc modules
!
!
!    5. mpe_dielectric_interfaces_common
!
!       Contains subroutines used in
!       multiple higher-level modules
!
!    To avoid circular dependencies modules exclusively use subroutines and
!    variables from lower level modules
!
!    Short description of public subroutines and functions (please refer to
!    headers for more specific details):
!
!     o update_neighbor_list_intra
!       TASKS
!       - Updates the neighbor list of a cavity point, based on distances
!         to all other points
!       CALLS
!         during constraint_dynamics and update_neighbor_lists_for_Voronoi
!
!     o distribute_loads
!       TASKS
!       - distributes cavity points between different cpus
!       CALL
!         during constraint_dynamics_parallel and get_new_grid
!
!     o scale_interaction_distance
!       TASKS
!       - scales interaction distance between cavity points (for
!         constraint dynamics) based on local curvature
!       CALL
!         during update_neighbor_list_intra and get_intra_distances
!
!  USES

   use types, only: dp
   use localorb_io, only: localorb_info
   use mpi_tasks, only: n_tasks
   use mpe_types, only: &
         InterfacePoint
   use mpe_dielectric_interfaces_common, only: &
         dot_d3

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
   public update_neighbor_list_intra
   public distribute_loads
   public scale_interaction_distance
   
   ! PUBLIC TYPES
   public   Neighbor, NeighborList

   real(dp), public :: charlensq ! (square of) characteristic length, important to account
                       ! for different system sizes and number of beads

   real(dp), public :: distsq_threshold  ! distance threshold between beads


   ! TYPE DEFINITIONS
   type Neighbor
      integer :: ind
      real(dp) :: relvec(3)
      real(dp) :: distsq
   end type

   type NeighborList
      type(Neighbor), allocatable :: nb(:)
   end type

contains

!******
!-------------------------------------------------------------------------------
!****s* isc_common/update_neighbor_list_intra
!  NAME
!    update_neighbor_list_intra
!  SYNOPSIS

subroutine update_neighbor_list_intra(points, i_point, &
      distsq_neighbor_threshold, nlist)

!  PURPOSE
!    This subroutine calculates the cartesian distance vectors and
!    squared cartesian distances between a certain point with index
!    i_point and all other points;
!    from this information a neighbor list is created
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(in) :: points(:)
   integer, intent(in) :: i_point
   real(dp), intent(in) :: distsq_neighbor_threshold

   type(NeighborList), intent(out) :: nlist

!  INPUTS
!   o points -- sampling points with assigned coordinate systems
!   o i_point -- index of point for which neighbor list is created
!   o distsq_neighbor_threshold -- threshold for squared cartesian distance
!                                  to consider two points neighboring
!  OUTPUT
!   o nlist -- index list of all neighboring points and the
!              corresponding distance information
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(132) :: info_str
   integer :: n_points
   type(Neighbor), allocatable :: tmp(:)
   integer :: i_nb, ncount

   n_points = size(points)

   !ALLOCATION
   allocate(tmp(n_points))

   ! reset neighbor count
   ncount = 1

   ! first loop: calculate distances and determine number of neighbors
   do i_nb = 1, n_points, 1

      ! distance calculation
      tmp(ncount)%relvec = points(i_point)%coord - points(i_nb)%coord
      tmp(ncount)%distsq = scale_interaction_distance( &
                     dot_d3(tmp(ncount)%relvec, tmp(ncount)%relvec), &
                     points(i_point)%normal, &
                     points(i_nb)%normal )

      ! when neighboring, keep in temporary array
      ! exclude the point itself
      if ( (tmp(ncount)%distsq.le.distsq_neighbor_threshold) .and. &
           (i_nb .ne. i_point) ) then
         tmp(ncount)%ind = i_nb
         ncount = ncount + 1
      endif

   enddo ! i_nb

   ! deallocate old neighbor list information
   if (allocated(nlist%nb)) &
      deallocate(nlist%nb)

   ! the counter of neighbors is 1 ahead, so
   ncount = ncount - 1

   ! check if neighbors have been found
   if ( ncount.eq.0 ) then
      write(info_str,'(2X,A)') '*** WARNING!'
      call localorb_info(info_str)
      write(info_str,'(2X,A,X,I6)') 'No neighboring bead found for bead no.', &
                                       i_point
      call localorb_info(info_str)
      write(info_str,'(2X,A)') 'Please choose larger distance cutoff.'
      call localorb_info(info_str)
   else
      ! allocate new neighbor list array
      allocate(nlist%nb(ncount))

      ! second loop: copy from temporary array
      do i_nb = 1, ncount, 1
         nlist%nb(i_nb) = tmp(i_nb)
      enddo ! i_nb
   endif

   ! DEALLOCATION
   deallocate(tmp)

end subroutine update_neighbor_list_intra


!******
!-------------------------------------------------------------------------------
!****s* isc_common/distribute_loads
!  NAME
!    distribute_loads
!  SYNOPSIS

subroutine distribute_loads(total_load, n_workers, load, offset)

!  PURPOSE
!    This subroutine calculates the distribution of a number of items
!    among a certain number of tasks in continuous chunks.
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: total_load
   integer, intent(in) :: n_workers

   integer, dimension(0:n_workers-1), intent(out) :: load
   integer, dimension(0:n_workers-1), intent(out) :: offset

!  INPUTS
!   o total_load -- number of items to be shared
!   o n_workers -- number of tasks that take a chunk
!  OUTPUT
!   o load -- number of items for each worker
!   o offset -- (index of first item -1) for each worker
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: remaining
   integer :: ind_task
   integer, parameter :: first_offset = 0

   ! first estimate
   load = total_load / n_workers

   ! now, some load might be remaining
   remaining = total_load - load(0)*n_workers

   ! distribute remaining items
   do ind_task = 0, remaining-1, 1
      load(ind_task) = load(ind_task) + 1
   enddo ! ind_task

   ! first task's offset
   offset(0) = first_offset

   ! calculate all other offsets
   do ind_task = 1, n_tasks-1, 1
      offset(ind_task) = offset(ind_task-1) + load(ind_task-1)
   enddo ! ind_task

end subroutine distribute_loads

!******
!-------------------------------------------------------------------------------
!****s* isc_common/scale_interaction_distance
!  NAME
!    scale_interaction_distance
!  SYNOPSIS

pure function scale_interaction_distance(distsq, n1, n2)

!  PURPOSE
!    This function calculates a scaling factor for the interaction distance
!    between two points; this correction is motivated by the following ideas:
!
!    the cavity will have different local curvature;
!    - at areas of large local curvature, the projection of the cartesian
!      distance on the normal vectors gets very small which leads to a
!      gathering of points in these areas
!    - simply re-directing the distance vector perpendicular to the normal
!      vector leads to an overestimation of the force as the cartesian
!      distance is smaller than the actual distance on the curved surface
!
!    Therefore, we re-direct the distance vector perpendicular to the normal
!    vector and approximately correct the distance by the factor calculated
!    by this function
!
!  USES
   implicit none

!  ARGUMENTS
   real(dp), intent(in) :: distsq
   real(dp), dimension(3), intent(in) :: n1
   real(dp), dimension(3), intent(in) :: n2

   real(dp) scale_interaction_distance

!  INPUTS
!   o distsq -- distance of the two points 1 and 2, squared
!   o n1 -- normal vector on point 1
!   o n2 -- normal vector on point 2
!  OUTPUT
!   o scale_interaction_distance -- scaling factor for the interaction distance
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! distance correction order
   integer, parameter :: isc_distance_correction_order = 2

   real(dp) :: one_minus_cosphi
   real(dp) :: scaling

   ! ratio between direct distance d on a secant vs distance on the
   ! (unit) circle c, where phi is the angle between the normal vectors
   ! (n1, n2) on the circle at the secant intersections (j, k):
   !
   !  c/d = (phi/2) / sin(phi/2)
   !
   ! thus
   !  c^2/d^2 = (phi/2)^2 / sin^2(phi/2)
   !
   ! alternate form (done by wolframalpha):
   !  c^2/d^2 = phi^2 / 2 / ( 1 - cos(phi) )
   !
   ! with
   !  cos(phi) = dot_product(n1, n2) / norm(n1) / norm(n2)
   !  phi^2 = acos^2( cos(phi) )
   !
   ! Taylor expansion of acos^2(x) at x=1 (done by wolframalpha):
   !
   !  acos^2(x) = 2(1-x) + 1/3*(1-x)^2 + 4/45*(1-x)^3 + O((1-x)^4)
   !
   ! and thus
   !
   !  c^2/d^2 = 1 + 1/6*(1-cos(phi)) + 2/45*(1-cos(phi))^2
   !

   ! normal vectors should be normalized already
   one_minus_cosphi = 1.e0_dp - dot_d3(n1, n2)

   if (isc_distance_correction_order.eq.2) then
      ! use third order expansion of acos^2
      scaling = 1.e0_dp + one_minus_cosphi/6.e0_dp + &
                  one_minus_cosphi*one_minus_cosphi/22.5e0_dp
   elseif (isc_distance_correction_order.eq.1) then
      ! use second order expansion of acos^2
      scaling = 1.e0_dp + one_minus_cosphi/6.e0_dp
   else
      ! no scaling (use first order expansion of acos^2)
      scaling = 1.e0_dp
   endif

   scale_interaction_distance = scaling * distsq
   !
   ! However, this does not seem to have much influence !

end function scale_interaction_distance

end module isc_common

