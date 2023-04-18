!****h* FHI-aims/isc_implicit_solvent_cavity
!  NAME
!    isc_implicit_solvent_cavity
!  SYNOPSIS

module isc_implicit_solvent_cavity

!  PURPOSE
!    This module contains routines to define an implicit solvent cavity (isc)
!    around the given structure in the QM calulation which is characterized
!    by a certain iso-density value {isc_isodensity_value};
!    The module is optimized for use in the module mpe_reaction_field
!
!    Subroutines shared with ifp_dielectric_interface_plane are contained in
!    mpe_dielectric_interfaces_common.f90
!    Subroutines to be called from the outside world are contained in
!    mpe_dielectric_interfaces.f90 and
!    mpe_dielectric_interfaces_common.f90
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
!              (this module)
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
!
!       Contains subroutines used in
!       multiple higher-level modules
!
!    To avoid circular dependencies modules exclusively use subroutines and
!    variables from lower level modules
!
!    Short description of public subroutines (please refer to subroutine
!    headers for more specific details):
!
!     o initialize_isc_keyword_defaults
!       TASKS
!       - initialize defaults for some keywords if not set by user
!       CALL
!         once during the initialization of the implicit solvent model
!
!     o get_new_grid
!       TASKS
!       - initialization of a new iso-density cavity from a spherical
!         Lebedev-Laikov grid
!       - allocation of public module variables
!       CALL
!         once during the initialization of the implicit solvent model
!
!     o update_implicit_solvent_cavity
!       TASKS
!       - optimization of cavity to current density value
!       - calculation of cavity's surface and volume for non-electrostatic
!         contributions
!       - redistribution of cavity points into subsets in different solvents
!       CALL
!         once if cavity is kept static, otherwise every step
!         in the scf cycle.
!         (provides necessary information for the calculation of the
!         reaction field)
!
!    Short description of public module variables:
!
!     o isc_surface_area, isc_cavity_volume
!       (approximate) surface area and volume of the cavity
!
!  USES

   use types, only: dp
   use constants, only: pi, bohr
   use localorb_io, only: localorb_info
   use timing, only: get_timestamps, get_times, output_times
   use mpi_tasks, only: myid, n_tasks, aims_stop, aims_stop_coll, &
                        mpi_comm_global
   use synchronize_mpi_basic, only: sync_vector_integer
   use quicksort_index, only: dquicksort_indexlist
   use mpe_constants, only: ISC_CONST
   use mpe_types, only: &
         DielectricInterface, &
         InterfacePoint, &
         InterfacePoint_vector_mpi_bcast
   use isc_projected_voronoi, only: &
         get_Voronoi_surface, get_volume_from_surface
   use runtime_choices, only: &
         isc_cavity_type, &
         isc_isodensity_value, &
         isc_calculate_surface_and_volume, &
         ifp, &
         ifp_dist, &
         ifp_normal, &
         isc_surface_curvature_correction, &
         isc_surface_sampling
   use mpe_dielectric_interfaces_common, only: &
         center, &
         R_atom, &
         grid_atoms, &
         charlen_for_if_gen, &
         get_perpendicular_d3, &
         cross_d3, &
         normalized_d3, &
         dot_d3, &
         iswap, &
         calculate_distance_to_plane, &
         merge_interfaces
   use isc_constraint_dynamics, only: constraint_dynamics, &
         initialize_cd_keyword_defaults
   use isc_common, only: &
         update_neighbor_list_intra, &
         distribute_loads, &
         distsq_threshold, &
         charlensq, &
         NeighborList


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
   public   initialize_isc_keyword_defaults
   public   get_new_grid
   public   update_implicit_solvent_cavity

   ! PUBLIC variables

   real(dp), public :: isc_surface_area, isc_cavity_volume, plane_origin(3)



   ! PRIVATE variables that are utilized in the creation and
   ! update process of the cavity, thus global

   !TODO make parameters adjustable
   real(dp), parameter :: exclude_points_near_overlap_of_spheres = 1.e-2_dp ! TODO: what does this do?
   real(dp), parameter :: exclude_points_on_interface_delta = 1.e-2_dp

contains


!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/distribute_cavity_interface
!  NAME
!    distribute_cavity_interface
!  SYNOPSIS

subroutine distribute_cavity_interface(interfaces, minimal_distance)

!  PURPOSE
!    This subroutine distributes the before gathered cavity sampling points
!    among the dielectric interfaces. Points too close to the interface
!    are discarded.
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(inout) :: interfaces(:)
   real(dp), intent(in) :: minimal_distance

!  INPUTS
!   o interfaces -- dielectric interfaces, gathered cavity
!   o minimal_distance -- all points closer to other interfaces are discarded
!  OUTPUT
!   o interfaces -- dielectric interfaces, distributed cavity
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'distribute_cavity_interface'
   character(132) :: info_str

   type(InterfacePoint), pointer :: pos_dist(:), neg_dist(:)
   type(InterfacePoint), target, allocatable :: tmp(:)
   integer :: n_interfaces

   n_interfaces = size(interfaces)
   select case (n_interfaces)
      case (1)
         continue ! everything's in place already

      case (4)
         ! prepare copy
         allocate(tmp(size(interfaces(1)%p)))
         tmp = interfaces(1)%p
         ! divide copy into parts
         call divide_points_with_plane(tmp, ifp_normal, ifp_dist, &
                                       minimal_distance, neg_dist, pos_dist )
         ! copy parts into respective interfaces
         deallocate(interfaces(1)%p)
         deallocate(interfaces(2)%p)
         allocate(interfaces(1)%p(size(pos_dist)))
         allocate(interfaces(2)%p(size(neg_dist)))
         interfaces(1)%p = pos_dist
         interfaces(2)%p = neg_dist
         ! clean-up
         neg_dist => null()
         pos_dist => null()
         deallocate(tmp)

      case default
         write(info_str,'(A)') '***Internal error: unexpected number of '//&
                            'interfaces'
         call localorb_info(info_str)
         call aims_stop('Internal error: unexpected number of '//&
                        'interfaces', func)
   end select

end subroutine distribute_cavity_interface



!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/divide_points_with_plane
!  NAME
!    divide_points_with_plane
!  SYNOPSIS

subroutine divide_points_with_plane(points, plane_normal, plane_dist, &
               distance_threshold, neg, pos)

!  PURPOSE
!    This subroutine sorts the array points such that all points on the
!    one side of the plane defined by its normal and distance from origin
!    end up as continuous block and the ones on the other side as well.
!    Points that are too close to the plane (distance < threshold) are
!    sorted to the ends of the array.
!
!    Method:
!    approach the array from both ends, negative distances should end up
!    in the left part and positive distances in the right part in the end.
!
!     A. calculate distances from the left until one is negative (l)
!        (close points are swapped with the first nonzero ones on-the-fly)
!     B. calculate distances from the right until one is positive (r)
!        (close points are swapped with the last nonzero ones on-the-fly)
!     C. if l<r, swap items and repeat from A
!     D. done, we have now something like (0+++++------00)
!        where n_neg and n_pos are pointers to the corresponding parts

!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), target, intent(inout) :: points(:)
   real(dp), intent(in) :: plane_normal(3), plane_dist, distance_threshold
   type(InterfacePoint), pointer, intent(out) :: neg(:), pos(:)

!  INPUTS
!   o points -- sampling points
!   o plane_normal -- normal vector of the bisecting plane
!   o plane_dist -- distance of the bisecting plane from the origin
!   o distance_threshold -- if point is closer to plane than this,
!                           it is discarded
!  OUTPUT
!   o points -- sampling points, resorted
!   o neg -- pointer to part with negative distances in points
!   o pos -- pointer to part with positive distances in points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: il, ir, l0, r0
   real(dp) :: dist
   real(dp), parameter :: ZERO = 0.e0_dp

   il = lbound(points,1)
   ir = ubound(points,1)
   l0 = il
   r0 = ir
   outer: do
      ! A
      fromleft: do while (il.le.ir)
         dist = calculate_distance_to_plane(points(il)%coord, &
                  plane_normal, plane_dist)
         if (abs(dist).lt.distance_threshold) then
            call swap_points(points, il, l0)
            l0 = l0 + 1
         elseif (dist.lt.ZERO) then
            exit fromleft
         endif
         il = il + 1
      enddo fromleft
      ! B
      fromright: do while (il.le.ir)
         dist = calculate_distance_to_plane(points(ir)%coord, &
                  plane_normal, plane_dist)
         if (abs(dist).lt.distance_threshold) then
            call swap_points(points, ir, r0)
            r0 = r0 - 1
         elseif (dist.gt.ZERO) then
            exit fromright
         endif
         ir = ir - 1
      enddo fromright
      ! C
      if (il.lt.ir) then
         call swap_points(points, il, ir)
         il = il + 1
         ir = ir - 1
      else
         exit outer
      endif
   enddo outer
   ! D
   neg => points(il:r0)
   pos => points(l0:ir)


   contains


   subroutine swap_points(points, ia, ib)
      implicit none
      type(InterfacePoint), intent(inout) :: points(:)
      integer, intent(in) :: ia, ib

      type(InterfacePoint) :: tmp
      tmp = points(ib)
      points(ib) = points(ia)
      points(ia) = tmp
   end subroutine swap_points

end subroutine divide_points_with_plane

!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/initialize_isc_keyword_defaults
!  NAME
!    initialize_isc_keyword_defaults
!  SYNOPSIS

subroutine initialize_isc_keyword_defaults()

!  PURPOSE
!    This subroutine initializes all implicit solvent cavity keywords
!    that are not of logical type, unless they have been specified by the user.
!
!    Please place sensible default values for all non-logical type keywords
!    here. Sanity checks for user input should go to read_control.f90
!
!  USES
   implicit none
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'initialize_isc_keyword_defaults'

   ! KEYWORDS: set defaults if not specified by user

   ! cavity type
   if ( isc_cavity_type .eq. ISC_CONST % CAVITY_UNDEF ) then
      call aims_stop_coll("Cavity type has not been specified")
   endif

   call initialize_cd_keyword_defaults()

   if (ifp) then
      isc_surface_sampling = ISC_CONST%CAV_SAMP_Even
   else
      isc_surface_sampling = ISC_CONST%CAV_SAMP_EqSph
   endif

end subroutine initialize_isc_keyword_defaults




!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/get_new_grid
!  NAME
!    get_new_grid
!  SYNOPSIS

subroutine get_new_grid(n_requested, interfaces, charlensq, distsq_threshold)

!  PURPOSE
!    This subroutine creates a new spherical Lebedev-Laikov grid with at least
!    {n_requested} points (max. 5810 per atom).
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_requested(:)

   type(DielectricInterface), intent(inout) :: interfaces(:)
   real(dp), intent(out) :: charlensq
   real(dp), intent(out) :: distsq_threshold

!  INPUTS
!   o n_requested -- number of points on cavity to generate
!  OUTPUT
!   o interfaces -- dielectric interfaces for this problem;
!                   only first interface is filled, sorting will be done later
!   o charlensq -- characteristic length, see comment below
!   o distsq_threshold -- distance threshold (squared) for neighboring points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   real(dp), parameter :: min_dist = 5.e-3_dp

   character(*), parameter :: func = 'get_new_grid'
   character(132) :: info_str

   integer :: n_grid_atoms

   integer, allocatable :: n_grid_points_atom(:)
   integer, allocatable :: n_grid_atoms_tasks(:), atoms_offset_tasks(:)
   real(dp), allocatable :: x_atom(:,:,:), R_pow3(:)
   integer, allocatable :: indices_keep(:,:), n_keep(:)
   integer :: n_points_task

   real(dp) :: Rk_inv, minaccept, R_max, sum_R_pow3
   real(dp) :: akjoverRk(3), lhsvec(3), RjoverRk

   integer :: j_atom, k_atom, n_points, i_point, i_keep, i_if, i_task

   intrinsic abs


   R_max = maxval(R_atom)
   n_grid_atoms = size(R_atom)

   allocate(n_grid_points_atom(n_grid_atoms))

   ! Calculate approximate number of grid points per atom:
   !
   ! Assume that the total volume filled up by small spheres
   ! around all atoms has spherical shape. Now, we can calculate
   ! the ratio of the big sphere's surface and the surface of
   ! a small sphere, which is equal to the ratio of the total
   ! number of points and the number of points on a small sphere
   ! for an equal density of points.

   select case (isc_surface_sampling)
      case (ISC_CONST%CAV_SAMP_EqSph)
         !TODO maxval or sum???
         n_grid_points_atom(1) = ceiling( maxval(n_requested(:min(2,ubound(n_requested,1)))) / &
                                       n_grid_atoms**(2.e0_dp/3.e0_dp) )

         ! get valid number of beads on initial grid
         n_grid_points_atom(:) = get_valid_lebedev(n_grid_points_atom(1))
      case (ISC_CONST%CAV_SAMP_Even)
         allocate(R_pow3(n_grid_atoms))
         R_pow3 = R_atom * R_atom * R_atom
         sum_R_pow3 = sum(R_pow3)
         do j_atom = 1, n_grid_atoms
            n_grid_points_atom(j_atom) = ceiling( maxval(n_requested(:min(2,ubound(n_requested,1)))) * &
                                       (R_pow3(j_atom) / sum_R_pow3)**(2.e0_dp/3.e0_dp) )
            n_grid_points_atom(j_atom) = get_valid_lebedev(n_grid_points_atom(j_atom))
         enddo
         deallocate(R_pow3)
      case (ISC_CONST%CAV_SAMP_Inval)
         call aims_stop('*** Internal error: Invalid cavity sampling type', func)
      case (ISC_CONST%CAV_SAMP_Undef)
         call aims_stop('*** Internal error: Cavity sampling undefined', func)
      case default
         call aims_stop('*** Internal error: Unknown cavity sampling type', func)
   end select

   ! PRE-ALLOCATION
   allocate(n_grid_atoms_tasks(0:n_tasks-1))
   allocate(atoms_offset_tasks(0:n_tasks-1))

   ! calculate center load-distribution
   call distribute_loads(n_grid_atoms, n_tasks, n_grid_atoms_tasks, &
                           atoms_offset_tasks)

   ! ALLOCATION
   allocate(x_atom(3,maxval(n_grid_points_atom),n_grid_atoms))
   x_atom(:,:,:) = 0.e0_dp
   allocate(n_keep(n_grid_atoms))
   allocate(indices_keep(maxval(n_grid_points_atom), &
               atoms_offset_tasks(myid)+1:&
               atoms_offset_tasks(myid)+n_grid_atoms_tasks(myid)))

   ! CREATE INITIAL SPHERICAL GRID
!   call create_unscaled_lebedev_mesh(n_grid_points_atom, x_atom)

   ! COUNT RELEVANT POINTS

   ! A grid point x_k at atom k at position a_k is given by the formula
   !
   !     x_k = R_k * x_0 + a_k                                (1)
   !
   ! where R_k is the radius of the atomic grid.
   !
   ! Exclusion:
   ! We can discard any point that is within the radius R_j of any other
   ! atom j/=k.
   !
   !     || x_k - a_j || <= R_j                               (2a)
   !
   ! Substituting (1) in (2a) we get
   !
   !     || R_k * x_0 + a_k - a_j || <= R_j                   (3a)
   !
   ! The triangle inequality yields
   !
   !     || R_k * x_0 + a_k - a_j || <= || R_k * x_0 || + || a_k - a_j ||
   !
   ! Thus if
   !
   !     || R_k * x_0 || + || a_j - a_k || <= R_j             (4a)
   !
   ! we can exclude all points x_k. Making use of the fact that the norm
   ! of x_0 is 1, we obtain
   !
   !     || a_j - a_k || <= R_j - R_k                         (5a)
   !
   ! which is sufficient but not necessary to exlude all points.
   !
   ! Inclusion:
   ! Any point that is not within the grid radius R_j of any other atom j/=k
   ! at position a_j has to fullfil:
   !
   !     || x_k - a_j || > R_j                                (2b)
   !
   ! Substituting (1) in (2b) we get:
   !
   !     || R_k * x_0 + a_k - a_j || > R_j                    (3b)
   !
   ! The triangle inequality yields
   !
   !     || R_k * x_0 + a_k - a_j || >= | || R_k * x_0 || - || a_j - a_k || |
   !
   ! thus it is sufficient (but not necessary) that Eq (5) holds.
   !
   !     | || R_k * x_0 || - || a_j - a_k || | > R_j          (5b)
   !
   ! further, R_k is a scalar and the norm of x_0 is 1,
   ! and we can also assume || a_j - a_k || > R_k, so we obtain
   !
   !     || a_j - a_k || > R_j + R_k                          (6b)
   !
   ! In case Eq. (6) holds we can directly accept all points wrt. atom j.
   ! Otherwise, Eq. (3) has to be evaluated individually.

   ! prepare for sync
   n_keep = 0

   ! loop over all atoms k on this task
   loop_k: do k_atom = atoms_offset_tasks(myid)+1, &
                        atoms_offset_tasks(myid)+n_grid_atoms_tasks(myid)

      call create_unscaled_lebedev_mesh(n_grid_points_atom(k_atom), &
              x_atom(:,1:n_grid_points_atom(k_atom),k_atom))

      Rk_inv = 1.e0_dp / R_atom(k_atom)

      ! initialize point count and indices
      n_keep(k_atom) = n_grid_points_atom(k_atom)
      do i_point = 1, n_grid_points_atom(k_atom), 1
         indices_keep(i_point,k_atom) = i_point
      enddo ! i_point

      ! loop over all other atoms j /= k
      loop_j: do j_atom = 1, n_grid_atoms, 1
         if (j_atom.eq.k_atom) &
            cycle loop_j

         ! get ratio of radii
         RjoverRk = R_atom(j_atom) * Rk_inv

         ! get difference vector between atoms divided by R_k
         akjoverRk = ( grid_atoms(:,j_atom) - grid_atoms(:,k_atom) ) * Rk_inv

         ! CASE 0: Eq. (5a)
         ! The distance between the two atoms is smaller than the
         ! absolute difference of the radii -> discard all points
         ! and exit loop_j
         ! remark: not possible for same radii, still unlikely otherwise
         if (dot_product(akjoverRk,akjoverRk).le.(RjoverRk - 1.e0_dp)**2 &
                 .and. RjoverRk .gt. 1.e0_dp) then
            ! discard all points, finished
            n_keep(k_atom) = 0
            exit loop_j          ! ATTENTION JUMP
         endif

         ! CASE 1: Eq. (6b)
         ! The distance between the two atoms is lager than the
         ! sum of their radii -> accept all points for this atom j
         ! remark: equality is excluded to avoid identical points
         if (dot_product(akjoverRk,akjoverRk).gt.(RjoverRk + 1.e0_dp)**2) then
            ! accept all points for this atom j
            cycle loop_j         ! ATTENTION JUMP
         endif

         ! get squared right hand side of Eq. (3)
         minaccept = (RjoverRk + min_dist*Rk_inv)**2

         ! CASE 2: Eq. (3b)
         ! neither CASE 0 nor CASE 1, so calculate every point individually
         ! loop over all remaining atomic grid points
         do i_point = n_keep(k_atom), 1, -1
            lhsvec = x_atom(:,indices_keep(i_point,k_atom),k_atom) - akjoverRk(:)
            ! remark: equality is included to avoid identical points
            if (dot_product(lhsvec,lhsvec).le.minaccept) then
               ! discard point:
               ! swap its index with the last (accepted) one
               call iswap(indices_keep(i_point,k_atom), &
                          indices_keep(n_keep(k_atom),k_atom))
               ! now reduce count
               n_keep(k_atom) = n_keep(k_atom) - 1
            endif
         enddo ! i_point

         if (n_keep(k_atom).eq.0) &
            exit loop_j          ! ATTENTION JUMP
      enddo loop_j ! j_atom
   enddo loop_k ! k_atom


   ! SYNCHRONIZE
   call sync_vector_integer(n_keep, n_grid_atoms)

   ! get total number of points
   n_points = sum(n_keep)

   ! ALLOCATION
   do i_if = 1, size(interfaces), 1
      if (allocated(interfaces(i_if)%p)) &
         deallocate(interfaces(i_if)%p)
   enddo ! i_if
   allocate(interfaces(1)%p(n_points))
   do i_if = 2, size(interfaces), 1
      allocate(interfaces(i_if)%p(0))
   enddo ! i_if

   ! prepare sync
   do i_point = 1, n_points, 1
      interfaces(1)%p(i_point)%coord = 0.e0_dp
      interfaces(1)%p(i_point)%normal = 0.e0_dp
      interfaces(1)%p(i_point)%tangents = 0.e0_dp
   enddo ! i_point

   ! FILL ARRAYS
   i_point = sum(n_keep(:atoms_offset_tasks(myid)))
   ! loop over all atoms k on this task
   do k_atom = atoms_offset_tasks(myid)+1, &
               atoms_offset_tasks(myid)+n_grid_atoms_tasks(myid)
      ! loop over points
      do i_keep = 1, n_keep(k_atom), 1
         i_point = i_point + 1
         interfaces(1)%p(i_point)%coord = R_atom(k_atom) * &
                        x_atom(:,indices_keep(i_keep,k_atom),k_atom) + &
                        grid_atoms(:,k_atom)
         interfaces(1)%p(i_point)%normal = &
                                       -x_atom(:,indices_keep(i_keep,k_atom),k_atom)
      enddo ! i_keep
   enddo ! k_atom

   ! SYNCHRONIZATION
   i_point = 0
   do i_task = 0, n_tasks-1, 1
      n_points_task = sum(n_keep(atoms_offset_tasks(i_task)+1:&
                     atoms_offset_tasks(i_task)+n_grid_atoms_tasks(i_task)))
      call InterfacePoint_vector_mpi_bcast( &
                  interfaces(1)%p(i_point+1:i_point+n_points_task), &
                  i_task, mpi_comm_global )
      i_point = i_point + n_points_task
   enddo ! i_task

   ! PARTIAL DEALLOCATION
   deallocate(x_atom)
   deallocate(n_keep)
   deallocate(indices_keep)
   deallocate(atoms_offset_tasks)
   deallocate(n_grid_atoms_tasks)

   ! CHARACTERISTIC LENGTH
      ! We define a characteristic length in terms of the area that can
      ! be assign to each bead on the initial sphere with radius R.
      !
      ! The surface area of the whole sphere is 4pi*R^2, so a single
      ! bead occupies a flat area of 4pi*R^2/N.
      ! The radius of a circle of corresponding area is 2*R/sqrt(N)
      ! we define the characteristic length "charlen" as the diameter of
      ! the circle
   charlensq = 4.e0_dp * ( 4.e0_dp * maxval(R_atom**2 / n_grid_points_atom ) )

   ! For interface generation: Estimate for smallest area element is 4pi*R^2/N
   ! calculate edge length of a square with the same area
   charlen_for_if_gen = 2.e0_dp * sqrt(pi) * minval(R_atom / sqrt(real(n_grid_points_atom)))

   distsq_threshold = 2.0e1_dp * charlensq

   ! Final report
   write(info_str,'(6X,A)') &
      "Initialized cavity with superposed spheres"
   call localorb_info(info_str)

   if (.not. ifp) then
           ! Atomic positions and radii are stored in a module variable for interface calculations.
           ! For interface calculations with isodensity cavity, deallocation will be done
           ! after interface initialization, for overlapping spheres after first update.
           ! Otherwise, deallocate here.
           deallocate(R_atom)
           deallocate(grid_atoms)
   endif

   deallocate(n_grid_points_atom)

   contains


   pure function get_valid_lebedev(N_min)

      implicit none

      integer, intent(in) :: N_min
      integer :: get_valid_lebedev

      integer, parameter :: max_i_lebedev = 32
      integer, dimension(max_i_lebedev), parameter :: valid_lebedev_N = (/ &
         6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974, &
         1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810 /)

      integer :: i_lebedev

      get_valid_lebedev = valid_lebedev_N(max_i_lebedev)

      do i_lebedev = 1, max_i_lebedev, 1
         if (valid_lebedev_N(i_lebedev).ge.N_min) then
            get_valid_lebedev = valid_lebedev_N(i_lebedev)
            exit
         endif
      enddo

   end function get_valid_lebedev


   subroutine create_unscaled_lebedev_mesh(N, x)

      implicit none

      integer, intent(in) :: N ! number of points of Lebedev-Laikov grid
                               ! only certain Ns are accepted!
      real(dp), dimension(3,N), intent(out) :: x


      character(*), parameter :: func = 'create_unscaled_lebedev_mesh'

      ! additional variables for Lebedev grids
      real(dp), dimension(:), allocatable :: w
      integer :: n_leb

      ! external Lebedev-Laikov subroutines
      external LD0006,LD0014,LD0026,LD0038,LD0050,LD0074,LD0086,LD0110,LD0146,&
               LD0170,LD0194,LD0230,LD0266,LD0302,LD0350,LD0434,LD0590,LD0770,&
               LD0974,LD1202,LD1454,LD1730,LD2030,LD2354,LD2702,LD3074,LD3470,&
               LD3890,LD4334,LD4802,LD5294,LD5810

      ! adapted from grids.f90
      if (N.eq.6) then
         allocate(w(N))
         call LD0006(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.14) then
         allocate(w(N))
         call LD0014(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.26) then
         allocate(w(N))
         call LD0026(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.38) then
         allocate(w(N))
         call LD0038(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.50) then
         allocate(w(N))
         call LD0050(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.74) then
         allocate(w(N))
         call LD0074(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.86) then
         allocate(w(N))
         call LD0086(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.110) then
         allocate(w(N))
         call LD0110(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.146) then
         allocate(w(N))
         call LD0146(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.170) then
         allocate(w(N))
         call LD0170(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.194) then
         allocate(w(N))
         call LD0194(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.230) then
         allocate(w(N))
         call LD0230(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.266) then
         allocate(w(N))
         call LD0266(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.302) then
         allocate(w(N))
         call LD0302(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.350) then
         allocate(w(N))
         call LD0350(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.434) then
         allocate(w(N))
         call LD0434(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.590) then
         allocate(w(N))
         call LD0590(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.770) then
         allocate(w(N))
         call LD0770(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.974) then
         allocate(w(N))
         call LD0974(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.1202) then
         allocate(w(N))
         call LD1202(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.1454) then
         allocate(w(N))
         call LD1454(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.1730) then
         allocate(w(N))
         call LD1730(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.2030) then
         allocate(w(N))
         call LD2030(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.2354) then
         allocate(w(N))
         call LD2354(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.2702) then
         allocate(w(N))
         call LD2702(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.3074) then
         allocate(w(N))
         call LD3074(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.3470) then
         allocate(w(N))
         call LD3470(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.3890) then
         allocate(w(N))
         call LD3890(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.4334) then
         allocate(w(N))
         call LD4334(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.4802) then
         allocate(w(N))
         call LD4802(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.5294) then
         allocate(w(N))
         call LD5294(x(1,:),x(2,:),x(3,:),w,n_leb)
      else if (N.eq.5810) then
         allocate(w(N))
         call LD5810(x(1,:),x(2,:),x(3,:),w,n_leb)
      else
         call aims_stop("No suitable Lebedev grid found.", func)
      endif

      ! DEALLOCATION
      if (allocated(w)) deallocate(w)

   end subroutine create_unscaled_lebedev_mesh


end subroutine get_new_grid


!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/calculate_tangent_vectors
!  NAME
!    calculate_tangent_vectors
!  SYNOPSIS

subroutine calculate_tangent_vectors(points, preferred_direction)

!  PURPOSE
!    This subroutine determines two normal vectors (t1, t2) that are both
!    orthogonal to the normal direction and also mutually orthogonal.
!    They are designed such that the vectors (n, t1, t2) form a right-handed
!    coordinate system.
!    If possible, the first tangent vector is chosen such that it is
!    perpendicular to the "preferred_direction".
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), dimension(:), intent(inout) :: points
   real(dp), dimension(3), intent(in) :: preferred_direction

!  INPUTS
!   o points -- points (coordinates and normal vectors)
!   o preferred_direction -- the first tangent vector is chosen such that it is
!                            perpendicular to the "preferred_direction"
!  OUTPUT
!   o points -- points (coordinates, normals and tangents)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: i_point

   do i_point = lbound(points,1), ubound(points,1), 1
      ! construct the perpendicular in-plane vectors
      ! first in-plane vector
      points(i_point)%tangents(:,1) = get_perpendicular_d3( &
                     points(i_point)%normal, preferred_direction)
      ! second in-plane vector
      ! (normalized by construction, right-handed system)
      points(i_point)%tangents(:,2) = cross_d3(points(i_point)%normal, &
                     points(i_point)%tangents(:,1))
   enddo ! i_point

end subroutine calculate_tangent_vectors



!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/calculate_cavity_measures
!  NAME
!    calculate_cavity_measures
!  SYNOPSIS

subroutine calculate_cavity_measures(points, distsq_threshold, area, volume, &
      curvature_correction, mpi_comm)

!  PURPOSE
!     calculate surface area and volume of cavity
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(inout) :: points(:)
   real(dp), intent(in) :: distsq_threshold
   logical, intent(in) :: curvature_correction
   integer, intent(in), optional :: mpi_comm

   real(dp), intent(out) :: area, volume

!  INPUTS
!   o points -- sampling points with assigned coordinate systems
!   o distsq_threshold -- squared distance threshold for neighboring points
!  OUTPUT
!   o area -- surface area of cavity
!   o volume -- volume of cavity
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = 'calculate_cavity_measures'
   character(132) :: info_str

   type(NeighborList), allocatable :: nlists(:)

   real(dp) :: cpu_pvoronoi, clock_pvoronoi
   real(dp) :: cpu_neighbor, clock_neighbor

   allocate(nlists(size(points)))

   write(info_str,'(4X,A)') 'Projected Voronoi algorithm'
   call localorb_info(info_str)

   call get_timestamps(cpu_neighbor, clock_neighbor)
   call update_neighbor_lists_for_Voronoi(points, distsq_threshold, nlists)
   call get_times(cpu_neighbor, clock_neighbor)
   call output_times('6X', 'Updating neighbor lists', &
             cpu_neighbor, clock_neighbor)

   call get_timestamps(cpu_pvoronoi, clock_pvoronoi)
   call get_Voronoi_surface(points, nlists, &
         area=area, & !TODO DELETE ME
         curvature_correction=curvature_correction, &
         mpi_comm=mpi_comm)
   call get_volume_from_surface(points, &
         volume=volume, &
         origin=plane_origin, &
         mpi_comm=mpi_comm)
   call get_times(cpu_pvoronoi, clock_pvoronoi)
   call output_times('6X', 'pVoronoi routine', &
             cpu_pvoronoi, clock_pvoronoi)

   write(info_str,'(4X,A,F18.3,2X,A)') 'Total surface area:  ', &
                                 area*bohr*bohr, 'AA^2'
   call localorb_info(info_str)
   write(info_str,'(4X,A,F18.3,2X,A)') 'Total cavity volume: ', &
                                 volume*bohr*bohr*bohr, 'AA^3'
   call localorb_info(info_str)

   deallocate(nlists)

end subroutine calculate_cavity_measures



!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/update_neighbor_lists_for_Voronoi
!  NAME
!    update_neighbor_lists_for_Voronoi
!  SYNOPSIS

subroutine update_neighbor_lists_for_Voronoi(points, distsq_threshold, nlists)

!  PURPOSE
!     create neighbor lists for points on cavity for given sitance threshold
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(in) :: points(:)
   real(dp), intent(in) :: distsq_threshold

   type(NeighborList), intent(out) :: nlists(:)

!  INPUTS
!   o points -- sampling points with assigned coordinate systems
!   o distsq_threshold -- squared distance threshold for neighboring points
!  OUTPUT
!   o nlists -- information about neighbors for each point
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = 'update_neighbor_lists_for_Voronoi'
   character(132) :: info_str

   integer :: errors, warnings
   integer :: n_p, n_q
   integer :: i_p

   write(info_str,'(4X,A)') 'creating neighbor list ...'
   call localorb_info(info_str)

   errors = 0
   warnings = 0
   do i_p = 1, size(points)
      call update_neighbor_list_intra(points, i_p, distsq_threshold, &
                                      nlists(i_p))
      n_q = size(nlists(i_p)%nb)
      if ( n_q .lt. 3 ) then
         errors = errors + 1
      elseif ( n_q .lt. 5 ) then
         warnings = warnings + 1
      endif
   enddo ! i_p

   if ( errors.gt.0 ) then
      write(info_str,'(4X,A)') 'Errors detected, calculation is not '//&
                            'possible. STOPPING.'
      call localorb_info(info_str)
      call aims_stop('*** Too few neighbors in neighbor list', func)
   else if ( warnings.gt.0 ) then
      write(info_str,'(4X,A)') 'Warnings detected, '//&
                               'calculation might be inaccurate!'
      call localorb_info(info_str)
      write(info_str,'(4X,A)') 'You should verify the results '//&
                               'with a larger distance cutoff.'
      call localorb_info(info_str)
   endif

end subroutine update_neighbor_lists_for_Voronoi



!******
!-------------------------------------------------------------------------------
!****s* isc_implicit_solvent_cavity/update_implicit_solvent_cavity
!  NAME
!    update_implicit_solvent_cavity
!  SYNOPSIS

subroutine update_implicit_solvent_cavity(interfaces, never_updated_before, &
                     cavity_has_changed, cd_exitcode, mpi_comm)

!  PURPOSE
!    This subroutine brings the cavity back in order after constraint dynamics.
!    It updates tangents, areas and volumes and distributes points into subsets
!    for different dielectrics
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(inout) :: interfaces(:)
   logical, intent(in)  :: never_updated_before
   logical, intent(out) :: cavity_has_changed

   integer, intent(in), optional :: mpi_comm
   integer, intent(out) :: cd_exitcode

!  INPUTS
!   o interfaces -- dielectric interfaces with cavity entirely in interfaces(1)
!   o need_tangents_and_measures_update -- indiciator whether update of tangents
!                                          and measures is necessary in first
!                                          update and after cavity has changed
!   o mpi_comm -- should MPI be used, the MPI communicator
!  OUTPUT
!   o interfaces -- dielectric interfaces; in interface case cavity will now
!                   be split into interfaces(1) and interfaces(2)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(132) :: info_str

   cavity_has_changed = .false.

   ! X. gather cavity in first interface
   if (ifp) &
     call merge_interfaces(interfaces(1),interfaces(2))

  ! A. optimize cavity
   if (.not. isc_cavity_type .eq. ISC_CONST%CAVITY_OvlpSph) then
      write(info_str,'(6X,A)') "Details on optimization with "//&
                                 "multipole density"
      call localorb_info(info_str)
      write(info_str,'(8X,A,I6)') 'number of points: ', &
                                 size(interfaces(1)%p)
      call localorb_info(info_str)

      cd_exitcode = constraint_dynamics(interfaces(1)%p, center)

      if ( cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOSTEP ) then
         write(info_str,'(6X,A,A)') "Cavity is already sufficiently ", &
                     "close to the actual isodensity surface."
         call localorb_info(info_str)
      else
         cavity_has_changed = .true.

      endif
   else
      cd_exitcode = ISC_CONST%CDYN_EXIT_CONV ! Needs to have a value for the
                                             ! following if statements
   endif ! .not. isc_cavity_type .eq. ISC_CONST%CAVITY_OvlpSph

   ! If constraint dynamics did not converge, skip the rest and go back to
   ! update_dielectric_interfaces, where restart will be written and aims will
   ! abort.
   if (.not. cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOTCONV) then

      if (cavity_has_changed.or.never_updated_before) then
         ! Need to delete points close to interface BEFORE updating neighbor list, area of deleted points
         ! will be missing otherwise
         ! TODO: solve this more elegantly
         if (ifp) then
            call distribute_cavity_interface(interfaces, &
                                    exclude_points_on_interface_delta)
            call merge_interfaces(interfaces(1),interfaces(2))
         endif

         ! B. calculate tangent vectors
         call calculate_tangent_vectors(interfaces(1)%p, ifp_normal)

         ! C. calculate surface area and cavity volume
         if (isc_calculate_surface_and_volume) then
            call calculate_cavity_measures(interfaces(1)%p, &
               distsq_threshold=distsq_threshold, &
               area=isc_surface_area, &
               volume=isc_cavity_volume, &
               curvature_correction=isc_surface_curvature_correction, &
               mpi_comm=mpi_comm)
         endif
      endif

      ! X. distribute the cavity points among the actual interfaces
      if (ifp) &
         call distribute_cavity_interface(interfaces, &
                                       exclude_points_on_interface_delta)
   endif ! .not. cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOTCONV

end subroutine update_implicit_solvent_cavity



end module isc_implicit_solvent_cavity

