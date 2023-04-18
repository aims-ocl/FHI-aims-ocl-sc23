!****h* FHI-aims/isc_implicit_solvent_cavity
!  NAME
!    isc_implicit_solvent_cavity
!  SYNOPSIS

module isc_projected_voronoi

!  PURPOSE
!     TODO
!  USES

   use types, only: dp
   use quicksort_index, only: dquicksort_indexlist
   use mpe_types, only: InterfacePoint
   use isc_common, only: NeighborList

   use mpi_tasks, only: aims_stop
   use synchronize_mpi_basic, only: sync_vector
   use localorb_io, only: localorb_info

   implicit none

   private

   public :: get_Voronoi_surface
   public :: get_volume_from_surface

contains

!-------------------------------------------------------------------------------
!****s* isc_projected_voronoi/get_volume_from_surface
!  NAME
!    get_volume_from_surface
!  SYNOPSIS

subroutine get_volume_from_surface( points, &
      volume, &
      mpi_comm, origin )

!  PURPOSE
!    This subroutine calculates signed incremental volume elements from 
!    the incremental surface areas assigned to all points.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(inout) :: points(:)

   real(dp), intent(out) :: volume

   integer, intent(in), optional :: mpi_comm
   real(dp), intent(in), optional :: origin(3)

!  INPUTS
!   o points -- point cloud with assigned coordinate systems and areas
!   o origin -- coordinate origin to define signed volumes
!               (default value: geometric mean of point coordinates)
!  OUTPUT
!   o points -- point cloud, now with assigned volume elements
!   o volume -- point cloud's volume
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   character(*), parameter :: func = 'get_volume_from_surface'
   character(132) :: info_str

   real(dp) :: volume_origin(3)
   real(dp) :: height

   real(dp), allocatable :: volume_elements(:)

   integer :: mpi_rank, mpi_size, mpi_ierr
   integer :: n_p
   integer :: i_p

   n_p = size(points)

   if (present(origin)) then
      volume_origin(:) = origin(:)
   else
      volume_origin(:) = 0.e0_dp 
      do i_p = 1, n_p
         volume_origin(:) = volume_origin(:) + points(i_p) % coord(:)
      enddo
      volume_origin(:) = volume_origin(:) / n_p
   endif

   if (present(mpi_comm)) then
      call MPI_Comm_rank(mpi_comm, mpi_rank, mpi_ierr)
      call MPI_Comm_size(mpi_comm, mpi_size, mpi_ierr)
   else
      mpi_rank = 0
      mpi_size = 1
   endif

   allocate(volume_elements(n_p))
   volume_elements = 0.e0_dp

   do i_p = mpi_rank+1, n_p, mpi_size
      ! The volume of a pyramid is
      ! V = g*h/3, where g is the area of the base and h the height.
      !
      ! the height is calculated from the Hesse normal form of the plane with
      ! normal vector n
      !  h = dot_product( n/||n||, o - p )
      ! with
      !  p: point on plane
      !  o: origin of coordinate system
      !
      ! The base area g is equal to the area of the Voronoi cell

      height = dot_product(volume_origin - points(i_p) % coord, &
                           points(i_p) % normal )
      volume_elements(i_p) = points(i_p) % area * height / 3.e0_dp
   enddo ! i_p

   ! SYNCHRONIZATION
   if (present(mpi_comm)) then
      call sync_vector(volume_elements, n_p, mpi_comm=mpi_comm)
   endif

   volume = sum(volume_elements)

   do i_p = 1, n_p
      points(i_p) % volume = volume_elements(i_p)
   enddo

   deallocate(volume_elements)

end subroutine get_volume_from_surface
!******

!-------------------------------------------------------------------------------
!****s* isc_projected_voronoi/get_Voronoi_surface
!  NAME
!    get_Voronoi_surface
!  SYNOPSIS

subroutine get_Voronoi_surface( points, nlists, &
      area, curvature_correction, mpi_comm )

!  PURPOSE
!    This subroutine approximates the surface area of a point cloud
!    via a local, two-dimensional Voronoi construction around all points.
!    From this, incremental surface elements are
!    defined and summed up to the total surface area.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(inout) :: points(:)
   type(NeighborList), intent(in) :: nlists(:)

   real(dp), intent(out), optional :: area

   logical, intent(in), optional :: curvature_correction
   integer, intent(in), optional :: mpi_comm

!  INPUTS
!   o points -- point cloud with assigned coordinate systems
!   o nlists -- information about neighbors in point cloud
!   o curvature_correction -- apply curvature correction (default: .true.)
!  OUTPUT
!   o points -- point cloud with assigned coordinate systems and areas
!   o area -- (optional) point cloud's surface area
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2017).
!  SOURCE

   logical, parameter :: default_curvature_correction = .true.

   real(dp), parameter :: collinear_threshold = 1.e-12_dp
   real(dp), parameter :: distance_threshold = 1.e-6_dp

   character(*), parameter :: func = 'get_Voronoi_surface'
   character(132) :: info_str

   logical :: use_curvature_correction

   real(dp), allocatable :: area_elements(:)
   real(dp), allocatable :: qtilde(:,:)
   real(dp), allocatable :: ztilde(:,:)
   real(dp), allocatable :: angles(:)
   real(dp), allocatable :: cos_opening_angles(:)
   integer, allocatable :: map_i_x_to_z(:)
   integer, allocatable :: ztilde_ind(:,:)

   real(dp) :: transform(2,3)
   real(dp) :: lambda, lambda_den, distance
   real(dp) :: point_area_doubled, curv_corr_factor, &
             mean_cos_opening_angle_half

   integer :: mpi_rank, mpi_size, mpi_ierr
   integer :: n_p, n_q
   integer :: i_p, i_q, i_l, i_x, i_z, i_z_pre
   integer :: max_q, max_z, n_z, n_x

   n_p = size(points)

   if (present(mpi_comm)) then
      call MPI_Comm_rank(mpi_comm, mpi_rank, mpi_ierr)
      call MPI_Comm_size(mpi_comm, mpi_size, mpi_ierr)
   else
      mpi_rank = 0
      mpi_size = 1
   endif

   if (present(curvature_correction)) then
      use_curvature_correction = curvature_correction
   else
      use_curvature_correction = default_curvature_correction
   endif
   curv_corr_factor = 1.e0_dp ! unchanged without curvature correction

   max_q = 0
   do i_p = 1, n_p
      max_q = max(max_q, size(nlists(i_p) % nb))
   enddo

   ! ALLOCATION
   allocate(area_elements(n_p))

   allocate(qtilde(2, max_q))
   allocate(cos_opening_angles(max_q))

   ! maximal number of intersections is M*(M-1)/2,
   ! where M is the number of neighbors
   max_z = (max_q*(max_q-1))/2
   allocate(ztilde(2, max_z))
   allocate(ztilde_ind(2, max_z))
   allocate(angles(max_z))
   allocate(map_i_x_to_z(max_z))


   ! INITIALIZATION
   area_elements = 0.e0_dp

   ! MPI parallelization implied
   do i_p = mpi_rank+1, n_p, mpi_size

      ! NOTE: we obtain an orthonormal basis from the in-plane vectors and
      ! the surface normal, therefore the basis transformation is a rotation
      ! and the transformation matrix "transform" is unitary.
      transform(1,:) = points(i_p) % tangents(:,1)
      transform(2,:) = points(i_p) % tangents(:,2)

      ! PART 1:
      ! determine all intersection points that possibly span the
      ! Voronoi cell

      n_z = 0
      ! loop over all neighbors
      n_q = size(nlists(i_p) % nb)
      do i_q = 1, n_q

         ! used for curvature correction, see comment below
         cos_opening_angles(i_q) = abs( dot_product( &
            points(i_p) % normal, &
            points(nlists(i_p)%nb(i_q)%ind) % normal ) )

         ! calculate coordinates in new basis, i.e.
         ! transform (= rotate) distance vectors
         qtilde(:,i_q) = matmul(transform, nlists(i_p) % nb(i_q) % relvec)

         ! The direction vectors of the perpendicular bisectors
         !   between point and neighbors:
         !   In 2D, (x, y) and (-y, x) form perpendicular vectors
         ! explicit construction via
         !   vtilde(1,i_q) = -qtilde(2,i_q)
         !   vtilde(2,i_q) =  qtilde(1,i_q)
         ! is not necessary

         ! calculate intersections of new bisector with all old bisectors
         do i_l = 1, i_q-1
            lambda_den = qtilde(1,i_q)*qtilde(2,i_l) &
                        - qtilde(1,i_l)*qtilde(2,i_q)
            ! lines do not intersect when this denominator of the full
            ! expression for lambda is (almost) zero
            if (abs(lambda_den) .gt. collinear_threshold) then
               lambda = ( qtilde(2,i_l)*(qtilde(2,i_l)-qtilde(2,i_q)) &
                        + qtilde(1,i_l)*(qtilde(1,i_l)-qtilde(1,i_q)) ) &
                     / lambda_den

               ! store intersection
               n_z = n_z + 1
               ztilde(1,n_z) = 0.5e0_dp*(qtilde(1,i_q) - lambda*qtilde(2,i_q))
               ztilde(2,n_z) = 0.5e0_dp*(qtilde(2,i_q) + lambda*qtilde(1,i_q))
               ztilde_ind(1,n_z) = i_l
               ztilde_ind(2,n_z) = i_q
            endif
         enddo ! i_l
      enddo ! i_q


      ! PART 2:
      ! Evaluate intersection points and keep only those, that are on the
      !   same side of all bisecting lines as the origin.
      !   This translates to: (signed) distance to all lines <= 0.
      ! These intersections span the so-called Voronoi cell.
      !
      ! For the subsequent triangulation, the remaining
      !   intersections are sorted according to their position vector angle
      !   The number of intersections is "n_z"
      !   The number of Voronoi cell corners is "n_x"
      n_x = 0
      loop_intersections: do i_z = 1, n_z
         ! determine distance of intersection to all lines
         do i_q = 1, n_q
            distance = dot_product( qtilde(:,i_q), &
                  ztilde(:,i_z) - 0.5e0_dp*qtilde(:,i_q) )
            ! "distance" > 0 means that intersection is not in Voronoi cell
            ! (see comment above). A small threshold is necessary to
            ! compensate for numerical inaccuracies
            if ( distance .gt. distance_threshold ) cycle loop_intersections
         enddo ! i_q

         ! corner identified
         n_x = n_x + 1
         map_i_x_to_z(n_x) = i_z
         ! order the intersection points
         angles(i_z) = atan2(ztilde(2,i_z),ztilde(1,i_z))
         !i_quadrant = ceiling(2*angles(n_x)/pi)+2
      enddo loop_intersections ! i_z

      ! check if there are sufficient spanning points of the Voronoi cell
      ! the minimum is three points that span a triangle
      ! however, usually there should be many more points and
      ! further checks might be helpful
      if ( n_x .lt. 3 ) then
         write(info_str,'(A,I3)') '*** ERROR: not enough spanning points '//&
               'for Voronoi cell of point ', i_p
         call localorb_info(info_str)
         write(info_str,'(4X,A,I3,A)') 'only ', n_x, ' corner points'
         call localorb_info(info_str)
         call aims_stop('Not enough spanning points.', func)
      end if

      ! sort Voronoi cell corners
      call dquicksort_indexlist(max_z, angles, &
               map_i_x_to_z, first=1, last=n_x)

      ! Now calculate area
      point_area_doubled = 0.e0_dp
      i_z_pre = map_i_x_to_z(n_x)
      do i_x = 1, n_x
         i_z = map_i_x_to_z(i_x)
         if (use_curvature_correction) then
            mean_cos_opening_angle_half = .125e0_dp * ( &
               cos_opening_angles(ztilde_ind(1,i_z)) + &
               cos_opening_angles(ztilde_ind(2,i_z)) + &
               cos_opening_angles(ztilde_ind(1,i_z_pre)) + &
               cos_opening_angles(ztilde_ind(2,i_z_pre)) )
            curv_corr_factor = 2.e0_dp/( 1.e0_dp + &
                              sqrt(.5e0_dp + mean_cos_opening_angle_half) )
         endif
         point_area_doubled = point_area_doubled + curv_corr_factor * &
                  ( ztilde(1,i_z_pre)*ztilde(2,i_z) &
                  - ztilde(2,i_z_pre)*ztilde(1,i_z) )
         i_z_pre = i_z
      enddo


      ! SOME DERIVATIONS/COMMENTS

      ! motivated by the ratio between the area S of a spherical calotte
      ! and a flat circle F of same segment radius, we formulate the following
      ! correction:
      !
      ! S / F = 2 / (1+cos(phi/2)) = 2 / ( 1 + sqrt( (1+cosphi)/2 ) )
      !
      ! where the opening angle phi is approximately given by the normal vector
      ! on point n and the normal vector on a neighbor m
      !
      ! cos(phi) = abs( dot_product( normals(n), normals(m) ) )
      !
      ! we simply average over the opening angles with all neighbors and
      ! calculate the correction term "curv_corr_factor" (see above)

      area_elements(i_p) = 0.5e0_dp * point_area_doubled

   enddo ! i_p

   ! SYNCHRONIZATION
   if (present(mpi_comm)) then
      call sync_vector(area_elements, n_p, mpi_comm=mpi_comm)
   endif

   do i_p = 1, n_p
      points(i_p) % area = area_elements(i_p)
   enddo

   if (present(area)) area = sum(area_elements)

   ! DEALLOCATION
   if ( allocated(area_elements) ) deallocate(area_elements)
   if ( allocated(map_i_x_to_z) ) deallocate(map_i_x_to_z)
   if ( allocated(angles) ) deallocate(angles)
   if ( allocated(ztilde_ind) ) deallocate(ztilde_ind)
   if ( allocated(ztilde) ) deallocate(ztilde)
   if ( allocated(cos_opening_angles) ) deallocate(cos_opening_angles)
   if ( allocated(qtilde) ) deallocate(qtilde)

end subroutine get_Voronoi_surface
!******


end module

