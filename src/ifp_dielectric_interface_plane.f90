!****h* FHI-aims/mpe_dielectric_interfaces_common
!  NAME
!    mpe_dielectric_interfaces_common
!  SYNOPSIS

module ifp_dielectric_interface_plane

!  PURPOSE
!    This module contains subroutines to define an interface plane between two
!    solvents in the MPE implicit solvation model. Subroutines shared with
!    isc_implicit_solvent_cavity are contained in
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
!                                                    (this module)
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
!     o create_interface_points_on_plane
!       TASKS
!       - initialization of new interface grid
!       CALL
!         once during the initialization of the implicit solvent model
!
!     o update_dielectric_interface_plane
!       TASKS
!       - redistribute interface points into one subset inside the cavity
!         and one subset outside of it
!       CALL
!         during update_dielectric_interfaces
!
!     o initialize_ifp_keyword_defaults
!       TASKS
!       - intialization of keywords if not set by the user
!       CALL
!         once during initialize_dielectric_interfaces_keyword_defaults
!
!  USES

   use constants, only: bohr, pi
   use mpi_tasks, only: aims_stop_coll, aims_stop
   use localorb_io, only: localorb_info
   use types, only: dp
   use mpe_types, only: &
         InterfacePoint, &
         InterfacePoint_vector_extract, &
         DielectricInterface
   use mpe_constants, only: ISC_CONST
   use runtime_choices, only: &
         ifp_rmin, &
         ifp_rmax, &
         ifp_manual, &
         ifp_rlog, &
         ifp_area_per_point, &
         ifp_n_shells_log, &
         ifp_n_angular_max, &
         ifp_min_dr, &
         ifp_n_shells_lin, &
         ifp_dist, &
         ifp_normal, &
         mpe_lmax_ep, &
         isc_cavity_type, &
         isc_isodensity_value
   use mpe_dielectric_interfaces_common, only: &
         get_perpendicular_d3, &
         cross_d3, &
         charlen_for_if_gen, &
         center, &
         R_atom, &
         grid_atoms, &
         iswap, &
         normalized_d3, &
         merge_interfaces, &
         calculate_distance_to_plane
   use rho_multipole_evaluation, only: &
         get_free_rho, &
         get_rho_multipole

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
   public create_interface_points_on_plane
   public initialize_ifp_keyword_defaults
   public update_dielectric_interface_plane

   contains


!******
!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/create_interface_points_on_plane
!  NAME
!    create_interface_points_on_plane
!  SYNOPSIS

subroutine create_interface_points_on_plane(n_desired, &
               plane_origin, normal_direction, points)

!  PURPOSE
!    This subroutine initializes a new representation of the interface
!      between two dielectrica from a simple 2D spherical mesh.
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_desired
   real(dp), intent(in) :: plane_origin(3), normal_direction(3)
   type(InterfacePoint), allocatable, intent(out) :: points(:)

!  INPUTS
!   o n_desired -- initial number of points on interface
!   o normal_direction -- normal direction of plane
!   o plane_origin -- origin of plane
!  OUTPUT
!   o points -- sampling points for dielectric interface
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'create_interface_points_on_plane'
   integer :: n_shells_log, n_shells_lin, n_angular_max, n_points_temp, n_points, i_point, i_atom, points_in_lin
   real(dp) :: rmin, rlog, rmax, dr, area_per_point_desired, dist, stretch_factor, center_to_plane
   real(dp), dimension(3) :: relvec
   real(dp), dimension(3,2) :: tangents
   type(InterfacePoint), allocatable :: points_temp(:)

   character(132) :: info_str

   if (ifp_manual) then

           ! Override automatic generation of grid densities
           if (ifp_min_dr .eq. -1.e0_dp .and. ifp_n_shells_lin .eq. -1) then
                   
                   ! Imported quantities
                   rlog = ifp_rlog
                   rmin = ifp_rmin
                   rmax = ifp_rmax
                   n_shells_log = ifp_n_shells_log
                   area_per_point_desired = ifp_area_per_point
                   n_angular_max = ifp_n_angular_max

                   ! Calculated quantities
                   dr = ( (rmax/rlog)**(1.0e0_dp/real(n_shells_log-1,dp)) &
                                 - 1.0e0_dp ) * rlog
                   n_shells_lin = int( (rlog-rmin)/dr )

           else if (ifp_rmax .eq. -1.0e0_dp .and. ifp_rlog .eq. -1.e0_dp) then
                   
                   ! Imported quantities
                   dr = ifp_min_dr
                   rmin = ifp_rmin
                   n_shells_lin = ifp_n_shells_lin
                   n_shells_log = ifp_n_shells_log
                   area_per_point_desired = ifp_area_per_point
                   n_angular_max = ifp_n_angular_max

                   ! Calculated quantities
                   rlog = dr * n_shells_lin + rmin
                   ! The following causes the linear regime to 'leak' into the nominally logarithmic region. 
                   ! However, it ensures the increasing ifp_n_shells_log only adds shells
                   ! and does not distort existing ones
                   rmax = (dr + rlog)**(n_shells_log-1) / rlog**(n_shells_log-2)

           else if (ifp_min_dr .eq. -1.e0_dp .and. ifp_rlog .eq. -1.e0_dp) then

                   ! TODO: implement this
                   call aims_stop_coll("Not implemented interface specification", func)

           else if (ifp_rmax .eq. -1.0e0_dp .and. ifp_n_shells_lin .eq. -1) then

                   ! Imported quantities
                   rlog = ifp_rlog
                   rmin = ifp_rmin
                   n_shells_log = ifp_n_shells_log
                   dr = ifp_min_dr
                   area_per_point_desired = ifp_area_per_point
                   n_angular_max = ifp_n_angular_max

                   ! Calculated quantities
                   n_shells_lin = int( (rlog-rmin)/dr )
                   ! The following causes the linear regime to 'leak' into the nominally logarithmic region. 
                   ! However, it ensures the increasing ifp_n_shells_log only adds shells
                   ! and does not distort existing ones
                   rmax = (dr + rlog)**(n_shells_log-1) / rlog**(n_shells_log-2)

           else if (ifp_n_shells_log .eq. -1 .and. ifp_n_shells_lin .eq. -1 ) then

                   ! Imported quantities
                   rlog = ifp_rlog
                   rmin = ifp_rmin
                   rmax = ifp_rmax
                   dr = ifp_min_dr
                   area_per_point_desired = ifp_area_per_point
                   n_angular_max = ifp_n_angular_max

                   ! Calculated quantities
                   n_shells_lin = nint( (rlog-rmin)/dr )
                   n_shells_log = nint( log(rmax/rlog)/(log(dr+rlog)-log(rlog)) ) + 1

           else
                   call aims_stop_coll("Not implemented interface specification", func)

           endif ! TODO: probably implement other combinations at some point...

           ! For array allocation
           n_points_temp = (n_shells_log + n_shells_lin) *n_angular_max

           write(info_str, '(4X,A)') 'ifp_manual was specified, overriding '//&
                   'mpe_degree_of_determination for the interface plane'
           call localorb_info(info_str)
   else
           center_to_plane = abs(calculate_distance_to_plane(center, &
                                        ifp_normal, ifp_dist))

           ! Estimate the area covered by the molecule from atomic radii
           rlog = 0.e0_dp
           do i_atom = 1, size(R_atom), 1
              relvec = grid_atoms(:,i_atom) - center

              ! Project out normal vector
              !relvec = relvec - dot_product(relvec, normal_direction)*normal_direction

              dist = R_atom(i_atom) + sqrt(dot_product(relvec, relvec))
              rlog = max(rlog, dist)
           enddo

           ! Approximately quadratic tesserae, which have approximately the same area elements
           ! as cavity surface points
           n_shells_lin  = int(2*rlog/charlen_for_if_gen) ! Not actually number of linear
                                                          ! shells at this point, v.i.

           ! rmin = dr/2, v.i.
           if (mod(n_shells_lin,2).eq.0) &
                   n_shells_lin = n_shells_lin + 1

           ! dr in consistency with other quantities
           dr = 2*rlog/n_shells_lin

           ! Equistant spacing at origin
           rmin = dr/2

           ! Until now was number of dr that fit in the diameter of a circle with radius rlog
           ! Will need actual number of shells later
           n_shells_lin = (n_shells_lin-1)/2

           area_per_point_desired = dr * dr

           ! Estimated number of points in equidistant region
           points_in_lin = nint(pi* (min(center_to_plane, rlog))**2 &
                   /area_per_point_desired)

           ! Reach n_angular max at innermost log shell
           n_angular_max = pi / area_per_point_desired * &
                           ((rlog+dr*.5e0_dp)**2 - (rlog-dr*.5e0_dp)**2)

           ! Assume molecule covers entire linear region
           n_shells_log = ceiling(real(n_desired-points_in_lin,dp) / real(n_angular_max,dp))
           rmax = (dr + rlog)**(n_shells_log-1) * rlog**(-n_shells_log+2)

           enforce_min_size: do
              if (rmax .ge. 4*rlog) exit enforce_min_size
              n_shells_log = n_shells_log + 1
              rmax = (dr + rlog)**(n_shells_log-1) * rlog**(-n_shells_log+2)
           enddo enforce_min_size
 
           ! For array allocation
           n_points_temp = (n_shells_log + n_shells_lin) *n_angular_max

           ! If cavity does not penetrate the interface, stretch the interface
           ! discretization

           ! rlog is estimated molecular radius
           stretch_factor = center_to_plane / rlog

           ! If molecule (in spherical approx.) does not penetrate the interface, stretch the grid
           if (stretch_factor .gt. 1.e0_dp) then
              write(info_str, '(A,E10.3)') 'Stretching interface grid by ', stretch_factor
              call localorb_info(info_str)
              rmin = stretch_factor * rmin
              rlog = stretch_factor * rlog
              rmax = stretch_factor * rmax
              dr   = stretch_factor * dr
              area_per_point_desired = area_per_point_desired * stretch_factor &
                                                              * stretch_factor
           endif

   endif

   if (.not. isc_cavity_type .eq. ISC_CONST % CAVITY_OvlpSph) then
           deallocate(R_atom)
           deallocate(grid_atoms)
   endif

   write(info_str, '(4X,A)') 'Generated interface plane with:'
   call localorb_info(info_str)
   write(info_str, '(6X,A,E14.7,A)') 'rmin          : ', rmin * bohr, ' AA'
   call localorb_info(info_str)
   write(info_str, '(6X,A,E14.7,A)') 'rmax          : ', rmax * bohr, ' AA'
   call localorb_info(info_str)
   write(info_str, '(6X,A,E14.7,A)') 'rlog          : ', rlog * bohr, ' AA'
   call localorb_info(info_str)
   write(info_str, '(6X,A,E14.7,A)') 'area_per_point: ', area_per_point_desired * bohr**2, ' AA^2'
   call localorb_info(info_str)
   write(info_str, '(6X,A,I10)')     'n_shells_log  : ', n_shells_log
   call localorb_info(info_str)
   write(info_str, '(6X,A,I10)')     'n_shells_lin  : ', n_shells_lin
   call localorb_info(info_str)
   write(info_str, '(6X,A,I10)')     'n_angular_max : ', n_angular_max
   call localorb_info(info_str)
   write(info_str, '(6X,A,E14.7,A)') 'min_dr        : ', dr * bohr, ' AA'
   call localorb_info(info_str)

   if (allocated(points)) &
      deallocate(points)
   allocate(points_temp(n_points_temp))

   ! define coordinate system on interface
   tangents(:,1) = get_perpendicular_d3(normal_direction)
   tangents(:,2) = cross_d3(normal_direction, tangents(:,1))

   call create_xy_loggrid( rmin, rlog, rmax, dr, n_shells_lin, n_shells_log, n_angular_max, &
          area_per_point_desired, points_temp, n_points )

   allocate(points(n_points))
   do i_point = 1, n_points, 1
       points(i_point) = points_temp(i_point)
   enddo
   deallocate(points_temp)

   call rotate_and_move_xy_grid(points, plane_origin, &
                                 normal_direction, tangents)

   contains


   subroutine create_xy_loggrid( rmin, rlog, rmax, dr, n_shells_lin, n_shells_log, n_angular_max, &
          area_per_point_desired, points_temp, n_points )

      implicit none

      real(dp), intent(in) :: rmin, rlog, rmax, area_per_point_desired
      integer, intent(in) :: n_shells_lin, n_shells_log, n_angular_max

      type(InterfacePoint), intent(out) :: points_temp(:)
      integer, intent(out) :: n_points

      real(dp), intent(inout) :: dr

      integer :: i_shell, i_angular, i_point, n_angular
      real(dp) :: angle, radius, radius_next, area, area_prev, area_next, area_per_point

      real(dp), parameter :: pi_quarter = pi*.25e0_dp

      character(132) :: info_str

      if (abs(dr - (rlog-rmin)/n_shells_lin) .ge. 1.e-15_dp) then ! TODO: numerical accuracy
              write(info_str, '(1X,A)') '* (MPE) Warning: Interface discretization parameters are '//&
                      'inconsistent. Overriding dr here.'
              call localorb_info(info_str)
              dr = (rlog-rmin)/(real(n_shells_lin,dp))
      endif

      n_angular = 2
      i_point             = 0
      radius_next         = rmin
      area_next           = 0.0e0_dp
      do i_shell = 0, n_shells_lin-1, 1
         area_prev      = area_next
         radius         = radius_next
         radius_next    = radius + dr
         ! A = pi * ( [1/2(r_{i+1}+r_{i})]^2 - [1/2(r_{i-1}+r_{i})]^2 ) / n
         area_next = pi_quarter*(radius+radius_next)**2
         area = area_next - area_prev
         n_angular = min(n_angular_max, ceiling(area/area_per_point_desired))
         area_per_point = area / n_angular
         do i_angular = 1, n_angular, 1
            i_point = i_point + 1
            angle = i_angular*2*pi/real(n_angular,dp)
            points_temp(i_point)%coord(1) = radius * cos(angle)
            points_temp(i_point)%coord(2) = radius * sin(angle)
            points_temp(i_point)%coord(3) = 0.e0_dp
            points_temp(i_point)%area     = area_per_point
            points_temp(i_point)%volume   = 0.e0_dp
         enddo ! i_angular
      enddo ! i_shell

      do i_shell = 0, n_shells_log-1, 1
         area_prev      = area_next
         radius         = radius_next
         ! Exponent is higher than formula in the paper by 1, because this is r_{i+1}
         radius_next         = (rmax/rlog)**((i_shell+1)/real(n_shells_log-1,dp))*rlog
         ! A = pi * ( [1/2(r_{i+1}+r_{i})]^2 - [1/2(r_{i-1}+r_{i})]^2 ) / n
         area_next = pi_quarter*(radius+radius_next)**2
         area = area_next - area_prev
         n_angular = min(n_angular_max, ceiling(area/area_per_point_desired))
         area_per_point = area / n_angular
         do i_angular = 1, n_angular, 1
            i_point = i_point + 1
            angle = i_angular*2*pi/real(n_angular,dp)
            points_temp(i_point)%coord(1) = radius * cos(angle)
            points_temp(i_point)%coord(2) = radius * sin(angle)
            points_temp(i_point)%coord(3) = 0.e0_dp
            points_temp(i_point)%area     = area_per_point
            points_temp(i_point)%volume   = 0.e0_dp ! Origin lies in plane
         enddo ! i_angular
      enddo ! i_shell

      n_points = i_point
   end subroutine


   subroutine rotate_and_move_xy_grid(points, dest, normal, tangents)
      implicit none
      type(InterfacePoint), intent(inout) :: points(:)
      real(dp), intent(in) :: dest(3), normal(3), tangents(3,2)

      integer :: i_ind
      real(dp)             :: rot_mat (3,3), quatnormsq, zplus1, zplus1sq, xsq, ysq
      real(dp), parameter  :: eps = 1.0e-7_dp
      character(132)       :: info_str

      intrinsic matmul

      zplus1 = normal(3) + 1.0e0_dp
      if ( zplus1 < eps ) then ! z- case
         write(info_str, '(4X, A)') 'Normal vector rotated into -z direction'
         call localorb_info(info_str)
         rot_mat(1,:)=(/ -1.0e0_dp,  0.0e0_dp,  0.0e0_dp /)
         rot_mat(2,:)=(/  0.0e0_dp, -1.0e0_dp,  0.0e0_dp /)
         rot_mat(3,:)=(/  0.0e0_dp,  0.0e0_dp, -1.0e0_dp /)
      else
         zplus1sq = zplus1**2
         xsq = normal(1)**2
         ysq = normal(2)**2
         quatnormsq = 1.0e0_dp / ( zplus1sq + xsq + ysq )

         rot_mat(1,:)=(/quatnormsq * (zplus1sq +ysq      -xsq      ), &
                        quatnormsq * (-2.0e0_dp*normal(1)*normal(2)), &
                        quatnormsq * ( 2.0e0_dp*zplus1   *normal(1)) /)
         rot_mat(2,:)=(/              rot_mat(1,2)                  , &
                        quatnormsq * (zplus1sq -ysq      +xsq      ), &
                        quatnormsq * ( 2.0e0_dp*zplus1   *normal(2)) /)
         rot_mat(3,:)=(/             -rot_mat(1,3)                  , &
                                     -rot_mat(2,3)                  , &
                        quatnormsq * (zplus1sq -ysq      -xsq      ) /)
      endif

      do i_ind = lbound(points,1), ubound(points,1), 1
         ! now comes the reduced multiplication of rotation matrix with grid
         points(i_ind)%coord = matmul(rot_mat, points(i_ind)%coord) + dest
         points(i_ind)%tangents = tangents
         points(i_ind)%normal = normal
      enddo ! i_ind
   end subroutine

end subroutine create_interface_points_on_plane


!******
!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/distribute_points
!  NAME
!    distribute_points
!  SYNOPSIS

subroutine distribute_points(points, remaining_points, &
                  values, upper_bound, lower_bound)

!  PURPOSE
!    This subroutine checks the array values and distributes all corresponding
!    coordinates into a set where the value is not within the specified bounds
!    and a remaining set.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), intent(inout), allocatable :: points(:)
   type(InterfacePoint), intent(out)  , allocatable :: remaining_points(:)
   real(dp), intent(in) :: values(:), upper_bound, lower_bound

!  INPUTS
!   o points -- sampling points
!   o values -- array of values for all points
!   o upper_bound -- point i is kept if values(i)<=upper_bound
!   o lower_bound --              .and. values(i)>=lower_bound
!  OUTPUT
!   o points -- points inside the boundaries
!   o remaining_points -- remaining points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: i_ind, i_tmp, i_tmprem, n_points, n_pointsrem
   type(InterfacePoint), dimension(:), allocatable :: tmp, tmprem

   n_points = size(points)
   allocate(tmp(n_points))
   allocate(tmprem(n_points))

   i_tmp = 0
   i_tmprem = 0
   do i_ind = 1, n_points, 1
      if ((values(i_ind).le.upper_bound).and. &
          (values(i_ind).ge.lower_bound)) then
         i_tmp = i_tmp + 1
         tmp(i_tmp) = points(i_ind)
      else
         i_tmprem = i_tmprem + 1
         tmprem(i_tmprem) = points(i_ind)
      endif
   enddo ! i_ind

   n_points = i_tmp
   n_pointsrem = i_tmprem
   if (allocated(remaining_points)) deallocate(remaining_points)
   deallocate(points)
   allocate(points(n_points))
   allocate(remaining_points(n_pointsrem))

   do i_ind = 1, n_points, 1
      points(i_ind) = tmp(i_ind)
   enddo ! i_ind
   do i_ind = 1, n_pointsrem, 1
      remaining_points(i_ind) = tmprem(i_ind)
   enddo

   deallocate(tmp)
   deallocate(tmprem)

end subroutine distribute_points



!******
!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/distribute_interface_points
!  NAME
!    distribute_interface_points
!  SYNOPSIS

subroutine distribute_interface_points(points_out, points_in)

!  PURPOSE
!    This subroutine distributes all sampling points on the interface
!    of two dielectrica into a set inside the cavity and a set outside.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), allocatable, intent(inout) :: points_out(:)
   type(InterfacePoint), allocatable, intent(out)   :: points_in(:)

!  INPUTS
!   o points_out -- sampling points on dielectric interface
!  OUTPUT
!   o points_out -- interface points outside the cavity
!   o points_in  -- 'ghost' interface points inside the cavity
!                   needed for interfacial energy
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   character(*), parameter :: func = 'distribute_interface_points'
   character(132) :: info_str
   integer :: n_points
   real(dp), allocatable :: rho(:), tmpx(:,:)

   n_points = size(points_out)
   allocate(rho(n_points))
   allocate(tmpx(3,n_points))

   call InterfacePoint_vector_extract(points_out, coords=tmpx)

   if ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPStat .or. &
          isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPDYN ) then
      call get_rho_multipole(n_points, tmpx, rho)
   elseif ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoFree ) then
      call get_free_rho(n_points, tmpx, rho)
   else
      ! raise error
      write(info_str,'(A)') '***Internal error: type of reference '//&
                            'density unknown'
      call localorb_info(info_str)
      call aims_stop('Internal error: type of reference density unknown', &
                           func)
   endif

   ! TODO: Density sometimes comes out negative. Should investigate if this is an actual bug!

   call distribute_points(points_out, points_in, rho, &
                  upper_bound=isc_isodensity_value, lower_bound=-1.e0_dp)

   deallocate(rho)
   deallocate(tmpx)

end subroutine distribute_interface_points


!******
!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/distribute_interface_points_OvlpSph
!  NAME
!    distribute_interface_points_OvlpSph
!  SYNOPSIS

subroutine distribute_interface_points_OvlpSph(points_out, points_in)

!  PURPOSE
!    This subroutine distributes all sampling points on the interface
!    of two dielectrica into a set inside the cavity and a set outside.
!    This is the version for the 'overlapping spheres' cavity type
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), allocatable, intent(inout) :: points_out(:)
   type(InterfacePoint), allocatable, intent(out)   :: points_in(:)

!  INPUTS
!   o points_out -- sampling points on dielectric interface
!  OUTPUT
!   o points_out -- interface points outside the cavity
!   o points_in  -- 'ghost' interface points inside the cavity
!                   needed for interfacial energy
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2018).
!  SOURCE

   character(*), parameter :: func = 'distribute_interface_points_OvlpSph'

   type(InterfacePoint), allocatable :: tmpp(:)

   integer, allocatable :: indices_keep(:)

   real(dp) :: minaccept

   real(dp) :: relvec(3)

   integer :: n_grid_atoms, n_keep, i_point, i_grid_atom, i_keep

   n_grid_atoms = size(R_atom)
   n_keep = size(points_out)
   allocate(indices_keep(n_keep))
   allocate(tmpp(n_keep), source=points_out)
   ! initialize count and indices
   do i_point = 1, n_keep, 1
      indices_keep(i_point) = i_point
   enddo ! i_point

   do i_grid_atom = 1, n_grid_atoms, 1
      ! check distance to plane
      if ( abs(dot_product(grid_atoms(:,i_grid_atom), ifp_normal) - &
              ifp_dist) .le. R_atom(i_grid_atom) ) then
         ! atom close enough, check distance to points
         minaccept = R_atom(i_grid_atom)**2
         do i_point = n_keep, 1, -1
            relvec = grid_atoms(:,i_grid_atom) - &
                                       points_out(indices_keep(i_point))%coord
            if ( dot_product(relvec,relvec).le.minaccept ) then
               call iswap(indices_keep(i_point),indices_keep(n_keep))
               n_keep = n_keep - 1
            endif
         enddo ! i_point
      endif
   enddo ! i_grid_atom

   ! copy remaining points to interface
   if (allocated(points_out)) &
      deallocate(points_out)
   allocate(points_out(n_keep))
   if (allocated(points_in)) &
      deallocate(points_in)
   allocate(points_in(size(tmpp)-n_keep))
   do i_keep = 1, n_keep, 1
      points_out(i_keep) = tmpp(indices_keep(i_keep))
   enddo
   do i_keep = 1, size(tmpp)-n_keep, 1
      points_in(i_keep) = tmpp(indices_keep(i_keep+n_keep))
   enddo

   deallocate(tmpp)
   deallocate(indices_keep)

   deallocate(R_atom)
   deallocate(grid_atoms)

end subroutine distribute_interface_points_OvlpSph


!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/initialize_ifp_keyword_defaults
!  NAME
!    initialize_ifp_keyword_defaults
!  SYNOPSIS

subroutine initialize_ifp_keyword_defaults()

!  PURPOSE
!    This subroutine initializes all dielectric interface plane keywords
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

   character(*), parameter :: func = 'initialize_ifp_keyword_defaults'

   ! KEYWORDS: set defaults if not specified by user

   ! dividing plane to simulate dielectric phase borders, right now only between
   ! a dielectric continuum and vacuum
   if ( ifp_rmin .le. 0.e0_dp ) &
           ifp_rmin = 1.e0_dp
   if ( ifp_rmax .le. 0.e0_dp .and. .not. ifp_manual) &
           ifp_rmax = 1.e2_dp

   ! normalize normal vector, just to be on the save side
   ifp_normal = normalized_d3(ifp_normal)

   if (ifp_manual) then
           if (ifp_area_per_point .eq. -1.e0_dp &
                   .or. ifp_n_angular_max .eq. -1) then
                   call aims_stop_coll("If ifp_manual is used, additional "//&
                          "parameters have to be given!", func)
          endif
          if (ifp_rmax .eq. -1.e0_dp &
                  .and. ifp_min_dr .eq. -1.e0_dp) then
                  call aims_stop_coll("If ifp_manual is used, min_dr "//&
                          "or rmax has to be given!", func)
          endif
          if (ifp_rlog .eq. -1.e0_dp &
                  .and. ifp_n_shells_lin .eq. -1) then
                  call aims_stop_coll("If ifp_manual is used, rlog "//&
                          "or n_shells_lin has to be given!", func)
          endif
  endif

end subroutine initialize_ifp_keyword_defaults



!******
!-------------------------------------------------------------------------------
!****s* ifp_dielectric_interface_plane/update_dielectric_interface_plane
!  NAME
!    update_dielectric_interface_plane
!  SYNOPSIS

subroutine update_dielectric_interface_plane(interface_out, interface_in)

!  PURPOSE
!    This subroutine redistributes interface points into a subset inside the
!    cavity and one outside
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(inout) :: interface_out, interface_in

!  INPUTS
!   o interface_out -- interface outside the cavity
!   o interface_in  -- interface inside the cavity
!  OUTPUT
!   o interface_out -- interface outside the cavity, updated
!   o interface_in  -- interface inside the cavity, updated
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   call merge_interfaces(interface_out,interface_in)
   if ( isc_cavity_type .eq. ISC_CONST % CAVITY_OvlpSph) then
           call distribute_interface_points_OvlpSph( &
                                        interface_out%p, interface_in%p)
   else
           call distribute_interface_points(interface_out%p, interface_in%p)
   endif

end subroutine update_dielectric_interface_plane



end module ifp_dielectric_interface_plane
