!****h* FHI-aims/mpe_dielectric_interfaces
!  NAME
!    mpe_dielectric_interfaces
!  SYNOPSIS

module mpe_dielectric_interfaces

!  PURPOSE
!    This module contains subroutines to handle generic dielectric interfaces
!    in the MPE implicit solvation model. Subroutines for specific types of
!    interfaces are contained in respective modules. Subroutines used in
!    multiple of the latter modules are contained in
!    mpe_dielectric_interfaces_common.f90
!
!    Schematically, the hierarchy of dielectric interface modules goes as
!    follows:
!
!    1.                      mpe_dielectric_interfaces
!                                  (this module)
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
!     o create_dielectric_interfaces
!       TASKS
!       - initialization of new interface grids either from restart file
!         or from scratch
!       - allocation of public module variables
!       CALL
!         once during the initialization of the implicit solvent model
!
!     o update_dielectric_interfaces
!       TASKS
!       - optimization of cavity to current density value
!       - calculation of cavity's surface and volume for non-electrostatic
!         contributions
!       - reassignment of points to correct subsets
!       - output of restart information (if desired)
!       CALL
!         once if cavity is kept static, otherwise every step
!         in the scf cycle.
!         (provides necessary information for the calculation of the
!         reaction field)
!
!     o cleanup_dielectric_interfaces
!       TASKS
!       - deallocation of public module variables
!       CALL
!         once during the final deallocations
!
!  USES

   use constants, only: bohr
   use mpi_tasks, only: aims_stop_coll, aims_stop, myid
   use types, only: dp
   use localorb_io, only: localorb_info
   use mpe_types, only: &
         DielectricInterface
   use dimensions, only: n_atoms
   use geometry, only: coords
   use runtime_choices, only: &
         ifp, &
         ifp_normal, &
         ifp_dist, &
         isc_cavity_restart_read, &
         isc_cavity_restart_read_file, &
         isc_cavity_type, &
         isc_cavity_restart_write, &
         isc_cavity_restart_write_file
   use mpe_constants, only: ISC_CONST
   use isc_implicit_solvent_cavity, only: &
         initialize_isc_keyword_defaults, &
         get_new_grid, &
         plane_origin, &
         update_implicit_solvent_cavity
   use isc_common, only: &
         charlensq, &
         distsq_threshold
   use ifp_dielectric_interface_plane, only: &
         create_interface_points_on_plane, &
         update_dielectric_interface_plane, &
         initialize_ifp_keyword_defaults
   use mpe_dielectric_interfaces_common, only: &
         dot_d3, &
         center
         
         
   
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

   ! PUBLIC subroutines
   public create_dielectric_interfaces
   public update_dielectric_interfaces
   public cleanup_dielectric_interfaces

   ! Note: module variables are initialized in set_mpe_defaults()
   logical, public :: never_updated_before


   
   contains


!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/create_dielectric_interfaces
!  NAME
!    create_dielectric_interfaces
!  SYNOPSIS

subroutine create_dielectric_interfaces(n_requested, interfaces)

!  PURPOSE
!    This subroutine initializes new dielectric interfaces either from a
!      restart file (if provided) or from scratch
!    In the latter case, the cavity will only be a superposition of spheres,
!      and in the interface case, points will not yet be assigned to the
!      correct subsets yet. Meaningful dielectric interfaces will only be
!      available after the first call to update_dielectric_interfaces
!
!  USES
   implicit none

!  ARGUMENTS
   integer, intent(in) :: n_requested(:)

   type(DielectricInterface), intent(inout) :: interfaces(:)
!  INPUTS
!   o n_requested -- desired number of points on cavity (index 1) and, in
!                    the interfaces case, on the dielectric interface (index 3)
!   o interfaces -- prepared data structure for dielectric interfaces
!  OUTPUT
!   o interfaces -- dielectric interfaces with sampling points added
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'create_dielectric_interfaces'

   character(132) :: info_str
   integer :: cd_exitcode, n_interfaces, n_cavity_points, i_center

   ! initialization of keywords (set defaults if not specified)
   call initialize_dielectric_interfaces_keyword_defaults()

   write(info_str,'(4X,A)') 'Implicit solvent cavity'
   call localorb_info(info_str)

   ! center of coordinates
   center = sum(coords, 2)/n_atoms
   if ( ifp ) then
           plane_origin = center + ifp_normal * &
                   ( ifp_dist - dot_d3(center,ifp_normal) )
   else
       plane_origin = center
   endif

   ! B. load or create sampling points
   if ( isc_cavity_restart_read ) then
      ! get restart information
      call import_interfaces_from_restart(isc_cavity_restart_read_file, &
                  interfaces, charlensq, distsq_threshold)
      write(info_str,'(6X,A,A)') "Read cavity from restart file ", &
         isc_cavity_restart_read_file
      call localorb_info(info_str)
   else
      ! no restart information, initialize from Lebedev spheres
      call get_new_grid(n_requested, interfaces, &
                        charlensq, distsq_threshold)
      if ( ifp ) then
         call create_interface_points_on_plane( n_requested(3), &
                 plane_origin, ifp_normal, interfaces(3)%p )
      endif
   endif

   n_interfaces = size(interfaces)
   n_cavity_points = size(interfaces(1)%p)

   if ( n_interfaces .eq. 4 ) &
      n_cavity_points = n_cavity_points + size(interfaces(2)%p)
   write(info_str,'(8X,A,I12)') &
      "points assigned to cavity: ", n_cavity_points
   call localorb_info(info_str)
   write(info_str,'(8X,A,E14.6)') &
      "characteristic length is ", sqrt(charlensq)
   call localorb_info(info_str)
   write(info_str,'(8X,A,E14.6)') &
      "neighbor distance cutoff is ", sqrt(distsq_threshold)
   call localorb_info(info_str)

end subroutine create_dielectric_interfaces



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/update_dielectric_interfaces
!  NAME
!    update_dielectric_interfaces
!  SYNOPSIS

subroutine update_dielectric_interfaces ( interfaces, &
            cavity_has_changed, mpi_comm )

!  PURPOSE
!    This subroutine should be called once before the first SCF cycle and,
!    in case the cavity should adapt to the actual current electron density (in
!    the multipole approximation), also every time before the reaction field
!    factors are calculated. It optimizes the position of all points to the
!    current density and, in the interface case, reassigns points to the correct
!    subsets.
!
!  USES
   implicit none

!  ARGUMENTS
   type(DielectricInterface), intent(inout) :: interfaces(:)
   logical, intent(out) :: cavity_has_changed

   integer, intent(in), optional :: mpi_comm
!  INPUTS
!   o interfaces -- dielectric interfaces
!   o mpi_comm -- should MPI be used, the MPI communicator
!  OUTPUT
!   o interfaces -- updated dielectric interfaces
!   o cavity_has_changed -- indicates whether the cavity has changed its shape
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'update_dielectric_interfaces'
   character(132) :: info_str
   integer :: cd_exitcode, i_if


   call update_implicit_solvent_cavity(interfaces, never_updated_before, &
             cavity_has_changed, cd_exitcode, mpi_comm=mpi_comm)

   if ( cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOTCONV ) then
      write(info_str,'(6X,A)') "constraint dynamics did not converge"
      call localorb_info(info_str)
      ! save last configuration
      if ( isc_cavity_restart_write ) then
         write(info_str,'(6X,A,A)') "writing cavity to restart file ", &
                                   isc_cavity_restart_write_file
         call localorb_info(info_str)
         call export_interfaces_to_restart(isc_cavity_restart_write_file, &
                        interfaces, charlensq, distsq_threshold)
     else
         write(info_str,'(6X,A)') "NOT writing cavity to restart file"
         call localorb_info(info_str)
      endif
      ! exit
      call aims_stop_coll('*** constraint dynamics did not converge!',&
                           func)
   endif

   ! A.1 remove points from cavity
   if ( ifp .and. (cavity_has_changed.or.never_updated_before)) then
      call update_dielectric_interface_plane(interfaces(3),interfaces(4))
   endif ! ifp

   !TODO integrate this somehow
   write(info_str,'(8X,A)') "cavity measures for interfaces"
   call localorb_info(info_str)
   write(info_str,'(8X,2X,2(1X,A18))') "area (AA^2)", "volume (AA^3)"
   call localorb_info(info_str)
   do i_if = 1, max(1,size(interfaces))
      write(info_str,'(8X,I2,2(1X,F18.7))') i_if, &
            sum(interfaces(i_if) % p(:) % area)*bohr*bohr, &
            sum(interfaces(i_if) % p(:) % volume)*bohr*bohr*bohr
      call localorb_info(info_str)
   enddo

   ! Z. write restart file
   if ( isc_cavity_restart_write ) then
      write(info_str,'(6X,A,A)') "writing cavity to restart file ", &
                                    isc_cavity_restart_write_file
      call localorb_info(info_str)
      call export_interfaces_to_restart(isc_cavity_restart_write_file, &
                     interfaces, charlensq, distsq_threshold)
   endif

   never_updated_before = .false.

end subroutine update_dielectric_interfaces


!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/cleanup_dielectric_interfaces
!  NAME
!    cleanup_dielectric_interfaces
!  SYNOPSIS

subroutine cleanup_dielectric_interfaces()

!  PURPOSE
!    This subroutine deallocates all module variables
!    TODO: I don't think it does!
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

   never_updated_before = .true.

end subroutine cleanup_dielectric_interfaces



!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/initialize_dielectric_interfaces_keyword_defaults
!  NAME
!    initialize_dielectric_interfaces_keyword_defaults
!  SYNOPSIS

subroutine initialize_dielectric_interfaces_keyword_defaults()

!  PURPOSE
!    This subroutine initializes all dielectric interface keywords
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

   character(*), parameter :: func = 'initialize_dielectric_interfaces_keyword_defaults'

   ! KEYWORDS: set defaults if not specified by user

   call initialize_isc_keyword_defaults()
   if (ifp) &
      call initialize_ifp_keyword_defaults()

end subroutine initialize_dielectric_interfaces_keyword_defaults

   
!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/export_interfaces_to_restart
!  NAME
!    export_interfaces_to_restart
!  SYNOPSIS

subroutine export_interfaces_to_restart(filename, interfaces, &
               charlensq, distsq_threshold)

!  PURPOSE
!    This subroutine writes all necessary restart information for the cavity
!    to a given file. The format of this file is very similar to the
!    established XYZ format and can be read by many viewers like
!    VMD or JMol
!
!  USES
   implicit none

!  ARGUMENTS
   character(*), intent(in) :: filename
   type(DielectricInterface), intent(inout) :: interfaces(:)
   real(dp), intent(in) :: charlensq
   real(dp), intent(in) :: distsq_threshold

!  INPUTS
!   o filename -- name of the output file, file will be overwritten
!   o interfaces -- dielectric interfaces, to be written to restart file
!   o charlensq -- characteristic length squared
!   o distsq_threshold -- neighbor distance threshold squared
!  OUTPUT
!    writes to file on disk
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   integer :: i_point, i_if, n_points
   real(dp), dimension(7) :: values

   n_points = 0
   do i_if = lbound(interfaces,1), ubound(interfaces,1), 1
      n_points = n_points + size(interfaces(i_if)%p)
   enddo

   if (myid.eq.0) then
      open(unit=147, file=trim(filename)//'.bin', status='unknown', action='write', form='unformatted')
      
      do i_if = lbound(interfaces,1), ubound(interfaces,1), 1
         do i_point = lbound(interfaces(i_if)%p,1), &
                      ubound(interfaces(i_if)%p,1), 1
            values(1:3)  = interfaces(i_if)%p(i_point)%coord
            values(4:6)  = interfaces(i_if)%p(i_point)%normal
            values(7)    = interfaces(i_if)%p(i_point)%area
            write(147) i_if, values
            values(1:3)  = interfaces(i_if)%p(i_point)%tangents(:,1)
            values(4:6)  = interfaces(i_if)%p(i_point)%tangents(:,2)
            values(7)    = interfaces(i_if)%p(i_point)%volume
            write(147) 0, values
         enddo ! i_point
      enddo ! i_if
      
      close(147)

      open(unit=20, file=filename, status='unknown', action='write')

      ! write number of points in first line (XYZ standard)
      write(20,*) n_points

      ! use second line (comments) to store additional information
      write(20,'(A,2(X,A,X,F18.6))') "#", "charlensq", charlensq, &
         "distsq_threshold", distsq_threshold

      ! remaining lines contain region, cartesian coordinates and normal vectors
      do i_if = lbound(interfaces,1), ubound(interfaces,1), 1
         do i_point = lbound(interfaces(i_if)%p,1), &
                      ubound(interfaces(i_if)%p,1), 1
            write(20,'(I2,X,6(1X,F12.6))') i_if, &
                              interfaces(i_if)%p(i_point)%coord*bohr, &
                              interfaces(i_if)%p(i_point)%normal

         enddo ! i_point
      enddo ! i_if

      close(20)
   endif

end subroutine export_interfaces_to_restart



!******
!-------------------------------------------------------------------------------
!****s* mpe_dielectric_interfaces/import_interfaces_from_restart
!  NAME
!    import_interfaces_from_restart
!  SYNOPSIS

subroutine import_interfaces_from_restart(filename, interfaces, &
                                       charlensq, distsq_threshold)

!  PURPOSE
!    This subroutine imports the restart information for the cavity
!    from a given (existing) file.
!    The output array interfaces is expected to be allocated,
!    but not the sampling points array. Allocation in proper size
!    is done in this routine.
!
!  USES
   implicit none

!  ARGUMENTS
   character(*), intent(in) :: filename

   type(DielectricInterface), intent(inout) :: interfaces(:)
   real(dp), intent(out) :: charlensq
   real(dp), intent(out) :: distsq_threshold

!  INPUTS
!   o filename -- name of the input file
!  OUTPUT
!   o interfaces -- dielectric interfaces contained in this restart file
!   o charlensq -- characteristic length squared
!   o distsq_threshold -- neighbor distance threshold squared
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'import_interfaces_from_restart'
   character(132) :: inputline
   character(132) :: info_str
   character(40) :: desc_str

   integer :: n_points, n_interfaces, i_if, i_point, i_dump
   integer, dimension(:), allocatable :: counters
   real(dp), dimension(7) :: values
   real(dp), dimension(6) :: values_humanreadable

   logical :: bin_exists

   n_interfaces = size(interfaces)

   allocate(counters(n_interfaces))

   ! this is done by all processes, race conditions shouldn't be a problem
   open(unit=20, file=filename, status='old', action='read')

   ! read number of total points
   read(20,'(A)') inputline
   read(inputline,*) n_points

   ! read characteristic length and distance threshold from comment line
   read(20,'(A)') inputline
   read(inputline,*) desc_str, desc_str, charlensq, desc_str, distsq_threshold

   ! read regions, check out-of-bounds
   counters = 0
   do i_point = 1, n_points, 1
      read(20,'(A)') inputline
      read(inputline,*) i_if, values_humanreadable
      if ( (i_if.lt.1) .or. (i_if.gt.n_interfaces) ) then
         write(info_str,'(2X,A)') 'Bad interface descriptor found in '//&
                                  'cavity restart file.'
         call localorb_info(info_str)
         write(info_str,'(2X,A)') 'Your restart file seems broken or '//&
                              'not suitable for this system. Aborting.'
         call localorb_info(info_str)
         call aims_stop('*** Bad interface descriptor found!', func)
      else
         counters(i_if) = counters(i_if) + 1
      endif
   enddo

   ! array allocation
   do i_if = 1, n_interfaces, 1
      if (allocated(interfaces(i_if)%p)) &
         deallocate(interfaces(i_if)%p)
      allocate(interfaces(i_if)%p(counters(i_if)))
   enddo ! i_if

   ! TODO: Now alloctes based on .xyz file and then reads from unformatted file
   !       Could be solved more elegantly

   close(20)

   inquire(file=trim(filename)//'.bin', exist=bin_exists)

   ! Import from unformatted file, if existent
   if (bin_exists) then
      open(unit=147, file=trim(filename)//'.bin', status='old', action='read', form='unformatted')

      counters = 0
      do i_point = 1, n_points, 1
         read(147) i_if, values
         counters(i_if) = counters(i_if) + 1
         interfaces(i_if)%p(counters(i_if))%coord         = values(1:3)
         interfaces(i_if)%p(counters(i_if))%normal        = values(4:6)
         interfaces(i_if)%p(counters(i_if))%area          = values(7)
         read(147) i_dump, values
         interfaces(i_if)%p(counters(i_if))%tangents(:,1) = values(1:3)
         interfaces(i_if)%p(counters(i_if))%tangents(:,2) = values(4:6)
         interfaces(i_if)%p(counters(i_if))%volume        = values(7)
      enddo

      close(147)

   ! Import from .xyz file (old version)
   else
      open(unit=20, file=filename, status='old', action='read')

      rewind(20)

      ! skip header lines
      read(20,'(A)') inputline
      read(20,'(A)') inputline

      counters = 0
      do i_point = 1, n_points, 1
         read(20,'(A)') inputline
         read(inputline,*) i_if, values_humanreadable
         counters(i_if) = counters(i_if) + 1
         interfaces(i_if)%p(counters(i_if))%coord         = values_humanreadable(1:3)/bohr
         interfaces(i_if)%p(counters(i_if))%normal        = values_humanreadable(4:6)
      enddo

      close(20)
   endif ! unformatted restart exists

   deallocate(counters)

end subroutine import_interfaces_from_restart



end module mpe_dielectric_interfaces
