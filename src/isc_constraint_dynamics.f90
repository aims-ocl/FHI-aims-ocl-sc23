!****h* FHI-aims/isc_constraint_dynamics
!  NAME
!    isc_constraint_dynamics
!  SYNOPSIS

module isc_constraint_dynamics

!  PURPOSE
!    This module contains subroutines and functions used for the MD-like
!    isodensity cavity optimization procedure.
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
!                                                    (this module)
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
!    Short description of public subroutines and functions (please refer to
!    headers for more specific details):
!
!     o constraint_dynamics
!       TASKS
!       - Optimizes moves the cavity points to a specified isodensity value
!         in and MD like algorithm
!       CALLS
!         during update_implicit_solvent_cavity
!
!     o initialize_cd_keyword_defaults
!       TASKS
!       - initializes keyword defaults specific to this algorithm
!       CALL
!         once during initialize_isc_keyword_defaults
!
!  USES

   use types, only: dp
   use constants, only: bohr
   use localorb_io, only: localorb_info, OL_norm
   use timing, only: get_timestamps, get_times, output_times
   use mpi_tasks, only: myid, n_tasks, aims_stop, &
                        mpi_comm_global
   use synchronize_mpi_basic, only: sync_vector, &
                        sync_find_max, &
                        get_max_double
   use quicksort_index, only: dquicksort_indexlist
   use rho_multipole_evaluation, only: &
                                       get_rho_multipole_and_gradient, &
                                       get_free_rho_and_gradient
   use mpe_constants, only: ISC_CONST
   use mpe_types, only: &
         InterfacePoint, &
         InterfacePoint_vector_extract, &
         InterfacePoint_vector_mpi_bcast
   use runtime_choices, only: &
         isc_cavity_type, &
         isc_isodensity_value, &
         isc_record_cavity_creation_file, &
         isc_record_every_n_steps, &
         isc_kill_ratio, &
         isc_dt, &
         isc_rho_k, &
         isc_rho_k_mod_constant, &
         isc_g_k, &
         isc_rep_k, &
         isc_dynamics_friction, &
         isc_rho_rel_deviation_threshold, &
         isc_gradient_threshold, &
         isc_update_nlist_interval, &
         isc_max_dyn_steps, &
         isc_try_restore_convergence
   use mpe_dielectric_interfaces_common, only: &
         normalized_d3, &
         dot_d3
   use isc_common, only: &
         update_neighbor_list_intra, &
         distribute_loads, &
         distsq_threshold, &
         charlensq, &
         scale_interaction_distance, &
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
   public constraint_dynamics
   public initialize_cd_keyword_defaults

contains

!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/initialize_cd_keyword_defaults
!  NAME
!    initialize_cd_keyword_defaults
!  SYNOPSIS

subroutine initialize_cd_keyword_defaults()

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

   character(*), parameter :: func = 'initialize_cd_keyword_defaults'

   ! KEYWORDS: set defaults if not specified by user

   ! timestep
   if ( isc_dt .le. 0.e0_dp ) isc_dt = 1.e-1_dp

   ! force constants
   if ( isc_rho_k .le. 0.e0_dp ) isc_rho_k = 1.e0_dp
   if ( isc_rho_k_mod_constant .le. 0.e0_dp ) isc_rho_k_mod_constant = 1.0e-1_dp
   if ( isc_g_k .le. 0.e0_dp )   isc_g_k = 2.0e0_dp
   if ( isc_rep_k .le. 0.e0_dp ) isc_rep_k = 1.0e-2_dp

   ! cavity convergence parameters
   if ( isc_rho_rel_deviation_threshold .le. 0.e0_dp ) &
      isc_rho_rel_deviation_threshold = 1.e-3_dp

   ! gradient threshold
   if ( isc_gradient_threshold .le. 0.e0_dp ) isc_gradient_threshold = 1.0e-8_dp

   ! neighbor list
   if ( isc_update_nlist_interval .le. 0 ) isc_update_nlist_interval = 50

   ! maximal number of steps in constraint dynamics simulation
   if ( isc_max_dyn_steps .le. 0 ) isc_max_dyn_steps = 300

   ! pseduo-friction term, scales velocities in the update step
   if ( isc_dynamics_friction .le. 0.e0_dp ) isc_dynamics_friction = 0.2e0_dp

end subroutine initialize_cd_keyword_defaults

!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/constraint_dynamics
!  NAME
!    constraint_dynamics
!  SYNOPSIS

function constraint_dynamics(points, gravity_center) &
   result(cd_exitcode)

!  PURPOSE
!    This function basically conducts a small molecular dynamics simulation.
!    There are, however, subtle differences as e.g. energy does not have
!    to be conserved.
!
!     o points
!       These points can move freely and their movement is driven by
!       four different forces (no more than three at a time):
!
!       - "density force"
!         drives the (free) point along the density gradient to the desired
!         iso-density value (either free atoms superposition or real density)
!
!       - "gravity force"
!         if density force is absent (too far off, no gradient), this pulls
!         the (free) point towards a defined center {gravity_center}
!
!       - "repulsive force"
!         acts between free points, but only perpendicular to the density
!         gradient so that the forces don't disturb each other;
!         this leads to a more regular distribution of the points on the
!         cavity surface
!
!    The coordinates of all free points and the corresponding normal vectors
!    are synchronized between the MPI threads. Each thread will only
!    take a certain subset of points, calculate the forces acting on it,
!    propagate it and broadcast the new positions to all other threads.
!
!    This is just a wrapper around the actual serial and parallel functions.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), allocatable, intent(inout) :: points(:)
   real(dp), intent(in) :: gravity_center(3)

   integer :: cd_exitcode

!  INPUTS
!   o points -- sampling points' coordinates and normals
!   o gravity_center -- coordinates of imaginary gravitational center
!                       (for gravity force)
!  OUTPUT
!   o points -- updated coordinates and normals of sampling points
!   o cd_exitcode -- indicates the exit status (i.e. did not reach convergence)
!                    use exit codes from ISC_CONST
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! LOCAL VARIABLES

   ! TODO: actually fix parallel cd with kill ratio (currently fails with GNU
   !       compiler). The following line is just a workaround!
   if (n_tasks .eq. 1 .or. isc_kill_ratio .gt. 0.e0_dp) then
      cd_exitcode = constraint_dynamics_serial(points, gravity_center)
   else
      cd_exitcode = constraint_dynamics_parallel(points, gravity_center)
   endif

end function constraint_dynamics


!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/constraint_dynamics_parallel
!  NAME
!    constraint_dynamics_parallel
!  SYNOPSIS

function constraint_dynamics_parallel(points, gravity_center) &
   result(cd_exitcode)

!  PURPOSE
!    This function basically conducts a small molecular dynamics simulation.
!    There are, however, subtle differences as e.g. energy does not have
!    to be conserved.
!
!     o points
!       These points can move freely and their movement is driven by
!       four different forces (no more than three at a time):
!
!       - "density force"
!         drives the (free) point along the density gradient to the desired
!         iso-density value (either free atoms superposition or real density)
!
!       - "gravity force"
!         if density force is absent (too far off, no gradient), this pulls
!         the (free) point towards a defined center {gravity_center}
!
!       - "repulsive force"
!         acts between free points, but only perpendicular to the density
!         gradient so that the forces don't disturb each other;
!         this leads to a more regular distribution of the points on the
!         cavity surface
!
!    The coordinates of all free points and the corresponding normal vectors
!    are synchronized between the MPI threads. Each thread will only
!    take a certain subset of points, calculate the forces acting on it,
!    propagate it and broadcast the new positions to all other threads.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), allocatable, intent(inout) :: points(:)
   real(dp), intent(in) :: gravity_center(3)

   integer :: cd_exitcode

!  INPUTS
!   o points -- sampling points' coordinates and normals
!   o gravity_center -- coordinates of imaginary gravitational center
!                       (for gravity force)
!  OUTPUT
!   o points -- updated coordinates and normals of sampling points
!   o cd_exitcode -- indicates the exit status (i.e. did not reach convergence)
!                    use exit codes from ISC_CONST
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! LOCAL VARIABLES

   character(*), parameter :: func = 'constraint_dynamics_parallel'
   character(132) :: info_str

   integer :: n_points, n_keep, n_kill, i_kill, i_loc, i_dev, i_keep, i_proc

   ! working arrays
   ! number of beads for all processes
   integer, allocatable :: load(:), new_load(:), offset(:)
   integer :: ind_task
   real(dp) :: ratio_remaining_points

   type(InterfacePoint), allocatable :: points_keep(:)
   type(NeighborList), allocatable :: mynblists(:)

   ! bead position, momenta and forces
   real(dp), allocatable :: tmpx(:,:), p(:,:), f(:,:)

   ! density and gradient at each point
   real(dp), allocatable :: myrho(:), myrho_gradient(:,:)

   ! tune "density" force
   real(dp), allocatable :: rho_rel_deviation(:)
   integer, allocatable :: dev_indices(:), nb_indices(:)
   real(dp) :: rho_max_rel_deviation, rho_k_mod

   ! square of maximal force component
   real(dp) :: f_max_sq, f_max_sq_task

   ! variables for timings
   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm
   real(dp) :: cpu_all, clock_all
   real(dp) :: cpu_dens, clock_dens, cpu_dens_tot, clock_dens_tot
   real(dp) :: cpu_intradist, clock_intradist, cpu_intradist_tot, &
             clock_intradist_tot
   real(dp) :: cpu_inter, clock_inter, cpu_inter_tot, clock_inter_tot
   real(dp) :: cpu_sync, clock_sync, cpu_sync_tot, clock_sync_tot
   real(dp) :: cpu_prop, clock_prop, cpu_prop_tot, clock_prop_tot

   integer :: ind_work, ind_n
   integer :: counter ! count number of dynamics steps

   logical :: flag_update_neighbor_list
   logical :: record_cavity_creation


   ! initialize with "not converged"
   cd_exitcode = ISC_CONST%CDYN_EXIT_NOTCONV

   n_points = size(points)

   ! get timestamp
   call get_timestamps(cpu_all, clock_all)
   cpu_dens_tot = 0
   clock_dens_tot = 0
   cpu_intradist_tot = 0
   clock_intradist_tot = 0
   cpu_inter_tot = 0
   clock_inter_tot = 0
   cpu_sync_tot = 0
   clock_sync_tot = 0
   cpu_prop_tot = 0
   clock_prop_tot = 0


   write(info_str,'(4X,A)') 'Starting constraint dynamics simulation'
   call localorb_info(info_str)

   ! state which type of density is used
   if ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPStat .or. &
          isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPDyn ) then
      write(info_str,'(6X,A)') 'using the multipole expansion density'
      call localorb_info(info_str)
   elseif ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoFree ) then
      write(info_str,'(6X,A)') 'using the free atoms superposition density'
      call localorb_info(info_str)
   else
      ! raise error
      write(info_str,'(A)') '***Internal error: type of reference '//&
                            'density unknown'
      call localorb_info(info_str)
      call aims_stop('Internal error: type of reference density unknown', &
                           func)
   endif

   record_cavity_creation = isc_record_every_n_steps .gt. 0
   if ( record_cavity_creation .and. myid.eq.0 ) then
      ! open output file
      open(unit=23, file=trim(isc_record_cavity_creation_file), &
         status='unknown', action='write')
      ! write initial configuration
      call write_xyz(23, points, write_normals=.false.)
   endif


   ! PRE-ALLOCATION
   allocate(load(0:n_tasks-1))
   allocate(new_load(0:n_tasks-1))
   allocate(offset(0:n_tasks-1))

   allocate(rho_rel_deviation(n_points))
   allocate(dev_indices(n_points))

   ! calculate load distribution
   call distribute_loads(n_points, n_tasks, load, offset)
   new_load = load

   ! ALLOCATION
   allocate(tmpx(3,load(myid)))
   allocate(p(3,load(myid)))
   allocate(f(3,load(myid)))

   allocate(mynblists(load(myid)))

   allocate(myrho(load(myid)))
   allocate(myrho_gradient(3,load(myid)))


   ! BEAD DYNAMICS SIMULATION
   p = 0.e0_dp
   rho_k_mod = 1.e0_dp

   ! careful when changing the initial value of the counter, as some
   ! counter mod interval checks have to trigger for the first cycle
   counter = 0
   mainloop: do
      ! updated kill number
      n_kill = ceiling(isc_kill_ratio*n_points)

      ! clear forces
      f = 0.e0_dp
      f_max_sq_task = 0.e0_dp

      ! check if neighbor lists need to be updated this cycle
      flag_update_neighbor_list = mod(counter, isc_update_nlist_interval).eq.0

      ! get timestamp
      call get_timestamps(cpu_dens, clock_dens)

      ! get density and density gradient at bead positions
      ! we distribute the load over all processes by "slicing" the positions
      ! array

      ! distinguish between free atoms sum and real density
      call InterfacePoint_vector_extract( &
               points(offset(myid)+1:offset(myid)+load(myid)), &
               coords=tmpx )
      if ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPStat .or. &
             isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPDyn ) then
         call get_rho_multipole_and_gradient(load(myid), tmpx, &
                                             myrho, myrho_gradient)
      elseif ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoFree ) then
         call get_free_rho_and_gradient(load(myid), tmpx, &
                                             myrho, myrho_gradient)
      else
         ! raise error
         call aims_stop('***Internal error: type of density unknown***')
      endif

      ! get relative density differences locally
      rho_rel_deviation = 0.e0_dp
      rho_rel_deviation(offset(myid)+1:offset(myid)+load(myid)) = &
                     1.e0_dp - myrho(:)/isc_isodensity_value
      ! sync to global array
      call sync_vector(rho_rel_deviation, n_points)

      if ( flag_update_neighbor_list .and. (counter.gt.0) &
                                    .and. (n_kill.gt.0) ) then
         write(info_str,'(6X,A,I8,A)') 'killing ', n_kill, ' points'
         call localorb_info(info_str)
         ! remove the "worst" points
         ! this only makes sense when the neighbor list is updated afterwards

         ! sort global array
         do i_dev = 1, n_points
            dev_indices(i_dev) = i_dev
         enddo
         call dquicksort_indexlist(n_points, -abs(rho_rel_deviation), &
                                    dev_indices)

         ! exchange removed points with end on certain task
         do i_kill = 1, n_kill
            ! find process
            i_dev = dev_indices(i_kill)
            i_loc = i_dev
            i_proc = 0
            do while (i_loc.gt.load(i_proc))
               i_loc = i_loc - load(i_proc)
               i_proc = i_proc + 1
            enddo
            ! move all the local quantities, restore local continuity
            if ( i_proc.eq.myid ) then
               p(:,i_loc) = p(:,new_load(myid))
               f(:,i_loc) = f(:,new_load(myid))
               mynblists(i_loc) = mynblists(new_load(myid))
               myrho(i_loc) = myrho(new_load(myid))
               myrho_gradient(:,i_loc) = myrho_gradient(:,new_load(myid))
            endif
            ! move global quantities, restore local continuity II
            points(i_dev) = points(offset(i_proc)+new_load(i_proc))
            rho_rel_deviation(i_dev) = rho_rel_deviation(offset(i_proc)+new_load(i_proc))
            ! update temporary count
            new_load(i_proc) = new_load(i_proc) - 1
         enddo ! i_kill

         load(0) = new_load(0)
         do i_proc = 1, n_tasks - 1
            ! update counts
            i_dev = offset(i_proc) + 1
            offset(i_proc) = offset(i_proc-1) + load(i_proc-1)
            load(i_proc) = new_load(i_proc)
            ! move global quantities, restore global continuity
            do i_keep = offset(i_proc)+1, offset(i_proc)+load(i_proc)
               points(i_keep) = points(i_dev)
               rho_rel_deviation(i_keep) = rho_rel_deviation(i_dev)
               i_dev = i_dev + 1
            enddo
         enddo
         n_points = sum(load)
         rho_rel_deviation(n_points+1:) = 0
      endif

      ! convergence criterion
      rho_max_rel_deviation = maxval(abs(rho_rel_deviation))

      ! diagonal elements of the interaction matrix
      do ind_work = 1, load(myid), 1

         ! avoid stalling points at density cutoff
         if ( (dot_product(myrho_gradient(:,ind_work), &
                          myrho_gradient(:,ind_work)) &
                  .le.isc_gradient_threshold) .and. &
              (rho_rel_deviation(offset(myid)+ind_work).gt.0.8e0_dp) ) then
            ! get normalized direction
            points(offset(myid)+ind_work)%normal = normalized_d3(&
                        gravity_center-points(offset(myid)+ind_work)%coord)

            ! add "gravity"
            f(:,ind_work) = f(:,ind_work) + isc_g_k * &
                        points(offset(myid)+ind_work)%normal
         else
            ! get normalized direction
            points(offset(myid)+ind_work)%normal = &
                        normalized_d3(myrho_gradient(:,ind_work))

            ! add "density force"
            f(:,ind_work) = f(:,ind_work) + isc_rho_k * rho_k_mod * &
                        rho_rel_deviation(offset(myid)+ind_work) * &
                        points(offset(myid)+ind_work)%normal
         endif

      enddo ! ind_work

      ! accumulate time
      call get_times(cpu_dens, clock_dens, cpu_dens_tot, clock_dens_tot)


      ! SYNCHRONIZE points

      ! get timestamp
      call get_timestamps(cpu_sync, clock_sync)

      ! points array
      do ind_task = 0, n_tasks-1, 1
         call InterfacePoint_vector_mpi_bcast(&
                  points(1+offset(ind_task):load(ind_task)+offset(ind_task)), &
                  ind_task, mpi_comm_global )
      enddo ! ind_task

      ! accumulate time
      call get_times(cpu_sync, clock_sync, cpu_sync_tot, clock_sync_tot)


      ! RECORD cavity generation "movie"
      if ( record_cavity_creation .and. (myid.eq.0) .and. &
           ( mod(counter, isc_record_every_n_steps).eq.0 ) ) &
         call write_xyz(23, points(:n_points), write_normals=.false.)


      ! use deviation from desired density value as convergence criterion
      if ( rho_max_rel_deviation .le. isc_rho_rel_deviation_threshold ) then
         write(info_str,'(4X,A,1X,I5,1X,A)') &
               'Constraint dynamics converged after', counter, 'steps.'
         call localorb_info(info_str)
         write(info_str,'(4X,A,E10.3,A,F6.3,A)') &
               'Largest remaining deviation from iso-density value is ', &
               rho_max_rel_deviation*isc_isodensity_value/(bohr**3), &
               ' e/AA^3 (', rho_max_rel_deviation*100, ' %)'
         call localorb_info(info_str)

         ! determine appropriate exit code
         if ( counter .eq. 0 ) then
            cd_exitcode = ISC_CONST%CDYN_EXIT_NOSTEP
         else
            cd_exitcode = ISC_CONST%CDYN_EXIT_CONV
         endif
         ! leave the main loop
         exit mainloop
      endif

      ! also exit loop after too many iterations
      if ( counter .ge. isc_max_dyn_steps ) then
         write(info_str,'(4X,A,1X,I5,1X,A)') &
            'Constraint dynamics did NOT converge within', counter, 'steps.'
         call localorb_info(info_str)
         write(info_str,'(4X,A,E10.3,A,F6.3,A)') &
               'Largest remaining deviation from iso-density value was ', &
               rho_max_rel_deviation*isc_isodensity_value/(bohr**3), &
               ' e/AA^3 (', rho_max_rel_deviation*100, ' %)'
         call localorb_info(info_str)

         cd_exitcode = ISC_CONST%CDYN_EXIT_NOTCONV
         ! leave main loop
         exit mainloop
      endif


      ! COUNTING the cycles happens here, after the exits
      counter = counter + 1

      ! calculate modification factor for density force
      ! this counteracts the linear decrease of the "density" force when
      ! approaching the desired density value
      rho_k_mod = max(1.e0_dp, isc_rho_k_mod_constant/rho_max_rel_deviation)


      ! get timestamp
      call get_timestamps(cpu_intradist, clock_intradist)

      if ( flag_update_neighbor_list ) then
         write(info_str,'(6X,A,I6)') &
            'updating neighbor list in cycle no. ', counter
         call localorb_info(info_str)

!         ! in the very beginning, charlensq has already been set, so don't
!         ! call the update function in the first cycle
!         if (counter.gt.1) then
!            call synced_determine_new_distsq_threshold(mynblists, &
!                                                      distsq_threshold)
!         endif

         ! update neighbor list
         do ind_work = 1, load(myid), 1
            call update_neighbor_list_intra(points(:n_points), &
                  ind_work+offset(myid), distsq_threshold, mynblists(ind_work))
         enddo ! ind_work

         call synced_print_neighbor_list_summary(mynblists) !TODO: delete me
      else
         ! get distance information between free beads.
         do ind_work = 1, load(myid), 1
            call get_intra_distances(points(:n_points), ind_work+offset(myid), &
                                       mynblists(ind_work))
         enddo ! ind_work
      endif

      ! accumulate time
      call get_times(cpu_intradist, clock_intradist, &
                     cpu_intradist_tot, clock_intradist_tot)


      ! get timestamp
      call get_timestamps(cpu_inter, clock_inter)

      ! get forces between beads, off-diagonal elements in interaction matrix
      ! MPI parallelization
      do ind_work = 1, load(myid), 1

         ! 1. forces between the free beads

         ! distances are known (see call above)
         ! loop over neighbors
         do ind_n = 1, size(mynblists(ind_work)%nb), 1
            ! add up interaction force
            ! relvecs_intra points in direction of bead ind_work, pos. relvec
            f(:,ind_work) = f(:,ind_work) + calculate_interaction_force( &
                                    mynblists(ind_work)%nb(ind_n)%relvec, &
                                    mynblists(ind_work)%nb(ind_n)%distsq, &
                                    points(offset(myid)+ind_work)%normal, &
                                    isc_rep_k )
         enddo ! ind_n
      enddo ! ind_work

      ! accumulate time
      call get_times(cpu_inter, clock_inter, cpu_inter_tot, clock_inter_tot)


      ! PROPAGATE

      ! get timestamp
      call get_timestamps(cpu_prop, clock_prop)

      ! now all the forces are known so we can propagate the system
      call leapfrog_propagator_w_friction( &
            points(offset(myid)+1:offset(myid)+load(myid)), &
            p, f, isc_dt, isc_dynamics_friction)

      ! accumulate time
      call get_times(cpu_prop, clock_prop, cpu_prop_tot, clock_prop_tot)

   enddo mainloop


   ! restore convergence by throwing away other points
   if ( (cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOTCONV) &
                     .and. isc_try_restore_convergence) then
      write(info_str,'(6X,A)') 'trying to find converged cavity by '//&
                              'excluding sampling points...'
      call localorb_info(info_str)
      ! sort global array
      do i_dev = 1, n_points
         dev_indices(i_dev) = i_dev
      enddo
      call dquicksort_indexlist(n_points, abs(rho_rel_deviation), &
                                    dev_indices)

      !TODO DELETE ME
      if ( myid.eq.0 ) then
         open(file="deviations.dat",action="write",status="unknown",unit=222)
            write(222,'(A,X,A6,X,A6,X,A16)') "#", "i_dev", "dev_index", "abs_dev"
         do i_dev = 1, n_points
            write(222,'(I6,X,I6,X,E16.8)') i_dev, dev_indices(i_dev), abs(rho_rel_deviation(dev_indices(i_dev)))
         enddo
         close(unit=222)
      endif

      allocate(points_keep(n_points))
      n_keep = 1
      do while (abs(rho_rel_deviation(dev_indices(n_keep))).lt.&
                  isc_rho_rel_deviation_threshold)
         points_keep(n_keep) = points(dev_indices(n_keep))
         n_keep = n_keep + 1
         if (n_keep.gt.n_points) exit
      enddo
      n_keep = n_keep - 1

      ! is the remaining set satisfactory?
      ratio_remaining_points = n_keep/real(n_points,dp)
      if ( ratio_remaining_points .ge. 0.5e0_dp ) then
         cd_exitcode = ISC_CONST%CDYN_EXIT_CONV
         write(info_str,'(8X,A)') 'succesful!'
         call localorb_info(info_str)
      else
         write(info_str,'(8X,A)') 'failed (too many points'//&
                                 ' have been deleted)!'
         call localorb_info(info_str)
      endif
      write(info_str,'(6X,I8,A,F6.2,A)') n_keep, ' (', &
                        100*ratio_remaining_points, '%) points remaining'
      call localorb_info(info_str)

      ! reallocate points array of proper size
      n_points = n_keep
      deallocate(points)
      allocate(points(n_points))
      ! copy back
      do i_dev = 1, n_points
         points(i_dev) = points_keep(i_dev)
      enddo
      deallocate(points_keep)

      if ( record_cavity_creation .and. myid.eq.0 ) &
         call write_xyz(23, points(:n_points), write_normals=.false.)
   endif

   ! finally, proper reallocation if necessary
   if (size(points).ne.n_points) then
      allocate(points_keep(n_points))
      do i_keep = 1, n_points
         points_keep(i_keep) = points(i_keep)
      enddo
      deallocate(points)
      allocate(points(n_points))
      do i_keep = 1, n_points
         points(i_keep) = points_keep(i_keep)
      enddo
      deallocate(points_keep)
   endif

   ! update distsq_threshold at the end of the run for Voronoi construction
   ! unless the chose geometry has been optimal and not even a single full
   ! cycle has been computed (e.g. due to restart)
   if ( counter.gt.1 ) &
      call synced_determine_new_distsq_threshold(mynblists, distsq_threshold)


   ! close output
   if ( record_cavity_creation .and. myid.eq.0 ) then
      close(unit=23)
   endif


   ! DEALLOCATION

   deallocate(myrho_gradient)
   deallocate(myrho)

   deallocate(mynblists)

   deallocate(tmpx)
   deallocate(p)
   deallocate(f)

   deallocate(rho_rel_deviation)
   deallocate(dev_indices)

   deallocate(offset)
   deallocate(load)


   ! accumulate time
   call get_times(cpu_all, clock_all)
   call output_times('6X', 'Constraint dynamics run', &
             cpu_all, clock_all, defprio)
   call output_times('8X', 'Evaluation of density force', &
             cpu_dens_tot, clock_dens_tot, defprio)
   call output_times('8X', 'Distance calculation between beads', &
             cpu_intradist_tot, clock_intradist_tot, defprio)
   call output_times('8X', 'Evaluation of interaction force', &
             cpu_inter_tot, clock_inter_tot, defprio)
   call output_times('8X', 'Synchronization', &
             cpu_sync_tot, clock_sync_tot, defprio)
   call output_times('8X', 'Propagation', &
             cpu_prop_tot, clock_prop_tot, defprio)

end function constraint_dynamics_parallel


!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/constraint_dynamics_serial
!  NAME
!    constraint_dynamics_serial
!  SYNOPSIS

function constraint_dynamics_serial(points, gravity_center) &
   result(cd_exitcode)

!  PURPOSE
!    This function basically conducts a small molecular dynamics simulation.
!    There are, however, subtle differences as e.g. energy does not have
!    to be conserved.
!
!     o points
!       These points can move freely and their movement is driven by
!       four different forces (no more than three at a time):
!
!       - "density force"
!         drives the (free) point along the density gradient to the desired
!         iso-density value (either free atoms superposition or real density)
!
!       - "gravity force"
!         if density force is absent (too far off, no gradient), this pulls
!         the (free) point towards a defined center {gravity_center}
!
!       - "repulsive force"
!         acts between free points, but only perpendicular to the density
!         gradient so that the forces don't disturb each other;
!         this leads to a more regular distribution of the points on the
!         cavity surface
!
!    The coordinates of all free points and the corresponding normal vectors
!    are synchronized between the MPI threads. Each thread will only
!    take a certain subset of points, calculate the forces acting on it,
!    propagate it and broadcast the new positions to all other threads.
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), allocatable, intent(inout) :: points(:)
   real(dp), intent(in) :: gravity_center(3)

   integer :: cd_exitcode

!  INPUTS
!   o points -- sampling points' coordinates and normals
!   o gravity_center -- coordinates of imaginary gravitational center
!                       (for gravity force)
!  OUTPUT
!   o points -- updated coordinates and normals of sampling points
!   o cd_exitcode -- indicates the exit status (i.e. did not reach convergence)
!                    use exit codes from ISC_CONST
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! LOCAL VARIABLES

   character(*), parameter :: func = 'constraint_dynamics_serial'
   character(132) :: info_str

   integer :: n_points, n_keep, n_kill, i_kill, i_dev, i_keep

   ! working arrays
   ! number of beads for all processes
   real(dp) :: ratio_remaining_points

   type(InterfacePoint), allocatable :: points_keep(:)
   type(NeighborList), allocatable :: mynblists(:)

   ! bead position, momenta and forces
   real(dp), allocatable :: tmpx(:,:), p(:,:), f(:,:)

   ! density and gradient at each point
   real(dp), allocatable :: myrho(:), myrho_gradient(:,:)

   ! tune "density" force
   real(dp), allocatable :: rho_rel_deviation(:)
   integer, allocatable :: dev_indices(:), nb_indices(:)
   real(dp) :: rho_max_rel_deviation, rho_k_mod

   ! square of maximal force component
   real(dp) :: f_max_sq, f_max_sq_task

   ! variables for timings
   character(*), parameter :: deffmt = '4X'
   integer, parameter :: defprio = OL_norm
   real(dp) :: cpu_all, clock_all
   real(dp) :: cpu_dens, clock_dens, cpu_dens_tot, clock_dens_tot
   real(dp) :: cpu_intradist, clock_intradist, cpu_intradist_tot, &
             clock_intradist_tot
   real(dp) :: cpu_inter, clock_inter, cpu_inter_tot, clock_inter_tot
   real(dp) :: cpu_sync, clock_sync, cpu_sync_tot, clock_sync_tot
   real(dp) :: cpu_prop, clock_prop, cpu_prop_tot, clock_prop_tot

   integer :: ind_work, ind_n
   integer :: counter ! count number of dynamics steps

   logical :: flag_update_neighbor_list
   logical :: record_cavity_creation


   ! initialize with "not converged"
   cd_exitcode = ISC_CONST%CDYN_EXIT_NOTCONV

   n_points = size(points)

   ! get timestamp
   call get_timestamps(cpu_all, clock_all)
   cpu_dens_tot = 0
   clock_dens_tot = 0
   cpu_intradist_tot = 0
   clock_intradist_tot = 0
   cpu_inter_tot = 0
   clock_inter_tot = 0
   cpu_sync_tot = 0
   clock_sync_tot = 0
   cpu_prop_tot = 0
   clock_prop_tot = 0


   write(info_str,'(4X,A)') 'Starting constraint dynamics simulation'
   call localorb_info(info_str)

   ! state which type of density is used
   if ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPStat .or. &
          isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPDyn ) then
      write(info_str,'(6X,A)') 'using the multipole expansion density'
      call localorb_info(info_str)
   elseif ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoFree ) then
      write(info_str,'(6X,A)') 'using the free atoms superposition density'
      call localorb_info(info_str)
   else
      ! raise error
      write(info_str,'(A)') '***Internal error: type of reference '//&
                            'density unknown'
      call localorb_info(info_str)
      call aims_stop('Internal error: type of reference density unknown', &
                           func)
   endif

   record_cavity_creation = isc_record_every_n_steps .gt. 0
   if ( record_cavity_creation ) then
      ! open output file
      open(unit=23, file=trim(isc_record_cavity_creation_file), &
         status='unknown', action='write')
      ! write initial configuration
      call write_xyz(23, points, write_normals=.false.)
   endif

   ! ALLOCATION
   allocate(rho_rel_deviation(n_points))
   allocate(dev_indices(n_points))

   allocate(tmpx(3,n_points))
   allocate(p(3,n_points))
   allocate(f(3,n_points))

   allocate(mynblists(n_points))

   allocate(myrho(n_points))
   allocate(myrho_gradient(3,n_points))


   ! BEAD DYNAMICS SIMULATION
   p = 0.e0_dp
   rho_k_mod = 1.e0_dp

   ! careful when changing the initial value of the counter, as some
   ! counter mod interval checks have to trigger for the first cycle
   counter = 0
   mainloop: do
      ! updated kill number
      n_kill = ceiling(isc_kill_ratio*n_points)

      ! clear forces
      f = 0.e0_dp
      f_max_sq_task = 0.e0_dp

      ! check if neighbor lists need to be updated this cycle
      flag_update_neighbor_list = mod(counter, isc_update_nlist_interval).eq.0

      ! get timestamp
      call get_timestamps(cpu_dens, clock_dens)

      ! get density and density gradient at bead positions

      ! distinguish between free atoms sum and real density
      call InterfacePoint_vector_extract( &
               points(:n_points), &
               coords=tmpx )
      if ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPStat .or. &
             isc_cavity_type.eq.ISC_CONST%CAVITY_RhoMPDyn ) then
         call get_rho_multipole_and_gradient(n_points, tmpx, &
                                             myrho, myrho_gradient)
      elseif ( isc_cavity_type.eq.ISC_CONST%CAVITY_RhoFree ) then
         call get_free_rho_and_gradient(n_points, tmpx, &
                                             myrho, myrho_gradient)
      else
         ! raise error
         call aims_stop('***Internal error: type of density unknown***')
      endif

      ! get relative density differences locally
      rho_rel_deviation = 0.e0_dp
      rho_rel_deviation(:n_points) = &
                     1.e0_dp - myrho(:n_points)/isc_isodensity_value

      if ( flag_update_neighbor_list .and. (counter.gt.0) &
                                    .and. (n_kill.gt.0) ) then
         write(info_str,'(6X,A,I8,A)') 'killing ', n_kill, ' points'
         call localorb_info(info_str)
         ! remove the "worst" points
         ! this only makes sense when the neighbor list is updated afterwards

         ! sort global array
         do i_dev = 1, n_points
            dev_indices(i_dev) = i_dev
         enddo
         call dquicksort_indexlist(n_points, -abs(rho_rel_deviation), &
                                    dev_indices)

         ! exchange removed points with end
         do i_kill = 1, n_kill
            i_dev = dev_indices(i_kill)
            ! move all the quantities, restore continuity
            p(:,i_dev) = p(:,n_points)
            f(:,i_dev) = f(:,n_points)
            mynblists(i_dev) = mynblists(n_points)
            myrho(i_dev) = myrho(n_points)
            myrho_gradient(:,i_dev) = myrho_gradient(:,n_points)
            points(i_dev) = points(n_points)
            rho_rel_deviation(i_dev) = rho_rel_deviation(n_points)
            ! update temporary count
            n_points = n_points - 1
         enddo ! i_kill
         rho_rel_deviation(n_points+1:) = 0
      endif

      ! convergence criterion
      rho_max_rel_deviation = maxval(abs(rho_rel_deviation))

      ! diagonal elements of the interaction matrix
      do ind_work = 1, n_points, 1

         ! avoid stalling points at density cutoff
         if ( (dot_product(myrho_gradient(:,ind_work), &
                          myrho_gradient(:,ind_work)) &
                  .le.isc_gradient_threshold) .and. &
              (rho_rel_deviation(ind_work).gt.0.8e0_dp) ) then
            ! get normalized direction
            points(ind_work)%normal = normalized_d3(&
                        gravity_center-points(ind_work)%coord)

            ! add "gravity"
            f(:,ind_work) = f(:,ind_work) + isc_g_k * &
                        points(ind_work)%normal
         else
            ! get normalized direction
            points(ind_work)%normal = &
                        normalized_d3(myrho_gradient(:,ind_work))

            ! add "density force"
            f(:,ind_work) = f(:,ind_work) + isc_rho_k * rho_k_mod * &
                        rho_rel_deviation(ind_work) * &
                        points(ind_work)%normal
         endif

      enddo ! ind_work

      ! accumulate time
      call get_times(cpu_dens, clock_dens, cpu_dens_tot, clock_dens_tot)


      ! SYNCHRONIZE points

      ! get timestamp
      call get_timestamps(cpu_sync, clock_sync)

      ! accumulate time
      call get_times(cpu_sync, clock_sync, cpu_sync_tot, clock_sync_tot)


      ! RECORD cavity generation "movie"
      if ( record_cavity_creation .and. &
           ( mod(counter, isc_record_every_n_steps).eq.0 ) ) &
         call write_xyz(23, points(:n_points), write_normals=.false.)


      ! use deviation from desired density value as convergence criterion
      if ( rho_max_rel_deviation .le. isc_rho_rel_deviation_threshold ) then
         write(info_str,'(4X,A,1X,I5,1X,A)') &
               'Constraint dynamics converged after', counter, 'steps.'
         call localorb_info(info_str)
         write(info_str,'(4X,A,E10.3,A,F6.3,A)') &
               'Largest remaining deviation from iso-density value is ', &
               rho_max_rel_deviation*isc_isodensity_value/(bohr**3), &
               ' e/AA^3 (', rho_max_rel_deviation*100, ' %)'
         call localorb_info(info_str)

         ! determine appropriate exit code
         if ( counter .eq. 0 ) then
            cd_exitcode = ISC_CONST%CDYN_EXIT_NOSTEP
         else
            cd_exitcode = ISC_CONST%CDYN_EXIT_CONV
         endif
         ! leave the main loop
         exit mainloop
      endif

      ! also exit loop after too many iterations
      if ( counter .ge. isc_max_dyn_steps ) then
         write(info_str,'(4X,A,1X,I5,1X,A)') &
            'Constraint dynamics did NOT converge within', counter, 'steps.'
         call localorb_info(info_str)
         write(info_str,'(4X,A,E10.3,A,F6.3,A)') &
               'Largest remaining deviation from iso-density value was ', &
               rho_max_rel_deviation*isc_isodensity_value/(bohr**3), &
               ' e/AA^3 (', rho_max_rel_deviation*100, ' %)'
         call localorb_info(info_str)

         cd_exitcode = ISC_CONST%CDYN_EXIT_NOTCONV
         ! leave main loop
         exit mainloop
      endif


      ! COUNTING the cycles happens here, after the exits
      counter = counter + 1

      ! calculate modification factor for density force
      ! this counteracts the linear decrease of the "density" force when
      ! approaching the desired density value
      rho_k_mod = max(1.e0_dp, isc_rho_k_mod_constant/rho_max_rel_deviation)


      ! get timestamp
      call get_timestamps(cpu_intradist, clock_intradist)

      if ( flag_update_neighbor_list ) then
         write(info_str,'(6X,A,I6)') &
            'updating neighbor list in cycle no. ', counter
         call localorb_info(info_str)

!         ! in the very beginning, charlensq has already been set, so don't
!         ! call the update function in the first cycle
!         if (counter.gt.1) then
!            call synced_determine_new_distsq_threshold(mynblists, &
!                                                      distsq_threshold)
!         endif

         ! update neighbor list
         do ind_work = 1, n_points, 1
            call update_neighbor_list_intra(points(:n_points), &
                  ind_work, distsq_threshold, mynblists(ind_work))
         enddo ! ind_work

         call synced_print_neighbor_list_summary(mynblists) !TODO: delete me
      else
         ! get distance information between free beads.
         do ind_work = 1, n_points, 1
            call get_intra_distances(points(:n_points), ind_work, &
                                       mynblists(ind_work))
         enddo ! ind_work
      endif

      ! accumulate time
      call get_times(cpu_intradist, clock_intradist, &
                     cpu_intradist_tot, clock_intradist_tot)


      ! get timestamp
      call get_timestamps(cpu_inter, clock_inter)

      ! get forces between beads, off-diagonal elements in interaction matrix
      do ind_work = 1, n_points, 1

         ! 1. forces between the free beads

         ! distances are known (see call above)
         ! loop over neighbors
         do ind_n = 1, size(mynblists(ind_work)%nb), 1
            ! add up interaction force
            ! relvecs_intra points in direction of bead ind_work, pos. relvec
            f(:,ind_work) = f(:,ind_work) + calculate_interaction_force( &
                                    mynblists(ind_work)%nb(ind_n)%relvec, &
                                    mynblists(ind_work)%nb(ind_n)%distsq, &
                                    points(ind_work)%normal, &
                                    isc_rep_k )
         enddo ! ind_n
      enddo ! ind_work

      ! accumulate time
      call get_times(cpu_inter, clock_inter, cpu_inter_tot, clock_inter_tot)


      ! PROPAGATE

      ! get timestamp
      call get_timestamps(cpu_prop, clock_prop)

      ! now all the forces are known so we can propagate the system
      call leapfrog_propagator_w_friction( &
            points(:n_points), &
            p, f, isc_dt, isc_dynamics_friction)

      ! accumulate time
      call get_times(cpu_prop, clock_prop, cpu_prop_tot, clock_prop_tot)

   enddo mainloop


   ! restore convergence by throwing away other points
   if ( (cd_exitcode.eq.ISC_CONST%CDYN_EXIT_NOTCONV) &
                     .and. isc_try_restore_convergence) then
      write(info_str,'(6X,A)') 'trying to find converged cavity by '//&
                              'excluding sampling points...'
      call localorb_info(info_str)
      ! sort global array
      do i_dev = 1, n_points
         dev_indices(i_dev) = i_dev
      enddo
      call dquicksort_indexlist(n_points, abs(rho_rel_deviation), &
                                    dev_indices)

      !TODO DELETE ME
      open(file="deviations.dat",action="write",status="unknown",unit=222)
         write(222,'(A,X,A6,X,A6,X,A16)') "#", "i_dev", "dev_index", "abs_dev"
      do i_dev = 1, n_points
         write(222,'(I6,X,I6,X,E16.8)') i_dev, dev_indices(i_dev), abs(rho_rel_deviation(dev_indices(i_dev)))
      enddo
      close(unit=222)

      allocate(points_keep(n_points))
      n_keep = 1
      do while (abs(rho_rel_deviation(dev_indices(n_keep))).lt.&
                  isc_rho_rel_deviation_threshold)
         points_keep(n_keep) = points(dev_indices(n_keep))
         n_keep = n_keep + 1
         if (n_keep.gt.n_points) exit
      enddo
      n_keep = n_keep - 1

      ! is the remaining set satisfactory?
      ratio_remaining_points = n_keep/real(n_points,dp)
      if ( ratio_remaining_points .ge. 0.5e0_dp ) then
         cd_exitcode = ISC_CONST%CDYN_EXIT_CONV
         write(info_str,'(8X,A)') 'succesful!'
         call localorb_info(info_str)
      else
         write(info_str,'(8X,A)') 'failed (too many points'//&
                                 ' have been deleted)!'
         call localorb_info(info_str)
      endif
      write(info_str,'(6X,I8,A,F6.2,A)') n_keep, ' (', &
                        100*ratio_remaining_points, '%) points remaining'
      call localorb_info(info_str)

      ! reallocate points array of proper size
      n_points = n_keep
      deallocate(points)
      allocate(points(n_points))
      ! copy back
      do i_dev = 1, n_points
         points(i_dev) = points_keep(i_dev)
      enddo
      deallocate(points_keep)

      if ( record_cavity_creation ) &
         call write_xyz(23, points(:n_points), write_normals=.false.)
   endif

   ! finally, proper reallocation if necessary
   if (size(points).ne.n_points) then
      allocate(points_keep(n_points))
      do i_keep = 1, n_points
         points_keep(i_keep) = points(i_keep)
      enddo
      deallocate(points)
      allocate(points(n_points))
      do i_keep = 1, n_points
         points(i_keep) = points_keep(i_keep)
      enddo
      deallocate(points_keep)
   endif

   ! update distsq_threshold at the end of the run for Voronoi construction
   ! unless the chose geometry has been optimal and not even a single full
   ! cycle has been computed (e.g. due to restart)
   if ( counter.gt.1 ) &
      call synced_determine_new_distsq_threshold(mynblists, distsq_threshold)


   ! close output
   if ( record_cavity_creation ) then
      close(unit=23)
   endif


   ! DEALLOCATION

   deallocate(myrho_gradient)
   deallocate(myrho)

   deallocate(mynblists)

   deallocate(tmpx)
   deallocate(p)
   deallocate(f)

   deallocate(rho_rel_deviation)
   deallocate(dev_indices)


   ! accumulate time
   call get_times(cpu_all, clock_all)
   call output_times('6X', 'Constraint dynamics run', &
             cpu_all, clock_all, defprio)
   call output_times('8X', 'Evaluation of density force', &
             cpu_dens_tot, clock_dens_tot, defprio)
   call output_times('8X', 'Distance calculation between beads', &
             cpu_intradist_tot, clock_intradist_tot, defprio)
   call output_times('8X', 'Evaluation of interaction force', &
             cpu_inter_tot, clock_inter_tot, defprio)
   call output_times('8X', 'Synchronization', &
             cpu_sync_tot, clock_sync_tot, defprio)
   call output_times('8X', 'Propagation', &
             cpu_prop_tot, clock_prop_tot, defprio)

end function constraint_dynamics_serial


pure function calculate_interaction_force(relvec, distsq, normal, fconstant)

   implicit none

   real(dp), dimension(3), intent(in) :: relvec
   real(dp), intent(in) :: distsq
   real(dp), dimension(3), intent(in) :: normal
   real(dp), intent(in) :: fconstant

   real(dp), dimension(3) :: calculate_interaction_force

   real(dp), dimension(3) :: force_dir
   real(dp) :: normsq

   ! get force direction
   ! we assume that relvec points in direction of the current bead, so we
   ! just project out the normal component
   force_dir = relvec - dot_d3(relvec, normal) * normal

   normsq = dot_d3(force_dir,force_dir)
   if (normsq .eq. 0.e0_dp) then
      ! force is zero
      calculate_interaction_force = 0.e0_dp
   else
      ! normalize force
      force_dir = force_dir / sqrt( normsq )
      calculate_interaction_force = fconstant * force_dir * charlensq/distsq
   endif

end function calculate_interaction_force


subroutine leapfrog_propagator_w_friction(points, p, f, dt, fric)

   implicit none

   type(InterfacePoint), intent(inout) :: points(:)
   real(dp), intent(inout) :: p(:,:)  ! momenta
   real(dp), intent(in) :: f(:,:)     ! forces
   real(dp), intent(in) :: dt         ! time step
   real(dp), intent(in) :: fric       ! friction factor

   integer :: i_point

   ! we enter with: x=x(t), p=p(t-dt/2), f=f(t)
   !! p(t+dt/2) = p(t-dt/2) + dt*f(t)
   p = ( p + dt * f ) * ( 1.e0_dp - fric )
   ! now: x=x(t), p=p(t+dt/2), f=f(t)
   !! x(t+dt) = x(t) + p(t+dt/2) * dt
   do i_point = lbound(points,1), ubound(points,1), 1
      points(i_point)%coord = points(i_point)%coord + &
                     p(:,i_point) * dt ! / points(i_point)%mass ! when masses should be considered
   enddo
   ! we exit with: x=x(t+dt), p=p(t+dt/2), f=f(t)

end subroutine leapfrog_propagator_w_friction


subroutine write_xyz(use_unit, points, write_normals)

   implicit none

   integer, intent(in) :: use_unit
   type(InterfacePoint), intent(in) :: points(:)
   logical, intent(in) :: write_normals

   integer :: i_point, n_points

   n_points = size(points)
   write(use_unit,'(I8)') n_points
   write(use_unit,*) ''

   if (write_normals) then
      do i_point = 1, n_points, 1
         write(use_unit,'(A2,6(2X,E16.6))') "X", &
                  points(i_point)%coord*bohr, points(i_point)%normal
      enddo ! i_point
   else
      do i_point = 1, n_points, 1
         write(use_unit,'(A2,3(2X,E16.6))') "X", &
                  points(i_point)%coord*bohr
      enddo ! i_point
   endif

end subroutine write_xyz



!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/synced_print_neighbor_list_summary
!  NAME
!    synced_print_neighbor_list_summary
!  SYNOPSIS

subroutine synced_print_neighbor_list_summary(nlists)

!  PURPOSE
!    This subroutine prints a very brief resumee of the number of neighbors
!    for all points (including necessary MPI synchronization)
!
!  USES
   implicit none

!  ARGUMENTS
   type(NeighborList), intent(in) :: nlists(:)

!  INPUTS
!   o nlists -- neighbor lists
!  OUTPUT
!    writes to standard output
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   character(*), parameter :: func = 'synced_print_neighbor_list_summary'
   character(132) :: info_str

   integer :: n_neighbors, i_list, max_neighbors, min_neighbors, min, max

   max = 0
   min = huge(n_neighbors) ! some very big number
   do i_list = lbound(nlists,1), ubound(nlists,1), 1
      n_neighbors = size(nlists(i_list)%nb)
      if (n_neighbors .gt. max) &
         max = n_neighbors
      if (n_neighbors .lt. min) &
         min = n_neighbors
   enddo ! i_list

   ! sync compare min and max with all other tasks
   min = -min ! invert to maximization problem
   call sync_find_max(max, max_neighbors)
   call sync_find_max(min, min_neighbors)
   min_neighbors = -min_neighbors

   write(info_str,'(8X,A,1X,I4,A,I4)') &
      'neighbor count statistics: min', min_neighbors, ', max', max_neighbors
   call localorb_info(info_str)

end subroutine synced_print_neighbor_list_summary


!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/synced_determine_new_distsq_threshold
!  NAME
!    synced_determine_new_distsq_threshold
!  SYNOPSIS

subroutine synced_determine_new_distsq_threshold(nlists, distsq_threshold)

!  PURPOSE
!    This subroutine determines a new, appropriate distance threshold
!    for considering two points neighboring.
!    This update is important if the cavity's shape changes significantly
!    during the optimization process in order to keep the number of neighbors
!    approximately constant.
!    Important: the distance threshold is synchronized over all MPI tasks
!    before returned
!
!  USES
   implicit none

!  ARGUMENTS
   type(NeighborList), dimension(:), intent(in) :: nlists

   real(dp), intent(out) :: distsq_threshold

!  INPUTS
!   o nlists -- neighbor lists
!  OUTPUT
!   o distsq_threshold -- squared distance threshold
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE

   ! do not allow more than the nearest 2*CAP_MEDIAN neighbors to be included
   integer, parameter :: CAP_MEDIAN = 20

   real(dp), parameter :: THRESHOLD_FACTOR = 1.7e0_dp

   integer :: n_lists, n_neighbors
   integer :: ind_j, ind_i, ind_m, max_neighbors
   integer, dimension(:), allocatable :: indices
   real(dp) :: max_median

   ! determine new characteristic length
   ! we define it as the maximal median of the distances found between
   ! all points and their old neighbors
   max_median = 0.e0_dp

   n_lists = size(nlists)

   max_neighbors = 0
   do ind_j = 1, n_lists, 1
      max_neighbors = max(max_neighbors, size(nlists(ind_j)%nb))
   enddo ! ind_j
   allocate(indices(max_neighbors))

   ! loop over all entries in neighbor list
   do ind_j = 1, n_lists, 1

      ! check if neighbor list is empty, should not happen!
      n_neighbors = size(nlists(ind_j)%nb)
      if (n_neighbors .eq. 0) then
         call aims_stop('neighbor list is empty, cannot determine'//&
                        ' new characteristic length')
      endif

      ! fill with standard values
      do ind_i = 1, n_neighbors, 1
         indices(ind_i) = ind_i
      enddo
      ! sort distances
      call dquicksort_indexlist(n_neighbors, nlists(ind_j)%nb(:)%distsq, &
                                 indices)
      ! get maximal median
      ind_m = indices( min( (n_neighbors+1)/2, CAP_MEDIAN ) )
      if ( nlists(ind_j)%nb(ind_m)%distsq .gt. max_median ) then
         max_median = nlists(ind_j)%nb(ind_m)%distsq
      endif

   enddo ! ind_j

   deallocate(indices)

   ! multiply by a threshold factor &
   ! SYNCHRONIZATION
   call get_max_double(distsq_threshold, max_median*THRESHOLD_FACTOR)

end subroutine synced_determine_new_distsq_threshold

!******
!-------------------------------------------------------------------------------
!****s* isc_constraint_dynamics/get_intra_distances
!  NAME
!    get_intra_distances
!  SYNOPSIS

subroutine get_intra_distances(points, i_point, nlist)

!  PURPOSE
!    This subroutine calculates the cartesian distance vectors and
!    squared cartesian distances (contained in {nlist}) between a certain
!    point with index i_point and all its neighboring points;
!    the latter information has to be provided in {nlist} as an input
!
!  USES
   implicit none

!  ARGUMENTS
   type(InterfacePoint), dimension(:), intent(in) :: points
   integer, intent(in) :: i_point

   type(NeighborList), intent(inout) :: nlist

!  INPUTS
!   o points -- points with assigned coordinate systems
!   o i_point -- index of point whose neighbor distances are calculated
!   o nlist -- neighbor list with indices of all neighboring points
!  OUTPUT
!   o nlist -- (updated) distance information to all neighboring points
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Development version, FHI-aims (2015).
!  SOURCE

   integer :: ind_n, ind_i

   do ind_n = 1, size(nlist%nb), 1
      ! auxiliary index variable
      ind_i = nlist%nb(ind_n)%ind

      ! distance calculation
      nlist%nb(ind_n)%relvec = points(i_point)%coord - points(ind_i)%coord
      nlist%nb(ind_n)%distsq = scale_interaction_distance( &
            dot_d3(nlist%nb(ind_n)%relvec, nlist%nb(ind_n)%relvec), &
            points(i_point)%normal, points(ind_i)%normal )
   enddo ! ind_n

end subroutine get_intra_distances



end module isc_constraint_dynamics

