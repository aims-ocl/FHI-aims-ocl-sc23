!****h* FHI-aims/mpi_utilities
!  NAME
!    mpi_utilities - provides the utilities for the parallel (MPI) environment
!  SYNOPSIS
      module mpi_utilities
!  PURPOSE
!    This module takes care of
!    * MPI_INIT and MPI_FINALIZE calls
!    * Host and processor accounting
!    * Task distribution
!  USES

      ! WPH: I have kept these global-level use statements here deliberately.
      !      I am quite comfortable with having this module be a "wrapper"
      !      module around the MPI infrastructure in aims, as long as this
      !      is relegated to ONLY modules strictly related to MPI.
      use mpi_tasks

!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE
!



!  Any other comments (not for documentation) follow below.

!  Declare global variables:

      implicit none


      integer, dimension(:), allocatable :: task_list
      integer, dimension(:,:), allocatable :: radial_task_list
      integer, dimension(:), allocatable :: batch_task_list

      integer, dimension(:), allocatable :: spline_atom_storage



      contains
!******
!------------------------------------------------------------------------------
!****s* mpi_utilities/initialize_mpi
!  NAME
!    initialize_mpi - initializes the MPI environment
!  SYNOPSIS
      subroutine initialize_mpi ()
!  PURPOSE
!    This subroutine:
!    * initializes the MPI environment,
!    * finds out the id of the thread, 
!    * finds out the total number of parallel tasks
!    If this is a serial run (linked to mpi_stubs instead of a real
!    MPI library), this information is passed on to the rest of the
!    program by setting use_mpi .false.
!  USES
      use localorb_io, only: localorb_info, OL_norm, use_unit
      implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE

!  Declare local variables:

      integer :: mpierr

      call MPI_INIT(mpierr)

      if (mpierr .eq. -1) then
         use_mpi = .false.
         n_tasks = 0

      else if (mpierr.ne.MPI_SUCCESS) then
         write(use_unit,*) "Error in MPI initialization."
         write(use_unit,*) "Exiting..."
         stop
      else
         use_mpi = .true.
      end if

      mpi_comm_global = MPI_COMM_WORLD

      ! This also gets the total number of tasks ...
      call get_my_task()

      end subroutine initialize_mpi
!******
!------------------------------------------------------------------------------
!****s* mpi_utilities/initial_mpi_report
!  NAME
!   initial_mpi_report
!  SYNOPSIS

      subroutine initial_mpi_report

!  PURPOSE
!    Gives out the initial MPI-status report
!  USES
      use localorb_io, only: localorb_allinfo, localorb_info, OL_norm, use_unit
      implicit none
!  ARGUMENTS 

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      character(LEN=MPI_MAX_PROCESSOR_NAME+100) :: info_str
      character(LEN=MPI_MAX_PROCESSOR_NAME) :: my_proc_name, tmp_name
      character(LEN=MPI_MAX_PROCESSOR_NAME) :: adjusted_proc_name
      integer :: my_proc_name_length

      if (use_mpi) then

        call get_my_processor( tmp_name, my_proc_name_length )
        write(my_proc_name, "(A)") adjustl(tmp_name(1:my_proc_name_length))

        ! Only task number zero reports number of tasks ...
        write(info_str,'(2X,A,I8,1X,A)') "Using ", n_tasks, &
          "parallel tasks."
        call localorb_info( info_str ,use_unit,'(A)',OL_norm)

        adjusted_proc_name = adjustl(my_proc_name)

        ! ... but this is the one and only time that each task reports to work.
        write(info_str,'(2X,A,I8,1X,A,A,1X,A,1X,A,I8,1X,A,I8,1X,A,I8,1X,A,I8)') "Task ", myid, & 
           "on host ", trim(adjusted_proc_name), "reporting.", &
           "shm_rank ", shm_rank, &
           ", shm_num_of_rank ", shm_num_of_rank, &
           ", shm_per0_rank ", shm_per0_rank, &
           ", shm_per0_num_of_rank ", shm_per0_num_of_rank
        call localorb_allinfo(info_str, use_unit, '(A)', OL_norm)

      end if

      end subroutine initial_mpi_report
!****** 
!------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_tasks
!  NAME
!    distribute_tasks
!  SYNOPSIS
      subroutine distribute_tasks(method)
!  PURPOSE
!    Distributes the work over atoms to different tasks.
!  USES
     use dimensions, only: n_atoms
     use grids, only: n_radial, n_angular
     use geometry, only: species

!  ARGUMENTS 
      integer :: method
!  INPUTS
!    o method -- the method for the parallel task distribution, zero gives a
!                round-robin distribution and one leads to a work-balanced distribution
!  OUTPUT
!    none -- at exit the array task_list is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!     local variables
      integer :: i_atom
      integer :: i_radial
      integer, dimension(:), allocatable :: work_list
      integer, dimension(:), allocatable :: task_work_table

      integer :: index
      integer :: index_2

      if (.not.allocated(task_list)) then
         allocate(task_list(n_atoms))
      end if

      if (.not.allocated(work_list)) then
         allocate(work_list(n_atoms))
      end if

      if (.not.allocated(task_work_table)) then
         allocate(task_work_table(n_tasks))
      end if

      select case (method)

      case(0)
         do i_atom = 1, n_atoms, 1
            task_list(i_atom) = MOD(i_atom, n_tasks)
         end do

      case(1)
         work_list = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               work_list(i_atom) = work_list(i_atom) + &
                    n_angular(i_radial, species(i_atom))
            enddo
         enddo

         task_list = 0
         task_work_table = 0

         do i_atom = 1, n_atoms, 1
            index = maxloc(work_list,1)
            index_2 = minloc(task_work_table,1)
            task_list(index) = index_2 - 1
            task_work_table(index_2) = task_work_table(index_2) &
                 + work_list(index)
            work_list(index) = 0
         enddo

      case default
         task_list = 0

      end select




      if (allocated(task_work_table)) then
         deallocate(task_work_table)
      end if
      if (allocated(work_list)) then
         deallocate(work_list)
      end if

      end subroutine distribute_tasks
!****** 

!----------------------------------------------------------
!****s* mpi_utilities/distribute_radial_tasks
!  NAME
!    distribute_radial_tasks
!  SYNOPSIS
      subroutine distribute_radial_tasks(method)
!  PURPOSE
!    Distributes the work over radial shells of the grid to different tasks.
!    This subroutine distributes the work over radial shells to different tasks
!    This is needed in initialize_integrals
!  USES
      use dimensions, only: n_atoms, n_full_points, n_max_radial
      use runtime_choices
      use localorb_io
      use grids, only: n_radial, n_angular
      use geometry, only: species
!
!  ARGUMENTS 
      integer :: method
!  INPUTS
!    o method -- the method for the parallel task distribution, minus one and zero give a
!                round-robin distribution and one and two lead to a work-balanced distributions
!  OUTPUT
!    none -- at exit the array radial_task_list is set (method=2 sets also task_list)
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     local variables
      integer :: i_atom
      integer :: i_radial
      integer :: index, index_2

      integer, dimension(:), allocatable :: radial_task_work_table

!     only for the atomic task distribution method
      integer, dimension(:), allocatable :: work_list
      integer, dimension(:), allocatable :: task_work_table
      character*100 :: info_str

!     begin work


      if (.not.allocated(radial_task_list)) then
         allocate(radial_task_list(n_max_radial,n_atoms))
      end if

      if (.not.allocated(radial_task_work_table)) then
         allocate(radial_task_work_table(n_tasks))
      end if

      select case (method)

!     a round-robin distribution with no output, used in the initialization stage
      case(-1)
         radial_task_list = 0
         index = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               index = index + 1
               radial_task_list(i_radial, i_atom) &
                    = MOD(index, n_tasks)
            enddo
         enddo

!     a simple round-robin distribution
      case(0)
         radial_task_list = 0
         radial_task_work_table = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               index = MOD(i_radial, n_tasks) + 1
               radial_task_list(i_radial, i_atom) = index - 1
               radial_task_work_table(index) = &
                    radial_task_work_table(index) + &
                    n_angular(i_radial,species(i_atom))
            enddo
         enddo

!     a distribution that tries to balance the load over threads by computing
!     the number of integration points assigned to each thread
      case(1)
         radial_task_list = 0
         radial_task_work_table = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               index = minloc(radial_task_work_table,1)
               radial_task_list(i_radial,i_atom) = index - 1
               radial_task_work_table(index) = &
                    radial_task_work_table(index) + &
                    n_angular(i_radial,species(i_atom))
            enddo
         enddo

!     An atomic balacing scheme, where the workload is first balanced over atoms and
!     then this information is extended to radial shells. This scheme sets also the
!     task_list variable for atoms, so that there is no need to reset that.
      case(2)

         if (.not.allocated(task_list)) then
            allocate(task_list(n_atoms))
         end if

         if (.not.allocated(work_list)) then
            allocate(work_list(n_atoms))
         end if

         if (.not.allocated(task_work_table)) then
            allocate(task_work_table(n_tasks))
         end if

!     first, balance the atomic tasks
         work_list = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               work_list(i_atom) = work_list(i_atom) + &
                    n_angular(i_radial, species(i_atom))
            enddo
         enddo

         task_list = 0
         task_work_table = 0

         do i_atom = 1, n_atoms, 1
            index = maxloc(work_list,1)
            index_2 = minloc(task_work_table,1)
            task_list(index) = index_2 - 1
            task_work_table(index_2) = task_work_table(index_2) &
                 + work_list(index)
            work_list(index) = 0
         enddo

!     then, create the radial tasks list based on the atomic distribution
         radial_task_list = 0
         radial_task_work_table = 0

         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               radial_task_list(i_radial,i_atom) = task_list(i_atom)
               radial_task_work_table(task_list(i_atom)+1) = &
                    radial_task_work_table(task_list(i_atom)+1) + &
                    n_angular(i_radial,species(i_atom))
            enddo
         enddo

!     the default case is that the master does it all
      case default
         radial_task_list = 0

      end select

      if ((use_mpi).and.(method.ne.-1)) then
!         if (myid.eq.0) then
            write(info_str,'(A)')''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A,A,I5,A)') &
            "Integration load balanced ", &
            "across ", n_tasks, " MPI tasks. "
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') &
                 "Work distribution over tasks is as follows:"
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            if (output_priority .le. OL_low) then 
               do index = 1, n_tasks, 1
                  write(info_str,'(2X,A,I5,A,I10,A)') &
                       "Task ", index-1, &
                       " has ", &
                       radial_task_work_table(index), &
                       " integration points."
                  call localorb_info(info_str,use_unit,'(A)',OL_low)
               enddo
            end if
      end if

      n_full_points = radial_task_work_table(myid+1)

      if (allocated(radial_task_work_table)) then
         deallocate(radial_task_work_table)
      end if
      if (allocated(work_list)) then
         deallocate(work_list)
      end if
      if (allocated(task_work_table)) then
         deallocate(task_work_table)
      end if

! FIXME - temporary error trap!!
      if ((method.ne.-1) .and. (n_full_points.eq.0)) then
        write(use_unit,'(1X,A,I5,A)') "* Thread number ", myid, &
        " was not assigned any integration points."
        write(use_unit,'(1X,A,A)') "* While this is not ", &
        "strictly impossible, the present code version"
        write(use_unit,'(1X,A,A)') "* is not set up to deal with", &
        " this choice. I will stop here, please modify "
        write(use_unit,'(1X,A,A)') "* the code (disable zero ", &
        "allocations etc.) before continuing."
        call aims_stop("No inegration points","distribute_radial_tasks")
      end if

      end subroutine distribute_radial_tasks
!****** 
!------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_batch_tasks
!  NAME
!    distribute_batch_tasks
!  SYNOPSIS
      subroutine distribute_batch_tasks( batch_sizes )
!  PURPOSE
!    Distributes the grid batches to different tasks by balancing the workload.
!  USES
      use dimensions, only: n_grid_batches, n_full_points, n_my_batches
      use runtime_choices
      use localorb_io
      implicit none
!
!  ARGUMENTS
      integer :: batch_sizes( n_grid_batches )
!  INPUTS
!    o batch_sizes -- number of grid points in each of the grid batches
!  OUTPUT
!    none -- at exit the array batch_task_list is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: i_batch, index
      integer :: batch_task_work_table(n_tasks)
      character*100 :: info_str

      if (.not.allocated(batch_task_list)) then
         allocate(batch_task_list(n_grid_batches))
      end if

      batch_task_list = 0
      n_full_points = 0

      batch_task_work_table = 0

      n_my_batches = 0

      do i_batch = 1, n_grid_batches, 1

         index = MINLOC(batch_task_work_table,1)
         batch_task_list(i_batch) = index - 1
         batch_task_work_table(index) = &
              batch_task_work_table(index) + batch_sizes(i_batch)

         if (myid.eq.batch_task_list(i_batch)) then
            n_my_batches = n_my_batches + 1
            n_full_points = n_full_points + batch_sizes(i_batch)
         end if

      end do

      if (use_mpi) then
!         if (myid.eq.0) then
            write(info_str,'(A)') ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A,A,I5,A)') &
            "Integration load balanced ", &
            "across ", n_tasks, " MPI tasks. "
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') &
                 "Work distribution over tasks is as follows:"
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            if (output_priority .le. OL_low) then 
               do index = 1, n_tasks, 1
                  write(info_str,'(2X,A,I5,A,I10,A)') &
                       "Task ", index-1, &
                       " has ", &
                       batch_task_work_table(index), &
                       " integration points."
                  call localorb_info(info_str,use_unit,'(A)',OL_low)
               end do
            end if
!         end if
      end if

      end subroutine distribute_batch_tasks
!****** 
!------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_batch_tasks_by_location
!  NAME
!    distribute_batch_tasks_by_location
!  SYNOPSIS
      subroutine distribute_batch_tasks_by_location(batch_sizes,coords)
!  PURPOSE
!    Distributes the grid batches to different tasks based on their location.
!  USES
      use dimensions, only: n_grid_batches, n_full_points, n_my_batches
      use runtime_choices
      use localorb_io
      use synchronize_mpi_basic, only: sync_integer_vector
!
!  ARGUMENTS 
      integer :: batch_sizes( n_grid_batches )
      real*8  :: coords(3,n_grid_batches)
!  INPUTS
!    o batch_sizes -- number of grid points in each of the grid batches
!    o coords -- coordinates of the grid batches
!  OUTPUT
!    none -- at exit the array batch_task_list is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE


      real*8, allocatable :: batch_desc(:,:)
      integer :: i_batch, my_off, my_len
      character*100 :: info_str

      if (.not.allocated(batch_task_list)) then
         allocate(batch_task_list(n_grid_batches))
      end if

      batch_task_list = 0

      ! Set up batch description

      allocate(batch_desc(5,n_grid_batches))

      do i_batch=1,n_grid_batches
         batch_desc(1:3,i_batch) = coords(1:3,i_batch)
         batch_desc(4,i_batch)   = batch_sizes(i_batch)
         batch_desc(5,i_batch)   = i_batch
      enddo

      ! Distribute batches
      my_off = 0
      my_len = n_grid_batches
      call distribute_batches_by_location(batch_desc, 0, n_tasks-1, my_off, my_len)

      ! Set my batches in batch_task_list

      n_full_points = 0
      n_my_batches = my_len

      batch_task_list(:) = 0

      do i_batch = my_off+1, my_off+my_len

         batch_task_list(nint(batch_desc(5,i_batch))) = myid
         n_full_points = n_full_points + batch_sizes(nint(batch_desc(5,i_batch)))

      end do

      call sync_integer_vector(batch_task_list, n_grid_batches)

      if (use_mpi) then
         write(info_str,'(A)') ''
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,A,I6,A)') &
              "Integration load balanced across ", n_tasks, &
              " MPI tasks."
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,2A)') "Work distribution over tasks ", &
            "is as follows:"
         call localorb_info(info_str,use_unit,'(A)',OL_norm)
         write(info_str,'(2X,A,I6,A,I10,A)') &
              "Task ", myid, " has ", n_full_points, &
              " integration points."
         call localorb_allinfo(info_str,use_unit,'(A)',OL_low)
      end if

      deallocate(batch_desc)

      end subroutine distribute_batch_tasks_by_location
!******
!----------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_batches_by_location
!  NAME
!    distribute_batches_by_location
!  SYNOPSIS
      recursive subroutine distribute_batches_by_location &
                           (batch_desc,cpu_a,cpu_b,my_off,my_len)
!  PURPOSE
!    Recursively distributes the batches to the tasks by their location
!  USES

      implicit none
!  ARGUMENTS 
      real*8, intent(inout) :: batch_desc(:,:)
      integer, intent(in) :: cpu_a, cpu_b
      integer, intent(inout) :: my_off, my_len
!  INPUTS
!    o batch_desc(1,.) -- x-coord of batch (typically coord of some point in batch)
!    o batch_desc(2,.) -- y-coord of batch
!    o batch_desc(3,.) -- z-coord of batch
!    o batch_desc(4,.) -- weight of batch
!    o batch_desc(5,.) -- number of batch, must be 1,2,3... on top level call
!    o cpu_a -- number of first CPU for distribution
!    o cpu_b -- number of last CPU for distribution
!    o my_off -- offset of batches on current CPU set (cpu_a to cpu_b) within batch_desc
!                must be 0 on top level call 
!    o my_len -- number of batches on current CPU set
!                must be the number of batches on top level call
!  OUTPUT
!    o my_off -- offset of my batches on my CPU within batch_desc
!    o my_len -- number of my batches
!    o batch_desc -- is sorted so that batch_desc(:,my_off+1:my_off+my_len) 
!                    contains the corresponding values for my batches.
!                    Normally only batch_desc(5,:) is relevant on output.
!                    Values outside my_off+1:my_off+my_len are undefined 
!                    and should not be used!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      real*8 total_weight, target_weight
      integer median, cpu_m, n, icoord
      real*8 xmax(3), xmin(3), d

      ! If there is only 1 CPU for distribution, we are done

      if(cpu_a==cpu_b) return

      ! Get geometric extensions and total number of points of all batches on current CPU set

      total_weight = 0
      xmax = -1.d99
      xmin =  1.d99
      do n=my_off+1,my_off+my_len
         xmax(:) = max(xmax(:),batch_desc(1:3,n))
         xmin(:) = min(xmin(:),batch_desc(1:3,n))
         total_weight = total_weight+batch_desc(4,n)
      enddo

      ! Get dimension with biggest distance from min to max

      d = 0
      icoord = 1
      do n=1,3
         if(xmax(n)-xmin(n) > d) then
            d = xmax(n)-xmin(n)
            icoord = n
         endif
      enddo

      ! CPU number where to split CPU set

      cpu_m = (cpu_a+cpu_b-1)/2

      ! If the number of CPUs is not even, we have to split the set
      ! of batches accordingly into to differently sized halfes
      ! in order to end with an equal weight on every CPU.

      target_weight = total_weight*(cpu_m-cpu_a+1) / (cpu_b-cpu_a+1)

      ! divide in two equal halfes with the weight of the left half approx. equal to target_weight

      call divide_values(batch_desc(:,my_off+1:my_off+my_len), 5, my_len, icoord, &
                         target_weight, median)

      ! Set my_off and my_len for next distribution step

      if(myid<=cpu_m) then
         ! my_off stays the same
         my_len = median
      else
         my_off = my_off+median
         my_len = my_len-median
      endif

      ! If there are only two CPUs, we are done

      if(cpu_b == cpu_a+1) return

      ! Further divide both halves recursively

      if(myid<=cpu_m) then
         call distribute_batches_by_location(batch_desc,cpu_a,  cpu_m,my_off,my_len)
      else
         call distribute_batches_by_location(batch_desc,cpu_m+1,cpu_b,my_off,my_len)
      endif


      end subroutine distribute_batches_by_location
!******
!----------------------------------------------------------------------------------
!****s* mpi_utilities/divide_values
!  NAME
!    divide_values
!  SYNOPSIS
   subroutine divide_values(x, ldx, num, icoord, target, ndiv)
!  PURPOSE
!    Divides an array x(:,num) so that
!    - all x(icoord,1:ndiv) are smaller than all x(icoord,ndiv+1:num)
!    - SUM(x(4,1:ndiv)) == target (approximatly)
!  USES
    use localorb_io, only: use_unit
    implicit none
!  ARGUMENTS 
    real*8, intent(inout) :: x(ldx,num)
    integer, intent(in)   :: ldx, num, icoord
    real*8, intent(in)    :: target
    integer, intent(out)  :: ndiv
!  INPUTS
!    o x -- on entry contains unsorted coords, weights, additional values
!    o ldx -- leading dimension of x
!    o num -- number of data points in x
!    o icoord -- number of coord (row in x) to be used for sorting
!    o target -- target weight for left part of sorted values
!  OUTPUT
!    o x -- on exit contains sorted values as described under PURPOSE
!    o ndiv -- number of values in left (lower) part
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2011).
!  SOURCE

    integer, parameter :: iwgt = 4 ! where to get the weight from

    integer ns, ne, nl, nr, i
    real*8 xmin, xmax, xmid, sl_tot, sr_tot, sl, sr, xmax_l, xmin_r
    real*8 xtmp(ldx)

    ! Safety check only
    if(num<=0) then
      ! In this case we must have run out of grid batches to be distributed to
      ! individual MPI tasks. This can happen when very many CPU cores are used
      ! for relatively small structures.
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': Error when trying to distribute batches of grid points'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': across different MPI tasks.'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': Zero or negative number of batches for one or more CPUs.'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': This could be a consequence of using very many MPI tasks'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': for a structure that contains relatively few atoms.'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': You may want to adjust the keywords "points_in_batch" or '
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': "batch_size_limit" to have fewer points per batch'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': and therefore more overall batches to distribute.'
      write(use_unit,'(1X,A,I8,A)') & 
        '* MPI process myid = ',myid,': Please read the manual for details.'
      write(use_unit,'(1X,A,I8)') '* The error was encountered in subroutine divide_values: num = ',num
      call aims_stop('Too many MPI tasks for current grid partitioning.','divide_values')
    endif

    xmin = minval(x(icoord,:))
    xmax = maxval(x(icoord,:))

    ns = 1
    ne = num

    sl_tot = 0
    sr_tot = 0

    do

      xmid = 0.5*(xmin+xmax)
      ! Make sure that at least 1 value is below xmid and 1 value is above xmid.
      ! Otherways xmin and xmax are equal or only 1 bit apart and there is no more
      ! separation possible
      if(xmid<=xmin .or. xmid>=xmax) exit

      ! divide values between ns..ne into values <= xmid and values > xmid
      ! with a quicksort like algorithm

      nl = ns
      nr = ne
      sl = 0
      sr = 0

      xmax_l = xmin
      xmin_r = xmax

      do while(nl<=nr)

        ! skip values in front <= xmid

        do while(nl <= nr)
          if(x(icoord,nl) > xmid) exit
          sl = sl + x(iwgt,nl)
          xmax_l = max(xmax_l,x(icoord,nl))
          nl = nl+1
        enddo

        ! skip values in back > xmid

        do while(nl <= nr)
          if(x(icoord,nr) <= xmid) exit
          sr = sr + x(iwgt,nr)
          xmin_r = min(xmin_r,x(icoord,nr))
          nr = nr-1
        enddo

        if(nl>nr) exit

        ! Exchange elements at nr/nl

        xtmp(:) = x(:,nl)

        x(:,nl) = x(:,nr)
        sl = sl + x(iwgt,nl)
        xmax_l = max(xmax_l,x(icoord,nl))
        nl = nl+1

        x(:,nr) = xtmp(:)
        sr = sr + x(iwgt,nr)
        xmin_r = min(xmin_r,x(icoord,nr))
        nr = nr-1

      enddo

      ! Safety check:
      ! Check that at least 1 value is in both halves, this must be always the case!
      ! Otherways there is something screwed up with the program logic and we better exit
      if(nl==ns .or. nr==ne) then
        print *,'INTERNAL error in program logic of divide_values'
        call aims_stop ('INTERNAL ERROR')
      endif

      ! we can keep one half of the sorted values as is whereas the other half has to be sorted again

      if(sl_tot+sl < target) then
        ! Left is ok, right must be sorted again
        ns = nl
        sl_tot = sl_tot+sl
        xmin = xmin_r
      else
        ! Right is ok, left must be sorted again
        ne = nr
        sr_tot = sr_tot+sr
        xmax = xmax_l
      endif

    enddo

    ! Safety check
    if(ns>ne) then
      print *,'INTERNAL error in divide_values: ns=',ns,' ne=',ne
      call aims_stop ('INTERNAL ERROR')
    endif

    ! Now the value searched must be somewhere between ns an ne
    ! Please note that the coords of x(icoord,ns:ne) are (nearly) the same

    do i = ns, ne
      sl_tot = sl_tot + x(iwgt,i)
      if(sl_tot >= target) then
        ndiv = i
        return
      endif
    enddo

    ! Well, we should never come here unless there is something wrong with the weights ...

    ndiv = ne

  end subroutine divide_values

!******
!------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_batch_tasks_metis
!  NAME
!    distribute_batch_tasks_metis
!  SYNOPSIS
      subroutine distribute_batch_tasks_metis( batch_sizes, &
           grid_partition, n_points_in_grid )
!  PURPOSE
!    Distributes the grid batches to different tasks based on their location using an
!    external graph partitioning library, Metis
!  USES
      use dimensions, only: n_grid_batches, n_atoms, n_full_points, n_my_batches
      use runtime_choices
      use localorb_io
      use pbc_lists, only: coords_center
      use grids, only: r_radial, n_radial, r_angular, n_angular
      use geometry, only: species
!
!  ARGUMENTS 
      integer :: n_points_in_grid
      integer :: batch_sizes( n_grid_batches )
      integer :: grid_partition( n_points_in_grid )
!  INPUTS
!    o n_points_in_grid -- number of points in the grid
!    o batch_sizes -- number of grid points in each of the grid batches
!    o grid_partition -- partition of the grid into batches
!  OUTPUT
!    none -- at exit the array batch_task_list is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


!     locals
      integer :: i_atom, i_radial, i_angular
      integer :: i_batch, i_point, i_coord
      integer :: index, mpierr
      integer :: batch_task_work_table(n_tasks)

      real*8, dimension(3) :: coord_current
      real*8, dimension(:,:), allocatable :: batch_com
      character*100 :: info_str

      if (.not.allocated(batch_task_list)) then
         allocate(batch_task_list(n_grid_batches))
      end if

      batch_task_list = 0

      if (myid.eq.0) then
         allocate(batch_com(3,n_grid_batches))
         batch_com = 0.0d0

         i_point = 0
         do i_atom = 1, n_atoms, 1
            do i_radial = 1, n_radial(species(i_atom)), 1
               do i_angular = 1, n_angular( i_radial,species(i_atom) )

                  i_point = i_point + 1

                  coord_current(:) = coords_center( :,i_atom ) + &
                       r_angular(:,i_angular,i_radial,species(i_atom)) * &
                       r_radial(i_radial, species(i_atom))

                  batch_com(:,grid_partition(i_point)) = &
                       batch_com(:,grid_partition(i_point)) + &
                       coord_current(:)

               end do
            end do
         end do

         do i_coord = 1, 3, 1
            batch_com(i_coord,:) = &
                 batch_com(i_coord,:)/dble(batch_sizes(:))
         end do

         call metis_qhull_batches_wrapper( batch_com(1,:), &
              batch_com(2,:), batch_com(3,:), batch_sizes, &
              batch_task_list, n_grid_batches, n_tasks )

         deallocate( batch_com )

      end if

      call MPI_BCAST( batch_task_list, n_grid_batches, MPI_INTEGER, &
           0, mpi_comm_global, mpierr)

      batch_task_work_table = 0
      n_my_batches = 0
      n_full_points = 0

      do i_batch = 1, n_grid_batches, 1

         index = batch_task_list(i_batch) + 1
         batch_task_work_table(index) = &
              batch_task_work_table(index) + batch_sizes(i_batch)

         if (myid.eq.batch_task_list(i_batch)) then
            n_my_batches = n_my_batches + 1
            n_full_points = n_full_points + batch_sizes(i_batch)
         end if

      end do

      if (use_mpi) then
!         if (myid.eq.0) then
            write(info_str,'(A)') ''
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A,A,I5,A)') &
            "Integration load balanced ", &
            "across ", n_tasks, " MPI tasks. "
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            write(info_str,'(2X,A)') &
                 "Work distribution over tasks is as follows:"
            call localorb_info(info_str,use_unit,'(A)',OL_norm)
            if (output_priority .le. OL_low) then 
               do index = 1, n_tasks, 1
                  write(info_str,'(2X,A,I5,A,I10,A)') &
                       "Task ", index-1, &
                       " has ", &
                       batch_task_work_table(index), &
                       " integration points."
                  call localorb_info(info_str,use_unit,'(A)',OL_low)
               enddo
            end if
      end if

      end subroutine distribute_batch_tasks_metis
!******
!------------------------------------------------------------------------------
!****s* mpi_utilities/distribute_spline_storage
!  NAME
!    distribute_spline_storage
!  SYNOPSIS
      subroutine distribute_spline_storage ()
!  PURPOSE
!    Distributes the storage of the spline arrays to different tasks if the flag
!    use_distributed_spline_storage is set.
!  USES
      use dimensions, only: n_atoms, n_spline_atoms, &
          use_distributed_spline_storage
      implicit none
!  ARGUMENTS 
!  INPUTS
!    none
!  OUTPUT
!    none -- at exit the array spline_atom_storage is set
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


      integer :: i_atom, i_spline_atom

      if (.not.allocated(spline_atom_storage)) then
         allocate(spline_atom_storage(n_atoms))
      end if

      n_spline_atoms = 0
      spline_atom_storage = 0

      if (.not.use_distributed_spline_storage) then

         i_spline_atom = 0
         do i_atom = 1, n_atoms, 1
            i_spline_atom = i_spline_atom + 1
            spline_atom_storage(i_atom) = i_spline_atom
         end do
         n_spline_atoms = n_atoms

      else

         i_spline_atom = 0
         do i_atom = 1, n_atoms, 1

            if (myid.eq.task_list(i_atom)) then
               i_spline_atom = i_spline_atom + 1
               spline_atom_storage(i_atom) = i_spline_atom
            end if
         end do
         n_spline_atoms = i_spline_atom

      end if

      end subroutine distribute_spline_storage
!******
!------------------------------------------------------------------------------
!****s* mpi_utilities/finalize_mpi
!  NAME
!    finalize_mpi
!  SYNOPSIS
      subroutine finalize_mpi()
!  PURPOSE
!    Finalizes the MPI-environment before the program terminates.
!  USES
!
!  ARGUMENTS 
!
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mpierr

! All moved to deallocate_mpi
!      if (allocated(task_list)) then
!         deallocate(task_list)
!      end if
!      if (allocated(radial_task_list)) then
!         deallocate(radial_task_list)
!      end if
!      if (allocated(batch_task_list)) then
!         deallocate(batch_task_list)
!      end if
!      if (allocated(spline_atom_storage)) then
!         deallocate(spline_atom_storage)
!      end if

      call MPI_FINALIZE(mpierr)

      end subroutine finalize_mpi
!******
!****s* mpi_utilities/deallocate_mpi
!  NAME
!    get_my_task
!  SYNOPSIS
      subroutine deallocate_mpi()
!  PURPOSE
!     Deallocate all MPI Arrays 
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    Andrew Logsdail, UCL
!  HISTORY
!    Release version, FHI-aims (2013).
!  SOURCE

      if (allocated(task_list)) then
         deallocate(task_list)
      end if
      if (allocated(radial_task_list)) then
         deallocate(radial_task_list)
      end if
      if (allocated(batch_task_list)) then
         deallocate(batch_task_list)
      end if
      if (allocated(spline_atom_storage)) then
         deallocate(spline_atom_storage)
      end if

      end subroutine deallocate_mpi
!******
!****s* mpi_utilities/get_my_task
!  NAME
!    get_my_task
!  SYNOPSIS
      subroutine get_my_task()
!  PURPOSE
!    Find the task index of this thread and set the corresponding values
!    in mpi_tasks.f90.  As this subroutine is called in initialize_mpi.f90,
!    it should not be needed to be called again.
!  USES
      use localorb_io, only: use_unit
      use opencl_util, only: mpi_per_node, mpi_task_per_gpu
      implicit none
!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    n_tasks and myid are set on exit
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mpierr
      integer :: tmp_int

      call MPI_COMM_SIZE(mpi_comm_global, n_tasks, mpierr)
      call MPI_COMM_RANK(mpi_comm_global, myid, mpierr)

      call MPI_COMM_SPLIT(mpi_comm_global, myid / mpi_per_node, myid, shm_comm, mpierr);
      call MPI_COMM_SIZE(shm_comm, shm_num_of_rank, mpierr)
      call MPI_COMM_RANK(shm_comm, shm_rank, mpierr)


      if (shm_rank .eq. 0) then
         tmp_int = shm_rank
      else
         tmp_int = MPI_UNDEFINED
      end if

      call MPI_COMM_SPLIT(mpi_comm_global, tmp_int, myid, shm_per0_comm, mpierr);

      if (shm_rank .eq. 0) then
         call MPI_COMM_SIZE(shm_per0_comm, shm_per0_num_of_rank, mpierr)
         call MPI_COMM_RANK(shm_per0_comm, shm_per0_rank, mpierr)
      end if

      if (mpierr.ne.MPI_SUCCESS) then
         write(use_unit,'(1X,A)') "* get_my_task() failed"
         write(use_unit,'(1X,A)') "* Exiting..."
         stop
      end if

      end subroutine get_my_task
!******
!-------------------------------------------------------------------------
!****s* mpi_utilities/get_my_processor
!  NAME
!    get_my_processor
!  SYNOPSIS
      subroutine get_my_processor(name, length)
!  PURPOSE
!    Find the hosta name of each MPI-task.
!  USES
      use localorb_io, only: use_unit
      implicit none
!  ARGUMENTS
      character(LEN=MPI_MAX_PROCESSOR_NAME) :: name
      integer :: length
!  INPUTS
!    none
!  OUTPUT
!    o name -- character array for the hostname
!    o length -- length of the hostname
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: mpierr

      call MPI_GET_PROCESSOR_NAME( name, length, mpierr )

      if (mpierr.ne.MPI_SUCCESS) then
         write(use_unit,'(1X,A)') "* get_my_processor() failed"
         write(use_unit,'(1X,A)') "* Exiting..."
         stop
      end if

      end subroutine get_my_processor

!******
!------------------------------------------------------------------------------

      end module mpi_utilities
