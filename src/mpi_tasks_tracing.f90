!****h* FHI-aims/mpi_tasks
!  NAME
!    mpi_tasks - provides basic utilities for MPI calculations
!  SYNOPSIS
      module mpi_tasks
!  PURPOSE
!    This module provides the very basic MPI-tools: finding the task number
!    and hostname.
!  USES
      use ifcore !this non-standard module is provided by Intel Fortran
      implicit none

      logical :: use_mpi    ! should we use MPI-based parallellism
      logical :: use_mpi_in_place
      integer :: n_tasks
      integer :: myid
      integer :: mpi_comm_global

      integer :: kpt_comm   ! One alternative communicator for possible
                            ! restriction of communication to k-point
                            ! carrying tasks only.
                            ! Only used in periodic Hartree Fock for now.

      integer, parameter :: STDERR = 0
      integer, parameter :: STDOUT = 6

      include 'mpif.h'
!  ARGUMENTS
!    none
!  INPUTS
!    none
!  OUTPUT
!    none
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
!
!******

      contains
!-------------------------------------------------------------------------
!****s* mpi_tasks/aims_stop
!  NAME
!    aims_stop
!  SYNOPSIS

      subroutine aims_stop(message, caller)

!  PURPOSE
!   If only one of processes is calling stop then in some computers
!   (e.g. BlueGene) the run is not stopping, but only get stuck to the state
!   where nothing is happening.  In this case it is better to call MPI_Abort
!   which kills all the processes.  So this is subroutine for replacing the
!   stop call.
!  USES
      use ipc, only: ipc_start_transaction, ipc_write_string
      implicit none

!  ARGUMENTS

      character(len=*), intent(IN), optional :: message
      character(len=*), intent(IN), optional :: caller


!  INPUTS
!    o message (optional) -- (error) message to be output. In a parallel run,
!                            the message might be output once per process
!                            in order to make really sure that it is output at least
!                            once.
!  OUTPUT
!    none
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
! SOURCE

        integer :: mpierr
        CHARACTER(LEN=4096) :: string_message

        if (present(message) .and. present(caller)) then
           write(string_message, "(1X,'*** Proc',I5,' in ',A,': ',A)") &
           & myid, trim(caller), trim(message)
        else if (present(message)) then
           write(string_message, "(1X,'*** Proc',I5,': ',A)") myid, &
           & trim(message)
        else if (present(caller)) then
           write(string_message, "(1X,'*** Proc',I5,' in ',A,': Error.')") &
           & myid, trim(caller)
        else
           write(string_message, "(1X,'*** Proc',I5,' : Error.')") myid
        end if

        write(STDOUT,'(A)') trim(string_message)
        write(STDERR,'(A)') trim(string_message)
        ! One extra instance. I do not think this is necessary.
        ! Commented.
        ! write(use_unit,'(A)') trim(string_message)

        if(ipc_start_transaction("AIMS_STOP")) then
           CALL ipc_write_string(trim(string_message))
        endif

        ! forcing a stack trace dump with Intel tools
        call TRACEBACKQQ("StackTrace of failing process", -1)
        if (n_tasks > 1) then
           call MPI_Abort(mpi_comm_global, 0, mpierr)
        end if
        stop

      end subroutine aims_stop
!******
!-------------------------------------------------------------------------
!****s* mpi_tasks/aims_warn
!  NAME
!    aims_warn
!  SYNOPSIS

      subroutine aims_warn(message, caller)

!  PURPOSE
!   The purpose for this routine is to give a warning if something strange
!   is happening on one process, but not stop the calculation.
!  USES

        implicit none

!  ARGUMENTS

      character(len=*), intent(IN), optional :: message
      character(len=*), intent(IN), optional :: caller


!  INPUTS
!    o message (optional) -- (warning) message to be output. In a parallel run,
!                            the message might be output once per process
!                            in order to make really sure that it is
!                            output at least once.
!  OUTPUT
!    none
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
! SOURCE


        integer :: mpierr

        if (present(message) .and. present(caller)) then
           write(STDOUT, "(1X,'*** Proc',I5,' in ',A,': ',A)") &
           & myid, trim(caller), trim(message)
        else if (present(message)) then
           write(STDOUT, "(1X,'*** Proc',I5,': ',A)") myid, trim(message)
        else if (present(caller)) then
           write(STDOUT, "(1X,'*** Proc',I5,' in ',A,': Warning.')") &
              myid, trim(caller)
        end if

      end subroutine aims_warn
!******
!****s* FHI-aims/check_allocation
!  NAME
!     check_allocation
!  SYNOPSIS

      subroutine check_allocation(info, name, caller, &
            dimension1, dimension2, dimension3, dimension4)

        !  PURPOSE
        !
        !     Check for memory allocation failures.
        !
        !     One might argue that this is complete nonsense because the
        !     default behavior after a failed allocation is to stop execution
        !     with a clear error message.  Unfortunately, some quite common
        !     platforms do *not* behave this way.
        !
        !     Paula Havu is more specific on this issue in the forum:
        !
        !         The reason because BlueGene needs the allocation checks (and
        !         before any ideas which could produce problems):
        !
        !         Sometimes there is problem with allocation specially with
        !         long runs. Without allocation checks, only the process with
        !         the problem collapses, and all the rest start waiting that
        !         one. And they wait before user or time limits kill the
        !         job. This takes computer resources without any result. This
        !         is why we need an allocation checks, which now calls
        !         aims_stop, which closes the total job. Note that normal stop
        !         is not a solution. This is also answer to question, why we
        !         need aims_stop (should always used instead of stop).
        !
        !         Note that allocation problem can appear in any variable,
        !         even small ones.
        !
        !
        !  USES

        implicit none

        ! ARGUMENTS

        integer, intent(IN) :: info
        character(*), intent(IN) :: name
        character(*), intent(IN), optional :: caller
        integer, intent(IN), optional :: dimension1
        integer, intent(IN), optional :: dimension2
        integer, intent(IN), optional :: dimension3
        integer, intent(IN), optional :: dimension4

        !  INPUTS
        ! o  info -  info integer from allocation
        ! o  name -  name of the variable. Name is printed out in the case of
        !            allocation did not work.
        ! o  caller - (optional) Additionally output caller name
        !  OUTPUT
        !   none
        !  AUTHOR
        !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
        !  SEE ALSO
        !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
        !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler, "Ab initio
        !    simulations with Numeric Atom-Centered Orbitals: FHI-aims",
        !    Computer Physics Communications (2008), submitted.
        !  COPYRIGHT
        !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
        !   e.V. Please note that any use of the "FHI-aims-Software" is
        !   subject to the terms and conditions of the respective license
        !   agreement."
        !  HISTORY
        !    Release version, FHI-aims (2008).
        !  SOURCE

        character*150 :: info_str

        if(info /= 0) then
           ! Write error message directly here, then we can use more than one line.

           write(info_str, '(1X,A,A,A)')  'Error in allocation of array ', trim(name), '.'

           if (present(caller)) then
              write(STDOUT, "(1X,'*** Proc',I10,' in ',A,': ',A)") &
              & myid, trim(caller), trim(info_str)
           else
              write(STDOUT, "(1X,'*** Proc',I10,': ',A)") &
              & myid, trim(info_str)
           end if

           if (present(dimension1)) write(STDOUT,"(1X,A,I10,A)") '*** The requested dimension (1) of the array was: ', dimension1, '.'
           if (present(dimension2)) write(STDOUT,"(1X,A,I10,A)") '*** The requested dimension (2) of the array was: ', dimension2, '.'
           if (present(dimension3)) write(STDOUT,"(1X,A,I10,A)") '*** The requested dimension (3) of the array was: ', dimension3, '.'
           if (present(dimension4)) write(STDOUT,"(1X,A,I10,A)") '*** The requested dimension (4) of the array was: ', dimension4, '.'

           write(STDOUT, '(1X,A)') &
             '*** This error may mean that you ran out of memory. If possible, check if more CPUs help.'
           write(STDOUT, '(1X,A)') &
             '*** In any case, reporting the problem on aimsclub would not hurt.'

           call aims_stop
        end if

      end subroutine check_allocation
!******
  !----------------------------------------------------------------------------
  !****s* mpi_tasks/aims_stop_coll
  !  NAME
  !    aims_stop_coll
  !  SYNOPSIS

  subroutine aims_stop_coll(message, caller)

    !  PURPOSE
    !
    !    This is a more graceful version of aims_stop() at the price of being
    !    collective (each process needs to call it).
    !
    !    It won't stop until all processors have arrived here; allowing all
    !    localorb_info() calls to be completed.
    !
    !    ATTENTION: Never, really NEVER call this procedure if you are not
    !    completely sure that all of the processors arrive here.  If one of
    !    them is missing, FHI-aims will hang, possibly without even having
    !    printed out the source of trouble.
    !
    !  USES

    implicit none

    !  ARGUMENTS

    character(len=*), intent(IN), optional :: message
    character(len=*), intent(IN), optional :: caller

    !  INPUTS
    !    o message -- Error message
    !    o func -- Function which calls this procedure
    !  OUTPUTS
    !    none
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    character(*), parameter :: func = 'aims_stop_coll'

    integer :: mpierr
    character*200 :: info_str
    if (present(message) .and. present(caller)) then
       write(info_str, "(1X,'*** ',A,': ',A)") trim(caller), trim(message)
    else if (present(message)) then
       write(info_str, "(1X,'*** ',A)") trim(message)
    else if (present(caller)) then
       write(info_str, "(1X,'*** ',A,': stop.')") trim(caller)
    else
       write(info_str, "(1X,'*** stop.')")
    end if
    if (myid == 0) then
       write(STDOUT, "(A)") trim(info_str)
       write(STDERR, "(A)") trim(info_str)
    end if
    if (use_mpi) call MPI_FINALIZE(mpierr)
    stop

  end subroutine aims_stop_coll
  !******

  !----------------------------------------------------------------------------
  !****s* mpi_tasks/setup_kpoint_communicator
  !  NAME
  !    setup_kpoint_communicator
  !  SYNOPSIS

  subroutine setup_kpoint_communicator (prev_n_k_points, n_k_points)

    !  PURPOSE
    !
    !  Handles the setup of a possible separate communicator for the group
    !  of processors among which
    !  k-points are distributed. Significant only if the number of CPUs
    !  is greater than the number of k-points, otherwise is the same as
    !  the global communicator
    !
    !  Currently, only used in periodic Hartree-Fock et al.
    !
    !  USES

    implicit none

    !  ARGUMENTS
    integer,intent(IN) :: prev_n_k_points
    integer,intent(IN) :: n_k_points
    !  INPUTS
    !    prev_n_k_points : This integer serves as a flag to see
    !                      (1) whether we have been here before at all (else it's 0)
    !                      (2) whether the communicator needs to be changed at all.
    !    n_k_points      : Number of k points
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2010).
    !  SOURCE

    integer:: mpierr
    integer:: world_group, kpt_group
    integer, dimension(3,1) ::  kpt_group_ranks

    if (prev_n_k_points.ne.0) then
      ! In this case, we have been here before, but the number of k-points has changed.
      ! Free old communicator first before creating a new one.
      ! ...but only do so if we are not using the global communicator.
      if (kpt_comm.ne.mpi_comm_global) then
         call mpi_comm_free(kpt_comm,mpierr)
      end if
    endif

    if(n_k_points.ge.n_tasks)then
       kpt_comm = mpi_comm_global
    else
       call MPI_comm_group(mpi_comm_global, world_group, mpierr)
       kpt_group_ranks(3,1) = 1
       kpt_group_ranks(1,1) = 1
       kpt_group_ranks(2,1) = n_k_points
       call MPI_Group_range_incl(world_group, 1, kpt_group_ranks, kpt_group, mpierr)
       call MPI_Comm_create(mpi_comm_global, kpt_group, kpt_comm, mpierr)
    endif

  end subroutine setup_kpoint_communicator

!******
!------------------------------------------------------------------------------

end module mpi_tasks
