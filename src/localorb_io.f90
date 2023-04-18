!****h* FHI-aims/localorb_io
!  NAME
!   localorb_io
!  SYNOPSIS

    module localorb_io

!  PURPOSE
!  Module contans routines for printing out information to output file.
!
!  USES

      implicit none

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
!******

      ! default priority for output: everything is printed.
      integer :: output_priority
      integer, parameter :: OL_low      = 0  ! Use for low priority messages
      integer, parameter :: OL_norm     = 1  ! Use for normal priority messages
      integer, parameter :: OL_high     = 2  ! Use for high priority messages
      character*100000, public  :: comm_string
      integer,public                  :: default_unit
      character*3, parameter, private :: default_format='(A)'
      integer,public                  :: use_unit
      character (len=80), public      :: SEPARATORLINE

      contains

      subroutine set_localorb_io_defaults ()
      implicit none

      output_priority = 1
      comm_string = ' '
      SEPARATORLINE = &
           "------------------------------------------------------------"

      end subroutine

!****s* localorb_io/localorb_info
!  NAME
!   localorb_info
!  SYNOPSIS

      subroutine localorb_info(message, unit, format, priority)


!  PURPOSE
!  The subroutine prints out line to output file
!
!  USES

      use ipc, only: ipc_start_transaction, ipc_write_string
      use mpi_tasks, only: myid
      implicit none

!  ARGUMENTS

      character(len=*) :: message
      integer, optional :: unit
      character(len=*), optional :: format
      integer, optional :: priority

!  INPUTS
!   o message -- printed out line
!   o unit -- file where the line is printed.
!   o format -- fomat of the printed line
!   o priority -- priority of printed message.
!
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE

      integer :: use_unit
      CHARACTER(LEN=4096) :: string_message

      if (present(priority)) then
         if (priority < output_priority) return
      end if
      if (present(unit)) then
         use_unit = unit
      else
         use_unit = default_unit
      end if

      if (myid == 0) then
         if (present(format)) then
             write(use_unit, format) trim(message)
             write(string_message, format) trim(message)
         else
             write(use_unit, default_format) trim(message)
             write(string_message, default_format) trim(message)
         end if
      end if

      ! Forward localorb message if localorb ipc is active.
      if(ipc_start_transaction("LOCALORB_MESSAGE")) then
         CALL ipc_write_string(string_message)
      end if


    end subroutine localorb_info
!******

subroutine localorb_info_alt(message)
    character(len=*), intent(in) :: message

    call localorb_info('  ' // message)
end subroutine

!****f* localorb_io/localorb_multi
! PURPOSE
!   Multiline version of localorb_info.
!******
subroutine localorb_multi( &
        msg1, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, msg10, &
        msg11, msg12, msg13, msg14, msg15, msg16, msg17, msg18, msg19, msg20, &
        unit, format, priority &
    )
    character(len=*), intent(in) :: msg1
    character(len=*), intent(in), optional :: &
        msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, msg10, &
        msg11, msg12, msg13, msg14, msg15, msg16, msg17, msg18, msg19, msg20
    integer, intent(in), optional :: unit
    character(len=*), intent(in), optional :: format
    integer, intent(in), optional :: priority

    call localorb_info(msg1, unit, format, priority)
    if (.not. present(msg2)) return
    call localorb_info(msg2, unit, format, priority)
    if (.not. present(msg3)) return
    call localorb_info(msg3, unit, format, priority)
    if (.not. present(msg4)) return
    call localorb_info(msg4, unit, format, priority)
    if (.not. present(msg5)) return
    call localorb_info(msg5, unit, format, priority)
    if (.not. present(msg6)) return
    call localorb_info(msg6, unit, format, priority)
    if (.not. present(msg7)) return
    call localorb_info(msg7, unit, format, priority)
    if (.not. present(msg8)) return
    call localorb_info(msg8, unit, format, priority)
    if (.not. present(msg9)) return
    call localorb_info(msg9, unit, format, priority)
    if (.not. present(msg10)) return
    call localorb_info(msg10, unit, format, priority)
    if (.not. present(msg11)) return
    call localorb_info(msg11, unit, format, priority)
    if (.not. present(msg12)) return
    call localorb_info(msg12, unit, format, priority)
    if (.not. present(msg13)) return
    call localorb_info(msg13, unit, format, priority)
    if (.not. present(msg14)) return
    call localorb_info(msg14, unit, format, priority)
    if (.not. present(msg15)) return
    call localorb_info(msg15, unit, format, priority)
    if (.not. present(msg16)) return
    call localorb_info(msg16, unit, format, priority)
    if (.not. present(msg17)) return
    call localorb_info(msg17, unit, format, priority)
    if (.not. present(msg18)) return
    call localorb_info(msg18, unit, format, priority)
    if (.not. present(msg19)) return
    call localorb_info(msg19, unit, format, priority)
    if (.not. present(msg20)) return
    call localorb_info(msg20, unit, format, priority)
end subroutine

!****s* localorb_io/localorb_allinfo
!  NAME
!   localorb_allinfo
!  SYNOPSIS

    subroutine localorb_allinfo(message, unit, format, priority)

      !  PURPOSE
      !     The subroutine prints out a line for each proc concertedly
      !  USES

      use mpi_tasks
      implicit none

      !  ARGUMENTS

      character(len=*) :: message
      integer, optional :: unit
      character(len=*), optional :: format
      integer, optional :: priority

      !  INPUTS
      !   o message -- printed out line
      !   o unit -- file where the line is printed.
      !   o format -- fomat of the printed line
      !   o priority -- priority of printed message.
      !  OUTPUT
      !    none
      !  AUTHOR
      !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
      !  HISTORY
      !    Release version, FHI-aims (2008).
      !  SOURCE

      integer, parameter :: outsize = 300
      integer :: use_unit
      integer :: i_task
      character(len=outsize) :: local_buffer
      character(len=outsize), allocatable :: global_buffer(:)
      integer :: mpierr, info
      character(*), parameter :: func = 'localorb_allinfo'

      if (present(priority)) then
         if (priority < output_priority) return
      end if

      if (present(format)) then
         write(local_buffer, format) trim(message)
      else
         write(local_buffer, default_format) trim(message)
      end if

      if (myid == 0) then
         allocate(global_buffer(n_tasks), stat=info)
      else
         allocate(global_buffer(1), stat=info)
      end if
      call check_allocation(info, 'global_buffer', func)
      if (use_mpi) then
         call MPI_Gather(local_buffer, outsize, MPI_CHARACTER, &
         &               global_buffer, outsize, MPI_CHARACTER, &
         &               0, mpi_comm_global, mpierr)
         if (mpierr /= 0) call aims_stop('MPI_Gather error',func)
      else
         write(global_buffer(1), "(A)") trim(local_buffer)
      end if

      if (myid == 0) then
         if (present(unit)) then
            use_unit = unit
         else
            use_unit = default_unit
         end if
         do i_task = 1, n_tasks
            if (global_buffer(i_task).ne."") then
               write(use_unit, "(A)") trim(global_buffer(i_task))
            end if
         end do
      end if
      deallocate(global_buffer)

    end subroutine localorb_allinfo
!******

!****f* localorb_io/str_lower
! PURPOSE
!   Converts a string to lower characters.
!******
function str_lower(str) result(lower)
    character(len=*), intent(in) :: str
    character(len=len(str)) :: lower

    integer :: i

    do i = 1, len(str)
        select case (str(i:i))
        case ('A':'Z')
            lower(i:i) = achar(iachar(str(i:i))+32)
        case default
            lower(i:i) = str(i:i)
        end select
    end do
end function

function str_replace(str, from, to) result(replace)
    character(len=*), intent(inout) :: str
    character(len=1), intent(in) :: from, to
    character(len=len(str)) :: replace

    integer :: i

    replace = str
    do i = 1, len(str)
        if (replace(i:i) == from) replace(i:i) = to
    end do
end function

      end module localorb_io
