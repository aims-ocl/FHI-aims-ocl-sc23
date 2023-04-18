!****f* FHI-aims/read_control_update
! PURPOSE
!   Reads control.update.in and updates settings accordingly.
! USAGE
!   It should be typically called after each SCF step.
! SEE ALSO
!   This was forked from read_control.f90.
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2015-07-20
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note 
!   that any use of the "FHI-aims-Software" is subject to the terms and 
!   conditions of the respective license agreement.
!******
subroutine read_control_update ()
    use localorb_io, only: use_unit, SEPARATORLINE, localorb_info
    use mpi_tasks, only: myid, aims_stop_coll
    use control_file, only: run_all_parsers

    implicit none

    character(len=*), parameter :: func = 'read_control_update'
    integer :: i_code, linecount, io_unit, ex_code
    character(len=132) inputline
    character(len=40) desc_str

    call localorb_info ('', use_unit)
    call localorb_info (SEPARATORLINE, use_unit, '(A)')
    call localorb_info ("Reading file control.update.in.", use_unit, '(10X,A)')
    call localorb_info (SEPARATORLINE, use_unit, '(A)')

    open (get_io_unit(io_unit), file="control.update.in")
    linecount = 0
    do
        read (io_unit, '(A)', iostat=i_code) inputline
        if (i_code < 0) exit 
        if (i_code > 0) then
            call aims_stop_coll("Unknown error reading file 'control.update.in'...", func)
        end if
        linecount = linecount+1
        read (inputline, *, iostat=i_code) desc_str
        if (i_code /= 0) cycle 
        if (desc_str(1:1) == '#') cycle 
        ex_code = run_all_parsers(desc_str, inputline)
        select case (ex_code) ! case (0) is normal exit, nothing is done
        case (88)
            go to 88
        case (99)
            go to 99
        case (1) ! no descriptor was found
            if (myid == 0) then
                write(use_unit, *) &
                    "Unknown descriptor ", desc_str, " in file control.in."
            end if
            return
        end select
    end do
    close (io_unit)
    if (myid == 0) then
        write (use_unit,* )
        write (use_unit, '(2X,A)') "Finished reading input file 'control.update.in'."
        write (use_unit, '(2X,A)') &
            "There are no checks done, make sure you know what you are doing."
        write (use_unit, *)
    end if
    call localorb_info (SEPARATORLINE, use_unit, '(A)')
    call localorb_info ('', use_unit)
    return
    88 continue
    write(use_unit, *) "Syntax error reading 'control.update.in' (missing arguments)"
    write(use_unit, *) "line: '" // trim(inputline) // "'"
    99 continue
    write(use_unit, *) "Syntax error reading 'control.update.in'"
    write(use_unit, *) "line: '" // trim(inputline) // "'"

contains

    integer function get_io_unit(io_unit_arg) result(io_unit)
        implicit none

        integer, intent(out), optional :: io_unit_arg

        logical :: unit_open
        integer :: i_io

        io_unit = -1
        do i_io = 10, 100
            inquire (unit=i_io, opened=unit_open)
            if (.not. unit_open) then
                io_unit = i_io
                exit
            end if
        end do
        if (present(io_unit_arg)) then
            io_unit_arg = io_unit
        end if
    end function get_io_unit
end subroutine read_control_update
