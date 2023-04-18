!****h* FHI-aims/control_file
! PURPOSE
!   Helper routines for processing control.in
! NOTES
!   This module was created to share code betewen read_control and 
!   read_control_update.
! AUTHOR
!   Jan Hermann
! CREATION DATE
!   2015-07-20
! COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften e.V. Please note 
!   that any use of the "FHI-aims-Software" is subject to the terms and 
!   conditions of the respective license agreement.
!******
module control_file
    use mpi_tasks, only: myid
    use localorb_io, only: use_unit, localorb_info

    private
    public :: run_all_parsers

    logical, public :: &
        flag_acc_forces, flag_acc_stress

contains

    !****f* FHI-aims/run_all_parsers
    ! PURPOSE
    !   Runs sequentially all parsers defined here on a given control.in line until 
    !   some parser actually parses it.
    ! NOTES
    !   The individual parser is handed the current descriptor and line and should 
    !   return:
    !   * 0: If it finds the descriptor that it deals with and ends normally.
    !   * 1: If it doesn't find anything.
    !   * 88: If it hits the end of file.
    !   * 99: If error occurs when reading.
    !******
    function run_all_parsers (desc_str, inputline) result(ex_code)
        implicit none

        ! inout so we can just copy-paste code from read_control.f90
        character(len=*), intent(inout) :: desc_str
        character(len=*), intent(in) :: inputline
        integer :: ex_code

        ex_code = parse_sc_accuracy_flags(desc_str, inputline)
        if (ex_code /= 1) return
        ! ex_code = parse...
        ! if (ex_code /= 1) return
        ! ex_code = parse...
        ! if (ex_code /= 1) return
        ! etc.
    end function run_all_parsers

    function parse_sc_accuracy_flags(desc_str, inputline) result(exit_code)
        use runtime_choices, only: &
            flag_acc_rho, &
            flag_acc_eev, &
            flag_acc_etot, &
            flag_acc_potjump, &
            sc_accuracy_eev, &
            sc_accuracy_rho, &
            sc_accuracy_etot, &
            sc_accuracy_potjump, &
            sc_accuracy_forces, &
            sc_accuracy_stress

        implicit none

        character(len=*), intent(inout) :: desc_str
        character(len=*), intent(in) :: inputline
        integer :: exit_code

        character(len=200) :: info_str

        exit_code = 0
        select case (desc_str)
        case ('sc_accuracy_rho')
            read (inputline, *, end=88, err=99) desc_str, sc_accuracy_rho
            if (myid == 0) then
                write (use_unit, '(2X,A,E11.4)') &
                    "Convergence accuracy of self-consistent charge density: ", &
                    sc_accuracy_rho
            end if
            flag_acc_rho = .true.
        case ('sc_accuracy_eev')
            read (inputline, *, end=88, err=99) desc_str, sc_accuracy_eev
            if (myid == 0) then
                write (use_unit, '(2X,A,E11.4)') &
                    "Convergence accuracy of sum of eigenvalues: ", &
                    sc_accuracy_eev
            end if
            flag_acc_eev = .true.
        case ('sc_accuracy_etot')
            read(inputline, *, end=88, err=99) desc_str, sc_accuracy_etot
            if (myid == 0) then
                write (use_unit, '(2X,A,E11.4)') &
                    "Convergence accuracy of total energy: ", &
                    sc_accuracy_etot
            end if
            flag_acc_etot = .true.
        case ('sc_accuracy_potjump')
            read(inputline, *, end=88, err=99) desc_str, sc_accuracy_potjump
            if (myid == 0) then
                write (use_unit, '(2X,A,E11.4)') &
                    "Convergence accuracy of potential jump (eV): ", &
                    sc_accuracy_potjump
            end if
            flag_acc_potjump = .true.
        case('sc_accuracy_forces')
            read(inputline, *, end=88, err=99) desc_str, desc_str
            if (desc_str == 'not_checked') then
                sc_accuracy_forces = -5.d0
                write (info_str, '(2X,A)') &
                    "| Force self-consistency will not be checked explicitly."
                call localorb_info (info_str)
            else
                read (inputline, *, end=88, err=99) desc_str, sc_accuracy_forces
                if (myid == 0) then
                    write (use_unit,'(2X,A,E11.4)') &
                        "Convergence accuracy of forces: ", sc_accuracy_forces
                    if (sc_accuracy_forces < 0.d0) then
                        write (use_unit, '(2X,A,A)') &
                            "| Zero or negative force accuracy criterion: ", &
                            "Force self-consistency will not be checked explicitly."
                    end if
                end if
            end if
            flag_acc_forces = .true.
        case('sc_accuracy_stress')
            read (inputline, *, end=88, err=99) desc_str, desc_str
            if (desc_str == 'not_checked') then
                sc_accuracy_stress = -5.d0
                write (info_str, '(2X,A)') &
                    "| Analytical stress self-consistency will not be checked explicitly."
                call localorb_info(info_str)
            else
                read(inputline, *, end=88, err=99) desc_str, sc_accuracy_stress
                if (myid == 0) then
                    write (use_unit,'(2X,A,E11.4)') &
                        "Convergence accuracy of analytical stress: ", &
                        sc_accuracy_stress
                    if (sc_accuracy_stress < 0.d0) then
                        write(use_unit, '(2X,A,A)') &
                            "| Zero or negative stress accuracy criterion: ", &
                            "Analytical stress self-consistency will not be checked explicitly."
                    end if
                end if
            end if
            flag_acc_stress = .true.
        case default
            exit_code = 1
        end select
        return
        88 continue
        exit_code = 88
        return
        99 continue
        exit_code = 99
        return
    end function parse_sc_accuracy_flags
end module control_file
