!****h* FHI-aims/restart_pt2
!  NAME
!    restart_pt2
!  SYNOPSIS
module restart_pt2
!  PURPOSE
!     module to house restart options for pt2 correlation calculations
!     implement the triple k_grid roops from a given restart point
!     to a given stop point
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
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
!    Release version, FHI-aims (2014).
!  SOURCE
  implicit none

  logical :: restart_zero_iteration
  logical :: restart_first_warnings = .true.

!******	
contains

!------------------------------------------------------------------------------
!****s* restart/write_restart_pt2_info
!  NAME
!    write_restart_pt2_info
!  SYNOPSIS
  subroutine write_restart_pt2_info
!  PURPOSE
!    write restart information to file
!  USES
    use gw_para
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks
    use mpi_utilities
    use physics
    use runtime_choices

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    integer :: i_basis, i_states, i_spin, i_k
    character*150 :: info_str

    if (myid .eq. 0) then
        if (restart_pt2_write) then
           ! read from file, but only in task 0
           write(info_str,'(2X,2A)') 'Writing cluster restart information to file ', restart_pt2_write_file
           call localorb_info(info_str,use_unit,'(A)')

           open(file = restart_pt2_write_file, unit = 7, status = 'unknown', form = 'formatted')
           write(7,'(2X,I8,I8,3f21.10)') n_tasks, pt2_finish, current_pt2, &
               current_pt2_os, current_pt2_ss

           close(unit = 7)
        end if
    endif
  end subroutine write_restart_pt2_info
!******

!------------------------------------------------------------------------------
!****s* restart/read_restart_pt2_info
!  NAME
!    read_restart_pt2_info
!  SYNOPSIS
  subroutine read_restart_pt2_info
!  PURPOSE
!    read scf information from last known scf cycle from file
!  USES
    use gw_para
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks
    use mpi_utilities
    use synchronize_mpi
    use synchronize_mpi_basic
    use runtime_choices

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none
    character*40 :: form1
    character*150 :: info_str

    integer :: tmp_n_tasks

    if (myid .eq. 0) then
        if (restart_pt2_read) then
           ! read from file, but only in task 0
           write(info_str,'(2X,A)') ''
           call localorb_info(info_str,use_unit,'(A)')
           write(info_str,'(2X,2A)') 'Reading PT2 restart information from file ', restart_pt2_read_file
           call localorb_info(info_str,use_unit,'(A)')

           open(file = restart_pt2_read_file, unit = 7, status = 'old', form = 'formatted')
           read(7,'(2X,I8,I8,3f21.10)')  tmp_n_tasks, pt2_finish, current_pt2, &
               current_pt2_os, current_pt2_ss
           write(info_str,'(2X,A,I6,A)') 'Restart the PT2 calculation from', &
               pt2_finish+1, ' round'
           call localorb_info(info_str,use_unit,'(A)')

           ! Check the parameters allocated
           if (tmp_n_tasks .ne. n_tasks) then
               write(info_str,'(2X,A)') '* WARNING: restart file cannot belong to the current system.'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') '*          Please check the CPU number'
               call localorb_info(info_str,use_unit,'(A)')
               write(info_str,'(2X,A)') '*          Aborting execution'
               call localorb_info(info_str,use_unit,'(A)')
               call aims_stop
           endif
           close(unit = 7)
        end if
    endif
  end subroutine read_restart_pt2_info
!******

!------------------------------------------------------------------------------
!****s* restart/write_restart_pt2_info
!  NAME
!    write_restart_pt2_info
!  SYNOPSIS
  subroutine write_restart_pt2_info_blacs
!  PURPOSE
!    write restart information to file
!  USES
    use gw_para
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks
    use mpi_utilities
    use physics
    use runtime_choices

    implicit none
!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    integer :: i_basis, i_states, i_spin, i_k
    character*150 :: info_str

    if (myid .eq. 0) then
        if (restart_pt2_write) then
           ! read from file, but only in task 0
           write(info_str,'(2X,2A)') 'Writing cluster restart information to file ', restart_pt2_write_file
           call localorb_info(info_str,use_unit,'(A)')

           open(file = restart_pt2_write_file, unit = 7, status = 'unknown', form = 'formatted')
           write(7,'(2X,I8,3f21.10)') pt2_finish, current_pt2, &
               current_pt2_os, current_pt2_ss

           close(unit = 7)
        end if
    endif
  end subroutine write_restart_pt2_info_blacs

!******
!------------------------------------------------------------------------------
!****s* restart/read_restart_pt2_info
!  NAME
!    read_restart_pt2_info
!  SYNOPSIS
  subroutine read_restart_pt2_info_blacs
!  PURPOSE
!    read scf information from last known scf cycle from file
!  USES
    use gw_para
    use localorb_io, only: localorb_info, use_unit
    use mpi_tasks
    use mpi_utilities
    use runtime_choices
    use synchronize_mpi
    use synchronize_mpi_basic

!  AUTHOR
!    FHI-aims team.
!  HISTORY
!    Release version, FHI-aims (2014).
!  INPUTS
!    none
!  OUTPUT
!    none
!  SEE ALSO
!    FHI-aims CPC publication (in copyright notice above)
!  SOURCE
    implicit none
    character*40 :: form1
    character*150 :: info_str

    integer :: tmp_n_tasks

    if (myid .eq. 0) then
        if (restart_pt2_read) then
           ! read from file, but only in task 0
           write(info_str,'(2X,A)') ''
           call localorb_info(info_str,use_unit,'(A)')
           write(info_str,'(2X,2A)') 'Reading PT2 restart information from file ', restart_pt2_read_file
           call localorb_info(info_str,use_unit,'(A)')

           open(file = restart_pt2_read_file, unit = 7, status = 'old', form = 'formatted')
           read(7,'(2X,I8,3f21.10)')  pt2_finish, current_pt2, &
               current_pt2_os, current_pt2_ss
           write(info_str,'(2X,A,I6,A)') 'Restart the PT2 calculation from', &
               pt2_finish+1, ' round'
           call localorb_info(info_str,use_unit,'(A)')

           close(unit = 7)
        end if
    endif
  end subroutine read_restart_pt2_info_blacs
!******

end module restart_pt2
