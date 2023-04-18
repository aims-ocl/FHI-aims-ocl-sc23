!****h* FHI-aims/restart_rpa
!  NAME
!    restart_rpa
!  SYNOPSIS
module restart_rpa
!  PURPOSE
!     module to house restart options for rpa correlation calculations
!     implement the frequence integration from a given restart point
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
!****s* restart/write_restart_rpa_info
!  NAME
!    write_restart_info
!  SYNOPSIS
  subroutine write_restart_rpa_info
!  PURPOSE
!    write restart information to file
!  USES
    use runtime_choices
    use gw_para
    use mpi_tasks
    use mpi_utilities
    use physics
    use localorb_io, only: localorb_info ,use_unit

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

    if (restart_rpa_write) then
       ! read from file, but only in task 0
       if (myid.eq.0) then              ! work only on CPU 0
          write(info_str,'(2X,2A)') 'Writing cluster restart information to file ', restart_rpa_write_file
          call localorb_info(info_str,use_unit,'(A)')

          open(file = restart_rpa_write_file, unit = 7, status = 'unknown', form = 'formatted')
          write(7,'(2X,I8,I8,f21.10)') n_full_freq, rpa_freq_start, current_crpa

          close(unit = 7)
       end if      ! (myid.eq.0)
    end if
  end subroutine write_restart_rpa_info
!******

!------------------------------------------------------------------------------
!****s* restart/read_restart_rpa_info
!  NAME
!    read_restart_rpa_info
!  SYNOPSIS
  subroutine read_restart_rpa_info
!  PURPOSE
!    read scf information from last known scf cycle from file
!  USES
    use runtime_choices
    use gw_para
    use mpi_tasks
    use mpi_utilities
    use synchronize_mpi
    use synchronize_mpi_basic
    use localorb_io, only: localorb_info, use_unit

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

    integer :: tmp_n_full_freq

    if (restart_rpa_read) then
       ! read from file, but only in task 0
       if (myid.eq.0) then
          write(info_str,'(2X,A)') ''
          call localorb_info(info_str,use_unit,'(A)')
          write(info_str,'(2X,2A)') 'Reading RPA restart information from file ', restart_rpa_read_file
          call localorb_info(info_str,use_unit,'(A)')

          open(file = restart_rpa_read_file, unit = 7, status = 'old', form = 'formatted')
          read(7,'(2X,I8,I8,f21.10)') tmp_n_full_freq, rpa_freq_start, current_crpa
          write(use_unit,'(2X,I8,I8,f21.10)') tmp_n_full_freq, rpa_freq_start, current_crpa

          ! Check the parameters allocated
          if (tmp_n_full_freq .ne. n_full_freq) then
              write(info_str,'(2X,A)') '* WARNING: restart file cannot belong to current system.'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') '*          Please check n_full_freq'
              call localorb_info(info_str,use_unit,'(A)')
              write(info_str,'(2X,A)') '*          Aborting execution'
              call localorb_info(info_str,use_unit,'(A)')
              call aims_stop
          endif


          close(unit = 7)
       end if      ! (myid.eq.0)

       ! MPI broadcast to all other threads
       call sync_integer(rpa_freq_start)
       call sync_real_number(current_crpa)

    end if
  end subroutine read_restart_rpa_info
!******

end module restart_rpa
