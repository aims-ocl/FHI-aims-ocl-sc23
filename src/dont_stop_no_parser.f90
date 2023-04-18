!****s* FHI-aims/stop_if_parser
!  NAME
!    stop_if_parser
!  SYNOPSIS

subroutine stop_if_parser()

  !  PURPOSE
  !
  !    Dummy only.
  !
  !  USES

  use runtime_choices
  use localorb_io
  use mpi_tasks
  use pbc_lists

  implicit none

  !  ARGUMENTS
  !    none
  !  INPUTS
  !    none
  !  OUTPUTS
  !    none
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2011).
  !  SOURCE

  character(*), parameter :: func = 'stop_if_parser(false)'

  if (dry_run) then

     ! These lists must still be checked for consistency.
     call initialize_bc_dependent_lists()

     call localorb_info("")
     call localorb_info("----------------------------------------")
     call localorb_info("  This was a dry run only, to check all input files.")
     call localorb_info("  So far, all input files must have looked valid.")
     call localorb_info("  This run will be stopped at this point.")
     call localorb_info("----------------------------------------")

     skip_SCF = .true.  ! The grids are not allocated yet ...
     call final_deallocations ( )

     ! Final timings are not too informative, either.
     ! call final_timings (.true.)

     call localorb_info("  Stopping.")

     call aims_stop_coll('End of the dry run was reached, apparently successfully so.','stop_if_parser')

  else
    ! Do nothing
    continue
  end if

end subroutine stop_if_parser
!******
