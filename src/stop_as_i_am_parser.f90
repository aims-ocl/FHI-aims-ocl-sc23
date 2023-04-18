!****s* FHI-aims/stop_if_parser
!  NAME
!    stop_if_parser
!  SYNOPSIS

subroutine stop_if_parser()

  !  PURPOSE
  !
  !    Called if this is a pure parser.  Deallocate, tidy up, output timings,
  !    and exit.
  !
  !  USES

  use timing
  use mpi_utilities
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

  character(*), parameter :: func = 'stop_if_parser(true)'

  call initialize_bc_dependent_lists()

  call localorb_info("")
  call localorb_info("----------------------------------------")
  call localorb_info("  Input files look valid.")

  skip_SCF = .true.  ! The grids are not allocated yet ...
  call final_deallocations ( )

  ! Final timings are not too informative, either.
  ! call final_timings (.true.)

  call localorb_info("  Stopping.")
  call finalize_mpi ( )
  stop

end subroutine stop_if_parser
!******
