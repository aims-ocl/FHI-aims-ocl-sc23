
!****s* FHI-aims/stop_illconditioning
!  NAME
!    stop_illconditioning
!  SYNOPSIS

subroutine stop_illconditioning

!  PURPOSE
!
!  USES

  use mpi_tasks
  implicit none

!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals",
!    Computer Physics Communications (2009).
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2009).
! SOURCE


  write(STDERR,'(1X,A,I8,A)') "* Attention: Ill-conditioned overlap matrix detected in thread number", myid, "!"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": This may mean that you are running with a very large basis set,"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": that you are running with an excessively large cutoff radius"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": in a periodic calculation, or that you are using diffuse Gaussian"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": basis_functions with an insufficient integration grid (try to "
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": increase radial_multiplier for all species in that case). You"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": should not usually encounter this behaviour in normal production calculations."
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": In any case, this behavior can often simply be overcome by setting"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": a large enough value of the basis_threshold parameter in control.in ."
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, &
    ": However, even then you MUST NOW ALSO SET the flag 'override_illconditioning .true.'"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": explicitly - else FHI-aims will stop (as in this case)."
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": Be ABSOLUTELY SURE to only ever set the 'override_illconditioning' flag"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": if you understand what you are doing, and always check convergence thoroughly."
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": Otherwise, running with a near-singular basis set may lead to absolutely"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": uncontrollable, completely wrong numerical results. When using basis_threshold"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": in an intelligent way, there is usually no problem - but by forcing you to set"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": the 'override_illconditioning' flag first, we are trying to make sure that you"
  write(STDERR,'(1X,A,I8,A)') "* Thread ", myid, ": are aware of the problem. Apologies for any inconvenience."

  call aims_stop('Ill-conditioned overlap matrix found. There may be a separate error file with details.','stop_illconditioning')

end subroutine stop_illconditioning
!******
