!****s* FHI-aims/mpi_init
!  NAME
!   mpi_init
!  SYNOPSIS
subroutine mpi_init(mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  implicit none

  integer :: mpierr

  mpierr = -1

end subroutine mpi_init
!******

!****s* FHI-aims/mpi_init
!  NAME
!   mpi_init
!  SYNOPSIS
subroutine mpi_initialized(flag,mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  implicit none
  logical :: flag
  integer :: mpierr

  mpierr = -1

end subroutine mpi_initialized
!******

!****s* FHI-aims/mpi_finalize
!  NAME
!   mpi_finalize
!  SYNOPSIS
subroutine mpi_finalize(mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  implicit none

  integer :: mpierr

  mpierr = 0

end subroutine mpi_finalize
!******

!****s* FHI-aims/mpi_comm_rank
!  NAME
!   mpi_comm_rank
!  SYNOPSIS
subroutine mpi_comm_rank(mpi_comm_global,myid, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  integer :: mpi_comm_global
  integer :: myid
  integer :: mpierr

  myid = 0
  mpierr = 0

end subroutine mpi_comm_rank
!******

!****s* FHI-aims/mpi_comm_size
!  NAME
!   mpi_comm_size
!  SYNOPSIS
subroutine mpi_comm_size(mpi_comm_global, n_tasks, mpierr)

  implicit none

  integer :: mpi_comm_global
  integer :: n_tasks
  integer :: mpierr

  n_tasks = 1
  mpierr = 0

end subroutine mpi_comm_size
!******

!****s* FHI-aims/mpi_comm_free
!  NAME
!   mpi_comm_free
!  SYNOPSIS
subroutine mpi_comm_free(comm, mpierr)

  use localorb_io, only: use_unit
  implicit none

  integer :: comm
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_ALLREDUCE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_comm_free
!******

!****s* FHI-aims/mpi_allreduce
!  NAME
!   mpi_allreduce
!  SYNOPSIS
subroutine mpi_allreduce(source, dest, size_mpi, type, method, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi

  real*8 :: source(size_mpi)
  real*8 :: dest(size_mpi)

  integer :: type
  integer :: method
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_ALLREDUCE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_allreduce
!******

!****s* FHI-aims/mpi_bcast
!  NAME
!   mpi_bcast
!  SYNOPSIS
subroutine mpi_bcast(source_and_dest, size_mpi, type, id, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi
  real*8 :: source_and_dest(size_mpi)

  integer :: type
  integer :: id
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_BCAST call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_bcast
!******

!****s* FHI-aims/mpi_barrier
!  NAME
!   mpi_barrier
!  SYNOPSIS
subroutine mpi_barrier(mpi_comm_global, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  implicit none

  integer :: mpi_comm_global
  integer :: mpierr

  mpierr = 0

end subroutine mpi_barrier
!******

!****s* FHI-aims/mpi_recv
!  NAME
!   mpi_recv
!  SYNOPSIS
subroutine mpi_recv(source_and_dest, size_mpi, type, id, tag, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi
  real*8 :: source_and_dest(size_mpi)

  integer :: type
  integer :: id
  integer :: tag
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_RECV call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_recv
!******

!****s* FHI-aims/mpi_send
!  NAME
!   mpi_send
!  SYNOPSIS
subroutine mpi_send(source_and_dest, size_mpi, type, id, tag, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi
  real*8 :: source_and_dest(size_mpi)

  integer :: type
  integer :: id
  integer :: tag
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_SEND call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_send
!******

!****s* FHI-aims/mpi_sendrecv
!  NAME
!   mpi_sendrecv
!  SYNOPSIS
subroutine  mpi_sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, &
     recvbuf, recvcount, recvtype, source, recvtag, comm, status, ierror)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: sendbuf, recvbuf
  integer :: sendcount, sendtype, dest, sendtag, recvcount, recvtype
  integer :: source, recvtag, comm, status(*), ierror

  write(use_unit,*) "This is a stub MPI_SENDRECV call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_sendrecv
!******

!****s* FHI-aims/mpi_sendrecv
!  NAME
!   mpi_sendrecv
!  SYNOPSIS
subroutine  mpi_sendrecv_replace(sendbuf, sendcount, sendtype, dest, sendtag, &
     source, recvtag, comm, status, ierror)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: sendbuf, recvbuf
  integer :: sendcount, sendtype, dest, sendtag, recvcount, recvtype
  integer :: source, recvtag, comm, status(*), ierror

  write(use_unit,*) "This is a stub MPI_SENDRECV call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_sendrecv_replace
!******
!****s* FHI-aims/mpi_alltoall
!  NAME
!   mpi_alltoall
!  SYNOPSIS
subroutine  mpi_alltoall(sendbuf, sendcount, sendtype, &
                         recvbuf, recvcount, recvtype, comm, ierror)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: sendbuf, recvbuf
  integer :: sendcount, sendtype, recvcount, recvtype
  integer :: comm, ierror

  write(use_unit,*) "This is a stub MPI_ALLTOALL call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_alltoall
!******
!****s* FHI-aims/mpi_alltoallv
!  NAME
!   mpi_alltoallv
!  SYNOPSIS
subroutine  mpi_alltoallv(sendbuf, sendcount, sdispl, sendtype, &
                          recvbuf, recvcount, rdispl, recvtype, comm, ierror)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: sendcount, sendtype, recvcount, recvtype
  integer :: sendbuf(1), recvbuf(1)
  integer :: sdispl(1), rdispl(1)
  integer :: comm, ierror

  write(use_unit,*) "This is a stub MPI_ALLTOALLV call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_alltoallv
!******

!****s* FHI-aims/mpi_abort
!  NAME
!   mpi_abort
!  SYNOPSIS
subroutine  mpi_abort(comm, errorcode, ierror)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: comm, errorcode, ierror
  write(use_unit,*) "This is a stub MPI_ABORT call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_abort
!******

!****s* FHI-aims/mpi_graph_create
!  NAME
!   mpi_graph_create
!  SYNOPSIS
subroutine mpi_graph_create( old_comm, nn, index, &
     edges, reorder, new_comm, mpierr )
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: old_comm, nn, new_comm, mpierr
  integer :: index, edges
  logical :: reorder

  write(use_unit,*) "This is a stub MPI_GRAPH_CREATE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_graph_create
!******

!****s* FHI-aims/mpi_get_processor_name
!  NAME
!   mpi_get_processor_name
!  SYNOPSIS
subroutine mpi_get_processor_name( name, length, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  character(*) :: name
  integer :: length
  integer :: mpierr

  name = 'non-mpi-localhost'
  length = 17
  mpierr = 0

end subroutine mpi_get_processor_name
!******

!****s* FHI-aims/mpi_comm_split
!  NAME
!   mpi_comm_split
!  SYNOPSIS
subroutine mpi_comm_split(mpi_comm_in, split_key, myid, my_comm, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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

  integer :: mpi_comm_in
  integer :: split_key
  integer :: myid
  integer :: my_comm
  integer :: mpierr

  myid = 0
  my_comm = mpi_comm_in

  mpierr = 0

end subroutine mpi_comm_split
!******

!****s* FHI-aims/mpi_comm_group
!  NAME
!   mpi_comm_group
!  SYNOPSIS
subroutine mpi_comm_group(mpi_comm_in, mpi_group, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: mpi_comm_in
  integer :: mpi_group
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_COMM_GROUP call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_comm_group
!******

!****s* FHI-aims/mpi_group_incl
!  NAME
!   mpi_group_incl
!  SYNOPSIS
subroutine mpi_group_incl(mpi_group_old, n, ranks, mpi_group_new, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: n
  integer, dimension(n) :: ranks
  integer :: mpi_group_old, mpi_group_new
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_COMM_INCL call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_group_incl
!******

!****s* FHI-aims/mpi_group_range_incl
!  NAME
!   mpi_group_range_incl
!  SYNOPSIS
subroutine mpi_group_range_incl(mpi_group_old, n, ranges, mpi_group_new, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: n
  integer, dimension(n,3) :: ranges
  integer :: mpi_group_old, mpi_group_new
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_GROUP_RANGE_INCL call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_group_range_incl
!******

!****s* FHI-aims/mpi_mod_comm_create
!  NAME
!   mpi_mod_comm_create
!  SYNOPSIS
subroutine mpi_comm_create(mpi_comm, mpi_group, mpi_comm_new, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: mpi_comm, mpi_comm_new
  integer :: mpi_group
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_COMM_CREATE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_comm_create
!******

!****s* FHI-aims/mpi_waitall
!  NAME
!   mpi_waitall
!  SYNOPSIS
subroutine mpi_waitall(num_reqs, mpi_requests, mpi_statuses, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none
  integer :: mpierr, num_reqs
  integer, dimension(num_reqs) :: mpi_requests
  integer, dimension(num_reqs,num_reqs) :: mpi_statuses

  write(use_unit,*) "This is a stub MPI_WAITALL call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_waitall
!******

!****s* FHI-aims/mpi_isend
!  NAME
!   mpi_isend
!  SYNOPSIS
subroutine mpi_isend(source_and_dest, size_mpi, type, id, tag, &
     comm, request, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi
  real*8 :: source_and_dest(size_mpi)

  integer :: type
  integer :: id
  integer :: tag
  integer :: comm
  integer :: err
  integer :: request

  write(use_unit,*) "This is a stub MPI_ISEND call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_isend
!******

!****s* FHI-aims/mpi_irecv
!  NAME
!   mpi_irecv
!  SYNOPSIS
subroutine mpi_irecv(source_and_dest, size_mpi, type, id, tag, &
     comm, request, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi
  real*8 :: source_and_dest(size_mpi)

  integer :: type
  integer :: id
  integer :: tag
  integer :: comm
  integer :: err
  integer :: request

  write(use_unit,*) "This is a stub MPI_IRECV call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_irecv
!******
!****s* FHI-aims/mpi_get_count
!  NAME
!   mpi_get_count
!  SYNOPSIS
subroutine mpi_get_count(mpi_status, type, count, mpierr)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
!    Release version, FHI-aims (2011).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: mpi_status(*)
  integer :: type
  integer :: count
  integer :: mpierr

  write(use_unit,*) "This is a stub MPI_GET_COUNT call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_get_count
!******
!****s* FHI-aims/mpi_wait
!  NAME
!   mpi_wait
!  SYNOPSIS
!call MPI_Wait(mpi_recv_requests(current_delta_v_recv_req),mpi_status,mpierr)
subroutine mpi_wait(request, mpi_status, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: mpi_status(*)
  integer :: err
  integer :: request

  write(use_unit,*) "This is a stub MPI_WAIT call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_wait
!******



!****f* FHI-aims/mpi_wtime
!  NAME
!   mpi_wtime
!  SYNOPSIS
real*8 function mpi_wtime()
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use timing, only: get_timestamps
  implicit none

  real*8 :: cpu_seconds, clock_seconds

  call get_timestamps(cpu_seconds, clock_seconds)
  mpi_wtime = clock_seconds

end function mpi_wtime
!******

!****s* FHI-aims/mpi_gather
!  NAME
!   mpi_gather
!  SYNOPSIS
subroutine mpi_gather( sndbuf, scount, datatype, recvbuf, rcount, rdatatype, &
                       root, comm, err &
                     )
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  logical :: sndbuf
  integer :: scount
  integer :: datatype
  logical :: recvbuf
  integer :: rcount
  integer :: rdatatype
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_gather call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_gather
!******

!****s* FHI-aims/mpi_allgather
!  NAME
!   mpi_allgather
!  SYNOPSIS
subroutine mpi_allgather( sndbuf, scount, datatype, recvbuf, rcount, rdatatype, &
                       comm, err &
                     )
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  logical :: sndbuf
  integer :: scount
  integer :: datatype
  logical :: recvbuf
  integer :: rcount
  integer :: rdatatype
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_gather call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_allgather
!******

!****s* FHI-aims/mpi_gatherv
!  NAME
!   mpi_gatherv
!  SYNOPSIS
subroutine mpi_gatherv( sndbuf, scount, datatype, recvbuf, rcounts, displs, rdatatype, &
                       root, comm, err &
                     )
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  logical :: sndbuf
  integer :: scount
  integer :: datatype
  logical :: recvbuf
  integer :: rcounts
  integer :: displs
  integer :: rdatatype
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_gatherv call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_gatherv
!******

!****s* FHI-aims/mpi_reduce
!  NAME
!   mpi_reduce
!  SYNOPSIS
subroutine mpi_reduce(source, dest, size_mpi, type, method, root, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: size_mpi

  real*8 :: source(size_mpi)
  real*8 :: dest(size_mpi)

  integer :: type
  integer :: method
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_REDUCE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_reduce
!******

!****s* FHI-aims/mpi_scatterv
!  NAME
!   mpi_scatterv
!  SYNOPSIS
subroutine mpi_scatter(sendbuf, sendcnts, sendtype, recvbuf, &
     recvcnt, recvtype, root, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: sendbuf
  integer :: sendcnts
  integer :: sendtype
  integer :: recvbuf
  integer :: recvcnt
  integer :: recvtype
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_REDUCE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_scatter
!******

!****s* FHI-aims/mpi_scatterv
!  NAME
!   mpi_scatterv
!  SYNOPSIS
subroutine mpi_scatterv(sendbuf, sendcnts, displs, sendtype, recvbuf, &
     recvcnt, recvtype, root, comm, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: sendbuf
  integer :: sendcnts
  integer :: displs
  integer :: sendtype
  integer :: recvbuf
  integer :: recvcnt
  integer :: recvtype
  integer :: root
  integer :: comm
  integer :: err

  write(use_unit,*) "This is a stub MPI_REDUCE call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_scatterv
!******

!****s* FHI-aims/mpi_op_create
!  NAME
!   mpi_op_create
!  SYNOPSIS
subroutine mpi_op_create(func, commute, op, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  external :: func
  logical :: commute
  integer :: op, err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_op_create
!******
!****s* FHI-aims/mpi_op_free
!  NAME
!   mpi_op_free
!  SYNOPSIS
subroutine mpi_op_free(op, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: op, err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_op_free
!******

!****s* FHI-aims/mpi_type_create_f90_integer
!  NAME
!   mpi_type_create_f90_integer
!  SYNOPSIS
subroutine mpi_type_create_f90_integer(r, newtype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: r
  integer :: newtype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_create_f90_integer
!******

!****s* FHI-aims/mpi_type_create_f90_real
!  NAME
!   mpi_type_create_f90_real
!  SYNOPSIS
subroutine mpi_type_create_f90_real(p, r, newtype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: p, r
  integer :: newtype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_create_f90_real
!******

!****s* FHI-aims/mpi_type_create_f90_complex
!  NAME
!   mpi_type_create_f90_complex
!  SYNOPSIS
subroutine mpi_type_create_f90_complex(p, r, newtype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: p, r
  integer :: newtype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_create_f90_complex
!******

!****s* FHI-aims/mpi_get_address
!  NAME
!   mpi_get_address
!  SYNOPSIS
subroutine mpi_get_address(location, address, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: location(*)
  integer :: address
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_get_address
!******

!****s* FHI-aims/mpi_type_create_struct
!  NAME
!   mpi_type_create_struct
!  SYNOPSIS
subroutine mpi_type_create_struct(count, array_of_blocklengths, &
                  array_of_displacements, array_of_types, newtype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: count
  integer :: array_of_blocklengths(*)
  integer :: array_of_displacements(*)
  integer :: array_of_types(*)
  integer :: newtype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_create_struct
!******

!****s* FHI-aims/mpi_type_create_resized
!  NAME
!   mpi_type_create_resized
!  SYNOPSIS
subroutine mpi_type_create_resized(oldtype, lb, extent, newtype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: oldtype, lb, extent, newtype, err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_create_resized
!******

!****s* FHI-aims/mpi_type_commit
!  NAME
!   mpi_type_commit
!  SYNOPSIS
subroutine mpi_type_commit(datatype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: datatype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_commit
!******

!****s* FHI-aims/mpi_type_free
!  NAME
!   mpi_type_free
!  SYNOPSIS
subroutine mpi_type_free(datatype, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Development version, FHI-aims (2014).
!  SOURCE
  use localorb_io, only: use_unit
  implicit none

  integer :: datatype
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_free
!******

!****s* FHI-aims/mpi_scatterv
!  NAME
!   mpi_scatterv
!  SYNOPSIS
subroutine mpi_test(request, flag, mpi_status, err)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: request
  logical :: flag
  integer :: mpi_status(*)
  integer :: err

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_test
!******

!****s* FHI-aims/mpi_get
!  NAME
!   mpi_get
!  SYNOPSIS
subroutine mpi_get( ORIGIN_ADDR, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK, TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE, WIN, IERROR )
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  complex*16 :: ORIGIN_ADDR
  INTEGER    :: TARGET_DISP, ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK, TARGET_COUNT, TARGET_DATATYPE, WIN, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_get
!******

!****s* FHI-aims/mpi_win_fence
!  NAME
!   mpi_win_fence
!  SYNOPSIS
subroutine MPI_WIN_FENCE(ASSERT, WIN, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  INTEGER ASSERT, WIN, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_win_fence
!******

!****s* FHI-aims/mpi_win_create
!  NAME
!   mpi_win_create
!  SYNOPSIS
subroutine MPI_WIN_CREATE(BASE, SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  complex*16 :: Base(*)
  INTEGER :: SIZE, DISP_UNIT, INFO, COMM, WIN, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_win_create
!******

!****s* FHI-aims/mpi_cart_sub
!  NAME
!   mpi_cart_sub
!  SYNOPSIS
subroutine MPI_CART_SUB(COMM, REMAIN_DIMS, COMM_NEW, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  LOGICAL :: REMAIN_DIMS(*)
  INTEGER :: COMM, COMM_NEW, IERROR


  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_cart_sub
!******


!****s* FHI-aims/mpi_cart_create
!  NAME
!   mpi_cart_create
!  SYNOPSIS
subroutine MPI_CART_CREATE(COMM_OLD, NDIMS, DIMS, PERIODS, REORDER, COMM_CART, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  LOGICAL :: PERIODS(*), REORDER
  INTEGER :: COMM_OLD, NDIMS, DIMS(*), COMM_CART, IERROR


  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_cart_create
!******

!****s* FHI-aims/mpi_type_size
!  NAME
!   mpi_type_size
!  SYNOPSIS
subroutine mpi_type_size(DATATYPE, SIZE, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  INTEGER :: DATATYPE, SIZE, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_size
!******

!****s* FHI-aims/mpi_win_free
!  NAME
!   mpi_win_free
!  SYNOPSIS
subroutine mpi_win_free(WIN, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  INTEGER :: WIN, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_win_free
!******

!****s* FHI-aims/mpi_comm_dup
!  NAME
!   mpi_comm_dup
!  SYNOPSIS
subroutine mpi_comm_dup(COMM, NEWCOMM, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  INTEGER :: COMM, NEWCOMM, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_comm_dup
!******

!****s* FHI-aims/mpi_allgatherv
!  NAME
!   mpi_allgatherv
!  SYNOPSIS
subroutine mpi_allgatherv(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, &
        RECVCOUNT, DISPLS, RECVTYPE, COMM, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: SENDBUF(*) ! I don't actually know if this is an integer or not, but code will be stopping shortly anyway
  complex*16 :: RECVBUF(*)
  INTEGER :: SENDCOUNT, SENDTYPE, RECVCOUNT(*)
  INTEGER :: DISPLS(*), RECVTYPE, COMM, IERROR

  write(use_unit,*) "This is a stub MPI_test call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_allgatherv
!******

!****s* FHI-aims/mpi_exscan
!  NAME
!   mpi_exscan
!  SYNOPSIS
subroutine mpi_exscan(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: SENDBUF(*)
  complex*16 :: RECVBUF(*)
  INTEGER :: COUNT, DATATYPE, OP, COMM, IERROR

  write(use_unit,*) "This is a stub MPI_Exscan call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_exscan
!******

!****s* FHI-aims/mpi_file_open
!  NAME
!   mpi_file_open
!  SYNOPSIS
subroutine mpi_file_open(COMM, FILENAME, MODE, INFO, FILEHANDLE, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  character(*) :: FILENAME
  integer :: COMM, MODE, INFO, FILEHANDLE, IERROR

  write(use_unit,*) "This is a stub MPI_File_open call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_open
!******

!****s* FHI-aims/mpi_file_close
!  NAME
!   mpi_file_close
!  SYNOPSIS
subroutine mpi_file_close(FILEHANDLE, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: FILEHANDLE, IERROR

  write(use_unit,*) "This is a stub MPI_File_close call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_close
!******

!****s* FHI-aims/mpi_file_read_at
!  NAME
!   mpi_file_read_at
!  SYNOPSIS
subroutine mpi_file_read_at(FILEHANDLE, OFFSET, BUF, COUNT, DATATYPE, STATUS, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: FILEHANDLE, BUF(*), COUNT, DATATYPE, STATUS(*), IERROR
  integer*8 :: OFFSET

  write(use_unit,*) "This is a stub MPI_File_read_at call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_read_at
!******

!****s* FHI-aims/mpi_file_read_at_all
!  NAME
!   mpi_file_read_at_all
!  SYNOPSIS
subroutine mpi_file_read_at_all(FILEHANDLE, OFFSET, BUF, COUNT, DATATYPE, STATUS, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: FILEHANDLE, BUF(*), COUNT, DATATYPE, STATUS(*), IERROR
  integer*8 :: OFFSET

  write(use_unit,*) "This is a stub MPI_File_read_at_all call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_read_at_all
!******

!****s* FHI-aims/mpi_file_write_at
!  NAME
!   mpi_file_write_at
!  SYNOPSIS
subroutine mpi_file_write_at(FILEHANDLE, OFFSET, BUF, COUNT, DATATYPE, STATUS, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: FILEHANDLE, BUF(*), COUNT, DATATYPE, STATUS(*), IERROR
  integer*8 :: OFFSET

  write(use_unit,*) "This is a stub MPI_File_write_at call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_write_at
!******

!****s* FHI-aims/mpi_file_write_at_all
!  NAME
!   mpi_file_write_at_all
!  SYNOPSIS
subroutine mpi_file_write_at_all(FILEHANDLE, OFFSET, BUF, COUNT, DATATYPE, STATUS, IERROR)
! PURPOSE
!   This is an MPI stub, i.e. a place holder. It ensures that AIMS can be
!   compiled on a single CPU computer without an external MPI installation.
!   As such it does nothing (except maybe to end the program, showing to a
!   programmer that something went wrong). Please refer to your MPI manual
!   or online documentation to find out about the arguments or what this
!   routine should actually be doing.
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
  use localorb_io, only: use_unit
  implicit none

  integer :: FILEHANDLE, BUF(*), COUNT, DATATYPE, STATUS(*), IERROR
  integer*8 :: OFFSET

  write(use_unit,*) "This is a stub MPI_File_write_at_all call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_write_at_all
!******
!****s* FHI-aims/mpi_file_get_size
!  NAME
!   mpi_file_get_size
!  SYNOPSIS
subroutine mpi_file_get_size(FH, SIZE, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: FH, IERROR
  integer*8 :: SIZE

  write(use_unit,*) "This is a stub MPI_File_get_size call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_get_size
!******
!****s* FHI-aims/mpi_file_set_size
!  NAME
!   mpi_file_set_size
!  SYNOPSIS
subroutine mpi_file_set_size(FH, SIZE, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: FH, IERROR
  integer*8 :: SIZE

  write(use_unit,*) "This is a stub MPI_File_set_size call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_set_size
!******
!****s* FHI-aims/mpi_file_set_view
!  NAME
!   mpi_file_set_view
!  SYNOPSIS
subroutine mpi_file_set_view(FH, DISP, ETYPE, FILETYPE, DATAREP, INFO, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: FH, ETYPE, FILETYPE, INFO, IERROR
  CHARACTER(*) :: DATAREP
  integer*8 :: DISP

  write(use_unit,*) "This is a stub MPI_File_set_view call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_set_view
!******
!****s* FHI-aims/mpi_error_string
!  NAME
!   mpi_error_string
!  SYNOPSIS
subroutine mpi_error_string(ERRORCODE, STRING, RESULTLEN, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: ERRORCODE, RESULTLEN, IERROR
  CHARACTER(*) :: STRING

  write(use_unit,*) "This is a stub MPI_Error_string call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_error_string
!******
!****s* FHI-aims/mpi_type_struct
!  NAME
!   mpi_type_struct
!  SYNOPSIS
subroutine mpi_type_struct(COUNT, ARRAY_OF_BLOCKLENGTHS, ARRAY_OF_DISPLACEMENTS, &
    ARRAY_OF_TYPES, NEWTYPE, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: COUNT, NEWTYPE, IERROR
  integer :: ARRAY_OF_BLOCKLENGTHS(*), ARRAY_OF_DISPLACEMENTS(*), ARRAY_OF_TYPES(*)

  write(use_unit,*) "This is a stub MPI_Type_struct call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_struct
!******
!****s* FHI-aims/mpi_type_extent
!  NAME
!   mpi_type_extent
!  SYNOPSIS
subroutine mpi_type_extent(DATATYPE, EXTENT, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: DATATYPE, EXTENT, IERROR

  write(use_unit,*) "This is a stub MPI_Type_extent call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_type_extent
!******
!****s* FHI-aims/mpi_file_read_all
!  NAME
!   mpi_file_read_all
!  SYNOPSIS
subroutine mpi_file_read_all(FH, BUF, COUNT, DATATYPE, STATUS, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: BUF(*)
  integer :: FH, COUNT, DATATYPE, STATUS(5), IERROR

  write(use_unit,*) "This is a stub MPI_File_read_all call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_read_all
!******
!****s* FHI-aims/mpi_file_write
!  NAME
!   mpi_file_write
!  SYNOPSIS
subroutine mpi_file_write(FH, BUF, COUNT, DATATYPE, STATUS, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: BUF(*)
  integer :: FH, COUNT, DATATYPE, STATUS(5), IERROR

  write(use_unit,*) "This is a stub MPI_File_write call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_file_write
!******
!****s* FHI-aims/mpi_iallgather
!  NAME
!   mpi_iallgather
!  SYNOPSIS
subroutine mpi_iallgather(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
    RECVTYPE, COMM, REQUEST, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: SENDBUF(*), RECVBUF(*)
  integer :: SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM, REQUEST, IERROR

  write(use_unit,*) "This is a stub MPI_Iallgather call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_iallgather
!******
!****s* FHI-aims/mpi_iallgatherv
!  NAME
!   mpi_iallgatherv
!  SYNOPSIS
subroutine mpi_iallgatherv(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
    DISPLS, RECVTYPE, COMM, REQUEST, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: SENDBUF(*), RECVBUF(*), RECVCOUNT(*), DISPLS(*)
  integer :: SENDCOUNT, SENDTYPE, RECVTYPE, COMM, REQUEST, IERROR

  write(use_unit,*) "This is a stub MPI_Iallgatherv call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_iallgatherv
!******

!******
!****s* FHI-aims/mpi_pack
!  NAME
!   mpi_pack
!  SYNOPSIS
subroutine mpi_pack(inbuf, incount, datatype, outbuf, outsize, &
    positions, COMM, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: inbuf(*), outbuf(*),  positions(*)
  integer :: incount, datatype, outsize, COMM, IERROR

  write(use_unit,*) "This is a stub MPI_pack call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_pack
!******

!******
!****s* FHI-aims/mpi_unpack
!  NAME
!   mpi_pack
!  SYNOPSIS
subroutine mpi_unpack(inbuf, insize, positions, datatype, outbuf, outcount, &
    COMM, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: inbuf(*), outbuf(*),  positions(*)
  integer :: insize, datatype, outcount, COMM, IERROR

  write(use_unit,*) "This is a stub MPI_pack call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_unpack
!******

!******
!****s* FHI-aims/mpi_dims_create
!  NAME
!   mpi_pack
!  SYNOPSIS
subroutine mpi_dims_create(nnodes, ndims, dims, IERROR)
  use localorb_io, only: use_unit
  implicit none

  integer :: dims(*)
  integer :: nnodes, ndims, IERROR

  write(use_unit,*) "This is a stub MPI_dims_create call."
  write(use_unit,*) "You should not be here. Check your"
  write(use_unit,*) "setting on the use_mpi variable, or"
  write(use_unit,*) "link against true MPI-libraries."
  stop

end subroutine mpi_dims_create
!******
