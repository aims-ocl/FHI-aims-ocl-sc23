!****s* FHI-aims/aims_shm_implemented
!  NAME
!    aims_shm_implemented
!  SYNOPSIS

subroutine aims_shm_implemented(flag)

  !  PURPOSE
  !     Stub returning 0 to signify absence of shm code.
  !  USES

  implicit none

  !  ARGUMENTS

  integer, intent(OUT) :: flag

  !  INPUTS
  !    none
  !  OUTPUTS
  !    o flag -- set to 0 because no shm code is compiled in
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
  !    Release version, FHI-aims (2010).
  !  SOURCE

  character(*), parameter :: func = 'aims_shm_implemented'

  flag = 0

end subroutine aims_shm_implemented
!******
!****s* FHI-aims/aims_shm_release
!  NAME
!   aims_shm_release
!  SYNOPSIS
subroutine aims_shm_release
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_release
!******
!****s* FHI-aims/aims_shm_init
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_init(shm_size, info)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE
  use mpi_tasks, only: STDERR

  implicit none
  
  integer :: shm_size
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_init
!******

!****s* FHI-aims/aims_shm_init64
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_init64(shm_size, info)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  integer*8 :: shm_size
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_init64
!******

!****s* FHI-aims/aims_shm_n_procs
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_n_procs(n_procs)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  integer :: n_procs
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_n_procs
!******

!****s* FHI-aims/aims_shm_myid
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_myid(myid)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  integer :: myid
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_myid
!******

!****s* FHI-aims/aims_shm_get
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_get( data, offset, size)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  real*8 :: data
  integer :: offset
  integer :: size
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_get
!******

!****s* FHI-aims/aims_shm_put
!  NAME
!   aims_shm_init
!  SYNOPSIS
subroutine aims_shm_put( data, offset, size)
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  real*8 :: data
  integer :: offset
  integer :: size
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_put
!******

!****s* FHI-aims/aims_shm_get_ptr64
!  NAME
!   aims_shm_get_ptr64
!  SYNOPSIS
subroutine aims_shm_get_ptr64
! PURPOSE
!   This is a stub, i.e. a place holder. It ensures that AIMS can be 
!   compiled without a C compiler if shared memory support is not requested.
!   As such it does nothing (except maybe to end the program, showing to a 
!   programmer that something went wrong). 
!   If you ended up here, you likely mixed up some input options.
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
!  SOURCE

  use mpi_tasks, only: STDERR
  implicit none
  
  integer :: info
  
  info=1

  write(STDERR,*) "* Attention. Your run has called a subroutine intended for "
  write(STDERR,*) "* shared-memory based communication in the Hartree potential. "
  write(STDERR,*) "* This functionality must be compiled into the code using a "
  write(STDERR,*) "* subroutine written in C. The code version that you are using "
  write(STDERR,*) "* was not compiled with support for shared memory. "
  write(STDERR,*) "* "
  write(STDERR,*) "* Either modify your control.in file, or build/use the code version "
  write(STDERR,*) "* including shared memory support."
  write(STDERR,*) "* Stopping the code for now. "
  stop
      
end subroutine aims_shm_get_ptr64
!******
