!****h* FHI-aims/mpi_routines
!  NAME
!    mpi_routines -- dummy module to resolve use statements in vibrations code
!  SYNOPSIS
module mpi_routines
!  PURPOSE
!    In the vibrations code, there is a dummy "use MPI_ROUTINES" statement
!    which is used as placeholder text, to be modified during the make process and
!    pointed at the *actual* MPI module.  This can confuse source code analysis
!    tools which are expecting an mpi_routines module to exist.
!
!    This, we've created this empty module to satisfy the use statement but have
!    deliberately kept it empty, i.e. not even any stub subroutines, to force the
!    linking stage to fail if this module is ever directly used.
!  USES
  implicit none
!  AUTHOR
!    William Huhn
!  HISTORY
!    June 2018 - Added.
!  COPYRIGHT
!    Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!    e.V. Please note that any use of the "FHI-aims-Software" is subject to
!    the terms and conditions of the respective license agreement."
end module mpi_routines
