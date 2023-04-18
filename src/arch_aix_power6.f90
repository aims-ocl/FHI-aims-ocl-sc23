!****h* FHI-aims/arch_aix_power6
!  NAME
!     arch_aix_power6
!  SYNOPSIS

module arch_specific

  !  PURPOSE
  !
  ! This module contains any architecture specific 
  ! calls for the "aix_power6" case.
  !
  !  PLEASE no dependencies on external modules if they can avoided.
  !  check_environment below is a special case (only an interface definition,
  !  has no dependecies of its own)
  !
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
  !
  !******

  implicit none

contains
  !
  !-------------------------------------------------------------------------
  !****f* arch_aix_power6/arch_erf
  !  NAME
  !    arch_erf
  !  SYNOPSIS

  real*8 function arch_erf(argument)

    !  PURPOSE
    !  The function calls Erf(x)  function 
    !
    implicit none
    !  ARGUMENTS

    real*8 :: argument

    !  INPUTS
    !    o argument -- argiment of erf-function
    !  OUTPUT
    !    o arch_erf -- value of erf-function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE



    arch_erf = erf(argument)

    return
  end function arch_erf
    !******
  !-------------------------------------------------------------------------
  !****f* arch_aix_power6/arch_erfc
  !  NAME
  !    arch_erfc
  !  SYNOPSIS

  real*8 function arch_erfc(argument)

    !  PURPOSE
    !  The function calls Erfc(x)  function 
    !
    implicit none
    !  ARGUMENTS

    real*8 :: argument

    !  INPUTS
    !    o argument -- argiment of erfc-function
    !  OUTPUT
    !    o arch_erf -- value of erfc-function
    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2008).
    !  SOURCE


    arch_erfc = erfc(argument)

    return
  end function arch_erfc
  !******
  !-------------------------------------------------------------------------

  subroutine set_openmp_threads ( myid )
    use check_environment
    implicit none

    ! current MPI process ID, only for information 
    ! (we do not want dependencies on module mpi_tasks or similar here)
    integer, intent(IN) :: myid

    ! local variables

    integer :: omp_set_num_threads, err, stat
    integer :: mkl_Set_num_threads
    integer, parameter :: openmp_threads = 1
    character(LEN=100) :: omp_variable


    if (.not.(check_environment_variable("OMP_NUM_THREADS",openmp_threads,myid))) then
       ! OMP_NUM_THREADS is not set to openmp_threads. Do this now (using a
       ! OpenMP directive)

       err = omp_set_num_threads(openmp_threads)
  
       if (myid .eq. 0) then
          print *,' '
          print *,' Setting OMP_NUM_THREADS to 1 '
          print *,' This is a precautionary measure! '
          print *,' You might want to change this behaviour'
          print *,' '
       endif

    endif


!    if (.not.(check_environment_variable("MP_USE_BULK_XFER","yes",myid))) then
!       ! MP_USE_BULK_XFER is not set to yes. We can not set it here
!       ! so just an information is printed in the function
!       
!       ! if possible: set variable here
!    endif

  
    if (.not.(check_environment_variable("MP_BULK_MIN_MSG_SIZE",32768,myid))) then
       ! MP_BULK_MIN_MSG_SIZE is not set to 32768. We can not set it here
       ! so just an information is printed in the function
       
       ! if possible: set variable here
    endif  

  

  
    if (.not.(check_environment_variable("MP_EUILIB","yes",myid))) then
       ! MP_EUILIB is not set to yes. We can not set it here
       ! so just an information is printed in the function
       
       ! if possible: set variable here
    endif  

  
    if (.not.(check_environment_variable("MP_EUIDEVICE","sn_all",myid))) then
       ! MP_EUIDEVICE is not set to sn_all. We can not set it here
       ! so just an information is printed in the function
       
       ! if possible: set variable here
    endif  

    if (.not.(check_environment_variable("MP_SHARED_MEMORY","yes",myid))) then
       ! MP_SHARED_MEMORY is not set to yes. We can not set it here
       ! so just an information is printed in the function
       
       ! if possible: set variable here
    endif  

    if (.not.(check_environment_variable("MEMORY_AFFINITY","MCM",myid))) then
       ! MEMORY_AFFINITY is not set to MCM. We can not set it here
       ! so just an information is printed in the function
       
       ! if possible: set variable here
    endif  

  end subroutine set_openmp_threads


end module arch_specific
