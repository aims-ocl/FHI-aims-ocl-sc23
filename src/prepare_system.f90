!****s* FHI-aims/prepare_system
!  NAME
!   prepare_system
!  SYNOPSIS

  subroutine prepare_system

!  PURPOSE
!  Driver for calling subroutines which prepare the system for good
!  performance and possible checks for known errors. In addition, this 
!  driver checks for functionality support such as scalapack or spglib.
!  Supported external libraries should be checked here to activate their
!  functionality in FHI-aims. Library related variables, if not otherwise
!  specified should be put in dimensions.f90 followed by 
!  sanity checks in readcontrol.f90
!  
!  Please add more tests here:
!
!  Preperations done:
!  * set OMP_NUM_THREADS to 1 (also MKL_NUM_THREADS)
!
!  Tests done:
!  * Determine the settings of the stacksize
!  * check the pdtran routine from LAPACK/SCALAPACK to prevent the code
!    running with a buggy MKL routine
!  * Check if FHI-aims was compiled with spglib
!
!  USES
    use mpi_tasks
    use localorb_io
    use arch_specific
    use dimensions, only : use_spglib
    use spglib_symmetry, only : check_spglib
    

    implicit none

!  ARGUMENTS
!  INPUTS
!    none
!  OUTPUT
!    none
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
!    Release version, FHI-aims (20011).
!  SOURCE



    ! local variables
    integer :: i, mypnum, nprocs
   
    call localorb_info('')
    call localorb_info('  Performing system and environment tests:')

    call set_openmp_threads(myid)

    call check_stacksize()

    ! pdtran test is only necessary with MPI usage
    if (use_mpi) then
       call localorb_info('  | Checking for scalapack...')
       call BLACS_Pinfo(mypnum, nprocs)
       if (nprocs.gt.0 ) then
          call localorb_info('  | Testing pdtran()...')
          do i=1, min(n_tasks, 8)
             call test_pdtran(i)
          enddo
          call localorb_info('  | All pdtran() tests passed.')
       else
          call localorb_info('  | No scalapack compiled in.')
       end if
    endif

    ! Check if spglib is available
    call check_spglib(use_spglib)

  end subroutine prepare_system
