!****p* FHI-aims/aims_real
!  NAME
!    aims_real -- Program to call the top-level FHI-aims subroutine
!  SYNOPSIS

    program aims_real

!  PURPOSE
!
!    This is FHI-aims code
!
!    In principle, FHI-aims can also be built as a library. In that case, aims
!    is called as a subroutine aims() (see main.f90). The present program embeds
!    this subroutine for the standalone executable.
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  HISTORY
!    Release version, FHI-aims (2008).
!  INPUT
!    none
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!     Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement.
!  SOURCE

! if an MPI-environment exists, initialize it - else, silently
! circumvent all MPI functionality using dummy subroutines
      use mpi_utilities

! Add this variables to asses BLACS environment
      integer :: me, nprocs

      call initialize_mpi ( )

      call aims ( mpi_comm_global, 6 , use_mpi )

      ! AJL/Feb2016
      call BLACS_Pinfo(me,nprocs)
      if (nprocs.gt.0) then
          ! Transfered from scalapack_wrapper.f90 as not suitable inside a library setup
          call BLACS_Exit(1)
      end if

      call finalize_mpi ( )

    end program aims_real
!******
