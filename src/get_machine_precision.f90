!****s* FHI-aims/get_machine_precision
!  NAME
!   get_machine_precision
!  SYNOPSIS
      subroutine get_machine_precision ( safe_minimum )
!  PURPOSE
!    Subroutine get_macine_precision determines the machine precision
!  USES
      use mpi_tasks
      use localorb_io

      implicit none

!  ARGUMENTS
      real*8 :: safe_minimum
!  INPUTS
!    none
!  OUTPUTS
!    o safe_minimum -- What is smallest  tolerable number of system.
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

!  functions
!  (from lapack)

      double precision DLAMCH

!  begin work

      call localorb_info( &
           "Determining machine precision: ", &
           use_unit,'(2X,A)')

      safe_minimum = dlamch ('S')

      if (myid.eq.0) then
         write(use_unit,*) " ", safe_minimum
      end if

      if (safe_minimum.eq.0.d0) then
         if (myid.eq.0) then
            write(use_unit,*) "The BLAS subroutine dlamch() claims that the", &
                 " minimum tolerable number on your system is ", &
                 safe_minimum, "."
            write(use_unit,*) "A value of zero is factually not correct, ", &
          "and may throw off any lapack calls later on. "
            write(use_unit,*) "Please re-examine your compiled version of ", &
                 "BLAS before continuing."
            write(use_unit,*) "If in doubt, compile dlamch() separately ", &
                 "without optimisation flags."
         end if
         stop
      end if

      return
    end subroutine get_machine_precision
!******
