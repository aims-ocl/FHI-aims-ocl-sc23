!****s* FHI-aims/reinitialize_noscf
!  NAME
!    reinitialize_scf
!  SYNOPSIS

    subroutine reinitialize_noscf &
    ( converged )

!  PURPOSE
!  High-level wrapper that avoids all reinitialization of the SCF
!
!  USES

      use localorb_io
      use dimensions
      use timing
      use runtime_choices
      use physics

      implicit none

!  ARGUMENTS

      logical :: converged

!  INPUTS
!    none
!  OUTPUT
!    o converged -- always .true.
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
!

      ! local variables
      character*100 :: info_str

      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      write(info_str,'(A)') '************************** W A R N I N G *******************************'
      call localorb_info ( info_str )
      write(info_str,'(A)') '* Skipping the SCF reinitialization.                                     *'
      call localorb_info ( info_str )
      write(info_str,'(A)') '************************************************************************'
      call localorb_info ( info_str )
      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      call allocate_physics ( )
      converged = .true.

    end subroutine reinitialize_noscf
!******
