!****s* FHI-aims/noscf_solver
!  NAME
!    scf_solver
!  SYNOPSIS

    subroutine noscf_solver &
    ( converged, enough_walltime_left )

!  PURPOSE
!  High-level wrapper to avoid the scf calculation:
!  The total energy, the entropy and the forces are set
!  to 0.0d0
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
      logical :: enough_walltime_left

!  INPUTS
!    none
!  OUTPUT
!    o  converged -- always .true.
!    o  enough_walltime_left -- always .true. 
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
!  SOURCE
!

! local variables
      character*100 :: info_str

      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      write(info_str,'(A)') '************************** W A R N I N G *******************************'
      call localorb_info ( info_str )
      write(info_str,'(A)') '* Skipping the usual SCF cycles                                        *'
      call localorb_info ( info_str )
      write(info_str,'(A)') '************************************************************************'
      call localorb_info ( info_str )
      write(info_str,'(A)') ''
      call localorb_info ( info_str )
      total_energy       = 0.0d0
      entropy_correction = 0.0d0
      total_forces(:,:)  = 0.0d0
      converged     = .true.
      enough_walltime_left = .true.

    end subroutine noscf_solver
!******
