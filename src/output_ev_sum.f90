!****s* FHI-aims/output_ev_sum
!  NAME
!   output_ev_sum
!  SYNOPSIS

      subroutine output_ev_sum &
      ( ev_sum, ev_sum_shifted &
      )

!  PURPOSE
!  The subroutine writes out the sum of eigenvalues.
!  This works only for the non-spinpolarised, non-charged case.
!  This is useful for the non-selfconsistent case.
!
!  USES

      use dimensions
      use localorb_io
      use constants

!  ARGUMENTS

      real*8 :: ev_sum
      real*8 :: ev_sum_shifted

!  INPUTS
!   o ev_sum -- sum of eigenvalues
!   o ev_sum_shifted -- shifted sum of eigenvalues
!      
!  OUTPUT
!   none
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


	


!  local variables

      character*100 :: info_str

!  begin work

      write(info_str,*)
      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A)') "Eigenvalue sums :"

      call localorb_info(info_str,use_unit,'(A)')
      write(info_str,'(2X,A,F16.5,A)')  "| Total sum of eigenvalues                       : ", &
        ev_sum*hartree, " eV"
      call localorb_info(info_str,use_unit,'(A)')

!      write(info_str,'(2X,A,F16.5,A)')
!     +  "| Core-level shifted sum of eigenvalues          : ",
!     +  ev_sum_shifted*hartree, " eV"
!      call localorb_info(info_str,use_unit,'(A)')

      write(info_str,'(2X,A,F16.5,A)') &
        "| Sum of eigenvalues per atom                    : ", ev_sum*hartree/dble(n_atoms), " eV"
      call localorb_info(info_str,use_unit,'(A)')

!      write(info_str,'(2X,A,F16.5,A)')
!     +  "| Core-level shifted sum of eigenvalues per atom : ",
!     +  ev_sum_shifted*hartree/dble(n_atoms), " eV"
!      call localorb_info(info_str,use_unit,'(A)')

!  that's all folks

    end subroutine output_ev_sum
!******	
