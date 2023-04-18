
!****s* FHI-aims/RRS-PBC/rrs_pbc()
!  NAME
!    rrs_pbc() 
!  SYNOPSIS

    subroutine rrs_pbc()

!  PURPOSE
!  High-level wrapper around the post RRS-PBC projection depending on the
!  cluster calculation.
!  This routine can only be called after a succeful scf procedure. Necessary
!  conditions for running are:
!  * an converged overlap matrix
!  * an converged hamiltonian matrix,
!  and all of that with the correct array dimensions
!
!  USES

      use localorb_io
      use dimensions
      use physics
      implicit none

!  ARGUMENTS

!  INPUTS
!    none
!  OUTPUT
!    none
!  COPYRIGHT
!   
!   
!   
!  HISTORY
!   
!  SOURCE
!

      ! imported variables


      ! local variables
      character*132 :: rrs_pbc_info

      ! print start info
      call localorb_info('',use_unit)
      call localorb_info( &
      "  ------------------------------------------------------------", &
            6,'(A)')
      write(rrs_pbc_info,'(2X,A)') &
         'Start RRS-PBC calculation based on converged FCM data'
      call localorb_info(rrs_pbc_info)
      call localorb_info( &
      "  ------------------------------------------------------------", &
            6,'(A)')
      call localorb_info('',use_unit)


      ! run RRS-PBC calculation
      !call store_rrs_pbc_data()
      write(use_unit,'(4X,A)') "------------------------------------ "
      call out_rrs_pbc_band()

      !call prepare_rrs_pbc_energy()
      write(use_unit,'(4X,A)') "------------------------------------ "


      ! print end info
      call localorb_info( &
      "  ------------------------------------------------------------", &
            6,'(A)')
      write(rrs_pbc_info,'(2X,A)') &
          'End RRS-PBC calculation'
      call localorb_info(rrs_pbc_info)
      call localorb_info( &
      "  ------------------------------------------------------------", &
            6,'(A)')
      call localorb_info('',use_unit)
      end subroutine

