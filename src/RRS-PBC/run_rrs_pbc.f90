!****s* FHI-aims/RRS-PBC/run_rrs_pbc()
!  NAME
!    run_rrs_pbc() 
!  SYNOPSIS

    subroutine run_rrs_pbc()

!  PURPOSE
!  This routine implement the actual RRS-PBC calculation, and can only be called
!  after parse_rrs_pbc()
!
!  USES

      use localorb_io
      use constants
      use dimensions
      use physics
      use geometry
      use numerical_utilities
      use basis
      use runtime_choices
      use mpi_tasks
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
      character*132 :: func = 'parse_rrs_pbc()'
      integer :: i, j, o, p
      integer :: info
      real*8, dimension(3)   :: tmp_k_vec
      real*8, dimension(3)   :: tmp_d_vec


      !call purge_rrs_pbc_hamiltonian()

      call out_rrs_pbc_band()

      !call prepare_rrs_pbc_energy()

      end subroutine

