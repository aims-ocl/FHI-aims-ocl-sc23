!****s* FHI-aims/RRS-PBC/debug_rrs_pbc()
!  NAME
!    debug_rrs_pbc() 
!  SYNOPSIS

    subroutine debug_rrs_pbc()

!  PURPOSE
!  This routine implement the actual RRS-PBC calculation, and can only be called
!  after parse_rrs_pbc() and scf_solver()
!
!  USES

      use localorb_io
      use dimensions
      use physics
      use pbc_lists
      use runtime_choices
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
    integer       :: i_cell_n


    do i_cell_n = 1, n_cells
        write(rrs_pbc_info,'(I3,3F8.6)') i_cell_n, cell_index(i_cell_n,1), cell_index(i_cell_n,2), cell_index(i_cell_n,3)
        write(rrs_pbc_info,'(3F8.6)') k_point_list(i_cell_n,1), k_point_list(i_cell_n,2), k_point_list(i_cell_n,3)
        call localorb_info(rrs_pbc_info,use_unit,'(A)')

    end subroutine

