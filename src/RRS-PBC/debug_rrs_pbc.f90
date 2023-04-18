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
      use physics
      use dimensions
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
    integer       :: i_cell_n, i_k_point


    if (allocated(cell_index)) then
        write(rrs_pbc_info,'(A)') 'no cell_index here'
        call localorb_info(rrs_pbc_info,use_unit,'(A)')
    else
        write(rrs_pbc_info,'(A)') 'cell_index still exists'
        call localorb_info(rrs_pbc_info,use_unit,'(A)')
    endif

    write(rrs_pbc_info,'(A)') 'number_of_super_cells still exists'
    call localorb_info(rrs_pbc_info,use_unit,'(A)')
    write(rrs_pbc_info,'(3I3)') number_of_super_cells
    call localorb_info(rrs_pbc_info,use_unit,'(A)')

    write(rrs_pbc_info,'(A,I3)') ' where n_cells = ', n_cells
    ! The last cell is not set in the above loop since it is a dummy cell
    call localorb_info(rrs_pbc_info,use_unit,'(A)')
    write(rrs_pbc_info,'(A,3F8.4)') 'In the K point of ',k_point_list(2,1),k_point_list(2,2),k_point_list(2,3)
    call localorb_info(rrs_pbc_info,use_unit,'(A)')
    do i_cell_n = 1, n_cells, 1
        !write(rrs_pbc_info,'(A,3F8.4)') '     cell index : ',cell_index(i_cell_n,1),cell_index(i_cell_n,1),cell_index(i_cell_n,1)
        !call localorb_info(rrs_pbc_info,use_unit,'(A)')
        write(rrs_pbc_info,'(A,2F8.4)') '        K phase : ',real(k_phase(i_cell_n,2)),aimag(k_phase(i_cell_n,2))
        call localorb_info(rrs_pbc_info,use_unit,'(A)')
    enddo
    write(rrs_pbc_info,'(A,3F8.4)') 'In the K point of ',k_point_list(1,1),k_point_list(1,2),k_point_list(1,3)
    call localorb_info(rrs_pbc_info,use_unit,'(A)')
    do i_cell_n = 1, n_cells, 1
        !write(rrs_pbc_info,'(A,3F8.4)') '     cell index : ',cell_index(i_cell_n,1),cell_index(i_cell_n,1),cell_index(i_cell_n,1)
        !call localorb_info(rrs_pbc_info,use_unit,'(A)')
        write(rrs_pbc_info,'(A,2F8.4)') '        K phase : ',real(k_phase(i_cell_n,1)),aimag(k_phase(i_cell_n,1))
        call localorb_info(rrs_pbc_info,use_unit,'(A)')
    enddo

    do i_k_point = 1, n_k_points, 1
        write(rrs_pbc_info,'(A,3F8.4)') 'The K point of ',&
            k_point_list(i_k_point,1), &
            k_point_list(i_k_point,2), &
            k_point_list(i_k_point,3)
        call localorb_info(rrs_pbc_info,use_unit,'(A)')
    enddo


    end subroutine

