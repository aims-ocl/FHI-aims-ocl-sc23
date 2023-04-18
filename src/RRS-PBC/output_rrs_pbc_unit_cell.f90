!****s* FHI-aims/RRS-PBC/output_rrs_pbc_unit_cell
!  NAME
!   output_rrs_pbc_unit_cell
!  SYNOPSIS

    subroutine output_rrs_pbc_unit_cell()
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
!  PURPOSE
!   The subroutine print out the cell element information of RRS-PBC method.
!   This subroutine is called from parse_rrs_pbc
!  INPUTS
!    none
!  OUTPUT
!    none
!  AUTHOR
!                                                                  
!  SEE ALSO
!    
!  COPYRIGHT
!   
!  HISTORY
!    
!  SOURCE

    integer :: i, j
    integer :: i_center
    integer :: start_cell
    integer :: n
    integer :: m
    integer :: i_code
    character*130 :: rrs_pbc_info
    character*30 desc_str
    character(*), parameter :: func = 'output_rrs_pbc_center_unit_cell'

    do i_center=1, rrs_pbc_n_center_atom, 1
        write(rrs_pbc_info,'(4X,A,I5,A,2(I5),A)') &
            '| Equal points associated with the center atom',&
            rrs_pbc_center_atom(1,i_center), ' (', &
            rrs_pbc_center_atom(2,i_center), &
            rrs_pbc_center_atom(3,i_center),'):'
        call localorb_info(rrs_pbc_info)
        start_cell = 1
        m = mod(rrs_pbc_n_equal,5)
        n = (rrs_pbc_n_equal-m)/5
        if (n > 0) then
            do i=1, n, 1
                write(rrs_pbc_info,'(4X,A,5(I8))') '| ',&
                    (rrs_pbc_equal_atom(1,j,i_center), &
                    j=start_cell,start_cell+4)
                start_cell = start_cell + 5
                call localorb_info(rrs_pbc_info)
            enddo
        endif
        if (m .eq. 1) then
            write(rrs_pbc_info,'(4X,A,I8)') '| ', &
                (rrs_pbc_equal_atom(1,j,i_center), &
                j=start_cell,start_cell+m-1)
            call localorb_info(rrs_pbc_info)
        elseif (m .eq. 2) then
            write(rrs_pbc_info,'(4X,A,2(I8))') '|',&
                (rrs_pbc_equal_atom(1,j,i_center), &
                j=start_cell,start_cell+m-1)
            call localorb_info(rrs_pbc_info)
        elseif (m .eq. 3) then
            write(rrs_pbc_info,'(4X,A,3(I8))') '|',&
                (rrs_pbc_equal_atom(1,j,i_center), &
                j=start_cell,start_cell+m-1)
            call localorb_info(rrs_pbc_info)
        elseif (m .eq. 4) then
            write(rrs_pbc_info,'(4X,A,4(I8))') '| ',&
                (rrs_pbc_equal_atom(1,j,i_center), &
                j=start_cell,start_cell+m-1)
            call localorb_info(rrs_pbc_info)
        endif
    enddo

    end subroutine

