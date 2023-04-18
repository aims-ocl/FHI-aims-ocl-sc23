!****s* FHI-aims/RRS-PBC/output_rrs_pbc_matrix_2D
!  NAME
!   output_rrs_pbc_matrix_2D
!  SYNOPSIS

    subroutine output_rrs_pbc_matrix_2D(dim_matrix,p_matrix)
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

    integer :: dim_matrix
    complex*16  :: p_matrix(dim_matrix,dim_matrix)
    integer :: i, j, o, p
    integer :: i_center
    integer :: start_cell
    integer :: n
    integer :: m
    integer :: i_code
    character*130 :: rrs_pbc_info
    character*30 desc_str
    character(*), parameter :: func = 'output_rrs_pbc_center_matrix'


    start_cell = 1
    m = mod(dim_matrix,5)
    n = (dim_matrix-m)/5
    if (n > 0) then
        do i=1, n, 1
            write(rrs_pbc_info,'(4X,A,6X,5(I16))') '| ',&
                start_cell, start_cell+1,start_cell+2, &
                start_cell+3, start_cell +4
            call localorb_info(rrs_pbc_info)
            do o=1,dim_matrix,1
                write(rrs_pbc_info,'(4X,A,I5,A,5(F16.8))') '| ',&
                    o,'R',(real(p_matrix(j,o)), &
                    j=start_cell,start_cell+4)
                call localorb_info(rrs_pbc_info)
                !write(rrs_pbc_info,'(4X,A,I5,A,5(F16.8))') '| ',&
                !    o,'I',(aimag(p_matrix(j,o)), &
                !    j=start_cell,start_cell+4)
                !call localorb_info(rrs_pbc_info)
            enddo
            start_cell = start_cell + 5
        enddo
    endif
    if (m .eq. 1) then
        write(rrs_pbc_info,'(4X,A,6X,I16)') '| ',&
            start_cell
        call localorb_info(rrs_pbc_info)
        do o=1,dim_matrix,1
            write(rrs_pbc_info,'(4X,A,I5,A,F16.8)') '| ',&
                o,'R',real(p_matrix(start_cell,o))
            call localorb_info(rrs_pbc_info)
            !write(rrs_pbc_info,'(4X,A,I5,A,F16.8)') '| ',&
            !    o,'I',aimag(p_matrix(start_cell,o))
            !call localorb_info(rrs_pbc_info)
        enddo
    elseif (m .eq. 2) then
        write(rrs_pbc_info,'(4X,A,6X,2(I16))') '| ',&
            start_cell, start_cell+1
        call localorb_info(rrs_pbc_info)
        do o=1,dim_matrix,1
            write(rrs_pbc_info,'(4X,A,I5,A,2(F16.8))') '| ',&
                o,'R',(real(p_matrix(j,o)), &
                j=start_cell,start_cell+1)
            call localorb_info(rrs_pbc_info)
            !write(rrs_pbc_info,'(4X,A,I5,A,2(F16.8))') '| ',&
            !    o,'I',(aimag(p_matrix(j,o)), &
            !    j=start_cell,start_cell+1)
            !call localorb_info(rrs_pbc_info)
        enddo
    elseif (m .eq. 3) then
        write(rrs_pbc_info,'(4X,A,6X,3(I16))') '| ',&
            start_cell, start_cell+1,start_cell+2
        call localorb_info(rrs_pbc_info)
        do o=1,dim_matrix,1
            write(rrs_pbc_info,'(4X,A,I5,A,3(F16.8))') '| ',&
                o,'R',(real(p_matrix(j,o)), &
                j=start_cell,start_cell+2)
            call localorb_info(rrs_pbc_info)
            !write(rrs_pbc_info,'(4X,A,I5,A,3(F16.8))') '| ',&
            !    o,'I',(aimag(p_matrix(j,o)), &
            !    j=start_cell,start_cell+2)
            !call localorb_info(rrs_pbc_info)
        enddo
    elseif (m .eq. 4) then
        write(rrs_pbc_info,'(4X,A,6X,4(I16))') '| ',&
            start_cell, start_cell+1,start_cell+2, &
            start_cell+3
        call localorb_info(rrs_pbc_info)
        do o=1,dim_matrix,1
            write(rrs_pbc_info,'(4X,A,I5,A,4(F16.8))') '| ',&
                o,'R',(real(p_matrix(j,o)), &
                j=start_cell,start_cell+3)
            call localorb_info(rrs_pbc_info)
            !write(rrs_pbc_info,'(4X,A,I5,A,4(F16.8))') '| ',&
            !    o,'I',(aimag(p_matrix(j,o)), &
            !    j=start_cell,start_cell+3)
            !call localorb_info(rrs_pbc_info)
        enddo
    endif
    end subroutine

