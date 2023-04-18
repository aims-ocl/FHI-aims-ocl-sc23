!****s* FHI-aims/RRS-PBC/output_array
!  NAME
!   output_array
!  SYNOPSIS

    subroutine output_array(dim_array,p_array)
!  USES
    use dimensions
    use runtime_choices
    use localorb_io
!  PURPOSE
!   The subroutine prints out the array '''p_array'''
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

    integer :: dim_array
    real*8  :: p_array(dim_array)
    integer :: i, j, o, p
    integer :: i_center
    integer :: start_cell
    integer :: n
    integer :: m
    character*130 :: rrs_pbc_info
    character*30 desc_str
    character(*), parameter :: func = 'output_array'


    start_cell = 1
    m = mod(dim_array,5)
    n = (dim_array-m)/5
    if (n > 0) then
        do i=1, n, 1
            write(rrs_pbc_info,'(4X,A,5(I16))') '| ',&
                start_cell, start_cell+1,start_cell+2, &
                start_cell+3, start_cell +4
            call localorb_info(rrs_pbc_info)
            write(rrs_pbc_info,'(4X,A,5(F16.8))') '| ',&
                (p_array(j), &
                j=start_cell,start_cell+4)
            call localorb_info(rrs_pbc_info)
            !write(rrs_pbc_info,'(4X,A,I5,A,5(F16.8))') '| ',&
            !    o,'I',(aimag(p_array(j,o)), &
            !    j=start_cell,start_cell+4)
            !call localorb_info(rrs_pbc_info)
            start_cell = start_cell + 5
        enddo
    endif
    if (m .eq. 1) then
        write(rrs_pbc_info,'(4X,A,I16)') '| ',&
            start_cell
        call localorb_info(rrs_pbc_info)
        write(rrs_pbc_info,'(4X,A,F16.8)') '| ',&
            p_array(start_cell)
        call localorb_info(rrs_pbc_info)
        !write(rrs_pbc_info,'(4X,A,I5,A,F16.8)') '| ',&
        !    o,'I',aimag(p_array(start_cell,o))
        !call localorb_info(rrs_pbc_info)
    elseif (m .eq. 2) then
        write(rrs_pbc_info,'(4X,A,2(I16))') '| ',&
            start_cell, start_cell+1
        call localorb_info(rrs_pbc_info)
        write(rrs_pbc_info,'(4X,A,2(F16.8))') '| ',&
            (p_array(j), &
            j=start_cell,start_cell+1)
        call localorb_info(rrs_pbc_info)
        !write(rrs_pbc_info,'(4X,A,I5,A,2(F16.8))') '| ',&
        !    o,'I',(aimag(p_array(j,o)), &
        !    j=start_cell,start_cell+1)
        !call localorb_info(rrs_pbc_info)
    elseif (m .eq. 3) then
        write(rrs_pbc_info,'(4X,A,3(I16))') '| ',&
            start_cell, start_cell+1,start_cell+2
        call localorb_info(rrs_pbc_info)
        write(rrs_pbc_info,'(4X,A,3(F16.8))') '| ',&
            (p_array(j), &
            j=start_cell,start_cell+2)
        call localorb_info(rrs_pbc_info)
        !write(rrs_pbc_info,'(4X,A,I5,A,3(F16.8))') '| ',&
        !    o,'I',(aimag(p_array(j,o)), &
        !    j=start_cell,start_cell+2)
        !call localorb_info(rrs_pbc_info)
    elseif (m .eq. 4) then
        write(rrs_pbc_info,'(4X,A,4(I16))') '| ',&
            start_cell, start_cell+1,start_cell+2, &
            start_cell+3
        call localorb_info(rrs_pbc_info)
        write(rrs_pbc_info,'(4X,A,4(F16.8))') '| ',&
            (p_array(j), &
            j=start_cell,start_cell+3)
        call localorb_info(rrs_pbc_info)
        !write(rrs_pbc_info,'(4X,A,I5,A,4(F16.8))') '| ',&
        !    o,'I',(aimag(p_array(j,o)), &
        !    j=start_cell,start_cell+3)
        !call localorb_info(rrs_pbc_info)
    endif
    end subroutine output_array

