!****s* FHI-aims/RRS-PBC/read_rrs_pbc_unit_cell
!  NAME
!   read_rrs_pbc_unit_cell
!  SYNOPSIS

    subroutine read_rrs_pbc_unit_cell(fn,iop)
!  USES
    use dimensions, only: rrs_pbc_center_atom, rrs_pbc_equal_atom, &
        rrs_pbc_n_equal, rrs_pbc_n_center_atom
    use mpi_tasks, only: aims_stop_coll
    implicit none
!  PURPOSE
!   The subroutine reads the center unit cell information from control.in for the 
!   RRS-PBC method.
!   This subroutine is called from read_control and dimension
!  INPUTS
!    fn   :: the file flow id of control.in
!    iop  :: 0 for actual cell_rrs_pbc loading
!            other value for check only
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

    integer :: fn, iop
    integer :: i, j, n, m
    integer :: start_cell
    integer :: tmp_equal
    integer :: i_code
    integer :: i_center
    character*130 :: inputline
    character*30 desc_str
    character*100 :: info_str
    character(*), parameter :: func = 'read_rrs_pbc_unit_cell'

    if (iop .eq. 0) then
        allocate(rrs_pbc_center_atom(3,rrs_pbc_n_center_atom))
        allocate(rrs_pbc_equal_atom(2,rrs_pbc_n_equal,rrs_pbc_n_center_atom))
    endif

    do i_center=1, rrs_pbc_n_center_atom, 1
        read(fn,'(A)',iostat=i_code) inputline
        if(i_code .ne. 0) then
           call aims_stop_coll("Unknown error reading file 'control.in'...", func)
        endif
        read(inputline,*,iostat=i_code) desc_str
        if(i_code/=0 .or. desc_str(1:1).eq.'#') then
           call aims_stop_coll( &
               "Error: empty and comment lines are forbidden following the keyword of cell_rrs_pbc" &
               , func)
        endif
        if (iop .eq. 0) then
            read(inputline,*) rrs_pbc_center_atom(1,i_center), tmp_equal
        else
            read(inputline,*) desc_str, tmp_equal
        endif

        start_cell = 1
        m = mod(tmp_equal,5)
        n = (tmp_equal-m)/5
        if (n > 0) then
            do i=1, n, 1
                read(fn,'(A)',iostat=i_code) inputline
                if(i_code .ne. 0) then
                   call aims_stop_coll("Unknown error reading file 'control.in'...", func)
                endif
                read(inputline,*,iostat=i_code) desc_str
                if(i_code/=0 .or. desc_str(1:1).eq.'#') then
                   call aims_stop_coll( &
                       "Error: empty and comment lines are forbidden following the keyword of cell_rrs_pbc" &
                       , func)
                endif

                if (iop.eq.0) then
                    read(inputline,*) (rrs_pbc_equal_atom(1,j,i_center), &
                            j=start_cell,start_cell+4)
                    start_cell = start_cell + 5
                endif
            enddo
        endif
        if (m >= 1) then
            read(fn,'(A)',iostat=i_code) inputline
            if(i_code .ne. 0) then
               call aims_stop_coll("Unknown error reading file 'control.in'...", func)
            endif
            read(inputline,*,iostat=i_code) desc_str
            if(i_code/=0 .or. desc_str(1:1).eq.'#') then
               call aims_stop_coll( &
                   "Error: empty and comment lines are forbidden following the keyword of cell_rrs_pbc" &
                   , func)
            endif
            
            if (iop.eq.0) then
                read(inputline,*) (rrs_pbc_equal_atom(1,j,i_center), &
                        j=start_cell,start_cell+m-1)
            endif
        endif
    enddo
    end subroutine

