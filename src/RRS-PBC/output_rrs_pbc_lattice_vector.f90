!****s* FHI-aims/RRS-PBC/output_rrs_pbc_lattice_vector
!  NAME
!   output_rrs_pbc_lattice_vector
!  SYNOPSIS

    subroutine output_rrs_pbc_lattice_vector()
!  USES
    use dimensions
    use geometry
    use runtime_choices
    use localorb_io
    use bravais
!  PURPOSE
!   The subroutine print out the lattice information of RRS-PBC method.
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

    character*130 :: rrs_pbc_info
    integer       :: i,j,n,m

    ! print out lattice vector
    write(rrs_pbc_info,'(4X,A)') &
        "| Lattice_vector                      :"
    call localorb_info( rrs_pbc_info )
    do i = 1, 3, 1
       write(rrs_pbc_info,'(4X,A,3(2X,F16.8),2X)') &
            "|", &
            (rrs_pbc_lattice_vector(j,i) * bohr, j=1,3,1)
       call localorb_info( rrs_pbc_info )
    enddo
    ! print out recip lattice vector
    write(rrs_pbc_info,'(4X,A)') &
        "| Recip_lattice_vector                :"
    call localorb_info( rrs_pbc_info )
    do i = 1, 3, 1
       write(rrs_pbc_info,'(4X,A,3(2X,F16.8),2X)') &
            "|", &
            (rrs_pbc_recip_lattice_vector(j,i) / bohr, j=1,3,1)
       call localorb_info( rrs_pbc_info )
    enddo
    ! print out inverse of the lattice vector
    write(rrs_pbc_info,'(4X,A)') &
        "| Inv_lattice_vector                  :"
    call localorb_info( rrs_pbc_info )
    do i = 1, 3, 1
        write(rrs_pbc_info,'(4X,A,3(2X, F16.8),2X)') &
            '|', &
            (rrs_pbc_inv_lattice_vector(j,i) /bohr, j=1,3,1)
        call localorb_info(rrs_pbc_info)
    enddo

    write(rrs_pbc_info,'(4X,A,E15.6,X,A)') & 
        '| Unit cell volume                     :', &
        rrs_pbc_cell_volume * bohr**3, '(A^3)'
    call localorb_info(rrs_pbc_info,use_unit,'(A)')

    end subroutine
