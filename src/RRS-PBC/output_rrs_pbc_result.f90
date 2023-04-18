!****s* FHI-aims/RRS-PBC/output_rrs_pbc_result
!  NAME
!   output_rrs_pbc_result
!  SYNOPSIS

    subroutine output_rrs_pbc_result(k_index,k_vector,iop)
!  USES
    use dimensions
    use physics
    use runtime_choices
    use localorb_io
!  PURPOSE
!   The subroutine print out the lattice information of RRS-PBC method.
!   This subroutine is called from parse_rrs_pbc
!  INPUTS
!     k_index  :: the index of the k point
!     k_vector :: the k vector in cartesian coordinate
!     iop      :: 1 to store KS_eigenvector
!     iop      :: 2 to print out KS_eigenvalue 
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

    integer             :: k_index, iop
    real*8,dimension(3) :: k_vector
    character*130       :: rrs_pbc_info
    integer             :: i,j

    ! print out lattice vector
    if (iop .eq. 1) then
        write(rrs_pbc_info,'(A,I5.5,A)') 'RRS_PBC_KS_eigenvector_',k_index,'.dat'
        open(unit=1982,file=rrs_pbc_info,FORM='UNFORMATTED', STATUS='REPLACE')
        write(1982) k_vector
        write(1982) rrs_pbc_KS_eigenvector
        close(1982)
    else if (iop .eq. 2) then
        write(rrs_pbc_info,'(A)') 'RRS_PBC_KS_eigenvalue.dat'
        open(unit=1982,file=rrs_pbc_info,FORM='UNFORMATTED', STATUS='REPLACE')
        write(1982) rrs_pbc_KS_eigenvalue
        close(1982)
    endif

    end subroutine
