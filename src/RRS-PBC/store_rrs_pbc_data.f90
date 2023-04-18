!****s* FHI-aims/RRS-PBC/store_rrs_pbc_date
!  NAME
!   store_rrs_pbc_data
!  SYNOPSIS

    subroutine store_rrs_pbc_data()
!  USES
    use dimensions
    use physics
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

    integer :: NAE, NBE


    Open(unit=1982,file='Parameter.dat',FORM='UNFORMATTED',STATUS='REPLACE')
    Write(1982) n_basis
    Write(1982) n_hamiltonian_matrix_size
    Write(1982) int(rrs_pbc_n_electron(1))
    Write(1982) int(rrs_pbc_n_electron(2))
    Close(1982)

!  1 : the Fock matrix for alpha electrons
    open(unit=1982,file='FAdata.dat',FORM='UNFORMATTED',STATUS='REPLACE')
    write(1982) hamiltonian(:,1)
! 2: the Fock matrix for beta eletrons
    If (NAE.ne.NBE) then
      open(unit=1982,file='FBdata.dat',FORM='UNFORMATTED',STATUS='REPLACE')
      write(1982) hamiltonian(:,2)
      close(1982)
    EndIf
! 3: the overlap matrix
    open(unit=1982,file='OLdata.dat',FORM='UNFORMATTED',STATUS='REPLACE')
    write(1982) overlap_matrix
    close(1982)
    end subroutine

