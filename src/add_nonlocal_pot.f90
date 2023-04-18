!****s* FHI-aims/add_nonlocal_pot
!  NAME
!    add_nonlocal_pot
!  SYNOPSIS
!
!

  subroutine add_nonlocal_pot(hamiltonian)

    !  PURPOSE
    !
    !    Adds the nonlocal potential to the hamiltonian. Hamiltonian matrix must not
    !    be of reduced size. Reason: even though two basisfctns may not overlap by themselves, 
    !    they are still may overlap with the same pseudo fctns and therefore contribute to the 
    !    nonlocal potentials.
    !          H(i,j) --> H(i,j) + \sum_{lm} <i|X_lm> V_nonlocal <X_lm|j>
    !  USES
    use dimensions
    use pseudodata

    implicit none

    !  ARGUMENTS
       real*8 :: hamiltonian(n_basis*(n_basis+1)/2,n_spin)

    !  AUTHOR
    !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
    !  HISTORY
    !    Release version, FHI-aims (2011).
    !  SOURCE


!     local variables


    ! counter
    
    integer :: i_basis_1, i_basis_2, i_spin, i_index

!    real*8, dimension((n_basis+1)*n_basis/2) :: nonlocal_matrix_2

!     nonlocal_matrix_2 = 0


!        do i_basis_2 = 1, n_basis
!          do i_basis_1 = 1, i_basis_2
!            i_index = i_basis_1 + (i_basis_2-1)*i_basis_2/2
!                  nonlocal_matrix_2(i_index) = nonlocal_matrix( i_basis_1, i_basis_2)
!           enddo 
!        enddo



      do i_spin = 1, n_spin

!        do i_index= 1, (n_basis+1)*n_basis/2

!            hamiltonian(i_index, i_spin) = &
!              hamiltonian(i_index, i_spin) &
!              + nonlocal_matrix_2(i_index)

        do i_basis_2 = 1, n_basis
          do i_basis_1 = 1, i_basis_2
            i_index = i_basis_1 + (i_basis_2-1)*i_basis_2/2
                 hamiltonian(i_index, i_spin) = &
                 hamiltonian(i_index, i_spin) &
                 + nonlocal_matrix( i_basis_1, i_basis_2)


          enddo 


        enddo
    
      enddo



  end subroutine add_nonlocal_pot

