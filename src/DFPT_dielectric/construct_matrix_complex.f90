!****s* FHI-aims/construct_matrix_complex
!  NAME
!   construct_matrix_complex
!  SYNOPSIS

subroutine construct_matrix_complex(matrix_sparse, matrix_complex,i_k_point)

  !  PURPOSE
  !   core code to contruct matrix_complex on all k_point
  !
  !   called by 
  !    construct_first_order_H_p1
  !
  !                                     shanghui 2015.11.16 


  ! USES

  use pbc_lists
  use dimensions
  implicit none

  !  ARGUMENTS

  real*8,intent(in) :: matrix_sparse(n_hamiltonian_matrix_size)
  complex*16,intent(inout) :: matrix_complex(n_basis,n_basis)
  integer, intent(in) :: i_k_point

  !  INPUTS
  !    o matrix_sparse -- matrix in PM_index
  !  OUTPUT
  !    o matrix_complex -- matrix in kpoint
  !  AUTHOR
  !    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !     Computer Physics Communications (2008), submitted.
  !  COPYRIGHT
  !   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
  !   e.V. Please note that any use of the "FHI-aims-Software" is subject to
  !   the terms and conditions of the respective license agreement."
  !  HISTORY
  !    Release version, FHI-aims (2008).
  !
  ! SOURCE

  integer :: i_cell, i_place, i_basis_1, i_basis_2

  ! initialize
  matrix_complex= (0.0d0,0.0d0)


     do i_cell = 1,n_cells_in_hamiltonian-1
     do i_basis_1 = 1, n_basis
    
        if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
        do i_place = index_hamiltonian(1,i_cell, i_basis_1), & 
                     index_hamiltonian(2,i_cell, i_basis_1)
    
           i_basis_2 =  column_index_hamiltonian(i_place)
    
         !----------(1) i_basis_1 > i_basis_2 case --------------
           matrix_complex(i_basis_1, i_basis_2) =  &
           matrix_complex(i_basis_1, i_basis_2) +  &
           dconjg(k_phase(i_cell,i_k_point))* &  !This is what did in construct_hamiltonian_and_ovl.f90 
           !k_phase(i_cell,i_k_point)* & 
           matrix_sparse( i_place) 
    
         !----------(2) i_basis_1 < i_basis_2--------------------------------------
          if(i_basis_1.ne.i_basis_2) then
           matrix_complex( i_basis_2, i_basis_1) =  &
           matrix_complex( i_basis_2, i_basis_1) +  &
           k_phase(i_cell,i_k_point)* & !This is what did in construct_hamiltonian_and_ovl.f90 for i < j
           !dconjg(k_phase(i_cell,i_k_point))* & 
           matrix_sparse( i_place) 
          endif
    
        enddo 
        endif
    
     enddo ! i_basis_1
     enddo ! i_cell
 

     
end subroutine construct_matrix_complex
!******
