!****s* FHI-aims/ construct_first_order_matrix_k_complex
!  NAME
!   construct_first_order_matrix_k_complex
!  SYNOPSIS

subroutine construct_first_order_matrix_k_complex(matrix_sparse, matrix_complex,i_k_point, i_coord)

  !  PURPOSE
  !   core code to contruct matrix_complex for dSuv(k)/dk_i
  !
  !
  !                                     shanghui 2017.01.10 


  ! USES

  use pbc_lists
  use dimensions
  use geometry, only: lattice_vector
  implicit none

  !  ARGUMENTS

  real*8,intent(in) :: matrix_sparse(n_hamiltonian_matrix_size)
  complex*16,intent(inout) :: matrix_complex(n_basis,n_basis)
  integer, intent(in) :: i_k_point
  integer, intent(in) :: i_coord

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
  real*8  :: R_i_cell


  !write(use_unit,*) 'lattice_vector(i_coord):',i_coord, lattice_vector(i_coord,1:3)


  ! initialize
  matrix_complex= (0.0d0,0.0d0)


     do i_cell = 1,n_cells_in_hamiltonian-1
  ! write(use_unit,*) 'cell_index:',cell_index(i_cell,1:3)
 
        R_i_cell = cell_index(i_cell, 1)*lattice_vector(i_coord,1) + &
                   cell_index(i_cell, 2)*lattice_vector(i_coord,2) + &
                   cell_index(i_cell, 3)*lattice_vector(i_coord,3)

     do i_basis_1 = 1, n_basis
    
        if( index_hamiltonian(1,i_cell, i_basis_1) > 0 )then
        do i_place = index_hamiltonian(1,i_cell, i_basis_1), & 
                     index_hamiltonian(2,i_cell, i_basis_1)
    
           i_basis_2 =  column_index_hamiltonian(i_place)
    
         !----------(1) i_basis_1 > i_basis_2 case --------------
           matrix_complex(i_basis_1, i_basis_2) =  &
           matrix_complex(i_basis_1, i_basis_2) +  &
           dconjg(k_phase(i_cell,i_k_point))* & 
           matrix_sparse( i_place) * (-R_i_cell)  
    
         !----------(2) i_basis_1 < i_basis_2--------------------------------------
          if(i_basis_1.ne.i_basis_2) then
           matrix_complex( i_basis_2, i_basis_1) =  &
           matrix_complex( i_basis_2, i_basis_1) +  &
           k_phase(i_cell,i_k_point)* & 
           matrix_sparse( i_place) * (R_i_cell) 
          endif
    
        enddo 
        endif
    
     enddo ! i_basis_1
     enddo ! i_cell
 

     
end subroutine construct_first_order_matrix_k_complex
!******
