!****s* FHI-aims/trans_sparse_to_matrix
!  NAME
!   trans_sparse_to_matrix
!  SYNOPSIS

subroutine trans_sparse_to_matrix( M_sparse, M_full )


  !  PURPOSE
  !   translate        M_sparse(n_hamiltonian_matrix_size) to 
  !   full matrix form M_full(n_basis_sc_DFPT,n_basis_sc_DFPT)
  !   shanghui 2015.07
  ! USES

  use pbc_lists
  use dimensions
  implicit none

  !  ARGUMENTS

  real*8,intent(in) :: M_sparse(n_hamiltonian_matrix_size)
  real*8,intent(out) :: M_full(n_basis_sc_DFPT,n_basis_sc_DFPT)  

  !  INPUTS
  !    o M_sparse -- sparse matrix using PM_index
  !  OUTPUT
  !    o M_full -- full matrix 
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


  integer :: iuo,io,io_trans,juo,jo,jo_trans 
  integer :: i_cell_in_hamiltonian, i_cell_in_sc_DFPT, i_index
  integer :: i_cell_trans, j_cell_trans
 
  M_full = 0.0d0
 
  do iuo = 1, n_basis
  do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1

     i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian) 

     if (index_hamiltonian(1,i_cell_in_hamiltonian, iuo) > 0) then
     do i_index = index_hamiltonian(1, i_cell_in_hamiltonian, iuo), &
            &     index_hamiltonian(2, i_cell_in_hamiltonian, iuo)

            juo =  column_index_hamiltonian(i_index)
            io  =  cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

            M_full(io,juo) =  M_sparse(i_index)
            M_full(juo,io) =  M_sparse(i_index)
            
         do j_cell_trans = 2,n_cells_in_sc_DFPT
            ! jo_trans = juo +j_cell 
            jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)  

            i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
            ! io_trans = iuo + i_cell_in_sc_DFPT + j_cell_in_sc_DFPT
            io_trans = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo)  
           
            M_full(io_trans,jo_trans) = M_sparse(i_index)
            M_full(jo_trans,io_trans) = M_sparse(i_index)
         enddo

      end do
      end if
   end do
   end do

  end subroutine trans_sparse_to_matrix
