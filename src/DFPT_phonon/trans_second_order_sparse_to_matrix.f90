!****s* FHI-aims/trans_second_order_sparse_to_matrix
!  NAME
!   trans_second_order_sparse_to_matrix
!  SYNOPSIS

subroutine trans_second_order_sparse_to_matrix(m_coord,m_center,n_coord, n_center, &  
                                              second_order_M_sparse, & 
                                              second_order_M_full )


  !  PURPOSE
  !   translate 
  ! second_order_M_sparse(3,n_centers_in_hamiltonian, 3,n_centers_in_hamiltonian, n_hamiltonian_matrix_size) to 
  !   full matrix form 
  ! second_order_M_full(n_basis_supercell,n_basis_supercell) 
  !                    @ m_coord m_center and n_coord n_center
  !   shanghui 2015.07
  ! USES

  use pbc_lists
  use dimensions
  implicit none

  !  ARGUMENTS

  integer, intent(in) :: m_coord, n_coord
  integer, intent(in) :: m_center, n_center
  real*8,intent(in) :: &  
  second_order_M_sparse(3,n_centers_in_hamiltonian,3,n_centers_in_hamiltonian, n_hamiltonian_matrix_size)
  real*8,intent(out) :: second_order_M_full(n_basis_supercell,n_basis_supercell)  

  !  INPUTS
  !    o second_order_M_sparse -- sparse matrix using PM_index
  !  OUTPUT
  !    o second_order_M_full   -- full matrix 
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
  integer :: i_cell,i_index,i_cell_trans,j_cell_trans
  integer :: m_cell,m_cell_trans,m_atom,m_center_trans
  integer :: n_cell,n_cell_trans,n_atom,n_center_trans
 
  m_cell = center_to_cell(centers_in_hamiltonian(m_center))
  m_atom = center_to_atom(centers_in_hamiltonian(m_center))

  n_cell = center_to_cell(centers_in_hamiltonian(n_center))
  n_atom = center_to_atom(centers_in_hamiltonian(n_center))

  second_order_M_full = 0.0d0
 
  do iuo = 1, n_basis
  do i_cell = 1, n_cells_in_hamiltonian-1
     if (index_hamiltonian(1,i_cell, iuo) > 0) then
     do i_index = index_hamiltonian(1, i_cell, iuo), &
            &     index_hamiltonian(2, i_cell, iuo)

            juo =  column_index_hamiltonian(i_index)
            io  =  cell_and_basis_to_basis_supercell(i_cell,iuo)

            second_order_M_full(io,juo) =  &  
            second_order_M_sparse(m_coord,m_center,n_coord,n_center,i_index)

            second_order_M_full(juo,io) =  &  
            second_order_M_sparse(m_coord,m_center,n_coord, n_center, i_index)
            
        do j_cell_trans = 2,n_cells_in_hamiltonian-1
            ! jo_trans = juo +j_cell 
           jo_trans = cell_and_basis_to_basis_supercell(j_cell_trans,juo)  

            ! io_trans = iuo + i_cell + j_cell
           i_cell_trans = position_in_hamiltonian_PBC_add(i_cell,j_cell_trans)
           io_trans = cell_and_basis_to_basis_supercell(i_cell_trans,iuo) 

            ! m_center_trans = m_center - j_cell
           m_cell_trans = position_in_hamiltonian_PBC(m_cell,j_cell_trans)
           m_center_trans = inv_centers_in_hamiltonian( & 
           cell_and_atom_to_center(m_cell_trans, m_atom))
 
            ! n_center_trans = n_center - j_cell
           n_cell_trans = position_in_hamiltonian_PBC(n_cell,j_cell_trans)
           n_center_trans = inv_centers_in_hamiltonian( & 
           cell_and_atom_to_center(n_cell_trans, n_atom))

           ! second_order_M_full(m_coord,m_center,n_coord,n_center,io_trans,jo_trans)
           second_order_M_full(io_trans,jo_trans) = &  
           second_order_M_sparse(m_coord,m_center_trans,n_coord,n_center_trans,i_index)

           
           second_order_M_full(jo_trans,io_trans) = &  
           second_order_M_sparse(m_coord,m_center_trans,n_coord,n_center_trans,i_index)

        enddo

      end do
      end if
   end do
   end do


  end subroutine trans_second_order_sparse_to_matrix
