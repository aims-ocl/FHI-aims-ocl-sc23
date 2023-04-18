!****s* FHI-aims/add_matrix_to_sparse
!  NAME
!   add_matrix_to_sparse
!  SYNOPSIS

subroutine add_matrix_to_sparse( k_atom, M_full, M_sparse)


  !  PURPOSE
  !   add full matrix form M_full(n_basis_sc_DFPT,n_basis_sc_DFPT)
  !   to M_sparse(n_hamiltonian_matrix_size) 
  !   shanghui 2015.07
  ! USES

  use pbc_lists
  use dimensions
  implicit none

  !  ARGUMENTS
  
  integer, intent(in) :: k_atom
  real*8,intent(in) :: M_full(n_basis_sc_DFPT,n_basis_sc_DFPT)  
  real*8,intent(out) :: M_sparse(n_centers_in_sc_DFPT, n_hamiltonian_matrix_size)

  !  INPUTS
  !    o M_full -- full matrix 
  !  OUTPUT
  !    o M_sparse -- sparse matrix using PM_index
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
  integer :: i_cell,i_cell_trans, i_index,j_cell,j_cell_trans, i_cell_in_sc_DFPT
  integer :: k_cell, k_center
 

 
  do iuo = 1, n_basis
  do i_cell = 1, n_cells_in_hamiltonian-1

     i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell)    

     if (index_hamiltonian(1,i_cell, iuo) > 0) then
     do i_index = index_hamiltonian(1, i_cell, iuo), &
            &     index_hamiltonian(2, i_cell, iuo)


      !------------------(1) j_cell = 1 ----------------------------   
      !                                                            ! 
      !      d M_sparse (u R1, v)           d M_full (u R1, v)     !
      !   -----------------------   +=   -----------------------   !
      !      d R_K                           d R_K                 !
      !                                                            !
 
          juo =  column_index_hamiltonian(i_index)
          io  =  cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

          M_sparse(k_atom,i_index) = M_sparse(k_atom,i_index) +  M_full(io,juo)


      !------------------(2) j_cell > 1 ----------------------------   
      !                                                            ! 
      !      d M_sparse (u R1, v)         d M_full (u R1+R2, v R2) !
      !   -----------------------   =   -----------------------    !
      !      d R_K-R2 [k_center]             d R_K [k_atom]        !
      !                                                            !
         do j_cell_trans = 2,n_cells_in_sc_DFPT                             ! j_cell_trans = R2 

            jo_trans     = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)
            i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans) ! i_cell_trans = R1 + R2
            io_trans     = cell_and_basis_to_basis_sc_DFPT( i_cell_trans,iuo) 
          
            k_cell       = cell_diff_sc_DFPT(1,j_cell_trans)                ! k_cell = 0 - R2
            k_center     = cell_and_atom_to_center_sc_DFPT( k_cell, k_atom)   
            M_sparse(k_center,i_index) = M_sparse(k_center,i_index) +  M_full(io_trans,jo_trans)

        enddo 

      end do
      end if



   end do
   end do

  end subroutine add_matrix_to_sparse
