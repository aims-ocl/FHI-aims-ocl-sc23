!****s* FHI-aims/trans_first_order_sparse_to_matrix
!  NAME
!   trans_first_order_sparse_to_matrix
!  SYNOPSIS

subroutine trans_first_order_sparse_to_matrix(k_coord,k_center_in_sc_DFPT, first_order_M_sparse, & 
                                              first_order_M_full )


  !  PURPOSE
  !   translate first_order_M_sparse(3,n_centers_in_sc_DFPT, n_hamiltonian_matrix_size) to 
  !   full matrix form first_order_M_full(n_basis_sc_DFPT,n_basis_sc_DFPT) 
  !                    @ k_coord and k_center
  !   shanghui 2015.07
  ! USES

  use pbc_lists
  use dimensions
  implicit none

  !  ARGUMENTS

  integer, intent(in) :: k_coord
  integer, intent(in) :: k_center_in_sc_DFPT
  real*8,intent(in) :: first_order_M_sparse(3,n_centers_in_sc_DFPT,n_hamiltonian_matrix_size)
  real*8,intent(out) :: first_order_M_full(n_basis_sc_DFPT,n_basis_sc_DFPT)  

  !  INPUTS
  !    o first_order_M_sparse -- sparse matrix using PM_index
  !  OUTPUT
  !    o first_order_M_full -- full matrix 
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
  integer :: k_cell_in_sc_DFPT, k_atom ,k_cell_trans, k_center_trans

 
  k_cell_in_sc_DFPT = center_in_sc_DFPT_to_cell_in_sc_DFPT(k_center_in_sc_DFPT)
  k_atom            = center_in_sc_DFPT_to_atom(k_center_in_sc_DFPT)

  first_order_M_full = 0.0d0
 
  do iuo = 1, n_basis
  do i_cell_in_hamiltonian = 1, n_cells_in_hamiltonian-1

     i_cell_in_sc_DFPT = cell_in_hamiltonian_to_cell_in_sc_DFPT(i_cell_in_hamiltonian)

     if (index_hamiltonian(1,i_cell_in_hamiltonian, iuo) > 0) then
     do i_index = index_hamiltonian(1, i_cell_in_hamiltonian, iuo), &
            &     index_hamiltonian(2, i_cell_in_hamiltonian, iuo)

      !------------------(1) j_cell = 1 ----------------------------   
      !                                                            ! 
      !       d M_full (u R1, v)          d M_sparse (u R1, v)     !
      !    ----------------------- =    -----------------------    !
      !        d R_K                      d R_K                    !
      !                                                            !


            juo =  column_index_hamiltonian(i_index)
            io  =  cell_and_basis_to_basis_sc_DFPT(i_cell_in_sc_DFPT,iuo)

            first_order_M_full(io,juo) =  &  
            first_order_M_sparse(k_coord,k_center_in_sc_DFPT,i_index)

            first_order_M_full(juo,io) =  &  
            first_order_M_sparse(k_coord,k_center_in_sc_DFPT,i_index)


      !------------------(2) j_cell > 1 ----------------------------   
      !                                                            ! 
      !    d M_full (u R1+R2, v+R2)      d M_sparse (u R1, v)      !
      !    ----------------------- =    -----------------------    !
      !        d R_K [k_center]       d R_K-R2 [k_center_trans]    !
      !                                                            !
            
        do j_cell_trans = 2, n_cells_in_sc_DFPT
            ! jo_trans = juo +j_cell 
           jo_trans = cell_and_basis_to_basis_sc_DFPT(j_cell_trans,juo)  

            ! io_trans = iuo + i_cell + j_cell
           i_cell_trans = cell_add_sc_DFPT(i_cell_in_sc_DFPT,j_cell_trans)
           io_trans = cell_and_basis_to_basis_sc_DFPT(i_cell_trans,iuo) 

            ! k_center_trans = k_center - j_cell
           k_cell_trans   = cell_diff_sc_DFPT(k_cell_in_sc_DFPT,j_cell_trans)
           k_center_trans = cell_and_atom_to_center_sc_DFPT(k_cell_trans,k_atom) 

 
           ! first_order_M_full(k_coord,k_center_in_sc_DFPT,io_trans,jo_trans)
           first_order_M_full(io_trans,jo_trans) = &  
           first_order_M_sparse(k_coord,k_center_trans,i_index)

           
           first_order_M_full(jo_trans,io_trans) = &  
           first_order_M_sparse(k_coord,k_center_trans,i_index)

        enddo

      end do
      end if
   end do
   end do


  end subroutine trans_first_order_sparse_to_matrix
