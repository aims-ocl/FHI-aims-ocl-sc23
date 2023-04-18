!****s* FHI-aims/prune_first_order_density_matrix_sparse
!  NAME
!    prune_first_order_density_matrix_sparse
!  SYNOPSIS

subroutine prune_first_order_density_matrix_sparse(first_order_density_matrix_sparse, & 
                                                   first_order_density_matrix_con, &
                                                   n_compute, i_basis_index)

!  PURPOSE
!    The subroutine saves the first-order-density matrix components belongs to non-zero basis functions
!    to first_order_density_matrix_con.
!
!  USES

  use dimensions
  use pbc_lists
  implicit none

!  ARGUMENTS

  real*8 :: first_order_density_matrix_sparse(3,n_centers_in_sc_DFPT,n_hamiltonian_matrix_size)
  real*8 :: first_order_density_matrix_con(3,n_centers_in_sc_DFPT,n_compute, n_compute)
  integer:: n_compute
  integer:: i_basis_index(n_compute)

!  INPUTS
!    o first_order_density_matrix_sparse -- total density matrix (packed matrix format)
!    o n_compute -- number of non-zero basis function in current grid batch 
!    o i_basis -- list of the non-zero basis functions in current grid batch
!
!  OUTPUT
!   o first_order_density_matrix_con -- values of density matrix balong to non-zero basis functions.
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
! SOURCE

  integer :: i_compute, j_compute
  integer :: i_basis, j_basis
  integer :: i_basis_uc, j_basis_uc 
  integer :: i_cell_in_hamiltonian, j_cell_in_hamiltonian 
  integer :: i_cell_in_sc_DFPT , j_cell_in_sc_DFPT
  integer :: i_coord, i_place
  integer :: k_center_in_sc_DFPT , k_cell_in_sc_DFPT, k_atom, k_cell_trans, k_center_trans 

  
  first_order_density_matrix_con = 0.d0

   do i_compute = 1,n_compute
      i_basis    = i_basis_index(i_compute)
      i_basis_uc = Cbasis_to_basis(i_basis)
      i_cell_in_hamiltonian     = center_to_cell(Cbasis_to_center(i_basis))
      i_cell_in_sc_DFPT         = center_in_sc_DFPT_to_cell_in_sc_DFPT( &
                                  center_to_center_in_sc_DFPT(Cbasis_to_center(i_basis)) )


   do j_compute = 1,n_compute
      j_basis    = i_basis_index(j_compute)
      j_basis_uc = Cbasis_to_basis(j_basis)

      if(j_basis_uc <= i_basis_uc) then
         j_cell_in_hamiltonian = center_to_cell(Cbasis_to_center(j_basis))
         j_cell_in_sc_DFPT     = center_in_sc_DFPT_to_cell_in_sc_DFPT( &
                                 center_to_center_in_sc_DFPT(Cbasis_to_center(j_basis)) )

      do i_place = &
      index_hamiltonian(1,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc),  &
      index_hamiltonian(2,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc)

         do i_coord =1 ,3
         do k_center_in_sc_DFPT =1 ,n_centers_in_sc_DFPT
            k_cell_in_sc_DFPT = center_in_sc_DFPT_to_cell_in_sc_DFPT(k_center_in_sc_DFPT)
            k_atom = center_in_sc_DFPT_to_atom(k_center_in_sc_DFPT)

            k_cell_trans = cell_diff_sc_DFPT(k_cell_in_sc_DFPT,j_cell_in_sc_DFPT)
            k_center_trans = cell_and_atom_to_center_sc_DFPT(k_cell_trans, k_atom)

         if( column_index_hamiltonian( i_place) == j_basis_uc)then

              first_order_density_matrix_con(i_coord,k_center_in_sc_DFPT,i_compute,j_compute) = & 
              first_order_density_matrix_sparse(i_coord,k_center_trans,i_place) 

              first_order_density_matrix_con(i_coord,k_center_in_sc_DFPT,j_compute,i_compute) = & 
              first_order_density_matrix_sparse(i_coord,k_center_trans,i_place) 

         endif

         enddo ! i_center
         enddo ! i_coord 

      enddo ! i_place
      endif ! j_basis_uc <= i_basis_uc

   enddo ! j_compute 
   enddo ! i_compute 


end subroutine prune_first_order_density_matrix_sparse
!******
