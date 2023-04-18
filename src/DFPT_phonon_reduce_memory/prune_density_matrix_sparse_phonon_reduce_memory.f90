!****s* FHI-aims/prune_density_matrix_sparse_phonon_reduce_memory
!  NAME
!    prune_density_matrix_sparse_phonon_reduce_memory
!  SYNOPSIS

subroutine prune_density_matrix_sparse_phonon_reduce_memory(density_matrix_sparse, & 
                                       density_matrix_con, &
                                       n_compute, i_basis_index)

!  PURPOSE
!    The subroutine change the density matrix components belongs to non-zero basis functions
!    to density_matrix_con.
!
!  USES

  use dimensions
  use pbc_lists
  implicit none

!  ARGUMENTS

  complex*16 :: density_matrix_sparse(n_hamiltonian_matrix_size)
  complex*16 :: density_matrix_con(n_compute, n_compute)
  integer:: n_compute
  integer:: i_basis_index(n_compute)

!  INPUTS
!    o density_matrix_sparse -- total density matrix (packed matrix format)
!    o n_compute -- number of non-zero basis function in current grid batch 
!    o i_basis -- list of the non-zero basis functions in current grid batch
!
!  OUTPUT
!   o density_matrix_con -- values of density matrix balong to non-zero basis functions.
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


  integer :: i_basis,j_basis,i_compute,j_compute,i_basis_uc,j_basis_uc,i_cell,j_cell
  integer :: i_coord, i_center,i_place
  integer :: k_cell,k_atom, k_cell_new,i_center_new

  !-------initial----------------------- 
  density_matrix_con = (0.d0,0.d0)

   do i_compute = 1,n_compute
      i_basis    = i_basis_index(i_compute)
      i_basis_uc = Cbasis_to_basis(i_basis)
      i_cell     = center_to_cell(Cbasis_to_center(i_basis))

   do j_compute = 1,n_compute
      j_basis    = i_basis_index(j_compute)
      j_basis_uc = Cbasis_to_basis(j_basis)

      if(j_basis_uc <= i_basis_uc) then
         j_cell = center_to_cell(Cbasis_to_center(j_basis))
 
      do i_place = &
         index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
         index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

         if( column_index_hamiltonian( i_place) == j_basis_uc)then

              density_matrix_con(i_compute,j_compute) = & 
              density_matrix_sparse(i_place) 

              density_matrix_con(j_compute,i_compute) = & 
              density_matrix_sparse(i_place) 

         endif

      enddo ! i_place
      endif ! j_basis_uc <= i_basis_uc

   enddo ! j_compute 
   enddo ! i_compute 



!-----------------shanghui begin to test DM1 and EDM1--------------
!  do i_basis = 1,n_Cbasis
!     i_basis_uc = Cbasis_to_basis(i_basis)
!     i_cell     = center_to_cell(Cbasis_to_center(i_basis))
!
!  do j_basis = 1,n_Cbasis
!     j_basis_uc = Cbasis_to_basis(j_basis)
!
!     if(j_basis_uc <= i_basis_uc) then
!        j_cell = center_to_cell(Cbasis_to_center(j_basis))
!
!     do i_place = &
!        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
!        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)
!
!        if( column_index_hamiltonian( i_place) == j_basis_uc)then
!
!
!          if(myid.eq.0)  write(use_unit,*) i_basis, j_basis, density_matrix_sparse(i_place)
!
!
!        endif
!
!     enddo ! i_place
!     endif ! j_basis_uc <= i_basis_uc
!
!  enddo ! j_compute 
!  enddo ! i_compute 
!

end subroutine prune_density_matrix_sparse_phonon_reduce_memory
!******
