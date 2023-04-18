!****s* FHI-aims/prune_density_matrix_sparse_polar_reduce_memeory
!  NAME
!    prune_density_matrix_sparse_polar_reduce_memory
!  SYNOPSIS

subroutine prune_density_matrix_sparse_polar_reduce_memory(density_matrix_sparse, & 
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

  real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
  real*8 :: density_matrix_con(n_compute, n_compute)
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


  integer :: i_basis,j_basis,i_compute,j_compute, i_start, i_end, i_place,  & 
             i_index_real

  !-------initial----------------------- 
  density_matrix_con = 0.0d0

  !-----(1) begin the slow version-----------------
  !do i_compute = 1,n_compute
  !   i_basis = i_basis_index(i_compute)
 
  !   i_start = index_hamiltonian(1,1, i_basis)
  !   i_end   = index_hamiltonian(2,1, i_basis)

  ! ! shanghui note: sparse matrix only have the i_basis > j_basis case  
  ! !                see prune_density_matrix_sparse.f90 
  !do j_compute = 1, i_compute
  !   j_basis = i_basis_index(j_compute)
 
  !   do i_place = i_start, i_end, 1

  !      if( column_index_hamiltonian(i_place) == j_basis)then

  !           density_matrix_con(i_compute,j_compute) = & 
  !           density_matrix_sparse(i_place) 

  !           density_matrix_con(j_compute,i_compute) = & 
  !           density_matrix_sparse(i_place) 

  !      endif

  !   enddo ! i_place

  !enddo ! j_compute 
  !enddo ! i_compute 
  !-----(1) end the slow version-----------------


  !-----(2) begin the faster version-----------------
   do i_compute = 1, n_compute, 1
      i_basis = i_basis_index(i_compute)
      i_start =  index_hamiltonian(1,1, i_basis)
      i_end   =  index_hamiltonian(2,1, i_basis)

      if(i_end<i_start) cycle ! otherways an undefined i_index_real is used!

      do j_compute = 1, i_compute, 1
         j_basis = i_basis_index(j_compute)

         place: do i_place = i_start, i_end, 1
            
            if( column_index_hamiltonian( i_place) == j_basis)then

               density_matrix_con(i_compute, j_compute) = density_matrix_sparse(i_place)
               density_matrix_con(j_compute, i_compute) = density_matrix_sparse(i_place)
               i_index_real = i_place
               exit place 
          
            else if(column_index_hamiltonian( i_place) > j_basis)then
               i_index_real = i_place
               exit place
               
            end if
         end do place
         i_start = i_index_real
      end do
   end do
  !-----(2) end the faster version-----------------


end subroutine prune_density_matrix_sparse_polar_reduce_memory
!******
