!****s* FHI-aims/unpack_matrices
!  NAME
!    unpack_matrices
!  SYNOPSIS

subroutine unpack_matrices(density_matrix_sparse, density_matrix_con, &
n_compute)

  !  PURPOSE
  !   The subroutine unpacks packed matrices
  !
  !  USES

  use dimensions
  use runtime_choices
  use pbc_lists
  implicit none

  !  ARGUMENTS

  real*8 :: density_matrix_sparse(n_hamiltonian_matrix_size)
  real*8 :: density_matrix_con(n_compute, n_compute)
  integer:: n_compute

  !  INPUTS
  !    o density_matrix_sparse -- any matrix (packed matrix format)
  !    o n_compute == n_centers_basis_T -- number of basis functions
  !
  !  OUTPUT
  !   o density_matrix_con -- values of density matrix belong to non-zero basis functions.
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



  integer :: i_compute_1, i_compute_2
  integer :: i_c1, i_c2
  integer :: i_cell_index, i_index_real
  integer :: i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: i_cell_old, i_cell_1
  integer :: offset(n_cells)
  integer :: offset_end(n_cells) 


  density_matrix_con = 0.d0

  if (packed_matrix_format /= PM_index) return

  if(n_periodic == 0)then

     do i_c2 = 1, n_compute, 1

        i_start =  index_hamiltonian(1,1, i_c2)
        i_end   =  index_hamiltonian(2,1, i_c2)


        do i_c1 = 1, i_c2, 1


           i_basis_2 = i_c1

           place: do i_place = i_start, i_end, 1

              if( column_index_hamiltonian( i_place) == i_basis_2)then

                 density_matrix_con(i_c1, i_c2) = density_matrix_sparse(i_place)
                 density_matrix_con(i_c2, i_c1) = density_matrix_sparse(i_place)
                 i_index_real = i_place
                 exit place 

              else if(column_index_hamiltonian( i_place) > i_basis_2)then
                 i_index_real = i_place
                 exit place

              end if
           end do place
           i_start = i_index_real
        end do
     end do


  else ! The periodic case


     ! Unfortunately the periodic systems can not use the searching routine used now in the clusters.
     ! This is because the peridic systems have extra packing for supercell information.

     do i_compute_1 = 1, n_compute, 1

        i_basis_1 = Cbasis_to_basis(i_compute_1)
        i_cell_old = center_to_cell(Cbasis_to_center(i_compute_1))

        do i_cell_1 = 1, n_cells 


           i_cell = position_in_hamiltonian( i_cell_old, i_cell_1) 

           offset(i_cell_1)     = index_hamiltonian(1,i_cell, i_basis_1)
           offset_end(i_cell_1) = index_hamiltonian(2,i_cell, i_basis_1)

        end do


        do i_compute_2 = 1, n_compute

           i_basis_2 = Cbasis_to_basis(i_compute_2)
           i_cell    =  center_to_cell(Cbasis_to_center(i_compute_2))


           if(i_basis_2 <= i_basis_1)then

              place_2: do i_place = offset(i_cell), offset_end(i_cell),1 

                 if( column_index_hamiltonian( i_place) == i_basis_2)then

                    density_matrix_con(i_compute_1, i_compute_2) = density_matrix_sparse(i_place)
                    density_matrix_con(i_compute_2, i_compute_1) = density_matrix_sparse(i_place)
                    exit place_2



                 else if(column_index_hamiltonian( i_place) > i_basis_2)then

                    exit place_2  

                 end if
              end do place_2
              offset(i_cell) = i_place
           end if
        end do
     end do


     ! The old routine just in case
     !
     !        do i_c2 = 1, n_compute, 1
     !           do i_c1 = 1, i_c2, 1
     !
     !
     !              if(Cbasis_to_basis(i_basis(i_c1)) <= Cbasis_to_basis(i_basis(i_c2)))then
     !
     !                 i_compute_2 = i_basis(i_c2)
     !                 i_compute_1 = i_basis(i_c1)
     !              else
     !                 i_compute_2 = i_basis(i_c1)
     !                 i_compute_1 = i_basis(i_c2)
     !              end if
     !
     !
     !           
     !              i_index_real =  find_position_in_hamiltonian( i_compute_2, i_compute_1) 
     !              density_matrix_con(i_c1, i_c2) = density_matrix_sparse(i_index_real)
     !
     !
     !           end do
     !        end do

  end if ! n_periodic == 0 - else


end subroutine unpack_matrices
!******
