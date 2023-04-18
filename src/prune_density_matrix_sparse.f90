!****s* FHI-aims/prune_density_matrix_sparse
!  NAME
!    prune_density_matrix_sparse
!  SYNOPSIS

subroutine prune_density_matrix_sparse(density_matrix_sparse, density_matrix_con, &
     n_compute, i_basis)

!  PURPOSE
!    The subroutine saves the density matrix components belongs to non-zero basis functions
!    to density_matrix_con.
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
  integer:: i_basis(n_compute)

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



  integer :: i_compute_2, i_compute_1, i_compute, i_center
  integer :: i_c1, i_c2
  integer :: i_index_real
  integer :: i_start, i_end, i_place, i_basis_2, i_basis_1, i_cell
  integer :: i_cell_2, i_cell_1, i_cell_2m1
  integer :: offset(n_cells)
  integer :: offset_end(n_cells) 

!NEC_CB
  integer::comp2basis(n_compute)
  integer::comp2cell(n_compute)

  
  density_matrix_con = 0.d0


  select case(packed_matrix_format)
     
  case(PM_index)

     if(n_periodic == 0)then

        do i_c2 = 1, n_compute, 1

           i_start =  index_hamiltonian(1,1, i_basis(i_c2))
           i_end   =  index_hamiltonian(2,1, i_basis(i_c2))

           if(i_end<i_start) cycle ! otherways an undefined i_index_real is used!


           do i_c1 = 1, i_c2, 1


              i_basis_2 = i_basis(i_c1)

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

!NEC_CB Do not read multi-indirectly addressed indices multiple times
        do i_compute = 1, n_compute, 1
          comp2basis(i_compute)=Cbasis_to_basis(i_basis(i_compute))
          i_center = Cbasis_to_center(i_basis(i_compute))
          comp2cell(i_compute)=center_to_cell(i_center)
        end do

        do i_compute_2 = 1, n_compute

           i_basis_2 = comp2basis(i_compute_2)
           i_cell_2 = comp2cell(i_compute_2)

           ! Now that we know i_cell_2, get for each i_cell_1 the cell-index
           ! corresponding to cell_2 - cell_1 (i_cell_2m1).
           ! As we also know i_basis_2, also get the corresponding part of
           ! density_matrix_sparse.
           do i_cell_1 = 1, n_cells 
              i_cell_2m1 = position_in_hamiltonian(i_cell_2, i_cell_1)
              offset(i_cell_1)     = index_hamiltonian(1,i_cell_2m1, i_basis_2)
              offset_end(i_cell_1) = index_hamiltonian(2,i_cell_2m1, i_basis_2)
           end do

           do i_compute_1 = 1, n_compute
             
               i_basis_1 = comp2basis(i_compute_1)
               i_cell_1    = comp2cell(i_compute_1)

               if(i_basis_1 <= i_basis_2)then

               place_2: do i_place = offset(i_cell_1), offset_end(i_cell_1)
                 
                  if (column_index_hamiltonian( i_place) == i_basis_1) then
                     !NEC_CB density_matrix_con is only needed as upper triangle
                     if (i_compute_2.le.i_compute_1) then
                       density_matrix_con(i_compute_2, i_compute_1) = density_matrix_sparse(i_place)
                     else
                       density_matrix_con(i_compute_1, i_compute_2) = density_matrix_sparse(i_place)
                     endif
                     exit place_2
                 
                  else if(column_index_hamiltonian( i_place) > i_basis_1)then

                     exit place_2  

                  end if
               end do place_2
               offset(i_cell_1) = i_place 
               end if

            end do ! i_compute_1
         end do ! i_compute_2
        
! The old routine just in case
!
!        do i_c2 = 1, n_compute, 1
!           do i_c1 = 1, i_c2, 1
!
!
!              if(Cbasis_to_basis(i_basis(i_c1)) <= Cbasis_to_basis(i_basis(i_c2)))then
!
!                 i_compute_1 = i_basis(i_c2)
!                 i_compute_2 = i_basis(i_c1)
!              else
!                 i_compute_1 = i_basis(i_c1)
!                 i_compute_2 = i_basis(i_c2)
!              end if
!
!
!           
!              i_index_real =  find_position_in_hamiltonian( i_compute_1, i_compute_2) 
!              density_matrix_con(i_c1, i_c2) = density_matrix_sparse(i_index_real)
!
!
!           end do
!        end do

     end if ! n_periodic == 0 - else


  end select ! packed matrix format










end subroutine prune_density_matrix_sparse
!******
