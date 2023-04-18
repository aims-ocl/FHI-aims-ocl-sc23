!****s* FHI-aims/evaluate_Omega_part_1_sparse
!  NAME
!    evaluate_Omega_part_1_sparse
!  SYNOPSIS

subroutine evaluate_Omega_part_1_sparse(n_points, partition_tab, &
           n_compute_c, i_basis_index,  & 
           wave, coords_npoints,  & 
           j_coord,  &
           Omega_part_1_sparse )


!  PURPOSE
!  calculate Omega_part_1_sparse
!  USES

  use dimensions
  use pbc_lists
  use geometry, only : lattice_vector

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(3,n_points), intent(in) :: coords_npoints
  integer , intent(in) :: j_coord

  real*8, dimension(n_hamiltonian_matrix_size_no_symmetry), intent(inout) :: Omega_part_1_sparse
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!
!  OUTPUT
!  Omega_part_1_sparse


  integer :: i_compute,j_compute, i_basis,j_basis, i_basis_uc, j_basis_uc
  integer :: i_cell, j_cell
  integer :: i_point, i_place
  real*8  :: R_j_cell 


  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)
     j_cell = center_to_cell(Cbasis_to_center(j_basis))
     R_j_cell = cell_index(j_cell, 1)*lattice_vector(j_coord,1) + &
                cell_index(j_cell, 2)*lattice_vector(j_coord,2) + &
                cell_index(j_cell, 3)*lattice_vector(j_coord,3) 

     do i_place = &
        index_hamiltonian_no_symmetry(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
        index_hamiltonian_no_symmetry(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

     if( column_index_hamiltonian_no_symmetry( i_place) == j_basis_uc)then

       do i_point = 1, n_points, 1

           Omega_part_1_sparse(i_place) = &
           Omega_part_1_sparse(i_place) - &
           partition_tab(i_point)* &                               
           wave(i_compute,i_point)*wave(j_compute,i_point)* & 
           (coords_npoints(j_coord,i_point) - R_j_cell )  

       enddo ! i_point

     endif ! column
     enddo ! i_place 

  enddo  ! j_compute
  enddo  ! i_compute 


end subroutine evaluate_Omega_part_1_sparse
!******
