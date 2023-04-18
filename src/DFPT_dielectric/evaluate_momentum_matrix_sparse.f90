!****s* FHI-aims/evaluate_momentum_matrix_sparse
!  NAME
!    evaluate_momentum_matrix_sparse
!  SYNOPSIS

subroutine evaluate_momentum_matrix_sparse(n_points, partition_tab, &
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave,   & 
           j_coord,  &
           momentum_matrix_sparse )

!  PURPOSE
!  calculate the first-order overlap matrix elements for phonon_gamma
!  shanghui,2013.12.30

!  first_order_S_sparse for phonon_reduce_memory 
!  shanghui, 2015.07.30

!  momentum_matrix_sparse for DFPT_dielectric
!  shanghui, 2015.11.27
!  USES

  use dimensions
  use basis  !basis_atom()
  use pbc_lists


!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave
  integer , intent(in) :: j_coord

  real*8, dimension(n_hamiltonian_matrix_size_no_symmetry), intent(inout) :: momentum_matrix_sparse
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!
!  OUTPUT
!  fire_order_S


  integer :: i_compute,j_compute, i_basis,j_basis, i_basis_uc, j_basis_uc
  integer :: i_cell, j_cell
  integer :: i_point, i_place
  !real*8, external :: ddot


  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)
     j_cell = center_to_cell(Cbasis_to_center(j_basis))

     do i_place = &
        index_hamiltonian_no_symmetry(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
        index_hamiltonian_no_symmetry(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

        if( column_index_hamiltonian_no_symmetry( i_place) == j_basis_uc)then
        do i_point = 1, n_points, 1

            momentum_matrix_sparse(i_place) = &
            momentum_matrix_sparse(i_place) + &
            partition_tab(i_point)* &
            wave(i_compute,i_point)*gradient_basis_wave(j_compute,j_coord,i_point)  

        enddo ! i_point
        endif ! column

     enddo ! i_place 

  enddo  ! j_compute
  enddo  ! i_compute 

end subroutine evaluate_momentum_matrix_sparse
!******
