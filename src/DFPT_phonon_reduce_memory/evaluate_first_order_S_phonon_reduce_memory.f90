!****s* FHI-aims/evaluate_first_order_S_phonon_reduce_memory
!  NAME
!    evaluate_first_order_S_phonon_reduce_memory
!  SYNOPSIS

subroutine evaluate_first_order_S_phonon_reduce_memory(n_points, partition_tab, &
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, coords_npoints,  & 
           i_q_point, j_atom, j_coord,  &
           first_order_S_sparse )



!  PURPOSE
!  calculate the first-order overlap matrix elements for phonon_gamma
!  shanghui,2013.12.30

!  first_order_S_sparse for phonon_reduce_memory 
!  shanghui, 2015.07.30
!  USES

  use dimensions
  use basis  !basis_atom()
  use geometry, only: recip_lattice_vector
  use pbc_lists


!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave
  real*8, dimension(3,n_points), intent(in) :: coords_npoints
  integer , intent(in) :: i_q_point
  integer , intent(in) :: j_atom
  integer , intent(in) :: j_coord

  complex*16, dimension(n_hamiltonian_matrix_size), intent(inout) :: first_order_S_sparse
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
  real*8  :: Gr(3)
  complex*16  ::  exp_iqr(3) 
  real*8, external :: ddot



  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)

  if(Cbasis_to_center(i_basis).ne.Cbasis_to_center(j_basis)) then ! here we use transition conservation


     if(j_basis_uc <= i_basis_uc) then
        j_cell = center_to_cell(Cbasis_to_center(j_basis))

     do i_place = &
        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

     if( column_index_hamiltonian( i_place) == j_basis_uc)then

       do i_point = 1, n_points, 1

           Gr(1) = ddot(3,recip_lattice_vector(1:3,1),1, coords_npoints(1:3,i_point),1)
           Gr(2) = ddot(3,recip_lattice_vector(1:3,2),1, coords_npoints(1:3,i_point),1)
           Gr(3) = ddot(3,recip_lattice_vector(1:3,3),1, coords_npoints(1:3,i_point),1)

           exp_iqr(1) = exp(-(0,1)*Gr(1)*k_point_list(i_q_point,1))
           exp_iqr(2) = exp(-(0,1)*Gr(2)*k_point_list(i_q_point,2))
           exp_iqr(3) = exp(-(0,1)*Gr(3)*k_point_list(i_q_point,3))

          if(j_atom.eq.Cbasis_to_atom(i_basis)) then
            first_order_S_sparse(i_place) = &
            first_order_S_sparse(i_place) - &
            partition_tab(i_point)* &
            gradient_basis_wave(i_compute,j_coord,i_point)*wave(j_compute,i_point)* &
            k_phase(i_cell,i_q_point)* &
            exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
          endif 

          if(j_atom.eq.Cbasis_to_atom(j_basis)) then
            first_order_S_sparse(i_place) = &
            first_order_S_sparse(i_place) - &
            partition_tab(i_point)* &                               
            gradient_basis_wave(j_compute,j_coord,i_point)*wave(i_compute,i_point)* & 
            k_phase(j_cell,i_q_point)* &
            exp_iqr(1)*exp_iqr(2)*exp_iqr(3)
          endif

       enddo ! i_point

     endif ! column
     enddo ! i_place 
     endif ! j_basis_uc <i_basis_uc


  endif  ! i_center .ne. j_center
  enddo  ! j_compute
  enddo  ! i_compute 



end subroutine evaluate_first_order_S_phonon_reduce_memory
!******
