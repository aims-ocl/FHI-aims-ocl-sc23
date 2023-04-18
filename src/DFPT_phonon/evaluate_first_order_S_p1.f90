!****s* FHI-aims/evaluate_first_order_S_p1
!  NAME
!    evaluate_first_order_S_p1
!  SYNOPSIS

subroutine evaluate_first_order_S_p1(  & 
           n_points, partition_tab, &
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           first_order_S_sparse)

!  PURPOSE
!    calculate the first-order overlap matrix elements for phonon_gamma

!  shanghui,2013.12.30
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


  real*8, dimension(3, n_centers_in_sc_DFPT, n_hamiltonian_matrix_size), intent(inout) :: & 
                                  first_order_S_sparse
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!  o  i_q_point -- the q point we works on
!
!  OUTPUT
!  fire_order_S


  integer :: i_point, i_place, i_coord,i_basis,j_basis,i_compute,j_compute 
  integer :: i_center_trans, j_center_trans
  integer :: i_basis_uc,j_basis_uc
  integer :: i_atom,j_atom
  integer :: i_cell_in_hamiltonian,j_cell_in_hamiltonian
  integer :: i_cell_in_sc_DFPT,j_cell_in_sc_DFPT



 
  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute) 
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_atom     = Cbasis_to_atom(i_basis)
     i_cell_in_hamiltonian     = center_to_cell(Cbasis_to_center(i_basis))
     i_cell_in_sc_DFPT         = center_in_sc_DFPT_to_cell_in_sc_DFPT( & 
                                 center_to_center_in_sc_DFPT(Cbasis_to_center(i_basis)) ) 

  do j_compute=1,n_compute_c,1
     j_basis=i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)


  if(j_basis_uc <= i_basis_uc) then
       j_atom = Cbasis_to_atom(j_basis) 
       j_cell_in_hamiltonian = center_to_cell(Cbasis_to_center(j_basis))
       j_cell_in_sc_DFPT     = center_in_sc_DFPT_to_cell_in_sc_DFPT( & 
                               center_to_center_in_sc_DFPT(Cbasis_to_center(j_basis)) )
 

  do i_place = &
        index_hamiltonian(1,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc),  &
        index_hamiltonian(2,position_in_hamiltonian(i_cell_in_hamiltonian,j_cell_in_hamiltonian), i_basis_uc)

  if( column_index_hamiltonian( i_place) == j_basis_uc)then

         !----------get the i_place related basis, and basis related center_trans-----
          i_center_trans = cell_and_atom_to_center_sc_DFPT( & 
                           cell_diff_sc_DFPT(i_cell_in_sc_DFPT,j_cell_in_sc_DFPT),i_atom) 

          j_center_trans = j_atom 
         !---------end get the i_place realated basis--------------------------------


       do i_point = 1, n_points, 1
  
        if(i_basis.ne.j_basis) then 
        do i_coord = 1, 3, 1

          first_order_S_sparse(i_coord, i_center_trans, i_place) = &
          first_order_S_sparse(i_coord, i_center_trans, i_place) - &
          partition_tab(i_point)*  &  
          gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)  


          first_order_S_sparse( i_coord, j_center_trans, i_place) = & 
          first_order_S_sparse( i_coord, j_center_trans, i_place) - &
          partition_tab(i_point)* &
          gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point)  
 
 
        enddo ! i_coord
        endif ! i_basis 

       enddo ! i_point

  endif ! column
  enddo ! i_place 
  endif ! j_basis_uc <i_basis_uc


  enddo ! j_compute 
  enddo ! i_compute



end subroutine evaluate_first_order_S_p1
!******
