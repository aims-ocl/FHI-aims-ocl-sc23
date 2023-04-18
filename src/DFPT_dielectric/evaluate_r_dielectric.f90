!****s* FHI-aims/evaluate_r_dielectric
!  NAME
!    evaluate_r_dielectric
!  SYNOPSIS

subroutine evaluate_r_dielectric& 
         ( n_points, partition_tab, & 
           n_compute_c, i_basis_index,  & 
           wave,            &
           j_coord,         &
           coords_npoints,  &
           r_sparse )

!  PURPOSE
!    calculate the r Hamiltion matrix elements for phonon_gamma
!    four terms:
!   (1) <X0| r   |X0> 

!
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 
  use pbc_lists

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in)    :: wave
  integer , intent(in) :: j_coord
  real*8, dimension(3, n_points),intent(in)    :: coords_npoints
  real*8, dimension(n_hamiltonian_matrix_size), intent(inout) :: r_sparse
!  INPUTS
!  o  n_points -- number of grid points in this grid batch
!  o  partition_tab -- values of partition function
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!  o  H_times_psi           -- H*psi  
!  o  first_order_rho       -- rho_tot(1)
!  o  first_order_potential -- V_free_hartree(1)+delta_V_hartree(1)
!  o  dVxc_drho             -- d Vxc/ drho
!
!  OUTPUT
!  first_order_H 


  integer :: i_compute,j_compute, i_basis,j_basis, i_basis_uc, j_basis_uc
  integer :: i_cell, j_cell
  integer :: i_point, i_place
  real*8  :: point_term

  real*8, external :: ddot


  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)

     if(j_basis_uc <= i_basis_uc) then
        j_cell = center_to_cell(Cbasis_to_center(j_basis))

     do i_place = &
        index_hamiltonian(1,position_in_hamiltonian(i_cell,j_cell), i_basis_uc),  &
        index_hamiltonian(2,position_in_hamiltonian(i_cell,j_cell), i_basis_uc)

     if( column_index_hamiltonian( i_place) == j_basis_uc)then

       do i_point = 1, n_points, 1

           point_term = partition_tab(i_point)    &
                    * wave(i_compute,i_point) * wave(j_compute,i_point)


!-------------------local term: <mu| -r |nu>---------------------
            r_sparse(i_place) = &
            r_sparse(i_place) - &
                  point_term*coords_npoints(j_coord,i_point)    


       enddo ! i_point

     endif ! column
     enddo ! i_place 
     endif ! j_basis_uc <i_basis_uc


  enddo  ! j_compute
  enddo  ! i_compute 


end subroutine evaluate_r_dielectric
!******
