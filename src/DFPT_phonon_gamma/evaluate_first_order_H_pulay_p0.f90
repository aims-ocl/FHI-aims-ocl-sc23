!****s* FHI-aims/evaluate_first_order_H_pulay_p0
!  NAME
!    evaluate_first_order_H_pulay_p0
!  SYNOPSIS

subroutine evaluate_first_order_H_pulay_p0(first_order_H_pulay, &
           n_points,partition_tab, n_compute_c, i_basis_index,  &
           H_times_psi, gradient_basis_wave)

!  PURPOSE
!    calculate the first-order Hamiltion pulay matrix elements for phonon_gamma
!    two terms:
!   (4) <X0| Hks(0)             |X1>
!   (5) <X1| Hks(0)             |X0>

!  shanghui,2013.12.30 
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
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi


  real*8, dimension(3, n_atoms, n_Cbasis,n_Cbasis), intent(inout) :: first_order_H_pulay
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  first_order_H_pulay 


  integer :: i_point,  i_coords,i_basis,j_basis,i_compute,j_compute
  

  
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

  do i_point = 1, n_points, 1

        do i_coords = 1, 3, 1
          first_order_H_pulay(i_coords, Cbasis_to_atom(i_basis),i_basis,j_basis) =    &
          first_order_H_pulay(i_coords, Cbasis_to_atom(i_basis),i_basis,j_basis) -    &
          partition_tab(i_point)* &
         gradient_basis_wave(i_compute,i_coords,i_point)*H_times_psi(j_compute,i_point)

          first_order_H_pulay(i_coords, Cbasis_to_atom(j_basis),i_basis,j_basis) =    &
          first_order_H_pulay(i_coords, Cbasis_to_atom(j_basis),i_basis,j_basis) -    &
          partition_tab(i_point)* &
         gradient_basis_wave(j_compute,i_coords,i_point)*H_times_psi(i_compute,i_point)
        enddo


  end do ! n_points

  enddo
  enddo
 
 
end subroutine evaluate_first_order_H_pulay_p0
!******
