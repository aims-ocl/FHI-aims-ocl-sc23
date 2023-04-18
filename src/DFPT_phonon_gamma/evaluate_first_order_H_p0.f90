!****s* FHI-aims/evaluate_first_order_H_p0
!  NAME
!    evaluate_first_order_H_p0
!  SYNOPSIS

subroutine evaluate_first_order_H_p0(first_order_H,n_points, partition_tab, & 
           n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           H_times_psi,first_order_rho,first_order_potential,dVxc_drho)

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements for phonon_gamma
!    four terms:
!   (1) <X0| free_V_hartree(1) |X0> 
!   (2) <X0| delta_V_hartree(1) |X0> 
!   (2) <X0| dVxc/drho * rho(1) |X0>   
!   (3) <X0| Hks(0)             |X1>
!   (4) <X1| Hks(0)             |X0>

!   shanghui  2013.12.30
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
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi

  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in) :: dVxc_drho


  real*8, dimension(3, n_atoms, n_Cbasis,n_Cbasis), intent(inout) :: first_order_H
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


  integer :: i_point, i_atom, i_coords,i_basis,j_basis,i_compute,j_compute
  real*8 :: point_term
 



  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

  do i_point = 1, n_points, 1
     point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 

     do i_atom = 1, n_atoms, 1
        do i_coords = 1, 3, 1

           first_order_H(i_coords, i_atom,i_basis,j_basis) =    &
           first_order_H(i_coords, i_atom,i_basis,j_basis) +    & 
           point_term*first_order_potential(i_coords, i_atom, i_point)   +    & !(1) + (2)
           point_term*dVxc_drho(i_point)*first_order_rho(i_coords,i_atom, i_point) !(3)

        end do
     end do

        do i_coords = 1, 3, 1
           first_order_H(i_coords, Cbasis_to_atom(i_basis),i_basis,j_basis) =    &
           first_order_H(i_coords, Cbasis_to_atom(i_basis),i_basis,j_basis) -    &
           partition_tab(i_point)* &                                                !(4)
           gradient_basis_wave(i_compute,i_coords,i_point)*H_times_psi(j_compute,i_point) 
  
           first_order_H(i_coords, Cbasis_to_atom(j_basis),i_basis,j_basis) =    &
           first_order_H(i_coords, Cbasis_to_atom(j_basis),i_basis,j_basis) -    &  !(5)
           partition_tab(i_point)* &  
           gradient_basis_wave(j_compute,i_coords,i_point)*H_times_psi(i_compute,i_point) 
        enddo
 

  end do ! n_points

  enddo
  enddo


end subroutine evaluate_first_order_H_p0
!******
