!****s* FHI-aims/evaluate_first_order_H
!  NAME
!    evaluate_first_order_H
!  SYNOPSIS

subroutine evaluate_first_order_H(first_order_H,n_points, &
           partition_tab, dist_tab_global, dir_tab_global,       &
           H_times_psi, n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           first_order_rho,v_hartree_gradient,dVxc_drho)

!  PURPOSE
!    calculate the first-order Hamiltion matrix elements.
!    five terms:
!   (1) <X0| Z_I(R-r)/|r-R_I|^3 |X0>   ----------
!   (2) <X0| Int{rho(1)/|r-r'|} |X0>     Vscf(1)
!   (3) <X0| dVxc/drho * rho(1) |X0>   ----------
!   (4) <X0| Hks(0)             |X1>
!   (5) <X1| Hks(0)             |X0>

!  shanghui,2012.03.05
!  shanghui,2012.04.23
!  shanghui,2012.05.29 : to complie 
!  USES

  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use basis 

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_global
  real*8, dimension(3,n_atoms, n_points), intent(in) :: dir_tab_global

  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham, n_points),intent(in) :: wave
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi

  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: v_hartree_gradient
  real*8, dimension(n_points), intent(in) :: dVxc_drho


  real*8, dimension(3, n_atoms, n_basis,n_basis), intent(inout) :: first_order_H
!  INPUTS
!  o  rho -- electron density
!  o  wave --
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_global -- (distance to atoms)**1
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  first_order_H 


  integer :: i_point, i_atom, i_coords,i_basis,j_basis,i_compute,j_compute
  real*8 :: point_term
  real*8 :: atomic_term

  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

     i_basis=i_basis_index(i_compute)
     j_basis=i_basis_index(j_compute)

  do i_point = 1, n_points, 1
     point_term = partition_tab(i_point) * wave(i_compute,i_point) * wave(j_compute,i_point) 
     do i_atom = 1, n_atoms, 1
        atomic_term = -species_z(species(i_atom)) * point_term &
                      / (dist_tab_global(i_atom, i_point))**3
        do i_coords = 1, 3, 1
           first_order_H(i_coords, i_atom,i_basis,j_basis) =    &
           first_order_H(i_coords, i_atom,i_basis,j_basis) +    & 
           atomic_term * dir_tab_global(i_coords, i_atom, i_point)   +    & !(1) 
          point_term*v_hartree_gradient(i_coords, i_atom, i_point)   +    & !(2)
          point_term*dVxc_drho(i_point)*first_order_rho(i_coords,i_atom, i_point) !(3)

    ! gradient_basis_wave(i_basis,i_coords,i_point)*H_times_psi(j_basis,i_point)+ & !(4)
    ! gradient_basis_wave(j_basis,i_coords,i_point)*H_times_psi(i_basis,i_point) ! (5)
        end do
     end do

        do i_coords = 1, 3, 1
           first_order_H(i_coords, basis_atom(i_basis),i_basis,j_basis) =    &
           first_order_H(i_coords, basis_atom(i_basis),i_basis,j_basis) -    &
           partition_tab(i_point)* &  
           gradient_basis_wave(i_compute,i_coords,i_point)*H_times_psi(j_compute,i_point) 
  
           first_order_H(i_coords, basis_atom(j_basis),i_basis,j_basis) =    &
           first_order_H(i_coords, basis_atom(j_basis),i_basis,j_basis) -    &
           partition_tab(i_point)* &  
           gradient_basis_wave(j_compute,i_coords,i_point)*H_times_psi(i_compute,i_point) 
        enddo


  end do ! n_points

  enddo
  enddo
  
end subroutine evaluate_first_order_H
!******
