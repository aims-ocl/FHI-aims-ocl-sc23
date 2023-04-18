!****s* FHI-aims/evaluate_first_order_S
!  NAME
!    evaluate_first_order_S
!  SYNOPSIS

subroutine evaluate_first_order_S_right(first_order_S, partition_tab, &
           n_points,n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave)

!  PURPOSE
!    calculate the first-order overlap matrix elements.

!  shanghui,2012.05.04
!  USES

  use dimensions
  !use species_data
  !use geometry
  !use runtime_choices
  use basis  !basis_atom()

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave

  real*8, dimension(3, n_atoms, n_basis,n_basis), intent(inout) :: first_order_S

!  INPUTS
!  o  rho -- electron density
!  o  wave --
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_sq_global -- (distance to atoms)**2
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  fire_order_S


  integer :: i_point, i_coords,i_basis,j_basis,i_compute,j_compute
  real*8 :: point_term

  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1
      
      i_basis=i_basis_index(i_compute)
      j_basis=i_basis_index(j_compute)

           if(i_basis.ne.j_basis) then ! here we use transition conservation

       do i_point = 1, n_points, 1
        do i_coords = 1, 3, 1

           first_order_S(i_coords, basis_atom(j_basis),i_basis,j_basis) =    &
           first_order_S(i_coords, basis_atom(j_basis),i_basis,j_basis) -    & 
           partition_tab(i_point)* & 
           gradient_basis_wave(j_compute,i_coords,i_point)*wave(i_compute,i_point) 

!           endif

        end do
       end do

          endif 
  enddo
  enddo

end subroutine evaluate_first_order_S_right
!******
