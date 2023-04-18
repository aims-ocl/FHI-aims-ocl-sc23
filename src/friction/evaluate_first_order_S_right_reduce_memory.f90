!****s* FHI-aims/evaluate_first_order_S
!  NAME
!    evaluate_first_order_S
!  SYNOPSIS

subroutine evaluate_first_order_S_right_reduce_memory(first_order_S, partition_tab, &
           n_points,n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave,j_atom,j_coord)

!  PURPOSE
!    calculate the first-order overlap matrix elements.

!  shanghui,2012.05.04
!  USES

  use dimensions
  use basis  !basis_atom()

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave

  real*8, dimension(n_basis,n_basis), intent(inout) :: first_order_S
  integer , intent(in) :: j_atom, j_coord
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


  integer :: i_point,i_basis,j_basis,i_compute,j_compute
  real*8 :: point_term

  
  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1
      
      i_basis=i_basis_index(i_compute)
      j_basis=i_basis_index(j_compute)

      if(i_basis.ne.j_basis) then ! here we use transition conservation

         do i_point = 1, n_points, 1

           !if(j_atom.eq.basis_atom(i_basis)) then
           !first_order_S(i_basis,j_basis) =    &
           !first_order_S(i_basis,j_basis) -    & 
           !partition_tab(i_point)*  &  
           !gradient_basis_wave(i_compute,j_coord,i_point)*wave(j_compute,i_point) 
           !endif                     

           if(j_atom.eq.basis_atom(j_basis)) then 
           first_order_S(i_basis,j_basis) =    &
           first_order_S(i_basis,j_basis) -    & 
           partition_tab(i_point)* & 
           gradient_basis_wave(j_compute,j_coord,i_point)*wave(i_compute,i_point) 
           endif 

         end do

      endif 
  enddo
  enddo




end subroutine evaluate_first_order_S_right_reduce_memory
!******
