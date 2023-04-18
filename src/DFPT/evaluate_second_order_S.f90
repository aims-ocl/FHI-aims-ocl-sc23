!****s* FHI-aims/evaluate_second_order_S
!  NAME
!    evaluate_second_order_S
!  SYNOPSIS

subroutine evaluate_second_order_S(second_order_S, partition_tab, &
           n_points,n_compute_c, i_basis_index, index_hessian, & 
           wave, gradient_basis_wave,hessian_basis_wave, & 
           atom_point, partition_deriv_local )

!  PURPOSE
!    calculate the second overlap matrix elements.
!  shanghui,2013.02.26

!    add  moving grid effect 
!  shanghui,2013.06.25
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
  integer , intent(in) :: index_hessian(3,3)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,6,n_points), intent(in) :: hessian_basis_wave
  integer, dimension(n_max_batch_size), intent(in) ::  atom_point
  real*8, dimension(3,n_atoms,n_points), intent(in) :: partition_deriv_local

  real*8, dimension(3, n_atoms,3,n_atoms, n_basis,n_basis), intent(inout) :: second_order_S
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  second_order_S


  integer :: i_point, i_coord,j_coord,i_basis,j_basis,i_compute,j_compute



   do i_point = 1, n_points, 1


    do i_compute=1,n_compute_c,1
    do j_compute=1,n_compute_c,1
      
      i_basis=i_basis_index(i_compute)
      j_basis=i_basis_index(j_compute)

       if(i_basis.ne.j_basis) then ! here we use transition conservation

        do i_coord = 1, 3, 1
        do j_coord = 1, 3, 1


           !------------(1)d<XI|XJ>/dRI dRI---------------------------------------
           if(atom_point(i_point).ne.basis_atom(i_basis) ) then !fixed  

           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(i_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(i_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(i_basis),i_point)* & 
           (-gradient_basis_wave(i_compute,i_coord,i_point))* & 
           wave(j_compute,i_point) + & 
           partition_tab(i_point)*  & 
           hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point)* & 
           wave(j_compute,i_point) 
         
           else !moving

           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(i_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(i_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(i_basis),i_point)* & 
           (-gradient_basis_wave(i_compute,i_coord,i_point))* & 
           wave(j_compute,i_point) + & 
           partition_tab(i_point)*  &
           gradient_basis_wave(i_compute,i_coord,i_point)*                &
           (-gradient_basis_wave(j_compute,j_coord,i_point))

           endif

           !------------(2)d<XI|XJ>/dRI dRJ---------------------------------------
           if(atom_point(i_point).ne.basis_atom(j_basis) ) then !fixed  

           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(j_basis),i_point)* & 
           (-gradient_basis_wave(i_compute,i_coord,i_point))* & 
           wave(j_compute,i_point) + & 
           partition_tab(i_point)*  &
           gradient_basis_wave(i_compute,i_coord,i_point)*                &
           gradient_basis_wave(j_compute,j_coord,i_point)

           else !moving 

           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(i_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(j_basis),i_point)* & 
           (-gradient_basis_wave(i_compute,i_coord,i_point))* & 
           wave(j_compute,i_point) + & 
           partition_tab(i_point)*  & 
           (-hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point))* & 
           wave(j_compute,i_point) 

           endif

           !------------(3)d<XI|XJ>/dRJ dRI---------------------------------------
           if(atom_point(i_point).ne.basis_atom(i_basis) ) then !fixed  
           second_order_S(j_coord,basis_atom(j_basis),i_coord,basis_atom(i_basis), & 
                          i_basis,j_basis) =    &
           second_order_S(j_coord,basis_atom(j_basis),i_coord,basis_atom(i_basis), & 
                          i_basis,j_basis) +    & 
           partition_deriv_local(i_coord,basis_atom(i_basis),i_point)* & 
           (-gradient_basis_wave(j_compute,j_coord,i_point))* & 
           wave(i_compute,i_point) + & 
           partition_tab(i_point)*  &  
           gradient_basis_wave(j_compute,j_coord,i_point)*                & 
           gradient_basis_wave(i_compute,i_coord,i_point)     

           else !moving 
           second_order_S(j_coord,basis_atom(j_basis),i_coord,basis_atom(i_basis), & 
                          i_basis,j_basis) =    &
           second_order_S(j_coord,basis_atom(j_basis),i_coord,basis_atom(i_basis), & 
                          i_basis,j_basis) +    & 
           partition_deriv_local(i_coord,basis_atom(i_basis),i_point)* & 
           (-gradient_basis_wave(j_compute,j_coord,i_point))* & 
           wave(i_compute,i_point) + & 
           partition_tab(i_point)*   & 
           (-hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point))* & 
           wave(i_compute,i_point)
           
           endif

           !------------(4)d<XI|XJ>/dRJ dRJ---------------------------------------
           if(atom_point(i_point).ne.basis_atom(j_basis) ) then !fixed 
 
           second_order_S(i_coord,basis_atom(j_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(j_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(j_basis),i_point)* & 
           (-gradient_basis_wave(j_compute,i_coord,i_point))* & 
           wave(i_compute,i_point) + & 
           partition_tab(i_point)*   & 
           hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point)* & 
           wave(i_compute,i_point)
 
           else !moving  
           second_order_S(i_coord,basis_atom(j_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) =    &
           second_order_S(i_coord,basis_atom(j_basis),j_coord,basis_atom(j_basis), &
                          i_basis,j_basis) +    &
           partition_deriv_local(j_coord,basis_atom(j_basis),i_point)* & 
           (-gradient_basis_wave(j_compute,i_coord,i_point))* & 
           wave(i_compute,i_point) + & 
           partition_tab(i_point)*  &  
           gradient_basis_wave(j_compute,j_coord,i_point)*                & 
           (-gradient_basis_wave(i_compute,i_coord,i_point))     
            
           endif

        enddo
        enddo

       endif

    end do
    end do

  enddo ! n_points


end subroutine evaluate_second_order_S
!******
