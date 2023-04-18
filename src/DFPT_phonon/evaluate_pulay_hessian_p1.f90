!****s* FHI-aims/evaluate_pulay_hessian_p1
!  NAME
!    evaluate_pulay_hessian_p1
!  SYNOPSIS

subroutine evaluate_pulay_hessian_p1(  & 
           partition_tab, &
           n_points, n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave,             & 
           index_hessian, hessian_basis_wave,  & 
           H_times_psi, &
           first_order_rho,first_order_potential, dVxc_drho,             &
           H_times_gradient_psi,  &
           density_matrix_compute, first_order_density_matrix_compute, &
           energy_density_matrix_compute, first_order_energy_density_matrix_compute, &
           pulay_hessian )



!  PURPOSE
!    calculate pulay_hessian

!  shanghui, 2015
!  USES

  use dimensions
  use basis  !basis_atom()
  use pbc_lists
  use mpi_tasks, only: check_allocation

!  ARGUMENTS

  implicit none
  !note: for grid related basis, use n_max_compute_ham, 
  !      for matrix, use n_compute_c
  !-------------(1) for_all-------------------------------------------------
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave
  !-------------(2) add for S2------------------------------------------------
  integer , intent(in) :: index_hessian(3,3)
  real*8, dimension(n_max_compute_ham,6,n_points), intent(in) :: hessian_basis_wave

  !-------------(3) add for H1_pulay and H2_pulay--------------------------------------------
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi
  real*8, dimension(3, n_centers_in_sc_DFPT, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_centers_in_sc_DFPT, n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in) :: dVxc_drho

  !-------------(4) add for H2_pulay-------------------------------------------------------
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: H_times_gradient_psi

  
  !-------------(5) input compute matrix-----------------
  real*8, dimension(n_compute_c,n_compute_c),intent(in) :: density_matrix_compute
  real*8, dimension(3,n_centers_in_sc_DFPT,n_compute_c,n_compute_c),intent(in) :: &
                   first_order_density_matrix_compute

  real*8, dimension(n_compute_c,n_compute_c),intent(in) :: energy_density_matrix_compute
  real*8, dimension(3,n_centers_in_sc_DFPT,n_compute_c,n_compute_c),intent(in) :: &
                   first_order_energy_density_matrix_compute

 
  !-------------(6) output pulay_hessian--------------------------------   
  real*8, dimension(3,n_centers_in_sc_DFPT,3,n_atoms),intent(out) :: pulay_hessian



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


  integer :: i_point, i_place
  integer :: i_coord , j_coord, k_coord
  integer :: i_compute, j_compute
  integer :: i_basis,j_basis
  integer :: i_basis_uc,j_basis_uc
  integer :: i_cell, j_cell, k_cell
  integer ::  k_center 
  integer :: i_center_trans, j_center_trans, k_center_trans

  integer :: k_cell_trans, k_atom
  integer :: j_center, i_center, i_atom, i_cell_trans

  real*8 :: point_term

  integer :: info

  real*8, dimension(:,:),  allocatable :: temp_first_order_S
  real*8, dimension(:,:),  allocatable :: temp_first_order_H_pulay

  real*8, dimension(:,:,:,:),  allocatable :: temp_second_order_S
  real*8, dimension(:,:,:,:),  allocatable :: temp_second_order_H_pulay


  allocate(temp_first_order_S(3,n_centers_in_sc_DFPT),stat=info)
  call check_allocation(info,'temp_first_order_S            ')
  allocate(temp_first_order_H_pulay(3,n_centers_in_sc_DFPT),stat=info)
  call check_allocation(info,'temp_first_order_H_pulay            ')

  allocate(temp_second_order_S(3,n_centers_in_sc_DFPT,3,n_centers_in_sc_DFPT),stat=info)
  call check_allocation(info,'temp_second_order_S            ')
  allocate(temp_second_order_H_pulay(3,n_centers_in_sc_DFPT,3,n_centers_in_sc_DFPT),stat=info)
  call check_allocation(info,'temp_second_order_H_pulay            ')



  do i_compute=1,n_compute_c,1
     i_basis=i_basis_index(i_compute) 
     i_basis_uc = Cbasis_to_basis(i_basis)

  do j_compute=1,n_compute_c,1
     j_basis=i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)

         !----------get the i_place related basis, and basis related center_trans-----
          i_center_trans = center_to_center_in_sc_DFPT(Cbasis_to_center(i_basis)) 
          j_center_trans = center_to_center_in_sc_DFPT(Cbasis_to_center(j_basis))
         !---------end get the i_place realated basis--------------------------------

            temp_first_order_S(1:3,1:n_centers_in_sc_DFPT) = 0.0d0 
      temp_first_order_H_pulay(1:3,1:n_centers_in_sc_DFPT) = 0.0d0 

           temp_second_order_S(1:3,1:n_centers_in_sc_DFPT,1:3,1:n_centers_in_sc_DFPT) = 0.0d0 
     temp_second_order_H_pulay(1:3,1:n_centers_in_sc_DFPT,1:3,1:n_centers_in_sc_DFPT) = 0.0d0 


    do i_point = 1, n_points, 1
!---------------(1) first_order_S---------------------------

       if(i_basis.ne.j_basis) then
       do i_coord = 1, 3, 1
 
          temp_first_order_S(i_coord, i_center_trans) = &
          temp_first_order_S(i_coord, i_center_trans) - &
          partition_tab(i_point)*  &  
          gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)  


!         temp_first_order_S( i_coord, j_center_trans) = & 
!         temp_first_order_S( i_coord, j_center_trans) - &
!         partition_tab(i_point)* &
!         gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point)  


       end do
       endif


!---------------(2) second_order_S---------------------------
       do i_coord = 1, 3, 1
       do j_coord = 1, 3, 1
           
          !------------[1] d<XI|XJ>/dRI dRI-----------
          temp_second_order_S(i_coord, i_center_trans, j_coord, i_center_trans) = &
          temp_second_order_S(i_coord, i_center_trans, j_coord, i_center_trans) + &
          partition_tab(i_point)*  &  
          hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point)*wave(j_compute,i_point)  
 
          !------------[2] d<XI|XJ>/dRI dRJ---------------------------------------
          temp_second_order_S(i_coord, i_center_trans, j_coord, j_center_trans) = &
          temp_second_order_S(i_coord, i_center_trans, j_coord, j_center_trans) + &
          partition_tab(i_point)*  &  
          gradient_basis_wave(i_compute,i_coord,i_point)*gradient_basis_wave(j_compute,j_coord,i_point)  

!         !------------[3] d<XI|XJ>/dRJ dRI--------------------------------------- 
!         temp_second_order_S(j_coord, j_center_trans, i_coord, i_center_trans) = &
!         temp_second_order_S(j_coord, j_center_trans, i_coord, i_center_trans) + &
!         partition_tab(i_point)*  &  
!         gradient_basis_wave(j_compute,j_coord,i_point)*gradient_basis_wave(i_compute,i_coord,i_point)  
!
!         !------------[4] d<XI|XJ>/dRJ dRJ---------------------------------------
!         temp_second_order_S(i_coord, j_center_trans, j_coord, j_center_trans) = &
!         temp_second_order_S(i_coord, j_center_trans, j_coord, j_center_trans) + &
!         partition_tab(i_point)*  &  
!         hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point)*wave(i_compute,i_point)  
        
       enddo 
       enddo


!---------------(3) first_order_H_pulay---------------------------
       do i_coord = 1,3
          temp_first_order_H_pulay(i_coord, i_center_trans) =  &
          temp_first_order_H_pulay(i_coord, i_center_trans) - &
                     partition_tab(i_point)                              &
                     *gradient_basis_wave(i_compute,i_coord,i_point)     &
                     *H_times_psi(j_compute,i_point)

!         temp_first_order_H_pulay( i_coord, j_center_trans) = &
!         temp_first_order_H_pulay( i_coord, j_center_trans) - &
!                    partition_tab(i_point)                              &
!                    *gradient_basis_wave(j_compute,i_coord,i_point)     &
!                    *H_times_psi(i_compute,i_point)
!
!
       enddo




!---------------(4) second_order_H_pulay---------------------------
       do i_coord = 1, 3, 1

          !---------------[1] begin <dXu/dRI| H |Xv> u=I_atom, v=J_atom 
          point_term= -partition_tab(i_point) * &
          gradient_basis_wave(i_compute,i_coord,i_point) * wave(j_compute,i_point)

       do k_center = 1, n_centers_in_sc_DFPT

            k_center_trans = k_center 


          !---------------[1.1]delta(I,K)*<d^2Xu/dRI dRK| H |Xv> 
          if(k_center_trans.eq.i_center_trans) then
            do k_coord = 1,3 
          temp_second_order_H_pulay(i_coord,i_center_trans, k_coord,k_center_trans) =    &
          temp_second_order_H_pulay(i_coord,i_center_trans, k_coord,k_center_trans) +    &
          partition_tab(i_point)*  &
          hessian_basis_wave(i_compute,index_hessian(i_coord,k_coord),i_point)* &
          H_times_psi(j_compute,i_point)
            enddo 
          endif
      

          !--------------[1.2] <dXu/dRI| dH/dRk |Xv>
          do k_coord = 1, 3          
          temp_second_order_H_pulay(i_coord,i_center_trans,k_coord,k_center_trans) =    &
          temp_second_order_H_pulay(i_coord,i_center_trans,k_coord,k_center_trans) +    &
          point_term*first_order_potential(k_coord, k_center, i_point)   +    &
          point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_center, i_point)
          enddo


          !---------------[1.3] delta(K,J) <d Xu/dRI | H |d Xv/ dRK>
          if(k_center_trans.eq.j_center_trans) then
             do k_coord = 1,3 
          temp_second_order_H_pulay(i_coord,i_center_trans,k_coord,k_center_trans) =    &
          temp_second_order_H_pulay(i_coord,i_center_trans,k_coord,k_center_trans) +    &
          partition_tab(i_point)*                                   &
          gradient_basis_wave(i_compute,i_coord,i_point) *          &
          H_times_gradient_psi(j_compute,k_coord,i_point)
             enddo 
          endif
       



       enddo  ! k_center
       enddo  ! i_coord 




!---------------[2] begin <Xu| H |dXv/dRJ> u=I_atom, v=J_atom
!      do j_coord =1 , 3 
!         point_term=-partition_tab(i_point) * &
!         gradient_basis_wave(j_compute,j_coord,i_point) * wave(i_compute,i_point)
!
!      do k_center = 1, n_centers_in_sc_DFPT
!
!           k_center_trans = k_center
!
!
!         !---------------[2.1] delta(I,K)*<d Xu/dRK | H |d Xv/ dRJ> 
!         if(k_center_trans.eq.i_center_trans) then
!         do k_coord = 1,3
!         temp_second_order_H_pulay(j_coord,j_center_trans,k_coord,k_center_trans) =    &
!         temp_second_order_H_pulay(j_coord,j_center_trans,k_coord,k_center_trans) +    &
!         partition_tab(i_point)*                                   &
!         gradient_basis_wave(j_compute,j_coord,i_point) *          &
!         H_times_gradient_psi(i_compute,k_coord,i_point)
!         enddo 
!         endif
!
!         !--------------[2.2] <Xu| dH/dRk |d Xv/ dRJ>
!         do k_coord = 1,3
!         temp_second_order_H_pulay(j_coord,j_center_trans,k_coord,k_center_trans) =    &
!         temp_second_order_H_pulay(j_coord,j_center_trans,k_coord,k_center_trans) +    &
!         point_term*first_order_potential(k_coord, k_center, i_point)   +    &
!         point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_center, i_point)
!
!         enddo  
!
!
!         !---------------[2.3] delta(J,K)*<Xu| H |d^2Xv/dRJ dRK> 
!         if(k_center_trans.eq.j_center_trans) then 
!         do k_coord = 1,3 
!         temp_second_order_H_pulay(j_coord,j_center_trans, k_coord,k_center_trans) =    &
!         temp_second_order_H_pulay(j_coord,j_center_trans, k_coord,k_center_trans) +    &
!         partition_tab(i_point)*  &
!         hessian_basis_wave(j_compute,index_hessian(j_coord,k_coord),i_point)* &
!         H_times_psi(i_compute,i_point)
!         enddo 
!         endif
!
!
!
!      enddo  ! k_center
!      enddo  ! j_coord 


    end do  ! i_point

   do j_coord = 1 ,3 
   do j_center = 1, n_centers_in_sc_DFPT

       j_cell = center_in_sc_DFPT_to_cell_in_sc_DFPT(j_center)
       j_center_trans  = center_in_sc_DFPT_to_atom(j_center)
 
    do i_coord =1, 3
    do i_center = 1, n_centers_in_sc_DFPT
 
       i_atom = center_in_sc_DFPT_to_atom(i_center) 
       i_cell = center_in_sc_DFPT_to_cell_in_sc_DFPT(i_center)
       i_cell_trans = cell_diff_sc_DFPT(i_cell,j_cell)
       i_center_trans = cell_and_atom_to_center_sc_DFPT(i_cell_trans, i_atom) 

       pulay_hessian(i_coord,i_center_trans,j_coord,j_center_trans)= &
       pulay_hessian(i_coord,i_center_trans,j_coord,j_center_trans)+ &
        !-------so we have pulay_hessian--->(1) 
        ! sum_{uv}{ dm(1)_{uv}*h_pulay(1)_{uv} }  
        ! -sum_{uv}{ edm(1)_{uv}*s(1)_{uv} }
        2.0d0*first_order_density_matrix_compute(j_coord, j_center, i_compute,j_compute)*  &
        temp_first_order_H_pulay(i_coord,i_center) -   &
        2.0d0*first_order_energy_density_matrix_compute(j_coord, j_center,i_compute,j_compute)*  &
        temp_first_order_S(i_coord,i_center)  +  &
        !-------and the pulay_hessian   --->(2) 
        ! sum_{uv}{ dm_{uv}*h_pulay(2)_{uv} }  
        ! -sum_{uv}{ edm_{uv}*s(2)_{uv} }
        2.0d0*density_matrix_compute(i_compute,j_compute)* &
        temp_second_order_H_pulay(i_coord,i_center,j_coord,j_center)- &
        2.0d0*energy_density_matrix_compute(i_compute,j_compute)* &
        temp_second_order_S(i_coord,i_center,j_coord,j_center)
    enddo  ! i_center
    enddo  ! i_coord

   enddo ! j_center 
   enddo ! j_coord



  enddo ! j_compute 
  enddo ! i_compute


 if (allocated( temp_first_order_S   )) deallocate( temp_first_order_S   )
 if (allocated( temp_first_order_H_pulay   )) deallocate( temp_first_order_H_pulay   )
 if (allocated( temp_second_order_S   )) deallocate( temp_second_order_S   )
 if (allocated( temp_second_order_H_pulay   )) deallocate( temp_second_order_H_pulay   )


end subroutine evaluate_pulay_hessian_p1
!******
