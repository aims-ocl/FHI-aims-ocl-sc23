!****s* FHI-aims/evaluate_second_order_H_pulay_p0
!  NAME
!    evaluate_second_order_H_pulay_p0
!  SYNOPSIS

subroutine evaluate_second_order_H_pulay_p0(second_order_H_pulay,  & 
           n_points,partition_tab,      &
           n_compute_c, i_basis_index, index_hessian,  & 
           wave, gradient_basis_wave, H_times_psi,              &
           H_times_gradient_psi,hessian_basis_wave,             &
           first_order_rho,first_order_potential,dVxc_drho) 

!  PURPOSE
!  for phonon_gamma 
!  shanghui, 2013.12.30

!  USES

  use dimensions
  use species_data ! species_z
  use geometry ! species
  use basis  !basis_atom()
  use pbc_lists

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  integer , intent(in) :: n_compute_c
  integer , intent(in) :: i_basis_index(n_compute_c)
  integer , intent(in) :: index_hessian(3,3)

  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: gradient_basis_wave
  real*8, dimension(n_max_compute_ham,n_points), intent(in) :: H_times_psi
  real*8, dimension(n_max_compute_ham,6,n_points), intent(in) :: hessian_basis_wave
  real*8, dimension(n_max_compute_ham,3,n_points), intent(in) :: H_times_gradient_psi  

  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in) :: dVxc_drho


  real*8, dimension(3, n_atoms,3,n_atoms, n_Cbasis,n_Cbasis), intent(inout) :: second_order_H_pulay
!  INPUTS
!  o  n_points -- number of grid points in this grid batch
!  o  partition_tab -- values of partition function
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  index_hessian -- (i_coord,j_coord) to a number
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!  o  H_times_psi           -- H*psi  
!  o  H_times_gradient_psi  -- H*grad(psi)
!  o  hessian_basis_wave    -- hessian of wave function in unit cell 
!  o  first_order_rho       -- rho_tot(1)
!  o  first_order_potential -- V_free_hartree(1)+delta_V_hartree(1)
!  o  dVxc_drho             -- d Vxc/ drho

!  OUTPUT
!  o second_order_H_pulay


  integer :: i_point, i_coord,j_coord,i_basis, j_basis,i_compute,j_compute, k_atom, & 
             k_coord,i_atom, k_center, current_center

  real*8  point_term



       ! basis_atom  change to  Cbasis_to_atom

   do i_point = 1, n_points, 1

    do i_compute=1,n_compute_c,1
    do j_compute=1,n_compute_c,1
      
       i_basis=i_basis_index(i_compute)
       j_basis=i_basis_index(j_compute)


!---------------(1) begin <dXu/dRI| H |Xv> u=I_atom, v=J_atom 
        do i_coord = 1, 3, 1
           !----point_term = dXu/dRI * Xv---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(i_compute,i_coord,i_point) * wave(j_compute,i_point)
           
         do k_center = 1,  n_centers_hartree_potential, 1
            current_center   = centers_hartree_potential(k_center)
            k_atom = center_to_atom(current_center)

           !--->(1.1)-fixed : delta(I,K)*<d^2Xu/dRI dRK| H |Xv> 
           if(current_center.eq.Cbasis_to_center(i_basis)) then 
  
              do k_coord=1,3 
          second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          hessian_basis_wave(i_compute,index_hessian(i_coord,k_coord),i_point)* & 
          H_times_psi(j_compute,i_point)   
              enddo  
 
           endif
 
           !--->(1.3)-fixed : delta(J,K)*<dXu/dRI | H |dXv/dRk> 
           if(current_center.eq.Cbasis_to_center(j_basis)) then 
 
              do k_coord=1,3 
          second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, & 
                          i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, &                                     
                          i_basis,j_basis) +    &
          partition_tab(i_point)*                                   &
          gradient_basis_wave(i_compute,i_coord,i_point) *          &
          H_times_gradient_psi(j_compute,k_coord,i_point) 
              enddo
 
           endif 
         
 
         enddo ! k_center


         do k_atom = 1, n_atoms
           !--->(1.2)-fixed : <dXu/dRI| dH/dRk |Xv>
           do k_coord = 1,3          
 
           second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) =    &
           second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) +    &
            point_term*first_order_potential(k_coord, k_atom, i_point)   +    & 
            point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_atom, i_point) 
        

           enddo 
         enddo  ! k_atom
 
          
        enddo  ! i_coord
 


!---------------(2) begin <Xu| H |dXv/dRJ> u=I_atom, v=J_atom 
        do j_coord = 1, 3, 1
           !----point_term = Xu * dXv/dRJ---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(j_compute,j_coord,i_point) * wave(i_compute,i_point)
 
 

         do k_center = 1,  n_centers_hartree_potential, 1
            current_center   = centers_hartree_potential(k_center)
            k_atom = center_to_atom(current_center)
 
           !--->(2.1)-fixed : delta(I,K)*<dXu/dRI| H |dXv/dRJ> 
           if(current_center.eq.Cbasis_to_center(i_basis)) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          gradient_basis_wave(j_compute,j_coord,i_point) *          &
          H_times_gradient_psi(i_compute,k_coord,i_point) 
              enddo  
 
           endif
 
           !--->(2.3)-fixed : delta(J,K)*<Xu| H |d^2Xv/dRJ dRK> 
           if(current_center.eq.Cbasis_to_center(j_basis) ) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          hessian_basis_wave(j_compute,index_hessian(j_coord,k_coord),i_point)* & 
          H_times_psi(i_compute,i_point)   
              enddo  
 
           endif
 
         enddo ! k_center


         do k_atom=1,n_atoms 


           !--->(2.2)-fixed : <Xu| dH/dRk |dXv/dRJ> 
           do k_coord=1,3          
  
            second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                           i_basis,j_basis) =    &
            second_order_H_pulay(j_coord,Cbasis_to_atom(j_basis),k_coord,k_atom, &
                           i_basis,j_basis) +    &
            point_term*first_order_potential(k_coord, k_atom, i_point)   +    & 
            point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_atom, i_point) 

           enddo ! k_coord 

         enddo ! k_atom

        enddo  ! j_coord
  

    end do !j_compute
    end do !i_comput

  enddo ! n_points

end subroutine evaluate_second_order_H_pulay_p0
!******
