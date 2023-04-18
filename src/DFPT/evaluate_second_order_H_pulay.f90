!****s* FHI-aims/evaluate_second_order_H_pulay
!  NAME
!    evaluate_second_order_H_pulay
!  SYNOPSIS

subroutine evaluate_second_order_H_pulay(second_order_H_pulay,  & 
           partition_tab,dist_tab_global, dir_tab_global,       &
           n_points,n_compute_c, i_basis_index, index_hessian,  & 
           wave, gradient_basis_wave, H_times_psi,              &
           H_times_gradient_psi,hessian_basis_wave,             &
           first_order_rho,first_order_potential,dVxc_drho,      & 
           atom_point, partition_deriv_local, &
           first_order_rho_moving_grid,first_order_potential_moving_grid)  

!  PURPOSE
!    calculate the second H_pulay elements.
!  shanghui,2013.02.26

!   add moving grid effect 
!  shanghui,2013.06.26

!  USES

  use dimensions
  use species_data ! species_z
  use geometry ! species
  !use runtime_choices
  use basis  !basis_atom()

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) ::  partition_tab
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_global
  real*8, dimension(3,n_atoms, n_points), intent(in) :: dir_tab_global
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


  integer, dimension(n_max_batch_size), intent(in) ::  atom_point
  real*8, dimension(3,n_atoms,n_points), intent(in) :: partition_deriv_local
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho_moving_grid
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_potential_moving_grid



  real*8, dimension(3, n_atoms,3,n_atoms, n_basis,n_basis), intent(inout) :: second_order_H_pulay
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch



!  OUTPUT
!  o second_order_H_pulay


  integer :: i_point, i_coord,j_coord,i_basis, j_basis,i_compute,j_compute, k_atom, & 
             k_coord,i_atom

  real*8  point_term, atomic_term
  real*8  atomic_term_moving_grid(3, n_atoms, n_points )



   do i_point = 1, n_points, 1

         !write(use_unit,*) i_point
         !write(use_unit,*) partition_tab(i_point),partition_deriv_local(1,1,i_point) 
         !write(use_unit,*) H_times_psi(1,i_point),H_times_gradient_psi(1,1,i_point)
         !write(use_unit,*) wave(1,i_point),gradient_basis_wave(1,1,i_point)
         !write(use_unit,*) Hessian_basis_wave(1,index_hessian(1,1),i_point)

     !-------------shanghui begin atomic_term_moving_grid------      
      do k_atom = 1, n_atoms 
      atomic_term_moving_grid(1:3, k_atom, i_point)=0.d0 

         do i_atom = 1, n_atoms 
   
          if(i_atom.ne.k_atom) then 
           do i_coord = 1, 3 
  
          atomic_term_moving_grid(i_coord, k_atom, i_point)= & 
          atomic_term_moving_grid(i_coord, k_atom, i_point)+ & 
            species_z(species(i_atom))  &
          * dir_tab_global(i_coord, i_atom, i_point) &
          / (dist_tab_global(i_atom, i_point))**3.0d0

           enddo ! i_coord
          endif  ! i.ne.k 
        
         enddo ! i_atom 
      enddo   ! k_atom 
     !-------------shanghui end atomic_term_moving_grid------      


    do i_compute=1,n_compute_c,1
    do j_compute=1,n_compute_c,1
      
       i_basis=i_basis_index(i_compute)
       j_basis=i_basis_index(j_compute)

!---------------(0) begin partition_deriv-------------------------
        do i_coord = 1, 3, 1
         do k_atom = 1,n_atoms 
         do k_coord = 1,3          
 
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                        i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                        i_basis,j_basis) +    &
          !------patrition_deriv(k_atom)---
          partition_deriv_local(k_coord,k_atom,i_point)* &
          (-gradient_basis_wave(i_compute,i_coord,i_point))* &
          H_times_psi(j_compute,i_point) 
 
         second_order_H_pulay(i_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
         second_order_H_pulay(i_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
         !------patrition_deriv(k_atom)---
         partition_deriv_local(k_coord,k_atom,i_point)* & 
         (-gradient_basis_wave(j_compute,i_coord,i_point))* & 
         H_times_psi(i_compute,i_point) 
          
          enddo 
          enddo 
        enddo  



!---------------(1) begin <dXu/dRI| H |Xv> u=I_atom, v=J_atom 
        do i_coord = 1, 3, 1
           !----point_term = dXu/dRI * Xv---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(i_compute,i_coord,i_point) * wave(j_compute,i_point)
           
         do k_atom = 1,n_atoms 
 
         if(atom_point(i_point).ne.k_atom ) then !fixed
 
           !--->(1.1)-fixed : delta(I,K)*<d^2Xu/dRI dRK| H |Xv> 
           if(k_atom.eq.basis_atom(i_basis)) then 
  
              do k_coord=1,3 
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          hessian_basis_wave(i_compute,index_hessian(i_coord,k_coord),i_point)* & 
          H_times_psi(j_compute,i_point)   
              enddo  
 
           endif
 
           !--->(1.2)-fixed : <dXu/dRI| dH/dRk |Xv>
           do k_coord = 1,3          
            atomic_term = -species_z(species(k_atom))  & 
                          * dir_tab_global(k_coord, k_atom, i_point) & 
                          / (dist_tab_global(k_atom, i_point))**3.0d0
 
           second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) =    &
           second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) +    &
            point_term*atomic_term   +    & 
            point_term*first_order_potential(k_coord, k_atom, i_point)   +    & 
            point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_atom, i_point) 
           enddo 
 
          
           !--->(1.3)-fixed : delta(J,K)*<dXu/dRI | H |dXv/dRk> 
           if(k_atom.eq.basis_atom(j_basis)) then  
 
              do k_coord=1,3 
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, & 
                          i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
          i_basis,j_basis) +    &
          partition_tab(i_point)*                                   &
          gradient_basis_wave(i_compute,i_coord,i_point) *          &
          H_times_gradient_psi(j_compute,k_coord,i_point) 
              enddo
 
           endif 
 

 
         
         else !moving 
 
           !--->(1.1)-moving : nodelta(I,K)*<d^2Xu/dRI dRK| H |Xv> 
           if(k_atom.ne.basis_atom(i_basis)) then 
  
              do k_coord=1,3 
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          (-hessian_basis_wave(i_compute,index_hessian(i_coord,k_coord),i_point))* & 
          H_times_psi(j_compute,i_point)   
              enddo  
 
           endif

           !--->(1.2)-moving : <dXu/dRI| dH/dRk |Xv>
           do k_coord = 1,3          
           second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) =    &
           second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
                          i_basis,j_basis) +    &
            point_term*atomic_term_moving_grid(k_coord, k_atom, i_point)  +    & 
            point_term*first_order_potential_moving_grid(k_coord, k_atom, i_point) + & 
            point_term*dVxc_drho(i_point)*first_order_rho_moving_grid(k_coord,k_atom, i_point) 
           enddo ! k_coord 
 
 
 
           !--->(1.3)-moving : nodelta(J,K)*<dXu/dRI | H |dXv/dRk> 
           if(k_atom.ne.basis_atom(j_basis))  then 
 
              do k_coord=1,3 
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, & 
                          i_basis,j_basis) =    &
          second_order_H_pulay(i_coord,basis_atom(i_basis),k_coord,k_atom, &
          i_basis,j_basis) +    &
          partition_tab(i_point)*                                   &
          gradient_basis_wave(i_compute,i_coord,i_point) *          &
          (-H_times_gradient_psi(j_compute,k_coord,i_point)) 
              enddo
 
           endif 
 
         endif ! fix or move 
 
         enddo ! k_atom 
        enddo  ! i_coord
 


!---------------(2) begin <Xu| H |dXv/dRJ> u=I_atom, v=J_atom 
        do j_coord = 1, 3, 1
           !----point_term = Xu * dXv/dRJ---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(j_compute,j_coord,i_point) * wave(i_compute,i_point)
 
         do k_atom=1,n_atoms 
 
         if(atom_point(i_point).ne.k_atom ) then !fixed
 
           !--->(2.1)-fixed : delta(I,K)*<dXu/dRI| H |dXv/dRJ> 
           if(k_atom.eq.basis_atom(i_basis)) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          gradient_basis_wave(j_compute,j_coord,i_point) *          &
          H_times_gradient_psi(i_compute,k_coord,i_point) 
              enddo  
 
           endif
 
           !--->(2.2)-fixed : <Xu| dH/dRk |dXv/dRJ> 
           do k_coord=1,3          
             atomic_term = -species_z(species(k_atom))     &
                     * dir_tab_global(k_coord, k_atom, i_point) & 
                     / (dist_tab_global(k_atom, i_point))**3.0d0  
  
            second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                           i_basis,j_basis) =    &
            second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                           i_basis,j_basis) +    &
            point_term*atomic_term    +    & 
            point_term*first_order_potential(k_coord, k_atom, i_point)   +    & 
            point_term*dVxc_drho(i_point)*first_order_rho(k_coord,k_atom, i_point) 
           enddo ! k_coord 
  

           !--->(2.3)-fixed : delta(J,K)*<Xu| H |d^2Xv/dRJ dRK> 
           if(k_atom.eq.basis_atom(j_basis)) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          hessian_basis_wave(j_compute,index_hessian(j_coord,k_coord),i_point)* & 
          H_times_psi(i_compute,i_point)   
              enddo  
 
           endif
 
 
         else  !moving 
 
           !--->(2.1)-moving : nodelta(I,K)*<dXu/dRI| H |dXv/dRJ> 
           if(k_atom.ne.basis_atom(i_basis)) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          gradient_basis_wave(j_compute,j_coord,i_point) *          &
          (-H_times_gradient_psi(i_compute,k_coord,i_point)) 
              enddo  
 
           endif
 
          !--->(2.2)-moving : <Xu| dH/dRk |dXv/dRJ> 
         do k_coord = 1,3          
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
           point_term*atomic_term_moving_grid(k_coord, k_atom, i_point)  +    & 
           point_term*first_order_potential_moving_grid(k_coord, k_atom, i_point) + & 
           point_term*dVxc_drho(i_point)*first_order_rho_moving_grid(k_coord,k_atom, i_point) 
         enddo ! k_coord 

         !--->(2.3)-moving : nodelta(J,K)*<Xu| H |d^2Xv/dRJ dRK> 
           if(k_atom.ne.basis_atom(j_basis)) then
  
              do k_coord=1,3 
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) =    &
          second_order_H_pulay(j_coord,basis_atom(j_basis),k_coord,k_atom, &
                         i_basis,j_basis) +    &
          partition_tab(i_point)*  & 
          (-hessian_basis_wave(j_compute,index_hessian(j_coord,k_coord),i_point))* & 
          H_times_psi(i_compute,i_point)   
              enddo  
 
           endif

         endif ! fix or move 
 
 
         enddo ! k_atom
        enddo  ! j_coord
  

    end do !j_compute
    end do !i_comput

  enddo ! n_points


end subroutine evaluate_second_order_H_pulay
!******
