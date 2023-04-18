!****s* FHI-aims/evaluate_pulay_hessian_fixed
!  NAME
!    evaluate_pulay_hessian_fixed
!  SYNOPSIS

subroutine evaluate_pulay_hessian_fixed_reduce_memory(pulay_hessian,  & 
           partition_tab,dist_tab_global, dir_tab_global,       &
           n_points,n_compute_c, i_basis_index, index_hessian,  & 
           wave, gradient_basis_wave, H_times_psi,              &
           H_times_gradient_psi,hessian_basis_wave,             &
           first_order_rho,first_order_potential,dVxc_drho,     &
           density_matrix,first_order_density_matrix, & 
           energy_density_matrix, first_order_energy_density_matrix, & 
           j_atom,j_coord,kinetic_gradient_basis_wave )

!  PURPOSE
!  calculate the second H_pulay elements.
!  shanghui,2013.02.26

!  add moving grid effect 
!  shanghui,2013.06.26

!  add pulay_hessian_fixed
!  shanghui, 2015.05.05
!  USES

  use dimensions
  use species_data ! species_z
  use geometry ! species
  use runtime_choices
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

  real*8, dimension( n_points), intent(in) :: first_order_rho
  real*8, dimension( n_points), intent(in) :: first_order_potential
  real*8, dimension(n_points), intent(in)  :: dVxc_drho

  real*8, dimension(n_basis,n_basis), intent(in) :: density_matrix
  real*8, dimension(n_basis,n_basis), intent(in) :: first_order_density_matrix
  real*8, dimension(n_basis,n_basis), intent(in) :: energy_density_matrix
  real*8, dimension(n_basis,n_basis), intent(in) :: first_order_energy_density_matrix
  integer, intent(in) :: j_atom,j_coord
!!!DSB atomic zora
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: kinetic_gradient_basis_wave
  real*8, dimension(3, n_atoms,3,n_atoms), intent(inout) :: pulay_hessian
!  INPUTS
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch



!  OUTPUT
!  o pulay_hessian


  integer :: i_point,i_atom, i_coord,i_basis, j_basis,i_compute,j_compute

  real*8 ::  point_term, atomic_term

  real*8 ::  temp_first_order_S(3,n_atoms),temp_first_order_H_pulay(3,n_atoms), & 
             temp_second_order_S(3,n_atoms),temp_second_order_H_pulay(3,n_atoms)

  integer :: i,j

 !------calculate--------

  do i_compute=1,n_compute_c,1
  do j_compute=1,n_compute_c,1

       i_basis=i_basis_index(i_compute)
       j_basis=i_basis_index(j_compute)

       temp_first_order_S(1:3,1:n_atoms) = 0.0d0
       temp_first_order_H_pulay(1:3,1:n_atoms) = 0.0d0

       temp_second_order_S(1:3,1:n_atoms) = 0.0d0
       temp_second_order_H_pulay(1:3,1:n_atoms) = 0.0d0

   do i_point = 1, n_points, 1
!---------------1. first_order_H_pulay---------------------------
       do i_coord = 1, 3, 1
          temp_first_order_H_pulay(i_coord,basis_atom(i_basis)) =    &
          temp_first_order_H_pulay(i_coord,basis_atom(i_basis)) -    &
            partition_tab(i_point)* &  
            gradient_basis_wave(i_compute,i_coord,i_point)*H_times_psi(j_compute,i_point) 

          temp_first_order_H_pulay(i_coord,basis_atom(j_basis)) =    &
          temp_first_order_H_pulay(i_coord,basis_atom(j_basis)) -    &
            partition_tab(i_point)* &  
            gradient_basis_wave(j_compute,i_coord,i_point)*H_times_psi(i_compute,i_point) 
          if(flag_rel .eq. REL_atomic_zora)then   !!!DSB
            temp_first_order_H_pulay(i_coord,basis_atom(i_basis)) =    &
            temp_first_order_H_pulay(i_coord,basis_atom(i_basis)) - &
              partition_tab(i_point)* &
              kinetic_gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)

            temp_first_order_H_pulay(i_coord,basis_atom(j_basis)) =    &
            temp_first_order_H_pulay(i_coord,basis_atom(j_basis)) - &
              partition_tab(i_point)* &
              kinetic_gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point)

          endif 
       enddo 
!---------------2. first_order_S---------------------------------
       if(i_basis.ne.j_basis) then
       do i_coord = 1, 3, 1
          temp_first_order_S(i_coord,basis_atom(i_basis)) =    &
          temp_first_order_S(i_coord,basis_atom(i_basis)) -    &
            partition_tab(i_point)* &  
            gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point) 

          temp_first_order_S(i_coord,basis_atom(j_basis)) =    &
          temp_first_order_S(i_coord,basis_atom(j_basis)) -    &
            partition_tab(i_point)* &  
            gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point) 
       enddo
       endif 




!---------------3. second_order_H_pulay--------------------------
        do i_coord = 1, 3, 1
           !----point_term = dXu/dRI * Xv---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(i_compute,i_coord,i_point) * wave(j_compute,i_point)
!---------------(1) begin <dXu/dRI| H |Xv> u=I_atom, v=J_atom 
            !--->(1.1)-fixed : delta(I,J)*<d^2Xu/dRI dRJ| H |Xv>
           if(j_atom.eq.basis_atom(i_basis)) then  
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) = & 
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) + & 
              partition_tab(i_point)*  &
              hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point)* &
              H_times_psi(j_compute,i_point)
           endif 

           !--->(1.2)-fixed : <dXu/dRI| dH/dRk |Xv>
            atomic_term = -species_z(species(j_atom))  & 
                          * dir_tab_global(j_coord, j_atom, i_point) & 
                          / (dist_tab_global(j_atom, i_point))**3.0d0
 
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) = &
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) + &
              point_term*atomic_term   +    & 
              point_term*first_order_potential(i_point)   +    & 
              point_term*dVxc_drho(i_point)*first_order_rho(i_point) 
          
           !--->(1.3)-fixed : delta(J,K)*<dXu/dRI | H |dXv/dRk> 
           if(j_atom.eq.basis_atom(j_basis)) then  
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) = &
            temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) + &
              partition_tab(i_point)*                                   &
              gradient_basis_wave(i_compute,i_coord,i_point) *          &
              H_times_gradient_psi(j_compute,j_coord,i_point) 

             if(flag_rel .eq. REL_atomic_zora)then   !!!DSB
               temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) = &
               temp_second_order_H_pulay(i_coord,basis_atom(i_basis)) + &
               partition_tab(i_point)*                                  &
                gradient_basis_wave(j_compute,i_coord,i_point) *          &
                kinetic_gradient_basis_wave(i_compute,j_coord,i_point)
             endif      

           endif 

!---------------(2) begin <Xu| H |dXv/dRI>--------------------------- 
           !----point_term = Xu * dXv/dRI---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(j_compute,i_coord,i_point) * wave(i_compute,i_point)
 
           !--->(2.1)-fixed : delta(i_basis,J)*<dXu/dRJ| H |dXv/dRI> 
           if(j_atom.eq.basis_atom(i_basis)) then
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) = & 
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) + &
              partition_tab(i_point)*  & 
              gradient_basis_wave(j_compute,i_coord,i_point) *          & ! |dXv/dRI> 
              H_times_gradient_psi(i_compute,j_coord,i_point)  ! <dXu/dRJ| H 
             if(flag_rel .eq. REL_atomic_zora)then   !!!DSB
              temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) = &
              temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) + &
                partition_tab(i_point)*  &
                gradient_basis_wave(i_compute,i_coord,i_point) *          & ! |dXv/dRI> 
                kinetic_gradient_basis_wave(j_compute,j_coord,i_point)  ! <dXu/dRJ| H 
             endif
           endif
 
           !--->(2.2)-fixed : <Xu| dH/dRJ |dXv/dRI> 
            atomic_term = -species_z(species(j_atom))     &
                     * dir_tab_global(j_coord, j_atom, i_point) & 
                     / (dist_tab_global(j_atom, i_point))**3.0d0  
  
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) = & 
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) + &
              point_term*atomic_term    +    & 
              point_term*first_order_potential(i_point)   +    & 
              point_term*dVxc_drho(i_point)*first_order_rho( i_point) 

           !--->(2.3)-fixed : delta(J,K)*<Xu| H |d^2Xv/dRJ dRK> 
           if(j_atom.eq.basis_atom(j_basis)) then
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) = & 
            temp_second_order_H_pulay(i_coord,basis_atom(j_basis)) + &
              partition_tab(i_point)*  & 
              hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point)* & 
              H_times_psi(i_compute,i_point)   
           endif
        enddo  ! i_coord
  

!---------------4.second_order_S----------------------------------------
        if(i_basis.ne.j_basis) then
        do i_coord = 1, 3, 1
           !----point_term = dXu/dRI * Xv---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(i_compute,i_coord,i_point) * wave(j_compute,i_point)
!---------------(1) begin <dXu/dRI|Xv> u=I_atom, v=J_atom 
            !--->(1.1)-fixed : delta(I,J)*<d^2Xu/dRI dRJ| Xv>
           if(j_atom.eq.basis_atom(i_basis)) then  
            temp_second_order_S(i_coord,basis_atom(i_basis)) = & 
            temp_second_order_S(i_coord,basis_atom(i_basis)) + & 
              partition_tab(i_point)*  &
              hessian_basis_wave(i_compute,index_hessian(i_coord,j_coord),i_point)* &
              wave(j_compute,i_point)
           endif 
          
           !--->(1.2)-fixed : delta(J,K)*<dXu/dRI  |dXv/dRk> 
           if(j_atom.eq.basis_atom(j_basis)) then  
            temp_second_order_S(i_coord,basis_atom(i_basis)) = &
            temp_second_order_S(i_coord,basis_atom(i_basis)) + &
              partition_tab(i_point)*                                   &
              gradient_basis_wave(i_compute,i_coord,i_point) *          &
              gradient_basis_wave(j_compute,j_coord,i_point) 
           endif 

!---------------(2) begin <Xu |dXv/dRI>--------------------------- 
           !----point_term = Xu * dXv/dRI---
           point_term=-partition_tab(i_point) * &
           gradient_basis_wave(j_compute,i_coord,i_point) * wave(i_compute,i_point)
 
           !--->(2.1)-fixed : delta(i_basis,J)*<dXu/dRJ |dXv/dRI> 
           if(j_atom.eq.basis_atom(i_basis)) then
            temp_second_order_S(i_coord,basis_atom(j_basis)) = & 
            temp_second_order_S(i_coord,basis_atom(j_basis)) + &
              partition_tab(i_point)*  & 
              gradient_basis_wave(j_compute,i_coord,i_point) *          & ! |dXv/dRI> 
              gradient_basis_wave(i_compute,j_coord,i_point)  ! <dXu/dRJ| 
           endif

           !--->(2.2)-fixed : delta(J,K)*<Xu |d^2Xv/dRJ dRK> 
           if(j_atom.eq.basis_atom(j_basis)) then
            temp_second_order_S(i_coord,basis_atom(j_basis)) = & 
            temp_second_order_S(i_coord,basis_atom(j_basis)) + &
              partition_tab(i_point)*  & 
              hessian_basis_wave(j_compute,index_hessian(i_coord,j_coord),i_point)* & 
              wave(i_compute,i_point)   
           endif
        enddo  ! i_coord
        endif

    enddo ! n_points
       
!       write(use_unit,*) 'i_basis,j_basis:',i_basis,j_basis
!       write(use_unit,*) '----------------------H_pulay---------------------------------------------------------'
!       write(use_unit,*) '************shanghui begain first_order_H_pulay(x 1)****************'
!        write(use_unit,'(f15.9)') temp_first_order_H_pulay(1,1 )
!       write(use_unit,*) '************shanghui end first_order_H_pulay****************'
!       write(use_unit,*) '************shanghui begain first_order_H_pulay(x 2)****************'
!        write(use_unit,'(f15.9)') temp_first_order_H_pulay(1,2 )
!       write(use_unit,*) '************shanghui end first_order_H_pulay****************'
!       write(use_unit,*) ''
!       write(use_unit,*) '************shanghui begain second_order_H_pulay(x 1 x 1)****************'
!        write(use_unit,'(f15.9)') temp_second_order_H_pulay(1,1)
!       write(use_unit,*) '************shanghui end second_order_H_pulay****************'
! 
!       write(use_unit,*) ''
!       write(use_unit,*) '----------------------S---------------------------------------------------------'
!       write(use_unit,*) '************shanghui begain first_order_S(x 1)****************'
!        write(use_unit,'(f15.9)') temp_first_order_S(1,1 )
!       write(use_unit,*) '************shanghui end first_order_S****************'
!       write(use_unit,*) ''
!       write(use_unit,*) '************shanghui begain first_order_S(x 2)****************'
!        write(use_unit,'(f15.9)') temp_first_order_S(1,2 )
!       write(use_unit,*) '************shanghui end first_order_S****************'
!       write(use_unit,*) '************shanghui begain second_order_S(x 1 x 1)****************'
!        write(use_unit,'(6f15.9)') temp_second_order_S(1,1)
!       write(use_unit,*) '************shanghui end second_order_S****************'


!--------------5. sum up to pulay_hessian------------------------------

    do i_coord =1, 3 
    do i_atom = 1, n_atoms 
       pulay_hessian(i_coord,i_atom,j_coord,j_atom)= &
       pulay_hessian(i_coord,i_atom,j_coord,j_atom)+ &
        !-------so we have pulay_hessian--->(1) 
        ! sum_{uv}{ dm(1)_{uv}*h_pulay(1)_{uv} }  
        ! -sum_{uv}{ edm(1)_{uv}*s(1)_{uv} }
        first_order_density_matrix(i_basis,j_basis)*  &
        temp_first_order_H_pulay(i_coord,i_atom)-   &
        first_order_energy_density_matrix(i_basis,j_basis)*  &
        temp_first_order_S(i_coord,i_atom)+  &
        !-------and the pulay_hessian   --->(2) 
        ! sum_{uv}{ dm_{uv}*h_pulay(2)_{uv} }  
        ! -sum_{uv}{ edm_{uv}*s(2)_{uv} }
        density_matrix(i_basis,j_basis)* &
        temp_second_order_H_pulay(i_coord,i_atom)- &
        energy_density_matrix(i_basis,j_basis)* &
        temp_second_order_S(i_coord,i_atom)
    enddo 
    enddo

  end do !j_compute
  end do !i_comput


end subroutine evaluate_pulay_hessian_fixed_reduce_memory
!******
