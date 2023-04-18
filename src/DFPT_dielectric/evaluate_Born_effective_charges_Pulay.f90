!****s* FHI-aims/evaluate_Born_effective_charges_Pulay
!  NAME
!    evaluate_Born_effective_charges_Pulay
!  SYNOPSIS

subroutine evaluate_Born_effective_charges_Pulay & 
         ( partition_tab, & 
           n_points, n_compute_c, i_basis_index,  & 
           wave, gradient_basis_wave, & 
           kinetic_gradient_basis_wave, & 
           H_times_psi, & 
           first_order_rho,first_order_potential, &
           dVxc_drho, &
           zero_order_density_matrix_compute, &
           first_order_density_matrix_compute, &
           first_order_energy_density_matrix_compute, &
           Born_effective_charges_Pulay, &
           test_Born_effective_charges_by_print_force )

!  PURPOSE
!
!   shanghui 2017.01.24, for Born effective charges
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
  real*8, dimension(n_max_compute_ham, 3, n_points),intent(in) :: kinetic_gradient_basis_wave 
  real*8, dimension(n_max_compute_ham, n_points), intent(in) :: H_times_psi

  real*8, dimension(n_points), intent(in) :: first_order_rho
  real*8, dimension(n_points), intent(in) :: first_order_potential
  !-------LDA---------------
  real*8, dimension(3,n_points), intent(in) :: dVxc_drho

  real*8, dimension(n_compute_c,n_compute_c), intent(in) :: zero_order_density_matrix_compute
  real*8, dimension(n_compute_c,n_compute_c), intent(in) :: first_order_density_matrix_compute
  real*8, dimension(n_compute_c,n_compute_c), intent(in) :: first_order_energy_density_matrix_compute

  real*8, dimension(3,n_atoms), intent(inout) :: Born_effective_charges_Pulay
  logical, intent(in) :: test_Born_effective_charges_by_print_force

!  INPUTS
!  o  n_points -- number of grid points in this grid batch
!  o  partition_tab -- values of partition function
!  o  n_compute_c -- basis set computed at this subroutine
!  o  i_basis_index -- i_compute to i_basis
!  o  wave -- wave function in unit cell
!  o  gradient_basis_wave -- gradient wave function in unit cell
!  o  H_times_psi           -- H*psi  
!
!  OUTPUT
!  Born_effect_charge_Pulay 


  integer :: i_compute,j_compute, i_basis,j_basis, i_basis_uc, j_basis_uc
  integer :: i_cell, j_cell
  integer :: i_point, i_place, i_coord,i_atom

  real*8 ::  temp_first_order_S(3,n_atoms)
  real*8 ::  temp_first_order_H_pulay(3,n_atoms), temp_second_order_H_pulay(3,n_atoms)


  do i_compute=1,n_compute_c,1
     i_basis    = i_basis_index(i_compute)
     i_basis_uc = Cbasis_to_basis(i_basis)
     i_cell     = center_to_cell(Cbasis_to_center(i_basis))

  do j_compute=1,n_compute_c,1
     j_basis    = i_basis_index(j_compute)
     j_basis_uc = Cbasis_to_basis(j_basis)
     j_cell     = center_to_cell(Cbasis_to_center(j_basis))

       temp_first_order_H_pulay(1:3,1:n_atoms) = 0.0d0
       temp_first_order_S(1:3,1:n_atoms) = 0.0d0
       temp_second_order_H_pulay(1:3,1:n_atoms) = 0.0d0

       do i_point = 1, n_points, 1

!---------------1. first_order_H_pulay---------------------------
           do i_coord = 1 ,3
            !delta(u',R3) * \int ( V(r)  d Xu'/dR3  Xv' ) 
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) = &
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) - &
            partition_tab(i_point)* &                                               
            gradient_basis_wave(i_compute,i_coord,i_point)*H_times_psi(j_compute,i_point)

            !delta(v',R3) * \int ( V(r) Xu' d Xv'/dR3 ) 
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) = &
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) - &
            partition_tab(i_point)* &                               
            gradient_basis_wave(j_compute,i_coord,i_point)*H_times_psi(i_compute,i_point) 

           enddo !i_coord


!---------------2. second_order_H_pulay---------------------------
           do i_coord = 1 ,3
            !delta(u',R3) * \int ( V(r)  d Xu'/dR3  Xv' ) 
            temp_second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) = &
            temp_second_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) - &
            partition_tab(i_point)* &                                               
            gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point) * &
            ( first_order_potential( i_point) + & 
              dVxc_drho(1,i_point)*first_order_rho(i_point) )

            !delta(v',R3) * \int ( V(r) Xu' d Xv'/dR3 ) 
            temp_second_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) = &
            temp_second_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) - &
            partition_tab(i_point)* &                               
            gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point)* &
            ( first_order_potential( i_point) + & 
              dVxc_drho(1,i_point)*first_order_rho(i_point) )

           enddo !i_coord



!---------------3. first_order_S = d Su'v'/ d ~u^*_I(q) ---------------------------
! shanghui note: By removing the following, I get exactly the same Pulay force as 
!       forces_densmat.f90 for Si system.
!           if(Cbasis_to_center(i_basis).ne.Cbasis_to_center(j_basis)) then ! here we use transition conservation
           do i_coord = 1 ,3
            !delta(u',R3) * \int ( d Xu'/dR3 e(iqr) e(-qR3)  Xv' ) 
            temp_first_order_S(i_coord,Cbasis_to_atom(i_basis)) = &
            temp_first_order_S(i_coord,Cbasis_to_atom(i_basis)) - &
            partition_tab(i_point)* &                                              
            gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)

            !delta(v',R3) * \int ( Xu' d Xv'/dR3 e(iqr) e(-qR3 ) 
            temp_first_order_S(i_coord,Cbasis_to_atom(j_basis)) = &
            temp_first_order_S(i_coord,Cbasis_to_atom(j_basis)) - &
            partition_tab(i_point)* &                               
            gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point) 
           enddo !i_coord
 !          endif !transition conservation
 
       enddo ! i_point


!---------------4. add atomic_zora to first_order_H_pulay---------------------------
       if(flag_rel == REL_atomic_zora)then
       do i_point = 1, n_points, 1
           do i_coord = 1 ,3
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) = &
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(i_basis)) - &
            partition_tab(i_point)* &                                               
            kinetic_gradient_basis_wave(i_compute,i_coord,i_point)*wave(j_compute,i_point)

            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) = &
            temp_first_order_H_pulay(i_coord,Cbasis_to_atom(j_basis)) - &
            partition_tab(i_point)* &                               
            kinetic_gradient_basis_wave(j_compute,i_coord,i_point)*wave(i_compute,i_point) 
           enddo !i_coord
       enddo ! i_point
       endif ! REL_atomic_zora


!--------------5. sum up to pulay_term------------------------------
   if(test_Born_effective_charges_by_print_force) then 

    do i_coord =1, 3
    do i_atom = 1, n_atoms
       Born_effective_charges_Pulay(i_coord,i_atom)= &
       Born_effective_charges_Pulay(i_coord,i_atom) &
        !-------so we have pulay_term for force--->(1) 
        ! -sum_{uv}{ dm(0)_{uv}*h_pulay(1)_{uv} }  
        ! +sum_{uv}{ edm(0)_{uv}*s(1)_{uv} }
        - first_order_density_matrix_compute(i_compute,j_compute)*  &
        temp_first_order_H_pulay(i_coord,i_atom)     &
        + first_order_energy_density_matrix_compute(i_compute,j_compute)*  &
        temp_first_order_S(i_coord,i_atom)  
    enddo
    enddo

   else 

    do i_coord =1, 3
    do i_atom = 1, n_atoms
       Born_effective_charges_Pulay(i_coord,i_atom)= &
       Born_effective_charges_Pulay(i_coord,i_atom) &
        !-------so we have pulay_term--->(1)
        ! -sum_{uv}{ dm(0)_{uv}*h_pulay(2)_{uv} } 
        ! -sum_{uv}{ dm(1)_{uv}*h_pulay(1)_{uv} }  
        ! +sum_{uv}{ edm(1)_{uv}*s(1)_{uv} }
        - zero_order_density_matrix_compute(i_compute,j_compute)*  &
        temp_second_order_H_pulay(i_coord,i_atom)     &
        - first_order_density_matrix_compute(i_compute,j_compute)*  &
        temp_first_order_H_pulay(i_coord,i_atom)     &
        + first_order_energy_density_matrix_compute(i_compute,j_compute)*  &
        temp_first_order_S(i_coord,i_atom)  
    enddo
    enddo

   endif ! test_force 

  enddo  ! j_compute
  enddo  ! i_compute 

end subroutine evaluate_Born_effective_charges_Pulay
!******
