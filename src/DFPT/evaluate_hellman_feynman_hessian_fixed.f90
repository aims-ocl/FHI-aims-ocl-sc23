!****s* FHI-aims/evaluate_hellman_feynman_hessian
!  NAME
!    evaluate_hellman_feynman_hessian
!  SYNOPSIS

subroutine evaluate_hellman_feynman_hessian_fixed(hellman_feynman_hessian,        &
           n_points,partition_tab, dist_tab_global, dir_tab_global,       &
           rho,first_order_rho)

!  PURPOSE
!    calculate Hellman_Feyman_hessian
!    
!   (1) Int{ first_order_rho(alpha) * Z_I *(R_I_beta-r)/|r-R_I|^3  }  
!   (2) Int{ rho * delta(alpha,beta) * (Z_I/|r-R_I|^3 - 3*Z_I*(R_alpha-r_alpha)^2/|r-R_I|^5  ) 

!  shanghui,2013.02.21 @ Berlin 

!  USES

  use dimensions
  use species_data
  use geometry
  use localorb_io, only: use_unit
!  use runtime_choices
!  use basis 

!  ARGUMENTS

  implicit none
  integer :: n_points
  real*8, dimension(n_points), intent(in) :: partition_tab
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_global
  real*8, dimension(3,n_atoms, n_points), intent(in) :: dir_tab_global

  real*8, dimension(n_points), intent(in) :: rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho

  real*8, dimension(3, n_atoms, 3, n_atoms), intent(out) :: hellman_feynman_hessian

!  INPUTS
!  o  rho -- electron density
!  o  dir_tab_global -- direction to atoms
!  o  dist_tab_global -- (distance to atoms)**1
!  o  partition_tab -- values of partition function
!  o  n_points -- number of grid points in this grid batch
!
!  OUTPUT
!  hellman_feyman_hessian


  integer :: i_point, i_atom, i_coord, j_atom, j_coord

  real*8 :: atomic_term(n_atoms)


 do i_point = 1, n_points, 1

  do i_atom=1,n_atoms
        atomic_term(i_atom) =  -species_z(species(i_atom))  &
                              / (dist_tab_global(i_atom, i_point))**3.0d0
  do i_coord=1,3   
     
   do j_atom=1,n_atoms
   do j_coord=1,3 


     enddo !j_coord



       if(j_atom.ne.i_atom) then ! fix R1,R2

        do j_coord=1,3 
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho(j_coord, j_atom,i_point) !&
        enddo    !j_coord
       
      else if(j_atom.eq.i_atom) then! fix R1,R1 
        do j_coord=1,3 
           if(j_coord.eq.i_coord) then
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho(j_coord, j_atom,i_point) &
       +partition_tab(i_point) * rho(i_point) *                     &
       (species_z(species(i_atom))/(dist_tab_global(i_atom, i_point))**3.0d0 -      & 
       3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  & 
    dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
          else 
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho(j_coord, j_atom,i_point) &
        +partition_tab(i_point) * rho(i_point) *                     &
       (-3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  & 
    dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
          endif
 
        enddo ! j_coord



      else 
        write(use_unit,*) 'error for j_atom in evaluate_hellman_feynman_hessian'
        stop 
      endif 

    enddo !j_atom

  enddo ! i_coord
  enddo ! i_atom

  
 end do ! n_points
end subroutine evaluate_hellman_feynman_hessian_fixed
!******
