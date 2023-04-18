!****s* FHI-aims/evaluate_hellman_feynman_hessian
!  NAME
!    evaluate_hellman_feynman_hessian
!  SYNOPSIS

subroutine evaluate_hellman_feynman_hessian(hellman_feynman_hessian,        &
           n_points,partition_tab,partition_deriv_local, dist_tab_global, dir_tab_global,       &
           rho,atom_point,first_order_rho, first_order_rho_moving_grid)

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
  real*8, dimension(3,n_atoms,n_points), intent(in) :: partition_deriv_local
  real*8, dimension(n_atoms, n_points), intent(in) :: dist_tab_global
  real*8, dimension(3,n_atoms, n_points), intent(in) :: dir_tab_global

  real*8, dimension(n_points), intent(in) :: rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho
  real*8, dimension(3, n_atoms, n_points), intent(in) :: first_order_rho_moving_grid
  integer, dimension(n_max_batch_size), intent(in) ::  atom_point

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

!!!debug!--------------------(1) begin middle point------------------
!!!debug         if(i_atom.eq.1.and.j_atom.eq.1.and. & 
!!!debug           i_coord.eq.2 .and. j_coord.eq.2 .and. & 
!!!debug        dabs(dir_tab_global(1,2,i_point)-0.6470719466891).lt.1.0d-4 .and. & 
!!!debug        dabs(dir_tab_global(2,2,i_point)).lt.1.0d-4 .and.  & 
!!!debug        dabs(dir_tab_global(3,2,i_point)).lt.1.0d-4) then 
!!!debug         
!!!debug           write(use_unit,*) '' 
!!!debug           write(use_unit,*) '=======in evaluate_hellman_hellman_hessian (middle-point)=======' 
!!!debug           write(use_unit,*) 'dir-to-atom1:',dir_tab_global(1:3, 1, i_point)
!!!debug           write(use_unit,*) 'dir-to-atom2:',dir_tab_global(1:3, 2, i_point)
!!!debug           write(use_unit,*) 'belong to atom--->',atom_point(i_point) 
!!!debug           write(use_unit,*) 'hessian:',& 
!!!debug                    ! partition_tab(i_point)  *   &       
!!!debug           dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
!!!debug           first_order_rho(j_coord, i_atom,i_point)+                      & 
!!!debug                    ! partition_tab(i_point) * & 
!!!debug                          rho(i_point)* & 
!!!debug        ( species_z(species(i_atom))/(dist_tab_global(i_atom, i_point))**3.0d0 -      &
!!!debug          3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  &
!!!debug   dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
!!!debug           write(use_unit,*) 'partition_tab:', partition_tab(i_point)
!!!debug
!!!debug         endif
!!!debug!--------------------(1) end middle point------------------
!!!debug
!!!debug!--------------------(2) begin near-R1 point------------------
!!!debug         if(i_atom.eq.1.and.j_atom.eq.1.and. & 
!!!debug            i_coord.eq.2 .and. j_coord.eq.2 .and. & 
!!!debug            dabs(dir_tab_global(1,1,i_point)+0.023828440).lt.1.0d-4 .and. & 
!!!debug            dabs(dir_tab_global(2,1,i_point)).lt.1.0d-4 .and. &
!!!debug            dabs(dir_tab_global(3,1,i_point)).lt.1.0d-4 ) then 
!!!debug
!!!debug           write(use_unit,*) '' 
!!!debug           write(use_unit,*) '===========in evaluate_hellman_hellman_hessian(near-R1)=======' 
!!!debug           write(use_unit,*) 'dir-to-atom1:',dir_tab_global(1:3, 1, i_point)
!!!debug           write(use_unit,*) 'dir-to-atom2:',dir_tab_global(1:3, 2, i_point)
!!!debug           write(use_unit,*) 'belong to atom--->',atom_point(i_point) 
!!!debug
!!!debug           write(use_unit,*) 'hessian:' ,                      & 
!!!debug                    ! partition_tab(i_point)  *   &       
!!!debug           dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
!!!debug           first_order_rho_moving_grid(j_coord, i_atom,i_point)                     
!!!debug
!!!debug           write(use_unit,*) 'partition_tab:', partition_tab(i_point)
!!!debug
!!!debug         endif
!!!debug!--------------------(2) end near-R1 point------------------
     enddo !j_coord



      if(j_atom.ne.i_atom.and.atom_point(i_point).ne.j_atom) then ! fix R1,R2

        do j_coord=1,3 
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho(j_coord, j_atom,i_point) &
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point)
        enddo    !j_coord
       
      else if(j_atom.ne.i_atom.and.atom_point(i_point).eq.j_atom) then !moving R1,R2
        do j_coord=1,3 
           
           if(j_coord.eq.i_coord) then  
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho_moving_grid(j_coord, j_atom,i_point) &
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point)  & 
        -partition_tab(i_point) * rho(i_point) *                     &
        (species_z(species(i_atom))/(dist_tab_global(i_atom, i_point))**3.0d0 -      &
         3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  &
    dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
           else 
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho_moving_grid(j_coord, j_atom,i_point) &
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point)  & 
        -partition_tab(i_point) * rho(i_point) *                     &
       (-3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  &
    dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
           endif 

        enddo !j_coord

      else if(j_atom.eq.i_atom.and.atom_point(i_point).ne.j_atom) then! fix R1,R1 
        do j_coord=1,3 
           if(j_coord.eq.i_coord) then
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho(j_coord, j_atom,i_point) &
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point) & 
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
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point) & 
        +partition_tab(i_point) * rho(i_point) *                     &
       (-3.0d0*species_z(species(i_atom))*dir_tab_global(i_coord, i_atom, i_point)*  & 
    dir_tab_global(j_coord, i_atom, i_point)/(dist_tab_global(i_atom, i_point))**5.0d0)
          endif
 
        enddo ! j_coord


      else if(j_atom.eq.i_atom.and.atom_point(i_point).eq.j_atom) then ! moving R1,R1 
        do j_coord=1,3 
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) =    &
         hellman_feynman_hessian(i_coord,i_atom,j_coord,j_atom) +    & 
         partition_tab(i_point)  *   &       
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * & 
         first_order_rho_moving_grid(j_coord, j_atom,i_point) &
        +partition_deriv_local(j_coord,j_atom,i_point)* & 
         dir_tab_global(i_coord, i_atom, i_point)*atomic_term(i_atom) * &
         rho(i_point)
        enddo  !j_coord 

      else 
        write(use_unit,*) 'error for j_atom in evaluate_hellman_feynman_hessian'
        stop 
      endif 

    enddo !j_atom

  enddo ! i_coord
  enddo ! i_atom

  
 end do ! n_points
end subroutine evaluate_hellman_feynman_hessian
!******
