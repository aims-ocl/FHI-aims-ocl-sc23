subroutine integrate_nlcorr(c_energy_pp)
  
  !calculate nonlocal correlation energy. SAG
  
  use octree_routines
  use physics, only: rho, rho_gradient, partition_tab
  use constants
  use dimensions
  use grids
  use runtime_choices
  
  
  implicit none 
  
  
  integer count2, i_full_points, i_my_batch, i_index, i_spin, i_coord
  real*8, dimension(3) :: coord_current
  real*8 epsilon, epsilon2, epsilon3, c_energy_pp, en_density_c_pp, en_nlcorr
  real*8 full_rho
  real*8, dimension(3) :: total_rho_gradient
  
  
  c_energy_pp = 0d0
  i_full_points = 0
  
  write(use_unit,*)" calculating nonlocal correlation energy"
  
  do i_my_batch = 1, n_my_batches, 1    !loop over all batches.
     
     do i_index = 1, batches(i_my_batch)%size, 1  !loop over one batch of integration pts. 
        
        epsilon = 0d0
        epsilon2 = 0d0
        epsilon3 = 0d0
        en_density_c_pp = 0d0
        full_rho = 0d0
        total_rho_gradient(:) = 0d0
        
        i_full_points = i_full_points + 1
        
        coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
                   
        
        if (spin_treatment.eq.1) then
           full_rho = rho(1,i_full_points) + rho(2,i_full_points)
           do i_coord = 1, 3, 1
              total_rho_gradient(i_coord) = &
                   rho_gradient(i_coord,1,i_full_points)+rho_gradient(i_coord,2,i_full_points)
           enddo
        else
           full_rho = rho(1,i_full_points)
           do i_coord = 1,3,1
              total_rho_gradient(i_coord) = rho_gradient(i_coord,1,i_full_points)
           enddo
        end if
                
        
        call nlcorr(coord_current,full_rho,total_rho_gradient(:),epsilon,epsilon2, epsilon3) 
        
        en_density_c_pp = epsilon
        
        
        do i_spin = 1,n_spin                  
           c_energy_pp = c_energy_pp + &
                rho(i_spin,i_full_points)* &
                en_density_c_pp * &
                partition_tab(i_full_points)
        enddo
        
     enddo
  enddo
  
  
  
  
  
  
  
end subroutine integrate_nlcorr
