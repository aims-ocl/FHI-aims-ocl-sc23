subroutine get_nlcorr(c_energy_pp)

!adds nonlocal correlation energy to total_energy in scf_solver
!XC energies and forces are NOT updated.   SAG

use octree_routines
use physics, only: rho, rho_gradient, partition_tab
use constants
use dimensions
use grids


implicit none 
 

integer count2, i_full_points, i_my_batch, i_index, i_spin
real*8, dimension(3) :: coord_current
real*8 epsilon, epsilon2, epsilon3, c_energy_pp, en_density_c_pp, en_nlcorr
!real*8 total_energy

character*100 info_str

             write(info_str,'(2X,A)') "VDW post-processed correction:"
             call localorb_info(info_str, use_unit,'(A)', OL_norm)


             if(associated(root_full%branch)) then  !Deallocate old tree if there is one.
                !write(use_unit,*)"Deallocating old tree" 
                count2 = 0
                call deallocate_tree(root_full,count2)
                !write(use_unit,*)"Number of branches deallocated (from root_full):",count2
             endif
             
             
             
             !initialize master cube, and go ahead and let it generate tree. 
             !treeflag = 3 !Set to 1 to use multipoles
             !write(use_unit,*)"Initializing master cube: (region of vdw interactions)"     

             epsfinal = 0.000001d0  !default
             call gen_master_cube()
             
             !calculate nlcorr in loop over all AIMS grid points.
             write(info_str,'(2X,A)') "calculating nonlocal correlation energy..."
             call localorb_info(info_str, use_unit,'(A)', OL_norm)
             c_energy_pp = 0d0
             i_full_points = 0
             do i_my_batch = 1, n_my_batches, 1    !loop over all batches.
                
                do i_index = 1, batches(i_my_batch)%size, 1  !loop over one batch of integration pts. 
                   
                   i_full_points = i_full_points + 1
                   
                   coord_current(:) = batches(i_my_batch) % points(i_index) % coords(:)
                              
                   epsilon = 0d0
                   epsilon2 = 0d0
                   epsilon3 = 0d0
                   en_density_c_pp = 0d0
                   
                   call nlcorr(coord_current,rho(1,i_full_points),   &
                        rho_gradient(:,:,i_full_points),epsilon,epsilon2, epsilon3) 
                  
                   en_density_c_pp = epsilon
                   
                   do i_spin = 1,1                  
                      c_energy_pp = c_energy_pp + &
                           rho(i_spin,i_full_points)* &
                           en_density_c_pp * &
                           partition_tab(i_full_points)
                   enddo
                   !write(21,*)"c_energy_pp", c_energy_pp
                   
                enddo
             enddo
                      
            
!!$             write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
!!$                  "Total energy, before vdw correction  ", total_energy, " Ha   ", total_energy*hartree, " eV"
!!$             total_energy = total_energy + c_energy_pp
!!$            
!!$             write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
!!$                  "Nonlocal correlation energy          ", c_energy_pp, " Ha   ", c_energy_pp*hartree, " eV"
!!$             write(use_unit,'(2X,A,F20.9,A,F20.9,A)') &
!!$                  "Total energy, now with vdw correction", total_energy, " Ha   ", total_energy*hartree, " eV"
!!$             write(use_unit,*)" "
             !write(use_unit,*)"Nonlocal correlation energy now included in total energy."
             !write(use_unit,*)"However, not included in XC analysis below. Nor in forces."
             !write(use_unit,*)" "
  
           




end subroutine get_nlcorr
