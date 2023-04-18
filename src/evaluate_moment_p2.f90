! ----------------------------------------------------------------------
!  Subroutine evaluate_moment_p2
!
!  Notes: * MPI implementation not checked
!  
!  Parameter:      
!     order:           order of moment
!     rho:             electron density
!     partition_tab:   partition tab for integration
!     electron_moment: (output) electrons' moment
!     ion_moment     : (output) ions' moment
!----------------------------------------------------------------------
!  This is a clone to evaluate_moment_p1, targeted towards periodic systems
!  The main difference is that is maps all grid points to a single cell, i.e, 
!  between the "upper" and the "lower" vacuum level
!  Note that the moments are calucated in all three directions, but
!  unless a vacuum-level in x and y are supplied (which FHI-aims doesn't 
!  support yet), only the z-component makes any sense. 
!----------------------------------------------------------------------

subroutine evaluate_moment_p2(order,rho,partition_tab,electron_moment,ion_moment)

  use runtime_choices, only: vacuum_z_level
  use dimensions
  use grids
  use geometry
  use spline
  use species_data
  use free_atoms
  use mpi_tasks
  use synchronize_mpi
  use localorb_io

  implicit none

!  input
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(n_full_points) :: partition_tab
  real*8 order

!  output
  real*8 :: electron_moment(3)
  real*8 :: ion_moment(3)

!  local variables
  real*8 :: grid_coord(3)
  integer :: i_spin
  integer :: i_coord
  integer :: i_full_points
  integer :: i_my_batch, i_atom, i_index
  real*8 :: coord_tmp(3)


! begin work
  i_full_points = 0
  grid_coord = 0.0d0
  electron_moment = 0.0d0
  ion_moment = 0.0d0

!    calculate dipole moment looping over the grid
  do i_my_batch = 1, n_my_batches, 1

     ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1               

!    coordinate in cartesian system 
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)

           !Make sure grid is between lower and upper vacuum level
           do while(grid_coord(3).lt.(vacuum_z_level-maxval(lattice_vector(3,:))))
                 grid_coord(3)=grid_coord(3)+maxval(lattice_vector(3,:))
           enddo
           do while(grid_coord(3).gt.vacuum_z_level)
                 grid_coord(3)=grid_coord(3)-maxval(lattice_vector(3,:))
           enddo

!    compute electron moment
           do i_spin = 1, n_spin, 1
              electron_moment(:)=electron_moment(:)+ &
                   rho(i_spin, i_full_points) * &
                   partition_tab(i_full_points)*(grid_coord(:))**order
!debug                   partition_tab(i_full_points)*(grid_coord(:))
           enddo

!    end loop over a batch
        enddo

!    end of work distribution
     ! endif
             
!    end loop over batches
  end do

!    compute ion moment
  do i_atom = 1, n_atoms, 1
     !Make sure atom is between lower and upper vacuum level
     coord_tmp(:)=coords(:,i_atom)
     do while(coord_tmp(3).lt.(vacuum_z_level-maxval(lattice_vector(3,:))))
           coord_tmp(3)=coord_tmp(3)+maxval(lattice_vector(3,:))
     enddo
     do while(coord_tmp(3).gt.vacuum_z_level)
           coord_tmp(3)=coord_tmp(3)-maxval(lattice_vector(3,:))
     enddo

     ion_moment(:) = ion_moment(:) + &
          dble(species_z(species(i_atom))) * (coord_tmp(:))**order
  end do
     
!     collective reduce operation, since rho is distributed
  call sync_moment(electron_moment)
  
  return
end subroutine evaluate_moment_p2


