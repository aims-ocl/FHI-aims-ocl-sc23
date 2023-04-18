! ----------------------------------------------------------------------
!  Subroutine evaluate_moment_p1
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
subroutine evaluate_moment_p1(order,rho,partition_tab,electron_moment,ion_moment)

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
     ion_moment(:) = ion_moment(:) + &
          dble(species_z(species(i_atom))) * (coords(:,i_atom))**order
  end do
     
!     collective reduce operation, since rho is distributed
  call sync_moment(electron_moment)
  
  return
end subroutine evaluate_moment_p1


