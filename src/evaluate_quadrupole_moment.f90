! ----------------------------------------------------------------------
!  Subroutine evaluate_quadrupole_moment
!  Quadrupole moment ist calculated following the definition in classical electrostatics, as you can find it e.g. in Jackson's textbook
!
!  Notes: * MPI implementation not checked
!  
!  Parameter:      
!     order:           order of moment, not used
!     rho:             electron density
!     partition_tab:   partition tab for integration
!     electron_quadrupole_moment: (output) electrons' moment
!     ion_quadrupole_moment     : (output) ions' moment
!----------------------------------------------------------------------
subroutine evaluate_quadrupole_moment(order,rho,partition_tab,electron_quadrupole_moment,ion_quadrupole_moment)

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
  real*8, dimension(3,3) :: electron_quadrupole_moment
  real*8, dimension(3,3) :: ion_quadrupole_moment

!  local variables
  real*8 :: grid_coord(3) 
  real*8 :: r_square
  integer :: i_spin
  integer :: i_coord
  integer :: j_coord
  integer :: i_full_points
  integer :: i_my_batch, i_atom, i_index


! begin work
  i_full_points = 0
  grid_coord = 0.0d0
  electron_quadrupole_moment = 0.0d0
  ion_quadrupole_moment = 0.0d0

!    calculate quadrupole moment looping over the grid
  do i_my_batch = 1, n_my_batches, 1

     ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1

           i_full_points = i_full_points + 1               

!    coordinate in cartesian system 
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)
	   r_square = grid_coord(1)**2+grid_coord(2)**2+grid_coord(3)**2
	   

!    compute electron quadrupole tensor
		do i_spin = 1, n_spin, 1
			do i_coord = 1, 3, 1
				do j_coord = 1, i_coord, 1
					if (i_coord==j_coord) then
						electron_quadrupole_moment(i_coord,j_coord) &
						=electron_quadrupole_moment(i_coord,j_coord)+& 
						rho(i_spin, i_full_points) * &
						partition_tab(i_full_points)* (   &
						3.0*grid_coord(i_coord)*grid_coord(j_coord)-&
						r_square)
					else 
						electron_quadrupole_moment(i_coord,j_coord) &
						=electron_quadrupole_moment(i_coord,j_coord)+&
						rho(i_spin, i_full_points) * &
						partition_tab(i_full_points)*   &
						3.0*grid_coord(i_coord)*grid_coord(j_coord)
					endif
				enddo
			enddo
		enddo


!    end loop over a batch
   enddo
!    end of work distribution
     ! endif
             
!    end loop over batches
  end do

     
!     collective reduce operation, since rho is distributed
  call sync_matrix(electron_quadrupole_moment,3,3)
  
!    compute  ion quadrupole tensor
	do i_atom= 1, n_atoms, 1
		r_square = coords(1,i_atom)**2+coords(2,i_atom)**2+coords(3,i_atom)**2
		do i_coord = 1, 3, 1
			do j_coord = 1, i_coord, 1
				if (i_coord==j_coord) then
					ion_quadrupole_moment(i_coord,j_coord) &
					=ion_quadrupole_moment(i_coord,j_coord)+&
						dble(species_z(species(i_atom))) *( &
					3.0*coords(i_coord,i_atom)*coords(j_coord,i_atom)-&
					r_square)
				else 
					ion_quadrupole_moment(i_coord,j_coord) &
					=ion_quadrupole_moment(i_coord,j_coord)+ &
						dble(species_z(species(i_atom))) *( &
					3.0*coords(i_coord,i_atom)*coords(j_coord,i_atom))
				endif
			enddo
		enddo
	enddo


  return
end subroutine evaluate_quadrupole_moment
