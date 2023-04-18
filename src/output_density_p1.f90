!****s* FHI-aims/output_density_p1
!  NAME
!   output_density_p1
!  SYNOPSIS

subroutine output_density_p1( rho, partition_tab )

!  PURPOSE
!
!  Subroutine output_density
!  Writes the total density and density minus free-atom density on the integration grid.   
!
!  USES
  use dimensions
  use grids
  use geometry
  use spline
  use free_atoms
  use mpi_tasks
  use localorb_io
  use constants
  use mpi_utilities

  implicit none

!  ARGUMENTS
  real*8, dimension(n_spin, n_full_points) :: rho
  real*8, dimension(n_full_points) :: partition_tab

!  INPUTS
!  o rho - electronic density
!  o partition_tab - the partition tab...
!
!  OUTPUTS
!  density on the integration grid
!
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals:
!    FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject
!   to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE


! local variables
  
  real*8 :: grid_coord(3)
  real*8 :: distance
  real*8 :: delta_rho
  real*8 :: i_r
  real*8 :: norm_density
  real*8 :: norm_diff_density
  real*8 :: norm_free_density

! counters

  integer :: i_atom
  integer :: i_atom_2
  integer :: i_radial
  integer :: i_angular
  integer :: i_coord
  integer :: i_point
  integer :: i_spin
  integer :: i_my_batch, i_index

!  functions

  real*8 :: get_distance

!  begin work

  
  norm_density = 0.
  norm_diff_density = 0.
  norm_free_density = 0.
  
  call localorb_info('',use_unit)
  call localorb_info( &
       "Writing the density and difference density at each integration", &
       6,'(2X,A)') 
  call localorb_info( &
       "point to file density.dat and diff_density.dat .", &
       6,'(2X,A)') 

  if (myid.eq.0) then
     open (50, file="diff_density.dat")
     open (60, file="density.dat")
  end if

  i_point = 0

!    calculate potential, atom by atom, and point by point 
  do i_my_batch = 1, n_my_batches, 1
     
        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_point = i_point + 1               

!           calculate grid point coordinate
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)

           delta_rho = 0.d0
           do i_spin = 1, n_spin, 1
              delta_rho = delta_rho + rho(i_spin, i_point)
           enddo
           
!     subtract the free atom density contribution                  
               
           do i_atom_2 = 1, n_atoms, 1
              
              distance = get_distance(grid_coord, &
                   coords(1, i_atom_2))
              i_r = &
                   invert_log_grid &
                   (distance, &
                   r_grid_min(species(i_atom_2)), &
                   r_grid_inc(species(i_atom_2)) )
              
              delta_rho = delta_rho - & 
                   pi4_inv * val_spline(i_r, free_rho_spl &
                   (1, 1, species(i_atom_2)), &
                   n_grid(species(i_atom_2)) )
              
              norm_free_density = norm_free_density + &
                   pi4_inv * val_spline(i_r, free_rho_spl &
                   (1, 1, species(i_atom_2)), &
                   n_grid(species(i_atom_2)) ) * &
                   partition_tab(i_point) 
              
           end do
           
           do i_spin = 1, n_spin, 1
              norm_density = norm_density + &
                   rho(i_spin, i_point) * partition_tab(i_point)
           enddo
           
           norm_diff_density = norm_diff_density &
                + delta_rho * partition_tab(i_point)
           
           !     write density and difference density

           if (myid.eq.0) then
              write(50,'(4(F20.10,1X))') &
                   (grid_coord(i_coord), i_coord=1,3,1), &
                   delta_rho
              write(60,'(5(F20.10,1X))') &
                   (grid_coord(i_coord), i_coord=1,3,1), &
                   (rho(i_spin, i_point), i_spin = 1, n_spin, 1)
           end if
           
           ! end loop over a batch
        end do
        
     !end if
     
     ! end loop over grid batches
  end do
      
  if (myid.eq.0) then 
     write(use_unit,'(2X,A, F10.4)')  "Norm of difference density: ", &
          norm_diff_density
     write(use_unit,'(2X, A, F10.4)') "Norm of full density      : ", &
          norm_density
     write(use_unit,'(2X, A, F10.4)') "Norm of free-atom density : ", &
          norm_free_density
     close(50)
     close(60)
  end if

  return
end subroutine output_density_p1

!----------------------------------------------------------------------
!******
