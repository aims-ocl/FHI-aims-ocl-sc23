!****s* FHI-aims/output_potential_p1
!  NAME
!    output_potential_p1
!  SYNOPSIS

subroutine output_rho( local_potential_parts)!, hartree_partition_tab )

!  PURPOSE
!  Subroutine output_potential
!
!  Writes the effective potential on the integration grid. In practise it prints out the Hartree potential part.
!
!  Note: This is not the effective potential for  a real fix, we must evaluate the 
!  xc potential right here so that we can actually write out the relevant bits.
!
!  USES

  use dimensions
  use runtime_choices
  use grids
  use geometry
  use mpi_tasks
  use mpi_utilities
  use species_data
  use constants
  use grids, only: w_angular
  use localorb_io, only: use_unit
  implicit none

!  ARGUMENTS

  real*8, dimension(n_full_points) :: local_potential_parts
 !real*8, dimension(n_full_points)      :: hartree_partition_tab

!  INPUTS
!   o local_potential_parts -- potential which is printed out.
!  OUTPUT
!    none
!  AUTHOR
!    FHI-aims team, Fritz-Haber Institute of the Max-Planck-Society
!  SEE ALSO
!    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
!    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
!    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
!    Computer Physics Communications (2008), submitted.
!  COPYRIGHT
!   Max-Planck-Gesellschaft zur Foerderung der Wissenschaften
!   e.V. Please note that any use of the "FHI-aims-Software" is subject to
!   the terms and conditions of the respective license agreement."
!  HISTORY
!    Release version, FHI-aims (2008).
!  SOURCE



!  local variables

  real*8 grid_coord(3)
  integer :: current_atom, current_radial,current_angular

!  counters

!  integer i_atom, i_atom_2
!  integer i_radial
!  integer i_angular
!  integer i_grid
  integer i_coord
  integer i_spin, i_grid
  
!  integer i_species

  integer :: i_point, i_my_batch, i_index, i_atom, i_grids

!  begin work



  if (myid.eq.0) then
     write(use_unit,'(2X,A,A)') &
          "Writing the linear Poisson-Boltzmann potential at each grid point", &
          " to file rho.dat ."
     open (50, file="rho.dat")
  end if



!   write(50,'(4(i7))') n_max_grid,n_max_radial,n_full_points, n_atoms
! 
!   do i_atom = 1, n_atoms, 1
!     write(50,'(4(F25.14,1X),2(i6,1X),3(F25.14,1X),(F25.14,1X))') r_grid_min(species(i_atom)),&
! 	 r_grid_inc(species(i_atom)), scale_radial(species(i_atom)), multipole_radius_free(species(i_atom)),&
! 	 n_radial(species(i_atom)), n_grid(species(i_atom)),  coords(1:3,i_atom), species_z(species(i_atom))
!     do i_grids = 1, n_max_grid, 1
!       write(50,'(F25.15)') r_grid(i_grids, species(i_atom))
!     end do
!   end do

  i_point = 0  
  
!     calculate potential, atom by atom, and point by point 
  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1

           i_point = i_point + 1

           !           calculate grid point coordinate
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)

           current_atom    = batches(i_my_batch) % points(i_index) % index_atom
	   current_radial  = batches(i_my_batch) % points(i_index) % index_radial
           !           write potential
	   
!            if (myid.eq.0) then
! 
!               write(50,'(2(i6,1X),6(F25.15,1X))') &
!                    current_atom, current_radial,r_radial(current_radial,species(current_atom)),&
! 		    hartree_partition_tab(i_point), (grid_coord(i_coord), i_coord=1,3,1), local_potential_parts(i_point)
!            end if

           if (myid.eq.0) then
              write(50,'(4(F20.10,1X))') &
                   (grid_coord(i_coord), i_coord=1,3,1), &
                   local_potential_parts(i_point)
                   
           end if
           
           ! end loop over a batch
        enddo
        
     !end if
     ! end loop over batche
  enddo

      
  if (myid.eq.0) then
     close(50)
  end if

  return
end subroutine output_rho
!----------------------------------------------------------------------
!******	
