!****s* FHI-aims/update_grid_coords
!  NAME
!   update_grid_coords
!  SYNOPSIS

subroutine update_grid_coords()

!  PURPOSE
!  The subroutine updates the grid coordinates. In every grid point the information
!  of atom, angular grid, and radial grid indexes are already saved. Now the actual
!  coordinates are calculated from these information.
!
!  USES

  use dimensions
  use grids
  use geometry
  use mpi_utilities
  implicit none

!  ARGUMENTS
!  INPUTS
!    none
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






  integer :: i_my_batch, i_index
  integer :: current_atom, current_radial, current_angular

  do i_my_batch = 1, n_my_batches, 1

        ! loop over one batch
        do i_index = 1, batches(i_my_batch)%size, 1
           
           current_atom = batches(i_my_batch)%points(i_index)%index_atom
           current_radial = batches(i_my_batch)%points(i_index)%index_radial
           current_angular = batches(i_my_batch)%points(i_index)%index_angular

           batches(i_my_batch)%points(i_index)% coords(:) = &
                coords(:, current_atom) + &
                r_angular(:, current_angular, current_radial,species(current_atom)) * &
                r_radial(current_radial, species(current_atom))

        end do
           
     ! end if
  end do

end subroutine update_grid_coords
!******	
