!****s* FHI-aims/output_potential_p1
!  NAME
!    output_potential_p1
!  SYNOPSIS

subroutine output_potential_p1( local_potential_parts )

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
  use localorb_io, only: use_unit
  use mpi_tasks
  use mpi_utilities
  implicit none

!  ARGUMENTS

  real*8, dimension(n_spin, n_full_points) :: local_potential_parts

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

!  counters

!  integer i_atom, i_atom_2
!  integer i_radial
!  integer i_angular
!  integer i_grid
  integer i_coord
  integer i_spin
  
!  integer i_species

  integer :: i_point, i_my_batch, i_index

!  begin work


  if (myid.eq.0) then
     write(use_unit,'(2X,A,A)') &
          "Writing the effective potential at each grid point", &
          " to file v_eff.dat ."
     open (50, file="v_eff.dat")
  end if

  i_point = 0

!     calculate potential, atom by atom, and point by point 
  do i_my_batch = 1, n_my_batches, 1

        do i_index = 1, batches(i_my_batch)%size, 1
           
           i_point = i_point + 1
           
           !           calculate grid point coordinate
           grid_coord(:) = batches(i_my_batch) % points(i_index) % coords(:)
           
           !           write potential
           
           if (myid.eq.0) then
              write(50,'(7(F20.10,1X))') &
                   (grid_coord(i_coord), i_coord=1,3,1), &
                   ( local_potential_parts(i_spin,i_point), &
                   i_spin = 1, n_spin, 1 )
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
end subroutine output_potential_p1
!----------------------------------------------------------------------
!******	
