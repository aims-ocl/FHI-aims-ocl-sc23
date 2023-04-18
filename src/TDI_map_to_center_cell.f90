!****s* thermodynamic_integration/TDI_map_to_center_cell
!  NAME
!    TDI_map_to_center_cell
!  SYNOPSIS

subroutine TDI_map_to_center_cell(point,i_atom)

!  PURPOSE
!   Maps the given coordinate POINTS to the center cell rount origin (periodic systems)
!   Additionally, it translates the equilibrium positions for the thermodynamic integration
!   by the EXACT same vector
!
!  USES

  use dimensions
  use geometry
  use pbc_lists
  use runtime_choices
  use mpi_tasks
  use synchronize_mpi
  use thermodynamic_integration
  implicit none

!  ARGUMENTS
  real*8  :: point(3)
  integer :: i_atom

!  INPUTS
!   o  coord  -- coordinate point before mapping
!   o  i_atom -- index of the atom
!  OUTPUT
!   o  coord -- coordinate point after mapping
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

  real*8:: old_point(3)
  integer :: segments

  if (n_periodic .eq. 0) return
  
  ! Save old position
  old_point = point

  ! Map to center cell
  call map_to_center_cell(point)
  
  ! Compute diff vector
  old_point =  point - old_point
   

  !Shift Equil. positions accoridngly
  !For actual coordinates
  TDI_Segment_atoms(:,i_atom) =  TDI_Segment_atoms(:,i_atom) + old_point

  !For the rest of the segments
  do segments=1,TDI_segments 
    TDI_atoms(segments,:,i_atom)  =   TDI_atoms(segments,:,i_atom) + old_point
  end do

end subroutine TDI_map_to_center_cell
!******
