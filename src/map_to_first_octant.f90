!****s* FHI-aims/map_to_first_octant
!  NAME
!    map_to_first_octant
!  SYNOPSIS

subroutine map_to_first_octant(coord)

!  PURPOSE
!   Maps the given coordinate point to the first octant (periodic systems).
!
!   ****************************************************************************************** 
!   **** If you make ANY changes here, please update "map_to_center_cell.f90" accordingly *** 
!   ****************************************************************************************** 
!
!  USES

  use dimensions
  use geometry
  implicit none

!  ARGUMENTS

  real*8:: coord(3)

!  INPUTS
!   o  coord -- coordinate point before mapping
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

  real*8:: length(3)

  if (n_periodic .eq. 0) return

! Atoms right on top of the cell limit will be moved to one 
! or the other side depending on a rounding error within 
! fortran. In order to have a more controlled shifting of 
! the atoms, a very small fraction is substracted from the 
! fractional coordinate and before converting it back 
! into real-space coordinates the small fraction is added 
! again, to ensure that atoms are not moved.

  length = matmul(map_to_center_cell_matrix, coord)
  length = length-1D-8
  length = length - dble(floor(length))+1D-8
  coord = matmul(lattice_vector, length)

end subroutine map_to_first_octant
!******
