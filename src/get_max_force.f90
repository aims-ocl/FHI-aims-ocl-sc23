!****f* FHI-aims/get_max_force
!  NAME
!   get_max_force
!  SYNOPSIS

real*8 function get_max_force(forces)

!  PURPOSE
!  The function finds out the maximum force component.
!
!  USES

  use dimensions
  implicit none

!  ARGUMENTS
    
  real*8, dimension(3, n_atoms) :: forces

!  INPUTS
!  o forces -- forces
!   
!  OUTPUT
!  o get_max_force -- maximum force component
!
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



  ! local variables
  real*8 :: max_force

  ! counter
  integer :: i_coord
  integer :: i_atom 

  max_force = 0.d0
  do i_atom = 1, n_atoms, 1
     do i_coord = 1, 3, 1
        if (max_force.lt.abs(forces(i_coord, i_atom))) then
           max_force = abs(forces(i_coord, i_atom))
        end if
     enddo
  enddo

  get_max_force = max_force

end function get_max_force
!******	
