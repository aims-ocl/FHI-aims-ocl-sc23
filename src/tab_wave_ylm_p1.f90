!****s* FHI-aims/tab_wave_ylm_p1
!  NAME
!   tab_wave_ylm_p1
!  SYNOPSIS

subroutine tab_wave_ylm_p1 &
     ( n_compute_atoms, atom_index, & 
     trigonom_tab, basis_l_max, l_ylm_max, &
     ylm_tab &
     )

!  PURPOSE
!  Subroutine tab_wave_ylm computes the Y_lm functions needed 
!  for wave functions between a 3D cartesian point and all atoms
!
!  USES

  use dimensions
  use geometry, only: species
  use pbc_lists
  implicit none

!  ARGUMENTS

  integer :: n_compute_atoms
  integer :: atom_index(n_compute_atoms)
  real*8 trigonom_tab ( 4, n_compute_atoms )
  integer basis_l_max (n_species)
  integer l_ylm_max
  real*8 ylm_tab( (l_ylm_max+1)**2, n_compute_atoms )

!  INPUTS
!    o n_compute_atoms -- number of relevant atoms
!    o atom_index -- list of relevant atoms
!    o trigonom_tab -- trigonometric functions from subroutine tab_trigonom_p0 
!    o basis_l_max -- maximum l index of basis functions
!    o l_ylm_max -- maximum l index of wave functions
!    
!  OUTPUT
!    o ylm_tab -- Y_lm functions
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
!



  !  counters

  integer i_atom

  !  begin work

  do i_atom = 1, n_compute_atoms, 1

     !       bootstrap table of ylm's
     call increment_ylm &
          ( trigonom_tab(1,i_atom), trigonom_tab(2,i_atom),  &
          trigonom_tab(3,i_atom), trigonom_tab(4,i_atom),  &
          0, basis_l_max(species(atom_index(i_atom))), &
          ylm_tab(1,i_atom) ) 

  enddo

  !  that's all folks

  return
end subroutine tab_wave_ylm_p1

!----------------------------------------------------------------------
!******
