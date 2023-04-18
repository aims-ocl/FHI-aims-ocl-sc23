!****s* FHI-aims/output_pi_molecular_dynamics
!  NAME
!    output_pi_molecular_dynamics
!  SYNOPSIS
subroutine output_pi_molecular_dynamics(free_energy)
!  PURPOSE
!    output routine for molecular dynamics
!    note that this routine outputs all the information for the LAST set of coordinates (not the freshly predicted ones!)
!    this is due to the fact that the generalized leap-frog algorithms don't determine the velocity at a given point until 
!    they actually predict the NEXT set of coordinates. 
!  USES

  use pi_molecular_dynamics
  use dimensions
  use species_data
  use geometry
  use runtime_choices
  use localorb_io
  use pbc_lists
  use timing
  use thermodynamic_integration
  use constants

  implicit none
!  ARGUMENTS
  real*8 :: free_energy
!  INPUTS
!    free_energy - energy from scf solver 
!          All other inputs are scavenged from the various modules in which they reside
!  OUTPUT
!    none
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
  integer :: i_atom, i_bead

  write(71,*) n_atoms*n_beads
  write(71,*)'comment'
  do i_atom = 1, n_atoms
    do i_bead = 1, n_beads
      write(71,900) species_name(species(i_atom)), coords_beads(:,i_bead,i_atom) * bohr
    end do
  end do

  write(72,*) n_atoms*n_beads
  write(72,*) 'comment'
  do i_atom = 1, n_atoms
    do i_bead = 1, n_beads
      write(72,900) species_name(species(i_atom)), v_beads_orig(:,i_bead,i_atom) * bohr
    end do
  end do
 
900   format(a20,3f16.9)
  
end subroutine output_pi_molecular_dynamics
!******
