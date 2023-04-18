!****s* FHI-aims/estimate_Ewald_radius
!  NAME
!    estimate_Ewald_radius
!  SYNOPSIS

subroutine estimate_Ewald_radius( n_atoms_ewald, volume, r0 )

  !  PURPOSE
  !
  !    Empirical estimate of the range separation parameter r0 in the 
  !    Ewald method for electrostatic energies in periodic systems.
  !
  !    r0 is the "hartree_convergence_parameter" keyword value in control.in .
  !
  !    The routine is deliberately not tied into the existing module
  !    infrastructure of FHI-aims since it is very short and thus more easily independent.
  !
  !    The optimum value of r0 could be determined much more rigorously and may be somewhat
  !    machine- and/or compiler- and optimization-dependent. 
  !
  !    However, here, our main objective is 
  !    not to determine an optimal choice for all times, but rather to limit some damage that
  !    arises if a too small r0 value is set for surface slab calculations with a large vacuum
  !    region. 
  !
  !    The empirical formula used here was determined for a specific system and is 
  !    NOT OPTIMAL for all systems. The test system was a dense slab of a carbon-nitride 
  !    material, with roughly 1000-atom per unit cell, varying the vacuum thickness.
  !    The result, expressed as a function of volume per atom, is
  !    as follows:
  !
  !    r0 = A0 * exp[ 1/3 * ln(v - A1) ]
  !
  !    i.e., an assumed cubic root behavior (no analytic reason), where 
  !
  !    A0 = 1.47941 bohr
  !    A1 = 1.85873 A^3
  !
  !    and v = volume per atom in A^3 .
  !
  !    This curve hits 3.0 bohr at somewhere between 9 and 10 Angstrom^3/atom .
  !
  !    The current minimum allowed r0 value is 2.5 bohr .
  !    The current maximum allowed r0 value is 5.0 bohr .
  !
  !    This range is chosen to prevent estimates that are accidentally out of bounds
  !    for other systems.
  !
  !    Specifically, volume per atom may not be the only relevant quantity. For instance,
  !    the choice of r0 according to this formula appears to be somewhat too large for
  !    a not-so-dense bulk material like GaAs, but it appears that the overall time for
  !    the Hartree potential does not become a dominant term in this case.
  !
  !    For larger vacuum regions (80 - 100 AA or even beyond), larger r0 values may well be
  !    beneficial, but it should be ensured explicitly that the resulting total energies 
  !    are identical (at sub-meV level) for a smaller, more expensive choice of r0.
  !
  !    If this is tested by anyone, please consider refining the default r0 estimate made here,
  !    especially 
  !    - towards larger vacuum spacings where the Ewald part of the Hartree potential
  !      may become computationally dominant.
  !    - for more bulk-like systems with reasonably low density, like GaAs or molecular crystals,
  !      where the simple volume formula may overestimate r0.
  !
  !    Note that 'output_level full' will, among other things, print out separate timings for
  !    the Ewald part of the Hartree potential.
  !
  !  USES

  use constants, only : bohr
  use localorb_io
  implicit none

  !  ARGUMENTS

  integer, intent(IN)   :: n_atoms_ewald
  real*8, intent(IN)    :: volume
  real*8, intent(OUT)   :: r0

  !  INPUTS
  !    o n_atoms_ewald -- number of atoms (usually in unit cell)
  !    o volume -- unit cell volume (in atomic units, i.e., bohr^3)
  !  OUTPUTS
  !    o r0 - Ewald range separation parameter, in atomic units, i.e., in bohr
  !  AUTHOR
  !    Volker Blum, Duke University
  !  SEE ALSO
  !    Volker Blum, Ralf Gehrke, Felix Hanke, Paula Havu, Ville Havu,
  !    Xinguo Ren, Karsten Reuter, and Matthias Scheffler,
  !    "Ab initio simulations with Numeric Atom-Centered Orbitals: FHI-aims",
  !    Computer Physics Communications 180, 2175 (2009).
  !  COPYRIGHT
  !    VB.
  !    Please note that any use of the "FHI-aims-Software" is subject to
  !    the terms and conditions of the respective license agreement.
  !  HISTORY
  !    Release version, FHI-aims (2016).
  !  SOURCE

  ! local variables

       real*8 :: A0, A1
       real*8 :: volume_per_atom
       character*120 :: info_str

       A0 = 1.47941  ! bohr
       A1 = 1.85873  ! A^3

       volume_per_atom = volume / dble(n_atoms_ewald) * bohr**3  ! in A^3

       if ((volume_per_atom-A1).gt.1.d-12) then
         r0 = A0 * (volume_per_atom - A1)**0.333333333333333d0
       else
         r0 = 3.d0 ! bohr
       end if

       ! boundaries
       r0 = max(r0, 2.5d0)
       r0 = min(r0, 5.0d0)

  ! Just before end of routine, write result.

       call localorb_info('  ', use_unit, '(A)')
       write(info_str,'(2X,A,F15.8,A)') & 
      "Range separation radius for Ewald summation (hartree_convergence_parameter): ", & 
       r0, " bohr. "
       call localorb_info(info_str,use_unit,'(A)')        

end subroutine estimate_Ewald_radius
!******
