"species_defaults" subdirectory LVL_fast

This is a directory of incomplete (experimental) species defaults that
documents some initial experiences with harmonic (finite difference)
molecular vibrations and ab initio molecular dynamics based on the
PBE0 functional and the "LVL_fast" method of expanding the 2-electron
Coulomb operator.

For hybrid functionals, LVL_fast comes very close (meV or sub-meV) to
the accuracy of RI-V in test calculations with light-element
molecules. RI_method LVL_fast is also the implementation which allows
exact exchange with gradients and/or periodic boundary conditions. 

However, the construction of the auxiliary basis set to expand the
Coulomb operator benefits from some modified choices in the species
defaults:

* The auxiliary basis must not be too small in order to retain
  accuracy; however, for "tier 1" or smaller (for light elements), some
  useful high-angular momentum components are not yet there and must
  be added by hand. (using the for_aux) keyword

* The auxiliary basis set must not be too large because a specific
  matrix inversion in the gradient calculation can introduce
  significant errors if applied to a nearly ill-conditioned
  matrix. This can happen already for "tier 2" basis sets with light
  elements, and is easily caught by limiting the number of auxiliary
  basis functions used per atom to a sensible value. 

* Finally, Hartree-Fock like exchange will scale (numerically) as
  O(N^4) as long as one adds basis functions on top of one
  another. Thus, while scalability with system size is much better,
  the step from tier 1 (light) to tier 2 (tight) is unusually
  expensive for hybrid functionals, and an intermediate level of
  accuracy may be warranted specifically for hybrid functionals. This
  level retains the grid, l_hartree and cut_pot settings at the
  "tight" level but reduces the overall number of basis functions
  included to just above the normal "light" standard. 

The species defaults provided here are absolutely unfinished at
present, in the sense that practically no elements exist except for
the ones which we tested explicitly.

Please take the present status as a the first step of documenting
optimized input parameters and contribute positive or negative
experiences. The goal is to provide a complete list of optimum
defaults in the future. 
