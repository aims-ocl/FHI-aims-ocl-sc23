# FHI-aims 2019 Release -- CHANGELOG
(Compilation in Progress; Describing new features, significant improvements)

FHI-aims, version 19----

### How to understand this CHANGELOG (and FHI-aims):

What we list here is a rough overview of recent changes in FHI-aims. A complete record of
what was done is available directly through the FHI-aims GitLab, to which every user is
welcome and encouraged to gain access. The GitLab repository is also the location where
bug fixes, improvements etc. will be made available.

To understand what this CHANGELOG promises and what it does not promise, it is important to
understand how FHI-aims is developed - that is, by a community of scientists in academia,
working on individual projects. There is no separation between "developers" and "users".

Thus, please understand the functionality listed below in this context. We
try to make things work. However, if any particular functionality does not work exactly as you would
expect, please be gentle on its original developer(s). Also research how the functionality
might perhaps work in an alternative way, or consider contributing a limited fix to a particular
issue back to FHI-aims. This will help the entire community. Thanks!

Changes since the release version 171221_1:


### FHI-aims at GitLab
(William Huhn, Victor Yu, ... who else?)

The source code has moved to: https://aims-git.rz-berlin.mpg.de

Please write to aims-coordinators@fhi-berlin.mpg.de to get access to the GitLab in case you have not so far.

Important features at the GitLab:

* Issue reporting
* Continuous Integration: Building + Regression test
* GUI for `git blame`, `git diff`, and other useful git features


### CMake for compiling FHI-aims
(Raul Laasner, Victor Yu, William Huhn, ... who else?)

Building of FHI-aims is now managed by CMake, which is a free and open-source build system generator.
The earlier (and alternative) "Makefile" based infrastructure will still be supported as a second option
for this release.

See the FHI-aims manual for more details and examples. Please also consider asking in our
Slack channel in case of any issues.


### Support of semilocal DFT simulations on GPU architectures
(William P. Huhn, Victor Yu, Raul Laasner ... who else?)

GPU implementation has been added to FHI-aims for select time-intensive operations via CUDA and the cuBLAS library. To use GPU acceleration with FHI-aims, one must compile the FHI-aims executable with GPU support, then enable GPU acceleration via keywords in the `control.in` file. More information may be found in the GPU appendix of the FHI-aims manual.

The following functionality has (apparently) stable GPU support:
* Electron density (when using density-matrix-based density update)
* Atomic forces and stress tensor (when using density-matrix-based density update)
* Hamiltonian and overlap matrices (semilocal DFT)
* NMR J-couplings and magnetic shielding tensors (done by Raul Laasner)

The following functionality has experimental GPU support:
* Solution of the Kohn-Sham eigenvalue equation via ELPA GPU (requires linkage to external/user-compiled ELPA library)

All other functionality currently does not have GPU support.

Note: It is highly recommended that users enable load balancing when running with GPU acceleration.


### Update of the electronic structure solvers in FHI-aims (including below-O(N^3) scaling)
(Victor Yu, William Huhn, others)

* Full support of density matrix solvers libOMM (cubic scaling), PEXSI (linear to quadratic scaling), and NTPoly (linear scaling) via ELSI
* New eigensolvers: SLEPc, EigenExa, MAGMA


### New and more stable restart infrastructure
(Victor Yu)

* MPI-IO based elsi_restart option for restarting a calculation - this is the recommended restart option in FHI-aims going forward


### JSON output format for FHI-aims
(William Huhn)

* JSON output framework via the FortJSON library.


### Improvement of Relaxation Algorithm (`relax_geometry trm`)
(Florian Knoop, Raul Laasner, others?)

* We updated the default algorithm for _periodic solids_ with unit cell relaxation in two ways:

1. We choose the initial guess for the hessian of lattice degrees of freedom such that is diagonal in _fractional coordinates_.
2. We optimize the cell in terms of strain transformations, i.e., we keep the fractional position of all atoms constant when moving the lattice in the geometry update.
3. The earlier algorithm (sometimes useful for deliberate symmetry breaking) is still available as trm_2012.
4. Without unit cell relaxation (i.e., if only atomic positions are changed) the earlier and new algorithms are identical.

* A long-standing issue with (very minor - 10^-11 or so) inconsistencies between geometries on different MPI tasks, for which
the code issued a warning, has been resolved.

* Output of final Hessian in `hessian.aims` file (instead of the end of `geometry.in`) (...)


### DFT+U extension
(Matthias Kick)

changed:
* occupation-matrix: default occupation matrix is now on-site (dual representation is still supported). turns out that the on-site occupation matrix tends to be more stable

added:
* occupation-matrix: on-site (default)
* forces: on-site occupation matrix has forces
* DFT+U matrix control: occupation-matrix can be initialized with predefined orbital occupancies, improves convergence for hard cases and enables convergence to distinct orbital configurations
* DFT+U ramping: the U value will be slightly increased during a calculation, improves convergence for hard cases
* enhanced localized projectors: arbitrary projector functions for calculating the occupation matrix can be used as long as they can be expressed as linear combination of basis functions

Please definitely see Matthias Kick's publication on how to apply DFT+U! (... ref ...)


### Frozen core approximation (experimental)
(Victor Yu)

This will make a difference for calculations that are (1) eigenvalue solver dominated
(solution of K.-S. equations) and (2) contain reasonably heavy elements, e.g., Ag, Au, Pt, ...
Large metallic slab systems with dense k-point grid are a particular application where
this may be helpful.

* free core orbitals when solving the Kohn-Sham eigenvalue problem
* Tested with single point total energy calculation, geometry relaxation, MD, periodic and non-periodic systems
* Usage: write `frozen_core_scf .true.` and `frozen_core_scf_cutoff -500` to freeze orbitals below -500 eV


### Real-space all-electron GW for periodic systems (experimental)
(Xinguo Ren)

* The code can in principle be applied to any gapped systems, but so far systematic benchmark tests have only been done for simple crystals.
* The current implementation is rather memory and CPU time intensive, but can be run in a massively parallel way.
* The FHI-aims-2009 tier 2 basis set is recommended to be used for production calculations, but one needs to be aware that, for certain systems (e.g., ZnO), the standard tier 2 basis alone might not be sufficient.
* In case there is any issue encountered in using the periodic GW code in FHI-aims, please contact Xinguo Ren (renxg@ustc.edu.cn).


### GW contour deformation treatment of the self-energy
(Dorothea Golze)

G0W0 with contour deformation (CD) for core states, but valence states with following features (JCTC paper from 2018)
* CD only for selected states, rest done with analytic continuation (AC)
* restart option for GW calculation
* print of self-energy matrix elements possible
* calculation of spectral function possible


### Parametrically-constrained relaxation
(Maja-Olivia Lenz)

* constrain relaxation to a user-specified set of parameters.
* relaxation is performed in this reduced space.
* allows for:
  - computational savings for high-symmetry structures
  - preserving crystal symmetry, e.g., to ensure that the relaxation remains within a particular structural prototype.
  - defining local symmetries and or local symmetry breakings, e.g., distortions.
* Usage: see the keywords `symmetry_n_params`, `symmetry_params`, `symmetry_lv` and `symmetry_frac` in the manual and the AFLOW Library of Crystallographic Prototypes, which is useful to generate this constraints.


### Other Changes and Improvements

* Intermediate Default Settings (accuracy and performance in between light and tight) (Volker Blum?)

* Performance Improvements DFPT polarizability/dielectric (Nathaniel Raimbault)

* Tetrahedron-Method for DOS - currently, serial only (Yi Yao)

* Improving performance and load balancing of exact exchange for heavy elements (`split_atoms`) (Sebastian Kokott)

* Simplified self-consistency for G0W0 using Hedin shift (CD and AC), see review paper (Dorothea Golze)

* Slater-Type orbitals (Raul Laasner)


# Waiting for input

* Magnetic Response – including GPU support (Raul Laasner)

* Libmbd (Integration of Many-Body-Dispersion Library, Jan Hermann)


### Work in progress :


### Improvements of Implicit Solvent Model for complex solute geometries (ongoing work)

* Cavity restart now produces both an `.xzy` file for visualization and an unformatted file used for actual restarts. Restarting an MPE calculation should be 100% exact now.
* Disabled forces for MPE. This was never implemented in the first place. Now aborts immediately after reading `control.in`.
* only allow `mpe_cavity_type rho_multipole_static` if one of the following is true:

  1. `mpe_skip_first_n_scf_steps` > 0 is given (edited) 
  2. A restart file (regular DFT restart) is given
  3. A cavity restart file is given (in this case `mpe_cavity_type` will have no effect unless it is `rho_multipole_dynamic`)

The reason is that people used to think that `rho_multipole_static` *without* restarting from a converged gas phase calculation was equivalent to `rho_free`, which it isn't  - the latter creates the isocavity from the superposition of free atom densities, the first from the density as it is in the first SCF cycle

* temporarily disabled `mpe_factorization_type qr`
* disable parallelization of the "walker dynamics" algorithm used in the cavity discretization **IF** `isc_kill_ratio`> 0 is given. With that keyword, parallelization led to inconsistent results with some builds (iirc with gnu compilers). Unfortunately, this slows down the algorithm, but I have not been able to trace down that bug yet.
(Jakob Filser)


### Fully-relativistic Implementation: self-consistent quasi-four-component relativity (experimental)

* New fully relativistic NAO basis set generated from solving the atomic Dirac equation
* New fully relativistic integrations. The code can perform fully relativistic calculations for band structures with a quasi four-component (Q4C) scheme at the cost of a two-component method level.
* Usage: write `spin collinear` and `relativistic 4c_dks` in control.in
(Rundong Zhao)

