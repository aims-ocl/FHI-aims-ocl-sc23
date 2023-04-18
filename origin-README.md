# The FHI-aims Electronic Structure Package #

<div align="center">
  <img src="aims-2010-11-01_800x800.png"><br><br>
</div>

FHI-aims is an efficient, accurate all-electron, full-potential electronic
structure code package for computational molecular and materials science
(non-periodic and periodic systems). The code supports DFT (semilocal and
hybrid) and many-body perturbation theory. FHI-aims is particularly efficient
for molecular systems and nanostructures, while maintaining high numerical
accuracy for all production tasks. Production calculations handle up to several
thousand atoms and can efficiently use (ten) thousands of cores.

FHI-aims is developed by an active, globally distributed community, including
significant developments at FHI, Duke University, TU Munich, USTC Hefei, Aalto
University, University of Luxembourg, TU Graz, Cardiff University and many other
institutions.

Specific functionality includes:

* Density functional theory (LDA, GGA, and mGGAs) for isolated molecules and
  periodic systems (solids, surfaces, ...)
* Preconstructed hierarchical basis sets across the periodic table (elements
  1-102) - from fast qualitative up to meV-converged accuracy
* Hartree-Fock and hybrid density functionals - tested up to 1,000 atoms in a
  linear scaling implementation (B3LYP, PBE0/PBEh, HSE)
* Range-separated long-range corrected hybrid density functionals (lc\_wPBEh)
* Approaches to van der Waals (Tkatchenko-Scheffler, vdW-DF, many-body
  dispersion (MBD))
* Many-body perturbation methods (currently non-periodic, for single-point
  geometries):
* MP2, RPA, renormalized second-order perturbation theory, G0W0, self-consistent
  GW, GW+SOSEX self-energy, doubly-hybrid functionals, and more under
  development
* Structure optimization, ab initio molecular dynamics, infrastructure for
  vibrations and phonons, ...
* Molecular transport (including "aitranss" maintained by Evers group in
  Regensburg)
* Seamlessly parallel from one up to (currently) (ten) thousands of CPUs
* FHI-aims' file formats are supported by the NOMAD repository, including data
  storage for 10 years, a citable Digital Object Indentifier (doi) for data
  sets, etc.

## Distribution of Source Code ##

FHI-aims is distributed in an "FHIaims" software package containing all files
needed to compile and run FHI-aims, including but not limited to:

* the source code needed to build an FHI-aims executable,
* test utilities to verify an FHI-aims executable,
* various utilities for plotting FHI-aims output,
* basis set definitions, and
* the FHI-aims manual.

We highly recommend that all users use the development version available on our
privately-hosted GitLab instance at
<https://aims-git.rz-berlin.mpg.de/aims/FHIaims>, as this version is the most
up-to-date.  Alternatively, the user may use pre-packaged release versions made
availably periodically.

## Installation ##

Instructions for installing FHI-aims may be found in the manual, which is a
standard LaTeX document that may be generated using the files found in the `doc`
folder of this package.

## Running FHI-aims ##

Instructions for running FHI-aims may be found in the manual, which is a
standard LaTeX document that may be generated using the files found in the `doc`
folder of this package.

## Contacting Us ##

We highly recommend that all FHI-aims users and developers join our Slack
channel at <https://fhi-aims.slack.com> .  The Slack channel is the fastest way
to get in contact with other FHI-aims users and developers.

* Slack account registration operates using a whitelist of approved e-mail
  domains.  If someone from your institution is already using the Slack channel,
  your institution is likely already on the whitelist, and you will be able to
  create an account on our Slack channel yourself.
* If you are unable to create an account because your institution's e-mail
  domain is not listed as one of our approved domains, contact us and we'll send
  you an invitation, as well as add your institution's e-mail domain to the
  approved list.

## Contributing to FHI-aims ##

Guidelines and suggested best practices for contributing to FHI-aims may be
found in the [CONTRIBUTING.md](CONTRIBUTING.md) file located in the root
directory of this repo.
