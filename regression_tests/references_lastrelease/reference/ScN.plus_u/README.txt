VB, 02.01.2011:

Added a regression test for periodic "LDA+U" in its present state, despite
some misgivings since the author (Norbert Nemec) is no longer actively
around. The express intention of this regression test is to ensure that the
working part of the LDA+U functionality remains in place and unbroken until
the development is picked up again. 

I reduced the original test case by Norbert (shown in the new "examples"
directory) in the following way:

k_grid 10 10 10 -> 6 6 6
cut_pot 3.0 1.0 1.0
sc_accuracy_rho 1e-5 -> 1e-4
sc_accuracy_etot 1e-6 -> 1e-5

This speeds up the test to approx 1 minute even on my desktop (virtual machine
on laptop). 

The caveat is that already now (code version 122110), I see that results
differ slightly but noticeably between 2 cores and 8 cores.

There is also a slight difference when using 2 cores but slightly different
versions of ifort, mkl, and environment (my desktop vs. theows17) but this
difference is much smaller.

All observed differences _could_ be due to the sensitive nature of projection
operations, which might act as a "noise enhancer" when coupled with the
nonlinear scf procedure. I presently have no opinion on the matter, but record
this information here for whoever picks up LDA+U again. 


*** Original comment by Norbert Nemec:

Compare e.g. with PhysRevB.62.13539:

Computed with plain LDA, the ScN system does not feature any band gap. LDA+U
allows to split the occupied and the unoccupied bands and reproduce the
experimentally observed gap. 

(The bands.svgz is produced with the utilities/aimsplot.py script.)
