Compilation instructions for NEB script + code:

In order to make this run both on serial and parallel systems, there are some
real and some fake MPI routines included in these modules, that have to be
linked into the main NEB program via the provided Makefile.

The program NEB.pl requires editing before it can be executed in parallel 
due to the missing poe statements in the serial version. Comment out the lines 
that are marked 'serial' and uncomment the ones that are marked 'parallel' and
you should be alright though.

*** We note that the "aimschain" tool (also documented and supplied) is newer
    and may be more robust than the older NEB infrastructure supplied here. We 
    recommend the "aimschain" tool since it is actively supported.
