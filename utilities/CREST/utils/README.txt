
Interface script and utilities to run AIMS with CREST.
See file headers for more details and usage.
* Remember to update the absolute paths specified in the "definitions" in the
* beginning of CREST_run_AIMS.sh and aims2crest.sh !!!

Ofer Sinai, 2015
For more info, see:
O. Sinai, O. T. Hofmann, P. Rinke, M. Scheffler. G. Heimel and L. Kronik,
Phys. Rev. B (2015)


FILES:
================================
- CREST_run_AIMS.sh
  Main interface script, called via command specified in crest input file.
  Runs FHI-AIMS with given sheet charge, extracts results and saves in
  CREST-readable format.

- aims2crest.sh
  Called by CREST_run_AIMS.sh. Extracts FHI-AIMS output and saves in
  CREST-readable format.

- sevcleanup.py
  Called by aims2crest.sh. Cleans up the electrostatic energy file.

- crest_in_check4aims.py
  Utility that goes over the CREST input file and checks some info against the
  FHI-AIMS input. Useful to make sure that the data in the input files match
  prior to beginning a CREST calculation.
