
Installation instructions for CREST program Matlab standalone application and
data on files included.

Ofer Sinai, 2015
For more info, see:
O. Sinai, O. T. Hofmann, P. Rinke, M. Scheffler. G. Heimel and L. Kronik,
Phys. Rev. B (2015)


1. Installing the Matlab runtime environment
(Based on the automatic instructions written by the MATLAB Compiler)
================================

You must have the appropriate MATLAB Compiler Runtime (MCR) installed on
your machine. You need the 64-bit version 7.17 (R2012a).

If the MCR is not installed, it can be freely downloaded from the
MathWorks website:

   http://www.mathworks.com/products/compiler/mcr/

where instructions on installing it on different machines can also be found.

If you have MATLAB installed, you can instead enter
  
      >>mcrinstaller
      
at MATLAB prompt. This MCR Installer command displays the location of the
MCR Installer. You can then run the installer.



2. Updating the crest shell script, crest_run.sh:
================================
  crest_run.sh sets up the correct environment and runs CREST_main. See the
  file header for instructions on use. You must manually modify the variables
  MCRROOT and CREST_BIN in this file:

  . MCRROOT is to be the directory where version 7.17 of MCR is installed or
    the directory where MATLAB is installed on the machine.
  . CREST_BIN should contain the path to the CREST_main executable.



3. Running CREST
================================
Simply run:
crest_run.sh <crest_input_file> <crest_output_file>




Files/folders in this directory
================================
- crest_run.sh
  Shell script to setup the correct environment and run CREST_main.

- CREST_main
  The compiled executable, to be run using crest_run.sh.

- MCR_readme.txt
  Readme generated automatically by MCR upon compilation.

- crest.in.full
  Input file containing all possible commands that CREST recognizes with short
  descriptions of their meanings. The closest we have to documentation at the
  moment.

- src/
  Subdirectory containing all Matlab source code for the program.

- license.txt
  The standard BSD 2-clause license
  (see http://opensource.org/licenses/BSD-2-Clause)

- This readme file.
