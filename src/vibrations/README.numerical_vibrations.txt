How to use the scripts numerical_vibrations.pl and numerical_vibrations_jmol.pl
===============================================================================

Preparation:
(1) compile with the ../Makefile target vibrations or vibrations.mkl. This requires
    the same Makefile settings as any other AIMS compilation.
(2) in ../../bin/aims.vibrations.${VERSION}.pl or ../../bin/aims.vibrations.${VERSION}.mpi.pl set 
   $EXE       = 'aims executable including path' ;

   if you are running MPI calculations, also set 
   $EXE_CALL, $HESSIAN_DIAGONALIZATION_CALL, $HESSIAN_DIAGONALIZATION_SUFFIX
   in order to make sure that the correct syntax is being used. 
   Working examples are provided for the poe environment with a specified number of processors.

----------------------------------------------------------------------------------------------

Running:

requires 
control.in & geometry.in 
files.

change into the directory where these files are and execute the script 
aims.vibrations.${VERSION}.pl [jobname] [delta]
after the changes described above. If no jobname is specified, "basic" is assumed. 

----------------------------------------------------------------------------------------------

Results : 

${jobname}.xyz     - molecular file that contains the structure, the eigenmodes, their 
		     frequencies, and their IR intensity 
		     to see a summary of all the modes, type
	    	     grep -e cm ${jobname}.xyz
 
${jobname}.vib.out - output of the Hessian diagonalizer, also contains information on the moments of 
		     inertia and the total mass, for molecular thermodynamics calculations. 

*.dat              - various files required to rerun the Hessian diagonalization.

There are a lot of intermediate files, which are cleaned at the end of a single script's run. 
