#!/usr/bin/perl
############################################################################################
#
#  Begin user-adjusted variables
#
# THESE FOUR VARIABLES NEED CHECKING:
# the AIMS bin directory 

$AIMS_BINDIR = "AIMS_BINARY_DIRECTORY" ;

# aims executable 

$EXE = "AIMS_BINARY";

# calling sequence: any prefixes for serial/parallel execution.
# serial example is default

$CALLING_PREFIX = "";
$CALLING_SUFFIX = "";  

#
# example for mpirun with two cores
#
#$CALLING_PREFIX = "mpirun -np 2 ";
#$CALLING_SUFFIX = "";
#
# example for poe with 16 cores:
#
#$CALLING_PREFIX = "poe ";
#$CALLING_SUFFIX = "-procs 16";
#
# End user-adjusted variables
#
############################################################################################

print "\n";
print "===================================================\n\n";
print "        Entering FHI-aims Phonon calculation       \n";
print "                 Version XXXXXX                    \n\n";
print "===================================================\n\n";

# JOBNAME - is the input parameter, the default name is "phonon"
$JOBNAME = @ARGV[0];if (@ARGV[0] eq ''){$JOBNAME = "phonon" ;}

$STARTIN       =  join '','startindex.',$JOBNAME,'.dat';
$STARTOUT      =  join '','>startindex.',$JOBNAME,'.dat';

# setup executable calls
$HESSIAN_DIAGONALIZATION_EXE  = "DYNAMICMATRIXDIAGONALIZER";
$SYMMETRY_EXE                 = "SYMMETRYCHECKER";
$EXE_CALL                     = "$CALLING_PREFIX $AIMS_BINDIR/$EXE $CALLING_SUFFIX";
$HESSIAN_DIAGONALIZATION_CALL = "$CALLING_PREFIX $AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE $CALLING_SUFFIX";
$SYMMETRY_CALL                = "$CALLING_PREFIX $AIMS_BINDIR/$SYMMETRY_EXE $CALLING_SUFFIX";

# default values
$force_cutoff     = 0.000;
$delta_default    = 0.01;
@supercell        = (0,0,0);
$startindex       = 0; 
$working_dir      = "phonon_workdir";
$force_file       = "phonon_displacement_forces.dat";
$band_output_file = "phonon_band_structure.dat";
$temp_force_file  = "phonon_temporary_force_file.dat";

# check executable files
if (! -e "$AIMS_BINDIR/$EXE"){print "Can't find AIMS executable $AIMS_BINDIR/$EXE. \n Aborting. \n";die ;}
if (! -e "$AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE"){print "Can't find Hessian diagonalizer $AIMS_BINDIR/$HESSIAN_DIAGONALIZATION_EXE. \nAborting. \n"; die ;}
if (! -e "$AIMS_BINDIR/$SYMMETRY_EXE"){print "Can't find symmetry checking routine $AIMS_BINDIR/$SYMMETRY_EXE. \nAborting. \n"; die ;}

# make temporary directory for all the stuff to go in, if it does not already exist
if (! -e "$working_dir/"){print "Creating working directory $working_dir/ for all necessary files. \n\n"; mkdir $working_dir; }
else{print "Found working directory $working_dir/ . \n will use its contents for further calculations.\n\n";}
# check input geometry and control files
if (! -e "control.in" ){print "No suitable control.in file found. \nAborting.\n" ;die ;}
if (! -e "geometry.in"){print "No suitable geometry.in file found.\nAborting.\n" ;die ;}

chdir $working_dir;

# check force constant file and delete if it already exists - will be set up by the complete run anyway.
if (-e "$force_file"){ system "rm $force_file"; }

# parse control.in and geometry.in 
# read ../control.in file (i.e. original!) for necessary parameters.
$n_bands = 0;
@bands = ();
open(CONTROL, "../control.in");
open(OUTPUT, ">control.in");
print "Parsing original control.in : \n";
while ($_ = <CONTROL>)
{
    @line = split " ", $_ ;
    if (!/\#/) {}    
    if($line[0] eq "phonon") 
    {
	if ($line[1] eq "supercell")
	{
	    @supercell = ($line[2],$line[3],$line[4]);
	    print " | Chosen supercell is : ", $supercell[0]," x ",$supercell[1]," x ",$supercell[2],"\n";
	}
	if ($line[1] eq "displacement")    
	{
	    $delta    = $line[2]; 
	    print " | The displacement is: ", $delta,"  Angstrom \n";
	}
	if ($line[1] eq "frequency_unit")
	{
	    $frequency_unit = $line[2];
	    if (($frequency_unit!~/cm\^-1/)&&($frequency_unit!~/THz/))    # wrong units, give error message & stop
	    {
		print " * WARNING: wrong frequency units, allowed values are cm^-1 and THz! \n";
		print " *          please correct and restart. \n";
		die; 
	    }
	    print " | The frequency output will be in units of ",$frequency_unit,"\n";
	}
	if ($line[1] eq "band")
	{
	    push @bands,[($line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$line[8],$line[9],$line[10])]; 
	    $n_bands++;		
	}	
	# keep those pieces of the control.in that really matter for the calculation of the phonons, in particular, avoid MD and relaxations
	$relaxing_structure = 0;
	if ((/restart/)||(/MD/)){}else{if (/relax_geometry/){$relaxing_structure = 1};print OUTPUT $_;}
    }
}
# these are things THAT MUST BE in the output control file. Otherwise, the phonon aren't happening ... 
# print OUTPUT "\nrestart_write_only restart.dat\n" ;
print OUTPUT "compute_forces .true. \n";
print OUTPUT "final_forces_cleaned .true. \n";
close(OUTPUT);
close(CONTROL);
print "\n";

# check if force and delta threshold was set, check if supercell was reasonable
if ($delta == 0){$delta = $delta_default;print STDOUT "Using default value for finite difference step size, delta = ",$delta,"\n";}
for ($icell = 0; $icell < 3; $icell++ ){if ($supercell[$icell] < 1){print "Nonsensical supercell dimension $icell of $supercell[$icell] in control.in!\n Please correct.\n Aborting.\n"; die}}
$use_periodic = 0;
@lattice_in = ();
open(GEOMETRY, "../geometry.in");
print "Parsing initial geometry.in : \n";
while (<GEOMETRY>)
{ 
    @line = split " ", $_ ; 
    if ($line[0]=~/\#/){}
    if ($line[0] eq "atom")
    {
	# add coord triple to the end of array coords
	push @coords_in, [($line[1], $line[2], $line[3])] ;
	$species_in[$n_atoms_in] = $line[4] ;
	print " | Found atom $species_in[$n_atoms_in] at $line[1] $line[2] $line[3]\n";
	$n_atoms_in++ ;
    }
    if ($line[0] eq "lattice_vector")   
    {
	# system is PERIODIC
	push @lattice_in, [($line[1], $line[2], $line[3])] ;
	$use_periodic++;
	print " | Found lattice vector $line[1] $line[2] $line[3] \n";
    }
}
close(GEOMETRY);
print "\n";

# check if there is a resonable number of atoms and lattice vectors
if (($use_periodic==0)||($n_atoms_in==0)){print "Not enough atoms and/or lattice vectors.\nPlease correct.\nAborting.\n";die;}

# check, if there is a startindex (if some jobs were already done) and if so read it
if (-e $STARTIN){open (INDEX, $STARTIN);$startindex=<INDEX>;close(INDEX);} 

# build central point geometry & run calculation
# --> build complete supercell, that means one only has to change the atoms in the first supercell at the beginning of all other calculations

# skip calculation of central lattice ... the restart option does not actually bring any advantages
# open(GEOMETRY,">geometry.in");
@coords_sc   = ( );
$n_supercell =  0 ;
$i_atom_sc   =  0 ;
@species     = ( );
for ($isc1 = 0; $isc1 < $supercell[0]; $isc1++)
{
    for ($isc2 = 0; $isc2 < $supercell[1]; $isc2++)
    {
	for ($isc3 = 0; $isc3 < $supercell[2]; $isc3++)
	{
	    for ($iatom = 0; $iatom < $n_atoms_in; $iatom++)
	    {
		$x = $coords_in[$iatom][0] + $isc1*$lattice_in[0][0] + $isc2*$lattice_in[1][0] + $isc3*$lattice_in[2][0];
		$y = $coords_in[$iatom][1] + $isc1*$lattice_in[0][1] + $isc2*$lattice_in[1][1] + $isc3*$lattice_in[2][1];
		$z = $coords_in[$iatom][2] + $isc1*$lattice_in[0][2] + $isc2*$lattice_in[1][2] + $isc3*$lattice_in[2][2];
		push @coords_sc, [($x,$y,$z)];
#		printf GEOMETRY "atom %15.8E %15.8E %15.8E $species_in[$iatom]\n",$x,$y,$z;
		$species[$i_atom_sc] = $species_in[$iatom];
		$i_atom_sc++;
	    }
	    $n_supercell++;
	}	
    }
}

# print lattice vectors
# printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[0]*$lattice_in[0][0], $supercell[0]*$lattice_in[0][1], $supercell[0]*$lattice_in[0][2];
# printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[1]*$lattice_in[1][0], $supercell[1]*$lattice_in[1][1], $supercell[1]*$lattice_in[1][2];
# printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[2]*$lattice_in[2][0], $supercell[2]*$lattice_in[2][1], $supercell[2]*$lattice_in[2][2];
# close(GEOMETRY);
# total number of atoms in DFT calculations, across all supercells:
$n_atoms = $n_supercell*$n_atoms_in;

# first do a single calculation at the central point (no displacement)
# create control-file

# $identifier = 'central' ;
# $out_file   = join "","phonon.",$identifier,".out" ;

# jump over jobs that are already done - but only jump over the aims calculation, do everything else
# $counter = 0;
# if (++$counter > $startindex) 
# {
#     print "Running SCF calculation on central point. \n";
#     system "$EXE_CALL > $out_file < /dev/null" ;
# }

# now get coordinates
#open (RESULT,"$out_file") ;		
#$conv = 1. ;
#$found = 0. ;
#while (<RESULT>)
#{

#    # extract final coordinates from output file and make sure that they are converged ... 
#    if (/Final\ atomic\ structure/) 
#    {
#	$_ = <RESULT>;
#	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
#	{
#	    $_ = <RESULT>;
#	    @line = split " ", $_ ;
#	    for ($i_coords = 0; $i_coords <3; $i_coords++)
#	    {
#		$coords_sc[$i_atom][$i_coords] = @line[$i_coords+1] ;
#	    }
#	}
#    }
#    if (/Have\ a\ nice\ day/){$found = 1 ;}
#    if (/WARNING\!\ SELF\-CONSISTENCY\ CYCLE\ DID\ NOT\ CONVERGE/){$conv = 0 ;}
#}
#close (RESULT) ;
#if (($found eq 0) || ($conv eq 0))
#{
#    print STDOUT " * WARNING: ",$identifier, " Not converged?\n" ;
#    print STDOUT " * Something is seriously wrong here. \n * Aborting.";
#    die ;
#}

# WARNING: if this is ever reinstated, the $counter check in the symmetry calculation must be changed
# write next start index: one new calculation is done, write number to file 'startindex.dat'
#open (INDEX, $STARTOUT) ;
#print INDEX $counter ;
#close(INDEX) ;

print "\nBeginning finite displacement calculations. \n\n";

# create control file from control.in.$JOBNAME with only those pieces that we particularly would like to have
open(CONTROL, "../control.in");
open(OUTPUT, ">control.in");
while ($_ = <CONTROL>){if ((/relax_geometry/)||(/MD/)||(/restart/)){}else{print OUTPUT $_;}}
# print OUTPUT "\nrestart_read_only restart.dat\n" ;
print OUTPUT "compute_forces .true. \n";
print OUTPUT "final_forces_cleaned .true. \n";
close (CONTROL) ;
close (OUTPUT)  ;

$force_number = 0;
open(HESSIAN, ">hessian.$JOBNAME.dat");
for ($i_atom = 0; $i_atom < $n_atoms_in; $i_atom++)
{
    for ($i_coord = 0; $i_coord < 3; $i_coord++)
    {
	# here we only calculate the hessian matrix row by row, the entire thing is only used in the 
	# actual phonon diagonalization script. 
	@hessian = ();
	for ($prefactor = -1; $prefactor < 2; $prefactor += 2)
	{
	    $identifier = join "","i_atom_",$i_atom,".i_coord_",$i_coord,".displ_",$prefactor*$delta;

	    $symmetry_found = 0;
	    $force_number ++ ;            # separate count for number of forces, which is not necessarily the same as number of DFT results

	    # Have already done more than one displacement
	    if ($counter > 0)
	    {
		print "Checking for symmetry equivalence with a previously calculated displacement for job $identifier ... ";
		# DEBUG:
		# print "\n$SYMMETRY_CALL new_displacement $force_file ../geometry.in ../control.in $force_number > symmetry_check.$identifier.out\n";
	        # call symmetry check
		system "$SYMMETRY_CALL new_displacement $force_file ../geometry.in ../control.in $force_number > symmetry_check.$identifier.out";
		# check output for any sign of useful forces
		open(SYMMETRY,"<symmetry_check.$identifier.out");
		$_=<SYMMETRY>;   # read first line in symmetry output
		if (/Found symmetry operation to match previously calculated displacement/)   # yes, have all necessary forces
		{
		    $symmetry_found = 1;
		    <SYMMETRY>;<SYMMETRY>;<SYMMETRY>;<SYMMETRY>;   # read (and forget) transformation matrix
		    @line = split " ", <SYMMETRY>;                 # this line contains the forces in the correct order 
		    $coord_count = 0;
		    for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2 ++ )
		    {
			for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2 ++)
			{
			    $forces[$i_atom_2][$i_coord_2] = $line[$coord_count];
			    $coord_count++
			}
		    }
		}
		close(SYMMETRY);
	    }

	    # if symmetry found, do hessian addition, if not do DFT calculation
	    if ($symmetry_found)
	    {
		# do hessian addition for this particular displacement & add line to force file for use in future symmetry checks
		$index = 0;
		open(FORCEFILE, ">>$force_file");
		for ($i_atom_2 = 0; $i_atom_2 < $n_atoms;  $i_atom_2++)
		{
		    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		    {
			@hessian[$index] +=  -$prefactor*$forces[$i_atom_2][$i_coord_2] / (2*$delta);
			$index++;
			printf FORCEFILE "%20.15E  ",$forces[$i_atom_2][$i_coord_2];
		    } 
		}
		print FORCEFILE "\n";
		close(FORCEFILE);
		print "found. \n";
	    }
	    else
	    {
		# have not found symmetry equivalent displacement among previous calculations: 
		#    need to do DFT calculations for that displacement now

		if ($counter > 0) {print "not found. \n";}
		# build finitely displaced geometries and run calculations
		open(GEOMETRY,">geometry.in.$identifier");
		for ($i_atom_print = 0; $i_atom_print < $n_atoms; $i_atom_print++)
		{
		    if ($i_atom_print == $i_atom)   # this one has a finite displacement
		    {
			@thiscoord = ($coords_sc[$i_atom_print][0],$coords_sc[$i_atom_print][1],$coords_sc[$i_atom_print][2]);
			
			$thiscoord[$i_coord] += $prefactor*$delta;
			# change relevant coordinate
			print GEOMETRY "atom $thiscoord[0] $thiscoord[1] $thiscoord[2] $species_in[$i_atom_print]\n";
			
		    }
		    else # this one is printed as is
		    {
			print GEOMETRY "atom $coords_sc[$i_atom_print][0] $coords_sc[$i_atom_print][1] $coords_sc[$i_atom_print][2] $species[$i_atom_print]\n";
		    }
		}
		# add lattice vectors for good measure
		printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[0]*$lattice_in[0][0], $supercell[0]*$lattice_in[0][1], $supercell[0]*$lattice_in[0][2];
		printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[1]*$lattice_in[1][0], $supercell[1]*$lattice_in[1][1], $supercell[1]*$lattice_in[1][2];
		printf GEOMETRY "lattice_vector %15.8E %15.8E %15.8E \n", $supercell[2]*$lattice_in[2][0], $supercell[2]*$lattice_in[2][1], $supercell[2]*$lattice_in[2][2];
		close(GEOMETRY);
		# actually make a geometry.in 
		system "cp geometry.in.$identifier geometry.in";
		# run aims on this setup
		$out_file = join "","phonon.",$identifier,".out" ;
		if (++$counter > $startindex)
		{
		    print "Computing DFT forces for job $out_file. \n";
		    system "$EXE_CALL > $out_file < /dev/null" ;
		}
		else 
		{
		    print "Already have DFT results for job $out_file . \n";
		}
		# check if converged and obtain DFT forces
		@forces = ();
		open (RESULT,"$out_file") ;		
		$conv = 1. ;
		$found = 0. ;
		while (<RESULT>)
		{
		    # extract all forces from output file and make sure that they are converged ... 
		    if (/Total\ atomic\ forces/)
		    {
			for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2++)
			{
			    # read line and split it according to the spaces
			    $_ = <RESULT> ;
			    @line = split " ", $_ ;
			    # extract single coordinate forces
			    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
			    {
				$forces[$i_atom_2][$i_coord_2] = $line[$i_coord_2 + 2];
			    }
			}
		    }
		    if (/Have\ a\ nice\ day/){$found=1;}
		    if (/WARNING\!\ SELF\-CONSISTENCY\ CYCLE\ DID\ NOT\ CONVERGE/){$conv=0;}
		}
		close (RESULT) ;
		# stop if not converged
		if (($found==0)||($conv==0)){print " * WARNING: ",$identifier, " Not converged?\n* Please check this problem before continuing.\n* Aborting.\n";die;}
		
		# write next start index: one new calculation is done, write number to file 'startindex.dat'
		open (INDEX, $STARTOUT) ;
		print INDEX $counter ;
		close (INDEX) ;
		
		print "Checking DFT results for applicable symmetry operations ... ";

		# print forces to temporary (unsymmetrized) force file 
		open(TEMPFORCE, ">$temp_force_file");
		for ($i_atom_2 = 0; $i_atom_2 < $n_atoms;  $i_atom_2++)
		{
		    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		    {
			printf TEMPFORCE "%20.15E  ",$forces[$i_atom_2][$i_coord_2];
		    } 
		}
		print TEMPFORCE "\n";
		close(TEMPFORCE);
		
		# call symmetry check
		system "$SYMMETRY_CALL DFT_symmetrization $temp_force_file ../geometry.in ../control.in $force_number > symmetry_check.$identifier.out";
		# remove one-line file used to store forces
		system "rm $temp_force_file";

		# obtain symmetrized forces from output of symmetry check
		open(SYMMETRY,"<symmetry_check.$identifier.out");
		$_ = <SYMMETRY>;   # first line contains number of symmetries 
		if (/Found identity to be the only applicable/)
		{
		    $number_of_symmetries = "none";
		}
		else
		{
		    @line = split " ", $_;
		    $number_of_symmetries = $line[1];		    
		}
		<SYMMETRY>;
		@line = split " ", <SYMMETRY>;   # this line contains the forces in the correct order 
		$coord_count = 0;
		for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2 ++ )
		{
		    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2 ++)
		    {
			$forces[$i_atom_2][$i_coord_2] = $line[$coord_count];
			$coord_count++;
		    }
		}
		close(SYMMETRY);

		print "$number_of_symmetries found.\n\n";

		$index = 0;
		# add contribution of this calculation to hessian matrix
		# print out line to force file, for use in the symmetry calculation. 
		open(FORCEFILE, ">>$force_file");
		for ($i_atom_2 = 0; $i_atom_2 < $n_atoms;  $i_atom_2++)
		{
		    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
		    {
			# This is inside the loop over all translations of one single atom, 
			# add the relevant component of the finite difference to the (already-present) hessian
			# note that the hessian is the actual second partial derivative of the energy wrt the
			# coordinates, i.e. H = -df/dx; hence the sign.
			@hessian[$index] +=  -$prefactor*$forces[$i_atom_2][$i_coord_2] / (2*$delta);
			$index++;
			printf FORCEFILE "%20.15E  ",$forces[$i_atom_2][$i_coord_2];
		    } 
		}
		print FORCEFILE "\n";
		close(FORCEFILE);
	    }   # end if symmetry found
	}  # end computing one of the Hessian rows.
	
	# print Hessian row to Hessian file
	$index = 0;
	for ($i_atom_2 = 0; $i_atom_2 < $n_atoms; $i_atom_2++)
	{
	    for ($i_coord_2 = 0; $i_coord_2 < 3; $i_coord_2++)
	    {
		printf HESSIAN "%20.10E  ", $hessian[$index];
		$index++;
	    }
	}
	print HESSIAN "\n";
    }
}
close(HESSIAN);

# run new diagonalization and dos code
print "\n Calling Hessian diagonalization routine ... \n";
chdir "..";
# DEBUG:
# print "$HESSIAN_DIAGONALIZATION_CALL $working_dir/hessian.$JOBNAME.dat geometry.in control.in $working_dir | tee output.$JOBNAME.dynamic_matrix_diagonalization\n\n";
system "$HESSIAN_DIAGONALIZATION_CALL $working_dir/hessian.$JOBNAME.dat geometry.in control.in $working_dir | tee output.$JOBNAME.dynamic_matrix_diagonalization";

# if any bands were computed, then write a continuous band file 
if ($n_bands > 0) 
{
    print "\n Writing completed phonon dispersion relation to file $band_output_file \n";
    
    @k_vectors = ();
    @point_spacings = ();
    # compute reciprocal lattice
    push @k_vectors, [($lattice_in[1][1]*$lattice_in[2][2]-$lattice_in[1][2]*$lattice_in[2][1],$lattice_in[1][2]*$lattice_in[2][0]-$lattice_in[1][0]*$lattice_in[2][2],$lattice_in[1][0]*$lattice_in[2][1]-$lattice_in[1][1]*$lattice_in[2][0])];
    push @k_vectors, [($lattice_in[2][1]*$lattice_in[0][2]-$lattice_in[2][2]*$lattice_in[0][1],$lattice_in[2][2]*$lattice_in[0][0]-$lattice_in[2][0]*$lattice_in[0][2],$lattice_in[2][0]*$lattice_in[0][1]-$lattice_in[2][1]*$lattice_in[0][0])];
    push @k_vectors, [($lattice_in[0][1]*$lattice_in[1][2]-$lattice_in[0][2]*$lattice_in[1][1],$lattice_in[0][2]*$lattice_in[1][0]-$lattice_in[0][0]*$lattice_in[1][2],$lattice_in[0][0]*$lattice_in[1][1]-$lattice_in[0][1]*$lattice_in[1][0])];
    $prefactor = 6.28318531/($lattice_in[0][0]*$k_vectors[0][0]+$lattice_in[0][1]*$k_vectors[0][1]+$lattice_in[0][2]*$k_vectors[0][2]);

    open(BANDOUTPUT,">$band_output_file");
    print BANDOUTPUT "#\n# plotting phonon bands from FHI-aims!\n#\n";
    for ($i=0;$i<3;$i++) 
    {
	for ($j=0;$j<3;$j++){$k_vectors[$i][$j] *= $prefactor;}
	printf BANDOUTPUT "# found reciprocal lattice vector %10.6f %10.6f %10.6f \n", $k_vectors[$i][0], $k_vectors[$i][1], $k_vectors[$i][2];
    }
    print BANDOUTPUT "#\n";
    $distance = 0.0;
    for ($i=0;$i<$n_bands;$i++)
    {
	$points     = $bands[$i][6];
	@start = ($k_vectors[0][0]*$bands[$i][0]+$k_vectors[1][0]*$bands[$i][1]+$k_vectors[2][0]*$bands[$i][2],
		  $k_vectors[0][1]*$bands[$i][0]+$k_vectors[1][1]*$bands[$i][1]+$k_vectors[2][1]*$bands[$i][2],
		  $k_vectors[0][2]*$bands[$i][0]+$k_vectors[1][2]*$bands[$i][1]+$k_vectors[2][2]*$bands[$i][2]);
	@end   = ($k_vectors[0][0]*$bands[$i][3]+$k_vectors[1][0]*$bands[$i][4]+$k_vectors[2][0]*$bands[$i][5],
		  $k_vectors[0][1]*$bands[$i][3]+$k_vectors[1][1]*$bands[$i][4]+$k_vectors[2][1]*$bands[$i][5],
		  $k_vectors[0][2]*$bands[$i][3]+$k_vectors[1][2]*$bands[$i][4]+$k_vectors[2][2]*$bands[$i][5]);
	$dist  = sqrt(($start[0]-$end[0])*($start[0]-$end[0])+($start[1]-$end[1])*($start[1]-$end[1])+($start[2]-$end[2])*($start[2]-$end[2]));
	$spacings = $dist/($points-1); 
	@point_spacings =  (@point_spacings,$spacings);
	printf BANDOUTPUT "# Starting point for band %02d, %6s = (%9.5f,%9.5f,%9.5f) will be at real distance = %11.5f \n", $i,$bands[$i][7],$start[0],$start[1],$start[2],$distance; 
	$distance += $dist;
	printf BANDOUTPUT "# Ending   point for band %02d, %6s = (%9.5f,%9.5f,%9.5f) will be at real distance = %11.5f \n#\n", $i,$bands[$i][8],$end[0],$end[1],$end[2],$distance; 
    }
    $distance = 0.0;
    for ($i=0;$i<$n_bands;$i++)
    {
	$points = $bands[$i][6];
	if ($i<10) {$file = join '', 'phonon_band100',$i+1,'.out';}
	else {if ($i<100) {$file = join '', 'phonon_band10',$i,'.out';}
	      else {$file = join '', 'phonon_band1',$i,'.out';}}
	if (-e "phonon_workdir/$file"){$file="phonon_workdir/$file";}
	open(BANDFILE,$file);
	$distance -= $point_spacings[$i];
	for ($j=0;$j<$points;$j++)
	{
	    $distance += $point_spacings[$i];
	    @line = split " ", <BANDFILE>;
	    $nbands = (@line-4)/2;
	    print BANDOUTPUT "$distance ";for ($k=0;$k<$nbands;$k++){$energy = $line[5+2*$k]; print BANDOUTPUT "$energy ";}print BANDOUTPUT "\n"; 
	}
	close(BANDFILE);
    }
    close(BANDOUTPUT);
    print "\n";
}

# all done, output quitting statement
print "\n Leaving FHI-aims phonon script. \n\n ==========================================================\n\n";
