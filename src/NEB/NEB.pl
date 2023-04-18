#!/usr/bin/perl

# set defaults for values that are absolutely necessary for running
$force_convergence_criterion = 1e-2;
$n_max_iteration             = 100;
$n_iteration_start           = 0;
$n_atoms_read                = 0;
$n_images_read               = 0;
$aims_geometry_file_template = "geometry.in.template";
$aims_control_file_template  = "control.in.template";

# variables for calling DFT code - in series:
 $EXE_PREFIX = "";
 $EXE_SUFFIX = "";
# variable for calling DFT code - in parallel; CHECK before using 
#$EXE_PREFIX = "poe ";
#$EXE_SUFFIX = "";

# set default control file, but read from input if it does exist!
$control_file = 'control_NEB.in';

# read input from control file specified in the command line, if it exists
if (@ARGV[0]) {$control_file = @ARGV[0];}

# check existence of control file, complain and stop if it ain't there:
if (!-e $control_file){ print "Control file $control_file does not exist!\nAborting execution."; exit;}

# read relevant data from control file
open(CONTROL, "$control_file");
while ($_ = <CONTROL>)
{
    if (!/\#/){    # only execute if no comment was involved!
	if (/dft_exe/)                       {@line = split " ", $_; $dft_exe                       = $line[1]; print "DFT executable:                ", $dft_exe,                       "\n";}
	if (/input_template/)                {@line = split " ", $_; $input_template                = $line[1]; print "xyz template:                  ", $input_template,                "\n";}
	if (/aims_control_file_template/)    {@line = split " ", $_; $aims_control_file_template    = $line[1]; print "AIMS input template:           ", $aims_control_file_template,    "\n";}
	if (/aims_geometry_file_template/)   {@line = split " ", $_; $aims_geometry_file_template   = $line[1]; print "AIMS geometry file template:   ", $aims_geometry_file_template,   "\n";}
	if (/aims_restart_string/)           {@line = split " ", $_; $aims_restart_string           = $line[1]; print "AIMS restart string:           ", $aims_restart_string,           "\n";}
	if (/n_iteration_start/)             {@line = split " ", $_; $n_iteration_start             = $line[1]; print "n_iteration_start:             ", $n_iteration_start,             "\n";}
	if (/n_max_iteration/)               {@line = split " ", $_; $n_max_iteration               = $line[1]; print "n_max_iteration:               ", $n_max_iteration,               "\n";}
	if (/force_convergence_criterion/)   {@line = split " ", $_; $force_convergence_criterion   = $line[1]; print "force_convergence_criterion:   ", $force_convergence_criterion,   "\n";}
	if (/start_BFGS/)                    {@line = split " ", $_; $start_BFGS                    = $line[1]; print "start_BFGS:                    ", $start_BFGS,                    "\n";}
	if (/save_data/)                     {@line = split " ", $_; $save_data                     = $line[1]; print "save_data:                     ", $save_data,                     "\n";}    
	if (/dft_code/)                      {@line = split " ", $_; $dft_code                      = $line[1]; print "DFT code:                      ", $dft_code,                      "\n";}
	if (/castep_control_file_template/)  {@line = split " ", $_; $castep_control_file_template  = $line[1]; print "castep control file template:  ", $castep_control_file_template,  "\n";}
	if (/castep_geometry_file_template/) {@line = split " ", $_; $castep_geometry_file_template = $line[1]; print "castep geometry file template: ", $castep_geometry_file_template, "\n";}
	if (/castep_restart_string/)         {@line = split " ", $_; $castep_restart_string         = $line[1]; print "castep restart string:         ", $castep_restart_string,         "\n";}
	if (/method/)                        {@line = split " ", $_; $method                        = $line[1]; print "method:                        ", $method,                        "\n";}
	if (/n_atoms/)                       {@line = split " ", $_; $n_atoms_read                  = $line[1];}
	if (/n_images/)                      {@line = split " ", $_; $n_images_read                 = $line[1];}
    }
}
close(CONTROL);
if (($dft_code!~'CASTEP')&&($dft_code!~'Castep')&&($dft_code!~'castep')){$dft_code='FHI-aims';}
$n_iteration = $n_iteration_start;
if (($n_iteration_start==0)&&(-e $save_data)){print "Someone has left a NEB save data file lying around: $save_data\nI am sure this was unintended, will delete this file before proceeding ... \n";system "rm -f $save_data";}

# check existence of necessary files and complain if they're not there, then stop ...
# dft_exe:
if (!-e $dft_exe){ print "Can't find DFT executable $dft_exe ...\nAborting execution.\n"; exit;}
if ($dft_code=~/FHI-aims/)
{
# aims_control_file_template exists:
     if (!-e $aims_control_file_template){ print "Can't find DFT control file template $aims_control_file_template ... \nAborting execution.\n"; exit;}
     else{&AIMS_check_periodic($aims_control_file_template);}
     if ((!-e $aims_geometry_file_template)&&($periodic)){print"Can't find DFT geometry template for periodic lattice vectors! \nPlease investigate.\n Aborting execution"}
}
else
{
    # Castep control/geometry templates exist?
    if (!-e $castep_control_file_template) { print "Can't find Castep parameter file template $castep_control_file_template ... \nAborting execution.\n"; exit;}
    if (!-e $castep_geometry_file_template){ print "Can't find Castep geometry file template $castep_control_file_template ... \nAborting execution.\n"; exit;}
}
# first iteration of the geometry:
$first_geometry_file = join '',$input_template,$n_iteration,'.xyz';
if (!-e $first_geometry_file){print "Could not find the initial NEB guess $first_geometry_file ... \nAborting execution.\n"; exit;}

# initial force convergence criterion should be larger than the maximal one ... 
$max_force_NEB = $force_convergence_criterion * 10;
# do iterations until either the maximum number is reached OR the forces have converged. 
while (($max_force_NEB > $force_convergence_criterion) && ($n_iteration < $n_max_iteration))
{

    print "opening file $input_template$n_iteration.xyz\n";
    open(INPUT, "$input_template$n_iteration.xyz");
    $n_images = 0;
    # read number of atoms, first entry in jmol movie ... 
    $n_atoms = <INPUT>;
    
    print "n_atoms = ${n_atoms}\n";

    # loop through images, copy geometry information to file ... 
    if ($dft_code=~/FHI-aims/)
    {
	while ($_ = <INPUT>)
	{
	    if (/START/){if (/CI-NEB/) {$method = "CI-NEB"}&print_geometry_file_AIMS("geometry.in.starting_configuration");}
	    if (/END/){&print_geometry_file_AIMS("geometry.in.final_configuration");}
	    if (/image/){$n_images++;&print_geometry_file_AIMS("geometry.in.image$n_images");}
	}
    } 
    else
    {
	while ($_ = <INPUT>)
	{
	    if (/START/){             system "cat $castep_geometry_file_template > $input_template.starting_configuration.cell";              open(GEOMETRY,">>$input_template.starting_configuration.cell");              &print_geometry_file_CASTEP;close(GEOMETRY);}
	    if (/END/)  {             system "cat $castep_geometry_file_template > $input_template.final_configuration.cell";                 open(GEOMETRY,">>$input_template.final_configuration.cell"   );              &print_geometry_file_CASTEP;close(GEOMETRY);}
	    if (/image/){$n_images++; system "cat $castep_geometry_file_template > $input_template.image$n_images.iteration$n_iteration.cell";open(GEOMETRY,">>$input_template.image$n_images.iteration$n_iteration.cell");&print_geometry_file_CASTEP;close(GEOMETRY);}
	}
    }
    close(INPUT);    

    # DFT calculations
    # starting point if necessary
    if ($n_iteration==0) 
    {
	print "DFT calculation of initial configuration ... ";
	if ($dft_code=~/FHI-aims/)
	{
	    system "cp geometry.in.starting_configuration geometry.in";
	    system "cat $aims_control_file_template | sed 's/$aims_restart_string/restart\.starting_configuration\.dat/g' > control.in";
	    system "echo 'compute_forces .true.' >> control.in";      
	    system "echo 'final_forces_cleaned .true.' >> control.in";
	    system "$EXE_PREFIX $dft_exe $EXE_SUFFIX > output.starting_configuration";
	    &check_AIMS_convergence('output.starting_configuration');
	}    
	else 
	{
	    system "cat $castep_control_file_template | sed 's/$castep_restart_string/castep_restart\.starting_configuration.dat/g' > $input_template.starting_configuration.param";
	    system "$EXE_PREFIX $dft_exe $input_template.starting_configuration $EXE_SUFFIX";
	    &check_CASTEP_convergence("$input_template.starting_configuration.castep");
	} 
    } 

    # all intermediate images
    for ($i_images = 1; $i_images <= $n_images; $i_images++)
    {
	# create proper restart files in the first iteration: copy from previous image ... 

	# FIXME: this does not work for periodic AIMS calculations at the moment, 
        # NEITHER does it work for CASTEP calculations!
	if ($n_iteration==0)
	{
	    if ($dft_code=~/FHI-aims/)
	    {
		if ($i_images == 1){&copy_AIMS_restart_files("restart.starting_configuration.dat","restart.image1.dat",$n_periodic, $n_tasks);}
		else
		{
		    $i_images_old = $i_images - 1;
		    &copy_AIMS_restart_files("restart.image$i_images_old.dat","restart.image$i_images.dat",$n_periodic,$n_tasks);
		}
		system "cat $aims_control_file_template | sed 's/$aims_restart_string/restart\.image$i_images\.dat/g' > control.in.image$i_images";
		system "echo 'compute_forces .true.' >> control.in.image$i_images";      
		system "echo 'final_forces_cleaned .true.' >> control.in.image$i_images";
	    }
	    else
	    {
		if ($i_images == 1)
		{
		    system "cp castep_restart.starting_configuration.dat castep_restart.image1.dat";
		}
		else
		{
		    $i_images_old = $i_images - 1;
		    system "cp castep_restart.image$i_images_old.dat castep_restart.image$i_images.dat";
		}
		system "cat $castep_control_file_template | sed 's/$castep_restart_string/castep_restart\.image$i_images.dat/g' > $input_template.image$i_images.iteration$n_iteration.param";
	    }
	}
	else 
	{
	    # create a param file for Castep information - simply move over old file ... 
	    if (($dft_code!~/FHI-aims/)&&($n_iteration>0))
	    {
		$n_last_iteration = $n_iteration - 1;
		system "mv $input_template.image$i_images.iteration$n_last_iteration.param $input_template.image$i_images.iteration$n_iteration.param";
	    }
	}

	print "DFT calculation of image $i_images ... ";
	if ($dft_code=~/FHI-aims/)
	{ system "cp geometry.in.image$i_images geometry.in"; system "cp control.in.image$i_images control.in"; system "$EXE_PREFIX $dft_exe $EXE_SUFFIX > output.image$i_images.iteration$n_iteration"; &check_AIMS_convergence("output.image$i_images.iteration$n_iteration");}
	else {system "$EXE_PREFIX $dft_exe $input_template.image$i_images.iteration$n_iteration $EXE_SUFFIX"; &check_CASTEP_convergence("$input_template.image$i_images.iteration$n_iteration.castep");}
    }

    # ending point, again only if necessary
    if ($n_iteration==0) 
    { 
 	print  "DFT calculation of final configuration ... ";
	if ($dft_code=~/FHI-aims/)
	{
	    &copy_AIMS_restart_files("restart.image$n_images.dat","restart.final_configuration.dat",$n_periodic,$n_tasks);
	    system "cat $aims_control_file_template | sed 's/$aims_restart_string/restart\.final_configuration\.dat/g' > control.in";
	    system "echo 'compute_forces .true.' >> control.in.final_configuration";      
	    system "echo 'final_forces_cleaned .true.' >> control.in.final_configuration";
	    system "cp geometry.in.final_configuration geometry.in"; 
	    system "$EXE_PREFIX $dft_exe $EXE_SUFFIX > output.final_configuration";
	    &check_AIMS_convergence('output.final_configuration'); 
	}
	else 
	{
	    system "cp castep_restart.image$n_images.dat castep_restart.final_configuration";
	    system "cat $castep_control_file_template | sed 's/$castep_restart_string/castep_restart\.final_configuration.dat/g' > $input_template.final_configuration.param";
	    system "$EXE_PREFIX $dft_exe $input_template.final_configuration $EXE_SUFFIX";
	    &check_CASTEP_convergence("$input_template.final_configuration.castep");
	}
    } 

    # make new movie with energies and all of the forces ... 
    print "Writing DFT results to jmol compatible file $input_template$n_iteration.with_forces.xyz \n";
    open(MOVIE,">$input_template$n_iteration.with_forces.xyz");
    print "   ... starting configuration \n";
    if ($dft_code=~/FHI-aims/)
    {
	&read_output_information_AIMS('output.starting_configuration');
    }
    else
    {
	&read_forces_CASTEP("$input_template.starting_configuration.castep");
	&read_geometry_CASTEP("$input_template.starting_configuration.cell");
    }
    $i_images = 0;
    &print_movie_iteration;
    for ($i_images = 1; $i_images <= $n_images; $i_images++)
    {
	print "   ... image $i_images \n";
	if ($dft_code=~/FHI-aims/)
	{
	    &read_output_information_AIMS("output.image$i_images.iteration$n_iteration");
	    system "rm -f output.image$i_images.iteration$n_iteration";
	}
	else
	{
	    &read_forces_CASTEP("$input_template.image$i_images.iteration$n_iteration.castep");
	    &read_geometry_CASTEP("$input_template.image$i_images.iteration$n_iteration.cell");
	    system "rm -f $input_template.image$i_images.iteration$n_iteration.castep";
	    system "rm -f $input_template.image$i_images.iteration$n_iteration.cell";
	    system "rm -f $input_template.image$i_images.iteration$n_iteration.check";
	    system "rm -f $input_template.image$i_images.iteration$n_iteration.bands";
	    system "rm -f $input_template.image$i_images.iteration$n_iteration.cst_esp";
	}
	&print_movie_iteration;
    }

    print "   ... final configuration \n";
    if ($dft_code=~/FHI-aims/)
    {
	&read_output_information_AIMS('output.final_configuration');
    }
    else
    {
	&read_forces_CASTEP("$input_template.final_configuration.castep");
	&read_geometry_CASTEP("$input_template.final_configuration.cell");
    }
    $i_images = $n_images+1;
    &print_movie_iteration;
    close(MOVIE);
    
    # make sure the NEB control file contains appropriate values of n_atoms and n_images:
    if (($n_atoms_read!=$n_atoms)||($n_images_read!=$n_images))
    {
	print "Control file contains inappropriate information on the number of atoms and images \nSetting proper values now.\n";
	system "mv $control_file $control_file.temp";
	open(CONTROL_BAD,"$control_file.temp");
	open(CONTROL_GOOD,">$control_file");
	while ($_ = <CONTROL_BAD>)
	{ 
	    if ($_=~/n_atoms/) {print CONTROL_GOOD "n_atoms   $n_atoms \n";}
	    else {if ($_=~/n_images/) {print CONTROL_GOOD "n_images   $n_images \n";}
		  else {print CONTROL_GOOD $_;}}
	}
	close(CONTROL_GOOD);
	close(CONTROL_BAD);
	system "rm $control_file.temp";
    }


    $inputname = "$input_template$n_iteration.with_forces.xyz";
    $n_iteration++;
    $outputname = "$input_template$n_iteration.xyz";
    print  "Calling NEB solver: \n";
    print  "$EXE_PREFIX ./NEB.x $control_file $inputname $outputname 2>&1 $EXE_SUFFIX\n\n";
    system "$EXE_PREFIX ./NEB.x $control_file $inputname $outputname 2>&1 $EXE_SUFFIX";

    # check for convergence - criterion was written by last NEB solver call
    open(INPUT, "$input_template$n_iteration.xyz");
    while ($_ = <INPUT>){if (/Maximum force component/){@line = split " ", $_;$max_force_NEB = $line[9];}}
    close(INPUT);
}
# END main program

#######################################################################################
#  sub print_geometry_file_AIMS
#  requires an open file handles GEOMETRY and INPUT
#  where INPUT points to an appropriate place in a jmol movie 
#  from which to extract the geometry
sub print_geometry_file_AIMS
{
    $geometry_name = $_[0];
    open(GEOMETRY,">$geometry_name");
    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
    {
	@line = split " ", <INPUT>;
	print  GEOMETRY "atom ",@line[1]," ",@line[2]," ",@line[3]," ",@line[0],"\n";
    }
    close(GEOMETRY);
    if (-e $aims_geometry_file_template)
    {
	print "cat $aims_geometry_file_template >> $geometry_name\n";
	system "cat $aims_geometry_file_template >> $geometry_name";
    }
}

#######################################################################################
#  sub check_AIMS_convergence
#  checks whether or not the file with name in the argument contains converged AIMS 
#  output information. Returns a variable $converged which contains a status-describing string. 
#  also sets the variables $n_tasks and $n_periodic
sub check_AIMS_convergence
{ 
    $converged_name = $_[0];
    $converged = "Output file not found";
    $n_tasks = 1;
    $n_periodic = 0;
    open(CONV_CHECK,"$converged_name");
    while ($_ = <CONV_CHECK>)
    {
	if (/Invoking FHI-aims/){$converged = "control.in or geometry.in not found";}
	if (/Begin self-consistency loop/){$converged = "has converged";}
	if (/Number of self-consistency cycles/){@line   = split " ", $_; $cycles = $line[6];}
	if (/SELF-CONSISTENCY CYCLE DID NOT CONVERGE/){$converged = "has NOT converged";}
	if (/parallel tasks\./){@line = split " ", $_; $n_tasks = $line[1]}
	if (/Unit cell/){$n_periodic = 3;}	   
    }
    close(CONV_CHECK);
    if ($converged=~/converged/){print "$converged in $cycles iterations \n"; if ($converged=~/NOT/){print "There is something wrong here, aborting. \n";exit;}}
    else{print "$converged.\n"; print "There is something wrong here, aborting!\n"; exit;}
}

#######################################################################################
#   sub AIMS_check_periodic
#   subroutine to check whether or not a supplied control file requires lattice vectors. 
#   if it does, it might be wise to add them to the geometry.in
sub AIMS_check_periodic
{
    $periodic = 0;
    open(CHECK_PERIODIC,$_[0]);
    while ($_ = <CHECK_PERIODIC>)
    {
	if (/k_grid/)           {$n_periodic = 3;}	   
	if (/k_points_external/){$n_periodic = 3;}
    }
    close(CHECK_PERIODIC)
}

#######################################################################################
#   sub copy_AIMS_restart_files
#   subroutine to copy AIMS checkpoint files - making sure that there are the right 
#   number of files for periodic and/or multi-processor calculations ... 
#
#   required arguments: input_base_name, output_base_name, $n_periodic, $n_tasks
sub copy_AIMS_restart_files
{
    $start_name = $_[0];
    $end_name   = $_[1];
    $periodic   = $_[2];
    $tasks      = $_[3];
    if (($periodic)&&($tasks>1))
    {
	$n_start = $tasks;
	if ($n_start>9) {$n_start = 9;}
	for ($i_start = 1; $i_start < $n_start; $i_start++){system "cp $start_name00$i_start $end_name00$i_start";}
	if ($n_tasks>9)
	{
	    $n_start=$n_tasks;
	    if ($n_start>99) {$n_start = 99;}
	    for ($i_start=10;$i_start<$n_start;$i_start++){system "cp $start_name0$i_start $end_name0$i_start";}	    
	}
	if ($n_tasks>99)
	{
	    $n_start=$n_tasks;
	    for ($i_start=100;$i_start<$n_start;$i_start++){system "cp $start_name$i_start $end_name$i_start";}	    
	}
    }
    else { system "cp $start_name $end_name"; }
}

#######################################################################################
#  sub read_output_information_AIMS
#  requires an open file handle OUTPUT, just the beginning.
#  obtains geometry, total energy and forces from a converged scf run
sub read_output_information_AIMS
{
    $converged_name = $_[0];
    $converged = "Output file not found";
    open(OUTPUT,"$converged_name");
    @geometry = ();
    @forces = ();
    $energy = 0;
    while ($_ = <OUTPUT>)
    {
	if (/Input geometry/)
	{
	    $_ = <OUTPUT>;
	    if ($_=~/\|\ Unit\ cell\:/) 
	    {
		<OUTPUT>;
		<OUTPUT>;
		<OUTPUT>;
	    }
	    <OUTPUT>;
	    <OUTPUT>;
	    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	    {
		@line = split " ", <OUTPUT>;
		push @geometry, [($line[4],$line[5],$line[6],@line[3])];
	    }
	}
	if (/Energy and forces in a compact form/)
	{
	    @line = split " ", <OUTPUT>;
	    @line = split " ", <OUTPUT>;
	    $energy = $line[5];
	}
	if (/Total atomic forces/)
	{
	    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	    {
		@line = split " ", <OUTPUT>;
		push @forces, [($line[2],$line[3],$line[4])];
	    }	    
	}
    }
    close(OUTPUT);
}

#######################################################################################
#  sub print_movie_iteration
#  writes jmol-compatible movie iteration from known geometry, energy and 
#  forces to file handle MOVIE
sub print_movie_iteration
{
    print MOVIE "$n_atoms";
    if ($i_images == 0){printf MOVIE "START - TOTAL ENERGY = %30.10f method: %s\n", $energy, $method;}
    else{if ($i_images == $n_images + 1){printf MOVIE "END - TOTAL ENERGY = %30.10f \n", $energy;}
    else{printf MOVIE "image %8i - ENERGY = %30.10f \n", $i_images, $energy;}}
    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
    {printf MOVIE "%s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n",$geometry[$i_atom][3],$geometry[$i_atom][0],$geometry[$i_atom][1],$geometry[$i_atom][2],$forces[$i_atom][0],$forces[$i_atom][1],$forces[$i_atom][2];}
}

sub print_geometry_file_CASTEP{print GEOMETRY "%BLOCK POSITIONS_ABS\n";for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++){@line = split " ", <INPUT>;print  GEOMETRY @line[0]," ",@line[1]," ",@line[2]," ",@line[3],"\n";}print GEOMETRY "%ENDBLOCK POSITIONS_ABS";}
sub read_forces_CASTEP {$forces_name = $_[0]; open(FORCES,"$forces_name");@forces=(); while ($_ = <FORCES>){if (/Symmetrised Forces/){@line = split " ", <FORCES>;@line = split " ", <FORCES>;@line = split " ", <FORCES>;@line = split " ", <FORCES>;@line = split " ", <FORCES>;for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++){@line = split " ", <FORCES>;push @forces, [(@line[3],@line[4],@line[5])];}}if (/Final free energy/){@line = split " ", $_;$energy = @line[5];}}close(FORCES);}
sub read_geometry_CASTEP{$geometry_name=$_[0];open(GEOMETRY,"$geometry_name");@geometry=(); while ($_ = <GEOMETRY>){if (/BLOCK POSITIONS_ABS/){for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++){@line = split " ", <GEOMETRY>;push @geometry, [(@line[1],@line[2],@line[3],@line[0])];}}}close(GEOMETRY);}
sub check_CASTEP_convergence{$converged_name = $_[0];$converged = "Output file not found";$cycles = -5;open(CONV_CHECK,"$converged_name");while ($_ = <CONV_CHECK>){if (/Symmetrised Forces/){$converged = "has converged";}if (/First principles methods using CASTEP/){$converged = "has NOT converged";}if (/--\ SCF/){$cycles = $cycles + 1;}}close(CONV_CHECK);$cycles = $cycles - 1;if ($converged=~/converged/){print "$converged in $cycles iterations. \n"; if ($converged=~/NOT/){print "There is something wrong here, aborting. \n";exit;}}else{print "$converged.\n"; print "There is something wrong here, aborting!\n"; exit;}}
