#! /usr/bin/perl
#
# RG 2006: small script to do a bh-optimization
#
# required sources:
# - this wrapper script
# - effernir
# - aims
#
# input files:
# - geometry.in.basic (structure to optimize)
# - control.in (control-file for aims)
# - control.in.bh (control-file for effernir)
# 
# created files:
# - dft_data.dat (all the so far calculated energies and forces
#   by localorb)
# - startindex.dat (contains just number of already calculated energies and forces)
#                  (read by get_cg_geometry)
# - temp.out (temporary logfile for localorb)
#
# temporary files are deleted afterwards

#$homedir = @ARGV[0] ;
$EXE_AIMS = @ARGV[0] ;
$EXE_BH = @ARGV[1] ;

@line = split "/", $EXE_AIMS ;
$exe_aims_local = @line[$#line];

@line = split "/", $EXE_BH ;
$exe_bh_local = @line[$#line];

# spin-polarized?
$spin_polarized = 0 ;

# file declarations
$out_file = "temp.out" ;
$bh_log_file = "bh_log.out" ;
$log_file = "log.out" ;
$dft_file = "dft_data.dat" ;

$STARTIN = "startindex.dat";
$STARTOUT = ">startindex.dat";

# seed for RNG
$seed_uni   = 93 ;
$seed_poiss = 16 ;

$startindex = 0 ;
$finished = 0 ;
$status = 0 ;

$n_atoms = 0 ;

# read startindex file
if (-e $STARTIN)
{
    open (INDEX, $STARTIN) ;
    $_ = <INDEX> ;
    @line = split " ", $_ ;
    $startindex = @line[0];
    $seed_uni = @line[1];
    $seed_poiss = @line[2];
    $relaxation_steps = @line[3];
    $self_consistency_cycles = @line[4];
    close (INDEX) ;
} else {
    open (INDEX, $STARTOUT) ;
    print INDEX "0 ", $seed_uni, " ", $seed_poiss, " 0  0\n" ;
    close (INDEX) ;
} 

system "rm -f $log_file" ;
system "rm -f $bh_log_file" ;

# first read geometry-file to get n_atoms
open (INPUT, "geometry.in.basic") ;
while (<INPUT>)
{ 
    if (/atom/)
    {
	$n_atoms++ ;
    }
}
close (INPUT) ;

until ($finished eq 1) {

    system "./$exe_bh_local > $out_file 2>&1" ;
    system "cat $out_file >> $bh_log_file" ;
#    system "cp $bh_log_file $homedir" ;

    if (-e "finished.tmp") {
	$finished = 1 ;
    }
    if ($finished eq 0) {

	$enough = 0 ;
	$abort = 0 ;

	until ($enough eq 1) {

	    system "poe ./$exe_aims_local > $out_file 2>&1" ;
	    system "cat $out_file >> $log_file " ;
    
	    $movie_file_name = join "", ">relax.", $startindex+1,".irc" ;
	    $restart_file_name = join "", "restart.", $startindex+1,".dat" ;

	    # check convergence
	    @line = () ;
	    $found = 0 ;
	    $conv  = 1 ;
	    $status = 0 ;
	    $updated_structure = 0 ;
	    $i_step = 0 ;
	    open (RESULT, "$out_file") ;
	    while (<RESULT>)
	    {
		if (/Total\ energy\ uncorrected/)
		{
		    @line = split " ", $_ ;
		    $energy[$i_step] = @line[5] ;
		}
		if (/Maximum\ force\ component/)
		{
		    @line = split " ", $_ ;
		    $max_force[$i_step] = @line[4] ;
		}
		if (/Total\ atomic\ forces/)
		{
		    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		    {
			$_ = <RESULT>;
			@line = split " ", $_ ;
			for ($i_coords = 1; $i_coords <=3; $i_coords++)
			{
			    $forces[$i_atom][$i_coords] = @line[$i_coords + 1] ;
			}
		    }
		}
		if (/Final\ atomic\ structure/) 
		{
		    $_ = <RESULT>;
		    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		    {
			$_ = <RESULT>;
			@line = split " ", $_ ;
			$species[$i_step][$i_atom] = $line[4] ;
			for ($i_coords = 1; $i_coords <=3; $i_coords++)
			{
			    $coords[$i_atom][$i_coords] = @line[$i_coords] ;
			    $coords_relax[$i_step][$i_atom][$i_coords] = @line[$i_coords] ;
			}
		    }
		    $i_step++ ;
		}
		if (/spin\ analysis/)
		{
		    $_ = <RESULT>;
		    $_ = <RESULT>;
		    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		    {
			$_ = <RESULT>;
			@line = split " ", $_ ;
			$spin_per_atom[$i_atom] = @line[2] ;
		    }
		}
		if (/Updated\ atomic\ structure/)
		{
		    $_ = <RESULT>;
		    $updated_structure = 1 ;
		    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		    {
			$_ = <RESULT>;
			@line = split " ", $_ ;
			$species[$i_step][$i_atom] = $line[4] ;
			for ($i_coords = 1; $i_coords <=3; $i_coords++)
			{
			    $coords[$i_atom][$i_coords] = @line[$i_coords] ;
			    $coords_relax[$i_step][$i_atom][$i_coords] = @line[$i_coords] ;
			}
		    }
		    $i_step++ ;
		}
		if (/Have\ a\ nice\ day/)
		{
		    $found = 1 ;
		}
		if (/WARNING\!\ SELF\-CONSISTENCY\ CYCLE\ DID\ NOT\ CONVERGE/)
		{
		    $conv = 0 ;
		}
		if (/Number\ of\ self\-consistency\ cycles/)
		{
		    @line = split " ", $_ ;
		    $self_consistency_cycles += $line[6] ;
		}
		if (/Number\ of\ relaxation\ steps/)
		{
		    @line = split " ", $_ ;
		    $relaxation_steps += $line[6] ;
		}
		if (/Maximum\ number\ of\ relaxation/)
		{
		    $status = 2 ;
		}
	    }
	    close (RESULT) ;
	    
	    # if not converged but some relaxation steps have been successful, restart with last geometry...
	    if (($found eq 0) || ($conv eq 0))
	    {
		if ($updated_structure eq 1) 
		{
		    print STDOUT "rewrite geometry_file...\n" ;
		    &rewrite_geometry_file() ;
		} else {
		    # cannot converge at all, try a new geometry...
		    print STDOUT "$startindex Not converged?\n" ;
#		$finished = 1 ;
		    $enough = 1 ;
		    $status = 2 ;
		}
            # if optimization was aborted because number of relaxation steps has been exceeded, try once more...
	    } elsif ($status eq 2) {
		if ($abort eq 0) {
		    print STDOUT "rewrite geometry_file...\n" ;
		    &rewrite_geometry_file() ;
		    $abort = 1 ;
		} else {
		    # give up...
		    $enough = 1 ;
		    $status = 2 ;
		}
	    } else {
		# yippiee...
		$enough = 1 ;
	    }
	    
	} # until $enough...
	$n_steps = $i_step ;
	&write_relax_movie_file() ;
	
	open (INDEX, $STARTOUT) ;
	print INDEX ++$startindex, " ", $seed_uni, " ", $seed_poiss, " " , $relaxation_steps, " ", $self_consistency_cycles, "\n" ;
	close (INDEX) ;
	
#	system "cp $STARTIN $homedir" ;
	system "cp restart.dat $restart_file_name" ;
#	system "cp $restart_file_name $homedir" ;
	
	open (DFT_OUT, ">>$dft_file") ;
	
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $coords[$i_atom][1], $coords[$i_atom][2], $coords[$i_atom][3] ;
	}
	printf DFT_OUT ("%30.15e %4d\n") , $energy[$n_steps-1], $status;
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $forces[$i_atom][1], $forces[$i_atom][2], $forces[$i_atom][3] ;
	}
	if ($spin_polarized == 1) 
	{
	    for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	    {
		printf DFT_OUT ("%30.15e\n") , $spin_per_atom[$i_atom] ; 
	    }
	}
	close(DFT_OUT) ;
	
#	system "cp $dft_file $homedir" ;
	
	close (RESULT) ;
	
	$average_relaxation_steps = ($relaxation_steps + 1) / $startindex ;
	$average_self_consistency_cycles = $self_consistency_cycles / ($startindex * $average_relaxation_steps) ;
	
	print STDOUT "---\n" ;
	printf STDOUT ("average # of relaxation steps per BH-loop             : %10.4f\n") , $average_relaxation_steps ;
	    printf STDOUT ("average # of self-consistency-cycles per configuration: %10.4f\n") , $average_self_consistency_cycles ;
	}
}

system "rm -f finished.tmp" ;
system "rm -f geometry.old" ;
system "rm -f opt_data.out" ;
system "rm -f $out_file" ;

#####################

sub rewrite_geometry_file() {
    system "cp geometry.in geometry.in.temp" ;

    open (OLD_GEOMETRY, "geometry.in.temp") ;
    open (NEW_GEOMETRY, ">geometry.in") ;

    $i_atom = 0 ;

    while (<OLD_GEOMETRY>) {
	if (/atom/) {
	    printf NEW_GEOMETRY ("atom %30.15e %30.15e %30.15e %s\n") , $coords[$i_atom][1], $coords[$i_atom][2], $coords[$i_atom][3], $species[1][$i_atom] ; 
	    $i_atom ++ ;
	} else {
	    print NEW_GEOMETRY $_ ;
	}
    }

    close (OLD_GEOMETRY) ;
    close (NEW_GEOMETRY) ;

#    system "rm geometry.in.temp" ;

}

################################

sub write_relax_movie_file() {
    
   open (RELAX, $movie_file_name) ;

   print RELAX "[Molden Format]\n" ;
   print RELAX "[GEOCONV]\n" ;

   print RELAX "energy\n" ;
   for ($i_step = 0; $i_step < $n_steps; $i_step++)
   {
       printf RELAX ("%30.15e\n"), $energy[$i_step] ;
   }
   
   print RELAX "max-force\n" ;
   for ($i_step = 0; $i_step < $n_steps; $i_step++)
   {
       printf RELAX ("%30.15e\n"), $max_force[$i_step] ;
   }
   
   print RELAX "[GEOMETRIES](XYZ)\n" ;
   
   for ($i_step = 0; $i_step < $n_steps; $i_step++)
   {
       printf RELAX ("  %4d\n"), $n_atoms ;
       printf RELAX "\n" ;
       for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
       {
	   $_ = <RESULT>;
	   @line = split " ", $_ ;
	   printf RELAX ("%s %30.15e %30.15e %30.15e\n") , $species[$i_step][$i_atom], $coords_relax[$i_step][$i_atom][1], 
	   $coords_relax[$i_step][$i_atom][2], $coords_relax[$i_step][$i_atom][3] ;
       }
   }
   
   close (RELAX) ;

}
