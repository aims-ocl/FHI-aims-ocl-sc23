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
# - geometry.old (last structure with dft-energy and -forces available)
# - geometry.opt (geometry.old plus opt_data.out from get_cg_geometry.x)
#                (so energy and max_force)
# - temp.out (temporary logfile for localorb)
#
# temporary files are deleted afterwards

$homedir = @ARGV[0] ;
$EXE_AIMS = @ARGV[1] ;
$EXE_BH = @ARGV[2] ;

@line = split "/", $EXE_AIMS ;
$exe_aims_local = @line[$#line];

@line = split "/", $EXE_BH ;
$exe_bh_local = @line[$#line];

# file declarations
$out_file = "temp.out" ;
$bh_log_file = "bh_log.out" ;
$log_file = "log.out" ;
$dft_file = "dft_data.dat" ;

$STARTIN = "startindex.dat";
$STARTOUT = ">startindex.dat";

# seed for RNG
$seed_uni   = 197 ;
$seed_poiss = 276 ;

$startindex = 0 ;
$finished = 0 ;
$status = 0 ;

$n_atoms = 0 ;

$relaxation_steps = 0 ;
$self_consistency_cycles = 0 ;

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
system "cp geometry.in.basic geometry.old" ;

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
    system "cp $bh_log_file $homedir" ;
    system "cat opt_data.out geometry.old > geometry.opt" ;

    if (-e "finished.tmp") {
	$finished = 1 ;
    }
    if ($finished eq 0) {
	system "poe ./$exe_aims_local > $out_file 2>&1" ;
	system "cat $out_file >> $log_file " ;
    
	$movie_file_name = join "", ">relax.", $startindex,".irc" ;

	# check convergence
	@line = () ;
	$found = 0 ;
	$conv  = 1 ;
	$status = 0 ;
	open (RESULT, "$out_file") ;
	open (MOVIE_FILE, "$movie_file_name") ;
	while (<RESULT>)
	{
	    if (/Total\ energy\ corrected/)
	    {
		@line = split " ", $_ ;
		$energy = @line[5] ;
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
		    for ($i_coords = 1; $i_coords <=3; $i_coords++)
		    {
			$coords[$i_atom][$i_coords] = @line[$i_coords] ;
		    }
		}
	    }
	    if (/Updated\ atomic\ structure/)
	    {
		$_ = <RESULT>;
		printf MOVIE_FILE ("  %4d\n"), $n_atoms ;
		printf MOVIE_FILE "\n" ;
		for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		{
		    $_ = <RESULT>;
		    @line = split " ", $_ ;
		    printf MOVIE_FILE ("%s %30.15e %30.15e %30.15e\n") , @line[4], @line[1], @line[2], @line[3] ;
		}
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
	    if (/Aborting\ optimization/)
	    {
		$status = 2 ;
	    }
	}
	close (RESULT) ;
	close (MOVIE_FILE) ;

	open (INDEX, $STARTOUT) ;
	print INDEX ++$startindex, " ", $seed_uni, " ", $seed_poiss, " " , $relaxation_steps, " ", $self_consistency_cycles, "\n" ;
	close (INDEX) ;

	system "cp $STARTIN $homedir" ;

	if (($found eq 0) || ($conv eq 0))
	{
	    print STDOUT "$startindex Not converged?\n" ;
	    $finished = 1 ;
	    $status = 2 ;
	}
	open (DFT_OUT, ">>$dft_file") ;
	
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $coords[$i_atom][1], $coords[$i_atom][2], $coords[$i_atom][3] ;
	}
	printf DFT_OUT ("%30.15e %4d\n") , $energy, $status ;
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $forces[$i_atom][1], $forces[$i_atom][2], $forces[$i_atom][3] ;
	}
	close(DFT_OUT) ;
	
	# store old geometry
	system "cp geometry.in geometry.old" ;
	system "cp $dft_file $homedir" ;

	close (RESULT) ;

	$average_relaxation_steps = $relaxation_steps / $startindex ;
	$average_self_consistency_cycles = $self_consistency_cycles / ($startindex * $average_relaxation_steps) ;
	
	print STDOUT "---\n" ;
	printf STDOUT ("average # of relaxation steps per BH-loop : %10.4f\n") , $average_relaxation_steps ;
	    printf STDOUT ("average # of self-consistency-cycles per configuration: %10.4f\n") , $average_self_consistency_cycles ;
	}
}


system "rm -f finished.tmp" ;
system "rm -f geometry.old" ;
system "rm -f opt_data.out" ;
system "rm -f $out_file" ;
