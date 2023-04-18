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

# localorb binary
$EXE_AIMS = '~/codes/fhi-aims/bin/aims.121106.mkl.x' ;

# get_cg_geometry binary
$EXE_BH = '~/bin/effernir.121106.x' ;

# file declarations
$out_file = "temp.out" ;
$cg_log_file = "bh_log.out" ;
$log_file = "log.out" ;
$dft_file = "dft_data.dat" ;

$STARTIN = "startindex.dat";
$STARTOUT = ">startindex.dat";

# seed for RNG
$seed = 11 ;

$startindex = 0 ;
$finished = 0 ;

$n_atoms = 0 ;

# read startindex file
if (-e $STARTIN)
{
    open (INDEX, $STARTIN) ;
    $startindex = <INDEX> ;
    close (INDEX) ;
} else {
    open (INDEX, $STARTOUT) ;
    print INDEX 0, $seed ;
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

    system "$EXE_BH > $cg_log_file 2>&1" ;

    system "cat opt_data.out geometry.old > geometry.opt" ;

    if (-e "finished.tmp") {
	$finished = 1 ;
    }
    if ($finished eq 0) {
	system "$EXE_AIMS > $out_file 2>&1" ;
	system "cat $out_file >> $log_file " ;
    
	open (INDEX, $STARTOUT) ;
	print INDEX ++$startindex ;
	close (INDEX) ;
	
	# check convergence
	@line = () ;
	$found = 0 ;
	$conv  = 1 ;
	open (RESULT, "$out_file") ;
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
			$forces[$i_atom][$i_coords] = @line[$i_coords] ;
		    }
		}
	    }
	    if (/\|\ Final\ atomic\ structure/) 
	    {
		$_ = <RESULT>;
		for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
		{
		    $_ = <RESULT>;
		    @line = split " ", $_ ;
		    for ($i_coords = 1; $i_coords <=3; $i_coords++)
		    {
			$coords[$i_atom][$i_coords] = @line[$i_coords + 3] ;
		    }
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
	}
	close (RESULT) ;
	if (($found eq 0) || ($conv eq 0))
	{
	    print STDOUT "$startindex Not converged?\n" ;
	    $finished = 1 ;
	}

	open (DFT_OUT, ">>$dft_file") ;
	
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $coords[$i_atom][1], $coords[$i_atom][2], $coords[$i_atom][3] ;
	}
	printf DFT_OUT ("%30.15e\n") , $energy ;
	for ($i_atom = 0; $i_atom < $n_atoms; $i_atom++)
	{
	    printf DFT_OUT ("%30.15e %30.15e %30.15e\n") , $forces[$i_atom][1], $forces[$i_atom][2], $forces[$i_atom][3] ;
	}
	close(DFT_OUT) ;
	
	# store old geometry
	system "cp geometry.in geometry.old" ;
	
	close (RESULT) ;
    }
}

system "rm -f finished.tmp" ;
system "rm -f geometry.old" ;
system "rm -f opt_data.out" ;
system "rm -f $out_file" ;
