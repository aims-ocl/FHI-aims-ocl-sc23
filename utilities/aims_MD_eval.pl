#!/usr/bin/perl

$infile = @ARGV[0] ;

open (INPUT,"<$infile") ;

#  always do printout of headline first ...
	$time = "#      Time [ps]";
	$T    = "Temperature [K]";
	$FE   = "E_tot (electronic) [eV]";
	$KE   = "E_kin (nuclei) [eV]";
	$TE   = "Total Energy [eV]";
	$TE_change   = "Total Energy Change [eV]";
	printf "%16s %16s %25s %25s %25s %25s\n", $time, $T, $FE, $KE, $TE, $TE_change;

$TE_zero = 0. ;

$iter = 0 ;
while (<INPUT>)
{

    if (!/\#/)
    {
      if (/Initial\ conditions\ for\ Born-Oppenheimer\ Molecular\ Dynamics:/)
      {
	&print_MD_status;
      }
      if (/Advancing\ structure\ using\ Born-Oppenheimer\ Molecular\ Dynamics/)
      {
  	<INPUT>;
	&print_MD_status;
      }
    }
	
}
close (INPUT) ;

sub print_MD_status
{
    <INPUT>;
    @line = split ' ', <INPUT> ;
    $time = $line[4] ;
    @line = split ' ', <INPUT> ;
    $FE   = $line[5] ;
    @line = split ' ', <INPUT> ;
    $T    = $line[4] ;
    @line = split ' ', <INPUT> ;
    $KE   = $line[5] ;
    @line = split ' ', <INPUT> ;
    $TE   = $line[5] ;
    if ($iter == 0 )
    {
	$TE_zero = $TE ;
        $iter = 1 ;
    }
    @line = split ' ', <INPUT>;
    if ($line[2]=~/value/ || $line[1] =~ /BDP/)
    {
	if ($line[2]=~/value/) {  # Nose-Hoover; skip two lines
	    $s    = $line[5];
	    @line = split ' ', <INPUT>;
	    $ps   = $line[5];
	    @line = split ' ', <INPUT>;
	}
        if ($line[1]=~/Nose-Hoover/ || $line[1] =~ /BDP/)
        {
            $H_NH = $line[4];
            printf "%18.8f %16.6f %25.6f %25.6f %25.6f %25.6f %25.6f\n", $time, $T, $FE, $KE, $TE, $H_NH, $TE-$TE_zero;
        }
        else
        {
            printf "%18.8f %16.6f %25.6f %25.6f %25.6f %25.6f\n", $time, $T, $FE, $KE, $TE, $TE-$TE_zero;
        }
    }
    else
    {
        printf "%16.6f %16.6f %25.6f %25.6f %25.6f %25.6f\n", $time, $T, $FE, $KE, $TE, $TE-$TE_zero;
    }
}
