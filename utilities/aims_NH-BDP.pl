#!/usr/bin/perl

$infile = @ARGV[0] ;

open (INPUT,"<$infile") ;

$time = "#      Time [ps]";
#$T    = "Temperature [K]";
$T    = "Temperature (nuclei)";
$FE   = "elec. Free Energy [eV]";
$KE   = "Nuclear kinetic energy";
$TE   = "Total Energy [eV]";
$NHH  = "Conserved Hamilt. [eV]";
printf "%15s %15s %23s %23s %23s %23s\n", $time, $T, $FE, $KE, $TE, $NHH;

$iter = 0 ;
while (<INPUT>)
{

    if (!/\#/)
    {
      if (/Advancing\ structure\ using\ Born-Oppenheimer\ Molecular\ Dynamics/)
      {
	<INPUT>;
	<INPUT>;
	&print_MD_status;
      }
    }
	
}
close (INPUT) ;

sub print_MD_status
{
    @line = split ' ', <INPUT> ;
    $time = @line[4] ;
    @line = split ' ', <INPUT> ;
    $FE   = @line[5] ;
    @line = split ' ', <INPUT> ;
    $T = @line[4] ;
    @line = split ' ', <INPUT> ;
    $KE   = @line[5] ;
    @line = split ' ', <INPUT> ;
    $TE   = @line[5] ;
    $NHH = 0.0;
    $_ = <INPUT> ;
    if(/BDP\ pseudo\-Hamiltonian/) {
        @line = split ' ', $_;
        $NHH =  @line[4] ;
    }
    else {
      $_ = <INPUT> ;
      $_ = <INPUT> ;
      if(/Nose-Hoover\ Hamiltonian/) {
    	@line = split ' ', $_;    
	$NHH =  @line[4] ; 
      }
    }
    printf "%15.6f %15.6f %23.6f %23.6f %23.6f %23.6f\n", $time, $T, $FE, $KE, $TE, $NHH;
}
