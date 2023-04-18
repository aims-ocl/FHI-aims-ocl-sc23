#!/usr/bin/perl

# USAGE:  ./get_relaxation_info.pl   aims_output_file.out
#
# Revision 2010/12: VA 
#
# Fixme: rejection for TRM method


use strict;
use warnings;

our $infile = $ARGV[0] ;

if ( ! -e $infile)
{
    print "\n*** Error: Input file '", $infile, "' does not appear to exist - please check.\n\n" ;
    exit; 
}

open INPUT, '<', "$infile" ;

our $optimizer=0 ;
our $geo_rejected  = 0;
our $disaster = 0 ;
our $stuck = 0 ;
our $first = 1 ;
our $n_rel = 0 ;
our $energy = 0 ;
our $first_energy = 0 ;
our $free_energy = 0 ;
our $first_free_energy = 0 ;
our $max_force   = 0 ;

print "\n\# Step Total energy [eV]   E-E(1) [meV]   Free energy [eV]   F-F(1) [meV]   max. force [eV/AA]\n \n" ;

while (<INPUT>)
{
    & check_optimizer ;
    & read_energies ;
    & read_forces ;

    if (/Advancing/)
    {
        if (!/\#/)
        {

            & check_rejection ;

            printf "%5i   %16.8f %14.6f   %16.8f %14.6f %10.6f ", 
                $n_rel, $energy, 
                1000*($energy - $first_energy), 
                $free_energy,
                1000*($free_energy - $first_free_energy),
                $max_force ;

            if ($geo_rejected == 1)
               { print " rejected. " ;}
            if ($disaster  == 1)
               { print " rejected: force <-> energy inconsistency? " ;}
            if ($stuck  == 1)
               { print " stuck. " ;}
            print "\n" ;
        }
    }
    elsif (/Present\ geometry\ is\ converged\./)
    {


        printf "%5i   %16.8f %14.6f   %16.8f %14.6f %10.6f ", 
                $n_rel, $energy, 
                1000*($energy - $first_energy), 
                $free_energy,
                1000*($free_energy - $first_free_energy),
                $max_force ;

        print " converged. \n" ;

    }
    elsif (/\*\ Aborting\ optimization./)
    {


        printf "%5i   %16.8f %14.6f   %16.8f %14.6f %10.6f ", 
                $n_rel, $energy, 
                1000*($energy - $first_energy), 
                $free_energy,
                1000*($free_energy - $first_free_energy),
                $max_force ;

        print " Aborted - too many steps. \n" ;

    }
}

#----------------------------------------------------------------------------------
# Subroutines
#----------------------------------------------------------------------------------

sub check_optimizer
{  
    if (/Geometry\ relaxation:\ Textbook BFGS/)
       {$optimizer = 1 } 
    elsif (/Geometry\ relaxation:\ Modified BFGS\ -\ TRM/)
       {$optimizer = 2} 
}

sub read_energies 
{
    if (/\|\ Total\ energy\ uncorrected/)
    {
        my @line = split " ", $_ ;
        $energy = $line[5] ;

        $_ = <INPUT>;
        $_ = <INPUT>;

        @line = split " ", $_ ;
        $free_energy = $line[5] ;

        if ($first==1)
        {
            $first_energy = $energy ;
            $first_free_energy = $free_energy ;
	    $first=0 ;
        }
    }
} 

sub read_forces 
{
   if (/Maximum\ force\ component\ is/)
   {
        my @line = split " ", $_ ;
        $max_force = $line[4] ;
        $n_rel++ ;
   }
}

sub check_rejection 
{
   $geo_rejected = 0 ;
   $disaster     = 0 ;
   $stuck        = 0 ;
   # Optimizer: BFGS
   if ($optimizer == 1)
   {
      $_ = <INPUT>;
      $_ = <INPUT>;

      my @line = split " ", $_ ;
      if ( $line[1] eq 'Perform'  )
         { $_ = <INPUT> ;}

      $_ = <INPUT>;

      @line = split " ", $_ ;
      if ( ! (($line[0] eq 'Predicted') || ($line[0] eq 'x')) )
         {$geo_rejected = 1 ;}
   }
   # Optimizer: TRM
   elsif ($optimizer == 2)
   {for (my $i = 0; $i<10 ; $i++ )
         {
         $_ = <INPUT>;
	 if (!$_) {last;}
         if ( (/\|\ Contraproductive\ step/) || (/\|\ Counterproductive\ step/) )
            {$geo_rejected = 1 ;}
         elsif( /\*\*/ )
            {$disaster = 1;}
            # This should only happen for BFGS->TRM switch
	 elsif(/Optimizer\ is\ stuck/)
	    {$stuck = 1;}
         }
   }
}
