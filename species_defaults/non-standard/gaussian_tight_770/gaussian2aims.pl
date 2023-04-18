#!/usr/bin/perl 
# use strict;
use warnings;

# This script converts gaussian format in FHI-AIMS
# To use the script, one put the gaussian basis in Gaussian94 format
# into a file called "gaussian.in", and then run the script
# ./gaussian2aims.pl <element> <basis>, e.g.,
# ./gaussian2aims.pl N cc-pVTZ
# The output file "N.cc-pVTZ" contains the same gausisan basis
# in FHI-aims format.

if($#ARGV<1) 
 {
  open (OUTPUT, ">gaussian.out") ;
 }
else
 {
  $filename= join "", $ARGV[0], ".", $ARGV[1];
  open (OUTPUT, ">$filename") ;
  print OUTPUT "# ", $ARGV[0], " ", $ARGV[1], "\n" ;
 }

@datagau = ();

open (INPUT, "gaussian.in") ;
$i_pool = 1 ;
$i_pool_1 = 1 ;
while (<INPUT>)
 {
     @line = split " ", $_ ;
     if($#line>0) 
     {
       for ($i_line=0; $i_line<=$#line ; $i_line++)
       {
         $datagau[$i_pool][$i_line] = $line[$i_line];
       }
       $i_pool ++;
     }
 }
 for ($i_pool_1=1; $i_pool_1<=$i_pool-1; $i_pool_1++)
 {
#    if ( $datagau[$i_pool_1][0] eq 'S') &&  ($datagau[$i_pool_1][1] >1) )
#    print $datagau[$i_pool_1][0], "\n";
#    print $i_pool_1, "\n";
    if ($datagau[$i_pool_1][0] eq "S" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '0 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'S' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '0 ', '1 ',  $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "P" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '1 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'P' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '1 ',  '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "D" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '2 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'D' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '2 ', '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "F" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '3 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'F' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '3 ', '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "G" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '4 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'G' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '4 ', '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "H" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '5 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'H' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '5 ', '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
    if ($datagau[$i_pool_1][0] eq "I" && $datagau[$i_pool_1][1] > 1 )
    {
     print OUTPUT ' gaussian ',  '6 ', $datagau[$i_pool_1][1], "\n";  
     for ($i_contra=1; $i_contra<=$datagau[$i_pool_1][1]; $i_contra++)
     {
      printf OUTPUT "%18.7f %20.7f  \n", $datagau[$i_pool_1+$i_contra][0], $datagau[$i_pool_1+$i_contra][1];
     }
    }
    elsif ($datagau[$i_pool_1][0] eq 'I' && $datagau[$i_pool_1][1] == 1 )
    {
     print OUTPUT ' gaussian ',  '6 ', '1 ', $datagau[$i_pool_1+1][0], "\n";  
    }
 }

close (INPUT);
close (OUTPUT);
