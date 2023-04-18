#!/usr/bin/perl

# VB: Really trivial script to grep for the number of self-consistency
#     cycles per iteration in a longer run.
#     Useful for debugging and optimization purposes

$iter = 0 ;

$now = 0 ;

      print "#   Step  s.c.f._iterations_taken \n" ;

while (<>) 
{

if (/Self\-consistency\ cycle\ converged/)
{
    $now = 1 ;
}

if ($now == 1)
{
  if (/iteration/)
  {
      @line = split " ", $_ ;
      $n_scf = @line[4] ;
      $iter++ ;
      $now = 0 ;

      printf ("%8i %8i \n",$iter, $n_scf) ;
  }  
}

}
