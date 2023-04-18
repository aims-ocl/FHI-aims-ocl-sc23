#!/usr/bin/perl

# This script outputs a .xyz movie when the input is a standard FHI-aims output
# The FHI-aims output can be the outcome of either 
# a geometry optimization or a molecular dynamics run

$infile = @ARGV[0] ;

open (INPUT,"<$infile") ;

$n_atoms = 0 ;

$iter = 0 ;
while (<INPUT>)
{

    if (/\|\ Number\ of\ atoms/)
    {
	@line = split ' ', $_ ;
        $n_atoms = @line[5] ;
    }

    # Only begin actual reading after control.in and geometry.in were
    # parsed and number of atoms is known.
    if ($n_atoms > 0) 
    {
      if (/Input\ geometry\:/)
      {
        $_ = <INPUT>;	
	if ($_=~/\|\ Unit\ cell\:/) 
	{
	    <INPUT>;
	    <INPUT>;
	    <INPUT>;
	}
        <INPUT>;	
        <INPUT>;	
	&print_geo ;
      }
      elsif (/\ Updated\ atomic\ structure\:/)
      {
        <INPUT>;
	&print_new_geo ;
      }
      elsif (/\ Final\ atomic\ structure\:/)
      {
        <INPUT>;
	&print_new_geo ;
      }
      elsif (/used in the preceding time step/)
      {
        <INPUT>;
        &print_new_geo ;
      }
    }

}

close (INPUT) ;

sub print_geo
{
	printf STDOUT "%12i\n",$n_atoms ;
        printf STDOUT "Iteration:%6i\n", $iter ;
        for ($i_atom=1;$i_atom<=$n_atoms;$i_atom++)
        {
	    $_ = <INPUT> ;
	    if (/velocity/) 
	    {
		$_ = <INPUT>;
		if (/atom/) 
		{
		    @line = split ' ', $_ ;
		    printf STDOUT " %s %16.6f %16.6f %16.6f\n",@line[3],@line[4],@line[5],@line[6] ;
		}
	    }
	    else
	    {
		@line = split ' ', $_ ;
		printf STDOUT " %s %16.6f %16.6f %16.6f\n",@line[3],@line[4],@line[5],@line[6] ;
	    }
        }
        $iter++ ;

}
sub print_new_geo
{
	printf STDOUT "%12i\n",$n_atoms ;
        printf STDOUT "Iteration:%6i\n", $iter ;
        for ($i_atom=1;$i_atom<=$n_atoms;$i_atom++)
        {
	    $_ = <INPUT> ;
	    if (/velocity/) 
	    {
		$_ = <INPUT>;
		if (/atom/) 
		{
		    @line = split ' ', $_ ;
		    printf STDOUT " %s %16.6f %16.6f %16.6f\n",@line[4],@line[1],@line[2],@line[3] ;
		}
	    }
       elsif (/Species/)
       {
         @line = split ' ', $_ ;
         printf STDOUT " %s %16.6f %16.6f %16.6f\n",@line[3],@line[4],@line[5],@line[6] ;
       }
	    else
	    {
		@line = split ' ', $_ ;
		printf STDOUT " %s %16.6f %16.6f %16.6f\n",@line[4],@line[1],@line[2],@line[3] ;
	    }
        }
        $iter++ ;

}
