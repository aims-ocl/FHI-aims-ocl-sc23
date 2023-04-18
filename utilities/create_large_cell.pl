#!/usr/bin/perl
#
# Takes small unit cell (fractional coordinates only) and creates larger cell
#

$factor[1] = 2 ;
$factor[2] = 2 ;
$factor[3] = 2 ;

# read input

$n_lattice_vectors = 0 ;

$n_atoms = 0 ;

if (! -e "geometry.in.basic" )
{
	 print "Error: File geometry.in.basic not found.\n" ;
	 print "Need file geometry.in.basic.\n" ;
         die;
}

open INPUT, '<', "geometry.in.basic" ;

while (<INPUT>)
{

  if (/lattice\_vector/)
  {
      $n_lattice_vector++ ;

      @line = split ' ', $_ ;

      $lattice_vector[1][$n_lattice_vector] = $line[1] ;
      $lattice_vector[2][$n_lattice_vector] = $line[2] ;
      $lattice_vector[3][$n_lattice_vector] = $line[3] ;
  }

  if (/atom/)
  {

     if (/atom\_frac/)
     {

	 $n_atoms++ ;

         @line = split ' ', $_ ;

         $coord_frac[1][$n_atoms] = $line[1] ;
         $coord_frac[2][$n_atoms] = $line[2] ;
         $coord_frac[3][$n_atoms] = $line[3] ;
         $species[$n_atoms] = $line[4] ;

     }
     else
     {
	 print "Error: atom encountered, but not atom\_frac.\n" ;
	 print "Only fractional coordinates are currently supported.\n" ;
         die;
     }
  }

}

close(INPUT) ;

if ( $n_lattice_vector ne 3 )
{
	 print "Error: Not 3 lattice vectors found.\n" ;
	 print "Need exactly 3 lattice vectors.\n" ;
         die;    
}

if ( $n_atoms eq 0 )
{
	 print "Error: No atoms found.\n" ;
	 print "Need atoms.\n" ;
         die;    
}

# write output

for ( $i_lv = 1; $i_lv <=3 ; $i_lv++ )
{
    printf "lattice\_vector  %20.8f %20.8f %20.8f\n", $lattice_vector[1][$i_lv]*$factor[$i_lv], $lattice_vector[2][$i_lv]*$factor[$i_lv], $lattice_vector[3][$i_lv]*$factor[$i_lv];

}

print "\n" ;

    for ( $i_cell_1 = 0 ; $i_cell_1 < $factor[1] ; $i_cell_1++ )
    {
        for ( $i_cell_2 = 0 ; $i_cell_2 < $factor[2] ; $i_cell_2++ )
        {
            for ( $i_cell_3 = 0 ; $i_cell_3 < $factor[3] ; $i_cell_3++ )
            {

                for ( $i_atom = 1 ; $i_atom <= $n_atoms ; $i_atom++ )
                { 

                    $i_cell[1] = $i_cell_1 ;
                    $i_cell[2] = $i_cell_2 ;
                    $i_cell[3] = $i_cell_3 ;

                    for ( $i_lv = 1; $i_lv <=3 ; $i_lv++)
                    {
                        $frac[$i_lv] = $coord_frac[$i_lv][$i_atom]/$factor[$i_lv] + $i_cell[$i_lv]/$factor[$i_lv] ;
                    }

                    printf "atom\_frac  %20.8f %20.8f %20.8f %s \n", $frac[1], $frac[2], $frac[3], $species[$i_atom] ;

                }
   
            }
        }
    }

