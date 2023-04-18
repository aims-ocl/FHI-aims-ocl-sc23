#!/usr/bin/perl

#  Simple script that takes two different geometry files and an scale factor as input,
#  and linearly interpolates both geometries according to the scale factor. The resulting 
#  intermediate geometry is written to STDOUT. The difference between the starting and
#  ending geometries is written to a file geometry_difference.out . Calling sequence:
#
#  > intermediate_geometries.pl <file_1> <file_2> <scale>
#
#  where <scale>=0.0 produces the same geometry as <file_1>, and <scale>=1.0 produces
#  the same geometry as <file_2> . This script has few safeguards but can be enormously
#  useful. 

$geometry_file[1] = "$ARGV[0]" ;
$geometry_file[2] = "$ARGV[1]" ;

$scale = $ARGV[2] ;

# learn first geometry

$i_file = 1 ;

open (IN,"<$geometry_file[1]") ;

&learn_geo;

close (IN) ;

# learn second geometry

$i_file = 2 ;

open (IN,"<$geometry_file[2]") ;

&learn_geo;

close (IN) ;

&average_geometries ;

&geometry_difference ;


# end main file


sub learn_geo
{

  $i_periodic = 0 ;
  $i_atom = 0 ;

  while(<IN>)
  {

      if (/lattice\_vector/)
      {
          $i_periodic = $i_periodic+1 ;

          @line = split " ", $_ ;

          $lattice[$i_file][$i_periodic][1] = $line[1] ;
          $lattice[$i_file][$i_periodic][2] = $line[2] ;
          $lattice[$i_file][$i_periodic][3] = $line[3] ;
      }
      if (/atom\_frac/) 
      {
          $i_atom++ ;

          @line = split " ", $_ ;

          $coord[$i_file][$i_atom][1] = 
             $line[1] * $lattice[$i_file][1][1]
            +$line[2] * $lattice[$i_file][2][1]
            +$line[3] * $lattice[$i_file][3][1];

          $coord[$i_file][$i_atom][2] = 
             $line[1] * $lattice[$i_file][1][2]
            +$line[2] * $lattice[$i_file][2][2]
            +$line[3] * $lattice[$i_file][3][2];
          
          $coord[$i_file][$i_atom][3] = 
             $line[1] * $lattice[$i_file][1][3]
            +$line[2] * $lattice[$i_file][2][3]
            +$line[3] * $lattice[$i_file][3][3];
          
          $species[$i_file][$i_atom] = $line[4] ;
      } 
      elsif (/atom/)
      {
          $i_atom++ ;

          @line = split " ", $_ ;

          $coord[$i_file][$i_atom][1] = $line[1] ;
          $coord[$i_file][$i_atom][2] = $line[2] ;
          $coord[$i_file][$i_atom][3] = $line[3] ;

          $species[$i_file][$i_atom] = $line[4] ;

      }


  }

  $n_periodic[$i_file] = $i_periodic ;
  $n_atoms[$i_file] = $i_atom ;

  printf STDOUT "# File %3i has %3i lattice vectors and %8i atoms.\n", $i_file, $n_periodic[$i_file], $n_atoms[$i_file] ; 

} # end subroutine learn_geo

sub average_geometries
{

    # check
    if ( $n_periodic[1] != $n_periodic[2] )
    {
          print STDERR "Error: Inconsistent number of lattice vectors in input files:\n" ;
          print STDERR "File 1: $n_periodic[1] , File 2: $n_periodic[2] .\n" ;
          die;       
    } 

    if ( $n_atoms[1] != $n_atoms[2])
    {
          print STDERR "Error: Inconsistent number of atoms in input files:\n" ;
          print STDERR "File 1: $n_atoms[1], File 2: $n_atoms[2] .\n" ;
          die;       
    } 

    print STDOUT "# Average geometry of $geometry_file[1] and $geometry_file[2], scale: $scale \n" ;
    print STDOUT "# \n" ; 

    for ($i_periodic = 1; $i_periodic <= $n_periodic[1]; $i_periodic++)
    {

       for ($i_coord = 1; $i_coord <= 3; $i_coord++ )
       {
	   $lattice_avg[$i_coord] = $lattice[1][$i_periodic][$i_coord] + $scale * ( $lattice[2][$i_periodic][$i_coord] - $lattice[1][$i_periodic][$i_coord] ) ;
       }

       printf STDOUT "  lattice_vector %18.10f  %18.10f  %18.10f  \n", $lattice_avg[1],$lattice_avg[2],$lattice_avg[3]; 
    }

    print STDOUT "# \n" ; 

    for ($i_atom = 1; $i_atom <= $n_atoms[1]; $i_atom++)
    {

       for ($i_coord = 1; $i_coord <= 3; $i_coord++ )
       {
	   $coord_avg[$i_coord] = $coord[1][$i_atom][$i_coord] + $scale * ( $coord[2][$i_atom][$i_coord] - $coord[1][$i_atom][$i_coord] ) ;
       }

       if ( $species[1][$i_atom] ne $species[2][$i_atom] )
       {
           print STDERR "Atom $i_atom : Species $species[1][$i_atom] and $species[2][$i_atom] do not match!\n" ;
	   die; 
       }

       printf STDOUT "  atom  %18.10f  %18.10f  %18.10f  ", $coord_avg[1],$coord_avg[2],$coord_avg[3]; 
       print STDOUT $species[1][$i_atom], "\n" ;
    }

} # end subroutine average_geometries

sub geometry_difference
{

    open (DIFF_FILE, ">geometry_difference.out") ;

    print DIFF_FILE "# Difference of $geometry_file[1] and $geometry_file[2] :" ;
    print DIFF_FILE "# \n" ; 

    for ($i_periodic = 1; $i_periodic <= $n_periodic[1]; $i_periodic++)
    {

       for ($i_coord = 1; $i_coord <= 3; $i_coord++ )
       {
	   $lattice_avg[$i_coord] = ( $lattice[2][$i_periodic][$i_coord] - $lattice[1][$i_periodic][$i_coord] ) ;
       }

       printf DIFF_FILE "  lattice_vector %18.10f  %18.10f  %18.10f  \n", $lattice_avg[1],$lattice_avg[2],$lattice_avg[3]; 
    }

    print DIFF_FILE "# \n" ; 

    for ($i_atom = 1; $i_atom <= $n_atoms[1]; $i_atom++)
    {

       for ($i_coord = 1; $i_coord <= 3; $i_coord++ )
       {
	   $coord_avg[$i_coord] = ( $coord[2][$i_atom][$i_coord] - $coord[1][$i_atom][$i_coord] ) ;
       }

       printf DIFF_FILE "  atom  %18.10f  %18.10f  %18.10f  ", $coord_avg[1],$coord_avg[2],$coord_avg[3]; 
       print DIFF_FILE $species[1][$i_atom], "\n" ;
    }

    close(DIFF_FILE) ;

} # end subroutine learn_geo
