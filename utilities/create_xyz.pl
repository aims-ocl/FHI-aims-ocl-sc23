#!/usr/bin/perl

# General settings

$n_per[1] = 2 ;
$n_per[2] = 2 ;
$n_per[3] = 2 ;

#

$infile = @ARGV[0] ;

open (INPUT,"<$infile") ;

$n_atoms = 0 ;
$dimensions = 0 ;

@atom = () ;
@lattice = () ;

while (<INPUT>)
{
  if (/atom/)
  {
      @line = split " ",$_ ;
      if ($#line==4)
      {
	  $n_atoms = $n_atoms+1 ;
          $atom[$n_atoms][1] = @line[1] ;
          $atom[$n_atoms][2] = @line[2] ;
          $atom[$n_atoms][3] = @line[3] ;
          $atom[$n_atoms][0] = @line[4] ;
#        print @line, "\n" ;
      }
  }
  elsif(/lattice\_vector/)
  {
      @line = split " ",$_ ;
      if ($#line==3)
      {
	  $dimensions = $dimensions+1 ;
          $lattice[$dimensions][1] = @line[1] ;
          $lattice[$dimensions][2] = @line[2] ;
          $lattice[$dimensions][3] = @line[3] ;
#        print @line, "\n" ;
      }
  }
}
close(INPUT) ;

$n_sites = $n_atoms ;
for ($i_dim=1; $i_dim<=$dimensions; $i_dim++)
{
$n_sites = $n_sites * $n_per[$i_dim]
}

print $n_sites, "\n" ;
print "# xyz file created automatically from geometry.in.\n" ;

# The following is a really clumsy hand-written recursive loop that accounts for 0, 1, 2, or 3 lattice_vectors inside geometry.in.

if ($dimensions==0) 
{

   for ($i_atom=1; $i_atom<=$n_atoms; $i_atom++)
   {
       printf "%s  %16.8f %16.8f %16.8f \n", $atom[$i_atom][0], $atom[$i_atom][1], $atom[$i_atom][2], $atom[$i_atom][3] ;
   }
}
else
{

   for ( $idim1=0; $idim1<$n_per[1]; $idim1++)
   {
       $offset_1[1] = $lattice[1][1] * $idim1 ;
       $offset_1[2] = $lattice[1][2] * $idim1 ;
       $offset_1[3] = $lattice[1][3] * $idim1 ;

       if ($dimensions==1)
       {
           for ($i_atom=1; $i_atom<=$n_atoms; $i_atom++)
           {
               printf "%s  %16.8f %16.8f %16.8f \n", $atom[$i_atom][0], 
                 $atom[$i_atom][1]+$offset_1[1], 
                 $atom[$i_atom][2]+$offset_1[2], 
                 $atom[$i_atom][3]+$offset_1[3] ;
            }

       }
       else
       {
           for ( $idim2=0; $idim2<$n_per[2]; $idim2++)
           {
              $offset_2[1] = $lattice[2][1] * $idim2 ;
              $offset_2[2] = $lattice[2][2] * $idim2 ;
              $offset_2[3] = $lattice[2][3] * $idim2 ;

              if ($dimensions==2)
              {
                   for ($i_atom=1; $i_atom<=$n_atoms; $i_atom++)
                   {
                       printf "%s  %16.8f %16.8f %16.8f \n", $atom[$i_atom][0], 
                          $atom[$i_atom][1]+$offset_1[1]+$offset_2[1], 
                          $atom[$i_atom][2]+$offset_1[2]+$offset_2[2], 
                          $atom[$i_atom][3]+$offset_1[3]+$offset_2[3] ;
                    }

              }
              else
              {

                  for ( $idim3=0; $idim3<$n_per[3]; $idim3++)
                  {
                      $offset_3[1] = $lattice[3][1] * $idim3 ;
                      $offset_3[2] = $lattice[3][2] * $idim3 ;
                      $offset_3[3] = $lattice[3][3] * $idim3 ;

                      for ($i_atom=1; $i_atom<=$n_atoms; $i_atom++)
                      {
                          printf "%s  %16.8f %16.8f %16.8f \n", $atom[$i_atom][0], 
                          $atom[$i_atom][1]+$offset_1[1]+$offset_2[1]+$offset_3[1], 
                          $atom[$i_atom][2]+$offset_1[2]+$offset_2[2]+$offset_3[2], 
                          $atom[$i_atom][3]+$offset_1[3]+$offset_2[3]+$offset_3[3] ;
                      } 
                 
		  }
              }

          }

       }

   }

}

