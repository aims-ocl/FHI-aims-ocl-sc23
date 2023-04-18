#!/usr/bin/perl

$infile =  "$ARGV[0]";

open (IN,"<$infile") ;

&get_data;

close (IN) ;

&output_data ;

&plot_data ;

sub get_data
{

    $iteration = 0 ;
    $read_now = 0 ;
    $spin = 0 ;
    $start = 0 ;

    $dens_conv = 0 ;
    $eev_conv = 0 ;
    $etot_conv = 0 ;

    while(<IN>)
    {

        if (/Spin\ treatment\:\ Spin\ density\ functional\ theory\ \-\ collinear\ spins\./)
        {
	    $spin = 1 ;
        }

        if (/Convergence\ accuracy\ of\ self\-consistent\ charge\ density\:/)
        {
	    @line = split " ", $_ ;
 	    $dens_conv = @line[6] ;
        }

        if (/Convergence\ accuracy\ of\ sum\ of\ eigenvalues\:/)
        {
	    @line = split " ", $_ ;
 	    $eev_conv = @line[6] ;
        }

        if (/Convergence\ accuracy\ of\ total\ energy\:/)
        {
	    @line = split " ", $_ ;
 	    $etot_conv = @line[5] ;
        }

        if (/Begin\ self\-consistency\ loop\:\ Initialization\./)
        {
            $start = 1 ;
        }

        if ($start == 1) 
        {
           if (/Self\-consistency\ convergence\ accuracy\:/)
           {
               $iteration++ ;
           }

           if ($spin == 1 )
           {
             if (/\|\ Change\ of\ charge\/spin\ density\ \:\ /)
             {
	       @line = split " ", $_ ;
               $density_change[$iteration] = @line[6] ;
               $spin_density_change[$iteration] = @line[7] ;
             }
           }
           else
           {
             if (/\|\ Change\ of\ charge\ density\ \ \ \ \ \ \:\ /)
             {
	       @line = split " ", $_ ;
               $density_change[$iteration] = @line[6] ;
             }
           }

           if (/\|\ Change\ of\ sum\ of\ eigenvalues\ \ \:\ /)
           {
	       @line = split " ", $_ ;
               $eev_change[$iteration] = @line[7] ;
           }

           if (/\|\ Change\ of\ total\ energy\ \ \ \ \ \ \ \ \:\ /)
           {
	       @line = split " ", $_ ;
               $etot_change[$iteration] = @line[6] ;
           }

        }

    }

}

sub output_data
{

    open (OUT,">convergence_data.dat") ;
 
    print OUT "# Iteration    Etot_change[eV]   Eigenvalues_change[eV]     Density_change  \n" ;

    for ( $i_iter = 1; $i_iter<=$iteration; $i_iter++ )
    {
        printf OUT "%8i  %20.10f  %20.10f  %20.10f  ", $i_iter, abs($etot_change[$i_iter]), abs($eev_change[$i_iter]),        abs($density_change[$i_iter]) ;

        if ($spin==1)
        {
	    printf OUT "%20.10f \n", abs($spin_density_change[$i_iter]);
        }
        else
        {
            print OUT "\n" ;
        }

    }

    close (OUT) ;

}

sub plot_data
{

    $it_min = 0 ;
    $it_max = $iteration * 1.1 ;

    open (TEMPLATE,"<scf_convergence_template.agr") ;
    open (OUT,">convergence_data.agr") ;

    while (<TEMPLATE>)
    {

        if (/\@\ \ \ \ \world/)
        {
	       @line = split " ", $_ ;
               print OUT "@    world $it_min, @line[3] $it_max, @line[5]\n" ;
        }
        elsif (/\?\?\?it\_min\?\?\?\ \?\?\?dens\_conv\?\?\?/)
        {
               print OUT "    $it_min    $dens_conv\n" ;
        }
        elsif (/\?\?\?it\_max\?\?\?\ \?\?\?dens\_conv\?\?\?/)
        {
               print OUT "    $it_max    $dens_conv\n" ;
        }
        elsif (/\?\?\?it\_min\?\?\?\ \?\?\?eev\_conv\?\?\?/)
        {
               print OUT "    $it_min    $eev_conv\n" ;
        }
        elsif (/\?\?\?it\_max\?\?\?\ \?\?\?eev\_conv\?\?\?/)
        {
               print OUT "    $it_max    $eev_conv\n" ;
        }
        elsif (/\?\?\?it\_min\?\?\?\ \?\?\?etot\_conv\?\?\?/)
        {
               print OUT "    $it_min    $etot_conv\n" ;
        }
        elsif (/\?\?\?it\_max\?\?\?\ \?\?\?etot\_conv\?\?\?/)
        {
               print OUT "    $it_max    $etot_conv\n" ;
        }
        elsif (/\?\?\?density\_data\?\?\?/)
        {

             for ( $i_iter = 1; $i_iter<=$iteration; $i_iter++ )
             {
                 printf OUT "%8i  %20.10f  \n", $i_iter, abs($density_change[$i_iter]) ;
            
             }

        }
        elsif (/\?\?\?spin\_density\_data\?\?\?/)
        {

             if ( $spin == 1)
             {
                 for ( $i_iter = 1; $i_iter<=$iteration; $i_iter++ )
                 {
                     printf OUT "%8i  %20.10f  \n", $i_iter, abs($spin_density_change[$i_iter]) ;
                 }
             }
             else
             {
                     printf OUT "%8i  %20.10f  \n", 0, 0.000000001  ;
             }

        }
        elsif (/\?\?\?eev\_data\?\?\?/)
        {
             for ( $i_iter = 1; $i_iter<=$iteration; $i_iter++ )
             {
                 printf OUT "%8i  %20.10f  \n", $i_iter, abs($eev_change[$i_iter]) ;            
             }

        }
        elsif (/\?\?\?etot\_data\?\?\?/)
        {
             for ( $i_iter = 1; $i_iter<=$iteration; $i_iter++ )
             {
                 printf OUT "%8i  %20.10f  \n", $i_iter, abs($etot_change[$i_iter]) ;            
             }

        }
        else
        {
	    print OUT $_ ;
        }

    }

    close (OUT) ;
    close (TEMPLATE) ;

    system "xmgrace convergence_data.agr" ;

}
