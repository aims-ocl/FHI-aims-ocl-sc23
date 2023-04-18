#!/usr/bin/perl

# turns cube file header units from Angstroms to bohr

$in_file = "$ARGV[0]" ;

#

open(INPUT,"$in_file") ;

$_ = <INPUT> ;
print $_ ;

$_ = <INPUT> ;
print $_ ;

$_ = <INPUT> ;
@line = split " ", $_ ;
$n_atoms = @line[0];
$origin[1] = @line[1]/0.52917721 ;
$origin[2] = @line[2]/0.52917721 ;
$origin[3] = @line[3]/0.52917721 ;

printf "%5i %10.6f %10.6f %10.6f\n", $n_atoms, $origin[1], $origin[2], $origin[3];

for ( $i_vec = 1; $i_vec<=3; $i_vec++ ) 
{
    $_ = <INPUT> ;
    @line = split " ", $_ ;
    $n_div = @line[0] ;
    $voxel[1] = @line[1]/0.52917721 ;
    $voxel[2] = @line[2]/0.52917721 ;
    $voxel[3] = @line[3]/0.52917721 ;

    printf "%5i %10.6f %10.6f %10.6f\n", $n_div, $voxel[1], $voxel[2], $voxel[3] ;
}

for ($i_atom = 1; $i_atom<=$n_atoms; $i_atom++)
{
    $_ = <INPUT> ;
    @line = split " ", $_ ;
    $n_at = @line[0] ;
    $chk = @line[1] ;
    $pos[1] = @line[2]/0.52917721 ;
    $pos[2] = @line[3]/0.52917721 ;
    $pos[3] = @line[4]/0.52917721 ;

    printf "%5i %10.6f %10.6f %10.6f %10.6f\n", $n_at, $chk, $pos[1], $pos[2], $pos[3] ;

}

while (<INPUT>)
{
    print $_ ;
}

close(INPUT) ;
