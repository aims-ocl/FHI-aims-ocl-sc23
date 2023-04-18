#!/usr/bin/perl

# create_geometry_zip.pl
# Joerg Meyer (jm), meyer@fhi-berlin.mpg.de 
#
# revision 1.0  2009/08 jm
# initial version
#
# revision 1.1  2009/08 nr
# order of geometry files corrected
# animation time delayed
#
# revision 1.2  2010/12 va 
# adapt to lattice optimization
#
# revision by VB 2014
# - Read multiple initial structures (e.g., in concatenated files from
#   different runs
#   Cat'ing different files in the same trajectory behind one another
#   may lead to duplicate geometries in the trajectory. However, there
#   is nothing the script can do about this, as the script has no way of
#   knowing why a specific FHI-aims output file was cat'ed behind another one.
#
# Writes "geometry-%06i.in" files for every geometry in FHI-aims standard output
# including lattice vector(s) into a subdirectory "geometries" which gets 
# compressed into "geometries.zip" (and is removed afterwards).
# (The geometry.in files can be of course be used as inputs for a new aims calculation.)
# Additionally, a Jmol script "geometries.spt" is written which animates 
# the geometries in "geometries.zip" including the information about 
# periodic boundary conditions (if available) when called e.g. via
#
#	<jmol> -s geometries.spt
#
# (loosely) based on create_xyz_movie.pl (with hopefully improved perl quality :-)


use strict ;				# enforce good perl :-)
use warnings;
use File::Spec::Functions ;		# to make hardcoded paths platform independent (see below)

# process aims output and write geometry files
our $infile = $ARGV[0] ;
open INPUT, '<', "$infile" ;

our $outdir = "geometries" ;
mkdir $outdir ;

our $n_iter = 0 ;
our $n_lattice_vectors = 0 ;
our $n_atoms = 0 ;
our @lattice_vectors = () ;
our @atoms = () ;
while (<INPUT>)
{
    if (/\|\ Number\ of\ atoms/)
    {
	&read_number_of_atoms ;
    }
    elsif (/\|\ Number\ of\ lattice\ vectors/)
    {
	&read_number_of_lattice_vectors ;
    }
    elsif (/Input\ geometry\:/)
    {
        # This case (imput geometry) will only happen more than once
        # if the file we are looking at was cat'ed based on
        # several output files from several separate runs.

        $_ = <INPUT>;	
	if ($_=~/\|\ Unit\ cell\:/) 
	{
		&read_lattice_vectors ;
	}
        <INPUT>;	
        <INPUT>;	
	&read_input_atoms ;
	&write_geometry ;
    }
    elsif (/\ Updated\ atomic\ structure\:/)
    {
        <INPUT>;
	if ($n_lattice_vectors > 0)
        {
		&read_lattice_vectors ;
                <INPUT>;
        }
	&read_atoms ;
	&write_geometry ;
    }
#    Reminder (VB):
#    The "final atomic structure" does not need to be read as a separate
#    geometry because it should be identical to the last "updated atomic structure"
#
#    elsif (/\ Final\ atomic\ structure\:/)
#    {
#	if ($n_lattice_vectors > 0)
#        {
#		&read_lattice_vectors ;
#                <INPUT>;
#        }
#        <INPUT>;
#	&read_atoms ;
#	&write_geometry ;
#    }
}
close(INPUT) ;

# prepare previous output for visualisation with Jmol
our $zipfile = "geometries.zip" ;
my $geomfiles = catfile $outdir,"geometry-*";
system "zip -m $zipfile $geomfiles" ;
rmdir $outdir ;
&write_jmol_script ;


### subroutines
# (copy & paste into other scripts is limited due to heavy use of global variables!)

sub read_number_of_atoms
{
	my @line = split ' ', $_ ;
        $n_atoms = $line[5] ;
}

sub read_number_of_lattice_vectors 
{
	my @line = split ' ', $_ ;
        $n_lattice_vectors = $line[6] ;
}

sub read_lattice_vectors
{
#	printf STDOUT "reading %3i lattice vectors \n", $n_lattice_vectors ;
	@lattice_vectors = () ;
	for ( my $i_lv=0; $i_lv<3; $i_lv++ )
	{
	    $_ = <INPUT>; 
	    my @line = split ' ', $_ ;
	    push @lattice_vectors, [($line[1],$line[2],$line[3])] ;
	}
}

sub write_lattice_vectors
{
	for ( my $i_lv=0; $i_lv<$n_lattice_vectors; $i_lv++ )
	{
	    my $x = $lattice_vectors[$i_lv][0] ;
	    my $y = $lattice_vectors[$i_lv][1] ;
	    my $z = $lattice_vectors[$i_lv][2] ;
	    printf OUTPUT "lattice_vector %16.6f %16.6f %16.6f\n", $x, $y, $z ;
	}
}

sub read_atoms
{
#	printf STDOUT "reading %5i atoms from iteration %10i \n", $n_atoms, $n_iter ;
	@atoms = () ;
        for ( my $i_atom=0; $i_atom<$n_atoms; $i_atom++ )
        {
	    $_ = <INPUT> ;
	    if (/velocity/) 
	    {
		$_ = <INPUT>;
		unless (/atom/) 
		{
			next ;
		}
	    }
	    my @line = split ' ', $_ ;
		push @atoms, [($line[4],$line[1],$line[2],$line[3])] ;
	}
	$n_iter++ ;
}

sub read_input_atoms
{
# Exact copy of read_atoms except for the push string order below which is different.
#	printf STDOUT "reading %5i atoms from iteration %10i \n", $n_atoms, $n_iter ;
	@atoms = () ;
        for ( my $i_atom=0; $i_atom<$n_atoms; $i_atom++ )
        {
	    $_ = <INPUT> ;
	    if (/velocity/) 
	    {
		$_ = <INPUT>;
		unless (/atom/) 
		{
			next ;
		}
	    }
	    my @line = split ' ', $_ ;
		push @atoms, [($line[3],$line[4],$line[5],$line[6])] ;
	}
	$n_iter++ ;
}

sub write_atoms
{
        for ( my $i_atom=0; $i_atom<$n_atoms; $i_atom++ )
	{
	    my $element = $atoms[$i_atom][0] ;
	    my $x = $atoms[$i_atom][1] ;
	    my $y = $atoms[$i_atom][2] ;
	    my $z = $atoms[$i_atom][3] ;
	    printf OUTPUT "atom %16.6f %16.6f %16.6f \t %s\n", $x, $y, $z, $element ;
	}
}

sub write_geometry
{
	# make hardcoded path platform independent
	my $geomfile = sprintf "geometry-%06i.in", $n_iter ;
	my $outfile = catfile $outdir, $geomfile ;
#	printf STDOUT "writing %s \n", $outfile ;
        open OUTPUT, '>', $outfile ;
	if ($n_lattice_vectors > 0)
	{
		&write_lattice_vectors ;
	}
	&write_atoms ;
	close OUTPUT ;
}

sub write_jmol_script
{	
	my $sptfile = "geometries.spt";
	open OUTPUT, '>', $sptfile ;
	print OUTPUT "load $zipfile \n" ;
	print OUTPUT "animation fps 1 \n" ;
	print OUTPUT "animation mode loop \n" ;
	print OUTPUT "animation on \n" ;
	close OUTPUT ;
}
