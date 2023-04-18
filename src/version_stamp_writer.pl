#!/usr/bin/perl

use strict;

# script to write a Fortran subroutine which writes the Version stamp and
# compile time.

# Obtain version stamp from input arg:
my $VERSION = $ARGV[0];

# obtain time and date from PERL function, and fix up to get proper format
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year = $year + 1900;
$mon  = $mon + 1;
my $date = sprintf("%4d/%02d/%02d", $year, $mon, $mday);
my $time = sprintf("%02d:%02d:%02d", $hour, $min, $sec);

# Where am I?
my $hostname = qx/hostname/;
chomp $hostname;
my $host = "host $hostname";

# Can we get the current Git ID?
my ($isgitrepos, $gitret, $ismod, $modstr, $gitrev, $gitstr);
my $perlret = system("git diff --raw --exit-code HEAD >/dev/null 2>&1");
if ($perlret == -1) {
    # execution failed quite early
    $isgitrepos = 0;
} else {
    print STDERR $gitret;
    $gitret = $perlret >> 8;
    $isgitrepos = $gitret == 0 || $gitret == 1;
}
if ($isgitrepos) {
    $ismod = $gitret == 1;
    if ($ismod) { $modstr = " (modified)"; } else { $modstr = ""; }
    $gitrev = qx/git show --no-color --pretty=oneline --abbrev-commit | head -n 1/;
    chomp $gitrev;
    # CC: clean git string from special chars that can break the compilation
    $gitrev =~ s/\@/[at]/g;
    $gitrev =~ s/[^A-Za-z 0-9\.,\[\]:\/-_]//g;
    if ( length($gitrev) > 50){
    	$gitrev = substr( $gitrev, 0, 50);
	$gitrev = $gitrev."[...]";
    };
    $gitstr = "Git rev.$modstr: $gitrev";
} elsif (-f "git-rev.txt") {
    my $catgitstr = qx/cat git-rev.txt/;
    chomp $catgitstr;
    $gitstr = "Based on $catgitstr";
}

if ($VERSION eq "--only-gitstr") {
    print "$gitstr\n";
} else {
    print "subroutine write_version_stamp()\n";
    print "  use localorb_io\n";
    print "  implicit none\n";
    print "  character*150 :: info_str\n";
    print "\n";
    print "  write(info_str,'(10X,A)') 'Version $VERSION'\n";
    print "  call localorb_info(info_str)\n";
    if ($gitstr) {
	print "  write(info_str,'(10X,A)') '$gitstr'\n";
	print "  call localorb_info(info_str)\n";
    }
    print "  write(info_str,'(10X,A)') 'Compiled on $date at $time on $host.'\n";
    print "  call localorb_info(info_str)\n";
    print "\n";
    print "end subroutine write_version_stamp\n";
}
