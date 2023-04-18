#!/usr/bin/perl

use strict;

# this PERL script generates a fortran subroutine which 
# outputs the version stamp and a compile time: 
# is based on extended script in the origial FHI-aims distribution

# obtain version stamp from input arg:
  my $VERSION = $ARGV[0];

# obtain time and date from PERL function, and fix up to get proper format
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);  
  $year = $year + 1900;
  $mon  = $mon + 1;
  my $date = sprintf("%4d/%02d/%02d", $year, $mon, $mday);
  my $time = sprintf("%02d:%02d:%02d", $hour, $min, $sec);

# where am i?
  my $hostname = qx/hostname/;
  chomp $hostname;
  my $host = "host $hostname";

  print "      subroutine output_version_stamp\n";
  print "       implicit none\n";
  print "\n";
  print "       write(6,fmt='(/,1x,a)') 'invoking aitranss (version: $VERSION) ... '\n";
  print "       write(6,fmt='(1x,a)') 'compiled on $date at $time on $host'\n";
  print "\n";
  print "      end subroutine output_version_stamp\n";

