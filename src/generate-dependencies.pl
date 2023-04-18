#!/usr/bin/perl

use strict;

# --- init

# global variables, read from *.d
my(@FILES, @MODS, @EXE, @USEFPP, %USES, %DEFS, @XPATTERNS);


# --------------------------------------------------------------------------- #
# ------------------------------- command line ------------------------------ #
# --------------------------------------------------------------------------- #

print "# Command line:\n";
print "# generate-dependencies.pl ", join(" ", @ARGV), "\n";
print "\n\n";

my @f90names = ();
my $modsuffix  = ".mod";
my $case;
while (my $opt = shift @ARGV) {
    if ($opt eq "--help" || $_ eq "-h") {
	print STDERR "Usage: $0 [options] f90names\n\n";
	print STDERR "  -u          uppercase module names\n";
	print STDERR "  -l          lowercase module names (default)\n";
	print STDERR "  -m suffix   suffix instead of \".mod\"\n";
	print STDERR "  -x pattern  exclude files matching pattern\n";
	print STDERR "\n";
	print STDERR "    The f90names are searched.\n";
	print STDERR "    Writes to stdout an .mk file to be included by\n";
	print STDERR "    GNU make.\n\n";
	exit 0;
    } elsif ($opt eq "-u") {
 	$case = "uc";
    } elsif ($opt eq "-l") {
 	$case = "lc";
    } elsif ($opt eq "-m") {
	$modsuffix = shift @ARGV;
    } elsif ($opt eq "-x") {
	push @XPATTERNS, shift @ARGV;
    } elsif ($opt =~ /^\-/) {
	die "Unknown option: $opt\n";
    } else {
	push @f90names, $opt;
    }
}
if (!@f90names) { die "Type \"$0 --help\" for help.\n"; }
if (!$case) { $case = "lc"; }


# --- parse files

foreach my $f90name (@f90names) {
    if (grep($f90name =~ /$_/, @XPATTERNS)) {
	next;
    }

    my $stem  = $f90name;
    $stem  =~ s/\.[^\.]+$//g;

    open(IN,"< $f90name") or die "Unable to open $f90name.\n";

    my(@defmods, @usemods, $exe, $fpp);
    while (<IN>) {
	if (/^\s*module\s+(\w+)/i)  {
	    my $defmod = $1;
	    if (grep $defmod eq $_, @usemods) {
		die "Module \"$defmod\" defined after usage.\n";
	    }
	    push @defmods, $defmod;
	}
	if (/^\s*use\s+(\w+)/i)     {
	    my $usemod = $1;
	    if (!grep $usemod eq $_, @defmods) {
		push @usemods, $usemod;
	    }
	}
	if (/^\s*program/i)         { $exe="yes"; }
	# line beginning with "#" means neccessity for preprocessing
	if (/^\#/)                  { $fpp="yes"; }
    }
    close(IN);

    # --- post edit

    if ( $case eq "uc" ) {
	@defmods = map(uc, @defmods);
	@usemods = map(uc, @usemods);
    } elsif ( $case eq "lc" ) {
	@defmods = map(lc, @defmods);
	@usemods = map(lc, @usemods);
    }	

    # ignore "module procedure" (out of interfaces)
    @defmods = grep(!/^procedure$/i, @defmods);

    # --- output

    push @FILES, $stem;
    if ($exe)     { push @EXE, $stem; }
    if ($fpp)     { push @USEFPP, $stem; }
    if (@defmods) { push @MODS, @defmods; }
    $DEFS{$stem} = \@defmods;
    $USES{$stem} = \@usemods;
}



# --------------------------------------------------------------------------- #
# ---------------------------------- parsing -------------------------------- #
# --------------------------------------------------------------------------- #

my @files;
my @usefpp;

my %srcuse = ();
my %srcdef = ();

@files = @FILES;
@usefpp  = @USEFPP;


# --- mod files for compiling

foreach my $file (@files) {

    my(@defmods, @usemods);
    # make modules right case
    @defmods = map($_ . $modsuffix, @{$DEFS{$file}});
    @usemods = map($_ . $modsuffix, @{$USES{$file}});

    $srcuse{$file} = [ @usemods ];
    $srcdef{$file} = [ @defmods ];

}

# --------------------------------------------------------------------------- #
# ---------------------------------- output --------------------------------- #
# --------------------------------------------------------------------------- #

# --- variables

if (@MODS) {
    print "# List of all .mod files\n";
    print "MODULE_HEADERS = ";
    foreach my $mod (@MODS) {
	print " \$(MODDIR)/$mod.mod";
    }
    print "\n\n";
}


# --- preprocess

if (@usefpp) {
    print "#  Now all f90 files which need a preprocessor are listed.\n";
    print "#  This could be used with a compile line like\n";
    print '#    $(FC) -c $($(@:.obj=)_FPPFLAGS) $(FFLAGS) $<', "\n";
    print "\n";
    foreach my $file (@usefpp) {
	print "${file}_FPPFLAGS = \$(FPPFLAGS)\n";
    }
    print "\n";
}

# --- <file>.mod: <file>.o ;

print "# Let all .mod files depend on the respective .o files "
    ,"with empty rules.\n";
foreach my $file (@files) {
    if (@{$srcdef{$file}}) {
	print join(" ", map("\$(MODDIR)/".$_,@{$srcdef{$file}})), " : ", "\$(OBJDIR)/${file}.o\n";
	foreach my $mod (@{$srcdef{$file}}) {
	    print "OBJECT_FILE_OF_$mod = $file\n";
	}
    }
}
print "\n";

print "# Provide variables that list an .o-file's contained modules\n";
foreach my $file (@files) {
    if (@{$srcdef{$file}}) {
	print "PROVIDED_MODFILES_$file = ", join(",", map("\$(MODDIR)/$_", @{$srcdef{$file}})), "\n";
    }
}

# --- <file>.o: <used_mods>.mod

print "# Let all .o files depend on the .mod files they use\n";
print "# either directly or indirectly.\n";
print "\n";
foreach my $file (@files) {
    if (@{$srcuse{$file}}) {
	print "\$(OBJDIR)/${file}.o  : ", join(" ", map("\$(MODDIR)/".$_,@{$srcuse{$file}})), "\n";
    }
}
print "\n";

