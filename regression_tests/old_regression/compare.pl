#!/usr/bin/perl

print "\n    Basic check of aims output consistency compared to reference file.\n" ;
print "    This check does not replace a line-by-line diff!!!\n" ;

$new_file = @ARGV[0] ;
$reference_file = @ARGV[1] ;

 print "\n    Comparing $new_file vs. $reference_file:\n" ;

# For now, compare only: 
# * Number of MPI threads
# * Number of scf cycles
# * Final total energy
# * Final scf timing

# Read new data

open(INPUT, "<$new_file") ;

&get_data ;

$threads_inp = $threads ;
$etot_inp = $etot ;
$scf_inp = $scf ;
$sing1_inp = $sing1 ;
$sing2_inp = $sing2 ;
$sing3_inp = $sing3 ;
$esp1_inp = $esp1 ;
$esp2_inp = $esp2 ;
$esp3_inp = $esp3 ;
$RRMS_inp = $RRMS ;
$homo_inp = $homo ;
$lumo_inp = $lumo ;
$time_inp = $time ;
$complete_inp = $complete ;

close(INPUT) ;

# read reference data

open(INPUT, "<$reference_file") ;

&get_data ;

$threads_ref = $threads ;
$etot_ref = $etot ;
$scf_ref = $scf ;
$sing1_ref = $sing1 ;
$sing2_ref = $sing2 ;
$sing3_ref = $sing3 ;
$esp1_ref = $esp1 ;
$esp2_ref = $esp2 ;
$esp3_ref = $esp3 ;
$RRMS_ref = $RRMS ;
$homo_ref = $homo ;
$lumo_ref = $lumo ;
$time_ref = $time ;
$complete_ref = $complete ;

close(INPUT) ;

# Compare output in table format

if ($reference_file eq "CO2_TDDFT.reference.out") {

  print  "    |                     New aims version          Reference aims version\n" ;
  printf "    | Threads:            %16i                %16i\n", $threads_inp, $threads_ref ;
  printf "    | Total energy:       %16.6f eV             %16.6f eV\n", $etot_inp, $etot_ref ;
  printf "    | No. of scf cycles:  %16i                %16i\n", $scf_inp, $scf_ref ;
  printf "    | 1st Singlet:        %16.4f eV             %16.4f eV\n", $sing1_inp, $sing1_ref ;
  printf "    | 2nd Singlet:        %16.4f eV             %16.4f eV\n", $sing2_inp, $sing2_ref ;
  printf "    | 3rd Singlet:        %16.4f eV             %16.4f eV\n", $sing3_inp, $sing3_ref ;
  printf "    | Total time:         %16.3f s              %16.3f s\n", $time_inp, $time_ref ;
  printf "    | SCF cycle converged:%16s                %16s\n", $complete_inp, $complete_ref ;

} elsif ($reference_file eq "H2O_esp.reference.out" or $reference_file eq "H2O_esp2.reference.out") {

  print  "    |                     New aims version          Reference aims version\n" ;
  printf "    | Threads:            %16i                %16i\n", $threads_inp, $threads_ref ;
  printf "    | Total energy:       %16.6f eV             %16.6f eV\n", $etot_inp, $etot_ref ;
  printf "    | No. of scf cycles:  %16i                %16i\n", $scf_inp, $scf_ref ;
  printf "    | Total time:         %16.3f s              %16.3f s\n", $time_inp, $time_ref ;
  printf "    | SCF cycle converged:%16s                %16s\n", $complete_inp, $complete_ref ;
  printf "    | \n" ;
  printf "    | ESP charges, transition (HOMO-1) -> (LUMO)\n" ;
  printf "    | ESP charge O:        %16.4f e_0          %16.4f e_0\n", $esp1_inp, $esp1_ref ;
  printf "    | ESP charge H:        %16.4f e_0          %16.4f e_0\n", $esp2_inp, $esp2_ref ;
  printf "    | ESP charge H:        %16.4f e_0          %16.4f e_0\n", $esp3_inp, $esp3_ref ;
  printf "    | RMMS        :        %16.4f              %16.4f \n", $RRMS_inp, $RRMS_ref ;

} elsif ($reference_file eq "GaAs_SOC.reference.out") {

  print  "    |                     New aims version          Reference aims version\n" ;
  printf "    | Threads:            %16i                %16i\n", $threads_inp, $threads_ref ;
  printf "    | Total energy:       %16.6f eV             %16.6f eV\n", $etot_inp, $etot_ref ;
  printf "    | No. of scf cycles:  %16i                %16i\n", $scf_inp, $scf_ref ;
  printf "    | Total time:         %16.3f s              %16.3f s\n", $time_inp, $time_ref ;
  printf "    | SCF cycle converged:%16s                %16s\n", $complete_inp, $complete_ref ;
  printf "    | \n" ;
  printf "    | Spin-orbit-coupling perturbed quantities\n" ;
  printf "    | VBM/HOMO:           %16.6f eV           %16.6f eV\n", $homo_inp, $homo_ref ;
  printf "    | CBM/LUMO:           %16.6f eV           %16.6f eV\n", $lumo_inp, $lumo_ref ;

}  
else {

  #print  "\n    Aims output comparison:\n" ;
  print  "    |                     New aims version          Reference aims version\n" ;
  printf "    | Threads:            %16i                %16i\n", $threads_inp, $threads_ref ;
  printf "    | Total energy:       %16.6f eV             %16.6f eV\n", $etot_inp, $etot_ref ;
  printf "    | No. of scf cycles:  %16i                %16i\n", $scf_inp, $scf_ref ;
  printf "    | Total time:         %16.3f s              %16.3f s\n", $time_inp, $time_ref ;
  printf "    | SCF cycle converged:%16s                %16s\n", $complete_inp, $complete_ref ;

}

# end main perl file

sub get_data
{

    $threads = 1 ; # default
    $etot = 0. ;
    $scf = -1 ;
    $time = -1. ;
    $homo = 0.;
    $lumo = 0.;
    $complete = "yes" ;
    $found = 0 ;

while (<INPUT>)
{

  #  number of threads (if any)  
  if (/parallel tasks/)
  {
      @line = split " ", $_ ;
      $threads = @line[1] ;
  }

  # final total energy
  if (/Total energy uncorrected/)
  {
      @line = split " ", $_ ;
      $etot = @line[5] ;
  }

  # Number of scf cycles
  if (/Number of self\-consistency cycles/)
  {
      @line = split " ", $_ ;
      $scf = @line[6] ;
  }

  if (/ 1. Singlet:/)
  {
      @line = split " ", $_ ;
      $sing1 = @line[3];
  }
  if (/ 2. Singlet:/)
  {
      @line = split " ", $_ ;
      $sing2 = @line[3];
  }
  if (/ 3. Singlet:/)
  {
      @line = split " ", $_ ;
      $sing3 = @line[3];
  }
  if (/ Atom         1: O, ESP charge:/)
  {
      @line = split " ", $_ ;
      $esp1 = @line[6];
  }
  if (/ Atom         2: H, ESP charge:/)
  {
      @line = split " ", $_ ;
      $esp2 = @line[6];
  }  
  if (/ Atom         3: H, ESP charge:/)
  {
      @line = split " ", $_ ;
      $esp3 = @line[6];
  }
  if (/ RRMS:/)
  {
      @line = split " ", $_ ;
      $RRMS = @line[1];
  }
  if (/Highest occupied state \(VBM\) at/)
  {
      @line = split " ", $_ ;
      $homo = @line[5];
  }
  if (/Lowest unoccupied state \(CBM\) at/)
  {
      @line = split " ", $_ ;
      $lumo = @line[5];
  }

  # Total execution time
  if (/\| Total time   /)
  {
      @line = split " ", $_ ;
      $time = @line[4] ;
  }

  # Did the scf cycle converge?
  if (/WARNING\! SELF\-CONSISTENCY CYCLE DID NOT CONVERGE/)
  {
      $complete = " no" ;
  }
  if (/Leaving FHI\-aims/)
  {
      if ($complete=="yes")
      {      $found = 1 ; }

  }

}

if ($found==0)
{
    $complete = " no" ;
}

}

