#!/bin/perl
# needed control.in.basic and geometry.in.basic
# assumption NVT_parrinello, tau 0.1

# IMPROVE declare all I/O files at the beginning

$BOLTZeVK = 8.617343e-5;
#parameters that might need adjustments
$INTEGRATOR = "NVT_parrinello";
$INTEGRATOR_TAU = 0.1;
$dirbase="rex_";
$plumed=0; # the default is PLUMED switched of; control.in.basic rules
srand(); # comment this line for reproducing the sequence

$status = @ARGV[0];
$logfile= @ARGV[1];
if (!$status) {
	print "Arg missing:, init, restart, postproc, stat\n";
  } 

if (!$logfile) {
	print "Missing log file specification (second argument)\n";
}

elsif ($status eq "init") {
  system "rm -rf $logfile ${dirbase}?? out.???? list_of_rex tt rex_par temp_list swap.info";
  # read control.in.rex 
  &read_rex;
  open(REXLLPAR,">rex_par");
  #print REXLLPAR "$n_rex\t$n_cpu\n";
  print REXLLPAR "$n_rex\t 0 \t $MAX_rex_steps\t";
  close(REXLLPAR);
  open(SW,">swap.info");
  print SW "1\n";
  close(SW);
  # create dirs and inputs
  if ( -e "geometry.in.basic" && -e "control.in.basic" ) {
        # check if PLUMED
	open(CONTRIN,"<control.in.basic");
	while(<CONTRIN>) {
		if (/plumed/) {
                @line = split ' ' , $_;
                if($line[0] ne '#' && $line[-1] eq '.true.') {
                        $plumed=1;      
                }
              }
	}
        close(CONTRIN);
   	for ($i=0; $i<$n_rex; $i++){
           $sfx=sprintf("%02d",$i);
           system "mkdir $dirbase$sfx || exit";
           open(CONTRIN,"<control.in.basic");
           open(CONTROUT,">$dirbase$sfx/control.in");
           while(<CONTRIN>) {
	      if (/MD\_run/ || /MD\_MB\_init/ || /MD\_schedule/ || /MD\_segment/ || /relax\_geometry/) {
		  #removing all unwanted lines
	      }
	      elsif (/MD_restart/ || /MD_clean_rotations/ || /output_level/) {
		  #removing all unwanted lines
	      }
	      else {
                 print CONTROUT $_;
              }
	   }
	   print CONTROUT <<ENDOFMOD;
\# All \"Physical model\", \"SCF convercenge\", \"Eigenvalue solution\", and basis set controls 
\# are left as in control.in.basic. Also \"MD_time_step\" is left.
\# Below are the MD controls imposed by the replica exchange scheme
  MD\_run $freq $INTEGRATOR $temp[$i+1] $INTEGRATOR_TAU
  MD\_MB\_init $temp[$i+1]
  MD\_restart\ \.true\.
  MD\_clean\_rotations\ \.true\.
  output\_level\ MD\_light
  MD_restart_binary .false.
ENDOFMOD
# The last option imposes an ascii restart file, in order to change velocities
#  MD\_restart\ \.false\.
           close(CONTRIN);
           close(CONTROUT);
           system "cp geometry.in.basic $dirbase$sfx/geometry.in";
  	}
	&write_target_T;
  } else {
        system "echo Exiting because of missing control\.in\.basic and\/or geometry\.in\.basic > EXIT";
        die;
  }

  if($plumed) {
	if(-e "plumed.basic") {
		for ($i=0; $i<$n_rex; $i++){
			$sfx=sprintf("%02d",$i);
			system "cp plumed.basic $dirbase$sfx/plumed.dat";
		}
	} else {
                system "echo Exiting because plumed\.basic is missing > EXIT";
                die;
        }
  }
  print REXLLPAR "$plumed\n";
  close(REXLLPAR);
  open(LG,">>$logfile");
  print LG "####### Starting a new Replica Exchange run  #################\n";
  if($plumed) { print LG "PLUMED is switched on\n \$plumed \= $plumed\n"; }
  close(LG);
}

elsif ($status eq "restart") {
  &read_rex;
  for ($i=0; $i<$n_rex; $i++) {
  	$sfx=sprintf("%02d",$i);
        #print "$dirbase$sfx\n";
        @neededfiles = ("control.in", "aims_MD_restart.dat", "geometry.in" );
	for ($j=0;$j<=$#neededfiles;$j++) {
		if (! -e "$dirbase$sfx/$neededfiles[$j]") {
			system "echo The existing $dirbase\?\? directories do not match control.in.rex > EXIT";
			die;
		}
	}
	if($plumed) {
		if (! -e "$dirbase$sfx/plumed.dat") {
                     	system "echo plumed.dat files are missing in some $dirbase\?\? directories > EXIT";
                        die;
		}
	}
   }
   open(LG,">>$logfile");
   print LG "####### Restarting a Replica Exchange run  #################\n";
   close(LG);
}

elsif ($status eq "postproc") {
  open(LG,">>$logfile");
  &read_rex;
  open(TEMPS,"<temp_list");
  $_ = <TEMPS> ;
  #while(<TEMPS>) {
  @line = split ' ', $_ ;
  for($i=1;$i<=$n_rex;$i++) {
  	$Tt[$i]=$line[$i];	# latest target Temperature, written in temp_list. Must agree with controls.in
  }
  $_ = <TEMPS> ;	
  @line = split ' ', $_ ;
  $plumed = $line[1];
  close(TEMPS);
  
  # list of E_final from outs
  for($i=0;$i<$n_rex;$i++) {	
  	$found = 0;
  	$conv = 1;
        $sfx=sprintf("%02d",$i);
        $tempout="$dirbase$sfx\/temp.out";
 	if (-e $tempout) {
  		open(AIMSOUT,"<$tempout");
		while(<AIMSOUT>) {
#               if(/Total energy \(el\.\+nuc\.\)/) {
                if(/Electronic\ free\ energy/) {
                        @line = split ' ', $_;
                        $TE[$i+1] = $line[5]; # Latest Total energy read from out
                }
           	if (/Have\ a\ nice\ day/) {
           		$found = 1 ;
           	}
           	if (/WARNING\!\ SELF\-CONSISTENCY\ CYCLE\ DID\ NOT\ CONVERGE/) {
           		$conv = 0 ;
           	}
		}
		close(AIMSOUT);
  		if (($found eq 0) || ($conv eq 0)) {
  			print LG " * WARNING: $tempout Not converged?\n" ;
       			print LG " Please check this problem before continuing. \n";
        		system "echo Exiting due to irregularities in $dirbase$sfx output > EXIT";
       			die;
		}
	}
	else {
  		print LG " * WARNING: $tempout Not converged?\n" ;
       		print LG " Please check this problem before continuing. \n";
        	system "echo Exiting due to irregularities in $dirbase$sfx output > EXIT";
       		die;
	}
  }
  # map $Tt[] wrt $temp[] from control.in.rex
  for($i=1;$i<=$n_rex;$i++) {	
  	for($j=1;$j<=$n_rex;$j++) {	
	   if($Tt[$i]==$temp[$j]) { $map[$j]=$i; } # chk for double entries?
	}
  }
  print LG "Tt\t @Tt\n","map\t @map\n";
  print LG "TE\t";
  for ($i=1;$i<=$n_rex;$i++) {printf LG "%8.4lf  ",$TE[$i];}
  print LG "\n";
  for ($i=1;$i<=$n_rex;$i++) { # the unchecked replicas (at boundaries) need initialization!
	$temp[$i] = $Tt[$i];
	$vfact[$i] = 1.;
  }
  open(SW,"<swap.info");
  $_ = <SW>;
  @line = split ' ',$_;
  $oneortwo = @line[0];
  close(SW);
  system "rm swap.info";
  open(SW,">swap.info");
  if ($oneortwo == 1) {
  	$start = 1;
	print LG "Start from one because swap.info contains $oneortwo\n";
	print SW "2\n";
  }
  else {
        $start = 2;
	print LG "Start from two because swap.info contains $oneortwo\n";
        print SW "1\n";

  }
  close(SW);
	
# useful for having rnd guess for starting T (1st or 2nd)
  #if(rand()<0.5) {
#	$start = 1;
#	#printf "Swap from 1st temperature \n";
#	} # neigh T starting from first
#  else {
#	$start = 2;
#	#printf "Swap from 2nd temperature \n";
#	} # neigh T starting from second
  for($i=$start;$i<=$n_rex-1;$i+=2) {
        $T2 = $Tt[$map[$i+1]];
        $T1 = $Tt[$map[$i]];
	$Boltz = ($TE[$map[$i+1]] - $TE[$map[$i]])*($T2 - $T1)/$T2/$T1/$BOLTZeVK; 
        print LG "swapping\t",$map[$i+1],"\t",$map[$i],"\t \@T \t",$temp[$map[$i+1]],"\t",$temp[$map[$i]],"\t";
	if ($Boltz < 0) { $accept=1; }
	elsif ($Boltz > 10) { $accept=0; }
  	elsif (rand() < exp(-$Boltz)) { $accept=1; }
	else { $accept=0; }
	if($accept) {
		$temp[$map[$i]] = $Tt[$map[$i+1]];
		$temp[$map[$i+1]] = $Tt[$map[$i]];
       		$vfact[$map[$i+1]] = sqrt ($temp[$map[$i+1]] / $temp[$map[$i]]);
		$vfact[$map[$i]] = 1./$vfact[$map[$i+1]];
                print LG "accepted\n";
	}
        else { 
                print LG "rejected\n";
        } 
   }

   &write_target_T;
   print LG "temp\t @temp\n","vfact\t @vfact\n";

   for ($i=0; $i<$n_rex ; $i++) {	
	$sfx=sprintf("%02d",$i);
        $tempout="$dirbase$sfx\/temp.out";
        open(CONTRIN,"<$dirbase$sfx\/control.in");
        open(CONTROUT,">$dirbase$sfx\/tmp");
        while(<CONTRIN>) {
        	if (/MD\_run/) {
        		@line = split ' ', $_ ;
        		$t = @line[1];
        		$t += $freq;
        		print CONTROUT "  MD\_run $t $line[2] $temp[$i+1] $line[4]\n";
        	}
        	elsif (/^  MD\_MB\_init/) {
        		print CONTROUT "\#$_";
        	}
                elsif (/MD\_restart\ \.false\./) {
			print CONTROUT "  MD\_restart\ \.true\.\n";
		}
        	else {
        		print CONTROUT $_;
        	}
        }
        close(CONTRIN);
        close(CONTROUT);
        system "mv $dirbase$sfx\/tmp $dirbase$sfx\/control.in";

	# scale velocities
	$k = $vfact[$i+1];
	# open and parse aims_MD_restart.dat (ASCII!!!)
	# rewrite all the lines except for velocities, v_half, v_last 
	open(RESTIN,"<$dirbase$sfx\/aims_MD_restart.dat");
	open(RESTOUT,">$dirbase$sfx\/tmp");
	while(<RESTIN>)	{
		if (/v_half/ || /v_last/) {
			@line = split ' ', $_ ;
			
			printf RESTOUT "%s               %30.20lf %30.20lf %30.20lf\n", $line[0],$line[1]*$k,$line[2]*$k,$line[3]*$k;
		}
		else {
			print RESTOUT $_;
		}
	} 
	close(RESTIN);
	close(RESTOUT);
	#system "mv $dirbase$sfx\/aims_MD_restart.dat $dirbase$sfx\/aims_MD_restart.dat.old";
	system "mv $dirbase$sfx\/tmp $dirbase$sfx\/aims_MD_restart.dat";

	$twrite = sprintf("%04d",$Tt[$i+1]);
        open(AIMSOUT,"<$tempout");
	open(ETRAJ,">>$dirbase$sfx/energy.trajectory");
	&energy_trajectory;
        close(ETRAJ);
        close(AIMSOUT);
        open(AIMSOUT,"<$tempout");
	open(XYZTRAJ,">>$dirbase$sfx/out.xyz");
        &jmolify;
	close(XYZTRAJ);
        close(AIMSOUT);
   	system "cat $tempout >> out.$twrite ";
	#print LG "Before copying COLVAR \$plumed \= $plumed\n"; 
	if($plumed) {
		system "cat $dirbase$sfx/COLVAR >> $dirbase$sfx/COLVAR.ark";
		system "cat $dirbase$sfx/COLVAR >> COLVAR.$twrite";
   		system "rm -f $dirbase$sfx/COLVAR";
	}
   	system "rm -f $tempout";
   }
   print LG "####### End of rex step #################\n";
   close(LG);
} 

elsif ($status eq "stat") {
  &read_rex;
  $n_steps = 0;
  $n_swap = 0;
  $n_accept = 0;
  for ($i=0;$i<$n_rex;$i++) {
   	$n_swapi[$i]=0;
   	$n_accepti[$i]=0;
  }
  if ( -e "$logfile") {
	open(LG,"<$logfile");  
        while(<LG>) {
  		$index1 = -1;
  		$index2 = -1;
		if (/End of rex step/) {
			$n_steps++;
                }
                elsif (/swapping/) {
			$n_swap++;
                        @line = split ' ', $_;
                        $T1 = $line[4];
                        $T2 = $line[5];
                        for ($i=1;$i<=$n_rex;$i++) {
                                if($temp[$i]==$T1) {$index1=$i;}
                                elsif($temp[$i]==$T2) {$index2=$i;}
                        }
 			if ($index1==-1 || $index2==-1 ) {
                                printf STDOUT "Error: one of the temperatures in $logfile does not correspond with control.in.rex\n";
   				exit;
                        }
                        $n_swapi[$index1]++;
                        $n_swapi[$index2]++;

			if (/accepted/) {
				$n_accept++;
				$n_accepti[$index1]++;
				$n_accepti[$index2]++;
			}
		}
	}
  }
  print STDOUT "########## Acceptance from $logfile #################\n";
  print STDOUT "Accepted swaps:\t",$n_accept,"\tAttempted swaps:\t",$n_swap,"\n";
  print STDOUT "Overall acceptance ratio:\t",$n_accept/$n_swap,"\n";
  for ($i=1;$i<=$n_rex;$i++) {
        if($n_swapi[$i]>0) {
		print STDOUT "Temperature\t",$temp[$i],"\tacceptance ratio\t",$n_accepti[$i]/$n_swapi[$i],"\n";
	}
  }
  print STDOUT "Total number of rex steps:\t",$n_steps,"\n";
}

else {
  print "$status is an unknown option\n";
  exit;
}

sub read_rex
{
  if ( -e "control.in.rex") { 	
        open(REX,"<control.in.rex");
  	while (<REX>) {
        	if (/n_rex/) {
           		@line = split ' ', $_;
           		$n_rex = $line[1];
        	}
        	elsif (/freq/) {
           		@line = split ' ', $_;
           		$freq = $line[1];
        	}
        	elsif (/temps/) {
           		@line = split ' ', $_;
           		if ($#line != $n_rex) {
              			print STDOUT "In control.in.rex there should be as many temperatures, after temps,\nas the numer of replicas declared in n_rex\n";
				system "echo Exiting due to irregularities in control.in.rex > EXIT";
              			exit;
           		}
           		for ($i=1; $i<=$n_rex ; $i++) {
              			$temp[$i] = $line[$i]; # temp array from 1 to n_rex
           		}
        	}
                elsif (/MAX_steps/) {
           		@line = split ' ', $_;
			$MAX_rex_steps = $line[1];
		}
  	}
	close(REX);
  }
  else {
       	system "echo Exiting because of missing control\.in\.rex > EXIT";
        die;
  }
}

sub write_target_T 
{
   open(TEMPS,">temp_list");
   print TEMPS "Target_temperatures\t";
   for ($i=1; $i<=$n_rex ; $i++) {
	print TEMPS $temp[$i],"\t"; 
   }
   print TEMPS "\nPLUMED $plumed\n";
   close(TEMPS);
}

sub energy_trajectory 
{
local($time,$TN,$FE,$KE,$TE,$NHH);
$time = "#      Time [ps]";
$TN    = "Temperature (nuclei)";
$FE   = "elec. Free Energy [eV]";
$KE   = "Nuclear kinetic energy";
$TE   = "Total Energy [eV]";
$NHH  = "Conserved Hamilt. [eV]";
printf ETRAJ "%15s %15s %23s %23s %23s %23s\n", $time, $TN, $FE, $KE, $TE, $NHH;
$iter = 0 ;
while (<AIMSOUT>)
{
    if (/Advancing\ structure\ using\ Born-Oppenheimer\ Molecular\ Dynamics/)
    {
	<AIMSOUT>;
	<AIMSOUT>;
	{
    		@line = split ' ', <AIMSOUT> ;
    		$time = @line[4] ;
    		@line = split ' ', <AIMSOUT> ;
    		$FE   = @line[5] ;
    		@line = split ' ', <AIMSOUT> ;
    		$T = @line[4] ;
    		@line = split ' ', <AIMSOUT> ;
    		$KE   = @line[5] ;
    		@line = split ' ', <AIMSOUT> ;
    		$TE   = @line[5] ;
    		$NHH = 0.0;
    		$_ = <AIMSOUT> ;
    		if(/BDP\ pseudo\-Hamiltonian/) {
        		@line = split ' ', $_;
        		$NHH =  @line[4] ;
    		}
    		else {
      			$_ = <AIMSOUT> ;
      			$_ = <AIMSOUT> ;
      			if(/Nose-Hoover\ Hamiltonian/) {
        			@line = split ' ', $_;
        			$NHH =  @line[4] ;
      			}
    		}
    		printf ETRAJ "%15.6f %15.6f %23.6f %23.6f %23.6f %23.6f\n", $time, $T, $FE, $KE, $TE, $NHH;
	}

    }
	
}
}

sub jmolify
{
$iter = 1 ;
while (<AIMSOUT>)
{

    if (/\|\ Number\ of\ atoms/)
    {
	@line = split ' ', $_ ;
        $n_atomsxyz = @line[5] ;
    }
    elsif (/\ Updated\ atomic\ structure\:/)
    {
        <AIMSOUT>;
	&print_new_geo ;
    }
    elsif (/\ Final\ atomic\ structure\:/)
    {
        <AIMSOUT>;
	&print_new_geo ;
    }
    elsif (/used\ in\ the\ preceding\ time\ step/)
    {
        <AIMSOUT>;
	&print_new_geo ;
    }
}

}

sub print_new_geo
{
	printf XYZTRAJ "%12i\n",$n_atomsxyz ;
        printf XYZTRAJ "Iteration:%6i\n", $iter ;
        for ($i_atom=1;$i_atom<=$n_atomsxyz;$i_atom++)
        {
	    $_ = <AIMSOUT> ;
	    if (/velocity/) 
	    {
		$_ = <AIMSOUT>;
		if (/atom/) 
		{
		    @line = split ' ', $_ ;
		    printf XYZTRAJ " %s %16.6f %16.6f %16.6f\n",@line[4],@line[1],@line[2],@line[3] ;
		}
	    }
	    else
	    {
		@line = split ' ', $_ ;
		printf XYZTRAJ " %s %16.6f %16.6f %16.6f\n",@line[4],@line[1],@line[2],@line[3] ;
	    }
        }
        $iter++ ;

}
