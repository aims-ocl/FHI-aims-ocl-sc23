#!/usr/bin/awk -f

function abs(x) {
  return (((x < 0.0) ? -x : x) + 0.0)
}

BEGIN {
  print "# FHI-aims SCF convergence (aims-SCF_convergence.awk)"
  print "# JM, 2011/10/23, j.meyer@tum.de"
  print "#"
}

/accuracy/ && /charge density/ { A_rho = $NF }
/accuracy/ && /sum of eigenvalues/ { A_EV = $NF }
/accuracy/ && /total energy/ { A_E = $NF }
/accuracy/ && /forces/ { A_F =$NF }

/Finished reading input file 'control.in'/ {
  print "# Target accuracies:"
  printf "%-5s  %11s  %-11s  %-11s  %-11s  %-11s \n", "#", "", "A_rho", "A_EV (eV)", "A_E (eV)", "A_F (eV/A)"
  printf "%-5s  %11s  %11.4e  %11.4e  %11.4e  %11.4e \n", "#", "", A_rho, A_EV, A_E, A_F
  print "#"
}

/Begin self-consistency loop: Initialization./{
  printf "%5s  %-11s  %-11s  %-11s  %-11s  %-11s   ", "# SCF", "D_t", "D_rho", "D_EV (eV)", "D_E (eV)", "D_F (eV/A)"
  printf "accuracies reached \n"
  D_F = "-"
}

/Change of charge/ {D_rho = $NF }
/Change of charge/ && /spin density/ {if (NF > 0) D_rho = $(NF-1) }
/Change of sum of eigenvalues/ { D_EV = $(NF-1) }
/Change of total energy/  { D_E = $(NF-1) }
/Change of forces/ { D_F = $(NF-1) }

/End self-consistency iteration/ { 
  Nscf = $5
  getline
  D_t = $7
  if (D_F == "-")
    printf "%5d  %11.3f  %11.4e  %+11.4e  %+11.4e  [not calc.]   ", Nscf, D_t, D_rho, D_EV, D_E
  else
    printf "%5d  %11.3f  %11.4e  %+11.4e  %+11.4e  %+11.4e   ", Nscf, D_t, D_rho, D_EV, D_E, D_F
  if ( abs(D_rho) < A_rho ) printf "A_rho  "
  if ( abs(D_EV) < A_EV ) printf "A_EV  "
  if ( abs(D_E) < A_E ) printf "A_E  "
  if ( ( D_F != "-") && ( abs(D_F) < A_F ) ) printf "A_F  "
  print ""
}
