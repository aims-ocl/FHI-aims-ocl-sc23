#!/usr/bin/ksh
############################################################
#      a i t r a n s s : ab initio transport simulations
#     (c) 2003-2012    : alexej bagrets,  andreas arnold
#                        florian weigend, ferdinand evers
#     institute of nanotechnology (int) &
#     institut fuer theorie der kondensierten materie (tkm)
#     karlsruhe institute of technology (kit)
#
#     author:         alexej.bagrets <at> kit.edu
#     date:           nov 2009
#     last revision:  jan 2012
############################################################

# -------------------------------------------------------------
# "tcontrol.aims.x"  script creates a mandatory file "tcontrol"
#                    which is required to run the "aitranss" 
#                    post-processig module after FHI-aims
# -------------------------------------------------------------

# info message

info_message() 
{
cat <<EOF1

 ===============================================================

      a i t r a n s s :  ab initio transport simulations

      (c)  2003-2013  :  alexei bagrets,  andreas arnold.
                         florian weigend, ferdinand evers

      institute of nanotechnology (int) &
      institut fuer theorie der kondensierten materie (tkm)
      karlsruhe institute of technology (kit)

 ===============================================================

 When using this software, please cite the following references:

 (a) A. Arnold, F. Weigend, and F. Evers,
     J. Chem. Phys. 126, 174101 (2007).
 (b) J. Wilhelm, M. Walz, M. Stendel, A. Bagrets, and F. Evers,
     Phys. Chem. Chem. Phys. 15, 6684 (2013).
 (c) A. Bagrets, J. Chem. Theory Comput. 9, 2801 (2013).

 Report bugs to: alexej.bagrets <at> kit.edu

 --------------------------------------------------------------
 "tcontrol.aims.x"   script creates a mandatory file "tcontrol"
                     which is required to run the "aitranss"
                     post-processing module after FHI-aims
 --------------------------------------------------------------
EOF1
}

# HELP message
helpout()
{
cat <<EOF2

 USAGE: tcontrol.aims.x [ -option <argument> ] ...

 where options & arguments are: 
 
			! electrodes geometry:
                        
  -lsurc <atom1>        three atoms which define an outermost 
  -lsurx <atom2>        LEFT surface layer of the extended 
  -lsury <atom3>        molecule

  -rsurc <atom4>        three atoms which define an outermost 
  -rsurx <atom5>        RIGHT surface layer of the extended 
  -rsury <atom6>        molecule

  -nlayers <number>     number of the atomic layers coupled  
                        to reservoirs via a self-energy
 
			! energy window, in Hartree [H], to 
			! output transmission function T(E) :

  -ener  <E1[H]>        initial energy point, E1   
  -estep <dE[H]>        energy step, dE    
  -eend  <E2[H]>        final energy point, E2     

			! output :
  -outfile <file_name>  output file name for T(E) [ default: TE.dat ] 

EOF2
}

help_line()
{
  echo ' for HELP, type: tcontrol.aims -h [or -help, --help] '
  echo 
}

# initialization
init_parameters()
{
 list="1 2 3"
 for i in $list ; do 
  latoms[$i]=0 
  ratoms[$i]=0
  latoms_found[$i]="no" ; ratoms_found[$i]="no"
 done 
 isigma[1]=0.1d0     # default values in [Hartree]
 isigma[2]=0.05d0    # for the imaginary piece 
 isigma[3]=0.025d0   # of the self-energy  
                                             
 nlayers_found="no" ;  out_file="TE.dat" 
 ener="undefined"   ;  estep="undefined"  ;  eend="undefined"
 ener_found="no"    ;  estep_found="no"   ;  eend_found="no"
 eta=1.0d-10
}

# prints out input paramaters to tcontrol file
print_tcontrol()
{
 echo '$lsurc   '${latoms[1]} >> $tfile
 echo '$lsurx   '${latoms[2]} >> $tfile
 echo '$lsury   '${latoms[3]} >> $tfile
 echo '$rsurc   '${ratoms[1]} >> $tfile
 echo '$rsurx   '${ratoms[2]} >> $tfile
 echo '$rsury   '${ratoms[3]} >> $tfile
 echo '$nlayers '$nlayers >> $tfile
 
 echo '$s1i     '${isigma[1]} >> $tfile
 echo '$s2i     '${isigma[2]} >> $tfile
 echo '$s3i     '${isigma[3]} >> $tfile

 echo '$ener    '$ener >> $tfile
 echo '$estep   '$estep >> $tfile
 echo '$eend    '$eend >> $tfile
 echo '$output  file='$out_file >> $tfile
 
 echo '$testing off' >> $tfile
}

# checks for geometry.in file; counts number of atoms
get_atoms_info()
{
  crd_file="geometry.in" 
  if [[ ! -r $crd_file ]] ; then
   echo
   echo ' cannot find/read file '$crd_file
   echo ' please, check your directory content or access options '
   echo
   exit 1
  fi
  
  echo '$coord   file='$crd_file >> $tfile
  
  print -n ' analysing geometry.in file ... '
  awk ' BEGIN { count=0; }
	/#/  { next; }
	/\$/ { next; }
	{ 
	  if ( $1 == "atom" ) { count++; }
	  else { next; }
        }
       END { print "$natoms ", count}
      ' $crd_file >> $tfile
  echo 'DONE'  
}

# checks for basis-indices.out file
check_basis()
{
  basis_file="basis-indices.out"

  echo
  print -n ' searching for basis-indices.out ... '

  if [[ ! -r $basis_file ]] ; then
   echo
   echo
   echo ' cannot find/read file '$basis_file
   echo ' please, check your directory content or access options '
   echo
   exit 1
  fi
  echo '$basis   file='$basis_file >> $tfile
  echo 'FOUND'
}

# checks for overlap integrals: omat.aims file
check_overlap()
{
  omat_file="omat.aims"

  echo
  print -n ' searching for omat.aims ... '

  if [[ ! -r $omat_file ]] ; then
   echo
   echo
   echo ' cannot find/read file '$omat_file
   echo ' please, check your directory content or access options '
   echo
   exit 1
  fi
  echo '$read_omat file='$omat_file >> $tfile
  echo 'FOUND'
}

# checks for mos/alpha/beta files; write down 'nsaos' info 
get_mos_info()
{
 mosfile="mos.aims"
 afile="alpha.aims"
 bfile="beta.aims" 

  echo
  print -n ' analysing mos.aims/alpha.aims/beta.aims files ... '

# check whether mos.aims or alpha.aims & beta.aims exist
  if [[ -r $mosfile ]] ; then
    echo '$scfmo   file='$mosfile >> $tfile
    ifile=$mosfile
  else 
    if [[ ( -r $afile ) && ( -r $bfile ) ]] ; then
      echo '$uhfmo_alpha file='$afile >> $tfile
      echo '$uhfmo_beta  file='$bfile >> $tfile
      ifile=$afile
    else
      echo
      echo
      echo ' cannot find/read either "'$mosfile'" file or "'$afile'/'$bfile'" files'
      echo ' please, check your directory content or access options '
      echo
      exit 1
    fi 
  fi
  
# get info on matix dimension
  awk '/nsaos/ {print $0 ; exit }' $ifile | awk -F= '{print $3}' > nsaos.tmp
  read nsaos < nsaos.tmp ; rm nsaos.tmp 
  echo '$nsaos  ' $nsaos >> $tfile
  echo 'DONE'
}

# main body
if [[ $1 = "" ]] ; then
  info_message ; helpout ;  exit 1

else   
  tfile=tcontrol
  info_message ; init_parameters
  argcheck=y
  while [[ "$argcheck" = "y" ]] ; do
   if [[ -n "$1" ]] ; then
     case $1 in
      "-h" | "-help" | "--help" ) helpout ;  exit 1 ;;      
      "-lsurc"   )  latoms[1]=$2 ; latoms_found[1]="yes" ; shift 2 ;;  
      "-lsurx"   )  latoms[2]=$2 ; latoms_found[2]="yes" ; shift 2 ;;  
      "-lsury"   )  latoms[3]=$2 ; latoms_found[3]="yes" ; shift 2 ;;  
      "-rsurc"   )  ratoms[1]=$2 ; ratoms_found[1]="yes" ; shift 2 ;;  
      "-rsurx"   )  ratoms[2]=$2 ; ratoms_found[2]="yes" ; shift 2 ;;  
      "-rsury"   )  ratoms[3]=$2 ; ratoms_found[3]="yes" ; shift 2 ;;  
      "-nlayers" )  nlayers=$2 ; nlayers_found="yes" ; shift 2 ;;
      "-ener"    )  ener=$2  ; ener_found="yes"  ;  shift 2 ;;
      "-estep"   )  estep=$2 ; estep_found="yes" ;  shift 2 ;;
      "-eend"    )  eend=$2  ; eend_found="yes"  ;  shift 2 ;;
      "-outfile" )  out_file=$2   ; shift 2 ;;   
      * ) 
       echo 
       echo ' wrong option "'$1'" found!' ; help_line ; 
       argcheck=n; exit 1 ;;
     esac
   else 
    argcheck=n
   fi
  done


# check for left/right electrodes
  list="1 2 3"
  for i in $list ; do
   if [[ ${latoms_found[$i]} == "no" ]] ; then
    echo
    echo ' you have not specified LEFT electrode properly ' 
    echo 
    help_line ; exit 1 
   elif [[ ${ratoms_found[$i]} == "no" ]] ; then
    echo
    echo ' you have not specified RIGHT electrode properly '
    echo 
    help_line ; exit 1
   fi 
  done 

# check for atoms fixing surfaces being different ...
   if [[ ${latoms[1]} == ${latoms[2]} || ${latoms[2]} == ${latoms[3]} || ${latoms[3]} == ${latoms[1]} ]] ; then
    echo
    echo ' you have not specified LEFT electrode properly: '
    echo ' atoms fixing the surface plane should be different! ' 
    echo
    help_line ; exit 1
   fi
   if [[ ${ratoms[1]} == ${ratoms[2]} || ${ratoms[2]} == ${ratoms[3]} || ${ratoms[3]} == ${ratoms[1]} ]] ; then
    echo
    echo ' you have not specified RIGHT electrode properly: '
    echo ' atoms fixing the surface plane should be different! ' 
    echo 
    help_line ; exit 1
   fi

# check for coupling region 
  if [[ $nlayers_found == "no" ]] ; then
    echo
    echo ' you have not specified a coupling region'
    echo ' option "-nlayers <number>" should be used '
    echo 
    help_line ; exit 1
  fi
  
# prepair tcontrol-file
  echo '#input data for the "aitranss" module' > $tfile

# keyword for the negf implementaion: not actually used from jan 2012 
# echo '$non-equilibrium on' >> $tfile   

# keyword to check for aims-input:
  echo '$aims_input on' >> $tfile

# default keyword for the landauer tranmission
  echo '$landauer on' >> $tfile

  get_atoms_info ;  check_basis
  check_overlap  ;  get_mos_info  

# default: core states are integrated out
  echo '$ecp     on' >> $tfile

  print_tcontrol    
  
  echo '$end' >> $tfile
  echo
  echo ' done with <tcontrol> : ready to run transport module '
  echo
  exit 0
fi
# all done
