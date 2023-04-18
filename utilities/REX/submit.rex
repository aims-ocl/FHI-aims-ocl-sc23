#!/bin/bash -l

#$ -S /bin/bash
### join stdout and stderr
#$ -j y
### change to current working dir
#$ -cwd
### send no mail
#$ -m n
### my email address
#$ -M my@email
### Parallel Environment
#$ -pe impi 96
### wallclock, e.g. 3600 seconds
#$ -l h_rt=24:00:00
### virtual memory (45G is max. on AIMS)
#$ -l h_vmem=22G

module load impi
module load mkl

################### attention of the user required here ###############

binary='aims.041010.meta.scalapack.shm.mpi.x'

# put  type='init', if initializing, 'restart' if restarting
type='init'
#type='restart'

# number of CPU per node (host)
ncpupn=12

#######################################################################


# some variables
hostfile='host.all'
list_of_hosts='hostlist'
perl_rex='rex.AIMS.pl'
logfile='log_rex'
infofile='rex_par'
touch $logfile

rm -f $hostfile rex_*/$list_of_hosts EXIT || exit

# thnec cluster
cat $TMPDIR/machines_$USER > $hostfile
# AIMS cluster
#cat $TMPDIR/machines > $hostfile

if [ $type != "init" -a $type != "restart" ]
then
  echo "-----------------------------------" >> $logfile
  echo "Unknown mode option, only init and restart are allowed" >> $logfile
  echo "-----------------------------------" >> $logfile
  exit 1
fi

perl $perl_rex $type $logfile
if [ -f EXIT ]
then
  cat EXIT >> $logfile
  exit
fi


# other variables after $perl_rex has run
n_of_replicas=`awk '{print $1}' $infofile`
n_of_sl4rep=$(($NSLOTS/$n_of_replicas))
rex_step=`awk '{print $2}' $infofile`
MAX_rex_step=`awk '{print $3}' $infofile`
n_of_hosts=`wc -l $hostfile | awk '{ print $1 }'`

# dummy check
# check if $NSLOTS is divisible by $ncpupn
if [ $(( $NSLOTS % $ncpupn )) -ne 0 ]; then
  echo "-----------------------------------" >> $logfile
  echo "Error! The requested number of slots (specified in job-submit-script) is not a multiple of the number of CPU per node (specified in job-submit-script)." >> $logfile
  echo "Aborting." >> $logfile
  echo "-----------------------------------" >> $logfile
  exit 1
fi

# check if $NSLOTS is greater (or equal) than $n_of_replicas
if [ $NSLOTS -lt $n_of_replicas ]; then
  echo "-----------------------------------" >> $logfile
  echo "Error! The requested number of slots (specified in job-submit-script) must be greater (or equal) than the number of replicas (also specified in control.in.rex)." >> $logfile
  echo "Aborting." >> $logfile
  echo "-----------------------------------" >> $logfile
  exit 1
fi

# check if $NSLOTS is divisible by $n_of_replicas
if [ $(( $NSLOTS % $n_of_replicas )) -eq 0 ]; then
  echo "-----------------------------------" >> $logfile
  echo "Info: The requested number of slots is a multiple of the number of replicas. That's perfect!" >> $logfile
  echo "-----------------------------------" >> $logfile
else
  echo "-----------------------------------" >> $logfile
  echo "Warning! The requested number of slots (specified in job-submit-script) is not a multiple of the number of replicas (specified in control.in.rex)." >> $logfile
  echo "That's not a problem per se. However, it might result in performance issues, because some slots (cores) will not be used." >> $logfile
  echo "In particular, $(( $NSLOTS % $n_of_replicas )) of $NSLOTS slots will not be used." >> $logfile
  echo "-----------------------------------" >> $logfile
fi

# check if $n_of_sl4rep is divisible by $ncpupn or vice versa
if [ $(( $n_of_sl4rep % $ncpupn )) -eq 0 ]; then
  echo "-----------------------------------" >> $logfile
  echo "Info: The number of slots per replica is a multiple of the number of CPU per node. That's perfect!" >> $logfile
  echo "-----------------------------------" >> $logfile
elif [ $(( $ncpupn % $n_of_sl4rep )) -eq 0 ]; then
  echo "-----------------------------------" >> $logfile
  echo "Info: The number of CPU per node is a multiple of the number of slots per replica. That's good." >> $logfile
  echo "-----------------------------------" >> $logfile
else
  echo "-----------------------------------" >> $logfile
  echo "Warning! The number of slots per replica is not a multiple of the number of CPU per node AND NEITHER is the number of CPU per node a multiple of the number of slots per replica." >> $logfile
  echo "That's not a problem per se. However, it might result in performance issues, because some parallel jobs will run across multiple nodes (hosts)." >> $logfile
  echo "-----------------------------------" >> $logfile
fi

# print some information to $logfile
echo "-----------------------------------" >> $logfile
echo "Number of slots: " $NSLOTS >> $logfile
echo "Number of replicas: " $n_of_replicas >> $logfile
echo "Number of slots per replica: " $n_of_sl4rep >> $logfile
echo "Number of unused slots: " $(( $NSLOTS % $n_of_replicas )) >> $logfile
echo "Number of hosts (nodes): " $n_of_hosts >> $logfile
echo "Number of CPU per node [hard-coded(!) in job-script (ncpupn)]: " $ncpupn >> $logfile
echo "-----------------------------------" >> $logfile

if [ -f list_of_geometries ]
then
  c=0
  cat list_of_geometries | while read j; do
        suffix=`printf "%02d" $c`
        cp $j rex_$suffix/geometry.in
        c=$(($c+1))
  done
fi

c=0
while [ $c -lt $n_of_hosts ]; do
  head -$(($c+1)) $hostfile | tail -1 > $hostfile.host`printf "%03d" $c`
  c=$(($c+1))
done

c=0
d=0
while [ $c -lt $NSLOTS ]; do
  if [ $(( $c % $ncpupn )) -eq 0 ] && [ $c -ne 0 ]; then
    d=$(($d+1))
  fi
  slot_suffix=`printf "%03d" $c`
  host_suffix=`printf "%03d" $d`
  cat $hostfile.host$host_suffix > $hostfile.host$host_suffix.slot$slot_suffix
  c=$(($c+1))
done

c=0
d=0
while [ $c -lt $NSLOTS ]; do
  if [ $(( $c % $n_of_sl4rep )) -eq 0 ] && [ $c -ne 0 ]; then
    d=$(($d+1))
  fi
  if [ $d -eq $n_of_replicas ]; then break; fi
  slot_suffix=`printf "%03d" $c`
  suffix=`printf "%02d" $d`
  if [ ! -f rex_$suffix/$list_of_hosts ]; then touch rex_$suffix/$list_of_hosts; fi
  cat $hostfile.host*.slot$slot_suffix >> rex_$suffix/$list_of_hosts
  c=$(($c+1))
done

rm $hostfile*

while true; do
  c=0
  while [ $c -lt ${n_of_replicas} ]; do
    suffix=`printf "%02d" $c`
    (cd rex_$suffix &&
      mpiexec -machinefile $list_of_hosts -n $n_of_sl4rep $binary >> temp.out) &
  c=$(($c+1))
  done
  wait
  date >> $logfile
  perl $perl_rex postproc $logfile || exit
  rex_step=$(($rex_step+1)) 
  echo $n_of_replicas $rex_step $MAX_rex_step > $infofile
  if [ $rex_step -ge $MAX_rex_step ]; then
   echo "The desired number of Replica Exchange steps has been performed, exiting" >> $logfile
   exit
  fi
  if [ -f EXIT ]; then 
    cat EXIT >> $logfile
    exit
  fi
done

