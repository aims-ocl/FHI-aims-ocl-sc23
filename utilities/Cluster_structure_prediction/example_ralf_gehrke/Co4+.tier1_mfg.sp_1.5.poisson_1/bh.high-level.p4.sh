#! /bin/csh -f

## @ job_name = Cu8.tier1.ini1
# @ job_type = parallel
# @ output  = $(job_name).out.$(jobid)
# @ error   = $(job_name).err.$(jobid)
# @ initialdir = ./
# @ class = lhuge
# @ environment= COPY_ALL
# @ node_usage= shared
# @ node = 1
# @ tasks_per_node = 32
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ stack_limit = 0.5gb
# @ data_limit = 1.25gb

# @ notification =  complete
# @ notify_user = gehrke@fhi-berlin.mpg.de

# @ queue

# get fresh AFS token
get-token

# settings
set EXE_AIMS = '/afs/ipp/home/r/rgehrke/bin/aims.061507.essl.mpi.x'
set EXE_BH   = '/afs/ipp/home/r/rgehrke/bin/effernir.110607.x'
# end of settings ...

set CURRENT_DIR = ${PWD}
set TMPDIR=/ptmp/${USER}/${LOADL_JOB_NAME}
mkdir -p $TMPDIR

cp * $TMPDIR
cp $EXE_AIMS $TMPDIR
cp $EXE_BH $TMPDIR

cd $TMPDIR

date
time ./bh_optimization.pl $CURRENT_DIR $EXE_AIMS $EXE_BH
date
