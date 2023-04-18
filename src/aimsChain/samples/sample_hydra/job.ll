#!/bin/bash
#
# Sample script for LoadLeveler
#
# @ job_name            = chain
# @ error   = job1.err.$(jobid)
# @ output  = job1.out.$(jobid)
# @ job_type = parallel
# @ environment= COPY_ALL
# @ node_usage= not_shared
# @ node = 1
# @ tasks_per_node = 16
# @ resources = ConsumableCpus(1)
# @ network.MPI = sn_all,not_shared,us
# @ wall_clock_limit    = 24:00:00
# @ notification = complete
# @ queue

#
# run the program
#
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_HOME/lib/intel64/

runchain.py
