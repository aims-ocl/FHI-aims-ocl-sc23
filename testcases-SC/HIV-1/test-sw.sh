#!/bin/bash

filename="../../build/aims.191127.scalapack.mpi.x"

node_num=1

export FHI_AIMS_OCL_RHO_MULTIPOLE_TILE_SIZE=512
export FHI_AIMS_OCL_RHO_MULTIPOLE_SHMEM=OFF
export FHI_AIMS_OCL=ON
export FHI_AIMS_OCL_DFPT_POLAR_RHO=ON
export FHI_AIMS_OCL_DFPT_DIEL_SUMUP=ON
export FHI_AIMS_OCL_DFPT_POLAR_H=ON
export MPI_PER_NODE=1
export GPU_PER_NODE=1
export FHI_AIMS_OCL_DEBUG_IO=ON

# cur_timestamp=$(date '+%Y-%m-%d-%H-%M-%S')
currentTime=`date "+%Y-%m-%d-%H-%M-%S"`

ofname="log.athread.sumup.${currentTime}.txt"

echo ${currentTime}
cat ./test.sh >> ${ofname}
echo "=======================================================" >> ${ofname}
echo "currentTime=${currentTime}" >> ${ofname}
echo "node_num=${node_num}" >> ${ofname}
echo "=======================================================" >> ${ofname}
ls ${filename} -alh >> ${ofname}

bsub -b -m 1 -J aims-c2h43h-athread-sumup -n ${node_num} -o ${ofname} -cgsp 64 -share_size 13000 -host_stack 1024 -priv_size 16 ${filename}
echo "=======================================================" >> ${ofname}

# watch tail -n 20 ${ofname}
# -------------------- end of script --------------------



