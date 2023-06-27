#!/bin/bash
set -u

node_num=128
mpi_task_each_node=32
gpu_each_node=4
# 448

export MPI_PER_NODE=mpi_task_each_node
export GPU_PER_NODE=gpu_each_node

for count in {1..4}
do
    ((ntask = node_num * mpi_task_each_node))
    ((gpu_num = node_num * gpu_each_node))
    echo "count=${count}, node_num=${node_num}, ntask=${ntask}, gpu_num=${gpu_num}"
    filename="sbatch_node${node_num}_ntask${ntask}_gpu${gpu_num}.sh"
    dirname="sbatch_node${node_num}_ntask${ntask}_gpu${gpu_num}"

    if [ ! -d ${dirname} ]; then
        mkdir ${dirname}
    fi
    cd ${dirname}

    cp ../sbatch_template.sh ./sbatch_use.sh
    sed "s/__NODE_NUM__/${node_num}/g" sbatch_use.sh > tmp1.txt
    sed "s/__mpi_task_each_node__/${mpi_task_each_node}/g" tmp1.txt > tmp2.txt
    sed "s/__gpu_each_node__/${gpu_each_node}/g" tmp2.txt > ${filename}

    rm tmp1.txt
    rm ./sbatch_use.sh

    cp ../control.in .

    if [ ! -f geometry.in ]; then
        cp ../geometry.in .
    fi

    sbatch ${filename}
    
    cd ..

    ((node_num = node_num * 2))
done
