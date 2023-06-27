#!/bin/bash

set -e

source /home/export/base/shisuan/swyjs/lxh/fhi-aims_CPE_O3_best-base-unopt/env.sh 

if [ "$1" = "1" ]
then
    cmake -C ./initial_cache.cmake -S . -B build2 -DCMAKE_PREFIX_PATH="" -DOpenCL_FOUND=True -DOpenCL_LIBRARY=/home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl -DOpenCL_INCLUDE_DIR=/home/export/base/shisuan/swyjs/cq/swcl/install/lib/swcl/host/include
fi
# cmake --build build -j $2