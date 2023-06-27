#!/bin/bash

set -e

if [ "$1" = "1" ]
then
    cmake -C ./initial_cache.cmake -S . -B build -DCMAKE_PREFIX_PATH="$HOME/project/install;$HOME/project/install" -DOpenCL_FOUND=True -DOpenCL_LIBRARY=/opt/intel/oneapi/compiler/2023.1.0/linux/lib/ -DOpenCL_INCLUDE_DIR=/opt/intel/oneapi/compiler/2023.1.0/linux/include/sycl/
fi
cmake --build build -j $2