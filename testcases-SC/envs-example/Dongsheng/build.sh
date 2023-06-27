#!/bin/bash

set -e

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <arg1> <arg2>"
  echo "<arg1> 0: build only, 1: cmake init and build"
  echo "<arg2> num of threads to build"
  exit 1
fi

source env.sh

if [ "$1" = "1" ]
then
    cmake -C ./initial_cache.cmake -S . -B build -DCMAKE_PREFIX_PATH="" -DOpenCL_FOUND=True -DOpenCL_LIBRARY=${MY_ROCM_PATH}/opencl/lib/ -DOpenCL_INCLUDE_DIR=${MY_ROCM_PATH}/opencl/include/
fi
cmake --build build -j $2
