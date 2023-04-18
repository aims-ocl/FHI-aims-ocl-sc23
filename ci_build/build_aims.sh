#!/bin/bash -eux

if [[ $# < 3 ]]
then
  echo "Please specify the CMake toolchain file location, source directory, root directory for intermediate files, and number of tasks to use when building aims."
  exit 1
fi

# The CMake toolchain file, defining compilation/install settings for FHI-aims
export CMAKE_TOOLCHAIN_FILE=$1
# The root directory for the current job
export JOB_HOME=$2
# Number of tasks to use when compiling FHI-aims
export N_TASKS_BUILD=$3
# The location where FHI-aims should be installed
export AIMS_HOME=${JOB_HOME}/install

echo "Building FHI-aims..."
cat "$CMAKE_TOOLCHAIN_FILE"
mkdir -p "$JOB_HOME" && cd "$JOB_HOME"
rm -rf build install && mkdir build && cd build
cmake -C "$CMAKE_TOOLCHAIN_FILE" "$JOB_HOME"
cat 'CMakeCache.txt'
make -j $N_TASKS_BUILD
make install
echo "Finished building FHI-aims."
