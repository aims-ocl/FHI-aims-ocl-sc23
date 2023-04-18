#!/bin/bash -eux

if [[ $# < 3 ]]
then
  echo "Please specify the install location, root directory for intermediate files, MPI launcher, and number of tasks to use for running tests."
  exit 1
fi

# The location where FHI-aims is installed
export AIMS_HOME=$1
# The root directory for the current job
export JOB_HOME=$2
# The MPI launcher to use
export MPIEXE=$3
# Number of tasks to use when running the FHI-aims regression tests
export N_TASKS_TEST=$4
# The terminal size
export TERMSIZE=200

echo "Running FHI-aims regression tests..."
cd "$JOB_HOME"
cd regression_tests
python3.3 ./regressiontools.py full --force --exclude="Symmetry-SPGlib" --exclude="Symmetry-RLSY" --mpiexe "$MPIEXE" --termsize=$TERMSIZE --cpus $N_TASKS_TEST --batch ./references_lastrelease/ "$AIMS_HOME/bin/aims.x"
echo "Finished running FHI-aims regression tests."
