#!/bin/bash

aimspath=../../bin
V=$1
R=$2
differ=kdiff3
diffopts="-p -s"

if [ -z "$V" -o -z "$R" -o "$V" = "-h" -o "-V" = "--help" ]; then
    cat <<EOF
Usage: ./regression_wrapper.sh NEWSTEM OLDSTEM

You might want to create a copy of this script and adjust it to your
needs.
EOF
exit
fi

./regression.py -n 2 run $aimspath/aims.$V.scalapack.mpi.x  -r $R.scalapack.mpi
./regression.py  -s  run $aimspath/aims.$V.serial.x         -r $R.serial

./regression.py diff $aimspath/aims.$V.scalapack.mpi.x \
                --differ "$differ" -O "$diffopts" -r $R.scalapack.mpi
./regression.py diff $aimspath/aims.$V.serial.x \
                --differ "$differ" -O "$diffopts" -r $R.serial
