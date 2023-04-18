#!/bin/bash

rm ../string.dat
k=-1
for i in paths/iteration*
do
    for j in $i/image*
    do
	cat $j | tail -n 1 >> ../string.dat
    done
    echo "" >> ../string.dat
    echo "" >> ../string.dat
    let k+=1
done
echo "j = $k" > ../plot.gnu
cat ../baseplot >> ../plot.gnu
