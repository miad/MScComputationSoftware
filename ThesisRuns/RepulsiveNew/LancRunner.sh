#!/bin/bash
num=$(awk 'BEGIN{for(i=-0.1;i<=0.1;i+=0.01)print i}')
for i in $num
do
	echo "Running  with " $i
	timeout 4000 ../../../Lanczos/lanc /net/data2/riklund/TwoRepulsive/lancANP1_g$i.obj 160 149 | tee /net/data2/riklund/TwoRepulsive/lanc$i.out
done
