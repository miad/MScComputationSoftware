#!/bin/bash
OIFS=$IFS
IFS=' '
read -a CB <<< "`cut -d ' ' -f 1 param.dat | awk '{printf("%s ",$0);}'`"
IFS=$OIFS

for i in ${CB[@]}
do 
    ./run.sh CONFDIR_$i run_$i.dat &
done