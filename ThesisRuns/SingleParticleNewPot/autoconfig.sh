#!/bin/bash
OIFS=$IFS
IFS=' '
read -a CB <<< "`cut -d ' ' -f 1 param.dat | awk '{printf("%s ",$0);}'`"
IFS=$OIFS

for i in ${CB[@]}
do 
	./config_generate.sh $i CONFDIR_$i &
done