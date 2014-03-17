#!/bin/bash
OIFS=$IFS
IFS=' '
read -a CB <<< "`cut -d ' ' -f 1 param.dat | awk '{printf("%s ",$0);}'`"
read -a XMIN <<< "`cut -d ' ' -f 2 param.dat | awk '{printf("%s ",$0);}'`"
read -a OMEGA <<< "`cut -d ' ' -f 3 param.dat | awk '{printf("%s ",$0);}'`"
read -a RCUT <<< "`cut -d ' ' -f 4 param.dat | awk '{printf("%s ",$0);}'`"
IFS=$OIFS
LEN=`expr ${#CB[*]} - 1`
for i in `seq 0 $LEN`
do
	cp "config_templ.conf" "config$i.conf"
	sed -i -e "s_!CB!_${CB[i]}_g" "config$i.conf"
	sed -i -e "s_!XMIN!_${XMIN[i]}_g" "config$i.conf"
	sed -i -e "s_!OMEGA!_${OMEGA[i]}_g" "config$i.conf"
	sed -i -e "s_!RCUT!_${RCUT[i]}_g" "config$i.conf"
done