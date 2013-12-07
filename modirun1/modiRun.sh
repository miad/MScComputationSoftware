#!/bin/bash
DEPTH=$1
TAB=`echo -e "\t"`

for (( DEPTH=0; DEPTH<=100; DEPTH++ ))
do
	for (( PRECISION=0; PRECISION<=9; PRECISION++))
	do
		cp config_modirun.conf config_modirun_run.conf
		sed -i "s/DEPTH/$DEPTH.$PRECISION/g" config_modirun_run.conf
		A=`./Compute --configFile config_modirun_run.conf | grep state`
		B=${A//'Bound state:'/$TAB}
		B=${B//'Resonant state:'/$TAB}
		echo $DEPTH.$PRECISION $TAB $B >> result.dat
		echo $DEPTH.$PRECISION $TAB $B
	done
done