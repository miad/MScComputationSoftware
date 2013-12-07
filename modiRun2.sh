#!/bin/bash
DEPTH=$1
TAB=`echo -e "\t"`

mkdir -p runData

for (( DEPTH=-300; DEPTH<=100; DEPTH++ ))
do
	for (( PRECISION=0; PRECISION<=9; PRECISION++))
	do
		cp config_modirun2.conf config_modirun_run2.conf
		sed -i "s/DEPTH/$DEPTH.$PRECISION/g" config_modirun_run2.conf
		A=`./Compute --configFile config_modirun_run2.conf | grep state`
		B=${A//'Bound state:'/$TAB}
		B=${B//'Resonant state:'/$TAB}
		echo $DEPTH.$PRECISION $TAB $B
		cp -r output runData/output_$DEPTH\_$PRECISION/

	done
done