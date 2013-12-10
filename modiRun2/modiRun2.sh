#!/bin/bash
DEPTH=$1
TAB=`echo -e "\t"`

mkdir -p runData

for (( DEPTH=-300; DEPTH<=100; DEPTH++ ))
do
	for (( PRECISION1=0; PRECISION1<=9; PRECISION1++))
	do
		for (( PRECISION2=0; PRECISION2<=9; PRECISION2+=10))
		do
			cp config_modirun2.conf config_modirun_run2.conf
			sed -i "s/DEPTH/$DEPTH.$PRECISION1$PRECISION2/g" config_modirun_run2.conf
			A=`./Compute --configFile config_modirun_run2.conf | grep state`
			B=${A//'Bound state:'/$TAB}
			B=${B//'Resonant state:'/$TAB}
			echo $DEPTH.$PRECISION1$PRECISION2 $TAB $B
			cp -r /net/data2/riklund/MR1 /net/data2/riklund/output/output_$DEPTH\_$PRECISION1$PRECISION2/
			
		done
	done
done