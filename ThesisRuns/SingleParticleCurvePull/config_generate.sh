#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 cB outdir"
	exit 1
fi

TMP_ARG=`mktemp`

CB=$1
OUTDIR=$2
./ArgsFinder.py $CB> $TMP_ARG
OIFS=$IFS
IFS=' '
read -a RLO <<< "`cut -d ' ' -f 1 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a RCUT <<< "`cut -d ' ' -f 2 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a XMIN <<< "`cut -d ' ' -f 3 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a POA <<< "`cut -d ' ' -f 4 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a POB <<< "`cut -d ' ' -f 5 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a POC <<< "`cut -d ' ' -f 6 $TMP_ARG | awk '{printf("%s ",$0);}'`"
read -a POD <<< "`cut -d ' ' -f 7 $TMP_ARG | awk '{printf("%s ",$0);}'`"


mkdir -p $OUTDIR

IFS=$OIFS
LEN=`expr ${#RLO[*]} - 1`
for i in `seq 0 $LEN`
do
	for j in `seq -f "%05g" 0 5000`
	do
		(
		name="$OUTDIR/config${i}_$j.conf";
		cp "config_templ.conf" $name;
		sed -i -e "s_!RLO!_${RLO[i]}_g" $name;
		sed -i -e "s_!RCUT!_${RCUT[i]}_g" $name;
		sed -i -e "s_!XMIN!_${XMIN[i]}_g" $name;
		sed -i -e "s_!POA!_${POA[i]}_g" $name;
		sed -i -e "s_!POB!_${POB[i]}_g" $name;
		sed -i -e "s_!POC!_${POC[i]}_g" $name;
		sed -i -e "s_!POD!_${POD[i]}_g" $name;
		sed -i -e "s@!CB!@$CB@g" $name;
		sed -i -e "s@!PULLED!@${j}.0E-6@g" $name;
		sed -i -e "s@!KCURVE!@${i}_${j}.dat@g" $name;
		sed -i -e "s@!KFOUND!@${i}_${j}.dat@g" $name;
		) &
	done
done


rm $TMP_ARG