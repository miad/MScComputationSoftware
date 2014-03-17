#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 configdir outfile"
fi
CONFIG_DIR=$1
OUTFILE=$2
rm -f $OUTFILE
for i in $(ls $CONFIG_DIR | grep -e "config[0-9]*.conf" | grep -v ~); 
do
	RLO=`grep 'Name = "RLoffset"' -A 1 $CONFIG_DIR/$i | sed -n 's_\(\s*Value = \)\([0-9]*.[0-9]*\)\(;\)_\2_p'`
	EE=`cd ../.. && ./Compute --configFile ThesisRuns/SingleParticleLift/$CONFIG_DIR/$i | grep "Resonant state"`
	if [ "$?" -ne "0" ]
	then
		echo $CONFIG_DIR/$i "gives nonzero return."
		exit 1
	fi
	echo "RLO = $RLO gives $EE" >> $OUTFILE
done