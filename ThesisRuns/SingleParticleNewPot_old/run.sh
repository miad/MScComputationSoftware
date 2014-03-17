#!/bin/bash
for i in $(ls | grep -e "config[0-9].conf" | grep -v ~); 
do
	CB=`grep 'Name = "cB"' -A 1 $i | sed -n 's_\(\s*Value = \[\)\([0-9].[0-9]*\), \([0-9].[0-9]*\)\(\]\;\)_\2_p'`
	EE=`cd ../.. && ./Compute --configFile ThesisRuns/SingleParticleNewPot/$i | grep "Resonant state"`
	echo "cB = $CB gives $EE"
done