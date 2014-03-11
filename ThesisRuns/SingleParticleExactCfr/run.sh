#!/bin/bash
for i in $(ls | grep conf | grep -v ~); 
do
	CB=`grep 'Name = "cB"' -A 1 $i | sed -n 's_\(\s*Value = \[\)\([0-9].[0-9]*\), \([0-9].[0-9]*\)\(\]\;\)_\2_p'`
	EE=`cd ../.. && ./Compute --configFile ThesisRuns/SingleParticleExactCfr/$i | grep "Resonant state"`
	echo "cB = $CB gives $EE"
done