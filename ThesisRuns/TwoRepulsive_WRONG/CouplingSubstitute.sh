#!/bin/bash
num=$(awk 'BEGIN{for(i=-0.1;i<=0.1;i+=0.01)print i}')
for i in $num
do
	cp curveANP1.conf tempCurve.conf
	sed -i -e "s:@COUPLINGCOEFF@:$i:g" tempCurve.conf
	sed -i -e "s:@OUTEXT@:$i:g" tempCurve.conf
	echo "Running  with " $i
	#../../Compute --configFile tempCurve.conf
	#rm tempCurve.conf
done
