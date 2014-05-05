#!/bin/bash
if [ -f "values.txt" ]
then
	rm "values.txt"
fi
for i in $(ls config_files | grep -e "config[0-9]*.conf" | grep -v ~); 
do
	VRC=`grep 'Name = "VRightCut"' -A 1 config_files/$i | sed -n 's_\(\s*Value = \[\)\([0-9]*.[0-9]*\), \([0-9]*.[0-9]*\)\(\]\;\)_\2_p'`
	RLO=`grep 'Name = "RLoffset"' -A 1 config_files/$i | sed -n 's_\(\s*Value = \)\([0-9]*.[0-9]*\)\(;\)_\2_p'`
	xmin=`sed -n 's_\(\s*Xmin = \[\)\([0-9]*.[0-9]*\), \([0-9]*.[0-9]*\)\(\]\;\)_\2_p' config_files/$i`

	echo -e "$i : \t VRC = $VRC \t RLO = $RLO \t xmin = $xmin"
	echo -e "$VRC \t $RLO \t $xmin" >> values.txt
done