#!/bin/bash
if [ -f "values.txt" ]
then
	rm "values.txt"
fi
for i in $(ls | grep -e "config[0-9].conf" | grep -v ~); 
do
	CB=`grep 'Name = "cB"' -A 1 $i | sed -n 's_\(\s*Value = \[\)\([0-9]*.[0-9]*\), \([0-9]*.[0-9]*\)\(\]\;\)_\2_p'`
	VRC=`grep 'Name = "VRightCut"' -A 1 $i | sed -n 's_\(\s*Value = \[\)\([0-9]*.[0-9]*\), \([0-9]*.[0-9]*\)\(\]\;\)_\2_p'`
	xmin=`sed -n 's_\(\s*Xmin = \[\)\([0-9]*.[0-9]*\), \([0-9]*.[0-9]*\)\(\]\;\)_\2_p' $i`
	omega=`sed -n 's_\(\s*AngularFrequency = \[\)\([0-9]*.[0-9]*\), *\([0-9]*.[0-9]*\)\(\]\;\)_\2_p' $i`

	echo -e "$i : \t cB = $CB \t VRightCut = $VRC \t xmin = $xmin \t omega = $omega"
	echo -e "$CB \t $VRC \t $xmin \t $omega" >> values.txt
done