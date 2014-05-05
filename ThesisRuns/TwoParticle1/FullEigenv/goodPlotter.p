#!/bin/bash
killall -9 gnuplot 2>/dev/null

xmax=20

values=("A" "C" "D" "E" "F")

LEN=`expr ${#values[@]} - 1`

plotStr="set xra [0.59273:0.59289]\\nset yra [-1.5E-4:-0.5E-4]\\nset title 'Eigenvalues, letter indicating curve'\\nset xlabel 'Real part'\\nset ylabel 'Imaginary part' \\nplot "
for i in `seq 0 $LEN`
do
	plotStr=$plotStr"'./eigenv${values[i]}.dat' u 1:2"
	if [ "$i" -ne "`expr ${#values[@]} - 1`" ]
	then
		plotStr=$plotStr",\\\\\\n"
	else
		plotStr=$plotStr";\\n\\nreplot"
	fi
done


echo -en $plotStr
echo -en $plotStr | gnuplot -persistent
