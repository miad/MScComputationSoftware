#!/bin/bash
#Quick & dirty
killall -9 gnuplot 2>/dev/null

xmax=20
dirName="/net/data2/riklund/ThesisTwoparticle1Output/"
#dirName="/net/data2/riklund/LplusG2/"

values=("Ag8_ok.out")
titles=("Real part" "Imaginary part" "Overlap")
use=('u 0:1:(abs($3-1)*$1) with errorbars' 'u 0:2:(100*abs($3-1)*$2) with errorbars' 'u 0:3')
yra=('[0.5927:0.59285]' '[-2E-4:0]' '[0.98:1.01]')
for v in "${values[@]}"
do
	echo $dirName"lanc$v"
	if [ -f $dirName"lanc$v" ]
	then
		./overlap.sh $dirName"lanc$v" > lanc$v.plot
	else
		echo "Could not find LANC file lanc$v." 
	fi
done

LEN=`expr ${#values[@]} - 1`

for j in `seq 0 2`
do
	plotStr[j]="set xra [0:$xmax]\\nset yra ${yra[j]}\\nset title '${titles[j]}'\\nplot"
	
	for i in `seq 0 $LEN`
	do
		plotStr[j]=${plotStr[j]}"'lanc${values[i]}.plot' ${use[j]}"
		if [ "$i" -ne "`expr ${#values[@]} - 1`" ]
		then
			plotStr[j]=${plotStr[j]}",\\\\\\n"
			else
			plotStr[j]=${plotStr[j]}";\\n\\n"
		fi
	done
done



for j in `seq 0 1`
do
echo -en ${plotStr[j]}
echo -en ${plotStr[j]} | gnuplot -persistent
done