#!/bin/bash
#Syntax self-descriptive from below:
InputFile=$1
LengthUnit=$2


read -a hArr < $InputFile


NumberOfWF=`awk '{print NF}' $InputFile | sort -nu | tail -n 1`
#echo $NumberOfWF
PlotterString="plot "
for ((WFline=2; WFline<=$NumberOfWF; WFline+=2))
do
	if [ "$WFline" -gt "3" ]
	then
		PlotterString="$PlotterString,"
	fi
	hdr=${hArr[$(($WFline-1))]}
	hdr=${hdr##"Real_wave_$(($WFline / 2 - 1))("}
	hdr=${hdr%%")"}
	PlotterString="$PlotterString '$InputFile' using 1:(\$$WFline*\$$WFline+\$$(($WFline + 1))*\$$(($WFline + 1))) with lines title '$hdr'"
done

#echo $PlotterString

gnuplot -persist << EOF
set term wxt 0
set xlabel 'Distance /$LengthUnit'
set ylabel 'No unit'
set title 'Wavefunctions Modulus Squared'
set key
set auto
$PlotterString
replot
EOF
