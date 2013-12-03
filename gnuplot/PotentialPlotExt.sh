#!/bin/bash
#Syntax self-descriptive from below:
InputFile=$1
WFfile=$2
A=`head $WFfile -n1 | wc -w`
NCOL=$(((A-1)/2))
echo $NCOL
#plot '$InputFile' using 1:2 with lines linecolor rgb "red" title "Potential",\
#   for [i=2:$NCOL] '$WFfile' using 1:7 with lines linecolor rgb "blue" title "G-L points
gnuplot -persist << EOF
set term pdf
set xlabel 'Distance /fm'
set ylabel 'Depth /MeV'
set title 'Potential and wavefunctions'
set nokey
set output "results.pdf"
plot '$InputFile' using 1:2 with lines linecolor rgb "red" title "Potential",\
    for [i=4:$A:3] '$WFfile' using 1:i with lines linecolor rgb "blue" title "Wavefunctions"
EOF
