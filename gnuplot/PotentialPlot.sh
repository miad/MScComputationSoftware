#!/bin/bash
#Syntax self-descriptive from below:
InputFile=$1
PrecFile=$2
LengthUnit=$3
DepthUnit=$4

gnuplot -persist << EOF
set term wxt 0
set xlabel 'Distance /$LengthUnit'
set ylabel 'Depth /$DepthUnit'
set title 'Potential function'
set nokey
plot '$InputFile' using 1:2 with lines linecolor rgb "red" title "Potential",\
   '$PrecFile' using 1:2 with points linecolor rgb "blue" title "G-L points"
EOF
