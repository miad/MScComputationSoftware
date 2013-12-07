#!/bin/bash
#Syntax self-descriptive from below:
InputFileBasis=$1
InputFileEigen=$2
LengthUnit=$3

gnuplot -persist << EOF

set term wxt 1
set xlabel 'Im(k) /$LengthUnit^(-1)'
set ylabel 'Re(k) /$LengthUnit^(-1)'
set title 'k-curves'
set key right bottom
plot '$InputFileBasis' using 1:2 with points linecolor rgb "red" title "k-values in basis",\
   '$InputFileEigen' using 1:2 with points linecolor rgb "green" title "Eigenstates"
EOF