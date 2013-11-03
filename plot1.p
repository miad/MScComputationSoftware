#!/usr/bin/gnuplot -persist
set term wxt 0
plot 'compute_output.dat' index 1 using 1:2 with lines title columnheader(1)

set term wxt 1
plot 'compute_output.dat' index 0 using 1:2 with points title columnheader(1),\
   'compute_output.dat' index 2 using 1:2 with points linecolor rgb "green" title columnheader(1)