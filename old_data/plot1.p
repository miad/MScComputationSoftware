#!/usr/bin/gnuplot -persist

set term wxt 0
set xlabel 'Distance /fm'
set ylabel 'Depth /MeV'
set title 'Potential function'
plot 'compute_output.dat' index 1 using 1:2 with lines title columnhead(1)

set term wxt 1
set xlabel 'Im(k) /fm^(-1)'
set ylabel 'Re(k) /fm^(-1)'
set title 'k-curves'
plot 'compute_output.dat' index 0 using 1:2 with points title columnheader(1),\
   'compute_output.dat' index 2 using 1:2 with points linecolor rgb "green" title columnheader(1)#,\
#	  'test4theory2.dat' using 1:2 with points linecolor rgb "blue" pointtype 6