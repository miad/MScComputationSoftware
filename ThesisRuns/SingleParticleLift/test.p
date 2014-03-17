#!/bin/bash
gnuplot -persist << EOF
plot './run_0.98989000000000000000.dat.parsed' u 1:2 lw 5 lt 1,\
'./run_0.99512000000000000000.dat.parsed' u 1:2 lw 5 lt 2,\
'./run_0.99806000000000000000.dat.parsed' u 1:2 lw 5 lt 3,\
'./run_0.99968000000000000000.dat.parsed' u 1:2 lw 5 lt 4,\
'./run_1.0000000000000000000.dat.parsed' u 1:2 lw 5 lt 5,\
'./run_1.0031100000000000000.dat.parsed' u 1:2 lw 5 lt 6,\
'./run_1.0035600000000000000.dat.parsed' u 1:2 lw 5 lt 7,\
'./run_1.0040700000000000000.dat.parsed' u 1:2 lw 5 lt 8,\
'./run_1.0045700000000000000.dat.parsed' u 1:2 lw 5 lt 9
EOF