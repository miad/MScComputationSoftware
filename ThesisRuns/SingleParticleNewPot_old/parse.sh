#!/bin/bash
sed -n -e 's_\(cB = \)\([0-9]\.[0-9]*\)\(.*\)\(E = \[ +\)\([0-9]\.[0-9]\)\([0-9]\{2\}\)\([0-9]*\)\(.* decayRate= \)\([0-9]*\.[0-9]*\)\(.*\)_\2 \6.\7 \9_p' $1 | awk '{printf("$%.5f$ & $%.7f$ & $%.7f$ & $$ & $$ \\\\\\hline \n", $1, $2, $3)}'