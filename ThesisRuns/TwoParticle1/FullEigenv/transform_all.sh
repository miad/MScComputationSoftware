#!/bin/bash
for i in `ls | grep 'eigenv[A-F].dat' | sed -n 's_\(eigenv\)\([A-F]\)\(\.dat\)_\2_p'`
do
	./transform.py "eigenv$i.dat" "kval$i.dat"
done