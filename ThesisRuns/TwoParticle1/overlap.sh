#!/bin/bash
if [ "$#" -ne "1" ]; then
	echo "Usage: $0 lanczos_outfile"
	exit 1
fi

sed -n -e 's_\(eigenvalues Lanczos (\)\([0-9\.]*\)\(,\)\([-0-9e\\.]*\)\() overlap (\)\(0.[5-9][-0-9e]*\|1[0-9\.]*\)\(,0)*\)_\2 \4 \6_p' $1