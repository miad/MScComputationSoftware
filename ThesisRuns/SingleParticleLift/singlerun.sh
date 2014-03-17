#!/bin/bash
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 configdir configfile outfile"
fi
CONFIG_DIR=$1
CONFIG_FILE=$2
OUTFILE=$3
HOSTDIR="/n/home/riklund/Dropbox/MastersProject_Drafts/ComputationSoftware/ThesisRuns/SingleParticleLift"
COMPUTEDIR="/n/home/riklund/Dropbox/MastersProject_Drafts/ComputationSoftware"

RLO=`grep 'Name = "RLoffset"' -A 1 "$HOSTDIR/$CONFIG_DIR/$CONFIG_FILE" | sed -n 's_\(\s*Value = \)\([0-9]*.[0-9]*\)\(;\)_\2_p'`


EE=`$COMPUTEDIR/./Compute --configFile $HOSTDIR/$CONFIG_DIR/$CONFIG_FILE | grep "Resonant state"`


echo "RLO = $RLO gives $EE" >> $HOSTDIR/$OUTFILE
