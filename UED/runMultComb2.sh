#!/bin/bash -x


#FILEPREFIX="runFullBend"
#RUNFILE="runScanLegendre.sh"

FILEPREFIX="runLists/MC/runCombineBend"
RUNFILE="combineScans.exe"
ind=6
NJOBS=10
until [ $ind -gt $NJOBS ]; do

  FILENAME=${FILEPREFIX}${ind}".txt"
  ./${RUNFILE} ${FILENAME}

  ind=$((ind + 1))
done
