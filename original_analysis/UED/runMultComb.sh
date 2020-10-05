#!/bin/bash -x


#FILEPREFIX="runFullBend"
#RUNFILE="runScanLegendre.sh"

FILEPREFIX="runList/MC/runCombineBend"
RUNFILE="combineScans.exe"
ind=1
NJOBS=40
until [ $ind -gt $NJOBS ]; do

  FILENAME=${FILEPREFIX}${ind}".txt"
  ./${RUNFILE} ${FILENAME}

  ind=$((ind + 1))
done
