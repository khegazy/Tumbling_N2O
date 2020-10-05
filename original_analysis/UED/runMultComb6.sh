#!/bin/bash -x


#FILEPREFIX="runFullBend"
#RUNFILE="runScanLegendre.sh"

FILEPREFIX="runLists/MC/runCombineBend"
RUNFILE="combineScans.exe"

ind=0
FILENAME=${FILEPREFIX}${ind}".txt"
./${RUNFILE} ${FILENAME}

ind=30
NJOBS=40
until [ $ind -gt $NJOBS ]; do

  FILENAME=${FILEPREFIX}${ind}".txt"
  ./${RUNFILE} ${FILENAME}

  ind=$((ind + 5))
done
