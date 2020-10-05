#!/bin/bash -x


FILEPREFIX="runFullBend"
RUNFILE="runScanLegendre.sh"

#FILEPREFIX="runCombineBend"

ind=0
NJOBS=40
until [ $ind -gt $NJOBS ]; do

  FILENAME=${FILEPREFIX}${ind}".txt"
  ./${RUNFILE} ${ind}

  ind=$((ind + 1))
done
