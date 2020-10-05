#!/bin/bash

ROOTFOLDER=/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/rootFiles/
NFILES=1

ROOTFILES=$(ls ${ROOTFOLDER})
RUNLIST=runList.txt

count=0
list=1
for filename in $ROOTFILES; do
  if [ "$((${count} % ${NFILES}))" -eq "0" ]; then
    RUNLIST=runList${list}.txt
    count=0
    list=$((list + 1))

    if [ -e "$RUNLIST" ]; then
      rm $RUNLIST
    fi
  fi

  echo ${ROOTFOLDER}${filename} >> ${RUNLIST}
  count=$((count+1))
done
