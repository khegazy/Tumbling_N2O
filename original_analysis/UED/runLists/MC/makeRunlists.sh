#!/bin/bash

ROOTFOLDER=/reg/neh/home5/khegazy/simulations/n2o/diffractionPatterns/output/diffraction/
ROOTSEARCH=${ROOTFOLDER}"*.root"
NFILES=1

ROOTFILES=`ls ${ROOTSEARCH}`
#ROOTFILES=$(ls ${ROOTFOLDER}"*.root")
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

  echo ${filename} >> ${RUNLIST}
  count=$((count+1))
done
