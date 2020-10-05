#!/bin/bash -x

OUTPUTFILE=legendreCoeffs

FILETORUN=scanLegendre.exe
INDEX=1
CLUSTER=0

DATAMC=1

if [ ! -z "$1" ]; then
  DATAMC=${1}
fi

if [ ! -z "$2" ]; then
  OUTPUTFILE=${2}
fi

if [ ! -z "$3" ]; then
  INDEX=${3}
fi

if [ ${DATAMC} -eq 0 ]; then
  RUNLIST=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/runLists/data/runList${INDEX}.txt
  OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/output/data/
elif [ ${DATAMC} -eq 1 ]; then
  RUNLIST=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/runLists/MC/runFull.txt
  OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/output/MC/
elif [ ${DATAMC} -eq 2 ]; then
  RUNLIST=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/runLists/MC/runFullBend${1}.txt
  OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/UED/output/MC/bend${1}/
else 
  echo "ERROR: MUST ENTER 0, 1, OR 2 FOR DATA/MC!!!"
  exit
fi

if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
fi 


#if [ ! -f ${FILETORUN} ]; then
#if [ "${OUTPUTFILE}" -eq "runoncluster" ]; then
if [ ${CLUSTER} -eq 1 ]; then
  echo "RUNNING ON CLUSTER"
  
  cp /reg/neh/home/khegazy/analysis/tumblingN2O/UED/${FILETORUN} .
  INDEX=${3}
  OUTPUTFILE="job"${INDEX}
  OUTPUTDIR="../"
fi


./${FILETORUN} ${RUNLIST} -Ofile ${OUTPUTFILE} -Odir ${OUTPUTDIR}


