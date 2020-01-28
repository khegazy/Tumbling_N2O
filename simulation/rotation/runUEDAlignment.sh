#!/bin/bash -x

OUTPUTFILE=anglePDFs
OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/alignmentPDFs/UED/

FILETORUN=simAlignment.exe
CLUSTER=1
INDEX=0
if [ ! -z "$1" ]; then
  INDEX=${1}
fi


#SAMPLESTEP=0.2575406
SAMPLESTEP=0.1
#STARTSAMPLETIME=315.914
#STARTSAMPLETIME=317.2726218
STARTSAMPLETIME=316.114
SAMPLETIME=`echo $STARTSAMPLETIME + $INDEX \* $SAMPLESTEP |bc`                    #ps


if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
fi 

if [ ! -z "$1" ]; then
  OUTPUTFILE=${1}
fi

#if [ ! -f ${FILETORUN} ]; then
#if [ "${OUTPUTFILE}" -eq "runoncluster" ]; then
if [ ${CLUSTER} -eq 1 ]; then
  echo "RUNNING ON CLUSTER"

  if [ ! -z "$1" ]; then
    INDEX=${1}
  else
    echo "ERROR: Must run cluster jobs with index argument!!!"
    exit
  fi

  cp /reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/${FILETORUN} .
  SAMPLETIME=`echo ${STARTSAMPLETIME} + $INDEX \* $SAMPLESTEP |bc`
  OUTPUTFILE="job"${INDEX}"_Temp-55_Intns-2.5_Time-"${SAMPLETIME}
  OUTPUTDIR="../"
fi


#./${FILETORUN} ${SAMPLETIME} -Ofile ${OUTPUTFILE} -Odir ${OUTPUTDIR} 
./${FILETORUN} -SampleTime ${SAMPLETIME}
