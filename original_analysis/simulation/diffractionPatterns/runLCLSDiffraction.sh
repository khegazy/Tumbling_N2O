#!/bin/bash -x

OUTPUTFILE=n2oDiff
OUTPUTDIR=/reg/d/psdm/amo/amoi0314/scratch/simThetaPhi/diffractionPatterns

FILETORUN=n2oDiffSim.exe
CLUSTER=0
INDEX=1
if [ ! -z "$1" ]; then
  INDEX=${1}
fi

NMOLS=150000

SAMPLESTEP=0.010						#ps
STARTSAMPLETIME=357.045
SAMPLETIME=`echo $STARTSAMPLETIME + $INDEX \* $SAMPLESTEP |bc`			#ps


if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
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

#  cp /reg/neh/home/khegazy/simulations/n2o/diffractionPatterns/${FILETORUN} .
#  SAMPLETIME=`echo $STARTSAMPLETIME - 8 \* 39.756 + $INDEX \* $SAMPLESTEP |bc`
  SAMPLETIME=`echo $STARTSAMPLETIME + $INDEX \* $SAMPLESTEP |bc`
  OUTPUTFILE="job"${INDEX}"_"${SAMPLETIME}
  OUTPUTDIR="../"

else
  if [ ! -z "$1" ]; then
    OUTPUTFILE=${OUTPUTFILE}${1}
  fi
fi


#PDFADDR=`ls /reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/alignmentPDFs/UED/job${INDEX}_*`
PDFADDR=`ls /reg/d/psdm/amo/amoi0314/scratch/simThetaPhi/job${INDEX}_*`

echo "index: "${INDEX}
echo ${FILETORUN}" "${NMOLS}" -PDF "${PDFADDR}" -Ofile "${OUTPUTFILE}" -Odir "${OUTPUTDIR}" -Index "${INDEX}

./${FILETORUN} ${NMOLS} -PDF ${PDFADDR} -Ofile ${OUTPUTFILE} -Odir ${OUTPUTDIR} -Index ${INDEX} 


