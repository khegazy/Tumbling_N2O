#!/bin/bash 

if [ -z "$2" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give the number of jobs and optional output directory"
  exit
fi

FILETORUN=run${2}Alignment.sh
EXECUTABLE=simAlignment.exe
MAINDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation

NJOBS=${1}
OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/clusterJobs/${EXP}/
#/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/alignmentPDFs/UED/plotJobs/

#make clean; make

if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
  mkdir ${OUTPUTDIR}/logs
  mkdir ${OUTPUTDIR}/run
else
  if [ ! -d "${OUTPUTDIR}/logs" ]; then
    echo "HERE"
    mkdir ${OUTPUTDIR}/logs
    mkdir ${OUTPUTDIR}/run
  fi
fi

cd ${OUTPUTDIR}/run
cp ${MAINDIR}/${FILETORUN} .
cp ${MAINDIR}/${EXECUTABLE} .
#sleep 5
#echo "runoncluster" > input.txt

sed -i 's/CLUSTER=0/CLUSTER=1/g' ${FILETORUN}
sleep 1
#stat ${FILETORUN}
#stat ${EXECUTABLE}

job=0
until [ $job -gt $NJOBS ]; do
  if ls ${OUTPUTDIR}/job${job}_* 1> /dev/null 2>&1; then
    echo "Output file already exists, now skipping!!!"
  else 
    sleep 2
    bsub -q psanaq -o "../logs/output"${job}".log" ./${FILETORUN} ${job}
  fi
  job=$((job + 1))
done


