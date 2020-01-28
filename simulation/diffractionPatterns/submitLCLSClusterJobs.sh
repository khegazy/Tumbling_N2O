#!/bin/bash -x

FILETORUN=runLCLSDiffraction.sh
EXECUTABLE=n2oDiffSim.exe
MAINDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/
COPYPREF=n2oDiffSim

if [ -z "$2" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give output directory and number of jobs"
  exit
fi

NJOBS=${2}
OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/output/${1}

make clean; make

if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
  mkdir ${OUTPUTDIR}/logs
  mkdir ${OUTPUTDIR}/run
else
  if [ ! -d "logs" ]; then
    mkdir ${OUTPUTDIR}/logs
    mkdir ${OUTPUTDIR}/run
  fi
fi

cd ${OUTPUTDIR}/run
cp ${MAINDIR}/${FILETORUN} .
cp ${MAINDIR}/${EXECUTABLE} .
cp ${MAINDIR}/${COPYPREF}* .
#sleep 5
#echo "runoncluster" > input.txt

sed -i 's/CLUSTER=0/CLUSTER=1/g' ${FILETORUN}
#sleep 1
#stat ${FILETORUN}
#stat ${EXECUTABLE}

job=0
until [ $job -gt $NJOBS ]; do
  if ls ${OUTPUTDIR}/job${job}_* 1> /dev/null 2>&1; then
    echo "Output file already exists, now skipping!!!"
  else
    sleep 5
    bsub -W 100:10 -q psanaq -o "../logs/output"${job}".log" ./${FILETORUN} ${job}
  fi
  job=$((job + 1))
done


