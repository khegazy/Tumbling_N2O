#!/bin/bash -x

FILETORUN=runDiffraction.sh
EXECUTABLE=n2oDiffSim.exe
MAINDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns
COPYPREF=n2oDiffSim
XYZFILE=N2O.xyz

if [ -z "$1" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give the number of jobs"
  exit
fi

NJOBS=${1}
OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/output/${1}
OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/clusterJobs

make clean; make

if [ ! -d "${OUTPUTDIR}" ]; then
  mkdir ${OUTPUTDIR}
  mkdir ${OUTPUTDIR}/logs
  mkdir ${OUTPUTDIR}/run
else
  if [ ! -d "${OUTPUTDIR}/logs" ]; then
    mkdir ${OUTPUTDIR}/logs
    mkdir ${OUTPUTDIR}/run
  fi
fi

cd ${OUTPUTDIR}/run
cp ${MAINDIR}/${FILETORUN} .
cp ${MAINDIR}/${EXECUTABLE} .
cp ${MAINDIR}/${COPYPREF}* .
cp ${MAINDIR}/${XYZFILE} .
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
    sleep 2
    #bsub -W 100:10 -q psanaq -o "../logs/output"${job}".log" ./${FILETORUN} ${job}
    bsub -q psanaq -W 6000 -o "../logs/output"${job}".log" ./${FILETORUN} ${job}
  fi
  job=$((job + 1))
done


