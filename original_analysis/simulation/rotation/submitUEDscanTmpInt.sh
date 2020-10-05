#!/bin/bash

FILETORUN=simAlignment.exe
EXECUTABLE=simAlignment.exe
MAINDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation


STARTTEMP=30.0
STEPTEMP=1.0
ENDTEMP=150.0
NTEMPS=3
#121
#25

STARTINTS=0.5
STEPINTS=0.1
ENDINTS=10.0
NINTS=96

STARTTIME=315
ENDTIME=321
SAMPLESTEP=0.01

OUTPUTDIR=/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/UEDscanTmpInts/

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

TEMP=${STARTTEMP}
TCNT=0
echo ${TEMP}
until [ $TCNT -gt $NTEMPS ]; do
  ICNT=0
  INTS=${STARTINTS}
  until [ $ICNT -gt $NINTS ]; do
    bsub -q psanaq -o "../logs/output_"${TEMP}"_"${INTS}".log" ./${FILETORUN} -Exp "UED" -Temp ${TEMP} -Ints ${INTS} -StartTime ${STARTTIME} -EndTime ${ENDTIME} -SampleStep ${SAMPLESTEP} -PDForBases 0
    echo $TEMP"  "$INTS
    INTS=`echo $INTS + $STEPINTS |bc`
    ICNT=$((ICNT + 1))
  done
  TEMP=`echo $TEMP + $STEPTEMP |bc`
  TCNT=$((TCNT + 1))
done

