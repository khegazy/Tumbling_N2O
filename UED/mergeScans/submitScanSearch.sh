#!/bin/bash 

MAINDIR=/reg/neh/home/khegazy/analysis/nitroBenzene/UED
MAINOUTPUTDIR=/reg/ued/ana/scratch/nitroBenzene/scanSearch/
FILETORUN=mergeScans.exe
REFERENCE=/reg/ued/ana/scratch/nitroBenzene/mergeScans/data-20180627_1551-referenceAzm[555].dat

if [ -z "$1" ]; then
  echo "ERROR SUBMITTING JOBS!!!   Must give the run name to search, or type ALL to search all runs!"
  exit
fi

RUNNAME=${1}"*"
if [ "${1}" = "ALL" ]; then
  RUNNAME="*"
fi


if [ -z "$2" ]; then
  NSCANS=3
else
  NSCANS=${2}
fi

if [ -z "$3" ]; then
  RESUBMIT=false
else
  RESUBMIT=true
fi

OUTPUTDIRSUFFIX="size"${NSCANS}
OUTPUTDIR=${MAINOUTPUTDIR}"/"${OUTPUTDIRSUFFIX}

echo ${OUTPUTDIR}
if [ ! -e "$OUTPUTDIR" ]; then
  mkdir ${OUTPUTDIR}
  cd ${OUTPUTDIR}
  mkdir logs
  mkdir plots
  cd -
fi

if [ ! -e "./runLists/scanSearch" ]; then
  mkdir ./runLists/scanSearch
fi

for file in runLists/run-${RUNNAME}
do

  NLINES=$(wc -l ${file})
  NLINES=${NLINES% run*}
  echo ${NLINES}

  CUTBEGINNING=${file#*run-}
  OUTPUTFILENAME=${CUTBEGINNING%.*}

  CUTBEGINNING=${file#*Lists/}
  NEWFILEPREFIX="./runLists/scanSearch/"${CUTBEGINNING%.*}

  CURLINE=1
  LASTLINE=$((NSCANS))
  until [ $CURLINE -gt $NLINES ]; do
    echo ${CURLINE}
    JOBNAME="scanLines"${CURLINE}"-"${LASTLINE}
    NEWRUNLIST=${NEWFILEPREFIX}"_"${JOBNAME}".txt"

    cp ${file} ${NEWRUNLIST}
    sed -i "${CURLINE},${LASTLINE}!d" ${NEWRUNLIST}

    echo ${OUTPUTDIR}"/logs/"${OUTPUTFILENAME}"-"${JOBNAME}".log"
    #bsub -q psanaq -o${OUTPUTDIR}"/logs/"${OUTPUTFILENAME}"-"${JOBNAME}".log" ./${FILETORUN} ${NEWRUNLIST} -ScanSearch ${JOBNAME} -OdirSuffix ${OUTPUTDIRSUFFIX} -CompareRef ${REFERENCE}
    bsub -q psanaq -o${OUTPUTDIR}"/logs/"${OUTPUTFILENAME}"-"${JOBNAME}".log" ./${FILETORUN} ${NEWRUNLIST} -ScanSearch ${JOBNAME} -OdirSuffix ${OUTPUTDIRSUFFIX}

    sleep 1

    CURLINE=$((CURLINE + 1))
    LASTLINE=$((LASTLINE + 1))
    if [ ${LASTLINE} -gt ${NLINES} ]; then
      LASTLINE=${NLINES}
    fi
  done
  #strPos=$(($strPos + $subDLength + 11))


  #CUTBEGINNING=${file#*runList_}
  #OUTPUTFILENAME=${CUTBEGINNING%.*} #${1}-${2}
  #ROOTFILENAME=${ROOTOUTDIR}${OUTPUTFILENAME}".root"

  #if [[ "$RESUBMIT" = true ]] && [[ ! -f ${ROOTFILENAME} ]]; then
  #  continue
  #fi

  #echo "Submitting "${file}
  #bsub -q psanaq -o${OUTPUTDIR}${OUTPUTFILENAME}".log" ./${FILETORUN} ${file}

done
