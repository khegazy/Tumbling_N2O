#!/bin/bash 

#LOWQ=-FillQ zeros 
#LOWQ=-lowQtheory ./results/sim-phenoxyRadical_LowQfill_scale-0.400000_Bins[555].dat

./pairCorr.exe 20180627_1551  
./pairCorr.exe 20180629_1630 
#./pairCorr.exe 20180630_1925 
#./pairCorr.exe 20180701_0746 

#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_

#SMOOTH='0.025000'
#./pairCorr.exe 20180627_1551 -Idir /reg/ued/ana/scratch/nitroBenzene/mergeScans/ -Fname data-20180627_1551-tSmeared-${SMOOTH}-sMsAzmAvgDiff[1152,555].dat -Osuf -tSmeared -FillQ zeros
#./pairCorr.exe 20180629_1630 -Idir /reg/ued/ana/scratch/nitroBenzene/mergeScans/ -Fname data-20180629_1630-tSmeared-${SMOOTH}-sMsAzmAvgDiff[402,555].dat -Osuf -tSmeared -lowQtheory ./results/sim-phenoxyRadicalLowQfill[555].dat
#./pairCorr.exe 20180630_1925 -Idir /reg/ued/ana/scratch/nitroBenzene/mergeScans/ -Fname data-20180630_1925-tSmeared-${SMOOTH}-sMsAzmAvgDiff[402,555].dat -Osuf -tSmeared -lowQtheory ./results/sim-phenoxyRadicalLowQfill[555].dat
#./pairCorr.exe 20180701_0746 -Idir /reg/ued/ana/scratch/nitroBenzene/mergeScans/ -Fname data-20180701_0746-tSmeared-${SMOOTH}-sMsAzmAvgDiff[402,555].dat -Osuf -tSmeared -lowQtheory ./results/sim-phenoxyRadicalLowQfill[555].dat

#######################################
#####  Time Dependend Simulation  #####
#######################################

NTIME=450
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname dissociation_phenyl-N2O_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[400,555].dat -Osuf -dissociation_phenyl-N2O -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname dissociation_phenyl-N2O_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[400,555].dat -Osuf -dissociation_phenyl-N2O_timeSmoothed -FillQ zeros
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname dissociation_nitrosobenzene-O_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[400,555].dat -Osuf -dissociation_phenyl-N2O -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname dissociation_nitrosobenzene-O_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[400,555].dat -Osuf -dissociation_nitrosobenzene-O_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname dissociation_nitrosobenzene-O_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[400,555].dat -Osuf -dissociation_nitrosobenzene-O -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname rotate90_nitrobenzene_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[80,555].dat -Osuf -rotate90_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname rotate90_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[80,555].dat -Osuf -rotate90 -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname NOflop_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -NOflop -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname NOflop_nitrobenzene_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -NOflop_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname diagBend_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -diagBend -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname diagBend_nitrobenzene_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -diagBend_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname axisBend_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -axisBend -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname axisBend_nitrobenzene_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -axisBend_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname C0flop_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -C0flop -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/neh/home/khegazy/analysis/nitrobenzene/simulation/diffractionPattern/output/timeDependent -Fname C0flop_nitrobenzene_azmAvgSMS_timeSmoothed_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[200,555].dat -Osuf -C0flop_timeSmoothed -FillQ zeros 
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/timeDependent/ -Fname rotation_nitrobenzene_azmAvgSMS_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[800,555].dat -Osuf -rotation-nitrobenzene -FillQ false 

#######################################################
#####  Difference PC with Simulated Final States  #####
#######################################################

#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -Osuf PhenoxyDiff -saveLowQtheory true -SubR true -FillQ false -FitPC true

##############################################
#####  Reference Image / Theory Compare  #####
##############################################

#RUNNAMES=( "20180627_1551" "20180629_1630" "20180630_1925" "20180701_0746" )

#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname nitrobenzene_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -FillQ false -Osuf Nitrobenzene -FitPC true

#for run in "${RUNNAMES[@]}"
#do
#  ./pairCorr.exe ${run} -Idir ../staticDiffraction/results/ -Fname staticDiffraction_${run}[555].dat -Osuf -reference -FillQ theory -FitPC true
#done

### change molecule to phenoxyRadical
#./pairCorr.exe simulateReference -Idir /reg/ued/ana/scratch/nitroBenzene/simulations/ -Fname phenoxyRadical_sMsPatternLineOut_Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[555].dat -FillQ false
