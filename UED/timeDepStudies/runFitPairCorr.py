import os
import sys
sys.path.append("/reg/neh/home/khegazy/baseTools/UEDanalysis/plots/scripts")
sys.path.append("/reg/neh/home/khegazy/baseTools/UEDanalysis/timeDepStudies")
from plotClass import plotCLASS
from fitPairCorr import fitPairCorr
import tensorflow as tf
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from enum import IntEnum

class molOpts(IntEnum):
  data          = 0
  dataStatic    = 1
  nitrobenzene  = 2
  phenoxyNO     = 3
  phenoxyNOdiff = 4

molNames = ["data", "dataStatic", "nitrobenzene", "phenoxyNO", "phenoxyNOdiff"]


if __name__ == "__main__":

  #####  Variables  ##### 
  #molecule = molOpts.nitrobenzene
  #molecule = molOpts.phenoxyNOdiff
  molecule = molOpts.data
  doFit   = True
  sineFit = True
  debug   = False

  parameters = {
      "maxTrain"      : 20000,
      "saveEvery"     : 20,
      "atoms"         : ["hydrogen", "carbon", "nitrogen", "oxygen"],
      "Natoms"        : [5, 6, 1, 2],
      "Rrange"        : [0,9],
      "Nbins"         : 555,
      "startBin"      : 50,
      "NfitFxns"      : 500,
      "doNormalEqn"   : False,
      "qPerPix"       : 0.0223,
      "elEnergy"      : 3.7e6,
      "sineTransFile" : None,
      "checkpointPath": "./results/fitting/",
      "smoothLoss"    : None,
      "nonNegLoss"    : False,
      "noiseLoss"     : None, #1,
      "timeDepLoss"   : False,
      "L1regularize"  : None,
      "L1NRregularize": None,
      "L2regularize"  : None}

  if parameters["smoothLoss"]:
    filterSize = {
      25  : 1.5,
      50  : 2.5,
      75  : 2,
      100 : 2,
      200 : 5,
      300 : 8,
      500 : 16,
      750 : 22}
  else:
    filterSize = {
        25  : 1,
        50  : 2.5,
        75  : 3,
        100 : 4,
        200 : 8,
        300 : 11,
        500 : 20,
        750 : 30}
  plc = plotCLASS()

  
  # Simulation/Data Paths
  nitrobenzene_bondFile =\
      "/reg/ued/ana/scratch/nitroBenzene/simulations/nitrobenzene_bonds_"
  phenoxyNO_bondFile =\
      "/reg/ued/ana/scratch/nitroBenzene/simulations/phenoxyRadical_bonds_"

  simDir = "/reg/ued/ana/scratch/nitroBenzene/simulations/"
  phenoxyNO_simPath       = "phenoxyRadical_sMsPatternLineOut_"\
                              + "Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_"\
                              + "elE-3700000.000000_Bins[555].dat"
  nitrobenzene_simPath    = "nitrobenzene_sMsPatternLineOut_"\
                              + "Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_"\
                              + "elE-3700000.000000_Bins[555].dat"
  nitrobenzeneATM_simPath = "nitrobenzene_atmDiffractionPatternLineOut_"\
                              + "Qmax-12.376500_Ieb-5.000000_scrnD-4.000000_"\
                              + "elE-3700000.000000_Bins[555].dat"


  #####  Get Data  #####
  if molecule is molOpts.data:
    mergeFolder = "/reg/ued/ana/scratch/nitroBenzene/mergeScans/"
    fxnToFit  = np.reshape(
                  np.fromfile(
                    mergeFolder
                      +"data-20180627_1551-sMsAzmAvgDiff[29,555].dat",
                    dtype=np.double),
                  (29,555))
    variance  = np.reshape(
                  np.fromfile(
                    mergeFolder
                      +"data-20180627_1551-sMsStandardDev[29,555].dat",
                    dtype=np.double),
                    (29,555))**2
    ###  Gaussian Damping  ###
    gDamp     = np.exp(-1*np.arange(555.)**2/(2*(555./4)**2))
    fxnToFit *= np.reshape(gDamp, (1,-1))

    plc.print1d(variance[:7,:], "testVariance", isFile=False)
    parameters["sineTransFile"] = None
    #parameters["smoothLoss"]      = 20
    #parameters["noiseLoss"]     = True
    #parameters["timeDepLoss"]   = True
    parameters["L1regularize"]    = 1e2
    #parameters["L1NRregularize"]  = 1e-1
    #parameters["L2regularize"]    = 5e-7
  elif molecule is molOpts.dataStatic:
    fxnToFit  = np.fromfile(
                  "../staticDiffraction/results/fitting/"
                    +"staticDiffraction_20180627_1551[555].dat",
                  #"../mergeScans/results/fitting/"
                  #+"referenceAzm-20180627_1551[555].dat",
                  dtype=np.double)
    variance  = np.fromfile(
                  "../mergeScans/results/fitting/"
                  +"referenceAzmsMsStandardDev-20180627_1551[555].dat",
                  dtype=np.double)**2
    parameters["sineTransFile"] =\
                  "./results/fitting/fitPairCorrNitrobenzene_maxR-9.862436[369].dat"
    parameters["smoothLoss"]      = 20
    parameters["nonNegLoss"]      = True
    parameters["L1regularize"]    = 1e-4
    parameters["L1NRregularize"]  = 1e-1
    parameters["L2regularize"]    = 5e-6
    #parameters["noiseLoss"]     = True
  elif molecule is molOpts.nitrobenzene:
    fxnToFit    = np.fromfile(
                    simDir+nitrobenzene_simPath, 
                    dtype=np.double)
    variance    = np.ones(fxnToFit.shape[0])
    bondPrefix  = nitrobenzene_bondFile
    parameters["sineTransFile"] =\
                  "./results/fitting/fitPairCorrNitrobenzene_maxR-9.862436[369].dat"
    parameters["smoothLoss"]      = 20
    parameters["nonNegLoss"]    = True
    #parameters["noiseLoss"]     = True
    parameters["L1regularize"]    = 1e-4
    parameters["L1NRregularize"]  = 1e-1
    parameters["L2regularize"]    = 5e-6
  elif molecule is molOpts.phenoxyNO:
    fxnToFit    = np.fromfile(
                    simDir+phenoxyNO_simPath, 
                    dtype=np.double)
    variance    = np.ones(fxnToFit.shape[0])
    bondPrefix  = nitrobenzene_bondFile
    parameters["sineTransFile"] =\
                  "./results/fitting/fitPairCorr"
    parameters["smoothLoss"]      = 20
    parameters["nonNegLoss"]    = True
    #parameters["noiseLoss"]     = True
    parameters["L1regularize"]    = 1e-4
    parameters["L1NRregularize"]  = 1e-1
    parameters["L2regularize"]    = 5e-6
  elif molecule is molOpts.phenoxyNOdiff:
    fxnToFit  = np.fromfile(
                  simDir+phenoxyNO_simPath, 
                  dtype=np.double)\
                - np.fromfile(
                  simDir+nitrobenzene_simPath, 
                  dtype=np.double)
    variance  = np.ones(fxnToFit.shape[0])
    parameters["sineTransFile"] =\
                  "./results/fitting/fitPairCorrPhenoxyDiff_maxR-9.862436[369].dat"
    parameters["smoothLoss"]      = 20
    #parameters["noiseLoss"]     = True
    parameters["L1regularize"]    = 1e-6
    parameters["L1NRregularize"]  = 1e-1
    parameters["L2regularize"]    = 5e-7
  else:
    raise RuntimeError("Do not recognize molecule %s" % (repr(molecule))) 

  if len(fxnToFit.shape) == 1:
    fxnToFit = np.reshape(fxnToFit, (1, parameters["Nbins"]))
    variance = np.reshape(variance, (1, parameters["Nbins"]))

  atomicScat = np.fromfile(simDir+nitrobenzeneATM_simPath, dtype=np.double)


  # Get Bonds
  bondTypes    = []
  for i in range(len(parameters["atoms"])):
    if parameters["atoms"][i] == "hydrogen":
      continue
    for j in range(i, len(parameters["atoms"])):
      if parameters["atoms"][j] == "hydrogen":
        continue
      if (i == j) and (parameters["Natoms"][i] == 1):
        continue
      bondTypes.append(parameters["atoms"][i]+"-"+parameters["atoms"][j])

  bondCoeffs  = np.zeros(
                  (len(bondTypes), parameters["NfitFxns"]), 
                  dtype=np.float32)
  if molecule is molOpts.nitrobenzene or molecule is molOpts.phenoxyNO:
    for i,bt in enumerate(bondTypes):
      bonds = np.fromfile(bondPrefix+bt+".dat", np.double)
      for bnd in bonds:
        bInd = int(parameters["NfitFxns"]*(bnd-parameters["Rrange"][0])
                  /(parameters["Rrange"][1]-parameters["Rrange"][0]))
        bondCoeffs[i,bInd] += 2
  elif molecule is molOpts.phenoxyNOdiff:
    oldBondTypes  = bondTypes[:]
    bondTypes     = []
    btInd         = 0
    for bt in oldBondTypes:
      if ("hydrogen" in bt) or ("carbon-carbon" == bt):
        continue
      bondTypes.append(bt)
      bonds = np.fromfile(nitrobenzene_bondFile+bt+".dat", np.double)
      for bnd in bonds:
        bInd = int(parameters["NfitFxns"]*(bnd-parameters["Rrange"][0])
                  /(parameters["Rrange"][1]-parameters["Rrange"][0]))
        bondCoeffs[btInd,bInd] -= 2

      if ("carbon-nitrogen" != bt) and ("oxygen-oxygen" != bt):
        bonds = np.fromfile(phenoxyNO_bondFile+bt+".dat", np.double)
        for bnd in bonds:
          bInd = int(parameters["NfitFxns"]*(bnd-parameters["Rrange"][0])
                    /(parameters["Rrange"][1]-parameters["Rrange"][0]))
          bondCoeffs[btInd,bInd] += 2
      btInd += 1




  #########################
  #####  Build Model  #####
  #########################

  fitCls = fitPairCorr(
      parameters,
      fxnToFit,
      variance,
      sineFit,
      debug)


  #####################
  #####  Fitting  #####
  #####################
  
  with tf.Session() as sess:

    if not parameters["doNormalEqn"] and doFit:
      sess.run(tf.global_variables_initializer())
      sess.run(tf.local_variables_initializer())
      fitCls.fit(sess)
    elif not parameters["doNormalEqn"]:
      fitCls.saver.restore(sess, 
        parameters["checkpointPath"]+"fitting")




    ###  Plot Fit  ###
    fits = fitCls.get_fit(sess)
    if "startBin" in parameters:
      fits[:,:parameters["startBin"]] = 0
    opts = {
        "xTitle" : r"Q [$\AA^{-1}$]"}
    print("shapes",fxnToFit.shape,fits.shape)
    for i in range(fitCls.Nfits):
      plotFit   = np.stack((fxnToFit[i,:], fits[i,:]))
      plotName  = "./plots/"+molNames[molecule]+"-fit"
      if fitCls.Nfits > 1:
        plotName += "-time-" + str(i)
      plc.print1d(
          plotFit, 
          plotName,
          xRange=[0, parameters["Nbins"]*parameters["qPerPix"]],
          isFile=False,
          options=opts)

    ###  Plot Coefficients  ###
    sineTrans = None
    if parameters["sineTransFile"] is not None:
      sineTrans = fitCls.evaluate(sess, [fitCls.sineTransform])
      sineTrans[0] /= 7.
      
      plc.print1d(sineTrans[0],
        "./plots/"+molNames[molecule]+"-fitSineTrans",
        xRange=parameters["Rrange"],
        isFile=False)


    fitCoeffs = fitCls.get_fitCoeff(sess)

    """
    print("FTC ",fitCoeffs.shape)
    for i in range(fitCls.Nfits):
      plotFit   = np.stack((fxnToFit, fits[i,:,:]))
      plotName  = "./plots/"+molNames[molecule]+"-fitCoeff"
      if fitCls.Nfits > 1:
        plotName += "-time-" + str(i)
      for i in range(fitCoeffs.shape[0]):
        plc.print1d(fitCoeffs[i,:], 
            "./plots/fitCoeffs_" + bondTypes[i], 
            xRange=parameters["Rrange"],
            isFile=False)
    """

    #fitCoeffs = fitCoeffs.sum(axis=0)
    opts = {
        "xTitle" : r"R [$\AA$]"}
    gaussCoeffs = np.zeros_like(fitCoeffs[:,0,:])
    for i in range(fitCls.Nfits):
      plotName  = "./plots/"+molNames[molecule]+"-fitCoeff"
      if fitCls.Nfits > 1:
        plotName += "-time-" + str(i)

      currentCoeffs = fitCoeffs[i,:,:]
      if sineTrans is not None:
        opts["labels"] = ["sine Transform", "Fit"]
        currentCoeffs  = np.vstack((np.array(sineTrans)*0.5, currentCoeffs))
      else:
        opts["labels"] = ["Fit"]
      if molecule is not molOpts.data:
        opts["labels"] += bondTypes
        currentCoeffs  = np.vstack((currentCoeffs, bondCoeffs/175.))

      plc.print1d(
          currentCoeffs,
          plotName, 
          xRange=parameters["Rrange"],
          isFile=False,
          options=opts)

      plotName  = "./plots/"+molNames[molecule]+"-fitCoeffFilt"
      if fitCls.Nfits > 1:
        plotName += "-time-" + str(i)
      gaussCoeffs[i,:]  = gaussian_filter1d(
                            fitCoeffs[i,:,:], 
                            filterSize[parameters["NfitFxns"]], 
                            axis=1)
      currentCoeffs     = gaussCoeffs[i,:]*6
      if sineTrans is not None:
        opts["labels"] = ["sine Transform", "Fit"]
        currentCoeffs  = np.vstack((np.array(sineTrans), currentCoeffs))
      else:
        opts["labels"] = ["Fit"]
      if molecule is not molOpts.data:
        opts["labels"] += bondTypes
        currentCoeffs  = np.vstack((currentCoeffs, bondCoeffs/120.))

      plc.print1d(
          currentCoeffs,
          plotName, 
          xRange=parameters["Rrange"],
          isFile=False)

    if fitCls.Nfits > 1:
      opts = {
          "xTitle"  : "Time [ps]",
          "yTitle"  : r"Q [$\AA^{-1}$]"}
      timeDelay = np.fromfile("../mergeScans/results/timeDelays["
            + str(fitCls.Nfits + 1) + "].dat", np.double)

      plotName = "./plots/"+molNames[molecule]+"-fitTdCoeff"
      plc.print2d(
          fitCoeffs[:,0,:],
          plotName,
          X=timeDelay,
          yRange=parameters["Rrange"],
          isFile=False,
          options=opts)

      opts["xSlice"] = [-0.3, 2]
      plc.print2d(
          fitCoeffs[:,0,:],
          plotName+"Tslice",
          X=timeDelay,
          yRange=parameters["Rrange"],
          isFile=False,
          options=opts)
      del opts["xSlice"]

      plotName = "./plots/"+molNames[molecule]+"-fitTdCoeffFilt"
      plc.print2d(
          gaussCoeffs,
          plotName,
          X=timeDelay,
          yRange=parameters["Rrange"],
          isFile=False,
          options=opts)

      opts["xSlice"] = [-0.3, 2]
      plc.print2d(
          gaussCoeffs,
          plotName+"Tslice",
          X=timeDelay,
          yRange=parameters["Rrange"],
          isFile=False,
          options=opts)
      del opts["xSlice"]




    ########################
    #####  Debuggings  #####
    ########################

    if debug:
      output = [
          fitCls.deBrogW,
          fitCls.molAtomicScat,
          fitCls.scatAmps,
          fitCls.qInp,
          fitCls.qEval,
          fitCls.bondScatAmps,
          fitCls.bondNormAmps,
          fitCls.interp,
          fitCls.Sargs,
          fitCls.sinusiods,
          fitCls.sinusiodsDR]

      deBrog, atmNorm, scatAmps, qInp, qEval, bondScatAmps,\
      bondNormAmps, interp, Sargs, sinusiods, sinusiodsDR\
          = sess.run(output)

      print("atmScat", atmNorm)
      opts = {"yLog" : True}
      print("shape qs", qInp.shape, qEval.shape, bondScatAmps.shape)
      print("DeBroglie Wavelength (angs): ", deBrog) # forgot to m->angs
      plc.print1d(np.vstack((atomicScat,atmNorm)), "./plots/testFit_atomicScat", options=opts, isFile=False)
      opts["labels"] = fitCls.atoms
      plc.print1d(scatAmps, "./plots/testFit_scatteringAmps", options=opts, isFile=False)
      plc.print1d(qInp[0,:,0], "./plots/testFit_qInp", isFile=False)
      plc.print1d(qEval[0,:,0], "./plots/testFit_qEval", isFile=False)
      opts["labels"] = fitCls.bondTypes
      plc.print1d(bondScatAmps, "./plots/testFit_bondScatAmps", 
          options=opts, isFile=False)
      opts["labels"] = fitCls.bondTypes
      plc.print1d(bondNormAmps, "./plots/testFit_bondNormAmps", isFile=False)
      opts["labels"] = fitCls.atoms
      plc.print1d(interp, "./plots/testFit_interp", options=opts, isFile=False)
      plc.print2d(Sargs, "./plots/testFit_Sargs", isFile=False)
      plc.print2d(sinusiods, "./plots/testFit_sinusiods", isFile=False)
      plc.print2d(sinusiodsDR, "./plots/testFit_sinusiodsDR", isFile=False)



