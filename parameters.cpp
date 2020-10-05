#include "/reg/neh/home5/khegazy/baseTools/tools/parameters.h"

parameterClass::parameterClass(std::string runName) {

  run = runName;

  // Molecule
  molecule = initialState;
  radicalNames.push_back("N2O");
  molName = radicalNames[molecule];


  // Image parameters
  QperPix = 0.028156*(137.0/112.0);
  hasI0 = false;
  Nlegendres = 1;
  NradLegBins = 50;
  NmaxRadBins = 750;
  NradAzmBins = 375;
  imgSize = 805;
  imgEdgeBuffer = 20;
  hasRef = true;
  refStagePosCut = 1054000;
  imgShutterTime = 60;
  imgNormRadMin = 0.16; //0.25; //0.4; //0.045; //0.06; //0.045;
  imgNormRadMax = 0.6; //0.721; //0.33; //0.721;


  // PV
  getPVs        = false;
  pvSampleTimes = 5;
  pressureSmear = 180;
  pvFolder      = "";

  // Power Scans
  /*
  range1Qbegin = 1.8;
  range1Qend   = 2.5;
  range2Qbegin = 2.8;
  range2Qend   = 3.75;
  */

  // Time Zero
  tZeroQranges.resize(0);
  tZeroRatio.resize(0);

  // Background removal
  XrayHighCut       = 5e4;
  XrayLowCut        = 4e4;
  XraySTDcut        = 3;
  XrayWindow        = 10;
  xRayHitDist       = false;

  refCorrection     = "NULL";

  hotPixel          = NAN;

  NshellOutlierLoops  = 3;
  bkgSTDcut           = 15;
  shellWidth          = 1;
  Npoly               = 3;
  stdCutLeft          = 3;
  stdChangeRatioLeft  = 0.39;
  stdAccRatioLeft     = 0.04;//0.070;
  fracShellSTDcutLeft = 0.0012;
  stdCutRight         = 4; //2.25; //2.75;
  stdChangeRatioRight = 0.17;
  stdAccRatioRight    = 0.01;//0.02;
  fracShellSTDcutRight= 0.002;
  outlierVerbose      = false;
  plotRadPixDist      = false;
  indicesPath         = "/reg/neh/home/khegazy/analysis/radialFitIndices/";

  /* For removeOutliers_stdInclude
  stdIncludeLeft      = 3; //1;
  distSTDratioLeft    = 0.5;
  stdIncludeRight     = 1;
  distSTDratioRight   = 0.75;
  */

  outlierSimpleSTDcut = -1; // for <0 do not use removeOutliersSimple

  outlierMapSTDcut        = 1.5;//75;
  outlierCoreValThresh    = 5000; //90; //65; //5e5;
  outlierCoreRad          = 2; 
  outlierClusterRad       = 3; 
  outlierMinClusterSize   = 50;
  outlierMinPixelSize     = outlierMinClusterSize;
  outlierMinDensity       = 0.0;//0.2;
  outlierShapeVarCut      = 400.5;
  outlierShapeEdgeCut     = 2.75;
  outlierBorderValThresh  = 10; //75; //1e5;
  outlierBorderDistLimit  = 2;
  outlierBorderRad        = 3;
  outlierPadRad           = 7;
  outlierrMaxScale        = 1;
  outlierrMinScale        = 0;
  outliercMaxScale        = 4;
  outliercMinScale        = 0;
 

  readoutStart         = 0.7; //0.94; // Use ratio < 1. Converts to bins at the end
  readoutEnd           = 1; // Use ratio < 1. Converts to bins at the end
 
  // Centering
  scanAvgCenter     = true;
  I0centers         = false;
  computeCenters    = false;
  centerFxnType     = 3;
  centerShellWidth  = 15; 
  centerSTDcut      = 3;

  meanInds.push_back(140); 
  meanInds.push_back(155); 
  meanInds.push_back(170); 
  meanInds.push_back(185); 
  meanInds.push_back(200); 

  meanInds.push_back(250); 

  centerDir = "/reg/ued/ana/scratch/N2O/preProcessing/centers/";

  // Q Filtering
  order  = 6;
  WnLow  = 0.005;
  WnHigh = 0.035;
  filterType = "lowpass";
  pltFilterVerbose = false;

  // Remove low order polynomial noise
  NlowOrderPoly         = 6;
  lowPolySubtractStudy  = false;


  // Merging scans
  mergeNormalizeImgs  = false;
  Qnormalize          = true;
  mergeSTDscale       = 3; //2.6; FIX ME CHANGE compare to Thomas
  mergeImageSTDScale  = 2.3;
  legImageNoiseCut    = 12;
  azmImageNoiseCut    = 105;

  mergeGaussSmoothRef = false;
  mergeGSmoothSTD     = 5;

  testMergeNbootStrap = false;
  useBootstrapSEM     = true;
  computeBootstrapSEM = false;
  mergeNbootstrap     = 10000;

  timeWnHigh = 0.8;
  timeFiltOrder  = 5;
  timeFilterType = "lowpass";
 
  smearTimeBinWindow  = 100;
  timeSmearSTD        = 0.025;
  scanImgAzmSTDcut    = 4;
  scanImgAzmRefSTDcut = 4;

  saveMergeIntermediates = false;
  saveMergeInterFolder = "/reg/ued/ana/scratch/N2O/mergeScans/intermediates";

  
  // Analysis Parameters
  signalRranges.resize(0);


  

  // Pair correlation parameters
  NautCpadding      = 10000;
  holeRat           = 0.15;
  rMaxLegRat        = 0.75;
  rMaxAzmRat        = 0.07; //0.052;
  padDecayRat       = 0.5;

  pCorrGaussFilter  = true;
  pCorrButterFilter = false;
  pCorrQcut         = 10;
  filterVar         = std::pow(NradAzmBins/4, 2);
  pCorrWnHigh       = 0.6;
  pCorrFilterOrder  = 2;
  pCorrFilterType   = "lowpass";
  lowQfillSimScale  = 0.3;
  fillLowQfile      = "NULL";
  fillLowQtheory    = false;
  fillLowQzeros     = true;
  fillLowQsine      = false;
  fillLowQfitTheory = false;
  useFilledHole     = false;


  subtractReference = false;

  // Simulation parameters
  compareSims     = false;
  simPairCorr     = true;
  getBonds        = true;
  simPltVerbose   = false;
  NradSimBins     = NradAzmBins;
  Iebeam          = 5;
  elEnergy        = 3.7e6;
  screenDist      = 4;
  xyzDir          = "/reg/neh/home/khegazy/analysis/2015/timeBasis_N2O_NO2/simulation/XYZfiles/";
  simOutputDir    = "/reg/ued/ana/scratch/N2O/simulations/";
  fsFitOffset     = false;
  fsFilterSMS     = false; 
  fsFilterVar     = std::pow(NradAzmBins/4.5, 2); 
  fitQbegin       = 1.2;
  fitQend         = 8;
  fsQfitBegin     = 1;
  fsQfitEnd       = 8;
  fsRfitBegin     = 1.1;
  fsRfitEnd       = 5;

  simHotFinalState  = false;
  hotFSrefVar       = filterVar*1.1;
  hotFStdepVar      = filterVar*0.9;

  preProcOutputDir    = "/reg/ued/ana/scratch/N2O/preProcessing/";
  preProcI0OutputDir  = "NULL";
  mergeScansOutputDir = "/reg/ued/ana/scratch/N2O/mergeScans/";
  scanSearchOutputDir = "/reg/ued/ana/scratch/N2O/scanSearch/";
  radialPixelDist     = "/reg/ued/ana/scratch/N2O/radialPixelDist/";
  backgroundImage     = "NULL"; // Use NULL if there is not background
  backgroundFolder    = "/reg/ued/ana/scratch/N2O/background/";
  indexPath           = "/reg/neh/home/khegazy/analysis/radialFitIndices/";

  pltCent     = false;
  verbose     = false; 
  pltVerbose  = false;


  int scaleStagePos = 1e4;
  if (runName.compare("fullRevival") == 0) {

    // Image parameters
    legStdCut = 3.0;
    NbinsSkip = 50; //31;
    imgMatType = "uint16";

    // Measurement parameters
    timeZero = 0.3;
    hasRef = false;

    // Filtering
    suppressBins = 75;
    padMaxHeight = 4;

    // Bad scans

    // Background
    backgroundImage = "NULL";
    hasLaserBkg = false;
    laserClusterRemoval = false;


    // Remove hole
    holeR = 540;
    holeC = 500;
    holeRad = 53;

    /////  Center Finding Parameters  /////
    // Rough center finding
    sigma = 8;
    doCOM_center = false;
    centR_estimate = 535;
    centC_estimate = 500;
    minRad = 100;
    maxRad = 325;
    meanInd = 350;
    COMstdScale = 0.075;

    cntrScale = 10;
    cntrMinScale = 1;
    cntrPowellTol = 1;
    cntrFracTol1d = 0.01;

    // Bad Regions //
    int row, col;
    int rad = 58;
    int hRad = 50;
    nanMap.clear();
    nanMap.resize(1024);
    for (uint ir=0; ir<nanMap.size(); ir++) {
      nanMap[ir].resize(1024, 0);
      for (uint ic=500; ic<nanMap[ir].size(); ic++) {
        nanMap[ir][ic] = 1;
      }
    }
    for (int ir=holeR-holeRad; ir<=holeR+holeRad; ir++) {
      for (int ic=holeC-holeRad; ic<=holeC+holeRad; ic++) {
        row = ir - holeR;
        col = ic - holeC;
        if (holeRad > std::sqrt(row*row + col*col)) {
          nanMap[ir][ic] = 1;
        }
      }
    }


  }
  else if (runName.compare("simulateReference") == 0) {
    
    NbinsSkip = 50; //28;
    switch (molecule) {
      case initialState: {
        xyzFiles.push_back("N2O.xyz");
        break;
      }
      default: {
        std::cerr << "ERROR: do not recognize molecule enum!!!\n";
        exit(0);
      }
    }
  }
  else if (runName.compare("doRunLists") == 0) {
  }
  else if (runName.compare("other") == 0) {
  }
  else {
    std::cerr << "ERROR: Cannot set values for run " + runName + "!!!\n";
    exit(0);
  }

  readoutAzmBinStart = readoutStart*NradAzmBins; //400
  readoutAzmBinEnd   = readoutEnd*NradAzmBins; //450
  readoutLegBinStart = readoutStart*NradLegBins;
  readoutLegBinEnd   = readoutEnd*NradLegBins;


  refStagePosCut *= scaleStagePos;
  maxQazm   = QperPix*(NradAzmBins - 1);
  maxRbins  = rMaxAzmRat*((NradAzmBins + NautCpadding)/2 + 1);
  maxR      = maxRbins*(2*PI/(QperPix*(NradAzmBins - 1 + NautCpadding)));

  R1Bin = 1/(2*PI/(QperPix*(NradAzmBins - 1 + NautCpadding)));
  R8Bin = 8/(2*PI/(QperPix*(NradAzmBins - 1 + NautCpadding)));

}


