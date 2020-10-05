#include "/reg/neh/home5/khegazy/baseTools/tools/parameters.h"

parameterClass::parameterClass(std::string runName) {

  // Molecule
  molecule = initialState;
  radicalNames.push_back("nitrousOxide");
  molName = radicalNames[molecule];


  // Image parameters
  Nlegendres      = 1;
  NradLegBins     = 50;
  NmaxRadBins     = 750;
  NradAzmBins     = 251;
  imgSize         = 895;
  imgEdgeBuffer   = 20;
  hasRef          = false;
  refStagePos     = -1;
  imgShutterTime  = 20;
  imgNormRadMin   = 0.045; 
  imgNormRadMax   = 0.721; 


  // PV
  getPVs        = false;
  pvSampleTimes = 5;
  pressureSmear = 180;
  pvFolder      = "";

  // Time Zero
  tZeroQranges.resize(3);
  tZeroQranges[0].push_back(3);     tZeroQranges[0].push_back(3.5);
  tZeroQranges[1].push_back(3.75);  tZeroQranges[1].push_back(5);
  tZeroQranges[2].push_back(2);     tZeroQranges[2].push_back(3);

  tZeroRatio.resize(2);
  tZeroRatio[0] = 1;
  tZeroRatio[1] = 0;

  // Merging scans
  normalizeImgs       = false;
  Qnormalize          = true;
  mergeSTDscale       = 3; 
  mergeImageSTDScale  = 2.3;
  legImageNoiseCut    = 12;
  azmImageNoiseCut    = 105;

  timeWnHigh          = 0.8;
  timeFiltOrder       = 5;
  timeFilterType      = "lowpass";
  smearTimeBinWindow  = 100;
  timeSmearSTD        = 0.025;
  scanImgAzmSTDcut    = 4;
  scanImgAzmRefSTDcut = 4;


  // Analysis Parameters
  signalRranges.resize(0);
  //signalRranges[0].push_back(1.1);  signalRranges[0].push_back(1.7);


  // Background removal
  XrayHighCut       = 3e4;
  XrayLowCut        = 1e4;
  XraySTDcut        = 3;
  XrayWindow        = 20;
  xRayHitDist       = false;

  hotPixel          = 1750;
  bkgSTDcut         = 15;
  shellWidth        = 1;
  Npoly             = 5;
  stdIncludeLeft    = 3; 
  distSTDratioLeft  = 0.5;
  stdCutLeft        = 3.5;
  meanBinSize       = 12;
  stdIncludeRight   = 1;
  distSTDratioRight = 0.75;
  stdChangeRatio    = 0.02;
  stdCutRight       = 7; 
  outlierSTDcut     = 7.;
  outlierVerbose    = false;
  radPixDist        = false;
  indicesPath       = "/reg/neh/home/khegazy/analysis/radialFitIndices/";

  outlierMapSTDcut        = 1.5;
  outlierCoreValThresh    = 5000;
  outlierCoreRad          = 2; 
  outlierClusterRad       = 3; 
  outlierMinClusterSize   = 50;
  outlierMinPixelSize     = outlierMinClusterSize;
  outlierMinDensity       = 0.0;
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
 

  readoutStart         = 0.7;   // Use ratio < 1. Converts to bins at the end
  readoutEnd           = 1;     // Use ratio < 1. Converts to bins at the end
   
  // Centering
  centerFxnType     = 3;
  centerShellWidth  = 15; 
  centerSTDcut      = 3;

  meanInds.push_back(140); 
  meanInds.push_back(155); 
  meanInds.push_back(170); 
  meanInds.push_back(185); 
  meanInds.push_back(200); 
  meanInds.push_back(250); 

  // Filtering
  order             = 5;
  WnLow             = 0.005;
  WnHigh            = 0.08;
  filterType        = "lowpass";
  pltFilterVerbose  = false;

  // Remove low order polynomial noise
  NlowOrderPoly         = 6;
  lowPolySubtractStudy  = false;

  // Pair correlation parameters
  NautCpadding      = 10000;
  holeRat           = 0.15;
  rMaxLegRat        = 0.75;
  rMaxAzmRat        = 0.07;
  padDecayRat       = 0.5;
  filterVar         = std::pow(NradAzmBins/4, 2);
  lowQfillSimScale  = 0.3;
  fillLowQfile      = "/reg/neh/home5/khegazy/analysis/pyridine/UED/timeDepStudies/results/sim-LowQfill[555].dat";
  fillLowQtheory    = true;


  subtractReference = true;

  // Simulation parameters
  compareSims     = false;
  simPairCorr     = true;
  getBonds        = true;
  simPltVerbose   = true;
  NradSimBins     = NradAzmBins;
  Iebeam          = 5;
  elEnergy        = 3.7e6;
  screenDist      = 4;
  xyzDir          = "../XYZfiles/";
  simOutputDir    = "./output/";
  fsFitOffset     = false;
  fsFilterSMS     = false; 
  fsFilterVar     = std::pow(NradAzmBins/4.5, 2); 
  fsQfitBegin     = 1;
  fsQfitEnd       = 4;
  fsRfitBegin     = 1.1;
  fsRfitEnd       = 8;
  //finalStates.push_back("");

  preProcOutputDir    = "/reg/ued/ana/scratch/pyridine/rootFiles/";
  mergeScansOutputDir = "/reg/ued/ana/scratch/pyridine/mergeScans/";
  scanSearchOutputDir = "/reg/ued/ana/scratch/pyridine/scanSearch/";
  radialPixelDist     = "/reg/ued/ana/scratch/pyridine/radialPixelDist/";
  backgroundImage     = "NULL";
  backgroundFolder    = "/reg/ued/ana/scratch/pyridine/background/";
  indexPath           = "/reg/neh/home/khegazy/analysis/radialFitIndices/";
  pltCent             = false;
  verbose             = true; 
  pltVerbose          = true;


  scaleStagePos = 1e4;
  if ((runName.compare("20180525_1914") == 0)
      || (runName.compare("e_off_laser_on") == 0)) {

    // Image parameters
    QperPix = 0.0223;
    legStdCut = 3.0;
    NbinsSkip = 30;

    // Filtering
    suppressBins = 75;
    padMaxHeight = 6.5;

    // Measurement parameters
    timeZero      = 0.125;
    hasRef        = true;
    refStagePos   = 154.3500;
    imgMatType    = "uint16";
    NfinalPoints  = 4;

    // Bad scans
    //badScans.push_back(scanNum);
    //badImages[scanNum].push_back(1530100);

    // Background
    backgroundImage = "backgroundImg-e_off_laser_on.dat";
    hasLaserBkg = true;
    laserClusterRemoval = false;

    int row, col;
    int rad = 52;
    int hRad = 55;
    nanMap.clear();
    nanMap.resize(1024);
    for (uint ir=0; ir<nanMap.size(); ir++) {
      nanMap[ir].resize(1024, 0);
    }
 
    // Remove hole
    holeR   = 573; 
    holeC   = 503; 
    holeRad = 45;

    /////  Center Finding Parameters  /////
    // Rough center finding
    sigma       = 8;
    blockCentR  = 573;
    blockCentC  = 500;
    minRad      = 70;
    maxRad      = 325;
    meanInd     = 350;
    COMstdScale = 0.075;

    cntrScale     = 10;
    cntrMinScale  = 1;
    cntrPowellTol = 1;
    cntrFracTol1d = 0.01;

    // PV files / variables
    //pvMap["pressure"]    = "pressureSampleChamber_06_27_2018_15_51_20-06_28_2018_13_27_15-15553.dat";
    //pvMap["UVcounts"]    = "UVsampleChamberCam_06_27_2018_15_51_20-06_28_2018_13_27_15-15553.dat";
    //pvMap["bunkerTemp"]  = "bunkerTemperature_06_27_2018_15_51_20-06_28_2018_13_27_15-15553.dat";
    //pvMap["highBayTemp"] = "highBayTemperature_06_27_2018_15_51_20-06_28_2018_13_27_15-15553.dat";

    throttle = 0; //uJ 

    // Variables for studies
    bkgStudyRanges.resize(0);
    //bkgStudyRanges[0].push_back(5.5);   bkgStudyRanges[0].push_back(6);

  }
  else if (runName.compare("simulateReference") == 0) {
    
    QperPix = 0.0223;
    NbinsSkip = 50; 
    switch (molecule) {
      case initialState: {
        xyzFiles.push_back("nitrousOxide.xyz");
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


  refStagePos *= scaleStagePos;
  maxQleg   = QperPix*NradLegBins;
  maxQazm   = 14.1207; //QperPix*NradAzmBins;
  maxRbins  = rMaxAzmRat*((NradAzmBins + NautCpadding)/2 + 1);
  maxR      = maxRbins*(2*PI/(QperPix*(NradAzmBins + NautCpadding)));
  std::cout<<"testing: "<<rMaxAzmRat<<"  "<<(NradAzmBins + NautCpadding)/2 + 1<<"  "<<QperPix<<"  "<<std::endl;

  R1Bin = 1/(2*PI/(QperPix*(NradAzmBins + NautCpadding)));
  R8Bin = 8/(2*PI/(QperPix*(NradAzmBins + NautCpadding)));
  std::cout<<"Rbins: "<<R1Bin<<"   "<<R8Bin<<"  "<<NradAzmBins + NautCpadding<<std::endl;

  std::cout<<"maxR: "<<maxR<<std::endl;
}


