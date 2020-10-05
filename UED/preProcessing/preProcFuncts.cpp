#include "/reg/neh/home5/khegazy/baseTools/UEDanalysis/preProcessing/preProcessing.h"



/////////////////////////////
/////  Days in a Month  /////
/////////////////////////////

std::map<int, int> ppFunct::monthLengths = {
  {1 , 31},
  {2 , 28},
  {3 , 31},
  {4 , 30},
  {5 , 31},
  {6 , 30},
  {7 , 31},
  {8 , 31},
  {9 , 30},
  {10 , 31},
  {11 , 30},
  {12 , 31}};


////////////////////////////////
////  Retrieving File Info  ////
////////////////////////////////

bool ppFunct::getScanRunInfo(std::vector<imgProc::imgInfoStruct> &imgINFO, std::string runListName, bool verbose) {


  ifstream fileNames(runListName.c_str());
  if (!fileNames.is_open()) {
    cerr << "ERROR: Cannot open runList file " << runListName << endl;
    exit(0);
  }

  size_t ipos, fpos;
  string fileName;
  std::map<float, imgProc::imgInfoStruct> imgInfoMap;

  std::string runType = "Run";
  if (runListName.find("Run") == string::npos) {
    cerr << "ERROR: Cannot find the runType!!!" <<endl;
    exit(0);
  }

  ipos = runListName.find(runType);
  ipos = runListName.find("-", ipos);
  fpos = runListName.find("_Scan", ipos);
  std::string curRun = runListName.substr(ipos+1, fpos-ipos-1);
  if (verbose) cout << "\tRUN: " << curRun;

  //ipos = runListName.find("converted");
  //fpos = runListName.find("/N2O", ipos);
  //std::string curDate = runListName.substr(ipos+10, fpos-ipos-10);
  std::string curDate = "20150904";
  if (verbose) std::cout << "\tDATE: " << curDate;

  ipos = runListName.find("Scan");
  fpos = runListName.find(".txt", ipos);
  int curScan = stoi(runListName.substr(ipos+5, fpos-ipos-5));
  if (verbose) cout << "\tSCAN: " << curScan;

  std::map<int, float> stagePositions;
  std::vector<int> months;
  std::vector<int> days;
  std::vector<int> times;
  int time = 0;

  /////  Import Measurement Parameters  /////
  std::string paramFile = "/reg/ued/ana/scratch/N2O/data/converted/";
  paramFile += curDate + "/N2O/";
  if (curScan == 321) {
    paramFile += "background/background.txt";
  }
  else {
    paramFile += curRun + "/" + curRun + ".txt";
  }

  ifstream fileParams(paramFile.c_str());
  if (!fileParams.is_open()) {
    cerr << "\nERROR: Cannot open param file " << paramFile << endl;
    exit(0);
  }

  std::string params, stgP1;
  getline(fileParams, params);
  cout<<endl;
  while (getline(fileParams, params)) {
   
    ipos = params.find(" ");
    ipos = params.find(" ", ipos+1);
    fpos = params.find(" ", ipos+1);
    int imgNum = stoi(params.substr(ipos+1, fpos-ipos-1));


    ipos = params.find(" ", fpos+1);
    fpos = params.find(".");
    stgP1 = params.substr(ipos+1, fpos-ipos-1);
    ipos = fpos;
    fpos = params.find("E");
    int32_t stgPos = stoi(stgP1+params.substr(ipos+1, fpos-ipos-1));
    if (params.find("E+1") != string::npos) {
      stgPos /= 10;
    }
    stagePositions[imgNum] = stgPos;
  }

  while (getline(fileNames, fileName)) {

    if (verbose) std::cout << "Now looking at file " << fileName << endl;

    imgProc::imgInfoStruct imgInfo;


    imgInfo.run = curRun;
    imgInfo.scan = curScan;
    imgInfo.date = curDate;
    imgInfo.runType = runType;
    
    // Finding file name
    ipos = fileName.rfind("/");
    imgInfo.fileName = fileName.substr(ipos, fileName.length()-ipos);
    imgInfo.path = fileName.substr(0, ipos);
    if (verbose) cout << "\tFILEPATH: " << imgInfo.path;
    if (verbose) cout << "\tFILENAME: " << imgInfo.fileName;

    // Finding image number
    ipos = imgInfo.fileName.find("image");
    fpos = imgInfo.fileName.find("_frame", ipos);
    imgInfo.imgNum    = stoi(imgInfo.fileName.substr(ipos+5, fpos-ipos-5));
    imgInfo.throttle  = -1;
    if (verbose) cout << "\tIMGNUM: " << imgInfo.imgNum;

    // Finding stage position
    imgInfo.stagePos = stagePositions[imgInfo.imgNum];
    if (verbose) cout << "\tSTAGE POSITION: " << imgInfo.stagePos;

    if (verbose) cout << endl;
    imgInfoMap[imgInfo.stagePos] = imgInfo;
  }

  for (auto itr: imgInfoMap) {
    itr.second.time = -1;//times[itr.second.imgNum-1];
    imgINFO.push_back(itr.second);
  }

  fileNames.close();
  if (verbose) cout << "\nINFO: Finished retrieving info from runList!\n\n";

  // Check stagePos is ordered
  if (verbose) {
    cout << "Checking images are in order of stage position\n";
    for (uint i=0; i<imgINFO.size(); i++) {
      cout << imgINFO[i].scan << "  " << imgINFO[i].run 
        << "  " << imgINFO[i].stagePos << endl;
    }
  }


  return true;
}

bool ppFunct::getI0RunInfo(std::map<int, std::string> &I0fileNames,
    std::string runListName, bool verbose) {

  return false;
}


void ppFunct::makeRunLists(std::vector<imgProc::imgInfoStruct> &imgINFO,
    std::string runName, std::string preProcFolder) {
  exit(1);
}
