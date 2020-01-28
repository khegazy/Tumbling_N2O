#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>

//From ROOT
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

//From OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>

//Home Grown
#include "/reg/neh/home/khegazy/baseScripts/tools.h"
#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/saveClass.h"



using namespace std;

//////////////////////////////////////
//   Declaration of useful classes  //
//////////////////////////////////////

  PLOTclass plt;


//////////////////////////////////////////
//   Declaration of analysis variables  //
//////////////////////////////////////////

   TChain         *fChain;   //!pointer to the analyzed TTree or TChain 
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           scanNum;
   Int_t           loopNum;
   Float_t         scanNorm;
   Int_t           imageNum;
   Int_t           frameNum;
   Float_t         stagePos;
   Int_t           imgTime;
   Float_t         temperature;
   Float_t         gain;
   Int_t           centerC;
   Int_t           centerR;
   Float_t         imgNorm;
   Float_t         imgDepth;
   vector<vector<double> > *imgOrig;
   vector<vector<double> > *img;
   vector<vector<double> > *imgBkg;
   vector<vector<double> > *imgSBkg;
   vector<vector<double> > *pimg;
   vector<vector<double> > *pimgSB;

   // List of branches
   TBranch        *b_scanNum;   //!
   TBranch        *b_loopNum;   //!
   TBranch        *b_scanNorm;   //!
   TBranch        *b_imageNum;   //!
   TBranch        *b_frameNum;   //!
   TBranch        *b_stagePos;   //!
   TBranch        *b_imgTime;   //!
   TBranch        *b_temperature;   //!
   TBranch        *b_gain;   //!
   TBranch        *b_centerC;   //!
   TBranch        *b_centerR;   //!
   TBranch        *b_imgNorm;   //!
   TBranch        *b_imgDepth;   //!
   TBranch        *b_imgOrig;   //!
   TBranch        *b_img;   //!
   TBranch        *b_imgBkg;   //!
   TBranch        *b_imgSBkg;   //!
   TBranch        *b_pimg;   //!
   TBranch        *b_pimgSB;   //!



//////////////////////////////////////
//  analysisClass for loading data  // 
//////////////////////////////////////

class analysisClass {

   public:
        analysisClass(string fileList);
        analysisClass(string fileList, string treeName);
        void initialize(string fileList, string treeName);
        ~analysisClass();
        uint64_t setupEnvironment();
        int loadEvent(uint64_t entry);

	clock_t start, stop;
};      


analysisClass::analysisClass(string fileList) {

  string treeName = "physics";
  initialize(fileList, treeName);
}


analysisClass::analysisClass(string fileList, string treeName) {

  initialize(fileList, treeName);
}
  

void analysisClass::initialize(string fileList, string treeName) {

  start=clock();

  ifstream files;
  files.open(fileList.c_str());
  if (!files.is_open()) {
    cerr<<"ERROR: Cannot open file list "<<fileList<<"!!!"<<endl;
    exit(0);
  } 
  
  ///// Initialize fChain (Chain of trees)
  fChain = new TChain(treeName.c_str());
  
  ///// Add trees from files in fileList to fChain
  string line;
  TFile* ftest;
  cout<<"!!!!!  Adding trees from files to fChain  !!!!!"<<endl;
  while (getline(files, line)) {
    cout<<"Adding tree '"<<treeName<<"' from file '"<<line<<"'"<<endl;
    ftest = TFile::Open(line.c_str());
    if (!ftest->Get(treeName.c_str())) cerr<<"WARNING: Cannot find tree "<<treeName<<" in file "<<line<<" !!!"<<endl;
    fChain->Add(line.c_str()); 
  }
  cout<<endl<<endl;

  ///// Set variable addresses to values

   fChain->SetBranchAddress("scanNum", &scanNum, &b_scanNum);
   fChain->SetBranchAddress("loopNum", &loopNum, &b_loopNum);
   fChain->SetBranchAddress("scanNorm", &scanNorm, &b_scanNorm);
   fChain->SetBranchAddress("imageNum", &imageNum, &b_imageNum);
   fChain->SetBranchAddress("frameNum", &frameNum, &b_frameNum);
   fChain->SetBranchAddress("stagePos", &stagePos, &b_stagePos);
   fChain->SetBranchAddress("imgTime", &imgTime, &b_imgTime);
   fChain->SetBranchAddress("temperature", &temperature, &b_temperature);
   fChain->SetBranchAddress("gain", &gain, &b_gain);
   fChain->SetBranchAddress("centerC", &centerC, &b_centerC);
   fChain->SetBranchAddress("centerR", &centerR, &b_centerR);
   fChain->SetBranchAddress("imgNorm", &imgNorm, &b_imgNorm);
   fChain->SetBranchAddress("imgDepth", &imgDepth, &b_imgDepth);
   fChain->SetBranchAddress("imgOrig", &imgOrig, &b_imgOrig);
   fChain->SetBranchAddress("img", &img, &b_img);
   fChain->SetBranchAddress("imgBkg", &imgBkg, &b_imgBkg);
   fChain->SetBranchAddress("imgSBkg", &imgSBkg, &b_imgSBkg);
   fChain->SetBranchAddress("pimg", &pimg, &b_pimg);
   fChain->SetBranchAddress("pimgSB", &pimgSB, &b_pimgSB);

}


analysisClass::~analysisClass() {

  delete fChain;
  stop = clock();
  cout<<endl<<"Total running time: "<<double(stop-start)/CLOCKS_PER_SEC<<endl<<endl;
}


uint64_t analysisClass::setupEnvironment() {

  /////////////////////////////////////////////////////////////////////
  //  Here goes anything needed to be setup before running the code  //
  /////////////////////////////////////////////////////////////////////

  loadEvent(0); 

  return fChain->GetEntries();
}


int analysisClass::loadEvent(uint64_t entry) {

  return fChain->GetEntry(entry);
}

  
