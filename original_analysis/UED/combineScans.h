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
#include <TF1.h>
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
#include "/reg/neh/home/khegazy/baseTools/tools/tools.h"
#include "/reg/neh/home/khegazy/baseTools/tools/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/plotClass.h"
#include "/reg/neh/home/khegazy/baseTools/tools/saveClass.h"



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
   Int_t           Nrows;
   Int_t           Ncols;
   Float_t         scanNorm;
   vector<double>          *imgNorms;
   vector<vector<double> > *img_flat;
   vector<vector<double> > *lgndr0;
   vector<vector<double> > *lgndr1;
   vector<vector<double> > *lgndr2;
   vector<vector<double> > *lgndr3;
   vector<vector<double> > *lgndr4;
   vector<vector<double> > *lgndr5;
   vector<vector<double> > *lgndr6;
   vector<vector<double> > *lgndr7;
   vector<vector<double> > *lgndr8;
   vector<vector<double> > *lgndr9;
   vector<vector<double> > *lgndr10;
   vector<vector<double> > *lgndr11;
   vector<vector<double> > *lgndr12;
   vector<vector<double> > *lgndr13;
   vector<vector<double> > *lgndr14;
   vector<vector<double> > *lgndr15;
   vector<vector<double> > *lgndr16;
   vector<vector<double> > *lgndr17;
   vector<vector<double> > *lgndr18;
   vector<vector<double> > *lgndr19;
   vector<vector<double> > *lgndr20;
   vector<vector<double> > *lgndr21;
   vector<vector<double> > *lgndr22;
   vector<vector<double> > *lgndr23;
   vector<vector<double> > *lgndr24;
   vector<vector<double> > *lgndr25;
   vector<vector<double> > *lgndr26;
   vector<vector<double> > *lgndr27;
   vector<vector<double> > *lgndr28;
   vector<vector<double> > *lgndr29;
   vector<vector<double> > *lgndr30;

   // List of branches
   TBranch        *b_scanNum;   //!
   TBranch        *b_loopNum;   //!
   TBranch        *b_scanNorm;   //!
   TBranch        *b_Nrows;   //!
   TBranch        *b_Ncols;   //!
   TBranch        *b_imgNorms;   //!
   TBranch        *b_img_flat;   //!
   TBranch        *b_lgndr0;   //!
   TBranch        *b_lgndr1;   //!
   TBranch        *b_lgndr2;   //!
   TBranch        *b_lgndr3;   //!
   TBranch        *b_lgndr4;   //!
   TBranch        *b_lgndr5;   //!
   TBranch        *b_lgndr6;   //!
   TBranch        *b_lgndr7;   //!
   TBranch        *b_lgndr8;   //!
   TBranch        *b_lgndr9;   //!
   TBranch        *b_lgndr10;   //!
   TBranch        *b_lgndr11;   //!
   TBranch        *b_lgndr12;   //!
   TBranch        *b_lgndr13;   //!
   TBranch        *b_lgndr14;   //!
   TBranch        *b_lgndr15;   //!
   TBranch        *b_lgndr16;   //!
   TBranch        *b_lgndr17;   //!
   TBranch        *b_lgndr18;   //!
   TBranch        *b_lgndr19;   //!
   TBranch        *b_lgndr20;   //!
   TBranch        *b_lgndr21;   //!
   TBranch        *b_lgndr22;   //!
   TBranch        *b_lgndr23;   //!
   TBranch        *b_lgndr24;   //!
   TBranch        *b_lgndr25;   //!
   TBranch        *b_lgndr26;   //!
   TBranch        *b_lgndr27;   //!
   TBranch        *b_lgndr28;   //!
   TBranch        *b_lgndr29;   //!
   TBranch        *b_lgndr30;   //!



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
   fChain->SetBranchAddress("Nrows", &Nrows, &b_Nrows);
   fChain->SetBranchAddress("Ncols", &Ncols, &b_Ncols);
   fChain->SetBranchAddress("img_flat", &img_flat, &b_img_flat);
   fChain->SetBranchAddress("imgNorms", &imgNorms, &b_imgNorms);
   fChain->SetBranchAddress("lgndr0", &lgndr0, &b_lgndr0);
   fChain->SetBranchAddress("lgndr1", &lgndr1, &b_lgndr1);
   fChain->SetBranchAddress("lgndr2", &lgndr2, &b_lgndr2);
   fChain->SetBranchAddress("lgndr3", &lgndr3, &b_lgndr3);
   fChain->SetBranchAddress("lgndr4", &lgndr4, &b_lgndr4);
   fChain->SetBranchAddress("lgndr5", &lgndr5, &b_lgndr5);
   fChain->SetBranchAddress("lgndr6", &lgndr6, &b_lgndr6);
   fChain->SetBranchAddress("lgndr7", &lgndr7, &b_lgndr7);
   fChain->SetBranchAddress("lgndr8", &lgndr8, &b_lgndr8);
   fChain->SetBranchAddress("lgndr9", &lgndr9, &b_lgndr9);
   fChain->SetBranchAddress("lgndr10", &lgndr10, &b_lgndr10);
   fChain->SetBranchAddress("lgndr11", &lgndr11, &b_lgndr11);
   fChain->SetBranchAddress("lgndr12", &lgndr12, &b_lgndr12);
   fChain->SetBranchAddress("lgndr13", &lgndr13, &b_lgndr13);
   fChain->SetBranchAddress("lgndr14", &lgndr14, &b_lgndr14);
   fChain->SetBranchAddress("lgndr15", &lgndr15, &b_lgndr15);
   fChain->SetBranchAddress("lgndr16", &lgndr16, &b_lgndr16);
   fChain->SetBranchAddress("lgndr17", &lgndr17, &b_lgndr17);
   fChain->SetBranchAddress("lgndr18", &lgndr18, &b_lgndr18);
   fChain->SetBranchAddress("lgndr19", &lgndr19, &b_lgndr19);
   fChain->SetBranchAddress("lgndr20", &lgndr20, &b_lgndr20);
   fChain->SetBranchAddress("lgndr21", &lgndr21, &b_lgndr21);
   fChain->SetBranchAddress("lgndr22", &lgndr22, &b_lgndr22);
   fChain->SetBranchAddress("lgndr23", &lgndr23, &b_lgndr23);
   fChain->SetBranchAddress("lgndr24", &lgndr24, &b_lgndr24);
   fChain->SetBranchAddress("lgndr25", &lgndr25, &b_lgndr25);
   fChain->SetBranchAddress("lgndr26", &lgndr26, &b_lgndr26);
   fChain->SetBranchAddress("lgndr27", &lgndr27, &b_lgndr27);
   fChain->SetBranchAddress("lgndr28", &lgndr28, &b_lgndr28);
   fChain->SetBranchAddress("lgndr29", &lgndr29, &b_lgndr29);
   fChain->SetBranchAddress("lgndr30", &lgndr30, &b_lgndr30);

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

  
