#include "preProcessing.h"

using namespace cv;
using namespace std;



class centerfnctr {

    public:
        vector< vector<double> >* img;
        int count;
        //centerfnctr(vector< vector<double> > &img_i) {img = &img_i;}
        //centerfnctr(); 
        double operator() (vector<double> vect) {
          count++;
          //cout<<"calling fxn          ";
          return imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 2, 80, 20);
          //return fabs(1-imgProc::centerSymXsqr(img, vect[0], vect[1], 8, 3, 100, 20));
        }  
};
 


int main(int argc, char* argv[]) {

  if (argc != 2) {
    cerr<<"ERROR: Must run program with a list of txt filenames containing information about the images!!!"<<endl;
    exit(0);
  }

  ifstream fileNames(argv[1]);
  if (!fileNames.is_open()) {
    cerr<<"ERROR: Cannot open file "<<argv[1]<<endl;
    exit(0);
  }


  PLOTclass plt;
  int roiw = 821;
  int roih = 821;
  int hotPixel = 50000;
  int maxBkg=8;
  double bAvg=0;
  vector< vector<double> > imgBkg;
  vector< vector<double> > imgBkgOrig;
  vector< vector<double> > imgOrig;
  vector< vector<double> > imgOrigSave;
  vector< vector<double> > imgSave(roih), imgSubBkgSave(roih);
  float imgNorm, bkgNorm;

  centerfnctr centfnctr;
  //centerfnctr centfnctr(imgOrigSave);
  vector<double> center(2);
  int centerC=0;
  int centerR=0;
  int pcentR, pcentC, Nslices, NmatBins, NhistBins;
  int Nrows;
  int Ncols;
  int imgDepth;

  for (int ir=0; ir<roih; ir++) {
    imgSave[ir].resize(roiw); 	
    imgSubBkgSave[ir].resize(roiw);
  }

  vector< vector<double> > radBins, radBkgBins, radSBkgBins;
  vector< vector<double> > radSubBkgBins, radSubBkgA;
  vector< vector<double> > legendreCoeff;

  vector< pair<int,int> > pmask, nmask;
  //makeMask( PI/3, 0.1, 0.25, roih, roiw, pmask, nmask);
  vector<double> aniParam;

  int curScan, curLoop;
  ifstream txtInput;
  string fileName;
  string line;
  size_t ipos, iipos, fpos;
  vector<string> inpVals(22);
  vector< vector<string> > inpValues;
  string date, type, folder, imageName, rFileName;

  int scan, loop, imageNum, frame, imgTime;
  float stagePos, temperature, gain;

  TFile* file=NULL;
  TTree* tree=NULL;


  /////////////////////////////////////
  /////  Making background image  /////
  /////////////////////////////////////

  double Nbkgs = 1;
  ifstream bkgFiles;
  string bkgFile;
  bkgFiles.open("backgroundFiles.txt");
  if (!bkgFiles.is_open()) {
    cerr << "ERROR: Cannot open file backgroundFiles.txt!!!\n\n";
    exit(0);
  }

  getline(bkgFiles, bkgFile);
  Mat imgMatBkg = imread(bkgFile.c_str() ,CV_LOAD_IMAGE_ANYDEPTH);
  imgProc::threshold(imgMatBkg, hotPixel);

  Nrows = imgMatBkg.rows;
  Ncols = imgMatBkg.cols;

  imgBkgOrig.resize(Nrows);
  imgOrig.resize(Nrows);
  for (int ir=0; ir<Nrows; ir++) {
    imgBkgOrig[ir].resize(Ncols,0);
    imgOrig[ir].resize(Ncols,0);
  }

  while (getline(bkgFiles, bkgFile)) {
    Mat imgMatBkg1 = imread(bkgFile.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
    imgProc::threshold(imgMatBkg1, hotPixel);

    int index=0;
    uint16_t* pbkg = imgMatBkg1.ptr<uint16_t>(0);
    for (int ir=Nrows-1; ir>=0; ir--) {
      for (int ic=0; ic<Ncols; ic++) {
        index = ((Nrows-1)-ir)*Ncols + ic;
        imgBkgOrig[ir][ic] += (double)pbkg[index];
      }
    }
    Nbkgs++;
  }

  // Rotate Background by PI/2 to change polarization
  //imgBkg = imgProc::imgRotatePId2(imgBkgOrig);
  imgBkg = imgBkgOrig;


  // Finding image center 
  center[0] = 500;
  center[1] = 500;
  centfnctr.count = 0;
  centfnctr.img = &imgBkg;

  tools::powellMin<centerfnctr> (centfnctr, center, 20, 1, 1.5, 0.1);


  // Averaging background and calculating norm
  bkgNorm = 0;
  for (int ir=0; ir<Nrows; ir++) {
    for (int ic=0; ic<Ncols; ic++) {
      imgBkg[ir][ic] /= Nbkgs;

      if (sqrt(pow(ir-center[0],2) + pow(ic-center[1],2)) >= 250) {
        bkgNorm += imgBkg[ir][ic];
      }
    }
  }



  while (getline(fileNames, fileName)) {

    cout<<"Now looking at file "<<fileName<<endl;

    ipos=0; fpos=fileName.find("/",0);
    while (fpos != string::npos) {
      iipos=ipos; ipos=fpos; 
      fpos=fileName.find("/",ipos+1);
    }
    fpos = fileName.find(".",ipos);
    date = fileName.substr(iipos-12, 8);
    type = fileName.substr(ipos+1, fpos-ipos-1);
    folder = fileName.substr(0, ipos+1);

    txtInput.open(fileName);
    if (!txtInput.is_open()) {
      cerr<<"ERROR: Cannot open image info file "<<fileName<<endl;
      exit(0);
    }
    curScan=-1; curLoop=-1;
    file=NULL; tree=NULL;
    getline(txtInput, line); 	// Need to skip first line (labels)
    while (getline(txtInput, line)) {

      // Reading input data from input txt
      ipos=0; fpos=line.find(" ",0);
      for (int inp=0; inp<22; inp++) {
        inpVals[inp] = line.substr(ipos, fpos-ipos);
      	ipos = fpos+1;
        fpos = line.find(" ",fpos+1);
      }

      inpValues.push_back(inpVals);
    }

    for (uint k=0; k<inpValues.size(); k++) {
     

      scan = atoi(inpValues[k][0].c_str());
      loop = atoi(inpValues[k][1].c_str());
      imageNum = atoi(inpValues[k][2].c_str());
      frame = atoi(inpValues[k][3].c_str());
      stagePos = atof(inpValues[k][4].c_str());
      imgTime = atoi(inpValues[k][5].c_str());
      temperature = atof(inpValues[k][9].c_str());
      gain = atof(inpValues[k][14].c_str());

      /////// Retrieving image ///////
      imageName = folder+"scan"+std::string(inpValues[k][0])+"_L"+std::string(inpValues[k][1])+"/image"+std::string(inpValues[k][2])+"_frame"+std::string(inpValues[k][3])+".png";
      cout<<"    Now looking at image "<<imageName<<endl;
      Mat imgOrigMat = imread(imageName.c_str() ,CV_LOAD_IMAGE_ANYDEPTH);

      Nrows = imgOrigMat.rows;
      Ncols = imgOrigMat.cols;
      imgDepth = imgOrigMat.depth();
 
      if (imgOrigMat.channels() != 1) {cerr<<"ERROR: Cannot deal with multichanneled images!!!"<<endl; exit(0);}

      /////// Removing hotpixels from original ///////
      imgProc::threshold(imgOrigMat, hotPixel);


// No long neccesary 
/*
      /////// Resizing Vectors ///////
      if (curScan == -1) {
        imgOrigSave.resize(Nrows);   
        imgBkg4.resize(Nrows);   	//imgBkg6.resize(Nrows);   
        for (int ir=0; ir<Nrows; ir++) {
          imgOrigSave[ir].resize(Ncols);   
          imgBkg4[ir].resize(Ncols);   	//imgBkg6[ir].resize(Ncols);   
        }
      }
*/

 
      /////// Create new files for each scan and loop ///////
      if (scan != curScan || loop != curLoop) {
     	curScan = scan;
	curLoop = loop;
	if (file) {
	  tree->Write();
	  file->Close();
 	}

	rFileName = "/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/rootFiles/"+date+"_Scan"+std::string(inpValues[k][0])+"_Loop"+std::string(inpValues[k][1])+".root";
        file = TFile::Open(rFileName.c_str(), "RECREATE");
        tree = new TTree("physics","physics");

	tree->Branch("scanNum", &scan, "scan/I");
	tree->Branch("loopNum", &loop, "loop/I");
	tree->Branch("imageNum", &imageNum, "imageNum/I");
	tree->Branch("frameNum", &frame, "frameNum/I");
	tree->Branch("stagePos", &stagePos, "stagePos/F");
	tree->Branch("imgTime", &imgTime, "imgTime/I");
	tree->Branch("temperature", &temperature, "temperature/F");
	tree->Branch("gain", &gain, "gain/F");
        tree->Branch("centerC", &centerC, "centerC/I");
        tree->Branch("centerR", &centerR, "centerR/I");
        //tree->Branch("Nrows", &Nrows, "Nrows/I");
        //tree->Branch("Ncols", &Ncols, "Ncols/I");
        tree->Branch("imgNorm", &imgNorm, "imgNorm/F");
        //tree->Branch("bkg4Norm", &bkg4Norm, "bkg4Norm/F");
        //tree->Branch("bkg6Norm", &bkg6Norm, "bkg6Norm/F");
        tree->Branch("imgDepth", &imgDepth, "imgDepth/I");
        tree->Branch("imgOrig", &imgOrigSave);
        tree->Branch("img", &imgSave);
        tree->Branch("imgBkg", &imgBkg);
        //tree->Branch("imgB4", &imgBkg4);
        //tree->Branch("normBkg6", &imgBkg6);
        tree->Branch("imgSBkg", &imgSubBkgSave);
        //tree->Branch("imgSB4", &imgSubBkg4Save);
        //tree->Branch("imgSubBkg6", &imgSubBkg6Save);
        tree->Branch("pimg", &radBins);
        tree->Branch("pimgSB", &radSBkgBins);
        //tree->Branch("pimgNB4", &radBkg4Bins); //check if actually normalized
        //tree->Branch("pimgNB6", &radBkg6Bins); //check if actually normalized
        //tree->Branch("pimgBA", &radBkgA);
        //tree->Branch("pimgB4A", &radBkg4A);
        //tree->Branch("pimgSBA", &radSubBkgA);
        //tree->Branch("pimgSB4", &radSubBkg4Bins);
        //tree->Branch("pimgSB4A", &radSubBkg4A);
        //tree->Branch("pimgSb6", &radSubBkg6Bins);
	//tree->Branch("aniParam", &aniParam);
        //tree->Branch("radLegendreCoeff", &legendreCoeff);

/*
        //////////////////////////////////
        //  Creating Background Images  //
        //////////////////////////////////

	string bkg1Name = folder+"scan"+std::string(inpValues[k][0])+"_L"+std::string(inpValues[k][1])+"/image1_frame0.png";
	Mat imgBkg1 = imread(bkg1Name.c_str() ,CV_LOAD_IMAGE_ANYDEPTH); 
        imgProc::threshold(imgBkg1, hotPixel);

        int index=0;
	uint16_t* pbkg = imgBkg1.ptr<uint16_t>(0);
	for (int ir=Nrows-1; ir>=0; ir--) {
          for (int ic=0; ic<Ncols; ic++) {
            index = ((Nrows-1)-ir)*Ncols + ic;
            imgBkg4[ir][ic] = (double)pbkg[index]/4.0;
            //imgBkg6[ir][ic] = (double)pbkg[index]/6.0;
	  }
   	}

	Mat imgTbkg;
        for (uint ibkg=0; ibkg<inpValues.size(); ibkg++) {
  	  if ((atoi(inpValues[ibkg][2].c_str()) > maxBkg) || (atoi(inpValues[ibkg][2].c_str()) == 1)) continue;
	  imgTbkg = imread((folder+"scan"+std::string(inpValues[ibkg][0])+"_L"+std::string(inpValues[ibkg][1])+"/image"+std::string(inpValues[ibkg][2])+"_frame"+std::string(inpValues[ibkg][3])+".png").c_str(), CV_LOAD_IMAGE_ANYDEPTH); 
	  imgProc::threshold(imgTbkg, hotPixel);
	  pbkg = imgTbkg.ptr<uint16_t>(0);
  	  if (atoi(inpValues[ibkg][2].c_str()) <=4) {
	    for (int ir=Nrows-1; ir>=0; ir--) {
              for (int ic=0; ic<Ncols; ic++) {
           	index = ((Nrows-1)-ir)*Ncols + ic;
            	imgBkg4[ir][ic] += (double)pbkg[index]/4.0;
	      }
   	    }
          }
  	  if (atoi(inpValues[ibkg][2].c_str()) <=6) {
 	    for (int ir=Nrows-1; ir>=0; ir--) {
              for (int ic=0; ic<Ncols; ic++) {
           	index = ((Nrows-1)-ir)*Ncols + ic;
            	//imgBkg6[ir][ic] += (double)pbkg[index]/6.0;
	      }
   	    }
          }
	}
*/
      }


      //////////////////////////////////////////////////
      //  Filling the Variables and filling the tree  //
      //////////////////////////////////////////////////


      /////// Filling image vector //////
      int index=0;
      uint16_t* pv = imgOrigMat.ptr<uint16_t>(0);
      for (int ir=Nrows-1; ir>=0; ir--) {
        for (int ic=0; ic<Ncols; ic++) {
          index = ((Nrows-1)-ir)*Ncols + ic;
          imgOrig[ir][ic] = pv[index];
	}
      }           

      // Rotate image by PI/2
      //imgOrigSave = imgProc::imgRotatePId2(imgOrig);
      imgOrigSave = imgOrig;

      /////// Finding image center ///////
      center[0] = 500;
      center[1] = 500;
      centfnctr.count = 0;
      centfnctr.img = &imgOrigSave;

      tools::powellMin<centerfnctr> (centfnctr, center, 20, 1, 1.5, 0.1);

      centerR = center[0];
      centerC = center[1];


      Nrows = roih;
      Ncols = roiw;
      ///////  Calculating Image Norm  ///////
      imgNorm=0;
      //imgNorm=bkg4Norm=bkg6Norm=0;
      for (int ir=-roih/2+centerR; ir<roih/2+centerR; ir++) {
	for (int ic=-roiw/2+centerC; ic<roiw/2+centerC; ic++) {
	  if (sqrt(pow(ir-centerR,2) + pow(ic-centerC,2)) < 250) continue;
	  //bkg4Norm += imgBkg4[ir][ic];
     	  //bkg6Norm += imgBkg6[ir][ic];
	  imgNorm += imgOrigSave[ir][ic];
	}
      }

      /////// Filling image vector //////
      int icol, irow;
      for (int ir=0; ir<roih; ir++) {
	for (int ic=0; ic<roiw; ic++) {
	  irow = ir - roih/2 + centerR;
	  icol = ic - roiw/2 + centerC;
	  imgSave[ir][ic] = imgOrigSave[irow][icol];
	  imgSubBkgSave[ir][ic] = imgOrigSave[irow][icol] - (imgNorm/bkgNorm)*imgBkg[irow][icol];	
	  //imgSubBkg4Save[ir][ic] = imgOrigSave[irow][icol] - (imgNorm/bkg4Norm)*imgBkg4[irow][icol];	
	  //imgSubBkg6Save[ir][ic] = imgOrigSave[irow][icol] - (imgNorm/bkg6Norm)*imgBkg6[irow][icol];	
	}
      }
         

      //////  Polar Plots //////
      pcentR = roih/2;		pcentC = roiw/2;
      Nslices = 60;	NmatBins = 410;	   NhistBins = 70;
      radBins = imgProc::polarBinning(imgSave, pcentR, pcentC, Nslices, NmatBins, NhistBins);
      radSBkgBins = imgProc::polarBinning(imgSubBkgSave, pcentR, pcentC, Nslices, NmatBins, NhistBins);
      radBkgBins = imgProc::polarBinning(imgBkg, pcentR, pcentC, Nslices, NmatBins, NhistBins);
      //radBkg2Bins = imgProc::polarBinning(imgBkg2, center, 16, 400, 20);
      //radBkg4Bins = imgProc::polarBinning(imgBkg4, pcentR, pcentC, Nslices, NmatBins, NhistBins);
      //radBkg6Bins = imgProc::polarBinning(imgBkg6, center, 16, 400, 20);
      //radBkg8Bins = imgProc::polarBinning(imgBkg8, center, 16, 400, 20);
      //radSubBkg2Bins = imgProc::polarBinning(imgSubBkg2, center, 16, 400, 20);
      //radSubBkg4Bins = imgProc::polarBinning(imgSubBkg4Save, pcentR, pcentC, Nslices, NmatBins, NhistBins);
      //radSubBkg6Bins = imgProc::polarBinning(imgSubBkg6Save, roih/2, roiw/2, 80, 410, 10, -1);
      //radSubBkg8Bins = imgProc::polarBinning(imgSubBkg8, center, 16, 400, 20);


      ///////////////////////////// NOT VALID ///////////////////////////////
/*
      //////  Subtract Average of Bkg4  //////
      radSubBkg4A.resize(radBkg4Bins.size());
      radBkg4A.resize(radBkg4Bins.size());
      for (uint ir=0; ir<radBkg4Bins.size(); ir++) {
   	bAvg=0;
        radSubBkg4A[ir].resize(radBkg4Bins[ir].size());
        radBkg4A[ir].resize(radBkg4Bins[ir].size());
        for (uint ic=0; ic<radBkg4Bins[ir].size(); ic++) {
	  radBkg4Bins[ir][ic] *= (imgNorm/bkg4Norm);
	  bAvg += radBkg4Bins[ir][ic];
        }
	bAvg /= (double)radBkg4Bins[ir].size();
   	for (uint ic=0; ic<radBkg4Bins[ir].size(); ic++) {
	  radSubBkg4A[ir][ic] = radBins[ir][ic] - bAvg;
	  radBkg4A[ir][ic] = bAvg;
  	}
      }


      //////  Subtract Average Background of Each Image  //////
      radSubBkgA.resize(radBins.size());
      radBkgA.resize(radBins.size());
      for (uint ir=0; ir<radBins.size(); ir++) {
	radSubBkgA[ir].resize(radBins[ir].size());
	radBkgA[ir].resize(radBins[ir].size());
	for (uint ic=0; ic<radBins[ir].size(); ic++) radSubBkgA[ir][ic] = radBins[ir][ic];
      }

      // Calculate average per radius
      for (uint ir=0; ir<radSubBkgA.size(); ir++) {
	bAvg = 0;
	for (uint ic=0; ic<radSubBkgA[ir].size(); ic++) bAvg += radSubBkgA[ir][ic];
   	bAvg /= (double) radSubBkgA[ir].size();
	for (uint ic=0; ic<radSubBkgA[ir].size(); ic++) {
	  radSubBkgA[ir][ic] -= bAvg;
	  radBkgA[ir][ic] = bAvg;
	}
      }	
*/


      //////  Legendre Coefficients  //////
      /*for (uint irh=0; irh<radHists.size(); irh++) {
        legendreCoeff.push_back(imgProc::legendre1dFit(radHists[irh], 10));
      }
 
      string lkj = "histLegCoeff"+inpVals[2];
      TH2F* histLegCoeff = new TH2F(lkj.c_str(),lkj.c_str(),9,2,11,radHists.size()-5,0,0.5);  
      for (uint i=0; i<legendreCoeff.size()-5; i++) {
        for (uint k=2; k<legendreCoeff[i].size(); k++) {
            histLegCoeff->SetBinContent(k-1, i+1, legendreCoeff[i][k]);
        }
      }
      histLegCoeff->GetYaxis()->SetTitle("Radius");
      histLegCoeff->GetXaxis()->SetTitle("Legendre Order");
      histLegCoeff->SetMaximum(300); histLegCoeff->SetMinimum(-700);
      plt.print2d(histLegCoeff);
      */

      tree->Fill();
      legendreCoeff.clear();
      aniParam.clear();
      //delete histLegCoeff;
    }
    txtInput.close();
    tree->Write();
    file->Close();
    cout<<endl<<endl<<endl;
  }
  fileNames.close();

return 1;
}

                                        
