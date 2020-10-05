#include "preProcessingMC.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }

  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName


  ///// Load environment and get the number of events /////
  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup


  FILE* refFile;
  refFile = fopen("/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/output/references/diffractionPattern_Bins-1024_Qmax-14.120700_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000.dat", "rb");
  double* bkgImg = new double[1024*1024];
  fread(bkgImg, sizeof(double), 1024*1024, refFile);
  fclose(refFile);
  
  refFile = fopen("/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/output/references/atmDiffractionPattern_Bins-1024_Qmax-14.120700_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000.dat", "rb");
  double* atmImg = new double[1024*1024];
  fread(atmImg, sizeof(double), 1024*1024, refFile);
  fclose(refFile);
 
  double imgRat = (*imgBkg)[0][0]/atmImg[0];
  vector< vector<double> > pimgSBproc, pimgproc;
  vector< vector<double> > imgSBproc;
  vector< vector<double> > imgOrigproc;
  vector< vector<double> > imgproc;
  vector< vector<double> > imgBkgproc;
  vector< vector<double> > imgSBprocOrig((*imgSBkg).size());
  vector< vector<double> > imgOrigprocOrig((*imgOrig).size());
  vector< vector<double> > imgprocOrig((*img).size());
  vector< vector<double> > imgBkgprocOrig((*imgSBkg).size());
  for (uint ir=0; ir<(*imgSBkg).size(); ir++) {
    imgSBprocOrig[ir].resize((*imgSBkg)[ir].size());
    imgprocOrig[ir].resize((*img)[ir].size());
    imgOrigprocOrig[ir].resize((*imgOrig)[ir].size());
    imgBkgprocOrig[ir].resize((*imgBkg)[ir].size());
  }

  int scanN = scanNum;
  int loopN = loopNum;
  float scanNrm = scanNorm;
  int imgN = imageNum + 1;
  int frmN = frameNum;
  float stgPos = stagePos;
  int imgTm = imgTime;
  float tmp = temperature;
  float gain = gain;
  int pcenterC = centerC;
  int pcenterR = centerR;
  float imgNrm = imgNorm;
  int imgDpth = imgDepth;


  TFile* outFile;
  size_t fpos = fileList.find("List");
  size_t ipos = fileList.find("Bend");
  if (ipos == string::npos) {
    ipos = fileList.find("bend");
  }
  if (ipos != string::npos) {   
    
    outFile = TFile::Open(("/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/rootFiles/procSimulation" 
		+ fileList.substr(ipos, fpos - ipos) + ".root").c_str(), "RECREATE");
  }
  else {
    outFile = TFile::Open("/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/rootFiles/procSimulation.root", "RECREATE");
  }
  TTree* tr = new TTree("physics", "physics");

   tr->Branch("scanNum", &scanN, "scanNum/I");
   tr->Branch("loopNum", &loopN, "loopNum/I");
   tr->Branch("scanNorm", &scanNrm, "scanNorm/F");
   tr->Branch("imageNum", &imgN, "imageNum/I");
   tr->Branch("frameNum", &frmN, "frameNum/I");
   tr->Branch("stagePos", &stgPos, "stagePos/F");
   tr->Branch("imgTime", &imgTm, "imgTime/I");
   tr->Branch("temperature", &tmp, "temperature/F");
   tr->Branch("gain", &gain, "gain/F");
   tr->Branch("centerC", &pcenterC, "centerC/I");
   tr->Branch("centerR", &pcenterR, "centerR/I");
   tr->Branch("imgNorm", &imgNrm, "imgNorm/F");
   tr->Branch("imgDepth", imgDpth, "imgDepth/F");
   tr->Branch("imgOrig", &imgOrigproc);
   tr->Branch("img", &imgproc);
   tr->Branch("imgBkg", &imgBkgproc);
   tr->Branch("imgSBkg", &imgSBproc);
   tr->Branch("pimg", &pimgproc);
   tr->Branch("pimgSB", &pimgSBproc);

  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

    if ((imageNum < 1) || (imageNum > 40)) continue;

    scanN = scanNum;
    loopN = loopNum;
    scanNrm = scanNorm;
    imgN = imageNum + 1 - 1;
    frmN = frameNum;
    stgPos = stagePos;
    imgTm = imgTime;
    tmp = temperature;
    gain = gain;
    pcenterC = centerC;
    pcenterR = centerR;
    imgNrm = imgNorm;
    imgDpth = imgDepth;

    for (uint ir=0; ir<(*imgSBkg).size(); ir++) {
      for (uint ic=0; ic<(*imgSBkg)[ir].size(); ic++) {
        imgSBprocOrig[ir][ic] = (*img)[ir][ic] - imgRat*bkgImg[ir*1024 + ic];
        imgprocOrig[ir][ic] = (*img)[ir][ic];
        imgBkgprocOrig[ir][ic] = (*imgBkg)[ir][ic];
        imgOrigprocOrig[ir][ic] = (*imgOrig)[ir][ic];
      }
    }

    // Rotate by PI/2
    //imgSBproc = imgProc::imgRotatePId2(imgSBprocOrig);
    //imgproc = imgProc::imgRotatePId2(imgprocOrig);
    //imgBkgproc = imgProc::imgRotatePId2(imgBkgprocOrig);
    //imgOrigproc = imgProc::imgRotatePId2(imgOrigprocOrig);
    imgSBproc = imgSBprocOrig;
    imgproc = imgprocOrig;
    imgBkgproc = imgBkgprocOrig;
    imgOrigproc = imgOrigprocOrig;

    int NhistBins = 70;
    int Nslices = 60;
    pimgSBproc = imgProc::polarBinning(imgSBproc,
                                 (int)centerR, (int)centerC, 
                                 Nslices, (int)centerR, NhistBins);
    pimgproc = imgProc::polarBinning(imgproc,
                                 (int)centerR, (int)centerC, 
                                 Nslices, (int)centerR, NhistBins);
    plt.printRC(pimgSBproc, "testImg"+to_string(imgN));

    tr->Fill();
  }

  tr->Write();
  outFile->Close();

  return 1;
}
