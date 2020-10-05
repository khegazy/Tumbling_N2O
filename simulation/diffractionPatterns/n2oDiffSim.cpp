#include <chrono>
#include "/reg/neh/home/khegazy/baseTools/simulation/diffractionSimulation/simulationTools/atomClass.h"
#include "/reg/neh/home/khegazy/baseTools/simulation/diffractionSimulation/simulationTools/moleculeClass.h"
#include "/reg/neh/home/khegazy/baseTools/simulation/diffractionSimulation/simulationTools/molEnsembleMC.h"
#include "/reg/neh/home/khegazy/baseTools/simulation/diffractionSimulation/simulationTools/diffractionClass.h"
#include "/reg/neh/home/khegazy/baseTools/tools/plotClass.h"
#include "/reg/neh/home/khegazy/baseTools/tools/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseTools/tools/saveClass.h"

using namespace std;

int main(int argc, char* argv[]) {

  PLOTclass plt;

  double seed = (double)clock();
  int Nmols = 500;
  int index = 0;
  std::string time="null";

  string pdfFile = "/reg/neh/home/khegazy/simulations/n2o/rotation/output/alignmentPDFs/job10_Temp-60.3_Intns-2.72_Time-356.726.root";
  //string pdfFile = "/reg/neh/home/khegazy/simulations/n2o/rotation/output/alignment/job39_357950.000000.root";
  string fileName = "n2oDiffSim";
  string outputDir = "/reg/d/psdm/amo/amoi0314/scratch/diffPatterns/";
  //string outputDir = "./";
  //string outputDir = "/reg/neh/home/khegazy/simulations/n2o/diffractionPatterns/output";
  if (argc > 1) {
    Nmols = atoi(argv[1]);
  }
  if (argc > 2) {
    for (int iarg=2; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg],"-Ofile")==0) {string str(argv[iarg+1]); fileName=str;}
      else if (strcmp(argv[iarg],"-PDF")==0) {string str(argv[iarg+1]); pdfFile=str;}
      else if (strcmp(argv[iarg],"-Odir")==0) {string str(argv[iarg+1]); outputDir=str;}
      else if (strcmp(argv[iarg],"-Index")==0) {index = atoi(argv[iarg+1]);}
      else if (strcmp(argv[iarg],"-Time")==0) {
        std::string time_(argv[iarg+1]);
        time = time_;
        //pdfFile = "/reg/d/psdm/amo/amoi0314/scratch/anglePDFs/UED/anglePDFs_Time-"
        pdfFile = "/reg/ued/ana/scratch/n2o/alignmentPDFs/UED/anglePDFs_Time-"
          + time + "_ThetaBins-200_PhiBins-200.root";
      }
      else {
        cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
        exit(0);
      }
    }
  }

  cout << "\n\nINFO: Number of molecules is now " << Nmols << "!!!\n";
  cout << "INFO: Using pdf: " << pdfFile << "!!!\n";
  
  double elEnergy = 3.7e6;
  double lambda = 2*PI*C_AU/sqrt(pow(elEnergy*eV_to_au + C_AU*C_AU,2) - pow(C_AU,4)); //au
  lambda /= angs_to_au;  // angs
  double k0 = 2*PI/lambda;
  cout<<"lambda/k0: "<<lambda<<"  "<<k0<<endl;

  MOLENSEMBLEMCclass n2oMC(seed, "N2O.xyz", pdfFile.c_str(),"thetaPhiPDF");

  n2oMC.Nmols = Nmols;
  cout<<"SET N MOLS"<<endl;
  //n2oMC.NmolAtoms = 3;
  //n2oMC.atomTypes.push_back(N);
  //n2oMC.atomTypes.push_back(O);


  n2oMC.useOrientationMC = true;
  n2oMC.orientationPDF = "thetaPhiPDF";

  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  n2oMC.makeMolEnsemble();


  uint ipix = 1024;
  //DIFFRACTIONclass diffP(&n2oMC, 5, 4, 3.7e6,
  //      "/reg/neh/home5/khegazy/simulations/scatteringAmplitudes/3.7MeV/",
  //      ipix, -0.024, 0.024, ipix, -0.024, 0.024);
  DIFFRACTIONclass diffP("name", &n2oMC, 14.1207, 5, 4, 3.7e6, ipix,
        "/reg/neh/home5/khegazy/baseTools/simulation/scatteringAmplitudes/3.7MeV/");

  cout<<"Starting diffraction pattern calculation."<<endl;
  diffP.diffPatternCalc();
  /*
  for (uint ir=0; ir<diffP.diffMolPattern.size(); ir++) {
    for (uint ic=0; ic<diffP.diffMolPattern[ir].size(); ic++) {

      if (true || abs(diffP.diffMolPattern[ir][ic]) > 1e10) {
        cout<<"TOOOO LARGE: "<<ir<<" / "<<ic<<"      "<<diffP.diffMolPattern[ir][ic]<<endl;
      }
    }
  }
  */

  std::cout << "Time difference = " 
      << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() 
      << "[s]" << std::endl;

  // Rotating the image PI/2
  diffP.diffPattern = imgProc::imgRotatePId2(diffP.diffPattern);
  diffP.diffMolPattern = imgProc::imgRotatePId2(diffP.diffMolPattern);
  diffP.diffAtmPattern = imgProc::imgRotatePId2(diffP.diffAtmPattern);

  if (time.compare("null") != 0) {
    fileName += "_time-" + time;
    cout<<"SAVING HERE"<<endl;
    save::saveDat<double>(
        diffP.diffMolPattern,
        //"/reg/d/psdm/amo/amoi0314/scratch/diffPatterns/molDiff_time-"
        "/reg/ued/ana/scratch/n2o/diffPatterns/molDiff_time-"
        + time + "_bins["
        + to_string(diffP.diffMolPattern.size()) + ","
        + to_string(diffP.diffMolPattern[0].size()) + "].dat");
    save::saveDat<double>(
        diffP.diffAtmPattern,
        //"/reg/d/psdm/amo/amoi0314/scratch/diffPatterns/atmDiff_time-"
        "/reg/ued/ana/scratch/n2o/diffPatterns/atmDiff_time-"
        + time + "_bins["
        + to_string(diffP.diffAtmPattern.size()) + ","
        + to_string(diffP.diffMolPattern[0].size()) + "].dat");
  }



  vector<PLOToptions> opts(2);
  vector<string> vals(2);
  //opts[0] = xSpan;      vals[0] = to_string(diffP.sMinX) + "," + to_string(diffP.sMaxX);
  //opts[0] = ySpan;      vals[0] = to_string(diffP.sMinZ) + "," + to_string(diffP.sMaxZ);
  opts[0] = xSpan;      vals[0] = to_string(-diffP.sMax) + "," + to_string(diffP.sMax);
  opts[0] = ySpan;      vals[0] = to_string(-diffP.sMax) + "," + to_string(diffP.sMax);

  //plt.printRC(diffP.diffPattern, fileName + "diffractionPattern", opts, vals);
  //plt.printRC(diffP.diffMolPattern, fileName + "molDiffractionPattern", opts, vals);
  //plt.printRC(diffP.diffAtmPattern, fileName + "atmDiffractionPattern", opts, vals);

  int scanN = 0;
  int loopN = 0;
  float scanNorm = 1;
  int frmN = 0;
  int imgTm = 0;
  float tmp = 0;
  float gain = 0;
  int centerC = ipix/2;
  int centerR = ipix/2;
  float imgNorm = 1;
  int imgDepth = 0;
  
  cout<<pdfFile<<endl;
  cout<<"H2 "<<endl;
  size_t ipos = pdfFile.find("_Time") + 6; //pdfFile.find("job") + 3;
  size_t fpos = pdfFile.find("_", ipos);
  float stgPs = stof(pdfFile.substr(ipos, fpos - ipos))*(C_SI*1e-15*1e3);
  //float stgPs = stof(pdfFile.substr(ipos, fpos - ipos))/(C_SI*1e-15*1e3);
  cout<<"F2 "<<stgPs<<endl;
  
  vector< vector<double> > pimg = imgProc::polarBinning(diffP.diffPattern,
					 centerR, centerC, 80, ipix/2, 40);
  vector< vector<double> > pimgSB = imgProc::polarBinning(diffP.diffMolPattern,
					 centerR, centerC, 80, ipix/2, 40);




  TFile* fl = TFile::Open((outputDir + "/" + fileName + ".root").c_str(), "RECREATE");
  TTree* tr = new TTree("physics", "physics");

  tr->Branch("scanNum", &scanN, "scanNum/I");
  tr->Branch("loopNum", &loopN, "loopNum/I");
  tr->Branch("scanNorm", &scanNorm, "scanNorm/F");
  tr->Branch("imageNum", &index, "imageNum/I");
  tr->Branch("frameNum", &frmN ,"frameNum/I");
  tr->Branch("stagePos", &stgPs ,"stagePos/F");
  tr->Branch("imgTime", &imgTm ,"imgTime/I");
  tr->Branch("temperature", &tmp ,"temperature/F");
  tr->Branch("gain", &gain ,"gain/F");
  tr->Branch("centerC", &centerC ,"centerC/I");
  tr->Branch("centerR", &centerR ,"centerR/I");
  tr->Branch("imgNorm", &imgNorm ,"imgNorm/F");
  tr->Branch("imgDepth", &imgDepth ,"imgDepth");
  tr->Branch("imgOrig", &diffP.diffPattern);
  tr->Branch("img", &diffP.diffPattern);
  tr->Branch("imgBkg", &diffP.diffAtmPattern);
  tr->Branch("imgSBkg", &diffP.diffMolPattern);
  tr->Branch("pimg", &pimg);
  tr->Branch("pimgSB", &pimgSB);
  tr->Branch("sPattern", &diffP.sPattern);

  tr->Fill();
  tr->Write();
  //n2oMC.orientGraph->Write();
  //n2oMC.orientDist->Write();
  fl->Close();


  // Make atomic scattering for analysis use

  /*
  for (uint ir=0; ir<diffP.diffMolPattern.size(); ir++) {
    for (uint ic=0; ic<diffP.diffMolPattern[ir].size(); ic++) {

      if (true || abs(diffP.diffMolPattern[ir][ic]) > 1e10) {
        cout<<"TOOOO LARGE: "<<ir<<" / "<<ic<<"      "<<diffP.diffMolPattern[ir][ic]<<endl;
      }
    }
  }
  */
  plt.printRC(diffP.diffMolPattern, "testMolDiff_"+time);

  return 0;
}
