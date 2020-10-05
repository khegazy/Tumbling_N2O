#include "scanLegendre.h" 


using namespace std;


int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }

cout<<"starting"<<endl;
  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName
cout<<"made analysis class"<<endl;


  ///// Load environment and get the number of events /////
  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup
cout<<"setup Env"<<endl;



  
  int order = 30;
  int scanN = scanNum;
  int loopN = loopNum;
  int Nimgs = 0;
  int Nrows = (*img).size();
  int Ncols = (*img)[0].size();
  float scanNorm = 0;
  vector< vector< vector<double> > > lgndrs(order+1);  		// [order][img#][rad]
  vector< vector< vector<double> > > lgndrsOrd(order+1);  	// [order][img#][rad]
  vector< vector<double> > img_flat;
  vector< vector<double> > img_flatOrd;
  vector<double> imgNorms;
  vector<double> tvect;
  vector<int> imgNs;
  vector<double> imgLgndr((*pimg).size());
  int index = 0;

cout<<"doing options"<<endl;
  string fileName = "legendreCoeffs.root";
  string outputDir = "";
  for (int iarg=2; iarg<argc; iarg+=2) {
cout<<argv[iarg]<<endl;
    if (strcmp(argv[iarg],"-Ofile")==0) {string str(argv[iarg+1]); fileName=str;}
    else if (strcmp(argv[iarg],"-Odir")==0) {string str(argv[iarg+1]); outputDir=str;}
    else {
      cerr<<"ERROR!!! Option "<<argv[iarg]<<" does not exist!"<<endl;
      exit(0);
    }
  }

cout<<"making output file: "<<outputDir+fileName+".root"<<endl;
  SAVEclass save(outputDir+fileName+".root", "physics");
  save.tree->Branch("scanNum", &scanN, "scanNum/I");
  save.tree->Branch("loopNum", &loopN, "loopNum/I");
  save.tree->Branch("scanNorm", &scanNorm, "scanNorm/F");
  save.tree->Branch("imgNorms", &imgNorms);
  save.tree->Branch("img_flat", &img_flatOrd);
  save.tree->Branch("Ncols", &Ncols, "Ncols/I");
  save.tree->Branch("Nrows", &Nrows, "Nrows/I");
  for (int io=0; io<=order; io++) save.tree->Branch(("lgndr"+to_string(io)).c_str(), &lgndrsOrd[io]);

  vector<double> flatImg((*img).size()*(*img)[0].size());
  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);

/*
    imgNs.push_back(imageNum);
    for (int io=0; io<=order; io++) {
      // Looping over radii
      for (uint ir=0; ir<(*pimg).size(); ir++) imgLgndr[ir] = imgProc::legendre1dCoeff((*pimg)[ir], io);
      lgndrs[io].push_back(tvect);
      lgndrs[io][ievt].assign(imgLgndr.begin(), imgLgndr.end());  
    }
*/
/*
    for (uint ir=0; ir<(*pimg).size(); ir++) imgLgndr[ir] = imgProc::legendre1dCoeff((*pimg)[ir], 0);
    lgndrs[0].push_back(tvect);
    lgndrs[0][ievt].assign(imgLgndr.begin(), imgLgndr.end());  
*/

    scanNorm += imgNorm;
    Nimgs++;

    if (imgNorm == 1) {
      double tempNorm = 0;
      for (int ir=0; ir<(*pimg).size(); ir++) {
 	for (int ic=0; ic<(*pimg)[ir].size(); ic++) {
	  tempNorm += (*pimg)[ir][ic];
	}
      }
      for (int ir=0; ir<(*pimg).size(); ir++) {
        for (int ic=0; ic<(*pimg)[ir].size(); ic++) {
          (*pimg)[ir][ic] /= tempNorm;
        }
      }
    }

    for (int ir=0; ir<(*img).size(); ir++) {
      for (int ic=0; ic<(*img)[ir].size(); ic++) {
        flatImg[ir*(*img)[ir].size() + ic] = (*img)[ir][ic];
      }
    }
    img_flat.push_back(flatImg);
    imgNorms.push_back(imgNorm);

    imgNs.push_back(imageNum);
    for (int io=0; io<=order; io++) {
      // Looping over radii
      for (uint ir=0; ir<(*pimg).size(); ir++) {
if (io==2 && ievt==3 && ir==30) { 
//for (int k=0; k<(*pimg)[ir].size(); k++) {cout<<(*pimg)[ir][k]<<endl;}
//cout<<"imgnorm: "<<imgNorm<<endl;
//cout<<imgProc::legendre1dCoeff((*pimg)[ir], io)/double(imgNorm)<<endl;
}
      imgLgndr[ir] = imgProc::legendre1dCosCoeff((*pimg)[ir], io)/imgNorm;
}
      lgndrs[io].push_back(tvect);
      lgndrs[io][ievt].assign(imgLgndr.begin(), imgLgndr.end());  
    }
   
    if (scanN != scanNum || loopN != loopNum || ievt == (Nentries-1)) {
      // Time ordering the images 
      img_flatOrd.resize(img_flat.size());
      cout<<"time ordering"<<endl;
      for (int io=0; io<=order; io++) {
cout<<"leg: "<<io<<endl;
	lgndrsOrd[io].resize(imgNs.size());
cout<<"resized"<<endl;
        for (uint im=0; im<imgNs.size(); im++) {
          cout<<"img "<<im; cout<<"  "<<imgNs[im]<<endl; 
          lgndrsOrd[io][imgNs[im]-1].assign(lgndrs[io][im].begin(), lgndrs[io][im].end()); 
          cout<<"ordered"<<endl;
          if (io == 0) {
            img_flatOrd[imgNs[im]-1].resize(img_flat[im].size(), 0);
            img_flatOrd[imgNs[im]-1].assign(img_flat[im].begin(), img_flat[im].end());
          }
        }
        //if (loopN%2 == 0) {reverse(lgndrsOrd[io].begin(), lgndrsOrd[io].end()); cout<<"reversed"<<endl;} // Even loops start from end
        lgndrs[io].clear();	// Need to clear this vector since using push_back
        cout<<"cleared"<<endl;
	//plt.printXY(lgndrsOrd[io], "RSLegendre_"+to_string(io));
        cout<<"plotted"<<endl;
      }
cout<<"ordered"<<endl;
      
      //plt.printXY(lgndrsOrd[2], "legendr2_"+to_string(index)); 
 
      scanNorm /= float(Nimgs);

      save.fillEvent();
      scanN = scanNum;
      loopN = loopNum;
      scanNorm = 0;
      Nimgs = 0;
      imgNorms.clear();
      img_flat.clear();
    }

  }
  save.writeFile();

  return 1;
}
