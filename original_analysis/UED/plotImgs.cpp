#include "plotImgs.h" 


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

  PLOTclass plt;
  TH2F* h;


  ///// Loop through events in the file /////
  for (uint64_t ievt=0; ievt<Nentries; ievt++) {
    analysis.loadEvent(ievt);


    h = plt.plotRC((*pimgBA), "pimgBA_"+to_string(ievt));
    //h->SetMaximum(2000);
    //h->SetMinimum(-2000);
    plt.print2d(h, "pimBAg_"+to_string(ievt));
    delete h;
  }

  return 1;
}
