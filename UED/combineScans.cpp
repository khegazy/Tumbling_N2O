#include "combineScans.h" 
#include <boost/math/special_functions/legendre.hpp>


using namespace std;


class peakfnctr {

    public:
	vector< vector<double> > *img;
	double sign=1;
	double operator() (vector<double> vect) {

          double sum=0;
          for (int itt=-3; itt<=3; itt++) {
            for (int irr=-3; irr<=3; irr++) sum += (*img)[vect[0]+irr][vect[1]+itt];
          }
	  return sum*sign;
	}
};



int main(int argc, char* argv[]) {

  if (argc<2) {
    cerr<<"ERROR: Missing input arguments, must run code ./analysis.exe 'fileList.txt' !!!"<<endl;
    cerr<<"         Can also run ./analysis.exe 'fileList.txt' 'treeName'"<<endl;
    exit(0);
  }

  
  /////  Flags  /////
  bool doMCbasis = true;

  string fileList(argv[1]);
  string treeName("physics");
  if (argc==3) string treeName(argv[2]);
  analysisClass analysis(fileList, treeName);  // Use this when specifying treeName


  ///// Load environment and get the number of events /////
  uint64_t Nentries;
  Nentries = analysis.setupEnvironment();  // Alter this function as needed for specific setup
  int dataMC;
  bool bendSearch = false;
  uint bendTimeItr = 21;
  string bendAng;
  TFile* output;
  ofstream outTxt;
  if (fileList.find("MC") != string::npos) {
    dataMC = 1;
    string tfname = "combinedResults";
    if (bendSearch) {
      tfname = "bendSearchResults";
    }
    size_t ipos = fileList.find("Bend");
    if (ipos != string::npos) {
      size_t fpos = fileList.find(".txt");
      bendAng = fileList.substr(ipos + 4, fpos - (ipos + 4));
      plt.setFolder("plots/MC/bend" + fileList.substr(ipos + 4, fpos - (ipos + 4)));
      output = TFile::Open(("output/MC/bend" 
			+ fileList.substr(ipos + 4, fpos - (ipos + 4)) 
			+ "/" + tfname + ".root").c_str(), "RECREATE");
      outTxt.open("output/MC/bend/combineOutput.txt");
    }
    else {
      plt.setFolder("plots/MC");
      output = TFile::Open("output/MC/combinedResults.root", "RECREATE");
      outTxt.open("output/MC/combineOutput.txt");
    }
  }
  else {
    dataMC = 0;
    bendSearch = false;
    plt.setFolder("plots/data");
    output = TFile::Open("output/data/combinedResults.root", "RECREATE");
    outTxt.open("output/data/combineOutput.txt");
  }


  // q calibration
  double maxQ;
  if (dataMC) {
    maxQ = 14.1207;
  }
  else {
    double qPer1024Pix = 0.028156*(137.0/112.0); //(26.5/22.0); //0.0358;
    uint NorigPix = 820;
    maxQ = NorigPix*qPer1024Pix/2.0;
  }


  int it, ir, io, ic, mIter, maxT, MAXT, MINT, shift, cnt;
  double sum, maxSum, minSum, slope, intcpt, scl, sbtr;

  int Nlg=31;
  int NradAvgs = 3;
  int irr, icc;
  vector< vector< vector<double> > > lgndrs(Nlg);	// io : q (scattering) : time
  vector< vector< vector< vector<double> > > > lgndr_vals(Nlg);	// io : q (scattering) : time
  vector< vector< vector< vector<double> > > > lgndr_scales(Nlg);	// io : q (scattering) : time
  vector< vector< vector<double> > > indVals(NradAvgs);	// Used to calculate the variance
  vector< vector<double> > indScls((*lgndr2).size());	// Used to calculate the variance
  vector< vector< vector<double> > > origDiffP;
  vector<double> count((*lgndr2).size()); 
  vector<double> origDiffPcount((*lgndr2).size(), 0); 
  vector<double> sumComb;
  vector<double> sumSing;
  vector< vector<double> > mask((*lgndr2)[0].size());
  for (ir=0; ir<(int)(*lgndr2)[0].size(); ir++) mask[ir].resize(21, 0.0); 
  for (ir=0; ir<(int)indVals.size(); ir++) indVals[ir].resize((*lgndr2).size());
  //vector< vector< vector<double> > > lgndrs(Nlg);

  /////  Aligning and combining scans by correlating with a mask  /////
  for (int itr=0; itr<3; itr++) {
    analysis.loadEvent(0);
    for (io=0; io<Nlg; io++) {
      lgndrs[io].clear();       lgndrs[io].resize((*lgndr2)[0].size());
      lgndr_vals[io].clear();   lgndr_vals[io].resize((*lgndr2)[0].size());
      lgndr_scales[io].clear(); lgndr_scales[io].resize((*lgndr2)[0].size());
      for (ir=0; ir<(int)(*lgndr2)[0].size(); ir++) {
        lgndrs[io][ir].clear();         lgndrs[io][ir].resize((*lgndr2).size(),0.0);
        lgndr_vals[io][ir].clear();     lgndr_vals[io][ir].resize((*lgndr2).size());
        lgndr_scales[io][ir].clear();   lgndr_scales[io][ir].resize((*lgndr2).size());
      }
    }
    count.clear();	count.resize((*lgndr2).size(),0.0);
    for (ir=0; ir<(int)indVals.size(); ir++) {
      indVals[ir].resize((*lgndr2).size());
    }
    indScls.resize((*lgndr2).size());
    if (itr == 2) {
      origDiffP.resize(Nrows);
      for (irr=0; irr<Nrows; irr++) {
        origDiffP[irr].resize(Ncols);
        for (icc=0; icc<Ncols; icc++) {
          origDiffP[irr][icc].resize((*lgndr2).size(), 0);
        }
      }
      origDiffPcount.resize((*lgndr2).size(), 0);
    }
    // Start with the first entry
    for (ir=0; ir<(int)(*lgndr2)[0].size(); ir++) {
      for (int k=0; k<Nlg; k++) {
        lgndrs[k][ir].clear();	      lgndrs[k][ir].resize((*lgndr2).size(),0.0);
        lgndr_vals[k][ir].clear();    lgndr_vals[k][ir].resize((*lgndr2).size());
        lgndr_scales[k][ir].clear();  lgndr_scales[k][ir].resize((*lgndr2).size());
      }
      for (it=0; it<(int)(*lgndr2).size(); it++) {
        count[it] = scanNorm;
        lgndr_vals[0][ir][it].push_back((*lgndr0)[it][ir]);
        lgndr_vals[1][ir][it].push_back((*lgndr1)[it][ir]);
        lgndr_vals[2][ir][it].push_back((*lgndr2)[it][ir]);
        lgndr_vals[3][ir][it].push_back((*lgndr3)[it][ir]);
        lgndr_vals[4][ir][it].push_back((*lgndr4)[it][ir]);
        lgndr_vals[5][ir][it].push_back((*lgndr5)[it][ir]);
        lgndr_vals[6][ir][it].push_back((*lgndr6)[it][ir]);
        lgndr_vals[7][ir][it].push_back((*lgndr7)[it][ir]);
        lgndr_vals[8][ir][it].push_back((*lgndr8)[it][ir]);
        lgndr_vals[9][ir][it].push_back((*lgndr9)[it][ir]);
        lgndr_vals[10][ir][it].push_back((*lgndr10)[it][ir]);
        lgndr_vals[11][ir][it].push_back((*lgndr11)[it][ir]);
        lgndr_vals[12][ir][it].push_back((*lgndr12)[it][ir]);
        lgndr_vals[13][ir][it].push_back((*lgndr13)[it][ir]);
        lgndr_vals[14][ir][it].push_back((*lgndr14)[it][ir]);
        lgndr_vals[15][ir][it].push_back((*lgndr15)[it][ir]);
        lgndr_vals[16][ir][it].push_back((*lgndr16)[it][ir]);
        lgndr_vals[17][ir][it].push_back((*lgndr17)[it][ir]);
        lgndr_vals[18][ir][it].push_back((*lgndr18)[it][ir]);
        lgndr_vals[19][ir][it].push_back((*lgndr19)[it][ir]);
        lgndr_vals[20][ir][it].push_back((*lgndr20)[it][ir]);
        lgndr_vals[21][ir][it].push_back((*lgndr21)[it][ir]);
        lgndr_vals[22][ir][it].push_back((*lgndr22)[it][ir]);
        lgndr_vals[23][ir][it].push_back((*lgndr23)[it][ir]);
        lgndr_vals[24][ir][it].push_back((*lgndr24)[it][ir]);
        lgndr_vals[25][ir][it].push_back((*lgndr25)[it][ir]);
        lgndr_vals[26][ir][it].push_back((*lgndr26)[it][ir]);
        lgndr_vals[27][ir][it].push_back((*lgndr27)[it][ir]);
        lgndr_vals[28][ir][it].push_back((*lgndr28)[it][ir]);
        lgndr_vals[29][ir][it].push_back((*lgndr29)[it][ir]);
        lgndr_vals[30][ir][it].push_back((*lgndr30)[it][ir]);
        lgndr_scales[0][ir][it].push_back(scanNorm);
        lgndrs[0][ir][it] = (*lgndr0)[it][ir];
        lgndrs[1][ir][it] = (*lgndr1)[it][ir];
        lgndrs[2][ir][it] = (*lgndr2)[it][ir];
        lgndrs[3][ir][it] = (*lgndr3)[it][ir];
        lgndrs[4][ir][it] = (*lgndr4)[it][ir];
        lgndrs[5][ir][it] = (*lgndr5)[it][ir];
        lgndrs[6][ir][it] = (*lgndr6)[it][ir];
        lgndrs[7][ir][it] = (*lgndr7)[it][ir];
        lgndrs[8][ir][it] = (*lgndr8)[it][ir];
        lgndrs[9][ir][it] = (*lgndr9)[it][ir];
        lgndrs[10][ir][it] = (*lgndr10)[it][ir];
        lgndrs[11][ir][it] = (*lgndr11)[it][ir];
        lgndrs[12][ir][it] = (*lgndr12)[it][ir];
        lgndrs[13][ir][it] = (*lgndr13)[it][ir];
        lgndrs[14][ir][it] = (*lgndr14)[it][ir];
        lgndrs[15][ir][it] = (*lgndr15)[it][ir];
        lgndrs[16][ir][it] = (*lgndr16)[it][ir];
        lgndrs[17][ir][it] = (*lgndr17)[it][ir];
        lgndrs[18][ir][it] = (*lgndr18)[it][ir];
        lgndrs[19][ir][it] = (*lgndr19)[it][ir];
        lgndrs[20][ir][it] = (*lgndr20)[it][ir];
        lgndrs[21][ir][it] = (*lgndr11)[it][ir];
        lgndrs[22][ir][it] = (*lgndr12)[it][ir];
        lgndrs[23][ir][it] = (*lgndr13)[it][ir];
        lgndrs[24][ir][it] = (*lgndr14)[it][ir];
        lgndrs[25][ir][it] = (*lgndr15)[it][ir];
        lgndrs[26][ir][it] = (*lgndr16)[it][ir];
        lgndrs[27][ir][it] = (*lgndr17)[it][ir];
        lgndrs[28][ir][it] = (*lgndr18)[it][ir];
        lgndrs[29][ir][it] = (*lgndr19)[it][ir];
        lgndrs[30][ir][it] = (*lgndr20)[it][ir];
	if (ir == 13) {
          for (irr=0; irr<(int)indVals.size(); irr++) {
	    indVals[irr][it].push_back((*lgndr2)[it][9 + irr]);
          }
	  indScls[it].push_back(scanNorm);
	}
      }
    }

    if (itr == 2 && (Nentries < 2)) {
      for (it=0; it<(int)(*lgndr2).size(); it++) {
        for (irr=0; irr<Nrows; irr++) {
          for (icc=0; icc<Ncols; icc++) {
            origDiffP[irr][icc][it] += (*img_flat)[it][irr*Ncols + icc]/(*imgNorms)[it];
          }
        }
      }
    }

    /// Loop through events to align and average together ///
    for (uint64_t ievt=1; ievt<Nentries; ievt++) {
      analysis.loadEvent(ievt);

      MAXT = MINT = -1;

      //  Maximizing the autocorrelation of legendres
      // If making the first mask drop pixels > 1000 
      if (itr==0) {
        sumComb.clear();  sumComb.resize(lgndrs[2][0].size(),0.0);
        sumSing.clear();  sumSing.resize((*lgndr2).size(),0.0);    
  
        for (it=0; it<(int)(*lgndr2).size(); it++) {
          sumSing[it]=0;
          for (ir=0.2*(*lgndr2)[it].size(); ir<=0.4*(*lgndr2)[it].size(); ir++) {
            if ((*lgndr2)[it][ir]<1000) {
              sumSing[it] += (*lgndr2)[it][ir];
            }
          }
        }
      
        for (it=0; it<(int)lgndrs[2][0].size(); it++) {
          sumComb[it]=0;
          for (ir=0.2*(lgndrs[2].size()); ir<=0.4*(lgndrs[2].size()); ir++) {
            if (lgndrs[2][ir][it]<1000) {
              sumComb[it] += lgndrs[2][ir][it];
            }
          }
        }
  
        // Maximize the correlation to find best alignment
        maxSum=-9999;    
        for (int k=-10; k<=10; k++) {
          sum=0;
          //for (it=0; it<(int)sumSing.size(); it++) 
          for (it=0; it<(int)sumComb.size(); it++) {
            //if (it>(int)sumComb.size()-1 || it<0 || k+it>=(int)sumSing.size() || k+it<=0) continue;
            if (k+it>=(int)sumSing.size() || k+it<0) continue;
            sum += sumComb[it]*sumSing[k+it];
          }
          if (sum > maxSum) {
            maxSum = sum;
            shift = k;
          }
        }
      }
      // Already made a mask
      else {
	// Maximizing correlation of mask with saved data
        maxSum=-9999;
        for (int k=10; k<=25; k++) {
          sum=0;
          for (it=0; it<(int)mask[0].size(); it++) {
            if (k+it>=(int)(*lgndr2).size() || k+it<=0) continue;
            for (ir=(int)(1.0+(*lgndr2)[0].size()/8.0); ir<(int)(*lgndr2)[0].size(); ir++) {
              sum += (*lgndr2)[k+it][ir]*mask[ir][it];
            }
          }
          if (sum > maxSum) {
            maxT=k; 
            maxSum=sum;
          }
        }
	// Maximizing correlation with aligned and averaged scans
        maxSum=-9999;
        for (int k=10; k<=25; k++) {
          sum=0;
          for (it=0; it<(int)mask[0].size(); it++) {
            if (k+it>=(int)lgndrs[2][0].size() || k+it<=0) continue;
            for (ir=(int)(1.0+(*lgndr2)[0].size()/8.0); ir<(int)(*lgndr2)[0].size(); ir++) {
              sum += lgndrs[2][ir][k+it]*mask[ir][it];
            }
          }
          if (sum > maxSum) {
            MAXT=k; 
            maxSum=sum;
          }
        }
        shift = MAXT-maxT;
      }  
 
 
      mIter = shift;
      // Add column to end if needed
      if ((mIter+(*lgndr2).size())>lgndrs[0][0].size()) {
        int size = (int)lgndrs[2][0].size();
        std::vector<double> eVec;
	//cout<<"add columns to end: "<<(mIter+(int)(*lgndr2).size() - size)<<endl;
        for (int j=0; j<(mIter+(int)(*lgndr2).size() - size); j++) {
          count.emplace_back(0);
	  indScls.emplace_back();
          for (irr=0; irr<(int)indVals.size(); irr++) {
	    indVals[irr].emplace_back();
          }
          for (int l=0; l<Nlg; l++) {
            for (ir=0; ir<(int)lgndrs[l].size(); ir++) {
              lgndrs[l][ir].emplace_back(0);
              lgndr_vals[l][ir].emplace_back(eVec);
              lgndr_scales[l][ir].emplace_back(eVec);
              
            }
          }
          if (itr == 2) {
            origDiffPcount.emplace_back(0);
            for (irr=0; irr<Nrows; irr++) {
              for (icc=0; icc<Ncols; icc++) {
                origDiffP[irr][icc].emplace_back(0);
              }
            }
          }

        }
      }
 
      // Add column to front if needed
      if (mIter<0) {
        std::vector<double> eVec;
        for (int j=0; j<fabs(mIter); j++) {
          count.emplace(count.begin(), 0);
	  indScls.emplace(indScls.begin());
          for (irr=0; irr<(int)indVals.size(); irr++) {
	    indVals[irr].emplace(indVals[irr].begin());
          }
          for (int l=0; l<Nlg; l++) {
            for (ir=0; ir<(int)lgndrs[l].size(); ir++) {
              lgndrs[l][ir].emplace(lgndrs[l][ir].begin(), 0);
              lgndr_vals[l][ir].emplace(lgndr_vals[l][ir].begin(), eVec);
              lgndr_scales[l][ir].emplace(lgndr_scales[l][ir].begin(), eVec);
            }
          }
          if (itr == 2) {
            origDiffPcount.emplace(origDiffPcount.begin(), 0);
            for (irr=0; irr<Nrows; irr++) {
              for (icc=0; icc<Ncols; icc++) {
                origDiffP[irr][icc].emplace(origDiffP[irr][icc].begin(), 0);
              }
            }
          }
        }
        MAXT += fabs(mIter);
        mIter=0;
      }
 
      // Averaging the scans together
      if (itr==0) {
        for (it=0; it<(int)(*lgndr2).size(); it++) {
          for (ir=0; ir<(int)(*lgndr2)[0].size(); ir++) {
            lgndrs[2][ir][mIter+it] = (*lgndr2)[it][ir] < 1000 ? lgndrs[2][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr2)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)) : lgndrs[2][ir][mIter+it]; 
          }
          count[mIter+it] += scanNorm;
        }
      }
      else {
        for (it=0; it<(int)(*lgndr2).size(); it++) {
 	  if (itr == 2) {
            for (irr=0; irr<(int)indVals.size(); irr++) {
	      indVals[irr][mIter+it].push_back((*lgndr2)[it][9 + irr]);
            }
	    indScls[mIter+it].push_back(scanNorm);
	  }
          for (ir=0; ir<(int)(*lgndr2)[0].size(); ir++) {
            lgndr_vals[0][ir][mIter+it].push_back((*lgndr0)[it][ir]);
            lgndr_vals[1][ir][mIter+it].push_back((*lgndr1)[it][ir]);
            lgndr_vals[2][ir][mIter+it].push_back((*lgndr2)[it][ir]);
            lgndr_vals[3][ir][mIter+it].push_back((*lgndr3)[it][ir]);
            lgndr_vals[4][ir][mIter+it].push_back((*lgndr4)[it][ir]);
            lgndr_vals[5][ir][mIter+it].push_back((*lgndr5)[it][ir]);
            lgndr_vals[6][ir][mIter+it].push_back((*lgndr6)[it][ir]);
            lgndr_vals[7][ir][mIter+it].push_back((*lgndr7)[it][ir]);
            lgndr_vals[8][ir][mIter+it].push_back((*lgndr8)[it][ir]);
            lgndr_vals[9][ir][mIter+it].push_back((*lgndr9)[it][ir]);
            lgndr_vals[10][ir][mIter+it].push_back((*lgndr10)[it][ir]);
            lgndr_vals[11][ir][mIter+it].push_back((*lgndr11)[it][ir]);
            lgndr_vals[12][ir][mIter+it].push_back((*lgndr12)[it][ir]);
            lgndr_vals[13][ir][mIter+it].push_back((*lgndr13)[it][ir]);
            lgndr_vals[14][ir][mIter+it].push_back((*lgndr14)[it][ir]);
            lgndr_vals[15][ir][mIter+it].push_back((*lgndr15)[it][ir]);
            lgndr_vals[16][ir][mIter+it].push_back((*lgndr16)[it][ir]);
            lgndr_vals[17][ir][mIter+it].push_back((*lgndr17)[it][ir]);
            lgndr_vals[18][ir][mIter+it].push_back((*lgndr18)[it][ir]);
            lgndr_vals[19][ir][mIter+it].push_back((*lgndr19)[it][ir]);
            lgndr_vals[20][ir][mIter+it].push_back((*lgndr20)[it][ir]);
            lgndr_vals[21][ir][mIter+it].push_back((*lgndr21)[it][ir]);
            lgndr_vals[22][ir][mIter+it].push_back((*lgndr22)[it][ir]);
            lgndr_vals[23][ir][mIter+it].push_back((*lgndr23)[it][ir]);
            lgndr_vals[24][ir][mIter+it].push_back((*lgndr24)[it][ir]);
            lgndr_vals[25][ir][mIter+it].push_back((*lgndr25)[it][ir]);
            lgndr_vals[26][ir][mIter+it].push_back((*lgndr26)[it][ir]);
            lgndr_vals[27][ir][mIter+it].push_back((*lgndr27)[it][ir]);
            lgndr_vals[28][ir][mIter+it].push_back((*lgndr28)[it][ir]);
            lgndr_vals[29][ir][mIter+it].push_back((*lgndr29)[it][ir]);
            lgndr_vals[30][ir][mIter+it].push_back((*lgndr30)[it][ir]);
            lgndr_scales[0][ir][mIter+it].push_back(scanNorm); 
            lgndrs[0][ir][mIter+it] = lgndrs[0][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr0)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[1][ir][mIter+it] = lgndrs[1][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr1)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[2][ir][mIter+it] = lgndrs[2][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr2)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[3][ir][mIter+it] = lgndrs[3][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr3)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[4][ir][mIter+it] = lgndrs[4][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr4)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[5][ir][mIter+it] = lgndrs[5][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr5)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[6][ir][mIter+it] = lgndrs[6][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr6)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[7][ir][mIter+it] = lgndrs[7][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr7)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[8][ir][mIter+it] = lgndrs[8][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr8)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[9][ir][mIter+it] = lgndrs[9][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr9)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm)); 
            lgndrs[10][ir][mIter+it] = lgndrs[10][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr10)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[11][ir][mIter+it] = lgndrs[11][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr11)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[12][ir][mIter+it] = lgndrs[12][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr12)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[13][ir][mIter+it] = lgndrs[13][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr13)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[14][ir][mIter+it] = lgndrs[14][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr14)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[15][ir][mIter+it] = lgndrs[15][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr15)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[16][ir][mIter+it] = lgndrs[16][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr16)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[17][ir][mIter+it] = lgndrs[17][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr17)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[18][ir][mIter+it] = lgndrs[18][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr18)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[19][ir][mIter+it] = lgndrs[19][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr19)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[20][ir][mIter+it] = lgndrs[20][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr20)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[21][ir][mIter+it] = lgndrs[21][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr21)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[22][ir][mIter+it] = lgndrs[22][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr22)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[23][ir][mIter+it] = lgndrs[23][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr23)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[24][ir][mIter+it] = lgndrs[24][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr24)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[25][ir][mIter+it] = lgndrs[25][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr25)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[26][ir][mIter+it] = lgndrs[26][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr26)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[27][ir][mIter+it] = lgndrs[27][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr27)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[28][ir][mIter+it] = lgndrs[28][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr28)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[29][ir][mIter+it] = lgndrs[29][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr29)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
            lgndrs[30][ir][mIter+it] = lgndrs[30][ir][mIter+it]*(count[mIter+it]/(count[mIter+it]+scanNorm)) + (*lgndr30)[it][ir]*(scanNorm/(count[mIter+it]+scanNorm));
          }
          if (itr == 2) {
            for (irr=0; irr<Nrows; irr++) {
              for (icc=0; icc<Ncols; icc++) {
                origDiffP[irr][icc][mIter+it] += (*img_flat)[it][irr*Ncols + icc]/(*imgNorms)[it];
              }
            }
            origDiffPcount[mIter+it] += 1;
          }
          count[mIter+it] += scanNorm;
        }

      }
    }

    if (itr == 2) {
      for (it=0; it<origDiffPcount.size(); it++) {
        for (irr=0; irr<Nrows; irr++) {
          for (icc=0; icc<Ncols; icc++) {
            origDiffP[irr][icc][it] /= origDiffPcount[it];
          }
        }
      }
    }


    // Create Mask
    ir = lgndrs[2].size()/4;
    sum=0; maxSum=-999999;
    for (it=10; it<(int)lgndrs[2][0].size()-10; it++) {
      sum=0;
      for (int itt=-3; itt<=3; itt++) {
        for (irr=-3; irr<=3; irr++) sum += lgndrs[2][ir+irr][it+itt];
      }
      if (sum > maxSum) {maxSum=sum; MAXT=it;}  
    }
    sum=0; minSum=999999;
    for (it=10; it<(int)lgndrs[2][0].size()-10; it++) {
      sum=0;
      for (int itt=-3; itt<=3; itt++) {
        for (irr=-3; irr<=3; irr++) sum += lgndrs[2][ir+irr][it+itt];
      }
      if (sum < minSum) {minSum=sum; MINT=it;}  
    }

    int itm=(MINT+MAXT)/2-10;
    for (it=0; it<(int)mask[0].size(); it++) {
      for (ir=(int)(1.0+(*lgndr2)[0].size()/8.0); ir<(int)(*lgndr2)[0].size(); ir++) {
        mask[ir][it]=lgndrs[2][ir][it+itm]; 
      }
    }
 
    plt.printRC(mask, "mask"+to_string(itr));
  }

  // Calculating Variance for specific rad (check)
  vector<double> varnc(lgndrs[2][10].size(), 0);
  for (it=0; it<lgndrs[2][10].size(); it++) {
    double anorm = 0;
    for (uint k=0; k<indVals[0][it].size(); k++) {
      for (int irr=0; irr<NradAvgs; irr++) {
        varnc[it] += pow(indVals[irr][it][k] - lgndrs[2][9 + irr][it], 2)
          *indScls[it][k]/((double)NradAvgs);
      }
      anorm += indScls[it][k];
    }
    varnc[it] /= (anorm*(indVals[0][it].size() - 1));
    //varnc[it] /= (indVals[0][it].size() - 1);
  }

  vector<double> sumlgndr2(lgndrs[2][10].size(), 0);
  for (ir=9; ir<9+NradAvgs; ir++) {
    for (it=0; it<sumlgndr2.size(); it++) {
      sumlgndr2[it] += lgndrs[2][ir][it]/((double)NradAvgs);
    }
  }

  FILE* dtFile =
      fopen(("../simulation/fitParameters/data/rawOrigDataP2_Size-"
      + to_string(sumlgndr2.size()) + ".dat").c_str(),
      "wb");
  fwrite(&sumlgndr2[0], sizeof(double), sumlgndr2.size(), dtFile);
  fclose(dtFile);


  // Subtracting DC
  for (int k=0; k<Nlg; k++) {
    for (ir=0; ir<(int)lgndrs[k].size(); ir++) {
      sum=0;
      for (it=0; it<(int)lgndrs[k][ir].size(); it++) sum += lgndrs[k][ir][it];
      sum /= lgndrs[k][ir].size();
      for (it=0; it<(int)lgndrs[k][ir].size(); it++) lgndrs[k][ir][it] -= sum;
    }
  }

/*
  vector<PLOToptions> vopts(2);
  vector<string> vvals(2);
  vopts[0] = maximum;
  vopts[1] = minimum;
  //vvals[0] = "5e-6";
  //vvals[1] = "-5e-6";
  vvals[0] = "0.004";
  vvals[1] = "-0.004";


  
  vector< vector<double> > tryvibe;
  vector<double> tvibe(lgndrs[0][0].size());
  for (ir=1; ir<(int)lgndrs[0].size(); ir+=2) {
    for (it=0; it<(int)lgndrs[0][ir].size(); it++) {
      tvibe[it] = lgndrs[0][ir][it] - lgndrs[0][ir-1][it];
    }
    tryvibe.push_back(tvibe);
    plt.print1d(tvibe, "vibLineOut" + to_string(ir), vopts, vvals);
  }

  plt.printRC(tryvibe, "vibrationStates", vopts, vvals);
  //plt.printRC(tryvibe, "vibrationStates");
*/


/*
  // Plotting beta coefficients
  for (il=1; il<=Nlg; il++) {
    for (ir=0; ir<(int)lgndrs[k].size(); ir++) {
      for (it=0; it<(int)lgndrs[k][ir].size(); it++) {
      betas[il][ir][it] = lgndrs[il][ir][it]/lgndrs[0][ir][it];
      }
    }
  }
*/

  // Saving Legendre Coefficients
  if (!dataMC) {
    save::saveDat<double>(lgndrs, 
        "./output/data/legendre_coefficients["
        + to_string(lgndrs.size()) + ","
        + to_string(lgndrs[0].size()) + ","
        + to_string(lgndrs[0][0].size()) + "].dat");
    cout<<"PRINTIN LGNDRS INFO\n\n";
    std::vector<double> tst;
    for (int i=0; i<40; i++) {
        cout<<lgndrs[2][15][i]<<endl;
        tst.push_back(lgndrs[2][15][i]);
    }
    plt.print1d(tst, "testL2plot");

    std::vector< std::vector< std::vector<double> > > lgndr_sem(Nlg);
    for (int ilg=0; ilg<Nlg; ilg++) {
      lgndr_sem[ilg].resize(lgndrs[ilg].size());
      for (ir=0; ir<(int)lgndrs[ilg].size(); ir++) {
        lgndr_sem[ilg][ir].resize(lgndrs[ilg][ir].size(), 0);
        for (it=0; it<(int)lgndrs[ilg][ir].size(); it++) {
          lgndr_sem[ilg][ir][it] = 0;
          double norm = 0;
          for (uint i=0; i<lgndr_vals[ilg][ir][it].size(); i++) {
            lgndr_sem[ilg][ir][it] += lgndr_scales[0][ir][it][i]
                *std::pow(lgndr_vals[ilg][ir][it][i] - lgndrs[ilg][ir][it], 2);
            norm += lgndr_scales[0][ir][it][i];
          }
          lgndr_sem[ilg][ir][it] /= norm;
          lgndr_sem[ilg][ir][it] /= lgndr_vals[ilg][ir][it].size();
          lgndr_sem[ilg][ir][it] = std::sqrt(lgndr_sem[ilg][ir][it]);
        }
      }
    }

    save::saveDat<double>(lgndr_sem, 
        "./output/data/legendre_coefficient_SEM["
        + to_string(lgndr_sem.size()) + ","
        + to_string(lgndr_sem[0].size()) + ","
        + to_string(lgndr_sem[0][0].size()) + "].dat");
 
  }

  cout<<"Number of time bins: "<<lgndrs[2][0].size()<<endl;

  // Smoothly send ends to 0
  double mult;
  for (int k=0; k<Nlg; k++) {
    for (int j=0; j<4; j++) {
      mult = pow(sin(((double)j/4.0)*PI/2.0),2);
      for (ir=0; ir<(int)lgndrs[k].size(); ir++)  {
   	lgndrs[k][ir][j] = mult*lgndrs[k][ir][4];
   	lgndrs[k][ir][lgndrs[k][ir].size()-1-j] = mult*lgndrs[k][ir][lgndrs[k][ir].size()-5-j];
      }
    }
  }

  /*
  // Scaling legendre coefficients for sms
  std::vector<double> atmDiffLgndr(lgndrs[0].size());
  save::importDat<double>(atmDiffLgndr,
      "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/newRefAnalysis/diffractionPattern/output/nitrousOxide_atmDiffractionPatternLineOut_Qmax-"
      + to_string(maxQ)
      + "_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins["
      + to_string(lgndrs[0].size()) + "].dat");
  
  for (int k=0; k<Nlg; k++) {
    for (ir=0; ir<(int)lgndrs[k].size(); ir++)  {
      for (int j=0; j<lgndrs[k][ir].size(); j++) {
        lgndrs[k][ir][j] *= maxQ*ir/(atmDiffLgndr[ir]*lgndrs[k].size());
      }
    }
  }
  */

  //// Plotting ////
  if (!bendSearch) {
    vector<double> Lmax(Nlg);
    vector<PLOToptions> opts(4);
    vector<string> valsl(4), valsr(4);
    vector<TH1*> hists(2);
    opts[0]=xLabel;		valsl[0]="Time [ps]";				valsr[0]="Time [ps]";
    opts[1]=yLabel;		valsl[1]="Q #left[ #AA^{-1} #right]";			valsr[1]="# Scans";
    opts[2]=xSpan;		valsl[2]="0,"+to_string(count.size()*0.1);  	valsr[2]="0,"+to_string(count.size()*0.1);
    opts[3]=ySpan;		valsl[3]="0,"+to_string(maxQ);  	valsr[3]="null";
  //  opts[4]=maximum;		valsl[4]="1e-7";					valsr[4]="50";
  //  opts[5]=minimum;		valsl[5]="-1e-7";					valsr[5]="-50";


    plt.doubleCanvas("double",600,600,'y',0.8);
    hists[1] = plt.plot1d(count, "scanAvg", opts, valsr);
    plt.print1d((TH1F*)hists[1], "AvgCount");
    for (int k=0; k<Nlg; k++) {
      hists[0] = plt.plotRC(lgndrs[k], "legendre"+to_string(k)+"Comb", opts, valsl);
      hists[0]->GetMaximum() > fabs(hists[0]->GetMinimum()) ? Lmax[k] = hists[0]->GetMaximum() : Lmax[k] = fabs(hists[0]->GetMinimum());
      plt.print2d((TH2F*)hists[0], string(hists[0]->GetName()))->Write();
      //plt.print2d((TH2F*)hists[0], string(hists[0]->GetName()), fileType, "ascii");
     // plt.print(hists, "legendre"+to_string(k)+"CombCount", pads, "double");

      hists[0]->Write();    
    }

    // Plotting the largest contributions of each Legendre
    plt.print1d(Lmax, "maxLegendreC", logy, "");
    
  }


  ///////////////////////////////////////////
  /////  Building the Rotational Basis  /////
  ///////////////////////////////////////////
  int npows = 1;
  int Nders = 2;  // 1 => no derivatives
  int NrotLg = 9;
  int id, ip;
  int N = lgndrs[0][0].size();
  // rotB[lgndr][der or not][pow][time]
  vector< vector< vector< vector<double> > > > rotB(NrotLg);

  for (uint k=0; k<rotB.size(); k++) {
    rotB[k].resize(Nders);
  }
  for (io=0; io<NrotLg; io++) {
    for (uint j=0; j<rotB[io].size(); j++) {
      rotB[io][j].resize(npows);
      for (int l=0; l<npows; l++) {
        rotB[io][j][l].resize(lgndrs[2*io][0].size(), 0.0);
      }
    }
  }


  ////  Finding the lowest order rotational basis  ////

  vector<double> anorms(rotB.size(), 0);

  for (io=0; io<NrotLg; io++) {
    cnt = 0;
    int strBin = (int)(lgndrs[2*io].size()/4.5);
    int endBin = (int)(lgndrs[2*io].size()/3.4);
    if (io >= 7) {
      strBin += io - 6;
      endBin += io - 6;
    }

    int shift = 2;
    for (ir=strBin; ir<endBin; ir++) {
      for (it=0; it<(int)lgndrs[shift + 2*io][ir].size(); it++) {
        anorms[io] += fabs(lgndrs[shift + 2*io][ir][it]);
      }
      for (it=0; it<(int)lgndrs[shift + 2*io][ir].size(); it++) {
        rotB[io][0][0][it] += lgndrs[shift + 2*io][ir][it];
      }
      cnt++;
    }
    for (it=0; it<(int)rotB[io][0][0].size(); it++) {
      rotB[io][0][0][it] /= (double)cnt; 
    }
  }

  /*
  // variance of rotB
  vector<double> varnc(lgndrs[2][ir].size(), 0);
  for (ir=11; ir<12; ir++) {
    double anorm = 0;
    for (it=0; it<(int)lgndrs[2][ir].size(); it++) {
      anorm += lgndrs[2][ir][it];
    }
    for (it=0; it<(int)lgndrs[2][ir].size(); it++) {
cout<<"test "<<it<<": "<<lgndrs[2][ir][it]/fabs(anorm) <<" "<<fabs(anorm)<<"  "<<rotB[0][0][it]<<endl;
      varnc[it] += pow(lgndrs[2][ir][it]/fabs(anorm) - rotB[0][0][it], 2);
    }
  }
  for (it=0; it<(int)rotB[0][0].size(); it++) varnc[it] /= (double)cnt;
exit(0);
*/

  //// Wiener Filter ////

  int FTsize = (int)(N/2+1);
  int noiseKthresh = (int)((N/2+1)*3/5);
  double* tSpace = (double*) fftw_malloc(sizeof(double)*N);
  fftw_complex* fSpace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*FTsize);
  vector< vector<complex<double> > > osnf(4);
  for (int j=0; j<4; j++) {
    osnf[j].resize(FTsize);
  }
  fftw_plan fftF;
  fftw_plan fftB;
  fftF = fftw_plan_dft_r2c_1d( N, tSpace, fSpace, FFTW_MEASURE);
  fftB = fftw_plan_dft_c2r_1d( N, fSpace, tSpace, FFTW_MEASURE);

  vector<double> linPart(FTsize-noiseKthresh, 0.0);
  TH1F* hist = new TH1F("fft", "fft", FTsize, 0, FTsize);
  TH1F* fitHist = new TH1F("lin", "lin", noiseKthresh, 0, noiseKthresh);
  TF1* linFit = new TF1("linFit", "[0] + [1]*x", 0, FTsize-noiseKthresh);


  for (io=0; io<rotB.size(); io++) {

    // Fill fftw time vector
    for (it=0; it<N; it++) {
      tSpace[it] = rotB[io][0][0][it];
    }

    plt.print1d(rotB[io][0][0], "rotBorig_" + to_string(io));
    // FT time vector
    fftw_execute(fftF);
    // Save values to go back
    for (int itt=0; itt<FTsize; itt++) {
      osnf[0][itt] = complex<double>(fSpace[itt][0], fSpace[itt][1]);
    }

    double bAvg=0;
    cnt=FTsize-noiseKthresh-1;
    for (int itt=FTsize-1; itt>noiseKthresh; itt--) {
      fitHist->SetBinContent(cnt+1, norm(osnf[0][itt]));
      bAvg += norm(osnf[0][itt]);
      cnt--;
    }
    bAvg /= FTsize-noiseKthresh;

    linFit->SetParameter(0, bAvg);
    linFit->SetParameter(1, 0);
    fitHist->Fit("linFit");
    intcpt = linFit->GetParameter(0);
    slope = linFit->GetParameter(1);

    for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, norm(osnf[0][ib]));
    //for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, pow(fSpace[ib][0],2)+pow(fSpace[ib][1],2));
    //hist->SetMaximum(10000);
    plt.print1d(hist,"rotBFT");
    fftw_execute(fftB);
    //for (it=0; it<(int)rotB[0][0].size(); it++) rotB[0][0][it] = tSpace[it];
    //plt.print1d(rotB[0][0], "rotBTest");


    // Subtracting the noise
    for (it=0; it<(int)osnf[0].size(); it++) {
      sbtr = linFit->GetMaximum(); //sbtr = slope*(it-noiseKthresh) + intcpt;
      scl = sqrt(sbtr/norm(osnf[0][it]));
      osnf[2][it] = complex<double>(osnf[0][it].real()*scl, osnf[0][it].imag()*scl);
      if (it+1 > noiseKthresh || sbtr > norm(osnf[0][it])) {
        osnf[1][it] = complex<double> (0,0);
        continue;
      }
      scl = sqrt((norm(osnf[0][it]) - sbtr)/norm(osnf[0][it]));
      osnf[1][it] = complex<double> (osnf[0][it].real()*scl, osnf[0][it].imag()*scl);
    }

    // Calculate Filter and fill fSpace
    for (it=0; it<(int)osnf[0].size(); it++) {
      osnf[3][it] = norm(osnf[1][it])/(norm(osnf[1][it]) + norm(osnf[2][it]));
    }

  /*
    for (int j=0; j<4; j++) {
      for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, norm(osnf[j][ib]));
      if (j==3) for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, sqrt(norm(osnf[j][ib])));
      plt.print1d(hist, "osnf_"+to_string(j));
    }
  */

    ///// Calculating non orthogonal rotation basis ////
    // Looping over the order of derivatives
    cnt=0;
    complex<double> der(0.0, 1.0);
    for (ir=0; ir<(int)rotB[io].size(); ir++) {
      for (it=0; it<FTsize; it++) {
        if (ir==0) {
          fSpace[it][0] = (osnf[0][it]*osnf[3][it]).real();
          fSpace[it][1] = (osnf[0][it]*osnf[3][it]).imag();
        }
        else {
          fSpace[it][0] = (osnf[0][it]*osnf[3][it]*pow(der*((double)it),ir)).real();
          fSpace[it][1] = (osnf[0][it]*osnf[3][it]*pow(der*((double)it),ir)).imag();
        }
      }

      fftw_execute(fftB);
      
      cnt = ir;
      for (id=0; id<npows; id++) {
        for (it=0; it<N; it++) {
          rotB[io][ir][id][it] = pow(tSpace[it],id+1)/((double)(1));
        }
        //gsOrth[cnt] = &(rotB[io][ir][id]);
        //cnt += rotB[io].size();
        plt.plot1d(rotB[io][ir][id], "rotBorig_L_" 
            + to_string(2*io) + "Der_" + to_string(ir) + "Pow_" + to_string(id))->Write();
      }
    }
  }

  // Orthogonalizing the rotation basis
  vector< vector<double>* > gsOrth(npows*rotB[0].size()*rotB.size());
  cnt = 0;
  for (io=0; io<NrotLg; io++) {
    for (ip=0; ip<npows; ip++) {
      for (int id=0; id<(int)rotB[io].size(); id++) {
        gsOrth[cnt] = &(rotB[io][id][ip]);
        cnt++;
      }
    }
  }

  tools::gramSchmidt(gsOrth, true);

  for (io=0; io<NrotLg; io++) {
    for (id=0; id<(int)rotB[io].size(); id++) {
      for (ip=0; ip<npows; ip++) {
        plt.plot1d(rotB[io][id][ip], "rotBasisOrth_L-" + to_string(2*io)
            + "_Der-" + to_string(id) 
            + "_Pow-" + to_string(ip))->Write();

        if (dataMC) {
          FILE* rbFile = 
            fopen(("output/MC/rotBasisOrth_L-" + to_string(2*io) 
                  + "_Der-" + to_string(id) 
                  + "_Pow-" + to_string(ip) + ".dat").c_str(),
                "wb");
          fwrite(&rotB[io][id][ip][0], sizeof(double), rotB[io][id][ip].size(), rbFile);
        }
      }
    }
  }

  cout<<"saved"<<endl;
  if (!dataMC) {
    io = 1;
    FILE* rbFile = 
        fopen(("../simulation/fitParameters/data/dataP" + to_string(io*2) + "_Size-" 
        + to_string(rotB[io][0][0].size()) + ".dat").c_str(),
        "wb");
    fwrite(&rotB[io][0][0][0], sizeof(double), rotB[io][0][0].size(), rbFile);
    fclose(rbFile); 
    plt.print1d(rotB[io][0][0], "testRotB");

    vector<double> tmplgndr2(lgndrs[2][10].size(), 0);
    int Navgs = 3;
    for (ir=9; ir<9+Navgs; ir++) {
      for (it=0; it<tmplgndr2.size(); it++) {
        tmplgndr2[ir] += lgndrs[2][ir][it]/((double)Navgs);
      }
    }

    rbFile =
        fopen(("../simulation/fitParameters/data/rawDataP" + to_string(io*2) + "_Size-"
        + to_string(rotB[io][0][0].size()) + ".dat").c_str(),
        "wb");
    fwrite(&lgndrs[2][10][0], sizeof(double), lgndrs[2][10].size(), rbFile);
    fclose(rbFile);

    rbFile =
        fopen(("../simulation/fitParameters/data/rawDataVarP" + to_string(io*2) + "_Size-"
        + to_string(varnc.size()) + ".dat").c_str(),
        "wb");
    fwrite(&varnc[0], sizeof(double), varnc.size(), rbFile);
    fclose(rbFile);
  }

  delete fitHist;
  delete hist; 

  fftw_destroy_plan(fftF);
  fftw_destroy_plan(fftB);
  fftw_free(tSpace);
  fftw_free(fSpace);


  ////  If data then retrieve simulated rotor basis  ////
  if (!dataMC && doMCbasis) {
    for (io=0; io<rotB.size(); io++) {
      for (id=0; id<rotB[id].size(); id++) {
        for (ip=0; ip<npows; ip++) {
          FILE* rbFile = 
            fopen(("output/MC/rotBasisOrth_L-" + to_string(2*io) 
                  + "_Der-" + to_string(id) 
                  + "_Pow-" + to_string(ip) + ".dat").c_str(),
                "rb");
          fread(&rotB[io][id][ip][0], sizeof(double), rotB[io][id][ip].size(), rbFile);
        }
      }
    }
  }

  cout<<"project"<<endl;
  /*
  for (io=0; io<rotB.size(); io++) {
    for (id=0; id<rotB[id].size(); id++) {
      for (ip=0; ip<npows; ip++) {
        for (int it=rotB[io][id][ip].size()-1; it>0; it--) {
          rotB[io][id][ip][it] = rotB[io][id][ip][it-1];
        }
        rotB[io][id][ip][0] = 0;
      }
    }
  }
  */

  ////  Project Onto the Rotational Basis  ////
  int ilg;
  int NlgProj = 5;
  double div=1;
  vector< vector< vector< vector< vector<double> > > > > rotBC(rotB.size());
  // Looping over rotB legendres
  for (io=0; io<NrotLg; io++) {
    rotBC[io].resize(NlgProj);
    // Looping over plot legendres
    for (ilg=0; ilg<NlgProj; ilg++) {
      rotBC[io][ilg].resize(Nders);
      // Looping over rotational basis first
      for (id=0; id<Nders; id++) {
        rotBC[io][ilg][id].resize(npows);
        for (ip=0; ip<npows; ip++) { 
          rotBC[io][ilg][id][ip].resize(lgndrs[io*2].size());
          for (ir=0; ir<(int)lgndrs[0].size(); ir++) {
            rotBC[io][ilg][id][ip][ir] = tools::vDotP(rotB[io][id][ip], lgndrs[ilg*2][ir])/div;
          }
	}
      }
    }
  }


  /*
  cout<<"open cos2"<<endl;
  ///// Project onto CosSq  /////
  int Ncspows = 6;
  gsOrth.resize(Ncspows);
  vector< vector<double> > cosP(Ncspows);
  vector< vector< vector<double> > > cosProj(Ncspows);
  for (uint i=0; i<Ncspows; i++) {
    cosP[i].resize(44);
    FILE* cosFile = fopen(("../simulation/rotation/output/expValCos" 
              + to_string((i+1)*2) + "_355.726000-360.026000.dat").c_str(), "rb");
    fread(&cosP[i][0], sizeof(double), cosP[i].size(), cosFile);
    fclose(cosFile);

    cosP[i].erase(cosP[i].begin()+41, cosP[i].end());
    cosP[i].erase(cosP[i].begin(), cosP[i].begin()+1);
    cout<<"IMP SIZE: "<<cosP.size();
    
    sum = 0;
    for (it=0; it<cosP[i].size(); it++) {
      sum += cosP[i][it];
    }
    sum /= (double)cosP[i].size();
    for (it=0; it<cosP[i].size(); it++) {
      cosP[i][it] -= sum;
    }
    cout<<"avg to subtract "<<i<<"   "<<sum<<endl;

    gsOrth[i] = &cosP[i];
    plt.print1d(cosP[i], "cos" + to_string((i+1)*2) + "BOrig");
  }

  vector<double> temp(cosP[0].size());
  for (uint i=0; i<cosP[0].size(); i++) temp[i] = cosP[0][i] - cosP[1][i]; 
  plt.print1d(temp, "cosDiff");
  tools::gramSchmidt(gsOrth, true);
   
  vector< vector<double> > cosPow(NlgProj);
  for (uint i=0; i<NlgProj; i++) {
    cosPow[i].resize(Ncspows, 0);
  }

  // Saving the coefficients
  if (!dataMC) {
    for (uint i=0; i<Ncspows; i++) {

      plt.print1d(cosP[i], "cos" + to_string((i+1)*2) + "BOrth");
      // Project onto cos2
      cosProj[i].resize(NlgProj);
      for (ilg=0; ilg<NlgProj; ilg++) {
        cosProj[i][ilg].resize(lgndrs[0].size());
        FILE* cosPfile = fopen(
            ("/reg/neh/home/khegazy/analysis/tumblingN2O/movie/data/UEDcos" 
            + to_string((i + 1)*2) + "CoefLeg" 
            + to_string(ilg) + "rad"
            + to_string(cosProj[i][ilg].size()) + ".dat").c_str(), "wb");
        for (ir=0; ir<(int)lgndrs[0].size(); ir++) {
          cosProj[i][ilg][ir] = tools::vDotP(cosP[i], lgndrs[ilg*2][ir]);
          cosPow[ilg][i] += pow(cosProj[i][ilg][ir], 1);
        }
        fwrite(&cosProj[i][ilg][0], sizeof(double), cosProj[i][ilg].size(), cosPfile);
        fclose(cosPfile);
      }
    }
  }

  // cos basis
  for (ilg=0; ilg<NlgProj; ilg++) {
    double sum = 0;
    for (io=0; io<Ncspows; io++) {
      sum += pow(cosPow[ilg][io], 2);
    }
    for (io=0; io<Ncspows; io++) {
      cosPow[ilg][io] = pow(cosPow[ilg][io], 2)/sum;
    }
  }


  h = (TH1F*)plt.print1d(cosPow[0], "leg_"+to_string(0));
  plt.MyC->cd();
  h->Draw();
  for (ilg=0; ilg<NlgProj; ilg++) {
    h = (TH1F*)plt.print1d(cosPow[ilg], "leg_"+to_string(ilg));
    h->SetLineColor(ilg*2);
    //h->Draw("same");
  }
  plt.MyC->Print("cospower.png");
  */
    


  cout<<"plotting"<<endl;
  vector< vector<double> > rotBmag(lgndrs[2].size());
  vector< vector<double> > rotBpow(NlgProj);
  for (ilg=0; ilg<NlgProj; ilg++) {
    rotBpow[ilg].resize(NrotLg*rotB[0].size()*npows, 0);
    for (ir=0; ir<(int)rotBmag.size(); ir++) {
      rotBmag[ir].clear(); //resize(rotBC.size()*rotBC[io].size()*npows, 0);
      for (ip=0; ip<npows; ip++) {
        for (io=0; io<(int)rotBC.size(); io++){
          for (id=0; id<Nders; id++) {
            int ind = ip + id*npows + io*rotBC[io][ilg].size()*npows;
            if (rotBpow[ilg].size() >= ind) {
              rotBpow[ilg].at(ind) = 0;
            }
            rotBpow[ilg][ind] += pow(rotBC[io][ilg][id][ip][ir], 2);
            rotBmag[ir].push_back(rotBC[io][ilg][id][ip][ir]);
          }
        }
      }
    }
    plt.plotRC(rotBmag, "RBcoeffs_" + to_string(ilg))->Write();
  }



  /// power spectrum ///
  // original basis
  for (ilg=0; ilg<NlgProj; ilg++) {
    double sum = 0;
    for (io=0; io<rotBpow[ilg].size(); io++) {
      sum += rotBpow[ilg][io];
    }
    for (io=0; io<rotBpow[ilg].size(); io++) {
      rotBpow[ilg][io] /= sum;
    }
  }


  /*
  TH1F* h = (TH1F*)plt.print1d(rotBpow[2], "leg_"+to_string(0));
  plt.MyC->cd();
  h->Draw();
  for (ilg=1; ilg<NlgProj; ilg++) {
    cout<<"ilg: "<<ilg<<endl;
    h = (TH1F*)plt.plot1d(rotBpow[ilg], "leg_"+to_string(ilg));
    cout<<"ksdf"<<endl;
    h->SetLineColor(ilg*2);
    h->Draw("same");
    cout<<"end"<<endl;
  }
  plt.MyC->Print("power.png");
  cout<<"finished plotting"<<endl;
  */

  cout<<"STARTING FFTBASIS"<<endl;
  /////  Fourier Basis  /////
  int NfftPows = 15;
  int NfftDers = 2;  // 1 => no derivatives
  int NfftBins = 40;
  // rotB[lgndr][der or not][pow][time]
  vector< vector< vector<double> > > fftBasis(NfftPows);
  vector<double> allBasis(NfftPows*NfftDers*NfftBins);
  FILE* cosFile = fopen("../rotorBasis/uedfftRotorbasis.dat", "rb");
  //FILE* cosFile = fopen("../rotorBasis/fftRotorbasis.dat", "rb");
  fread(&allBasis[0], sizeof(double), allBasis.size(), cosFile);
  fclose(cosFile);

  for (ip=0; ip<NfftPows; ip++) {
    fftBasis[ip].resize(NfftDers);
    for (id=0; id<NfftDers; id++) {
      fftBasis[ip][id].resize(NfftBins);
      for (it=0; it<NfftBins; it++) {
        fftBasis[ip][id][it] = allBasis[ip*NfftDers*NfftBins + id*NfftBins + it];
      }
    }
  }

  ///  Project onto fourier basis  ///
  vector< vector< vector< vector<double> > > > fftProj(NfftPows);
  vector< vector<double> > fftPow(NlgProj);
  for (int i=0; i<NlgProj; i++) {
    fftPow[i].resize(NfftPows*NfftDers, 0);
  }
  for (ip=0; ip<NfftPows; ip++) {
    fftProj[ip].resize(NfftDers);
    for (id=0; id<NfftDers; id++) {
      plt.print1d(fftBasis[ip][id], "fft" + to_string(ip) + "_" + to_string(id) + "BOrth");
      fftProj[ip][id].resize(NlgProj);
      for (ilg=0; ilg<NlgProj; ilg++) {
        fftProj[ip][id][ilg].resize(lgndrs[0].size());
//        FILE* fftPfile = fopen(
//            ("/reg/neh/home/khegazy/analysis/tumblingN2O/movie/data/UEDfft_" 
//            + to_string(ip) + "-" + to_string(id) + "_CoefLeg" 
//            + to_string(ilg*2) + "rad"
//            + to_string(fftProj[ip][id][ilg].size()) + ".dat").c_str(), "wb");
        for (ir=0; ir<(int)lgndrs[0].size(); ir++) {
          fftProj[ip][id][ilg][ir] = tools::vDotP(fftBasis[ip][id], lgndrs[ilg*2][ir]);
          fftPow[ilg][ip*NfftDers + id] += pow(fftProj[ip][id][ilg][ir], 2);
        }
        if (ilg == 2) {
          cout<<"pow: "<<ip<<"  "<<id<<"  "<<fftPow[ilg][ip*NfftDers + id]<<endl;
        }
//        fwrite(&fftProj[ip][id][ilg][0], sizeof(double), fftProj[ip][id][ilg].size(), fftPfile);
//        fclose(fftPfile);
      }
    }
  }

  if (!dataMC) {
    for (ip=0; ip<NfftPows; ip++) {
      for (id=0; id<NfftDers; id++) {
        for (ilg=0; ilg<NlgProj; ilg++) {
          FILE* fftPfile = fopen(
              ("/reg/neh/home/khegazy/analysis/tumblingN2O/movie/data/UEDfft_" 
              + to_string(ip) + "-" + to_string(id) + "_CoefLeg" 
              + to_string(ilg*2) + "rad"
              + to_string(fftProj[ip][id][ilg].size()) + ".dat").c_str(), "wb");
          fwrite(&fftProj[ip][id][ilg][0], sizeof(double), fftProj[ip][id][ilg].size(), fftPfile);
          fclose(fftPfile);
        }
      }
    }
  }


  for (ilg=0; ilg<NlgProj; ilg++) {
    plt.print1d(fftPow[ilg], "fftPowerSpect" + to_string(ilg*2));
  }

  cout<<"ENDING FFT BASIS"<<endl;




  //////////////////////////////////////////////////////////////////
  /////  Build Diffraction Patterns from Legendre Projections  /////
  //////////////////////////////////////////////////////////////////
  
  double Npix = 501; // Need odd number of pixels for autocorrelation
  int centdP = Npix/2;
  int ipix = (int)Npix;
  if (ipix%2 == 0) cerr<<endl<<endl<<"WARNING: Npix is even, ifft for autocoralation will have a phase!!!"<<endl<<endl;
  int Nrad = lgndrs[2].size();
  int irad;
  double rconv = Nrad/(Npix/2);
  double ang, rad;
  vector< vector< vector<double> > > diffP(N);
  vector< vector< vector<double> > > diffPpol(N);
  vector< vector<double> > diffPVar(ipix);
  vector< vector<double> > diffPVarCos2(ipix);
  vector< vector<double> > diffPpolVar(Nrad);
  for (it=0; it<N; it++) {
    diffP[it].resize(ipix);
    diffPpol[it].resize(Nrad);
    for (ir=0; ir<(int)diffP[it].size(); ir++) diffP[it][ir].resize(ipix, 0.0);
    for (ir=0; ir<Nrad; ir++) diffPpol[it][ir].resize(ipix, 0.0);
  }
  for (ir=0; ir<ipix; ir++) {
    diffPVar[ir].resize(ipix, 0.0);
    diffPVarCos2[ir].resize(ipix, 0.0);
  }
  for (ir=0; ir<Nrad; ir++) diffPpolVar[ir].resize(ipix, 0.0);

  double NdivBin = 4; // must be 1, 4, 16, ...
  double binArea, binHstep, binStep;
  bool hitBorder;

  uint stLg = 2;
  uint fnLg = 11;
  if (bendSearch) {
    stLg = 4;
    fnLg = 11;
  }

  std::vector<double> atmDiff(ipix/2 + 1);
  save::importDat<double>(atmDiff,
      "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/newRefAnalysis/diffractionPattern/output/nitrousOxide_atmDiffractionPatternLineOut_Qmax-"
      + to_string(maxQ)
      + "_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins["
      + to_string(251) + "].dat");
  cout<<"/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/newRefAnalysis/diffractionPattern/output/nitrousOxide_atmDiffractionPatternLineOut_Qmax-" + to_string(maxQ) + "_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins[" + to_string(251) + "].dat"<<endl;
  std::vector<double> molDiff(ipix/2 + 1);
  save::importDat<double>(molDiff,
      "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/diffractionPatterns/newRefAnalysis/diffractionPattern/output/nitrousOxide_molDiffractionPatternLineOut_Qmax-"
      + to_string(maxQ)
      + "_Ieb-5.000000_scrnD-4.000000_elE-3700000.000000_Bins["
      + to_string(251) + "].dat");


  // Calculate error from including 
  //    X numbers of rotational basis
  vector<double> rotBin(NrotLg*rotB[0].size()*npows, 0);
  vector<double> rotBinChiSq(rotBin.size(), 0);
  int chiSqTime = 21;

  // Normalize the diffraction pattern
  double lgNorm = 0;
  for (it=0; it<N; it++) {
    lgNorm += fabs(lgndrs[2][((int)Nrad/15)][it]);
  }
  cout<<"start diffp"<<endl;
  for (it=0; it<N; it++) {
//if (it!=26) continue;
    // Creating diffraction pattern in cartesian coordinates
    binArea = 1/NdivBin;
    binStep = 1/NdivBin;
    binHstep = binStep/2;
    for (ir=0; ir<=(int)(Npix/2); ir++) {
      for (ic=0; ic<=(int)(Npix/2); ic++) {
	diffP[it][ir][ic] = 0;
        hitBorder = false;
        if (it == chiSqTime) {
          for (int i=0; i<rotBin.size(); i++) {
            rotBin[i] = 0;
          }
        }

	// Looping through sub bins
	for (irr=0; irr<(int)NdivBin/2; irr++) {
	  for (icc=0; icc<(int)NdivBin/2; icc++) {

 	    rad = sqrt(pow((ir+binHstep+binStep*irr)-Npix/2, 2) + pow((ic+binHstep+binStep*icc)-Npix/2, 2));
            irad = (int)(rad*rconv)-1 < 0 ? 0 : (int)(rad*rconv)-1;
	    ang = atan(((ir+binHstep+binStep*irr)-Npix/2)/((ic+binHstep+binStep*icc)-Npix/2));  //FIX ANGLE
	    if (rad <= Npix/2) {
              for (io=stLg; io<fnLg; io+=2) { // Ignoring legendre 0
	        //diffP[it][ir][ic] += rotBC[id][ip][io][irad]*rotB[id][ip][it]*cos(ang);  // Consider area of pixel???
	        //diffP[it][ir][ic] += binArea*rotBC[id][ip][io][irad]*rotB[id][ip][it]*boost::math::legendre_p(io, cos(ang));
	        diffP[it][ir][ic] += binArea*lgndrs[io][irad][it]
                                     *boost::math::legendre_p(io, cos(ang))/lgNorm;

                // Calculate error from including 
                //    X numbers of rotational basis
                if (it == chiSqTime) {
                  if (io/2 >= NlgProj) {
                    continue;
                  }

                  int ind;
                  for (int i=0; i<NrotLg; i++) {
                    for (id=0; id<rotBC[i][io/2].size(); id++) {
                      for (ip=0; ip<npows; ip++) {
                        ind = ip + id*npows + i*rotBC[i][io/2].size()*npows;
                        //cout<<ir<<"/"<<ic<<" "<<ind<<" :  "<<binArea<<" "<<diffP[it][ir][ic]<<"/"<<rotBin[ind]<<" "<<rotBC[i][io/2][id][ip][irad]<<" "<<rotB[i][id][ip][it]<<" "<<boost::math::legendre_p(io, cos(ang))<<endl;
                        rotBin[ind] += binArea
                                    *rotBC[i][io/2][id][ip][irad]*rotB[i][id][ip][it]
                                    *boost::math::legendre_p(io, cos(ang))/lgNorm;
                      }
                    }
                  }
                }


	      }

	    }
	    else {
              diffP[it][ir][ic] = 0;
              hitBorder = true;
              break;
            }
	  }

          if (hitBorder) {
            break;
          }
	}

//cout<<ir<<"/"<<ic<<"  "<<ipix-ir-1<<"/"<<ipix-ic-1<<"  "<<ipix-ir-1<<"/"<<ic<<"  "<<ir<<"/"<<ipix-ic-1<<endl;
	diffP[it][ipix-ir-1][ipix-ic-1] = diffP[it][ir][ic];	// Diffraction pattern is symmetric
	diffP[it][ipix-ir-1][ic] = diffP[it][ir][ic];	// Diffraction pattern is symmetric
	diffP[it][ir][ipix-ic-1] = diffP[it][ir][ic];	// Diffraction pattern is symmetric

        // Calculate error from including 
        //    X numbers of rotational basis
        if ((it == chiSqTime) && !hitBorder) {

          bool verb = false;
            if ((ir==(int)(Npix/2)) && (ic==(int)(Npix/3))) {
              verb = true;
              cout<<"target: "<<diffP[it][ir][ic]<<endl;
            }
          for (uint i=0; i<rotBin.size(); i++) {
            if (i != 0) {
              rotBin[i] += rotBin[i-1];
            }
            if (verb) cout<<" new RB: "<<rotBin[i]<<endl;

            rotBinChiSq[i] += pow((diffP[it][ir][ic] - rotBin[i])/diffP[it][ir][ic], 2);

            /*
            if (sqrt(pow(ir-Npix/2, 2) + pow(ic-Npix/2, 2))/(Npix/2) <= 1) {
              int irad = (int)(*sqrt(pow(ir-Npix/2, 2) + pow(ic-Npix/2, 2))/(Npix/2));
              rotBinRadChiSq[i][irad] += pow((diffP[it][ir][ic] - rotBin[i])/diffP[it][ir][ic], 2);
            }
            */
            if (verb) cout<<" new CS: "<<rotBinChiSq[i]<<endl;
          }
        }

      }
    }

/*
    // Creating diffraction pattern in polar coordinates
    for (ir=0; ir<Nrad; ir++) {
      binArea = (pow(ir+1,2)-pow(ir,2))*PI/(sqrt(NdivBin)*Npix); 
      for (ic=0; ic<=(int)(Npix/4); ic++) {
	// Looping through sub bins
	for (icc=0; icc<(int)NdivBin/2; icc++) {
	  for (id=0; id<1; id++) { // only non div
	    for (ip=0; ip<1; ip++) { // only 0th normally ip<npows
	      for (io=0; io<(int)rotB.size(); io++) {	 
                for (ilg=0; ilg<(int)rotB[io].size(); ilg++) {
	          diffPpol[it][ir][ic] += binArea*rotBC[io][ilg][id][ip][ir]
                      *rotB[io][id][ip][it]
                      *boost::math::legendre_p(ilg*2, cos((ic+binHstep+binStep*icc)*2.0*PI/Npix));  	
                }
	      }
	    }
	  }
	}
	diffPpol[it][ir][ic+(int)(Npix/4)] = diffPpol[it][ir][ic];
	diffPpol[it][ir][ic+(int)(2*Npix/4)] = diffPpol[it][ir][ic];
	diffPpol[it][ir][ic+(int)(3*Npix/4)] = diffPpol[it][ir][ic];
      }
    }
*/
  }

  for (uint i=0; i<rotBinChiSq.size(); i++) {
    rotBinChiSq[i] = sqrt(rotBinChiSq[i])/pow(Npix/2 + 1, 2);
  }

  plt.print1d(rotBinChiSq, "rotBinChiSq");
              
  
  string dpName, dpPName, dpSName;
  vector<PLOToptions> opts(4);
  vector<PLOToptions> optspol(2);
  vector<string> dfpvars(4);
  vector<string> dfppvars(2);
  optspol[0]=maximum;		dfppvars[0] = "1500";
  optspol[1]=minimum;		dfppvars[1] = "-1500";
  opts[0] = xLabel;		dfpvars[0] = "Q #left[ #AA^{-1} #right]";
  opts[1] = yLabel;		dfpvars[1] = "Q #left[ #AA^{-1} #right]";
  opts[2] = xSpan;		dfpvars[2] = to_string(-maxQ)+","+to_string(maxQ);
  opts[3] = ySpan;		dfpvars[3] = dfpvars[2];
  //opts[4]=maximum;		/*dfpvars[4] = "1500";*/		dfpvars[4] = to_string(5*0.2*(1 + 0.8));	
  //opts[5]=minimum;		/*dfpvars[5] = "-1500";*/		dfpvars[5] = to_string(-5*0.2*(1 + 0.8));;
  for (it=0; it<N; it++) {
    if (bendSearch && (it != bendTimeItr)) {
      continue;
    }
    dpName = "diffP";
    dpSName = "diffSMS";
    dpPName = "diffPpol";
    if (it<10) {
      dpName += "0";
      dpPName += "0";
      dpSName += "0";
    }
    dpName += to_string(it);
    dpPName += to_string(it);
    dpSName += to_string(it);
    plt.plotRC(diffP[it], dpName, opts, dfpvars)->Write();
    //plt.plotRC(diffPpol[it], dpPName, optspol, dfppvars)->Write();
    std::vector< std::vector<double> > sms(diffP[it].size());
    for (int ir=0; ir<sms.size(); ir++) {
      sms[ir].resize(diffP[it][ir].size());
      for (int ic=0; ic<sms[ir].size(); ic++) {
        float hyp = std::sqrt(std::pow(ir-centdP, 2) + std::pow(ic-centdP, 2));
        if ((uint)hyp < atmDiff.size()) {
          sms[ir][ic] = diffP[it][ir][ic]/atmDiff[(int)hyp];
          //sms[ir][ic] = diffP[it][ir][ic]*((maxQ*hyp)/(atmDiff.size()))/atmDiff[(int)hyp];
        }
        else {
          sms[ir][ic] = 0;
        }
      }
    }
    plt.plotRC(sms, dpSName, opts, dfpvars)->Write();

  }
  cout<<"PASSED DIFFP"<<endl;

  // Line out of diffP at different angles
  if (false && !dataMC) {
    FILE* dpFile;
    int bLength = Npix/4;
    vector<double> dpout(bLength);
    vector<double> angs = {0, PI/6., PI/3., PI/2.};
    vector<int> times = {21, 26};
    for (it=0; it<times.size(); it++) {
      for (uint ia=0; ia<angs.size(); ia++) {
        for (int ird=Npix/16; ird<bLength; ird++) {
          dpout[ird] = 0;
          irad = (int)(ird*rconv)-1 < 0 ? 0 : (int)(ird*rconv)-1;
          for (io=2; io<Nlg; io+=2) {   // Consider only even Legendres (skip 0)
            dpout[ird] += lgndrs[io][irad][times[it]]*boost::math::legendre_p(io, cos(angs[ia]));
          }
        }
        dpFile = fopen(("output/data/diffPLineOut_T-"
          	+ to_string(times[it]) + "Ang-" + to_string(angs[ia]) 
          	+ "_Sr-" + to_string(maxQ*((double)Npix)/16.0) 
          	+ "_Er-" + to_string(maxQ*((double)Npix)/4.0) + ".dat").c_str(), "wb");
        fwrite(&dpout[0], sizeof(double), dpout.size(), dpFile);
        fclose(dpFile);
      }
    }
  } 


  /*
  // Plotting regions of variation
  for (id=0; id<1; id++) { // only non div
    for (ip=0; ip<1; ip++) { // only 0th normally ip<npows
      //for (int irr=0; irr<pix; irr++) diffPVar[irr].resize(pix,0.0); //Use only if one wants each rotB seperately
      //for (int irr=0; irr<Nrad; irr++) diffPpolVar[irr].resize(pix,0.0); //Use only if one wants each rotB seperately
      for (io=0; io<1(int)rotBC.size(); io++) {	// Only looking at the first vector
        for (ilg=1; ilg<(int)rotBC[io].size(); ilg++) {

          // Cartesian coordinate pattern
          for (ir=0; ir<ipix; ir++) {
            for (ic=0; ic<ipix; ic++) {
              rad = sqrt(pow(ir-Npix/2, 2) + pow(ic-Npix/2, 2));
              irad = (int)(rad*rconv)-1 < 0 ? 0 : (int)(rad*rconv)-1;
              ang = atan(((double)(ir-Npix/2))/((double)(ic-Npix/2)));  //FIX ANGLE
              if (rad <= Npix/2) {
                diffPVar[ir][ic] += rotBC[io][ilg][id][ip][irad]*boost::math::legendre_p(2*ilg, cos(ang));  // Consider area of pixel???  
                if ((id == 0) && (ip == 0) && (io == 0)) {
                  diffPVarCos2[ir][ic] += cosProj[0][ilg][irad]*boost::math::legendre_p(2*ilg, cos(ang));  // Consider area of pixel???  
                }
              }
              else {
                diffPVar[ir][ic] += 0;
                diffPVarCos2[ir][ic] += 0;
              }
            }
          }

          // Polar coordinate pattern
          for (ir=0; ir<Nrad; ir++) {
            for (ic=0; ic<ipix; ic++) {
              diffPpolVar[ir][ic] += rotBC[io][ilg][id][ip][ir]*boost::math::legendre_p(2*ilg, cos((double)ic*2.0*PI/ipix));  // Consider area of pixel???  
            }
          }
        }
      }
      //plt.printRC(diffPVar, "diffP"+to_string(io)+"Var"+to_string(id)+to_string(ip), opts, valsdv); 
    }
  }


  plt.plotRC(diffPVar, "diffPVar")->Write();
  plt.printRC(diffPVarCos2, "diffPVarCos2")->Write();
  plt.plotRC(diffPpolVar, "diffPpolVar")->Write();
  */



  //////////////////////////////////////////////////////////
  /////  Reverse FFT to get molecular autocorrelation  /////
  //////////////////////////////////////////////////////////

  int fftPix = 15*Npix + 1;	// Must be an odd number for ifft
  double atcRat = 0.12*201/Npix;  //(int)(fftPix*Npix/((double)fftPix + Npix)); 
  int maxAtcItr = atcRat*fftPix; //maxAtcBinVal*fftPix;
  opts.resize(4);
  vector<string> acvars(4);
  opts[0] = xLabel;		acvars[0] = "#left[ #AA #right]";
  opts[1] = yLabel;		acvars[1] = "#left[ #AA #right]";
  opts[2] = xSpan;		acvars[2] = to_string(-atcRat*2*PI*(double)((int)(Npix/2))/(2*maxQ))+","+to_string(atcRat*2*PI*(double)((int)(Npix/2))/(2*maxQ));
  //opts[2] = xSpan;		acvars[2] = to_string(-((Npix))/(2*410*0.0358))+","+to_string(((Npix))/(2*410*0.0358));
  opts[3] = ySpan;		acvars[3] = acvars[2];

  vector<string> acDvars(6);
  acDvars[0] = acvars[0];
  acDvars[1] = acvars[1];
  acDvars[2] = acvars[2];
  acDvars[3] = acvars[3];

  vector< vector< vector<double> > > autCRel(N);
  vector< vector< vector<double> > > autCImg(N);
  for (it=0; it<N; it++) {
    autCRel[it].resize(maxAtcItr);
    autCImg[it].resize(maxAtcItr);
    for (ir=0; ir<(int)autCRel[it].size(); ir++) {
      autCRel[it][ir].resize(maxAtcItr, 0.0);
      autCImg[it][ir].resize(maxAtcItr, 0.0);
    }
  }


  fftw_complex* autCfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pow(fftPix,2));
  fftw_complex* difPfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pow(fftPix,2));
  fftw_plan fftAC;
  fftAC = fftw_plan_dft_2d(fftPix, fftPix, difPfft, autCfft, FFTW_BACKWARD, FFTW_MEASURE);


  int iref = Npix/13;
  int itr, itc;
  double hyp;
  vector< vector<double> > tfft(fftPix);
  double decLength = (fftPix - Npix)/(2*1.1);
  double omega = PI/(2*decLength);
  string atcName;
  vector< vector<double> > atcRad;
  vector< vector<double> > autLgndr(N);
  int atcNHbins = 50;
  double maxAtcBinVal = 1.0/((double)Npix*2);
 

  /////  Searching for Bending (altering diffP)  /////
  //double DFPthresh = 3.5e-4;
  double DFPthresh = 7e-4;
  vector< vector<double> > bendOrigDFP(diffP[bendTimeItr].size());
  vector< vector<double> > DFPchiSqMap(diffP[bendTimeItr].size());
  vector< vector<double> > chiSqMap(maxAtcItr);
  if (bendSearch) {
    it = bendTimeItr;
    if (!tools::fileExists("references/diffDataPattern_Bins-"
                + to_string(diffP[it].size()) + "_time-"
                + to_string(it) + ".dat")) {
        cerr << "ERROR: Cannot find diffraction reference file "
	  	+ ("references/diffDataPattern_Bins-"
                + to_string(diffP[it].size()) + "_time-"
                + to_string(it) + ".dat") + "\n\n\n"; 
        exit(0);
    }
                      
    FILE* inpDFPfile = fopen(("references/diffDataPattern_Bins-"
                + to_string(diffP[it].size()) + "_time-"
                + to_string(it) + ".dat").c_str(), "rb");
    double* inpDFP = new double[diffP[it].size()*diffP[it].size()];
    fread(inpDFP, sizeof(double), diffP[it].size()*diffP[it].size(), inpDFPfile);

    int ind;
    double norm = 1.0;
    int cent = diffP[it].size()/2;
    uint matSize = pow(diffP[it].size(), 2);
    //Eigen::MatrixXd X(matSize, 2);
    //Eigen::MatrixXd W = Eigen::MatrixXd::Zero(matSize, matSize);
    //Eigen::VectorXd Y(matSize);
    Eigen::Vector2d Th;

    for (int irr=0; irr<diffP[it].size(); irr++) {
      for (int icc=0; icc<diffP[it][irr].size(); icc++) {
        ind = irr*diffP[it].size() + icc;
        //X(ind, 0) = 1;
        //X(ind, 1) = diffP[it][irr][icc];
        //Y(ind) = inpDFP[ind];
        rad = sqrt((pow(irr - cent, 2) + pow(icc - cent, 2)));
        if ((rad > iref) && (rad < cent)) {
          if (fabs(diffP[it][irr][icc]) > DFPthresh) {
            norm = 1.0/fabs(diffP[it][irr][icc]);
          }
          //W(ind, ind) = norm;
        }
        else {
          //W(ind, ind) = 0;
        }
      }
    }
    
cout<<"calculating theta"<<endl;
    //Th = ((X.transpose()*W*X).inverse())*(X.transpose()*W*Y);

    Th(0) = 0.000301894;
    Th(1) = 1.10089;
cout<<"THETAS: "<<Th(0)<<"   "<<Th(1)<<endl;
cout<<"filling containers"<<endl;
    double DFPchiSq = 0;
    double NdfpChiSqBins = 0;
    for (int irr=0; irr<diffP[bendTimeItr].size(); irr++) {
      bendOrigDFP[irr].resize(diffP[bendTimeItr][irr].size(), 0);
      DFPchiSqMap[irr].resize(diffP[bendTimeItr][irr].size(), 0);
      for (int icc=0; icc<diffP[bendTimeItr][irr].size(); icc++) {
	ind = irr*diffP[it].size() + icc;
        bendOrigDFP[irr][icc] = Th(1)*diffP[it][irr][icc] + Th(0);
        diffP[it][irr][icc] = inpDFP[ind] - bendOrigDFP[irr][icc];
        rad = sqrt((pow(irr - cent, 2) + pow(icc - cent, 2)));
        if ((rad > iref) && (rad < cent) && (rad < cent*8.5/14.) && (fabs(bendOrigDFP[irr][icc]) > DFPthresh)) {
 	  DFPchiSq += pow(diffP[it][irr][icc], 2)/fabs(bendOrigDFP[irr][icc]);
          NdfpChiSqBins++;
          DFPchiSqMap[irr][icc] = pow(diffP[it][irr][icc], 2)/fabs(bendOrigDFP[irr][icc]);
	}
      }
    }

    double stdev = 0;
    double mean = DFPchiSq/NdfpChiSqBins;
    for (int irr=0; irr<diffP[bendTimeItr].size(); irr++) {
      for (int icc=0; icc<diffP[bendTimeItr][irr].size(); icc++) {
	ind = irr*diffP[it].size() + icc;
        rad = sqrt((pow(irr - cent, 2) + pow(icc - cent, 2)));
        if ((rad > iref) && (rad < cent) && (rad < cent*8.5/14.) && (fabs(bendOrigDFP[irr][icc]) > DFPthresh)) {
          stdev += pow(DFPchiSqMap[irr][icc] - mean, 2);
        }
      }
    }

    DFPchiSq = 0;
    stdev = sqrt(stdev/NdfpChiSqBins);
    for (int irr=0; irr<diffP[bendTimeItr].size(); irr++) {
      for (int icc=0; icc<diffP[bendTimeItr][irr].size(); icc++) {
	ind = irr*diffP[it].size() + icc;
        rad = sqrt((pow(irr - cent, 2) + pow(icc - cent, 2)));
        if ((rad > iref) && (rad < cent) && (rad < cent*8.5/14.) && (fabs(bendOrigDFP[irr][icc]) > DFPthresh)) {
          if (DFPchiSqMap[irr][icc] < mean + 3.5*stdev) {
            DFPchiSq += DFPchiSqMap[irr][icc];
          }
          else {
            DFPchiSqMap[irr][icc] = 0;
          }
        }
      }
    }

cout<<"DFP CHI SQ: "<<DFPchiSq<<endl;
    plt.printRC(DFPchiSqMap, "DFPchiSq", xSpan, dfpvars[2])->Write();
  } 
  uint bendItr = 0;

  // End of bending segment


  for (it=0; it<N; it++) {
    if (bendSearch && (it != bendTimeItr)) {
      continue;
    }
    if (bendItr > 1) {
      break;
    }

    autLgndr[it].resize(atcNHbins);
    for (ir=0; ir<=fftPix/2; ir++) {
      for (ic=0; ic<=fftPix/2; ic++) {

        hyp = sqrt(pow(ir, 2) + pow(ic, 2));
        int ind = ir*fftPix + ic;
        if (hyp < Npix/2) {
	  difPfft[ind][0] = diffP[it][ir + centdP][ic + centdP];
        }
       	else if ((hyp < iref)) {
	  //difPfft[ind][0] = diffP[it][(centdP+(int)(iref*((double)ir)/hyp))%ipix][(centdP+(int)(iref*((double)ic)/hyp))%ipix]*pow(cos((PI/(2*iref))*(iref-hyp)), 1);
	  difPfft[ind][0] = diffP[it][(centdP+(int)(iref*((double)ir)/hyp))%ipix][(centdP+(int)(iref*((double)ic)/hyp))%ipix]*molDiff[ind]/molDiff[iref];
	}
	else {
	  difPfft[ind][0] = 0;
 	}

        if (hyp <= Npix/2) {
          difPfft[ind][0] = difPfft[ind][0]
            /atmDiff[(int)hyp]
            //*(maxQ*hyp/(atmDiff.size()))/atmDiff[(int)hyp]
            *std::exp(-0.5*std::pow(hyp/110, 2));
            //*std::exp(-0.5*std::pow(hyp/75, 2));

            //*std::pow(cos((hyp-Npix/4)*PI/(Npix/1.5)), 2);
            //*std::pow(sin(hyp*PI/(Npix/2)), 2);
          cout<<ir<<"  "<<ic<<"  "<<difPfft[ind][0]<<endl;
        }
        difPfft[ind][1] = 0.0;

	if (ic != 0) difPfft[ir*(fftPix)+fftPix-ic][0] = difPfft[ir*fftPix+ic][0];
	if (ir != 0) difPfft[(fftPix-ir)*(fftPix)+ic][0] = difPfft[ir*fftPix+ic][0];
        if (ic != 0 && ir != 0) difPfft[(fftPix-ir)*(fftPix)+fftPix-ic][0] = difPfft[ir*fftPix+ic][0]; 


	difPfft[ir*(fftPix)+fftPix-1-ic][1] = difPfft[(fftPix-1-ir)*(fftPix)+ic][1] = difPfft[(fftPix-1-ir)*(fftPix)+fftPix-1-ic][1] = 0;
      }
    }

    fftw_execute(fftAC);

    for (ir=0; ir<maxAtcItr; ir++) { 
      for (ic=0; ic<maxAtcItr; ic++) { 
        if (ic > maxAtcItr/2) itc = fftPix-1+(ic-maxAtcItr);
	else itc = ic;
	if (ir > maxAtcItr/2) itr = fftPix-1+(ir-maxAtcItr);
	else itr = ir;
	autCRel[it][((int)(maxAtcItr/2)+ir)%maxAtcItr][((int)(maxAtcItr/2)+ic)%maxAtcItr] = autCfft[itr*fftPix+itc][0]/(fftPix*fftPix);
	//autCImg[it][((int)(atcBins/2)+ir)%atcBins][((int)(atcBins/2)+ic)%atcBins] = autCfft[itr*fftPix+itc][1]/(fftPix*fftPix);
      }
    }

    // Plotting autocorrelation
    if (!bendSearch) {
      atcName = "autCRel";
      if (it < 10) atcName += "0";
      atcName += to_string(it);
      //plt.printRC(autCRel[it], atcName)->Write();
      plt.plotRC(autCRel[it], atcName, opts, acvars)->Write();
      plt.plotRC(imgProc::polarBinning(autCRel[it], 
                  autCRel[it].size()/2+0.5, autCRel[it].size()/2+0.5,
                  80, autCRel[it].size()/2-1, autCRel[it].size()/4),
                  atcName + "Polar")->Write();


      atcName = "autCImg";
      if (it < 10) atcName += "0";
      atcName += to_string(it);
      plt.plotRC(autCImg[it], atcName, opts, acvars)->Write();
    }


    // Compare data and MC autocorrelations with least squares 
    if (!dataMC && (it == bendTimeItr)) {
      FILE* outATCfile = fopen(("references/autCorr_Bins-" 
		+ to_string(autCRel[it].size()) + "_time-" 
		+ to_string(it) + ".dat").c_str(), "wb");
      for (uint ii=0; ii<autCRel[it].size(); ii++) {
        fwrite(&autCRel[it][ii][0], sizeof(double), autCRel[it][ii].size(), outATCfile);
      }
      fclose(outATCfile);
 
      FILE* outDFPfile = fopen(("references/diffDataPattern_Bins-" 
		+ to_string(diffP[it].size()) + "_time-" 
		+ to_string(it) + ".dat").c_str(), "wb");
      for (uint ii=0; ii<diffP[it].size(); ii++) {
        fwrite(&diffP[it][ii][0], sizeof(double), diffP[it][ii].size(), outDFPfile);
      }
      fclose(outDFPfile);
    }
    else if (bendSearch) {

      if (!bendItr) {
        for (int irr=0; irr<(int)autCRel[it].size(); irr++) {
	  chiSqMap[irr].resize(maxAtcItr, 0);
 	  for (int icc=0; icc<(int)autCRel[it][irr].size(); icc++) {
            chiSqMap[irr][icc] = pow(autCRel[it][irr][icc], 2);
	  }
        }
        for (int irr=0; irr<(int)diffP[it].size(); irr++) {
          for (int icc=0; icc<(int)diffP[it][irr].size(); icc++) {
            diffP[it][irr][icc] = bendOrigDFP[irr][icc];
          }
        }
        it--;
      }
      else {
	for (int irr=0; irr<(int)autCRel[it].size(); irr++) {
          for (int icc=0; icc<(int)autCRel[it][irr].size(); icc++) {
            chiSqMap[irr][icc] /= fabs(autCRel[it][irr][icc]);
          }
        }

      }
      bendItr++;


/*
      if (!tools::fileExists("references/autCorr_Bins-"
                + to_string(autCRel[it].size()) + "_time-"
                + to_string(it) + ".dat")) {
	cerr << "ERROR: Cannot find autocorrelation reference file " 
	 	<< "references/autCorr_Bins-"
                + to_string(autCRel[it].size()) + "_time-"
                + to_string(it) + ".dat!!!! \nPlese run this script with data first\n\n";
	//exit(0);
	continue;
      }
      
      string diffAtcName = "autChiSq";
      if (it < 10) diffAtcName += "0";
      diffAtcName += to_string(it);

      FILE* inpATCfile = fopen(("references/autCorr_Bins-"
                + to_string(autCRel[it].size()) + "_time-"
                + to_string(it) + ".dat").c_str(), "rb");
      double* inpATC = new double[autCRel[it].size()*autCRel[it].size()];
      fread(inpATC, sizeof(double), autCRel[it].size()*autCRel[it].size(), inpATCfile);

      // Fitting scaling and offset parameter
      int ind;
      double norm = 1.0;
      double NchiSqBins = 0;
      int cent = autCRel[it].size()/2;
      int iref = cent/10;
      vector< vector<double> > testThresh(autCRel[it].size());
      double thresh = 5.5e-6;
      uint matSize = pow(autCRel[it].size(), 2);
      Eigen::MatrixXd X(matSize, 2);
      Eigen::MatrixXd W = Eigen::MatrixXd::Zero(matSize, matSize);
      Eigen::VectorXd Y(matSize);
      Eigen::Vector2d Th;

      for (int irr=0; irr<autCRel[it].size(); irr++) {
        testThresh[irr].resize(autCRel[it][irr].size(),1);
 	for (int icc=0; icc<autCRel[it][irr].size(); icc++) {
	  ind = irr*autCRel[it].size() + icc;
	  X(ind, 0) = 1;
	  X(ind, 1) = autCRel[it][irr][icc];
	  Y(ind) = inpATC[ind];
	  rad = sqrt((pow(irr - cent, 2) + pow(icc - cent, 2)));
	  if ((rad > iref) && (rad < cent) && (fabs(Y(ind)) > thresh)) {
	    if (autCRel[it][irr][icc]) {
	      norm = 1.0/pow(autCRel[it][irr][icc], 2);
	    }
	    W(ind, ind) = norm;
	    NchiSqBins ++; 
	  }
	  else {
	    W(ind, ind) = 0;
	    testThresh[irr][icc] = 0;
	  }
  	}
      }
      
      Th = ((X.transpose()*W*X).inverse())*(X.transpose()*W*Y);

   
      // Calculating percent error
      vector< vector<double> > autChiSq(autCRel[it].size());
      double sum = 0;
      double prevVal;
      norm = 1;
      for (int irr=0; irr<autCRel[it].size(); irr++) {
  	autChiSq[irr].resize(autCRel[it].size(), 0);
 	for (int icc=0; icc<autCRel[it][irr].size(); icc++) {
	  rad = sqrt(pow(irr - cent, 2) + pow(icc - cent, 2));
	  if ((rad > iref) && (rad < cent)) { 
            if (testThresh[irr][icc]) {
	      norm = inpATC[irr*autCRel[it].size() + icc];
	      autChiSq[irr][icc] = pow((Th(1)*autCRel[it][irr][icc]  + Th(0)
			- inpATC[irr*autCRel[it].size() + icc])
			/norm, 2)/NchiSqBins;
	    //if (autChiSq[irr][icc] >= 0.001 ) {
	    //  autChiSq[irr][icc] = 0;
	    //}
            sum += autChiSq[irr][icc];
testThresh[irr][icc] = 0;
autCRel[it][irr][icc] /= norm;
	    }
	    else {
autCRel[it][irr][icc] = 0;
testThresh[irr][icc] = autCRel[it][irr][icc];
	    }
	  }
	  else {
	    autChiSq[irr][icc] = 0;
	  }
        }
      }

delete plt.printRC(testThresh, "testThresh");
plt.printRC(autCRel[it], "autCAlnNorm")->Write();
      //plt.printRC(autChiSq, diffAtcName, maximum, "0.00007")->Write(); 
      plt.printRC(autChiSq, diffAtcName)->Write(); 
      plt.printRC(imgProc::polarBinning(autChiSq, 
		autCRel[it].size()/2+0.5, autCRel[it].size()/2+0.5,
		80, autCRel[it].size()/2-1, autCRel[it].size()/4),
		diffAtcName + "Pol")->Write();

      outTxt << "Autocorrelation Chi Squared " << it << " : " << sum << endl;
      cout << "Autocorrelation Chi Squared " << it << " : " << sum << endl;
*/

 
    }
 

    //atcRad = imgProc::polarBinning(autCRel[it], Npix/2+1, Npix/2+1, 80, Npix/2 - 1, atcNHbins);
    //for (uint iar=0; iar<atcRad.size(); iar++) {
    //  autLgndr[it][iar] = imgProc::legendre1dCoeff(atcRad[iar], 2);
   // }
  }


  if (bendSearch) {
    it = bendTimeItr;
      if (!tools::fileExists("references/autCorr_Bins-"
                + to_string(autCRel[it].size()) + "_time-"
                + to_string(it) + ".dat")) {
        cerr << "ERROR: Cannot find autocorrelation reference file " 
                << "references/autCorr_Bins-"
                + to_string(autCRel[it].size()) + "_time-"
                + to_string(it) + ".dat!!!! \nPlese run this script with data first\n\n";
        exit(0);
       }
//delete plt.printRC(testThresh, "testThresh");
//plt.printRC(autCRel[it], "autCAlnNorm")->Write();
      //plt.printRC(autChiSq, diffAtcName, maximum, "0.00007")->Write(); 
    FILE* inpATCfile = fopen(("references/autCorr_Bins-" 
		+ to_string(autCRel[it].size()) + "_time-" 
		+ to_string(it) + ".dat").c_str(), "rb");
    double* inpATC = new double[autCRel[it].size()*autCRel[it].size()];
    fread(inpATC, sizeof(double), autCRel[it].size()*autCRel[it].size(), inpATCfile);
 

cout<<"plotting"<<endl;
    string diffAtcName = "autChiSq";
    if (it < 10) diffAtcName += "0";
    diffAtcName += to_string(it);

cout<<"summing"<<endl;
    //double thresh = 1.2e-6;
    double thresh = 3.5e-6;
    //double ATCthresh = 1.05e-6;

//TH1F* thist = new TH1F();
    int ind = 0;
    int cent = autCRel[it].size()/2;
    int iref = cent/5.;
    //int tmin = cent/2.;
    //int tmax = cent*11./15.;
    int tmin = cent/2.2;
    int tmax = cent*10./13.;
cout<<"min/max: "<<tmin<<"  "<<tmax<<endl;
    double NchiSqBins = 0;
    double totChiSq = 0;
    for (int irr=0; irr<autCRel[bendTimeItr].size(); irr++) {
      for (int icc=0; icc<autCRel[bendTimeItr][irr].size(); icc++) {
 	rad = sqrt(pow(irr - cent, 2) + pow(icc - cent, 2));
        ind = irr*autCRel[bendTimeItr].size() + icc;
        //chiSqMap[irr][icc] /= fabs(inpATC[ind]);
        if ((fabs(inpATC[ind]) > thresh) && (rad > iref) && (rad < cent)) {
        //if ((fabs(autCRel[it][irr][icc]) > ATCthresh) 
	//	&& (rad > iref) && (rad < cent)
	//	&& (((rad < tmin) || (rad > tmax)))
	//	&& (chiSqMap[irr][icc] < 1e-5)) {
//thist->Fill(chiSqMap[irr][icc]);
          totChiSq += chiSqMap[irr][icc];
	  NchiSqBins++;
	}
	else {
	  chiSqMap[irr][icc] = 0;
	}
      }
    }
   
    if (!NchiSqBins) {
      cerr << "\n\nERROR: There were no bins added in the Chi Sq!!!\n\n";
    }
    totChiSq /= NchiSqBins;    

cout<<"plotting"<<endl;
    plt.printRC(chiSqMap, diffAtcName)->Write(); 
    plt.printRC(imgProc::polarBinning(chiSqMap, 
      	autCRel[bendTimeItr].size()/2+0.5, autCRel[bendTimeItr].size()/2+0.5,
      	80, autCRel[bendTimeItr].size()/2-1, autCRel[bendTimeItr].size()/4),
      	diffAtcName + "Pol")->Write();

cout<<"plotted"<<endl;
    outTxt << "Autocorrelation Chi Squared " << it << " : " << totChiSq << endl;
    cout << "Autocorrelation Chi Squared " << it << " : " << totChiSq << endl;
  }



/*
  vector<PLOToptions> vopts(2);
  vector<string> vvals(2);
  vopts[0] = maximum;
  vopts[1] = minimum;
  //vvals[0] = "5e-6";
  //vvals[1] = "-5e-6";
  vvals[0] = "0.15e-9";
  vvals[1] = "-0.15e-9";


  plt.printXY(autLgndr, "autLgndr2", vopts, vvals);
*/

/*
  //// Wiener Filter Time Array ////
  
  vector< vector<double> > timeD(Npix/2);
  for (uint j=0; j<timeD.size(); j++) {
    timeD[j].resize(N);
  }

  for (int ib=Npix/2; ib<Npix; ib++) {

    // Fill fftw time vector
    for (it=0; it<N; it++) {
      tSpace[it] = diffP[it][Npix/2][ib];
if (ib==121) cout<<diffP[it][Npix/2][ib]<<"   ";
    }
 cout<<endl<<endl;

    // FT time vector
    fftw_execute(fftF);
    // Save values to go back
    for (int itt=0; itt<FTsize; itt++) {
      osnf[0][itt] = complex<double>(fSpace[itt][0], fSpace[itt][1]);
if (ib==121) cout<<osnf[0][itt]<<"   ";
    }
cout<<endl<<endl;
 
    //for (int ibb=0; ibb<FTsize; ibb++) hist->SetBinContent(ibb+1, norm(osnf[0][ibb]));
    //for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, pow(fSpace[ib][0],2)+pow(fSpace[ib][1],2));
    //hist->SetMaximum(10000);
    //plt.print1d(hist,"rotBFT");
    //fftw_execute(fftB);
    //for (it=0; it<(int)rotB[0][0].size(); it++) rotB[0][0][it] = tSpace[it];
    //plt.print1d(rotB[0][0], "rotBTest");
 
 
    // Subtracting the noise
    sbtr = 0;
    for (int itt=FTsize-1; itt>noiseKthresh; itt--) {
      //if (norm(osnf[0][itt]) > sbtr) sbtr = norm(osnf[0][itt]);
      sbtr += norm(osnf[0][itt]);
    }
    sbtr /= (double)((FTsize - 1) - noiseKthresh); 

    for (it=0; it<(int)osnf[0].size(); it++) {
      scl = sqrt(sbtr/norm(osnf[0][it]));
      osnf[2][it] = complex<double>(osnf[0][it].real()*scl, osnf[0][it].imag()*scl);
      if (it+1 > noiseKthresh || sbtr > norm(osnf[0][it])) {
        osnf[1][it] = complex<double> (0,0);
        continue;
      }
      scl = sqrt((norm(osnf[0][it]) - sbtr)/norm(osnf[0][it]));
      osnf[1][it] = complex<double> (osnf[0][it].real()*scl, osnf[0][it].imag()*scl);
    }

    // Calculate Filter and fill fSpace
    for (it=0; it<(int)osnf[0].size(); it++) {
      osnf[3][it] = norm(osnf[1][it])/(norm(osnf[1][it]) + norm(osnf[2][it]));
if (ib==121) cout<<osnf[3][it]<<"   ";
    }
cout<<endl<<endl;

    for (it=0; it<FTsize; it++) {
      fSpace[it][0] = (osnf[0][it]*osnf[3][it]).real();
      fSpace[it][1] = (osnf[0][it]*osnf[3][it]).imag();
    }

    fftw_execute(fftB);
    
    for (it=0; it<N; it++) {
      timeD[ib-Npix/2][it] = tSpace[it];
if (ib==121) cout<<osnf[3][it]<<"   ";
}
cout<<endl<<endl;
    plt.print1d(timeD[ib-Npix/2], "qTimeDep"+to_string(ib));
  }
*/


  ////  FFT time domain for J state information  ////
  vector< vector< vector<double> > > cntrfgl(Nlg);
  for (int il=0; il<Nlg; il++) {
    cntrfgl[il].resize(Nrad);
    for (ir=0; ir<Nrad; ir++) {
      cntrfgl[il][ir].resize(N/2);
    }
  }

  fftw_complex* inpTimefft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fftw_complex* outTimefft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*N);
  fftw_plan fftTM;
  fftTM = fftw_plan_dft_1d(N, inpTimefft, outTimefft, FFTW_FORWARD, FFTW_MEASURE);

  opts.resize(5);
  vector<string> cfvars(5);
  opts[0] = xLabel;		cfvars[0] = "cm^{-1}";
  opts[1] = yLabel;		cfvars[1] = "Q #left[ #AA^{-1} #right]";
  opts[2] = xSpan;		cfvars[2] = "0,"+to_string((N/2+1)/(C_CMpS*N*50e-15)); //cfvars[2] = "0,"+to_string((N+1)*33.35641/(N*50));
  opts[3] = ySpan;		cfvars[3] = "0,"+to_string(maxQ);
  opts[4] = maximum;		cfvars[4] = "0.3";

  double tmNorm;
  for (int il=0; il<Nlg; il++) {
    for (ir=0; ir<Nrad; ir++) {
 
      // Filling fft array
      for (it=0; it<N; it++) {
        inpTimefft[it][0] = lgndrs[il][ir][it];
        inpTimefft[it][1] = 0.0;
      }

      fftw_execute(fftTM);
 
      tmNorm = 0;
      for (it=0; it<N/2; it++) tmNorm += sqrt(pow(outTimefft[int(N/2+it)][0],2) + pow(outTimefft[int(N/2+it)][1],2));
      for (it=0; it<N/2; it++) cntrfgl[il][ir][it] = sqrt(pow(outTimefft[it][0],2) + pow(outTimefft[it][1],2))/tmNorm;
    }

    plt.printRC(cntrfgl[il], "centrifugal_"+to_string(il), opts, cfvars);
    //if (it == 21) {
    //  for (ir=0; ir<Nrad+addBin; ir++) plt.print1d(autCPow[it][ir], "autCAng_"+to_string(ir));
    //}
  }

/*
  fftw_complex* autCfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pix*pix);
  fftw_complex* difPfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*pix*pix);
  fftw_plan fftAC;
  fftAC = fftw_plan_dft_2d(pix, pix, difPfft, autCfft, FFTW_BACKWARD, FFTW_MEASURE);

vector< vector<double> > test(pix);
for (int j=0; j<pix; j++) test[j].resize(pix, 0.0);
  int ac, ar, centr, centc;
  centr = (int)(pix/2);
  centc = (int)(pix/2);
  for (ar=0; ar<pix; ar++) {
    for (ac=0; ac<pix; ac++) {
      difPfft[ar*pix+ac][0] = diffP[24][(centr+ar)%pix][(centc+ac)%pix];
      difPfft[ar*pix+ac][1] = 0.0;
      test[ar][ac] = diffP[24][(centr+ar)%pix][(centc+ac)%pix];
    }
  }

cout<<"start fft"<<endl;
  fftw_execute(fftAC);
cout<<"finish fft"<<endl;

  vector< vector<double> > autoC(pix);
  for (ar=0; ar<pix; ar++) autoC[ar].resize(pix, 0.0);
  for (ar=0; ar<pix; ar++) {
    for (ac=0; ac<pix; ac++) {
//cout<<ar<<"/"<<ac<<" ";
//cout<<(centr+ar)%pix<<"/"<<(centc+ac)%pix<<" ";
      autoC[(centr+ar)%pix][(centc+ac)%pix] = autCfft[ar*pix + ac][0];
      //cout<<"autoC: "<<autoC[(centr+ar)%pix][(centc+ac)%pix]<<endl;
      //cout<<autCfft[ar*pix + ac][0]<<"/"<<autCfft[ar*pix + ac][1]<<"  ";
    }
cout<<endl<<endl;
  }

  plt.printRC(test, "autoTest");
  plt.printRC(autoC, "autoCorr");

*/


  fftw_destroy_plan(fftAC);    
  fftw_free(difPfft);
  fftw_free(autCfft);
//  for (int k=0; k<11; k++) fftw_free(fSpaces[k]);
cout<<"done"<<endl;
  output->Close();
  outTxt.close();
  return 1;
}

/*
  ////// Wiener Filter the entire image //////

  int N = lgndrs[2][0].size();
  int FTsize = (int)(N/2+1);
  int noiseKthresh = (int)((N/2+1)*3/5);
  double* tSpace = (double*) fftw_malloc(sizeof(double)*N);
  fftw_complex* fSpace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*FTsize);
  vector< vector< vector<complex<double> > > > fSpaces(11);
  vector< vector< vector<double> > > tSpaces(11);
  fftw_plan fftF;
  fftw_plan fftB;
  fftF = fftw_plan_dft_r2c_1d( N, tSpace, fSpace, FFTW_MEASURE);
  fftB = fftw_plan_dft_c2r_1d( N, fSpace, tSpace, FFTW_MEASURE);

  vector<double> linPart(FTsize-noiseKthresh, 0.0);
  TH1F* hist = new TH1F("fft", "fft", FTsize, 0, FTsize);
  TH1F* fitHist = new TH1F("lin", "lin", noiseKthresh, 0, noiseKthresh);
  TF1* linFit = new TF1("linFit", "[0] + [1]*x", 0, FTsize-noiseKthresh);

  for (int k=2; k<3; k++) {
    linPart.assign(linPart.size(), 0);
    fSpaces[k].resize(lgndrs[k].size());
    tSpaces[k].resize(lgndrs[k].size());
    for (ir=0; ir<(int)lgndrs[k].size(); ir++) {
      fSpaces[k][ir].resize(FTsize);
      // Fill fftw time vector
      for (it=0; it<(int)lgndrs[k][ir].size(); it++) {
        tSpace[it] = lgndrs[k][ir][it];
      }

      // FT time vector
      fftw_execute(fftF);
      // Save values to go back
      for (int itt=0; itt<FTsize; itt++) {
	fSpaces[k][ir][itt] = complex<double>(fSpace[itt][0], fSpace[itt][1]);
      }

//    double bAvg=0;
      // Find linear fit of noise from last 2/5 bins
      int cnt=linPart.size()-1;
      for (int itt=FTsize-1; itt>noiseKthresh; itt--) {
	linPart[cnt] += sqrt(norm(fSpaces[k][ir][itt]));  
//	bAvg += sqrt(norm(fSpaces[k][ir][itt]));
//        fitHist->SetBinContent(cnt+1, sqrt(norm(fSpaces[k][ir][itt])));
	cnt--;
      }

//    //for (uint il=0; il<linPart.size(); il++) bAvg += linPart[il];
//    linFit->SetParameter(0, bAvg/((double)linPart.size()));
//    linFit->SetParameter(1, 0);
//    fitHist->Fit("linFit", "q");
//cout<<"row slope intcpt: "<<ir<<"   "<<linFit->GetParameter(1)<<"   "<<linFit->GetParameter(0)<<"   "<<bAvg/((double)linPart.size())<<endl;     

      for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, sqrt(pow(fSpace[ib][0],2)+pow(fSpace[ib][1],2)));
      hist->SetMaximum(10000);
      //plt.print1d(hist,"fftRad"+to_string(ir));
    }


    // Fitting the noise
    for (uint ih=0; ih<linPart.size(); ih++) fitHist->SetBinContent(ih+1, linPart[ih]/((double)lgndrs[k].size()));
    linFit->SetParameter(0, linPart[0]);
    linFit->SetParameter(1, 0);
    fitHist->Fit("linFit", "q");
    intcpt = linFit->GetParameter(0);
    slope = linFit->GetParameter(1);
cout<<"lgndr slope intcpt: "<<k<<"   "<<slope<<"   "<<intcpt<<endl;

//    intcpt = 2989.33;
//    slope = 1130.94; 

    for (ir=0; ir<(int)fSpaces[k].size(); ir++) {
      // Subtracting the noise
      for (int ic=0; ic<FTsize; ic++) {
   	sbtr = slope*(ic-noiseKthresh) + intcpt;
        if (ic+1 > noiseKthresh || sbtr > sqrt(norm(fSpaces[k][ir][ic]))) {
          fSpaces[k][ir][ic] = complex<double> (0,0);
          continue;
        }
        scl = (sqrt(norm(fSpaces[k][ir][ic])) - sbtr)/sqrt(norm(fSpaces[k][ir][ic]));
        fSpaces[k][ir][ic] = complex<double> (fSpaces[k][ir][ic].real()*scl, fSpaces[k][ir][ic].imag()*scl);
      }

      for (int ib=0; ib<FTsize; ib++) hist->SetBinContent(ib+1, sqrt(norm(fSpaces[k][ir][ib])));
      hist->SetMaximum(10000);
      plt.print1d(hist,"fftWF"+to_string(ir));
   
      
      // Inverse FT to get cleaner image
      for (uint ic=0; ic<fSpaces[k][ir].size(); ic++) {
	fSpace[ic][0] = fSpaces[k][ir][ic].real();
	fSpace[ic][1] = fSpaces[k][ir][ic].imag();
      }
	
      fftw_execute(fftB);

      // Save values to go back
      tSpaces[k][ir].resize(N);
      for (int itt=0; itt<N; itt++) {
	tSpaces[k][ir][itt] = tSpace[itt];
      }
    }
    
 
    // Plotting filtered legendres
    plt.printRC(tSpaces[k], "WFlegendre_"+to_string(k));
*/

/*
          for (it=0; it<(int)(*lgndr2).size(); it++) {
	    lgndrs[0][k][it]=(*lgndr0)[ir][it];    lgndrs[1][k][it]=(*lgndr1)[it][ir];    lgndrs[2][k][it]=(*lgndr2)[it][ir];
	    lgndrs[3][k][it]=(*lgndr3)[it][ir];    lgndrs[4][k][it]=(*lgndr4)[it][ir];    lgndrs[5][k][it]=(*lgndr5)[it][ir];
	    lgndrs[6][k][it]=(*lgndr6)[it][ir];    lgndrs[7][k][it]=(*lgndr7)[it][ir];    lgndrs[8][k][it]=(*lgndr8)[it][ir];
	    lgndrs[9][k][it]=(*lgndr9)[it][ir];    lgndrs[10][k][it]=(*lgndr10)[it][ir];
          }
*/



