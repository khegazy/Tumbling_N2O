#include "/reg/neh/home/khegazy/baseTools/simulation/rigidRotor/RIGROTSIMclass.h"
#include <dirent.h>

enum expENUM {UED, LCLS, exp_1, exp_2, exp_3, exp_4, Nexps};
std::vector<std::string> expNames{"UED","LCLS","exp1","exp2","exp3","exp4"};


int main(int argc, char* argv[]) {


  RIGROTSIMclass rrs;
  PLOTclass plt;

  expENUM exp = UED;
  if (argc > 2) {
    for (int iarg=1; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg], "-Exp") == 0) {
        for (int i=0; i<(int)expNames.size(); i++) {
          std::string expN(argv[iarg+1]);
          if (expN.compare(expNames[i]) == 0) {
            exp = static_cast<expENUM>(i);
          }
        }
      }
    }
  }

  cout << "Running Experiment " << expNames[exp] <<endl;

  //// Variables to Change ////
                                                 
  rrs.startTime  = 0;               /* ps */    
  rrs.endTime    = 400;             /* ps */   
                                                 
  rrs.temperature    = 55; //60.30;         /* K */             
  rrs.laserIntensity = 2.5e12; //2.72e12;        /* W/cm^2 */  
                                                     
  rrs.rotConstB  = 0.419011;      /* cm^(-1) */      
  rrs.rotConstD  = 1.76e-7;       /* cm^(-1) */     
  rrs.rotConstH  = 0.16e-13;      /* cm^(-1) */    
  rrs.deltaAlpha = 2.994e-30;     /* m^3 */          
                                                    
  //double startSampleTime = 357.804-1.5+0.01;     /* ps */              
  //double endSampleTime = 357.804+1.5+0.01;       /* ps */              
  rrs.startSampleTime  = 316.114; //315.914; //355.726; //354.726;      /* ps */    
  rrs.endSampleTime    = 320.014; //320.214;//360.026;//362;       /* ps */
  rrs.sampleStep       = 1e-1;      /* ps */       
  rrs.dTimeEvolStep    = 1e-3;      /* ps */       
  rrs.pulseTlength     = 0.1;       /* ps */       
  rrs.Npulses          = 8;                       
  rrs.pulseTspacing    = 39.794;    /* ps */        
  //rrs.pulseTspacing = 39.757;     /* ps */         
                                                   
                                                   
  rrs.doCosSqEV     = true;                         
  rrs.doSinCosEV    = false;                         
  rrs.cosOrder      = 2;                            
  rrs.sinCosOrder   = 1;                            
  rrs.makePDFs      = false;                         
  rrs.NYlmEVs       = 10;
  rrs.evenOnlyAxisDist  = true;
  rrs.startPDFTime  = 316.114;          /* ps */
  rrs.endPDFTime    = 320.014;          /* ps */           
  //double startPDFTime = 357.804-1.5+0.01;
  //double endPDFTime = 357.804+1.5+0.01;          /* ps */      
                                                   
  rrs.MAXj           = 70;                        
  rrs.hasQuarterRev  = false;                      
  rrs.indexRefr      = 1;                        
  rrs.savePDFformat  = "ROOT";                   

  rrs.fileName  = "anglePDFs";
  rrs.outputDir = "./";
  rrs.PDFfileName  = "anglePDFs";
  rrs.PDFoutputDir = "/reg/d/psdm/amo/amoi0314/scratch/anglePDFs/";

  if (exp == UED) {
    rrs.PDFoutputDir = rrs.PDFoutputDir + "UED/";
  }

  /////  LCLS Parameters  /////
  if (exp == LCLS) {
    rrs.temperature      = 300;              /* K */           
    rrs.laserIntensity   = 2.5e12;           /* W/cm^2 */      
    rrs.pulseTspacing    = 39.840;           /* ps */          
    rrs.startSampleTime  = 318.384; //317.0; //357.045;     /* ps */     
    rrs.endSampleTime    = 319.39588; //320.50; //360.025;       /* ps */    
    rrs.sampleStep       = 4.936e-3; //4.2305e-03;      /* ps */           
    //rrs.sampleStep = 5e-03;      /* ps */            
    rrs.MAXj             = 100;                         
    rrs.savePDFformat    = "binary";                    
    rrs.vibKey           = "n2o";                       
    rrs.outputDir = "/reg/d/psdm/amo/amoi0314/scratch/simThetaPhi/";
    rrs.PDFoutputDir = rrs.PDFoutputDir + "LCLS/";
  }


  /////  Other Experiments  /////
  if (exp == exp_1) {
    rrs.temperature       = 300;              /* K */           
    rrs.laserIntensity    = 2.5e12;           /* W/cm^2 */      
    rrs.pulseTspacing     = 39.840;           /* ps */          
    rrs.startSampleTime   = 315.0;            /* ps */     
    rrs.endSampleTime     = 321.0;            /* ps */    
    rrs.sampleStep        = 2.5e-2;             /* ps */       
    rrs.NYlmEVs           = 10;
    rrs.makePDFs          = true;                         

    rrs.Npulses           = 8;                       
    rrs.MAXj              = 100;
    rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_1/";
  }

  if (exp == exp_2) {
    rrs.temperature       = 300;              /* K */           
    rrs.laserIntensity    = 2.5e13;           /* W/cm^2 */      
    rrs.pulseTspacing     = 39.840;           /* ps */          
    rrs.startSampleTime   = 315.0;            /* ps */     
    rrs.endSampleTime     = 321.0;            /* ps */    
    rrs.sampleStep        = 2.5e-2;             /* ps */       
    rrs.NYlmEVs           = 10;
    rrs.makePDFs          = true;                         

    rrs.Npulses           = 8;                       
    rrs.MAXj              = 100;
    rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_2/";
  }

  if (exp == exp_3) {
    rrs.temperature       = 60;               /* K */           
    rrs.laserIntensity    = 2.5e13;           /* W/cm^2 */      
    rrs.pulseTspacing     = 39.794;           /* ps */        
    rrs.startSampleTime   = 315.0;            /* ps */     
    rrs.endSampleTime     = 321.0;            /* ps */    
    rrs.sampleStep        = 2.5e-2;             /* ps */       
    rrs.NYlmEVs           = 10;
    rrs.makePDFs          = true;                         

    rrs.Npulses           = 8;                       
    rrs.MAXj              = 100;
    rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_3/";
  }

  if (exp == exp_4) {
    rrs.temperature       = 60;               /* K */           
    rrs.laserIntensity    = 2.5e12;           /* W/cm^2 */      
    rrs.pulseTspacing     = 39.794;           /* ps */        
    rrs.startSampleTime   = 315.0;            /* ps */     
    rrs.endSampleTime     = 321.0;            /* ps */    
    rrs.sampleStep        = 2.5e-2;             /* ps */       
    rrs.NYlmEVs           = 10;
    rrs.makePDFs          = true;                         

    rrs.Npulses           = 8;                       
    rrs.MAXj              = 100;
    rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_4/";
  }



  /////  Plotting Parameters  /////
  vector<PLOToptions> opts(3);
  vector<string> optVals(3);
  opts[0] = xSpan;              optVals[0] = to_string(rrs.startSampleTime)
						+","+to_string(rrs.endSampleTime);
  opts[1] = xLabel;             optVals[1] = "Time [ps]";
  opts[2] = yLabel;             

  /////  Command Line Parameters  /////
  if (argc > 2) {
    for (int iarg=1; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg],"-Ofile")==0) {
	string str(argv[iarg+1]); 
	rrs.fileName=str;
	cout << "INFO: Replacing fileName with " << str << "!!!\n";
      }
      else if (strcmp(argv[iarg],"-Odir")==0) {
	string str(argv[iarg+1]); 
	rrs.outputDir=str;
	cout << "INFO: Replacing output dir with " << str << "!!!\n";
      }
      else if (strcmp(argv[iarg],"-Temp")==0) {
        rrs.temperature = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg],"-Ints")==0) {
        rrs.laserIntensity = atof(argv[iarg+1])*1e12;
      }
      else if (strcmp(argv[iarg],"-SampleTime")==0) {
        rrs.startPDFTime  = atof(argv[iarg+1]);
        rrs.endPDFTime    = rrs.startPDFTime;
      }
      else if (strcmp(argv[iarg],"-StartTime")==0) {
        rrs.startSampleTime = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg],"-EndTime")==0) {
        rrs.endSampleTime = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg],"-SampleStep")==0) {
        rrs.sampleStep = atof(argv[iarg+1])*1e-3;
      }
      else if (strcmp(argv[iarg], "-Exp") == 0) {
        continue;
      }
      else {
        cerr << "ERROR!!! Option " << argv[iarg] << " does not exist!" << endl;
        exit(0);
      }
    }
  }



  ////////////////////////////
  /////  Run Simulation  /////
  ////////////////////////////

  cout<<"starting simulation"<<endl;
  rrs.runSimulation();

  cout<<"start plotting pop"<<endl;
  /////  Plot Population  /////
  std::vector<PLOToptions> Popts(3);
  std::vector<std::string> Pvals(3);
  Popts[0] = logy;    Pvals[0] = "true";
  Popts[1] = xLabel;  Pvals[1] = "J";
  Popts[2] = yLabel;  Pvals[2] = "Population";
  std::vector<double> population1 = rrs.getPopulationDistribution();

  /////  Saving Projection Coeffs onto Ylm  /////
  if (rrs.NYlmEVs) {
    int j;
    for (int j_=0; j_<=rrs.NYlmEVs; j_++) {
      if (rrs.evenOnlyAxisDist) {
        j = 2*j_;
      }
      else {
        j = j_;
      }
      std::string folderName = 
          "/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/YlmExpVals/"
          + expNames[exp] + "_"
          + to_string(rrs.temperature) + "_"
          + to_string(rrs.laserIntensity*1e-12) + "_"
          + to_string(rrs.startSampleTime) + "_"
          + to_string(rrs.endSampleTime);

      if (!opendir(folderName.c_str())) {
        mkdir(folderName.c_str(), 0777);
      }

      save::saveDat<double>(rrs.YlmEVals[j_],
          folderName + "/"
          + "expValYlm_L-" + to_string(j)
          + "_M-0_time-" + to_string(rrs.startSampleTime)
          + "-" + to_string(rrs.endSampleTime)
          + "_bins[" + to_string(rrs.YlmEVals[j_].size())
          + "].dat");
      cout<<"saved"<<endl;
    }
    cout<<"done saving wave fxn"<<endl;
  }


  cout<<"saving cos2"<<endl;
  /////  Saving Expectation Values  ///
  if (rrs.doCosSqEV) {
    std::ofstream asciiFile;
    std::string outputFolder = "/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/cosExpVals/"
          + expNames[exp] + "/";
    //std::string outputFolder = "/reg/neh/home5/khegazy/analysis/tumblingN2O/movie/data/";
    //if (exp == LCLS) {
    //  outputFolder = "/reg/d/psdm/amo/amoi0314/scratch/cosExpVals/";
    //}
    for (uint i=0; i<rrs.cosSqEVals.size(); i++) {
      FILE* file = fopen((outputFolder + "expValCos" + to_string(2 + i*2) 
            + "_" + to_string(rrs.temperature)
            + "_" + to_string(rrs.laserIntensity*1e-12)
            + "_" + to_string(rrs.startSampleTime) 
            + "-" + to_string(rrs.endSampleTime) 
            + "-" + to_string(rrs.sampleStep) + ".dat").c_str()
                         , "wb");
      fwrite(&rrs.cosSqEVals[i][0], sizeof(double), 
          rrs.cosSqEVals[i].size(), file);

      if (exp = LCLS) {
        asciiFile.open((outputFolder + "expValCos" + to_string(2 + i*2)
              + "_" + to_string(rrs.startSampleTime) 
              + "-" + to_string(rrs.endSampleTime) + ".ascii").c_str());
        for (uint j=0; j<rrs.cosSqEVals[i].size(); j++) {
          asciiFile << rrs.cosSqEVals[i][j];
          if (j != rrs.cosSqEVals[i].size() - 1) {
            asciiFile << "\t";
          }
        }
        asciiFile.close();
      }
    }
  }
 
  cout<<"done"<<endl;
  cout<<"sincos"<<endl;
  if (rrs.doSinCosEV) {
    std::ofstream asciiFile;
    std::string outputFolder = "output/";
    //std::string outputFolder = "/reg/neh/home5/khegazy/analysis/tumblingN2O/movie/data/";
    if (exp == LCLS) {
      outputFolder = "/reg/d/psdm/amo/amoi0314/scratch/sinCosExpVals/";
    }
    for (uint i=0; i<rrs.sinCosEVals.size(); i++) {
      FILE* file = fopen((outputFolder + "expValSinCos" + to_string(2 + i*2) 
            + "_" + to_string(rrs.startSampleTime) 
            + "-" + to_string(rrs.endSampleTime) 
            + "-" + to_string(rrs.sampleStep) + ".dat").c_str()
                         , "wb");
      fwrite(&rrs.sinCosEVals[i][0], sizeof(double), 
          rrs.sinCosEVals[i].size(), file);

      if (exp == LCLS) {
        asciiFile.open((outputFolder + "expValSinCos" + to_string(2 + i*2)
              + "_" + to_string(rrs.startSampleTime) 
              + "-" + to_string(rrs.endSampleTime) + ".ascii").c_str());
        for (uint j=0; j<rrs.sinCosEVals[i].size(); j++) {
          asciiFile << rrs.sinCosEVals[i][j];
          if (j != rrs.sinCosEVals[i].size() - 1) {
            asciiFile << "\t";
          }
        }
        asciiFile.close();
      }
    }
  }
  cout<<"done"<<endl;

  cout<<"plotting"<<endl;
  /////  Plotting Results  /////
  if (rrs.doCosSqEV && false) {
    optVals[2] = "<Cos^{2}#theta>";
    for (uint i=0; i<rrs.cosSqEVals.size(); i++) {
      delete plt.print1d(rrs.cosSqEVals[i], 
            "./plots/cos" + to_string(2 + 2*i) 
            + "_" + to_string(rrs.startSampleTime) 
            + "-" + to_string(rrs.endSampleTime),
            opts, optVals);
    }
  }
  if (rrs.doSinCosEV) {
    optVals[2] = "<Sin(#theta)Cos(#theta)>";
    for (uint i=0; i<rrs.sinCosEVals.size(); i++) {
      delete plt.print1d(rrs.sinCosEVals[i], 
            "./plots/sinCos" + to_string(2 + 2*i) 
            + "_" + to_string(rrs.startSampleTime) 
            + "-" + to_string(rrs.endSampleTime),
            opts, optVals);
    }
  }
  cout<<"done"<<endl;


  return 0;
}
