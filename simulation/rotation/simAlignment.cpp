#include "/reg/neh/home/khegazy/baseTools/simulation/rigidRotor/RIGROTSIMclass.h"
#include <dirent.h>

//enum expENUM {UED, LCLS, exp_1, exp_2, exp_3, exp_4, Nexps};
//std::vector<std::string> expNames{"UED","LCLS","exp1","exp2","exp3","exp4"};


int main(int argc, char* argv[]) {


  RIGROTSIMclass rrs;
  PLOTclass plt;

  //expENUM exp = UED;
  std::string exp = "UED";
  if (argc > 2) {
    for (int iarg=1; iarg<argc; iarg+=2) {
      if (strcmp(argv[iarg], "-Exp") == 0) {
        std::string temp(argv[iarg+1]);
        exp = temp;

        //for (int i=0; i<(int)expNames.size(); i++) {
        //  std::string expN(argv[iarg+1]);
        //  if (expN.compare(expNames[i]) == 0) {
        //    exp = static_cast<expENUM>(i);
        //  }
        //}
      }
    }
  }

  //cout << "Running Experiment " << expNames[exp] <<endl;
  cout << "Running Experiment " << exp <<endl;

  //// Variables to Change ////
                                                 
  rrs.startTime  = 0;               /* ps */    
  rrs.endTime    = 400;             /* ps */   
                                                 
  rrs.rotConstB  = 0.419011;      /* cm^(-1) */      
  rrs.rotConstD  = 1.76e-7;       /* cm^(-1) */     
  rrs.rotConstH  = 0.16e-13;      /* cm^(-1) */    
  rrs.deltaAlpha = 2.994e-30;     /* m^3 */          
  rrs.vibKey     = "n2o";                       
                                                    
  //double startSampleTime = 357.804-1.5+0.01;     /* ps */              
  //double endSampleTime = 357.804+1.5+0.01;       /* ps */              
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
  rrs.NYlmEVs       = 5;
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
  //rrs.PDFoutputDir = "/reg/d/psdm/amo/amoi0314/scratch/anglePDFs/";
  rrs.PDFoutputDir = "/reg/ued/ana/scratch/N2O/alignmentPDFs/";
  //rrs.PDFoutputDir = "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/temp/alignmentPDFs/";

  std::vector<double> exp_temp{ 50,   100,  250,  500,  1000};
  std::vector<double> exp_ints{1e11,  5e11, 1e12, 5e12, 1e13};
  std::vector<double> exp_strT{ 316,    316.00, 316.00, 316.00,  316.00,
                                316,    316,    316,    316.00,  316.25,
                                316.25, 316.25, 316.25, 316.25,  316.50,
                                316.75, 316.75, 316.75, 316.75,  316.75,
                                317.00, 317.00, 317.00, 317.00,  317.00};
  double delta;
  std::vector<double> exp_endT;
  for (uint i=0; i<exp_strT.size(); i++) {
    delta = (exp_strT[i] - 316)/8;
    exp_strT[i] = exp_strT[i]/8 - 2 + delta;//-= 276.5;
    exp_endT.push_back(exp_strT[i] + 5);
  }


  if (exp.compare("UED") == 0) {
    rrs.temperature    = 84; //55; //60.30;         /* K */             
    rrs.laserIntensity = 3.5e12; //2.5e12; //2.72e12;        /* W/cm^2 */  
    rrs.pulseTspacing    = 39.794;    /* ps */        
    rrs.startSampleTime  = 316.000; //316.114; //315.914; //355.726; //354.726;      /* ps */    
    rrs.endSampleTime    = 319.900; //320.214;//360.026;//362;       /* ps */
    rrs.sampleStep       = 1e-1;      /* ps */       
    rrs.PDFoutputDir = rrs.PDFoutputDir + "UED/";
  }
  else if (exp.compare("LCLS") == 0) {
    rrs.temperature      = 300;              /* K */           
    rrs.laserIntensity   = 2.5e12;           /* W/cm^2 */      
    rrs.pulseTspacing    = 39.840;           /* ps */          
    rrs.startSampleTime  = 317; //318.384; //317.0; //357.045;     /* ps */     
    rrs.endSampleTime    = 321; //319.39588; //320.50; //360.025;       /* ps */    
    rrs.sampleStep       = 4.936e-3; //4.2305e-03;      /* ps */           
    //rrs.sampleStep = 5e-03;      /* ps */            
    rrs.MAXj             = 100;                         
    rrs.savePDFformat    = "binary";                    
    rrs.vibKey           = "n2o";                       
    //rrs.outputDir = "/reg/d/psdm/amo/amoi0314/scratch/simThetaPhi/";
    rrs.PDFoutputDir = rrs.PDFoutputDir + "LCLS/";
  }
  else if (exp.find("exp") != std::string::npos) {
    int num = stoi(exp.substr(3, exp.length() - 3));
    int intsInd = num % exp_ints.size();
    int tempInd = (num - intsInd)/exp_temp.size();

    rrs.Npulses         = 1;
    rrs.MAXj            = 100;                        
    rrs.temperature     = exp_temp[tempInd];
    rrs.laserIntensity  = exp_ints[intsInd];
    rrs.sampleStep      = 2.5e-2;      /* ps */       
    rrs.startSampleTime = exp_strT[num];
    rrs.endSampleTime   = exp_endT[num];
    rrs.PDFoutputDir    = rrs.PDFoutputDir + exp + "/";
    cout<<"EXP INFO: "<<exp<<" "<<num<<" "<<intsInd<<" "<<tempInd<<" "<<rrs.temperature<<" "<<rrs.laserIntensity<<endl;
  }
  else {
    cerr << "Cannot handle experiment named " << exp <<endl;
    std::exit(0);
  }

//  /////  Other Experiments  /////
//  if (exp == exp_1) {
//    rrs.temperature       = 300;              /* K */           
//    rrs.laserIntensity    = 2.5e12;           /* W/cm^2 */      
//    rrs.pulseTspacing     = 39.840;           /* ps */          
//    rrs.startSampleTime   = 315.0;            /* ps */     
//    rrs.endSampleTime     = 321.0;            /* ps */    
//    rrs.sampleStep        = 2.5e-2;             /* ps */       
//    rrs.NYlmEVs           = 10;
 //   rrs.makePDFs          = true;                         
 //
 //   rrs.Npulses           = 8;                       
 //   rrs.MAXj              = 100;
 //   rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_1/";
 // }
 //
 // if (exp == exp_2) {
 //   rrs.temperature       = 300;              /* K */           
 //   rrs.laserIntensity    = 2.5e13;           /* W/cm^2 */      
 //   rrs.pulseTspacing     = 39.840;           /* ps */          
 //   rrs.startSampleTime   = 315.0;            /* ps */     
 //   rrs.endSampleTime     = 321.0;            /* ps */    
 //   rrs.sampleStep        = 2.5e-2;             /* ps */       
 //   rrs.NYlmEVs           = 10;
 //   rrs.makePDFs          = true;                         
 //
 //   rrs.Npulses           = 8;                       
 //   rrs.MAXj              = 100;
 //   rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_2/";
 // }
 //
 // if (exp == exp_3) {
 //   rrs.temperature       = 60;               /* K */           
 //   rrs.laserIntensity    = 2.5e13;           /* W/cm^2 */      
 //   rrs.pulseTspacing     = 39.794;           /* ps */        
 //   rrs.startSampleTime   = 315.0;            /* ps */     
 //   rrs.endSampleTime     = 321.0;            /* ps */    
 //   rrs.sampleStep        = 2.5e-2;             /* ps */       
 //   rrs.NYlmEVs           = 10;
 //   rrs.makePDFs          = true;                         
 //
 //   rrs.Npulses           = 8;                       
 //   rrs.MAXj              = 100;
 //   rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_3/";
 // }
//
//  if (exp == exp_4) {
//    rrs.temperature       = 60;               /* K */           
//    rrs.laserIntensity    = 2.5e12;           /* W/cm^2 */      
//    rrs.pulseTspacing     = 39.794;           /* ps */        
//    rrs.startSampleTime   = 315.0;            /* ps */     
//    rrs.endSampleTime     = 321.0;            /* ps */    
//    rrs.sampleStep        = 2.5e-2;             /* ps */       
//    rrs.NYlmEVs           = 10;
 //   rrs.makePDFs          = true;                         
 //
 //   rrs.Npulses           = 8;                       
 //   rrs.MAXj              = 100;
 //   rrs.PDFoutputDir = rrs.PDFoutputDir + "exp_4/";
 // }
 
  //std::string baseDir = "/reg/d/psdm/amo/amoi0314/scratch/UED_tumbling/";
  std::string baseDir = "/reg/ued/ana/scratch/N2O/";
  //std::string baseDir = "/reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/temp/";


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
        cout<<"Sampletime: "<<rrs.endPDFTime<<endl;
      }
      else if (strcmp(argv[iarg], "-SampleTimeOffset") == 0) {
        rrs.startPDFTime  = rrs.startSampleTime + atof(argv[iarg+1]);
        rrs.endPDFTime    = rrs.startPDFTime;
        cout<<"Sampletime: "<<rrs.endPDFTime<<endl;
      }
      else if (strcmp(argv[iarg],"-StartTime")==0) {
        rrs.startSampleTime = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg],"-EndTime")==0) {
        rrs.endSampleTime = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg],"-SampleStep")==0) {
        rrs.sampleStep = atof(argv[iarg+1]);
      }
      else if (strcmp(argv[iarg], "-PDForBases") == 0) {
        cout<<"TESTING: "<<argv[iarg+1]<<"  "<<atoi(argv[iarg+1])<<"  "<<(bool)atoi(argv[iarg+1])<<endl;
        rrs.makePDFs = (bool)atoi(argv[iarg+1]);
        if (rrs.makePDFs) {
          cout<<" IT IS TRUE"<<endl;
        }
        rrs.doCosSqEV = !rrs.makePDFs;
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
  if (rrs.makePDFs) {
    rrs.NYlmEVs = 0;
    rrs.doCosSqEV = false;
  }

  // Exit if YlmEVs exist
  if (rrs.NYlmEVs) {
    std::string folderName = 
        baseDir + "YlmExpVals/"
        //+ expNames[exp] + "_"
        + exp + "/" + exp + "_"
        + to_string(rrs.temperature) + "_"
        + to_string(rrs.laserIntensity*1e-12) + "_"
        + to_string(rrs.startSampleTime) + "_"
        + to_string(rrs.endSampleTime);

    if (!opendir(folderName.c_str())) {
      mkdir(folderName.c_str(), 0777);
    }

    bool exists = true;
    int j;
    for (int j_=0; j_<=rrs.NYlmEVs; j_++) {
      if (rrs.evenOnlyAxisDist) {
        j = 2*j_;
      }
      else {
        j = j_;
      }
      std::string fileName = folderName + "/"
          + "expValYlm_L-" + to_string(j)
          + "_M-0_time-" + to_string(rrs.startSampleTime)
          + "-" + to_string(rrs.endSampleTime)
          + "_bins[" + to_string(1+int((rrs.endSampleTime-rrs.startSampleTime)/rrs.sampleStep))
          + "].dat";


      exists = exists && (tools::fileExists(fileName));
    }
    if (exists && false) {
      cout << "<Ylm> already exist!" << endl;
      exit(0);
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
    cout << "Starting to save <Ylm>" << endl;
    std::string folderName = 
        baseDir + "YlmExpVals/"
        //+ expNames[exp] + "_"
        + exp;
   
    cout << folderName << endl;
    if (!opendir(folderName.c_str())) {
      mkdir(folderName.c_str(), 0777);
      sleep(10);
    }

    folderName = folderName + "/" + exp + "_"
        + to_string(rrs.temperature) + "_"
        + to_string(rrs.laserIntensity*1e-12) + "_"
        + to_string(rrs.startSampleTime) + "_"
        + to_string(rrs.endSampleTime);

    if (!opendir(folderName.c_str())) {
      mkdir(folderName.c_str(), 0777);
    }
     cout << "FOLDER: "<<folderName<<endl;

    int j;
    for (int j_=0; j_<=rrs.NYlmEVs; j_++) {
      if (rrs.evenOnlyAxisDist) {
        j = 2*j_;
      }
      else {
        j = j_;
      }

      save::saveDat<double>(rrs.YlmEVals[j_],
          folderName + "/"
          + "expValYlm_L-" + to_string(j)
          + "_M-0_time-" + to_string(rrs.startSampleTime)
          + "-" + to_string(rrs.endSampleTime)
          + "_bins[" + to_string(rrs.YlmEVals[j_].size())
          + "].dat");
      cout<<"saved"<<endl;
      /*
      delete plt.print1d(rrs.YlmEVals[j_], 
            "./plots/cos" + to_string(j) 
            + "_exp-" + exp
            + "_testing",
            opts, optVals);
      */

    }
    cout<<"done saving wave fxn"<<endl;
  }


  cout<<"saving cos2"<<endl;
  /////  Saving Expectation Values  ///
  if (false && rrs.doCosSqEV) {
    std::ofstream asciiFile;
    std::string outputFolder = baseDir +  "cosExpVals/"
          + exp + "/";
          //+ expNames[exp] + "/";
    if (!opendir(outputFolder.c_str())) {
      mkdir(outputFolder.c_str(), 0777);
    }
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

      if (exp.compare("LCLS") == 0) {
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
    if (exp.compare("LCLS")) {
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

      if (exp.compare("LCLS") == 0) {
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
  if (false && rrs.doCosSqEV) {
    optVals[2] = "<Cos^{2}#theta>";
    for (uint i=0; i<rrs.cosSqEVals.size(); i++) {
      delete plt.print1d(rrs.cosSqEVals[i], 
            "./plots/cos" + to_string(2 + 2*i) 
            + "_exp-" + exp
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
