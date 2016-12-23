/*******************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: Presently in the github
 *
 *
 * Authors:
 *
 *   Osipenkov, Ilya, Texas A&M - ilyao@fnal.gov
 *
 * Description: Use to constuct, fit and analyze toy datasets
 *
 ********************************************************************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TString.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TH1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom3.h>

using namespace std;
const bool Verbose=true;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////---------- Abstract Base Class used to hold information and perform operations on the individual fit results ------------////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GeneralFitResult{
public:
  GeneralFitResult(int nProc, int nVars): NProc(nProc), NVars(nVars) {
    theoryYield.clear();
    yieldFromGen.clear();
    ProcIdx.clear();
    for (int i=0; i<NProc; i++) {
      theoryYield.push_back(-1.0);
      yieldFromGen.push_back(-1.0);
      yieldToGen.push_back(-1.0);
      ProcIdx.push_back(-1);
    }
    Convergence=-1;
  }

  ///NProc and NVar reperesent external information which must be set (via a constructor) for each derived class
  const int NProc; //Number of processes (i.e. Signal + Backgrounds) in the fit
  const int NVars; //Number of variables (i.e. processes + free parameters) in the fit
  vector<TString> ProcName; //Names for each process used in the log file. Must be assigned for each derived class
  virtual void assignProcNames() = 0;
  
  ///Fit Input and Output
  vector<double> theoryYield, yieldFromGen, yieldToGen;
  vector<int> ProcIdx;

  TMatrixDSym* fitFullCov;//The Full Covariance Matrix
  TVectorD* fitPar;//Fitted Parameter Values
  TVectorD* fitErr;//Their Errors
  double fitTotal, fitTotalErr, dataTotal, chi2dof, chi2prob;
  int Convergence; //0: Full Convergence, 1: Did Not Converge, 2: Converged, but the Covariance Matrix (from HESSE) is not positive definite, -1: Initialization Value/Did not record convergence

  void printSummary();
  void computeYieldsToGen(int InitRand=333, bool includeNonPoissonFluctuations = true, bool skipConfigurationsWithNegativeYields = true);
  /// Delete the pointer content
  virtual ~GeneralFitResult() {
    delete fitFullCov;
    delete fitPar;
    delete fitErr;
  }

protected:
  void readLogInfo(TString);

};

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void GeneralFitResult::printSummary() {
  for (int i=0; i<NProc; i++) {
    cout << "Process " << i << " : Name=" << ProcName[i] << ", PredictedEventYield=" << theoryYield[i] << ", GeneratedEventYield=" << yieldFromGen[i] << ", ProcessIndex=" << ProcIdx[i] <<endl;
    cout << "FitFraction=" << (*fitPar)(ProcIdx[i]) << " +/- " << (*fitErr)(ProcIdx[i]) << endl;
    cout << "YieldToBeUsedInToyGeneration=" << yieldToGen[i] << endl;
  }
  cout << "Covariance Matrix = " << endl;
  for (int j=0; j<NVars;j++) {
    for (int k=0; k<NVars;k++) {
      cout << (*fitFullCov)(j,k) << " ";
    }
    cout << endl;
  }
  cout << "fitTotal=" << fitTotal << " +/- " << fitTotalErr << ", dataTotal=" << dataTotal << ", chi2dof=" << chi2dof << ", chi2prob=" << chi2prob << ", Convergence=" << Convergence << endl;
  cout << endl;

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void GeneralFitResult::computeYieldsToGen(int InitRand, bool includeNonPoissonFluctuations, bool skipConfigurationsWithNegativeYields)
//// Recomputes the yields to be used for toy generation in appropriate amounts, randomized based on corresponding errors (when includeNonPoissonFluctuations = true) and then according to a Poisson distribution.
{
  TRandom PRand(InitRand);
  TRandom GRand(InitRand*137+2015);
  double smearedProcNorm[NProc];
  bool repeatGeneration=true;

  TMatrixDSymEigen CM(*fitFullCov);

  if ( Verbose ) {
    cout << "Recomputing YieldsToGen for toy generation (InitRand=" << InitRand << ")" << endl;
    if ( Convergence!=0 ) {
      cout << "WARNING : Convergence problems in the original fit" << endl;
    }
    cout << "Covariance Matrix = " << endl;
    for (int j=0; j<NVars;j++) {
      for (int k=0; k<NVars;k++) {
	cout << (*fitFullCov)(j,k) << " ";
      }
      cout << endl;
    }

    cout << "Expected Variable Values = " << endl;
    for (int j=0; j<NVars;j++) {
      cout << (*fitPar)(j) << " ";
    }
    cout << endl;
  }
  TVectorD sigsqVal(NVars);
  sigsqVal = CM.GetEigenValues();
  TVectorD errVal(NVars);
  TMatrixD EVec = CM.GetEigenVectors();
  for (int j=0; j<NVars;j++) {
    errVal(j)=TMath::Sqrt(sigsqVal(j));
  }
  TVectorD RErr(NVars);//Errors (to be generated) in the Rotated System
  TVectorD UnRErr(NVars);//Errors in the physical coordinate system

  while ( repeatGeneration ) {
    repeatGeneration=false;
    for (int j=0; j<NVars;j++) {
      RErr(j)=GRand.Gaus(0,errVal(j));//Fluctuation in the rotated coordinate system
    }
    UnRErr=EVec*RErr;//Fluctuation in the physical coordinate system
    for (int i=0; i<NProc; i++) {
      //Get the generated yields prior to adding poisson fluctuations
      smearedProcNorm[i]=(*fitPar)(ProcIdx[i])+UnRErr(ProcIdx[i]);
      if ( includeNonPoissonFluctuations ) {
	yieldToGen[i]=smearedProcNorm[i]*theoryYield[i];
      } else {
	yieldToGen[i]=(*fitPar)(ProcIdx[i])*theoryYield[i];//Take the predicted event count from the fit
      }
    }

    for (int i=0; i<NProc; i++) {
      //Either skip the configuration or apply poisson smearing
      if ( (yieldToGen[i]<0)&&(skipConfigurationsWithNegativeYields) ) {
	repeatGeneration=true;
	break;
      } else {
	yieldToGen[i]=PRand.Poisson(yieldToGen[i]);
      }
    }
  }

  if ( Verbose ) {
    /// Display the process information before and after smearing
    for (int i=0; i<NProc; i++) {
      cout << ProcName[i] << ": Ntheory=" << theoryYield[i] << " initial fit_nrm=" << (*fitPar)(ProcIdx[i]) << ", smeared fit_nrm=" << smearedProcNorm[i] << ",  Ngen=" << yieldToGen[i] << endl;
    }
  }

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void GeneralFitResult::readLogInfo(TString inLogFileName) {
  /// Go through the log file, identify key strings and record subsequent information
  TString sPar1, sPar2, sPar3, sPar4, sPar5;
  double dParV[NVars];
  double dPar1;
  char logline[2000];
  cout << "Processing Log File = " << inLogFileName << endl;
  ifstream inLogFile(inLogFileName);
  bool isCovMPositiveDefinite=true;
  bool genParsUsed=false;
  TMatrixDSym* ptFullCov = new TMatrixDSym(NVars);//The Full Covariance Matrix 
  TVectorD* ptFitPar = new TVectorD(NVars);//Fitted Parameter Values
  TVectorD* ptFitErr = new TVectorD(NVars);//Their Errors
  while ( inLogFile.good() ) {
    inLogFile.getline(logline,2000);
    istringstream str(logline);
    str >> sPar1 >> sPar2 >> sPar3 >> sPar4;
    if ( (sPar2=="fit")&&(sPar4=="null") ) { break; }//Null Fit should be performed after the default (if at all). All information from it is ignored; as will chi2 info, which is printed after the Null fit results.

    ///Get the initial values
    if ( (sPar1=="load")&&(sPar2=="data") ) {
      inLogFile.getline(logline,2000);//skip a line
      bool getadditionalyields=true;
      while ( getadditionalyields ) {
	///Parse the string beginning with 'RooRealVar::n_', either with or without skipping the string starting with 'explicitly'. Otherwise discontinue.
	inLogFile.getline(logline,2000);
	istringstream strA(logline);
	strA >> sPar1 >> sPar2 >> dPar1;
	if ( sPar1=="explicitly" ) {
     	  inLogFile.getline(logline,2000);
	  istringstream strB(logline);
	  strB >> sPar1 >> sPar2 >> dPar1;
	}
	for (int i=0; i<NProc; i++) {
	  if ( sPar1.Contains("n_"+ProcName[i]) ) { theoryYield[i]=dPar1; }
	}
	if ( !sPar1.Contains("RooRealVar::") ) { getadditionalyields=false; }
      }
    }

    ///Get the number of generated events
    if ( sPar1=="genPars:" ) { genParsUsed=true; }
    if ( (sPar3=="RooRealVar::")&&genParsUsed ) {
      str >> sPar5 >> dPar1;
      for (int i=0; i<NProc; i++) {
	if ( sPar4==("n_"+ProcName[i]) ) { yieldFromGen[i]=dPar1; }
      }
    }


    /// Get the Covariance Matrix from the fit
    if ( sPar3=="**HESSE" ) {
      ///Step 1: Check the Convergence Status from HESSE
      for (int dummy=0; dummy<3; dummy++) {
    	inLogFile.getline(logline,2000);
    	istringstream stra(logline);
    	stra >> sPar1 >> sPar2 >> sPar3 >> sPar4;
    	if ( dummy==1 ) {
    	  if ( (sPar3=="CALCULATED")&&(sPar4=="SUCCESSFULLY") ) {
    	    isCovMPositiveDefinite=true;
    	    //cout << "Covariance Matrix is Positive Definite" << endl;
    	  } else {
    	    isCovMPositiveDefinite=false;
    	    cout << "Covariance Matrix is NOT Positive Definite" << endl;
    	  }
    	}
    	if ( sPar2=="ERROR:InputArguments" ) {
    	  cout << "Skipping Error Messages" << endl;
    	  bool skiperrors=true;
    	  while (skiperrors) {
    	    inLogFile.getline(logline,2000);
    	    istringstream strb(logline);
    	    strb >> sPar1 >> sPar2 >> sPar3 >> sPar4;
    	    if ( sPar2!="ERROR:InputArguments" ) { skiperrors=false;}
    	  }
    	}
      }
      if ( isCovMPositiveDefinite ) {
    	istringstream strA(logline);
    	strA >> sPar1 >> sPar2 >> sPar3 >> sPar4;
    	cout << " HESSE_Convergence = " << sPar4 << endl;
      }

      ///Step 2: Get the Covariance Matrix
      int nSkip=5;
      if ( !isCovMPositiveDefinite ) { nSkip=4; }//no "COVARIANCE MATRIX CALCULATED SUCCESSFULLY" line
      for (int j=0; j<nSkip+NVars; j++) {
    	inLogFile.getline(logline,2000);
    	istringstream strc(logline);
    	strc >> sPar1 >> sPar2 >> sPar3 >> sPar4;
    	if ( sPar1=="WARNING" ) {
    	  cout << "Skipping WARNING - AT LIMIT message" << endl;
    	  inLogFile.getline(logline,2000);
    	}
      }
      for (int j=0; j<NVars; j++) {
    	inLogFile.getline(logline,2000);
    	istringstream strd(logline);
    	strd >> sPar1 >> sPar2;
    	if ( (sPar1=="ELEMENTS")&&(sPar2=="ABOVE") ) {
    	  inLogFile.getline(logline,2000);
    	}
    	istringstream strB(logline);
	//cout << "NVars=" << NVars << endl;
	for (int k=0; k<NVars; k++) {
	  //cout << "k=" << k << endl;
    	  strB >> dParV[k];
    	  (*ptFullCov)(j,k)=dParV[k];
    	}
      }
      /// Make sure the values below the diagonal are the ones kept (the ones above aren't always recorded in the log file)
      for (int j=0; j<NVars; j++) {
    	for (int k=0; k<NVars; k++) {
    	  if ( (*ptFullCov)(j,k)!=(*ptFullCov)(k,j) ) {
    	    cout << "Symmetrizing the covariance matrix" << endl;
    	    (*ptFullCov)(j,k)=(*ptFullCov)(k,j);
    	  }
    	}
      }
    }

    ///Get the Convergence Status from MIGRAD and HESSE
    if ( (sPar1=="Status")&&(sPar2=":") ) {
      if ( (sPar3=="MIGRAD=0")&&(sPar4=="HESSE=0") ) {
    	Convergence=0;
    	if ( !isCovMPositiveDefinite ) {
    	  ///Note: Covariance Matrix gets processed before the convergence status
    	  Convergence=2;
    	  cout << "Fit converged but did not generate a positive definite error matrix" << endl;
    	}
      } else {
    	Convergence=1;
    	cout << "Fit did not converge: " << sPar2 << ", " << sPar3 << endl;
      }
    }
    fitFullCov=ptFullCov;

    ///Get the Total Yield from the fit and in the data
    if ( (sPar1=="total")&&(sPar2=="expected:") ) {
	istringstream strA(logline);
	strA >> sPar1 >> sPar2 >> dPar1;
	fitTotal=dPar1;
	inLogFile.getline(logline,2000);
	istringstream strB(logline);
	strB >> sPar1 >> sPar2 >> dPar1;
	dataTotal=dPar1;
    }

    ///Get the parameter values from the fit, their errors and indices (denoting processes are listed among the floating parameters)
    if ( (sPar1=="Floating")&&(sPar3=="InitialValue") ) {
      inLogFile.getline(logline,2000);//Skip a line
      for (int j=0; j<NVars; j++) {
	inLogFile.getline(logline,2000);
	istringstream strA(logline);
	strA >> sPar1 >> dPar1 >> (*ptFitPar)(j) >> sPar2 >> (*ptFitErr)(j); //dPar1 is the initial value
	// (*ptFitPar)(j) = dParV[1]; //dParV[0] is the initial value
	// (*ptFitErr)(j)=dParV[2];
	for (int i=0; i<NProc; i++) {
	  if ( sPar1==(ProcName[i]+"_nrm") ) { ProcIdx[i]=j; }
	}
      }
      fitPar=ptFitPar;
      fitErr=ptFitErr;
    }

    ///Get the chi2dof and chi2prob
    if ( (sPar1=="chi2:") ) {
      str >> sPar5 >> dPar1;
      chi2dof=dPar1;
      inLogFile.getline(logline,2000);
      istringstream strA(logline);
      strA >> sPar1 >> sPar2 >> dPar1;
      chi2prob=dPar1;
    }

  }

  ///Comptute the error on the total from the covariance matrix
  fitTotalErr=0.0;
  for (int i=0; i<NProc;i++) {
    for (int j=0; j<NProc;j++) {
      fitTotalErr = fitTotalErr + (*fitFullCov)(ProcIdx[i],ProcIdx[j])*theoryYield[i]*theoryYield[j];
    }
  }
  fitTotalErr=TMath::Sqrt(fitTotalErr);

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////---------- Derived Classes: Create one for each fit type ------------////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
class DibosonAntiBtagMuFit: public GeneralFitResult{
public:
  DibosonAntiBtagMuFit(TString inLogFileName):GeneralFitResult(3,7) {
    assignProcNames();
    readLogInfo(inLogFileName);
  }
  void assignProcNames();
};

void DibosonAntiBtagMuFit::assignProcNames() {
  ProcName.clear();
  ProcName.push_back("diboson");
  ProcName.push_back("WpJ");
  ProcName.push_back("top");
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
class DibosonAntiBtagElFit: public GeneralFitResult{
public:
  DibosonAntiBtagElFit(TString inLogFileName):GeneralFitResult(4,8) {
    assignProcNames();
    readLogInfo(inLogFileName);
  }
  void assignProcNames();
};

void DibosonAntiBtagElFit::assignProcNames() {
  ProcName.clear();
  ProcName.push_back("diboson");
  ProcName.push_back("WpJ");
  ProcName.push_back("top");
  ProcName.push_back("QCD");
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
class DibosonBtagFit: public GeneralFitResult{
public:
  DibosonBtagFit(TString inLogFileName):GeneralFitResult(4,7) {
    assignProcNames();
    readLogInfo(inLogFileName);
  }
  void assignProcNames();
};

void DibosonBtagFit::assignProcNames() {
  ProcName.clear();
  ProcName.push_back("diboson");
  ProcName.push_back("WpJ");
  ProcName.push_back("top");
  ProcName.push_back("WHbb");
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
class DibosonBoostedFit: public GeneralFitResult{
public:
  DibosonBoostedFit(TString inLogFileName):GeneralFitResult(3,6) {
    assignProcNames();
    readLogInfo(inLogFileName);
  }
  void assignProcNames();
};

void DibosonBoostedFit::assignProcNames() {
  ProcName.clear();
  ProcName.push_back("diboson");
  ProcName.push_back("WpJ");
  ProcName.push_back("top");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////---------- General Functions used to implement the toy fits and collect information from the logs ------------//////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
bool logfileContains(const char* logfileName, const char* stringOfInterest) 
//// Check that a particular string is present in the log file (or any other text file) 
{
  bool thefileContains=false;
  char logline[2000];
  ifstream inLogFile(logfileName);
  while ( inLogFile.good() ) {
    inLogFile.getline(logline,2000);
    TString str(logline);
    if ( str.Contains(stringOfInterest) ) {
      thefileContains=true;
    }
  }

  return thefileContains;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
GeneralFitResult* createFitResultInMemory(TString logFileName, int FitTypeFlag = -1) 
//// Create the appropriate derived class in memory and return a pointer to it.
//// The FitTypeFlag is used to indicate which fit type will be created: 1 - DibosonAntiBtagMu, 2 - DibosonAntiBtagEl, 3 - DibosonBtag, 4 - DibosonBoosted
{
  GeneralFitResult* ftResult;

  switch ( FitTypeFlag ) {
  case 1:
    ftResult = new DibosonAntiBtagMuFit(logFileName);
    break;
  case 2:
    ftResult = new DibosonAntiBtagElFit(logFileName);
    break;
  case 3:
    ftResult= new DibosonBtagFit(logFileName);
    break;
  case 4:
    ftResult = new DibosonBoostedFit(logFileName);
    break;
  default:
    cerr << "Error: Unable to create fit result class" << endl;
    break;
  }

  return ftResult;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void runToyMCFits(int NStart, int NEnd, int fitTypeFlag = -1, int InitRand = 2653, const char* inputLogName = "Fit_DibosonStandard_mu.txt", const char* CmdStringHead = "python -i runDiboson8TevFit.py -j 2 -m Diboson8TeVConfig --runPdfGenToySim --genParamFiles DibosonStandardFitOutPars_mu.txt --noNull --nosig", const char* CmdStringArgs = "topDibosonParameters.txt dibosonDibosonParameters.txt ZpJDibosonParameters.txt WpJDibosonParameters.txt", const char* logDir="./", const char* logPrefix="StandardFit_WpJDefault_Muon_", bool fullySmearAroundCenter = true, bool skipConfigurationsWithNegativeYields = true ) 
//// Construct a command string from input commands and generated parameters and implement the fit.
//// One time check the output and if covariance matrix "forced positive-definite" rerun the fit using the parameters from the first fit (TempFitPars.txt) as a starting point
{
  TString Command,Command_body,Command_logend,LogName,ConfigName,ParamfileName;
  char I_char[7];
  GeneralFitResult* fitResult = createFitResultInMemory(inputLogName,fitTypeFlag);

  for (int nFit=NStart; nFit<(NEnd+1); nFit++) {
    cout << "Running Toy Fit " << nFit << endl;
    ///Construct the command string
    //Log file:
    sprintf(I_char,"%i",nFit);
    LogName=".log";
    LogName=I_char+LogName;
    LogName=logPrefix+LogName;
    LogName=logDir+LogName;
    Command_logend=" > "+LogName;
    //Random Seeding:
    Command_body=" ";
    sprintf(I_char,"%i",InitRand+nFit*3+89);
    Command_body=I_char+Command_body;
    Command_body=" --seed "+Command_body;
    //Yields to use in generation:
    fitResult->computeYieldsToGen(InitRand+587+nFit*2,fullySmearAroundCenter,skipConfigurationsWithNegativeYields);

    for (int i=0; i<fitResult->NProc; i++) {
      sprintf(I_char,"%i",int(fitResult->yieldToGen[i]+0.5));
      Command_body=I_char+Command_body;
      Command_body=" --ext" + fitResult->ProcName[i] + " " + Command_body;//Note the flag with which the external yields should be passed to the .py file
    }
    // create a file with temporary fit output parameters
    Command_body=" --fitparfn TempFitPars.txt" + Command_body;

    //Combine to form the command string
    Command = CmdStringHead + Command_body + CmdStringArgs + Command_logend;
    cout << "Command=" << Command << endl;
    system(Command);

    //Check if the fit produced a positive-definite covariance matrix and if not rerun it using the output parameters as a starting point
    cout << "Checking the status of the covariance matrix ... ";
    if ( logfileContains(LogName,"covariance matrix quality: Full matrix, but forced positive-definite") ) {
      cout << "not positive-definite. Rerun the fit." << endl;
      Command = "mv " + LogName + " " + LogName + "_unconverged";
      cout << "Command=" << Command << endl;
      system(Command);
      Command = CmdStringHead  + Command_body + "TempFitPars.txt" + Command_logend;
      cout << "Command=" << Command << endl;
      system(Command);
    } else {
      cout << "is positive-definite" << endl;
    }

  }
  delete fitResult;
}


// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void processLogs(int fitTypeFlag = -1, int NLogs=1, const char* processlogPrefix="./temp/BtagFitMuon_ShapeBased_VariableYield_GenWpJDataDefaultFitDefault_", const char* outFileName="./test_LogSummary.root")
//// Processes the logs and produces an output file. Log files are of the form processlogPrefix*.log, where *=0-(NLogs-1)
{
  TString inLogFileName;
  char Lg_char[5];
  GeneralFitResult* fitResult;
  TString processName;

  TFile *outfile = new TFile(outFileName, "RECREATE");
  TTree *OutTree = new TTree("OutTree","OutTree");

  cout << "Initializing with Log 0" << endl;
  inLogFileName="0.log";
  inLogFileName=processlogPrefix+inLogFileName;
  fitResult = createFitResultInMemory(inLogFileName,fitTypeFlag);

  const int NProc = fitResult->NProc;
  double theoryValProcess[NProc], genProcess[NProc], numProcess[NProc], errProcess[NProc];

  for (int i=0; i<NProc; i++) {
    processName=fitResult->ProcName[i];
    OutTree->Branch("TheoryYield_"+processName,&(theoryValProcess[i]),"TheoryYield_"+processName+"/D");
    OutTree->Branch("GenYield_"+processName,&(genProcess[i]),"GenYield_"+processName+"/D");
    OutTree->Branch("FitYield_"+processName,&(numProcess[i]),"FitYield_"+processName+"/D");
    OutTree->Branch("FitErr_"+processName,&(errProcess[i]),"FitErr_"+processName+"/D");
  }
  OutTree->Branch("fitTotal",&(fitResult->fitTotal),"fitTotal/D");
  OutTree->Branch("fitTotalErr",&(fitResult->fitTotalErr),"fitTotalErr/D");
  OutTree->Branch("dataTotal",&(fitResult->fitTotal),"dataTotal/D");
  OutTree->Branch("chi2dof",&(fitResult->chi2dof),"chi2dof/D");
  OutTree->Branch("chi2prob",&(fitResult->chi2prob),"chi2prob/D");
  OutTree->Branch("Convergence",&(fitResult->Convergence),"Convergence/I");

  /// Process the logs
  for (int lg=0; lg<NLogs; lg++) {
    cout << "Processing Log " << lg << " : " << endl;
    if ( lg!=0 ) {//Don't reload the first log
      sprintf(Lg_char,"%i",lg);
      inLogFileName=".log";
      inLogFileName=Lg_char+inLogFileName;
      inLogFileName=processlogPrefix+inLogFileName;
      fitResult = createFitResultInMemory(inLogFileName,fitTypeFlag);
    }

    for (int i=0; i<NProc; i++) {
      theoryValProcess[i]=fitResult->theoryYield[i];
      genProcess[i] = fitResult->yieldFromGen[i];
      numProcess[i] = ((*fitResult->fitPar)(fitResult->ProcIdx[i]))*(theoryValProcess[i]);
      errProcess[i] = ((*fitResult->fitErr)(fitResult->ProcIdx[i]))*(theoryValProcess[i]);
      cout << "Recording " << fitResult->ProcName[i] << ": generated=" << genProcess[i] << ", fitted = " << numProcess[i] << " +/- " << errProcess[i] << endl;
    }
    OutTree->Fill();

    delete fitResult;
  }

  outfile->Write();
  outfile->Close();

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void printLogSummary(int fitTypeFlag = -1, const char* inLogFileName="./MyLog.log")
//// Print the summary of a particular fit log file
{
  GeneralFitResult* fitResult = createFitResultInMemory(inLogFileName,fitTypeFlag);
  fitResult->printSummary();
}
