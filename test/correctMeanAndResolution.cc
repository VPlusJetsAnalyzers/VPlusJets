/*******************************************************************
 * Project: CMS detector at the CERN, diboson analysis in the semileptonic channel
 *
 *
 * Authors:
 *
 *   Osipenkov, Ilya, Texas A&M - ilyao@fnal.gov
 *
 * Description: Use to correct mean and resolution parameters based on the control sample
 *
 ********************************************************************/

#include <iostream>
#include <fstream>
#include <strstream>
#include <vector>
#include <string>

#include <stdio.h>      /* printf, fgets */
#include <stdlib.h>     /* atof */

#include <TFile.h>
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

void correctMeanAndResolution(TString paramfilePrefix="dibosonDibosonParameters",TString paramfileSuffix="_mu.txt", int correctionCode = 0)
//// Takes the input parameter file (name=paramfilePrefix+parmfileSuffix) and creates an output parameter file (name=paramfilePrefix+"_SidebandModified"+corrSuffix+paramfileSuffix) with the corrected mean and resolution/sigma
{
  ifstream inparFile(paramfilePrefix+paramfileSuffix);
  TString outparFileName;
  ofstream outparFile;
  double meanCorr=-1.0;
  double sigmaCorr=0.0;

  TString sPar[9];
  char inputline[500];
  char tempVal_char[10];
  double val=-1;
  double uncorrval=-1;
  TString nameVar, meanVar, sigmaVar;
  TString corrSuffix="";
  int processCode=-1;
  switch ( correctionCode ) {
  case 1:
    //diboson central
    processCode = 0;
    meanCorr=0.93;
    sigmaCorr=1.091;
    break;
  case 2:
    //diboson +1 sigma
    corrSuffix="p1Sigma";
    processCode = 0;
    meanCorr=1.27;
    sigmaCorr=1.127;
    break;
  case 3:
    //diboson -1 sigma
    corrSuffix="m1Sigma";
    processCode = 0;
    meanCorr=0.59;
    sigmaCorr=1.055;
    break;
  case 11:
    //top central
    processCode = 1;
    meanCorr=0.93;
    sigmaCorr=1.091;
    break;
  case 12:
    //top +1 sigma
    corrSuffix="p1Sigma";
    processCode = 1;
    meanCorr=1.27;
    sigmaCorr=1.127;
    break;
  case 13:
    //top -1 sigma
    corrSuffix="m1Sigma";
    processCode = 1;
    meanCorr=0.59;
    sigmaCorr=1.055;
    break;
  default:
    cerr <<"Error, input the correct correctionCode" << endl;
    break;
  }

  switch ( processCode ) {
  case 0:
    nameVar="diboson";
    meanVar="mean_"+nameVar;
    sigmaVar="resolution_"+nameVar;
    break;
  case 1:
    nameVar="top";
    meanVar="mean_"+nameVar;
    sigmaVar="sigma_"+nameVar;
    break;
  default:
    cerr <<"Error, input the correct processCode" << endl;
    break;
  }


  outparFile.open(paramfilePrefix+"_SidebandModified"+corrSuffix+paramfileSuffix,ios::out);


  while ( inparFile.good() ) {
    inparFile.getline(inputline,500);
    istrstream str(inputline);
    for (int i=0; i<9;i++) { str >> sPar[i]; }

    if ( !sPar[0].IsNull() ) {
      //cout << sPar[0] << endl;
      val=atof(sPar[2]);

      if ( sPar[0].Contains(meanVar) ) {
	uncorrval = val;
	val=meanCorr+val;
	cout << "Changed " << sPar[0] << " from " << uncorrval << " to " << val << endl;
      }
      if ( sPar[0].Contains(sigmaVar) ) {
	uncorrval = val;
	val=sigmaCorr*val;
	cout << "Changed " << sPar[0] << " from " << uncorrval << " to " << val << endl;
      }

      sprintf(tempVal_char,"%f",val);
      sPar[2]=tempVal_char;

      for (int i=0; i<9;i++) { outparFile << sPar[i] << " " ; }
      outparFile << endl;
    }

  }


}



