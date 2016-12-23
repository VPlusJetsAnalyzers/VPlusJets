////////////////////////////////////////////////////////////////////////////////////////////////
////                    CMS - Semileptonic Channel (jj+l+nu final state)                    ////
////                                    Diboson Analysis                                    ////
////    Created by Osipenkov, Ilya (Texas A&M) : ilyao@fnal.gov, ilyao@physics.tamu.edu     ////
////////////////////////////////////////////////////////////////////////////////////////////////
////                        Unfolding Tools utilizing RooUnfold by                          ////
////       Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>       ////
////////////////////////////////////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iomanip>
#include <fstream>
#include <set>
#include <algorithm>

#include "TROOT.h"
#include "TEventList.h"

#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldSvd.h"
#include "src/RooUnfoldTUnfold.h"
#include "src/RooUnfoldBinByBin.h"
#include "src/TSVDUnfold_local.h"

#include "/uscms_data/d3/ilyao/CMSSW_5_3_15_patch1/src/ElectroWeakAnalysis/VPlusJets/test/PlotTools.cc"

#endif


////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////*********      Class Definitions     *********/////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
struct MCEvent{
  MCEvent(int evtNumberIn, double mcVarValIn): evtNumber(evtNumberIn), mcVarVal(mcVarValIn), isDetectorEvt(false) {};
  MCEvent(int evtNumberIn, double mcVarValIn, bool isDetectorEvtIn): evtNumber(evtNumberIn), mcVarVal(mcVarValIn), isDetectorEvt(isDetectorEvtIn) {};
  int evtNumber;
  double mcVarVal;
  double evtWt;
  bool isDetectorEvt;
  bool addAsMissing;
  bool addToHistMCTruth;
  ///MCEvents are used in sets: elements a and b are considered equal iff !(a < b) && !(b < a)
  bool operator<(const MCEvent &other) const {
    if ( evtNumber<other.evtNumber ) { return true; }
    if ( evtNumber==other.evtNumber ) {
      if ( (mcVarVal<other.mcVarVal)&&(abs((mcVarVal-other.mcVarVal)/(mcVarVal+other.mcVarVal))>0.0001) ) {
	return true;
      }
    }
    return false;
  }
};


int getNegExponent(double numIn){
  numIn=abs(numIn);
  if ( numIn >=1 ) {
    cout << "Error: can not compute negative exponent, absolute value of the input is not less than 1" << endl;
    return -1;
  }
  int exp=0;
  int expMax=30;
  while ( (numIn<1)&&(exp<expMax) ) {
    exp++;
    numIn*=10;
  }
  if ( exp==expMax ) { cout << "Warning: maximum allowed exponent computed. Input is likely 0" << endl; }
  return exp;
}



///////////////////////////////           Main Class:            ///////////////////////////////
class UnfoldDibosonInPt: public SingleProcess<processes> {
/// Inputs: a) Diboson process, either explicitly as a single process or extracted from a SimulationStack, b) Data-MC object containing histograms and uncertainties.
/// Build In: MC Selection - variable, cuts, binning.
public:
  UnfoldDibosonInPt(const SingleProcess<processes>& dibosonProcessIn, bool isSimulationTestIn, CanvasProperties cnvPropertiesIn, TH1F* histToUnfoldIn = NULL): SingleProcess<processes>(dibosonProcessIn), isSimulationTest(isSimulationTestIn), cnvProperties(cnvPropertiesIn) {
    setAllHistPointersToNULL();
    setVariableProperties();
    if ( !isSimulationTest ) {
      if ( histToUnfoldIn!=NULL ) {
	histToUnfold=new TH1F(*histToUnfoldIn);
      } else {
	cout << "Error: no histogram to unfold is specified nor is the unfolding is not based on simulation" << endl;
      }
    }
    fillResponseAndHists();
  }
  ~UnfoldDibosonInPt() { 
//     if ( histToUnfold!=NULL ) { 
//       delete histToUnfold;
//     }
//     if ( histMCTruth!=NULL ) { 
//       delete histMCTruth;
//     }
//     if ( histReco!=NULL ) { 
//       delete histReco;
//     }
//     if ( histReco_ErrUpper!=NULL ) { 
//       delete histReco_ErrUpper;
//     }
//     if ( histReco_ErrLower!=NULL ) { 
//       delete histReco_ErrLower;
//     }
  }

  ///Various tools
  void setAllHistPointersToNULL () {
    histToUnfold=NULL;
    histMCTruth=NULL;
    histReco=NULL;
    histReco_ErrUpper=NULL;
    histReco_ErrLower=NULL;
  }
  void setVariableProperties();
  void fillResponseAndHists();
  void printRawResponseProperties();
  void printResponsePropertiesToTableFile(const char*);

  void runSVDUnfolding();
  void runBayesUnfolding();
  void augmentHistRecoUncertaintiesBasedOnAlternateMCTableInFile(const char*);
  void printRecoHistToRawTableFile(const char*);
  void printRecoHistToFormattedTableFile(const char*);
  void fillRecoErrorHistograms();
  bool nullHistograms() {
    if ( (histToUnfold==NULL)||(histMCTruth==NULL)||(histReco==NULL)||(histReco_ErrUpper==NULL)||(histReco_ErrLower==NULL) ) { return true; }
    return false;
  }

  void formatAllHistograms();
  void rescaleAllHistogramsByLumi(double=19250.0);
  void printHistogramEntriesAndIntegrals();
  void drawHistogramsSeparatelyAndSaveWithPrefix(TString);
  TCanvas* drawOverlayedHistsWithTitle(TString);
  void placeCMSLumiFrame(TPad* pad, int iPeriod=3, int iPosX=10) { cnvProperties.setVarsAndRunCMSLumi(pad,iPeriod,iPosX); }

private:
  TString getCuts_hadlep(TString inFilePath) {
    if ( inFilePath.Contains("_mu") ) { 
      return cuts_MCTruthhadlep_mu;
    } 
    return cuts_MCTruthhadlep_el;
  }
  TString getCuts_lephad(TString inFilePath) {
    if ( inFilePath.Contains("_mu") ) { 
      return cuts_MCTruthlephad_mu;
    }
    return cuts_MCTruthlephad_el;
  }

  const bool isSimulationTest;//Fill the response based on first half of events and use to test on the second half of simulation if true
  static TString cuts_MCTruthhadlep_base, cuts_MCTruthlephad_base, cuts_MCTruthhadlep_mu, cuts_MCTruthlephad_mu, cuts_MCTruthhadlep_el, cuts_MCTruthlephad_el;

  TH1F* histToUnfold;
  TH1F* histMCTruth;
  RooUnfoldResponse response;
//   double nDetectorTotal;
//   double nFakeTotal;
//   double nMissingTotal;
  TH1F* histReco;
  TH1F *histReco_ErrUpper, *histReco_ErrLower;

  /// Binning for detector level and MCTruth variables (pTs)
  int nBinsDetectorVar;
  double minDetectorVar;
  double maxDetectorVar;
  int nBinsTrueVar; 
  double minTrueVar;
  double maxTrueVar;

  /// Plot Formatting
  CanvasProperties cnvProperties;

};

// TString UnfoldDibosonInPt::cuts_MCTruthhadlep_base = "(W_V_hadronic_gen[0]==1)&&(W_V_hadronic_gen[1]==0)&&(W_V_pt_gen[0]>70.0)&&(W_V_pt_gen[0]<240.0)&&(W_V_dau0_pt_gen[0]>30)&&(W_V_dau1_pt_gen[0]>30)&&(abs(W_V_dau0_eta_gen[0])<2.4)&&(abs(W_V_dau1_eta_gen[0])<2.4)&&(sqrt(W_V_et_gen[1]*W_V_et_gen[1]-W_V_pt_gen[1]*W_V_pt_gen[1])>25)&&(W_neutrino_pt_gen>20)";
// TString UnfoldDibosonInPt::cuts_MCTruthhadlep_mu = UnfoldDibosonInPt::cuts_MCTruthhadlep_base + "&&(W_muon_pt_gen>20)&&(abs(W_muon_eta_gen)<2.1)";
// TString UnfoldDibosonInPt::cuts_MCTruthhadlep_el = UnfoldDibosonInPt::cuts_MCTruthhadlep_base + "&&(W_electron_pt_gen>25)&&(abs(W_electron_eta_gen)<2.5)";

// TString UnfoldDibosonInPt::cuts_MCTruthlephad_base = "(W_V_hadronic_gen[0]==0)&&(W_V_hadronic_gen[1]==1)&&(W_V_pt_gen[1]>70.0)&&(W_V_pt_gen[1]<240.0)&&(W_V_dau0_pt_gen[1]>30)&&(W_V_dau1_pt_gen[1]>30)&&(abs(W_V_dau0_eta_gen[1])<2.4)&&(abs(W_V_dau1_eta_gen[1])<2.4)&&(sqrt(W_V_et_gen[0]*W_V_et_gen[0]-W_V_pt_gen[0]*W_V_pt_gen[0])>25)&&(W_neutrino_pt_gen>20)";
// TString UnfoldDibosonInPt::cuts_MCTruthlephad_mu = UnfoldDibosonInPt::cuts_MCTruthlephad_base + "&&(W_muon_pt_gen>20)&&(abs(W_muon_eta_gen)<2.1)";
// TString UnfoldDibosonInPt::cuts_MCTruthlephad_el = UnfoldDibosonInPt::cuts_MCTruthlephad_base + "&&(W_electron_pt_gen>25)&&(abs(W_electron_eta_gen)<2.5)";

TString UnfoldDibosonInPt::cuts_MCTruthhadlep_base = "(W_V_hadronic_gen[0]==1)&&(W_V_hadronic_gen[1]==0)&&(W_V_pt_gen[0]>70.0)&&(W_V_pt_gen[0]<240.0)";
TString UnfoldDibosonInPt::cuts_MCTruthhadlep_mu = UnfoldDibosonInPt::cuts_MCTruthhadlep_base + "&&(W_muon_pt_gen>20)&&(abs(W_muon_eta_gen)<2.1)";
TString UnfoldDibosonInPt::cuts_MCTruthhadlep_el = UnfoldDibosonInPt::cuts_MCTruthhadlep_base + "&&(W_electron_pt_gen>25)&&(abs(W_electron_eta_gen)<2.5)";

TString UnfoldDibosonInPt::cuts_MCTruthlephad_base = "(W_V_hadronic_gen[0]==0)&&(W_V_hadronic_gen[1]==1)&&(W_V_pt_gen[1]>70.0)&&(W_V_pt_gen[1]<240.0)";
TString UnfoldDibosonInPt::cuts_MCTruthlephad_mu = UnfoldDibosonInPt::cuts_MCTruthlephad_base + "&&(W_muon_pt_gen>20)&&(abs(W_muon_eta_gen)<2.1)";
TString UnfoldDibosonInPt::cuts_MCTruthlephad_el = UnfoldDibosonInPt::cuts_MCTruthlephad_base + "&&(W_electron_pt_gen>25)&&(abs(W_electron_eta_gen)<2.5)";


void UnfoldDibosonInPt::setVariableProperties() {
  ///Get the binning based on the diboson process
  SingleFileInput firstDibosonFile = allProcessFiles[0];
  AnalysisVariable detectorVar = firstDibosonFile.getAnalysisVariable();
  nBinsDetectorVar = detectorVar.getNumberOfBins();
  minDetectorVar = detectorVar.getMin();
  maxDetectorVar = detectorVar.getMax();
  //Unfold to the same binning as the detector pT
  nBinsTrueVar = nBinsDetectorVar; 
  minTrueVar = minDetectorVar;
  maxTrueVar = maxDetectorVar;
}

void UnfoldDibosonInPt::fillResponseAndHists() {
  /// If isSimulationTest==true fill the response matrix based on the first half of MC, while histMCTruth and histToUnfold are filled based on the second half of MC (file by file)
  /// If isSimulationTest==false fill the response matrix and histMCTruth based on all of MC, while histToUnfold should be provided externally
  if ( isEmpty() ) {
    cout << "Unable to make the Response Matrices, no diboson files provided" << endl;
    return;
  }

  cout << "Initializing histMCTruth" << endl;
  histMCTruth = new TH1F("histMCTruth","histMCTruth",nBinsTrueVar,minTrueVar,maxTrueVar);
  histMCTruth->Sumw2();
  if ( isSimulationTest ) {
    cout << "Initializing histToUnfold" << endl;
    histToUnfold = new TH1F("histToUnfold","histToUnfold",nBinsDetectorVar,minDetectorVar,maxDetectorVar);
    histToUnfold->Sumw2();
  }

  cout << "Initializing the Response Matrix" << endl;
  response = RooUnfoldResponse(nBinsDetectorVar,minDetectorVar,maxDetectorVar,nBinsTrueVar,minTrueVar,maxTrueVar);

  ///Iterate over the diboson files and record the relevant info
  TFile* currentFile;
  TTree* currentTree;
  TString filePath, treeName, detectorCuts;
  TString cuts_MCTruthhadlep, cuts_MCTruthlephad;
  int nEntries = -1;
  double externalWeight = -1.0;
  TEventList *list_All, *list_Detector, *list_MCTruthhadlep, *list_MCTruthlephad;
  Long_t ent;
  double evtWt=-1.0;
  int evtNumber;
  float effwt, puwt;
  float JetPFCor_Pt[6];
  float JetPFCor_Phi[6];
  float W_V_pt_gen[6];
  bool detectorEvt=false, mcEvt=false;
  float detectorVarVal=-1.0, mcVarVal=-1.0;

  cout << "Filling the Response Matrix and histograms ..." << endl;
  int nFillTotal=0;///*Unweighted* event counts
  int nFakeTotal=0;
  int nMissingTotal=0;
  int nDuplicateMissingEvents = 0;
  std::set<MCEvent> mcEvents;
  //map<int,MissingEvent> missingEvents; //There should not be any overlap between files when it comes to events passing detector cuts; however the same events may be considered missing for multiple files
  for ( vector<SingleFileInput>::iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    SingleFileInput currentDibosonFile = *fileIter;
    filePath = currentDibosonFile.getFilePath();
    currentFile = new TFile(filePath);
    treeName = currentDibosonFile.getTreeName();
    currentTree = (TTree*)currentFile->Get(treeName);
    nEntries=currentTree->GetEntries();
    detectorCuts=currentDibosonFile.produceCutStringWithoutWeight();
    externalWeight=currentDibosonFile.getKFactorXLumiDivByNumGenerated();

    cout << "Processing file in path: " << filePath << " with " << nEntries << " entries" << endl;
    /// Store All Events
    currentTree->Draw(">>templist_All");
    list_All = (TEventList*)gDirectory->Get("templist_All");
    nEntries=list_All->GetN();
    cout << nEntries << " total entries in the file; externalWeight = " << externalWeight << endl;
    /// Record Detector Events
    currentTree->Draw(">>templist_Detector",detectorCuts);
    list_Detector = (TEventList*)gDirectory->Get("templist_Detector");
    cout << list_Detector->GetN() << " entries passing default (detector level) cuts" << endl;
    /// Record Events Passing MC Truth Cuts
    cuts_MCTruthhadlep=getCuts_hadlep(filePath);
    cuts_MCTruthlephad=getCuts_lephad(filePath);
    currentTree->Draw(">>templist_MCTruthhadlep",cuts_MCTruthhadlep);
    list_MCTruthhadlep = (TEventList*)gDirectory->Get("templist_MCTruthhadlep");
    cout << list_MCTruthhadlep->GetN() << " entries passing MCTruth hadlep cuts" << endl;
    currentTree->Draw(">>templist_MCTruthlephad",cuts_MCTruthlephad);
    list_MCTruthlephad = (TEventList*)gDirectory->Get("templist_MCTruthlephad");
    cout << list_MCTruthlephad->GetN() << " entries passing MCTruth lephad cuts" << endl;

    //Set the branches
    currentTree->SetBranchStatus("*",0); //disable all branches
    currentTree->SetBranchStatus("event_evtNo",1);
    currentTree->SetBranchAddress("event_evtNo",&evtNumber);
    currentTree->SetBranchStatus("effwt",1);
    currentTree->SetBranchAddress("effwt",&effwt);
    currentTree->SetBranchStatus("puwt",1);
    currentTree->SetBranchAddress("puwt",&puwt);
    currentTree->SetBranchStatus("JetPFCor_Pt",1);
    currentTree->SetBranchAddress("JetPFCor_Pt",&JetPFCor_Pt);
    currentTree->SetBranchStatus("JetPFCor_Phi",1);
    currentTree->SetBranchAddress("JetPFCor_Phi",&JetPFCor_Phi);
    currentTree->SetBranchStatus("W_V_pt_gen[6]",1);
    currentTree->SetBranchAddress("W_V_pt_gen[6]",&W_V_pt_gen);

    /// Go through the events and either fill a response or record the fakes.
    int nFill = 0;
    int nFake = 0;
    int nMissing = 0;
    for (int i=0; i<nEntries ; i++) {
      detectorEvt=false;
      mcEvt=false;
      detectorVarVal=-1.0;
      mcVarVal=-1.0;
      ent = list_All->GetEntry(i);
      currentTree->GetEntry(ent);
      evtWt=externalWeight*effwt*puwt;

      if ( list_Detector->Contains(ent) ) {
	detectorEvt=true;
	detectorVarVal=sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]));
      }
      if ( list_MCTruthhadlep->Contains(ent) ) {
	mcEvt=true;
	mcVarVal=W_V_pt_gen[0];
      }
      if ( list_MCTruthlephad->Contains(ent) ) {
	mcEvt=true;
	mcVarVal=W_V_pt_gen[1];
      }

      /// Fill the responses
      if ( detectorEvt&&mcEvt ) {
	MCEvent evt(evtNumber,mcVarVal,true);
	if ( mcEvents.count(evt) ) {
	  mcEvents.erase(evt);//delete the missing event (which should have isDetectorEvt=false)
	  nMissing--;//what we previously thought was a missing event is actually a fill event in this channel
	  nDuplicateMissingEvents++;
	}
	mcEvents.insert(evt);
	if ( isSimulationTest ) {
	  if ( (i<=(nEntries/2)) ) {//***Comment out for closure test
	    response.Fill(detectorVarVal,mcVarVal,evtWt);
	    nFill++;
	  } else {
	    histMCTruth->Fill(mcVarVal,evtWt);
	    histToUnfold->Fill(detectorVarVal,evtWt);
	  }
	} else {
	  response.Fill(detectorVarVal,mcVarVal,evtWt);
	  nFill++;
	  histMCTruth->Fill(mcVarVal,evtWt);
	}
      }

      if ( detectorEvt&&(!mcEvt) ) {
	if ( isSimulationTest ) {
	  if ( (i<=(nEntries/2)) ) {//***Comment out for closure test
	    response.Fake(detectorVarVal,evtWt);
	    nFake++;
	  } else {
	    histToUnfold->Fill(detectorVarVal,evtWt);
	  }
	} else {
	  response.Fake(detectorVarVal,evtWt);
	  nFake++;
	}
      }

      if ( (!detectorEvt)&&mcEvt ) {
	MCEvent evt(evtNumber,mcVarVal);
	if ( mcEvents.count(evt) ) {
	  nDuplicateMissingEvents++;
	} else {
	  nMissing++;
	  evt.evtWt=evtWt;
	  if ( isSimulationTest ) {
	    if ( (i<=(nEntries/2)) ) {//***Comment out for closure test
	      evt.addAsMissing=true;
	      evt.addToHistMCTruth=false;
	    } else {
	      evt.addAsMissing=false;
	      evt.addToHistMCTruth=true;
	    }
	  } else {
	    evt.addAsMissing=true;
	    evt.addToHistMCTruth=true;
	  }
	  mcEvents.insert(evt);
	}
      }
    }

    nFillTotal+=nFill;
    nFakeTotal+=nFake;
    nMissingTotal+=nMissing;
    cout << "Added " << nFill << " Fill, " << nFake << " Fake and " << nMissing << " Missing events (to be added at the end)" << endl;
    cout << "Total " << nFillTotal << " FillTotal, " << nFakeTotal << " FakeTotal and " << nMissingTotal << " MissingTotal events at this point" << endl;
    cout << nDuplicateMissingEvents << " duplicates at this point (not added to the response)" << endl;
    cout << endl;
  }

  int cntColl=0;
  int cntAddedAsMiss=0;
  int cntAddedToHist=0;
  for ( set<MCEvent>::iterator setIter = mcEvents.begin(); setIter != mcEvents.end(); ++setIter ) {
    MCEvent evt=*setIter;
    cntColl++;
    if ( !evt.isDetectorEvt ) {
      if ( evt.addAsMissing ) { 
	response.Miss(evt.mcVarVal,evt.evtWt); 
	cntAddedAsMiss++;
      }
      if ( evt.addToHistMCTruth ) { 
	histMCTruth->Fill(evt.mcVarVal,evt.evtWt); 
	cntAddedToHist++;
      }
    }
  }

  cout << "MC Collection Summary: cntColl=" << cntColl << ", cntAddedAsMiss=" << cntAddedAsMiss << ", cntAddedToHist(in addition to Fills)=" << cntAddedToHist << endl;

}

void UnfoldDibosonInPt::printRawResponseProperties() {
  cout << "Print The Response Matrix (row,column)=(measured,truth)" << endl;
  response.Mresponse().Print("f=%8.4f");
  cout << endl;
  cout << "Print The Errors on the Unfolding Matrix (row,column)=(measured,truth)" << endl;
  response.Eresponse().Print("%6.4d");
  cout << endl;
  cout << "Print The Fakes" << endl;
  response.Vfakes().Print("f=%8.4f");
}

void UnfoldDibosonInPt::printResponsePropertiesToTableFile(const char* outTableFilePath) {
  char buffer[100];
  TString tableLine;
  cout << "Printing Response Matrix and Bin Summary to file: " << outTableFilePath << endl;
  ofstream outTableFile(outTableFilePath);

  outTableFile << "Response Matrix With Errors (In TableFormat)" << endl;
  double val=-1.0, err=-1.0;
  int exponent=-1;
  int nBinsTrue_Min=0, nBinsTrue_Max=0;
  const int nBinTrueColumnsToDisplay=5;
  while ( nBinsTrue_Max<nBinsTrueVar ) {
    nBinsTrue_Min=nBinsTrue_Max;
    nBinsTrue_Max=std::min(nBinsTrue_Max+nBinTrueColumnsToDisplay,nBinsTrueVar);
    outTableFile << "\\hline" << endl;
    tableLine="Bin$_{Detector}$";//16 characters
    for (int j=nBinsTrue_Min; j<nBinsTrue_Max;j++) { 
      sprintf(buffer," & Bin$_{MCTruth}$=%-13i",j);//Bin$_{MCTruth}$= -> 16characters
      tableLine+=buffer;
    }
    outTableFile << tableLine << " \\\\"<< endl;
    outTableFile << "\\hline" << endl;
    for (int i=0; i<nBinsDetectorVar;i++) {
      sprintf(buffer,"%-16i",i);
      tableLine=buffer;
      for (int j=nBinsTrue_Min; j<nBinsTrue_Max;j++) {
	val=response.Mresponse().GetMatrixArray()[i*nBinsTrueVar+j];
	err=response.Eresponse().GetMatrixArray()[i*nBinsTrueVar+j];
	if ( val<1e-10 ) {
	  if ( err>1e-10 ) { cout << "Warning: response value is near 0, but the error is not. Recoring both as 0..." << endl; }
	  sprintf(buffer," & %-29s","$0\\pm0$");
	} else {
	  exponent=getNegExponent(val);
	  val/=pow(10,-exponent);
	  err/=pow(10,-exponent);
	  sprintf(buffer," & $(%4.2f\\pm%-5.2f)\\times10^{-%i}$",val,err,exponent);//$\\times10^{-%i}$ -> 15characters
	}
	tableLine+=buffer;
      }
      outTableFile << tableLine << " \\\\"<< endl;
    }
    outTableFile << "\\hline" << endl;
  }
  outTableFile << "\v" << endl;
  outTableFile << endl;

  outTableFile << "Weighted Event Counts Bin By Bin as A Fraction of Total (based on MC input)" << endl;
  sprintf(buffer,"%-20s & %-19s & %-19s","Detector Bin","Frac$_{Detector}$","Frac$_{Fakes}$");
  outTableFile << buffer << " \\\\"<< endl;
  outTableFile << "\\hline" << endl;
  double nDetectorTotal=response.Vmeasured().Norm1();
  double fracDetector=-1.0, fracFakes=-1.0;
  double fracDetectorTot=0.0, fracFakesTot=0.0;
  int expDetector=0, expFakes=0;
  for (int i=0; i<nBinsDetectorVar;i++) {
    fracDetector=response.Vmeasured()[i]/nDetectorTotal;
    fracDetectorTot+=fracDetector;
    expDetector=getNegExponent(fracDetector);
    fracDetector/=pow(10,-expDetector);
    fracFakes=response.Vfakes()[i]/nDetectorTotal;
    fracFakesTot+=fracFakes;
    expFakes=getNegExponent(fracFakes);
    fracFakes/=pow(10,-expFakes);
    sprintf(buffer,"%-20i & $%4.2f\\times10^{-%i}$ & $%4.2f\\times10^{-%i}$",i,fracDetector,expDetector,fracFakes,expFakes);
    outTableFile << buffer << " \\\\" << endl;
  }
  outTableFile << "\\hline" << endl;
  sprintf(buffer,"%-20s & $%-17.2f$ & $%-17.2f$","Total",fracDetectorTot,fracFakesTot);
  outTableFile << buffer << " \\\\" << endl;
  outTableFile << "\\hline" << endl;
  outTableFile.close();
  cout << "... done" << endl;
}

void UnfoldDibosonInPt::runSVDUnfolding() {
  cout << endl;
  cout << "Performing Svd Unfolding" << endl;
  RooUnfoldSvd unfoldSvd (&response,histToUnfold,4);//default=4
  histReco = (TH1F*) unfoldSvd.Hreco();
}

void UnfoldDibosonInPt::runBayesUnfolding() {
  cout << endl;
  cout << "Performing Bayes Unfolding" << endl;
  RooUnfoldBayes   unfoldBayes (&response,histToUnfold,4,false);//default=4
  histReco = (TH1F*) unfoldBayes.Hreco();
}

void UnfoldDibosonInPt::augmentHistRecoUncertaintiesBasedOnAlternateMCTableInFile(const char* inTableFilePath) {
  /// The difference between the bin values and values in the table is added in quadrature to the current error
  /// The table lines should be of the form "BinMin-BinMax & Value ..."
  ifstream inTableFile(inTableFilePath);
  if ( !inTableFile.good() ) {
    cout << "Error: unable to augment the histogram uncertainties, no in file table at: " << inTableFilePath << endl;
    return;
  }
  char linein[300];
  int nBins = histReco->GetNbinsX();
  double recoVal = -1.0, recoErr = -1.0, alternateMCVal = -1.0, alternateMCErr = -1.0;
  TString temp1, temp2;
  inTableFile.getline(linein,300); //skip the first line
  for (int n=1; n<(nBins+1);n++) {
    recoVal=histReco->GetBinContent(n);
    recoErr=histReco->GetBinError(n);
    inTableFile.getline(linein,300);
    istringstream tableLine(linein);
    tableLine >> temp1 >> temp2 >> alternateMCVal;
    alternateMCErr = recoVal - alternateMCVal;
    histReco->SetBinError(n,sqrt(recoErr*recoErr+alternateMCErr*alternateMCErr));
  }
}

void UnfoldDibosonInPt::printRecoHistToRawTableFile(const char* outTableFilePath) {
  /// The table lines will be of the form used as an input to augmentHistRecoUncertaintiesBasedOnAlternateMCTableInFile: "BinMin-BinMax & Value $\\pm$ Error \\"
  cout << "Printing histReco to RawTableFile: " << outTableFilePath << endl;
  ofstream outTableFile(outTableFilePath);
  int nBins = histReco->GetNbinsX();
  double binXMin = -1.0, binXWidth = -1.0, binYVal = -1.0, binYErr = -1.0;

  outTableFile << "Table Listing bin (center) locations, along with corresponding values and errors" << std::endl;
  for (int n=1; n<(nBins+1);n++) {
    binXMin=histReco->GetBinLowEdge(n);
    binXWidth=histReco->GetBinWidth(n);
    binYVal=histReco->GetBinContent(n);
    binYErr=histReco->GetBinError(n);
    outTableFile << std::setprecision(0) << fixed << binXMin << "-" << binXMin+binXWidth << " & " << std::setprecision(2) << scientific << binYVal << " $\\pm$ " << binYErr << "  \\\\";
    outTableFile << std::endl;
  }
  outTableFile.close();
  cout << "...done" << endl;
}

void UnfoldDibosonInPt::printRecoHistToFormattedTableFile(const char* outTableFilePath) {
  /// The table lines will be of the form used in the documentation: "BinMin-BinMax & $(Value\\pmError)\times10^{-exp}$ \\"
  cout << "Printing histReco to formatted table file: " << outTableFilePath << endl;
  char buffer[100];
  //TString tableLine;
  ofstream outTableFile(outTableFilePath);
  int nBins = histReco->GetNbinsX();
  double binXMin = -1.0, binXWidth = -1.0, binYVal = -1.0, binYErr = -1.0;
  int binYExp=0;
  outTableFile << "Table Listing bin (center) locations, along with corresponding values and errors" << std::endl;
  outTableFile << "\\hline" << endl;
  outTableFile << "\\hline" << endl;
  sprintf(buffer,"%-20s & %-26s","V$_{had}$ \\pt (GeV)","Cross Section (pb)");
  outTableFile << buffer << " \\\\" << endl;
  outTableFile << "\\hline" << endl;
  for (int n=1; n<(nBins+1);n++) {
    binXMin=histReco->GetBinLowEdge(n);
    binXWidth=histReco->GetBinWidth(n);
    binYVal=histReco->GetBinContent(n);
    binYErr=histReco->GetBinError(n);
    if ( binYVal>1.0 ) {
      sprintf(buffer,"%3.0f-%-16.0f & $%4.2f\\pm%-19.2f$",binXMin,binXMin+binXWidth,binYVal,binYErr);
      outTableFile << buffer << " \\\\" << endl;
      continue;
    }
    binYExp=getNegExponent(binYVal);
    binYVal/=pow(10,-binYExp);
    binYErr/=pow(10,-binYExp);
    sprintf(buffer,"%3.0f-%-16.0f & $(%4.2f\\pm%4.2f)\\times10^{-%i}$",binXMin,binXMin+binXWidth,binYVal,binYErr,binYExp);//$\\times10^{-%i}$ -> 15characters
    outTableFile << buffer << " \\\\" << endl;
  }
  outTableFile << "\\hline" << endl;
  outTableFile.close();
  cout << "...done" << endl;
}


void UnfoldDibosonInPt::fillRecoErrorHistograms() {
  if ( histReco==NULL ) {
    cout << "Error: unable to construct error histograms. The Reco hisogram is NULL." << endl;
    return;
  }
  histReco_ErrUpper = new TH1F(*histReco);
  histReco_ErrUpper->SetTitle("Upper Error Envelope");
  histReco_ErrLower = new TH1F(*histReco);
  histReco_ErrLower->SetTitle("Upper Lower Envelope");

  int nBins = histReco->GetNbinsX();
  double recoVal = -1.0, recoErr = -1.0;
  for (Int_t n=1; n<(nBins+1);n++) {
    recoVal=histReco->GetBinContent(n);
    recoErr=histReco->GetBinError(n);
    histReco_ErrUpper->SetBinContent(n,recoVal+recoErr);
    histReco_ErrLower->SetBinContent(n,recoVal-recoErr);
  }
}

void UnfoldDibosonInPt::formatAllHistograms() {
  if ( nullHistograms() ) { 
    cout << "Error: unable to format all histograms. At least one of the histograms is NULL." << endl;
    return;
  }

  //histToUnfold->Sumw2();
  histToUnfold->SetMarkerColor(kBlue);
  histToUnfold->SetMarkerStyle(8);
  histToUnfold->SetLineColor(kBlue);
  histToUnfold->SetLineStyle(1);
  histToUnfold->SetLineWidth(2);

  //histMCTruth->Sumw2();
  histMCTruth->SetMarkerStyle(8);
  histMCTruth->SetLineStyle(7);

  //histReco->Sumw2();
  histReco->SetMarkerColor(kRed);
  histReco->SetMarkerStyle(8);
  histReco->SetLineColor(kRed);
  histReco->SetLineWidth(2);
  histReco->SetMaximum(1.2*histReco_ErrUpper->GetMaximum());
  histReco->SetMinimum(0.0);

  //histReco_ErrUpper->Sumw2();
  //histReco_ErrLower->Sumw2();
  histReco_ErrUpper->SetLineStyle(2);
  histReco_ErrUpper->SetLineColor(kRed);
  histReco_ErrLower->SetLineStyle(2);
  histReco_ErrLower->SetLineColor(kRed);
  histReco_ErrUpper->SetFillColor(kRed);
  histReco_ErrUpper->SetFillStyle(3253);
  histReco_ErrLower->SetFillColor(926);
}

void UnfoldDibosonInPt::rescaleAllHistogramsByLumi(double lumi) {
  if ( nullHistograms() ) { 
    cout << "Error: unable to rescale all histograms by lumi. At least one of the histograms is NULL." << endl;
    return;
  }
  histToUnfold->Scale(1.0/lumi);
  histMCTruth->Scale(1.0/lumi);
  histReco->Scale(1.0/lumi);
  histReco_ErrUpper->Scale(1.0/lumi);
  histReco_ErrLower->Scale(1.0/lumi);
}

void UnfoldDibosonInPt::printHistogramEntriesAndIntegrals() {
  if ( nullHistograms() ) { 
    cout << "Error: unable to print histogram entries and integrals. At least one of the histograms is NULL." << endl;
    return;
  }
  cout << "histToUnfold : nEntries=" << histToUnfold->GetEntries() << ", Integral=" << histToUnfold->Integral() << endl;
  cout << "histMCTruth : nEntries=" << histMCTruth->GetEntries() << ", Integral=" << histMCTruth->Integral() << endl;
  cout << "histReco : nEntries=" << histReco->GetEntries() << ", Integral=" << histReco->Integral() << endl;
  cout << "histReco_ErrUpper : nEntries=" << histReco_ErrUpper->GetEntries() << ", Integral=" << histReco_ErrUpper->Integral() << endl;
  cout << "histReco_ErrLower : nEntries=" << histReco_ErrLower->GetEntries() << ", Integral=" << histReco_ErrLower->Integral() << endl;
}

void UnfoldDibosonInPt::drawHistogramsSeparatelyAndSaveWithPrefix(TString savePrefix) {
  if ( nullHistograms() ) { 
    cout << "Error: unable to draw the histograms separately. At least one of the histograms is NULL." << endl;
    return;
  }
  TCanvas* cnvToUnfold = cnvProperties.makeCnv("cnvToUnfold");
  cnvToUnfold->cd();
  histToUnfold->SetTitle("Histogram To Unfold");
  histToUnfold->Draw("Ep");
  cnvToUnfold->Update();
  saveCanvasAs(cnvToUnfold,savePrefix+"toUnfold.png"); 
  TCanvas* cnvMCTruth = cnvProperties.makeCnv("cnvMCTruth");
  cnvMCTruth->cd();
  histMCTruth->SetTitle("True Value");
  histMCTruth->Draw("Ep");
  cnvMCTruth->Update();
  saveCanvasAs(cnvMCTruth,savePrefix+"MCTruth.png"); 
  TCanvas* cnvUnfolded = cnvProperties.makeCnv("cnvUnfolded");
  cnvUnfolded->cd();
  histReco->SetTitle("Unfolded Events");
  histReco->Draw("hist");
  histReco_ErrUpper->Draw("same hist");
  histReco_ErrLower->Draw("same hist");
  histReco->Draw("same hist");
  histReco->Draw("axis same");
  cnvUnfolded->Update();
  saveCanvasAs(cnvUnfolded,savePrefix+"Unfolded.png"); 
}

TCanvas* UnfoldDibosonInPt::drawOverlayedHistsWithTitle(TString title) {
  TCanvas* cnvTemp = cnvProperties.makeCnv(title);
  TLegend *lgndTemp = new TLegend(0.7,0.65,0.90,0.9);
  cnvTemp->cd();
  PlotProperties unfoldedPlotProperties(title,0.25,0.5,0.97,"dijet/V_{had} candidate p_{T}",0.9,"#sigma_{Fiducial} pb",1.2);
  histReco = unfoldedPlotProperties.labelHistogram(histReco);
  histReco->Draw("hist");
  lgndTemp->AddEntry(histReco,"Unfolded Distribution","l");
  histReco_ErrUpper->Draw("same hist");
  histReco_ErrLower->Draw("same hist");
  lgndTemp->AddEntry(histReco_ErrUpper,"Unfolded Error","f");
  histToUnfold->Draw("Ep:same");
  lgndTemp->AddEntry(histToUnfold,"Original Distribution","l");
  if ( isSimulationTest ) {
    histMCTruth->SetLineStyle(1);
    histMCTruth->Draw("Ep:same");
    lgndTemp->AddEntry(histMCTruth,"True Distribution","p");
  } else {
    histMCTruth->Draw("same hist");
    lgndTemp->AddEntry(histMCTruth,"MCTrue","l");
  }
  histReco->Draw("same hist");
  histReco->Draw("axis same");
  lgndTemp->Draw();
  cnvTemp->Update();
  return cnvTemp;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////********* Main Implementation*********////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void UnfoldDiboson(bool performSVDUnfolding=true, bool performBayesUnfolding=false)
{
  /// Setup:
#ifdef __CINT__
  gSystem->Load("libRooUnfold");
#endif
  gStyle->SetOptStat(0);
  gROOT->ProcessLine(".L tdrstyle.C");
  initializeCombinedInputs();

  TString outputSuffix="";
  bool displayHistograms=true;
  bool saveHistograms=true;
  TString plotDir = "plots/";
  if ( saveHistograms ) { displayHistograms=true; }

  CanvasProperties canvas900x600("cnv900x600",10,10,900,600);
  AnalysisVariable pTjj("sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",70,240,17);//Default
  //AnalysisVariable pTjj("sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",70,240,10);

  PlotProperties pTStackPlotProperties("StackPlot: Dijet p_{T} - All Resolved Channels",0.25,0.5,0.97,"dijet/V_{had} candidate p_{T}",1.0,"Event Count - Normalized to Fit Output",1.2);
  StackPlotProperties<processes> defaultChanHistProperties(canvas900x600,pTStackPlotProperties,make_pair(0.7,0.6),make_pair(0.90,0.90),legendColor,orderedNameColorPairs);

  SingleProcess<processes> processCollection[nResolvedChannels][nMCProcesses];
  vector < SingleProcess<processes> > dataByChannel(nResolvedChannels, SingleProcess<processes>(data));


  /// 1. Load the process info
  cout << "Loading process info ..." << endl;
  map<processes,vector< pair<TString,double> > > currentChannel;
  vector< pair<TString,double> > currentProcess;
  bool emptyProcess=false;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    cout << "Loading channel " << chanTitleLabel[chan] << " channel" << endl;
    currentChannel = simFilesAndCoeffs[chan];
    for (int prc=0; prc<nMCProcesses;prc++) { 
      processes proc = (processes)prc;
      if ( currentChannel.count(proc) ) {
	emptyProcess=false;
	currentProcess = currentChannel[proc];
	for ( vector< pair<TString,double> >::iterator fileIter = currentProcess.begin(); fileIter != currentProcess.end(); ++fileIter ) {
	  if ( proc==qcd ) {
	    processCollection[chan][proc]+=SingleFileInput(fileIter->first,pTjj,fileIter->second,qcdCutString[chan],dataEventWeightString);
	  } else {
	    processCollection[chan][proc]+=SingleFileInput(fileIter->first,pTjj,fileIter->second,defaultCutString[chan],mcEventWeightString);
	  }
	}
	processCollection[chan][proc].setDesignation(proc);
	processCollection[chan][proc].remakeProcessAndFileHistograms();
      } else {
	emptyProcess=true;
      }
      cout << procLabel[proc] << " :" << endl;
      if ( emptyProcess ) {
	cout << "No files contributing to the process vs. " << expectedProcessYields[chan][proc] << " events expected for mjj cuts" << endl;
      } else {
	cout << processCollection[chan][proc].getProcessEvents() << " events expected based on contributing files vs. " << expectedProcessYields[chan][proc] << " events expected for mjj cuts" << endl;
      }
    }
    cout << endl;
  }
  cout << "Filling the data processes ..." << endl;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    dataByChannel[chan]+=SingleFileInput(dataFiles[chan],pTjj,1.0,defaultCutString[chan],dataEventWeightString);
  }
  cout << endl;

  /// 2. Reweigh to fit results with corrected frac used for the diboson; combine the individual processes; compare to data and reweigh VJets+QCD to match the data event count
  cout << "Reweighting to the fit results (diboson is bias corrected) ..." << endl;
  double scalingFactor=-1.0;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    for (int prc=0; prc<nMCProcesses;prc++) { 
      processes proc = (processes)prc;
      scalingFactor=fitFracVal[chan][proc];
      if ( proc==diboson ) { scalingFactor=signalCorrectedFitFracVal[chan]; }
      if ( !processCollection[chan][proc].isEmpty() ) { processCollection[chan][proc].rescaleProcessAndFilesBy(scalingFactor); }
    }
  }

  cout << "Combining data and MC for each channel ..." << endl;
  SimulationStack<processes> emptyChanStack(defaultChanHistProperties);
  vector < DataAndMCStack<processes> > dataPlusMC(nResolvedChannels,DataAndMCStack<processes>(emptyChanStack));
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    cout << chanTitleLabel[chan] << " channel" << endl;
    dataPlusMC[chan].setStackHistTitle("StackPlot: Dijet p_{T} - " + chanTitleLabel[chan]);
    dataPlusMC[chan].addDataProcess(dataByChannel[chan]);
    for (int prc=0; prc<nMCProcesses;prc++) { 
      processes proc = (processes)prc;
      if ( !processCollection[chan][proc].isEmpty() ) { 
	dataPlusMC[chan].addSimulationProcess(processCollection[chan][proc]); 
      }
    }
  }

  cout << "Rescaling V+Jets and QCD MC to match the event count to data..." << endl;
  vector<processes> backgroundsToScale;
  backgroundsToScale.clear();
  backgroundsToScale.push_back(vjets);
  backgroundsToScale.push_back(qcd);
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    dataPlusMC[chan].remakeIndividualFileProcessStackHistAndLgnd();
    dataPlusMC[chan].rescaleDesignatedProcessesToMatchMCToData(backgroundsToScale);
    dataPlusMC[chan].remakeIndividualFileProcessStackHistAndLgnd();
    cout << "Channel : " << chanTitleLabel[chan] << endl;
    cout << dataPlusMC[chan].eventsInStackHist() << " rescaled MC events vs " << dataPlusMC[chan].getDataEvents() << " Data events with " << expectedDataYield[chan] << " expected for mjj cuts"  << endl;
  }

  if ( displayHistograms ) {
    TCanvas* datavsmccnv[nResolvedChannels];
    for (int ch=0; ch<nResolvedChannels;ch++) {
      channels chan = (channels)ch;
      datavsmccnv[chan] = dataPlusMC[chan].drawDataAndMCOnCanvas("datavsmccnv"+chanPlotLabel[chan]);
      datavsmccnv[chan]->Update();
      if ( saveHistograms ) { saveCanvasAs(datavsmccnv[chan],plotDir + "dijetpT_SvdUnfolding_DataUnfold_" + chanPlotLabel[chan] + "_Stacked"+outputSuffix+".png"); }//keep the same prefix, regardless of the unfolding procedure, to avoid duplication
    }
  }
  cout << endl;


  /// 3. Make the subtracted histograms for each channel 
  cout << "Making the subtracted histograms ..." << endl;
  vector<processes> backgroundList;
  backgroundList.clear();
  backgroundList.push_back(vjets);
  backgroundList.push_back(top);
  backgroundList.push_back(whbb);
  backgroundList.push_back(qcd);

  SubtractedHist dataMinusBkg[nResolvedChannels];
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    dataMinusBkg[chan] = SubtractedHist(dataPlusMC[chan],backgroundList,signalCorrectedCombinedFitFracErr[chan]);
  }

  if ( displayHistograms ) {
    TCanvas *dataminusbkgcnv[nResolvedChannels];
    for (int ch=0; ch<nResolvedChannels;ch++) {
      channels chan = (channels)ch;
      dataminusbkgcnv[chan]=dataMinusBkg[chan].drawSubtractedHistsWithTitle("Subtracted Data vs Diboson MC : "+chanTitleLabel[chan]);
      if ( saveHistograms ) { saveCanvasAs(dataminusbkgcnv[chan],plotDir + "dijetpT_SvdUnfolding_DataUnfold_" + chanPlotLabel[chan] + "_Subtracted"+outputSuffix+".png"); }//keep the same prefix, regardless of the unfolding procedure, to avoid duplication
    }
  }
  cout << endl;

  /// 4. Combine stacked and subtracted histograms for the four channels into one
  cout << "Combining resolved channels ..." << endl;
  DataAndMCStack<processes> resolvedDataPlusMC(emptyChanStack);
  SubtractedHist dataMinusBkgTotal;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    resolvedDataPlusMC+=dataPlusMC[chan];
    dataMinusBkgTotal+=dataMinusBkg[chan];
  }
  resolvedDataPlusMC.remakeIndividualFileProcessStackHistAndLgnd();
  cout << "Four resolved channels combined : " << endl;
  cout << resolvedDataPlusMC.eventsInStackHist() << " rescaled MC events vs " << resolvedDataPlusMC.getDataEvents() << " Data events"  << endl;

  if ( displayHistograms ) {
    TCanvas *datavsmctotcnv;
    resolvedDataPlusMC.setStackHistTitle("");
    datavsmctotcnv = resolvedDataPlusMC.drawDataAndMCOnCanvas("datavsmctotcnv");
    resolvedDataPlusMC.placeCMSLumiFrame(datavsmctotcnv,2,0);
    datavsmctotcnv->Update();
    if ( saveHistograms ) { saveCanvasAs(datavsmctotcnv,plotDir + "dijetpT_SvdUnfolding_DataUnfold_AllChannels_Stacked"+outputSuffix+".png"); }

    TCanvas *dataminusbkgtotcnv;
    dataminusbkgtotcnv=dataMinusBkgTotal.drawSubtractedHistsWithTitle("");
    dataMinusBkgTotal.placeCMSLumiFrame(dataminusbkgtotcnv,2,0);
    dataminusbkgtotcnv->Update();
    if ( saveHistograms ) { saveCanvasAs(dataminusbkgtotcnv,plotDir + "dijetpT_SvdUnfolding_DataUnfold_AllChannels_Subtracted"+outputSuffix+".png"); }
  }
  cout << endl;
  cout << endl;
  cout << endl;

//////////////////////////////       Perform The Unfolding      ///////////////////////////////
  cout << "UNFOLDING ..." << endl;

  /// 5a) Unfold using the response matrix computed based on PYTHIA diboson MC
  cout << "Unfolding with PYTHIA ..." << endl;
  cout << "Loading process info ..." << endl;
  SingleProcess<processes> pythiaDiboson(diboson);
  vector< pair<TString,double> > currentPythiaProcess;
  double expectedYieldTotal=0;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    cout << "Loading channel " << chanTitleLabel[chan] << " channel" << endl;
    currentPythiaProcess = dibosonPythiaFilesAndCoeffs[chan];
    for ( vector< pair<TString,double> >::iterator fileIter = currentPythiaProcess.begin(); fileIter != currentPythiaProcess.end(); ++fileIter ) {
      pythiaDiboson+=SingleFileInput(fileIter->first,pTjj,fileIter->second,defaultCutString[chan],mcEventWeightString);
    }
    expectedYieldTotal+=expectedProcessYields[chan][diboson];
  }
  pythiaDiboson.remakeProcessAndFileHistograms();
  cout << "PYTHIA based diboson :" << endl;
  cout << pythiaDiboson.getProcessEvents() << " events expected based on contributing files vs. " << expectedYieldTotal << " events expected for aMC@NLO with mjj cuts" << endl;
  cout << endl;

  cout << "Unfolding MC ..." << endl;
  if ( performSVDUnfolding ) {
    cout << "Run SVD Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingPythiaMCSVD(pythiaDiboson,true,canvas900x600);
    pTUnfoldingPythiaMCSVD.runSVDUnfolding();
    pTUnfoldingPythiaMCSVD.fillRecoErrorHistograms();
    //   cout << "Hist Info Before Rescaling:" << endl;
    //   pTUnfoldingPythiaMCSVD.printHistogramEntriesAndIntegrals();
    //   cout << endl;
    //pTUnfoldingPythiaMCSVD.printRawResponseProperties();
    pTUnfoldingPythiaMCSVD.rescaleAllHistogramsByLumi();
    //   cout << "Hist Info After Rescaling:" << endl;
    //   pTUnfoldingPythiaMCSVD.printHistogramEntriesAndIntegrals();
    //   cout << endl;
    pTUnfoldingPythiaMCSVD.formatAllHistograms();
    //  pTUnfoldingPythiaMCSVD.drawHistogramsSeparatelyAndSaveWithPrefix(plotDir+"TestPYTHIA_");
    //pTUnfoldingPythiaMCSVD.printRecoHistToRawTableFile(plotDir+"dijetpT_SvdUnfolding_MCUnfold_PYTHIA_UnfoldingTable_Raw"+outputSuffix+".txt");
    //pTUnfoldingPythiaMCSVD.printRecoHistToFormattedTableFile(plotDir+"dijetpT_SvdUnfolding_MCUnfold_PYTHIA_UnfoldingTable"+outputSuffix+".txt");
    TCanvas *unfoldedPythiaMCSVDCnv = pTUnfoldingPythiaMCSVD.drawOverlayedHistsWithTitle("");
    pTUnfoldingPythiaMCSVD.placeCMSLumiFrame(unfoldedPythiaMCSVDCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedPythiaMCSVDCnv,plotDir+"dijetpT_SvdUnfolding_MCUnfold_PYTHIA_AllChannels_Unfolded"+outputSuffix+".png"); }
  }

  if ( performBayesUnfolding ) {
    cout << "Run Bayes Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingPythiaMCBayes(pythiaDiboson,true,canvas900x600);
    pTUnfoldingPythiaMCBayes.runBayesUnfolding();
    pTUnfoldingPythiaMCBayes.fillRecoErrorHistograms();
    pTUnfoldingPythiaMCBayes.rescaleAllHistogramsByLumi();
    pTUnfoldingPythiaMCBayes.formatAllHistograms();
    TCanvas *unfoldedPythiaMCBayesCnv = pTUnfoldingPythiaMCBayes.drawOverlayedHistsWithTitle("");
    pTUnfoldingPythiaMCBayes.placeCMSLumiFrame(unfoldedPythiaMCBayesCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedPythiaMCBayesCnv,plotDir + "dijetpT_BayesUnfolding_MCUnfold_PYTHIA_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }


  cout << "Unfolding Data (based on a response computed with PYTHIA) ..." << endl;
  if ( performSVDUnfolding ) {
    cout << "Run SVD Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingPythiaDataSVD(pythiaDiboson,false,canvas900x600,dataMinusBkgTotal.makeDataPlusFullErrHist());
    pTUnfoldingPythiaDataSVD.runSVDUnfolding();
    pTUnfoldingPythiaDataSVD.fillRecoErrorHistograms();
    pTUnfoldingPythiaDataSVD.printResponsePropertiesToTableFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_PYTHIA_ResponseAndSummaryTable"+outputSuffix+".txt");
    pTUnfoldingPythiaDataSVD.rescaleAllHistogramsByLumi();
    pTUnfoldingPythiaDataSVD.formatAllHistograms();
    //   pTUnfoldingPythiaDataSVD.printHistogramEntriesAndIntegrals();
    //   pTUnfoldingPythiaDataSVD.drawHistogramsSeparatelyAndSave();
    pTUnfoldingPythiaDataSVD.printRecoHistToRawTableFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_PYTHIA_UnfoldingTable_Raw"+outputSuffix+".txt");
    pTUnfoldingPythiaDataSVD.printRecoHistToFormattedTableFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_PYTHIA_UnfoldingTable"+outputSuffix+".txt");
    TCanvas *unfoldedPythiaDataSVDCnv = pTUnfoldingPythiaDataSVD.drawOverlayedHistsWithTitle("");
    pTUnfoldingPythiaDataSVD.placeCMSLumiFrame(unfoldedPythiaDataSVDCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedPythiaDataSVDCnv,plotDir+"dijetpT_SvdUnfolding_DataUnfold_PYTHIA_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }

  if ( performBayesUnfolding ) {
    cout << "Run Bayes Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingPythiaDataBayes(pythiaDiboson,false,canvas900x600,dataMinusBkgTotal.makeDataPlusFullErrHist());
    pTUnfoldingPythiaDataBayes.runBayesUnfolding();
    pTUnfoldingPythiaDataBayes.fillRecoErrorHistograms();
    pTUnfoldingPythiaDataBayes.printResponsePropertiesToTableFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_PYTHIA_ResponseAndSummaryTable"+outputSuffix+".txt");
    pTUnfoldingPythiaDataBayes.rescaleAllHistogramsByLumi();
    pTUnfoldingPythiaDataBayes.formatAllHistograms();
    pTUnfoldingPythiaDataBayes.printRecoHistToRawTableFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_PYTHIA_UnfoldingTable_Raw"+outputSuffix+".txt");
    pTUnfoldingPythiaDataBayes.printRecoHistToFormattedTableFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_PYTHIA_UnfoldingTable"+outputSuffix+".txt");
    TCanvas *unfoldedPythiaDataBayesCnv = pTUnfoldingPythiaDataBayes.drawOverlayedHistsWithTitle("");
    pTUnfoldingPythiaDataBayes.placeCMSLumiFrame(unfoldedPythiaDataBayesCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedPythiaDataBayesCnv,plotDir+"dijetpT_BayesUnfolding_DataUnfold_PYTHIA_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
    cout << endl;
    cout << endl;
  }


  /// 5b) Unfold using the response matrix computed based on aMC@NLO diboson MC
  cout << "Unfolding with aMC@NLO ..." << endl;
  cout << "Unfolding MC ..." << endl;
//   cout << "Reloading process info directly from MC ..." << endl;
//   SingleProcess<processes> amcnloDiboson(diboson);
//   vector< pair<TString,double> > currentAMCNLOProcess;
//   double expectedYieldTotalforAMCNLO=0;
//   for (int ch=0; ch<nResolvedChannels;ch++) {
//     channels chan = (channels)ch;
//     cout << "Loading channel " << chanTitleLabel[chan] << " channel" << endl;
//     currentAMCNLOProcess = simFilesAndCoeffs[chan][diboson];
//     for ( vector< pair<TString,double> >::iterator fileIter = currentAMCNLOProcess.begin(); fileIter != currentAMCNLOProcess.end(); ++fileIter ) {
//       amcnloDiboson+=SingleFileInput(fileIter->first,pTjj,fileIter->second,defaultCutString[chan],mcEventWeightString);
//     }
//     expectedYieldTotalforAMCNLO+=expectedProcessYields[chan][diboson];
//   }
//   amcnloDiboson.remakeProcessAndFileHistograms();
//   cout << "aMC@NLO based diboson :" << endl;
//   cout << amcnloDiboson.getProcessEvents() << " events expected based on contributing files vs. " << expectedYieldTotalforAMCNLO  << " events expected for aMC@NLO with mjj cuts" << endl;
//   cout << endl;
  if ( performSVDUnfolding ) {
    cout << "Run SVD Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingAMCNLOMCSVD(resolvedDataPlusMC.getStackMap()[diboson],true,canvas900x600);
    //UnfoldDibosonInPt pTUnfoldingAMCNLOMCSVD(amcnloDiboson,true,canvas900x600);
    pTUnfoldingAMCNLOMCSVD.runSVDUnfolding();
    pTUnfoldingAMCNLOMCSVD.fillRecoErrorHistograms();
    //   cout << "Hist Info Before Rescaling:" << endl;
    //   pTUnfoldingAMCNLOMCSVD.printHistogramEntriesAndIntegrals();
    pTUnfoldingAMCNLOMCSVD.rescaleAllHistogramsByLumi();
    //   cout << "Hist Info After Rescaling:" << endl;
    //   pTUnfoldingAMCNLOMCSVD.printHistogramEntriesAndIntegrals();
    pTUnfoldingAMCNLOMCSVD.formatAllHistograms();
    //  pTUnfoldingAMCNLOMCSVD.drawHistogramsSeparatelyAndSaveWithPrefix(plotDir+"TestAMCNLO_");
    TCanvas *unfoldedAMCNLOMCSVDCnv = pTUnfoldingAMCNLOMCSVD.drawOverlayedHistsWithTitle("");
    pTUnfoldingAMCNLOMCSVD.placeCMSLumiFrame(unfoldedAMCNLOMCSVDCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedAMCNLOMCSVDCnv,plotDir + "dijetpT_SvdUnfolding_MCUnfold_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }

  if ( performBayesUnfolding ) {
    cout << "Run Bayes Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingAMCNLOMCBayes(resolvedDataPlusMC.getStackMap()[diboson],true,canvas900x600);
    pTUnfoldingAMCNLOMCBayes.runBayesUnfolding();
    pTUnfoldingAMCNLOMCBayes.fillRecoErrorHistograms();
    pTUnfoldingAMCNLOMCBayes.rescaleAllHistogramsByLumi();
    pTUnfoldingAMCNLOMCBayes.formatAllHistograms();
    TCanvas *unfoldedAMCNLOMCBayesCnv = pTUnfoldingAMCNLOMCBayes.drawOverlayedHistsWithTitle("");
    pTUnfoldingAMCNLOMCBayes.placeCMSLumiFrame(unfoldedAMCNLOMCBayesCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedAMCNLOMCBayesCnv,plotDir + "dijetpT_BayesUnfolding_MCUnfold_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }

  cout << "Unfolding Data ..." << endl;
  if ( performSVDUnfolding ) {
    cout << "Run SVD Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingAMCNLODataSVD(resolvedDataPlusMC.getStackMap()[diboson],false,canvas900x600,dataMinusBkgTotal.makeDataPlusFullErrHist());
    //cout << "Response Properties:" << endl;
    //pTUnfoldingAMCNLODataSVD.printResponseProperties();
    pTUnfoldingAMCNLODataSVD.runSVDUnfolding();
    pTUnfoldingAMCNLODataSVD.fillRecoErrorHistograms();
    pTUnfoldingAMCNLODataSVD.printResponsePropertiesToTableFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_ResponseAndSummaryTable"+outputSuffix+".txt");
    pTUnfoldingAMCNLODataSVD.rescaleAllHistogramsByLumi();
    pTUnfoldingAMCNLODataSVD.augmentHistRecoUncertaintiesBasedOnAlternateMCTableInFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_PYTHIA_UnfoldingTable_Raw"+outputSuffix+".txt");
    pTUnfoldingAMCNLODataSVD.formatAllHistograms();
    pTUnfoldingAMCNLODataSVD.printRecoHistToFormattedTableFile(plotDir+"dijetpT_SvdUnfolding_DataUnfold_UnfoldingTable"+outputSuffix+".txt");
    TCanvas *unfoldedAMCNLODataSVDCnv = pTUnfoldingAMCNLODataSVD.drawOverlayedHistsWithTitle("");
    pTUnfoldingAMCNLODataSVD.placeCMSLumiFrame(unfoldedAMCNLODataSVDCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedAMCNLODataSVDCnv,plotDir + "dijetpT_SvdUnfolding_DataUnfold_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }

  if ( performBayesUnfolding ) {
    cout << "Run Bayes Unfolding ..." << endl;
    UnfoldDibosonInPt pTUnfoldingAMCNLODataBayes(resolvedDataPlusMC.getStackMap()[diboson],false,canvas900x600,dataMinusBkgTotal.makeDataPlusFullErrHist());
    pTUnfoldingAMCNLODataBayes.runBayesUnfolding();
    pTUnfoldingAMCNLODataBayes.fillRecoErrorHistograms();
    pTUnfoldingAMCNLODataBayes.printResponsePropertiesToTableFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_ResponseAndSummaryTable"+outputSuffix+".txt");
    pTUnfoldingAMCNLODataBayes.rescaleAllHistogramsByLumi();
    pTUnfoldingAMCNLODataBayes.augmentHistRecoUncertaintiesBasedOnAlternateMCTableInFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_PYTHIA_UnfoldingTable_Raw"+outputSuffix+".txt");
    pTUnfoldingAMCNLODataBayes.formatAllHistograms();
    pTUnfoldingAMCNLODataBayes.printRecoHistToFormattedTableFile(plotDir+"dijetpT_BayesUnfolding_DataUnfold_UnfoldingTable"+outputSuffix+".txt");
    TCanvas *unfoldedAMCNLODataBayesCnv = pTUnfoldingAMCNLODataBayes.drawOverlayedHistsWithTitle("");
    pTUnfoldingAMCNLODataBayes.placeCMSLumiFrame(unfoldedAMCNLODataBayesCnv,2,0);
    if ( saveHistograms ) { saveCanvasAs(unfoldedAMCNLODataBayesCnv,plotDir + "dijetpT_BayesUnfolding_DataUnfold_AllChannels_Unfolded"+outputSuffix+".png"); }
    cout << endl;
  }

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//


#ifndef __CINT__
int main () { UnfoldDiboson(); return 0; }  // Main program when run stand-alone
#endif
