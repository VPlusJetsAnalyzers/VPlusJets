//// Namespace of inputs to the diboson analysis
//// Created 01/19/2016
#include <vector>
#include <map>
#include <utility>
#include <cstdio>
#include <iostream>

#include <TString.h>

using std::vector;
using std::map;
using std::pair;
using std::make_pair;

using std::cout;
using std::endl;

namespace dibosonAnalysis {

/////////////////////////////////////  Basic Definitions:  /////////////////////////////////////
  const int nChannels=6;
  const int nResolvedChannels=4;
  enum channels{mu_antibtagged,el_antibtagged,mu_btagged,el_btagged,mu_boosted,el_boosted};
  const int nProcesses=6;
  const int nMCProcesses=nProcesses-1;
  enum processes{diboson,vjets,top,qcd,whbb,data};

////////////////////////////////////  Graphics And Labels:  ////////////////////////////////////
  TString chanTitleLabel[nChannels] = {"#mu - resolved, anti-btagged","el - resolved, anti-btagged","#mu - resolved, btagged","el - resolved, btagged","#mu - boosted","el - boosted"};
  TString chanPlotLabel[nChannels] = {"muResolvedAntibtagged","elResolvedAntibtagged","muResolvedBtagged","elResolvedBtagged","muBoosted","elBoosted"};
  TString procLabel[nProcesses-1] = {"WW+WZ", "V+jets", "top", "multijet", "WHbb" };
  Color_t procColor[nProcesses-1] = {kAzure+8,kRed,kGreen+2,kGray,kBlue};
  processes processesInStackOrder[nProcesses-1] = {diboson,whbb,top,qcd,vjets};
  vector< pair< processes,pair<TString,Color_t> > > orderedNameColorPairs;
  Color_t legendColor = kWhite;

//////////////////////////  Fit/Validation Input Values and Results:  //////////////////////////
  double expectedProcessYields[nChannels][nProcesses] = { 
    {2261.8,118257,4568.2,-1.0,-1.0},//mu_antibtagged
    {2046.1,102230,3872.8,8519.6,-1.0},//el_antibtagged
    {172.2,6884.5,7491.3,-1.0,22.7},//mu_btagged
    {124.7,6247.1,6180.2,-1.0,20.1},//el_btagged
    {205.3,2636.3,706.4,-1.0,-1.0},//mu_boosted
    {171.6,2149.1,578.6,-1.0,-1.0}//el_boosted
  };
  double expectedDataYield[nChannels] = {126764,121709,14001,12586,4194,3481};
  double signalAxEff[nChannels] = {0.00142,0.00129,0.000108,0.0000786,0.000132,0.000111};
  double lumiMuon = 19300.0;
  double lumiElectron = 19200.0;
  double lumi[nChannels] = {lumiMuon,lumiElectron,lumiMuon,lumiElectron,lumiMuon,lumiElectron};

  double fitFracVal[nChannels][nProcesses] = { 
    {1.7720,0.9994,1.0004,-1.0,-1.0},//mu_antibtagged
    {1.9635,1.0461,1.0023,0.8033,-1.0},//el_antibtagged
    {1.4861,0.9777,0.9367,-1.0,1.0},//mu_btagged
    {1.7065,1.0068,0.9815,-1.0,1.0},//el_btagged
    {1.3474,1.2187,0.9965,-1.0,-1.0},//mu_boosted
    {2.2266,1.1724,0.9997,-1.0,-1.0}//el_boosted
  };
  double fitFracErr[nChannels][nProcesses] = { 
    {0.359,0.0075,0.060,-1.0,-1.0},//mu_antibtagged
    {0.390,0.0271,0.060,0.340,-1.0},//el_antibtagged
    {0.792,0.0372,0.034,-1.0,0.0},//mu_btagged
    {0.881,0.0360,0.036,-1.0,0.0},//el_btagged
    {0.769,0.0616,0.080,-1.0,-1.0},//mu_boosted
    {0.840,0.0694,0.100,-1.0,-1.0}//el_boosted
  };

  double signalYieldBias[nChannels] = {-263.0,131.5,-6.0,23.5,17.4,12.8};
  double signalErrorInflation[nChannels] = {0.877,0.842,0.909,0.926,0.850,0.986};
  double signalQuadError[nChannels] = {331.2,139.7,43.1,18.2,77.3,55.8};//From studies with alternate W+Jets parameterizations
  double signalCorrelatedError[nChannels] = {796.0,224.1,61.1,56.6,0.0,0.0};//MC Syst errors 
  double signalCorrectedFitFracVal[nChannels];
  double signalCorrectedCombinedFitFracErr[nChannels];

///////////////////////////////  Input Files And Coefficients:  ////////////////////////////////
  /// Coefficients and corrections:
  ///Diboson: realign cross-section to the theory value and generation with the full phase space similar to the fit configuration
  double crossXWWaMCNLO=57.52;
  double fracPhaseSpaceWWmain=0.444408;
  double fracPhaseSpaceWWtauleptonic=0.0246772*2;
  //double fracPhaseSpaceWW=fracPhaseSpaceWWmain+fracPhaseSpaceWWtauleptonic;
  double crossXWZaMCNLO=24.20;
  double fracPhaseSpaceWZmain=0.230996;
  double fracPhaseSpaceWZtauleptonic=0.0114347;
  //double fracPhaseSpaceWZ=fracPhaseSpaceWZmain+fracPhaseSpaceWZtauleptonic;
  //Take the relative fractions of WW to WZ from the same MC (aMCNLO), but the ovarall cross section from theory
  double crossXDibosonTheory=59.8+22.9;
  double crossXWW=crossXDibosonTheory*crossXWWaMCNLO/(crossXWWaMCNLO+crossXWZaMCNLO);
  double crossXWZ=crossXDibosonTheory*crossXWZaMCNLO/(crossXWWaMCNLO+crossXWZaMCNLO);
  /// Vjets: different kfactors for different channels
  double vjetsKfactor_antibtagged=1.16;
  double vjetsKfactor_btagged=2.05;
  /// QCD: explicitly accounted for electrons_antibtagged
  double qcdDirectCorrection=8519.6/1728;//Corrects QCD to give the right event count in the mjj window (which is calculated by the fit machinery as a fraction of the data)

  /// Combined Inputs:
  TString inputFileDir = "/eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/";
  //  TString signalFileDir = "/uscms_data/d3/ilyao/Unfolding/RooUnfold-1.1.1/inputfiles/";
  TString signalFileDir = "/eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/UnfoldingInput";
  map<channels, map<processes,vector< pair<TString,double> > > > simFilesAndCoeffs;//Typically Coeff=kFactor*Lumi*CrossX/NGenerated
  map<channels,TString> dataFiles;
  map<channels,vector< pair<TString,double> > > dibosonPythiaFilesAndCoeffs;

////////////////////////////////  Cuts And Analysis Variables:  ////////////////////////////////
  TString defaultCutString[nChannels] = { "( (sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)&&(Mass2j_PFCor>48.000)&&(Mass2j_PFCor<160.000) )", //mu_antibtagged
					  "( ((sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(W_electron_pt>30))&&(Mass2j_PFCor>48.000)&&(Mass2j_PFCor<160.000) )", //el_antibtagged
					  "( ((sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&(JetPFCor_bDiscriminatorCSV[0]>0.244)&&(JetPFCor_bDiscriminatorCSV[1]>0.244)&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.))&&(Mass2j_PFCor>40.000)&&(Mass2j_PFCor<160.000) )", //mu_btagged
					  "( ((sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&(JetPFCor_bDiscriminatorCSV[0]>0.244)&&(JetPFCor_bDiscriminatorCSV[1]>0.244)&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(W_electron_pt>30))&&(Mass2j_PFCor>40.000)&&(Mass2j_PFCor<160.000) )", //el_btagged
					  "( ((W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&(GroomedJet_CA8_pt[1]<80)&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet[0]>2.0)&&(GroomedJet_CA8_deltaR_lca8jet[0]>1.57)&&(numPFCorJetBTags<1)&&(vbf_event==0)&&(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.))&&(GroomedJet_CA8_mass_pr[0]>40.000)&&(GroomedJet_CA8_mass_pr[0]<140.000) )", //mu_boosted
					  "( ((W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&(GroomedJet_CA8_pt[1]<80)&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet[0]>2.0)&&(GroomedJet_CA8_deltaR_lca8jet[0]>1.57)&&(numPFCorJetBTags<1)&&(vbf_event==0)&&(event_met_pfmet >70)&&(W_electron_pt>35))&&(GroomedJet_CA8_mass_pr[0]>40.000)&&(GroomedJet_CA8_mass_pr[0]<140.000) )" //el_boosted
  };
  TString qcdCutString[nChannels] = { "N/A",
				      "( ((sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>20)&&(W_electron_pfIsoEA>0.3)&&(W_electron_pt>30))&&(Mass2j_PFCor>48.000)&&(Mass2j_PFCor<160.000) )",
				      "N/A","N/A","N/A","N/A" 
  };				 
  TString mcEventWeightString = "effwt*puwt";
  TString dataEventWeightString = "";
  TString defaultVarName="Mass2j_PFCor";// AnalysisVariable mjj("Mass2j_PFCor",48,160,14);
  double defaultVarMin=48, defaultVarBtagMin = 40, defaultVarMax=160;
  int defaultVarNBins=14, defaultVarBtagNBins = 15;

};

void initializeCombinedInputs () {
  using namespace dibosonAnalysis;

  simFilesAndCoeffs[mu_antibtagged][diboson].clear();
  simFilesAndCoeffs[mu_antibtagged][diboson].push_back( make_pair(signalFileDir+"FT_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NominalWWpTwt.root",lumiMuon*crossXWW*fracPhaseSpaceWWmain/799899.) );
  //simFilesAndCoeffs[mu_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532.root",lumiMuon*crossXWW*fracPhaseSpaceWWmain/799899.) );//Use if FT file is not available for all cuts
  simFilesAndCoeffs[mu_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_mu_WZtoLNuQQ_amcnlo_100k_CMSSW532.root",lumiMuon*crossXWZ*fracPhaseSpaceWZmain/99791.) );
  simFilesAndCoeffs[mu_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_mu_WW_lvtauvtau_CMSSW532.root",lumiMuon*crossXWW*fracPhaseSpaceWWtauleptonic/(24898.+24898.)) );
  simFilesAndCoeffs[mu_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_mu_WZ_lvtautau_CMSSW532.root",lumiMuon*crossXWZ*fracPhaseSpaceWZtauleptonic/(24897.+24899.)) );
  simFilesAndCoeffs[el_antibtagged][diboson].clear();
  simFilesAndCoeffs[el_antibtagged][diboson].push_back( make_pair(signalFileDir+"FT_el_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NominalWWpTwt.root",lumiElectron*crossXWW*fracPhaseSpaceWWmain/799899.) );
  //simFilesAndCoeffs[el_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_el_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532.root",lumiElectron*crossXWW*fracPhaseSpaceWWmain/799899.) );//Use if FT file is not available for all cuts
  simFilesAndCoeffs[el_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_el_WZtoLNuQQ_amcnlo_100k_CMSSW532.root",lumiElectron*crossXWZ*fracPhaseSpaceWZmain/99791.) );
  simFilesAndCoeffs[el_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_el_WW_lvtauvtau_CMSSW532.root",lumiElectron*crossXWW*fracPhaseSpaceWWtauleptonic/(24898.+24898.)) );
  simFilesAndCoeffs[el_antibtagged][diboson].push_back( make_pair(signalFileDir+"RD_el_WZ_lvtautau_CMSSW532.root",lumiElectron*crossXWZ*fracPhaseSpaceWZtauleptonic/(24897.+24899.)) );
  simFilesAndCoeffs[mu_btagged][diboson] = simFilesAndCoeffs[mu_antibtagged][diboson];
  simFilesAndCoeffs[el_btagged][diboson] = simFilesAndCoeffs[el_antibtagged][diboson];

  simFilesAndCoeffs[mu_antibtagged][vjets].clear();
  simFilesAndCoeffs[mu_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_mu_W1Jets_CMSSW532.root",lumiMuon*vjetsKfactor_antibtagged*5400.0/19871598.0) );
  simFilesAndCoeffs[mu_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_mu_W2Jets_CMSSW532.root",lumiMuon*vjetsKfactor_antibtagged*1750.0/33004921.0) );
  simFilesAndCoeffs[mu_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_mu_W3Jets_CMSSW532.root",lumiMuon*vjetsKfactor_antibtagged*519.0/15059503.0) );
  simFilesAndCoeffs[mu_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_mu_W4Jets_CMSSW532.root",lumiMuon*vjetsKfactor_antibtagged*214.0/12842803.0) );
  simFilesAndCoeffs[mu_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_mu_ZpJ_CMSSW532.root",lumiMuon*3503.71/30209426.0) );
  simFilesAndCoeffs[el_antibtagged][vjets].clear();
  simFilesAndCoeffs[el_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_el_W1Jets_CMSSW532.root",lumiElectron*vjetsKfactor_antibtagged*5400.0/19871598.0) );
  simFilesAndCoeffs[el_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_el_W2Jets_CMSSW532.root",lumiElectron*vjetsKfactor_antibtagged*1750.0/33004921.0) );
  simFilesAndCoeffs[el_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_el_W3Jets_CMSSW532.root",lumiElectron*vjetsKfactor_antibtagged*519.0/15059503.0) );
  simFilesAndCoeffs[el_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_el_W4Jets_CMSSW532.root",lumiElectron*vjetsKfactor_antibtagged*214.0/12842803.0) );
  simFilesAndCoeffs[el_antibtagged][vjets].push_back( make_pair(inputFileDir+"RD_el_ZpJ_CMSSW532.root",lumiElectron*3503.71/30209426.0) );
  simFilesAndCoeffs[mu_btagged][vjets] = simFilesAndCoeffs[mu_antibtagged][vjets];
  simFilesAndCoeffs[el_btagged][vjets] = simFilesAndCoeffs[el_antibtagged][vjets];
  for (int i=0; i<4;i++) { //Correct the btagged kfactors for W+Jets 
    // simFilesAndCoeffs[mu_btagged][vjets][i].second=(simFilesAndCoeffs[mu_btagged][vjets][i].second)*vjetsKfactor_btagged/vjetsKfactor_antibtagged;
    // simFilesAndCoeffs[el_btagged][vjets][i].second=(simFilesAndCoeffs[el_btagged][vjets][i].second)*vjetsKfactor_btagged/vjetsKfactor_antibtagged;
    simFilesAndCoeffs[mu_btagged][vjets][i].second*=vjetsKfactor_btagged/vjetsKfactor_antibtagged;
    simFilesAndCoeffs[el_btagged][vjets][i].second*=vjetsKfactor_btagged/vjetsKfactor_antibtagged;
  }

  simFilesAndCoeffs[mu_antibtagged][top].clear();
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_TTbar_CMSSW532.root",lumiMuon*246.73/6893735.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopT_T_CMSSW532.root",lumiMuon*55.531/3758221.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopT_Tbar_CMSSW532.root",lumiMuon*30.0042/1935066.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopS_T_CMSSW532.root",lumiMuon*3.89394/259960.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopS_Tbar_CMSSW532.root",lumiMuon*1.75776/139974.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopTW_T_CMSSW532.root",lumiMuon*11.1773/497657.0) );
  simFilesAndCoeffs[mu_antibtagged][top].push_back( make_pair(inputFileDir+"RD_mu_STopTW_Tbar_CMSSW532.root",lumiMuon*11.1773/493458.0) );
  simFilesAndCoeffs[el_antibtagged][top].clear();
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_TTbar_CMSSW532.root",lumiElectron*246.73/6893735.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopT_T_CMSSW532.root",lumiElectron*55.531/3758221.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopT_Tbar_CMSSW532.root",lumiElectron*30.0042/1935066.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopS_T_CMSSW532.root",lumiElectron*3.89394/259960.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopS_Tbar_CMSSW532.root",lumiElectron*1.75776/139974.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopTW_T_CMSSW532.root",lumiElectron*11.1773/497657.0) );
  simFilesAndCoeffs[el_antibtagged][top].push_back( make_pair(inputFileDir+"RD_el_STopTW_Tbar_CMSSW532.root",lumiElectron*11.1773/493458.0) );
  simFilesAndCoeffs[mu_btagged][top] = simFilesAndCoeffs[mu_antibtagged][top];
  simFilesAndCoeffs[el_btagged][top] = simFilesAndCoeffs[el_antibtagged][top];

  simFilesAndCoeffs[el_antibtagged][qcd].clear();
  simFilesAndCoeffs[el_antibtagged][qcd].push_back( make_pair(inputFileDir+"RDQCD_WenuJets_Isog0p3NoElMVA_19p2invfb.root",qcdDirectCorrection) );

  simFilesAndCoeffs[mu_btagged][whbb].clear();
  simFilesAndCoeffs[mu_btagged][whbb].push_back( make_pair(inputFileDir+"RD_mu_WH_WToLNu_HToBB_M-125_CMSSW532.root",lumiMuon*0.6966*0.577*(0.1075+0.1057+0.1125)/999998.) );
  simFilesAndCoeffs[el_btagged][whbb].clear();
  simFilesAndCoeffs[el_btagged][whbb].push_back( make_pair(inputFileDir+"RD_el_WH_WToLNu_HToBB_M-125_CMSSW532.root",lumiElectron*0.6966*0.577*(0.1075+0.1057+0.1125)/999998.) );

  dataFiles[mu_antibtagged] = inputFileDir + "RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root";
  dataFiles[el_antibtagged] = inputFileDir + "RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root";
  dataFiles[mu_btagged] = dataFiles[mu_antibtagged];
  dataFiles[el_btagged] = dataFiles[el_antibtagged];

  dibosonPythiaFilesAndCoeffs[mu_antibtagged].clear();
  dibosonPythiaFilesAndCoeffs[mu_antibtagged].push_back( make_pair(signalFileDir+"RD_mu_WW_CMSSW532.root",lumiMuon*57.25/9450414.) );
  dibosonPythiaFilesAndCoeffs[mu_antibtagged].push_back( make_pair(signalFileDir+"RD_mu_WZ_CMSSW532.root",lumiMuon*22.88/10000267.) );
  dibosonPythiaFilesAndCoeffs[el_antibtagged].clear();
  dibosonPythiaFilesAndCoeffs[el_antibtagged].push_back( make_pair(signalFileDir+"RD_el_WW_CMSSW532.root",lumiElectron*57.25/9450414.) );
  dibosonPythiaFilesAndCoeffs[el_antibtagged].push_back( make_pair(signalFileDir+"RD_el_WZ_CMSSW532.root",lumiElectron*22.88/10000267.) );
  dibosonPythiaFilesAndCoeffs[mu_btagged] = dibosonPythiaFilesAndCoeffs[mu_antibtagged];
  dibosonPythiaFilesAndCoeffs[el_btagged] = dibosonPythiaFilesAndCoeffs[el_antibtagged];

  ///Obtain compute the signal Yield and Error after the bias corrections & systematics studies
  double corrEvtErr=-1.0;
  for (int ch=0; ch<nChannels;ch++) {
    signalCorrectedFitFracVal[ch]=(fitFracVal[ch][diboson]*expectedProcessYields[ch][diboson]-signalYieldBias[ch])/expectedProcessYields[ch][diboson];
    corrEvtErr=expectedProcessYields[ch][diboson]*fitFracErr[ch][diboson]*signalErrorInflation[ch];
    signalCorrectedCombinedFitFracErr[ch]=(sqrt(corrEvtErr*corrEvtErr+signalQuadError[ch]*signalQuadError[ch]+signalCorrelatedError[ch]*signalCorrelatedError[ch]))/expectedProcessYields[ch][diboson];
  }

  ///Initialize the graphics options
  orderedNameColorPairs.clear();
  for (int i=0; i<nProcesses-1;i++) {
    processes proc = processesInStackOrder[i];
    orderedNameColorPairs.push_back(  make_pair( proc, make_pair(procLabel[proc],procColor[proc]) )  );
  }
}

void printSummaryTables() 
/// Print the tables for the documentation
{
  initializeCombinedInputs();
  using namespace dibosonAnalysis;
  char buffer[100];
  TString tableLine;

  std::cout << "SMP-13-008 - tab:yields - Expected event yields and ..." << std::endl;
  // TString predictedFormatString="%-30s & %-18.0f & %-18.0f & %-18.0f & %-18.0f & %-18.0f & %-18.0f \\\\ \n";
  // TString fitResultFormatString="%-30s & %-4.2f$\\pm$%-9.2f & %-4.2f$\\pm$%-9.2f & %-4.2f$\\pm$%-9.2f & %-4.2f$\\pm$%-9.2f & %-4.2f$\\pm$%-9.2f & %-4.2f$\\pm$%-9.2f \\\\ \n";
  char bkgPredictedLabel[nProcesses][30]={"NA","V+Jets Predicted","Top Predicted","Multijet Predicted","$WHbb(125)$ Predicted","NA"};
  char fitResultLabel[nProcesses][30]={"NA","$Frac_{V+Jets}$","$Frac_{Top}$","$Frac_{Multijet}$","$Frac_{WHbb(125)}$","NA"};
  char chanTitle[nChannels][30]={"$\\Pgm$ anti-b-tag","$\\Pe$ anti-b-tag","$\\Pgm$ b-tag","$\\Pe$ b-tag","$\\Pgm$ merged","$\\Pe$ merged"};

  sprintf(buffer,"%-30s","Quantity");
  tableLine=buffer;
  for (int ch=0; ch<nChannels;ch++) {
    sprintf(buffer," & %-18s",chanTitle[ch]);
    tableLine+=buffer;
  }
  cout << tableLine << " \\\\" << endl;
  cout << "\\hline" << endl;

  sprintf(buffer,"%-30s","Data Events");
  tableLine=buffer;
  for (int ch=0; ch<nChannels;ch++) {
    sprintf(buffer," & %-18.0f",expectedDataYield[ch]);
    tableLine+=buffer;
  }
  cout << tableLine << " \\\\" << endl;
  for (int prc=1; prc<(nProcesses-1);prc++) {//skip diboson(for now) and data(above)
    sprintf(buffer,"%-30s",bkgPredictedLabel[prc]);
    tableLine=buffer;
    for (int ch=0; ch<nChannels;ch++) {
      if ( expectedProcessYields[ch][prc]<-0.99 ) {
	sprintf(buffer," & %-18s","---");
      } else {
	sprintf(buffer," & %-18.0f",expectedProcessYields[ch][prc]);
      }
      tableLine+=buffer;
    }
    cout << tableLine << " \\\\" << endl;
    sprintf(buffer,"%-30s",fitResultLabel[prc]);
    tableLine=buffer;
    for (int ch=0; ch<nChannels;ch++) {
      if ( fitFracVal[ch][prc]<-0.99 ) {
	sprintf(buffer," & %-18s","---");
      } else {
	if ( fitFracErr[ch][prc]<0.000001 ) {
	  sprintf(buffer," & %-18.2f",fitFracVal[ch][prc]);
	} else {
	  sprintf(buffer," & %-4.2f$\\pm$%-9.2f",fitFracVal[ch][prc],fitFracErr[ch][prc]);
	}
      }
      tableLine+=buffer;
    }
    cout << tableLine << " \\\\" << endl;
  }
  cout << "\\hline" << endl;

  sprintf(buffer,"%-30s","$N_{WW+WZ}$ Predicted");
  tableLine=buffer;
  for (int ch=0; ch<nChannels;ch++) {
    sprintf(buffer," & %-18.0f",expectedProcessYields[ch][diboson]);
    tableLine+=buffer;
  }
  cout << tableLine << " \\\\" << endl;
  sprintf(buffer,"%-30s","$Frac_{WW+WZ}$");
  tableLine=buffer;
  for (int ch=0; ch<nChannels;ch++) {
    sprintf(buffer," & %-4.2f$\\pm$%-9.2f",fitFracVal[ch][diboson],fitFracErr[ch][diboson]);
    tableLine+=buffer;
  }
  cout << tableLine << " \\\\" << endl;
  sprintf(buffer,"%-30s","$Frac_{WW+WZ}$ Corrected");
  tableLine=buffer;
  for (int ch=0; ch<nChannels;ch++) {
    sprintf(buffer," & %-4.2f$\\pm$%-9.2f",signalCorrectedFitFracVal[ch],signalCorrectedCombinedFitFracErr[ch]);
    tableLine+=buffer;
  }
  cout << tableLine << " \\\\" << endl;
  cout << "NOTE: The last three lines of the table should come from running 'processDibosonFitOutput_BLUE.cc'" << endl;
  cout << endl;
  cout << "\v" << endl;

  ///Second Table:
  std::cout << "AN-12-464 - table:FitTotalsAndComparisons - Fractions of the expected yields ..." << std::endl;
  char tableProcessLabel[nProcesses][30]={"Diboson","V+Jets","Top","Multijet","$WH(150)\\to bb$","NA"};
  double predicted_mu=-1.0, predicted_el=-1.0, extracted_mu=-1.0, extracted_el=-1.0;
  int ch_mu=-1, ch_el=-1;

  for (int ch=0; ch<nChannels;ch+=2) {
    ch_mu=ch;
    ch_el=ch+1;
    predicted_mu=0;
    predicted_el=0;
    extracted_mu=0;
    extracted_el=0;
    cout << "Making the table for " << chanTitleLabel[ch_mu] << " and " <<  chanTitleLabel[ch_el] << endl;
    sprintf(buffer,"%-20s & %-20s & %-20s & %-20s & %-20s","Process","Predicted","Extracted Fraction","Predicted","Extracted Fraction");
    cout << buffer << " \\\\" << endl;
    for (int prc=0; prc<(nProcesses-1);prc++) {//data (&corrected diboson are listed separately)
      sprintf(buffer,"%-20s",tableProcessLabel[prc]);
      tableLine=buffer;
      if ( expectedProcessYields[ch_mu][prc]<-0.99 ) {//mu
	sprintf(buffer," & %-20s","---");
      } else {
	predicted_mu+=expectedProcessYields[ch_mu][prc];
	sprintf(buffer," & %-20.0f",expectedProcessYields[ch_mu][prc]);
      }
      tableLine+=buffer;
      if ( fitFracVal[ch_mu][prc]<-0.99 ) {
	sprintf(buffer," & %-20s","---");
      } else {
	extracted_mu+=(fitFracVal[ch_mu][prc]*expectedProcessYields[ch_mu][prc]);
	if ( fitFracErr[ch_mu][prc]<0.000001 ) {
	  sprintf(buffer," & %-20.2f",fitFracVal[ch_mu][prc]);
	} else {
	  sprintf(buffer," & %-4.2f$\\pm$%-11.2f",fitFracVal[ch_mu][prc],fitFracErr[ch_mu][prc]);
	}
      }
      tableLine+=buffer;
      if ( expectedProcessYields[ch_el][prc]<-0.99 ) {//el
	sprintf(buffer," & %-20s","---");
      } else {
	predicted_el+=expectedProcessYields[ch_el][prc];
	sprintf(buffer," & %-20.0f",expectedProcessYields[ch_el][prc]);
      }
      tableLine+=buffer;
      if ( fitFracVal[ch_el][prc]<-0.99 ) {
	sprintf(buffer," & %-20s","---");
      } else {
	extracted_el+=(fitFracVal[ch_el][prc]*expectedProcessYields[ch_el][prc]);
	if ( fitFracErr[ch_el][prc]<0.000001 ) {
	  sprintf(buffer," & %-20.2f",fitFracVal[ch_el][prc]);
	} else {
	  sprintf(buffer," & %-4.2f$\\pm$%-11.2f",fitFracVal[ch_el][prc],fitFracErr[ch_el][prc]);
	}
      }
      tableLine+=buffer;
      cout << tableLine << " \\\\" << endl;
    }
    sprintf(buffer,"%-20s & %-20.0f & %-20.0f & %-20.0f & %-20.0f","Total Yields",predicted_mu,extracted_mu,predicted_el,extracted_el);
    cout << buffer << " \\\\" << endl;
    cout << "\\hline" << endl;
    sprintf(buffer,"%-20s & %-20.0f & %-4.2f$\\pm$%-11.2f & %-20.0f & %-4.2f$\\pm$%-11.2f","Corrected Diboson",expectedProcessYields[ch_mu][diboson],signalCorrectedFitFracVal[ch_mu],signalCorrectedCombinedFitFracErr[ch_mu],expectedProcessYields[ch_el][diboson],signalCorrectedFitFracVal[ch_el],signalCorrectedCombinedFitFracErr[ch_el]);
    cout << buffer << " \\\\" << endl;
    cout << "\\hline" << endl;
    sprintf(buffer,"%-20s & %-20.0f & %-20s & %-20.0f & %-20s","Data",expectedDataYield[ch_mu],"---",expectedDataYield[ch_el],"---");
    cout << buffer << " \\\\" << endl;
    cout << "\\hline" << endl;
  }

}
