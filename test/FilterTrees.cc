/*******************************************************************
 * Project: CMS detector at the CERN, diboson analysis in the semileptonic channel
 *
 *
 * Authors:
 *
 *   Osipenkov, Ilya, Texas A&M - ilyao@fnal.gov
 *
 * Description: Use to filter the trees with the option of making additional cuts, adding corrections to the weight, changing the Jet Resolution, etc. in the process
 *
 ********************************************************************/

#include <iostream>
#include <iomanip>
#include <strstream>
//#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

using namespace std;

#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TBranch.h>
#include <TROOT.h>
#include "TLorentzVector.h"


// #include "DataFormats/FWLite/interface/Handle.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
// #include "FWCore/ParameterSet/interface/FileInPath.h"



void FilterTrees()
{

}


void FilterSingleFile(int reweighWWpT=0, TString outFileName="FT_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NominalWWpTwt.root", TString inFileName="RD_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532.root", TString outDirName ="/uscms_data/d3/ilyao/Diboson8TeV/FilteredTrees/", TString inDirName = "/eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/")
//// Generate the n-tuples with extra cuts applied
{
  bool passAllCuts=true;
//   vector < TLorentzVector > jLV, lLV, METLV;
//   double Mjj, WmT, DRlj1, DRlj2, DEtajj, Ptjj;
  float JetPFCor_Pt[8], JetPFCor_Eta[8], JetPFCor_Phi[8], JetPFCor_dphiMET[8], JetPFCor_bDiscriminatorCSV[8];
  float W_mt, TopWm, W_pt, event_met_pfmet;
  float W_electron_pt, W_electron_eta, W_muon_pt, W_muon_eta, ptlvjj;
  float genWWptresumwtnominal, genWWptresumwtresummup, genWWptresumwtresummdn, genWWptresumwtrenormup, genWWptresumwtrenormdn;
  int event_evtNo, vbf_event, W_electron_charge, W_muon_charge ;



  float W_V_pt_gen[6];
  bool W_V_hadronic_gen[6];
  float W_electron_pt_gen, W_electron_eta_gen, W_muon_pt_gen, W_muon_eta_gen;

  float effwt, puwt;
  float Mass2j_PFCor;

  //Read the input tree
  TFile* infile = new TFile(inDirName+inFileName);
  TTree* InTree = (TTree*)infile->Get("WJet");


//   EventNtuple* InEvtTree = new EventNtuple();
//   InTree->SetBranchAddress("EvtNtuple", &InEvtTree);
  //Create the output tree
  TFile *outfile = new TFile(outDirName+outFileName, "RECREATE");
  TTree* OutTree = new TTree("WJet", "Filtered Output Tree");
//   EventNtuple * OutEvtTree = new EventNtuple();
//   OutTree->Branch("EvtNtuple", "EventNtuple", &OutEvtTree);
//Set The Branches
  InTree->SetBranchAddress("JetPFCor_Pt",&JetPFCor_Pt);
  OutTree->Branch("JetPFCor_Pt",&JetPFCor_Pt,"JetPFCor_Pt[8]/F");
  InTree->SetBranchAddress("JetPFCor_Eta",&JetPFCor_Eta);
  OutTree->Branch("JetPFCor_Eta",&JetPFCor_Eta,"JetPFCor_Eta[8]/F");
  InTree->SetBranchAddress("JetPFCor_Phi",&JetPFCor_Phi);
  OutTree->Branch("JetPFCor_Phi",&JetPFCor_Phi,"JetPFCor_Phi[8]/F");
  InTree->SetBranchAddress("JetPFCor_dphiMET",&JetPFCor_dphiMET);
  OutTree->Branch("JetPFCor_dphiMET",&JetPFCor_dphiMET,"JetPFCor_dphiMET[8]/F");
  InTree->SetBranchAddress("JetPFCor_bDiscriminatorCSV",&JetPFCor_bDiscriminatorCSV);
  OutTree->Branch("JetPFCor_bDiscriminatorCSV",&JetPFCor_bDiscriminatorCSV,"JetPFCor_bDiscriminatorCSV[8]/F");

  InTree->SetBranchAddress("W_mt",&W_mt);
  OutTree->Branch("W_mt",&W_mt,"W_mt/F");
  InTree->SetBranchAddress("TopWm",&TopWm);
  OutTree->Branch("TopWm",&TopWm,"TopWm/F");
  InTree->SetBranchAddress("W_pt",&W_pt);
  OutTree->Branch("W_pt",&W_pt,"W_pt/F");
  InTree->SetBranchAddress("event_met_pfmet",&event_met_pfmet);
  OutTree->Branch("event_met_pfmet",&event_met_pfmet,"event_met_pfmet/F");
  InTree->SetBranchAddress("W_electron_pt",&W_electron_pt);
  OutTree->Branch("W_electron_pt",&W_electron_pt,"W_electron_pt/F");
  InTree->SetBranchAddress("W_electron_charge",&W_electron_charge);
  OutTree->Branch("W_electron_charge",&W_electron_charge,"W_electron_charge/I");
  InTree->SetBranchAddress("W_electron_eta",&W_electron_eta);
  OutTree->Branch("W_electron_eta",&W_electron_eta,"W_electron_eta/F");
  InTree->SetBranchAddress("W_muon_pt",&W_muon_pt);
  OutTree->Branch("W_muon_pt",&W_muon_pt,"W_muon_pt/F");
  InTree->SetBranchAddress("W_muon_charge",&W_muon_charge);
  OutTree->Branch("W_muon_charge",&W_muon_charge,"W_muon_charge/I");
  InTree->SetBranchAddress("W_muon_eta",&W_muon_eta);
  OutTree->Branch("W_muon_eta",&W_muon_eta,"W_muon_eta/F");
  InTree->SetBranchAddress("ptlvjj",&ptlvjj);
  OutTree->Branch("ptlvjj",&ptlvjj,"ptlvjj/F");

  InTree->SetBranchAddress("event_evtNo",&event_evtNo);
  OutTree->Branch("event_evtNo",&event_evtNo,"event_evtNo/I");
  InTree->SetBranchAddress("vbf_event",&vbf_event);
  OutTree->Branch("vbf_event",&vbf_event,"vbf_event/I");

  InTree->SetBranchAddress("effwt",&effwt);
  OutTree->Branch("effwt",&effwt,"effwt/F");
  InTree->SetBranchAddress("puwt",&puwt);
  OutTree->Branch("puwt",&puwt,"puwt/F");
  InTree->SetBranchAddress("Mass2j_PFCor",&Mass2j_PFCor);
  OutTree->Branch("Mass2j_PFCor",&Mass2j_PFCor,"Mass2j_PFCor/F");

  //MCTruth:
  InTree->SetBranchAddress("W_V_pt_gen[6]",&W_V_pt_gen);
  OutTree->Branch("W_V_pt_gen[6]",&W_V_pt_gen,"W_V_pt_gen[6]/F");
  InTree->SetBranchAddress("W_V_hadronic_gen[6]",&W_V_hadronic_gen);
  OutTree->Branch("W_V_hadronic_gen[6]",&W_V_hadronic_gen,"W_V_hadronic_gen[6]/O");

  InTree->SetBranchAddress("W_electron_pt_gen",&W_electron_pt_gen);
  OutTree->Branch("W_electron_pt_gen",&W_electron_pt_gen,"W_electron_pt_gen/F");
  InTree->SetBranchAddress("W_electron_eta_gen",&W_electron_eta_gen);
  OutTree->Branch("W_electron_eta_gen",&W_electron_eta_gen,"W_electron_eta_gen/F");
  InTree->SetBranchAddress("W_muon_pt_gen",&W_muon_pt_gen);
  OutTree->Branch("W_muon_pt_gen",&W_muon_pt_gen,"W_muon_pt_gen/F");
  InTree->SetBranchAddress("W_muon_eta_gen",&W_muon_eta_gen);
  OutTree->Branch("W_muon_eta_gen",&W_muon_eta_gen,"W_muon_eta_gen/F");

  //Reweighing:
  InTree->SetBranchAddress("genWWptresumwtnominal",&genWWptresumwtnominal);
  InTree->SetBranchAddress("genWWptresumwtresummup",&genWWptresumwtresummup);
  InTree->SetBranchAddress("genWWptresumwtresummdn",&genWWptresumwtresummdn);
  InTree->SetBranchAddress("genWWptresumwtrenormup",&genWWptresumwtrenormup);
  InTree->SetBranchAddress("genWWptresumwtrenormdn",&genWWptresumwtrenormdn);



  int nEntries=InTree->GetEntries();
  cout << "nEntries=" << nEntries << endl;


  for (Int_t i=0; i<nEntries; i++) {
    InTree->GetEntry(i);
    //cout << "W_V_pt_gen[0]=" << W_V_pt_gen[0] << endl;

    switch (reweighWWpT) {
    case 0:
      //Nominal
      effwt*=genWWptresumwtnominal;
      break;
    case +1:
      //ME-PS Scale Up
      effwt*=genWWptresumwtresummup;
      break;
    case -1:
      //ME-PS Scale Down
      effwt*=genWWptresumwtresummdn;
      break;
    case +2:
      //Factorization-Renormalization Scale Up
      effwt*=genWWptresumwtrenormup;
      break;
    case -2:
      //Factorization-Renormalization Scale Down
      effwt*=genWWptresumwtrenormdn;
      break;
    default:
      //No Reweighting
      break;
    }

    if ( passAllCuts ) {
      OutTree->Fill();
    }

  }

  outfile->Write();
  outfile->Close();

}


// *Br 1645 :genWWpt   : genWWpt/F                                              *
// *Entries :    69733 : Total  Size=     281866 bytes  File Size  =     256992 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.09     *
// *............................................................................*
// *Br 1646 :genWWptresumwtnominal : genWWptresumwtnominal/F                    *
// *Entries :    69733 : Total  Size=     282314 bytes  File Size  =     232821 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.21     *
// *............................................................................*
// *Br 1647 :genWWptresumwtresummdn : genWWptresumwtresummdn/F                  *
// *Entries :    69733 : Total  Size=     282346 bytes  File Size  =     234458 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.20     *
// *............................................................................*
// *Br 1648 :genWWptresumwtresummup : genWWptresumwtresummup/F                  *
// *Entries :    69733 : Total  Size=     282346 bytes  File Size  =     190304 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.48     *
// *............................................................................*
// *Br 1649 :genWWptresumwtrenormdn : genWWptresumwtrenormdn/F                  *
// *Entries :    69733 : Total  Size=     282346 bytes  File Size  =     229638 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.23     *
// *............................................................................*
// *Br 1650 :genWWptresumwtrenormup : genWWptresumwtrenormup/F                  *
// *Entries :    69733 : Total  Size=     282346 bytes  File Size  =     235312 *
// *Baskets :       28 : Basket Size=       8704 bytes  Compression=   1.20     *
// *............................................................................*


































// ///////////////////////////////////////////////////////////////////////////////////////
// ////// Functions to shift the pT (and associated values) by +/- JES Uncertainty
// ///////////////////////////////////////////////////////////////////////////////////////
// void ShiftDueToJESUncertainty(const char* inFileName, const char* outFileName, bool ShiftUp)
// //// Shift the pT by the uncertainty JES corrections; either up 1sigma (ShiftUp=true) or down 1sigma (ShiftUp=false).
// //// Load gSystem->Load("~ilyao/MATRIXELEMENT/CMSSW_4_2_8/lib/slc5_amd64_gcc434/libCondFormatsJetMETObjects.so") before compiling
// {

//   //double in_METEt, out_METEt;
//   vector < TLorentzVector > in_jLV, out_jLV, in_lLV;
//   TLorentzVector cjLV;
//   double pTUncorr, pTCorr, pXCorr, pYCorr, ECorr;
//   double out_Mjj, out_DRlj1, out_DRlj2;
//   double uncert, j_eta;

//   //read the input tree
//   TFile* infile = new TFile(inFileName);
//   TTree* InTree = (TTree*)infile->Get("EvtTree");
//   EventNtuple* InEvtTree = new EventNtuple();
//   InTree->SetBranchAddress("EvtNtuple", &InEvtTree);
//   // create the output tree
//   TFile *outfile = new TFile(outFileName, "RECREATE");
//   TTree* OutTree = new TTree("EvtTree", "Output tree for matrix element");
//   EventNtuple * OutEvtTree = new EventNtuple();
//   OutTree->Branch("EvtNtuple", "EventNtuple", &OutEvtTree);


//   int nEntries;
//   nEntries=InTree->GetEntries();
//   cout << "nEntries=" << nEntries << endl;

//   edm::FileInPath fip("mytestV13_AK5PF_Uncertainty.txt");
//   JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fip.fullPath());

//   for (Int_t i=0; i<nEntries; i++) {
//     InTree->GetEntry(i);
//     /// Record the branches being altered
//     //in_METEt=InEvtTree->METEt;
// //     cout << "in_METEt=" << in_METEt << endl;
// //     out_METEt=in_METEt+10.0;
// //     OutEvtTree->METEt=out_METEt;

//     in_jLV=InEvtTree->jLV;
//     in_lLV=InEvtTree->lLV;
//     out_jLV.clear();
//     OutEvtTree->lLV.clear();

//     for (Int_t j=0; j<2; j++) {
//       //cout << "jLV[" << j << "].Pt()=" << in_jLV[j].Pt() << endl;
//       pTUncorr=in_jLV[j].Pt();
//       j_eta=in_jLV[j].Eta();
//       jecUnc->setJetEta(j_eta); // Give rapidity of jet you want uncertainty on
//       jecUnc->setJetPt(pTUncorr);// Also give the corrected pt of the jet you want the uncertainty on
//       uncert = jecUnc->getUncertainty(true); // In principle, boolean controls if uncertainty on +ve or -ve side is returned (asymmetric errors) but not yet implemented.
//       //cout << "pt=" << pTUncorr << " eta=" << j_eta << " : uncert=" << uncert << endl;
//       if ( ShiftUp ) {
// 	pTCorr=pTUncorr*(1.0+uncert);
//       } else {
// 	pTCorr=pTUncorr*(1.0-uncert);
//       }
//       pXCorr=(pTCorr/pTUncorr)*in_jLV[j].Px();
//       pYCorr=(pTCorr/pTUncorr)*in_jLV[j].Py();
//       ECorr=sqrt(pTCorr*pTCorr+in_jLV[j].Pz()*in_jLV[j].Pz()+in_jLV[j].M()*in_jLV[j].M());
//       cjLV.SetPxPyPzE(pXCorr,pYCorr,in_jLV[j].Pz(),ECorr);
//       //OutEvtNtuple->lLV.push_back(cjLV);
//       out_jLV.push_back(cjLV);
//     }
//     OutEvtTree->jLV=out_jLV;
//     out_Mjj=(out_jLV[0]+out_jLV[1]).M();
//     OutEvtTree->Mjj=out_Mjj;
//     out_DRlj1=in_lLV[0].DeltaR(out_jLV[0]);
//     out_DRlj2=in_lLV[0].DeltaR(out_jLV[1]);
//     OutEvtTree->DRlj1=out_DRlj1;
//     OutEvtTree->DRlj2=out_DRlj2;

//     /// Record the unaltered branches
//     //OutEvtTree->jLV=InEvtTree->jLV;
//     OutEvtTree->METLV=InEvtTree->METLV; 
//     OutEvtTree->lLV=InEvtTree->lLV;
//     OutEvtTree->matchedGenParticles=InEvtTree->matchedGenParticles;
//     OutEvtTree->matchedpdgId=InEvtTree->matchedpdgId;
//     OutEvtTree->matchedDeltaR=InEvtTree->matchedDeltaR;
//     OutEvtTree->jBtag=InEvtTree->jBtag;
//     OutEvtTree->lQ=InEvtTree->lQ;
//     OutEvtTree->ldetComp=InEvtTree->ldetComp;
//     OutEvtTree->run=InEvtTree->run;
//     OutEvtTree->event=InEvtTree->event;
//     //OutEvtTree->Mjj=InEvtTree->Mjj;
//     OutEvtTree->leptonCat=InEvtTree->leptonCat;
//     OutEvtTree->leptonCat_passAll=InEvtTree->leptonCat_passAll;
//     //OutEvtTree->DRlj1=InEvtTree->DRlj1;
//     //OutEvtTree->DRlj2=InEvtTree->DRlj2;
//     OutEvtTree->Thetalj1pj2=InEvtTree->Thetalj1pj2;
//     OutEvtTree->lTotIso=InEvtTree->lTotIso;
//     OutEvtTree->lecalIso=InEvtTree->lecalIso;
//     OutEvtTree->lhcalIso=InEvtTree->lhcalIso;
//     OutEvtTree->ltrkIso=InEvtTree->ltrkIso;
//     OutEvtTree->METEt=InEvtTree->METEt;
//     OutEvtTree->lPhi=InEvtTree->lPhi;
//     //Fill The Output Tree 
//     OutTree->Fill();

//   }

// //   outfile->cd();
// //   OutTree->Write("EvtTree");
//   outfile->Write();
//   outfile->Close();

// }

// void ShiftJESMultiplePATSets(const char* inputmode, const char* targetDir, int StartProcess, int EndProcess) 
// //// Generate the shifted due to JES up and down files. The input the form targetDir+PLabel+inputmode+"_outfile.root", and the output is of the form targetDir+PLabel+inputmode+"_JESp1s.root", targetDir+PLabel+inputmode+"_JESm1s.root".
// {
//   TString infilename, outfilename;
//   InitializeLabels(PLabel,CLabel);

//   for (Int_t np=StartProcess; np<(EndProcess+1);np++) {
//     cout << "Process=" << PLabel[np] << endl;
//     infilename="_outfile.root";
//     infilename=inputmode+infilename;
//     infilename=PLabel[np]+infilename;
//     infilename=targetDir+infilename;
//     outfilename=inputmode;
//     outfilename=PLabel[np]+outfilename;
//     outfilename=targetDir+outfilename;

//     ShiftDueToJESUncertainty(infilename,outfilename+"_JESp1s.root",true);
//     ShiftDueToJESUncertainty(infilename,outfilename+"_JESm1s.root",false);
//   }
  
// }
