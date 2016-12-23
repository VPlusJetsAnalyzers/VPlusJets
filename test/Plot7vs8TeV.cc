//=====================================================================-*-C++-*-
// File and Version Information:
//      $Id: RooUnfoldExample.cxx 279 2011-02-11 18:23:44Z T.J.Adye $
//
// Description:
//      Simple example usage of the RooUnfold package using toy MC.
//
// Authors: Tim Adye <T.J.Adye@rl.ac.uk> and Fergus Wilson <fwilson@slac.stanford.edu>
//
//==============================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TString.h"
#include "TColor.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEventList.h"
#include "TLatex.h"
#include "TObjArray.h"
#include "TFractionFitter.h"

#endif

//==============================================================================
// Global definitions
//==============================================================================
//bool VerboseOutput = false;
const char* TreeName = "WJet";
const char* VarLabel = "dijet p_{T}";
const char* VaribaleName = "sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))";
const double varMin = 75;
const double varMax = 305-3*10;
const int numBins = 23-3;//10GeV bins
TString Cuts_8TeVData = "( (sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)&&(Mass2j_PFCor>48.000)&&(Mass2j_PFCor<160.000) )";

TString datafilename = "root://cmseos:1094//eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root";
TString PlotDir = "./";
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////*********Helper Functions*********//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void labelAHistogram(TH1D* &h, TString plotTitle, TString xaxisTitle, TString yaxisTitle, double plotTitleX = -1.0, double yaxisOffset = -1.0)
//// Add a title and axis labels to a histogram (to be drawn)
{
  h->SetTitle(plotTitle);
  gStyle->SetTitleY(0.97);
  if ( plotTitleX > 0 ) {
    gStyle->SetTitleX(plotTitleX);
  }
  h->SetXTitle(xaxisTitle);
  h->SetYTitle(yaxisTitle);
  if ( yaxisOffset > 0 ) {
    h->GetYaxis()->SetTitleOffset(yaxisOffset);
  }

}


// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void fill8TeVHist(TH1D* &h, TString evtWt = "effwt*puwt", const char* VarName = VaribaleName, const char* treeName = TreeName)
////Fill the data histogram
{

  TString Restrictions, histload_str;
  char Max_char[30], Min_char[30];

  /// Make a string defining the restrictions (Range + Additional Ones).
  
  sprintf(Max_char,"%.2e",varMax);
  sprintf(Min_char,"%.2e",varMin);
  Restrictions = ") )*"+evtWt;
  Restrictions = Cuts_8TeVData + Restrictions;
  Restrictions=")&&("+Restrictions;
  Restrictions=Max_char+Restrictions;
  Restrictions="<"+Restrictions;
  Restrictions=VarName+Restrictions;
  Restrictions=")&&("+Restrictions;
  Restrictions=Min_char+Restrictions;
  Restrictions=">"+Restrictions;
  Restrictions=VarName+Restrictions;
  Restrictions="( ("+Restrictions;

  cout << "Restrictions=" << Restrictions << endl;

  //Open the files
  TFile* f;
  TTree* InTree;
  f = new TFile(datafilename);
  InTree = (TTree*)f->Get(treeName);
  
  //Fill the histogram
  histload_str=">>htemp";
  histload_str=VarName+histload_str;
//   if ( VerboseOutput ) {
//     cout << "hist input =" << histload_str << endl;
//   }
    
  TH1D* htemp = new TH1D("htemp","htemp",numBins,varMin,varMax);
  InTree->Draw(histload_str,Restrictions);
  h = (TH1D*)htemp->Clone();
//   if ( dosumW2 ) {
//     h->Sumw2();
//   }
  //  h->Scale(scale);
  cout << "NEntries=" << h->GetEntries() << " , Integral=" << h->Integral() << endl; 
  delete htemp;
  
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void fill7TeVDataHist(TH1D* &h)
//// Fill the 7TeV Histogram
{

  double binVal[23] = {10300,8700,7100,5600,4400,3400,2650,2050,1550,1250,1000,800,600, 550,450,380,300,250,200,160,150,100,70}; //80-200 GeV
  double binErr = 118;// Set to the width of the point on the control plots

  ///Add the error for each bin in quadrature
  for (Int_t n=1; n<(numBins+1);++n) {
    h->SetBinContent(n,binVal[n-1]);
    h->SetBinError(n,binErr);
  }
    //    binErr=(h_FitErr->GetBinContent(n))*(h_FitErr->GetBinContent(n));
//     binErr=binErr+(fitFracErr_VJets[nChan]*h_VJets->GetBinContent(n))*(fitFracErr_VJets[nChan]*h_VJets->GetBinContent(n));
//     binErr=binErr+(fitFracErr_Top[nChan]*h_Top->GetBinContent(n))*(fitFracErr_Top[nChan]*h_Top->GetBinContent(n));
//     if ( gc_QCD[nChan]>0 ) { binErr=binErr+(fitFracErr_QCD[nChan]*h_QCD->GetBinContent(n))*(fitFracErr_QCD[nChan]*h_QCD->GetBinContent(n)); }

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////********* Main Implementation*********////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void Plot7vs8TeV(const char* SavePrefix = 0)
{
  double val_Ratio, val_8TeV, val_7TeV;
  double err_Ratio = 1000;
  /// Graphics Options:
  gStyle->SetOptStat(0);
  //  gROOT->ProcessLine(".L tdrstyle.C");

  /// Define the histograms
  TH1D *h_7TeV, *h_8TeV, *h_Ratio;
  h_7TeV = new TH1D("h_7TeV","7TeV Data",numBins,varMin,varMax);
  //  h_7TeV->Sumw2();
  //  h_Diboson_Comb->SetFillColor(color_Diboson);
  h_8TeV = new TH1D("h_8TeV","8TeV Data",numBins,varMin,varMax);
  h_Ratio = new TH1D("h_Ratio","8/7 TeV Data",numBins,varMin,varMax);
  //  h_Subtracted_Comb->SetLineWidth(2);
  //  h_Subtracted_Comb->SetMarkerStyle(8);

  /// Fill the histograms
  fill8TeVHist(h_8TeV);
  fill7TeVDataHist(h_7TeV);

  for (Int_t n=1; n<(numBins+1);++n) {
    val_8TeV=h_8TeV->GetBinContent(n);
    cout << "bin " << n << ", number of 8 TeV events=" << val_8TeV << endl;
    val_7TeV=h_7TeV->GetBinContent(n);
    val_Ratio=val_8TeV/val_7TeV;
    err_Ratio=(sqrt(val_7TeV*val_7TeV*( h_8TeV->GetBinError(n))*(h_8TeV->GetBinError(n))+val_8TeV*val_8TeV*(h_7TeV->GetBinError(n))*(h_7TeV->GetBinError(n)) ))/(val_7TeV*val_7TeV);

    h_Ratio->SetBinContent(n,val_Ratio);
    h_Ratio->SetBinError(n,err_Ratio);
  }

  /// Plot the histograms
  TCanvas *cnv_8TeV, *cnv_7TeV, *cnv_Ratio;
  cnv_8TeV = new TCanvas("cnv_8TeV","8TeV",10,10,900,600);
  cnv_7TeV = new TCanvas("cnv_7TeV","7TeV",10,10,900,600);
  cnv_Ratio = new TCanvas("cnv_Ratio","Ratio",10,10,900,600);

  ///Make Comparison Plots for the current channel
  cout << "Plotting... " << endl;
  cnv_8TeV->cd();
  cnv_8TeV->Clear();
    // TLegend *lgnd1 = new TLegend(0.7,0.65,0.90,0.9);
    // lgnd1->AddEntry(h_Diboson,"WW+WZ Truth","l");
  labelAHistogram(h_8TeV,"8 TeV Data",VarLabel,"8TeV Event Count",0.3);
  h_8TeV->SetMarkerStyle(8);
  h_8TeV->Draw("ep");
    // labelAHistogram(h_Diboson_MCTruth,"Expected Diboson MC With vs Without Truth : " + chanTitleLabel[chan],VarLabel,"Expected Event Count",0.12);
    // lgnd1->AddEntry(h_Diboson,"WW+WZ Default","f");
    // h_Diboson->Draw("hist same");
    // h_Diboson_MCTruth->SetMaximum(1.1*h_Diboson->GetMaximum());
    // h_Diboson_MCTruth->Draw("hist same");
    // lgnd1->Draw();
  cnv_8TeV->Update();

  cnv_7TeV->cd();
  cnv_7TeV->Clear();
  labelAHistogram(h_7TeV,"7 TeV Data",VarLabel,"7TeV Event Count",0.3);
  h_7TeV->SetMarkerStyle(8);
  h_7TeV->Draw("ep");
  cnv_7TeV->Update();

  cnv_Ratio->cd();
  cnv_Ratio->Clear();
  labelAHistogram(h_Ratio,"Ratio of 8 to 7 TeV Data",VarLabel,"8TeV/7TeV Event Count",0.3);
  h_Ratio->SetMarkerStyle(8);
  h_Ratio->Draw("ep");
  cnv_Ratio->Update();



  /// Save histograms for this channel
  if ( SavePrefix!=0 ) {
    cnv_8TeV->SaveAs(PlotDir+SavePrefix+"_8TeV.png");
    cnv_7TeV->SaveAs(PlotDir+SavePrefix+"_7TeV.png");
    cnv_Ratio->SaveAs(PlotDir+SavePrefix+"_Ratio.png");
  }


}

