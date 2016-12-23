#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixT.h>
#include <TMatrixTSym.h>

#include "TMath.h"
#include "TF1.h"
#include "TMinuit.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TPaveText.h"

#include "CMS_lumi.h"
#include "Blue.h"

//#include "/uscms_data/d3/ilyao/CMSSW_5_3_15_patch1/src/ElectroWeakAnalysis/VPlusJets/test/inputDibosonMCInfoFitAndValidationResults.cpp"
#include "../inputDibosonMCInfoFitAndValidationResults.cpp"
using namespace dibosonAnalysis;
using namespace std;


////////// ******* Global Parameters ********* ////////////
const int NCHANNELS=6;
const double LumiUnc = 0.026; // Included as an external systematic error (fully correlated between the channels)
const bool RunBLUEComb = true;
const bool VerboseOutput = true;


// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void cmspre()
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.9,0.93,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.65,0.93,Form("#scale[0.5]{#lower[-0.15]{#it{#int}}}#it{L} dt = %0.0f#kern[0.2]{fb}^{-1}", 19.3));
  latex.SetTextAlign(11); // align left
//  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");
  latex.DrawLatex(0.15,0.93,"CMS");

}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void 
CMS_lumi( TPad* pad, int iPeriod, int iPosX )
{            
  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 ) 
    {
      lumiText = lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }
   
  std::cout << lumiText << std::endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont(cmsTextFont);
	  latex.SetTextSize(cmsTextSize*t);
	  latex.SetTextAlign(align_);
	  latex.DrawLatex(posX_, posY_, cmsText);
	  if( writeExtraText ) 
	    {
	      latex.SetTextFont(extraTextFont);
	      latex.SetTextAlign(align_);
	      latex.SetTextSize(extraTextSize*t);
	      latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	    }
	}
    }
  else if( writeExtraText )
    {
      if( iPosX==0) 
	{
	  posX_ =   l +  relPosX*(1-l-r);
	  posY_ =   1-t+lumiTextOffset*t;
	}
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);      
    }
  return;
}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//

int plot6Chan_WW_WZ_crossX_8TeV(double cY_Val[6], double cY_Err[6], double cY_ErrStat[6], double cYTot_Val, double cYTot_Err, double cYTot_ErrStat) {

// ------- Setup the canvas ------- 
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetEndErrorSize(5);

  // gStyle->SetPadRightMargin(0.2);
  // gStyle->SetPadBottomMargin(0.3);
  // gStyle->SetErrorX(0.5);


  TCanvas *c1 = new TCanvas( "c1", " Cross Section WW+WZ ", 200, 10, 1000, 1000);
  //  TPad*    c1_1 =(TPad*)c1->GetPrimitive("c1_1");
  c1->cd(1);
  c1->SetTicks(1,1);
  
  /// Define the y values and arrays for y values
  double cYTheory_Val = 82.7; // 59.8+22.9 from arXiv:1408.5243v1 [hep-ph]
  double cYTheory_Err = 2.5; // SMP-13-008 (+2.51/-1.74)

  double y[8] = {1., 2., 3., 4., 5., 6., 7., 8.};
  double yExtra[10] = {1., 2., 3., 4., 5., 6., 7., 8. ,9. ,10.};
  double ey[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
  double eyBand[8] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  double eyBandExtra[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

  double frameVal[10], frameErr[10];
  double theoryVal[8], theoryErr[8];
  double combineVal[8], combineErr[8];
  double Result[8], EResult[8], EResultStat[8];

  /// Fill the y value arrays
  double plotMin=10.0;
  double plotMax=350;
  for (Int_t i=0; i<10; i++) {
    frameVal[i] = (plotMin+plotMax)/2.0;
    frameErr[i] = (plotMax-plotMin)/2.0;
  }
  for (Int_t i=0; i<8; i++) {
    theoryVal[i] = cYTheory_Val;
    theoryErr[i] = cYTheory_Err;
  }
  for (Int_t i=0; i<8; i++) {
    combineVal[i] = cYTot_Val;
    combineErr[i] = cYTot_Err;
  }

  for (Int_t i=0; i<6; i++) {
    Result[7-i]=cY_Val[i];
    EResult[7-i]=cY_Err[i];
    EResultStat[7-i]=cY_ErrStat[i];
  }
  Result[1] = cYTot_Val;    
  EResult[1] = cYTot_Err; 
  EResultStat[1] = cYTot_ErrStat; 
  Result[0] = cYTheory_Val;
  EResult[0] = cYTheory_Err;
  EResultStat[0] = cYTheory_Err;


  /// Construct the individual graphs and combine
  ///Global features:
  TGraphErrors  *frameBand = new TGraphErrors(10, frameVal, yExtra, frameErr, eyBandExtra);
  frameBand->SetTitle();
  frameBand->GetXaxis()->SetLabelSize(0.04);
  frameBand->GetXaxis()->SetTitleSize(0.05);
  frameBand->GetXaxis()->SetTitle("#sigma(pb)");
  frameBand->GetXaxis()->SetTitleOffset(0.85);
  frameBand->GetXaxis()->SetNdivisions(510);//10 primary, 5 secondary

  frameBand->GetYaxis()->SetLabelSize(0.0);//No Labels
  frameBand->GetYaxis()->SetTitle("channel");
  frameBand->GetYaxis()->SetTitleOffset(0.5);
  frameBand->GetYaxis()->SetTitleSize(0.05);
  frameBand->GetYaxis()->SetNdivisions(20);
  //frameBand -> GetYaxis() -> CenterTitle();
  frameBand->SetFillColor(0);
  frameBand->SetFillStyle(1001);

  ///Individual Graphs:
  TGraphErrors  *theoryBand = new TGraphErrors(8, theoryVal, y, theoryErr, eyBand);
  theoryBand->SetFillColor(kGreen+2);
  theoryBand->SetFillStyle(1001);

  TGraphErrors  *combineBand = new TGraphErrors(8, combineVal, y, combineErr, eyBand);
  combineBand->SetFillColor(kCyan+2);
  combineBand->SetFillStyle(1001);

  TGraphErrors  *resultGraphStat = new TGraphErrors(8, Result, y, EResultStat, ey);
  resultGraphStat->SetLineWidth(2);

  TGraphErrors  *resultGraphBar = new TGraphErrors(8, Result, y, EResult, ey);//Redundant graph used to ensure that the vertical ends of the error bars aren't plotted as dashed
  resultGraphBar->SetLineWidth(2);

  TGraphErrors  *resultGraph = new TGraphErrors(8, Result, y, EResult, ey);
  resultGraph->SetMarkerColor(kRed);
  resultGraph->SetMarkerStyle(20);
  resultGraph->SetMarkerSize(1.2);
  resultGraph->SetLineWidth(2);
  resultGraph->SetLineStyle(2);

  /// Combined Graph
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(frameBand,"a2");
  mg->Add(theoryBand,"2");
  mg->Add(combineBand,"2");
  mg->Add(resultGraphStat,"P");
  mg->Add(resultGraphBar,"||");
  mg->Add(resultGraph,"P");
  mg->Draw();

  /// Labeling
  TLatex l_Title;
  l_Title.SetTextSize(0.042);
  l_Title.DrawLatex(60,10.4, "#sigma [pp #rightarrow W(#rightarrowl#nu) W/Z(#rightarrowq#bar{q})]:");
  TLatex l_Combine;
  char text1[40];
  sprintf(text1,"Combined #sigma = %5.1f #pm %4.1f pb (Cyan)",cYTot_Val,cYTot_Err);
  l_Combine.SetTextSize(0.031);
  l_Combine.DrawLatex(75,9.8,text1);
  TLatex l_Theory;
  char text2[40];
  sprintf(text2,"Theory (NLO) = %4.1f #pm %3.1f pb (Green)",Result[0],EResult[0]);
  l_Theory.SetTextSize(0.031);
  l_Theory.DrawLatex(75,9.3,text2);


  double textXStart=290.0;
  double textSize=0.030;
  TLatex l_c1;
  l_c1.SetTextSize(textSize);
  l_c1.DrawLatex(textXStart,7.9, "anti-btag #mu");
  TLatex l_c2;
  l_c2.SetTextSize(textSize);
  l_c2.DrawLatex(textXStart,6.9, "anti-btag el"); 
  TLatex l_c3;
  l_c3.SetTextSize(textSize);
  l_c3.DrawLatex(textXStart,5.9, "btag #mu");
  TLatex l_c4;
  l_c4.SetTextSize(textSize);
  l_c4.DrawLatex(textXStart,4.9, "btag el");
  TLatex l_c5;
  l_c5.SetTextSize(textSize);
  l_c5.DrawLatex(textXStart,3.9, "merged #mu");
  TLatex l_c6;
  l_c6.SetTextSize(textSize);
  l_c6.DrawLatex(textXStart,2.9, "merged el"); 
  TLatex l_c7;
  l_c7.SetTextSize(textSize);
  l_c7.DrawLatex(textXStart,1.9, "combined");
  TLatex l_c8;
  l_c8.SetTextSize(textSize);
  l_c8.DrawLatex(textXStart,0.9, "theory"); 

  //cmspre(); 
  ///Load the lumi parameters
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19 fb^{-1}"; // default is "19.7 fb^{-1}"
  //  relPosX = 0.7; // Adjust the label postion
  cmsTextSize = 0.4;

  CMS_lumi(c1,2,0);
  c1->Update();
  c1->SaveAs("WW_WZ_crossX_6chan_8TeV.png");
  return 0;

}


// -----------------------------------------------------------------------------------------------------------------------------------------------------------//
void crossXcomb(bool plotCrossX, const char* vals, const char* errs_stat, const char* errs_syst, const char* errs_correlated, double errSystExt_All = sqrt(0.05*0.05+LumiUnc*LumiUnc+0.03*0.03), double errSystExt_Resolved = sqrt(0.08*0.08+0.05*0.05+0.12*0.12), double errSystExt_Merged = sqrt(0.03*0.03+0.32*0.32), const int nChan = NCHANNELS) {
//// Currently the total syst error is a combination of documented sources + WWpt reweighting error + aMC@NLO vs PowHeg AxEff
//// Input the cross-section values, statistical errors and systematic errors
//// The resolved and boosted channels are handled separately and the external systematics are then included assuming 100% correlations between the channels. 
  bool printCorrMatrices=true;

  double Val[nChan];
  double Err[nChan];
  double Weight[nChan], WeightUncorr[nChan];
  double Err_stat[nChan];
  double Err_syst[nChan];
  double Err_correlated[nChan];
  TMatrixDSym Cov(nChan), CovInv(nChan);
  TMatrixTSym<double> CovUncorr(nChan);
  TMatrixTSym<double> CovUncorrInv(nChan);

  double weight_norm=0;
  double weightUncorr_norm=0;
  double errTot_stat=0;
  double errTot_syst=0;
  double errTot_lumi=0;
  double errTot_comb=0;
  double errTotUncorr_stat=0;
  double errTotUncorr_syst=0;
  double errTotUncorr_lumi=0;
  double errTotUncorr_comb=0;
  double crossXTot=0;
  double crossXTotUncorr=0;

  double weightTot=0;
  double weightTotUncorr=0;
  istrstream inVals(vals);
  istrstream inErr_stat(errs_stat);
  istrstream inErr_syst(errs_syst);
  istrstream inErr_correlated(errs_correlated);

  //Read in the values & errors
  for (Int_t i=0; i<nChan; i++) {
    inVals >> Val[i];
    inErr_stat >> Err_stat[i];
    inErr_syst >> Err_syst[i];
    inErr_correlated >> Err_correlated[i];
    Err[i]=sqrt(Err_stat[i]*Err_stat[i]+Err_syst[i]*Err_syst[i]);
  }
  ///Input the Covariance Matrix
  for (Int_t i=0; i<nChan; i++) {
    for (Int_t j=0; j<nChan; j++) {
      CovUncorr(i,j)=0;
      Cov(i,j)=errSystExt_All*errSystExt_All*Val[i]*Val[j]+Err_correlated[i]*Err_correlated[j];//Include the total errors (100% correlation among all channels)
      if ( (i<4)&&(j<4) ) { Cov(i,j)=Cov(i,j)+errSystExt_Resolved*errSystExt_Resolved*Val[i]*Val[j]; } //Add the resolved channel errors (100% correlated)
      if ( (i>=4)&&(j>=4) ) { Cov(i,j)=Cov(i,j)+errSystExt_Merged*errSystExt_Merged*Val[i]*Val[j]; } //Add the boosted channel errors (100% correlated)
      if ( i==j ) { 
	Cov(i,j)=Cov(i,j)+Err[i]*Err[i]; //Add the individual channel errors (uncorrelated)
	CovUncorr(i,j)=errSystExt_All*errSystExt_All*Val[i]*Val[j]+Err[i]*Err[i];
	if ( (i<4)&&(j<4) ) { CovUncorr(i,j)=CovUncorr(i,j)+errSystExt_Resolved*errSystExt_Resolved*Val[i]*Val[j]; }
	if ( (i>=4)&&(j>=4) ) { CovUncorr(i,j)=CovUncorr(i,j)+errSystExt_Merged*errSystExt_Merged*Val[i]*Val[j]; }
      }
    }
  }
  ///Invert the matrices and compute the weights
  CovInv=Cov;
  CovInv=CovInv.Invert();
  //CovInv.Print();
  CovUncorrInv = CovUncorr;
  CovUncorrInv = CovUncorrInv.Invert();
  //CovUncorrInv.Print();

  if ( printCorrMatrices ) {
    cout << "Covariance Matrix:"  << setprecision(1) << endl;
    for (Int_t i=0; i<nChan; i++) {
      for (Int_t j=0; j<nChan; j++) {
	cout << Cov(i,j) << " ";
      }
      cout << endl;
    }
    cout << "Inverse of the Covariance Matrix:" << setprecision(7) << endl;
    for (Int_t i=0; i<nChan; i++) {
      for (Int_t j=0; j<nChan; j++) {
	cout << CovInv(i,j) << " ";
      }
      cout << endl;
    }

    cout << "Uncorrelated Covariance Matrix:" << setprecision(1) << endl;
    for (Int_t i=0; i<nChan; i++) {
      for (Int_t j=0; j<nChan; j++) {
	cout << CovUncorr(i,j) << " ";
      }
      cout << endl;
    }
    cout << "Inverse of the Uncorrelated Covariance Matrix:" << setprecision(7) << endl;
    for (Int_t i=0; i<nChan; i++) {
      for (Int_t j=0; j<nChan; j++) {
	cout << CovUncorrInv(i,j) << " ";
      }
      cout << endl;
    }

  }

  for (Int_t i=0; i<nChan; i++) {
    for (Int_t j=0; j<nChan; j++) {
      weight_norm=weight_norm+CovInv(i,j);
      weightUncorr_norm=weightUncorr_norm+CovUncorrInv(i,j);
    }
  }

  for (Int_t i=0; i<nChan; i++) {
    Weight[i]=0;
    WeightUncorr[i]=0;
    for (Int_t j=0; j<nChan; j++) {
      Weight[i]=Weight[i]+CovInv(i,j)/weight_norm;
      WeightUncorr[i]=WeightUncorr[i]+CovUncorrInv(i,j)/weightUncorr_norm;
    }
  }

  ///Compute and display the cross-sections with errors
  for (Int_t i=0; i<nChan; i++) {
    cout  << setprecision(1) << fixed << "Channel " << i << " : " << endl;
    cout << "crossX = " << Val[i] << " : Error=" << sqrt(Cov(i,i)) << "=" << Err_stat[i] << "(stat)\\pm " << sqrt(Cov(i,i)-Err_stat[i]*Err_stat[i]-LumiUnc*LumiUnc*Val[i]*Val[i]) << "(syst)  :  UncorrErr=" << sqrt(CovUncorr(i,i)) << setprecision(3) << ", Weight=" << Weight[i] << ", UncorrWeight=" << WeightUncorr[i] << endl;
    crossXTot=crossXTot+Weight[i]*Val[i];
    crossXTotUncorr=crossXTotUncorr+WeightUncorr[i]*Val[i];
    errTot_stat=errTot_stat+Weight[i]*Weight[i]*Err_stat[i]*Err_stat[i];
    errTotUncorr_stat=errTotUncorr_stat+WeightUncorr[i]*WeightUncorr[i]*Err_stat[i]*Err_stat[i];
    for (Int_t j=0; j<nChan; j++) {
      errTot_comb=errTot_comb+Weight[i]*Weight[j]*Cov(i,j);
      errTotUncorr_comb=errTotUncorr_comb+WeightUncorr[i]*WeightUncorr[j]*CovUncorr(i,j);
    }
    //errTot_syst=errTot_syst+Weight[i]*Weight[i]*Err_syst[i]*Err_syst[i];
    weightTot=weightTot+Weight[i];
    weightTotUncorr=weightTotUncorr+WeightUncorr[i];
  }



  if (RunBLUEComb) { 
    //gSystem->Load("libBlue.so");
    const Int_t NUnc=3;
    Blue blc(nChan,NUnc);//NumEst=nChan, NumUnc=3 - Uncorrelated, FullyCorrelated, Correlated either among resolved or among boosted
    // Create arrays containing values and errors
    vector <Double_t *> chanPt; //runs over the number of channels
    for (Int_t s=0; s<nChan; s++) {
      chanPt.push_back( new Double_t[NUnc+1]() );
      chanPt[s][0]=Val[s];
      chanPt[s][1]=Err[s];
      chanPt[s][2]=errSystExt_All*Val[s]+Err_correlated[s];
      if ( s<4 ) {
	chanPt[s][3]=errSystExt_Resolved*Val[s];
      }
      if ( (s==4)||(s==5) ) {
	chanPt[s][3]=errSystExt_Merged*Val[s];
      }
    }

    //Add the correlations
    double MCorr[NUnc][nChan][nChan];
    // Initialize as identity
    for (Int_t r=0; r<NUnc; r++) {
      for (Int_t s=0; s<nChan; s++) {
	for (Int_t t=0; t<nChan; t++) {
	  MCorr[r][s][t]=0.0;
	  if ( s==t ) {  MCorr[r][s][t]=1.0; }
	}
      }
    }
    /// Add correlation terms
    // Uncertainty 0 corresponds to uncorrelated terms
    // Uncertainty 1 corresponds to fully correlated terms
    for (Int_t s=0; s<nChan; s++) {
      for (Int_t t=0; t<nChan; t++) {
	MCorr[1][s][t]=1.0;
      }
    }
    // Uncertainty 2 corresponds to terms which are correlated either among resolved or among boosted (but not between the two topologies)
    for (Int_t s=0; s<nChan; s++) {
      for (Int_t t=0; t<nChan; t++) {
	if ( (s<4)&&(t<4) ) {
	  //Resolved
	  MCorr[2][s][t]=1.0;
	}
	if ( ((s==4)||(s==5)) && ((t==4)||(t==5)) ) {
	  //Boosted
	  MCorr[2][s][t]=1.0;
	}
      }
    }


    /// Create a vector of flat arrays for input to BLUE
    vector <Double_t *> flatCorrsPt; //runs over the number of uncertainties
    for (Int_t r=0; r<NUnc; r++) {
      flatCorrsPt.push_back( new Double_t[nChan*nChan]() );
      for (Int_t s=0; s<nChan; s++) {
	for (Int_t t=0; t<nChan; t++) {
	  flatCorrsPt[r][s*nChan+t]=MCorr[r][s][t];
	  //cout << "(r,s,t)=(" << r << "," << s << "," << t << "), MCorr[r][s][t]=" << MCorr[r][s][t] << ", flatCorrsPt[r][s*nChan+t]=" << flatCorrsPt[r][s*nChan+t] << endl;
	}
      }
    }

    /// Examine the flattened channel parameters and correlations
    if (  VerboseOutput ) {
      cout << "Channel Parameters: " << endl;
      for (Int_t s=0; s<nChan; s++) {
	cout << "Channel " << s << ":";
	for (Int_t r=0; r<NUnc+1; r++) {
	  cout << " " << chanPt[s][r];
	}
	cout << endl;
      }

      cout << "Correlations: " << endl;
      for (Int_t r=0; r<NUnc; r++) {
	cout << "Uncertainty Type " << r << ":";
	for (Int_t st=0; st<nChan*nChan; st++) {
	  cout << " " << flatCorrsPt[r][st];
	}
	cout << endl;
      }
    }

    /// Fill the Blue
    for (Int_t s=0; s<nChan; s++) {
      blc.FillEst(s,chanPt[s]);
    }
    for (Int_t r=0; r<NUnc; r++) {
      blc.FillCor(r,flatCorrsPt[r]);
    }

    /// Solve
    blc.SetPrintLevel(1);
    blc.FixInp();

    //blc.Solve();
    // blc.PrintEst();
    // blc.PrintResult();

    blc.SolvePosWei();

    blc.PrintEst();
    blc.PrintResult();

    Double_t * tmp = new Double_t[NUnc+1]();
    blc.GetResult(tmp);
    //cout << "tmp[0]=" << tmp[0] << ", tmp[1]=" << tmp[1]  << ", tmp[3]=" << tmp[3]<<endl;
    crossXTot=tmp[0];
    errTot_comb=0;
    for (Int_t r=1; r<NUnc+1; r++) { errTot_comb+=tmp[r]*tmp[r]; }

  } 

  errTot_lumi = LumiUnc*LumiUnc*crossXTot*crossXTot;
  errTotUncorr_lumi = LumiUnc*LumiUnc*crossXTotUncorr*crossXTotUncorr;
  errTot_syst=sqrt(errTot_comb-errTot_stat-errTot_lumi);
  errTotUncorr_syst=sqrt(errTotUncorr_comb-errTotUncorr_stat-errTotUncorr_lumi);
  errTot_comb=sqrt(errTot_comb);
  errTotUncorr_comb=sqrt(errTotUncorr_comb);
  errTot_stat=sqrt(errTot_stat);
  errTotUncorr_stat=sqrt(errTotUncorr_stat);
  errTot_lumi=sqrt(errTot_lumi);
  errTotUncorr_lumi=sqrt(errTotUncorr_lumi);

  //  cout << "FitCrossX = " << crossXTot << "\\pm " << errTot_stat << "(stat.) \\pm " << errTot_syst << "(syst.)" << endl;
  cout << "WeightTot = " << weightTot << ", WeightTotUncorr = " << weightTotUncorr << " (should be 1)" << endl;
  cout << "TotalCrossX = " << setprecision(1) << fixed << crossXTot << " \\pm " << errTot_stat << "\\text{(stat)} \\pm " << errTot_syst << "\\text{(syst)} \\pm " << errTot_lumi << "\\text{(lumi)}" << endl;
  cout << "Without Correlations : TotalCrossX = " << crossXTotUncorr << " \\pm " << errTotUncorr_stat << "\\text{(stat)} \\pm " << errTotUncorr_syst << "\\text{(syst)} \\pm " << errTotUncorr_lumi << "\\text{(lumi)}" << endl;

  //  Add the external systematics as well as luminosity uncertainty to individual channel measurements
  for (Int_t i=0; i<nChan; i++) {
    Err[i]=sqrt(Cov(i,i));
  }

  if ( plotCrossX ) {
    plot6Chan_WW_WZ_crossX_8TeV(Val,Err,Err_stat,crossXTot,errTot_comb,errTot_stat);
  }
  

} 

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//


void processDibosonFitOutput_BLUE(){
//// Compute the various output summaries when given the expected number of signal (Diboson) events and fitted fractions & corresponding errors for each channel. The total fitted event counts in the data are also given in order to compute the stat vs syst uncertainties.
//// Compute the corrected diboson fraction, corresponding event yields with stat (=sqrt{NTotal}) and (fit) systematic errors, individual cross sections and the total cross-section
  
  ///Input Quantities
  TString ChannelName[NCHANNELS];
  double ExpectedYield[NCHANNELS];
  double FitFracVal[NCHANNELS];
  double FitFracErr[NCHANNELS];
  double DataEvts[NCHANNELS];
  double YieldBias[NCHANNELS];//Subtract from the fitted yield in order to correct
  double ErrorInflation[NCHANNELS];//Multiply the fit error by this factor in order to correct
  double QuadError[NCHANNELS];//Add in quadrature to the corrected fit error (in events) in order to correct
  double CorrelatedError[NCHANNELS];//Take to be correlated among the channels
  double Ae[NCHANNELS];
  double Lumi[NCHANNELS];//in pb^-1
  ///Computed Quantities
  double FitYieldVal[NCHANNELS];
  double FitYieldErr[NCHANNELS];
  double CorrYieldVal[NCHANNELS];
  double CorrYieldErr[NCHANNELS];
  double CorrFracVal[NCHANNELS];
  double CombFracErr[NCHANNELS];
  double YieldStatErr[NCHANNELS];
  double CorrYieldSystErr[NCHANNELS];//Internal/Individual channel systematic error defined as the difference between the corrected fit error and the statistical error (=sqrt(NTotal)) for that channel.
  double CrossXVal[NCHANNELS];
  double CrossXErr[NCHANNELS];
  double CrossXStatErr[NCHANNELS];
  double CrossXSystErr[NCHANNELS];
  double CrossXCorrelatedErr[NCHANNELS];
  char crossX_char[7];
  TString str_CrossXVal="";
  TString str_CrossXStatErr="";
  TString str_CrossXSystErr="";
  TString str_CrossXCorrelatedErr="";

  double FitDibosonYield=0;
  double FitDibosonErr=0;
  double CorrDibosonYield=0;
  double CorrDibosonErr=0;

  ///Input the parameters from the fit
  ChannelName[0]="muons_anti-btagged";
  ExpectedYield[0]=expectedProcessYields[mu_antibtagged][diboson];
  FitFracVal[0]=fitFracVal[mu_antibtagged][diboson];
  FitFracErr[0]=fitFracErr[mu_antibtagged][diboson];
  DataEvts[0]=expectedDataYield[mu_antibtagged];
  YieldBias[0]=signalYieldBias[mu_antibtagged];
  ErrorInflation[0]=signalErrorInflation[mu_antibtagged];
  CorrelatedError[0]=signalCorrelatedError[mu_antibtagged]; //Alternate W+Jets 
  QuadError[0]=signalQuadError[mu_antibtagged]; //MC Syst errors
  Ae[0]=signalAxEff[mu_antibtagged];
  Lumi[0]=lumiMuon;

  ChannelName[1]="electrons_anti-btagged";
  ExpectedYield[1]=expectedProcessYields[el_antibtagged][diboson];
  FitFracVal[1]=fitFracVal[el_antibtagged][diboson];
  FitFracErr[1]=fitFracErr[el_antibtagged][diboson];
  DataEvts[1]=expectedDataYield[el_antibtagged];
  YieldBias[1]=signalYieldBias[el_antibtagged];
  ErrorInflation[1]=signalErrorInflation[el_antibtagged];
  CorrelatedError[1]=signalCorrelatedError[el_antibtagged];
  QuadError[1]=signalQuadError[el_antibtagged];
  Ae[1]=signalAxEff[el_antibtagged];
  Lumi[1]=lumiElectron;

  ChannelName[2]="muons_btagged";
  ExpectedYield[2]=expectedProcessYields[mu_btagged][diboson];
  FitFracVal[2]=fitFracVal[mu_btagged][diboson];
  FitFracErr[2]=fitFracErr[mu_btagged][diboson];
  DataEvts[2]=expectedDataYield[mu_btagged];
  YieldBias[2]=signalYieldBias[mu_btagged];
  ErrorInflation[2]=signalErrorInflation[mu_btagged];
  CorrelatedError[2]=signalCorrelatedError[mu_btagged];
  QuadError[2]=signalQuadError[mu_btagged];
  Ae[2]=signalAxEff[mu_btagged];
  Lumi[2]=lumiMuon;

  ChannelName[3]="electrons_btagged";
  ExpectedYield[3]=expectedProcessYields[el_btagged][diboson];
  FitFracVal[3]=fitFracVal[el_btagged][diboson];
  FitFracErr[3]=fitFracErr[el_btagged][diboson];
  DataEvts[3]=expectedDataYield[el_btagged];
  YieldBias[3]=signalYieldBias[el_btagged];
  ErrorInflation[3]=signalErrorInflation[el_btagged];
  CorrelatedError[3]=signalCorrelatedError[el_btagged];
  QuadError[3]=signalQuadError[el_btagged];
  Ae[3]=signalAxEff[el_btagged];
  Lumi[3]=lumiElectron;

  ChannelName[4]="muons_boosted";
  ExpectedYield[4]=expectedProcessYields[mu_boosted][diboson];
  FitFracVal[4]=fitFracVal[mu_boosted][diboson];
  FitFracErr[4]=fitFracErr[mu_boosted][diboson];
  DataEvts[4]=expectedDataYield[mu_boosted];
  YieldBias[4]=signalYieldBias[mu_boosted];
  ErrorInflation[4]=signalErrorInflation[mu_boosted];
  CorrelatedError[4]=signalCorrelatedError[mu_boosted];
  QuadError[4]=signalQuadError[mu_boosted];
  Ae[4]=signalAxEff[mu_boosted];
  Lumi[4]=lumiMuon;

  ChannelName[5]="electrons_boosted";
  ExpectedYield[5]=expectedProcessYields[el_boosted][diboson];
  FitFracVal[5]=fitFracVal[el_boosted][diboson];
  FitFracErr[5]=fitFracErr[el_boosted][diboson];
  DataEvts[5]=expectedDataYield[el_boosted];
  YieldBias[5]=signalYieldBias[el_boosted];
  ErrorInflation[5]=signalErrorInflation[el_boosted];
  CorrelatedError[5]=signalCorrelatedError[el_boosted];
  QuadError[5]=signalQuadError[el_boosted];
  Ae[5]=signalAxEff[el_boosted];
  Lumi[5]=lumiElectron;

  bool statErrIsgtPartialTotErr = false;

  ///Compute Output Quantities
  for (Int_t i=0; i<NCHANNELS; i++) {
    FitYieldVal[i] = FitFracVal[i]*ExpectedYield[i];
    FitDibosonYield = FitDibosonYield+FitYieldVal[i];

    FitYieldErr[i] = FitFracErr[i]*ExpectedYield[i];
    FitDibosonErr = FitDibosonErr+FitYieldErr[i]*FitYieldErr[i];

    CorrYieldVal[i] = FitYieldVal[i]-YieldBias[i];
    CorrDibosonYield = CorrDibosonYield+CorrYieldVal[i];

    CorrYieldErr[i] = FitYieldErr[i]*ErrorInflation[i];
    CorrYieldErr[i] = sqrt(CorrYieldErr[i]*CorrYieldErr[i]+QuadError[i]*QuadError[i]+CorrelatedError[i]*CorrelatedError[i]);//Contains all of the errors and is not passed on to the crossXcomb
    CorrDibosonErr = CorrDibosonErr+CorrYieldErr[i]*CorrYieldErr[i];

    CorrFracVal[i] = CorrYieldVal[i]/ExpectedYield[i];
    CombFracErr[i] = CorrYieldErr[i]/ExpectedYield[i];
    YieldStatErr[i] = sqrt(DataEvts[i]);
    statErrIsgtPartialTotErr = false;
    if ( YieldStatErr[i]*YieldStatErr[i]>(CorrYieldErr[i]*CorrYieldErr[i]-CorrelatedError[i]*CorrelatedError[i]) ) {
      statErrIsgtPartialTotErr = true;
      // Temporary alter the values since only a combination of Stat & Uncorrelated Syst enters in the combination
      YieldStatErr[i] = sqrt(CorrYieldErr[i]*CorrYieldErr[i]-CorrelatedError[i]*CorrelatedError[i]);
      CorrYieldSystErr[i] = 0;
    } else {
      CorrYieldSystErr[i] = sqrt(CorrYieldErr[i]*CorrYieldErr[i]-CorrelatedError[i]*CorrelatedError[i]-YieldStatErr[i]*YieldStatErr[i]);
    }
    //cout << "i=" << i << ", CorrYieldErr[i]=" <<  CorrYieldErr[i] << ", YieldStatErr[i]=" <<  YieldStatErr[i] << ", CorrYieldSystErr[i]=" <<  CorrYieldSystErr[i] << endl;

    //Cross Section Computation Inputs
    CrossXVal[i] = CorrYieldVal[i]/(Ae[i]*Lumi[i]);
    sprintf(crossX_char,"%f",CrossXVal[i]);
    str_CrossXVal = str_CrossXVal + crossX_char;
    str_CrossXVal = str_CrossXVal + " ";

    CrossXErr[i] = CorrYieldErr[i]/(Ae[i]*Lumi[i]);

    CrossXStatErr[i] = YieldStatErr[i]/(Ae[i]*Lumi[i]);
    sprintf(crossX_char,"%f",CrossXStatErr[i]);
    str_CrossXStatErr = str_CrossXStatErr + crossX_char;
    str_CrossXStatErr = str_CrossXStatErr + " ";

    CrossXSystErr[i] = CorrYieldSystErr[i]/(Ae[i]*Lumi[i]);
    sprintf(crossX_char,"%f",CrossXSystErr[i]);
    str_CrossXSystErr = str_CrossXSystErr + crossX_char;
    str_CrossXSystErr = str_CrossXSystErr + " ";

    CrossXCorrelatedErr[i] = CorrelatedError[i]/(Ae[i]*Lumi[i]);
    sprintf(crossX_char,"%f",CrossXCorrelatedErr[i]);
    str_CrossXCorrelatedErr = str_CrossXCorrelatedErr + crossX_char;
    str_CrossXCorrelatedErr = str_CrossXCorrelatedErr + " ";
//     if ( statErrIsgtPartialTotErr ) {
//       //Set the YieldStatError back to sqrt{NData], etc.
//       YieldStatErr[i] = sqrt(DataEvts[i]);
//       CrossXStatErr[i] = YieldStatErr[i]/(Ae[i]*Lumi[i]);
//     }
  }
  FitDibosonErr=sqrt(FitDibosonErr);
  CorrDibosonErr=sqrt(CorrDibosonErr);

  cout << "Extracted Diboson Events : " << FitDibosonYield << "\\pm " << FitDibosonErr << endl;
  cout << "Corrected Diboson Events : " << CorrDibosonYield << "\\pm " << CorrDibosonErr << endl;

  int cC = 0;//Channel Counter
  bool produceSummaryTable = true;
  while ( cC<NCHANNELS ) {

    cout << ChannelName[cC] << "   ||   " << ChannelName[cC+1] << "   Results Summary" << endl;
    cout << "------------------------------- table:FitTotalsAndComparisons --------------------------" << endl;
    //    cout << setprecision(2) << scientific;
    cout << " Predicted Yield & Fitted Frac || Predicted Yield & Fitted Frac  " << endl;
    cout << "Diboson & " << setprecision(0) << fixed << ExpectedYield[cC] << " & " << setprecision(2) << fixed << FitFracVal[cC] << "$\\pm$" << FitFracErr[cC];
    cout  << "   &   "   << setprecision(0) << fixed << ExpectedYield[cC+1] << " & " << setprecision(2) << fixed << FitFracVal[cC+1] << "$\\pm$" << FitFracErr[cC+1] << endl;
    cout << "Corrected Diboson & " << setprecision(0) << fixed << ExpectedYield[cC] << " & " << setprecision(2) << fixed << CorrFracVal[cC] << "$\\pm$" << CombFracErr[cC];
    cout  << "   &   "   << setprecision(0) << fixed << ExpectedYield[cC+1] << " & " << setprecision(2) << fixed << CorrFracVal[cC+1] << "$\\pm$" << CombFracErr[cC+1] << endl;
    cout << "Data & " << setprecision(0) << fixed << DataEvts[cC] << " & ---        & " << DataEvts[cC+1] << " & --- " << endl;
    cout << endl;

    cout << "------------------------------- tab:dibosonYield --------------------------" << endl;
    cout << " Corrected Val || Corrected Val  " << endl;
    cout << "N Diboson & " << setprecision(0) << fixed << CorrYieldVal[cC] << "$\\pm$" << YieldStatErr[cC] << "$\\pm$" << sqrt(CorrYieldErr[cC]*CorrYieldErr[cC]-YieldStatErr[cC]*YieldStatErr[cC]);
    cout  << "   &   "   << CorrYieldVal[cC+1] << "$\\pm$" << YieldStatErr[cC+1] << "$\\pm$" << sqrt(CorrYieldErr[cC+1]*CorrYieldErr[cC+1]-YieldStatErr[cC+1]*YieldStatErr[cC+1]) << endl;
    cout << "$\\sigma^{fid}_{WW+WZ}$(pb) & " << setprecision(3) << fixed << Ae[cC]*CrossXVal[cC] << "$\\pm$" << Ae[cC]*CrossXStatErr[cC] << "$\\pm$" << Ae[cC]*sqrt(CrossXErr[cC]*CrossXErr[cC]-CrossXStatErr[cC]*CrossXStatErr[cC]);
    cout  << "   &   " << Ae[cC+1]*CrossXVal[cC+1] << "$\\pm$" << Ae[cC+1]*CrossXStatErr[cC+1] << "$\\pm$" << Ae[cC+1]*sqrt(CrossXErr[cC+1]*CrossXErr[cC+1]-CrossXStatErr[cC+1]*CrossXStatErr[cC+1]) << endl;
    cout << "$\\sigma_{WW+WZ}$(pb) & " << setprecision(2) << fixed << CrossXVal[cC] << "$\\pm$" << CrossXStatErr[cC] << "$\\pm$" << sqrt(CrossXErr[cC]*CrossXErr[cC]-CrossXStatErr[cC]*CrossXStatErr[cC]);
    cout  << "   &   " << CrossXVal[cC+1] << "$\\pm$" << CrossXStatErr[cC+1] << "$\\pm$" << sqrt(CrossXErr[cC+1]*CrossXErr[cC+1]-CrossXStatErr[cC+1]*CrossXStatErr[cC+1]) << endl;
    cout << endl;

    cout << "------------------------------- tab:SystematicCorrections --------------------------" << endl;
    cout << " Corrected Val || Corrected Val  " << endl;
    cout << "Fit Yield & " << setprecision(0) << fixed << CorrYieldVal[cC] << "$\\pm$" << CorrYieldErr[cC] << "   &   " << CorrYieldVal[cC+1] << "$\\pm$" << CorrYieldErr[cC+1] << endl;
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << endl;

    cC=cC+2;
  }

  if ( produceSummaryTable ) {
    cout << "------------------------------ Summary Table --------------------------" << endl;
    cout << "Channel " <<  ChannelName[0] << " | " <<  ChannelName[1] << " | " <<  ChannelName[2] << " | " <<  ChannelName[3] << " | " <<  ChannelName[4] << " | " <<  ChannelName[5] << endl;
    cout << "Data Events & " << setprecision(0) << fixed << DataEvts[0] << " & " << DataEvts[1] << " & " << DataEvts[2] << " & " << DataEvts[3] << " & " << DataEvts[4] << " & " << DataEvts[5] << "  \\\\" << endl;
    cout << "$N_{WW+WZ}$ Predicted & " << ExpectedYield[0] << " & " << ExpectedYield[1] << " & " << ExpectedYield[2] << " & " << ExpectedYield[3] << " & " << ExpectedYield[4] << " & " << ExpectedYield[5] << "  \\\\" << endl;

    cout << "$Frac_{WW+WZ}$ & " << setprecision(2) << FitYieldVal[0]/ExpectedYield[0] << "$\\pm$" << FitYieldErr[0]/ExpectedYield[0] << " & " << FitYieldVal[1]/ExpectedYield[1] << "$\\pm$" << FitYieldErr[1]/ExpectedYield[1] << " & " << FitYieldVal[2]/ExpectedYield[2] << "$\\pm$" << FitYieldErr[2]/ExpectedYield[2] << " & " << FitYieldVal[3]/ExpectedYield[3] << "$\\pm$" << FitYieldErr[3]/ExpectedYield[3] << " & " << FitYieldVal[4]/ExpectedYield[4] << "$\\pm$" << FitYieldErr[4]/ExpectedYield[4] << " & " << FitYieldVal[5]/ExpectedYield[5] << "$\\pm$" << FitYieldErr[5]/ExpectedYield[5] << "  \\\\" << endl;

    cout << "$Frac_{WW+WZ}$ Corrected & " << CorrYieldVal[0]/ExpectedYield[0] << "$\\pm$" << CorrYieldErr[0]/ExpectedYield[0] << " & " << CorrYieldVal[1]/ExpectedYield[1] << "$\\pm$" << CorrYieldErr[1]/ExpectedYield[1] << " & " << CorrYieldVal[2]/ExpectedYield[2] << "$\\pm$" << CorrYieldErr[2]/ExpectedYield[2] << " & " << CorrYieldVal[3]/ExpectedYield[3] << "$\\pm$" << CorrYieldErr[3]/ExpectedYield[3] << " & " << CorrYieldVal[4]/ExpectedYield[4] << "$\\pm$" << CorrYieldErr[4]/ExpectedYield[4] << " & " << CorrYieldVal[5]/ExpectedYield[5] << "$\\pm$" << CorrYieldErr[5]/ExpectedYield[5] << "  \\\\" << setprecision(0) << endl;

    cout << "$\\sigma^{fid}_{WW+WZ}$(\\unit{pb}) & " << setprecision(3) << fixed;
    for ( int i=0; i<NCHANNELS; i++) {
      cout << CrossXVal[i]*Ae[i] << "$\\pm$" << CrossXErr[i]*Ae[i]  << " & ";
    }
    cout  << "  \\\\" << endl;

    cout << "$\\mathcal{A}\\varepsilon$ & " << setprecision(2) << scientific << Ae[0] << " & " <<  Ae[1] << " & " <<  Ae[2] << " & " <<  Ae[3] << " & " <<  Ae[4] << " & " << Ae[5] << "  \\\\" << endl;

    cout << "$\\sigma^{Fit}_{WW+WZ}$(\\unit{pb}) & " << setprecision(1) << fixed << CrossXVal[0] << "$\\pm$" << CrossXErr[0]  << " & " <<  CrossXVal[1] << "$\\pm$" << CrossXErr[1]  << " & " <<  CrossXVal[2] << "$\\pm$" << CrossXErr[2]  << " & " <<  CrossXVal[3] << "$\\pm$" << CrossXErr[3]  << " & " <<  CrossXVal[4] << "$\\pm$" << CrossXErr[4]  << " & " <<  CrossXVal[5] << "$\\pm$" << CrossXErr[5]  << "  \\\\" << endl;
    cout << endl;
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << endl;
  }


  bool produceComparisonTable = true;
  if ( produceComparisonTable ) {
    cout << "------------------------------ Comparison Table --------------------------" << endl;
    cout << "Channel " <<  ChannelName[0] << " | " <<  ChannelName[1] << " | " <<  ChannelName[2] << " | " <<  ChannelName[3] << " | " <<  ChannelName[4] << " | " <<  ChannelName[5] << endl;

    cout << "$N_{WW+WZ}$ Fitted & " << setprecision(1) << FitYieldVal[0] << "$\\pm$" << FitYieldErr[0] << " & " << FitYieldVal[1] << "$\\pm$" << FitYieldErr[1] << " & " << FitYieldVal[2] << "$\\pm$" << FitYieldErr[2] << " & " << FitYieldVal[3] << "$\\pm$" << FitYieldErr[3] << " & " << FitYieldVal[4] << "$\\pm$" << FitYieldErr[4] << " & " << FitYieldVal[5] << "$\\pm$" << FitYieldErr[5] << "  \\\\" << endl;

    cout << "Yield Bias" << setprecision(1) << fixed;
    for ( int i=0; i<NCHANNELS; i++) {
      cout << " & " << YieldBias[i] ;
    }
    cout  << "  \\\\" << endl;

    cout << "Fractional Error Bias" << setprecision(2) << fixed;
    for ( int i=0; i<NCHANNELS; i++) {
      cout << " & " << ErrorInflation[i] ;
    }
    cout  << "  \\\\" << endl;

    cout << "V+Jets MC Systematic" << setprecision(0) << fixed; //Alternate W+Jets 
    for ( int i=0; i<NCHANNELS; i++) {
      cout << " & " << CorrelatedError[i] ;
    }
    cout  << "  \\\\" << endl;

    cout << "Shape Systematic" << setprecision(0) << fixed;
    for ( int i=0; i<NCHANNELS; i++) {
      cout << " & " << QuadError[i] ;
    }
    cout  << "  \\\\" << endl;

    cout << "$N_{WW+WZ}$ Corrected & " << setprecision(1) << CorrYieldVal[0] << "$\\pm$" << CorrYieldErr[0] << " & " << CorrYieldVal[1] << "$\\pm$" << CorrYieldErr[1] << " & " << CorrYieldVal[2] << "$\\pm$" << CorrYieldErr[2] << " & " << CorrYieldVal[3] << "$\\pm$" << CorrYieldErr[3] << " & " << CorrYieldVal[4] << "$\\pm$" << CorrYieldErr[4] << " & " << CorrYieldVal[5] << "$\\pm$" << CorrYieldErr[5] << "  \\\\" << setprecision(0) << endl;

    cout << "$\\mathcal{A}\\varepsilon$ & " << setprecision(5) << fixed << Ae[0] << " & " <<  Ae[1] << " & " <<  Ae[2] << " & " <<  Ae[3] << " & " <<  Ae[4] << " & " << Ae[5] << " \\\\" << endl;

    cout << "$\\sigma^{Fit}_{WW+WZ}$(\\unit{pb}) & " << setprecision(1) << fixed << CrossXVal[0] << "$\\pm$" << CrossXErr[0]  << " & " <<  CrossXVal[1] << "$\\pm$" << CrossXErr[1]  << " & " <<  CrossXVal[2] << "$\\pm$" << CrossXErr[2]  << " & " <<  CrossXVal[3] << "$\\pm$" << CrossXErr[3]  << " & " <<  CrossXVal[4] << "$\\pm$" << CrossXErr[4]  << " & " <<  CrossXVal[5] << "$\\pm$" << CrossXErr[5]  << "  \\\\" << endl;
    cout << endl;
    cout << "-------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << endl;
  }


  ///Compute the total cross section
  cout << "Computing the total cross section:" << endl;
  bool plotCrossX=true; // set to true in order to plot the cross sections
  crossXcomb(plotCrossX,str_CrossXVal,str_CrossXStatErr,str_CrossXSystErr,str_CrossXCorrelatedErr);
  ///Compute the total cross section


}

// -----------------------------------------------------------------------------------------------------------------------------------------------------------//


void computeTopControlSampleSystematic(bool isMuon=true){
//// Compute the systematic associated with the conversion from the top control region to the signal region by comparing the ratio to the one computed with alternate Matching/Scale samples

  TFile* f_default;
  TFile* f[4];
  TTree* InTree_default;
  TTree* InTree[4];
  TString cutstr_baseline = "(W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet>2.0)&&(GroomedJet_CA8_deltaR_lca8jet>1.57)&&(vbf_event==0)&&(GroomedJet_CA8_mass_pr[0]>40.000)&&(GroomedJet_CA8_mass_pr[0]<140.000)";//&&(numPFCorJetBTags<1)
  double NSigRegion, NContRegion, Rdefault, Ralt;
  


  if ( isMuon ) {
    cout << "Processing Systematics for the Muon channel" << endl;
    cutstr_baseline = "(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.)&&"+cutstr_baseline;
    f_default = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_mu_TTJetsPoheg_CMSSW532.root");
    f[0] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_mu_TTbar_scaleup_CMSSW532.root");
    f[1] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_mu_TTbar_scaledown_CMSSW532.root");
    f[2] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_mu_TTbar_matchup_CMSSW532.root");
    f[3] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_mu_TTbar_matchdown_CMSSW532.root");
  } else {
    cout << "Processing Systematics for the Electron channel" << endl;
    cutstr_baseline = "(event_met_pfmet >70)&&(W_electron_pt>35)&&"+cutstr_baseline;
    f_default = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_el_TTJetsPoheg_CMSSW532.root");
    f[0] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_el_TTbar_scaleup_CMSSW532.root");
    f[1] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_el_TTbar_scaledown_CMSSW532.root");
    f[2] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_el_TTbar_matchup_CMSSW532.root");
    f[3] = new TFile("/eos/uscms/store/user/lnujj/BoostedDibosonFitPostMoriond2013/RD_el_TTbar_matchdown_CMSSW532.root");
  }



  TString cutstr_signal = "(numPFCorJetBTags<1)&&" + cutstr_baseline;
  TString cutstr_control = "(numPFCorJetBTags>=1)&&" + cutstr_baseline;

  cout << "Signal Region Cuts = " << cutstr_signal << endl;
  cout << "Control Region Cuts = " << cutstr_control << endl;



  //Fill the histogram
  TString histload_str=">>htemp";
  histload_str="GroomedJet_CA8_mass_pr[0]"+histload_str;
  //  cout << "hist input =" << histload_str << endl;
  TH1F* htemp = new TH1F("htemp","htemp",20,40,140);
//   TH1F* hs_default, hc_default;
//   TH1F* hs[4], hc[4];

//   InTree->Draw(histload_str,Restrictions);
//   h = (TH1F*)htemp->Clone();
//   if ( dosumW2 ) {
//     h->Sumw2();
//   }



  InTree_default = (TTree*)f_default->Get("WJet");
  InTree_default->Draw(histload_str,cutstr_signal);
  NSigRegion=htemp->GetEntries();
  InTree_default->Draw(histload_str,cutstr_control);
  NContRegion=htemp->GetEntries();
  Rdefault=NSigRegion/NContRegion;
  cout << "Default Sample : NSigRegion=" << NSigRegion << ", NContRegion=" << NContRegion << ", Ratio=" << Rdefault << endl;

  for (Int_t n=0; n<4; n++) {
    InTree[n] = (TTree*)f[n]->Get("WJet");
    InTree[n]->Draw(histload_str,cutstr_signal);
    NSigRegion=htemp->GetEntries();
    InTree[n]->Draw(histload_str,cutstr_control);
    NContRegion=htemp->GetEntries();
    Ralt=NSigRegion/NContRegion;
    cout << "Alternate Sample " << n << endl;
    cout << "NSigRegion=" << NSigRegion << ", NContRegion=" << NContRegion << ", Ratio=" << Ralt << endl;
    cout << "Systematic Error: abs(Ralt/Rdefault - 1) = " << abs(Ralt/Rdefault -1) << endl;
  }




}



