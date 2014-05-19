#include <iostream>
#include "TLatex.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TCanvas.h"
//=====================================================================================
// SYNOPSIS:
//   1. Prepare "InData" and "OutDir" directories; e.g., "ln -s . OutDir" to go to current dir
//   2. root [0] .x mkControlPlots.C(0) for electron data, or
//      root [0] .x mkControlPlots.C(1) for muon data
//
// ====================================================================================
// Self Function
// ====================================================================================
/*

*/
void cmspre(double intlumifbinv)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.04);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.85,0.93,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.65,0.93,Form("#int #font[12]{L} dt = %.1f fb^{-1}", (float)intlumifbinv));
//  latex.DrawLatex(0.65,0.93,Form("2012B", (float)intlumifbinv));

  latex.SetTextAlign(11); // align left
//  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");
  latex.DrawLatex(0.15,0.96,"CMS preliminary");

}
#if 0
  PlotVar_t(char *inpv,double inmaxr,double inminr,int innbin,int inslog,char *inxl,char *inoutf,char *inoutf2,double inamax,double inamin,int inanb,int inhp,int indl) :
    plotvar(inpv),
    MAXRange(inmaxr),MINRange(inminr),NBINS(innbin),
    slog(inslog),xlabel(inxl),outfile(inoutf),outfile2(inoutf2),
    AMAXRange(inamax),AMINRange(inamin),ANBINS(inanb),
    hplot(inhp,drawleg(indl) {}
#endif
struct plotVar_t {
  char* plotvar;
  double MINRange;
  double MAXRange;
  int    NBINS;
  int    slog;
  char* xlabel;
  char* outfile;
  double AMINRange;
  double AMAXRange;
  int    ANBINS;
  int    mva_in;
  int    hplot;
  int  drawleg;
};

void mkControlPlotsDibosonBoosted(bool domu=true, const char* titlePrefix = "", const char* supplementalCutString="", bool scaleMCtoData=true,  bool domva=false,bool dovbf=false)
{
  gROOT->ProcessLine(".L tdrstyle.C");
  
  bool fitcorrYields=false;
  double fitcorr_WJets=1.0;
  double fitcorr_diboson=1.0;
  double fitcorr_top=1.0;
  if ( fitcorrYields ) {
    if ( domu ) {
      fitcorr_WJets=3299.9/3269.5;
      fitcorr_diboson=486.6/245.0;
      fitcorr_top=395.0/396.3;
    } else {
      fitcorr_WJets=2628.2/2768.8;
      fitcorr_diboson=509.8/218.4;
      fitcorr_top=322.9/322.1;
    }
  }
  double ttkfactor=0.95;
  if ( !domu ) { ttkfactor=0.92; }


  const double intLUMI = 19300.;
  const double H190_scale   = 8.782*0.115*1.5*2.* intLUMI/194504; // V
  const double H300_scale   = 3.606*0.101*1.5*2.* intLUMI/197284 ; // V
  const double H500_scale   = 1.439*0.080*1.5*2.* intLUMI/198470 ; // V

  //const double WJets_scale   = fitcorr_WJets*1.3*228.9* intLUMI/8955318; // pT>100 sample
  const double WJets_scale   = fitcorr_WJets*1.3*23.5* intLUMI/9492452; //pT>180 sample
  const double W4Jets_scale  = 172.6 * intLUMI/5000700;
  const double WW_scale      = fitcorr_diboson * 54.838   * intLUMI/9450414; // V
  const double WZ_scale      = fitcorr_diboson * 32.3161   * intLUMI/10000267; // V
  const double ZZ_scale      = 0.0*8.059   * intLUMI/9702850; // V

  const double QCD_scale     = 364000000 *    0.00037 * intLUMI/7529312 ; // V
  const double ZJets_scale   = 0.*3503.71  * intLUMI/30209426; // V
  const double ttbar_scale   = fitcorr_top * ttkfactor * 225.197 * intLUMI/20975917; // V
  const double SToppS_scale  = fitcorr_top * ttkfactor * 1.75776  * intLUMI/139974; // V anti-top
  const double SToppT_scale  = fitcorr_top * ttkfactor * 30.0042  * intLUMI/1935066; // V
  const double SToppTW_scale = fitcorr_top * ttkfactor * 11.1773  * intLUMI/493458; // V
  const double STopS_scale   = fitcorr_top * ttkfactor * 3.89394  * intLUMI/259960; // top
  const double STopT_scale   = fitcorr_top * ttkfactor * 55.531  * intLUMI/3758221; //V
  const double STopTW_scale  = fitcorr_top * ttkfactor * 11.1773  * intLUMI/497657; // V


  const plotVar_t plotvars[] = {

  //    plotvar	MINRange  MAXRange  NBINS  slog xlabel outfile AMINRange  AMAXRange ANBINS mva_in hplot drawleg
       { "GroomedJet_CA8_pt_pr[0]",   200,   800, 30, 3, "pruned jet p_{T}",        "GroomedJet_pt_pr",      200,  800, 30, 0, 0, 1 },
       { "GroomedJet_CA8_mass_pr[0]",   0,   140, 28, 3, "pruned jet mass",        "GroomedJet_mass",      0,  140, 28, 0, 0, 0 },
       // { "GroomedJet_CA8_mass_pr[0]",   0,   140, 28, 3, "pruned jet mass",        "GroomedJet_mass",      65,  105, 28, 0, 0, 1 }, // signal region
       { "(GroomedJet_CA8_pt[0]>150.0)+(GroomedJet_CA8_pt[1]>150.0)+(GroomedJet_CA8_pt[2]>150.0)+(GroomedJet_CA8_pt[3]>150.0)+(GroomedJet_CA8_pt[4]>150.0)+(GroomedJet_CA8_pt[5]>150.0)",     0, 5, 5, 3,   "Jet Multiplicity",            "CA8_JetMultiplicityg150",       0,  5, 5, 0, 0, 1 }, //jet multiplicity
       { "(W_muon_pfiso_sumChargedHadronPt+max(0,W_muon_pfiso_sumNeutralHadronEt+W_muon_pfiso_sumPhotonEt-0.5*W_muon_pfiso_sumPUPt))/W_muon_pt",   0,   0.12, 24, 3, "RelIso_PF",        "RelIso",      0,  0.12, 24, 0, 0, 0 },
       { "W_pt",   200,   800, 30, 3, "Leptonic W  p_{T}", "W_pt",      200,  800, 30, 0, 0, 1 },
       { "W_muon_pt",         0, 520, 26, 3,  "Muon p_{T} (GeV)",     "W_muon_pt",       0,  520, 26, 0, 0, 1 },
       { "W_muon_eta",      -2.7, 2.7, 18, 1,  "Muon #eta",            "W_muon_eta",    -2.4,  2.4, 16, 0, 0, 0 },
       { "W_electron_pt",         0, 520, 26, 3,  "Electron p_{T} (GeV)",     "W_elec_pt",       0,  520, 26, 0, 0, 1 },
       { "W_electron_eta",      -2.7, 2.7, 18, 1,  "Electron #eta",            "W_elec_eta",    -2.7,  2.7, 18, 0, 0, 0 },
       { "event_nPV",      0, 50., 50, 1,  "Num PV",            "event_nPV",    0.,  50., 50, 0, 0, 0 },
       { "event_met_pfmet",  0, 600, 24, 3,   "PF MET (GeV)",  "event_met_pfmet",     0, 600, 24, 0, 0, 0 },
       { "GroomedJet_numberbjets_csvm",     0, 5, 5, 3,   "nbjets_csmv",            "GroomedJet_numberbjets_csvm",       0,  5, 5, 0, 1, 1 },
       { "GroomedJet_CA8_pt[0]",   200,   800, 40, 3, "ungroomed jet p_{T}",        "GroomedJet_pt",      200,  800, 40, 0, 0, 0 },
       { "", 0.0,0.0,0,0,"","",0.,0.,0.,0,0,0 }
  };


  //  const char* the_cut = "1";
  //  double BINWIDTH = ((MAXRange-MINRange)/NBINS);

  //  Signal Region
  TString base_cutString = "((W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&( (GroomedJet_CA8_deltaR_lca8jet[1]<-900.0)||(GroomedJet_CA8_deltaR_lca8jet[1]>7.0) )&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet[0]>2.0)&&(GroomedJet_CA8_deltaR_lca8jet[0]>1.57)&&(numPFCorJetBTags<1)&&(vbf_event==0)&&(GroomedJet_CA8_mass_pr[0]>40.000)&&(GroomedJet_CA8_mass_pr[0]<140.000)&&";

//   // Top Control Sample Region
//   TString base_cutString = "effwt*puwt*((W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(GroomedJet_numberbjets_csvm>=1)&&(vbf_event==0)&&";

  TString lepton_cutString, data_cutString, mc_cutString;
  if ( domu ) {
    lepton_cutString="(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.)";
  } else {
    lepton_cutString="(event_met_pfmet >70)&&(W_electron_pt>35)";
  }
  TString cutString=")";
  if ( supplementalCutString!="" ) {
    cutString=supplementalCutString+cutString;
    cutString="&&"+cutString;
  }
  cutString=lepton_cutString+cutString;
  cutString= base_cutString+cutString;

  
  data_cutString = cutString;
  mc_cutString = "effwt*puwt*"+cutString;
  cout << "mc_cutString = " << mc_cutString << endl;

  //cutString="effwt*puwt*((W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&( (GroomedJet_CA8_deltaR_lca8jet[1]<-900.0)||(GroomedJet_CA8_deltaR_lca8jet[1]>7.0) )&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet[0]>2.0)&&(GroomedJet_CA8_deltaR_lca8jet[0]>1.57)&&(numPFCorJetBTags<1)&&(vbf_event==0)&&(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.)&&(GroomedJet_CA8_mass_pr[0]>40.000)&&(GroomedJet_CA8_mass_pr[0]<140.000))";


  //TCut the_cut("effwt*puwt*((event_nPV>15)&&(event_nPV<100)&&(W_pt>200.)&&(GroomedJet_CA8_pt[0]>200)&&(abs(GroomedJet_CA8_eta[0])<2.4)&&(GroomedJet_CA8_mass_pr[0]>40)&&(GroomedJet_CA8_tau2tau1[0]<0.55)&&(ggdboostedWevt==1)&&(GroomedJet_CA8_deltaphi_METca8jet>2.0)&&(GroomedJet_CA8_deltaR_lca8jet>1.57)&&(numPFCorJetBTags<1)&&(vbf_event==0)&&(event_met_pfmet >50)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>30.))");


  TCut thecut_mc(mc_cutString);
  TCut thecut_data(data_cutString);
//  TCut the_cut_data("effwt*puwt*((ggdevt==2||ggdevt==3) && fit_status==0  &&abs(JetPFCor_dphiMET[0])>0.4 && abs(W_muon_eta)<2.1&&event_runNo>193700)");

  if (dovbf)
    thecut_mc = TCut("vbf_event");

  TCut the_cut2("effwt*puwt*((ggdevt==2) && fit_status==0  &&abs(JetPFCor_dphiMET[0])>0.4 && abs(W_muon_eta)<2.1)");
  TCut the_cut3("effwt*puwt*((ggdevt==3) && fit_status==0  &&abs(JetPFCor_dphiMET[0])>0.4 && abs(W_muon_eta)<2.1)");
  //TCut the_cut("(fit_status==0 && TopWm!=0 )*effwt*puwt");

  TFile f("plotvar_histo.root", "RECREATE");

 // Get the input trees:

 // Data

  TFile *fin2,*H190_file,*H500_file,*H300_file,*wwShape_file,*wzShape_file,*zzShape_file,*wjetsShape_file,*w4jetShape_file,*ttbar_file,*qcd_file1,*zjets_file,*stops_file,*stopt_file,*stoptW_file,*stopbars_file,*stopbart_file,*stopbartW_file;


  if (domu) {
    fin2            = new TFile("InData/RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root", "read");
    wwShape_file    = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
    wzShape_file    = new TFile("InData/RD_mu_WZ_CMSSW532.root", "READ");
    zzShape_file    = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
    if (dovbf) {
      wjetsShape_file = new TFile("InData/RD_mu_W4Jets_CMSSW428.root","READ");
    } else {
      //wjetsShape_file = new TFile("InData/RD_mu_WJets_madgraph_CMSSW532.root", "READ");// pT>100 sample
      wjetsShape_file = new TFile("InData/RD_mu_WpJ_PT180_Madgraph_CMSSW532.root", "READ");// pT>180 sample
    }
    ttbar_file      = new TFile("InData/RD_mu_TTJets_poheg_CMSSW532_v2.root", "READ");
    //ttbar_file      = new TFile("InData/RD_mu_TTJetsPoheg_CMSSW532.root", "READ");
    zjets_file      = new TFile("InData/RD_mu_WW_CMSSW532.root", "READ");
    stops_file      = new TFile("InData/RD_mu_STopS_T_CMSSW532.root", "READ");
    stopt_file      = new TFile("InData/RD_mu_STopT_T_CMSSW532.root", "READ");
    stoptW_file     = new TFile("InData/RD_mu_STopTW_T_CMSSW532.root", "READ");
    stopbars_file =  new TFile("InData/RD_mu_STopS_Tbar_CMSSW532.root", "READ");
    stopbart_file =  new TFile("InData/RD_mu_STopT_Tbar_CMSSW532.root", "READ");
    stopbartW_file =  new TFile("InData/RD_mu_STopTW_Tbar_CMSSW532.root", "READ");


  } else { // electrons
    fin2            = new TFile("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root", "read");
    wwShape_file    = new TFile("InData/RD_el_WW_CMSSW532.root", "READ");
    wzShape_file    = new TFile("InData/RD_el_WZ_CMSSW532.root", "READ");
    zzShape_file    = new TFile("InData/RD_el_WW_CMSSW532.root", "READ");
    if (dovbf) {
      wjetsShape_file = new TFile("InData/RD_el_W4Jets_CMSSW428.root","READ"); 
    } else {
      //wjetsShape_file = new TFile("InData/RD_el_WJets_madgraph_CMSSW532.root", "READ");// pT>100 sample
      wjetsShape_file = new TFile("InData/RD_el_WpJ_PT180_Madgraph_CMSSW532.root", "READ");// pT>180 sample
    } 
    ttbar_file      = new TFile("InData/RD_el_TTJets_poheg_CMSSW532_v2.root", "READ");
    zjets_file      = new TFile("InData/RD_el_WW_CMSSW532.root", "READ");
    stops_file      = new TFile("InData/RD_el_STopS_T_CMSSW532.root", "READ");
    stopt_file      = new TFile("InData/RD_el_STopT_T_CMSSW532.root", "READ");
    stoptW_file     = new TFile("InData/RD_el_STopTW_T_CMSSW532.root", "READ");
    stopbars_file =  new TFile("InData/RD_el_STopS_Tbar_CMSSW532.root", "READ");
    stopbart_file =  new TFile("InData/RD_el_STopT_Tbar_CMSSW532.root", "READ");
    stopbartW_file =  new TFile("InData/RD_el_STopTW_Tbar_CMSSW532.root", "READ");

/*     fin2            = new TFile("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_3p5invfb.root", "READ"); */
/* //    H190_file       = new TFile("InData/RD_el_HWWMH190_CMSSW428.root", "READ"); */
/* //    H500_file       = new TFile("InData/RD_el_HWWMH500_CMSSW428.root", "READ"); */
/*     wwShape_file    = new TFile("InData/RD_el_WW_CMSSW532.root", "READ"); */
/*     wzShape_file    = new TFile("InData/RD_el_WZ_CMSSW532.root", "READ"); */
/*     if (dovbf) */
/*       wjetsShape_file = new TFile("InData/RD_el_W4Jets_CMSSW428.root","READ"); */
/*     else */
/*       wjetsShape_file = new TFile("InData/RD_el_WpJ_CMSSW532.root", "READ"); */
/*     ttbar_file      = new TFile("InData/RD_el_TTbar_CMSSW532.root", "READ"); */
/*     qcd_file1       = new TFile("InData/RDQCD_WenuJets_Isog0p3NoElMVA_1p6invfb.root", "READ"); */
/* //    qcd_file1       = new TFile("InData/RDQCD_WenuJets_DataAll_GoldenJSON_2p1invfb2011.root", "READ"); */
/*     zjets_file      = new TFile("InData/RD_el_ZpJ_CMSSW532.root", "READ"); */
/*     stops_file      = new TFile("InData/RD_el_STopS_T_CMSSW532.root", "READ"); */
/*     stopt_file      = new TFile("InData/RD_el_STopT_T_CMSSW532.root", "READ"); */
/*     stoptW_file     = new TFile("InData/RD_el_STopTW_T_CMSSW532.root", "READ"); */
  }

  TTree* treedata = (TTree*) fin2->Get("WJet");
  double nData = treedata->GetEntries();
  std::cout << "ndata =" << nData <<std::endl;

//  TTree* treeh190  = (TTree*)       H190_file->Get("WJet");
//  TTree* treeh500  = (TTree*)       H500_file->Get("WJet");
//  TTree* treeh300  = (TTree*)       H300_file->Get("WJet");
  TTree* treeww    = (TTree*)    wwShape_file->Get("WJet");
  TTree* treewz    = (TTree*)    wzShape_file->Get("WJet");
  TTree* treezz    = (TTree*)    zzShape_file->Get("WJet");
  TTree* treewj    = (TTree*) wjetsShape_file->Get("WJet");
  TTree* treettb   = (TTree*)      ttbar_file->Get("WJet");
//  TTree* treeqcd   = (TTree*)       qcd_file1->Get("WJet");
  TTree* treezj    = (TTree*)      zjets_file->Get("WJet");
  TTree* treests   = (TTree*)      stops_file->Get("WJet");
  TTree* treestt   = (TTree*)      stopt_file->Get("WJet");
  TTree* treestw   = (TTree*)     stoptW_file->Get("WJet");
  TTree* treestps   = (TTree*)      stopbars_file->Get("WJet");
  TTree* treestpt   = (TTree*)      stopbart_file->Get("WJet");
  TTree* treestpw   = (TTree*)     stopbartW_file->Get("WJet");


  for (int ivar=0; ; ivar++) {

    plotVar_t pv;
    if (dovbf)
      pv = vbfplotvars[ivar];
    else
      if(domva) 
        pv = higgsplotvars[ivar];
      else
        pv = plotvars[ivar];
 
    if ( !strlen(pv.plotvar) ) break;

    std::cout << TString(pv.plotvar) << "\t"<<pv.MINRange<<"\t" << pv.MAXRange<<"\t" << pv.NBINS<<"\tTHE CUT " << endl;

    if (domu) {
      if (strstr(pv.plotvar,"el")) continue;
    } else {
      if (strstr(pv.plotvar,"mu")) continue;
    }

    if (dovbf && pv.mva_in)
      thecut_mc = TCut("effwt*puwt*(vbf_event && vbf_wjj_m > 65.0 && vbf_wjj_m < 95.0)"); // plot only events in the signal region

    const double BINWIDTH = ((pv.MAXRange-pv.MINRange)/pv.NBINS);

    TH1* th1data  = new TH1D("th1data",  "th1data",  pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1data1 = new TH1D("th1data1", "th1data1", pv.NBINS, pv.MINRange, pv.MAXRange);
    TBox *errbox = new TBox(pv.AMINRange,0.974,pv.AMAXRange,1.026);
    errbox->SetFillColor(kYellow);
    treedata->Draw(TString(pv.plotvar)+TString(">>th1data"), thecut_data, "goff");
    treedata->Draw(TString(pv.plotvar)+TString(">>th1data1"), thecut_data, "goff");
 
    // Get WW/WZ/ZZ

    TH1* th1ww = new TH1D("th1ww", "th1ww", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1wz = new TH1D("th1wz", "th1wz", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1zz = new TH1D("th1zz", "th1zz", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1ww->Sumw2();
    th1wz->Sumw2();
    th1zz->Sumw2();

    treeww->Draw(TString(pv.plotvar)+TString(">>th1ww"), thecut_mc, "goff");
    treewz->Draw(TString(pv.plotvar)+TString(">>th1wz"), thecut_mc, "goff");
    treezz->Draw(TString(pv.plotvar)+TString(">>th1zz"), thecut_mc, "goff");

    // Get WJets

    TH1* th1wjets  = new TH1D("th1wjets",  "th1wjets",  pv.NBINS ,pv.MINRange,pv.MAXRange);

    treewj->Draw(TString(pv.plotvar)+TString(">>th1wjets"), thecut_mc, "goff");

    th1wjets->Sumw2();

    // Get ttbar
 
    TH1* th1Top = new TH1D("th1Top", "th1Top", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1Top->Sumw2();
    // cross section: 157.5 pb, events_gen = 3701947 (These are summer11 TTJets sample
  
    treettb->Draw(TString(pv.plotvar)+TString(">>th1Top"), thecut_mc, "goff");
/*
    // Get QCD

    TH1* th1qcd = new TH1D("th1qcd", "th1qcd", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1qcd->Sumw2();
 
    treeqcd->Draw(TString(pv.plotvar)+TString(">>th1qcd"), thecut_mc, "goff");
    int n2 = treeqcd->Draw(TString(pv.plotvar),  thecut_mc2, "goff");
    int n3 = treeqcd->Draw(TString(pv.plotvar),  thecut_mc3 , "goff");

    std::cout << "got qcd " << " n2 " << n2 <<  " n3  " << n3 <<std::endl; 

*/
    // Get Z+Jets

    TH1* th1zjets = new TH1D("th1zjets", "th1zjets", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1zjets->Sumw2();
    treezj->Draw(TString(pv.plotvar)+TString(">>th1zjets"), thecut_mc, "goff");

    // Get Single top

    TH1* th1stops = new TH1D("th1stops", "th1stops", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1stopt = new TH1D("th1stopt", "th1stopt", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1stoptw = new TH1D("th1stoptw", "th1stoptw", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1stops->Sumw2();
    th1stopt->Sumw2();
    th1stoptw->Sumw2();

    treests->Draw(TString(pv.plotvar)+TString(">>th1stops"), thecut_mc, "goff");
    treestt->Draw(TString(pv.plotvar)+TString(">>th1stopt"), thecut_mc, "goff");
    treestw->Draw(TString(pv.plotvar)+TString(">>th1stoptw"), thecut_mc, "goff");

    TH1* th1stopps = new TH1D("th1stopps", "th1stopps", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1stoppt = new TH1D("th1stoppt", "th1stoppt", pv.NBINS, pv.MINRange, pv.MAXRange);
    TH1* th1stopptw = new TH1D("th1stopptw", "th1stopptw", pv.NBINS, pv.MINRange, pv.MAXRange);
    th1stopps->Sumw2();
    th1stoppt->Sumw2();
    th1stopptw->Sumw2();
    treestps->Draw(TString(pv.plotvar)+TString(">>th1stopps"), thecut_mc, "goff");
    treestpt->Draw(TString(pv.plotvar)+TString(">>th1stoppt"), thecut_mc, "goff");
    treestpw->Draw(TString(pv.plotvar)+TString(">>th1stopptw"), thecut_mc, "goff");

    /* TFile* stopps_file =  new TFile("InData/RD_mu_STopS_Tbar_CMSSW532.root", "READ"); */
    /* TTree* tree64 = (TTree*) stopps_file->Get("WJet"); */
    /* TFile* stoppt_file =  new TFile("InData/RD_mu_STopT_Tbar_CMSSW532.root", "READ"); */
    /* TTree* tree65 = (TTree*) stoppt_file->Get("WJet"); */
    /* TFile* stopptW_file =  new TFile("InData/RD_mu_STopTW_Tbar_CMSSW532.root", "READ"); */
    /* TTree* tree66 = (TTree*) stopptW_file->Get("WJet"); */

    /* TH1* th1stopps = new TH1D("th1stopps", "th1stopps", pv.NBINS, pv.MINRange, pv.MAXRange); */
    /* TH1* th1stoppt = new TH1D("th1stoppt", "th1stoppt", pv.NBINS, pv.MINRange, pv.MAXRange); */
    /* TH1* th1stopptw = new TH1D("th1stopptw", "th1stopptw", pv.NBINS, pv.MINRange, pv.MAXRange); */
    /* th1stopps->Sumw2(); */
    /* th1stoppt->Sumw2(); */
    /* th1stopptw->Sumw2(); */
    /* tree64->Draw(TString(pv.plotvar)+TString(">>th1stopps"), thecut_mc, "goff"); */
    /* tree65->Draw(TString(pv.plotvar)+TString(">>th1stoppt"), thecut_mc, "goff"); */
    /* tree66->Draw(TString(pv.plotvar)+TString(">>th1stopptw"), thecut_mc, "goff"); */

    // Setup the canvas

//    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();
    tdrStyle->SetErrorX(0.5);
    tdrStyle->SetPadRightMargin(0.05);

    tdrStyle->SetLegendBorderSize(0);

    th1data->Sumw2();

    TCanvas* c1 = new TCanvas(pv.plotvar,pv.plotvar,10,10, 800, 800);
    TPad *d1, *d2;

    c1->Divide(1,2,0,0);
    d1 = (TPad*)c1->GetPad(1);
    d1->SetPad(0.01,0.30,0.95,0.99);
    d2 = (TPad*)c1->GetPad(2);
    d2->SetPad(0.01,0.02,0.95,0.30);

    // Setup the stack, scale the histos

    THStack* hs = new THStack("hs","MC contribution");
    th1Top->Scale(ttbar_scale);
    th1Top->SetFillColor(kGreen+2);
    th1Top->SetLineColor(kGreen+2);
    th1Top->SetLineWidth(0);

    th1stops->Scale(STopS_scale);
    th1stops->SetFillColor(7);
    th1stops->SetLineWidth(2);
    th1stopt->Scale(STopT_scale);
    th1stopt->SetFillColor(13);
    th1stopt->SetLineWidth(2);
    th1stoptw->Scale(STopTW_scale);
    th1stoptw->SetFillColor(9);
    th1stoptw->SetLineWidth(2);

    th1stopps->Scale(SToppS_scale);
    th1stopps->SetFillColor(7);
    th1stopps->SetLineWidth(2);
    th1stoppt->Scale(SToppT_scale);
    th1stoppt->SetFillColor(13);
    th1stoppt->SetLineWidth(2);
    th1stopptw->Scale(SToppTW_scale);
    th1stopptw->SetFillColor(9);
    th1stopptw->SetLineWidth(2);

    th1wjets->Scale(dovbf? W4Jets_scale : WJets_scale);
    th1wjets->SetFillColor(kRed);
    th1wjets->SetLineColor(kRed);
    th1wjets->SetLineWidth(0);
    th1ww->Scale(WW_scale);
    th1ww->SetFillColor(kAzure+8);
    th1ww->SetLineColor(kAzure+8);
    th1ww->SetLineWidth(0);
    th1wz->Scale(WZ_scale);
    th1wz->SetFillColor(11);
    th1wz->SetLineWidth(0);
    th1zz->Scale(ZZ_scale);
    th1zz->SetFillColor(11);
    th1zz->SetLineWidth(0);

    // th1qcd->Scale(QCD_scale);
  
//    th1qcd->SetFillColor(kGray+1);
//    th1qcd->SetLineColor(kGray+1);
//    th1qcd->SetLineWidth(0);
    th1zjets->Scale(ZJets_scale);
    th1zjets->SetFillColor(kYellow);
    th1zjets->SetLineColor(kYellow);
    th1zjets->SetLineWidth(0);
//    std::cout << " qcd " << th1qcd->Integral()   << std::endl;
    std::cout << "wjets "   << th1wjets->Integral()   << std::endl;
    std::cout << "tt "   << th1Top->Integral()   << std::endl;
    std::cout << "ww "   << th1ww->Integral()    << std::endl;
    std::cout << "wz "   << th1wz->Integral()    << std::endl;
    std::cout << "zz "   << th1zz->Integral()    << std::endl;
    std::cout << "zjets "    << th1zjets->Integral() << std::endl;
  
    double den_qcd = 
      th1Top->Integral()+
      th1stops->Integral()+
      th1stopt->Integral()+
      th1stoptw->Integral()+
      th1stopps->Integral()+
      th1stoppt->Integral()+
      th1stopptw->Integral()+
      th1wjets->Integral()+
      th1ww->Integral()+
      th1wz->Integral()+
      th1zz->Integral()+
      th1zjets->Integral();
/*
    double qcd_scale;

    if (domu)
      qcd_scale = (n2*0.002 + n3*0.000) / (n2+n3) ;//muon
    else
      qcd_scale = (n2*0.0637 + n3*0.02) / (n2+n3); //electron


*/
//    std::cout << " qcd_scale  " << qcd_scale <<std::endl;
//    th1qcd->Scale(qcd_scale*den_qcd/th1qcd->Integral()); 

    double den = 
      th1Top->Integral()+
      th1stops->Integral()+
      th1stopt->Integral()+
      th1stoptw->Integral()+
      th1stopps->Integral()+
      th1stoppt->Integral()+
      th1stopptw->Integral()+
      th1wjets->Integral()+
      th1ww->Integral()+
      th1wz->Integral()+
      th1zz->Integral()+
      th1zjets->Integral();
//      th1qcd->Integral();

    std::cout << "N MC Events after cuts = " << den <<std::endl;
    std::cout << "N Data Events after the cuts = " <<  th1data->Integral() << std::endl;

    if ( scaleMCtoData ) {
      cout << "Rescaling MC to Data" << endl;
      //    th1qcd->Scale   (th1data->Integral()/den); std::cout <<"qcd "   << th1qcd->Integral()   << std::endl;
      th1Top->Scale   (th1data->Integral()/den); std::cout <<"tt "    << th1Top->Integral()   << std::endl;
      th1stops->Scale (th1data->Integral()/den); std::cout <<"stops " << th1stops->Integral() << std::endl;
      th1stopt->Scale (th1data->Integral()/den); std::cout <<"stopt " << th1stopt->Integral() << std::endl;
      th1stoptw->Scale(th1data->Integral()/den); std::cout <<"stoptw "<< th1stoptw->Integral()<< std::endl;
      th1stopps->Scale (th1data->Integral()/den); std::cout <<"stops " << th1stopps->Integral() << std::endl;
      th1stoppt->Scale (th1data->Integral()/den); std::cout <<"stopt " << th1stoppt->Integral() << std::endl;
      th1stopptw->Scale(th1data->Integral()/den); std::cout <<"stoptw "<< th1stopptw->Integral()<< std::endl;
      th1wjets->Scale (th1data->Integral()/den); std::cout <<"wjets " << th1wjets->Integral() << std::endl;
      th1ww->Scale    (th1data->Integral()/den); std::cout <<"ww "    << th1ww->Integral()    << std::endl;
      th1wz->Scale    (th1data->Integral()/den); std::cout << "wz "   << th1wz->Integral()    << std::endl;
      th1zz->Scale    (th1data->Integral()/den); std::cout << "zz "   << th1zz->Integral()    << std::endl;
      th1zjets->Scale (th1data->Integral()/den); std::cout << "z "    << th1zjets->Integral() << std::endl;
      //    th1H500->Scale(th1data->Integral()/den);
      //    th1H300->Scale(th1data->Integral()/den);
      //    th1H190->Scale(th1data->Integral()/den);
    }
    cout<<"(th1data->Integral()/den) = "<< (th1data->Integral()/den) <<endl;
    double den2 =
      th1Top->Integral()+
      th1stops->Integral()+
      th1stopt->Integral()+
      th1stoptw->Integral()+
      th1stopps->Integral()+
      th1stoppt->Integral()+
      th1stopptw->Integral()+
      th1wjets->Integral()+
      th1ww->Integral()+
      th1wz->Integral()+
      th1zz->Integral()+
      th1zjets->Integral();
//      th1qcd->Integral();

    std::cout << "den2 " << den2 << std::endl;

    th1Top->Add(th1stopptw,1);
    th1Top->Add(th1stoppt,1);
    th1Top->Add(th1stopps,1);
    th1Top->Add(th1stoptw,1);
    th1Top->Add(th1stopt,1);
    th1Top->Add(th1stops,1);
    th1ww->Add(th1wz,1);
    th1ww->Add(th1zz,1);

    // Sum all the backgrounds

    TH1D *th1tot = (TH1D*)th1wjets->Clone("th1tot");
    th1tot->Reset();
    th1tot->Add(th1ww,1);
//    th1tot->Add(th1qcd,1);
    th1tot->Add(th1Top,1);
    th1tot->Add(th1wjets,1);
    th1tot->Add(th1zjets,1);
    TH1D* th1totClone = ( TH1D*) th1tot->Clone("th1totClone");
    th1totClone->SetMarkerStyle(0);
    th1totClone->SetFillStyle(3003);
    th1totClone->SetFillColor(11);
    th1totClone->SetLineColor(0);
    double binErr(0.0);
    for(int i=0; i<th1totClone->GetNbinsX(); ++i) {
      binErr = sqrt(
		    (th1ww->GetBinError(i))**2 +
//		    (th1qcd->GetBinError(i))**2 +
		    (th1Top->GetBinError(i))**2 +
		    (th1wjets->GetBinError(i))**2 +
		    (th1zjets->GetBinError(i))**2);
      th1totClone->SetBinError(i, binErr);
    }

    // Compose the stack

    hs->Add(th1zjets);
    d1->cd();
    gPad->SetBottomMargin(0.0);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);

//    hs->Add(th1qcd);
    hs->Add(th1Top);
    hs->Add(th1ww);
    hs->Add(th1wjets);

    // Set up the legend

    float  legX0=0.65, legX1=0.99, legY0=0.4, legY1=0.88;
    // float  legX0=0.35, legX1=0.85, legY0=0.4, legY1=0.88;
    // float  legX0=0.18, legX1=0.52, legY0=0.4, legY1=0.88;
    TLegend * Leg = new TLegend( legX0, legY0, legX1, legY1);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(0);
    Leg->SetTextSize(0.04);
    if (domu)
      Leg->AddEntry(th1data,  "Muon Data",  "PLE");
    else
      Leg->AddEntry(th1data,  "Electron Data",  "PLE");
    Leg->AddEntry(th1wjets,  "V+jets",  "f");
    Leg->AddEntry(th1ww,  "WV ",  "f");
    Leg->AddEntry(th1Top,  "top",  "f");
//    Leg->AddEntry(th1qcd,  "QCD",  "f");

//    Leg->AddEntry(th1zjets,  "Z+Jets",  "f");
    Leg->AddEntry(th1tot,  "MC Uncertainty",  "f");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")) Leg->AddEntry(th1H500,  "H (500) x 200",  "L");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")&&!strstr(pv.plotvar,"mva2j500mu")) Leg->AddEntry(th1H300,  "H (300) x 200",  "L");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j500mu")) Leg->AddEntry(th1H190,  "H (190) x 100",  "L");
    Leg->SetFillColor(0);

    ///Get the mean and RMS values
    cout << "WJets Mean = " << th1wjets->GetMean() << " +/- " << th1wjets->GetMeanError() << " || WJets RMS = " << th1wjets->GetRMS() << " +/- " << th1wjets->GetRMSError() << endl;
    TH1* th1comb  = new TH1D("th1comb",  "th1comb",  pv.NBINS, pv.MINRange, pv.MAXRange);
    th1comb->Add(th1Top);
    th1comb->Add(th1ww);
    th1comb->Add(th1wjets);
    cout << "All MC Mean = " << th1comb->GetMean() << " +/- " << th1comb->GetMeanError() << " || All MC RMS = " << th1comb->GetRMS() << " +/- " << th1comb->GetRMSError() << endl;
    cout << "Data Mean = " << th1data->GetMean() << " +/- " << th1data->GetMeanError() << " || Data RMS = " << th1data->GetRMS() << " +/- " << th1data->GetRMSError() << endl;
    delete th1comb;


/*
    th1H500->SetLineColor(kBlack);
    th1H500->SetLineWidth(3);
    th1H500->Scale(200);

    th1H300->SetLineColor(kBlack);
    th1H300->SetLineWidth(3);
    th1H300->SetLineStyle(2);
    th1H300->Scale(200);

    th1H190->SetLineColor(kBlue);
    th1H190->SetLineWidth(3);
    th1H190->Scale(100);
*/

    TH1* th1totempty = new TH1D("th1totempty", "th1totempty", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
    th1data->SetMarkerStyle(20);
    th1data->SetMarkerSize(1.25);
    th1data->SetLineWidth(2);

    th1tot->SetFillStyle(3001);
    th1tot->SetFillColor(1);
    th1tot->SetLineColor(1);
    th1tot->SetMarkerStyle(0);

    char tmpc[100];    sprintf(tmpc,"Events / %.1f GeV",BINWIDTH);
    if (pv.slog==1)    sprintf(tmpc,"Events/ %.1f",BINWIDTH);
    if (pv.slog==2)    sprintf(tmpc,"Events/ %.2f",BINWIDTH);
    if (pv.slog==3)    sprintf(tmpc,"Events/ %.0f GeV",BINWIDTH);
    if (pv.slog==6)    sprintf(tmpc,"Events/ %.1f rad",BINWIDTH);
    th1totempty->SetYTitle(tmpc);
    //  th1totempty->GetYaxis()->SetTitleSize(0.1);
    th1totempty->GetYaxis()->SetTitleOffset(1.2);
    th1totempty->GetYaxis()->SetLabelOffset(0.01);
    //  th1totempty->GetYaxis()->CenterTitle(true);
    th1totempty->GetYaxis()->SetLabelSize(0.04);
    // th1totClone->Draw("e3");   

    th1tot->SetMinimum(0.0);
    int maxbin = th1data->GetMaximumBin();
    float maxval = th1data->GetBinContent(maxbin);
    std::cout << "maxval " <<maxval <<std::endl;
//    th1totempty->SetMaximum(2.5*maxval);
    th1totempty->SetMaximum(1.6*maxval);
    th1totempty->SetMinimum(0.0);
    if(pv.slog==1) th1totempty->SetMaximum(1.6*maxval);
    th1data->SetMinimum(0.0);

    // Draw it all

    th1totempty->Draw();
    //th1tot->Draw("e2same");
    th1data->Draw("esame");
    hs->Draw("samehist");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")) th1H500->Draw("same");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j190mu")&&!strstr(pv.plotvar,"mva2j500mu")) th1H300->Draw("same");
//    if (pv.hplot ==1&&!strstr(pv.plotvar,"mva2j500mu")) th1H190->Draw("same");
    th1tot->Draw("e2same");

    th1data->Draw("esame");
    cmspre(intLUMI/1000.0);    
    if (pv.drawleg ==1)  Leg->Draw();  
    // th1data->Draw("Axissame");
    gPad->RedrawAxis();
    d2->cd();

    TH1F    * hhratio    = (TH1F*) th1data->Clone("hhratio")  ;
    hhratio->Sumw2();
    hhratio->SetStats(0);

    gPad->SetLeftMargin(0.14);
    gPad->SetTopMargin(0);
    gPad->SetRightMargin(0.05);
    gPad->SetFrameBorderSize(0);
    gPad->SetBottomMargin(0.3);
    gPad->SetTickx();

    hhratio->SetMarkerSize(1.25);
    //  hhratio->GetYaxis()->SetRangeUser(0.48,1.52);
    hhratio->GetYaxis()->SetRangeUser(0.3,1.7);
    hhratio->GetXaxis()->SetTitle(pv.xlabel);
    hhratio->GetXaxis()->SetTitleOffset(0.9);
    hhratio->GetXaxis()->SetTitleSize(0.15);
    hhratio->GetXaxis()->SetLabelSize(0.15);
    hhratio->SetYTitle("Ratio Data/MC");
    hhratio->GetYaxis()->SetTitleSize(0.1);
    hhratio->GetYaxis()->SetTitleOffset(0.5);
    hhratio->GetYaxis()->CenterTitle(true);
    hhratio->GetYaxis()->SetLabelSize(0.1);
    std::cout << hhratio->GetNbinsX() << std::endl;
    std::cout << th1tot->GetNbinsX() << std::endl;
    hhratio->Divide(th1tot);
    double binError(0.0), mcbinentry(0.0), mcerror(0.0);
    for(int i=0; i<hhratio->GetNbinsX(); ++i) {
      binError = hhratio->GetBinError(i);
      mcerror = th1tot->GetBinError(i);
      mcbinentry = th1tot->GetBinContent(i);
      if(mcbinentry>0.) mcerror /= mcbinentry;
      else mcerror = 0.0;
      binError = sqrt(binError**2 + mcerror**2);
      hhratio->SetBinError(i, binError);
    }
    TH1D *th1emptyclone = new TH1D("th1emptyclone", "th1emptyclone", pv.ANBINS, pv.AMINRange, pv.AMAXRange);
    th1emptyclone->GetYaxis()->SetRangeUser(0.6,1.3999);
    th1emptyclone->GetXaxis()->SetTitle(pv.xlabel);
    th1emptyclone->GetXaxis()->SetTitleOffset(0.9);
    th1emptyclone->GetXaxis()->SetTitleSize(0.15);
    th1emptyclone->GetXaxis()->SetLabelSize(0.15);
    th1emptyclone->SetYTitle("Ratio Data/MC");
    th1emptyclone->GetYaxis()->SetTitleSize(0.1);
    th1emptyclone->GetXaxis()->SetNdivisions(505);
    th1emptyclone->GetYaxis()->SetNdivisions(505);
    th1emptyclone->GetYaxis()->SetTitleOffset(0.5);
    th1emptyclone->GetYaxis()->CenterTitle(true);
    th1emptyclone->GetYaxis()->SetLabelSize(0.1);
    th1emptyclone->Draw();
    errbox->Draw();
    hhratio->Draw("esame");
    TLine *line; line = new TLine(pv.AMINRange,1.0,pv.AMAXRange,1.0);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    line->SetLineColor(1);
    line->Draw();

    TString outfile = TString("OutDir/")+(domu?TString("mu_"):TString("el_"))+TString(titlePrefix)+TString(pv.outfile);

    //c1->Print(outfile+"_topControlSample.png");
    c1->Print(outfile+".png");
    //c1->Print(outfile+".C");
    //gPad->WaitPrimitive();
    c1->Modified();
    c1->Update();
    //c1->SaveAs(outfile+".pdf"); 

  } // var loop

  // f.Write();

}

