#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TCut.h"
#include "TString.h"
#include "TF1.h"
#include "TPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TProfile2D.h"
#include "THStack.h"


#undef ATGCISSIGNAL

const double intLUMIinvpb_mu = 19300.; // for normalization of signal
const double intLUMIinvpb_el = 19200.; // (background norm comes from data)

//double WJets_scale   = 37509.0    * intLUMIinvpb/(18353019+50768992);
//double WJets_scale   = 228.9*1.3  * intLUMIinvpb/8955318;
//double WJets_scale   =  26.4*1.3  * intLUMIinvpb/9492452;
//double ZJets_scale   =  3503.71   * intLUMIinvpb/30209426;

double WJ_scale_lo   =    23.5    / 9492452; // pt180 samp xsec from PREP
double WJ_scale_hi   =    0.0587  / 155501;  // pt600 samp

//double WW_scale_lo   =    33.61     / 9450414; // xsec from PREP
//double WW_scale_hi   =    0.005235  / 1000139; // xsec from PREP

// aMC@NLO samples produced privately:
double WW_scale_lo   =    0.2222038 * 2 * 1.636 * 3      / 119692;
double WW_scale_hi   =    0.2222038 * 2 * 1.079e-02 * 3  / 118889;

double WW_scale_NLO  =    57.2  / 33.615;   // xsec from MCFM paper (arxiv:1105.0020v1)
double WZ_scale_NLO  =    22.88 / 10000267; // xsec from MCFM paper - on-shell Z, no Wgam* component

double stitchvalgev = 640;
double stitchvalgev = 640;

// For WV, take NLO prediction as opposed to our fit results (smaller uncertainty)
//double WWplusWZscale =    80.1 * intLUMIinvpb;

double ttbar_scale   =   225.197 /20975917;
double mtt700_scale  =    15.38  / 3082808;
double mtt1k_scale   =     3.063 / 1249109;
double SToppS_scale  =     1.76  /  139974;
double SToppT_scale  =    30.7   / 1935066;
double SToppTW_scale =    11.1   /  493458;
double STopS_scale   =     3.79  /  259960;
double STopT_scale   =    56.4   / 3758221;
double STopTW_scale  =    11.1   /  497657;

double  mc2data_scale;

// NORMALIZE TO THE CROSS-SECTION FIT YIELD RESULTS IN THE SIGNAL REGION:
//

// //Muons (boosted):
// double    WV_fityield_muboosted =  317.822509443; //  +/-  132.893054439
// double   top_fityield_muboosted =  659.894037959; //  +/-   52.869716158
// double WJets_fityield_muboosted = 1827.11376406;  //  +/-   66.7367004534

// //Electrons (boosted):
// double    WV_fityield_elboosted =  387.089737128; // +/-  106.493654322
// double   top_fityield_elboosted =  539.602494805; // +/-   53.1969479612
// double WJets_fityield_elboosted = 1380.22945171;  // +/-   52.9665002837

// YIELDS WITH SECOND JET VETO (Pt<150):

// //Muons (boosted):
// double    WV_fityield_muboosted =    281.733229243; //  ±  120.611950138
// double   top_fityield_muboosted =    593.98965759;  //  ±   47.5785980599 
// double WJets_fityield_muboosted =   1701.45698625;  //  ±   61.5501720859 
				    								    
// //Electrons (boosted):		    							    
// double    WV_fityield_elboosted =    298.816791102; //  ±  125.156825962
// double   top_fityield_elboosted =    488.108124474; //  ±   47.9910346173
// double WJets_fityield_elboosted =   1309.77948103;  //  ±   62.8213214703 

// YIELDS WITH SECOND JET/lepton sep. reqrmnt. (deltaR(l,ca8jet)>7.0):

//Muons (boosted):
double WV_fityield_muboosted    =  369.97472832;  //  +/-  128.653534669
double top_fityield_muboosted   =  249.405360996; //  +/-  20.0052840402
double WJets_fityield_muboosted = 1349.50888792;  //  +/-  72.6138461002

//Electrons (boosted):
double WV_fityield_elboosted    =  387.529582723; //  +/-  93.0619917601
double top_fityield_elboosted   =  202.06545164;  //  +/-  20.1262682472
double WJets_fityield_elboosted = 1075.33103014;  //  +/-  53.2322420384

// Yield uncertainties for sum(W+jets,top)
// computed with correlations from the fit,
// and divided by the total background

double mu_bkgd_norm_error = (77.65)/(WJets_fityield_muboosted+top_fityield_muboosted);
double el_bkgd_norm_error = (70.76)/(WJets_fityield_elboosted+top_fityield_elboosted);

// Add absolute errors for WW, W+Z, and W-Z NLO from MCFM paper in quadrature
// (arxiv:1105.0020v1)
//
double mu_sig_norm_error =  sqrt((.041 * .041 * 57.25 *57.25)+
				 (.052 * .052 * 14.48 *14.48)+
				 (.054 * .054 * 8.4   *  8.4)
				 )/(57.25+14.48+8.4); // (132.90)/317.82;

double el_sig_norm_error =  mu_sig_norm_error; // (106.49)/387.09;

bool domu=true;
float yRatioMin = 0.2;
float yRatioMax = 2.2;

//TFile* sigLambdaZ;

////////// ALL input trees ///////////
TTree* treedata;
TTree* treewwlo;
TTree* treewwhi;
TTree* treewz;
TTree* treewjlo;
TTree* treewjhi;
TTree* treettblo,*treettbmd,*treettbhi;
//TTree* treeqcd;
TTree* treezj;
TTree* treests;
TTree* treestt;
TTree* treestw;
TTree* tree64;
TTree* tree65;
TTree* tree66;



////////// ALL histograms ///////////
TH1* th1data;
TH1* th1wwlo;
TH1* th1wwhi;
TH1* th1ww;
TH1* th1wv;
TH1* th1wz;
TH1* th1wjets;
TH1* th1wjlo;
TH1* th1wjhi;
TH1* systUp;
TH1* systDown;
TH1* th1toplo;
TH1* th1topmd;
TH1* th1tophi;
TH1* th1Top;
//TH1* th1qcd;
TH1* th1zjets;
TH1* th1stops;
TH1* th1stopt;
TH1* th1stoptw;
TH1* th1stopps;
TH1* th1stoppt;
TH1* th1stopptw;
TH1D *th1tot;
TH1D* th1totClone;
TH1* th1totempty;
TH1D* th1emptyclone;
TH1F* hhratio;
TH1F* hhratioUp;
TH1F* hhratioDown;
TH1* signalForDisplay;
TH1* signalRatioForDisplay;

bool saveDataCards_ = true;
//bool saveDataCards_ = false;
TF1 *gaus2;

//======================================================================


void InstantiateTrees() {

  ////////// ALL input files ///////////
  TFile* fin2;
  TFile* wwShapeLo_file;
  TFile* wwShapeHi_file;
  TFile* wzShape_file;
  TFile* wjShapeLo_file;
  TFile* wjShapeHi_file;
  TFile* ttbar_file;
  TFile* ttb_mtt700to1k;
  TFile* ttb_mtt1ktoinf;
  //  TFile* qcd_file1;
  TFile* zjets_file;
  TFile* stops_file;
  TFile* stopt_file;
  TFile* stoptW_file;
  TFile* stopps_file;
  TFile* stoppt_file;
  TFile* stopptW_file;

 if (domu) {
   fin2            = new TFile("InData/RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root", "read");
   wwShapeLo_file  = new TFile("InMC/RD_mu_WW_minPt150_amcnlo_CMSSW532.root", "READ");
   wwShapeHi_file  = new TFile("InMC/RD_mu_WW_minPt500_amcnlo_CMSSW532.root", "READ");
   wzShape_file    = new TFile("InMC/RD_mu_WZ_CMSSW532.root", "READ");
// wjShape_file    = new TFile("InMC/RD_mu_WpJ_CMSSW532.root", "READ");
   wjShapeLo_file  = new TFile("InMC/RD_mu_WpJ_PT180_Madgraph_CMSSW532.root", "READ");
   wjShapeHi_file  = new TFile("InMC/RD_mu_WpJ_minPt600_CMSSW532.root", "READ");
// ttbar_file      = new TFile("InMC/RD_mu_TTbar_CMSSW532.root", "READ");
   ttbar_file      = new TFile("InMC/RD_mu_TTJets_poheg_CMSSW532.root", "READ");
   ttb_mtt700to1k  = new TFile("InMC/RD_mu_Mtt700to1000_CMSSW532.root", "READ");
   ttb_mtt1ktoinf  = new TFile("InMC/RD_mu_Mtt1000toinf_CMSSW532.root", "READ");
// qcd_file1       = new TFile("InMC/RD_mu_WpJ_PT180_Madgraph_CMSSW532.root", "READ"); // use WpJ, shapes are the same
// zjets_file      = new TFile("InMC/RD_mu_ZpJ_CMSSW532.root", "READ");
   stops_file      = new TFile("InMC/RD_mu_STopS_T_CMSSW532.root", "READ");
   stopt_file      = new TFile("InMC/RD_mu_STopT_T_CMSSW532.root", "READ");
   stoptW_file     = new TFile("InMC/RD_mu_STopTW_T_CMSSW532.root", "READ");
   stopps_file     = new TFile("InMC/RD_mu_STopS_Tbar_CMSSW532.root", "READ");
   stoppt_file     = new TFile("InMC/RD_mu_STopT_Tbar_CMSSW532.root", "READ");
   stopptW_file    = new TFile("InMC/RD_mu_STopTW_Tbar_CMSSW532.root", "READ");
 } else { // electrons

   fin2            = new TFile("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_19p2invfb.root", "READ");
   wwShapeLo_file  = new TFile("InMC/RD_el_WW_minPt150_amcnlo_CMSSW532.root", "READ");
   wwShapeHi_file  = new TFile("InMC/RD_el_WW_minPt500_amcnlo_CMSSW532.root", "READ");
   wzShape_file    = new TFile("InMC/RD_el_WZ_CMSSW532.root", "READ");
// wjShape_file    = new TFile("InMC/RD_el_WpJ_CMSSW532.root", "READ");
   wjShapeLo_file  = new TFile("InMC/RD_el_WpJ_PT180_Madgraph_CMSSW532.root", "READ");
   wjShapeHi_file  = new TFile("InMC/RD_el_WpJ_minPt600_CMSSW532.root", "READ");
// ttbar_file      = new TFile("InMC/RD_el_TTbar_CMSSW532.root", "READ");
   ttbar_file      = new TFile("InMC/RD_el_TTJets_poheg_CMSSW532.root", "READ"); // boosted sample
   ttb_mtt700to1k  = new TFile("InMC/RD_el_Mtt700to1000_CMSSW532.root", "READ"); // boosted sample
   ttb_mtt1ktoinf  = new TFile("InMC/RD_el_Mtt1000toinf_CMSSW532.root", "READ"); // boosted sample
// qcd_file1       = new TFile("InData/RDQCD_WenuJets_Isog0p3NoElMVA_19p2invfb.root", "READ");
// zjets_file      = new TFile("InMC/RD_el_ZpJ_CMSSW532.root", "READ");
   stops_file      = new TFile("InMC/RD_el_STopS_T_CMSSW532.root", "READ");
   stopt_file      = new TFile("InMC/RD_el_STopT_T_CMSSW532.root", "READ");
   stoptW_file     = new TFile("InMC/RD_el_STopTW_T_CMSSW532.root", "READ");
   stopps_file     = new TFile("InMC/RD_el_STopS_Tbar_CMSSW532.root", "READ");
   stoppt_file     = new TFile("InMC/RD_el_STopT_Tbar_CMSSW532.root", "READ");
   stopptW_file    = new TFile("InMC/RD_el_STopTW_Tbar_CMSSW532.root", "READ");
 }


  treedata = (TTree*) fin2->Get("WJet");
  double nData = treedata->GetEntries();
  std::cout << "ndata =" << nData <<std::endl;

  //// ------------ Get all trees
  treewwlo  = (TTree*)    wwShapeLo_file->Get("WJet");
  treewwhi  = (TTree*)    wwShapeHi_file->Get("WJet");
  treewz    = (TTree*)    wzShape_file->Get("WJet");
  treewjlo  = (TTree*)    wjShapeLo_file->Get("WJet");
  treewjhi  = (TTree*)    wjShapeHi_file->Get("WJet");
  treettblo = (TTree*)    ttbar_file->Get("WJet");
  treettbmd = (TTree*)    ttb_mtt700to1k->Get("WJet");
  treettbhi = (TTree*)    ttb_mtt1ktoinf->Get("WJet");
  //treeqcd   = (TTree*)    qcd_file1->Get("WJet");
  //treezj    = (TTree*)    zjets_file->Get("WJet");
  treests   = (TTree*)    stops_file->Get("WJet");
  treestt   = (TTree*)    stopt_file->Get("WJet");
  treestw   = (TTree*)    stoptW_file->Get("WJet");
  tree64    = (TTree*)    stopps_file->Get("WJet");
  tree65    = (TTree*)    stoppt_file->Get("WJet");
  tree66    = (TTree*)    stopptW_file->Get("WJet");

  //// ------------ Create a tree branch for dijet pt
  const char* dijetPt = "sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))";
  treedata->SetAlias("dijetPt", dijetPt);
  treewwlo->SetAlias("dijetPt", dijetPt);
  treewwhi->SetAlias("dijetPt", dijetPt);
  treewz->SetAlias("dijetPt", dijetPt);
  treewjlo->SetAlias("dijetPt", dijetPt);
  treewjhi->SetAlias("dijetPt", dijetPt);
  treettblo->SetAlias("dijetPt", dijetPt);
  //treeqcd->SetAlias("dijetPt", dijetPt);
  //treezj->SetAlias("dijetPt", dijetPt);
  treests->SetAlias("dijetPt", dijetPt);
  treestt->SetAlias("dijetPt", dijetPt);
  treestw->SetAlias("dijetPt", dijetPt);
  tree64->SetAlias("dijetPt", dijetPt);
  tree65->SetAlias("dijetPt", dijetPt);
  tree66->SetAlias("dijetPt", dijetPt);
}                                                    // InstantiateTrees

//======================================================================

void ScaleHistos(int channel)
{

  bool domu = true;
  if(channel==1 || channel==3) domu = false;

#if 0
  double WJets_norm = 1.;
  double VV_norm    = 1.;
  double Top_norm   = 1.;
  WJets_scale   *= WJets_norm;
  WW_scale      *= VV_norm;
  WZ_scale      *= VV_norm;
  ttbar_scale   *= Top_norm;
  SToppS_scale  *= Top_norm;
  SToppT_scale  *= Top_norm;
  SToppTW_scale *= Top_norm;
  STopS_scale   *= Top_norm;
  STopT_scale   *= Top_norm;
  STopTW_scale  *= Top_norm;
#endif

  // Print all bin contents prior to scaling:
  printf ("%-5s %7s %10s %10s %9s %11s %9s\n","bin","Pt(GeV)","TTbar","WJets","WW","WZ","Data");
  for (int ibin=1; ibin<=th1wz->GetNbinsX()+1; ibin++) {
    if (ibin>th1wz->GetNbinsX())
      printf ("%-4d ovrflow", ibin);
    else
      printf ("%-4d %3.0f-%3.0f", ibin,
	      th1wz->GetXaxis()->GetBinLowEdge(ibin),
	      th1wz->GetXaxis()->GetBinUpEdge(ibin));
    printf (" %5.4g+-%4.1f",  th1toplo->GetBinContent(ibin)+
	                      th1topmd->GetBinContent(ibin)+
	                      th1tophi->GetBinContent(ibin),
	                      th1toplo->GetBinError(ibin)+
	                      th1topmd->GetBinError(ibin)+
	                      th1toplo->GetBinError(ibin));
    if (th1wz->GetBinLowEdge(ibin) < stitchvalgev)
      printf (" %5.4g+-%4.1f",th1wjlo->GetBinContent(ibin),th1wjlo->GetBinError(ibin));
    else
      printf (" %5.4g+-%4.1f",th1wjhi->GetBinContent(ibin),th1wjhi->GetBinError(ibin));

    if (th1wz->GetBinLowEdge(ibin) < stitchvalgev)
      printf (" %7.1f+-%4.1f",   th1wwlo->GetBinContent(ibin),   th1wwlo->GetBinError(ibin));
    else
      printf (" %7.1f+-%4.1f",   th1wwhi->GetBinContent(ibin),   th1wwhi->GetBinError(ibin));

    printf (" %5.3g+-%4.1f",   th1wz->GetBinContent(ibin),   th1wz->GetBinError(ibin));
    printf (" %5.3g+-%4.1f", th1data->GetBinContent(ibin), th1data->GetBinError(ibin));
    printf ("\n");
  }

  double WV_fityield;
  double WJets_fityield;
  double top_fityield;
  switch(channel) {
  case 2:                                    //muon boosted
    WJets_fityield = WJets_fityield_muboosted;
    WV_fityield    =    WV_fityield_muboosted;
    top_fityield   =   top_fityield_muboosted;
    break;
  case 3:                                    //electron boosted
    WJets_fityield = WJets_fityield_elboosted;
    WV_fityield    =    WV_fityield_elboosted;
    top_fityield   =   top_fityield_elboosted;
    break;
  default:
    cerr << "only boosted channels, please" << endl;
    exit(-1);
  }

  th1toplo->Scale(ttbar_scale);
  th1topmd->Scale(mtt700_scale);
  th1tophi->Scale(mtt1k_scale);

  th1stops->Scale(STopS_scale);
  th1stops->SetFillColor(7);
  th1stops->SetLineWidth(2);
  th1stopt->Scale(STopT_scale);
  th1stopt->SetFillColor(13);
  th1stopt->SetLineWidth(2);
  th1stoptw->Scale(STopTW_scale);
  th1stoptw->SetFillColor(9);
  th1stoptw->SetLineWidth(2);

  // Add all the top together:
  th1topmd->Add(th1tophi,1);
  th1Top->Add(th1toplo,th1topmd,1.0,0.5);
  th1Top->Add(th1stoptw,1);
  th1Top->Add(th1stopt,1);
  th1Top->Add(th1stops,1);
  th1Top->SetFillColor(kGreen+2);
  th1Top->SetLineColor(kGreen+2);
  th1Top->SetLineWidth(0);

  th1wjlo->Scale(WJ_scale_lo);
  th1wjhi->Scale(WJ_scale_hi);
#if 1
  // stitch the two histograms together into one, post-smearing.
  for (int ibin=0;ibin<=th1wjets->GetNbinsX()+1;ibin++) {
    if (th1wjets->GetBinLowEdge(ibin)<stitchvalgev) {
      th1wjets->SetBinContent(ibin,th1wjlo->GetBinContent(ibin));
      th1wjets->SetBinError  (ibin,th1wjlo->GetBinError  (ibin));
    } else {
      th1wjets->SetBinContent(ibin,th1wjhi->GetBinContent(ibin));
      th1wjets->SetBinError  (ibin,th1wjhi->GetBinError  (ibin));
    }
  }
#else
  th1wjets = (TH1 *)th1wjlo->Clone();
#endif

  th1wjets->SetFillColor(kRed);
  th1wjets->SetLineColor(kRed);
  th1wjets->SetLineWidth(0);

  th1wwlo->Scale(WW_scale_lo);
  th1wwhi->Scale(WW_scale_hi);
#if 1
  // stitch the two histograms together into one, post-smearing.
  for (int ibin=0;ibin<=th1ww->GetNbinsX()+1;ibin++) {
    if (th1ww->GetBinLowEdge(ibin)<560) {
      th1ww->SetBinContent(ibin,th1wwlo->GetBinContent(ibin));
      th1ww->SetBinError  (ibin,th1wwlo->GetBinError  (ibin));
    } else {
      th1ww->SetBinContent(ibin,th1wwhi->GetBinContent(ibin));
      th1ww->SetBinError  (ibin,th1wwhi->GetBinError  (ibin));
    }
  }
#else
  th1ww = (TH1 *)th1wwlo->Clone();
#endif


  th1ww->Scale(WW_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));
  th1ww->SetFillColor(kAzure+8);
  th1ww->SetLineColor(kAzure+8);
  th1ww->SetLineWidth(0);

  th1wz->Scale(WZ_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));
  th1wz->SetFillColor(11);
  th1wz->SetLineWidth(0);

  //th1qcd->SetFillColor(kGray+1);
  //th1qcd->SetLineColor(kGray+1);
  //th1qcd->SetLineWidth(0);
  //th1zjets->Scale(ZJets_scale);
  //th1zjets->SetFillColor(kYellow);
  //th1zjets->SetLineColor(kYellow);
  //th1zjets->SetLineWidth(0);
  //std::cout << " qcd " << th1qcd->Integral()   << std::endl;
  std::cout << "wjets "  << th1wjets->Integral() << std::endl;
  std::cout << "tt "     << th1Top->Integral()   << std::endl;
  std::cout << "ww "     << th1ww->Integral()    << std::endl;
  std::cout << "wz "     << th1wz->Integral()    << std::endl;
  //std::cout << "z "    << th1zjets->Integral() << std::endl;

#if 0
  double den_qcd = 
    th1Top->Integral()+
    //th1stops->Integral()+
    //th1stopt->Integral()+
    //th1stoptw->Integral()+
    th1wjets->Integral()+
    th1ww->Integral()+
    th1wz->Integral();
    //th1zjets->Integral();

  double qcd_scale = 0.0;

  if (channel == 1) // electron, resolved dijet
    qcd_scale = 0.03;

  //std::cout << " qcd_scale  " << qcd_scale <<std::endl;
  th1qcd->Scale(qcd_scale*den_qcd/th1qcd->Integral()); 
#endif

  std::cout <<" data " <<  th1data->Integral() << std::endl;

#if 0  // no longer scaling total to data...

  // include SM WW in the backgrounds to determine MC/data normalization
  double den = 
    th1Top->Integral()+
    th1stops->Integral()+
    th1stopt->Integral()+
    th1stoptw->Integral()+
    th1wjets->Integral()+
    th1ww->Integral()+
    th1wz->Integral();
    //th1zjets->Integral()+
    //th1qcd->Integral();

  std::cout << "den = " <<den <<std::endl;
  std::cout <<" data " <<  th1data->Integral() << std::endl;

  mc2data_scale = th1data->Integral()/den;

  //th1qcd->Scale   (mc2data_scale); std::cout <<"qcd "   << th1qcd->Integral()   << std::endl;
  th1Top->Scale   (mc2data_scale); std::cout <<"tt "    << th1Top->Integral()   << std::endl;
  th1stops->Scale (mc2data_scale); std::cout <<"stops " << th1stops->Integral() << std::endl;
  th1stopt->Scale (mc2data_scale); std::cout <<"stopt " << th1stopt->Integral() << std::endl;
  th1stoptw->Scale(mc2data_scale); std::cout <<"stoptw "<< th1stoptw->Integral()<< std::endl;
  th1wjets->Scale (mc2data_scale); std::cout <<"wjets " << th1wjets->Integral() << std::endl;
  th1ww->Scale    (mc2data_scale); std::cout <<"ww "    << th1ww->Integral()    << std::endl;
  th1wz->Scale    (mc2data_scale); std::cout <<"wz "    << th1wz->Integral()    << std::endl;
  //th1zjets->Scale (mc2data_scale); std::cout <<"z "     << th1zjets->Integral() << std::endl;

  double den2 =
    th1Top->Integral()+
    th1stops->Integral()+
    th1stopt->Integral()+
    th1stoptw->Integral()+
    th1wjets->Integral()+
    th1ww->Integral()+
    th1wz->Integral()+
    //th1zjets->Integral()+
    //th1qcd->Integral();

  std::cout << "den2 " << den2 << std::endl;
#endif

  th1wv = (TH1 *)th1ww->Clone("wv");
  th1wv->Add(th1wz,1);

#if 0
  double  wvscale = WV_fityield/ th1wv->Integral();
  th1ww ->Scale( wvscale );
  th1wz ->Scale( wvscale );
  th1wv ->Scale( wvscale );
#endif

  // rescale combined histos to the fit yield
  //
  th1wjets->Scale(WJets_fityield/th1wjets->Integral());

  th1Top->Scale(top_fityield/th1Top->Integral());

  // Print totals after scaling
  printf ("scaled:\n");

  // Print all bin contents prior after scaling:
  printf ("%-5s %7s %10s %10s %9s %9s\n","bin","Pt(GeV)","Top","WJets","WW+WZ","Data");
  for (int ibin=1; ibin<=th1wv->GetNbinsX()+1; ibin++) {
    if (ibin>th1wv->GetNbinsX())
      printf ("%-4d ovrflow", ibin);
    else
      printf ("%-4d %3.0f-%3.0f", ibin,
	      th1wv->GetXaxis()->GetBinLowEdge(ibin),
	      th1wv->GetXaxis()->GetBinUpEdge(ibin));
    printf (" %5.4g+-%4.2f",  th1Top->GetBinContent(ibin),  th1Top->GetBinError(ibin));
    printf (" %5.3g+-%4.1f",th1wjets->GetBinContent(ibin),th1wjets->GetBinError(ibin));
    printf (" %5.2f+-%4.2f",   th1wv->GetBinContent(ibin),   th1wv->GetBinError(ibin));
    printf (" %5.3g+-%4.1f", th1data->GetBinContent(ibin), th1data->GetBinError(ibin));
    printf ("\n");
  }
  
  //mc2data_scale = th1data->Integral()/(top_fityield+WV_fityield+WJets_fityield);

}                                                         // ScaleHistos


//======================================================================

void AddOverflowBin(TH1* hist) {
  int nBinsTot = hist->GetNbinsX();
  double lastbin = hist->GetBinContent(nBinsTot);
  double overflow = hist->GetBinContent(nBinsTot+1);
  hist->SetBinContent(nBinsTot, lastbin+overflow);
  hist->SetBinError(nBinsTot, sqrt(lastbin+overflow));
}

//======================================================================

void AddOverflowBin(TH1* hist, TF1 *f, double dm_max) {
  int nBinsTot = hist->GetNbinsX();
  double lastbin = hist->GetBinContent(nBinsTot);
  double overflow = f->Integral(dm_max, dm_max + 400.);
  hist->SetBinContent(nBinsTot, lastbin+overflow);
  hist->SetBinError(nBinsTot, sqrt(lastbin+overflow));
}

//======================================================================

// Sum all the backgrounds
void SumAllBackgrounds() {

  //-------- First add overflow bin ----------------
  AddOverflowBin(th1ww);
  AddOverflowBin(th1wz);
  AddOverflowBin(th1wv);

  AddOverflowBin(th1wjets);
  //AddOverflowBin(th1zjets);
  AddOverflowBin(th1Top);
  //AddOverflowBin(th1qcd);
  AddOverflowBin(th1data);

  //-------- Now sum of all bkg histograms ----------
  th1tot = (TH1D*)th1wjets->Clone("th1tot");
  th1tot->Reset();
#ifdef ATGCISSIGNAL
  cout << "I'M ADDING WV TO THE BACKGROUND!!!!" << endl;
  th1tot->Add(th1wv,1);
#endif
  //th1tot->Add(th1qcd,1);
  th1tot->Add(th1Top,1);
  th1tot->Add(th1wjets,1);
  //th1tot->Add(th1zjets,1);

  th1tot->SetFillStyle(3001);
  th1tot->SetFillColor(1);
  th1tot->SetLineColor(1);
  th1tot->SetMarkerStyle(0);
  th1tot->SetMinimum(0.0);


  //-------- Needed for plotting ----------
  th1totClone = ( TH1D*) th1tot->Clone("th1totClone");

#ifndef ATGCISSIGNAL // wasn't added before, so add now for plotting
  th1totClone->Add(th1wv,1);
#endif
  th1totClone->SetMarkerStyle(0);
  th1totClone->SetFillStyle(3003);
  th1totClone->SetFillColor(11);
  th1totClone->SetLineColor(0);

  double binErr(0.0);
  for(int i=0; i<th1totClone->GetNbinsX(); ++i) {
    binErr = sqrt(
		  (th1wv->GetBinError(i))**2 +
		  //(th1qcd->GetBinError(i))**2 +
		  (th1Top->GetBinError(i))**2 +
		  (th1wjets->GetBinError(i))**2);
		  //(th1zjets->GetBinError(i))**2);
    th1totClone->SetBinError(i, binErr);
  }
  

  //-------- Ratio histogram ----------
  hhratio    = (TH1F*) th1data->Clone("hhratio")  ;
  hhratio->Sumw2();
  hhratio->SetMarkerSize(1.25);
  hhratio->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
  hhratio->Divide(th1totClone);
  double binError(0.0), mcbinentry(0.0), mcerror(0.0);
  for(int i=0; i<hhratio->GetNbinsX(); ++i) {
    binError = hhratio->GetBinError(i);
    mcerror = th1tot->GetBinError(i);
    mcbinentry = th1tot->GetBinContent(i);
    if(mcbinentry>0.) mcerror /= mcbinentry;
    else mcerror = 0.0;
    binError = sqrt(binError**2 + mcerror**2);
    hhratio->SetBinError(i, binError);
    if(th1data->GetBinContent(i)<0.1) hhratio->SetBinError(i,0.0);
  }
}                                                   // SumAllBackgrounds

//======================================================================


TLegend* GetLegend(int channel)
{
  // float  legX0=0.5, legX1=0.89, legY0=0.41, legY1=0.86;
  float  legX0=0.6, legX1=0.93, legY0=0.45, legY1=0.88;
  if(channel > 1) { legX0=0.6; legY0=0.48; }

  bool domu = true;
  if(channel==1 || channel==3) domu = false;

  TLegend * Leg = new TLegend( legX0, legY0, legX1, legY1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(0);
  Leg->SetTextSize(0.05);
  Leg->AddEntry(th1data, domu ?  "Muon Data" : "Electron Data",  "PLE");
  Leg->AddEntry(th1wv,  "SM WW+WZ ",  "f");
  Leg->AddEntry(th1wjets,  "W+jets",  "f");
  Leg->AddEntry(th1Top,  "top",  "f");
  //if(channel == 1) Leg->AddEntry(th1qcd,  "Multijet",  "f");
  //Leg->AddEntry(th1zjets,  "Z+Jets",  "f");
  Leg->AddEntry(th1tot,  "MC error",  "f");
  //Leg->AddEntry(systUp,  "Shape error",  "f");
  //Leg->AddEntry(signalForDisplay,  "#lambda_{Z}=0.03",  "l");
  Leg->SetFillColor(0);

  return Leg;
}


//======================================================================

void SetupEmptyHistogram(int bins, double dm_min, double dm_max, char* xtitle)
//void SetupEmptyHistogram(int bins, double* ptbins, char* xtitle)
{
  //th1totempty = new TH1D("th1totempty", "th1totempty", bins, ptbins); // dm_min, dm_max);
  th1totempty = new TH1D("th1totempty", "th1totempty", bins,dm_min, dm_max);
  char tmpc[100];    sprintf(tmpc,"Events / %d GeV", (int) (dm_max-dm_min)/bins);
  th1totempty->SetYTitle(tmpc);
  th1totempty->GetYaxis()->SetTitleOffset(1);
  th1totempty->GetYaxis()->SetLabelOffset(0.01);
  th1totempty->GetYaxis()->SetLabelSize(0.08);
  th1totempty->GetYaxis()->SetTitleSize(0.08);

  int maxbin = th1data->GetMaximumBin();
  float maxval = th1data->GetBinContent(maxbin);
  std::cout << "maxval " <<maxval <<std::endl;
  th1totempty->SetMaximum(1.2*maxval);
  th1totempty->SetMinimum(0.0);


  //th1emptyclone = new TH1D("th1emptyclone", "th1emptyclone", bins, ptbins); // dm_min, dm_max);
  th1emptyclone = new TH1D("th1emptyclone", "th1emptyclone", bins, dm_min, dm_max);
  th1emptyclone->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
  th1emptyclone->GetXaxis()->SetTitle(xtitle);
  th1emptyclone->GetXaxis()->SetTitleOffset(0.9);
  th1emptyclone->GetXaxis()->SetTitleSize(0.15);
  th1emptyclone->GetXaxis()->SetLabelSize(0.15);
  th1emptyclone->SetYTitle("Data/MC  ");
  th1emptyclone->GetYaxis()->SetTitleSize(0.2);
  th1emptyclone->GetXaxis()->SetNdivisions(505);
  th1emptyclone->GetYaxis()->SetNdivisions(503);
  th1emptyclone->GetYaxis()->SetTitleOffset(0.4);
  th1emptyclone->GetYaxis()->CenterTitle(true);
  th1emptyclone->GetYaxis()->SetLabelSize(0.2);
  th1emptyclone->GetXaxis()->SetLabelSize(0.2);
  th1emptyclone->GetXaxis()->SetTitleSize(0.2);
  th1emptyclone->GetYaxis()->SetMoreLogLabels();
  th1emptyclone->GetYaxis()->SetNoExponent();
}                                                //  SetupEmptyHistogram

//======================================================================

//------- Get signal histogram -------

TCut
GetSigRatioFunction(float lambdaZ, float dkappaGamma, float deltaG1,const char *obsrvbl, const char *wworwz="ww")
{
  //printf("%g\t%g\t%g\r",lambdaZ,dkappaGamma,deltaG1);

  TFile f("ATGC_shape_coefficients.root");
  TProfile2D* p0_prof; 
  TProfile2D* p1_prof; 
  TProfile2D* p2_prof; 
  TProfile2D* p3_prof; 
  TProfile2D* p4_prof; 
  float var1, var2;

  if(fabs(deltaG1)<0.000001) {
    p0_prof = (TProfile2D*) f.Get(Form("%s_p0_lambda_dkg",wworwz));
    p1_prof = (TProfile2D*) f.Get(Form("%s_p1_lambda_dkg",wworwz));
    p2_prof = (TProfile2D*) f.Get(Form("%s_p2_lambda_dkg",wworwz));
    p3_prof = (TProfile2D*) f.Get(Form("%s_p3_lambda_dkg",wworwz));
    p4_prof = (TProfile2D*) f.Get(Form("%s_p4_lambda_dkg",wworwz));
    var1 = lambdaZ; 
    var2 = dkappaGamma; 
    //printf("Looking up coefficients for %s lambdaZ=%g, dkappaGamma=%g\n",wworwz,var1,var2);
  }
  else if(fabs(dkappaGamma)<0.000001) { 
    p0_prof = (TProfile2D*) f.Get(Form("%s_p0_lambda_dg1",wworwz));
    p1_prof = (TProfile2D*) f.Get(Form("%s_p1_lambda_dg1",wworwz));
    p2_prof = (TProfile2D*) f.Get(Form("%s_p2_lambda_dg1",wworwz));
    p3_prof = (TProfile2D*) f.Get(Form("%s_p3_lambda_dg1",wworwz));
    p4_prof = (TProfile2D*) f.Get(Form("%s_p4_lambda_dg1",wworwz));
    var1 = lambdaZ; 
    var2 = deltaG1;
    //printf("Looking up coefficients for %s lambdaZ=%g, deltaG1=%g\n",wworwz,var1,var2);
  }
  else if(fabs(lambdaZ)<0.000001) {
    p0_prof = (TProfile2D*) f.Get(Form("%s_p0_dkg_dg1",wworwz));
    p1_prof = (TProfile2D*) f.Get(Form("%s_p1_dkg_dg1",wworwz));
    p2_prof = (TProfile2D*) f.Get(Form("%s_p2_dkg_dg1",wworwz));
    p3_prof = (TProfile2D*) f.Get(Form("%s_p3_dkg_dg1",wworwz));
    p4_prof = (TProfile2D*) f.Get(Form("%s_p4_dkg_dg1",wworwz));
    var1 = dkappaGamma; 
    var2 = deltaG1;
    //printf("Looking up coefficients for %s dkappaGamma=%g, deltaG1=%g\n",wworwz,var1,var2);
  }
  else {
    assert(0);
  }

  TCut sigratio(
#ifdef ATGCISSIGNAL
    TString::Format("(%g-1.0+",p0_prof->Interpolate(var1,var2))+
#else
    TString::Format("(%g+",    p0_prof->Interpolate(var1,var2))+
#endif
    TString::Format("%s*(%g+",obsrvbl,p1_prof->Interpolate(var1,var2))+
    TString::Format("%s*(%g+",obsrvbl,p2_prof->Interpolate(var1,var2))+
    TString::Format("%s*(%g+",obsrvbl,p3_prof->Interpolate(var1,var2))+
    TString::Format("%s*%g))))",obsrvbl,p4_prof->Interpolate(var1,var2))
    );

  return sigratio;
}                                                 // GetSigRatioFunction


//======================================================================

void cmspre()
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.05);

  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.85,0.93,"#sqrt{s} = 8 TeV");
  latex.SetTextAlign(31); // align right
  latex.DrawLatex(0.65,0.93,Form("#scale[0.5]{#lower[-0.15]{#it{#int}}}#it{L} dt = %0.1f#kern[0.2]{fb}^{-1}", 19.3));
  latex.SetTextAlign(11); // align left
//  latex.DrawLatex(0.15,0.93,"CMS,  #sqrt{s} = 7 TeV");//preliminary 2011");
  latex.DrawLatex(0.15,0.93,"CMS Preliminary");

}

//======================================================================

void fillMChisto(TH1 *thehisto,
		 TTree *thetree,
		 const TString& theobservable,
		 const TCut& thecut)
{
  TString drawstr = theobservable+Form(">>%s",thehisto->GetName());

  cout << "filling histo " << thehisto->GetName() << endl;

  cout <<"tree->Draw("<<drawstr<<","<<TString((const char *)thecut)<<",goff);"<<endl;
  thetree->Draw(drawstr,TCut(thecut),"goff");
}                                                         // fillMChisto

//======================================================================

const double ptbins_boosted[15] = {
  200,225,250,275,300,325,350,375,400,425,450,500,550,625,925
};

#define NUMETABINS 5
const double etabins[NUMETABINS+1] = {
  0.0,0.5,1.1,1.7,2.3,5.0
};

//======================================================================

//////---------- channel: 0==muon dijet, 1== electron dijet
/////                     2==muon  boosted,   3== electron boosted
void makeATGCLimitDataCards(int channel) {

//   const Int_t bins = 8; 
//   const Float_t dm_min = 200.; 
//   const Float_t dm_max = 600.;

  Int_t bins = 7; 
  Float_t dm_min = 100.; 
  Float_t dm_max = 275.;
  if(channel>1) { bins = 15; dm_min = 200.; dm_max = 800; }

  Int_t binsmc = (int)(10*(dm_max - dm_min)); // units of 0.1GeV, set up for smearing post-binning!

  domu = true;
  if(channel==1 || channel==3) domu = false;

 
  TString outfile = (domu?TString("mu_"):TString("el_"))+ 
    (channel<2?TString("dijet"):TString("boosted"));
  TFile* outputForLimit = TFile::Open(outfile+".root", "recreate");

  h_smearfactor = new TH1D("h_smearfactor","h_smearfactor",100,0.5,1.5);
  h_sigma_MC    = new TH1D("h_sigma_MC","h_sigma_MC",50,0,1);


  TString cutsDijet("(W_pt<200.) && (dijetPt>70.) && (abs(JetPFCor_Eta[0])<2.4) && (abs(JetPFCor_Eta[1])<2.4) && (abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5) &&(abs(JetPFCor_dphiMET[0])>0.4) &&(W_mt>30.) &&(JetPFCor_Pt[0]>40.) &&(JetPFCor_Pt[1]>35.) &&(JetPFCor_Pt[2]<30.) &&(JetPFCor_bDiscriminatorCSV[0]<0.244) &&(JetPFCor_bDiscriminatorCSV[1]<0.244) && (Mass2j_PFCor>70. && Mass2j_PFCor<100.)");



  // Do not put jet pt in the cut string here, since it is going to be smeared
  TString cutsMerged("(vbf_event==0) && (W_pt>200.) &&(abs(GroomedJet_CA8_eta[0])<2.4)&&(ggdboostedWevt==1) && (GroomedJet_CA8_deltaphi_METca8jet[0]>2.0) && (GroomedJet_CA8_deltaR_lca8jet[0]>1.57) && (GroomedJet_CA8_deltaR_lca8jet[1]<-900 || GroomedJet_CA8_deltaR_lca8jet[1]>7.0) && (numPFCorJetBTags<1) && (GroomedJet_CA8_tau2tau1[0]<0.55) && (GroomedJet_CA8_mass_pr[0]>70. && GroomedJet_CA8_mass_pr[0]<100.)");



  TString        lepton_cut = "(event_met_pfmet >30) && (W_electron_pt>35.)";
  if(channel==0) lepton_cut = "(event_met_pfmet >25) &&(abs(W_muon_eta)<2.1) && (W_muon_pt>25.)";
  if(channel==1) lepton_cut = "(event_met_pfmet >30) && (W_electron_pt>30.)";
  if(channel==2) lepton_cut = "(event_met_pfmet >50) &&(abs(W_muon_eta)<2.1) && (W_muon_pt>30.)";
  if(channel==3) lepton_cut = "(event_met_pfmet >70) && (W_electron_pt>35.)";

  TString and(" && ");

  TString jet_cut = cutsDijet;
  TString jetptcut,mcjetptcut;

  if(channel>1) jet_cut = cutsMerged;

  char* observable = "dijetPt";
  char* mcobservable = "dijetPt";
  char* xtitle = "p_{T}^{jj} [GeV]"; 
  double jetptcutval;
  if(channel>1) {
    double jetthresh = 80;

    /*observable = 
      "(GroomedJet_CA8_pt[0]>jetthresh)+\
      (GroomedJet_CA8_pt[1]>jetthresh)+\
      (GroomedJet_CA8_pt[2]>jetthresh)+\
      (GroomedJet_CA8_pt[3]>jetthresh)+\
      (GroomedJet_CA8_pt[4]>jetthresh)+\
      (GroomedJet_CA8_pt[5]>jetthresh)";*/

    mcobservable = "GroomedJet_CA8_pt_smeared[0]";
    //mcobservable = "GroomedJet_CA8_pt[0]";
    observable = "GroomedJet_CA8_pt[0]";
    jetptcutval = 200.;
    jetptcut   = Form("(%s > %f)",  observable,jetptcutval);
    mcjetptcut = Form("(%s > %f)",mcobservable,jetptcutval);
    xtitle = "p_{T}^{j} [GeV]";
  }


  TCut mccut( TString("(effwt*puwt)*(")+ lepton_cut+and+jet_cut+and+mcjetptcut + TString(")") );
  TCut datacut( TString("(") + lepton_cut+and+jet_cut+and+jetptcut + TString(")") );
  //TCut the_cut( TString("(effwt*puwt)*(")+ lepton_cut+and+jet_cut+and+jetptcut + TString(")") );

  // for combining ttbar files
  TString mttstr
    ("sqrt((W_top_E+W_atop_E)^2 -(W_top_px+W_atop_px)^2 -(W_top_py+W_atop_py)^2 -(W_top_pz+W_atop_pz)^2)");

  // weight default sample by 1 for mtt<700, by half above for adding the high mtt samples
  TCut mttwt(TString("(")+mttstr+TString("<700)?1.0:0.5"));

  InstantiateTrees();


  //th1data  = new TH1D("th1data",  "th1data",  bins, ptbins_boosted); // bins, dm_min, dm_max);
  th1data  = new TH1D("th1data",  "th1data",  bins, dm_min, dm_max);
  th1data->Sumw2();
  th1data->SetMarkerStyle(20);
  th1data->SetMarkerSize(1.25);
  th1data->SetLineWidth(2);
  th1data->SetMinimum(0.0);

  TString drawstr = TString(observable)+TString(">>th1data");

  cout <<
    TString("treedata->Draw(\"")+drawstr+TString("\", \"")+
    TString((const char*)datacut)+TString("\", \"goff\")") << endl;

  treedata->Draw(drawstr, datacut, "goff");

  // ------- Get WW/WZ ------- 
  th1wwlo = new TH1D("th1wwlo", "th1wwlo", bins, dm_min, dm_max);
  th1wwhi = new TH1D("th1wwhi", "th1wwhi", bins, dm_min, dm_max);
  th1ww   = new TH1D("th1ww", "th1ww", bins, dm_min, dm_max);
  th1wz   = new TH1D("th1wz", "th1wz", bins, dm_min, dm_max);
  //th1wz = new TH1D("th1wz", "th1wz", bins, ptbins_boosted);
  th1wwlo->Sumw2();
  th1wwhi->Sumw2();
  th1wz->Sumw2();

  fillMChisto(th1wwlo,treewwlo,mcobservable,mccut);
  fillMChisto(th1wwhi,treewwhi,mcobservable,mccut);

  fillMChisto(th1wz,treewz,mcobservable,mccut);

  // ------- Get ttbar ------- 
  th1Top   = new TH1D("th1Top", "th1Top", bins, dm_min, dm_max);
  th1toplo = new TH1D("th1toplo", "th1toplo", bins, dm_min, dm_max);
  th1topmd = new TH1D("th1topmd", "th1topmd", bins, dm_min, dm_max);
  th1tophi = new TH1D("th1tophi", "th1tophi", bins, dm_min, dm_max);

  //th1Top = new TH1D("th1Top", "th1Top", bins, ptbins_boosted);
  th1Top->Sumw2();
  th1toplo->Sumw2();
  th1topmd->Sumw2();
  th1tophi->Sumw2();

  fillMChisto(th1toplo,treettblo,mcobservable,mttwt*mccut);
  fillMChisto(th1topmd,treettbmd,mcobservable,mccut);
  fillMChisto(th1tophi,treettbhi,mcobservable,mccut);

    // ------- Get WJets ------- 
  th1wjlo  = new TH1D("th1wjlo",  "th1wjlo",  bins, dm_min, dm_max);
  th1wjhi  = new TH1D("th1wjhi",  "th1wjhi",  bins, dm_min, dm_max);
  th1wjets = new TH1D("th1wjets", "th1wjets", bins, dm_min, dm_max);
  th1wjlo->Sumw2();
  th1wjhi->Sumw2();

  fillMChisto(th1wjlo,treewjlo,mcobservable,mccut);
  fillMChisto(th1wjhi,treewjhi,mcobservable,mccut);

  //th1wjets  = new TH1D("th1wjets",  "th1wjets", bins, ptbins_boosted);

  // ------- Get QCD ------- 
  //th1qcd = new TH1D("th1qcd", "th1qcd", bins, dm_min, dm_max);
  //th1qcd = new TH1D("th1qcd", "th1qcd", bins, ptbins_boosted);
  //th1qcd->Sumw2();
  //treeqcd->Draw(TString(observable)+TString(">>th1qcd"), mccut, "goff");

  // ------- Get Z+Jets ------- 
  //th1zjets = new TH1D("th1zjets", "th1zjets", bins, dm_min, dm_max);
  //th1zjets->Sumw2();
  //treezj->Draw(TString(observable)+TString(">>th1zjets"), mccut, "goff");


  // ------- Get Single top ------- 
  th1stops = new TH1D("th1stops", "th1stops", bins, dm_min, dm_max);
  th1stopt = new TH1D("th1stopt", "th1stopt", bins, dm_min, dm_max);
  th1stoptw = new TH1D("th1stoptw", "th1stoptw", bins, dm_min, dm_max);
  //th1stops = new TH1D("th1stops", "th1stops", bins, ptbins_boosted);
  //th1stopt = new TH1D("th1stopt", "th1stopt", bins, ptbins_boosted);
  //th1stoptw = new TH1D("th1stoptw", "th1stoptw", bins, ptbins_boosted);
  th1stops->Sumw2();
  th1stopt->Sumw2();
  th1stoptw->Sumw2();
  
  fillMChisto(th1stops,treests,mcobservable,mccut);
  fillMChisto(th1stopt,treestt,mcobservable,mccut);
  fillMChisto(th1stoptw,treestw,mcobservable,mccut);
 
  th1stopps = new TH1D("th1stopps", "th1stopps", bins, dm_min, dm_max);
  th1stoppt = new TH1D("th1stoppt", "th1stoppt", bins, dm_min, dm_max);
  th1stopptw = new TH1D("th1stopptw", "th1stopptw", bins, dm_min, dm_max);
  //th1stopps = new TH1D("th1stopps", "th1stopps", bins, ptbins_boosted);
  //th1stoppt = new TH1D("th1stoppt", "th1stoppt", bins, ptbins_boosted);
  //th1stopptw = new TH1D("th1stopptw", "th1stopptw", bins, ptbins_boosted);
  th1stopps->Sumw2();
  th1stoppt->Sumw2();
  th1stopptw->Sumw2();

  fillMChisto(th1stopps,tree64,mcobservable,mccut);
  fillMChisto(th1stoppt,tree65,mcobservable,mccut);
  fillMChisto(th1stopptw,tree66,mcobservable,mccut);

  // ---- Scale the histos ---- 
  ScaleHistos(channel);

#if 0    
  // ---- Make smooth diboson shape ----------
  TH1D* th1wvclone = (TH1D *)th1wv->Clone("th1wvclone");
  float tmin = 200.0;
  if(channel>1) tmin = 300.0;
  gaus2 = new TF1("gaus2","gaus", tmin, 1000000000.);
  th1wvclone->Fit(gaus2,"I0","");
#endif
    
  // ---- Empty histograms for display/plotting ---- 
  SetupEmptyHistogram(bins, dm_min, dm_max, xtitle);
  //SetupEmptyHistogram(bins, ptbins_boosted, xtitle);
  
  // ---- Sum all backgrounds ----------
  TH1D* th1wv_no_overflow = (TH1D *)th1wv->Clone("th1wv_no_overflow");
  SumAllBackgrounds();


#if 0
  // ---- Get signal histogram ----------
  TCut wwsigratio= GetSigRatioFunction(0.03, 0.0, 0.0,mcobservable);
  TCut wzsigratio= GetSigRatioFunction(0.03, 0.0, 0.0,mcobservable,"wz");
  
  TCut wwsigcut = mccut*wwsigratio;
  TCut wzsigcut = mccut*wzsigratio;
  
  TH1D *signalForDisplay = new TH1D("signalForDisplay","signalForDisplay",bins,dm_min,dm_max);
  //TH1D *signalForDisplay = new TH1D("signalForDisplay","signalForDisplay",bins,ptbins_boosted);
  
  drawstr = TString(mcobservable)+TString(">>signalForDisplay");

  cout <<
    TString("wwtree->Draw(\"")+drawstr+TString("\", \"")+
    TString((const char*)wwsigcut)+TString("\", \"goff\")") << endl;

  th1wwlo->Reset();th1wwhi->Reset();
  fillMChisto(th1wwlo,treewwlo,mcobservable,wwsigcut);
  fillMChisto(th1wwhi,treewwhi,mcobservable,wwsigcut);

  // stitch the two histograms together into one, post-smearing.
  for (int ibin=1;ibin<=signalForDisplay->GetNbinsX();ibin++) {
    if (signalForDisplay->GetBinLowEdge(ibin)<560) {
      signalForDisplay->SetBinContent(ibin,th1wwlo->GetBinContent(ibin));
      signalForDisplay->SetBinError  (ibin,th1wwlo->GetBinError  (ibin));
    } else {
      signalForDisplay->SetBinContent(ibin,th1wwhi->GetBinContent(ibin));
      signalForDisplay->SetBinError  (ibin,th1wwhi->GetBinError  (ibin));
    }
  }

  cout << "signalForDisplay nentries = " << signalForDisplay->GetEntries() << endl;

  // ----- need to subtract the diboson contribution 
  signalForDisplay->SetLineWidth(2);
  signalForDisplay->SetLineColor(1);
  signalForDisplay->SetFillColor(0);

  // not precise, but close enough for gubmint work...
  signalForDisplay->Scale(mc2data_scale*WW_scale);

  //-------- Add overflow bin ----------------
  AddOverflowBin(signalForDisplay);
#endif

  // ---- Compose the stack ----------
  THStack* hs = new THStack("hs","MC contribution");
  //hs->Add(th1zjets);
  //hs->Add(th1qcd);
  hs->Add(th1Top);
  hs->Add(th1wjets);
  hs->Add(th1wv);  // add WW in with backgrounds for display only
  //hs->Add(signalForDisplay);



  // ---- Stack for shape systematics Up ----------
  double bkgd_norm_fracerror = domu ? mu_bkgd_norm_error : el_bkgd_norm_error;
  double sig_norm_fracerror = domu ? mu_sig_norm_error : el_sig_norm_error;

  cout << "Background normalization fractional systematic = " << bkgd_norm_fracerror << endl;
  cout << "Signal     normalization fractional systematic = " << sig_norm_fracerror << endl;

  //TF1* formScaleUp = new TF1("formScaleUp", "1.0+0.4*log(x/5)", dm_min, dm_max);
  //TF1* formScaleDn = new TF1("formScaleDn", "1.0-0.2*log(x/5)", dm_min, dm_max);
  TF1* formScaleUp = new TF1("formScaleUp", Form("1.0+%f",bkgd_norm_fracerror), dm_min, dm_max);
  TF1* formScaleDn = new TF1("formScaleDn", Form("1.0-%f",bkgd_norm_fracerror), dm_min, dm_max);

  systUp = (TH1D*) th1wjets->Clone("systUp");
  systUp->Multiply(formScaleUp);
  //systUp->Add(th1zjets);
  //systUp->Add(th1qcd);
  systUp->Add(th1Top);
  systUp->SetFillColor(0);
  systUp->SetLineStyle(2);
  systUp->SetLineColor(2);
  systUp->SetLineWidth(3);

  // ---- Stack for shape systematics Down ----------
  systDown = (TH1D*) th1wjets->Clone("systDown");
  systDown->Multiply(formScaleDn);
  //systDown->Add(th1zjets);
  //systDown->Add(th1qcd);
  systDown->Add(th1Top);
  systDown->SetFillColor(0);
  systDown->SetLineWidth(3);
  systDown->SetLineStyle(2);
  systDown->SetLineColor(2);

  /////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  

  // ------- Setup the canvas ------- 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadBottomMargin(0.3);
  // gStyle->SetErrorX(0.5);

  TCanvas* c1 = new TCanvas("dijetPt", "", 10,10, 500, 500);
  TPad *d1, *d2;
  c1->Divide(1,2,0,0);
  d1 = (TPad*)c1->GetPad(1);
  d1->SetPad(0.01,0.30,0.95,0.99);
  d2 = (TPad*)c1->GetPad(2);
  d2->SetPad(0.01,0.02,0.95,0.30);

  d1->cd();
  gPad->SetBottomMargin(0.005);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.04);
  // gPad->SetLeftMargin(0.14);


  // Draw it all
  double ymax= 5000000.;
  double ymin= 7.0;
  if(channel>1) { 
    ymax= 3000.;
    ymin= 0.50;
  }

  th1totempty->GetYaxis()->SetRangeUser(ymin, ymax);
  th1data->GetYaxis()->SetRangeUser(ymin, ymax);
  th1totempty->Draw();
  hs->Draw("samehist");
  for (int i=1;i<=th1tot->GetNbinsX();i++)
    {
      double val = th1tot->GetBinContent(i); // / (ptbins_boosted[i]-ptbins_boosted[i-1]);
      double err = fabs(th1tot->GetBinError(i)); // / (ptbins_boosted[i]-ptbins_boosted[i-1]);
      TBox *b = new TBox(th1tot->GetBinLowEdge(i),
			 val-err,th1tot->GetBinLowEdge(i+1),val+err);
      b->SetLineColor(1);
      b->SetFillColor(1);
      b->SetFillStyle(3001);
      b->SetLineStyle(3001);	 
      b->Draw();
    }
  //data2draw->Draw("esame");
  th1data->Draw("esame");
  cmspre(); 
  // Set up the legend
  TLegend* Leg = GetLegend(channel);   
  Leg->Draw();  
  d1->SetLogy();
  gPad->RedrawAxis();


  d2->cd();
  gPad->SetTopMargin(0.02);
  gPad->SetRightMargin(0.04);
  gPad->SetFrameBorderSize(0);
  gPad->SetBottomMargin(0.45);
  gPad->SetTickx();

  th1emptyclone->Draw();
  hhratio->Draw("esame");
  //hhratioUp->Draw("hist lsame");
  //hhratioDown->Draw("hist lsame");
  TLine *line; line = new TLine(dm_min,1.0,dm_max,1.0);
  line->SetLineStyle(1);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->Draw();
  c1->Print(TString("OutDir/")+outfile+TString("_fatjetPt.png"));
  //gPad->WaitPrimitive();
  c1->Modified();
  c1->Update();
  c1->SaveAs(TString("OutDir/")+outfile+TString("_fatjetPt.pdf")); 
  c1->SaveAs(TString("OutDir/")+outfile+TString("_fatjetPt.root")); 

   

  ///// -------------------------------//////

  if(saveDataCards_) {
    outputForLimit->cd();
    th1data->SetName("data_obs");     th1data->Write("data_obs");
    th1tot->SetName("background");    th1tot->Write("background");
    th1ww->SetName("ww");             th1ww->Write("ww");
    th1ww->SetName("wz");             th1wz->Write("wz");
    th1wv->SetName("diboson");        th1wv->Write("diboson");

    th1wv_no_overflow->Write("th1wv_no_overflow");
    char* tempname = "background_backshapeUp";
    if(channel==0) tempname = "background_mudijet_backshapeUp";
    if(channel==1) tempname = "background_eldijet_backshapeUp";
    if(channel==2) tempname = "background_muboosted_backshapeUp";
    if(channel==3) tempname = "background_elboosted_backshapeUp";
    systUp->SetName(tempname);
    systUp->Write(tempname);
    tempname = "background_backshapeDown";
    if(channel==0) tempname = "background_mudijet_backshapeDown";
    if(channel==1) tempname = "background_eldijet_backshapeDown";
    if(channel==2) tempname = "background_muboosted_backshapeDown";
    if(channel==3) tempname = "background_elboosted_backshapeDown";
    systDown->SetName(tempname);
    systDown->Write(tempname);

    h_smearfactor->Write();
    h_sigma_MC->Write();

    outputForLimit->Close();

  } ///// close if saveDataCards_

    //delete th1wvclone;
}

//======================================================================

void makeATGCLimitDataCards(void)
{
  makeATGCLimitDataCards(2);
  //makeATGCLimitDataCards(3);
}

