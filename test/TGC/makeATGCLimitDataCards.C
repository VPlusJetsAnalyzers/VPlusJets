#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TChain.h"
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
#include "TRandom.h"
#include "TMath.h"
#include "TEntryList.h"

#include "TFitResult.h"
#include "TFitResultPtr.h"

#undef ATGCISSIGNAL

//#define SIGNORM_FROM_FIT
#undef SIGNORM_FROM_FIT

const bool blinded      = false;

const float lz4display=0.01;
const float dkg4display=0.0;

// for  asym. binning...
//const double ptbins_boosted[15] = {
//  200,225,250,275,300,325,350,375,400,425,450,500,550,625,925
//};

const double intLUMIinvpb_mu = 19300.; // for normalization of signal
const double intLUMIinvpb_el = 19200.; // (background norm comes from data)

//const double WJets_scale   = 37509.0    * intLUMIinvpb/(18353019+50768992);
//const double WJets_scale   = 228.9*1.3  * intLUMIinvpb/8955318;
//const double WJets_scale   =  26.4*1.3  * intLUMIinvpb/9492452;
//const double ZJets_scale   =  3503.71   * intLUMIinvpb/30209426;

const double WJ_scale_lo   =    23.5    / 9492452; // pt180 samp xsec from PREP
const double WJ_scale_hi   =    0.0587  / 155501;  // pt600 samp

//double WW_scale_lo   =    33.61     / 9450414; // xsec from PREP
//double WW_scale_hi   =    0.005235  / 1000139; // xsec from PREP

// aMC@NLO samples produced privately:
const double WW_scale_lo   =    0.2222038 * 2 * 1.883      / 119692;
const double WW_scale_hi   =    0.2222038 * 2 * 1.269e-02  / 118889;
const double WZ_scale_lo   =    0.230996  *     1.00125782 / (58005+95230);
const double WZ_scale_hi   =    0.230996  *     0.00754    / 118889;

const double WW_scale_NLO  =    1.0; // 57.2;             // xsec from MCFM paper (arxiv:1105.0020v1)
const double WZ_scale_NLO  =    1.0; // 22.88 / 10000267; // xsec from MCFM paper - on-shell Z, no Wgam* component

const double stitchwwgev = 600;
const double stitchwzgev = 520;
const double stitchwjgev = 680;

// For WV, take NLO prediction as opposed to our fit results (smaller uncertainty)
//double WWplusWZscale =    80.1 * intLUMIinvpb;

const double ttbar_scale   =   225.197 /20975917;
const double mtt700_scale  =    15.38  / 3082808;
const double mtt1k_scale   =     3.063 / 1249109;
// don't bother with these, only tW contributes in this phase space
//const double SToppS_scale  =     1.76  /  139974;
//const double SToppT_scale  =    30.7   / 1935066;
//const double STopS_scale   =     3.79  /  259960;
//const double STopT_scale   =    56.4   / 3758221;
const double STopTbarW_scale =    11.1   /  493458;
const double STopTW_scale    =    11.1   /  497657;


// NORMALIZE TO THE CROSS-SECTION FIT YIELD RESULTS IN THE SIGNAL REGION:
//

//Muons (boosted):
const double WV_fityield_muboosted    =  204.118392963; // +/- 116.455352775
const double top_fityield_muboosted   =  450.117983309; // +/-  36.0567046084
const double WJets_fityield_muboosted = 1318.48178491;  // +/-  66.6457513819

//Electrons (boosted):
const double WV_fityield_elboosted    =  285.019244526; // +/- 107.560258373
const double top_fityield_elboosted   =  363.719433478; // +/-  36.2755682802
const double WJets_fityield_elboosted = 1023.34933926;  // +/-  60.6034935714

// Yield uncertainties for sum(W+jets,top)
// computed with correlations from the fit,
// and divided by the total background

const double mu_bkgd_norm_error = (77.65)/(WJets_fityield_muboosted+top_fityield_muboosted);
const double el_bkgd_norm_error = (70.76)/(WJets_fityield_elboosted+top_fityield_elboosted);

#ifdef SIGNORM_FROM_FIT
// Taking signal normalization and norm unc. from the fit:
const double mu_sig_norm_error = 112.25/376.65;
const double el_sig_norm_error = 95.265/410.23;
#else
// Add absolute errors for WW, W+Z, and W-Z NLO from MCFM paper in quadrature
// (arxiv:1105.0020v1)
//
const double mu_sig_norm_error =  sqrt((.041 * .041 * 57.25 *57.25)+
				       (.052 * .052 * 14.48 *14.48)+
				       (.054 * .054 * 8.4   *  8.4)
				       )/(57.25+14.48+8.4);
const double el_sig_norm_error =  mu_sig_norm_error;
#endif

const float yRatioMin = 0.2;
const float yRatioMax = 2.2;

bool domu=true;
bool injectSignal=false;
double  mc2data_scale;

//TFile* sigLambdaZ;

////////// ALL input trees ///////////
TChain* treedata;
TChain* treewwlo;
TChain* treewwhi;
TChain* treewzlo;
TChain* treewzhi;
TChain* treewjlo;
TChain* treewjhi;
TChain* treewjsherpa;
TChain* treettblo,*treettbmd,*treettbhi;
//TChain* treeqcd;
TChain* treezj;
//TChain* treestops;
//TChain* treestopt;
//TChain* treestopstbar;
//TChain* treestopttbar;
TChain* treestoptw;
TChain* treestoptbarw;


const int    bins   = 15;
const double dm_min = 200.; 
const double dm_max = 800.;
const double binwidth = (dm_max-dm_min)/bins;

class histos_t {
public:
  TH1D *hobs;
  TH1D *hcum;

  explicit histos_t() {}
  explicit histos_t(const TString& name,bool ensumw2=true);
  ~histos_t() {}
  inline void Scale(double scale)  { hobs->Scale(scale); hcum->Scale(scale); }
  inline void SetFillColor(int fc) { hobs->SetFillColor(fc); hcum->SetFillColor(fc); }
  inline void SetFillStyle(int fs) { hobs->SetFillStyle(fs); hcum->SetFillStyle(fs); }
  inline void SetLineColor(int lc) { hobs->SetLineColor(lc); hcum->SetLineColor(lc); }
  inline void SetLineWidth(int lw) { hobs->SetLineWidth(lw); hcum->SetLineWidth(lw); }
  inline void SetLineStyle(int ls) { hobs->SetLineStyle(ls); hcum->SetLineStyle(ls); }
  inline double Integral(void)     { return hobs->Integral(); }
  inline void SetMarkerStyle(int ms)   { hobs->SetMarkerStyle(ms); hcum->SetMarkerStyle(ms); }
  inline void SetMarkerSize(double ms) { hobs->SetMarkerSize(ms); hcum->SetMarkerSize(ms); }
  inline void SetMinimum(double min)   { hobs->SetMinimum(min); hcum->SetMinimum(min); }

  TAxis *GetXaxis(void)            { return hobs->GetXaxis(); }

  inline void Add(const histos_t *ad, double c1=1) {
    hobs->Add(ad->hobs, c1); hcum->Add(ad->hcum, c1);
  }
  inline void Add(const histos_t *ad1, const histos_t *ad2,
		  double c1=1, double c2=1) {
    hobs->Add(ad1->hobs, ad2->hobs, c1, c2);
    hcum->Add(ad1->hcum, ad2->hcum, c1, c2);
  }
  inline void Divide(const histos_t *denom) {
    hobs->Divide(denom->hobs);
    hcum->Divide(denom->hcum);
  }
  histos_t *Clone(const TString& name) {
    histos_t *clone = new histos_t();
    clone->hobs = (TH1D *)this->hobs->Clone(name);
    clone->hcum = (TH1D *)this->hcum->Clone(name+"_cum");
    return clone;
  }

  void Sumw2() { hobs->Sumw2(); hcum->Sumw2(); }

  void fillFromTree(TChain *thetree,
		    const TString& theobservable,
		    const TCut& thecut,
		    const TString& wtstr="",
		    double scalefactor=1.,
		    double locutoff = 0,
		    double hicutoff = 9.99e99);
  
private:
};


//======================================================================

histos_t::histos_t(const TString& name,bool ensumw2)
{
  hobs = new TH1D(name,name,bins,dm_min,dm_max);
  hcum = new TH1D(name+"_cum",name+"_cum",bins,dm_min,dm_max);

  if (ensumw2) {
    hobs->Sumw2();
    hcum->Sumw2();
  }
}

//======================================================================

void
histos_t::fillFromTree(TChain *thetree,
		       const TString& theobservable,
		       const TCut& thecut,
		       const TString& wtstr,
		       double scalefactor,
		       double locutoff,
		       double hicutoff)
{
  cout << "filling  "<<hobs->GetName()<<", "<< hcum->GetName()<<endl;

  TCut newcut = thecut;
  if (locutoff > 0)   newcut = TCut(TString("(") + newcut.GetTitle() + Form("&& (%s > %g))",theobservable.Data(),locutoff));
  if (hicutoff < 1e4) newcut = TCut(TString("(") + newcut.GetTitle() + Form("&& (%s < %g))",theobservable.Data(),hicutoff));

  cout << "Determining event list from cut " << newcut.GetTitle() << endl;

  thetree->Draw(">>elist",newcut,"entrylist");
  TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
  thetree->SetEntryList(elist);

  elist->Print();

  TString newwtstr = Form("%g",scalefactor);
  if (wtstr.Length())
    newwtstr += Form("*%s",wtstr.Data());

  cout << 
    thetree->GetName()+
    TString("->Draw(")+
    theobservable +
    Form(">>+%s",hobs->GetName()) +
    TString(",") +
    newwtstr +  
    TString(",goff);")
       << endl;

  thetree->Draw(theobservable+Form(">>+%s",hobs->GetName()),newwtstr,"goff");

  cout << hobs->GetName() << " now has " << hobs->GetEntries() << " entries." << endl;

#if 1
  // for cumulative histogram,
  // ttree::Draw fill histo for appropriate bin for obs value + all bins below, down to the min
  //
  for (int ibin=0; ibin<=bins; ibin++) {
    TString cumobsstr(Form("TMath::Min(%g,",dm_max)+theobservable+Form("-%g)",ibin*binwidth));
    TString cumcut = newwtstr +TString("*") + TString("((") + theobservable + Form("-%g)>= %g)",ibin*binwidth,dm_min);
    //if (ibin<2)
    cout <<
      thetree->GetName() +
      TString("->Draw(") + 
      cumobsstr +
      Form(">>%s",hcum->GetName()) + 
      TString(",") +
      cumcut + 
      TString(",goff);")
	 <<endl;
    thetree->Draw(cumobsstr+Form(">>+%s",hcum->GetName()),cumcut,"goff");
  }

  cout << hcum->GetName() << " now has " << hcum->GetEntries() << " entries." << endl;
#endif

  thetree->SetEntryList(0);

}                                              // histos_t::fillFromTree

//======================================================================

////////// ALL histograms ///////////

histos_t *th1data;
histos_t *th1ww;
histos_t *th1wz;
histos_t *th1wv;
histos_t *th1wjets;
histos_t *th1ttbar;
//histos_t *th1stops;
//histos_t *th1stopt;
//histos_t *th1stopps;
//histos_t *th1stoppt;
histos_t *th1stoptw;
histos_t *th1stoptbarw;
histos_t *th1Top;
histos_t *th1tot;
histos_t *hhratio;
histos_t *wwatgc4Display;
histos_t *wzatgc4Display;

TH1D *th1totempty;
TH1D* th1emptyclone;
//TH1F* hhratioUp;
//TH1F* hhratioDown;

bool saveDataCards_ = !blinded;

TF1 *gaus2;

//======================================================================


void InstantiateTrees() {

  treedata  = new TChain("WJet");

  treewwlo  = new TChain("WJet");
  treewwhi  = new TChain("WJet");
  treewzlo  = new TChain("WJet");
  treewzhi  = new TChain("WJet");
  treewjlo  = new TChain("WJet");
  treewjhi  = new TChain("WJet");
  treettblo = new TChain("WJet");
  treettbmd = new TChain("WJet");
  treettbhi = new TChain("WJet");
  //treeqcd = new TChain("WJet");
  //treezj  = new TChain("WJet");
  //treests   = new TChain("WJet");
  //treestt   = new TChain("WJet");
  //tree64    = new TChain("WJet");
  //tree65    = new TChain("WJet");
  treestoptw    = new TChain("WJet");
  treestoptbarw = new TChain("WJet");

  treewjsherpa = new TChain("WJet");

  //// ------------ Get all trees
 if (domu) {
   treedata->Add("InData/RD_WmunuJets_DataAll_GoldenJSON_v1.root");
   treedata->Add("InData/RD_WmunuJets_DataAll_GoldenJSON_v2.root");
   treewwlo->Add("InMC/RD_mu_WW_minPt150_amcnlo_CMSSW532.root");
   treewwhi->Add("InMC/RD_mu_WW_minPt500_amcnlo_CMSSW532.root");
   treewzlo->Add("InMC/RD_mu_WZ_minPt150_amcnlo_CMSSW532.root");
   treewzlo->Add("InMC/RD_mu_WZ_minPt150_amcnlo_add_CMSSW532.root"); // N.B.
   treewzhi->Add("InMC/RD_mu_WZ_minPt500_amcnlo_CMSSW532.root");
// treewj->Add("InMC/RD_mu_WpJ_CMSSW532.root");
   treewjlo->Add("InMC/RD_mu_WpJ_PT180_Madgraph_CMSSW532.root");
   treewjhi->Add("InMC/RD_mu_WpJ_minPt600_CMSSW532.root");

   treewjsherpa->Add("InMC/RD_mu_WJets_sherpa_CMSSW532.root");

// treettbar->Add("InMC/RD_mu_TTbar_CMSSW532.root");
   treettblo->Add("InMC/RD_mu_TTJets_poheg_CMSSW532.root");
   treettbmd->Add("InMC/RD_mu_Mtt700to1000_CMSSW532.root");
   treettbhi->Add("InMC/RD_mu_Mtt1000toinf_CMSSW532.root");
// treeqcd->Add("InMC/RD_mu_WpJ_PT180_Madgraph_CMSSW532.root"); // use WpJ, shapes are the same
// treezj->Add("InMC/RD_mu_ZpJ_CMSSW532.root");
// treests->Add("InMC/RD_mu_STopS_T_CMSSW532.root");
// treestt->Add("InMC/RD_mu_STopT_T_CMSSW532.root");
// tree64->Add("InMC/RD_mu_STopS_Tbar_CMSSW532.root");
// tree65->Add("InMC/RD_mu_STopT_Tbar_CMSSW532.root");
   treestoptw->Add("InMC/RD_mu_STopTW_T_CMSSW532.root");
   treestoptbarw->Add("InMC/RD_mu_STopTW_Tbar_CMSSW532.root");

 } else { // electrons

   treedata->Add("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_v1.root");
   treedata->Add("InData/RD_WenuJets_DataAllSingleElectronTrigger_GoldenJSON_v2.root");
   treewwlo->Add("InMC/RD_el_WW_minPt150_amcnlo_CMSSW532.root");
   treewwhi->Add("InMC/RD_el_WW_minPt500_amcnlo_CMSSW532.root");
   treewzlo->Add("InMC/RD_el_WZ_minPt150_amcnlo_CMSSW532.root");
   treewzlo->Add("InMC/RD_el_WZ_minPt150_amcnlo_add_CMSSW532.root"); // N.B.
   treewzhi->Add("InMC/RD_el_WZ_minPt500_amcnlo_CMSSW532.root");
// treewj->Add("InMC/RD_el_WpJ_CMSSW532.root");
   treewjlo->Add("InMC/RD_el_WpJ_PT180_Madgraph_CMSSW532.root");
   treewjhi->Add("InMC/RD_el_WpJ_minPt600_CMSSW532.root");

   treewjsherpa->Add("InMC/RD_el_WJets_sherpa_CMSSW532.root");

// treettbar->Add("InMC/RD_el_TTbar_CMSSW532.root");
   treettblo->Add("InMC/RD_el_TTJets_poheg_CMSSW532.root");
   treettbmd->Add("InMC/RD_el_Mtt700to1000_CMSSW532.root");
   treettbhi->Add("InMC/RD_el_Mtt1000toinf_CMSSW532.root");
// treeqcd->Add("InMC/"InData/RDQCD_WenuJets_Isog0p3NoElMVA_19p2invfb.root");
// treezj->Add("InMC/RD_el_ZpJ_CMSSW532.root");
// treests->Add("InMC/RD_el_STopS_T_CMSSW532.root");
// treestt->Add("InMC/RD_el_STopT_T_CMSSW532.root");
// tree64->Add("InMC/RD_el_STopS_Tbar_CMSSW532.root");
// tree65->Add("InMC/RD_el_STopT_Tbar_CMSSW532.root");
   treestoptw->Add("InMC/RD_el_STopTW_T_CMSSW532.root");
   treestoptbarw->Add("InMC/RD_el_STopTW_Tbar_CMSSW532.root");
 }

  double nData = treedata->GetEntries();
  std::cout << "ndata =" << nData <<std::endl;

#if 0
  //// ------------ Create a tree branch for dijet pt
  const char* dijetPt = "sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))";
  treedata->SetAlias("dijetPt", dijetPt);
  treewwlo->SetAlias("dijetPt", dijetPt);
  treewwhi->SetAlias("dijetPt", dijetPt);
  treewzlo->SetAlias("dijetPt", dijetPt);
  treewzhi->SetAlias("dijetPt", dijetPt);
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
#endif
}                                                    // InstantiateTrees

//======================================================================
// This is a TH1 method in ROOT 6.x

TH1 *GetCumulative(TH1 *src)
{
  double sum=0;
  TH1* hintegrated = (TH1*)src->Clone(TString(src->GetName()) + "_cum");
  hintegrated->Reset();

  const Int_t nbinsx = src->GetNbinsX();
  for (Int_t binx = nbinsx; binx >= 1; --binx) {
    const Int_t bin = hintegrated->GetBin(binx);
    sum += src->GetBinContent(bin);
    hintegrated->SetBinContent(bin, sum);
  }
  return hintegrated;
}

//======================================================================

void ScaleHistos(int channel)
{
#if 0
  // Print all bin contents prior to scaling:
  printf ("%-5s %7s %10s %10s %9s %11s %9s\n","bin","Pt(GeV)","TTbar","WJets","WW","WZ","Data");
  for (int ibin=1; ibin<=bins+1; ibin++) {
    TAxis *xax = thwzlo->GetXaxis();

    if (ibin > bins)
      printf ("%-4d ovrflow", ibin);
    else
      printf ("%-4d %3.0f-%3.0f", ibin,
	      xax->GetBinLowEdge(ibin),
	      xax->GetBinUpEdge(ibin));
    printf (" %5.4g+-%4.1f",  th1toplo->hobs->GetBinContent(ibin)+
	                      th1topmd->hobs->GetBinContent(ibin)+
	                      th1tophi->hobs->GetBinContent(ibin),
	                      th1toplo->hobs->GetBinError(ibin)+
	                      th1topmd->hobs->GetBinError(ibin)+
	                      th1tophi->hobs->GetBinError(ibin));
    if (xax->GetBinLowEdge(ibin) < stitchwjgev)
      printf (" %5.4g+-%4.1f", th1wjlo->hobs->GetBinContent(ibin), th1wjlo->hobs->GetBinError(ibin));
    else
      printf (" %5.4g+-%4.1f", th1wjhi->hobs->GetBinContent(ibin), th1wjhi->hobs->GetBinError(ibin));

    if (xax->GetBinLowEdge(ibin) < stitchwwgev)
      printf (" %7.1f+-%4.1f", th1wwlo->hobs->GetBinContent(ibin), th1wwlo->hobs->GetBinError(ibin));
    else
      printf (" %7.1f+-%4.1f", th1wwhi->hobs->GetBinContent(ibin), th1wwhi->hobs->GetBinError(ibin));

    if (xax->GetBinLowEdge(ibin) < stitchwzgev)
      printf (" %7.1f+-%4.1f", th1wzlo->hobs->GetBinContent(ibin), th1wzlo->hobs->GetBinError(ibin));
    else
      printf (" %7.1f+-%4.1f", th1wzhi->hobs->GetBinContent(ibin), th1wzhi->hobs->GetBinError(ibin));

    printf   (" %5.3g+-%4.1f", th1data->hobs->GetBinContent(ibin), th1data->hobs->GetBinError(ibin));
    printf   ("\n");
  }
#endif

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

  //th1stops->SetFillColor(7);    th1stops->SetLineWidth(2);
  //th1stopt->SetFillColor(13);   th1stopt->SetLineWidth(2);
  //th1stopps->SetFillColor(7);   th1stopps->SetLineWidth(2);
  //th1stoppt->SetFillColor(13);  th1stoppt->SetLineWidth(2);
  th1stoptw->SetFillColor(9);     th1stoptw->SetLineWidth(2);
  th1stoptbarw->SetFillColor(9);  th1stoptbarw->SetLineWidth(2);

  // Add all the top together:
  th1Top->Add(th1ttbar,1);
  th1Top->Add(th1stoptw,1);
  th1Top->Add(th1stoptbarw,1);
  //th1Top->Add(th1stopt,1);
  //th1Top->Add(th1stops,1);
  //th1Top->Add(th1stoppt,1);
  //th1Top->Add(th1stopps,1);

  th1Top->SetFillColor(kGreen+2);
  th1Top->SetLineColor(kGreen+2);
  th1Top->SetLineWidth(0);

  // for one-time test:
  // TF1* formScaleUp = new TF1("formScaleUp", "1.0-0.0002*x", dm_min, dm_max);
  // th1wjets->Multiply(formScaleUp);

  //th1wjets->Scale(37509.0 * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el)/93158078); // sherpa sample

  th1wjets->SetFillColor(kRed);
  th1wjets->SetLineColor(kRed);
  th1wjets->SetLineWidth(0);

  th1ww->Scale(WW_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));
  th1wz->Scale(WZ_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));

  std::cout << "wjets " << th1wjets->hobs->Integral() << std::endl;
  std::cout << "top "   << th1Top->hobs->Integral()   << std::endl;
  std::cout << "ww "    << th1ww->hobs->Integral()    << std::endl;
  std::cout << "wz "    << th1wz->hobs->Integral()    << std::endl;
  std::cout <<" data "  << th1data->hobs->Integral()  << std::endl;

  th1wv->Add(th1ww,1);
  th1wv->Add(th1wz,1);
  th1wv->SetFillColor(kAzure+8);
  th1wv->SetLineColor(kAzure+8);
  th1wv->SetLineWidth(0);

#ifdef SIGNORM_FROM_FIT
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
  for (int ibin=1; ibin<=bins+1; ibin++) {
    if (ibin>bins)
      printf ("%-4d ovrflow", ibin);
    else
      printf ("%-4d %3.0f-%3.0f", ibin,
	      th1wv->GetXaxis()->GetBinLowEdge(ibin),
	      th1wv->GetXaxis()->GetBinUpEdge(ibin));
    printf (" %5.4g+-%4.2f",  th1Top->hobs->GetBinContent(ibin),  th1Top->hobs->GetBinError(ibin));
    printf (" %5.3g+-%4.1f",th1wjets->hobs->GetBinContent(ibin),th1wjets->hobs->GetBinError(ibin));
    printf (" %5.2f+-%4.2f",   th1wv->hobs->GetBinContent(ibin),   th1wv->hobs->GetBinError(ibin));
    printf (" %5.3g+-%4.1f", th1data->hobs->GetBinContent(ibin), th1data->hobs->GetBinError(ibin));
    printf ("\n");
  }
 
  //mc2data_scale = th1data->Integral()/(top_fityield+WV_fityield+WJets_fityield);

}                                                         // ScaleHistos


//======================================================================

void AddOverflowBin(TH1 *hist,bool calcerror=true) {
  int nBinsTot = hist->GetNbinsX();
  double lastbin = hist->GetBinContent(nBinsTot);
  double overflow = hist->GetBinContent(nBinsTot+1);
  hist->SetBinContent(nBinsTot, lastbin+overflow);
  if (calcerror) hist->SetBinError(nBinsTot, sqrt(lastbin+overflow));
}

//======================================================================

// Sum all the backgrounds
void SumAllBackgrounds() {

  //-------- First add overflow bin ----------------
  AddOverflowBin(th1ww->hobs);
  AddOverflowBin(th1wz->hobs);
  AddOverflowBin(th1wv->hobs);

  AddOverflowBin(th1wjets->hobs);
  //AddOverflowBin(th1zjets);
  AddOverflowBin(th1Top->hobs);
  //AddOverflowBin(th1qcd);
  AddOverflowBin(th1data->hobs,false);

  //-------- Now sum of all bkg histograms ----------
  th1tot->hobs = (TH1D*)th1wjets->hobs->Clone("th1tot");     th1tot->hobs->Reset();
  th1tot->hcum = (TH1D*)th1wjets->hcum->Clone("th1totcum");  th1tot->hcum->Reset();

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

}                                                   // SumAllBackgrounds

//======================================================================

void makeRatioHisto()
{
#if 0
  //-------- Needed for plotting ----------
  TH1D* th1totClone = ( TH1D*) th1tot->Clone("th1totClone");

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
#endif

  //-------- Ratio histogram ----------
  hhratio    = th1data->Clone("hhratio")  ;
  hhratio->Sumw2();
  hhratio->SetMarkerStyle(20);
  hhratio->SetMarkerSize(1.25);
  hhratio->SetLineWidth(2);
  hhratio->hobs->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
  hhratio->hcum->GetYaxis()->SetRangeUser(yRatioMin, yRatioMax);
  hhratio->Divide(th1tot);
#if 0
  double binError(0.0), mcbinentry(0.0), mcerror(0.0);
  for(int i=0; i<hhratio->GetNbinsX(); ++i) {
    binError   = hhratio->GetBinError(i);
    mcerror    = th1tot->GetBinError(i);
    mcbinentry = th1tot->GetBinContent(i);
    if(mcbinentry>0.) mcerror /= mcbinentry;
    else mcerror = 0.0;
    binError = sqrt(binError**2 + mcerror**2);
    hhratio->SetBinError(i, binError);
    if(th1data->GetBinContent(i)<0.1) hhratio->SetBinError(i,0.0);
  }
#endif
}                                                      // makeRatioHisto

//======================================================================

TLegend* GetLegend(int channel)
{
  // float  legX0=0.5, legX1=0.89, legY0=0.41, legY1=0.86;
  float  legX0=0.65, legX1=0.96, legY0=0.50, legY1=0.90;

  TLegend * Leg = new TLegend( legX0, legY0, legX1, legY1);
  Leg->SetFillColor(0);
  Leg->SetFillStyle(0);
  Leg->SetTextSize(0.045);
  if (injectSignal)
    Leg->AddEntry(th1data->hobs, domu ?  "Sig. inject toy mu data" : "Sig. inject toy el. data",  "PLE");
  else
    Leg->AddEntry(th1data->hobs, domu ?  "Muon Data" : "Electron Data",  "PLE");

  Leg->AddEntry(th1wv->hobs,  "SM WW+WZ ",  "f");
  Leg->AddEntry(th1wjets->hobs,  "W+jets",  "f");
  Leg->AddEntry(th1Top->hobs,  "top",  "f");
  //if(channel == 1) Leg->AddEntry(th1qcd,  "Multijet",  "f");
  //Leg->AddEntry(th1zjets,  "Z+Jets",  "f");
  Leg->AddEntry(th1tot->hobs,  "MC error",  "f");
  //Leg->AddEntry(systUp,  "Shape error",  "f");
#if 0
  Leg->AddEntry(wwatgc4Display,
		Form("WW, (#lambda_{Z}, #Delta#kappa_{#gamma})=(%g,%g)",
		     lz4display,dkg4display), "l");
  Leg->AddEntry(wzatgc4Display,
		Form("WZ, (#lambda_{Z}, #Delta#kappa_{#gamma})=(%g,%g)",
		     lz4display,dkg4display), "l");
#else
  Leg->AddEntry(wwatgc4Display->hobs,Form("WW, #lambda_{Z} = %g",lz4display), "l");
  Leg->AddEntry(wzatgc4Display->hobs,Form("WZ, #lambda_{Z} = %g",lz4display), "l");
#endif
  Leg->SetFillColor(0);

  return Leg;
}


//======================================================================

void SetupEmptyHistogram(char* xtitle)
//void SetupEmptyHistogram(char* xtitle)
{
  //th1totempty = new TH1D("th1totempty", "th1totempty", bins, ptbins); // dm_min, dm_max);
  th1totempty = new TH1D("th1totempty", "th1totempty", bins,dm_min, dm_max);
  char tmpc[100];    sprintf(tmpc,"Events / %d GeV", (int) (dm_max-dm_min)/bins);
  th1totempty->SetYTitle(tmpc);
  th1totempty->GetYaxis()->SetTitleOffset(1);
  th1totempty->GetYaxis()->SetLabelOffset(0.01);
  th1totempty->GetYaxis()->SetLabelSize(0.08);
  th1totempty->GetYaxis()->SetTitleSize(0.08);

  int maxbin = th1data->hobs->GetMaximumBin();
  float maxval = th1data->hobs->GetBinContent(maxbin);
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

TString
GetSigRatioFunction(const char *obsrvbl, const char *wworwz="ww")
{
  //printf("%g\t%g\t%g\r",lambdaZ,dkappaGamma,deltaG1);

  TFile f("ATGC_shape_coefficients.root");

  TString sigratiostr;
  TString closestr;

  float lambdaZ = lz4display;
  float dkg     = dkg4display;

  for (int i=0; i<=6; i++) {
    TString pname(Form("%s_p%d_lambda_dkg",wworwz,i));
    TProfile2D *prof = (TProfile2D*) f.Get(pname);
    if (!prof) continue;
    if (i) {
      sigratiostr += TString::Format("%s*(%g+",obsrvbl,prof->Interpolate(lambdaZ,dkg));
      closestr += TString(")");
    } else {
      //#ifdef ATGCSIGNAL // the ATGC specific portion is stacked on top of the SM contribution
      sigratiostr += TString::Format("%g-1.0+",prof->Interpolate(lambdaZ,dkg));
      //#else
      //sigratiostr += TString::Format("%g+",    prof->Interpolate(lambdaZ,dkg));
      //#endif
    }
  }
  sigratiostr += TString("0")+closestr;
  //printf("Looking up coefficients for %s lambdaZ=%g, dkappaGamma=%g\n",wworwz,lambda,dkg);

  cout << "sigratio function: " << sigratiostr << endl;

  return sigratiostr;
}                                                 // GetSigRatioFunction

//======================================================================

void cmspre()
{
  TLatex latex;
  latex.SetNDC();

  latex.SetTextSize(0.05);
  latex.SetTextFont(61); // helvetica bold
  latex.SetTextAlign(11); // align left
  latex.DrawLatex(0.15,0.93,"CMS");

  latex.SetTextSize(0.04);
  latex.SetTextFont(52); // helvetica italic
  latex.DrawLatex(0.25,0.93,"Preliminary");

  latex.SetTextAlign(31); // align right
  latex.SetTextSize(0.037);
  latex.SetTextFont(42); // helvetica italic
  latex.DrawLatex(0.94,0.93,Form("%0.1f#kern[0.2]{fb}^{-1} (8 TeV)",
				 (domu ? intLUMIinvpb_mu : intLUMIinvpb_el)/1000.));
}

//======================================================================

void makeSignalHistos4Display(const TString& mcobservable,
			      const TCut&    mccut,
			      const TString& wtstr)
{
  // ---- Get signal histogram ----------
  TString wwsigratio= GetSigRatioFunction(mcobservable);
  TString wzsigratio= GetSigRatioFunction(mcobservable,"wz");

  wwatgc4Display = new histos_t("th1wwatgc");
  wzatgc4Display = new histos_t("th1wzatgc");

  TString wwwtstr = TString("(")+wwsigratio+TString(")");
  TString wzwtstr = TString("(")+wzsigratio+TString(")");
  if (wtstr.Length()) {
    wwwtstr += TString("*")+wtstr;
    wzwtstr += TString("*")+wtstr;
  }

  wwatgc4Display->fillFromTree(treewwlo,mcobservable,mccut,wwwtstr,WW_scale_lo,0,stitchwwgev);
  wzatgc4Display->fillFromTree(treewzlo,mcobservable,mccut,wzwtstr,WZ_scale_lo,0,stitchwzgev);

  wwatgc4Display->fillFromTree(treewwhi,mcobservable,mccut,wwwtstr,WW_scale_hi,stitchwwgev);
  wzatgc4Display->fillFromTree(treewzhi,mcobservable,mccut,wzwtstr,WZ_scale_hi,stitchwzgev);

  wwatgc4Display->Scale(WW_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));
  wzatgc4Display->Scale(WZ_scale_NLO * (domu ? intLUMIinvpb_mu : intLUMIinvpb_el));

  cout << "wwatgc4Display nentries = " << wwatgc4Display->hobs->GetEntries() << endl;

  wwatgc4Display->SetLineWidth(2);
  wwatgc4Display->SetLineColor(1);
  wwatgc4Display->SetFillColor(0);

  wzatgc4Display->SetLineWidth(2);
  wzatgc4Display->SetLineColor(1);
  wzatgc4Display->SetFillColor(0);
  wzatgc4Display->SetLineStyle(2);

  //-------- Add overflow bin ----------------
  AddOverflowBin(wwatgc4Display->hobs);
  AddOverflowBin(wzatgc4Display->hobs);
}                                            // makeSignalHistos4Display

//======================================================================

void DrawItAll(int channel, 
	       THStack *hs, 
	       TH1 *htot,
	       TH1 *hdata,
	       TH1 *hratio,
	       const TString& outfile)
{
  // ------- Setup the canvas ------- 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  // gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadBottomMargin(0.3);
  // gStyle->SetErrorX(0.5);

  TCanvas* c1 = new TCanvas(hs->GetName(), "", 10,10, 500, 500);

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
    ymin= 0.10;
  }

  th1totempty->GetYaxis()->SetRangeUser(ymin, ymax);

  th1totempty->Draw();
  hs->Draw("samehist");
  for (int i=1;i<=bins;i++)
    {
      double val = htot->GetBinContent(i); // / (ptbins_boosted[i]-ptbins_boosted[i-1]);
      double err = fabs(htot->GetBinError(i)); // / (ptbins_boosted[i]-ptbins_boosted[i-1]);
      TBox *b = new TBox(htot->GetBinLowEdge(i),
			 val-err,htot->GetBinLowEdge(i+1),val+err);
      b->SetLineColor(1);
      b->SetFillColor(1);
      b->SetFillStyle(3001);
      b->SetLineStyle(3001);	 
      b->Draw();
    }
  //data2draw->Draw("esame");

  hdata->SetMarkerStyle(20);
  hdata->SetMarkerSize(1.25);
  hdata->SetLineWidth(2);
  hdata->SetMinimum(0.0);
  hdata->GetYaxis()->SetRangeUser(ymin, ymax);
  hdata->SetBinErrorOption(TH1::kPoisson);

  cout << "Bin error option = " << hdata->GetBinErrorOption() << endl;
  cout << "fSumw2.fN = " << hdata->GetSumw2N() << endl;

  hdata->Draw("E0 same");

  cout << "Bin error option = " << hdata->GetBinErrorOption() << endl;
  cout << "fSumw2.fN = " << hdata->GetSumw2N() << endl;

  cmspre(); 
  // Set up the legend
  TLegend* Leg = GetLegend(channel);   
  Leg->Draw();  
  gPad->SetLogy();
  gPad->RedrawAxis();

  d2->cd();
  gPad->SetTopMargin(0.02);
  gPad->SetRightMargin(0.04);
  gPad->SetFrameBorderSize(0);
  gPad->SetBottomMargin(0.45);
  gPad->SetTickx();

  th1emptyclone->Draw();

  hratio->Draw("E0 same");
  //hhratioUp->Draw("hist lsame");
  //hhratioDown->Draw("hist lsame");

#if 1
  {
    TF1 *f1 = new TF1("f1", "pol1", dm_min, dm_max);
    f1->SetParameters(1.0,0.0);
    f1->SetParameters(1,0.0);
    //f1->FixParameter(0,1);
    //f1->FixParameter(1,0);
    //cout<<" par1   "f1->GetParameter(0)<<endl;
    //cout<<" par2   "f1->GetParameter(1)<<endl;

    //TFitResultPtr r =  hratio->Fit("f1", "RBS");
    //TFitResultPtr r = hratio->Fit(myFunc,"S");
    //r->Print("V");     // print full information of fit including covariance matrix
  }
#endif

  TLine *line; line = new TLine(dm_min,1.0,dm_max,1.0);
  line->SetLineStyle(1);
  line->SetLineWidth(1);
  line->SetLineColor(1);
  line->Draw();

  //gPad->WaitPrimitive();
  c1->Modified();
  c1->Update();

  TString outfilename = TString("OutDir/")+outfile;

  if (injectSignal)
    outfilename += TString("_fatjetPt");
  else if (blinded)
    outfilename += TString("_fatjetPt_blinded");
  else
    outfilename += TString("_fatjetPt_unblinded");

  if (TString(hs->GetName()).Contains("cum"))
    outfilename += TString("_cum");

  c1->Print (outfilename + TString(".png"));
  c1->SaveAs(outfilename + TString(".pdf"));
  c1->SaveAs(outfilename + TString(".C"));
}                                                           // DrawItAll   

//======================================================================

//////---------- channel: 0==muon dijet, 1== electron dijet
/////                     2==muon  boosted,   3== electron boosted
void makeATGCLimitDataCards(int channel) {

//   const Int_t bins = 8; 
//   const Float_t dm_min = 200.; 
//   const Float_t dm_max = 600.;

//  Int_t bins = 7; 
//  Float_t dm_min = 100.; 
//  Float_t dm_max = 275.;

  if (channel < 0) { injectSignal = true; channel = abs(channel); }

  if (channel < 2) {
    cout << "Script does not currently handle unboosted channels - fix me!" << endl;
    exit(-1);
  }

  if (injectSignal)
    cout << "*** SIGNAL INJECTION MODE, NOT USING REAL DATA ***" << endl;

  domu = true;
  if(channel==1 || channel==3) domu = false;

  //histos_t th1qcd;
  //histos_t th1zjets;

  TString outfile = (domu?TString("mu"):TString("el"))+ 
                    (channel<2?TString("dijet"):TString("boosted")) +
                    (injectSignal? TString("_siginject"):TString(""));

  TFile* outputForLimit=0;
  if (saveDataCards_) {
    outputForLimit = TFile::Open(outfile+".root", "recreate");
  }

  // histos with different norm. scale factors
  th1data    = new histos_t("th1data",false);
  th1ww      = new histos_t("th1ww");
  th1wz      = new histos_t("th1wz");
  th1wv      = new histos_t("th1wv");
  th1wjets   = new histos_t("th1wjets");
  th1ttbar   = new histos_t("th1ttbar");
  //th1stops   = new histos_t("th1stops");
  //th1stopt   = new histos_t("th1stopt");
  //th1stopstbar  = new histos_t("th1stopstbar");
  //th1stopttbar  = new histos_t("th1stopttbar");
  th1stoptw    = new histos_t("th1stoptw");
  th1stoptbarw = new histos_t("th1stoptbarw");
  th1Top       = new histos_t("th1Top");

  th1tot     = new histos_t();

  TString cutsDijet("(W_pt<200.) && (dijetPt>70.) && (abs(JetPFCor_Eta[0])<2.4) && (abs(JetPFCor_Eta[1])<2.4) && (abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5) &&(abs(JetPFCor_dphiMET[0])>0.4) &&(JetPFCor_Pt[0]>40.) &&(JetPFCor_Pt[1]>35.) &&(JetPFCor_Pt[2]<30.) &&(JetPFCor_bDiscriminatorCSV[0]<0.244) &&(JetPFCor_bDiscriminatorCSV[1]<0.244) && (Mass2j_PFCor>70. && Mass2j_PFCor<100.)");


  // Do not put jet pt in the cut string here, since it is going to be smeared
  TString cutsMerged("(vbf_event==0) && (W_pt>200.) &&(abs(GroomedJet_CA8_eta[0])<2.4)&&(ggdboostedWevt==1) && (GroomedJet_CA8_deltaphi_METca8jet[0]>2.0) && (GroomedJet_CA8_deltaR_lca8jet[0]>1.57) && (numPFCorJetBTags<1) && (GroomedJet_CA8_tau2tau1[0]<0.55) && (GroomedJet_CA8_mass_pr[0]>70. && GroomedJet_CA8_mass_pr[0]<100.)");

  TString wtstr = "(puwt*effwt)";

  TString        lepton_cut = "(event_met_pfmet >30) && (W_mt>30.) && (W_electron_pt>35.)";
  if(channel==0) lepton_cut = "(event_met_pfmet >25) && (W_mt>30.) && (W_muon_pt>25.) && (abs(W_muon_eta)<2.1)";
  if(channel==1) lepton_cut = "(event_met_pfmet >30) && (W_mt>30.) && (W_electron_pt>30.)";
  if(channel==2) lepton_cut = "(event_met_pfmet >50) && (W_mt>30.) && (W_muon_pt>30.) && (abs(W_muon_eta)<2.1)";
  if(channel==3) lepton_cut = "(event_met_pfmet >70) && (W_mt>30.) && (W_electron_pt>35.)";

  TString And = " && ";
  TString Open = "(";
  TString Close = ")";

  TString jet_cut = cutsDijet;
  TString jetptcuts,mcjetptcuts,wwhimcjetptcuts,wzhimcjetptcuts;

  if(channel>1) jet_cut = cutsMerged;

  char* observable = "dijetPt";
  char* mcobservable = "dijetPt";
  char* xtitle = "p_{T}^{jj} [GeV]"; 
  double jetptcutval,wwjetptcutoff,wzjetptcutoff;
  if(channel>1) {
    double jetthresh = 80;
    /*observable = 
      "(GroomedJet_CA8_pt[0]>jetthresh)+\
      (GroomedJet_CA8_pt[1]>jetthresh)+\
      (GroomedJet_CA8_pt[2]>jetthresh)+\
      (GroomedJet_CA8_pt[3]>jetthresh)+\
      (GroomedJet_CA8_pt[4]>jetthresh)+\
      (GroomedJet_CA8_pt[5]>jetthresh)";*/

#if 0
    mcobservable = "W_pt";
    observable = "W_pt";
#else
    mcobservable = "GroomedJet_CA8_pt_smeared[0]";
    observable = "GroomedJet_CA8_pt[0]";
#endif

    char *jet1var= "GroomedJet_CA8_pt[0]";
    char *jet2var= "GroomedJet_CA8_pt[1]";
    char *mcjet1var="GroomedJet_CA8_pt_smeared[0]";
    char *mcjet2var="GroomedJet_CA8_pt_smeared[1]";

    jetptcutval = 200.;
    wwjetptcutoff = 1300.; // value where the WW hi pt MC starts to run out of statistics
    wzjetptcutoff = 1100.; // value where the WZ hi pt MC starts to run out of statistics
    jetptcuts   = Form("(%s > %f) && (%s < %f)",  jet1var,jetptcutval,jet2var,jetthresh);
    //mcjetptcut = Form("(%s > %f)",jet1var,jetptcutval);
    mcjetptcuts     = Form("((%s > %f) && (%s < %f))&&(%s < %f)",mcjet1var,jetptcutval,mcjet1var,wwjetptcutoff,mcjet2var,jetthresh);
    wwhimcjetptcuts = Form("((%s > %f) && (%s < %f))&&(%s < %f)",mcjet1var,jetptcutval,mcjet1var,wwjetptcutoff,mcjet2var,jetthresh);
    wzhimcjetptcuts = Form("((%s > %f) && (%s < %f))&&(%s < %f)",mcjet1var,jetptcutval,mcjet1var,wzjetptcutoff,mcjet2var,jetthresh);
    xtitle = "p_{T}^{j} [GeV]";
  }


  TCut mccut( Open+lepton_cut+And+jet_cut+And+mcjetptcuts+Close );

  TCut datacut;
  if (blinded)
    datacut = TCut( Open + lepton_cut + And + jet_cut + And + jetptcuts + And +
		    Form("(%s < 520.)", observable) + Close );
  else
    datacut = TCut( Open + lepton_cut + And + jet_cut + And + jetptcuts + Close );

  InstantiateTrees();

  if (!injectSignal) {
    th1data->fillFromTree(treedata,observable,datacut);
    cout << "th1data->hobs fSumw2.fN = " << th1data->hobs->GetSumw2N() << endl;
  }

  // ------- Get WW/WZ ------- 

  th1ww->fillFromTree(treewwlo,mcobservable,mccut,wtstr,WW_scale_lo,0,stitchwwgev);
  th1wz->fillFromTree(treewzlo,mcobservable,mccut,wtstr,WZ_scale_lo,0,stitchwzgev);

  // cut off integration of high tail when statistics become low and the parametrization is invalid: 
  TCut wwhimccut( Open + lepton_cut+And+jet_cut+And+wwhimcjetptcuts + Close );
  TCut wzhimccut( Open + lepton_cut+And+jet_cut+And+wzhimcjetptcuts + Close );

  th1ww->fillFromTree(treewwhi,mcobservable,wwhimccut,wtstr,WW_scale_hi,stitchwwgev);
  th1wz->fillFromTree(treewzhi,mcobservable,wzhimccut,wtstr,WZ_scale_hi,stitchwzgev);

  // ------- Get ttbar ------- 

  // for combining ttbar files
  TString mttstr = 
    TString("sqrt((W_top_E+W_atop_E)^2 -") +
    TString("     (W_top_px+W_atop_px)^2 -") +
    TString("     (W_top_py+W_atop_py)^2 -") +
    TString("     (W_top_pz+W_atop_pz)^2)");

  // weight default sample by 1 for mtt<700, by half above for adding the high mtt samples
  TString mttwt = Open+Open+mttstr+TString("<700)?1.0:0.5")+Close + TString("*") + wtstr;

  th1ttbar->fillFromTree(treettblo,mcobservable,mccut,mttwt,ttbar_scale);
  th1ttbar->fillFromTree(treettbmd,mcobservable,mccut,wtstr,0.5*mtt700_scale);
  th1ttbar->fillFromTree(treettbhi,mcobservable,mccut,wtstr,0.5*mtt1k_scale);

    // ------- Get WJets ------- 
  //th1wjets = new TH1D("th1wjets", "th1wjets", bins, dm_min, dm_max);

  th1wjets->fillFromTree(treewjlo,mcobservable,mccut,wtstr,WJ_scale_lo,0,stitchwjgev);
  th1wjets->fillFromTree(treewjhi,mcobservable,mccut,wtstr,WJ_scale_hi,stitchwjgev);

  // sherpa sample
  //th1wjets->fillFromTree("th1wjets",treewjsherpa,mcobservable,mccut,wtstr);

  // ------- Get QCD ------- 
  //treeqcd->Draw(TString(observable)+TString(">>th1qcd"), mccut, "goff");

  // ------- Get Z+Jets ------- 
  //treezj->Draw(TString(observable)+TString(">>th1zjets"), mccut, "goff");


  // ------- Get Single top ------- 
  
  //th1stops->fillFromTree  (treests,mcobservable,mccut,wtstr,STopS_scale);
  //th1stopt->fillFromTree  (treestt,mcobservable,mccut,wtstr,STopT_scale);
  //th1stopps->fillFromTree (tree64, mcobservable,mccut,wtstr,SToppS_scale);
  //th1stoppt->fillFromTree (tree65, mcobservable,mccut,wtstr,SToppT_scale);
  th1stoptw->fillFromTree   (treestoptw,mcobservable,mccut,wtstr,STopTW_scale);
  th1stoptbarw->fillFromTree(treestoptbarw, mcobservable,mccut,wtstr,STopTbarW_scale);

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
  SetupEmptyHistogram(xtitle);
  
  // ---- Sum all backgrounds ----------
  TH1D* th1wv_no_overflow = (TH1D *)th1wv->hobs->Clone("th1wv_no_overflow");
  SumAllBackgrounds();

  makeSignalHistos4Display(mcobservable, mccut, wtstr);

  // ---- Compose the stacks ----------
  THStack* hs = new THStack("hs","MC contribution");

  THStack* hscum = new THStack("hscum","MC contribution");

  //hs->Add(th1zjets); 
 //hs->Add(th1qcd);
  hs->Add(th1Top->hobs);         hscum->Add(th1Top->hcum);
  hs->Add(th1wjets->hobs);       hscum->Add(th1wjets->hcum);
  hs->Add(th1wv->hobs);          hscum->Add(th1wv->hcum); // add WW in with backgrounds for display only
  hs->Add(wwatgc4Display->hobs); hscum->Add(wwatgc4Display->hcum);
  hs->Add(wzatgc4Display->hobs); hscum->Add(wzatgc4Display->hcum);

#if 0
  if (injectSignal) {
    TH1D* th1splusb = ( TH1D*) th1tot->Clone("th1splusb");
    th1splusb->Add(th1wv);
    th1splusb->Add(wwatgc4Display);
    th1splusb->Add(wzatgc4Display);

    int scale=1;
    for (int i=0; i<gRandom->Poisson(scale*th1splusb->Integral()); ++i)
      th1data->Fill(th1splusb->GetRandom());

    th1data->Scale(1./scale);
  }
#endif

  makeRatioHisto();

  // ---- Stack for shape systematics Up ----------
  double bkgd_norm_fracerror = domu ? mu_bkgd_norm_error : el_bkgd_norm_error;
  double sig_norm_fracerror = domu ? mu_sig_norm_error : el_sig_norm_error;

  cout << "Background normalization fractional systematic = " << bkgd_norm_fracerror << endl;
  cout << "Signal     normalization fractional systematic = " << sig_norm_fracerror << endl;

  //TF1* formScaleUp = new TF1("formScaleUp", "1.0+0.4*log(x/5)", dm_min, dm_max);
  //TF1* formScaleDn = new TF1("formScaleDn", "1.0-0.2*log(x/5)", dm_min, dm_max);
  //TF1* formScaleUp = new TF1("formScaleUp", Form("1.0+%f",bkgd_norm_fracerror), dm_min, dm_max);
  //TF1* formScaleDn = new TF1("formScaleDn", Form("1.0-%f",bkgd_norm_fracerror), dm_min, dm_max);
  TF1* formScaleUp = new TF1("formScaleUp", "1.0-0.0002*x", dm_min, dm_max);
  TF1* formScaleDn = new TF1("formScaleDn", "1.0+0.0002*x", dm_min, dm_max);

  TH1D *systUp = (TH1D*) th1wjets->hobs->Clone("systUp");
  systUp->Multiply(formScaleUp);
  systUp->Scale(th1wjets->hobs->Integral()/systUp->Integral());
  //systUp->Add(th1zjets);
  //systUp->Add(th1qcd);
  systUp->Add(th1Top->hobs);
  systUp->SetFillColor(0);
  systUp->SetLineStyle(2);
  systUp->SetLineColor(2);
  systUp->SetLineWidth(3);

  // ---- Stack for shape systematics Down ----------
  TH1D *systDown = (TH1D*) th1wjets->hobs->Clone("systDown");
  systDown->Multiply(formScaleDn);
  systDown->Scale(th1wjets->hobs->Integral()/systDown->Integral());
  //systDown->Add(th1zjets);
  //systDown->Add(th1qcd);
  systDown->Add(th1Top->hobs);
  systDown->SetFillColor(0);
  systDown->SetLineWidth(3);
  systDown->SetLineStyle(2);
  systDown->SetLineColor(2);

  /////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  
  DrawItAll(channel, hs, th1tot->hobs, th1data->hobs, hhratio->hobs, outfile);
  DrawItAll(channel, hscum, th1tot->hcum, th1data->hcum, hhratio->hcum, outfile);

  ///// -------------------------------//////

  if(saveDataCards_) {
    outputForLimit->cd();
    th1data->hobs->SetName("data_obs");     th1data->hobs->Write("data_obs");
    th1tot->hobs->SetName("background");    th1tot->hobs->Write("background");
    th1ww->hobs->SetName("ww");             th1ww->hobs->Write("ww");
    th1ww->hobs->SetName("wz");             th1wz->hobs->Write("wz");
    th1wv->hobs->SetName("diboson");        th1wv->hobs->Write("diboson");

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

