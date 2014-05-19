#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TTree.h"
#include "TString.h"
#include "TStyle.h"
#include "TRandom.h"


TH1D *h_smearfactor;
TH1D *h_sigma_MC;


TRandom *rand;

//======================================================================

#define NUMETABINS 5
const double etabins[NUMETABINS+1] = {
  0.0,0.5,1.1,1.7,2.3,5.0
};

typedef struct {
  double etastart;
  double etastop;
  double datamcratio;
  double staterror;
  double plussyst;
  double minussyst;
} jertablentry_t;

const jertablentry_t jertable[] = {
  {0.0,	0.5, 	1.052, 0.012, 0.062, 0.061},
  {0.5,	1.1, 	1.057, 0.012, 0.056, 0.055},
  {1.1,	1.7, 	1.096, 0.017, 0.063, 0.062},
  {1.7,	2.3, 	1.134, 0.035, 0.087, 0.085},
  {2.3,	5.0, 	1.288, 0.127, 0.155, 0.153}, 
};

double smearJetPt(double jetptgev, double absjeteta)
{
  double ratio=1;
  double fractional_sigma_MC =
    sqrt(16/(jetptgev*jetptgev) + 0.04/jetptgev + 0.0049);
  //sqrt((4/pt)^2 + (0.2/sqrt(pt))^2 + (0.07)^2)

  h_sigma_MC->Fill(fractional_sigma_MC);

  for (int ietabin=0; ietabin<NUMETABINS+1; ietabin++)
    if (absjeteta<jertable[ietabin].etastop) break;

  if (ietabin < NUMETABINS+1)
    ratio = jertable[ietabin].datamcratio;
  else
    ratio = jertable[NUMETABINS].datamcratio;

  double smearfactor = 1;
  if (ratio > 1)
    smearfactor = rand->Gaus(1,fractional_sigma_MC*sqrt(ratio*ratio - 1));

  h_smearfactor->Fill(smearfactor);

  return (jetptgev*smearfactor);
}

//======================================================================

void mkSmearedJetPtTree(const char *indir,
			const char *outdir,
			const char *filename,
			const char *jetptstr = "GroomedJet_CA8_pt",
			const char *jetetastr = "GroomedJet_CA8_eta")
{

  rand = new TRandom();

  TString inpath  = Form("%s/%s",indir,filename);
  TString outpath = Form("%s/%s",outdir,filename);

  cout << inpath << " --> " << outpath << endl;

  TFile *infile = new TFile(inpath,"READ");  assert (infile);
  TTree *intree = (TTree*)infile->Get("WJet"); assert (intree);

  TFile *outfile = new TFile(outpath,"RECREATE");
  TTree *smearedtree = intree->CloneTree(0);

  h_smearfactor = new TH1D("h_smearfactor","h_smearfactor",100,0.5,1.5);
  h_sigma_MC    = new TH1D("h_sigma_MC","h_sigma_MC",300,0,0.15);

  Float_t    smeared_pts[6];
  Float_t          jetpt[6];
  Float_t         jeteta[6];

  TString newjetptname = TString(jetptstr)+TString("_smeared");
  smearedtree->Branch(newjetptname,
		      &smeared_pts,
		      Form("%s[6]/F",newjetptname.Data()));

  infile->cd();

  intree->SetBranchAddress(jetptstr, &jetpt);
  intree->SetBranchAddress(jetetastr, &jeteta);


  Long64_t nentries = intree->GetEntriesFast();
  for (Long64_t ient = 0; ient<nentries; ient++) {
    intree->GetEntry(ient);
    for (int i=0; i<6; i++)
      smeared_pts[i] = smearJetPt(jetpt[i],TMath::Abs(jeteta[i]));
    smearedtree->Fill();
    if (!(ient%1000)) { printf("%ld\r",ient); fflush(stdout); }
  }

  outfile->cd();
  smearedtree->Write();
  h_smearfactor->Write();
  h_sigma_MC->Write();
  infile->Close();

  delete smearedtree;
  delete infile;
  delete outfile;
  delete rand;
}                                                 //  mkSmearedJetPtTree

//======================================================================
#if 0
void mkSmearedJetPtTree() {

  const char *jetptname = "GroomedJet_CA8_pt";
  const char *jetetaname = "GroomedJet_CA8_eta";

  //mkSmearedJetPtTree("InData","/tmp","RD_el_WW_CMSSW532.root",           jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_WWToAnything_CMSSW532.root", jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_WZ_CMSSW532.root",           jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_WpJ_PT180_Madgraph_CMSSW532.root", jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_TTJets_poheg_CMSSW532.root", jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_Mtt700to1000_CMSSW532.root", jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_Mtt1000toinf_CMSSW532.root", jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_STopS_T_CMSSW532.root",      jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_STopT_T_CMSSW532.root",      jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_STopTW_T_CMSSW532.root",     jetptname, jetetaname);
  mkSmearedJetPtTree("InData","/tmp",   "RD_el_STopS_Tbar_CMSSW532.root",   jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_STopT_Tbar_CMSSW532.root",   jetptname, jetetaname);
  //mkSmearedJetPtTree("InData",".",   "RD_el_STopTW_Tbar_CMSSW532.root",  jetptname, jetetaname);

  mkSmearedJetPtTree("InData","/tmp","RD_mu_WW_CMSSW532.root",           jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_WWToAnything_CMSSW532.root", jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_WZ_CMSSW532.root",           jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_WpJ_PT180_Madgraph_CMSSW532.root", jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_TTJets_poheg_CMSSW532.root", jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_Mtt700to1000_CMSSW532.root", jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_Mtt1000toinf_CMSSW532.root", jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopS_T_CMSSW532.root",      jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopT_T_CMSSW532.root",      jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopTW_T_CMSSW532.root",     jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopS_Tbar_CMSSW532.root",   jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopT_Tbar_CMSSW532.root",   jetptname, jetetaname);
  mkSmearedJetPtTree("InData",".",   "RD_mu_STopTW_Tbar_CMSSW532.root",  jetptname, jetetaname);
}

#endif
