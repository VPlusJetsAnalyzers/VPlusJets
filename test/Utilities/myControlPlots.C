#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <assert.h>
#include "TROOT.h"
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
#include "TCut.h"
#include "TSystem.h"

#include "tdrstyle.C"
#include "utils.C" // Tokenize

typedef struct {
  int     index;
  TString samplename;
  TString treefilename;
  double xsecpblumi;
  double otherscale;
  int    nMCevents;
  int    colorcode;
  int    stackit;
}
SampleInfo_t;

#include "controlplotvars_boosted.h"
#include "controlplotvars_common.h"
#include "controlplotvars_higgs.h"
#include "controlplotvars_vbf.h"

using namespace std;

double intLUMIinvpb;

//=====================================================================================
// SYNOPSIS:
//   1. Prepare "InData" and "OutDir" directories; e.g., "ln -s . OutDir" to go to current dir
//   2. Prepare "cuttable.txt" of cut names and cut strings
//   3. root [0] .L myControlPlots.C+
//      root [1] myControlPlots("cuttable.txt")
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

//======================================================================

void loadCutString(const char *filename, TString& cutstring)
{
  FILE *fp = fopen(filename,"r");
  if (!fp) {
    cout << "Error, file " << TString(filename) << " not found." << endl;
    exit(-1);
  }

  char line[512];

  for (int i=0; !feof(fp) && fgets(line,512,fp); i++) {
    if (!strlen(line) || line[0]=='#') continue; // comments are welcome

    if (cutstring.Length()) cutstring += " && ";

    string strline(line);
    strline.erase(strline.size()-1,1);     // i.e., "pop_back" in c++11. shed the \n
    vector<string> fields;

    // expect columns with fields cutname, cutvalue, possible embedded spaces both
    // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
    // instead.
    //
    Tokenize(strline,fields, "\"");

    //for (size_t j=0; j<fields.size(); j++)
    //cout << j << ": \"" << fields[j] << "\"" << endl;

    assert (fields.size()==3);
    cutstring += TString(fields.at(2));
  }
}                                                       // loadCutString

//======================================================================

class Sample {
public:
  Sample(const SampleInfo_t& sinfo) {
    info_ = sinfo;
    //cout << "sample = " << name_ << endl;
    TFile *f = new TFile (sinfo.treefilename, "READ"); if (!f) { cerr << "Couldn't find file " << sinfo.treefilename << endl; exit(-1); }
    tree_ =  (TTree *)f->Get("WJet"); if (!tree_) { cerr << "Couldn't find tree WJet in file " << sinfo.treefilename << endl; exit(-1); }
  }
  ~Sample() { delete tree_; }
  TTree *Tree() const { return tree_; }
  TString name() const { return info_.samplename; }
  TString filename() const { return info_.treefilename; }
  bool stackit() const { return (info_.stackit != 0); }
  int colorcode() const { return info_.colorcode; }
  double otherscale() const { return info_.otherscale; }
  TH1 *Draw(const plotVar_t& pv, const TCut& cut, const TCut& cutSQ ) {
    cout << "\tDrawing " << pv.plotvar << " for sample = " << info_.samplename << " ... ";
    TString hname = TString("th1")+ pv.outfile + Form("%d",info_.index);
    TH1 *histo = new TH1D(hname, hname, pv.NBINS, pv.MINRange, pv.MAXRange);
    assert(histo);
    histo->Sumw2();
    assert(tree_);
    cout << tree_->Draw(pv.plotvar+TString(">>")+hname, cut, "goff") << " events" << endl;

    if (strlen((const char *)cutSQ)) {
      hname = TString("th1") + pv.outfile + Form("%d",info_.index) + TString("SQ");
      TH1 *histoSQ = new TH1D(hname, hname, pv.NBINS, pv.MINRange, pv.MAXRange);
      tree_->Draw(pv.plotvar+TString(">>")+hname, cutSQ, "goff");
      for(int hi=1;hi<=pv.NBINS;hi++) {
	histo->SetBinError(hi,sqrt(histoSQ->GetBinContent(hi)));
      }
      delete histoSQ;
    }
    if (info_.nMCevents)
      histo->Scale(info_.xsecpblumi*info_.otherscale/(double)info_.nMCevents);

    return histo;
  }
private:
  SampleInfo_t info_;
  TTree *tree_;
};

//======================================================================

void loadSamples(const char *filename,vector<Sample *>& samples)
{
  FILE *fp = fopen(filename,"r");
  if (!fp) {
    cout << "Error, file " << TString(filename) << " not found." << endl;
    exit(-1);
  }

  char line[512];

  intLUMIinvpb=-1; // obvious error condition

  for (int i=0; !feof(fp) && fgets(line,512,fp); i++) {
    if (!strlen(line) || line[0]=='#') continue; // comments are welcome

    string strline(line);
    strline.erase(strline.size()-1,1);     // i.e., "pop_back" in c++11. shed the \n
    vector<string> fields;

    // expect columns with fields cutname, cutvalue, possible embedded spaces both
    // within and between, so " " or "\t" cannot be used as delimiters. Require quotes
    // instead.
    //
    Tokenize(strline,fields, " \t");

    //for (size_t j=0; j<fields.size(); j++)
    //cout << j << ": \"" << fields[j] << "\"" << endl;

    assert (fields.size()==7);

    SampleInfo_t s;
    s.index        = i;
    s.samplename   = fields[0];
    s.treefilename = fields[1];
    s.xsecpblumi   = str2dbl(fields[2]);
    s.otherscale   = str2dbl(fields[3]);
    s.nMCevents    = str2int(fields[4]);
    s.colorcode    = str2int(fields[5]);
    s.stackit      = str2int(fields[6]);
    
    cout << "Loading sample " << s.samplename << " -> " << s.treefilename << endl;

    if (!samples.size()) {
      if (s.samplename.EqualTo("data")) {
	intLUMIinvpb = s.xsecpblumi;
	s.xsecpblumi = 1;
	cout << "intLUMI = " << intLUMIinvpb << " pb^-1" << endl;
      } else {
	cerr << "First sample in the table must be 'data'" << endl;
	exit(-1);
      }
    } else if (!s.treefilename.Contains("QCD")) {
      s.otherscale *= intLUMIinvpb;
    }

    samples.push_back(new Sample(s) );
  }
}                                                         // loadSamples

//======================================================================

void myControlPlots(const char *cuttablefilename,
		    const char *qcdcuttablefilename,
		    const char *samplefilename,
		    const plotVar_t plotvars[])
{
  //gROOT->ProcessLine(".L tdrstyle.C");

  TString unwtcutstring, qcdcutstring;

  loadCutString(cuttablefilename, unwtcutstring);
  if (strlen(qcdcuttablefilename))
    loadCutString(qcdcuttablefilename, qcdcutstring);

  //  const char* the_cut = "1";
  //  double BINWIDTH = ((MAXRange-MINRange)/NBINS);

 // Get the input trees:

  vector<Sample *> samples;

  loadSamples(samplefilename,samples);

 // Data

  Sample *sdata = samples[0];
  cout << "ndata =" << sdata->Tree()->GetEntries() <<endl;

  TFile f("plotvar_histo.root", "RECREATE");

  //============================================================
  //  VARIABLE LOOP
  //============================================================

  for (int ivar=0; ; ivar++) {

    plotVar_t pv = plotvars[ivar];
 
    if ( !pv.plotvar.Length() ) break;

    cout << pv.plotvar << "\t"<<pv.MINRange<<"\t" << pv.MAXRange<<"\t" << pv.NBINS<<"\tTHE CUT " << endl;

    if ( sdata->Tree()->Draw(pv.plotvar,"","goff",1) == -1 ) { // check if the variable exists in the tree
      cout << "\t...can't be plotted!" << endl;
      continue;
    }

    TCut the_cut(TString("effwt*puwt*")+unwtcutstring);
    TCut the_cutE(TString("effwt*puwt*puwt*")+unwtcutstring);
    TCut qcd_cut;
    if (qcdcutstring.Length())
      qcd_cut = TCut(TString("effwt*puwt*")+qcdcutstring);

    TCut nullcut("");

    const double BINWIDTH = ((pv.MAXRange-pv.MINRange)/pv.NBINS);

    map<TString, TH1 *> m_histos;
    map<TString, bool> m_stacked;

    double totevents = 0.;
    TH1 * th1qcd = NULL;
    double qcdfrac = 0.;

    //============================================================
    // DRAW THE VARIABLE FOR ALL SAMPLES, CREATE HISTOS
    //============================================================

    for (size_t isamp=0; isamp<samples.size(); isamp++) {
      Sample *s = samples[isamp];

      m_stacked[s->name()] = false;

      TH1 *h;

      if (s->name().EqualTo("data"))            h = s->Draw(pv, the_cut, nullcut); // effwt*puwt==1 for data!
      else if (s->filename().Contains("QCD")) { h = s->Draw(pv, qcd_cut, nullcut);
	th1qcd = h;
	qcdfrac = s->otherscale();
      } else {                                  h = s->Draw(pv, the_cut, the_cutE);
	if (s->stackit()) {
	  totevents += h->Integral();
	}
      }

      map<TString, TH1 *>::iterator mit = m_histos.find(s->name());
      if (mit == m_histos.end()) {
	if (s->stackit()) {
	  h->SetFillColor(s->colorcode());
	  h->SetLineColor(s->colorcode());
	  h->SetLineWidth(0);
	}
	m_histos[s->name()] = h;
      } else {
	mit->second->Add(h);
      }
    }

    //============================================================
    // COUNT EVENTS, RENORM TO DATA, CONSTRUCT THE TSTACK & LEGEND
    //============================================================

    TH1 *th1data = m_histos["data"];

    assert(th1data);

    double ndata = th1data->Integral();

    if (th1qcd) {
      // QCD = qcdfrac * data, QCD + Sum(MC) = data, ergo...
      //th1qcd->Scale(qcdfrac*totevents/((1-qcdfrac)*th1qcd->Integral()));

      th1qcd->Scale(qcdfrac*ndata/th1qcd->Integral()); // matches previous script

      totevents += th1qcd->Integral();
    }

    double renorm = ndata/totevents;

    cout << "den = " << totevents << endl;
    cout << "data = " << ndata  <<endl;
    cout  << "data/den = " << renorm << endl;


    // Setup the stack and total
    THStack* hs = new THStack("hs","MC contribution");
    TH1D *th1tot = new TH1D("th1tot", "th1tot", pv.NBINS, pv.MINRange, pv.MAXRange);

    // Set up the legend

    float  legX0=0.65, legX1=0.99, legY0=0.4, legY1=0.88;
    // float  legX0=0.35, legX1=0.85, legY0=0.4, legY1=0.88;
    // float  legX0=0.18, legX1=0.52, legY0=0.4, legY1=0.88;
    TLegend * Leg = new TLegend( legX0, legY0, legX1, legY1);
    Leg->SetFillColor(0);
    Leg->SetFillStyle(0);
    Leg->SetTextSize(0.04);

    if (TString(cuttablefilename).Contains("Mu"))
      Leg->AddEntry(th1data,  "Muon Data",  "PLE");
    else
      Leg->AddEntry(th1data,  "Electron Data",  "PLE");

    vector<double> binErrSQ(pv.NBINS,0.);

    vector<pair<TString, TH1 *> > v_legentries;

    for (size_t isamp=1; isamp<samples.size(); isamp++) {
      Sample *s = samples[isamp];
      if (m_stacked[s->name()]) continue;

      map<TString, TH1 *>::iterator mit = m_histos.find(s->name());
      TH1 *h = mit->second;
      h->Scale(renorm);

      cout << s->name() << " = " << h->Integral() << endl;

      if(s->stackit()) {
	hs->Add(h);
	th1tot->Add(h);
	m_stacked[s->name()] = true;
	v_legentries.push_back(*mit);
	for (int ibin=1; ibin <= pv.NBINS; ibin++)
	  binErrSQ[ibin-1] += h->GetBinError(ibin)*h->GetBinError(ibin);
      }
    }

    // Reverse the order for the legend
    for (vector<pair<TString, TH1 *> >::reverse_iterator
	   rit = v_legentries.rbegin();
	 rit != v_legentries.rend();
	 rit++)
      Leg->AddEntry(rit->second, rit->first, "F");

    TH1D* th1totClone = ( TH1D*) th1tot->Clone("th1totClone");
    th1totClone->SetMarkerStyle(0);
    th1totClone->SetFillStyle(3003);
    th1totClone->SetFillColor(11);
    th1totClone->SetLineColor(0);

    for(int ibin=1; ibin<=th1totClone->GetNbinsX(); ++ibin) {
      th1totClone->SetBinError(ibin, sqrt(binErrSQ[ibin-1]));
    }

    //============================================================
    // SETUP THE CANVAS
    //============================================================

//    gROOT->ProcessLine(".L tdrstyle.C");
    setTDRStyle();
    tdrStyle->SetErrorX(0.5);
    tdrStyle->SetPadRightMargin(0.05);

    tdrStyle->SetLegendBorderSize(0);

    TCanvas* c1 = new TCanvas(pv.plotvar,pv.plotvar,10,10, 800, 800);
    TPad *d1, *d2;

    c1->Divide(1,2,0,0);
    d1 = (TPad*)c1->GetPad(1);
    d1->SetPad(0.01,0.30,0.95,0.99);
    d2 = (TPad*)c1->GetPad(2);
    d2->SetPad(0.01,0.02,0.95,0.30);

    // Compose the stack

    d1->cd();
    gPad->SetBottomMargin(0.0);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.14);

    Leg->AddEntry(th1tot,  "MC Uncertainty",  "f");

    Leg->SetFillColor(0);

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

    th1tot->SetMinimum(0.01);
    int maxbin = th1data->GetMaximumBin();
    float maxval = th1data->GetBinContent(maxbin);
    cout << "maxval " <<maxval <<endl;
//    th1totempty->SetMaximum(2.5*maxval);
    th1totempty->SetMaximum(1.6*maxval);
    th1totempty->SetMinimum(0.01);
    if(pv.slog==1) th1totempty->SetMaximum(1.6*maxval);
    th1data->SetMinimum(0.01);

    // Draw it all

    th1totempty->Draw();
    //th1tot->Draw("e2same");
    th1data->Draw("esame");
    hs->Draw("samehist");

    th1tot->Draw("e2same");

    th1data->Draw("esame");
    cmspre(intLUMIinvpb/1000.0);    
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
    hhratio->GetYaxis()->SetTitleSize(0.1);
    hhratio->GetYaxis()->SetTitleOffset(0.5);
    hhratio->GetYaxis()->CenterTitle(true);
    hhratio->GetYaxis()->SetLabelSize(0.1);
    cout << hhratio->GetNbinsX() << endl;
    cout << th1tot->GetNbinsX() << endl;
    hhratio->Divide(th1tot);
    double binError(0.0), mcbinentry(0.0), mcerror(0.0);
    for(int i=0; i<hhratio->GetNbinsX(); ++i) {
      binError = hhratio->GetBinError(i);
      mcerror = th1tot->GetBinError(i);
      mcbinentry = th1tot->GetBinContent(i);
      if(mcbinentry>0.) mcerror /= mcbinentry;
      else mcerror = 0.0;
      binError = sqrt(binError*binError + mcerror*mcerror);
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

    TBox *errbox = new TBox(pv.AMINRange,0.974,pv.AMAXRange,1.026); // lumi systematic uncertainty
    errbox->SetFillColor(kGray);
    errbox->Draw();

#if 0
 	TF1 *f1 = new TF1("f1", "pol1",  pv.AMINRange, pv.AMAXRange);
	//f1->SetParameters(1.0,0.0);
       // f1->SetParameters(1,0.0);
	f1->FixParameter(0,1);
	f1->FixParameter(1,0);
	//cout<<" par1   "f1->GetParameter(0)<<endl;
	//cout<<" par2   "f1->GetParameter(1)<<endl;

        TFitResultPtr r =  hhratio->Fit("f1", "RBS");
	//TFitResultPtr r = hhratio->Fit(myFunc,"S");
       	    r->Print("V");     // print full information of fit including covariance matrix
#endif

    hhratio->Draw("esame");
    TLine *line; line = new TLine(pv.AMINRange,1.0,pv.AMAXRange,1.0);
    line->SetLineStyle(1);
    line->SetLineWidth(1);
    line->SetLineColor(1);
    line->Draw();

    TString outfile = TString("OutDir/")+TString(gSystem->BaseName(cuttablefilename)).ReplaceAll(".txt","")+TString("_")+pv.outfile;

    c1->Print(outfile+".png");
    c1->Print(outfile+".C");
    //gPad->WaitPrimitive();
    c1->Modified();
    c1->Update();
    c1->SaveAs(outfile+".pdf"); 

  } // var loop

  f.Write();

}                                                                // myControlPlots

//================================================================================

void dibresNobtagElplots()
{
  myControlPlots("DibosonResolvedElCutsNoBtag.txt",
		 "QCDcuts8TeV.txt",
		 "DibosonResolvedElAbtagSamples8TeV.txt",
		 commonplotvars);
}

//================================================================================


void dibresabtagElplots()
{
  myControlPlots("DibosonResolvedElCutsAntibtag.txt",
		 "QCDcuts8TeV.txt",
		 "DibosonResolvedElAbtagSamples8TeV.txt",
		 commonplotvars);
}

void dibresabtagMuplots()
{
  myControlPlots("DibosonResolvedMuCutsAntibtag.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars);
}

//================================================================================


void dibresbtagElplots()
{
  myControlPlots("DibosonResolvedElCutsBtag.txt",
		 "",
		 "DibosonResolvedElBtagSamples8TeV.txt",
		 commonplotvars);
}

void dibresbtagMuplots()
{
  myControlPlots("DibosonResolvedMuCutsBtag.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars);
}

//================================================================================


void dibbooElplots()
{
  myControlPlots("DibosonBoostedElCuts8TeV.txt",
		 "",
		 "DibosonBoostedElSamples8TeV.txt",
		 boostedplotvars);
}

void dibbooMuplots()
{
  myControlPlots("DibosonBoostedMuCuts8TeV.txt",
		 "",
		 "DibosonBoostedMuSamples8TeV.txt",
		 boostedplotvars);
}

//================================================================================

void dibresabtagElplotsByPU()
{
  gSystem->Exec("cat DibosonResolvedElCutsAntibtag.txt loPUcuts.txt >>/tmp/DibosonResolvedElCutsAntibtag_loPU.txt");
  gSystem->Exec("cat QCDcuts8TeV.txt loPUcuts.txt >>/tmp/QCDcuts8TeV_loPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsAntibtag_loPU.txt",
		 "/tmp/QCDcuts8TeV_loPU.txt",
		 "DibosonResolvedElAbtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedElCutsAntibtag.txt medPUcuts.txt >>/tmp/DibosonResolvedElCutsAntibtag_medPU.txt");
  gSystem->Exec("cat QCDcuts8TeV.txt medPUcuts.txt >>/tmp/QCDcuts8TeV_medPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsAntibtag_medPU.txt",
		 "/tmp/QCDcuts8TeV_medPU.txt",
		 "DibosonResolvedElAbtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedElCutsAntibtag.txt hiPUcuts.txt >>/tmp/DibosonResolvedElCutsAntibtag_hiPU.txt");
  gSystem->Exec("cat QCDcuts8TeV.txt hiPUcuts.txt >>/tmp/QCDcuts8TeV_hiPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsAntibtag_hiPU.txt",
		 "/tmp/QCDcuts8TeV_hiPU.txt",
		 "DibosonResolvedElAbtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj
}

//================================================================================

void dibresabtagMuplotsByPU()
{
  gSystem->Exec("cat DibosonResolvedMuCutsAntibtag.txt loPUcuts.txt >>/tmp/DibosonResolvedMuCutsAntibtag_loPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsAntibtag_loPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedMuCutsAntibtag.txt medPUcuts.txt >>/tmp/DibosonResolvedMuCutsAntibtag_medPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsAntibtag_medPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedMuCutsAntibtag.txt hiPUcuts.txt >>/tmp/DibosonResolvedMuCutsAntibtag_hiPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsAntibtag_hiPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj
}

//================================================================================

void dibresbtagElplotsByPU()
{
  gSystem->Exec("cat DibosonResolvedElCutsBtag.txt loPUcuts.txt >>/tmp/DibosonResolvedElCutsBtag_loPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsBtag_loPU.txt",
		 "",
		 "DibosonResolvedElBtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedElCutsBtag.txt medPUcuts.txt >>/tmp/DibosonResolvedElCutsBtag_medPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsBtag_medPU.txt",
		 "",
		 "DibosonResolvedElBtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedElCutsBtag.txt hiPUcuts.txt >>/tmp/DibosonResolvedElCutsBtag_hiPU.txt");
  myControlPlots("/tmp/DibosonResolvedElCutsBtag_hiPU.txt",
		 "",
		 "DibosonResolvedElBtagSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj
}

//================================================================================

void dibresbtagMuplotsByPU()
{
  gSystem->Exec("cat DibosonResolvedMuCutsBtag.txt loPUcuts.txt >>/tmp/DibosonResolvedMuCutsBtag_loPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsBtag_loPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedMuCutsBtag.txt medPUcuts.txt >>/tmp/DibosonResolvedMuCutsBtag_medPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsBtag_medPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj

  gSystem->Exec("cat DibosonResolvedMuCutsBtag.txt hiPUcuts.txt >>/tmp/DibosonResolvedMuCutsBtag_hiPU.txt");
  myControlPlots("/tmp/DibosonResolvedMuCutsBtag_hiPU.txt",
		 "",
		 "DibosonResolvedMuSamples8TeV.txt",
		 commonplotvars); // cut the variables down to jet0 pt and mjj
}

//================================================================================

void dibbooElplotsByPU()
{
  gSystem->Exec("cat DibosonBoostedElCuts8TeV.txt loPUcuts.txt >>/tmp/DibosonBoostedElCuts8TeV_loPU.txt");
  myControlPlots("/tmp/DibosonBoostedElCuts8TeV_loPU.txt",
		 "",
		 "DibosonBoostedElSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ

  gSystem->Exec("cat DibosonBoostedElCuts8TeV.txt medPUcuts.txt >>/tmp/DibosonBoostedElCuts8TeV_medPU.txt");
  myControlPlots("/tmp/DibosonBoostedElCuts8TeV_medPU.txt",
		 "",
		 "DibosonBoostedElSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ

  gSystem->Exec("cat DibosonBoostedElCuts8TeV.txt hiPUcuts.txt >>/tmp/DibosonBoostedElCuts8TeV_hiPU.txt");
  myControlPlots("/tmp/DibosonBoostedElCuts8TeV_hiPU.txt",
		 "",
		 "DibosonBoostedElSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ
}

//================================================================================

void dibbooMuplotsByPU()
{
  gSystem->Exec("cat DibosonBoostedMuCuts8TeV.txt loPUcuts.txt >>/tmp/DibosonBoostedMuCuts8TeV_loPU.txt");
  myControlPlots("/tmp/DibosonBoostedMuCuts8TeV_loPU.txt",
		 "",
		 "DibosonBoostedMuSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ

  gSystem->Exec("cat DibosonBoostedMuCuts8TeV.txt medPUcuts.txt >>/tmp/DibosonBoostedMuCuts8TeV_medPU.txt");
  myControlPlots("/tmp/DibosonBoostedMuCuts8TeV_medPU.txt",
		 "",
		 "DibosonBoostedMuSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ

  gSystem->Exec("cat DibosonBoostedMuCuts8TeV.txt hiPUcuts.txt >>/tmp/DibosonBoostedMuCuts8TeV_hiPU.txt");
  myControlPlots("/tmp/DibosonBoostedMuCuts8TeV_hiPU.txt",
		 "",
		 "DibosonBoostedMuSamples8TeV.txt",
		 boostedplotvars); // cut the variables down to pt(J) and mJ
}

//================================================================================
