
void cmsLabel(TCanvas *canvas,
	      double intlumi_mu,
	      double intlumi_el)
{
  TLatex latex;

  latex.SetNDC();

  canvas->cd();

  latex.SetTextSize(0.05);
  latex.SetTextFont(61); // helvetica bold
  latex.SetTextAlign(11); // align left
  latex.DrawLatex(0.18,0.93,"CMS");

  latex.SetTextSize(0.04);
  latex.SetTextFont(52); // helvetica italic
  latex.DrawLatex(0.29,0.93,"Preliminary");

  latex.SetTextAlign(31); // align right
  latex.SetTextSize(0.037);
  latex.SetTextFont(42); // helvetica italic

  if (intlumi_el > 0 && intlumi_mu > 0)
    latex.DrawLatex(0.94,0.93,Form("%4.1lf fb^{-1} (#mu)+%4.1lf fb^{-1} (e) (8 TeV)",intlumi_mu,intlumi_el));
  else if (intlumi_mu > 0)
    latex.DrawLatex(0.94,0.93,Form("%4.1lf fb^{-1} (#mu) (8 TeV)",intlumi_mu));
  else if (intlumi_el > 0)
    latex.DrawLatex(0.94,0.93,Form("%4.1lf fb^{-1} (e) (8 TeV)",intlumi_el));

  canvas->Update();
}

void atgcstyle()
{
  // import os
  // macroPath = gROOT.GetMacroPath()
  // macroPath += os.environ["CMSSW_BASE"] + "/src/ElectroWeakAnalysis/VPlusJets/test:"
  // gROOT.SetMacroPath(macroPath)
  // del os

  gROOT->SetStyle("Plain");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat("iouRMe");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(1112);
  gStyle->SetOptTitle(0);
  
  gStyle->SetCanvasDefH(600); // Height of canvas
  gStyle->SetCanvasDefW(600); // Width of canvas
  gStyle->SetErrorX(0.);
  
  gStyle->SetMarkerStyle(20);
  
  // For the fit/function:
  gStyle->SetFuncColor(2);
  gStyle->SetFuncStyle(1);
  gStyle->SetFuncWidth(1);
  
  //  Margins:
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.18); // was 0.16
  gStyle->SetPadRightMargin(0.06);// was 0.02
  
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.07, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.3); // was 1.25
  
  //  For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.06, "XYZ");
  gStyle->SetNdivisions(505, "XYZ");
  
  // if (gSystem->DynamicPathName("libFWCoreFWLite.so",True););:
  //     print "adding RooFit ...",
  //     scramCmd = ["scram","tool","info","roofitcore"]
  //     grepCmd = ["grep", "INCLUDE"]
  //     pscram = subprocess.Popen(scramCmd, stdout = subprocess.PIPE);
  //     pgrep = subprocess.Popen(grepCmd, stdin=pscram.stdout,
  //                              stdout=subprocess.PIPE);
  //     pscram.stdout.close();
  //     output = pgrep.communicate();[0]
  //     if (pgrep.returncode == 0);:
  //         roofitinc = output.split("=");[1].rstrip();
  //         //// print roofitinc
  //         gROOT.GetInterpreter();.AddIncludePath(roofitinc);
  //         roofitinc = "-I"" + roofitinc + """
  //         gSystem.AddIncludePath(roofitinc);
  //         print "done"
  //     else:
  //         print "failed"
  //         print "scram returned:",pscram.returncode,"grep:",pgrep.returncode
}

