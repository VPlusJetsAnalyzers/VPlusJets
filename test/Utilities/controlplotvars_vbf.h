#ifndef _VBFCONTROLPLOTVARS_H
#define _VBFCONTROLPLOTVARS_H

#include "plotvar_t.h"

////////////////////////////////////////////VBF
const plotVar_t vbfplotvars[] = {
  //    plotvar MINRange  MAXRange  NBINS  slog xlabel    outfile  AMINRange  AMAXRange ANBINS mva_in hplot drawleg
#if 0
  { "vbf_wbj_pt/vbf_wjj_m",
                        0.3, 0.7,  20, 2, "Jet 2 p_{T}/m_{jj}", "vbfjet2pt_ov_mjj", 0.25, 0.75, 25, 0, 0, 0 },
  
  { "vbf_waj_pt",       40,   300, 27, 3, "Hadronic W Jet1 p_{T}", "vbfwjeta_pt",      20,  300, 28, 0, 0, 1 },
  { "vbf_wbj_pt",       30,   200, 17, 3, "Hadronic W Jet2 p_{T}", "vbfwjetb_pt",      20,  200, 18, 0, 0, 0 },
  { "vbf_waj_eta",     -2.5 , 2.5, 25, 1, "Hadronic W Jet1 #eta",  "vbfwjeta_eta",   -2.7,  2.7, 27, 0, 0, 0 },
  { "vbf_wbj_eta",     -2.5 , 2.5, 25, 1, "Hadronic W Jet2 #eta",  "vbfwjetb_eta",   -2.7,  2.7, 27, 0, 0, 0 },
#endif
  { "vbf_aj_pt",        40,   300, 27, 3, "Tagjet1  p_{T}",        "vbftagjeta_pt",    20,  300, 28, 0, 0, 1 },
#if 0
  { "vbf_bj_pt",        30,   200, 17, 3, "Tagjet2  p_{T}",        "vbftagjetb_pt",    20,  200, 18, 0, 0, 0 },
  { "vbf_aj_eta",      -4.5 , 4.5, 45, 1, "Tagjet1 #eta",          "vbftagjeta_eta", -5.0,  5.0, 50, 0, 0, 0 },
  { "vbf_bj_eta",      -4.5 , 4.5, 45, 1, "Tagjet2 #eta",          "vbftagjetb_eta", -5.0,  5.0, 50, 0, 0, 0 },
  { "vbf_wjj_m",        40,   400, 36, 3, "m_{jj} (GeV)",          "vbfwjjmass",       30,  400, 37, 0, 0, 0 },
  { "W_mt",             50,   170, 24, 3, "W Transverse Mass (GeV)", "vbf_W_mt",       40,  170, 26, 0, 0, 1 },
  { "event_met_pfmet",  25,   205, 18, 3, "pf MET (GeV)",     "event_met_pfmet"  ,     15, 205, 19,  0, 0, 0 },

  { "vbf_lvjj_m",        140, 700, 56, 3, "m_{l#nujj} (GeV)",  "vbfmlvjj",        110, 700, 59, 0, 0, 1 },

  // MVA training variables:
  { "vbf_lvjj_pt",         0, 150, 15, 3, "p_{T} of WW (GeV)", "vbfptlvjj",         0, 150, 15, 1, 0, 0 },
  { "vbf_lvjj_y",       -2.5, 2.5, 25, 1, "WW rapidity" ,      "vbfetalvjj"   ,  -3.0, 3.0, 30, 1, 0, 0 },
  { "vbf_wjj_ang_ha",      0,   1, 10, 1, "Cos(#theta_{1})" ,  "vbfangha",          0,   1, 10, 1, 0, 1 },
  { "vbf_wjj_ang_hb",      0,   1, 10, 1, "Cos(#theta_{2})" ,  "vbfanghb",          0,   1, 10, 1, 0, 0 },
  { "vbf_wjj_ang_hs",   -1.0,   1, 20, 1, "Cos(#theta*)" ,     "vbfanghs",       -1.0,   1, 20, 1, 0, 0 },
  { "vbf_wjj_ang_phi",  -3.2, 3.2, 32, 6, "#Phi (rad)" ,       "vbfangphi",      -3.5, 3.5, 35, 1, 0, 0 },
  { "vbf_wjj_ang_phib", -3.2, 3.2, 32, 6, "#Phi_{1} (rad)",    "vbfangphib",     -3.5, 3.5, 35, 1, 0, 0 },
  { "W_electron_charge",-1.2, 1.2, 24, 1, "Electron Charge" ,  "vbfWelcharge",   -1.2, 1.2, 24, 1, 0, 0 },
  { "W_muon_charge",    -1.2, 1.2, 24, 1, "Muon Charge" ,      "vbfWmucharge",   -1.2, 1.2, 24, 1, 0, 0 },
  { "vbf_jj_deta",         3, 9.0, 30, 1, "Tagjet  #Delta#eta",    "vbftagjet_deta",      3,  9.0, 30, 1, 0, 0 },
  { "vbf_jj_m",         300, 1200, 18, 3, "Tagjet Invariant Mass (GeV)", "vbftagjet_mass", 200,  1200, 20, 1, 0, 0 },

  { "mvavbf500mu",       -0.2, 1.2, 28, 2, "Likelihood Output ",       "vbfmva500",  -0.2, 1.2, 28, 0, 0, 1 },
  { "mvavbf350mu",          0,   1, 25, 2, "Likelihood Output mu350 ", "vbfmva350",     0,   1, 10, 0, 0, 1 },
  { "mvavbf180el",          0,   1, 25, 2, "Likelihood Output el180 ", "vbfmva180",     0,   1, 10, 0, 0, 1 },
#endif
  { "", 0.0,0.0,0,0,"","",0.,0.,0,0,0,0 }
};

#endif // _VBFCONTROLPLOTVARS_H
