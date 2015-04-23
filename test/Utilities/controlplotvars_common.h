#ifndef _COMMONCONTROLPLOTVARS_H
#define _COMMONCONTROLPLOTVARS_H

#include "plotvar_t.h"

const plotVar_t commonplotvars[] = {

//    plotvar	MINRange  MAXRange  NBINS  slog xlabel outfile AMINRange  AMAXRange ANBINS mva_in hplot drawleg

//    { "JetPFCor_Pt[1]/Mass2j_PFCor",
//                          0.3, 0.7,  20, 2, "Jet 2 p_{T}/m_{jj}", "jet2pt_ov_mjj", 0.25, 0.75, 25, 0, 0, 0 },
  { "JetPFCor_Pt[0]",   30,   200, 17, 3, "Leading Jet  p_{T}", "jetld_pt",      20,  200, 18, 0, 0, 1 },
  { "JetPFCor_Pt[1]",   30,   150, 12, 3, "Jet 2 p_{T}",        "jetnt_pt",      20,  150, 13, 0, 0, 0 },
  //{ "JetPFCor_Pt[2]",   30,   150, 12, 3, "Jet 3 p_{T}",        "jet3_pt",       20,  150, 13, 0, 0, 0 },
  { "JetPFCor_Eta[0]", -2.4 , 2.4, 12, 1, "Leading Jet  #eta",  "jetld_eta",   -2.6 , 2.6, 13, 0, 0, 0 },
  { "JetPFCor_Eta[1]", -2.4 , 2.4, 12, 1, "Jet 2 #eta",         "jetnt_eta",   -2.6 , 2.6, 13, 0, 0, 0 },
  //{ "JetPFCor_Eta[2]", -2.4 , 2.4, 12, 1, "Jet 3 #eta",         "jet3_eta",    -2.6 , 2.6, 13, 0, 0, 0 },
  { "Mass2j_PFCor",     40, 400, 18, 3,   "m_{jj} (GeV)",            "mjj",       20,  400, 19, 0, 0, 1 },
  { "W_mt",             20, 140, 12, 3,   "W Transverse Mass (GeV)", "W_mt",      20,  140, 12, 0, 0, 1 },

  { "event_met_pfmet",  25, 155, 13, 3,   "pf MET (GeV)",  "event_met_pfmet",     15,  155, 14, 0, 0, 0 },

  //{ "W_mtMVA",             20, 140, 12, 3,   "W Transverse Mass (GeV) withMVA MET", "W_mtMVA",      20,  140, 12, 0, 1, 1 },
  //{ "abs(event_metMVA_metPhi-W_muon_phi)",  0, 6.28, 20, 3,   " #Delta#phi(lepton,MVA MET) ",  "event_metMVA_lep_dphi",     0, 6.28, 20, 0, 0, 0 },
  //{ "abs(event_met_pfmetPhi-W_muon_phi)",  0, 6.28, 20, 3,   " #Delta#phi(lepton,pf MET) ",  "event_pfmet_lep_dphi",     0, 6.28, 20, 0, 0, 0 },

  //{ "event_met_pfsumet",  0, 2000, 40, 3,   "pf SET (GeV)",  "event_met_pfsumet",     0, 2000, 40, 0, 0, 0 },
  //{ "event_met_pfmetPhi",  -3.14, 3.14, 26, 3,   "pf MET #phi",  "event_met_pfmetPhi",     -3.14, 3.14, 26, 0, 0, 0 },
  //{ "event_metMVA_metPhi",  -3.14, 3.14, 26, 3,   "MVA MET #phi ",  "event_metMVA_metPhi",     -3.14, 3.14, 26, 0, 0, 0 },
  //{ "event_metMVA_met",  25, 155, 13, 3,   "MVA MET (GeV)",  "event_metMVA_met",     15, 155, 14, 0, 0, 0 },

  { "sqrt(JetPFCor_Pt[0]*JetPFCor_Pt[0]+JetPFCor_Pt[1]*JetPFCor_Pt[1]+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",
                          45, 205, 16, 3,   "Dijet System p_{T} (GeV)", "dijet_pt", 35, 205, 17, 0, 0, 0},  

  { "(JetPFCor_Eta[0]-JetPFCor_Eta[1])",
                         -3.0, 3.0, 15, 1,  "#Delta #eta (j,j)",    "deltaeta_jj",   -3.2,  3.2, 16, 0, 0, 0 },
  { "W_muon_pt",         15, 155, 14, 3,  "Muon p_{T} (GeV)",     "W_muon_pt",       25,  155, 13, 0, 0, 0 },
  { "W_muon_eta",      -2.7, 2.7, 18, 1,  "Muon #eta",            "W_muon_eta",    -2.4,  2.4, 16, 0, 0, 0 },

  { "W_electron_et",         30, 160, 13, 3,  "Electron E_{T} (GeV)",     "W_electron_et",       30,  160, 13, 0, 0, 0 },
  { "W_electron_eta",      -2.7, 2.7, 18, 1,  "Electron #eta",            "W_electron_eta",    -2.4,  2.4, 16, 0, 0, 0 },

//    { "event_nPV",      0, 50., 50, 1,  "Num PV",            "event_nPV",    0.,  50., 50, 0, 0, 0 },

  { "JetPFCor_dphiMET[0]", -3.14, 3.14, 16, 1, "#Delta #phi (Leading Jet, MET)", "deltaphi_jetldmet", -3.4, 3.4, 17, 0, 0, 0 },
#if 0
  { "sqrt((JetPFCor_Eta[0]-JetPFCor_Eta[1])**2+(abs(abs(abs(JetPFCor_Phi[0]-JetPFCor_Phi[1])-TMath::Pi())-TMath::Pi()))**2)",
                             0.4, 5.0, 23, 1, "#Delta R_{jj}",    "deltaRjj", 0.0, 5.4, 27, 0, 0, 1},
#endif
  { "", 0.0,0.0,0,0,"","",0.,0.,0,0,0,0 }
};


#endif // _COMMONCONTROLPLOTVARS_H
