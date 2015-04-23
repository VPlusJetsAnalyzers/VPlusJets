#ifndef _BOOSTEDCONTROLPLOTVARS_H
#define _BOOSTEDCONTROLPLOTVARS_H

#include "plotvar_t.h"

const plotVar_t boostedplotvars[] = {

  //    plotvar	MINRange  MAXRange  NBINS  slog xlabel outfile AMINRange  AMAXRange ANBINS mva_in hplot drawleg
  { "GroomedJet_CA8_pt_pr[0]",   200,   800, 30, 3, "pruned jet p_{T}",        "GroomedJet_pt_pr",      200,  800, 30, 0, 0, 1 },
  { "GroomedJet_CA8_mass_pr[0]",   0,   140, 28, 3, "pruned jet mass",        "GroomedJet_mass",      0,  140, 28, 0, 0, 0 },
//{ "GroomedJet_CA8_mass_pr[0]",   0,   140, 28, 3, "pruned jet mass",        "GroomedJet_mass",      65,  105, 28, 0, 0, 1 }, // signal region
//{ "(GroomedJet_CA8_pt[0]>150.0)+(GroomedJet_CA8_pt[1]>150.0)+(GroomedJet_CA8_pt[2]>150.0)+(GroomedJet_CA8_pt[3]>150.0)+(GroomedJet_CA8_pt[4]>150.0)+(GroomedJet_CA8_pt[5]>150.0)",     0, 5, 5, 3,   "Jet Multiplicity",            "CA8_JetMultiplicityg150",       0,  5, 5, 0, 0, 1 }, //jet multiplicity
//{ "(W_muon_pfiso_sumChargedHadronPt+max(0,W_muon_pfiso_sumNeutralHadronEt+W_muon_pfiso_sumPhotonEt-0.5*W_muon_pfiso_sumPUPt))/W_muon_pt",   0,   0.12, 24, 3, "RelIso_PF",        "RelIso",      0,  0.12, 24, 0, 0, 0 },

  { "W_pt",   200,   800, 30, 3, "Leptonic W  p_{T}", "W_pt",      200,  800, 30, 0, 0, 1 },
  { "W_muon_pt",         0, 520, 26, 3,  "Muon p_{T} (GeV)",     "W_muon_pt",       0,  520, 26, 0, 0, 1 },
  { "W_muon_eta",      -2.7, 2.7, 18, 1,  "Muon #eta",            "W_muon_eta",    -2.4,  2.4, 16, 0, 0, 0 },
  { "W_electron_pt",         0, 520, 26, 3,  "Electron p_{T} (GeV)",     "W_elec_pt",       0,  520, 26, 0, 0, 1 },
  { "W_electron_eta",      -2.7, 2.7, 18, 1,  "Electron #eta",            "W_elec_eta",    -2.7,  2.7, 18, 0, 0, 0 },
//{ "event_nPV",      0, 50., 50, 1,  "Num PV",            "event_nPV",    0.,  50., 50, 0, 0, 0 },
  { "event_met_pfmet",  0, 600, 24, 3,   "PF MET (GeV)",  "event_met_pfmet",     0, 600, 24, 0, 0, 0 },
//{ "GroomedJet_numberbjets_csvm",     0, 5, 5, 3,   "nbjets_csmv",            "GroomedJet_numberbjets_csvm",       0,  5, 5, 0, 1, 1 },
//{ "GroomedJet_CA8_pt[0]",   200,   800, 40, 3, "ungroomed jet p_{T}",        "GroomedJet_pt",      200,  800, 40, 0, 0, 0 },

  { "", 0.0,0.0,0,0,"","",0.,0.,0,0,0,0 }
};


#endif // _BOOSTEDCONTROLPLOTVARS_H
