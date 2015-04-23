#ifndef _HIGGSCONTROLPLOTVARS_H
#define _HIGGSCONTROLPLOTVARS_H

#include "plotvar_t.h"

////////////////////////////////////////////higgs
const plotVar_t higgsplotvars[] = {
  //    plotvar MINRange  MAXRange  NBINS  slog xlabel     outfile  AMINRange  AMAXRange ANBINS mva_in hplot drawleg

  //{ "MassV2j_PFCor",     140, 900, 19, 3, "m_{l#nujj} (GeV)",  "mlvjj",     100, 900, 20, 0, 1, 1 },
  //{ "MassV2j_PFCor_MVAMET",     140, 700, 14, 3, "m_{l#nujj} (GeV) MVA MET",  "mlvjjMVA",     100, 700, 15, 0, 1, 1 },
#if 0
  { "ptlvjj",              0, 150, 15, 3, "p_{T} of WW (GeV)", "ptlvjj",      0, 150, 15, 1, 0, 0 },
  { "ylvjj",            -2.5, 2.5, 25, 1, "WW rapidity" ,      "etalvjj",  -3.0, 3.0, 30, 1, 0, 0 },
  { "ang_ha",              0,   1,  5, 1, "Cos(#theta_{1})" ,  "ha",          0,   1,  5, 1, 0, 1 },
#endif
  { "ang_hb",              0,   1,  5, 1, "Cos(#theta_{2})" ,  "hb",          0,   1,  5, 1, 0, 0 },
  { "ang_hs",           -1.0,   1, 20, 1, "Cos(#theta*)" ,     "hs",       -1.0,   1, 20, 1, 0, 0 },
#if 0
  { "ang_phi",          -3.2, 3.2, 16, 6, "#Phi (rad)" ,       "phi",      -3.6, 3.6, 18, 1, 0, 0 },
  { "ang_phib",         -3.2, 3.2, 16, 6, "#Phi_{1} (rad)",    "phib",     -3.6, 3.6, 18, 1, 0, 0 },
  { "W_muon_charge",    -1.2, 1.2, 24, 1, "Muon Charge" ,      "charge",   -1.2, 1.2, 24, 1, 0, 0 },

  //    { "W_electron_charge",-1.2, 1.2, 24, 1, "Electron Charge" ,  "charge",   -1.2, 1.2, 24, 1, 0, 0 },
  { "mva2j500mu",       -0.2, 1.2, 28, 2, "Likelihood discriminant ( M_{H}=500 GeV )",  "mva2j500", -0.2, 1.2, 28, 0, 1, 1 },
  { "mva2j200mu",          0,   1, 25, 2, "Likelihood discriminant ( M_{H}=200 GeV ) ", "mva2j200", -0.1, 1.4, 30, 0, 1, 1 },
  //{ "mva2j350mu",          0,   1, 25, 2, "Likelihood discriminant ( M_{H}=350 GeV )",  "mva2j350_top", 0,   1, 10, 0, 0, 1 },

  { "mva2j180el",          0,   1, 25, 2, "Likelihood Output 2j el180 ",  "mva2j180_top", 0,   1, 10, 0, 0, 1 },
  { "qgld_Summer11CHS[0]", 0,   1, 25, 2, "Q/G Likelihood of Leading Jet","jetld_qgl",    0,   1, 25, 0, 0, 0 },
  { "qgld_Summer11CHS[1]", 0,   1, 25, 2, "Q/G Likelihood of Second Jet", "jetnt_qgl",    0,   1, 25, 0, 0, 0 },
#endif
  { "", 0.0,0.0,0,0,"","",0.,0.,0,0,0,0 }
};

#endif // _HIGGSCONTROLPLOTVARS_H
