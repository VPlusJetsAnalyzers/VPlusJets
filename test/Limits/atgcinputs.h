#ifndef _ATGCINPUTS_H
#define _ATGCINPUTS_H

// never uncomment this line!
//#define ISHWW // now you've done it!

#define NUMMASSPTS 1
const int masspts[NUMMASSPTS] = { 150 };

#define NUMCHAN 2
#define ELORMUCHAR 15

const char *channames[NUMCHAN] = {
  //"WVsemilepElboosted",
  "WVsemileptonic_elboosted",
  //"eldijet",
  "WVsemileptonic_muboosted",
  //"mudijet",
};

const TString AnalStr       = "ATGC";

const double mutrigeff      = 1.;
const double eltrigeff      = 1.;
const double sigtrigeffunc  = 0.01;
const double siglepteffunc  = 0.02;
const double siglumiunc     = 0.026;
const double signal_xs_unc  = 0.034;

// for cut-and-count limits:
const double dijetptmingev  = 350.;

const float LAMBDAZ_MIN = -0.06;
const float LAMBDAZ_MAX =  0.06;
const float LAMBDAZ_INC = 0.001; //  61 pts

const float dKG_MIN =  -0.12;
const float dKG_MAX =   0.12;
const float dKG_INC =   0.01;    // 25 pts

const float dg1_MIN =  -0.05;
const float dg1_MAX =   0.07;
const float dg1_INC =   0.002;   // x61 pts

// FILE/HISTO STRUCTURE: assumed same name objects contained in different files for the different inputs

//const char *dir = "/uscms_data/d2/kalanand/junk/vplusjets/CMSSW_4_2_8/src/ElectroWeakAnalysis/VPlusJets/test/TGC/";
const char *dir = ".";

const char *signalfmtstr_lzvsdkg  = "lambdaZ_%.3f_deltaKappaGamma_%.3f";
const char *signalfmtstr_lzvsdg1  = "lambdaZ_%.3f_deltaG1_%.3f";
const char *signalfmtstr_dkgvsdg1 = "deltaKappaGamma_%.3f_deltaG1_%.3f";

const char *dataobjname   = "data_obs";
const char *bkgdobjprefix = "bkg";
const char *signame       = "WWgammaZ";

#endif // _ATGCINPUTS_H
