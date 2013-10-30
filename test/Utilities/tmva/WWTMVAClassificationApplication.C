/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "/uscms/home/kukarzev/nobackup/root/root_v5.32.00/tmva/test/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;

void WWTMVAClassificationApplication( TString methodName = "Likelihood", double mH=400., int njets=2, TString chan="el" ) 
{   
#ifdef __CINT__
    gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif
    
    //---------------------------------------------------------------
    
    // This loads the library
    TMVA::Tools::Instance();
    
    // --------------------------------------------------------------------------------------------------
    
    // --- Create the Reader object
    
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    
    
    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    Float_t _WWpt, _Wpt, _MET, _LepCharge, _J1QGL, _J2QGL;
    Float_t _costheta1, _costheta2, _costhetaS, _Phi, _Phi2;
    reader->AddVariable( "WWpt := ptlvjj", &_WWpt );
    reader->AddVariable( "Wpt := W_pt", &_Wpt );
    reader->AddVariable( "MET := event_met_pfmet", &_MET );
    if (chan = "mu"){
        reader->AddVariable( "LepCharge := W_muon_charge", &_LepCharge );
    }
    else if (chan = "el"){
        reader->AddVariable( "LepCharge := W_electron_charge", &_LepCharge );        
    }
    else{
        std::cout << "Invalid channel!" << std::endl;
        return;
    }
    reader->AddVariable( "J1QGL := JetPFCor_QGLikelihood[0]", &_J1QGL );
    reader->AddVariable( "J2QGL := JetPFCor_QGLikelihood[1]", &_J2QGL );
    reader->AddVariable( "costheta1 := ang_ha", &_costheta1 );
    reader->AddVariable( "costheta2 := ang_hb", &_costheta2 );
    reader->AddVariable( "costhetaS := ang_hs", &_costhetaS );
    reader->AddVariable( "Phi := ang_phi", &_Phi );
    reader->AddVariable( "Phi2 := ang_phib", &_Phi2 );
    
    Int_t _run,_lumi,_event;
    Float_t _mjj,_mlvjj,_masslvjj;
    reader->AddSpectator("run := event_runNo", &_run);
    reader->AddSpectator("lumi := event_lumi", &_lumi);
    reader->AddSpectator("event := event_evtNo", &_event);
    reader->AddSpectator("mjj := Mass2j_PFCor", &_mjj);
    reader->AddSpectator("mlvjj := MassV2j_PFCor", &_mlvjj);
    reader->AddSpectator("masslvjj := masslvjj", &_masslvjj);

    //-------------------------------------------
    
    // --- Book the MVA methods
    
    TString dir    = "/uscms_data/d2/ntran/physics/VJets/WWlvjj/output_TMVAdiscriminant/WWlvjj_TMVAoutput/";
    TString prefix = "TMVAClassification";
    
    //TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
    char oweightfile[192];
    sprintf(oweightfile, "%s/m%3.0f/TMVAClassification_%3.0f_nJ%i_%s_%s.weights.xml",dir.Data(),mH,mH,njets,chan.Data(),methodName.Data());
    reader->BookMVA( methodName, oweightfile ); 
    
    /*
     // Prepare input tree (this must be replaced by your data source)
     // in this example, there is a toy tree with signal and one with background events
     // we'll later on use only the "signal" events for the test in this example.
     //   
     TFile *input(0);
     TString fname = "./tmva_example.root";   
     if (!gSystem->AccessPathName( fname )) 
     input = TFile::Open( fname ); // check if file in local directory exists
     else    
     input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
     
     if (!input) {
     std::cout << "ERROR: could not open data file" << std::endl;
     exit(1);
     }
     std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
     
     // --- Event loop
     
     // Prepare the event tree
     // - here the variable names have to corresponds to your tree
     // - you can use the same variables as above which is slightly faster,
     //   but of course you can use different ones and copy the values inside the event loop
     //
     std::cout << "--- Select signal sample" << std::endl;
     TTree* theTree = (TTree*)input->Get("TreeS");
     Float_t userVar1, userVar2;
     theTree->SetBranchAddress( "var1", &userVar1 );
     theTree->SetBranchAddress( "var2", &userVar2 );
     theTree->SetBranchAddress( "var3", &var3 );
     theTree->SetBranchAddress( "var4", &var4 );
     
     // Efficiency calculator for cut method
     Int_t    nSelCutsGA = 0;
     Double_t effS       = 0.7;
     
     std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests
     
     std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
     TStopwatch sw;
     sw.Start();
     for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {
     
     if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
     
     theTree->GetEntry(ievt);
     
     var1 = userVar1 + userVar2;
     var2 = userVar1 - userVar2;
     
     // --- Return the MVA outputs and fill into histograms
     
     if (Use["CutsGA"]) {
     // Cuts is a special case: give the desired signal efficienciy
     Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
     if (passed) nSelCutsGA++;
     }
     
     if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
     if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
     if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
     if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
     if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
     if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
     if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
     if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
     if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
     if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
     if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
     if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
     if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
     if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
     if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
     if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
     if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
     if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
     if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
     if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
     if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
     if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
     if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
     if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
     if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
     if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
     if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
     if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
     if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
     if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );
     
     // Retrieve also per-event error
     if (Use["PDEFoam"]) {
     Double_t val = reader->EvaluateMVA( "PDEFoam method" );
     Double_t err = reader->GetMVAError();
     histPDEFoam   ->Fill( val );
     histPDEFoamErr->Fill( err );         
     if (err>1.e-50) histPDEFoamSig->Fill( val/err );
     }         
     
     // Retrieve probability instead of MVA output
     if (Use["Fisher"])   {
     probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
     rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
     }
     }
     
     // Get elapsed time
     sw.Stop();
     std::cout << "--- End of event loop: "; sw.Print();
     
     // Get efficiency for cuts classifier
     if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
     << " (for a required signal efficiency of " << effS << ")" << std::endl;
     
     if (Use["CutsGA"]) {
     
     // test: retrieve cuts for particular signal efficiency
     // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
     TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;
     
     if (mcuts) {      
     std::vector<Double_t> cutsMin;
     std::vector<Double_t> cutsMax;
     mcuts->GetCuts( 0.7, cutsMin, cutsMax );
     std::cout << "--- -------------------------------------------------------------" << std::endl;
     std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
     for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
     std::cout << "... Cut: " 
     << cutsMin[ivar] 
     << " < \"" 
     << mcuts->GetInputVar(ivar)
     << "\" <= " 
     << cutsMax[ivar] << std::endl;
     }
     std::cout << "--- -------------------------------------------------------------" << std::endl;
     }
     }
     
     // --- Write histograms
     
     TFile *target  = new TFile( "TMVApp.root","RECREATE" );
     if (Use["Likelihood"   ])   histLk     ->Write();
     if (Use["LikelihoodD"  ])   histLkD    ->Write();
     if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
     if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
     if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
     if (Use["PDERS"        ])   histPD     ->Write();
     if (Use["PDERSD"       ])   histPDD    ->Write();
     if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
     if (Use["KNN"          ])   histKNN    ->Write();
     if (Use["HMatrix"      ])   histHm     ->Write();
     if (Use["Fisher"       ])   histFi     ->Write();
     if (Use["FisherG"      ])   histFiG    ->Write();
     if (Use["BoostedFisher"])   histFiB    ->Write();
     if (Use["LD"           ])   histLD     ->Write();
     if (Use["MLP"          ])   histNn     ->Write();
     if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
     if (Use["MLPBNN"       ])   histNnbnn  ->Write();
     if (Use["CFMlpANN"     ])   histNnC    ->Write();
     if (Use["TMlpANN"      ])   histNnT    ->Write();
     if (Use["BDT"          ])   histBdt    ->Write();
     if (Use["BDTD"         ])   histBdtD   ->Write();
     if (Use["BDTG"         ])   histBdtG   ->Write(); 
     if (Use["RuleFit"      ])   histRf     ->Write();
     if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
     if (Use["SVM_Poly"     ])   histSVMP   ->Write();
     if (Use["SVM_Lin"      ])   histSVML   ->Write();
     if (Use["FDA_MT"       ])   histFDAMT  ->Write();
     if (Use["FDA_GA"       ])   histFDAGA  ->Write();
     if (Use["Category"     ])   histCat    ->Write();
     if (Use["Plugin"       ])   histPBdt   ->Write();
     
     // Write also error and significance histos
     if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }
     
     // Write also probability hists
     if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
     target->Close();
     
     std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
     
     delete reader;
     
     std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
     */
} 
