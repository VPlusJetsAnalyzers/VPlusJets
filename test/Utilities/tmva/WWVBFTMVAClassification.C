// @(#)root/tmva $Id: WWVBFTMVAClassification.C,v 1.1 2011/12/16 14:26:35 kalanand Exp $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l TMVAClassification.C\(\"Fisher,Likelihood\"\)                       *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

// read input data file with ascii format (otherwise ROOT) ?
Bool_t ReadDataFromAsciiIFormat = kFALSE;

void WWVBFTMVAClassification( TString myMethodList = "", double mH=125., int njets = 2, TString chan="mu" ) 
{
	// The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
	// if you use your private .rootrc, or run from a different directory, please copy the 
	// corresponding lines from .rootrc

	// methods to be processed can be given as an argument; use format:
	//
	// mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
	//
	// if you like to use a method via the plugin mechanism, we recommend using
	// 
	// mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
	// (an example is given for using the BDT as plugin (see below),
	// but of course the real application is when you write your own
	// method based)

	// this loads the library
	TMVA::Tools::Instance();

	//---------------------------------------------------------------
	// default MVA methods to be trained + tested
	std::map<std::string,int> Use;

	Use["Cuts"]            = 0;
	Use["CutsD"]           = 0;
	Use["CutsPCA"]         = 0;
	Use["CutsGA"]          = 0;
	Use["CutsSA"]          = 0;
	// ---
	Use["Likelihood"]      = 0;
	Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
	Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
	Use["LikelihoodKDE"]   = 0;
	Use["LikelihoodMIX"]   = 0;
	// ---
	Use["PDERS"]           = 0;
	Use["PDERSD"]          = 0;
	Use["PDERSPCA"]        = 0;
	Use["PDERSkNN"]        = 0; // depreciated until further notice
	Use["PDEFoam"]         = 0;
	// --
	Use["KNN"]             = 0;
	// ---
	Use["HMatrix"]         = 0;
	Use["Fisher"]          = 0;
	Use["FisherG"]         = 0;
	Use["BoostedFisher"]   = 0;
	Use["LD"]              = 0;
	// ---
	Use["FDA_GA"]          = 0;
	Use["FDA_SA"]          = 0;
	Use["FDA_MC"]          = 0;
	Use["FDA_MT"]          = 0;
	Use["FDA_GAMT"]        = 0;
	Use["FDA_MCMT"]        = 0;
	// ---
	Use["MLP"]             = 0; // this is the recommended ANN
	Use["MLPBFGS"]         = 0; // recommended ANN with optional training method
	Use["CFMlpANN"]        = 0; // *** missing
	Use["TMlpANN"]         = 0; 
	// ---
	Use["SVM"]             = 0;
	// ---
	Use["BDT"]             = 1;
	Use["BDTD"]            = 0;
	Use["BDTG"]            = 0;
	Use["BDTB"]            = 0;
	// ---
	Use["RuleFit"]         = 0;
	// ---
	Use["Plugin"]          = 0;
	// ---------------------------------------------------------------

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassification" << std::endl;

	if (myMethodList != "") {
		for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

		std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
		for (UInt_t i=0; i<mlist.size(); i++) {
			std::string regMethod(mlist[i]);

			if (Use.find(regMethod) == Use.end()) {
				std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
				for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
				std::cout << std::endl;
				return;
			}
			Use[regMethod] = 1;
		}
	}

	// Create a new root output file.
	///////TString outfileName( "TMVA.root" );
	char outfileName[192];
	sprintf(outfileName,"TMVA_%3.0f_nJ%i_%s.root",mH,njets,chan.Data());
	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

	// Create the factory object. Later you can choose the methods
	// whose performance you'd like to investigate. The factory will
	// then run the performance analysis for you.
	//
	// The first argument is the base of the name of all the
	// weightfiles in the directory weight/ 
	//
	// The second argument is the output file for the training results
	// All TMVA output can be suppressed by removing the "!" (not) in 
	// front of the "Silent" argument in the option string
	char classifierName[192];
	sprintf(classifierName,"TMVAClassification_%3.0f_nJ%i_%s",mH,njets,chan.Data());
	TMVA::Factory *factory = new TMVA::Factory( classifierName, outputFile, 
			"!V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D" );
	//TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, 
	//                                           "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );

	// If you wish to modify default settings 
	// (please check "src/Config.h" to see all available global options)
	//    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
	//    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

	// Define the input variables that shall be used for the MVA training
       if (chan == "mu"){
                factory->AddVariable("LepCharge := W_muon_charge", 'F');
        }
        else if (chan == "el"){
                factory->AddVariable("LepCharge := W_electron_charge", 'F');
        }
        else{
                std::cout << "Invalid channel!" << std::endl;
                return;
        }

   	factory->AddVariable("Tagjet1pt := hvbf_aj_pt", 'F');
	factory->AddVariable("Tagjet2pt := hvbf_bj_pt", 'F');
    	//factory->AddVariable("Wjet1pt := hvbf_waj_pt", 'F');
	//factory->AddVariable("Wjet2pt := hvbf_wbj_pt", 'F');

       // factory->AddVariable("Tagjet1eta := hvbf_aj_eta", 'F');
       // factory->AddVariable("Tagjet2eta := hvbf_bj_eta", 'F');
       // factory->AddVariable("Wjet1eta := hvbf_waj_eta", 'F');
        //factory->AddVariable("Wjet2eta := hvbf_wbj_eta", 'F');

       //factory->AddVariable("Tagjet1phi := hvbf_aj_phi", 'F');
       // factory->AddVariable("Tagjet2phi := hvbf_bj_phi", 'F');
       // factory->AddVariable("Wjet1phi := hvbf_waj_phi", 'F');
       // factory->AddVariable("Wjet2phi := hvbf_wbj_phi", 'F')
//	factory->AddVariable("MET := event_met_pfmet", 'F');
/*
    factory->AddVariable("VBFcostheta_a := hvbf_wjj_ang_ha", 'F');
    factory->AddVariable("VBFcostheta_b := hvbf_wjj_ang_hb", 'F');
    factory->AddVariable("VBFcosthetaS := hvbf_wjj_ang_hs", 'F');
    factory->AddVariable("VBFPhi := hvbf_wjj_ang_phi", 'F');
    factory->AddVariable("VBFPhi2 := hvbf_wjj_ang_phib", 'F');
*/
// W jets varialbes

    factory->AddVariable("WJDeleaEta := hvbf_wjj_deta",'F');
    //factory->AddVariable("WJDeltaPhi := hvbf_wjj_dphi",'F');
    factory->AddVariable("WJMass := hvbf_wjj_m",'F'); 
//tag jets variables

    factory->AddVariable("VBFDeltaPhi := hvbf_jj_dphi",'F');
    factory->AddVariable("VBFTagJetMass := hvbf_jj_m",'F'); // preserved for fitting
  //  factory->AddVariable("WWrapidity := hvbf_lvjj_Rapidity",'F');
   // factory->AddVariable("WWpt := hvbf_lvjj_pt",'F');

// eta diff between ws and tag jets
//    factory->AddVariable("hW_tag1_deta := hvbf_hW_tag1_deta",'F');
//    factory->AddVariable("hW_tag2_deta := hvbf_hW_tag2_deta",'F');
//    factory->AddVariable("lW_tag1_deta := hvbf_lW_tag1_deta",'F');
//    factory->AddVariable("lW_tag2_deta := hvbf_lW_tag2_deta",'F');

//phi diff
//    factory->AddVariable("lW_hW_dphi := hvbf_lW_hW_deltaphi",'F');

    //factory->AddVariable("VBF_wjet1_QGd := hvbf_wjet1_QGd",'F');
    //factory->AddVariable("VBF_wjet2_QGd := hvbf_wjet2_QGd",'F');
  //  factory->AddVariable("VBF_tagjet1_QGd := hvbf_tagjet1_QGd",'F');
  //  factory->AddVariable("VBF_tagjet2_QGd := hvbf_tagjet2_QGd",'F');

   /* factory->AddVariable("VBF_wjet1_btagCSV := hvbf_wjet1_btagCSV",'F');
    factory->AddVariable("VBF_wjet2_btagCSV := hvbf_wjet2_btagCSV",'F');
    factory->AddVariable("VBF_tagjet1_btagCSV := hvbf_tagjet1_btagCSV",'F');
    factory->AddVariable("VBF_tagjet1_btagCSV := hvbf_tagjet1_btagCSV",'F');
*/
    factory->AddVariable("VBF_hvbf_lvjj_ZeppenField := hvbf_lvjj_ZeppenField",'F');

    //factory->AddVariable("VBFTagJet1Eta := hvbf_aj_eta",'F');
    //factory->AddVariable("VBFTagJet2Eta := hvbf_bj_eta",'F');


	// You can add so-called "Spectator variables", which are not used in the MVA training, 
	// but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
	// input variables, the response values of all trained MVAs, and the spectator variables
	factory->AddSpectator("run := event_runNo", "I");
	factory->AddSpectator("lumi := event_lumi", "I");
	factory->AddSpectator("event := event_evtNo", "I");
	factory->AddSpectator("mjj := Mass2j_PFCor", "F");
	factory->AddSpectator("mlvjj := MassV2j_PFCor", "F");
	factory->AddSpectator("hvbfevt := hvbf_event", "F");


	// read training and test data
	char signalOutputName[250];
	sprintf(signalOutputName,"/eos/uscms/store/user/lnujj/RDtrees_11Feb14/RD_%s_mh%3.0f_CMSSW532.root",chan.Data(),mH);

	TFile *input1 = TFile::Open( signalOutputName );

	char backgroundOutputName1[250];
	sprintf(backgroundOutputName1,"/eos/uscms/store/user/lnujj/RDtrees_11Feb14/RD_%s_W4Jets_CMSSW532_old.root",chan.Data());
	    char backgroundOutputName2[192];
	  sprintf(backgroundOutputName2,"/eos/uscms/store/user/lnujj/RDtrees_11Feb14/RD_%s_TTbar_CMSSW532.root",chan.Data());

	TFile *input2 = TFile::Open( backgroundOutputName1 );
	TFile *input3 = TFile::Open( backgroundOutputName2 );

	std::cout << "--- TMVAClassification : Using input file: " << input1->GetName() << std::endl;
	std::cout << "--- TMVAClassification : Using input file: " << input2->GetName() << std::endl;

	TTree *old_signal     = (TTree*)input1->Get("WJet");
	TTree *background1 = (TTree*)input2->Get("WJet");
	TFile newfile("small.root","recreate");

	 TTree * signal = old_signal->CloneTree(old_signal->GetEntries()*0.66);
	//TTree * signal = old_signal;
	 
	    TTree *background2 = (TTree*)input3->Get("WJet");
	// global event weights per tree (see below for setting event-wise weights)
	//    Double_t signalWeight     = 19.52*2.*.0314*10386/299968.;
	//    Double_t signalWeight     = 19.52*2.*.0314*10386/6893735.;
	//    Double_t signalWeight     = 225.197*10386./6893735.;
	//  Double_t backgroundWeight1 = 36257.2*10386./18353019.;

        Double_t signalWeight     = ((0.03907/2.0)* 193000)/414963;
	Double_t backgroundWeight1 = (214.0*193000)/4369420; 
        Double_t backgroundWeight2 = (225.197*193000)/6893735;


/*	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight1 = 1.0;
        Double_t backgroundWeight2 = 1.0;
*/

	//    Double_t backgroundWeight2 = 225.197*10386./6893735.;

	// ====== register trees ====================================================
	//
	// the following method is the prefered one:
	// you can add an arbitrary number of signal or background trees
	factory->AddSignalTree    ( signal,     signalWeight     );
	factory->AddBackgroundTree( background1, backgroundWeight1 );
	    factory->AddBackgroundTree( background2, backgroundWeight2 );

	// To give different trees for training and testing, do as follows:
	//    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
	//    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

	// Use the following code instead of the above two or four lines to add signal and background 
	// training and test events "by hand"
	// NOTE that in this case one should not give expressions (such as "var1+var2") in the input 
	//      variable definition, but simply compute the expression before adding the event
	// 
	//    // --- begin ----------------------------------------------------------
	//    std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
	//    Float_t  treevars[4];
	//    for (Int_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
	//    for (Int_t i=0; i<signal->GetEntries(); i++) {
	//       signal->GetEntry(i);
	//       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
	//       // add training and test events; here: first half is training, second is testing
	//       // note that the weight can also be event-wise	
	//       if (i < signal->GetEntries()/2) factory->AddSignalTrainingEvent( vars, signalWeight ); 
	//       else                            factory->AddSignalTestEvent    ( vars, signalWeight ); 
	//    }
	//
	//    for (Int_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
	//    for (Int_t i=0; i<background->GetEntries(); i++) {
	//       background->GetEntry(i); 
	//       for (Int_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
	//       // add training and test events; here: first half is training, second is testing
	//       // note that the weight can also be event-wise	
	//       if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight ); 
	//       else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight ); 
	//    }
	//    // --- end ------------------------------------------------------------
	//
	// ====== end of register trees ==============================================


	// This would set individual event weights (the variables defined in the 
	// expression need to exist in the original TTree)
	//    for signal    :
        char sig_w[200];
	if((mH==170 )||(mH==126 || (mH==125)) )
        {sprintf(sig_w,"effwt*puwt");}
        if((mH== 180) || (mH==190 )|| (mH==200) || (mH==250) ||(mH== 300) ||(mH== 350) )
        {
        sprintf(sig_w,"effwt*complexpolewtggH%3.0f",mH);
        }
        if((mH== 400) || (mH==450 )|| (mH==500) || (mH==550) ||(mH== 600)  )

        {
        sprintf(sig_w,"effwt*interferencewtggH%3.0f*complexpolewtggH%3.0f",mH,mH);
        }

//	cout<<"sig w exp "<<sig_w<<endl; 

  //      if(mH== (170 || 180 || 190 || 200 || 250 || 300 || 350 ))
//	factory->SetSignalWeightExpression("effwt*complexpolewtggh%3.0f",mH);
  //     else if (mH== (400 || 450 || 500 || 550 || 600 ))
//        factory->SetSignalWeightExpression("effwt*complexpolewtggH180");
        factory->SetSignalWeightExpression(sig_w);
	//    for background: 
	factory->SetBackgroundWeightExpression("effwt*puwt");
	// factory->SetBackgroundWeightExpression("weight");

	// Apply additional cuts on the signal and background samples (can be different)
	//   TCut mycuts = "abs(eta)>1.5"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	//    TCut mycutb = "abs(eta)>1.5"; // for example: TCut mycutb = "abs(var1)<0.5";


	TString mass4bodycut( "" );
	if(njets==2 ) {
		if(chan.Contains("mu")) {
			//if(mH==125.) mass4bodycut = "(ggh125_Mlvjj >100 && ggh125_Mlvjj<250)"; // 2j170mu
			//if(mH==160.) mass4bodycut = "(fit_mlvjj<200)"; // 2j170mu


                        if((mH==126.) || (mH==125)) mass4bodycut = "( hvbf_lvjj_m > 0 )"; // 2j170mu
			if(mH==170.) mass4bodycut = "((fit_mlvjj>165 )&& (fit_mlvjj<245) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j170mu =====
			if(mH==180.) mass4bodycut = "((fit_mlvjj>165) &&( fit_mlvjj<245)&& (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1)) "; // 2j180mu
			if(mH==190.) mass4bodycut = "((fit_mlvjj>165 )&& (fit_mlvjj<245) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j190mu
			if(mH==200.) mass4bodycut = "((fit_mlvjj>165 )&& (fit_mlvjj<245) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j200mu
			if(mH==250.) mass4bodycut = "((fit_mlvjj>200) && (fit_mlvjj<400) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j250mu
			if(mH==300.) mass4bodycut = "((fit_mlvjj>240) && (fit_mlvjj<400) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j300mu
			if(mH==350.) mass4bodycut = "((fit_mlvjj>250) && (fit_mlvjj<475) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j350mu
			if(mH==400.) mass4bodycut = "((fit_mlvjj>300) && (fit_mlvjj<600) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j400mu
			if(mH==450.) mass4bodycut = "((fit_mlvjj>340 ) && (fit_mlvjj<780) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j450mu
			if(mH==500.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j500mu
			if(mH==550.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j550mu
			if(mH==600.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) && (abs(JetPFCor_dphiMET[0])>0.4) && (abs(W_muon_eta)<2.1))"; // 2j600mu

		}

		if(chan.Contains("el")) {
			//if(mH==125.) mass4bodycut = "(ggh125_Mlvjj >100 && ggh125_Mlvjj<250)"; // 2j125el
		//	if(mH==160.) mass4bodycut = "(fit_mlvjj<200)"; // 2j160el
			if(mH==170.) mass4bodycut = "((fit_mlvjj>165)&&(fit_mlvjj<245)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j170el =====
			if(mH==180.) mass4bodycut = "((fit_mlvjj>165) && (fit_mlvjj<245) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j180el
			if(mH==190.) mass4bodycut = "((fit_mlvjj>165) && (fit_mlvjj<245)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)) "; // 2j190el
			if(mH==200.) mass4bodycut = "((fit_mlvjj>165) && (fit_mlvjj<245)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)) "; // 2j200el
			if(mH==250.) mass4bodycut = "((fit_mlvjj>200) && (fit_mlvjj<400) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j250el
			if(mH==300.) mass4bodycut = "((fit_mlvjj>240) && (fit_mlvjj<400)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)) "; // 2j300el
			if(mH==350.) mass4bodycut = "((fit_mlvjj>250) && (fit_mlvjj<475)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)) "; // 2j350el
			if(mH==400.) mass4bodycut = "((fit_mlvjj>300) && (fit_mlvjj<600)&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)) "; // 2j400el
			if(mH==450.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j450el
			if(mH==500.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j500el
			if(mH==550.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j550el
			if(mH==600.) mass4bodycut = "((fit_mlvjj>340) && (fit_mlvjj<780) &&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35))"; // 2j600el

		}
	}

/*	if(njets==3) {
		if(chan.Contains("mu")) {
			if(mH==125.) mass4bodycut = "(ggh125_Mlvjj >100 && ggh125_Mlvjj<250)"; // 2j170mu
			if(mH==160.) mass4bodycut = "(fit_mlvjj<200)"; // 2j170mu
			if(mH==170.) mass4bodycut = "(fit_mlvjj>150 && fit_mlvjj<271)"; // 3j170mu =====
			if(mH==180.) mass4bodycut = "(fit_mlvjj>175 && fit_mlvjj<284)"; // 3j180mu
			if(mH==190.) mass4bodycut = "(fit_mlvjj>185 && fit_mlvjj<290)"; // 3j190mu
			if(mH==200.) mass4bodycut = "(fit_mlvjj>188 && fit_mlvjj<293)"; // 3j200mu
			if(mH==250.) mass4bodycut = "(fit_mlvjj>216 && fit_mlvjj<300)"; // 3j250mu
			if(mH==300.) mass4bodycut = "(fit_mlvjj>241 && fit_mlvjj<355)"; // 3j300mu
			if(mH==350.) mass4bodycut = "(fit_mlvjj>269 && fit_mlvjj<407)"; // 3j350mu
			if(mH==400.) mass4bodycut = "(fit_mlvjj>300 && fit_mlvjj<465)"; // 3j400mu
			if(mH==450.) mass4bodycut = "(fit_mlvjj>332 && fit_mlvjj<518)"; // 3j450mu
			if(mH==500.) mass4bodycut = "(fit_mlvjj>362 && fit_mlvjj<569)"; // 3j500mu
			if(mH==550.) mass4bodycut = "(fit_mlvjj>398 && fit_mlvjj<616)"; // 3j550mu
			if(mH==600.) mass4bodycut = "(fit_mlvjj>419 && fit_mlvjj<660)"; // 3j600mu
		}
		if(chan.Contains("el")) {
			if(mH==125.) mass4bodycut = "(ggh125_Mlvjj >110 && ggh125_Mlvjj<250)"; // 2j125el
			if(mH==160.) mass4bodycut = "(fit_mlvjj<200)"; // 2j160el
			if(mH==170.) mass4bodycut = "(fit_mlvjj>150 && fit_mlvjj<271)"; // 3j170el =====
			if(mH==180.) mass4bodycut = "(fit_mlvjj>175 && fit_mlvjj<284)"; // 3j180el
			if(mH==190.) mass4bodycut = "(fit_mlvjj>185 && fit_mlvjj<290)"; // 3j190el
			if(mH==200.) mass4bodycut = "(fit_mlvjj>188 && fit_mlvjj<293)"; // 3j200el
			if(mH==250.) mass4bodycut = "(fit_mlvjj>216 && fit_mlvjj<300)"; // 3j250el
			if(mH==300.) mass4bodycut = "(fit_mlvjj>241 && fit_mlvjj<355)"; // 3j300el
			if(mH==350.) mass4bodycut = "(fit_mlvjj>269 && fit_mlvjj<407)"; // 3j350el
			if(mH==400.) mass4bodycut = "(fit_mlvjj>300 && fit_mlvjj<465)"; // 3j400el
			if(mH==450.) mass4bodycut = "(fit_mlvjj>332 && fit_mlvjj<518)"; // 3j450el
			if(mH==500.) mass4bodycut = "(fit_mlvjj>362 && fit_mlvjj<569)"; // 3j500el
			if(mH==550.) mass4bodycut = "(fit_mlvjj>398 && fit_mlvjj<616)"; // 3j550el
			if(mH==600.) mass4bodycut = "(fit_mlvjj>419 && fit_mlvjj<660)";  // 3j600el
		}

	}
*/
	char mycutschar[5000];
	/*if(mH != 125)
	  {sprintf(mycutschar,"ggdevt == %i &&(Mass2j_PFCor>65 && Mass2j_PFCor<95) && %s", njets, mass4bodycut);}
	  if(mH == 125)
	  {sprintf(mycutschar,"ggdevt == %i", njets);}

	  if(mH == 125 || mH == 160)
	  {sprintf(mycutschar,"(ggdevt == %i) && ( Mass2j_PFCor<95 ) && %s",njets,mass4bodycut);}
	// {sprintf(mycutschar,"(gg125evt == 1) && ( ggh125_wjj_m<80 ) && %s",mass4bodycut);}
	else
	//    {sprintf(mycutschar,"(gg125evt == 1) && (Mass2j_PFCor>65 && Mass2j_PFCor<95) && %s",mass4bodycut);}
*/
//	{sprintf(mycutschar,"(fit_status==0)&&(W_mt>30) &&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244)) && ((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244)) && ((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(Mass2j_PFCor>66)&&(Mass2j_PFCor<98)&&%s",mass4bodycut.Data());}
    {sprintf(mycutschar,"( (hvbf_event == 1) &&  (hvbf_wjj_m > 65.0 && hvbf_wjj_m < 95.0) && (hvbf_waj_pt>30.) && (hvbf_wbj_pt >30.) && (hvbf_aj_pt >30.) && (hvbf_bj_pt >30.)) && %s",mass4bodycut.Data());}
	

 
//	if(njets==2) 
/*	{

	//	if(mH == 125 || mH == 160)
	//	{sprintf(mycutschar,"(ggdevt == %i) && ( Mass2j_PFCor<95 ) && %s",njets,mass4bodycut);}
	//	else
	//	{
			if(chan.Contains("mu")) 
			{sprintf(mycutschar," ((fit_status==0) &&(W_mt>30) && ((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244)) && ((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244)) && ((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(abs(JetPFCor_dphiMET[0])>0.4)&&(abs(W_muon_eta)<2.1)&&(Mass2j_PFCor>66)&&(Mass2j_PFCor<98)) &&  %s",njets,mass4bodycut);}

			else if (chan.Contains("el"))
			{sprintf(mycutschar,"((fit_status==0)&&(W_mt>30)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(abs(JetPFCor_dphiMET[0])>0.8)&&(W_electron_pt>35)&&(Mass2j_PFCor>66)&&(Mass2j_PFCor<98))  &&  %s",njets,mass4bodycut);}
		//}
	}
*/

	TCut mycuts (mycutschar);

	// tell the factory to use all remaining events in the trees after training for testing:
	factory->PrepareTrainingAndTestTree( mycuts, mycuts,
			"nTrain_Signal=0:nTrain_Background=0:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=None:!V" );
        //          "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );


	//            "SplitMode=Random::SplitSeed=10:NormMode=EqualNumEvents" ); // Ruchika's settings

	//          "nTrain_Signal=0:nTrain_Background=6000:nTest_Signal=0:nTest_Background=6000:SplitMode=Random:NormMode=NumEvents:!V" );

	// If no numbers of events are given, half of the events in the tree are used for training, and 
	// the other half for testing:
	//factory->PrepareTrainingAndTestTree( mycuts, "SplitMode=random:!V" );  
	// To also specify the number of testing events, use:
	// factory->PrepareTrainingAndTestTree( mycuts, 
	//                                    "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );  

	// ---- Book MVA methods
	//
	// please lookup the various method configuration options in the corresponding cxx files, eg:
	// src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
	// it is possible to preset ranges in the option string in which the cut optimisation should be done:
	// "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

	// Cut optimisation
	if (Use["Cuts"])
		factory->BookMethod( TMVA::Types::kCuts, "Cuts", 
				"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

	if (Use["CutsD"])
		factory->BookMethod( TMVA::Types::kCuts, "CutsD", 
				"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

	if (Use["CutsPCA"])
		factory->BookMethod( TMVA::Types::kCuts, "CutsPCA", 
				"!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

	if (Use["CutsGA"])
		factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
				"H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

	if (Use["CutsSA"])
		factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
				"!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

	// Likelihood
	if (Use["Likelihood"])
	  factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", 
			       "!H:!V:!TransformOutput:PDFInterpol=Spline2"
			       ":NSmoothSig[0]=0:NSmoothBkg[0]=0"
			       ":NBinsSig[0]=5:NBinsBkg[0]=5"
			       ":PDFInterpolSig[0]=Spline0:PDFInterpolBkg[0]=Spline0"
			       ":NSmoothSig[1]=10:NSmoothBkg[1]=10"
			       ":NBinsSig[1]=20:NBinsBkg[1]=20"
			       ":NSmoothSig[2]=10:NSmoothBkg[2]=10"
			       ":NBinsSig[2]=20:NBinsBkg[2]=20"
			       ":NSmoothSig[3]=10:NSmoothBkg[3]=10"
			       ":NBinsSig[3]=20:NBinsBkg[3]=20"
			       ":NSmoothSig[4]=10:NSmoothBkg[4]=10"
			       ":NBinsSig[4]=20:NBinsBkg[4]=20"
                               ":NSmoothSig[5]=10:NSmoothBkg[5]=10"
                               ":NBinsSig[5]=20:NBinsBkg[5]=20"
                               ":NSmoothSig[6]=10:NSmoothBkg[6]=10"
                               ":NBinsSig[6]=20:NBinsBkg[6]=20"
                               ":NSmoothSig[7]=10:NSmoothBkg[7]=10"
                               ":NBinsSig[7]=20:NBinsBkg[7]=20"
/*                               ":NSmoothSig[8]=10:NSmoothBkg[8]=10"
                               ":NBinsSig[8]=20:NBinsBkg[8]=20"
                               ":NSmoothSig[9]=10:NSmoothBkg[9]=10"
                               ":NBinsSig[9]=20:NBinsBkg[9]=20"
                               ":NSmoothSig[10]=10:NSmoothBkg[10]=10"
                               ":NBinsSig[10]=20:NBinsBkg[10]=20"
                               ":NSmoothSig[11]=10:NSmoothBkg[11]=10"
                               ":NBinsSig[11]=20:NBinsBkg[11]=20"
                               ":NSmoothSig[12]=10:NSmoothBkg[12]=10"
                               ":NBinsSig[12]=20:NBinsBkg[12]=20"
                               ":NSmoothSig[13]=10:NSmoothBkg[13]=10"
                               ":NBinsSig[13]=20:NBinsBkg[13]=20"
                               ":NSmoothSig[14]=10:NSmoothBkg[14]=10"
                               ":NBinsSig[14]=20:NBinsBkg[14]=20"*/
			       ":NSmooth=10"
			       ":NAvEvtPerBin=10:CheckHist"
			      //  "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50"
			       ); 

	// test the decorrelated likelihood
	if (Use["LikelihoodD"])
		factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", 
				"!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ); 

	if (Use["LikelihoodPCA"])
		factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA", 
				"!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=5:NSmoothBkg[0]=5:NSmooth=5:NAvEvtPerBin=100:VarTransform=PCA" ); 

	// test the new kernel density estimator
	if (Use["LikelihoodKDE"])
		factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE", 
				"!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

	// test the mixed splines and kernel density estimator (depending on which variable)
	if (Use["LikelihoodMIX"])
		factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX", 
				"!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

	// test the multi-dimensional probability density estimator
	// here are the options strings for the MinMax and RMS methods, respectively:
	//      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );   
	//      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );   
	if (Use["PDERS"])
		factory->BookMethod( TMVA::Types::kPDERS, "PDERS", 
				"!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

	if (Use["PDERSkNN"])
		factory->BookMethod( TMVA::Types::kPDERS, "PDERSkNN", 
				"!H:!V:VolumeRangeMode=kNN:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

	if (Use["PDERSD"])
		factory->BookMethod( TMVA::Types::kPDERS, "PDERSD", 
				"!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

	if (Use["PDERSPCA"])
		factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA", 
				"!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

	// Multi-dimensional likelihood estimator using self-adapting phase-space binning
	if (Use["PDEFoam"])
		factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam", 
				"H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:CutNmin=T:Nmin=100:Kernel=None:Compress=T" );

	// K-Nearest Neighbour classifier (KNN)
	if (Use["KNN"])
		factory->BookMethod( TMVA::Types::kKNN, "KNN", 
				"H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
	// H-Matrix (chi2-squared) method
	if (Use["HMatrix"])
		factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" ); 

	// Fisher discriminant   
	if (Use["Fisher"])
		factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=60:NsmoothMVAPdf=10" );

	// Fisher with Gauss-transformed input variables
	if (Use["FisherG"])
		factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

	// Composite classifier: ensemble (tree) of boosted Fisher classifiers
	if (Use["BoostedFisher"])
		factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2");

	// Linear discriminant (same as Fisher)
	if (Use["LD"])
		factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None" );

	// Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
	if (Use["FDA_MC"])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

	if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

	if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

	if (Use["FDA_MT"])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

	if (Use["FDA_GAMT"])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

	if (Use["FDA_MCMT"])
		factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
				"H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

	// TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
	if (Use["MLP"])
		factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=10:EpochMonitoring" );

	if (Use["MLPBFGS"])
		factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=10:TrainingMethod=BFGS:!EpochMonitoring" );


	// CF(Clermont-Ferrand)ANN
	if (Use["CFMlpANN"])
		factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

	// Tmlp(Root)ANN
	if (Use["TMlpANN"])
		factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

	// Support Vector Machine
	if (Use["SVM"])
		factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

	// Boosted Decision Trees
	if (Use["BDTG"]) // Gradient Boost
		factory->BookMethod( TMVA::Types::kBDT, "BDTG", 
				"!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.30:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5" );

	if (Use["BDT"])  // Adaptive Boost
		factory->BookMethod( TMVA::Types::kBDT, "BDT", 
				"!H:!V:NTrees=400:nEventsMin=20:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
	//		"!H:!V:NTrees=200:nEventsMin=20:MaxDepth=9:BoostType=Bagging:UseYesNoLeaf=false:UseRandomisedTrees=true:SeparationType=GiniIndex:UseNvars=9:PruneMethod=NoPruning"); //Ruchika's settings

	if (Use["BDTB"]) // Bagging
		factory->BookMethod( TMVA::Types::kBDT, "BDTB", 
				"!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

	if (Use["BDTD"]) // Decorrelation + Adaptive Boost
		factory->BookMethod( TMVA::Types::kBDT, "BDTD", 
				"!H:!V:NTrees=400:nEventsMin=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:VarTransform=Decorrelate" );

	// RuleFit -- TMVA implementation of Friedman's method
	if (Use["RuleFit"])
		factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
				"H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

	// For an example of the category classifier, see: TMVAClassificationCategory

	// --------------------------------------------------------------------------------------------------

	// As an example how to use the ROOT plugin mechanism, book BDT via
	// plugin mechanism
	if (Use["Plugin"]) {
		//
		// first the plugin has to be defined, which can happen either through the following line in the local or global .rootrc:
		//
		// # plugin handler          plugin name(regexp) class to be instanciated library        constructor format
		// Plugin.TMVA@@MethodBase:  ^BDT                TMVA::MethodBDT          TMVA.1         "MethodBDT(TString,TString,DataSet&,TString)"
		// 
		// or by telling the global plugin manager directly
		gPluginMgr->AddHandler("TMVA@@MethodBase", "BDT", "TMVA::MethodBDT", "TMVA.1", "MethodBDT(TString,TString,DataSet&,TString)");
		factory->BookMethod( TMVA::Types::kPlugins, "BDT",
				"!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=CostComplexity:PruneStrength=50" );
	}

	// --------------------------------------------------------------------------------------------------

	// ---- Now you can tell the factory to train, test, and evaluate the MVAs

	// Train MVAs using the set of training events
	factory->TrainAllMethods();

	// ---- Evaluate all MVAs using the set of test events
	factory->TestAllMethods();

	// ----- Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();    

	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;      

	delete factory;

	// Launch the GUI for the root macros
	if (!gROOT->IsBatch()) TMVAGui( outfileName );

}
