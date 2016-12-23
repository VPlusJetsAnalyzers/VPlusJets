{
 gSystem->Load("libFWCoreFWLite.so");
 AutoLibraryLoader::enable();
 gSystem->Load("libDataFormatsFWLite.so");
 gROOT->SetStyle ("Plain");
 gSystem->Load("libRooFit") ;
 ////TAMUWW libraries:
//  /// CMSSW4XX
//  gSystem->Load("../lib/slc5_amd64_gcc434/libTAMUWWMEPATNtuple.so");
//  gSystem->Load("../lib/slc5_amd64_gcc434/libTAMUWWSpecialTools.so");
 /// CMSSW5XX
//  gSystem->Load("../lib/slc5_amd64_gcc462/libTAMUWWMEPATNtuple.so");
//  gSystem->Load("../lib/slc5_amd64_gcc462/libTAMUWWSpecialTools.so");
 using namespace RooFit ;
 //cout << "loaded TAMU stuff" << endl;

 if (gSystem->DynamicPathName("libFWCoreFWLite.so",true)) {
   cout << "loading FWlite.\n";
   gSystem->Load("libFWCoreFWLite.so");
   AutoLibraryLoader::enable();
   gSystem->Load("libDataFormatsFWLite.so");
   gSystem->Load("libDataFormatsPatCandidates.so");
   cout << "adding RooFit.\n";
   FILE *pipe = gSystem->OpenPipe("scram tool info roofitcore | grep INCLUDE | awk -F '[=]' '{print $2}'", "r");
   TString roofitinc;
   roofitinc.Gets(pipe);
   gSystem->ClosePipe(pipe);

   gROOT->GetInterpreter()->AddIncludePath(roofitinc);
   roofitinc.Prepend("-I\"");
   roofitinc += "\"";
   gSystem->AddIncludePath(roofitinc);
 }

//  TStyle *plain  = new TStyle("plain","Plain Style (no colors/fill areas)");
//  plain->SetCanvasBorderMode(0);
//  plain->SetPadBorderMode(0);
//  plain->SetPadColor(0);
//  plain->SetCanvasColor(0);
//  plain->SetCanvasDefH(350);
//  plain->SetCanvasDefW(450);
//  plain->SetTitleColor(1); //The X title of histograms
//  //plain->SetStatColor(0); 
//  plain->cd();//  this becomes now the current style gStyle
}
