////////////////////////////////////////////////////////////////////////////////////////////////
////                    CMS - Semileptonic Channel (jj+l+nu final state)                    ////
////                       Diboson, Higgs and  Matrix Element Analyses                      ////
////    Created by Osipenkov, Ilya (Texas A&M) : ilyao@fnal.gov, ilyao@physics.tamu.edu     ////
////////////////////////////////////////////////////////////////////////////////////////////////
////                    General Tools for plotting Monte Carlo and Data                     ////
////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <sstream>
#include <numeric>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
//#include <TH1.h>
#include <TH1F.h>
//#include <THStack.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "CMS_lumi.h"
#include "TMath.h"
//#include "TF1.h"
//#include "TMinuit.h"
//#include "TLatex.h"

#include "inputDibosonMCInfoFitAndValidationResults.cpp"
using namespace dibosonAnalysis;
// using std::pair;
// using std::make_pair;
// using std::vector;
// using std::map;

///Standard-library names:
using std::cout;
using std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////*********      Global Function     *********//////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////  Histogram And Canvas Manipulation:  /////////////////////////////
TH1F* addSecondHistToFirst(TH1F* histIn1, TH1F* histIn2){
  /// Updates histIn1 and returns a pointer to it
  if ( histIn2==NULL ) { return histIn1; } //nothing to add
  if ( histIn1==NULL ) { 
    histIn1 = new TH1F(*histIn2); 
  } else {
    histIn1->Add(histIn2);
  }
  return histIn1;
}

void saveCanvasAs(TCanvas *cnvIn, TString cnvName) { cnvIn->SaveAs(cnvName); }

void CMS_lumi( TPad* pad, int iPeriod, int iPosX )
///Standard CMS Lumi & Tdr style format
{            
  bool outOfFrame    = false;
  if( iPosX/10==0 ) 
    {
      outOfFrame = true;
    }
  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float H = pad->GetWh();
  float W = pad->GetWw();
  float l = pad->GetLeftMargin();
  float t = pad->GetTopMargin();
  float r = pad->GetRightMargin();
  float b = pad->GetBottomMargin();
  //  float e = 0.025;

  pad->cd();

  TString lumiText;
  if( iPeriod==1 )
    {
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==2 )
    {
      lumiText += lumi_8TeV;
      lumiText += " (8 TeV)";
    }
  else if( iPeriod==3 ) 
    {
      lumiText = lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
    }
  else if ( iPeriod==4 )
    {
      lumiText += lumi_13TeV;
      lumiText += " (13 TeV)";
    }
  else if ( iPeriod==7 )
    { 
      if( outOfFrame ) lumiText += "#scale[0.85]{";
      lumiText += lumi_13TeV; 
      lumiText += " (13 TeV)";
      lumiText += " + ";
      lumiText += lumi_8TeV; 
      lumiText += " (8 TeV)";
      lumiText += " + ";
      lumiText += lumi_7TeV;
      lumiText += " (7 TeV)";
      if( outOfFrame) lumiText += "}";
    }
  else if ( iPeriod==12 )
    {
      lumiText += "8 TeV";
    }
   
  std::cout << lumiText << std::endl;

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

  if( outOfFrame )
    {
      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);
    }
  
  pad->cd();

  float posX_=0;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  if( !outOfFrame )
    {
      if( drawLogo )
	{
	  posX_ =   l + 0.045*(1-l-r)*W/H;
	  posY_ = 1-t - 0.045*(1-t-b);
	  float xl_0 = posX_;
	  float yl_0 = posY_ - 0.15;
	  float xl_1 = posX_ + 0.15*H/W;
	  float yl_1 = posY_;
	  TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
	  TPad* pad_logo = new TPad("logo","logo", xl_0, yl_0, xl_1, yl_1 );
	  pad_logo->Draw();
	  pad_logo->cd();
	  CMS_logo->Draw("X");
	  pad_logo->Modified();
	  pad->cd();
	}
      else
	{
	  latex.SetTextFont(cmsTextFont);
	  latex.SetTextSize(cmsTextSize*t);
	  latex.SetTextAlign(align_);
	  latex.DrawLatex(posX_, posY_, cmsText);
	  if( writeExtraText ) 
	    {
	      latex.SetTextFont(extraTextFont);
	      latex.SetTextAlign(align_);
	      latex.SetTextSize(extraTextSize*t);
	      latex.DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, extraText);
	    }
	}
    }
  else if( writeExtraText )
    {
      if( iPosX==0) 
	{
	  posX_ =   l +  relPosX*(1-l-r);
	  posY_ =   1-t+lumiTextOffset*t;
	}
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);      
    }
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////*********      Class Definitions     *********/////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////  Graphics Manipulation Classes:  ///////////////////////////////
class CanvasProperties{
public:
  CanvasProperties() {}
  CanvasProperties(TString cnvNameIn, int cnvTopXIn, int cnvTopYIn, int cnvWidthIn, int cnvHeightIn, bool cmslumi_writeExtraTextIn=true, TString cmslumi_extraTextIn="Preliminary", TString cmslumi_lumi_8TeVIn="19 fb^{-1}", float cmslumi_relPosXIn=0.3): cnvName(cnvNameIn),cnvTopX(cnvTopXIn),cnvTopY(cnvTopYIn),cnvWidth(cnvWidthIn),cnvHeight(cnvHeightIn),cmslumi_writeExtraText(cmslumi_writeExtraTextIn),cmslumi_extraText(cmslumi_extraTextIn),cmslumi_lumi_8TeV(cmslumi_lumi_8TeVIn),cmslumi_relPosX(cmslumi_relPosXIn) {}

  TString getCnvName() {return cnvName;}
  int getCnvTopX() {return cnvTopX;}
  int getCnvTopY() {return cnvTopY;}
  int getCnvWidth() {return cnvWidth;}
  int getCnvHeight() {return cnvHeight;}

  TCanvas* makeCnv(TString="");
  void setCMSLumiGlobalVars() {
    writeExtraText = cmslumi_writeExtraText;
    extraText = cmslumi_extraText;
    lumi_8TeV = cmslumi_lumi_8TeV;
    relPosX = cmslumi_relPosX;
  }
  void setVarsAndRunCMSLumi(TPad* pad, int iPeriod=3, int iPosX=10) {
    setCMSLumiGlobalVars();
    CMS_lumi(pad,iPeriod,iPosX);
  }

private:
  TString cnvName;
  int cnvTopX;
  int cnvTopY;
  int cnvWidth;
  int cnvHeight;
  ///Set prior to running CMS_lumi
  bool cmslumi_writeExtraText; // if extra text
  TString cmslumi_extraText;  // default extra text is "Preliminary"
  TString cmslumi_lumi_8TeV; // default is "19.7 fb^{-1}"
  float cmslumi_relPosX; // Adjust the label postion
};

TCanvas* CanvasProperties::makeCnv(TString cnvNameInput) {
  ///Make a canvas with the properties described in the class and an optional name
  TString canvasName=cnvNameInput;
  if ( canvasName=="" ) { canvasName=cnvName; }
  TCanvas* cnvTemp = new TCanvas(canvasName,canvasName,cnvTopX,cnvTopY,cnvWidth,cnvHeight);
  return cnvTemp;
}


class PlotProperties{
public:
  PlotProperties() {}
  PlotProperties(TString plotTitleIn, double plotTitleXIn, double plotTitleWIn, double plotTitleYIn, TString xAxisTitleIn, double xAxisTitleOffsetIn, TString yAxisTitleIn, double yAxisTitleOffsetIn): plotTitle(plotTitleIn), plotTitleX(plotTitleXIn), plotTitleW(plotTitleWIn), plotTitleY(plotTitleYIn), xAxisTitle(xAxisTitleIn), xAxisTitleOffset(xAxisTitleOffsetIn), yAxisTitle(yAxisTitleIn), yAxisTitleOffset(yAxisTitleOffsetIn) {}

  void setPlotProperties(PlotProperties plotPropIn) { *this=PlotProperties(plotPropIn); }
  TString getPlotTitle() {return plotTitle;}
  void setPlotTitle(TString plotTitleIn) { plotTitle=plotTitleIn; }
  double getPlotTitleX() {return plotTitleX;}
  double getPlotTitleY() {return plotTitleY;}
  TString getXAxisTitle() {return xAxisTitle;}
  double getXAxisTitleOffset() {return xAxisTitleOffset;}
  TString getYAxisTitle() {return yAxisTitle;}
  double getYAxisTitleOffset() {return yAxisTitleOffset;}

  TH1F* labelHistogram(TH1F*);


private:
  TString plotTitle;
  double plotTitleX;
  double plotTitleW;
  double plotTitleY;
  TString xAxisTitle;
  double xAxisTitleOffset;
  TString yAxisTitle;
  double yAxisTitleOffset;
};

TH1F* PlotProperties::labelHistogram(TH1F* histToLabel)
//// Add a title and axis labels to a histogram and return a pointer to it
{
  histToLabel->SetTitle(plotTitle);
  if ( plotTitleX > 0 ) {
    gStyle->SetTitleX(plotTitleX);
  }
  if ( plotTitleW > 0 ) {
    gStyle->SetTitleW(plotTitleW);
  }
  if ( plotTitleY > 0 ) {
    gStyle->SetTitleY(plotTitleY);
  }
  histToLabel->SetXTitle(xAxisTitle);
  if ( xAxisTitleOffset > 0.0 ) {
    histToLabel->GetXaxis()->SetTitleOffset(xAxisTitleOffset);
  }
  histToLabel->SetYTitle(yAxisTitle);
  if ( yAxisTitleOffset > 0.0 ) {
    histToLabel->GetYaxis()->SetTitleOffset(yAxisTitleOffset);
  }
  return histToLabel;
}

template<class DesignationType> class StackPlotProperties: public CanvasProperties, public PlotProperties{
public:
  StackPlotProperties(): CanvasProperties(), PlotProperties() {
    plotMinimum=0.0;
    plotMaximumIsUsed=false;
  }
  StackPlotProperties(CanvasProperties cnvPropIn, PlotProperties plotPropIn, pair<double,double> lgndBottomLeftIn, pair<double,double> lgndTopRightIn, Color_t lgndFillColorIn, vector<  pair< DesignationType,pair<TString,Color_t> >  > nameColorPairInPlotOrder ): CanvasProperties(cnvPropIn), PlotProperties(plotPropIn) {
    lgndFillColor=lgndFillColorIn;
    lgndX1=lgndBottomLeftIn.first;
    lgndY1=lgndBottomLeftIn.second;
    lgndX2=lgndTopRightIn.first;
    lgndY2=lgndTopRightIn.second;
    designationsTopToBottom.clear();
    for ( typename vector<  pair< DesignationType,pair<TString,Color_t> >  > ::iterator inputProcessIter = nameColorPairInPlotOrder.begin();  inputProcessIter!= nameColorPairInPlotOrder.end(); ++inputProcessIter ) {
      DesignationType procDesignation = inputProcessIter->first;
      TString procName = (inputProcessIter->second).first;
      Color_t procColor = (inputProcessIter->second).second;
      designationsTopToBottom.push_back(procDesignation);
      processNames[procDesignation]=procName;
      processColors[procDesignation]=procColor;
    }
    plotMinimum=0.0;
    plotMaximumIsUsed=false;
  }

  double getLgndX1() {return lgndX1;}
  double getLgndY1() {return lgndY1;}
  double getLgndX2() {return lgndX2;}
  double getLgndY2() {return lgndY2;}
  Color_t getLgndFillColor() {return lgndFillColor;}

  vector<DesignationType> getDesignationsTopToBottom() {return designationsTopToBottom;}
  TString getProcessName(DesignationType Designation) {
    if ( processNames.count(Designation) ) {
      return processNames[Designation];
    } else {
      cout << "Unable to find the process name for the specified designation" << endl;
      return "";
    }
  }
  Color_t getProcessColor(DesignationType Designation) {
    if ( processColors.count(Designation) ) {
      return processColors[Designation];
    } else {
      cout << "Unable to find the process color for the specified designation" << endl;
      return kYellow;
    }
  }

  void setPlotMinimum(double plotMinimumIn) { plotMinimum=plotMinimumIn; }
  double getPlotMinimum() { return plotMinimum; }
  void setPlotMaximum(double plotMaximumIn) {
    plotMaximumIsUsed=true;
    plotMaximum=plotMaximumIn;
  }
  bool usePlotMaximum() { return plotMaximumIsUsed; }
  double getPlotMaximum() { return plotMaximum; }


private:
  double lgndX1;
  double lgndY1;
  double lgndX2;
  double lgndY2;
  Color_t lgndFillColor;

  vector< DesignationType > designationsTopToBottom;//List the order of processes appearing in the stack plot (assuming they are included) top to bottom
  map< DesignationType,TString > processNames;
  map< DesignationType,Color_t > processColors;

  double plotMinimum;
  bool plotMaximumIsUsed;
  double plotMaximum;

};

////////////////////////  Classes For Reading And Combining Processes:  ////////////////////////
///////////////  Step 1 - set up the infrastructure for working with a single file.
class AnalysisVariable{
  ///The variable to draw the plots & perform analysis in
public:
  AnalysisVariable() {};
  AnalysisVariable(TString nameIn, double minIn, double maxIn, double numberOfBinsIn, TString optionalCuts=""): name(nameIn), min(minIn), max(maxIn), numberOfBins(numberOfBinsIn) {
    ///Construct the cut string to include minimum, maximum as well as optional cuts
    std::ostringstream cutStream;
    cutStream << "(" << name << ">" << min << ")&&(" << name << "<" << max << ")";
    if ( optionalCuts!="" ) { cutStream << "&&" << optionalCuts; }
    associatedCutString = cutStream.str();
  }

  TString getName() { return name; }
  double getMin() { return min; }
  double getMax() { return max; }
  double getNumberOfBins() { return numberOfBins; }
  TString getAssociatedCutString() { return associatedCutString; }

  void printSummary() const ;

private:
  TString name;
  double min;
  double max;
  double numberOfBins;
  TString associatedCutString;
};

void AnalysisVariable::printSummary() const{
  cout << "Variable=" << name << ", range=(" << min << "," << max << ") in " << numberOfBins << " bins, with associated cuts = " << associatedCutString << endl;
}

class SingleFileInput{
public:
  SingleFileInput() { fileHist=NULL; }
  SingleFileInput(TString filePathIn, AnalysisVariable variableIn, double kFactorXLumiXCrossSecDivByNumGeneratedIn, TString baseCutsIn="",  TString eventWeightIn="", TString treeNameIn="WJet"): filePath(filePathIn), treeName(treeNameIn), kFactorXLumiXCrossSecDivByNumGenerated(kFactorXLumiXCrossSecDivByNumGeneratedIn), variable(variableIn), baseCuts(baseCutsIn), eventWeight(eventWeightIn) { 
    additionalScalingFactor=1.0;
    fileHist=NULL;
  }
  ///Need to define the copy constructor, assignment operator and destructor to make sure the histogram (and not the pointer) is properly copied/deleted.
  SingleFileInput(const SingleFileInput& fileInputIn) {
    filePath = fileInputIn.filePath;
    treeName = fileInputIn.treeName;
    kFactorXLumiXCrossSecDivByNumGenerated = fileInputIn.kFactorXLumiXCrossSecDivByNumGenerated;
    variable = fileInputIn.variable;
    baseCuts = fileInputIn.baseCuts;
    eventWeight = fileInputIn.eventWeight;
    additionalScalingFactor = fileInputIn.additionalScalingFactor;
    TH1F* fileHistIn = fileInputIn.fileHist;
    if ( fileHistIn!=NULL ) {
      fileHist=new TH1F(*fileHistIn);
    } else {
      fileHist=NULL;
    }
  }
  SingleFileInput& operator=(const SingleFileInput &rhs) {
    filePath = rhs.filePath;
    treeName = rhs.treeName;
    kFactorXLumiXCrossSecDivByNumGenerated = rhs.kFactorXLumiXCrossSecDivByNumGenerated;
    variable = rhs.variable;
    baseCuts = rhs.baseCuts;
    eventWeight = rhs.eventWeight;
    additionalScalingFactor = rhs.additionalScalingFactor;
    if ( fileHist!=NULL ) {delete fileHist;}
    TH1F* fileHistIn = rhs.fileHist;
    if ( fileHistIn!=NULL ) {
      fileHist=new TH1F(*fileHistIn);
    } else {
      fileHist=NULL;
    }
    return *this;
  }
  ~SingleFileInput() { 
    //cout << "initating the single file destructor" << endl;
    if ( fileHist!=NULL ) { 
      //cout << "destroy" << endl;
      delete fileHist;
    }
    //cout << "finished" << endl;
  }

  void setBaseCuts(TString& baseCutsIn) { baseCuts=baseCutsIn; }
  void setVariable(const AnalysisVariable& variableIn) { variable=variableIn; }
  void setEventWeight(TString& eventWeightIn) { eventWeight=eventWeightIn; }
  void multiplyAdditionalScalingFactor(double multiple) { additionalScalingFactor=additionalScalingFactor*multiple; }
  void resetAddidionalScalingFactor() { additionalScalingFactor=1.0; }

  AnalysisVariable getAnalysisVariable() { return variable; }
  TString getFilePath() { return filePath; }
  TString getTreeName() { return treeName; }
  double getKFactorXLumiDivByNumGenerated() { return kFactorXLumiXCrossSecDivByNumGenerated; }
  TString getBaseCuts() { return baseCuts; }
  // TString getEventWeight() { return eventWeight; }

  void printSummary() const;

  TString produceCutStringWithoutWeight();
  TString produceCutString();
  void remakeFileHist();
  void makeFileHistIfNull();
  void scaleBy(double scalingFactor) {
    multiplyAdditionalScalingFactor(scalingFactor);
    if ( fileHist ) { fileHist->Scale(scalingFactor); }
  }
  TH1F* getFileHistPtr() const { return fileHist; }
  int getNTreeEntries() const { return nTreeEntries; }

private:
  TString filePath; //Should be set when the object is created and not change
  TString treeName;
  int nTreeEntries;
  double kFactorXLumiXCrossSecDivByNumGenerated;//Standard weight - the number expected will just be NumberOfEventsAfterCuts * kFactorXLumiXCrossSecDivByNumGenerated.
  AnalysisVariable variable;
  TString baseCuts;
  TString eventWeight;
  double additionalScalingFactor;//Keeps track of any additional scaling modifications (e.g. based on fit results) and is applied when making the histogram
  TH1F* fileHist;
};

TH1F* addToHistAHistFromFile(TH1F* histIn, const SingleFileInput& fileIn) {
  /// Updates histIn and returns a pointer to it
  TH1F* fileHist = fileIn.getFileHistPtr();
  histIn = addSecondHistToFirst(histIn,fileHist);
  return histIn;
}

void SingleFileInput::printSummary() const {
  cout << "Summary for file: " << filePath << endl;
  cout << "kFactorXLumiXCrossSecDivByNumGenerated = " << kFactorXLumiXCrossSecDivByNumGenerated << endl;
  cout << "baseCuts = " << baseCuts << endl;
  cout << "eventWeight = " << eventWeight << endl;
  cout << "additionalScalingFactor = " << additionalScalingFactor << endl;
  variable.printSummary();
  cout << endl;
}

TString SingleFileInput::produceCutStringWithoutWeight() {
  ///Construct the full cut string taking into account baseCuts and the range of the variable
  TString cutString = baseCuts;
  cutString = cutString + "&&";
  cutString = cutString + variable.getAssociatedCutString();
  return cutString;
}

TString SingleFileInput::produceCutString() {
  ///Construct the full cut string taking into account baseCuts, the range of the variable and eventWeight when making the histogram
  TString cutString = produceCutStringWithoutWeight();
  if ( eventWeight!="" ) {
    cutString = "( " + cutString;
    cutString = cutString+ " )*";
    cutString = cutString+eventWeight;
  }
  return cutString;
}

void SingleFileInput::remakeFileHist() {
  //Open the files
  TFile* eventFile = new TFile(filePath);
  TTree* eventTree = (TTree*)eventFile->Get(treeName);
  nTreeEntries = eventTree->GetEntries();

  //Fill a temporary histogram
  TH1F* htemp = new TH1F("htemp","htemp",variable.getNumberOfBins(),variable.getMin(),variable.getMax());
  TString htempFillString=variable.getName() + ">>htemp";
  eventTree->Draw(htempFillString,produceCutString());
  htemp->Sumw2();
  htemp->Scale(kFactorXLumiXCrossSecDivByNumGenerated*additionalScalingFactor);
  if ( fileHist!=NULL ) { delete fileHist; }//delete the previous
  fileHist = htemp;
}

void SingleFileInput::makeFileHistIfNull() {
  if ( fileHist==NULL ) { remakeFileHist(); }
}

///////////////  Step 2 - combine the files into a single process.
template<class DesignationType> class SingleProcess{
/// Manages a signle physics process, which may be composed from a number of files with different weights, cuts, etc.
/// Is set up as a template class in order to processes with different designation type, which will usually be a string (e.g. "diboson") or an enum created for the analysis. The DesignationType should have != and < operators defined.
public:
  SingleProcess() { processHist=NULL; }
  SingleProcess(DesignationType processDesignationIn): processDesignation(processDesignationIn) { processHist=NULL; }
  ///Need to define the copy constructor, assignment operator and destructor to make sure the histogram (and not the pointer) is properly copied/deleted.
  SingleProcess(const SingleProcess& processIn) {
    processDesignation=processIn.processDesignation;
    allProcessFiles=processIn.allProcessFiles;
    TH1F* processHistIn = processIn.processHist;
    if ( processHistIn!=NULL ) {
      processHist=new TH1F(*processHistIn);
    } else {
      processHist=NULL;
    }
  }
  SingleProcess& operator=(const SingleProcess &rhs) {
    processDesignation=rhs.processDesignation;
    allProcessFiles=rhs.allProcessFiles;
    if ( processHist!=NULL ) {delete processHist;}
    addProcessHistAndHist(rhs.processHist);
    return *this;
  }
  ~SingleProcess() { 
    //cout << "initating the single process destructor" << endl;
    if ( processHist!=NULL ) { 
      //cout << "destroy" << endl;
      delete processHist;
    }
    //cout << "finished" << endl;
  }

  ///Overloaded +=,+ 
  SingleProcess& operator+=(const SingleFileInput& );
  SingleProcess& operator+=(const SingleProcess& );

  ///Various tools
  void setDesignation(DesignationType designationIn) { processDesignation=designationIn; }
  void remakeIndividualFileHistograms();
  void fillNullIndividualFileHistograms();
  void combineFileHistogramsToANewProcessHistogram();
  void remakeProcessAndFileHistograms();
  void rescaleProcessAndFilesBy(double);
  void setVariableInAllFilesTo(const AnalysisVariable& );

  vector<SingleFileInput> getAllProcessFiles() { return allProcessFiles; }
  TH1F* getProcessHistPtr() { return processHist; }
  TH1F* getProcessHistPtr() const { return processHist; }
  double getProcessEvents() const { return processHist->Integral(); }
  DesignationType getProcessDesignation() const { return processDesignation; }

  bool isEmpty() const { return allProcessFiles.empty(); }
  void printSummary() const;

protected:
  vector<SingleFileInput> allProcessFiles;

private:
  void addProcessHistAndHist(TH1F* histIn) { 
    processHist = addSecondHistToFirst(processHist,histIn);
  }

  ///Class members
  DesignationType processDesignation;
  TH1F* processHist;
};

/// Overload the += and + to do most of the work manipulating SingleProcess objects
template<class DesignationType> SingleProcess<DesignationType>& SingleProcess<DesignationType>::operator+=(const SingleFileInput& fileInput) {
  allProcessFiles.push_back(fileInput);
  processHist=addToHistAHistFromFile(processHist,fileInput);
  return *this;
}

template<class DesignationType> SingleProcess<DesignationType>& SingleProcess<DesignationType>::operator+=(const SingleProcess<DesignationType>& processIn) {
  //cout << "1, designation=" << processDesignation << endl;
  if ( allProcessFiles.empty() ) {
    processDesignation=processIn.processDesignation;
  } else {
    if ( processDesignation!=processIn.processDesignation ) { cout << "Warning: adding processes with different designations" << endl; }
  }
  //cout << "2, designation=" << processDesignation << endl;
  std::copy(processIn.allProcessFiles.begin(), processIn.allProcessFiles.end(), std::back_inserter(allProcessFiles));
  //cout << "3, designation=" << processDesignation << endl;
  addProcessHistAndHist(processIn.processHist);
  //cout << "4, designation=" << processDesignation << endl;
  return *this;
}

template<class DesignationType> SingleProcess<DesignationType> operator+(const SingleProcess<DesignationType>& lhsSingleProcess, const SingleFileInput& rhsFileInput) {
  SingleProcess<DesignationType> processTemp(lhsSingleProcess);
  processTemp+=rhsFileInput;
  return processTemp;
}

template<class DesignationType> SingleProcess<DesignationType> operator+(const SingleFileInput& lhsFileInput, const SingleProcess<DesignationType>& rhsSingleProcess) {
  SingleProcess<DesignationType> processTemp(rhsSingleProcess);
  processTemp+=lhsFileInput;
  return processTemp;
}

template<class DesignationType> SingleProcess<DesignationType> operator+(const SingleProcess<DesignationType>& lhsSingleProcess, const SingleProcess<DesignationType>& rhsSingleProcess) {
  SingleProcess<DesignationType> processTemp(lhsSingleProcess);
  processTemp+=rhsSingleProcess;
  return processTemp;
}

/// Define the utility functions
template<class DesignationType> void SingleProcess<DesignationType>::remakeIndividualFileHistograms() {
  for ( vector<SingleFileInput>::iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    fileIter->remakeFileHist();
  }
}

template<class DesignationType> void SingleProcess<DesignationType>::fillNullIndividualFileHistograms() {
  for ( vector<SingleFileInput>::iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    fileIter->makeFileHistIfNull();
  }
}

template<class DesignationType> void SingleProcess<DesignationType>::combineFileHistogramsToANewProcessHistogram() {
  if (processHist!=NULL) { delete processHist; }
  TH1F* tempHist = NULL;
  processHist = std::accumulate(allProcessFiles.begin(),allProcessFiles.end(),tempHist,addToHistAHistFromFile);
}

template<class DesignationType> void SingleProcess<DesignationType>::remakeProcessAndFileHistograms() {
  if ( isEmpty() ) {
    cout << "Error: No files present for the process" << endl;
    return;
  }
  remakeIndividualFileHistograms();
  combineFileHistogramsToANewProcessHistogram();
}

template<class DesignationType> void SingleProcess<DesignationType>::rescaleProcessAndFilesBy(double scaleFactor) {
  processHist->Scale(scaleFactor);
  for ( vector<SingleFileInput>::iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    fileIter->scaleBy(scaleFactor);
  }
}

template<class DesignationType> void SingleProcess<DesignationType>::setVariableInAllFilesTo(const AnalysisVariable& variableIn) {
  for ( vector<SingleFileInput>::iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    fileIter->setVariable(variableIn);
  }
}

template<class DesignationType> void SingleProcess<DesignationType>::printSummary() const {
  for ( vector<SingleFileInput>::const_iterator fileIter = allProcessFiles.begin(); fileIter != allProcessFiles.end(); ++fileIter ) {
    fileIter->printSummary();
    cout << endl;
  }
}

///////////////  Step 3 - construct a stack histogram from individual processes
template<class DesignationType> class SimulationStack {
/// The stackHist is automatically updated every time a new process/component to a process is added
/// There is nothing preventing the instances from having data histograms, but the class is designed for dealing with a stack of simulation processes
/// Is set up as a template class in order to processes with different designation type, which will usually be a string (e.g. "diboson") or an enum created for the analysis. The DesignationType should have != and < operators defined.
public:
  SimulationStack(StackPlotProperties<DesignationType> stackHistPropertiesIn): stackHistProperties(stackHistPropertiesIn) { 
    stackLgnd=NULL;
  }
  ///Need to define the copy constructor, assignment operator and destructor to make sure the histogram (and not the pointer) is properly copied/deleted.
  SimulationStack(const SimulationStack& stackIn) {
    stackHistProperties=stackIn.stackHistProperties;
    stackMap=stackIn.stackMap;
    stackLgnd = NULL;
    //cout << "x" << endl;
    remakeProcessHistHistToDrawAndLgnd();
    //cout << "y" << endl;
  }
  SimulationStack& operator=(const SimulationStack& rhs) {
    stackHistProperties=rhs.stackHistProperties;
    stackMap=rhs.stackMap;
    stackLgnd = NULL;
    remakeProcessHistHistToDrawAndLgnd();
    return *this;
  }
  ~SimulationStack() {
    ///Do not delete the histograms & legend after the code exits (keep them displayed on the canvas). 
    //clearStackHist();
    //if (stackLgnd!=NULL) { delete stackLgnd; }
  }

  ///Overloaded +=
  SimulationStack& operator+=(const SingleProcess<DesignationType>& );
  SimulationStack& operator+=(const SimulationStack& );

  ///Various tools
  void setPlotProperties(PlotProperties plotPropIn) { stackHistProperties.setPlotProperties(plotPropIn); }
  void setStackHistTitle(TString histTitleIn) { stackHistProperties.setPlotTitle(histTitleIn); }
  void fillNullIndividualFileHistogramsForAllProcesses();
  void remakeIndividualFileHistogramsForAllProcesses();
  void combineFileHistogramsToANewProcessHistogramForAllProcesses();

  void clearStackHist();
  void reinitializeLgnd();
  virtual void populateHistsToDrawAndLgnd();
  virtual void remakeProcessHistHistToDrawAndLgnd();
  void remakeIndividualFileProcessStackHistAndLgnd();
  TCanvas* drawStackHistOnCanvas(TString="");
  void placeCMSLumiFrame(TPad* pad, int iPeriod=3, int iPosX=10) { stackHistProperties.setVarsAndRunCMSLumi(pad,iPeriod,iPosX); }

  void rescaleProcessBy(DesignationType,double);
  void rescaleDesignatedProcessesBy(vector<DesignationType>,double);
  void rescaleStackBy(double);
  void setVariableInAllStackProcessesTo(const AnalysisVariable& );
  StackPlotProperties<DesignationType> getStackHistProperties() { return stackHistProperties; }
  map< DesignationType,SingleProcess<DesignationType> > getStackMap() { return stackMap; }
  map< DesignationType,SingleProcess<DesignationType> > getStackMap() const { return stackMap; }
  bool isEmpty() const { return stackMap.empty(); }
  double eventsInStackHist() const;
  void printStackSummary() const;

private:
  StackPlotProperties<DesignationType> stackHistProperties;
  map< DesignationType,SingleProcess<DesignationType> > stackMap;
  vector<TH1F*> stackHistComponents;
  TLegend* stackLgnd;
};

/// Overload the += and + operators
template<class DesignationType> SimulationStack<DesignationType>& SimulationStack<DesignationType>::operator+=(const SingleProcess<DesignationType>& processIn) {
  stackMap[processIn.getProcessDesignation()] += processIn;
  //remakeProcessHistHistToDrawAndLgnd();
  return *this;
}

template<class DesignationType> SimulationStack<DesignationType>& SimulationStack<DesignationType>::operator+=(const SimulationStack<DesignationType>& stackIn) {
  map<DesignationType,SingleProcess<DesignationType> > mapIn = stackIn.stackMap;
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = mapIn.begin(); mapIter != mapIn.end(); ++mapIter ) {
    stackMap[mapIter->first]+=(mapIter->second);
  }
  //remakeProcessHistHistToDrawAndLgnd();
  return *this;
}

template<class DesignationType> SimulationStack<DesignationType> operator+(const SimulationStack<DesignationType>& lhsSimulationStack, const SingleProcess<DesignationType>& rhsSingleProcess ) {
  SimulationStack<DesignationType> stackTemp(lhsSimulationStack);
  stackTemp+=rhsSingleProcess;
  //stackTemp.remakeProcessHistHistToDrawAndLgnd();
  return stackTemp;
}

template<class DesignationType> SimulationStack<DesignationType> operator+(const SingleProcess<DesignationType>& lhsSingleProcess, const SimulationStack<DesignationType>& rhsSimulationStack) {
  SimulationStack<DesignationType> stackTemp(rhsSimulationStack);
  stackTemp+=lhsSingleProcess;
  //stackTemp.remakeProcessHistHistToDrawAndLgnd();
  return stackTemp;
}

template<class DesignationType> SimulationStack<DesignationType> operator+(const SimulationStack<DesignationType>& lhsSimulationStack, const SimulationStack<DesignationType>& rhsSimulationStack) {
  SimulationStack<DesignationType> stackTemp(lhsSimulationStack);
  stackTemp+=rhsSimulationStack;
  //stackTemp.remakeProcessHistHistToDrawAndLgnd();
  return stackTemp;
}

/// Define the utility functions
template<class DesignationType> void SimulationStack<DesignationType>::clearStackHist() {
  if ( !stackHistComponents.empty() ) {
    for ( vector<TH1F*>::iterator histIter = stackHistComponents.begin(); histIter != stackHistComponents.end(); ++histIter ) {
      if ( (*histIter!=NULL) ) { delete *histIter; }
    }
    stackHistComponents.clear();
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::reinitializeLgnd() {
  //cout << "r1" << endl;
  //clearStackHist();
  //cout << "r2" << endl;
  //cout << "Reinitializing the legend using the base class version" << endl;
  if (stackLgnd!=NULL) { delete stackLgnd; }
  //cout << "r3" << endl;
  stackLgnd = new TLegend(stackHistProperties.getLgndX1(),stackHistProperties.getLgndY1(),stackHistProperties.getLgndX2(),stackHistProperties.getLgndY2()); //e.g. (0.6,0.6,0.90,0.9)
  stackLgnd->SetFillColor(stackHistProperties.getLgndFillColor());
  //cout << "r5" << endl;
}


template<class DesignationType> void SimulationStack<DesignationType>::fillNullIndividualFileHistogramsForAllProcesses() {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    (mapIter->second).fillNullIndividualFileHistograms();
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::remakeIndividualFileHistogramsForAllProcesses() {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    (mapIter->second).remakeIndividualFileHistograms();
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::combineFileHistogramsToANewProcessHistogramForAllProcesses() {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    (mapIter->second).combineFileHistogramsToANewProcessHistogram();
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::populateHistsToDrawAndLgnd() {
  /// Produce a stack histogram and corresponding legend with the appropriate properties, arranging the individual histograms in the desired order
  //cout << "produceStackHistAndLgnd:" << endl;
  //reinitializeStackHistAndLgnd();

  //cout << "B" << endl;
  vector<DesignationType> orderedDesignations = stackHistProperties.getDesignationsTopToBottom();
  for ( typename vector<DesignationType>::iterator orderedDesigIter = orderedDesignations.begin(); orderedDesigIter != orderedDesignations.end(); ++orderedDesigIter ) {
    //cout << "C" << endl;
    /// If the process is a part of the stackMap add it to the stackHist and Legend
    DesignationType currentDesignation=*orderedDesigIter;
    //cout << "process = " << procLabel[currentDesignation] << ", stackMap.count(currentDesignation) = " << stackMap.count(currentDesignation) << endl;
    if ( stackMap.count(currentDesignation) ) {
      //cout << "adding the process " << procLabel[currentDesignation] << endl;
      //stackMap[currentDesignation].printSummary();
      TH1F* currentProcessHist = stackMap[currentDesignation].getProcessHistPtr();
      if ( currentProcessHist!=NULL ) {
      // The histogram for each process is a sum of the the histogram for that process plus all below it. These can then drawn on top of each other showing the right process.
	for ( vector<TH1F*>::iterator histIter = stackHistComponents.begin(); histIter != stackHistComponents.end(); ++histIter ) {
	  (*histIter)->Add(currentProcessHist);
	}
	stackHistComponents.push_back(new TH1F(*currentProcessHist));
	TH1F* lastHist = stackHistComponents.back();
	// 	cout << "C2" << endl;
// 	Color_t cc = stackHistProperties.getProcessColor(currentDesignation);
// 	cout << "C2.5, Nentries = " << lastHist->GetEntries() << endl;
	lastHist->SetFillColor(stackHistProperties.getProcessColor(currentDesignation));
	stackLgnd->AddEntry(lastHist,stackHistProperties.getProcessName(currentDesignation),"f");
      }
    }
    //    cout << "D" << endl;
  }//for
}

template<class DesignationType> void SimulationStack<DesignationType>::remakeProcessHistHistToDrawAndLgnd() {
  if ( isEmpty() ) { return; }
  combineFileHistogramsToANewProcessHistogramForAllProcesses();
  clearStackHist();
  reinitializeLgnd();
  populateHistsToDrawAndLgnd();
}

template<class DesignationType> void SimulationStack<DesignationType>::remakeIndividualFileProcessStackHistAndLgnd() {
  remakeIndividualFileHistogramsForAllProcesses();
  remakeProcessHistHistToDrawAndLgnd();
}

template<class DesignationType> TCanvas* SimulationStack<DesignationType>::drawStackHistOnCanvas(TString cnvNameIn) {
  TCanvas* cnvTemp = stackHistProperties.makeCnv(cnvNameIn);
  cnvTemp->cd();
  cnvTemp->Clear();
  TH1F* firstHist;
  if ( !stackHistComponents.empty() ) {
    firstHist = stackHistComponents[0];
    firstHist = stackHistProperties.labelHistogram(firstHist);
    firstHist->Draw("hist");
    firstHist->SetMinimum(stackHistProperties.getPlotMinimum());
    if ( stackHistProperties.usePlotMaximum() ) {
      firstHist->SetMaximum(stackHistProperties.getPlotMaximum());
    }
    for ( vector<TH1F*>::iterator histIter = ++stackHistComponents.begin(); histIter != stackHistComponents.end(); ++histIter ) {
      TH1F* nextHist=*histIter;
      if ( (nextHist!=NULL) ) { 
	nextHist->Draw("hist:same");
      }
    }
    firstHist->Draw("axis same");
  } else {
    cout << "Warning: unable to draw stack histograms: the stackHistComponents vector is empty" << endl;
  }
  stackLgnd->Draw();
  cnvTemp->Update();
  return cnvTemp;
}

template<class DesignationType> void SimulationStack<DesignationType>::rescaleProcessBy(DesignationType processDesignation, double scaleFactor) {
  if ( stackMap.count(processDesignation) ) {
    stackMap[processDesignation].rescaleProcessAndFilesBy(scaleFactor);
  } else {
    cout << "Warning: unable to rescale the process. No matching designation found in stackMap" << endl;
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::rescaleDesignatedProcessesBy(vector<DesignationType> processDesignations, double scaleFactor) {
  for ( typename vector<DesignationType>::iterator procIter = processDesignations.begin(); procIter != processDesignations.end(); ++procIter ) {
    rescaleProcessBy(*procIter,scaleFactor);
  }
  remakeProcessHistHistToDrawAndLgnd();
}

template<class DesignationType> void SimulationStack<DesignationType>::rescaleStackBy(double scaleFactor) {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    DesignationType currentDesignation = mapIter->first;
    rescaleProcessBy(currentDesignation,scaleFactor);
  }
  remakeProcessHistHistToDrawAndLgnd();
}

template<class DesignationType> void SimulationStack<DesignationType>::setVariableInAllStackProcessesTo(const AnalysisVariable& variableIn) {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    (mapIter->second).setVariableInAllFilesTo(variableIn);
  }
}

template<class DesignationType> void SimulationStack<DesignationType>::printStackSummary() const {
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::const_iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    (mapIter->second).printSummary();
  }
}

template<class DesignationType> double SimulationStack<DesignationType>::eventsInStackHist() const {
  double nEvents = -1.0;
  if ( !stackHistComponents.empty() ) {
    nEvents = stackHistComponents[0]->Integral();
  } else {
    cout << "Warning: the stackHist is empty - unable to compute event count (returning -1.0)" << endl;
  }
  return nEvents;
}

///////////////  Step 4 - Combine Data with Stacked MC
template<class DesignationType> class DataAndMCStack: public SimulationStack<DesignationType> {
/// Combines Data and Monte Carlo
/// Is set up as a template class in order to processes with different designation type, which will usually be a string (e.g. "diboson") or an enum created for the analysis. The DesignationType should have != and < operators defined.
public:
  //DataAndMCStack(): SimulationStack<DesignationType>() { dataHist=NULL; }
  DataAndMCStack(SimulationStack<DesignationType>& simulationStackIn): SimulationStack<DesignationType>(simulationStackIn) { dataHist=NULL; }
  DataAndMCStack(SimulationStack<DesignationType>& simulationStackIn, SingleProcess<DesignationType>& dataProcessIn): SimulationStack<DesignationType>(simulationStackIn), dataProcess(dataProcessIn) { dataHist=NULL; }

  DataAndMCStack(const DataAndMCStack& stackIn): SimulationStack<DesignationType>(stackIn) {
    //cout << "r" << endl;
    dataProcess=stackIn.dataProcess;
    //cout << "rr" << endl;
    remakeAllDataHists();
    //cout << "rrr" << endl;
  }
  DataAndMCStack& operator=(const DataAndMCStack& rhs) {
    dataProcess=rhs.dataProcess;
    dataHist = NULL;
    //cout << "t" << endl;
    SimulationStack<DesignationType>::operator=(rhs);
    //cout << "tt" << endl;
    return *this;
  }
  ~DataAndMCStack() {
    ///Do not delete the histograms & legend after the code exits (keep them displayed on the canvas). 
    //clearStackHist();
    //if (stackLgnd!=NULL) { delete stackLgnd; }
  }


  ///Add Data, MC or another DataAndMCStack
  DataAndMCStack& addDataProcess(const SingleProcess<DesignationType>& );
  DataAndMCStack& addSimulationProcess(const SingleProcess<DesignationType>& );
  DataAndMCStack& addSimulationStack(const SimulationStack<DesignationType>& );
  DataAndMCStack& operator+=(const DataAndMCStack& );

  ///Various tools
  void remakeDataProcessHists() { dataProcess.remakeProcessAndFileHistograms();}
  void remakeAllDataHists();
  void populateHistsToDrawAndLgnd();
  void remakeProcessHistHistToDrawAndLgnd();

  TCanvas* drawDataAndMCOnCanvas(TString="");
  // TCanvas* drawDataMinusDesignatedMCsOnCanvas(vector<DesignationType>&,TString="",TString="Data-Bkg",TString="Signal");

  double getDataEvents() { return dataProcess.getProcessEvents(); }
  TH1F* getDataHistPtr() { return dataHist; }
  void rescaleMCtoData();
  void rescaleDesignatedProcessesToMatchMCToData(vector<DesignationType>);
  void setVariableInAllProcessesTo(const AnalysisVariable& );

  void printDataSummary() const;
  void printSummary() const;

private:
  SingleProcess<DesignationType> dataProcess;
  TH1F* dataHist;
};

/// Add Data, MC or another DataAndMCStack
template<class DesignationType> DataAndMCStack<DesignationType>& DataAndMCStack<DesignationType>::addDataProcess(const SingleProcess<DesignationType>& dataIn) {
  dataProcess+=dataIn;
  return *this;
}

template<class DesignationType> DataAndMCStack<DesignationType>& DataAndMCStack<DesignationType>::addSimulationProcess(const SingleProcess<DesignationType>& simulationIn) {
  SimulationStack<DesignationType>::operator+=(simulationIn);
  return *this;
}

template<class DesignationType> DataAndMCStack<DesignationType>& DataAndMCStack<DesignationType>::addSimulationStack(const SimulationStack<DesignationType>& simulationStackIn) {
  SimulationStack<DesignationType>::operator+=(simulationStackIn);
  return *this;
}

template<class DesignationType> DataAndMCStack<DesignationType>& DataAndMCStack<DesignationType>::operator+=(const DataAndMCStack<DesignationType>& dataAndMCIn) {
  cout << "Using DataAndMCStack +=" << endl;
  SimulationStack<DesignationType>::operator+=(dataAndMCIn);//add the simulation components
  dataProcess+=dataAndMCIn.dataProcess;
  return *this;
}

template<class DesignationType> DataAndMCStack<DesignationType> operator+(const DataAndMCStack<DesignationType>& lhsDataAndMCStack, const SimulationStack<DesignationType>& rhsSimulationStack ) {
  DataAndMCStack<DesignationType> dataAndMCStackTemp(lhsDataAndMCStack);
  dataAndMCStackTemp+=rhsSimulationStack;
  return dataAndMCStackTemp;
}

template<class DesignationType> DataAndMCStack<DesignationType> operator+(const SimulationStack<DesignationType>& lhsSimulationStack, const DataAndMCStack<DesignationType>& rhsDataAndMCStack) {
  SimulationStack<DesignationType> dataAndMCStackTemp(rhsDataAndMCStack);
  dataAndMCStackTemp+=lhsSimulationStack;
  return dataAndMCStackTemp;
}

template<class DesignationType> DataAndMCStack<DesignationType> operator+(const DataAndMCStack<DesignationType>& lhsDataAndMCStack, const DataAndMCStack<DesignationType>& rhsDataAndMCStack) {
  SimulationStack<DesignationType> dataAndMCStackTemp(lhsDataAndMCStack);
  dataAndMCStackTemp+=rhsDataAndMCStack;
  return dataAndMCStackTemp;
}

/// Define the utility functions
template<class DesignationType> void DataAndMCStack<DesignationType>::remakeAllDataHists() {
  if ( dataProcess.isEmpty() ) { return; }
  remakeDataProcessHists();
  dataHist = new TH1F(*dataProcess.getProcessHistPtr());
}


template<class DesignationType> void DataAndMCStack<DesignationType>::populateHistsToDrawAndLgnd() {
  cout << "Running DataAndMCStack version of populateHistsToDrawAndLgnd()" << endl;
  if ( dataHist!=NULL ) { delete dataHist; }
  dataHist = new TH1F(*dataProcess.getProcessHistPtr());
  dataHist->SetLineWidth(2);
  dataHist->SetMarkerStyle(8);
  SimulationStack<DesignationType>::stackLgnd->AddEntry(dataHist,"Data","l");
  SimulationStack<DesignationType>::populateHistsToDrawAndLgnd();
}

template<class DesignationType> void DataAndMCStack<DesignationType>::remakeProcessHistHistToDrawAndLgnd() {
  cout << "Running DataAndMCStack version of remakeProcessHistHistToDrawAndLgnd()" << endl;
  this->remakeAllDataHists();
  SimulationStack<DesignationType>::remakeProcessHistHistToDrawAndLgnd();
}

template<class DesignationType> TCanvas* DataAndMCStack<DesignationType>::drawDataAndMCOnCanvas(TString cnvNameIn) {
  TCanvas* cnvTemp = this->drawStackHistOnCanvas(cnvNameIn);
  cnvTemp->cd();
  dataHist->Draw("Ep:same");
  cnvTemp->Update();
  return cnvTemp;
}

// template<class DesignationType>TCanvas* DataAndMCStack<DesignationType>::drawDataMinusDesignatedMCsOnCanvas(vector<DesignationType>& processDesignationsToSubtract, TString cnvNameIn, TString dataMinusSubtractedLabel, TString remainingProcessLabel) {
//   TH1F* dataMinusSubtracted = new TH1F(*dataHist);
//   dataMinusSubtracted->SetLineWidth(2);
//   dataMinusSubtracted->SetMarkerStyle(8);
//   TLegend* subtractedLgnd = new TLegend(0.7,0.75,0.90,0.9);
//   subtractedLgnd->SetFillColor(this->stackHistProperties.getLgndFillColor());
//   subtractedLgnd->AddEntry(dataMinusSubtracted,dataMinusSubtractedLabel,"l");
//   TH1F* remaining;
//   bool firstRemainingProcess=true;
//   for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = (this->stackMap).begin(); mapIter != (this->stackMap).end(); ++mapIter ) {
//     DesignationType currentDesignation = mapIter->first;
//     if ( std::find(processDesignationsToSubtract.begin(), processDesignationsToSubtract.end(), currentDesignation) != processDesignationsToSubtract.end() ) {
//       dataMinusSubtracted->Add((mapIter->second).getProcessHistPtr(),-1.0);
//     } else {
//       if (firstRemainingProcess) {
// 	firstRemainingProcess=false;
// 	remaining = new TH1F(*(mapIter->second).getProcessHistPtr());
// 	remaining->SetFillColor(this->stackHistProperties.getProcessColor(currentDesignation));
// 	subtractedLgnd->AddEntry(remaining,remainingProcessLabel,"f");
//       } else {
// 	remaining->Add((mapIter->second).getProcessHistPtr());
//       }
//     }
//   }
//   if ( firstRemainingProcess ) {
//     cout << "Warning: no remaining processes after subtraction" << endl;
//     return NULL;
//   }

//   TCanvas* cnvTemp = this->stackHistProperties.makeCnv(cnvNameIn);
//   cnvTemp->cd();
//   dataMinusSubtracted->Draw("Ep");
//   this->stackHistProperties.labelHistogram(dataMinusSubtracted);
//   remaining->Draw("hist:same");
//   dataMinusSubtracted->Draw("Ep:same");
//   subtractedLgnd->Draw();
//   cnvTemp->Update();
//   return cnvTemp;
// }

template<class DesignationType> void DataAndMCStack<DesignationType>::rescaleMCtoData() {
  double dataEvts=getDataEvents();
  double mcEvts=this->eventsInStackHist();
  if ( mcEvts>0.001 ) {
    this->rescaleStackBy(dataEvts/mcEvts);
  } else {
    cout << "Warning: unable to rescale MC events to data, MC event cout <=0.001" << endl;
  }
}

template<class DesignationType> void DataAndMCStack<DesignationType>::rescaleDesignatedProcessesToMatchMCToData(vector<DesignationType> processDesignationsToScale) {
  double dataEvts=getDataEvents();
  double mcEvtsToScale=0, mcEvtsNotToScale=0; //Before scaling mcEvts=mcEvtsNotToScale+mcEvtsToScale. After scaling mcEvts=mcEvtsNotToScale+mcEvtsToScale*(dataEvts-mcEvtsNotToScale)/mcEvtsToScale=dataEvts;
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = this->stackMap.begin(); mapIter != this->stackMap.end(); ++mapIter ) {
    DesignationType currentDesignation = mapIter->first;
    if ( std::find(processDesignationsToScale.begin(), processDesignationsToScale.end(), currentDesignation) != processDesignationsToScale.end() ) {
      mcEvtsToScale+=(mapIter->second).getProcessEvents();
    } else {
      mcEvtsNotToScale+=(mapIter->second).getProcessEvents();
    }
  }
  if ( mcEvtsToScale>0.001 ) {
    rescaleDesignatedProcessesBy(processDesignationsToScale,(dataEvts-mcEvtsNotToScale)/mcEvtsToScale);
  } else {
    cout << "Warning: unable to rescale MC events to data, MC event cout <=0.001" << endl;
  }
}

template<class DesignationType> void DataAndMCStack<DesignationType>::setVariableInAllProcessesTo(const AnalysisVariable& variableIn) {
  this->setVariableInAllStackProcessesTo(variableIn);
  dataProcess.setVariableInAllFilesTo(variableIn);
}

template<class DesignationType> void DataAndMCStack<DesignationType>::printDataSummary() const {
  cout << "Data:" << endl;
  dataProcess.printSummary();
  cout << endl;
}

template<class DesignationType> void DataAndMCStack<DesignationType>::printSummary() const {
  printDataSummary();
  this->printStackSummary();
}

///////////////  Step 5 - Make a data-mc histogram and accound for additional errors
class SubtractedHist{
/// Most of the work is performed by the constructor
/// Histograms with full errors need to be constructed in the end due to the way fit error & additional errors factor in
public:
  SubtractedHist() {
    dataMinusSubtracted=NULL;
    remainingSimulation=NULL;
    fitErrHist=NULL;
    subtractedLgnd=NULL;
  }
  template<class DesignationType> SubtractedHist(DataAndMCStack<DesignationType>&, vector<DesignationType>&, double=1.0, TString="Data-Bkg", TString="Signal");
  ///Need to define the copy constructor, assignment operator and destructor to make sure the histograms/legend (and not the pointers) are properly copied/deleted.
  SubtractedHist(const SubtractedHist& subtractedHistIn) {
    if ( (subtractedHistIn.dataMinusSubtracted==NULL)||(subtractedHistIn.remainingSimulation==NULL)||(subtractedHistIn.fitErrHist==NULL)||(subtractedHistIn.subtractedLgnd==NULL) ) {
      dataMinusSubtracted=NULL;
      remainingSimulation=NULL;
      fitErrHist=NULL;
      subtractedLgnd=NULL;
    } else {
      dataMinusSubtracted=new TH1F(*subtractedHistIn.dataMinusSubtracted);
      remainingSimulation=new TH1F(*subtractedHistIn.remainingSimulation);
      fitErrHist=new TH1F(*subtractedHistIn.fitErrHist);
      subtractedLgnd=new TLegend(*subtractedHistIn.subtractedLgnd);
      cnvProperties=subtractedHistIn.cnvProperties;
    }
  }
  SubtractedHist& operator=(const SubtractedHist& rhs) {
    if ( (rhs.dataMinusSubtracted==NULL)||(rhs.remainingSimulation==NULL)||(rhs.fitErrHist==NULL)||(rhs.subtractedLgnd==NULL) ) {
      dataMinusSubtracted=NULL;
      remainingSimulation=NULL;
      fitErrHist=NULL;
      subtractedLgnd=NULL;
    } else {
      dataMinusSubtracted=new TH1F(*rhs.dataMinusSubtracted);
      remainingSimulation=new TH1F(*rhs.remainingSimulation);
      fitErrHist=new TH1F(*rhs.fitErrHist);
      subtractedLgnd=new TLegend(*rhs.subtractedLgnd);
      cnvProperties=rhs.cnvProperties;
    }
    return *this;
  }
  ~SubtractedHist() { 
    ///Keep the histograms in memory & displayed on the canvas
    // if ( dataMinusSubtracted!=NULL ) { delete dataMinusSubtracted; }
    // if ( remainingSimulation!=NULL ) { delete remainingSimulation; }
    // if ( fitErrHist!=NULL ) { delete fitErrHist; }
    // if ( subtractedLgnd!=NULL ) { delete subtractedLgnd; }
  }

  TH1F* makeDataPlusFullErrHist();
  TCanvas* drawSubtractedHistsWithTitle(TString="");
  void placeCMSLumiFrame(TPad* pad, int iPeriod=3, int iPosX=10) { cnvProperties.setVarsAndRunCMSLumi(pad,iPeriod,iPosX); }

  ///Add multiple SubtractedHists
  SubtractedHist& operator+=(const SubtractedHist& );

private:
  TH1F* dataMinusSubtracted;
  TH1F* remainingSimulation;
  TH1F* fitErrHist;
  TLegend* subtractedLgnd;
  CanvasProperties cnvProperties;
};

template<class DesignationType> SubtractedHist::SubtractedHist(DataAndMCStack<DesignationType>& dataAndMCStackIn, vector<DesignationType>& processDesignationsToSubtract, double fitErrFrac, TString dataMinusSubtractedLabel, TString remainingProcessLabel): dataMinusSubtracted(NULL), remainingSimulation(NULL), fitErrHist(NULL), subtractedLgnd(NULL)  { 
  if ( dataAndMCStackIn.getDataHistPtr()==NULL ) { return; }
  dataMinusSubtracted = new TH1F(*(dataAndMCStackIn.getDataHistPtr()));
  fitErrHist = new TH1F(*dataMinusSubtracted);
  StackPlotProperties<DesignationType> stackHistProperties = dataAndMCStackIn.getStackHistProperties();
  cnvProperties = CanvasProperties(stackHistProperties);
  dataMinusSubtracted->SetLineWidth(2);
  dataMinusSubtracted->SetMarkerStyle(8);
  subtractedLgnd = new TLegend(0.7,0.75,0.90,0.9);
  subtractedLgnd->SetFillColor(stackHistProperties.getLgndFillColor());
  subtractedLgnd->AddEntry(dataMinusSubtracted,dataMinusSubtractedLabel,"l");
  fitErrHist->SetLineStyle(2);
  fitErrHist->SetLineWidth(2);
  subtractedLgnd->AddEntry(fitErrHist,"Full Error","l");

  bool firstRemainingProcess=true;
  map< DesignationType,SingleProcess<DesignationType> > stackMap = dataAndMCStackIn.getStackMap();
  for ( typename map< DesignationType,SingleProcess<DesignationType> >::iterator mapIter = stackMap.begin(); mapIter != stackMap.end(); ++mapIter ) {
    DesignationType currentDesignation = mapIter->first;
    if ( std::find(processDesignationsToSubtract.begin(), processDesignationsToSubtract.end(), currentDesignation) != processDesignationsToSubtract.end() ) {
      dataMinusSubtracted->Add((mapIter->second).getProcessHistPtr(),-1.0);
    } else {
      if (firstRemainingProcess) {
	firstRemainingProcess=false;
	remainingSimulation = new TH1F(*(mapIter->second).getProcessHistPtr());
	remainingSimulation->SetFillColor(stackHistProperties.getProcessColor(currentDesignation));
	subtractedLgnd->AddEntry(remainingSimulation,remainingProcessLabel,"f");
      } else {
	remainingSimulation->Add((mapIter->second).getProcessHistPtr());
      }
    }
  }
  if ( firstRemainingProcess ) {
    cout << "Warning: no remaining processes after subtraction" << endl;
    return;
  }
  stackHistProperties.labelHistogram(dataMinusSubtracted);

  int NBins = fitErrHist->GetNbinsX();
  for (int n=1;n<(NBins+1);n++) {
    fitErrHist->SetBinContent(n,sqrt( (fitErrFrac*dataMinusSubtracted->GetBinContent(n))*(fitErrFrac*dataMinusSubtracted->GetBinContent(n)) ));
  }
}

TH1F* SubtractedHist::makeDataPlusFullErrHist()
//// Create the histogram which includes error due to fit as well as additional systematics
//// All histograms should have the same binning
{
  double systScaleErr = 0.070;//Including Lumi]
  if ( dataMinusSubtracted==NULL ) { return NULL; }
  TH1F* fullErrHist = new TH1F(*dataMinusSubtracted);
  int NBins = fullErrHist->GetNbinsX();
  double binVal = -1, binErrSquared = -1;

  for (int n=1;n<(NBins+1);n++) {
  ///Add the relevant errors for each bin in quadrature
    binVal=fullErrHist->GetBinContent(n);
    binErrSquared=(fullErrHist->GetBinError(n))*(fullErrHist->GetBinError(n));
    binErrSquared=binErrSquared+(systScaleErr*binVal)*(systScaleErr*binVal);
    binErrSquared=binErrSquared+(fitErrHist->GetBinContent(n))*(fitErrHist->GetBinContent(n));
    fullErrHist->SetBinError(n,sqrt(binErrSquared));
  }

  fullErrHist->SetLineStyle(2);
  fullErrHist->SetLineWidth(2);
  return fullErrHist;
}

TCanvas* SubtractedHist::drawSubtractedHistsWithTitle(TString histTitle) {
  TCanvas* cnvTemp = cnvProperties.makeCnv(histTitle);
  cnvTemp->cd();
  TH1F* fullErrHist = makeDataPlusFullErrHist();
  //dataMinusSubtracted->SetMaximum(1.5*dataMinusSubtracted->GetMaximum());
  //Set maximum to 115% of max bin value+error
  int binMax = fullErrHist->GetMaximumBin();
  dataMinusSubtracted->SetMaximum(1.15*(fullErrHist->GetBinContent(binMax)+fullErrHist->GetBinError(binMax)));
  //Draw:
  dataMinusSubtracted->Draw("Ep");
  dataMinusSubtracted->SetTitle(histTitle);
  remainingSimulation->Draw("hist:same");
  fullErrHist->Draw("Ep:same");
  dataMinusSubtracted->Draw("Ep:same");
  dataMinusSubtracted->Draw("axis same");
  subtractedLgnd->Draw();
  cnvTemp->Update();
  return cnvTemp;
}

///Add multiple SubtractedHists
SubtractedHist& SubtractedHist::operator+=(const SubtractedHist& subtractedHistIn) {
  if ( (subtractedHistIn.dataMinusSubtracted==NULL)||(subtractedHistIn.remainingSimulation==NULL)||(subtractedHistIn.fitErrHist==NULL)||(subtractedHistIn.subtractedLgnd==NULL) ) {
    cout << "Warning: unable to add the subtracted hist, at least one of the members is empty" << endl;
    return *this;
  } else {
    if ( (dataMinusSubtracted==NULL)||(remainingSimulation==NULL)||(fitErrHist==NULL)||(subtractedLgnd==NULL) ) {
      cout << "Note: one of the subtracted hist members is empty, addition resets them to be the same as the input members" << endl;
      dataMinusSubtracted=new TH1F(*subtractedHistIn.dataMinusSubtracted);
      remainingSimulation=new TH1F(*subtractedHistIn.remainingSimulation);
      fitErrHist=new TH1F(*subtractedHistIn.fitErrHist);
      subtractedLgnd=new TLegend(*subtractedHistIn.subtractedLgnd);
      cnvProperties=subtractedHistIn.cnvProperties;
    } else {
      /// Add the histograms (the legends & cnv properties should be the same)
      dataMinusSubtracted->Add(subtractedHistIn.dataMinusSubtracted);
      remainingSimulation->Add(subtractedHistIn.remainingSimulation);
      ///Add the fit errors in quadrature
      int NBins = fitErrHist->GetNbinsX();
      double binErrSqared;
      for (int n=1;n<(NBins+1);n++) {
	binErrSqared=(fitErrHist->GetBinContent(n))*(fitErrHist->GetBinContent(n));
	binErrSqared=binErrSqared+(subtractedHistIn.fitErrHist->GetBinContent(n))*(subtractedHistIn.fitErrHist->GetBinContent(n));
	fitErrHist->SetBinContent(n,sqrt(binErrSqared));
      }
    }
  }
  return *this;
}

SubtractedHist operator+(const SubtractedHist& lhsSubtractedHist, const SubtractedHist& rhsSubtractedHist) {
  SubtractedHist subtractedHistTemp(lhsSubtractedHist);
  subtractedHistTemp+=rhsSubtractedHist;
  return subtractedHistTemp;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////*********    Run Some Tests To Validate    *********//////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
void PlotTools () {
  ////UnitTests:
  gStyle->SetOptStat(0);

  /*
//   ///CMS Lumi Stuff
//   writeExtraText = true;       // if extra text
//   extraText  = "Preliminary";  // default extra text is "Preliminary"
//   lumi_8TeV  = "19 fb^{-1}"; // default is "19.7 fb^{-1}"
//   relPosX = 0.3; // Adjust the label postion

  ///Draw histograms and compare event counts to expectations for the WW file and the WW+WZ combination
  bool displayHistograms=true;
  bool displayAndSaveHistograms=false;
  if ( displayAndSaveHistograms ) { displayHistograms=true; }
  bool printSummaries=false;

  ///Input Files & Parameters
  double lumi=19300.0;
  TString inputFileDir = "/eos/uscms/store/user/lnujj/DibosonFitPostMoriond2013/";
  TString wwFilePath = inputFileDir + "FT_mu_WWtoLNuQQ_amcnlo_800k_WWptResumWts_CMSSW532_NominalWWpTwt.root";
  TString wzFilePath = inputFileDir + "RD_mu_WZtoLNuQQ_amcnlo_100k_CMSSW532.root";
  const int NVjetsFiles = 5;
  TString vjetsFilePaths[NVjetsFiles] = { inputFileDir+"RD_mu_W1Jets_CMSSW532.root", inputFileDir+"RD_mu_W2Jets_CMSSW532.root", inputFileDir+"RD_mu_W3Jets_CMSSW532.root", inputFileDir+"RD_mu_W4Jets_CMSSW532.root", inputFileDir+"RD_mu_ZpJ_CMSSW532.root"};
  double wjets_kfactor = 1.16;
  double vjetsCoeffs[NVjetsFiles] = { lumi*wjets_kfactor*5400.0/19871598.0, lumi*wjets_kfactor*1750.0/33004921.0, lumi*wjets_kfactor*519.0/15059503.0, lumi*wjets_kfactor*214.0/12842803.0, lumi*3503.71/30209426.0 }; 
  const int NTopFiles = 7;
  TString topFilePaths[NTopFiles] = { inputFileDir+"RD_mu_TTbar_CMSSW532.root", inputFileDir+"RD_mu_STopT_T_CMSSW532.root", inputFileDir+"RD_mu_STopT_Tbar_CMSSW532.root", inputFileDir+"RD_mu_STopS_T_CMSSW532.root", inputFileDir+"RD_mu_STopS_Tbar_CMSSW532.root", inputFileDir+"RD_mu_STopTW_T_CMSSW532.root", inputFileDir+"RD_mu_STopTW_Tbar_CMSSW532.root" };
  double topCoeffs[NTopFiles] = { lumi*225.197/6893735.0, lumi*55.531/3758221.0, lumi* 30.0042/1935066.0, lumi*3.89394/259960.0, lumi*1.75776/139974.0, lumi*11.1773/497657.0, lumi*11.1773/493458.0 };
  double crossX_WW_aMCNLO=57.52*0.444408;
  double crossX_WZ_aMCNLO=24.20*0.230996;
  double coeffWW = lumi*crossX_WW_aMCNLO/799899.0;
  double coeffWZ = lumi*crossX_WZ_aMCNLO/99791.0;
  TString cutString = "( (sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))>70.)&&(abs(JetPFCor_Eta[0]-JetPFCor_Eta[1])<1.5)&&(abs(JetPFCor_dphiMET[0])>0.4)&&(W_mt>30.)&&(JetPFCor_Pt[1]>35.)&&(JetPFCor_Pt[0]>40.)&&(JetPFCor_Pt[2]<30.)&&((abs(JetPFCor_Eta[0])>2.4)||(JetPFCor_Pt[0]<30.)||(JetPFCor_bDiscriminatorCSV[0]<0.244))&&((abs(JetPFCor_Eta[1])>2.4)||(JetPFCor_Pt[1]<30.)||(JetPFCor_bDiscriminatorCSV[1]<0.244))&&((abs(JetPFCor_Eta[2])>2.4)||(JetPFCor_Pt[2]<30.)||(JetPFCor_bDiscriminatorCSV[2]<0.244))&&((abs(JetPFCor_Eta[3])>2.4)||(JetPFCor_Pt[3]<30.)||(JetPFCor_bDiscriminatorCSV[3]<0.244))&&((abs(JetPFCor_Eta[4])>2.4)||(JetPFCor_Pt[4]<30.)||(JetPFCor_bDiscriminatorCSV[4]<0.244))&&((abs(JetPFCor_Eta[5])>2.4)||(JetPFCor_Pt[5]<30.)||(JetPFCor_bDiscriminatorCSV[5]<0.244))&&(W_pt<200.)&&(vbf_event==0)&&(event_met_pfmet>25)&&(abs(W_muon_eta)<2.1)&&(W_muon_pt>25.)&&(Mass2j_PFCor>48.000)&&(Mass2j_PFCor<160.000) )";
  TString eventWeightString = "effwt*puwt";
  TString dataFilePath = inputFileDir + "RD_WmunuJets_DataAll_GoldenJSON_19p3invfb.root";
  TString dataEventWeightString = "";
  AnalysisVariable mjj("Mass2j_PFCor",48,160,14);


  /// diboson tests
  double numWWExpected = 1803.8;
  double numWZExpected = 354.3;
  double numDataExpected = 126764;

  SingleFileInput ww(wwFilePath,mjj,coeffWW,cutString,eventWeightString);
  //ww.printSummary();
  ww.remakeFileHist();
  TH1F* wwHist = new TH1F(*ww.getFileHistPtr());
  cout << "wwHist->Integral()=" << wwHist->Integral() << " vs. expected " << numWWExpected << endl;

  if ( displayHistograms ) {
    TCanvas *cnvww = new TCanvas("cnvww","cnvww",10,10,900,600);
    cnvww->cd();
    cnvww->Clear();
    wwHist->Draw("ep");
    cnvww->Update();
  }
  cout << "Finished Testing WW" << endl;
  cout << endl;
  cout << endl;
  cout << endl;


  cout << "Testing wwPlusWz Process" << endl;
  cout << endl;
  SingleFileInput wz(wzFilePath,mjj,coeffWZ,cutString,eventWeightString);
  SingleProcess<TString> diboson("diboson"), diboson_ww("wwcomponent"), diboson_wz("wzcomponent"), diboson_Added("diboson");
  diboson+=ww;
  diboson+=wz;
  diboson_ww+=ww;
  diboson_wz+=wz;
  diboson_Added=diboson_ww+diboson_wz;

  cout << "WW and WZ added as files summary" << endl;
  if ( printSummaries ) { diboson.printSummary(); }
  cout << endl;

  cout << "WW and WZ added as processes summary" << endl;
  if ( printSummaries ) { diboson_Added.printSummary(); }
  cout << endl;

  diboson.remakeProcessAndFileHistograms();
  TH1F* dibosonHist = new TH1F(*diboson.getProcessHistPtr());
  diboson_Added.remakeProcessAndFileHistograms();
  TH1F* diboson_AddedHist = new TH1F(*diboson_Added.getProcessHistPtr());
  cout << "dibosonHist->Integral()=" << dibosonHist->Integral();
  cout << " vs. diboson_AddedHist->Integral()=" << diboson_AddedHist->Integral();
  cout << " vs. expected " << numWWExpected+numWZExpected << endl;
  cout << endl;

  if ( displayHistograms ) {
    TCanvas *cnvdiboson = new TCanvas("cnvdiboson","cnvdiboson",10,10,900,600);
    cnvdiboson->cd();
    cnvdiboson->Clear();
    dibosonHist->Draw("Ep");
    cnvdiboson->Update();

    TCanvas *cnvdiboson_Added = new TCanvas("cnvdiboson_Added","cnvdiboson_Added",10,10,900,600);
    cnvdiboson_Added->cd();
    cnvdiboson_Added->Clear();
    diboson_AddedHist->Draw("hist");
    cnvdiboson_Added->Update();
  }

  cout << "Finished The Diboson Tests" << endl;
  cout << endl;
  cout << endl;


  cout << "Starting Stack Tests" << endl;
  SingleProcess<TString> vjets("vjets");
  for (int i=0; i<NVjetsFiles; i++) {
    SingleFileInput vjetsFileTemp(vjetsFilePaths[i],mjj,vjetsCoeffs[i],cutString,eventWeightString);
    vjets+=vjetsFileTemp;
  }

  SingleProcess<TString> top("top");
  for (int i=0; i<NTopFiles; i++) {
    SingleFileInput topFileTemp(topFilePaths[i],mjj,topCoeffs[i],cutString,eventWeightString);
    top+=topFileTemp;
  }

  CanvasProperties canvas900x600("cnv900x600",10,10,900,600);
  PlotProperties mjjStackPlotProperties("StackPlot: DijetMass - #mu Channel",0.25,0.5,0.97,"M_{jj}",0.9,"Expected Event Count",1.2);

  vector< pair< TString,pair<TString,Color_t> > > orderedNameColorPairs;
  orderedNameColorPairs.clear();
  orderedNameColorPairs.push_back(  make_pair( "diboson", make_pair("WW+WZ",kAzure+8) )  );
  orderedNameColorPairs.push_back(  make_pair( "top", make_pair("Top",kGreen+2) )  );
  orderedNameColorPairs.push_back(  make_pair( "qcd", make_pair("QCD",kGray) )  );
  orderedNameColorPairs.push_back(  make_pair( "vjets", make_pair("VJets",kRed) )  );

  StackPlotProperties<TString> muonChanHistProperties(canvas900x600,mjjStackPlotProperties,make_pair(0.6,0.6),make_pair(0.90,0.90),kWhite,orderedNameColorPairs);

  SimulationStack<TString> muonChanStack(muonChanHistProperties);
  cout << "Adding to the stack..." << endl;
  muonChanStack+=diboson;
  muonChanStack+=vjets;
  muonChanStack+=top;
  SimulationStack<TString> muonChanStackpT(muonChanStack);

  cout << "Muon Channel Stack Summary" << endl;
  if ( printSummaries ) { muonChanStack.printStackSummary(); }
  cout << endl;
  cout << "Remake all the mjj histograms" << endl;
  muonChanStack.remakeIndividualFileProcessStackHistAndLgnd();

  cout << "Replace mjj with dijet pT, remake the histograms and draw on cms_lumi template" << endl;
  AnalysisVariable pTjj("sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",70,240,17);
  muonChanStackpT.setVariableInAllStackProcessesTo(pTjj);
  PlotProperties pTStackPlotProperties("StackPlot: Dijet p_{T} - #mu Channel",0.25,0.5,0.97,"dijet/V_{had} candidate p_{T}",1.0,"Expected Event Count",1.2);
  muonChanStackpT.setPlotProperties(pTStackPlotProperties);
  muonChanStackpT.remakeIndividualFileProcessStackHistAndLgnd();

  cout << muonChanStack.eventsInStackHist() << " events in the mjj stack plot" << endl;
  cout << muonChanStackpT.eventsInStackHist() << " events in the dijet pT stack plot" << endl;

  if ( displayHistograms ) {
    TCanvas *mjjcnv, *ptcnv;
    mjjcnv = muonChanStack.drawStackHistOnCanvas("mjjcnv");
    mjjcnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(mjjcnv,"./Test_mjjStackHist.png"); }
    muonChanStackpT.setStackHistTitle("");
    ptcnv = muonChanStackpT.drawStackHistOnCanvas("ptcnv");
    muonChanStackpT.placeCMSLumiFrame(ptcnv,2,0);
    ptcnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(ptcnv,"./Test_ptStackHistOnCMSLumi.png"); }
  }


  cout << "Finished Stack Tests" << endl;
  cout << endl;
  cout << endl;

  cout << "Starting Tests With Data" << endl;

  SingleFileInput dataFile(dataFilePath,mjj,1.0,cutString,dataEventWeightString);
  SingleProcess<TString> data("data");
  data+=dataFile;
  DataAndMCStack<TString> dataPlusMCMuonChan(muonChanStack,data);
  if ( printSummaries ) { dataPlusMCMuonChan.printSummary(); }
  dataPlusMCMuonChan.remakeIndividualFileProcessStackHistAndLgnd();
  cout << dataPlusMCMuonChan.getDataEvents() << " Data events vs " << numDataExpected  << " expected and " << dataPlusMCMuonChan.eventsInStackHist() << " MC events" << endl;

//   cout << "Rescaling MC to data..." << endl;
//   dataPlusMCMuonChan.rescaleMCtoData();
//   dataPlusMCMuonChan.remakeIndividualFileProcessStackHistAndLgnd();
//   cout << dataPlusMCMuonChan.getDataEvents() << " Data events vs " << numDataExpected  << " expected and " << dataPlusMCMuonChan.eventsInStackHist() << " MC events" << endl;

  cout << "Rescaling V+Jets and QCD MC to match the event count to data..." << endl;
  vector<TString> backgroundsToScale;
  backgroundsToScale.clear();
  backgroundsToScale.push_back("qcd");
  backgroundsToScale.push_back("vjets");
  dataPlusMCMuonChan.rescaleDesignatedProcessesToMatchMCToData(backgroundsToScale);
  dataPlusMCMuonChan.remakeIndividualFileProcessStackHistAndLgnd();
  cout << dataPlusMCMuonChan.getDataEvents() << " Data events vs " << numDataExpected  << " expected and " << dataPlusMCMuonChan.eventsInStackHist() << " MC events" << endl;

  if ( displayHistograms ) {
    TCanvas *datavsmccnv, *dataminusnondibosoncnv;
    dataPlusMCMuonChan.setStackHistTitle("");
    datavsmccnv = dataPlusMCMuonChan.drawDataAndMCOnCanvas("datavsmccnv");
    dataPlusMCMuonChan.placeCMSLumiFrame(datavsmccnv,2,0);
    datavsmccnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(datavsmccnv,"./Test_mjjRescaledVJetsVsData.png"); }
    vector<TString> backgroundList;
    backgroundList.clear();
    //backgroundList.push_back("diboson");
    backgroundList.push_back("top");
    backgroundList.push_back("qcd");
    backgroundList.push_back("vjets");
    dataPlusMCMuonChan.setStackHistTitle("Data-Backgrounds: DijetMass - #mu Channel");
    dataminusnondibosoncnv = dataPlusMCMuonChan.drawDataMinusDesignatedMCsOnCanvas(backgroundList,"dataminusnondibosoncnv");
    dataminusnondibosoncnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(dataminusnondibosoncnv,"./Test_mjjDataMinusBackgroundsWithOnlyVJetsRescaled.png"); }
  }

  cout << "Finished Tests With Data" << endl;
  cout << endl;
  cout << endl;


  cout << "Finished All Tests" << endl;
*/

}

void RunTestsOnAnalysisInputs() {
  /// Options:
  bool displayHistograms=true;
  bool displayAndSaveHistograms=true;
  if ( displayAndSaveHistograms ) { displayHistograms=true; }
  bool printSummaries=false;
  initializeCombinedInputs();
  gStyle->SetOptStat(0);

  /// 1. Compare the expected numbers vs obtained those obtained by directly applying the cuts to the input files
  AnalysisVariable defaultVarAntiBtag(defaultVarName,defaultVarMin,defaultVarMax,defaultVarNBins);
  AnalysisVariable defaultVarBtag(defaultVarName,defaultVarBtagMin,defaultVarMax,defaultVarBtagNBins);
  SingleProcess<processes> processCollection[nResolvedChannels][nProcesses];
  map<processes,vector< pair<TString,double> > > currentChannel;
  vector< pair<TString,double> > currentProcess;
  bool emptyProcess=false;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    AnalysisVariable defaultVar = defaultVarAntiBtag;
    if ( (chan==mu_btagged)||(chan==el_btagged) ) { defaultVar = defaultVarBtag; }
    cout << "Running comparisons for " << chanTitleLabel[chan] << " channel" << endl;
    currentChannel = simFilesAndCoeffs[chan];
    for (int prc=0; prc<nProcesses;prc++) { 
      processes proc = (processes)prc;
      if ( currentChannel.count(proc) ) {
	emptyProcess=false;
	currentProcess = currentChannel[proc];
	for ( vector< pair<TString,double> >::iterator fileIter = currentProcess.begin(); fileIter != currentProcess.end(); ++fileIter ) {
	  if ( proc==qcd ) {
	    processCollection[chan][proc]+=SingleFileInput(fileIter->first,defaultVar,fileIter->second,qcdCutString[chan],dataEventWeightString);
	  } else {
	    processCollection[chan][proc]+=SingleFileInput(fileIter->first,defaultVar,fileIter->second,defaultCutString[chan],mcEventWeightString);
	  }
	}
	processCollection[chan][proc].setDesignation(proc);
	processCollection[chan][proc].remakeProcessAndFileHistograms();
      } else {
	emptyProcess=true;
      }
      cout << procLabel[proc] << " :" << endl;
      if ( emptyProcess ) {
	cout << "No files contributing to the process vs. " << expectedProcessYields[chan][proc] << " events expected" << endl;
      } else {
	cout << processCollection[chan][proc].getProcessEvents() << " events expected based on contributing files vs. " << expectedProcessYields[chan][proc] << " events expected" << endl;
      }
    }
    cout << endl;
  }

  cout << endl;
  cout << endl;

  /// 2. Reweigh to fit results, combine the individual processes and compare to data - small inconsistencies are possible since the Extended Maximum Likelihood Fit does not require an exact match.
  cout << "Reweighting to the fit results ..." << endl;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    for (int prc=0; prc<nProcesses;prc++) { 
      processes proc = (processes)prc;
      if ( !processCollection[chan][proc].isEmpty() ) { processCollection[chan][proc].rescaleProcessAndFilesBy(fitFracVal[chan][proc]); }
    }
  }
  cout << "Filling the data processes ..." << endl;
  vector < SingleProcess<processes> > dataByChannel(nResolvedChannels, SingleProcess<processes>(data));
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    AnalysisVariable defaultVar = defaultVarAntiBtag;
    if ( (chan==mu_btagged)||(chan==el_btagged) ) { defaultVar = defaultVarBtag; }
    dataByChannel[chan]+=SingleFileInput(dataFiles[chan],defaultVar,1.0,defaultCutString[chan],dataEventWeightString);
  }

  if ( printSummaries ) {
    cout << "Printing Data Summaries" << endl;
    for (int ch=0; ch<nResolvedChannels;ch++) {
      channels chan = (channels)ch;
      dataByChannel[chan].printSummary();
    }
  }

  cout << "Combining data and MC for each channel ..." << endl;
  CanvasProperties canvas900x600("cnv900x600",10,10,900,600);
  PlotProperties mjjStackPlotProperties("StackPlot: DijetMass - Insert Channel Title",0.25,0.5,0.97,"M_{jj}",0.9,"Event Count",1.2);
  StackPlotProperties<processes> defaultChanHistProperties(canvas900x600,mjjStackPlotProperties,make_pair(0.7,0.6),make_pair(0.90,0.90),legendColor,orderedNameColorPairs);
  SimulationStack<processes> emptyChanStack(defaultChanHistProperties);
  vector < DataAndMCStack<processes> > dataPlusMC(nResolvedChannels,DataAndMCStack<processes>(emptyChanStack));
  TCanvas* datavsmccnv[nResolvedChannels];
  for (int ch=0; ch<nResolvedChannels;ch++) {
    //cout << "ch=" << ch << endl;
    channels chan = (channels)ch;
    cout << chanTitleLabel[chan] << " channel" << endl;
    dataPlusMC[chan].setStackHistTitle("StackPlot: DijetMass - " + chanTitleLabel[chan]);
    dataPlusMC[chan].addDataProcess(dataByChannel[chan]);
    //cout << "a, ch=" << ch << endl;
    for (int prc=0; prc<nProcesses;prc++) { 
      //cout << "prc=" << prc << endl;
      processes proc = (processes)prc;
      if ( !processCollection[chan][proc].isEmpty() ) { 
	//cout << "a, prc=" << prc << endl;
	dataPlusMC[chan].addSimulationProcess(processCollection[chan][proc]); 
	//cout << "b, prc=" << prc << endl;
      }
    }
    if ( printSummaries ) { dataPlusMC[chan].printSummary(); }
    dataPlusMC[chan].remakeIndividualFileProcessStackHistAndLgnd();
    cout << dataPlusMC[chan].getDataEvents() << " Data events vs " << expectedDataYield[chan] << " expected and " << dataPlusMC[chan].eventsInStackHist() << " MC events reweighted to fit results" << endl;

    if ( displayHistograms ) {
      datavsmccnv[chan] = dataPlusMC[chan].drawDataAndMCOnCanvas("datavsmccnv"+chanPlotLabel[chan]);
      datavsmccnv[chan]->Update();
      if ( displayAndSaveHistograms ) { saveCanvasAs(datavsmccnv[chan],"./Test_mjjScaledToFit" + chanPlotLabel[chan] + ".png"); }
    }
  }
  cout << endl;
  cout << endl;

  /// 3. Combine the four channels to one
  cout << "Resetting the variable to pT and combining resolved channels ..." << endl;
  AnalysisVariable pTjj("sqrt(JetPFCor_Pt[0]**2+JetPFCor_Pt[1]**2+2*JetPFCor_Pt[0]*JetPFCor_Pt[1]*cos(JetPFCor_Phi[0]-JetPFCor_Phi[1]))",70,240,17);
  PlotProperties pTStackPlotProperties("StackPlot: Dijet p_{T} - All Resolved Channels",0.25,0.5,0.97,"dijet/V_{had} candidate p_{T}",1.0,"Expected Event Count",1.2);
  DataAndMCStack<processes> resolvedDataPlusMC(emptyChanStack);
  //resolvedDataPlusMC.setStackHistTitle("StackPlot: DijetPt - Resolved Channels");
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    dataPlusMC[chan].setVariableInAllProcessesTo(pTjj);
    dataPlusMC[chan].setPlotProperties(pTStackPlotProperties);
    dataPlusMC[chan].remakeIndividualFileProcessStackHistAndLgnd();
    resolvedDataPlusMC+=dataPlusMC[chan];
  }
  
//   resolvedDataPlusMC.setVariableInAllProcessesTo(pTjj);
//   resolvedDataPlusMC.setPlotProperties(pTStackPlotProperties);
  resolvedDataPlusMC.remakeIndividualFileProcessStackHistAndLgnd();

  cout << resolvedDataPlusMC.getDataEvents() << " total data events vs " << resolvedDataPlusMC.eventsInStackHist() << " total MC events reweighted to fit results" << endl;
  if ( displayHistograms ) {
    TCanvas *datavsmctotcnv;
    datavsmctotcnv = resolvedDataPlusMC.drawDataAndMCOnCanvas("datavsmctotcnv");
    datavsmctotcnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(datavsmctotcnv,"./Test_ptScaledToFit_TotalResolved.png"); }
    // resolvedDataPlusMC.setStackHistTitle("");
    // dataminusbkgtotcnv = resolvedDataPlusMC.drawDataMinusDesignatedMCsOnCanvas(backgroundList,"dataminusbkgtotcnv");
    // resolvedDataPlusMC.placeCMSLumiFrame(dataminusbkgtotcnv,2,0);
    // dataminusbkgtotcnv->Update();
    // if ( displayAndSaveHistograms ) { saveCanvasAs(dataminusbkgtotcnv,"./Test_ptDataMinusBackgroundsScaledToFit_TotalResolved.png"); }
  }
  cout << endl;
  cout << endl;

  /// 4. Make the subtracted histograms for each channel and combine
  cout << "Making the subtracted histograms ..." << endl;
  vector<processes> backgroundList;
  backgroundList.clear();
  backgroundList.push_back(vjets);
  backgroundList.push_back(top);
  backgroundList.push_back(whbb);
  backgroundList.push_back(qcd);

  SubtractedHist dataMinusBkg[nResolvedChannels], dataMinusBkgTotal;
  for (int ch=0; ch<nResolvedChannels;ch++) {
    channels chan = (channels)ch;
    dataMinusBkg[chan] = SubtractedHist(dataPlusMC[chan],backgroundList,fitFracErr[chan][diboson]);
    dataMinusBkgTotal+=dataMinusBkg[chan];
  }

  if ( displayHistograms ) {
    TCanvas *dataminusbkgcnv[nResolvedChannels], *dataminusbkgtotcnv;
    for (int ch=0; ch<nResolvedChannels;ch++) {
      channels chan = (channels)ch;
      dataPlusMC[chan].remakeIndividualFileProcessStackHistAndLgnd();
      dataminusbkgcnv[chan]=dataMinusBkg[chan].drawSubtractedHistsWithTitle("Subtracted Data vs Diboson MC : "+chanTitleLabel[chan]);
      if ( displayAndSaveHistograms ) { saveCanvasAs(dataminusbkgcnv[chan],"./Test_ptDataMinusBackgroundsScaledToFit"+chanPlotLabel[chan]+".png"); }
    }

    dataminusbkgtotcnv=dataMinusBkgTotal.drawSubtractedHistsWithTitle("");
    dataMinusBkgTotal.placeCMSLumiFrame(dataminusbkgtotcnv,2,0);
    dataminusbkgtotcnv->Update();
    if ( displayAndSaveHistograms ) { saveCanvasAs(dataminusbkgtotcnv,"./Test_ptDataMinusBackgroundsScaledToFit_TotalResolved.png"); }
  }

}

