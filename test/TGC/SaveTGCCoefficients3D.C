#include <vector>
#include "atgcgrid.h"

void drawCoeffProfiles1pair(char *infile,
			    char *proc,  // WW or WZ
			    char *varx,  // what we want to call the variables
			    char *vary,
			    char *trvarx, // what the tree calls the variables
			    char *trvary,
			    int nptsx,float xmin,float xmax,
			    int nptsy,float ymin,float ymax,
			    std::vector<TProfile2D*>& outgraphs)
{
  TFile * fin = new TFile(infile);
  TTree* tr = (TTree*) fin->Get("tree"); 
  for (int i=0; i<=6; i++) {
    TString hname(Form("%s_p%d_%s_%s",proc,i,varx,vary));
    cout << "Drawing " << TString(hname)<<endl;
    TProfile2D* prof2d = new TProfile2D(hname,"",nptsx,xmin,xmax,nptsy,ymin,ymax);
    prof2d->SetTitle(Form("%s;%s;%s;p%d",hname.Data(),varx,vary,i));
    tr->Draw(Form("p%d:%s:%s>>%s",i,trvary,trvarx,hname.Data()),"","goff");
    outgraphs.push_back(prof2d);
  }
}

void SaveTGCCoefficients3D(){
   
   std::vector<TProfile2D*> outgraphs;

   TFile fout( "ATGC_shape_coefficients.root", "recreate");

   ////// ---- first process lambda:dkappa combination --------

   drawCoeffProfiles1pair("ww_root_kappa-lambda.root","ww",
			  "lambda","dkg",
			  "lambda","dkappa",
			  nLZ_pts,   LZ_MIN-( LZ_INC/2.),  LZ_MAX+( LZ_INC/2.),
			  ndKG_pts, dKG_MIN-(dKG_INC/2.), dKG_MAX+(dKG_INC/2.),
			  outgraphs);


   ////// ---- second process lambda:dg1 combination --------

   drawCoeffProfiles1pair("ww_root_lambda-g1.root","ww",
			  "lambda","dg1",
			  "lambda","dg1",
			  nLZ_pts,   LZ_MIN-( LZ_INC/2.),  LZ_MAX+( LZ_INC/2.),
			  ndg1_pts, dg1_MIN-(dg1_INC/2.), dg1_MAX+(dg1_INC/2.),
			  outgraphs);

   ////// ---- third process dkg:dg1 combination --------


   drawCoeffProfiles1pair("ww_root_kappa-g1.root", "ww",
			  "dkg","dg1",
			  "dkappa","dg1",
			  ndKG_pts, dKG_MIN-(dKG_INC/2.), dKG_MAX+(dKG_INC/2.),
			  ndg1_pts, dg1_MIN-(dg1_INC/2.), dg1_MAX+(dg1_INC/2.),
			  outgraphs);

   drawCoeffProfiles1pair("wz_root_kappa-lambda.root","wz",
			  "lambda","dkg",
			  "lambda","dkappa",
			  nLZ_pts,   LZ_MIN-( LZ_INC/2.),  LZ_MAX+( LZ_INC/2.),
			  ndKG_pts, dKG_MIN-(dKG_INC/2.), dKG_MAX+(dKG_INC/2.),
			  outgraphs);


   ////// ---- second process lambda:dg1 combination --------

   drawCoeffProfiles1pair("wz_root_lambda-g1.root","wz",
			  "lambda","dg1",
			  "lambda","dg1",
			  nLZ_pts,   LZ_MIN-( LZ_INC/2.),  LZ_MAX+( LZ_INC/2.),
			  ndg1_pts, dg1_MIN-(dg1_INC/2.), dg1_MAX+(dg1_INC/2.),
			  outgraphs);

   ////// ---- third process dkg:dg1 combination --------


   drawCoeffProfiles1pair("wz_root_kappa-g1.root", "wz",
			  "dkg","dg1",
			  "dkappa","dg1",
			  ndKG_pts, dKG_MIN-(dKG_INC/2.), dKG_MAX+(dKG_INC/2.),
			  ndg1_pts, dg1_MIN-(dg1_INC/2.), dg1_MAX+(dg1_INC/2.),
			  outgraphs);


   /////----- Let's save all histograms in the output file ----
   fout.cd();

   for (int i=0; i<outgraphs.size(); i++)
     outgraphs[i]->Write();

   fout.Close();
}
