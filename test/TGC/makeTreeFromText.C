
void makeTreeFromText() {

   makeTreeFromText("params-ww-kappa-lambda-cteq6l1-fitted-withCuts.list", 
		    "ww_root_kappa-lambda.root",
		    "lambda/F:dkappa:p0:p1:p2:p3:p4:p5:chi2");

   makeTreeFromText("params-ww-lambda-g1z-cteq6l1-fitted-withCuts.list", 
		    "ww_root_lambda-g1.root",
		    "lambda/F:dg1:p0:p1:p2:p3:p4:p5:chi2");

   makeTreeFromText("params-ww-kappa-g1z-cteq6l1-fitted-withCuts.list",
		    "ww_root_kappa-g1.root",
		    "dkappa/F:dg1:p0:p1:p2:p3:p4:p5:chi2");

   // Note: p4 missing for WZ functions

   makeTreeFromText("params-wz-kappa-lambda-cteq6l1-fitted-withCuts.list", 
		    "wz_root_kappa-lambda.root",
		    "lambda/F:dkappa:p0:p1:p2:p3:p5:p6:chi2");

   makeTreeFromText("params-wz-lambda-g1z-cteq6l1-fitted-withCuts.list", 
		    "wz_root_lambda-g1.root",
		    "lambda/F:dg1:p0:p1:p2:p3:p5:p6:chi2");

   makeTreeFromText("params-wz-kappa-g1z-cteq6l1-fitted-withCuts.list",
		    "wz_root_kappa-g1.root",
		    "dkappa/F:dg1:p0:p1:p2:p3:p5:p6:chi2");
}

void makeTreeFromText(const char* inputFile,
		      const char* outputFile,
		      const char *fmt) {
   TFile fout( outputFile, "recreate");
   TTree tree("tree","tree");

   tree.ReadFile( inputFile, fmt);

   tree.Write();
   // fout.Close();
}
