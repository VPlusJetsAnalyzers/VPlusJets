void plotsqgl_eff(){

  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",800,800); 

 TFile f1 ("WWlvjj_TMVAoutput_21jan14/m300/TMVA_300_nJ2_mu.root");

	f1.cd("Method_Likelihood/Likelihood");
       TH1F *newqg = (TH1F*)MVA_Likelihood_rejBvsS->Clone("newqg");
	newqg->SetLineColor(kBlack);
	newqg->SetLineWidth(2);
	newqg->SetStats(kFALSE);



TFile f2 ("WWlvjj_TMVAoutput_21jan14_old/m300/TMVA_300_nJ2_mu.root");



        f2.cd("Method_Likelihood/Likelihood");
       TH1F *oldqg = (TH1F*)MVA_Likelihood_rejBvsS->Clone("oldqg");
        oldqg->SetLineColor(kGreen);
        oldqg->SetLineWidth(2);
        oldqg->SetStats(kFALSE);

	newqg->Draw();
        oldqg->Draw("sames");

	gPad->Modified(); gPad->Update();
/* TLegend *tleg4_1 = new TLegend(0.2444481,0.665852,0.5550148,0.8811685,NULL,"brNDC");
tleg4_1 ->SetTextSize(0.043);
tleg4_1->SetBorderSize(1);
tleg4_1->SetLineColor(1);
tleg4_1->SetLineStyle(1);
tleg4_1->SetLineWidth(2);
tleg4_1->SetFillColor(19);
tleg4_1->SetFillStyle(1001);
  tleg4_1->AddEntry(newqg," New QG");
  tleg4_1->AddEntry(oldqg," Old QG");
   tleg4_1->Draw();
*/
//******************************

 c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/21Jan_ggh/mu_compare_300_MVA_Likelihood_rejBvsS.png");


}

