void plotsqgl_ov(){

  gROOT->Reset();
  TCanvas *c1 = new TCanvas("c1","c1",800,800); 

// TFile f1 ("WWlvjj_TMVAoutput_final_newQGL/m300/TMVA_300_nJ2_el.root");
// TFile f1 ("WWlvjj_TMVAoutput_final_newQGL/m500/TMVA_500_nJ2_el.root");
 TFile f1 ("WWlvjj_TMVAoutput_final_newQGL/m500/TMVA_500_nJ2_mu.root");

// TFile f1 ("WWlvjj_TMVAoutput_final_newQGL/m300/TMVA_300_nJ2_mu.root");

	f1.cd("Method_Likelihood/Likelihood");
       TH1F *Train_S = (TH1F*)MVA_Likelihood_Train_S->Clone("Train_S");
	Train_S->SetLineColor(kBlack);
	Train_S->SetLineWidth(2);
        Train_S->SetStats(kFALSE);
	Train_S->Rebin(4);

      TH1F *Train_B = (TH1F*)MVA_Likelihood_Train_B->Clone("Train_B");
        Train_B->SetLineColor(kRed);
        Train_B->SetLineWidth(2);
        Train_B->SetStats(kFALSE);
        Train_B->Rebin(4);

       TH1F *Test_S = (TH1F*)MVA_Likelihood_S->Clone("Test_S");
        Test_S->SetLineColor(kBlue);
        Test_S->SetLineWidth(2);
        Test_S->SetStats(kFALSE);
        Test_S->Rebin(4);

      TH1F *Test_B = (TH1F*)MVA_Likelihood_B->Clone("Test_B");
        Test_B->SetLineColor(kGreen);
        Test_B->SetLineWidth(2);
        Test_B->SetStats(kFALSE);
        Test_B->Rebin(4);



//TFile f2 ("oldQG/WWlvjj_TMVAoutput_final_w/m300/TMVA_300_nJ2_el.root");
//TFile f2 ("oldQG/WWlvjj_TMVAoutput_final_w/m500/TMVA_500_nJ2_el.root");
//TFile f2 ("oldQG/WWlvjj_TMVAoutput_final_w/m500/TMVA_500_nJ2_mu.root");
//TFile f2 ("oldQG/WWlvjj_TMVAoutput_final_w/m300/TMVA_300_nJ2_mu.root");

//TFile f2 ("WWlvjj_TMVAoutput_withoutQGL/m500/TMVA_500_nJ2_mu.root");

//TFile f2 ("WWlvjj_TMVAoutput_withoutQGL/m500/TMVA_500_nJ2_el.root");

//TFile f2 ("WWlvjj_TMVAoutput_withoutQGL/m300/TMVA_300_nJ2_el.root");
//TFile f2 ("WWlvjj_TMVAoutput_withoutQGL/m300/TMVA_300_nJ2_mu.root");


/*        f2.cd("Method_Likelihood/Likelihood");
       TH1F *oldqg = (TH1F*)MVA_Likelihood_rejBvsS->Clone("oldqg");
        oldqg->SetLineColor(kGreen);
        oldqg->SetLineWidth(2);
        oldqg->SetStats(kFALSE);
*/
        Train_S->Draw();
        Train_B->Draw("sames");
        Test_S->Draw("sames");
	Test_B->Draw("sames");

	gPad->Modified(); gPad->Update();
/* TLegend *tleg4_1 = new TLegend(0.2444481,0.665852,0.5550148,0.8811685,NULL,"brNDC");
tleg4_1 ->SetTextSize(0.043);
tleg4_1->SetBorderSize(1);
tleg4_1->SetLineColor(1);
tleg4_1->SetLineStyle(1);
tleg4_1->SetLineWidth(2);
tleg4_1->SetFillColor(19);
tleg4_1->SetFillStyle(1001);
  tleg4_1->AddEntry(Train_S," New QG");
  tleg4_1->AddEntry(oldqg," Old QG");
   tleg4_1->Draw();
*/
//******************************
// c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/QGL/ROC/W4Jets_newaNormalizedtoArea_QGl_dis[0]30.png");
// c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/QGL/9dec/RECO/W4Jets_newaNormalizedtoArea_QGl_dis_nochs_anti_btag_[0]30.png");
// c1->SaveAs("HiggsMass300_iMVA_Likelihood_rejBvsS.png");
// c1->SaveAs("el_500_iMVA_Likelihood_rejBvsS.png");
// c1->SaveAs("mu_500_iMVA_Likelihood_rejBvsS.png");
// c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/ggh15Dec2013_Including_newQGl/Compare/mu_300_MVA_Likelihood_rejBvsS.png");

// c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/ggh15Dec2013_Including_newQGl/Compare/el_300_MVA_Likelihood_rejBvsS.png");

 c1->SaveAs("/afs/fnal.gov/files/home/room1/ajay/public_html/pub_ac/plots/ggh15Dec2013_Including_newQGl/Compare/mu_500_MVA_Likelihood_TestTrain.png");


}

