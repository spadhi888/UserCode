{
  gStyle->SetOptStat(0);

  TFile *fLM0 = new TFile("baslineeffunmatched.root");
  //  TFile *fLM0 = new TFile("baslineeffmatched.root");
  TH2F* nom = LM0x_hm12m0_all;
  TH2F* sig1 = LM0x_hm12m01sig_all;
  TH2F* sig2 = LM0x_hm12m02sig_all;
  TH2F* sig3 = LM0x_hm12m03sig_all;
  TH2F* sig4 = LM0x_hm12m04sig_all;


  nom->SetDirectory(0);
  sig1->SetDirectory(0);
  sig2->SetDirectory(0);
  sig3->SetDirectory(0);
  sig4->SetDirectory(0);
  fLM0->Close();


  TFile *fLM1 = new TFile("baslineeffmatched.root");
  TH2F* matched_nom = LM0x_hm12m0_all;
  TH2F* matched_sig1 = LM0x_hm12m01sig_all;
  TH2F* matched_sig2 = LM0x_hm12m02sig_all;
  TH2F* matched_sig3 = LM0x_hm12m03sig_all;
  TH2F* matched_sig4 = LM0x_hm12m04sig_all;


  matched_nom->SetDirectory(0);
  matched_sig1->SetDirectory(0);
  matched_sig2->SetDirectory(0);
  matched_sig3->SetDirectory(0);
  matched_sig4->SetDirectory(0);
  fLM1->Close();

  TFile *f = new TFile("stopeff.root","RECREATE");

  TH2F* hexbase = new TH2F("hexbase","hexbase", 60, 0.,1500., 40, 0. , 1000.);
  TH2F* hexsig1 = new TH2F("hexsig1","hexsig1", 60, 0.,1500., 40, 0. , 1000.);
  TH2F* hexsig2 = new TH2F("hexsig2","hexsig2", 60, 0.,1500., 40, 0. , 1000.);
  TH2F* hexsig3 = new TH2F("hexsig3","hexsig3", 60, 0.,1500., 40, 0. , 1000.);
  TH2F* hexsig4 = new TH2F("hexsig4","hexsig4", 60, 0.,1500., 40, 0. , 1000.);

  int binno = 0;
  
 for (int m0=0; m0 < 1500; m0=m0+25) {
     for (int m12=0; m12 < 1000; m12=m12+25) {
       float  outbase = nom->GetBinContent(nom->FindBin(m0,m12) );
       float  outsig1 = sig1->GetBinContent(sig1->FindBin(m0,m12) );
       float  outsig2 = sig2->GetBinContent(sig2->FindBin(m0,m12) );
       float  outsig3 = sig3->GetBinContent(sig3->FindBin(m0,m12) );
       float  outsig4 = sig4->GetBinContent(sig4->FindBin(m0,m12) );

       float  outmatched_base = matched_nom->GetBinContent(matched_nom->FindBin(m0,m12) );
       float  outmatched_sig1 = matched_sig1->GetBinContent(matched_sig1->FindBin(m0,m12) );
       float  outmatched_sig2 = matched_sig2->GetBinContent(matched_sig2->FindBin(m0,m12) );
       float  outmatched_sig3 = matched_sig3->GetBinContent(matched_sig3->FindBin(m0,m12) );
       float  outmatched_sig4 = matched_sig4->GetBinContent(matched_sig4->FindBin(m0,m12) );


       if (m0 == 1000 && nom->GetBinContent(nom->FindBin(m0+25,m12)) > 0 && nom->GetBinContent(nom->FindBin(m0-25,m12)) > 0 && nom->GetBinContent(nom->FindBin(m0,m12)) == 0) outbase = 0.5*(nom->GetBinContent(nom->FindBin(m0+25,m12)) + nom->GetBinContent(nom->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig1->GetBinContent(sig1->FindBin(m0+25,m12)) > 0 && sig1->GetBinContent(sig1->FindBin(m0-25,m12)) > 0 && sig1->GetBinContent(sig1->FindBin(m0,m12)) == 0) outsig1 = 0.5*(sig1->GetBinContent(sig1->FindBin(m0+25,m12)) + sig1->GetBinContent(sig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig2->GetBinContent(sig2->FindBin(m0+25,m12)) > 0 && sig2->GetBinContent(sig2->FindBin(m0-25,m12)) > 0 && sig2->GetBinContent(sig2->FindBin(m0,m12)) == 0) outsig2 = 0.5*(sig2->GetBinContent(sig2->FindBin(m0+25,m12)) + sig2->GetBinContent(sig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig3->GetBinContent(sig3->FindBin(m0+25,m12)) > 0 && sig3->GetBinContent(sig3->FindBin(m0-25,m12)) > 0 && sig3->GetBinContent(sig3->FindBin(m0,m12)) == 0) outsig3 = 0.5*(sig3->GetBinContent(sig3->FindBin(m0+25,m12)) + sig3->GetBinContent(sig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig4->GetBinContent(sig4->FindBin(m0+25,m12)) > 0 && sig4->GetBinContent(sig4->FindBin(m0-25,m12)) > 0 && sig4->GetBinContent(sig4->FindBin(m0,m12)) == 0) outsig4 = 0.5*(sig4->GetBinContent(sig4->FindBin(m0+25,m12)) + sig4->GetBinContent(sig4->FindBin(m0-25,m12)));



       if (m0 == 1000 && matched_nom->GetBinContent(matched_nom->FindBin(m0+25,m12)) > 0 && matched_nom->GetBinContent(matched_nom->FindBin(m0-25,m12)) > 0 && matched_nom->GetBinContent(matched_nom->FindBin(m0,m12)) == 0) outmatched_base = 0.5*(matched_nom->GetBinContent(matched_nom->FindBin(m0+25,m12)) + matched_nom->GetBinContent(matched_nom->FindBin(m0-25,m12)));

       if (m0 == 1000 && matched_sig1->GetBinContent(matched_sig1->FindBin(m0+25,m12)) > 0 && matched_sig1->GetBinContent(matched_sig1->FindBin(m0-25,m12)) > 0 && matched_sig1->GetBinContent(matched_sig1->FindBin(m0,m12)) == 0) outmatched_sig1 = 0.5*(matched_sig1->GetBinContent(matched_sig1->FindBin(m0+25,m12)) + matched_sig1->GetBinContent(matched_sig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && matched_sig2->GetBinContent(matched_sig2->FindBin(m0+25,m12)) > 0 && matched_sig2->GetBinContent(matched_sig2->FindBin(m0-25,m12)) > 0 && matched_sig2->GetBinContent(matched_sig2->FindBin(m0,m12)) == 0) outmatched_sig2 = 0.5*(matched_sig2->GetBinContent(matched_sig2->FindBin(m0+25,m12)) + matched_sig2->GetBinContent(matched_sig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && matched_sig3->GetBinContent(matched_sig3->FindBin(m0+25,m12)) > 0 && matched_sig3->GetBinContent(matched_sig3->FindBin(m0-25,m12)) > 0 && matched_sig3->GetBinContent(matched_sig3->FindBin(m0,m12)) == 0) outmatched_sig3 = 0.5*(matched_sig3->GetBinContent(matched_sig3->FindBin(m0+25,m12)) + matched_sig3->GetBinContent(matched_sig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && matched_sig4->GetBinContent(matched_sig4->FindBin(m0+25,m12)) > 0 && matched_sig4->GetBinContent(matched_sig4->FindBin(m0-25,m12)) > 0 && matched_sig4->GetBinContent(matched_sig4->FindBin(m0,m12)) == 0) outmatched_sig4 = 0.5*(matched_sig4->GetBinContent(matched_sig4->FindBin(m0+25,m12)) + matched_sig4->GetBinContent(matched_sig4->FindBin(m0-25,m12)));



       binno = nom->FindBin(m0,m12);
       if (outbase > 0) hexbase->SetBinContent(binno, outmatched_base/outbase);
       if (outsig1 > 0) hexsig1->SetBinContent(binno, outmatched_sig1/outsig1);
       if (outsig2 > 0) hexsig2->SetBinContent(binno, outmatched_sig2/outsig2);
       if (outsig3 > 0) hexsig3->SetBinContent(binno, outmatched_sig3/outsig3);
       if (outsig4 > 0) hexsig4->SetBinContent(binno, outmatched_sig4/outsig4);

       //       cout << m0 << " " << m12 << " " << outsig1 << endl;
     }
 }

 f->Write();

/*
  gStyle->SetOptStat(0);
  TLatex* mytext = new TLatex(500.,800.,"CMS, L_{int} = 0.98 fb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  // Plots 
  hexbase->Draw("COLZ");
  hexbase->SetXTitle("m_{#tilde{g}} (GeV)");
  hexbase->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext0 = new TLatex(500.,750.,"H_{T} > 80 GeV, MET > 30 GeV");
  mytext0->SetTextSize(0.03);
  mytext0->Draw("same");
  mytext->Draw("same");


  hexsig1->Draw("COLZ");
  hexsig1->SetXTitle("m_{#tilde{g}} (GeV)");
  hexsig1->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext1 = new TLatex(500.,750.,"H_{T} > 400 GeV, MET > 120 GeV");
  mytext1->SetTextSize(0.03);
  mytext1->Draw("same");
  mytext->Draw("same");


  hexsig2->Draw("COLZ");
  hexsig2->SetXTitle("m_{#tilde{g}} (GeV)");
  hexsig2->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext2 = new TLatex(500.,750.,"H_{T} > 400 GeV, MET > 50 GeV");
  mytext2->SetTextSize(0.03);
  mytext2->Draw("same");
  mytext->Draw("same");


  hexsig3->Draw("COLZ");
  hexsig3->SetXTitle("m_{#tilde{g}} (GeV)");
  hexsig3->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext3 = new TLatex(500.,750.,"H_{T} > 200 GeV, MET > 120 GeV");
  mytext3->SetTextSize(0.03);
  mytext3->Draw("same");
  mytext->Draw("same");


  hexsig4->Draw("COLZ");
  hexsig4->SetXTitle("m_{#tilde{g}} (GeV)");
  hexsig4->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext4 = new TLatex(500.,750.,"H_{T} > 80 GeV, MET > 100 GeV");
  mytext4->SetTextSize(0.03);
  mytext4->Draw("same");
  mytext->Draw("same");

*/




}
