{
  gStyle->SetOptStat(0);

// stopeff.root - baseline
// stopjeceff.root - JES
// stopfakeeffmatched.root - Signal contamination
// Add 2% for PDF

  TFile *fLM0 = new TFile("stopeff.root");
  TH2F* nom = hexbase;
  TH2F* sig1 = hexsig1;
  TH2F* sig2 = hexsig2;
  TH2F* sig3 = hexsig3;
  TH2F* sig4 = hexsig4;
  nom->SetDirectory(0);
  sig1->SetDirectory(0);
  sig2->SetDirectory(0);
  sig3->SetDirectory(0);
  sig4->SetDirectory(0);
  fLM0->Close();

  TFile *fLM1 = new TFile("stopjeceff.root");
  TH2F* jecnom = hexjecbase;
  TH2F* jecsig1 = hexjecsig1;
  TH2F* jecsig2 = hexjecsig2;
  TH2F* jecsig3 =  hexjecsig3;
  TH2F* jecsig4 = hexjecsig4;
  jecnom->SetDirectory(0);
  jecsig1->SetDirectory(0);
  jecsig2->SetDirectory(0);
  jecsig3->SetDirectory(0);
  jecsig4->SetDirectory(0);
  fLM1->Close();


  TFile *fLM2 = new TFile("stopfakeeffmatched.root");
  TH2F* fakenom = hexfakebase;
  TH2F* fakesig1 = hexfakesig1;
  TH2F* fakesig2 = hexfakesig2;
  TH2F* fakesig3 =  hexfakesig3;
  TH2F* fakesig4 = hexfakesig4;
  fakenom->SetDirectory(0);
  fakesig1->SetDirectory(0);
  fakesig2->SetDirectory(0);
  fakesig3->SetDirectory(0);
  fakesig4->SetDirectory(0);
  fLM2->Close();


  TFile *f = new TFile("stoptoteffmatched.root","RECREATE"); 

   TH2F* hextotbase = new TH2F("hextotbase","hextotbase", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hextotig1 = new TH2F("hextotsig1","hextotsig1", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hextotsig2 = new TH2F("hextotsig2","hextotsig2", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hextotsig3 = new TH2F("hextotsig3","hextotsig3", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hextotsig4 = new TH2F("hextotsig4","hextotsig4", 60, 0.,1500., 40, 0. , 1000.);

   int binno = 0;


 for (int m0=0; m0 < 1500; m0=m0+25) {
     for (int m12=0; m12 < 1000; m12=m12+25) {

       float  outbase = nom->GetBinContent(nom->FindBin(m0,m12) );
       float  outsig1 = sig1->GetBinContent(sig1->FindBin(m0,m12) );
       float  outsig2 = sig2->GetBinContent(sig2->FindBin(m0,m12) );
       float  outsig3 = sig3->GetBinContent(sig3->FindBin(m0,m12) );
       float  outsig4 = sig4->GetBinContent(sig4->FindBin(m0,m12) );

       float  outjecbase = jecnom->GetBinContent(jecnom->FindBin(m0,m12) );
       float  outjecsig1 = jecsig1->GetBinContent(jecsig1->FindBin(m0,m12) );
       float  outjecsig2 = jecsig2->GetBinContent(jecsig2->FindBin(m0,m12) );
       float  outjecsig3 = jecsig3->GetBinContent(jecsig3->FindBin(m0,m12) );
       float  outjecsig4 = jecsig4->GetBinContent(jecsig4->FindBin(m0,m12) );


       float  outfakebase = fakenom->GetBinContent(fakenom->FindBin(m0,m12) );
       float  outfakesig1 = fakesig1->GetBinContent(fakesig1->FindBin(m0,m12) );
       float  outfakesig2 = fakesig2->GetBinContent(fakesig2->FindBin(m0,m12) );
       float  outfakesig3 = fakesig3->GetBinContent(fakesig3->FindBin(m0,m12) );
       float  outfakesig4 = fakesig4->GetBinContent(fakesig4->FindBin(m0,m12) );

       //       if (outjecbase > 0) cout << m0 << " " << m12 << " " << outbase << " " << 1+(sqrt(outjecbase*outjecbase + outfakebase*outfakebase + 0.02*0.02)) << endl;
       // if (outjecsig1 > 0) cout << m0 << " " << m12 << " " << outsig1 << " " << 1+(sqrt(outjecsig1*outjecsig1 + outfakesig1*outfakesig1 + 0.02*0.02)) << endl;
       // if (outjecsig2 > 0) cout << m0 << " " << m12 << " " << outsig2 << " " << 1+(sqrt(outjecsig2*outjecsig2 + outfakesig2*outfakesig2 + 0.02*0.02)) << endl;
       // if (outjecsig3 > 0) cout << m0 << " " << m12 << " " << outsig3 << " " << 1+(sqrt(outjecsig3*outjecsig3 + outfakesig3*outfakesig3 + 0.02*0.02)) << endl;
       if (outjecsig4 > 0) cout << m0 << " " << m12 << " " << outsig4 << " " << 1+(sqrt(outjecsig4*outjecsig4 + outfakesig4*outfakesig4 + 0.02*0.02)) << endl;
       
       binno = nom->FindBin(m0,m12);

       hextotbase->SetBinContent(binno, sqrt(outjecbase*outjecbase + outfakebase*outfakebase + 0.02*0.02));
       hextotsig1->SetBinContent(binno, sqrt(outjecsig1*outjecsig1 + outfakesig1*outfakesig1 + 0.02*0.02));
       hextotsig2->SetBinContent(binno, sqrt(outjecsig2*outjecsig2 + outfakesig2*outfakesig2 + 0.02*0.02));
       hextotsig3->SetBinContent(binno, sqrt(outjecsig3*outjecsig3 + outfakesig3*outfakesig3 + 0.02*0.02));
       hextotsig4->SetBinContent(binno, sqrt(outjecsig4*outjecsig4 + outfakesig4*outfakesig4 + 0.02*0.02));

     }
 }
f->Write();
/*
gStyle->SetOptStat(0);
  TLatex* mytext = new TLatex(500.,800.,"CMS, L_{int} = 0.98 fb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  // Plots 
  hextotbase->Draw("COLZ");
  hextotbase->SetXTitle("m_{#tilde{g}} (GeV)");
  hextotbase->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext0 = new TLatex(500.,700.,"Total Unc., H_{T} > 80 GeV, MET > 30 GeV");
  mytext0->SetTextSize(0.03);
  mytext0->Draw("same");
  mytext->Draw("same");


  hextotsig1->Draw("COLZ");
  hextotsig1->SetXTitle("m_{#tilde{g}} (GeV)");
  hextotsig1->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext1 = new TLatex(500.,750.,"Total Unc., H_{T} > 400 GeV, MET > 120 GeV");
  mytext1->SetTextSize(0.03);
  mytext1->Draw("same");
  mytext->Draw("same");


  hextotsig2->Draw("COLZ");
  hextotsig2->SetXTitle("m_{#tilde{g}} (GeV)");
  hextotsig2->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext2 = new TLatex(500.,750.,"Total Unc., H_{T} > 400 GeV, MET > 50 GeV");
  mytext2->SetTextSize(0.03);
  mytext2->Draw("same");
  mytext->Draw("same");

  hextotsig3->Draw("COLZ");
  hextotsig3->SetXTitle("m_{#tilde{g}} (GeV)");
  hextotsig3->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext3 = new TLatex(500.,750.,"Total Unc., H_{T} > 200 GeV, MET > 120 GeV");
  mytext3->SetTextSize(0.03);
  mytext3->Draw("same");
  mytext->Draw("same");


  hextotsig4->Draw("COLZ");
  hextotsig4->SetXTitle("m_{#tilde{g}} (GeV)");
  hextotsig4->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext4 = new TLatex(500.,750.,"Total Unc., H_{T} > 80 GeV, MET > 100 GeV");
  mytext4->SetTextSize(0.03);
  mytext4->Draw("same");
  mytext->Draw("same");


*/

}
