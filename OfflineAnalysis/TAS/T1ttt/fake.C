{
  gStyle->SetOptStat(0);

  //  TFile *fLM0 = new TFile("baslineeffmatched.root");
  TFile *fLM0 = new TFile("baslineeffunmatched.root");
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

  TFile *fLM1 = new TFile("histo_yieldSFeffNOTmatched.root");
  //  TFile *fLM1 = new TFile("histo_yieldSFeffmatchedv1.root");

  TH2F* jecUpnom = LM0x_hm12m0_all;
  TH2F* jecUpsig1 = LM0x_hm12m01sig_all;
  TH2F* jecUpsig2 = LM0x_hm12m02sig_all;
  TH2F* jecUpsig3 = LM0x_hm12m03sig_all;
  TH2F* jecUpsig4 = LM0x_hm12m04sig_all;
  jecUpnom->SetDirectory(0);
  jecUpsig1->SetDirectory(0);
  jecUpsig2->SetDirectory(0);
  jecUpsig3->SetDirectory(0);
  jecUpsig4->SetDirectory(0);
  fLM1->Close();

  //  TFile *fLM2 = new TFile("histo_estimateDFeffmatchedv1.root");
  TFile *fLM2 = new TFile("histo_yieldDFeffNOTmatched.root");
  TH2F* jecDnnom = LM0x_hm12m0_all;
  TH2F* jecDnsig1 = LM0x_hm12m01sig_all;
  TH2F* jecDnsig2 = LM0x_hm12m02sig_all;
  TH2F* jecDnsig3 = LM0x_hm12m03sig_all;
  TH2F* jecDnsig4 = LM0x_hm12m04sig_all;
  jecDnnom->SetDirectory(0);
  jecDnsig1->SetDirectory(0);
  jecDnsig2->SetDirectory(0);
  jecDnsig3->SetDirectory(0);
  jecDnsig4->SetDirectory(0);
  fLM2->Close();


  TFile *f = new TFile("stopfakeeffmatched.root","RECREATE"); 

   TH2F* hexfakebase = new TH2F("hexfakebase","hexfakebase", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexfakeig1 = new TH2F("hexfakesig1","hexfakesig1", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexfakesig2 = new TH2F("hexfakesig2","hexfakesig2", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexfakesig3 = new TH2F("hexfakesig3","hexfakesig3", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexfakesig4 = new TH2F("hexfakesig4","hexfakesig4", 60, 0.,1500., 40, 0. , 1000.);

   int binno = 0;

 for (int m0=0; m0 < 1500; m0=m0+25) {
     for (int m12=0; m12 < 1000; m12=m12+25) {
       float  outbase = nom->GetBinContent(nom->FindBin(m0,m12) );
       float  outsig1 = sig1->GetBinContent(sig1->FindBin(m0,m12) );
       float  outsig2 = sig2->GetBinContent(sig2->FindBin(m0,m12) );
       float  outsig3 = sig3->GetBinContent(sig3->FindBin(m0,m12) );
       float  outsig4 = sig4->GetBinContent(sig4->FindBin(m0,m12) );

       if (m0 == 1000 && nom->GetBinContent(nom->FindBin(m0+25,m12)) > 0 && nom->GetBinContent(nom->FindBin(m0-25,m12)) > 0 && nom->GetBinContent(nom->FindBin(m0,m12)) == 0) outbase = 0.5*(nom->GetBinContent(nom->FindBin(m0+25,m12)) + nom->GetBinContent(nom->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig1->GetBinContent(sig1->FindBin(m0+25,m12)) > 0 && sig1->GetBinContent(sig1->FindBin(m0-25,m12)) > 0 && sig1->GetBinContent(sig1->FindBin(m0,m12)) == 0) outsig1 = 0.5*(sig1->GetBinContent(sig1->FindBin(m0+25,m12)) + sig1->GetBinContent(sig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig2->GetBinContent(sig2->FindBin(m0+25,m12)) > 0 && sig2->GetBinContent(sig2->FindBin(m0-25,m12)) > 0 && sig2->GetBinContent(sig2->FindBin(m0,m12)) == 0) outsig2 = 0.5*(sig2->GetBinContent(sig2->FindBin(m0+25,m12)) + sig2->GetBinContent(sig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig3->GetBinContent(sig3->FindBin(m0+25,m12)) > 0 && sig3->GetBinContent(sig3->FindBin(m0-25,m12)) > 0 && sig3->GetBinContent(sig3->FindBin(m0,m12)) == 0) outsig3 = 0.5*(sig3->GetBinContent(sig3->FindBin(m0+25,m12)) + sig3->GetBinContent(sig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && sig4->GetBinContent(sig4->FindBin(m0+25,m12)) > 0 && sig4->GetBinContent(sig4->FindBin(m0-25,m12)) > 0 && sig4->GetBinContent(sig4->FindBin(m0,m12)) == 0) outsig4 = 0.5*(sig4->GetBinContent(sig4->FindBin(m0+25,m12)) + sig4->GetBinContent(sig4->FindBin(m0-25,m12)));


       // JEC Up

       float  outjecUpbase = jecUpnom->GetBinContent(jecUpnom->FindBin(m0,m12) );
       float  outjecUpsig1 = jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0,m12) );
       float  outjecUpsig2 = jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0,m12) );
       float  outjecUpsig3 = jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0,m12) );
       float  outjecUpsig4 = jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0,m12) );

       if (m0 == 1000 && jecUpnom->GetBinContent(jecUpnom->FindBin(m0+25,m12)) > 0 && jecUpnom->GetBinContent(jecUpnom->FindBin(m0-25,m12)) > 0 && jecUpnom->GetBinContent(jecUpnom->FindBin(m0,m12)) == 0) outjecUpbase = 0.5*(jecUpnom->GetBinContent(jecUpnom->FindBin(m0+25,m12)) + jecUpnom->GetBinContent(jecUpnom->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0+25,m12)) > 0 && jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0-25,m12)) > 0 && jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0,m12)) == 0) outjecUpsig1 = 0.5*(jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0+25,m12)) + jecUpsig1->GetBinContent(jecUpsig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0+25,m12)) > 0 && jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0-25,m12)) > 0 && jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0,m12)) == 0) outjecUpsig2 = 0.5*(jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0+25,m12)) + jecUpsig2->GetBinContent(jecUpsig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0+25,m12)) > 0 && jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0-25,m12)) > 0 && jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0,m12)) == 0) outjecUpsig3 = 0.5*(jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0+25,m12)) + jecUpsig3->GetBinContent(jecUpsig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0+25,m12)) > 0 && jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0-25,m12)) > 0 && jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0,m12)) == 0) outjecUpsig4 = 0.5*(jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0+25,m12)) + jecUpsig4->GetBinContent(jecUpsig4->FindBin(m0-25,m12)));

       // JEC Dn

       float  outjecDnbase = jecDnnom->GetBinContent(jecDnnom->FindBin(m0,m12) );
       float  outjecDnsig1 = jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0,m12) );
       float  outjecDnsig2 = jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0,m12) );
       float  outjecDnsig3 = jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0,m12) );
       float  outjecDnsig4 = jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0,m12) );

       if (m0 == 1000 && jecDnnom->GetBinContent(jecDnnom->FindBin(m0+25,m12)) > 0 && jecDnnom->GetBinContent(jecDnnom->FindBin(m0-25,m12)) > 0 && jecDnnom->GetBinContent(jecDnnom->FindBin(m0,m12)) == 0) outjecDnbase = 0.5*(jecDnnom->GetBinContent(jecDnnom->FindBin(m0+25,m12)) + jecDnnom->GetBinContent(jecDnnom->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0+25,m12)) > 0 && jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0-25,m12)) > 0 && jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0,m12)) == 0) outjecDnsig1 = 0.5*(jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0+25,m12)) + jecDnsig1->GetBinContent(jecDnsig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0+25,m12)) > 0 && jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0-25,m12)) > 0 && jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0,m12)) == 0) outjecDnsig2 = 0.5*(jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0+25,m12)) + jecDnsig2->GetBinContent(jecDnsig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0+25,m12)) > 0 && jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0-25,m12)) > 0 && jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0,m12)) == 0) outjecDnsig3 = 0.5*(jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0+25,m12)) + jecDnsig3->GetBinContent(jecDnsig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0+25,m12)) > 0 && jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0-25,m12)) > 0 && jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0,m12)) == 0) outjecDnsig4 = 0.5*(jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0+25,m12)) + jecDnsig4->GetBinContent(jecDnsig4->FindBin(m0-25,m12)));

       float errbase = 0.;
       float errsig1 = 0.;
       float errsig2 = 0.;
       float errsig3 = 0.;
       float errsig4 = 0.;

       if (outbase > 0)  errbase = fabs(outjecUpbase - outjecDnbase)/outbase;
       if (outsig1 > 0)  errsig1 = fabs(outjecUpsig1 - outjecDnsig1)/outsig1;
       if (outsig2 > 0)  errsig2 = fabs(outjecUpsig2 - outjecDnsig2)/outsig2;
       if (outsig3 > 0)  errsig3 = fabs(outjecUpsig3 - outjecDnsig3)/outsig3;
       if (outsig4 > 0)  errsig4 = fabs(outjecUpsig4 - outjecDnsig4)/outsig4;


       //       if (outbase > 0)  errbase = fabs(outjecUpbase - outjecDnbase);
       //       if (outsig1 > 0)  errsig1 = fabs(outjecUpsig1 - outjecDnsig1);
       //       if (outsig2 > 0)  errsig2 = fabs(outjecUpsig2 - outjecDnsig2);
       //       if (outsig3 > 0)  errsig3 = fabs(outjecUpsig3 - outjecDnsig3);
       //       if (outsig4 > 0)  errsig4 = fabs(outjecUpsig4 - outjecDnsig4);

       binno = nom->FindBin(m0,m12);
       //       if (outbase > 0) cout << m0 << " " << m12 << " " << errbase << endl;

       hexfakebase->SetBinContent(binno, errbase);
       hexfakesig1->SetBinContent(binno, errsig1);
       hexfakesig2->SetBinContent(binno, errsig2);
       hexfakesig3->SetBinContent(binno, errsig3);
       hexfakesig4->SetBinContent(binno, errsig4);

       //       cout << m0 << " " << m12 << " " << outsig1 << endl;
     }
 }

f->Write();
/*

  gStyle->SetOptStat(0);
  TLatex* mytext = new TLatex(500.,800.,"CMS, L_{int} = 0.98 fb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  // Plots 
  hexfakebase->Draw("COLZ");
  hexfakebase->SetXTitle("m_{#tilde{g}} (GeV)");
  hexfakebase->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext0 = new TLatex(500.,700.,"Eff. Fakes, H_{T} > 80 GeV, MET > 30 GeV");
  mytext0->SetTextSize(0.03);
  mytext0->Draw("same");
  mytext->Draw("same");


  hexfakesig1->Draw("COLZ");
  hexfakesig1->SetXTitle("m_{#tilde{g}} (GeV)");
  hexfakesig1->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext1 = new TLatex(500.,750.,"Eff SB, H_{T} > 400 GeV, MET > 120 GeV");
  mytext1->SetTextSize(0.03);
  mytext1->Draw("same");
  mytext->Draw("same");


  hexfakesig2->Draw("COLZ");
  hexfakesig2->SetXTitle("m_{#tilde{g}} (GeV)");
  hexfakesig2->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext2 = new TLatex(500.,750.,"Eff. SB, H_{T} > 400 GeV, MET > 50 GeV");
  mytext2->SetTextSize(0.03);
  mytext2->Draw("same");
  mytext->Draw("same");


  hexfakesig3->Draw("COLZ");
  hexfakesig3->SetXTitle("m_{#tilde{g}} (GeV)");
  hexfakesig3->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext3 = new TLatex(500.,750.,"Eff. Fakes, H_{T} > 200 GeV, MET > 120 GeV");
  mytext3->SetTextSize(0.03);
  mytext3->Draw("same");
  mytext->Draw("same");


  hexfakesig4->Draw("COLZ");
  hexfakesig4->SetXTitle("m_{#tilde{g}} (GeV)");
  hexfakesig4->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext4 = new TLatex(500.,750.,"Eff. Fakes, H_{T} > 80 GeV, MET > 100 GeV");
  mytext4->SetTextSize(0.03);
  mytext4->Draw("same");
  mytext->Draw("same");

*/




}
