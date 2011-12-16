{
  gStyle->SetOptStat(0);

  TFile *fLM0 = new TFile("histoFastNominal.root");
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

  TFile *fLM1 = new TFile("histo_yield_JECup.root");
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

  TFile *fLM2 = new TFile("histo_yield_JECdown.root");
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


  TFile *f = new TFile("stopjeceff.root","RECREATE"); 

   TH2F* hexjecbase = new TH2F("hexjecbase","hexjecbase", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexjecsig1 = new TH2F("hexjecsig1","hexjecsig1", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexjecsig2 = new TH2F("hexjecsig2","hexjecsig2", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexjecsig3 = new TH2F("hexjecsig3","hexjecsig3", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexjecsig4 = new TH2F("hexjecsig4","hexjecsig4", 60, 0.,1500., 40, 0. , 1000.);

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

       if (outbase > 0)  errbase = TMath::Max(fabs(outbase - outjecUpbase)/outbase, fabs(outbase - outjecDnbase)/outbase);
       if (outsig1 > 0)  errsig1 = TMath::Max(fabs(outsig1 - outjecUpsig1)/outsig1, fabs(outsig1 - outjecDnsig1)/outsig1);
       if (outsig2 > 0)  errsig2 = TMath::Max(fabs(outsig2 - outjecUpsig2)/outsig2, fabs(outsig2 - outjecDnsig2)/outsig2);
       if (outsig3 > 0)  errsig3 = TMath::Max(fabs(outsig3 - outjecUpsig3)/outsig3, fabs(outsig3 - outjecDnsig3)/outsig3);
       if (outsig4 > 0)  errsig4 = TMath::Max(fabs(outsig4 - outjecUpsig4)/outsig4, fabs(outsig4 - outjecDnsig4)/outsig4);

       binno = nom->FindBin(m0,m12);
       //       if (outbase > 0) cout << m0 << " " << m12 << " " << errbase << endl;

       hexjecbase->SetBinContent(binno, errbase);
       hexjecsig1->SetBinContent(binno, errsig1);
       hexjecsig2->SetBinContent(binno, errsig2);
       hexjecsig3->SetBinContent(binno, errsig3);
       hexjecsig4->SetBinContent(binno, errsig4);

       //       cout << m0 << " " << m12 << " " << outsig1 << endl;
     }
 }

f->Write();

/*
  gStyle->SetOptStat(0);
  TLatex* mytext = new TLatex(500.,800.,"CMS, L_{int} = 0.98 fb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  // Plots 
  hexjecbase->Draw("COLZ");
  hexjecbase->SetXTitle("m_{#tilde{g}} (GeV)");
  hexjecbase->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext0 = new TLatex(500.,700.,"JES Uncert., H_{T} > 80 GeV, MET > 30 GeV");
  mytext0->SetTextSize(0.03);
  mytext0->Draw("same");
  mytext->Draw("same");


  hexjecsig1->Draw("COLZ");
  hexjecsig1->SetXTitle("m_{#tilde{g}} (GeV)");
  hexjecsig1->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext1 = new TLatex(500.,750.,"JES Uncert., H_{T} > 400 GeV, MET > 120 GeV");
  mytext1->SetTextSize(0.03);
  mytext1->Draw("same");
  mytext->Draw("same");


  hexjecsig2->Draw("COLZ");
  hexjecsig2->SetXTitle("m_{#tilde{g}} (GeV)");
  hexjecsig2->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext2 = new TLatex(500.,750.,"JES Uncert., H_{T} > 400 GeV, MET > 50 GeV");
  mytext2->SetTextSize(0.03);
  mytext2->Draw("same");
  mytext->Draw("same");


  hexjecsig3->Draw("COLZ");
  hexjecsig3->SetXTitle("m_{#tilde{g}} (GeV)");
  hexjecsig3->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext3 = new TLatex(500.,750.,"JES Uncert., H_{T} > 200 GeV, MET > 120 GeV");
  mytext3->SetTextSize(0.03);
  mytext3->Draw("same");
  mytext->Draw("same");


  hexjecsig4->Draw("COLZ");
  hexjecsig4->SetXTitle("m_{#tilde{g}} (GeV)");
  hexjecsig4->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext4 = new TLatex(500.,750.,"JES Uncet., H_{T} > 80 GeV, MET > 100 GeV");
  mytext4->SetTextSize(0.03);
  mytext4->Draw("same");
  mytext->Draw("same");

*/




}
