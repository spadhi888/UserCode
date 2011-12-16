{
  gStyle->SetOptStat(0);

  // Normalization


  TFile *fLM0x = new TFile("stopeff.root");
  TH2F* normal = hexbase;
  normal->SetDirectory(0);
  fLM0x->Close();


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

  TFile *fLM1 = new TFile("histo_yield_PDFMAX.root");
  TH2F* pdfUpnom = LM0x_hm12m0_all;
  TH2F* pdfUpsig1 = LM0x_hm12m01sig_all;
  TH2F* pdfUpsig2 = LM0x_hm12m02sig_all;
  TH2F* pdfUpsig3 = LM0x_hm12m03sig_all;
  TH2F* pdfUpsig4 = LM0x_hm12m04sig_all;
  pdfUpnom->SetDirectory(0);
  pdfUpsig1->SetDirectory(0);
  pdfUpsig2->SetDirectory(0);
  pdfUpsig3->SetDirectory(0);
  pdfUpsig4->SetDirectory(0);
  fLM1->Close();

  TFile *fLM2 = new TFile("histo_yield_PDFMIN.root");
  TH2F* pdfDnnom = LM0x_hm12m0_all;
  TH2F* pdfDnsig1 = LM0x_hm12m01sig_all;
  TH2F* pdfDnsig2 = LM0x_hm12m02sig_all;
  TH2F* pdfDnsig3 = LM0x_hm12m03sig_all;
  TH2F* pdfDnsig4 = LM0x_hm12m04sig_all;
  pdfDnnom->SetDirectory(0);
  pdfDnsig1->SetDirectory(0);
  pdfDnsig2->SetDirectory(0);
  pdfDnsig3->SetDirectory(0);
  pdfDnsig4->SetDirectory(0);
  fLM2->Close();


  TFile *fLM3 = new TFile("histo_yield_PDFOther.root");
  TH2F* pdfOthnom = LM0x_hm12m0_all;
  TH2F* pdfOthsig1 = LM0x_hm12m01sig_all;
  TH2F* pdfOthsig2 = LM0x_hm12m02sig_all;
  TH2F* pdfOthsig3 = LM0x_hm12m03sig_all;
  TH2F* pdfOthsig4 = LM0x_hm12m04sig_all;
  pdfOthnom->SetDirectory(0);
  pdfOthsig1->SetDirectory(0);
  pdfOthsig2->SetDirectory(0);
  pdfOthsig3->SetDirectory(0);
  pdfOthsig4->SetDirectory(0);
  fLM3->Close();



  TFile *f = new TFile("stoppdfeff.root","RECREATE"); 

   TH2F* hexpdfbase = new TH2F("hexpdfbase","hexpdfbase", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfsig1 = new TH2F("hexpdfsig1","hexpdfsig1", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfsig2 = new TH2F("hexpdfsig2","hexpdfsig2", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfsig3 = new TH2F("hexpdfsig3","hexpdfsig3", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfsig4 = new TH2F("hexpdfsig4","hexpdfsig4", 60, 0.,1500., 40, 0. , 1000.);


   TH2F* hexpdfothbase = new TH2F("hexpdfothbase","hexpdfothbase", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfothsig1 = new TH2F("hexpdfothsig1","hexpdfothsig1", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfothsig2 = new TH2F("hexpdfothsig2","hexpdfothsig2", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfothsig3 = new TH2F("hexpdfothsig3","hexpdfothsig3", 60, 0.,1500., 40, 0. , 1000.);
   TH2F* hexpdfothsig4 = new TH2F("hexpdfothsig4","hexpdfothsig4", 60, 0.,1500., 40, 0. , 1000.);

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


       // PDF Up

       float  outpdfUpbase = pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0,m12) );
       float  outpdfUpsig1 = pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0,m12) );
       float  outpdfUpsig2 = pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0,m12) );
       float  outpdfUpsig3 = pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0,m12) );
       float  outpdfUpsig4 = pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0,m12) );

       if (m0 == 1000 && pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0+25,m12)) > 0 && pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0-25,m12)) > 0 && pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0,m12)) == 0) outpdfUpbase = 0.5*(pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0+25,m12)) + pdfUpnom->GetBinContent(pdfUpnom->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0+25,m12)) > 0 && pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0-25,m12)) > 0 && pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0,m12)) == 0) outpdfUpsig1 = 0.5*(pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0+25,m12)) + pdfUpsig1->GetBinContent(pdfUpsig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0+25,m12)) > 0 && pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0-25,m12)) > 0 && pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0,m12)) == 0) outpdfUpsig2 = 0.5*(pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0+25,m12)) + pdfUpsig2->GetBinContent(pdfUpsig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0+25,m12)) > 0 && pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0-25,m12)) > 0 && pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0,m12)) == 0) outpdfUpsig3 = 0.5*(pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0+25,m12)) + pdfUpsig3->GetBinContent(pdfUpsig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0+25,m12)) > 0 && pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0-25,m12)) > 0 && pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0,m12)) == 0) outpdfUpsig4 = 0.5*(pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0+25,m12)) + pdfUpsig4->GetBinContent(pdfUpsig4->FindBin(m0-25,m12)));

       // PDF Dn

       float  outpdfDnbase = pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0,m12) );
       float  outpdfDnsig1 = pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0,m12) );
       float  outpdfDnsig2 = pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0,m12) );
       float  outpdfDnsig3 = pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0,m12) );
       float  outpdfDnsig4 = pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0,m12) );

       if (m0 == 1000 && pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0+25,m12)) > 0 && pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0-25,m12)) > 0 && pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0,m12)) == 0) outpdfDnbase = 0.5*(pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0+25,m12)) + pdfDnnom->GetBinContent(pdfDnnom->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0+25,m12)) > 0 && pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0-25,m12)) > 0 && pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0,m12)) == 0) outpdfDnsig1 = 0.5*(pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0+25,m12)) + pdfDnsig1->GetBinContent(pdfDnsig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0+25,m12)) > 0 && pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0-25,m12)) > 0 && pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0,m12)) == 0) outpdfDnsig2 = 0.5*(pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0+25,m12)) + pdfDnsig2->GetBinContent(pdfDnsig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0+25,m12)) > 0 && pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0-25,m12)) > 0 && pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0,m12)) == 0) outpdfDnsig3 = 0.5*(pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0+25,m12)) + pdfDnsig3->GetBinContent(pdfDnsig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0+25,m12)) > 0 && pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0-25,m12)) > 0 && pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0,m12)) == 0) outpdfDnsig4 = 0.5*(pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0+25,m12)) + pdfDnsig4->GetBinContent(pdfDnsig4->FindBin(m0-25,m12)));


       // PDF Dn

       float  outpdfOthbase = pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0,m12) );
       float  outpdfOthsig1 = pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0,m12) );
       float  outpdfOthsig2 = pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0,m12) );
       float  outpdfOthsig3 = pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0,m12) );
       float  outpdfOthsig4 = pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0,m12) );

       if (m0 == 1000 && pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0+25,m12)) > 0 && pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0-25,m12)) > 0 && pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0,m12)) == 0) outpdfOthbase = 0.5*(pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0+25,m12)) + pdfOthnom->GetBinContent(pdfOthnom->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0+25,m12)) > 0 && pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0-25,m12)) > 0 && pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0,m12)) == 0) outpdfOthsig1 = 0.5*(pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0+25,m12)) + pdfOthsig1->GetBinContent(pdfOthsig1->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0+25,m12)) > 0 && pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0-25,m12)) > 0 && pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0,m12)) == 0) outpdfOthsig2 = 0.5*(pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0+25,m12)) + pdfOthsig2->GetBinContent(pdfOthsig2->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0+25,m12)) > 0 && pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0-25,m12)) > 0 && pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0,m12)) == 0) outpdfOthsig3 = 0.5*(pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0+25,m12)) + pdfOthsig3->GetBinContent(pdfOthsig3->FindBin(m0-25,m12)));

       if (m0 == 1000 && pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0+25,m12)) > 0 && pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0-25,m12)) > 0 && pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0,m12)) == 0) outpdfOthsig4 = 0.5*(pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0+25,m12)) + pdfOthsig4->GetBinContent(pdfOthsig4->FindBin(m0-25,m12)));



       float errbase = 0.;
       float errsig1 = 0.;
       float errsig2 = 0.;
       float errsig3 = 0.;
       float errsig4 = 0.;


       float errmstwbase = 0.;
       float errmstwsig1 = 0.;
       float errmstwsig2 = 0.;
       float errmstwsig3 = 0.;
       float errmstwsig4 = 0.;

       if (outbase > 0)  errbase = TMath::Max(fabs(outbase - outpdfUpbase)/outbase, fabs(outbase - outpdfDnbase)/outbase);
       if (outsig1 > 0)  errsig1 = TMath::Max(fabs(outsig1 - outpdfUpsig1)/outsig1, fabs(outsig1 - outpdfDnsig1)/outsig1);
       if (outsig2 > 0)  errsig2 = TMath::Max(fabs(outsig2 - outpdfUpsig2)/outsig2, fabs(outsig2 - outpdfDnsig2)/outsig2);
       if (outsig3 > 0)  errsig3 = TMath::Max(fabs(outsig3 - outpdfUpsig3)/outsig3, fabs(outsig3 - outpdfDnsig3)/outsig3);
       if (outsig4 > 0)  errsig4 = TMath::Max(fabs(outsig4 - outpdfUpsig4)/outsig4, fabs(outsig4 - outpdfDnsig4)/outsig4);

       if (outbase > 0)  errmstwbase = fabs(outbase - outpdfOthbase)/outbase;
       if (outsig1 > 0)  errmstwsig1 = fabs(outsig1 - outpdfOthsig1)/outsig1;
       if (outsig2 > 0)  errmstwsig2 = fabs(outsig2 - outpdfOthsig2)/outsig2;
       if (outsig3 > 0)  errmstwsig3 = fabs(outsig3 - outpdfOthsig3)/outsig3;
       if (outsig4 > 0)  errmstwsig4 = fabs(outsig4 - outpdfOthsig4)/outsig4;

       binno = nom->FindBin(m0,m12);
       //       if (outbase > 0) cout << m0 << " " << m12 << " " << errbase << " " << normal->GetBinContent(normal->FindBin(m0,m12) ) << endl;
       
       float normerr = normal->GetBinContent(normal->FindBin(m0,m12) );

       hexpdfbase->SetBinContent(binno, errbase*normerr);
       hexpdfsig1->SetBinContent(binno, errsig1*normerr);
       hexpdfsig2->SetBinContent(binno, errsig2*normerr);
       hexpdfsig3->SetBinContent(binno, errsig3*normerr);
       hexpdfsig4->SetBinContent(binno, errsig4*normerr);


       hexpdfothbase->SetBinContent(binno, errmstwbase*normerr);
       hexpdfothsig1->SetBinContent(binno, errmstwsig1*normerr);
       hexpdfothsig2->SetBinContent(binno, errmstwsig2*normerr);
       hexpdfothsig3->SetBinContent(binno, errmstwsig3*normerr);
       hexpdfothsig4->SetBinContent(binno, errmstwsig4*normerr);

       //       cout << m0 << " " << m12 << " " << outsig1 << endl;
     }
 }

f->Write();

/*
  gStyle->SetOptStat(0);
  TLatex* mytext = new TLatex(500.,800.,"CMS, L_{int} = 0.98 fb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  // Plots 
  hexpdfbase->Draw("COLZ");
  hexpdfbase->SetXTitle("m_{#tilde{g}} (GeV)");
  hexpdfbase->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext0 = new TLatex(500.,700.,"PDF Uncert., H_{T} > 80 GeV, MET > 30 GeV");
  mytext0->SetTextSize(0.03);
  mytext0->Draw("same");
  mytext->Draw("same");


  hexpdfsig1->Draw("COLZ");
  hexpdfsig1->SetXTitle("m_{#tilde{g}} (GeV)");
  hexpdfsig1->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext1 = new TLatex(500.,750.,"PDF Uncert., H_{T} > 400 GeV, MET > 120 GeV");
  mytext1->SetTextSize(0.03);
  mytext1->Draw("same");
  mytext->Draw("same");


  hexpdfsig2->Draw("COLZ");
  hexpdfsig2->SetXTitle("m_{#tilde{g}} (GeV)");
  hexpdfsig2->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext2 = new TLatex(500.,750.,"PDF Uncert., H_{T} > 400 GeV, MET > 50 GeV");
  mytext2->SetTextSize(0.03);
  mytext2->Draw("same");
  mytext->Draw("same");


  hexpdfsig3->Draw("COLZ");
  hexpdfsig3->SetXTitle("m_{#tilde{g}} (GeV)");
  hexpdfsig3->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext3 = new TLatex(500.,750.,"PDF Uncert., H_{T} > 200 GeV, MET > 120 GeV");
  mytext3->SetTextSize(0.03);
  mytext3->Draw("same");
  mytext->Draw("same");


  hexpdfsig4->Draw("COLZ");
  hexpdfsig4->SetXTitle("m_{#tilde{g}} (GeV)");
  hexpdfsig4->SetYTitle("m_{#tilde{#chi^{0}_{1}}} (GeV)");
  TLatex* mytext4 = new TLatex(500.,750.,"PDF Uncet., H_{T} > 80 GeV, MET > 100 GeV");
  mytext4->SetTextSize(0.03);
  mytext4->Draw("same");
  mytext->Draw("same");

*/




}
