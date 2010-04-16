{
  // Exclusion based on 100 pb-1
//---LEP curves
  double ch_m0[] = {0,0,600,1000};
  double ch_m12[] = {0,160,140,120};
  TGraph* ch_gr = new TGraph(4,ch_m0,ch_m12);

  double sl_m0[] = {0,0,30,50,60,75,80,90,100};
  double sl_m12[] = {0,245,240,220,200,150,100,50,0};

  TGraph* sl_gr = new TGraph(9,sl_m0,sl_m12);

  ch_gr->SetFillStyle(3325);
  sl_gr->SetFillStyle(3325);

  //--tevatron curves

//First Blob
Double_t x1[15] =  { 0.1,   3,   20,  34,  50,  58,  69,  80,  91, 100, 106, 110, 105,  82., 60.};
Double_t y1[15] = {170, 184,  194, 207, 225, 234, 240, 244, 246, 246, 240, 236, 222, 190, 160.1};

//Second Blob
Double_t x2[10] = { 70,       84,    92.5, 107,  120  ,130 , 140  ,150 , 160, 170 };
Double_t y2[10] = {160,   176.47,  186.47, 204,  217  ,211 , 200  , 190 , 179,160 };


TGraph *cdfn1 = new TGraph(15,x1,y1);
TGraph *cdfn2 = new TGraph(10,x2,y2);

cdfn1->SetFillColor(4);
cdfn1->SetLineColor(4);
cdfn1->SetLineWidth(2);
cdfn1->SetMarkerColor(4);
cdfn1->SetMarkerStyle(21);

cdfn2->SetFillColor(4);
cdfn2->SetLineColor(4);
cdfn2->SetLineWidth(2);
cdfn2->SetMarkerColor(4);
cdfn2->SetMarkerStyle(21);


int nss0 = 9;
 double xss0[nss0]={0,50,100,150,200,250,300,350,400};
 double yss0[nss0]={240.00,240.00,260.00,260.00,200.00,180.00,180.00,200.00,140.00};
TGraph *g_ss0 = new TGraph(nss0,xss0, yss0);


g_ss0->SetLineColor(6);
g_ss0->SetLineStyle(1);
g_ss0->SetLineWidth(4);
g_ss0->SetMarkerColor(kRed);
g_ss0->SetMarkerStyle(20);
g_ss0->SetMarkerSize(0.7);
g_ss0->SetMarkerSize(0.7);
g_ss0->SetFillStyle(3004);


 int nss1 = 11;
 double xss1[nss1]={0,50,100,150,200,250,300,350,400,450,500};
 double yss1[nss1]={240.00,220.00,240.00,180.00,160.00,160.00,160.00,140.00,140.00,140.00,100.00};
 TGraph *g_ss1 = new TGraph(nss1,xss1, yss1);

g_ss1->SetLineColor(2);
g_ss1->SetLineStyle(2);
g_ss1->SetLineWidth(4);
g_ss1->SetMarkerColor(kBlue);
g_ss1->SetMarkerStyle(21);
g_ss1->SetMarkerSize(0.7);


 vC1 = new TCanvas("vC1","square",0,10,700,700);
 vC1->cd(1);
 TH1F *vFrame = gPad->DrawFrame(0,100,400,400);
 vFrame->SetXTitle("m_{0} (GeV)");
 vFrame->SetYTitle("m_{1/2} (GeV)");
 vFrame->SetTitle("Exclusion using msugra scan: tan\\beta = 3, A0 = 0, \\mu > 0 ");
 vFrame->SetLabelSize(0.03, "Y");
 vFrame->SetLabelSize(0.03, "X");
 
 ch_gr->Draw("F SAME");
 sl_gr->Draw("F SAME");
 cdfn1->Draw("L SAME");
 cdfn1->Draw("F1 SAME");
 
 cdfn2->Draw("L SAME");
 cdfn2->Draw("F1 SAME");
 
 g_ss0->Draw("C SAME");
 g_ss1->Draw("C SAME");

 //g_cdf->Draw("P");
 //g_cdf1->Draw("P");
 // TText *t1 = new TLatex(100,360,"CMS Preliminary @ 7 TeV, 100 pb^{-1} ");
 // t1->SetTextSize(0.028);
 // t1->Draw();
 
 TLegend *leg = new TLegend(0.4,0.7, 0.9, 0.97, "CMS Preliminary, \\sqrt{s} = 7 TeV, 100 pb^{-1}");
 leg.SetFillStyle(0);
 leg.SetBorderSize(0);
 leg.SetTextSize(0.03);

 leg->AddEntry(g_ss0,"Observe 0 event","l");
 leg->AddEntry(g_ss1,"Observe 1 event","l");
 leg->AddEntry(cdfn1,"CDF 3l 95% CL (3.2 fb^{-1})","f");
 leg->AddEntry(ch_gr,"LEP excluded regions","f");
 leg->Draw();
// vC1->SaveAs("exclusion100ss.eps");
 vC1->SaveAs("exclusion100ss.png");

}
