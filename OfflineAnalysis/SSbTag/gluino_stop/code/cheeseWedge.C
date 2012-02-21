{
 // Makes basic histograms from the root tree
  
  gStyle->SetPadRightMargin(0.12);   // default 0.1
  gStyle->SetTitleOffset(1.20, "Y");  // default 1
  gStyle->SetOptStat(0);

  float xmax = 1000.;
  float xmin = 400.;
  float ymax = 1000.;
  float ymin = 100.;

  TH2F* empty = new TH2F("empty","T2 model",100,xmin,xmax,100,ymin,ymax);
  empty->GetYaxis()->SetTitle("m(#tilde{t}) GeV");
  empty->GetXaxis()->SetTitle("m(#tilde{g}) GeV");


  TLatex gg;
  TLatex gg2;
  TLatex latexLabel;
  latexLabel.SetTextSize(0.04);
  char * selection ="Same Sign dilpetons with btag selection";
  char * masses = Form("m(#chi_{1}^{0}) = 50 GeV");
  gg.SetTextSize(0.04);
  gg2.SetTextSize(0.04);

  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  empty->Draw();


  // The kinematical limits
  float mtop = 175.;
  float mb   = 5;
  float x50  = 50. + 2. * mtop ;
  float y50c = 50. + mtop;
  float y50m = xmax - mtop;
  TLine l1_50 = TLine(x50, y50c, xmax, y50c);
  TLine l2_50 = TLine(x50, y50c, xmax, y50m);

  float x150  = 150. + 2. * mtop;
  float y150c = 150. + mtop;
  float y150m = xmax - mtop;
  TLine l1_150 = TLine(x150, y150c, xmax, y150c);
  TLine l2_150 = TLine(x150, y150c, xmax, y150m);

  l1_150->SetLineColor(kBlue);
  l2_150->SetLineColor(kBlue);
  l1_50->SetLineColor(kRed);
  l2_50->SetLineColor(kRed);
  l1_150->Draw();
  l2_150->Draw();
  l1_50->Draw();
  l2_50->Draw();


  // A polyline with the 150 smoothed limit
  TPolyLine *p50 = new TPolyLine();
  p50->SetLineColor(kRed);
  p50->SetFillStyle(3244);
  p50->SetFillColor(kRed);
  p50->SetNextPoint(820, y50c);
  p50->SetNextPoint(820, 255);
  p50->SetNextPoint(850, 380);
  p50->SetNextPoint(845, 500);
  p50->SetNextPoint(820, 640);
  p50->SetNextPoint(795, 615);
  p50->SetNextPoint(795, 400);
  p50->SetNextPoint(805, 350);
  p50->SetNextPoint(805, 325);
  p50->SetNextPoint(795, y50c);

  // A polyline with the 300 smoothed limit
  TPolyLine *p150 = new TPolyLine();
  p150->SetLineColor(kBlue);
  p150->SetFillStyle(3344);
  p150->SetFillColor(kBlue);
  p150->SetNextPoint(845, y150c);
  p150->SetNextPoint(825, 475);
  p150->SetNextPoint(825, 560);
  p150->SetNextPoint(815, 640);
  p150->SetNextPoint(775, 595);
  p150->SetNextPoint(785, 400);
  p150->SetNextPoint(805, y150c);

  p50->Draw("fl");
  p150->Draw("fl");

  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.08*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.16*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.24*(ymax-ymin),"#sqrt{s} = 7 TeV, L=4.7 fb^{-1} ");
  gg2.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.32*(ymax-ymin), "Exclusion #sigma^{prod} = #sigma^{NLO+NLL} #pm 1 #sigma");

  latexLabel.DrawLatex(575, 250, "#color[2]{m(#chi_{1}^{0}) = 50 GeV}");
  latexLabel.DrawLatex(575, 350, "#color[4]{m(#chi_{1}^{0}) = 150 GeV}");

}  
