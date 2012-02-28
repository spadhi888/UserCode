{
    // Makes basic histograms from the root tree
  
    gStyle->SetPadRightMargin(0.12);   // default 0.1
    gStyle->SetTitleOffset(1.20, "Y");  // default 1
    gStyle->SetOptStat(0);

    float xmax = 1000.;
    float xmin = 300.;
    float ymax = 1000.;
    float ymin = 200.;

    TH2F* empty = new TH2F("empty","Gluino Sbottom model",100,xmin,xmax,100,ymin,ymax);
    empty->GetYaxis()->SetTitle("m(#tilde{b}_{1}) GeV");
    empty->GetXaxis()->SetTitle("m(#tilde{g}) GeV");

    TLatex gg;
    gg.SetTextSize(0.035);

    TLatex gg2;
    gg2.SetTextSize(0.035);

    TLatex latexLabel;
    latexLabel.SetTextSize(0.035);
    const char *selection       = "Same Sign dileptons with btag selection";
    const char *obligatory_text = "CMS Preliminary, #sqrt{s} = 7 TeV, L_{int} = 5.0 fb^{-1}";
    const char *central_text    = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL}";
    const char *bands_text      = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL} #pm 1 #sigma";
    char  *masses               = Form("m(#tilde{#chi}^{0}_{1}) = 50 GeV");

    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    empty->Draw();

    // The kinematical limits
    float mtop   = 175.;
    float mb     = 5;
    float x150   = 150. + mtop + mb;
    float y150c  = 150. + mtop;
    float y150m  = xmax - mb;
    TLine l1_150 = TLine(x150, y150c, xmax, y150c);
    TLine l2_150 = TLine(x150, y150c, xmax, y150m);

    float x300   = 300. + mtop + mb;
    float y300c  = 300. + mtop;
    float y300m  = xmax - mb;
    TLine l1_300 = TLine(x300, y300c, xmax, y300c);
    TLine l2_300 = TLine(x300, y300c, xmax, y300m);

    l1_300->SetLineColor(kBlue);
    l2_300->SetLineColor(kBlue);
    l1_150->SetLineColor(kRed);
    l2_150->SetLineColor(kRed);
    l1_300->Draw();
    l2_300->Draw();
    l1_150->Draw();
    l2_150->Draw();

    // A polyline with the 150 smoothed limit
    TPolyLine *p150 = new TPolyLine();
    p150->SetLineColor(kRed);
    p150->SetFillStyle(3244);
    p150->SetFillColor(kRed);
    p150->SetNextPoint(780.,y150c);
    p150->SetNextPoint(780.,350.);
    p150->SetNextPoint(845.,600.);
    p150->SetNextPoint(850.,700.);
    p150->SetNextPoint(760.,760.-mb);
    p150->SetNextPoint(730.,730.-mb);
    p150->SetNextPoint(795.,600.);
    p150->SetNextPoint(760.,350.);
    p150->SetNextPoint(760.,y150c);
    p150->SetNextPoint(780.,y150c);


    // A polyline with the 300 smoothed limit
    TPolyLine *p300 = new TPolyLine();
    p300->SetLineColor(kBlue);
    p300->SetFillStyle(3344);
    p300->SetFillColor(kBlue);
    p300->SetNextPoint(830.,y300c);
    p300->SetNextPoint(830.,475.);
    p300->SetNextPoint(850.,550.);
    p300->SetNextPoint(850.,675.);
    p300->SetNextPoint(840.,750.);
    p300->SetNextPoint(780.,780.-mb);
    p300->SetNextPoint(740.,740.-mb);
    p300->SetNextPoint(800.,700.);
    p300->SetNextPoint(800.,475.);
    p300->SetNextPoint(800.,y300c);
    p150->Draw("fl");
    p300->Draw("fl");


    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin), selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
    gg2.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);

    latexLabel.DrawLatex(500, 345, "#color[2]{m(#tilde{#chi}_{1}^{+}) = 150 GeV}");
    latexLabel.DrawLatex(560, 495, "#color[4]{m(#tilde{#chi}_{1}^{+}) = 300 GeV}");
    c1->Print("B2_CheeseWedge.pdf");
}  
