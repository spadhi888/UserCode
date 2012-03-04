#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"

// Stupid macro to make plots of the 7 events (see below)
// run     ls   event type id1 id2  pt1     pt2    nj(nb)  HT      met 
// 166888 425 507242758 3  -11 -11 87.9777 62.2746 6 (2) 563.705 83.8312 
// 175887 86  92144536  2  -11 -13 45.3336 142.789 3 (2) 171.995 31.4517 
// 180076 233 397121374 3  -11 -11 25.1588 70.3593 3 (2) 227.588 146.563 
// 166408 105 100231594 1   13  13 24.3261 62.1673 3 (2) 372.221 51.0459 
// 166565 815 726094730 2   11  13 23.916 46.4423 3 (2)  296.215 150.137 
// 179497 183 229894176 1  -13 -13 231.253 21.1616 2 (2) 252.454 34.5121 
// 176929 171 264222933 2  -11 -13 24.8385 43.6281 3 (2) 230.402 100.221 
void plotMetHet (float intLumi) {    
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    TCanvas *c1 = new TCanvas();
    c1->SetFillColor(0);
    c1->GetPad(0)->SetRightMargin(0.07);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->GetPad(0)->SetBorderSize(2);
    c1->GetPad(0)->SetLeftMargin(0.1407035);
    c1->GetPad(0)->SetTopMargin(0.08);
    c1->GetPad(0)->SetBottomMargin(0.13);

    TLatex latexText;
    latexText.SetTextSize(0.045);

    // fill the 2D histograms of MET vs HT
    TH2F* hee = new TH2F("hee","Met vs HT",600,0.,600.,200,0.,200.);
    TH2F* hem = new TH2F("hem","Met vs HT",600,0.,600.,200,0.,200.);
    TH2F* hmm = new TH2F("hmm","Met vs HT",600,0.,600.,200,0.,200.);

    hee->Fill(563.705, 83.8312);
    hem->Fill(171.995, 31.4517);
    hee->Fill(227.588, 146.563);
    hmm->Fill(372.221, 51.0459);
    hem->Fill(296.215, 150.137);
    hmm->Fill(252.454, 34.5121);
    hem->Fill(230.402, 100.221);

    hee->SetMarkerStyle(20); // circles
    hem->SetMarkerStyle(21); // squares
    hmm->SetMarkerStyle(22); // triangles

    hee->SetMarkerSize(2);
    hem->SetMarkerSize(2);
    hmm->SetMarkerSize(2);

    gStyle->SetOptStat(0);
    hee->GetXaxis()->SetTitle("H_{T} (GeV)");
    hee->GetYaxis()->SetTitle("#slash{E}_{T} (GeV)");
    hee->Draw();
    hem->Draw("same");
    hmm->Draw("same");
    //  TLine l1(80,0,80,200);
    // TLine l2(0,30,600,30);
    // l1.Draw();
    // l2.Draw();

    TBox b1(0,0,80,200);
    TBox b2(0,0,600,30);
    b1.SetFillColor(1);
    b1.SetFillStyle(3001);
    b2.SetFillColor(1);
    b2.SetFillStyle(3001);
    b1.Draw();
    b2.Draw();

    float xmin = hee->GetXaxis()->GetXmin();
    float xmax = hee->GetXaxis()->GetXmax();
    float ymax = hee->GetYaxis()->GetXmax();
    float x = xmin + 0.1 * (xmax-xmin);
    float y = ymax + 0.025 * ymax;
    latexText.SetTextSize(0.045);
    latexText.DrawLatex(x, y, Form("CMS Preliminary, #sqrt{s} = 7 TeV, L_{int} = %2.1f fb^{-1}", intLumi));

    TLegend *leg = new TLegend(440, 120, 560, 190, "", "br");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hee, "ee", "P");
    leg->AddEntry(hem, "e#mu", "P");
    leg->AddEntry(hmm, "#mu#mu", "P");
    leg->Draw();

    c1->Print("figs/MetVsHt.root");
    c1->Print("figs/MetVsHt.pdf");

    //============================================

    // Now the 1D histograms... this is the option
    // to plot things in bins of "events-per-10-gev" or not"
    bool per10 = true;
    float scale[3] = {1., 1., 1.};

    //-------------------------------------------------------
    // HT first
    TCanvas* cht = new TCanvas();
    cht->SetFillColor(0);
    cht->GetPad(0)->SetRightMargin(0.07);
    cht->SetFillColor(0);
    cht->SetBorderMode(0);
    cht->GetPad(0)->SetBorderSize(2);
    cht->GetPad(0)->SetLeftMargin(0.1407035);
    cht->GetPad(0)->SetTopMargin(0.08);
    cht->GetPad(0)->SetBottomMargin(0.13);
    int nht = 3;
    float ht[4]   = {80., 200., 320., 600.}; // the bins
    float htbg[3] = {2.224, 2.742, 2.630};   // the bg
    float hter[3] = {0.899, 1.234, 1.132};   // the bg err
    float htd[3]  = {1., 4., 2.};            // the data counts

    // the scale factor if we do events-per-10-gev
    if (per10) {
        for (int i=0; i<nht; i++) {
            scale[i] = 10./(ht[i+1]-ht[i]);
        }
    }

    // the histograms (bg, data, and a copy of bg)
    TH1F* htbghist   = new TH1F("htbg",  "htbg", nht, ht);
    TH1F* htdatahist = new TH1F("htdata","htdata", nht, ht);
    TH1F* htdummy    = new TH1F("htbg",  "htbg", nht, ht);

    // fill the histograms
    for (int i=0; i<nht; i++) {
        htbghist->SetBinContent(i+1, htbg[i]*scale[i]);
        htbghist->SetBinError(i+1, hter[i]*scale[i]);
        htdummy->SetBinContent(i+1, htbg[i]*scale[i]);
        htdatahist->SetBinContent(i+1,htd[i]*scale[i]);
        htdatahist->SetBinError(i+1,sqrt(htd[i])*scale[i]);
    }


    // now plot things
    htbghist->SetLineColor(kBlue);
    htbghist->SetMarkerColor(kBlue);
    htbghist->SetFillStyle(3002);
    htbghist->SetFillColor(kBlue); 
    htbghist->SetMinimum(0.);
    htbghist->GetXaxis()->SetTitle("H_{T} (GeV)");
    float ytext;
    if (per10) {
        htbghist->SetMaximum(0.6);
        htbghist->GetYaxis()->SetTitle("Events / 10 GeV");
        ytext = 0.5;
    } else {
        htbghist->SetMaximum(8.);
        ytext = 6.5;
        htbghist->GetYaxis()->SetTitle("Events");  
    }
    htdatahist->SetMarkerStyle(23);
    htdatahist->SetMarkerSize(2);
    htbghist->Draw("E2");
    htdummy->Draw("same");
    htdatahist->Draw("esamex0");
    xmin = htbghist->GetXaxis()->GetXmin();
    xmax = htbghist->GetXaxis()->GetXmax();
    ymax = htbghist->GetYaxis()->GetXmax();
    x = xmin + 0.1 * (xmax-xmin);
    y = 0.62;
    latexText.DrawLatex(x, y, Form("CMS Preliminary, #sqrt{s} = 7 TeV, L_{int} = %2.1f fb^{-1}", intLumi));
    latexText.DrawLatex(450, 0.55, "#slash{E}_{T} > 30 GeV");
    cht->Print("figs/Ht.root");
    cht->Print("figs/Ht.pdf");

    //-------------------------------------------------------
    // MET second
    TCanvas* cmet = new TCanvas();
    cmet->SetFillColor(0);
    cmet->GetPad(0)->SetRightMargin(0.07);
    cmet->SetFillColor(0);
    cmet->SetBorderMode(0);
    cmet->GetPad(0)->SetBorderSize(2);
    cmet->GetPad(0)->SetLeftMargin(0.1407035);
    cmet->GetPad(0)->SetTopMargin(0.08);
    cmet->GetPad(0)->SetBottomMargin(0.13);
    int nmet = 3;
    float met[4]   = {30., 50., 120., 200.}; // the bins
    float metbg[3] = {2.589, 4.102, 1.033};   // the bg
    float meter[3] = {1.209, 1.465, 0.592};   // the bg err
    float metd[3]  = {2., 3., 2.};            // the data counts

    // the scale factor if we do events-per-10-gev
    if (per10) {
        for (int i=0; i<nmet; i++) {
            scale[i] = 10./(met[i+1]-met[i]);
        }
    }

    // the histograms (bg, data, and a copy of bg)
    TH1F* metbghist   = new TH1F("metbg",  "metbg", nmet, met);
    TH1F* metdatahist = new TH1F("metdata","metdata", nmet, met);
    TH1F* metdummy    = new TH1F("metbg",  "metbg", nmet, met);

    // fill the histograms
    for (int i=0; i<nmet; i++) {
        metbghist->SetBinContent(i+1, metbg[i]*scale[i]);
        metbghist->SetBinError(i+1, meter[i]*scale[i]);
        metdummy->SetBinContent(i+1, metbg[i]*scale[i]);
        metdatahist->SetBinContent(i+1,metd[i]*scale[i]);
        metdatahist->SetBinError(i+1,sqrt(metd[i])*scale[i]);
    }

    // now plot things
    metbghist->SetLineColor(kBlue);
    metbghist->SetMarkerColor(kBlue);
    metbghist->SetFillStyle(3002);
    metbghist->SetFillColor(kBlue); 
    metbghist->SetMinimum(0.);
    metbghist->GetXaxis()->SetTitle("#slash{E}_{T} (GeV)");
    float ytext2;
    if (per10) {
        metbghist->SetMaximum(2.);
        metbghist->GetYaxis()->SetTitle("Events / 10 GeV");
        ytext2 = 1.7;
    } else {
        metbghist->SetMaximum(8.);
        ytext2 = 6.5;
        metbghist->GetYaxis()->SetTitle("Events");  
    }
    metdatahist->SetMarkerStyle(23);
    metdatahist->SetMarkerSize(2);
    metbghist->Draw("E2");
    metdummy->Draw("same");
    metdatahist->Draw("esamex0");
    xmin = metbghist->GetXaxis()->GetXmin();
    xmax = metbghist->GetXaxis()->GetXmax();
    ymax = metbghist->GetYaxis()->GetXmax();
    x = xmin + 0.1 * (xmax-xmin);
    y = 2.05;
    latexText.DrawLatex(x, y, Form("CMS Preliminary, #sqrt{s} = 7 TeV, L_{int} = %2.1f fb^{-1}", intLumi));
    latexText.DrawLatex(150, 1.8, "H_{T} > 80 GeV");
    cmet->Print("figs/Met.root");
    cmet->Print("figs/Met.pdf");
}
