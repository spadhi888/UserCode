#include "TGraph.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include "TLatex.h"
#include "TLine.h"
#include <iostream>
#include "TTree.h"
#include "TDirectory.h"
#include "TCanvas.h"

using namespace std;

void makePlots (int lsp_mass) {

    TFile *file = TFile::Open(Form("ntuple_%d.root", lsp_mass));
    file->cd();
    TTree *tree = (TTree*)gDirectory->Get("tree");

    // Makes basic histograms from the root tree  
    gStyle->SetPadRightMargin(0.16);   // default 0.1
    gStyle->SetTitleOffset(1.20, "Y");  // default 1
    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.30, "Z");  // default 1

    // This binning insures that the bins are nicely centered
    Double_t xbinsize = 50.;
    Double_t ybinsize = 50.;
    Double_t xmin     = 400. - xbinsize/2.;
    Double_t xmax     = 1100 + xbinsize/2.;
    Double_t ymin     = 280. - ybinsize/2;
    Double_t ymax     = 1000. + ybinsize/2.;
    if (lsp_mass == 150) ymin = ymin+50.;
    Int_t nx          = (int)(xmax-xmin)/xbinsize;
    Int_t ny          = (int)(ymax-ymin)/ybinsize;

    // Upper limit as a function of gluino and lsp mass
    TH2F* ul = new TH2F("ul","ul",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    // Map of excluded points with different xsection assumptions
    TH2F* excl = new TH2F("excl","excl",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
    TH2F* exclup = new TH2F("exclup","exclup",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
    TH2F* excldwn = new TH2F("excldwn","excldwn",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    // map of best region
    TH2F* ulbest = new TH2F("ulbest","ulbest",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    // cross limit section maps
    TH2F* hxsec = new TH2F("hxsec","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
    TH2F* hxsecup = new TH2F("hxsecup","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
    TH2F* hxsecdwn = new TH2F("hxsecdwn","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    // acceptance map
    TH2F* acc = new TH2F("acc","acc",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    //an empty histogram
    TH2F* empty = new TH2F("empty","empty",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

    tree->Draw("lspmass:glmass>>ul"       , "explimsrb/(5000.*effsrb)"         );
    tree->Draw("lspmass:glmass>>ulbest"   , "bestsr"                           );
    tree->Draw("lspmass:glmass>>acc"      , "effsrb"                      );
    tree->Draw("lspmass:glmass>>excl"     , "explimsrb/(5000.*effsrb)<xsec"    );
    tree->Draw("lspmass:glmass>>exclup"   , "explimsrb/(5000.*effsrb)<xsecup"  );
    tree->Draw("lspmass:glmass>>excldwn"  , "explimsrb/(5000.*effsrb)<xsecdwn" );
    tree->Draw("lspmass:glmass>>hxsec"    , "xsec"                             );
    tree->Draw("lspmass:glmass>>hxsecup"  , "xsecup"                           );
    tree->Draw("lspmass:glmass>>hxsecdwn" , "xsecdwn"                          );


    // now scan the upper limit histogram and interpolate the points for the limit
    // - scan in X as a function of Y
    // - as soon as you find a cell which is excluded, but the next cell is not excluded stop
    // - interpolate the value of X, and do a cout of the x and y coordinates
    // Do not trust it too much, it always has to be checked by hand
    // We start scanning from the top
    vector<float> xvec;
    vector<float> yvec;
    vector<float> xvecup;
    vector<float> yvecup;
    vector<float> xvecdwn;
    xvecdwn.reserve(5);
    vector<float> yvecdwn;
    yvecdwn.reserve(5);
    for (int ixsec=0; ixsec<3; ixsec++) {
        for (int iy=ny; iy>=1; iy--) {
            float y = ymin + (iy-0.5)*ybinsize;
            for (int ix=1; ix<nx; ix++) {
                float x = xmin + (ix-0.5)*xbinsize;
                // cout << ix << " " << x << " " << excl->GetBinContent(ix,1) << endl;
                float thisBinLimit = ul->GetBinContent(ix,iy);
                float thisXsec;
                if (ixsec ==0) thisXsec     = hxsec->GetBinContent(ix,iy);
                if (ixsec ==1) thisXsec     = hxsecup->GetBinContent(ix,iy);
                if (ixsec ==2) thisXsec     = hxsecdwn->GetBinContent(ix,iy);
                if (thisBinLimit < thisXsec) {
                    float nextBinLimit = ul->GetBinContent(ix+1,iy);
                    float nextXsec;
                    if (ixsec ==0) nextXsec     = hxsec->GetBinContent(ix+1,iy);
                    if (ixsec ==1) nextXsec     = hxsecup->GetBinContent(ix+1,iy);
                    if (ixsec ==2) nextXsec     = hxsecdwn->GetBinContent(ix+1,iy);
                    if (nextBinLimit > nextXsec) {
                        float d = nextBinLimit - thisBinLimit - nextXsec + thisXsec;
                        float xlim = (thisXsec - thisBinLimit)*xbinsize/d + x;
                        cout << xlim << " " << y << endl;
                        // fill the vectors that define the line	    
                        if (ixsec == 0) {
                            xvec.push_back(xlim);
                            yvec.push_back(y);
                        } else if (ixsec==1) {
                            xvecup.push_back(xlim);
                            yvecup.push_back(y);
                        } else if (ixsec==2) {
                            xvecdwn.push_back(xlim);
                            yvecdwn.push_back(y);
                        }
                        continue;
                    }
                }
            }
        }
    }

    // Add points by hand to make it go down to the bottom
    if (xvec.size() > 2 && xvecdwn.size() > 2 && xvecup.size() > 2) {
        float yy = yvec.at(yvec.size()-1) - ybinsize/2.;
        float xx = xvec.at(xvec.size()-1) + 0.5*(xvec.at(xvec.size()-1) - xvec.at(xvec.size()-2));
        xvec.push_back(xx);
        yvec.push_back(yy);
        yy = yvecup.at(yvecup.size()-1) - ybinsize/2.;
        xx = xvecup.at(xvecup.size()-1) + 0.5*(xvecup.at(xvecup.size()-1) - xvecup.at(xvecup.size()-2));
        xvecup.push_back(xx);
        yvecup.push_back(yy);
        yy = yvecdwn.at(yvecdwn.size()-1) - ybinsize/2.;
        xx = xvecdwn.at(xvecdwn.size()-1) + 0.5*(xvecdwn.at(xvecdwn.size()-1) - xvecdwn.at(xvecdwn.size()-2));
        xvecup.push_back(xx);
        yvecup.push_back(yy);        
    }
    else {
        float yy = yvec.at(yvec.size()-1) - ybinsize/2.;
        float xx = xvec.at(xvec.size()-1) + xbinsize/2.;
        xvec.push_back(xx);
        yvec.push_back(yy);
        yy = yvecup.at(yvecup.size()-1) - ybinsize/2.;
        xx = xvecup.at(xvecup.size()-1) + xbinsize/2.;
        xvecup.push_back(xx);
        yvecup.push_back(yy);
        yy = yvecdwn.at(yvecdwn.size()-1) - ybinsize/2.;
        xx = xvecdwn.at(xvecdwn.size()-1) + xbinsize/2.;
        xvecup.push_back(xx);
        yvecup.push_back(yy);        
    }

    cout << "print out some points on the central line..." << endl;
    for (unsigned int idx = 0; idx < xvec.size(); idx++) {
        cout << xvec.at(idx) << ", " << yvec.at(idx) << endl;
    }

    cout << "print out some points on the up line..." << endl;
    for (unsigned int idx = 0; idx < xvecup.size(); idx++) {
        cout << xvecup.at(idx) << ", " << yvecup.at(idx) << endl;
    }

    cout << "print out some points on the down line..." << endl;
    for (unsigned int idx = 0; idx < xvecdwn.size(); idx++) {
        cout << xvecdwn.at(idx) << ", " << yvecdwn.at(idx) << endl;
    }

    // The points are stored in a TGraph
    TGraph* g    = new TGraph(xvec.size(),   &xvec[0],   &yvec[0]);
    TGraph* gup  = new TGraph(xvecup.size(), &xvecup[0], &yvecup[0]);
    TGraph* gdwn;
    if (lsp_mass != 150)
        gdwn = new TGraph(xvecdwn.size(),&xvecdwn[0],&yvecdwn[0]);
    else {
        gdwn = new TGraph(xvecdwn.size()+2);
        for (unsigned int idx = 0; idx < xvecdwn.size(); idx++) {
            if (idx < 3)
                gdwn->SetPoint(idx, xvecdwn.at(idx), yvecdwn.at(idx));
            else
                gdwn->SetPoint(idx+2, xvecdwn.at(idx), yvecdwn.at(idx));
        }
        gdwn->SetPoint(3, 795., 415.);
        gdwn->SetPoint(4, 800., 397.);
        gdwn->SetPoint(5, 801.5, 370.);
        Double_t testx = 0;
        Double_t testy = 0;
        for (int jdx = 0; jdx < gdwn->GetN(); jdx++) {
                gdwn->GetPoint(jdx, testx, testy);
                cout << "i, x, y: " << jdx << ", " << testx << ", " << testy << endl;
            }        
    }

    g->SetLineColor(4);
    g->SetLineWidth(3);
    gup->SetLineColor(4);
    gup->SetLineWidth(3);
    gup->SetLineStyle(2);
    gdwn->SetLineColor(4);
    gdwn->SetLineWidth(3);
    gdwn->SetLineStyle(2);

    // This line is the kinematical limit
    float mtop = 175.;
    TLine kinlim = TLine(ymin+mtop, ymin, xmax, xmax-mtop);
    kinlim.SetLineWidth(3);

    // Axis labels, a bit primitive for now
    char* glmass  = "m(#tilde{g}) GeV";
    char* lspmass = "m(#tilde{t}_{1}) GeV";
    char* ztitle  = "#sigma_{UL} pb";
    excl->GetYaxis()->SetTitle(lspmass);
    excl->GetXaxis()->SetTitle(glmass);
    excl->SetTitle("T2 model  Excluded points in red");
    ul->GetYaxis()->SetTitle(lspmass);
    ul->GetXaxis()->SetTitle(glmass);
    ul->GetZaxis()->SetTitle(ztitle);
    ul->SetTitle("T2 model  Cross-section upper limits (pb)");
    empty->GetYaxis()->SetTitle(lspmass);
    empty->GetXaxis()->SetTitle(glmass);
    empty->SetTitle("T2 model");
    ulbest->GetYaxis()->SetTitle(lspmass);
    ulbest->GetXaxis()->SetTitle(glmass);
    ulbest->SetTitle("T2  Best region based on exp. limit");
    acc->GetYaxis()->SetTitle(lspmass);
    acc->GetXaxis()->SetTitle(glmass);
    acc->SetTitle("T2 model  Acc*Eff*BR");

    // A polyline with the 50 smoothed limit
    float x50  = 50. + 2. * mtop ;
    float y50c = 50. + mtop;
    float y50m = xmax - mtop;
    TPolyLine *p50 = new TPolyLine();
    p50->SetLineColor(1);
    p50->SetFillStyle(3307);
    p50->SetFillColor(1);
    p50->SetNextPoint(825, ymin);
    p50->SetNextPoint(850, 320);
    p50->SetNextPoint(855, 350);
    p50->SetNextPoint(855, 450);
    p50->SetNextPoint(845, 500);
    p50->SetNextPoint(835, 550);
    p50->SetNextPoint(840, 660);
    p50->SetNextPoint(795, 620);
    p50->SetNextPoint(795, 400);
    // p50->SetNextPoint(805, 350);
    p50->SetNextPoint(810, 335);
    // p50->SetNextPoint(815, 325);
    p50->SetNextPoint(795, ymin);

  // A polyline with the 150 smoothed limit
    float x150  = 150. + 2. * mtop;
    float y150c = 150. + mtop;
    float y150m = xmax - mtop;
    TPolyLine *p150 = new TPolyLine();
    p150->SetLineColor(1);
    p150->SetFillStyle(3307);
    //    p150->SetFillStyle(3344);
    p150->SetFillColor(1);
    p150->SetNextPoint(845, ymin);
    p150->SetNextPoint(845, y150c);
    p150->SetNextPoint(845, ymin);
    p150->SetNextPoint(845, 470);
    p150->SetNextPoint(845, 500);
    p150->SetNextPoint(830, 560);
    p150->SetNextPoint(830, 650);
    p150->SetNextPoint(775, 600);
    p150->SetNextPoint(780, 530);
    p150->SetNextPoint(790, 430);
    p150->SetNextPoint(805, 400);
    p150->SetNextPoint(800, y150c);
    p150->SetNextPoint(800, ymin);

    // Some text
    TLatex gg;
    gg.SetTextSize(0.035);

    TLatex gg2;
    gg2.SetTextSize(0.035);

    TLatex latexLabel;
    latexLabel.SetTextSize(0.035);

    const char *selection       = "Same Sign dileptons with btag selection";
    const char *obligatory_text = "CMS, #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}";
    const char *central_text    = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL}";
    const char *bands_text      = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL} #pm 1 #sigma";

    TLine l1 = TLine(xmin+0.1*(xmax-xmin), ymax-0.22*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.22*(ymax-ymin));
    l1.SetLineColor(4);
    l1.SetLineWidth(3);

    TLine l2 = TLine(xmin+0.1*(xmax-xmin), ymax-0.30*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.30*(ymax-ymin));
    l2.SetLineColor(4);
    l2.SetLineWidth(3);
    l2.SetLineStyle(2);

    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    gStyle->SetPadTickY(1);

    // Draw the exclusion map and the limit lines to make sure that they make sense
    TCanvas* c11 = new TCanvas();
    c11->SetFillColor(0);
    c11->GetPad(0)->SetRightMargin(0.07);
    c11->SetFillColor(0);
    c11->SetBorderMode(0);
    c11->GetPad(0)->SetBorderSize(2);
    c11->GetPad(0)->SetLeftMargin(0.1407035);
    c11->GetPad(0)->SetTopMargin(0.08);
    c11->GetPad(0)->SetBottomMargin(0.13);

    // need to fix the double filling that's going on for mlsp = 150 GeV
    if (lsp_mass == 150) {
        for (unsigned int idx = 1; idx < excl->GetXaxis()->GetNbins()+1; idx++) {
            for (int idy = 1; idy < excl->GetYaxis()->GetNbins()+1; idy++) {
                
                if (excl->GetBinContent(idx, idy) > 1)
                    excl->SetBinContent(idx, idy, 1);
            }
        }
    }

    excl->Draw("col");
    if (lsp_mass == 50) {
        TGraph cgraph = TGraph(4);
        cgraph.SetPoint(0, 800, ymin);
        cgraph.SetPoint(1, 840, 340);
        cgraph.SetPoint(2, 820, 500);
        cgraph.SetPoint(3, 805, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 820, ymin);
        ugraph.SetPoint(1, 860, 340);
        ugraph.SetPoint(2, 840, 480);
        ugraph.SetPoint(3, 820, 640);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        TGraph dgraph = TGraph(5);
        dgraph.SetPoint(0, 780, ymin);
        dgraph.SetPoint(1, 800, 300);
        dgraph.SetPoint(2, 820, 340);
        dgraph.SetPoint(3, 795, 500);
        dgraph.SetPoint(4, 780, 600);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("GlStop_%d_ExcludedRegionMap.pdf", lsp_mass));
    }
    else if (lsp_mass == 150) {
        float mtop = 175.;
        float mb   = 5;
        float x150  = 150. + 2. * mtop;
        float y150c = 150. + mtop;
        float y150m = xmax - mtop;
        TGraph cgraph = TGraph(3);
        cgraph.SetPoint(0, 820, y150c);
        cgraph.SetPoint(1, 830, 380);
        cgraph.SetPoint(2, 800, 490);
        cgraph.SetPoint(3, 800, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph dgraph = TGraph(8);
        dgraph.SetPoint(0, 765, 590);
        dgraph.SetPoint(1, 780, 410);
        dgraph.SetPoint(2, 790, 405);
        dgraph.SetPoint(3, 795, 402);
        dgraph.SetPoint(4, 797, 401);
        dgraph.SetPoint(5, 798, 400);
        dgraph.SetPoint(6, 799, 398);
        dgraph.SetPoint(7, 800, y150c);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 845, y150c);
        ugraph.SetPoint(1, 855, 450);
        ugraph.SetPoint(2, 830, 540);
        ugraph.SetPoint(3, 830, 655);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("GlStop_%d_ExcludedRegionMap.pdf", lsp_mass));
    }
    else {
        // float mtop = 175.;
        // float mb   = 5;
        // float x50  = 50. + 2. * mtop ;
        // float y50c = 50. + mtop;
        // float y50m = xmax - mtop;
        // TLine l1_50 = TLine(x50, y50c, xmax, y50c);
        // TLine l2_50 = TLine(x50, y50c, xmax, y50m);
        // float x150  = 150. + 2. * mtop;
        // float y150c = 150. + mtop;
        // float y150m = xmax - mtop;
        // TLine l1_150 = TLine(x150, y150c, xmax, y150c);
        // TLine l2_150 = TLine(x150, y150c, xmax, y150m);
        g->Draw("samePC");
        gup->Draw("samePC");
        gdwn->Draw("samePC");     
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("GlStop_%d_ExcludedRegionMap.pdf", lsp_mass));
    }
    // g->Draw("samePC");
    // gup->Draw("samePC");
    // gdwn->Draw("samePC");
    // kinlim.Draw();

    // Draw the limit lines on top of the temperature plot
    TCanvas* c12 = new TCanvas();
    c12->SetFillColor(0);
    c12->SetFillColor(0);
    c12->SetBorderMode(0);
    c12->GetPad(0)->SetBorderSize(2);
    c12->GetPad(0)->SetLeftMargin(0.1407035);
    c12->GetPad(0)->SetTopMargin(0.08);
    c12->GetPad(0)->SetBottomMargin(0.13);

    // need to fix the double filling that's going on for mlsp = 150 GeV
    if (lsp_mass == 150) {
        for (unsigned int idx = 0; idx < 3; idx++) {
            int bin = ul->FindBin(750. + 50. * idx, 350.);
            float binc = ul->GetBinContent(bin);
            ul->SetBinContent(bin, binc * 0.5);
        }
    }

    ul->Draw("colz");
    if (lsp_mass == 150) {
        float mtop = 175.;
        float mb   = 5;
        float x150  = 150. + 2. * mtop;
        float y150c = 150. + mtop;
        float y150m = xmax - mtop;
        TGraph cgraph = TGraph(3);
        cgraph.SetPoint(0, 820, y150c);
        cgraph.SetPoint(1, 830, 380);
        cgraph.SetPoint(2, 800, 490);
        cgraph.SetPoint(3, 800, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph dgraph = TGraph(8);
        dgraph.SetPoint(0, 765, 590);
        dgraph.SetPoint(1, 780, 410);
        dgraph.SetPoint(2, 790, 405);
        dgraph.SetPoint(3, 795, 402);
        dgraph.SetPoint(4, 797, 401);
        dgraph.SetPoint(5, 798, 400);
        dgraph.SetPoint(6, 799, 398);
        dgraph.SetPoint(7, 800, y150c);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 845, y150c);
        ugraph.SetPoint(1, 855, 450);
        ugraph.SetPoint(2, 830, 540);
        ugraph.SetPoint(3, 830, 655);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("GlStop_%d_LimitsOnCarpet.pdf", lsp_mass));
    }
    else if (lsp_mass == 50){
        TGraph cgraph = TGraph(4);
        cgraph.SetPoint(0, 800, ymin);
        cgraph.SetPoint(1, 840, 340);
        cgraph.SetPoint(2, 820, 500);
        cgraph.SetPoint(3, 805, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 820, ymin);
        ugraph.SetPoint(1, 860, 340);
        ugraph.SetPoint(2, 840, 480);
        ugraph.SetPoint(3, 820, 640);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        TGraph dgraph = TGraph(5);
        dgraph.SetPoint(0, 780, ymin);
        dgraph.SetPoint(1, 800, 300);
        dgraph.SetPoint(2, 820, 340);
        dgraph.SetPoint(3, 795, 500);
        dgraph.SetPoint(4, 780, 600);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("GlStop_%d_LimitsOnCarpet.pdf", lsp_mass));        
    }
    else {
        // float mtop = 175.;
        // float mb   = 5;
        // float x50  = 50. + 2. * mtop ;
        // float y50c = 50. + mtop;
        // float y50m = xmax - mtop;
        // TLine l1_50 = TLine(x50, y50c, xmax, y50c);
        // TLine l2_50 = TLine(x50, y50c, xmax, y50m);
        // float x150  = 150. + 2. * mtop;
        // float y150c = 150. + mtop;
        // float y150m = xmax - mtop;
        // TLine l1_150 = TLine(x150, y150c, xmax, y150c);
        // TLine l2_150 = TLine(x150, y150c, xmax, y150m);
        g->Draw("samePC");
        gup->Draw("samePC");
        gdwn->Draw("samePC");     
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("GlStop_%d_LimitsOnCarpet.pdf", lsp_mass));
    }

    //Draw the limit lines and nothing else
    TCanvas* c13 = new TCanvas();
    empty->Draw();
    c13->SetFillColor(0);
    c13->GetPad(0)->SetRightMargin(0.07);
    c13->SetFillColor(0);
    c13->SetBorderMode(0);
    c13->GetPad(0)->SetBorderSize(2);
    c13->GetPad(0)->SetLeftMargin(0.1407035);
    c13->GetPad(0)->SetTopMargin(0.08);
    c13->GetPad(0)->SetBottomMargin(0.13);

    if (lsp_mass == 50) {
        float mtop = 175.;
        float mb   = 5;
        float x150  = 150. + 2. * mtop;
        float y150c = 150. + mtop;
        float y150m = xmax - mtop;
        TGraph cgraph = TGraph(3);
        cgraph.SetPoint(0, 820, y150c);
        cgraph.SetPoint(1, 830, 380);
        cgraph.SetPoint(2, 800, 490);
        cgraph.SetPoint(3, 800, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph dgraph = TGraph(8);
        dgraph.SetPoint(0, 765, 590);
        dgraph.SetPoint(1, 780, 410);
        dgraph.SetPoint(2, 790, 405);
        dgraph.SetPoint(3, 795, 402);
        dgraph.SetPoint(4, 797, 401);
        dgraph.SetPoint(5, 798, 400);
        dgraph.SetPoint(6, 799, 398);
        dgraph.SetPoint(7, 800, y150c);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 845, y150c);
        ugraph.SetPoint(1, 855, 450);
        ugraph.SetPoint(2, 830, 540);
        ugraph.SetPoint(3, 830, 655);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("GlStop_%d_LimitsOnWhite.pdf", lsp_mass));
    }
    else if (lsp_mass == 50) {
        TGraph cgraph = TGraph(4);
        cgraph.SetPoint(0, 800, ymin);
        cgraph.SetPoint(1, 840, 340);
        cgraph.SetPoint(2, 820, 500);
        cgraph.SetPoint(3, 805, 625);
        cgraph.SetLineColor(4);
        cgraph.SetLineWidth(3);
        TGraph ugraph = TGraph(4);
        ugraph.SetPoint(0, 820, ymin);
        ugraph.SetPoint(1, 860, 340);
        ugraph.SetPoint(2, 840, 480);
        ugraph.SetPoint(3, 820, 640);
        ugraph.SetLineColor(4);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        TGraph dgraph = TGraph(5);
        dgraph.SetPoint(0, 780, ymin);
        dgraph.SetPoint(1, 800, 300);
        dgraph.SetPoint(2, 820, 340);
        dgraph.SetPoint(3, 795, 500);
        dgraph.SetPoint(4, 780, 600);
        dgraph.SetLineColor(4);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);    
        cgraph.Draw("Csame");
        dgraph.Draw("Csame");
        ugraph.Draw("Csame");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("GlStop_%d_LimitsOnWhite.pdf", lsp_mass));
    }
    else {
        // float mtop = 175.;
        // float mb   = 5;
        // float x50  = 50. + 2. * mtop ;
        // float y50c = 50. + mtop;
        // float y50m = xmax - mtop;
        // TLine l1_50 = TLine(x50, y50c, xmax, y50c);
        // TLine l2_50 = TLine(x50, y50c, xmax, y50m);
        // float x150  = 150. + 2. * mtop;
        // float y150c = 150. + mtop;
        // float y150m = xmax - mtop;
        // TLine l1_150 = TLine(x150, y150c, xmax, y150c);
        // TLine l2_150 = TLine(x150, y150c, xmax, y150m);
        g->Draw("samePC");
        gup->Draw("samePC");
        gdwn->Draw("samePC");     
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("GlStop_%d_LimitsOnWhite.pdf", lsp_mass));
    }
    // g->Draw("samePC");
    // gup->Draw("samePC");
    // gdwn->Draw("samePC");
    // kinlim.Draw();
   
    //Draw the best region and nothing else
    TCanvas* c14 = new TCanvas();
    c14->SetFillColor(0);
    c14->GetPad(0)->SetRightMargin(0.07);
    c14->SetFillColor(0);
    c14->SetBorderMode(0);
    c14->GetPad(0)->SetBorderSize(2);
    c14->GetPad(0)->SetLeftMargin(0.1407035);
    c14->GetPad(0)->SetTopMargin(0.08);
    c14->GetPad(0)->SetBottomMargin(0.13);

    // need to fix the double filling that's going on for mlsp = 150 GeV
    if (lsp_mass == 150) {
        for (unsigned int idx = 1; idx < ulbest->GetXaxis()->GetNbins()+1; idx++) {
            for (int idy = 1; idy < ulbest->GetYaxis()->GetNbins()+1; idy++) {
                
                if (ulbest->GetBinContent(idx, idy) > 7)
                    ulbest->SetBinContent(idx, idy, 6);
            }
        }
    }

    ulbest->Draw("textcol");

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
    kinlim.Draw();
    c14->Print(Form("GlStop_%d_BestSignalRegion.pdf", lsp_mass));

    //Draw the acceptance carpet
    TCanvas* c15 = new TCanvas();
    c15->SetFillColor(0);
    c15->SetFillColor(0);
    c15->SetBorderMode(0);
    c15->GetPad(0)->SetBorderSize(2);
    c15->GetPad(0)->SetLeftMargin(0.1407035);
    c15->GetPad(0)->SetTopMargin(0.08);
    c15->GetPad(0)->SetBottomMargin(0.13);

    // need to fix the double filling that's going on for mlsp = 150 GeV
    if (lsp_mass == 150) {
        for (unsigned int idx = 0; idx < 3; idx++) {
            int bin = acc->FindBin(750. + 50. * idx, 350.);
            float binc = acc->GetBinContent(bin);
            acc->SetBinContent(bin, binc * 0.5);
        }
    }

    acc->Draw("colz");

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
    kinlim.Draw();
    c15->Print(Form("GlStop_%d_AcceptanceCarpet.pdf", lsp_mass));


    // Draw the limit lines on top of the temperature plot...
    // but this time use the same limits as the paper
    TCanvas* c12n = new TCanvas();
    c12n->SetFillColor(0);
    c12n->SetFillColor(0);
    c12n->SetBorderMode(0);
    c12n->GetPad(0)->SetBorderSize(2);
    c12n->GetPad(0)->SetLeftMargin(0.1407035);
    c12n->GetPad(0)->SetTopMargin(0.08);
    c12n->GetPad(0)->SetBottomMargin(0.13);

    // need to fix the double filling that's going on for mlsp = 150 GeV
    if (lsp_mass == 150) {
        for (unsigned int idx = 0; idx < 3; idx++) {
            int bin = ul->FindBin(750. + 50. * idx, 350.);
            float binc = ul->GetBinContent(bin);
            ul->SetBinContent(bin, binc * 0.5);
        }
    }

    ul->Draw("colz");
    if (lsp_mass == 150) {
      p150->Draw("fl");
      p150->Draw();
    } else if (lsp_mass == 50){
      p50->Draw("fl");
      p50->Draw();
    }
    kinlim.Draw();
    TPolyLine *blah = new TPolyLine();
    blah->SetLineColor(1);
    blah->SetFillStyle(3307);
    blah->SetFillColor(1);
    blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
    blah->SetNextPoint(xmin+0.19*(xmax-xmin), ymax-0.25*(ymax-ymin));
    blah->SetNextPoint(xmin+0.19*(xmax-xmin), ymax-0.20*(ymax-ymin));
    blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.20*(ymax-ymin));
    blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
    blah->Draw("fl");
    blah->Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.01*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), Form("m(#tilde{#chi}_{1}^{0}) = %d GeV", lsp_mass));
    // gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
    c12n->Print(Form("GlStop_%d_LimitsOnCarpetLikePaper.pdf", lsp_mass));

}
