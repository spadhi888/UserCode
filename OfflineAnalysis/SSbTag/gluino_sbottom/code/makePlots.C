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
#include "TPolyLine.h"

using namespace std;

void makePlots (int chi_mass) {

    TFile *file = TFile::Open(Form("ntuple_%d.root", chi_mass));
    file->cd();
    TTree *tree = (TTree*)gDirectory->Get("tree");

    Double_t mtop = 175.;
    Double_t mb   = 5.;

    // Makes basic histograms from the root tree  
    gStyle->SetPadRightMargin(0.16);   // default 0.1
    gStyle->SetTitleOffset(1.20, "Y");  // default 1
    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.20, "Z");  // default 1

    // This binning insures that the bins are nicely centered
    Double_t xbinsize = 25.;
    Double_t ybinsize = 50.;
    Double_t xmin     = 300. - xbinsize/2.;
    Double_t xmax     = 1000 + xbinsize/2.;
    Double_t ymin     = 350. - ybinsize/2;
    Double_t ymax     = 1000. + ybinsize/2.;
    if (chi_mass > 151.) {
        xbinsize = 50.;
        ybinsize = 25.;
        xmin     = 300. - xbinsize/2.;
        xmax     = 1000 + xbinsize/2.;
        ymin     = 350. - ybinsize/2;
        ymax     = 1000. + ybinsize/2.;
    }
    if (chi_mass > 201.) ymin = 475. - ybinsize/2;
    Int_t nx = (int)(xmax-xmin)/xbinsize;
    Int_t ny = (int)(ymax-ymin)/ybinsize;
  

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


    tree->Draw("lspmass:glmass>>ul"       ,"explimsrb/(5000.*effsrb)"         );
    tree->Draw("lspmass:glmass>>ulbest"   ,"bestsr"                           );
    tree->Draw("lspmass:glmass>>acc"      ,"1.*effsrb"                      );
    tree->Draw("lspmass:glmass>>excl"     ,"explimsrb/(5000.*effsrb)<xsec"    );
    tree->Draw("lspmass:glmass>>exclup"   ,"explimsrb/(5000.*effsrb)<xsecup"  );
    tree->Draw("lspmass:glmass>>excldwn"  ,"explimsrb/(5000.*effsrb)<xsecdwn" );
    tree->Draw("lspmass:glmass>>hxsec"    ,   "xsec"                          );
    tree->Draw("lspmass:glmass>>hxsecup"  , "xsecup"                          );
    tree->Draw("lspmass:glmass>>hxsecdwn" ,"xsecdwn"                          );

    // now scan the upper limit histogram and interpolate the points for the limit
    // - scan in X as a function of Y
    // - as soon as you find a cell which is excluded, but the next cell is not excluded stop
    // - interpolate the value of X, and do a cout of the x and y coordinates
    // Do not trust it too much, it always has to be cheked by hand
    // We start scanning from the top
    vector<float> xvec;
    vector<float> yvec;
    vector<float> xvecup;
    vector<float> yvecup;
    vector<float> xvecdwn;
    vector<float> yvecdwn;
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
    TGraph* gdwn = new TGraph(xvecdwn.size(),&xvecdwn[0],&yvecdwn[0]);
    g->SetLineColor(1);
    g->SetLineWidth(3);
    gup->SetLineColor(1);
    gup->SetLineWidth(3);
    gup->SetLineStyle(2);
    gdwn->SetLineColor(1);
    gdwn->SetLineWidth(3);
    gdwn->SetLineStyle(2);

    // Axis labels, a bit primitive for now
    char* glmass = "m(#tilde{g}) GeV";
    char* lspmass = "m(#tilde{b}_{1}) GeV";
    char* ztitle  = "#sigma_{UL} pb";

    excl->GetYaxis()->SetTitle(lspmass);
    excl->GetXaxis()->SetTitle(glmass);
    excl->SetTitle("B2 model  Excluded points in red");
    ul->GetYaxis()->SetTitle(lspmass);
    ul->GetXaxis()->SetTitle(glmass);
    ul->GetZaxis()->SetTitle(ztitle);
    ul->SetTitle("B2 model Cross-section upper limits (pb)");
    empty->GetYaxis()->SetTitle(lspmass);
    empty->GetXaxis()->SetTitle(glmass);
    empty->SetTitle("B2 model");
    ulbest->GetYaxis()->SetTitle(lspmass);
    ulbest->GetXaxis()->SetTitle(glmass);
    ulbest->SetTitle("B2 model  Best region based on exp. limit");
    acc->GetYaxis()->SetTitle(lspmass);
    acc->GetXaxis()->SetTitle(glmass);
    acc->SetTitle("B2 model  Acc*Eff*BR for best signal region in percent");

    // This line is the kinematical limit....make it exactly as the paper
    TLine kinlim;
    if (chi_mass == 150) {
      kinlim = TLine(150.+mtop+mb, 150.+mtop, xmax, xmax-mb);

    } else if (chi_mass == 300) {
      kinlim = TLine(300.+mtop+mb, 300.+mtop, xmax, xmax-mb);
    } else {
      kinlim = TLine(xmin, max(xmin-mb,mtop+chi_mass), xmax, max(xmax-mb,mtop+chi_mass-mb-mb-mb));
    }
    kinlim.SetLineWidth(3);


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
    const char *masses          = Form("m(#tilde{#chi}^{0}_{1}) = 50 GeV and m(#tilde{#chi}^{+}_{1}) = %d GeV", chi_mass);

    TLine l1 = TLine(xmin+0.1*(xmax-xmin), ymax-0.22*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.22*(ymax-ymin));
    l1.SetLineColor(1);
    l1.SetLineWidth(3);

    TLine l2 = TLine(xmin+0.1*(xmax-xmin), ymax-0.30*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.30*(ymax-ymin));
    l2.SetLineColor(1);
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

    excl->Draw("col");
    if (chi_mass == 150) {
        TGraph cgraph = TGraph(17);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 800, 400);
        cgraph.SetPoint(2, 820, 460);
        cgraph.SetPoint(3, 835, 580);
        cgraph.SetPoint(4, 837, 590);
        cgraph.SetPoint(5, 838, 595);
        cgraph.SetPoint(6, 839, 600);                
        cgraph.SetPoint(7, 840, 605);
        cgraph.SetPoint(8, 841, 608);
        cgraph.SetPoint(9, 842, 615);
        cgraph.SetPoint(10, 843, 620);
        cgraph.SetPoint(11, 844, 625);
        cgraph.SetPoint(12, 845, 630);
        cgraph.SetPoint(13, 846, 650);
        cgraph.SetPoint(14, 848, 680);
        cgraph.SetPoint(15, 850, 720);
        cgraph.SetPoint(16, 840, 775);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(17);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 780, 400);
        dgraph.SetPoint(2, 800, 460);
        dgraph.SetPoint(3, 815, 580);
        dgraph.SetPoint(4, 817, 590);
        dgraph.SetPoint(5, 818, 595);
        dgraph.SetPoint(6, 819, 600);                
        dgraph.SetPoint(7, 820, 605);
        dgraph.SetPoint(8, 821, 608);
        dgraph.SetPoint(9, 822, 615);
        dgraph.SetPoint(10, 823, 620);
        dgraph.SetPoint(11, 824, 625);
        dgraph.SetPoint(12, 825, 630);
        dgraph.SetPoint(13, 826, 650);
        dgraph.SetPoint(14, 828, 680);
        dgraph.SetPoint(15, 830, 710);
        dgraph.SetPoint(16, 825, 775);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(17);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 770, ymin);
        ugraph.SetPoint(1, 810, 400);
        ugraph.SetPoint(2, 835, 460);
        ugraph.SetPoint(3, 855, 580);
        ugraph.SetPoint(4, 857, 590);
        ugraph.SetPoint(5, 858, 595);
        ugraph.SetPoint(6, 859, 600);                
        ugraph.SetPoint(7, 860, 605);
        ugraph.SetPoint(8, 861, 608);
        ugraph.SetPoint(9, 862, 615);
        ugraph.SetPoint(10, 863, 620);
        ugraph.SetPoint(11, 864, 625);
        ugraph.SetPoint(12, 865, 630);
        ugraph.SetPoint(13, 866, 650);
        ugraph.SetPoint(14, 868, 680);
        ugraph.SetPoint(15, 870, 730);
        ugraph.SetPoint(16, 860, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("B2_ExcludedRegionMap_%d.pdf", chi_mass));
    }
    else if (chi_mass == 200) {
        TGraph cgraph = TGraph(18);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 805, 420);
        cgraph.SetPoint(2, 820, 450);
        cgraph.SetPoint(3, 840, 540);
        cgraph.SetPoint(4, 845, 600);
        cgraph.SetPoint(5, 846, 620);
        cgraph.SetPoint(6, 847, 640);
        cgraph.SetPoint(7, 848, 660);
        cgraph.SetPoint(8, 849, 680);
        cgraph.SetPoint(9, 849, 735);
        cgraph.SetPoint(10, 848.75, 736.25);
        cgraph.SetPoint(11, 848.5, 737.5);
        cgraph.SetPoint(12, 848, 740);
        cgraph.SetPoint(13, 847.5, 742.5);
        cgraph.SetPoint(14, 847, 745);
        cgraph.SetPoint(15, 846, 750);
        cgraph.SetPoint(16, 845, 755);
        cgraph.SetPoint(17, 840, 780);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(12);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 770, 420);
        dgraph.SetPoint(2, 800, 450);
        dgraph.SetPoint(3, 820, 540);
        dgraph.SetPoint(4, 825, 600);
        dgraph.SetPoint(5, 826, 620);
        dgraph.SetPoint(6, 827, 640);
        dgraph.SetPoint(7, 828, 660);
        dgraph.SetPoint(8, 829, 680);
        dgraph.SetPoint(9, 829, 710);
        dgraph.SetPoint(10, 830, 725);
        dgraph.SetPoint(11, 831, 740);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(19);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 780, ymin);
        ugraph.SetPoint(1, 825, 415);
        ugraph.SetPoint(2, 835, 425);
        ugraph.SetPoint(3, 845, 455);
        ugraph.SetPoint(4, 850, 470);
        ugraph.SetPoint(5, 855, 485);
        ugraph.SetPoint(6, 857, 505);
        ugraph.SetPoint(7, 860, 530);
        ugraph.SetPoint(8, 861.5, 550);
        ugraph.SetPoint(9, 863, 570);
        ugraph.SetPoint(10, 864, 585);
        ugraph.SetPoint(11, 865, 600);
        ugraph.SetPoint(12, 865, 740);
        ugraph.SetPoint(13, 864.5, 742);
        ugraph.SetPoint(14, 863.75, 743.75);
        ugraph.SetPoint(15, 862.5, 747.5);
        ugraph.SetPoint(16, 860, 755);
        ugraph.SetPoint(17, 855, 770);
        ugraph.SetPoint(18, 850, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("B2_ExcludedRegionMap_%d.pdf", chi_mass));        
    }
    else if (chi_mass == 300) {
        TGraph cgraph = TGraph(15);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 820, ymin);
        cgraph.SetPoint(1, 820, 635);
        cgraph.SetPoint(2, 819.5, 645);
        cgraph.SetPoint(3, 819, 655);
        cgraph.SetPoint(4, 818.5, 665);
        cgraph.SetPoint(5, 818, 675);
        cgraph.SetPoint(6, 816, 700);
        cgraph.SetPoint(7, 815, 710);
        cgraph.SetPoint(8, 813.75, 712.5);
        cgraph.SetPoint(9, 812.5, 715);
        cgraph.SetPoint(10, 810, 720);
        cgraph.SetPoint(11, 805, 725);
        cgraph.SetPoint(12, 800, 730);
        cgraph.SetPoint(13, 790, 740);
        cgraph.SetPoint(14, 780, 750);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(9);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 800, ymin);
        dgraph.SetPoint(1, 800, 660);
        dgraph.SetPoint(2, 799.75, 661.25);
        dgraph.SetPoint(3, 799.5, 662.5);
        dgraph.SetPoint(4, 799, 665);
        dgraph.SetPoint(5, 797.5, 670);
        dgraph.SetPoint(6, 795, 680);
        dgraph.SetPoint(7, 790, 700);
        dgraph.SetPoint(8, 780, 730);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(14);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 835, ymin);
        ugraph.SetPoint(1, 835, 660);
        ugraph.SetPoint(2, 834.5, 667.5);
        ugraph.SetPoint(3, 834, 675);
        ugraph.SetPoint(4, 833, 690);
        ugraph.SetPoint(5, 832.5, 697.5);
        ugraph.SetPoint(6, 832, 705);
        ugraph.SetPoint(7, 831.5, 712.5);
        ugraph.SetPoint(8, 831, 720);
        ugraph.SetPoint(9, 830, 730);
        ugraph.SetPoint(10, 825, 740);
        ugraph.SetPoint(11, 815, 750);
        ugraph.SetPoint(12, 805, 760);
        ugraph.SetPoint(13, 795, 770);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c11->Print(Form("B2_ExcludedRegionMap_%d.pdf", chi_mass));        
    }

    // Draw the limit lines on top of the temperature plot
    TCanvas* c12 = new TCanvas();
    c12->SetFillColor(0);
    c12->SetFillColor(0);
    c12->SetBorderMode(0);
    c12->GetPad(0)->SetBorderSize(2);
    c12->GetPad(0)->SetLeftMargin(0.1407035);
    c12->GetPad(0)->SetTopMargin(0.08);
    c12->GetPad(0)->SetBottomMargin(0.13);

    ul->Draw("colz");
    if (chi_mass == 150) {
        TGraph cgraph = TGraph(17);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 800, 400);
        cgraph.SetPoint(2, 820, 460);
        cgraph.SetPoint(3, 835, 580);
        cgraph.SetPoint(4, 837, 590);
        cgraph.SetPoint(5, 838, 595);
        cgraph.SetPoint(6, 839, 600);                
        cgraph.SetPoint(7, 840, 605);
        cgraph.SetPoint(8, 841, 608);
        cgraph.SetPoint(9, 842, 615);
        cgraph.SetPoint(10, 843, 620);
        cgraph.SetPoint(11, 844, 625);
        cgraph.SetPoint(12, 845, 630);
        cgraph.SetPoint(13, 846, 650);
        cgraph.SetPoint(14, 848, 680);
        cgraph.SetPoint(15, 850, 720);
        cgraph.SetPoint(16, 840, 775);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(17);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 780, 400);
        dgraph.SetPoint(2, 800, 460);
        dgraph.SetPoint(3, 815, 580);
        dgraph.SetPoint(4, 817, 590);
        dgraph.SetPoint(5, 818, 595);
        dgraph.SetPoint(6, 819, 600);                
        dgraph.SetPoint(7, 820, 605);
        dgraph.SetPoint(8, 821, 608);
        dgraph.SetPoint(9, 822, 615);
        dgraph.SetPoint(10, 823, 620);
        dgraph.SetPoint(11, 824, 625);
        dgraph.SetPoint(12, 825, 630);
        dgraph.SetPoint(13, 826, 650);
        dgraph.SetPoint(14, 828, 680);
        dgraph.SetPoint(15, 830, 710);
        dgraph.SetPoint(16, 825, 775);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(17);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 770, ymin);
        ugraph.SetPoint(1, 810, 400);
        ugraph.SetPoint(2, 835, 460);
        ugraph.SetPoint(3, 855, 580);
        ugraph.SetPoint(4, 857, 590);
        ugraph.SetPoint(5, 858, 595);
        ugraph.SetPoint(6, 859, 600);                
        ugraph.SetPoint(7, 860, 605);
        ugraph.SetPoint(8, 861, 608);
        ugraph.SetPoint(9, 862, 615);
        ugraph.SetPoint(10, 863, 620);
        ugraph.SetPoint(11, 864, 625);
        ugraph.SetPoint(12, 865, 630);
        ugraph.SetPoint(13, 866, 650);
        ugraph.SetPoint(14, 868, 680);
        ugraph.SetPoint(15, 870, 730);
        ugraph.SetPoint(16, 860, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("B2_LimitsOnCarpet_%d.pdf", chi_mass));

	// Do it again but using the same limits that we used for the paper
	float x150   = 150. + mtop + mb;
	float y150c  = 150. + mtop;
	float y150m  = xmax - mb;
	TPolyLine *p150 = new TPolyLine();
	p150->SetLineColor(1);
	p150->SetLineWidth(1);
	p150->SetLineStyle(1);
	// p150->SetFillStyle(3395);
	p150->SetFillStyle(3305);
	p150->SetFillColor(1);
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
	TPolyLine *blah = new TPolyLine();
	blah->SetLineColor(1);
	blah->SetLineWidth(1);
	blah->SetLineStyle(1);
	blah->SetFillStyle(3305);
	blah->SetFillColor(1);
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.20*(ymax-ymin));
	blah->SetNextPoint(xmin+0.18*(xmax-xmin), ymax-0.20*(ymax-ymin));
	blah->SetNextPoint(xmin+0.18*(xmax-xmin), ymax-0.25*(ymax-ymin));
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
	TCanvas* c12n = new TCanvas();
	c12n->SetFillColor(0);
	c12n->SetFillColor(0);
	c12n->SetBorderMode(0);
	c12n->GetPad(0)->SetBorderSize(2);
	c12n->GetPad(0)->SetLeftMargin(0.1407035);
	c12n->GetPad(0)->SetTopMargin(0.08);
	c12n->GetPad(0)->SetBottomMargin(0.13);	
	ul->Draw("colz");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
        //l2.Draw();
        //l1.Draw();
	p150->Draw("fl");
	p150->Draw();
	blah->Draw("fl");
	blah->Draw();
        c12n->Print(Form("B2_LimitsOnCarpet_LikePaper_%d.pdf", chi_mass));
    }
    else if (chi_mass == 200) {
        TGraph cgraph = TGraph(18);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 805, 420);
        cgraph.SetPoint(2, 820, 450);
        cgraph.SetPoint(3, 840, 540);
        cgraph.SetPoint(4, 845, 600);
        cgraph.SetPoint(5, 846, 620);
        cgraph.SetPoint(6, 847, 640);
        cgraph.SetPoint(7, 848, 660);
        cgraph.SetPoint(8, 849, 680);
        cgraph.SetPoint(9, 849, 735);
        cgraph.SetPoint(10, 848.75, 736.25);
        cgraph.SetPoint(11, 848.5, 737.5);
        cgraph.SetPoint(12, 848, 740);
        cgraph.SetPoint(13, 847.5, 742.5);
        cgraph.SetPoint(14, 847, 745);
        cgraph.SetPoint(15, 846, 750);
        cgraph.SetPoint(16, 845, 755);
        cgraph.SetPoint(17, 840, 780);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(13);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 770, 420);
        dgraph.SetPoint(2, 800, 450);
        dgraph.SetPoint(3, 820, 540);
        dgraph.SetPoint(4, 825, 600);
        dgraph.SetPoint(5, 826, 620);
        dgraph.SetPoint(6, 827, 640);
        dgraph.SetPoint(7, 828, 660);
        dgraph.SetPoint(8, 829, 680);
        dgraph.SetPoint(9, 829, 710);
        dgraph.SetPoint(10, 830, 725);
        dgraph.SetPoint(11, 831, 740);
        dgraph.SetPoint(12, 820, 775);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(19);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 780, ymin);
        ugraph.SetPoint(1, 825, 415);
        ugraph.SetPoint(2, 835, 425);
        ugraph.SetPoint(3, 845, 455);
        ugraph.SetPoint(4, 850, 470);
        ugraph.SetPoint(5, 855, 485);
        ugraph.SetPoint(6, 857, 505);
        ugraph.SetPoint(7, 860, 530);
        ugraph.SetPoint(8, 861.5, 550);
        ugraph.SetPoint(9, 863, 570);
        ugraph.SetPoint(10, 864, 585);
        ugraph.SetPoint(11, 865, 600);
        ugraph.SetPoint(12, 865, 740);
        ugraph.SetPoint(13, 864.5, 742);
        ugraph.SetPoint(14, 863.75, 743.75);
        ugraph.SetPoint(15, 862.5, 747.5);
        ugraph.SetPoint(16, 860, 755);
        ugraph.SetPoint(17, 855, 770);
        ugraph.SetPoint(18, 850, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("B2_LimitsOnCarpet_%d.pdf", chi_mass));
    }
    else if (chi_mass == 300) {
        TGraph cgraph = TGraph(15);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 820, ymin);
        cgraph.SetPoint(1, 820, 635);
        cgraph.SetPoint(2, 819.5, 645);
        cgraph.SetPoint(3, 819, 655);
        cgraph.SetPoint(4, 818.5, 665);
        cgraph.SetPoint(5, 818, 675);
        cgraph.SetPoint(6, 816, 700);
        cgraph.SetPoint(7, 815, 710);
        cgraph.SetPoint(8, 813.75, 712.5);
        cgraph.SetPoint(9, 812.5, 715);
        cgraph.SetPoint(10, 810, 720);
        cgraph.SetPoint(11, 805, 725);
        cgraph.SetPoint(12, 800, 730);
        cgraph.SetPoint(13, 790, 740);
        cgraph.SetPoint(14, 780, 750);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(9);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 800, ymin);
        dgraph.SetPoint(1, 800, 660);
        dgraph.SetPoint(2, 799.75, 661.25);
        dgraph.SetPoint(3, 799.5, 662.5);
        dgraph.SetPoint(4, 799, 665);
        dgraph.SetPoint(5, 797.5, 670);
        dgraph.SetPoint(6, 795, 680);
        dgraph.SetPoint(7, 790, 700);
        dgraph.SetPoint(8, 780, 730);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(14);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 835, ymin);
        ugraph.SetPoint(1, 835, 660);
        ugraph.SetPoint(2, 834.5, 667.5);
        ugraph.SetPoint(3, 834, 675);
        ugraph.SetPoint(4, 833, 690);
        ugraph.SetPoint(5, 832.5, 697.5);
        ugraph.SetPoint(6, 832, 705);
        ugraph.SetPoint(7, 831.5, 712.5);
        ugraph.SetPoint(8, 831, 720);
        ugraph.SetPoint(9, 830, 730);
        ugraph.SetPoint(10, 825, 740);
        ugraph.SetPoint(11, 815, 750);
        ugraph.SetPoint(12, 805, 760);
        ugraph.SetPoint(13, 795, 770);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c12->Print(Form("B2_LimitsOnCarpet_%d.pdf", chi_mass));

	// Do it again but using the same limits that we used for the paper
	float x300   = 300. + mtop + mb;
	float y300c  = 300. + mtop;
	float y300m  = xmax - mb;

	TPolyLine *p300 = new TPolyLine();
	p300->SetLineColor(1);
	p300->SetFillStyle(3305);
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
	TPolyLine *blah = new TPolyLine();
	blah->SetLineColor(1);
	blah->SetLineWidth(1);
	blah->SetLineStyle(1);
	blah->SetFillStyle(3305);
	blah->SetFillColor(1);
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.20*(ymax-ymin));
	blah->SetNextPoint(xmin+0.18*(xmax-xmin), ymax-0.20*(ymax-ymin));
	blah->SetNextPoint(xmin+0.18*(xmax-xmin), ymax-0.25*(ymax-ymin));
	blah->SetNextPoint(xmin+0.1*(xmax-xmin), ymax-0.25*(ymax-ymin));
	TCanvas* c12nn = new TCanvas();
	c12nn->SetFillColor(0);
	c12nn->SetFillColor(0);
	c12nn->SetBorderMode(0);
	c12nn->GetPad(0)->SetBorderSize(2);
	c12nn->GetPad(0)->SetLeftMargin(0.1407035);
	c12nn->GetPad(0)->SetTopMargin(0.08);
	c12nn->GetPad(0)->SetBottomMargin(0.13);	
	ul->Draw("colz");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
        //l2.Draw();
        //l1.Draw();
	p300->Draw("fl");
	p300->Draw();
	blah->Draw("fl");
	blah->Draw();
        c12nn->Print(Form("B2_LimitsOnCarpet_LikePaper_%d.pdf", chi_mass));
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


    if (chi_mass == 150) {
        TGraph cgraph = TGraph(18);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 805, 420);
        cgraph.SetPoint(2, 820, 450);
        cgraph.SetPoint(3, 840, 540);
        cgraph.SetPoint(4, 845, 600);
        cgraph.SetPoint(5, 846, 620);
        cgraph.SetPoint(6, 847, 640);
        cgraph.SetPoint(7, 848, 660);
        cgraph.SetPoint(8, 849, 680);
        cgraph.SetPoint(9, 849, 735);
        cgraph.SetPoint(10, 848.75, 736.25);
        cgraph.SetPoint(11, 848.5, 737.5);
        cgraph.SetPoint(12, 848, 740);
        cgraph.SetPoint(13, 847.5, 742.5);
        cgraph.SetPoint(14, 847, 745);
        cgraph.SetPoint(15, 846, 750);
        cgraph.SetPoint(16, 845, 755);
        cgraph.SetPoint(17, 840, 780);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(17);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 780, 400);
        dgraph.SetPoint(2, 800, 460);
        dgraph.SetPoint(3, 815, 580);
        dgraph.SetPoint(4, 817, 590);
        dgraph.SetPoint(5, 818, 595);
        dgraph.SetPoint(6, 819, 600);                
        dgraph.SetPoint(7, 820, 605);
        dgraph.SetPoint(8, 821, 608);
        dgraph.SetPoint(9, 822, 615);
        dgraph.SetPoint(10, 823, 620);
        dgraph.SetPoint(11, 824, 625);
        dgraph.SetPoint(12, 825, 630);
        dgraph.SetPoint(13, 826, 650);
        dgraph.SetPoint(14, 828, 680);
        dgraph.SetPoint(15, 830, 710);
        dgraph.SetPoint(16, 825, 775);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(17);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 770, ymin);
        ugraph.SetPoint(1, 810, 400);
        ugraph.SetPoint(2, 835, 460);
        ugraph.SetPoint(3, 855, 580);
        ugraph.SetPoint(4, 857, 590);
        ugraph.SetPoint(5, 858, 595);
        ugraph.SetPoint(6, 859, 600);                
        ugraph.SetPoint(7, 860, 605);
        ugraph.SetPoint(8, 861, 608);
        ugraph.SetPoint(9, 862, 615);
        ugraph.SetPoint(10, 863, 620);
        ugraph.SetPoint(11, 864, 625);
        ugraph.SetPoint(12, 865, 630);
        ugraph.SetPoint(13, 866, 650);
        ugraph.SetPoint(14, 868, 680);
        ugraph.SetPoint(15, 870, 730);
        ugraph.SetPoint(16, 860, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("B2_LimitsOnWhite_%d.pdf", chi_mass));
    }
    else if (chi_mass == 200) {
        TGraph cgraph = TGraph(10);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 760, ymin);
        cgraph.SetPoint(1, 805, 420);
        cgraph.SetPoint(2, 820, 450);
        cgraph.SetPoint(3, 840, 540);
        cgraph.SetPoint(4, 845, 600);
        cgraph.SetPoint(5, 846, 620);
        cgraph.SetPoint(6, 847, 640);
        cgraph.SetPoint(7, 848, 660);
        cgraph.SetPoint(8, 849, 680);
        cgraph.SetPoint(9, 849, 735);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(12);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 740, ymin);
        dgraph.SetPoint(1, 770, 420);
        dgraph.SetPoint(2, 800, 450);
        dgraph.SetPoint(3, 820, 540);
        dgraph.SetPoint(4, 825, 600);
        dgraph.SetPoint(5, 826, 620);
        dgraph.SetPoint(6, 827, 640);
        dgraph.SetPoint(7, 828, 660);
        dgraph.SetPoint(8, 829, 680);
        dgraph.SetPoint(9, 829, 710);
        dgraph.SetPoint(10, 830, 725);
        dgraph.SetPoint(11, 831, 740);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(19);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 780, ymin);
        ugraph.SetPoint(1, 825, 415);
        ugraph.SetPoint(2, 835, 425);
        ugraph.SetPoint(3, 845, 455);
        ugraph.SetPoint(4, 850, 470);
        ugraph.SetPoint(5, 855, 485);
        ugraph.SetPoint(6, 857, 505);
        ugraph.SetPoint(7, 860, 530);
        ugraph.SetPoint(8, 861.5, 550);
        ugraph.SetPoint(9, 863, 570);
        ugraph.SetPoint(10, 864, 585);
        ugraph.SetPoint(11, 865, 600);
        ugraph.SetPoint(12, 865, 740);
        ugraph.SetPoint(13, 864.5, 742);
        ugraph.SetPoint(14, 863.75, 743.75);
        ugraph.SetPoint(15, 862.5, 747.5);
        ugraph.SetPoint(16, 860, 755);
        ugraph.SetPoint(17, 855, 770);
        ugraph.SetPoint(18, 850, 795);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("B2_LimitsOnWhite_%d.pdf", chi_mass));
    }
    else if (chi_mass == 300) {
        TGraph cgraph = TGraph(15);
        cgraph.SetLineColor(1);
        cgraph.SetLineWidth(3);
        cgraph.SetPoint(0, 820, ymin);
        cgraph.SetPoint(1, 820, 635);
        cgraph.SetPoint(2, 819.5, 645);
        cgraph.SetPoint(3, 819, 655);
        cgraph.SetPoint(4, 818.5, 665);
        cgraph.SetPoint(5, 818, 675);
        cgraph.SetPoint(6, 816, 700);
        cgraph.SetPoint(7, 815, 710);
        cgraph.SetPoint(8, 813.75, 712.5);
        cgraph.SetPoint(9, 812.5, 715);
        cgraph.SetPoint(10, 810, 720);
        cgraph.SetPoint(11, 805, 725);
        cgraph.SetPoint(12, 800, 730);
        cgraph.SetPoint(13, 790, 740);
        cgraph.SetPoint(14, 780, 750);
        cgraph.Draw("sameC");
        TGraph dgraph = TGraph(9);
        dgraph.SetLineColor(1);
        dgraph.SetLineWidth(3);
        dgraph.SetLineStyle(2);
        dgraph.SetPoint(0, 800, ymin);
        dgraph.SetPoint(1, 800, 660);
        dgraph.SetPoint(2, 799.75, 661.25);
        dgraph.SetPoint(3, 799.5, 662.5);
        dgraph.SetPoint(4, 799, 665);
        dgraph.SetPoint(5, 797.5, 670);
        dgraph.SetPoint(6, 795, 680);
        dgraph.SetPoint(7, 790, 700);
        dgraph.SetPoint(8, 780, 730);
        dgraph.Draw("sameC");
        TGraph ugraph = TGraph(14);
        ugraph.SetLineColor(1);
        ugraph.SetLineWidth(3);
        ugraph.SetLineStyle(2);
        ugraph.SetPoint(0, 835, ymin);
        ugraph.SetPoint(1, 835, 660);
        ugraph.SetPoint(2, 834.5, 667.5);
        ugraph.SetPoint(3, 834, 675);
        ugraph.SetPoint(4, 833, 690);
        ugraph.SetPoint(5, 832.5, 697.5);
        ugraph.SetPoint(6, 832, 705);
        ugraph.SetPoint(7, 831.5, 712.5);
        ugraph.SetPoint(8, 831, 720);
        ugraph.SetPoint(9, 830, 730);
        ugraph.SetPoint(10, 825, 740);
        ugraph.SetPoint(11, 815, 750);
        ugraph.SetPoint(12, 805, 760);
        ugraph.SetPoint(13, 795, 770);
        ugraph.Draw("sameC");
        kinlim.Draw();
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
        latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
        gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), central_text);
        gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.32*(ymax-ymin), bands_text);
        l2.Draw();
        l1.Draw();
        c13->Print(Form("B2_LimitsOnWhite_%d.pdf", chi_mass));
    }
   
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

    ulbest->Draw("textcol");
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
    c14->Print(Form("B2_BestSignalRegion_%d.pdf", chi_mass));

    //Draw the acceptance carpet
    TCanvas* c15 = new TCanvas();
    c15->SetFillColor(0);
    c15->SetFillColor(0);
    c15->SetBorderMode(0);
    c15->GetPad(0)->SetBorderSize(2);
    c15->GetPad(0)->SetLeftMargin(0.1407035);
    c15->GetPad(0)->SetTopMargin(0.08);
    c15->GetPad(0)->SetBottomMargin(0.13);

    acc->Draw("colz");

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.16*(ymax-ymin), masses);
    kinlim.Draw();
    c15->Print(Form("B2_AcceptanceCarpet_%d.pdf", chi_mass));
}
