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

void makePlots () {

    TFile *file = TFile::Open("ntuple.root");
    file->cd();
    TTree *tree = (TTree*)gDirectory->Get("tree");

    // Makes basic histograms from the root tree
    gStyle->SetPadRightMargin(0.16);   // default 0.1
    gStyle->SetTitleOffset(1.20, "Y");  // default 1
    gStyle->SetOptStat(0);
    gStyle->SetTitleOffset(1.20, "Z");  // default 1

    // This binning insures that the bins are nicely centered
    Double_t xbinsize = 25.;
    Double_t ybinsize = 25.;
    Double_t xmin     = 400. - xbinsize/2.;
    Double_t xmax     = 1100 + xbinsize/2.;
    Double_t ymin     = 50. - ybinsize/2;
    Double_t ymax     = 800. + ybinsize/2.;
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


    tree->Draw("lspmass:glmass>>ul","explimsrb/(5000.*effsrb)");
    tree->Draw("lspmass:glmass>>ulbest","bestsr");
    tree->Draw("lspmass:glmass>>acc","effsrb");
    tree->Draw("lspmass:glmass>>excl","explimsrb/(5000.*effsrb)<xsec");
    tree->Draw("lspmass:glmass>>exclup","explimsrb/(5000.*effsrb)<xsecup");
    tree->Draw("lspmass:glmass>>excldwn","explimsrb/(5000.*effsrb)<xsecdwn");
    tree->Draw("lspmass:glmass>>hxsec",   "xsec");
    tree->Draw("lspmass:glmass>>hxsecup", "xsecup");
    tree->Draw("lspmass:glmass>>hxsecdwn","xsecdwn");

    // There is a vertical slice of missing points.  Fill the acceptance and 
    // cross-section limits by interpolation
    for (int iy=1; iy<=ny; iy++) {
        for (int ix=2; ix<ny; ix++) {
            if (hxsec->GetBinContent(ix,iy) == 0 && hxsec->GetBinContent(ix-1,iy) > 0) {
                float x = xmin + (ix-0.5)*xbinsize;
                float y = ymin + (iy-0.5)*xbinsize;
                cout << "Interpolating at x,y = " << x << " " << y << endl;
                acc->Fill(x,y,0.5*(acc->GetBinContent(ix-1,iy) + acc->GetBinContent(ix+1,iy)));
                ul->Fill(x,y,0.5*(ul->GetBinContent(ix-1,iy) + ul->GetBinContent(ix+1,iy)));
            }
        }
    }

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
                        // cout << xlim << " " << y << endl;
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
    TLine kinlim = TLine(2.*mtop+ymin, ymin, xmax, xmax-2*mtop);
    kinlim.SetLineWidth(3);

    char * selection ="Same Sign dileptons with btag selection";
    const char *obligatory_text = "CMS, #sqrt{s} = 7 TeV, L_{int} = 4.98 fb^{-1}";
    const char *central_text    = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL}";
    const char *bands_text      = "Exclusion #sigma^{prod} = #sigma^{NLO+NLL} #pm 1 #sigma";
    char* ztitle  = "#sigma_{UL} pb";

    // Axis labels, a bit primitive for now
    char* glmass = "m(#tilde{g}) GeV";
    char* lspmass = "m(#tilde{#chi}_{1}^{0}) GeV";
    excl->GetYaxis()->SetTitle(lspmass);
    excl->GetXaxis()->SetTitle(glmass);
    excl->SetTitle("T1 model  Excluded points in red");
    ul->GetYaxis()->SetTitle(lspmass);
    ul->GetXaxis()->SetTitle(glmass);
    ul->GetZaxis()->SetTitle(ztitle);
    ul->SetTitle("T1 model  Cross-section upper limits (pb)");
    empty->GetYaxis()->SetTitle(lspmass);
    empty->GetXaxis()->SetTitle(glmass);
    empty->SetTitle("T1 model");
    ulbest->GetYaxis()->SetTitle(lspmass);
    ulbest->GetXaxis()->SetTitle(glmass);
    ulbest->SetTitle("T1 model  Best region based on exp. limit");
    acc->GetYaxis()->SetTitle(lspmass);
    acc->GetXaxis()->SetTitle(glmass);
    acc->SetTitle("T1 model  Acc*Eff*BR for best signal region in percent");

    // Some text
    TLatex gg;
    gg.SetTextSize(0.035);

    TLatex gg2;
    gg2.SetTextSize(0.035);

    TLatex latexLabel;
    latexLabel.SetTextSize(0.035);

    TLine l1 = TLine(xmin+0.1*(xmax-xmin), ymax-0.14*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.14*(ymax-ymin));
    l1.SetLineColor(4);
    l1.SetLineWidth(3);

    TLine l2 = TLine(xmin+0.1*(xmax-xmin), ymax-0.22*(ymax-ymin), xmin+0.2*(xmax-xmin), ymax-0.22*(ymax-ymin));
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

    excl->Draw("col");
    TGraph cgraph = TGraph(6);
    cgraph.SetLineColor(4);
    cgraph.SetLineWidth(3);
    cgraph.SetPoint(0, 830, ymin);
    cgraph.SetPoint(1, 825, 130);
    cgraph.SetPoint(2, 820, 160);
    cgraph.SetPoint(3, 810, 240);
    cgraph.SetPoint(4, 800, 280);
    cgraph.SetPoint(5, 750, 345);
    cgraph.Draw("sameC");
    TGraph dgraph = TGraph(7);
    dgraph.SetLineColor(4);
    dgraph.SetLineWidth(3);
    dgraph.SetLineStyle(2);
    dgraph.SetPoint(0, 810, ymin);
    dgraph.SetPoint(1, 805, 130);
    dgraph.SetPoint(2, 800, 160);
    dgraph.SetPoint(3, 780, 235);
    dgraph.SetPoint(4, 770, 270);
    dgraph.SetPoint(5, 750, 315);
    dgraph.SetPoint(6, 730, 325);
    dgraph.Draw("sameC");
    TGraph ugraph = TGraph(7);
    ugraph.SetLineColor(4);
    ugraph.SetLineWidth(3);
    ugraph.SetLineStyle(2);
    ugraph.SetPoint(0, 860, ymin);
    ugraph.SetPoint(1, 855, 130);
    ugraph.SetPoint(2, 850, 160);
    ugraph.SetPoint(3, 830, 250);
    ugraph.SetPoint(4, 820, 295);
    ugraph.SetPoint(5, 800, 340);
    ugraph.SetPoint(6, 775, 370);
    ugraph.Draw("sameC");
    
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.16*(ymax-ymin), central_text);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
    l2.Draw();
    l1.Draw();
    c11->Print("T1tttt_ExcludedRegionMap.pdf");


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
    cgraph = TGraph(6);
    cgraph.SetLineColor(4);
    cgraph.SetLineWidth(3);
    cgraph.SetPoint(0, 830, ymin);
    cgraph.SetPoint(1, 825, 130);
    cgraph.SetPoint(2, 820, 160);
    cgraph.SetPoint(3, 810, 240);
    cgraph.SetPoint(4, 800, 280);
    cgraph.SetPoint(5, 750, 345);
    cgraph.Draw("sameC");
    dgraph = TGraph(7);
    dgraph.SetLineColor(4);
    dgraph.SetLineWidth(3);
    dgraph.SetLineStyle(2);
    dgraph.SetPoint(0, 810, ymin);
    dgraph.SetPoint(1, 805, 130);
    dgraph.SetPoint(2, 800, 160);
    dgraph.SetPoint(3, 780, 235);
    dgraph.SetPoint(4, 770, 270);
    dgraph.SetPoint(5, 750, 315);
    dgraph.SetPoint(6, 730, 325);
    dgraph.Draw("sameC");
    ugraph = TGraph(7);
    ugraph.SetLineColor(4);
    ugraph.SetLineWidth(3);
    ugraph.SetLineStyle(2);
    ugraph.SetPoint(0, 860, ymin);
    ugraph.SetPoint(1, 855, 130);
    ugraph.SetPoint(2, 850, 160);
    ugraph.SetPoint(3, 830, 250);
    ugraph.SetPoint(4, 820, 295);
    ugraph.SetPoint(5, 800, 340);
    ugraph.SetPoint(6, 775, 370);
    ugraph.Draw("sameC");
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.16*(ymax-ymin), central_text);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
    l2.Draw();
    l1.Draw();
    c12->Print("T1tttt_LimitsOnCarpet.pdf");

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

    cgraph = TGraph(6);
    cgraph.SetLineColor(4);
    cgraph.SetLineWidth(3);
    cgraph.SetPoint(0, 830, ymin);
    cgraph.SetPoint(1, 825, 130);
    cgraph.SetPoint(2, 820, 160);
    cgraph.SetPoint(3, 810, 240);
    cgraph.SetPoint(4, 800, 280);
    cgraph.SetPoint(5, 750, 345);
    cgraph.Draw("sameC");
    dgraph = TGraph(7);
    dgraph.SetLineColor(4);
    dgraph.SetLineWidth(3);
    dgraph.SetLineStyle(2);
    dgraph.SetPoint(0, 810, ymin);
    dgraph.SetPoint(1, 805, 130);
    dgraph.SetPoint(2, 800, 160);
    dgraph.SetPoint(3, 780, 235);
    dgraph.SetPoint(4, 770, 270);
    dgraph.SetPoint(5, 750, 315);
    dgraph.SetPoint(6, 730, 325);
    dgraph.Draw("sameC");
    ugraph = TGraph(7);
    ugraph.SetLineColor(4);
    ugraph.SetLineWidth(3);
    ugraph.SetLineStyle(2);
    ugraph.SetPoint(0, 860, ymin);
    ugraph.SetPoint(1, 855, 130);
    ugraph.SetPoint(2, 850, 160);
    ugraph.SetPoint(3, 830, 250);
    ugraph.SetPoint(4, 820, 295);
    ugraph.SetPoint(5, 800, 340);
    ugraph.SetPoint(6, 775, 370);
    ugraph.Draw("sameC");
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    gg.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.16*(ymax-ymin), central_text);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.24*(ymax-ymin), bands_text);
    l2.Draw();
    l1.Draw();
    c13->Print("T1tttt_LimitsOnWhite.pdf");
   
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

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    kinlim.Draw();
    c14->Print("T1tttt_BestSignalRegion.pdf");

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
    kinlim.Draw();
    c15->Print("T1tttt_AcceptanceCarpet.pdf");

    //Draw the smoothed limit lines and nothing else
    TCanvas* c16 = new TCanvas();
    c16->SetFillColor(0);
    c16->GetPad(0)->SetRightMargin(0.07);
    c16->SetFillColor(0);
    c16->SetBorderMode(0);
    c16->GetPad(0)->SetBorderSize(2);
    c16->GetPad(0)->SetLeftMargin(0.1407035);
    c16->GetPad(0)->SetTopMargin(0.08);
    c16->GetPad(0)->SetBottomMargin(0.13);

    empty->Draw();
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.16*(ymax-ymin), bands_text);  

    // A polyline for the legend
    TPolyLine *lg = new TPolyLine();
    lg->SetLineColor(kBlue);
    lg->SetFillStyle(3244);
    lg->SetFillColor(kBlue);
    lg->SetNextPoint(xmin+0.1*(xmax-xmin),ymax-0.16*(ymax-ymin));
    lg->SetNextPoint(xmin+0.2*(xmax-xmin),ymax-0.16*(ymax-ymin));
    lg->SetNextPoint(xmin+0.2*(xmax-xmin),ymax-0.13*(ymax-ymin));
    lg->SetNextPoint(xmin+0.1*(xmax-xmin),ymax-0.13*(ymax-ymin));
    lg->SetNextPoint(xmin+0.1*(xmax-xmin),ymax-0.16*(ymax-ymin));
    lg->Draw("fl");

    // A polyline with the smoothed limit
    TPolyLine *psm = new TPolyLine();
    psm->SetLineColor(kBlue);
    psm->SetFillStyle(3244);
    psm->SetFillColor(kBlue);
    psm->SetNextPoint(805.,ymin);
    psm->SetNextPoint(805.,50.);
    psm->SetNextPoint(790.,180.);
    psm->SetNextPoint(780.,250.);
    psm->SetNextPoint(740.,740.-2*mtop);
    psm->SetNextPoint(770.,770.-2*mtop);
    psm->SetNextPoint(815.,300.);
    psm->SetNextPoint(840.,200.);
    psm->SetNextPoint(855.,50.);
    psm->SetNextPoint(860.,ymin);
    psm->SetNextPoint(805.,ymin);
    psm->Draw("fl");
    c16->Print("T1tttt_SmoothLimitsOnWhite.pdf");


    //Draw the smoothed limit lines with the upper limit carpet
    TCanvas* c17 = new TCanvas();

    c17->SetFillColor(0);
    c17->SetFillColor(0);
    c17->SetBorderMode(0);
    c17->GetPad(0)->SetBorderSize(2);
    c17->GetPad(0)->SetLeftMargin(0.1407035);
    c17->GetPad(0)->SetTopMargin(0.08);
    c17->GetPad(0)->SetBottomMargin(0.13);

    ul->Draw("colz");

    // empty->Draw();
    kinlim.Draw();

    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax+0.04*(ymax-ymin), obligatory_text);
    latexLabel.DrawLatex(xmin+0.1*(xmax-xmin), ymax-0.08*(ymax-ymin),selection);
    gg2.DrawLatex(xmin+0.21*(xmax-xmin), ymax-0.16*(ymax-ymin), bands_text);  

    // A polyline for the legend
    lg->SetLineColor(1);
    lg->SetFillStyle(3344);
    lg->SetFillColor(1);
    lg->Draw("fl");
    lg->Draw();

    // A polyline with the smoothed limit
    psm->SetLineColor(1);
    psm->SetFillStyle(3244);
    psm->SetFillColor(1);
    psm->Draw("fl");
    psm->Draw();


    c17->Print("T1tttt_LimitsOnCarpetLikePaper.pdf");


}
