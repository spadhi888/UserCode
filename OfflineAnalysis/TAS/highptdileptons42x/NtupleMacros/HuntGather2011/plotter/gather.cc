#include "BabySample.h"
#include "cuts.h"
#include "gather.h"
#include "../../Tools/goodrun.h"
#include "BabyDorkIdentifier.h"

#include "TPad.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TROOT.h"
#include "TTreePlayer.h"

#include <algorithm>
#include <vector>
#include <iostream>

float GetAllLumi(TChain *c, float zPerPb)
{   

    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));
    TCut c_remove_end2010bad("remove_end2010bad", "run <= 149294 || run >= 160325");

    reset_babydorkidentifier();
    int n_new = c->GetEntries(c_goodrunplus+c_remove_end2010bad+inclusivez_dilep);
    reset_babydorkidentifier();
    float lumi = ((float)(n_new))/zPerPb;
    return lumi;

}


float GetNewLumi(TChain *c, float zPerPb)
{

    int erun = max_run();
    int els  = max_run_max_lumi();
    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_newrun(Form("(run>%i||(run==%i&&ls>%i))", erun, erun, els));
    TCut c_remove_end2010bad("remove_end2010bad", "run <= 149294 || run >= 160325");

    reset_babydorkidentifier();
    int n_new = c->GetEntries(c_newrun+c_remove_end2010bad+inclusivez_dilep);
    reset_babydorkidentifier();
    float lumi = ((float)(n_new))/zPerPb;
    return lumi;

}

float GetZPerPb(TChain *c, float lumi)
{

    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    TCut c_goodrun(Form("((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls)", brun, erun, brun, bls, erun, els));

    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized

    // brun:bls -> erun:els
    reset_babydorkidentifier();
    int n_goodrun = c->GetEntries(c_goodrun+inclusivez_dilep);
    // total
    reset_babydorkidentifier();
    // that which is new

    float zPerPb = ((float)(n_goodrun))/lumi;
    std::cout << "in the good run list for " << lumi << " there are " << n_goodrun/lumi << " Z/pb" << std::endl;

    return zPerPb;
}


float GetIntLumi(TChain *c, float lumi, int brun, int bls, int erun, int els)
{

    TCut c_goodrun(Form("((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls)", brun, erun, brun, bls, erun, els));

    // goodrun plus events beyond range of goodrun
    // which are not goodrun penalized
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun_json(run,ls))||(run>%i||(run==%i&&ls>%i))", brun, erun, brun, bls, erun, els, erun, erun, els));
    TCut c_remove_end2010bad("remove_end2010bad", "run <= 149294 || run >= 160325");

    // brun:bls -> erun:els
    reset_babydorkidentifier();
    int n_goodrun = c->GetEntries(c_goodrun+inclusivez_dilep); 
    // total
    reset_babydorkidentifier();
    int n_total   = c->GetEntries(c_goodrunplus+c_remove_end2010bad+inclusivez_dilep);
    // that which is new
    int n_new     = n_total-n_goodrun;

    float newlumi = ((float)(n_new*lumi))/(float)n_goodrun;

    std::cout << "in the good run list for " << lumi << " there are " << n_goodrun << std::endl;
    std::cout << "out of the good run list there are " << n_new << std::endl;
    std::cout << "that means a new lumi of " << newlumi << std::endl;

    return lumi+newlumi;
}

float GetIntLumi(BabySample *bs, float lumi)
{
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    std::cout << brun << ", " << bls << " : " << erun << ", " << els << std::endl;
    return GetIntLumi(bs->chain(), lumi, brun, bls, erun, els);
}

TCanvas* TagAndProbe(const char *savename, TCut evsel, TCut var1, TCut var2,
        TCut tag1, TCut tag2, TCut probe1, TCut probe2, TCut sel1, TCut sel2, 
        float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated, std::vector<BabySample*> bss)
{

    std::vector<TH1F*> vh_mc_2TT;
    std::vector<TH1F*> vh_mc_TP;
    std::vector<TH1F*> vh_mc_TF;
    std::vector<TH1F*> vh_data_2TT;
    std::vector<TH1F*> vh_data_TP;
    std::vector<TH1F*> vh_data_TF;
    TH1F* h1_2TT_1;
    TH1F* h1_TP_1;
    TH1F* h1_TF_1;
    TH1F* h1_2TT_2;
    TH1F* h1_TP_2;
    TH1F* h1_TF_2;
    //TH1F* h1_2TT;
    //TH1F* h1_TP;
    //TH1F* h1_TF;

    //
    // define the categories
    //

    // remember we need to plot the variable for the probe
    // when the other leg is the tag... so this means a few
    // permutations have to be defined

    // 2TT means both legs are tags
    TCut cut_2TT = tag1 && tag2 && evsel;
    cut_2TT.SetName(TString("2TT"));
    // TP means one leg is a tag and the other is a probe that passes the selection
    // Note: TP is exclusive from 2TT
    TCut cut_TP_1 = (tag2 && probe1 && sel1 && evsel) && !(tag1) && evsel;
    TCut cut_TP_2 = (tag1 && probe2 && sel2 && evsel) && !(tag2) && evsel;
    cut_TP_1.SetName(TString("TP_1"));
    cut_TP_2.SetName(TString("TP_2"));
    // TF means one lef is a tag and the other is a probe that fails the selection
    TCut cut_TF_1 = (tag2 && probe1 && !sel1 && evsel) && !(tag1) && evsel;
    TCut cut_TF_2 = (tag1 && probe2 && !sel2 && evsel) && !(tag2) && evsel;
    cut_TF_1.SetName(TString("TF_1"));
    cut_TF_2.SetName(TString("TF_2"));

    //
    // loop on baby samples
    //

    for(unsigned int i = 0; i < bss.size(); ++i) {

        // construct the selection for this sample
        // and get the plot of that selection

        // fill the right leg variable depending 
        // which leg is the tag

        h1_2TT_1   = Plot(bss[i], var1, cut_2TT, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
        h1_TP_1    = Plot(bss[i], var1, cut_TP_1, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);   
        h1_TF_1    = Plot(bss[i], var1, cut_TF_1, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
        h1_2TT_2   = Plot(bss[i], var2, cut_2TT, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
        h1_TP_2    = Plot(bss[i], var2, cut_TP_2, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
        h1_TF_2    = Plot(bss[i], var2, cut_TF_2, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);

        // now combine them

        //h1_2TT = (TH1F*)h1_2TT_1->Clone();
        //h1_TP = (TH1F*)h1_TP_1->Clone();
        //h1_TF = (TH1F*)h1_TF_1->Clone();
        h1_2TT_1->Add(h1_2TT_2);
        h1_TP_1->Add(h1_TP_2);
        h1_TF_1->Add(h1_TF_2);

        // if the plot doesn't already exist then
        // add it to the appropriate vector of plots
        if (bss[i]->type() == DATA && (find(vh_data_2TT.begin(), vh_data_2TT.end(), h1_2TT_1) == vh_data_2TT.end())) {
            vh_data_2TT.push_back(h1_2TT_1);
            vh_data_TP.push_back(h1_TP_1);
            vh_data_TF.push_back(h1_TF_1);
        } else if (bss[i]->type() == BACKGROUND && find(vh_mc_2TT.begin(), vh_mc_2TT.end(), h1_2TT_1) == vh_mc_2TT.end()) {
            vh_mc_2TT.push_back(h1_2TT_1);
            vh_mc_TP.push_back(h1_TP_1);
            vh_mc_TF.push_back(h1_TF_1);
        }

    }

    // MUST have at least one data and at least one MC
    if (vh_mc_2TT.size() == 0 || vh_data_2TT.size() == 0) {
        std::cout << "[DrawAll] ERROR - MUST have at least one MC and at least one data" << std::endl;
        return 0;
    }

    //
    // sort histograms by their total contributions
    // makes for prettier stack plots
    //

    sort(vh_mc_2TT.begin(), vh_mc_2TT.end(), sortHistsByIntegral);
    sort(vh_mc_TP.begin(), vh_mc_TP.end(), sortHistsByIntegral);
    sort(vh_mc_TF.begin(), vh_mc_TF.end(), sortHistsByIntegral);
    sort(vh_data_2TT.begin(), vh_data_2TT.end(), sortHistsByIntegral);
    sort(vh_data_TP.begin(), vh_data_TP.end(), sortHistsByIntegral);
    sort(vh_data_TF.begin(), vh_data_TF.end(), sortHistsByIntegral);

    //
    // do the stacking
    // NOTE - signals are overlaid NOT stacked
    //

    makeStack(vh_mc_2TT);
    makeStack(vh_mc_TP);
    makeStack(vh_mc_TF);
    makeStack(vh_data_2TT);
    makeStack(vh_data_TP);
    makeStack(vh_data_TF);

    // data denominator and numerator for eff
    TH1F *h1_data_denom = (TH1F*)vh_data_2TT[0]->Clone("data_denom");
    h1_data_denom->Add(vh_data_TP[0]);
    h1_data_denom->Add(vh_data_TF[0]);
    TH1F *h1_data_numer = (TH1F*)vh_data_2TT[0]->Clone("data_numer");
    h1_data_numer->Add(vh_data_TP[0]);

    // mc denominator and numerator for eff
    TH1F *h1_mc_denom = (TH1F*)vh_mc_2TT[0]->Clone("mc_denom");
    h1_mc_denom->Add(vh_mc_TP[0]);
    h1_mc_denom->Add(vh_mc_TF[0]);
    TH1F *h1_mc_numer = (TH1F*)vh_mc_2TT[0]->Clone("mc_numer");
    h1_mc_numer->Add(vh_mc_TF[0]);

    //TGraphAsymmErrors* gr_eff_data = new TGraphAsymmErrors(nbins);
    TH1F *gr_eff_data = (TH1F*)h1_data_numer->Clone();
    gr_eff_data->SetName(TString("gr_") + h1_data_denom->GetName());
    gr_eff_data->SetTitle(TString(savename));
    //gr_eff_data->BayesDivide(h1_data_numer, h1_data_denom);
    ComputeEfficiency(gr_eff_data, h1_data_denom, h1_data_numer);
    gr_eff_data->SetMarkerColor(kRed);

    TGraphAsymmErrors* gr_eff_mc = new TGraphAsymmErrors(nbins);
    gr_eff_mc->SetName(TString("gr_") + h1_mc_denom->GetName());
    gr_eff_mc->SetTitle(TString(savename));
    //gr_eff_mc->Divide(h1_mc_numer, h1_mc_denom);
    gr_eff_mc->SetMarkerColor(kBlue);    

    TLegend* leg = new TLegend(0.7,0.25,0.95,0.45);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);
    //leg->AddEntry(gr_eff_mc, "MC", "lp");
    leg->AddEntry(gr_eff_data, "Data", "lp");

    //
    // do the drawing
    //

    TCanvas *c1 = new TCanvas(savename);
    c1->SetTopMargin(0.08);
    c1->Divide(2, 2);

/*   
    // do the background MC histograms
    for(unsigned int i = 0; i < vh_mc_2TT.size(); ++i)
    {

        // 2TT
        c1->cd(1);
        if (i == 0) {
            vh_mc_2TT[0]->Draw("hist");
            vh_mc_2TT[0]->GetXaxis()->SetNdivisions(504);
        } else vh_mc_2TT[i]->Draw("histsame");
        // TP
        c1->cd(2);
        if (i == 0) {
            vh_mc_TP[0]->Draw("hist");
            vh_mc_TP[0]->GetXaxis()->SetNdivisions(504);
        } else vh_mc_TP[i]->Draw("histsame");
        // TF
        c1->cd(3);
        if (i == 0) {
            vh_mc_TF[0]->Draw("hist");
            vh_mc_TF[0]->GetXaxis()->SetNdivisions(504);
        } else vh_mc_TF[i]->Draw("histsame");

    }
*/

    // do the data histogram
    // 2TT
    c1->cd(1);
    vh_data_2TT[0]->Draw("e1");
    //h1_data_denom->Draw();
    //float ymax = vh_data_2TT[0]->GetMaximum() > vh_mc_2TT[0]->GetMaximum()
    //    ? vh_data_2TT[0]->GetMaximum()+2*sqrt(vh_data_2TT[0]->GetMaximum()) : 
    //    vh_mc_2TT[0]->GetMaximum() + 2*sqrt(vh_mc_2TT[0]->GetMaximum());
    //vh_mc_2TT[0]->SetMaximum(ymax);

    // TP
    c1->cd(2);
    vh_data_TP[0]->Draw("e1");
    //h1_data_numer->Draw();
    //ymax = vh_data_TP[0]->GetMaximum() > vh_mc_TP[0]->GetMaximum()
    //    ? vh_data_TP[0]->GetMaximum()+2*sqrt(vh_data_TP[0]->GetMaximum()) :
    //    vh_mc_TP[0]->GetMaximum() + 2*sqrt(vh_mc_TP[0]->GetMaximum());
    //vh_mc_TP[0]->SetMaximum(ymax);

    // TF
    c1->cd(3);
    vh_data_TF[0]->Draw("e1");
    //ymax = vh_data_TF[0]->GetMaximum() > vh_mc_TF[0]->GetMaximum()
    //    ? vh_data_TF[0]->GetMaximum()+2*sqrt(vh_data_TF[0]->GetMaximum()) :
    //    vh_mc_TF[0]->GetMaximum() + 2*sqrt(vh_mc_TF[0]->GetMaximum());
    //vh_mc_TF[0]->SetMaximum(ymax);

    // now the efficiency
    c1->cd(4);
    //gr_eff_data->Draw("AP");
    gr_eff_data->Draw("HIST E1");
    gr_eff_data->GetYaxis()->SetRangeUser(0.0, 1.1);
    //gr_eff_mc->Draw("P");
    leg->Draw("SAME");

    PrintBins(h1_data_denom);
    PrintBins(h1_data_numer);
    std::cout << "---" << std::endl;
    PrintBins(gr_eff_data);

    // draw the legend and tidy up
    c1->RedrawAxis();
    gDrawAllCount++;
    reset_babydorkidentifier();

    return c1;
}



TCanvas* TriggerMonitor(const char *savename, TCut var, TCut sel, TCut trig, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated, BabySample *bs)
{

    TH1F* h1_pass;
    TH1F* h1_total;

    //
    // compute the efficiency
    //

    TCut cut_pass = sel + trig;
    cut_pass.SetName(TString(sel.GetName())+"_"+TString(trig.GetName()));
    TCut cut_fail = sel + !trig;
    cut_fail.SetName(TString(sel.GetName())+"_not_"+TString(trig.GetName()));

    h1_total    = Plot(bs, var, sel, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
    h1_pass     = Plot(bs, var, cut_pass, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);

    // make a failing histogram in order to get list of
    // failing events
    TH1F *h1_fail = Plot(bs, var, cut_fail, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);
    //
    //

    TGraphAsymmErrors* gr_eff = new TGraphAsymmErrors();
    gr_eff->SetName(TString("gr_") + h1_pass->GetName());
    gr_eff->SetTitle(TString(savename));
    gr_eff->BayesDivide(h1_pass, h1_total);
    //PrintBins(gr_eff);

    //
    // do the drawing
    //

    TCanvas *c1 = new TCanvas(savename);
    c1->SetTopMargin(0.08);
    c1->cd();

    // do the data histogram
    gr_eff->Draw("AP");
    gr_eff->GetXaxis()->SetTitle(var.GetTitle());
    gr_eff->GetYaxis()->SetTitle("Efficiency");
    gr_eff->GetXaxis()->SetNdivisions(504);
    //gr_eff->GetYaxis()->SetRangeUser(0.0, 1.10);
 
    // draw the legend and tidy up
    c1->RedrawAxis();
    gDrawAllCount++;
    reset_babydorkidentifier();
    delete h1_pass;
    delete h1_total;

    return c1;
}

TCanvas* DrawAll(TCut var, const char *savename, TCut sel, float intlumipb, unsigned int nbins, float xlo, float xhi, bool integrated,
        std::vector<BabySample*> bss)
{

    std::vector<TH1F*> vh_background;
    std::vector<TH1F*> vh_signal;
    std::vector<TH1F*> vh_data;
    TH1F* buffer;

    //
    // loop on baby samples
    //

    for(unsigned int i = 0; i < bss.size(); ++i) {

        // construct the selection for this sample
        // and get the plot of that selection
        TCut selection = sel;// + bss[i]->presel();
        buffer = Plot (bss[i], var, selection, intlumipb, nbins, xlo, xhi, integrated, gDrawAllCount);

        // if the plot doesn't already exist then
        // add it to the appropriate vector of plots
        if (bss[i]->type() == DATA && (find(vh_data.begin(), vh_data.end(), buffer) == vh_data.end())) {
            vh_data.push_back(buffer);
        } else if (bss[i]->type() == BACKGROUND && find(vh_background.begin(), vh_background.end(), buffer) == vh_background.end()) {
            vh_background.push_back(buffer);
        } else if (bss[i]->type() == SIGNAL && find(vh_signal.begin(), vh_signal.end(), buffer) == vh_signal.end()) {
            vh_signal.push_back(buffer);
        }


    }

    // MUST have at least one data and at least one MC
    if (vh_background.size() == 0 || vh_data.size() == 0) {
        std::cout << "[DrawAll] ERROR - MUST have at least one MC and at least one data" << std::endl;
        return 0;
    }

    //
    // sort histograms by their total contributions
    // makes for prettier stack plots
    //

    sort(vh_background.begin(), vh_background.end(), sortHistsByIntegral);
    sort(vh_signal.begin(), vh_signal.end(), sortHistsByIntegral);
    sort(vh_data.begin(), vh_data.end(), sortHistsByIntegral);

    //
    // do the stacking
    // NOTE - signals are overlaid NOT stacked
    //

    makeStack(vh_background);
    makeStack(vh_data);

    //
    // set up the legend
    //

    TLegend* leg = new TLegend(0.0, 0.4, 1.0, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);
    leg->SetTextSize(0.08);

    //
    // do the drawing
    //

    TCanvas *c1 = new TCanvas(savename, savename, 800, 600);
    c1->cd();
    TPad *pad1 = new TPad("p_main", "p_main", 0.0, 0.0, 0.75, 1.0);
    pad1->SetBottomMargin(0.13);
    pad1->SetRightMargin(0.07);
    pad1->Draw();
    c1->cd();
    TPad *pad2 = new TPad("p_leg", "p_leg", 0.75, 0.0, 1.0, 1.0);
    pad2->SetTopMargin(0.01);
    pad2->SetRightMargin(0.01);
    pad2->SetBottomMargin(0.13);
    pad2->Draw();

    pad1->cd();
    // do the background MC histograms
    for(unsigned int i = 0; i < vh_background.size(); ++i) 
    { 
        // the zeroth sets the axis
        if (i == 0) {
            vh_background[0]->Draw("hist");
            vh_background[0]->GetXaxis()->SetNdivisions(504);
        }
        // the others are drawn the same
        else vh_background[i]->Draw("histsame");
        // add each histogram to the legend
        addToLegend(leg, vh_background.at(i), "f");
    }

    // do the signal MC histograms
    for(unsigned int i = 0; i < vh_signal.size(); ++i)
    {
        vh_signal[i]->Draw("histsame");
        addToLegend(leg, vh_signal.at(i), "l");
    }

    // do the data histogram
    vh_data[0]->Draw("samee1");
    addToLegend(leg, vh_data[0], "lp");

    // set the optimum y axis scale
    float ymax = vh_data[0]->GetMaximum() > vh_background[0]->GetMaximum() 
        ? vh_data[0]->GetMaximum()+2*sqrt(vh_data[0]->GetMaximum()) : 
        vh_background[0]->GetMaximum() + 2*sqrt(vh_background[0]->GetMaximum());
    if (vh_signal.size() > 0) {
        ymax = ymax > vh_signal[0]->GetMaximum() ? ymax : 
            vh_signal[0]->GetMaximum() + 2*sqrt(vh_signal[0]->GetMaximum());
    }
    vh_background[0]->SetMaximum(ymax);

    // draw the legend and tidy up
    pad2->cd();
    leg->Draw();

    c1->cd();
    c1->RedrawAxis();
    gDrawAllCount++;
    reset_babydorkidentifier();

    return c1;
}

TH1F* Plot(BabySample *bs, TCut var, TCut selection, float intlumipb, 
        unsigned int nbins, float xlo, float xhi, bool integrated, unsigned int isfx)
{

    //
    // Construct the name and title
    // and command to draw
    // 

    if (!strcmp("CUT", var.GetName()))
        var.SetName(var.GetTitle());
   
    char *name = 0;
    if (integrated) name = Form("%s_%s_%s_%i_int", bs->pfx(), selection.GetName(), var.GetName(), gDrawAllCount);
    else name = Form("%s_%s_%s_%i", bs->pfx(), selection.GetName(), var.GetName(), gDrawAllCount);
    char *title = Form("%s_%s_%s, ~%.2f/pb", bs->pfx(), selection.GetName(), var.GetName(), 1e3*intlumipb);
    char *drawcommand = Form("%s>>+%s", var.GetTitle(), name);

    //
    // set the normalisation scale
    //
    TCut scale = Form("scale1fb*(%f/1000.0)*%f", intlumipb, bs->kfactor());
    if (bs->type() == DATA || bs->type() == TPMC) scale = "1.0";

    //
    // If the histogram does not already exist
    // then create it
    //

    TH1F *h = 0;
    if (!(h = (TH1F*)gROOT->FindObjectAny(name))) {
        h = new TH1F(name, title, nbins, xlo, xhi);
        h->Sumw2();
    }

    //
    // Draw the histogram
    //

    // fill an event list for the plot selection
    // and then get this event list
    bs->chain()->Draw(">>evtlist", selection, "goff");
    TEventList* elist = (TEventList*)gDirectory->Get("evtlist");

    // set the event list for events passing
    // the plot selection and draw the plot
    bs->chain()->SetEventList(elist);
    bs->chain()->Draw(drawcommand, scale, "goff");

    // if the events are data then
    // use the event list to print a scan
    if (bs->type() == DATA) {
        TTreePlayer *tp = (TTreePlayer*)bs->chain()->GetPlayer();
        tp->SetScanRedirect(kTRUE);
        tp->SetScanFileName(Form("../output/%s_%s_%s_%s.out", selection.GetName(), var.GetName(), bs->pfx(), bs->pfx2()));
        // trileptons
        if (bs->chain()->GetBranch("pt3"))
            bs->chain()->Scan("dataset:run:ls:evt:pt1:pt2:pt3", "", "colsize=20");
        // dileptons
        else if (bs->chain()->GetBranch("pt2")) // dilep
            bs->chain()->Scan("dataset:run:ls:evt:pt1:pt2", "", "colsize=20");
        // single leptons
        else
            bs->chain()->Scan("dataset:run:ls:evt:pt1", "", "colsize=20");
    }

    // reset the event list for this baby
    // to the one internally stored
    // usually corresponding to a common preselection
    bs->resetEventList();
    elist->Reset();

    //
    // Set style options
    //

    if (bs->type() == DATA) {
        h->SetMarkerColor(bs->color());
        h->SetMarkerStyle(bs->style());
    } else if (bs->type() == BACKGROUND || bs->type() == TPMC) {
        h->SetLineColor(bs->color());
        h->SetFillColor(bs->color());
        h->SetFillStyle(bs->style());
    } else if (bs->type() == SIGNAL) {
        h->SetLineColor(bs->color());
        h->SetLineStyle(bs->style());
        h->SetLineWidth(2);
    }

    h->SetTitle(Form("%s, ~%.2f/pb", selection.GetName(), intlumipb));
    if (integrated) {
        h = slideIntegrated(h);
        h->GetXaxis()->SetTitle(Form("integrated %s", var.GetName()));
    }
    else
        h->GetXaxis()->SetTitle(var.GetName());

    //
    // Move overflow to the last bin
    //

    float overflowerr = h->GetBinError(nbins+1);
    float overflow = h->GetBinContent(nbins+1);
    float endbinerr = h->GetBinError(nbins);
    float endbin = h->GetBinContent(nbins);
    h->SetBinContent(nbins,overflow+endbin);
    h->SetBinError(nbins, sqrt(overflow+endbin));
    h->SetBinContent(nbins+1,0.);
    h->SetEntries(h->GetEntries() - overflow);

    return h;

}

bool sortHistsByIntegral(TH1* h1, TH1* h2)
{
    return h1->Integral() > h2->Integral();
}

void addToLegend(TLegend *leg, TH1F *hist, TString opt)
{

    // For extracting sample prefix
    TPRegexp preg("^([^_]+)_.*$");
    TString  s_pfx("");

    s_pfx   = ((TObjString*)(preg.MatchS(TString(hist->GetName()))->At(1)))->GetString();
    leg->AddEntry(hist, s_pfx.Data(), opt);

}

void makeStack(std::vector<TH1F*> &v_hists)
{

    for(unsigned int i = 0; i < v_hists.size(); ++i) {
        for(unsigned int j = i+1; j < v_hists.size(); ++j) {
            (v_hists[i])->Add(v_hists[j], 1);
        }
    }

}


TH1F* slideIntegrated(TH1F* h1)
{
    int NbinsX = h1->GetNbinsX();
    for(int i = 1; i <= NbinsX; ++i) {
        // don't forget the overflow!
        float integral = h1->Integral(i,NbinsX+1);
        h1->SetBinContent(i,integral);
        h1->SetBinError(i,TMath::Sqrt(integral));
    }

    return h1;
}

void PreselectBabies(std::vector<BabySample*> bss, TCut cut)
{

    // set up the good run list cut
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    TCut c_goodrunplus(Form("((((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun(run,ls))||(run>%i||(run==%i&&ls>%i)))||!isdata",
                brun, erun, brun, bls, erun, els, erun, erun, els));

    std::cout << "Preselecting Babies..." << std::endl;
    for (unsigned int i = 0; i < bss.size(); ++i) {
        bss[i]->setEventList(cut+c_goodrunplus);
    }

}

void PrintBins(const TH1F *h1)
{

    unsigned int nbins = h1->GetNbinsX();
    for (unsigned int i = 1; i <= nbins; ++i)
    {
        float low = h1->GetBinLowEdge(i);
        float high = h1->GetBinWidth(i);
        std::cout << "Bin: " << low << " - " << high << ": " << h1->GetBinContent(i) << " \\pm " << h1->GetBinError(i) << std::endl; 
    }

}

void PrintBins(const TGraphAsymmErrors *h1)
{

    unsigned int nbins = h1->GetN();
    for (unsigned int i = 0; i < nbins; ++i)
    {
        float low = h1->GetErrorXlow(i);
        float high = h1->GetErrorXhigh(i);
        Double_t x, y;
        h1->GetPoint(i, x, y);
        std::cout << "Bin: " << x - low << " - " << x + high << ": " << y << " - " << h1->GetErrorYlow(i) << " + " << h1->GetErrorYhigh(i) << std::endl;
    }

}

void ZeroBinError(TH1F *h1)
{

    unsigned int nbins = h1->GetNbinsX();
    for (unsigned int i = 1; i <= nbins; ++i)
    {
        h1->SetBinContent(i, int(h1->GetBinContent(i)));
        h1->SetBinError(i, sqrt(int(h1->GetBinContent(i))));
    }

}

void ComputeEfficiency(TH1F *gr_eff_data, const TH1F *h1_data_denom, const TH1F *h1_data_numer)
{

    unsigned int nbins = h1_data_denom->GetNbinsX();
    for (unsigned int i = 1; i <= nbins; ++i)
    {

        float eff = 0.0;
        float err = 0.0;
        if (h1_data_denom->GetBinContent(i) != 0) {
            eff = h1_data_numer->GetBinContent(i) / h1_data_denom->GetBinContent(i);
            err = sqrt(eff*(1-eff)/h1_data_denom->GetBinContent(i));
        }
        gr_eff_data->SetBinContent(i, eff);
        gr_eff_data->SetBinError(i, err);
        std::cout << "setting (" << i << "): " << eff << " \\pm " << err << std::endl;
    }

}




