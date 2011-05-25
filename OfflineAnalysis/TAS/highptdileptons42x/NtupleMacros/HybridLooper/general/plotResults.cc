
#include "plotResults.h"
#include "../../Tools/HistogramUtilities.h"
#include "../../Tools/AllDataSources.h"
#include "../../Tools/Utilities.h"
#include "../../Tools/histtools.cc"

#include "TROOT.h"

#include "TGraphAsymmErrors.h"
#include "TArrow.h"

	// for 2_2_X
	const static sources_t theSignal = 
		(1<<H_WJETS);

	const static sources_t theBackground = 
		(1<<H_QCD30);
	const static sources_t theSources = 
		(1<<H_WJETS) |
		(1<<H_PHOTONJETS) |
		(1<<H_QCD30);

TArrow *getArrow(THStack *st, TString det, float cutValEB, float cutValEE)
{

        float cutVal; 
        if (det == "eb") cutVal = cutValEB;
        else cutVal = cutValEE;
        TArrow *arr_cut = new TArrow(cutVal, st->GetMaximum()/2.0, cutVal, 0, 0.05, "|>");
        arr_cut->SetLineWidth(2.0);
	return arr_cut;
}

void plot2DSB(HistogramUtilities &h1, TString name, TString xTitle, TString yTitle, TString saveName, TString det)
{

        TH2F *h2_signal = h1.get2dHistogram(theSignal, name, "", det, 1);
	h2_signal->SetTitle(";" + xTitle + ";" + yTitle);
        TH2F *h2_background = h1.get2dHistogram(theBackground, name, "", det, 1);
        h2_background->SetTitle(";" + xTitle + ";" + yTitle);
	TCanvas *c1 = new TCanvas();
	c1->SetCanvasSize(600, 300);
	c1->Divide(2, 1);
	c1->cd(1);
	h2_signal->Draw("COLZ");
	c1->cd(2);
	h2_background->Draw("COLZ");
        Utilities::saveCanvas(c1, "results/" + saveName + "_2D_linz_" + name + "_" + det);

        c1->cd(1);
	c1->cd(1)->SetLogz();
        h2_signal->Draw("COLZ");
        c1->cd(2);
	c1->cd(2)->SetLogz();
        h2_background->Draw("COLZ");
        Utilities::saveCanvas(c1, "results/" + saveName + "_2D_logz_" + name + "_" + det);

	delete c1;
	delete h2_signal;
	delete h2_background;

}

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin, bool legendOnRight, float cutValEB, float cutValEE)
{

	std::cout << "[plotEff] " << name << std::endl;
	TH1F *h1_signal = h1.getHistogram(theSignal, name, "", det, rebin, "s_");
        TH1F *h1_background = h1.getHistogram(theBackground, name, "", det, rebin, "b_");
	TH1F *h1_both = h1.getHistogram(theBackground | theSignal, name, "", det, rebin, "");

	//
	// S/B plots
	//

	TH1F *h1_signal_cumulated = (TH1F*)(cumulate(*h1_signal, ascending).Clone());
        TH1F *h1_background_cumulated = (TH1F*)(cumulate(*h1_background, ascending).Clone());
	TH1F *h1_sob = (TH1F*)h1_signal_cumulated->Clone();
	float sTotal = h1_signal->Integral(0, h1_signal_cumulated->GetNbinsX() + 1);
        float bTotal = h1_background->Integral(0, h1_background_cumulated->GetNbinsX() + 1);

	//int bin_bg50 = 0;
	int bin_eff99 = 0;
	int bin_eff98 = 0;
	int bin_eff95 = 0;
	int bin_eff80 = 0;
	for (Int_t bin = 0; bin < h1_signal_cumulated->GetNbinsX() + 1; ++bin)
	{
		float s = h1_signal_cumulated->GetBinContent(bin);
		float b = h1_background_cumulated->GetBinContent(bin);
		float sob = 0;
		if (b > 0.0) { 
			sob = s/b;
		}
		h1_sob->SetBinContent(bin, sob);
		//std::cout << s/sTotal << std::endl;
		if (ascending) {
			if (s/sTotal >= 0.99 && bin_eff99 == 0) bin_eff99 = bin;
                	if (s/sTotal >= 0.98 && bin_eff98 == 0) bin_eff98 = bin;
                	if (s/sTotal >= 0.95 && bin_eff95 == 0) bin_eff95 = bin;
                	if (s/sTotal >= 0.80 && bin_eff80 == 0) bin_eff80 = bin;
		}
		else {
                        if (s/sTotal <= 0.99 && bin_eff99 == 0) bin_eff99 = bin;
                        if (s/sTotal <= 0.98 && bin_eff98 == 0) bin_eff98 = bin;
                        if (s/sTotal <= 0.95 && bin_eff95 == 0) bin_eff95 = bin;
                        if (s/sTotal <= 0.80 && bin_eff80 == 0) bin_eff80 = bin;
		}

		//if (b/bTotal >= 0.50 && bin_bg50 == 0) bin_bg50 = bin;
	}

        TArrow *arr_eff99 = new TArrow(h1_signal->GetBinCenter(bin_eff99), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff99), 0, 0.05, "|>");
	arr_eff99->SetLineColor(kGreen);
	arr_eff99->SetFillColor(kGreen);
	arr_eff99->SetLineWidth(2);
        TArrow *arr_99 = new TArrow(0.99, 1.0, 0.99, 0, 0.05, "|>");
        arr_99->SetLineColor(kGreen);
        arr_99->SetFillColor(kGreen);
        arr_99->SetLineWidth(2);

        TArrow *arr_eff98 = new TArrow(h1_signal->GetBinCenter(bin_eff98), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff98), 0, 0.05, "|>");
        arr_eff98->SetLineColor(kBlue);
        arr_eff98->SetFillColor(kBlue);
        arr_eff98->SetLineWidth(2);
        TArrow *arr_98 = new TArrow(0.98, 1.0, 0.98, 0, 0.05, "|>");
        arr_98->SetLineColor(kBlue);
        arr_98->SetFillColor(kBlue);
        arr_98->SetLineWidth(2);

        TArrow *arr_eff95 = new TArrow(h1_signal->GetBinCenter(bin_eff95), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff95), 0, 0.05, "|>");
        arr_eff95->SetLineColor(kRed);
        arr_eff95->SetFillColor(kRed);
        arr_eff95->SetLineWidth(2);
        TArrow *arr_95 = new TArrow(0.95, 1.0, 0.95, 0, 0.05, "|>");
        arr_95->SetLineColor(kRed);
        arr_95->SetFillColor(kRed);
        arr_95->SetLineWidth(2);

        TArrow *arr_eff80 = new TArrow(h1_signal->GetBinCenter(bin_eff80), h1_both->GetMaximum()/2.0, h1_signal->GetBinCenter(bin_eff80), 0, 0.05, "|>");
        arr_eff80->SetLineColor(kMagenta);
        arr_eff80->SetFillColor(kMagenta);
        arr_eff80->SetLineWidth(2);
        TArrow *arr_80 = new TArrow(0.80, 1.0, 0.80, 0, 0.05, "|>");
        arr_80->SetLineColor(kMagenta);
        arr_80->SetFillColor(kMagenta);
        arr_80->SetLineWidth(2);

        //TArrow *arr_bg50 = new TArrow(h1_background->GetBinCenter(bin_bg50), h1_both->GetMaximum()/2.0, h1_background->GetBinCenter(bin_bg50), 0, 0.05, "|>");
        //arr_bg50->SetLineColor(kBlack);
        //arr_bg50->SetFillColor(kBlack);
        //arr_bg50->SetLineWidth(2);

        TCanvas *c = new TCanvas();
	c->cd();	
	h1_sob->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_sob_" + name + "_" + det);	

        //
        // Acceptance rejection curves
        //

        // latter bool is "ascending"
        bool normalise = true;
        TGraph *gr = (TGraph*)(eff_rej(*h1_signal, *h1_background, normalise, ascending).Clone());
        gr->SetTitle(name + ";Signal;Background");
        gr->SetMarkerStyle(23);

        c->cd();
        gr->Draw("AP");
        gr->GetXaxis()->SetRangeUser(0.95, 1.1);
        gr->GetYaxis()->SetRangeUser(0.00, 1.1); 
        arr_99->Draw();
        arr_98->Draw();
        arr_95->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_eff_" + name + "_" + det);

        gr->GetXaxis()->SetRangeUser(0.00, 1.1);
        gr->GetYaxis()->SetRangeUser(0.00, 1.1);
        arr_80->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_effzoomout_" + name + "_" + det);


	//
	// S and B overlays
	//

        h1_signal->SetLineWidth(2);
	h1_signal->SetLineColor(kBlue);
	h1_signal->SetLineStyle(kDashed);
        h1_background->SetLineWidth(2);
	h1_background->SetLineColor(kGreen);
	TLegend *lg = 0;
	if (! legendOnRight) lg = new TLegend(0.2, 0.75, 0.5, 0.9);
	else lg = new TLegend(0.55, 0.7, 0.9, 0.9);

	float cutVal;
	if (det == "eb") cutVal = cutValEB;
	else cutVal = cutValEE;
	TArrow *arr_cut = new TArrow(cutVal, h1_both->GetMaximum()/2.0, cutVal, 0, 0.05, "|>");
	arr_cut->SetLineWidth(2.0);

  	lg->SetFillColor(kWhite);
  	lg->SetLineColor(kWhite);
  	lg->SetShadowColor(kWhite);
	TString upperDet = det;
	upperDet.ToUpper();
	lg->AddEntry(h1_signal, "Signal (" + upperDet + ")", "l");
	lg->AddEntry(h1_background, "Background (" + upperDet + ")", "l");
	TCanvas *c_sb = new TCanvas();
	c_sb->cd();
	if (h1_signal->GetMaximum() > h1_background->GetMaximum()) {
		h1_signal->Draw();
		h1_signal->GetYaxis()->SetRangeUser(0, h1_signal->GetMaximum() 
							+ h1_signal->GetMaximum()*0.1);
		h1_background->Draw("SAME");
	}
	else {
        	h1_background->Draw();
                h1_background->GetYaxis()->SetRangeUser(0, h1_background->GetMaximum() 
								+ h1_background->GetMaximum()*0.1);
        	h1_signal->Draw("SAME");
	}
	if (cutVal > 0.0) arr_cut->Draw();
	if (bin_eff99 != 0) arr_eff99->Draw();
        if (bin_eff98 != 0) arr_eff98->Draw();
        if (bin_eff95 != 0) arr_eff95->Draw();
	if (bin_eff80 != 0) arr_eff80->Draw();
	//if (bin_bg50 != 0) arr_bg50->Draw();

	lg->Draw();
        Utilities::saveCanvas(c_sb, "results/" + saveName + "_sb_lin" + name + "_" + det);

	h1_signal->Draw();
        Utilities::saveCanvas(c_sb, "results/" + saveName + "_s_lin" + name + "_" + det);


	TCanvas *c_sb_log = new TCanvas();
	c_sb_log->cd();
	c_sb_log->SetLogy();
        if (h1_signal->GetMaximum() > h1_background->GetMaximum()) {
                h1_signal->Draw();
                h1_signal->GetYaxis()->SetRangeUser(0.1, h1_signal->GetMaximum() 
                                                        + h1_signal->GetMaximum());
                h1_background->Draw("SAME");
        }
        else {
                h1_background->Draw();
                h1_background->GetYaxis()->SetRangeUser(0.1, h1_background->GetMaximum() 
                                                                + h1_background->GetMaximum());
                h1_signal->Draw("SAME");
        }
        if (cutVal > 0.0) arr_cut->Draw();
        if (bin_eff99 != 0) arr_eff99->Draw();
        if (bin_eff98 != 0) arr_eff98->Draw();
        if (bin_eff95 != 0) arr_eff95->Draw();
	if (bin_eff80 != 0) arr_eff80->Draw();
        //if (bin_bg50 != 0) arr_bg50->Draw();

        lg->Draw();
	Utilities::saveCanvas(c_sb_log, "results/" + saveName + "_sb_log" + name + "_" + det);

        h1_signal->Draw();
        Utilities::saveCanvas(c_sb_log, "results/" + saveName + "_s_log" + name + "_" + det);

	delete gr;
	delete h1_signal;
	delete h1_background;
	delete c;
	delete c_sb;
	delete c_sb_log;
	delete lg;
	delete arr_cut;
	delete h1_both;
	delete arr_eff99;
	delete arr_eff98;
	delete arr_eff95;
	//delete arr_bg50;
	std::cout << "[plotEff] Done" << std::endl;
}



void plotEffVar(HistogramUtilities &h1, TString name, TString det, TString saveName, Int_t rebin = 1)
{

        TCanvas *c = new TCanvas();
        TH1F *h1_eff = 0;
        TH1F *h1_total = 0;
        TH1F *h1_pass = 0;
	TGraphAsymmErrors *gr_eff_s = new TGraphAsymmErrors();
        TGraphAsymmErrors *gr_eff_bg = new TGraphAsymmErrors();

        h1_total = h1.getHistogram(theSignal, name, "", det + "_denom", rebin, "");
        h1_pass = h1.getHistogram(theSignal, name, "", det + "_numer", rebin, "");
        //h1_total->Rebin(rebin);
        //h1_pass->Rebin(rebin);
	gr_eff_s->BayesDivide(h1_pass, h1_total);
        gr_eff_s->SetMarkerColor(kBlue);
        gr_eff_s->SetLineColor(kBlue);
        //h1_eff = (TH1F*)h1_pass->Clone();
        //h1_eff->Reset();
        //h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        //h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_s_") + h1_pass->GetName());
        c->cd();
        //h1_eff->SetLineWidth(2);
        //h1_eff->SetMarkerStyle(20);
	//h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
	//h1_eff->GetXaxis()->SetRangeUser(0, 100.0);
        //h1_eff->Draw("E1");
	gr_eff_s->Draw("AP");
	//gr_eff_s->GetXaxis()->SetRangeUser(0, 100.0);
        gr_eff_s->GetYaxis()->SetRangeUser(0, 1.1);
	gr_eff_s->GetXaxis()->SetTitle(h1_total->GetXaxis()->GetTitle());
        gr_eff_s->GetYaxis()->SetTitle("Efficiency");
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_s_" + name + "_" + det);

        h1_total = h1.getHistogram(theBackground, name, "", det + "_denom", rebin, "");
        h1_pass = h1.getHistogram(theBackground, name, "", det + "_numer", rebin, "");
	//h1_total->Rebin(rebin);
	//h1_pass->Rebin(rebin);
        gr_eff_bg->BayesDivide(h1_pass, h1_total);
	gr_eff_bg->SetMarkerColor(kGreen);
	gr_eff_bg->SetLineColor(kGreen);
        //h1_eff->Reset();
        //h1_eff->Divide(h1_pass, h1_total, 1.0, 1.0, "B");
        //h1_eff->SetName(h1_pass->GetName());
        c->SetName(TString("c_bg_") + h1_pass->GetName());
        c->cd();
	//h1_eff->SetLineWidth(2);
	//h1_eff->SetMarkerStyle(20);
        //h1_eff->GetYaxis()->SetRangeUser(0, 1.1);
        //h1_eff->GetXaxis()->SetRangeUser(0, 100.0);
        //h1_eff->Draw("E1");
	gr_eff_bg->Draw("AP");
        //gr_eff_bg->GetXaxis()->SetRangeUser(0, 100.0);
        gr_eff_bg->GetYaxis()->SetRangeUser(0, 1.1);
        gr_eff_bg->GetXaxis()->SetTitle(h1_total->GetXaxis()->GetTitle());
	gr_eff_bg->GetYaxis()->SetTitle("Efficiency");
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_bg_" + name + "_" + det);

        TLegend *lg = new TLegend(0.5, 0.2, 0.9, 0.4);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
	lg->SetFillStyle(0);
        lg->SetShadowColor(kWhite);
        TString upperDet = det;
        upperDet.ToUpper();
        lg->AddEntry(gr_eff_s, "Signal (" + upperDet + ")", "lp");
        lg->AddEntry(gr_eff_bg, "Background (" + upperDet + ")", "lp");

        c->SetName(TString("c_sb_") + h1_pass->GetName());
        c->cd();
	c->Clear();
        gr_eff_s->Draw("AP");
        gr_eff_bg->Draw("P");
	lg->Draw();
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_sb_" + name + "_" + det);

	//
	// like "data" (use all sources)
	//

        h1_total = h1.getHistogram(theSources, name, "", det + "_denom", rebin, "");
        h1_pass = h1.getHistogram(theSources, name, "", det + "_numer", rebin, "");
        gr_eff_s->BayesDivide(h1_pass, h1_total);
        gr_eff_s->SetMarkerColor(kBlue);
        gr_eff_s->SetLineColor(kBlue);
        c->SetName(TString("c_data_") + h1_pass->GetName());
        c->cd();
        gr_eff_s->Draw("AP");
        gr_eff_s->GetYaxis()->SetRangeUser(0, 1.1);
        gr_eff_s->GetXaxis()->SetTitle(h1_total->GetXaxis()->GetTitle());
        gr_eff_s->GetYaxis()->SetTitle("Efficiency");
        Utilities::saveCanvas(c, "results/" + saveName + "_effVar_data_" + name + "_" + det);

        delete c;
	delete lg;
        //delete h1_eff;
        delete h1_pass;
        delete h1_total;
	delete gr_eff_s;
	delete gr_eff_bg;
}

void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin, float cutValEB, float cutValEE)
{

        THStack *st = h1.getStack(theSources, name, "", det, rebin);
        TLegend *lg_all = h1.getLegend(theSources, name, "", det);
	lg_all->SetX1(0.4);

        TCanvas *c1 = new TCanvas();
        c1->cd();
        st->Draw();
	st->GetXaxis()->SetTitle(titleX); 
        lg_all->Draw();
	if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        Utilities::saveCanvas(c1, "results/" + saveName  + "_lin_" + name + "_" + det);

	c1->SetLogy();
	st->SetMinimum(1.0);
        st->Draw();	
	lg_all->Draw();
        if (det == "ee" && cutValEE != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        if (det == "eb" && cutValEB != -1) getArrow(st, det, cutValEB, cutValEE)->Draw();
        Utilities::saveCanvas(c1, "results/" + saveName  + "_log_" + name + "_" + det);

	delete c1;
	delete st;
	delete lg_all;

}

void plotResultsW(TString det, TString fileStamp)
{

        gROOT->ProcessLine(".L tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");

        // luminorm for 1pb-1
        // luminosity is already normalised to 1pb-1 in the looper
        std::vector<DataSource> sources;
        sources.push_back( fH_WJETS() );
        sources.push_back( fH_QCD30() );
        sources.push_back( fH_PHOTONJETS() );
        HistogramUtilities h1(fileStamp + ".root", sources, 1.0);

        // W studies related
        //

        plotStack(h1, "pt", "p_{T} (GeV)", fileStamp, det, 2);

	// N-1
	plotStack(h1, "nm1_met", "tcMet (GeV)", fileStamp, det, 2);

}

