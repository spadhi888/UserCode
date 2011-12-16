#include "Utils/SMS_utils.C"

#include "Utils/macro_utils.C"
#include "TObjArray.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <sstream>
#include <TGraph2D.h>


using namespace std;

TStyle * gStyle;


void drawGraph() {

   TCanvas *tempGraph = new TCanvas();
   tempGraph->cd();
   
   gStyle->SetCanvasDefH(400);
   gStyle->SetCanvasDefW(400);
   //myFileGraph.root
   gStyle->SetPadLeftMargin(0.4);
   gStyle->SetPadRightMargin(1.);

   TH2D *empty = new TH2D("empty","",60,0.,1500.,60,0.,1500.);
   empty->SetMaximum(1);
   empty->SetTitle("");
   empty->SetXTitle("m_{#tilde{g}} (GeV)");
   empty->SetYTitle("m_{#chi^{0}} (GeV)");   

   empty->Draw("");

   TLine *line;
   line = new TLine(50., 50., 1200., 1200.);
   line->SetLineStyle(3);   

   TFile *file = TFile::Open("myFileGraph.root");
   
   TH2F * h_alphaT = (TH2F*) file->Get("graph_h_T1_alphaT_limit_1.000000");
   TH2F * h_RA2b = (TH2F*) file->Get("graph_h_RA2b_best_T1bbbb_cloned_1.000000");
   TH2F * h_MT2b = (TH2F*) file->Get("graph_h_MT2blow_T1bbbb_cloned_1.000000");
   
   h_alphaT->SetLineColor(1);
   h_RA2b->SetLineColor(2);
   h_MT2b->SetLineColor(4);

   h_alphaT->SetLineWidth(5);
   h_RA2b->SetLineWidth(5);
   h_MT2b->SetLineWidth(5);

   if(h_alphaT) h_alphaT->Draw("L same");
   if(h_RA2b) h_RA2b->Draw("L same");
   if(h_MT2b) h_MT2b->Draw("L same");
   line ->Draw("same");    

   char * legT1bbbb="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2b + LSP; m(#tilde{q})>>m(#tilde{g})";

   TLatex latexLabel;
   latexLabel.SetTextSize(0.035);
   latexLabel.SetNDC();
   
   latexLabel.DrawLatex(0.2, 0.75, legT1bbbb);   

  
   TLegend* this_leg = new TLegend(0.7,0.15,0.9,0.45);
   this_leg->SetFillColor(0);
   this_leg->SetBorderSize(0);
   this_leg->SetTextSize(0.035);
   
   if(h_alphaT) this_leg->AddEntry(h_alphaT,"alphaT" , "l");
   if(h_RA2b) this_leg->AddEntry(h_RA2b,"RA2b" , "l");
   if(h_MT2b) this_leg->AddEntry(h_MT2b,"MT2b" , "l");
  
   this_leg->Draw();

   tempGraph->SaveAs("graph.gif");
   tempGraph->Print("graph.eps");
   
   return;


}

/////

TH1 * get1dHisto(char * fileName, char * histoName, double LSPmass){

  cout << "fileName " << fileName << endl;
  cout << "histoName " << histoName << endl;

  TFile *RA4 = TFile::Open(fileName);
  
  TH2F * h_RA4_loose_2d = (TH2F*) RA4->Get(histoName);

  TString name = TString(h_RA4_loose_2d->GetName());

  //  TH1 * h_RA4_loose = (TH2F*) RA4->Get("h_RA4_T2tt_Loose");
  //  TH1 * h_RA4_loose= new TH1D(histoName,histoName,200, 0., 2000.);
  TH1 * h_RA4_loose= new TH1D(histoName,histoName,80, 0., 2000.);

  double bin_Y=h_RA4_loose_2d->GetYaxis()->FindBin(LSPmass);
  double bin_X=h_RA4_loose_2d->GetXaxis()->FindBin(225);

  //  if(name.Contains("pl"))  cout << "totBinX " << h_RA4_loose_2d->GetNbinsX() << " totBinY " << h_RA4_loose_2d->GetNbinsY() << " bin_X " << bin_X  << " bin_Y " << bin_Y << " content "<< h_RA4_loose_2d->GetBinContent(bin_X,bin_Y) << endl;

  if(name.Contains("pl")) bin_Y=5;
    //cout << "found a non emty bin " << h_RA4_loose_2d->GetBinContent(5,10) << endl;

  /*
  //   for(double i=0; i<(h_RA4_loose_2d->GetNbinsX()+1); i++) {
  for(double i=0; i<(9+1); i++) {
    //      for(double j=0; j<(h_RA4_loose_2d->GetNbinsX()+1); j++) {
      for(double j=0; j<(30+1); j++) {
	
	bool empty=(h_RA4_loose_2d->GetBinContent(i,j)==0);

	if(name.Contains("pl") && (!empty)) cout << "found non empty bin X " << i << " Y "<< j << " "<< h_RA4_loose_2d->GetBinContent(i,j) <<  endl;
	
      }
    }
  */

  //  Int_t  nBinX= limit->GetXaxis()->GetNbins();
  //  double xLow=limit->GetXaxis()->GetBinLowEdge(1); 
  //  double xHigh=limit->GetXaxis()->GetBinLowEdge(nBinX+1);

  for(double j=0; j<(h_RA4_loose_2d->GetNbinsX()+1); j++) {

    /*
    if(name.Contains("pl")) {
      cout << "bin_X " << bin_X << " bin content " << h_RA4_loose_2d->GetBinContent(bin_X,bin_Y) << endl;      
      cout << "bin_Y " << bin_Y << " bin content " << h_RA4_loose_2d->GetBinContent(j,bin_Y) << endl;
    }
    */

    //    if(!(h_RA4_loose_2d->GetBinContent(j,bin_Y)==0)) h_RA4_loose->Fill(h_RA4_loose_2d->GetXaxis()->GetBinCenter(j),h_RA4_loose_2d->GetBinContent(j,bin_Y)); 
    //    if(!(h_RA4_loose_2d->GetBinContent(j,bin_Y)==0)) h_RA4_loose->Fill(h_RA4_loose_2d->GetXaxis()->GetBinCenter(j),h_RA4_loose_2d->GetBinContent(j,bin_Y)); 
    if(!(h_RA4_loose_2d->GetBinContent(j,bin_Y)==0) && (name.Contains("pl"))) h_RA4_loose->Fill(h_RA4_loose_2d->GetXaxis()->GetBinCenter(j)+25,h_RA4_loose_2d->GetBinContent(j,bin_Y)); 
    if(!(h_RA4_loose_2d->GetBinContent(j,bin_Y)==0) && (!name.Contains("pl"))) h_RA4_loose->Fill(h_RA4_loose_2d->GetXaxis()->GetBinCenter(j),h_RA4_loose_2d->GetBinContent(j,bin_Y)); 

  }

  /*
  TCanvas *cPlot1 = new TCanvas("cPlot1","cPlot1");
  
  h_RA4_loose_2d->Draw("same p");
  
  if(name.Contains("pl")) cPlot1->SaveAs("test.C");
  */

  return h_RA4_loose;

}


void gluinoStopSummary() {

  TH1 *hRef=getHisto("reference_xSec.root", "gluino" , "0", 1,1);
  
  TH1 * h_RA2b_best=get1dHisto("RA2bT1tttt/T1tttt_results.root","best_xsUL_all", 50.);
  TH1 * h_SS_best=get1dHisto("RESULT/SS_T1tttt_bestLimit_matched.root","crossSec_cloned", 50.);
  
  //  cout << " file2 " << endl;
  
  TCanvas *cPlot1 = new TCanvas("cPlot1","cPlot1");
  cPlot1->SetLogy(1);
  cPlot1->SetGridy(1);
  cPlot1->SetGridx(1);

  //  hRef->SetRangeUser(300.,800);
  hRef->GetYaxis()->SetTitle("#sigma x BR (#tilde{g}#tilde{g} #rightarrow 4t+2LSP) (pb)");
  hRef->GetXaxis()->SetTitle("m_{#tilde{g}} GeV");
  //  hRef->GetXaxis()->SetRangeUser(300.,600.);
  hRef->GetXaxis()->SetRangeUser(400.,800.);

  hRef->SetTitle("");
  hRef->SetMaximum(100);
  hRef->SetMinimum(0.001);
  hRef->Draw("hist");

  h_RA2b_best->SetMarkerColor(46);
  h_RA2b_best->SetMarkerStyle(20);
  h_RA2b_best->Draw("p same");

  h_SS_best->SetMarkerColor(38);
  h_SS_best->SetMarkerStyle(20);
  h_SS_best->Draw("p same");

  TLegend* this_leg = new TLegend(0.2,0.15,0.9,0.4);
  this_leg->SetFillColor(0);
  this_leg->SetBorderSize(0);
  this_leg->SetTextSize(0.035);

  if(hRef) this_leg->AddEntry(hRef,"(NLO Prospino)" , "l");

  if(h_RA2b_best) this_leg->AddEntry(h_RA2b_best,"CMS SUS-11-006 (Inclusive jets with b-tag)" , "p");
  if(h_SS_best) this_leg->AddEntry(h_SS_best,"CMS SUS-11-010 (Same Sign dileptons)" , "p");

  this_leg->Draw();

  TLatex  t;
  //  t.SetTextAngle(45.);
  t.SetTextSize(0.03);
  t.SetTextAlign(33);

  t.DrawText(650, 50,"L = 1.1 / fb,  m(LSP)=50 (GeV)");
  //  t.DrawText(550,0.02,"L=1.1 /fb");
  //  t.DrawText(300, 0.85," CMS Preliminary");
  //  t.DrawText(300, 0.8, "#sqrt{s} = 7 TeV L=1.1 fb^{-1} ");

  cPlot1->SaveAs("/tmp/gluinoStopResult.gif");
  cPlot1->Print("/tmp/gluinoStopResult.eps");


}


void stopSummary() {

  TH1 *hRef=getHisto("reference_xSec_stop.root", "stop" , "0", 1,1);

  cout << " file1 " << endl;

  TH1 * h_RA4_loose=get1dHisto("T2tt_RA4_LS.root","h_RA4_T2tt_Loose_cloned", 50.);
  TH1 * h_RA4_tight=get1dHisto("T2tt_RA4_LS.root","h_RA4_T2tt_Tight_cloned", 50.);

  TH1 * h_RA4_LP=get1dHisto("RA4georgia/limit.root","limit_sms/hLimit_pl", 50.);

  TH1 * h_RA2b_1loose=get1dHisto("RA2bT2tt/T2tt_xsULs_finerScan.root","h2d_1bloose_xsUL", 50.);
  TH1 * h_RA2b_1tight=get1dHisto("RA2bT2tt/T2tt_xsULs_finerScan.root","h2d_1btight_xsUL", 50.);
  TH1 * h_RA2b_2loose=get1dHisto("RA2bT2tt/T2tt_xsULs_finerScan.root","h2d_2bloose_xsUL", 50.);
  TH1 * h_RA2b_2tight=get1dHisto("RA2bT2tt/T2tt_xsULs_finerScan.root","h2d_2btight_xsUL", 50.);

  TH1 * h_RA2_highHT=get1dHisto("RA2Seema/CLs_SMS_T2tt.root","highHT_obsLimit", 50.);
  TH1 * h_RA2_highMHT=get1dHisto("RA2Seema/CLs_SMS_T2tt.root","mediumHTMHT_obsLimit", 50.);
  TH1 * h_RA2_highHTMHT=get1dHisto("RA2Seema/CLs_SMS_T2tt.root","highHTMHT_obsLimit", 50.);

  TH1 * h_MT2b_T2tt=get1dHisto("T2tt_Leo/T2tt_MT2b.root","h_MT2_T2tt_cloned", 50.);

  TH1 * h_OS_T2tt_ht=get1dHisto("OSBen/OS_T2tt.root","hxsec_highht", 50.);
  TH1 * h_OS_T2tt_met=get1dHisto("OSBen/OS_T2tt.root","hxsec_highmet", 50.);

  TH1 * h_RA1=get1dHisto("RA1Ted/RA1_T2tt_xsLimit.root","UpperLimit_2D", 50.);

  cout << " file2 " << endl;
  
  TCanvas *cPlot1 = new TCanvas("cPlot1","cPlot1");
  cPlot1->SetLogy(1);
  cPlot1->SetGridy(1);
  cPlot1->SetGridx(1);

  //  hRef->SetRangeUser(300.,800);
  hRef->GetYaxis()->SetTitle("#sigma x BR (#tilde{t}#tilde{t} #rightarrow 2t+2LSP) (pb)");
  hRef->GetXaxis()->SetTitle("m_{#tilde{t}}(GeV)");
  //  hRef->GetXaxis()->SetRangeUser(300.,600.);
  hRef->GetXaxis()->SetRangeUser(200.,600.);

  hRef->SetMaximum(500);
  hRef->SetMinimum(0.0001);
  hRef->Draw("hist");

  h_RA4_loose->SetMarkerSize(1.5);
  h_RA4_loose->SetMarkerColor(30);
  h_RA4_loose->SetMarkerStyle(20);
  h_RA4_loose->Draw("p same");

  h_RA4_tight->SetMarkerSize(1.5);
  h_RA4_tight->SetMarkerColor(38);
  h_RA4_tight->SetMarkerStyle(20);
  h_RA4_tight->Draw("p same");

  h_RA4_LP->SetMarkerSize(1.5);
  h_RA4_LP->SetMarkerColor(40);
  h_RA4_LP->SetMarkerStyle(21);
  h_RA4_LP->Draw("p same");

  h_OS_T2tt_met->SetMarkerSize(1.5);
  h_OS_T2tt_met->SetMarkerColor(kMagenta);
  h_OS_T2tt_met->SetMarkerStyle(22);
  h_OS_T2tt_met->Draw("p same");

  h_OS_T2tt_ht->SetMarkerSize(1.5);
  h_OS_T2tt_ht->SetMarkerColor(kViolet-6);
  h_OS_T2tt_ht->SetMarkerStyle(22);
  h_OS_T2tt_ht->Draw("p same");

  h_RA2b_1loose->SetMarkerSize(1.5);
  h_RA2b_1loose->SetMarkerColor(42);
  h_RA2b_1loose->Draw("p same");

  h_RA2b_1loose->SetMarkerSize(1.5);
  h_RA2b_2loose->SetMarkerColor(46);
  h_RA2b_2loose->Draw("p same");

  h_RA2b_1tight->SetMarkerSize(1.5);
  h_RA2b_1tight->SetMarkerColor(48);
  h_RA2b_1tight->Draw("p same");

  h_RA2b_2tight->SetMarkerSize(1.5);
  h_RA2b_2tight->SetMarkerColor(44);
  h_RA2b_2tight->Draw("p same");

  h_MT2b_T2tt->SetMarkerSize(1.5);
  h_MT2b_T2tt->SetMarkerStyle(29);
  h_MT2b_T2tt->SetMarkerColor(12);
  h_MT2b_T2tt->Draw("p same");

  h_RA2_highHT->SetMarkerSize(1.5);
  h_RA2_highHT->SetMarkerColor(kRed);
  h_RA2_highHT->SetMarkerStyle(23);
  h_RA2_highHT->Draw("p same");

  h_RA2_highMHT->SetMarkerSize(1.5);
  h_RA2_highMHT->SetMarkerColor(kPink-1);
  h_RA2_highMHT->SetMarkerStyle(23);
  h_RA2_highMHT->Draw("p same");

  h_RA2_highHTMHT->SetMarkerSize(1.5);
  h_RA2_highHTMHT->SetMarkerColor(kRed-1);
  h_RA2_highHTMHT->SetMarkerStyle(23);
  h_RA2_highHTMHT->Draw("p same");

  h_RA1->SetMarkerSize(1.5);
  h_RA1->SetMarkerStyle(29);
  h_RA1->SetMarkerColor(1);
  h_RA1->Draw("p same");

  TLegend* this_leg = new TLegend(0.5,0.12,0.9,0.45);
  this_leg->SetFillColor(0);
  this_leg->SetBorderSize(0);
  this_leg->SetTextSize(0.035);

  if(hRef) this_leg->AddEntry(hRef,"(NLO Prospino)" , "l");

  if(h_RA4_loose) this_leg->AddEntry(h_RA4_loose,"1l+MET, LS Loose" , "p");
  if(h_RA4_tight) this_leg->AddEntry(h_RA4_tight,"1l+MET, LS Tight" , "p");

  if(h_RA4_LP) this_leg->AddEntry(h_RA4_LP,"1l+MET, LP" , "p");

  if(h_OS_T2tt_met) this_leg->AddEntry(h_OS_T2tt_met,"OS+MET" , "p");
  if(h_OS_T2tt_ht) this_leg->AddEntry(h_OS_T2tt_ht,"OS+HT" , "p");

  if(h_RA2b_1tight) this_leg->AddEntry(h_RA2b_1tight,"#geq 1b, Tight" , "p");
  if(h_RA2b_2tight) this_leg->AddEntry(h_RA2b_2tight,"#geq 2b, Tight" , "p");
  if(h_RA2b_1loose) this_leg->AddEntry(h_RA2b_1loose,"#geq 1b, Loose" , "p");
  if(h_RA2b_2loose) this_leg->AddEntry(h_RA2b_2loose,"#geq 2b, Loose" , "p");

  if(h_MT2b_T2tt) this_leg->AddEntry(h_MT2b_T2tt,"MT2" , "p");
  
  if(h_RA2_highHT) this_leg->AddEntry(h_RA2_highHT,"high HT" , "p");
  if(h_RA2_highMHT) this_leg->AddEntry(h_RA2_highMHT,"high MHT" , "p");
  if(h_RA2_highHTMHT) this_leg->AddEntry(h_RA2_highHTMHT,"high HT + high MHT" , "p");
  
  if(h_RA1) this_leg->AddEntry(h_RA1,"alphaT" , "p");

  this_leg->Draw();

  TText t;
  //  t.SetTextAngle(45.);
  t.SetTextSize(0.05);
  t.SetTextAlign(33);

  t.DrawText(350,0.01,"m(LSP)=50");
  t.DrawText(350,0.02,"L=1.1 /fb");
  //  t.DrawText(300, 0.85," CMS Preliminary");
  //  t.DrawText(300, 0.8, "#sqrt{s} = 7 TeV L=1.1 fb^{-1} ");

  cPlot1->SaveAs("/tmp/dalfonso/StopResult.gif");
  cPlot1->Print("/tmp/dalfonso/StopResult.eps");

  return;

}


void xSecParse() {
  
  //sg  0  0     0.0    0.0    1.0  300.0  150.0  0.000  579.     0.841E-03  825.     0.108E-02 1.4258  173.      246. 

  TGraph2D *gr2= new TGraph2D("/afs/cern.ch/user/d/dalfonso/xSec_gq.txt", "%lg %lg %lg");

  TH2D*	gr= gr2->GetHistogram();

  TH2D *g=clone2D(gr);

  TH2D *empty = new TH2D("empty","",60,0.,1500.,60,0.,1500.);
  //  empty->SetMaximum(1);
  //  empty->SetYTitle("Arbitrary Units");
  empty->SetXTitle("m(Gluino)");
  empty->SetYTitle("m(SQ)");
  empty->SetZTitle("xSec (pb-1)");
  empty->SetTitle("");
  
  TCanvas *cPlot = new TCanvas("cPlot","cPlot");
  cPlot->SetLogz(1);

  empty->Draw("hist");
  if(gr2) gr2->Draw("colz same");

  TFile f1("TGQ_xSec.root","UPDATE");
  //  Stop10->Write();
  g->Write();
  f1.Close();

}

void LHEparse() {
  
  TGraph *gr1= new TGraph("/data/dalfonso/LHE/list_TGQ_0.0.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr2= new TGraph("/data/dalfonso/LHE/list_TGQ_0.2.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr3= new TGraph("/data/dalfonso/LHE/list_TGQ_0.4.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr4= new TGraph("/data/dalfonso/LHE/list_TGQ_0.6.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr5= new TGraph("/data/dalfonso/LHE/list_TGQ_0.8.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr6= new TGraph("/data/dalfonso/LHE/list_T1tttt.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr7= new TGraph("/data/dalfonso/LHE/list_T2bw_T2bw_0.25.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr8= new TGraph("/data/dalfonso/LHE/list_T2bw_T2bw_0.75.txt", "%*s %*s %*s %*s %*s %lg %lg");
  TGraph *gr9= new TGraph("/data/dalfonso/LHE/list_T2bw_T2bw_0.5.txt", "%*s %*s %*s %*s %*s %lg %lg");

  //  TCanvas *c1 = new TCanvas("cPlot1","cPlot1");
  
  if(gr1) gr1->SetMarkerColor(1);
  if(gr2) gr2->SetMarkerColor(2);
  if(gr3) gr3->SetMarkerColor(3);
  if(gr4) gr4->SetMarkerColor(4);
  if(gr5) gr5->SetMarkerColor(5);
  if(gr6) gr5->SetMarkerColor(6);

  if(gr1) gr1->SetMarkerStyle(1);
  if(gr2) gr2->SetMarkerStyle(2);
  if(gr3) gr3->SetMarkerStyle(3);
  if(gr4) gr4->SetMarkerStyle(4);
  if(gr5) gr5->SetMarkerStyle(5);
  if(gr6) gr5->SetMarkerStyle(6);
  
  /*
    int nBinY= histo->GetNbinsY();
    double yLow=histo->GetBinLowEdge(1); 
    double yHigh=histo->GetBinLowEdge(nBinY+1); 
  */

  TH2D *empty = new TH2D("empty","",60,0.,1500.,60,0.,1500.);
  empty->SetMaximum(1);
  //  empty->SetYTitle("Arbitrary Units");
  empty->SetXTitle("m(Gluino)");
  empty->SetYTitle("m(Gluino)");
  empty->SetTitle("");
  
  TCanvas *cPlot = new TCanvas("cPlot","cPlot");
  //  cPlot->SetLogy(1);

  empty->Draw("hist");
  if(gr9) gr9->Draw("A* same");

  //  if(gr2) gr2->Draw("A* same");
  //  if(gr3) gr3->Draw("A* same");
  //  if(gr4) gr4->Draw("A* same");
  //  if(gr5) gr5->Draw("A* same");

  /*
    TLegend* this_leg = new TLegend(0.2,0.65,0.4,0.73);
    this_leg->SetFillColor(0);
    this_leg->SetBorderSize(0);
    this_leg->SetTextSize(0.035);
    if(gr1) this_leg->AddEntry(gr1,"x=0.0" , "p");
    if(gr2) this_leg->AddEntry(gr2,"x=0.2" , "p");
    if(gr3) this_leg->AddEntry(gr3,"x=0.4" , "p");
    if(gr4) this_leg->AddEntry(gr4,"x=0.6" , "p");
    if(gr5) this_leg->AddEntry(gr5,"x=0.8" , "p");
    this_leg->Draw();
  */
  
  //  cPlot->SaveAs("T2bw_0.25.gif");
  // cPlot->SaveAs("T2bw_0.75.gif");
  cPlot->SaveAs("T2bw_0.5.gif");
  
}

void xSec() {

 Double_t xSt[186]={

    150.0,
    160.0,
    170.0,
    180.0,
    190.0,
    200.0,
    210.0,
    220.0,
    230.0,
    240.0,
    250.0,
    260.0,
    270.0,
    280.0,
    290.0,
    300.0,
    310.0,
    320.0,
    330.0,
    340.0,
    350.0,
    360.0,
    370.0,
    380.0,
    390.0,
    400.0,
    410.0,
    420.0,
    430.0,
    440.0,
    450.0,
    460.0,
    470.0,
    480.0,
    490.0,
    500.0,
    510.0,
    520.0,
    530.0,
    540.0,
    550.0,
    560.0,
    570.0,
    580.0,
    590.0,
    600.0,
    610.0,
    620.0,
    630.0,
    640.0,
    650.0,
    660.0,
    670.0,
    680.0,
    690.0,
    700.0,
    710.0,
    720.0,
    730.0,
    740.0,
    750.0,
    760.0,
    770.0,
    780.0,
    790.0,
    800.0,
    810.0,
    820.0,
    830.0,
    840.0,
    850.0,
    860.0,
    870.0,
    880.0,
    890.0,
    900.0,
    910.0,
    920.0,
    930.0,
    940.0,
    950.0,
    960.0,
    970.0,
    980.0,
    990.0,
    1000.0,
    1010.0,
    1020.0,
    1030.0,
    1040.0,
    1050.0,
    1060.0,
    1070.0,
    1080.0,
    1090.0,
    1100.0,
    1110.0,
    1120.0,
    1130.0,
    1140.0,
    1150.0,
    1160.0,
    1170.0,
    1180.0,
    1190.0,
    1200.0,
    1210.0,
    1220.0,
    1230.0,
    1240.0,
    1250.0,
    1260.0,
    1270.0,
    1280.0,
    1290.0,
    1300.0,
    1310.0,
    1320.0,
    1330.0,
    1340.0,
    1350.0,
    1360.0,
    1370.0,
    1380.0,
    1390.0,
    1400.0,
    1410.0,
    1420.0,
    1430.0,
    1440.0,
    1450.0,
    1460.0,
    1470.0,
    1480.0,
    1490.0,
    1500.0,
    1510.0,
    1520.0,
    1530.0,
    1540.0,
    1550.0,
    1560.0,
    1570.0,
    1580.0,
    1590.0,
    1600.0,
    1610.0,
    1620.0,
    1630.0,
    1640.0,
    1650.0,
    1660.0,
    1670.0,
    1680.0,
    1690.0,
    1700.0,
    1710.0,
    1720.0,
    1730.0,
    1740.0,
    1750.0,
    1760.0,
    1770.0,
    1780.0,
    1790.0,
    1800.0,
    1810.0,
    1820.0,
    1830.0,
    1840.0,
    1850.0,
    1860.0,
    1870.0,
    1880.0,
    1890.0,
    1900.0,
    1910.0,
    1920.0,
    1930.0,
    1940.0,
    1950.0,
    1960.0,
    1970.0,
    1980.0,
    1990.0,
    2000.0
 };


 Double_t ySt[186]={

    35.0,
    25.0,
    18.1,
    13.4,
    9.98,
    7.55,
    5.76,
    4.45,
    3.46,
    2.72,
    2.15,
    1.71,
    1.37,
    1.11,
    0.897,
    0.731,
    0.599,
    0.493,
    0.408,
    0.338,
    0.282,
    0.236,
    0.198,
    0.166,
    0.140,
    0.119,
    0.101,
    0.860E-01,
    0.734E-01,
    0.628E-01,
    0.538E-01,
    0.462E-01,
    0.398E-01,
    0.343E-01,
    0.297E-01,
    0.257E-01,
    0.223E-01,
    0.194E-01,
    0.169E-01,
    0.147E-01,
    0.128E-01,
    0.112E-01,
    0.980E-02,
    0.859E-02,
    0.753E-02,
    0.662E-02,
    0.582E-02,
    0.512E-02,
    0.451E-02,
    0.398E-02,
    0.351E-02,
    0.311E-02,
    0.275E-02,
    0.243E-02,
    0.216E-02,
    0.191E-02,
    0.170E-02,
    0.151E-02,
    0.134E-02,
    0.119E-02,
    0.106E-02,
    0.944E-03,
    0.841E-03,
    0.750E-03,
    0.669E-03,
    0.597E-03,
    0.534E-03,
    0.477E-03,
    0.426E-03,
    0.381E-03,
    0.341E-03,
    0.306E-03,
    0.274E-03,
    0.245E-03,
    0.220E-03,
    0.197E-03,
    0.177E-03,
    0.159E-03,
    0.142E-03,
    0.128E-03,
    0.115E-03,
    0.103E-03,
    0.928E-04,
    0.835E-04,
    0.751E-04,
    0.675E-04,
    0.608E-04,
    0.547E-04,
    0.492E-04,
    0.443E-04,
    0.399E-04,
    0.359E-04,
    0.324E-04,
    0.292E-04,
    0.263E-04,
    0.237E-04,
    0.214E-04,
    0.193E-04,
    0.174E-04,
    0.157E-04,
    0.141E-04,
    0.127E-04,
    0.115E-04,
    0.104E-04,
    0.935E-05,
    0.844E-05,
    0.762E-05,
    0.687E-05,
    0.620E-05,
    0.560E-05,
    0.505E-05,
    0.456E-05,
    0.411E-05,
    0.371E-05,
    0.335E-05,
    0.302E-05,
    0.273E-05,
    0.246E-05,
    0.222E-05,
    0.201E-05,
    0.181E-05,
    0.163E-05,
    0.147E-05,
    0.133E-05,
    0.120E-05,
    0.108E-05,
    0.975E-06,
    0.880E-06,
    0.793E-06,
    0.715E-06,
    0.645E-06,
    0.582E-06,
    0.524E-06,
    0.473E-06,
    0.426E-06,
    0.384E-06,
    0.346E-06,
    0.313E-06,
    0.282E-06,
    0.255E-06,
    0.230E-06,
    0.207E-06,
    0.187E-06,
    0.168E-06,
    0.152E-06,
    0.137E-06,
    0.123E-06,
    0.111E-06,
    0.100E-06,
    0.901E-07,
    0.811E-07,
    0.730E-07,
    0.657E-07,
    0.591E-07,
    0.532E-07,
    0.478E-07,
    0.430E-07,
    0.386E-07,
    0.347E-07,
    0.312E-07,
    0.280E-07,
    0.251E-07,
    0.226E-07,
    0.202E-07,
    0.181E-07,
    0.163E-07,
    0.146E-07,
    0.131E-07,
    0.117E-07,
    0.105E-07,
    0.938E-08,
    0.839E-08,
    0.750E-08,
    0.671E-08,
    0.600E-08,
    0.536E-08,
    0.478E-08,
    0.427E-08,
    0.381E-08,
    0.340E-08,
    0.303E-08,
    0.270E-08,
    0.241E-08,
    0.214E-08,
    0.191E-08,
    0.170E-08,
 };

  Double_t ySt2[186]={
    53.1,
    38.2,
    27.9,
    20.7,
    15.6,
    11.8,
    9.10,
    7.06,
    5.53,
    4.36,
    3.47,
    2.78,
    2.24,
    1.81,
    1.47,
    1.20,
    0.995,
    0.823,
    0.683,
    0.566,
    0.474,
    0.398,
    0.335,
    0.285,
    0.241,
    0.205,
    0.176,
    0.150,
    0.129,
    0.111,
    0.952E-01,
    0.821E-01,
    0.710E-01,
    0.615E-01,
    0.534E-01,
    0.464E-01,
    0.404E-01,
    0.353E-01,
    0.308E-01,
    0.270E-01,
    0.236E-01,
    0.208E-01,
    0.182E-01,
    0.159E-01,
    0.141E-01,
    0.124E-01,
    0.110E-01,
    0.968E-02,
    0.857E-02,
    0.761E-02,
    0.674E-02,
    0.598E-02,
    0.532E-02,
    0.473E-02,
    0.425E-02,
    0.379E-02,
    0.338E-02,
    0.302E-02,
    0.270E-02,
    0.241E-02,
    0.216E-02,
    0.193E-02,
    0.173E-02,
    0.155E-02,
    0.139E-02,
    0.125E-02,
    0.112E-02,
    0.100E-02,
    0.906E-03,
    0.811E-03,
    0.730E-03,
    0.657E-03,
    0.593E-03,
    0.534E-03,
    0.482E-03,
    0.435E-03,
    0.393E-03,
    0.355E-03,
    0.321E-03,
    0.290E-03,
    0.262E-03,
    0.238E-03,
    0.215E-03,
    0.195E-03,
    0.176E-03,
    0.160E-03,
    0.145E-03,
    0.131E-03,
    0.119E-03,
    0.108E-03,
    0.982E-04,
    0.890E-04,
    0.808E-04,
    0.733E-04,
    0.666E-04,
    0.604E-04,
    0.550E-04,
    0.500E-04,
    0.454E-04,
    0.413E-04,
    0.376E-04,
    0.342E-04,
    0.311E-04,
    0.283E-04,
    0.258E-04,
    0.235E-04,
    0.214E-04,
    0.195E-04,
    0.177E-04,
    0.161E-04,
    0.147E-04,
    0.134E-04,
    0.122E-04,
    0.111E-04,
    0.101E-04,
    0.921E-05,
    0.839E-05,
    0.766E-05,
    0.696E-05,
    0.635E-05,
    0.579E-05,
    0.528E-05,
    0.481E-05,
    0.438E-05,
    0.399E-05,
    0.364E-05,
    0.332E-05,
    0.303E-05,
    0.276E-05,
    0.252E-05,
    0.230E-05,
    0.209E-05,
    0.191E-05,
    0.173E-05,
    0.159E-05,
    0.144E-05,
    0.132E-05,
    0.120E-05,
    0.110E-05,
    0.100E-05,
    0.914E-06,
    0.837E-06,
    0.764E-06,
    0.697E-06,
    0.635E-06,
    0.579E-06,
    0.529E-06,
    0.483E-06,
    0.441E-06,
    0.402E-06,
    0.366E-06,
    0.335E-06,
    0.307E-06,
    0.280E-06,
    0.259E-06,
    0.232E-06,
    0.216E-06,
    0.196E-06,
    0.175E-06,
    0.159E-06,
    0.145E-06,
    0.132E-06,
    0.120E-06,
    0.110E-06,
    0.100E-06,
    0.911E-07,
    0.830E-07,
    0.754E-07,
    0.686E-07,
    0.624E-07,
    0.567E-07,
    0.516E-07,
    0.485E-07,
    0.425E-07,
    0.388E-07,
    0.352E-07,
    0.320E-07,
    0.290E-07,
    0.264E-07,
    0.239E-07,
    0.217E-07,
    0.197E-07,
    0.184E-07,
    0.162E-07,
    0.147E-07,
    0.134E-07

  };

  Double_t xG[38]={
    100.0,
    125.0,
    150.0,
    175.0,
    200.0,
    225.0,
    250.0,
    275.0,
    300.0,
    325.0,
    350.0,
    375.0,
    400.0,
    425.0,
    450.0,
    475.0,
    50.0, // this is the little one
    500.0,
    525.0,
    550.0,
    575.0,
    600.0,
    625.0,
    650.0,
    675.0,
    700.0,
    725.0,
    75.0,
    750.0,
    775.0,
    800.0,
    825.0,
    850.0,
    875.0,
    900.0,
    925.0,
    950.0, // this is the highest
  };

  Double_t yG[38]={
    0.212E+05,
    0.717E+04,
    0.286E+04,
    0.128E+04,
    625.,
    326.,
    180.,
    104.,
    62.1,
    38.3,
    24.2,
    15.7,
    10.4,
    6.97,
    4.76,
    3.30,
    0.484E+06,
    2.31,
    1.64,
    1.17,
    0.847,
    0.617,
    0.453,
    0.335,
    0.249,
    0.186,
    0.140,
    0.810E+05,
    0.106,
    0.801E-01,
    0.610E-01,
    0.466E-01,
    0.358E-01,
    0.275E-01,
    0.212E-01,
    0.164E-01,
    0.128E-01
  };


  Double_t xSq[18]={
    100.0,
    150.0,
    200.0,
    250.0,
    300.0,
    350.0,
    400.0,
    450.0,
    500.0,
    550.0,
    600.0,
    650.0,
    700.0,
    750.0,
    800.0,
    850.0,
    900.0,
  };


  Double_t ySq[18]={
    0.393E+04,
    547.,
    123.,
    36.8,
    13.3,
    5.48,
    2.49,
    1.22,
    0.630,
    0.342,
    0.193,
    0.112,
    0.666E-01,
    0.404E-01,
    0.249E-01,
    0.156E-01,
    0.981E-02,
  };


  TH1D * Gluino= new TH1D("gluino","gluino pair production",40, 0., 1000.);
  TH1D * Squark= new TH1D("squark","squark pair production",40, 0., 1000.);
  TH1D * StopLO= new TH1D("stopLO","stop pair production",200, 0., 2000.);
  TH1D * Stop= new TH1D("stop","stop pair production",200, 0., 2000.);

  /*
  for (Int_t i=1;i<38+1;i++) {
    Gluino->SetBinContent(Gluino->FindBin(xG[i-1]),yG[i-1]);
  }

  for (Int_t i=1;i<18+1;i++) {
    Squark->SetBinContent(Squark->FindBin(xSq[i-1]),ySq[i-1]);
  }
  */

  for (Int_t i=1;i<186+1;i++) {
    StopLO->SetBinContent(StopLO->FindBin(xSt[i-1]),ySt[i-1]);
  }

  for (Int_t i=1;i<186+1;i++) {
    Stop->SetBinContent(Stop->FindBin(xSt[i-1]),ySt2[i-1]);
  }

  Squark->SetYTitle("pb");
  Gluino->SetYTitle("pb");
  StopLO->SetYTitle("pb");
  Stop->SetYTitle("pb");

  StopLO->SetXTitle("m_{Stop}");
  Stop->SetXTitle("m_{Stop}");
  Squark->SetXTitle("m_{Squark}");
  Gluino->SetXTitle("m_{Gluino}");

  Squark->SetTitle("");
  Gluino->SetTitle("");
  StopLO->SetTitle("");
  Stop->SetTitle("");

  TFile f("reference_xSec_stop.root","UPDATE");
  //  Stop10->Write();
  Stop->Write();
  f.Close();

  //  return ;

}

TH2F * restyleHisto(TH2F * h1, bool doLimit, bool doError, string topology) {

  cout << "getting Name " << h1->GetName() << endl;  

  TH2F *h=clone2F(h1);
  TString name = TString(h->GetName());

  /*
  if(name.Contains("RA2b") && name.Contains("T1bbbb") && name.Contains("best")) {
    TFile f("test_RA2b.root","UPDATE");
    h1->Write();
    h->Write();
    f.Close();
  }
  */

  //  cout << "out of clone" << endl;

  if(doLimit)  h->SetMaximum(4);
  if(doLimit)  h->SetMinimum(0.01);
  if(!doLimit && doError)  h->SetMaximum(0.5);
  if(!doLimit && doError)  h->SetMinimum(0.01);
  if(!doLimit && !doError && (topology=="T1lnu"))  h->SetMaximum(0.20);
  if(!doLimit && !doError && (topology=="T1lh"))  h->SetMaximum(0.50);
  if(!doLimit && !doError && ((topology=="T1lnu") || (topology=="T5zz") || (topology=="T1lh")))  h->SetMinimum(0.0);
  if(!doLimit && !doError && ((topology=="T1") || (topology=="T2") || (topology=="T1bbbb")))  h->SetMaximum(0.50);
  if(!doLimit && !doError && ((topology=="T1") || (topology=="T2") || (topology=="T1bbbb")))  h->SetMinimum(0.0);
  if((topology=="T5zz") && (name.Contains("warren")) && (!doLimit) && (!doError)) h->SetMaximum(0.40);
  if(!doLimit && !doError && ((topology=="T1tttt")))  h->SetMaximum(0.04);

  h->SetTitle("");
  h->SetXTitle("m_{#tilde{g}} (GeV)");
  if(topology=="T2") h->SetXTitle("m_{#tilde{q}} (GeV)");
  if(topology=="T2bb") h->SetXTitle("m_{#tilde{b}} (GeV)");
  if(topology=="T2tt") h->SetXTitle("m_{#tilde{t}} (GeV)");

  h->SetYTitle("m_{#chi^{0}} (GeV)");
  h->SetTitle("");
  if(doLimit) h->SetZTitle("95% CL upper limit on #sigma (pb) (CL_{s})");
  if(!doLimit && doError) h->SetZTitle("Error(A #times #varepsilon)/(A #times #varepsilon)");
  if(!doLimit && !doError) h->SetZTitle("A #times #varepsilon");
  if((topology=="T5zz") && !doLimit && !doError) h->SetZTitle("A #times #varepsilon (#geq 1 Z(ll))");

  //  if(((topology=="T1") || (topology=="T2") || (topology=="T1bbbb"))) 
  h->GetXaxis()->SetRangeUser(50.,1200.);
  h->GetYaxis()->SetRangeUser(50.,1200.);
  if((topology=="T1bbbb") || (topology=="T1") || (topology=="T2")) h->GetXaxis()->SetRangeUser(300.,1200.);
  //  if((topology=="T1bbbb") || (topology=="T1") || (topology=="T2")) h->GetXaxis()->SetLabelSize(0.038);
  if((topology=="T1bbbb") || (topology=="T1tttt") || (topology=="T1") || (topology=="T2") || (topology=="T1lh") || (topology=="T2tt") ) h->GetXaxis()->SetNdivisions(505.);

  if((topology=="T1tttt")) h->GetXaxis()->SetRangeUser(400.,1200.);

  if((topology=="TGQ")) h->GetXaxis()->SetRangeUser(400.,1500.);
  if((topology=="TGQ")) h->GetYaxis()->SetRangeUser(450.,1500.);

  if(name.Contains("OS_2010_HT")) h->GetXaxis()->SetRangeUser(300.,1200.);
  if(name.Contains("OS_2010_MET")) h->GetXaxis()->SetRangeUser(350.,1200.);



  //comp_T1lnu_OS_2010_HT

  h->GetZaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetTitleOffset(1.2);

  //  gROOT->ForceStyle();

  return h;

}

//void labelling(TPad * pad , TH2F *histo,char * selection,char * leg, char * topology, bool doLimit){
void labelling(TH2F *histo,char * selection,char * leg, string topology, bool doLimit){

  TString name = TString(histo->GetName());

  cout << "name " << name.Data() << endl;

  TLatex latexLabel;
  latexLabel.SetTextSize(0.035);
  latexLabel.SetNDC();

  latexLabel.DrawLatex(0.2, 0.85," CMS Preliminary");
  if((topology=="T1") || (topology=="T2")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=1.1 fb^{-1} ");
  if((topology=="T1bbbb") || (topology=="T1tttt") || (topology=="T2tt") || (topology=="T2bb")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=1.1 fb^{-1} ");
  if((topology=="T1lnu")|| (topology=="T1lh")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=0.98 fb^{-1} ");
  if(name.Contains("ANN")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=1.15 fb^{-1} ");
  if((topology=="T5zz") && name.Contains("jzb")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=1.6 fb^{-1} ");
  //  if((topology=="T5zz") && name.Contains("jzb")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=0.191 fb^{-1} ");
  if((topology=="T5zz") && name.Contains("warren")) latexLabel.DrawLatex(0.2, 0.8, "#sqrt{s} = 7 TeV L=0.98 fb^{-1} ");
  latexLabel.DrawLatex(0.2, 0.75, selection);

  //  if((topology=="T5zz") && name.Contains("warren"))  cout << "ciao ciao warren" << endl;
  //  if((topology=="T5zz") && name.Contains("jzb"))  cout << "ciao ciao jzb " << endl;

  //  latexLabel.SetTextSize(0.03);
  latexLabel.DrawLatex(0.1, 0.91, leg);

  TLine *line;
  TLine *line2;
  line = new TLine(50., 50., 1200., 1200.);
  if(topology=="T5zz") line = new TLine(50.+92., 50., 1200., 1200-92.);
  if(topology=="T2tt") line = new TLine(50.+173., 50., 1200., 1200-173.);
  if(topology=="T2tt") line2 = new TLine(50., 50., 1200., 1200.);

  if((topology=="T1tttt")) line = new TLine(400., 400., 1200., 1200.);
  if((topology=="T1tttt")) line2 = new TLine(50.+350., 50., 1200., 1200-350.);

  if((topology=="T1bbbb") || (topology=="T1") || (topology=="T2")) line = new TLine(300., 300., 1200., 1200.);

  if((topology=="T1lh") && name.Contains("OS_2010_HT")) line = new TLine(300., 300., 1200., 1200.);
  if((topology=="T1lh") && name.Contains("OS_2010_MET")) line = new TLine(350., 350., 1200., 1200.);

  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("same");

  if(topology=="T2tt") line2->SetLineStyle(2);
  if(topology=="T2tt") line2->SetLineWidth(2);
  if(topology=="T2tt") line2->Draw("same");
 
  if(topology=="T1tttt") line2->SetLineStyle(2);
  if(topology=="T1tttt") line2->SetLineWidth(2);
  if(topology=="T1tttt") line2->Draw("same"); 

  TLine *lineGQ02;
  TLine *lineGQ04;
  TLine *lineGQ08;

  lineGQ02 = new TLine(50., 0.2*50., 1200., 0.2*1200.);
  lineGQ04 = new TLine(50., 0.4*50., 1200., 0.4*1200.);
  lineGQ08 = new TLine(50., 0.8*50., 1200., 0.8*1200.);

  lineGQ02->SetLineColor(2);
  lineGQ02->SetLineStyle(2);
  lineGQ02->SetLineWidth(5);
  //  lineGQ02->Draw("same");

  lineGQ04->SetLineColor(2);
  lineGQ04->SetLineStyle(2);
  lineGQ04->SetLineWidth(5);
  //  lineGQ04->Draw("same");

  lineGQ08->SetLineColor(2);
  lineGQ08->SetLineStyle(2);
  lineGQ08->SetLineWidth(5);
  //  lineGQ08->Draw("same");

  /*
  if(name.Contains("T1") && name.Contains("RA2_"))  lineGQ02->Draw("same");
  if(name.Contains("T2") && name.Contains("RA2_"))  lineGQ02->Draw("same");

  if(name.Contains("T1") && name.Contains("RA2_"))  lineGQ04->Draw("same");
  if(name.Contains("T2") && name.Contains("RA2_"))  lineGQ04->Draw("same");

  if(name.Contains("T1") && name.Contains("RA2_"))  lineGQ08->Draw("same");
  if(name.Contains("T2") && name.Contains("RA2_"))  lineGQ08->Draw("same");  
  
  */
  /*
    TText t;
    t.SetTextAngle(45.);
    t.SetTextSize(0.03);
    t.SetTextAlign(33);
    double x1=histo->GetXaxis()->GetBinCenter(histo->GetXaxis()->FindBin(500))+300;
    double y1=histo->GetYaxis()->GetBinCenter(histo->GetYaxis()->FindBin(500))+320;
    
    if(!(topology=="T5zz"))t.DrawText(x1,y1,"M(G) > mLSP");
    if(topology=="T5zz") t.DrawText(x1,y1,"M(Q) > mLSP+92 GeV");
  */
  
  //  return;
  if(!doLimit) return;
  gStyle->SetOptStat(0);
  cout << "before here " << endl;

  TGraph * refGraph_histo_1=getRefXsecGraph(histo,topology,1.);

  cout << "I'm here 1"<< endl;

  TGraph * refGraph_histo_13=getRefXsecGraph(histo,topology,1./3);
  cout << "I'm here 13 "<< endl;

  TGraph * refGraph_histo_3=getRefXsecGraph(histo,topology,3.);

  cout << "I'm here 3 "<< endl;

  TGraph * refGraph_histo_2=getRefXsecGraph(histo,topology,2.);

  cout << "I'm here 2 "<< endl;

  /*
  if(name.Contains("RA2b_best") || name.Contains("T1bbbb") || (name.Contains("alphaT") && name.Contains("T1") ) ) {
    
    TFile f("myFileGraph.root","UPDATE");
    if(refGraph_histo_1) refGraph_histo_1->Write();
    //    if(refGraph_histo_13) refGraph_histo_13->Write();
    //    if(refGraph_histo_3) refGraph_histo_3->Write();
    f.Close();

  }
  */
 
  if(refGraph_histo_1) refGraph_histo_1->SetLineWidth(3);
  if(refGraph_histo_13) refGraph_histo_13->SetLineWidth(3);
  if(refGraph_histo_3)  refGraph_histo_3->SetLineWidth(3);

  if(refGraph_histo_1) refGraph_histo_1->SetLineStyle(1);
  if(refGraph_histo_13) refGraph_histo_13->SetLineStyle(2);
  if(refGraph_histo_3) refGraph_histo_3->SetLineStyle(3);

  if(refGraph_histo_1) refGraph_histo_1->Draw("L same");
  if(refGraph_histo_13) refGraph_histo_13->Draw("L same");
  if(refGraph_histo_3) refGraph_histo_3->Draw("L same");

  //  TGraph * refGraph_histo_2;
  
  /*
  if((topology=="T2tt") || (topology=="T1tttt")) {
    
    refGraph_histo_2=getRefXsecGraph(histo,topology,2.);
    
    if(refGraph_histo_2)  refGraph_histo_2->SetLineWidth(3);
    if(refGraph_histo_2) refGraph_histo_2->SetLineStyle(4);
    if(refGraph_histo_2) refGraph_histo_2->Draw("L same");
    
  }
  
*/
  /*
    TGraph * refGraph_histo_10;
    
    if((topology=="T2tt") || (topology=="T1tttt")) {
    
    refGraph_histo_10=getRefXsecGraph(histo,topology,10.);
    
    if(refGraph_histo_10)  refGraph_histo_10->SetLineWidth(3);
    if(refGraph_histo_10) refGraph_histo_10->SetLineStyle(4);
    if(refGraph_histo_10) refGraph_histo_10->Draw("L same");    

    }

  */

  //  cout << "here 1"<< endl;
  
  //  if(refGraph_histo_1) refGraph_histo_1->SetLineStyle(1);
  
  //  cout << "here 2"<< endl;
  
  TLegend* this_leg = new TLegend(0.2,0.65,0.4,0.73);
  this_leg->SetFillColor(0);
  this_leg->SetBorderSize(0);
  this_leg->SetTextSize(0.035);
  if(refGraph_histo_1) this_leg->AddEntry(refGraph_histo_1,"#sigma^{prod} = #sigma^{NLO-QCD}" , "l");
  if(refGraph_histo_3) this_leg->AddEntry(refGraph_histo_3,"#sigma^{prod} = 3 #times #sigma^{NLO-QCD}" , "l");
  //  if((topology!="T1tttt") && (refGraph_histo_13)) this_leg->AddEntry(refGraph_histo_13,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}" , "l");
  if((refGraph_histo_13)) this_leg->AddEntry(refGraph_histo_13,"#sigma^{prod} = 1/3 #times #sigma^{NLO-QCD}" , "l");
  if(refGraph_histo_2 && (topology=="T1tttt")) this_leg->AddEntry(refGraph_histo_2,"#sigma^{prod} = 2 #times #sigma^{NLO-QCD}" , "l");
  //  if(refGraph_histo_10) this_leg->AddEntry(refGraph_histo_10,"#sigma^{prod} = 10 #times #sigma^{NLO-QCD}" , "l");
  this_leg->Draw();

  return;

  if(topology=="T1lnu" || topology=="T1") {

    //    char * filename="";
    //    if(topology=="T1lnu") filename="myFileT1lnu.root";
    //    if(topology=="T1bbbb") filename="myFileT1bbbb.root";

    TCanvas *tempGraph = new TCanvas();
    tempGraph->cd();
    
    if(refGraph_histo_1) refGraph_histo_1->Draw("L");
    
    TFile f("myFileGraph.root","UPDATE");
    if(refGraph_histo_1) refGraph_histo_1->Write();
    f.Close();

    delete tempGraph;
    
  }

  return;

}

struct lesserFirst {
 bool operator()(const pair<double,int>& t1, const pair<double,int>& t2) const { return t1.first < t2.first; }
};

vector<TH2F*> takeMinimum(vector<TH2F*> objVecExp,vector<TH2F*> objVecObs) {

  /////

  int nBinX= objVecExp.at(0)->GetNbinsX();
  double xLow=objVecExp.at(0)->GetBinLowEdge(1);
  double xHigh=objVecExp.at(0)->GetBinLowEdge(nBinX+1);
  
  int nBinY= objVecExp.at(0)->GetNbinsY();
  double yLow=objVecExp.at(0)->GetBinLowEdge(1); 
  double yHigh=objVecExp.at(0)->GetBinLowEdge(nBinY+1); 
  
  TH2F *crossSec = new TH2F("crossSec","crossSec",nBinX,xLow,xHigh,nBinY,yLow,yHigh);
  TH2F *analysis = new TH2F("crossSec","crossSec",nBinX,xLow,xHigh,nBinY,yLow,yHigh);

  for(int i=1; i <= objVecExp.at(0)->GetNbinsX(); ++i) {
    for(int j=1; j <= objVecExp.at(0)->GetNbinsY(); ++j) {
      
      if(objVecExp.at(0)->GetBinContent(i,j)==0) continue;

      int size=objVecExp.size();
      vector<pair<double,int> > myvector;
      pair<double,int> mypair;

      for(int k=0; k < size; ++k) {
	if(objVecExp.at(k)->GetBinContent(i,j)==0) continue;

	mypair.first=objVecExp.at(k)->GetBinContent(i,j);
	mypair.second=k;

	myvector.push_back(mypair);
      }
     
      //      myvector.resize(size);

      std::sort(myvector.begin(), myvector.end(), lesserFirst());
     
      pair<double,int> myFirstPair=myvector.at(0);       

      // this will return the analyiss that got the better expected limit
      int indexOfMin=myFirstPair.second;

      //
      double minExp=myFirstPair.first;

      // this will return the minimum of the Observed limit based on the expected
    
      TString name = TString(objVecExp.at(0)->GetName());
      
      double mGL=objVecObs.at(indexOfMin)->GetBinCenter(i);
      double mLSP=objVecObs.at(indexOfMin)->GetBinCenter(j);
      
      if(name.Contains("T2") && name.Contains("RA2_") && fabs(mGL-mLSP)<250) indexOfMin=2;
      if(name.Contains("T1") && name.Contains("RA2_") && fabs(mGL-mLSP)<300 && mLSP>300) indexOfMin=2;
      if(name.Contains("T1") && name.Contains("RA2_") && fabs(mGL)>900 && indexOfMin==0) indexOfMin=1;

      double min=objVecObs.at(indexOfMin)->GetBinContent(i,j);

      crossSec->Fill(objVecExp.at(0)->GetXaxis()->GetBinCenter(i),objVecExp.at(0)->GetYaxis()->GetBinCenter(j),min);
      analysis->Fill(objVecExp.at(0)->GetXaxis()->GetBinCenter(i),objVecExp.at(0)->GetYaxis()->GetBinCenter(j),indexOfMin+1); // +1 to remove the white
      
    }
  }

  vector<TH2F*> objvec;
  objvec.push_back(crossSec);
  objvec.push_back(analysis);
  //  cout<<"objvec.size() = "<<objvec.size()<<endl;

  return objvec;
  
}


void printCanvas(TH2F * histo,char * histoName, char * selection,char *legT5zz,char * label,bool doLimit,bool doError,bool doLogy){

  TCanvas *temp = new TCanvas();
  temp->cd();

  if(doLogy)  gPad->SetLogz(1);

  gStyle->SetCanvasDefH(400);
  gStyle->SetCanvasDefW(400);

  histo->Draw("colz");
  labelling(histo,selection,legT5zz,label,doLimit);

  char limitPNG[1000];
  char limitEPS[1000];
  sprintf(limitPNG,"RESULT/h_limit_%s.png",histoName);
  sprintf(limitEPS,"RESULT/h_limit_%s.eps",histoName);

  char effPNG[1000];
  char effEPS[1000];
  sprintf(effPNG,"RESULT/h_eff_%s.png",histoName);
  sprintf(effEPS,"RESULT/h_eff_%s.eps",histoName);

  char effStatErrPNG[1000];
  char effStatErrEPS[1000];
  sprintf(effStatErrPNG,"RESULT/h_effStatErr_%s.png",histoName);
  sprintf(effStatErrEPS,"RESULT/h_effStatErr_%s.eps",histoName);

  if(doLimit) temp->SaveAs(limitPNG);
  if(doLimit) temp->SaveAs(limitEPS);

  if(!doLimit && !doError) temp->SaveAs(effPNG);
  if(!doLimit && !doError) temp->SaveAs(effEPS);

  if(!doLimit && doError) temp->SaveAs(effStatErrPNG);
  if(!doLimit && doError) temp->SaveAs(effStatErrEPS);
  
  delete temp;

}

void SMSplot(bool doLimit,bool doLogy,bool doError) {

  bool writeRootFile=false;

  //  gStyle->SetPalette(1);
  //  TH2F *computeEfficiencyError(const TH2F *hEff, const TH1F *denom){

  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetPadLeftMargin(0.15);

  gStyle->SetPalette(1);

  gStyle->SetTitleYOffset(1.3);
  gStyle->SetTitleXOffset(1.1);
  //  gStyle->SetTitleZOffset(1);

  // SS-Ronny T1tttt
  TFile *_file_T1tttt_SS_2010_HTmet = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_highPt_400_120.root");
  TFile *_file_T1tttt_SS_2010_HT = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_highPt_400_50.root");
  TFile *_file_T1tttt_SS_2010_met = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_highPt_200_120.root");

  TFile *_file_T1tttt_SS_lowPt_HTmet = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_lowPt_400_120.root");
  TFile *_file_T1tttt_SS_lowPt_HT = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_lowPt_400_50.root");
  TFile *_file_T1tttt_SS_lowPt_met = TFile::Open("SSRonny/SMSResults_SMS_T1tttt_lowPt_200_120.root");

  TFile *_file_T1tttt_SS_Sanjay = TFile::Open("SSSanjay/RA5SSlimitMatched.root");

  TH2F * h_T1tttt_SS_2010_sig1_exp;
  if(doLimit) h_T1tttt_SS_2010_sig1_exp = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig1expectedUL");
  h_T1tttt_SS_2010_sig1_exp->SetName("h_T1tttt_SS_2010_sig1_exp");
  h_T1tttt_SS_2010_sig1_exp=restyleHisto(h_T1tttt_SS_2010_sig1_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig1_exp->SetName("h_T1tttt_SS_2010_sig1_limit");

  TH2F * h_T1tttt_SS_2010_sig2_exp;
  if(doLimit) h_T1tttt_SS_2010_sig2_exp = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig2expectedUL");
  h_T1tttt_SS_2010_sig2_exp->SetName("h_T1tttt_SS_2010_sig2_exp");
  h_T1tttt_SS_2010_sig2_exp=restyleHisto(h_T1tttt_SS_2010_sig2_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig2_exp->SetName("h_T1tttt_SS_2010_sig2_limit");

  TH2F * h_T1tttt_SS_2010_sig3_exp;
  if(doLimit) h_T1tttt_SS_2010_sig3_exp = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig3expectedUL");
  h_T1tttt_SS_2010_sig3_exp->SetName("h_T1tttt_SS_2010_sig3_exp");
  h_T1tttt_SS_2010_sig3_exp=restyleHisto(h_T1tttt_SS_2010_sig3_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig3_exp->SetName("h_T1tttt_SS_2010_sig3_limit");

  TH2F * h_T1tttt_SS_2010_sig4_exp;
  if(doLimit) h_T1tttt_SS_2010_sig4_exp = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig4expectedUL");
  h_T1tttt_SS_2010_sig4_exp->SetName("h_T1tttt_SS_2010_sig4_exp");
  h_T1tttt_SS_2010_sig4_exp=restyleHisto(h_T1tttt_SS_2010_sig4_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig4_exp->SetName("h_T1tttt_SS_2010_sig4_limit");


  TH2F * h_T1tttt_SS_2010_HTmet_exp;
  if(doLimit) h_T1tttt_SS_2010_HTmet_exp = (TH2F*) _file_T1tttt_SS_2010_HTmet->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_2010_HTmet_exp->SetName("h_T1tttt_SS_2010_HTmet_exp");
  h_T1tttt_SS_2010_HTmet_exp=restyleHisto(h_T1tttt_SS_2010_HTmet_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_HTmet_exp->SetName("h_T1tttt_SS_2010_HTmet_limit");

  TH2F * h_T1tttt_SS_2010_HT_exp;
  if(doLimit) h_T1tttt_SS_2010_HT_exp = (TH2F*) _file_T1tttt_SS_2010_HT->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_2010_HT_exp->SetName("h_T1tttt_SS_2010_HT_exp");
  h_T1tttt_SS_2010_HT_exp=restyleHisto(h_T1tttt_SS_2010_HT_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_HT_exp->SetName("h_T1tttt_SS_2010_HT_limit");

  TH2F * h_T1tttt_SS_2010_met_exp;
  if(doLimit) h_T1tttt_SS_2010_met_exp = (TH2F*) _file_T1tttt_SS_2010_met->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_2010_met_exp->SetName("h_T1tttt_SS_2010_met_exp");
  h_T1tttt_SS_2010_met_exp=restyleHisto(h_T1tttt_SS_2010_met_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_met_exp->SetName("h_T1tttt_SS_2010_met_limit");

  TH2F * h_T1tttt_SS_lowPt_HTmet_exp;
  if(doLimit) h_T1tttt_SS_lowPt_HTmet_exp = (TH2F*) _file_T1tttt_SS_lowPt_HTmet->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_HTmet_exp->SetName("h_T1tttt_SS_lowPt_HTmet_exp");
  h_T1tttt_SS_lowPt_HTmet_exp=restyleHisto(h_T1tttt_SS_lowPt_HTmet_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_HTmet_exp->SetName("h_T1tttt_SS_lowPt_HTmet_limit");

  TH2F * h_T1tttt_SS_lowPt_HT_exp;
  if(doLimit) h_T1tttt_SS_lowPt_HT_exp = (TH2F*) _file_T1tttt_SS_lowPt_HT->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_HT_exp->SetName("h_T1tttt_SS_lowPt_HT_exp");
  h_T1tttt_SS_lowPt_HT_exp=restyleHisto(h_T1tttt_SS_lowPt_HT_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_HT_exp->SetName("h_T1tttt_SS_lowPt_HT_limit");

  TH2F * h_T1tttt_SS_lowPt_met_exp;
  if(doLimit) h_T1tttt_SS_lowPt_met_exp = (TH2F*) _file_T1tttt_SS_lowPt_met->Get("Expected_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_met_exp->SetName("h_T1tttt_SS_lowPt_met_exp");
  h_T1tttt_SS_lowPt_met_exp=restyleHisto(h_T1tttt_SS_lowPt_met_exp,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_met_exp->SetName("h_T1tttt_SS_lowPt_met_limit");

  //////


  TH2F * h_T1tttt_SS_2010_sig1;
  if(doLimit) h_T1tttt_SS_2010_sig1 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig1UL");
  if(!doLimit)  h_T1tttt_SS_2010_sig1 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig1Eff");
  h_T1tttt_SS_2010_sig1->SetName("h_T1tttt_SS_2010_sig1");
  h_T1tttt_SS_2010_sig1=restyleHisto(h_T1tttt_SS_2010_sig1,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig1->SetName("h_T1tttt_SS_2010_sig1_limit");
  if(!doLimit) h_T1tttt_SS_2010_sig1->SetName("h_T1tttt_SS_2010_sig1_efficiency");

  TH2F * h_T1tttt_SS_2010_sig2;
  if(doLimit) h_T1tttt_SS_2010_sig2 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig2UL");
  if(!doLimit)  h_T1tttt_SS_2010_sig2 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig2Eff");
  h_T1tttt_SS_2010_sig2->SetName("h_T1tttt_SS_2010_sig2");
  h_T1tttt_SS_2010_sig2=restyleHisto(h_T1tttt_SS_2010_sig2,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig2->SetName("h_T1tttt_SS_2010_sig2_limit");
  if(!doLimit) h_T1tttt_SS_2010_sig2->SetName("h_T1tttt_SS_2010_sig2_efficiency");

  TH2F * h_T1tttt_SS_2010_sig3;
  if(doLimit) h_T1tttt_SS_2010_sig3 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig3UL");
  if(!doLimit)  h_T1tttt_SS_2010_sig3 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig3Eff");
  h_T1tttt_SS_2010_sig3->SetName("h_T1tttt_SS_2010_sig3");
  h_T1tttt_SS_2010_sig3=restyleHisto(h_T1tttt_SS_2010_sig3,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig3->SetName("h_T1tttt_SS_2010_sig3_limit");
  if(!doLimit) h_T1tttt_SS_2010_sig3->SetName("h_T1tttt_SS_2010_sig3_efficiency");

  TH2F * h_T1tttt_SS_2010_sig4;
  if(doLimit) h_T1tttt_SS_2010_sig4 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig4UL");
  if(!doLimit)  h_T1tttt_SS_2010_sig4 = (TH2F*) _file_T1tttt_SS_Sanjay->Get("h2sig4Eff");
  h_T1tttt_SS_2010_sig4->SetName("h_T1tttt_SS_2010_sig4");
  h_T1tttt_SS_2010_sig4=restyleHisto(h_T1tttt_SS_2010_sig4,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_sig4->SetName("h_T1tttt_SS_2010_sig4_limit");
  if(!doLimit) h_T1tttt_SS_2010_sig4->SetName("h_T1tttt_SS_2010_sig4_efficiency");


  TH2F * h_T1tttt_SS_2010_HTmet;
  if(doLimit) h_T1tttt_SS_2010_HTmet = (TH2F*) _file_T1tttt_SS_2010_HTmet->Get("Observed_limit_SMS_T1tttt");
  if(!doLimit)  h_T1tttt_SS_2010_HTmet= (TH2F*) _file_T1tttt_SS_2010_HTmet->Get("efficiency_SMS_T1tttt");
  h_T1tttt_SS_2010_HTmet->SetName("h_T1tttt_SS_2010_HTmet");
  h_T1tttt_SS_2010_HTmet=restyleHisto(h_T1tttt_SS_2010_HTmet,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_HTmet->SetName("h_T1tttt_SS_2010_HTmet_limit");
  if(!doLimit) h_T1tttt_SS_2010_HTmet->SetName("h_T1tttt_SS_2010_HTmet_efficiency");

  TH2F * h_T1tttt_SS_2010_HT;
  if(doLimit) h_T1tttt_SS_2010_HT = (TH2F*) _file_T1tttt_SS_2010_HT->Get("Observed_limit_SMS_T1tttt");
  if(!doLimit)  h_T1tttt_SS_2010_HT= (TH2F*) _file_T1tttt_SS_2010_HT->Get("efficiency_SMS_T1tttt");
  h_T1tttt_SS_2010_HT->SetName("h_T1tttt_SS_2010_HT");
  h_T1tttt_SS_2010_HT=restyleHisto(h_T1tttt_SS_2010_HT,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_HT->SetName("h_T1tttt_SS_2010_HT_limit");
  if(!doLimit) h_T1tttt_SS_2010_HT->SetName("h_T1tttt_SS_2010_HT_efficiency");

  TH2F * h_T1tttt_SS_2010_met;
  if(doLimit) h_T1tttt_SS_2010_met = (TH2F*) _file_T1tttt_SS_2010_met->Get("Observed_limit_SMS_T1tttt");
  if(!doLimit)  h_T1tttt_SS_2010_met= (TH2F*) _file_T1tttt_SS_2010_met->Get("efficiency_SMS_T1tttt");
  h_T1tttt_SS_2010_met->SetName("h_T1tttt_SS_2010_met");
  h_T1tttt_SS_2010_met=restyleHisto(h_T1tttt_SS_2010_met,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_2010_met->SetName("h_T1tttt_SS_2010_met_limit");
  if(!doLimit) h_T1tttt_SS_2010_met->SetName("h_T1tttt_SS_2010_met_efficiency");

  TH2F * h_T1tttt_SS_lowPt_HTmet;
  if(doLimit) h_T1tttt_SS_lowPt_HTmet = (TH2F*) _file_T1tttt_SS_lowPt_HTmet->Get("Observed_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_HTmet->SetName("h_T1tttt_SS_lowPt_HTmet");
  h_T1tttt_SS_lowPt_HTmet=restyleHisto(h_T1tttt_SS_lowPt_HTmet,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_HTmet->SetName("h_T1tttt_SS_lowPt_HTmet_limit");

  TH2F * h_T1tttt_SS_lowPt_HT;
  if(doLimit) h_T1tttt_SS_lowPt_HT = (TH2F*) _file_T1tttt_SS_lowPt_HT->Get("Observed_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_HT->SetName("h_T1tttt_SS_lowPt_HT");
  h_T1tttt_SS_lowPt_HT=restyleHisto(h_T1tttt_SS_lowPt_HT,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_HT->SetName("h_T1tttt_SS_lowPt_HT_limit");

  TH2F * h_T1tttt_SS_lowPt_met;
  if(doLimit) h_T1tttt_SS_lowPt_met = (TH2F*) _file_T1tttt_SS_lowPt_met->Get("Observed_limit_SMS_T1tttt");
  h_T1tttt_SS_lowPt_met->SetName("h_T1tttt_SS_lowPt_met");
  h_T1tttt_SS_lowPt_met=restyleHisto(h_T1tttt_SS_lowPt_met,doLimit,doError,"T1tttt");
  if(doLimit) h_T1tttt_SS_lowPt_met->SetName("h_T1tttt_SS_lowPt_met_limit");

  ////



  //$$$$$$$
  //$$$$$$$
  //$$$$$$$
  
  //  char * leg="GL GL -> 4jets + LSPs";
  //  char * leg="pp #rightarrow #tilde{G} #tilde{G} #rightarrow 4jets + 2 lept + 2 #nu + LSPs ";
  char * legTGQ="pp #rightarrow #tilde{g} #tilde{q}, #tilde{g} #rightarrow 2q + LSP; #tilde{q} #rightarrow q + LSP; m(#tilde{q})>>m(#tilde{g})";
  char * legT2="pp #rightarrow #tilde{q} #tilde{q}, #tilde{q} #rightarrow q + LSP; m(#tilde{g})>>m(#tilde{q})";
  char * legT1="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2q + LSP; m(#tilde{q})>>m(#tilde{g})";
  char * legT1bbbb="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2b + LSP; m(#tilde{b})>>m(#tilde{g})";
  char * legT1tttt="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2t + LSP; m(#tilde{t})>>m(#tilde{g})";

  char * legT1lnu="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2q + #chi^{#pm}, #chi^{#pm} #rightarrow l^{#pm}#nu + LSP;m(#tilde{q})>>m(#tilde{g})";  
  char * legT1lh="pp #rightarrow #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2q + LSP, #tilde{g} #rightarrow 2j + #chi^{0},#chi^{0} #rightarrow l^{+}l^{-}LSP;m(#tilde{q})>>m(#tilde{g})";  
  char * legT5zz="pp #rightarrow  #tilde{g} #tilde{g}, #tilde{g} #rightarrow 2q + #chi^{0}, #chi^{0} #rightarrow Z + LSP;m(#tilde{q})>>m(#tilde{g})";

  char * legT2tt="pp #rightarrow #tilde{t} #tilde{t}, #tilde{t} #rightarrow t + LSP; m(#tilde{g})>>m(#tilde{t})";

  //m(#chi^{0})=0.5*(m(#tilde{g})+m(#tilde{LSP})

  //  char * exp="Experimental uncertainty A #times #varepsilon";
  //  char * teo="Theoretical uncertainty A #times #varepsilon";
  char * selectionHTOS="OS: High H_{T} selection";
  char * selectionMETOS="OS: High #slash{E}_{T} selection";
  ////  char * selectionHTmetOS="OS: High #slash{E}_{T} High H_{T} selection";

  char * selectionANNOS="OS: ANN selection";

  char * selectionHTSS="SS: High H_{T} selection";
  char * selectionMHTSS="SS: High #slash{E}_{T} selection";
  char * selectionHTmetSS="SS: High #slash{E}_{T} High H_{T} selection";
  char * selectionHTmetSS_sb="SS (SB): High #slash{E}_{T} High H_{T} selection";

  char * selectionHTSS_lowPt="SS lowPt: High H_{T} selection";
  char * selectionmetSS_lowPt="SS lowPt : High #slash{E}_{T} selection";
  char * selectionHTmetSS_lowPt="SS lowPt : High #slash{E}_{T} High H_{T} selection";

  char * selectionSS_OR="SS : OR all selection";
  char * selectionOS_OR="OS : OR all selection";
  char * selectionRA2b_OR="#slash{E}_{T} + b : OR all selection";
  char * selectionRA2_OR="#slash{E}_{T} + jets: OR all selection";
  char * selectionWarren_OR="Z+MET : OR all selection";
  char * selectionJZB_OR="ZJB : OR all selection";

  char * selection_RA2_HT="H_{T} #geq 800 #slash{H}_{T} #geq 200";
  char * selection_RA2_MHT="H_{T} #geq 500 #slash{H}_{T} #geq 350";
  char * selection_RA2_lowHTMHT="H_{T} #geq 350 #slash{H}_{T} #geq 200";
  char * selection_RA2_highHTMHT="H_{T} #geq 800 #slash{H}_{T} #geq 500";

  char * selection_RA4_Tight="LS: tight";
  char * selection_RA4_Loose="LS: loose";

  char * selectionJZB50="LowJZB selection";  
  char * selectionJZB100="HighJZB selection";  

  char * selection_2b_Loose="#geq 2b, Loose";  
  char * selection_2b_Tight="#geq 2b, Tight";  

  char * selection_1b_Loose="#geq 1b, Loose";  
  char * selection_1b_Tight="#geq 1b, Tight";  

  char * selectionWarrenTight="High #slash{E}_{T} selection";
  char * selectionWarrenLoose="Low #slash{E}_{T} selection";

  char *selectionTed="#alpha_{T}";

  //$$$$$$
  //$$$$$$
  //$$$$$$




  printCanvas(h_T1tttt_SS_2010_HTmet,"T1tttt_SS_2010_HTmet", selectionHTmetSS, legT1tttt,"T1tttt",doLimit,doError,doLogy);
  //  printCanvas(h_T1tttt_SS_2010_HTmet,"T1tttt_SS_2010_HTmet", selectionHTmetSS_sb, legT1tttt,"T1tttt",doLimit,doError,doLogy);

  if(doLimit) {

  char * labelSS[6]={"1=2010HTmet","2=2010HT","3=2010met","4=lowPtHTmet","5=lowPtHT","6=lowPtmet"};

  vector<TH2F*> objvecExpT1tttt;
  objvecExpT1tttt.push_back(h_T1tttt_SS_2010_HTmet_exp);
  objvecExpT1tttt.push_back(h_T1tttt_SS_2010_HT_exp);
  objvecExpT1tttt.push_back(h_T1tttt_SS_2010_met_exp);
  objvecExpT1tttt.push_back(h_T1tttt_SS_lowPt_HTmet_exp);
  objvecExpT1tttt.push_back(h_T1tttt_SS_lowPt_HT_exp);
  objvecExpT1tttt.push_back(h_T1tttt_SS_lowPt_met_exp);

  vector<TH2F*> objvecObsT1tttt;
  objvecObsT1tttt.push_back(h_T1tttt_SS_2010_HTmet);
  objvecObsT1tttt.push_back(h_T1tttt_SS_2010_HT);
  objvecObsT1tttt.push_back(h_T1tttt_SS_2010_met);
  objvecObsT1tttt.push_back(h_T1tttt_SS_lowPt_HTmet);
  objvecObsT1tttt.push_back(h_T1tttt_SS_lowPt_HT);
  objvecObsT1tttt.push_back(h_T1tttt_SS_lowPt_met);


  vector<TH2F*> objvecExpSPT1tttt;
  objvecExpSPT1tttt.push_back(h_T1tttt_SS_2010_sig1_exp);
  objvecExpSPT1tttt.push_back(h_T1tttt_SS_2010_sig2_exp);
  objvecExpSPT1tttt.push_back(h_T1tttt_SS_2010_sig3_exp);
  objvecExpSPT1tttt.push_back(h_T1tttt_SS_2010_sig4_exp);

  vector<TH2F*> objvecObsSPT1tttt;
  objvecObsSPT1tttt.push_back(h_T1tttt_SS_2010_sig1);
  objvecObsSPT1tttt.push_back(h_T1tttt_SS_2010_sig2);
  objvecObsSPT1tttt.push_back(h_T1tttt_SS_2010_sig3);
  objvecObsSPT1tttt.push_back(h_T1tttt_SS_2010_sig4);

  //  vector<TH2F*> manySST1tttt=takeMinimum(objvecExpT1tttt,objvecObsT1tttt);

  vector<TH2F*> manySST1tttt=takeMinimum(objvecExpSPT1tttt,objvecObsSPT1tttt);

  TH2F* crossSecT1tttt=manySST1tttt.at(0);
  TH2F* analysisT1tttt=manySST1tttt.at(1);
  //  if(doLimit) crossSecT1tttt->Scale(1./1000);

  crossSecT1tttt=restyleHisto(crossSecT1tttt,doLimit,doError,"T1tttt");
  analysisT1tttt=restyleHisto(analysisT1tttt,doLimit,doError,"T1tttt");

  TCanvas *temp3t = new TCanvas();
  temp3t->cd();
  if(doLogy)  gPad->SetLogz(1);

  crossSecT1tttt->Draw("colz");
  labelling(crossSecT1tttt,selectionSS_OR,legT1tttt,"T1tttt",doLimit);

  if(doLimit) temp3t->SaveAs("RESULT/comparisonSSMatched_xSec_T1tttt.png");
  if(doLimit) temp3t->SaveAs("RESULT/comparisonSSMatched_xSec_T1tttt.eps");

  //  if(!doLimit) temp3t->SaveAs("RESULT/comparisonSS_eff_T1tttt.png");
  //  if(!doLimit) temp3t->SaveAs("RESULT/comparisonSS_eff_T1tttt.eps");

 
  if(doLimit) {
  // TFile fSS("SS_T1tttt_bestLimit.root","RECREATE");
  TFile fSS("SS_T1tttt_bestLimit.root","UPDATE");
  //  Stop10->Write();
  crossSecT1tttt->Write();
  fSS.Close();
  }

  delete temp3t;

  TCanvas *temp4t = new TCanvas();
  temp4t->cd();
  
  gStyle->SetNumberContours(6);
  analysisT1tttt->SetMinimum(1);
  analysisT1tttt->SetMaximum(6);
  analysisT1tttt->SetZTitle("");
  //  analysisT1tttt->GetZaxis()->SetLabelOffset(99);
  analysisT1tttt->Draw("colz text");
 
  TText t1;
  t1.SetTextAngle(0);
  t1.SetTextSize(0.04);
  t1.SetTextAlign(33);
  Float_t x1, y1;
  x1 = analysisT1tttt->GetXaxis()->GetBinCenter(analysisT1tttt->GetXaxis()->FindBin(700));

  for (int i=0;i<6;i++) {
    if(i==0) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(600));
    if(i==1) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(700));
    if(i==2) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(800));
    if(i==3) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(900));
    if(i==4) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(1000));
    if(i==5) y1 = analysisT1tttt->GetYaxis()->GetBinCenter(analysisT1tttt->GetYaxis()->FindBin(1100));
     t1.DrawText(x1,y1,labelSS[i]);
  }
  
  temp4t->SaveAs("RESULT/comparisonSS_analysis_T1tttt.png");
  temp4t->SaveAs("RESULT/comparisonSS_analysis_T1tttt.eps");

  delete temp4t;
  gStyle->SetNumberContours(20);

  ///////


  }


}


void SMS() {

  // for limit: 
  //  bool doLimit=true; bool doLogy=true; bool doError=false;
  SMSplot(true,true,false);

  //for Efficiency: 
  //  doLimit=false; doLogy=true; doError=false;
  SMSplot(false,false,false); 

  // for Efficiency Error: 
  //  doLimit=false; doLogy=false; doError=true;
  //////  SMSplot(false,false,true); 

}

