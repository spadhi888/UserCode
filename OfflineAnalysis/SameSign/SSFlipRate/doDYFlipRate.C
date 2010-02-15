#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include <iostream>
#include "DYFlipRate.h"
#endif __CINT__
using namespace std;


TString  bbFile("Zee_BBhistos.root");
TString beFile("Zee_BEhistos.root");
TString allFile("Zee_mcallhistos.root");


void TestPrediction(const char* barrelBarrelFile, const char* barrelEndcapFile, const char* allRegionsFile) {

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
//  TFile *fbb = TFile::Open(barrelBarrelFile.c_str(),"READ");
  TFile *fbb = TFile::Open(barrelBarrelFile,"READ");
  
  
  TH2F *bb_FLRPtvsEta = (TH2F*)fbb->Get("bb_FLRPtvsEta");
  TH2F *bb_mcFLRPtvsEta = (TH2F*)fbb->Get("bb_mcFLRPtvsEta");
  TH1F *bb_FLREta  = (TH1F*)fbb->Get("bb_FLREta");
  TH1F *bb_mcFLREta  = (TH1F*)fbb->Get("bb_mcFLREta");
  TH1F *bb_numPt = (TH1F*)fbb->Get("bb_diffPt"); //predicted wrong sign distribution
  TH1F *bb_denomPt = (TH1F*)fbb->Get("bb_denomPt"); //distribution of all bb electrons
  bb_FLRPtvsEta->SetDirectory(rootdir);
  bb_mcFLRPtvsEta->SetDirectory(rootdir);
  bb_FLREta->SetDirectory(rootdir);
  bb_mcFLREta->SetDirectory(rootdir);
  bb_numPt->SetDirectory(rootdir);
  bb_denomPt->SetDirectory(rootdir);
  fbb->Close();
  
  
//  TFile *fbe = TFile::Open(barrelEndcapFile.c_str(),"READ");
  TFile *fbe = TFile::Open(barrelEndcapFile,"READ");
  
  TH2F *be_FLRPtvsEta = (TH2F*)fbe->Get("be_FLRPtvsEta");
  TH2F *be_mcFLRPtvsEta = (TH2F*)fbe->Get("be_mcFLRPtvsEta");
  TH1F *be_FLREta  = (TH1F*)fbe->Get("be_FLREta");
  TH1F *be_mcFLREta  = (TH1F*)fbe->Get("be_mcFLREta");
  TH1F *be_FLRPt  = (TH1F*)fbe->Get("be_FLRPt");
  TH1F *be_mcFLRPt  = (TH1F*)fbe->Get("be_mcFLRPt");
  //get the ss and os Pt distributions
  TH1F *be_ssPt      = (TH1F*)fbe->Get("be_ssPt");
  TH1F *be_osPt      = (TH1F*)fbe->Get("be_osPt");
  
  //this is the MC Pt distribution of the incorrectly assigned electrons
  TH1F *be_mcPtW = (TH1F*)fbe->Get("be_mcPtW");
  be_mcPtW->SetName("be_mcPtW");

  //add the same sign and opposite sign distributions to get the total distribution
  TH1F *be_denomPt = (TH1F*)be_ssPt->Clone();
  be_denomPt->SetName("be_denomPt");
  be_denomPt->Add(be_osPt);
  be_denomPt->Sumw2();

  TH1F *be_numPt = (TH1F*)be_denomPt->Clone();
  be_numPt->SetName("be_numPt");
  be_numPt->Sumw2();

  for(int i = 0; i < be_denomPt->GetNbinsX()+1; i++) {
    //scale the total by the probability to get a fake
    be_numPt->SetBinContent(i, be_denomPt->GetBinContent(i)*be_FLRPt->GetBinContent(i));
    double denomPtErr2 = pow(be_denomPt->GetBinError(i),2);
    double FLRPtErr2 = pow(be_FLRPt->GetBinError(i),2);
    double err = sqrt(pow(be_denomPt->GetBinContent(i), 2)*FLRPtErr2 + pow(be_FLRPt->GetBinContent(i),2)*denomPtErr2);
    be_numPt->SetBinError(i, err);
  }

  //now multiply the Pt distribution with the FLRPt distribution
  be_FLRPtvsEta->SetDirectory(rootdir);
  be_mcFLRPtvsEta->SetDirectory(rootdir);
  be_FLREta->SetDirectory(rootdir);
  be_mcFLREta->SetDirectory(rootdir);
  be_FLRPt->SetDirectory(rootdir);
  be_mcFLRPt->SetDirectory(rootdir);
  be_ssPt->SetDirectory(rootdir);
  be_osPt->SetDirectory(rootdir);
  be_mcPtW->SetDirectory(rootdir);
  be_denomPt->SetDirectory(rootdir);
  be_numPt->SetDirectory(rootdir);
  
  
  fbe->Close();

//  TFile *fall = TFile::Open(allRegionsFile.c_str(), "READ");
  TFile *fall = TFile::Open(allRegionsFile, "READ");
  
  TH2F *all_mcFLRPtvsEta = (TH2F*)fall->Get("all_mcFLRPtvsEta");
  all_mcFLRPtvsEta->SetName("all_mcFLRPtvsEta");

  TH1F *all_mcFLREta     = (TH1F*)fall->Get("all_mcFLREta");
  all_mcFLREta->SetName("all_mcFLREta");



  all_mcFLRPtvsEta->SetDirectory(rootdir);
  all_mcFLREta->SetDirectory(rootdir);
  fall->Close();

  TH1F *all_FLREta = (TH1F*)be_FLREta->Clone();
  all_FLREta->SetName("all_FLREta");
  all_FLREta->Add(bb_FLREta);

  TH2F *all_FLRPtvsEta = (TH2F*)be_FLRPtvsEta->Clone();
  all_FLRPtvsEta->SetName("all_FLRPtvsEta");
  all_FLRPtvsEta->Add(bb_FLRPtvsEta);
  //get the Pt distributions ->all electrons
  TH1F *all_denomPt = (TH1F*)bb_denomPt->Clone();
  all_denomPt->SetName("all_denomPt");
  all_denomPt->Add(be_denomPt);
  
  //get the wrong sign electrons
  TH1F *all_numPt = (TH1F*)bb_numPt->Clone();
  all_numPt->SetName("all_numPt");
  all_numPt->Add(be_numPt);

  TH1F* all_FLRPt = (TH1F*)all_numPt->Clone();
  all_FLRPt->SetName("all_FLRPt");
  DivideHists(all_FLRPt, all_denomPt);

  all_mcFLREta->SetLineColor(2);
  all_mcFLREta->Draw();
  all_FLREta->Draw("same");

}

void DivideHists(TH1F *num, TH1F *denom) {

  for(int i = 0; i <= num->GetNbinsX(); i++) {
    double valnum   = num->GetBinContent(i);
    double valdenom = denom->GetBinContent(i);
    double sigmanum2   = pow(num->GetBinError(i),2);
    double sigmadenom2 = pow(denom->GetBinError(i),2);
    double err = 0.0;
    if(valdenom < 0.00000001) 
      num->SetBinContent(i,0);
    else {
      err = sqrt(sigmanum2/(valdenom*valdenom) + sigmadenom2*pow(valnum,2)/pow(valdenom,4));
      num->SetBinContent(i,valnum/valdenom);
      num->SetBinError(i,err);
    }
      
  }
}


void doAll(){

  //  gROOT->ProcessLine(Form(".x setup.C(%d)", 1));
  gROOT->ProcessLine(".L histtools.C++");
  gSystem->CompileMacro("DYFlipRate.C","++k", "libDYFlipRate");
  
  TChain *ch_dy = new TChain("Events");
  ch_dy->Add("/store/disk01/cms2-V01-03-01/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root");

  DYFlipRate *loop = new DYFlipRate();

  loop->ScanChainDY_BB(ch_dy, -1, "bb");
  saveHist(bbFile);
  deleteHistos();
  cout << "Done with Barrel" << endl;
 
  loop->ScanChainDY_BB(ch_dy, -1, "all");
  saveHist(allFile, "*mc*");
  deleteHistos();
  cout << "Done with MC" << endl;
  
  loop->ScanChainDY_BE(ch_dy, -1, beFile);
  saveHist(beFile);
  deleteHistos();

  cout << "Done with BE" << endl;
  
  delete loop;
  
  TestPrediction(bbFile, beFile, allFile);

}
