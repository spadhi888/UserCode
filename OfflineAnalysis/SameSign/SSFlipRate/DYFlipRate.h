// -*- C++ -*-
#ifndef DYFlipRate_H
#define DYFlipRate_H
#include "CMS2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

class DYFlipRate {

public: 
 
  int ScanChainDY_BB( TChain* chain, int nEventsMax = -1, std::string region="bb", 
		      bool useSingleElectronWeights = false);
  
  int ScanChainDY_BE( TChain* chain, int nEvents = -1, std::string BBfname="");

  int getClosestStatus3ParticleIdx(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4);
  int bookHistos(std::string prefix);
  void SubtractHistsCorrectly(TH1F* ss, TH1F *os);
  void SubtractHistsCorrectly(TH2F* ss, TH2F *os);
  void DivideHists(TH1F *num, TH1F *denom);
  void DivideHists(TH2F *num, TH2F *denom);
  void DivideHistsBinomial(TH1F *num, TH1F *denom);
  void DivideHistsBinomial(TH2F *num, TH2F *denom);

  // BB histos
  
  TH1F *h_ssPt[3]; //0 is Plus, 1 is Minus, 3 is the sum
  TH1F *h_osPt[3];
  
  TH1F *h_ssEta[3];
  TH1F *h_osEta[3];
  
  //2d distributions
  TH2F *h_ssPtvsEta[3];
  TH2F *h_osPtvsEta[3];
  
  
  TH1F *h_mcssPt[3];
  TH1F *h_mcosPt[3];
  
  TH1F *h_mcssEta[3];
  TH1F *h_mcosEta[3];
  
  //2d distributions
  TH2F *h_mcssPtvsEta[3];
  TH2F *h_mcosPtvsEta[3];
  
  //seg = single electron gun applied weight
  TH1F *h_allPt_seg;
  TH1F *h_allEta_seg;
  
  TH1F *h_allFLRPt_seg;
  TH1F *h_allFLREta_seg;
  TH2F *h_allPtvsEta_seg;
  
  TH1F *h_allPt;
  TH1F *h_allEta;
  TH2F *h_allPtvsEta;
  
  //predicted wrong sign distributions
  TH1F *h_diffPt[3];
  TH1F *h_diffEta[3];
  TH2F *h_diffPtvsEta[3];
  
  TH1F *h_FLRPt;
  TH1F *h_FLREta;
  TH2F *h_FLRPtvsEta;
  
  TH1F *h_mcFLRPt;
  TH1F *h_mcFLREta;
  TH2F *h_mcFLRPtvsEta;
  
  
  // BE histos
  TH1F *h_ssPtv;
  TH1F *h_ssEtav;
  TH2F *h_ssPtvsEtav;
  
  TH1F *h_osPtv;
  TH1F *h_osEtav;
  TH2F *h_osPtvsEtav;
  
  //right charge
  TH1F  *h_mcPtC;
  TH1F  *h_mcEtaC;
  TH2F  *h_mcPtvsEtaC;
  
  //wrong charge
  TH1F  *h_mcPtW;
  TH1F  *h_mcEtaW;
  TH2F  *h_mcPtvsEtaW;
  
  //Fake Rates
//  TH1F *h_FLRPt;
//  TH1F *h_FLREta;
//  TH2F *h_FLRPtvsEta;
  
  //Monte Carlo Fake Rates
//  TH1F *h_mcFLRPt;
//  TH1F *h_mcFLREta;
//  TH2F *h_mcFLRPtvsEta;
  
  
  //temp histos
  TH1F *h_APt;
  TH1F *h_BPt;
  
  TH1F *h_AEta;
  TH1F *h_BEta;
  
  TH2F *h_APtvsEta;
  TH2F *h_BPtvsEta;

};
#endif

