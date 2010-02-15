#include <iostream>
#include <vector>
#include <map>

#include "TChain.h"
#include "TFile.h"
#include "TChainElement.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "CMS2.h"

#include "selections.cc"
#include "utilities.cc"

#include "TH1F.h"
#include "TH2F.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include <algorithm>
#include "TRandom2.h"
#include <fstream>
#include "TChain.h"

#include "DYFlipRate.h"
#include "fliprate_egun.cc"

CMS2 cms2;
using namespace tas;


//////////////////////////////////////////////////////////////////////////////

int DYFlipRate::bookHistos(std::string prefix) {

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  std::vector<string> ch;
  ch.push_back("Plus");
  ch.push_back("Minus");
  ch.push_back("");

  Double_t eta[8]  = {0.0, 0.5,1.0,1.28, 1.56, 1.84, 2.12, 2.4};
  int numberBins_eta = 7;
  Double_t pt[5] = {10.0, 30.0, 40.0, 50.0, 200.0};
  int numberBins_pt = 4;

  h_allPt_seg = new TH1F(Form("%s_allPt_seg",prefix.c_str()), "Pt distribution of all electrons, OS and SS with the electron FLR weight from the electron gun", numberBins_pt, pt);

  h_allEta_seg = new TH1F(Form("%s_allEta_seg",prefix.c_str()), "Eta distribution of all electrons, OS and SS with the electron FLR weight from the electron gun", numberBins_eta, eta);

  h_allPtvsEta_seg = new TH2F(Form("%s_allPtvsEta_seg",prefix.c_str()), "Pt vs Eta distribution for all electrons with the electron FLR from the electron gun", numberBins_eta, eta, numberBins_pt, pt);


  h_allPt = new TH1F(Form("%s_allPt",prefix.c_str()), "Pt distribution of all electrons, OS and SS", numberBins_pt, pt);

  h_allEta = new TH1F(Form("%s_allEta",prefix.c_str()), "Eta distribution of all electrons, OS and SS", numberBins_eta, eta);

  h_allPtvsEta = new TH2F(Form("%s_allPtvsEta",prefix.c_str()), "Pt vs Eta distribution for all electrons", numberBins_eta, eta, numberBins_pt, pt);

  //flip rate
  h_allFLRPt_seg = new TH1F(Form("%s_allFLRPt_seg",prefix.c_str()), "FLR distribution of all electrons vs Pt, OS and SS with the electron FLR weight from the single electron gun", numberBins_pt, pt);

  h_allFLREta_seg = new TH1F(Form("%s_allFLREta_seg",prefix.c_str()), "FLR distribution of all electrons, OS and SS with the electron FLR weight from the single electron gun", numberBins_eta, eta);


  h_allPt_seg->SetDirectory(rootdir);
  h_allEta_seg->SetDirectory(rootdir);
  h_allPtvsEta_seg->SetDirectory(rootdir);
  
  h_allPt_seg->Sumw2();
  h_allEta_seg->Sumw2();
  h_allPtvsEta_seg->Sumw2();
  
  h_allPt->SetDirectory(rootdir);
  h_allEta->SetDirectory(rootdir);
  h_allPtvsEta->SetDirectory(rootdir);
  
  h_allPt->Sumw2();
  h_allEta->Sumw2();
  h_allPtvsEta->Sumw2();
  
  h_allFLRPt_seg->SetDirectory(rootdir);
  h_allFLREta_seg->SetDirectory(rootdir);

  h_allFLRPt_seg->Sumw2();
  h_allFLREta_seg->Sumw2();

  // BE

  h_ssPtv      = new TH1F((prefix+"_ssPtv").c_str(), "Pt distribution, SS", numberBins_pt, pt);
  h_ssEtav     = new TH1F((prefix+"_ssEtav").c_str(), "SS", numberBins_eta, eta);
  
  h_osPtv      = new TH1F((prefix+"_osPtv").c_str(), "Pt distribution, SS", numberBins_pt, pt);
  h_osEtav     = new TH1F((prefix+"_osEtav").c_str(), "SS", numberBins_eta, eta);

  
  h_ssPtvsEtav = new TH2F((prefix+"_ssPtvsEtav").c_str(), "Pt vs Eta distribution, SS",
			  numberBins_eta, eta, numberBins_pt, pt);
  h_osPtvsEtav = new TH2F((prefix+"_osPtvsEtav").c_str(), "Pt vs Eta distribution, OS",
			  numberBins_eta, eta, numberBins_pt, pt);
  

  //Monte Carlo Histograms
  h_mcPtC      = new TH1F((prefix+"_mcPtC").c_str(),  "Correctly assigned charge Monte Carlo Pt distribution", numberBins_pt, pt);
  h_mcEtaC    = new TH1F((prefix+"_mcEtaC").c_str(), "Correctly assigned charge Monte Carlo Eta distribution", numberBins_eta, eta);
  h_mcPtvsEtaC = new TH2F((prefix+"_mcPtvsEtaC").c_str(), "Correctly assigned charge Monte Carlo Pt vs Eta distribution",
			  numberBins_eta, eta, numberBins_pt, pt);
  
  h_mcPtW      = new TH1F((prefix+"_mcPtW").c_str(),  "Wrongly assigned charge Monte Carlo Pt distribution", numberBins_pt, pt);
  h_mcEtaW    = new TH1F((prefix+"_mcEtaW").c_str(), "Wrongly assigned charge Monte Carlo Eta distribution", numberBins_eta, eta);
  h_mcPtvsEtaW = new TH2F((prefix+"_mcPtvsEtaW").c_str(), "Wrongly assigned charge Monte Carlo Pt vs Eta distribution",
                          numberBins_eta, eta, numberBins_pt, pt);



  h_APt = (TH1F*)h_ssPtv->Clone();
  h_APt->SetName("h_APt");
  h_BPt = (TH1F*)h_ssPtv->Clone();
  h_BPt->SetName("h_BPt");


  h_AEta = (TH1F*)h_ssEtav->Clone();
  h_AEta->SetName("h_AEta");
  h_BEta = (TH1F*)h_ssEtav->Clone();
  h_BEta->SetName("h_BEta");


  h_APtvsEta = (TH2F*)h_ssPtvsEtav->Clone();
  h_APtvsEta->SetName("h_APtvsEta");

  h_BPtvsEta = (TH2F*)h_ssPtvsEtav->Clone();
  h_BPtvsEta->SetName("h_APtvsEta");

  h_ssPtv->Sumw2();
  h_osPtv->Sumw2();

  h_ssEtav->Sumw2();
  h_osEtav->Sumw2();

  h_ssPtvsEtav->Sumw2();
  h_osPtvsEtav->Sumw2();

  h_ssPtv->SetDirectory(rootdir);
  h_osPtv->SetDirectory(rootdir);

  h_ssEtav->SetDirectory(rootdir);
  h_osEtav->SetDirectory(rootdir);

  h_ssPtvsEtav->SetDirectory(rootdir);
  h_osPtvsEtav->SetDirectory(rootdir);

  h_mcPtC->Sumw2();
  h_mcEtaC->Sumw2();
  h_mcPtvsEtaC->Sumw2();

  h_mcPtW->Sumw2();
  h_mcEtaW->Sumw2();
  h_mcPtvsEtaW->Sumw2();

  h_mcPtC->SetDirectory(rootdir);
  h_mcEtaC->SetDirectory(rootdir);
  h_mcPtvsEtaC->SetDirectory(rootdir);

  h_mcPtW->SetDirectory(rootdir);
  h_mcEtaW->SetDirectory(rootdir);
  h_mcPtvsEtaW->SetDirectory(rootdir);

  h_APt->Sumw2();
  h_BPt->Sumw2();

  h_AEta->Sumw2();
  h_BEta->Sumw2();

  h_APtvsEta->Sumw2();
  h_BPtvsEta->Sumw2();


  // BB rest

  for(unsigned int i = 0; i < 3; i++) {

    //Double_t pt[5] = {10.,40.,70.,100.,200.}; //works ok ->first try
    //this is to try and see where the stats are. The pt bins will be 20 bins, 0 to 200
    //Double_t eta[10]  = {0.0,0.25, 0.5,0.75,1.0,1.28, 1.56, 1.84, 2.12, 2.4};

    //combine some bins to reduce the existance of 0 bins

    h_ssPt[i] = new TH1F(string(prefix+"_ss"+ch[i]+"Pt").c_str(),
                         string("Pt distribution of "+ch[i]+" electrons from same sign events").c_str(),
                         numberBins_pt, pt);
    h_osPt[i] = new TH1F(string(prefix+"_os"+ch[i]+"Pt").c_str(),
                         string("Pt distribution of "+ch[i]+" electrons from OS events").c_str(),
                         numberBins_pt, pt);


    h_mcssPt[i] = new TH1F(string(prefix+"_mcss"+ch[i]+"Pt").c_str(),
                           string("Pt distribution (reco) of "+ch[i]+" electrons matched to MC els with wrong sign").c_str(),
                           numberBins_pt, pt);
    h_mcosPt[i] = new TH1F(string(prefix+"_mcos"+ch[i]+"Pt").c_str(),
                           string("Pt distribution (reco) of "+ch[i]+" electrons matched to MC els with the right  sign").c_str(),
                           numberBins_pt, pt);


    h_ssEta[i] = new TH1F(string(prefix+"_ss"+ch[i]+"Eta").c_str(),
                          string("Eta distribution of "+ch[i]+" electrons from SS events").c_str(),
                          numberBins_eta, eta);
    h_osEta[i] = new TH1F(string(prefix+"_os"+ch[i]+"Eta").c_str(),
                          string("Eta distribution of "+ch[i]+" electrons from OS events").c_str(),
                          numberBins_eta, eta);
    h_mcssEta[i] = new TH1F(string(prefix+"_mcss"+ch[i]+"Eta").c_str(),
                            string("Eta dist (reco) of "+ch[i]+" electrons matched to MC els with wrong sign").c_str(),
                            numberBins_eta, eta);
    h_mcosEta[i] = new TH1F(string(prefix+"_mcos"+ch[i]+"Eta").c_str(),
                            string("Eta dist (reco) of "+ch[i]+" electrons matched to MC els with the right sign").c_str(),
                            numberBins_eta, eta);


    //TH2Fs
    h_ssPtvsEta[i] = new TH2F(string(prefix+"_ss"+ch[i]+"PtvsEta").c_str(),
                              string("Pt vs #eta distribution of "+ch[i]+" electrons from same sign events").c_str(),
                              numberBins_eta, eta, numberBins_pt, pt);
    h_osPtvsEta[i] = new TH2F(string(prefix+"_os"+ch[i]+"PtvsEta").c_str(),
                              string("Pt vs #eta distribution of "+ch[i]+" electrons from OS events").c_str(),
                              numberBins_eta, eta, numberBins_pt, pt);
    h_mcssPtvsEta[i] = new TH2F(string(prefix+"_mcss"+ch[i]+"PtvsEta").c_str(),
                                string("Pt vs #eta distriution (reco) of "+ch[i]+" electrons matched to MC els with wrong sign").c_str(),
                                numberBins_eta, eta, numberBins_pt, pt);
    h_mcosPtvsEta[i] = new TH2F(string(prefix+"_mcos"+ch[i]+"PtvsEta").c_str(),
                                string("Pt vs #eta dist (reco) of "+ch[i]+" electrons matched to MC els with the right  sign").c_str(),
                                numberBins_eta, eta, numberBins_pt, pt);


    h_ssPt[i]->SetDirectory(rootdir);
    h_osPt[i]->SetDirectory(rootdir);

    h_mcssPt[i]->SetDirectory(rootdir);
    h_mcssEta[i]->SetDirectory(rootdir);

    h_mcosPt[i]->SetDirectory(rootdir);
    h_mcosEta[i]->SetDirectory(rootdir);

    h_ssEta[i]->SetDirectory(rootdir);
    h_osEta[i]->SetDirectory(rootdir);

    h_ssPtvsEta[i]->SetDirectory(rootdir);
    h_osPtvsEta[i]->SetDirectory(rootdir);


    h_mcssPtvsEta[i]->SetDirectory(rootdir);
    h_mcosPtvsEta[i]->SetDirectory(rootdir);


    h_ssPt[i]->Sumw2();
    h_osPt[i]->Sumw2();
    h_mcssPt[i]->Sumw2();
    h_mcssEta[i]->Sumw2();

    h_ssEta[i]->Sumw2();
    h_osEta[i]->Sumw2();

    h_ssPtvsEta[i]->Sumw2();
    h_osPtvsEta[i]->Sumw2();
    h_mcssPtvsEta[i]->Sumw2();
    h_mcosPtvsEta[i]->Sumw2();

  }
  return 0;
}
//////////////////////////////////////////////////////////////////////////////

int DYFlipRate::ScanChainDY_BB( TChain* chain, int nEventsMax, std::string region, bool useSingleElectronWeights) {

  if(region == "") {
    cout << "Need to specify an eta region" << endl;
    return 0;
  }
  
  if(region != "all" && region != "be" && region != "bb" && region != "ee") {
    cout << region << " is not a valid region" << endl;
    return 0;
  }

  TObjArray *listOfFiles = chain->GetListOfFiles();
  
  unsigned int nEventsChain=0;
  unsigned int nEventsTotal = 0;
  bookHistos(region);
  int ll_matchedToPhotonOS = 0;
  int ll_matchedToPhotonSS = 0;
  
  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    nEventsChain += nEvents;
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
      
      if(nEventsTotal%50000 == 0)
        std::cout << "Event processed = " << nEventsTotal << std::endl;
      
      if(!cms2.passHLTTrigger("HLT_Ele15_SW_L1R"))
        continue;
      
      for(unsigned int i = 0; i < cms2.hyp_p4().size(); i++) {
	
        if(cms2.hyp_type()[i] != 3 )
          continue;
	
        if(cms2.hyp_p4()[i].M() < 76. || cms2.hyp_p4()[i].M() > 106.)
          continue;
	
        float ltPt  = TMath::Min(cms2.hyp_lt_p4()[i].Pt(), 199.9);
        float ltEta = TMath::Min(fabs(cms2.hyp_lt_p4()[i].Eta()), 2.399);
	
        float llPt  = TMath::Min(cms2.hyp_ll_p4()[i].Pt(), 199.9);
        float llEta = TMath::Min(fabs(cms2.hyp_ll_p4()[i].Eta()), 2.399);
	
        unsigned int ltidx = cms2.hyp_lt_index()[i];
        unsigned int llidx = cms2.hyp_ll_index()[i];
	
        if(region == "bb") {
          //do only barrel-barrel events
          if(ltEta > 1.0 || llEta > 1.0)
            continue;
        }

        if(region == "ee" ) {
	  //do only endcap-endcap events
	  if(ltEta < 1.0 || llEta < 1.0)
	    continue;
        }
        if(region == "be") {
	  //ask that it is an barrel-endcap event
          if(ltEta > 1.0 && llEta > 1.0)
            continue;
          if(ltEta < 1.0 && llEta < 1.0)
            continue;
        }
	
	
        if(!GoodSusyLeptonID(cms2.hyp_lt_id()[i], cms2.hyp_lt_index()[i]))
          continue;
        if(!GoodSusyLeptonID(cms2.hyp_ll_id()[i], cms2.hyp_ll_index()[i]))
          continue;
        if(!GoodSusyLeptonWithIsolation(cms2.hyp_lt_id()[i], cms2.hyp_lt_index()[i]))
          continue;
        if(!GoodSusyLeptonWithIsolation(cms2.hyp_ll_id()[i], cms2.hyp_ll_index()[i]))
          continue;
	
	
        if(conversionElectron(ltidx))
          continue;
        if(conversionElectron(llidx))
          continue;
        if(isChargeFlip(ltidx))
          continue;
        if(isChargeFlip(llidx))
          continue;
	
	
        float ltweight = 1.0;
        float llweight = 1.0;
        if(useSingleElectronWeights) {
          ltweight = getSingleEleFlipRate(ltPt, ltEta);
          llweight = getSingleEleFlipRate(llPt, llEta);
        }
	
        h_allPt_seg->Fill(ltPt, ltweight);
        h_allEta_seg->Fill(ltEta, ltweight);
        h_allPtvsEta_seg->Fill(ltEta, ltPt, ltweight);
	
        h_allPt_seg->Fill(llPt, llweight);
        h_allEta_seg->Fill(llEta, llweight);
        h_allPtvsEta_seg->Fill(llEta, llPt, llweight);
	
        h_allPt->Fill(ltPt);
        h_allEta->Fill(ltEta);
        h_allPtvsEta->Fill(ltEta, ltPt);

        h_allPt->Fill(llPt);
        h_allEta->Fill(llEta);
        h_allPtvsEta->Fill(llEta, llPt);


        float weightlt_Pt      = 1.;
        float weightlt_Eta     = 1.;
        float weightlt_PtvsEta = 1.;

        float weightll_Pt      = 1.;
        float weightll_Eta     = 1.;
        float weightll_PtvsEta = 1.;

        //OS case
        if(cms2.hyp_lt_id()[i]*cms2.hyp_ll_id()[i] < 0) {
          if(cms2.hyp_lt_id()[i] > 0 ) {
            //+
            h_osPt[0] ->Fill(ltPt, weightlt_Pt);
            h_osEta[0]->Fill(ltEta, weightlt_Eta);
            h_osPtvsEta[0]->Fill(ltEta, ltPt, weightlt_PtvsEta);
            //-
            h_osPt[1] ->Fill(llPt, weightll_Pt);
            h_osEta[1]->Fill(llEta, weightll_Eta);
            h_osPtvsEta[1]->Fill(llEta, llPt, weightll_PtvsEta);
          }
          if(cms2.hyp_lt_id()[i] < 0 ) {//tight lepton is -
            //+
            h_osPt[0] ->Fill(llPt, weightll_Pt);
            h_osEta[0]->Fill(llEta, weightll_Eta);
            h_osPtvsEta[0]->Fill(llEta, llPt, weightll_PtvsEta);
            //-
            h_osPt[1] ->Fill(ltPt, weightlt_Pt);
            h_osEta[1]->Fill(ltEta, weightlt_Eta);
            h_osPtvsEta[1]->Fill(ltEta, ltPt, weightlt_PtvsEta);
          }// if tight lepton is -
          h_osPt[2]->Fill(ltPt, weightlt_Pt);
          h_osPt[2]->Fill(llPt, weightll_Pt);
          h_osEta[2]->Fill(ltEta, weightlt_Eta);
          h_osEta[2]->Fill(llEta, weightll_Eta);
          h_osPtvsEta[2]->Fill(ltEta, ltPt, weightlt_PtvsEta);
          h_osPtvsEta[2]->Fill(llEta, llPt, weightll_PtvsEta);
	  
        } else {//SS sign
          if(cms2.hyp_lt_id()[i] > 0) {
            //++
            h_ssPt[0] ->Fill(ltPt, weightlt_Pt);
            h_ssEta[0]->Fill(ltEta, weightlt_Eta);
            h_ssPtvsEta[0]->Fill(ltEta, ltPt, weightlt_PtvsEta);
            h_ssPt[0] ->Fill(llPt, weightll_Pt);
            h_ssEta[0]->Fill(llEta, weightll_Eta);
            h_ssPtvsEta[0]->Fill(llEta, llPt, weightll_PtvsEta);
          } else {
            //--
            h_ssPt[1] ->Fill(ltPt, weightlt_Pt);
            h_ssEta[1]->Fill(ltEta, weightlt_Eta);
            h_ssPtvsEta[1]->Fill(ltEta, ltPt, weightlt_PtvsEta);
            h_ssPt[1] ->Fill(llPt,weightll_Pt);
            h_ssEta[1]->Fill(llEta, weightll_Eta);
            h_ssPtvsEta[1]->Fill(llEta, llPt, weightll_PtvsEta);
          }

          //cout << cms2.hyp_p4()[i].M() << endl;
          //dumpDocLines

          h_ssPt[2]->Fill(ltPt, weightlt_Pt);
          h_ssPt[2]->Fill(llPt, weightll_Pt);
          h_ssEta[2]->Fill(ltEta, weightlt_Eta);
          h_ssEta[2]->Fill(llEta,weightll_Eta);
          h_ssPtvsEta[2]->Fill(ltEta, ltPt, weightlt_PtvsEta);
          h_ssPtvsEta[2]->Fill(llEta, llPt, weightll_PtvsEta);
        }//SS
        //////////////////////////////////////////////////////////////////////////////
        //now do the MC matching
        //////////////////////////////////////////////////////////////////////////////
        int lt_genid  = cms2.els_mc_id()[cms2.hyp_lt_index()[i]];
        int ll_genid  = cms2.els_mc_id()[cms2.hyp_ll_index()[i]];
        int lt_mom_genid = cms2.els_mc_motherid()[cms2.hyp_lt_index()[i]];
        int ll_mom_genid = cms2.els_mc_motherid()[cms2.hyp_ll_index()[i]];


        int lt_genid3 = cms2.els_mc3_id()[cms2.hyp_lt_index()[i]];
        int ll_genid3 = cms2.els_mc3_id()[cms2.hyp_ll_index()[i]];
        int lt_mom_genid3 = cms2.els_mc3_motherid()[cms2.hyp_lt_index()[i]];
        int ll_mom_genid3 = cms2.els_mc3_motherid()[cms2.hyp_ll_index()[i]];


        if(abs(lt_genid) == 11 && lt_mom_genid == 23) { //if matched to electron from Z
          if(cms2.hyp_lt_id()[i]*lt_genid < 0) { //is opposite sign from the mc particle, so the wrong sign (ss)
            if(cms2.hyp_lt_id()[i] > 0) { //lt Plus
              h_mcssPt[0]  ->Fill(ltPt);
              h_mcssEta[0] ->Fill(ltEta);
              h_mcssPtvsEta[0]->Fill(ltEta, ltPt);
            }
            if(cms2.hyp_lt_id()[i] < 0) {//lt Minus
              h_mcssPt[1]  ->Fill(ltPt);
              h_mcssEta[1] ->Fill(ltEta);
              h_mcssPtvsEta[1]->Fill(ltEta, ltPt);
            }
            h_mcssPt[2]  ->Fill(ltPt);
            h_mcssEta[2] ->Fill(ltEta);
            h_mcssPtvsEta[2]->Fill(ltEta, ltPt);
          } else {
            if(cms2.hyp_lt_id()[i] > 0) { //is the same sign as the mc particle
              h_mcosPt[0]  ->Fill(ltPt);
              h_mcosEta[0] ->Fill(ltEta);
              h_mcosPtvsEta[0]->Fill(ltEta, ltPt);
            }
            if(cms2.hyp_lt_id()[i] < 0) {
              h_mcosPt[1]  ->Fill(ltPt);
              h_mcosEta[1] ->Fill(ltEta);
              h_mcosPtvsEta[1]->Fill(ltEta, ltPt);
            }
            h_mcosPt[2]  ->Fill(ltPt);
            h_mcosEta[2] ->Fill(ltEta);
            h_mcosPtvsEta[2]->Fill(ltEta, ltPt);
          }// } else {
        }//if matched to electron from Z


        if(abs(lt_genid) == 22 && abs(lt_mom_genid) == 11) {
          //cout << "tight from gamma" << endl; //lots here!!!!->look at mc3
          if(abs(lt_genid3) == 11 && lt_mom_genid3 == 23) {
            if(cms2.hyp_lt_id()[i]*lt_genid3 < 0) { //the gen particle and reco are the opposite sign, (wrong sign, so ss)
              if(cms2.hyp_lt_id()[i] > 0) {//lt Plus
                h_mcssPt[0]  ->Fill(ltPt);
                h_mcssEta[0] ->Fill(ltEta);
                h_mcssPtvsEta[0]->Fill(ltEta, ltPt);
              }
              if(cms2.hyp_lt_id()[i] < 0) {//lt Minus
                h_mcssPt[1] ->Fill(ltPt);
                h_mcssEta[1]->Fill(ltEta);
                h_mcssPtvsEta[1]->Fill(ltEta, ltPt);
              }//ltMinus
              h_mcssPt[2] ->Fill(ltPt);
              h_mcssEta[2]->Fill(ltEta);
              h_mcssPtvsEta[1]->Fill(ltEta, ltPt);
            } else {//opposite sign
              if(cms2.hyp_lt_id()[i] > 0) {//lt Plus
                h_mcosPt[0]  ->Fill(ltPt);
                h_mcosEta[0] ->Fill(ltEta);
                h_mcosPtvsEta[0]->Fill(ltEta, ltPt);
              }
              if(cms2.hyp_lt_id()[i] < 0) { //lt Minus
                h_mcosPt[1]  ->Fill(ltPt);
                h_mcosEta[1] ->Fill(ltEta);
                h_mcosPtvsEta[1]->Fill(ltEta, ltPt);
              }
              h_mcosPt[2]  ->Fill(ltPt);
              h_mcosEta[2] ->Fill(ltEta);
              h_mcosPtvsEta[2]->Fill(ltEta, ltPt);
            }//} else {//opposite sign
          }
        }



        if(abs(ll_genid) == 11 && ll_mom_genid == 23) { //if matched to electron from Z
          if(cms2.hyp_ll_id()[i]*ll_genid < 0) { //is opposite sign from the mc particle, so ss
            if(cms2.hyp_ll_id()[i] > 0) {//ll Plus
              h_mcssPt[0] ->Fill(llPt);
              h_mcssEta[0]->Fill(llEta);
              h_mcssPtvsEta[0]->Fill(llEta, llPt);
            }
            if(cms2.hyp_ll_id()[i] < 0) {//ll Minus
              h_mcssPt[1]  ->Fill(llPt);
              h_mcssEta[1] ->Fill(llEta);
              h_mcssPtvsEta[1]->Fill(llEta, llPt);
            }
            h_mcssPt[2]  ->Fill(llPt);
            h_mcssEta[2] ->Fill(llEta);
            h_mcssPtvsEta[2]->Fill(llEta, llPt);
          } else {// is same sign as the mc particle
            if(cms2.hyp_ll_id()[i] > 0) {//ll Plus
              h_mcosPt[0] ->Fill(llPt);
              h_mcosEta[0]->Fill(llEta);
              h_mcosPtvsEta[0]->Fill(llEta, llPt);
            }
            if(cms2.hyp_ll_id()[i] < 0) {//ll Minus
              h_mcosPt[1]  ->Fill(llPt);
              h_mcosEta[1] ->Fill(llEta);
              h_mcosPtvsEta[1]->Fill(llEta, llPt);
            }
            h_mcosPt[2]  ->Fill(llPt);
            h_mcosEta[2] ->Fill(llEta);
            h_mcosPtvsEta[2]->Fill(llEta, llPt);
          }//} else {
        }

        if(abs(ll_genid) == 22 && abs(ll_mom_genid) == 11) {
          //cout << "loose from gamma" << endl;
          if(abs(ll_genid3) == 11 && ll_mom_genid3 == 23){
            if(cms2.hyp_ll_id()[i]*ll_genid3 < 0) {
              ll_matchedToPhotonSS++;

              if(cms2.hyp_ll_id()[i] > 0) {//ll Plus
                h_mcssPt[0] ->Fill(llPt);
                h_mcssEta[0]->Fill(llEta);
                h_mcssPtvsEta[0]->Fill(llEta, llPt);
              }
              if(cms2.hyp_ll_id()[i] < 0) {//ll Minus
                h_mcssPt[1]  ->Fill(llPt);
                h_mcssEta[1] ->Fill(llEta);
                h_mcssPtvsEta[1]->Fill(llEta, llPt);
              }
              h_mcssPt[2]  ->Fill(llPt);
              h_mcssEta[2] ->Fill(llEta);
              h_mcssPtvsEta[2]->Fill(llEta, llPt);
            } else {
              ll_matchedToPhotonOS++;
              if(cms2.hyp_ll_id()[i] > 0) {//ll Plus
                h_mcosPt[0] ->Fill(llPt);
                h_mcosEta[0]->Fill(llEta);
                h_mcosPtvsEta[0]->Fill(llEta, llPt);
              }
              if(cms2.hyp_ll_id()[i] < 0) {//ll Minus
                h_mcosPt[1]  ->Fill(llPt);
                h_mcosEta[1] ->Fill(llEta);
                h_mcosPtvsEta[1]->Fill(llEta, llPt);
              }
              h_mcosPt[2]  ->Fill(llPt);
              h_mcosEta[2] ->Fill(llEta);
              h_mcosPtvsEta[2]->Fill(llEta, llPt);
            }//} else {
          }//if(abs(ll_genid3) == 11 && ll_mom_genid3 == 23){
        }//if(abs(ll_genid) == 22 && abs(ll_mom_genid) == 11) {


      }//cms2.hyp loop
      
    }//events loop
  }//file loop


  //now do the scaling and the subtracting of the SS and OS
  //distributions to get the estimate of the pure wrong sign dists
   vector<string> v_ch;
   v_ch.push_back("Plus");
   v_ch.push_back("Minus");
   v_ch.push_back("");
   TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
   
   for(unsigned int i = 0; i < 2; i++ ) {
     double nOS = h_osPt[i]->Integral();
     double nSS = h_ssPt[i]->Integral();
     h_osPt[i]  ->Scale(0.5*nSS/nOS); //ignore the error on this scale factor as its very small
     h_osEta[i] ->Scale(0.5*nSS/nOS); //ignore the error on this scale factor as its very small

     
     //2D histo stuff
     h_osPtvsEta[i]->Scale(0.5*nSS/nOS);
     
     h_diffPt[i]      = (TH1F*)h_ssPt[i]->Clone();
     h_diffEta[i]     = (TH1F*)h_ssEta[i]->Clone();
     h_diffPtvsEta[i] = (TH2F*)h_ssPtvsEta[i]->Clone();
     h_diffPt[i]      ->SetName((region+"_diff"+v_ch[i]+"Pt").c_str());
     h_diffEta[i]     ->SetName((region+"_diff"+v_ch[i]+"Eta").c_str());
     h_diffPtvsEta[i] ->SetName((region+"_diff"+v_ch[i]+"PtvsEta").c_str());
     h_diffPt[i]      ->SetDirectory(rootdir);
     h_diffEta[i]      ->SetDirectory(rootdir);
     h_diffPtvsEta[i]   ->SetDirectory(rootdir);
     SubtractHistsCorrectly(h_diffPt[i], h_osPt[i]);
     SubtractHistsCorrectly(h_diffEta[i], h_osEta[i]);
     SubtractHistsCorrectly(h_diffPtvsEta[i], h_osPtvsEta[i]);
     
   }
   
   //now get my prediction and combine the Plus and Minus
   h_diffPt[2] = (TH1F*)h_diffPt[0]->Clone();
   h_diffPt[2]->SetName((region+"_diffPt").c_str());
   h_diffPt[2]->SetTitle("Predicted Pt distribution of Wrong Sign Electrons");
   //  //   h_diffPt[2]->Sumw2();
   h_diffPt[2]->Add(h_diffPt[1]);
   
   //do the same for Eta
   h_diffEta[2] = (TH1F*)h_diffEta[0]->Clone();
   h_diffEta[2]->SetName((region+"_diffEta").c_str());
   h_diffEta[2]->SetTitle("Predicted Eta distribution of Wrong Sign Electrons");
   // //   h_diffEta[2]->Sumw2();
   h_diffEta[2]->Add(h_diffEta[1]);

   //do it for the 2d histo
   h_diffPtvsEta[2] = (TH2F*)h_diffPtvsEta[0]->Clone();
   h_diffPtvsEta[2]->SetName((region+"_diffPtvsEta").c_str());
   h_diffPtvsEta[2]->SetTitle("Predicted Pt vs Eta distribution of Wrong Sign Electrons");
   // //   h_diffPtvsEta[2]->Sumw2();
   h_diffPtvsEta[2]->Add(h_diffPtvsEta[1]);

   //now make the predicted FLR hists
   h_FLRPt= (TH1F*)h_diffPt[2]->Clone();
   h_FLRPt->SetTitle("predicted FLR as a function of electron Pt");
   h_FLRPt->SetName((region+"_FLRPt").c_str());
   TH1F *h_denom = (TH1F*)h_osPt[2]->Clone();
   h_denom->SetName((region+"_denomPt").c_str());
   h_denom->Add(h_ssPt[2]);
   h_FLRPt->SetName((region+"_FLRPt").c_str());
   h_FLRPt->SetTitle("FLR extracted from Z->ee events as a function of Pt");
   DivideHists(h_FLRPt, h_denom);
   
   //Eta
   h_FLREta= (TH1F*)h_diffEta[2]->Clone();
   h_FLREta->SetTitle("predicted FLR as a function of electron Eta");
   h_FLREta->SetName((region+"_FLREta").c_str());
   h_denom = (TH1F*)h_osEta[2]->Clone();
   h_denom->SetName((region+"_denomEta").c_str());
   h_denom->Add(h_ssEta[2]);
   h_FLREta->SetName((region+"_FLREta").c_str());
   h_FLREta->SetTitle("FLR extracted from Z->ee events as a function of Eta");
   DivideHists(h_FLREta, h_denom);
   
   //PtvsEta
   h_FLRPtvsEta = (TH2F*)h_diffPtvsEta[2]->Clone();
   h_FLRPtvsEta->SetTitle("Predicted FLR as a function of the pt and eta of the electron");
   h_FLRPtvsEta->SetName((region+"_FLRPtvsEta").c_str());
   TH2F *h2_denom = (TH2F*)h_osPtvsEta[2]->Clone();
   h2_denom->SetName((region+"_denomPtvsEta").c_str());
   h2_denom->Add(h_ssPtvsEta[2]);
   h_FLRPtvsEta->SetName((region+"_FLRPtvsEta").c_str());
   h_FLRPtvsEta->SetTitle("FLR extracted from Z->ee events as a function of Pt and \eta");
   DivideHists(h_FLRPtvsEta, h2_denom);
   
   
   //now make the MC truth FLR hists
   h_denom = (TH1F*)h_mcssPt[2]->Clone();
   h_denom->SetName((region+"_mcdenomPt").c_str());
   h_denom->Add(h_mcosPt[2]);
   h_mcFLRPt  = (TH1F*)h_mcssPt[2]->Clone();
   h_mcFLRPt->SetName((region+"_mcFLRPt").c_str());
   h_mcFLRPt->SetTitle("Truth matched FLR extracted from Z->ee events as a function of Pt");
   DivideHistsBinomial(h_mcFLRPt, h_denom);
   //Eta
   h_denom = (TH1F*)h_mcssEta[2]->Clone();
   h_denom->SetName((region+"_mcdenomEta").c_str());
   h_denom->Add(h_mcosEta[2]);
   h_mcFLREta  = (TH1F*)h_mcssEta[2]->Clone();
   h_mcFLREta ->SetName((region+"_mcFLREta").c_str());
   h_mcFLREta ->SetTitle("Truth matched FLR extracted from Z->ee events as a function of Eta");
   DivideHistsBinomial(h_mcFLREta, h_denom);
   
   //PtvsEta
   h2_denom = (TH2F*)h_mcssPtvsEta[2]->Clone();
   h2_denom->SetName((region+"_mcdenomPtvsEta").c_str());
   h2_denom->Add(h_mcosPtvsEta[2]);
   h_mcFLRPtvsEta  = (TH2F*)h_mcssPtvsEta[2]->Clone();
   h_mcFLRPtvsEta ->SetName((region+"_mcFLRPtvsEta").c_str());
   h_mcFLRPtvsEta ->SetTitle("Truth matched FLR extracted from Z->ee events as a function of Pt and Eta");
   DivideHistsBinomial(h_mcFLRPtvsEta, h2_denom);
   
   
   //calculate the errors on the numerator and denominator of the weighted histograms using the 2d histograms
   for(int binx = 1; binx < h_allPtvsEta->GetNbinsX()+1; binx++) {
     double sumNentries = 0.;
     double sumNentries2 = 0.;
     double avgFLR = 0.;
     double avgFLRErr = 0.;
     for(int biny = 1; biny < h_allPtvsEta->GetNbinsY()+1; biny++) {
       sumNentries += h_allPtvsEta->GetBinContent(binx, biny);
       sumNentries2 += pow(sumNentries,2);
     }
     
     cout << sumNentries << "  " << h_allEta->GetBinContent(binx) << endl;
     
     avgFLR = getSingleEleFlipRate(20, h_allEta->GetBinCenter(binx) )
       + getSingleEleFlipRate(30, h_allEta->GetBinCenter(binx))
       + getSingleEleFlipRate(40, h_allEta->GetBinCenter(binx))
       + getSingleEleFlipRate(50, h_allEta->GetBinCenter(binx));
     avgFLR = avgFLR/4.;
     
     avgFLRErr = getSingleEleFlipRateError(20, h_allEta->GetBinCenter(binx) )
       + getSingleEleFlipRateError(30, h_allEta->GetBinCenter(binx))
       + getSingleEleFlipRateError(40, h_allEta->GetBinCenter(binx))
       + getSingleEleFlipRateError(50, h_allEta->GetBinCenter(binx));
     avgFLRErr = avgFLRErr/4.;
     
     

     double numErr2 = avgFLRErr*avgFLRErr*sumNentries2 + avgFLR*avgFLR*sumNentries;
     double denomErr2 = sumNentries;
     

     cout << "avgFLR " << avgFLR
	  << " avgFLRErr: " << avgFLRErr
	  << " sumNentries: " << sumNentries
	  << " sumNentries2: " << sumNentries2
	  << " numErr2: " << numErr2
	  << " denomErr2 " << denomErr2 << endl;
     
     double err = sqrt(numErr2/pow(h_allEta->GetBinContent(binx),2) + denomErr2*sumNentries2/pow(h_allEta->GetBinContent(binx),4));
     double temp = h_allEta_seg->GetBinContent(binx);
     temp = temp/(h_allEta->GetBinContent(binx));
     h_allFLREta_seg->SetBinContent(binx, temp);
     h_allFLREta_seg->SetBinError(binx, err);
     
     
   }
   
   if ( nEventsChain != nEventsTotal ) {
     std::cout << "ERROR: number of events from files is not equal to total number of events" << nEventsChain << "  " << nEventsTotal << std::endl;
   }
   return 0;
}
//////////////////////////////////////////////////////////////////////////////

int DYFlipRate::ScanChainDY_BE( TChain* chain, int nEvents, const char* BBfname) {


  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=0;
  if(nEvents==-1)
    nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;

  //prefix is hardwired....this code is only used for events
  //where one electron is in the barrel and the other is
  //in the endcap anyway
  bookHistos("be");

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
//  TFile *fbb = TFile::Open(BBfname.c_str(), "READ");
  TFile *fbb = TFile::Open(BBfname, "READ");
  TH1F *h_FLRPtbb = (TH1F*)fbb->Get("bb_FLRPt");
  h_FLRPtbb->SetName("bb_FLRPt");
  h_FLRPtbb->SetDirectory(rootdir);
  TH1F *h_FLREtabb = (TH1F*)fbb->Get("bb_FLREta");
  h_FLREtabb->SetName("bb_FLREta");
  // //  h_FLREtabb->Sumw2();
  h_FLREtabb->SetDirectory(rootdir);
  TH2F *h_FLRPtvsEtabb = (TH2F*)fbb->Get("bb_FLRPtvsEta");
  h_FLRPtvsEtabb->SetName("h_FLRPtvsEtabb");
  // //  h_FLRPtvsEtabb->Sumw2();
  h_FLRPtvsEtabb->SetDirectory(rootdir);
  fbb->Close();


  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;

      if(nEventsTotal%50000 == 0)
        std::cout << "Event processed = " << nEventsTotal << std::endl;

      if(!passHLTTrigger("HLT_Ele15_SW_L1R"))
        continue;

      for(unsigned int i = 0; i < cms2.hyp_p4().size(); i++) {

        if(cms2.hyp_type()[i] != 3 )
          continue;


        if(cms2.hyp_p4()[i].M() < 76. || cms2.hyp_p4()[i].M() > 106.)
          continue;

        float ltPt  = TMath::Min(cms2.hyp_lt_p4()[i].Pt(), 199.9);
        float ltEta = TMath::Min(fabs(cms2.hyp_lt_p4()[i].Eta()), 2.399);

        float llPt  = TMath::Min(cms2.hyp_ll_p4()[i].Pt(), 199.9);
        float llEta = TMath::Min(fabs(cms2.hyp_ll_p4()[i].Eta()), 2.399);

        unsigned int ltidx = cms2.hyp_lt_index()[i];
        unsigned int llidx = cms2.hyp_ll_index()[i];

        //ask that it is an barrel-endcap event
        if(ltEta > 1.0 && llEta > 1.0)
          continue;

        if(ltEta < 1.0 && llEta < 1.0)
          continue;


        if(!GoodSusyLeptonID(cms2.hyp_lt_id()[i], cms2.hyp_lt_index()[i]))
          continue;
        if(!GoodSusyLeptonID(cms2.hyp_ll_id()[i], cms2.hyp_ll_index()[i]))
          continue;
        if(!GoodSusyLeptonWithIsolation(cms2.hyp_lt_id()[i], cms2.hyp_lt_index()[i]))
          continue;
        if(!GoodSusyLeptonWithIsolation(cms2.hyp_ll_id()[i], cms2.hyp_ll_index()[i]))
          continue;

        if(conversionElectron(ltidx))
            continue;
          if(conversionElectron(llidx))
            continue;
          if(isChargeFlip(ltidx))
            continue;

          if(isChargeFlip(llidx))
            continue;


        //eta and pt of barrel electron
        float bEta = ltEta < 1.0 ? ltEta : llEta;
        float bPt  = ltEta < 1.0 ? ltPt  : llPt;
        float eEta = ltEta > 1.0 ? ltEta : llEta;
        float ePt  = ltEta > 1.0 ? ltPt  : llPt;
        float weight = h_FLRPtbb->GetBinContent(h_FLRPtbb->FindBin(bPt));
        h_APt->Fill(ePt, weight);
        h_BPt->Fill(ePt, 1-2*weight);
        weight = h_FLREtabb->GetBinContent(h_FLREtabb->FindBin(bEta));
        h_AEta->Fill(eEta, weight);
        h_BEta->Fill(eEta, 1-2*weight);
        weight = h_FLRPtvsEtabb->GetBinContent(h_FLRPtvsEtabb->FindBin(bEta, bPt));
        h_APtvsEta->Fill(eEta, ePt, weight);
        h_BPtvsEta->Fill(eEta, ePt, 1-2*weight);


        //OS case
        if(cms2.hyp_lt_id()[i]*cms2.hyp_ll_id()[i] < 0) {
          if(ltEta > 1.0) {
            h_osPtv->Fill(ltPt);
            h_osEtav->Fill(ltEta);
            h_osPtvsEtav->Fill(ltEta, ltPt);
          }
          if(llEta > 1.0) {
            h_osPtv->Fill(llPt);
            h_osEtav->Fill(llEta);
            h_osPtvsEtav->Fill(llEta, llPt);
          }
        } else {//SS case
          if(ltEta > 1.0) {
            h_ssPtv->Fill(ltPt);
            h_ssEtav->Fill(ltEta);
            h_ssPtvsEtav->Fill(ltEta, ltPt);
          }
          if(llEta > 1.0) {
            h_ssPtv->Fill(llPt);
            h_ssEtav->Fill(llEta);
            h_ssPtvsEtav->Fill(llEta, llPt);
          }

        }//SS case

        //////////////////////////////////////////////////////////////////////////////
        //now do the MC matching
        //////////////////////////////////////////////////////////////////////////////
        int lt_genid  = cms2.els_mc_id()[cms2.hyp_lt_index()[i]];
        int ll_genid  = cms2.els_mc_id()[cms2.hyp_ll_index()[i]];
        int lt_genid3 = cms2.els_mc3_id()[cms2.hyp_lt_index()[i]];
        int ll_genid3 = cms2.els_mc3_id()[cms2.hyp_ll_index()[i]];
        int lt_mom_genid = cms2.els_mc_motherid()[cms2.hyp_lt_index()[i]];
        int ll_mom_genid = cms2.els_mc_motherid()[cms2.hyp_ll_index()[i]];
        int lt_mom_genid3 = cms2.els_mc3_motherid()[cms2.hyp_lt_index()[i]];
        int ll_mom_genid3 = cms2.els_mc3_motherid()[cms2.hyp_ll_index()[i]];


        if(abs(lt_genid) == 11 && lt_mom_genid == 23) { //if matched to electron from Z
          if(cms2.hyp_lt_id()[i]*lt_genid < 0) { //is opposite sign from the mc particle, so the wrong sign (ss)
            if(ltEta > 1.0) {
              h_mcPtW  ->Fill(ltPt);
              h_mcEtaW ->Fill(ltEta);
              h_mcPtvsEtaW->Fill(ltEta, ltPt);
            }
          } else {
            if(ltEta > 1.0) {
              h_mcPtC ->Fill(ltPt);
              h_mcEtaC ->Fill(ltEta);
              h_mcPtvsEtaC->Fill(ltEta, ltPt);
            }
          }// } else {
        }//if matched to electron from Z


        if(abs(lt_genid) == 22 && abs(lt_mom_genid) == 11) {
          //cout << "tight from gamma" << endl; //lots here!!!!->look at mc3
          if(abs(lt_genid3) == 11 && lt_mom_genid3 == 23) {
            if(cms2.hyp_lt_id()[i]*lt_genid3 < 0) { //the gen particle and reco are the opposite sign, so Wrong sign
              if(ltEta > 1.0) {
                h_mcPtW  ->Fill(ltPt);
                h_mcEtaW ->Fill(ltEta);
                h_mcPtvsEtaW->Fill(ltEta, ltPt);
              }
            } else {//opposite sign, so Correct
              if(ltEta > 1.0) {
                h_mcPtC  ->Fill(ltPt);
                h_mcEtaC ->Fill(ltEta);
                h_mcPtvsEtaC->Fill(ltEta, ltPt);
              }
            }//} else {//opposite sign
          }
        }



        if(abs(ll_genid) == 11 && ll_mom_genid == 23) { //if matched to electron from Z
          if(cms2.hyp_ll_id()[i]*ll_genid < 0) { //is opposite sign from the mc particle, so Wrong Sign
            if(llEta > 1.0) {
              h_mcPtW  ->Fill(llPt);
              h_mcEtaW ->Fill(llEta);
              h_mcPtvsEtaW->Fill(llEta, llPt);
            }
          } else {// is same sign as the mc particle
            if(llEta > 1.0) {
              h_mcPtC  ->Fill(llPt);
              h_mcEtaC ->Fill(llEta);
              h_mcPtvsEtaC->Fill(llEta, llPt);
            }
          }//} else {
        }

        if(abs(ll_genid) == 22 && abs(ll_mom_genid) == 11) {
          //cout << "loose from gamma" << endl;
          if(abs(ll_genid3) == 11 && ll_mom_genid3 == 23){
            if(cms2.hyp_ll_id()[i]*ll_genid3 < 0) { //is oppsite sign from the MC particle, so Wrong Sign
              if(llEta > 1.0) {
                h_mcPtW  ->Fill(llPt);
                h_mcEtaW ->Fill(llEta);
                h_mcPtvsEtaW->Fill(llEta, llPt);
              }
            } else {
              if(llEta > 1.0) {
                h_mcPtC  ->Fill(llPt);
                h_mcEtaC ->Fill(llEta);
                h_mcPtvsEtaC->Fill(llEta, llPt);
              }
            }//} else {
          }//if(abs(ll_genid3) == 11 && ll_mom_genid3 == 23){
        }//if(abs(ll_genid) == 22 && abs(ll_mom_genid) == 11) {

      }//hyp loop
    }//events loop
  }//file loop

  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }


  h_mcFLRPt = (TH1F*)h_mcPtW->Clone();
  h_mcFLRPt->SetName("be_mcFLRPt");
  TH1F *h_temp = (TH1F*)h_mcPtW->Clone();
  h_temp->SetName("be_denommcPt");
  h_temp->Add(h_mcPtC);
  DivideHistsBinomial(h_mcFLRPt, h_temp);
  h_mcFLRPt->SetDirectory(rootdir);


  h_mcFLREta = (TH1F*)h_mcEtaW->Clone();
  h_mcFLREta->SetName("be_mcFLREta");
  h_temp = (TH1F*)h_mcEtaW->Clone();
  h_temp->SetName("h_denommcEta");
  h_temp->Add(h_mcEtaC);
  DivideHistsBinomial(h_mcFLREta, h_temp);
  h_mcFLREta->SetDirectory(rootdir);


  h_mcFLRPtvsEta = (TH2F*)h_mcPtvsEtaW->Clone();
  h_mcFLRPtvsEta->SetName("be_mcFLRPtvsEta");
  TH2F *h2_temp = (TH2F*)h_mcPtvsEtaW->Clone();
  h2_temp->SetName("be_denommcPtvsEta");
  h2_temp->Add(h_mcPtvsEtaC);
  DivideHistsBinomial(h_mcFLRPtvsEta, h2_temp);
  h_mcFLRPtvsEta->SetDirectory(rootdir);

  h_mcFLRPt->SetLineColor(kRed);
  h_mcFLREta->SetLineColor(kRed);

  delete h_temp;
  delete h2_temp;


  // /////////////////////////////////////////////////
  // New shit
  // /////////////////////////////////////////////////
  h_FLRPt = (TH1F*)h_ssPtv->Clone();
  h_FLRPt->SetName("be_FLRPt");
  SubtractHistsCorrectly(h_FLRPt, h_APt);
  DivideHists(h_FLRPt, h_BPt);
  //h_FLRPt->Add(h_APt, -1);
  //h_FLRPt->Divide(h_BPt);

  h_FLREta = (TH1F*)h_ssEtav->Clone();
  h_FLREta->SetName("be_FLREta");
  SubtractHistsCorrectly(h_FLREta, h_AEta);
  DivideHists(h_FLREta, h_BEta);
  //h_FLREta->Add(h_AEta, -1);
  //h_FLREta->Divide(h_BEta);

  h_FLRPtvsEta = (TH2F*)h_ssPtvsEtav->Clone();
  h_FLRPtvsEta->SetName("be_FLRPtvsEta");
  SubtractHistsCorrectly(h_FLRPtvsEta, h_APtvsEta);
  DivideHists(h_FLRPtvsEta, h_BPtvsEta);
  //h_FLRPtvsEta->Add(h_APtvsEta, -1);
  //h_FLRPtvsEta->Divide(h_BPtvsEta);


  return 0;

}

//////////////////////////////////////////////////////////////////////////////
int DYFlipRate::getClosestStatus3ParticleIdx(LorentzVector p4) {

  float dR = 0.3;
  int closestIdx = -999;
  for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {
    double temp = dRbetweenVectors(p4, cms2.genps_p4()[i]);
    if(temp < dR) {
      dR = temp;
      closestIdx = i;
    }
  }
  return closestIdx;
}

//////////////////////////////////////////////////////////////////////////////
void DYFlipRate::SubtractHistsCorrectly(TH1F* ss, TH1F *os) {


  if(ss->GetNbinsX() != os->GetNbinsX() ) {
    cout << "Something is wrong!!! Exiting" << endl;
    return;
  }

  for(int i = 0; i <= ss->GetNbinsX(); i++) {
    double temp = ss->GetBinContent(i) - os->GetBinContent(i);
    double err = sqrt(pow(ss->GetBinError(i),2) + pow(os->GetBinError(i), 2));
    ss->SetBinContent(i, temp);
    ss->SetBinError(i, err);
  }
}
//////////////////////////////////////////////////////////////////////////////
//doesn't take into correlations.....
void DYFlipRate::SubtractHistsCorrectly(TH2F* ss, TH2F *os) {


  if(ss->GetNbinsX() != os->GetNbinsX() ) {
    cout << "Something is wrong!!! Exiting" << endl;
    return;
  }

  for(int ix = 0; ix <= ss->GetNbinsX(); ix++) {
    for(int iy = 0; iy <= ss->GetNbinsY(); iy++){
      double temp = ss->GetBinContent(ix, iy) - os->GetBinContent(ix, iy);
      double err = sqrt(pow(ss->GetBinError(ix,iy),2) + pow(os->GetBinError(ix,iy), 2));
      ss->SetBinContent(ix, iy, temp);
      ss->SetBinError(ix, iy, err);
    }//y loop
  }//x loop
}//void SubtractHistsCorrectly(TH2F* ss, TH2F *os) {
//////////////////////////////////////////////////////////////////////////////
void DYFlipRate::DivideHists(TH1F *num, TH1F *denom) {

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
//////////////////////////////////////////////////////////////////////////////
void DYFlipRate::DivideHists(TH2F *num, TH2F *denom) {

  for(int ix = 0; ix <= num->GetNbinsX(); ix++) {
    for(int iy = 0; iy <= num->GetNbinsY(); iy++) {
      double valnum   = num->GetBinContent(ix, iy);
      double valdenom = denom->GetBinContent(ix, iy);
      double sigmanum2   = pow(num->GetBinError(ix,iy),2);
      double sigmadenom2 = pow(denom->GetBinError(ix,iy),2);
      double err = 0.0;
      if(valdenom < 0.00000001)
        num->SetBinContent(ix,iy,0);
      else {
        err = sqrt(sigmanum2/(valdenom*valdenom) + sigmadenom2*pow(valnum,2)/pow(valdenom,4));
        num->SetBinContent(ix,iy,valnum/valdenom);
        num->SetBinError(ix,iy,err);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
void DYFlipRate::DivideHistsBinomial(TH1F *num, TH1F *denom) {

  for(int i = 0; i <= num->GetNbinsX(); i++) {
    double valnum   = num->GetBinContent(i);
    double valdenom = denom->GetBinContent(i);
    double err = 0.0;
    if(valdenom < 0.00000001)
      num->SetBinContent(i,0);
    else {
      double p = valnum/valdenom;
      err = sqrt(p*(1-p)/valdenom);
      num->SetBinContent(i,valnum/valdenom);
      num->SetBinError(i,err);
    }

  }
}
//////////////////////////////////////////////////////////////////////////////
void DYFlipRate::DivideHistsBinomial(TH2F *num, TH2F *denom) {

  for(int ix = 0; ix <= num->GetNbinsX(); ix++) {
    for(int iy = 0; iy <= num->GetNbinsY(); iy++) {
      double valnum   = num->GetBinContent(ix, iy);
      double valdenom = denom->GetBinContent(ix, iy);
      double err = 0.0;
      if(valdenom < 0.00000001)
        num->SetBinContent(ix,iy,0);
      else {
        double p = valnum/valdenom;
        err = sqrt(p*(1-p)/valdenom);
        num->SetBinContent(ix,iy,valnum/valdenom);
        num->SetBinError(ix,iy,err);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////
