/* Usage:
   root [0] .L ScanChain.C++
   root [1] TChain *chain = new TChain("Events")
   root [2] chain->Add("merged_ntuple.root")
   root [3] ScanChain(chain)
*/

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "../CORE/mcSelections.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/triggerUtils.cc"
#include "../CORE/jetSelections.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/trackSelections.cc"

/*
#include "../CORE_Nov05_11pb/muonSelections.cc"
#include "../CORE_Nov05_11pb/electronSelections.cc"
#include "../CORE_Nov05_11pb/electronSelectionsParameters.cc"
*/

// CMS2
#include "CMS2.cc"
using namespace tas;

void progress( int nEventsTotal, int nEventsChain ){
  int period = 1000;
  if(nEventsTotal%1000 == 0) {
    // xterm magic from L. Vacavant and A. Cerri
    if (isatty(1)) {
      if( ( nEventsChain - nEventsTotal ) > period ){
        float frac = (float)nEventsTotal/(nEventsChain*0.01);
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
             "\033[0m\033[32m <---\033[0m\015", frac);
        fflush(stdout);
      }
      else {
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", 100.);
        cout << endl;
      }
    }
  }
}

double effmodel(int id, float pt) {
 double eff = 0.0;
 if (abs(id) == 11) {
  if (pt > 80) eff = 0.818;
//  if (pt > 40 && pt <= 80) eff = 0.68 + 0.0018*pt;
//  if (pt > 10 && pt <= 40) eff = 0.33 + 0.011*pt;

 if (pt > 10 && pt <=80) eff = 0.125 + 0.0317*pt - 0.00056*pt*pt + 4.34e-06*pt*pt*pt - 1.216e-08*pt*pt*pt*pt;
 }

 if (abs(id) == 13) {
  if (pt > 70) eff = 0.852;
//  if (pt > 40 && pt <= 70) eff = 0.772 + 0.001135*pt;
//  if (pt > 10 && pt <= 40) eff = 0.642 + 0.0042*pt;
  if (pt > 10 && pt <=70) eff = 0.605 + 0.0076*pt -7.8e-05*pt*pt + 3.155e-07*pt*pt*pt - 4.25e-10*pt*pt*pt*pt;
 }

 return eff;
}

bool isGoodJet(LorentzVector jetp4, double ptCut, double absEtaCut, double dRCut, bool muJetClean, bool applyAlignmentCorrection, bool removedEtaCutInEndcap){

  if(jetp4.Pt() < ptCut) return false;
  if(fabs(jetp4.Eta()) > absEtaCut) return false;

  for (unsigned int i=0; i < els_p4().size(); i++) {
    if (els_p4().at(i).pt() < 10.)     continue;
    if( !pass_electronSelection(i, electronSelection_ssV3, false, false) ) continue;
    double dR_ell = ROOT::Math::VectorUtil::DeltaR(els_p4().at(i),jetp4);
    if (dR_ell < dRCut) return false;
  }

  
  if (muJetClean){
      for (unsigned int i=0; i < mus_p4().size(); i++) {
        if (mus_p4().at(i).pt() < 10.)     continue;
        if( !muonId(i, NominalSSv3) ) continue;
        if (muonIsoValue(i, false) > 0.15) continue;
        double dR_mll = ROOT::Math::VectorUtil::DeltaR(mus_p4().at(i),jetp4);
        if (dR_mll < dRCut) return false;
    }
  }
  return true;
}



int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

  // Example Histograms
  TFile *file = new TFile("test.root","RECREATE");

  int pt_bins   = 48;
  double pt_min = 10.0;
  double pt_max = 250.0;
 
  float ptbins[] = {5, 10, 20, 30, 40, 50, 60, 200};
  unsigned ptnbins = sizeof(ptbins)/sizeof(ptbins[0])-1; 

 
  int eta_bins    = 50;
  double eta_min  = -2.5;
  double eta_max  = 2.5;

  // el pt
  TH1F *h_el_pt_den       = new TH1F("el_pt_den",       "Electron P_{T} denominator",    ptnbins, ptbins);
  TH1F *h_el_pt_numId     = new TH1F("el_pt_numId",     "Electron P_{T} Id",             ptnbins, ptbins);
  TH1F *h_el_pt_numIdIso  = new TH1F("el_pt_numIdIso",  "Electron P_{T} Id & Iso",       ptnbins, ptbins);
  TH1F *h_el_pt_effId     = new TH1F("el_pt_effId",     "Electron Efficiency P_{T}",     ptnbins, ptbins);
  TH1F *h_el_pt_effIdIso  = new TH1F("el_pt_effIdIso",  "Electron Efficiency P_{T}",     ptnbins, ptbins);

// barrel

  TH1F *h_el_barrelpt_den       = new TH1F("el_barrelpt_den",       "Electron P_{T} denominator",    ptnbins, ptbins);
  TH1F *h_el_barrelpt_numId     = new TH1F("el_barrelpt_numId",     "Electron P_{T} Id",             ptnbins, ptbins);
  TH1F *h_el_barrelpt_numIdIso  = new TH1F("el_barrelpt_numIdIso",  "Electron P_{T} Id & Iso",       ptnbins, ptbins);
  TH1F *h_el_barrelpt_effId     = new TH1F("el_barrelpt_effId",     "Electron Efficiency P_{T}",     ptnbins, ptbins);
  TH1F *h_el_barrelpt_effIdIso  = new TH1F("el_barrelpt_effIdIso",  "Electron Efficiency P_{T}",     ptnbins, ptbins);
//endcap

  TH1F *h_el_endcappt_den       = new TH1F("el_endcappt_den",       "Electron P_{T} denominator",    ptnbins, ptbins);
  TH1F *h_el_endcappt_numId     = new TH1F("el_endcappt_numId",     "Electron P_{T} Id",             ptnbins, ptbins);
  TH1F *h_el_endcappt_numIdIso  = new TH1F("el_endcappt_numIdIso",  "Electron P_{T} Id & Iso",       ptnbins, ptbins);
  TH1F *h_el_endcappt_effId     = new TH1F("el_endcappt_effId",     "Electron Efficiency P_{T}",     ptnbins, ptbins);
  TH1F *h_el_endcappt_effIdIso  = new TH1F("el_endcappt_effIdIso",  "Electron Efficiency P_{T}",     ptnbins, ptbins);


  // el eta
  TH1F *h_el_eta_den      = new TH1F("el_eta_den",      "Electron #eta denominator",  eta_bins, eta_min, eta_max);
  TH1F *h_el_eta_numId    = new TH1F("el_eta_numId",    "Electron #eta Id",           eta_bins, eta_min, eta_max);
  TH1F *h_el_eta_numIdIso = new TH1F("el_eta_numIdIso", "Electron #eta Id & Iso",     eta_bins, eta_min, eta_max);
  TH1F *h_el_eta_effId    = new TH1F("el_eta_effId",    "Electron Efficiency #eta",   eta_bins, eta_min, eta_max);
  TH1F *h_el_eta_effIdIso = new TH1F("el_eta_effIdIso", "Electron Efficiency #eta",   eta_bins, eta_min, eta_max);

  // mu pt
  TH1F *h_mu_pt_den       = new TH1F("mu_pt_den",       "Muon P_{T} denominator",        ptnbins, ptbins);
  TH1F *h_mu_pt_numId     = new TH1F("mu_pt_numId",     "Muon P_{T} Id",                 ptnbins, ptbins);
  TH1F *h_mu_pt_numIdIso  = new TH1F("mu_pt_numIdIso",  "Muon P_{T} Id & Iso",           ptnbins, ptbins);
  TH1F *h_mu_pt_effId     = new TH1F("mu_pt_effId",     "Muon Efficiency P_{T}",         ptnbins, ptbins);
  TH1F *h_mu_pt_effIdIso  = new TH1F("mu_pt_effIdIso",  "Muon Efficiency P_{T}",   ptnbins, ptbins);

// barrel

  TH1F *h_mu_barrelpt_den       = new TH1F("mu_barrelpt_den",       "Muon P_{T} denominator",        ptnbins, ptbins);
  TH1F *h_mu_barrelpt_numId     = new TH1F("mu_barrelpt_numId",     "Muon P_{T} Id",                 ptnbins, ptbins);
  TH1F *h_mu_barrelpt_numIdIso  = new TH1F("mu_barrelpt_numIdIso",  "Muon P_{T} Id & Iso",           ptnbins, ptbins);
  TH1F *h_mu_barrelpt_effId     = new TH1F("mu_barrelpt_effId",     "Muon Efficiency P_{T}",         ptnbins, ptbins);
  TH1F *h_mu_barrelpt_effIdIso  = new TH1F("mu_barrelpt_effIdIso",  "Muon Efficiency P_{T}",   ptnbins, ptbins);

// encap

  TH1F *h_mu_endcappt_den       = new TH1F("mu_endcappt_den",       "Muon P_{T} denominator",        ptnbins, ptbins);
  TH1F *h_mu_endcappt_numId     = new TH1F("mu_endcappt_numId",     "Muon P_{T} Id",                 ptnbins, ptbins);
  TH1F *h_mu_endcappt_numIdIso  = new TH1F("mu_endcappt_numIdIso",  "Muon P_{T} Id & Iso",           ptnbins, ptbins);
  TH1F *h_mu_endcappt_effId     = new TH1F("mu_endcappt_effId",     "Muon Efficiency P_{T}",         ptnbins, ptbins);
  TH1F *h_mu_endcappt_effIdIso  = new TH1F("mu_endcappt_effIdIso",  "Muon Efficiency P_{T}",   ptnbins, ptbins);


  // mu eta
  TH1F *h_mu_eta_den      = new TH1F("mu_eta_den",      "Muon #eta denominator",      eta_bins, eta_min, eta_max);
  TH1F *h_mu_eta_numId    = new TH1F("mu_eta_numId",    "Muon #eta Id",               eta_bins, eta_min, eta_max);
  TH1F *h_mu_eta_numIdIso = new TH1F("mu_eta_numIdIso", "Muon #eta Id & Iso",         eta_bins, eta_min, eta_max);
  TH1F *h_mu_eta_effId    = new TH1F("mu_eta_effId",    "Muon Efficieny #eta",        eta_bins, eta_min, eta_max);
  TH1F *h_mu_eta_effIdIso = new TH1F("mu_eta_effIdIso", "Muon Efficiency #eta",       eta_bins, eta_min, eta_max);

  // Jets

  TH1F *h_genSumJet     = new TH1F("genSumJet",      "genSumJet",      40, 0., 1000.);
  TH1F *h_recoSumJet    = new TH1F("recoSumJet",    "recoSumJet",     40, 0., 1000.);
  TH1F *h_effSumJet     = new TH1F("effSumJet",     "effSumJet",     40, 0., 1000.);

// MET

  TH1F *h_genMET     = new TH1F("genMET",      "genMET",      100, 0., 2000.);
  TH1F *h_recoMET    = new TH1F("recoMET",    "recoMET",     100, 0., 2000.);
  TH1F *h_effMET     = new TH1F("effMET",     "effMET",     100, 0., 2000.);

  TH1F *h_resoMET     = new TH1F("resoMET",     "resoMET",     200, 0., 2.);
  TH1F *h_resoSumJetpt     = new TH1F("resoSumJetpt",     "resoSumJetpt",     200, 0., 2.);

  TH1F *Njet_ee       = new TH1F("Njet_ee",       "Njet_ee",    0, 0, 5);
  TH1F *Njet_mm       = new TH1F("Njet_mm",       "Njet_mm",    0, 0, 5);
  TH1F *Njet_em       = new TH1F("Njet_em",       "Njet_em",    0, 0, 5);
  TH1F *Njet_all       = new TH1F("Njet_all",       "Njet_all",    0, 0, 5);

  std::vector<std::string> jetcorr_filenames;
  jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt");
  jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");

  FactorizedJetCorrector *jet_corrector = makeJetCorrector(jetcorr_filenames);

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
    
      // Get Event Content
      cms2.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      progress( nEventsTotal, nEventsChain );
      if ( nEventsTotal%10000 == 0 ) {
        std::cout << "Event: " << nEventsTotal << endl;
      }

/*
      // Gen particles only

      float weight = 1.;
      bool data = false;
      if(!data) {
//          weight = evt_scale1fb()*0.035*157.5/(evt_xsec_incl()*evt_kfactor());
          weight = evt_scale1fb()*0.035;   
      }

     // Basic trigger

//     bool passMu = passHLTTrigger("HLT_Mu9");
//     bool passEl = passHLTTrigger("HLT_Ele10_LW_L1R") || passHLTTrigger("HLT_Ele10_LW_EleId_L1R") || passHLTTrigger("HLT_Ele15_LW_L1R");

     vector< pair <float, int > > hypothese;

      for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {//status 3 loop
       if (cms2.genps_p4().at(i).pt() > 10 && fabs(cms2.genps_p4().at(i).eta()) < 2.4 ) {
       if(abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13) {  
         pair<float,int> lep(genps_p4().at(i).pt(), genps_id()[i]); 
         hypothese.push_back(lep);
       }

       if (abs(cms2.genps_id()[i]) == 15 ) {  // taus
          cms2.genps_lepdaughter_id()[i].size(); 
          bool neutrino = false;

            for(unsigned int kk = 0; kk < cms2.genps_lepdaughter_id()[i].size(); kk++) {
               int daughter = abs(cms2.genps_lepdaughter_id()[i][kk]);
               if( daughter == 12 || daughter == 14) neutrino = true;
             }//daughter loop

        // Make sure they pass eta cuts

            for(unsigned int j = 0; j < cms2.genps_lepdaughter_id()[i].size(); j++) { //loop over the tau's status1 daughters
               if((abs(cms2.genps_lepdaughter_id()[i][j]) == 11 || abs(cms2.genps_lepdaughter_id()[i][j]) == 13 )&& fabs(genps_lepdaughter_p4()[i][j].eta()) < 2.4 && neutrino) {

                 pair<float,int> lep(genps_lepdaughter_p4()[i][j].pt(), genps_lepdaughter_id()[i][j]);
                 hypothese.push_back(lep);
                 continue;
//             dumpDocLines();
           } 
         }
        } // tau
       }
      }//status 3 loop

//     if (nele > 1 && nele20 > 0 && (nelchp > 1 || nelchm > 1)) dumpDocLines();
//     if (nmu > 1 && nmu20 > 0 && (nmuchp > 1 || nmuchm > 1)) dumpDocLines();

// Same Sign

  for(unsigned int i = 0; i < hypothese.size(); i++) {
      pair<float,int> lep1 = hypothese[i];
    for(unsigned int j = i+1; j < hypothese.size(); j++) {
      pair<float,int> lep2 = hypothese[j];
      if(TMath::Max(lep1.first,lep2.first) < 20) continue;
      if(TMath::Min(lep1.first,lep2.first) < 10) continue;
      if (lep1.second * lep2.second < 0) continue;
// OS
//      if (lep1.second * lep2.second > 0) continue;
      float wt = weight;   
//        float wt = effmodel(lep1.second, lep1.first)*effmodel(lep2.second, lep2.first)*weight;
         Njet_all->Fill(4, wt);   
         if (abs(lep1.second) == 11 && abs(lep2.second) ==11 ) Njet_ee->Fill(1, wt);
         else if (abs(lep1.second) == 13 && abs(lep2.second) ==13 ) Njet_mm->Fill(2, wt);
         else Njet_em->Fill(3, wt);
//       cout << lep1.first << "  " << lep2.first << endl;
//       cout << lep1.second << "  " << lep2.second << endl;
   }
  }

// Get the jets
     vector<LorentzVector> v_jetP4s;
     double sumet_calo = 0.0;
     double sumecalo = 0.0;

        for (unsigned int i = 0; i < pfjets_p4().size(); i++)
         {
            LorentzVector jp4 = pfjets_p4()[i];
            float jet_cor = jetCorrection(jp4, jet_corrector);
            v_jetP4s.push_back(jp4 * jet_cor);
         }

        for (unsigned int j = 0; j < v_jetP4s.size(); ++j) {
         if ( !passesPFJetID(j)) continue;
         if (isGoodJet(v_jetP4s.at(j), 30, 2.5, 0.4, true, false, false)) {
            sumet_calo += v_jetP4s[j].Pt();
            sumecalo += v_jetP4s[j].E();
         }
       }


// MET and HT Turn ON curve

// evt_pfmet

     float genMETpx = 0.0;
     float genMETpy = 0.0;
     float gensumpt = 0.0;
     float gensume = 0.0;

//     dumpDocLines();
     
      for(unsigned int i = 0; i < cms2.genps_id().size(); i++) {//status 3 loop
        if(abs(cms2.genps_id()[i]) == 1000022 || abs(cms2.genps_id()[i]) == 12 || abs(cms2.genps_id()[i]) == 14 || abs(cms2.genps_id()[i]) == 16) {
           genMETpx += genps_p4().at(i).px();
           genMETpy += genps_p4().at(i).py();
        }
       if(abs(cms2.genps_id()[i]) < 6 || abs(cms2.genps_id()[i]) == 21 ) {
          if (genps_p4().at(i).pt() > 30 && fabs(genps_p4().at(i).eta()) < 2.5) {
          gensumpt += genps_p4().at(i).pt();
          gensume += genps_p4().at(i).E();
         }
        }
      } 

//    dumpDocLines();     

//  RECO pfJet/MET

     h_genMET->Fill(sqrt((genMETpx*genMETpx) + (genMETpy*genMETpy)), weight);
     h_genSumJet->Fill(gensumpt,  weight);

     if (evt_pfmet() > 80) {
         h_recoMET->Fill(sqrt((genMETpx*genMETpx) + (genMETpy*genMETpy)), weight);
         h_resoMET->Fill(evt_pfmet()/(sqrt((genMETpx*genMETpx) + (genMETpy*genMETpy))), weight);
       }
     if (sumet_calo > 350) {
         h_recoSumJet->Fill(gensumpt, weight);
         if (gensumpt > 80) h_resoSumJetpt->Fill(sumet_calo/gensumpt, weight);
      }

 //     if (sumet_calo < 300) continue;
      if (evt_pfmet() < 80) continue;

*/

      // loop over electrons
      for(unsigned int iEl=0; iEl<els_p4().size(); iEl++){
        // reco Pt
        double ptEl1   = els_p4().at(iEl).pt();
//        double etaEl  = els_p4().at(iEl).eta();

        double ptEl  = els_mc_p4().at(iEl).pt();
        double etaEl = els_mc_p4().at(iEl).eta();


        // denominator
        if( ptEl < 10. ) continue;
        if( ptEl1 < 10. ) continue;
        if( fabs(etaEl) > 2.4 ) continue;
        if( leptonIsFromW(iEl,11, true) < 1 ) continue;
        h_el_pt_den->Fill(ptEl);
        h_el_eta_den->Fill(etaEl);
        if (fabs(etaEl) < 1.479 ) h_el_barrelpt_den->Fill(ptEl);
        if (fabs(etaEl) > 1.479 ) h_el_endcappt_den->Fill(ptEl);

        // numerator Id
        if( !pass_electronSelection(iEl, electronSelection_ssV3_noIso, false, false) ) continue;
        h_el_pt_numId->Fill(ptEl);
        h_el_eta_numId->Fill(etaEl);
        if (fabs(etaEl) < 1.479 ) h_el_barrelpt_numId->Fill(ptEl);
        if (fabs(etaEl) > 1.479 ) h_el_endcappt_numId->Fill(ptEl);


        // numerator Id & Iso
        if( !pass_electronSelection(iEl, electronSelection_ssV3, false, false) ) continue;
        h_el_pt_numIdIso->Fill(ptEl);
        h_el_eta_numIdIso->Fill(etaEl);
        if (fabs(etaEl) < 1.479 ) h_el_barrelpt_numIdIso->Fill(ptEl);
        if (fabs(etaEl) > 1.479 ) h_el_endcappt_numIdIso->Fill(ptEl);

      }
      
      // loop over muons
      for(unsigned int iMu=0; iMu<mus_p4().size(); iMu++){

        double ptMu1   = mus_p4().at(iMu).pt();
//        double etaMu  = mus_p4().at(iMu).eta();

        double ptMu   = mus_mc_p4().at(iMu).pt();
        double etaMu  = mus_mc_p4().at(iMu).eta();


        // denominator
        if( ptMu < 10 ) continue;
        if( ptMu1 < 10 ) continue;
        if( fabs(etaMu) > 2.4 ) continue;
        if( leptonIsFromW(iMu,13,true) < 1 ) continue;
        h_mu_pt_den->Fill(ptMu);
        h_mu_eta_den->Fill(etaMu);
        if (fabs(etaMu) < 1.479 ) h_mu_barrelpt_den->Fill(ptMu);
        if (fabs(etaMu) > 1.479 ) h_mu_endcappt_den->Fill(ptMu);

        // numerator Id
        //if( !muonIdNotIsolated(iMu, NominalTTbarV2) ) continue;
        if( !muonId(iMu, NominalSSv3) ) continue;
        h_mu_pt_numId->Fill(ptMu);
        h_mu_eta_numId->Fill(etaMu);
        if (fabs(etaMu) < 1.479 ) h_mu_barrelpt_numId->Fill(ptMu);
        if (fabs(etaMu) > 1.479 ) h_mu_endcappt_numId->Fill(ptMu);

        // numerator Id & Iso
        //if( !muonId(iMu, NominalTTbarV2) ) continue;
        if( !muonId(iMu, NominalSSv3) ) continue;
        if (muonIsoValue(iMu, false) > 0.15) continue;
        // if( muonIsoValueFlorida(iMu) > 0.15 ) continue;

        h_mu_pt_numIdIso->Fill(ptMu);
        h_mu_eta_numIdIso->Fill(etaMu);

        if (fabs(etaMu) < 1.479 ) h_mu_barrelpt_numIdIso->Fill(ptMu);
        if (fabs(etaMu) > 1.479 ) h_mu_endcappt_numIdIso->Fill(ptMu);


      }

    }
  
    delete tree;
    f.Close();
  }
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  // errors
  h_el_pt_effId->Sumw2();
  h_el_pt_effIdIso->Sumw2();

  h_el_barrelpt_effId->Sumw2();
  h_el_barrelpt_effIdIso->Sumw2();
  h_el_endcappt_effId->Sumw2();
  h_el_endcappt_effIdIso->Sumw2();

  h_el_eta_effId->Sumw2();
  h_el_eta_effIdIso->Sumw2();
  h_mu_pt_effId->Sumw2();
  h_mu_pt_effIdIso->Sumw2();

  h_mu_barrelpt_effId->Sumw2();
  h_mu_barrelpt_effIdIso->Sumw2();
  h_mu_endcappt_effId->Sumw2();
  h_mu_endcappt_effIdIso->Sumw2();


  h_mu_eta_effId->Sumw2();
  h_mu_eta_effIdIso->Sumw2();
  h_effMET->Sumw2();
  h_effSumJet->Sumw2();

  // efficiencies
  h_el_pt_effId->Divide( h_el_pt_numId, h_el_pt_den, 1.0, 1.0, "B" ); 
  h_el_pt_effIdIso->Divide( h_el_pt_numIdIso, h_el_pt_den, 1.0, 1.0, "B" ); 

  h_el_barrelpt_effId->Divide( h_el_barrelpt_numId, h_el_barrelpt_den, 1.0, 1.0, "B" ); 
  h_el_barrelpt_effIdIso->Divide( h_el_barrelpt_numIdIso, h_el_barrelpt_den, 1.0, 1.0, "B" );

  h_el_endcappt_effId->Divide( h_el_endcappt_numId, h_el_endcappt_den, 1.0, 1.0, "B" );
  h_el_endcappt_effIdIso->Divide( h_el_endcappt_numIdIso, h_el_endcappt_den, 1.0, 1.0, "B" );


  h_el_eta_effId->Divide( h_el_eta_numId, h_el_eta_den, 1.0, 1.0, "B" ); 
  h_el_eta_effIdIso->Divide( h_el_eta_numIdIso, h_el_eta_den, 1.0, 1.0, "B" ); 

  h_mu_pt_effId->Divide( h_mu_pt_numId, h_mu_pt_den, 1.0, 1.0, "B" ); 
  h_mu_pt_effIdIso->Divide( h_mu_pt_numIdIso, h_mu_pt_den, 1.0, 1.0, "B" ); 

  h_mu_barrelpt_effId->Divide( h_mu_barrelpt_numId, h_mu_barrelpt_den, 1.0, 1.0, "B" ); 
  h_mu_barrelpt_effIdIso->Divide( h_mu_barrelpt_numIdIso, h_mu_barrelpt_den, 1.0, 1.0, "B" );

  h_mu_endcappt_effId->Divide( h_mu_endcappt_numId, h_mu_endcappt_den, 1.0, 1.0, "B" );
  h_mu_endcappt_effIdIso->Divide( h_mu_endcappt_numIdIso, h_mu_endcappt_den, 1.0, 1.0, "B" );


  h_mu_eta_effId->Divide( h_mu_eta_numId, h_mu_eta_den, 1.0, 1.0, "B" ); 
  h_mu_eta_effIdIso->Divide( h_mu_eta_numIdIso, h_mu_eta_den, 1.0, 1.0, "B" ); 

// Jets and MET

  h_effMET->Divide( h_recoMET, h_genMET, 1.0, 1.0, "B" );
  h_effSumJet->Divide( h_recoSumJet, h_genSumJet, 1.0, 1.0, "B" );

  // return
  file->Write();
  return 0;
}
