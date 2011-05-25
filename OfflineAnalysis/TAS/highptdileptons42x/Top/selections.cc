#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "TSystem.h"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/metSelections.cc"
#include "CMS2.cc"

using namespace tas;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);//, bool used0wrtPV = false);
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);//, bool used0wrtPV  = false);
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);//, bool used0wrtPV = false);
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);//, bool used0wrtPV = false);
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);
// std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix);
std::pair<float,float> getMet(const string algo, unsigned int hypIdx);
bool inZmassWindow (float mass);
bool passTriggersMu9orLisoE15(int dilType);
int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug, bool usePtOnlyForWeighting);
bool isFakeableMuon(int index);
double getd0wrtPV(LorentzVector p4, float d0);
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi);
bool isFakeableElectron(int index, string prefix, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);
bool additionalZvetoSUSY2010 (int i_hyp, bool applyAlignmentCorrection, bool removedEtaCutInEndcap);
bool passEGTrigger(bool mc, int type);
bool passMuTrigger(bool mc, int type);
int nHLTObjects(string arg);
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps);
int getNbtags(const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator);
std::pair<int,float> getMinbltags(unsigned int hypIdx, const vector<LorentzVector> v_jetP4s, const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator);
bool sortLVByPt(LorentzVector lv1, LorentzVector lv2);
LorentzVector p4HLTObject(string arg, int) ;


Bool_t sortLVByPt(LorentzVector lv1, LorentzVector lv2) {
  return lv1.pt() > lv2.pt();
}


bool isFakeableElectron (int index, string prefix, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {
  TPMERegexp re1("v1", "g");
  TPMERegexp re2("v2", "g");
  TPMERegexp re3("v3", "g");

  if (re1.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v1, applyAlignmentCorrection, removedEtaCutInEndcap);
  if (re2.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v2, applyAlignmentCorrection, removedEtaCutInEndcap);
  if (re3.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v3, applyAlignmentCorrection, removedEtaCutInEndcap);

  return false;
}

bool isFakeableMuon (int index) {
     cout << "This is for ttbar check the ISO values" << endl;
     return muonId(index, muonSelectionFO_mu_ttbar);
}


/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {//, bool used0wrtPV) {

  if(abs(id) == 11) {
    return (pass_electronSelection(lepIdx, electronSelection_ttbarV2_noiso, applyAlignmentCorrection, removedEtaCutInEndcap));
  }


  if(abs(id) == 13) {
    if ( mus_p4()[lepIdx].pt() < 5.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    } 
    return muonIdNotIsolated(lepIdx, NominalTTbarV2);
  }

  return true;
}

/******************************************************************************************/
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) { //, bool used0wrtPV) {

// Covers both the noiso part

  if(!isGoodLeptonNoIso(id, lepIdx, applyAlignmentCorrection, removedEtaCutInEndcap)) return false;

  if(abs(id)== 11) {
       if (!pass_electronSelection(lepIdx, electronSelection_ttbarV2, applyAlignmentCorrection, removedEtaCutInEndcap)) return false;
  }

  // 13 is a muon
  if(abs(id) == 13)
    if(muonIsoValue(lepIdx) > 0.15)   return false;
  return true;
}


bool additionalZvetoSUSY2010(int i_hyp, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i < mus_p4().size(); i++) {
    bool hypLep1 = false;
    if (mus_p4().at(i).pt() < 10.)     continue;

    if (!isGoodLeptonNoIso(13, i, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;

    if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 13 && hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 13 && hyp_ll_index()[i_hyp] == i ) hypLep1 = true;
    
    for (unsigned int j=i+1; j < mus_p4().size(); j++) {
      bool hypLep2 = false;
      if (mus_p4().at(j).pt() < 10.) continue;

      if (!isGoodLeptonNoIso(13, j, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;

      if (mus_charge().at(i) == mus_charge().at(j)) continue;
      if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 13 && hyp_lt_index()[i_hyp] == j ) hypLep2 = true;
      if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 13 && hyp_ll_index()[i_hyp] == j ) hypLep2 = true;

      if ((muonIsoValue(i) > 0.10) && (muonIsoValue(j) > 0.10))   continue;	
      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;
      // Make the invariant mass
      LorentzVector vec = mus_p4().at(i) + mus_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;
    }
  }

  // now, look for Z->ee
  for (unsigned int i=0; i < els_p4().size(); i++) {
    bool hypLep1 = false;
    if (els_p4().at(i).pt() < 10.)     continue;

    if (!isGoodLeptonNoIso(11, i, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;

    if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 11 && hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 11 && hyp_ll_index()[i_hyp] == i ) hypLep1 = true;

    for (unsigned int j=i+1; j<els_p4().size(); j++) {
      bool hypLep2 = false;
      if (els_p4().at(j).pt() < 10.) continue;

      if (!isGoodLeptonNoIso(11, j, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;
      if (els_charge().at(i) == els_charge().at(j)) continue;

      if ((!pass_electronSelection(i, electronSelection_ss_Iso, applyAlignmentCorrection, removedEtaCutInEndcap)) && (!pass_electronSelection(j, electronSelection_ss_Iso, applyAlignmentCorrection, removedEtaCutInEndcap))) continue;

      if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 11 && hyp_lt_index()[i_hyp] == j ) hypLep2 = true;
      if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 11 && hyp_ll_index()[i_hyp] == j ) hypLep2 = true;

      if ( hypLep1 && hypLep2 ) continue;
      if ( !hypLep1 && !hypLep2 ) continue;

      // Make the invariant mass
      LorentzVector vec = els_p4().at(i) + els_p4().at(j);
      if ( inZmassWindow(vec.mass()) ) return true;
      
    }
  }

  return veto;
}


/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {//, bool used0wrtPV) {
  
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
    return false;

  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap) {//, bool used0wrtPV) {


  if(!isGoodLeptonwIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
    return false;
  if(!isGoodLeptonwIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap))//, used0wrtPV)
    return false;


  return true;
}

/******************************************************************************************/     
// is it a good jet?
/******************************************************************************************/     
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean){
		     
  if(jetp4.Pt() < ptCut)
    return false;
  if(fabs(jetp4.Eta()) > absEtaCut)
    return false;
  
  double dR_ll = ROOT::Math::VectorUtil::DeltaR(hyp_ll_p4()[hypIdx],jetp4);
  double dR_lt = ROOT::Math::VectorUtil::DeltaR(hyp_lt_p4()[hypIdx],jetp4);
  
  if (abs(hyp_ll_id()[hypIdx]) == 11){
    if (dR_ll < dRCut) return false;
  }
  if (abs(hyp_lt_id()[hypIdx]) == 11){
    if (dR_lt < dRCut) return false;
  }

  if (muJetClean){
    if (abs(hyp_ll_id()[hypIdx]) == 13){
      if (dR_ll < dRCut) return false;
    }
    if (abs(hyp_lt_id()[hypIdx]) == 13){
      if (dR_lt < dRCut) return false;
    }
  }

  return true;

}
/******************************************************************************************/     
//return the MET and the MET phi instead of a bool because the MT2 needs it
/******************************************************************************************/     
/* Old one
std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix) {
 
 if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET") {
    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
    return make_pair(-99999., -99999.);
  }
  if(algo == "tcMET") {
    double tcmet = evt_tcmet();
    double tcmetPhi = evt_tcmetPhi();
    correctTcMETForHypMus(hypIdx, tcmet, tcmetPhi);
    return make_pair(tcmet, tcmetPhi);
  }
  if(algo == "muCorMET")
    return make_pair(evt_metMuonCorr(), evt_metMuonCorrPhi());
  if(algo == "pfMET")
    return make_pair(evt_pfmet(), evt_pfmetPhi());


  return make_pair(-99999., -99999);

}
*/

std::pair<float,float> getMet(const string algo, unsigned int hypIdx) {
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET" && algo != "tcMET35X" && algo != "tcMET_looper") {
    cout << algo << "IS NOT A RECOGNIZED MET ALGORITHM!!!!! PLEASE CHECK YOUR CODE!!!";
    return make_pair(-99999., -99999.);
  }

  
  if(algo == "tcMET") {

    float tcmetX = evt_tcmet()*cos(evt_tcmetPhi());
    float tcmetY = evt_tcmet()*sin(evt_tcmetPhi());
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), tcmetX, tcmetY, usingTcMet);
    
    return make_pair(sqrt(tcmetX*tcmetX + tcmetY*tcmetY), atan2(tcmetY, tcmetX));
  }
  
  if(algo == "muCorMET") {

    float metX = evt_metMuonCorr()*cos(evt_metMuonCorrPhi());
    float metY = evt_metMuonCorr()*sin(evt_metMuonCorrPhi());
    
    if(abs(hyp_lt_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_lt_index().at(hypIdx), metX, metY, usingCaloMet);
    if(abs(hyp_ll_id()[hypIdx]) == 13)
      fixMetForThisMuon(hyp_ll_index().at(hypIdx), metX, metY, usingCaloMet);

    return make_pair(sqrt(metX*metX + metY*metY), atan2(metY, metX));
  }
  
  //nothing to do here because they're perfect
  if(algo == "pfMET") 
    return make_pair(evt_pfmet(), evt_pfmetPhi());
  
  
  return make_pair(-99999., -99999);
  
}


/******************************************************************************************/     
//trigger requirement
/******************************************************************************************/         
bool passTriggersMu9orLisoE15(int dilType) {
  
  //TString method
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_SW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  
  if (dilType == 0 && ! (hltMu9) ) return false;
  if ((dilType == 1 || dilType == 2) && ! (hltMu9 || hlt_ele15_lw_l1r)) return false;
  if (dilType == 3 && ! hlt_ele15_lw_l1r) return false;     

  return true;
}

/******************************************************************************************/     
//hypothesis disabmiguation
/******************************************************************************************/     
int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug, bool usePtOnlyForWeighting){
  
  int result = -1;
  unsigned int nGoodHyps = goodHyps.size();
  if ( nGoodHyps == 0 ) return result;

  float maxWeight = -1;
  unsigned int maxWeightIndex = 9999;
  
  for (unsigned int hypIdxL=0; hypIdxL < nGoodHyps; ++hypIdxL){
    unsigned int hypIdx = goodHyps[hypIdxL];
    float hypWeight_lt = 0;
    float hypWeight_ll = 0;
    float hypWeight_iso = 0;
    float hypWeight = 0;
    unsigned int i_lt = cms2.hyp_lt_index().at(hypIdx);
    unsigned int i_ll = cms2.hyp_ll_index().at(hypIdx);

    int id_lt = cms2.hyp_lt_id().at(hypIdx);
    int id_ll = cms2.hyp_ll_id().at(hypIdx);

    //float isoTk_lt = leptonTrkIsolationTTDil08(id_lt, i_lt);
    //float isoTk_ll = leptonTrkIsolationTTDil08(id_ll, i_ll);
    float isoTk_lt, isoTk_ll;
    float isoCal_lt, isoCal_ll;
    if(abs(id_lt) == 11) {
      isoTk_lt = els_p4()[i_lt].Pt()/(els_tkJuraIso()[i_lt]+els_p4()[i_lt].Pt());
      isoCal_lt = els_p4()[i_lt].Pt()/(els_hcalIso()[i_lt]+els_ecalJuraIso()[i_lt]+els_p4()[i_lt].Pt());
    } else {
      isoTk_lt = mus_p4()[i_lt].Pt()/(mus_iso03_sumPt()[i_lt]+mus_p4()[i_lt].Pt());
      isoCal_lt = mus_p4()[i_lt].Pt()/(mus_iso03_emEt()[i_lt]+mus_iso03_hadEt()[i_lt]+mus_p4()[i_lt].Pt());
    }
    if(abs(id_ll) == 11) {
      isoTk_ll = els_p4()[i_ll].Pt()/(els_tkJuraIso()[i_ll]+els_p4()[i_ll].Pt());
      isoCal_ll = els_p4()[i_ll].Pt()/(els_hcalIso()[i_ll]+els_ecalJuraIso()[i_ll]+els_p4()[i_ll].Pt());
    } else {
      isoTk_ll = mus_p4()[i_ll].Pt()/(mus_iso03_sumPt()[i_ll]+mus_p4()[i_ll].Pt());
      isoCal_ll = mus_p4()[i_ll].Pt()/(mus_iso03_emEt()[i_ll]+mus_iso03_hadEt()[i_ll]+mus_p4()[i_ll].Pt());
    }

  
    
    //ad-hoc selection of weights
    if (abs(id_lt) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25); //shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting && cms2.els_egamma_tightId().at(i_lt)) hypWeight_lt += 0.2;
    }
    if (abs(id_lt) == 13){
      //I want to select "trk & cal"-isolated ones	    
      hypWeight_iso += (isoTk_lt*isoCal_lt - 0.25);//shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting) hypWeight_lt += 0.4;
    }
    if (abs(id_ll) == 11){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting && cms2.els_egamma_tightId().at(i_ll)) hypWeight_ll += 0.2;
    }
    if (abs(id_ll) == 13){
      //I want to select "trk & cal"-isolated ones
      hypWeight_iso *= (isoTk_ll*isoCal_ll - 0.25); //shift by 0.25 to be positive-definite
      if (! usePtOnlyForWeighting) hypWeight_ll += 0.4;
    }
    float pt_lt = cms2.hyp_lt_p4().at(hypIdx).pt();
    float pt_ll = cms2.hyp_ll_p4().at(hypIdx).pt();
    if(pt_lt > 20.)
      hypWeight_lt += (1. - 20./pt_lt*20./pt_lt);
    else
      hypWeight_lt += (1. - 10./pt_lt*10./pt_lt);
    if(pt_ll > 20.)
      hypWeight_ll += (1. - 20./pt_ll*20./pt_ll);
    else
      hypWeight_ll += (1. - 10./pt_ll*10./pt_ll);
    
    if (usePtOnlyForWeighting){
      hypWeight = hypWeight_ll*hypWeight_lt; //again, desire to have both good
    } else {
      hypWeight = hypWeight_ll*hypWeight_lt*hypWeight_iso; //again, desire to have both good
    }

    if (hypWeight > maxWeight){
      maxWeight = hypWeight;
      maxWeightIndex = hypIdx;
    }
  }


  //Now let's implement the Strasbourg-type disambiguation/dispatch
  //ee
//   std::cout << "right before here" << endl;
//   {
//     std::cout << "here" << endl;
//     std::vector<unsigned int> looseEls(0);
//     std::vector<unsigned int> looseMus(0);
//     //for (unsigned int iEl =0; iEl < cms2.els_p4().size(); ++iEl){
//     //if (looseElectronSelectionTTDil08(iEl)){
//     //looseEls.push_back(iEl);
//     //}
//     //}
//     for (unsigned int iMu =0; iMu < cms2.mus_p4().size(); ++iMu){
//       //if (looseMuonSelectionTTDil08(iMu)){
//       looseMus.push_back(iMu);
//     }
//     //}
    
//     bool pass_elec = false;
//     if (looseEls.size()>1){
//       if (cms2.els_charge().at(looseEls[0]) != cms2.els_charge().at(looseEls[1])){
// 	pass_elec = true;
//       }
//     //   if (looseMus.size()>0 && cms2.mus_p4().at(looseMus[0]).pt() > cms2.els_p4().at(looseEls[1]).pt()) pass_elec = false;
// //       if (looseMus.size()>0 && 
// // 	  ( ( muonTrkIsolationTTDil08(looseMus[0]) > electronTrkIsolationTTDil08(looseEls[0]) 
// // 	      && cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0]) )
// // 	    || ( muonTrkIsolationTTDil08(looseMus[0]) > electronTrkIsolationTTDil08(looseEls[1])
// // 		 && cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0]))
// // 	    )
// // 	  ) pass_elec = false; 
//     }
//     bool pass_muon = false;
//     if (looseMus.size()>1){
//       for (unsigned int iMu=0; iMu < looseMus.size(); ++iMu){
// 	for (unsigned int jMu=iMu+1; jMu < looseMus.size(); ++jMu){
// 	  if (cms2.mus_charge().at(looseMus[iMu]) != cms2.mus_charge().at(looseMus[jMu])) pass_muon = true;
// 	}
//       }
//       if (looseEls.size()>0 && cms2.els_p4().at(looseEls[0]).pt() > cms2.mus_p4().at(looseMus[1]).pt()
// 	  && cms2.mus_charge().at(looseMus[1]) != cms2.els_charge().at(looseEls[0])) pass_muon = false;
//     }
//     bool pass_elecmuon = false;
//     if (looseMus.size() > 0 && looseEls.size() > 0){
//       if (! pass_elec && ! pass_muon ){
// 	if (cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0])) pass_elecmuon = true;
// 	if (! pass_elecmuon && looseEls.size()>1){
// 	  if (cms2.mus_charge().at(looseMus[0]) != cms2.els_charge().at(looseEls[0])) pass_elecmuon = true;
// 	}
//       }
//     }

//     unsigned int passStatus = 0;
//     if (pass_muon) passStatus++;
//     if (pass_elecmuon) passStatus++;
//     if (pass_elec) passStatus++;
//     if (passStatus > 1) std::cout<<"ERROR: inconsistent assignment"<<std::endl;
//     if (passStatus == 1){
//       if (pass_muon) strasbourgDilType = 0;
//       if (pass_elecmuon) strasbourgDilType = 1;
//       if (pass_elec) strasbourgDilType = 2;
//     }
//   }

 //  if (printDebug){
//     int genpDilType = genpDileptonType();
//     if (genpDilType>=0 ){ std::cout<<"Dil type "<<genpDilType<<std::endl;
//       if (nGoodHyps > 1){
// 	int maxWeightType = cms2.hyp_type().at(maxWeightIndex);
// 	if ((maxWeightType == 0 && genpDilType == 0)
// 	    || ( (maxWeightType == 1 || maxWeightType == 2) && genpDilType == 1)
// 	    || (maxWeightType == 3 && genpDilType == 2)){
// 	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
// 		   <<" assigned correctly by maxWeight method";
// 	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type().at(goodHyps[iih]);
// 	  std::cout<<std::endl;
// 	} else {
// 	  std::cout<<"Dil type "<<genpDilType<<" ; Strasbourg dil type "<<strasbourgDilType 
// 		   <<" assigned incorrectly by maxWeight method";
// 	  std::cout<<" out of"; for(unsigned int iih=0;iih<nGoodHyps;++iih)std::cout<<" "<<cms2.hyp_type().at(goodHyps[iih]);
// 	  std::cout<<std::endl;	    
// 	}
//       }
//     }
//     int nMCTruth = 0;
//     for(unsigned int iih=0;iih<nGoodHyps;++iih) if (matchesMCTruthDilExtended(goodHyps[iih])) nMCTruth++;
//     std::cout<<"Ne: "<<genpCountPDGId_Pt20h24(11)<<" nmu: "<<genpCountPDGId_Pt20h24(13)<<" ntau: "<<genpCountPDGId_Pt20h24(15)
// 	     <<" ngood "<<nGoodHyps<<" SBtype "<<strasbourgDilType
// 	     <<" hyp_typeM: "<<cms2.hyp_type()[maxWeightIndex]<<" matchMC "<<(matchesMCTruthDilExtended(maxWeightIndex)? 1 : 0)
// 	     <<" nMatches "<<nMCTruth
// 	     <<" ltid "<< cms2.hyp_lt_id()[maxWeightIndex]
// 	     <<" ltmcid "<< cms2.hyp_lt_mc_id()[maxWeightIndex]<<" ltmcmid "<< cms2.hyp_lt_mc_motherid()[maxWeightIndex]
// 	     <<" llid "<< cms2.hyp_ll_id()[maxWeightIndex]
// 	     <<" llmcid "<< cms2.hyp_ll_mc_id()[maxWeightIndex]<<" llmcmid "<< cms2.hyp_ll_mc_motherid()[maxWeightIndex]
// 	     <<std::endl;    
//   }

  result = maxWeightIndex;
  return result;
}
  

/******************************************************************************************/     
//electron impact parameter with respect to the highest sumpt pv
/******************************************************************************************/     
double getd0wrtPV(LorentzVector p4, float d0) {

   
  double max_sumpt = -1;
  int i_max = -1;
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_isFake().size());
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_position().size());
  assert(cms2.vtxs_sumpt().size() == cms2.vtxs_covMatrix().size());
  for (unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i) {
    if (cms2.vtxs_isFake().at(i))
	       continue;
    if (cms2.vtxs_sumpt().at(i) > max_sumpt) {
      max_sumpt = cms2.vtxs_sumpt().at(i);
      i_max = i;
    }
  }

   if (i_max != -1) {
     const double bx = vtxs_position().at(i_max).x();
     const double by = vtxs_position().at(i_max).y();
     double phi = p4.phi();
     double d0vtx = d0 - bx * sin(phi) + by * cos(phi);
     return d0vtx;
   }

   
   cout << "did not find a PV!!!" << endl;
   return 99999;


 }


//*****************************************************************************************
//correct MET for hyp mus that are not used in MET correction
//*****************************************************************************************
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi){
  if (cms2.hyp_type()[hypIdx] ==3) return;
  double lmetx = met*cos(metPhi);
  double lmety = met*sin(metPhi);

  unsigned int i_lt = cms2.hyp_lt_index()[hypIdx];
  unsigned int i_ll = cms2.hyp_ll_index()[hypIdx];
  if (abs(cms2.hyp_lt_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_lt] == 0){
      lmetx+= - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x();
      lmety+= - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y();
    } else if (cms2.mus_tcmet_flag()[i_lt] == 4){
      lmetx+= - cms2.mus_tcmet_deltax()[i_lt] - cms2.mus_met_deltax()[i_lt] - cms2.mus_p4()[i_lt].x(); 
      lmety+= - cms2.mus_tcmet_deltay()[i_lt] - cms2.mus_met_deltay()[i_lt] - cms2.mus_p4()[i_lt].y(); 
    }
  }
  if (abs(cms2.hyp_ll_id()[hypIdx])==13){
    if(cms2.mus_tcmet_flag()[i_ll] == 0){ 
      lmetx+= - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x(); 
      lmety+= - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y(); 
    } else if (cms2.mus_tcmet_flag()[i_ll] == 4){ 
      lmetx+= - cms2.mus_tcmet_deltax()[i_ll] - cms2.mus_met_deltax()[i_ll] - cms2.mus_p4()[i_ll].x();  
      lmety+= - cms2.mus_tcmet_deltay()[i_ll] - cms2.mus_met_deltay()[i_ll] - cms2.mus_p4()[i_ll].y();  
    } 
  }
  met = sqrt(lmetx*lmetx+lmety*lmety);
  metPhi = atan2(lmety,lmetx);

  return;
}

// Returns the number of objects passing a given trigger
// Returns zero if the trigger failed
// Returns -1 if the trigger passed but no onjects were found
//--------------------------------------------------------
int nHLTObjects( string arg ){

  // put the trigger name into a string
  TString HLTTrigger( arg );

  // Did the trigger pass?
  if ( !(cms2.passHLTTrigger(HLTTrigger)) ) return 0;

  // The trigger passed, see how many associated objects there are
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "nHLTObjects: Found Trigger: " << arg << endl;
  }
  else {
    cout << "nHLTObjects: Cannot find Trigger " << arg << endl;
    return 0;
  }

  int nobj = 0;
  for( unsigned int i=0; i < cms2.hlt_trigObjs_p4().at(trigIndx).size(); i++ ){
    nobj++;
    //cout << "\t" << i << ", (pt, eta, phi): " << hlt_trigObjs_p4().at(trigIndx).at(i).pt() << " "
    //              << hlt_trigObjs_p4().at(trigIndx).at(i).eta() << " " << hlt_trigObjs_p4().at(trigIndx).at(i).phi() << endl;
  }

  // cout << " Number of jets = " << njets << endl;

  if (nobj == 0) return -1;
  return nobj;
}

LorentzVector p4HLTObject( string arg, int objNumber){

  TString HLTTrigger( arg );
  int trigIndx = -1;
  vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
  vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();
  vector<TString>::const_iterator found_it = find(begin_it, end_it, HLTTrigger );
  if( (found_it != end_it) ){
    trigIndx = found_it - begin_it;
    //cout << "p4HLTObject: Found Trigger: " << arg << endl;
  }
  else {
    cout << "p4HLTObject: Cannot find Trigger: " << arg << endl;
    gSystem->Exit(1);
  }

  int nobj = cms2.hlt_trigObjs_p4().at(trigIndx).size();
  if (nobj == 0 ) {
    cout << "ERROR: nobj == 0" << endl;
    gSystem->Exit(1);
  }

  if (objNumber > (nobj-1)) {
    cout << "ERROR: requested object number " << objNumber << " but we only have " << nobj <<endl;
    gSystem->Exit(1);
  }

  return cms2.hlt_trigObjs_p4().at(trigIndx).at(objNumber);

}

bool passMuTrigger(bool mc, int type) {

  if(mc) {
    if(passHLTTrigger("HLT_Mu9"))
      return true;

    
  } else {

    if(cms2.evt_run() <= 145000) {
      if(nHLTObjects("HLT_Mu9") != 0)
        return true;
    }

    if(cms2.evt_run() > 145000 && cms2.evt_run() <= 147120) {
      if(nHLTObjects("HLT_Mu11") != 0)
        return true;
    }
  
    if(cms2.evt_run() > 147120) {
      if(nHLTObjects("HLT_Mu15_v1") != 0)
        return true;
    }

  }

 // keep bool mc for future data partitions 

 return false;
}

bool passEGTrigger(bool mc, int type) {

  if (mc) {

    int e10 = nHLTObjects("HLT_Ele10_SW_L1R");
    for (int i=0; i<e10; i++) {
      LorentzVector p4 = p4HLTObject("HLT_Ele10_SW_L1R", i);
      if (p4.Pt() > 15.) return true;
    }
 

  } else {  // data now

    if (cms2.evt_run() < 138000) {
      int e10 = nHLTObjects("HLT_Ele10_LW_L1R");
      for (int i=0; i<e10; i++) {
        LorentzVector p4 = p4HLTObject("HLT_Ele10_LW_L1R", i);
        if (p4.Pt() > 15.) return true;
      }
    }

    if (cms2.evt_run() >= 138000 && cms2.evt_run() < 141900) {
      int e15 = nHLTObjects("HLT_Ele15_LW_L1R");
      if (e15 != 0) 
        return true;
    }

    if (cms2.evt_run() >= 141900 && cms2.evt_run() <= 144000) {
      int e15 = nHLTObjects("HLT_Ele15_SW_L1R");
      if (e15 != 0) 
        return true;
    }   
    
    if (cms2.evt_run() > 144000 && cms2.evt_run() <= 144114) {

      int e15caloId = nHLTObjects("HLT_Ele15_SW_CaloEleId_L1R");
      if (e15caloId != 0)
        return true;
      
      int e20 = nHLTObjects("HLT_Ele20_SW_L1R");
      if (e20 != 0)
        return true;
      
      int ed10 = nHLTObjects("HLT_DoubleEle10_SW_L1R");
      if(ed10 != 0)
        return true;
    }

    if(cms2.evt_run() > 146000 && cms2.evt_run() <= 147120) {

      if(nHLTObjects("HLT_DoubleEle10_SW_L1R") != 0 )
        return true;
      
      if(nHLTObjects("HLT_Ele17_SW_CaloEleId_L1R") != 0)
        return true;
    }

    if(cms2.evt_run() > 147120 && cms2.evt_run() <= 148100 ) {
    
      if(nHLTObjects("HLT_DoubleEle15_SW_L1R_v1") != 0)
        return true;

      if(nHLTObjects("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1") != 0)
        return true;
      
      if(nHLTObjects("HLT_Ele17_SW_TightEleId_L1R") != 0)
        return true;

    }
    if(cms2.evt_run()  > 148100) {

      vector<TString>::const_iterator begin_it = cms2.hlt_trigNames().begin();
      vector<TString>::const_iterator end_it = cms2.hlt_trigNames().end();

      std::string HLTTrigger = "HLT_DoubleEle17_SW_L1R_v1";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }
      
      HLTTrigger = "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }

      HLTTrigger = "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }

      
      HLTTrigger = "HLT_Ele22_SW_TighterEleId_L1R_v3";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }


      HLTTrigger = "HLT_Ele22_SW_TighterEleId_L1R_v2";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
          return true;
      }

      
      HLTTrigger = "HLT_Ele17_SW_TighterEleIdIsol_L1R_v3";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }


      HLTTrigger = "HLT_Ele17_SW_TighterEleIdIsol_L1R_v2";
      if(find(begin_it, end_it, HLTTrigger ) != end_it )  {
        if(nHLTObjects(HLTTrigger) != 0)
        return true;
      }

    }
  }
  return false;
}

unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps) {
  
  float maxSumPt = 0;
  unsigned int bestHypIdx = 0;
  for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
    
    unsigned int index = v_goodHyps.at(i);
    float sumPt = hyp_lt_p4()[index].Pt() + hyp_ll_p4()[index].Pt();
    if( sumPt > maxSumPt) {
      maxSumPt = sumPt;
      bestHypIdx = index;
    }
  }

  return bestHypIdx;

}

int getNbtags(const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator) {

  if(jetAlgo != "jptJets" && jetAlgo != "caloJets" && jetAlgo != "pfJets") {
    cout << "Unknown jet algorithm. Returning spurious value" << endl;
    return -9999;
  }

  if(bTagDiscriminator != "trackCountingHighEffBJetTag" &&
     bTagDiscriminator != "simpleSecondaryVertexHighEffBJetTag" && 
     bTagDiscriminator != "simpleSecondaryVertexHighPurBJetTag") {
    cout << "Unknown bTag Discriminator. Returning spurious value" << endl;
    return -9999;
  }

  if(jetAlgo == "jptJets") {
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jpts_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
          ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jpts_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
          ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jpts_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
          ntags++;
      }
      return ntags;
    }
  }


  if(jetAlgo == "caloJets") {           
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
          ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
          ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(jets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
          ntags++;
      }
      return ntags;
    }
  }
  
  if(jetAlgo == "pfJets") {             
    int ntags = 0;
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(pfjets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7)
          ntags++;
      }
      return ntags;
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(pfjets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74)
          ntags++;
      }
      return ntags;
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        if(pfjets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2)
          ntags++;
      }
      return ntags;
    }
  }
 return -9999;
}

std::pair<int,float> getMinbltags(unsigned int hypIdx, const vector<LorentzVector> v_jetP4s, const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator) {

   LorentzVector ll_p4 = hyp_ll_p4()[hypIdx];
   LorentzVector lt_p4 = hyp_lt_p4()[hypIdx];
   vector<LorentzVector> v_bjets_p4;
   LorentzVector llbjet; 
   LorentzVector ltbjet; 

  if(jetAlgo != "jptJets" && jetAlgo != "caloJets" && jetAlgo != "pfJets") {
    cout << "Unknown jet algorithm. Returning spurious value" << endl;
    return make_pair(-999, -999.);
  }

  if(bTagDiscriminator != "trackCountingHighEffBJetTag" &&
     bTagDiscriminator != "simpleSecondaryVertexHighEffBJetTag" &&
     bTagDiscriminator != "simpleSecondaryVertexHighPurBJetTag") {
    cout << "Unknown bTag Discriminator. Returning spurious value" << endl;
    return make_pair(-999, -999.);
  }

  if(jetAlgo == "jptJets") {
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jpts_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7) v_bjets_p4.push_back(p4);
      }
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jpts_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74) v_bjets_p4.push_back(p4);
      }
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jpts_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2) v_bjets_p4.push_back(p4);
      }
    }
  }


  if(jetAlgo == "caloJets") {
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7) v_bjets_p4.push_back(p4);
      }
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74) v_bjets_p4.push_back(p4);
      }
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(jets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2) v_bjets_p4.push_back(p4);
      }
    }
  }
  if(jetAlgo == "pfJets") {
    if(bTagDiscriminator == "trackCountingHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) { 
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(pfjets_trackCountingHighEffBJetTag()[v_jetsIdx.at(i)] > 1.7) v_bjets_p4.push_back(p4);
      }
    }else if(bTagDiscriminator == "simpleSecondaryVertexHighEffBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(pfjets_simpleSecondaryVertexHighEffBJetTag()[v_jetsIdx.at(i)] > 1.74) v_bjets_p4.push_back(p4);
      }
    } else if(bTagDiscriminator == "simpleSecondaryVertexHighPurBJetTag") {
      for(unsigned int i = 0; i < v_jetsIdx.size(); i++) {
        LorentzVector p4 = v_jetP4s[v_jetsIdx.at(i)];
        if(pfjets_simpleSecondaryVertexHighPurBJetTags()[v_jetsIdx.at(i)] > 2) v_bjets_p4.push_back(p4);
      }
    }
  }
// Sort the jets
 std::sort(v_bjets_p4.begin(), v_bjets_p4.end(), sortLVByPt);

 if (v_bjets_p4.size() > 0) {
     llbjet = v_bjets_p4[0] + ll_p4;
     ltbjet = v_bjets_p4[0] + lt_p4;
     return make_pair(v_bjets_p4.size(), min(llbjet.mass(), ltbjet.mass()));
  } else {
  return make_pair(0, -999.);
 }
  return make_pair(0, -999.);
}





