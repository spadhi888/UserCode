#include <assert.h>
#include <algorithm>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "TPRegexp.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/muonSelections.cc"
#include "../CORE/metSelections.cc"
#include "CMS2.cc"

using namespace tas;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
bool isGoodLeptonNoIso(int id, int lepIdx);//, bool used0wrtPV = false);
bool isGoodLeptonwIso(int id, int lepIdx);//, bool used0wrtPV  = false);
bool isGoodHypNoIso(int hypIdx);//, bool used0wrtPV = false);
bool isGoodHypwIso(int hypIdx);//, bool used0wrtPV = false);
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);
std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix);
bool inZmassWindow (float mass);
bool passTriggersMu9orLisoE15(int dilType);
int eventDilIndexByWeightTTDil08(const std::vector<unsigned int>& goodHyps, int& strasbourgDilType, bool printDebug, bool usePtOnlyForWeighting);
bool isFakeDenominatorElectron_v1(unsigned int lepIdx);
bool isFakeDenominatorElectron_v2(unsigned int lepIdx);
bool isFakeableMuon(int index);
double getd0wrtPV(LorentzVector p4, float d0);
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi);
bool isFakeableElectron(int index, string prefix);


bool isFakeableElectron (int index, string prefix) {
  TPMERegexp re1("v1", "g");
  TPMERegexp re2("v2", "g");
  TPMERegexp re3("v3", "g");

  if (re1.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v1);
  if (re2.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v2);
  if (re3.Match(prefix)) return pass_electronSelection(index, electronSelectionFO_el_ttbarV1_v3);

  return false;
}

bool isFakeableMuon (int index) {
     return muonId(index, muonSelectionFO_mu_ttbar);
}


bool isGoodLeptonLooseID(int id, int lepIdx) {

  if(abs(id) == 11) {
//    if (!pass_electronSelection(lepIdx, electronSelection_ttbarV2_noiso)) return false;
//              cuts_t cuts_passedele = electronSelection(lepIdx);
//              std::cout << bool((cuts_passedele & (1ll<<ELEISO_REL010))) << " ";
//              std::cout << bool((cuts_passedele & (1ll<<ELEID_CAND01))) << " ";
//              std::cout << bool((cuts_passedele & (1ll<<ELEIP_400))) << " ";
//              std::cout << bool((cuts_passedele & (1ll<<ELESEED_ECAL))) << " ";
//              std::cout << bool((cuts_passedele & (1ll<<ELEETA_250))) << " ";
//             std::cout << bool((cuts_passedele & (1ll<<ELENOTCONV_DISTDCOT002))) << " ";
//              std::cout << bool((cuts_passedele & (1ll<<ELENOMUON_010))) << std::endl;
//    if (!pass_electronSelection(lepIdx, electronSelection_ttbarV2_noiso)) return false;
//              std::cout << "Passed Sanjay " << std::endl;
  }

  if(abs(id) == 13) {
    if ( TMath::Abs(mus_p4()[lepIdx].eta()) > 2.5)  return false; // eta cut
    if (mus_gfit_chi2().at(lepIdx)/mus_gfit_ndof().at(lepIdx) >= 10) return false; //glb fit chisq
    if (((mus_type().at(lepIdx)) & (1<<1)) == 0)    return false; // global muon
    if (((mus_type().at(lepIdx)) & (1<<2)) == 0)    return false; // tracker muon
    if (mus_validHits().at(lepIdx) < 11)            return false; // # of tracker hits
    if (cms2.mus_gfit_validSTAHits().at(lepIdx) == 0)    return false; // Glb fit must have hits in mu chambers

  }
  return true;
}

/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx) {//, bool used0wrtPV) {

  if(abs(id) == 11) {
    const cuts_t elIDcuts =   
	  (1ll<<ELEID_VBTF_35X_90) |
	  (1ll<<ELEIP_400) |
	  (1ll<<ELENOMUON_010) |
	  (1ll<<ELENOTCONV_HITPATTERN) |
	  (1ll<<ELENOTCONV_DISTDCOT002) |
	  (1ll<<ELESCET_010) |
	  (1ll<<ELEPT_010) |
	  (1ll<<ELEETA_250) |
	  (1ll<<ELESEED_ECAL) |
	  (1ll<<ELENOSPIKE_SWISS005);


	//    bool isSpike = isSpikeElectron(lepIdx);
	//    unsigned int answer_vbtf = electronId_VBTF(lepIdx, VBTF_35X_90);
	//    bool elsvbtf90_ = ( ( answer_vbtf & (1ll<<ELEID_ID) ) == (1ll<<ELEID_ID) );
    return (pass_electronSelection(lepIdx, elIDcuts));

  }


  if(abs(id) == 13) {
    if ( mus_p4()[lepIdx].pt() < 5.) {
      std::cout << "muonID ERROR: requested muon is too low pt,  Abort." << std::endl;
      return false;
    }
    if ( TMath::Abs(mus_p4()[lepIdx].eta()) > 2.5)  return false; // eta cut
    if (mus_gfit_chi2().at(lepIdx)/mus_gfit_ndof().at(lepIdx) >= 10) return false; //glb fit chisq
    if (((mus_type().at(lepIdx)) & (1<<1)) == 0)    return false; // global muon
    if (((mus_type().at(lepIdx)) & (1<<2)) == 0)    return false; // tracker muon
    if (mus_validHits().at(lepIdx) < 11)            return false; // # of tracker hits
    //if(used0wrtPV) {
    //if(fabs(getd0wrtPV(mus_p4()[lepIdx], mus_d0()[lepIdx])) > 0.04)   return false;
    //return false;
    //} else 
    if (TMath::Abs(mus_d0corr().at(lepIdx)) > 0.02) return false; // d0 from beamspot
    if (cms2.mus_gfit_validSTAHits().at(lepIdx) == 0)    return false; // Glb fit must have hits in mu chambers
  }
  

  return true;
}

/******************************************************************************************/     
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx) { //, bool used0wrtPV) {

// Covers both the noiso part
 
  if(!isGoodLeptonNoIso(id, lepIdx))
    return false;

  // 11 is a electron
  if(abs(id)== 11) {
       const cuts_t elISOcuts =   (1ll<<ELEISO_REL015);
       if (!pass_electronSelection(lepIdx, elISOcuts))
            return false;
  }

  // 13 is a muon
  if(abs(id) == 13)
    if(muonIsoValue(lepIdx) > 0.15)   return false;

  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis good (all cuts but isolation?)
/******************************************************************************************/
bool isGoodHypNoIso(int hypIdx) {//, bool used0wrtPV) {
  
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx]))//, used0wrtPV)
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx]))//, used0wrtPV)
    return false;

  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx) {//, bool used0wrtPV) {


  if(cms2.hyp_type()[hypIdx] == 3) {
    /*
    cout << "evt_event: " << evt_event() << endl;
    cout << "lt, ll pt: " << hyp_lt_p4()[hypIdx].Pt() << "," << hyp_ll_p4()[hypIdx].Pt() << endl;
    cout << "lt, ll eta: " << hyp_lt_p4()[hypIdx].Eta() << "," << hyp_ll_p4()[hypIdx].Eta() << endl;
    cout << "lt iso: " << electronIsolation_relsusy_cand1(hyp_lt_index()[hypIdx], true) << endl;
    cout << "ll iso: " << electronIsolation_relsusy_cand1(hyp_ll_index()[hypIdx], true) << endl;
    */
  }
  
  
  if(!isGoodLeptonwIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx]))//, used0wrtPV)
    return false;
  if(!isGoodLeptonwIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx]))//, used0wrtPV)
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

  
  //the cut is 30 for ee/mm hyps to reject DY
  //20 for emu
  //if(hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3) {
  //if(met < 30.) 
  //return false;
  //  }
  //if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
  //    if(met < 20.)
  //    return false;
  //}
 
  //return true;
  
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
//electron FO v1
/******************************************************************************************/     
bool isFakeDenominatorElectron_v1(unsigned int lepIdx) {
//  if (fabs(els_p4()[lepIdx].Eta()) > 2.5)    return false;
//  if (els_p4()[lepIdx].Pt() < 20.)           return false;
//  if (!electronId_noMuon(lepIdx))            return false;
//  if (isFromConversionPartnerTrack(lepIdx))  return false;
// //  if (electronIsolation_relsusy_cand1(lepIdx, true) > 0.40) return false;

  return true;
  
}

/******************************************************************************************/     
//electron FO v1
/******************************************************************************************/     
bool isFakeDenominatorElectron_v2(unsigned int lepIdx) {

//  if (fabs(els_p4()[lepIdx].Eta()) > 2.5)    return false;
//  if (els_p4()[lepIdx].Pt() < 20.)           return false;
//  if (!electronId_noMuon(lepIdx))            return false;
//  if (isFromConversionPartnerTrack(lepIdx))  return false;
//  if (electronIsolation_relsusy_cand1(lepIdx, true) > 0.10) return false;

  return true;

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
