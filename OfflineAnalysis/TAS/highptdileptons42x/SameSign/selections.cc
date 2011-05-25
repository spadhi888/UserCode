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
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);//, bool used0wrtPV = false);
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);//, bool used0wrtPV  = false);
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);//, bool used0wrtPV = false);
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);//, bool used0wrtPV = false);
bool isGoodDilHypJet(LorentzVector jetp4, unsigned int& hypIdx, double ptCut, double absEtaCut, double dRCut, bool muJetClean);
bool isGoodJet(LorentzVector jetp4, double ptCut, double absEtaCut, double dRCut, bool muJetClean, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int hypIdx);
// std::pair<float,float> getMet(string& algo, unsigned int hypIdx, std::string prefix);
std::pair<float,float> getMet(const string algo, unsigned int hypIdx);
bool inZmassWindow (float mass);
bool passTriggersMu9orLisoE15(int dilType);
bool isFakeableMuon(int index, int vtx2011);
double getd0wrtPV(LorentzVector p4, float d0);
void correctTcMETForHypMus(unsigned int hypIdx, double& met, double& metPhi);
bool isFakeableElectron(int index, string prefix, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);
bool additionalZvetoSUSY2010 (int i_hyp, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011);
bool passEGTrigger(unsigned int hypIdx, bool mc);
bool passMuTrigger(unsigned int hypIdx);
int nHLTObjects(string arg);
int getNbtags(const vector<unsigned int> v_jetsIdx, const string jetAlgo, const string bTagDiscriminator);
unsigned int selectHypByHighestSumPt(const vector<unsigned int> &v_goodHyps);

LorentzVector p4HLTObject(string arg, int) ;


bool isFakeableElectron (int index, string prefix, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) {
  TPMERegexp re1("V1", "g");
  TPMERegexp re2("V2", "g");
  TPMERegexp re3("V3", "g");

  if (re1.Match(prefix)) return pass_electronSelection(index, electronSelectionFOV3_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011);
  if (re2.Match(prefix)) return pass_electronSelection(index, electronSelectionFOV3_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011);
  if (re3.Match(prefix)) return pass_electronSelection(index, electronSelectionFOV3_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011);

  return false;
}

bool isFakeableMuon (int index, int vtx2011) {
     return muonId(index, muonSelectionFO_ssV3);
}


/******************************************************************************************/     
// good lepton (either mu or electron, no isolation cuts)
/******************************************************************************************/
bool isGoodLeptonNoIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) {//, bool used0wrtPV) {

  if(abs(id) == 11) {
    return (pass_electronSelection(lepIdx, electronSelection_ssV3_noIso, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011));
  }


  if(abs(id) == 13) {
    return muonId(lepIdx, NominalSSv3, vtx2011);
  }

  return true;
}

/******************************************************************************************/
// isolated lepton (either mu or electron)
/******************************************************************************************/
bool isGoodLeptonwIso(int id, int lepIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) { //, bool used0wrtPV) {

// Covers both the noiso part

  if(!isGoodLeptonNoIso(id, lepIdx, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) return false;

  if(abs(id)== 11) {
       if (!pass_electronSelection(lepIdx, electronSelection_ssV3_iso, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) return false;
  }

  // 13 is a muon
  if(abs(id) == 13)
    if(muonIsoValue(lepIdx, false) > 0.15)   return false;
  return true;
}


bool additionalZvetoSUSY2010(int i_hyp, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) {
  bool veto=false;

  // first, look for Z->mumu
  for (unsigned int i=0; i < mus_p4().size(); i++) {
    bool hypLep1 = false;
    if (mus_p4().at(i).pt() < 5.)     continue;

    if (!isGoodLeptonNoIso(13, i, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;

    if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 13 && hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 13 && hyp_ll_index()[i_hyp] == i ) hypLep1 = true;
    
    for (unsigned int j=i+1; j < mus_p4().size(); j++) {
      bool hypLep2 = false;
      if (mus_p4().at(j).pt() < 5.)     continue;

      if (!isGoodLeptonNoIso(13, j, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;

      if (mus_charge().at(i) == mus_charge().at(j)) continue;
      if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 13 && hyp_lt_index()[i_hyp] == j ) hypLep2 = true;
      if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 13 && hyp_ll_index()[i_hyp] == j ) hypLep2 = true;
      if ((muonIsoValue(i, false) > 0.15) && (muonIsoValue(j, false) > 0.15)) continue;
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

    if (!isGoodLeptonNoIso(11, i, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;

    if ( TMath::Abs(hyp_lt_id()[i_hyp]) == 11 && hyp_lt_index()[i_hyp] == i ) hypLep1 = true;
    if ( TMath::Abs(hyp_ll_id()[i_hyp]) == 11 && hyp_ll_index()[i_hyp] == i ) hypLep1 = true;

    for (unsigned int j=i+1; j<els_p4().size(); j++) {
      bool hypLep2 = false;
      if (els_p4().at(j).pt() < 10.) continue;

      if (!isGoodLeptonNoIso(11, j, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;
      if (els_charge().at(i) == els_charge().at(j)) continue;

      if ((!pass_electronSelection(i, electronSelection_ssV3_iso, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) && (!pass_electronSelection(j, electronSelection_ssV3_iso, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011))) continue;
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
bool isGoodHypNoIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) {//, bool used0wrtPV) {
  
  if(!isGoodLeptonNoIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011))//, used0wrtPV)
     return false;
  if(!isGoodLeptonNoIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011))//, used0wrtPV)
    return false;

  return true;
}

/******************************************************************************************/     
// are the leptons in the hypothesis isolated?
/******************************************************************************************/     
bool isGoodHypwIso(int hypIdx, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int vtx2011) {//, bool used0wrtPV) {


  if(!isGoodLeptonwIso(hyp_lt_id()[hypIdx], hyp_lt_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011))//, used0wrtPV)
    return false;
  if(!isGoodLeptonwIso(hyp_ll_id()[hypIdx], hyp_ll_index()[hypIdx], applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011))//, used0wrtPV)
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

bool isGoodJet(LorentzVector jetp4, double ptCut, double absEtaCut, double dRCut, bool muJetClean, bool applyAlignmentCorrection, bool removedEtaCutInEndcap, int hypIdx){


  if(jetp4.Pt() < ptCut) return false;
  if(fabs(jetp4.Eta()) > absEtaCut) return false;

  for (unsigned int i=0; i < els_p4().size(); i++) {
    if (els_p4().at(i).pt() < 10.)  continue;
    int vtx2011 = hypsFromSameVtx2011Int(hypIdx, 1.0, true , false);
    if (!isGoodLeptonwIso(11, i, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;
    double dR_ell = ROOT::Math::VectorUtil::DeltaR(els_p4().at(i),jetp4);
    if (dR_ell < dRCut) return false;
  }


  if (muJetClean){
      for (unsigned int i=0; i < mus_p4().size(); i++) {
        if (mus_p4().at(i).pt() < 10.)     continue;
        int vtx2011 = hypsFromSameVtx2011Int(hypIdx, 1.0, true , false);
        if (!isGoodLeptonwIso(13, i, applyAlignmentCorrection, removedEtaCutInEndcap, vtx2011)) continue;
        double dR_mll = ROOT::Math::VectorUtil::DeltaR(mus_p4().at(i),jetp4);
        if (dR_mll < dRCut) return false;
    }
  }
  return true;
}

std::pair<float,float> getMet(const string algo, unsigned int hypIdx) {
  
  if(algo != "tcMET" && algo != "muCorMET" && algo != "pfMET" && algo != "tcMET35X") {
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

bool passMuTrigger(unsigned int hypIdx) {

         if (passHLTTrigger("HLT_Mu9"))
                  return true;

         if (passHLTTrigger("HLT_Mu15_v1"))
                  return true;

         if (cms2.hyp_type()[hypIdx] == 0)
         {
                  if (passHLTTrigger("HLT_DoubleMu3"))
                           return true;
                  if (passHLTTrigger("HLT_DoubleMu5_v1"))
                           return true;
         }

         if (cms2.hyp_type()[hypIdx] == 1 || cms2.hyp_type()[hypIdx] == 2)
         {
                  if (passHLTTrigger("HLT_Mu5_Ele9_v1"))
                           return true;
                  if (passHLTTrigger("HLT_Mu8_Ele8_v1"))
                           return true;
         }
 // keep bool mc for future data partitions 

 return false;
}

bool passEGTrigger(unsigned int hypIdx, bool mc) {

         unsigned int hypType = cms2.hyp_type()[hypIdx];

         if (mc)
         {
                  if (passHLTTrigger("HLT_Ele10_LW_L1R"))
                           return true;

                  if (passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
                           return true;

                  if (passHLTTrigger("HLT_Ele15_LW_L1R"))
                           return true;

                  if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
                           return true;
         }
         else // data now
         {
                  if(evt_run() < 138000)
                  {
                           if (passHLTTrigger("HLT_Ele10_LW_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele15_LW_L1R"))
                                        return true;

                           if (hypType == 3)
                                        if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
                                                 return true;
                  }

                  if(evt_run() >= 138000 && evt_run() < 141900)
                  {
                           if (passHLTTrigger("HLT_Ele10_LW_EleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele15_LW_L1R"))
                                        return true;

                           if (hypType == 3)
                                        if (passHLTTrigger("HLT_DoubleEle5_SW_L1R"))
                                                 return true;
                  }
                  if (cms2.evt_run() >= 141900)
                  {
                           if (passHLTTrigger("HLT_Ele10_SW_EleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele15_SW_CaloEleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele15_SW_EleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele17_SW_LooseEleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele17_SW_CaloEleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele17_SW_EleId_L1R"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1"))
                                        return true;

                           if (passHLTTrigger("HLT_Ele17_SW_TightEleIdIsol_L1R_v1"))
                                        return true;

                           if (hypType == 3)
                           {
                                        if (passHLTTrigger("HLT_DoubleEle10_SW_L1R"))
                                                 return true;

                                        if (passHLTTrigger("HLT_DoubleEle15_SW_L1R_v1"))
                                                 return true;
                           }
                  }
         }


  return false;
}

/*****************************************************************************************/
//get the number of jets passing btag discriminator cuts
// takes as arguments a vector of indices, the jet algorithm
// and the btag discriminator
// the working points are hard coded
/*****************************************************************************************/
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


