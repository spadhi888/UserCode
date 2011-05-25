#include "babymakercommon.h"
#include "trilepbabymaker.h" 

#include <algorithm>
#include <iostream>

#include "TChain.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"

#include "CORE/CMS2.h"
#include "CORE/mcSelections.h"
#include "CORE/electronSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/metSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/trackSelections.h"

void trilepbabymaker::ScanChain (const char *inputFilename, const char *babyFilename, int nEvents)
{
    TChain *chain = new TChain("Events");
    chain->Add(inputFilename);
    TObjArray *listOfFiles = chain->GetListOfFiles();    

    unsigned int nEventsChain=0;
    if (nEvents==-1) 
        nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;

    // make a baby ntuple
    MakeBabyNtuple(babyFilename);

    // file loop
    TIter fileIter(listOfFiles);
    TFile *currentFile = 0;
    while ((currentFile = (TFile*)fileIter.Next()))
    {
        TFile f(currentFile->GetTitle());
        TTree *tree = (TTree*)f.Get("Events");
        cms2.Init(tree);

        //Event Loop
        unsigned int nEvents = tree->GetEntries();
        for(unsigned int event = 0; event < nEvents; ++event)
        {
            cms2.GetEntry(event);
            ++nEventsTotal;

            // Progress feedback to the user
            CMS2::progress(nEventsTotal, nEventsChain );

            // trilepton hypothesis stuff
            for (unsigned hypi = 0; hypi < cms2.hyp_trilep_bucket().size(); ++hypi)
            {
                int index1 = cms2.hyp_trilep_first_index()[hypi];
                int index2 = cms2.hyp_trilep_second_index()[hypi];
                int index3 = cms2.hyp_trilep_third_index()[hypi];

                int type1 = abs(cms2.hyp_trilep_first_type()[hypi]);
                int type2 = abs(cms2.hyp_trilep_second_type()[hypi]);
                int type3 = abs(cms2.hyp_trilep_third_type()[hypi]);

                // require any muons to be global OR tracker (i.e. not SA only)
                if (type1 == 1 && !(cms2.mus_type()[index1] & 6))
                    continue;
                if (type2 == 1 && !(cms2.mus_type()[index2] & 6))
                    continue;
                if (type3 == 1 && !(cms2.mus_type()[index3] & 6))
                    continue;

                // initialize baby quantities
                InitBabyNtuple();

                // event stuff
                strcpy(dataset_, cms2.evt_dataset().Data());
                run_        = cms2.evt_run();
                ls_         = cms2.evt_lumiBlock();
                evt_        = cms2.evt_event();
                isdata_     = cms2.evt_isRealData();
                nvtx_ = 0;
                for( unsigned int i=0; i<cms2.vtxs_isFake().size(); i++ ) {
                    if( !cms2.vtxs_isFake()[i] && cms2.vtxs_isValid()[i] )
                        ++nvtx_;
                }


                if (!isdata_) {
                    int nlep  = leptonGenpCount_lepTauDecays(ngenels_, ngenmus_, ngentaus_);
                    scale1fb_ = cms2.evt_scale1fb();
                    pthat_    = cms2.genps_pthat();
                }

                pfmet_      = cms2.evt_pfmet();
                tcmet_      = cms2.evt_pf_tcmet();
                calotcmet_  = cms2.evt_tcmet();
                ntrks_      = cms2.trks_trk_p4().size();

                float thePFMetPhi     = cms2.evt_pfmetPhi();
                float theTCMetPhi     = cms2.evt_pf_tcmetPhi();
                float theCaloTCMetPhi = cms2.evt_tcmetPhi();

                float metx  = tcmet_ * cos(theTCMetPhi);
                float mety  = tcmet_ * sin(theTCMetPhi);
                float cmetx = calotcmet_ * cos(theCaloTCMetPhi);
                float cmety = calotcmet_ * sin(theCaloTCMetPhi);
                for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                    if (!wasMetCorrectedForThisMuon(muj, usingTcMet) && muonIdNotIsolated(muj, NominalTTbarV2)) {
                        fixMetForThisMuon(muj, metx, mety, usingTcMet);
                        fixMetForThisMuon(muj, cmetx, cmety, usingTcMet);
                    }
                }
                tcmet_ = sqrt(metx * metx + mety * mety);
                theTCMetPhi = atan2(mety, metx);
                calotcmet_ = sqrt(cmetx * cmetx + cmety * cmety);
                theCaloTCMetPhi = atan2(cmety, cmetx);

                // loop over muons and electrons to get ngoodlep
                ngoodlep_ = 0;
                ngoodmus_ = 0;
                ngoodels_ = 0;
                for (unsigned muii = 0; muii < cms2.mus_p4().size(); ++muii) {
                    if (cms2.mus_p4()[muii].pt() > 20. && muonId(muii, NominalTTbarV2)) {
                        ++ngoodlep_;
                        ++ngoodmus_;
                    }
                }
                for (unsigned eli = 0; eli < cms2.els_p4().size(); ++eli) {
                    if (cms2.els_p4()[eli].pt() > 20. && pass_electronSelection(eli, electronSelection_ttbarV2)) {
                        ++ngoodlep_;
                        ++ngoodels_;
                    }
                }

                // initialize meff to 0
                pfmeff_ = 0.;
                tcmeff_ = 0.;

                hyp_type_    = cms2.hyp_trilep_bucket()[hypi];
                pt1_         = type1 == 1 ? cms2.mus_p4()[index1].pt() : cms2.els_p4()[index1].pt();
                pt2_         = type2 == 1 ? cms2.mus_p4()[index2].pt() : cms2.els_p4()[index2].pt();
                pt3_         = type3 == 1 ? cms2.mus_p4()[index3].pt() : cms2.els_p4()[index3].pt();
                eta1_        = type1 == 1 ? cms2.mus_p4()[index1].eta() : cms2.els_p4()[index1].eta();
                eta2_        = type2 == 1 ? cms2.mus_p4()[index2].eta() : cms2.els_p4()[index2].eta();
                eta3_        = type3 == 1 ? cms2.mus_p4()[index3].eta() : cms2.els_p4()[index3].eta();
                phi1_        = type1 == 1 ? cms2.mus_p4()[index1].phi() : cms2.els_p4()[index1].phi();
                phi2_        = type2 == 1 ? cms2.mus_p4()[index2].phi() : cms2.els_p4()[index2].phi();
                phi3_        = type3 == 1 ? cms2.mus_p4()[index3].phi() : cms2.els_p4()[index3].phi();
                d0corr1_     = type1 == 1 ? cms2.mus_d0corr()[index1] : cms2.els_d0corr()[index1];
                d0corr2_     = type2 == 1 ? cms2.mus_d0corr()[index2] : cms2.els_d0corr()[index2];
                d0corr3_     = type3 == 1 ? cms2.mus_d0corr()[index3] : cms2.els_d0corr()[index3];
                eormu1_      = type1 == 1 ? -13 * cms2.mus_charge()[index1] : -11 * cms2.els_charge()[index1];
                eormu2_      = type2 == 1 ? -13 * cms2.mus_charge()[index2] : -11 * cms2.els_charge()[index2];
                eormu3_      = type3 == 1 ? -13 * cms2.mus_charge()[index3] : -11 * cms2.els_charge()[index3];
                dphipfmet1_  = type1 == 1 ? deltaPhi(cms2.mus_p4()[index1].phi(), thePFMetPhi) : deltaPhi(cms2.els_p4()[index1].phi(), thePFMetPhi);
                dphipfmet2_  = type2 == 1 ? deltaPhi(cms2.mus_p4()[index2].phi(), thePFMetPhi) : deltaPhi(cms2.els_p4()[index2].phi(), thePFMetPhi);
                dphipfmet3_  = type3 == 1 ? deltaPhi(cms2.mus_p4()[index3].phi(), thePFMetPhi) : deltaPhi(cms2.els_p4()[index3].phi(), thePFMetPhi);
                dphitcmet1_  = type1 == 1 ? deltaPhi(cms2.mus_p4()[index1].phi(), theTCMetPhi) : deltaPhi(cms2.els_p4()[index1].phi(), theTCMetPhi);
                dphitcmet2_  = type2 == 1 ? deltaPhi(cms2.mus_p4()[index2].phi(), theTCMetPhi) : deltaPhi(cms2.els_p4()[index2].phi(), theTCMetPhi);
                dphitcmet3_  = type3 == 1 ? deltaPhi(cms2.mus_p4()[index3].phi(), theTCMetPhi) : deltaPhi(cms2.els_p4()[index3].phi(), theTCMetPhi);

                // add MET to meff
                pfmeff_ += pfmet_;
                tcmeff_ += tcmet_;

                // add lepton pt to meff
                pfmeff_ += (pt1_ + pt2_ + pt3_);
                tcmeff_ += (pt1_ + pt2_ + pt3_);					

                // clean jets for hyp leptons
                VofP4 theJets;
                std::vector<unsigned int> theJetIndices;
                sumjetpt_ = 0.0;
                for (unsigned int jeti = 0; jeti < cms2.pfjets_p4().size(); ++jeti)
                {
                    LorentzVector vjet = cms2.pfjets_p4()[jeti];
                    bool jetIsLep = false;
                    for (unsigned int muj = 0; muj < cms2.mus_p4().size(); ++muj) {
                        LorentzVector vlep = cms2.mus_p4()[muj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && cms2.mus_p4()[muj].pt() > 20. && muonId(muj, NominalTTbarV2))
                            jetIsLep = true;
                    }
                    for (unsigned int elj = 0; elj < cms2.els_p4().size(); ++elj) {
                        LorentzVector vlep = cms2.els_p4()[elj];
                        if (dRbetweenVectors(vjet, vlep) < 0.4 && cms2.els_p4()[elj].pt() > 20. && pass_electronSelection(elj, electronSelection_ttbarV2))
                            jetIsLep = true;
                    }
                    if (jetIsLep) continue;

                    if (cms2.pfjets_p4()[jeti].pt() > 30. && fabs(cms2.pfjets_p4()[jeti].eta()) < 2.5 && isGoodPFJet(jeti)) {
                        theJets.push_back(cms2.pfjets_p4()[jeti]);
                        theJetIndices.push_back(jeti);
                        sumjetpt_ += vjet.Pt();
                    }
                }
                std::sort(theJets.begin(), theJets.end(), sortByPt);
                std::sort(theJetIndices.begin(), theJetIndices.end(), sortByPFJetPt);

                njets_        = theJetIndices.size();
                jet1pt_       = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].pt()  : -999999.;
                jet1eta_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].eta() : -999999.;
                jet1phi_      = theJetIndices.size() > 0 ? cms2.pfjets_p4()[theJetIndices[0]].phi() : -999999.;
                jet2pt_       = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].pt()  : -999999.;
                jet2eta_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].eta() : -999999.;
                jet2phi_      = theJetIndices.size() > 1 ? cms2.pfjets_p4()[theJetIndices[1]].phi() : -999999.;
                jet3pt_       = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].pt()  : -999999.;
                jet3eta_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].eta() : -999999.;
                jet3phi_      = theJetIndices.size() > 2 ? cms2.pfjets_p4()[theJetIndices[2]].phi() : -999999.;

                LorentzVector dijetP4;
                jetmass_ = theJetIndices.size() > 1 ? sqrt((cms2.pfjets_p4()[theJetIndices[0]]+cms2.pfjets_p4()[theJetIndices[1]]).M2()) : -999999.; 

                double mindphipfmet = 999999.;
                double mindphitcmet = 999999.;
                neffbtags_   = 0;
                npurbtags_   = 0;
                ntceffbtags_ = 0;
                ntcpurbtags_ = 0;
                jet1isBtag_  = 0;
                jet2isBtag_  = 0;
                jet3isBtag_  = 0;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    if (cms2.pfjets_simpleSecondaryVertexHighEffBJetTag()[theJetIndices[jeti]] > 1.74)
                    {
                        ++neffbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
                    if (cms2.pfjets_simpleSecondaryVertexHighPurBJetTags()[theJetIndices[jeti]] > 2.)
                    {
                        ++npurbtags_;

                        if (jeti == 0)
                            jet1isBtag_ = 1;
                        else if (jeti == 1)
                            jet2isBtag_ = 1;
                        else if (jeti == 2)
                            jet3isBtag_ = 1;
                    }
                    if (cms2.pfjets_trackCountingHighEffBJetTag()[jeti] > 1.7)
                        ++ntceffbtags_;
                    if (cms2.pfjets_trackCountingHighPurBJetTag()[jeti] > 1.19)
                        ++ntcpurbtags_;

                    float currdphipfmet = deltaPhi(thePFMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphipfmet < mindphipfmet)
                        mindphipfmet = currdphipfmet;

                    float currdphitcmet = deltaPhi(theTCMetPhi, cms2.pfjets_p4()[theJetIndices[jeti]].phi());
                    if (currdphitcmet < mindphitcmet)
                        mindphitcmet = currdphitcmet;

                    // add jet pt to meff
                    pfmeff_ += cms2.pfjets_p4()[theJetIndices[jeti]].pt();
                    tcmeff_ += cms2.pfjets_p4()[theJetIndices[jeti]].pt();
                }

                dphipfmetjet_ = mindphipfmet;
                dphitcmetjet_ = mindphitcmet;

                // now, find jet closest to each hyp lepton
                float mindrjet1 = 999999.;
                float mindrjet2 = 999999.;
                float mindrjet3 = 999999.;
                for (unsigned int jeti = 0; jeti < theJetIndices.size(); ++jeti)
                {
                    // for the first lepton
                    float deta1 = type1 == 1 ? cms2.mus_p4()[index1].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta() : cms2.els_p4()[index1].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi1 = type1 == 1 ? cms2.mus_p4()[index1].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi() : cms2.els_p4()[index1].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi();
                    float currdrjet1 = sqrt(deta1*deta1+dphi1*dphi1);
                    if (currdrjet1 < mindrjet1)
                        mindrjet1 = currdrjet1;

                    // for the second lepton
                    float deta2 = type2 == 1 ? cms2.mus_p4()[index2].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta() : cms2.els_p4()[index2].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi2 = type2 == 1 ? cms2.mus_p4()[index2].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi() : cms2.els_p4()[index2].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi();
                    float currdrjet2 = sqrt(deta2*deta2+dphi2*dphi2);
                    if (currdrjet2 < mindrjet2)
                        mindrjet2 = currdrjet2;

                    // for the third lepton
                    float deta3 = type3 == 1 ? cms2.mus_p4()[index3].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta() : cms2.els_p4()[index3].eta()-cms2.pfjets_p4()[theJetIndices[jeti]].eta();
                    float dphi3 = type3 == 1 ? cms2.mus_p4()[index3].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi() : cms2.els_p4()[index3].phi()-cms2.pfjets_p4()[theJetIndices[jeti]].phi();
                    float currdrjet3 = sqrt(deta3*deta3+dphi3*dphi3);
                    if (currdrjet3 < mindrjet3)
                        mindrjet3 = currdrjet3;
                }

                drjet1_ = mindrjet1;
                drjet2_ = mindrjet2;
                drjet3_ = mindrjet3;

                if (type1 == 1)
                {
                    iso1_   = muonIsoValue(index1);
                    type1_  = cms2.mus_type()[index1];					 
                    mu1_muonidfull_   = muonId(index1, NominalTTbarV2); 
                    mu1_muonid_       = muonIdNotIsolated(index1, NominalTTbarV2); 
                    mu1_muonidfullV1_ = muonId(index1, NominalTTbar); 
                    mu1_muonidV1_     = muonIdNotIsolated(index1, NominalTTbar); 
                    mu1_goodmask_     = cms2.mus_goodmask()[index1];
                    mu1_gfitchi2_     = cms2.mus_gfit_chi2()[index1] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index1]/cms2.mus_gfit_ndof()[index1];
                    mu1_cosmic_       = isCosmics(index1);
                    mu1_siHits_       = cms2.mus_validHits()[index1];
                    mu1_saHits_       = cms2.mus_gfit_validSTAHits()[index1];
                    mu1_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[index1];
                    mu1_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[index1];
                    if (cms2.mus_pfmusidx()[index1] > -1)
                        mu1_isPFmuon_ = 1;

                    int trkidx1 = cms2.mus_trkidx()[index1];
                    d0vtx1_ = cms2.trks_d0vtx()[trkidx1];
                }
                else if (type1 == 2)
                {
                    iso1_   = electronIsolation_rel(index1, true);
                    type1_  = cms2.els_type()[index1];

                    e1_cand01full_  = pass_electronSelection(index1, electronSelection_ttbar);
                    e1_cand01_      = electronId_cand(index1, CAND_01);
                    e1_vbtf90full_  = pass_electronSelection(index1, electronSelection_ttbarV2);
                    electronIdComponent_t answer_vbtf90 = electronId_VBTF(index1, VBTF_35X_90);
                    e1_vbtf90_      = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf85 = electronId_VBTF(index1, VBTF_35X_85);
                    e1_vbtf85_      = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf80 = electronId_VBTF(index1, VBTF_35X_80);
                    e1_vbtf80_      = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf70 = electronId_VBTF(index1, VBTF_35X_70);
                    e1_vbtf70_      = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    e1_scet_        = cms2.els_eSC()[index1] / cosh(cms2.els_etaSC()[index1]);
                    e1_eopin_       = cms2.els_eOverPIn()[index1];
                    e1_hoe_         = cms2.els_hOverE()[index1];
                    e1_dphiin_      = cms2.els_dPhiIn()[index1];
                    e1_detain_      = cms2.els_dEtaIn()[index1];
                    e1_e25Me55_     = cms2.els_e2x5Max()[index1] / cms2.els_e5x5()[index1];
                    e1_sigieie_     = cms2.els_sigmaIEtaIEta()[index1];
                    e1_eMe55_       = cms2.els_eMax()[index1] / cms2.els_e5x5()[index1];
                    e1_nmHits_      = cms2.els_exp_innerlayers()[index1];
                    e1_dcot_        = cms2.els_conv_dcot()[index1];
                    e1_dist_        = cms2.els_conv_dist()[index1];
                    e1_drmu_        = cms2.els_closestMuon()[index1] < 0 ? -999999. : cms2.els_musdr()[index1];
                    e1_isspike_     = isSpikeElectron(index1);
                    e1_scCharge_    = cms2.els_sccharge()[index1];
                    e1_gsfCharge_   = cms2.els_trk_charge()[index1];
                    e1_ctfCharge_   = cms2.els_trkidx()[index1] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[index1]] : -999999;

                    int trkidx1 = cms2.els_trkidx()[index1];
                    if (trkidx1 >= 0)
                        d0vtx1_ = cms2.trks_d0vtx()[trkidx1];
                }

                if (type2 == 1)
                {
                    iso2_   = muonIsoValue(index2);
                    type2_  = cms2.mus_type()[index2];					 
                    mu2_muonidfull_   = muonId(index2, NominalTTbarV2); 
                    mu2_muonid_       = muonIdNotIsolated(index2, NominalTTbarV2); 
                    mu2_muonidfullV1_ = muonId(index2, NominalTTbar); 
                    mu2_muonidV1_     = muonIdNotIsolated(index2, NominalTTbar); 
                    mu2_goodmask_     = cms2.mus_goodmask()[index2];
                    mu2_gfitchi2_     = cms2.mus_gfit_chi2()[index2] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index2]/cms2.mus_gfit_ndof()[index2];
                    mu2_cosmic_       = isCosmics(index2);
                    mu2_siHits_       = cms2.mus_validHits()[index2];
                    mu2_saHits_       = cms2.mus_gfit_validSTAHits()[index2];
                    mu2_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[index2];
                    mu2_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[index2];
                    if (cms2.mus_pfmusidx()[index2] > -1)
                        mu2_isPFmuon_ = 1;

                    int trkidx2 = cms2.mus_trkidx()[index2];
                    d0vtx2_ = cms2.trks_d0vtx()[trkidx2];
                }
                else if (type2 == 2)
                {
                    iso2_   = electronIsolation_rel(index2, true);
                    type2_  = cms2.els_type()[index2];

                    e2_cand01full_  = pass_electronSelection(index2, electronSelection_ttbar);
                    e2_cand01_      = electronId_cand(index2, CAND_01);
                    e2_vbtf90full_  = pass_electronSelection(index2, electronSelection_ttbarV2);
                    electronIdComponent_t answer_vbtf90 = electronId_VBTF(index2, VBTF_35X_90);
                    e2_vbtf90_      = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf85 = electronId_VBTF(index2, VBTF_35X_85);
                    e2_vbtf85_      = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf80 = electronId_VBTF(index2, VBTF_35X_80);
                    e2_vbtf80_      = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf70 = electronId_VBTF(index2, VBTF_35X_70);
                    e2_vbtf70_      = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    e2_scet_        = cms2.els_eSC()[index2] / cosh(cms2.els_etaSC()[index2]);
                    e2_eopin_       = cms2.els_eOverPIn()[index2];
                    e2_hoe_         = cms2.els_hOverE()[index2];
                    e2_dphiin_      = cms2.els_dPhiIn()[index2];
                    e2_detain_      = cms2.els_dEtaIn()[index2];
                    e2_e25Me55_     = cms2.els_e2x5Max()[index2] / cms2.els_e5x5()[index2];
                    e2_sigieie_     = cms2.els_sigmaIEtaIEta()[index2];
                    e2_eMe55_       = cms2.els_eMax()[index2] / cms2.els_e5x5()[index2];
                    e2_nmHits_      = cms2.els_exp_innerlayers()[index2];
                    e2_dcot_        = cms2.els_conv_dcot()[index2];
                    e2_dist_        = cms2.els_conv_dist()[index2];
                    e2_drmu_        = cms2.els_closestMuon()[index2] < 0 ? -999999. : cms2.els_musdr()[index2];
                    e2_isspike_     = isSpikeElectron(index2);
                    e2_scCharge_    = cms2.els_sccharge()[index2];
                    e2_gsfCharge_   = cms2.els_trk_charge()[index2];
                    e2_ctfCharge_   = cms2.els_trkidx()[index2] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[index2]] : -999999;

                    int trkidx2 = cms2.els_trkidx()[index2];
                    if (trkidx2 >= 0)
                        d0vtx2_ = cms2.trks_d0vtx()[trkidx2];
                }

                if (type3 == 1)
                {
                    iso3_   = muonIsoValue(index3);
                    type3_  = cms2.mus_type()[index3];					 
                    mu3_muonidfull_   = muonId(index3, NominalTTbarV2); 
                    mu3_muonid_       = muonIdNotIsolated(index3, NominalTTbarV2); 
                    mu3_muonidfullV1_ = muonId(index3, NominalTTbar); 
                    mu3_muonidV1_     = muonIdNotIsolated(index3, NominalTTbar); 
                    mu3_goodmask_     = cms2.mus_goodmask()[index3];
                    mu3_gfitchi2_     = cms2.mus_gfit_chi2()[index3] < -9000. ? -999999. : cms2.mus_gfit_chi2()[index3]/cms2.mus_gfit_ndof()[index3];
                    mu3_cosmic_       = isCosmics(index3);
                    mu3_siHits_       = cms2.mus_validHits()[index3];
                    mu3_saHits_       = cms2.mus_gfit_validSTAHits()[index3];
                    mu3_emVetoDep_    = cms2.mus_iso_ecalvetoDep()[index3];
                    mu3_hadVetoDep_   = cms2.mus_iso_hcalvetoDep()[index3];
                    if (cms2.mus_pfmusidx()[index3] > -1)
                        mu3_isPFmuon_ = 1;

                    int trkidx3 = cms2.mus_trkidx()[index3];
                    d0vtx3_ = cms2.trks_d0vtx()[trkidx3];
                }
                else if (type3 == 2)
                {
                    iso3_   = electronIsolation_rel(index3, true);
                    type3_  = cms2.els_type()[index3];

                    e3_cand01full_  = pass_electronSelection(index3, electronSelection_ttbar);
                    e3_cand01_      = electronId_cand(index3, CAND_01);
                    e3_vbtf90full_  = pass_electronSelection(index3, electronSelection_ttbarV2);
                    electronIdComponent_t answer_vbtf90 = electronId_VBTF(index3, VBTF_35X_90);
                    e3_vbtf90_      = (answer_vbtf90 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf85 = electronId_VBTF(index3, VBTF_35X_85);
                    e3_vbtf85_      = (answer_vbtf85 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf80 = electronId_VBTF(index3, VBTF_35X_80);
                    e3_vbtf80_      = (answer_vbtf80 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    electronIdComponent_t answer_vbtf70 = electronId_VBTF(index3, VBTF_35X_70);
                    e3_vbtf70_      = (answer_vbtf70 & (1ll<<ELEID_ID)) == (1ll<<ELEID_ID);
                    e3_scet_        = cms2.els_eSC()[index3] / cosh(cms2.els_etaSC()[index3]);
                    e3_eopin_       = cms2.els_eOverPIn()[index3];
                    e3_hoe_         = cms2.els_hOverE()[index3];
                    e3_dphiin_      = cms2.els_dPhiIn()[index3];
                    e3_detain_      = cms2.els_dEtaIn()[index3];
                    e3_e25Me55_     = cms2.els_e2x5Max()[index3] / cms2.els_e5x5()[index3];
                    e3_sigieie_     = cms2.els_sigmaIEtaIEta()[index3];
                    e3_eMe55_       = cms2.els_eMax()[index3] / cms2.els_e5x5()[index3];
                    e3_nmHits_      = cms2.els_exp_innerlayers()[index3];
                    e3_dcot_        = cms2.els_conv_dcot()[index3];
                    e3_dist_        = cms2.els_conv_dist()[index3];
                    e3_drmu_        = cms2.els_closestMuon()[index3] < 0 ? -999999. : cms2.els_musdr()[index3];
                    e3_isspike_     = isSpikeElectron(index3);
                    e3_scCharge_    = cms2.els_sccharge()[index3];
                    e3_gsfCharge_   = cms2.els_trk_charge()[index3];
                    e3_ctfCharge_   = cms2.els_trkidx()[index3] > -1 ? cms2.trks_charge()[cms2.els_trkidx()[index3]] : -999999;

                    int trkidx3 = cms2.els_trkidx()[index3];
                    if (trkidx3 >= 0)
                        d0vtx3_ = cms2.trks_d0vtx()[trkidx3];
                }

                FillBabyNtuple();
            }
        }

        if (nEventsChain != nEventsTotal)
            std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
    }

    CloseBabyNtuple();
}

void trilepbabymaker::InitBabyNtuple ()
{
    // event stuff
    memset(dataset_, '\0', 200);
    run_          = -999999;
    ls_           = -999999;
    evt_          = -999999;
    nvtx_         = -999999;
    isdata_       = 1;
    scale1fb_     = -999999.;
    pthat_        = -999999.;
    hyp_type_     = -999999;
    pfmet_        = -999999.;
    tcmet_        = -999999.;
    calotcmet_    = -999999.;
    ntrks_        = -999999;
    njets_        = -999999;
    jet1pt_       = -999999.;
    jet2pt_       = -999999.;
    jet3pt_       = -999999.;
    sumjetpt_     = -999999.;
    jet1eta_      = -999999.;
    jet2eta_      = -999999.;
    jet3eta_      = -999999.;
    jet1phi_      = -999999.;
    jet2phi_      = -999999.;
    jet3phi_      = -999999.;
    jetmass_      = -999999.;
    jet1isBtag_   = 0;
    jet2isBtag_   = 0;
    jet3isBtag_   = 0;
    dphipfmetjet_ = -999999.;
    dphitcmetjet_ = -999999.;
    neffbtags_    = -999999;
    npurbtags_    = -999999;
    pfmeff_       = -999999.;
    tcmeff_       = -999999.;

    // lepton stuff
    ngoodlep_     = -999999;
    ngoodmus_     = -999999;
    ngoodels_     = -999999;
    ngenels_      = -999999;
    ngenmus_      = -999999;
    ngentaus_     = -999999;

    eormu1_       = -999999;
    type1_        = -999999;
    pt1_          = -999999.;
    eta1_         = -999999.;
    phi1_         = -999999.;
    iso1_         = -999999.;
    d0corr1_      = -999999.;
    d0vtx1_       = -999999.;
    dphipfmet1_   = -999999.;
    dphitcmet1_   = -999999.;
    drjet1_       = -999999.;

    eormu2_       = -999999;
    type2_        = -999999;
    pt2_          = -999999.;
    eta2_         = -999999.;
    phi2_         = -999999.;
    iso2_         = -999999.;
    d0corr2_      = -999999.;
    d0vtx2_       = -999999.;
    dphipfmet2_   = -999999.;
    dphitcmet2_   = -999999.;
    drjet2_       = -999999.;

    eormu3_       = -999999;
    type3_        = -999999;
    pt3_          = -999999.;
    eta3_         = -999999.;
    phi3_         = -999999.;
    iso3_         = -999999.;
    d0corr3_      = -999999.;
    d0vtx3_       = -999999.;
    dphipfmet3_   = -999999.;
    dphitcmet3_   = -999999.;
    drjet3_       = -999999.;

    // muon stuff
    mu1_muonidfull_   = 0;
    mu1_muonid_       = 0;
    mu1_muonidfullV1_ = 0;
    mu1_muonidV1_     = 0;
    mu1_goodmask_     = -999999;
    mu1_gfitchi2_     = -999999.;
    mu1_cosmic_       = 0;
    mu1_siHits_       = -999999;
    mu1_saHits_       = -999999;
    mu1_emVetoDep_    = -999999.;
    mu1_hadVetoDep_   = -999999.;

    mu2_muonidfull_   = 0;
    mu2_muonid_       = 0;
    mu2_muonidfullV1_ = 0;
    mu2_muonidV1_     = 0;
    mu2_goodmask_     = -999999;
    mu2_gfitchi2_     = -999999.;
    mu2_cosmic_       = 0;
    mu2_siHits_       = -999999;
    mu2_saHits_       = -999999;
    mu2_emVetoDep_    = -999999.;
    mu2_hadVetoDep_   = -999999.;

    mu3_muonidfull_   = 0;
    mu3_muonid_       = 0;
    mu3_muonidfullV1_ = 0;
    mu3_muonidV1_     = 0;
    mu3_goodmask_     = -999999;
    mu3_gfitchi2_     = -999999.;
    mu3_cosmic_       = 0;
    mu3_siHits_       = -999999;
    mu3_saHits_       = -999999;
    mu3_emVetoDep_    = -999999.;
    mu3_hadVetoDep_   = -999999.;

    // electron stuff
    e1_cand01full_  = 0;
    e1_cand01_      = 0;
    e1_vbtf90full_  = 0;
    e1_vbtf90_      = 0;
    e1_vbtf85_      = 0;
    e1_vbtf80_      = 0;
    e1_vbtf70_      = 0;
    e1_scet_        = -999999.;
    e1_eopin_       = -999999.;
    e1_hoe_         = -999999.;
    e1_dphiin_      = -999999.;
    e1_detain_      = -999999.;
    e1_e25Me55_     = -999999.;
    e1_sigieie_     = -999999.;
    e1_eMe55_       = -999999.;
    e1_nmHits_      = -999999;
    e1_dcot_        = -999999.;
    e1_dist_        = -999999.;
    e1_drmu_        = -999999.;
    e1_isspike_     = 0;
    e1_ctfCharge_   = -999999;
    e1_gsfCharge_   = -999999;
    e1_scCharge_    = -999999;

    e2_cand01full_  = 0;
    e2_cand01_      = 0;
    e2_vbtf90full_  = 0;
    e2_vbtf90_      = 0;
    e2_vbtf85_      = 0;
    e2_vbtf80_      = 0;
    e2_vbtf70_      = 0;
    e2_scet_        = -999999.;
    e2_eopin_       = -999999.;
    e2_hoe_         = -999999.;
    e2_dphiin_      = -999999.;
    e2_detain_      = -999999.;
    e2_e25Me55_     = -999999.;
    e2_sigieie_     = -999999.;
    e2_eMe55_       = -999999.;
    e2_nmHits_      = -999999;
    e2_dcot_        = -999999.;
    e2_dist_        = -999999.;
    e2_drmu_        = -999999.;
    e2_isspike_     = 0;
    e2_ctfCharge_   = -999999;
    e2_gsfCharge_   = -999999;
    e2_scCharge_    = -999999;

    e3_cand01full_  = 0;
    e3_cand01_      = 0;
    e3_vbtf90full_  = 0;
    e3_vbtf90_      = 0;
    e3_vbtf85_      = 0;
    e3_vbtf80_      = 0;
    e3_vbtf70_      = 0;
    e3_scet_        = -999999.;
    e3_eopin_       = -999999.;
    e3_hoe_         = -999999.;
    e3_dphiin_      = -999999.;
    e3_detain_      = -999999.;
    e3_e25Me55_     = -999999.;
    e3_sigieie_     = -999999.;
    e3_eMe55_       = -999999.;
    e3_nmHits_      = -999999;
    e3_dcot_        = -999999.;
    e3_dist_        = -999999.;
    e3_drmu_        = -999999.;
    e3_isspike_     = 0;
    e3_ctfCharge_   = -999999;
    e3_gsfCharge_   = -999999;
    e3_scCharge_    = -999999;
}

void trilepbabymaker::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();

    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    // event stuff
    babyTree_->Branch("dataset",      &dataset_,     "dataset[200]/C");
    babyTree_->Branch("run",          &run_,         "run/I"         );
    babyTree_->Branch("ls",           &ls_,          "ls/I"          );
    babyTree_->Branch("evt",          &evt_,         "evt/I"         );
    babyTree_->Branch("nvtx",         &nvtx_,        "nvtx/I"        );
    babyTree_->Branch("isdata",       &isdata_,      "isdata/I"      );
    babyTree_->Branch("scale1fb",     &scale1fb_,    "scale1fb/F"    );
    babyTree_->Branch("pthat",        &pthat_,       "pthat/F"       );
    babyTree_->Branch("hyp_type",     &hyp_type_,    "hyp_type/I"    );
    babyTree_->Branch("pfmet",        &pfmet_,       "pfmet/F"       );
    babyTree_->Branch("tcmet",        &tcmet_,       "tcmet/F"       );
    babyTree_->Branch("calotcmet",    &calotcmet_,   "calotcmet/F"   );
    babyTree_->Branch("ntrks",        &ntrks_,       "ntrks/I"       );
    babyTree_->Branch("njets",        &njets_,       "njets/I"       ); // uncorrected pt > 20
    babyTree_->Branch("jet1pt",       &jet1pt_,      "jet1pt/F"      );
    babyTree_->Branch("jet2pt",       &jet2pt_,      "jet2pt/F"      );
    babyTree_->Branch("jet3pt",       &jet3pt_,      "jet3pt/F"      );
    babyTree_->Branch("sumjetpt",     &sumjetpt_,    "sumjetpt/F"    );      
    babyTree_->Branch("jet1eta",      &jet1eta_,     "jet1eta/F"     );
    babyTree_->Branch("jet2eta",      &jet2eta_,     "jet2eta/F"     );
    babyTree_->Branch("jet3eta",      &jet3eta_,     "jet3eta/F"     );
    babyTree_->Branch("jet1phi",      &jet1phi_,     "jet1phi/F"     );
    babyTree_->Branch("jet2phi",      &jet2phi_,     "jet2phi/F"     );
    babyTree_->Branch("jet3phi",      &jet3phi_,     "jet3phi/F"     );
    babyTree_->Branch("jetmass",      &jetmass_,     "jetmass/F"     );
    babyTree_->Branch("jet1isBtag",   &jet1isBtag_,  "jet1isBtag/O"  );
    babyTree_->Branch("jet2isBtag",   &jet2isBtag_,  "jet2isBtag/O"  );
    babyTree_->Branch("jet3isBtag",   &jet3isBtag_,  "jet3isBtag/O"  );
    babyTree_->Branch("dphipfmetjet", &dphipfmetjet_,"dphipfmetjet/F");
    babyTree_->Branch("dphitcmetjet", &dphitcmetjet_,"dphitcmetjet/F");
    babyTree_->Branch("neffbtags",    &neffbtags_,   "neffbtags/I"   );
    babyTree_->Branch("npurbtags",    &npurbtags_,   "npurbtags/I"   );
    babyTree_->Branch("ntceffbtags",  &ntceffbtags_, "ntceffbtags/I" );
    babyTree_->Branch("ntcpurbtags",  &ntcpurbtags_, "ntcpurbtags/I" );
    babyTree_->Branch("pfmeff",       &pfmeff_,      "pfmeff/F"      );
    babyTree_->Branch("tcmeff",       &tcmeff_,      "tcmeff/F"      );


    // lepton stuff
    babyTree_->Branch("ngoodlep",   &ngoodlep_,   "ngoodlep/I"  );
    babyTree_->Branch("ngoodmus",   &ngoodmus_,   "ngoodmus/I"  );
    babyTree_->Branch("ngoodels",   &ngoodels_,   "ngoodels/I"  );
    babyTree_->Branch("ngenels",    &ngenels_,    "ngenels/I"   );
    babyTree_->Branch("ngenmus",    &ngenmus_,    "ngenmus/I"   );
    babyTree_->Branch("ngentaus",   &ngentaus_,   "ngentaus/I"  );

    babyTree_->Branch("eormu1",     &eormu1_,     "eormu1/I"    );
    babyTree_->Branch("type1",      &type1_,      "type1/I"     );
    babyTree_->Branch("pt1",        &pt1_,        "pt1/F"       );
    babyTree_->Branch("eta1",       &eta1_,       "eta1/F"      );
    babyTree_->Branch("phi1",       &phi1_,       "phi1/F"      );
    babyTree_->Branch("iso1",       &iso1_,       "iso1/F"      );
    babyTree_->Branch("d0corr1",    &d0corr1_,    "d0corr1/F"   );
    babyTree_->Branch("d0vtx1",     &d0vtx1_,     "d0vtx1/F"    );
    babyTree_->Branch("dphipfmet1", &dphipfmet1_, "dphipfmet1/F");
    babyTree_->Branch("dphitcmet1", &dphitcmet1_, "dphitcmet1/F");
    babyTree_->Branch("drjet1",     &drjet1_,     "drjet1/F"    );

    babyTree_->Branch("eormu2",     &eormu2_,     "eormu2/I"    );
    babyTree_->Branch("type2",      &type2_,      "type2/I"     );
    babyTree_->Branch("pt2",        &pt2_,        "pt2/F"       );
    babyTree_->Branch("eta2",       &eta2_,       "eta2/F"      );
    babyTree_->Branch("phi2",       &phi2_,       "phi2/F"      );
    babyTree_->Branch("iso2",       &iso2_,       "iso2/F"      );
    babyTree_->Branch("d0corr2",    &d0corr2_,    "d0corr2/F"   );
    babyTree_->Branch("d0vtx2",     &d0vtx2_,     "d0vtx2/F"    );
    babyTree_->Branch("dphipfmet2", &dphipfmet2_, "dphipfmet2/F");
    babyTree_->Branch("dphitcmet2", &dphitcmet2_, "dphitcmet2/F");
    babyTree_->Branch("drjet2",     &drjet2_,     "drjet2/F"    );

    babyTree_->Branch("eormu3",     &eormu3_,     "eormu3/I"    );
    babyTree_->Branch("type3",      &type3_,      "type3/I"     );
    babyTree_->Branch("pt3",        &pt3_,        "pt3/F"       );
    babyTree_->Branch("eta3",       &eta3_,       "eta3/F"      );
    babyTree_->Branch("phi3",       &phi3_,       "phi3/F"      );
    babyTree_->Branch("iso3",       &iso3_,       "iso3/F"      );
    babyTree_->Branch("d0corr3",    &d0corr3_,    "d0corr3/F"   );
    babyTree_->Branch("d0vtx3",     &d0vtx3_,     "d0vtx3/F"    );
    babyTree_->Branch("dphipfmet3", &dphipfmet3_, "dphipfmet3/F");
    babyTree_->Branch("dphitcmet3", &dphitcmet3_, "dphitcmet3/F");
    babyTree_->Branch("drjet3",     &drjet3_,     "drjet3/F"    );

    // muon stuff
    babyTree_->Branch("mu1_muonidfull",   &mu1_muonidfull_,   "mu1_muonidfull/O"  );
    babyTree_->Branch("mu1_muonid",       &mu1_muonid_,       "mu1_muonid/O"      );
    babyTree_->Branch("mu1_muonidfullV1", &mu1_muonidfullV1_, "mu1_muonidfullV1/O");
    babyTree_->Branch("mu1_muonidV1",     &mu1_muonidV1_,     "mu1_muonidV1/O"    );
    babyTree_->Branch("mu1_goodmask",     &mu1_goodmask_,     "mu1_goodmask/I"    );
    babyTree_->Branch("mu1_gfitchi2",     &mu1_gfitchi2_,     "mu1_gfitchi2/F"    );
    babyTree_->Branch("mu1_cosmic",       &mu1_cosmic_,       "mu1_cosmic/O"      );
    babyTree_->Branch("mu1_siHits",       &mu1_siHits_,       "mu1_siHits/I"      );
    babyTree_->Branch("mu1_saHits",       &mu1_saHits_,       "mu1_saHits/I"      );
    babyTree_->Branch("mu1_emVetoDep",    &mu1_emVetoDep_,    "mu1_emVetoDep/F"   );
    babyTree_->Branch("mu1_hadVetoDep",   &mu1_hadVetoDep_,   "mu1_hadVetoDep/F"  );
    babyTree_->Branch("mu1_isPFmuon",     &mu1_isPFmuon_,     "mu1_isPFmuon/O"    );

    babyTree_->Branch("mu2_muonidfull",   &mu2_muonidfull_,   "mu2_muonidfull/O"  );
    babyTree_->Branch("mu2_muonid",       &mu2_muonid_,       "mu2_muonid/O"      );
    babyTree_->Branch("mu2_muonidfullV1", &mu2_muonidfullV1_, "mu2_muonidfullV1/O");
    babyTree_->Branch("mu2_muonidV1",     &mu2_muonidV1_,     "mu2_muonidV1/O"    );
    babyTree_->Branch("mu2_goodmask",     &mu2_goodmask_,     "mu2_goodmask/I"    );
    babyTree_->Branch("mu2_gfitchi2",     &mu2_gfitchi2_,     "mu2_gfitchi2/F"    );
    babyTree_->Branch("mu2_cosmic",       &mu2_cosmic_,       "mu2_cosmic/O"      );
    babyTree_->Branch("mu2_siHits",       &mu2_siHits_,       "mu2_siHits/I"      );
    babyTree_->Branch("mu2_saHits",       &mu2_saHits_,       "mu2_saHits/I"      );
    babyTree_->Branch("mu2_emVetoDep",    &mu2_emVetoDep_,    "mu2_emVetoDep/F"   );
    babyTree_->Branch("mu2_hadVetoDep",   &mu2_hadVetoDep_,   "mu2_hadVetoDep/F"  );
    babyTree_->Branch("mu2_isPFmuon",     &mu2_isPFmuon_,     "mu2_isPFmuon/O"    );

    babyTree_->Branch("mu3_muonidfull",   &mu3_muonidfull_,   "mu3_muonidfull/O"  );
    babyTree_->Branch("mu3_muonid",       &mu3_muonid_,       "mu3_muonid/O"      );
    babyTree_->Branch("mu3_muonidfullV1", &mu3_muonidfullV1_, "mu3_muonidfullV1/O");
    babyTree_->Branch("mu3_muonidV1",     &mu3_muonidV1_,     "mu3_muonidV1/O"    );
    babyTree_->Branch("mu3_goodmask",     &mu3_goodmask_,     "mu3_goodmask/I"    );
    babyTree_->Branch("mu3_gfitchi2",     &mu3_gfitchi2_,     "mu3_gfitchi2/F"    );
    babyTree_->Branch("mu3_cosmic",       &mu3_cosmic_,       "mu3_cosmic/O"      );
    babyTree_->Branch("mu3_siHits",       &mu3_siHits_,       "mu3_siHits/I"      );
    babyTree_->Branch("mu3_saHits",       &mu3_saHits_,       "mu3_saHits/I"      );
    babyTree_->Branch("mu3_emVetoDep",    &mu3_emVetoDep_,    "mu3_emVetoDep/F"   );
    babyTree_->Branch("mu3_hadVetoDep",   &mu3_hadVetoDep_,   "mu3_hadVetoDep/F"  );
    babyTree_->Branch("mu3_isPFmuon",     &mu3_isPFmuon_,     "mu3_isPFmuon/O"    );

    // electron stuff
    babyTree_->Branch("e1_cand01full", &e1_cand01full_, "e1_cand01full/O");
    babyTree_->Branch("e1_cand01",     &e1_cand01_,     "e1_cand01/O"    );
    babyTree_->Branch("e1_vbtf90full", &e1_vbtf90full_, "e1_vbtf90full/O");
    babyTree_->Branch("e1_vbtf90",     &e1_vbtf90_,     "e1_vbtf90/O"    );
    babyTree_->Branch("e1_vbtf85",     &e1_vbtf85_,     "e1_vbtf85/O"    );
    babyTree_->Branch("e1_vbtf80",     &e1_vbtf80_,     "e1_vbtf80/O"    );
    babyTree_->Branch("e1_vbtf70",     &e1_vbtf70_,     "e1_vbtf70/O"    );
    babyTree_->Branch("e1_scet",       &e1_scet_,       "e1_scet/F"      );
    babyTree_->Branch("e1_eopin",      &e1_eopin_,      "e1_eopin/F"     );
    babyTree_->Branch("e1_hoe",        &e1_hoe_,        "e1_hoe/F"       );
    babyTree_->Branch("e1_dphiin",     &e1_dphiin_,     "e1_dphiin/F"    );
    babyTree_->Branch("e1_detain",     &e1_detain_,     "e1_detain/F"    );
    babyTree_->Branch("e1_e25Me55",    &e1_e25Me55_,    "e1_e25Me55/F"   );
    babyTree_->Branch("e1_sigieie",    &e1_sigieie_,    "e1_sigieie/F"   );
    babyTree_->Branch("e1_eMe55",      &e1_eMe55_,      "e1_eMe55/F"     ); // for spikes
    babyTree_->Branch("e1_nmHits",     &e1_nmHits_,     "e1_nmHits/I"    );
    babyTree_->Branch("e1_dcot",       &e1_dcot_,       "e1_dcot/F"      );
    babyTree_->Branch("e1_dist",       &e1_dist_,       "e1_dist/F"      );
    babyTree_->Branch("e1_drmu",       &e1_drmu_,       "e1_drmu/F"      );
    babyTree_->Branch("e1_isspike",    &e1_isspike_,    "e1_isspike/O"   );
    babyTree_->Branch("e1_ctfCharge",  &e1_ctfCharge_,  "e1_ctfCharge/I" );
    babyTree_->Branch("e1_gsfCharge",  &e1_gsfCharge_,  "e1_gsfCharge/I" );
    babyTree_->Branch("e1_scCharge",   &e1_scCharge_,   "e1_scCharge/I"  );

    babyTree_->Branch("e2_cand01full", &e2_cand01full_, "e2_cand01full/O");
    babyTree_->Branch("e2_cand01",     &e2_cand01_,     "e2_cand01/O"    );
    babyTree_->Branch("e2_vbtf90full", &e2_vbtf90full_, "e2_vbtf90full/O");
    babyTree_->Branch("e2_vbtf90",     &e2_vbtf90_,     "e2_vbtf90/O"    );
    babyTree_->Branch("e2_vbtf85",     &e2_vbtf85_,     "e2_vbtf85/O"    );
    babyTree_->Branch("e2_vbtf80",     &e2_vbtf80_,     "e2_vbtf80/O"    );
    babyTree_->Branch("e2_vbtf70",     &e2_vbtf70_,     "e2_vbtf70/O"    );
    babyTree_->Branch("e2_scet",       &e2_scet_,       "e2_scet/F"      );
    babyTree_->Branch("e2_eopin",      &e2_eopin_,      "e2_eopin/F"     );
    babyTree_->Branch("e2_hoe",        &e2_hoe_,        "e2_hoe/F"       );
    babyTree_->Branch("e2_dphiin",     &e2_dphiin_,     "e2_dphiin/F"    );
    babyTree_->Branch("e2_detain",     &e2_detain_,     "e2_detain/F"    );
    babyTree_->Branch("e2_e25Me55",    &e2_e25Me55_,    "e2_e25Me55/F"   );
    babyTree_->Branch("e2_sigieie",    &e2_sigieie_,    "e2_sigieie/F"   );
    babyTree_->Branch("e2_eMe55",      &e2_eMe55_,      "e2_eMe55/F"     ); // for spikes
    babyTree_->Branch("e2_nmHits",     &e2_nmHits_,     "e2_nmHits/I"    );
    babyTree_->Branch("e2_dcot",       &e2_dcot_,       "e2_dcot/F"      );
    babyTree_->Branch("e2_dist",       &e2_dist_,       "e2_dist/F"      );
    babyTree_->Branch("e2_drmu",       &e2_drmu_,       "e2_drmu/F"      );
    babyTree_->Branch("e2_isspike",    &e2_isspike_,    "e2_isspike/O"   );
    babyTree_->Branch("e2_ctfCharge",  &e2_ctfCharge_,  "e2_ctfCharge/I" );
    babyTree_->Branch("e2_gsfCharge",  &e2_gsfCharge_,  "e2_gsfCharge/I" );
    babyTree_->Branch("e2_scCharge",   &e2_scCharge_,   "e2_scCharge/I"  );

    babyTree_->Branch("e3_cand01full", &e3_cand01full_, "e3_cand01full/O");
    babyTree_->Branch("e3_cand01",     &e3_cand01_,     "e3_cand01/O"    );
    babyTree_->Branch("e3_vbtf90full", &e3_vbtf90full_, "e3_vbtf90full/O");
    babyTree_->Branch("e3_vbtf90",     &e3_vbtf90_,     "e3_vbtf90/O"    );
    babyTree_->Branch("e3_vbtf85",     &e3_vbtf85_,     "e3_vbtf85/O"    );
    babyTree_->Branch("e3_vbtf80",     &e3_vbtf80_,     "e3_vbtf80/O"    );
    babyTree_->Branch("e3_vbtf70",     &e3_vbtf70_,     "e3_vbtf70/O"    );
    babyTree_->Branch("e3_scet",       &e3_scet_,       "e3_scet/F"      );
    babyTree_->Branch("e3_eopin",      &e3_eopin_,      "e3_eopin/F"     );
    babyTree_->Branch("e3_hoe",        &e3_hoe_,        "e3_hoe/F"       );
    babyTree_->Branch("e3_dphiin",     &e3_dphiin_,     "e3_dphiin/F"    );
    babyTree_->Branch("e3_detain",     &e3_detain_,     "e3_detain/F"    );
    babyTree_->Branch("e3_e25Me55",    &e3_e25Me55_,    "e3_e25Me55/F"   );
    babyTree_->Branch("e3_sigieie",    &e3_sigieie_,    "e3_sigieie/F"   );
    babyTree_->Branch("e3_eMe55",      &e3_eMe55_,      "e3_eMe55/F"     ); // for spikes
    babyTree_->Branch("e3_nmHits",     &e3_nmHits_,     "e3_nmHits/I"    );
    babyTree_->Branch("e3_dcot",       &e3_dcot_,       "e3_dcot/F"      );
    babyTree_->Branch("e3_dist",       &e3_dist_,       "e3_dist/F"      );
    babyTree_->Branch("e3_drmu",       &e3_drmu_,       "e3_drmu/F"      );
    babyTree_->Branch("e3_isspike",    &e3_isspike_,    "e3_isspike/O"   );
    babyTree_->Branch("e3_ctfCharge",  &e3_ctfCharge_,  "e3_ctfCharge/I" );
    babyTree_->Branch("e3_gsfCharge",  &e3_gsfCharge_,  "e3_gsfCharge/I" );
    babyTree_->Branch("e3_scCharge",   &e3_scCharge_,   "e3_scCharge/I"  );
}

void trilepbabymaker::FillBabyNtuple()
{
    babyTree_->Fill();
}

void trilepbabymaker::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}
