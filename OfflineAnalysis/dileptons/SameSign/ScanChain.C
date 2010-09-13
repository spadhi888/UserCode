/* 
   Usage:
   root [0] .L ScanChain.C++
   root [1] TFile *_file0 = TFile::Open("merged_ntuple.root")
   root [2] TChain *chain = new TChain("Events")
   root [3] chain->Add("merged_ntuple.root")

   There are several places where one may create CMS2 cms2
   It can be done here (in a doAll.C script), i.e.:

   root [4] CMS2 cms2 
   It can be done in the source as is done below, or it can be
   ascertained by including CORE/CMS2.cc as is commented out
   below.  They are all the same, and everything will work so
   long as it is created somewhere globally.

   root [5] ScanChain(chain)
*/
#include <stdexcept>

#include <iostream>
#include <vector>
#include <algorithm>

#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TString.h"
#include "selections.cc"
#include "Histograms.cc"
#include "CMS2.cc"
#include "goodrun.cc"
#include "../CORE/mcSelections.cc"
#include "../CORE/MT2/MT2.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/SimpleFakeRate.cc"
#include "../CORE/triggerUtils.cc"
#include "../NtupleMacros/Tools/fliprate_egun.cc"
#include "../CORE/jetSelections.cc"


#define JETPTCUT 30.0

using namespace tas;
unsigned int numtightLeps;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool usePtGt2020         = false;
bool usePtGt2010         = false;
bool usePtGt1010         = false;
bool excludePtGt2020     = false;
bool used0corrPV         = false;
bool applylepIDCuts      = false;
bool applyFOv1Cuts       = false;
bool applyFOv2Cuts       = false;
bool applyFOv3Cuts       = false;
bool applylepIsoCuts     = false;
bool applyTriggers       = false;
bool vetoZmass           = false;
bool requireZmass        = false;
bool hypDisamb           = false;
bool useCorMET           = false;
bool useOSleptons        = false;
bool useSSleptons        = false;
bool usetcMET            = false;
bool usepfMET            = false;
bool vetoMET             = false;
bool vetoProjectedMET    = false;
bool usejptJets          = false;
bool usecaloJets         = false;
bool usepfJets           = false;
bool vetoJets            = false;
bool requireEcalEls      = false;
bool estimateDoubleFakes = false;
bool estimateSingleFakes = false;
bool useFlipRateEstimation = false;
bool applyAlignmentCorrection  = false;
bool removedEtaCutInEndcap = false;
bool chargeFlip        = false;

double sumfr=0;

Bool_t sortByPt(LorentzVector lv1, LorentzVector lv2) {
  return lv1.pt() > lv2.pt();
}

// Resurrect the DorkyEventIdentifier  
struct DorkyEventIdentifier {
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
       return run < other.run;
     if (event != other.event)
       return event < other.event;
     if(lumi != other.lumi)
       return lumi < other.lumi;
     return false;
}

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
     if (run != other.run)
          return false;
     if (event != other.event)
          return false;
     return true;
}

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
     std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
          already_seen.insert(id);
     return !ret.second;
}


void FillHistograms(const unsigned int hypIdx, const vector<unsigned int> v_jets, const vector<unsigned int> v_jetsNoEtaCut,
		    const pair<float, float>, const float weight, std::string prefix);
double getFRWeight(const int hypIdx, string  elFRversion);
float GetValueTH2F(Float_t x, Float_t y, const TH2F* h);
void EndJob();



void ScanChain( TChain* chain, vector<TString> v_Cuts, string prefix="", 
		bool doFRestimation = false, float jetTriggerPt = -9999., float lumi = 0.01, float NLOCS=-9999.) { //lumi in fb-1


  //deal with the cuts
  useOSleptons         = find(v_Cuts.begin(), v_Cuts.end(), "useOSleptons"         ) != v_Cuts.end();
  useSSleptons         = find(v_Cuts.begin(), v_Cuts.end(), "useSSleptons"         ) != v_Cuts.end();
  usePtGt2020          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2020"          ) != v_Cuts.end();
  usePtGt2010          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2010"          ) != v_Cuts.end();
  usePtGt1010          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt1010"          ) != v_Cuts.end();      
  used0corrPV          = find(v_Cuts.begin(), v_Cuts.end(), "used0corrPV"          ) != v_Cuts.end();
  excludePtGt2020      = find(v_Cuts.begin(), v_Cuts.end(), "excludePtGt2020"      ) != v_Cuts.end();
  applylepIDCuts       = find(v_Cuts.begin(), v_Cuts.end(), "applylepIDCuts"       ) != v_Cuts.end(); 
  applyFOv1Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv1Cuts"        ) != v_Cuts.end(); 
  applyFOv2Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv2Cuts"        ) != v_Cuts.end(); 
  applyFOv3Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv3Cuts"        ) != v_Cuts.end(); 
  applylepIsoCuts      = find(v_Cuts.begin(), v_Cuts.end(), "applylepIsoCuts"      ) != v_Cuts.end(); 
  applyTriggers        = find(v_Cuts.begin(), v_Cuts.end(), "applyTriggers"        ) != v_Cuts.end();
  vetoZmass            = find(v_Cuts.begin(), v_Cuts.end(), "vetoZmass"            ) != v_Cuts.end();
  requireZmass         = find(v_Cuts.begin(), v_Cuts.end(), "requireZmass"         ) != v_Cuts.end();
  hypDisamb            = find(v_Cuts.begin(), v_Cuts.end(), "hypDisamb"            ) != v_Cuts.end();
  useCorMET            = find(v_Cuts.begin(), v_Cuts.end(), "useCorMET"            ) != v_Cuts.end();
  usetcMET             = find(v_Cuts.begin(), v_Cuts.end(), "usetcMET"             ) != v_Cuts.end();   
  usepfMET             = find(v_Cuts.begin(), v_Cuts.end(), "usepfMET"             ) != v_Cuts.end();
  vetoMET              = find(v_Cuts.begin(), v_Cuts.end(), "vetoMET"              ) != v_Cuts.end();
  vetoProjectedMET     = find(v_Cuts.begin(), v_Cuts.end(), "vetoProjectedMET"     ) != v_Cuts.end();
  usecaloJets          = find(v_Cuts.begin(), v_Cuts.end(), "usecaloJets"          ) != v_Cuts.end();
  usejptJets           = find(v_Cuts.begin(), v_Cuts.end(), "usejptJets"           ) != v_Cuts.end();
  usepfJets            = find(v_Cuts.begin(), v_Cuts.end(), "usepfJets"            ) != v_Cuts.end();
  chargeFlip           = find(v_Cuts.begin(), v_Cuts.end(), "chargeFlip"           ) != v_Cuts.end();
  vetoJets             = find(v_Cuts.begin(), v_Cuts.end(), "vetoJets"             ) != v_Cuts.end();
  requireEcalEls       = find(v_Cuts.begin(), v_Cuts.end(), "requireEcalEls"       ) != v_Cuts.end();
  estimateDoubleFakes  = find(v_Cuts.begin(), v_Cuts.end(), "estimateDoubleFakes"  ) != v_Cuts.end();
  estimateSingleFakes  = find(v_Cuts.begin(), v_Cuts.end(), "estimateSingleFakes"  ) != v_Cuts.end();
  useFlipRateEstimation  = find(v_Cuts.begin(), v_Cuts.end(), "useFlipRateEstimation"  ) != v_Cuts.end();
  applyAlignmentCorrection  = find(v_Cuts.begin(), v_Cuts.end(), "applyAlignmentCorrection"  ) != v_Cuts.end();
  removedEtaCutInEndcap  = find(v_Cuts.begin(), v_Cuts.end(), "removedEtaCutInEndcap"  ) != v_Cuts.end();


  cout << "REMEMBER THAT JET COUNTING NOW USES " << JETPTCUT << " GeV JETS" << endl;


  if(doFRestimation && jetTriggerPt < 0) {
    cout << "Need to specify the jetTrigger to do the FRs. Exiting" << endl;
    return;
  }

  if (doFRestimation) {
    if(estimateDoubleFakes + estimateSingleFakes == 0) {
      cout << "Use either estimateDoubleFakes for QCD or estimateSingleFakes for ttbar and WJets" << endl;
      return;
    }

    if(estimateDoubleFakes + estimateSingleFakes > 1) {
      cout << "Pick One: estimateDoubleFakes or estimateSingleFakes" << endl;
      return;
    }
  }

  if (useFlipRateEstimation + useSSleptons > 1)
    cout << " Pick One: Cannot estimate FlipRate using SS leptons " << endl;


  if(usePtGt2020 && usePtGt1010) {
    cout << "You want to use pt > 20,20, and >10,10. While this is ok, its a bit redundant. Please pick the cut you really want" << endl;
    return;
  }

  if(!usePtGt2020 && excludePtGt2020) {
    cout << "Cannot exlcude higer pt region and not use the lower Dilepton pt region. Will end up with no events passing in the greater than 2 jet bin" << endl;
    return;
  }

  if(applylepIDCuts + applyFOv1Cuts + applyFOv2Cuts +applyFOv3Cuts > 1 ) {
    cout << "You have selected too maky ID choices. Please pick one" << endl;
    return;
  }

  if(useCorMET + usetcMET + usepfMET > 1) {
    cout << "Can't use two or more MET algoritms at the same time! Make up your mind!" << endl;
    return;
  }

  if(useCorMET + usetcMET + usepfMET == 0) {
    cout << "Need to use at least one MET algo to make the MT2 and MT2J plots. Please pick either calo, tc or pf MET" << endl;
    return;
  }

  if(usejptJets + usecaloJets +usepfJets > 1) {
    cout << "Can't use two jet algoritms at the same time! Make up your mind!" << endl;
    return;
  }

  if(usejptJets + usecaloJets + usepfJets== 0) {
    cout << "Need to pick at least one jet algoritm for jet counting. Pick either calo or jpt jets!" << endl;
    return;
  }

  if(requireZmass + vetoZmass > 1) {
    cout << "Can't both require and reject events in the Z mass windown at the same time!" << endl;
    return;
  }

  if(vetoMET + vetoProjectedMET > 1) {
    cout << "Projected MET includes the MET veto. If you want projected MET, unselect \"vetoMET\"" << endl;
    return;
  }

  if(usetcMET)
    cout << "REMEMBER THAT WERE COMBINING THE SPRING10 AND SUMMER09 SAMPLES. THERE ARE NO CALOTOWERS STORED IN THE SPRING10 SAMPLES. SO WE CANT CORRECT THE TCMET ON THE FLY AND THERE ARE SOME IF STATEMENTS THAT NEED TO BE CLEANED UP WHEN WE MOVE STRAIGHT TO THE SPRING SAMPLES" << endl;

// JEC 

  std::vector<std::string> jetcorr_filenames;
  if (usecaloJets)
  {
         jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5Calo.txt");
         jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5Calo.txt");
  } 
  else if (usejptJets)
  {
         jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5JPT.txt");
         jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5JPT.txt");
  }
  else if (usepfJets)
  {
        jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt");
        jetcorr_filenames.push_back("../CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");
  }

  FactorizedJetCorrector *jet_corrector = makeJetCorrector(jetcorr_filenames);


  TObjArray *listOfFiles = chain->GetListOfFiles();
  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  numtightLeps = 0;
  unsigned int nGoodHyps[4] = {0,0,0,0};
  unsigned int nGoodEvts = 0;
  
  
  //fix kfactor for DYmm, tautau events
  //if(prefix == "DYmm") {
  //cout << "SCALING DYMM AND TAUTAU EVENTS BY 1.14 BECAUSE WE FORGOT THE KFACTORS. BE CAREFUL!!!!!!" << endl;
  //lumi = lumi * 1.14;
  //}
  //bool isDYee;
  //TObjArray *objArray = chain->GetListOfFiles();
  //if(


  //book Histograms
  bookHistos(prefix.c_str());

  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  //set the good run file
  set_goodrun_file("goodruns.txt");
  

  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);

// deal with MC later
// set the alignement corrections for data to be true

      if(prefix == "data") {
        applyAlignmentCorrection = true;
        DorkyEventIdentifier id = { evt_run(),evt_event(),evt_lumiBlock() };
        if (is_duplicate(id)) continue;
      } else {
        applyAlignmentCorrection = false;
      }

      ++nEventsTotal;
      // Progress feedback to the user
      if(nEventsTotal%2000 == 0) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
	  fflush(stdout);
	}
      }//if(nEventsTotal%20000 == 0) {
      //select good runs, lumis
      if(prefix == "data" && !goodrun(evt_run(), evt_lumiBlock())) continue;

      // Cleaning cuts.

      if (!cleaning_standardAugust2010()) continue;

      //get the channels correct
      int nels = 0;
      int nmus = 0;
      int ntaus = 0;
      int nleps = 0;
      if(prefix != "data")
	nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
      if (prefix == "ttdil"   &&  nleps  != 2) continue;
      if (prefix == "ttotr"    &&  nleps == 2) continue;
      if (prefix == "DYee"     &&  nels != 2) continue;
      if (prefix == "DYmm"     &&  nmus != 2) continue;
      if (prefix == "DYtautau" &&  ntaus != 2) continue;

      float pthat_cutoff = 30.;
      if (prefix == "qcdpt15" && genps_pthat() > pthat_cutoff) continue;

/*
      //splice together the DY samples - if its madgraph, then we do nothing
      if(TString(prefix).Contains("DY")) {
	   bool doNotContinue = false;
	   if (TString(evt_dataset()).Contains("madgraph") == true) {
	     for(unsigned int i = 0; i < genps_p4().size(); i++){
	       if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() < 50.) doNotContinue = true;
	     }
	   } else if (TString(evt_dataset()).Contains("M10to20") == true) {
	     // Do 10-20
	     for(unsigned int i = 0; i < genps_p4().size(); i++){
	       if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 20.) doNotContinue = true;
	     }
	   } else {
	     // do 20-50
             for(unsigned int i = 0; i < genps_p4().size(); i++){
               if(abs(genps_id()[i]) == 23 && (genps_p4()[i].M() < 20. || genps_p4()[i].M() > 50.)) doNotContinue = true;
             }
	   }
	   if(doNotContinue) continue;
	 }
*/

      //splice together the DY samples - if its madgraph, then we do nothing
      if(TString(prefix).Contains("DY") && TString(evt_dataset()).Contains("madgraph") == false) {      
        bool doNotContinue = false;
        for(unsigned int i = 0; i < genps_p4().size(); i++){
          if(abs(genps_id()[i]) == 23 && genps_p4()[i].M() > 50.)
            doNotContinue = true;
        }
        if(doNotContinue)
          continue;     
      }

      //10pb-1 
      float weight = 1.;
      if(prefix != "data") {
	if(NLOCS > 0)
	  weight = evt_scale1fb()*lumi*NLOCS/(evt_xsec_incl()*evt_kfactor());
	else 
	  weight = evt_scale1fb()*lumi;	  
      }

      if(TString(prefix).Contains("DY") ) {
        if(TString(evt_dataset()).Contains("madgraph") == true) { //mll > 50
          //weight = weight*2800./2400./kFactor;
          weight = weight*3048./2400.;
        } else if(TString(evt_dataset()).Contains("M10to20") == true) { //10 < mll < 20 
          weight = weight*3457./2659.;
        } else { // 20 < mll < 50
          weight = weight * 1666./1300.;
        }
      }

      vector<unsigned int> v_goodHyps;
      v_goodHyps.clear();
      vector<float> v_weights;
      v_weights.clear();
      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); hypIdx++) {//hyploop
	
	int id_lt = hyp_lt_id()[hypIdx]; 
	int id_ll = hyp_ll_id()[hypIdx];
	int idx_lt = hyp_lt_index()[hypIdx];
	int idx_ll = hyp_ll_index()[hypIdx];
	LorentzVector lt_p4 = hyp_lt_p4()[hypIdx]; 
	LorentzVector ll_p4 = hyp_ll_p4()[hypIdx];
	int type = hyp_type()[hypIdx];	

//        cout << "run,lumi,number: " << evt_run() << "," << evt_lumiBlock() << "," << evt_event() << " hyp Idx: " << hypIdx << endl;
//        cout << "lt Pt: " << lt_p4.Pt() << " ll Pt" << ll_p4.Pt() << endl;

	// opposite charge
	if (useOSleptons)  {
	  if (id_lt * id_ll > 0) continue;
//          if ((els_trk_charge().at(idx_lt))*(els_trk_charge().at(idx_ll)) > 0 ) continue; // GSF
//          if (cms2.els_trkidx().at(idx_lt) < 0) continue; //CTF
//          if (cms2.els_trkidx().at(idx_ll) < 0) continue; //CTF
//          if ((cms2.trks_charge().at(cms2.els_trkidx().at(idx_lt)))*(cms2.trks_charge().at(cms2.els_trkidx().at(idx_ll))) > 0) continue; //CTF
//          if ((cms2.els_sccharge().at(idx_lt))*(cms2.els_sccharge().at(idx_ll)) > 0 ) continue; //SCCharge
        }

	// same charge
        if (useSSleptons) {
          if (id_lt * id_ll < 0) continue;
  //        if ((cms2.els_sccharge().at(idx_lt))*(cms2.els_sccharge().at(idx_ll)) < 0 ) continue; //SCCharge
  //        if ((els_trk_charge().at(idx_lt))*(els_trk_charge().at(idx_ll)) < 0 ) continue; // GSF
  //        if (cms2.els_trkidx().at(idx_lt) < 0) continue; //CTF
  //        if (cms2.els_trkidx().at(idx_ll) < 0) continue; //CTF
  //        if ((cms2.trks_charge().at(cms2.els_trkidx().at(idx_lt)))*(cms2.trks_charge().at(cms2.els_trkidx().at(idx_ll))) < 0) continue; //CTF
        }

	//if a muon, always require global and tracker
	if(abs(id_lt)==13) {
	  if (((mus_type()[idx_lt]) & (1<<1)) == 0)    continue; // global muon
	  if (((mus_type()[idx_lt]) & (1<<2)) == 0)    continue; // tracker muon
	  
	}
	if(abs(id_ll)==13) {
	  if (((mus_type()[idx_ll]) & (1<<1)) == 0)    continue; // global muon
	  if (((mus_type()[idx_ll]) & (1<<2)) == 0)    continue; // tracker muon
	}

	if(requireEcalEls) {
	  //ask that the electron is ecal driven
	  if(abs(id_lt) == 11) {
	    if (!(els_type()[idx_lt] & (1<<2)))
	      continue;
	  }
	  if(abs(id_ll) == 11) {
	    if (!(els_type()[idx_ll] & (1<<2)))
	      continue;
	  }
	}

	//regardless of jet bin
	//cut at Pt > 20, 20
	if(usePtGt2020) {
	  if(lt_p4.Pt() < 20. || ll_p4.Pt() < 20.)
	    continue;
	}

	if(usePtGt1010) {
	  if(lt_p4.Pt() < 10. || ll_p4.Pt() < 10.)
	    continue;
	}


	//cut at tight Pt > 20, loose Pt > 10
	if(usePtGt2010) {
	  if(TMath::Max(lt_p4.Pt(),ll_p4.Pt()) < 20)
	    continue;
	  if(TMath::Min(lt_p4.Pt(),ll_p4.Pt()) < 10)
	    continue;
	}
	

	//only look at events where loose lepton Pt < 20, tight > 20
	if(excludePtGt2020) {
	  if(lt_p4.Pt() > 20. && ll_p4.Pt() > 20.)
	    continue; 
	  if(lt_p4.Pt() < 10 || ll_p4.Pt() < 10.)
	    continue; 
	}

        //require trigger?
        if(applyTriggers) {
        //  bool passMu = passHLTTrigger("HLT_Mu9");
          bool runningOnMC = true;
          if(prefix == "data") runningOnMC = false;

          bool passEl = passEGTrigger(runningOnMC, type);
          bool passMu = passMuTrigger(runningOnMC, type);

          if(type == 0 && !passMu) continue;
          if(type == 3 && !passEl) continue;
          if((type == 1 || type == 2) && !passMu && !passEl) continue;

        }//applyTriggers

        //apply the vertex requirement

        if(!hypsFromSameVtx(hypIdx)) continue;
	/*
	//require trigger? - From Claudio
	if(applyTriggers) {
	  bool passMu = passHLTTrigger("HLT_Mu9");
	  bool passEl = false;
          if(prefix == "data") {
            passEl = goodEGTrigger5July2010(0); 
          } else {
            passEl = goodEGTrigger5July2010(1); 
          }

	  if(type == 0 && !passMu)
	    continue;
	  if(type == 3 && !passEl)
	    continue;
	  if((type == 1 || type == 2) && !passMu && !passEl)
	    continue;
//         if (!passTriggersMu9orLisoE15(type)) continue;
 	}//applyTriggers

	*/

	// Additional Z Veto to supress WZ and ZZ

//        if ( additionalZvetoSUSY2010(hypIdx, applyAlignmentCorrection, removedEtaCutInEndcap)) continue;

//       z mass window
	if(vetoZmass) {
//	  if(type == 0 || type == 3) {
	  if(type == 3) {
	    if (inZmassWindow(hyp_p4()[hypIdx].mass())) 
	      continue;
	  }
	}//vetoZmass


// Require Z mass in ee only

	if(requireZmass) {
//	  if(type == 0 || type == 3) {
	  if(type == 3) {
	    if (!inZmassWindow(hyp_p4()[hypIdx].mass())) 
	      continue;
	  } else 
              continue;
	}//requireZmass


        // Lepton from W
        bool fromw = false;
        if (abs(id_ll) == 13 && leptonIsFromW(hyp_ll_index()[hypIdx], hyp_ll_id()[hypIdx]) == 1) fromw = true;
        if (abs(id_lt) == 13 && leptonIsFromW(hyp_lt_index()[hypIdx], hyp_lt_id()[hypIdx]) == 1) fromw = true;
        if (!fromw) continue;


	if(doFRestimation) {
	  //unsigned int elFRversion = 9999; FR versions and string
	  string elFRversion;
	  if(applyFOv1Cuts)
	    elFRversion = Form("eFRv1%du", (int)jetTriggerPt);
	  else if(applyFOv2Cuts)
	    elFRversion = Form("eFRv2%du", (int)jetTriggerPt);
	  else if(applyFOv3Cuts)
	    elFRversion = Form("eFRv3%du", (int)jetTriggerPt);
          else
            elFRversion = "undefined"; 

	  if(elFRversion == "undefined") {
	    cout << "asking for FR version that is not supported. quitting" << endl;
	    return;
	  }
	  float FRweight = getFRWeight(hypIdx, elFRversion); 
	  if(FRweight < -1.) 
	    continue;
	  v_goodHyps.push_back(hypIdx);
	  v_weights.push_back(FRweight); 
	  continue;
	}
	
	//require lepton ID cuts?
	if(applylepIDCuts) {
	  if(!isGoodHypNoIso(hypIdx, applyAlignmentCorrection, removedEtaCutInEndcap))//, used0corrPV) 
	    continue;
	}
	
	//common for muons, so do all here
	if(applyFOv1Cuts || applyFOv2Cuts || applyFOv3Cuts) {

	  if (estimateDoubleFakes) {
 
	    if(abs(id_lt) == 13) 	
	      if(!isFakeableMuon(idx_lt)) continue;
	    if(abs(id_ll) == 13) 
	      if(!isFakeableMuon(idx_ll)) continue;
	  
	    if(abs(id_lt) == 11) {
	      if(applyFOv1Cuts && !pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	      if(applyFOv2Cuts && !pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	      if(applyFOv3Cuts && !pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	    }
	    if(abs(id_ll) == 11) {
	      if(applyFOv1Cuts && !pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	      if(applyFOv2Cuts && !pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	      if(applyFOv3Cuts && !pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap) ) continue;
	    }
	  }

	  if (estimateSingleFakes){
	    // mumu only
	    if(hyp_type()[hypIdx] == 0) {

	      bool isGoodMlt = false;
	      bool isGoodMll = false;
	      bool isFOMlt   = false;
	      bool isFOMll   = false;
	      bool evtFRmm   = false;
	      
	      if (isGoodLeptonwIso(13, idx_lt, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodMlt = true;
	      if (isGoodLeptonwIso(13, idx_ll, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodMll = true;
	      if (isFakeableMuon(idx_lt)) isFOMlt = true;
	      if (isFakeableMuon(idx_ll)) isFOMll = true;
	      if (isGoodMlt && isGoodMll) continue;
	      if (isGoodMlt && !isGoodMll && isFOMll) evtFRmm = true;
	      if (isGoodMll && !isGoodMlt && isFOMlt) evtFRmm = true;
	      if ( isGoodMlt && isGoodMll) continue;
	      if (!evtFRmm) continue;
	    }
	    // ee case
	    if(hyp_type()[hypIdx] == 3) {

	      bool isGoodElt = false;
	      bool isGoodEll = false;
	      bool isFOElt   = false;
	      bool isFOEll   = false;
	      bool evtFRee   = false;

	      if (isGoodLeptonwIso(11, idx_lt, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodElt = true;
              if (isGoodLeptonwIso(11, idx_ll, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodEll = true;

	      if(applyFOv1Cuts && pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOElt = true;
              if(applyFOv2Cuts && pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOElt = true;
              if(applyFOv3Cuts && pass_electronSelection(idx_lt, electronSelectionFO_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOElt = true;

	      if(applyFOv1Cuts && pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEll = true;
              if(applyFOv2Cuts && pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEll = true;
              if(applyFOv3Cuts && pass_electronSelection(idx_ll, electronSelectionFO_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEll = true;

	      if( isGoodElt && !isGoodEll && isFOEll) evtFRee = true;
	      if( isGoodEll && !isGoodElt && isFOElt) evtFRee = true;
	      if ( isGoodElt && isGoodEll) continue;
	      if (!evtFRee) continue;
	    }

	    // emu case

	    if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
	      int iEl = 0;
	      int iMu = 0;
	      if(hyp_type()[hypIdx] == 2) {
		iEl = hyp_lt_index()[hypIdx];
		iMu = hyp_ll_index()[hypIdx];
	      }
	      if (hyp_type()[hypIdx] == 1) {
		iEl = hyp_ll_index()[hypIdx];
		iMu = hyp_lt_index()[hypIdx];
	      }
	      bool isGoodEl = false;
	      bool isFOEl   = false;
	      bool isGoodMu = false;
	      bool isFOMu   = false;
	      bool evtFRemu = false;

              if (isGoodLeptonwIso(13, iMu, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodMu = true;
              if (isGoodLeptonwIso(11, iEl, applyAlignmentCorrection, removedEtaCutInEndcap) ) isGoodEl = true;
	      if (isFakeableMuon(iMu)) isFOMu = true;
              if(applyFOv1Cuts && pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v1, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEl = true;
              if(applyFOv2Cuts && pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v2, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEl = true;
              if(applyFOv3Cuts && pass_electronSelection(iEl, electronSelectionFO_ssVBTF80_v3, applyAlignmentCorrection, removedEtaCutInEndcap) ) isFOEl = true;
	      if (isGoodMu && !isGoodEl && isFOEl) evtFRemu = true;
	      if (isGoodEl && !isGoodMu && isFOMu) evtFRemu = true;
	      if (isGoodEl && isGoodMu) continue;
	      if (!evtFRemu) continue;

	    }
	  }
	}//	if(applyFOv1Cuts || applyFOv2Cuts || applyFOv3Cuts) 
	
	// require lepton isolation cuts?
	if(applylepIsoCuts) {
	  if(!isGoodHypwIso(hypIdx, applyAlignmentCorrection, removedEtaCutInEndcap))
	    continue;
	}

     if (chargeFlip) {
       bool mischarge = false;

       if (TMath::Abs(hyp_ll_id()[hypIdx]) == 11) {
         int elIndex = hyp_ll_index()[hypIdx];
         if ( isChargeFlip3agree(elIndex)) mischarge = true;
       }

       if (TMath::Abs(hyp_lt_id()[hypIdx]) == 11) {
         int elIndex = hyp_lt_index()[hypIdx];
         if ( isChargeFlip3agree(elIndex)) mischarge = true;
       }
       if (mischarge) continue;
      }


       if (useFlipRateEstimation) { 
         if (hyp_type()[hypIdx] == 3 && useOSleptons) {
           float flip_wt = -99;
           float flip_ll = 0.0;
           float flip_lt = 0.0;
           if( TMath::Abs(hyp_lt_id()[idx_lt])==11 ) { flip_lt = getSingleEleFlipRate(els_p4()[hyp_lt_index()[idx_lt]].pt(), els_p4()[hyp_lt_index()[idx_lt]].eta());}
           if( TMath::Abs(hyp_ll_id()[idx_ll])==11 ) { flip_ll = getSingleEleFlipRate(els_p4()[hyp_ll_index()[idx_ll]].pt(), els_p4()[hyp_ll_index()[idx_ll]].eta());}
           flip_wt = (flip_ll/(1-flip_ll)+flip_lt/(1-flip_lt));
           if (flip_wt <= 0) continue;
           v_goodHyps.push_back(hypIdx);
           v_weights.push_back(flip_wt);
           continue;
        } else {
          continue;
        }
       }
	
	v_goodHyps.push_back(hypIdx);	
	v_weights.push_back(1);	
      }//hypothesis loop
      
      //
      // perform hypothesis disambiguation
      //
      if(v_goodHyps.size() == 0) continue;
      
      if(hypDisamb) {
	int strasbourgDilType = -1;
	unsigned int goodHyp = (unsigned int)eventDilIndexByWeightTTDil08(v_goodHyps, strasbourgDilType, false, false);
        vector<unsigned int>::const_iterator goodHyp_it = find(v_goodHyps.begin(), v_goodHyps.end(), goodHyp);
	if(goodHyp_it == v_goodHyps.end()) {
	  cout << "The weight index does not correspond to the index of the best hypothesis!!!! Something is wrong" 
               << "We will quit" << endl;
          return;
	}

	//clear this vector and put in the goodHyp vector in here so we can then save some space
	//and loop over this vector below. Useful when we're not using the hypDisambiguation
	//get the index of the goodHyp in the vector of goodHyps
	unsigned int goodHyp_idx = goodHyp_it - v_goodHyps.begin();
	//get the weight of the corresponding goodHyp
	float goodHyp_weight = v_weights[goodHyp_idx];
	v_goodHyps.clear();
        v_weights.clear();
	v_goodHyps.push_back(goodHyp);
	v_weights.push_back(goodHyp_weight);
      }//if(hypDisamb)


      //now loop over the good hypotheses. If we require hypothesis disambiguation,
      //we will only have one entry in the vector of good hypothesis
      for(unsigned int i = 0; i < v_goodHyps.size(); i++) {
	
	unsigned int hypIdx = v_goodHyps[i];
	weight = weight*v_weights[i];
	int type = hyp_type()[hypIdx];

	//get the jets passing cuts 
	vector<unsigned int> v_goodJets;
	vector<unsigned int> v_goodJetsNoEtaCut;
	vector<LorentzVector> v_jetP4s;
        double sumet_calo = 0.0;
	if(usecaloJets) {
	  for(unsigned int i = 0; i < jets_p4().size(); i++) {
//	    v_jetP4s.push_back(jets_p4()[i]*jets_cor()[i]); //jets are uncorrected in our ntuplesa
                 LorentzVector jp4 = jets_p4()[i];
                 float jet_cor = jetCorrection(jp4, jet_corrector);
                 v_jetP4s.push_back(jp4 * jet_cor);
           }
	}
	if(usejptJets) {
          for (unsigned int i = 0; i < jpts_p4().size(); i++) {
//            v_jetP4s.push_back(jpts_p4()[i] * jpts_cor()[i]);a
                 LorentzVector jp4 = jpts_p4()[i];
                 float jet_cor = jetCorrection(jp4, jet_corrector);
                 v_jetP4s.push_back(jp4 * jet_cor);
          }
	}
	if(usepfJets) {
//	  v_jetP4s = pfjets_p4();

          for (unsigned int i = 0; i < pfjets_p4().size(); i++)
           {
                LorentzVector jp4 = pfjets_p4()[i];
                float jet_cor = jetCorrection(jp4, jet_corrector);
                v_jetP4s.push_back(jp4 * jet_cor);
           }
        }


	for (unsigned int j = 0; j < v_jetP4s.size(); ++j) {
	  if (isGoodDilHypJet(v_jetP4s.at(j), hypIdx, JETPTCUT, 2.5, 0.4, true)) {
              if (usecaloJets && !passesCaloJetID(v_jetP4s.at(j))) continue;
              if (usejptJets && !passesCaloJetID(v_jetP4s.at(j))) continue;
              if (usepfJets && !passesPFJetID(j)) continue;

	    v_goodJets.push_back(j);
            sumet_calo += v_jetP4s[j].Pt();
           }
	  if(isGoodDilHypJet(v_jetP4s.at(j), hypIdx, JETPTCUT, 9999, 0.4, true)) {

            if (usecaloJets && !passesCaloJetID(v_jetP4s.at(j))) continue;
            if (usejptJets && !passesCaloJetID(v_jetP4s.at(j))) continue;
            if (usepfJets && !passesPFJetID(j)) continue;

	    v_goodJetsNoEtaCut.push_back(j);
          }
	}

	//if we want to veto on nJets, do it here
	if(vetoJets) {
	  if(v_goodJets.size() < 2) continue;
          if (sumet_calo < 60) continue;
	}
	
	//if we want to lower the pt cut, do the selection here
	//this is njet dependent. Only want to lower the pt cut to 20,10 
	//in the nJet > 1 bin
	if(usePtGt2010) {
	  if(max(hyp_lt_p4()[hypIdx].Pt(), hyp_ll_p4()[hypIdx].Pt()) < 20.)
	    continue;
	  if(min(hyp_lt_p4()[hypIdx].Pt(), hyp_ll_p4()[hypIdx].Pt()) < 10.)
	    continue;
	}//if(usePtGt2010)

	
	// MET cut
	string metAlgo;
	if(useCorMET   ) metAlgo  = "CorMET";
	if(usetcMET    ) metAlgo  = "tcMET";
	if(usepfMET    ) metAlgo  = "pfMET"; 	
	pair<float, float> p_met; //met and met phi
	if(usetcMET || usepfMET) {
           if (usetcMET && prefix == "data")
               p_met = getMet("tcMET", hypIdx);
           else
               p_met = getMet(metAlgo, hypIdx);
	  
	  if(p_met.first < 0) {
	    cout << "Something is wrong with the MET. Exiting" << endl;
	    return;
	  }


	  if(vetoMET) {
//            if(p_met.first < 50.)   continue;
            if(hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3) {
              if(p_met.first < 30.) continue;
            }
            if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
              if(p_met.first < 20.)   continue;
            }

	  }

	  if(vetoProjectedMET) {
	    if(hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3) {
	      if(p_met.first < 30.) continue;
	    }
	    if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
	      if(p_met.first < 20.)   continue;
	    }
	    if(v_goodJets.size() < 2 && hyp_p4()[hypIdx].M() < 80.) {
	      if(projectedMET(p_met.first, p_met.second, hypIdx) < 10)
		continue;
	    }
	  }
	  
	} else if(useCorMET) {

	  cout << "THIS HAS BEEN COMMENTED OUT FOR NOW. DOES NOTHING" << endl;
	  /*
	  //this is to agree with Slava
	  float globalJESscaleRescale = 1.;
	  bool muJetClean = true;
	  float metx = met_pat_metCor_hyp(hypIdx)*cos(met_pat_metPhiCor_hyp(hypIdx));
	  float mety = met_pat_metCor_hyp(hypIdx)*sin(met_pat_metPhiCor_hyp(hypIdx));
	  //DOUBLE CHECK THIS !!!!!!!
	  if(usecaloJets) {
	  unsigned int nJused = 0;
	  unsigned int nJ = jets_p4().size();
	  for (unsigned int iJ = 0; iJ < nJ; ++iJ){
	  if (isGoodDilHypJet(jets_cor()[iJ]*jets_p4()[iJ], hypIdx, 0, 2.4, 0.4,muJetClean)){
	  metx -= cms2.jets_cor()[iJ]*cms2.jets_p4()[iJ].x()*(globalJESscaleRescale - 1.); 
	  mety -= cms2.jets_cor()[iJ]*cms2.jets_p4()[iJ].y()*(globalJESscaleRescale - 1.); 
	  nJused++;
	  }
	  }//jet loop
	  }//if use caloJETS
	  
	  p_met = make_pair(sqrt(metx*metx + mety*mety), atan2(mety, metx));
	  if(vetoMET) {
	  if(hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3) {
	  if(p_met.first < 30.) continue;
	  }
	  if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
	  if(p_met.first < 20.)   continue;
	  }
	  }

	  if(vetoProjectedMET) {
	  if(v_goodJets.size() < 2 && hyp_p4()[hypIdx].M() < 80.) {
	  if(projectedMET(p_met.first, p_met.second, hypIdx) < 10)
	  continue;
	  }
	  if(hyp_type()[hypIdx] == 0 || hyp_type()[hypIdx] == 3) {
	  if(p_met.first < 30.) continue;
	  }
	  if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
	  if(p_met.first < 20.)   continue;
	  }
	  }
	  */
	}//if vetoing on corrected caloMet
       
//        float isolt = -99.;
//        float isoll = -99.;
        sumet_calo += hyp_lt_p4()[hypIdx].Pt();
        sumet_calo += hyp_ll_p4()[hypIdx].Pt();
	nGoodHyps[hyp_type()[hypIdx]]++;	
	if(prefix == "data") {
	  cout <<"Final_" << evt_run() << "_" << evt_lumiBlock() << "_" << evt_event() << "_Type_" << hyp_type()[hypIdx] << "_Njet_" << v_goodJets.size() <<"_MET_" <<p_met.first << "_SumJetPt_" << sumet_calo << "_Mll_" << hyp_p4()[hypIdx].mass() << "_Pt1_" << hyp_lt_p4()[hypIdx].Pt() << "_Eta1_" << hyp_lt_p4()[hypIdx].Eta() << "_Pt2_" << hyp_ll_p4()[hypIdx].Pt() << "_Eta2_" << hyp_ll_p4()[hypIdx].Eta() << "_ID1_" <<  hyp_lt_id()[hypIdx] << "_ID2_" << hyp_ll_id()[hypIdx] << endl;
/*
	  string type = "";
	  if(hyp_type()[hypIdx] == 0) {
	    type = "mumu";
            int idx_lt = hyp_lt_index()[hypIdx];
            int idx_ll = hyp_ll_index()[hypIdx];
            cout << " mumu_Iso_" << muonIsoValue(idx_lt) << "_" << muonIsoValue(idx_ll) << "_pt_" << TMath::Max(hyp_lt_p4()[hypIdx].Pt(),hyp_ll_p4()[hypIdx].Pt()) << "_" << TMath::Min(hyp_lt_p4()[hypIdx].Pt(),hyp_ll_p4()[hypIdx].Pt()) << "_Mll_" << hyp_p4()[hypIdx].mass();
	  } else if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
	    type = "emu";
            int iEl = 0;
            int iMu = 0;
            if(hyp_type()[hypIdx] == 2) {
              iEl = hyp_lt_index()[hypIdx];
              iMu = hyp_ll_index()[hypIdx];
            }
            if (hyp_type()[hypIdx] == 1) {
              iEl = hyp_ll_index()[hypIdx];
              iMu = hyp_lt_index()[hypIdx];
            }
           cout << " emu_eIso_" <<  electronIsolation_rel(iEl, true) << "_muIsoi_" << muonIsoValue(iMu) << "_ePt_" << els_p4()[iEl].pt() << "_mPt_" << mus_p4()[iMu].pt() << "_Mll_" << hyp_p4()[hypIdx].mass();
           cout << " emu_dEtaIn_" <<  els_dEtaIn()[iEl] << "_dPhiIn_" << els_dPhiIn()[iEl] << "_hOverE_" << els_hOverE()[iEl] << "_sigmaIEtaIEta_" << els_sigmaIEtaIEta()[iEl];
           cout << " emu_EleCharge_" << els_trk_charge().at(iEl) << "_GSF_" << trks_charge().at(els_trkidx().at(iEl)) << "_Pix_" << els_sccharge().at(iEl);
          } else {
	    type = "ee";
            int idx_lt = hyp_lt_index()[hypIdx];
            int idx_ll = hyp_ll_index()[hypIdx];
            cout << " ee_e1Iso_" << electronIsolation_rel(idx_lt, true) << "_e2Iso_" << electronIsolation_rel(idx_ll, true) << "_e1Pt1_" << els_p4()[idx_lt].pt() <<"_e1Pt2_" <<  els_p4()[idx_ll].pt() << "_Mll_" << hyp_p4()[hypIdx].mass();
            cout << " ee1_dEtaIn_" <<  els_dEtaIn()[idx_lt] << "_dPhiIn_" << els_dPhiIn()[idx_lt] << "_hOverE_" << els_hOverE()[idx_lt] << "_sigmaIEtaIEta_" << els_sigmaIEtaIEta()[idx_lt];
            cout << " ee2_dEtaIn_" <<  els_dEtaIn()[idx_ll] << "_dPhiIn_" << els_dPhiIn()[idx_ll] << "_hOverE_" << els_hOverE()[idx_ll] << "_sigmaIEtaIEta_" << els_sigmaIEtaIEta()[idx_ll];
            cout << " ee1_Charge_" << els_trk_charge().at(idx_lt) << "_GSF_" << trks_charge().at(els_trkidx().at(idx_lt)) << "_Pix_" << els_sccharge().at(idx_lt);
            cout << " ee2_Charge_" << els_trk_charge().at(idx_ll) << "_GSF_" << trks_charge().at(els_trkidx().at(idx_ll)) << "_Pix_" << els_sccharge().at(idx_ll);
          }
          cout <<"_"<<endl;
*/
	}
	FillHistograms(hypIdx, v_goodJets, v_goodJetsNoEtaCut, p_met, weight, prefix);	
	
      }//hypothesis loop
      if(v_goodHyps.size() > 0)
	nGoodEvts++;
    }//event loop
  }//file loop
  cout << "Number of good ee Hyps:  " << nGoodHyps[3] << endl;
  cout << "Number of good mumu Hyps: " << nGoodHyps[0] << endl;
  cout << "Number of good emu Hyps: " << nGoodHyps[1] + nGoodHyps[2] << endl;
  cout << "Number of good evts before overlap: " << nGoodEvts << endl; 
  cout << "Number of good evts after overlap: " << nGoodHyps[0] + nGoodHyps[1] + nGoodHyps[2] + nGoodHyps[3] << endl; 
  cout << "sum fr: " << sumfr << endl;
  EndJob();
  cout << "Total number of events after overlap " << nEventsTotal << endl; 
  cout << "Total number of events before overlap " << nEventsChain << endl; 
  if ( nEventsChain != nEventsTotal ) {
//    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }
    
    
}

void FillHistograms(const unsigned int hypIdx, const vector<unsigned int> v_jets, const vector<unsigned int> v_jetsNoEtaCut,
		    const pair<float, float> p_met, const float weight, std::string prefix) {
  
  int type   = hyp_type()[hypIdx];
  int lt_id  = hyp_lt_id()[hypIdx];
  int lt_idx = hyp_lt_index()[hypIdx];
  LorentzVector lt_p4 = hyp_lt_p4()[hypIdx];	
  
  int ll_id  = hyp_ll_id()[hypIdx];
  int ll_idx = hyp_ll_index()[hypIdx];
  LorentzVector ll_p4 = hyp_ll_p4()[hypIdx];
  LorentzVector hypp4 = hyp_p4()[hypIdx];	
  

  //make a vector of corrected jets
  vector<LorentzVector> v_jets_p4;
  for(unsigned int i = 0; i < v_jets.size(); i++) {
    if(usecaloJets)
      v_jets_p4.push_back(jets_p4()[v_jets[i]]*jets_cor()[v_jets[i]]);
    else if(usejptJets)
      v_jets_p4.push_back(jpts_p4()[v_jets[i]]);
    else if(usepfJets)
      v_jets_p4.push_back(pfjets_p4()[v_jets[i]]);
  }
			  

  //sort the jets
  std::sort(v_jets_p4.begin(), v_jets_p4.end(), sortByPt);
  
  //instead of duplicating all the fills twice (once for the correct type and then for the 
  //histogram that combines them all, we can loop over a vector of ints of size 2 whose first 
  //entry is for index of the histogram of all the channels. From Histograms.cc (bookhistos):
  //suffixall[0] = "ee";
  //suffixall[1] = "mm";
  //suffixall[2] = "em";
  //suffixall[3] = "all";
  
  vector<int> v_type;
  v_type.push_back(3);	
  if(type == 1 || type == 2) //emu
    v_type.push_back(2);
  if(type == 0) //mumu
    v_type.push_back(1);
  if(type == 3) //ee
    v_type.push_back(0);	


  unsigned int jetBin = v_jets.size();
  if(jetBin > 3)
    jetBin = 4;

  int tttype = -99;
  int lttype = -99;
  int lltype = -99;

  // Classify the different types based on Truth  
  if (prefix != "data") {
    //To distinguish WW (=1), WO (=2), and OO (=3) 
    tttype = ttbarconstituents(hypIdx);
    // Semileptonic
    lttype = leptonIsFromW(hyp_lt_index()[hypIdx], hyp_lt_id()[hypIdx] );
    lltype = leptonIsFromW(hyp_ll_index()[hypIdx], hyp_ll_id()[hypIdx] );
  }

  
  for(unsigned int i = 0; i < v_type.size(); i++) {
    
    //which channel?
    unsigned int ch  = v_type.at(i);
      
    
    hnJet[ch]                         ->Fill(jetBin,        weight);

    if (tttype == 1) {
      hnJetWW[ch]->Fill(jetBin, weight);
    } else if (tttype == 2) {
      hnJetWO[ch]->Fill(jetBin, weight);
      if ( lttype == -1 || lttype == -2 || lltype == -1 || lltype == -2) {
	hnJetWOSemilep[ch]->Fill(jetBin, weight);
      } else {
	hnJetWOOther[ch]->Fill(jetBin, weight);
      }
    } else {
      hnJetOO[ch]->Fill(jetBin, weight);
    }

    if(inZmassWindow(hypp4.mass())) {
      hnJetinZwindow[ch]             ->Fill(jetBin,        weight);
      hdilMassTightWindow[ch][jetBin]->Fill(hypp4.mass(),  weight); 
    }
    else
      hnJetoutZwindow[ch]             ->Fill(jetBin,        weight);	
    
    if(abs(lt_id)==11) {
      helePt[ch][jetBin]              ->Fill(lt_p4.Pt(),    weight);
      helePhi[ch][jetBin]             ->Fill(lt_p4.Phi(),   weight); 
      heleEta[ch][jetBin]             ->Fill(lt_p4.Eta(),   weight);
// -s      float relIso = electronIsolation_relsusy_cand1(lt_idx, true);
      float relIso = 1;
      float sumIso = relIso*max(els_p4()[lt_idx].Pt(), (float)20.);
      helSumIso[ch][jetBin]           ->Fill(sumIso, weight); 
      helRelIso[ch][jetBin]           ->Fill(relIso, weight); 	
      helRelIsoTrack[ch][jetBin]      ->Fill(els_tkJuraIso()[lt_idx], weight);
      double caloIso = els_hcalIso()[lt_idx];
      if (fabs(els_etaSC()[lt_idx]) > 1.479) caloIso += els_ecalIso()[lt_idx];
      if (fabs(els_etaSC()[lt_idx]) <= 1.479) caloIso += max(0., (els_ecalIso()[lt_idx] -1.));
      helRelIsoCalo[ch][jetBin]       ->Fill(caloIso, weight                    );
      heled0BS[ch][jetBin]            ->Fill(els_d0corr()[lt_idx], weight       );
      heled0PV[ch][jetBin]            ->Fill(getd0wrtPV(els_trk_p4()[lt_idx], els_d0()[lt_idx]), weight);
      heleEmaxOE5x5[ch][jetBin]       ->Fill(els_eMax()[lt_idx]/els_e5x5()[lt_idx], weight);

      helIsoTrack[ch][jetBin]         ->Fill(els_tkJuraIso()[lt_idx], weight);
      helIsoEcal[ch][jetBin]          ->Fill(els_ecalIso()[lt_idx],   weight);
      helIsoHcal[ch][jetBin]          ->Fill(els_hcalIso()[lt_idx],   weight);

      if (fabs(els_etaSC()[lt_idx]) < 1.479) {
	helIsoTrackb[ch][jetBin]         ->Fill(els_tkJuraIso()[lt_idx], weight);
	helIsoEcalb[ch][jetBin]          ->Fill(els_ecalIso()[lt_idx],   weight);
	helIsoHcalb[ch][jetBin]          ->Fill(els_hcalIso()[lt_idx],   weight);
      } else {
	helIsoTracke[ch][jetBin]         ->Fill(els_tkJuraIso()[lt_idx], weight);
	helIsoEcale[ch][jetBin]          ->Fill(els_ecalIso()[lt_idx],   weight);
	helIsoHcale[ch][jetBin]          ->Fill(els_hcalIso()[lt_idx],   weight);
      }
      
    }//if(abs(lt_id)==11) {
    if(abs(ll_id)==11) {
      helePt[ch][jetBin]              ->Fill(ll_p4.Pt(),    weight);	
      helePhi[ch][jetBin]             ->Fill(ll_p4.Phi(),   weight); 
      heleEta[ch][jetBin]             ->Fill(ll_p4.Eta(),   weight);
//-s      float relIso = electronIsolation_relsusy_cand1(ll_idx, true);
      float relIso = 1;
      float sumIso = relIso*max(els_p4()[ll_idx].Pt(), (float)20.);
      helSumIso[ch][jetBin]           ->Fill(sumIso, weight); 
      helRelIso[ch][jetBin]           ->Fill(relIso, weight); 	
      helRelIsoTrack[ch][jetBin]      ->Fill(els_tkJuraIso()[ll_idx], weight);
      double caloIso = els_hcalIso()[lt_idx];
      if (fabs(els_etaSC()[ll_idx]) > 1.479) caloIso += els_ecalIso()[ll_idx];
      if (fabs(els_etaSC()[ll_idx]) <= 1.479) caloIso += max(0., (els_ecalIso()[ll_idx] -1.));
      helRelIsoCalo[ch][jetBin]       ->Fill(caloIso, weight);	
      heled0BS[ch][jetBin]            ->Fill(els_d0corr()[ll_idx], weight);
      heled0PV[ch][jetBin]            ->Fill(getd0wrtPV(els_trk_p4()[ll_idx], els_d0()[ll_idx]), weight);
      heleEmaxOE5x5[ch][jetBin]       ->Fill(els_eMax()[ll_idx]/els_e5x5()[ll_idx], weight);

      if (fabs(els_etaSC()[ll_idx]) < 1.479) {
	helIsoTrackb[ch][jetBin]         ->Fill(els_tkJuraIso()[ll_idx], weight);
	helIsoEcalb[ch][jetBin]          ->Fill(els_ecalIso()[ll_idx],   weight);
	helIsoHcalb[ch][jetBin]          ->Fill(els_hcalIso()[ll_idx],   weight);
      } else {
	helIsoTracke[ch][jetBin]         ->Fill(els_tkJuraIso()[ll_idx], weight);
	helIsoEcale[ch][jetBin]          ->Fill(els_ecalIso()[ll_idx],   weight);
	helIsoHcale[ch][jetBin]          ->Fill(els_hcalIso()[ll_idx],   weight);
      }

    }//if(abs(ll_id)==11) {
    if(abs(lt_id) == 13) {
      hmuPt[ch][jetBin]               ->Fill(lt_p4.Pt(),    weight); 
      hmuPtFromSilicon[ch][jetBin]    ->Fill(mus_trk_p4()[lt_idx].Pt(), weight);
      hmuPhi[ch][jetBin]              ->Fill(lt_p4.Phi(),   weight);	
      hmuEta[ch][jetBin]              ->Fill(lt_p4.Eta(),   weight);
      float relIso = muonIsoValue(lt_idx);
      float sumIso = relIso*max(mus_p4()[lt_idx].Pt(), (float)20);
      hmuSumIso[ch][jetBin]           ->Fill(min(sumIso, (float)24.99), weight);	
      hmuRelIso[ch][jetBin]           ->Fill(relIso, weight);
      hmuRelIsoTrack[ch][jetBin]      ->Fill(mus_iso03_sumPt()[lt_idx]/max(mus_p4()[lt_idx].Pt(), (float)20.), weight);
      hmuRelIsoCalo[ch][jetBin]       ->Fill((mus_iso03_hadEt()[lt_idx]+mus_iso03_emEt()[lt_idx])/max(mus_p4()[lt_idx].Pt(), (float)20.), weight);	
      hmud0BS[ch][jetBin]             ->Fill(mus_d0corr()[lt_idx], weight);
      hmud0PV[ch][jetBin]             ->Fill(getd0wrtPV(mus_trk_p4()[lt_idx], mus_d0()[lt_idx]), weight);
    }
    if(abs(ll_id) == 13) {
      hmuPt[ch][jetBin]               ->Fill(ll_p4.Pt(),    weight); 
      hmuPtFromSilicon[ch][jetBin]    ->Fill(mus_trk_p4()[ll_idx].Pt(), weight);
      hmuPhi[ch][jetBin]              ->Fill(ll_p4.Phi(), weight);	
      hmuEta[ch][jetBin]              ->Fill(ll_p4.Eta(),   weight);
      float relIso = muonIsoValue(ll_idx);
      float sumIso = relIso*max(mus_p4()[ll_idx].Pt(), (float)20);
      hmuSumIso[ch][jetBin]           ->Fill(min(sumIso, (float)24.99), weight);	
      hmuRelIso[ch][jetBin]           ->Fill(relIso, weight);
      hmuRelIsoTrack[ch][jetBin]      ->Fill(mus_iso03_sumPt()[ll_idx]/max(mus_p4()[ll_idx].Pt(), (float)20), weight);
      hmuRelIsoCalo[ch][jetBin]       ->Fill((mus_iso03_hadEt()[ll_idx]+mus_iso03_emEt()[ll_idx])/max(mus_p4()[ll_idx].Pt(), (float)20.), weight);	
      hmud0BS[ch][jetBin]             ->Fill(mus_d0corr()[ll_idx], weight);
      hmud0PV[ch][jetBin]             ->Fill(getd0wrtPV(mus_trk_p4()[ll_idx], mus_d0()[ll_idx]), weight);
    }
    
    hminLepPt[ch][jetBin]             ->Fill(min(lt_p4.Pt(), ll_p4.Pt()), weight); 
    hmaxLepPt[ch][jetBin]             ->Fill(max(lt_p4.Pt(), ll_p4.Pt()), weight);
    
    double dphi = fabs(lt_p4.Phi() - ll_p4.Phi());
    if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
    
    hdphiLep[ch][jetBin]              ->Fill(dphi,                    weight);
    hdilMass[ch][jetBin]              ->Fill(hypp4.mass(),            weight);	
    hdilPt[ch][jetBin]                ->Fill(hypp4.Pt(),              weight);
    hmet[ch][jetBin]                  ->Fill(evt_metMuonJESCorr(),    weight); 
    hmetPhi[ch][jetBin]               ->Fill(evt_metMuonJESCorrPhi(), weight); 
    hpfmet[ch][jetBin]                ->Fill(evt_pfmet(),             weight);
    hpfmetPhi[ch][jetBin]             ->Fill(evt_pfmetPhi(),          weight); 
/*
//    double tcmet = evt_tcmet();
//    double tcmetPhi = evt_tcmetPhi();
    std::pair<float, float> p_tcmet;
    if(prefix == "data") 
         p_tcmet = getMet("tcMET", hypIdx);
    else 
         p_tcmet = getMet("tcMET35X", hypIdx);
    
    htcmet[ch][jetBin]                ->Fill(p_tcmet.first,           weight); 
    htcmetPhi[ch][jetBin]             ->Fill(p_tcmet.second,        weight); 
    
    //projected MET
    hprojmet[ch][jetBin]              ->Fill(projectedMET(evt_metMuonJESCorr(),
							  evt_metMuonJESCorrPhi(),
							  hypIdx),    weight);  
    hprojpfmet[ch][jetBin]            ->Fill(projectedMET(evt_pfmet(), evt_pfmetPhi(),
							  hypIdx), weight); 
    hprojtcmet[ch][jetBin]            ->Fill(projectedMET(p_tcmet.first, p_tcmet.second,
							  hypIdx), weight);


    hmetVsDilepPt[ch][jetBin]         ->Fill(hypp4.Pt(), evt_metMuonCorr(), weight); 
    hmetOverPtVsDphi[ch][jetBin]      ->Fill(dphi, evt_metMuonCorr()/hypp4.Pt(), weight); 

    hpfmetVsDilepPt[ch][jetBin]       ->Fill(hypp4.Pt(), evt_pfmet(),       weight); 
    hpfmetOverPtVsDphi[ch][jetBin]    ->Fill(dphi, evt_pfmet()/hypp4.Pt(),  weight);  

    htcmetVsDilepPt[ch][jetBin]       ->Fill(hypp4.Pt(), p_tcmet.first,       weight);
    htcmetOverPtVsDphi[ch][jetBin]    ->Fill(dphi, p_tcmet.first/hypp4.Pt(),  weight);

*/
     hdphillvsmll[ch][jetBin]          ->Fill(hypp4.mass(), dphi, weight);

// Meff
    float meff = 0.;
    if(v_jets_p4.size() > 0) {
      for(unsigned int i = 0; i < v_jets_p4.size(); i++) {
        meff += v_jets_p4[i].Pt();
      }
    }
    meff += ll_p4.Pt();
    meff += lt_p4.Pt();

    hnMeff[ch]->Fill(meff, weight);
    
    if(v_jets_p4.size() > 0) {
      hptJet1[ch][jetBin]             ->Fill(v_jets_p4[0].Pt(),   weight);
      hetaJet1[ch][jetBin]            ->Fill(v_jets_p4[0].Eta(),  weight); 
    }
    if(v_jets_p4.size() > 1) {
      hptJet2[ch][jetBin]             ->Fill(v_jets_p4[1].Pt(),   weight); 
      hetaJet2[ch][jetBin]            ->Fill(v_jets_p4[1].Eta(),  weight); 
    }
    if(v_jets_p4.size() > 2) {
      hptJet3[ch][jetBin]             ->Fill(v_jets_p4[2].Pt(),   weight); 
      hetaJet3[ch][jetBin]            ->Fill(v_jets_p4[2].Eta(),  weight); 
    }
    if(v_jets_p4.size() > 3) {
      hptJet4[ch][jetBin]             ->Fill(v_jets_p4[3].Pt(),   weight); 
      hetaJet4[ch][jetBin]            ->Fill(v_jets_p4[3].Eta(),  weight); 
    }
    
    //fill MT2 stuff - use jets with no eta cut
    v_jets_p4.clear();	
    for(unsigned int i = 0; i < v_jetsNoEtaCut.size(); i++) {
      if(usecaloJets)
	v_jets_p4.push_back(jets_p4()[v_jetsNoEtaCut[i]]*jets_cor()[v_jetsNoEtaCut[i]]);
      else if(usejptJets) {
	v_jets_p4.push_back(jpts_p4()[v_jetsNoEtaCut[i]]);	
      } else if(usepfJets) {
	v_jets_p4.push_back(pfjets_p4()[v_jetsNoEtaCut[i]]);	
      }
    }
    hmt2[ch][jetBin]                   ->Fill(MT2(p_met.first, p_met.second, lt_p4, ll_p4, 0.0, false), weight);	
    if(v_jets.size() > 1) 
      hmt2J[ch][jetBin]                ->Fill(MT2J(p_met.first, p_met.second, lt_p4, ll_p4, v_jets_p4, 0.0, BISECT, false),weight);	
  }
  
    
}

// *****************************************************************
//get the FR weight
// *****************************************************************
double getFRWeight(const int hypIdx, string elFRversion) {

  bool isGoodMut = false;
  bool isGoodMul = false;
  bool isFOMut   = false;
  bool isFOMul   = false;

  if(hyp_type()[hypIdx] == 0) {
    
    unsigned int iMut = hyp_lt_index()[hypIdx];
    unsigned int iMul = hyp_ll_index()[hypIdx];
    
    if(isGoodLeptonwIso(13, iMut, applyAlignmentCorrection, removedEtaCutInEndcap) )
      isGoodMut = true;
    
    if(isGoodLeptonwIso(13, iMul, applyAlignmentCorrection, removedEtaCutInEndcap) )
      isGoodMul = true;
    
    if(isFakeableMuon(iMut))
      isFOMut = true;
       
    if(isFakeableMuon(iMul))
      isFOMul = true;
    
    if(!isFOMut || !isFOMul)
      return -9999.;
    
 //   SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", "muFR15u");
  //  SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FakeRates31May.root", "muFR15u");
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FakeRates7July.root", "muFR15u");
 //   SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", "muFR15u");
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30FR_SSAug12.root", "iso10_muFR15u");
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/jmtFR_SSAug12.root", "iso10_muFR15u");

// Aug31st
//       SimpleFakeRate fake("../data/SSFakeRates31August.root", "iso10_muFR15u");
     SimpleFakeRate fake("../data/qcd30_SSFakeRates31August.root", "iso10_muFR15u");

    if (estimateDoubleFakes) {

      if( isGoodMut || isGoodMul)
	return -9999.;
      
      double FRMut = fake.getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
      double FRMul = fake.getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
      cout << "mm, FRlt, FRll and FR/(1-FR)FR/(1-FR) " << FRMut << "  " << FRMul << "  " << (FRMut/(1-FRMut))*(FRMul/(1-FRMul)) <<  endl;
      return (FRMut/(1-FRMut))*(FRMul/(1-FRMul));
    } else if (estimateSingleFakes) {
      
      //need one to be a Numerator lepton, and the other to be FO but not num
      if( isGoodMut && !isGoodMul && isFOMul) {
        double FR = fake.getFR(mus_p4()[iMul].pt(), mus_p4()[iMul].eta());
        cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }

      //check the other muon
      if( isGoodMul && !isGoodMut && isFOMut) {
        double FR = fake.getFR(mus_p4()[iMut].pt(), mus_p4()[iMut].eta());
        cout << "mm, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
    }
    return -9999.;
  }//mumu case

  //now we do the ee case
  if(hyp_type()[hypIdx] == 3) {
	  
    unsigned int iElt = hyp_lt_index()[hypIdx];
    unsigned int iEll = hyp_ll_index()[hypIdx];
	  
    bool isGoodElt = false;
    bool isGoodEll = false;
    bool isFOElt   = false;
    bool isFOEll   = false;

    if(isGoodLeptonwIso(11, iElt, applyAlignmentCorrection, removedEtaCutInEndcap))
      isGoodElt = true;
       
    if(isGoodLeptonwIso(11, iEll, applyAlignmentCorrection, removedEtaCutInEndcap))
      isGoodEll = true;

    if(isFakeableElectron(iElt,elFRversion, applyAlignmentCorrection, removedEtaCutInEndcap))
      isFOElt = true;
    
    if(isFakeableElectron(iEll,elFRversion, applyAlignmentCorrection, removedEtaCutInEndcap))
      isFOEll = true;

    if( !isFOElt || !isFOEll)
      return -9999.;
    
    //double FRElt = elFakeProb(iElt, elFRversion);
    //double FREll = elFakeProb(iEll, elFRversion);
    char* url;
    url = new char[elFRversion.length() + 1];
    strcpy(url, elFRversion.c_str()); 
  //  SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FakeRates31May.root", url);
  //  SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FakeRates7July.root", url);
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", url);
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", url);
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30FR_SSAug12.root", url);
//    SimpleFakeRate fake("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/jmtFR_SSAug12.root", url);

// Aug31 

//      SimpleFakeRate fake("../data/SSFakeRates31August.root", url);
      SimpleFakeRate fake("../data/qcd30_SSFakeRates31August.root", url);

    delete [] url;
    if(estimateDoubleFakes) {

      if( isGoodElt || isGoodEll)
	return -9999.;
      
      double FRElt = fake.getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
      double FREll = fake.getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
      sumfr = sumfr +  (FRElt/(1-FRElt))*(FREll/(1-FREll));

      cout << "ee, FRlt, FRll, and FR/(1-FR)xFR/(1-FR) " << FRElt << "  " << FREll << "  " << (FRElt/(1-FRElt))*(FREll/(1-FREll)) << endl;
      return (FRElt/(1-FRElt))*(FREll/(1-FREll));

    } else if (estimateSingleFakes) {
      
      if(isGoodElt && !isGoodEll && isFOEll) {
        double FR = fake.getFR(els_p4()[iEll].pt(), els_p4()[iEll].eta());
        cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      //check the other electron
      if(isGoodEll && !isGoodElt && isFOElt) {
        double FR = fake.getFR(els_p4()[iElt].pt(), els_p4()[iElt].eta());
        cout << "ee, FR and FR/(1-FR) " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
    }
    return -9999.;
  }//ee case

  if(hyp_type()[hypIdx] == 1 || hyp_type()[hypIdx] == 2) {
     int iEl = 0;
    int iMu = 0;
    if(hyp_type()[hypIdx] == 2) {
      iEl = hyp_lt_index()[hypIdx];
      iMu = hyp_ll_index()[hypIdx];
    } 
    if (hyp_type()[hypIdx] == 1) {
      iEl = hyp_ll_index()[hypIdx];
      iMu = hyp_lt_index()[hypIdx];
    } 
    bool isGoodEl = false;
    bool isFOEl   = false;
    bool isGoodMu = false;
    bool isFOMu   = false;
    if(isGoodLeptonwIso(11, iEl, applyAlignmentCorrection, removedEtaCutInEndcap))
      isGoodEl = true;
    
    if(isGoodLeptonwIso(13, iMu, applyAlignmentCorrection, removedEtaCutInEndcap))
      isGoodMu = true;

    if(isFakeableElectron(iEl,elFRversion, applyAlignmentCorrection, removedEtaCutInEndcap))
      isFOEl = true;
    
    if(isFakeableMuon(iMu))
      isFOMu = true;
    

    if(!isFOMu || !isFOEl)
      return -9999.;

    char* url;
    url = new char[elFRversion.length() + 1];
    strcpy(url, elFRversion.c_str());

//    SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", "muFR15u");
 //   SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", "muFR15u");
//    SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30FR_SSAug12.root", "iso10_muFR15u");
//    SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/jmtFR_SSAug12.root", "iso10_muFR15u");
//    SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/SSFakeRates31August.root", "iso10_muFR15u");
//    SimpleFakeRate mufr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30_SSFakeRates31August.root", "iso10_muFR15u");
//    SimpleFakeRate elfr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/FR_qcd30_SSJul26.root", url);
//    SimpleFakeRate elfr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30FR_SSAug12.root", url);
//    SimpleFakeRate elfr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/jmtFR_SSAug12.root", url);
//    SimpleFakeRate elfr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/SSFakeRates31August.root", url);
//    SimpleFakeRate elfr("/home/users/spadhi/CMS/TAS/Jul2010/dileptons/data/qcd30_SSFakeRates31August.root", url);

// Aug31st

//     SimpleFakeRate mufr("../data/SSFakeRates31August.root", "iso10_muFR15u");  
//     SimpleFakeRate elfr("../data/SSFakeRates31August.root", url); 

     SimpleFakeRate mufr("../data/qcd30_SSFakeRates31August.root", "iso10_muFR15u");
     SimpleFakeRate elfr("../data/qcd30_SSFakeRates31August.root", url);


    delete [] url;

    if (estimateDoubleFakes) {
      if(isGoodMu || isGoodEl)
	return -9999.;
      double FRMu = mufr.getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
      double FREl = elfr.getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());

      cout << "emu, FRMu, FREl, FR/(1-FR)xFR/(1-FR) " << FRMu << "  " << FREl << "  " << (FRMu/(1-FRMu))*(FREl/(1-FREl)) << endl;  
      return (FRMu/(1-FRMu))*(FREl/(1-FREl));
      
    } else if (estimateSingleFakes) {

      //need one to be a numerator lepton and the other to be a FO
      if(isGoodMu && !isGoodEl && isFOEl) {
        double FR = elfr.getFR(els_p4()[iEl].pt(), els_p4()[iEl].eta());
        cout << "emu, el FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
      if(isGoodEl && !isGoodMu && isFOMu) {
        double FR = mufr.getFR(mus_p4()[iMu].pt(), mus_p4()[iMu].eta());
        cout << "emu, mu FR, FR/(1-FR): " << FR << ", " << FR/(1-FR) << endl;
        return FR/(1-FR);
      }
    }
    return -9999.;
  } //emu case
  return -9999.;
}

// *****************************************************************
//get the bin value
// *****************************************************************
float GetValueTH2F(Float_t x, Float_t y, const TH2F* h) {
  Int_t binx = h->GetXaxis()->FindBin(x);
  Int_t biny = h->GetYaxis()->FindBin(y);
  return h->GetBinContent(binx, biny);
}


void EndJob() {

  for(unsigned int i = 0; i < 4; i++) {
    for(unsigned int j = 2; j < 5; j++) { //loop to fill the geq2 histos
      
      helePt[i][5]->Add(helePt[i][j]);
      hmuPt[i][5]->Add(hmuPt[i][j]);
      hmuPtFromSilicon[i][5]->Add(hmuPtFromSilicon[i][j]);
      hminLepPt[i][5]->Add(hminLepPt[i][j]);
      hmaxLepPt[i][5]->Add(hmaxLepPt[i][j]);
      helePhi[i][5]->Add(helePhi[i][j]);
      hmuPhi[i][5]->Add(hmuPhi[i][j]);
      hdphiLep[i][5]->Add(hdphiLep[i][j]);
      heleEta[i][5]->Add(heleEta[i][j]);
      hmuEta[i][5]->Add(hmuEta[i][j]);
      heled0BS[i][5]->Add(heled0BS[i][j]);
      heled0PV[i][5]->Add(heled0PV[i][j]);
      hmud0BS[i][5]->Add(heled0PV[i][j]);
      hmud0PV[i][5]->Add(heled0PV[i][j]);
      heleEmaxOE5x5[i][5]->Add(heleEmaxOE5x5[i][j]);


      helIsoTrack[i][5]->Add(helIsoTrack[i][j]);
      helIsoEcal[i][5]->Add(helIsoEcal[i][j]);
      helIsoHcal[i][5]->Add(helIsoHcal[i][j]);

      helIsoTrackb[i][5]->Add(helIsoTrackb[i][j]);
      helIsoEcalb[i][5]->Add(helIsoEcalb[i][j]);
      helIsoHcalb[i][5]->Add(helIsoHcalb[i][j]);

      helIsoTracke[i][5]->Add(helIsoTracke[i][j]);
      helIsoEcale[i][5]->Add(helIsoEcale[i][j]);
      helIsoHcale[i][5]->Add(helIsoHcale[i][j]);
           


      hdilMass[i][5]->Add(hdilMass[i][j]);
      hdilMassTightWindow[i][5]->Add(hdilMassTightWindow[i][j]);
      hdilPt[i][5]->Add(hdilPt[i][j]);
      hmet[i][5]->Add(hmet[i][j]);
      hmetPhi[i][5]->Add(hmetPhi[i][j]);
      hpfmet[i][5]->Add(hpfmet[i][j]);
      hpfmetPhi[i][5]->Add(hpfmetPhi[i][j]); 
      htcmet[i][5]->Add(htcmet[i][j]);
      htcmetPhi[i][5]->Add(htcmetPhi[i][j]);
      hptJet1[i][5]->Add(hptJet1[i][j]);
      hptJet2[i][5]->Add(hptJet2[i][j]);
      hptJet3[i][5]->Add(hptJet3[i][j]);
      hptJet4[i][5]->Add(hptJet4[i][j]);
      hetaJet1[i][5]->Add(hetaJet1[i][j]);
      hetaJet2[i][5]->Add(hetaJet2[i][j]);
      hetaJet3[i][5]->Add(hetaJet3[i][j]);
      hetaJet4[i][5]->Add(hetaJet4[i][j]);
      heleSumPt[i][5]->Add(heleSumPt[i][j]);
      hmuSumPt[i][5]->Add(hmuSumPt[i][j]);
      
      hmuSumIso[i][5]->Add(hmuSumIso[i][j]);
      helSumIso[i][5]->Add(helSumIso[i][j]);
      hmuRelIso[i][5]->Add(hmuRelIso[i][j]);
      helRelIso[i][5]->Add(helRelIso[i][j]);
      hmuRelIsoTrack[i][5]->Add(hmuRelIsoTrack[i][j]);
      helRelIsoTrack[i][5]->Add(helRelIsoTrack[i][j]);
      hmuRelIsoCalo[i][5]->Add(hmuRelIsoCalo[i][j]);
      helRelIsoCalo[i][5]->Add(helRelIsoCalo[i][j]);
      
      hmt2[i][5]->Add(hmt2[i][j]);
      hmt2J[i][5]->Add(hmt2J[i][j]);
      
    }
   
    for(unsigned int j = 0; j < 5; j++) { //loop to fill the all jet bins histos
       
      helePt[i][6]->Add(helePt[i][j]);
      hmuPt[i][6]->Add(hmuPt[i][j]);
      hmuPtFromSilicon[i][6]->Add(hmuPtFromSilicon[i][j]);
      hminLepPt[i][6]->Add(hminLepPt[i][j]);
      hmaxLepPt[i][6]->Add(hmaxLepPt[i][j]);
      helePhi[i][6]->Add(helePhi[i][j]);
      hmuPhi[i][6]->Add(hmuPhi[i][j]);
      hdphiLep[i][6]->Add(hdphiLep[i][j]);
      heleEta[i][6]->Add(heleEta[i][j]);
      hmuEta[i][6]->Add(hmuEta[i][j]);
      heled0BS[i][6]->Add(heled0BS[i][j]);
      heled0PV[i][6]->Add(heled0PV[i][j]);
      hmud0BS[i][6]->Add(heled0PV[i][j]);
      hmud0PV[i][6]->Add(heled0PV[i][j]);
      heleEmaxOE5x5[i][6]->Add(heleEmaxOE5x5[i][j]);


      helIsoTrack[i][5]->Add(helIsoTrack[i][j]);
      helIsoEcal[i][5]->Add(helIsoEcal[i][j]);
      helIsoHcal[i][5]->Add(helIsoHcal[i][j]);

      helIsoTrackb[i][5]->Add(helIsoTrackb[i][j]);
      helIsoEcalb[i][5]->Add(helIsoEcalb[i][j]);
      helIsoHcalb[i][5]->Add(helIsoHcalb[i][j]);

      helIsoTracke[i][5]->Add(helIsoTracke[i][j]);
      helIsoEcale[i][5]->Add(helIsoEcale[i][j]);
      helIsoHcale[i][5]->Add(helIsoHcale[i][j]);
           
      hdilMass[i][6]->Add(hdilMass[i][j]);
      hdilMassTightWindow[i][6]->Add(hdilMassTightWindow[i][j]);
      hdilPt[i][6]->Add(hdilPt[i][j]);
      hmet[i][6]->Add(hmet[i][j]);
      hmetPhi[i][6]->Add(hmetPhi[i][j]);
      hpfmet[i][6]->Add(hpfmet[i][j]);
      hpfmetPhi[i][6]->Add(hpfmetPhi[i][j]); 
      htcmet[i][6]->Add(htcmet[i][j]);
      htcmetPhi[i][6]->Add(htcmetPhi[i][j]);
      hptJet1[i][6]->Add(hptJet1[i][j]);
      hptJet2[i][6]->Add(hptJet2[i][j]);
      hptJet3[i][6]->Add(hptJet3[i][j]);
      hptJet4[i][6]->Add(hptJet4[i][j]);
      hetaJet1[i][6]->Add(hetaJet1[i][j]);
      hetaJet2[i][6]->Add(hetaJet2[i][j]);
      hetaJet3[i][6]->Add(hetaJet3[i][j]);
      hetaJet4[i][6]->Add(hetaJet4[i][j]);
      heleSumPt[i][6]->Add(heleSumPt[i][j]);
      hmuSumPt[i][6]->Add(hmuSumPt[i][j]);
      
      hmuSumIso[i][6]->Add(hmuSumIso[i][j]);
      helSumIso[i][6]->Add(helSumIso[i][j]);
      hmuRelIso[i][6]->Add(hmuRelIso[i][j]);
      helRelIso[i][6]->Add(helRelIso[i][j]);
      hmuRelIsoTrack[i][6]->Add(hmuRelIsoTrack[i][j]);
      helRelIsoTrack[i][6]->Add(helRelIsoTrack[i][j]);
      hmuRelIsoCalo[i][6]->Add(hmuRelIsoCalo[i][j]);
      helRelIsoCalo[i][6]->Add(helRelIsoCalo[i][j]);
      
      hmt2[i][6]->Add(hmt2[i][j]);
      hmt2J[i][6]->Add(hmt2J[i][j]);

    }
 
  }

}

