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
#include "../CORE/CMS2.cc"
#include "goodrun.cc"
// #include "../CORE/fakerates.cc"
#include "../CORE/mcSelections.cc"
#include "../CORE/MT2/MT2.cc"
#include "../CORE/trackSelections.cc"
#include "../CORE/eventSelections.cc"
#include "../CORE/SimpleFakeRate.cc"

#define JETPTCUT 20.0

using namespace tas;
unsigned int numtightLeps;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

bool usePtGt2020;
bool usePtGt2010;
bool usePtGt1010;
bool excludePtGt2020;
bool used0corrPV;
bool applylepIDCuts;
bool applyFOv1Cuts;
bool applyFOv2Cuts;
bool applyFOv3Cuts;
bool applyLooseIDCuts;
bool applylepIsoCuts;
bool applylepLooseIsoCuts;
bool applyTriggers;
bool vetoZmass;
bool requireZmass;
bool hypDisamb;
bool useCorMET;
bool useOSleptons;
bool useSSleptons;
bool usetcMET;
bool usepfMET;
bool vetoMET;
bool vetoProjectedMET;
bool usejptJets;
bool usecaloJets;
bool usepfJets;
bool vetoJets;
bool requireEcalEls;
bool chargeFlip;
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
void EndJob();


void ScanChain( TChain* chain, vector<TString> v_Cuts, string prefix="", 
		bool doFRestimation = false, float jetTriggerPt = -9999., float lumi = 0.01, float NLOCS=-9999.) { //lumi in fb-1


  //deal with the cuts
  useOSleptons          = find(v_Cuts.begin(), v_Cuts.end(), "useOSleptons"          ) != v_Cuts.end();
  useSSleptons          = find(v_Cuts.begin(), v_Cuts.end(), "useSSleptons"          ) != v_Cuts.end();
  usePtGt2020          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2020"          ) != v_Cuts.end();
  usePtGt2010          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt2010"          ) != v_Cuts.end();
  usePtGt1010          = find(v_Cuts.begin(), v_Cuts.end(), "usePtGt1010"          ) != v_Cuts.end();      
  used0corrPV          = find(v_Cuts.begin(), v_Cuts.end(), "used0corrPV"          ) != v_Cuts.end();
  excludePtGt2020      = find(v_Cuts.begin(), v_Cuts.end(), "excludePtGt2020"      ) != v_Cuts.end();
  applylepIDCuts       = find(v_Cuts.begin(), v_Cuts.end(), "applylepIDCuts"       ) != v_Cuts.end(); 
  applyFOv1Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv1Cuts"        ) != v_Cuts.end(); 
  applyFOv2Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv2Cuts"        ) != v_Cuts.end(); 
  applyFOv3Cuts        = find(v_Cuts.begin(), v_Cuts.end(), "applyFOv3Cuts"        ) != v_Cuts.end(); 
  applyLooseIDCuts     = find(v_Cuts.begin(), v_Cuts.end(), "applyLooseIDCuts"     ) != v_Cuts.end(); 	  
  applylepIsoCuts      = find(v_Cuts.begin(), v_Cuts.end(), "applylepIsoCuts"      ) != v_Cuts.end(); 
  applylepLooseIsoCuts = find(v_Cuts.begin(), v_Cuts.end(), "applylepLooseIsoCuts" ) != v_Cuts.end();
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
  chargeFlip            = find(v_Cuts.begin(), v_Cuts.end(), "chargeFlip"            ) != v_Cuts.end();
  vetoJets             = find(v_Cuts.begin(), v_Cuts.end(), "vetoJets"             ) != v_Cuts.end();
  requireEcalEls       = find(v_Cuts.begin(), v_Cuts.end(), "requireEcalEls"       ) != v_Cuts.end();

  cout << "REMEMBER THAT JET COUNTING NOW USES 20 GeV JETS" << endl;


  if(doFRestimation && jetTriggerPt < 0) {
    cout << "Need to specify the jetTrigger to do the FRs. Exiting" << endl;
    return;
  }

  //preliminary checks
  if(usePtGt2020 && usePtGt1010) {
    cout << "You want to use pt > 20,20, and >10,10. While this is ok, its a bit redundant. Please pick the cut you really want" << endl;
    return;
  }

  if(!usePtGt2020 && excludePtGt2020) {
    cout << "Cannot exlcude higer pt region and not use the lower Dilepton pt region. Will end up with no events passing in the greater than 2 jet bin" << endl;
    return;
  }

  if(applylepIDCuts  + applyLooseIDCuts + applyFOv1Cuts + applyFOv2Cuts +applyFOv3Cuts > 1 ) {
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
      
      if(prefix == "data") {
        DorkyEventIdentifier id = { evt_run(),evt_event(),evt_lumiBlock() };
        if (is_duplicate(id)) continue;
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
      if(prefix == "data" && !goodrun(evt_run(), evt_lumiBlock()))
	continue;

      // Cleaning cuts.
      if(!cleaning_standard(prefix == "data")) continue;
      //      if (!cleaning_goodVertex()) continue;
      //      if (!cleaning_goodTracks()) continue;

      //get the channels correct
      int nels, nmus, ntaus;
      int nleps = 0;
      if(prefix != "data")
	nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);
      if (prefix == "ttdil"   &&  nleps  != 2) continue;
      if (prefix == "ttotr"    &&  nleps == 2) continue;
      if (prefix == "DYee"     &&  nleps != 2) continue;
      if (prefix == "DYmm"     &&  nleps != 2) continue;
      if (prefix == "DYtautau" &&  nleps != 2) continue;
      
         
      //10pb-1 
      float weight = 1.;

      if(prefix != "data") {
	if(NLOCS > 0)
	  weight = evt_scale1fb()*lumi*NLOCS/(evt_xsec_incl()*evt_kfactor());
	else 
	  weight = evt_scale1fb()*lumi;	  
      }

      
      vector<unsigned int> v_goodHyps;
      v_goodHyps.clear();
      vector<float> v_weights;
      v_weights.clear();
      

      // Electrons
      std::vector<LorentzVector> eleCollection;
      eleCollection.clear();

      for(int iEl = 0 ; iEl < els_p4().size(); iEl++) {
	
	Double_t pt = els_p4()[iEl].Pt();
	Double_t eta = els_p4()[iEl].Eta();
	//	Double_t phi = els_p4()[iEl].Phi();

	if (!(els_type()[iEl] & (1<<2))) continue; 
	if (pt < 10 ) continue;  
	if (fabs(eta) > 2.5)  continue;   

 	// ID cuts.
	if (!isGoodLeptonwIso(11, iEl)) continue;  
	//	cout << els_charge()[iEl] << endl;
	eleCollection.push_back(els_p4()[iEl]);
      }

      if (eleCollection.size() > 1) {
	sort(eleCollection.begin(), eleCollection.end(),  sortByPt);
      }


      // Muons
      std::vector<LorentzVector> musCollection;
      musCollection.clear();

      for(int iMu = 0 ; iMu < mus_p4().size(); iMu++) {
	
	Double_t pt = mus_p4()[iMu].Pt();
	Double_t eta = mus_p4()[iMu].Eta();
	//	Double_t phi = mus_p4()[iMu].Phi();

        if (pt < 10 ) continue;   
        if (fabs(eta) > 2.5)  continue;    
	if (((mus_type()[iMu]) & (1<<1)) == 0)    continue; // global muon
	if (((mus_type()[iMu]) & (1<<2)) == 0)    continue; // tracker muon
	
	// ID cuts.
	if (!isGoodLeptonwIso(13, iMu)) continue;  
	// cout << mus_charge()[iMu] << endl;
        musCollection.push_back(mus_p4()[iMu]); 
      }
      
      if (musCollection.size() > 1) {
	sort(musCollection.begin(), musCollection.end(),  sortByPt);
      }

      // Jets

	vector<unsigned int> v_goodJets;
	vector<LorentzVector> v_jetP4s;

	if(usecaloJets) {
	  for(unsigned int i = 0; i < jets_p4().size(); i++) 
	    v_jetP4s.push_back(jets_p4()[i]*jets_cor()[i]); //jets are uncorrected in our ntuples
	}

	if(usejptJets)
	  v_jetP4s = jpts_p4();
	if(usepfJets) 
	  v_jetP4s = pfjets_p4();
	for (unsigned int j = 0; j < v_jetP4s.size(); ++j) {
	  if (v_jetP4s.at(j).Pt() < JETPTCUT) continue; 
	  if (fabs(v_jetP4s.at(j).Eta()) > 2.5) continue;
	  
	  bool overlap = false;
	  for(unsigned int k = 0; k < eleCollection().size(); k++) {
	    if (dRbetweenVectors(eleCollection()[k],v_jetP4s()[j]) < 0.4) overlap = true; 
	  }
	  for(unsigned int k = 0; k < eleCollection().size(); k++) {
	    if (dRbetweenVectors(musCollection()[k],v_jetP4s()[j]) < 0.4) overlap = true; 
	  }
	  if (overlap) continue;
	  v_goodJets.push_back(j);
	}


	// MET cut
	string metAlgo;
	if(useCorMET   ) metAlgo  = "CorMET";
	if(usetcMET    ) metAlgo  = "tcMET";
	if(usepfMET    ) metAlgo  = "pfMET"; 	
	pair<float, float> p_met; //met and met phi
	if(usetcMET || usepfMET) {
	  p_met = make_pair(evt_pfmet(), evt_pfmetPhi());
	  
	  if(p_met.first < 0) {
	    cout << "Something is wrong with the MET. Exiting" << endl;
	    return;
	  }
	  
	  cout << p_met.first << endl; 
	}

	// Hypothesis based analysis
	
      for(unsigned int hypIdx = 0; hypIdx < hyp_p4().size(); hypIdx++) {//hyploop
	// UNUSED 
	
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
  
  for(unsigned int i = 0; i < v_type.size(); i++) {
    
    //which channel?
    unsigned int ch  = v_type.at(i);
      
    
    hnJet[ch]                         ->Fill(jetBin,        weight);
    
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
    double tcmet = evt_tcmet();
    double tcmetPhi = evt_tcmetPhi();
    if(prefix == "data") {
      metStruct cortcmet = correctedTCMET();
      tcmet = cortcmet.met;
      tcmetPhi = cortcmet.metphi;
    }
    correctTcMETForHypMus(hypIdx, tcmet, tcmetPhi);
    htcmet[ch][jetBin]                ->Fill(tcmet,           weight); 
    htcmetPhi[ch][jetBin]             ->Fill(tcmetPhi,        weight); 
    
    //projected MET
    hprojmet[ch][jetBin]              ->Fill(projectedMET(evt_metMuonJESCorr(),
							  evt_metMuonJESCorrPhi(),
							  hypIdx),    weight);  
    hprojpfmet[ch][jetBin]            ->Fill(projectedMET(evt_pfmet(), evt_pfmetPhi(),
							  hypIdx), weight); 
    hprojtcmet[ch][jetBin]            ->Fill(projectedMET(tcmet, tcmetPhi,
							  hypIdx), weight);


    hmetVsDilepPt[ch][jetBin]         ->Fill(hypp4.Pt(), evt_metMuonCorr(), weight); 
    hmetOverPtVsDphi[ch][jetBin]      ->Fill(dphi, evt_metMuonCorr()/hypp4.Pt(), weight); 

    hpfmetVsDilepPt[ch][jetBin]       ->Fill(hypp4.Pt(), evt_pfmet(),       weight); 
    hpfmetOverPtVsDphi[ch][jetBin]    ->Fill(dphi, evt_pfmet()/hypp4.Pt(),  weight);  

    htcmetVsDilepPt[ch][jetBin]       ->Fill(hypp4.Pt(), tcmet,       weight);
    htcmetOverPtVsDphi[ch][jetBin]    ->Fill(dphi, tcmet/hypp4.Pt(),  weight);
    */
    hdphillvsmll[ch][jetBin]          ->Fill(hypp4.mass(), dphi, weight);

    
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
    
    //numTightLep[ch][jetBin]; 
    //heleSumPt[ch][jetBin]; 
    //hmuSumPt[ch][jetBin]; 
   
    
        
    //hnJetLepVeto[ch]; 

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
//    hmt2[ch][jetBin]                   ->Fill(MT2(p_met.first, p_met.second, lt_p4, ll_p4, 0.0, false), weight);	
//    if(v_jets.size() > 1) 
//      hmt2J[ch][jetBin]                ->Fill(MT2J(p_met.first, p_met.second, lt_p4, ll_p4, v_jets_p4, 0.0, BISECT, false),weight);	
  }
    
}

// ****************************************************************
// Does the event pass the monster event cuts??
// ****************************************************************

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

