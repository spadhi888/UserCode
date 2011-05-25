#include <iostream>
#include <fstream>
#include <vector>

#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"

//#include "CMS2.h"
#include "CMS2.cc"

//#include "../CORE/selections.cc"
//#include "../CORE/utilities.cc"
//#include "../Tools/tools.cc"
#include "../CORE/electronSelections.cc"
#include "../CORE/electronSelectionsParameters.cc"
#include "../CORE/MITConversionUtilities.cc"
#include "../CORE/trackSelections.cc"

using namespace tas;
//CMS2 cms2;

TH1F* book1DHist(const char* name, const char* title, unsigned int nbins, float low, float high, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,nbins,low,high);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH1F* book1DHist(const char* name, const char* title, unsigned int xnbins, const float *xbin, const char* xtitle, const char* ytitle, int color) {
  // return histogram instance with called Sumw2
  TH1F *hist = new TH1F(name,title,xnbins,xbin);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->Sumw2();
  hist->SetFillColor(color);
  hist->SetLineColor(color);
   
  return hist;   
}

TH2F* book2DHist(const char* name, const char* title, 
                 unsigned int xnbins, float xlow, float xhigh, 
                 unsigned int ynbins, float ylow, float yhigh, 
                 const char* xtitle, const char* ytitle, const char* ztitle, 
                 int color) {
  // return histogram instance with called Sumw2
  TH2F *hist = new TH2F(name,title,xnbins,xlow,xhigh,ynbins,ylow,yhigh);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2(); 
  hist->SetFillColor(color);
  hist->SetStats(kFALSE);
  //hist->SetLineColor(color);
   
  return hist;   
}

TH2F* book2DHistVar(const char* name, const char* title, 
                 unsigned int xnbins, const float *xbin, 
                 unsigned int ynbins, const float *ybin, 
                 const char* xtitle, const char* ytitle, const char* ztitle, 
                 int color) {
  // return histogram instance with called Sumw2
  TH2F *hist = new TH2F(name,title,xnbins,xbin,ynbins,ybin);
  hist->SetXTitle(xtitle);
  hist->SetYTitle(ytitle);
  hist->SetZTitle(ztitle);
  hist->Sumw2(); 
  hist->SetFillColor(color);
  hist->SetStats(kFALSE);
  //hist->SetLineColor(color);
   
  return hist;   
}

int ScanChain( TChain* chain) {

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain=chain->GetEntries();;
  unsigned int nEventsTotal = 0;

  // book histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  unsigned int nBinsPt 	= 55;
  //unsigned int nBinsPt= 5;
  float lowBinsPt 	= 0.;
  float highBinsPt 	= 110.;

  unsigned int nBinsEta = 52;
  //unsigned int nBinsEta = 5;
  float lowBinsEta      = -2.6;
  float highBinsEta     =  2.6;

  // unfixed bin size
  //float ptbins[10]= {10, 20, 30, 40, 50, 60, 70, 80, 90, 150};
  //float etabins[8] = {0, 0.5, 1.0, 1.28, 1.56, 1.84, 2.12, 2.5};
  //float etabins[] = {0, 0.5, 1.0, 1.28, 1.56, 1.84, 2.12, 2.4};
  //float etabins[] = {0, 1.28, 1.56, 1.84, 2.12, 2.5};
//   float ptbins[] = {10, 30, 50, 70, 100};
//  float etabins[] = {0, 0.5, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5};
//  float ptbins[] = {10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110};
//   float etabins[] = {0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};
//  float etabins[] = {-2.5, -2.2, -2.1, -2.0, -1.8, -1.479, -1.0, -0.5, 0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};
// - Last year
  float ptbins[] = {10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
//  float etabins[] = {0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};
  float etabins[] = {-2.5, -2.2, -2.1, -2.0, -1.8, -1.479, -1.0, -0.5, 0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};


  unsigned etanbins = sizeof(etabins)/sizeof(etabins[0])-1;
  unsigned ptnbins = sizeof(ptbins)/sizeof(ptbins[0])-1; 

  TH1F *els_pt_sim 			= book1DHist("els_pt_sim", 
						     "true Electron: p_{T}^{true}",
                                                     ptnbins,
                                                     ptbins,
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_sim 			= book1DHist("els_eta_sim", 
						     "true Electron: #eta^{true}",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
                                                     etanbins,
                                                     etabins,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim 			= book1DHist("els_pt_recosim", 
						     "true Electron: p_{T}^{true}",
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
                                                     ptnbins,
                                                     ptbins,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim 		= book1DHist("els_eta_recosim", 
						     "true Electron: #eta^{true}",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
                                                     etanbins,
                                                     etabins,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco 			= book1DHist("els_pt_reco", 
						     "reco Electron: p_{T}^{true}",
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
                                                     ptnbins,
                                                     ptbins,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco 			= book1DHist("els_eta_reco", 
						     "reco Electron: #eta^{true}",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
                                                     etanbins,
                                                     etabins,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_incorCharge 	= book1DHist("els_pt_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: p_{T}^{true}",
						     ptnbins,
						     ptbins,
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_incorCharge 	= book1DHist("els_eta_recosim_incorCharge", 
						     "true Electron with incorrect reconstructed Charge: #eta^{true}",
						     etanbins,
						     etabins,
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_recosim_corCharge 	= book1DHist("els_pt_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: p_{T}^{true}",
						     ptnbins,
						     ptbins,
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_recosim_corCharge 	= book1DHist("els_eta_recosim_corCharge", 
						     "true Electron with correct reconstructed Charge: #eta^{true}",
						     etanbins,
						     etabins,
						     //						     nBinsEta,
						     // 						     lowBinsEta,
						     //						     highBinsEta,
						     "#eta^{true}",
						     "Electrons",2);
  TH1F *els_pt_reco_corCharge 		= book1DHist("els_pt_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: p_{T}^{true}",
						     //						     nBinsPt,
						     //						     lowBinsPt,
						     //						     highBinsPt,
                                                     ptnbins,
                                                     ptbins,
						     "p_{T}^{true} [GeV]",
						     "Electrons",2);
  TH1F *els_eta_reco_corCharge 		= book1DHist("els_eta_reco_corCharge", 
						     "reco Electron with correct reconstructed Charge: #eta^{true}",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
                                                     etanbins,
                                                     etabins,
						     "#eta^{true}",
						     "Electrons",2);

  TH1F *els_trkId         		= book1DHist("els_trkId", 
						     "Track ID of associated to reco Electron",
						     30,
						     -1050.,
						     450.,
						     "Track ID",
						     "Electrons",2);

  TH1F *els_eta_before3agree 		= book1DHist("els_eta_before3agree", 
						     "",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
						     etanbins,
						     etabins,
						     "#eta^{true}",
						     "Electrons",2);

  TH1F *els_eta_after3agree 		= book1DHist("els_eta_after3agree", 
						     "",
						     //						     nBinsEta,
						     //						     lowBinsEta,
						     //						     highBinsEta,
						     etanbins,
						     etabins,
						     "#eta^{true}",
						     "Electrons",2);

  //
  // 2D histograms 
  //
  TH2F *els_2d_eta_Pt_corCharge		= book2DHistVar("els_2d_eta_Pt_corCharge",
						     "2D histogram of eta vs Pt corCharge",
						     etanbins,
						     etabins,
						     ptnbins,
						     ptbins,
						     "#eta",
						     "p_{T} [GeV]",
						     " ",2);

  TH2F *els_2d_eta_Pt_incorCharge	= book2DHistVar("els_2d_eta_Pt_incorCharge",
						     "2D histogram of eta vs Pt incorCharge",
						     etanbins,
						     etabins,
						     ptnbins,
						     ptbins,
						     "#eta",
						     "p_{T} [GeV]",
						     " ",2);

  TH2F *els_2d_eta_Pt_ratio 		= book2DHistVar("els_2d_eta_Pt_ratio",
						     "2D histogram of eta vs Pt ratio",
						     etanbins,
						     etabins,
						     ptnbins,
						     ptbins,
						     "#eta",
						     "p_{T} [GeV]",
						     //"MisID rate",2);
						     " ",2);

  int num_beforecut=0;
  int num_ELEID_CAND02=0;
  int num_ELEID_EXTRA=0;
  int num_ELENOTCONV_DISTDCOT002=0;
  int num_ELENOTCONV_HITPATTERN=0;
  int num_ELECHARGE_NOFLIP=0;

  // individual cuts in ELEID_CAND02
  int num_ELEID_CAND02_ELESEED_ECAL=0;
  int num_ELEID_CAND02_ELEETA_250=0;
  int num_ELEID_CAND02_ELENOMUON_010=0;
  int num_ELEID_CAND02_ELEID_CAND02=0;
  int num_ELEID_CAND02_ELEISO_REL010=0;
  int num_ELEID_CAND02_ELENOTCONV_DISTDCOT002=0;
  int num_ELEID_CAND02_ELEIP_200=0;


  int num_cand02=0;
  int num_extra=0;
  int num_hitpattern=0;
  int num_partnertrack=0;
  int num_chargeflip=0;
  int num_cand02_ecal=0;
  int num_cand02_eta250=0;
  int num_cand02_noMuon=0;
  int num_cand02_cand02=0;
  int num_cand02_impact=0;
  int num_cand02_relsusy010=0;


  // file loop
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( currentFile = (TFile*)fileIter.Next() ) {
    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    //Event Loop
    unsigned int nEvents = tree->GetEntries();
    // nEvents = 100;
    for( unsigned int event = 0; event < nEvents; ++event) {
      cms2.GetEntry(event);
      ++nEventsTotal;
     
      if ( nEventsTotal%100000 == 0 ) {
	std::cout << "Event: " << nEventsTotal << endl;
      }


      // loop over true electrons
      for ( unsigned int els = 0;
	    els < genps_p4().size();
	    ++els ) {

	// pt cut
	if( (genps_p4()[els].pt()) < 10 ) continue;
	if( fabs(genps_p4()[els].eta()) > 2.4 ) continue;


	// check that electron is final state electron
	if ( genps_status()[els] != 1 ) continue;
	
	// check for true electron
	if ( TMath::Abs(genps_id()[els]) != 11 ) continue;

	// fill true histrograms
	els_pt_sim->Fill(genps_p4()[els].pt());
	els_eta_sim->Fill(genps_p4()[els].eta());

      }

// down to here, no problem.

      // loop over reco electrons
      for (unsigned int els = 0; 
           els < els_p4().size(); 
	   ++els) {

	num_beforecut++;

	//
	// cuts
	//

	// no missing hit
//	if(els_exp_innerlayers()[els] != 0) continue;

	// pt
	if( (els_p4().at(els).Pt()) < 10 ) continue;
	if( fabs(els_p4().at(els).eta()) > 2.4 ) continue;

	if (!pass_electronSelection(els, electronSelection_ssV3, false, false)) continue;
	els_eta_before3agree->Fill((els_p4().at(els).eta()));

	// 3 charges agree
//	if (!pass_electronSelection(els, electronSelection_ss_Flip)) continue;
	els_eta_after3agree->Fill((els_p4().at(els).eta())); 

	// check how many electrons don't have an associated track
	els_trkId->Fill(els_trkidx().at(els));
	// tmp charge variable
	//double charge = els_charge().at(els); //
	double charge = els_trk_charge().at(els);
	
	// fill reco
	els_pt_reco->Fill(els_p4().at(els).Pt());
	els_eta_reco->Fill((els_p4().at(els).eta()));

	// fill reco_corCharge
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {
	  els_pt_reco_corCharge->Fill(els_p4().at(els).Pt());
	  els_eta_reco_corCharge->Fill((els_p4().at(els).eta()));
	}

	// exclude reco which has no true electron match
	if ( abs(els_mc_id()[els]) != 11 ) continue;

	// fill recosim
	els_pt_recosim->Fill(els_mc_p4().at(els).Pt());
	els_eta_recosim->Fill((els_mc_p4().at(els).eta()));

	// correct charge identified
	if ( (charge == -1 && els_mc_id().at(els) == 11) || (charge == 1 && els_mc_id().at(els) == -11) ) {

	  els_pt_recosim_corCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_corCharge->Fill((els_mc_p4().at(els).eta()));
	  els_2d_eta_Pt_corCharge->Fill((els_mc_p4().at(els).eta()), els_mc_p4().at(els).Pt()); //2d

	  // incorrect charge identified
	} else {

	  els_pt_recosim_incorCharge->Fill(els_mc_p4().at(els).Pt());
	  els_eta_recosim_incorCharge->Fill((els_mc_p4().at(els).eta()));
	  els_2d_eta_Pt_incorCharge->Fill((els_mc_p4().at(els).eta()), els_mc_p4().at(els).Pt()); //2d
	}

      }
      
    }
  }




  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  return 0;
}
