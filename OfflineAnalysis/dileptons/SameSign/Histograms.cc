#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TMath.h"
#include "TDirectory.h"

#define MAXPT 1000
#define MINPT 0
#define NBINSPT 50
#define NBINSETA 60

#define RELISOMIN 0
#define RELISOMAX 10
#define NBINSRELISO 100

//histos from:http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/JRibnik/CMS2/NtupleMacros/TTDil/slava/ttDilCounts_looper.C?revision=1.20.4.1&view=markup&pathrev=PAS-TOP-09-002_vFinal_072009


//in the 2d arrays, the first index is the different hypothesis types (ee, mm, em, all)
//the second index is the jet bin (0 jet, 1 jet, greater than 2 jets


TH1F* hnJet[4];                   // Njet distributions
TH1F* hnJetWW[4];
TH1F* hnJetWO[4];
TH1F* hnJetWOSemilep[4];
TH1F* hnJetWOOther[4];
TH1F* hnJetOO[4];
TH1F* hnJetinZwindow[4];          //usefull for DY estimate
TH1F* hnJetoutZwindow[4];         //usefull for DY estimate
TH1F* hnTrueFVFit[4];                // True FVFit prob
TH1F* hnFVFit[4];                // Other FVFit prob
TH1F* hnd0Corr[4];                // d0Corr of the other lepton
TH1F* hnd0TrueCorr[4];                // d0Corr of the other lepton
TH1F* hnEleIsolation[4];          //Isolation of bad electron
TH1F* hnEleTrueIsolation[4];          //Isolation of bad electron
TH1F* hnMuIsolation[4];          //Isolation of bad muon
TH1F* hnMuTrueIsolation[4];          //Isolation of bad muon
TH1F* hnMeff[4];                   // Njet distributions
TH1F* hnSumptj[4];
TH1F* hnMET[4];
TH1F* helePt[4][7];               // electron Pt
TH1F* hmuPt[4][7];                // muon Pt
TH1F* hmuPtFromSilicon[4][7];     // muon Pt (from tracker)
TH1F* hminLepPt[4][7];            // minimum lepton Pt
TH1F* hmaxLepPt[4][7];            // maximum lepton Pt
TH1F* helePhi[4][7];              // electron phi
TH1F* hmuPhi[4][7];               // muon phi
TH1F* hdphiLep[4][7];             // delta phi between leptons
TH1F* heleEta[4][7];              // electron eta
TH1F* hmuEta[4][7];               // muon eta
TH1F* heled0BS[4][7];             // electron d0 corrected for the beamspot
TH1F* heled0PV[4][7];             // electron d0 corrected for the primary vertex
TH1F* hmud0BS[4][7];              // muon d0 corrected for the beamspot
TH1F* hmud0PV[4][7];              // muon d0 corrected for the primary vertex
TH1F *heleEmaxOE5x5[4][7];        // Emax / E5x5 for electrons


TH1F* hdilMass[4][7];             // dilepton mass
TH1F* hdilMassTightWindow[4][7];  // dilepton mass, but zooming around Z
TH1F* hdilPt[4][7];               // dilepton Pt
TH1F* hmet[4][7];                 // MET
TH1F* hmetPhi[4][7];              // MET phi
TH1F* hpfmet[4][7];               // PF MET
TH1F* hpfmetPhi[4][7];            // PF MET phi
TH1F* htcmet[4][7];               // tc MET
TH1F* htcmetPhi[4][7];            // tc MET phi
TH1F* hprojmet[4][7];             // projected caloMET
TH1F* hprojpfmet[4][7];           // PF MET
TH1F* hprojtcmet[4][7];           // tc MET



TH2F* hmetVsDilepPt[4][7];        // MET vs dilepton Pt
TH2F* hmetOverPtVsDphi[4][7];     // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

TH2F* hpfmetVsDilepPt[4][7];      // PAT MET vs dilepton Pt
TH2F* hpfmetOverPtVsDphi[4][7];   // PAT MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

TH2F* htcmetVsDilepPt[4][7];      // tc MET vs dilepton Pt
TH2F* htcmetOverPtVsDphi[4][7];   // tc MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

TH2F* hdphillvsmll[4][7];         // delta phi between leptons vs dilepton mass
TH1F* hptJet1[4][7];              // Pt of 1st jet
TH1F* hptJet2[4][7];              // Pt of 2nd jet
TH1F* hptJet3[4][7];              // Pt of 3rd jet
TH1F* hptJet4[4][7];              // Pt of 4th jet
TH1F* hetaJet1[4][7];             // eta of 1st jet
TH1F* hetaJet2[4][7];             // eta of 2nd jet
TH1F* hetaJet3[4][7];             // eta of 3rd jet
TH1F* hetaJet4[4][7];             // eta of 4th jet
TH1F* numTightLep[4][7];          // number of tight leptons per event.
TH1F* heleSumPt[4][7];            // sumPt for electron isolation
TH1F* hmuSumPt[4][7];             // sumPt for muon isolation
TH1F* hmuSumIso[4][7];            // sum of trk pt, em et, had et in cone of 0.3
TH1F* helSumIso[4][7];            // sum of trk pt, em et, had et in cone of 0.3
TH1F* helRelIso[4][7];            //  Iso variable defined as pt/(pt+sum) for electron
TH1F* hmuRelIso[4][7];            //  Iso variable defined as pt/(pt+sum) for muons
TH1F* helRelIsoTrack[4][7];       //  Iso variable defined as pt/(pt+sum) for electron
TH1F* hmuRelIsoTrack[4][7];       //  Iso variable defined as pt/(pt+sum) for muons
TH1F* helRelIsoCalo[4][7];        //  Iso variable defined as pt/(pt+sum) for electron
TH1F* hmuRelIsoCalo[4][7];        //  Iso variable defined as pt/(pt+sum) for muons

TH1F* helIsoTrack[4][7];           // Jurassic Track Isolation
TH1F* helIsoEcal[4][7];            // Jurassic Ecal Isolation
TH1F* helIsoHcal[4][7];            // Hcal Isolation

TH1F* helIsoTrackb[4][7];         // Jurassic Track Isolation, barrel
TH1F* helIsoEcalb[4][7];          // Jurassic Ecal Isolation, barrel
TH1F* helIsoHcalb[4][7];          // Hcal Isolation, barrel

TH1F* helIsoTracke[4][7];         // Jurassic Track Isolation, endcap
TH1F* helIsoEcale[4][7];          // Jurassic Ecal Isolation, endcap
TH1F* helIsoHcale[4][7];          // Hcal Isolation, endcap

TH1F* hnJetLepVeto[4];            //njet distribution after requiring numTightLep < 3.


//MT2 stuff
TH1F* hmt2[4][7];
TH1F* hmt2J[4][7];//MT2 with jets only makes sense if you have >1 jet


void bookHistos(const char *prefix) {
  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  
  char *jetbins[5] = {"0", "1", "2", "3", "#geq 4"};
  char *jets[7]    = {"0", "1", "2", "3", "4", "geq2", "all"};
  char *suffixall[4];
  suffixall[0] = "ee";
  suffixall[1] = "mm";
  suffixall[2] = "em";
  suffixall[3] = "all";
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  for (int i=0; i<4; i++) {
    for (int j=0; j<7; j++) {
      
      if (j == 0){
	hnJet[i] = new TH1F((string(prefix)+"_hnJet_" + string(suffixall[i])).c_str(),(string(prefix)+"_nJet_" + string(suffixall[i])).c_str(),
			    5,0.,5.);	
        hnJetWW[i] = new TH1F((string(prefix)+"_hnJetWW_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetWW_" + string(suffixall[i])).c_str(),
                            5,0.,5.); 
        hnMeff[i] = new TH1F((string(prefix)+"_hnMeff_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnMeff_" + string(suffixall[i])).c_str(),
                            100,0.,2000.);
        hnSumptj[i] = new TH1F((string(prefix)+"_hnSumptj_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnSumptj_" + string(suffixall[i])).c_str(),
                            100,0.,2000.);
        hnMET[i] = new TH1F((string(prefix)+"_hnMET_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnMET_" + string(suffixall[i])).c_str(),
                            NBINSPT, MINPT, MAXPT);
        hnJetWO[i] = new TH1F((string(prefix)+"_hnJetWO_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetWO_" + string(suffixall[i])).c_str(),
                            5,0.,5.);
        hnJetWOSemilep[i] = new TH1F((string(prefix)+"_hnJetWOSemilep_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetWOSemilep_" + string(suffixall[i])).c_str(),
                            5,0.,5.);
        hnJetWOOther[i] = new TH1F((string(prefix)+"_hnJetWOOther_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetWOOther_" + string(suffixall[i])).c_str(),
                            5,0.,5.);
        hnJetOO[i] = new TH1F((string(prefix)+"_hnJetOO_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetOO_" + string(suffixall[i])).c_str(),
                            5,0.,5.);
        hnd0Corr[i] = new TH1F((string(prefix)+"_hnd0Corr_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnd0Corr_" + string(suffixall[i])).c_str(), 100,-0.05,0.05);
        hnd0TrueCorr[i] = new TH1F((string(prefix)+"_hnd0TrueCorr_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnd0TrueCorr_" + string(suffixall[i])).c_str(), 100,-0.05,0.05);
        hnTrueFVFit[i] = new TH1F((string(prefix)+"_hnTrueFVFit_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnTrueFVFit_" + string(suffixall[i])).c_str(), 100,0.,1);
        hnFVFit[i] = new TH1F((string(prefix)+"_hnFVFit_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnFVFit_" + string(suffixall[i])).c_str(), 100,0.,1);
        hnEleIsolation[i] = new TH1F((string(prefix)+"_hnEleIsolation_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnEleIsolation_" + string(suffixall[i])).c_str(), 100,0.,1.);
        hnEleTrueIsolation[i] = new TH1F((string(prefix)+"_hnEleTrueIsolation_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnEleTrueIsolation_" + string(suffixall[i])).c_str(), 100,0.,1.);
        hnMuIsolation[i] = new TH1F((string(prefix)+"_hnMuIsolation_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnMuIsolation_" + string(suffixall[i])).c_str(), 100,0.,1.);
        hnMuTrueIsolation[i] = new TH1F((string(prefix)+"_hnMuTrueIsolation_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnMuTrueIsolation_" + string(suffixall[i])).c_str(), 100,0.,1.);

	hnJetinZwindow[i] = new TH1F((string(prefix)+"_hnJetinZwindow_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetinZwindow_" + string(suffixall[i])).c_str(),
				     5,0.,5.);	
	hnJetoutZwindow[i] = new TH1F((string(prefix)+"_hnJetoutZwindow_" + string(suffixall[i])).c_str(),(string(prefix)+"_hnJetoutZwindow_" + string(suffixall[i])).c_str(),
				      5,0.,5.);	
	
	hnJet[i]->SetDirectory(rootdir);
	hnJet[i]->GetXaxis()->SetTitle("nJets");

	hnJetinZwindow[i]->SetDirectory(rootdir);
	hnJetinZwindow[i]->GetXaxis()->SetTitle("nJets");

	hnJetoutZwindow[i]->SetDirectory(rootdir);
	hnJetoutZwindow[i]->GetXaxis()->SetTitle("nJets");

	for(int k = 0; k<5; k++) {
	  hnJet[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJet[i]->GetXaxis()->SetLabelSize(0.07);
	  
	  hnJetinZwindow[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJetinZwindow[i]->GetXaxis()->SetLabelSize(0.07);
	  
	  hnJetoutZwindow[i]->GetXaxis()->SetBinLabel(k+1, jetbins[k]);
	  hnJetoutZwindow[i]->GetXaxis()->SetLabelSize(0.07);
	  
	}
      }
    
      char *suffix[4];
      suffix[0] = Form("%sj_ee",  jets[j]);
      suffix[1] = Form("%sj_mm",  jets[j]);
      suffix[2] = Form("%sj_em",  jets[j]);
      suffix[3] = Form("%sj_all", jets[j]);
      
      helePt[i][j] = new TH1F((string(prefix)+"_helePt_" + string(suffix[i])).c_str(),(string(prefix)+"_elePt_" + string(suffix[i])).c_str(),
			      NBINSPT, MINPT, MAXPT);
      helePt[i][j]->SetDirectory(rootdir);
      helePt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmuPt[i][j]  = new TH1F((string(prefix)+"_hmuPt_" + string(suffix[i])).c_str(),(string(prefix)+"_muPt_" + string(suffix[i])).c_str(),
			      NBINSPT, MINPT, MAXPT);
      hmuPt[i][j]->SetDirectory(rootdir);
      hmuPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmuPtFromSilicon[i][j]  = new TH1F((string(prefix)+"_hmuPtFromSilicon_" + string(suffix[i])).c_str(),
					 (string(prefix)+"_muPtFromSilicon_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      hmuPtFromSilicon[i][j]->SetDirectory(rootdir);
      hmuPtFromSilicon[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hminLepPt[i][j]  = new TH1F((string(prefix)+"_hminLepPt_" + string(suffix[i])).c_str(),
				  (string(prefix)+"_minLepPt_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      hminLepPt[i][j]->SetDirectory(rootdir);
      hminLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      hmaxLepPt[i][j]  = new TH1F((string(prefix)+"_hmaxLepPt_" + string(suffix[i])).c_str(),
				  (string(prefix)+"_maxLepPt_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      hmaxLepPt[i][j]->SetDirectory(rootdir);
      hmaxLepPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
      
    
      helePhi[i][j] = new TH1F((string(prefix)+"_helePhi_" + string(suffix[i])).c_str(),(string(prefix)+"_elePhi_" + string(suffix[i])).c_str(),
			       50,-1*TMath::Pi(), TMath::Pi());
      helePhi[i][j]->SetDirectory(rootdir);
      helePhi[i][j]->GetXaxis()->SetTitle("#phi");
      

      hmuPhi[i][j]  = new TH1F((string(prefix)+"_hmuPhi_" + string(suffix[i])).c_str(),(string(prefix)+"_muPhi_" + string(suffix[i])).c_str(),
			       50,-1*TMath::Pi(), TMath::Pi());
      hmuPhi[i][j]->SetDirectory(rootdir);
      hmuPhi[i][j]->GetXaxis()->SetTitle("#phi");
      
    
      hdphiLep[i][j]  = new TH1F((string(prefix)+"_hdphiLep_" + string(suffix[i])).c_str(),(string(prefix)+"_dphiLep_" + string(suffix[i])).c_str(),
				 50,0., TMath::Pi());
      hdphiLep[i][j]->SetDirectory(rootdir);
      hdphiLep[i][j]->GetXaxis()->SetTitle("#delta#phi_{ll}");
      
      
      heleEta[i][j] = new TH1F((string(prefix)+"_heleEta_" + string(suffix[i])).c_str(),(string(prefix)+"_eleEta_" + string(suffix[i])).c_str(),
			       NBINSETA, -3., 3.);
      heleEta[i][j]->SetDirectory(rootdir);
      heleEta[i][j]->GetXaxis()->SetTitle("#eta");
      
	
      hmuEta[i][j]  = new TH1F((string(prefix)+"_hmuEta_" + string(suffix[i])).c_str(),(string(prefix)+"_muEta_" + string(suffix[i])).c_str(),
			       NBINSETA, -3., 3.);
      hmuEta[i][j]->SetDirectory(rootdir);
      hmuEta[i][j]->GetXaxis()->SetTitle("#eta");



      heled0BS[i][j]  = new TH1F((string(prefix)+"_heled0BS_" + string(suffix[i])).c_str(),(string(prefix)+"_heled0BS_" + string(suffix[i])).c_str(),
				 50,-0.1,0.1);
      heled0BS[i][j]->SetDirectory(rootdir);
      heled0BS[i][j]->GetXaxis()->SetTitle("d0 BScorr");
      heled0PV[i][j]  = new TH1F((string(prefix)+"_heled0PV_" + string(suffix[i])).c_str(),(string(prefix)+"_heled0PV_" + string(suffix[i])).c_str(),
				 50,-0.1,0.1);
      heled0PV[i][j]->SetDirectory(rootdir);
      heled0PV[i][j]->GetXaxis()->SetTitle("d0 PVcorr");


      hmud0BS[i][j]  = new TH1F((string(prefix)+"_hmud0BS_" + string(suffix[i])).c_str(),(string(prefix)+"_hmud0BS_" + string(suffix[i])).c_str(),
				50,-0.1,0.1);
      hmud0BS[i][j]->SetDirectory(rootdir);
      hmud0BS[i][j]->GetXaxis()->SetTitle("d0 BScorr");
      hmud0PV[i][j]  = new TH1F((string(prefix)+"_hmud0PV_" + string(suffix[i])).c_str(),(string(prefix)+"_hmud0PV_" + string(suffix[i])).c_str(),
				 50,-0.1,0.1);
      hmud0PV[i][j]->SetDirectory(rootdir);
      hmud0PV[i][j]->GetXaxis()->SetTitle("d0 PVcorr");

      heleEmaxOE5x5[i][j] = new TH1F((string(prefix)+"_heleEmaxOE5x5_" + string(suffix[i])).c_str(), (string(prefix)+"_heleEmaxOE5x5_" + string(suffix[i])).c_str(),
				     50,0,1);

	
      hdilMass[i][j] = new TH1F((string(prefix)+"_hdilMass_" + string(suffix[i])).c_str(),(string(prefix)+"_dilMass_" + string(suffix[i])).c_str(),
				300, 0., 300.);
      hdilMass[i][j]->SetDirectory(rootdir);
      hdilMass[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      

      hdilMassTightWindow[i][j] = new TH1F((string(prefix)+"_hdilMassTightWindow_" + string(suffix[i])).c_str(),
					   (string(prefix)+"_dilMassTightWindow_" + string(suffix[i])).c_str(),
					   120, 60., 120.);
      hdilMassTightWindow[i][j]->SetDirectory(rootdir);
      hdilMassTightWindow[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      
    
      hdilPt[i][j] = new TH1F((string(prefix)+"_hdilPt_" + string(suffix[i])).c_str(),(string(prefix)+"_dilPt_" + string(suffix[i])).c_str(),
			      100, 0., 300.);
      hdilPt[i][j]->SetDirectory(rootdir);
      hdilPt[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      //changed binning from 2 GeV to 10 GeV
      hmet[i][j] = new TH1F((string(prefix)+"_hmet_" + string(suffix[i])).c_str(),(string(prefix)+"_met_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      hmet[i][j]->SetDirectory(rootdir);
      hmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hmetPhi[i][j] = new TH1F((string(prefix)+"_hmetPhi_" + string(suffix[i])).c_str(),(string(prefix)+"_metPhi_" + string(suffix[i])).c_str(),
			       50,-1*TMath::Pi(), TMath::Pi());
      hmetPhi[i][j]->SetDirectory(rootdir);
      hmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hmetVsDilepPt[i][j] = new TH2F((string(prefix)+"_hmetVsDilepPt_" + string(suffix[i])).c_str(),
				     (string(prefix)+"_metVsDilepPt_" + string(suffix[i])).c_str(),
				     NBINSPT, MINPT, MAXPT,100,0.,200.);
      hmetVsDilepPt[i][j]->SetDirectory(rootdir);
      hmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      hmetOverPtVsDphi[i][j] = new TH2F((string(prefix)+"_hmetOverPtVsDphi_" + string(suffix[i])).c_str(),
					(string(prefix)+"_metOverPtVsDphi_" + string(suffix[i])).c_str(),
					30,0.,3.,25,0.,TMath::Pi());
      hmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      hmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");


      //pf
      hpfmet[i][j] = new TH1F((string(prefix)+"_hpfmet_" + string(suffix[i])).c_str(),(string(prefix)+"_pfmet_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      hpfmet[i][j]->SetDirectory(rootdir);
      hpfmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      hpfmetPhi[i][j] = new TH1F((string(prefix)+"_hpfmetPhi_" + string(suffix[i])).c_str(),(string(prefix)+"_pfmetPhi_" + string(suffix[i])).c_str(),
				 50,-1*TMath::Pi(), TMath::Pi());
      hpfmetPhi[i][j]->SetDirectory(rootdir);
      hpfmetPhi[i][j]->GetXaxis()->SetTitle("#phi");

      hpfmetVsDilepPt[i][j] = new TH2F((string(prefix)+"_hpfmetVsDilepPt_" + string(suffix[i])).c_str(),
				       (string(prefix)+"_pfmetVsDilepPt_" + string(suffix[i])).c_str(),
				       NBINSPT, MINPT, MAXPT,100,0.,200.);
      hpfmetVsDilepPt[i][j]->SetDirectory(rootdir);
      hpfmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      hpfmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      hpfmetOverPtVsDphi[i][j] = new TH2F((string(prefix)+"_hpfmetOverPtVsDphi_" + string(suffix[i])).c_str(),
					  (string(prefix)+"_pfmetOverPtVsDphi_" + string(suffix[i])).c_str(),
					  30,0.,3.,25,0.,TMath::Pi());
      hpfmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      hpfmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      hpfmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");
    
      //tc
      htcmet[i][j] = new TH1F((string(prefix)+"_htcmet_" + string(suffix[i])).c_str(),(string(prefix)+"_tcmet_" + string(suffix[i])).c_str(),NBINSPT, MINPT, MAXPT);
      htcmet[i][j]->SetDirectory(rootdir);
      htcmet[i][j]->GetXaxis()->SetTitle("MET (GeV)");

      htcmetPhi[i][j] = new TH1F((string(prefix)+"_htcmetPhi_" + string(suffix[i])).c_str(),(string(prefix)+"_tcmetPhi_" + string(suffix[i])).c_str(),
				 50,-1*TMath::Pi(), TMath::Pi());
      htcmetPhi[i][j]->SetDirectory(rootdir);
      htcmetPhi[i][j]->GetXaxis()->SetTitle("#phi");


      hprojmet[i][j]   = new TH1F((string(prefix)+"_hprojmet_" + string(suffix[i])).c_str(),(string(prefix)+"_hprojmet_" + string(suffix[i])).c_str(),
				  40,0.,200.);
      hprojpfmet[i][j] = new TH1F((string(prefix)+"_hprojpfmet_" + string(suffix[i])).c_str(),(string(prefix)+"_hprojpfmet_" + string(suffix[i])).c_str(),
				  40,0.,200.);;
      hprojtcmet[i][j] = new TH1F((string(prefix)+"_hprojtcmet_" + string(suffix[i])).c_str(),(string(prefix)+"_hprojtcmet_" + string(suffix[i])).c_str(),
				  40,0.,200.);;




      htcmetVsDilepPt[i][j] = new TH2F((string(prefix)+"_htcmetVsDilepPt_" + string(suffix[i])).c_str(),
				       (string(prefix)+"_tcmetVsDilepPt_" + string(suffix[i])).c_str(),
				       NBINSPT, MINPT, MAXPT,100,0.,200.);
      htcmetVsDilepPt[i][j]->SetDirectory(rootdir);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("Pt_{ll} (GeV)");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("Met (GeV)");
    
      htcmetOverPtVsDphi[i][j] = new TH2F((string(prefix)+"_htcmetOverPtVsDphi_" + string(suffix[i])).c_str(),
					  (string(prefix)+"_tcmetOverPtVsDphi_" + string(suffix[i])).c_str(),
					  30,0.,3.,25,0.,TMath::Pi());
      htcmetOverPtVsDphi[i][j]->SetDirectory(rootdir);
      htcmetVsDilepPt[i][j]->GetXaxis()->SetTitle("#Delta#Phi");
      htcmetVsDilepPt[i][j]->GetYaxis()->SetTitle("MET/Pt_{ll}");
    
    

      hdphillvsmll[i][j] = new TH2F((string(prefix)+"_dphillvsmll_" + string(suffix[i])).c_str(),
				    (string(prefix)+"_dphillvsmll_" + string(suffix[i])).c_str(),
				    100,10.,210.,50,0., TMath::Pi());
      hdphillvsmll[i][j]->SetDirectory(rootdir);
      hdphillvsmll[i][j]->GetXaxis()->SetTitle("Mass_{ll} (GeV)");
      hdphillvsmll[i][j]->GetYaxis()->SetTitle("#delta#phi_{ll}");

      hptJet1[i][j] = new TH1F((string(prefix)+"_hptJet1_" + string(suffix[i])).c_str(),(string(prefix)+"_ptJet1_" + string(suffix[i])).c_str(),
			       100, 0., 300.);
      hptJet1[i][j]->SetDirectory(rootdir);
      hptJet1[i][j]->GetXaxis()->SetTitle("Pt (GeV)");

      hptJet2[i][j] = new TH1F((string(prefix)+"_hptJet2_" + string(suffix[i])).c_str(),(string(prefix)+"_ptJet2_" + string(suffix[i])).c_str(),
			       100, 0., 300.);
      hptJet2[i][j]->SetDirectory(rootdir);
      hptJet2[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
  
      hptJet3[i][j] = new TH1F((string(prefix)+"_hptJet3_" + string(suffix[i])).c_str(),(string(prefix)+"_ptJet3_" + string(suffix[i])).c_str(),
			       100, 0., 300.);
      hptJet3[i][j]->SetDirectory(rootdir);
      hptJet3[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
    
      hptJet4[i][j] = new TH1F((string(prefix)+"_hptJet4_" + string(suffix[i])).c_str(),(string(prefix)+"_ptJet4_" + string(suffix[i])).c_str(),
			       100, 0., 300.);
      hptJet4[i][j]->SetDirectory(rootdir);
      hptJet4[i][j]->GetXaxis()->SetTitle("Pt (GeV)");
    
      hetaJet1[i][j] = new TH1F((string(prefix)+"_hetaJet1_" + string(suffix[i])).c_str(),(string(prefix)+"_etaJet1_" + string(suffix[i])).c_str(),
				50, -4., 4.);
      hetaJet1[i][j]->SetDirectory(rootdir);
      hetaJet1[i][j]->GetXaxis()->SetTitle("#eta");

      hetaJet2[i][j] = new TH1F((string(prefix)+"_hetaJet2_" + string(suffix[i])).c_str(),(string(prefix)+"_etaJet2_" + string(suffix[i])).c_str(),
				50, -4., 4.);
      hetaJet2[i][j]->SetDirectory(rootdir);
      hetaJet2[i][j]->GetXaxis()->SetTitle("#eta");
 
      hetaJet3[i][j] = new TH1F((string(prefix)+"_hetaJet3_" + string(suffix[i])).c_str(),(string(prefix)+"_etaJet3_" + string(suffix[i])).c_str(),
				50, -4., 4.);
      hetaJet3[i][j]->SetDirectory(rootdir);
      hetaJet3[i][j]->GetXaxis()->SetTitle("#eta");
    
      hetaJet4[i][j] = new TH1F((string(prefix)+"_hetaJet4_" + string(suffix[i])).c_str(),(string(prefix)+"_etaJet4_" + string(suffix[i])).c_str(),
				50, -4., 4.);
      hetaJet4[i][j]->SetDirectory(rootdir);
      hetaJet4[i][j]->GetXaxis()->SetTitle("#eta");
    
      heleSumPt[i][j] = new TH1F((string(prefix)+"_heleSumPt_" + string(suffix[i])).c_str(),(string(prefix)+"_heleSumPt_" + string(suffix[i])).c_str(),
				 100, 0., 25.);
      heleSumPt[i][j]->SetDirectory(rootdir);
      heleSumPt[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuSumPt[i][j] = new TH1F((string(prefix)+"_hmuSumPt_" + string(suffix[i])).c_str(),(string(prefix)+"_hmuSumPt_" + string(suffix[i])).c_str(),
				100, 0., 25.);
      hmuSumPt[i][j]->SetDirectory(rootdir);
      hmuSumPt[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuSumIso[i][j] = new TH1F((string(prefix)+"_hmuSumIso_" + string(suffix[i])).c_str(),(string(prefix)+"_hmuSumIso_" + string(suffix[i])).c_str(),
				 100, 0., 25.);
      hmuSumIso[i][j]->SetDirectory(rootdir);
      hmuSumIso[i][j]->GetXaxis()->SetTitle("#SigmaPt");
      helSumIso[i][j] = new TH1F((string(prefix)+"_helSumIso_" + string(suffix[i])).c_str(),(string(prefix)+"_helSumIso_" + string(suffix[i])).c_str(),
				 100, 0., 25.);
      helSumIso[i][j]->SetDirectory(rootdir);
      helSumIso[i][j]->GetXaxis()->SetTitle("#SigmaPt");
    
      hmuRelIso[i][j] = new TH1F((string(prefix)+"_hmuRelIso_" + string(suffix[i])).c_str(),(string(prefix)+"_hmuRelIso_" + string(suffix[i])).c_str(),
				 NBINSRELISO, RELISOMIN, RELISOMAX);
      hmuRelIso[i][j]->SetDirectory(rootdir);
      helRelIso[i][j] = new TH1F((string(prefix)+"_helRelIso_" + string(suffix[i])).c_str(),(string(prefix)+"_helRelIso_" + string(suffix[i])).c_str(),
				 NBINSRELISO, RELISOMIN, RELISOMAX);
      helRelIso[i][j]->SetDirectory(rootdir);
      // tracker
      hmuRelIsoTrack[i][j] = new TH1F((string(prefix)+"_hmuRelIsoTrack_" + string(suffix[i])).c_str(),(string(prefix)+"_hmuRelIsoTrack_" + string(suffix[i])).c_str(),
				      NBINSRELISO, RELISOMIN, RELISOMAX);
      hmuRelIsoTrack[i][j]->SetDirectory(rootdir);
      helRelIsoTrack[i][j] = new TH1F((string(prefix)+"_helRelIsoTrack_" + string(suffix[i])).c_str(),(string(prefix)+"_helRelIsoTrack_" + string(suffix[i])).c_str(),
				      NBINSRELISO, RELISOMIN, RELISOMAX);
      helRelIsoTrack[i][j]->SetDirectory(rootdir);
      // calorimeter
      hmuRelIsoCalo[i][j] = new TH1F((string(prefix)+"_hmuRelIsoCalo_" + string(suffix[i])).c_str(),(string(prefix)+"_hmuRelIsoCalo_" + string(suffix[i])).c_str(),
				     NBINSRELISO, RELISOMIN, RELISOMAX);
      hmuRelIsoCalo[i][j]->SetDirectory(rootdir);
      helRelIsoCalo[i][j] = new TH1F((string(prefix)+"_helRelIsoCalo_" + string(suffix[i])).c_str(),(string(prefix)+"_helRelIsoCalo_" + string(suffix[i])).c_str(),
				     NBINSRELISO, RELISOMIN, RELISOMAX);
      helRelIsoCalo[i][j]->SetDirectory(rootdir);


      helIsoTrack[i][j] = new TH1F((string(prefix)+"_helIsoTrack_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoTrack_" + string(suffix[i])).c_str(),
					 40,0,40);
      helIsoEcal[i][j]  = new TH1F((string(prefix)+"_helIsoEcal_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoEcal_" + string(suffix[i])).c_str(),
					 30,0,30);
      helIsoHcal[i][j]  = new TH1F((string(prefix)+"_helIsoHcal_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoHcal_" + string(suffix[i])).c_str(),
					 20,0,20);

      helIsoTrackb[i][j] = new TH1F((string(prefix)+"_helIsoTrackb_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoTrackb_" + string(suffix[i])).c_str(),
					 40,0,40);
      helIsoEcalb[i][j]  = new TH1F((string(prefix)+"_helIsoEcalb_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoEcalb_" + string(suffix[i])).c_str(),
					 30,0,30);
      helIsoHcalb[i][j]  = new TH1F((string(prefix)+"_helIsoHcalb_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoHcalb_" + string(suffix[i])).c_str(),
					 20,0,20);

      helIsoTracke[i][j] = new TH1F((string(prefix)+"_helIsoTracke_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoTracke_" + string(suffix[i])).c_str(),
					 40,0,40);
      helIsoEcale[i][j]  = new TH1F((string(prefix)+"_helIsoEcale_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoEcale_" + string(suffix[i])).c_str(),
					 30,0,30);
      helIsoHcale[i][j]  = new TH1F((string(prefix)+"_helIsoHcale_" + string(suffix[i])).c_str(),(string(prefix)+"_helIsoHcale_" + string(suffix[i])).c_str(),
					 20,0,20);
      
    
      hmt2[i][j] = new TH1F((string(prefix) + "_mt2_" + string(suffix[i])).c_str(),(string(prefix) + "_mt2_", string(suffix[i])).c_str(),
			    10,0.,500.);
      hmt2[i][j]->SetDirectory(rootdir);
      
      hmt2J[i][j] = new TH1F((string(prefix) + "_mt2J_" + string(suffix[i])).c_str(),(string(prefix) + "_mt2J_", string(suffix[i])).c_str(),
			     100,0.,500.);
      hmt2J[i][j]->SetDirectory(rootdir);
       

      if (j==0){
	hnJet[i]->Sumw2();
	hnJetWW[i]->Sumw2();
	hnJetWO[i]->Sumw2();
	hnJetWOSemilep[i]->Sumw2();
	hnJetWOOther[i]->Sumw2();
	hnJetOO[i]->Sumw2();
	hnd0Corr[i]->Sumw2();
	hnTrueFVFit[i]->Sumw2();
	hnFVFit[i]->Sumw2();
	hnd0TrueCorr[i]->Sumw2();
	hnEleIsolation[i]->Sumw2();
	hnEleTrueIsolation[i]->Sumw2();
	hnMuIsolation[i]->Sumw2();
	hnMuTrueIsolation[i]->Sumw2();
        hnMeff[i]->Sumw2();
	hnJetinZwindow[i]->Sumw2();
	hnJetoutZwindow[i]->Sumw2();
        hnSumptj[i]->Sumw2();
        hnMET[i]->Sumw2();
      }
      helePt[i][j]->Sumw2();
      hmuPt[i][j]->Sumw2();
      hmuPtFromSilicon[i][j]->Sumw2();
      hminLepPt[i][j]->Sumw2();
      hmaxLepPt[i][j]->Sumw2();
      helePhi[i][j]->Sumw2();
      hmuPhi[i][j]->Sumw2();
      hdphiLep[i][j]->Sumw2();
      heleEta[i][j]->Sumw2();
      hmuEta[i][j]->Sumw2();
      heled0BS[i][j]->Sumw2();
      heled0PV[i][j]->Sumw2();
      hmud0BS[i][j]->Sumw2();
      hmud0PV[i][j]->Sumw2();
      heleEmaxOE5x5[i][j]->Sumw2();

      helIsoTrack[i][j]->Sumw2();
      helIsoEcal[i][j]->Sumw2(); 
      helIsoHcal[i][j]->Sumw2(); 

      helIsoTrackb[i][j]->Sumw2();
      helIsoEcalb[i][j]->Sumw2(); 
      helIsoHcalb[i][j]->Sumw2(); 
      
      helIsoTracke[i][j]->Sumw2();
      helIsoEcale[i][j]->Sumw2(); 
      helIsoHcale[i][j]->Sumw2(); 


      hdilMass[i][j]->Sumw2();
      hdilMassTightWindow[i][j]->Sumw2();
      hdilPt[i][j]->Sumw2();
      hmet[i][j]->Sumw2();
      hmetPhi[i][j]->Sumw2();
      hpfmet[i][j]->Sumw2();
      hpfmetPhi[i][j]->Sumw2();
      htcmet[i][j]->Sumw2();
      htcmetPhi[i][j]->Sumw2();
      hptJet1[i][j]->Sumw2();
      hptJet2[i][j]->Sumw2();
      hptJet3[i][j]->Sumw2();
      hptJet4[i][j]->Sumw2();
      hetaJet1[i][j]->Sumw2();
      hetaJet2[i][j]->Sumw2();
      hetaJet3[i][j]->Sumw2();
      hetaJet4[i][j]->Sumw2();
      heleSumPt[i][j]->Sumw2();
      hmuSumPt[i][j]->Sumw2();

      hmuSumIso[i][j]->Sumw2();
      helSumIso[i][j]->Sumw2();
      hmuRelIso[i][j]->Sumw2();
      helRelIso[i][j]->Sumw2();
      hmuRelIsoTrack[i][j]->Sumw2();
      helRelIsoTrack[i][j]->Sumw2();
      hmuRelIsoCalo[i][j]->Sumw2();
      helRelIsoCalo[i][j]->Sumw2();
      
      hmt2[i][j]->Sumw2();
      hmt2J[i][j]->Sumw2();

    }
  }//channel loop
}//CMS2::bookHistos()


