//-----------------------------------------------------------------------------
// File:        pmssmanalyzer.cc
// Description: Analyzer for ntuples created by Mkntuple
// Created:     Fri Dec  3 11:22:14 2010 by mkntanalyzer.py
// Author:      Sanjay Padhi
// $Revision: 1.19 $
//-----------------------------------------------------------------------------
#include "pmssmanalyzerosleptons.h"

#ifdef PROJECT_NAME
#include "PhysicsTools/Mkntuple/interface/pdg.h"
#else
#include "pdg.h"
#endif

using namespace std;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line
  
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);
  
  // Get names of ntuple files to be processed and open chain of ntuples
  
  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Analysis GEN");
  if ( !stream.good() ) error("unable to open ntuple file(s)");
  
  // Get number of events to be read

  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;
  
  // Select variables to be read

  selectVariables(stream);

  //---------------------------------------------------------------------------
  // Book histograms etc.
  //---------------------------------------------------------------------------
  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the following
  // line

  TApplication app("analyzer", &argc, argv);
  
  histogramFile hfile(cmdline.histfilename);
  
  // Histograms
  // Jets
  TH1F *h_recoSumJet    = new TH1F("recoSumJet",    "recoSumJet",     100, 0., 2000.);
  TH1F *Njet_ee       = new TH1F("Njet_ee",       "Njet_ee",    0, 0, 5);
  TH1F *Njet_mm       = new TH1F("Njet_mm",       "Njet_mm",    0, 0, 5);
  TH1F *Njet_em       = new TH1F("Njet_em",       "Njet_em",    0, 0, 5);
  TH1F *Njet_all       = new TH1F("Njet_all",       "Njet_all",    0, 0, 5);

  // MET
  TH1F *h_recoMET    = new TH1F("recoMET",    "recoMET",     100, 0., 2000.);

  // Leptons

 // Generic variables:
//  double xsect = 4.888; // LM1
//  double xsect = 38.93; // LM0

  double xsect = 1.0;
  double lumi = 35;
  double weight = (xsect*lumi)/nevents;
  double wt = 1;

  // Cut counts and names:
  int ncut = 5; // number of cuts
  int cutc[ncut]; // cut count
  for (int i=0; i<ncut; i++) {
    cutc[i] = 0;
  }
  
  string cutn[ncut];
  
  
  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------
  
  for (int entry=0; entry < nevents; ++entry)
    {
      if (entry % 1000 == 0 ) cout << "Event: " << entry << endl;
      // Read event into memory
      stream.read(entry);
      // Fill the objects
      fillObjects();
      
      // Leptons
      vector< pair <TLorentzVector, int > > hypothesep4;

      /*
      // Gen Level
      for (unsigned int i = 0; i < Particle.size(); i++) {
	if (Particle[i].PT < 10) continue;
	if (!(fabs(Particle[i].Eta) < 2.4)) continue;
	if (Particle[i].Status != 3) continue;
	if (abs(Particle[i].PID) == 11 || abs(Particle[i].PID) == 13) {
	  pair<float,int> lep(Particle[i].PT, Particle[i].PID); 
	  hypothese.push_back(lep);
	}
      }

      
      for (unsigned int i = 0; i < Particle.size(); i++) {
	if (Particle[i].PT < 10) continue;
	if (!(fabs(Particle[i].Eta) < 2.4)) continue;
	if (Particle[i].Status != 1) continue;
	if (abs(Particle[i].PID) == 15) {
	  if (abs(Particle[i].D1) == 11 || abs(Particle[i].D2) == 13) cout << Particle[i].D1 << " particles " << Particle[i].D2 << endl;
	}	
      }
      */

      // reco level-Electrons
      std::vector<Electron_s> vElectron;
      
      for (unsigned int i=0; i<Electron.size(); i++) {
	if (!(Electron[i].PT > 10)) continue;
	if (!(fabs(Electron[i].Eta) < 2.4)) continue;
	if (!(Electron[i].IsolFlag == 1) ) continue;
	vElectron.push_back(Electron[i]);
	
	// Hypothesis
        TLorentzVector ele(Electron[i].Px, Electron[i].Py, Electron[i].Pz, Electron[i].E);
        pair<TLorentzVector,int> lepp4(ele, (Electron[i].Charge)*11);
        hypothesep4.push_back(lepp4);
      }
      
      // reco level Muons
      std::vector<Muon_s> vMuon;
      
      for (unsigned int i=0; i<Muon.size(); i++) {
	if (!(Muon[i].PT > 10)) continue;
	if (!(fabs(Muon[i].Eta) < 2.4)) continue;
	if (!(Muon[i].IsolFlag == 1) ) continue;
	vMuon.push_back(Muon[i]);
	
	// Hypothesis
        TLorentzVector mu(Muon[i].Px, Muon[i].Py, Muon[i].Pz, Muon[i].E);
        pair<TLorentzVector,int> lepp4(mu, (Muon[i].Charge)*13);
        hypothesep4.push_back(lepp4);

      }

	  
      // Jets 
      float sumpt = 0.0;
      int njet = 0;
      
      std::vector<Jet_s> sJet;
      for (unsigned int i=0; i<Jet.size(); i++) {
	if (!(Jet[i].PT > 30)) continue;
	if (!(fabs(Jet[i].Eta) < 2.5)) continue;
	sumpt += Jet[i].PT;
	njet++;
	sJet.push_back(Jet[i]);
      }
      
      // HT and MHT:
      double HT = 0;
      double MHx = 0;
      double MHy = 0;
      for (unsigned int i=0; i<sJet.size(); i++) {
	HT += sJet[i].PT;
	MHx += sJet[i].Px;
	MHy += sJet[i].Py;
      }
//      double MHT = sqrt(MHx*MHx + MHy*MHy);
      double MET = ETmis[0].ET;
      double ysig = MET/sqrt(sumpt);
      
      // ---------------------
      // -- Event Selection --
      // ---------------------
      
      // Hypotheses
      bool osdil = false;
      double invmass = 0.;
      int index1 = -9;
      int index2 = -9;

      for (unsigned int i = 0; i < hypothesep4.size(); i++) {
        pair<TLorentzVector,int> lep1 = hypothesep4[i];
        for(unsigned int j = i+1; j < hypothesep4.size(); j++) {
          pair<TLorentzVector,int> lep2 = hypothesep4[j];
          if(TMath::Max(lep1.first.Perp(),lep2.first.Perp()) < 20) continue;
          if(TMath::Min(lep1.first.Perp(),lep2.first.Perp()) < 10) continue;
          if (lep1.second * lep2.second > 0) continue;
          TLorentzVector ptot = lep1.first + lep2.first;
          // Z Veto
          if (abs(lep1.second) == 11 && abs(lep2.second) == 11 && ptot.M() > 76 && ptot.M() < 106) continue;
          if (abs(lep1.second) == 13 && abs(lep2.second) == 13 && ptot.M() > 76 && ptot.M() < 106) continue;
          float leppt = lep1.first.Perp() + lep2.first.Perp();
          if (leppt > invmass) {
            index1 = i;
            index2 = j;
            invmass = leppt;
          }
        }
      }


     if (invmass > 10) {
          pair<TLorentzVector,int> lep1 = hypothesep4[index1];
          pair<TLorentzVector,int> lep2 = hypothesep4[index2];
          TLorentzVector ptot = lep1.first + lep2.first;
          if (ptot.M() < 10) continue;
          if (ysig < 8.5) continue;
          osdil = true;
          // Use efficiency model for the generator leptons 
          // wt = effmodel(lep1.second, lep1.first.Perp())*effmodel(lep2.second, lep2.first.Perp());
          wt = 1.0;
          if (MET < 50 ) continue;
          if (sumpt < 300) continue;
          if (njet < 2) continue;

          h_recoMET->Fill(MET, weight*wt);
          h_recoSumJet->Fill(sumpt, weight*wt);
          Njet_all->Fill(TMath::Min(njet, 4), weight*wt);
          if (abs(lep1.second) == 11 && abs(lep2.second) ==11 ) Njet_ee->Fill(TMath::Min(njet, 4), weight*wt);
          else if (abs(lep1.second) == 13 && abs(lep2.second) ==13 ) Njet_mm->Fill(TMath::Min(njet, 4), weight*wt);
          else Njet_em->Fill(TMath::Min(njet, 4), weight*wt);
     }
      // if ( !SUSY ) continue;
      // Cut flows
      cutc[0]++;
      if (MET < 50 ) continue;  
      cutc[1]++;
      if (sumpt < 300) continue;
      cutc[2]++;
      if (njet < 2) continue;
      cutc[3]++;
      if (!osdil) continue;
      cutc[4]++;
    }
      
      // Cut names:
      cutn[0] = "Total";
      cutn[1] = "MET>50";
      cutn[2] = "HT>300";
      cutn[3] = "nj>=2";
      cutn[4] = "osnlep";

      

      // Print the cutflow:
      printf("\n +++ Selection results for %.0f pb-1:\n\n", lumi);
      printf("%4s %20s   %10s %10s %10s\n",
	     " ", "Cut", "Nevents", "Abs eff", "Rel eff");
      printf("%4s %20s   %10s %10s %10s\n",
	     " ", " ", " ", "Ni/N0", "Ni/Ni-1");
      printf("%-4s %20s : %10.2f %10.2f %10.2f\n",
	     "N0", "Total", cutc[0]*weight, (cutc[0]*100)/double(cutc[0]), 
	     (cutc[0]*100)/double(cutc[0]));
      for (int i=1; i<ncut; i++) {
	printf("N%-3d %20s : %10.2f %10.2f %10.2f\n",
	       i, cutn[i].c_str(), cutc[i]*weight*wt, (cutc[i]*100)/double(cutc[0]), 
	       (cutc[i]*100)/double(cutc[i-1]) );
      }
    
     std::cout << "---- pMSSM Results with xsec=1pb ---- " << filenames[0] << "  " << cutc[ncut-1]*weight << std::endl; 
  stream.close();
  hfile.close();
  return 0;
}
