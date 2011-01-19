//-----------------------------------------------------------------------------
// File:        pmssmanalyzercms.cc
// Description: Analyzer for ntuples created by Mkntuple
// Created:     Sun Jan  9 11:27:57 2011 by mkntanalyzer.py
// Author:      Sezen Sekmen
// $Revision: 1.19 $
//-----------------------------------------------------------------------------
#include "pmssmanalyzerdijet.h"

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

  TH1F *h_HT = new TH1F("h_HT", "h_HT", 50, 350, 2500);
  TH1F *h_aT = new TH1F("h_aT", "h_aT", 50, 0.2, 1.5);
  TH1F *h_dphis = new TH1F("h_dphis", "h_dphis", 50, 0, 3.142);

  // Generic variables:
//  double xsect = 4.888;
  double xsect = 1.0;
  double lumi = 35;
  double weight = (xsect*lumi)/nevents;

  cout << "Weight: " << weight << endl;

  // Cut counts and names:
  int ncut = 15; // number of cuts
  int cutc[ncut]; // cut count
  for (int i=0; i<ncut; i++) {
    cutc[i] = 0;
  }

  string cutn[ncut];


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry)
    {
      if (entry % 1000 == 0 ) cout << "Event: " << entry << endl;

      // Read event into memory
      stream.read(entry);

      // Fill in the objects:
      fillObjects();

      // ----------------------
      // -- Object Selection --
      // ----------------------

      // Vertices:
      /*
      std::vector<vertex_s> svertex;
      for (unsigned int i=0; i<vertex.size(); i++) {
        if (!(vertex[i].isFake==0)) continue;
        if (!(vertex[i].ndof > 4)) continue;
        if (!(fabs(vertex[i].z) < 24)) continue;
        if (!(vertex[i].position_Rho < 2.0)) continue;
        svertex.push_back(vertex[i]);
      }
      */

      // Jets - selected:
      std::vector<Jet_s> sJet;
      std::vector<Jet_s> sJet30;
      for (unsigned int i=0; i<Jet.size(); i++) {
        if (!(Jet[i].PT > 30)) continue;
        sJet30.push_back(Jet[i]);
        if (!(Jet[i].PT > 50)) continue;
        if (!(fabs(Jet[i].Eta) < 3.0)) continue;
	/*
	// In CMSSW, but not in Delphes:
	if (fabs(Jet[i].eta) < 2.55)
	  if (!(Jet[i].emEnergyFraction > 0.01)) continue;
	if (fabs(Jet[i].eta) >= 2.55)
	  if (!(Jet[i].emEnergyFraction > -0.9)) continue;
	if (Jet[i].pt > 80)
	  if (!(Jet[i].emEnergyFraction < 1)) continue;
        if (!(Jet[i].JetID_n90Hits >= 2)) continue;
        if (!(Jet[i].JetID_fHPD < 0.98)) continue;
	*/
        sJet.push_back(Jet[i]);
      }

      // Jets - vetoed:
      std::vector<Jet_s> vJet;
      for (unsigned int i=0; i<Jet.size(); i++) {
        if (!(Jet[i].PT > 50)) continue;
        if (!(fabs(Jet[i].Eta) >= 3.0)) continue;
	// In CMSSW, but not in Delphes:
        //if (!(Jet[i].JetID_n90Hits < 2)) continue;
        //if (!(Jet[i].JetID_fHPD >= 0.98)) continue;
        vJet.push_back(Jet[i]);
      }

      // Electrons - vetoed:
      std::vector<Electron_s> vElectron;
      for (unsigned int i=0; i<Electron.size(); i++) {
        if (!(Electron[i].PT > 10)) continue;
	if (!(fabs(Electron[i].Eta) < 2.5)) continue;
	// In CMSSW, but not in Delphes:
	//double reliso = (Electron[i].dr03TkSumPt +
	//               Electron[i].dr03EcalRecHitSumEt +
	//               Electron[i].dr03HcalTowerSumEt) / Electron[i].et;
        //if (!(reliso < 0.15)) continue;
	if (!(Electron[i].IsolFlag == 1) ) continue;
	vElectron.push_back(Electron[i]);
      }

      // Muons - vetoed:
      std::vector<Muon_s> vMuon;
      for (unsigned int i=0; i<Muon.size(); i++) {
        if (!(Muon[i].PT > 10)) continue;
	if (!(fabs(Muon[i].Eta) < 2.5)) continue;
	// In CMSSW, but not in Delphes:
	//if (!(Muon[i].GlobalMuonPromptTight)) continue;
	//double reliso = (Muon[i].trackIso + Muon[i].ecalIso + Muon[i].hcalIso) / Muon[i].pt;
	//if (!(reliso < 0.15)) continue;
	if (!(Muon[i].IsolFlag == 1) ) continue;
	vMuon.push_back(Muon[i]);
      }

      /*
      // In CMSSW, but not in Delphes:
      // Odd Muons - vetoed:
      std::vector<Muon_s> voMuon;
      for (unsigned int i=0; i<Muon.size(); i++) {
        if (!(Muon[i].pt > 10) ) continue;
	if (!(fabs(Muon[i].eta) < 2.5) ) continue;
	if (!(Muon[i].isGlobalMuon) ) continue;
	if (!(Muon[i].GlobalMuonPromptTight == false) ) continue;
	voMuon.push_back(Muon[i]);
      }
      */

      // Photons - vetoed:
      std::vector<Photon_s> vPhoton;
      for (unsigned int i=0; i<Photon.size(); i++) {
        if (!(Photon[i].PT > 25)) continue;
	if (!(fabs(Photon[i].Eta) < 2.5)) continue;
	if (!(Photon[i].EHoverEE < 0.05)) continue;
	// In CMSSW, but not in Delphes:
	//if (!(Photon[i].ecalRecHitSumEtConeDR04 < 4.2+0.006*Photon[i].pt)) continue;
	//if (!(Photon[i].hcalTowerSumEtConeDR04 < 2.2+0.0025*Photon[i].pt)) continue;
	//if (!(Photon[i].hadronicOverEm < 0.05)) continue;
	//if (!(Photon[i].trkSumPtHollowConeDR04 < 2+0.001*Photon[i].pt)) continue;
	//if (fabs(Photon[i].eta) < 1.479)
	//if (!(Photon[i].sigmaIetaIeta < 0.013)) continue;
	//if (fabs(Photon[i].eta) >= 1.479)
	//if (!(Photon[i].sigmaIetaIeta < 0.030)) continue;
	// Says optional in the twiki for these cuts.  Numbers suggest that
	// RA1 guys did not use this.
	//if (!(Photon[i].hasPixelSeed == false)) continue;
	vPhoton.push_back(Photon[i]);
      }


      /*
      // Photon veto from genparticles:
      std::vector<Particle_s> vPhoton;
      for (unsigned int i=0; i<Particle.size(); i++) {
	if (!(Particle[i].Status == 3) ) continue;
	if (!(Particle[i].PID == 22) ) continue;
        if (!(Particle[i].PT > 25)) continue;
	if (!(fabs(Particle[i].Eta) < 2.5)) continue;
	vPhoton.push_back(Particle[i]);
      }
      */


      // ---------------
      // -- Variables --
      // ---------------

      // HT and MHT:
      double HT = 0;
      double MHx = 0;
      double MHy = 0;
      for (unsigned int i=0; i<sJet.size(); i++) {
        HT += sJet[i].PT;
	MHx += sJet[i].Px;
	MHy += sJet[i].Py;
      }
      double MHT = sqrt(MHx*MHx + MHy*MHy);


      // ---------------------
      // -- Event Selection --
      // ---------------------

      // if ( !SUSY ) continue;

      //if (!(svertex.size() >= 1) ) continue;
      cutc[0]++;
      if (!(sJet.size() >= 2) ) continue;
      cutc[1]++;
      if (!(HT > 250) ) continue;
      cutc[2]++;
      if (!(vPhoton.size() == 0) ) continue;
      cutc[3]++;
      if (!(vElectron.size() == 0) ) continue;
      cutc[4]++;
      //if (!(voMuon.size() == 0) ) continue;
      cutc[5]++;
      if (!(vMuon.size() == 0) ) continue;
      cutc[6]++;
      cutc[7]++;
      if (!(vJet.size() == 0) ) continue;
      cutc[8]++;
      if (!(fabs(sJet[0].Eta) < 2.5) ) continue;
      cutc[9]++;
      if (!(fabs(sJet[1].PT) > 100) ) continue;
      cutc[10]++;
      if (!(HT > 350) ) continue;
      cutc[11]++;
      h_HT->Fill(HT, weight);
      std::vector<Jet_s> pj = makepseudojets(sJet);
      double aT = alphaT(pj);
      h_aT->Fill(aT, weight);
      if (!(aT > 0.55) ) continue;
      cutc[12]++;
      double dphistar = deltaphistar(sJet30);
      h_dphis->Fill(dphistar, weight);
      //if (!(dphistar > 0.5) ) continue;
      cutc[13]++;
      if (!(MHT/ETmis[0].ET < 1.25) ) continue;
      cutc[14]++;
    }


  // Cut names:
  cutn[0] = "Total";
  cutn[1] = "nj>=2";
  cutn[2] = "HT>250";
  cutn[3] = "ngamma==0";
  cutn[4] = "ne==0";
  cutn[5] = "noddmu==0";
  cutn[6] = "nmu==0";
  cutn[7] = "nbadmuinj";
  cutn[8] = "noddj";
  cutn[9] = "j1eta<2.5";
  cutn[10] = "j2pt>100";
  cutn[11] = "HT>350";
  cutn[12] = "alphaT>0.55";
  cutn[13] = "Dead ECAL";
  cutn[14] = "MHT/MET";

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
	   i, cutn[i].c_str(), cutc[i]*weight, (cutc[i]*100)/double(cutc[0]), 
	   (cutc[i]*100)/double(cutc[i-1]) );
  }

   std::cout << "---- pMSSM Results with xsec=1pb ---- " << filenames[0] << "  " << cutc[ncut-1]*weight << std::endl;



  stream.close();
  hfile.close();
  return 0;
}
