#ifndef PMSSMANALYZER_H
#define PMSSMANALYZER_H
//-----------------------------------------------------------------------------
// File:        pmssmanalyzer.h
// Description: Analyzer header for ntuples created by Mkntuple
// Created:     Fri Dec  3 11:22:14 2010 by mkntanalyzer.py
// Author:      Sezen Sekmen
// $Revision: 1.19 $
//-----------------------------------------------------------------------------

// -- System

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <set>

#ifdef PROJECT_NAME

// --- CMSSW

#include "PhysicsTools/Mkntuple/interface/treestream.h"
#include "PhysicsTools/Mkntuple/interface/pdg.h"

#else

#include "treestream.h"
#include "pdg.h"

#endif

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
//-----------------------------------------------------------------------------
// -- Utilities
//-----------------------------------------------------------------------------
void
error(std::string message)
{
  std::cout << "** error ** " << message << std::endl;
  exit(0);
}

std::string 
strip(std::string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

std::string
nameonly(std::string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}
//-----------------------------------------------------------------------------
struct skimFile
{
  skimFile(itreestream& stream, std::string filename)
   : filename_(filename),
	 file_(new TFile(filename.c_str(), "recreate")),
	 tree_(stream.tree()->CloneTree(0)),
	 b_weight_(0)
  {
	std::cout << "events will be skimmed to file "
			  << filename_ << std::endl;
	file_->cd();
	hist_ = new TH1F("counts", "", 1,0,1);
	hist_->SetBit(TH1::kCanRebin);
	hist_->SetStats(0);
  }

  void addEvent(double weight=1)
  {
    weight_ = weight;
	file_   = tree_->GetCurrentFile();
	file_->cd();
	
	// dynamically add a weight branch to cloned tree
	if ( b_weight_ == 0 )
      {
	    b_weight_ = tree_->Branch("eventWeight", &weight_, "eventWeight/D");
	  }
	tree_->Fill();
  }

  void count(std::string cond)
  {
    hist_->Fill(cond.c_str(), 1);
  }
  
  void close()
  {
	std::cout << "==> events skimmed to file " << filename_ << std::endl;
	file_ = tree_->GetCurrentFile();
	file_->cd();
	file_->Write("", TObject::kOverwrite);
	hist_->Write("", TObject::kOverwrite);
	file_->Close();
  }

  std::string filename_;  
  TFile* file_;
  TTree* tree_;
  TH1F*  hist_;
  TBranch* b_weight_;
  double     weight_;
};

struct commandLine
{
  std::string progname;
  std::string filelist;
  std::string histfilename;
  std::string skimfilename;
};

struct histogramFile
{
  histogramFile(std::string histfilename)
   : filename_(histfilename),
	 file_(new TFile(filename_.c_str(), "recreate")) 
  {}

  void close()
  {
	std::cout << "==> histograms saved to file " << filename_ << std::endl;
	file_->cd();
	file_->Write("", TObject::kOverwrite);
	file_->ls();
	file_->Close();
  }

  std::string filename_;
  TFile* file_;
};

void
decodeCommandLine(int argc, char** argv, commandLine& cl)
{
  cl.progname = std::string(argv[0]);

  // 1st (optional) argument
  if ( argc > 1 )
	cl.filelist = std::string(argv[1]);
  else
	cl.filelist = std::string("filelist.txt");

  // 2nd (optional) command line argument
  if ( argc > 2 ) 
	cl.histfilename = std::string(argv[2]);
  else
	cl.histfilename = cl.progname + std::string("_histograms");

  // 3rd (optional) command line argument
  if ( argc > 3 ) 
	cl.skimfilename = std::string(argv[3]);
  else
	cl.skimfilename = cl.progname + std::string("_skim");

  // Make sure extension is ".root"
  cl.histfilename = nameonly(cl.histfilename);
  cl.histfilename += ".root";

  cl.skimfilename = nameonly(cl.skimfilename);
  cl.skimfilename += ".root";
}

// Read ntuple filenames from file list

std::vector<std::string>
getFilenames(std::string filelist)
{
  std::ifstream stream(filelist.c_str());
  if ( !stream.good() ) error("unable to open file: " + filelist);

  // Get list of ntuple files to be processed

  std::vector<std::string> v;
  std::string filename;
  while ( stream >> filename )
	if ( strip(filename) != "" ) v.push_back(filename);
  return v;
}
//-----------------------------------------------------------------------------
// -- Declare variables to be read
//-----------------------------------------------------------------------------
std::vector<float>	CaloTower_E(989,0);
std::vector<float>	CaloTower_ET(989,0);
std::vector<float>	CaloTower_E_em(989,0);
std::vector<float>	CaloTower_E_had(989,0);
std::vector<float>	CaloTower_Eta(989,0);
std::vector<float>	CaloTower_Phi(989,0);
int	CaloTower_size;
std::vector<float>	ETmis_ET(3,0);
std::vector<float>	ETmis_Phi(3,0);
std::vector<float>	ETmis_Px(3,0);
std::vector<float>	ETmis_Py(3,0);
int	ETmis_size;
std::vector<int>	Electron_Charge(9,0);
std::vector<float>	Electron_E(9,0);
std::vector<float>	Electron_EHoverEE(9,0);
std::vector<float>	Electron_EtRatio(9,0);
std::vector<float>	Electron_Eta(9,0);
std::vector<float>	Electron_EtaCalo(9,0);
std::vector<int>	Electron_IsolFlag(9,0);
std::vector<float>	Electron_PT(9,0);
std::vector<float>	Electron_Phi(9,0);
std::vector<float>	Electron_PhiCalo(9,0);
std::vector<float>	Electron_Px(9,0);
std::vector<float>	Electron_Py(9,0);
std::vector<float>	Electron_Pz(9,0);
std::vector<float>	Electron_SumEt(9,0);
std::vector<float>	Electron_SumPt(9,0);
int	Electron_size;
std::vector<double>	Event_CouplingQCD(3,0);
std::vector<double>	Event_CouplingQED(3,0);
std::vector<int>	Event_Nparticles(3,0);
std::vector<long>	Event_Number(3,0);
std::vector<int>	Event_ProcessID(3,0);
std::vector<double>	Event_ScalePDF(3,0);
std::vector<double>	Event_Weight(3,0);
int	Event_size;
int	FP420hits;
float	FP420hits_E;
float	FP420hits_S;
float	FP420hits_T;
float	FP420hits_Tx;
float	FP420hits_Ty;
float	FP420hits_X;
float	FP420hits_Y;
float	FP420hits_genE;
float	FP420hits_genEta;
float	FP420hits_genPT;
float	FP420hits_genPhi;
float	FP420hits_genPx;
float	FP420hits_genPy;
float	FP420hits_genPz;
int	FP420hits_pid;
float	FP420hits_q2;
int	FP420hits_side;
int	FP420hits_size;
std::vector<int>	Jet_Btag(27,0);
std::vector<float>	Jet_E(27,0);
std::vector<float>	Jet_EHoverEE(27,0);
std::vector<float>	Jet_Eta(27,0);
std::vector<int>	Jet_NCalo(27,0);
std::vector<int>	Jet_NTracks(27,0);
std::vector<float>	Jet_PT(27,0);
std::vector<float>	Jet_Phi(27,0);
std::vector<float>	Jet_Px(27,0);
std::vector<float>	Jet_Py(27,0);
std::vector<float>	Jet_Pz(27,0);
int	Jet_size;
std::vector<int>	Muon_Charge(9,0);
std::vector<float>	Muon_E(9,0);
std::vector<float>	Muon_EHoverEE(9,0);
std::vector<float>	Muon_EtRatio(9,0);
std::vector<float>	Muon_Eta(9,0);
std::vector<float>	Muon_EtaCalo(9,0);
std::vector<int>	Muon_IsolFlag(9,0);
std::vector<float>	Muon_PT(9,0);
std::vector<float>	Muon_Phi(9,0);
std::vector<float>	Muon_PhiCalo(9,0);
std::vector<float>	Muon_Px(9,0);
std::vector<float>	Muon_Py(9,0);
std::vector<float>	Muon_Pz(9,0);
std::vector<float>	Muon_SumEt(9,0);
std::vector<float>	Muon_SumPt(9,0);
int	Muon_size;
std::vector<float>	Particle_Charge(2537,0);
std::vector<int>	Particle_D1(2537,0);
std::vector<int>	Particle_D2(2537,0);
std::vector<float>	Particle_E(2537,0);
std::vector<float>	Particle_Eta(2537,0);
std::vector<float>	Particle_M(2537,0);
std::vector<int>	Particle_M1(2537,0);
std::vector<int>	Particle_M2(2537,0);
std::vector<int>	Particle_PID(2537,0);
std::vector<float>	Particle_PT(2537,0);
std::vector<float>	Particle_Phi(2537,0);
std::vector<float>	Particle_Px(2537,0);
std::vector<float>	Particle_Py(2537,0);
std::vector<float>	Particle_Pz(2537,0);
std::vector<int>	Particle_Status(2537,0);
std::vector<float>	Particle_T(2537,0);
std::vector<float>	Particle_X(2537,0);
std::vector<float>	Particle_Y(2537,0);
std::vector<float>	Particle_Z(2537,0);
int	Particle_size;
std::vector<float>	Photon_E(49,0);
std::vector<float>	Photon_EHoverEE(49,0);
std::vector<float>	Photon_Eta(49,0);
std::vector<float>	Photon_PT(49,0);
std::vector<float>	Photon_Phi(49,0);
std::vector<float>	Photon_Px(49,0);
std::vector<float>	Photon_Py(49,0);
std::vector<float>	Photon_Pz(49,0);
int	Photon_size;
int	RP220hits;
float	RP220hits_E;
float	RP220hits_S;
float	RP220hits_T;
float	RP220hits_Tx;
float	RP220hits_Ty;
float	RP220hits_X;
float	RP220hits_Y;
float	RP220hits_genE;
float	RP220hits_genEta;
float	RP220hits_genPT;
float	RP220hits_genPhi;
float	RP220hits_genPx;
float	RP220hits_genPy;
float	RP220hits_genPz;
int	RP220hits_pid;
float	RP220hits_q2;
int	RP220hits_side;
int	RP220hits_size;
std::vector<float>	TauJet_Charge(7,0);
std::vector<float>	TauJet_E(7,0);
std::vector<float>	TauJet_EHoverEE(7,0);
std::vector<float>	TauJet_Eta(7,0);
std::vector<int>	TauJet_NCalo(7,0);
std::vector<int>	TauJet_NTracks(7,0);
std::vector<float>	TauJet_PT(7,0);
std::vector<float>	TauJet_Phi(7,0);
std::vector<float>	TauJet_Px(7,0);
std::vector<float>	TauJet_Py(7,0);
std::vector<float>	TauJet_Pz(7,0);
int	TauJet_size;
std::vector<float>	Tracks_Charge(357,0);
std::vector<float>	Tracks_E(357,0);
std::vector<float>	Tracks_Eta(357,0);
std::vector<float>	Tracks_EtaOuter(357,0);
std::vector<float>	Tracks_PT(357,0);
std::vector<float>	Tracks_Phi(357,0);
std::vector<float>	Tracks_PhiOuter(357,0);
std::vector<float>	Tracks_Px(357,0);
std::vector<float>	Tracks_Py(357,0);
std::vector<float>	Tracks_Pz(357,0);
std::vector<float>	Tracks_Vx(357,0);
std::vector<float>	Tracks_Vy(357,0);
std::vector<float>	Tracks_Vz(357,0);
int	Tracks_size;
std::vector<float>	ZDChits_E(9,0);
std::vector<float>	ZDChits_T(9,0);
std::vector<float>	ZDChits_genE(9,0);
std::vector<float>	ZDChits_genEta(9,0);
std::vector<float>	ZDChits_genPT(9,0);
std::vector<float>	ZDChits_genPhi(9,0);
std::vector<float>	ZDChits_genPx(9,0);
std::vector<float>	ZDChits_genPy(9,0);
std::vector<float>	ZDChits_genPz(9,0);
std::vector<int>	ZDChits_hadronic_hit(9,0);
std::vector<int>	ZDChits_pid(9,0);
std::vector<int>	ZDChits_side(9,0);
int	ZDChits_size;


//-----------------------------------------------------------------------------
// -- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("CaloTower.E", CaloTower_E);
  stream.select("CaloTower.ET", CaloTower_ET);
  stream.select("CaloTower.E_em", CaloTower_E_em);
  stream.select("CaloTower.E_had", CaloTower_E_had);
  stream.select("CaloTower.Eta", CaloTower_Eta);
  stream.select("CaloTower.Phi", CaloTower_Phi);
  stream.select("CaloTower_size", CaloTower_size);
  stream.select("ETmis.ET", ETmis_ET);
  stream.select("ETmis.Phi", ETmis_Phi);
  stream.select("ETmis.Px", ETmis_Px);
  stream.select("ETmis.Py", ETmis_Py);
  stream.select("ETmis_size", ETmis_size);
  stream.select("Electron.Charge", Electron_Charge);
  stream.select("Electron.E", Electron_E);
  stream.select("Electron.EHoverEE", Electron_EHoverEE);
  stream.select("Electron.EtRatio", Electron_EtRatio);
  stream.select("Electron.Eta", Electron_Eta);
  stream.select("Electron.EtaCalo", Electron_EtaCalo);
  stream.select("Electron.IsolFlag", Electron_IsolFlag);
  stream.select("Electron.PT", Electron_PT);
  stream.select("Electron.Phi", Electron_Phi);
  stream.select("Electron.PhiCalo", Electron_PhiCalo);
  stream.select("Electron.Px", Electron_Px);
  stream.select("Electron.Py", Electron_Py);
  stream.select("Electron.Pz", Electron_Pz);
  stream.select("Electron.SumEt", Electron_SumEt);
  stream.select("Electron.SumPt", Electron_SumPt);
  stream.select("Electron_size", Electron_size);
  stream.select("Event.CouplingQCD", Event_CouplingQCD);
  stream.select("Event.CouplingQED", Event_CouplingQED);
  stream.select("Event.Nparticles", Event_Nparticles);
  stream.select("Event.Number", Event_Number);
  stream.select("Event.ProcessID", Event_ProcessID);
  stream.select("Event.ScalePDF", Event_ScalePDF);
  stream.select("Event.Weight", Event_Weight);
  stream.select("Event_size", Event_size);
  stream.select("FP420hits", FP420hits);
  stream.select("FP420hits.E", FP420hits_E);
  stream.select("FP420hits.S", FP420hits_S);
  stream.select("FP420hits.T", FP420hits_T);
  stream.select("FP420hits.Tx", FP420hits_Tx);
  stream.select("FP420hits.Ty", FP420hits_Ty);
  stream.select("FP420hits.X", FP420hits_X);
  stream.select("FP420hits.Y", FP420hits_Y);
  stream.select("FP420hits.genE", FP420hits_genE);
  stream.select("FP420hits.genEta", FP420hits_genEta);
  stream.select("FP420hits.genPT", FP420hits_genPT);
  stream.select("FP420hits.genPhi", FP420hits_genPhi);
  stream.select("FP420hits.genPx", FP420hits_genPx);
  stream.select("FP420hits.genPy", FP420hits_genPy);
  stream.select("FP420hits.genPz", FP420hits_genPz);
  stream.select("FP420hits.pid", FP420hits_pid);
  stream.select("FP420hits.q2", FP420hits_q2);
  stream.select("FP420hits.side", FP420hits_side);
  stream.select("FP420hits_size", FP420hits_size);
  stream.select("Jet.Btag", Jet_Btag);
  stream.select("Jet.E", Jet_E);
  stream.select("Jet.EHoverEE", Jet_EHoverEE);
  stream.select("Jet.Eta", Jet_Eta);
  stream.select("Jet.NCalo", Jet_NCalo);
  stream.select("Jet.NTracks", Jet_NTracks);
  stream.select("Jet.PT", Jet_PT);
  stream.select("Jet.Phi", Jet_Phi);
  stream.select("Jet.Px", Jet_Px);
  stream.select("Jet.Py", Jet_Py);
  stream.select("Jet.Pz", Jet_Pz);
  stream.select("Jet_size", Jet_size);
  stream.select("Muon.Charge", Muon_Charge);
  stream.select("Muon.E", Muon_E);
  stream.select("Muon.EHoverEE", Muon_EHoverEE);
  stream.select("Muon.EtRatio", Muon_EtRatio);
  stream.select("Muon.Eta", Muon_Eta);
  stream.select("Muon.EtaCalo", Muon_EtaCalo);
  stream.select("Muon.IsolFlag", Muon_IsolFlag);
  stream.select("Muon.PT", Muon_PT);
  stream.select("Muon.Phi", Muon_Phi);
  stream.select("Muon.PhiCalo", Muon_PhiCalo);
  stream.select("Muon.Px", Muon_Px);
  stream.select("Muon.Py", Muon_Py);
  stream.select("Muon.Pz", Muon_Pz);
  stream.select("Muon.SumEt", Muon_SumEt);
  stream.select("Muon.SumPt", Muon_SumPt);
  stream.select("Muon_size", Muon_size);
  stream.select("Particle.Charge", Particle_Charge);
  stream.select("Particle.D1", Particle_D1);
  stream.select("Particle.D2", Particle_D2);
  stream.select("Particle.E", Particle_E);
  stream.select("Particle.Eta", Particle_Eta);
  stream.select("Particle.M", Particle_M);
  stream.select("Particle.M1", Particle_M1);
  stream.select("Particle.M2", Particle_M2);
  stream.select("Particle.PID", Particle_PID);
  stream.select("Particle.PT", Particle_PT);
  stream.select("Particle.Phi", Particle_Phi);
  stream.select("Particle.Px", Particle_Px);
  stream.select("Particle.Py", Particle_Py);
  stream.select("Particle.Pz", Particle_Pz);
  stream.select("Particle.Status", Particle_Status);
  stream.select("Particle.T", Particle_T);
  stream.select("Particle.X", Particle_X);
  stream.select("Particle.Y", Particle_Y);
  stream.select("Particle.Z", Particle_Z);
  stream.select("Particle_size", Particle_size);
  stream.select("Photon.E", Photon_E);
  stream.select("Photon.EHoverEE", Photon_EHoverEE);
  stream.select("Photon.Eta", Photon_Eta);
  stream.select("Photon.PT", Photon_PT);
  stream.select("Photon.Phi", Photon_Phi);
  stream.select("Photon.Px", Photon_Px);
  stream.select("Photon.Py", Photon_Py);
  stream.select("Photon.Pz", Photon_Pz);
  stream.select("Photon_size", Photon_size);
  stream.select("RP220hits", RP220hits);
  stream.select("RP220hits.E", RP220hits_E);
  stream.select("RP220hits.S", RP220hits_S);
  stream.select("RP220hits.T", RP220hits_T);
  stream.select("RP220hits.Tx", RP220hits_Tx);
  stream.select("RP220hits.Ty", RP220hits_Ty);
  stream.select("RP220hits.X", RP220hits_X);
  stream.select("RP220hits.Y", RP220hits_Y);
  stream.select("RP220hits.genE", RP220hits_genE);
  stream.select("RP220hits.genEta", RP220hits_genEta);
  stream.select("RP220hits.genPT", RP220hits_genPT);
  stream.select("RP220hits.genPhi", RP220hits_genPhi);
  stream.select("RP220hits.genPx", RP220hits_genPx);
  stream.select("RP220hits.genPy", RP220hits_genPy);
  stream.select("RP220hits.genPz", RP220hits_genPz);
  stream.select("RP220hits.pid", RP220hits_pid);
  stream.select("RP220hits.q2", RP220hits_q2);
  stream.select("RP220hits.side", RP220hits_side);
  stream.select("RP220hits_size", RP220hits_size);
  stream.select("TauJet.Charge", TauJet_Charge);
  stream.select("TauJet.E", TauJet_E);
  stream.select("TauJet.EHoverEE", TauJet_EHoverEE);
  stream.select("TauJet.Eta", TauJet_Eta);
  stream.select("TauJet.NCalo", TauJet_NCalo);
  stream.select("TauJet.NTracks", TauJet_NTracks);
  stream.select("TauJet.PT", TauJet_PT);
  stream.select("TauJet.Phi", TauJet_Phi);
  stream.select("TauJet.Px", TauJet_Px);
  stream.select("TauJet.Py", TauJet_Py);
  stream.select("TauJet.Pz", TauJet_Pz);
  stream.select("TauJet_size", TauJet_size);
  stream.select("Tracks.Charge", Tracks_Charge);
  stream.select("Tracks.E", Tracks_E);
  stream.select("Tracks.Eta", Tracks_Eta);
  stream.select("Tracks.EtaOuter", Tracks_EtaOuter);
  stream.select("Tracks.PT", Tracks_PT);
  stream.select("Tracks.Phi", Tracks_Phi);
  stream.select("Tracks.PhiOuter", Tracks_PhiOuter);
  stream.select("Tracks.Px", Tracks_Px);
  stream.select("Tracks.Py", Tracks_Py);
  stream.select("Tracks.Pz", Tracks_Pz);
  stream.select("Tracks.Vx", Tracks_Vx);
  stream.select("Tracks.Vy", Tracks_Vy);
  stream.select("Tracks.Vz", Tracks_Vz);
  stream.select("Tracks_size", Tracks_size);
  stream.select("ZDChits.E", ZDChits_E);
  stream.select("ZDChits.T", ZDChits_T);
  stream.select("ZDChits.genE", ZDChits_genE);
  stream.select("ZDChits.genEta", ZDChits_genEta);
  stream.select("ZDChits.genPT", ZDChits_genPT);
  stream.select("ZDChits.genPhi", ZDChits_genPhi);
  stream.select("ZDChits.genPx", ZDChits_genPx);
  stream.select("ZDChits.genPy", ZDChits_genPy);
  stream.select("ZDChits.genPz", ZDChits_genPz);
  stream.select("ZDChits.hadronic_hit", ZDChits_hadronic_hit);
  stream.select("ZDChits.pid", ZDChits_pid);
  stream.select("ZDChits.side", ZDChits_side);
  stream.select("ZDChits_size", ZDChits_size);

}

//-----------------------------------------------------------------------------
// --- These structs can be filled by calling fillObjects()
// --- after the call to stream.read(...)
//-----------------------------------------------------------------------------
struct CaloTower_s
{
  float	E;
  float	ET;
  float	E_em;
  float	E_had;
  float	Eta;
  float	Phi;
};

std::vector<CaloTower_s> CaloTower(989);

std::ostream& operator<<(std::ostream& os, const CaloTower_s& o)
{
  char r[1024];
  os << "CaloTower" << std::endl;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "ET", (double)o.ET); os << r;
  sprintf(r, "  %-32s: %f\n", "E_em", (double)o.E_em); os << r;
  sprintf(r, "  %-32s: %f\n", "E_had", (double)o.E_had); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  return os;
}

struct ETmis_s
{
  float	ET;
  float	Phi;
  float	Px;
  float	Py;
};

std::vector<ETmis_s> ETmis(3);

std::ostream& operator<<(std::ostream& os, const ETmis_s& o)
{
  char r[1024];
  os << "ETmis" << std::endl;
  sprintf(r, "  %-32s: %f\n", "ET", (double)o.ET); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  return os;
}

struct Electron_s
{
  int	Charge;
  float	E;
  float	EHoverEE;
  float	EtRatio;
  float	Eta;
  float	EtaCalo;
  int	IsolFlag;
  float	PT;
  float	Phi;
  float	PhiCalo;
  float	Px;
  float	Py;
  float	Pz;
  float	SumEt;
  float	SumPt;
};

std::vector<Electron_s> Electron(9);

std::ostream& operator<<(std::ostream& os, const Electron_s& o)
{
  char r[1024];
  os << "Electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Charge", (double)o.Charge); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "EHoverEE", (double)o.EHoverEE); os << r;
  sprintf(r, "  %-32s: %f\n", "EtRatio", (double)o.EtRatio); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "EtaCalo", (double)o.EtaCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "IsolFlag", (double)o.IsolFlag); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "PhiCalo", (double)o.PhiCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  sprintf(r, "  %-32s: %f\n", "SumEt", (double)o.SumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "SumPt", (double)o.SumPt); os << r;
  return os;
}

struct Event_s
{
  double	CouplingQCD;
  double	CouplingQED;
  int	Nparticles;
  long	Number;
  int	ProcessID;
  double	ScalePDF;
  double	Weight;
};

std::vector<Event_s> Event(3);

std::ostream& operator<<(std::ostream& os, const Event_s& o)
{
  char r[1024];
  os << "Event" << std::endl;
  sprintf(r, "  %-32s: %f\n", "CouplingQCD", (double)o.CouplingQCD); os << r;
  sprintf(r, "  %-32s: %f\n", "CouplingQED", (double)o.CouplingQED); os << r;
  sprintf(r, "  %-32s: %f\n", "Nparticles", (double)o.Nparticles); os << r;
  sprintf(r, "  %-32s: %f\n", "Number", (double)o.Number); os << r;
  sprintf(r, "  %-32s: %f\n", "ProcessID", (double)o.ProcessID); os << r;
  sprintf(r, "  %-32s: %f\n", "ScalePDF", (double)o.ScalePDF); os << r;
  sprintf(r, "  %-32s: %f\n", "Weight", (double)o.Weight); os << r;
  return os;
}

struct Jet_s
{
  int	Btag;
  float	E;
  float	EHoverEE;
  float	Eta;
  int	NCalo;
  int	NTracks;
  float	PT;
  float	Phi;
  float	Px;
  float	Py;
  float	Pz;
};

std::vector<Jet_s> Jet(27);

std::ostream& operator<<(std::ostream& os, const Jet_s& o)
{
  char r[1024];
  os << "Jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Btag", (double)o.Btag); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "EHoverEE", (double)o.EHoverEE); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "NCalo", (double)o.NCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "NTracks", (double)o.NTracks); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  return os;
}

struct Muon_s
{
  int	Charge;
  float	E;
  float	EHoverEE;
  float	EtRatio;
  float	Eta;
  float	EtaCalo;
  int	IsolFlag;
  float	PT;
  float	Phi;
  float	PhiCalo;
  float	Px;
  float	Py;
  float	Pz;
  float	SumEt;
  float	SumPt;
};

std::vector<Muon_s> Muon(9);

std::ostream& operator<<(std::ostream& os, const Muon_s& o)
{
  char r[1024];
  os << "Muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Charge", (double)o.Charge); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "EHoverEE", (double)o.EHoverEE); os << r;
  sprintf(r, "  %-32s: %f\n", "EtRatio", (double)o.EtRatio); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "EtaCalo", (double)o.EtaCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "IsolFlag", (double)o.IsolFlag); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "PhiCalo", (double)o.PhiCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  sprintf(r, "  %-32s: %f\n", "SumEt", (double)o.SumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "SumPt", (double)o.SumPt); os << r;
  return os;
}

struct Particle_s
{
  float	Charge;
  int	D1;
  int	D2;
  float	E;
  float	Eta;
  float	M;
  int	M1;
  int	M2;
  int	PID;
  float	PT;
  float	Phi;
  float	Px;
  float	Py;
  float	Pz;
  int	Status;
  float	T;
  float	X;
  float	Y;
  float	Z;
};

std::vector<Particle_s> Particle(2537);

std::ostream& operator<<(std::ostream& os, const Particle_s& o)
{
  char r[1024];
  os << "Particle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Charge", (double)o.Charge); os << r;
  sprintf(r, "  %-32s: %f\n", "D1", (double)o.D1); os << r;
  sprintf(r, "  %-32s: %f\n", "D2", (double)o.D2); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "M", (double)o.M); os << r;
  sprintf(r, "  %-32s: %f\n", "M1", (double)o.M1); os << r;
  sprintf(r, "  %-32s: %f\n", "M2", (double)o.M2); os << r;
  sprintf(r, "  %-32s: %f\n", "PID", (double)o.PID); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  sprintf(r, "  %-32s: %f\n", "Status", (double)o.Status); os << r;
  sprintf(r, "  %-32s: %f\n", "T", (double)o.T); os << r;
  sprintf(r, "  %-32s: %f\n", "X", (double)o.X); os << r;
  sprintf(r, "  %-32s: %f\n", "Y", (double)o.Y); os << r;
  sprintf(r, "  %-32s: %f\n", "Z", (double)o.Z); os << r;
  return os;
}

struct Photon_s
{
  float	E;
  float	EHoverEE;
  float	Eta;
  float	PT;
  float	Phi;
  float	Px;
  float	Py;
  float	Pz;
};

std::vector<Photon_s> Photon(49);

std::ostream& operator<<(std::ostream& os, const Photon_s& o)
{
  char r[1024];
  os << "Photon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "EHoverEE", (double)o.EHoverEE); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  return os;
}

struct TauJet_s
{
  float	Charge;
  float	E;
  float	EHoverEE;
  float	Eta;
  int	NCalo;
  int	NTracks;
  float	PT;
  float	Phi;
  float	Px;
  float	Py;
  float	Pz;
};

std::vector<TauJet_s> TauJet(7);

std::ostream& operator<<(std::ostream& os, const TauJet_s& o)
{
  char r[1024];
  os << "TauJet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Charge", (double)o.Charge); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "EHoverEE", (double)o.EHoverEE); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "NCalo", (double)o.NCalo); os << r;
  sprintf(r, "  %-32s: %f\n", "NTracks", (double)o.NTracks); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  return os;
}

struct Tracks_s
{
  float	Charge;
  float	E;
  float	Eta;
  float	EtaOuter;
  float	PT;
  float	Phi;
  float	PhiOuter;
  float	Px;
  float	Py;
  float	Pz;
  float	Vx;
  float	Vy;
  float	Vz;
};

std::vector<Tracks_s> Tracks(357);

std::ostream& operator<<(std::ostream& os, const Tracks_s& o)
{
  char r[1024];
  os << "Tracks" << std::endl;
  sprintf(r, "  %-32s: %f\n", "Charge", (double)o.Charge); os << r;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "Eta", (double)o.Eta); os << r;
  sprintf(r, "  %-32s: %f\n", "EtaOuter", (double)o.EtaOuter); os << r;
  sprintf(r, "  %-32s: %f\n", "PT", (double)o.PT); os << r;
  sprintf(r, "  %-32s: %f\n", "Phi", (double)o.Phi); os << r;
  sprintf(r, "  %-32s: %f\n", "PhiOuter", (double)o.PhiOuter); os << r;
  sprintf(r, "  %-32s: %f\n", "Px", (double)o.Px); os << r;
  sprintf(r, "  %-32s: %f\n", "Py", (double)o.Py); os << r;
  sprintf(r, "  %-32s: %f\n", "Pz", (double)o.Pz); os << r;
  sprintf(r, "  %-32s: %f\n", "Vx", (double)o.Vx); os << r;
  sprintf(r, "  %-32s: %f\n", "Vy", (double)o.Vy); os << r;
  sprintf(r, "  %-32s: %f\n", "Vz", (double)o.Vz); os << r;
  return os;
}

struct ZDChits_s
{
  float	E;
  float	T;
  float	genE;
  float	genEta;
  float	genPT;
  float	genPhi;
  float	genPx;
  float	genPy;
  float	genPz;
  int	hadronic_hit;
  int	pid;
  int	side;
};

std::vector<ZDChits_s> ZDChits(9);

std::ostream& operator<<(std::ostream& os, const ZDChits_s& o)
{
  char r[1024];
  os << "ZDChits" << std::endl;
  sprintf(r, "  %-32s: %f\n", "E", (double)o.E); os << r;
  sprintf(r, "  %-32s: %f\n", "T", (double)o.T); os << r;
  sprintf(r, "  %-32s: %f\n", "genE", (double)o.genE); os << r;
  sprintf(r, "  %-32s: %f\n", "genEta", (double)o.genEta); os << r;
  sprintf(r, "  %-32s: %f\n", "genPT", (double)o.genPT); os << r;
  sprintf(r, "  %-32s: %f\n", "genPhi", (double)o.genPhi); os << r;
  sprintf(r, "  %-32s: %f\n", "genPx", (double)o.genPx); os << r;
  sprintf(r, "  %-32s: %f\n", "genPy", (double)o.genPy); os << r;
  sprintf(r, "  %-32s: %f\n", "genPz", (double)o.genPz); os << r;
  sprintf(r, "  %-32s: %f\n", "hadronic_hit", (double)o.hadronic_hit); os << r;
  sprintf(r, "  %-32s: %f\n", "pid", (double)o.pid); os << r;
  sprintf(r, "  %-32s: %f\n", "side", (double)o.side); os << r;
  return os;
}


void fillObjects()
{

  CaloTower.resize(CaloTower_E.size());
  for(unsigned int i=0; i < CaloTower_E.size(); ++i)
    {
      CaloTower[i].E	= CaloTower_E[i];
      CaloTower[i].ET	= CaloTower_ET[i];
      CaloTower[i].E_em	= CaloTower_E_em[i];
      CaloTower[i].E_had	= CaloTower_E_had[i];
      CaloTower[i].Eta	= CaloTower_Eta[i];
      CaloTower[i].Phi	= CaloTower_Phi[i];
    }

  ETmis.resize(ETmis_ET.size());
  for(unsigned int i=0; i < ETmis_ET.size(); ++i)
    {
      ETmis[i].ET	= ETmis_ET[i];
      ETmis[i].Phi	= ETmis_Phi[i];
      ETmis[i].Px	= ETmis_Px[i];
      ETmis[i].Py	= ETmis_Py[i];
    }

  Electron.resize(Electron_Charge.size());
  for(unsigned int i=0; i < Electron_Charge.size(); ++i)
    {
      Electron[i].Charge	= Electron_Charge[i];
      Electron[i].E	= Electron_E[i];
      Electron[i].EHoverEE	= Electron_EHoverEE[i];
      Electron[i].EtRatio	= Electron_EtRatio[i];
      Electron[i].Eta	= Electron_Eta[i];
      Electron[i].EtaCalo	= Electron_EtaCalo[i];
      Electron[i].IsolFlag	= Electron_IsolFlag[i];
      Electron[i].PT	= Electron_PT[i];
      Electron[i].Phi	= Electron_Phi[i];
      Electron[i].PhiCalo	= Electron_PhiCalo[i];
      Electron[i].Px	= Electron_Px[i];
      Electron[i].Py	= Electron_Py[i];
      Electron[i].Pz	= Electron_Pz[i];
      Electron[i].SumEt	= Electron_SumEt[i];
      Electron[i].SumPt	= Electron_SumPt[i];
    }

  Event.resize(Event_CouplingQCD.size());
  for(unsigned int i=0; i < Event_CouplingQCD.size(); ++i)
    {
      Event[i].CouplingQCD	= Event_CouplingQCD[i];
      Event[i].CouplingQED	= Event_CouplingQED[i];
      Event[i].Nparticles	= Event_Nparticles[i];
      Event[i].Number	= Event_Number[i];
      Event[i].ProcessID	= Event_ProcessID[i];
      Event[i].ScalePDF	= Event_ScalePDF[i];
      Event[i].Weight	= Event_Weight[i];
    }

  Jet.resize(Jet_Btag.size());
  for(unsigned int i=0; i < Jet_Btag.size(); ++i)
    {
      Jet[i].Btag	= Jet_Btag[i];
      Jet[i].E	= Jet_E[i];
      Jet[i].EHoverEE	= Jet_EHoverEE[i];
      Jet[i].Eta	= Jet_Eta[i];
      Jet[i].NCalo	= Jet_NCalo[i];
      Jet[i].NTracks	= Jet_NTracks[i];
      Jet[i].PT	= Jet_PT[i];
      Jet[i].Phi	= Jet_Phi[i];
      Jet[i].Px	= Jet_Px[i];
      Jet[i].Py	= Jet_Py[i];
      Jet[i].Pz	= Jet_Pz[i];
    }

  Muon.resize(Muon_Charge.size());
  for(unsigned int i=0; i < Muon_Charge.size(); ++i)
    {
      Muon[i].Charge	= Muon_Charge[i];
      Muon[i].E	= Muon_E[i];
      Muon[i].EHoverEE	= Muon_EHoverEE[i];
      Muon[i].EtRatio	= Muon_EtRatio[i];
      Muon[i].Eta	= Muon_Eta[i];
      Muon[i].EtaCalo	= Muon_EtaCalo[i];
      Muon[i].IsolFlag	= Muon_IsolFlag[i];
      Muon[i].PT	= Muon_PT[i];
      Muon[i].Phi	= Muon_Phi[i];
      Muon[i].PhiCalo	= Muon_PhiCalo[i];
      Muon[i].Px	= Muon_Px[i];
      Muon[i].Py	= Muon_Py[i];
      Muon[i].Pz	= Muon_Pz[i];
      Muon[i].SumEt	= Muon_SumEt[i];
      Muon[i].SumPt	= Muon_SumPt[i];
    }

  Particle.resize(Particle_Charge.size());
  for(unsigned int i=0; i < Particle_Charge.size(); ++i)
    {
      Particle[i].Charge	= Particle_Charge[i];
      Particle[i].D1	= Particle_D1[i];
      Particle[i].D2	= Particle_D2[i];
      Particle[i].E	= Particle_E[i];
      Particle[i].Eta	= Particle_Eta[i];
      Particle[i].M	= Particle_M[i];
      Particle[i].M1	= Particle_M1[i];
      Particle[i].M2	= Particle_M2[i];
      Particle[i].PID	= Particle_PID[i];
      Particle[i].PT	= Particle_PT[i];
      Particle[i].Phi	= Particle_Phi[i];
      Particle[i].Px	= Particle_Px[i];
      Particle[i].Py	= Particle_Py[i];
      Particle[i].Pz	= Particle_Pz[i];
      Particle[i].Status	= Particle_Status[i];
      Particle[i].T	= Particle_T[i];
      Particle[i].X	= Particle_X[i];
      Particle[i].Y	= Particle_Y[i];
      Particle[i].Z	= Particle_Z[i];
    }

  Photon.resize(Photon_E.size());
  for(unsigned int i=0; i < Photon_E.size(); ++i)
    {
      Photon[i].E	= Photon_E[i];
      Photon[i].EHoverEE	= Photon_EHoverEE[i];
      Photon[i].Eta	= Photon_Eta[i];
      Photon[i].PT	= Photon_PT[i];
      Photon[i].Phi	= Photon_Phi[i];
      Photon[i].Px	= Photon_Px[i];
      Photon[i].Py	= Photon_Py[i];
      Photon[i].Pz	= Photon_Pz[i];
    }

  TauJet.resize(TauJet_Charge.size());
  for(unsigned int i=0; i < TauJet_Charge.size(); ++i)
    {
      TauJet[i].Charge	= TauJet_Charge[i];
      TauJet[i].E	= TauJet_E[i];
      TauJet[i].EHoverEE	= TauJet_EHoverEE[i];
      TauJet[i].Eta	= TauJet_Eta[i];
      TauJet[i].NCalo	= TauJet_NCalo[i];
      TauJet[i].NTracks	= TauJet_NTracks[i];
      TauJet[i].PT	= TauJet_PT[i];
      TauJet[i].Phi	= TauJet_Phi[i];
      TauJet[i].Px	= TauJet_Px[i];
      TauJet[i].Py	= TauJet_Py[i];
      TauJet[i].Pz	= TauJet_Pz[i];
    }

  Tracks.resize(Tracks_Charge.size());
  for(unsigned int i=0; i < Tracks_Charge.size(); ++i)
    {
      Tracks[i].Charge	= Tracks_Charge[i];
      Tracks[i].E	= Tracks_E[i];
      Tracks[i].Eta	= Tracks_Eta[i];
      Tracks[i].EtaOuter	= Tracks_EtaOuter[i];
      Tracks[i].PT	= Tracks_PT[i];
      Tracks[i].Phi	= Tracks_Phi[i];
      Tracks[i].PhiOuter	= Tracks_PhiOuter[i];
      Tracks[i].Px	= Tracks_Px[i];
      Tracks[i].Py	= Tracks_Py[i];
      Tracks[i].Pz	= Tracks_Pz[i];
      Tracks[i].Vx	= Tracks_Vx[i];
      Tracks[i].Vy	= Tracks_Vy[i];
      Tracks[i].Vz	= Tracks_Vz[i];
    }

  ZDChits.resize(ZDChits_E.size());
  for(unsigned int i=0; i < ZDChits_E.size(); ++i)
    {
      ZDChits[i].E	= ZDChits_E[i];
      ZDChits[i].T	= ZDChits_T[i];
      ZDChits[i].genE	= ZDChits_genE[i];
      ZDChits[i].genEta	= ZDChits_genEta[i];
      ZDChits[i].genPT	= ZDChits_genPT[i];
      ZDChits[i].genPhi	= ZDChits_genPhi[i];
      ZDChits[i].genPx	= ZDChits_genPx[i];
      ZDChits[i].genPy	= ZDChits_genPy[i];
      ZDChits[i].genPz	= ZDChits_genPz[i];
      ZDChits[i].hadronic_hit	= ZDChits_hadronic_hit[i];
      ZDChits[i].pid	= ZDChits_pid[i];
      ZDChits[i].side	= ZDChits_side[i];
    }
}

using namespace std;

void getcombin(int n, int r, vector<vector<int> >& vcomb);

std::vector<Jet_s> makepseudojets(std::vector<Jet_s> sJet)
{
  double PT1MPT2min = 10000;
  std::vector<int> npj1min;
  std::vector<int> npj2min;

  //cout << sJet.size() << endl;
  for (unsigned int i=0; i<sJet.size(); i++) {
    int i1 = i+1;
    int i2 = sJet.size() - (i+1);
    if (i1>i2) continue;
    //cout << "i1, i2: " << i1 << " " << i2 << endl;
    vector<vector<int> > vcomb;
    getcombin(sJet.size(), i2, vcomb);
    for (unsigned int icomb=0; icomb<vcomb.size(); icomb++) {
      std::vector<int> npj1, npj2;
      //cout << "combination: " << icomb << endl;
      set<int> setcomb(vcomb[icomb].begin(), vcomb[icomb].end());
      //cout << "Complements: " << icomb+1 << "\t";
      double pj1PT = 0;
      double pj2PT = 0;
      for (unsigned int j=0; j<sJet.size(); j++) {
	if (setcomb.find(j)!=setcomb.end()) {
	  pj1PT += sJet[j].PT;
	  npj1.push_back(j);
	} else {
	  pj2PT += sJet[j].PT;
	  npj2.push_back(j);
	  //cout << j << "\t"; 
	}
      }
      double PT1MPT2 = fabs(pj1PT - pj2PT);
      //cout << "diff: " << PT1MPT2 << endl;
      //cout << PT1MPT2 << " " << PT1MPT2min << endl;
      if (PT1MPT2 <= PT1MPT2min) {
	PT1MPT2min = PT1MPT2;
	npj1min = npj1;
	npj2min = npj2; 
      }

      //cout << endl;
    }
  }

  /*
  cout << "ptdiff: " << PT1MPT2min << endl;
  cout << "Results:" << endl;
  cout << "ptdiff: " << PT1MPT2min << endl;
  cout << "pj1: " ;
  for (unsigned int i=0; i<npj1min.size(); i++) {
    cout << npj1min[i] << " ";
  }
  cout << endl;
  cout << "pj2: " ;
  for (unsigned int i=0; i<npj2min.size(); i++) {
    cout << npj2min[i] << " ";
  }
  cout << endl;
  
  for (unsigned int i=0; i<sJet.size(); i++) {
    cout << "PT " << i << " " << sJet[i].PT << endl;
  }
  */

  Jet_s pj1, pj2;
  pj1.Px = 0;
  pj1.Py = 0;
  pj1.PT = 0;
  for (unsigned int i=0; i<npj1min.size(); i++) {
    pj1.Px += sJet[npj1min[i]].Px;
    pj1.Py += sJet[npj1min[i]].Py;
    pj1.PT += sJet[npj1min[i]].PT;
    //cout << npj1min[i] << " " << pj1.PT << endl;
  }

  pj2.Px = 0;
  pj2.Py = 0;
  pj2.PT = 0;
  for (unsigned int i=0; i<npj2min.size(); i++) {
    pj2.Px += sJet[npj2min[i]].Px;
    pj2.Py += sJet[npj2min[i]].Py;
    pj2.PT += sJet[npj2min[i]].PT;
    //cout << npj2min[i] << " " << pj2.PT << endl;
  }
  
  std::vector<Jet_s> pJet;
  pJet.push_back(pj1);
  pJet.push_back(pj2);

  return pJet;
}


double MT(std::vector<Jet_s> sJet)
{
  double PT = 0;
  double Px = 0;
  double Py = 0;
  for (unsigned int i=0; i<sJet.size(); i++) {
    PT += sJet[i].PT;
    Px += sJet[i].Px;
    Py += sJet[i].Py;
  }
  return sqrt(PT*PT - Px*Px - Py*Py);
}

double alphaT(std::vector<Jet_s> pJet)
{
  if (pJet.size()!=2) return -1;
  double PT = pJet[0].PT + pJet[1].PT;
  double Px = pJet[0].Px + pJet[1].Px;
  double Py = pJet[0].Py + pJet[1].Py;

  double mT = sqrt(PT*PT - Px*Px - Py*Py);
  double aT;
  if (pJet[0].PT > pJet[1].PT) {
    aT = pJet[1].PT / mT;
  } else {
    aT = pJet[0].PT / mT;
  }

  return aT;
}

double deltaphistar(std::vector<Jet_s> sJet)
{
  double jitotpx = 0;
  double jitotpy = 0;
  double jitotpz = 0;
  double jitotE = 0;
  for (unsigned int i=0; i<sJet.size(); i++) {
    jitotpx -= sJet[i].Px;
    jitotpy -= sJet[i].Py;
    jitotpz -= sJet[i].Pz;
    jitotE -= sJet[i].E;
  }

  TLorentzVector lv_sum;
  TLorentzVector lv_jk;
  double deltaphistar = 100;
  for (unsigned int k=0; k<sJet.size(); k++) {
    lv_sum.SetPxPyPzE(jitotpx+sJet[k].Px, jitotpy+sJet[k].Py, jitotpz+sJet[k].Pz,
                      jitotE+sJet[k].E);
    lv_jk.SetPxPyPzE(sJet[k].Px, sJet[k].Py, sJet[k].Pz, sJet[k].E);
    double dphi = fabs(lv_sum.DeltaPhi(lv_jk));
    if (dphi < deltaphistar)
      deltaphistar = dphi;
  }

  return deltaphistar;

}



#endif
