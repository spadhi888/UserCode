{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("Events"); 
  ch->Add("/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML01_7TeV_v0/V03-04-13-01-gmsb/merged_ntuple.root");
  ScanChain(ch); 
}