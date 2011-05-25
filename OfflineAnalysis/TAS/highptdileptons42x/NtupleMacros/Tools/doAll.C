{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("Events"); 
  ch->Add("/nfs-4/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v2_AOD/V04-01-03/DoubleElectronTriggerSkim/skimmed_ntuple_163630_2.root");
  ScanChain(ch); 
}