{
  gSystem->Load("../NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->ProcessLine(".L ScanChain.C++");
  TChain *ch = new TChain("Events"); 
//  ch->Add("/nfs-3/userdata/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");
//  ch->Add("/data/tmp/spadhi/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/*.root");
//  ch->Add("/nfs-3/userdata/cms2/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/*.root");
//  ch->Add("/nfs-3/userdata/cms2/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-START38_V12-v2/V03-06-17/*.root");
//  ch->Add("/data/tmp/spadhi/LM0_Spring10-START3X_V26_S09-v1/*.root");
  ch->Add("/nfs-3/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
//  ch->Add("/nfs-3/userdata/cms2/LM1_Spring10-START3X_V26_S09-v1/V03-04-13-01/*.root");
//  ch->Add("/nfs-3/userdata/cms2/LM0_Spring10-START3X_V26_S09-v1/V03-04-13-01/*.root");
  ScanChain(ch); 

  cout << "All done " << endl;
}
