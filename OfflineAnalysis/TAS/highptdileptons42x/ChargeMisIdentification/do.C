{

  using namespace std;

  // Load various tools  
  gROOT->ProcessLine(".x setup.C(true)");
  gSystem->Load("../NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");
  gROOT->LoadMacro("ScanChain.C+");

  TChain *chain = new TChain("Events");
  //chain->Add("/nfs-3/userdata/cms2/SingleElectronPt7to100_CMSSW_3_8_4_patch3/V03-06-14/merged_ntuple_65.root");
//  chain->Add("/nfs-3/userdata/cms2/SingleElectronPt7to100_CMSSW_3_8_4_patch3/V03-06-14/merged_ntuple*.root");
//  chain->Add("/nfs-4/userdata/spadhi/TAS/tas_singleelectronAprl30/*.root");
//  chain->Add("/nfs-4/userdata/spadhi/TAS/tas_singleelectronMay05/*.root");
  chain->Add("/nfs-4/userdata/spadhi/TAS/tas_singleelectronMay13/*.root");

  ScanChain(chain);

  //save all the histograms

 cout << "all done" << endl;
    
  const char* outFile = "myHist.root";
  hist::saveHist(outFile);
  hist::deleteHistos();
}
