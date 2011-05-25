{

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
// TCHEL - 1.7
// TCHEM - 3.3 *
// TCHET - 10.2
// SSVHEM - 1.74 *
// SSVHET - 3.05

  vector<TString> a;
  a.clear();

  // Analysis Dirs
  a.push_back("Common");
//  a.push_back("Top");
//  a.push_back("QCDLike");
  a.push_back("SameSign");

  // Compilation  
  TString dir = gSystem->pwd();
  TString run = gSystem->pwd()+string("/run");
  gSystem->Load("NtupleMacros/Tools/MiniFWLite/libMiniFWLite.so");

  int kk = a.size()-1;
  for ( size_t j=0; j<a.size(); j++ ){
    gSystem->AddIncludePath("-I$ROOTSYS/include -I+dir/+a[j]");
    void *thedir = gSystem->OpenDirectory(a[j]);
    const char *dirEntry;
    gSystem->ChangeDirectory(a[j]);
    gROOT->ProcessLine(".! rm *.so *.d ");
    if (j == kk) {
      cout << "Processing within the analysis dir " << a[j] << endl;
      gROOT->ProcessLine(".L ScanChain.C+");
    } else {
      TPMERegexp re(".C", "g");
      while ((dirEntry = gSystem->GetDirEntry(thedir))) {
	const TString s1(dirEntry);
	if(*dirEntry =='.') continue;
	if (re.Match(s1)) gSystem->CompileMacro(dirEntry, "++k");
      }
    }
    gSystem->ChangeDirectory(dir);
  }
  gSystem->ChangeDirectory(run);

  // Analysis Flags
  // Flags for files to run over
  bool    rundata        = false;
  bool    rundata2010    = false;
  bool    runttbar       = true;
  bool    runttotr       = false;
  bool    runWjets       = false;
  bool    runDYee        = false;
  bool    runDYmm        = false;
  bool    runDYtautau    = false;
  bool    runQCDPt15     = false; //2010
  bool    runQCDPt30     = false; //2010
  bool    runGamma15     = false; //2010
  bool    runWW          = false; //2010
  bool    runtW          = false; 
  bool    runLM0         = false; //2010
  bool    runLM1         = false; //2010
  bool    runVqq         = false; //2010
  bool    runWc          = false; //2010
  bool    runWZ          = false; //2010
  bool    runZZ          = false; //2010
  bool    runVV          = false;
  bool    runVgamma      = false; //2010
  bool    runWgamma      = false; //2010


  float   jetTriggerPt   = 40;    
  
  //NLO cross-sections
//  float   kttdil    = 157.5;
//  float   kWjets    = 31314.;
  float   kttdil    = 1.;
  float   kWjets    = 1.;
  float   kDYee     = 1.;
  float   kDYmm     = 1.;
  float   kDYtautau = 1.;
  float   kqcd15    = 1.;
  float   minbias   = 1.;

  bool    doFRestimation = false;    


vector<TString> v_Cuts;
// v_Cuts.push_back("useOSleptons");         // OS leptons
v_Cuts.push_back("useSSleptons");         // SS leptons
// v_Cuts.push_back("usePtGt2020");         // use leptons with pt > 20
v_Cuts.push_back("usePtGt2010");         // one lepton > 20, other > 10
// v_Cuts.push_back("usePtGt1010");           // both leptons > 10
// v_Cuts.push_back("excludePtGt2020");     // one lepton > 10, < 20 other >20
// v_Cuts.push_back("used0corrPV");           // use the d0 corrected for the hihest PV 
v_Cuts.push_back("applylepIDCuts");      // apply tight ID cuts
v_Cuts.push_back("applylepIsoCuts");     // tight iso cuts 
// v_Cuts.push_back("applyAlignmentCorrection"); // apply alignment corrections
// v_Cuts.push_back("removedEtaCutInEndcap"); // apply alignment corrections
// v_Cuts.push_back("applyFOv1Cuts");
// v_Cuts.push_back("applyFOv2Cuts");
// v_Cuts.push_back("applyFOv3Cuts");
v_Cuts.push_back("applyTriggers");       // apply triggers
//v_Cuts.push_back("vetoZmass");           // no leptons in zmass
// v_Cuts.push_back("requireZmass");        // leptons only in zmass
//v_Cuts.push_back("hypDisamb");           // do hyp. disambiguation
//v_Cuts.push_back("useCorMET");           // use corrected calo MET ---> NOT SUPPORTED RIGHT NOW
// v_Cuts.push_back("usetcMET");            // use tcMET
v_Cuts.push_back("usepfMET");   //use PFMET 
// v_Cuts.push_back("chargeFlip");   //use chargeFlip
v_Cuts.push_back("vetoMET");             // cut on MET  
//v_Cuts.push_back("vetoProjectedMET");    // cut on projected MET
// v_Cuts.push_back("usecaloJets");         // use caloJETs for jet counting
// v_Cuts.push_back("usejptJets");          // use jpt jets for jet counting
v_Cuts.push_back("usepfJets");  // use pf jets for jet counting
v_Cuts.push_back("vetoJets");
v_Cuts.push_back("requireEcalEls");
//
// v_Cuts.push_back("requireBTag");
// v_Cuts.push_back("requireTCHEL");
// v_Cuts.push_back("requiredouble");

// v_Cuts.push_back("useFlipRateEstimation");
// v_Cuts.push_back("estimateSingleFakes");
// v_Cuts.push_back("estimateDoubleFakes");

  TChain  *ch_data    = new TChain("Events");  
  TChain  *ch_ttbar   = new TChain("Events");
  TChain  *ch_ttbaro   = new TChain("Events");
  TChain  *ch_wjets   = new TChain("Events");
  TChain  *ch_vgamma   = new TChain("Events");
  TChain  *ch_wgamma   = new TChain("Events");
  TChain  *ch_dyee    = new TChain("Events");
  TChain  *ch_dymm    = new TChain("Events");
  TChain  *ch_dytt    = new TChain("Events");
  TChain  *ch_qcdpt15 = new TChain("Events");
  TChain  *ch_qcdpt30 = new TChain("Events");
  TChain  *ch_gamma15 = new TChain("Events");
  TChain  *ch_ww      = new TChain("Events");
  TChain  *ch_tw      = new TChain("Events");
  TChain  *ch_lm0      = new TChain("Events");
  TChain  *ch_lm1      = new TChain("Events");
  TChain  *ch_vqq      = new TChain("Events");
  TChain  *ch_wc       = new TChain("Events");
  TChain  *ch_wz      = new TChain("Events");
  TChain  *ch_zz      = new TChain("Events");
  TChain  *ch_vv      = new TChain("Events");


// Luminosity
//  const float LUMINORM = 0.00006511; // 78 nb-1
//   const float LUMINORM = 0.00030355; // 303.55 nb-1
//   const float LUMINORM = 0.00083828; // 838.28 nb-1
//   const float LUMINORM = 0.0011; // 1.1 pb-1
//   const float LUMINORM = 0.001; // 1 pb-1
//     const float LUMINORM = 0.01; // 10 pb-1
//   const float LUMINORM = 0.00288; // 2.88 pb-1
//   const float LUMINORM = 0.0361; // 36.1 pb-1
//   const float LUMINORM = 0.04341; // 43.41 pb-1
//   const float LUMINORM = 0.153; // 43.41 pb-1
   const float LUMINORM = 0.1911; // 191.1 pb-1

//     const float LUMINORM = 1; // 1 fb-1

  if(rundata2010) {
    cout << "Doing 2010 data" << endl;

      ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple*.root");
      ch_data->Add("/nfs-3/userdata/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed_ntuple*.root");
      ch_data->Add("/nfs-3/userdata/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple*.root");

      ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple*.root");
      ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14-00/diLepPt1020Skim/skimmed_ntuple*.root");
      ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple*.root");

    ScanChain(ch_data, v_Cuts, "data", doFRestimation, jetTriggerPt);
    hist::color("data", kRed);
  }

  if(rundata) {
    cout << "Doing 2011 data" << endl;

//DoubleElectron re-reco
      ch_data->Add("/nfs-4/userdata/cms2/DoubleElectron_Run2011A-Apr22ReReco-v2_AOD/V04-01-05/DoubleElectronTriggerSkim/*.root");
//v1 dataset
      ch_data->Add("/nfs-4/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v1_AOD/V04-00-13/DoubleElectronTriggerSkim_merged/*.root");
      ch_data->Add("/nfs-4/userdata/cms2/DoubleMu_Run2011A-PromptReco-v1_AOD/V04-00-13/DoubleMuTriggerSkim_merged/*.root");
      ch_data->Add("/nfs-4/userdata/cms2/MuEG_Run2011A-PromptReco-v1_AOD/V04-00-13/*.root");
      ch_data->Add("/nfs-4/userdata/cms2/SingleMu_Run2011A-PromptReco-v1_AOD/V04-00-13/*.root");

      ch_data->Add("/nfs-4/userdata/cms2/DoubleElectron_Run2011A-PromptReco-v2_AOD/V04-01-03/DoubleElectronTriggerSkim/*.root");
      ch_data->Add("/nfs-4/userdata/cms2/DoubleMu_Run2011A-PromptReco-v2_AOD/V04-01-03/DoubleMuTriggerSkim/*.root");
      ch_data->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-01-03/MuEG_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/*.root");
      ch_data->Add("/hadoop/cms/store/user/jaehyeok/CMSSW_4_1_2_patch1_V04-01-03/SingleMu_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03 /*.root");

    ScanChain(ch_data, v_Cuts, "data", doFRestimation, jetTriggerPt);
    hist::color("data", kRed);
  }


  if(runttbar) {
    cout << "Doing the ttbar sample" << endl;
    ch_ttbar->Add("/nfs-4/userdata/cms2/TT_TuneZ2_7TeV-pythia6-tauola_Summer11-PU_S3_START42_V11-v2/V04-02-09/merged_ntuple*root"); 
    ScanChain(ch_ttbar, v_Cuts, "ttbar",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ttbar", kRed);
  }
  
  if(runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    ch_ttbaro->Add("/nfs-3/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*root");
    ScanChain(ch_ttbaro, v_Cuts,"ttotr", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ttotr", 30);
  }

  if(runVV) { // 
    cout << "Processing VV" << endl;
    ch_vv->Add("/nfs-3/userdata/cms2/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ch_vv->Add("/nfs-3/userdata/cms2/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ch_vv->Add("/nfs-3/userdata/cms2/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ScanChain(ch_vv, v_Cuts, "vv", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("vv", kRed-5);
  } 

  if(runtW) {
    cout << "Processing tW" << endl;
    ch_tw->Add("/nfs-3/userdata/cms2/TToBLNu_TuneZ2_t-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ch_tw->Add("/nfs-3/userdata/cms2/TToBLNu_TuneZ2_s-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ch_tw->Add("/nfs-3/userdata/cms2/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ScanChain(ch_tw, v_Cuts, "tw", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("tw", kRed-3);
  } 

  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
      ch_wjets->Add("/nfs-3/userdata/cms2/WToENu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
      ch_wjets->Add("/nfs-3/userdata/cms2/WToMuNu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
      ch_wjets->Add("/nfs-3/userdata/cms2/WToMuNu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");

    ScanChain(ch_wjets,v_Cuts, "wj", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wjets", 40);
  }

  if (runDYee) {
    cout << "Processing DY->ee" << endl;
//        ch_dyee->Add("/nfs-3/userdata/cms2/DYToEE_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
//        ch_dyee->Add("/nfs-3/userdata/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
        ch_dyee->Add("/nfs-3/userdata/cms2/DYToEE_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/dilep_ZMassLessThan50Skim/*.root");
        ch_dyee->Add("/nfs-3/userdata/cms2/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ScanChain(ch_dyee, v_Cuts, "DYee",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
//      ch_dymm->Add("/nfs-3/userdata/cms2/DYToMuMu_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
//      ch_dymm->Add("/nfs-3/userdata/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
      ch_dymm->Add("/nfs-3/userdata/cms2/DYToMuMu_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/dilep_ZMassLessThan50Skim/*.root");
      ch_dymm->Add("/nfs-3/userdata/cms2/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ScanChain(ch_dymm, v_Cuts, "DYmm",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYmm", kCyan);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
//      ch_dytt->Add("/nfs-3/userdata/cms2/DYToTauTau_M-10To20_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
//      ch_dytt->Add("/nfs-3/userdata/cms2/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
      ch_dytt->Add("/nfs-3/userdata/cms2/DYToTauTau_M-20_CT10_TuneZ2_7TeV-powheg-pythia-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/dilep_ZMassLessThan50Skim/*.root");
      ch_dytt->Add("/nfs-3/userdata/cms2/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/merged*.root");
    ScanChain(ch_dytt,v_Cuts, "DYtautau",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYtautau", kBlack);
  }


  if (runVgamma) { // 2010
    cout << "Processing Vgamma.."<<endl;
    ch_vgamma->Add("/nfs-3/userdata/cms2/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
    ScanChain(ch_vgamma,v_Cuts, "vgamma", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("vgamma", 40);
  } 

  if (runWgamma) { // 2010
    cout << "Processing Wgamma.."<<endl;
    ch_vgamma->Add("/nfs-3/userdata/cms2/Wgamma_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_wgamma,v_Cuts, "wgamma", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wgamma", 40);
  }

  if(runQCDPt15) { // 2010
    cout << "Processing QCDPt15" << endl;
    ch_qcdpt15->Add("/nfs-3/userdata/cms2/EarlyDataSamples/QCD_Pt15_Spring10-START3X_V26_S09-v1_Single5GeV/V03-04-08-01/*.root");
    ScanChain(ch_qcdpt15, v_Cuts, "qcd15", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("qcd15", 28);
  }

  if(runQCDPt30) { //2010
    cout << "Processing QCDPt30" << endl;
    ch_qcdpt30->Add("/nfs-3/userdata/cms2/EarlyDataSamples/QCD_Pt30_Spring10-START3X_V26_S09-v1_Single5GeV/V03-04-08-01/merged_ntuple*.root"); 
    ScanChain(ch_qcdpt30, v_Cuts, "qcd30", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("qcd30", 27);
  }

  if(runGamma15) { //2010
    cout << "Processing runGamma15" << endl;
    ch_gamma15->Add("/nfs-3/userdata/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*.root");
    ScanChain(ch_gamma15, v_Cuts, "ph15", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ph15", 26);
  } 

  if(runWW) { //2010
    cout << "Processing WW" << endl;
    ch_ww->Add("/nfs-3/userdata/cms2/WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/*.root");
    ScanChain(ch_ww, v_Cuts, "WW" ,doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ww", 31);
  }

  if(runLM0) { // 2010
    cout << "Processing LM0" << endl;
    ch_lm0->Add("/nfs-3/userdata/cms2/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/*.root");
    ScanChain(ch_lm0, v_Cuts, "LM0x", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("LM0x", kRed-3);
  }

  if(runLM1) { // 2010
    cout << "Processing LM1" << endl;
    ch_lm1->Add("/nfs-3/userdata/cms2/LM1_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/*.root");
    ScanChain(ch_lm1, v_Cuts, "LM1x", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("LM1x", kRed-3);
  }


  if(runVqq) { //2010
    cout << "Processing Vqq" << endl;
    ch_vqq->Add("/nfs-3/userdata/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_vqq, v_Cuts, "Vqq", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("Vqq", kRed-3);
  }

  if(runWc) {// 2010
    cout << "Processing Wc" << endl;
    ch_wc->Add("/nfs-3/userdata/cms2/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/*.root");
    ScanChain(ch_wc, v_Cuts, "wc", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wc", kRed-3);
  }

  if(runWZ) {//2010
    cout << "Processing WZ" << endl;
    ch_wz->Add("/nfs-3/userdata/cms2/WZ_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_wz, v_Cuts, "wz", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wz", kRed-4);
  }

  if(runZZ) {//2010
    cout << "Processing ZZ" << endl;
    ch_zz->Add("/nfs-3/userdata/cms2/ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/*.root");
    ScanChain(ch_zz, v_Cuts, "zz", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("zz", kRed-5);
  }


  TString cutstring = "";
  for(unsigned int i = 0; i < v_Cuts.size(); i++) 
    cutstring       = cutstring + "_" + v_Cuts.at(i);
  
  if(doFRestimation)
    cutstring = "FRhist" + cutstring + ".root";
  else cutstring   = "hist" + cutstring + ".root";
  
  cout << "Saving histograms to: " << cutstring << endl;
  hist::saveHist(cutstring.Data());
  hist::deleteHistos();
  hist::loadHist(cutstring.Data(),0,"*_hnJet_*");
  hist::loadHist(cutstring.Data(),0,"*_hnbJet_*");

  std::vector<TString> v_prfxsToCombine;
  v_prfxsToCombine.push_back("DYee");
  v_prfxsToCombine.push_back("DYmm");
  hist::combineHists(v_prfxsToCombine, "DYeemm");



  //hist::loadHist(cutstring.Data());
  //print table
  //browseStacks( true, false , "qcd15", cutstring); //need to supply a reference channel. Should not be "data" because "data" should not be in the stacks
  if (doFRestimation) {
 //   printNJets(false, "%6.2f","ttdil", true,true,true,false, true);
    printNJets(false, "%6.2f","ttdil", false,true,true,false, true);
  } else {
//   printNJets(false, "%6.2f","ttdil",false,true,true);
   printNJets(false, "%6.2f","ttdil", true,false,true, false, false); 
  }
}
