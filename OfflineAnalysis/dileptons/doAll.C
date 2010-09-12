{
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
  bool    runttbar       = true;
  bool    runttotr       = false;
  bool    runWjets       = true;
  bool    runDYee        = true;
  bool    runDYmm        = true;
  bool    runDYtautau    = true;
  bool    runQCDPt15     = false;
  bool    runQCDPt30     = false;
  bool    runGamma15     = false;
  bool    runWW          = true;
  bool    runtW          = true;
  bool    runLM0         = false;
  bool    runVqq         = false;
  bool    runWc          = false;
  bool    runWZ          = true;
  bool    runZZ          = true;
  bool    runVV          = false;
  bool    runVgamma      = true;
  bool    runWgamma      = false;


  float   jetTriggerPt   = 15;    
  
  //NLO cross-sections
  float   kttdil    = 157.5;
  float   kWjets    = 31314.;
  float   kDYee     = 1.;
  float   kDYmm     = 1.;
  float   kDYtautau = 1.;
  float   kqcd15    = 1.;
  float   minbias   = 1.;

  bool    doFRestimation = true;    


vector<TString> v_Cuts;
// v_Cuts.push_back("useOSleptons");         // OS leptons
v_Cuts.push_back("useSSleptons");         // SS leptons
// v_Cuts.push_back("usePtGt2020");         // use leptons with pt > 20
v_Cuts.push_back("usePtGt2010");         // one lepton > 20, other > 10
// v_Cuts.push_back("usePtGt1010");           // both leptons > 10
// v_Cuts.push_back("excludePtGt2020");     // one lepton > 10, < 20 other >20
// v_Cuts.push_back("used0corrPV");           // use the d0 corrected for the hihest PV 
// v_Cuts.push_back("applylepIDCuts");      // apply tight ID cuts
// v_Cuts.push_back("applylepIsoCuts");     // tight iso cuts 
v_Cuts.push_back("applyAlignmentCorrection"); // apply alignment corrections
v_Cuts.push_back("removedEtaCutInEndcap"); // apply alignment corrections
v_Cuts.push_back("applyFOv1Cuts");
// v_Cuts.push_back("applyFOv2Cuts");
// v_Cuts.push_back("applyFOv3Cuts");
v_Cuts.push_back("applyTriggers");       // apply triggers
// v_Cuts.push_back("vetoZmass");           // no leptons in zmass
// v_Cuts.push_back("requireZmass");        // leptons only in zmass
//v_Cuts.push_back("hypDisamb");           // do hyp. disambiguation
//v_Cuts.push_back("useCorMET");           // use corrected calo MET ---> NOT SUPPORTED RIGHT NOW
v_Cuts.push_back("usetcMET");            // use tcMET
// v_Cuts.push_back("usepfMET");   //use PFMET 
// v_Cuts.push_back("chargeFlip");   //use chargeFlip
// v_Cuts.push_back("vetoMET");             // cut on MET  
//v_Cuts.push_back("vetoProjectedMET");    // cut on projected MET
// v_Cuts.push_back("usecaloJets");         // use caloJETs for jet counting
//v_Cuts.push_back("usejptJets");          // use jpt jets for jet counting
v_Cuts.push_back("usepfJets");  // use pf jets for jet counting
//v_Cuts.push_back("vetoJets");
v_Cuts.push_back("requireEcalEls");
// v_Cuts.push_back("useFlipRateEstimation");
v_Cuts.push_back("estimateSingleFakes");
// v_Cuts.push_back("estimateDoubleFakes");

  TChain  *ch_data    = new TChain("Events");  
  TChain  *ch_ttbar   = new TChain("Events");
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
//     const float LUMINORM = 0.01; // 10 pb-1
   const float LUMINORM = 0.00288; // 2.88 pb-1

//     const float LUMINORM = 1; // 1 fb-1

  if(rundata) {
    cout << "Doing data" << endl;

 //   ch_data->Add("/hadoop/cms/store/user/spadhi/CMS2_V03-04-26-02/Commissioning10-SD_EG-Jun14thSkim_v1/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/MinimumBias_Commissioning10-SD_EG-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-Jul16thReReco-v2_RECO/V03-04-26-07/singleLepPt5Skim//*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-07/singleLepPt10Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/*.root");

//    ch_data->Add("/hadoop/cms/store/user/spadhi/CMS2_V03-04-26-02/Commissioning10-SD_Mu-Jun14thSkim_v1/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/MinimumBias_Commissioning10-SD_Mu-Jun14thSkim_v1_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-Jun14thReReco_v1_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-25/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-01/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-Jul16thReReco-v1_RECO/V03-04-26-07/singleLepPt5Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-07/singleLepPt10Skim/*.root");
    ch_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/*.root");

    ScanChain(ch_data, v_Cuts, "data", doFRestimation, jetTriggerPt);
    hist::color("data", kRed);
  }

  if(runttbar) {
    cout << "Doing the ttbar sample" << endl;
    // ch_ttbar->Add("/nfs-3/userdata/cms2/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08/*.root"); 
    ch_ttbar->Add("/nfs-3/userdata/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/*.root"); 
    ScanChain(ch_ttbar, v_Cuts, "ttbar",doFRestimation, jetTriggerPt, LUMINORM, kttdil);
    hist::color("ttbar", kYellow);
  }
  
  if(runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
  //  ch_ttbar->Add("/nfs-3/userdata/cms2/TTbar_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ch_ttbar->Add("/nfs-3/userdata/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/*.root");
    ScanChain(ch_ttbar, v_Cuts,"ttotr", doFRestimation, jetTriggerPt, LUMINORM, kttotr);
    hist::color("ttotr", 30);
  }

  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
//    ch_wjets->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Wenu_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
//    ch_wjets->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Wmunu_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
//    ch_wjets->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Wtaunu_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
    ch_wjets->Add("/nfs-3/userdata/cms2/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/*.root");
    ScanChain(ch_wjets,v_Cuts, "wj", doFRestimation, jetTriggerPt, LUMINORM, kWjets);
    hist::color("wjets", 40);
  }

  if (runVgamma) {
    cout << "Processing Vgamma.."<<endl;
    ch_vgamma->Add("/nfs-3/userdata/cms2/PhotonVJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
    ScanChain(ch_vgamma,v_Cuts, "vgamma", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("vgamma", 40);
  } 

  if (runWgamma) {
    cout << "Processing Wgamma.."<<endl;
    ch_vgamma->Add("/nfs-3/userdata/cms2/Wgamma_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_wgamma,v_Cuts, "wgamma", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wgamma", 40);
  }

  if(runQCDPt15) {
    cout << "Processing QCDPt15" << endl;
    ch_qcdpt15->Add("/nfs-3/userdata/cms2/EarlyDataSamples/QCD_Pt15_Spring10-START3X_V26_S09-v1_Single5GeV/V03-04-08-01/*.root");
    ScanChain(ch_qcdpt15, v_Cuts, "qcd15", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("qcd15", 28);
  }

  if(runQCDPt30) {
    cout << "Processing QCDPt30" << endl;
    ch_qcdpt30->Add("/nfs-3/userdata/cms2/EarlyDataSamples/QCD_Pt30_Spring10-START3X_V26_S09-v1_Single5GeV/V03-04-08-01/merged_ntuple*.root"); 
    ScanChain(ch_qcdpt30, v_Cuts, "qcd30", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("qcd30", 27);
  }

  if(runGamma15) {
    cout << "Processing runGamma15" << endl;
    ch_gamma15->Add("/nfs-3/userdata/cms2/PhotonJet_Pt15_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*.root");
    ScanChain(ch_gamma15, v_Cuts, "ph15", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ph15", 26);
  } 

  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    ch_dyee->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Zee_Spring10-START3X_V26_S09-v1/*.root"                          );
    ch_dyee->Add("/nfs-3/userdata/cms2/EarlyDataSamples/DYee_M10to20_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root"    );
    ch_dyee->Add("/nfs-3/userdata/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root"         );
    ScanChain(ch_dyee, v_Cuts, "DYee",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    ch_dymm->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Zmumu_Spring10-START3X_V26_S09-v1/*.root"                        );
    ch_dymm->Add("/nfs-3/userdata/cms2/EarlyDataSamples/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/*.root"               );
    ch_dymm->Add("/nfs-3/userdata/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root"         );
    ScanChain(ch_dymm, v_Cuts, "DYmm",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYmm", kCyan);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    ch_dytt->Add("/nfs-3/userdata/cms2/EarlyDataSamples/Ztautau_Spring10-START3X_V26_S09-v1/*.root"                      );
    ch_dytt->Add("/nfs-3/userdata/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root"         );
    ScanChain(ch_dytt,v_Cuts, "DYtautau",doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("DYtautau", kBlack);
  }

  if(runWW) {
    cout << "Processing WW" << endl;
    ch_ww->Add("/nfs-3/userdata/cms2/WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/*.root");
    ScanChain(ch_ww, v_Cuts, "WW" ,doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("ww", 31);
  }

  if(runtW) {
    cout << "Processing tW" << endl;
    ch_tw->Add("/nfs-3/userdata/cms2/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/*.root");
    ScanChain(ch_tw, v_Cuts, "tw", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("tw", kRed-3);
  }

  if(runLM0) {
    cout << "Processing LM0" << endl;
    ch_lm0->Add("/nfs-3/userdata/cms2/LM0_Spring10-START3X_V26_S09-v1/V03-04-13-01/*.root");
    ScanChain(ch_lm0, v_Cuts, "LM0x", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("LM0x", kRed-3);
  }

  if(runVqq) {
    cout << "Processing Vqq" << endl;
    ch_vqq->Add("/nfs-3/userdata/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_vqq, v_Cuts, "Vqq", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("Vqq", kRed-3);
  }

  if(runWc) {
    cout << "Processing Wc" << endl;
    ch_wc->Add("/nfs-3/userdata/cms2/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/*.root");
    ScanChain(ch_wc, v_Cuts, "wc", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wc", kRed-3);
  }

  if(runWZ) {
    cout << "Processing WZ" << endl;
    ch_wz->Add("/nfs-3/userdata/cms2/WZ_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");
    ScanChain(ch_wz, v_Cuts, "wz", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("wz", kRed-4);
  }

  if(runZZ) {
    cout << "Processing ZZ" << endl;
    ch_zz->Add("/nfs-3/userdata/cms2/ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/*.root");
    ScanChain(ch_zz, v_Cuts, "zz", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("zz", kRed-5);
  }

  if(runVV) {
    cout << "Processing VV" << endl;
    ch_zz->Add("/nfs-3/userdata/cms2/VVJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-01/*.root");
    ScanChain(ch_vv, v_Cuts, "vv", doFRestimation, jetTriggerPt, LUMINORM);
    hist::color("vv", kRed-5);
  }


  TString cutstring = "";
  for(unsigned int i = 0; i < v_Cuts.size(); i++) 
    cutstring       = cutstring + "_" + v_Cuts.at(i);
  
  if(doFRestimation)
    cutstring = "FRhist" + cutstring + ".root";
  else cutstring   = "hist" + cutstring + ".root";
  
  cout << "Saving histograms to: " << cutstring << endl;
  hist::saveHist(cutstring.Data());
  //hist::delteHistos();
  //hist::loadHist(cutstring.Data());
  //print table
  //browseStacks( true, false , "qcd15", cutstring); //need to supply a reference channel. Should not be "data" because "data" should not be in the stacks
  if (doFRestimation) {
 //   printNJets(false, "%6.2f","ttdil", true,true,true,false, true);
    printNJets(false, "%6.2f","ttdil", false,true,true,false, true);
  } else {
   printNJets(false, "%6.2f","ttdil",false,true,true);
  }
}
