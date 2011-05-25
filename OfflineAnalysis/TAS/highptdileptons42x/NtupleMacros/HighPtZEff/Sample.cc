#include "TChain.h"
#include "Sample.h"
#include "tools.h"
#include "../CORE/selections.h"
#include "../CORE/CMS2.h"
#include <cstdlib>
#include <string>
#include <iostream>

bool filterByProcess (enum Process sample)
{
     switch (sample) {
     case DYee: 
          return isDYee();
     case DYmm:
          return isDYmm();
     case DYtt:
          return isDYtt();
     default:
	  return true;
     }
     return true;
}

static const std::string prefix = (getenv("CMS2_NTUPLE_LOCATION") != 0) ? 
     std::string(getenv("CMS2_NTUPLE_LOCATION")) + "/" : "/data/tmp/";
     
//WW file
Sample fWW ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/WW_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true, 0. };
     return ret;
}

//WW file
Sample fWW_excl ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WW, kRed, 1, "ww", true, 0. };
     return ret;
}

//WZ file
Sample fWZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/WZ_incl_Summer08_IDEAL_V9_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, WZ, kBlue, 1, "wz", true, 0. };
     return ret;
}

//ZZ file
Sample fZZ ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/ZZ_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ZZ, kGreen, 1, "zz", true, 0. };
     return ret;
}

//Wjets file
Sample fWjets ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wjets", true, 0. };
     return ret;
}

//Wjets file
Sample fWc ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Wc-madgraph_Fall08_IDEAL_V9_reco-v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wjets, 40, 1, "wc", true };
     return ret;
}

// "vlqq" sample
Sample fVlqq ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/VQQ-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DY, kBlack, 1, "vlqq", true };
     return ret;
}

//DYee file
Sample fDYee ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1, "dyee", true, 0. };
     return ret;
}

//DYmm file
Sample fDYmm ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1, "dymm", true, 0. };
     return ret;
}

//DYtt file
Sample fDYtt ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1, "dytt", true, 0. };
     return ret;
}


// Pythia DY with no filtering when making ntuples
// These samples also contain pdf_info and genlepdaughters
//
// DYee
Sample fDYee_nofilter    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-06-NoFilter/Zee_M20/merged_ntuple*.root";
     //std::string sample = "/data/fio-1/tmp/wandrews/cms2-V01-02-06-NoFilter/Zee_M20/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1, "dyee_nofilter", true, 0. };
     return ret;
}
// DYmm
Sample fDYmm_nofilter    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-06-NoFilter/Zmumu_M20/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1, "dymm_nofilter", true, 0. };
     return ret;
}
// DYtt
Sample fDYtt_nofilter    ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-06-NoFilter/Ztautau_M20/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1, "dytt_nofilter", true, 0. };
     return ret;
}

//High Pt Z samples -- electrons
Sample fZeeJet80to120_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt80to120_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet80to120_nofilter", true, 0. };
     return ret;
}

Sample fZeeJet120to170_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt120to170_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet120to170_nofilter", true, 0. };
     return ret;
}

Sample fZeeJet170to230_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt170to230_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet170to230_nofilter", true, 0. };
     return ret;
}

Sample fZeeJet230to300_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt230to300_Summer08_IDEAL_V9_v3/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet230to300_nofilter", true, 0. };
     return ret;
}

Sample fZeeJet300toInf_nofilter()  
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt300toInf_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet300toInf_nofilter", true, 0. };
     return ret;
}

Sample fZeeJet300toInf_nofilter_sngl() //just the file used for new branches, but in original form
{
     TChain *c = new TChain("Events");
     std::string sample = "/home/users/wandrews/second/CMSSW_2_2_10/src/CMS2/NtupleMaker/ntuple_sngl.root"; //CAREFUL HERE!!!!
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJet300toInf_nofilter_sngl", true, 0. };
     return ret;
}

//Sample fZeeJet300toInf_nofilter_eleiso() //new branch(es) for ele iso
//{
//     TChain *c = new TChain("Events");
//     std::string sample = "/home/users/wandrews/second/CMSSW_2_2_10/src/CMS2/NtupleMaker/ntuple_eleiso.root";
//     c->Add(sample.c_str());
//     Sample ret = { c, DYee, kBlack, 1, "ZeeJet300toInf_nofilter_eleiso", true, 0. };
//     return ret;
//}

//All the high pt samples together
Sample fZeeJetALL80toInf_nofilter()  
{
     TChain *c = new TChain("Events");
	 std::string sample;
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt80to120_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt120to170_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt170to230_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());	 
	 sample= prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt230to300_Summer08_IDEAL_V9_v3/merged_ntuple*.root";
     c->Add(sample.c_str());
     sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZeeJet_Pt300toInf_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kBlack, 1, "ZeeJetALL80toInf_nofilter", true, 0. };
     return ret;
}

//High Pt Z samples -- muons
Sample fZmmJet80to120_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt80to120_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet80to120_nofilter", true, 0. };
     return ret;
}

Sample fZmmJet120to170_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt120to170_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet120to170_nofilter", true, 0. };
     return ret;
}

Sample fZmmJet170to230_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt170to230_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet170to230_nofilter", true, 0. };
     return ret;
}

Sample fZmmJet230to300_nofilter()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt230to300_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet230to300_nofilter", true, 0. };
     return ret;
}

Sample fZmmJet300toInf_nofilter()  
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt300toInf_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet300toInf_nofilter", true, 0. };
     return ret;
}

//Sample fZmmJet300toInf_nofilter_sngl() //just the file used for new branches, but in original form
//{
//     TChain *c = new TChain("Events");
//     std::string sample = "/home/users/wandrews/second/CMSSW_2_2_10/src/CMS2/NtupleMaker/ntuple_sngl.root"; //CAREFUL HERE!!!!
//     c->Add(sample.c_str());
//     Sample ret = { c, DYmm, kBlack, 1, "ZmmJet300toInf_nofilter_sngl", true, 0. };
//     return ret;
//}

//All the high pt samples together
Sample fZmmJetALL80toInf_nofilter()  
{
     TChain *c = new TChain("Events");
	 std::string sample;
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt80to120_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt120to170_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt170to230_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());	 
	 sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt230to300_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     sample = prefix + "dlevans/cms2-V01-02-11-NoFilter/ZmumuJet_Pt300toInf_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kBlack, 1, "ZmmJetALL80toInf_nofilter", true, 0. };
     return ret;
}


// low-mass DY sample
Sample fAstar ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/AstarJets-madgraph_Fall08_IDEAL_V9_v2/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DY, kBlack, 1, "astar", true };
     return ret;
}

// low-mass DY --> tau tau sample
Sample fDY20tt ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Ztautau_M20_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYtt, kBlack, 1, "dy20tt", true };
     return ret;
}

Sample fDY20mm ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Zmumu_M20_Summer08_IDEAL_V9_reco-v2/merged_ntuple**.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYmm, kCyan, 1, "dy20mm", true };
     return ret;
}

Sample fDY20ee ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Zee_M20_Summer08_IDEAL_V9_reco-v3/merged_ntuple**.root";
     c->Add(sample.c_str());
     Sample ret = { c, DYee, kMagenta, 1, "dy20ee", true };
     return ret;
}

// Wgamma
Sample fWgamma ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Wgamma_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Wgamma, kBlack, 1, "wgamma", true };
     return ret;
}

// Zgamma
Sample fZgamma ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/Zgamma_Summer08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, Zgamma, kBlack, 1, "zgamma", true };
     return ret;
}

//ttbar file
Sample fttbar ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbar", true, 0. };
     return ret;
}

//ttbar file
Sample fttbar_taula ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/TauolaTTbar-Pythia/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, ttbar, kYellow, 1, "ttbartauola", true, 0. };
     return ret;
}

Sample ftW ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "/cms2-V01-02-06/SingleTop_tWChannel-madgraph-LHE/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "tw", true, 0. };
     return ret;
}

Sample fSingleTop_tChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SingleTop_tChannel-madgraph-LHE/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletopt", true, 0. };
     return ret;
}

Sample fSingleTop_sChannel ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SingleTop_sChannel-madgraph-LHE/merged_ntuple*.root";
     c->Add(sample.c_str());
     Sample ret = { c, tW, 63, 1, "singletops", true, 0. };
     return ret;
}


Sample fLM1 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM1-sftsht/skim/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM1, 37, 1, "LM1", false, 0. };
     return ret;
}

Sample fLM2 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM2-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM2, 37, 1, "LM2", false, 0. };
     return ret;
}

Sample fLM3 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM3-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM3, 37, 1, "LM3", false, 0. };
     return ret;
}

Sample fLM4 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM4-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM4, 37, 1, "LM4", false, 0. };
     return ret;
}

Sample fLM5 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM5-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM5, 37, 1, "LM5", false, 0. };
     return ret;
}

Sample fLM6 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM6-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM6, 37, 1, "LM6", false, 0. };
     return ret;
}

Sample fLM7 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM7-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM7, 37, 1, "LM7", false, 0. };
     return ret;
}

Sample fLM8 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM8-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM8, 37, 1, "LM8", false, 0. };
     return ret;
}

Sample fLM9 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM9-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM9, 37, 1, "LM9", false, 0. };
     return ret;
}

Sample fLM10 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM10-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM10, 37, 1, "LM10", false, 0. };
     return ret;
}

Sample fLM11 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/SUSY_LM11-sftsht/merge*.root";
     c->Add(sample.c_str());
     Sample ret = { c, LM11, 37, 1, "LM11", false, 0. };
     return ret;
}



// QCD samples
Sample fInclusiveMu5Pt50 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/InclusiveMu5Pt50/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMu5Pt50, 28, 1, "InclusiveMu5Pt50", true, 0. };
     return ret;
}

Sample fInclusiveMuPt15 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/InclusiveMuPt15/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, InclusiveMuPt15, 28, 1, "InclusiveMuPt15", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_BCtoE_Pt20to30/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt20to30, 28, 1, "QCDBCtoEPt20to30", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_BCtoE_Pt30to80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt30to80, 28, 1, "QCDBCtoEPt30to80", true, 0. };
     return ret;
}

Sample fQCDBCtoEPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_BCtoE_Pt80to170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDBCtoEPt80to170, 28, 1, "QCDBCtoEPt80to170", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt20to30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_EMenriched_Pt20to30/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt20to30, 28, 1, "QCDEMenrichedPt20to30", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_EMenriched_Pt30to80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt30to80, 28, 1, "QCDEMenrichedPt30to80", true, 0. };
     return ret;
}

Sample fQCDEMenrichedPt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCD_EMenriched_Pt80to170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDEMenrichedPt80to170, 28, 1, "QCDEMenrichedPt80to170", true, 0. };
     return ret;
}

Sample fQCDpt30 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt30_v2/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt30, 28, 1, "QCDpt30", true, 0. };
     return ret;
}

Sample fQCDpt30to80 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt30_v2/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt30to80, 28, 1, "QCDpt30to80", true, 80. };
     return ret;
}

Sample fQCDpt80to170 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt80/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt80to170, 28, 1, "QCDpt80to170", true, 170 };
     return ret;
}

Sample fQCDpt170to300 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt170/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt170to300, 28, 1, "QCDpt170to300", true, 300 };
     return ret;
}
Sample fQCDpt300to470 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt300/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt300to470, 28, 1, "QCDpt300to470", true, 470 };
     return ret;
}
Sample fQCDpt470to800 ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt470/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt470to800, 28, 1, "QCDpt470to800", true, 800 };
     return ret;
}
Sample fQCDpt800toInf ()
{
     TChain *c = new TChain("Events");
     std::string sample = prefix + "cms2-V01-02-06/QCDpt800/*.root";
     c->Add(sample.c_str());
     Sample ret = { c, QCDpt800toInf, 28, 1, "QCDpt800toInf", true, 999999999 };
     return ret;
}
