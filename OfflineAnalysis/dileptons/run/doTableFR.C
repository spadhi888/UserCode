#include <iostream>
#include <vector>
#include <iomanip>

#include "TH1F.h"
#include "TH1F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"

using namespace std;

void doTableFR() {

//  gROOT->ProcessLine(".x Tools/setup.C");

  TFile *ftt = TFile::Open("ntuple.root");

  pair<Float_t, Float_t> tt_mm = doErrors(ftt, "ttbar_hnJet_mm");
  pair<Float_t, Float_t> tt_ee = doErrors(ftt, "ttbar_hnJet_ee");
  pair<Float_t, Float_t> tt_em = doErrors(ftt, "ttbar_hnJet_em");
  pair<Float_t, Float_t> tt_all = doErrors(ftt, "ttbar_hnJet_all");

  pair<Float_t, Float_t> wj_mm = doErrors(ftt, "wj_hnJet_mm");
  pair<Float_t, Float_t> wj_ee = doErrors(ftt, "wj_hnJet_ee");
  pair<Float_t, Float_t> wj_em = doErrors(ftt, "wj_hnJet_em");
  pair<Float_t, Float_t> wj_all = doErrors(ftt, "wj_hnJet_all");

  pair<Float_t, Float_t> tw_mm = doErrors(ftt, "tw_hnJet_mm");
  pair<Float_t, Float_t> tw_ee = doErrors(ftt, "tw_hnJet_ee");
  pair<Float_t, Float_t> tw_em = doErrors(ftt, "tw_hnJet_em");
  pair<Float_t, Float_t> tw_all = doErrors(ftt, "tw_hnJet_all");

  pair<Float_t, Float_t> qcd_mm = doErrors(ftt, "vgamma_hnJet_mm");
  pair<Float_t, Float_t> qcd_ee = doErrors(ftt, "vgamma_hnJet_ee");
  pair<Float_t, Float_t> qcd_em = doErrors(ftt, "vgamma_hnJet_em");
  pair<Float_t, Float_t> qcd_all = doErrors(ftt, "vgamma_hnJet_all");


//  pair<Float_t, Float_t> qcd_mm = doErrors(ftt, "qcd30_hnJet_mm");
//  pair<Float_t, Float_t> qcd_ee = doErrors(ftt, "qcd30_hnJet_ee");
//  pair<Float_t, Float_t> qcd_em = doErrors(ftt, "qcd30_hnJet_em");
//  pair<Float_t, Float_t> qcd_all = doErrors(ftt, "qcd30_hnJet_all");

//  pair<Float_t, Float_t> qcd1_mm = doErrors(ftt, "qcd15_hnJet_mm");
//  pair<Float_t, Float_t> qcd1_ee = doErrors(ftt, "qcd15_hnJet_ee");
//  pair<Float_t, Float_t> qcd1_em = doErrors(ftt, "qcd15_hnJet_em");
//  pair<Float_t, Float_t> qcd1_all = doErrors(ftt, "qcd15_hnJet_all");

  pair<Float_t, Float_t> DYeemm = doErrors(ftt, "DYee_hnJet_mm");
  pair<Float_t, Float_t> DYeeee = doErrors(ftt, "DYee_hnJet_ee");
  pair<Float_t, Float_t> DYeeem = doErrors(ftt, "DYee_hnJet_em");
  pair<Float_t, Float_t> DYeeall = doErrors(ftt, "DYee_hnJet_all");

  pair<Float_t, Float_t> DYmmmm = doErrors(ftt, "DYmm_hnJet_mm");
  pair<Float_t, Float_t> DYmmee = doErrors(ftt, "DYmm_hnJet_ee");
  pair<Float_t, Float_t> DYmmem = doErrors(ftt, "DYmm_hnJet_em");
  pair<Float_t, Float_t> DYmmall = doErrors(ftt, "DYmm_hnJet_all");

  pair<Float_t, Float_t> DYtautaumm = doErrors(ftt, "DYtautau_hnJet_mm");
  pair<Float_t, Float_t> DYtautauee = doErrors(ftt, "DYtautau_hnJet_ee");
  pair<Float_t, Float_t> DYtautauem = doErrors(ftt, "DYtautau_hnJet_em");
  pair<Float_t, Float_t> DYtautauall = doErrors(ftt, "DYtautau_hnJet_all");

  pair<Float_t, Float_t> ww_mm = doErrors(ftt, "WW_hnJet_mm");
  pair<Float_t, Float_t> ww_ee = doErrors(ftt, "WW_hnJet_ee");
  pair<Float_t, Float_t> ww_em = doErrors(ftt, "WW_hnJet_em");
  pair<Float_t, Float_t> ww_all = doErrors(ftt, "WW_hnJet_all");

//  pair<Float_t, Float_t> Vqq_mm = doErrors(ftt, "Vqq_hnJet_mm");
//  pair<Float_t, Float_t> Vqq_ee = doErrors(ftt, "Vqq_hnJet_ee");
//  pair<Float_t, Float_t> Vqq_em = doErrors(ftt, "Vqq_hnJet_em");
//  pair<Float_t, Float_t> Vqq_all = doErrors(ftt, "Vqq_hnJet_all");

//  pair<Float_t, Float_t> wc_mm = doErrors(ftt, "wc_hnJet_mm");
//  pair<Float_t, Float_t> wc_ee = doErrors(ftt, "wc_hnJet_ee");
//  pair<Float_t, Float_t> wc_em = doErrors(ftt, "wc_hnJet_em");
//  pair<Float_t, Float_t> wc_all = doErrors(ftt, "wc_hnJet_all");

  pair<Float_t, Float_t> wz_mm = doErrors(ftt, "wz_hnJet_mm");
  pair<Float_t, Float_t> wz_ee = doErrors(ftt, "wz_hnJet_ee");
  pair<Float_t, Float_t> wz_em = doErrors(ftt, "wz_hnJet_em");
  pair<Float_t, Float_t> wz_all = doErrors(ftt, "wz_hnJet_all");

  pair<Float_t, Float_t> zz_mm = doErrors(ftt, "zz_hnJet_mm");
  pair<Float_t, Float_t> zz_ee = doErrors(ftt, "zz_hnJet_ee");
  pair<Float_t, Float_t> zz_em = doErrors(ftt, "zz_hnJet_em");
  pair<Float_t, Float_t> zz_all = doErrors(ftt, "zz_hnJet_all");


//  float QCD_ee = qcd_ee.first + qcd1_ee.first;
//  float QCD_mm = qcd_mm.first + qcd1_mm.first;
//  float QCD_em = qcd_em.first + qcd1_em.first;
//  float QCD_all = QCD_ee + QCD_mm + QCD_em;

//  float QCDE_ee = sqrt(pow(qcd_ee.second,2)+pow(qcd1_ee.second,2)); 
//  float QCDE_mm = sqrt(pow(qcd_mm.second,2)+pow(qcd1_mm.second,2)); 
//  float QCDE_em = sqrt(pow(qcd_em.second,2)+pow(qcd1_em.second,2)); 
//  float QCDE_all = sqrt(pow(qcd_all.second,2)+pow(qcd1_all.second,2)); 

  float DY_ee = DYeeee.first + DYmmee.first + DYtautauee.first;
  float DY_mm = DYeemm.first + DYmmmm.first + DYtautaumm.first;
  float DY_em = DYeeem.first + DYmmem.first + DYtautauem.first;
  float DY_all = DY_ee + DY_mm + DY_em;

  float DY_eeErr = sqrt(pow(DYeeee.second,2)+pow(DYmmee.second,2)+pow(DYtautauee.second,2));
  float DY_mmErr = sqrt(pow(DYeemm.second,2)+pow(DYmmmm.second,2)+pow(DYtautaumm.second,2));
  float DY_emErr = sqrt(pow(DYeeem.second,2)+pow(DYmmem.second,2)+pow(DYtautauem.second,2));
  float DY_allErr = sqrt(pow(DYeeall.second,2)+pow(DYmmall.second,2)+pow(DYtautauall.second,2));

  // string pm = "+/-";
  string pm = " &plusmn; ";
  //  string pm = " \\pm ";


  cout << "| SS Leptons | ttbar | SingleTop | Wjets | DY | WW | WZ | ZZ |" <<  endl;
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(2);
  cout << "| ee | " << tt_ee.first << pm << tt_ee.second << " | " 
       << tw_ee.first << pm << tw_ee.second << " | " 
       << wj_ee.first << pm << wj_ee.second << " | "
       << DY_ee << pm << DY_eeErr << " | "
//       << QCD_ee << pm << QCDE_ee  << " | "
       << ww_ee.first << pm << ww_ee.second << " | " 
//       << Vqq_ee.first << pm << Vqq_ee.second << " | " 
//       << wc_ee.first << pm << wc_ee.second << " | " 
       << wz_ee.first << pm << wz_ee.second << " | " 
       << zz_ee.first << pm << zz_ee.second << " | " 
       << vgamma_ee.first << pm << vgamma_ee.second << " | " 
       << endl;
  cout << "| &mu;&mu; | " << tt_mm.first << pm << tt_mm.second << " | "
       << tw_mm.first << pm << tw_mm.second << " | "
       << wj_mm.first << pm << wj_mm.second << " | "
       << DY_mm << pm << DY_mmErr << " | "
//       << QCD_mm << pm << QCDE_mm << " | "
       << ww_mm.first << pm << ww_mm.second << " | "
//       << Vqq_mm.first << pm << Vqq_mm.second << " | "
//       << wc_mm.first << pm << wc_mm.second << " | "
       << wz_mm.first << pm << wz_mm.second << " | "
       << zz_mm.first << pm << zz_mm.second << " | "
       << vgamma_mm.first << pm << vgamma_mm.second << " | "
       << endl;
  cout << "| e&mu; | " << tt_em.first << pm << tt_em.second << " | "
       << tw_em.first << pm << tw_em.second << " | "
       << wj_em.first << pm << wj_em.second << " | "
       << DY_em << pm << DY_emErr << " | "
//       << QCD_em << pm << QCDE_em << " | "
       << ww_em.first << pm << ww_em.second << " | "
//       << Vqq_em.first << pm << Vqq_em.second << " | "
//       << wc_em.first << pm << wc_em.second << " | "
       << wz_em.first << pm << wz_em.second << " | "
       << zz_em.first << pm << zz_em.second << " | "
       << vgamma_em.first << pm << vgamma_em.second << " | "
       << endl;
  cout << "| total | " << tt_ee.first+tt_mm.first+tt_em.first << pm << tt_all.second << " | "
       << tw_ee.first+tw_mm.first+tw_em.first << pm << tw_all.second << " | "
       << wj_ee.first+wj_mm.first+wj_em.first << pm << wj_all.second << " | "
       << DY_all << pm << DY_allErr << " | "
//       << QCD_all << pm << QCDE_all << " | "
       << ww_ee.first+ww_mm.first+ww_em.first << pm << ww_all.second << " | "
//       << Vqq_ee.first+Vqq_mm.first+Vqq_em.first << pm << Vqq_all.second << " | "
//       << wc_ee.first+wc_mm.first+wc_em.first << pm << wc_all.second << " | "
       << wz_ee.first+wz_mm.first+wz_em.first << pm << wz_all.second << " | "
       << zz_ee.first+zz_mm.first+zz_em.first << pm << zz_all.second << " | "
       << vgamma_ee.first+vgamma_mm.first+vgamma_em.first << pm << vgamma_all.second << " | "
       << endl;
}


pair<Float_t, Float_t> doErrors(TFile* ftt, const char* histname) {
//  cout << histname << endl;
  TH1F* h1 = dynamic_cast<TH1F*>(ftt->Get(histname));

  float total = 0.0;
  float  err=  0.0;

  for(int j = 1; j < h1->GetNbinsX() + 1; j++) {
    total = h1->GetBinContent(j) + total;
    err = pow(h1->GetBinError(j),2) + err;
  }
  int aa = h1->GetEntries();
  if (aa > 0 ) { 
    total =  h1->Integral();
    err = total/(sqrt(aa));
//    cout << "  " << histname << "  " << total << " " << err << "  " << aa << endl;
    return make_pair(total, err);
  }
  return make_pair(total, sqrt(err));
}

