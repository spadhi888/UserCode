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

void doTable() {

//  gROOT->ProcessLine(".x Tools/setup.C");

  TFile *ftt = TFile::Open("ntuple.root");

  pair<Float_t, Double_t> tt_mm = doErrors(ftt, "ttbar_hnJet_mm");
  pair<Float_t, Double_t> tt_ee = doErrors(ftt, "ttbar_hnJet_ee");
  pair<Float_t, Double_t> tt_em = doErrors(ftt, "ttbar_hnJet_em");
  pair<Float_t, Double_t> tt_all = doErrors(ftt, "ttbar_hnJet_all");

  pair<Float_t, Double_t> wj_mm = doErrors(ftt, "wj_hnJet_mm");
  pair<Float_t, Double_t> wj_ee = doErrors(ftt, "wj_hnJet_ee");
  pair<Float_t, Double_t> wj_em = doErrors(ftt, "wj_hnJet_em");
  pair<Float_t, Double_t> wj_all = doErrors(ftt, "wj_hnJet_all");

  pair<Float_t, Double_t> tw_mm = doErrors(ftt, "tw_hnJet_mm");
  pair<Float_t, Double_t> tw_ee = doErrors(ftt, "tw_hnJet_ee");
  pair<Float_t, Double_t> tw_em = doErrors(ftt, "tw_hnJet_em");
  pair<Float_t, Double_t> tw_all = doErrors(ftt, "tw_hnJet_all");

//  pair<Float_t, Double_t> qcd_mm = doErrors(ftt, "qcd30_hnJet_mm");
//  pair<Float_t, Double_t> qcd_ee = doErrors(ftt, "qcd30_hnJet_ee");
//  pair<Float_t, Double_t> qcd_em = doErrors(ftt, "qcd30_hnJet_em");
//  pair<Float_t, Double_t> qcd_all = doErrors(ftt, "qcd30_hnJet_all");

//  pair<Float_t, Double_t> qcd1_mm = doErrors(ftt, "qcd15_hnJet_mm");
//  pair<Float_t, Double_t> qcd1_ee = doErrors(ftt, "qcd15_hnJet_ee");
//  pair<Float_t, Double_t> qcd1_em = doErrors(ftt, "qcd15_hnJet_em");
//  pair<Float_t, Double_t> qcd1_all = doErrors(ftt, "qcd15_hnJet_all");

  pair<Float_t, Double_t> DYeemm = doErrors(ftt, "DYee_hnJet_mm");
  pair<Float_t, Double_t> DYeeee = doErrors(ftt, "DYee_hnJet_ee");
  pair<Float_t, Double_t> DYeeem = doErrors(ftt, "DYee_hnJet_em");
  pair<Float_t, Double_t> DYeeall = doErrors(ftt, "DYee_hnJet_all");

  pair<Float_t, Double_t> DYmmmm = doErrors(ftt, "DYmm_hnJet_mm");
  pair<Float_t, Double_t> DYmmee = doErrors(ftt, "DYmm_hnJet_ee");
  pair<Float_t, Double_t> DYmmem = doErrors(ftt, "DYmm_hnJet_em");
  pair<Float_t, Double_t> DYmmall = doErrors(ftt, "DYmm_hnJet_all");

  pair<Float_t, Double_t> DYtautaumm = doErrors(ftt, "DYtautau_hnJet_mm");
  pair<Float_t, Double_t> DYtautauee = doErrors(ftt, "DYtautau_hnJet_ee");
  pair<Float_t, Double_t> DYtautauem = doErrors(ftt, "DYtautau_hnJet_em");
  pair<Float_t, Double_t> DYtautauall = doErrors(ftt, "DYtautau_hnJet_all");

  pair<Float_t, Double_t> ww_mm = doErrors(ftt, "WW_hnJet_mm");
  pair<Float_t, Double_t> ww_ee = doErrors(ftt, "WW_hnJet_ee");
  pair<Float_t, Double_t> ww_em = doErrors(ftt, "WW_hnJet_em");
  pair<Float_t, Double_t> ww_all = doErrors(ftt, "WW_hnJet_all");

//  pair<Float_t, Double_t> Vqq_mm = doErrors(ftt, "Vqq_hnJet_mm");
//  pair<Float_t, Double_t> Vqq_ee = doErrors(ftt, "Vqq_hnJet_ee");
//  pair<Float_t, Double_t> Vqq_em = doErrors(ftt, "Vqq_hnJet_em");
//  pair<Float_t, Double_t> Vqq_all = doErrors(ftt, "Vqq_hnJet_all");

//  pair<Float_t, Double_t> wc_mm = doErrors(ftt, "wc_hnJet_mm");
//  pair<Float_t, Double_t> wc_ee = doErrors(ftt, "wc_hnJet_ee");
//  pair<Float_t, Double_t> wc_em = doErrors(ftt, "wc_hnJet_em");
//  pair<Float_t, Double_t> wc_all = doErrors(ftt, "wc_hnJet_all");

  pair<Float_t, Double_t> wz_mm = doErrors(ftt, "wz_hnJet_mm");
  pair<Float_t, Double_t> wz_ee = doErrors(ftt, "wz_hnJet_ee");
  pair<Float_t, Double_t> wz_em = doErrors(ftt, "wz_hnJet_em");
  pair<Float_t, Double_t> wz_all = doErrors(ftt, "wz_hnJet_all");

  pair<Float_t, Double_t> zz_mm = doErrors(ftt, "zz_hnJet_mm");
  pair<Float_t, Double_t> zz_ee = doErrors(ftt, "zz_hnJet_ee");
  pair<Float_t, Double_t> zz_em = doErrors(ftt, "zz_hnJet_em");
  pair<Float_t, Double_t> zz_all = doErrors(ftt, "zz_hnJet_all");

  pair<Float_t, Double_t> vgamma_mm = doErrors(ftt, "vgamma_hnJet_mm");
  pair<Float_t, Double_t> vgamma_ee = doErrors(ftt, "vgamma_hnJet_ee");
  pair<Float_t, Double_t> vgamma_em = doErrors(ftt, "vgamma_hnJet_em");
  pair<Float_t, Double_t> vgamma_all = doErrors(ftt, "vgamma_hnJet_all");



  float DY_ee = DYeeee.first + DYmmee.first + DYtautauee.first;
  float DY_mm = DYeemm.first + DYmmmm.first + DYtautaumm.first;
  float DY_em = DYeeem.first + DYmmem.first + DYtautauem.first;
  float DY_all = DYeeall.first + DYmmall.first + DYtautauall.first;

  float DY_eeErr = sqrt(pow(DYeeee.second,2)+pow(DYmmee.second,2)+pow(DYtautauee.second,2));
  float DY_mmErr = sqrt(pow(DYeemm.second,2)+pow(DYmmmm.second,2)+pow(DYtautaumm.second,2));
  float DY_emErr = sqrt(pow(DYeeem.second,2)+pow(DYmmem.second,2)+pow(DYtautauem.second,2));
  float DY_allErr = sqrt(pow(DYeeall.second,2)+pow(DYmmall.second,2)+pow(DYtautauall.second,2));


//  float QCD_ee = qcd_ee.first + qcd1_ee.first;
//  float QCD_mm = qcd_mm.first + qcd1_mm.first;
//  float QCD_em = qcd_em.first + qcd1_em.first;
//  float QCD_all = qcd_all.first + qcd1_all.first;

//  float QCDE_ee = sqrt(pow(qcd_ee.second,2)+pow(qcd1_ee.second,2));
//  float QCDE_mm = sqrt(pow(qcd_mm.second,2)+pow(qcd1_mm.second,2));
//  float QCDE_em = sqrt(pow(qcd_em.second,2)+pow(qcd1_em.second,2));
//  float QCDE_all = sqrt(pow(qcd_all.second,2)+pow(qcd1_all.second,2));

  float total_ee = tt_ee.first + tw_ee.first + wj_ee.first + DY_ee + ww_ee.first + wz_ee.first + zz_ee.first;
  float totalE_ee = sqrt(pow(tt_ee.second,2)+ pow(tw_ee.second,2) + pow(wj_ee.second,2) + pow(DY_eeErr,2) + pow(ww_ee.second,2) + pow(wz_ee.second,2) + pow(zz_ee.second,2));

  float total_mm = tt_mm.first + tw_mm.first + wj_mm.first + DY_mm + ww_mm.first + wz_mm.first + zz_mm.first;
  float totalE_mm = sqrt(pow(tt_mm.second,2)+ pow(tw_mm.second,2) + pow(wj_mm.second,2) + pow(DY_mmErr,2) + pow(ww_mm.second,2) + pow(wz_mm.second,2) + pow(zz_mm.second,2));

  float total_em = tt_em.first + tw_em.first + wj_em.first + DY_em + ww_em.first + wz_em.first + zz_em.first;
  float totalE_em = sqrt(pow(tt_em.second,2)+ pow(tw_em.second,2) + pow(wj_em.second,2) + pow(DY_emErr,2) + pow(ww_em.second,2) + pow(wz_em.second,2) + pow(zz_em.second,2));
  
  float total_all = tt_all.first + tw_all.first + wj_all.first + DY_all + ww_all.first + wz_all.first + zz_all.first;
  float totalE_all = sqrt(pow(tt_all.second,2)+ pow(tw_all.second,2) + pow(wj_all.second,2) + pow(DY_allErr,2) + pow(ww_all.second,2) + pow(wz_all.second,2) + pow(zz_all.second,2));

  // string pm = "+/-";
  string pm = " &plusmn; ";
  //  string pm = " \\pm ";


  cout << "| SS Leptons | Total MC | ttbar | SingleTop | Wjets | DY | WW | WZ | ZZ | Vgamma (not in total) |" <<  endl;
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(2);
  cout << "| ee | " << total_ee << pm << totalE_ee << " | "
       << tt_ee.first << pm << tt_ee.second << " | " 
       << tw_ee.first << pm << tw_ee.second << " | " 
       << wj_ee.first << pm << wj_ee.second << " | "
       << DY_ee << pm << DY_eeErr << " | "
//       << QCD_ee << pm << QCDE_ee << " | "
       << ww_ee.first << pm << ww_ee.second << " | " 
//       << Vqq_ee.first << pm << Vqq_ee.second << " | " 
//       << wc_ee.first << pm << wc_ee.second << " | " 
       << wz_ee.first << pm << wz_ee.second << " | " 
       << zz_ee.first << pm << zz_ee.second << " | " 
       << vgamma_ee.first << pm << vgamma_ee.second << " | " 
       << endl;
  cout << "| &mu;&mu; | "  << total_mm << pm << totalE_mm << " | "
       << tt_mm.first << pm << tt_mm.second << " | "
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
  cout << "| e&mu; | "  << total_em << pm << totalE_em << " | "
       << tt_em.first << pm << tt_em.second << " | "
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
  cout << "| total | "  << total_all << pm << totalE_all << " | "
       << tt_all.first << pm << tt_all.second << " | "
       << tw_all.first << pm << tw_all.second << " | "
       << wj_all.first << pm << wj_all.second << " | "
       << DY_all << pm << DY_allErr << " | "
//       << QCD_all << pm << QCDE_all << " | "
       << ww_all.first << pm << ww_all.second << " | "
//       << Vqq_all.first << pm << Vqq_all.second << " | "
//       << wc_all.first << pm << wc_all.second << " | "
       << wz_all.first << pm << wz_all.second << " | "
       << zz_all.first << pm << zz_all.second << " | "
       << vgamma_all.first << pm << vgamma_all.second << " | "
       << endl;

}


pair<Float_t, Double_t> doErrors(TFile* ftt, const char* histname) {
//  cout << histname << endl;
  TH1F* h1 = dynamic_cast<TH1F*>(ftt->Get(histname));

  float total = 0.0;
  float  err=  0.0;

  for(int j = 1; j < h1->GetNbinsX() + 1; j++) {
    total = h1->GetBinContent(j) + total;
    err = pow(h1->GetBinError(j),2) + err;
  }
  return make_pair(total, sqrt(err));
}

