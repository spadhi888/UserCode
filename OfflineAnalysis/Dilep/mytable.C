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


void mytable(){
// Load and compile something to allow proper treatment of vectors
// Not clear that it is needed
//gROOT->LoadMacro("loader.C+");

// Load various tools
//gROOT->ProcessLine(".x ../Tools/setup.C");

//hist file:
TFile *ftt = TFile::Open("./processed_data_tag.root");

 pair<Float_t, Double_t> tt_em = doErrors(ftt, "ttbar_hnJet_em");
 pair<Float_t, Double_t> tt_mm = doErrors(ftt, "ttbar_hnJet_mm");
 pair<Float_t, Double_t> tt_ee = doErrors(ftt, "ttbar_hnJet_ee");
 pair<Float_t, Double_t> tt_all = doErrors(ftt, "ttbar_hnJet_all");
 
 pair<Float_t, Double_t> dy_em = doErrors(ftt, "dy_hnJet_em");
 pair<Float_t, Double_t> dy_mm = doErrors(ftt, "dy_hnJet_mm");
 pair<Float_t, Double_t> dy_ee = doErrors(ftt, "dy_hnJet_ee");
 pair<Float_t, Double_t> dy_all = doErrors(ftt, "dy_hnJet_all");
 
 pair<Float_t, Double_t> wjets_em = doErrors(ftt, "wjets_hnJet_em");
 pair<Float_t, Double_t> wjets_mm = doErrors(ftt, "wjets_hnJet_mm");
 pair<Float_t, Double_t> wjets_ee = doErrors(ftt, "wjets_hnJet_ee");
 pair<Float_t, Double_t> wjets_all = doErrors(ftt, "wjets_hnJet_all");
 
 pair<Float_t, Double_t> wz_em = doErrors(ftt, "wz_hnJet_em");
 pair<Float_t, Double_t> wz_mm = doErrors(ftt, "wz_hnJet_mm");
 pair<Float_t, Double_t> wz_ee = doErrors(ftt, "wz_hnJet_ee");
 pair<Float_t, Double_t> wz_all = doErrors(ftt, "wz_hnJet_all");
 
 pair<Float_t, Double_t> ww_em = doErrors(ftt, "ww_hnJet_em");
 pair<Float_t, Double_t> ww_mm = doErrors(ftt, "ww_hnJet_mm");
 pair<Float_t, Double_t> ww_ee = doErrors(ftt, "ww_hnJet_ee");
 pair<Float_t, Double_t> ww_all = doErrors(ftt, "ww_hnJet_all");
 
 pair<Float_t, Double_t> zz_em = doErrors(ftt, "zz_hnJet_em");
 pair<Float_t, Double_t> zz_mm = doErrors(ftt, "zz_hnJet_mm");
 pair<Float_t, Double_t> zz_ee = doErrors(ftt, "zz_hnJet_ee");
 pair<Float_t, Double_t> zz_all = doErrors(ftt, "zz_hnJet_all");
 
 pair<Float_t, Double_t> tw_em = doErrors(ftt, "tw_hnJet_em");
 pair<Float_t, Double_t> tw_mm = doErrors(ftt, "tw_hnJet_mm");
 pair<Float_t, Double_t> tw_ee = doErrors(ftt, "tw_hnJet_ee");
 pair<Float_t, Double_t> tw_all = doErrors(ftt, "tw_hnJet_all");
 


 float total_ee = tt_ee.first+dy_ee.first+wjets_ee.first+wz_ee.first+ww_ee.first+zz_ee.first+tw_ee.first;
 float total_mm = tt_mm.first+dy_mm.first+wjets_mm.first+wz_mm.first+ww_mm.first+zz_mm.first+tw_mm.first;
 float total_em = tt_em.first+dy_em.first+wjets_em.first+wz_em.first+ww_em.first+zz_em.first+tw_em.first; 
 float total_all = tt_all.first+dy_all.first+wjets_all.first+wz_all.first+ww_all.first+zz_all.first+tw_all.first;


 float total_eeErr = sqrt(pow(tt_ee.second,2)+pow(dy_ee.second,2)+pow(wjets_ee.second,2)+pow(wz_ee.second,2)+pow(ww_ee.second,2)+pow(zz_ee.second,2)+pow(tw_ee.second,2));
 float total_mmErr = sqrt(pow(tt_mm.second,2)+pow(dy_mm.second,2)+pow(wjets_mm.second,2)+pow(wz_mm.second,2)+pow(ww_mm.second,2)+pow(zz_mm.second,2)+pow(tw_mm.second,2));
 float total_emErr = sqrt(pow(tt_em.second,2)+pow(dy_em.second,2)+pow(wjets_em.second,2)+pow(wz_em.second,2)+pow(ww_em.second,2)+pow(zz_em.second,2)+pow(tw_em.second,2));
 float total_allErr = sqrt(pow(tt_all.second,2)+pow(dy_all.second,2)+pow(wjets_all.second,2)+pow(wz_all.second,2)+pow(ww_all.second,2)+pow(zz_all.second,2)+pow(tw_all.second,2));

 // string pm = "+/-";
 string pm = " &plusmn; ";
 //  string pm = " \\pm ";

 cout.setf(ios::fixed, ios::floatfield);
 cout.precision(2);


 cout << "| OS Leptons | Total SM | ttbar | SingleTop | WZ | ZZ | WW | DY | Wjets | " << endl;
 cout << "| ee | " << total_ee << pm << total_eeErr << " | " << tt_ee.first << pm << tt_ee.second << " | " << tw_ee.first << pm << tw_ee.second << " | " << wz_ee.first << pm  << wz_ee.second << " | " << zz_ee.first << pm << zz_ee.second << " | " << ww_ee.first << pm << ww_ee.second << " | " << dy_ee.first << pm  << dy_ee.second << " | " << wjets_ee.first << pm << wjets_ee.second << " | " << endl;
 cout << "| &mu;&mu; | " << total_mm << pm << total_mmErr << " | " << tt_mm.first << pm << tt_mm.second << " | " << tw_mm.first << pm << tw_mm.second << " | " << wz_mm.first << pm  << wz_mm.second << " | " << zz_mm.first << pm << zz_mm.second << " | " << ww_mm.first << pm << ww_mm.second << " | " << dy_mm.first << pm  << dy_mm.second << " | " << wjets_mm.first << pm << wjets_mm.second << " | " << endl;
 cout << "| e&mu; | " << total_em << pm << total_emErr << " | " << tt_em.first << pm << tt_em.second << " | " << tw_em.first << pm << tw_em.second << " | " << wz_em.first << pm  << wz_em.second << " | " << zz_em.first << pm << zz_em.second << " | " << ww_em.first << pm << ww_em.second << " | " << dy_em.first << pm  << dy_em.second << " | " << wjets_em.first << pm << wjets_em.second << " | " << endl;
 cout << "| total | " << total_all << pm << total_allErr << " | "<< tt_all.first << pm << tt_all.second << " | " << tw_all.first << pm << tw_all.second << " | " << wz_all.first << pm  << wz_all.second << " | " << zz_ee.first+zz_mm.first+zz_em.first << pm << zz_all.second << " | " << ww_all.first << pm << ww_all.second << " | " << dy_all.first << pm  << dy_all.second << " | " << wjets_all.first << pm << wjets_all.second << " | " << endl;


 // Now for the LM points

 pair<Float_t, Double_t> lm0x_em = doErrors(ftt, "lm0x_hnJet_em");
 pair<Float_t, Double_t> lm0x_mm = doErrors(ftt, "lm0x_hnJet_mm");
 pair<Float_t, Double_t> lm0x_ee = doErrors(ftt, "lm0x_hnJet_ee");
 pair<Float_t, Double_t> lm0x_all = doErrors(ftt, "lm0x_hnJet_all");


 pair<Float_t, Double_t> lm1x_em = doErrors(ftt, "lm1x_hnJet_em");
 pair<Float_t, Double_t> lm1x_mm = doErrors(ftt, "lm1x_hnJet_mm");
 pair<Float_t, Double_t> lm1x_ee = doErrors(ftt, "lm1x_hnJet_ee");
 pair<Float_t, Double_t> lm1x_all = doErrors(ftt, "lm1x_hnJet_all");

 pair<Float_t, Double_t> lm2x_em = doErrors(ftt, "lm2x_hnJet_em");
 pair<Float_t, Double_t> lm2x_mm = doErrors(ftt, "lm2x_hnJet_mm");
 pair<Float_t, Double_t> lm2x_ee = doErrors(ftt, "lm2x_hnJet_ee");
 pair<Float_t, Double_t> lm2x_all = doErrors(ftt, "lm2x_hnJet_all");

 pair<Float_t, Double_t> lm3x_em = doErrors(ftt, "lm3x_hnJet_em");
 pair<Float_t, Double_t> lm3x_mm = doErrors(ftt, "lm3x_hnJet_mm");
 pair<Float_t, Double_t> lm3x_ee = doErrors(ftt, "lm3x_hnJet_ee");
 pair<Float_t, Double_t> lm3x_all = doErrors(ftt, "lm3x_hnJet_all");


 pair<Float_t, Double_t> lm4x_em = doErrors(ftt, "lm4x_hnJet_em");
 pair<Float_t, Double_t> lm4x_mm = doErrors(ftt, "lm4x_hnJet_mm");
 pair<Float_t, Double_t> lm4x_ee = doErrors(ftt, "lm4x_hnJet_ee");
 pair<Float_t, Double_t> lm4x_all = doErrors(ftt, "lm4x_hnJet_all");


 pair<Float_t, Double_t> lm5x_em = doErrors(ftt, "lm5x_hnJet_em");
 pair<Float_t, Double_t> lm5x_mm = doErrors(ftt, "lm5x_hnJet_mm");
 pair<Float_t, Double_t> lm5x_ee = doErrors(ftt, "lm5x_hnJet_ee");
 pair<Float_t, Double_t> lm5x_all = doErrors(ftt, "lm5x_hnJet_all");


 pair<Float_t, Double_t> lm6x_em = doErrors(ftt, "lm6x_hnJet_em");
 pair<Float_t, Double_t> lm6x_mm = doErrors(ftt, "lm6x_hnJet_mm");
 pair<Float_t, Double_t> lm6x_ee = doErrors(ftt, "lm6x_hnJet_ee");
 pair<Float_t, Double_t> lm6x_all = doErrors(ftt, "lm6x_hnJet_all");


 pair<Float_t, Double_t> lm7x_em = doErrors(ftt, "lm7x_hnJet_em");
 pair<Float_t, Double_t> lm7x_mm = doErrors(ftt, "lm7x_hnJet_mm");
 pair<Float_t, Double_t> lm7x_ee = doErrors(ftt, "lm7x_hnJet_ee");
 pair<Float_t, Double_t> lm7x_all = doErrors(ftt, "lm7x_hnJet_all");


 pair<Float_t, Double_t> lm8x_em = doErrors(ftt, "lm8x_hnJet_em");
 pair<Float_t, Double_t> lm8x_mm = doErrors(ftt, "lm8x_hnJet_mm");
 pair<Float_t, Double_t> lm8x_ee = doErrors(ftt, "lm8x_hnJet_ee");
 pair<Float_t, Double_t> lm8x_all = doErrors(ftt, "lm8x_hnJet_all");


 pair<Float_t, Double_t> lm9x_em = doErrors(ftt, "lm9x_hnJet_em");
 pair<Float_t, Double_t> lm9x_mm = doErrors(ftt, "lm9x_hnJet_mm");
 pair<Float_t, Double_t> lm9x_ee = doErrors(ftt, "lm9x_hnJet_ee");
 pair<Float_t, Double_t> lm9x_all = doErrors(ftt, "lm9x_hnJet_all");

 cout << "| OS Leptons | LM0 | LM1 | LM2 | LM3 | LM4 | LM5 | LM6 | LM7 | LM8 | LM9 |" << endl;

 cout << "| ee | " << lm0x_ee.first << pm << lm0x_ee.second << " | " << lm1x_ee.first << pm << lm1x_ee.second << " | " << lm2x_ee.first << pm << lm2x_ee.second << " | " << lm3x_ee.first << pm << lm3x_ee.second << " | " << lm4x_ee.first << pm << lm4x_ee.second << " | " << lm5x_ee.first << pm << lm5x_ee.second << " | " << lm6x_ee.first << pm << lm6x_ee.second << " | " << lm7x_ee.first << pm << lm7x_ee.second << " | " << lm8x_ee.first << pm << lm8x_ee.second << " | " << lm9x_ee.first  << pm << lm9x_ee.second << " | " << endl;

 cout << "| &mu;&mu; | " << lm0x_mm.first << pm << lm0x_mm.second << " | " << lm1x_mm.first << pm << lm1x_mm.second << " | " << lm2x_mm.first << pm << lm2x_mm.second << " | " << lm3x_mm.first << pm << lm3x_mm.second << " | " << lm4x_mm.first << pm << lm4x_mm.second << " | " << lm5x_mm.first << pm << lm5x_mm.second << " | " << lm6x_mm.first << pm << lm6x_mm.second << " | " << lm7x_mm.first << pm << lm7x_mm.second << " | " << lm8x_mm.first << pm << lm8x_mm.second << " | " << lm9x_mm.first  << pm << lm9x_mm.second << " | " << endl;

 cout << "| e&mu; | " << lm0x_em.first << pm << lm0x_em.second << " | " << lm1x_em.first << pm << lm1x_em.second << " | " << lm2x_em.first << pm << lm2x_em.second << " | " << lm3x_em.first << pm << lm3x_em.second << " | " << lm4x_em.first << pm << lm4x_em.second << " | " << lm5x_em.first << pm << lm5x_em.second << " | " << lm6x_em.first << pm << lm6x_em.second << " | " << lm7x_em.first << pm << lm7x_em.second << " | " << lm8x_em.first << pm << lm8x_em.second << " | " << lm9x_em.first  << pm << lm9x_em.second << " | " << endl;

 cout << "| total | " << lm0x_all.first << pm << lm0x_all.second << " | " << lm1x_all.first << pm << lm1x_all.second << " | " << lm2x_all.first << pm << lm2x_all.second << " | " << lm3x_all.first << pm << lm3x_all.second << " | " << lm4x_all.first << pm << lm4x_all.second << " | " << lm5x_all.first << pm << lm5x_all.second << " | " << lm6x_all.first << pm << lm6x_all.second << " | " << lm7x_all.first << pm << lm7x_all.second << " | " << lm8x_all.first << pm << lm8x_all.second << " | " << lm9x_all.first  << pm << lm9x_all.second << " | " << endl;





 cout << "|Same Sign | SM+LM0 | SM+LM1 | SM+LM2 | SM+LM3 | SM+LM4 | SM+LM5 | SM+LM6 | SM+LM7 | SM+LM8 | SM+LM9 |" <<endl;
 cout << "|Observed |" << lm0x_all.first+total_all << pm << sqrt(pow(lm0x_all.second,2)+pow(total_allErr, 2)) << " | " 
      << lm1x_all.first+total_all << pm << sqrt(pow(lm1x_all.second,2)+pow(total_allErr,2)) << " | " 
      << lm2x_all.first+total_all << pm << sqrt(pow(lm2x_all.second,2)+pow(total_allErr,2)) << " | " 
      << lm3x_all.first+total_all << pm << sqrt(pow(lm3x_all.second,2)+pow(total_allErr,2)) << " | " 
      << lm4x_all.first+total_all << pm << sqrt(pow(lm4x_all.second,2)+pow(total_allErr,2)) << " | "
      << lm5x_all.first+total_all << pm << sqrt(pow(lm5x_all.second,2)+pow(total_allErr,2)) << " | "
      << lm6x_all.first+total_all << pm << sqrt(pow(lm6x_all.second,2)+pow(total_allErr,2)) << " | "
      << lm7x_all.first+total_all << pm << sqrt(pow(lm7x_all.second,2)+pow(total_allErr,2)) << " | "
      << lm8x_all.first+total_all << pm << sqrt(pow(lm8x_all.second,2)+pow(total_allErr,2)) << " | "
      << lm9x_all.first+total_all << pm << sqrt(pow(lm9x_all.second,2)+pow(total_allErr,2)) << " | "
      << endl;
 cout << "| Predicted | 4.10 | 1.26 | 0.83 | 1.24 | 0.91 | 0.81 | 0.84 | 0.83 | 1.03 | 1.00 |" << endl;



}

pair<Float_t, Double_t> doErrors(TFile* ftt, const char* histname) {

  TH1F* h1 = dynamic_cast<TH1F*>(ftt->Get(histname));

  float total = 0.0;
  float  err=  0.0;

  for(int j = 1; j < h1->GetNbinsX() + 1; j++) {
    total = h1->GetBinContent(j) + total;
    err = pow(h1->GetBinError(j),2) + err;
  }
  return make_pair(total, sqrt(err));
}

//  LocalWords:  plusmn
