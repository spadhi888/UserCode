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

void doTTmytableFR() {

//  gROOT->ProcessLine(".x Tools/setup.C");

  TFile *ftt = TFile::Open("ntuple.root");

  float tt[4];
  float ttww[4];
  float ttwo[4];
  float ttwosemi[4];
  float ttwoother[4];
  float other[4];

  for (int i = 0; i < 4; i++){
    tt[i] = 0.0;
    ttww[i] = 0.0;
    ttwosemi[i] = 0.0;
    ttwoother[i] = 0.0;
    other[i] = 0.0;
  }


  pair<Float_t, Float_t> ttww_ee = doErrors(ftt, "ttbar_hnJetWW_ee");
  pair<Float_t, Float_t> ttww_mm = doErrors(ftt, "ttbar_hnJetWW_mm");
  pair<Float_t, Float_t> ttww_em = doErrors(ftt, "ttbar_hnJetWW_em");
  pair<Float_t, Float_t> ttww_all = doErrors(ftt, "ttbar_hnJetWW_all");


  pair<Float_t, Float_t> ttwo_ee = doErrors(ftt, "ttbar_hnJetWO_ee");
  pair<Float_t, Float_t> ttwo_mm = doErrors(ftt, "ttbar_hnJetWO_mm");
  pair<Float_t, Float_t> ttwo_em = doErrors(ftt, "ttbar_hnJetWO_em");
  pair<Float_t, Float_t> ttwo_all = doErrors(ftt, "ttbar_hnJetWO_all");

  pair<Float_t, Float_t> ttwosemi_ee = doErrors(ftt, "ttbar_hnJetWOSemilep_ee");
  pair<Float_t, Float_t> ttwosemi_mm = doErrors(ftt, "ttbar_hnJetWOSemilep_mm");
  pair<Float_t, Float_t> ttwosemi_em = doErrors(ftt, "ttbar_hnJetWOSemilep_em");
  pair<Float_t, Float_t> ttwosemi_all = doErrors(ftt, "ttbar_hnJetWOSemilep_all");

  pair<Float_t, Float_t> ttwosemioth_ee = doErrors(ftt, "ttbar_hnJetWOOther_ee");
  pair<Float_t, Float_t> ttwosemioth_mm = doErrors(ftt, "ttbar_hnJetWOOther_mm");
  pair<Float_t, Float_t> ttwosemioth_em = doErrors(ftt, "ttbar_hnJetWOOther_em");
  pair<Float_t, Float_t> ttwosemioth_all = doErrors(ftt, "ttbar_hnJetWOOther_all");

  pair<Float_t, Float_t> other_ee = doErrors(ftt, "ttbar_hnJetOO_ee");
  pair<Float_t, Float_t> other_mm = doErrors(ftt, "ttbar_hnJetOO_mm");
  pair<Float_t, Float_t> other_em = doErrors(ftt, "ttbar_hnJetOO_em");
  pair<Float_t, Float_t> other_all = doErrors(ftt, "ttbar_hnJetOO_all");

  pair<Float_t, Float_t> tt_mm = doErrors(ftt, "ttbar_hnJet_mm");
  pair<Float_t, Float_t> tt_ee = doErrors(ftt, "ttbar_hnJet_ee");
  pair<Float_t, Float_t> tt_all = doErrors(ftt, "ttbar_hnJet_all");


  float total_ee = ttww_ee.first + ttwo_ee.first + other_ee.first;
  float total_mm = ttww_mm.first + ttwo_mm.first + other_mm.first;
  float total_em = ttww_em.first + ttwo_em.first + other_em.first;
  float total_all = ttww_all.first + ttwo_all.first + other_all.first;

  float total_eeErr = sqrt(pow(ttww_ee.second,2)+pow(ttwo_ee.second,2)+pow(other_ee.second,2));
  float total_mmErr = sqrt(pow(ttww_mm.second,2)+pow(ttwo_mm.second,2)+pow(other_mm.second,2));
  float total_emErr = sqrt(pow(ttww_em.second,2)+pow(ttwo_em.second,2)+pow(other_em.second,2));
  float total_allErr = sqrt(pow(ttww_all.second,2)+pow(ttwo_all.second,2)+pow(other_all.second,2));


  // string pm = "+/-";
  string pm = " &plusmn; ";
  //  string pm = " \\pm ";


  cout << "| SS Leptons | Total | Type-I | Type-II | Type-II a) | Type-II b) | Type-III |" <<  endl;
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(2);
  cout << "| ee | " << total_ee << pm << total_eeErr << " | " 
       << ttww_ee.first << pm << ttww_ee.second << " | " 
       << ttwo_ee.first << pm << ttwo_ee.second << " | "
       << ttwosemi_ee.first << pm << ttwosemi_ee.second << " | "
       << ttwosemioth_ee.first << pm << ttwosemioth_ee.second << " | "
       << other_ee.first << pm << other_ee.second << " | " 
       << endl;
  cout << "| &mu;&mu; | " << total_mm << pm << total_mmErr << " | " 
       << ttww_mm.first << pm << ttww_mm.second << " | "
       << ttwo_mm.first << pm << ttwo_mm.second << " | "
       << ttwosemi_mm.first << pm << ttwosemi_mm.second << " | "
       << ttwosemioth_mm.first << pm << ttwosemioth_mm.second << " | "
       << other_mm.first << pm << other_mm.second << " | "
       << endl;
  cout << "| e&mu; | " << total_em << pm << total_emErr << " | "
       << ttww_em.first << pm << ttww_em.second << " | "
       << ttwo_em.first << pm << ttwo_em.second << " | "
       << ttwosemi_em.first << pm << ttwosemi_em.second << " | "
       << ttwosemioth_em.first << pm << ttwosemioth_em.second << " | "
       << other_em.first << pm << other_em.second << " | "
       << endl;
  
  cout << "| total | " << total_all << pm << total_allErr << " | "

       << ttww_all.first << pm << ttww_all.second << " | "
       << ttwo_all.first << pm << ttwo_all.second << " | "
       << ttwosemi_all.first << pm << ttwosemi_all.second << " | "
       << ttwosemioth_all.first << pm << ttwosemioth_all.second << " | "
       << other_all.first << pm << other_all.second << " | "
       << endl;

  
  tt[0]=ttbar_hnJet_ee->Integral();
  ttww[0] = ttbar_hnJetWW_ee->Integral();
  ttwo[0]=ttbar_hnJetWO_ee->Integral();
  ttwosemi[0]=ttbar_hnJetWOSemilep_ee->Integral();
  ttwoother[0]= ttbar_hnJetWOOther_ee->Integral();
  other[0]=ttbar_hnJetOO_ee->Integral();

  tt[1]=ttbar_hnJet_mm->Integral();
  ttww[1] = ttbar_hnJetWW_mm->Integral();
  ttwo[1]=ttbar_hnJetWO_mm->Integral();
  ttwosemi[1]=ttbar_hnJetWOSemilep_mm->Integral();
  ttwoother[1]= ttbar_hnJetWOOther_mm->Integral();
  other[1]=ttbar_hnJetOO_mm->Integral();

  tt[2]=ttbar_hnJet_em->Integral();
  ttww[2] = ttbar_hnJetWW_em->Integral();
  ttwo[2]=ttbar_hnJetWO_em->Integral();
  ttwosemi[2]=ttbar_hnJetWOSemilep_em->Integral();
  ttwoother[2]= ttbar_hnJetWOOther_em->Integral();
  other[2]=ttbar_hnJetOO_em->Integral();

  tt[3]=ttbar_hnJet_all->Integral();
  ttww[3] = ttbar_hnJetWW_all->Integral();
  ttwo[3]=ttbar_hnJetWO_all->Integral();
  ttwosemi[3]=ttbar_hnJetWOSemilep_all->Integral();
  ttwoother[3]= ttbar_hnJetWOOther_all->Integral();
  other[3]=ttbar_hnJetOO_all->Integral();


  cout << "| SS Leptons | Total | Type-I | Type-II | Type-II a) | Type-II b) | Type-III |" <<  endl;
  cout << "| ee | " << ttbar_hnJet_ee->GetEntries() << " | " << ttbar_hnJetWW_ee->GetEntries() << " | " << ttbar_hnJetWO_ee->GetEntries() << " | " << ttbar_hnJetWOSemilep_ee->GetEntries() << " | " << ttbar_hnJetWOOther_ee->GetEntries() << " | " << ttbar_hnJetOO_ee->GetEntries() << " | " << endl;
  cout << "| &mu;&mu; | " << ttbar_hnJet_mm->GetEntries() << " | " << ttbar_hnJetWW_mm->GetEntries() << " | " << ttbar_hnJetWO_mm->GetEntries() << " | " << ttbar_hnJetWOSemilep_mm->GetEntries() << " | " << ttbar_hnJetWOOther_mm->GetEntries() << " | " << ttbar_hnJetOO_mm->GetEntries() << " | " << endl;
  cout << "| e&mu; | " << ttbar_hnJet_em->GetEntries() << " | " << ttbar_hnJetWW_em->GetEntries() << " | " << ttbar_hnJetWO_em->GetEntries() << " | " << ttbar_hnJetWOSemilep_em->GetEntries() << " | " << ttbar_hnJetWOOther_em->GetEntries() << " | " << ttbar_hnJetOO_em->GetEntries() << " | " << endl;
  cout << "| Total | " << ttbar_hnJet_all->GetEntries() << " | " << ttbar_hnJetWW_all->GetEntries() << " | " << ttbar_hnJetWO_all->GetEntries() << " | " << ttbar_hnJetWOSemilep_all->GetEntries() << " | " << ttbar_hnJetWOOther_all->GetEntries() << " | " << ttbar_hnJetOO_all->GetEntries() << " | " << endl;


  cout << "| Same Sign | Type I+II+III | Type-I | Type-II | Type-III | Type-II-SemiLep | Type-II-Fakes | " <<  endl;
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(2);
  
  for (int i = 0; i < 4 ; i++) {
    cout << "| | " << ttww[i]+ttwo[i]+other[i] << " | " << ttww[i] << " | " << ttwo[i] << " | " << other[i] << " | " << ttwosemi[i] << " | " << ttwoother[i] << " | " << endl;
  }


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


