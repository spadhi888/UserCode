{

TChain *data = new TChain("Events");
TChain *mc = new TChain("Events");
mc->Add("/tas03/disk01/dilephunt_mc/cms2-V03-03-09_QCD_Pt15_Summer09-MC_31X_V3_7TeV-v1/skim/*.root");
data->Add("/tas03/disk01/kalavase/MinimumBias_Commissioning10-Apr20Skim_GOODCOLL-v1/V03-03-19/lepskim/*.root");
data->Add("/tas03/disk01/dilephunt/skim/dilepskim_1338*.root");
gROOT->ProcessLine(".L goodrun.cc++");
gROOT->ProcessLine(".L histtools.C++");

mc->Draw(">>goodmc", "hyp_lt_p4.Pt()>7.&&hyp_ll_p4.Pt()>7.&&hyp_lt_id*hyp_ll_id<0&&hyp_type==3");
data->Draw(">>gooddata", "hyp_lt_p4.Pt()>7.&&hyp_ll_p4.Pt()>7.&&hyp_lt_id*hyp_ll_id<0&&goodrun(evt_run,evt_lumiBlock)&&hyp_type==3");
data->SetEventList(gooddata);
mc->SetEventList(goodmc);


data->Draw("(els_tkJuraIso + els_ecalIso + els_hcalIso)/TMath::Max(els_p4.Pt(),20.)>>d_relIso_e(50,-1,4)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("(els_tkJuraIso + els_ecalIso + els_hcalIso)/TMath::Max(els_p4.Pt(),20.)>>mc_relIso_e(50,-1,4)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");

data->Draw("(els_tkJuraIso + TMath::Max(els_ecalIso-1.,0.) + els_hcalIso)/TMath::Max(els_p4.Pt(),20.)>>d_relIso_b(50,-1,4)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("(els_tkJuraIso + TMath::Max(els_ecalIso-1.,0.) + els_hcalIso)/TMath::Max(els_p4.Pt(),20.)>>mc_relIso_b(50,-1,4)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");


data->Draw("els_tkJuraIso>>d_tkIso_b(50,0,50)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_tkJuraIso>>mc_tkIso_b(50,0,50)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");

data->Draw("els_tkJuraIso>>d_tkIso_e(50,0,50)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_tkJuraIso>>mc_tkIso_e(50,0,50)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");


data->Draw("els_ecalIso>>d_ecalIso_b(35,0,35)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_ecalIso>>mc_ecalIso_b(35,0,35)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");

data->Draw("els_ecalIso>>d_ecalIso_e(35,0,35)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_ecalIso>>mc_ecalIso_e(35,0,35)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");


data->Draw("els_hcalIso>>d_hcalIso_b(20,0,20)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_hcalIso>>mc_hcalIso_b(20,0,20)", "abs(els_etaSC) < 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");

data->Draw("els_hcalIso>>d_hcalIso_e(20,0,20)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");
mc->Draw("els_hcalIso>>mc_hcalIso_e(20,0,20)", "abs(els_etaSC) > 1.479&&els_p4.Pt()>7.&&(els_type&(1<<2))");




TH1F *mc_relIso = (TH1F*)mc_relIso_b->Clone();
mc_relIso->SetName("mc_relIso");
mc_relIso->Add(mc_relIso_e);

TH1F *d_relIso = (TH1F*)d_relIso_b->Clone();
d_relIso->SetName("d_relIso");
d_relIso->Add(d_relIso_e);


TH1F *mc_tkIso = (TH1F*)mc_tkIso_b->Clone();
mc_tkIso->SetName("mc_tkIso");
mc_tkIso->Add(mc_tkIso_e);

TH1F *d_tkIso = (TH1F*)d_tkIso_b->Clone();
d_tkIso->SetName("d_tkIso");
d_tkIso->Add(d_tkIso_e);


TH1F *mc_ecalIso = (TH1F*)mc_ecalIso_b->Clone();
mc_ecalIso->SetName("mc_ecalIso");
mc_ecalIso->Add(mc_ecalIso_e);

TH1F *d_ecalIso = (TH1F*)d_ecalIso_b->Clone();
d_ecalIso->SetName("d_ecalIso");
d_ecalIso->Add(d_ecalIso_e);


TH1F *mc_hcalIso = (TH1F*)mc_hcalIso_b->Clone();
mc_hcalIso->SetName("mc_hcalIso");
mc_hcalIso->Add(mc_hcalIso_e);

TH1F *d_hcalIso = (TH1F*)d_hcalIso_b->Clone();
d_hcalIso->SetName("d_hcalIso");
d_hcalIso->Add(d_hcalIso_e);


hist::saveHist("iso.root");

float scale = d_relIso->GetEntries()/mc_relIso->GetEntries();
mc_relIso->Scale(scale);
d_relIso->SetLineWidth(2);
mc_relIso->SetLineWidth(2);
d_relIso->SetLineColor(kRed);
d_relIso->Draw("e");
mc_relIso->Draw("sames");


scale = d_tkIso->GetEntries()/mc_tkIso->GetEntries();
mc_tkIso->Scale(scale);
d_tkIso->SetLineWidth(2);
mc_tkIso->SetLineWidth(2);
d_tkIso->SetLineColor(kRed);
d_tkIso->Draw("e");
mc_tkIso->Draw("sames");


float scale = d_ecalIso->GetEntries()/mc_ecalIso->GetEntries();
mc_ecalIso->Scale(scale);
d_ecalIso->SetLineWidth(2);
mc_ecalIso->SetLineWidth(2);
d_ecalIso->SetLineColor(kRed);
d_ecalIso->Draw("e");
mc_ecalIso->Draw("sames");


float scale = d_hcalIso->GetEntries()/mc_hcalIso->GetEntries();
mc_hcalIso->Scale(scale);
d_hcalIso->SetLineWidth(2);
mc_hcalIso->SetLineWidth(2);
d_hcalIso->SetLineColor(kRed);
d_hcalIso->Draw("e");
mc_hcalIso->Draw("sames");



}




