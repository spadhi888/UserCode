{

  TFile *f = new TFile("test.root");
  //TFile *f = new TFile("eff_11.root");

  el_pt_effId->SetLineColor(9);
  el_eta_effId->SetLineColor(9);
  mu_pt_effId->SetLineColor(9);
  mu_eta_effId->SetLineColor(9);

  el_pt_effIdIso->SetLineColor(kRed);
  el_eta_effIdIso->SetLineColor(kRed);
  mu_pt_effIdIso->SetLineColor(kRed);
  mu_eta_effIdIso->SetLineColor(kRed);

  double x1 = 0.7;
  double x2 = 0.95;
  double y1 = 0.35;
  double y2 = 0.2;

  TLegend *leg1 = new TLegend(x1,y1,x2,y2);
  leg1->AddEntry(el_pt_effId, "ID");
  leg1->AddEntry(el_pt_effIdIso, "ID & Iso");

  TLegend *leg2 = new TLegend(x1,y1,x2,y2);
  leg2->AddEntry(el_eta_effId, "ID");
  leg2->AddEntry(el_eta_effIdIso, "ID & Iso");

  TLegend *leg3 = new TLegend(x1,y1,x2,y2);
  leg3->AddEntry(mu_pt_effId, "ID");
  leg3->AddEntry(mu_pt_effIdIso, "ID & Iso");

  TLegend *leg4 = new TLegend(x1,y1,x2,y2);
  leg4->AddEntry(mu_eta_effId, "ID");
  leg4->AddEntry(mu_eta_effIdIso, "ID & Iso");

  gStyle->SetOptStat("");

  TCanvas *c = new TCanvas();
  c->SetWindowSize(1100,850);
  c->Divide(2,2);
  c->cd(1);
  el_pt_effId->Draw();
  el_pt_effIdIso->Draw("SAMES");
  leg1->Draw();
  c->cd(2);
  el_eta_effId->Draw();
  el_eta_effIdIso->Draw("SAMES");
  leg2->Draw();
  c->cd(3);
  mu_pt_effId->Draw();
  mu_pt_effIdIso->Draw("SAMES");
  leg3->Draw();
  c->cd(4);
  mu_eta_effId->Draw();
  mu_eta_effIdIso->Draw("SAMES");
  leg4->Draw();

}
