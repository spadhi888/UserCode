{
  // Makes basic histograms from the root tree
  // No prettifications yet

  gStyle->SetOptStat(0);

  // This binning insures that the bins are nicely centered
  Double_t xbinsize = 25.;
  Double_t ybinsize = 25.;
  Double_t xmin    = 400. - xbinsize/2.;
  Double_t xmax    = 1100 + xbinsize/2.;
  Double_t ymin    = 50. - ybinsize/2;
  Double_t ymax    = 800. + ybinsize/2.;
  Int_t nx = (xmax-xmin)/xbinsize;
  Int_t ny = (ymax-ymin)/ybinsize;

  // Upper limit as a function of gluino and lsp mass
  TH2F* ul = new TH2F("ul","ul",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

  // Map of excluded points with different xsection assumptions
  TH2F* excl = new TH2F("excl","excl",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
  TH2F* exclup = new TH2F("exclup","exclup",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
  TH2F* excldwn = new TH2F("excldwn","excldwn",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

  // map of best region
  TH2F* ulbest = new TH2F("ulbest","ulbest",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

  tree->Draw("lspmass:glmass>>ul","explimsrb/(4680.*effsrb)");
  tree->Draw("lspmass:glmass>>ulbest","bestsr");
  tree->Draw("lspmass:glmass>>excl","explimsrb/(4680.*effsrb)<xsec");
  tree->Draw("lspmass:glmass>>exclup","explimsrb/(4680.*effsrb)<xsecup");
  tree->Draw("lspmass:glmass>>excldwn","explimsrb/(4680.*effsrb)<xsecdwn");

}
