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

  // cross limit section maps
  TH2F* hxsec = new TH2F("hxsec","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
  TH2F* hxsecup = new TH2F("hxsecup","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);
  TH2F* hxsecdwn = new TH2F("hxsecdwn","hxsec",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

  // acceptance map
  TH2F* acc = new TH2F("acc","acc",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);

  //an empty histogram
  TH2F* empyt = new TH2F("empty","empty",nx,xmin,xmin+nx*xbinsize,ny,ymin,ymin+ny*ybinsize);


  tree->Draw("lspmass:glmass>>ul","explimsrb/(4680.*effsrb)");
  tree->Draw("lspmass:glmass>>ulbest","bestsr");
  tree->Draw("lspmass:glmass>>acc","100.*effsrb");
  tree->Draw("lspmass:glmass>>excl","explimsrb/(4680.*effsrb)<xsec");
  tree->Draw("lspmass:glmass>>exclup","explimsrb/(4680.*effsrb)<xsecup");
  tree->Draw("lspmass:glmass>>excldwn","explimsrb/(4680.*effsrb)<xsecdwn");
  tree->Draw("lspmass:glmass>>hxsec",   "xsec");
  tree->Draw("lspmass:glmass>>hxsecup", "xsecup");
  tree->Draw("lspmass:glmass>>hxsecdwn","xsecdwn");


  // now scan the upper limit histogram and interpolate the points for the limit
  // - scan in X as a function of Y
  // - as soon as you find a cell which is excluded, but the next cell is not excluded stop
  // - interpolate the value of X, and do a cout of the x and y coordinates
  // Do not trust it too much, it always has to be cheked by hand
  // We start scanning from the top
  vector<float> xvec;
  vector<float> yvec;
  vector<float> xvecup;
  vector<float> yvecup;
  vector<float> xvecdwn;
  vector<float> yvecdwn;
  for (int ixsec=0; ixsec<3; ixsec++) {
    for (int iy=ny; iy>=1; iy--) {
      float y = ymin + (iy-0.5)*ybinsize;
      for (int ix=1; ix<nx; ix++) {
	float x = xmin + (ix-0.5)*xbinsize;
	// cout << ix << " " << x << " " << excl->GetBinContent(ix,1) << endl;
	float thisBinLimit = ul->GetBinContent(ix,iy);
	float thisXsec;
	if (ixsec ==0) thisXsec     = hxsec->GetBinContent(ix,iy);
	if (ixsec ==1) thisXsec     = hxsecup->GetBinContent(ix,iy);
	if (ixsec ==2) thisXsec     = hxsecdwn->GetBinContent(ix,iy);
	if (thisBinLimit < thisXsec) {
	  float nextBinLimit = ul->GetBinContent(ix+1,iy);
	  float nextXsec;
	  if (ixsec ==0) nextXsec     = hxsec->GetBinContent(ix+1,iy);
	  if (ixsec ==1) nextXsec     = hxsecup->GetBinContent(ix+1,iy);
	  if (ixsec ==2) nextXsec     = hxsecdwn->GetBinContent(ix+1,iy);
	  if (nextBinLimit > nextXsec) {
	    float d = nextBinLimit - thisBinLimit - nextXsec + thisXsec;
	    float xlim = (thisXsec - thisBinLimit)*xbinsize/d + x;
	    cout << xlim << " " << y << endl;
	    // fill the vectors that define the line	    
	    if (ixsec == 0) {
	      xvec.push_back(xlim);
	      yvec.push_back(y);
	    } else if (ixsec==1) {
	      xvecup.push_back(xlim);
	      yvecup.push_back(y);
	    } else if (ixsec==2) {
	      xvecdwn.push_back(xlim);
	      yvecdwn.push_back(y);
	    }
	    continue;
	  }
	}
      }
    }
  }
  // Add points by hand to make it go down to the bottom
  float yy = yvec.at(yvec.size()-1) - ybinsize/2.;
  float xx = xvec.at(xvec.size()-1) + 0.5*(xvec.at(xvec.size()-1) - xvec.at(xvec.size()-2));
  xvec.push_back(xx);
  yvec.push_back(yy);
  yy = yvecup.at(yvecup.size()-1) - ybinsize/2.;
  xx = xvecup.at(xvecup.size()-1) + 0.5*(xvecup.at(xvecup.size()-1) - xvecup.at(xvecup.size()-2));
  xvecup.push_back(xx);
  yvecup.push_back(yy);
  yy = yvecdwn.at(yvecdwn.size()-1) - ybinsize/2.;
  xx = xvecdwn.at(xvecdwn.size()-1) + 0.5*(xvecdwn.at(xvecdwn.size()-1) - xvecdwn.at(xvecdwn.size()-2));
  xvecup.push_back(xx);
  yvecup.push_back(yy);


  // The points are stored in a TGraph
  TGraph* g    = new TGraph(xvec.size(),   &xvec[0],   &yvec[0]);
  TGraph* gup  = new TGraph(xvecup.size(), &xvecup[0], &yvecup[0]);
  TGraph* gdwn = new TGraph(xvecdwn.size(),&xvecdwn[0],&yvecdwn[0]);
  g->SetLineColor(4);
  g->SetLineWidth(3);
  gup->SetLineColor(4);
  gup->SetLineWidth(3);
  gup->SetLineStyle(2);
  gdwn->SetLineColor(4);
  gdwn->SetLineWidth(3);
  gdwn->SetLineStyle(2);

  // This line is the kinematical limit
  float mtop = 180.;
  TLine kinlim = TLine(2.*mtop+ymin, ymin, xmax, xmax-2*mtop);

  // Axis labels, a bit primitive for now
  excl->GetYaxis()->SetTitle("LSP Mass (GeV)");
  excl->GetXaxis()->SetTitle("Gluino Mass (GeV)");
  excl->SetTitle("T1tttt model  Excluded points in red");
  ul->GetYaxis()->SetTitle("LSP Mass (GeV)");
  ul->GetXaxis()->SetTitle("Gluino Mass (GeV)");
  ul->SetTitle("T1tttt model  Cross-section upper limits (pb)");
  empty->GetYaxis()->SetTitle("LSP Mass (GeV)");
  empty->GetXaxis()->SetTitle("Gluino Mass (GeV)");
  empty->SetTitle("T1tttt model");
  ulbest->GetYaxis()->SetTitle("LSP Mass (GeV)");
  ulbest->GetXaxis()->SetTitle("Gluino Mass (GeV)");
  ulbest->SetTitle("T1tttt model  Best region based on exp. limit");
  acc->GetYaxis()->SetTitle("LSP Mass (GeV)");
  acc->GetXaxis()->SetTitle("Gluino Mass (GeV)");
  acc->SetTitle("T1tttt model  Acc*Eff*BR for best signal region in percent");





  // Some text
  TLatex gg;
  TLatex gg2;
  TLatex latexLabel;
  latexLabel.SetTextSize(0.035);
  char * selection ="Same Sign with btag selection";
  gg.SetTextSize(0.035);
  gg2.SetTextSize(0.035);
  TLine l1 = TLine(xmin+0.05*(xmax-xmin), ymax-0.175*(ymax-ymin), xmin+0.14*(xmax-xmin), 
                   ymax-0.175*(ymax-ymin));
  l1.SetLineColor(4);
  l1.SetLineWidth(3);
  TLine l2 = TLine(xmin+0.05*(xmax-xmin), ymax-0.225*(ymax-ymin), xmin+0.14*(xmax-xmin), 
                   ymax-0.225*(ymax-ymin));
  l2.SetLineColor(4);
  l2.SetLineWidth(3);
  l2.SetLineStyle(2);




  // Draw the exclusion map and the limit lines to make sure that they make sense
  excl->Draw("col");
  g->Draw("samePC");
  gup->Draw("samePC");
  gdwn->Draw("samePC");
  kinlim.Draw();
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.05*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.10*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.15*(ymax-ymin),"#sqrt{s} = 7 TeV L=4.7 fb^{-1} ");
  gg.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.2*(ymax-ymin), "Exclusion (NLO+NLL xsection)");
  gg2.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.25*(ymax-ymin), "Exclusion (NLO+NLL xsection #pm 1 #sigma)");
  l2.Draw();
  l1.Draw();


  // Draw the limit lines on top of the temperature plot
  TCanvas* c12 = new TCanvas();
  ul->Draw("colz");
  g->Draw("samePC");
  gup->Draw("samePC");
  gdwn->Draw("samePC");
  kinlim.Draw();
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.05*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.10*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.15*(ymax-ymin),"#sqrt{s} = 7 TeV L=4.7 fb^{-1} ");
  gg.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.2*(ymax-ymin), "Exclusion (NLO+NLL xsection)");
  gg2.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.25*(ymax-ymin), "Exclusion (NLO+NLL xsection #pm 1 #sigma)");
  l2.Draw();
  l1.Draw();

  //Draw the limit lines and nothing else
  TCanvas* c13 = new TCanvas();
  empty->Draw();
  g->Draw("samePC");
  gup->Draw("samePC");
  gdwn->Draw("samePC");
  kinlim.Draw();
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.05*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.10*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.15*(ymax-ymin),"#sqrt{s} = 7 TeV L=4.7 fb^{-1} ");
  gg.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.2*(ymax-ymin), "Exclusion (NLO+NLL xsection)");
  gg2.DrawLatex(xmin+0.15*(xmax-xmin), ymax-0.25*(ymax-ymin), "Exclusion (NLO+NLL xsection #pm 1 #sigma)");
  l2.Draw();
  l1.Draw();
   
  //Draw the best region and nothing else
  TCanvas* c14 = new TCanvas();
  ulbest->Draw("text");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.05*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.10*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.15*(ymax-ymin),"#sqrt{s} = 7 TeV L=4.7 fb^{-1} ");

  //Draw the best region and nothing else
  TCanvas* c15 = new TCanvas();
  acc->Draw("colz");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.05*(ymax-ymin),"CMS Preliminary");
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.10*(ymax-ymin),selection);
  latexLabel.DrawLatex(xmin+0.05*(xmax-xmin), ymax-0.15*(ymax-ymin),"#sqrt{s} = 7 TeV L=4.7 fb^{-1} ");
  kinlim.Draw();
}
