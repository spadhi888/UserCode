
//#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
//#include "Tools/Utilities.h"
//#include "Tools/histtools.cc"

#include "TROOT.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "/home/users/wandrews/macros/comparison.C"
//#include "comparison.C"

bool verbose = false;

void saveStack(THStack* st, TLegend* leg, bool cpylog, double ymin, double ymax, TString name="", THStack* st2=0) {

  double oldymax = st->GetMaximum();
  //double oldymin = st->GetMinimum();
  TCanvas *c = new TCanvas();
  TString usename;
  if( name != "" ) {
	usename = name;//+".png";
	//st->SetName(name);
	st->SetTitle(name);
  }
  else
	usename = (TString)st->GetName();//+".png";

  //for the y axis won't give me a max...just gives 1--don't use GetYaxis for stacks, use GetMaximum
  //cout << "ymin " << ymin << "  ymax " << ymax << "  st y max " << st->GetYaxis()->GetXmax() << "  max " << st->GetMaximum() << endl;
  if( ymax != -999 ) {
	//st->SetMinimum( ymin );
	st->SetMaximum( ymax );
  }
  //else
  //st->SetMinimum( ymin );
  //st->GetYaxis()->SetRangeUser( ymin, st->GetMaximum()+0.01*st->GetMaximum() );

  st->Draw("hist"); //must draw before the axis range is defined...(uh, don't ask me)
  if( st2 != 0 ) {
	st->SetMinimum( st2->GetMinimum() );
	st2->Draw("same,hist");
  }

  if( verbose )
	cout << "stack " << usename << "  max " << st->GetMaximum() << "  " << oldymax << "  min " << st->GetMinimum() << endl;
  //c->Update();
  //getc(stdin);
  leg->Draw();
  c->SaveAs(usename+".png");

  if(cpylog) {
	gPad->SetLogy();

	//st->SetMaximum( oldymax );
	//if( st2 == 0 ) { //don't reset min if st2 exists--this is wrong
	st->SetMinimum( ymin ); //need to reset regardless
	//cout << "set min to ymin " << ymin << endl;
	//}
  
	st->Draw("hist");
 	leg->Draw();
	c->Update();
	//c->SaveAs((TString)st->GetName()+"_log.png");
	c->SaveAs(usename+"_log.png");
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog=false, double ymin=1., double ymax=-999, TString name="") {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, bool cpylog=false, TString title2="", double scale=1.0, TString name="", double ymin=1., double ymax=-999) {
  THStack *st = h->getSumStack(theSources, title, subtitle, suff1, suff2, 1, title2, scale);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void make2fileStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, bool cpylog=false, double ymin=1., double ymax=-999, TString name="") {
  THStack *st = h->get2fileStack(theSources, title, subtitle, suff1);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}

void makeSumDifStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suff1, TString suff2, TString suff3, bool cpylog=false, TString name="", double ymin=1., double ymax=-999) {
  //second to last argument is scale--leave 1
  THStack *stpos = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3, 1, 1);
  THStack *stneg = h->getSumDifStack(theSources, title, subtitle, suff1, suff2, suff3, 1, -1); //negative
  saveStack( stpos, leg, cpylog, ymin, ymax, name, stneg );
}

void makeHistOverlay( HistogramUtilities* h, sources_t theSources1,  sources_t theSources2, TString title1, TString title2, TString subtitle, TString suff1, TString suff2, bool cpylog=false, bool cpyscale=false, bool cpydiff=false, TString nameprefix="", TString namesuffix="", double ymin=1., double ymax=-999) {
  TH1F* h1 = h->getHistogram(theSources1, title1, "", suff1, 1, nameprefix);
  TH1F* h2 = h->getHistogram(theSources2, title2, "", suff2, 1, nameprefix);
  if( namesuffix != "" ) {
	h1->SetName( (TString)h1->GetName() + namesuffix );
	h2->SetName( (TString)h2->GetName() + namesuffix );
  }
  over_save(h2, h1);
  if( cpylog )
	over_save(h2, h1, true, true); //two bools set log, and put log in name
  //same two scaled
  if( cpyscale ) {
	h2->Scale( h1->Integral()/h2->Integral() );
	h2->SetName((TString)h2->GetName() + "_scale");
	over_save(h2, h1);
	if( cpylog )
	  over_save(h2, h1, true, true);
  }
  //same two difference--keep the scaling if cpyscale
  if( cpydiff ) {
	h1->Add( h2, -1. );
	h2->SetName( (TString)h2->GetName() + "_diff" ); //change name on h2 for consistency with above
	h2->SetTitle( (TString)h2->GetName() + "_diff" ); //for display purposes
	TCanvas *c1 = new TCanvas();
	h1->Draw();
	c1->SaveAs( (TString)h2->GetName() + ".png" );
  }
}

//fit functions

//if you use either low or high, you must use both
double* getHistSigma( TH1F* h, double low=0, double high=0 ) {
  //TF1 f = new TF1("f", "gaus");
  //TF1 *f;
  if( low == 0 && high == 0 ) //don't use range--full range of hist is default
	//f = new TF1("f", "gaus");
	h->Fit("gaus", "Q"); //Q for quiet
  else
	//f = new TF1("f", "gaus", low, high);
	h->Fit("gaus", "Q", "", low, high); //Q for quiet
	
  //h->Fit("gaus", "Q"); //Q for quiet
  //h->Fit( "f", "RQ"); //Q for quiet
  TF1 *f2 = h->GetFunction("gaus");
  double* ret = new double[2]; //need new b'c it's returned
  ret[0] = f2->GetParameter(2);
  ret[1] = f2->GetParError(2);
  return ret;
  //return f2->GetParameter(2); // 0 is a scale constant, 1 is mean, 2 is sigma
  //return {f->GetParameter(2), f->GetParError(2)}; 
	//f->GetChisquare(); //this is useless for now
}

//for gaussians
double getTailRatio( TH1F* h ) {
  h->Fit("gaus");
  TF1 *f = h->GetFunction("gaus");
  double sigma = f->GetParameter(2);
  double mean  = f->GetParameter(1);
  //int binlow = h->GetBin( mean - sigma );
  //int binhgh = h->GetBin( mean + sigma );
  int binlow = h->GetXaxis()->FindFixBin( mean - sigma );
  int binhgh = h->GetXaxis()->FindFixBin( mean + sigma );
  cout << "hist " << h->GetName() << "  mean " << mean << "  sigma " << sigma << "  binlow " << binlow << "  binhigh " << binhgh << endl;
  double intsigma = h->Integral( binlow, binhgh );
  double inttaillow = h->Integral( 0, binlow-1 ); //include underflow
  double inttailhgh = h->Integral( binhgh+1, h->GetNbinsX()+1 ); //include overflow
  return intsigma/(inttaillow + inttailhgh);
}

//I need this function bc root doesn't do integrals with x-axis values, only with bins...fucking root...
double getIntegral( TH1F* h, double low1, double hgh1=0 ) {
  int binlow1 = h->GetXaxis()->FindFixBin( low1 );
  if( hgh1 != 0 ) {
	int binhgh1 = h->GetXaxis()->FindFixBin( hgh1 );
	return h->Integral( binlow1, binhgh1 );
  }
  else 
	return h->Integral( binlow1, h->GetNbinsX()+1 );
}

//ratio function : return ratio of integral(low2,hgh2)/integral(low1,hgh1)
//hgh2=0 is for infinity as up range
double* getIntegralRatio( TH1F* h, double low1, double hgh1, double low2, double hgh2=0 ) { 
  int binlow1 = h->GetXaxis()->FindFixBin( low1 );
  int binhgh1 = h->GetXaxis()->FindFixBin( hgh1 );
  int binlow2 = h->GetXaxis()->FindFixBin( low2 );
  int binhgh2 = 0;
  double int1 = h->Integral( binlow1, binhgh1 );
  double int2 = 0;
  if( hgh2 == 0 ) //default is infinity as upper limit
	int2 = h->Integral( binlow2, h->GetNbinsX()+1 );
  else {
	binhgh2 = h->GetXaxis()->FindFixBin( hgh2 );
	int2 = h->Integral( binlow2, binhgh2 );
  }

  if( int1 == 0 )
	return 0;
  else {
	// f = x/y, err_f/f = sqrt( err_x^2/x^2 + err_y^2/y^2 ) -->> err_f^2 = x/y^2 + x^2/y^3 -->> (f/y)(1+f)
	double itg = int2/int1;
	double *err = new double[2];
	//err = {itg, (itg/int1)*(1 + itg) };
	err[0] = itg;
	err[1] = sqrt( (itg/int1)*(1 + itg) );
	return err;
  }
}


//HSqrt: takes sq root of each bin--since takes ref, changes argument
void HSqrt( TH1F*& h ) {
  for( int i=0; i<h->GetNbinsX(); i++ ) 
	h->SetBinContent( i, sqrt( h->GetBinContent(i) ) );
}


//results function

void plotResults() {

  sources_t theSources = sources_all;
  
  //HistogramUtilities h1("Results.root");
  //HistogramUtilities h1("Susy_Results.root");
  HistogramUtilities* h1 = new HistogramUtilities("Results.root");
  HistogramUtilities* h2 = new HistogramUtilities("Results_Nar.root");

  TString all = "all";
  TString ee = "ee"; //these correspond to defintion in DileptonHypType.h
  TString mm = "mm";
  TString em = "em"; 

  //Need this in order to not store clones, which could have same names in 2 files, thus breaking file_->Get()
  TH1::AddDirectory(false); 

  //makeStack(h1, theSources, "tcMet", "", all);
  
  TLegend *lg_all = h1->getLegend(theSources, "sumJetPt_outz", "", all);

  makeStack(h1, lg_all, theSources, "sumJetPt", "", all, true);
  makeStack(h1, lg_all, theSources, "sumJetPt_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "sumJetPt_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet", "", all, true); //magnitude of met
  makeStack(h1, lg_all, theSources, "tcMet_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_x", "", all, true); //x comp
  makeStack(h1, lg_all, theSources, "tcMet_x_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_x_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_y", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_y_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_xy", "", all, true); //each x and y comps seperately
  makeStack(h1, lg_all, theSources, "tcMet_xy_outz", "", all, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_inz", "", all, true);

  makeStack(h1, lg_all, theSources, "tcMet_xy", "", em, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_outz", "", em, true);
  makeStack(h1, lg_all, theSources, "tcMet_xy_inz", "", em, true);

  makeSumStack(h1, lg_all, theSources, "tcMet_xy", "", ee, mm, true);
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_outz", "", ee, mm, true);
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_inz", "", ee, mm, true);

  makeSumStack(h1, lg_all, theSources, "tcMet_xy"     , "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy");
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_outz", "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy_outz");
  makeSumStack(h1, lg_all, theSources, "tcMet_xy_inz" , "", all, all, true, "genMet_xy", -1., "tcMet-genMet_xy_inz"); 

  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy", "", ee, mm, em, true);
  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy_outz", "", ee, mm, em, true);
  makeSumDifStack(h1, lg_all, theSources, "tcMet_xy_inz", "", ee, mm, em, true);

  //verbose = true;
  //cout << "setting verbose = true" << endl;
  //h1->setVerbose( true );
  //makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy", "", ee, mm, em, true, "nody_tcMet_xy_eemm-em", -20, 370); //just for this one, set range manually
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy", "", ee, mm, em, true, "nody_tcMet_xy_eemm-em");
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy_outz", "", ee, mm, em, true, "nody_tcMet_xy_outz_eemm-em");
  makeSumDifStack(h1, lg_all, sources_nody, "tcMet_xy_inz", "", ee, mm, em, true, "nody_tcMet_xy_inz_eemm-em");
  //verbose = false;
  //cout << "setting verbose = false" << endl;
  //h1->setVerbose( false );

  //integral hists
  TH1F* h_int_ratio[nsjpbins];
  //this value actually comes from looper.h via preprocessor directive...kind of scary, kind of cool
  //int nsjpbins = 6; //nsjpbins is set as 6 here (was 4)
  //nsjpbins = 6;
  double* res_ratio1[nsjpbins];
  TString name_ratio1 = "Ratio_200+_over_50-100_tcMet_inz_eemm-em";
  TH1F* hres_ratio1 = new TH1F(name_ratio1 , name_ratio1, nsjpbins, 0, nsjpbins );
  hres_ratio1->Sumw2();
  
  double* res_ratio2[nsjpbins];
  TString name_ratio2 = "Ratio_100+_over_0-50_tcMet_inz_eemm-em";
  TH1F* hres_ratio2 = new TH1F( name_ratio2, name_ratio2, nsjpbins, 0, nsjpbins );
  hres_ratio2->Sumw2();
  
  double* res_ratio3[nsjpbins];
  TString name_ratio3 = "Ratio_150+_over_50-100_tcMet_inz_eemm-em";
  TH1F* hres_ratio3 = new TH1F( name_ratio3, name_ratio3, nsjpbins, 0, nsjpbins );
  hres_ratio3->Sumw2();

  //table of integral heading
  cout << "|*sjp bin*|*0-50*|*0-100*|*50-100*|*100-150*|*100+*|*150+*|*200+*|" << endl;

  //sigma vs sjp -- sigma from fit
  double* sigma1[nsjpbins];
  TH1F* hsigma1 = new TH1F( "Sigma_fit_full", "Sigma_fit_full", nsjpbins, 0, nsjpbins );
  hsigma1->Sumw2();
  double* sigma2[nsjpbins];
  TH1F* hsigma2 = new TH1F( "Sigma_fit_pm1sigma", "Sigma_fit_pm1sigma", nsjpbins, 0, nsjpbins );
  hsigma2->Sumw2();
  //double* sigma3[nsjpbins];
  TH1F* hsigma3 = new TH1F( "Sigma_fit_itr", "Sigma_fit_itr", nsjpbins, 0, nsjpbins );
  hsigma3->Sumw2();

  //stat error hists -- sigma from yield only, no fit
  TCanvas *c3 = new TCanvas();
  //TH1F* hsigma_stat = new TH1F( "Sigma_stat_tcMet", "Sigma_stat_tcMet",
  //TH1F* hsigma_stat = h1->getHistogramSum(theSources, "tcMet_inz", "", ee, mm);
  //TH1F* hsigma_stat->Add( h1->getHistogram(theSources, "tcMet_inz", "", em) );
  TH1F* hsigma_stat = h1->getHistogram(theSources, "tcMet_inz", "", all);
  hsigma_stat->Draw();
  c3->SaveAs( (TString)hsigma_stat->GetName()+"_hist.png" );
  HSqrt( hsigma_stat );
  c3->SaveAs( (TString)hsigma_stat->GetName()+"_sqrt.png" );

  TH1F* h_sigma_stat_den = h1->getHistogramSum(theSources, "tcMet_inz", "", ee, mm);
  h_sigma_stat_den->Add( h1->getHistogram(theSources, "tcMet_inz", "", em), -1. ); //subtract em
  hsigma_stat->Divide( h_sigma_stat_den );
  hsigma_stat->Draw();
  c3->SaveAs( "Sigma_stat_tcMet_inz.png" );

  //syst error hists -- sigma from yield only, no fit
  TH1F* hsigma_syst = h1->getHistogramSum(sources_nody, "tcMet_inz", "", ee, mm);
  hsigma_syst->Add( h1->getHistogram(sources_nody, "tcMet_inz", "", em), -1. );
  hsigma_syst->Divide( h_sigma_stat_den ); //stat and syst have same denominator
  hsigma_syst->Draw();
  c3->SaveAs( "Sigma_syst_tcMet_inz.png" );

  //stat, syst error from yield vs sjp
  //double sigmastat[nsjpbins];
  TH1F* hsigma_stat_vsjp = new TH1F( "Sigma_stat_tcMet_inz_vsjp", "Sigma_stat_tcMet_inz_vsjp", nsjpbins, 0, nsjpbins );
  hsigma_stat_vsjp->Sumw2();
  //double sigmasyst[nsjpbins];
  TH1F* hsigma_syst_vsjp = new TH1F( "Sigma_syst_tcMet_inz_vsjp", "Sigma_syst_tcMet_inz_vsjp", nsjpbins, 0, nsjpbins );
  hsigma_syst_vsjp->Sumw2();

  for(int i=0; i<nsjpbins; i++) {

	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_outz_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_inz_sjp", i), "", all, true);

	TString nametcsjp = Form("%s%i", "tcMet_xy_sjp", i);
	TString nmtcousjp = Form("%s%i", "tcMet_xy_outz_sjp", i);
	TString nmtcinsjp = Form("%s%i", "tcMet_xy_inz_sjp", i);

	TString namegnsjp = Form("%s%i", "genMet_xy_sjp", i);
	TString nmgnousjp = Form("%s%i", "genMet_xy_outz_sjp", i);
	TString nmgninsjp = Form("%s%i", "genMet_xy_inz_sjp", i);

	makeStack(h1, lg_all, theSources, nametcsjp, "", all, true);
	makeStack(h1, lg_all, theSources, nmtcousjp, "", all, true);
	makeStack(h1, lg_all, theSources, nmtcinsjp, "", all, true);

	makeStack(h1, lg_all, theSources, nametcsjp, "", em, true);
	makeStack(h1, lg_all, theSources, nmtcousjp, "", em, true);
	makeStack(h1, lg_all, theSources, nmtcinsjp, "", em, true);

	makeSumStack(h1, lg_all, theSources, nametcsjp, "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, nmtcousjp, "", ee, mm, true);
	makeSumStack(h1, lg_all, theSources, nmtcinsjp, "", ee, mm, true);

	makeSumStack(h1, lg_all, theSources, nametcsjp, "", all, all, true, namegnsjp, -1., "tcMet-" + namegnsjp);
	makeSumStack(h1, lg_all, theSources, nmtcousjp, "", all, all, true, nmgnousjp, -1., "tcMet-" + nmgnousjp);
	makeSumStack(h1, lg_all, theSources, nmtcinsjp, "", all, all, true, nmgninsjp, -1., "tcMet-" + nmgninsjp);

	makeSumDifStack(h1, lg_all, theSources, nametcsjp, "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, nmtcousjp, "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, nmtcinsjp, "", ee, mm, em, true);

	makeSumDifStack(h1, lg_all, sources_nody, nametcsjp, "", ee, mm, em, true, "nody_" + nametcsjp + "_eemm-em");
	makeSumDifStack(h1, lg_all, sources_nody, nmtcousjp, "", ee, mm, em, true, "nody_" + nmtcousjp + "_eemm-em");
	makeSumDifStack(h1, lg_all, sources_nody, nmtcinsjp, "", ee, mm, em, true, "nody_" + nmtcinsjp + "_eemm-em");

	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_sjp", i)     , "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_outz_sjp", i), "", ee, mm, em, true);
	makeSumDifStack(h1, lg_all, theSources, Form("%s%i", "tcMet_inz_sjp", i) , "", ee, mm, em, true);

	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_outz_sjp", i), "", all, true);
	makeStack(h1, lg_all, theSources, Form("%s%i", "dphi_tcMetl_inz_sjp", i), "", all, true);

	//for this one, the one without 'den' has weight as weight*tcmet, so don't look at, just use for numerator in divide
	makeStack(h1, lg_all, theSources, Form("%s%i", "tcMet_mllden_sjp", i), "", all, true);

	//integral ratio
	//TH1F* h_int_ratio = h1->getHistogram(theSources, Form("%s%i", "tcMet_sjp", i), "", all);
	//new ratio plots use inz, after subtraction only. since it's integral, just use hists, not stack
	//TH1F* h_int_ratio = h1->getHistogramSum(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", ee, mm);
	h_int_ratio[i] = h1->getHistogramSum(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", ee, mm);
	TH1F* h_int_ratio_em = h1->getHistogram(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", em);
	h_int_ratio[i]->Add( h_int_ratio_em, -1. ); //subtract em
	TCanvas *c2 = new TCanvas();
	//h_int_ratio->Draw();
	//c2->SaveAs( (TString)h_int_ratio->GetName()+".png" );
	//c2->SaveAs( (TString)Form("%s%i", "tcMet_inz_sjp", i, "_eemm-em_hist")+".png" );

	//res_ratio1[i] = 0; //initialize here
	res_ratio1[i] = getIntegralRatio( h_int_ratio[i], 50, 100, 200 );
	hres_ratio1->Fill( i, res_ratio1[i][0] );
	hres_ratio1->SetBinError( i, res_ratio1[i][1] ); //bin # = axis value

	//res_ratio2[i] = 0; //initialize here
	res_ratio2[i] = getIntegralRatio( h_int_ratio[i], 0, 50, 100 );
	hres_ratio2->Fill( i, res_ratio2[i][0] );
	hres_ratio2->SetBinError( i, res_ratio1[i][1] ); //bin # = axis value

	//res_ratio3[i] = 0; //initialize here
	res_ratio3[i] = getIntegralRatio( h_int_ratio[i], 50, 100, 150 );
	hres_ratio3->Fill( i, res_ratio3[i][0] );
	hres_ratio3->SetBinError( i, res_ratio1[i][1] ); //bin # = axis value

	//table of integrals
	//cout << "|*sjp bin*|0-50|0-100|50-100|100-150|100+|150+|200+|" << endl;
	cout << "|*" << i << "*|"
		 << getIntegral( h_int_ratio[i], 0, 50 ) << "|"
		 << getIntegral( h_int_ratio[i], 0, 100 ) << "|"
		 << getIntegral( h_int_ratio[i], 50, 100 ) << "|"
		 << getIntegral( h_int_ratio[i], 100, 150 ) << "|"
		 << getIntegral( h_int_ratio[i], 100 ) << "|"
		 << getIntegral( h_int_ratio[i], 150 ) << "|"
		 << getIntegral( h_int_ratio[i], 200 ) << "|"
		 << endl;

	//sigma fit hists
	TH1F* h_sigma = h1->getHistogramSum(theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", ee, mm);
	h_sigma->Add( h1->getHistogram(theSources, Form("%s%i", "tcMet_xy_inz_sjp", i), "", em), -1. ); //subtract em
	//TH1F* h_sigma = h2->getHistogram(theSources, Form("%s%i", "tcMet_xy_sjp", i), "", all);

	//sigma1[i] = 0; //initialize here
	sigma1[i] = getHistSigma( h_sigma );
	hsigma1->Fill( i, sigma1[i][0] );
	hsigma1->SetBinError( i, sigma1[i][1] );

	//sigma2[i] = 0; //initialize here
	sigma2[i] = getHistSigma( h_sigma, -sigma1[i][0], sigma1[i][0] ); //range of fit is +/- 1 sigma using sigma of full range
	hsigma2->Fill( i, sigma2[i][0] );
	hsigma2->SetBinError( i, sigma2[i][1] );

	//iterative fitting
	//double* sigma_itr[] = {sigma2[i][0],sigma2[i][1]}; //error is irrelevant...
	double* sigma_itr = sigma1[i];
	for( int j=0; j<2; j++ ) { //two iterations are enough
	  sigma_itr = getHistSigma( h_sigma, -sigma_itr[0], sigma_itr[0] ); //range of fit is +/- 1 sigma using sigma of full range
	  //cout << sigma_itr[0] << "  " ;
	}
	//cout << endl;

	//put iterative fitting results in sigma3
	//sigma3[i] = 0; //initialize here
	//sigma3[i] = getHistSigma( h_sigma, -2*sigma1[i][0], 2*sigma1[i][0] ); //range of fit is +/- 2 sigma using sigma of full range
	hsigma3->Fill( i, sigma_itr[0] );
	hsigma3->SetBinError( i, sigma_itr[1] );
	//cout << sigma1[i] << "  " << sigma2[i] << "  " << sigma3[i] << "  " << endl; //put in hist
	//cout << sigma1[i][0] << "  ";// << sigma2[i][0] << "  " ; //put in hist
	
	//sigma stat hists
	TH1F* hsigma_stat_sjp = h1->getHistogram(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", all);
	double stat_num = sqrt( hsigma_stat_sjp->Integral() );
	//TCanvas *c3 = new TCanvas();
	//hsigma_stat->Draw();
	//c3->SaveAs( (TString)hsigma_stat->GetName()+"_hist.png" );
	HSqrt( hsigma_stat_sjp );
	//c3->SaveAs( (TString)hsigma_stat->GetName()+"_sqrt.png" );

	TH1F* hsigma_stat_sjp_den = h1->getHistogramSum(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", ee, mm);
	hsigma_stat_sjp_den->Add( h1->getHistogram(theSources, Form("%s%i", "tcMet_inz_sjp", i), "", em), -1. ); //subtract em
	double stat_den = hsigma_stat_sjp_den->Integral();
	hsigma_stat_sjp->Divide( hsigma_stat_sjp_den );
	hsigma_stat_sjp->Draw("hist"); //"hist" option turns off error bars (because that's the most sensible thing to call it)
	c2->SaveAs( (TString)Form("%s%i", "Sigma_stat_tcMet_inz_sjp", i)+".png" ); 

	//syst error hists -- sigma from yield only, no fit
	TH1F* hsigma_syst_sjp = h1->getHistogramSum(sources_nody, Form("%s%i", "tcMet_inz_sjp", i), "", ee, mm);
	hsigma_syst_sjp->Add( h1->getHistogram(sources_nody, Form("%s%i", "tcMet_inz_sjp", i), "", em), -1. );
	double syst_num = hsigma_syst_sjp->Integral();
	hsigma_syst_sjp->Divide( hsigma_stat_sjp_den ); //stat and syst have same denominator
	hsigma_syst_sjp->Draw();
	c2->SaveAs( (TString)Form("%s%i", "Sigma_syst_tcMet_inz_sjp", i)+".png" ); 

	//sigma stat, syst from yield vs sjp
	// f = x/y, err_f/f = sqrt( err_x^2/x^2 + err_y^2/y^2 ) -->> err_f^2 = x/y^2 + x^2/y^3 = (f/y)(1+f)
	double ratio = stat_num/stat_den;
	hsigma_stat_vsjp->Fill( i, ratio );
	hsigma_stat_vsjp->SetBinError( i, sqrt( (ratio/stat_den)*(1+ratio) ) );
	ratio = syst_num/stat_den;
	hsigma_syst_vsjp->Fill( i, ratio );
	hsigma_syst_vsjp->SetBinError( i, sqrt( (ratio/stat_den)*(1+ratio) ) );
  }


  //save integral ratio
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptStat(0);

  h_int_ratio[0]->GetYaxis()->SetRangeUser( -1.0, 1.0 );
  for( int i=0; i<nsjpbins; i++ ) {
	h_int_ratio[i]->SetLineColor( i );
	h_int_ratio[i]->Draw("same");
  }
  //c1->SaveAs( (TString)Form("%s%i", "tcMet_inz_sjp", i, "_eemm-em_hist")+".png" );
  //gPad->SetLogy();
  c1->SaveAs( "tcMet_inz_sjpall_eemm-em.png" );
  //gPad->SetLogy(0);
  hres_ratio1->Draw();
  c1->SaveAs((TString)hres_ratio1->GetName()+".png");
  hres_ratio2->Draw();
  c1->SaveAs((TString)hres_ratio2->GetName()+".png");
  hres_ratio3->Draw();
  c1->SaveAs((TString)hres_ratio3->GetName()+".png");
  gStyle->SetOptStat(1110111);

  //save sigma_fit
  hsigma1->Draw();
  c1->SaveAs((TString)hsigma1->GetName()+".png");
  hsigma2->Draw();
  c1->SaveAs((TString)hsigma2->GetName()+".png");
  hsigma3->Draw();
  c1->SaveAs((TString)hsigma3->GetName()+".png");
  
  //save sigma stat, syst vs sjp
  hsigma_stat_vsjp->Draw();
  c1->SaveAs((TString)hsigma_stat_vsjp->GetName()+".png");
  hsigma_syst_vsjp->Draw();
  c1->SaveAs((TString)hsigma_syst_vsjp->GetName()+".png");


  //do 2d hist
  TH2F* metxy = h1->get2dHistogram(theSources, "tcMet_xvy", "", all);
  //TString filename = "Results_plot.root";
  //TFile outf(filename,"RECREATE") ;
  metxy->SetName("tcMet_xvy");
  metxy->Draw("box");
  c1->SaveAs((TString)metxy->GetName()+".png");

  //do sum hists

  //subtraction check
  THStack *tcMet_xy_nody = h1->getSumStack(sources_nody, "tcMet_xy", "", ee, mm);
  TH1F* tcMet_xy_em = h1->getHistogram(theSources, "tcMet_xy", "", em);
  tcMet_xy_nody->Draw("hist");
  tcMet_xy_em->Draw("same");
  c1->SaveAs("tcMet_xy_sf_nody_of_comp.png");
  gPad->SetLogy();
  c1->SaveAs("tcMet_xy_sf_nody_of_comp_log.png");
  gPad->SetLogy(0);

  //second subtraction check
  THStack *tcMet_xy_subpos = h1->getSumDifStack(theSources, "tcMet_xy", "", ee, mm, em, 1, 1);
  THStack *tcMet_xy_subneg = h1->getSumDifStack(theSources, "tcMet_xy", "", ee, mm, em, 1, -1);
  TH1F* tcMet_xy_dy = h1->getHistogram(sources_dy, "tcMet_xy", "", all);
  tcMet_xy_dy->SetLineColor( 1 ); //0=white, 1=black
  tcMet_xy_dy->SetMarkerColor( 1 );
  tcMet_xy_subpos->SetMinimum( tcMet_xy_subneg->GetMinimum() );
  tcMet_xy_subpos->Draw("hist");
  tcMet_xy_subneg->Draw("same,hist");
  tcMet_xy_dy->Draw("same");
  c1->SaveAs("tcMet_xy_eemm-em_dy_comp.png");
  gPad->SetLogy();
  c1->SaveAs("tcMet_xy_eemm-em_dy_comp_log.png");
  gPad->SetLogy(0);

  //tcMet_x v tcMet_y
  //bools are copylog, copyscale, copydiff
  makeHistOverlay( h1, theSources, theSources, "tcMet_x", "tcMet_y", "", all, all, true, false, true);
  //makeHistOverlay( h1, theSources, theSources, "tcMet_x_nar", "tcMet_y_nar", "", all, all, true, false, true); //replaced with below
  makeHistOverlay( h2, theSources, theSources, "tcMet_x", "tcMet_y", "", all, all, true, false, true,"", "_nar"); //add suffix for h2 
  // * / 
  //xy_in vs xy_out
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true );
  makeHistOverlay( h1, theSources, theSources, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true );
  sources_t sources_dyee = (1ll << H_DYEE);
  sources_t sources_dymm = (1ll << H_DYMM);
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dyee_");
  makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dyee_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dymm_");
  makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dymm_");
  //same in narrow
  makeHistOverlay( h2, sources_dyee, sources_dyee, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dyee_", "_nar");
  makeHistOverlay( h2, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dyee_", "_nar");
  makeHistOverlay( h2, sources_dymm, sources_dymm, "tcMet_xy_inz", "tcMet_xy_outz", "", all, all, true, true, true, "dymm_", "_nar");
  makeHistOverlay( h2, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, "dymm_", "_nar");

  //these are wrong--ignore for now
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dyee, sources_dyee, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dyee_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_nar", "tcMet_xy_outz_nar", "", all, all, true, true, true, (TString)"dymm_");
  //makeHistOverlay( h1, sources_dymm, sources_dymm, "tcMet_xy_inz_sjp0", "tcMet_xy_outz_sjp0", "", all, all, true, true, true, (TString)"dymm_");


  //fitting -- use h2 b'c finer binning

  //tcMet_xy
  TH1F* tcMet_xy_all = h2->getHistogram(theSources, "tcMet_xy", "", all);
  double* sigma_xy_all = getHistSigma( tcMet_xy_all );
  double trato_xy_all = getTailRatio( tcMet_xy_all );
  cout << "name " << tcMet_xy_all->GetName() << "  sigma " << sigma_xy_all[0] << "  ratio core/tail " << trato_xy_all << endl;

  //tcMet_x
  TH1F* tcMet_x_all = h2->getHistogram(theSources, "tcMet_x", "", all);
  double* sigma_x_all = getHistSigma( tcMet_x_all );
  double trato_x_all = getTailRatio( tcMet_x_all );
  cout << "name " << tcMet_x_all->GetName() << "  sigma " << sigma_x_all[0] << "  ratio core/tail " << trato_x_all << endl;

  //tcMet_y
  TH1F* tcMet_y_all = h2->getHistogram(theSources, "tcMet_y", "", all);
  double* sigma_y_all = getHistSigma( tcMet_y_all );
  double trato_y_all = getTailRatio( tcMet_y_all );
  cout << "name " << tcMet_y_all->GetName() << "  sigma " << sigma_y_all[0] << "  ratio core/tail " << trato_y_all << endl;

  //for dyee+mm
  TH1F* dy_tcMet_xy_all = h2->getHistogram(sources_dy, "tcMet_xy", "", all);
  double* dy_sigma_xy_all = getHistSigma( dy_tcMet_xy_all );
  double dy_trato_xy_all = getTailRatio( dy_tcMet_xy_all );
  cout << "name " << dy_tcMet_xy_all->GetName() << "  sigma " << dy_sigma_xy_all[0] << "  ratio core/tail " << dy_trato_xy_all << endl;

  cout << endl << endl;

  //twiki table
  cout << "Fit results for met xy distribution\n";
  cout << "||*DYee+DYmm*|*All Samples*|\n";
  cout << "|*sigma*|" << dy_sigma_xy_all[0] << "|" << sigma_xy_all[0] << "|\n";
  cout << "|*ratio*|" << dy_trato_xy_all << "|" << trato_xy_all << "|\n";
  //add note on twiki on definition of all samples, and ratio

}

