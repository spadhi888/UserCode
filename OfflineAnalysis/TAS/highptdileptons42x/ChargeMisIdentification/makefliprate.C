{
	
  TCanvas *canvas = new TCanvas;
  TPaveText *box = 0;

  //TH1F *els_MT_noflip = (TH1F*)gFile->Get("els_Mt_noflip");
  //els_MT_noflip->Drax();
  TH2F *els_2d_eta_Pt_corCharge= (TH2F*)gFile->Get("els_2d_eta_Pt_corCharge");
  TH2F *els_2d_eta_Pt_incorCharge= (TH2F*)gFile->Get("els_2d_eta_Pt_incorCharge");


  //float etabins[] = {0, 1.28, 1.56, 1.84, 2.12, 2.5};
  //float etabins[] = {0, 0.5, 1.0, 1.28, 1.56, 1.84, 2.12, 2.4};
  //float ptbins[] = {10, 30, 40, 50, 70, 200};

//  float etabins[] = {0, 0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4};
//  float ptbins[] = {10, 20, 30, 40, 50, 60, 70, 80, 100};
//  float etabins[] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
//  float ptbins[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};

//  float etabins[] = {0, 0.5, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7,
//                     1.8, 1.9, 2.0, 2.2, 2.5};
//  float ptbins[] = {10, 30, 40, 50, 60, 70, 80, 100};


// - flip rate
//  float ptbins[] = {10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
//  float etabins[] = {0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};

  float ptbins[] = {10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100};
  float etabins[] = {0, 0.5, 1.0, 1.479, 1.8, 2.0, 2.1, 2.2, 2.5};


  unsigned etanbins = sizeof(etabins)/sizeof(etabins[0])-1;
  unsigned ptnbins = sizeof(ptbins)/sizeof(ptbins[0])-1; 

  // open a file
  ofstream outfile;
  outfile.open ("fliprate_egun.cc");


  // write one file 
  outfile << "//-------------------------------------------------------------------" << endl;
  outfile << "// Returns the sign flip rate from the electron gun" << endl;
  outfile << "// Note that the electron gun only goes up to 100 GeV" << endl;
  outfile << "// For electron pt>100 returns the flip rate in the 90-100 Gev bin" << endl;
  outfile << "// For electron pt<10  returns zero" << endl;
  outfile << "// For electron abs(eta)>2.4 returns zero & prints out an error" << endl;
  outfile << "//" << endl;
  outfile << "// Usage:" << endl;
  outfile << "// double flipRate   = getSingleEleFlipRate(el_pt, el_eta)" << endl;
  outfile << "// double fliRateErr = getSingleEleFlipRateError(el_pt, el_eta)" << endl;
  outfile << "//" << endl;
  outfile << "// Claudio & Derek 23 July 2009" << endl;
  outfile << "//-------------------------------------------------------------" << endl;
  outfile << "#include \"fliprate_egun.h\"" << endl;
  outfile << "#include <iostream>" << endl;
  outfile << "#include <stdio.h>" << endl;
  outfile << "#include <math.h>" << endl;
  outfile << endl;
  outfile << "using namespace std;" << endl;
  outfile << endl;

  outfile << "double getSingleEleNum(double el_pt, double el_eta) {" << endl;
  outfile << endl;

  outfile << "  el_eta = fabs(el_eta);" << endl;
  outfile << endl;


  // 2d ratio
  Float_t deno=0, num=0;
  //for(Int_t i=0; i<7; i++)
  for(Int_t i=1; i<=etanbins; i++)
  {
    outfile << "  if( el_eta < " ; 
    if(i==etanbins) outfile << "= 2.4 ){ \n" ;
    else outfile << etabins[i] << " ){ \n" ;

    //for(Int_t j=9; j>0; j--)
    for(Int_t j=ptnbins; j>0; j--)
    {
       //deno = els_2d_eta_Pt_corCharge->GetBinContent(i, j)	
       //     +els_2d_eta_Pt_incorCharge->GetBinContent(i, j);
       num = els_2d_eta_Pt_incorCharge->GetBinContent(i, j);
       outfile << "    if( el_pt > " ; 
       if(j==1) outfile << ptbins[j-1] << " ) return " << num << ";\n" 
                        << "    return 0.0;\n"<< "  }" << endl;
       else outfile << ptbins[j-1] << " ) return " << num << ";" << endl;

       //els_2d_eta_Pt_ratio->SetBinContent(i+1, j+1, num/deno);
    }
  } 

  outfile << "  std::cout << \"Error: eta > 2.4 value found\" << endl;" << endl; 
  outfile << "  return 0.0;" << endl;
  outfile << "}" << endl;
  outfile << endl;
  outfile << endl;

  outfile << "double getSingleEleDenom(double el_pt, double el_eta) {" << endl;
  outfile << endl;
  outfile << "  el_eta = std::fabs(el_eta);" <<endl;
  outfile << endl;


  deno=0, num=0;
  //for(Int_t i=0; i<7; i++)
  for(Int_t i=1; i<=etanbins; i++)
  {
    outfile << "  if( el_eta < " ; 
    if(i==etanbins) outfile << "= 2.4 ){ \n" ;
    else outfile << etabins[i] << " ){ \n" ;

    //for(Int_t j=9; j>0; j--)
    for(Int_t j=ptnbins; j>0; j--)
    {
       deno = els_2d_eta_Pt_corCharge->GetBinContent(i, j)	
            +els_2d_eta_Pt_incorCharge->GetBinContent(i, j);
       //num = els_2d_eta_Pt_incorCharge->GetBinContent(i, j);
       outfile << "    if( el_pt > " ; 
       if(j==1) outfile << ptbins[j-1] << " ) return " << deno << ";\n" 
                        << "    return 0.0;\n"<< "  }" << endl;
       else outfile << ptbins[j-1] << " ) return " << deno << ";" << endl;

       //els_2d_eta_Pt_ratio->SetBinContent(i+1, j+1, num/deno);
    }
  }

  outfile << "  std::cout << \"Error: eta > 2.5 value found\" << endl;" << endl;
  outfile << "  return 0.0;" << endl;
  outfile << "}" << endl;
  outfile << endl;
  outfile << endl;

  outfile << "double getSingleEleFlipRate(double el_pt, double el_eta) {" << endl;
  outfile << "  if( el_pt < 10.0 || fabs(el_eta) > 2.5 ){" << endl;
  outfile << "    std::cout << \"Error in 'getSingleEleFlipRate': pt or eta value found out of range\" << endl;" << endl; 
  outfile << "    return 0.0;" << endl;
  outfile << "  }" << endl;
  outfile << "  return getSingleEleNum(el_pt, fabs(el_eta))/getSingleEleDenom(el_pt, fabs(el_eta));" << endl;
  outfile << "}" << endl;
  outfile << endl;

  outfile << "double getSingleEleFlipRateError(double el_pt, double el_eta) {" << endl;
  outfile << "  //the binomial error" << endl;
  outfile << "  if( el_pt < 10.0 || fabs(el_eta) > 2.5 ){" << endl;
  outfile << "    std::cout << \"Error in 'getSingleEleFlipRate': pt or eta value found out of range\" << endl;" << endl;
  outfile << "    return 0.0;" << endl;
  outfile << "  }" << endl;
  outfile << "  double num   = getSingleEleNum(el_pt,   fabs(el_eta));" << endl;
  outfile << "  double denom = getSingleEleDenom(el_pt, fabs(el_eta));" << endl;
  outfile << "  double p = num/denom;" << endl;
  outfile << "  return sqrt(p*(1-p)/denom);" << endl;
  outfile << "}" << endl;

  outfile.close();



  //els_2d_eta_Pt_corCharge->Draw();
  //els_2d_eta_Pt_incorCharge->Draw();

}
