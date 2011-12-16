{


  gSystem->Load("libFWCoreFWLite");

  cout << "Loading CMS TDR style..." << flush;
  #include "/afs/cern.ch/user/d/dalfonso/scratch0/UserCode/dalfonso/Utils/setTDRStyle.C"
//  #include "setTDRStyle.C"
  setTDRStyle();
  tdrStyle->SetOptStat(0000000);
///  tdrStyle->SetOptStat(111111111);
  tdrStyle->SetTitleSize(0.05, "XYZ");
//  tdrStyle->SetTitleSize(0.1, "XYZ");
//  tdrStyle->SetOptTitle(0);
  tdrStyle->SetPadBottomMargin(0.14);
//  tdrStyle->SetPadTopMargin(0.2);
  tdrStyle->SetPadTopMargin(0.1);
  tdrStyle->SetPadLeftMargin(0.18);
//  tdrStyle->SetPadLeftMargin(0.1);
//  tdrStyle->SetPadRightMargin(0.1);
  tdrStyle->SetPadRightMargin(0.2);
//  tdrStyle->SetTitleXOffset(1.0);
  tdrStyle->SetTitleXOffset(0.6);
  tdrStyle->SetTitleYOffset(1.);
//  tdrStyle->SetTitleOffset(0.1,"Z");
//  tdrStyle->SetTitleYOffset(1.);
  tdrStyle->SetNdivisions(505, "X");
  tdrStyle->SetErrorX(0.5);
  tdrStyle->SetPalette(1,0);
  gROOT->ForceStyle();
  cout << " done." << endl;

  setTDRStyle();

//   cout << "Loading Parang, a plot-making suite..." << flush;
//   if (gSystem->Load("libUtilitiesParang") == 0) {
//     gSystem->AddIncludePath(" -I$CMSSW_BASE/src ");
//     cout << " OK. " << endl;
//   }
//   else  cout << " Uh oh! " << endl;

//    AutoLibraryLoader::enable();

}
