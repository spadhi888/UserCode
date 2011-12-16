#include "Riostream.h"
void upperlimit() {
   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("upperlimit.C","");
   dir.ReplaceAll("/./","/");
   ifstream in1;
   ifstream in2;
   ifstream in3;
   ifstream in4;



   Float_t x1, x2, x3, x4, x5, x6;
   Int_t nlines = 0;
   double lumi = 980.;
   //   TFile *f = new TFile("RA5SSlimitMatched.root","RECREATE");
   TFile *f = new TFile("RA5SSlimitNOTMatched.root","RECREATE");
   TH2F *h2sig1UL = new TH2F("h2sig1UL","h2sig1UL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig2UL = new TH2F("h2sig2UL","h2sig2UL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig3UL = new TH2F("h2sig3UL","h2sig3UL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig4UL = new TH2F("h2sig4UL","h2sig4UL",60, 0.,1500., 40, 0. , 1000.);

   TH2F *h2sig1expectedUL = new TH2F("h2sig1expectedUL","h2sig1expectedUL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig2expectedUL = new TH2F("h2sig2expectedUL","h2sig2expectedUL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig3expectedUL = new TH2F("h2sig3expectedUL","h2sig3expectedUL",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig4expectedUL = new TH2F("h2sig4expectedUL","h2sig4expectedUL",60, 0.,1500., 40, 0. , 1000.);

   TH2F *h2sig1Eff = new TH2F("h2sig1Eff","h2sig1Eff",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig2Eff = new TH2F("h2sig2Eff","h2sig2Eff",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig3Eff = new TH2F("h2sig3Eff","h2sig3Eff",60, 0.,1500., 40, 0. , 1000.);
   TH2F *h2sig4Eff = new TH2F("h2sig4Eff","h2sig4Eff",60, 0.,1500., 40, 0. , 1000.);

   //   in1.open(Form("upperlimit_sig1_matched.txt",dir.Data()));
   in1.open(Form("upperlimit_sig1_NOTmatched.txt",dir.Data()));
   while (1) {
     in1 >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
      if (!in1.good()) break;
      if (nlines < 5) printf("x1=%8f, x2=%8f, x3=%8f, x4=%8f, x5=%8f, x6=%8f\n",x1,x2,x3,x4,x5,x6);
      //      cout << x <<  "  " << y << "  " << z << endl; 
      h2sig1UL->Fill(x1, x2, x5/(lumi*x3));
      h2sig1expectedUL->Fill(x1, x2, x6/(lumi*x3));
      h2sig1Eff->Fill(x1, x2, x3);
      nlines++;
   }
   printf(" found %d points\n",nlines);
   in1.close();


   nlines = 0;
   //   in2.open(Form("upperlimit_sig2_matched.txt",dir.Data()));
   in2.open(Form("upperlimit_sig2_NOTmatched.txt",dir.Data()));
   while (1) {
     in2 >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
      if (!in2.good()) break;
      if (nlines < 5) printf("x1=%8f, x2=%8f, x3=%8f, x4=%8f, x5=%8f, x6=%8f\n",x1,x2,x3,x4,x5,x6);
      //      cout << x <<  "  " << y << "  " << z << endl; 
      h2sig2UL->Fill(x1, x2, x5/(lumi*x3));
      h2sig2expectedUL->Fill(x1, x2, x6/(lumi*x3));
      h2sig2Eff->Fill(x1, x2, x3);
      nlines++;
   }
   printf(" found %d points\n",nlines);
   in2.close();


   nlines = 0;
   //   in3.open(Form("upperlimit_sig3_matched.txt",dir.Data()));
   in3.open(Form("upperlimit_sig3_NOTmatched.txt",dir.Data()));
   while (1) {
     in3 >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
      if (!in3.good()) break;
      if (nlines < 5) printf("x1=%8f, x2=%8f, x3=%8f, x4=%8f, x5=%8f, x6=%8f\n",x1,x2,x3,x4,x5,x6);
      //      cout << x <<  "  " << y << "  " << z << endl; 
      h2sig3UL->Fill(x1, x2, x5/(lumi*x3));
      h2sig3expectedUL->Fill(x1, x2, x6/(lumi*x3));
      h2sig3Eff->Fill(x1, x2, x3);
      nlines++;
   }
   printf(" found %d points\n",nlines);
   in3.close();

   nlines = 0;
   //   in4.open(Form("upperlimit_sig3_matched.txt",dir.Data()));
   in4.open(Form("upperlimit_sig4_NOTmatched.txt",dir.Data()));
   while (1) {
     in4 >> x1 >> x2 >> x3 >> x4 >> x5 >> x6;
      if (!in4.good()) break;
      if (nlines < 5) printf("x1=%8f, x2=%8f, x3=%8f, x4=%8f, x5=%8f, x6=%8f\n",x1,x2,x3,x4,x5,x6);
      //      cout << x <<  "  " << y << "  " << z << endl; 
      h2sig4UL->Fill(x1, x2, x5/(lumi*x3));
      h2sig4expectedUL->Fill(x1, x2, x6/(lumi*x3));
      h2sig4Eff->Fill(x1, x2, x3);
      nlines++;
   }
   printf(" found %d points\n",nlines);
   in4.close();

   f->Write();

}
