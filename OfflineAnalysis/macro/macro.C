#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include <fstream>
#include "TLegendEntry.h"
#include "TMath.h"


using namespace std;


void macro() {

  TFile* thefile=new TFile("processed_data_tag_scanNEW.root" ,"read");
  TH2F *all  = (TH2F*)gDirectory->Get("susyscan_hm12m0_all");

 
  double prob = 0.;
  double probm0 = 0.;
  double probm1 = 0.;
  
  double confidence0 = 3.2;
  double confidence1 = 4.7;


  TH2F* hex0 = new TH2F("hex0","hex0",80, 0.,4000., 26, 100. , 620.);
  TH2F* hex1 = new TH2F("hex1","hex1",80, 0.,4000., 26, 100. , 620.);
  TH2F* hexx =  new TH2F("hexx","hexx",80, 0.,4000., 26, 100. , 620.);


  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(2);

  for (int m0=0; m0 < 1450; m0=m0+50) {
    double max0 = 0;
    double max1 = 0;
    probm0 = 0;
    probm1 = 0;
    for (int m12=100; m12 < 620; m12=m12+20) {
      prob       = all->GetBinContent(all->FindBin(m0,m12) );
      if (prob > 0) hexx->Fill(m0,m12,prob);
      if (prob > confidence0) {
        hex0->Fill(m0,m12,prob);
        if (m12 > max0) {
	  max0 = m12;
	  probm0 = prob;
	}
      } 
      if (prob > confidence1) {
	hex1->Fill(m0,m12,prob);
        if (m12 > max1) {
	  max1 = m12;
	  probm1 = prob;
	}
      }
      
    }

    if (max0 > 0 ) cout << "| Exclude_0 | " << m0 << " |  " << max0 << " | " << probm0 << " | " << endl;

    if (max1 > 0 ) cout << "| Exclude_1 | " << m0 << " |  " << max1 << " | " << probm1 << " | " <<  endl;


  } 


}
