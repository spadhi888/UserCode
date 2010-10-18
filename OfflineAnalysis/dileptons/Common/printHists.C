#include "THStack.h"
//#include "specFun.C"
#include <iostream>
#include "TMath.h"
#include "histtools.C"
#include "TDirectory.h"
#include "TROOT.h"
using namespace std;

/*
void getSoverRootN(double& rat, double& ratE, double s, double n, double sE, double nE){
  rat = s; rat /= n > 0 ? sqrt(n) : 1.;
  ratE = ( sE*sE* ( 1. - 0.5*rat)*( 1. - 0.5*rat)  + (nE*nE - sE*sE)*0.25*rat*rat);
  ratE /= n > 0 ? n : 1.;
  ratE = sqrt(ratE);
}

double bOnlyProb(double s, double b, double bE){
  //at some point need to put some protections here or find an appropriate code
  unsigned int lowExp = floor(s+b);
  double pSum = 0;
  for (int i=lowExp; i>=0; --i){
    if (b> 0. && bE/b>0.03 && s>0.1*b && bE<0.5*s){
      pSum+= poisson_smeared_prob(i,b, bE);
    } else {//use regular Poisson here
      pSum+= TMath::Poisson(i,b);
    }
  }
  return pSum;
}
*/
std::string formatFloat(double x, const char* formatS){
  std::string xS = Form(Form("%s", formatS),x);
  double xB = atof(xS.c_str());
  if (x>0 && xB==0){
    xS = Form(" %6.1g",x);
  }
  return xS;
}

void printNJets( bool latex=false, const char* formatS = "%6.1f", const char* signalS= "ttdil", 
		 bool combineVVsamples = true, bool combineDYsamples = true, bool combineJetBins = false, bool printProbs=false, bool printErrorsForData = false){
  char* suffix[4];
  suffix[0] = "ee";
  suffix[1] = "mm";
  suffix[2] = "em";
  suffix[3] = "all";


  std::string pmSign  = latex ? " \\pm " : " &plusmn; ";
  std::string colSep  = latex ? " & " : " | ";
  std::string beginL  = latex ? ""   : " | ";
  std::string endL    = latex ? " \\\\ " : " | ";
  std::string mathSep = latex ? "$" : "";


  if (latex) {
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\begin{center}" << std::endl;
    std::cout << "{\\small" << std::endl;
  }

  bool haveVVsamples = false;
  if(gROOT->FindObjectAny("ww_hnJet_ee") != NULL || 
     gROOT->FindObjectAny("zz_hnJet_ee") != NULL || 
     gROOT->FindObjectAny("wz_hnJet_ee") != NULL)
    haveVVsamples = true;

  bool haveData = false;
  if(gROOT->FindObjectAny("data_hnJet_ee") != NULL)
    haveData = true;
  
  for (int sample=0; sample<4; sample++) {
    
    hist::stack(Form("st_hnJet_%s", suffix[sample]), Form("^[^_]+_hnJet_%s$", suffix[sample]));
    THStack* thisStack = (THStack*) 
      gROOT->FindObjectAny(Form("st_hnJet_%s", suffix[sample]));
    //    std::cout<<"Found "<<thisStack->GetName()<<std::endl;

    int nHists = thisStack->GetHists()->GetSize();
    if (latex) {
      std::cout << "\\begin{minipage}{0.48\\textwidth}" << std::endl;
      if (sample == 0 || sample == 2){
	std::cout << "\\begin{tabular}{|l|c|c|c|} \\hline" << std::endl;
	std::cout << "\\hline \\hline" << std::endl;
	std::cout << "\\multicolumn{4}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
	std::cout << "   & $N_{jets}=0$       & $N_{jets}=1$      & $N_{jets} \\geq 2$ \\"<<"\\ \\hline" <<std::endl;
      } else {
	std::cout << "\\begin{tabular}{||c|c|c|} \\hline" << std::endl;
	std::cout << "\\hline \\hline" << std::endl;
	std::cout << "\\multicolumn{3}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
	std::cout << "   $N_{jets}=0$       & $N_{jets}=1$      & $N_{jets} \\geq 2$ \\"<<"\\ \\hline" <<std::endl;
      }
    } else {      
      if(!combineJetBins) {
	std::cout<<"====================================================";
	std::cout<<"Table for "<< suffix[sample]<<std::endl;
	std::cout<<" |   *sample* |       *nJet = 0*       |      *nJet = 1*        |      *nJet >= 2*      |"<<std::endl;
      } else {
	std::cout<<"===============================";
	std::cout<<"Table for "<< suffix[sample]<<std::endl;
	std::cout<<" |   *sample* |       *Yield*       |" << std::endl;
      }
    }
    //CHANGE LATER
    bool showFirstCol = true;
  
    //print signal first
    for(int iH=0; iH< nHists; ++iH){
      
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TString process(h1F->GetName());
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
      bool isSig = std::string(sampleNs->At(0)->GetName()) == std::string(signalS);
      
      if(!isSig)
	continue;
      
      double n0 = h1F->GetBinContent(1);
      double n0E = h1F->GetBinError(1);
      
      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);     

      double n2 = 0;
      double n2E = 0;
      
      for(unsigned int i = 3; i <= h1F->GetNbinsX() + 1; i++) {
	n2 = n2 + h1F->GetBinContent(i);
	double temp = h1F->GetBinError(i);
	n2E = sqrt(n2E*n2E + temp*temp);
      }

      
      std::string firstCol= Form("%9s ",sampleNs->At(0)->GetName());
      std::cout << beginL;
      if (showFirstCol) std::cout<< firstCol << colSep;
      if(!combineJetBins) {
	std::cout << mathSep << formatFloat(n0,formatS) <<pmSign<< formatFloat(n0E,formatS)<<mathSep<<colSep
		  << mathSep << formatFloat(n1,formatS) <<pmSign<< formatFloat(n1E,formatS)<<mathSep<<colSep
		  << mathSep << formatFloat(n2,formatS) <<pmSign<< formatFloat(n2E,formatS)<<mathSep
		  << endL
		  << std::endl;
      } else {
	double n = n0 + n1 + n2;
	double nE = sqrt(n0E*n0E + n1E*n1E + n2E*n2E);
	std::cout << mathSep << formatFloat(n,formatS) << pmSign << formatFloat(nE,formatS) << mathSep 
		  << endL << std::endl;
      }
    }//signal only


    //when we want to print the MC in the end
    double n0MCbkg = 0;
    double n0MCbkgE = 0;
    double n1MCbkg = 0;
    double n1MCbkgE = 0;
    double n2MCbkg = 0;
    double n2MCbkgE = 0;    
    
    //if you want to combine the VV samples
    vector<double> v_VV(3,0.0);
    vector<double> v_VVE(3,0.0); 
    //if you want to combine the DY samples
    vector<double> v_DY(3,0.0);
    vector<double> v_DYE(3,0.0);

    //combine the DY and VV samples
    bool foundVVsamples = false;
    bool foundDYsamples = false;
    for(int iH=0; iH<nHists; ++iH){
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TString process(h1F->GetName());
      if(combineVVsamples) {
	if(process.Contains("ww") || 
	   process.Contains("wz") ||
	   process.Contains("zz") ) {
	  foundVVsamples = true;
	  for(unsigned int j = 0; j < 2; j++) {
	    v_VV.at(j)  = v_VV.at(j)+h1F->GetBinContent(j+1);
	    v_VVE.at(j) = sqrt(pow(v_VVE.at(j),2)+pow(h1F->GetBinError(j+1),2));
	  }
	  for(unsigned int j = 2; j <= h1F->GetNbinsX(); j++) {
	    v_VV.at(2)  = v_VV.at(2)+h1F->GetBinContent(j+1);
	    v_VVE.at(2) = sqrt(pow(v_VVE.at(2),2)+pow(h1F->GetBinError(j+1),2));
	  }
	  continue;
	}
      }
      if(combineDYsamples) {
	if(process.Contains("DY")) {
	  foundDYsamples = true;
	  for(unsigned int j = 0; j < 2; j++) {
	    v_DY.at(j)  = v_DY.at(j)+h1F->GetBinContent(j+1);
	    v_DYE.at(j) = sqrt(pow(v_DYE.at(j),2)+pow(h1F->GetBinError(j+1),2));
	  }
	  for(unsigned int j = 2; j <= h1F->GetNbinsX(); j++) {
	    v_DY.at(2)  = v_DY.at(2)+h1F->GetBinContent(j+1);
	    v_DYE.at(2) = sqrt(pow(v_DYE.at(2),2)+pow(h1F->GetBinError(j+1),2));
	  }
	  continue;
	}
      }//combineDYsamples
    }
    if(combineVVsamples && foundVVsamples) {
      std::cout << beginL;
      if (showFirstCol) std::cout<<  Form("%9s ","VV") << colSep;
      if(!combineJetBins) {
      std::cout << mathSep << formatFloat(v_VV[0],formatS) <<pmSign<< formatFloat(v_VVE[0],formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(v_VV[1],formatS) <<pmSign<< formatFloat(v_VVE[1],formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(v_VV[2],formatS) <<pmSign<< formatFloat(v_VVE[2],formatS)<<mathSep
		<< endL
		<< std::endl;
      } else {   
	double n = v_VV[0] + v_VV[1] + v_VV[2];
	double nE = sqrt(pow(v_VVE[0],2) + pow(v_VVE[1], 2) + pow(v_VVE[2], 2));
	std::cout << mathSep << formatFloat(n,formatS) << pmSign << formatFloat(nE,formatS) << mathSep 
		  << endL << std::endl;
      }
    }
    if(combineDYsamples && foundDYsamples) {
       std::cout << beginL;
       if (showFirstCol) std::cout<<  Form("%9s ","DY") << colSep;
       if(!combineJetBins) {
       std::cout << mathSep << formatFloat(v_DY[0],formatS) <<pmSign<< formatFloat(v_DYE[0],formatS)<<mathSep<<colSep
		 << mathSep << formatFloat(v_DY[1],formatS) <<pmSign<< formatFloat(v_DYE[1],formatS)<<mathSep<<colSep
		 << mathSep << formatFloat(v_DY[2],formatS) <<pmSign<< formatFloat(v_DYE[2],formatS)<<mathSep
		 << endL
		 << std::endl;
       } else {
	 double n = v_DY[0] + v_DY[1] + v_DY[2];
	 double nE = sqrt(pow(v_DYE[0],2) + pow(v_DYE[1], 2) + pow(v_DYE[2], 2));
	 std::cout << mathSep << formatFloat(n,formatS) << pmSign << formatFloat(nE,formatS) << mathSep 
		   << endL << std::endl;
       }
    }
    

    //need to add the total DY, VV if we decided to combine the samples
    n0MCbkg  = n0MCbkg + v_VV[0];
    n0MCbkgE = sqrt(pow(n0MCbkgE,2) + pow(v_VVE[0],2));
    n1MCbkg  = n1MCbkg + v_VV[1];
    n1MCbkgE = sqrt(pow(n1MCbkgE,2) + pow(v_VVE[1],2));
    n2MCbkg = n2MCbkg + v_VV[2]; 
    n2MCbkgE = sqrt(pow(n2MCbkgE,2) + pow(v_VVE[2],2));

    n0MCbkg  = n0MCbkg + v_DY[0];
    n0MCbkgE = sqrt(pow(n0MCbkgE,2) + pow(v_DYE[0],2));
    n1MCbkg  = n1MCbkg + v_DY[1];
    n1MCbkgE = sqrt(pow(n1MCbkgE,2) + pow(v_DYE[1],2));
    n2MCbkg = n2MCbkg + v_DY[2]; 
    n2MCbkgE = sqrt(pow(n2MCbkgE,2) + pow(v_DYE[2],2));


    //now do the MC samples that are not data and not signal
    for(int iH=0; iH< nHists; ++iH){
      
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TString process(h1F->GetName());
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
      bool isSig = std::string(sampleNs->At(0)->GetName()) == std::string(signalS);
      
      if(isSig) //skip signal
	continue;
      if(TString(sampleNs->At(0)->GetName()).Contains("data")) //skip data
	continue;

      if(combineDYsamples && process.Contains("DY"))
	continue;
      if(combineVVsamples && process.Contains("ww"))
	continue;
      if(combineVVsamples && process.Contains("wz"))
	continue;
      if(combineVVsamples && process.Contains("zz"))
	continue;

      double n0 = h1F->GetBinContent(1);
      double n0E = h1F->GetBinError(1);
      n0MCbkg = n0MCbkg + n0;
      n0MCbkgE = sqrt(n0MCbkgE*n0MCbkgE + n0E*n0E);
      
      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);
      n1MCbkg = n1MCbkg + n1;
      n1MCbkgE = sqrt(n1MCbkgE*n1MCbkgE + n1E*n1E);


      double n2 = 0;
      double n2E = 0;      
      for(unsigned int i = 3; i <= h1F->GetNbinsX() + 1; i++) {
	n2 = n2 + h1F->GetBinContent(i);
	double temp = h1F->GetBinError(i);
	n2E = sqrt(n2E*n2E + temp*temp);
      }
      n2MCbkg = n2MCbkg + n2;
      n2MCbkgE = sqrt(n2MCbkgE*n2MCbkgE + n2E*n2E);

      
      std::string firstCol= Form("%9s ",sampleNs->At(0)->GetName());
      std::cout << beginL;
      if (showFirstCol) std::cout<< firstCol << colSep;
      if(!combineJetBins) {
	std::cout << mathSep << formatFloat(n0,formatS) <<pmSign<< formatFloat(n0E,formatS)<<mathSep<<colSep
		  << mathSep << formatFloat(n1,formatS) <<pmSign<< formatFloat(n1E,formatS)<<mathSep<<colSep
		  << mathSep << formatFloat(n2,formatS) <<pmSign<< formatFloat(n2E,formatS)<<mathSep
		  << endL
		  << std::endl;
      } else {
	double n = n0 + n1 + n2;
	double nE = sqrt(n0E*n0E + n1E*n1E + n2E*n2E);
	std::cout << mathSep << formatFloat(n,formatS) << pmSign << formatFloat(nE,formatS) << mathSep 
		  << endL << std::endl;
      }
    }//MC non signal only

    //now print the total MC
    std::string firstCol= Form("%9s ","Total MC");
    std::cout << beginL;
    if (showFirstCol) std::cout<< firstCol << colSep;
    if(!combineJetBins) {
      std::cout << mathSep << formatFloat(n0MCbkg,formatS) <<pmSign<< formatFloat(n0MCbkgE,formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(n1MCbkg,formatS) <<pmSign<< formatFloat(n1MCbkgE,formatS)<<mathSep<<colSep
		<< mathSep << formatFloat(n2MCbkg,formatS) <<pmSign<< formatFloat(n2MCbkgE,formatS)<<mathSep
		  << endL
		<< std::endl;
    } else {
      double n = n0MCbkg + n1MCbkg + n2MCbkg;
      double nE = sqrt(n0MCbkgE*n0MCbkgE + n1MCbkgE*n1MCbkgE + n2MCbkgE*n2MCbkgE);
      std::cout << mathSep << formatFloat(n,formatS) << pmSign << formatFloat(nE,formatS) << mathSep 
		<< endL << std::endl;
      }



    //now do the data, if we find it
    for(int iH=0; iH< nHists; ++iH){
      
      TH1F* h1F = (TH1F*)(thisStack->GetHists()->At(iH));
      TString process(h1F->GetName());
      TObjArray* sampleNs = TString(h1F->GetName()).Tokenize("_");
            
      if(!TString(sampleNs->At(0)->GetName()).Contains("data")) //do data only
	continue;
      
      double n0 = h1F->GetBinContent(1);
      double n0E = h1F->GetBinError(1);
      
      double n1 = h1F->GetBinContent(2);
      double n1E = h1F->GetBinError(2);

      double n2 = 0;
      double n2E = 0;      
      for(unsigned int i = 3; i <= h1F->GetNbinsX() + 1; i++) {
	n2 = n2 + h1F->GetBinContent(i);
	double temp = h1F->GetBinError(i);
	n2E = sqrt(n2E*n2E + temp*temp);
      }
      
      std::string firstCol= Form("%9s ",sampleNs->At(0)->GetName());
      std::cout << beginL;
      if (showFirstCol) std::cout<< firstCol << colSep;
      if(!combineJetBins) {
	//0 jet bin
	std::cout << mathSep << formatFloat(n0,formatS);
	if(printErrorsForData)
	  cout <<pmSign<< formatFloat(n0E, formatS);
	cout << mathSep<<colSep;

	//1jet bin
	cout << mathSep << formatFloat(n1,formatS);
	if(printErrorsForData)
	  cout <<pmSign<< formatFloat(n1E,formatS);
	cout <<mathSep<<colSep;

	//2nd jet bin
	cout << mathSep << formatFloat(n2,formatS);
	if(printErrorsForData)
	  cout <<pmSign<< formatFloat(n2E,formatS);
	cout << mathSep << colSep;
	cout << endL<< std::endl;
      } else {
	double n = n0 + n1 + n2;
	double nE = sqrt(n0E*n0E + n1E*n1E + n2E*n2E);
	std::cout << mathSep << formatFloat(n,formatS);
	if(printErrorsForData)
	  cout <<pmSign<< formatFloat(nE,formatS);
	cout << mathSep;
	cout << endL << std::endl;
      }
    }//data only

   

  }//over samples

}

/*
void printN2JetsColumns( const std::vector<std::string>& pfxs, int nPfx, bool latex = false, bool noErrs = false){
char* suffix[4];
suffix[0] = "ee";
suffix[1] = "mm";
suffix[2] = "em";
suffix[3] = "all";


if (latex) {
std::cout << "\\begin{table}" << std::endl;
std::cout << "\\begin{center}" << std::endl;
std::cout << "{\\small" << std::endl;
}

     
for (int sample=0; sample<4; sample++) {
    
THStack* thisStack[10]; memset(thisStack, 0, sizeof(thisStack));
for (int iP = 0; iP< nPfx; ++iP){
hist::stack(Form("st_%s_hnJet_%s", pfxs[iP].c_str(), suffix[sample]), 
  Form("%s_[a-zA-Z0-9]*_hnJet_%s$", pfxs[iP].c_str(), suffix[sample]));
thisStack[iP] = (THStack*) 
  gROOT->FindObjectAny(Form("st_%s_hnJet_%s", pfxs[iP].c_str(), suffix[sample]));
}
//    std::cout<<"Found "<<thisStack->GetName()<<std::endl;

int nHists = thisStack[0]->GetHists()->GetSize();
if (latex) {
std::cout << "\\begin{minipage}{0.48\\textwidth}" << std::endl;
std::cout << "\\begin{tabular}{|l|";
for(int iP=0;iP<nPfx; ++iP) std::cout<<"c|";
std::cout<<"} \\hline" << std::endl;
std::cout << "\\hline \\hline" << std::endl;
std::cout << "\\multicolumn{"<< nPfx+1 <<"}{|l|}{" << suffix[sample] << " final state} \\"<<"\\  " << std::endl;
std::cout<<" & ";
for (int iP = 0; iP < nPfx; ++iP) {
std::cout<< pfxs[iP].c_str();
if (iP != nPfx - 1) std::cout<<" & ";
}
std::cout <<"\\\\ \\hline" <<std::endl;
} else {
std::cout<<"===================================================="<<std::endl;
std::cout<<suffix[sample]<<std::endl;
std::cout<<" |   sample  |";
for (int iP = 0; iP < nPfx; ++iP){
std::cout<< "        "<<pfxs[iP].c_str()<<"        |";
}
std::cout<<std::endl;
}
double n2all[10]; memset(n2all, 0, sizeof(n2all));
double n2allE[10]; memset(n2allE, 0, sizeof(n2allE));
for(int iH=0; iH< nHists; ++iH){
TH1F* h1F[10];
for (int iP = 0; iP < nPfx; ++iP){
h1F[iP] = (TH1F*)(thisStack[iP]->GetHists()->At(iH));
}
TObjArray* sampleNs = TString(h1F[0]->GetName()).Tokenize("_");

//fisrt column first
if (latex) {
std::cout<<Form("%9s",sampleNs->At(1)->GetName()) << " & ";
} else {
std::cout<<" | "<<Form("%9s",sampleNs->At(1)->GetName())<<" | ";
}

int nBins = h1F[0]->GetNbinsX();
for (int iP = 0; iP < nPfx; ++iP){
double n2 = 0; 
for (int i=3; i<= nBins+1; ++i) n2+= h1F[iP]->GetBinContent(i);
double n2E = 0; 
for (int i=3; i<= nBins+1; ++i) n2E += h1F[iP]->GetBinError(i)*h1F[iP]->GetBinError(i);
n2E = sqrt(n2E);
n2all[iP] += n2;
n2allE[iP] += n2E*n2E;
if (latex) {
std::cout <<" $"<<Form("%6.1f",n2);
if (! noErrs){
std::cout <<" \\pm "<< Form("%6.1f",n2E);
}
std::cout <<" $";
if (iP != nPfx - 1){
std::cout<< " & ";
} else {
std::cout <<" \\"<<"\\"<<std::endl;
}
} else {
std::cout <<Form("%6.1f",n2) <<" &plusmn; "<<Form("%6.1f",n2E) <<" | ";
if (iP == nPfx - 1) std::cout<<std::endl;
}
}
}
if (latex) {
std::cout<<"\\hline"<<std::endl;
std::cout<<Form("%9s","Total")<<" & ";
for (int iP = 0; iP < nPfx; ++iP){
std::cout<<" $"<<Form("%6.1f",n2all[iP]);
if (! noErrs) {
std::cout<<" \\pm "<< Form("%6.1f",n2allE[iP]);
}
std::cout<<" $";
if (iP != nPfx - 1){
std::cout<< " & ";
} else {
std::cout <<" \\"<<"\\"<<std::endl;
}
}
std::cout << "\\hline " << std::endl;
std::cout << "\\end{tabular}" << std::endl;
std::cout << "\\end{minipage}" << std::endl;
} else {
std::cout<<" | "<<Form("%9s","Total")<<" | ";
for (int iP = 0; iP < nPfx; ++iP){
std::cout<<Form("%6.1f",n2all[iP]) <<" &plusmn; "<<Form("%6.1f",n2allE[iP]) <<" | ";
if (iP == nPfx - 1) std::cout<<std::endl;
}
}
if(!latex)
  std::cout << " " << std::endl;
}
 
if (latex) {
std::cout << "}" <<std::endl; //close \small fontsize
std::cout << "\\caption{\\label{tab:dummyLabel} Put your caption here}" << std::endl;
std::cout << "\\end{center}" << std::endl;
std::cout << "\\end{table}" << std::endl;
}
}


void getJESSyst(const char* formatS = "% 4.f"){
//input the names of the files here
hist::loadHist("sk_tdilgp032309/myHist_1957888__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15.root", 
  "mid", "*_hnJet_*");
hist::loadHist("sk_tdilgp032309/myHist_10346496__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15_jets10Up.root", 
  "jup", "*_hnJet_*");
hist::loadHist("sk_tdilgp032309/myHist_18735104__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15_jets10Dn.root", 
  "jdn", "*_hnJet_*");

char* suffix[4];
suffix[0] = "ee";
suffix[1] = "mm";
suffix[2] = "em";
suffix[3] = "all";
char* suffixLatex[4] = { "$\\eepm$", "$\\mmpm$", "$\\empm$", "All"};

TH1F* ttdil_mid[4];
TH1F* ttdil_jup[4];
TH1F* ttdil_jdn[4];

TH1F* tw_mid[4];
TH1F* tw_jup[4];
TH1F* tw_jdn[4];

TH1F* DYtautau_mid[4];
TH1F* DYtautau_jup[4];
TH1F* DYtautau_jdn[4];

TH1F* ww_mid[4];
TH1F* ww_jup[4];
TH1F* ww_jdn[4];

TH1F* wz_mid[4];
TH1F* wz_jup[4];
TH1F* wz_jdn[4];

TH1F* zz_mid[4];
TH1F* zz_jup[4];
TH1F* zz_jdn[4];

for (int iSam = 0; iSam < 4; ++iSam){
ttdil_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_ttdil_hnJet_%s",suffix[iSam]));
ttdil_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_ttdil_hnJet_%s",suffix[iSam]));
ttdil_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_ttdil_hnJet_%s",suffix[iSam]));

tw_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_tW_hnJet_%s",suffix[iSam]));
tw_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_tW_hnJet_%s",suffix[iSam]));
tw_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_tW_hnJet_%s",suffix[iSam]));

DYtautau_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_DYtautau_hnJet_%s",suffix[iSam]));
DYtautau_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_DYtautau_hnJet_%s",suffix[iSam]));
DYtautau_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_DYtautau_hnJet_%s",suffix[iSam]));

ww_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_ww_hnJet_%s",suffix[iSam]));
ww_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_ww_hnJet_%s",suffix[iSam]));
ww_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_ww_hnJet_%s",suffix[iSam]));

wz_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_wz_hnJet_%s",suffix[iSam]));
wz_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_wz_hnJet_%s",suffix[iSam]));
wz_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_wz_hnJet_%s",suffix[iSam]));

zz_mid[iSam] = (TH1F*)gDirectory->Get(Form("mid_zz_hnJet_%s",suffix[iSam]));
zz_jup[iSam] = (TH1F*)gDirectory->Get(Form("jup_zz_hnJet_%s",suffix[iSam]));
zz_jdn[iSam] = (TH1F*)gDirectory->Get(Form("jdn_zz_hnJet_%s",suffix[iSam]));

}

//  if(latex)
// just latex for now
std::cout<<  "\\begin{table}"                                              <<std::endl;
std::cout<<  "\\begin{center}"                                             <<std::endl;
std::cout<<  "\\begin{tabular}{|l|cc|cc|cc|cc|cc|cc|}"                              <<std::endl;
std::cout<<  "\\hline "                                                    <<std::endl;
std::cout<<  " & \\multicolumn{4}{|c|}{\\ttbar} & ";
std::cout<<  "   \\multicolumn{4}{|c|}{single-top} & ";
std::cout<<  "   \\multicolumn{4}{|c|}{$VV+\\dytt$} \\\\"<<std::endl;

std::cout<<  " & \\multicolumn{2}{|c|}{$N_\\jets=0,1$} & ";
std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} & "<<std::endl;
std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets=0,1$} & ";
std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} & "<<std::endl;
std::cout<<  " \t  \\multicolumn{2}{|c|}{$N_\\jets=0,1$} &";
std::cout<<  "   \\multicolumn{2}{|c|}{$N_\\jets\\geq 2$} \\\\"<<std::endl;
std::cout<<  "Channel & $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ & ";
std::cout<<  "  $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ & " <<std::endl;
std::cout<<  "  $-10\\%$ & $+10\\%$ & $-10\\%$ & $+10\\%$ \\\\ " <<std::endl;
std::cout<<  "\\hline "                                                    <<std::endl;
for (int iSam=0; iSam< 4; ++iSam){
if (iSam == 3) std::cout<<"\\hline"<<std::endl;
std::cout<<suffixLatex[iSam]<<" & ";
std::cout<< Form("$% 4.2f$ & ", (ttdil_jdn[iSam]->Integral(0,2)/ttdil_mid[iSam]->Integral(0,2) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (ttdil_jup[iSam]->Integral(0,2)/ttdil_mid[iSam]->Integral(0,2) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (ttdil_jdn[iSam]->Integral(3,9)/ttdil_mid[iSam]->Integral(3,9) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (ttdil_jup[iSam]->Integral(3,9)/ttdil_mid[iSam]->Integral(3,9) - 1.)*100.);

std::cout<< Form("$% 4.2f$ & ", (tw_jdn[iSam]->Integral(0,2)/tw_mid[iSam]->Integral(0,2) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (tw_jup[iSam]->Integral(0,2)/tw_mid[iSam]->Integral(0,2) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (tw_jdn[iSam]->Integral(3,9)/tw_mid[iSam]->Integral(3,9) - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (tw_jup[iSam]->Integral(3,9)/tw_mid[iSam]->Integral(3,9) - 1.)*100.);

//0,1
double dyttVV_mid = DYtautau_mid[iSam]->Integral(0,2);
dyttVV_mid += ww_mid[iSam]->Integral(0,2);
dyttVV_mid += wz_mid[iSam]->Integral(0,2);
dyttVV_mid += zz_mid[iSam]->Integral(0,2);
double dyttVV_jup = DYtautau_jup[iSam]->Integral(0,2);
dyttVV_jup += ww_jup[iSam]->Integral(0,2);
dyttVV_jup += wz_jup[iSam]->Integral(0,2);
dyttVV_jup += zz_jup[iSam]->Integral(0,2);
double dyttVV_jdn = DYtautau_jdn[iSam]->Integral(0,2);
dyttVV_jdn += ww_jdn[iSam]->Integral(0,2);
dyttVV_jdn += wz_jdn[iSam]->Integral(0,2);
dyttVV_jdn += zz_jdn[iSam]->Integral(0,2);
std::cout<< Form("$% 4.2f$ & ", (dyttVV_jdn/dyttVV_mid - 1.)*100.);
std::cout<< Form("$% 4.2f$ & ", (dyttVV_jup/dyttVV_mid - 1.)*100.);
//gt2
dyttVV_mid = DYtautau_mid[iSam]->Integral(3,9);
dyttVV_mid += ww_mid[iSam]->Integral(3,9);
dyttVV_mid += wz_mid[iSam]->Integral(3,9);
dyttVV_mid += zz_mid[iSam]->Integral(3,9);
dyttVV_jup = DYtautau_jup[iSam]->Integral(3,9);
dyttVV_jup += ww_jup[iSam]->Integral(3,9);
dyttVV_jup += wz_jup[iSam]->Integral(3,9);
dyttVV_jup += zz_jup[iSam]->Integral(3,9);
dyttVV_jdn = DYtautau_jdn[iSam]->Integral(3,9);
dyttVV_jdn += ww_jdn[iSam]->Integral(3,9);
dyttVV_jdn += wz_jdn[iSam]->Integral(3,9);
dyttVV_jdn += zz_jdn[iSam]->Integral(3,9);
std::cout<< Form("$% 4.2f$ & ", (dyttVV_jdn/dyttVV_mid - 1.)*100.);
std::cout<< Form("$% 4.2f$ \\\\ ", (dyttVV_jup/dyttVV_mid - 1.)*100.);

std::cout<<std::endl;
}
std::cout<<  "\\hline "                                                    <<std::endl;
std::cout<<  "\\end{tabular}"                                              <<std::endl;
std::cout<<  "\\caption{\\label{tab:dummyLabel} Relative changes in the number of selected events"<<std::endl;
std::cout<<  "expressed in $\\%$ for changes in the jet energy scale down and up by $10\\%$.}"    <<std::endl;
std::cout<<  "\\end{center}"                                               <<std::endl;
std::cout<<  "\\end{table}"                                                <<std::endl;
}
*/
