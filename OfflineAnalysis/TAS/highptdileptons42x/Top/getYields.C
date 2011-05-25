#include <iostream>
#include "histtools.C"
#include "printHists.C"


void getYields(TString dataFName, TString mcFName, float scaleMC, bool combineJetBins = false, const char* formatS = "%6.3f", bool latex = true) {


  /*
    gROOT->ProcessLine(".L getMyHistosNames.C+");
    gROOT->ProcessLine(".L histtools.C++");
    gROOT->ProcessLine(".L browseStacks.C++");
    gROOT->ProcessLine(".L printHists.C++");
  */

    
    
    hist::loadHist(mcFName.Data(),0,"*_hnJet_*");
    hist::loadHist(mcFName.Data(),0,"*_htcmet_allj_*");
    hist::loadHist(mcFName.Data(),0,"*hdilMass_allj_*");
    hist::loadHist(mcFName.Data(),0,"*bTagtk*_allj_*");    
    hist::scale("*_*", scaleMC);
  
    hist::loadHist(dataFName.Data(),0,"data_hnJet_*");
    hist::loadHist(dataFName.Data(),0,"data_htcmet_allj_*");
    hist::loadHist(dataFName.Data(),0,"data_*hdilMass_allj_*");
    hist::loadHist(dataFName.Data(),0,"data_*bTagtk*_allj_*");
    
    
    std::vector<TString> v_prfxsToCombine;
    v_prfxsToCombine.push_back("qcd15");
    v_prfxsToCombine.push_back("qcd30");
    hist::combineHists(v_prfxsToCombine,"QCD");
    
    
    v_prfxsToCombine.clear();
    v_prfxsToCombine.push_back("DYee");
    v_prfxsToCombine.push_back("DYmm");
    hist::combineHists(v_prfxsToCombine, "DYeemm");
    
    
    
    /*
      void browseStacks( bool makePictures=false, bool wait=true , string referencePrefix="ttdil", 
		   TString outfilename = "stacks.ps",int addHistName = 1, 
		   Double_t maxYScaleF = 1.,  bool logScale = false, bool setMinZero = true, 
		   int colorScheme = 2, bool noLegend = false, int orderScheme = 0,
		   bool saveAsPNGs = false, char* bsmName = "") {
		   
		   
		   void printNJets( bool latex=false, const char* formatS = "%6.1f", const char* signalS= "ttdil", 
		   bool combineVVsamples = true, bool combineDYsamples = true, bool combineJetBins = false, bool printProbs=false, bool printErrorsForData = false){
    */
    
   
    
    printNJets(latex, formatS,"ttdil", true,false,combineJetBins, false); 
    //browseStacks( true, false , "ttdil", dataFName, 4, 7, true, false, 3, false, 0, true);
    //browseStacks( true, false , "DYeemm", dataFName, 4, 27, true, false, 3, false, 0, true);
    hist::deleteHistos();
    gDirectory->Clear();
    

    
  
}
