/* Usage:
   root [0] .L ScanChain.C++
   root [1] TChain *chain = new TChain("Events")
   root [2] chain->Add("merged_ntuple.root")
   root [3] ScanChain(chain)
*/

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

// CMS2
#include "CMS2.cc"
using namespace tas;


int ScanChain( TChain* chain, int nEvents = -1, std::string skimFilePrefix="") {

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  TH1F *samplehisto = new TH1F("samplehisto", "Example histogram", 200,0,200);
  samplehisto->SetDirectory(rootdir);

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {
    // Get File Content
    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEvents = tree->GetEntries();
    for( unsigned int event = 0; event < nEvents; ++event) {
    
      // Get Event Content
      cms2.GetEntry(event);
      ++nEventsTotal;
    
      // Progress
      progress( nEventsTotal, nEventsChain );
    }
  
    delete tree;
    f.Close();
  }
  if ( nEventsChain != nEventsTotal ) {
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }

  
  // Example Histograms
  samplehisto->Draw();
  
  // return
  return 0;
}
