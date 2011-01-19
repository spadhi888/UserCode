//-----------------------------------------------------------------------------
// File:        pmssmanalyzer.cc
// Description: Analyzer for ntuples created by Mkntuple
// Created:     Fri Dec  3 11:22:14 2010 by mkntanalyzer.py
// Author:      Sezen Sekmen
// $Revision: 1.19 $
//-----------------------------------------------------------------------------
#include "pmssmanalyzer.h"

#ifdef PROJECT_NAME
#include "PhysicsTools/Mkntuple/interface/pdg.h"
#else
#include "pdg.h"
#endif

using namespace std;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line

  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples

  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Analysis GEN");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

  // Select variables to be read

  selectVariables(stream);

  //---------------------------------------------------------------------------
  // Book histograms etc.
  //---------------------------------------------------------------------------
  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the following
  // line

  TApplication app("analyzer", &argc, argv);

  histogramFile hfile(cmdline.histfilename);

  // Histograms


  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry < nevents; ++entry)
	{
	  // Read event into memory
	  stream.read(entry);
	 
	  // ---------------------
	  // -- Event Selection --
	  // ---------------------

	  // if ( !SUSY ) continue;
	}

  stream.close();
  hfile.close();
  return 0;
}
