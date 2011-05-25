#include <unistd.h>
#include <string>
#include "Looper.h"
#include "Tools/Sample.h"
#include "Tools/tools.h"

using std::string;

// this enum says which samples should actually be used (to shorten
// looping time if you only care about the yields for one or two
// samples)
enum {
     LOOP_WW	,
     LOOP_WZ	,
     LOOP_ZZ	,
     LOOP_WJETS	,
     LOOP_DYEE	,
     LOOP_DYMM	,
     LOOP_DYTT	,
     LOOP_TTBAR	,
     LOOP_TW	,
};

// helper function used to print yield tables
void printTable (const Looper **hists, int n, const char *fname, 
		 uint32 which_ones)
{
     FILE *f = 0;
     if (fname == 0 || strlen(fname) == 0)
	  f = stdin;
     else f = fopen(fname, "w");
     if (f == 0) {
	  perror("printing table");
	  return;
     }
     fprintf(f, "|%4s", "");
     for (int j = 0; j < n; ++j) {
	  fprintf(f, "|*%13s*", hists[j]->SampleName().c_str());
     }
     //fprintf(f, "|%30s  |\n", "total");
	 fprintf(f, "|%13s  |\n", "total");
     for (int i = 0; i < 4; ++i) {
	   fprintf(f, "|%3s ", dilepton_hypo_names[i]);
	   double cands = 0;
	   double w2 = 0;
	   for (int j = 0; j < n; ++j) {
		 //fprintf(f, "|%10.1f &plusmn;%10.1f", 
		 fprintf(f, "| %5.1f +- %5.1f", 
				 hists[j]->CandsPassing(DileptonHypType(i)),
				 hists[j]->RMS(DileptonHypType(i)));
		 cands += hists[j]->CandsPassing(DileptonHypType(i));
		 w2 += hists[j]->RMS(DileptonHypType(i)) * 
		   hists[j]->RMS(DileptonHypType(i));
	  }
	  //fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
	  fprintf(f, "| %5.1f +- %5.1f|\n", cands, sqrt(w2));
     }
     if (f != stdin) 
	  fclose(f);
}

// run a looper on each sample and produce a yield table; arguments:
//
// class Looper: which type of looper to run (usually: Looper)
// cuts: cut definition from Looper.h (usually: baseline_cuts)
// name: name for the output files (usually: "Results", which produces Results.tbl, Results.root, Results.log)
// which_ones: which samples to run (usually: all); to run only WW and ttbar, use: (1 << LOOP_WW) | (1 << LOOP_TTBAR)
//
// examples:
// run<Looper>(baseline_cuts, "Results", 1 << LOOP_WW)				// produce table with default cuts, WW only
// run<Looper>(baseline_cuts, "Results", 1 << LOOP_WW | 1 << LOOP_WJETS)	// produce table with default cuts, WW and Wjets only
// run<Looper>(baseline_cuts, "Results")					// produce table with default cuts, all samples
//last argument added by me (Warren) for using narrow bins, axis range
template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff, bool usenar=false)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     Looper looper_ww	(fWW()		, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     Looper looper_wz	(fWZ()		, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     Looper looper_zz	(fZZ()		, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_wjets(fWjets()	, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     Looper looper_dyee	(fDYee()	, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm	(fDYmm()	, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dytt	(fDYtt()	, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_ttbar(fttbar()	, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw	(ftW()		, cuts, log.c_str(), usenar);	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     // then we collect them all and print a table
     const Looper *loopers[] = { 
	  &looper_ww          ,
	  &looper_wz          ,
	  &looper_zz          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results ()
{
  cuts_t baseline_cuts = 0;
  baseline_cuts = baseline_metres;
  //uint32 which_ones = 0;
  //which_ones |= (1 << LOOP_DYEE);
  //which_ones |= (1 << LOOP_DYMM);
  //which_ones |= (1 << LOOP_DYTT);
  //which_ones |= (1 << LOOP_TTBAR);
  //return run<Looper>(baseline_cuts, "Results", which_ones);
  return run<Looper>(baseline_cuts, "Results");
}

int Results_Nar()
{
  cuts_t baseline_cuts = 0;
  baseline_cuts = baseline_metres;
  return run<Looper>(baseline_cuts, "Results_Nar", 0xffffffff, true);
}
  
//don't use below two
int Susy_Results ()
{
  cuts_t susy_baseline_cuts = baseline_susy;
  uint32 which_ones = 0;
  which_ones |= (1 << LOOP_TTBAR);
  //return run<Looper>(susy_baseline_cuts, "Susy_Results", which_ones);
  return run<Looper>(susy_baseline_cuts, "Susy_Results");
}

int In_Zwindow ()
{
  cuts_t baseline_cuts_zwindow = 0;
     return run<Looper>(baseline_cuts_zwindow, "In_Zwindow");
}
