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
  LOOP_LM1	,
  LOOP_LM2	,
  LOOP_LM3	,
  LOOP_LM4	,
  LOOP_LM5	,
  LOOP_LM6	,
  LOOP_LM7	,
  LOOP_LM8	,
  LOOP_LM9	,
  LOOP_LM10	,
  LOOP_LM11	,
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
  fprintf(f, "| %10s", "");
  for (int j = 0; j < n; ++j) {
    fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
  }
  fprintf(f, "|%30s  |\n", "SM total");
  for (int i = 0; i < 4; ++i) {
    fprintf(f, "|%10s  ", dilepton_hypo_names[i]);
    double cands = 0;
    double w2 = 0;
    for (int j = 0; j < n; ++j) {
      fprintf(f, "|  %10.1f &plusmn; %10.1f", 
	      hists[j]->CandsPassing(DileptonHypType(i)),
	      hists[j]->RMS(DileptonHypType(i)));
      // only sum up SM samples
      if ( hists[j]->SampleSM() ) {
	cands += hists[j]->CandsPassing(DileptonHypType(i));
	w2 += hists[j]->RMS(DileptonHypType(i)) * 
	  hists[j]->RMS(DileptonHypType(i));
      }
      if (not (which_ones & 1 << j))
	continue;
//       const DSGFakeRateLooper *looper = 
// 	dynamic_cast<const DSGFakeRateLooper *>(hists[j]);
//       if (looper != 0) {
// 	fprintf(f, " + %5.1f &minus; %5.1f", 
// 		looper->CandsPassingSystHi(DileptonHypType(i)) 
// 		- looper->CandsPassing(DileptonHypType(i)),
// 		looper->CandsPassing(DileptonHypType(i)) 
// 		- looper->CandsPassingSystLo(DileptonHypType(i)));
// 	fprintf(f, "(stat) &plusmn; %5.1f (fake)", 
// 		looper->FakeSyst(DileptonHypType(i)));
//       }
    }
    fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
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
template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
  const string hist = name + ".root";
  const string tbl = name + ".tbl";
  const string log = name + ".log";
  // by default, we run this list of samples; if we're told by the
  // which_ones bit field to skip a sample, we skip it
  Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
  Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
  Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
  Looper looper_wjets   	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
  Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
  Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
  Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
  Looper looper_ttbar	        (fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
  Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
  Looper looper_lm1		(fLM1()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM1   )) looper_lm1         .Loop();
  Looper looper_lm2		(fLM2()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM2   )) looper_lm2         .Loop();
  Looper looper_lm3		(fLM3()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM3   )) looper_lm3         .Loop();
  Looper looper_lm4		(fLM4()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM4   )) looper_lm4         .Loop();
  Looper looper_lm5		(fLM5()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM5   )) looper_lm5         .Loop();
  Looper looper_lm6		(fLM6()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM6   )) looper_lm6         .Loop();
  Looper looper_lm7		(fLM7()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM7   )) looper_lm7         .Loop();
  Looper looper_lm8		(fLM8()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM8   )) looper_lm8         .Loop();
  Looper looper_lm9		(fLM9()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM9   )) looper_lm9         .Loop();
  Looper looper_lm10		(fLM10()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM10  )) looper_lm10        .Loop();
  Looper looper_lm11		(fLM11()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM11  )) looper_lm11        .Loop();
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
    &looper_lm1         ,
    &looper_lm2         ,
    &looper_lm3         ,
    &looper_lm4         ,
    &looper_lm5         ,
    &looper_lm6         ,
    &looper_lm7         ,
    &looper_lm8         ,
    &looper_lm9         ,
    &looper_lm10        ,
    &looper_lm11        ,
  };
  printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
  return 0;
}

int BaseLine ()
{
  return run<Looper>(dsg_baseline_cuts, "BaseLine");
}

int MET_10 ()
{
  return run<Looper>(dsg_met_10_cuts, "MET_10");
}

int MET_1 ()
{
  return run<Looper>(dsg_met_1_cuts, "MET_1");
}

int SUMET_10 ()
{
  return run<Looper>(dsg_sumet_10_cuts, "SUMET_10");
}

int SUMET_1 ()
{
  return run<Looper>(dsg_sumet_1_cuts, "SUMET_1");
}
