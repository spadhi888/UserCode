
#ifndef PLOTRESULTS_H
#define PLOTRESULTS_H

#include "TString.h"

class TCanvas;
class HistogramUtilities;
class THStack;
class TArrow;

void plotEff(HistogramUtilities &h1, TString name, TString saveName, TString det, bool ascending, int rebin = 1, bool legendOnRight = true, float cutValEB = -1.0, float cutValEE = -1.0);
void plotStack(HistogramUtilities &h1, TString name, TString titleX, TString saveName, TString det, int rebin = 1, float cutValEB = -1.0, float cutValEE = -1.0);


TArrow *getArrow(TString det, THStack *st,float cutValEB, float cutValEE);


void plot2DSB(HistogramUtilities &h1, TString name, TString xTitle, TString yTitle, TString saveName, TString det);

void plotAllResultsID();
void plotResultsID(TString det, TString fileStamp);

void plotAllResultsAN2009_098();
void plotResultsAN2009_098(TString det, TString fileStamp);

void plotResultsW(TString det, TString fileStamp);
void plotAllResultsW();


void plotAllResultsQCDVal();
void plotResultsQCDVal(TString det, TString fileStamp);


void test();

#endif

