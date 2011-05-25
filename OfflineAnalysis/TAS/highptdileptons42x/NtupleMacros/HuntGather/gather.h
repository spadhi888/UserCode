#ifndef GATHER_H
#define GATHER_H

class BabySample;
class TCanvas;
class TChain;
class TCut;
class TH1F;

float GetIntLumi(float lumi, int brun, int bls, int erun, int els);
float GetIntLumi(float lumi);
TH1F* Plot(const char *pfx, const char *pfx2, TChain *chain, TCut field, TCut sel, TCut presel, float intlumifb, float kfactor,
           unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx = 0);
TH1F* Plot(TCut field, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
           BabySample *bs, unsigned int isfx = 0);
// These should only be used with data, where intlumifb need not be specified and 
// most likely you aren't scaling
TH1F* Plot(const char *pfx, const char *pfx2,TChain *chain, TCut field, TCut sel, TCut presel, float kfactor,
           unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx = 0);
TH1F* Plot(const char *pfx, const char *pfx2,TChain *chain, TCut field, TCut sel, TCut presel,
           unsigned int nbins, float xlo, float xhi, bool integrated, bool isdata, unsigned int isfx = 0);
TH1F* Plot(TCut field, TCut sel, TCut presel, unsigned int nbins, float xlo, float xhi, bool integrated,
           BabySample *bs, unsigned int isfx = 0);

bool sortHistsByIntegral(TH1* h1, TH1* h2);
TH1F* slideIntegrated(TH1F*);
TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
           std::vector<BabySample*> bss);
TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated,
		 BabySample *bs1,    BabySample *bs2 =0, BabySample *bs3 =0, BabySample *bs4 =0, BabySample *bs5 =0,
		 BabySample *bs6 =0, BabySample *bs7 =0, BabySample *bs8 =0, BabySample *bs9 =0, BabySample *bs10=0,
		 BabySample *bs11=0, BabySample *bs12=0, BabySample *bs13=0, BabySample *bs14=0, BabySample *bs15=0,
		 BabySample *bs16=0, BabySample *bs17=0, BabySample *bs18=0, BabySample *bs19=0, BabySample *bs20=0);

// Predefines and uses what are most likley the only BabySamples
// one needs for gathering
TCanvas* DrawAll(TCut field, const char *savename, TCut sel, TCut presel, float intlumifb, unsigned int nbins, float xlo, float xhi, bool integrated);

unsigned int gDrawAllCount = 0;

#endif
