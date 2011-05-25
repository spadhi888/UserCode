{
//---------------------------
// Load some useful tools
//----------------------------
gROOT->LoadMacro("eff2.C");
gSystem->CompileMacro("histtools.C", "++k", "libhisttools");
gStyle->SetOptStat(0);

//-------------------------------------------
// Here you load the lepton data that
// you want to use to make your fake rate.
// It should be a baby ntuple.  
// Make sure you have one available
//-------------------------------------------
TChain *ch2 = new TChain("tree");
//ch2->Add("QCDpt30_pfjet.root");
ch2->Add("QCDpt30.root");

//----------------------------------------------------
// Here we define cuts for numerator and denominator
//-----------------------------------------------------
// A cut against Ws
TCut notWCut = "tcmet<20 && mt<25"; 

// A pt cut...
// Remember, we use 35 (at least for muons) to
// minimize the impact of Ws
TCut ptCut   = "pt>10 && pt<35";

// Only consider events with >=1 jet of uncorrected
// calojet pt above a threshold separated by at least
// dR...here the threshold is set to 15 GeV
TCut jetCut20  = "ptpfj1cor>20";
TCut jetCut40  = "ptpfj1cor>40";
TCut jetCut60  = "ptpfj1cor>60";

// The trigger selection
TCut trgCutEle = "(el10_lw>1||el10_sw>1)";
TCut trgCutMu = "mu9>1";

// The 3-charge consistency 
TCut q3Cut = "q3";



// ------------------------------------------------
// electron
// ------------------------------------------------
// The numerator selection
TCut isNum20 = "numSSAug9&&abs(id)==11"+trgCutEle+jetCut20+notWCut+q3Cut;
TCut isNum40 = "numSSAug9&&abs(id)==11"+trgCutEle+jetCut40+notWCut+q3Cut;
TCut isNum60 = "numSSAug9&&abs(id)==11"+trgCutEle+jetCut60+notWCut+q3Cut;

//The denominator selection
TCut isDenomv120 = "v1SSAug9&&abs(id)==11"+trgCutEle+jetCut20+notWCut;
TCut isDenomv220 = "v2SSAug9&&abs(id)==11"+trgCutEle+jetCut20+notWCut;
TCut isDenomv320 = "v3SSAug9&&abs(id)==11"+trgCutEle+jetCut20+notWCut;
TCut isDenomv140 = "v1SSAug9&&abs(id)==11"+trgCutEle+jetCut40+notWCut;
TCut isDenomv240 = "v2SSAug9&&abs(id)==11"+trgCutEle+jetCut40+notWCut;
TCut isDenomv340 = "v3SSAug9&&abs(id)==11"+trgCutEle+jetCut40+notWCut;
TCut isDenomv160 = "v1SSAug9&&abs(id)==11"+trgCutEle+jetCut60+notWCut;
TCut isDenomv260 = "v2SSAug9&&abs(id)==11"+trgCutEle+jetCut60+notWCut;
TCut isDenomv360 = "v3SSAug9&&abs(id)==11"+trgCutEle+jetCut60+notWCut;

//-------------------------------------------------------
// Now you define the pt and eta bins for your fake rate
//-------------------------------------------------------
double ybinel[5]={10.,15.,20.,25.,35.};
int nbinsyel = 4;
double xbinel[5]={0.0, 1.0, 1.479, 2.0, 2.5};
int nbinsxel = 4;

// create an array of Histograms 
TObjArray Hlist(0);      

//--------------------------------------------------------
// Book your numerator and denominator histograms
//--------------------------------------------------------
TH2F* num20 = new TH2F("num20","num20",  nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov120  = new TH2F("fov120", "fov120",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov220  = new TH2F("fov220", "fov220",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov320  = new TH2F("fov320", "fov320",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* num40 = new TH2F("num40","num40",  nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov140  = new TH2F("fov140", "fov140",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov240  = new TH2F("fov240", "fov240",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov340  = new TH2F("fov340", "fov340",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* num60 = new TH2F("num60","num60",  nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov160  = new TH2F("fov160", "fov160",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov260  = new TH2F("fov260", "fov260",   nbinsxel, xbinel, nbinsyel, ybinel);
TH2F* fov360  = new TH2F("fov360", "fov360",   nbinsxel, xbinel, nbinsyel, ybinel);

//------------------------------------------
// Fill the Histograms
//-------------------------------------------
ch2->Draw("pt:abs(eta)>>num20",isNum20);
ch2->Draw("pt:abs(eta)>>fov120", isDenomv120);
ch2->Draw("pt:abs(eta)>>fov220", isDenomv220);
ch2->Draw("pt:abs(eta)>>fov320", isDenomv320);
ch2->Draw("pt:abs(eta)>>num40",isNum40);
ch2->Draw("pt:abs(eta)>>fov140", isDenomv140);
ch2->Draw("pt:abs(eta)>>fov240", isDenomv240);
ch2->Draw("pt:abs(eta)>>fov340", isDenomv340);
ch2->Draw("pt:abs(eta)>>num60",isNum60);
ch2->Draw("pt:abs(eta)>>fov160", isDenomv160);
ch2->Draw("pt:abs(eta)>>fov260", isDenomv260);
ch2->Draw("pt:abs(eta)>>fov360", isDenomv360);

//------------------------------------------
// Get the fake rate
// The output histogram name is "fr"
//------------------------------------------
TH2F* frv120 = eff2(fov120,num20,"frv120");
TH2F* frv140 = eff2(fov140,num40,"frv140");
TH2F* frv160 = eff2(fov160,num60,"frv160");
TH2F* frv220 = eff2(fov220,num20,"frv220");
TH2F* frv240 = eff2(fov240,num40,"frv240");
TH2F* frv260 = eff2(fov260,num60,"frv260");
TH2F* frv320 = eff2(fov320,num20,"frv320");
TH2F* frv340 = eff2(fov340,num40,"frv340");
TH2F* frv360 = eff2(fov360,num60,"frv360");



//------------------------------------------
// Muons
//------------------------------------------
// The numerator selection
TCut isNum20Mu = "numSS&&abs(id)==13"+trgCutMu+jetCut20+notWCut+ptCut;
TCut isNum40Mu = "numSS&&abs(id)==13"+trgCutMu+jetCut40+notWCut+ptCut;
TCut isNum60Mu = "numSS&&abs(id)==13"+trgCutMu+jetCut60+notWCut+ptCut;

//The denominator selection
TCut isDenom20Mu04 = "fo_04&&abs(id)==13"+trgCutMu+jetCut20+notWCut+ptCut;
TCut isDenom20Mu10 = "fo_10&&abs(id)==13"+trgCutMu+jetCut20+notWCut+ptCut;
TCut isDenom40Mu04 = "fo_04&&abs(id)==13"+trgCutMu+jetCut40+notWCut+ptCut;
TCut isDenom40Mu10 = "fo_10&&abs(id)==13"+trgCutMu+jetCut40+notWCut+ptCut;
TCut isDenom60Mu04 = "fo_04&&abs(id)==13"+trgCutMu+jetCut60+notWCut+ptCut;
TCut isDenom60Mu10 = "fo_10&&abs(id)==13"+trgCutMu+jetCut60+notWCut+ptCut;


//-------------------------------------------------------
// Now you define the pt and eta bins for your fake rate
//-------------------------------------------------------
double ybinmu[5]={10.,15.,20.,25.,35.};
int nbinsymu = 4;
//double xbinmu[4]={0.0, 1.0, 1.5, 2.5};
double xbinmu[4]={0.0, 1.0, 1.479, 2.5};
int nbinsxmu = 3;

// create an array of Histograms 
//TObjArray Hlist(0);      

//--------------------------------------------------------
// Book your numerator and denominator histograms
//--------------------------------------------------------
TH2F* num20Mu = new TH2F("num20Mu","num20Mu",  nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo20Mu04  = new TH2F("fo20Mu04", "fo20Mu04",   nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo20Mu10  = new TH2F("fo20Mu10", "fo20Mu10",   nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* num40Mu = new TH2F("num40Mu","num40Mu",  nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo40Mu04  = new TH2F("fo40Mu04", "fo40Mu04",   nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo40Mu10  = new TH2F("fo40Mu10", "fo40Mu10",   nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* num60Mu = new TH2F("num60Mu","num60Mu",  nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo60Mu04  = new TH2F("fo60Mu04", "fo60Mu04",   nbinsxmu, xbinmu, nbinsymu, ybinmu);
TH2F* fo60Mu10  = new TH2F("fo60Mu10", "fo60Mu10",   nbinsxmu, xbinmu, nbinsymu, ybinmu);

//------------------------------------------
// Fill the Histograms
//-------------------------------------------
ch2->Draw("pt:abs(eta)>>num20Mu",isNum20Mu);
ch2->Draw("pt:abs(eta)>>fo20Mu04", isDenom20Mu04);
ch2->Draw("pt:abs(eta)>>fo20Mu10", isDenom20Mu10);
ch2->Draw("pt:abs(eta)>>num40Mu",isNum40Mu);
ch2->Draw("pt:abs(eta)>>fo40Mu04", isDenom40Mu04);
ch2->Draw("pt:abs(eta)>>fo40Mu10", isDenom40Mu10);
ch2->Draw("pt:abs(eta)>>num60Mu",isNum60Mu);
ch2->Draw("pt:abs(eta)>>fo60Mu04", isDenom60Mu04);
ch2->Draw("pt:abs(eta)>>fo60Mu10", isDenom60Mu10);

//------------------------------------------
// Get the fake rate
// The output histogram name is "fr"
//------------------------------------------
TH2F* fr20Mu04 = eff2(fo20Mu04,num20Mu,"fr20Mu04");
TH2F* fr20Mu10 = eff2(fo20Mu10,num20Mu,"fr20Mu10");
TH2F* fr40Mu04 = eff2(fo40Mu04,num40Mu,"fr40Mu04");
TH2F* fr40Mu10 = eff2(fo40Mu10,num40Mu,"fr40Mu10");
TH2F* fr60Mu04 = eff2(fo60Mu04,num60Mu,"fr60Mu04");
TH2F* fr60Mu10 = eff2(fo60Mu10,num60Mu,"fr60Mu10");


// open a file and save histograms using histtools.C
const char* outFile = "SSFakeRate.root";
hist::saveHist(outFile);
hist::deleteHistos();

}
