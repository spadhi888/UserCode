#ifndef zmet_looper_h
#define zmet_looper_h

#include <vector>
#include <map>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
typedef map<unsigned int, unsigned int> m_uiui;

class  TChain;
class  TH1F;
class  TH2F;
class  TRandom3;
class  TTree;
struct metStruct;

class zmet_looper
{
    public: 
        zmet_looper();
        ~zmet_looper() {}

        enum JetTypeEnum { e_JPT = 0, e_calo };
        // e_JPT     :   jpt jets
        // e_calo    :   l1 and l2 corrected calo jets
        enum MetTypeEnum { e_tcmet = 0, e_muon, e_muonjes };
        // e_tcmet   :   track corrected met
        // e_muon    :   calo met with muon corrections
        // e_muonjes :   calo met with muon and jet energy scale corrections
        enum ZVetoEnum   { e_standard = 0, e_allzveto, e_nozveto, e_selectz };
        // e_standard:   apply Z-veto to same-flavor pairs
        // e_allzveto:   apply Z-veto regardless of lepton flavor
        // e_nozveto :   no Z-veto
        // e_selectz :   select Z by requiring SF OS pair in Z mass window

        int  ScanChain(TChain *chain, char *prefix = "", float kFactor = 1., int prescale = 1., 
                       JetTypeEnum jetType = e_JPT, 
                       MetTypeEnum metType = e_tcmet,
                       ZVetoEnum zveto = e_standard,
                       bool doFakeApp = false,
                       bool calculateTCMET = false);
        void BookHistos (char *prefix);
        bool passZSelection (int hypIdx);
        bool passTrigger (int dilType);
      
        // Set globals
        void set_susybaseline (bool b) { g_susybaseline = b; }
        void set_createTree   (bool b) { g_createTree   = b; }
        void set_useBitMask   (bool b) { g_useBitMask   = b; }

        // Baby ntuple methods
        void makeTree (char *prefix);
        void closeTree ();



    private:

        // Globals
        bool g_susybaseline;
        bool g_createTree;
        bool g_useBitMask;
        TRandom3 *random3_;

        // Baby ntuple variables
        TFile  *outFile;
        TTree  *outTree;
        Float_t weight_;
        Int_t   proc_;
        Int_t   leptype_;
        Int_t   njets_;
        Float_t dilmass_;
        Float_t tcmet_;
        Float_t tcsumet_;
        Float_t tcmetphi_;
        Float_t mt2_;
        Float_t mt2j_;
        Float_t sumjetpt_;
        Float_t dileta_;
        Float_t dilpt_;
        Float_t dildphi_;
        Float_t vecjetpt_;
        Int_t   pass_;
        Int_t   passz_;
        Float_t m0_;
        Float_t m12_;
        Float_t ptl1_;
        Float_t ptl2_;
        Float_t ptj1_;
        Float_t ptj2_;
        Float_t meff_;
        Float_t mt_;

	// the fake rate accessors

/* 	SimpleFakeRate fr_el */
/* 	SimpleFakeRate fr_mu; */

        // Lots and lots of histograms
        TH1F* hmt2j_signal[4][4];
        TH1F* hmt2j_control[4][4];
        TH1F* hmt2j_all[4][4];
        TH2F* hmet_dilpt_signal[4][4];
        TH2F* hmet_dilpt_control[4][4];
        TH2F* hmet_dilpt_all[4][4];
              
        TH1F* hmt[4][4];
        TH1F* hetaz[4][4];
        TProfile* htcsumet_tcmet_prof[4][4]; 
        TProfile* hsumJetPt_tcmetsqrtsumet_prof[4][4]; 
        TProfile* hgensumet_genmet_prof[4][4]; 
        //TProfile* hsumJetPt_tcmetpowtcsumet_prof[4][4][101]; 
        //TH2F* hsumJetPt_tcmetpowtcsumet_th2[4][4][101];
        //TH2F*     hsumJetPt_tcmetpowtcsumet_th2[101];
        //TProfile* hsumJetPt_tcmetpowtcsumet_prof[101];

        TH2F* hsumJetPt_tcmet[4][4];
        TH2F* hsumJetPt_tcmetsqrtsumet[4][4]; 
        TH2F* hsumJetPt_tcmetsumet[4][4]; 
        TH2F* hetaZ_tcmet[4][4]; 
        TH2F* hetaZ_tcmetsqrtsumet[4][4]; 
        TH2F* hetaZ_tcmetsumet[4][4]; 

        TH2F* hdilMass_tcmet[4][4];
        TH1F* hmt2core[4][4];               // MT2 from CORE
        TH1F* hmt2jcore[4][4];               // MT2J from CORE
        TH1F* hmt2j[4][4];                   // potentially custom MT2J
        TH1F* hsumJetPt[4][4];               // scalar sum jet Et
        TH2F* hDtcmetgenmetVsumJetPt[4][4];
        TH2F* hDmetmuonjesgenmetVsumJetPt[4][4];
        TH1F* hmeffJet[4][4];                // scalar sum jet pt + scalar sum dil pt + (tc)met

        TH1F* hsumJptPt[4][4];               // scalar sum JPT jet Et
        TH1F* hmeffJPT[4][4];                // scalar sum JPT jet pt + scalar sum dil pt + (tc)met

        TH1F* hsumHypPt[4][4];               // scalar sum JPT jet Et
        TH1F* hmeffHyp[4][4];                // scalar sum JPT jet pt + scalar sum dil pt + (tc)met

        TH1F* hnJet[4];                      // Njet distributions
        TH1F* hnJpt[4];                      // Njpt distributions
        TH1F* hnHypJet[4];                   // Hyp Njet distributions
        TH1F* helePt[4][4];                  // electron Pt
        TH1F* hmuPt[4][4];                   // muon Pt
        TH1F* hminLepPt[4][4];               // minimum lepton Pt
        TH1F* hmaxLepPt[4][4];               // maximum lepton Pt
        TH1F* helePhi[4][4];                 // electron phi
        TH1F* hmuPhi[4][4];                  // muon phi
        TH1F* hdphiLep[4][4];                // delta phi between leptons
        TH1F* heleEta[4][4];                 // electron eta
        TH1F* hmuEta[4][4];                  // muon eta
        TH1F* hdilMass[4][4];                // dilepton mass
        TH1F* hdilPt[4][4];                  // dilepton Pt
        TH1F* hdilPt_zveto[4][4];            // dilepton Pt with z-veto applied
        TH1F* hdilPtSmeared[4][4];           // dilepton Pt with Gaussian smearing

        TH1F* hgenmet[4][4];                 // MET corrected for muons and JES
        TH1F* hgenmetPhi[4][4];              // MET corrected for muons and JES phi
        TH1F* hmetmuon[4][4];                // MET corrected for muons and JES
        TH1F* hmetmuonPhi[4][4];             // MET corrected for muons and JES phi
        TH1F* hmetmuonjes[4][4];             // MET corrected for muons and JES
        TH1F* hmetmuonjesPhi[4][4];          // MET corrected for muons and JES phi
        TH1F* htcmet[4][4];                  // tc MET
        TH1F* htcmetPhi[4][4];               // tc MET phi
        TH1F* hpfmet[4][4];                  // tc MET
        TH1F* hpfmetPhi[4][4];               // tc MET phi

        TH2F* hmetmuonVsDilepPt[4][4];       // MET vs dilepton Pt
        TH2F* hmetmuonOverPtVsDphi[4][4];    // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
        TH2F* hmetmuonjesVsDilepPt[4][4];    // MET vs dilepton Pt
        TH2F* hmetmuonjesOverPtVsDphi[4][4]; // MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt
        TH2F* htcmetVsDilepPt[4][4];         // tc MET vs dilepton Pt
        TH2F* htcmetOverPtVsDphi[4][4];      // tc MET/Lepton Pt vs DeltaPhi between MET and Lepton Pt

        TH2F* hdphillvsmll[4][4];            // delta phi between leptons vs dilepton mass
        TH1F* hptJet1[4][4];                 // Pt of 1st jet
        TH1F* hptJet2[4][4];                 // Pt of 2nd jet
        TH1F* hptJet3[4][4];                 // Pt of 3rd jet
        TH1F* hptJet4[4][4];                 // Pt of 4th jet
        TH1F* hetaJet1[4][4];                // eta of 1st jet
        TH1F* hetaJet2[4][4];                // eta of 2nd jet
        TH1F* hetaJet3[4][4];                // eta of 3rd jet
        TH1F* hetaJet4[4][4];                // eta of 4th jet
        TH1F* hptJpt1[4][4];                 // Pt of 1st JPT jet
        TH1F* hptJpt2[4][4];                 // Pt of 2nd JPT jet
        TH1F* hptJpt3[4][4];                 // Pt of 3rd JPT jet
        TH1F* hptJpt4[4][4];                 // Pt of 4th JPT jet
        TH1F* hetaJpt1[4][4];                // eta of 1st JPT jet
        TH1F* hetaJpt2[4][4];                // eta of 2nd JPT jet
        TH1F* hetaJpt3[4][4];                // eta of 3rd JPT jet
        TH1F* hetaJpt4[4][4];                // eta of 4th JPT jet=
        TH1F* hptHypJet1[4][4];              // Pt of 1st JPT jet
        TH1F* hptHypJet2[4][4];              // Pt of 2nd JPT jet
        TH1F* hptHypJet3[4][4];              // Pt of 3rd JPT jet
        TH1F* hptHypJet4[4][4];              // Pt of 4th JPT jet
        TH1F* hetaHypJet1[4][4];             // eta of 1st JPT jet
        TH1F* hetaHypJet2[4][4];             // eta of 2nd JPT jet
        TH1F* hetaHypJet3[4][4];             // eta of 3rd JPT jet
        TH1F* hetaHypJet4[4][4];             // eta of 4th JPT jet
};

#endif
