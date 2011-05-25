// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "TDatabasePDG.h"
#include <vector>
#include <string>

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
  CUT_LT_PT,
  CUT_LL_PT,
  CUT_SAME_SIGN,
  CUT_OPP_SIGN,
  CUT_PASS2_MET,
  CUT_PASS4_MET,
  CUT_PASS2_TCMET,
  CUT_PASS4_TCMET,
  CUT_LT_GOOD,
  CUT_LL_GOOD,
  CUT_EL_GOOD,
  CUT_MU_GOOD,
  CUT_LT_ISO,
  CUT_LL_ISO,
  CUT_ONE_ISO,
  CUT_TWO_ISO,
  CUT_EL_ISO,
  CUT_MU_ISO,
  CUT_LT_CALOISO,
  CUT_LL_CALOISO,
  CUT_ONE_CALOISO,
  CUT_TWO_CALOISO,
  CUT_EL_CALOISO,
  CUT_PASS_ZVETO,
  CUT_IN_Z_WINDOW,
  CUT_MUON_TAGGED,
  CUT_ELFAKE_FAKEABLE_OBJECT,
  CUT_ELFAKE_NUMERATOR,
  CUT_ELFAKE_NOT_NUMERATOR,
  CUT_MUFAKE_FAKEABLE_OBJECT,
  CUT_MUFAKE_NUMERATOR,
  CUT_MUFAKE_NOT_NUMERATOR,
  CUT_MUFAKE_LT_FAKEABLE_OBJECT,
  CUT_MUFAKE_LT_NUMERATOR,
  CUT_MUFAKE_LT_NOT_NUMERATOR,
  CUT_MUFAKE_LL_FAKEABLE_OBJECT,
  CUT_MUFAKE_LL_NUMERATOR,
  CUT_MUFAKE_LL_NOT_NUMERATOR,
  CUT_ELFAKE_LT_FAKEABLE_OBJECT,
  CUT_ELFAKE_LT_NUMERATOR,
  CUT_ELFAKE_LT_NOT_NUMERATOR,
  CUT_ELFAKE_LL_FAKEABLE_OBJECT,
  CUT_ELFAKE_LL_NUMERATOR,
  CUT_ELFAKE_LL_NOT_NUMERATOR,
  CUT_MORE_THAN_TWO_TRACKS,
  CUT_TRUE_MU_FROM_W,
  CUT_TRUE_MU_FROM_W_WJETS,
  CUT_TRUE_LT_MU_FROM_W,
  CUT_TRUE_LT_MU_FROM_W_WJETS,
  CUT_TRUE_LL_MU_FROM_W,
  CUT_TRUE_LL_MU_FROM_W_WJETS,
  CUT_TRUE_EL_FROM_W,
  CUT_TRUE_EL_FROM_W_WJETS,
  CUT_TRUE_LT_EL_FROM_W,
  CUT_TRUE_LT_EL_FROM_W_WJETS,
  CUT_TRUE_LL_EL_FROM_W,
  CUT_TRUE_LL_EL_FROM_W_WJETS,
  CUT_NOT_TRUE_GAMMA_FROM_MUON,
  CUT_NOT_TRUE_LT_GAMMA_FROM_MUON,
  CUT_NOT_TRUE_LL_GAMMA_FROM_MUON,
};


//----------------------------------------------------------------------
// Cut combinations for selections.  These are examples that are used
// for various tables in the WW analysis.
// ----------------------------------------------------------------------

// Note: some quick reminders for bitwise operations.  
//
// - To turn the enum into a cut bit, use the CUT_BIT(x) function.
//   For example, to set the bit for CUT_LT_PT, use CUT_BIT(CUT_LT_PT).
//
// - Combine bits with the bitwise OR operator (|).  For example, to
//   turn on the bits for CUT_LT_PT and CUT_LL_PT, use
//   CUT_BIT(CUT_LT_PT) | CUT_BIT(CUT_LL_PT).  For another example,
//   making a new set of cuts that is the same as an existing set with
//   one cut added (the hypothetical CUT_EXTRA_CUT), use the
//   following: 
//   cuts_t baseline_cuts_with_extra_cut = baseline_cuts | CUT_BIT(CUT_EXTRA_CUT);
//
// - To turn off a bit (useful when defining a set of cuts that is the
//   same as another set of cuts with one cut removed), use the
//   following: 
//   cuts_t baseline_cuts_without_lt_pt = baseline_cuts & ~CUT_BIT(CUT_LT_PT);
 
// define useful cut combinations here

// this is the current baseline set of cuts
const static cuts_t baseline_cuts = 
  (CUT_BIT(CUT_TRUE_MU_FROM_W_WJETS) ) | // back in 090311, IBL, why was this out?? Also added 20 GeV cut to muon leg.
  (CUT_BIT(CUT_LT_PT)		) | 
  (CUT_BIT(CUT_LL_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)		) | 
  //  (CUT_BIT(CUT_PASS4_MET)		) |  
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_CALOISO)	) |  
  (CUT_BIT(CUT_LL_CALOISO)	)
   |
    (CUT_BIT(CUT_NOT_TRUE_GAMMA_FROM_MUON))
;   

// denominator object cuts for the fake rate prediction 
const static cuts_t fakerate_denominator_cuts = (baseline_cuts & 
						 ~(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
						   CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) |
						   CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)
//                                                       |
//                                                       (CUT_BIT(CUT_NOT_TRUE_GAMMA_FROM_MUON))
                                                   ) ) |
  CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT) 
//  | CUT_BIT(CUT_EL_GOOD) |  CUT_BIT(CUT_EL_ISO) // removed 090226 IBL
;
						 
// numerator object cuts for the fake rate prediction 
const static cuts_t fakerate_numerator_cuts = 
  fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NUMERATOR);

// denominator and not numerator (this is the yield that should be
// multiplied by FR / (1 - FR))
const static cuts_t fakerate_denominator_not_numerator_cuts = 
  fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

// WW looper: make all the standard histograms and yields
class Looper : public LooperBase {

public:
  // constructor; tell the looper what sample to loop on (see
  // Tools/Sample.h), what cuts candidates need to pass, and a file
  // name for dumping log messages
  Looper (Sample, cuts_t cuts, const char *logfilename = 0);
  virtual ~Looper () { }

protected:
  // this is where we book our histograms
  virtual void	BookHistos ();

  // filter out this event.  If FilterEvent returns true, no
  // further processing is done on this event
  virtual bool	FilterEvent();

  // we define an analysis-specific EventSelect(), DilepSelect(),
  // TrilepSelect() and QuadlepSelect() that check which cuts the
  // event, dilepton/trilepton/quadlepton candidate passes
  virtual cuts_t	EventSelect	();
  virtual cuts_t	DilepSelect 	(int idx);
  virtual cuts_t	TrilepSelect 	(int idx);
  virtual cuts_t	QuadlepSelect 	(int idx);
  // we define an analysis-specific set of FillEventHistos(),
  // FillDilepHistos(), FillTrilepHistos() and FillQuadlepHistos()
  // that fill our histograms.  
  // 
  // the framework calls our FillEventHistos() function for every event
  virtual void	FillEventHistos ();
  // the framework calls our FillDilepHistos() function for every
  // dilepton candidate; the argument is the index of the candidate
  // in the dilepton block
  virtual void	FillDilepHistos (int idx);
  // the framework calls our FillTrilepHistos() function for every
  // trilepton candidate; the argument is the index of the candidate
  // in the trilepton block
  virtual void	FillTrilepHistos (int idx);
  // the framework calls our FillQuadlepHistos() function for every
  // quadlepton candidate; the argument is the index of the candidate
  // in the quadlepton block
  virtual void	FillQuadlepHistos (int idx);
  // at the end of the loop, we get a callback to do things like
  // printing a status message
  virtual void	End		();

  // weight for dilepton candidate
  virtual double 	Weight (int);

public:
  // these functions are called by the table-printing code
  virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
  virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
  virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------


  TH1F	*hnJet_[4];
  TH1F	*helPt_[4];
  TH1F	*helEta_[4];
  TH1F	*hmet_[4];
  TH1F	*helEtaPtBelow20_[4];
  TH1F	*helEtaPtAbove20_[4];
  

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double	cands_passing_[4];
  double	cands_passing_w2_[4];
  unsigned int	cands_count_[4];
//   std::vector<std::string> mybits;
  

};

#endif
