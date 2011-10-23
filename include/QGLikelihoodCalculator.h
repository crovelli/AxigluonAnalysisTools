// ------------------------------------------------------------
//  
//    QGLikelihoodCalculator - Class
//    for the computation of the QG likelihood.
//    Needs files provided by having run the
//    Ntp1Finalizer_QG on QCD samples.
//
// ------------------------------------------------------------

#include <string>

#include "TFile.h"
#include "TH1F.h"

class QGLikelihoodCalculator {

 public:

  QGLikelihoodCalculator( const std::string& fileName="./data/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root", unsigned nPtBins=20, unsigned int nRhoBins=17 );
  virtual ~QGLikelihoodCalculator() {};

  float computeQGLikelihood( float pt, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );
  float computeQGLikelihoodPU( float pt, float rhoPF, int nCharged, int nNeutral, float ptD, float rmsCand=-1. );

  float likelihoodProduct( float nCharged, float nNeutral, float ptD, float rmsCand, TH1F* h1_nCharged, TH1F* h1_nNeutral, TH1F* h1_ptD, TH1F* h1_rmsCand);

 private:

  TFile* histoFile_;

  unsigned int nPtBins_;
  unsigned int nRhoBins_;

};

