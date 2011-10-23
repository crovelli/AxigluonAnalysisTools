//-------------------------------------------------------

#ifndef AxigluonSelection_h
#define AxigluonSelection_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "AxigluonAnalysisTools/include/Axigluon.hh"
#include "AxigluonAnalysisTools/include/CutBasedAxiSelector.hh" 
#include "AxigluonAnalysisTools/include/QGLikelihoodCalculator.h"
#include "AxigluonAnalysisTools/include/RedAxiTree.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class AxigluonSelection : public Axigluon{
public:
  
  //! constructor
  AxigluonSelection(TTree *tree=0);
  //! destructor
  virtual ~AxigluonSelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies(std::string filename);
  //! set the required triggers masks (one per channel)
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  //! set the not-required triggers masks (one per channel)
  void setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  
private:

  //! get the hardest electron after the different steps
  int getBestElectronPair_acceptance();
  int getBestElectronPair_id( std::vector<int> acceptEle );
  int getBestElectronPair_isol( std::vector<int> idEle );
  int getBestElectronPair_conv( std::vector<int> isolEle );
  int getBestElectronPair_ip( std::vector<int> convEle );

  //! get the hardest muon after the different steps
  int getBestMuonPair_acceptance();
  int getBestMuonPair_id( std::vector<int> acceptMu ); 
  int getBestMuonPair_isol( std::vector<int> idMu );
  int getBestMuonPair_ip( std::vector<int> isoMu ); 

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsEle(int myEle);
  void setKinematicsMu(int myMuon);
  //! reset the kinematic quantities at the beginning of event and after the selection if needed
  void resetKinematicsStart();
  void resetKinematics();

  //! count jet multiplicity
  int numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) ;
  //! count the soft muons
  int numSoftMuons(std::vector<int> muonToRemove);
  //! count the extra leptons (id, iso, d0,acceptance etc) 
  int numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove );
  //! di-jet system proprierties
  float getDiJetInvMass(int theLJ, int theSJ);
  float getDiJetPt(int theLJ, int theSJ);
  float getDiJetDeta(int theLJ, int theSJ);
  //! set the electron ID variables to dump
  void setEleIdVariables(int ele);
  //! reload the trigger mask_s_ (one per channel)
  bool reloadTriggerMask(int runN);
  //! get the trigger answer depending on the channel
  bool hasPassedHLT(int channel);
  //! muonID
  bool isAxigluonMuonID(int muonIndex);
  //! to correct charged met
  TVector3 oneLepPFChargedMet(TVector3 lept);

  //! to evaluate eleID
  CutBasedEleIDSelector EgammaCutBasedID;

  //! to evaluate full selection efficiency
  Selection *_selectionEle, *_selectionMu;
  CutBasedAxiSelector CutBasedAxiSelectionEle;  
  CutBasedAxiSelector CutBasedAxiSelectionMu;   

  //! be verbose during runtime
  bool _verbose;

  //! mass hypotesis
  int _massVal;

  //! an integer defining the sub-channel
  enum { ee=0, mm=1 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[2];
  
  //! trigger masks
  std::vector<int>  m_requiredTriggersEle, m_requiredTriggersMu;
  std::vector<int>  m_notRequiredTriggersEle, m_notRequiredTriggersMu;
  std::vector<std::string> requiredTriggersEle, requiredTriggersMu;
  std::vector<std::string> notRequiredTriggersEle, notRequiredTriggersMu;

  //! leptons for the analysis
  int theElectron, theMuon;
  int thePreElectron, thePreMuon;

  //! kinematics of the event: W
  std::vector<int> eleCands[2], muCands[2];
  TLorentzVector *m_p4Lepton[2];
  TVector3 *m_p3Lepton[2];
  float leptonPt[2];
  float WmT_ChPFMet[2], WmT_PFMet[2];

  //! kinematics of the event: jets
  int theLeadingJet[2];
  int theSecondJet[2];
  float leadJetBtag[2], subleadJetBtag[2];

  //! kinematics of the event: met
  TVector3 *m_p3PFMet; 
  TVector3 m_p3ChPFMet[2];
  float modPFMet, modChPFMet[2];

  //! quark-gluon likelihood
  QGLikelihoodCalculator *qglikeli;
  
  //! vectors to store indices of best candidates
  std::vector<int> _acceptEleAll, _idEleAll, _isolEleAll, _convEleAll, _ipEleAll;
  std::vector<int> _acceptMuonsAll, _idMuonsAll, _isolMuonsAll, _ipMuonsAll;

  //! reduced tree for event selection (on at least 2leptons events)
  RedAxiTree *myOutTreeEle;  
  RedAxiTree *myOutTreeMu;

  //! variables for EleID
  int myRecoflag;
  float myPt, myEta, myPhi;
  float myDeta, myDphi, myHoe, mySee, myLh;
  float myTrackerIso, myHcalIso, myEcalIso, myCombinedIso;
  int myCharge;

  //! name of rootfile with dataset
  std::string _datasetName;  
};
#endif
