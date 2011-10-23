#ifndef CutBasedAxiSelector_h
#define CutBasedAxiSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CutBasedAxiSelector {

public:

  //! constructor
  CutBasedAxiSelector();

  //! copy constructor
  CutBasedAxiSelector( const CutBasedAxiSelector& selector );

  //! destructor
  virtual ~CutBasedAxiSelector();   

  //! configure from files
  void Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle);

  //! get the applied selection
  Selection* GetSelection() { return _selection; }

  //! set event by event observables
  void SetProcessID(int processID)    { m_processID = processID; }
  void SetMcTruth(bool mctruth)       { m_foundMcTree = mctruth; }
  void SetWeight(float weight)        { m_weight = weight; }
  void SetHLT(bool passedHLT)         { m_passedHLT = passedHLT; }
  void SetIsChannel(bool channelsel)  { m_isThisChannel = channelsel; }
  void SetElectronId(int isEleId)     { m_isElectronId  = isEleId; }
  void SetElectronIsolation(int isEleIsol)  { m_isElectronIsol = isEleIsol; }
  void SetElectronConvRejection(int isEleConvRej) { m_isElectronConvRej = isEleConvRej; }
  void SetElectronIp(int isEleIp) { m_isElectronIp = isEleIp; }
  void SetLeptPt(float leptPt) { m_leptPt = leptPt; }
  void SetNExtraLeptons(int nextralep) { m_nExtraLeptons = nextralep; }
  void SetMet(float met) { m_met = met;}
  void SetWMt(float wmt) { m_wmt = wmt;}
  void SetNJets(int njets) { m_nJets = njets;}
  void SetDijetMass(float dijInvMass)  { m_dijetInvMass = dijInvMass;}
  void SetDijetPt(float dijPt)  { m_dijetPt = dijPt;}          
  void SetDijetDeltaEta(float dijDe)  { m_dijetDeta = dijDe;}  
  void SetBTagLeadJet(float btlead) { m_leadJetBtag = btlead;}  
  void SetBTagSubleadJet(float btslead) { m_subleadJetBtag = btslead;}  
  void SetNSoftMuons(int nsoftmu) { m_nSoftMuons    = nsoftmu; }

  //! get output of the selector
  bool output();

  //! steps 
  bool outputStep0() { return m_step0; }
  bool outputStep1() { return m_step1; }
  bool outputStep2() { return m_step2; }
  bool outputStep3() { return m_step3; }
  bool outputStep4() { return m_step4; }
  bool outputStep5() { return m_step5; }
  bool outputStep6() { return m_step6; }
  bool outputStep7() { return m_step7; }
  bool outputStep8() { return m_step8; }
  bool outputStep9() { return m_step9; }
  bool outputStep10() { return m_step10; }
  bool outputStep11() { return m_step11; }
  bool outputStep12() { return m_step12; }
  bool outputStep13() { return m_step13; }
  bool outputStep14() { return m_step14; }
  bool outputStep15() { return m_step15; }
  bool outputStep16() { return m_step16; }

  //! display the electron efficiency
  void displayEfficiencies(std::string datasetName);
  
private:
  
  float m_weight;
  bool m_foundMcTree;
  bool m_passedHLT;
  bool m_isThisChannel;
  int m_isElectronId, m_isElectronIsol, m_isElectronConvRej, m_isElectronIp;
  int m_nJets, m_nSoftMuons, m_nExtraLeptons;
  float m_leptPt, m_met, m_wmt;
  float m_dijetInvMass, m_dijetPt, m_dijetDeta;
  float m_leadJetBtag, m_subleadJetBtag;
  int m_processID;

  //! contains the preselection cuts
  Selection* _selection;

  //! counters for the efficiencies display, based on electron candidates
  Counters* globalCounter;

  //! different selection steps
  bool m_step0, m_step1, m_step2, m_step3, m_step4, m_step5, m_step6, m_step7, m_step8, m_step9;  
  bool m_step10, m_step11, m_step12, m_step13, m_step14, m_step15, m_step16;  

  //! this is to do an efficiency for each process in the sample 
  //! (if more than one is present)
  //! to turn on it, use SetProcessID(int processID) with processID=!-1
  std::map<int, Counters*> multiProcessCounter;

};

#endif
