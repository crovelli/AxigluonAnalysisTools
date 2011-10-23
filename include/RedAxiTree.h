#ifndef RedAxiTree_h
#define RedAxiTree_h

class TFile;
class TTree;

class G3EventProxy;

class RedAxiTree {
public:
   RedAxiTree(const char * filename = "axi.root");
  ~RedAxiTree();

  //! add the electron ID+iso variables for the selected best electrons
  void addElectronInfos();

  //! add steps
  void addSteps();
  
  //! add kinematics
  void addKinematics();

  //! event by event final dataset fill
  void fillAll(int nvtx, float lpt, int njets, float ljbt, float sjbt, int nsoftmu, int numExtraLep, float dijim, float dijpt, float djdeta, float pfmet, float chmet, float pfmt, float chmt, float ljl, float sjl, float pjl);
    
  // kinematics
  void fillKinematics(float pxPFMet, float pyPFMet, float pzPFMet, float pxChMet, float pyChMet, float pzChMet, float pxLeadJet, float pyLeadJet, float pzLeadJet, float pxSecJet, float pySecJet, float pzSecJet, float pxL, float pyL, float pzL);

  //! fill electron ID variables
  void fillElectrons(float pt, float eta, float phi, float deta, float dphi, float hoe, float see, float trackerIso, float hcalIso, float ecalIso, float combinedIso, int charge, float lh);

  //! fill the run,lumi, event number
  void fillRunInfos(int run, int lumi, int event, float puweight, bool HLT);   

  //! steps
  void fillSteps(bool s0, bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16);
    
  //! effectively store the events in the tree
  void store();

  //! save in the ROOT file
  void save();

private:

  // general
  bool myHLT;
  int myRun, myLS, myEvent;
  int myNVtx, myNjets, myNSoftMu, myNumExtraLep;
  float myPUWeight;
  float myLeptPt;
  float myLeadingJetBTag, mySecondJetBTag;
  float myDijetInvMass, myDijetPt, myDijetDeta;
  float myPFMet, myChPFMet, myPFWMt, myChPFWMt;
  float myLeadingJetLike, mySecondJetLike, myProductJetLike; 

  // steps
  bool mySteps[17]; 

  // electron variables
  float myPt, myEta, myPhi;
  float myDeta, myDphi, myHoe, mySee, myLh;
  float myTrackerIso, myHcalIso, myEcalIso, myCombinedIso;
  int myCharge;

  // kinematics
  float myPxPFMet,     myPyPFMet,      myPzPFMet;
  float myPxChPFMet,   myPyChPFMet,    myPzChPFMet;
  float myPxLeadJet,   myPyLeadJet,    myPzLeadJet;
  float myPxSecondJet, myPySecondJet,  myPzSecondJet;
  float myPxLept,      myPyLept,       myPzLept;

  TFile* myFile;
  TTree* myTree; 
};

#endif // RedAxiTree_h
