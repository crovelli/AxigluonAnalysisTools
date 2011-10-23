#include "AxigluonAnalysisTools/include/RedAxiTree.h"

// C++
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

// Root
#include "TFile.h"
#include "TTree.h"

RedAxiTree::RedAxiTree(const char * filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","axi tree");

  // GENERAL block
  myTree->Branch("run",                 &myRun,                 "run/I");
  myTree->Branch("ls",                  &myLS,                  "ls/I");
  myTree->Branch("event",               &myEvent,               "event/I");
  myTree->Branch("puweight",            &myPUWeight,            "puweight/F");
  myTree->Branch("hlt",                 &myHLT,                 "hlt/O");
  myTree->Branch("nVtx",                &myNVtx,                "nVtx/I");
  myTree->Branch("leptPt",              &myLeptPt,              "leptPt/F");   
  myTree->Branch("njets",               &myNjets,               "njets/I");
  myTree->Branch("leadingJetBTag",      &myLeadingJetBTag,      "leadingJetBTag/F");  
  myTree->Branch("secondJetBTag",       &mySecondJetBTag,       "secondJetBTag/F");   
  myTree->Branch("nSoftMu",             &myNSoftMu,             "nSoftMu/I");
  myTree->Branch("numExtraLep",         &myNumExtraLep,         "numExtraLep/I");   
  myTree->Branch("dijetInvMass",        &myDijetInvMass,        "dijetInvMass/F");  
  myTree->Branch("dijetPt",             &myDijetPt,             "dijetPt/F");  
  myTree->Branch("dijetDeta",           &myDijetDeta,           "dijetDeta/F");  
  myTree->Branch("PFMet",               &myPFMet,               "PFMet/F");    
  myTree->Branch("chPFMet",             &myChPFMet,             "chPFMet/F");
  myTree->Branch("PFWMt",               &myPFWMt,               "PFWMt/F");    
  myTree->Branch("chPFWMt",             &myChPFWMt,             "chPFWMt/F");  
  myTree->Branch("leadingJetLike",      &myLeadingJetLike,      "leadingJetLike/F");  
  myTree->Branch("secondJetLike",       &mySecondJetLike,       "secondJetLike/F");   
  myTree->Branch("productJetLike",      &myProductJetLike,      "productJetLike/F");   
}

RedAxiTree::~RedAxiTree() {

  delete myFile;
}

void RedAxiTree::addElectronInfos() {
  
  myTree->Branch("pt",   &myPt,   "pt/F");
  myTree->Branch("eta",  &myEta,  "eta/F");
  myTree->Branch("phi",  &myPhi,  "phi/F");
  myTree->Branch("deta", &myDeta, "deta/F");
  myTree->Branch("dphi", &myDphi, "dphi/F");
  myTree->Branch("hoe",  &myHoe,  "hoe/F");
  myTree->Branch("see",  &mySee,  "see/F");
  myTree->Branch("trackerIso",  &myTrackerIso,  "trackerIso/F");
  myTree->Branch("hcalIso",     &myHcalIso,     "hcalIso/F");
  myTree->Branch("ecalIso",     &myEcalIso,     "ecalIso/F");
  myTree->Branch("combinedIso", &myCombinedIso, "combinedIso/F");
  myTree->Branch("charge",  &myCharge, "charge/I");
  myTree->Branch("lh",      &myLh,     "lh/F");
}

void RedAxiTree::addSteps() {

  myTree->Branch("step", mySteps, "step[17]/O"); 
}

void RedAxiTree::addKinematics() {

  myTree->Branch("pxPFMet",     &myPxPFMet,    "pxPFMet/F");
  myTree->Branch("pyPFMet",     &myPyPFMet,    "pyPFMet/F");
  myTree->Branch("pzPFMet",     &myPzPFMet,    "pzPFMet/F");
  myTree->Branch("pxChPFMet",   &myPxChPFMet,  "pxChPFMet/F");
  myTree->Branch("pyChPFMet",   &myPyChPFMet,  "pyChPFMet/F");
  myTree->Branch("pzChPFMet",   &myPzChPFMet,  "pzChPFMet/F");
  myTree->Branch("pxLeadJet",   &myPxLeadJet,  "pxLeadJet/F"); 
  myTree->Branch("pyLeadJet",   &myPyLeadJet,  "pyLeadJet/F");
  myTree->Branch("pzLeadJet",   &myPzLeadJet,  "pzLeadJet/F");
  myTree->Branch("pxSecondJet", &myPxSecondJet,"pxSecondJet/F");
  myTree->Branch("pySecondJet", &myPySecondJet,"pySecondJet/F");
  myTree->Branch("pzSecondJet", &myPzSecondJet,"pzSecondJet/F");
  myTree->Branch("pxLept",      &myPxLept,     "pxLept/F");
  myTree->Branch("pyLept",      &myPyLept,     "pyLept/F");
  myTree->Branch("pzLept",      &myPzLept,     "pzLept/F");
}

void RedAxiTree::store() {

  myTree->Fill();
}


void RedAxiTree::save() {

  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void RedAxiTree::fillAll(int nvtx, float lpt, int njets, float ljbt, float sjbt, int nsmu, int nextral, float dijim, float dijpt, float dijdeta, float pfm, float chm, float pfmt, float chmt, float ljl, float sjl, float pjl) {

  myNVtx = nvtx;
  myLeptPt = lpt; 
  myNjets  = njets;
  myLeadingJetBTag = ljbt;
  mySecondJetBTag  = sjbt;
  myNSoftMu      = nsmu;
  myNumExtraLep  = nextral;
  myDijetInvMass = dijim;
  myDijetPt      = dijpt;
  myDijetDeta    = dijdeta;
  myPFMet   = pfm;
  myChPFMet = chm;          
  myPFWMt   = pfmt; 
  myChPFWMt = chmt;   
  myLeadingJetLike = ljl;
  mySecondJetLike  = sjl;
  myProductJetLike = pjl;
}

void RedAxiTree::fillSteps(bool s0, bool s1, bool s2, bool s3, bool s4, bool s5, bool s6, bool s7, bool s8, bool s9, bool s10, bool s11, bool s12, bool s13, bool s14, bool s15, bool s16) {

  mySteps[0]  = s0;
  mySteps[1]  = s1;
  mySteps[2]  = s2;
  mySteps[3]  = s3;
  mySteps[4]  = s4;
  mySteps[5]  = s5;
  mySteps[6]  = s6;
  mySteps[7]  = s7;
  mySteps[8]  = s8;
  mySteps[9]  = s9;
  mySteps[10] = s10;
  mySteps[11] = s11;
  mySteps[12] = s12;
  mySteps[13] = s13;
  mySteps[14] = s14;
  mySteps[15] = s15;
  mySteps[16] = s16;
}

void RedAxiTree::fillKinematics(float pxPFMet, float pyPFMet, float pzPFMet,
				float pxChMet, float pyChMet, float pzChMet,
				float pxLeadJet, float pyLeadJet, float pzLeadJet, 
				float pxSecJet,  float pySecJet, float pzSecJet,    
				float pxL, float pyL, float pzL) {

  myPxPFMet   = pxPFMet;
  myPyPFMet   = pyPFMet;
  myPzPFMet   = pzPFMet;
  myPxChPFMet = pxChMet;
  myPyChPFMet = pyChMet;
  myPzChPFMet = pzChMet;
  myPxLeadJet   = pxLeadJet;
  myPyLeadJet   = pyLeadJet;
  myPzLeadJet   = pzLeadJet;
  myPxSecondJet = pxSecJet;
  myPySecondJet = pySecJet;
  myPzSecondJet = pzSecJet;
  myPxLept = pxL;
  myPyLept = pyL;
  myPzLept = pzL;
}

void RedAxiTree::fillElectrons(float pt, float eta, float phi, float deta, float dphi, float hoe, float see, float trackerIso, float hcalIso, float ecalIso, float combinedIso, int charge,float lh) {

  myPt   = pt;
  myEta  = eta;
  myPhi  = phi;
  myDeta = deta;
  myDphi = dphi;
  myHoe  = hoe;
  mySee  = see;
  myTrackerIso = trackerIso;
  myHcalIso = hcalIso;
  myEcalIso = ecalIso;
  myCombinedIso = combinedIso;
  myCharge = charge;
  myLh = lh;
}

void RedAxiTree::fillRunInfos(int run, int lumi, int event, float puweight, bool HLT) {

  myRun   = run;
  myLS    = lumi;
  myEvent = event;
  myPUWeight = puweight;
  myHLT = HLT;
}

