#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2F.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <iostream>

using namespace std;

int FRWeights = 0;

void addFRWeights() { FRWeights = 1;};
float getOfflineEff(float pT, float eta, TH2F *myH);

void addWeights(const char* filename, float baseW, int processId, int finalstate) {
  
  cout << "Adding weight branch to file " << filename << " with weight " << baseW << endl;

  TFile *fileOrig = 0;
  TTree *treeOrig = 0;

  fileOrig = TFile::Open(filename);
  if( fileOrig ) {
    fileOrig->cd();
    treeOrig = (TTree*)fileOrig->Get("T1");
  } else {
    cout << "File " << filename << " not existing !" << endl;
    return;
  }

  // reading root files with electrons and muons efficiencies
  TFile fileSFmuons42("/cmsrm/pc23_2/crovelli/data/muonTeP_LP/Muons_vpvPlusExpo_OutputScaleFactorMap_Summer11_42X.root");
  TH2F *histoSFmuons42 = (TH2F*)fileSFmuons42.Get("hScaleFactorMap");
  //
  TFile fileSFEle42("/cmsrm/pc23_2/crovelli/data/muonTeP_LP/Electrons_vpvPlusExpo_OutputScaleFactorMap_Summer11_42X.root");
  TH2F *histoSFele42 = (TH2F*)fileSFEle42.Get("hScaleFactorMap");
  
  if ( treeOrig ) {
    int nentriesOrig = treeOrig->GetEntries();

    TFile *fileNew = TFile::Open(filename,"recreate");
    TTree *treeNew = new TTree("t1new","tree with only selected events");

    std::vector<TTree*> trees; 
    trees.push_back(treeNew);
    
    // add also a branch with jet category (1 for njets=0, -1 for njets=1: useful for the fit)
    // and a branch with float final selection bool (for roofit)
    bool            hlt;
    Bool_t          step[17];
    Int_t           run;
    Int_t           ls;
    Int_t           event;
    Int_t           nVtx;
    Int_t           njets;
    Int_t           numExtraLep;
    Float_t         puweight;
    Float_t         leptPt1;
    Float_t         leptPt2;
    Float_t         leadingJetBTag;
    Float_t         secondJetBTag;
    Float_t         dijetInvMass;
    Float_t         dijetPt;
    Float_t         dijetDeta;
    Float_t         PFMet;
    Float_t         chPFMet;
    Float_t         PFWMt;
    Float_t         chPFWMt;
    /*
    Float_t         pt;
    Float_t         eta;
    Float_t         phi;
    Float_t         deta;
    Float_t         dphi;
    Float_t         hoe;
    Float_t         see;
    Float_t         trackerIso;
    Float_t         hcalIso;
    Float_t         ecalIso;
    Float_t         combinedIso;
    Int_t           charge;
    Float_t         lh;
    */
    Float_t         pxPFMet;
    Float_t         pyPFMet;
    Float_t         pzPFMet;
    Float_t         pxChPFMet;
    Float_t         pyChPFMet;
    Float_t         pzChPFMet;
    Float_t         pxLeadJet;
    Float_t         pyLeadJet;
    Float_t         pzLeadJet;
    Float_t         pxSecondJet;
    Float_t         pySecondJet;
    Float_t         pzSecondJet;
    Float_t         pxLept1;
    Float_t         pyLept1;
    Float_t         pzLept1;
    Float_t         pxLept2;
    Float_t         pyLept2;
    Float_t         pzLept2;
    Float_t         leadingJetLike;
    Float_t         secondJetLike;
    Float_t         productJetLike;

    treeOrig->SetBranchAddress("run", &run);
    treeOrig->SetBranchAddress("ls", &ls);
    treeOrig->SetBranchAddress("event", &event);
    treeOrig->SetBranchAddress("puweight", &puweight);
    treeOrig->SetBranchAddress("hlt", &hlt);
    treeOrig->SetBranchAddress("nVtx", &nVtx);
    treeOrig->SetBranchAddress("leptPt1", &leptPt1);
    treeOrig->SetBranchAddress("leptPt2", &leptPt2);
    treeOrig->SetBranchAddress("njets", &njets);
    treeOrig->SetBranchAddress("leadingJetBTag", &leadingJetBTag);
    treeOrig->SetBranchAddress("secondJetBTag", &secondJetBTag);
    treeOrig->SetBranchAddress("numExtraLep", &numExtraLep);
    treeOrig->SetBranchAddress("dijetInvMass", &dijetInvMass);
    treeOrig->SetBranchAddress("dijetPt", &dijetPt);
    treeOrig->SetBranchAddress("dijetDeta", &dijetDeta);
    treeOrig->SetBranchAddress("PFMet", &PFMet);
    treeOrig->SetBranchAddress("chPFMet", &chPFMet);
    treeOrig->SetBranchAddress("PFWMt", &PFWMt);
    treeOrig->SetBranchAddress("chPFWMt", &chPFWMt);
    treeOrig->SetBranchAddress("step", step);
    /*
    treeOrig->SetBranchAddress("pt", &pt);
    treeOrig->SetBranchAddress("eta", &eta);
    treeOrig->SetBranchAddress("phi", &phi);
    treeOrig->SetBranchAddress("deta", &deta);
    treeOrig->SetBranchAddress("dphi", &dphi);
    treeOrig->SetBranchAddress("hoe", &hoe);
    treeOrig->SetBranchAddress("see", &see);
    treeOrig->SetBranchAddress("trackerIso", &trackerIso);
    treeOrig->SetBranchAddress("ecalIso", &ecalIso);
    treeOrig->SetBranchAddress("hcalIso", &hcalIso);
    treeOrig->SetBranchAddress("combinedIso", &combinedIso);
    treeOrig->SetBranchAddress("charge", &charge);    
    treeOrig->SetBranchAddress("lh", &lh);
    */
    treeOrig->SetBranchAddress("pxPFMet", &pxPFMet);
    treeOrig->SetBranchAddress("pyPFMet", &pyPFMet);
    treeOrig->SetBranchAddress("pzPFMet", &pzPFMet);
    treeOrig->SetBranchAddress("pxChPFMet", &pxChPFMet);
    treeOrig->SetBranchAddress("pyChPFMet", &pyChPFMet);
    treeOrig->SetBranchAddress("pzChPFMet", &pzChPFMet);
    treeOrig->SetBranchAddress("pxLeadJet", &pxLeadJet);
    treeOrig->SetBranchAddress("pyLeadJet", &pyLeadJet);
    treeOrig->SetBranchAddress("pzLeadJet", &pzLeadJet);
    treeOrig->SetBranchAddress("pxSecondJet", &pxSecondJet);
    treeOrig->SetBranchAddress("pySecondJet", &pySecondJet);
    treeOrig->SetBranchAddress("pzSecondJet", &pzSecondJet);
    treeOrig->SetBranchAddress("pxLept1", &pxLept1);
    treeOrig->SetBranchAddress("pyLept1", &pyLept1);
    treeOrig->SetBranchAddress("pzLept1", &pzLept1);
    treeOrig->SetBranchAddress("pxLept2", &pxLept2);
    treeOrig->SetBranchAddress("pyLept2", &pyLept2);
    treeOrig->SetBranchAddress("pzLept2", &pzLept2);
    treeOrig->SetBranchAddress("leadingJetLike", &leadingJetLike);
    treeOrig->SetBranchAddress("secondJetLike",  &secondJetLike);
    treeOrig->SetBranchAddress("productJetLike", &productJetLike);

    // additional
    Float_t effW   = 1.0;   
    
    for(int i=0; i<(int)trees.size();i++) {
      TTree *theTreeNew = trees[i];

      // the selected final state: ele=0, mu=1
      theTreeNew->Branch("channel", &finalstate, "channel/I");

      // one integer containing the process identifier (for MC, 0 for data)
      theTreeNew->Branch("dataset", &processId, "dataset/I");

      // Copy branches
      theTreeNew->Branch("run", &run, "run/I");
      theTreeNew->Branch("ls", &ls, "ls/I");
      theTreeNew->Branch("event", &event, "event/I");
      theTreeNew->Branch("puW", &puweight, "puW/F");
      theTreeNew->Branch("effW", &effW, "effW/F");
      theTreeNew->Branch("pfmet", &PFMet, "pfmet/F");
      theTreeNew->Branch("chmet", &chPFMet, "chmet/F");  
      theTreeNew->Branch("pfwmt", &PFWMt, "pfwmt/F");
      theTreeNew->Branch("chwmt", &chPFWMt, "chwmt/F");  
      theTreeNew->Branch("trigger", &hlt, "trigger/O");
      theTreeNew->Branch("nvtx", &nVtx, "nvtx/I");
      theTreeNew->Branch("leptPt1", &leptPt1, "leptPt1/F");
      theTreeNew->Branch("leptPt2", &leptPt2, "leptPt2/F");
      // theTreeNew->Branch("leptEta", &eta, "leptEta/F");
      // theTreeNew->Branch("leptPhi", &phi, "leptphi/F");
      theTreeNew->Branch("njets", &njets, "njets/I");
      theTreeNew->Branch("leadingJetBTag", &leadingJetBTag, "leadingJetBTag/F");
      theTreeNew->Branch("secondJetBTag", &secondJetBTag, "secondJetBTag/F");
      theTreeNew->Branch("nextra", &numExtraLep, "nextra/I");
      theTreeNew->Branch("dijetInvMass", &dijetInvMass, "dijetInvMass/F");
      theTreeNew->Branch("dijetPt", &dijetPt, "dijetPt/F");
      theTreeNew->Branch("dijetDeta", &dijetDeta, "dijetDeta/F");
      theTreeNew->Branch("step", step, "step[25]/O");
      /*
      theTreeNew->Branch("deta_eleid", &deta, "deta_eleid/F");
      theTreeNew->Branch("dphi_eleid", &dphi, "dphi_eleid/F");
      theTreeNew->Branch("hoe_eleid", &hoe, "hoe_eleid/F");
      theTreeNew->Branch("see_eleid", &see, "see_eleid/F");
      theTreeNew->Branch("lh_eleid", &lh, "lh_eleid/F");
      theTreeNew->Branch("charge_eleid", &charge, "charge_eleid/I");
      theTreeNew->Branch("trackerIso", &trackerIso, "trackerIso/F");
      theTreeNew->Branch("ecalIso", &ecalIso, "ecalIso/F");
      theTreeNew->Branch("hcalIso", &hcalIso, "hcalIso/F");
      theTreeNew->Branch("combinedIso", &combinedIso, "combinedIso/F");
      */
      theTreeNew->Branch("pxPFMet", &pxPFMet, "pxPFMet/F");
      theTreeNew->Branch("pyPFMet", &pyPFMet, "pyPFMet/F");
      theTreeNew->Branch("pzPFMet", &pzPFMet, "pzPFMet/F");
      theTreeNew->Branch("pxChPFMet", &pxChPFMet, "pxChPFMet/F");
      theTreeNew->Branch("pyChPFMet", &pyChPFMet, "pyChPFMet/F");
      theTreeNew->Branch("pzChPFMet", &pzChPFMet, "pzChPFMet/F");
      theTreeNew->Branch("pxLeadJet", &pxLeadJet, "pxLeadJet[3]/F");
      theTreeNew->Branch("pyLeadJet", &pyLeadJet, "pyLeadJet[3]/F");
      theTreeNew->Branch("pzLeadJet", &pzLeadJet, "pzLeadJet[3]/F");
      theTreeNew->Branch("pxSecondJet", &pxSecondJet, "pxSecondJet[3]/F");
      theTreeNew->Branch("pySecondJet", &pySecondJet, "pySecondJet[3]/F");
      theTreeNew->Branch("pzSecondJet", &pzSecondJet, "pzSecondJet[3]/F");
      theTreeNew->Branch("pxLept1", &pxLept1, "pxLept1/F");
      theTreeNew->Branch("pyLept1", &pyLept1, "pyLept1/F");
      theTreeNew->Branch("pzLept1", &pzLept1, "pzLept1/F");
      theTreeNew->Branch("pxLept2", &pxLept2, "pxLept2/F");
      theTreeNew->Branch("pyLept2", &pyLept2, "pyLept2/F");
      theTreeNew->Branch("pzLept2", &pzLept2, "pzLept2/F");
      theTreeNew->Branch("leadingJetLike", &leadingJetLike, "leadingJetLike/F");
      theTreeNew->Branch("secondJetLike",  &secondJetLike,  "secondJetLike/F");
      theTreeNew->Branch("productJetLike", &productJetLike, "productJetLike/F");

      theTreeNew->Branch("baseW", &baseW,  "baseW/F");
    }

    int j =0;

    for(int i=0; i<nentriesOrig; i++) {
      if (i%10000 == 0) std::cout << ">>> Weighting event # " << i << " / " << nentriesOrig << " entries" << std::endl;
      treeOrig->GetEntry(i);
      
      if (processId>0) { // MC => apply scale factors
	if (finalstate==0) {   // mu
	  // effW = getOfflineEff(leptPt, eta, histoSFmuons42);    
	  effW = 1.;
	}
	else if (finalstate==1) { // ele
	  // effW = getOfflineEff(leptPt, eta, histoSFele42);
	  effW = 1.;
	}
      } else { // data
	effW = 1.;
      }
      
      if(processId < 100 || processId >= 1000 ) { // MC
	treeNew->Fill();
      } else { // data: apply the trigger 
	if(hlt) {
	  treeNew->Fill();
	}
      }
      j++;
    }
  
    fileNew->cd();
    treeNew->Write();
    fileNew->Close();

    fileOrig->cd();
    fileOrig->Close();

  } else {
    cout << "Tree T1 not present in the file " << filename << endl;
    return;
  }
}

float getOfflineEff(float pT, float eta, TH2F *myH) {

  float theEff=-1.;
  
  int   xBins = myH->GetXaxis()->GetNbins();
  float xMin  = myH->GetXaxis()->GetBinLowEdge(1);
  float xMax  = myH->GetXaxis()->GetBinUpEdge(xBins);
  int   yBins = myH->GetYaxis()->GetNbins();
  float yMin  = myH->GetYaxis()->GetBinLowEdge(1);
  float yMax  = myH->GetYaxis()->GetBinUpEdge(yBins);
  int theBin = myH->FindBin(pT, fabs(eta));
  if (pT>xMin && pT<xMax && fabs(eta)>yMin && fabs(eta)<yMax) {
    theEff = myH->GetBinContent(theBin);
  } else {
    theEff = 1.;
  }

  return theEff;
}


