#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#define NSAMPLES 8 

using namespace std;

double weight(double ngen, double xsec, double filtereff, double lumi = 1);

void countEvents() {

  char nametree[200];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EE");

  cout << "nametree = " << nametree << endl;

  // signals
  TChain *signalChains[4];
  int mAxi[4] = {150,500,1000,1500};
  for(int imass=0; imass<4;imass++) {
    char mass[5];
    sprintf(mass,"%d",mAxi[imass]);
    signalChains[imass] = new TChain(nametree);

    TString hSample("results/Summer11_V1/AxigluonW_M-");
    hSample += TString(mass)+TString("_TuneZ2_7TeV-calchep-pythia/*Counters.root");
    signalChains[imass]->Add(hSample.Data());
  }

  // backgrounds
  TChain *chains[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    chains[isample] = new TChain(nametree);
  }

  // nominal sample first, then the systematics ones
  chains[0]->Add("results/Summer11_V1/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root");

  chains[1]->Add("results/Summer11_V1/T_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root");
  chains[2]->Add("results/Summer11_V1/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root");
  chains[3]->Add("results/Summer11_V1/T_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root");
  chains[4]->Add("results/Summer11_V1/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root");
  chains[5]->Add("results/Summer11_V1/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root");
  chains[6]->Add("results/Summer11_V1/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root");

  chains[7]->Add("results/Summer11_V1/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root");

  cout << "chains added. " << endl;

  std::vector<TString>  signalSampleName;
  for(int imass=0; imass<4;imass++) {
    char mass[5];
    sprintf(mass,"%d",mAxi[imass]);
    signalSampleName.push_back(TString("results/merged/axiW")+TString(mass)+TString("_Ele.root"));  
  }

  std::vector<TString> sampleName;  
  sampleName.push_back("results/merged/Wjets_Ele.root");                // 0

  sampleName.push_back("results/merged/SingleT_sChannel_Ele.root");     // 1 
  sampleName.push_back("results/merged/SingleTbar_sChannel_Ele.root");  // 2
  sampleName.push_back("results/merged/SingleT_tChannel_Ele.root");     // 3
  sampleName.push_back("results/merged/SingleTbar_tChannel_Ele.root");  // 4
  sampleName.push_back("results/merged/SingleT_tWChannel_Ele.root");    // 5
  sampleName.push_back("results/merged/SingleTbar_tWChannel_Ele.root"); // 6

  sampleName.push_back("results/merged/TTbar_Ele.root");                // 7


  std::map<int,float> axigluon_xsec;  
  axigluon_xsec.insert(std::make_pair(150,77.3440));
  axigluon_xsec.insert(std::make_pair(500,3.9928));
  axigluon_xsec.insert(std::make_pair(1000,0.3164));
  axigluon_xsec.insert(std::make_pair(1500,0.0382));

  std::vector<double> signalXSec;
  for(int imass=0; imass<4;imass++) {
    double xsec = axigluon_xsec[mAxi[imass]]; 
    signalXSec.push_back(xsec);
  }

  std::vector<float> sampleXsec;
  // wjets
  sampleXsec.push_back(31314.); // madgraph // 0
  // powheg samples. Xsecs taken from PREP: http://cms.cern.ch/iCMS/prep/requestmanagement?pwg=TOP&campid=Summer11
  sampleXsec.push_back(2.341); // 1
  sampleXsec.push_back(1.265); // 2
  sampleXsec.push_back(3.572); // 3
  sampleXsec.push_back(1.843); // 4
  // prep says 7.46 pb. 7.87 is from Guillelmo
  sampleXsec.push_back(7.87); // 5
  sampleXsec.push_back(7.87); // 6
  // ttbar
  sampleXsec.push_back(157.5); // 7

  std::vector<double> signalProcId;
  for(int imass=0; imass<4;imass++) {
    signalProcId.push_back(1000+mAxi[imass]);
  }

  std::vector<int> sampleProcessId; 
  sampleProcessId.push_back(80); // Wjets
  sampleProcessId.push_back(15); // t(s-cha) 
  sampleProcessId.push_back(16); // tbar(s-cha) 
  sampleProcessId.push_back(13); // t(t-cha)
  sampleProcessId.push_back(14); // tbar(t-cha)
  sampleProcessId.push_back(19); // tW
  sampleProcessId.push_back(20); // tbarW
  sampleProcessId.push_back(10); // ttbar

  float nEvH[4];
  for(int imass=0; imass<4; imass++) nEvH[imass] = 0.0;
  for(int imass=0; imass<4; imass++) {

    cout << "\tProcessing signal sample mass # " << mAxi[imass] << "..." << endl;
      
    Int_t   nCuts;
    Float_t nSel[17];   //[nCuts]   
      
    // List of branches
    TBranch *b_nCuts;   //!
    TBranch *b_nSel;   //!
    signalChains[imass]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    signalChains[imass]->SetBranchAddress("nSel",  nSel,   &b_nSel);
    
    Long64_t nentries = signalChains[imass]->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = signalChains[imass]->GetEntry(jentry);   nbytes += nb;
      nEvH[imass] += nSel[0];
    }
  }

  // backgrounds
  float nEv[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) nEv[isample] = 0.0;

  for(int isample=0; isample<NSAMPLES; isample++) {

    Int_t   nCuts;
    Float_t nSel[25];   //[nCuts]
    
    // List of branches
    TBranch *b_nCuts;   //!
    TBranch *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;
      nEv[isample] += nSel[0];
    }
  }


  if (sampleXsec.size() != sampleName.size() ) cout << "nasty error! check sizes..." << endl;

  std::ofstream weightsFile;
  weightsFile.open("weightTrees.sh");
  weightsFile << "#! /bin/sh\n\n" << std::endl;
  weightsFile << "lumiEE=$1" << std::endl;
  weightsFile << "lumiMM=$2" << std::endl;
  
  // now write ele
  weightsFile << "echo \"Adding weights for ele datasets for \" $lumiEE \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<4; imass++) {
    double massXsec = signalXSec[imass];
    TString massSampleName = signalSampleName[imass];
    double massId = signalProcId[imass];
    float w = weight(nEvH[imass], massXsec, 1., 1.);
    weightsFile << "addWeights(\"" << massSampleName.Data() << "\", " << w << "*$lumiEE, " << massId << " ,1);" << std::endl;
  }

  for(int isample=0; isample<NSAMPLES; isample++) {
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    weightsFile << "addWeights(\"" << sampleName[isample].Data() << "\", " << w << "*$lumiEE, " << sampleProcessId[isample] << " ,1);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;
  

  // now write mu
  weightsFile << "echo \"Adding weights for muon datasets for \" $lumiMM \" pb-1...\"" << std::endl;
  weightsFile << "root -l -b <<EOF" << std::endl;
  weightsFile << ".L addWeightsToTree.cc+" << std::endl;
  for(int imass=0; imass<4; imass++) {
    double massXsec = signalXSec[imass];
    TString massSampleName = signalSampleName[imass];
    double massId = signalProcId[imass];
    cout << "Events processed for sample: " << massSampleName << " = " << nEvH[imass] << endl;
    float w = weight(nEvH[imass], massXsec, 1., 1.);
    TString massSampleNameMM = massSampleName.ReplaceAll("_Ele","_Mu");
    weightsFile << "addWeights(\"" << massSampleNameMM.Data() << "\", " << w << "*$lumiMM, " << massId << " ,0);" << std::endl;
  }
  for(int isample=0; isample<NSAMPLES; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMM = sampleName[isample].ReplaceAll("_Ele","_Mu");
    weightsFile << "addWeights(\"" << sampleNameMM.Data() << "\", " << w << "*$lumiMM, " << sampleProcessId[isample] << " ,0);" << std::endl;
  }
  weightsFile << ".q\n\nEOF\n" << std::endl;
  
  weightsFile << "echo \"done weighting.\"" << std::endl;
}

double weight(double ngen, double xsec, double filtereff, double lumi) {

  if(ngen==0) return 0;
  return xsec * filtereff * lumi / ngen;

}

