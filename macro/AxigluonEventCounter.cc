#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#define NSAMPLES 13 

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

  chains[8]->Add("results/Summer11_V1/WW_TuneZ2_7TeV_pythia6_tauola/*Counters.root");
  chains[9]->Add("results/Summer11_V1/WZ_TuneZ2_7TeV_pythia6_tauola/*Counters.root");

  chains[10]->Add("results/Summer11_V1/GJets_TuneZ2_40_HT_100_7TeV-madgraph/*Counters.root");
  chains[11]->Add("results/Summer11_V1/GJets_TuneZ2_100_HT_200_7TeV-madgraph/*Counters.root");

  chains[12]->Add("results/Summer11_V1/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/*Counters.root");
  
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

  sampleName.push_back("results/merged/WW_Ele.root");                   // 8
  sampleName.push_back("results/merged/WZ_Ele.root");                   // 9

  sampleName.push_back("results/merged/GJ_40-100_Ele.root");           // 10
  sampleName.push_back("results/merged/GJ_100-200_Ele.root");          // 11

  sampleName.push_back("results/merged/DY_Ele.root");               // 12

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
  // single-top, powheg samples. Xsecs taken from PREP: http://cms.cern.ch/iCMS/prep/requestmanagement?pwg=TOP&campid=Summer11
  sampleXsec.push_back(2.341);  // 1
  sampleXsec.push_back(1.265);  // 2
  sampleXsec.push_back(35.72);  // 3
  sampleXsec.push_back(18.43);  // 4
  sampleXsec.push_back(7.46);   // 5
  sampleXsec.push_back(7.46);   // 6
  // ttbar
  sampleXsec.push_back(163.);   // 7  - by A.Giammanco, private mail
  // dibosons
  sampleXsec.push_back(47.);    // 8 - s(NLO qqWW+ggWW) = 47 pb, gg/Tot = 0.0305 [K. Ellis]                             
  sampleXsec.push_back(18.2);   // 9 - from ?? (same as HWW)
  // gamma+jets
  sampleXsec.push_back(25690.);  // 10 - from prep
  sampleXsec.push_back(5213.);   // 11 - from prep
  // DY, m>50                        
  sampleXsec.push_back(3048.);   // 12 - from PGiulio

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
  sampleProcessId.push_back(10000); // WW
  sampleProcessId.push_back(10001); // WZ
  sampleProcessId.push_back(20000); // GJets 40-100
  sampleProcessId.push_back(20001); // GJets 100-200 
  sampleProcessId.push_back(30000); // DY, m>50

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

