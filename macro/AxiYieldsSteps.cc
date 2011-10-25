#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

float yieldErrPoisson(float nEst1, float n1, float nEst2=0, float n2=0, float nEst3=0, float n3=0, float nEst4=0, float n4=0, float nEst5=0, float n5=0, float nEst6=0, float n6=0);

using namespace std;

string fullSelCuts[19];

float Axi_fullSel[19];
float Wj_fullSel[19];
float SingleTop_fullSel[19];
float ttbar_fullSel[19];
float WW_fullSel[19];
float WZ_fullSel[19];
float ZZ_fullSel[19];
float gjets_fullSel[19];
float DY_fullSel[19];
float data_fullSel[19];

float Axi_fullSel_err[19];
float Wj_fullSel_err[19];
float SingleTop_fullSel_err[19];
float ttbar_fullSel_err[19];
float WW_fullSel_err[19];
float WZ_fullSel_err[19];
float ZZ_fullSel_err[19];
float gjets_fullSel_err[19];
float DY_fullSel_err[19];
float data_fullSel_err[19];

float Axi_eff_fullSel[19];
float Wj_eff_fullSel[19];
float SingleTop_eff_fullSel[19];
float ttbar_eff_fullSel[19];
float WW_eff_fullSel[19];
float WZ_eff_fullSel[19];
float ZZ_eff_fullSel[19];
float gjets_eff_fullSel[19];
float DY_eff_fullSel[19];

float Axi_finaleff_fullSel;    
float Wj_finaleff_fullSel;
float SingleTop_finaleff_fullSel;
float ttbar_finaleff_fullSel;
float WW_finaleff_fullSel;
float WZ_finaleff_fullSel;
float ZZ_finaleff_fullSel;
float gjets_finaleff_fullSel;
float DY_finaleff_fullSel;

float Axi_finaleff;    
float Wj_finaleff;
float SingleTop_finaleff;
float ttbar_finaleff;
float WW_finaleff;
float WZ_finaleff;
float ZZ_finaleff;
float gjets_finaleff;
float DY_finaleff;

float Axi_final[4][3];
float Wj_final[4][3];
float SingleTop_final[4][3];
float ttbar_final[4][3];
float WW_final[4][3];
float WZ_final[4][3];
float ZZ_final[4][3];
float gjets_final[4][3];
float DY_final[4][3];
float bkg_final[4][3];

float Axi_final_err[4][3];
float Wj_final_err[4][3];
float SingleTop_final_err[4][3];
float ttbar_final_err[4][3];
float WW_final_err[4][3];
float WZ_final_err[4][3];
float ZZ_final_err[4][3];
float gjets_final_err[4][3];
float DY_final_err[4][3];
float bkg_final_err[4][3];
float data_final[4][3];

// x-sections
std::map<int,float> axigluon_xsec_masses;

// backgrounds
float Wjet_xsec = 31314.;
float TTjets_xsec = 163.;
float SingleTopS_xsec   = 2.341;
float SingleTbarS_xsec  = 1.265;
float SingleTopT_xsec   = 35.72;
float SingleTbarT_xsec  = 18.43;
float SingleTopTW_xsec  = 7.46;
float SingleTbarTW_xsec = 7.46;
float WW_xsec = 47.;
float WZ_xsec = 18.2;
float ZZ_xsec = 7.41;   
float GJets_40_100_xsec = 25690.;
float GJets_100_200_xsec = 5213.;
float DY_xsec = 3048.;

string sampleNames[16];

void computeYields(float lumi=1, const char* finalstate="EE", int mass=0) {

  TChain *chains_fullSel[16];

  for(int isample=0; isample<16; isample++) {
    char fullsel_treename[200];
    sprintf(fullsel_treename,"FULL_SELECTION_EVENT_COUNTER_%s",finalstate);
    chains_fullSel[isample] = new TChain(fullsel_treename);
  }

  // signa x-secs
  axigluon_xsec_masses.insert(std::make_pair(150,77.3440));
  axigluon_xsec_masses.insert(std::make_pair(500,3.9928));
  axigluon_xsec_masses.insert(std::make_pair(1000,0.3164));
  axigluon_xsec_masses.insert(std::make_pair(1500,0.0382));
  
  // data
  sampleNames[0]  = "data";
  // backgrounds
  sampleNames[1]  = "W+jets";
  sampleNames[2]  = "TTbar+jets";
  sampleNames[3]  = "SingleTop_sChannel";
  sampleNames[4]  = "SingleTbar_sChannel";
  sampleNames[5]  = "SingleTop_tChannel";
  sampleNames[6]  = "SingleTbar_tChannel";
  sampleNames[7]  = "SingleTop_tWChannel";
  sampleNames[8]  = "SingleTbar_tWChannel";
  sampleNames[9]  = "WW";
  sampleNames[10] = "WZ";
  sampleNames[11] = "ZZ";
  sampleNames[12] = "GJets_40_100";
  sampleNames[13] = "GJets_100_200";
  sampleNames[14] = "DY";
  // signal
  sampleNames[15] = "AG+W->jj";

  float Axi_xsec;
  
  // mc
  char dir_mc[1000];                  
  sprintf(dir_mc,"results/Summer11_V1/");
  cout << "dir_mc = " << dir_mc << endl;

  // data
  char dir_data[1000]; 
  sprintf(dir_data,"results/Summer11_V1");
  cout << "dir_data = " << dir_data << endl;

  if(mass==0) { // use the default mass (150) to print the cut-by cut table

    int exampleMass=150;
    Axi_xsec = axigluon_xsec_masses[0];
    
    char AxiSample[500];
    sprintf(AxiSample,"AxigluonW_M-%d_TuneZ2_7TeV-calchep-pythia/*Counters.root",exampleMass);

    // backgrounds
    chains_fullSel[1]->Add(TString(dir_mc)+TString("/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir_mc)+TString("/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir_mc)+TString("/T_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir_mc)+TString("/T_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir_mc)+TString("/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root"));         // chiara" default?
    chains_fullSel[8]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root"));      // chiara" default?
    chains_fullSel[9]->Add(TString(dir_mc)+TString("/WW_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[10]->Add(TString(dir_mc)+TString("/WZ_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[11]->Add(TString(dir_mc)+TString("/ZZ_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[12]->Add(TString(dir_mc)+TString("/GJets_TuneZ2_40_HT_100_7TeV-madgraph/*Counters.root"));       
    chains_fullSel[13]->Add(TString(dir_mc)+TString("/GJets_TuneZ2_100_HT_200_7TeV-madgraph/*Counters.root"));       
    chains_fullSel[14]->Add(TString(dir_mc)+TString("/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/*Counters.root"));       
    // signal
    chains_fullSel[15]->Add(TString(dir_mc)+TString(AxiSample));       
    
    // data
    // if (strcmp(finalstate,"EE")==0) chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleElectron/*Counters.root"));
    // if (strcmp(finalstate,"MM")==0) {
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleMu/*Counters.root"));
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/SingleMu/*Counters.root"));
    // }
  } else {

    Axi_xsec = axigluon_xsec_masses[mass];    

    char AxiSample[500];
    sprintf(AxiSample,"AxigluonW_M-%d_TuneZ2_7TeV-calchep-pythia/*Counters.root",mass);

    // backgrounds
    chains_fullSel[1]->Add(TString(dir_mc)+TString("/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[2]->Add(TString(dir_mc)+TString("/TTJets_TuneZ2_7TeV-madgraph-tauola/*Counters.root"));       
    chains_fullSel[3]->Add(TString(dir_mc)+TString("/T_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[4]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_s-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[5]->Add(TString(dir_mc)+TString("/T_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[6]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_t-channel_7TeV-powheg-tauola/*Counters.root"));       
    chains_fullSel[7]->Add(TString(dir_mc)+TString("/T_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root"));         // chiara" default?
    chains_fullSel[8]->Add(TString(dir_mc)+TString("/Tbar_TuneZ2_tW-channel-DR_7TeV-powheg-tauola/*Counters.root"));      // chiara" default?
    chains_fullSel[9]->Add(TString(dir_mc)+TString("/WW_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[10]->Add(TString(dir_mc)+TString("/WZ_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[11]->Add(TString(dir_mc)+TString("/ZZ_TuneZ2_7TeV_pythia6_tauola/*Counters.root"));       
    chains_fullSel[12]->Add(TString(dir_mc)+TString("/GJets_TuneZ2_40_HT_100_7TeV-madgraph/*Counters.root"));       
    chains_fullSel[13]->Add(TString(dir_mc)+TString("/GJets_TuneZ2_100_HT_200_7TeV-madgraph/*Counters.root"));       
    chains_fullSel[14]->Add(TString(dir_mc)+TString("/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/*Counters.root"));       
    // signal
    chains_fullSel[15]->Add(TString(dir_mc)+TString(AxiSample));       

    // data
    // if (strcmp(finalstate,"EE")==0) chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleElectron/*Counters.root"));
    // if (strcmp(finalstate,"MM")==0) {
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleMu/*Counters.root"));
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/SingleMu/*Counters.root"));
    // }
  }

  float nFullSelTot[19][16];
  for(int isample=0; isample<16; isample++) {
    for(int icut=0; icut<19; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }
  
  // full selection
  int nCutsAnaFull = 19;
  for(int isample=0; isample<16; isample++) {

    // List of branches    
    Int_t           nCutsFull;
    Float_t         nSelFull[19];   //[nCuts]
    TBranch        *b_nCutsFull;    //!
    TBranch        *b_nSelFull;     //!
    chains_fullSel[isample]->SetBranchAddress("nCuts", &nCutsFull, &b_nCutsFull);
    chains_fullSel[isample]->SetBranchAddress("nSel",  nSelFull,   &b_nSelFull);
    
    Long64_t nentriesFull = chains_fullSel[isample]->GetEntries();

    // full selection
    for (Long64_t jentry=0; jentry<nentriesFull;jentry++) {
      Long64_t nb2;
      nb2 = chains_fullSel[isample]->GetEntry(jentry);   
      for(int icut=0; icut<nCutsFull; icut++) nFullSelTot[icut][isample] += nSelFull[icut];
    }
  }
  
  // eff at full selection level
  for(int icut=0; icut<nCutsAnaFull; icut++) {

    Axi_fullSel[icut] = lumi * Axi_xsec * nFullSelTot[icut][15]/nFullSelTot[0][15];
    Axi_fullSel_err[icut] = yieldErrPoisson(lumi * Axi_xsec * nFullSelTot[icut][15]/nFullSelTot[0][15], nFullSelTot[icut][15]);

    Wj_fullSel[icut] = lumi * Wjet_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1];
    Wj_fullSel_err[icut] = yieldErrPoisson(lumi * Wjet_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1], nFullSelTot[icut][1]);
    
    ttbar_fullSel[icut] = lumi * TTjets_xsec * nFullSelTot[icut][2]/nFullSelTot[0][2];
    ttbar_fullSel_err[icut] = yieldErrPoisson(lumi * TTjets_xsec * nFullSelTot[icut][2]/nFullSelTot[0][2], nFullSelTot[icut][2]);

    float singletop_tmp=0.;
    singletop_tmp += lumi * SingleTopS_xsec   * nFullSelTot[icut][3]/nFullSelTot[0][3];
    singletop_tmp += lumi * SingleTbarS_xsec  * nFullSelTot[icut][4]/nFullSelTot[0][4];
    singletop_tmp += lumi * SingleTopT_xsec   * nFullSelTot[icut][5]/nFullSelTot[0][5];
    singletop_tmp += lumi * SingleTbarT_xsec  * nFullSelTot[icut][6]/nFullSelTot[0][6];
    singletop_tmp += lumi * SingleTopTW_xsec  * nFullSelTot[icut][7]/nFullSelTot[0][7];
    singletop_tmp += lumi * SingleTbarTW_xsec * nFullSelTot[icut][8]/nFullSelTot[0][8];
    SingleTop_fullSel[icut] = singletop_tmp;

    SingleTop_fullSel_err[icut] = yieldErrPoisson(lumi * SingleTopS_xsec  * nFullSelTot[icut][3]/nFullSelTot[0][3], nFullSelTot[icut][3],
						  lumi * SingleTbarS_xsec * nFullSelTot[icut][4]/nFullSelTot[0][4], nFullSelTot[icut][4],
                                                  lumi * SingleTopT_xsec  * nFullSelTot[icut][5]/nFullSelTot[0][5], nFullSelTot[icut][5],
                                                  lumi * SingleTbarT_xsec * nFullSelTot[icut][6]/nFullSelTot[0][6], nFullSelTot[icut][6],
                                                  lumi * SingleTopTW_xsec * nFullSelTot[icut][7]/nFullSelTot[0][7], nFullSelTot[icut][7],
                                                  lumi * SingleTbarTW_xsec * nFullSelTot[icut][8]/nFullSelTot[0][8], nFullSelTot[icut][8]);

    WW_fullSel[icut] = lumi * WW_xsec * nFullSelTot[icut][9]/nFullSelTot[0][9];
    WW_fullSel_err[icut] = yieldErrPoisson(lumi * WW_xsec * nFullSelTot[icut][9]/nFullSelTot[0][9], nFullSelTot[icut][9]);
    
    WZ_fullSel[icut] = lumi * WZ_xsec * nFullSelTot[icut][10]/nFullSelTot[0][10];
    WZ_fullSel_err[icut] = yieldErrPoisson(lumi * WZ_xsec * nFullSelTot[icut][10]/nFullSelTot[0][10], nFullSelTot[icut][10]);

    ZZ_fullSel[icut] = lumi * ZZ_xsec * nFullSelTot[icut][11]/nFullSelTot[0][11];
    ZZ_fullSel_err[icut] = yieldErrPoisson(lumi * ZZ_xsec * nFullSelTot[icut][11]/nFullSelTot[0][11], nFullSelTot[icut][11]);
    
    float gjets_tmp=0.;
    gjets_tmp += lumi * GJets_40_100_xsec   * nFullSelTot[icut][12]/nFullSelTot[0][12];
    gjets_tmp += lumi * GJets_100_200_xsec  * nFullSelTot[icut][13]/nFullSelTot[0][13];
    gjets_fullSel[icut] = gjets_tmp;

    DY_fullSel[icut] = lumi * DY_xsec * nFullSelTot[icut][14]/nFullSelTot[0][14];
    DY_fullSel_err[icut] = yieldErrPoisson(lumi * DY_xsec * nFullSelTot[icut][14]/nFullSelTot[0][14], nFullSelTot[icut][14]);

    // data
    data_fullSel[icut] = nFullSelTot[icut][0];
    
    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][15]>0) Axi_eff_fullSel[icut] = Axi_fullSel[icut]/Axi_fullSel[icut-1];
    else Axi_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][1]>0) Wj_eff_fullSel[icut] = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else Wj_eff_fullSel[icut] = 0.0;
    
    if(icut>0 && nFullSelTot[icut-1][2]>0) ttbar_eff_fullSel[icut] = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else ttbar_eff_fullSel[icut] = 0.0;

    if(icut>0 && SingleTop_fullSel[icut-1]>0) SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[icut-1];
    else SingleTop_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][9]>0) WW_eff_fullSel[icut] = nFullSelTot[icut][9]/nFullSelTot[icut-1][9];
    else WW_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][10]>0) WZ_eff_fullSel[icut] = nFullSelTot[icut][10]/nFullSelTot[icut-1][10];
    else WZ_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][11]>0) ZZ_eff_fullSel[icut] = nFullSelTot[icut][11]/nFullSelTot[icut-1][11];
    else ZZ_eff_fullSel[icut] = 0.0;

    if(icut>0 && gjets_fullSel[icut-1]>0) gjets_eff_fullSel[icut] = gjets_fullSel[icut]/gjets_fullSel[icut-1];
    else gjets_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][14]>0) DY_eff_fullSel[icut] = nFullSelTot[icut][14]/nFullSelTot[icut-1][14];
    else DY_eff_fullSel[icut] = 0.0;

    if(icut==0) { 
      Axi_eff_fullSel[icut]       = Axi_fullSel[icut]/Axi_fullSel[0];
      Wj_eff_fullSel[icut]        = Wj_fullSel[icut]/Wj_fullSel[0];
      ttbar_eff_fullSel[icut]     = nFullSelTot[icut][2]/nFullSelTot[0][2];
      SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[0];
      WW_eff_fullSel[icut]        = WW_fullSel[icut]/WW_fullSel[0];      
      WZ_eff_fullSel[icut]        = WZ_fullSel[icut]/WZ_fullSel[0];      
      ZZ_eff_fullSel[icut]        = ZZ_fullSel[icut]/ZZ_fullSel[0];      
      gjets_eff_fullSel[icut]     = gjets_fullSel[icut]/gjets_fullSel[0];
      DY_eff_fullSel[icut]        = DY_fullSel[icut]/DY_fullSel[0];      
    }
  }
  
  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(Axi_fullSel[0]>0)  Axi_finaleff_fullSel  = Axi_fullSel[nCutsAnaFull-1]/Axi_fullSel[0];
  else Axi_finaleff_fullSel = 0.0;

  if(Wj_fullSel[0]>0) Wj_finaleff_fullSel = Wj_fullSel[nCutsAnaFull-1]/Wj_fullSel[0];
  else Wj_finaleff_fullSel = 0.0;

  if(nFullSelTot[0][2]>0) ttbar_finaleff_fullSel = nFullSelTot[nCutsAnaFull-1][2]/nFullSelTot[0][2];
  else ttbar_finaleff_fullSel = 0.0;

  if(SingleTop_fullSel[0]>0) SingleTop_finaleff_fullSel = SingleTop_fullSel[nCutsAnaFull-1]/SingleTop_fullSel[0];
  else SingleTop_finaleff_fullSel = 0.0;

  if(WW_fullSel[0]>0) WW_finaleff_fullSel = WW_fullSel[nCutsAnaFull-1]/WW_fullSel[0];
  else WW_finaleff_fullSel = 0.0;

  if(WZ_fullSel[0]>0) WZ_finaleff_fullSel = WZ_fullSel[nCutsAnaFull-1]/WZ_fullSel[0];
  else WZ_finaleff_fullSel = 0.0;

  if(ZZ_fullSel[0]>0) ZZ_finaleff_fullSel = ZZ_fullSel[nCutsAnaFull-1]/ZZ_fullSel[0];
  else ZZ_finaleff_fullSel = 0.0;

  if(gjets_fullSel[0]>0) gjets_finaleff_fullSel = gjets_fullSel[nCutsAnaFull-1]/gjets_fullSel[0];
  else gjets_finaleff_fullSel = 0.0;

  if(DY_fullSel[0]>0) DY_finaleff_fullSel = DY_fullSel[nCutsAnaFull-1]/DY_fullSel[0];
  else DY_finaleff_fullSel = 0.0;
}

void setupCuts() {

  fullSelCuts[0]="event";
  fullSelCuts[1]="MCtruth";
  fullSelCuts[2]="trigger";
  fullSelCuts[3]="channel preSel.";
  fullSelCuts[4]="lepton Id";
  fullSelCuts[5]="lepton isolation";
  fullSelCuts[6]="conv. rej.";
  fullSelCuts[7]="lepton d0";
  fullSelCuts[8]="lepton pT";
  fullSelCuts[9]="met";
  fullSelCuts[10]="W mT";
  fullSelCuts[11]="extra lepton veto";
  fullSelCuts[12]=">=2 jets";
  fullSelCuts[13]="dijetInvMass";
  fullSelCuts[14]="dijetPt";
  fullSelCuts[15]="dijetDeta";
  fullSelCuts[16]="leadJetBtag";
  fullSelCuts[17]="subleadJetBtag";
  fullSelCuts[18]="nSoftMuons";
}

void printLatex(float lumi, const char* finalstate,int mass) {

  setupCuts();

  if(strcmp(finalstate,"EE") && strcmp(finalstate,"MM")) {
    cout << "ERROR! finalstate must be one among EE/MM. Exiting..." << endl;
    return;
  } else {
    cout << " === NOW COMPUTING YIELDS FOR FINAL STATE: " << finalstate << " ===" << endl;
  }
  
  computeYields(lumi,finalstate,mass);
  
  char namefile[200];
  sprintf(namefile,"yields_byCut.tex");
  ofstream textfile;
  textfile.open(namefile, ios_base::app);
  textfile.precision(2);

  textfile << "\\documentclass{article}" << endl;
  textfile << "\\setlength\\textheight{9.8in}" << endl;
  textfile << "\\usepackage{rotating}" << endl;
  textfile << "\\begin{document}" << endl;
  textfile << "\\begin{table}[p]" << endl
	   << "\\begin{tiny}" << endl
	   << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  textfile << "selection & data & AG & W$(l\\nu)$+jets & $t\\bar{t}$ & single t \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  textfile << "\\hline" << endl;
  textfile << "\\hline" << endl;
  
  for(int icut=0; icut<19; icut++) {
    
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << data_fullSel[icut]       << "\t&\t"
	     << Axi_fullSel[icut]        << " (" << 100. * Axi_eff_fullSel[icut]        << "\\%)" << "\t&\t"
	     << Wj_fullSel[icut]         << " (" << 100. * Wj_eff_fullSel[icut]         << "\\%)" << "\t&\t"
	     << ttbar_fullSel[icut]      << " (" << 100. * ttbar_eff_fullSel[icut]      << "\\%)" << "\t&\t"
	     << SingleTop_fullSel[icut]  << " (" << 100. * SingleTop_eff_fullSel[icut]  << "\\%)" << "\t";
    textfile << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total fullselection "  << "\t&\t"
	   << data_fullSel[18]        << "\t&\t"
	   << Axi_fullSel[18]         << " (" << 100. * Axi_finaleff_fullSel       << "\\%)"    << "\t&\t"
	   << Wj_fullSel[18]          << " (" << 100. * Wj_finaleff_fullSel        << "\\%)"    << "\t&\t"
	   << ttbar_fullSel[18]       << " (" << 100. * ttbar_finaleff_fullSel     << "\\%)"    << "\t&\t"
	   << SingleTop_fullSel[18]   << " (" << 100. * SingleTop_finaleff_fullSel << "\\%)"    << "\t";
  textfile << "\t\\\\" << endl;
  
  textfile << "\\hline" << endl;
  textfile << "\\hline" << endl
	   << "\\end{tabular}" << endl
	   << "\\end{center}" << endl
	   << "\\end{tiny}" << endl
	   << "\\caption{Axigluon with $m_{AG}$ = "  << mass << " GeV/c$^2$. Breakdown of signal and backgrounds events in "
	   << lumi << " $pb^{-1}$ for " << finalstate << " final state.} " << endl 
	   << "\\end{table}" << endl;

  textfile << "\\begin{table}[p]" << endl
	   << "\\begin{tiny}" << endl
	   << "\\begin{center}" << endl;
  textfile << "\\begin{tabular}{|c|c|c|c|c|c|}" << endl;
  textfile << "\\hline" << endl;
  textfile << "selection & WW & WZ & ZZ & gamma+jets & DY \t\\\\" << endl;
  textfile << "\\hline" << endl; 
  textfile << "\\hline" << endl;
  textfile << "\\hline" << endl;
  
  for(int icut=0; icut<19; icut++) {
    
    textfile << fullSelCuts[icut] << "\t&\t";
    
    textfile << fixed
	     << WW_fullSel[icut]         << " (" << 100. * WW_eff_fullSel[icut]         << "\\%)" << "\t&\t"
	     << WZ_fullSel[icut]         << " (" << 100. * WZ_eff_fullSel[icut]         << "\\%)" << "\t&\t"
	     << ZZ_fullSel[icut]         << " (" << 100. * ZZ_eff_fullSel[icut]         << "\\%)" << "\t&\t"
	     << gjets_fullSel[icut]      << " (" << 100. * gjets_eff_fullSel[icut]      << "\\%)" << "\t&\t"
	     << DY_fullSel[icut]         << " (" << 100. * DY_eff_fullSel[icut]         << "\\%)" << "\t";
    textfile << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total fullselection "  << "\t&\t"
	   << WW_fullSel[18]          << " (" << 100. * WW_finaleff_fullSel        << "\\%)"    << "\t&\t"
	   << WZ_fullSel[18]          << " (" << 100. * WZ_finaleff_fullSel        << "\\%)"    << "\t&\t"
	   << ZZ_fullSel[18]          << " (" << 100. * ZZ_finaleff_fullSel        << "\\%)"    << "\t&\t"
	   << gjets_fullSel[18]       << " (" << 100. * gjets_finaleff_fullSel     << "\\%)"    << "\t&\t"
	   << DY_fullSel[18]          << " (" << 100. * DY_finaleff_fullSel        << "\\%)"    << "\t";
  textfile << "\t\\\\" << endl;
  
  textfile << "\\hline" << endl;
  textfile << "\\hline" << endl
	   << "\\end{tabular}" << endl
	   << "\\end{center}" << endl
	   << "\\end{tiny}" << endl
	   << "\\caption{Axigluon with $m_{AG}$ = "  << mass << " GeV/c$^2$. Breakdown of signal and backgrounds events in "
	   << lumi << " $pb^{-1}$ for " << finalstate << " final state.} " << endl 
	   << "\\end{table}" << endl;

  textfile << "\\end{document}" << endl;
}

float yieldErrPoisson(float nEst1, float n1, float nEst2, float n2, float nEst3, float n3, float nEst4, float n4, float nEst5, float n5, float nEst6, float n6) {

  float sum=0;
  if(n1>0) sum += pow(nEst1,2)/n1;
  if(n2>0) sum += pow(nEst2,2)/n2;
  if(n3>0) sum += pow(nEst3,2)/n3;
  if(n4>0) sum += pow(nEst4,2)/n4;
  if(n5>0) sum += pow(nEst5,2)/n5;
  if(n6>0) sum += pow(nEst6,2)/n6;
  
  return sqrt(sum);
}
