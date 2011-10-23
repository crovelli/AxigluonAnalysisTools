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
float data_fullSel[19];

float Axi_fullSel_err[19];
float Wj_fullSel_err[19];
float SingleTop_fullSel_err[19];
float ttbar_fullSel_err[19];
float data_fullSel_err[19];

float Axi_eff_fullSel[19];
float Wj_eff_fullSel[19];
float SingleTop_eff_fullSel[19];
float ttbar_eff_fullSel[19];

float Axi_finaleff_fullSel;    
float Wj_finaleff_fullSel;
float SingleTop_finaleff_fullSel;
float ttbar_finaleff_fullSel;

float Axi_finaleff;    
float Wj_finaleff;
float SingleTop_finaleff;
float ttbar_finaleff;

float Axi_final[4][3];
float Wj_final[4][3];
float SingleTop_final[4][3];
float ttbar_final[4][3];
float bkg_final[4][3];

float Axi_final_err[4][3];
float Wj_final_err[4][3];
float SingleTop_final_err[4][3];
float ttbar_final_err[4][3];
float bkg_final_err[4][3];
float data_final[4][3];

// x-sections
std::map<int,float> axigluon_xsec_masses;

// backgrounds
float Wjet_xsec = 31314.;
float TTjets_xsec = 157.5;
float SingleTopS_xsec   = 2.341;
float SingleTbarS_xsec  = 1.265;
float SingleTopT_xsec   = 3.572;
float SingleTbarT_xsec  = 1.843;
float SingleTopTW_xsec  = 7.87;
float SingleTbarTW_xsec = 7.87;

string sampleNames[10];

void computeYields(float lumi=1, const char* finalstate="EE", int mass=0) {

  TChain *chains_fullSel[10];

  for(int isample=0; isample<10; isample++) {
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
  sampleNames[0] = "data";
  // backgrounds
  sampleNames[1] = "W+jets";
  sampleNames[2] = "TTbar+jets";
  sampleNames[3] = "SingleTop_sChannel";
  sampleNames[4] = "SingleTbar_sChannel";
  sampleNames[5] = "SingleTop_tChannel";
  sampleNames[6] = "SingleTbar_tChannel";
  sampleNames[7] = "SingleTop_tWChannel";
  sampleNames[8] = "SingleTbar_tWChannel";
  // signal
  sampleNames[9] = "AG+W->jj";

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
    
    // signal
    chains_fullSel[9]->Add(TString(dir_mc)+TString(AxiSample));       
    
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
    
    // signal
    chains_fullSel[9]->Add(TString(dir_mc)+TString(AxiSample));       

    // data
    // if (strcmp(finalstate,"EE")==0) chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleElectron/*Counters.root"));
    // if (strcmp(finalstate,"MM")==0) {
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/DoubleMu/*Counters.root"));
    //  chains_fullSel[0]->Add(TString(dir_data)+TString("/SingleMu/*Counters.root"));
    // }
  }

  float nFullSelTot[19][10];
  for(int isample=0; isample<10; isample++) {
    for(int icut=0; icut<19; icut++) { nFullSelTot[icut][isample] = 0.0; }
  }
  
  // full selection
  int nCutsAnaFull = 19;
  for(int isample=0; isample<10; isample++) {

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

    Axi_fullSel[icut] = lumi * Axi_xsec * nFullSelTot[icut][9]/nFullSelTot[0][9];
    Axi_fullSel_err[icut] = yieldErrPoisson(lumi * Axi_xsec * nFullSelTot[icut][9]/nFullSelTot[0][9], nFullSelTot[icut][9]);

    Wj_fullSel[icut] = lumi * Wjet_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1];
    Wj_fullSel_err[icut] = yieldErrPoisson(lumi * Wjet_xsec * nFullSelTot[icut][1]/nFullSelTot[0][1], nFullSelTot[icut][1]);
    
    ttbar_fullSel[icut] = lumi * TTjets_xsec * nFullSelTot[icut][2]/nFullSelTot[0][2];
    ttbar_fullSel_err[icut] = yieldErrPoisson(lumi * TTjets_xsec * nFullSelTot[icut][2]/nFullSelTot[0][2], nFullSelTot[icut][2]);

    float singletop_tmp=0.;
    singletop_tmp += lumi * SingleTopS_xsec   * nFullSelTot[icut][2]/nFullSelTot[0][2];
    singletop_tmp += lumi * SingleTbarS_xsec  * nFullSelTot[icut][3]/nFullSelTot[0][3];
    singletop_tmp += lumi * SingleTopT_xsec   * nFullSelTot[icut][4]/nFullSelTot[0][4];
    singletop_tmp += lumi * SingleTbarT_xsec  * nFullSelTot[icut][5]/nFullSelTot[0][5];
    singletop_tmp += lumi * SingleTopTW_xsec  * nFullSelTot[icut][6]/nFullSelTot[0][6];
    singletop_tmp += lumi * SingleTbarTW_xsec * nFullSelTot[icut][7]/nFullSelTot[0][7];
    SingleTop_fullSel[icut] = singletop_tmp;

    SingleTop_fullSel_err[icut] = yieldErrPoisson(lumi * SingleTopS_xsec  * nFullSelTot[icut][2]/nFullSelTot[0][2], nFullSelTot[icut][2],
						  lumi * SingleTbarS_xsec * nFullSelTot[icut][3]/nFullSelTot[0][3], nFullSelTot[icut][3],
                                                  lumi * SingleTopT_xsec  * nFullSelTot[icut][4]/nFullSelTot[0][4], nFullSelTot[icut][4],
                                                  lumi * SingleTbarT_xsec * nFullSelTot[icut][5]/nFullSelTot[0][5], nFullSelTot[icut][5],
                                                  lumi * SingleTopTW_xsec * nFullSelTot[icut][6]/nFullSelTot[0][6], nFullSelTot[icut][6],
                                                  lumi * SingleTbarTW_xsec * nFullSelTot[icut][7]/nFullSelTot[0][7], nFullSelTot[icut][7]);

    // data
    data_fullSel[icut] = nFullSelTot[icut][0];
    
    // efficiencies
    if(icut>0 && nFullSelTot[icut-1][9]>0) Axi_eff_fullSel[icut] = Axi_fullSel[icut]/Axi_fullSel[icut-1];
    else Axi_eff_fullSel[icut] = 0.0;

    if(icut>0 && nFullSelTot[icut-1][1]>0) Wj_eff_fullSel[icut] = nFullSelTot[icut][1]/nFullSelTot[icut-1][1];
    else Wj_eff_fullSel[icut] = 0.0;
    
    if(icut>0 && nFullSelTot[icut-1][2]>0) ttbar_eff_fullSel[icut] = nFullSelTot[icut][2]/nFullSelTot[icut-1][2];
    else ttbar_eff_fullSel[icut] = 0.0;

    if(icut>0 && SingleTop_fullSel[icut-1]>0) SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[icut-1];
    else SingleTop_eff_fullSel[icut] = 0.0;

    if(icut==0) { 
      Axi_eff_fullSel[icut]       = Axi_fullSel[icut]/Axi_fullSel[0];
      Wj_eff_fullSel[icut]        = Wj_fullSel[icut]/Wj_fullSel[0];
      ttbar_eff_fullSel[icut]     = nFullSelTot[icut][2]/nFullSelTot[0][2];
      SingleTop_eff_fullSel[icut] = SingleTop_fullSel[icut]/SingleTop_fullSel[0];
    }
  }
  

  // final efficiency after full selections (-4 = 3 x jets + 1=final)
  if(nFullSelTot[0][9]>0)  Axi_finaleff_fullSel  = Axi_fullSel[nCutsAnaFull]/Axi_fullSel[0];
  else Axi_finaleff_fullSel = 0.0;

  if(Wj_fullSel[0]>0) Wj_finaleff_fullSel = Wj_fullSel[nCutsAnaFull]/Wj_fullSel[0];
  else Wj_finaleff_fullSel = 0.0;

  if(nFullSelTot[0][2]>0) ttbar_finaleff_fullSel = nFullSelTot[nCutsAnaFull][2]/nFullSelTot[0][2];
  else ttbar_finaleff_fullSel = 0.0;

  if(SingleTop_fullSel[0]>0) SingleTop_finaleff_fullSel = SingleTop_fullSel[nCutsAnaFull]/SingleTop_fullSel[0];
  else SingleTop_finaleff_fullSel = 0.0;
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
  fullSelCuts[12]=">2 jets";
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
	     << data_fullSel[icut]   << "\t&\t"
	     << Axi_fullSel[icut]    << " (" << 100. * Axi_eff_fullSel[icut]   << "\\%)" << "\t&\t"
	     << Wj_fullSel[icut]     << " (" << 100. * Wj_eff_fullSel[icut]  << "\\%)" << "\t&\t"
	     << ttbar_fullSel[icut]    << " (" << 100. * ttbar_eff_fullSel[icut] << "\\%)" << "\t&\t"
	     << SingleTop_fullSel[icut]    << " (" << 100. * SingleTop_eff_fullSel[icut] << "\\%)" << "\t";
    textfile << "\t\\\\" << endl;
  }
  
  textfile << "\\hline" << endl;
  
  textfile << "total fullselection " << "\t&\t"
	   << data_fullSel[18]   << "\t&\t"
	   << Axi_fullSel[18]    << " (" << 100. * Axi_finaleff_fullSel  << "\\%)"     << "\t&\t"
	   << Wj_fullSel[18]     << " (" << 100. * Wj_finaleff_fullSel << "\\%)"     << "\t&\t"
	   << ttbar_fullSel[18]    << " (" << 100. * ttbar_finaleff_fullSel << "\\%)"    << "\t&\t"
	   << SingleTop_fullSel[18]    << " (" << 100. * SingleTop_finaleff_fullSel << "\\%)"    << "\t";
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
