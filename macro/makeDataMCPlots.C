#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"

#include <iostream>

#define NSPECIES 4
#define NVARIABLES 10
#define NCUTS 1

void makeDataMCPlots()
{
  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x tdrstyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetOptTitle(0); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);

  TString suffix="";

  TString species[NSPECIES];
  species[0]="Data";
  species[1]="W+jets";
  species[2]="top";
  species[3]="AG+W";

  Color_t colors[NSPECIES];
  colors[0]=kBlack;
  colors[1]=kViolet;
  colors[2]=kOrange+8;
  colors[3]=kOrange;

  Color_t lineColors[NSPECIES];
  lineColors[0]=kBlack;
  lineColors[1]=kViolet+3;
  lineColors[2]=kOrange+4;
  lineColors[3]=kOrange+7;

  int legendOrder[NSPECIES];
  legendOrder[0]=0;
  legendOrder[1]=3;
  legendOrder[2]=2;
  legendOrder[3]=1;

  TString files[NSPECIES];
  // files[0]="results_data/datasets_trees/data_Wenu.root";
  files[1]="results/datasets_trees/Wjets_ll.root";
  files[2]="results/datasets_trees/top_ll.root";
  files[3]="results/datasets_trees/axiW500_ll.root";

  TString plotsDir="./axi/";
  TFile* fOut=new TFile("histos_"+suffix+".root","RECREATE");
  
  char icut[NCUTS][100];
  TH1F* histos[NSPECIES][NCUTS][NVARIABLES];   // 4 species, 1 cut levels, 10 variables
  
  TString variables[NVARIABLES];
  variables[0]="pfmet";
  variables[1]="chmet";
  variables[2]="pfwmt";
  variables[3]="chwmt";
  variables[4]="leptPt";
  variables[5]="dijetInvMass";
  variables[6]="dijetDeta";
  variables[7]="dijetPt";
  variables[8]="njets";
  variables[9]="nvtx";

  TString units[NVARIABLES];
  units[0]="GeV";
  units[1]="GeV";
  units[2]="GeV";
  units[3]="GeV";
  units[4]="GeV";
  units[5]="GeV";
  units[6]="";
  units[7]="GeV";
  units[8]="";
  units[9]="";

  int nbins[NVARIABLES];
  nbins[0]=40;
  nbins[1]=40;
  nbins[2]=50;  
  nbins[3]=50;  
  nbins[4]=75;  
  nbins[5]=75;  
  nbins[6]=50;  
  nbins[7]=50;
  nbins[8]=10;
  nbins[9]=25;

  float range[NVARIABLES][2]; // 8 variables, min, max
  // pf met
  range[0][0]=0.;
  range[0][1]=160.;
  // ch met
  range[1][0]=0.;
  range[1][1]=160.;
  // pf mT
  range[2][0]=0.;    
  range[2][1]=200.;   
  // ch mT
  range[3][0]=0.;    
  range[3][1]=200.;   
  // lepton pt
  range[4][0]=0.;
  range[4][1]=150.;   
  // m(jj)
  range[5][0]=200.;
  range[5][1]=800.;   
  // etaj1 - etaj2
  range[6][0]=-10.;
  range[6][1]=10.;
  // pT(jj)
  range[7][0]=0.;
  range[7][1]=300.;   
  // njets
  range[8][0]=0.;
  range[8][1]=10.;
  // nvtx
  range[9][0]=1.;
  range[9][1]=26.;

  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="PF met";
  xaxisLabel[1]="charged met";
  xaxisLabel[2]="PF W m_{T}";
  xaxisLabel[3]="charged W m_{T}";
  xaxisLabel[4]="lepton pT";
  xaxisLabel[5]="m(jj)";
  xaxisLabel[6]="etaj1 - etaj2";
  xaxisLabel[7]="pT (jj)";
  xaxisLabel[8]="n jets";
  xaxisLabel[9]="n vtx";

  TString binSize[NVARIABLES];
  for (int z=0;z<NVARIABLES;++z) {
    for (int j=0;j<NCUTS;++j) {
      sprintf(icut[j],"icut%d",j);
      // for (int i=0;i<NSPECIES;++i) {   // chiara
      for (int i=1;i<NSPECIES;++i) {
	histos[i][j][z]=new TH1F(variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),variables[z]+"_W_"+species[i]+"_"+TString(icut[j]),nbins[z],range[z][0],range[z][1]);
	char binsiz[10];
	sprintf(binsiz,"%2.0f",(range[z][1]-range[z][0])/nbins[z]);
	binSize[z]=TString(binsiz);
      }
    }
  }

  TString cut[NCUTS];
  cut[0]="(step[7] && njets>=2)*";

  TString intLumi="1000";
  TFile *_file[NSPECIES];
  TTree *T1[NSPECIES];
  
  TCanvas* c1= new TCanvas("test","test",800,800);
  
  //  for (int i=0;i<NSPECIES;++i) {    // chiara
  for (int i=1;i<NSPECIES;++i) {
    _file[i]=TFile::Open(files[i]);
    T1[i] = (TTree*)_file[i]->Get("t1new");
  }
  
  int nspeciesToRun=NSPECIES;
  for (int z=0;z<NVARIABLES;++z) {
    for (int j=0;j<NCUTS;++j) {
      // for (int i=0;i<nspeciesToRun;++i) {   // chiara
      for (int i=1;i<nspeciesToRun;++i) {
	fOut->cd();
	TString histoName=variables[z]+"_W_"+species[i]+"_"+TString(icut[j]);
	std::cout << "Producing " << histoName << std::endl;
	if (T1[i]==0) {
	  std::cout << "Tree not found" << std::endl;
	  return;
	}
	if (i!=0)
	  T1[i]->Project(histoName,variables[z],cut[j]+"baseW");
	  // T1[i]->Project(histoName,variables[z],cut[j]+"baseW*puW*effW");   // chiara
	else
	  T1[i]->Project(histoName,variables[z],cut[j]+"1.");
	std::cout << "Done " << histoName << std::endl;
      }
      
      THStack histo_MC(variables[z]+"_MC",variables[z]+"_MC");
      for (int i=1;i<nspeciesToRun;++i) {
	histos[i][j][z]->SetFillColor(colors[i]);
	histos[i][j][z]->SetLineColor(lineColors[i]);
	histo_MC.Add(histos[i][j][z]);
      }
	  
      // float maximum=TMath::Max(histo_MC.GetMaximum(),histos[0][j][z]->GetMaximum());  // chiara
      float maximum=histo_MC.GetMaximum();
      histo_MC.SetMinimum(0.01);
      histo_MC.SetMaximum(maximum*1.6);
      histo_MC.Draw("");
      histo_MC.GetXaxis()->SetTitle(xaxisLabel[z]+" ["+units[z]+"]");
      histo_MC.GetYaxis()->SetTitle("Events/"+binSize[z]+" "+units[z]);
      
      // histos[0][j][z]->SetMarkerStyle(20);   // chiara
      // histos[0][j][z]->SetMarkerSize(1.3);
      // histos[0][j][z]->Draw("EP1SAME");

      TPaveText pt1(0.6,0.83,0.8,0.9,"NDC");
      pt1.SetTextSize(0.028);
      pt1.SetTextAlign(12);
      pt1.SetFillColor(0);
      pt1.SetBorderSize(0);
      pt1.AddText("CMS Preliminary 2011");
      pt1.AddText("");
      pt1.AddText("");
      pt1.AddText("");
      pt1.AddText("#sqrt{s}=7 TeV L_{int}="+intLumi+" pb^{-1}");
      pt1.Draw();
      
      c1->Update();
      TLegendEntry *legge;
      TLegend *leg;
      leg = new TLegend(0.6,0.65,0.93,0.8,cut[j]);
      leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.025);
      leg->SetFillColor(0);
      // for (int i=0;i<nspeciesToRun;++i) {   // chiara
      for (int i=1;i<nspeciesToRun;++i) {
	if (i == 0)
	  legge = leg->AddEntry(histos[legendOrder[i]][j][z],"Data", "lpe");
	else
	  legge = leg->AddEntry(histos[legendOrder[i]][j][z],species[legendOrder[i]],"f");
      }
      leg->Draw();
      c1->SetLogy(0);
      c1->Update();
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+"_"+suffix+".png");
      c1->SaveAs(plotsDir+variables[z]+"DataMC_"+TString(icut[j])+"_"+suffix+".root");
    }
  }
  
  fOut->Write();
  fOut->Close();
  
}
