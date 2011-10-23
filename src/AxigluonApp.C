
// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "AxigluonAnalysisTools/include/Application.hh"
#include "CommonTools/include/TriggerMask.hh"
#if Application == 1
#include "AxigluonAnalysisTools/src/AxigluonSelection.cc"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[150];
  char outputFileName[150];
  char dataset[150];
  if ( argc < 2 ){
    std::cout << "missing argument: insert at least inputFile with list of root files" << std::endl; 
    std::cout << "AxiApp inputFile [outputFile] [1=MC,0=data] [dataset]" << std::endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  if (argc < 3 ) strcpy(outputFileName,argv[1]);
  else strcpy(outputFileName,argv[2]);
  int isMC=1;
  if(argc==5) {
    isMC=atoi(argv[3]);
    strcpy(dataset,argv[4]);
  }
  
  // -------------------------
  // Loading the file
  TChain *theChain = new TChain("ntp1");
  char Buffer[500];
  char MyRootFile[2000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  char tmpFileName[256];
  vector<string> filesToRemove;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
        sscanf(Buffer,"%s",MyRootFile);
        theChain->Add(TString(MyRootFile));
        std::cout << "chaining " << MyRootFile << std::endl;
      }
  }
  inputFile->close();
  delete inputFile;

#if Application == 1

  AxigluonSelection axi(theChain);
  axi.SetDatasetName(outputFileName);

  std::vector<std::string> maskEle, maskMu;
  std::vector<std::string> maskNotEle, maskNotMu;
  
  if(isMC) {
    maskEle.push_back("1-1:HLT_Ele17_SW_TighterEleIdIsol_L1R_v3");
    maskMu.push_back("1-1:HLT_DoubleMu5_v1");
  } else {
    TString DatasetName(dataset);
    if(DatasetName.Contains("DoubleElectron")) {
      maskEle.push_back("1-170052:HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v");
    } else if(DatasetName.Contains("DoubleMu")) {
      maskMu.push_back("1-164237:HLT_DoubleMu7_v");
    }
  }

  axi.setRequiredTriggers(maskEle,0);
  axi.setRequiredTriggers(maskMu,1);
  axi.setNotRequiredTriggers(maskNotEle,0);
  axi.setNotRequiredTriggers(maskNotMu,1);

  axi.Loop();
  axi.displayEfficiencies(outputFileName);

#endif

  return 0;

}
