#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "AxigluonAnalysisTools/include/AxigluonSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"

#include <iostream>
#include <string>
#include <algorithm>

#include <TTree.h>

using namespace bits;

AxigluonSelection::AxigluonSelection(TTree *tree) 
  : Axigluon(tree) {
  
  // choose the Axigluon Mass
  std::string axiConfigDir;
  std::string axiConfigDirMass;
  std::ifstream setfile("config/axi/axiMass.txt");
  axiConfigDir="config/axi/";
  if(!setfile.good()) {
    std::cout << "Cannot read the axigluon mass file to choose the selection: config/axi/axiMass.txt" << std::endl
	      << "using default (500 GeV)" << std::endl;
    axiConfigDirMass="config/axi/h500/";
  }
  else {
    std::string var, massVal; 
    bool found=false;
    while (1) {
      setfile >> var >> massVal;
      _massVal = atoi(massVal.c_str());
      if(!setfile.good()) break;
      if(var.compare("axiMass")==0) { 
	found=true;
	axiConfigDirMass="config/axi/h" + massVal + "/";
	std::cout << "Reading configuration for Axi mass = " << massVal << " GeV/c^2" << std::endl;
	break;
      }
    }
  }
  
  // selection efficiencies
  std::string fileCuts     = axiConfigDirMass + "lnuCuts.txt";  
  std::string fileSwitches = axiConfigDir + "lnuSwitches.txt";  
  cout << "debug: axiConfigDir = " << axiConfigDir << ", fileSwitches = " << fileSwitches.c_str() << endl;
  cout << "debug: axiConfigDirMass = " << axiConfigDirMass << ", fileCuts = " << fileCuts.c_str() << endl;
  CutBasedAxiSelectionEle.Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER EE"); 
  CutBasedAxiSelectionMu.Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER MM"); 
  _selectionEle = CutBasedAxiSelectionEle.GetSelection();  
  _selectionMu  = CutBasedAxiSelectionMu.GetSelection();

  //  extra selection efficiencies  - to be put here not to pass the full list of leptons to the preselection class
  _selectionEle->addCut("etaElectronAcc");    
  _selectionEle->addCut("electronBigCrack");
  _selectionEle->addCut("ptElectronAcc");
  _selectionEle->addCut("etaMuonAcc");
  _selectionEle->addCut("ptMuonAcc");
  _selectionEle->addCut("muGlobalIso");
  _selectionEle->addCut("electronIP");
  _selectionEle->addCut("electronDz");
  _selectionEle->addCut("muonIP");
  _selectionEle->addCut("muonDz");
  _selectionEle->addCut("softMuPt");  
  _selectionEle->addCut("etJetAcc");  
  _selectionEle->addCut("etaJetAcc");  
  _selectionEle->addCut("jetConeWidth");  
  _selectionEle->addSwitch("isData"); 
  _selectionEle->addSwitch("goodRunLS");
  _selectionEle->addSwitch("ecalDrivenOnly");
  _selectionEle->addSwitch("MCtruth");
  _selectionEle->addSwitch("trigger");
  _selectionEle->addSwitch("leptonId");
  _selectionEle->addSwitch("leptonIso");
  _selectionEle->addSwitch("leptonD0");
  _selectionEle->addSwitch("convRej");
  _selectionEle->addStringParameter("electronIDType");   

  // configuring electron id
  TString selectionString(_selectionEle->getStringParameter("electronIDType"));
  cout << "=== CONFIGURING " << selectionString << " ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/axi/electronId/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/axi/electronId/");

  // configuring the quark-gluon likelihood
  qglikeli = new QGLikelihoodCalculator();
  
  // Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _selectionEle->getSwitch("goodRunLS") << " isData is " <<  _selectionEle->getSwitch("isData") << std::endl;

  // To read good run list!
  if (_selectionEle->getSwitch("goodRunLS") && _selectionEle->getSwitch("isData")) {
    std::string goodRunJsonFile = "config/json/goodCollisions2011_axi.json";
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics
  m_p3PFMet = new TVector3(0.,0.,0.);   
  for(int theChannel=0; theChannel<2; theChannel++) { 
    m_p4Lepton[theChannel]  = new TLorentzVector(0.,0.,0.,0.);
    m_p3Lepton[theChannel]  = new TVector3(0.,0.,0.);
  }
}

AxigluonSelection::~AxigluonSelection(){

  delete m_p3PFMet;    
  for(int theChannel=0; theChannel<2; theChannel++) { 
    delete m_p4Lepton[theChannel];
    delete m_p3Lepton[theChannel];
  }
  delete _selectionEle;
  delete _selectionMu;
  
  delete qglikeli;

  myOutTreeEle -> save();
  myOutTreeMu  -> save();
}

void AxigluonSelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;
  
  // kinematics reduced tree
  std::string reducedTreeNameEle = _datasetName+"-datasetEle.root";
  std::string reducedTreeNameMu  = _datasetName+"-datasetMu.root";
  myOutTreeEle = new RedAxiTree(reducedTreeNameEle.c_str());
  myOutTreeMu  = new RedAxiTree(reducedTreeNameMu.c_str());
  myOutTreeEle->addElectronInfos();
  myOutTreeEle->addKinematics();
  myOutTreeMu ->addKinematics();
  myOutTreeEle->addSteps();
  myOutTreeMu ->addSteps();

  // running over all events
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // PU reweighting
  PUWeight* fPUWeight = new PUWeight();  // Isidro's function

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    if(nEle>100)  { cout << "in this event there are " << nEle  << " electrons: skip it" << endl; continue; }
    if(nMuon>100) { cout << "in this event there are " << nMuon << " muons: skip it" << endl; continue; } 
    
    resetKinematicsStart();

    // event weight
    float weight = 1;
    
    // weight for the PU observed in 2011 data = chiara, per il momento off
    // if ( !_selectionEle->getSwitch("isData") ) weight *= fPUWeight->GetWeight(nPU[1]);    // Isidro's function
 
    // Good Run selection
    if (_selectionEle->getSwitch("isData") && _selectionEle->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_selectionEle->getSwitch("isData") && _selectionEle->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    
    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask(runNumber);
    bool passedHLT[2];
    passedHLT[ee] = hasPassedHLT(ee);
    passedHLT[mm] = hasPassedHLT(mm);

    // -------------------------------------------------------------
    // vertex selection - we only consider the first vertex of the list ( = highest sumPT^2)
    bool isGoodVertex = true;
    if (nPV<1) isGoodVertex = false;
    float rhoVtx = sqrt(PVxPV[0]*PVxPV[0] + PVyPV[0]*PVyPV[0]);
    if ( isFakePV[0] )       isGoodVertex = false;
    if ( ndofPV[0]<4 )       isGoodVertex = false;
    if ( fabs(PVzPV[0])>24.) isGoodVertex = false;
    if ( rhoVtx>2 )          isGoodVertex = false; 
    
    // -------------------------------------------------------------
    
    // MC truth infos
    bool promptEle = false;
    bool promptMu  = false;
    if( _selectionEle->getSwitch("MCtruth") ) {
      for (int iMC=0; iMC<30; iMC++) {
	int thisId = idMc[iMC];
	int mother = mothMc[iMC];
	int mothId = idMc[mother];
	if (abs(thisId)==11 && abs(mothId)==24) promptEle = true;
	if (abs(thisId)==13 && abs(mothId)==24) promptMu  = true;
	if (promptEle || promptMu) continue; 
      }
    }

    // get the best electron and best muon ==> tu be used to select ALL the possible channels at the beginning only
    int thePreElectron = getBestElectronPair_acceptance();
    int thePreMuon     = getBestMuonPair_acceptance();

    // reconstructed channel
    m_channel[ee] = false;     
    m_channel[mm] = false;
    if ( (thePreElectron>-1) && isGoodVertex ) m_channel[ee] = true;    
    if ( (thePreMuon>-1) && isGoodVertex )     m_channel[mm] = true;  
    
    if (_verbose) {
      std::cout << "nEle = "      << nEle << "\tnMuon = "  << nMuon << std::endl;
      std::cout << "indices: "    << thePreElectron << " " << thePreMuon << std::endl;
      std::cout << "chargeEle = " << chargeEle[thePreElectron] 
		<< "\t chargeMuon = " << chargeMuon[thePreMuon] << std::endl;
      std::cout << "channel: ee = " << m_channel[ee] << "\tmm = " << m_channel[mm] << std::endl;
    }


    // -------------------------------------------------------------
    // W->enu candidates: preparing vectors of candidates and selecting the highest pT ele after each step

    // eleID, for electrons in acceptance
    int theBestIdEle = getBestElectronPair_id(_acceptEleAll);   

    // isolation, for identified electrons
    int theBestIsolEle = getBestElectronPair_isol(_idEleAll); 

    // conversion rejection, for isolated electrons
    int theBestConvEle = getBestElectronPair_conv(_isolEleAll);     

    // transverse impact parameter, for electrons passing conversion rejection
    int theBestIpEle = getBestElectronPair_ip(_convEleAll);     

    // the highest pT electron at this point is the one I use for my analysis since it passed the full lepton selection
    theElectron = theBestIpEle;
    
    // to be used in the following
    int theIdElectron(theBestIdEle);
    int theIsolElectron(theBestIsolEle);
    int theConvElectron(theBestConvEle);
    int theIpElectron(theBestIpEle);
    

    // -------------------------------------------------------------
    // W->munu candidates: preparing vectors of candidates and selecting the highest pT mu after each step

    // muID, for muons in acceptance
    int theBestIdMuon = getBestMuonPair_id(_acceptMuonsAll); 

    // isolation, for identified muons
    int theBestIsolMuon = getBestMuonPair_isol(_idMuonsAll); 

    // transverse impact parameter, for isolated muons
    int theBestIpMuon = getBestMuonPair_ip(_isolMuonsAll);     
    
    // the highest pT muon at this point is the one I use for my analysis since it passed the full lepton selection
    theMuon = theBestIpMuon;

    // to be used in the following
    int theIdMuonMinus(theBestIdMuon);
    int theIsolMuonMinus(theBestIsolMuon);
    int theIpMuonMinus(theBestIpMuon);


    // -------------------------------------------------------------
    // set of kinematics: : now I've the final lepton 
    resetKinematics();
    
    // MET is an event variable. Independent on the channel
    m_p3PFMet -> SetXYZ(pxPFMet[0],pyPFMet[0],pzPFMet[0]);  
    modPFMet = m_p3PFMet->Pt();   
    
    // setting all the channel dependent variables
    setKinematicsEle(theElectron);
    setKinematicsMu(theMuon);

    // -------------------------------------------------------------    
    int njets[2], nsoftmu[2], nextraleptons[2];
    float dijetInvMass[2], dijetPt[2], dijetDeta[2];
    float QGLikeLead[2], QGLikeSublead[2], QGLikeProd[2];

    // initialize the btags for the leading and subleading jets to unphysical value
    for(int ichan=0; ichan<2; ichan++) {
      leadJetBtag[ichan]    = -2000.;   
      subleadJetBtag[ichan] = -2000.;   
    }

    // variables useful for the selection
    for(int ichan=0; ichan<2; ichan++) {

      int lead    = theLeadingJet[ichan];
      int sublead = theSecondJet[ichan];

      // jet counter 
      njets[ichan] = numJets(eleCands[ichan],muCands[ichan],ichan);

      // soft muon counter 
      nsoftmu[ichan] = numSoftMuons(muCands[ichan]);
      
      // extra lepton counter
      nextraleptons[ichan] = numExtraLeptons(eleCands[ichan],muCands[ichan]);

      // di-jet invariant mass
      dijetInvMass[ichan] = getDiJetInvMass(lead, sublead);

      // di-jet pT
      dijetPt[ichan] = getDiJetPt(lead, sublead);
      
      // deltaEta between the two jets
      dijetDeta[ichan] = getDiJetDeta(lead, sublead);

      // quark-gluon likelihood 
      float rhoPF = rhoFastjet; 
      float leadPt    = GetPt(pxAK5PFPUcorrJet[lead],pyAK5PFPUcorrJet[lead]);
      float subleadPt = GetPt(pxAK5PFPUcorrJet[sublead],pyAK5PFPUcorrJet[sublead]);
      int leadNcharged    = chargedHadronMultiplicityAK5PFPUcorrJet[lead];
      int subleadNcharged = chargedHadronMultiplicityAK5PFPUcorrJet[sublead];
      int leadNneutral    = neutralHadronMultiplicityAK5PFPUcorrJet[lead] + photonMultiplicityAK5PFPUcorrJet[lead];
      int subleadNneutral = neutralHadronMultiplicityAK5PFPUcorrJet[sublead] + photonMultiplicityAK5PFPUcorrJet[sublead];
      float leadPtD    = ptDAK5PFPUcorrJet[lead];
      float subleadPtD = ptDAK5PFPUcorrJet[sublead];
      QGLikeLead[ichan]    = qglikeli->computeQGLikelihoodPU( leadPt,    rhoPF, leadNcharged,    leadNneutral,    leadPtD );
      QGLikeSublead[ichan] = qglikeli->computeQGLikelihoodPU( subleadPt, rhoPF, subleadNcharged, subleadNneutral, subleadPtD );
      QGLikeProd[ichan]    = QGLikeLead[ichan]*QGLikeSublead[ichan];
    }

    // ---------------------------------------
    // filling counters for the different final states

    // W->enu
    CutBasedAxiSelectionEle.SetWeight(weight);               
    CutBasedAxiSelectionEle.SetMcTruth(promptEle);               
    CutBasedAxiSelectionEle.SetHLT(passedHLT[ee]);               
    CutBasedAxiSelectionEle.SetIsChannel(m_channel[ee]);     
    CutBasedAxiSelectionEle.SetElectronId(theIdElectron);                 
    CutBasedAxiSelectionEle.SetElectronIsolation(theIsolElectron);        
    CutBasedAxiSelectionEle.SetElectronConvRejection(theConvElectron);    
    CutBasedAxiSelectionEle.SetElectronIp(theIpElectron);                 
    CutBasedAxiSelectionEle.SetLeptPt(leptonPt[ee]);     
    CutBasedAxiSelectionEle.SetNJets(njets[ee]);
    CutBasedAxiSelectionEle.SetBTagLeadJet(leadJetBtag[ee]);
    CutBasedAxiSelectionEle.SetBTagSubleadJet(subleadJetBtag[ee]);
    CutBasedAxiSelectionEle.SetNSoftMuons(nsoftmu[ee]);
    CutBasedAxiSelectionEle.SetNExtraLeptons(nextraleptons[ee]);
    CutBasedAxiSelectionEle.SetDijetMass(dijetInvMass[ee]);
    CutBasedAxiSelectionEle.SetDijetPt(dijetPt[ee]);
    CutBasedAxiSelectionEle.SetDijetDeltaEta(dijetDeta[ee]);
    CutBasedAxiSelectionEle.SetMet(modPFMet);         // chiara: qui dovrai avere gia' deciso quale met usare	    
    CutBasedAxiSelectionEle.SetWMt(WmT_PFMet[ee]);    // chiara: qui dovrai avere gia' deciso quale met usare	    

    // after each step
    bool isSelected    = CutBasedAxiSelectionEle.output();
    bool outputStep0   = CutBasedAxiSelectionEle.outputStep0();
    bool outputStep1   = CutBasedAxiSelectionEle.outputStep1();
    bool outputStep2   = CutBasedAxiSelectionEle.outputStep2();
    bool outputStep3   = CutBasedAxiSelectionEle.outputStep3();
    bool outputStep4   = CutBasedAxiSelectionEle.outputStep4();
    bool outputStep5   = CutBasedAxiSelectionEle.outputStep5();
    bool outputStep6   = CutBasedAxiSelectionEle.outputStep6();
    bool outputStep7   = CutBasedAxiSelectionEle.outputStep7();
    bool outputStep8   = CutBasedAxiSelectionEle.outputStep8();
    bool outputStep9   = CutBasedAxiSelectionEle.outputStep9();
    bool outputStep10  = CutBasedAxiSelectionEle.outputStep10();
    bool outputStep11  = CutBasedAxiSelectionEle.outputStep11();
    bool outputStep12  = CutBasedAxiSelectionEle.outputStep12();
    bool outputStep13  = CutBasedAxiSelectionEle.outputStep13();
    bool outputStep14  = CutBasedAxiSelectionEle.outputStep14();
    bool outputStep15  = CutBasedAxiSelectionEle.outputStep15();
    bool outputStep16  = CutBasedAxiSelectionEle.outputStep16();

    myOutTreeEle -> fillRunInfos(runNumber, lumiBlock, eventNumber, weight, passedHLT[ee]);

    myOutTreeEle -> fillAll(nPV, leptonPt[ee], njets[ee], leadJetBtag[ee], subleadJetBtag[ee], nsoftmu[ee], nextraleptons[ee], 
			    dijetInvMass[ee], dijetPt[ee], dijetDeta[ee], modPFMet, modChPFMet[ee], WmT_PFMet[ee], WmT_ChPFMet[ee],
			    QGLikeLead[ee], QGLikeSublead[ee], QGLikeProd[ee]);
    
    myOutTreeEle -> fillSteps(outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16);

    setEleIdVariables(theElectron);
    myOutTreeEle -> fillElectrons(myPt, myEta, myPhi, myDeta, myDphi, myHoe, mySee, myTrackerIso, myHcalIso, myEcalIso, myCombinedIso, myCharge, myLh);

    int theLJ  = theLeadingJet[ee];
    int theSJ  = theSecondJet[ee];
    float pxLJ = pxAK5PFPUcorrJet[theLJ];
    float pyLJ = pyAK5PFPUcorrJet[theLJ];
    float pzLJ = pzAK5PFPUcorrJet[theLJ];
    float pxSJ = pxAK5PFPUcorrJet[theSJ];
    float pySJ = pyAK5PFPUcorrJet[theSJ];
    float pzSJ = pzAK5PFPUcorrJet[theSJ];

    myOutTreeEle -> fillKinematics(m_p3PFMet[ee].Px(), m_p3PFMet[ee].Py(), m_p3PFMet[ee].Pz(), 
				   m_p3ChPFMet[ee].Px(), m_p3ChPFMet[ee].Py(), m_p3ChPFMet[ee].Pz(),
				   pxLJ, pyLJ, pzLJ, pxSJ, pySJ, pzSJ,
				   m_p4Lepton[ee]->Px(), m_p4Lepton[ee]->Py(), m_p4Lepton[ee]->Pz()) ;

    // dumping final tree, only if there is 1 lepton in the acceptance
    if(outputStep1) myOutTreeEle -> store();


    
    // ---------------------------------------
    // W -> munu
    CutBasedAxiSelectionMu.SetWeight(weight);               
    CutBasedAxiSelectionMu.SetMcTruth(promptMu);               
    CutBasedAxiSelectionMu.SetHLT(passedHLT[mm]);               
    CutBasedAxiSelectionMu.SetIsChannel(m_channel[mm]);     
    CutBasedAxiSelectionMu.SetElectronId(theIdMuonMinus);
    CutBasedAxiSelectionMu.SetElectronIsolation(theIsolMuonMinus);
    CutBasedAxiSelectionMu.SetElectronConvRejection(theIsolMuonMinus);
    CutBasedAxiSelectionMu.SetElectronIp(theIpMuonMinus);
    CutBasedAxiSelectionMu.SetLeptPt(leptonPt[mm]);     
    CutBasedAxiSelectionMu.SetNJets(njets[mm]);
    CutBasedAxiSelectionMu.SetBTagLeadJet(leadJetBtag[mm]);
    CutBasedAxiSelectionMu.SetBTagSubleadJet(subleadJetBtag[mm]);
    CutBasedAxiSelectionMu.SetNSoftMuons(nsoftmu[mm]);
    CutBasedAxiSelectionMu.SetNExtraLeptons(nextraleptons[mm]);
    CutBasedAxiSelectionMu.SetDijetMass(dijetInvMass[mm]);
    CutBasedAxiSelectionMu.SetDijetPt(dijetPt[mm]);
    CutBasedAxiSelectionMu.SetDijetDeltaEta(dijetDeta[mm]);
    CutBasedAxiSelectionMu.SetMet(modPFMet);         // chiara: qui dovrai avere gia' deciso quale met usare	    
    CutBasedAxiSelectionMu.SetWMt(WmT_PFMet[mm]);    // chiara: qui dovrai avere gia' deciso quale met usare	    

    // after each step
    isSelected    = CutBasedAxiSelectionMu.output();
    outputStep0   = CutBasedAxiSelectionMu.outputStep0();
    outputStep1   = CutBasedAxiSelectionMu.outputStep1();
    outputStep2   = CutBasedAxiSelectionMu.outputStep2();
    outputStep3   = CutBasedAxiSelectionMu.outputStep3();
    outputStep4   = CutBasedAxiSelectionMu.outputStep4();
    outputStep5   = CutBasedAxiSelectionMu.outputStep5();
    outputStep6   = CutBasedAxiSelectionMu.outputStep6();
    outputStep7   = CutBasedAxiSelectionMu.outputStep7();
    outputStep8   = CutBasedAxiSelectionMu.outputStep8();
    outputStep9   = CutBasedAxiSelectionMu.outputStep9();
    outputStep10  = CutBasedAxiSelectionMu.outputStep10();
    outputStep11  = CutBasedAxiSelectionMu.outputStep11();
    outputStep12  = CutBasedAxiSelectionMu.outputStep12();
    outputStep13  = CutBasedAxiSelectionMu.outputStep13();
    outputStep14  = CutBasedAxiSelectionMu.outputStep14();
    outputStep15  = CutBasedAxiSelectionMu.outputStep15();
    outputStep16  = CutBasedAxiSelectionMu.outputStep16();

    myOutTreeMu -> fillRunInfos(runNumber, lumiBlock, eventNumber, weight, passedHLT[mm]);

    myOutTreeMu -> fillAll(nPV, leptonPt[mm], njets[mm], leadJetBtag[mm], subleadJetBtag[mm], nsoftmu[mm], nextraleptons[mm], 
			   dijetInvMass[mm], dijetPt[mm], dijetDeta[mm], modPFMet, modChPFMet[mm], WmT_PFMet[mm], WmT_ChPFMet[mm],
			   QGLikeLead[mm], QGLikeSublead[mm], QGLikeProd[mm]);
    
    myOutTreeMu -> fillSteps(outputStep0, outputStep1, outputStep2, outputStep3, outputStep4, outputStep5, outputStep6, outputStep7, outputStep8, outputStep9, outputStep10, outputStep11, outputStep12, outputStep13, outputStep14, outputStep15, outputStep16);
    
    int theLJm  = theLeadingJet[mm];
    int theSJm  = theSecondJet[mm];
    float pxLJm = pxAK5PFPUcorrJet[theLJm];
    float pyLJm = pyAK5PFPUcorrJet[theLJm];
    float pzLJm = pzAK5PFPUcorrJet[theLJm];
    float pxSJm = pxAK5PFPUcorrJet[theSJm];
    float pySJm = pyAK5PFPUcorrJet[theSJm];
    float pzSJm = pzAK5PFPUcorrJet[theSJm];

    myOutTreeMu -> fillKinematics(m_p3PFMet[mm].Px(), m_p3PFMet[mm].Py(), m_p3PFMet[mm].Pz(), 
				  m_p3ChPFMet[mm].Px(), m_p3ChPFMet[mm].Py(), m_p3ChPFMet[mm].Pz(),
				  pxLJm, pyLJm, pzLJm, pxSJm, pySJm, pzSJm,
				  m_p4Lepton[mm]->Px(), m_p4Lepton[mm]->Py(), m_p4Lepton[mm]->Pz()) ;

    // dumping final tree, only if there is 1 lepton in the acceptance
    if(outputStep1) myOutTreeMu -> store();
  }
}


void AxigluonSelection::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full W->enu selections: " << std::endl;
  CutBasedAxiSelectionEle.displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full W->munu selections: " << std::endl;
  CutBasedAxiSelectionMu.displayEfficiencies(datasetName);

  // simple cuts based or like based ele id
  std::cout << "electron ID: " << std::endl;
  EgammaCutBasedID.displayEfficiencies();
}

int AxigluonSelection::getBestElectronPair_acceptance() {   
  
  int theLep = -1;
  float maxPtLep = -1000.;
  
  _acceptEleAll.clear();

  for(int i=0;i<nEle;i++) {

    TVector3 pLepton(pxEle[i],pyEle[i],pzEle[i]);
    float thisPt=pLepton.Pt();

    // only ecal driven
    Utils anaUtils;
    bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[i], isEcalDriven);
    if(_selectionEle->getSwitch("ecalDrivenOnly") && !ecaldriven ) continue;    

    if(_selectionEle->getSwitch("etaElectronAcc") && !_selectionEle->passCut("etaElectronAcc",etaEle[i]) ) continue;

    if(_selectionEle->getSwitch("electronBigCrack") && _selectionEle->passCut("electronBigCrack",fabs(etaEle[i])) ) continue;

    if(_selectionEle->getSwitch("ptElectronAcc") && !_selectionEle->passCut("ptElectronAcc",thisPt) ) continue;
    
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = i; }

    _acceptEleAll.push_back(i);   
  }
  _acceptEleAll = sortElectronsByPt(_acceptEleAll);

  return theLep;
}

int AxigluonSelection::getBestElectronPair_id( std::vector<int> acceptEle ) {  

  int theLep = -1;
  float maxPtLep = -1000.;

  _idEleAll.clear();

  for (int iEle=0; iEle<acceptEle.size(); iEle++) {
    int thisEle = acceptEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    float thisPt = GetPt(pxEle[thisEle],pyEle[thisEle]);
    isEleIDAndDenom(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    if (!theElectronID) continue;
    
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisEle; }

    _idEleAll.push_back(thisEle);  
  }
  _idEleAll = sortElectronsByPt(_idEleAll);

  return theLep;
}

int AxigluonSelection::getBestElectronPair_isol( std::vector<int> idEle ) {  

  int theLep = -1;
  float maxPtLep = -1000.;

  _isolEleAll.clear();

  for (int iEle=0; iEle<idEle.size(); iEle++) {
    int thisEle = idEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    isEleIDAndDenom(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);

    if (!theElectronIsol) continue;
    
    float thisPt = GetPt(pxEle[thisEle],pyEle[thisEle]);
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisEle; }

    _isolEleAll.push_back(thisEle);  
  }
  _isolEleAll = sortElectronsByPt(_isolEleAll);

  return theLep;
}

int AxigluonSelection::getBestElectronPair_conv( std::vector<int> isolEle ) { 

  int theLep = -1;
  float maxPtLep = -1000.;
  
  _convEleAll.clear();

  for (int iEle=0; iEle<isolEle.size(); iEle++) {
    int thisEle = isolEle[iEle];
    
    bool theElectronID, theElectronIsol, theElectronConvRej;
    theElectronID = theElectronIsol = theElectronConvRej = true;
    
    isEleIDAndDenom(thisEle,&theElectronID,&theElectronIsol,&theElectronConvRej,&EgammaCutBasedID);
    
    if (!theElectronConvRej) continue;

    float thisPt = GetPt(pxEle[thisEle],pyEle[thisEle]);
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisEle; }

    _convEleAll.push_back(thisEle);      
  }
  _convEleAll = sortElectronsByPt(_convEleAll);

  return theLep;
}

int AxigluonSelection::getBestElectronPair_ip( std::vector<int> convEle ) {  // chiara ok, ma i tagli da definire con francesco

  int theLep = -1;
  float maxPtLep = -1000.;

  _ipEleAll.clear();

  for (int iEle=0; iEle<convEle.size(); iEle++) {
    int thisEle = convEle[iEle];

    int gsfTrack = gsfTrackIndexEle[thisEle]; 
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[gsfTrack];   
    if (_selectionEle->getSwitch("electronIP") && (!_selectionEle->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionEle->getSwitch("electronDz") && (!_selectionEle->passCut("electronDz",dzEle)) )  continue;

    float thisPt = GetPt(pxEle[thisEle],pyEle[thisEle]);
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisEle; }

    _ipEleAll.push_back(thisEle);  
  }
  _ipEleAll = sortElectronsByPt(_ipEleAll);

  return theLep;
}

int AxigluonSelection::getBestMuonPair_acceptance() {  
  
  int theLep = -1;
  float maxPtLep = -1000.;
  
  _acceptMuonsAll.clear();

  for(int i=0;i<nMuon;i++) {

    if(_selectionEle->getSwitch("etaMuonAcc") && !_selectionEle->passCut("etaMuonAcc",etaMuon[i]) ) continue;

    float thisPt = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEle->getSwitch("ptMuonAcc") && !_selectionEle->passCut("ptMuonAcc",thisPt) ) continue;
    
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = i; }
    
    _acceptMuonsAll.push_back(i);  
  }
  _acceptMuonsAll = sortMuonsByPt(_acceptMuonsAll);

  return theLep;
}

int AxigluonSelection::getBestMuonPair_id( std::vector<int> acceptMu ) {   
   
  int theLep = -1;
  float maxPtLep = -1000.;
  
  _idMuonsAll.clear();

  for(int iMu=0; iMu<acceptMu.size(); iMu++) {
    
    int thisMu = acceptMu[iMu];

    bool theMuonID = isAxigluonMuonID(thisMu);   // "tight" muonID by muPOG
    if (!theMuonID) continue;

    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);
    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisMu; }
    
    _idMuonsAll.push_back(thisMu);   
  }
  _idMuonsAll = sortMuonsByPt(_idMuonsAll);
  
  return theLep;
}

int AxigluonSelection::getBestMuonPair_isol( std::vector<int> idMu ) { 
  
  int theLep = -1;
  float maxPtLep = -1000.;
   
  _isolMuonsAll.clear();

  for(int iMu=0; iMu<idMu.size(); iMu++) {

    int thisMu   = idMu[iMu];
    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);
    
    float muonTrackerForGlobal = sumPt03Muon[thisMu];
    float muonEcalForGlobal    = emEt03Muon[thisMu];
    float muonHcalForGlobal    = hadEt03Muon[thisMu]; 
    float theMuonGlobalSum     = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal - rhoFastjet*TMath::Pi()*0.3*0.3;
    float theRelMuonIso        = theMuonGlobalSum/thisPt; 
    if(_selectionEle->getSwitch("muGlobalIso") && !_selectionEle->passCut("muGlobalIso",theRelMuonIso)) continue;  

    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisMu; }

    _isolMuonsAll.push_back(thisMu);   
  }
  _isolMuonsAll = sortMuonsByPt(_isolMuonsAll);

  return theLep;
}

int AxigluonSelection::getBestMuonPair_ip( std::vector<int> isoMu ) {  

  int theLep = -1;
  float maxPtLep = -1000.;

  _ipMuonsAll.clear();  

  for(int iMu=0; iMu<isoMu.size(); iMu++) {

    int thisMu = isoMu[iMu];

    float thisPt = GetPt(pxMuon[thisMu],pyMuon[thisMu]);    
    
    int ctfMuon   = trackIndexMuon[thisMu]; 
    float dxyMuon = transvImpactParTrack[ctfMuon];
    float dzMuon  = PVzPV[0] - trackVzTrack[ctfMuon];   
    if (_selectionEle->getSwitch("muonIP") && (!_selectionEle->passCut("muonIP",dxyMuon)) ) continue;   
    if (_selectionEle->getSwitch("muonDz") && (!_selectionEle->passCut("muonDz",dzMuon)) )  continue;   

    if (thisPt> maxPtLep){ maxPtLep = thisPt; theLep = thisMu; }

    _ipMuonsAll.push_back(thisMu);   
  }
  _ipMuonsAll = sortMuonsByPt(_ipMuonsAll);

  return theLep;
}

void AxigluonSelection::setKinematicsEle(int myEle) {

  if (myEle > -1) {
    leptonPt[ee] = GetPt(pxEle[myEle],pyEle[myEle]);
    eleCands[ee].push_back(myEle);
    m_p3Lepton[ee] -> SetXYZ(pxEle[myEle], pyEle[myEle], pzEle[myEle]);
    m_p4Lepton[ee] -> SetXYZT(pxEle[myEle], pyEle[myEle], pzEle[myEle], energyEle[myEle]);
    m_p3ChPFMet[ee] = oneLepPFChargedMet(m_p4Lepton[ee]->Vect());  
    modChPFMet[ee]  = (m_p3ChPFMet[ee]).Pt();     
    WmT_ChPFMet[ee] = sqrt(2 * m_p4Lepton[ee]->Pt() * modChPFMet[ee] * (1-cos(m_p3Lepton[ee]->Angle(m_p3ChPFMet[ee]))));
    WmT_PFMet[ee]   = sqrt(2 * m_p4Lepton[ee]->Pt() * modPFMet * (1-cos(m_p3Lepton[ee]->Angle(*m_p3PFMet))) );
  }
}

void AxigluonSelection::setKinematicsMu(int myMuon) {
  
  if (myMuon > -1) {
    leptonPt[mm] = GetPt(pxMuon[myMuon],pyMuon[myMuon]);
    muCands[mm].push_back(myMuon);
    m_p3Lepton[mm] -> SetXYZ(pxMuon[myMuon], pyMuon[myMuon], pzMuon[myMuon]);
    m_p4Lepton[mm] -> SetXYZT(pxMuon[myMuon], pyMuon[myMuon], pzMuon[myMuon], energyMuon[myMuon]);
    m_p3ChPFMet[mm] = oneLepPFChargedMet(m_p4Lepton[mm]->Vect());
    modChPFMet[mm]  = (m_p3ChPFMet[mm]).Pt();
    WmT_ChPFMet[mm] = sqrt(2 * m_p4Lepton[mm]->Pt() * modChPFMet[mm] * (1-cos(m_p3Lepton[mm]->Angle(m_p3ChPFMet[mm]))));
    WmT_PFMet[mm]   = sqrt(2 * m_p4Lepton[mm]->Pt() * modPFMet * (1-cos(m_p3Lepton[mm]->Angle(*m_p3PFMet))) );
  }
}

void AxigluonSelection::resetKinematicsStart() {

  theElectron    = -1;
  theMuon        = -1;
  thePreElectron = -1;
  thePreMuon     = -1;
}

void AxigluonSelection::resetKinematics() {

  m_p3PFMet -> SetXYZ(0,0,0);
  modPFMet = 0.;

  for(int theChannel=0; theChannel<2; theChannel++) {
    leptonPt[theChannel] = 0.;
    eleCands[theChannel].clear();
    muCands[theChannel] .clear();
    m_p3Lepton[theChannel] -> SetXYZ(0,0,0);                                                        
    m_p4Lepton[theChannel] -> SetXYZT(0,0,0,0);                                                        
    m_p3ChPFMet[theChannel].SetXYZ(0,0,0);
    modChPFMet[theChannel]  = 0.;
    WmT_ChPFMet[theChannel] = 0.;
    WmT_PFMet[ee] = 0.;
  }
}

int AxigluonSelection::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  float ETMax=0.;
  float ETMax2=0.;

  theLeadingJet[theChannel]=-1;   
  theSecondJet[theChannel] =-1;   

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);
    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    
    // PF jet ID variables - senti francesco....
    /*
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrEnergyAK5PFPUcorrJet[j];
    //     if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
    //                   chargedHadFraction,chargedMultiplicity,chargedEmFraction, Axi::loose)) continue;
    */
    
    bool foundMatch = false;

    // check if the electron falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {
      int ele = eleToRemove[i];
      if ( ele>-1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR = fabs( p3Jet.DeltaR(p3Ele) );
        if(_selectionEle->getSwitch("jetConeWidth") && _selectionEle->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    // check if the muon falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu>-1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR(p3Muon) );
        if(_selectionEle->getSwitch("jetConeWidth") && _selectionEle->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    if(_selectionEle->getSwitch("etaJetAcc") && !_selectionEle->passCut("etaJetAcc", fabs(etaAK5PFPUcorrJet[j]))) continue;

    if ( pt>ETMax2 && pt>ETMax ) {
      int theSecond = theLeadingJet[theChannel];
      int theLead   = j;
      theSecondJet[theChannel]  = theSecond;   
      theLeadingJet[theChannel] = theLead;    
      ETMax2 = ETMax;
      ETMax  = pt;
      subleadJetBtag[theChannel] = trackCountingHighEffBJetTagsAK5PFPUcorrJet[theSecond];   
      leadJetBtag[theChannel]    = trackCountingHighEffBJetTagsAK5PFPUcorrJet[theLead];     
    } else if ( pt>ETMax2 && pt<ETMax ) {
      theSecondJet[theChannel] = j;
      subleadJetBtag[theChannel] = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];
      ETMax2 = pt;
    }

    if(_selectionEle->getSwitch("etJetAcc") && !_selectionEle->passCut("etJetAcc", pt)) continue;

    num++;   
  }

  return num;
}

int AxigluonSelection::numSoftMuons(std::vector<int> muonToRemove) {

  int num = 0;
  for(int i=0; i<nMuon; ++i) {

    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) { 
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;

    float pt = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEle->getSwitch("softMuPt") && !_selectionEle->passCut("softMuPt",pt)) continue;  

    bool theMuonID = isAxigluonMuonID(i);   
    if (!theMuonID) continue;

    if (pt>20) {  // hardcoded
      float muonTrackerForGlobal = sumPt03Muon[i];
      float muonEcalForGlobal    = emEt03Muon[i];
      float muonHcalForGlobal    = hadEt03Muon[i];
      float theMuonGlobalSum     = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal - rhoFastjet*TMath::Pi()*0.3*0.3;
      float theRelMuonIso        = theMuonGlobalSum/pt; 
      if(_selectionEle->getSwitch("muGlobalIso") && !_selectionEle->passCut("muGlobalIso",theRelMuonIso)) continue;  
    }

    int ctfMuon = trackIndexMuon[i]; 
    float dxyMuon = transvImpactParTrack[ctfMuon];
    float dzMuon  = PVzPV[0] - trackVzTrack[ctfMuon];   
    if (_selectionEle->getSwitch("muonIP") && (!_selectionEle->passCut("muonIP",dxyMuon)) ) continue;   
    if (_selectionEle->getSwitch("muonDz") && (!_selectionEle->passCut("muonDz",dzMuon)) )  continue;   

    num++;
  }
  return num;
}

int AxigluonSelection::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    
    bool isSelEle=false;
    for(int eleSel=0; eleSel<(int)eleToRemove.size(); eleSel++) {
      if(i==eleToRemove[eleSel]) isSelEle=true;
    }
    if(isSelEle) continue;

    Utils anaUtils;
    bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[i], isEcalDriven);
    if(_selectionEle->getSwitch("ecalDrivenOnly") && !ecaldriven ) continue;    
    if(_selectionEle->getSwitch("etaElectronAcc") && !_selectionEle->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_selectionEle->getSwitch("electronBigCrack") && _selectionEle->passCut("electronBigCrack",fabs(etaEle[i])) ) continue;
    if(_selectionEle->getSwitch("ptElectronAcc")  && !_selectionEle->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])) ) continue;

    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;
    isEleIDAndDenom(i,&theId,&theIso,&theConvRej,&EgammaCutBasedID);
    if(!theId || !theIso || !theConvRej) continue;
    
    int track = gsfTrackIndexEle[i];
    float dxyEle = transvImpactParGsfTrack[track];
    float dzEle  = PVzPV[0] - trackVzGsfTrack[track];   
    if (_selectionEle->getSwitch("electronIP") && (!_selectionEle->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionEle->getSwitch("electronDz") && (!_selectionEle->passCut("electronDz",dzEle)) ) continue;

    numEle++;
  }

  int numMu = 0;
  for(int i=0; i<nMuon; ++i) {
    
    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) {
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;
    
    float ptMu = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEle->getSwitch("etaMuonAcc") && !_selectionEle->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_selectionEle->getSwitch("ptMuonAcc") && !_selectionEle->passCut("ptMuonAcc",ptMu) ) continue;

    bool theMuonID = isAxigluonMuonID(i);   
    if (!theMuonID) continue;

    float muonTrackerForGlobal = sumPt03Muon[i];
    float muonEcalForGlobal    = emEt03Muon[i];
    float muonHcalForGlobal    = hadEt03Muon[i];
    float theMuonGlobalSum     = muonTrackerForGlobal + muonEcalForGlobal + muonHcalForGlobal - rhoFastjet*TMath::Pi()*0.3*0.3;
    float theRelMuonIso        = theMuonGlobalSum/ptMu; 
    if(_selectionEle->getSwitch("muGlobalIso") && !_selectionEle->passCut("muGlobalIso",theRelMuonIso)) continue;  

    int ctfMuon = trackIndexMuon[i]; 
    float dxyMuon = transvImpactParTrack[ctfMuon];
    float dzMuon  = PVzPV[0] - trackVzTrack[ctfMuon];   
    if (_selectionEle->getSwitch("muonIP") && (!_selectionEle->passCut("muonIP",dxyMuon)) ) continue;   
    if (_selectionEle->getSwitch("muonDz") && (!_selectionEle->passCut("muonDz",dzMuon)) )  continue;   

    numMu++;
  }
  
  return numEle + numMu;
}

float AxigluonSelection::getDiJetInvMass(int theLJ, int theSJ) {

  TVector3 p3LJet(pxAK5PFPUcorrJet[theLJ],pyAK5PFPUcorrJet[theLJ],pzAK5PFPUcorrJet[theLJ]);
  TVector3 p3SJet(pxAK5PFPUcorrJet[theSJ],pyAK5PFPUcorrJet[theSJ],pzAK5PFPUcorrJet[theSJ]);
  TLorentzVector p4LJet(p3LJet, energyAK5PFPUcorrJet[theLJ]);
  TLorentzVector p4SJet(p3SJet, energyAK5PFPUcorrJet[theSJ]);
  return (p4LJet+p4SJet).M();
}

float AxigluonSelection::getDiJetPt(int theLJ, int theSJ) {

  TVector3 p3LJet(pxAK5PFPUcorrJet[theLJ],pyAK5PFPUcorrJet[theLJ],pzAK5PFPUcorrJet[theLJ]);
  TVector3 p3SJet(pxAK5PFPUcorrJet[theSJ],pyAK5PFPUcorrJet[theSJ],pzAK5PFPUcorrJet[theSJ]);
  TLorentzVector p4LJet(p3LJet, energyAK5PFPUcorrJet[theLJ]);
  TLorentzVector p4SJet(p3SJet, energyAK5PFPUcorrJet[theSJ]);
  return (p4LJet+p4SJet).Pt();
}

float AxigluonSelection::getDiJetDeta(int theLJ, int theSJ) {

  TVector3 p3LJet(pxAK5PFPUcorrJet[theLJ],pyAK5PFPUcorrJet[theLJ],pzAK5PFPUcorrJet[theLJ]);
  TVector3 p3SJet(pxAK5PFPUcorrJet[theSJ],pyAK5PFPUcorrJet[theSJ],pzAK5PFPUcorrJet[theSJ]);
  float eta1 = p3LJet.Eta();
  float eta2 = p3SJet.Eta();
  float deltaEta = eta1-eta2;
  return deltaEta;
}

void AxigluonSelection::setEleIdVariables(int eleIndex) {   

  myPt   = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
  myEta  = etaEle[eleIndex];
  myPhi  = phiEle[eleIndex];
  myDeta = deltaEtaAtVtxEle[eleIndex];
  myDphi = deltaPhiAtVtxEle[eleIndex];
  myHoe  = hOverEEle[eleIndex];
  mySee  = SigmaiEiE(eleIndex);
  myTrackerIso = dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
  myHcalIso    = dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
  myEcalIso    = dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3;
  float combinedIso = 0.0;
  if (fabs(myEta)<1.476) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  myCombinedIso = ( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / myPt );
  myCharge = chargeEle[eleIndex];
  myLh = eleIdLikelihoodEle[eleIndex];
}

// specific for HWW that has multiple channels with different HLT requirements
bool AxigluonSelection::reloadTriggerMask(int runN) {

  std::vector<int> triggerMask;
  
  // load the triggers required for Wenu
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersEle.begin();fIter!=requiredTriggersEle.end();++fIter) {  
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_requiredTriggersEle = triggerMask;

  // load the triggers NOT required for Wenu
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersEle.begin();fIter!=notRequiredTriggersEle.end();++fIter) {
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_notRequiredTriggersEle = triggerMask;

  // load the triggers required for Wmunu
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=requiredTriggersMu.begin();fIter!=requiredTriggersMu.end();++fIter) {
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_requiredTriggersMu = triggerMask;

  // load the triggers NOT required for Wmunu
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersMu.begin();fIter!=notRequiredTriggersMu.end();++fIter) {
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos) {
	triggerMask.push_back( indexHLT[i] ) ;
	break;
      }
    }
  }
  m_notRequiredTriggersMu = triggerMask;
}

bool AxigluonSelection::hasPassedHLT(int channel) {
  Utils anaUtils;
  if(channel==ee) {
    bool required    = anaUtils.getTriggersOR(m_requiredTriggersEle, firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersEle, firedTrg);
    return (required && !notRequired);
  } else if(channel==mm) {
    bool required    = anaUtils.getTriggersOR(m_requiredTriggersMu, firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersMu, firedTrg);
    return (required && !notRequired);
  }
  return true;
}

void AxigluonSelection::setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) requiredTriggersEle=reqTriggers;
  else if(channel==mm) requiredTriggersMu=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

void AxigluonSelection::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {
  if(channel==ee) notRequiredTriggersEle=reqTriggers;
  else if(channel==mm) notRequiredTriggersMu=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;
}

bool AxigluonSelection::isAxigluonMuonID(int muonIndex) {    

  bool isGood = true;

  bool isGlobal = (muonIdMuon[muonIndex] >> 13)%2;

  int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
  int normChi2 = trackNormalizedChi2GlobalMuonTrack[globalMuonTrack];
  int gTrackH = trackValidHitsGlobalMuonTrack[globalMuonTrack];
  int matchMu = numberOfMatchesMuon[muonIndex];
  if (!isGlobal)    isGood = false;
  if (normChi2>=10) isGood = false;
  if (gTrackH<=0)   isGood = false;
  if (matchMu<=1)   isGood = false;

  int track  = trackIndexMuon[muonIndex];
  int trackH = trackValidHitsTrack[track];
  int pixels = numberOfValidPixelBarrelHitsTrack[track] + numberOfValidPixelEndcapHitsTrack[track];
  if (trackH<=10 ) isGood = false;
  if (pixels<1 )   isGood = false;

  return isGood;
}

TVector3 AxigluonSelection::oneLepPFChargedMet(TVector3 lept) {   

  float chMetP3x = pxPFChMet[0];
  float chMetP3y = pyPFChMet[0];

  // charged PF MET has been computed with all the PF cands (inverted -p)                                                          
  // first remove the contribution in dR = 0.1 to avoid double counting                                                            
  for(int i=0; i<nReducedPFCand; i++) {
    TVector3 pfCandP3(pxReducedPFCand[i],pyReducedPFCand[i],pzReducedPFCand[i]);
    if(pfCandP3.DeltaR(lept)<=0.1) {
      chMetP3x += pxReducedPFCand[i];
      chMetP3y += pyReducedPFCand[i];
    }
  }

  // then add back the RECO leptons                                                                                                
  chMetP3x -= (lept.Px());
  chMetP3y -= (lept.Py());

  return TVector3(chMetP3x,chMetP3y,0.0);
}
