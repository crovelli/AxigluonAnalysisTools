#include "AxigluonAnalysisTools/include/CutBasedAxiSelector.hh"
#include <iostream>
#include <math.h>

CutBasedAxiSelector::CutBasedAxiSelector() {

  // all steps 
  m_step0 = false;
  m_step1 = false;
  m_step2 = false;
  m_step3 = false;
  m_step4 = false;
  m_step5 = false;
  m_step6 = false;
  m_step7 = false;
  m_step8 = false;
  m_step9 = false;
  m_step10 = false;
  m_step11 = false;
  m_step12 = false;
  m_step13 = false;
  m_step14 = false;
  m_step15 = false;
  m_step16 = false;

  m_processID = -1;
}

CutBasedAxiSelector::CutBasedAxiSelector( const CutBasedAxiSelector& selector ) {

  m_weight = selector.m_weight;
  m_passedHLT = selector.m_passedHLT;
  m_isThisChannel = selector.m_isThisChannel;
  m_isElectronId = selector.m_isElectronId;
  m_isElectronIsol = selector.m_isElectronIsol;
  m_isElectronConvRej = selector.m_isElectronConvRej;
  m_isElectronIp = selector.m_isElectronIp;
  m_nExtraLeptons = selector.m_nExtraLeptons;
  m_leptPt = selector.m_leptPt;  
  m_met = selector.m_met;
  m_wmt = selector.m_wmt;       
  m_nJets = selector.m_nJets;
  m_dijetInvMass = selector.m_dijetInvMass;  
  m_dijetPt = selector.m_dijetPt;        
  m_dijetDeta = selector.m_dijetDeta;       
  m_leadJetBtag = selector.m_leadJetBtag;   
  m_subleadJetBtag = selector.m_subleadJetBtag; 
  m_nSoftMuons = selector.m_nSoftMuons;
  multiProcessCounter = selector.multiProcessCounter;  

  // steps
  m_step0  = selector.m_step0;
  m_step1  = selector.m_step1;
  m_step2  = selector.m_step2;
  m_step3  = selector.m_step3;
  m_step4  = selector.m_step4;
  m_step5  = selector.m_step5;
  m_step6  = selector.m_step6;
  m_step7  = selector.m_step7;
  m_step8  = selector.m_step8;
  m_step9  = selector.m_step9;
  m_step10 = selector.m_step10;
  m_step11 = selector.m_step11;
  m_step12 = selector.m_step12;
  m_step13 = selector.m_step13;
  m_step14 = selector.m_step14;
  m_step15 = selector.m_step15;
  m_step16 = selector.m_step16;
}

CutBasedAxiSelector::~CutBasedAxiSelector() {}

void CutBasedAxiSelector::Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle) {

  _selection = new Selection(std::string(fileCuts),std::string(fileSwitches));

  // these cuts are applied in the AxiSelection class, but are configured here
  _selection->addSwitch("MCtruth"); 
  _selection->addSwitch("trigger"); 
  _selection->addSwitch("leptonId");  
  _selection->addSwitch("leptonIso"); 
  _selection->addSwitch("leptonD0");  
  _selection->addSwitch("convRej");   
  _selection->addCut("muGlobalIso");  
  _selection->addCut("electronIP");
  _selection->addCut("electronDz");
  _selection->addCut("muonIP");
  _selection->addCut("muonDz");
  _selection->addCut("ptLepton");   
  _selection->addCut("met");   
  _selection->addCut("wmt");   
  _selection->addCut("nExtraLeptons"); 
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("jetConeWidth");
  _selection->addCut("nJets");
  _selection->addCut("dijetInvMass"); 
  _selection->addCut("dijetPt");      
  _selection->addCut("dijetDeta");    
  _selection->addCut("leadJetBtag"); 
  _selection->addCut("subleadJetBtag");
  _selection->addCut("nSoftMuons"); 
  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle(theTitle);
  globalCounter->AddVar("event"); 
  globalCounter->AddVar("MCtruth"); 
  globalCounter->AddVar("trigger"); // 0
  globalCounter->AddVar("preselected"); // 1
  globalCounter->AddVar("leptonId"); // 2
  globalCounter->AddVar("leptonIso"); // 3
  globalCounter->AddVar("convRej"); // 4
  globalCounter->AddVar("leptonD0"); // 5   
  globalCounter->AddVar("ptLepton"); // 6
  globalCounter->AddVar("met"); // 7
  globalCounter->AddVar("wmt"); // 8
  globalCounter->AddVar("nExtraLeptons"); // 9
  globalCounter->AddVar("nJets"); // 10
  globalCounter->AddVar("dijetInvMass"); // 11
  globalCounter->AddVar("dijetPt"); // 12
  globalCounter->AddVar("dijetDeta"); // 13
  globalCounter->AddVar("leadJetBtag"); // 14
  globalCounter->AddVar("subleadJetBtag"); // 15
  globalCounter->AddVar("nSoftMuons"); // 16
}

bool CutBasedAxiSelector::output() {

  Counters *theCounter=0;

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter = multiProcessCounter.find(m_processID);

    if ( iter == multiProcessCounter.end() ) {
      
      std::cout << "First time I get process " << m_processID 
		<< ": adding a counter" << std::endl;

      char buffer[200];
      sprintf(buffer,"Event counter for process %d", m_processID);
      
      Counters *processCounter = new Counters();
      processCounter->SetTitle(buffer);
      processCounter->AddVar("event");
      processCounter->AddVar("MCtruth");
      processCounter->AddVar("trigger");
      processCounter->AddVar("preselected");
      processCounter->AddVar("leptonId");
      processCounter->AddVar("leptonIso");      
      processCounter->AddVar("convRej");
      processCounter->AddVar("leptonD0");
      processCounter->AddVar("ptLepton");
      processCounter->AddVar("met");
      processCounter->AddVar("wmt");
      processCounter->AddVar("nExtraLeptons");
      processCounter->AddVar("nJets");
      processCounter->AddVar("dijetInvMass");
      processCounter->AddVar("dijetPt");
      processCounter->AddVar("dijetDeta");
      processCounter->AddVar("nSoftMuons");
      processCounter->AddVar("leadJetBtag");
      processCounter->AddVar("subleadJetBtag");
      multiProcessCounter.insert( std::make_pair(m_processID,processCounter) );      
    }

    theCounter = multiProcessCounter[m_processID];

  }
  
  else theCounter = globalCounter;
  
  // steps
  m_step0 = false;
  m_step1 = false;
  m_step2 = false;
  m_step3 = false;
  m_step4 = false;
  m_step5 = false;
  m_step6 = false;
  m_step7 = false;
  m_step8 = false;
  m_step9 = false;
  m_step10 = false;
  m_step11 = false;
  m_step12 = false;
  m_step13 = false;
  m_step14 = false;
  m_step15 = false;
  m_step16 = false;

  theCounter->IncrVar("event",m_weight);
  
  if (_selection->getSwitch("MCtruth") && !m_foundMcTree) return false;
  theCounter->IncrVar("MCtruth",m_weight);

  if(_selection->getSwitch("trigger") && !m_passedHLT ) return false;
  theCounter->IncrVar("trigger",m_weight); 
  m_step0 = true;

  if(!m_isThisChannel) return false;
  theCounter->IncrVar("preselected",m_weight);
  m_step1 = true;
  
  if (_selection->getSwitch("leptonId") && (m_isElectronId<0) ) return false; 
  theCounter->IncrVar("leptonId",m_weight);
  m_step2 = true;
  
  if (_selection->getSwitch("leptonIso") && (m_isElectronIsol<0)) return false; 
  theCounter->IncrVar("leptonIso",m_weight);
  m_step3 = true;

  if (_selection->getSwitch("convRej") && (m_isElectronConvRej<0) ) return false; 
  theCounter->IncrVar("convRej",m_weight);
  m_step4 = true;

  if (_selection->getSwitch("leptonD0") && (m_isElectronIp<0) ) return false; 
  theCounter->IncrVar("leptonD0",m_weight);
  m_step5 = true;
  
  if (_selection->getSwitch("ptLepton") && !_selection->passCut("ptLepton", m_leptPt)) return false;
  theCounter->IncrVar("ptLepton",m_weight);
  m_step6 = true;

  if (_selection->getSwitch("met") && !_selection->passCut("met",m_met)) return false; 
  theCounter->IncrVar("met",m_weight);
  m_step7 = true;

  if (_selection->getSwitch("wmt") && !_selection->passCut("wmt",m_wmt)) return false; 
  theCounter->IncrVar("wmt",m_weight);
  m_step8 = true;

  if (_selection->getSwitch("nExtraLeptons") && !_selection->passCut("nExtraLeptons",m_nExtraLeptons)) return false;
  theCounter->IncrVar("nExtraLeptons",m_weight);
  m_step9 = true;

  if (_selection->getSwitch("nJets") && !_selection->passCut("nJets",m_nJets)) return false;
  theCounter->IncrVar("nJets",m_weight);
  m_step10 = true;

  if (_selection->getSwitch("dijetInvMass") && !_selection->passCut("dijetInvMass",m_dijetInvMass)) return false;
  theCounter->IncrVar("dijetInvMass",m_weight);
  m_step11 = true;

  if (_selection->getSwitch("dijetPt") && !_selection->passCut("dijetPt",m_dijetPt)) return false;
  theCounter->IncrVar("dijetPt",m_weight);
  m_step12 = true;

  if (_selection->getSwitch("dijetDeta") && !_selection->passCut("dijetDeta",m_dijetDeta)) return false;
  theCounter->IncrVar("dijetDeta",m_weight);
  m_step13 = true;

  if (_selection->getSwitch("nSoftMuons") && !_selection->passCut("nSoftMuons",m_nSoftMuons)) return false;
  theCounter->IncrVar("nSoftMuons",m_weight);
  m_step14 = true;

  if (_selection->getSwitch("leadJetBtags") && !_selection->passCut("leadJetBtags",m_leadJetBtag)) return false;
  theCounter->IncrVar("leadJetBtags",m_weight);
  m_step15 = true;

  if (_selection->getSwitch("subleadJetBtags") && !_selection->passCut("subleadJetBtags",m_subleadJetBtag)) return false;
  theCounter->IncrVar("subleadJetBtags",m_weight);
  m_step16 = true;

  return true;
}


void CutBasedAxiSelector::displayEfficiencies(std::string datasetName) {

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter;
    for( iter=multiProcessCounter.begin(); iter!=multiProcessCounter.end(); ++iter ) {

      Counters *theCounter = iter->second;

      theCounter->Draw();
      theCounter->Draw("trigger","MCtruth");
      theCounter->Draw("preselected","trigger");
      theCounter->Draw("leptonId","preselected");     
      theCounter->Draw("leptonIso","leptonId");
      theCounter->Draw("convRej","leptonIso");
      theCounter->Draw("leptonD0","convRej");
      theCounter->Draw("ptLepton","leptonD0");
      theCounter->Draw("met","ptLepton");
      theCounter->Draw("wmt","met");
      theCounter->Draw("nExtraLeptons","wmt");
      theCounter->Draw("nJets","nExtraLeptons");
      theCounter->Draw("dijetInvMass","nJets");
      theCounter->Draw("dijetInvPt","dijetInvMass");
      theCounter->Draw("dijetInvDeta","dijetPt");
      theCounter->Draw("nSoftMuons","dijetInvDeta");
      theCounter->Draw("leadJetBtag","nSoftMuons");
      theCounter->Draw("subleadJetBtag","leadJetBtag");
      theCounter->Draw("subleadJetBtag","preselected");
    }
  }

  else {

    char namefile[500];
    sprintf(namefile,"%s-Counters.root",datasetName.c_str());
    
    globalCounter->Draw();
    globalCounter->Draw("trigger","MCtruth");
    globalCounter->Draw("preselected","trigger");
    globalCounter->Draw("leptonId","preselected");     
    globalCounter->Draw("leptonIso","leptonId");
    globalCounter->Draw("convRej","leptonIso");
    globalCounter->Draw("leptonD0","convRej");
    globalCounter->Draw("ptLepton","leptonD0");
    globalCounter->Draw("met","ptLepton");
    globalCounter->Draw("wmt","met");
    globalCounter->Draw("nExtraLeptons","wmt");
    globalCounter->Draw("nJets","nExtraLeptons");
    globalCounter->Draw("dijetInvMass","nJets");
    globalCounter->Draw("dijetInvPt","dijetInvMass");
    globalCounter->Draw("dijetInvDeta","dijetPt");
    globalCounter->Draw("nSoftMuons","dijetInvDeta");
    globalCounter->Draw("leadJetBtag","nSoftMuons");
    globalCounter->Draw("subleadJetBtag","leadJetBtag");
    globalCounter->Draw("subleadJetBtag","preselected");
    
    globalCounter->Save(namefile,"update");
  }

}
