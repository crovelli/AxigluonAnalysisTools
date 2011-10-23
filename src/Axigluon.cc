#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <fstream>
#include <sstream>

#include "TString.h"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "AxigluonAnalysisTools/include/Axigluon.hh"

using namespace bits;

Axigluon::Axigluon(TTree *tree) : AxigluonBase(tree)
{
  jsonFile = "";
  lastFile = "";

  //initialize the JES objects for calo and PF
  jecUnc_calo = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/GR_R_42_V19_AK5Calo_Uncertainty.txt");
  jecUnc_PF = 
    (JetCorrectionUncertainty*) new JetCorrectionUncertainty("data/JES/GR_R_42_V19_AK5PF_Uncertainty.txt");

}

Axigluon::~Axigluon()
{
  // By this time, the destructor of AxigluonBase has not yet been called.
  // This means that the tree has not yet been deleted.
  // So, we do nothing here.
}

void Axigluon::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool Axigluon::hasPassedHLT() {
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

void Axigluon::setJsonGoodRunList(const std::string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void Axigluon::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool Axigluon::isGoodRunLS()
{
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //  std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //      std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
	return true;
    }
  return false;
}

bool Axigluon::reloadTriggerMask(bool newVersion)
{
  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {   
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
            if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
              {
                triggerMask.push_back( indexHLT[i] ) ;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {

        std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
        TTree *treeCond;
        std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
        std::vector<std::string>  *nameHLT_;
        std::vector<unsigned int> *indexHLT_;

        //To get the pointers for the vectors
        nameHLT_=0;
        indexHLT_=0;

        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);

        std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++) 
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
                  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (int i=0;i<m_requiredTriggers.size();++i)
          std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
}

float Axigluon::mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met) {
  float pTll = (pl1.Vect() + pl2.Vect()).Pt();
  float mll = (pl1 + pl2).M();
  float El = sqrt(pTll*pTll + mll*mll);
  float pTnu = met.Pt();
  float Enu = sqrt(pTnu*pTnu + mll*mll);
  float Ex = (pl1+pl2).X() + met.X();
  float Ey = (pl1+pl2).Y() + met.Y();
  float mnu = mll;

  return sqrt(mll*mll + mnu*mnu + 2*(El*Enu-Ex*Ex-Ey*Ey));
}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float Axigluon::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float Axigluon::SigmaiPiP(int electron) {
  float spp;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    spp = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      spp = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      spp = 999.;
    }
  }
  return spp;
}

bool Axigluon::isPFJetID(float eta, float neutralHadFrac, float neutralEmFraction, int nConstituents, float chargedHadFraction, 
                      float chargedMultiplicity, float chargedEmFraction, int WP) {
  switch(WP) {
  case none:
    return true;
    break;
  case loose:
    if(neutralHadFrac>=0.99 || neutralEmFraction>=0.99 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case medium:
    if(neutralHadFrac>=0.95 || neutralEmFraction>=0.95 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  case tight:
    if(neutralHadFrac>=0.90 || neutralEmFraction>=0.90 || nConstituents<=1) return false;
    if(fabs(eta)<2.4 && (chargedHadFraction==0 || chargedMultiplicity==0 || chargedEmFraction>=0.99) ) return false;
    break;
  default:
    std::cout << "Jet::isPFJetID(nt WP). Requested wrong Working point. Available are loose, medium, tight." << std::endl;
    return false;
  }
  return true;
}

/// *****************************************
/// From C. Rogan, M. Pierini, M. Spiropulu 
/// *****************************************

//
// this is the 'new' MRstar
//
double Axigluon::CalcMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  return temp;
}



//
// this is the 'new' MRstar, times 'gamma_{R*}' - I would recommend making 'R*' with this as 
// the denominator and 'M_{T}^{R}' as the numerator (the next function in this file)
//
double Axigluon::CalcGammaMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar
  temp *= mygamma;

  return temp;
}

//
// This is 'M_{T}^{R}', the guy that should be used in the numerator of 'R' or 'R*'
//
double Axigluon::CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){

  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;

  temp = sqrt(temp);

  return temp;
}

std::vector<int> Axigluon::sortElectronsByPt(std::vector<int> electrons) {
  int tmp;
  int max;
  for(int i=0;i<(int)electrons.size();i++) {
      max = i;
      float maxelePt = GetPt(pxEle[i],pyEle[i]);
      for(int x=i; x<(int)electrons.size(); x++) {
        float xelePt = GetPt(pxEle[x],pyEle[x]);
        if(xelePt > maxelePt) {
          max = x;
        }
      }
      tmp = electrons[i];
      electrons[i] = electrons[max];
      electrons[max] = tmp;
  }
  return electrons;
}

std::vector<int> Axigluon::sortMuonsByPt(std::vector<int> muons) {
  int tmp;
  int max;
  for(int i=0;i<(int)muons.size();i++) {
      max = i;
      float maxmuonPt = GetPt(pxMuon[i],pyMuon[i]);
      for(int x=i; x<(int)muons.size(); x++) {
        float xmuonPt = GetPt(pxMuon[x],pyMuon[x]);
        if(xmuonPt > maxmuonPt) {
          max = x;
        }
      }
      tmp = muons[i];
      muons[i] = muons[max];
      muons[max] = tmp;
  }
  return muons;
}

TLorentzVector Axigluon::GetJESCorrected(TLorentzVector p4jet, const char *ScaleDirection) {

  float mass = p4jet.M();
  float ptUnscaled = p4jet.Pt();
  
  // estimate the uncertainty
  jecUnc_PF->setJetEta(p4jet.Eta());
  jecUnc_PF->setJetPt(ptUnscaled);

  int scaleEnergy = 0;
  if(TString(ScaleDirection).Contains("Up")) scaleEnergy = 1.0;
  if(TString(ScaleDirection).Contains("Down")) scaleEnergy = -1.0;

  // apply the uncertainty
  float pt = ptUnscaled + scaleEnergy*jecUnc_PF->getUncertainty(true)*ptUnscaled;
  float p = pt/fabs(sin(p4jet.Theta()));
  float energy = sqrt(p*p+mass*mass);

  TLorentzVector p4Scaled;
  p4Scaled.SetPtEtaPhiE(pt,p4jet.Eta(),p4jet.Phi(),energy);

  return p4Scaled;
}

TVector3 Axigluon::pfChargedMet(TVector3 lep1, TVector3 lep2) {

  float chMetP3x = pxPFChMet[0];
  float chMetP3y = pyPFChMet[0];

  // charged PF MET has been computed with all the PF cands (inverted -p)
  // first remove the contribution in dR = 0.1 to avoid double counting
  for(int i=0; i<nReducedPFCand; i++) {
    TVector3 pfCandP3(pxReducedPFCand[i],pyReducedPFCand[i],pzReducedPFCand[i]);
    if(pfCandP3.DeltaR(lep1)<=0.1 || pfCandP3.DeltaR(lep2)<=0.1) {
      chMetP3x += pxReducedPFCand[i];
      chMetP3y += pyReducedPFCand[i];
    }
  }

  // then add back the RECO leptons
  chMetP3x -= (lep1.Px() + lep2.Px());
  chMetP3y -= (lep1.Py() + lep2.Py());
  
  return TVector3(chMetP3x,chMetP3y,0.0);

}

TVector3 Axigluon::corrSaclayMet(TVector3 lep1, TVector3 lep2) {

  float sMetP3x = pxsaclayPFMet[0];
  float sMetP3y = pysaclayPFMet[0];

  // saclay MET has been computed with all the PF cands (inverted -p)
  // first remove the contribution in dR = 0.1 to avoid double counting
  for(int i=0; i<nReducedPFCand; i++) {
    TVector3 pfCandP3(pxReducedPFCand[i],pyReducedPFCand[i],pzReducedPFCand[i]);
    if(pfCandP3.DeltaR(lep1)<=0.1 || pfCandP3.DeltaR(lep2)<=0.1) {
      sMetP3x += pxReducedPFCand[i];
      sMetP3y += pyReducedPFCand[i];
    }
  }

  // then add back the RECO leptons
  sMetP3x -= (lep1.Px() + lep2.Px());
  sMetP3y -= (lep1.Py() + lep2.Py());
  
  return TVector3(sMetP3x,sMetP3y,0.0);
}

std::string Axigluon::getHLTPathForRun(int runN, std::string fullname) {
  TString fullName = TString(fullname.c_str());
  TObjArray* selectionTokens = fullName.Tokenize(":");
  if (selectionTokens->GetEntries()!=2) {
    std::cout << "Wrong trigger strings " << selectionTokens->GetEntries() << std::endl;
    return std::string("NOPATH");
  }
  TString RunRange =((TObjString*)(*selectionTokens)[0])->GetString();
  TString HLTPathName =((TObjString*)(*selectionTokens)[1])->GetString();
  
  TObjArray* runs = RunRange.Tokenize("-");
  if (runs->GetEntries()!=2) {
    std::cout << "Wrong trigger run range strings " << runs->GetEntries() << std::endl;
    return std::string("NOPATH");    
  }
  
  const char *minStr = (((TObjString*)(*runs)[0])->GetString()).Data();
  const char *maxStr = (((TObjString*)(*runs)[1])->GetString()).Data();

  int min = atoi(minStr);
  int max = atoi(maxStr);

  if(runN>=min && runN<=max) return std::string(HLTPathName.Data());
  else return std::string("NOPATH");
}

bool Axigluon::isEleDenomFake(int theEle, bool *isDenomEleID, bool *isDenomEleIso, CutBasedEleIDSelector *thisCutBasedID) {

  Utils anaUtils;
  bool isGoodDenom = true;
  bool isGoodDenomID, isGoodDenomIso;
  isGoodDenomID = isGoodDenomIso = true;
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance	         
  if( fabs(p3Ele.Eta()) > 2.5 ) { isGoodDenom = false; isGoodDenomID = false; isGoodDenomIso = false; }  
  if( p3Ele.Pt() < 10.  )       { isGoodDenom = false; isGoodDenomID = false; isGoodDenomIso = false; }

  // barrel or endcap 
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);

  // taking shower shape
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  if ( ecalDriven) { 
    sc = superClusterIndexEle[theEle]; 
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]); 
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle]; 
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]); 
  }
  if ( sc < 0 ) { isGoodDenom = false; isGoodDenomID = false; isGoodDenomIso = false; }

  // sigmaIetaIeta
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; isGoodDenomID = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; isGoodDenomID = false; }

  // isolation
  // float ecalIsol    = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
  float ecalIsolAbs = 0.0;
  if ( isEleEB ) ecalIsolAbs = max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEtEle[theEle];
  float ecalIsol = ecalIsolAbs/p3Ele.Pt();

  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if (ecalIsol>0.2)    { isGoodDenom = false; isGoodDenomIso = false; }
  if (hcalIsol>0.2)    { isGoodDenom = false; isGoodDenomIso = false; }
  if (trackerIsol>0.2) { isGoodDenom = false; isGoodDenomIso = false; }

  // H/E 
  if ( isEleEB && hOverEEle[theEle]>0.12) { isGoodDenom = false; isGoodDenomID = false; }
  if (!isEleEB && hOverEEle[theEle]>0.10) { isGoodDenom = false; isGoodDenomID = false; }

  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) { isGoodDenom = false; isGoodDenomID = false; }
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) { isGoodDenom = false; isGoodDenomID = false; }

  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) { isGoodDenom = false; isGoodDenomID = false; }
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) { isGoodDenom = false; isGoodDenomID = false; }

  // full conversion rejection                                                                                                    
  bool dummyId, dummyIso, trueConvRej;
  dummyId = dummyIso = trueConvRej = true;
  isEleID(theEle,&dummyId,&dummyIso,&trueConvRej,thisCutBasedID);  
  if (!trueConvRej) isGoodDenomID = false;

  *isDenomEleID = isGoodDenomID;
  *isDenomEleIso = isGoodDenomIso;

  return isGoodDenom;
}

bool Axigluon::isMuonDenomFake(int theMuon, bool *isDenomMuonID, bool *isDenomMuonIso) {

  bool isGoodDenom = true;
  bool isGoodDenomID, isGoodDenomIso;
  isGoodDenomID = isGoodDenomIso = true;
  TVector3 p3Muon(pxMuon[theMuon], pyMuon[theMuon], pzMuon[theMuon]);
  
  // acceptance	   
  if( fabs(p3Muon.Eta()) > 2.5 ) { isGoodDenom = false; isGoodDenomID = false; isGoodDenomIso = false; }
  if( p3Muon.Pt() < 10. )        { isGoodDenom = false; isGoodDenomID = false; isGoodDenomIso = false; }

  // muonID
  bool isTight = true;
  isMuonID(theMuon, &isTight);
  if (!isTight) { isGoodDenom = false; isGoodDenomID = false; }
  
  // isolation
  float thePFMuonIso = pfCombinedIsoMuon[theMuon]/p3Muon.Pt();
  if ( thePFMuonIso > 0.4 ) { isGoodDenom = false; isGoodDenomIso = false; }   // this is LM2; LM1 < 1.0
  
  // IP
  int ctfMuon   = trackIndexMuon[theMuon]; 
  float dxyMuon = transvImpactParTrack[ctfMuon];
  float dzMuon  = PVzPV[0] - trackVzTrack[ctfMuon];  
  if (fabs(dxyMuon)>0.1 ) { isGoodDenom = false; isGoodDenomID = false; }     // this is LM2; LM1 < 0.1
  if (fabs(dzMuon)>0.1  ) { isGoodDenom = false; isGoodDenomID = false; }     // this is LM2; LM1 < 0.1
  
  *isDenomMuonID = isGoodDenomID;
  *isDenomMuonIso = isGoodDenomIso;

  return isGoodDenom;
}

void Axigluon::isMuonID(int muonIndex, bool *muonIdOutput) {

  *muonIdOutput = true;

  Utils anaUtils; 
  bool flagGlobalMu = false;
  if(anaUtils.muonIdVal(muonIdMuon[muonIndex],AllGlobalMuons)) {
    int globalMuonTrack = combinedTrackIndexMuon[muonIndex];
    if(trackNormalizedChi2GlobalMuonTrack[globalMuonTrack] < 10 && 
       trackValidHitsGlobalMuonTrack[globalMuonTrack] > 0 &&
       numberOfMatchesMuon[muonIndex] > 1 ) flagGlobalMu = true; // to be used when new trees are available
  }

  bool flagTrackerMu = false;
  if( (anaUtils.muonIdVal(muonIdMuon[muonIndex],AllTrackerMuons) &&
       anaUtils.muonIdVal(muonIdMuon[muonIndex],TMLastStationTight)) ) flagTrackerMu  = true;

  if(!(flagGlobalMu || flagTrackerMu)) {
    *muonIdOutput = false;
    return;
  }
    
  int track = trackIndexMuon[muonIndex];

  if(trackValidHitsTrack[track]<=10) *muonIdOutput = false;

  if( (numberOfValidPixelBarrelHitsTrack[track]+numberOfValidPixelEndcapHitsTrack[track])<1 ) *muonIdOutput = false; 

  float ptTrack = sqrt( pxTrack[track]*pxTrack[track] + pyTrack[track]*pyTrack[track] );
  float sign = fabs(ptErrorTrack[track]/ptTrack);
  if (sign>=0.1) *muonIdOutput = false;
}

void Axigluon::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {
 
  *eleIdOutput = *isolOutput = *convRejOutput = false; 

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  float pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxSC[sc];
      e4SwissCross = e4SwissCrossSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagSC[sc];
      seedTime = timeSC[sc];
      seedChi2 = chi2SC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  
  thisCutBasedID->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  thisCutBasedID->SetRecoFlag(recoFlagsEle[eleIndex]);
  thisCutBasedID->applyElectronIDOnPFlowElectrons(true);
  thisCutBasedID->SetHOverE( HoE );
  thisCutBasedID->SetS9S25( s9s25 );
  thisCutBasedID->SetDEta( deta );
  thisCutBasedID->SetDPhiIn( dphiin );
  thisCutBasedID->SetDPhiOut( dphiout );
  thisCutBasedID->SetBremFraction( fbrem );
  thisCutBasedID->SetSigmaEtaEta( see );
  thisCutBasedID->SetSigmaPhiPhi( spp );
  thisCutBasedID->SetEOverPout( eopout );
  thisCutBasedID->SetEOverPin( eop );
  thisCutBasedID->SetElectronClass ( classificationEle[eleIndex] );
  //  thisCutBasedID->SetLikelihood( likelihoodRatio(eleIndex,*LH) );
  thisCutBasedID->SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  thisCutBasedID->SetNBrem( nbremsEle[eleIndex] );
  thisCutBasedID->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                
  thisCutBasedID->SetTrkIsolation ( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );                        
  thisCutBasedID->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pt );         
  float iso = 0.0;
  if ( anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB) ) iso = dr03TkSumPtEle[eleIndex] + max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else iso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  thisCutBasedID->SetCombinedIsolation( (iso - rhoFastjet*TMath::Pi()*0.3*0.3) / pt );
  thisCutBasedID->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pt );
  thisCutBasedID->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  thisCutBasedID->SetConvDist( fabs(convDistEle[eleIndex]) );
  thisCutBasedID->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  thisCutBasedID->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );

  // ECAL cleaning variables
  thisCutBasedID->m_cleaner->SetE1(e1);
  thisCutBasedID->m_cleaner->SetE4SwissCross(e4SwissCross);
  thisCutBasedID->m_cleaner->SetFiducialFlag(fidFlagSC);
  thisCutBasedID->m_cleaner->SetSeedFlag(seedRecHitFlag);
  thisCutBasedID->m_cleaner->SetSeedTime(seedTime);
  thisCutBasedID->m_cleaner->SetSeedChi2(seedChi2);

  //  return egammaCutBasedID.output(); // class dependent result
  *eleIdOutput = thisCutBasedID->outputNoClassEleId();
  *isolOutput = thisCutBasedID->outputNoClassIso();
  *convRejOutput = thisCutBasedID->outputNoClassConv();
}

void Axigluon::isEleIDAndDenom(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, CutBasedEleIDSelector *thisCutBasedID) {

  bool tightId, tightIso, tightConvRej;
  tightId = tightIso = tightConvRej = true;
  isEleID(eleIndex,&tightId,&tightIso,&tightConvRej,thisCutBasedID);

  bool denomId, denomIso;
  denomId = denomIso = true;
  isEleDenomFake(eleIndex,&denomId,&denomIso,thisCutBasedID); 
  
  // denominator definition is only applied on ID (because different algorithm is used offline and in HLT)
  *eleIdOutput = (tightId && denomId);
  *isolOutput = tightIso;
  *convRejOutput = tightConvRej;

}
