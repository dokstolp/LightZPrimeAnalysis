// -*- C++ -*-
//
// Package:    LightZPrimeAnalysis/JetAnalyzerMC
// Class:      JetAnalyzerMC
// 
/**\class JetAnalyzerMC JetAnalyzerMC.cc LightZPrimeAnalysis/JetAnalyzerMC/plugins/JetAnalyzerMC.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Tue, 23 Feb 2016 04:57:10 GMT
// Second Author: Usama Hussain
//


// system include files
#include <memory>
#include <vector>
#include <list>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "LightZPrimeAnalysis/JetWidthCalculator/interface/JetWidthCalculator.hh"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class JetAnalyzerMC : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JetAnalyzerMC(const edm::ParameterSet&);
      ~JetAnalyzerMC();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      // ----------member data ---------------------------

  edm::EDGetTokenT< vector<reco::PFCandidate> > pfCandsToken;
  edm::EDGetTokenT< vector<reco::PFJet> > pfJetsToken;
  edm::EDGetTokenT< vector<reco::PFMET> > pfMETsToken;
  edm::EDGetToken electronCollection_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >pileupSummaryInfoToken_;

  // elecontr ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;
 
  //genlevelInfo Weights
  float genWeight_;
  int   numWeights_;
  vector<float>  genWeights_;
  float nTrueVertices_;
  int   nup_;
  int   numGenJets_; 
  float genHT_; 

  //some must have variables for tuples
  int     run_;
  int  event_;
  int     lumis_; 
  
  //jet variables
  vector<float> jetPt_;
  vector<float> jetEn_;
  vector<float> jetEta_;
  vector<float> jetPhi_;
  vector<float> jetCHF_;
  vector<float> jetNHF_;
  vector<float> jetCEF_;
  vector<float> jetNEF_;
  vector<int>   jetNCH_;
  vector<float> jetHFHAE_;
  vector<float> jetHFEME_;
  vector<int>   jetNConstituents_;
  vector<float> jetEtaWidth_;
  vector<float> jetPhiWidth_;
  vector<float> jetEtaWidthInECal_;
  vector<float> jetEtaWidthInHCal_;
  vector<float> jetPhiWidthInECal_;
  vector<float> jetPhiWidthInHCal_;
  
  vector<bool>  jetPFLooseId_;
 
  //electron variables
  Int_t          nEle_;
  vector<float>  elePt_;
  vector<float>  eleEta_;
  vector<float>  elePhi_;
  int PassVeto_;
  int PassLoose_;
  int PassMedium_;
  int PassTight_;
  int PassHEEP_;
 
  double totalET;
  double HT;
  float METValue;
  float  METPhi;
  float  METsumEt_;
  float  METmEtSig_;
  float  METSig_;
  double j1PT;
  double j1Eta;
  double j1Phi;
  double j1CHdFr;
  double j1NHdFr;
  double j1CEmFr;
  double j1NEmFr;
  double j1etaWidth;
  double j1phiWidth;
  double j1etaWidthInECal;
  double j1phiWidthInECal;
  double j1etaWidthInHCal;
  double j1phiWidthInHCal;
  double j2PT;
  double j2Eta;
  double j2Phi;
  double j2CHdFr;
  double j2NHdFr;
  double j2CEmFr;
  double j2NEmFr;
  double j2etaWidth;
  double j2phiWidth;
  double j2etaWidthInECal;
  double j2phiWidthInECal;
  double j2etaWidthInHCal;
  double j2phiWidthInHCal;

  uint32_t nJets;
  uint32_t nGoodJets;

  uint32_t j1nCons;
  uint32_t j2nCons;
  
  TTree* tree;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JetAnalyzerMC::JetAnalyzerMC(const edm::ParameterSet& iConfig)
{

  genEventInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  pfCandsToken = consumes< vector<reco::PFCandidate> >(edm::InputTag("particleFlow"));
  pfJetsToken = consumes< vector<reco::PFJet> >(edm::InputTag("ak4PFJetsCHS"));
  pfMETsToken = consumes< vector<reco::PFMET> >(edm::InputTag("pfMet"));
  electronCollection_ = mayConsume<edm::View<reco::GsfElectron> >(edm::InputTag("gedGsfElectrons"));
  lheEventProductToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  genEventInfoToken_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));  
  pileupSummaryInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));

  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-veto"));
  eleLooseIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-loose"));
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-medium"));
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-tight")); 
  eleHEEPIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"));  
  
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("JetTree", "Jet data for analysis");

  tree->Branch("run",     &run_);
  tree->Branch("event",   &event_);
  tree->Branch("lumis",   &lumis_);  
  //lhe and genlevel Info  
  tree->Branch("genEventWeight", &genWeight_, "genWeight/F"); //genEventInfo->weight();
  tree->Branch("genlheWeights", &genWeights_); // genWeights_.push_back(lheInfo->weights()[i].wgt);
  tree->Branch("numWeights", &numWeights_, "numWeights/I"); //genWeights_.size();
  tree->Branch("nTrueVertices", &nTrueVertices_, "nTrueVertices/F");
  tree->Branch("NUP", &nup_, "NUP/I");
  tree->Branch("numGenJets", &numGenJets_, "numGenJets/I");
  tree->Branch("genHT", &genHT_, "genHT/F"); 
 
  tree->Branch("totalET", &totalET);
  tree->Branch("HT", &HT);
  tree->Branch("METValue", &METValue);
  tree->Branch("METPhi", &METPhi);
  tree->Branch("METsumEt",       &METsumEt_);
  tree->Branch("METmEtSig",      &METmEtSig_);
  tree->Branch("METSig",         &METSig_);
  tree->Branch("nJets", &nJets);
  tree->Branch("jetPt",&jetPt_);
  tree->Branch("jetEta",          &jetEta_);
  tree->Branch("jetEn",           &jetEn_);
  tree->Branch("jetPhi",          &jetPhi_);
  tree->Branch("jetPFLooseId", &jetPFLooseId_);
  tree->Branch("jetCHF",       &jetCHF_);
  tree->Branch("jetNHF",       &jetNHF_);
  tree->Branch("jetCEF",       &jetCEF_);
  tree->Branch("jetNEF",       &jetNEF_);
  tree->Branch("jetNCH",       &jetNCH_); 
  tree->Branch("jetHFHAE",         &jetHFHAE_);
  tree->Branch("jetHFEME",         &jetHFEME_);
  tree->Branch("jetNConstituents", &jetNConstituents_);
  tree->Branch("jetEtaWidth", &jetEtaWidth_);
  tree->Branch("jetPhiWidth", &jetPhiWidth_);
  tree->Branch("jetEtaWidthInECal", &jetEtaWidthInECal_);
  tree->Branch("jetEtaWidthInHCal", &jetEtaWidthInHCal_);
  tree->Branch("jetPhiWidthInECal", &jetPhiWidthInECal_);
  tree->Branch("jetPhiWidthInHCal", &jetPhiWidthInHCal_); 
  tree->Branch("nEle",                    &nEle_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_); 
  tree->Branch("ElectronPassVetoID", &PassVeto_);
  tree->Branch("ElectronPassLooseID", &PassLoose_);
  tree->Branch("ElectronPassMediumID", &PassMedium_);
  tree->Branch("ElectronPassTightID", &PassTight_);
  tree->Branch("ElectronPassHEEPID", &PassHEEP_);
  
  tree->Branch("nGoodJets", &nGoodJets);
  tree->Branch("j1PT", &j1PT);
  tree->Branch("j1Eta", &j1Eta);
  tree->Branch("j1Phi", &j1Phi);
  tree->Branch("j1CHdFr", &j1CHdFr);
  tree->Branch("j1NHdFr", &j1NHdFr);
  tree->Branch("j1CEmFr", &j1CEmFr);
  tree->Branch("j1NEmFr", &j1NEmFr);
  tree->Branch("j1nCons", &j1nCons);
  tree->Branch("j1etaWidth", &j1etaWidth);
  tree->Branch("j1phiWidth", &j1phiWidth);
  tree->Branch("j1etaWidthInECal", &j1etaWidthInECal);
  tree->Branch("j1phiWidthInECal", &j1phiWidthInECal);
  tree->Branch("j1etaWidthInHCal", &j1etaWidthInHCal);
  tree->Branch("j1phiWidthInHCal", &j1phiWidthInHCal);
  tree->Branch("j2PT", &j2PT);
  tree->Branch("j2Eta", &j2Eta);
  tree->Branch("j2Phi", &j2Phi);
  tree->Branch("j2CHdFr", &j2CHdFr);
  tree->Branch("j2NHdFr", &j2NHdFr);
  tree->Branch("j2CEmFr", &j2CEmFr);
  tree->Branch("j2NEmFr", &j2NEmFr);
  tree->Branch("j2nCons", &j2nCons);
  tree->Branch("j2etaWidth", &j2etaWidth);
  tree->Branch("j2phiWidth", &j2phiWidth);
  tree->Branch("j2etaWidthInECal", &j2etaWidthInECal);
  tree->Branch("j2phiWidthInECal", &j2phiWidthInECal);
  tree->Branch("j2etaWidthInHCal", &j2etaWidthInHCal);
  tree->Branch("j2phiWidthInHCal", &j2phiWidthInHCal);
}


JetAnalyzerMC::~JetAnalyzerMC()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
// ------------ method called for each event  ------------
void
JetAnalyzerMC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle< vector<reco::PFCandidate> > pfCands;
   iEvent.getByToken(pfCandsToken, pfCands);

   Handle< vector<reco::PFJet> > pfJets;
   iEvent.getByToken(pfJetsToken, pfJets);

   Handle< vector<reco::PFMET> > pfMETs;
   iEvent.getByToken(pfMETsToken, pfMETs);
   
   run_    = iEvent.id().run();
   event_  = iEvent.id().event();
   lumis_  = iEvent.luminosityBlock();   

   nEle_ = 0;
   
   Handle<edm::View<reco::GsfElectron> > electronHandle;
   iEvent.getByToken(electronCollection_, electronHandle);

   Handle<edm::ValueMap<bool> >  veto_id_decisions;
   iEvent.getByToken(eleVetoIdMapToken_ ,         veto_id_decisions); 
   Handle<edm::ValueMap<bool> >  loose_id_decisions;
   iEvent.getByToken(eleLooseIdMapToken_ ,        loose_id_decisions);
   Handle<edm::ValueMap<bool> >  medium_id_decisions;
   iEvent.getByToken(eleMediumIdMapToken_,        medium_id_decisions); 
   Handle<edm::ValueMap<bool> >  tight_id_decisions; 
   iEvent.getByToken(eleTightIdMapToken_,         tight_id_decisions);
   Handle<edm::ValueMap<bool> >  heep_id_decisions;
   iEvent.getByToken(eleHEEPIdMapToken_ ,         heep_id_decisions);
   
   //Clear previous events
   jetPt_.clear();
   jetEn_.clear();
   jetEta_.clear();
   jetPhi_.clear();
   jetPFLooseId_                           .clear();
   jetCHF_                                 .clear();
   jetNHF_                                 .clear();
   jetCEF_                                 .clear();
   jetNEF_                                 .clear();
   jetNCH_                                 .clear();
   jetHFHAE_                               .clear();
   jetHFEME_                               .clear();
   jetNConstituents_                       .clear();   
   jetEtaWidth_.clear();
   jetPhiWidth_.clear();
   jetEtaWidthInECal_.clear();
   jetPhiWidthInECal_.clear();
   jetEtaWidthInHCal_.clear();
   jetPhiWidthInHCal_.clear();
   
   elePt_                      .clear();
   eleEta_                     .clear();
   elePhi_                     .clear();
  
  //Treatment of the negative event weight (especially in MG5_aMCNLO)
  //weight can be retrieved by GenEventInfoProduct. In each event, we can take it as followings
  
  Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(genEventInfoToken_, genEventInfo);
  
  genWeight_ = 0.;
  if (genEventInfo.isValid()) {
      genWeight_ = genEventInfo->weight();
  }

  //PileupSummaryInfo
  Handle<std::vector<PileupSummaryInfo> > pileupSummaryInfo;
  iEvent.getByToken(pileupSummaryInfoToken_, pileupSummaryInfo);   

  nTrueVertices_ = 0;
  if (pileupSummaryInfo.isValid() && pileupSummaryInfo->size()>0) {
      nTrueVertices_ = pileupSummaryInfo->at(1).getTrueNumInteractions();
  }
  
  //calculate lhe information for weighting and accessing gen level HT 
  nup_ = 0;
  numGenJets_ = 0; //number of outgoing partons
  genHT_ = 0;
  genWeights_.clear();

  Handle<LHEEventProduct> lheEventProduct;  
  iEvent.getByToken( lheEventProductToken_, lheEventProduct);
 
  if (lheEventProduct.isValid()){ 
  	const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
	std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
  	nup_ = lheEvent.NUP;
  	for ( int idxParticle = 0; idxParticle < nup_; ++idxParticle ) {
   		int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
   		int status = lheEvent.ISTUP[idxParticle];
                if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
       			genHT_+= TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
       			++numGenJets_;
   		} 
  	}
	for (size_t i=0; i<lheEventProduct->weights().size(); ++i) {
		genWeights_.push_back(lheEventProduct->weights()[i].wgt);
	}
  }  
  numWeights_ = genWeights_.size();
   // Set event level quantities

   totalET = 0;
   HT = 0;
   nJets = 0;
   nGoodJets = 0;

   // Compute HT
   
   for(uint32_t j = 0; j < pfJets->size(); j++) {
     const pat::Jet &jet = (*pfJets)[j];
     if(j == 0) {
       j1PT = jet.pt();
       j1Eta = jet.eta();
       j1Phi = jet.phi();
       j1CHdFr = jet.chargedHadronEnergyFraction();
       j1NHdFr = jet.neutralHadronEnergyFraction();
       j1CEmFr = jet.chargedEmEnergyFraction();
       j1NEmFr = jet.neutralEmEnergyFraction();
       j1nCons = jet.nConstituents();
       JetWidthCalculator jwc(jet);
       j1etaWidth = jwc.getEtaWidth();
       j1phiWidth = jwc.getPhiWidth();
       j1etaWidthInECal = jwc.getEtaWidthInECal();
       j1phiWidthInECal = jwc.getPhiWidthInECal();
       j1etaWidthInHCal = jwc.getEtaWidthInHCal();
       j1phiWidthInHCal = jwc.getPhiWidthInHCal();
     }
     else if(j == 1) {
       j2PT = jet.pt();
       j2Eta = jet.eta();
       j2Phi = jet.phi();
       j2CHdFr = jet.chargedHadronEnergyFraction();
       j2NHdFr = jet.neutralHadronEnergyFraction();
       j2CEmFr = jet.chargedEmEnergyFraction();
       j2NEmFr = jet.neutralEmEnergyFraction();
       j2nCons = jet.nConstituents();
       JetWidthCalculator jwc(jet);
       j2etaWidth = jwc.getEtaWidth();
       j2phiWidth = jwc.getPhiWidth();
       j2etaWidthInECal = jwc.getEtaWidthInECal();
       j2phiWidthInECal = jwc.getPhiWidthInECal();
       j2etaWidthInHCal = jwc.getEtaWidthInHCal();
       j2phiWidthInHCal = jwc.getPhiWidthInHCal();
     }
     if(jet.pt() > 30.) { // Good jets have high pt
       HT += jet.pt();
       nGoodJets++;
     }
     //storing jetpT,jetEta,jetPhi and all other variables
     jetPt_.push_back(jet.pt());
     jetEn_.push_back(jet.energy());
     jetEta_.push_back(jet.eta());  
     jetPhi_.push_back(jet.phi());
     jetCEF_.push_back(   jet.chargedEmEnergyFraction());
     jetNEF_.push_back(   jet.neutralEmEnergyFraction());
     jetCHF_.push_back(   jet.chargedHadronEnergyFraction());
     jetNHF_.push_back(   jet.neutralHadronEnergyFraction());
     jetNCH_.push_back(   jet.chargedMultiplicity());    
     jetHFHAE_.push_back( jet.HFHadronEnergy());
     jetHFEME_.push_back( jet.HFEMEnergy());
     jetNConstituents_.push_back(jet.numberOfDaughters());
     JetWidthCalculator jwc(jet);
     jetEtaWidth_.push_back(jwc.getEtaWidth());
     jetPhiWidth_.push_back(jwc.getPhiWidth());
     jetEtaWidthInECal_.push_back(jwc.getEtaWidthInECal());
     jetPhiWidthInECal_.push_back(jwc.getPhiWidthInECal());
     jetEtaWidthInHCal_.push_back(jwc.getEtaWidthInHCal());
     jetPhiWidthInHCal_.push_back(jwc.getPhiWidthInHCal());     
     totalET += jet.pt(); // Use all jets
     nJets++;
   
     //jet PF Loose ID
     bool jetID = true;
     if (fabs(jet.eta()) <= 3.0) {
      if (!(jet.neutralHadronEnergyFraction() < 0.99))                       jetID = false;
      if (!(jet.neutralEmEnergyFraction() < 0.99))                           jetID = false;
      if (!((jet.chargedMultiplicity() + jet.neutralMultiplicity()) > 1))  jetID = false;
      if (fabs(jet.eta()) <= 2.4) {
        if (!(jet.chargedHadronEnergyFraction() > 0))  jetID = false;
        if (!(jet.chargedMultiplicity() > 0))          jetID = false;
        if (!(jet.chargedEmEnergyFraction() < 0.99))   jetID = false;
      }
    }
    if (fabs(jet.eta()) > 3.0) {
      if (!(jet.neutralEmEnergyFraction() < 0.90))  jetID = false;
      if (!(jet.neutralMultiplicity() > 10))        jetID = false;
    }
    jetPFLooseId_.push_back(jetID);

  }
   
   //ID cuts stored in simple histos
   PassVeto_ = 0;
   PassLoose_ = 0;
   PassMedium_ =0;
   PassTight_ = 0;
   PassHEEP_ = 0; 
   for(uint32_t i = 0; i < electronHandle->size(); i++) {
     //const pat::Electron &electron = (*electronHandle)[i];
     const auto el = electronHandle->ptrAt(i);
     elePt_              .push_back(el->pt());
     eleEta_             .push_back(el->eta());
     elePhi_             .push_back(el->phi());
     nEle_++;
   

     // const auto el = electronHandle->ptrAt(nEle_);
     bool isPassVeto  = (*veto_id_decisions)[el];
     if(isPassVeto) PassVeto_ = isPassVeto ? 1 : 0;

     bool isPassLoose  = (*loose_id_decisions)[el];
     if(isPassLoose) PassLoose_ = isPassLoose ? 1 : 0;

     bool isPassMedium = (*medium_id_decisions)[el];
     if(isPassMedium) PassMedium_ = isPassMedium ? 1 : 0;

     bool isPassTight  = (*tight_id_decisions)[el];
     if(isPassTight) PassTight_ = isPassTight ? 1 : 0;

     bool isPassHEEP = (*heep_id_decisions)[el];
     if(isPassHEEP) PassHEEP_ = isPassHEEP ? 1 : 0;
      
     }  
   METValue = -99, METPhi = -99, METsumEt_ = -99, METmEtSig_ = -99, METSig_ = -99;
 
   METValue = (*pfMETs)[0].pt();
   METPhi = (*pfMETs)[0].phi();
   METsumEt_ = (*pfMETs)[0].sumEt();
   METmEtSig_ = ((*pfMETs)[0].mEtSig() < 1.e10) ? (*pfMETs)[0].mEtSig() : 0;
   METSig_ = ((*pfMETs)[0].significance() < 1.e10) ? (*pfMETs)[0].significance() : 0;    
   
   std::cout << "TotalET = " << totalET << "; nJets = " << nJets 
	     << "; HT = " << HT << "; nGoodJets = " << nGoodJets 
	     << "; MET = (" << METValue << ", " << METPhi << ")" 
	     << "; nPFCands = " << pfCands->size() << std::endl;

   tree->Fill();
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetAnalyzerMC::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetAnalyzerMC::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnalyzerMC::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzerMC);
