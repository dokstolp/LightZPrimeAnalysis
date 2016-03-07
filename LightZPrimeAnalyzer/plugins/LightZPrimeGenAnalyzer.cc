// -*- C++ -*-
//
// Package:    LightZPrimeAnalysis/LightZPrimeGenAnalyzer
// Class:      LightZPrimeGenAnalyzer
// 
/**\class LightZPrimeGenAnalyzer LightZPrimeGenAnalyzer.cc LightZPrimeAnalysis/LightZPrimeGenAnalyzer/plugins/LightZPrimeGenAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Sridhara Rao Dasu
//         Created:  Tue, 23 Feb 2016 04:57:10 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <list>
#include <cmath>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "TLorentzVector.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class LightZPrimeGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit LightZPrimeGenAnalyzer(const edm::ParameterSet&);
  ~LightZPrimeGenAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------

  edm::EDGetTokenT< vector<reco::GenJet> > genJetsToken;
  edm::EDGetTokenT< vector<reco::GenMET> > genMETsToken;
  edm::EDGetTokenT< vector<reco::GenParticle> > genParticlesToken;

  double totalET;
  double HT;
  double METValue;
  double METPhi;
  double j1PT;
  double j1Eta;
  double j1Phi;
  double j1EMFr;
  double j1CHdFr;
  double j2PT;
  double j2Eta;
  double j2Phi;
  double j2EMFr;
  double j2CHdFr;
  double j12PT;
  double j12Eta;
  double j12Phi;
  double j12Mass;
  double dRzPRemnants;

  uint32_t nGenPs;
  uint32_t nJets;
  uint32_t j1nCons;
  uint32_t j2nCons;

  TLorentzVector MET;
  TLorentzVector goodPion;
  TLorentzVector goodOSPartner;
  TLorentzVector goodSSPartner;
  TLorentzVector goodOSPair;
  TLorentzVector goodSSPair;

  TTree* tree;

};

//
// constants, enums and typedefs
//

class GoodGenPion {
public:
  const reco::GenParticle &candidate;
  double neutralEMIsolation;
  double neutralHDIsolation;
  double chargedHDIsolation;
  int nNeutralEMNeighbors;
  int nNeutralHDNeighbors;
  int nChargedHDNeighbors;
  GoodGenPion(const reco::GenParticle c, double nEI, double nHI, double cHI, int nNE, int nNH, int nCH) :
    candidate(c), neutralEMIsolation(nEI), neutralHDIsolation(nHI), chargedHDIsolation(cHI),
    nNeutralEMNeighbors(nNE), nNeutralHDNeighbors(nNH), nChargedHDNeighbors(nCH) {;}
};

// comparison, based on pt of the candidates
bool compareGPCands(const GoodGenPion& first, const GoodGenPion& second)
{
  return ( first.candidate.pt() < second.candidate.pt() );
}

double getDeltaR(double eta1, double phi1, double eta2, double phi2) {
  double dR = sqrt((eta2 - eta1) * (eta2 - eta1) + (phi2 - phi1) * (phi2 - phi1));
  return dR;
}

//
// static data member definitions
//

//
// constructors and destructor
//
LightZPrimeGenAnalyzer::LightZPrimeGenAnalyzer(const edm::ParameterSet& iConfig)
{
  genParticlesToken = consumes<vector<reco::GenParticle> > (edm::InputTag("genParticles"));
  genJetsToken = consumes< vector<reco::GenJet> >(edm::InputTag("ak4GenJetsNoNu"));
  genMETsToken = consumes< vector<reco::GenMET> >(edm::InputTag("genMetTrue"));
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("GenEventTree", "Generator level event data");

  tree->Branch("totalET", &totalET);
  tree->Branch("HT", &HT);
  tree->Branch("METValue", &METValue);
  tree->Branch("METPhi", &METPhi);
  tree->Branch("j1PT", &j1PT);
  tree->Branch("j1Eta", &j1Eta);
  tree->Branch("j1Phi", &j1Phi);
  tree->Branch("j1EMFr", &j1EMFr);
  tree->Branch("j1CHdFr", &j1CHdFr);
  tree->Branch("j2PT", &j2PT);
  tree->Branch("j2Eta", &j2Eta);
  tree->Branch("j2Phi", &j2Phi);
  tree->Branch("j2EMFr", &j2EMFr);
  tree->Branch("j2CHdFr", &j2CHdFr);
  tree->Branch("j12PT", &j12PT);
  tree->Branch("j12Eta", &j12Eta);
  tree->Branch("j12Phi", &j12Phi);
  tree->Branch("j12Mass", &j12Mass);
  tree->Branch("dRzPRemnants", &dRzPRemnants);

  tree->Branch("nGenPs", &nGenPs);
  tree->Branch("nJets", &nJets);
  tree->Branch("j1nCons", &j1nCons);
  tree->Branch("j2nCons", &j2nCons);

}

LightZPrimeGenAnalyzer::~LightZPrimeGenAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LightZPrimeGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  Handle< vector<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken, genParticles);

  Handle< vector<reco::GenJet> > genJets;
  iEvent.getByToken(genJetsToken, genJets);

  Handle< vector<reco::GenMET> > genMETs;
  iEvent.getByToken(genMETsToken, genMETs);

  // Get the generated ZPrime

  nGenPs = genParticles->size();
  TLorentzVector zPrime;
  for(uint32_t p = 0; p < genParticles->size(); p++) {
    const reco::GenParticle &particle = (*genParticles)[p];
    if(particle.pdgId() == 600001) {
      zPrime.SetPtEtaPhiE(particle.pt(), particle.eta(), particle.phi(), particle.energy());
      break;
    }
  }

  // Initialize defaults

  totalET = 0;
  HT = 0;
  nJets = 0;
  j1PT = j1Eta = j1Phi = j1EMFr = j1CHdFr = 0;
  j1nCons = 0;
  j2PT = j2Eta = j2Phi = j2EMFr = j2CHdFr = 0;
  j2nCons = 0;
  j12PT = j12Eta = j12Phi = 0;

  // Set event level quantities

  MET.SetPtEtaPhiE((*genMETs)[0].pt(), 0., (*genMETs)[0].phi(), 0.);
  METValue = MET.Pt();
  METPhi = MET.Phi();

  // Compute HT
   
  double dRzPj1 = -999.;
  double dRzPj2 = -999.;
  double dRzPj12 = -999.;
  TLorentzVector j12;
  for(uint32_t j = 0; j < genJets->size(); j++) {
    const reco::GenJet &jet = (*genJets)[j];
    if(jet.pt() > 30.) {
      HT += jet.pt();
      nJets++;
    }
    totalET += jet.pt();
    if(j == 0) {
      j1PT = jet.pt();
      j1Eta = jet.eta();
      j1Phi = jet.phi();
      j1EMFr = jet.emEnergy() / jet.energy();
      j1CHdFr = jet.hadEnergy() / jet.energy();
      j1nCons = jet.nConstituents();
      dRzPj1 = getDeltaR(zPrime.Eta(), zPrime.Phi(), jet.eta(), jet.phi());
      dRzPRemnants = dRzPj1;
    }
    if(j == 1) {
      j2PT = jet.pt();
      j2Eta = jet.eta();
      j2Phi = jet.phi();
      j2EMFr = jet.emEnergy() / jet.energy();
      j2CHdFr = jet.hadEnergy() / jet.energy();
      j2nCons = jet.nConstituents();
      dRzPj2 = getDeltaR(zPrime.Eta(), zPrime.Phi(), jet.eta(), jet.phi());
      dRzPRemnants = min(dRzPj1, dRzPj2);
      j12 = 
	TLorentzVector((*genJets)[0].px(), (*genJets)[0].py(), (*genJets)[0].pz(), (*genJets)[0].energy()) +
	TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy());
      j12PT = j12.Pt();
      j12Eta = j12.Eta();
      j12Phi = j12.Phi();
      j12Mass = j12.M();
      dRzPj12 = getDeltaR(zPrime.Eta(), zPrime.Phi(), j12.Eta(), j12.Phi());
      dRzPRemnants = min(dRzPRemnants, dRzPj12);
    }
  }
  
  // Fill the TTree
  tree->Fill();
  
}

// ------------ method called once each job just before starting event loop  ------------
void 
LightZPrimeGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LightZPrimeGenAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LightZPrimeGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LightZPrimeGenAnalyzer);
