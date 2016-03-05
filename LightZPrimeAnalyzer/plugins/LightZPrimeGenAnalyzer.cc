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
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "TLorentzVector.h"

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

  uint32_t nJets;

  TLorentzVector MET;
  TLorentzVector goodPion;
  TLorentzVector goodOSPartner;
  TLorentzVector goodSSPartner;
  TLorentzVector goodOSPair;
  TLorentzVector goodSSPair;

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

//
// static data member definitions
//

//
// constructors and destructor
//
LightZPrimeGenAnalyzer::LightZPrimeGenAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   genParticlesToken = consumes<vector<reco::GenParticle> > (edm::InputTag("genParticles"));
   genJetsToken = consumes< vector<reco::GenJet> >(edm::InputTag("ak4GenJetsNoNu"));
   genMETsToken = consumes< vector<reco::GenMET> >(edm::InputTag("genMetTrue"));
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

   // Header

   static bool first = true;
   if(first) {
     first = false;
     std::cout 
       << "  nGenPs TotalET nJets30      HT     MET  METPhi    j1PT   j1Eta   j1Phi j1nCons  j1EMFr j1CHdFr    j2PT   j2Eta   j2Phi j2nCons  j2EMFr j2CHdFr"
       << std::endl;
   }

   // Set event level quantities

   totalET = 0; // Computed in pfCands loop
   HT = 0; // Computed in pfJets loop
   nJets = 0; // Computed in pfJets loop
   MET.SetPtEtaPhiE((*genMETs)[0].pt(), 0., (*genMETs)[0].phi(), 0.);

   // Compute HT
   
   for(uint32_t j = 0; j < genJets->size(); j++) {
     const reco::GenJet &jet = (*genJets)[j];
     if(jet.pt() > 30.) {
       HT += jet.pt();
       nJets++;
     }
     totalET += jet.pt();
   }

   std::cout.precision(2);

   std::cout
     << std::fixed
     << setw(8)
     << genParticles->size()
     << setw(8)
     << totalET
     << setw(8)
     << nJets
     << setw(8)
     << HT
     << setw(8)
     << MET.Pt()
     << setw(8)
     << MET.Phi();

   if(genJets->size() > 0)
     std::cout
       << std::fixed
       << setw(8)
       << (*genJets)[0].pt()
       << setw(8)
       << (*genJets)[0].eta()
       << setw(8)
       << (*genJets)[0].phi()
       << setw(8)
       << (*genJets)[0].nConstituents()
       << setw(8)
       << (*genJets)[0].emEnergy() / (*genJets)[0].energy()
       << setw(8)
       << (*genJets)[0].hadEnergy() / (*genJets)[0].energy();
   else
     std::cout
       << std::fixed
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0;
   
   if(genJets->size() > 1)
     std::cout
       << std::fixed
       << setw(8)
       << (*genJets)[1].pt()
       << setw(8)
       << (*genJets)[1].eta()
       << setw(8)
       << (*genJets)[1].phi()
       << setw(8)
       << (*genJets)[1].nConstituents()
       << setw(8)
       << (*genJets)[1].emEnergy() / (*genJets)[1].energy()
       << setw(8)
       << (*genJets)[1].hadEnergy() / (*genJets)[1].energy();
   else
     std::cout
       << std::fixed
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0
       << setw(8)
       << 0;
   
   std::cout << std::endl;

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
