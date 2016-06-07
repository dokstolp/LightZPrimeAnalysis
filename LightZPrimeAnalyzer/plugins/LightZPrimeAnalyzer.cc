// -*- C++ -*-
//
// Package:    LightZPrimeAnalysis/LightZPrimeAnalyzer
// Class:      LightZPrimeAnalyzer
// 
/**\class LightZPrimeAnalyzer LightZPrimeAnalyzer.cc LightZPrimeAnalysis/LightZPrimeAnalyzer/plugins/LightZPrimeAnalyzer.cc

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

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "LightZPrimeAnalysis/JetWidthCalculator/interface/JetWidthCalculator.hh"

#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class LightZPrimeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit LightZPrimeAnalyzer(const edm::ParameterSet&);
      ~LightZPrimeAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  edm::EDGetTokenT< vector<pat::PackedCandidate> > packedPFCandsToken;
  edm::EDGetTokenT< vector<pat::Jet> > slimmedJetsToken;
  edm::EDGetTokenT< vector<pat::MET> > slimmedMETsToken;

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

class GoodPionCandidate {
public:
  const pat::PackedCandidate &candidate;
  double neutralEMIsolation;
  double neutralHDIsolation;
  double chargedHDIsolation;
  int nNeutralEMNeighbors;
  int nNeutralHDNeighbors;
  int nChargedHDNeighbors;
  GoodPionCandidate(const pat::PackedCandidate c, double nEI, double nHI, double cHI, int nNE, int nNH, int nCH) :
    candidate(c), neutralEMIsolation(nEI), neutralHDIsolation(nHI), chargedHDIsolation(cHI),
    nNeutralEMNeighbors(nNE), nNeutralHDNeighbors(nNH), nChargedHDNeighbors(nCH) {;}
};

// comparison, based on pt of the candidates
bool compareGPCands(const GoodPionCandidate& first, const GoodPionCandidate& second)
{
  return ( first.candidate.pt() < second.candidate.pt() );
}

//
// static data member definitions
//

//
// constructors and destructor
//
LightZPrimeAnalyzer::LightZPrimeAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   packedPFCandsToken = consumes< vector<pat::PackedCandidate> >(edm::InputTag("packedPFCandidates"));
   slimmedJetsToken = consumes< vector<pat::Jet> >(edm::InputTag("slimmedJets"));
   slimmedMETsToken = consumes< vector<pat::MET> >(edm::InputTag("slimmedMETs"));
}


LightZPrimeAnalyzer::~LightZPrimeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LightZPrimeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle< vector<pat::PackedCandidate> > pfCands;
   iEvent.getByToken(packedPFCandsToken, pfCands);

   Handle< vector<pat::Jet> > pfJets;
   iEvent.getByToken(slimmedJetsToken, pfJets);

   Handle< vector<pat::MET> > pfMETs;
   iEvent.getByToken(slimmedMETsToken, pfMETs);

   // Set event level quantities

   totalET = 0; // Computed in pfCands loop
   HT = 0; // Computed in pfJets loop
   nJets = 0; // Computed in pfJets loop
   MET.SetPtEtaPhiE((*pfMETs)[0].pt(), 0., (*pfMETs)[0].phi(), 0.);

   // Compute HT
   
   for(uint32_t j = 0; j < pfJets->size(); j++) {
     const pat::Jet &jet = (*pfJets)[j];
     JetWidthCalculator jwc(jet);
     std::cout << jwc.getEtaWidth() << std::endl;
     if(jet.pt() > 30.) {
       HT += jet.pt();
       nJets++;
     }
   }
   std::cout << "HT = " << HT << "; nJets = " << nJets 
	     << "; MET = (" << MET.Pt() << ", " << MET.Phi() << ")" 
	     << "; nPFCands = " << pfCands->size() << std::endl;

   // Select good pions, meaning those which are isolated except for another 
   // "higher PT" charged pion in the immediate neighborhood

   list<GoodPionCandidate> goodPions;
   for(uint32_t i = 0; i < pfCands->size(); i++) {
     // Select only charged hadrons
     const pat::PackedCandidate &iCand = (*pfCands)[i];
     totalET += iCand.pt();
     if(abs(iCand.pdgId()) == 211) {
       // Find good charged hadrons
       double neutralEMIsolation = 0;
       double neutralHDIsolation = 0;
       double chargedHDIsolation = 0;
       int nNeutralEMNeighbors = 0;
       int nNeutralHDNeighbors = 0;
       int nChargedHDNeighbors = 0;
       double highestChargedNeighborPT = 0;
       for(uint32_t j = 0; j < pfCands->size(); j++) {
	 if(i != j) {	   
	   const pat::PackedCandidate &jCand = (*pfCands)[j];
	   // FIXME: Make sure this is from the same primary vertex as original
	   double ijDeltaR = deltaR(iCand.eta(), iCand.phi(), jCand.eta(), jCand.phi());
	   if(ijDeltaR < 0.5) {
	     if(jCand.pdgId() == 22) {
	       neutralEMIsolation += jCand.pt();
	       nNeutralEMNeighbors++;
	     }
	     else if(jCand.pdgId() == 130) {
	       neutralHDIsolation += jCand.pt();
	       nNeutralHDNeighbors++;
	     }
	     else if(abs(jCand.pdgId()) == 211) {
	       if(jCand.pt() > highestChargedNeighborPT) highestChargedNeighborPT = jCand.pt();
	       chargedHDIsolation += jCand.pt();
	       nChargedHDNeighbors++;
	     }
	     else {
	       //	       std::cout << "Unexpected pdgId() ignored; Id = " << jCand.pdgId() << std::endl; 
	     }
	   }
	   // To speed things up a bit
	   double combinedIsolation = neutralEMIsolation + neutralHDIsolation;
	   if(nChargedHDNeighbors > 1) {
	     combinedIsolation += (neutralEMIsolation - highestChargedNeighborPT);
	   }
	   if(combinedIsolation > iCand.pt() * 0.50) break;
	 }
       }
       neutralEMIsolation = neutralEMIsolation / iCand.pt();
       neutralHDIsolation = neutralHDIsolation / iCand.pt();
       chargedHDIsolation = chargedHDIsolation / iCand.pt();
       if(neutralEMIsolation < 0.10 &&
	  neutralHDIsolation < 0.10 &&
	  chargedHDIsolation < 0.10) {
	 goodPions.push_back(GoodPionCandidate(iCand, neutralEMIsolation, neutralHDIsolation, chargedHDIsolation, 
					       nNeutralEMNeighbors, nNeutralHDNeighbors, nChargedHDNeighbors));
       }
     }
   }
   goodPion = TLorentzVector();
   goodOSPartner = TLorentzVector();
   goodOSPair = TLorentzVector();
   goodSSPartner = TLorentzVector();
   goodSSPair = TLorentzVector();
   std::cout << "totalET = " << totalET 
	     << "; nGoodPions = " << goodPions.size()
	     << std::endl;
   if(goodPions.size() > 0) {
     goodPions.sort(compareGPCands);
     // Select highest PT good charged hadron as goodPion
     goodPion.SetPtEtaPhiE(goodPions.begin()->candidate.pt(), 
			   goodPions.begin()->candidate.eta(), 
			   goodPions.begin()->candidate.phi(), 
			   goodPions.begin()->candidate.energy());
     if(goodPions.size() > 1) {
       // Select next highest oppositely charged hadron as goodOSPartner, which should be available in the signal case
       // Select next highest identically charged hadron as goodSSPartner, which is background
       // Form OS and SS pair candidates
       // The SS pair can serve as a combinatorial background estimate
       bool osPairFound = false;
       bool ssPairFound = false;
       std::list<GoodPionCandidate>::const_iterator iterator = goodPions.begin();
       for(iterator++; iterator != goodPions.end(); iterator++) {
	 GoodPionCandidate partner = *iterator;
	 if(partner.candidate.pdgId() == (-goodPions.begin()->candidate.pdgId())) {
	   goodOSPartner.SetPtEtaPhiE(partner.candidate.pt(), partner.candidate.eta(), partner.candidate.phi(), partner.candidate.energy());
	   // Form good opposite sign pair candidate
	   goodOSPair = goodPion + goodOSPartner;
	   osPairFound = true;
	 }
	 else if(partner.candidate.pdgId() == goodPions.begin()->candidate.pdgId()) {
	   goodSSPartner.SetPtEtaPhiE(partner.candidate.pt(), partner.candidate.eta(), partner.candidate.phi(), partner.candidate.energy());
	   // Form good same sign pair candidate
	   goodSSPair = goodPion + goodSSPartner;
	   ssPairFound = true;
	 }
	 if(osPairFound && ssPairFound) break;
       }
     }
   }
   std::cout << "goodPion: (Pt, Eta, Phi) = (" 
	     << goodPion.Pt() << ", " << goodPion.Eta() << ", " << goodPion.Phi() << ")" << std::endl;
   std::cout << "goodOSPartner: (Pt, Eta, Phi) = (" 
	     << goodOSPartner.Pt() << ", " << goodOSPartner.Eta() << ", " << goodOSPartner.Phi() << ")" << std::endl;
   std::cout << "goodOSPair: (M, Pt) = (" << goodOSPair.M() << ", " << goodOSPair.Pt() << std::endl;
   std::cout << "goodSSPartner: (Pt, Eta, Phi) = (" 
	     << goodSSPartner.Pt() << ", " << goodSSPartner.Eta() << ", " << goodSSPartner.Phi() << ")" << std::endl;
   std::cout << "goodSSPair: (M, Pt) = (" << goodOSPair.M() << ", " << goodOSPair.Pt() << std::endl;
}

// ------------ method called once each job just before starting event loop  ------------
void 
LightZPrimeAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
LightZPrimeAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LightZPrimeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LightZPrimeAnalyzer);
