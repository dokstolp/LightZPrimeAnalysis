#include "LightZPrimeAnalysis/JetWidthCalculator/interface/JetWidthCalculator.hh"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <math.h>
JetWidthCalculator::JetWidthCalculator(const reco::PFJet& jet) {
  jet.pt();
  jet.eta();
  jet.phi();
  //float jetet = jet.et();
  double pfCand1pt = 0;
  double pfCand2pt = 0;  
  double etSum = 0;
  double etaSum = 0;
  double etaSqSum = 0;
  double phiSum = 0;
  double phiSqSum = 0;
  double etSumInECal = 0;
  double etaSumInECal = 0;
  double etaSqSumInECal = 0;
  double phiSumInECal = 0;
  double phiSqSumInECal = 0;
  double etSumInHCal = 0;
  double etaSumInHCal = 0;
  double etaSqSumInHCal = 0;
  double phiSumInHCal = 0;
  double phiSqSumInHCal = 0;


  minEt = 9999.0;
  minPt = 9999.0;

  for(int p=0;p<9;p++){
	  NConstituents.push_back(0);
	  Nunknown.push_back(0);
	  NchargedHad.push_back(0);
	  Nelectrons.push_back(0);
	  Nmuons.push_back(0);
	  Ngamma.push_back(0);
	  NneutralHad.push_back(0);
  }

  //0,1,2,5,7,10,15,20,25
  //0,0.25,0.5,0.75,1.0,2.0,5.0,7.0,10.0

  const std::vector<reco::PFCandidatePtr> pfCands = jet.getPFConstituents();
 // std::cout<<"No. of pfCands/Constituents of a jet: "<<pfCands.size()<<std::endl;
 // double nPhotons = 0;
 // double nCHPions = 0;
 // double nMiscParticles=0;
  for(uint32_t i = 0; i < pfCands.size(); i++) {
    const edm::Ptr<reco::PFCandidate> pfCand = pfCands[i];
    if(pfCand->et()<minEt) minEt = pfCand->et();
    if(pfCand->pt()<minPt) minPt = pfCand->pt();
    NConstituents[0]++;
    if(pfCand->particleId() == 0) Nunknown[0]++;
    if(pfCand->particleId() == 1) NchargedHad[0]++;
    if(pfCand->particleId() == 2) Nelectrons[0]++;
    if(pfCand->particleId() == 3) Nmuons[0]++;
    if(pfCand->particleId() == 4) Ngamma[0]++;
    if(pfCand->particleId() == 5) NneutralHad[0]++;

    //std::cout<<i+1<<")"<<"PfCand constituent of Jet: "<< pfCand->pdgId()<<std::endl;
  /*  if (pfCand->pdgId() == 211 || pfCand->pdgId() == -211){
    nCHPions++;
    }
    else if (pfCand->pdgId() == 22){
    nPhotons++;
    }
    else{
    nMiscParticles++;
    }*/
    ptSum +=pfCand->pt();
    etSum += pfCand->et();
    //float frac = (etSum/jetet);
    //if (frac >= 0.90){
    //std::cout<<"nPFCand carrying "<< frac*100 <<" of jet energy: "<<i+1<<std::endl;
    //}
   // std::cout<<"particle is type: "<<pfCand->particleId()<<std::endl;
    etaSum += (pfCand->eta() * pfCand->et());
    etaSqSum += (pfCand->eta() * pfCand->eta() * pfCand->et());
    phiSum += (pfCand->phi() * pfCand->et());
    phiSqSum += (pfCand->phi() * pfCand->phi() * pfCand->et());
    etSumInECal += pfCand->ecalEnergy();
    etaSumInECal += (pfCand->eta() * pfCand->ecalEnergy());
    etaSqSumInECal += (pfCand->eta() * pfCand->eta() * pfCand->ecalEnergy());
    phiSumInECal += (pfCand->phi() * pfCand->ecalEnergy());
    phiSqSumInECal += (pfCand->phi() * pfCand->phi() * pfCand->ecalEnergy());
    etSumInHCal += pfCand->hcalEnergy();
    etaSumInHCal += (pfCand->eta() * pfCand->hcalEnergy());
    etaSqSumInHCal += (pfCand->eta() * pfCand->eta() * pfCand->hcalEnergy());
    phiSumInHCal += (pfCand->phi() * pfCand->hcalEnergy());
    phiSqSumInHCal += (pfCand->phi() * pfCand->phi() * pfCand->hcalEnergy());
    if(pfCand->et() <= 0.25) continue;
    NConstituents[1]++;
    if(pfCand->particleId() == 0) Nunknown[1]++;
    if(pfCand->particleId() == 1) NchargedHad[1]++;
    if(pfCand->particleId() == 2) Nelectrons[1]++;
    if(pfCand->particleId() == 3) Nmuons[1]++;
    if(pfCand->particleId() == 4) Ngamma[1]++;
    if(pfCand->particleId() == 5) NneutralHad[1]++;
    if(pfCand->et() <= 0.5) continue;
    NConstituents[2]++;
    if(pfCand->particleId() == 0) Nunknown[2]++;
    if(pfCand->particleId() == 1) NchargedHad[2]++;
    if(pfCand->particleId() == 2) Nelectrons[2]++;
    if(pfCand->particleId() == 3) Nmuons[2]++;
    if(pfCand->particleId() == 4) Ngamma[2]++;
    if(pfCand->particleId() == 5) NneutralHad[2]++;
    if(pfCand->et() <= 0.75) continue;
    NConstituents[3]++;
    if(pfCand->particleId() == 0) Nunknown[3]++;
    if(pfCand->particleId() == 1) NchargedHad[3]++;
    if(pfCand->particleId() == 2) Nelectrons[3]++;
    if(pfCand->particleId() == 3) Nmuons[3]++;
    if(pfCand->particleId() == 4) Ngamma[3]++;
    if(pfCand->particleId() == 5) NneutralHad[3]++;
    if(pfCand->et() <= 1.0) continue;
    NConstituents[4]++;
    if(pfCand->particleId() == 0) Nunknown[4]++;
    if(pfCand->particleId() == 1) NchargedHad[4]++;
    if(pfCand->particleId() == 2) Nelectrons[4]++;
    if(pfCand->particleId() == 3) Nmuons[4]++;
    if(pfCand->particleId() == 4) Ngamma[4]++;
    if(pfCand->particleId() == 5) NneutralHad[4]++;
    if(pfCand->et() <= 2.0) continue;
    NConstituents[5]++;
    if(pfCand->particleId() == 0) Nunknown[5]++;
    if(pfCand->particleId() == 1) NchargedHad[5]++;
    if(pfCand->particleId() == 2) Nelectrons[5]++;
    if(pfCand->particleId() == 3) Nmuons[5]++;
    if(pfCand->particleId() == 4) Ngamma[5]++;
    if(pfCand->particleId() == 5) NneutralHad[5]++;
    if(pfCand->et() <= 5.0) continue;
    NConstituents[6]++;
    if(pfCand->particleId() == 0) Nunknown[6]++;
    if(pfCand->particleId() == 1) NchargedHad[6]++;
    if(pfCand->particleId() == 2) Nelectrons[6]++;
    if(pfCand->particleId() == 3) Nmuons[6]++;
    if(pfCand->particleId() == 4) Ngamma[6]++;
    if(pfCand->particleId() == 5) NneutralHad[6]++;
    if(pfCand->et() <= 7.0) continue;
    NConstituents[7]++;
    if(pfCand->particleId() == 0) Nunknown[7]++;
    if(pfCand->particleId() == 1) NchargedHad[7]++;
    if(pfCand->particleId() == 2) Nelectrons[7]++;
    if(pfCand->particleId() == 3) Nmuons[7]++;
    if(pfCand->particleId() == 4) Ngamma[7]++;
    if(pfCand->particleId() == 5) NneutralHad[7]++;
    if(pfCand->et() <= 10.0) continue;
    NConstituents[8]++;
    if(pfCand->particleId() == 0) Nunknown[8]++;
    if(pfCand->particleId() == 1) NchargedHad[8]++;
    if(pfCand->particleId() == 2) Nelectrons[8]++;
    if(pfCand->particleId() == 3) Nmuons[8]++;
    if(pfCand->particleId() == 4) Ngamma[8]++;
    if(pfCand->particleId() == 5) NneutralHad[8]++;
  }
  if (pfCands.size()>1){
      pfCand1pt = pfCands[0]->pt();
      pfCand2pt = pfCands[1]->pt();
//      std::cout<<"pfCand1pt("<<pfCands[0]->pdgId()<<"): " <<pfCand1pt<<std::endl;
//      std::cout<<"pfCand2pt("<<pfCands[1]->pdgId()<<"): " <<pfCand2pt<<std::endl;
      }
  else{
      pfCand1pt = pfCands[0]->pt();
//      std::cout<<"pfCand1pt("<<pfCands[0]->pdgId()<<"): " <<pfCand1pt<<std::endl;
      }
  pfCand12PtSum = pfCand1pt + pfCand2pt;
//  std::cout<<"ptSum: " <<ptSum<<std::endl;
//  std::cout<<"pfCand12PtSum: "<<pfCand12PtSum<<std::endl;
  pt12ratio = (pfCand12PtSum/ptSum);
//  std::cout<<"pt12ratio: "<<pt12ratio<<std::endl; 
  if(etSum < 0.000001) etSum = 0.000001; // To avoid NaNs
  double etaAve = etaSum / etSum;
  double etaSqAve = etaSqSum / etSum;
  etaWidth = sqrt(etaSqAve - etaAve * etaAve);
  double phiAve = phiSum / etSum;
  double phiSqAve = phiSqSum / etSum;
  phiWidth = sqrt(phiSqAve - phiAve * phiAve);
  if(etSumInECal < 0.000001) etSumInECal = 0.000001; // To avoid NaNs
  double etaAveInECal = etaSumInECal / etSumInECal;
  double etaSqAveInECal = etaSqSumInECal / etSumInECal;
  etaWidthInECal = sqrt(etaSqAveInECal - etaAveInECal * etaAveInECal);
  double phiAveInECal = phiSumInECal / etSumInECal;
  double phiSqAveInECal = phiSqSumInECal / etSumInECal;
  phiWidthInECal = sqrt(phiSqAveInECal - phiAveInECal * phiAveInECal);
  if(etSumInECal < 0.000001) etSumInHCal = 0.000001; // To avoid NaNs
  double etaAveInHCal = etaSumInHCal / etSumInHCal;
  double etaSqAveInHCal = etaSqSumInHCal / etSumInHCal;
  etaWidthInHCal = sqrt(etaSqAveInHCal - etaAveInHCal * etaAveInHCal);
  double phiAveInHCal = phiSumInHCal / etSumInHCal;
  double phiSqAveInHCal = phiSqSumInHCal / etSumInHCal;
  phiWidthInHCal = sqrt(phiSqAveInHCal - phiAveInHCal * phiAveInHCal);
  //std::cout<<"Type of pfCands/Constituents:"<<std::endl;
  //std::cout <<"nPhotons: " <<nPhotons<<","<<" nCHPions: "<<nCHPions<<","<<" OtherMesons: "<<nMiscParticles<<std::endl;
 // std::cout<<"etaSum: "<<etaSum<<std::endl;
  //std::cout<<"etSum: "<<etSum<<std::endl;
  //std::cout<<"jet Et: "<<jetet<<std::endl;
  //std::cout<<"etaSqSum: "<<etaSqSum<<std::endl;
 // std::cout<<"phiSum: "<<phiSum<<std::endl;
 // std::cout<<"phiSqSum: "<<phiSqSum<<std::endl;
}

JetWidthCalculator::~JetWidthCalculator() {;}
