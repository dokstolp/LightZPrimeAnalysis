#include "LightZPrimeAnalysis/JetWidthCalculator/interface/JetWidthCalculator.hh"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include <math.h>

JetWidthCalculator::JetWidthCalculator(const pat::Jet& jet) {
  jet.pt();
  jet.eta();
  jet.phi();
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
  const std::vector<reco::PFCandidatePtr> pfCands = jet.getPFConstituents();
  for(uint32_t i = 0; i < pfCands.size(); i++) {
    const edm::Ptr<reco::PFCandidate> pfCand = pfCands[i];
    etSum += pfCand->et();
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
  }
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
}

JetWidthCalculator::~JetWidthCalculator() {;}
