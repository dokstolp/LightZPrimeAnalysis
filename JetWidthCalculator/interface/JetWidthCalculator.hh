#ifndef JetWidthCalculator_hh
#define JetWidthCalculator_hh

// This is a simple helper class to compute energy weighted eta and phi 
// widths of a jet.  Width in ECal and HCal are available 

#include "DataFormats/PatCandidates/interface/Jet.h"
#include <vector>
using namespace std;
class JetWidthCalculator {

public:

  JetWidthCalculator(const reco::PFJet&);

  virtual ~JetWidthCalculator();

  double getEtaWidth() {return etaWidth;}
  double getEtaWidthInECal() {return etaWidthInECal;}
  double getEtaWidthInHCal() {return etaWidthInHCal;}

  double getPhiWidth() {return phiWidth;}
  double getPhiWidthInECal() {return phiWidthInECal;}
  double getPhiWidthInHCal() {return phiWidthInHCal;}

  double getPFCand12PtSum() {return pfCand12PtSum;}
  double getPFCandsPtSum() {return ptSum;}
  double getPFCand12Ratio() {return pt12ratio;}

  double getMinPt() {return minPt;}
  double getMinEt() {return minEt;}
  std::vector<int> getNConstituents() {return NConstituents;}
  std::vector<int> getNunknown() {return Nunknown;}
  std::vector<int> getNchargedHad() {return NchargedHad;}
  std::vector<int> getNelectrons() {return Nelectrons;}
  std::vector<int> getNmuons() {return Nmuons;}
  std::vector<int> getNgamma() {return Ngamma;}
  std::vector<int> getNneutralHad() {return NneutralHad;}
   
private:

  // No default constructor is possible

  JetWidthCalculator();

  // No copy constructor is needed

  JetWidthCalculator(const JetWidthCalculator&);

  // No equality operator is needed

  const JetWidthCalculator& operator=(const JetWidthCalculator&);

  double etaWidth;
  double etaWidthInECal;
  double etaWidthInHCal;

  double phiWidth;
  double phiWidthInECal;
  double phiWidthInHCal;
  
  double pfCand12PtSum;
  double ptSum;
  double pt12ratio;

  double minPt;
  double minEt;

  std::vector<int> NConstituents;
  std::vector<int> Nunknown;
  std::vector<int> NchargedHad;
  std::vector<int> Nelectrons;
  std::vector<int> Nmuons;
  std::vector<int> Ngamma;
  std::vector<int> NneutralHad;
};

#endif
