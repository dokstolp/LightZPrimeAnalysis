//For use with Ntuples made from JetAnalyzer
////Required arguments: 1 is folder containing input files, 2 is output file path, 3 is maxEvents (-1 to run over all events), 4 is reportEvery
////
////To compile using rootcom to an executable named 'analyze':
////$ ./rootcom ZprimeJetsClass analyze
////
////To run, assuming this is compiled to an executable named 'analyze':
////$ ./analyze /hdfs/store/user/uhussain/Zprime_Ntuples/ /cms/uhussain/MonoZprimeJet/CMSSW_8_0_8/src/LightZPrimeAnalysis/JetAnalyzer/test/output.root -1 10000
////Runs over every event in the folder Zprime_Ntuples, reporting progress every 10000 events
////and storing the resulting histograms in the file output.root.
////
//
#define ZprimeJetsClass_cxx
#include "ZprimeJetsClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1F.h"
#include <iostream>
#include <bitset>
#include <climits>
#include <cstring>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TStopwatch.h"
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <set>
using namespace std;
using std::vector;
int main(int argc, const char* argv[])
{ 
  Long64_t maxEvents = atof(argv[3]);
  if (maxEvents < -1LL)
  {
    std::cout<<"Please enter a valid value for maxEvents (parameter 3)."<<std::endl;
    return 1;
  }
  int reportEvery = atof(argv[4]);
  if (reportEvery < 1)
  {
    std::cout<<"Please enter a valid value for reportEvery (parameter 4)."<<std::endl;
    return 1;
  }
  //const char* file2 = argv[2];

  ZprimeJetsClass t(argv[1],argv[2]);
  t.Loop(maxEvents,reportEvery);
  return 0;
}

void ZprimeJetsClass::Loop(Long64_t maxEvents, int reportEvery)
{
  if (fChain == 0) return;
  int nTotal;
   nTotal = 0;   
   // Book histograms for recording properties of leading jet that passes dR and MET cut
  // TFile* histFile = new TFile(file2, "RECREATE");
  // h_deltar = new TH1F("j1deltaR","j1deltaR; #DeltaR of Leading Jet",50,0,0.51);h_deltar->Sumw2();
  Long64_t nentries = fChain->GetEntries();
  std::cout<<"Coming in: "<<std::endl;
  std::cout<<"nentries:"<<nentries<<std::endl;
  Long64_t nentriesToCheck = nentries;   

  std::vector<int> jetCand;
  jetCand.clear();

  //jetCandidate that passes dPhiJetMET cut out of the above jetCand
  std::vector<int> jetCand1;
  jetCand1.clear();

  float dphimin=-99;
  //Event is rejected if it contains a HighPtMuon faking MET
  bool muonVeto;

  if (maxEvents != -1LL && nentries > maxEvents)
    nentriesToCheck = maxEvents;
  nTotal = nentriesToCheck;
  Long64_t nbytes = 0, nb = 0;
  std::cout<<"Running over "<<nTotal<<" events."<<std::endl;
  for (Long64_t jentry=0; jentry<nentriesToCheck;jentry++) {
    event_.clear();
    event_info.clear();
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    jetCand   = getJetCand(200,2.4,0.8,0.1);


    muonVeto = HighPtMuonVeto(50.0); //muonPt cut of 50 GeV
    double deltar = 0.0 ;
    deltar= dR(j1etaWidth,j1phiWidth);
    fillHistos(0);
    float metcut= 0.0;
    metcut = (fabs(caloMET-pfMET))/pfMET;
    //std::cout<<"|caloMET-pfMET|/pfMET: "<<metcut<<std::endl;
    if (metFilters==0)
      {     
        fillHistos(1);
	if (HLTPFMET300==1) 
      	  {
	    fillHistos(2);
	    if(jetCand.size()>0 && j1NEmFr<0.6 && j1CHdFr > 0.4)
	      {
		fillHistos(3);
	    	if (pfMET>320)
	          {
		    fillHistos(4);
		    if (metcut<0.5)
		      {
		        fillHistos(5);
			double minDPhiJetMET = TMath::Pi();
			double minDPhiJetMET_first4 = TMath::Pi();
			for(int j = 0; j < jetCand.size(); j++)
			  {
			    if(DeltaPhi((*jetPhi)[j],pfMETPhi)<minDPhiJetMET)
			      {
				minDPhiJetMET = DeltaPhi((*jetPhi)[j],pfMETPhi);
				if(j < 4)
				  minDPhiJetMET_first4 = DeltaPhi((*jetPhi)[j],pfMETPhi);
			      }
			  }
			h_dphimin->Fill(minDPhiJetMET_first4);
			if(dPhiJetMETcut(jetCand))
			  {
			    fillHistos(6);
			    h_deltar->Fill(deltar);
			    if (deltar<0.1)
			      {
				h_deltar_cut->Fill(deltar);
				fillHistos(7);
			      }
			  }
		      }
		  }
	      }
	  }
      }
  
    tree->Fill();

    if (jentry%reportEvery == 0)
       {
       std::cout<<"Finished entry "<<jentry<<"/"<<(nentriesToCheck-1)<<std::endl;
       }
  
  }
   
   
   //save the histograms
//   histFile->Write();
//   histFile->Close();
   
}//Closing the Loop function

void ZprimeJetsClass::BookHistos(const char* file2)
{
  fileName = new TFile(file2, "RECREATE");
  tree = new TTree("ZprimeJet","ZprimeJet");
  tree->Branch("event_","std::vector<unsigned int>",&event_);
  tree->Branch("event_info","std::vector<double>",&event_info);
  fileName->cd();

  h_dphimin = new TH1F("h_dphimin","h_dphimin; Minimum dPhiJetMET",50,0,3.2);h_dphimin->Sumw2();
  h_deltar  = new TH1F("h_deltar","h_deltar; #DeltaR for Leading Jet", 10,0,0.1);h_deltar->Sumw2();
  h_deltar_cut = new TH1F("h_deltar_cut","h_deltar_cut; #DeltaR for Leading Pencil Jet",10,0,0.1);h_deltar_cut->Sumw2();
  for(int i=0; i<8; i++){

     char ptbins[100];
     sprintf(ptbins, "_%d", i);
     std::string histname(ptbins);
     h_metFilters[i] = new TH1F(("metFilters"+histname).c_str(),"metFilters",50,0,40); h_metFilters[i]->Sumw2();
     h_HLTPFMET300[i] = new TH1F(("HLTPFMET300"+histname).c_str(),"HLTJet",50,0,2);h_HLTPFMET300[i]->Sumw2();
     h_HLTPFMET170_HBHECleaned[i] = new TH1F(("HLTPFMET170"+histname).c_str(),"HLTJet",50,0,2);h_HLTPFMET170_HBHECleaned[i]->Sumw2();     
     h_HT[i] = new TH1F(("HT"+histname).c_str(), "HT",100,0,2500);h_HT[i]->Sumw2();
     h_nJets[i]   = new TH1F(("nJets"+histname).c_str(), "nJets;Number of Jets", 130, 0, 130);h_nJets[i]->Sumw2();
     h_nGoodJets[i]   = new TH1F(("nGoodJets"+histname).c_str(), "nGoodJets;Number of Jets that pass cut(30 GeV)", 10, 0, 10);h_nGoodJets[i]->Sumw2();
     h_pfMET[i] = new TH1F(("pfMET"+histname).c_str(), "pfMET",50,280,1000);h_pfMET[i] ->Sumw2();
     h_pfMETPhi[i] = new TH1F(("pfMETPhi"+histname).c_str(), "pfMETPhi",50,-4,4);h_pfMETPhi[i]->Sumw2();
     h_j1Pt[i]  = new TH1F(("j1pT"+histname).c_str(), "j1pT;p_{T} of Leading Jet [GeV]", 50, 100, 1000);h_j1Pt[i]->Sumw2();
     h_j1Eta[i] = new TH1F(("j1Eta"+histname).c_str(), "j1Eta; #eta of Leading Jet", 100, -3.0, 3.0);h_j1Eta[i]->Sumw2();
     h_j1Phi[i] = new TH1F(("j1Phi"+histname).c_str(), "j1Phi; #phi of Leading Jet", 100, -4, 4);h_j1Phi[i]->Sumw2();
     h_j1etaWidth[i] = new TH1F(("j1etaWidth"+histname).c_str(),"j1etaWidh; #eta width of Leading Jet", 50,0,0.1);h_j1etaWidth[i] ->Sumw2();
     h_j1phiWidth[i] = new TH1F(("j1phiWidth"+histname).c_str(),"j1phiWidth; #phi width of Leading Jet", 50, 0,0.1);h_j1phiWidth[i]->Sumw2();
     h_j1etaWidthInECal[i] = new TH1F(("j1etaWidthInECal"+histname).c_str(),"j1etaWidthInECal; #eta width in ECAL of Leading Jet", 40,0,0.2); h_j1etaWidthInECal[i]->Sumw2();
     h_j1etaWidthInHCal[i] = new TH1F(("j1etaWidthInHCal"+histname).c_str(),"j1etaWidthInHCal; #eta width in HCAL of Leading Jet", 40,0,0.2);h_j1etaWidthInHCal[i]->Sumw2();
     h_j1phiWidthInECal[i] = new TH1F(("j1phiWidthInECal"+histname).c_str(),"j1phiWidthInECal; #phi width in ECAL of Leading Jet", 40,0,0.2);h_j1phiWidthInECal[i]->Sumw2(); 
     h_j1phiWidthInHCal[i] = new TH1F(("j1phiWidthInHCal"+histname).c_str(),"j1phiWidthInHCal; #phi width in HCAL of Leading Jet", 40,0,0.2);h_j1phiWidthInHCal[i]->Sumw2(); 
     h_j1nCons[i] = new TH1F (("j1nCons"+histname).c_str(),"j1NConstituents; Number of Constituents of Leading Jet",100, 0, 100);h_j1nCons[i]->Sumw2();
     h_j1CEF[i] = new TH1F(("j1CEF"+histname).c_str(), "j1CEF; Charged Electromagnetic Energy Fraction of leading Jet",50, 0, 1.0);h_j1CEF[i]->Sumw2();
     h_j1NEF[i] = new TH1F(("j1NEF"+histname).c_str(), "j1NEF; Neutral Electromagnetic Energy Fraction of leading Jet",30, 0, 0.6);h_j1NEF[i]->Sumw2();
     h_j1CHF[i] = new TH1F(("j1CHF"+histname).c_str(), "j1CHF; Charged Hadron Energy Fraction of leading Jet",30, 0.4, 1.0);h_j1CHF[i]->Sumw2();
     h_j1NHF[i] = new TH1F(("j1NHF"+histname).c_str(), "j1NHF; Neutral Hadron Energy Fraction of leading Jet",40, 0.0, 0.8);h_j1NHF[i]->Sumw2();
   }
}

double ZprimeJetsClass::dR(double jetetaWidth, double jetphiWidth)
{
  double deltar = sqrt(jetetaWidth*jetetaWidth + jetphiWidth*jetphiWidth);
  return deltar;
}


void ZprimeJetsClass::fillHistos(int histoNumber)
{
  h_metFilters[histoNumber]->Fill(metFilters);
  h_HLTPFMET170_HBHECleaned[histoNumber]->Fill(HLTPFMET170_HBHECleaned);
  h_HLTPFMET300[histoNumber]->Fill(HLTPFMET300);
  h_HT[histoNumber]->Fill(HT);
  h_nJets[histoNumber]->Fill(nJets);
  h_nGoodJets[histoNumber]->Fill(nGoodJets);
  h_pfMET[histoNumber]->Fill(pfMET);
  h_pfMETPhi[histoNumber]->Fill(pfMETPhi);
  h_j1Pt[histoNumber]->Fill(j1PT);
  h_j1Eta[histoNumber]->Fill(j1Eta);
  h_j1Phi[histoNumber]->Fill(j1Phi);
  h_j1etaWidth[histoNumber]->Fill(j1etaWidth);
  h_j1phiWidth[histoNumber]->Fill(j1phiWidth);
  h_j1etaWidthInECal[histoNumber]->Fill(j1etaWidthInECal);
  h_j1etaWidthInHCal[histoNumber]->Fill(j1etaWidthInHCal);
  h_j1phiWidthInECal[histoNumber]->Fill(j1phiWidthInECal);
  h_j1phiWidthInHCal[histoNumber]->Fill(j1phiWidthInHCal);
  h_j1nCons[histoNumber]->Fill(j1nCons);
  h_j1CEF[histoNumber]->Fill(j1CEmFr);
  h_j1NEF[histoNumber]->Fill(j1NEmFr);
  h_j1CHF[histoNumber]->Fill(j1CHdFr);
  h_j1NHF[histoNumber]->Fill(j1NHdFr);

  
}

//Gives the (minimum) separation in phi between the specified phi values
//Must return a positive value
float ZprimeJetsClass::DeltaPhi(float phi1, float phi2)
{
  float pi = TMath::Pi();
  float dphi = fabs(phi1-phi2);
  if(dphi>pi)
    dphi = 2.0*pi - dphi;
  return dphi;
}
float ZprimeJetsClass::dPhiJetMETmin(std::vector<int> jets)
{
  float dPhimin=TMath::Pi();
  int njetsMax = jets.size();
  if(njetsMax > 4)
    njetsMax = 4;
  for(int j=0;j< njetsMax; j++)
    {
      float dPhi = DeltaPhi((*jetPhi)[j],pfMETPhi);
      std::cout<<"DeltaPhi: "<<dPhi<<std::endl;
      if(dPhi < dPhimin){
        dPhimin = dPhi;
     }
   }
  return dPhimin;
}
bool ZprimeJetsClass::HighPtMuonVeto(double muonPtCut)
{
  //pass veto if no muon in the event is found with Pt > muonPtCut 
  bool passes = true;
  for (int i=0;i<nMu;i++)
    {
      if((*muPt)[i] > muonPtCut){
        passes = false;
    }
   }
  return passes;
}

std::vector<int> ZprimeJetsClass::getJetCand(double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut){

  std::vector<int> tmpCand;
  tmpCand.clear();

  for(int p=0;p<nJets;p++)
    {

      bool kinematic = (*jetPt)[p] > jetPtCut && (*jetNHF)[p] < jetNHFCut && (*jetCHF)[p] > jetCHFCut && fabs((*jetEta)[p])<jetEtaCut;

      if((*jetPFLooseId)[p]==1 && kinematic){
        tmpCand.push_back(p);
      }
    }

  return tmpCand;
}

bool ZprimeJetsClass::dPhiJetMETcut(std::vector<int> jets)
{
  //reject jet if it is found within DeltaPhi(jet,MET) < 0.5 
  bool passes = false;
  
  int njetsMax = jets.size();
  //Only look at first four jets (because that's what monojet analysis do)
  if(njetsMax > 4)
    njetsMax = 4;
  int j=0;
  for(;j< njetsMax; j++){

    if(DeltaPhi((*jetPhi)[j],pfMETPhi) < 0.5)
      break;
  }

  if(j==njetsMax)
    passes = true;

  return passes;
  
} 



//bool ZprimeJetsClass::electron_veto(int jet_index)
//{
//  bool veto_passed = true;
//  for(int i = 0; i < nEle; i++)
//    {
//      if (if(elePt->at(i)> 10 && eleEta->at(i) < 2.4 && )
//    }
//
//
//}
