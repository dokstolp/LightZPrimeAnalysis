//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 17 06:52:45 2016 by ROOT version 6.06/01
// from TTree JetTree/Jet data for analysis
// found on file: /hdfs/store/user/uhussain/Zprime_Ntuples/MET/crab_JetAnalyzer_MET_2016B/160816_111908/0000/JetAnalyzer_487.root
//////////////////////////////////////////////////////////

#ifndef ZprimeJetsClass_h
#define ZprimeJetsClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <list>
#include <vector>
#include <bitset>
#include <TCanvas.h>
#include <TSystem.h>
#include <TPostScript.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TRef.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TDCacheFile.h>
#include <TLorentzVector.h>
#include "TIterator.h"
// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
using namespace std;
class ZprimeJetsClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::vector<unsigned int> event_;
   std::vector<double> event_info;
   TFile *fileName;
   TTree *tree;

   TH1F *h_dphimin,*h_deltar_cut,*h_deltar, *h_metFilters[9], *h_HLTPFMET300[9], *h_HLTPFMET170_HBHECleaned[9], *h_HT[9], *h_nJets[9], *h_nGoodJets[9], *h_pfMET[9], *h_pfMETPhi[9], *h_j1Pt[9], *h_j1Eta[9], *h_j1Phi[9], *h_j1etaWidth[9], *h_j1phiWidth[9], *h_j1etaWidthInECal[9], *h_j1etaWidthInHCal[9], *h_j1phiWidthInECal[9], *h_j1phiWidthInHCal[9], *h_j1nCons[9], *h_j1CEF[9], *h_j1NEF[9], *h_j1CHF[9], *h_j1NHF[9];

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
                
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   UInt_t          metFilters;
   Int_t           npv;
   Int_t           nTrksPV;
   Float_t         rho;
   Double_t        totalET;
   Double_t        HT;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         caloMET;
   UInt_t          nJets;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetEn;
   vector<float>   *jetPhi;
   vector<bool>    *jetPFLooseId;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<int>     *jetNCH;
   vector<float>   *jetHFHAE;
   vector<float>   *jetHFEME;
   vector<int>     *jetNConstituents;
   vector<float>   *jetEtaWidth;
   vector<float>   *jetPhiWidth;
   vector<float>   *jetEtaWidthInECal;
   vector<float>   *jetEtaWidthInHCal;
   vector<float>   *jetPhiWidthInECal;
   vector<float>   *jetPhiWidthInHCal;
   Int_t           nEle;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<int>     *ElectronPassVetoID;
   vector<int>     *ElectronPassLooseID;
   vector<int>     *ElectronPassMediumID;
   vector<int>     *ElectronPassTightID;
   vector<int>     *ElectronPassHEEPID;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<bool>    *muIsLooseID;
   vector<bool>    *muIsMediumID;
   vector<bool>    *muIsTightID;
   vector<bool>    *muIsSoftID;
   vector<bool>    *muIsHighPtID;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaAtCalo;
   vector<float>   *eleSigmaIEtaIEta;
   vector<float>   *eleSigmaIPhiIPhi;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *elePFMiniIso;
   vector<float>   *eleIDMVANonTrg;
   vector<float>   *eleIDMVATrg;
   vector<float>   *eledEtaseedAtVtx;
   vector<float>   *eleE1x5;
   vector<float>   *eleE2x5;
   vector<float>   *eleE5x5;
   vector<float>   *eleE1x5Full5x5;
   vector<float>   *eleE2x5Full5x5;
   vector<float>   *eleE5x5Full5x5;
   vector<float>   *eleR9Full5x5;
   vector<float>   *eleisoChargedHadrons;
   vector<float>   *eleisoNeutralHadrons;
   vector<float>   *eleisoPhotons;
   vector<float>   *eleisoChargedFromPU;
   Int_t           nPhotons;
   vector<float>   *phopt;
   vector<float>   *phoeta;
   vector<float>   *phophi;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phohOverE;
   vector<int>     *phohasPixelSeed;
   vector<float>   *phoisoChargedHadrons;
   vector<float>   *phoWorstisoChargedHadrons;
   vector<float>   *phoisoNeutralHadrons;
   vector<float>   *phoisoPhotons;
   vector<float>   *trackPt;
   vector<float>   *trackEta;
   vector<float>   *trackPhi;
   UInt_t          j1nTracks;
   Double_t        j1trk12PT;
   Double_t        j1trk1PT;
   Double_t        j1trk1Eta;
   Double_t        j1trk1Phi;
   Double_t        j1trk2PT;
   Double_t        j1trk2Eta;
   Double_t        j1trk2Phi;
   UInt_t          j2nTracks;
   Double_t        j2trk12PT;
   Double_t        j2trk1PT;
   Double_t        j2trk1Eta;
   Double_t        j2trk1Phi;
   Double_t        j2trk2PT;
   Double_t        j2trk2Eta;
   Double_t        j2trk2Phi;
   UInt_t          nGoodJets;
   Double_t        j1PT;
   Double_t        j1Eta;
   Double_t        j1Phi;
   Double_t        j1CHdFr;
   Double_t        j1NHdFr;
   Double_t        j1CEmFr;
   Double_t        j1NEmFr;
   Double_t        j1PhoEFr;
   Double_t        j1EleEFr;
   Double_t        j1MuEFr;
   Double_t        j1CMuEFr;
   UInt_t          j1nCons;
   Int_t           j1CMty;
   Int_t           j1NMty;
   Int_t           j1CHdMty;
   Int_t           j1NHdMty;
   Int_t           j1PhoMty;
   Int_t           j1EleMty;
   Int_t           j1MuMty;
   Double_t        j1etaWidth;
   Double_t        j1phiWidth;
   Double_t        j1etaWidthInECal;
   Double_t        j1phiWidthInECal;
   Double_t        j1etaWidthInHCal;
   Double_t        j1phiWidthInHCal;
   Double_t        j2PT;
   Double_t        j2Eta;
   Double_t        j2Phi;
   Double_t        j2CHdFr;
   Double_t        j2NHdFr;
   Double_t        j2CEmFr;
   Double_t        j2NEmFr;
   Double_t        j2PhoEFr;
   Double_t        j2EleEFr;
   Double_t        j2MuEFr;
   Double_t        j2CMuEFr;
   UInt_t          j2nCons;
   Int_t           j2CMty;
   Int_t           j2NMty;
   Int_t           j2CHdMty;
   Int_t           j2NHdMty;
   Int_t           j2PhoMty;
   Int_t           j2EleMty;
   Int_t           j2MuMty;
   Double_t        j2etaWidth;
   Double_t        j2phiWidth;
   Double_t        j2etaWidthInECal;
   Double_t        j2phiWidthInECal;
   Double_t        j2etaWidthInHCal;
   Double_t        j2phiWidthInHCal;
   Int_t           HLTPFMET300;
   Int_t           HLTPFMET170_HBHECleaned;
   Int_t           HLTPhoton165_HE10;
   Int_t           HLTPhoton175;
   Int_t           HLTPhoton75;
   Int_t           HLTPhoton90;
   Int_t           HLTPhoton120;
   Int_t           HLTPFJet40;
   Int_t           HLTPFJet60;
   Int_t           HLTPFJet80;
   // List of branches

   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_totalET;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_caloMET;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetEtaWidth;   //!
   TBranch        *b_jetPhiWidth;   //!
   TBranch        *b_jetEtaWidthInECal;   //!
   TBranch        *b_jetEtaWidthInHCal;   //!
   TBranch        *b_jetPhiWidthInECal;   //!
   TBranch        *b_jetPhiWidthInHCal;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_ElectronPassVetoID;   //!
   TBranch        *b_ElectronPassLooseID;   //!
   TBranch        *b_ElectronPassMediumID;   //!
   TBranch        *b_ElectronPassTightID;   //!
   TBranch        *b_ElectronPassHEEPID;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIsLooseID;   //!
   TBranch        *b_muIsMediumID;   //!
   TBranch        *b_muIsTightID;   //!
   TBranch        *b_muIsSoftID;   //!
   TBranch        *b_muIsHighPtID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaAtCalo;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_elePFMiniIso;   //!
   TBranch        *b_eleIDMVANonTrg;   //!
   TBranch        *b_eleIDMVATrg;   //!
   TBranch        *b_eledEtaseedAtVtx;   //!
   TBranch        *b_eleE1x5;   //!
   TBranch        *b_eleE2x5;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE1x5Full5x5;   //!
   TBranch        *b_eleE2x5Full5x5;   //!
   TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleisoChargedHadrons;   //!
   TBranch        *b_eleisoNeutralHadrons;   //!
   TBranch        *b_eleisoPhotons;   //!
   TBranch        *b_eleisoChargedFromPU;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_phopt;   //!
   TBranch        *b_phoeta;   //!
   TBranch        *b_phophi;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phohOverE;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoisoChargedHadrons;   //!
   TBranch        *b_phoWorstisoChargedHadrons;   //!
   TBranch        *b_phoisoNeutralHadrons;   //!
   TBranch        *b_phoisoPhotons;   //!
   TBranch        *b_trackPt;   //!
   TBranch        *b_trackEta;   //!
   TBranch        *b_trackPhi;   //!
   TBranch        *b_j1nTracks;   //!
   TBranch        *b_j1trk12PT;   //!
   TBranch        *b_j1trk1PT;   //!
   TBranch        *b_j1trk1Eta;   //!
   TBranch        *b_j1trk1Phi;   //!
   TBranch        *b_j1trk2PT;   //!
   TBranch        *b_j1trk2Eta;   //!
   TBranch        *b_j1trk2Phi;   //!
   TBranch        *b_j2nTracks;   //!
   TBranch        *b_j2trk12PT;   //!
   TBranch        *b_j2trk1PT;   //!
   TBranch        *b_j2trk1Eta;   //!
   TBranch        *b_j2trk1Phi;   //!
   TBranch        *b_j2trk2PT;   //!
   TBranch        *b_j2trk2Eta;   //!
   TBranch        *b_j2trk2Phi;   //!
   TBranch        *b_nGoodJets;   //!
   TBranch        *b_j1PT;   //!
   TBranch        *b_j1Eta;   //!
   TBranch        *b_j1Phi;   //!
   TBranch        *b_j1CHdFr;   //!
   TBranch        *b_j1NHdFr;   //!
   TBranch        *b_j1CEmFr;   //!
   TBranch        *b_j1NEmFr;   //!
   TBranch        *b_j1PhoEFr;   //!
   TBranch        *b_j1EleEFr;   //!
   TBranch        *b_j1MuEFr;   //!
   TBranch        *b_j1CMuEFr;   //!
   TBranch        *b_j1nCons;   //!
   TBranch        *b_j1CMty;   //!
   TBranch        *b_j1NMty;   //!
   TBranch        *b_j1CHdMty;   //!
   TBranch        *b_j1NHdMty;   //!
   TBranch        *b_j1PhoMty;   //!
   TBranch        *b_j1EleMty;   //!
   TBranch        *b_j1MuMty;   //!
   TBranch        *b_j1etaWidth;   //!
   TBranch        *b_j1phiWidth;   //!
   TBranch        *b_j1etaWidthInECal;   //!
   TBranch        *b_j1phiWidthInECal;   //!
   TBranch        *b_j1etaWidthInHCal;   //!
   TBranch        *b_j1phiWidthInHCal;   //!
   TBranch        *b_j2PT;   //!
   TBranch        *b_j2Eta;   //!
   TBranch        *b_j2Phi;   //!
   TBranch        *b_j2CHdFr;   //!
   TBranch        *b_j2NHdFr;   //!
   TBranch        *b_j2CEmFr;   //!
   TBranch        *b_j2NEmFr;   //!
   TBranch        *b_j2PhoEFr;   //!
   TBranch        *b_j2EleEFr;   //!
   TBranch        *b_j2MuEFr;   //!
   TBranch        *b_j2CMuEFr;   //!
   TBranch        *b_j2nCons;   //!
   TBranch        *b_j2CMty;   //!
   TBranch        *b_j2NMty;   //!
   TBranch        *b_j2CHdMty;   //!
   TBranch        *b_j2NHdMty;   //!
   TBranch        *b_j2PhoMty;   //!
   TBranch        *b_j2EleMty;   //!
   TBranch        *b_j2MuMty;   //!
   TBranch        *b_j2etaWidth;   //!
   TBranch        *b_j2phiWidth;   //!
   TBranch        *b_j2etaWidthInECal;   //!
   TBranch        *b_j2phiWidthInECal;   //!
   TBranch        *b_j2etaWidthInHCal;   //!
   TBranch        *b_j2phiWidthInHCal;   //!
   TBranch        *b_HLTPFMET300;   //!
   TBranch        *b_HLTPFMET170_HBHECleaned;   //!
   TBranch        *b_HLTPhoton165_HE10;   //!
   TBranch        *b_HLTPhoton175;   //!
   TBranch        *b_HLTPhoton75;   //!
   TBranch        *b_HLTPhoton90;   //!
   TBranch        *b_HLTPhoton120;   //!
   TBranch        *b_HLTPFJet40;   //!
   TBranch        *b_HLTPFJet60;   //!
   TBranch        *b_HLTPFJet80;   //!
   ZprimeJetsClass(const char* file1,const char* file2);
   virtual ~ZprimeJetsClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(Long64_t maxEvents, int reportEvery);
   virtual Bool_t   Notify();   
   virtual void     Show(Long64_t entry = -1);
   virtual void BookHistos(const char* file2);
   virtual double dR(double jetetaWidth, double jetphiWidth);
   virtual double deltaR(double eta1, double phi1, double eta2, double phi2);
   virtual void fillHistos(int histoNumber);
   virtual float DeltaPhi(float phi1, float phi2);
   virtual bool HighPtMuonVeto(double muonPtCut);
   virtual vector<int> getJetCand(double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut);
   virtual bool dPhiJetMETcut(std::vector<int> jets);
   virtual float dPhiJetMETmin(std::vector<int> jets);
   virtual bool electron_veto_looseID(int jet_index, float elePtCut);
   virtual bool muon_veto_looseID(int jet_index, float muPtCut);
};

#endif

#ifdef ZprimeJetsClass_cxx
ZprimeJetsClass::ZprimeJetsClass(const char* file1,const char* file2) 
{
  TChain *chain = new TChain("JetAnalyzer/JetTree");
  TString path = file1;
  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter nextlist(fileList);
  TSystemFile* filename;
  int fileNumber = 0;
  int maxFiles = -1;
  while ((filename = (TSystemFile*)nextlist()) && fileNumber >  maxFiles)
    {
      //Debug
    std::cout<<"file path found: "<<(path+filename->GetName())<<std::endl;
    std::cout<<"name: "<<(filename->GetName())<<std::endl;
    std::cout<<"fileNumber: "<<fileNumber<<std::endl;

     TString dataset = "JetAnalyzer_";
     TString  FullPathInputFile = (path+filename->GetName());
     TString name = filename->GetName();
     if(name.Contains(dataset))
       {
         std::cout<<"Adding FullPathInputFile to chain:"<<FullPathInputFile<<std::endl;
         std::cout<<std::endl;
         chain->Add(FullPathInputFile);
       }

      fileNumber++;
    }
  std::cout<<"All files added."<<std::endl;
  std::cout<<"Initializing chain."<<std::endl;
  Init(chain);
  BookHistos(file2);
}

ZprimeJetsClass::~ZprimeJetsClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   tree->Write();
   fileName->Close();
}

Int_t ZprimeJetsClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeJetsClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ZprimeJetsClass::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   jetPt = 0;
   jetEta = 0;
   jetEn = 0;
   jetPhi = 0;
   jetPFLooseId = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetHFHAE = 0;
   jetHFEME = 0;
   jetNConstituents = 0;
   jetEtaWidth = 0;
   jetPhiWidth = 0;
   jetEtaWidthInECal = 0;
   jetEtaWidthInHCal = 0;
   jetPhiWidthInECal = 0;
   jetPhiWidthInHCal = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   ElectronPassVetoID = 0;
   ElectronPassLooseID = 0;
   ElectronPassMediumID = 0;
   ElectronPassTightID = 0;
   ElectronPassHEEPID = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIsLooseID = 0;
   muIsMediumID = 0;
   muIsTightID = 0;
   muIsSoftID = 0;
   muIsHighPtID = 0;
   muD0 = 0;
   muDz = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleD0 = 0;
   eleDz = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaAtCalo = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIPhiIPhi = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   elePFMiniIso = 0;
   eleIDMVANonTrg = 0;
   eleIDMVATrg = 0;
   eledEtaseedAtVtx = 0;
   eleE1x5 = 0;
   eleE2x5 = 0;
   eleE5x5 = 0;
   eleE1x5Full5x5 = 0;
   eleE2x5Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   eleR9Full5x5 = 0;
   eleisoChargedHadrons = 0;
   eleisoNeutralHadrons = 0;
   eleisoPhotons = 0;
   eleisoChargedFromPU = 0;
   phopt = 0;
   phoeta = 0;
   phophi = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phohOverE = 0;
   phohasPixelSeed = 0;
   phoisoChargedHadrons = 0;
   phoWorstisoChargedHadrons = 0;
   phoisoNeutralHadrons = 0;
   phoisoPhotons = 0;
   trackPt = 0;
   trackEta = 0;
   trackPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("metFilters", &metFilters, &b_metFilters);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("nTrksPV", &nTrksPV, &b_nTrksPV);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("totalET", &totalET, &b_totalET);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("caloMET", &caloMET, &b_caloMET);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", &jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", &jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", &jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetEtaWidth", &jetEtaWidth, &b_jetEtaWidth);
   fChain->SetBranchAddress("jetPhiWidth", &jetPhiWidth, &b_jetPhiWidth);
   fChain->SetBranchAddress("jetEtaWidthInECal", &jetEtaWidthInECal, &b_jetEtaWidthInECal);
   fChain->SetBranchAddress("jetEtaWidthInHCal", &jetEtaWidthInHCal, &b_jetEtaWidthInHCal);
   fChain->SetBranchAddress("jetPhiWidthInECal", &jetPhiWidthInECal, &b_jetPhiWidthInECal);
   fChain->SetBranchAddress("jetPhiWidthInHCal", &jetPhiWidthInHCal, &b_jetPhiWidthInHCal);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("ElectronPassVetoID", &ElectronPassVetoID, &b_ElectronPassVetoID);
   fChain->SetBranchAddress("ElectronPassLooseID", &ElectronPassLooseID, &b_ElectronPassLooseID);
   fChain->SetBranchAddress("ElectronPassMediumID", &ElectronPassMediumID, &b_ElectronPassMediumID);
   fChain->SetBranchAddress("ElectronPassTightID", &ElectronPassTightID, &b_ElectronPassTightID);
   fChain->SetBranchAddress("ElectronPassHEEPID", &ElectronPassHEEPID, &b_ElectronPassHEEPID);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIsLooseID", &muIsLooseID, &b_muIsLooseID);
   fChain->SetBranchAddress("muIsMediumID", &muIsMediumID, &b_muIsMediumID);
   fChain->SetBranchAddress("muIsTightID", &muIsTightID, &b_muIsTightID);
   fChain->SetBranchAddress("muIsSoftID", &muIsSoftID, &b_muIsSoftID);
   fChain->SetBranchAddress("muIsHighPtID", &muIsHighPtID, &b_muIsHighPtID);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaAtCalo", &eledEtaAtCalo, &b_eledEtaAtCalo);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("elePFMiniIso", &elePFMiniIso, &b_elePFMiniIso);
   fChain->SetBranchAddress("eleIDMVANonTrg", &eleIDMVANonTrg, &b_eleIDMVANonTrg);
   fChain->SetBranchAddress("eleIDMVATrg", &eleIDMVATrg, &b_eleIDMVATrg);
   fChain->SetBranchAddress("eledEtaseedAtVtx", &eledEtaseedAtVtx, &b_eledEtaseedAtVtx);
   fChain->SetBranchAddress("eleE1x5", &eleE1x5, &b_eleE1x5);
   fChain->SetBranchAddress("eleE2x5", &eleE2x5, &b_eleE2x5);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE1x5Full5x5", &eleE1x5Full5x5, &b_eleE1x5Full5x5);
   fChain->SetBranchAddress("eleE2x5Full5x5", &eleE2x5Full5x5, &b_eleE2x5Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleisoChargedHadrons", &eleisoChargedHadrons, &b_eleisoChargedHadrons);
   fChain->SetBranchAddress("eleisoNeutralHadrons", &eleisoNeutralHadrons, &b_eleisoNeutralHadrons);
   fChain->SetBranchAddress("eleisoPhotons", &eleisoPhotons, &b_eleisoPhotons);
   fChain->SetBranchAddress("eleisoChargedFromPU", &eleisoChargedFromPU, &b_eleisoChargedFromPU);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("phopt", &phopt, &b_phopt);
   fChain->SetBranchAddress("phoeta", &phoeta, &b_phoeta);
   fChain->SetBranchAddress("phophi", &phophi, &b_phophi);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phohOverE", &phohOverE, &b_phohOverE);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoisoChargedHadrons", &phoisoChargedHadrons, &b_phoisoChargedHadrons);
   fChain->SetBranchAddress("phoWorstisoChargedHadrons", &phoWorstisoChargedHadrons, &b_phoWorstisoChargedHadrons);
   fChain->SetBranchAddress("phoisoNeutralHadrons", &phoisoNeutralHadrons, &b_phoisoNeutralHadrons);
   fChain->SetBranchAddress("phoisoPhotons", &phoisoPhotons, &b_phoisoPhotons);
   fChain->SetBranchAddress("trackPt", &trackPt, &b_trackPt);
   fChain->SetBranchAddress("trackEta", &trackEta, &b_trackEta);
   fChain->SetBranchAddress("trackPhi", &trackPhi, &b_trackPhi);
   fChain->SetBranchAddress("j1nTracks", &j1nTracks, &b_j1nTracks);
   fChain->SetBranchAddress("j1trk12PT", &j1trk12PT, &b_j1trk12PT);
   fChain->SetBranchAddress("j1trk1PT", &j1trk1PT, &b_j1trk1PT);
   fChain->SetBranchAddress("j1trk1Eta", &j1trk1Eta, &b_j1trk1Eta);
   fChain->SetBranchAddress("j1trk1Phi", &j1trk1Phi, &b_j1trk1Phi);
   fChain->SetBranchAddress("j1trk2PT", &j1trk2PT, &b_j1trk2PT);
   fChain->SetBranchAddress("j1trk2Eta", &j1trk2Eta, &b_j1trk2Eta);
   fChain->SetBranchAddress("j1trk2Phi", &j1trk2Phi, &b_j1trk2Phi);
   fChain->SetBranchAddress("j2nTracks", &j2nTracks, &b_j2nTracks);
   fChain->SetBranchAddress("j2trk12PT", &j2trk12PT, &b_j2trk12PT);
   fChain->SetBranchAddress("j2trk1PT", &j2trk1PT, &b_j2trk1PT);
   fChain->SetBranchAddress("j2trk1Eta", &j2trk1Eta, &b_j2trk1Eta);
   fChain->SetBranchAddress("j2trk1Phi", &j2trk1Phi, &b_j2trk1Phi);
   fChain->SetBranchAddress("j2trk2PT", &j2trk2PT, &b_j2trk2PT);
   fChain->SetBranchAddress("j2trk2Eta", &j2trk2Eta, &b_j2trk2Eta);
   fChain->SetBranchAddress("j2trk2Phi", &j2trk2Phi, &b_j2trk2Phi);
   fChain->SetBranchAddress("nGoodJets", &nGoodJets, &b_nGoodJets);
   fChain->SetBranchAddress("j1PT", &j1PT, &b_j1PT);
   fChain->SetBranchAddress("j1Eta", &j1Eta, &b_j1Eta);
   fChain->SetBranchAddress("j1Phi", &j1Phi, &b_j1Phi);
   fChain->SetBranchAddress("j1CHdFr", &j1CHdFr, &b_j1CHdFr);
   fChain->SetBranchAddress("j1NHdFr", &j1NHdFr, &b_j1NHdFr);
   fChain->SetBranchAddress("j1CEmFr", &j1CEmFr, &b_j1CEmFr);
   fChain->SetBranchAddress("j1NEmFr", &j1NEmFr, &b_j1NEmFr);
   fChain->SetBranchAddress("j1PhoEFr", &j1PhoEFr, &b_j1PhoEFr);
   fChain->SetBranchAddress("j1EleEFr", &j1EleEFr, &b_j1EleEFr);
   fChain->SetBranchAddress("j1MuEFr", &j1MuEFr, &b_j1MuEFr);
   fChain->SetBranchAddress("j1CMuEFr", &j1CMuEFr, &b_j1CMuEFr);
   fChain->SetBranchAddress("j1nCons", &j1nCons, &b_j1nCons);
   fChain->SetBranchAddress("j1CMty", &j1CMty, &b_j1CMty);
   fChain->SetBranchAddress("j1NMty", &j1NMty, &b_j1NMty);
   fChain->SetBranchAddress("j1CHdMty", &j1CHdMty, &b_j1CHdMty);
   fChain->SetBranchAddress("j1NHdMty", &j1NHdMty, &b_j1NHdMty);
   fChain->SetBranchAddress("j1PhoMty", &j1PhoMty, &b_j1PhoMty);
   fChain->SetBranchAddress("j1EleMty", &j1EleMty, &b_j1EleMty);
   fChain->SetBranchAddress("j1MuMty", &j1MuMty, &b_j1MuMty);
   fChain->SetBranchAddress("j1etaWidth", &j1etaWidth, &b_j1etaWidth);
   fChain->SetBranchAddress("j1phiWidth", &j1phiWidth, &b_j1phiWidth);
   fChain->SetBranchAddress("j1etaWidthInECal", &j1etaWidthInECal, &b_j1etaWidthInECal);
   fChain->SetBranchAddress("j1phiWidthInECal", &j1phiWidthInECal, &b_j1phiWidthInECal);
   fChain->SetBranchAddress("j1etaWidthInHCal", &j1etaWidthInHCal, &b_j1etaWidthInHCal);
   fChain->SetBranchAddress("j1phiWidthInHCal", &j1phiWidthInHCal, &b_j1phiWidthInHCal);
   fChain->SetBranchAddress("j2PT", &j2PT, &b_j2PT);
   fChain->SetBranchAddress("j2Eta", &j2Eta, &b_j2Eta);
   fChain->SetBranchAddress("j2Phi", &j2Phi, &b_j2Phi);
   fChain->SetBranchAddress("j2CHdFr", &j2CHdFr, &b_j2CHdFr);
   fChain->SetBranchAddress("j2NHdFr", &j2NHdFr, &b_j2NHdFr);
   fChain->SetBranchAddress("j2CEmFr", &j2CEmFr, &b_j2CEmFr);
   fChain->SetBranchAddress("j2NEmFr", &j2NEmFr, &b_j2NEmFr);
   fChain->SetBranchAddress("j2PhoEFr", &j2PhoEFr, &b_j2PhoEFr);
   fChain->SetBranchAddress("j2EleEFr", &j2EleEFr, &b_j2EleEFr);
   fChain->SetBranchAddress("j2MuEFr", &j2MuEFr, &b_j2MuEFr);
   fChain->SetBranchAddress("j2CMuEFr", &j2CMuEFr, &b_j2CMuEFr);
   fChain->SetBranchAddress("j2nCons", &j2nCons, &b_j2nCons);
   fChain->SetBranchAddress("j2CMty", &j2CMty, &b_j2CMty);
   fChain->SetBranchAddress("j2NMty", &j2NMty, &b_j2NMty);
   fChain->SetBranchAddress("j2CHdMty", &j2CHdMty, &b_j2CHdMty);
   fChain->SetBranchAddress("j2NHdMty", &j2NHdMty, &b_j2NHdMty);
   fChain->SetBranchAddress("j2PhoMty", &j2PhoMty, &b_j2PhoMty);
   fChain->SetBranchAddress("j2EleMty", &j2EleMty, &b_j2EleMty);
   fChain->SetBranchAddress("j2MuMty", &j2MuMty, &b_j2MuMty);
   fChain->SetBranchAddress("j2etaWidth", &j2etaWidth, &b_j2etaWidth);
   fChain->SetBranchAddress("j2phiWidth", &j2phiWidth, &b_j2phiWidth);
   fChain->SetBranchAddress("j2etaWidthInECal", &j2etaWidthInECal, &b_j2etaWidthInECal);
   fChain->SetBranchAddress("j2phiWidthInECal", &j2phiWidthInECal, &b_j2phiWidthInECal);
   fChain->SetBranchAddress("j2etaWidthInHCal", &j2etaWidthInHCal, &b_j2etaWidthInHCal);
   fChain->SetBranchAddress("j2phiWidthInHCal", &j2phiWidthInHCal, &b_j2phiWidthInHCal);
   fChain->SetBranchAddress("HLTPFMET300", &HLTPFMET300, &b_HLTPFMET300);
   fChain->SetBranchAddress("HLTPFMET170_HBHECleaned", &HLTPFMET170_HBHECleaned, &b_HLTPFMET170_HBHECleaned);
   fChain->SetBranchAddress("HLTPhoton165_HE10", &HLTPhoton165_HE10, &b_HLTPhoton165_HE10);
   fChain->SetBranchAddress("HLTPhoton175", &HLTPhoton175, &b_HLTPhoton175);
   fChain->SetBranchAddress("HLTPhoton75", &HLTPhoton75, &b_HLTPhoton75);
   fChain->SetBranchAddress("HLTPhoton90", &HLTPhoton90, &b_HLTPhoton90);
   fChain->SetBranchAddress("HLTPhoton120", &HLTPhoton120, &b_HLTPhoton120);
   fChain->SetBranchAddress("HLTPFJet40", &HLTPFJet40, &b_HLTPFJet40);
   fChain->SetBranchAddress("HLTPFJet60", &HLTPFJet60, &b_HLTPFJet60);
   fChain->SetBranchAddress("HLTPFJet80", &HLTPFJet80, &b_HLTPFJet80);
   Notify();
}

Bool_t ZprimeJetsClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeJetsClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeJetsClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeJetsClass_cxx
