//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug 28 11:27:23 2016 by ROOT version 6.06/01
// from TTree JetTree/Jet data for analysis
// found on file: /hdfs/store/user/uhussain/Zprime_Ntuples_Aug30/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_QCD_HT100to200/160830_160139/0000/JetAnalyzerMC_999.root
//////////////////////////////////////////////////////////

#ifndef ZprimeJetsClass_MC_h
#define ZprimeJetsClass_MC_h

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
class ZprimeJetsClass_MC {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::vector<unsigned int> event_;
   std::vector<double> event_info;
   TFile *fileName;
   TTree *tree;

   TH1F *h_dphimin, *h_deltar_cut, *h_deltar, *h_HT[8], *h_nJets[8], *h_nGoodJets[8], *h_pfMET[8], *h_pfMETPhi[8], *h_j1Pt[8], *h_j1Eta[8], *h_j1Phi[8], *h_j1etaWidth[8], *h_j1phiWidth[8], *h_j1etaWidthInECal[8], *h_j1etaWidthInHCal[8], *h_j1phiWidthInECal[8], *h_j1phiWidthInHCal[8], *h_j1nCons[8], *h_j1CEF[8], *h_j1NEF[8], *h_j1CHF[8], *h_j1NHF[8];
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumis;
   UInt_t          metFilters;
   Int_t           npv;
   Int_t           nTrksPV;
   Float_t         genEventWeight;
   vector<float>   *genlheWeights;
   Int_t           numWeights;
   Float_t         nTrueVertices;
   Int_t           NUP;
   Int_t           numGenJets;
   Float_t         genHT;
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
   Int_t           ElectronPassVetoID;
   Int_t           ElectronPassLooseID;
   Int_t           ElectronPassMediumID;
   Int_t           ElectronPassTightID;
   Int_t           ElectronPassHEEPID;
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
   UInt_t          nGoodJets;
   Double_t        j1PT;
   Double_t        j1Eta;
   Double_t        j1Phi;
   Double_t        j1CHdFr;
   Double_t        j1NHdFr;
   Double_t        j1CEmFr;
   Double_t        j1NEmFr;
   UInt_t          j1nCons;
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
   UInt_t          j2nCons;
   Double_t        j2etaWidth;
   Double_t        j2phiWidth;
   Double_t        j2etaWidthInECal;
   Double_t        j2phiWidthInECal;
   Double_t        j2etaWidthInHCal;
   Double_t        j2phiWidthInHCal;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_metFilters;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_nTrksPV;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genlheWeights;   //!
   TBranch        *b_numWeights;   //!
   TBranch        *b_nTrueVertices;   //!
   TBranch        *b_NUP;   //!
   TBranch        *b_numGenJets;   //!
   TBranch        *b_genHT;   //!
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
   TBranch        *b_nGoodJets;   //!
   TBranch        *b_j1PT;   //!
   TBranch        *b_j1Eta;   //!
   TBranch        *b_j1Phi;   //!
   TBranch        *b_j1CHdFr;   //!
   TBranch        *b_j1NHdFr;   //!
   TBranch        *b_j1CEmFr;   //!
   TBranch        *b_j1NEmFr;   //!
   TBranch        *b_j1nCons;   //!
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
   TBranch        *b_j2nCons;   //!
   TBranch        *b_j2etaWidth;   //!
   TBranch        *b_j2phiWidth;   //!
   TBranch        *b_j2etaWidthInECal;   //!
   TBranch        *b_j2phiWidthInECal;   //!
   TBranch        *b_j2etaWidthInHCal;   //!
   TBranch        *b_j2phiWidthInHCal;   //!

   ZprimeJetsClass_MC(const char* file1,const char* file2);
   virtual ~ZprimeJetsClass_MC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop(Long64_t maxEvents, int reportEvery);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void BookHistos(const char* file2);
   virtual double dR(double jetetaWidth, double jetphiWidth);
   virtual void fillHistos(int histoNumber);
   virtual float DeltaPhi(float phi1, float phi2);
   virtual bool HighPtMuonVeto(double muonPtCut);
   virtual vector<int> getJetCand(double jetPtCut, double jetEtaCut, double jetNHFCut, double jetCHFCut);
   virtual vector<int> dPhiJetMETcut(std::vector<int> jets);
   virtual float dPhiJetMETmin(std::vector<int> jets);
};

#endif

#ifdef ZprimeJetsClass_MC_cxx
ZprimeJetsClass_MC::ZprimeJetsClass_MC(const char* file1,const char* file2)
{
  TChain *chain = new TChain("JetAnalyzerMC/JetTree");
  TString path = file1;
  TSystemDirectory sourceDir("hi",path);
  TList* fileList = sourceDir.GetListOfFiles();
  TIter nextlist(fileList);
  TSystemFile* filename;
  int fileNumber = 0;
  int maxFiles = -1;
  while ((filename = (TSystemFile*)nextlist()) && fileNumber >  maxFiles)
    {
    std::cout<<"file path found: "<<(path+filename->GetName())<<std::endl;
    std::cout<<"name: "<<(filename->GetName())<<std::endl;
    std::cout<<"fileNumber: "<<fileNumber<<std::endl;

     TString dataset = "JetAnalyzerMC_";
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

ZprimeJetsClass_MC::~ZprimeJetsClass_MC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
   fileName->cd();
   fileName->Write();
   tree->Write();
   fileName->Close();
}

Int_t ZprimeJetsClass_MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ZprimeJetsClass_MC::LoadTree(Long64_t entry)
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

void ZprimeJetsClass_MC::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   genlheWeights = 0;
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
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genWeight);
   fChain->SetBranchAddress("genlheWeights", &genlheWeights, &b_genlheWeights);
   fChain->SetBranchAddress("numWeights", &numWeights, &b_numWeights);
   fChain->SetBranchAddress("nTrueVertices", &nTrueVertices, &b_nTrueVertices);
   fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
   fChain->SetBranchAddress("numGenJets", &numGenJets, &b_numGenJets);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
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
   fChain->SetBranchAddress("nGoodJets", &nGoodJets, &b_nGoodJets);
   fChain->SetBranchAddress("j1PT", &j1PT, &b_j1PT);
   fChain->SetBranchAddress("j1Eta", &j1Eta, &b_j1Eta);
   fChain->SetBranchAddress("j1Phi", &j1Phi, &b_j1Phi);
   fChain->SetBranchAddress("j1CHdFr", &j1CHdFr, &b_j1CHdFr);
   fChain->SetBranchAddress("j1NHdFr", &j1NHdFr, &b_j1NHdFr);
   fChain->SetBranchAddress("j1CEmFr", &j1CEmFr, &b_j1CEmFr);
   fChain->SetBranchAddress("j1NEmFr", &j1NEmFr, &b_j1NEmFr);
   fChain->SetBranchAddress("j1nCons", &j1nCons, &b_j1nCons);
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
   fChain->SetBranchAddress("j2nCons", &j2nCons, &b_j2nCons);
   fChain->SetBranchAddress("j2etaWidth", &j2etaWidth, &b_j2etaWidth);
   fChain->SetBranchAddress("j2phiWidth", &j2phiWidth, &b_j2phiWidth);
   fChain->SetBranchAddress("j2etaWidthInECal", &j2etaWidthInECal, &b_j2etaWidthInECal);
   fChain->SetBranchAddress("j2phiWidthInECal", &j2phiWidthInECal, &b_j2phiWidthInECal);
   fChain->SetBranchAddress("j2etaWidthInHCal", &j2etaWidthInHCal, &b_j2etaWidthInHCal);
   fChain->SetBranchAddress("j2phiWidthInHCal", &j2phiWidthInHCal, &b_j2phiWidthInHCal);
   Notify();
}

Bool_t ZprimeJetsClass_MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ZprimeJetsClass_MC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ZprimeJetsClass_MC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ZprimeJetsClass_MC_cxx
