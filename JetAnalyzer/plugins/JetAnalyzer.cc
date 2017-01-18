// -*- C++ -*-
//
// Package:    LightZPrimeAnalysis/JetAnalyzer
// Class:      JetAnalyzer
// 
//cclass JetAnalyzer JetAnalyzer.cc LightZPrimeAnalysis/JetAnalyzer/plugins/JetAnalyzer.cc

// Description: [one line class summary]



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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "LightZPrimeAnalysis/JetWidthCalculator/interface/JetWidthCalculator.hh"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class JetAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit JetAnalyzer(const edm::ParameterSet&);
    ~JetAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    // ----------member data ---------------------------

    edm::EDGetTokenT< vector<reco::PFCandidate> > pfCandsToken;
    edm::EDGetTokenT< vector<reco::PFJet> > pfJetsToken;
    edm::EDGetTokenT< vector<pat::Jet> > JetsToken;
    edm::EDGetTokenT< vector<reco::CaloJet> > caloJetsToken;
    edm::EDGetTokenT< vector<reco::PFMET> > pfMETsToken;
    edm::EDGetTokenT< vector<reco::CaloMET> > caloMETToken;  
    edm::EDGetToken electronCollection_;
    edm::EDGetTokenT<edm::View<pat::Muon> > muonToken;  
    //edm::EDGetTokenT< vector<reco::Muon> > muonToken;  
    edm::EDGetTokenT<vector<reco::Track> > trackToken;  
    edm::EDGetTokenT<bool> globalHandle_;
    edm::EDGetTokenT<bool> hcalNoiseHandle_;
    edm::EDGetTokenT<bool> hcalIsoNoiseHandle_;
    edm::EDGetTokenT<bool> eCALTPHandle_;
    edm::EDGetTokenT<bool> bADSCHandle_;
    edm::EDGetTokenT<bool> BadChCandFilterToken_;
    edm::EDGetTokenT<bool> BadPFMuonFilterToken_;
    edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;//we need this to store genHT which is required to stitch inclusive samples with HT-binned samples
    edm::EDGetTokenT<GenEventInfoProduct> generatorLabel_;//we need this to store neagtive weights information
    edm::EDGetTokenT<vector<PileupSummaryInfo> > puCollection_; //Required for Pileup reweighting
    edm::EDGetTokenT<vector<reco::GenParticle> > genParticlesCollection_; //genParticles info required for boson Pt reweighting
    edm::EDGetToken photonsToken_;

    edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
    edm::EDGetTokenT<vector<reco::Vertex> > vtxToken_;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::EDGetTokenT<double> rhoToken_; 
    edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    // elecontr ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken_;

    edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
    edm::EDGetTokenT<edm::ValueMap<float> > phoWorstChargedIsolationToken_; 
    edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
    edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 

    void setbit(UShort_t& x, UShort_t bit) {
      UShort_t a = 1;
      x |= (a << bit);
    }

    //some genlevel variables
    int              nup_;  
    int              numGenJets_; 

    vector<float>    pdf_;
    float            pthat_;
    float            processID_;
    float            genWeight_;
    float            genHT_;

    Int_t            nPUInfo_;
    vector<int>      nPU_;
    vector<int>      puBX_;
    vector<float>    puTrue_;

    Int_t            nMC_;
    vector<int>      mcPID;
    vector<float>    mcVtx;
    vector<float>    mcVty;
    vector<float>    mcVtz;
    vector<float>    mcPt;
    vector<float>    mcMass;
    vector<float>    mcEta;
    vector<float>    mcPhi;
    vector<float>    mcE;
    vector<float>    mcEt;
    vector<int>      mcIndex;
    vector<UShort_t> mcStatusFlag;
    vector<int>      mcStatus;

    double bosonMass;
    double bosonPt;
    int neutrinos;
    vector<int>      BosonPID;
    vector<float>    BosonVtx;
    vector<float>    BosonVty;
    vector<float>    BosonVtz;
    vector<float>    BosonPt;
    vector<float>    BosonMass;
    vector<float>    BosonEta;
    vector<float>    BosonPhi;
    vector<float>    BosonE;
    vector<float>    BosonEt;
    vector<int>      BosonStatus;
    //some must have variables for tuples
    int              run_;
    int              event_;
    int              lumis_;
    int              npv_;
    int              nTrksPV_;
    int              nVtx_;
    Float_t          rho_;  
    //jet variables

    vector<float> calojetPt_;
    vector<float> calojetEn_;
    vector<float> calojetEta_;
    vector<float> calojetPhi_;
    vector<float> calojetEEF_;
    vector<float> calojetHEF_;
    vector<float> calojetTowersArea_;
    vector<float> calojetMaxEInEmTowers_;
    vector<float> calojetMaxEInHadTowers_;

    //double jet1Ptdiff; //variable to test JEC
    vector<float> jetPx_;
    vector<float> jetPy_;
    vector<float> jetPz_;
    vector<float> jetPt_;
    vector<float> jetEn_;
    vector<float> jetEta_;
    vector<float> jetPhi_;
    vector<float> jetCHF_;
    vector<float> jetNHF_;
    vector<float> jetCEF_;
    vector<float> jetNEF_;
    vector<int>   jetNCH_;
    vector<int>   jetNNP_;
    vector<float> jetMUF_;
    vector<float> jetHFHAE_;
    vector<float> jetHFEME_;
    vector<int>   jetNConstituents_;
    vector<float> jetCSV2BJetTags_; // recommended
    vector<float> jetJetProbabilityBJetTags_;
    vector<float> jetpfCombinedMVAV2BJetTags_;  
    vector<float> jetEtaWidth_;
    vector<float> jetPhiWidth_;
    vector<float> jetEtaWidthInECal_;
    vector<float> jetEtaWidthInHCal_;
    vector<float> jetPhiWidthInECal_;
    vector<float> jetPhiWidthInHCal_;

    vector<float> jetPhotonEnergyFraction_;
    vector<float> jetJetArea_;
    vector<float> jetMt_;
    vector<float> jetMass_;
    vector<float> jetMaxDistance_;
    vector<float> jetPhiPhiMoment_;
    vector<float> jetConstituentEtaPhiSpread_;
    vector<float> jetConstituentPtDistribution_;
    vector<float> jetPileup_;
    vector<int>   jetN60_;
    vector<int>   jetN90_;

    vector<bool>  jetPFLooseId_;

    //electron variables
    Int_t          nEle_;
    vector<float>  elePt_;
    vector<float>  eleEta_;
    vector<float>  elePhi_;
    vector<int>    PassVeto_;
    vector<int>    PassLoose_;
    vector<int>    PassMedium_;
    vector<int>    PassTight_;
    vector<int>    PassHEEP_; 


    vector<int>    eleCharge_;
    vector<int>    eleChargeConsistent_;
    vector<float>  eleEn_;
    vector<float>  eleSCEn_;
    vector<float>  eleESEn_;
    vector<float>  eleESEnP1_;
    vector<float>  eleESEnP2_;
    vector<float>  eleESEnP1Raw_;
    vector<float>  eleESEnP2Raw_;
    vector<float>  eleD0_;
    vector<float>  eleDz_;
    vector<float>  eleR9_;
    vector<float>  eleCalibPt_;
    vector<float>  eleCalibEn_;
    vector<float>  eleSCEta_;
    vector<float>  eleSCPhi_;
    vector<float>  eleSCRawEn_;
    vector<float>  eleSCEtaWidth_;
    vector<float>  eleSCPhiWidth_;
    vector<float>  eleHoverE_;
    vector<float>  eleEoverP_;
    vector<float>  eleEoverPout_;
    vector<float>  eleEoverPInv_;
    vector<float>  eleBrem_;
    vector<float>  eledEtaAtVtx_;
    vector<float>  eledPhiAtVtx_;
    vector<float>  eledEtaAtCalo_;
    vector<float>  eleSigmaIEtaIEta_;
    vector<float>  eleSigmaIPhiIPhi_;
    vector<float>  eleSigmaIEtaIEtaFull5x5_;
    vector<float>  eleSigmaIPhiIPhiFull5x5_;
    vector<int>    eleConvVeto_;
    vector<int>    eleMissHits_;
    vector<float>  eleESEffSigmaRR_;
    vector<float>  elePFClusEcalIso_;
    vector<float>  elePFClusHcalIso_;
    vector<float>  elePFMiniIso_;
    vector<float>  eleIDMVANonTrg_;
    vector<float>  eleIDMVATrg_;
    vector<float>  eledEtaseedAtVtx_;
    vector<float>  eleE1x5_;
    vector<float>  eleE2x5_;
    vector<float>  eleE5x5_;
    vector<float>  eleE1x5Full5x5_;
    vector<float>  eleE2x5Full5x5_;
    vector<float>  eleE5x5Full5x5_;
    vector<float>  eleR9Full5x5_;

    std::vector<Float_t> eleisoChargedHadrons_;
    std::vector<Float_t> eleisoNeutralHadrons_;
    std::vector<Float_t> eleisoPhotons_;
    std::vector<Float_t> eleisoChargedFromPU_;


    Int_t nPhotons_;

    std::vector<Float_t> phopt_;
    std::vector<Float_t> phoeta_;
    std::vector<Float_t> phophi_;

    std::vector<Float_t> phopx_;
    std::vector<Float_t> phopy_;
    std::vector<Float_t> phopz_;
    std::vector<Float_t> phoE_;
    std::vector<Float_t> phor9_;

    // Variables typically used for cut based photon ID
    std::vector<Float_t>  phoSigmaIEtaIEtaFull5x5_ ;
    std::vector<Float_t>  phohOverE_;
    std::vector<Int_t>    phohasPixelSeed_;

    std::vector<Float_t> phoisoChargedHadrons_;
    std::vector<Float_t> phoWorstisoChargedHadrons_;
    std::vector<Float_t> phoisoNeutralHadrons_;
    std::vector<Float_t> phoisoPhotons_;





    //muon variables 
    Int_t          nMu_;
    vector<float>  muPt_;
    vector<float>  muEn_;
    vector<float>  muEta_;
    vector<float>  muPhi_;
    vector<int>    muCharge_;
    vector<int>    muType_;

    vector<Bool_t> muIsLooseID_;
    vector<Bool_t> muIsMediumID_;
    vector<Bool_t> muIsTightID_;
    vector<Bool_t> muIsSoftID_;
    vector<Bool_t> muIsHighPtID_;
    vector<float>  muD0_;
    vector<float>  muDz_;
    vector<float>  muChi2NDF_;
    vector<float>  muInnerD0_;
    vector<float>  muInnerDz_;
    vector<int>    muTrkLayers_;
    vector<int>    muPixelLayers_;
    vector<int>    muPixelHits_;
    vector<int>    muMuonHits_;
    vector<int>    muStations_;
    vector<int>    muMatches_;
    vector<int>    muTrkQuality_;
    vector<float>  muIsoTrk_;
    vector<float>  muPFChIso_;
    vector<float>  muPFPhoIso_;
    vector<float>  muPFNeuIso_;
    vector<float>  muPFPUIso_;
    vector<float>  muInnervalidFraction_;
    vector<float>  musegmentCompatibility_;
    vector<float>  muchi2LocalPosition_;
    vector<float>  mutrkKink_;
    vector<float>  muBestTrkPtError_;
    vector<float>  muBestTrkPt_;



    //track variables
    vector<float>  trkPt_;
    vector<float>  trkEta_;
    vector<float>  trkPhi_;

    uint32_t metFilters_;
    double totalET;
    double HT;
    float pfMET;
    float  pfMETPhi;
    float  pfMETsumEt_;
    float  pfMETmEtSig_;
    float  pfMETSig_;
    float caloMET;
    double j1PT;
    double j1Eta;
    double j1Phi;
    double j1CHdFr;
    double j1NHdFr;
    double j1CEmFr;
    double j1NEmFr;
    double j1PhoEFr;
    double j1EleEFr;
    double j1MuEFr;
    double j1CMuEFr;
    double j1etaWidth;
    double j1phiWidth;
    double j1etaWidthInECal;
    double j1phiWidthInECal;
    double j1etaWidthInHCal;
    double j1phiWidthInHCal;
    double j1PFCand12PtSum;
    double j1PFCandsPtSum;
    double j1PFCand12Ratio;
    //track variables associated with leading Jet
    uint32_t j1nTracks;
    bool j1trk1Quality;
    int j1trk1Charge;
    double j1dRCand1Trk1;
    int j1Trk1PFCand1Match;
    bool j1trk2Quality;
    int j1trk2Charge;
    double j1dRCand2Trk2;
    int j1Trk2PFCand2Match;
    double j1trk1PtError;
    double j1trk2PtError;
    double j1ptTrk1Ratio;//j1trk1PtError/j1trk1PT
    double j1ptTrk2Ratio;
    double j1trk12PT;
    double j1trk12Ratio;//(j1trk12PT/j1pT)
    double j1trk1PT;
    double j1trk1Eta;
    double j1trk1Phi;
    double j1trk2PT;
    double j1trk2Eta;
    double j1trk2Phi;
    //track variables associated with second leading Jet
    uint32_t j2nTracks;
    bool j2trk1Quality;
    int j2trk1Charge;
    double j2dRCand1Trk1;
    int j2Trk1PFCand1Match;
    bool j2trk2Quality;
    int j2trk2Charge;
    double j2dRCand2Trk2;
    int j2Trk2PFCand2Match;
    double j2trk1PtError;
    double j2trk2PtError;
    double j2ptTrk1Ratio;//j2trk1PtError/j2trk1PT
    double j2ptTrk2Ratio;
    double j2trk12PT;
    double j2trk12Ratio;//(j2trk12PT/j1pT)
    double j2trk1PT;
    double j2trk1Eta;
    double j2trk1Phi;
    double j2trk2PT;
    double j2trk2Eta;
    double j2trk2Phi;

    double j2PT;
    double j2Eta;
    double j2Phi;
    double j2CHdFr;
    double j2NHdFr;
    double j2CEmFr;
    double j2NEmFr;
    double j2PhoEFr;
    double j2EleEFr;
    double j2MuEFr;
    double j2CMuEFr;
    double j2etaWidth;
    double j2phiWidth;
    double j2etaWidthInECal;
    double j2phiWidthInECal;
    double j2etaWidthInHCal;
    double j2phiWidthInHCal;
    double j2PFCand12PtSum;
    double j2PFCandsPtSum;
    double j2PFCand12Ratio;

    uint32_t nJets;
    uint32_t nGoodJets;

    uint32_t j1nCons;
    int j1nCarrying60;
    int j1nCarrying90;
    int j1nCarrying92;
    int j1nCarrying94;
    int j1nCarrying96;
    int j1nCarrying98;
    uint32_t j2nCons;
    int j2nCarrying60;
    int j2nCarrying90;
    int j2nCarrying92;
    int j2nCarrying94;
    int j2nCarrying96;
    int j2nCarrying98;
    int j1CMty;
    int j1NMty;
    int j1CHdMty;
    int j1NHdMty;
    int j1PhoMty;
    int j1EleMty;
    int j1MuMty;
    int j2CMty;
    int j2NMty;
    int j2CHdMty;
    int j2NHdMty;
    int j2PhoMty;
    int j2EleMty;
    int j2MuMty;
    int HLTMET300_;
    int HLTMET170_HBHE;
    int HLTPhoton165_HE10_;
    int HLTPhoton175_;
    int HLTPhoton75_;
    int HLTPhoton90_;
    int HLTPhoton120_;
    int HLTPFJet40_;
    int HLTPFJet60_;
    int HLTPFJet80_;

    TTree* tree;
    TH1F    *hPU_;
    TH1F    *hPUTrue_;
    TH1F    *hGenWeight_;
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
JetAnalyzer::JetAnalyzer(const edm::ParameterSet& iConfig):

  pfJetsToken(consumes< vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jetsTag")))
{

  pfCandsToken = consumes< vector<reco::PFCandidate> >(edm::InputTag("particleFlow"));
  //pfJetsToken = consumes< vector<reco::PFJet> >(edm::InputTag("jetsTag"));
  JetsToken = consumes< vector<pat::Jet> >(edm::InputTag("selectedPatJets"));
  caloJetsToken = consumes< vector<reco::CaloJet> >(edm::InputTag("ak4CaloJets"));
  pfMETsToken = consumes< vector<reco::PFMET> >(edm::InputTag("pfMet"));
  caloMETToken = consumes< vector<reco::CaloMET> >(edm::InputTag("caloMet"));  
  electronCollection_ = mayConsume<edm::View<reco::GsfElectron> >(edm::InputTag("gedGsfElectrons"));
  photonsToken_ = mayConsume<edm::View<reco::Photon> >(edm::InputTag("gedPhotons"));
  muonToken = consumes<edm::View<pat::Muon> >(edm::InputTag("selectedPatMuons"));
  //muonToken = consumes< vector<reco::Muon> >(edm::InputTag("muons"));
  trackToken = consumes< vector<reco::Track> >(edm::InputTag("generalTracks"));  
  trgResultsLabel_ = consumes<edm::TriggerResults>  (edm::InputTag("TriggerResults", "", "HLT"));
  globalHandle_= consumes<bool>(edm::InputTag("globalTightHalo2016Filter", ""));  
  hcalNoiseHandle_ = consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult"));
  hcalIsoNoiseHandle_ = consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"));
  eCALTPHandle_ =  consumes<bool>(edm::InputTag("EcalDeadCellTriggerPrimitiveFilter", "")); 
  bADSCHandle_ = consumes<bool>(edm::InputTag("eeBadScFilter", ""));
  BadChCandFilterToken_= consumes<bool>(edm::InputTag("BadChargedCandidateFilter"));
  BadPFMuonFilterToken_= consumes<bool>(edm::InputTag("BadPFMuonFilter"));  
  lheEventProductToken_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  generatorLabel_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  puCollection_ = consumes<vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));
  genParticlesCollection_ = consumes<vector<reco::GenParticle> >(edm::InputTag("genParticles")); 
  vtxToken_ = consumes< vector<reco::Vertex> >(edm::InputTag("offlinePrimaryVertices"));
  beamSpotToken_    = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
  rhoToken_    = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
  conversionsToken_ = mayConsume< reco::ConversionCollection>(edm::InputTag("allConversions"));
  eleVetoIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-veto"));
  eleLooseIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-loose"));
  eleMediumIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-medium"));
  eleTightIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-tight")); 
  eleHEEPIdMapToken_ = consumes<edm::ValueMap<bool> >(edm::InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"));  
  phoChargedIsolationToken_ =consumes <edm::ValueMap<float> >(edm::InputTag("photonIDValueMapProducer:phoChargedIsolation"));
  phoWorstChargedIsolationToken_ =consumes <edm::ValueMap<float> >(edm::InputTag("photonIDValueMapProducer:phoWorstChargedIsolation"));
  phoNeutralHadronIsolationToken_ =consumes <edm::ValueMap<float> >(edm::InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"));
  phoPhotonIsolationToken_= consumes <edm::ValueMap<float> >(edm::InputTag("photonIDValueMapProducer:phoPhotonIsolation"));

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("JetTree", "Jet data for analysis");

  tree->Branch("run",     &run_);
  tree->Branch("event",   &event_);
  tree->Branch("lumis",   &lumis_);  
  tree->Branch("metFilters", &metFilters_);
  tree->Branch("npv",&npv_,"npv/I");
  tree->Branch("nTrksPV",&nTrksPV_);
  tree->Branch("rho",&rho_);
  tree->Branch("NUP", &nup_, "NUP/I");
  tree->Branch("numGenJets", &numGenJets_, "numGenJets/I");

  tree->Branch("pdf",           &pdf_);
  tree->Branch("pthat",         &pthat_);
  tree->Branch("processID",     &processID_);
  tree->Branch("genWeight",     &genWeight_);
  tree->Branch("genHT",         &genHT_);

  tree->Branch("nPUInfo",       &nPUInfo_);
  tree->Branch("nPU",           &nPU_);
  tree->Branch("puBX",          &puBX_);
  tree->Branch("puTrue",        &puTrue_);

  tree->Branch("nMC",          &nMC_);
  tree->Branch("mcPID",        &mcPID);
  tree->Branch("mcVtx",        &mcVtx);
  tree->Branch("mcVty",        &mcVty);
  tree->Branch("mcVtz",        &mcVtz);
  tree->Branch("mcPt",         &mcPt);
  tree->Branch("mcMass",       &mcMass);
  tree->Branch("mcEta",        &mcEta);
  tree->Branch("mcPhi",        &mcPhi);
  tree->Branch("mcE",          &mcE);
  tree->Branch("mcEt",         &mcEt);
  tree->Branch("mcIndex",      &mcIndex);
  tree->Branch("mcStatusFlag", &mcStatusFlag); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
  tree->Branch("mcStatus",     &mcStatus);     // status of the particle

  hPU_        = fs->make<TH1F>("hPU",        "number of pileup",      200,  0, 200);
  hPUTrue_    = fs->make<TH1F>("hPUTrue",    "number of true pilepu", 1000, 0, 200);
  hGenWeight_ = fs->make<TH1F>("hGenWeight", "Gen weights",           2,    0, 2);


  //genlevel info about W/Z boson
  tree->Branch("bosonMass", &bosonMass);
  tree->Branch("bosonPt",&bosonPt);
  tree->Branch("neutrinos",&neutrinos);
  tree->Branch("BosonPID",        &BosonPID);
  tree->Branch("BosonVtx",        &BosonVtx);
  tree->Branch("BosonVty",        &BosonVty);
  tree->Branch("BosonVtz",        &BosonVtz);
  tree->Branch("BosonPt",         &BosonPt);
  tree->Branch("BosonMass",       &BosonMass);
  tree->Branch("BosonEta",        &BosonEta);
  tree->Branch("BosonPhi",        &BosonPhi);
  tree->Branch("BosonE",          &BosonE);
  tree->Branch("BosonEt",         &BosonEt);
  tree->Branch("BosonStatus",         &BosonStatus);

  tree->Branch("totalET", &totalET);
  tree->Branch("HT", &HT);
  tree->Branch("pfMET", &pfMET);
  tree->Branch("pfMETPhi", &pfMETPhi);
  tree->Branch("pfMETsumEt",       &pfMETsumEt_);
  tree->Branch("pfMETmEtSig",      &pfMETmEtSig_);
  tree->Branch("pfMETSig",         &pfMETSig_);
  tree->Branch("caloMET", &caloMET);

  tree->Branch("calojetPt", &calojetPt_);
  tree->Branch("calojetEn", &calojetEn_);
  tree->Branch("calojetEta",&calojetEta_);
  tree->Branch("calojetPhi",&calojetPhi_);
  tree->Branch("calojetEEF",&calojetEEF_);
  tree->Branch("calojetHEF",&calojetHEF_);
  tree->Branch("calojetTowersArea",&calojetTowersArea_);
  tree->Branch("calojetMaxEInEmTowers",&calojetMaxEInEmTowers_);
  tree->Branch("calojetMaxEInHadTowers",&calojetMaxEInHadTowers_);

  tree->Branch("nJets", &nJets);
  //tree->Branch("jet1Ptdiff", &jet1Ptdiff); //checking that JEC are applied correctly. It should be zero
  tree->Branch("jetPx",           &jetPx_);
  tree->Branch("jetPy",           &jetPy_);
  tree->Branch("jetPz",           &jetPz_);
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
  tree->Branch("jetNNP",       &jetNNP_);
  tree->Branch("jetMUF",       &jetMUF_);
  tree->Branch("jetHFHAE",         &jetHFHAE_);
  tree->Branch("jetHFEME",         &jetHFEME_);
  tree->Branch("jetNConstituents", &jetNConstituents_);
  tree->Branch("jetCSV2BJetTags", &jetCSV2BJetTags_);
  tree->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
  tree->Branch("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags_);
  tree->Branch("jetEtaWidth", &jetEtaWidth_);
  tree->Branch("jetPhiWidth", &jetPhiWidth_);
  tree->Branch("jetEtaWidthInECal", &jetEtaWidthInECal_);
  tree->Branch("jetEtaWidthInHCal", &jetEtaWidthInHCal_);
  tree->Branch("jetPhiWidthInECal", &jetPhiWidthInECal_);  
  tree->Branch("jetPhiWidthInHCal", &jetPhiWidthInHCal_);

  tree->Branch("jetPhotonEnergyFraction",&jetPhotonEnergyFraction_);
  tree->Branch("jetJetArea",&jetJetArea_);
  tree->Branch("jetMt",           &jetMt_);
  tree->Branch("jetMass",         &jetMass_);
  tree->Branch("jetMaxDistance",&jetMaxDistance_);
  tree->Branch("jetPhiPhiMoment",&jetPhiPhiMoment_);
  tree->Branch("jetConstituentEtaPhiSpread",&jetConstituentEtaPhiSpread_);
  tree->Branch("jetConstituentPtDistribution",&jetConstituentPtDistribution_);
  tree->Branch("jetPileup",&jetPileup_);
  tree->Branch("jetN60",&jetN60_);
  tree->Branch("jetN90",&jetN90_);

  tree->Branch("nMu",           &nMu_);
  tree->Branch("muPt",          &muPt_);
  tree->Branch("muEn",          &muEn_);
  tree->Branch("muEta",         &muEta_);
  tree->Branch("muPhi",         &muPhi_);
  tree->Branch("muCharge",      &muCharge_);
  tree->Branch("muType",        &muType_);
  tree->Branch("muIsLooseID",   &muIsLooseID_);
  tree->Branch("muIsMediumID",  &muIsMediumID_);
  tree->Branch("muIsTightID",   &muIsTightID_);
  tree->Branch("muIsSoftID",    &muIsSoftID_);
  tree->Branch("muIsHighPtID",  &muIsHighPtID_);
  tree->Branch("muD0",          &muD0_);
  tree->Branch("muDz",          &muDz_);
  tree->Branch("muChi2NDF",     &muChi2NDF_);
  tree->Branch("muInnerD0",     &muInnerD0_);
  tree->Branch("muInnerDz",     &muInnerDz_);
  tree->Branch("muTrkLayers",   &muTrkLayers_);
  tree->Branch("muPixelLayers", &muPixelLayers_);
  tree->Branch("muPixelHits",   &muPixelHits_);
  tree->Branch("muMuonHits",    &muMuonHits_);
  tree->Branch("muStations",    &muStations_);
  tree->Branch("muMatches",     &muMatches_);
  tree->Branch("muTrkQuality",  &muTrkQuality_);
  tree->Branch("muIsoTrk",      &muIsoTrk_);
  tree->Branch("muPFChIso",     &muPFChIso_);
  tree->Branch("muPFPhoIso",    &muPFPhoIso_);
  tree->Branch("muPFNeuIso",    &muPFNeuIso_);
  tree->Branch("muPFPUIso",     &muPFPUIso_);
  tree->Branch("muInnervalidFraction",   &muInnervalidFraction_);
  tree->Branch("musegmentCompatibility", &musegmentCompatibility_);
  tree->Branch("muchi2LocalPosition",    &muchi2LocalPosition_);
  tree->Branch("mutrkKink",              &mutrkKink_);
  tree->Branch("muBestTrkPtError",       &muBestTrkPtError_);
  tree->Branch("muBestTrkPt",            &muBestTrkPt_);


  tree->Branch("nEle",                    &nEle_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_); 
  tree->Branch("ElectronPassVetoID", &PassVeto_);
  tree->Branch("ElectronPassLooseID", &PassLoose_);
  tree->Branch("ElectronPassMediumID", &PassMedium_); 
  tree->Branch("ElectronPassTightID", &PassTight_);
  tree->Branch("ElectronPassHEEPID", &PassHEEP_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("eleR9",                   &eleR9_);
  tree->Branch("eleCalibPt",              &eleCalibPt_);
  tree->Branch("eleCalibEn",              &eleCalibEn_);
  tree->Branch("eleSCEta",                &eleSCEta_);
  tree->Branch("eleSCPhi",                &eleSCPhi_);
  tree->Branch("eleSCRawEn",              &eleSCRawEn_);
  tree->Branch("eleSCEtaWidth",           &eleSCEtaWidth_);
  tree->Branch("eleSCPhiWidth",           &eleSCPhiWidth_);
  tree->Branch("eleHoverE",               &eleHoverE_);
  tree->Branch("eleEoverP",               &eleEoverP_);
  tree->Branch("eleEoverPout",            &eleEoverPout_);
  tree->Branch("eleEoverPInv",            &eleEoverPInv_);
  tree->Branch("eleBrem",                 &eleBrem_);
  tree->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  tree->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  tree->Branch("eledEtaAtCalo",           &eledEtaAtCalo_);
  tree->Branch("eleSigmaIEtaIEta",        &eleSigmaIEtaIEta_);
  tree->Branch("eleSigmaIPhiIPhi",        &eleSigmaIPhiIPhi_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("eleConvVeto",             &eleConvVeto_);
  tree->Branch("eleMissHits",             &eleMissHits_);
  tree->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  tree->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);
  tree->Branch("elePFMiniIso",            &elePFMiniIso_);
  tree->Branch("eleIDMVANonTrg",          &eleIDMVANonTrg_);
  tree->Branch("eleIDMVATrg",             &eleIDMVATrg_);
  tree->Branch("eledEtaseedAtVtx",        &eledEtaseedAtVtx_);
  tree->Branch("eleE1x5",                 &eleE1x5_);
  tree->Branch("eleE2x5",                 &eleE2x5_);
  tree->Branch("eleE5x5",                 &eleE5x5_);
  tree->Branch("eleE1x5Full5x5",          &eleE1x5Full5x5_);
  tree->Branch("eleE2x5Full5x5",          &eleE2x5Full5x5_);
  tree->Branch("eleE5x5Full5x5",          &eleE5x5Full5x5_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleisoChargedHadrons"      , &eleisoChargedHadrons_);
  tree->Branch("eleisoNeutralHadrons"      , &eleisoNeutralHadrons_);
  tree->Branch("eleisoPhotons"             , &eleisoPhotons_);
  tree->Branch("eleisoChargedFromPU"       , &eleisoChargedFromPU_);

  tree->Branch("nPhotons",                    &nPhotons_);
  tree->Branch("phopx",&phopx_);
  tree->Branch("phopy",&phopy_);
  tree->Branch("phopz",&phopz_);
  tree->Branch("phoE",&phoE_);
  tree->Branch("phopt"  ,  &phopt_    );
  tree->Branch("phoeta" ,  &phoeta_ );
  tree->Branch("phophi" ,  &phophi_ );
  tree->Branch("phor9",&phor9_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5"   , &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phohOverE"                 ,  &phohOverE_);
  tree->Branch("phohasPixelSeed"           ,  &phohasPixelSeed_);

  tree->Branch("phoisoChargedHadrons"      , &phoisoChargedHadrons_);
  tree->Branch("phoWorstisoChargedHadrons" , &phoWorstisoChargedHadrons_);
  tree->Branch("phoisoNeutralHadrons"      , &phoisoNeutralHadrons_);
  tree->Branch("phoisoPhotons"             , &phoisoPhotons_);


  tree->Branch("trackPt",&trkPt_);
  tree->Branch("trackEta",&trkEta_);
  tree->Branch("trackPhi",&trkPhi_);
  //Two Leading Tracks Associated with the leading Jet
  //
  tree->Branch("j1nTracks",&j1nTracks);
  tree->Branch("j1trk1Quality",&j1trk1Quality);
  tree->Branch("j1trk1Charge",&j1trk1Charge);
  tree->Branch("j1dRCand1Trk1",&j1dRCand1Trk1);
  tree->Branch("j1Trk1PFCand1Match",&j1Trk1PFCand1Match);
  tree->Branch("j1trk2Quality",&j1trk2Quality);
  tree->Branch("j1trk2Charge",&j1trk2Charge);
  tree->Branch("j1dRCand2Trk2",&j1dRCand2Trk2);
  tree->Branch("j1Trk2PFCand2Match",&j1Trk2PFCand2Match);
  tree->Branch("j1trk1PtError",&j1trk1PtError);
  tree->Branch("j1trk2PtError",&j1trk2PtError);
  tree->Branch("j1ptTrk1Ratio",&j1ptTrk1Ratio);
  tree->Branch("j1ptTrk2Ratio",&j1ptTrk2Ratio);
  tree->Branch("j1trk12PT", &j1trk12PT);
  tree->Branch("j1trk12Ratio",&j1trk12Ratio);
  tree->Branch("j1trk1PT", &j1trk1PT);
  tree->Branch("j1trk1Eta", &j1trk1Eta);
  tree->Branch("j1trk1Phi", &j1trk1Phi);
  tree->Branch("j1trk2PT", &j1trk2PT);
  tree->Branch("j1trk2Eta", &j1trk2Eta);
  tree->Branch("j1trk2Phi", &j1trk2Phi);

  //Two leading Tracks associated with the second leading Jet
  tree->Branch("j2nTracks",&j2nTracks);
  tree->Branch("j2trk1Quality",&j2trk1Quality);
  tree->Branch("j2trk1Charge",&j2trk1Charge);
  tree->Branch("j2dRCand1Trk1",&j2dRCand1Trk1);
  tree->Branch("j2Trk1PFCand1Match",&j2Trk1PFCand1Match);
  tree->Branch("j2trk2Quality",&j2trk2Quality);
  tree->Branch("j2trk2Charge",&j2trk2Charge);
  tree->Branch("j2dRCand2Trk2",&j2dRCand2Trk2);
  tree->Branch("j2Trk2PFCand2Match",&j2Trk2PFCand2Match);
  tree->Branch("j2trk1PtError",&j2trk1PtError);
  tree->Branch("j2trk2PtError",&j2trk2PtError); 
  tree->Branch("j2ptTrk1Ratio",&j2ptTrk1Ratio);
  tree->Branch("j2ptTrk2Ratio",&j2ptTrk2Ratio);
  tree->Branch("j2trk12PT", &j2trk12PT);
  tree->Branch("j2trk12Ratio",&j2trk12Ratio);
  tree->Branch("j2trk1PT", &j2trk1PT);
  tree->Branch("j2trk1Eta", &j2trk1Eta);
  tree->Branch("j2trk1Phi", &j2trk1Phi);
  tree->Branch("j2trk2PT", &j2trk2PT);
  tree->Branch("j2trk2Eta", &j2trk2Eta);
  tree->Branch("j2trk2Phi", &j2trk2Phi);

  tree->Branch("nGoodJets", &nGoodJets);
  tree->Branch("j1nCarrying60",&j1nCarrying60);
  tree->Branch("j1nCarrying90",&j1nCarrying90); //nConstituents of Leading Jet carrying 90% of jet energy
  tree->Branch("j1nCarrying92",&j1nCarrying92); //nConstituents of Leading Jet carrying 92%
  tree->Branch("j1nCarrying94",&j1nCarrying94); //nConstituents of Leading Jet carrying 94%
  tree->Branch("j1nCarrying96",&j1nCarrying96); //nConstituents of Leading Jet carrying 96%
  tree->Branch("j1nCarrying98",&j1nCarrying98); //nConstituents of Leading Jet carrying 98%  
  tree->Branch("j1PT", &j1PT);
  tree->Branch("j1Eta", &j1Eta);
  tree->Branch("j1Phi", &j1Phi);
  tree->Branch("j1CHdFr", &j1CHdFr);
  tree->Branch("j1NHdFr", &j1NHdFr);
  tree->Branch("j1CEmFr", &j1CEmFr);
  tree->Branch("j1NEmFr", &j1NEmFr);
  tree->Branch("j1PhoEFr", &j1PhoEFr);
  tree->Branch("j1EleEFr", &j1EleEFr);
  tree->Branch("j1MuEFr", &j1MuEFr);
  tree->Branch("j1CMuEFr", &j1CMuEFr);
  tree->Branch("j1nCons", &j1nCons);
  tree->Branch("j1CMty", &j1CMty);
  tree->Branch("j1NMty", &j1NMty);
  tree->Branch("j1CHdMty", &j1CHdMty);
  tree->Branch("j1NHdMty", &j1NHdMty);
  tree->Branch("j1PhoMty", &j1PhoMty);
  tree->Branch("j1EleMty", &j1EleMty);
  tree->Branch("j1MuMty", &j1MuMty);
  tree->Branch("j1etaWidth", &j1etaWidth);
  tree->Branch("j1phiWidth", &j1phiWidth);
  tree->Branch("j1etaWidthInECal", &j1etaWidthInECal);
  tree->Branch("j1phiWidthInECal", &j1phiWidthInECal);
  tree->Branch("j1etaWidthInHCal", &j1etaWidthInHCal);
  tree->Branch("j1phiWidthInHCal", &j1phiWidthInHCal);

  //PFCandidate Info associated with leading Jet 
  tree->Branch("j1PFCand12PtSum", &j1PFCand12PtSum);
  tree->Branch("j1PFCandsPtSum", &j1PFCandsPtSum);   
  tree->Branch("j1PFCand12Ratio",&j1PFCand12Ratio);

  tree->Branch("j2nCarrying60",&j2nCarrying60);
  tree->Branch("j2nCarrying90",&j2nCarrying90); //nConstituents of Leading Jet carrying 90% of jet energy
  tree->Branch("j2nCarrying92",&j2nCarrying92); //nConstituents of Leading Jet carrying 92%
  tree->Branch("j2nCarrying94",&j2nCarrying94); //nConstituents of Leading Jet carrying 94%
  tree->Branch("j2nCarrying96",&j2nCarrying96); //nConstituents of Leading Jet carrying 96%
  tree->Branch("j2nCarrying98",&j2nCarrying98); //nConstituents of Leading Jet carrying 98%  
  tree->Branch("j2PT", &j2PT);
  tree->Branch("j2Eta", &j2Eta);
  tree->Branch("j2Phi", &j2Phi);
  tree->Branch("j2CHdFr", &j2CHdFr);
  tree->Branch("j2NHdFr", &j2NHdFr);
  tree->Branch("j2CEmFr", &j2CEmFr);
  tree->Branch("j2NEmFr", &j2NEmFr);
  tree->Branch("j2PhoEFr", &j2PhoEFr);
  tree->Branch("j2EleEFr", &j2EleEFr);
  tree->Branch("j2MuEFr", &j2MuEFr);
  tree->Branch("j2CMuEFr", &j2CMuEFr);
  tree->Branch("j2nCons", &j2nCons);
  tree->Branch("j2CMty", &j2CMty);
  tree->Branch("j2NMty", &j2NMty);
  tree->Branch("j2CHdMty", &j2CHdMty);
  tree->Branch("j2NHdMty", &j2NHdMty);
  tree->Branch("j2PhoMty", &j2PhoMty);
  tree->Branch("j2EleMty", &j2EleMty);
  tree->Branch("j2MuMty", &j2MuMty);
  tree->Branch("j2etaWidth", &j2etaWidth);
  tree->Branch("j2phiWidth", &j2phiWidth);
  tree->Branch("j2etaWidthInECal", &j2etaWidthInECal);
  tree->Branch("j2phiWidthInECal", &j2phiWidthInECal);
  tree->Branch("j2etaWidthInHCal", &j2etaWidthInHCal);
  tree->Branch("j2phiWidthInHCal", &j2phiWidthInHCal);
  //PFCandidate Info associated with Second leading Jet
  tree->Branch("j2PFCand12PtSum", &j2PFCand12PtSum);
  tree->Branch("j2PFCandsPtSum", &j2PFCandsPtSum);
  tree->Branch("j2PFCand12Ratio",&j2PFCand12Ratio);

  tree->Branch("HLTPFMET300",               &HLTMET300_);
  tree->Branch("HLTPFMET170_HBHECleaned", &HLTMET170_HBHE);
  tree->Branch("HLTPhoton165_HE10",               &HLTPhoton165_HE10_);
  tree->Branch("HLTPhoton175",               &HLTPhoton175_);
  tree->Branch("HLTPhoton75",               &HLTPhoton75_);
  tree->Branch("HLTPhoton90",               &HLTPhoton90_);
  tree->Branch("HLTPhoton120",               &HLTPhoton120_);
  tree->Branch("HLTPFJet40",               &HLTPFJet40_);
  tree->Branch("HLTPFJet60",               &HLTPFJet60_);
  tree->Branch("HLTPFJet80",               &HLTPFJet80_);
}


JetAnalyzer::~JetAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
  void
JetAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  metFilters_ = 0;
  if (iEvent.isRealData()){

    Handle<bool> globalHandle;
    iEvent.getByToken(globalHandle_, globalHandle);
    bool GlobalHaloResult_ = *globalHandle;

    Handle<bool> hcalNoiseHandle;
    iEvent.getByToken(hcalNoiseHandle_, hcalNoiseHandle);
    bool HBHENoiseResult_ = *hcalNoiseHandle;

    Handle<bool> hcalIsoNoiseHandle;
    iEvent.getByToken(hcalIsoNoiseHandle_, hcalIsoNoiseHandle);
    bool HBHEIsoNoiseResult_ = *hcalIsoNoiseHandle;

    Handle<bool> eCALTPHandle;
    iEvent.getByToken(eCALTPHandle_, eCALTPHandle);
    bool EcalDeadCellTFResult_ = *eCALTPHandle;

    Handle<bool> bADSCHandle;
    iEvent.getByToken(bADSCHandle_, bADSCHandle);
    bool EEBadSCResult_ = *bADSCHandle;

    Handle<bool> ifilterbadChCand;
    iEvent.getByToken(BadChCandFilterToken_, ifilterbadChCand);
    bool filterbadChCandidate = *ifilterbadChCand;

    Handle<bool> ifilterbadPFMuon;
    iEvent.getByToken(BadPFMuonFilterToken_, ifilterbadPFMuon);
    bool filterbadPFMuon = *ifilterbadPFMuon;

    // std::cout<<"event: "<<event_<<std::endl;
    // std::cout<<"metFilters: "<<metFilters_<<std::endl;

    if ( !HBHENoiseResult_      ) metFilters_ += 1;
    if ( !HBHEIsoNoiseResult_   ) metFilters_ += 2; 
    if ( !GlobalHaloResult_        ) metFilters_ += 4; 
    if ( !EEBadSCResult_        ) metFilters_ += 8;
    if ( !EcalDeadCellTFResult_ ) metFilters_ += 16;
    if ( !filterbadChCandidate) metFilters_ += 32;
    //std::cout<<"metFilters(badChCandidateFilter): "<<metFilters_<<std::endl;
    if ( !filterbadPFMuon) metFilters_ += 64;
    //std::cout<<"metFilters(badPFMuonFilter): "<<metFilters_<<std::endl;   
  }

  Handle< vector<reco::PFCandidate> > pfCands;
  iEvent.getByToken(pfCandsToken, pfCands);

  Handle< vector<reco::PFJet> > pfJets;
  iEvent.getByToken(pfJetsToken, pfJets);

  Handle< vector<reco::CaloJet> > caloJets;
  iEvent.getByToken(caloJetsToken, caloJets);

  // Handle< vector<pat::Jet> > Jets;
  // iEvent.getByToken(JetsToken, Jets);

  Handle< vector<reco::PFMET> > pfMETs;
  iEvent.getByToken(pfMETsToken, pfMETs);


  Handle< vector<reco::CaloMET> > caloMETs;
  iEvent.getByToken(caloMETToken, caloMETs);

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();

  trkPt_                      .clear();
  trkEta_                     .clear();
  trkPhi_                     .clear();

  Handle< vector<reco::Track> > recoTracks;
  iEvent.getByToken(trackToken,recoTracks);

  for(uint32_t i = 0; i < recoTracks->size(); i++) {
    const reco::Track& track = (*recoTracks)[i];
    trkPt_    .push_back(track.pt());
    trkEta_   .push_back(track.eta());
    trkPhi_   .push_back(track.phi());
  }   

  Handle< vector<reco::Vertex> > vtxHandle;
  iEvent.getByToken(vtxToken_, vtxHandle);

  //number of primary vertices 

  if (vtxHandle->empty()) return;
  const reco::Vertex &pv = vtxHandle->front();
  npv_ = vtxHandle->size();    
  nTrksPV_ = pv.nTracks();




  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionsToken_, conversions);

  // Get rho value
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  rho_ = *rhoH;

  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  

  //gen level variables cleared from last execution
  nup_ = 0;
  numGenJets_ = 0; //number of outgoing partons
  genHT_ = 0;
  nPUInfo_ =0;
  nPU_          .clear();
  puBX_         .clear();
  puTrue_       .clear();
  BosonPID       .clear();
  BosonVtx       .clear();
  BosonVty       .clear();
  BosonVtz       .clear();
  BosonPt        .clear();
  BosonMass      .clear();
  BosonEta       .clear();
  BosonPhi       .clear();
  BosonE         .clear();
  BosonEt        .clear();
  BosonStatus    .clear();
  pdf_          .clear();

  mcPID       .clear();
  mcVtx       .clear();
  mcVty       .clear();
  mcVtz       .clear();
  mcPt        .clear();
  mcMass      .clear();
  mcEta       .clear();
  mcPhi       .clear();
  mcE         .clear();
  mcEt        .clear();
  mcIndex     .clear();
  mcStatusFlag.clear();
  mcStatus    .clear();

  nMC_ = 0;

  pthat_     = -99;
  processID_ = -99;
  genWeight_ = -99;


  Handle<GenEventInfoProduct> genEventInfoHandle;
  iEvent.getByToken(generatorLabel_, genEventInfoHandle);

  if (genEventInfoHandle.isValid()) {

    if (genEventInfoHandle->pdf()) {
      pdf_.push_back(genEventInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->x.first);     // x value of parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->x.second);    // x value of parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->scalePDF);    // scale of the hard interaction
    }

    if (genEventInfoHandle->hasBinningValues())
      pthat_ = genEventInfoHandle->binningValues()[0];
    processID_ = genEventInfoHandle->signalProcessID();
    genWeight_ = genEventInfoHandle->weight();
    //std::cout<<"Gen weight: "<<genWeight_<<std::endl;
    if (genWeight_ >= 0) hGenWeight_->Fill(0.5);    
    else hGenWeight_->Fill(1.5);
  } else
    edm::LogWarning("ggNtuplizer") << "no GenEventInfoProduct in event";


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
  }
  Handle<vector<PileupSummaryInfo> > genPileupHandle;
  iEvent.getByToken(puCollection_, genPileupHandle);

  if (genPileupHandle.isValid()) {
    for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
	hPU_->Fill(pu->getPU_NumInteractions());
	hPUTrue_->Fill(pu->getTrueNumInteractions());
      }

      nPU_   .push_back(pu->getPU_NumInteractions());
      puTrue_.push_back(pu->getTrueNumInteractions());
      puBX_  .push_back(pu->getBunchCrossing());

      nPUInfo_++;
    }
  }
  else
    edm::LogWarning("JetAnalyzer") << "no PileupSummaryInfo in event";

  Handle<vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  reco::Candidate::LorentzVector p4(0,0,0,0);
  if (genParticlesHandle.isValid()) {
    int genIndex = 0;
    neutrinos = 0; //count the no.of neutrinos in an event.

    for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
      genIndex++;
      int pdg = abs(ip->pdgId());
      
      const reco::Candidate *p = (const reco::Candidate*)&(*ip);
      if (!p->mother()) continue;

      mcPID    .push_back(p->pdgId());
      mcVtx    .push_back(p->vx());
      mcVty    .push_back(p->vy());
      mcVtz    .push_back(p->vz());
      mcPt     .push_back(p->pt());
      mcMass   .push_back(p->mass());
      mcEta    .push_back(p->eta());
      mcPhi    .push_back(p->phi());
      mcE      .push_back(p->energy());
      mcEt     .push_back(p->et());
      mcStatus .push_back(p->status());

      UShort_t tmpStatusFlag = 0;
      if (ip->fromHardProcessFinalState()) setbit(tmpStatusFlag, 0);
      if (ip->isPromptFinalState())        setbit(tmpStatusFlag, 1);
      if (ip->isHardProcess())  setbit(tmpStatusFlag, 2);

      // if genParticle is W or Z, check its decay type
      if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
	for (size_t k=0; k < p->numberOfDaughters(); ++k) {
	  const reco::Candidate *dp = p->daughter(k);
	  if (abs(dp->pdgId())<=6)                               setbit(tmpStatusFlag, 4);
	  else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12) setbit(tmpStatusFlag, 5);
	  else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14) setbit(tmpStatusFlag, 6);
	  else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16) setbit(tmpStatusFlag, 7);
	}
      } 
      mcStatusFlag.push_back(tmpStatusFlag);
      mcIndex.push_back(genIndex-1);

      if ((ip->isHardProcess()) &&  (pdg==11 || pdg==13 || pdg==15 || pdg==12||pdg==14 ||pdg==16)   ){
	p4+= ip->p4();
      }
      if((ip->isHardProcess()) && (pdg==12||pdg==14 ||pdg==16)){
	neutrinos++;
      }

      //select Z/W boson
      bool boson = ((    ip->pdgId()  == 23 && ip->isHardProcess()) || 
	  (abs(ip->pdgId()) == 24 && ip->isHardProcess()));

      if (boson) {

	BosonPID    .push_back(ip->pdgId());
	BosonVtx    .push_back(ip->vx());
	BosonVty    .push_back(ip->vy());
	BosonVtz    .push_back(ip->vz());
	BosonPt     .push_back(ip->pt());
	BosonMass   .push_back(ip->mass());
	BosonEta    .push_back(ip->eta());
	BosonPhi    .push_back(ip->phi());
	BosonE      .push_back(ip->energy());
	BosonEt     .push_back(ip->et());
	BosonStatus .push_back(ip->status());    
      }
    }
    bosonMass = p4.M();
    bosonPt = p4.Pt(); 
  }

  else
    edm::LogWarning("JetAnalyzer") << "no reco::GenParticles in event";
  //   if (vtxHandle.isValid()) {
  //   nVtx_ = 0;
  //   for (uint32_t v = 0; v < vtxHandle->size(); v++) {
  //   const reco::Vertex& vertex = (*vtxHandle)[v];
  //   if (v == 0) {
  //       vtx_0 = vertex;
  //   } 
  //   if (nVtx_ == 0) {
  //       nTrksPV_ = vertex.nTracks();
  //    }
  //   }
  //   }
  //  

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

  nPhotons_=0;

  edm::Handle<edm::View<reco::Photon> > photons;
  iEvent.getByToken(photonsToken_, photons);
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  iEvent.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;
  iEvent.getByToken(phoWorstChargedIsolationToken_, phoWorstChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  iEvent.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  iEvent.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);


  nMu_ = 0;


  edm::Handle<edm::View<pat::Muon> > recoMuons;
  //edm::Handle<vector<reco::Muon> > recoMuons;
  iEvent.getByToken(muonToken, recoMuons);



  //   Handle< vector<pat::Muon> > recoMuons;
  //   iEvent.getByToken(muonToken, recoMuons);

  //HLT treatment
  HLTMET300_               = -99;
  HLTMET170_HBHE =-99;
  HLTPhoton165_HE10_ =-99;
  HLTPhoton175_ =-99;
  HLTPhoton75_ =-99;
  HLTPhoton90_ =-99;
  HLTPhoton120_ =-99;
  HLTPFJet40_=-99;
  HLTPFJet60_=-99;
  HLTPFJet80_=-99;

  Handle<edm::TriggerResults> trgResultsHandle;
  iEvent.getByToken(trgResultsLabel_, trgResultsHandle);


  const edm::TriggerNames &trgNames = iEvent.triggerNames(*trgResultsHandle);

  for (size_t i = 0; i < trgNames.size(); ++i) {
    const string &name = trgNames.triggerName(i);
    //Jet triggers

    if (name.find("HLT_PFMET300_v") != string::npos) {HLTMET300_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_PFMET170_HBHECleaned_v") != string::npos) {HLTMET170_HBHE = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_Photon165_HE10_v") != string::npos) {HLTPhoton165_HE10_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_Photon175_v") != string::npos) {HLTPhoton175_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_Photon75_v") != string::npos) {HLTPhoton75_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_Photon90_v") != string::npos) {HLTPhoton90_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_Photon120_v") != string::npos) {HLTPhoton120_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_PFJet40_v") != string::npos) {HLTPFJet40_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_PFJet60_v") != string::npos) {HLTPFJet60_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
    if (name.find("HLT_PFJet80_v") != string::npos) {HLTPFJet80_ = (trgResultsHandle->accept(i)) ? 1 : 0;}
  }
  //Clear previous events



  calojetPt_ .clear();
  calojetEn_ .clear();
  calojetEta_.clear();
  calojetPhi_.clear();
  calojetEEF_.clear();
  calojetHEF_.clear();
  calojetTowersArea_.clear();
  calojetMaxEInEmTowers_.clear();
  calojetMaxEInHadTowers_.clear();


  jetPx_                                  .clear();
  jetPy_                                  .clear();
  jetPz_                                  .clear();
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
  jetNNP_                                 .clear();
  jetMUF_                                 .clear();
  jetHFHAE_                               .clear();
  jetHFEME_                               .clear();
  jetNConstituents_                       .clear();
  jetCSV2BJetTags_                        .clear();
  jetJetProbabilityBJetTags_              .clear();
  jetpfCombinedMVAV2BJetTags_             .clear();
  jetEtaWidth_.clear();
  jetPhiWidth_.clear();   
  jetEtaWidthInECal_.clear();
  jetPhiWidthInECal_.clear();
  jetEtaWidthInHCal_.clear();
  jetPhiWidthInHCal_.clear(); 

  jetPhotonEnergyFraction_.clear();
  jetJetArea_.clear();
  jetMt_                                  .clear();
  jetMass_                                .clear();
  jetMaxDistance_.clear();
  jetPhiPhiMoment_.clear();
  jetConstituentEtaPhiSpread_.clear();
  jetConstituentPtDistribution_.clear();
  jetPileup_.clear();
  jetN60_.clear();
  jetN90_.clear();

  
  PassVeto_.clear();
  PassLoose_.clear();
  PassMedium_.clear();
  PassTight_.clear();
  PassHEEP_.clear();
  eleCharge_                  .clear();
  eleChargeConsistent_        .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleD0_                      .clear();
  eleDz_                      .clear();
  elePt_                      .clear();
  eleEta_                     .clear();
  elePhi_                     .clear();
  eleR9_                      .clear();
  eleCalibPt_                 .clear();
  eleCalibEn_                 .clear();
  eleSCEta_                   .clear();
  eleSCPhi_                   .clear();
  eleSCRawEn_                 .clear();
  eleSCEtaWidth_              .clear();
  eleSCPhiWidth_              .clear();
  eleHoverE_                  .clear();
  eleEoverP_                  .clear();
  eleEoverPout_               .clear();
  eleEoverPInv_               .clear();
  eleBrem_                    .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eledEtaAtCalo_              .clear();
  eleSigmaIEtaIEta_           .clear();
  eleSigmaIPhiIPhi_           .clear();
  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  eleConvVeto_                .clear();
  eleMissHits_                .clear();
  elePFClusEcalIso_           .clear();
  elePFClusHcalIso_           .clear();
  elePFMiniIso_               .clear();
  eleIDMVANonTrg_             .clear();
  eleIDMVATrg_                .clear();
  eledEtaseedAtVtx_           .clear();
  eleE1x5_                    .clear();
  eleE2x5_                    .clear();
  eleE5x5_                    .clear();
  eleE1x5Full5x5_             .clear();
  eleE2x5Full5x5_             .clear();
  eleE5x5Full5x5_             .clear();
  eleR9Full5x5_               .clear();
  eleisoChargedHadrons_.clear();
  eleisoNeutralHadrons_.clear();
  eleisoPhotons_.clear();
  eleisoChargedFromPU_.clear();


  muPt_         .clear();
  muEn_         .clear();
  muEta_        .clear();
  muPhi_        .clear();
  muCharge_     .clear();
  muType_       .clear();
  muPt_         .clear();
  muEn_         .clear();
  muEta_        .clear();
  muPhi_        .clear();
  muCharge_     .clear();
  muType_       .clear();
  muIsLooseID_  .clear();
  muIsMediumID_ .clear();
  muIsTightID_  .clear();
  muIsSoftID_   .clear();
  muIsHighPtID_ .clear();
  muD0_         .clear();
  muDz_         .clear();
  muChi2NDF_    .clear();
  muInnerD0_    .clear();
  muInnerDz_    .clear();
  muTrkLayers_  .clear();
  muPixelLayers_.clear();
  muPixelHits_  .clear();
  muMuonHits_   .clear();
  muStations_   .clear();
  muMatches_    .clear();
  muTrkQuality_ .clear();
  muIsoTrk_     .clear();
  muPFChIso_    .clear();
  muPFPhoIso_   .clear();
  muPFNeuIso_   .clear();
  muPFPUIso_    .clear();
  muInnervalidFraction_  .clear();
  musegmentCompatibility_.clear();
  muchi2LocalPosition_   .clear();
  mutrkKink_             .clear();
  muBestTrkPtError_      .clear();
  muBestTrkPt_           .clear();



  // Set event level quantities

  phopt_.clear();
  phoeta_.clear();
  phophi_.clear();
  //
  phoSigmaIEtaIEtaFull5x5_.clear();
  phohOverE_.clear();
  phohasPixelSeed_.clear();
  //
  phoisoChargedHadrons_.clear();
  phoWorstisoChargedHadrons_.clear();
  phoisoNeutralHadrons_.clear();
  phoisoPhotons_.clear();

  phopx_.clear();
  phopy_.clear();
  phopz_.clear();
  phoE_ .clear();
  phor9_.clear();


  totalET = 0;
  HT = 0;
  nJets = 0;
  nGoodJets = 0;
  j1Trk1PFCand1Match = 0;
  j1Trk2PFCand2Match = 0;
  j2Trk1PFCand1Match = 0;
  j2Trk2PFCand2Match = 0;

  for(uint32_t j = 0; j < caloJets->size(); j++) {
    const reco::CaloJet &jet = (*caloJets)[j];
    calojetPt_.push_back( jet.pt());
    calojetEn_.push_back( jet.energy());
    calojetEta_.push_back(jet.eta());
    calojetPhi_.push_back(jet.phi());
    calojetEEF_.push_back(   jet.emEnergyFraction());
    calojetHEF_.push_back(   jet.energyFractionHadronic());
    calojetTowersArea_.push_back(     jet.towersArea());
    calojetMaxEInEmTowers_.push_back( jet.maxEInEmTowers());
    calojetMaxEInHadTowers_.push_back(jet.maxEInHadTowers());
  }
  //jet1Ptdiff = 0;
  Handle< vector<pat::Jet> > Jets;
  iEvent.getByToken(JetsToken, Jets);
  for(uint32_t J = 0; J < Jets->size(); J++) {
    const pat::Jet &j = (*Jets)[J];
    //if(j.pt() < 20) continue;
    //b-tagging
    jetCSV2BJetTags_           .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetJetProbabilityBJetTags_ .push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
    jetpfCombinedMVAV2BJetTags_.push_back(j.bDiscriminator("pfCombinedMVAV2BJetTags"));
  }
  //std::cout<<"pat jetPt[0]: "<<(*Jets)[0].pt()<<std::endl;  
  //std::cout<<"pat Uncorrected jetPt[0]: "<<(*Jets)[0].correctedJet("Uncorrected").pt()<<std::endl;

  //test the JEC by filling the following variable
  //jet1Ptdiff = ((*Jets)[0].pt() - (*pfJets)[0].pt());
  // Compute HT
  //std::cout<<run_<<":"<<lumis_<<":"<<event_<<std::endl;
  for(uint32_t j = 0; j < pfJets->size(); j++) {
    const reco::PFJet &jet = (*pfJets)[j];
    //std::cout<<"reco jetPt[0]: "<<(*pfJets)[0].pt()<<std::endl;
    if(j == 0) {
      //std::cout<<"j1PT: " <<jet.pt()<<std::endl;
      j1PT = jet.pt();
      j1Eta = jet.eta();
      j1Phi = jet.phi();
      j1CHdFr = jet.chargedHadronEnergyFraction();
      j1NHdFr = jet.neutralHadronEnergyFraction();
      j1CEmFr = jet.chargedEmEnergyFraction();
      j1NEmFr = jet.neutralEmEnergyFraction();
      j1PhoEFr = jet.photonEnergyFraction();
      j1EleEFr = jet.electronEnergyFraction();
      j1MuEFr = jet.muonEnergyFraction();
      j1CMuEFr = jet.chargedMuEnergyFraction();
      j1nCons = jet.nConstituents();
      j1nCarrying60 = jet.nCarrying(0.60);
      j1nCarrying90 = jet.nCarrying(0.90); //return # of consitutents carrying 90% of energy
      j1nCarrying92 = jet.nCarrying(0.92);
      j1nCarrying94 = jet.nCarrying(0.94);
      j1nCarrying96 = jet.nCarrying(0.96);
      j1nCarrying98 = jet.nCarrying(0.98);
      //std::cout<<"Info starts about Leading Jet"<<std::endl;
      //get all tracks in the jets. All PFCandidates hold a reference to a track.
      const reco::TrackRefVector j1tracks = jet.getTrackRefs();
      std::vector<reco::TrackRef> j1tracksRef;//make a copy of RefVector to sort it
      j1nTracks = j1tracks.size();
      if(j1nTracks >0){
	//std::cout<<"Number of tracks(Leading Jet): "<< j1nTracks<<std::endl;
	for (const auto &trkRef : j1tracks){
	  reco::Track trk = *trkRef;
	  j1tracksRef.push_back(trkRef);
	}
	const std::vector<reco::PFCandidatePtr> pfCands = jet.getPFConstituents(); //get pfCands to match with tracks 
	std::sort(j1tracksRef.begin(),j1tracksRef.end(),[](const reco::TrackRef& track1, const reco::TrackRef& track2)
	    {return track1->pt() > track2->pt();}); //sorting the copy of the TrackRefVector by pT
	//std::cout<<"TrackRef Copy begin"<<std::endl;
	if(j1tracksRef.size()>1){
	  j1trk1Quality = j1tracksRef.at(0)->quality(reco::TrackBase::TrackQuality::highPurity);//determining the quality of track
	  j1trk1Charge = j1tracksRef.at(0)->charge();
	  j1trk1PT = j1tracksRef.at(0)->pt();
	  if(j1trk1Charge!=0){
	    j1trk1PtError=j1tracksRef.at(0)->ptError(); //error on Pt (set to 1000 TeV if charge==0 for safety) 
	    j1ptTrk1Ratio = (j1trk1PtError/j1trk1PT);}
	  j1trk1Eta = j1tracksRef.at(0)->eta();
	  j1trk1Phi = j1tracksRef.at(0)->phi();
	  //std::cout<<"Leading track Pt: "<< j1trk1PT<<std::endl;
	  j1trk2Quality = j1tracksRef.at(1)->quality(reco::TrackBase::TrackQuality::highPurity);
	  //std::cout<<"DPt/Pt: "<<j1ptTrk1Ratio<<std::endl;
	  //matching leading track candidate with PF-Candidate
	  j1dRCand1Trk1 = deltaR((*j1tracksRef.at(0)),(*pfCands[0]));
	  if (j1dRCand1Trk1 < 0.1){j1Trk1PFCand1Match = 1;}
	  j1trk2PT = j1tracksRef.at(1)->pt();
	  j1trk2Charge = j1tracksRef.at(1)->charge();
	  if(j1trk2Charge!=0){
	    j1trk2PtError = j1tracksRef.at(1)->ptError();//ptError of second leading track
	    j1ptTrk2Ratio = (j1trk2PtError/j1trk2PT);}
	  j1trk2Eta = j1tracksRef.at(1)->eta();
	  j1trk2Phi = j1tracksRef.at(1)->phi();
	  // std::cout<<"Second leading track Pt: "<< j1trk2PT<<std::endl;
	  //std::cout<<"DPt/Pt: "<<j1ptTrk2Ratio<<std::endl;
	  j1dRCand2Trk2 = deltaR((*j1tracksRef.at(1)),(*pfCands[1]));
	  if (j1dRCand2Trk2 < 0.1){j1Trk2PFCand2Match = 1;}         
	  j1trk12PT = j1trk1PT+j1trk2PT;//access pt of the first and second track in list
	}
	else{
	  j1trk1Quality = j1tracksRef.at(0)->quality(reco::TrackBase::TrackQuality::highPurity);
	  j1trk1Charge = j1tracksRef.at(0)->charge();
	  j1trk1PT = j1tracksRef.at(0)->pt();
	  if(j1trk1Charge!=0){
	    j1trk1PtError = j1tracksRef.at(0)->ptError(); //error on Pt (set to 1000 TeV if charge==0 for safety)
	    j1ptTrk1Ratio = (j1trk1PtError/j1trk1PT);}
	  //std::cout<<"DPt/Pt: "<<j1ptTrk1Ratio<<std::endl;
	  j1dRCand1Trk1 = deltaR((*j1tracksRef.at(0)),(*pfCands[0]));
	  if (j1dRCand1Trk1 < 0.1){j1Trk1PFCand1Match = 1;}
	  j1trk1Eta = j1tracksRef.at(0)->eta();
	  j1trk1Phi = j1tracksRef.at(0)->phi();
	  j1trk12PT = j1trk1PT;
	}
	//std::cout<<"j1trk1PtError: "<<j1trk1PtError<<std::endl;
	//std::cout<<"SumPT of two leading Tracks: " << j1trk12PT << std::endl;
	//std::cout<<"PT of the Leading Jet: " << j1PT << std::endl;
	j1trk12Ratio = (j1trk12PT/j1PT);
	//std::cout<<"Trk12PT/j1PT: "<< j1trk12Ratio <<std::endl;
      }
      j1CMty = jet.chargedMultiplicity();
      j1NMty = jet.neutralMultiplicity();
      j1CHdMty = jet.chargedHadronMultiplicity();
      j1NHdMty = jet.neutralHadronMultiplicity();
      j1PhoMty = jet.photonMultiplicity();
      j1EleMty = jet.electronMultiplicity();
      j1MuMty = jet.muonMultiplicity();
      JetWidthCalculator jwc(jet);
      j1etaWidth = jwc.getEtaWidth();
      j1phiWidth = jwc.getPhiWidth();
      j1etaWidthInECal = jwc.getEtaWidthInECal();
      j1phiWidthInECal = jwc.getPhiWidthInECal();
      j1etaWidthInHCal = jwc.getEtaWidthInHCal();
      j1phiWidthInHCal = jwc.getPhiWidthInHCal();
      j1PFCand12PtSum = jwc.getPFCand12PtSum();
      j1PFCandsPtSum = jwc.getPFCandsPtSum();
      j1PFCand12Ratio = jwc.getPFCand12Ratio(); //j1PFCand12PtSum/j1PFCandsPtSum
      //std::cout<<"j1PFCand12Ratio: "<<j1PFCand12Ratio<<std::endl;
      //std::cout<<"end of info about leading jet"<<std::endl;
    }
    else if(j == 1) {
      //std::cout<<"start of info about second leading jet"<<std::endl;	
      j2PT = jet.pt();
      j2Eta = jet.eta();
      j2Phi = jet.phi();
      j2CHdFr = jet.chargedHadronEnergyFraction();
      j2NHdFr = jet.neutralHadronEnergyFraction();
      j2CEmFr = jet.chargedEmEnergyFraction();
      j2NEmFr = jet.neutralEmEnergyFraction();
      j2PhoEFr = jet.photonEnergyFraction();
      j2EleEFr = jet.electronEnergyFraction();
      j2MuEFr = jet.muonEnergyFraction();
      j2CMuEFr = jet.chargedMuEnergyFraction();
      j2nCons = jet.nConstituents();
      j2nCarrying60 = jet.nCarrying(0.60);
      j2nCarrying90 = jet.nCarrying(0.90); //return # of consitutents carrying 90% of energy
      j2nCarrying92 = jet.nCarrying(0.92);
      j2nCarrying94 = jet.nCarrying(0.94);
      j2nCarrying96 = jet.nCarrying(0.96);
      j2nCarrying98 = jet.nCarrying(0.98);
      //std::cout<<"Info starts about Leading Jet"<<std::endl;
      //get all tracks in the jets. All PFCandidates hold a reference to a track.
      const reco::TrackRefVector j2tracks = jet.getTrackRefs();
      std::vector<reco::TrackRef> j2tracksRef;//make a copy of RefVector to sort it
      j2nTracks = j2tracks.size();
      if(j2nTracks >0){
	//std::cout<<"Number of tracks(Second Leading Jet): "<< j2nTracks<<std::endl;
	for (const auto &trkRef : j2tracks){
	  reco::Track trk = *trkRef;
	  j2tracksRef.push_back(trkRef);
	}
	const std::vector<reco::PFCandidatePtr> pfCands = jet.getPFConstituents(); //get pfCands to match with tracks 
	std::sort(j2tracksRef.begin(),j2tracksRef.end(),[](const reco::TrackRef& track1, const reco::TrackRef& track2)
	    {return track1->pt() > track2->pt();}); //sorting the copy of the TrackRefVector by pT
	//std::cout<<"TrackRef Copy begin"<<std::endl;
	if(j2tracksRef.size()>1){
	  j2trk1Quality = j2tracksRef.at(0)->quality(reco::TrackBase::TrackQuality::highPurity);//determining the quality of track
	  j2trk1Charge = j2tracksRef.at(0)->charge();
	  j2trk1PT = j2tracksRef.at(0)->pt();
	  if(j2trk1Charge!=0){
	    j2trk1PtError = j2tracksRef.at(0)->ptError(); //ptError for leading track
	    j2ptTrk1Ratio = (j2trk1PtError/j2trk1PT);}
	  j2trk1Eta = j2tracksRef.at(0)->eta();
	  j2trk1Phi = j2tracksRef.at(0)->phi();
	  //std::cout<<"Leading track Pt: "<< j2trk1PT<<std::endl;
	  j2trk2Quality = j2tracksRef.at(1)->quality(reco::TrackBase::TrackQuality::highPurity);
	  //std::cout<<"DPt/Pt: "<<j2ptTrk1Ratio<<std::endl;
	  //matching leading track candidate with PF-Candidate
	  j2dRCand1Trk1 = deltaR((*j2tracksRef.at(0)),(*pfCands[0]));
	  if (j2dRCand1Trk1 < 0.1){j2Trk1PFCand1Match = 1;}
	  j2trk2PT = j2tracksRef.at(1)->pt();
	  j2trk2Charge = j2tracksRef.at(1)->charge();
	  if(j2trk2Charge!=0){
	    j2trk2PtError = j2tracksRef.at(1)->ptError();//ptError of second leading track
	    j2ptTrk2Ratio = (j2trk2PtError/j2trk2PT);}
	  j2trk2Eta = j2tracksRef.at(1)->eta();
	  j2trk2Phi = j2tracksRef.at(1)->phi();
	  //std::cout<<"Second leading track Pt: "<< j2trk2PT<<std::endl;
	  //std::cout<<"DPt/Pt: "<<j2ptTrk2Ratio<<std::endl;
	  j2dRCand2Trk2 = deltaR((*j2tracksRef.at(1)),(*pfCands[1]));
	  if (j2dRCand2Trk2 < 0.1){j2Trk2PFCand2Match = 1;}
	  j2trk12PT = j2trk1PT+j2trk2PT;//access pt of the first and second track in list
	}
	else{
	  j2trk1Quality = j2tracksRef.at(0)->quality(reco::TrackBase::TrackQuality::highPurity);
	  j2trk1Charge = j2tracksRef.at(0)->charge();
	  j2trk1PT = j2tracksRef.at(0)->pt();
	  if(j2trk1Charge!=0){
	    j2trk1PtError = j2tracksRef.at(0)->ptError(); //ptError for leading track         
	    j2ptTrk1Ratio = (j2trk1PtError/j2trk1PT);}
	  //std::cout<<"DPt/Pt: "<<j2ptTrk1Ratio<<std::endl;
	  j2dRCand1Trk1 = deltaR((*j2tracksRef.at(0)),(*pfCands[0]));
	  if (j2dRCand1Trk1 < 0.1){j2Trk1PFCand1Match = 1;}
	  j2trk1Eta = j2tracksRef.at(0)->eta();
	  j2trk1Phi = j2tracksRef.at(0)->phi();
	  j2trk12PT = j2trk1PT;
	}
	//std::cout<<"SumPT of two leading Tracks: " << j2trk12PT << std::endl;
	//std::cout<<"PT of the Second Leading Jet: " << j2PT << std::endl;
	j2trk12Ratio = (j2trk12PT/j2PT);
	//std::cout<<"Trk12PT/j2PT: "<< j2trk12Ratio<<std::endl;
      }
      j2CMty = jet.chargedMultiplicity();
      j2NMty = jet.neutralMultiplicity();
      j2CHdMty = jet.chargedHadronMultiplicity();
      j2NHdMty = jet.neutralHadronMultiplicity();
      j2PhoMty = jet.photonMultiplicity();
      j2EleMty = jet.electronMultiplicity();
      j2MuMty = jet.muonMultiplicity();
      JetWidthCalculator jwc(jet);
      j2etaWidth = jwc.getEtaWidth();
      j2phiWidth = jwc.getPhiWidth();
      j2etaWidthInECal = jwc.getEtaWidthInECal();
      j2phiWidthInECal = jwc.getPhiWidthInECal();
      j2etaWidthInHCal = jwc.getEtaWidthInHCal();
      j2phiWidthInHCal = jwc.getPhiWidthInHCal();
      j2PFCand12PtSum = jwc.getPFCand12PtSum();
      j2PFCandsPtSum = jwc.getPFCandsPtSum();
      j2PFCand12Ratio = jwc.getPFCand12Ratio(); //j2PFCand12PtSum/j2PFCandsPtSum
      //std::cout<<"j2PFCand12Ratio: "<<j2PFCand12Ratio<<std::endl;

    }
    // std::cout<<"end of info about second leading jet"<<std::endl; 

    if(jet.pt() > 30.) { // Good jets have high pt
      HT += jet.pt();
      nGoodJets++;
    }
    //storing jetpT,jetEta,jetPhi and all other variables
    jetPx_.push_back(jet.px());
    jetPy_.push_back(jet.py());
    jetPz_.push_back(jet.pz());
    jetPt_.push_back(jet.pt());
    jetEn_.push_back(jet.energy());
    jetEta_.push_back(jet.eta());  
    jetPhi_.push_back(jet.phi());
    jetCEF_.push_back(jet.chargedEmEnergyFraction());
    jetNEF_.push_back(jet.neutralEmEnergyFraction());
    jetCHF_.push_back(jet.chargedHadronEnergyFraction());
    jetNHF_.push_back(jet.neutralHadronEnergyFraction());
    jetNCH_.push_back(jet.chargedMultiplicity());    
    jetNNP_.push_back(jet.neutralMultiplicity());
    jetMUF_.push_back(jet.muonEnergyFraction());
    jetHFHAE_.push_back( jet.HFHadronEnergy());
    jetHFEME_.push_back( jet.HFEMEnergy());
    jetNConstituents_.push_back(jet.numberOfDaughters());
    /*JetWidthCalculator jwc(jet);
      jetEtaWidth_.push_back(jwc.getEtaWidth());     
      jetPhiWidth_.push_back(jwc.getPhiWidth());
      jetEtaWidthInECal_.push_back(jwc.getEtaWidthInECal());
      jetPhiWidthInECal_.push_back(jwc.getPhiWidthInECal());
      jetEtaWidthInHCal_.push_back(jwc.getEtaWidthInHCal());
      jetPhiWidthInHCal_.push_back(jwc.getPhiWidthInHCal());*/

    jetPhotonEnergyFraction_.push_back(jet.photonEnergyFraction());
    jetJetArea_.push_back(jet.jetArea());
    jetMt_.push_back(jet.mt());
    jetMass_.push_back(jet.mass() );
    jetMaxDistance_.push_back(jet.maxDistance());
    jetPhiPhiMoment_.push_back(jet.phiphiMoment());
    jetConstituentEtaPhiSpread_.push_back(jet.constituentEtaPhiSpread());
    jetConstituentPtDistribution_.push_back(jet.constituentPtDistribution());
    jetPileup_.push_back(jet.pileup());
    jetN60_.push_back(jet.nCarrying(0.6));
    jetN90_.push_back(jet.nCarrying(0.9));

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

  for(uint32_t i = 0; i < electronHandle->size(); i++) {
    //const pat::Electron &electron = (*electronHandle)[i];
    const auto el = electronHandle->ptrAt(i);
    elePt_              .push_back(el->pt());
    eleEta_             .push_back(el->eta());
    elePhi_             .push_back(el->phi());

    eleCharge_          .push_back(el->charge());
    eleChargeConsistent_.push_back((Int_t)el->isGsfCtfScPixChargeConsistent());
    eleEn_              .push_back(el->energy());

    eleR9_              .push_back(el->r9());
    eleSCEn_            .push_back(el->superCluster()->energy());
    eleSCEta_           .push_back(el->superCluster()->eta());
    eleSCPhi_           .push_back(el->superCluster()->phi());
    eleSCRawEn_         .push_back(el->superCluster()->rawEnergy());
    eleSCEtaWidth_      .push_back(el->superCluster()->etaWidth());
    eleSCPhiWidth_      .push_back(el->superCluster()->phiWidth());
    eleHoverE_          .push_back(el->hcalOverEcal());
    eleEoverP_          .push_back(el->eSuperClusterOverP());
    eleEoverPout_       .push_back(el->eEleClusterOverPout());
    eleBrem_            .push_back(el->fbrem());
    eledEtaAtVtx_       .push_back(el->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_       .push_back(el->deltaPhiSuperClusterTrackAtVtx());
    eledEtaAtCalo_      .push_back(el->deltaEtaSeedClusterTrackAtCalo());
    eleSigmaIEtaIEta_   .push_back(el->sigmaIetaIeta()); ///new sigmaietaieta
    eleSigmaIPhiIPhi_   .push_back(el->sigmaIphiIphi());
    eleMissHits_        .push_back(el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));

    bool passConvVeto = !ConversionTools::hasMatchedConversion(*el,conversions,theBeamSpot->position());
    eleConvVeto_.push_back( (int) passConvVeto );

    reco::GsfTrackRef theTrack = el->gsfTrack();
    eleD0_.push_back( (-1) * theTrack->dxy(pv.position() ) );
    eleDz_.push_back( theTrack->dz( pv.position() ) );


    // VID calculation of (1/E - 1/p)
    if (el->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
    else if (!std::isfinite(el->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
    else  eleEoverPInv_.push_back((1.0 - el->eSuperClusterOverP())/el->ecalEnergy());

    float dEtaSeedValue = 
      el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() 
      ? 
      el->deltaEtaSuperClusterTrackAtVtx() 
      - el->superCluster()->eta() 
      + el->superCluster()->seed()->eta() 
      : std::numeric_limits<float>::max();
    eledEtaseedAtVtx_   .push_back(dEtaSeedValue);

    eleE1x5_            .push_back(el->e1x5());
    eleE2x5_            .push_back(el->e2x5Max());
    eleE5x5_            .push_back(el->e5x5());

    eleSigmaIEtaIEtaFull5x5_    .push_back(el->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_    .push_back(el->full5x5_sigmaIphiIphi());
    eleE1x5Full5x5_             .push_back(el->full5x5_e1x5());
    eleE2x5Full5x5_             .push_back(el->full5x5_e2x5Max());
    eleE5x5Full5x5_             .push_back(el->full5x5_e5x5());
    eleR9Full5x5_               .push_back(el->full5x5_r9());

    PassVeto_.push_back((*veto_id_decisions)[el]);
    PassLoose_.push_back((*loose_id_decisions)[el]);
    PassMedium_.push_back((*medium_id_decisions)[el]);
    PassTight_.push_back((*tight_id_decisions)[el]);
    PassHEEP_.push_back((*heep_id_decisions)[el]);

    reco::GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
    eleisoChargedHadrons_.push_back( pfIso.sumChargedHadronPt );
    eleisoNeutralHadrons_.push_back( pfIso.sumNeutralHadronEt );
    eleisoPhotons_.push_back( pfIso.sumPhotonEt );
    eleisoChargedFromPU_.push_back( pfIso.sumPUPt );



    nEle_++;

  }

  //Loop over photon collection
  for (size_t i = 0; i < photons->size(); ++i){
    const auto pho = photons->ptrAt(i);

    if( pho->pt() < 15 ) 
      continue;

    nPhotons_++;

    phopt_  .push_back( pho->pt() );
    phoeta_ .push_back( pho->superCluster()->eta() );
    phophi_ .push_back( pho->superCluster()->phi() );

    phohOverE_                .push_back( pho->hadTowOverEm() );
    phohasPixelSeed_          .push_back( (Int_t)pho->hasPixelSeed() );

    float chIso =  (*phoChargedIsolationMap)[pho];
    float worstchIso =  (*phoWorstChargedIsolationMap)[pho];
    float nhIso =  (*phoNeutralHadronIsolationMap)[pho];
    float phIso = (*phoPhotonIsolationMap)[pho];
    phoisoChargedHadrons_ .push_back( chIso );
    phoWorstisoChargedHadrons_ .push_back( worstchIso );
    phoisoNeutralHadrons_ .push_back( nhIso );
    phoisoPhotons_        .push_back( phIso );

    phoSigmaIEtaIEtaFull5x5_ .push_back(pho->full5x5_sigmaIetaIeta());

    phopx_.push_back(pho->px());
    phopy_.push_back(pho->py());
    phopz_.push_back(pho->pz());
    phoE_ .push_back(pho->energy());
    phor9_.push_back(pho->r9());

  }




  for (edm::View<pat::Muon>::const_iterator muon = recoMuons->begin(); muon != recoMuons->end(); ++muon) {
    if (muon->pt() < 3) continue;
    if (! (muon->isPFMuon() || muon->isGlobalMuon() || muon->isTrackerMuon())) continue;

    muPt_    .push_back(muon->pt());
    muEn_    .push_back(muon->energy());
    muEta_   .push_back(muon->eta());
    muPhi_   .push_back(muon->phi());
    muCharge_.push_back(muon->charge());
    muType_.push_back(muon->type());
    muD0_    .push_back(muon->muonBestTrack()->dxy(pv.position()));
    muDz_    .push_back(muon->muonBestTrack()->dz(pv.position()));
    muIsLooseID_ .push_back(muon->isLooseMuon());;
    muIsMediumID_.push_back(muon->isMediumMuon());
    muIsTightID_ .push_back(muon->isTightMuon(pv));
    muIsSoftID_  .push_back(muon->isSoftMuon(pv));
    muIsHighPtID_.push_back(muon->isHighPtMuon(pv));
    muBestTrkPtError_        .push_back(muon->muonBestTrack()->ptError());
    muBestTrkPt_             .push_back(muon->muonBestTrack()->pt());
    musegmentCompatibility_  .push_back(muon->segmentCompatibility());
    muchi2LocalPosition_     .push_back(muon->combinedQuality().chi2LocalPosition);
    mutrkKink_               .push_back(muon->combinedQuality().trkKink);
    const reco::TrackRef glbmu = muon->globalTrack();
    const reco::TrackRef innmu = muon->innerTrack();


    if (glbmu.isNull()) {
      muChi2NDF_ .push_back(-99.);
      muMuonHits_.push_back(-99);
    } else {
      muChi2NDF_.push_back(glbmu->normalizedChi2());
      muMuonHits_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
    }


    if (innmu.isNull()) {
      muInnerD0_     .push_back(-99.);
      muInnerDz_     .push_back(-99.);
      muTrkLayers_   .push_back(-99);
      muPixelLayers_ .push_back(-99);
      muPixelHits_   .push_back(-99);
      muTrkQuality_  .push_back(-99);

      muInnervalidFraction_ .push_back(-99);
    } else {
      muInnerD0_     .push_back(innmu->dxy(pv.position()));
      muInnerDz_     .push_back(innmu->dz(pv.position()));
      muTrkLayers_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
      muPixelLayers_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
      muPixelHits_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
      muTrkQuality_  .push_back(innmu->quality(reco::TrackBase::highPurity));

      muInnervalidFraction_ .push_back(innmu->validFraction());
    }

    muStations_ .push_back(muon->numberOfMatchedStations());
    muMatches_  .push_back(muon->numberOfMatches());
    muIsoTrk_   .push_back(muon->trackIso());
    muPFChIso_  .push_back(muon->pfIsolationR04().sumChargedHadronPt);
    muPFPhoIso_ .push_back(muon->pfIsolationR04().sumPhotonEt);
    muPFNeuIso_ .push_back(muon->pfIsolationR04().sumNeutralHadronEt);
    muPFPUIso_  .push_back(muon->pfIsolationR04().sumPUPt);


    nMu_++;
  }


  pfMET = -99, pfMETPhi = -99, pfMETsumEt_ = -99, pfMETmEtSig_ = -99, pfMETSig_ = -99;

  pfMET = (*pfMETs)[0].pt();
  pfMETPhi = (*pfMETs)[0].phi();
  pfMETsumEt_ = (*pfMETs)[0].sumEt();
  pfMETmEtSig_ = ((*pfMETs)[0].mEtSig() < 1.e10) ? (*pfMETs)[0].mEtSig() : 0;
  pfMETSig_ = ((*pfMETs)[0].significance() < 1.e10) ? (*pfMETs)[0].significance() : 0;    

  caloMET = -99;
  caloMET = (*caloMETs)[0].pt();

  //std::cout << "TotalET = " << totalET << "; nJets = " << nJets 
  //    << "; HT = " << HT << "; nGoodJets = " << nGoodJets 
  //  << "; pfMET = (" << pfMET << ", " << pfMETPhi << ")" 
  //  << "; nPFCands = " << pfCands->size() << std::endl;

  tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
  void 
JetAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
  void 
JetAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JetAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetAnalyzer);
