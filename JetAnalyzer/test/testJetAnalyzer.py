import FWCore.ParameterSet.Config as cms

process = cms.Process('JetAnalyzer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/data/Run2016B/MET/AOD/PromptReco-v2/000/274/968/00000/FABA80F3-2F32-E611-88CB-02163E013556.root'))
#'file:/cms/uhussain/ZprimeDM/CMSSW_7_6_1/src/Zprime_UH-Mzp1_MET300-1000_3b.root'))
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/Asympt25ns_MCRUN2_74_V9-v2/50000/1E44B2CF-C127-E511-945D-008CFA152144.root'))
#'/store/mc/RunIIFall15MiniAODv1/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext1-v1/00000/50E50F96-E308-E611-9688-0CC47A6C0716.root'))
#'/store/mc/RunIISpring16MiniAODv2/TTToSemiLeptonic_13TeV_ScaleUp-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/3E18D980-1F20-E611-9CE6-0CC47A6C183A.root'))
#'file:/cms/uhussain/ZprimeDM/CMSSW_7_6_1/src/Zprime_SUS-RunIIWinter15pLHEMiniAODv2-1000_4.root'))

process.TFileService = cms.Service("TFileService", fileName = cms.string('JetAnalyzer_3triggers.root'))
#addFilterInfoAOD_ = True
process.load("LightZPrimeAnalysis.JetAnalyzer.JetAnalyzer_cfi")
process.load("LightZPrimeAnalysis.JetAnalyzer.ggMETFilters_cff")
#process.JetAnalyzer.addFilterInfo=cms.bool(True)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#    src = cms.InputTag("prunedGenParticles"),
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex  = cms.untracked.bool(True)
)

#process.p = cms.Path(process.printTree+process.JetAnalyzer)

process.p = cms.Path(process.JetAnalyzer)
