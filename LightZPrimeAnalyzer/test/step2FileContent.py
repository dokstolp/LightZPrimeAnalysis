import FWCore.ParameterSet.Config as cms

process = cms.Process('LightZPrimeAnalyzer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
#'file:/hdfs/store/user/uhussain/Step2-GENSIM-Zprime_SUS-RunIIWinter15pLHE-10000_2/Zprime_SUS-RunIIWinter15pLHE-10000_2-Zprime_SUS-RunIIWinter15pLHE-10000_18_numEvent500.root'
'file:/hdfs/store/user/tuanqui/ZprimeM10_STEP2-step2/step2-step1.root'
))

process.load("LightZPrimeAnalysis.LightZPrimeAnalyzer.LightZPrimeAnalyzer_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex  = cms.untracked.bool(True)
)

process.p = cms.Path(process.printTree+process.LightZPrimeAnalyzer)
