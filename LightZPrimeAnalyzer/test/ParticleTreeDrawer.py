import FWCore.ParameterSet.Config as cms

process = cms.Process('ParticleTreeDrawer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/hdfs/store/user/tuanqui/ZprimeM10_STEP2-step2/step2-step1.root'))

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag("genParticles"),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex  = cms.untracked.bool(True)
)

process.p = cms.Path(process.printTree)
