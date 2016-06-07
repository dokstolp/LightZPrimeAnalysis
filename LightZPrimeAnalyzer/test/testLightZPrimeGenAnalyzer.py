import FWCore.ParameterSet.Config as cms

process = cms.Process('LightZPrimeGenAnalyzer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/nfs_scratch/dasu/LightZPrimeDecay/GEN-1GeV.root'))

process.TFileService = cms.Service("TFileService", fileName = cms.string('LightZPrimeGenAnalyzer-1GeV.root'))

process.load("LightZPrimeAnalysis.LightZPrimeAnalyzer.LightZPrimeGenAnalyzer_cfi")

process.p = cms.Path(process.LightZPrimeGenAnalyzer)
