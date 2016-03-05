import FWCore.ParameterSet.Config as cms

process = cms.Process('LightZPrimeGenAnalyzer')

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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:/nfs_scratch/dasu/LightZPrimeDecay/GEN.root'))

process.load("LightZPrimeAnalysis.LightZPrimeAnalyzer.LightZPrimeGenAnalyzer_cfi")

process.p = cms.Path(process.LightZPrimeGenAnalyzer)
