import FWCore.ParameterSet.Config as cms

process = cms.Process('ParticleTreeDrawer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('inputFile', 'file:EDMFileWithGenParticles.root', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Input file')
options.register('sourceTag', 'genParticles', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'genParticle input tag name')
options.register('maxEvents', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Maximum number of events')
options.parseArguments()

process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(options.inputFile))

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
    src = cms.InputTag(options.sourceTag),
    printP4 = cms.untracked.bool(False),
    printPtEtaPhi = cms.untracked.bool(True),
    printVertex = cms.untracked.bool(False),
    printStatus = cms.untracked.bool(True),
    printIndex  = cms.untracked.bool(True)
)

process.p = cms.Path(process.printTree)
