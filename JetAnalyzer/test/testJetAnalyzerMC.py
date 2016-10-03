import FWCore.ParameterSet.Config as cms

process = cms.Process('JetAnalyzerMC')

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

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:0036B88E-5009-E611-9C57-0CC47A009148.root'))
#'/store/mc/RunIISpring16DR80/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v3/00000/0036B88E-5009-E611-9C57-0CC47A009148.root'))
#
# Set up electron ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.AOD

process.TFileService = cms.Service("TFileService", fileName = cms.string('JetAnalyzerMC.root'))
#addFilterInfoAOD_ = True
process.load("LightZPrimeAnalysis.JetAnalyzer.JetAnalyzerMC_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#turn on VID producer,
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#    src = cms.InputTag("prunedGenParticles"),
#    src = cms.InputTag("genParticles"),
#    printP4 = cms.untracked.bool(False),
#    printPtEtaPhi = cms.untracked.bool(True),
#    printVertex = cms.untracked.bool(False),
#    printStatus = cms.untracked.bool(True),
#    printIndex  = cms.untracked.bool(True)
#)

#process.p = cms.Path(process.printTree+process.JetAnalyzerMC)

process.p = cms.Path(process.egmGsfElectronIDSequence*process.JetAnalyzerMC)
#dump_file = open("dump_file_withfilter.py", "w")
#dump_file.write(process.dumpPython())
