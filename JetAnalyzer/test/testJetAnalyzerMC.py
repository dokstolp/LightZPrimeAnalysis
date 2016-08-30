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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:0036B88E-5009-E611-9C57-0CC47A009148.root'))
#'/store/data/Run2016B/MET/AOD/PromptReco-v2/000/274/968/00000/FABA80F3-2F32-E611-88CB-02163E013556.root'))

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
#process.JetAnalyzer.addFilterInfo=cms.bool(True)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#turn on VID producer,
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V2_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.p = cms.Path(process.egmGsfElectronIDSequence*process.JetAnalyzerMC)