import FWCore.ParameterSet.Config as cms

process = cms.Process('JetAnalyzer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.Geometry.GeometryRecoDB_cff")


process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/data/Run2016C/SinglePhoton/AOD/23Sep2016-v1/90000/009435BC-1F86-E611-9241-0CC47A6C1034.root'))#
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/data/Run2016D/MET/AOD/PromptReco-v2/000/276/315/00000/286DBB5B-F444-E611-AA93-02163E011B53.root'))#
# Set up electron ID (VID framework)



process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process, outputModules = [] )

#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.AOD

process.TFileService = cms.Service("TFileService", fileName = cms.string('JetAnalyzer.root'))
#addFilterInfoAOD_ = True
process.load("LightZPrimeAnalysis.JetAnalyzer.JetAnalyzer_cfi")
process.load("LightZPrimeAnalysis.JetAnalyzer.ggMETFilters_cff")
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


process.p = cms.Path(process.egmGsfElectronIDSequence*process.photonIDValueMapProducer*process.ggMETFiltersSequence*process.JetAnalyzer)
#dump_file = open("dump_file_withfilter.py", "w")
#dump_file.write(process.dumpPython())
