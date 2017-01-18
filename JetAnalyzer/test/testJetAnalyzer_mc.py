import FWCore.ParameterSet.Config as cms

from LightZPrimeAnalysis.JetAnalyzer.JetAnalyzer_cfi import *

#process = cms.Process('JetAnalyzer')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'

process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryRecoDB_cff")


process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
#process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/data/Run2016D/MET/AOD/PromptReco-v2/000/276/315/00000/286DBB5B-F444-E611-AA93-02163E011B53.root'))#
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_1.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_2.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_3.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_4.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_5.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_6.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_7.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_8.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_9.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_10.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_11.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_12.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_13.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_14.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_15.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_16.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_17.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_18.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_19.root','file:root://cmsxrootd.hep.wisc.edu//store/user/uhussain/Zprime_3a-Zprime_UH-Mzp1_MET300_pLHE-10000_3b/Zprime_3b_20.root'))#'/store/mc/RunIISpring16DR80/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/00000/0062C1C4-8419-E611-9214-002590A83308.root'))## Set up electron ID (VID framework)
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/mc/RunIISpring16DR80/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/00000/00773447-C000-E611-B4E6-0CC47A4D761A.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('/store/mc/RunIISpring16DR80/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_v3_ext1-v1/00000/04D527BE-0206-E611-BD65-002590FD5A72.root'))
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring('file:pickevents_tracks.root'))

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

process.TFileService = cms.Service("TFileService", fileName = cms.string('JetAnalyzer_ZJets.root'))
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


process.p = cms.Path(process.jetSequence*process.egmGsfElectronIDSequence*process.photonIDValueMapProducer*process.ggMETFiltersSequence*process.JetAnalyzer)
#dump_file = open("dump_file_withfilter.py", "w")
#dump_file.write(process.dumpPython())
