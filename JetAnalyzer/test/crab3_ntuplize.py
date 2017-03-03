from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'ntuples_df1en0'
config.General.workArea = 'darkPhoton'

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'testJetAnalyzer_mc.py'
config.JobType.outputFiles = ['darkPhoton_df1en0.root']
config.JobType.pyCfgParams = ['isRealData=0']

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
NJOBS = 500
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/group/lpcgg/dokstolp'
config.Data.outputDatasetTag = 'df1en0'
config.Data.outputPrimaryDataset = 'ntuples'
config.Data.userInputFiles = list(open('RecoOuts/reco_file_df1en0.txt'))

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'cmseos.fnal.gov'
