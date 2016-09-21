#multicrab

dataset = {
  'dataset1': '/MET/Run2016B-PromptReco-v2/AOD',
  'dataset2': '/MET/Run2016C-PromptReco-v2/AOD',
  'dataset3': '/MET/Run2016D-PromptReco-v2/AOD'
}

if __name__ == '__main__':
 from CRABAPI.RawCommand import crabCommand

def submit(config):
 res = crabCommand('submit', config = config)

from CRABClient.UserUtilities import config
config = config()
name = 'Zprime_Ntuples_Sep20'
config.General.workArea = 'crab_'+name
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.Data.publication = False
config.Site.storageSite = 'T2_US_Wisconsin'
config.JobType.psetName = 'testJetAnalyzer.py'
config.JobType.outputFiles = ['JetAnalyzer.root']
config.section_('Data') 
#config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt'
#Already submitted:
listOfSamples = ['dataset1','dataset2','dataset3']
for sample in listOfSamples:  
  config.General.requestName = sample
  config.Data.inputDataset = dataset[sample]
  config.Data.unitsPerJob = 15
  config.Data.totalUnits = -1
  config.Data.outLFNDirBase = '/store/user/uhussain/'+name
  submit(config)
