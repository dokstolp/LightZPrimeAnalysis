#multicrab

dataset = {
  'GJets_HT-40To100': '/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM',
  'GJets_HT-100To200': '/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v2/AODSIM',
  'GJets_HT-200To400': '/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'GJets_HT-400To600': '/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'GJets_HT-600ToInf': '/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'TTJets': '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'QCD_Pt-80to120': '/QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'QCD_Pt-120to170': '/QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'QCD_Pt-170to300': '/QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'QCD_Pt-300toInf': '/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'W1Jets': '/W1JetsToLNu_NuPt-200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'W2Jets': '/W2JetsToLNu_NuPt-200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'W3Jets': '/W3JetsToLNu_NuPt-200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM',
  'W4Jets': '/W3JetsToLNu_NuPt-200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16DR80-PUSpring16_80X_mcRun2_asymptotic_2016_v3-v1/AODSIM'
 }

if __name__ == '__main__':
 from CRABAPI.RawCommand import crabCommand

def submit(config):
 res = crabCommand('submit', config = config)

from CRABClient.UserUtilities import config
config = config()
name = 'Zprime_Ntuples'
config.General.workArea = 'crab_workArea_'+name
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.pluginName = 'Analysis'
config.Data.publication = False
config.Site.storageSite = 'T2_US_Wisconsin'
config.JobType.psetName = 'testJetAnalyzer.py'
#config.JobType.inputFiles = ['Spring16_25nsV1_MC_L2Relative_AK8PFchs.txt','Spring16_25nsV1_MC_L3Absolute_AK8PFchs.txt']
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#Already submitted:
#listOfSamples = ['ZNuNuGJets','ZLLGJets','WGJets','WToMuNu','WToTauNu','GJets_HT-40To100','GJets_HT-100To200','GJets_HT-200To400','GJets_HT-400To600','GJets_HT-600ToInf','TTGJets','QCD_Pt-80to120','QCD_Pt-120to170','QCD_Pt-170to300','QCD_Pt-300toInf','DYJetsToLL_M-50_HT-100to200','DYJetsToLL_M-50_HT-200to400','DYJetsToLL_M-50_HT-400to600','DYJetsToLL_M-50_HT-600toInf','DMAVMx-1000Mv-10000','DMAVMx-1000Mv-10','DMAVMx-1000Mv-1995','DMAVMx-10Mv-10000','DMAVMx-10Mv-100','DMAVMx-10Mv-10','DMAVMx-10Mv-15','DMAVMx-150Mv-10','DMAVMx-150Mv-295','DMAVMx-150Mv-500','DMAVMx-1Mv-1000','DMAVMx-1Mv-100','DMAVMx-1Mv-10','DMAVMx-1Mv-2000','DMAVMx-1Mv-20','DMAVMx-1Mv-300','DMAVMx-1Mv-50','DMAVMx-500Mv-10000','DMAVMx-500Mv-10','DMAVMx-500Mv-2000','DMAVMx-500Mv-995','DMAVMx-50Mv-200','DMAVMx-50Mv-300','DMAVMx-50Mv-50','DMAVMx-50Mv-95','DMEWKscalarMx-100','DMEWKscalarMx-10','DMEWKscalarMx-1','DMEWKscalarMx-200','DMEWKscalarMx-400','DMEWKscalarMx-50','DMEWKscalarMx-800','DMVMx-1000Mv-10000','DMVMx-1000Mv-10','DMVMx-1000Mv-1995','DMVMx-10Mv-100','DMVMx-10Mv-15','DMVMx-10Mv-50','DMVMx-150Mv-10000','DMVMx-150Mv-1000','DMVMx-150Mv-10','DMVMx-150Mv-200','DMVMx-150Mv-295','DMVMx-1Mv-1000','DMVMx-1Mv-10','DMVMx-1Mv-2000','DMVMx-1Mv-200','DMVMx-1Mv-20','DMVMx-1Mv-300','DMVMx-1Mv-500','DMVMx-1Mv-50','DMVMx-500Mv-10000','DMVMx-500Mv-10','DMVMx-500Mv-2000','DMVMx-500Mv-500','DMVMx-50Mv-10000','DMVMx-50Mv-200','DMVMx-50Mv-300','DMVMx-50Mv-50','DMVMx-50Mv-95']
listOfSamples = ['W1Jets','W2Jets','W3Jets','W4Jets']
for sample in listOfSamples:
  config.General.requestName = sample
  config.Data.inputDataset = dataset[sample]
  config.Data.unitsPerJob = 1
  config.Data.totalUnits = -1
  config.Data.outLFNDirBase = '/store/user/gomber/'+name
  submit(config)
