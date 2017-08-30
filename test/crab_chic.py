from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

myname='rootuple-chic-2017'
mydata='/Charmonium/Run2017B-PromptReco-v1/MINIAOD'
myjson='Cert_294927-297723_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt'

config.General.requestName = myname
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run-chic-miniaod.py'
config.JobType.outputFiles = ['rootuple.root']

config.Data.inputDataset = mydata
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.lumiMask = myjson

config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.publishDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outputDatasetTag  = myname
config.Site.storageSite = 'T3_US_FNALLPC'
