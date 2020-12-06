from CRABClient.UserUtilities import config
config = config()

import datetime, time
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d-%H-%M')

ii=1

myrun='run_chic_hi.py'
mydata='/PADoubleMuon/PARun2016C-PromptReco-v1/AOD'

if ii==1:
  myname='chic_hi_pPb'
  mylumi='Cert_285479-285832_HI8TeV_PromptReco_pPb_Collisions16_JSON_NoL1T_MuonPhys.txt'
else :
  myname='chic_hi_Pbp'
  mylumi='Cert_285952-286496_HI8TeV_PromptReco_Pbp_Collisions16_JSON_NoL1T_MuonPhys.txt'

config.General.requestName = myname+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = myrun
config.JobType.outputFiles = ['rootuple_chic_hi.root']
config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = mydata
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = mylumi

config.Data.outLFNDirBase = '/store/user/asanchez'
config.Data.publication = False
config.Data.publishDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outputDatasetTag  = myname
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Data.outLFNDirBase = '/store/group/phys_bphys/asanchez'
