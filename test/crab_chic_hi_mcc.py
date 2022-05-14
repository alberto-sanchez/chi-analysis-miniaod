from CRABClient.UserUtilities import config
config = config()

import datetime, time
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%y%m%d_%H%M')

ii=6

mydbs='phys03'
if ii<=1:
  mydata='/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v9-4e07c3c67d0ff0e1e9aecbbcb6a514fc/USER'
  myname='chic_hi_pPb_mcc_v9'
if ii==2:
  mydata='/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_RECO_v8-61416874db099c53202c8cb2d81ec4a3/USER'
  myname='chic_hi_pPb_mcc_v8'
if ii==3 :
  mydata='/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_Pbp_RECO_v8-61416874db099c53202c8cb2d81ec4a3/USER'
  myname='chic_hi_Pbp_mcc_v8'
if ii==4 :
  mydata='/Chi_c_pPb8TeV_privateMC_GEN/okukral-Chi_c_pPb8TeV_MC_Pbp_RECO_v9-4e07c3c67d0ff0e1e9aecbbcb6a514fc/USER'
  myname='chic_hi_Pbp_mcc_v9'
if ii==5 :
  mydata ='/Chic1Chic2_JpsiTogg_MuMuTogg_pThat4p5_Pbp-EmbEPOS_8p16_pythia8_evtgen/pPb816Summer16DR-80X_mcRun2_pA_v4-v7/AODSIM'
  myname='chic_hi_mc_v12'
  mydbs ='global'
if ii>=6 :
  mydata ='/Chic1Chic2_JpsiTogg_MuMuTogg_pThat4p5_pPb-EmbEPOS_8p16_pythia8_evtgen/pPb816Summer16DR-80X_mcRun2_pA_v4-v5/AODSIM'
  myname='chic_hi_pPb_mc_v12'
  mydbs ='global'


config.General.workArea = 'crab_jobs'
config.General.requestName = myname+'-'+st
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_chic_hi_mcc.py'
config.JobType.outputFiles = ['rootuple_chic_hi_mc.root']
#config.JobType.allowUndistributedCMSSW = True

config.Data.inputDataset = mydata

config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10

config.Data.outLFNDirBase = '/store/user/asanchez'
config.Data.publication = False
config.Data.publishDBS  = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
#config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSReader/'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.outputDatasetTag  = myname
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T2_CH_CERN'
#config.Data.outLFNDirBase = '/store/group/phys_bphys/asanchez'
