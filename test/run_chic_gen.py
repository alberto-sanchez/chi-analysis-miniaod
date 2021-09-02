# cfg for running on AOD

input_filename = '/store/hidata/PARun2016C/PADoubleMuon/AOD/PromptReco-v1/000/285/505/00000/18609D6A-4CAF-E611-AA50-02163E01184F.root'
ouput_filename = 'rootuple_chic_hi_mc2.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source       = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))

process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options      = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.rootuple = cms.EDAnalyzer('GenChicRootupler')

process.p = cms.Path(process.rootuple)
