# cfg for running on AOD

input_filename = '/store/hidata/PARun2016C/PADoubleMuon/AOD/PromptReco-v1/000/285/505/00000/18609D6A-4CAF-E611-AA50-02163E01184F.root'
ouput_filename = 'rootuple_chic_hi.root'

import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_v15', '')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents    = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source       = cms.Source("PoolSource",fileNames = cms.untracked.vstring(input_filename))

process.TFileService = cms.Service("TFileService",fileName = cms.string(ouput_filename))
process.options      = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.load('HeavyIonsAnalysis.Configuration.hfCoincFilter_cff')
process.primaryVertexFilterPA = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && abs(z) <= 25 && position.Rho <= 2 && tracksSize >= 2"),
    filter = cms.bool(True),
)
process.noScraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

## ==== Trigger matching
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import *
## with some customization  O: values copied over from the HI onia trees (or T&P - same)
process.muonL1Info.maxDeltaR = 0.3
process.muonL1Info.maxDeltaEta = 0.2
process.muonL1Info.fallbackToME1 = True
process.muonMatchHLTL1.maxDeltaR = 0.3
process.muonMatchHLTL1.maxDeltaEta = 0.2
process.muonMatchHLTL1.fallbackToME1 = True
process.muonMatchHLTL2.maxDeltaR = 0.3
process.muonMatchHLTL2.maxDPtRel = 10.0
process.muonMatchHLTL3.maxDeltaR = 0.1
process.muonMatchHLTL3.maxDPtRel = 10.0

## For L1 muons
addHLTL1Passthrough(process)
useL1Stage2Candidates(process)
process.patTrigger.collections.remove("hltL1extraParticles")
process.patTrigger.collections.append("hltGmtStage2Digis:Muon")
process.muonMatchHLTL1.matchedCuts = cms.string('coll("hltGmtStage2Digis:Muon")')
process.muonMatchHLTL1.useStage2L1 = cms.bool(True)
process.muonMatchHLTL1.useMB2InOverlap = cms.bool(True)
process.muonMatchHLTL1.preselection = cms.string("")
appendL1MatchingAlgo(process)

process.selectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('patMuonsWithTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
            ' && abs(innerTrack.dxy) < 0.3'
            ' && abs(innerTrack.dz)  < 20.'
            ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
            ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
            ' && innerTrack.quality(\"highPurity\")'
            ' && ((abs(eta) <= 0.9 && pt > 2.5) || (0.9 < abs(eta) <= 2.4 && pt > 1.0))'
   ),
   filter = cms.bool(True)
)

### Trigger selection
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.triggerResultsFilter.triggerConditions = cms.vstring('HLT_PAL1DoubleMu*_v*')
process.triggerResultsFilter.hltResults = cms.InputTag("TriggerResults","","HLT")
process.triggerResultsFilter.l1tResults = cms.InputTag("gtStage2Digis") #O: needs to be this
process.triggerResultsFilter.throw = False

### Filter sequence
process.fastFilter = cms.Sequence(process.hfCoincFilter + process.primaryVertexFilterPA + process.noScraping + process.triggerResultsFilter)

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('selectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isTrackerMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isTrackerMuon")
process.onia2MuMuPAT.dimuonSelection=cms.string("mass > 2.0 && abs(daughter('muon1').innerTrack.dz - daughter('muon2').innerTrack.dz) < 25 && pt>5.0")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)
process.onia2MuMuPAT.resolvePileUpAmbiguity = cms.bool(True)

process.Onia2MuMuFiltered = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("onia2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.0 < mass && pt > 3. && charge==0 && userFloat('vProb') > 0.0"),
      do_trigger_match    = cms.bool(True),
      HLTFilters          = cms.vstring(
   'HLT_PAL1DoubleMu0_v*',
   'HLT_PAL1DoubleMuOpen_v*',
   'HLT_PAL2DoubleMu0_v*',
   'HLT_PAL3DoubleMu0_v*'
      ),
)

process.DiMuonCounter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("Onia2MuMuFiltered"),
    minNumber = cms.uint32(1),
    filter    = cms.bool(True)
)

# The low energy photons are reconstructed here.
import HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi
process.oniaPhotonCandidates = HeavyFlavorAnalysis.Onia2MuMu.OniaPhotonConversionProducer_cfi.PhotonCandidates.clone()

process.chiProducer = cms.EDProducer('OniaPhotonProducer',
    conversions     = cms.InputTag("oniaPhotonCandidates","conversions"),
    dimuons         = cms.InputTag("Onia2MuMuFiltered"),
    pi0OnlineSwitch = cms.bool(False),
    deltaMass       = cms.vdouble(0.0,2.0),
    dzmax           = cms.double(99.),
    triggerMatch    = cms.bool(False)  # trigger match is performed in Onia2MuMuFiltered
)

process.chiFitter1S = cms.EDProducer('OniaPhotonKinematicFit',
                          chi_cand = cms.InputTag("chiProducer"),
                          upsilon_mass = cms.double(3.0969), # GeV   1S = 9.46030   2S = 10.02326    3S = 10.35520  J/psi=3.0969
                          product_name = cms.string("y1S"),
			  beamSpotTag = cms.InputTag('offlineBeamSpot')
                         )

process.chiSequence = cms.Sequence(
   process.fastFilter *
   process.patMuonsWithTriggerSequence *
   process.selectedMuons *
   process.onia2MuMuPAT *
   process.Onia2MuMuFiltered *
   process.DiMuonCounter *
   process.oniaPhotonCandidates *
   process.chiProducer *
   process.chiFitter1S
)

process.rootuple = cms.EDAnalyzer('HIchicRootupler',
                          chic_cand = cms.InputTag("chiProducer"),
                          psi_cand = cms.InputTag("Onia2MuMuFiltered"),
                          refit1S  = cms.InputTag("chiFitter1S","y1S"),
                          primaryVertices = cms.InputTag("offlinePrimaryVertices"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
			  beamSpotTag = cms.InputTag('offlineBeamSpot'),
                          isMC = cms.bool(False),
                          only_best =  cms.bool(False),
                          FilterNames = cms.vstring(
   'HLT_PAL1DoubleMu0',
   'HLT_PAL1DoubleMuOpen',
   'HLT_PAL2DoubleMu0',
   'HLT_PAL3DoubleMu0',
   'HLT_PAFullTracks_Multiplicity120',
   'HLT_PAFullTracks_Multiplicity150',
   'HLT_PAFullTracks_Multiplicity185_part1',
   'HLT_PAFullTracks_Multiplicity185_part2',
   'HLT_PAFullTracks_Multiplicity185_part3',
   'HLT_PAFullTracks_Multiplicity185_part4',
   'HLT_PAFullTracks_Multiplicity185_part5',
   'HLT_PAFullTracks_Multiplicity185_part6',
   'HLT_PAFullTracks_Multiplicity220',
   'HLT_PAFullTracks_Multiplicity250',
   'HLT_PAFullTracks_Multiplicity280'
                          )
)

process.p = cms.Path(process.chiSequence*process.rootuple)
