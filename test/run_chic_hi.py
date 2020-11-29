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

# get patmuons, no trigger info
import PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi
process.oniaPATMuonsWithoutTrigger = PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi.patMuons.clone(
    muonSource = 'muons',
    embedTrack          = True,
    embedCombinedMuon   = True,
    embedStandAloneMuon = True,
    embedPFCandidate    = False,
    embedCaloMETMuonCorrs = cms.bool(False),
    embedTcMETMuonCorrs   = cms.bool(False),
    embedPfEcalEnergy     = cms.bool(False),
    embedPickyMuon = False,
    embedTpfmsMuon = False,
    userIsolation = cms.PSet(),   # no extra isolation beyond what's in reco::Muon itself
    isoDeposits = cms.PSet(),     # no heavy isodeposits
    addGenMatch = False,          # no mc
)

# select soft-patmuons
process.selectedMuons = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('oniaPATMuonsWithoutTrigger'),
   cut = cms.string('muonID(\"TMOneStationTight\")'
            ' && abs(innerTrack.dxy) < 0.3'
            ' && abs(innerTrack.dz)  < 20.'
            ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
            ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
            ' && innerTrack.quality(\"highPurity\")'
            ' && (abs(eta) <= 2.4 && pt > 1.3)'
   ),
   filter = cms.bool(True)
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
   triggerConditions = cms.vstring(
   'HLT_PAL1DoubleMu0_v*',
   'HLT_PAL1DoubleMuOpen_v*',
   'HLT_PAL2DoubleMu0_v*',
   'HLT_PAL3DoubleMu0_v*',
   'HLT_PAFullTracks_Multiplicity120_v*',
   'HLT_PAFullTracks_Multiplicity150_v*',
   'HLT_PAFullTracks_Multiplicity185_part1_v*',
   'HLT_PAFullTracks_Multiplicity185_part2_v*',
   'HLT_PAFullTracks_Multiplicity185_part3_v*',
   'HLT_PAFullTracks_Multiplicity185_part4_v*',
   'HLT_PAFullTracks_Multiplicity185_part5_v*',
   'HLT_PAFullTracks_Multiplicity185_part6_v*',
   'HLT_PAFullTracks_Multiplicity220_v*',
   'HLT_PAFullTracks_Multiplicity250_v*',
   'HLT_PAFullTracks_Multiplicity280_v*'
   ),
   hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
   l1tResults = cms.InputTag( "" ),
   throw = cms.bool(False)
)

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('selectedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlinePrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.7 < mass && mass < 3.5")
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

process.Onia2MuMuFiltered = cms.EDProducer('DiMuonFilter',
      OniaTag             = cms.InputTag("onia2MuMuPAT"),
      singlemuonSelection = cms.string(""),
      dimuonSelection     = cms.string("2.7 < mass && mass < 3.5 && pt > 0. && charge==0 && userFloat('vProb') > 0.01"),
      do_trigger_match    = cms.bool(False),
      HLTFilters          = cms.vstring(
   'HLT_PAL1DoubleMu0_v*',
   'HLT_PAL1DoubleMuOpen_v*',
   'HLT_PAL2DoubleMu0_v*',
   'HLT_PAL3DoubleMu0_v*',
   'HLT_PAFullTracks_Multiplicity120_v*',
   'HLT_PAFullTracks_Multiplicity150_v*',
   'HLT_PAFullTracks_Multiplicity185_part1_v*',
   'HLT_PAFullTracks_Multiplicity185_part2_v*',
   'HLT_PAFullTracks_Multiplicity185_part3_v*',
   'HLT_PAFullTracks_Multiplicity185_part4_v*',
   'HLT_PAFullTracks_Multiplicity185_part5_v*',
   'HLT_PAFullTracks_Multiplicity185_part6_v*',
   'HLT_PAFullTracks_Multiplicity220_v*',
   'HLT_PAFullTracks_Multiplicity250_v*',
   'HLT_PAFullTracks_Multiplicity280_v*'
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
    dzmax           = cms.double(0.5),
    triggerMatch    = cms.bool(False)  # trigger match is performed in Onia2MuMuFiltered
)

process.chiFitter1S = cms.EDProducer('OniaPhotonKinematicFit',
                          chi_cand = cms.InputTag("chiProducer"),
                          upsilon_mass = cms.double(3.0969), # GeV   1S = 9.46030   2S = 10.02326    3S = 10.35520  J/psi=3.0969
                          product_name = cms.string("y1S")
                         )

process.chiSequence = cms.Sequence(
   process.oniaPATMuonsWithoutTrigger *
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
                          isMC = cms.bool(False),
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
