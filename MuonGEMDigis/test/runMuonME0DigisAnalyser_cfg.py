import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("me0", eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
#process.load('Configuration.StandardSequences.Generator_cff')
#process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step3.root')
)

process.me0Digis = cms.EDAnalyzer("MuonME0DigisAnalyser",
    simTrackTag = cms.InputTag('g4SimHits'),
    simHitTag = cms.InputTag('g4SimHits', "MuonME0Hits"),
    me0DigiTag = cms.InputTag("simMuonME0Digis"), 
    me0DigiSimLinkTag = cms.InputTag("simMuonME0Digis", "ME0"),
    me0RecHitTag = cms.InputTag("me0RecHits"),
    me0SegmentTag = cms.InputTag("me0Segments"),
    min_pt = cms.double(5.0), # minimum pt for muon SimTrack
    min_quality = cms.double(0.6),
    min_num_layers = cms.uint32(4),
    min_digis = cms.uint32(3),
    max_muons = cms.uint32(3),
    # genParticleTag = cms.InputTag("genParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuonME0DigisAnalysis.root')
)

process.p = cms.Path(process.me0Digis)
