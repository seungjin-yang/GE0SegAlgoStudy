import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11M9_cff import Phase2C11M9
process = cms.Process('reRECO',Phase2C11M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D67Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/users/jlee/CMSSW_11_3_0_pre1/src/31021.0_TenMuE_0_200+2026D67+TenMuE_0_200_pythia8_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal/step3.root')
)

process.gemDigis = cms.EDAnalyzer("MuonGEMDigisAnalyser",
    simTrackTag = cms.InputTag('g4SimHits'),
    simHitTag = cms.InputTag('g4SimHits',"MuonGEMHits"),
    gemDigiTag = cms.InputTag("simMuonGEMDigis"), 
    gemDigiSimLinkTag = cms.InputTag("simMuonGEMDigis","GEM"),
    gemSegmentTag = cms.InputTag("gemSegments"),
    min_pt = cms.double(5.0), # minimum pt for muon SimTrack
    min_quality = cms.double(0.6),
    min_num_layers = cms.uint32(4),
    # genParticleTag = cms.InputTag("genParticles"),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuonGEMDigisAnalysis.root')
)

process.p = cms.Path(process.gemDigis)
