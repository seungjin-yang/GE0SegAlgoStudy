import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("me0", eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D28Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D28_cff')
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
    #fileNames = cms.untracked.vstring('file:/scratch/me0/10mu/step2_10mu.root')
    fileNames = cms.untracked.vstring('file:/scratch/me0/minbias/step2_minbias.root')
)

process.me0Digis = cms.EDAnalyzer("MuonME0DigisAnalyser",
    me0DigiToken = cms.InputTag("simMuonME0Digis"), 
    me0DigiSimLinkToken = cms.InputTag("simMuonME0Digis","ME0"),
    simTrackCollection = cms.InputTag('g4SimHits'),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuonME0DigisAnalysis.root')
)

process.p = cms.Path(process.me0Digis)
