import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Dump", eras.Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D24Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D24_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag


# process.GlobalTag = GlobalTag(process.GlobalTag, '90X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

# fileNames = cms.untracked.vstring('file:/xrootd_UOS/store/user/yekang/CRAB_PrivateMC/me0_tenMu_GEN_SIM_DIGI_RECO/190315_092951/0000/tenMu_GEN-SIM-DIGI-RECO_1.root')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:tenMu_GEN-SIM-DIGI-RECO_1.root')
)

process.dumper = cms.EDAnalyzer("MuonME0DigisAnalyser",
    me0DigiToken = cms.InputTag("simMuonME0Digis"), 
    me0DigiSimLinkToken = cms.InputTag("simMuonME0Digis","ME0"),
    simTrackCollection = cms.InputTag('g4SimHits'),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('MuonME0DigisAnalysis.root')
)

process.p = cms.Path(process.dumper)
