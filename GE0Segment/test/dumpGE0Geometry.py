import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9
process = cms.Process('ANALYSIS',Phase2C11I13M9)

GEOMETRY = "2026D76"

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended{}Reco_cff'.format(GEOMETRY))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('MuonTriggering.GE0Segment.GE0GeometryDumper_cfi')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source = cms.Source("EmptySource")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('GE0-{}.root'.format(GEOMETRY))
)

process.GE0 = process.GE0GeometryDumper.clone(
    treename=cms.untracked.string(GEOMETRY)
)
process.p = cms.Path(process.GE0)
