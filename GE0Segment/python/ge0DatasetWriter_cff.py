import FWCore.ParameterSet.Config as cms
from MuonTriggering.GE0Segment.GE0DatasetWriter_cfi import GE0DatasetWriter

ge0DatasetWriter = cms.Sequence(GE0DatasetWriter)
