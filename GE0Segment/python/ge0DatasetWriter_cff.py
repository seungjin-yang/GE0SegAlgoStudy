import FWCore.ParameterSet.Config as cms
from GE0SegAlgoStudy.GE0Segment.GE0_cfi import GE0

ge0DatasetWriter = cms.Sequence(GE0)
