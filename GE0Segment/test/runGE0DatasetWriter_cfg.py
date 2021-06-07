# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step4 --conditions auto:phase2_realistic_T21 --era Phase2C11I13M9 --geometry Extended2026D76 --filein file:step3.root --fileout file:step4.root --python_filename step_4_cfg.py --step NONE --no_exec --number -1
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C11I13M9_cff import Phase2C11I13M9

process = cms.Process('NONE',Phase2C11I13M9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D76Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('MuonTriggering.GE0Segment.ge0DatasetWriter_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/T2_KR_KISTI/store/user/seyang/GE0_EightMu_Pt_5_200_Eta_1p7_3p2_noPU/CMSSW_12_0_0_pre1_step_3_GEN-SIM-RECO/210515_065103/0001/step3_1300.root'),
    secondaryFileNames = cms.untracked.vstring()
)

# Production Info

# Output definition

# Additional output definition
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('step4.root')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# Path and EndPath definitions
process.GE0 = cms.Path(process.ge0DatasetWriter)

# Schedule definition
process.schedule = cms.Schedule(process.GE0)
