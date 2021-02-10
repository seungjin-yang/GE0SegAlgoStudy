# Auto generated configuration file
# NOTE I modify a script created by: 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: TenMuExtendedE_0_100_pythia8_cfi --conditions auto:phase2_realistic -n 10 --era Phase2 --eventcontent FEVTDEBUG --relval 10000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --geometry Extended2023D17 --fileout file:step1.root
# 
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2_cff import Phase2

process = cms.Process('SIM',Phase2)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D17_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedHLLHC_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('TenMuExtendedPt_0_100_pythia8_cfi nevts:10'), # NOTE E to Pt
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step1.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

# NOTE I can't find TenMu with flat pt distribution using `runTheMatrix.py -w upgrade -n`
# So I made this cfi based on `TenMuExtendedE_0_100_pythia8_cfi_GEN_SIM.py`
# You can find the way to create original cfi on the top of this script.
# https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_1/GeneratorInterface/Pythia8Interface/plugins/Py8PtGun.cc
process.generator = cms.EDFilter("Pythia8PtGun", # NOTE Pythia8EGun to Pythia8PtGun
    PGunParameters = cms.PSet(
        AddAntiParticle = cms.bool(True),
        MaxPt = cms.double(100.0), # NOTE MaxE to MaxPt
        MaxEta = cms.double(2.7), # NOTE target to ME0
        MaxPhi = cms.double(3.14159265359),
        MinPt = cms.double(5.0), # NOTE MaxE to MaxPt, we want muons with energy greater than 5 GeV
        MinEta = cms.double(2.1), # NOTE target to ME0
        MinPhi = cms.double(-3.14159265359),
        ParticleID = cms.vint32(-13, -13, -13, -13, -13)
    ),
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring()
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    psethack = cms.string('Ten mu pt 0 to 100') # NOTE e to pt
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path).insert(0, process.generator)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
