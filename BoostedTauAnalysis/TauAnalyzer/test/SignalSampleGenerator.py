# Auto generated configuration file
# using: 
# Revision: 1.381.2.6 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: test --conditions START53_V7A::All -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:7E33v2,RAW2DIGI,L1Reco,RECO --eventcontent AODSIM --fileout=Summer12_DR53X_NMSSMHiggs.root --number=1000 --mc --no_exec --beamspot Realistic8TeVCollision --datamix NODATAMIXER --pileup 2012_Summer_50ns_PoissonOOTPU
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2012_Summer_50ns_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_7E33v2_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#so the seed is randomized for each job
from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()
process.RandomNumberGeneratorService.saveFileName =  cms.untracked.string(
	"RandomEngineState_JOBNUM.log"
	)

from Configuration.Generator.PythiaUEZ2starSettings_cfi import *

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(300)
)

# Input source
process.source = cms.Source("EmptySource")
## process.source.firstRun = cms.untracked.uint32(JOBNUM)
process.source.firstRun = cms.untracked.uint32(1)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.6 $'),
    annotation = cms.untracked.string('test nevts:1000'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Summer12_DR53X_NMSSMHiggs_JOBNUM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
process.GlobalTag.globaltag = 'START53_V26::All' #latest tag for CMSSWv5.3.3 as of 9-Aug-13

process.generator = cms.EDFilter(
	"Pythia6GeneratorFilter",
	ExternalDecays = cms.PSet(
	Tauola = cms.untracked.PSet(
	UseTauolaPolarization = cms.bool(True),
	InputCards = cms.PSet(
	mdtau = cms.int32(0),
	pjak2 = cms.int32(0),
	pjak1 = cms.int32(0)
	)
        ),
	parameterSets = cms.vstring('Tauola')
	),
	maxEventsToPrint = cms.untracked.int32(0),
	pythiaPylistVerbosity = cms.untracked.int32(0),
	filterEfficiency = cms.untracked.double(1.0),
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	comEnergy = cms.double(8000.0),
	crossSection = cms.untracked.double(0.1913),
	UseExternalGenerators = cms.untracked.bool(True),
	PythiaParameters = cms.PSet(
	pythiaUESettingsBlock,
	processParameters = cms.vstring(

	#the full NMSSM can only be simulated in Pythia if an SLHA file is taken as input
	#here, we use the 2HDM as an approximation

	#setup
	'MSEL        = 0    !User defined processes',
## 	'MSUB(24)=1 !qqbar-->Z(h/H1(h/H2))',
## 	'MSUB(171)=1 !qqbar-->Z(H/H1(H/H2))',
	'MSUB(172)=1 !qqbar-->W+(H/H1(H/H2))',
##  	'MSUB(152)=1 !gg-->H/H1(H/H2)',
	'MSTP(4)=1 !modify h/H1 (h/H2) couplings',
## 	'IMSS(13)=1 !recognize NMSSM particles', #for reading in from an SLHA file

	#masses of the 3 neutral scalar Higgses in the NMSSM, lightest to heaviest
## 	'PMAS(C25,1)=115.0 !h/H1 mass',
## 	'PMAS(C35,1)=125.0 !H/H2 mass',
	'PMAS(C25,1)=125.0 !h/H2 mass',
	'PMAS(C35,1)=115.0 !H/H1 mass',	
	'PMAS(C45,1)=500.0 !H3 mass',
	
	#charged Higgs mass?

	#mass of the light pseudoscalar Higgs in the NMSSM
	'PMAS(C36,1)=9.0 !a mass',

	#mass of the CP-odd Higgs in the NMSSM
	'PMAS(C46,1)=500.0 !A2 mass',

	#mass of the extra singlino neutralino in the NMSSM
	#appears not to work unless an SLHA file is supplied and corresponding switches are set
## 	'PMAS(1000045,1)=100.0 !chi5 mass',

	#h/H1 (h/H2) couplings (these play little role for ZH production, mH > mh)
## 	'PARU(161)=10. !h/H1 (h/H2) coupling to d-type quarks',
## 	'PARU(162)=0.1 !h/H1 (h/H2) coupling to u-type quarks',
## 	'PARU(163)=10. !h/H1 (h/H2) coupling to leptons',
	'PARU(161)=1. !h/H1 (h/H2) coupling to d-type quarks', #=1 to be more SM-like
	'PARU(162)=1. !h/H1 (h/H2) coupling to u-type quarks', #=1 to be more SM-like
	'PARU(163)=1. !h/H1 (h/H2) coupling to leptons', #=1 to be more SM-like
	'PARU(164)=0.71 !h/H1 (h/H2) coupling to Z', #set to 1/sqrt(2) so PARU(164)^2 + PARU(174)^2 = 1
	'PARU(165)=0.71 !h/H1 (h/H2) coupling to W', #set to 1/sqrt(2) so PARU(165)^2 + PARU(175)^2 = 1
	'PARU(168)=0.0 !h/H1 (h/H2) coupling to H+/- in AA-->h/H1 (h/H2) loops',

	#H/H2 (H/H1) couplings
## 	'PARU(171)=10. !H/H2 (H/H1) coupling to d-type quarks',
## 	'PARU(172)=0.1 !H/H2 (H/H1) coupling to u-type quarks',
## 	'PARU(173)=10. !H/H2 (H/H1) coupling to leptons',
## 	'PARU(171)=0. !H/H2 (H/H1) coupling to d-type quarks', #force H-->aa
## 	'PARU(172)=0. !H/H2 (H/H1) coupling to u-type quarks', #force H-->aa
## 	'PARU(173)=0. !H/H2 (H/H1) coupling to leptons', #force H-->aa
	'PARU(171)=0.25 !H/H2 (H/H1) coupling to d-type quarks', #prefer H-->aa
	'PARU(172)=0.25 !H/H2 (H/H1) coupling to u-type quarks', #prefer H-->aa
	'PARU(173)=0.25 !H/H2 (H/H1) coupling to leptons', #prefer H-->aa
	'PARU(174)=0.71 !H/H2 (H/H1) coupling to Z', #set to 1/sqrt(2) so PARU(164)^2 + PARU(174)^2 = 1
	'PARU(175)=0.71 !H/H2 (H/H1) coupling to W', #set to 1/sqrt(2) so PARU(165)^2 + PARU(175)^2 = 1
	'PARU(176)=0.0 !H/H2 (H/H1) coupling to (h/H1 (h/H2))(h/H1 (h/H2))',
	'PARU(177)=4.0 !H/H2 (H/H1) coupling to aa', #what should this be?  want it dominant
	'PARU(178)=0.0 !H/H2 (H/H1) coupling to H+/- in AA-->H/H2 (H/H1) loops',

	#H3 couplings?

	#a couplings
	'PARU(181)=10. !a coupling to d-type quarks',
	'PARU(182)=0.1 !a coupling to u-type quarks',
	'PARU(183)=10. !a coupling to leptons',
	'PARU(184)=0. !a coupling to Z', #0 by default?
	'PARU(185)=0. !a coupling to W', #0 by default?
	'PARU(186)=0. !a coupling to Z(h/H1 (h/H2))',
	'PARU(187)=0. !a coupling to Z(H/H2 (H/H1))',
	'PARU(188)=0. !a coupling to Zprime(h/H1 (h/H2))',
	'PARU(189)=0. !a coupling to Zprime(H/H2 (H/H1))',
	'PARU(190)=0. !a coupling to H+/- in AA-->a loops',

	#tan beta (PARU(141))?

	#only allow leptonic W decays (this should affect the LO cross section reported by Pythia)
	#does this mess with status 2 Ws in hadronic tau decays?
	'24:ALLOFF',
	'24:ONIFANY 11 13 15',

	#only allow leptonic Z decays (this should affect the LO cross section reported by Pythia)
## 	'23:ALLOFF !turn on all Z decays',
## 	'23:ONIFANY 11 13 15 !only turn on Z-->ll decays',

	#default stuff that has to do with fragmentation
	'!MSTJ(11)=3 !something to do with fragmentation',
	'!MSTJ(27)=2 !something to do with fragmentation',
	'!PARJ(55)=-0.006 !something to do with fragmentation'
	),
	parameterSets = cms.vstring('pythiaUESettings', 'processParameters')
	)
	)
				 

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
