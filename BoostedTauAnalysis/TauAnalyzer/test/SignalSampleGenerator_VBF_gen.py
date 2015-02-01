# Auto generated configuration file
# using: 
# Revision: 1.381.2.6 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: test --conditions START53_V7A::All -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT:7E33v2,RAW2DIGI,L1Reco,RECO --eventcontent AODSIM --fileout=Summer12_DR53X_NMSSMHiggs.root --number=1000 --mc --no_exec --beamspot Realistic8TeVCollision --datamix NODATAMIXER --pileup 2012_Summer_50ns_PoissonOOTPU
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

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

#change the PU dataset for mix
process.mix.input.fileNames = cms.untracked.vstring(
	'/store/relval/CMSSW_5_3_6-START53_V14/RelValProdMinBias/GEN-SIM-RAW/v2/00000/4677049F-042A-E211-8525-0026189438E8.root',
	'/store/relval/CMSSW_5_3_6-START53_V14/RelValProdMinBias/GEN-SIM-RAW/v2/00000/52000D8A-032A-E211-BC94-00304867BFA8.root'
	)

#so the seed is randomized for each job
from IOMC.RandomEngine.RandomServiceHelper import  RandomNumberServiceHelper
randHelper =  RandomNumberServiceHelper(process.RandomNumberGeneratorService)
randHelper.populate()
process.RandomNumberGeneratorService.saveFileName =  cms.untracked.string(
	"RandomEngineState_JOBNUM.log"
	)

from Configuration.Generator.PythiaUEZ2starSettings_cfi import *

process.maxEvents = cms.untracked.PSet(
## 	input = cms.untracked.int32(1000)
	input = cms.untracked.int32(10)
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
##     eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
##     outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Summer12_DR53X_NMSSMHiggs_gen_JOBNUM.root'),
##     dataset = cms.untracked.PSet(
##         filterName = cms.untracked.string(''),
##         dataTier = cms.untracked.string('')
##     ),
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
	mdtau = cms.int32(212), #>=1 tau_mu
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
##	'MSUB(172)=1 !qqbar-->W+(H/H1(H/H2))',
##  	'MSUB(152)=1 !gg-->H/H1(H/H2)',
## 	'MSUB(174)=1 !WW fusion',
## 	'MSUB(173)=1 !ZZ fusion',
	'MSUB(124)=1 !WW fusion',
	'MSUB(123)=1 !ZZ fusion',
	'MSTP(4)=1 !modify h/H1 (h/H2) couplings',
## 	'IMSS(13)=1 !recognize NMSSM particles', #for reading in from an SLHA file

	#masses of the 3 neutral scalar Higgses in the NMSSM, lightest to heaviest
## 	'PMAS(C25,1)=115.0 !h/H1 mass',
## 	'PMAS(C35,1)=125.0 !H/H2 mass',
	'PMAS(C25,1)=500.0 !h/H2 mass',
	'PMAS(C35,1)=125.0 !H/H1 mass',	
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

	#allowed a1 decays
	'MDME(420,1)=0 !d               dbar',
	'MDME(421,1)=0 !u               ubar',
	'MDME(422,1)=0 !s               sbar',
	'MDME(423,1)=0 !c               cbar',
	'MDME(424,1)=0 !b               bbar',
	'MDME(425,1)=0 !t               tbar',
	'MDME(426,1)=0 !bprime          bprimebar',
	'MDME(427,1)=0 !tprime          tprimebar',
	'MDME(428,1)=0 !e-              e+',
	'MDME(429,1)=0 !mu-             mu+',
	'MDME(430,1)=1 !tau-            tau+',
	'MDME(431,1)=0 !tauprime-       tauprime+',
	'MDME(432,1)=0 !g               g',
	'MDME(433,1)=0 !gamma           gamma',
	'MDME(434,1)=0 !gamma           Z0',
	'MDME(435,1)=0 !Z0              Z0',
	'MDME(436,1)=0 !W+              W-',
	'MDME(437,1)=0 !Z0              h0',
	'MDME(438,1)=0 !h0              h0',
	'MDME(439,1)=0 !W+              H-',
	'MDME(440,1)=0 !H+              W-',
	'MDME(441,1)=0 !~chi_10         ~chi_10',
	'MDME(442,1)=0 !~chi_20         ~chi_10',
	'MDME(443,1)=0 !~chi_20         ~chi_20',
	'MDME(444,1)=0 !~chi_30         ~chi_10',
	'MDME(445,1)=0 !~chi_30         ~chi_20',
	'MDME(446,1)=0 !~chi_30         ~chi_30',
	'MDME(447,1)=0 !~chi_40         ~chi_10',
	'MDME(448,1)=0 !~chi_40         ~chi_20',
	'MDME(449,1)=0 !~chi_40         ~chi_30',
	'MDME(450,1)=0 !~chi_40         ~chi_40',
	'MDME(451,1)=0 !~chi_1+         ~chi_1-',
	'MDME(452,1)=0 !~chi_1+         ~chi_2-',
	'MDME(453,1)=0 !~chi_2+         ~chi_1-',
	'MDME(454,1)=0 !~chi_2+         ~chi_2-',
	'MDME(455,1)=0 !~d_L            ~d_Lbar',
	'MDME(456,1)=0 !~d_R            ~d_Rbar',
	'MDME(457,1)=0 !~d_L            ~d_Rbar',
	'MDME(458,1)=0 !~d_Lbar         ~d_R',
	'MDME(459,1)=0 !~u_L            ~u_Lbar',
	'MDME(460,1)=0 !~u_R            ~u_Rbar',
	'MDME(461,1)=0 !~u_L            ~u_Rbar',
	'MDME(462,1)=0 !~u_Lbar         ~u_R',
	'MDME(463,1)=0 !~s_L            ~s_Lbar',
	'MDME(464,1)=0 !~s_R            ~s_Rbar',
	'MDME(465,1)=0 !~s_L            ~s_Rbar',
	'MDME(466,1)=0 !~s_Lbar         ~s_R',
	'MDME(467,1)=0 !~c_L            ~c_Lbar',
	'MDME(468,1)=0 !~c_R            ~c_Rbar',
	'MDME(469,1)=0 !~c_L            ~c_Rbar',
	'MDME(470,1)=0 !~c_Lbar         ~c_R',
	'MDME(471,1)=0 !~b_1            ~b_1bar',
	'MDME(472,1)=0 !~b_2            ~b_2bar',
	'MDME(473,1)=0 !~b_1            ~b_2bar',
	'MDME(474,1)=0 !~b_1bar         ~b_2',
	'MDME(475,1)=0 !~t_1            ~t_1bar',
	'MDME(476,1)=0 !~t_2            ~t_2bar',
	'MDME(477,1)=0 !~t_1            ~t_2bar',
	'MDME(478,1)=0 !~t_1bar         ~t_2',
	'MDME(479,1)=0 !~e_L-           ~e_L+',
	'MDME(480,1)=0 !~e_R-           ~e_R+',
	'MDME(481,1)=0 !~e_L-           ~e_R+',
	'MDME(482,1)=0 !~e_L+           ~e_R-',
	'MDME(483,1)=0 !~nu_eL          ~nu_eLbar',
	'MDME(484,1)=0 !~nu_eR          ~nu_eRbar',
	'MDME(485,1)=0 !~nu_eL          ~nu_eRbar',
	'MDME(486,1)=0 !~nu_eLbar       ~nu_eR',
	'MDME(487,1)=0 !~mu_L-          ~mu_L+',
	'MDME(488,1)=0 !~mu_R-          ~mu_R+',
	'MDME(489,1)=0 !~mu_L-          ~mu_R+',
	'MDME(490,1)=0 !~mu_L+          ~mu_R-',
	'MDME(491,1)=0 !~nu_muL         ~nu_muLbar',
	'MDME(492,1)=0 !~nu_muR         ~nu_muRbar',
	'MDME(493,1)=0 !~nu_muL         ~nu_muRbar',
	'MDME(494,1)=0 !~nu_muLbar      ~nu_muR',
	'MDME(495,1)=0 !~tau_1-         ~tau_1+',
	'MDME(496,1)=0 !~tau_2-         ~tau_2+',
	'MDME(497,1)=0 !~tau_1-         ~tau_2+',
	'MDME(498,1)=0 !~tau_1+         ~tau_2-',
	'MDME(499,1)=0 !~nu_tauL        ~nu_tauLbar',
	'MDME(500,1)=0 !~nu_tauR        ~nu_tauRbar',
	'MDME(501,1)=0 !~nu_tauL        ~nu_tauRbar',
	'MDME(502,1)=0 !~nu_tauLbar     ~nu_tauR',

	#tan beta (PARU(141))?

##	#only allow leptonic W decays (this should affect the LO cross section reported by Pythia)
##	#does this mess with status 2 Ws in hadronic tau decays?
##	'24:ALLOFF',
##	'24:ONIFANY 11 13 15',

##	#only allow leptonic Z decays (this should affect the LO cross section reported by Pythia)
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

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
W_PDGID = 24
TAU_PDGID = 15
MU_PDGID = 13
NUMU_PDGID = 14
D_PDGID = 1
U_PDGID = 2
S_PDGID = 3
C_PDGID = 4
B_PDGID = 5
T_PDGID = 6
G_PDGID = 21
ANY_PDGID = 0

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

#tau hadronic decay types
TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

#no consideration of pT rank
ANY_PT_RANK = -1

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

#di-mu+X filter
process.genDiTauMuTauXSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    primaryTauDecayType = cms.uint32(TAU_MU),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(2),
    makeAllCollections = cms.bool(False)
    )

process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
## process.generation_step = cms.Path(process.pgen*process.genDiTauMuTauXSelector)
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
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step)
process.schedule.extend([process.endjob_step,process.AODSIMoutput_step])

# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
