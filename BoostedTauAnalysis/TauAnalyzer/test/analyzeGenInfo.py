import FWCore.ParameterSet.Config as cms

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15
MU_PDGID = 13

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.int32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_WH_files1-500.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files501-1000.root'
    )
)

#analyze gen infor for Wh sample
process.analyzeGenInfo = cms.EDAnalyzer(
    'GenAnalyzer',
    outFileName = cms.string('/data1/yohay/Wh_gen_analysis.root'),
    genParticleTag = cms.InputTag('genParticles'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet
    )

process.p = cms.Path(process.analyzeGenInfo)
