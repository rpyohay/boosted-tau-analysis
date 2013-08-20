import FWCore.ParameterSet.Config as cms

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15
MU_PDGID = 13

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'
    )
)

#analyze gen infor for Wh sample
process.analyzeGenInfo = cms.EDAnalyzer(
    'GenAnalyzer',
    outFileName = cms.string('/data1/yohay/Wh_gen_analysis_v2.root'),
    genParticleTag = cms.InputTag('genParticles'),
    PUTag = cms.InputTag('addPileupInfo'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet
    )

process.p = cms.Path(process.analyzeGenInfo)
