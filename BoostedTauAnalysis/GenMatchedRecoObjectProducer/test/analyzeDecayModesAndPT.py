import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15
MU_PDGID = 13

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

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_gg_files1-500.root',
    'file:/data1/yohay/NMSSMHiggs_gg_files501-1000.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.int32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

#produce gen tau collections
process.genTauSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    primaryTauDecayType = cms.uint32(TAU_ALL),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(0),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    countSister = cms.bool(True),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(0),
    makeAllCollections = cms.bool(True)
    )

#analyze
process.genTauAnalyzer = cms.EDAnalyzer(
    'DecayModePTRankAnalyzer',
    outFileName = cms.string('/data1/yohay/NMSSMHiggs_gg_decayModePTRank_analysis.root'),
    tauMuInputTags = cms.VInputTag(
    cms.InputTag('genTauSelector', 'decayModeMuPTRank0', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayModeMuPTRank1', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayModeMuPTRank2', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayModeMuPTRank3', 'ANALYSIS')),
    tau1ProngInputTags = cms.VInputTag(
    cms.InputTag('genTauSelector', 'decayMode1ProngPTRank0', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1ProngPTRank1', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1ProngPTRank2', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1ProngPTRank3', 'ANALYSIS')),
    tau1Prong1Pi0InputTags = cms.VInputTag(
    cms.InputTag('genTauSelector', 'decayMode1Prong1Pi0PTRank0', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong1Pi0PTRank1', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong1Pi0PTRank2', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong1Pi0PTRank3', 'ANALYSIS')),
    tau1Prong2Pi0InputTags = cms.VInputTag(
    cms.InputTag('genTauSelector', 'decayMode1Prong2Pi0PTRank0', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong2Pi0PTRank1', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong2Pi0PTRank2', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode1Prong2Pi0PTRank3', 'ANALYSIS')),
    tau3ProngInputTags = cms.VInputTag(
    cms.InputTag('genTauSelector', 'decayMode3ProngPTRank0', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode3ProngPTRank1', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode3ProngPTRank2', 'ANALYSIS'),
    cms.InputTag('genTauSelector', 'decayMode3ProngPTRank3', 'ANALYSIS')),
    genParticleTag = cms.InputTag('genParticles'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    decayModeColors = cms.vuint32(1, 2, 4, 6, 8),
    decayModeStyles = cms.vuint32(20, 21, 22, 23, 24),
    decayModeEntries = cms.vstring('#tau_{#mu}', '#tau_{had}, 1 prong',
                                   '#tau_{had}, 1 prong + 1 #pi^{0}',
                                   '#tau_{had}, 1 prong + 2 #pi^{0}', '#tau_{had}, 3 prong')
    )

#path
process.p = cms.Path(process.genTauSelector*process.genTauAnalyzer)
