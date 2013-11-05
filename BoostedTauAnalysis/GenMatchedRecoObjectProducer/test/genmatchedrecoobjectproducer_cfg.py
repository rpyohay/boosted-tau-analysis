import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15

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
    'file:/data1/yohay/NMSSMHiggs_gg_skim.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.int32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0)
                                   )

#add a collection of gen muons from gen muonic tau decays, split by pT rank
process.genTauMuSelector0 = cms.EDFilter('GenObjectProducer',
                                         genParticleTag = cms.InputTag('genParticles'),
                                         absMatchPDGID = cms.uint32(TAU_PDGID),
                                         genTauDecayIDPSet = commonGenTauDecayIDPSet,
                                         primaryTauDecayType = cms.uint32(TAU_MU),
                                         sisterTauDecayType = cms.uint32(TAU_ALL),
                                         primaryTauPTRank = cms.int32(0),
                                         primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
                                         sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
                                         primaryTauAbsEtaMax = cms.double(-1.0),
                                         countSister = cms.bool(False),
                                         applyPTCuts = cms.bool(False),
                                         countKShort = cms.bool(True),
                                         minNumGenObjectsToPassFilter = cms.uint32(0),
                                         makeAllCollections = cms.bool(True)
                                         )
process.genTauMuSelector1 = process.genTauMuSelector0.clone()
process.genTauMuSelector1.primaryTauPTRank = cms.int32(1)
process.genTauMuSelector2 = process.genTauMuSelector0.clone()
process.genTauMuSelector2.primaryTauPTRank = cms.int32(2)
process.genTauMuSelector3 = process.genTauMuSelector0.clone()
process.genTauMuSelector3.primaryTauPTRank = cms.int32(3)

#add a collection of gen 1-prong taus, split by pT rank
process.genTau1ProngSelector0 = process.genTauMuSelector0.clone()
process.genTau1ProngSelector0.primaryTauDecayType = cms.uint32(TAU_HAD)
process.genTau1ProngSelector0.primaryTauHadronicDecayType = cms.int32(TAU_1PRONG_0NEUTRAL)
process.genTau1ProngSelector1 = process.genTau1ProngSelector0.clone()
process.genTau1ProngSelector1.primaryTauPTRank = cms.int32(1)
process.genTau1ProngSelector2 = process.genTau1ProngSelector0.clone()
process.genTau1ProngSelector2.primaryTauPTRank = cms.int32(2)
process.genTau1ProngSelector3 = process.genTau1ProngSelector0.clone()
process.genTau1ProngSelector3.primaryTauPTRank = cms.int32(3)

#add a collection of gen 1-prong + 1 pi0 taus, split by pT rank
process.genTau1Prong1Pi0Selector0 = process.genTau1ProngSelector0.clone()
process.genTau1Prong1Pi0Selector0.primaryTauHadronicDecayType = cms.int32(TAU_1PRONG_1NEUTRAL)
process.genTau1Prong1Pi0Selector1 = process.genTau1Prong1Pi0Selector0.clone()
process.genTau1Prong1Pi0Selector1.primaryTauPTRank = cms.int32(1)
process.genTau1Prong1Pi0Selector2 = process.genTau1Prong1Pi0Selector0.clone()
process.genTau1Prong1Pi0Selector2.primaryTauPTRank = cms.int32(2)
process.genTau1Prong1Pi0Selector3 = process.genTau1Prong1Pi0Selector0.clone()
process.genTau1Prong1Pi0Selector3.primaryTauPTRank = cms.int32(3)

#add a collection of gen 1-prong + 2 pi0 taus, split by pT rank
process.genTau1Prong2Pi0Selector0 = process.genTau1ProngSelector0.clone()
process.genTau1Prong2Pi0Selector0.primaryTauHadronicDecayType = cms.int32(TAU_1PRONG_2NEUTRAL)
process.genTau1Prong2Pi0Selector1 = process.genTau1Prong2Pi0Selector0.clone()
process.genTau1Prong2Pi0Selector1.primaryTauPTRank = cms.int32(1)
process.genTau1Prong2Pi0Selector2 = process.genTau1Prong2Pi0Selector0.clone()
process.genTau1Prong2Pi0Selector2.primaryTauPTRank = cms.int32(2)
process.genTau1Prong2Pi0Selector3 = process.genTau1Prong2Pi0Selector0.clone()
process.genTau1Prong2Pi0Selector3.primaryTauPTRank = cms.int32(3)

#add a collection of gen 3-prong taus, split by pT rank
process.genTau3ProngSelector0 = process.genTau1ProngSelector0.clone()
process.genTau3ProngSelector0.primaryTauHadronicDecayType = cms.int32(TAU_3PRONG_0NEUTRAL)
process.genTau3ProngSelector1 = process.genTau3ProngSelector0.clone()
process.genTau3ProngSelector1.primaryTauPTRank = cms.int32(1)
process.genTau3ProngSelector2 = process.genTau3ProngSelector0.clone()
process.genTau3ProngSelector2.primaryTauPTRank = cms.int32(2)
process.genTau3ProngSelector3 = process.genTau3ProngSelector0.clone()
process.genTau3ProngSelector3.primaryTauPTRank = cms.int32(3)

#add a collection of any gen tau->mu from NMSSM a decay
process.genTauMuSelectorAny = cms.EDFilter('GenObjectProducer',
                                           genParticleTag = cms.InputTag('genParticles'),
                                           absMatchPDGID = cms.uint32(TAU_PDGID),
                                           genTauDecayIDPSet = commonGenTauDecayIDPSet,
                                           primaryTauDecayType = cms.uint32(TAU_MU),
                                           sisterTauDecayType = cms.uint32(TAU_MU),
                                           primaryTauPTRank = cms.int32(ANY_PT_RANK),
                                           primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
                                           sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
                                           primaryTauAbsEtaMax = cms.double(-1.0),
                                           countSister = cms.bool(True),
                                           applyPTCuts = cms.bool(False),
                                           countKShort = cms.bool(True),
                                           minNumGenObjectsToPassFilter = cms.uint32(0),
                                           makeAllCollections = cms.bool(False)
                                           )

#add a collection of any gen tau from NMSSM a decay
process.genTauSelector = process.genTauMuSelectorAny.clone()
process.genTauSelector.primaryTauDecayType = cms.uint32(TAU_ALL)
process.genTauSelector.sisterTauDecayType = cms.uint32(TAU_ALL)

#add a collection of any gen tau from NMSSM a decay with a hadronic sister
process.genBoostedDiTauHadSelector = process.genTauMuSelectorAny.clone()
process.genBoostedDiTauHadSelector.primaryTauDecayType = cms.uint32(TAU_ALL)
process.genBoostedDiTauHadSelector.countSister = cms.bool(False)
process.genBoostedDiTauHadSelector.sisterTauDecayType = cms.uint32(TAU_HAD)
process.genBoostedDiTauHadSelector.minNumGenObjectsToPassFilter = cms.uint32(2)

## #add a collection of tight muons
## #(cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId) to the event
## process.tightMuonSelector = cms.EDFilter('CustomMuonSelector',
##                                          muonTag = cms.InputTag('muons'),
##                                          vtxTag = cms.InputTag('offlinePrimaryVertices'),
##                                          PFIsoMax = cms.double(-1.0),
##                                          etaMax = cms.double(-1.0)
##                                          )

#add collections of gen muonic tau decays matched to trigger muons
## process.Mu17Mu8MatchedGenMuonSelector = cms.EDProducer(
##     'trgMatchedGenParticleProducer',
##     InputProducer = cms.InputTag('genMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_Mu17_Mu8_v16', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8")
##     )
## process.IsoMu20eta2p1CentralPFNoPUJet30MatchedGenMuonSelector = cms.EDProducer(
##     'trgMatchedGenParticleProducer',
##     InputProducer = cms.InputTag('genMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_IsoMu20_eta2p1_CentralPFNoPUJet30_v2', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag(
##     "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f20L3crIsoFiltered10"
##     )
##     )
## process.DoubleMu14Mass8PFMET40v3MatchedGenMuonSelector = cms.EDProducer(
##     'trgMatchedGenParticleProducer',
##     InputProducer = cms.InputTag('genMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_DoubleMu14_Mass8_PFMET40_v3', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag(
##     "hltL1DoubleMu10MuOpenORDoubleMu103p5L3DiMu14Mass8Filtered"
##     )
##     )
## process.IsoMu24eta2p1v3MatchedGenMuonSelector = cms.EDProducer(
##     'trgMatchedGenParticleProducer',
##     InputProducer = cms.InputTag('genMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_IsoMu24_eta2p1_v11', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag(
##     "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10"
##     )
##     )
## process.Mu50eta2p1v6MatchedGenMuonSelector = cms.EDProducer(
##     'trgMatchedGenParticleProducer',
##     InputProducer = cms.InputTag('genMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_Mu50_eta2p1_v6', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag(
##     "hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered50Q"
##     )
##     )
process.IsoMu15eta2p1L1ETM20v3MatchedGenMuonSelector = cms.EDProducer(
    'trgMatchedGenParticleProducer',
    InputProducer = cms.InputTag('genTauMuSelectorAny'),
    hltTags = cms.VInputTag(cms.InputTag('HLT_IsoMu15_eta2p1_L1ETM20_v3', '', 'HLT')),
    isTriggerFilter = cms.untracked.bool(True),
    HLTSubFilters = cms.untracked.VInputTag(
    "hltL3crIsoL1sMu12Eta2p1L1f0L2f12QL3f15QL3crIsoFiltered10"
    )
    )

## #add a collection of tight muons matched to 8 GeV trigger muons
## process.triggerMatchedMuonSelector = cms.EDProducer(
##     'trgMatchedMuonProducer',
##     InputProducer = cms.InputTag('tightMuonSelector'),
##     hltTags = cms.VInputTag(cms.InputTag('HLT_Mu17_Mu8_v16', '', 'HLT')),
##     isTriggerFilter = cms.untracked.bool(True),
##     HLTSubFilters = cms.untracked.VInputTag("hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8")
##     )

#add a collection of jets matched to gen boosted di-tau pairs where one of the pair is a hadronic
#tau
process.genMatchedJetSelector = cms.EDFilter(
    'GenMatchedJetProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genBoostedDiTauHadSelector'),
    recoObjTag = cms.InputTag('ak5PFJets'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(True),
    useGenObjPTRank = cms.bool(False),
    nOutputColls = cms.uint32(2),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(0)
    )

#analyze
process.genTriggerAnalyzer = cms.EDAnalyzer(
    'DecayModePTRankAnalyzer',
    outFileName = cms.string('/data1/yohay/NMSSMHiggs_gg_trigger_analysis_v2.root'),
    tauMuInputTags = cms.VInputTag(
    cms.InputTag('genTauMuSelector0', 'decayModeMuPTRank0', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayModeMuPTRank1', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayModeMuPTRank2', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayModeMuPTRank3', 'OWNPARTICLES')),
    tau1ProngInputTags = cms.VInputTag(
    cms.InputTag('genTauMuSelector0', 'decayMode1ProngPTRank0', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1ProngPTRank1', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1ProngPTRank2', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1ProngPTRank3', 'OWNPARTICLES')),
    tau1Prong1Pi0InputTags = cms.VInputTag(
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong1Pi0PTRank0', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong1Pi0PTRank1', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong1Pi0PTRank2', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong1Pi0PTRank3', 'OWNPARTICLES')),
    tau1Prong2Pi0InputTags = cms.VInputTag(
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong2Pi0PTRank0', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong2Pi0PTRank1', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong2Pi0PTRank2', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode1Prong2Pi0PTRank3', 'OWNPARTICLES')),
    tau3ProngInputTags = cms.VInputTag(
    cms.InputTag('genTauMuSelector0', 'decayMode3ProngPTRank0', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode3ProngPTRank1', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode3ProngPTRank2', 'OWNPARTICLES'),
    cms.InputTag('genTauMuSelector0', 'decayMode3ProngPTRank3', 'OWNPARTICLES')),
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

#analyze
process.jetAnalyzer = cms.EDAnalyzer(
    'JetAnalyzer',
    outFileName = cms.string('/data1/yohay/jet_analysis.root'),
    jetTags = cms.VInputTag(cms.InputTag('genMatchedJetSelector', 'coll0', 'OWNPARTICLES'),
                            cms.InputTag('genMatchedJetSelector', 'coll1', 'OWNPARTICLES')),
    METTag = cms.InputTag('pfMet'),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}')
    )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/data1/yohay/debug_EDM.root'),
                               SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
                               )

  
## process.p = cms.Path(process.genMuonSelector*process.Mu17Mu8MatchedGenMuonSelector*
##                      process.IsoMu20eta2p1CentralPFNoPUJet30MatchedGenMuonSelector*
##                      process.DoubleMu14Mass8PFMET40v3MatchedGenMuonSelector*
##                      process.IsoMu24eta2p1v3MatchedGenMuonSelector*
##                      process.Mu50eta2p1v6MatchedGenMuonSelector*
##                      process.genMuonTriggerAnalyzer)

## process.p = cms.Path(process.genTauMuSelector0 + process.genTauMuSelector1 +
##                      process.genTauMuSelector2 + process.genTauMuSelector3 +
##                      process.genTau1ProngSelector0 + process.genTau1ProngSelector1 +
##                      process.genTau1ProngSelector2 + process.genTau1ProngSelector3 +
##                      process.genTau1Prong1Pi0Selector0 + process.genTau1Prong1Pi0Selector1 +
##                      process.genTau1Prong1Pi0Selector2 + process.genTau1Prong1Pi0Selector3 +
##                      process.genTau1Prong2Pi0Selector0 + process.genTau1Prong2Pi0Selector1 +
##                      process.genTau1Prong2Pi0Selector2 + process.genTau1Prong2Pi0Selector3 +
##                      process.genTau3ProngSelector0 + process.genTau3ProngSelector1 +
##                      process.genTau3ProngSelector2 + process.genTau3ProngSelector3 +
##                      process.genTriggerAnalyzer)

## process.p = cms.Path(process.genTauSelector + process.genTauMuSelectorAny +
##                      process.genTauMuSelector0 + process.genBoostedDiTauHadSelector + 
##                      process.IsoMu15eta2p1L1ETM20v3MatchedGenMuonSelector +
##                      process.genMatchedJetSelector + process.genTriggerAnalyzer)

process.p = cms.Path(process.genBoostedDiTauHadSelector)

## process.e = cms.EndPath(process.out)
