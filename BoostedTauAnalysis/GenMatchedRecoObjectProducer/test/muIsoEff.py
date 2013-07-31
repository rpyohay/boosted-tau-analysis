import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

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
    'file:/data1/yohay/FILE_NAME.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#define a parameter set to be passed to all modules that utilize GenTauDecayID for A1 decays
commonGenTauDecayIDPSetSignal = cms.PSet(momPDGID = cms.int32(A_PDGID),
                                         chargedHadronPTMin = cms.double(0.0),
                                         neutralHadronPTMin = cms.double(0.0),
                                         chargedLeptonPTMin = cms.double(0.0),
                                         totalPTMin = cms.double(0.0)
                                         )

#define a parameter set to be passed to all modules that utilize GenTauDecayID for Z decays
commonGenTauDecayIDPSetZMuMu = commonGenTauDecayIDPSetSignal.clone()
## commonGenTauDecayIDPSetZMuMu.momPDGID = cms.int32(Z_PDGID)
commonGenTauDecayIDPSetZMuMu.momPDGID = cms.int32(24)

#produce gen tau-->mu + tau-->had collection in A1 decays, only saving the muon
process.genTauMuTauHadSelectorSignal = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSetSignal,
    primaryTauDecayType = cms.uint32(TAU_MU),
    sisterTauDecayType = cms.uint32(TAU_HAD),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#produce gen Z-->mumu collection in A1 decays, saving both muons
process.genMuMuSelectorZMuMu = process.genTauMuTauHadSelectorSignal.clone()
process.genMuMuSelectorZMuMu.absMatchPDGID = cms.uint32(MU_PDGID)
process.genMuMuSelectorZMuMu.genTauDecayIDPSet = commonGenTauDecayIDPSetZMuMu
process.genMuMuSelectorZMuMu.primaryTauDecayType = cms.uint32(TAU_ALL)
process.genMuMuSelectorZMuMu.sisterTauDecayType = cms.uint32(TAU_ALL)
## process.genMuMuSelectorZMuMu.countSister = cms.bool(True)
process.genMuMuSelectorZMuMu.countSister = cms.bool(False)
## process.genMuMuSelectorZMuMu.minNumGenObjectsToPassFilter = cms.uint32(2)
process.genMuMuSelectorZMuMu.minNumGenObjectsToPassFilter = cms.uint32(1)

#add a collection of tight muons in |eta| < 2.1 failing tight PF isolation
process.tightNonIsoMuonSelector = cms.EDFilter('CustomMuonSelector',
                                               muonTag = cms.InputTag('muons'),
                                               vtxTag = cms.InputTag('offlinePrimaryVertices'),
##                                                PFIsoMax = cms.double(0.12),
                                               PFIsoMax = cms.double(-1.0),
                                               detectorIsoMax = cms.double(0.1),
                                               PUSubtractionCoeff = cms.double(0.5),
                                               usePFIso = cms.bool(True),
                                               passIso = cms.bool(False),
                                               etaMax = cms.double(2.1),
                                               minNumObjsToPassFilter = cms.uint32(0)
                                               )

#add a collection of tight muons in |eta| < 2.1 passing tight detector isolation
process.tightDetectorIsoMuonSelector = process.tightNonIsoMuonSelector.clone()
process.tightDetectorIsoMuonSelector.usePFIso = cms.bool(False)
process.tightDetectorIsoMuonSelector.passIso = cms.bool(True)

#add a collection of tight muons in |eta| < 2.1 passing tight PF isolation w/o PU subtraction
process.tightPFIsoNoPUSubtractionMuonSelector = process.tightNonIsoMuonSelector.clone()
process.tightPFIsoNoPUSubtractionMuonSelector.PUSubtractionCoeff = cms.double(0.0)
process.tightPFIsoNoPUSubtractionMuonSelector.passIso = cms.bool(True)

#add a collection of tight muons in |eta| < 2.1
process.tightMuonSelector = process.tightNonIsoMuonSelector.clone()
process.tightMuonSelector.PFIsoMax = cms.double(-1.0)

#add a collection of tight muons in |eta| < 2.1 failing tight PF isolation matched to gen muons
#from tau-->mu + tau-->had A1 decays
process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal = cms.EDFilter(
    'GenMatchedMuonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauMuTauHadSelectorSignal'),
    recoObjTag = cms.InputTag('tightNonIsoMuonSelector'),
    genTauDecayIDPSet = commonGenTauDecayIDPSetSignal,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(0)
    )

#add a collection of tight muons in |eta| < 2.1 failing tight PF isolation matched to gen muons
#from Z-->mumu decays
process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu = process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal.clone()
process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu.selectedGenParticleTag = cms.InputTag(
    'genMuMuSelectorZMuMu'
    )
process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu.genTauDecayIDPSet = commonGenTauDecayIDPSetZMuMu

#add a collection of tight muons in |eta| < 2.1 passing tight detector isolation matched to gen
#muons from tau-->mu + tau-->had A1 decays
process.genTauMuTauHadMatchedTightDetectorIsoMuonSelectorSignal = process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal.clone()
process.genTauMuTauHadMatchedTightDetectorIsoMuonSelectorSignal.recoObjTag = cms.InputTag(
    'tightDetectorIsoMuonSelector'
    )

#add a collection of tight muons in |eta| < 2.1 passing tight detector isolation matched to gen
#muons from Z-->mumu decays
process.genMuMuMatchedTightDetectorIsoMuonSelectorZMuMu = process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu.clone()
process.genMuMuMatchedTightDetectorIsoMuonSelectorZMuMu.recoObjTag = cms.InputTag(
    'tightDetectorIsoMuonSelector'
    )

#add a collection of tight muons in |eta| < 2.1 passing tight PF isolation w/o PU subtraction
#matched to gen muons from tau-->mu + tau-->had A1 decays
process.genTauMuTauHadMatchedTightPFIsoNoPUSubtractionMuonSelectorSignal = process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal.clone()
process.genTauMuTauHadMatchedTightPFIsoNoPUSubtractionMuonSelectorSignal.recoObjTag = cms.InputTag(
    'tightPFIsoNoPUSubtractionMuonSelector'
    )

#add a collection of tight muons in |eta| < 2.1 passing tight PF isolation w/o PU subtraction
#matched to gen muons from Z-->mumu decays
process.genMuMuMatchedTightPFIsoNoPUSubtractionMuonSelectorZMuMu = process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu.clone()
process.genMuMuMatchedTightPFIsoNoPUSubtractionMuonSelectorZMuMu.recoObjTag = cms.InputTag(
    'tightPFIsoNoPUSubtractionMuonSelector'
    )

#add a collection of tight muons in |eta| < 2.1 from tau-->mu + tau-->had A1 decays
process.genTauMuTauHadMatchedTightMuonSelectorSignal = process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal.clone()
process.genTauMuTauHadMatchedTightMuonSelectorSignal.recoObjTag = cms.InputTag('tightMuonSelector')

#add a collection of tight muons in |eta| < 2.1 from Z-->mumu decays
process.genMuMuMatchedTightMuonSelectorZMuMu = process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu.clone()
process.genMuMuMatchedTightMuonSelectorZMuMu.recoObjTag = cms.InputTag('tightMuonSelector')

#produce reco muons failing tight PF isolation matched to gen tau-->mu decays matched to
#IsoMu24_eta2p1 trigger objects from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal = cms.EDProducer(
    'trgMatchedMuonProducer',
    InputProducer = cms.InputTag('genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal'),
    hltTags = cms.VInputTag(cms.InputTag('HLT_IsoMu24_eta2p1_v11', '', 'HLT')),
    isTriggerFilter = cms.untracked.bool(True),
    HLTSubFilters = cms.untracked.VInputTag(
    'hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10'
    )
    )

#produce reco muons failing tight PF isolation matched to gen muons matched to IsoMu24_eta2p1
#trigger objects from Z-->mumu decays
process.IsoMu24eta2p1TightNonIsoMatchedGenMuMuSelectorZMuMu = process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal.clone()
process.IsoMu24eta2p1TightNonIsoMatchedGenMuMuSelectorZMuMu.InputProducer = cms.InputTag(
    'genMuMuMatchedTightNonIsoMuonSelectorZMuMu'
    )

#produce reco muons passing tight detector isolation matched to gen tau-->mu decays matched to
#IsoMu24_eta2p1 trigger objects from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightDetectorIsoMatchedGenTauMuTauHadSelectorSignal = process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal.clone()
process.IsoMu24eta2p1TightDetectorIsoMatchedGenTauMuTauHadSelectorSignal.InputProducer = cms.InputTag('genTauMuTauHadMatchedTightDetectorIsoMuonSelectorSignal')

#produce reco muons passing tight detector isolation matched to gen muons matched to IsoMu24_eta2p1
#trigger objects from Z-->mumu decays
process.IsoMu24eta2p1TightDetectorIsoMatchedGenMuMuSelectorZMuMu = process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal.clone()
process.IsoMu24eta2p1TightDetectorIsoMatchedGenMuMuSelectorZMuMu.InputProducer = cms.InputTag(
    'genMuMuMatchedTightDetectorIsoMuonSelectorZMuMu'
    )

#produce reco muons passing tight PF isolation w/o PU subtraction matched to gen tau-->mu decays
#matched to IsoMu24_eta2p1 trigger objects from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenTauMuTauHadSelectorSignal = process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal.clone()
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenTauMuTauHadSelectorSignal.InputProducer = cms.InputTag('genTauMuTauHadMatchedTightPFIsoNoPUSubtractionMuonSelectorSignal')

#produce reco muons passing tight PF isolation w/o PU subtraction matched to gen muons matched to
#IsoMu24_eta2p1 trigger objects from Z-->mumu decays
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenMuMuSelectorZMuMu = process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal.clone()
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenMuMuSelectorZMuMu.InputProducer = cms.InputTag('genMuMuMatchedTightPFIsoNoPUSubtractionMuonSelectorZMuMu')

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 failing tight PF isolation
#matched to gen muons from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightNonIsoEffSignal = cms.EDAnalyzer(
    'MuonEfficiencyAnalyzer',
    outFileName = cms.string('/data1/yohay/NMSSMHiggs_gg_muLeg_tauMuTauHad_IsoMu24eta2p1.root'),
    denominatorTag = cms.InputTag('genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal'),
    numeratorTag = cms.InputTag('IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal'),
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

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 failing tight PF isolation
#matched to gen muons from Z-->mumu decays
process.IsoMu24eta2p1TightNonIsoEffZMuMu = cms.EDAnalyzer(
##     'MuonTagAndProbeEfficiencyAnalyzer',
    'MuonEfficiencyAnalyzer',
##     outFileName = cms.string('/data1/yohay/ZMuMu_IsoMu24eta2p1.root'),
    outFileName = cms.string('/data1/yohay/NMSSMHiggs_WH_IsoMu24eta2p1.root'),
##     tagTagTag = cms.InputTag('IsoMu24eta2p1TightNonIsoMatchedGenMuMuSelectorZMuMu'),
    numeratorTag = cms.InputTag('IsoMu24eta2p1TightNonIsoMatchedGenMuMuSelectorZMuMu'),
##     tagFailTag = cms.InputTag('genMuMuMatchedTightNonIsoMuonSelectorZMuMu'),
    denominatorTag = cms.InputTag('genMuMuMatchedTightNonIsoMuonSelectorZMuMu'),
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

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 passing tight detector
#isolation matched to gen muons from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightDetectorIsoEffSignal = process.IsoMu24eta2p1TightNonIsoEffSignal.clone()
process.IsoMu24eta2p1TightDetectorIsoEffSignal.outFileName = cms.string(
    '/data1/yohay/NMSSMHiggs_gg_muLeg_tightDetectorIso_tauMuTauHad_IsoMu24eta2p1.root'
    )
process.IsoMu24eta2p1TightDetectorIsoEffSignal.denominatorTag = cms.InputTag(
    'genTauMuTauHadMatchedTightDetectorIsoMuonSelectorSignal'
    )
process.IsoMu24eta2p1TightDetectorIsoEffSignal.numeratorTag = cms.InputTag(
    'IsoMu24eta2p1TightDetectorIsoMatchedGenTauMuTauHadSelectorSignal'
    )

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 passing tight detector
#isolation matched to gen muons from Z-->mumu decays
process.IsoMu24eta2p1TightDetectorIsoEffZMuMu = process.IsoMu24eta2p1TightNonIsoEffZMuMu.clone()
process.IsoMu24eta2p1TightDetectorIsoEffZMuMu.outFileName = cms.string(
    '/data1/yohay/ZMuMu_tightDetectorIso_IsoMu24eta2p1.root'
    )
process.IsoMu24eta2p1TightDetectorIsoEffZMuMu.tagTagTag = cms.InputTag(
    'IsoMu24eta2p1TightDetectorIsoMatchedGenMuMuSelectorZMuMu'
    )
process.IsoMu24eta2p1TightDetectorIsoEffZMuMu.tagFailTag = cms.InputTag(
    'genMuMuMatchedTightDetectorIsoMuonSelectorZMuMu'
    )

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 passing tight PF isolation
#w/o PU subtraction matched to gen muons from tau-->mu + tau-->had A1 decays
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffSignal = process.IsoMu24eta2p1TightNonIsoEffSignal.clone()
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffSignal.outFileName = cms.string(
    '/data1/yohay/NMSSMHiggs_gg_muLeg_tightPFIsoNoPUSubtraction_tauMuTauHad_IsoMu24eta2p1.root'
    )
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffSignal.denominatorTag = cms.InputTag(
    'genTauMuTauHadMatchedTightPFIsoNoPUSubtractionMuonSelectorSignal'
    )
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffSignal.numeratorTag = cms.InputTag(
    'IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenTauMuTauHadSelectorSignal'
    )

#compute efficiency of IsoMu24_eta2p1 for tight muons in |eta| < 2.1 passing tight PF isolation
#w/o PU subtraction matched to gen muons from Z-->mumu decays
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffZMuMu = process.IsoMu24eta2p1TightNonIsoEffZMuMu.clone()
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffZMuMu.outFileName = cms.string(
    '/data1/yohay/ZMuMu_tightPFIsoNoPUSubtraction_IsoMu24eta2p1.root'
    )
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffZMuMu.tagTagTag = cms.InputTag(
    'IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenMuMuSelectorZMuMu'
    )
process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffZMuMu.tagFailTag = cms.InputTag(
    'genMuMuMatchedTightPFIsoNoPUSubtractionMuonSelectorZMuMu'
    )

#analyze tight muon isolation for signal A1-->tau(-->mu)tau(-->had) events
process.tightMuonIsoAnalyzerSignal = cms.EDAnalyzer(
    'IsolationAnalyzer',
    outFileName = cms.string('/data1/yohay/isoAnalysis_tightMuons_signal.root'),
    muonTag = cms.InputTag('genTauMuTauHadMatchedTightMuonSelectorSignal'),
    muonPFIsoPUSubtractionCoeff = cms.double(0.5)
    )

#analyze tight muon isolation for Z-->mumu events
process.tightMuonIsoAnalyzerZMuMu = process.tightMuonIsoAnalyzerSignal.clone()
process.tightMuonIsoAnalyzerZMuMu.outFileName = cms.string(
    '/data1/yohay/isoAnalysis_tightMuons_ZMuMu.root')
process.tightMuonIsoAnalyzerZMuMu.muonTag = cms.InputTag('genMuMuMatchedTightMuonSelectorZMuMu')

#sequence of modules to be run for the signal A1-->tau(-->mu)tau(-->had) analysis on non-isolated
#muons
process.tightNonIsoSignalSequence = cms.Sequence(
    process.genTauMuTauHadSelectorSignal*process.tightNonIsoMuonSelector*
    process.genTauMuTauHadMatchedTightNonIsoMuonSelectorSignal*
    process.IsoMu24eta2p1TightNonIsoMatchedGenTauMuTauHadSelectorSignal*
    process.IsoMu24eta2p1TightNonIsoEffSignal
    )

#sequence of modules to be run for the Z-->mumu analysis on non-isolated muons
process.tightNonIsoZMuMuSequence = cms.Sequence(
    process.genMuMuSelectorZMuMu*process.tightNonIsoMuonSelector*
    process.genMuMuMatchedTightNonIsoMuonSelectorZMuMu*
    process.IsoMu24eta2p1TightNonIsoMatchedGenMuMuSelectorZMuMu*
    process.IsoMu24eta2p1TightNonIsoEffZMuMu
    )

#sequence of modules to be run for the signal A1-->tau(-->mu)tau(-->had) analysis on detector
#isolated muons
process.tightDetectorIsoSignalSequence = cms.Sequence(
    process.genTauMuTauHadSelectorSignal*process.tightDetectorIsoMuonSelector*
    process.genTauMuTauHadMatchedTightDetectorIsoMuonSelectorSignal*
    process.IsoMu24eta2p1TightDetectorIsoMatchedGenTauMuTauHadSelectorSignal*
    process.IsoMu24eta2p1TightDetectorIsoEffSignal
    )

#sequence of modules to be run for the Z-->mumu analysis on detector isolated muons
process.tightDetectorIsoZMuMuSequence = cms.Sequence(
    process.genMuMuSelectorZMuMu*process.tightDetectorIsoMuonSelector*
    process.genMuMuMatchedTightDetectorIsoMuonSelectorZMuMu*
    process.IsoMu24eta2p1TightDetectorIsoMatchedGenMuMuSelectorZMuMu*
    process.IsoMu24eta2p1TightDetectorIsoEffZMuMu
    )

#sequence of modules to be run for the signal A1-->tau(-->mu)tau(-->had) analysis on PF isolated
#muons w/o PU subtraction
process.tightPFIsoNoPUSubtractionSignalSequence = cms.Sequence(
    process.genTauMuTauHadSelectorSignal*process.tightPFIsoNoPUSubtractionMuonSelector*
    process.genTauMuTauHadMatchedTightPFIsoNoPUSubtractionMuonSelectorSignal*
    process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenTauMuTauHadSelectorSignal*
    process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffSignal
    )

#sequence of modules to be run for the Z-->mumu analysis on PF isolated muons w/o PU subtraction
process.tightPFIsoNoPUSubtractionZMuMuSequence = cms.Sequence(
    process.genMuMuSelectorZMuMu*process.tightPFIsoNoPUSubtractionMuonSelector*
    process.genMuMuMatchedTightPFIsoNoPUSubtractionMuonSelectorZMuMu*
    process.IsoMu24eta2p1TightPFIsoNoPUSubtractionMatchedGenMuMuSelectorZMuMu*
    process.IsoMu24eta2p1TightPFIsoNoPUSubtractionEffZMuMu
    )

#sequence of modules to be run for the signal A1-->tau(-->mu)tau(-->had) tight muon isolation
#analysis
process.tightMuonIsolationSignalSequence = cms.Sequence(
    process.genTauMuTauHadSelectorSignal*process.tightMuonSelector*
    process.genTauMuTauHadMatchedTightMuonSelectorSignal*process.tightMuonIsoAnalyzerSignal
    )

#sequence of modules to be run for the Z-->mumu tight muon isolation analysis
process.tightMuonIsolationZMuMuSequence = cms.Sequence(
    process.genMuMuSelectorZMuMu*process.tightMuonSelector*
    process.genMuMuMatchedTightMuonSelectorZMuMu*process.tightMuonIsoAnalyzerZMuMu
    )

#path
process.p = cms.Path(process.SEQUENCE)
