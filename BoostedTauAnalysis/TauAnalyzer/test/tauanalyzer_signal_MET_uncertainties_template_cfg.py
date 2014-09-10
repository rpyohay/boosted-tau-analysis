import FWCore.ParameterSet.Config as cms

process = cms.Process("MUHADANALYSIS")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
W_PDGID = 24
TAU_PDGID = 15
MU_PDGID = 13
NUMU_PDGID = 14
NUTAU_PDGID = 16
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

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )
FILES

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck') #speed?

#for L1GtStableParametersRcd and jet corrections
#START52_V9B is recommended for JEC in Summer12 CMSSWv5.2 MC
#START52_V9 is what the Summer12 CMSSWv5.2 MC was produced with
#START53_V15 is the latest recommended JEC tag
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V15::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#for mu-less jets
#process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load('BoostedTauAnalysis/TauAnalyzer/tauanalyzer_cfi')

#for jet corrections
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

#define a parameter set to be passed to all modules that utilize GenTauDecayID for signal taus
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

#define a parameter set for the W-->munu selector
WMuNuPSet = commonGenTauDecayIDPSet.clone()
WMuNuPSet.momPDGID = cms.vint32(W_PDGID)

#only proceed if event is a true W-->munu event
process.genWMuNuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(MU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(NUMU_PDGID),
    genTauDecayIDPSet = WMuNuPSet,
    primaryTauDecayType = cms.uint32(TAU_ALL),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(0),
    makeAllCollections = cms.bool(False)
    )

#produce a collection of muons from W-->tau(-->mu)nu
process.genWTauNuSelector = process.genWMuNuSelector.clone()
process.genWTauNuSelector.absMatchPDGIDs = cms.vuint32(TAU_PDGID)
process.genWTauNuSelector.sisterAbsMatchPDGID = cms.uint32(NUTAU_PDGID)
process.genWTauNuSelector.primaryTauDecayType = cms.uint32(TAU_MU)

#produce gen parton collection
process.genPartonSelector = cms.EDFilter(
    'GenPartonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    partonAbsEtaMax = cms.double(-1.0),
    partonPTMin = cms.double(-1.0),
    minNumGenObjectsToPassFilter = cms.uint32(0)
    )

#produce gen muon collection
process.genMuSelector = cms.EDFilter('PdgIdCandViewSelector',
                                     src = cms.InputTag('genParticles'),
                                     pdgId = cms.vint32(-13, 13))

#produce collection of gen muons from a1-->tau(-->mu)tau(-->had) decay
process.genTauMuSelector = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGIDs = cms.vuint32(TAU_PDGID),
    sisterAbsMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    primaryTauDecayType = cms.uint32(TAU_MU),
    sisterTauDecayType = cms.uint32(TAU_HAD),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
    primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(0),
    makeAllCollections = cms.bool(False)
    )

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets that pass the isolation
#discriminator
#this will produce a ref to the cleaned tau collection
process.muHadIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    tauTag = cms.InputTag('muHadTauSelector', '', 'SKIM'),
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr', '',
                                'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'SKIM'), 
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets', '', 'SKIM'),
    overlapCandTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(10.0),
    etaMax = cms.double(2.4),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#b-tag filter
process.IsoBVetoFilter = cms.EDFilter(
    'BVetoFilter',
    tauTag = cms.InputTag('muHadIsoTauSelector'),
    oldJetTag = cms.InputTag('ak5PFJets'),
    jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
    bTagInfoTag = cms.InputTag('combinedSecondaryVertexBJetTags'),
    CSVMax = cms.double(0.679),
    passFilter = cms.bool(True),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#create a collection of corrected jets with pT > 20 GeV and |eta| < 2.4 distinct from the W muon
#and isolated tau
#this collection has no memory of the uncorrected jets
process.corrJetDistinctIsoTauSelector = cms.EDFilter(
    'CustomJetSelector',
    tauTag = cms.InputTag('IsoBVetoFilter'),
    overlapCandTag = cms.InputTag('WIsoMuonSelector'),
    oldJetTag = cms.InputTag('ak5PFJets'),
    jetMuonMapTag = cms.InputTag('CleanJets'),
    pTMin = cms.double(30.0),
    absEtaMax = cms.double(2.4),
    dR = cms.double(0.3),
    minNumObjsToPassFilter = cms.uint32(0),
    maxNumObjsToPassFilter = cms.int32(-1)
    )

#analyze isolated taus, e/g scale shifted
process.highMTEGScaleDownMuHadIsoTauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaEGScale_SAMPLE_VERSION.root'
    ),
    tauTag = cms.InputTag('IsoBVetoFilter'),
    METTag = cms.InputTag('patType1CorrectedPFMetElectronEnDownSmeared'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    muonPFIsoPUSubtractionCoeff = cms.double(0.5),
    genMatchedMuonTag = cms.InputTag('WIsoMuonSelector'),
    oldJetTag = cms.InputTag('ak5PFJets'),
    newJetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
    oldNewJetMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
    genParticleTag = cms.InputTag('genWMuNuSelector'),
    genTauMuTag = cms.InputTag('genTauMuSelector'),
    genWTauMuTag = cms.InputTag('genWTauNuSelector'),
    tauHadIsoTag = cms.InputTag('hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr', '',
                                'SKIM'),
    allMuonTag = cms.InputTag('muons'),
    muonGenParticleTag = cms.InputTag('genMuSelector'),
    PUTag = cms.InputTag('addPileupInfo'),
    vtxTag = cms.InputTag('offlinePrimaryVertices'),
    allGenParticleTag = cms.InputTag('genParticles'),
    corrJetTag = cms.InputTag('corrJetDistinctIsoTauSelector'),
    bJetTag = cms.InputTag('combinedSecondaryVertexBJetTags'),
    dR = cms.double(0.3),
    tauPTMin = cms.double(10.0), #GeV
    tauDecayMode = cms.int32(TAU_ALL_HAD),
    uncorrJetPTMin = cms.double(0.0), #GeV
    tauArbitrationMethod = cms.string("m"),
    PUScenario = cms.string("PUSCENARIO"),
    zCut = cms.double(0.1),
    RcutFactor = cms.double(0.5),
    CSVMax = cms.double(0.679),
    MC = cms.bool(True),
    reweight = cms.bool(False),
    bTagScaleShift = cms.string("mean"),
    sample = cms.string("SAMPLE"),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    hltTags = cms.VInputTag(cms.InputTag("HLT_IsoMu24_eta2p1_v1", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v2", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v3", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v4", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v5", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v6", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v7", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v8", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v9", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v10", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v11", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v12", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v13", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v14", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v15", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_IsoMu24_eta2p1"),
    theRightHLTSubFilter = cms.InputTag("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3cr"),
    HLTSubFilters = cms.untracked.VInputTag("")
    )
process.lowMTEGScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTEGScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaEGScale_SAMPLE_VERSION.root'
    )
process.highMTEGScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTEGScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaEGScale_SAMPLE_VERSION.root'
    )
process.highMTEGScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetElectronEnUpSmeared'
    )
process.lowMTEGScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTEGScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaEGScale_SAMPLE_VERSION.root'
    )

#analyze isolated taus, JER shifted
process.highMTJERDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJERDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaJER_SAMPLE_VERSION.root'
    )
process.highMTJERDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetResDownSmeared'
    )
process.lowMTJERDownMuHadIsoTauAnalyzer = process.highMTJERDownMuHadIsoTauAnalyzer.clone()
process.lowMTJERDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaJER_SAMPLE_VERSION.root'
    )
process.highMTJERUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJERUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaJER_SAMPLE_VERSION.root'
    )
process.highMTJERUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetResUpSmeared'
    )
process.lowMTJERUpMuHadIsoTauAnalyzer = process.highMTJERUpMuHadIsoTauAnalyzer.clone()
process.lowMTJERUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaJER_SAMPLE_VERSION.root'
    )

#analyze isolated taus, JES shifted
process.highMTJESDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJESDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaJES_SAMPLE_VERSION.root'
    )
process.highMTJESDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetEnDownSmeared'
    )
process.lowMTJESDownMuHadIsoTauAnalyzer = process.highMTJESDownMuHadIsoTauAnalyzer.clone()
process.lowMTJESDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaJES_SAMPLE_VERSION.root'
    )
process.highMTJESUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJESUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaJES_SAMPLE_VERSION.root'
    )
process.highMTJESUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetEnUpSmeared'
    )
process.lowMTJESUpMuHadIsoTauAnalyzer = process.highMTJESUpMuHadIsoTauAnalyzer.clone()
process.lowMTJESUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaJES_SAMPLE_VERSION.root'
    )

#analyze isolated taus, muon energy scale shifted
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnDownSmeared'
    )
process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnUpSmeared'
    )
process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )

#analyze isolated taus, tau energy scale shifted
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnDownSmeared'
    )
process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnUpSmeared'
    )
process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )

#analyze isolated taus, unclustered energy scale shifted
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnDownSmeared'
    )
process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )

process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnUpSmeared'
    )
process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )

#analyze isolated taus, b veto scale shifted
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetSmeared'
    )
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.reweight = cms.bool(True)
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.bTagScaleShift = cms.string("min")

process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer = process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.highMTBVetoScaleUpMuHadIsoTauAnalyzer = process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.bTagScaleShift = cms.string("max")
process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer = process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )

#output
process.output = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    fileName = cms.untracked.string(
    'EDMOUTFILE'
    )
    )

#MT filter
process.MTFilter.minMT = cms.double(50.)
process.highMTFilter = process.MTFilter.clone()
process.highMTFilter.METTag = cms.InputTag('patType1CorrectedPFMetSmeared')
process.lowMTFilter = process.highMTFilter.clone()
process.lowMTFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, e/g scale shifted
process.highMTEGScaleDownFilter = process.MTFilter.clone()
process.highMTEGScaleDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetElectronEnDownSmeared')
process.lowMTEGScaleDownFilter = process.highMTEGScaleDownFilter.clone()
process.lowMTEGScaleDownFilter.passFilter = cms.bool(False)
process.highMTEGScaleUpFilter = process.MTFilter.clone()
process.highMTEGScaleUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetElectronEnUpSmeared')
process.lowMTEGScaleUpFilter = process.highMTEGScaleUpFilter.clone()
process.lowMTEGScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, JER shifted
process.highMTJERDownFilter = process.MTFilter.clone()
process.highMTJERDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetResDownSmeared')
process.lowMTJERDownFilter = process.highMTJERDownFilter.clone()
process.lowMTJERDownFilter.passFilter = cms.bool(False)
process.highMTJERUpFilter = process.MTFilter.clone()
process.highMTJERUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetResUpSmeared')
process.lowMTJERUpFilter = process.highMTJERUpFilter.clone()
process.lowMTJERUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, JES shifted
process.highMTJESDownFilter = process.MTFilter.clone()
process.highMTJESDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetEnDownSmeared')
process.lowMTJESDownFilter = process.highMTJESDownFilter.clone()
process.lowMTJESDownFilter.passFilter = cms.bool(False)
process.highMTJESUpFilter = process.MTFilter.clone()
process.highMTJESUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetEnUpSmeared')
process.lowMTJESUpFilter = process.highMTJESUpFilter.clone()
process.lowMTJESUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, muon energy scale shifted
process.highMTMuEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTMuEnergyScaleDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetMuonEnDownSmeared')
process.lowMTMuEnergyScaleDownFilter = process.highMTMuEnergyScaleDownFilter.clone()
process.lowMTMuEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTMuEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTMuEnergyScaleUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetMuonEnUpSmeared')
process.lowMTMuEnergyScaleUpFilter = process.highMTMuEnergyScaleUpFilter.clone()
process.lowMTMuEnergyScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, tau energy scale shifted
process.highMTTauEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTTauEnergyScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnDownSmeared'
    )
process.lowMTTauEnergyScaleDownFilter = process.highMTTauEnergyScaleDownFilter.clone()
process.lowMTTauEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTTauEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTTauEnergyScaleUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetTauEnUpSmeared')
process.lowMTTauEnergyScaleUpFilter = process.highMTTauEnergyScaleUpFilter.clone()
process.lowMTTauEnergyScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, unclustered energy scale shifted
process.highMTUnclusteredEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTUnclusteredEnergyScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnDownSmeared'
    )
process.lowMTUnclusteredEnergyScaleDownFilter = process.highMTUnclusteredEnergyScaleDownFilter.clone()
process.lowMTUnclusteredEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTUnclusteredEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTUnclusteredEnergyScaleUpFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnUpSmeared'
    )
process.lowMTUnclusteredEnergyScaleUpFilter = process.highMTUnclusteredEnergyScaleUpFilter.clone()
process.lowMTUnclusteredEnergyScaleUpFilter.passFilter = cms.bool(False)

#OS filter for tau_mu W_mu charge product
process.OSSFFilterIso = cms.EDFilter('OSSFFilter',
                                     WMuonTag = cms.InputTag('WIsoMuonSelector'),
                                     tauTag = cms.InputTag('IsoBVetoFilter'),
                                     jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
                                     passFilter = cms.bool(True)
                                     )

#SS filter for tau_mu tau_had charge product
process.SSSFFilterIso = cms.EDFilter('SSSFFilter',
                                     tauTag = cms.InputTag('IsoBVetoFilter'),
                                     jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
                                     passFilter = cms.bool(True)
                                     )

#muon trigger object filter
process.muonTriggerObjectFilter = cms.EDFilter(
    'MuonTriggerObjectFilter',
    recoObjTag = cms.InputTag("WIsoMuonSelector"),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "HLT"),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
    triggerDelRMatch = cms.untracked.double(0.1),
    hltTags = cms.VInputTag(cms.InputTag("HLT_IsoMu24_eta2p1_v1", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v2", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v3", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v4", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v5", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v6", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v7", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v8", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v9", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v10", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v11", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v12", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v13", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v14", "", "HLT"),
                            cms.InputTag("HLT_IsoMu24_eta2p1_v15", "", "HLT")
                            ),
    theRightHLTTag = cms.InputTag("HLT_IsoMu24_eta2p1"),
    theRightHLTSubFilter = cms.InputTag("hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3cr"),
    HLTSubFilters = cms.untracked.VInputTag("")
    )

#sequences
process.begin = cms.Path(
    process.genWMuNuSelector*
    process.genWTauNuSelector*
    process.genPartonSelector*
    process.genMuSelector*
    process.genTauMuSelector
    )
process.baseIsoTauAnalysisSequence = cms.Sequence(
    process.muHadIsoTauSelector*
    process.IsoBVetoFilter*
    process.corrJetDistinctIsoTauSelector*
    process.muonTriggerObjectFilter*
    process.OSSFFilterIso*
    process.SSSFFilterIso
    )
process.highMTEGScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTEGScaleDownFilter*
    process.highMTEGScaleDownMuHadIsoTauAnalyzer
    )
process.lowMTEGScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTEGScaleDownFilter*
    process.lowMTEGScaleDownMuHadIsoTauAnalyzer
    )
process.highMTJERDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTJERDownFilter*
    process.highMTJERDownMuHadIsoTauAnalyzer
    )
process.lowMTJERDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTJERDownFilter*
    process.lowMTJERDownMuHadIsoTauAnalyzer
    )
process.highMTJESDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTJESDownFilter*
    process.highMTJESDownMuHadIsoTauAnalyzer
    )
process.lowMTJESDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTJESDownFilter*
    process.lowMTJESDownMuHadIsoTauAnalyzer
    )
process.highMTMuEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTMuEnergyScaleDownFilter*
    process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.lowMTMuEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTMuEnergyScaleDownFilter*
    process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.highMTTauEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTTauEnergyScaleDownFilter*
    process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.lowMTTauEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTTauEnergyScaleDownFilter*
    process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.highMTUnclusteredEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTUnclusteredEnergyScaleDownFilter*
    process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.lowMTUnclusteredEnergyScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTUnclusteredEnergyScaleDownFilter*
    process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer
    )
process.highMTEGScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTEGScaleUpFilter*
    process.highMTEGScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTEGScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTEGScaleUpFilter*
    process.lowMTEGScaleUpMuHadIsoTauAnalyzer
    )
process.highMTJERUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTJERUpFilter*
    process.highMTJERUpMuHadIsoTauAnalyzer
    )
process.lowMTJERUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTJERUpFilter*
    process.lowMTJERUpMuHadIsoTauAnalyzer
    )
process.highMTJESUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTJESUpFilter*
    process.highMTJESUpMuHadIsoTauAnalyzer
    )
process.lowMTJESUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTJESUpFilter*
    process.lowMTJESUpMuHadIsoTauAnalyzer
    )
process.highMTMuEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTMuEnergyScaleUpFilter*
    process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTMuEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTMuEnergyScaleUpFilter*
    process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.highMTTauEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTTauEnergyScaleUpFilter*
    process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTTauEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTTauEnergyScaleUpFilter*
    process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.highMTUnclusteredEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTUnclusteredEnergyScaleUpFilter*
    process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTUnclusteredEnergyScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTUnclusteredEnergyScaleUpFilter*
    process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer
    )
process.highMTBVetoScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTFilter*
    process.highMTBVetoScaleDownMuHadIsoTauAnalyzer
    )
process.lowMTBVetoScaleDownIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTFilter*
    process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer
    )
process.highMTBVetoScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTFilter*
    process.highMTBVetoScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTBVetoScaleUpIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTFilter*
    process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer
    )

#path
process.schedule = cms.Schedule(
    process.begin,
    process.highMTEGScaleDownIsoTauAnalysis,
    process.lowMTEGScaleDownIsoTauAnalysis,
    process.highMTJERDownIsoTauAnalysis,
    process.lowMTJERDownIsoTauAnalysis,
    process.highMTJESDownIsoTauAnalysis,
    process.lowMTJESDownIsoTauAnalysis,
    process.highMTMuEnergyScaleDownIsoTauAnalysis,
    process.lowMTMuEnergyScaleDownIsoTauAnalysis,
    process.highMTTauEnergyScaleDownIsoTauAnalysis,
    process.lowMTTauEnergyScaleDownIsoTauAnalysis,
    process.highMTUnclusteredEnergyScaleDownIsoTauAnalysis,
    process.lowMTUnclusteredEnergyScaleDownIsoTauAnalysis,
    process.highMTEGScaleUpIsoTauAnalysis,
    process.lowMTEGScaleUpIsoTauAnalysis,
    process.highMTJERUpIsoTauAnalysis,
    process.lowMTJERUpIsoTauAnalysis,
    process.highMTJESUpIsoTauAnalysis,
    process.lowMTJESUpIsoTauAnalysis,
    process.highMTMuEnergyScaleUpIsoTauAnalysis,
    process.lowMTMuEnergyScaleUpIsoTauAnalysis,
    process.highMTTauEnergyScaleUpIsoTauAnalysis,
    process.lowMTTauEnergyScaleUpIsoTauAnalysis,
    process.highMTUnclusteredEnergyScaleUpIsoTauAnalysis,
    process.lowMTUnclusteredEnergyScaleUpIsoTauAnalysis,
    process.highMTBVetoScaleDownIsoTauAnalysis,
    process.lowMTBVetoScaleDownIsoTauAnalysis,
    process.highMTBVetoScaleUpIsoTauAnalysis,
    process.lowMTBVetoScaleUpIsoTauAnalysis
    )
## process.e = cms.EndPath(process.output)
