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
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('BoostedTauAnalysis/CleanJets/cleanjets_cfi')
process.load('BoostedTauAnalysis/TauAnalyzer/tauanalyzer_cfi')
process.load('JetMETCorrections.Type1MET.pfMETCorrections_cff')

#for jet corrections
process.load('JetMETCorrections.Configuration.JetCorrectionServices_cff')

process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("RecoJets.JetAssociationProducers.ak5JTA_cff")
from RecoBTag.SoftLepton.softLepton_cff import *
from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoBTau.JetTagComputer.combinedMVA_cff import *
process.impactParameterTagInfos.jetTracks = cms.InputTag("ak5JetTracksAssociatorAtVertex")
process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM')
process.ak5JetTracksAssociatorAtVertex.tracks = cms.InputTag("generalTracks")
process.btagging = cms.Sequence(
    process.ak5JetTracksAssociatorAtVertex*
    # impact parameters and IP-only algorithms
    process.impactParameterTagInfos*
    (process.trackCountingHighEffBJetTags +
     process.trackCountingHighPurBJetTags +
     process.jetProbabilityBJetTags +
     process.jetBProbabilityBJetTags +
     # SV tag infos depending on IP tag infos, and SV (+IP) based algos
     process.secondaryVertexTagInfos*
     (process.simpleSecondaryVertexHighEffBJetTags +
      process.simpleSecondaryVertexHighPurBJetTags +
      process.combinedSecondaryVertexBJetTags +
      process.combinedSecondaryVertexMVABJetTags) +
     process.ghostTrackVertexTagInfos*
     process.ghostTrackBJetTags)##  +
##     process.softPFMuonsTagInfos*
##     process.softPFMuonBJetTags *
##     process.softPFElectronsTagInfos*
##     process.softPFElectronBJetTags
    )

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

#produce a collection with the highest pT W muon in the event
process.highestPTWMuonSelector = cms.EDFilter('HighestPTMuonRefSelector',
                                              objRefTag = cms.InputTag('WIsoMuonSelector')
                                              )

#clean the jets of trigger muons, then rebuild the taus
process.CleanJetsTrig = process.CleanJets.clone()
process.CleanJetsTrig.muonSrc = cms.InputTag('highestPTWMuonSelector')
process.CleanJetsTrig.PFCandSrc = cms.InputTag('particleFlow')
process.CleanJetsTrig.cutOnGenMatches = cms.bool(False)
process.CleanJetsTrig.outFileName = cms.string(
    'CLEANJETSOUTFILE'
    )
process.recoTauAK5PFJets08Region.src = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                    "MUHADANALYSIS")
process.recoTauAK5PFJets08Region.jetMuonMapTag = cms.InputTag("CleanJetsTrig", "", "MUHADANALYSIS")
process.ak5PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                      "MUHADANALYSIS")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu", "MUHADANALYSIS")
process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                             "MUHADANALYSIS")
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                        "MUHADANALYSIS")
process.recoTauCommonSequence = cms.Sequence(process.CleanJetsTrig*
                                             process.ak5PFJetTracksAssociatorAtVertex*
                                             process.recoTauAK5PFJets08Region*
                                             process.recoTauPileUpVertices*
                                             process.pfRecoTauTagInfoProducer
                                             )
process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence)


#find taus in |eta| < 2.3 matched to muon-tagged cleaned jets that pass the isolation
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
    pTMin = cms.double(20.0),
    etaMax = cms.double(2.3),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#find taus in |eta| < 2.3 matched to the trigger muon-cleaned jets that pass
#decay mode finding and isolation discriminator
process.trigMuHadIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'MUHADANALYSIS'),
    tauHadIsoTag = cms.InputTag(''),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'MUHADANALYSIS'), 
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'MUHADANALYSIS')
    ),
    jetTag = cms.InputTag('CleanJetsTrig', 'ak5PFJetsNoMu', 'MUHADANALYSIS'),
    muonRemovalDecisionTag = cms.InputTag('CleanJetsTrig', '', 'MUHADANALYSIS'),
    overlapCandTag = cms.InputTag(''),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(10.0),
    etaMax = cms.double(2.3),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(0)
    )

#b-tag filter
process.IsoBVetoFilter = cms.EDFilter(
    'BVetoFilter',
    tauTag = cms.InputTag('muHadIsoTauSelector'),
    oldJetTag = cms.InputTag('ak5PFJets'),
    jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
    bTagInfoTag = cms.InputTag('combinedSecondaryVertexBJetTags', '', 'MUHADANALYSIS'),
    CSVMax = cms.double(0.679),
    passFilter = cms.bool(True),
    minNumObjsToPassFilter = cms.uint32(1)
    )


#trigger muon neighbouring lepton filters
process.electronSelector = cms.EDFilter(
    'CustomElectronSelector',
    baseCandidateTag = cms.InputTag("particleFlow"),
    pTMin = cms.double(7.0),
    etaMax = cms.double(2.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )
process.trigMuonEFilter = cms.EDFilter(
    'TriggerMuonElectronFilter',
    muonTag = cms.InputTag("highestPTWMuonSelector"),
    recoObjTag = cms.InputTag("electronSelector"),
    delRMin = cms.double(0.4)
    )
process.trigMuonMuFilter = cms.EDFilter(
    'TriggerMuonMuonFilter',
    muonTag = cms.InputTag("highestPTWMuonSelector"),
    recoObjTag = cms.InputTag("tauMuonSelector"),
    delRMin = cms.double(0.4)
    )
process.trigMuonTauFilter = cms.EDFilter(
    'TriggerMuonTauFilter',
    muonTag = cms.InputTag("highestPTWMuonSelector"),
    recoObjTag = cms.InputTag("trigMuHadIsoTauSelector"),
    delRMin = cms.double(0.4)
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

#analyze isolated taus, *NotSmeared MET (i.e. central value for MET uncertainty)
process.highMTMuHadIsoTauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    'PREFIX_highMT_central_SAMPLE_VERSION.root'
    ),
    tauTag = cms.InputTag('IsoBVetoFilter'),
    METTag = cms.InputTag('patType1CorrectedPFMetNotSmeared'),
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
    bJetTag = cms.InputTag('combinedSecondaryVertexBJetTags', '', 'MUHADANALYSIS'),
    dR = cms.double(0.3),
    tauPTMin = cms.double(20.0), #GeV
    tauDecayMode = cms.int32(TAU_ALL_HAD),
    uncorrJetPTMin = cms.double(0.0), #GeV
    tauArbitrationMethod = cms.string("m"),
    PUScenario = cms.string("PUSCENARIO"),
    zCut = cms.double(0.1),
    RcutFactor = cms.double(0.5),
    CSVMax = cms.double(0.679),
    MC = cms.bool(True),
    higgsReweight = cms.bool(HIGGSREW),
    reweight = cms.bool(True),
    isHighMT = cms.bool(True),
    bTagScaleShift = cms.string("mean"),
    sample = cms.string("SAMPLE"),
    muHadMassBins = cms.vdouble(0.0, 1.0, 2.0, 3.0, 4.0, 11.0),
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
process.lowMTMuHadIsoTauAnalyzer = process.highMTMuHadIsoTauAnalyzer.clone()
process.lowMTMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_central_SAMPLE_VERSION.root'
    )
process.lowMTMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, e/g scale shifted
process.highMTEGScaleDownMuHadIsoTauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaEGScale_SAMPLE_VERSION.root'
    ),
    tauTag = cms.InputTag('IsoBVetoFilter'),
    METTag = cms.InputTag('patType1CorrectedPFMetElectronEnDownNotSmeared'),
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
    bJetTag = cms.InputTag('combinedSecondaryVertexBJetTags', '', 'MUHADANALYSIS'),
    dR = cms.double(0.3),
    tauPTMin = cms.double(20.0), #GeV
    tauDecayMode = cms.int32(TAU_ALL_HAD),
    uncorrJetPTMin = cms.double(0.0), #GeV
    tauArbitrationMethod = cms.string("m"),
    PUScenario = cms.string("PUSCENARIO"),
    zCut = cms.double(0.1),
    RcutFactor = cms.double(0.5),
    CSVMax = cms.double(0.679),
    MC = cms.bool(True),
    higgsReweight = cms.bool(HIGGSREW),
    reweight = cms.bool(True),
    isHighMT = cms.bool(True),
    bTagScaleShift = cms.string("mean"),
    sample = cms.string("SAMPLE"),
    muHadMassBins = cms.vdouble(0.0, 1.0, 2.0, 3.0, 4.0, 11.0),
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
process.lowMTEGScaleDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)
process.highMTEGScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTEGScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaEGScale_SAMPLE_VERSION.root'
    )
process.highMTEGScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetElectronEnUpNotSmeared'
    )
process.lowMTEGScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTEGScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaEGScale_SAMPLE_VERSION.root'
    )
process.lowMTEGScaleUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, JER shifted
process.highMTJERDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJERDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaJER_SAMPLE_VERSION.root'
    )
process.highMTJERDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetResDownNotSmeared'
    )
process.lowMTJERDownMuHadIsoTauAnalyzer = process.highMTJERDownMuHadIsoTauAnalyzer.clone()
process.lowMTJERDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaJER_SAMPLE_VERSION.root'
    )
process.lowMTJERDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)
process.highMTJERUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJERUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaJER_SAMPLE_VERSION.root'
    )
process.highMTJERUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetResUpNotSmeared'
    )
process.lowMTJERUpMuHadIsoTauAnalyzer = process.highMTJERUpMuHadIsoTauAnalyzer.clone()
process.lowMTJERUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaJER_SAMPLE_VERSION.root'
    )
process.lowMTJERUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, JES shifted
process.highMTJESDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJESDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaJES_SAMPLE_VERSION.root'
    )
process.highMTJESDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetEnDownNotSmeared'
    )
process.lowMTJESDownMuHadIsoTauAnalyzer = process.highMTJESDownMuHadIsoTauAnalyzer.clone()
process.lowMTJESDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaJES_SAMPLE_VERSION.root'
    )
process.lowMTJESDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)
process.highMTJESUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTJESUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaJES_SAMPLE_VERSION.root'
    )
process.highMTJESUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetJetEnUpNotSmeared'
    )
process.lowMTJESUpMuHadIsoTauAnalyzer = process.highMTJESUpMuHadIsoTauAnalyzer.clone()
process.lowMTJESUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaJES_SAMPLE_VERSION.root'
    )
process.lowMTJESUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, muon energy scale shifted
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnDownNotSmeared'
    )
process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTMuEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTMuEnergyScaleDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnUpNotSmeared'
    )
process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTMuEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaMuEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTMuEnergyScaleUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, tau energy scale shifted
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnDownNotSmeared'
    )
process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTTauEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTTauEnergyScaleDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnUpNotSmeared'
    )
process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTTauEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaTauEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTTauEnergyScaleUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, unclustered energy scale shifted
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnDownNotSmeared'
    )
process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer = process.highMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTUnclusteredEnergyScaleDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnUpNotSmeared'
    )
process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer = process.highMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaUnclusteredEnergyScale_SAMPLE_VERSION.root'
    )
process.lowMTUnclusteredEnergyScaleUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

#analyze isolated taus, b veto scale shifted
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer = process.highMTEGScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_minus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.METTag = cms.InputTag(
    'patType1CorrectedPFMetPFlow'
    )
process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.bTagScaleShift = cms.string("min")

process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer = process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.clone()
process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_minus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

process.highMTBVetoScaleUpMuHadIsoTauAnalyzer = process.highMTBVetoScaleDownMuHadIsoTauAnalyzer.clone()
process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_highMT_plus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.bTagScaleShift = cms.string("max")
process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer = process.highMTBVetoScaleUpMuHadIsoTauAnalyzer.clone()
process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer.outFileName = cms.string(
    'PREFIX_lowMT_plus1SigmaBVetoScale_SAMPLE_VERSION.root'
    )
process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer.isHighMT = cms.bool(False)

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
process.MTFilter.objTag = cms.InputTag('highestPTWMuonSelector')
process.highMTFilter = process.MTFilter.clone()
process.highMTFilter.METTag = cms.InputTag('patType1CorrectedPFMetNotSmeared')
process.lowMTFilter = process.highMTFilter.clone()
process.lowMTFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, e/g scale shifted
process.highMTEGScaleDownFilter = process.MTFilter.clone()
process.highMTEGScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetElectronEnDownNotSmeared'
    )
process.lowMTEGScaleDownFilter = process.highMTEGScaleDownFilter.clone()
process.lowMTEGScaleDownFilter.passFilter = cms.bool(False)
process.highMTEGScaleUpFilter = process.MTFilter.clone()
process.highMTEGScaleUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetElectronEnUpNotSmeared')
process.lowMTEGScaleUpFilter = process.highMTEGScaleUpFilter.clone()
process.lowMTEGScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, JER shifted
process.highMTJERDownFilter = process.MTFilter.clone()
process.highMTJERDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetResDownNotSmeared')
process.lowMTJERDownFilter = process.highMTJERDownFilter.clone()
process.lowMTJERDownFilter.passFilter = cms.bool(False)
process.highMTJERUpFilter = process.MTFilter.clone()
process.highMTJERUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetResUpNotSmeared')
process.lowMTJERUpFilter = process.highMTJERUpFilter.clone()
process.lowMTJERUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, JES shifted
process.highMTJESDownFilter = process.MTFilter.clone()
process.highMTJESDownFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetEnDownNotSmeared')
process.lowMTJESDownFilter = process.highMTJESDownFilter.clone()
process.lowMTJESDownFilter.passFilter = cms.bool(False)
process.highMTJESUpFilter = process.MTFilter.clone()
process.highMTJESUpFilter.METTag = cms.InputTag('patType1CorrectedPFMetJetEnUpNotSmeared')
process.lowMTJESUpFilter = process.highMTJESUpFilter.clone()
process.lowMTJESUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, muon energy scale shifted
process.highMTMuEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTMuEnergyScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnDownNotSmeared'
    )
process.lowMTMuEnergyScaleDownFilter = process.highMTMuEnergyScaleDownFilter.clone()
process.lowMTMuEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTMuEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTMuEnergyScaleUpFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetMuonEnUpNotSmeared'
    )
process.lowMTMuEnergyScaleUpFilter = process.highMTMuEnergyScaleUpFilter.clone()
process.lowMTMuEnergyScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, tau energy scale shifted
process.highMTTauEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTTauEnergyScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnDownNotSmeared'
    )
process.lowMTTauEnergyScaleDownFilter = process.highMTTauEnergyScaleDownFilter.clone()
process.lowMTTauEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTTauEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTTauEnergyScaleUpFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetTauEnUpNotSmeared'
    )
process.lowMTTauEnergyScaleUpFilter = process.highMTTauEnergyScaleUpFilter.clone()
process.lowMTTauEnergyScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, unclustered energy scale shifted
process.highMTUnclusteredEnergyScaleDownFilter = process.MTFilter.clone()
process.highMTUnclusteredEnergyScaleDownFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnDownNotSmeared'
    )
process.lowMTUnclusteredEnergyScaleDownFilter = process.highMTUnclusteredEnergyScaleDownFilter.clone()
process.lowMTUnclusteredEnergyScaleDownFilter.passFilter = cms.bool(False)
process.highMTUnclusteredEnergyScaleUpFilter = process.MTFilter.clone()
process.highMTUnclusteredEnergyScaleUpFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetUnclusteredEnUpNotSmeared'
    )
process.lowMTUnclusteredEnergyScaleUpFilter = process.highMTUnclusteredEnergyScaleUpFilter.clone()
process.lowMTUnclusteredEnergyScaleUpFilter.passFilter = cms.bool(False)

#MT filter for isolated taus, b veto scale shifted
process.highMTBVetoScaleFilter = process.MTFilter.clone()
process.highMTBVetoScaleFilter.METTag = cms.InputTag(
    'patType1CorrectedPFMetPFlow'
    )
process.lowMTBVetoScaleFilter = process.highMTBVetoScaleFilter.clone()
process.lowMTBVetoScaleFilter.passFilter = cms.bool(False)

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
    process.genTauMuSelector*
    process.highestPTWMuonSelector
    )
process.baseIsoTauAnalysisSequence = cms.Sequence(
    process.muHadIsoTauSelector*
    process.btagging*
    process.IsoBVetoFilter*
    process.corrJetDistinctIsoTauSelector*
    process.muonTriggerObjectFilter*
    process.OSSFFilterIso*
    process.SSSFFilterIso*
    process.PFTau*
    process.trigMuHadIsoTauSelector*
    process.electronSelector*
    process.trigMuonEFilter*
    process.trigMuonMuFilter*
    process.trigMuonTauFilter
    )
process.highMTIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTFilter*
    process.highMTMuHadIsoTauAnalyzer
    )
process.lowMTIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTFilter*
    process.lowMTMuHadIsoTauAnalyzer
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
process.highMTBVetoScaleIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.highMTBVetoScaleFilter*
    process.highMTBVetoScaleDownMuHadIsoTauAnalyzer*
    process.highMTBVetoScaleUpMuHadIsoTauAnalyzer
    )
process.lowMTBVetoScaleIsoTauAnalysis = cms.Path(
    process.baseIsoTauAnalysisSequence*
    process.lowMTBVetoScaleFilter*
    process.lowMTBVetoScaleDownMuHadIsoTauAnalyzer*
    process.lowMTBVetoScaleUpMuHadIsoTauAnalyzer
    )

#path
process.schedule = cms.Schedule(
    process.begin,
    process.highMTIsoTauAnalysis,
    process.lowMTIsoTauAnalysis,
    process.highMTEGScaleDownIsoTauAnalysis,
    process.lowMTEGScaleDownIsoTauAnalysis,
#    process.highMTJERDownIsoTauAnalysis, #v6 WH and v7 ggH skims don't apply JER smearing
#    process.lowMTJERDownIsoTauAnalysis, #v6 WH and v7 ggH skims don't apply JER smearing
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
#    process.highMTJERUpIsoTauAnalysis, #v6 WH and v7 ggH skims don't apply JER smearing
#    process.lowMTJERUpIsoTauAnalysis, #v6 WH and v7 ggH skims don't apply JER smearing
    process.highMTJESUpIsoTauAnalysis,
    process.lowMTJESUpIsoTauAnalysis,
    process.highMTMuEnergyScaleUpIsoTauAnalysis,
    process.lowMTMuEnergyScaleUpIsoTauAnalysis,
    process.highMTTauEnergyScaleUpIsoTauAnalysis,
    process.lowMTTauEnergyScaleUpIsoTauAnalysis,
    process.highMTUnclusteredEnergyScaleUpIsoTauAnalysis,
    process.lowMTUnclusteredEnergyScaleUpIsoTauAnalysis,
    process.highMTBVetoScaleIsoTauAnalysis,
    process.lowMTBVetoScaleIsoTauAnalysis
    )
## process.e = cms.EndPath(process.output)
