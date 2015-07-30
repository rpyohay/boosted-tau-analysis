import FWCore.ParameterSet.Config as cms
from subprocess import *

process = cms.Process("EFFANALYSIS")

#user
#USER = Popen("whoami", stdout=PIPE, shell=True).stdout.read().strip('\n')
USER = "yohay"

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

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/' + USER + '/DIR/EDM_files/data_no_selection_aMASS_VERSION.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
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
#process.ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJets")
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

#define a parameter set to be passed to all modules that utilize GenTauDecayID
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

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
                                                    "EFFANALYSIS")
process.recoTauAK5PFJets08Region.jetMuonMapTag = cms.InputTag("CleanJetsTrig", "", "EFFANALYSIS")
process.ak5PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                      "EFFANALYSIS")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu", "EFFANALYSIS")
process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                             "EFFANALYSIS")
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJetsTrig", "ak5PFJetsNoMu",
                                                        "EFFANALYSIS")
process.recoTauCommonSequence = cms.Sequence(process.CleanJetsTrig*
                                             process.ak5PFJetTracksAssociatorAtVertex*
                                             process.recoTauAK5PFJets08Region*
                                             process.recoTauPileUpVertices*
                                             process.pfRecoTauTagInfoProducer
                                             )
process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence)

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
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'EFFANALYSIS'),
    tauHadIsoTag = cms.InputTag(''),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS'), 
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'EFFANALYSIS')
    ),
    jetTag = cms.InputTag('CleanJetsTrig', 'ak5PFJetsNoMu', 'EFFANALYSIS'),
    muonRemovalDecisionTag = cms.InputTag('CleanJetsTrig', '', 'EFFANALYSIS'),
    overlapCandTag = cms.InputTag(''),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(10.0),
    etaMax = cms.double(2.3),
    isoMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(0)
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

#OS filter for tau_mu W_mu charge product
process.OSSFFilterIso = cms.EDFilter('OSSFFilter',
                                     WMuonTag = cms.InputTag('WIsoMuonSelector'),
                                     tauTag = cms.InputTag('muHadIsoTauSelector'),
                                     jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
                                     passFilter = cms.bool(True)
                                     )

#SS filter for tau_mu tau_had charge product
process.SSSFFilterIso = cms.EDFilter('SSSFFilter',
                                     tauTag = cms.InputTag('muHadIsoTauSelector'),
                                     jetMuonMapTag = cms.InputTag('CleanJets', '', 'SKIM'),
                                     passFilter = cms.bool(True)
                                     )

#MT filter
process.load('BoostedTauAnalysis/TauAnalyzer/tauanalyzer_cfi')
process.MTFilter.minMT = cms.double(50.)
process.MTFilter.objTag = cms.InputTag('highestPTWMuonSelector')
process.MTFilter.METTag = cms.InputTag("patType1CorrectedPFMetPFlow")
process.highMTFilter = process.MTFilter.clone()
process.lowMTFilter = process.MTFilter.clone()
process.lowMTFilter.passFilter = cms.bool(False)

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

#produce gen tau mu collection
process.genTauMuSelector = cms.EDFilter(
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
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#produce gen-tau-mu-matched mu+had objects (using the HPS tau)
process.genTauMuMatchedRecoTauSelector = cms.EDFilter(
    'GenMatchedTauProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauMuSelector'),
    recoObjTag = cms.InputTag('muHadIsoTauSelector'),
    baseRecoObjTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(1)
    )

#produce gen-tau-mu-matched mu+had objects (using the HPS tau) passing a CSVM veto
process.genTauMuMatchedRecoTauBVetoSelector = cms.EDFilter(
    'BVetoFilter',
    tauTag = cms.InputTag('genTauMuMatchedRecoTauSelector'),
    #oldJetTag = cms.InputTag('ak5PFJets'),
    oldJetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    jetMuonMapTag = cms.InputTag('CleanJets'),
    bTagInfoTag = cms.InputTag('combinedSecondaryVertexBJetTags', '', 'EFFANALYSIS'),
    CSVMax = cms.double(0.679),
    passFilter = cms.bool(False),
    minNumObjsToPassFilter = cms.uint32(0)
    )

#compute CSVM veto efficiency for gen-tau-mu-matched mu+had objects (using the HPS tau)
process.bVetoEffAnalyzer = cms.EDAnalyzer(
    'TauEfficiencyAnalyzer',
    outFileName = cms.string('b_veto_eff_SAMPLE.root'),
    denominatorTag = cms.InputTag('genTauMuMatchedRecoTauSelector'),
    numeratorTag = cms.InputTag('genTauMuMatchedRecoTauBVetoSelector'),
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

process.bVetoEffHighMTAnalyzer = process.bVetoEffAnalyzer.clone()
process.bVetoEffHighMTAnalyzer.outFileName = cms.string('b_veto_eff_highMT_SAMPLE.root')
process.bVetoEffLowMTAnalyzer = process.bVetoEffAnalyzer.clone()
process.bVetoEffLowMTAnalyzer.outFileName = cms.string('b_veto_eff_lowMT_SAMPLE.root')

#sequences
process.selectionSeqHighMT = cms.Sequence(process.highestPTWMuonSelector*
                                          process.btagging*
                                          process.muHadIsoTauSelector*
                                          process.muonTriggerObjectFilter*
                                          process.OSSFFilterIso*
                                          process.SSSFFilterIso*
                                          process.PFTau*
                                          process.trigMuHadIsoTauSelector*
                                          process.electronSelector*
                                          process.trigMuonEFilter*
                                          process.trigMuonMuFilter*
                                          process.trigMuonTauFilter*
                                          process.highMTFilter*
                                          process.genTauMuSelector*
                                          process.genTauMuMatchedRecoTauSelector*
                                          process.genTauMuMatchedRecoTauBVetoSelector)
process.selectionSeqLowMT = cms.Sequence(process.highestPTWMuonSelector*
                                         process.btagging*
                                         process.muHadIsoTauSelector*
                                         process.muonTriggerObjectFilter*
                                         process.OSSFFilterIso*
                                         process.SSSFFilterIso*
                                         process.PFTau*
                                         process.trigMuHadIsoTauSelector*
                                         process.electronSelector*
                                         process.trigMuonEFilter*
                                         process.trigMuonMuFilter*
                                         process.trigMuonTauFilter*
                                         process.lowMTFilter*
                                         process.genTauMuSelector*
                                         process.genTauMuMatchedRecoTauSelector*
                                         process.genTauMuMatchedRecoTauBVetoSelector)
#path
process.high = cms.Path(process.selectionSeqHighMT*
                     process.bVetoEffHighMTAnalyzer)
process.low = cms.Path(process.selectionSeqLowMT*
                     process.bVetoEffLowMTAnalyzer)
process.schedule = cms.Schedule(process.high,
                                process.low)
