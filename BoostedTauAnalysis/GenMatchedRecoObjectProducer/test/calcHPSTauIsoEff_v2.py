import FWCore.ParameterSet.Config as cms

process = cms.Process("EFFANALYSIS")

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
    'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All') #?

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#for mu-less jets
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('BoostedTauAnalysis/CleanJets/cleanjets_cfi')

#define a parameter set to be passed to all modules that utilize GenTauDecayID
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
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#require event to fire IsoMu24_eta2p1
process.IsoMu24eta2p1Selector = process.hltHighLevel.clone()
process.IsoMu24eta2p1Selector.HLTPaths = cms.vstring('HLT_IsoMu24_eta2p1_v11')

#search for a muon with pT > 25 GeV as in WHbb CMS AN-2012/349 and proceed if one can be found
#this will produce a ref to the original muon collection
process.WMuonPTSelector = cms.EDFilter('MuonRefSelector',
                                       src = cms.InputTag('muons'),
                                       cut = cms.string('pt > 25.0'),
                                       filter = cms.bool(True)
                                       )

#search for a loose PF isolated tight muon in |eta| < 2.1 with pT > 25 GeV
#(see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1 for
#isolation definition;CMS AN-2012/349 uses loose isolation working point for WHbb muon selection)
#this will produce a ref to the original muon collection
process.WIsoMuonSelector = cms.EDFilter('CustomMuonSelector',
                                        baseMuonTag = cms.InputTag('muons'),
                                        muonTag = cms.InputTag('WMuonPTSelector'),
                                        vtxTag = cms.InputTag('offlinePrimaryVertices'),
                                        muonID = cms.string('tight'),
                                        PFIsoMax = cms.double(0.2),
                                        detectorIsoMax = cms.double(-1.0),
                                        PUSubtractionCoeff = cms.double(0.5),
                                        usePFIso = cms.bool(True),
                                        passIso = cms.bool(True),
                                        etaMax = cms.double(2.1),
                                        minNumObjsToPassFilter = cms.uint32(1)
                                        )

#produce gen tau mu collection
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
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#produce muons in |eta| < 2.1
process.recoMuonEtaSelector = cms.EDFilter('MuonRefSelector',
                                           src = cms.InputTag('muons'),
                                           cut = cms.string('abs(eta) < 2.1')
                                           )

#produce muons with pT > 5 GeV
process.recoMuonPTSelector = process.recoMuonEtaSelector.clone()
process.recoMuonPTSelector.cut = cms.string('pt > 5.0')

#produce gen-tau-mu-matched muons in |eta| < 2.1
process.genTauMuMatchedRecoMuonEtaSelector = cms.EDFilter(
    'GenMatchedMuonProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauMuSelector'),
    recoObjTag = cms.InputTag('recoMuonEtaSelector'),
    baseRecoObjTag = cms.InputTag('muons'),
    genTauDecayIDPSet = commonGenTauDecayIDPSet,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
    minNumGenObjectsToPassFilter = cms.uint32(0)
    )

#produce gen-tau-mu-matched muons with pT > 5 GeV
process.genTauMuMatchedRecoMuonPTSelector = process.genTauMuMatchedRecoMuonEtaSelector.clone()
process.genTauMuMatchedRecoMuonPTSelector.recoObjTag = cms.InputTag('recoMuonPTSelector')

#produce muons in |eta| < 2.1 passing pT cut
process.recoMuonPTEtaSelector = cms.EDFilter(
    'MuonRefSelector',
    src = cms.InputTag('muons'),
    cut = cms.string('(pt > 5.0) && (abs(eta) < 2.1)'),
    )

#produce gen-tau-mu-matched muons passing pT and eta cut
process.genTauMuMatchedRecoMuonPTEtaSelector = process.genTauMuMatchedRecoMuonEtaSelector.clone()
process.genTauMuMatchedRecoMuonPTEtaSelector.recoObjTag = cms.InputTag('recoMuonPTEtaSelector')

#produce gen-tau-mu-matched muons passing pT cut, soft ID, and eta cut
process.genTauMuMatchedRecoSoftMuonSelector = cms.EDFilter(
    'CustomMuonSelector',
    baseMuonTag = cms.InputTag('muons'),
    muonTag = cms.InputTag('genTauMuMatchedRecoMuonPTEtaSelector'),
    vtxTag = cms.InputTag('offlinePrimaryVertices'),
    vetoMuonTag = cms.InputTag('WIsoMuonSelector'),
    muonID = cms.string('soft'),
    PFIsoMax = cms.double(0.2),
    detectorIsoMax = cms.double(-1.0),
    PUSubtractionCoeff = cms.double(0.5),
    usePFIso = cms.bool(True),
    passIso = cms.bool(True),
    etaMax = cms.double(-1.0),
    minNumObjsToPassFilter = cms.uint32(0)
    )

#produce gen-tau-mu-matched muons passing pT cut and soft ID
process.genTauMuMatchedRecoPTSoftMuonSelector = process.genTauMuMatchedRecoSoftMuonSelector.clone()
process.genTauMuMatchedRecoPTSoftMuonSelector.muonTag = cms.InputTag(
    'genTauMuMatchedRecoMuonPTSelector'
    )

#produce gen-tau-mu-matched muons passing eta cut and soft ID
process.genTauMuMatchedRecoEtaSoftMuonSelector = process.genTauMuMatchedRecoSoftMuonSelector.clone()
process.genTauMuMatchedRecoEtaSoftMuonSelector.muonTag = cms.InputTag(
    'genTauMuMatchedRecoMuonEtaSelector'
    )

#clean the jets of soft muons passing pT and eta cut, then rebuild the taus
process.CleanJets.muonSrc = cms.InputTag('genTauMuMatchedRecoSoftMuonSelector')
process.CleanJets.PFCandSrc = cms.InputTag('particleFlow')
process.CleanJets.cutOnGenMatches = cms.bool(False)
process.CleanJets.outFileName = cms.string('NMSSMSignal_MuProperties.root')
process.recoTauAK5PFJets08Region.src = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.ak5PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.recoTauAK5PFJets08Region.jetMuonMapTag = cms.InputTag("CleanJets", "", "EFFANALYSIS")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "EFFANALYSIS")
process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJets", "ak5PFJetsNoMu",
                                                             "EFFANALYSIS")
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu",
                                                        "EFFANALYSIS")
process.recoTauCommonSequence = cms.Sequence(process.CleanJets*
                                             process.ak5PFJetTracksAssociatorAtVertex*
                                             process.recoTauAK5PFJets08Region*
                                             process.recoTauPileUpVertices*
                                             process.pfRecoTauTagInfoProducer
                                             )
process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence)

#produced gen-tau-mu-matched HPS taus accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauSelector = cms.EDFilter(
    'CustomTauSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'EFFANALYSIS'),
    tauDiscriminatorTags = cms.VInputTag(),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'EFFANALYSIS'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets', '', 'EFFANALYSIS'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(-1.0),
    etaMax = cms.double(-1.0),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(0)
    )

#produce gen-tau-mu-matched muons passing pT, eta, and soft ID with accompanying HPS taus
process.genTauMuMatchedSoftMuTauSelector = cms.EDFilter(
    'TauMatchedMuonSelector',
    tauTag = cms.InputTag('genTauMuMatchedRecoTauSelector'),
    muonTag = cms.InputTag('genTauMuMatchedRecoSoftMuonSelector'),
    jetMuonMapTag = cms.InputTag('CleanJets'),
    minNumObjsToPassFilter = cms.uint32(0)
    )

#produced gen-tau-mu-matched HPS taus passing pT accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauPTSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTSelector.pTMin = cms.double(10.0)

#produced gen-tau-mu-matched HPS taus passing eta accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauEtaSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauEtaSelector.etaMax = cms.double(2.4)

#produced gen-tau-mu-matched HPS taus passing pT and eta accompanying muons passing pT, eta, and
#soft ID
process.genTauMuMatchedRecoTauPTEtaSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTEtaSelector.pTMin = cms.double(10.0)
process.genTauMuMatchedRecoTauPTEtaSelector.etaMax = cms.double(2.4)

#produced gen-tau-mu-matched HPS taus passing pT and decay mode finding accompanying muons passing
#pT, eta, and soft ID
process.genTauMuMatchedRecoTauPTDMFSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTDMFSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauPTDMFSelector.pTMin = cms.double(10.0)

#produced gen-tau-mu-matched HPS taus passing eta and decay mode finding accompanying muons
#passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauEtaDMFSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauEtaDMFSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauEtaDMFSelector.etaMax = cms.double(2.4)

#produced gen-tau-mu-matched HPS taus passing pT, eta, and decay mode finding accompanying muons
#passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauPTEtaDMFSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTEtaDMFSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauPTEtaDMFSelector.pTMin = cms.double(10.0)
process.genTauMuMatchedRecoTauPTEtaDMFSelector.etaMax = cms.double(2.4)

#produced gen-tau-mu-matched HPS taus passing pT, decay mode finding, and medium isolation
#accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauPTDMFIsoSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTDMFIsoSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS'),
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauPTDMFIsoSelector.pTMin = cms.double(10.0)

#produced gen-tau-mu-matched HPS taus passing eta, decay mode finding, and medium isolation
#accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauEtaDMFIsoSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauEtaDMFIsoSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS'),
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauEtaDMFIsoSelector.etaMax = cms.double(2.4)

#produced gen-tau-mu-matched HPS taus passing pT, eta, decay mode finding, and medium isolation
#accompanying muons passing pT, eta, and soft ID
process.genTauMuMatchedRecoTauPTEtaDMFIsoSelector = process.genTauMuMatchedRecoTauSelector.clone()
process.genTauMuMatchedRecoTauPTEtaDMFIsoSelector.tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'EFFANALYSIS'),
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'EFFANALYSIS')
    )
process.genTauMuMatchedRecoTauPTEtaDMFIsoSelector.pTMin = cms.double(10.0)
process.genTauMuMatchedRecoTauPTEtaDMFIsoSelector.etaMax = cms.double(2.4)

#compute pT efficiency of soft ID for gen-tau-mu-matched muons
process.softMuPTEffAnalyzer = cms.EDAnalyzer(
    'MuonEfficiencyAnalyzer',
    outFileName = cms.string('softMuPTEff.root'),
    denominatorTag = cms.InputTag('genTauMuMatchedRecoMuonEtaSelector'),
    numeratorTag = cms.InputTag('genTauMuMatchedRecoEtaSoftMuonSelector'),
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

#compute eta efficiency of soft ID for gen-tau-mu-matched muons
process.softMuEtaEffAnalyzer = process.softMuPTEffAnalyzer.clone()
process.softMuEtaEffAnalyzer.outFileName = cms.string(
    'softMuEtaEff.root'
    )
process.softMuEtaEffAnalyzer.denominatorTag = cms.InputTag('genTauMuMatchedRecoMuonPTSelector')
process.softMuEtaEffAnalyzer.numeratorTag = cms.InputTag('genTauMuMatchedRecoPTSoftMuonSelector')

#compute efficiency of HPS taus accompanying muons passing pT, eta, and soft ID
process.HPSTauEffAnalyzer = cms.EDAnalyzer(
    'CandidateEfficiencyAnalyzer',
    outFileName = cms.string('HPSTauEff.root'),
    denominatorTag = cms.InputTag('genTauMuMatchedRecoSoftMuonSelector'),
    numeratorTag = cms.InputTag('genTauMuMatchedSoftMuTauSelector'),
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

#compute efficiency of HPS taus passing pT accompanying muons passing pT, eta, and soft ID
process.HPSTauPTEffAnalyzer = process.HPSTauEffAnalyzer.clone()
process.HPSTauPTEffAnalyzer.outFileName = cms.string('HPSTauPTEff.root')
process.HPSTauPTEffAnalyzer.denominatorTag = cms.InputTag('genTauMuMatchedRecoTauSelector')
process.HPSTauPTEffAnalyzer.numeratorTag = cms.InputTag('genTauMuMatchedRecoTauPTSelector')

#compute efficiency of HPS taus passing pT and eta accompanying muons passing pT, eta, and soft ID
process.HPSTauPTEtaEffAnalyzer = process.HPSTauEffAnalyzer.clone()
process.HPSTauPTEtaEffAnalyzer.outFileName = cms.string('HPSTauPTEtaEff.root')
process.HPSTauPTEtaEffAnalyzer.denominatorTag = cms.InputTag('genTauMuMatchedRecoTauPTSelector')
process.HPSTauPTEtaEffAnalyzer.numeratorTag = cms.InputTag('genTauMuMatchedRecoTauPTEtaSelector')

#compute efficiency of HPS taus passing pT, eta, and decay mode finding accompanying muons passing
#pT, eta, and soft ID
process.HPSTauPTEtaDMFEffAnalyzer = process.HPSTauEffAnalyzer.clone()
process.HPSTauPTEtaDMFEffAnalyzer.outFileName = cms.string('HPSTauPTEtaDMFEff.root')
process.HPSTauPTEtaDMFEffAnalyzer.denominatorTag = cms.InputTag(
    'genTauMuMatchedRecoTauPTEtaSelector'
    )
process.HPSTauPTEtaDMFEffAnalyzer.numeratorTag = cms.InputTag(
    'genTauMuMatchedRecoTauPTEtaDMFSelector'
    )

#compute efficiency of HPS taus passing pT, eta, decay mode, and medium isolation finding
#accompanying muons passing pT, eta, and soft ID
process.HPSTauPTEtaDMFIsoEffAnalyzer = process.HPSTauEffAnalyzer.clone()
process.HPSTauPTEtaDMFIsoEffAnalyzer.outFileName = cms.string('HPSTauPTEtaDMFIsoEff.root')
process.HPSTauPTEtaDMFIsoEffAnalyzer.denominatorTag = cms.InputTag(
    'genTauMuMatchedRecoTauPTEtaDMFSelector'
    )
process.HPSTauPTEtaDMFIsoEffAnalyzer.numeratorTag = cms.InputTag(
    'genTauMuMatchedRecoTauPTEtaDMFIsoSelector'
    )

#sequences
process.baseSelectionSeq = cms.Sequence(process.IsoMu24eta2p1Selector*process.WMuonPTSelector*
                                        process.WIsoMuonSelector)
process.softMuCommonSeq = cms.Sequence(process.genTauMuSelector*process.recoMuonEtaSelector*
                                       process.recoMuonPTSelector*
                                       process.recoMuonPTEtaSelector*
                                       process.genTauMuMatchedRecoMuonPTSelector*
                                       process.genTauMuMatchedRecoMuonEtaSelector*
                                       process.genTauMuMatchedRecoMuonPTEtaSelector*
                                       process.genTauMuMatchedRecoSoftMuonSelector*
                                       process.genTauMuMatchedRecoPTSoftMuonSelector*
                                       process.genTauMuMatchedRecoEtaSoftMuonSelector)
process.softMuPTEffSeq = cms.Sequence(process.softMuPTEffAnalyzer)
process.softMuEtaEffSeq = cms.Sequence(process.softMuEtaEffAnalyzer)
process.tauRelEffSeq = cms.Sequence(process.PFTau*process.genTauMuMatchedRecoTauSelector*
                                    process.genTauMuMatchedSoftMuTauSelector*
                                    process.genTauMuMatchedRecoTauPTSelector*
                                    process.genTauMuMatchedRecoTauEtaSelector*
                                    process.genTauMuMatchedRecoTauPTEtaSelector*
                                    process.genTauMuMatchedRecoTauPTDMFSelector*
                                    process.genTauMuMatchedRecoTauEtaDMFSelector*
                                    process.genTauMuMatchedRecoTauPTEtaDMFSelector*
                                    process.genTauMuMatchedRecoTauPTDMFIsoSelector*
                                    process.genTauMuMatchedRecoTauEtaDMFIsoSelector*
                                    process.genTauMuMatchedRecoTauPTEtaDMFIsoSelector*
                                    process.HPSTauEffAnalyzer*process.HPSTauPTEffAnalyzer*
                                    process.HPSTauPTEtaEffAnalyzer*
                                    process.HPSTauPTEtaDMFEffAnalyzer*
                                    process.HPSTauPTEtaDMFIsoEffAnalyzer)

#path
process.p = cms.Path(process.baseSelectionSeq*process.softMuCommonSeq*process.softMuPTEffSeq*
                     process.softMuEtaEffSeq*process.tauRelEffSeq)
