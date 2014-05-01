import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

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
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    ),
    skipEvents = cms.untracked.uint32(0)
    )

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V7F::All') #MC
## process.GlobalTag.globaltag = cms.string('FT_53_V21_AN4::All') #data

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#for mu-less jets
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("RecoTauTag.RecoTau.RecoTauPiZeroProducer_cfi")
process.load('BoostedTauAnalysis/CleanJets/cleanjets_cfi')

#define a parameter set to be passed to all modules that utilize GenTauDecayID for signal taus
commonGenTauDecayIDPSet = cms.PSet(momPDGID = cms.vint32(A_PDGID),
                                   chargedHadronPTMin = cms.double(0.0),
                                   neutralHadronPTMin = cms.double(0.0),
                                   chargedLeptonPTMin = cms.double(0.0),
                                   totalPTMin = cms.double(0.0))

#define a parameter set for the W-->munu selector
WMuNuPSet = commonGenTauDecayIDPSet.clone()
WMuNuPSet.momPDGID = cms.vint32(W_PDGID)

#define a parameter set for background W+jet jet-parton matching
WRecoilJetPSet = commonGenTauDecayIDPSet.clone()
WRecoilJetPSet.momPDGID = cms.vint32(ANY_PDGID)

#output commands
skimEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_kt6PFJets_rho_*",
    "keep *_TriggerResults_*_*",
    "keep *_addPileupInfo_*_*",
    "keep recoGenParticles_genParticles_*_*",
    "keep recoMuons_muons_*_*",
    "keep recoPFJets_ak5PFJets_*_*",
    "keep *_pfMet_*_*",
    "keep *_WIsoMuonSelector_*_*",
    "keep recoPFJets_CleanJets_ak5PFJetsNoMu_*",
    "keep recoPFJetsrecoPFJetrecoPFJetsrecoPFJetedmrefhelperFindUsingAdvanceedmRefedmValueMap_CleanJets_*_*",
    "keep recoMuonsRefsedmValueMap_CleanJets_*_*",
    "keep booledmValueMap_CleanJets_*_*",
    "keep *_hpsPFTauDiscriminationByDecayModeFinding_*_SKIM",
    "keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr_*_SKIM",
    "keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr_*_SKIM",
    "keep recoPFTaus_hpsPFTauProducer_*_SKIM",
    "keep *_offlinePrimaryVertices_*_*"
    )
    )

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
process.IsoMu24eta2p1Selector.HLTPaths = cms.vstring('HLT_IsoMu24_eta2p1_v1',
                                                     'HLT_IsoMu24_eta2p1_v2',
                                                     'HLT_IsoMu24_eta2p1_v3',
                                                     'HLT_IsoMu24_eta2p1_v4',
                                                     'HLT_IsoMu24_eta2p1_v5',
                                                     'HLT_IsoMu24_eta2p1_v6',
                                                     'HLT_IsoMu24_eta2p1_v7',
                                                     'HLT_IsoMu24_eta2p1_v8',
                                                     'HLT_IsoMu24_eta2p1_v9',
                                                     'HLT_IsoMu24_eta2p1_v10',
                                                     'HLT_IsoMu24_eta2p1_v11',
                                                     'HLT_IsoMu24_eta2p1_v12',
                                                     'HLT_IsoMu24_eta2p1_v13',
                                                     'HLT_IsoMu24_eta2p1_v14',
                                                     'HLT_IsoMu24_eta2p1_v15')
process.IsoMu24eta2p1Selector.throw = cms.bool(False)

#search for a muon with pT > 25 GeV as in WHbb CMS AN-2012/349 and proceed if one can be found
#this will produce a ref to the original muon collection
process.WMuonPTSelector = cms.EDFilter('MuonRefSelector',
                                       src = cms.InputTag('muons'),
                                       cut = cms.string('pt > 25.0'),
                                       filter = cms.bool(True)
                                       )

#search for a loose PF isolated tight muon in |eta| < 2.1 with pT > 25 GeV
#(see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1 for
#isolation definition; CMS AN-2012/349 uses loose isolation working point for WHbb muon selection)
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
                                        minNumObjsToPassFilter = cms.uint32(2)
                                        )

#search for a muon with pT > 5 GeV as in HZZ4l analysis and proceed if one can be found
#this will produce a ref to the original muon collection
process.tauMuonPTSelector = cms.EDFilter('MuonRefSelector',
                                         src = cms.InputTag('muons'),
                                         cut = cms.string('pt > 5.0'),
                                         filter = cms.bool(True)
                                         )

#search for soft muons
#(see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Soft_Muon) not overlapping with
#the W muon in |eta| < 2.4
#this will produce a ref to the original muon collection
process.tauMuonSelector = cms.EDFilter('CustomMuonSelector',
                                       baseMuonTag = cms.InputTag('muons'),
                                       muonTag = cms.InputTag('tauMuonPTSelector'),
                                       vtxTag = cms.InputTag('offlinePrimaryVertices'),
                                       vetoMuonTag = cms.InputTag('WIsoMuonSelector'),
                                       muonID = cms.string('soft'),
                                       PFIsoMax = cms.double(0.2),
                                       detectorIsoMax = cms.double(-1.0),
                                       PUSubtractionCoeff = cms.double(0.5),
                                       usePFIso = cms.bool(True),
                                       passIso = cms.bool(True),
                                       etaMax = cms.double(2.4),
                                       minNumObjsToPassFilter = cms.uint32(1)
                                       )

#clean the jets of soft muons, then rebuild the taus
process.CleanJets.muonSrc = cms.InputTag('tauMuonSelector')
process.CleanJets.PFCandSrc = cms.InputTag('particleFlow')
process.CleanJets.cutOnGenMatches = cms.bool(False)
process.CleanJets.outFileName = cms.string('NMSSMSignal_MuProperties.root')
process.recoTauAK5PFJets08Region.src = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "SKIM")
process.ak5PFJetsRecoTauPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "SKIM")
process.combinatoricRecoTaus.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "SKIM")
process.ak5PFJetTracksAssociatorAtVertex.jets = cms.InputTag("CleanJets", "ak5PFJetsNoMu",
                                                             "SKIM")
process.ak5PFJetsLegacyHPSPiZeros.jetSrc = cms.InputTag("CleanJets", "ak5PFJetsNoMu", "SKIM")
process.recoTauCommonSequence = cms.Sequence(process.CleanJets*
                                             process.ak5PFJetTracksAssociatorAtVertex*
                                             process.recoTauAK5PFJets08Region*
                                             process.recoTauPileUpVertices*
                                             process.pfRecoTauTagInfoProducer
                                             )
process.PFTau = cms.Sequence(process.recoTauCommonSequence*process.recoTauClassicHPSSequence)

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets that pass the medium isolation
#discriminator
#this will produce a ref to the cleaned tau collection
process.muHadIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'SKIM'), 
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    overlapCandTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets
#this will produce a ref to the cleaned tau collection
process.muHadTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    overlapCandTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    pTMin = cms.double(10.0),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets that fail the medium isolation
#discriminator
#this will produce a ref to the cleaned tau collection
process.muHadNonIsoTauSelector = cms.EDFilter(
    'CustomTauSepFromMuonSelector',
    tauTag = cms.InputTag('muHadTauSelector'),
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    overlapCandTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(False),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#output
process.selectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('data_selected.root')
    )
process.antiSelectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('data_anti-selected.root')
    )
process.noSelectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('data_no_selection.root')
    )

#sequences
process.antiSelectionSequence = cms.Sequence(process.IsoMu24eta2p1Selector*process.WMuonPTSelector*
                                             process.WIsoMuonSelector*process.tauMuonPTSelector*
                                             process.tauMuonSelector*process.PFTau*
                                             process.muHadTauSelector*
                                             process.muHadNonIsoTauSelector)
process.selectionSequence = cms.Sequence(process.IsoMu24eta2p1Selector*process.WMuonPTSelector*
                                         process.WIsoMuonSelector*process.tauMuonPTSelector*
                                         process.tauMuonSelector*process.PFTau*
                                         process.muHadIsoTauSelector)
process.noSelectionSequence = cms.Sequence(process.IsoMu24eta2p1Selector*process.WMuonPTSelector*
                                           process.WIsoMuonSelector)

## #selection path
## process.p = cms.Path(process.selectionSequence)
## process.e = cms.EndPath(process.selectedOutput)

#anti-selection path
## process.p = cms.Path(process.antiSelectionSequence)
## process.e = cms.EndPath(process.antiSelectedOutput)

#no selection path
process.p = cms.Path(process.noSelectionSequence)
process.e = cms.EndPath(process.noSelectedOutput)
