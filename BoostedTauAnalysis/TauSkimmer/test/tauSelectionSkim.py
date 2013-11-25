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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring('file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root'),
    skipEvents = cms.untracked.uint32(0)
    )

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V7F::All') #START53_V7A used for production

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
    "keep *",
    "drop *_*ak7*_*_*",
    "drop *_*kt4*_*_*",
    "drop *_kt6GenJets_*_*",
    "drop *_kt6CaloJets*_*_*",
    "drop *_kt6PFJetsCentral*_*_*",
    "drop *_kt6PFJets_sigma*_*",
    "drop *_kt6JetID_*_*",
    "drop *_kt6PFJets_rhos_*",
    "drop recoPFJets_kt6PFJets_*_*",
    "drop *_fixedGridRho*_*_*",
    "drop *_hfEMClusters_*_*",
    "drop *_eid*_*_*",
    "drop *_muonMETValueMapProducer_muCorrData_*",
    "drop *_muons_muonShowerInformation_*",
    "drop *_muons_combined_*",
    "drop *_muons_csc_*",
    "drop *_muons_dt_*",
    "drop l1extraL1HFRingss_*_*_*",
    "drop *_muonsFromCosmics*_*_*",
    "drop recoCaloClusters_*_*_*",
    "drop recoPreshowerCluster*_*_*_*",
    "drop *_hfRecoEcalCandidate_*_*",
    "drop *_generalV0Candidates_*_*",
    "drop *_selectDigi_*_*",
    "drop *_*BJetTags*_*_*",
    "drop *_castorreco_*_*",
    "drop *_reduced*RecHits*_*_*",
    "drop *_PhotonIDProd_*_*",
    "drop *_*_*photons*_*",
    "drop *_dedx*_*_*",
    "drop *_*_cosmicsVeto_*",
    "drop *_muonTCMETValueMapProducer_*_*",
    "drop *_BeamHaloSummary_*_*",
    "drop *_caloRecoTau*_*_*",
    "drop *_GlobalHaloData_*_*",
    "drop *_hpsTancTau*_*_*",
    "drop *_shrinkingConePFTau*_*_*",
    "drop *_ak5CaloJets_*_*",
    "drop *_ak5TrackJets_*_*",
    "drop *_*_uncleanOnly*_*",
    "drop recoCaloMETs_*_*_*",
    "drop recoConversions_*_*_*",
    "drop *_CastorTowerReco_*_*",
    "drop *_uncleanedOnlyGsfElectron*_*_*",
    "drop recoJPTJets_*_*_*",
    "drop recoMETs_*_*_*",
    "drop *_photons_*_*",
    "drop *_photonCore_*_*",
    "drop *_ak5PFJetsRecoTauPiZeros_*_RECO",
    "drop *_ak5PFJetsRecoTauPiZeros_*_HLT",
    "drop *_hpsPFTauDiscrimination*_*_RECO",
    "drop *_hpsPFTauDiscrimination*_*_HLT",
    "drop *_hpsPFTauProducer_*_RECO",
    "drop *_hpsPFTauProducer_*_HLT",
    "drop *_recoTauAK5PFJets08Region_*_SKIM",
    "drop *_ak5PFJetTracksAssociatorAtVertex_*_SKIM",
    "drop *_ak5PFJetsLegacyHPSPiZeros_*_SKIM",
    "drop *_combinatoricRecoTausDiscriminationByLeadingPionPtCut_*_SKIM",
    "drop *_combinatoricRecoTausHPSSelector_*_SKIM",
    "drop *_hpsSelectionDiscriminator_*_SKIM",
    "drop *_combinatoricRecoTaus_*_SKIM",
    "drop *_hpsPFTauProducerSansRefs_*_SKIM",
    "drop *_pfRecoTauTagInfoProducer_*_SKIM",
    "drop *_recoTauPileUpVertices_*_SKIM",
    "drop *_towerMaker_*_*",
    "drop *_particleFlowClusterHF*_*_*",
    "drop *_pfPhotonTranslator_*_*",
    "drop *_pfElectronTranslator_*_*",
    "drop *_correctedHybridSuperClusters_*_*",
    "drop *_correctedMulti5x5SuperClustersWithPreshower_*_*",
    "drop *_*Voronoi*_*_*"
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
process.IsoMu24eta2p1Selector.HLTPaths = cms.vstring('HLT_IsoMu24_eta2p1_v13')

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

#produce AK5 PF jets in |eta| < 2.4
process.jetSelector = cms.EDFilter('PFJetRefSelector',
                                   src = cms.InputTag('ak5PFJets'),
                                   cut = cms.string('abs(eta) < 2.4'),
                                   filter = cms.bool(True)
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
    'CustomTauSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'SKIM'), 
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets
#this will produce a ref to the cleaned tau collection
process.muHadTauSelector = cms.EDFilter(
    'CustomTauSelector',
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    passDiscriminator = cms.bool(True),
    etaMax = cms.double(2.4),
    dR = cms.double(0.5),
    minNumObjsToPassFilter = cms.uint32(1)
    )

#find taus in |eta| < 2.4 matched to muon-tagged cleaned jets that fail the medium isolation
#discriminator
#this will produce a ref to the cleaned tau collection
process.muHadNonIsoTauSelector = cms.EDFilter(
    'CustomTauSelector',
    tauTag = cms.InputTag('muHadTauSelector'),
    baseTauTag = cms.InputTag('hpsPFTauProducer', '', 'SKIM'),
    tauDiscriminatorTags = cms.VInputTag(
    cms.InputTag('hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr', '', 'SKIM')
    ),
    jetTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'SKIM'),
    muonRemovalDecisionTag = cms.InputTag('CleanJets'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
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
    fileName = cms.untracked.string('WNJetsToLNu_selected.root')
    )
process.antiSelectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('WNJetsToLNu_anti-selected.root')
    )
process.noSelectedOutput = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    outputCommands = skimEventContent.outputCommands,
    fileName = cms.untracked.string('WNJetsToLNu_no_selection.root')
    )

#sequences
process.antiSelectionSequence = cms.Sequence(process.genWMuNuSelector*
                                             process.IsoMu24eta2p1Selector*process.WMuonPTSelector*
                                             process.WIsoMuonSelector*process.jetSelector*
                                             process.tauMuonPTSelector*process.tauMuonSelector*
                                             process.PFTau*process.muHadTauSelector*
                                             process.muHadNonIsoTauSelector)
process.selectionSequence = cms.Sequence(process.genWMuNuSelector*process.IsoMu24eta2p1Selector*
                                         process.WMuonPTSelector*process.WIsoMuonSelector*
                                         process.jetSelector*process.tauMuonPTSelector*
                                         process.tauMuonSelector*process.PFTau*
                                         process.muHadIsoTauSelector)
process.noSelectionSequence = cms.Sequence(process.genWMuNuSelector*process.IsoMu24eta2p1Selector*
                                           process.WMuonPTSelector*process.WIsoMuonSelector*
                                           process.jetSelector*process.tauMuonPTSelector*
                                           process.tauMuonSelector*process.PFTau)

## #selection path
## process.p = cms.Path(process.selectionSequence)
## process.e = cms.EndPath(process.selectedOutput)

#anti-selection path
## process.p = cms.Path(process.antiSelectionSequence)
## process.e = cms.EndPath(process.antiSelectedOutput)

#no selection path
process.p = cms.Path(process.noSelectionSequence)
process.e = cms.EndPath(process.noSelectedOutput)
