import FWCore.ParameterSet.Config as cms

process = cms.Process("MERGE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.GlobalTag.globaltag = "START53_V7F::All"

readFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    FILES
    )
    )

process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

process.output = cms.OutputModule(
    "PoolOutputModule",
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
    "drop *_recoTauAK5PFJets08Region_*_MERGE",
    "drop *_ak5PFJetTracksAssociatorAtVertex_*_MERGE",
    "drop *_ak5PFJetsLegacyHPSPiZeros_*_MERGE",
    "drop *_combinatoricRecoTausDiscriminationByLeadingPionPtCut_*_MERGE",
    "drop *_combinatoricRecoTausHPSSelector_*_MERGE",
    "drop *_hpsSelectionDiscriminator_*_MERGE",
    "drop *_combinatoricRecoTaus_*_MERGE",
    "drop *_hpsPFTauProducerSansRefs_*_MERGE",
    "drop *_pfRecoTauTagInfoProducer_*_MERGE",
    "drop *_recoTauPileUpVertices_*_MERGE",
    "drop *_towerMaker_*_*",
    "drop *_particleFlowClusterHF*_*_*",
    "drop *_pfPhotonTranslator_*_*",
    "drop *_pfElectronTranslator_*_*",
    "drop *_correctedHybridSuperClusters_*_*",
    "drop *_correctedMulti5x5SuperClustersWithPreshower_*_*",
    "drop *_*Voronoi*_*_*"
    ),
    fileName = cms.untracked.string(
    'OUTPUT_FILE_NAME'
    )
    )

## process.p = cms.Path(process.PFTau)

process.end = cms.EndPath(process.output)
