import FWCore.ParameterSet.Config as cms

TauAnalyzer = cms.EDAnalyzer(
    'TauAnalyzer',
    outFileName = cms.string(
    '/data1/yohay/NMSSMHiggs_SM-like_analysis_v2.root'
    ),
    genParticleTag = cms.InputTag("genParticles"),
    tauTag = cms.InputTag("hpsPFTauProducer"),
    muonTag = cms.InputTag("muons"),
    vtxTag = cms.InputTag("offlinePrimaryVertices"),
    jetTag = cms.InputTag("ak5PFJets"),
    HPSDiscriminatorTags = cms.VInputTag(
    cms.InputTag("hpsPFTauDiscriminationByDecayModeFinding"),
    cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA")
    ),
##     momPDGID = cms.int32(23), #Z for Drell-Yan
    momPDGID = cms.int32(36), #a for signal
    genMuTauPTMin = cms.double(0.0), #GeV
    genMuPTMin = cms.double(0.0), #GeV
    effVsEtaPTMin = cms.double(10.0), #GeV
    zCut = cms.double(0.1),
    RcutFactor = cms.double(0.5)
    )

OSSFFilter = cms.EDFilter(
    'OSSFFilter',
    WMuonTag = cms.InputTag("WIsoMuonSelector"),
    tauMuonTag = cms.InputTag("tauMuonSelector")
    )

METFilter = cms.EDFilter(
    'METFilter',
    minMET = cms.double(40.),
    METTag = cms.InputTag("pfMet")
    )

MTFilter = cms.EDFilter(
    'PATMTFilter',
    minMT = cms.double(50.),
    passFilter = cms.bool(True),
    METTag = cms.InputTag("patType1CorrectedPFMetPFlow"),
    objTag = cms.InputTag("WIsoMuonSelector")
    )

DelPhiFilter = cms.EDFilter(
    'DelPhiFilter',
    delPhiCut = cms.double(0.5),
    WMuonTag = cms.InputTag("WIsoMuonSelector"),
    tauMuonTag = cms.InputTag("tauMuonSelector")
    )
