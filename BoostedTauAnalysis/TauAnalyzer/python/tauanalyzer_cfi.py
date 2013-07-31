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
    effVsEtaPTMin = cms.double(10.0) #GeV
    )
