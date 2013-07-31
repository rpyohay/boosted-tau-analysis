import FWCore.ParameterSet.Config as cms

CleanJets = cms.EDProducer(
    'CleanJets',
    jetSrc = cms.InputTag("ak5PFJets"),
    outFileName = cms.string('/data1/friccita/NMSSMSignalWH_MuProperties.root'),
    momPDGID = cms.uint32(36),
    genMuTauPTMin = cms.double(0.0), #GeV
    genMuPTMin = cms.double(20.0), #GeV                           
    cutOnGenMatches = cms.bool(True),
    thisTag = cms.InputTag('CleanJets', 'ak5PFJetsNoMu', 'ANALYSIS')
    )
