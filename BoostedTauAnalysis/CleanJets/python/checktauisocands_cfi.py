import FWCore.ParameterSet.Config as cms

CheckTauIsoCands = cms.EDAnalyzer('CheckTauIsoCands',
                                  tauSrc = cms.InputTag("hpsPFTauProducer", "", "OWNPARTICLES"),
                                  outFileName = cms.string('/data1/friccita/NMSSM_WH_TauIsoCandidateHistos_alltaus.root'),
                                  momPDGID = cms.uint32(36),
                                  genMuTauPTMin = cms.double(0.0), #GeV
                                  genMuPTMin = cms.double(20.0) #GeV

)
