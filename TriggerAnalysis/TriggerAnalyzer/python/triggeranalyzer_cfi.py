import FWCore.ParameterSet.Config as cms

TriggerAnalyzer = cms.EDAnalyzer('TriggerAnalyzer',
                                 #HLTProcessName = cms.untracked.string("HLT"),
                                 #triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", "HLT"),
#                                 unprescaledHLTPaths = cms.untracked.vstring(
#    "HLT_Photon75_CaloIdVL_IsoL_v1", "HLT_DoublePhoton5_IsoVL_CEP_v1", "HLT_DoublePhoton33_v1",
#    "HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", "HLT_Photon26_IsoVL_Photon18_v1",
#    "HLT_Photon26_CaloIdL_IsoVL_Photon18_v1",
#    "HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1",
#    "HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1", "HLT_Photon60_CaloIdL_HT200_v1",
#    "HLT_Photon70_CaloIdL_HT200_v1", "HLT_Photon70_CaloIdL_HT300_v1",
#    "HLT_Photon70_CaloIdL_MHT30_v1", "HLT_Photon70_CaloIdL_MHT50_v1"
#    ),
#                                 outputFile = cms.untracked.string(
#    "/data2/yohay/trigger_analysis_35GeV-30GeV.root"
#    )
)
