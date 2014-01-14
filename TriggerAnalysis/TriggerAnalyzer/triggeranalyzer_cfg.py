import FWCore.ParameterSet.Config as cms

process = cms.Process("TRIGGERANALYZER")

process.load("FWCore.MessageService.MessageLogger_cfi")

#stuff you need for the trigger for some unknown reason
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')

#global tag: needed for trigger stuff
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START53_V7F::All')

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'root://eoscms//eos/cms/store/user/yohay/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-Summer12_DR53X-PU_S10_START53_V7A-v1-AODSIM_skim_DYEnriched/data_no_selection_100_1_tAb.root'
    )
                            )
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.load('TriggerAnalysis/TriggerAnalyzer/triggeranalyzer_cfi')
process.TriggerAnalyzer.outputFile = cms.untracked.string(
    'dummy.root'
    )
process.TriggerAnalyzer.unprescaledHLTPaths = cms.untracked.vstring(
    'HLT_Photon32_CaloIdL_Photon26_CaloIdL'
    )

process.p = cms.Path(process.TriggerAnalyzer)
