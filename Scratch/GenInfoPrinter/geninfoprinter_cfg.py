import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'root://eoscms//eos/cms/store/user/friccita/WToMuNu_skim/Summer12_WJetsToLNu_WMuNuSkim_1000_3_TDu.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

process.demo = cms.EDAnalyzer('GenInfoPrinter')

process.p = cms.Path(process.demo)
