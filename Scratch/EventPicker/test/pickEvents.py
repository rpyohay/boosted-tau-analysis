import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/data1/yohay/NMSSMHiggs_WH_files1-250_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files251-500_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files501-750_24Sep12.root',
    'file:/data1/yohay/NMSSMHiggs_WH_files751-1000_24Sep12.root'
    ),
    skipEvents = cms.untracked.uint32(0)
    )

process.demo = cms.EDFilter(
    'EventPicker',
    evts = cms.VEventRange(
    EVENTS
    )
    )

process.out = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p')),
    fileName = cms.untracked.string('/data1/yohay/histHIST_xBinXBIN_yBinYBIN.root')
    )

process.p = cms.Path(process.demo)

process.e = cms.EndPath(process.out)
