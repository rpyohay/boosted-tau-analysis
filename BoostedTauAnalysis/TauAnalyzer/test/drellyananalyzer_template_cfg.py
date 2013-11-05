import FWCore.ParameterSet.Config as cms

process = cms.Process("DRELLYANANALYSIS")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

readFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    skipEvents = cms.untracked.uint32(0)
    )
FILES

## #for L1GtStableParametersRcd and jet corrections
## #START52_V9B is recommended for JEC in Summer12 CMSSWv5.2 MC
## #START52_V9 is what the Summer12 CMSSWv5.2 MC was produced with
## #START53_V15 is the latest recommended JEC tag
## process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
## process.GlobalTag.globaltag = cms.string('START53_V15::All')

#analyze muons
process.DrellYanAnalyzer = cms.EDAnalyzer(
    'DrellYanAnalyzer',
    outFileName = cms.string(
    'OUTFILE'
    ),
    METTag = cms.InputTag('pfMet'),
    muonTag = cms.InputTag('WIsoMuonSelector'),
    PUTag = cms.InputTag('addPileupInfo'),
    vtxTag = cms.InputTag('offlinePrimaryVertices'),
    PUSubtractionCoeff = cms.double(0.5),
    PUScenario = cms.string("PUSCENARIO"),
    MC = cms.bool(MCFLAG),
    minMDimuon = cms.double(0.0),
    minPTMuon = cms.double(0.0)
    )

#path
process.p = cms.Path(process.DrellYanAnalyzer)
