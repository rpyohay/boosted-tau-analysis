import FWCore.ParameterSet.Config as cms

process = cms.Process("TAUISOCANDS")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')
process.load("BoostedTauAnalysis.CleanJets.checktauisocands_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data1/friccita/NMSSMHiggs_WH_PFTausNoMu.root'    
    )
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)



#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('Anything.root')
#)
  
process.p = cms.Path(process.CheckTauIsoCands)

#process.e = cms.EndPath(process.out)

#processDumpFile = open('PFTausNoMu_isocands.dump', 'w'
#)

#print >> processDumpFile, process.dumpPython()
