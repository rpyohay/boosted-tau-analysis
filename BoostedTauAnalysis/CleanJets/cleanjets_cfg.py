import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("BoostedTauAnalysis.CleanJets.cleanjets_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/data1/yohay/NMSSMHiggs_gg_skim_1000Files.root'    
    )
)
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
#                                   src = cms.InputTag("genParticles"),                                                                 
#                                   printP4 = cms.untracked.bool(False),
#                                   printPtEtaPhi = cms.untracked.bool(False),
#                                   printVertex = cms.untracked.bool(False),
#                                   printStatus = cms.untracked.bool(True),
#                                   printIndex = cms.untracked.bool(False),
#                                   status = cms.untracked.vint32( 3 )
#                                   )
process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('NMSSMHiggs_ZH_skim_PFJetsNoMu.root')
#)
  
#process.p = cms.Path(process.printTree*process.CleanJets)
process.p = cms.Path(process.CleanJets)
#process.e = cms.EndPath(process.out)
