import FWCore.ParameterSet.Config as cms

process = cms.Process("OWNPARTICLES")

#PDG IDs
A_PDGID = 36
Z_PDGID = 23
TAU_PDGID = 15
MU_PDGID = 13

#tau decay types
TAU_HAD = 0
TAU_MU = 1
TAU_E = 2
TAU_ALL = 3

#tau hadronic decay types
TAU_ALL_HAD = -1
TAU_1PRONG_0NEUTRAL = 0
TAU_1PRONG_1NEUTRAL = 1
TAU_1PRONG_2NEUTRAL = 2
TAU_1PRONG_3NEUTRAL = 3
TAU_1PRONG_NNEUTRAL = 4
TAU_2PRONG_0NEUTRAL = 5
TAU_2PRONG_1NEUTRAL = 6
TAU_2PRONG_2NEUTRAL = 7
TAU_2PRONG_3NEUTRAL = 8
TAU_2PRONG_NNEUTRAL = 9
TAU_3PRONG_0NEUTRAL = 10
TAU_3PRONG_1NEUTRAL = 11
TAU_3PRONG_2NEUTRAL = 12
TAU_3PRONG_3NEUTRAL = 13
TAU_3PRONG_NNEUTRAL = 14
TAU_RARE = 15

#no consideration of pT rank
ANY_PT_RANK = -1

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    INPUT_FILES
    ),
    skipEvents = cms.untracked.uint32(0)
    )

#for L1GtStableParametersRcd
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V9::All')

#for HLT selection
process.load('HLTrigger/HLTfilters/hltHighLevel_cfi')

#for rerunning tau reconstruction
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

#define a parameter set to be passed to all modules that utilize GenTauDecayID for A1 decays
commonGenTauDecayIDPSetSignal = cms.PSet(momPDGID = cms.int32(A_PDGID),
                                         chargedHadronPTMin = cms.double(0.0),
                                         neutralHadronPTMin = cms.double(0.0),
                                         chargedLeptonPTMin = cms.double(0.0),
                                         totalPTMin = cms.double(0.0)
                                         )

#define a parameter set to be passed to all modules that utilize GenTauDecayID for Z decays
commonGenTauDecayIDPSetZTauTau = commonGenTauDecayIDPSetSignal.clone()
commonGenTauDecayIDPSetZTauTau.momPDGID = cms.int32(Z_PDGID)

#produce gen tau-->had + tau-->anything collection in A1 decays, only saving the had tau
process.genTauHadXSelectorSignal = cms.EDFilter(
    'GenObjectProducer',
    genParticleTag = cms.InputTag('genParticles'),
    absMatchPDGID = cms.uint32(TAU_PDGID),
    genTauDecayIDPSet = commonGenTauDecayIDPSetSignal,
    primaryTauDecayType = cms.uint32(TAU_HAD),
    sisterTauDecayType = cms.uint32(TAU_ALL),
    primaryTauPTRank = cms.int32(ANY_PT_RANK),
##     primaryTauHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauHadronicDecayType = cms.int32(TAU_1PRONG_0NEUTRAL),
    sisterHadronicDecayType = cms.int32(TAU_ALL_HAD),
    primaryTauAbsEtaMax = cms.double(-1.0),
    primaryTauPTMin = cms.double(-1.0),
    countSister = cms.bool(False),
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    minNumGenObjectsToPassFilter = cms.uint32(1),
    makeAllCollections = cms.bool(False)
    )

#produce gen Z-->tau(-->had)tau(-->had) collection, saving both taus
process.genTauHadTauHadSelectorZTauTau = process.genTauHadXSelectorSignal.clone()
process.genTauHadTauHadSelectorZTauTau.genTauDecayIDPSet = commonGenTauDecayIDPSetZTauTau
process.genTauHadTauHadSelectorZTauTau.sisterTauDecayType = cms.uint32(TAU_HAD)
## process.genTauHadTauHadSelectorZTauTau.primaryTauDecayType = cms.uint32(TAU_ALL)
process.genTauHadTauHadSelectorZTauTau.sisterHadronicDecayType = cms.int32(TAU_1PRONG_0NEUTRAL)
process.genTauHadTauHadSelectorZTauTau.countSister = cms.bool(True)
process.genTauHadTauHadSelectorZTauTau.minNumGenObjectsToPassFilter = cms.uint32(2)

#add a collection of HPS taus in |eta| < 2.1 passing medium isolation discriminator and anti-e/mu
#discriminators (cf.
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2012#Object_ID)
process.medIsoTauHadSelector = cms.EDFilter(
    'CustomTauSelector',
##     tauTag = cms.InputTag('hpsPFTauProducer'),
    baseTauTag = cms.InputTag('hpsPFTauProducer'),
    tauDiscriminatorTags = cms.VInputTag(
    'hpsPFTauDiscriminationByDecayModeFinding',
    'hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr'
    ),
    etaMax = cms.double(2.1),
##     minNumObjsToPassFilter = cms.uint32(0)
    minNumObjsToPassFilter = cms.uint32(2)
    )

#add a collection of HPS taus passing anti-e/mu discriminators in |eta| < 2.1
process.tauHadSelector = process.medIsoTauHadSelector.clone()
process.tauHadSelector.tauDiscriminatorTags = cms.VInputTag(
    'hpsPFTauDiscriminationByDecayModeFinding'
    )

#add a collection of HPS taus in |eta| < 2.1 passing medium isolation discriminator and anti-e/mu
#discriminators matched to gen tau-->had decays from tau-->had + tau-->anything A1 decays
process.genTauHadXMatchedMedIsoTauHadSelectorSignal = cms.EDFilter(
    'GenMatchedTauProducer',
    genParticleTag = cms.InputTag('genParticles'),
    selectedGenParticleTag = cms.InputTag('genTauHadXSelectorSignal'),
    recoObjTag = cms.InputTag('medIsoTauHadSelector'),
    baseRecoObjTag = cms.InputTag('hpsPFTauProducer'),
    genTauDecayIDPSet = commonGenTauDecayIDPSetSignal,
    applyPTCuts = cms.bool(False),
    countKShort = cms.bool(True),
    pTRank = cms.int32(ANY_PT_RANK),
    makeAllCollections = cms.bool(False),
    useGenObjPTRank = cms.bool(True),
    nOutputColls = cms.uint32(1),
    dR = cms.double(0.3),
##     minNumGenObjectsToPassFilter = cms.uint32(0)
    minNumGenObjectsToPassFilter = cms.uint32(1)
    )

#add a collection of HPS taus in |eta| < 2.1 passing medium isolation discriminator and anti-e/mu
#discriminators matched to gen tau-->had decays from Z-->tau(-->had)tau(-->had) decays
process.genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau = process.genTauHadXMatchedMedIsoTauHadSelectorSignal.clone()
process.genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau.selectedGenParticleTag = cms.InputTag(
    'genTauHadTauHadSelectorZTauTau'
    )
process.genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau.genTauDecayIDPSet = commonGenTauDecayIDPSetZTauTau

#add a collection of HPS taus passing anti-e/mu discriminators in |eta| < 2.1 from tau-->had +
#tau-->anything A1 decays
process.genTauHadXMatchedTauHadSelectorSignal = process.genTauHadXMatchedMedIsoTauHadSelectorSignal.clone()
process.genTauHadXMatchedTauHadSelectorSignal.recoObjTag = cms.InputTag('tauHadSelector')

#add a collection of HPS taus passing anti-e/mu discriminators in |eta| < 2.1 from
#Z-->tau(-->had)tau(-->had) decays
process.genTauHadTauHadMatchedTauHadSelectorZTauTau = process.genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau.clone()
process.genTauHadTauHadMatchedTauHadSelectorZTauTau.recoObjTag = cms.InputTag('tauHadSelector')

#produce HPS taus in |eta| < 2.1 passing medium isolation discriminator and anti-e/mu
#discriminators matched to gen tau-->had decays matched to 
#DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 trigger objects from tau-->had + tau-->anything A1
#decays
## process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal = cms.EDProducer('trgMatchedTauProducer', InputProducer = cms.InputTag('genTauHadXMatchedMedIsoTauHadSelectorSignal'), hltTags = cms.VInputTag(cms.InputTag('HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1_v2', '', 'HLT')), isTriggerFilter = cms.untracked.bool(True), HLTSubFilters = cms.untracked.VInputTag('hltDoublePFTau35TrackPt5MediumIsolationProng2Dz02'))
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal = cms.EDProducer('trgMatchedTauProducer', InputProducer = cms.InputTag('genTauHadXMatchedMedIsoTauHadSelectorSignal'), hltTags = cms.VInputTag(cms.InputTag('HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1_v2', '', 'HLT')), isTriggerFilter = cms.untracked.bool(False), HLTSubFilters = cms.untracked.VInputTag('hltDoublePFTau35TrackPt5MediumIsolationProng2Dz02'), checkObjMatching = cms.bool(False))

#produce HPS taus in |eta| < 2.1 passing medium isolation discriminator and anti-e/mu
#discriminators matched to gen tau-->had decays matched to 
#DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 trigger objects from Z-->tau(-->had)tau(-->had) decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadTauHadSelectorZTauTau = process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal.clone()
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadTauHadSelectorZTauTau.InputProducer = cms.InputTag(
    'genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau'
    )

#produce HPS taus in |eta| < 2.1 passing anti-e/mu discriminators matched to gen tau-->had decays
#matched to DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 trigger objects from tau-->had +
#tau-->anything A1 decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadXSelectorSignal = process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal.clone()
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadXSelectorSignal.InputProducer = cms.InputTag('genTauHadXMatchedTauHadSelectorSignal')

#produce HPS taus in |eta| < 2.1 passing anti-e/mu discriminators matched to gen tau-->had decays
#matched to DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 trigger objects from
#Z-->tau(-->had)tau(-->had) decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadTauHadSelectorZTauTau = process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal.clone()
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadTauHadSelectorZTauTau.InputProducer = cms.InputTag('genTauHadTauHadMatchedTauHadSelectorZTauTau')

#compute efficiency of DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 for HPS taus in |eta| < 2.1
#passing medium isolation discriminator and anti-e/mu discriminators matched to gen tau-->had
#decays from tau-->had + tau-->anything A1 decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffSignal = cms.EDAnalyzer(
    'TauEfficiencyAnalyzer',
    outFileName = cms.string('/data1/yohay/NMSSMHiggs_gg_medIso_tauLeg_tauHadX_DoubleMediumIsoPFTau35Trk5eta2p1Prong1.root'),
    denominatorTag = cms.InputTag('genTauHadXMatchedMedIsoTauHadSelectorSignal'),
    numeratorTag = cms.InputTag(
    'DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal'
    ),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    decayModeColors = cms.vuint32(1, 2, 4, 6, 8),
    decayModeStyles = cms.vuint32(20, 21, 22, 23, 24),
    decayModeEntries = cms.vstring('#tau_{#mu}', '#tau_{had}, 1 prong',
                                   '#tau_{had}, 1 prong + 1 #pi^{0}',
                                   '#tau_{had}, 1 prong + 2 #pi^{0}', '#tau_{had}, 3 prong')
    )

#compute efficiency of DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 for HPS taus in |eta| < 2.1
#passing medium isolation discriminator and anti-e/mu discriminators matched to gen
#Z-->tau(-->had)tau(-->had) decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffZTauTau = cms.EDAnalyzer(
    'TauTagAndProbeEfficiencyAnalyzer',
    outFileName = cms.string(
    '/data1/yohay/ZTauTau_medIso_tauHadTauHad_DoubleMediumIsoPFTau35Trk5eta2p1Prong1.root'
    ),
    tagTagTag = cms.InputTag(
    'DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadTauHadSelectorZTauTau'
    ),
    tagFailTag = cms.InputTag('genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau'),
    pTRankColors = cms.vuint32(1, 2, 4, 6),
    pTRankStyles = cms.vuint32(20, 21, 22, 23),
    pTRankEntries = cms.vstring('Highest p_{T}', 'Second highest p_{T}', 'Third highest p_{T}',
                                'Lowest p_{T}'),
    decayModeColors = cms.vuint32(1, 2, 4, 6, 8),
    decayModeStyles = cms.vuint32(20, 21, 22, 23, 24),
    decayModeEntries = cms.vstring('#tau_{#mu}', '#tau_{had}, 1 prong',
                                   '#tau_{had}, 1 prong + 1 #pi^{0}',
                                   '#tau_{had}, 1 prong + 2 #pi^{0}', '#tau_{had}, 3 prong')
    )

#compute efficiency of DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 for HPS taus in |eta| < 2.1
#passing anti-e/mu discriminators matched to gen tau-->had decays from tau--> + tau-->anything A1
#decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffSignal = process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffSignal.clone()
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffSignal.outFileName = cms.string(
    '/data1/yohay/NMSSMHiggs_gg_tauLeg_tauHadX_DoubleMediumIsoPFTau35Trk5eta2p1Prong1.root'
    )
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffSignal.denominatorTag = cms.InputTag(
    'genTauHadXMatchedTauHadSelectorSignal'
    )
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffSignal.numeratorTag = cms.InputTag(
    'DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadXSelectorSignal'
    )

#compute efficiency of DoubleMediumIsoPFTau35_Trk5_eta2p1_Prong1 for HPS taus in |eta| < 2.1
#passing anti-e/mu discriminators matched to gen Z-->tau(-->had)tau(-->had) decays
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffZTauTau = process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffZTauTau.clone()
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffZTauTau.outFileName = cms.string(
    '/data1/yohay/ZTauTau_tauHadTauHad_DoubleMediumIsoPFTau35Trk5eta2p1Prong1.root'
    )
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffZTauTau.tagTagTag = cms.InputTag(
    'DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadTauHadSelectorZTauTau'
    )
process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffZTauTau.tagFailTag = cms.InputTag(
    'genTauHadTauHadMatchedTauHadSelectorZTauTau'
    )

## #analyze tight muon isolation for signal A1-->tau(-->mu)tau(-->had) events
## process.tauHadIsoAnalyzerSignal = cms.EDAnalyzer(
##     'IsolationAnalyzer',
##     outFileName = cms.string('/data1/yohay/isoAnalysis_tightMuons_signal.root'),
##     muonTag = cms.InputTag('genTauHadXMatchedTauHadSelectorSignal'),
##     muonPFIsoPUSubtractionCoeff = cms.double(0.5)
##     )

## #analyze tight muon isolation for Z-->mumu events
## process.tauHadIsoAnalyzerZTauTau = process.tauHadIsoAnalyzerSignal.clone()
## process.tauHadIsoAnalyzerZTauTau.outFileName = cms.string(
##     '/data1/yohay/isoAnalysis_tightMuons_ZTauTau.root')
## process.tauHadIsoAnalyzerZTauTau.muonTag = cms.InputTag('genTauHadTauHadMatchedTauHadSelectorZTauTau')

#sequence of modules to be run for the signal A1-->tau(-->had)tau(-->X) analysis on isolated taus
process.medIsoSignalSequence = cms.Sequence(
    process.genTauHadXSelectorSignal*process.medIsoTauHadSelector*
    process.genTauHadXMatchedMedIsoTauHadSelectorSignal*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadXSelectorSignal*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffSignal
    )

#sequence of modules to be run for the Z-->tau(-->had)tau(-->had) analysis on isolated taus
process.medIsoZTauTauSequence = cms.Sequence(
    process.genTauHadTauHadSelectorZTauTau*process.medIsoTauHadSelector*
    process.genTauHadTauHadMatchedMedIsoTauHadSelectorZTauTau*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadMatchedGenTauHadTauHadSelectorZTauTau*process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1MedIsoTauHadEffZTauTau
    )

#sequence of modules to be run for the signal A1-->tau(-->had)tau(-->X) analysis on all taus
process.signalSequence = cms.Sequence(
    process.genTauHadXSelectorSignal*process.tauHadSelector*
    process.genTauHadXMatchedTauHadSelectorSignal*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadXSelectorSignal*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffSignal
    )

#sequence of modules to be run for the Z-->tau(-->had)tau(-->had) analysis on all taus
process.ZTauTauSequence = cms.Sequence(
    process.genTauHadTauHadSelectorZTauTau*process.tauHadSelector*
    process.genTauHadTauHadMatchedTauHadSelectorZTauTau*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadMatchedGenTauHadTauHadSelectorZTauTau*
    process.DoubleMediumIsoPFTau35Trk5eta2p1Prong1TauHadEffZTauTau
    )

## #sequence of modules to be run for the signal A1-->tau(-->mu)tau(-->had) tight muon isolation
## #analysis
## process.tauHadIsolationSignalSequence = cms.Sequence(
##     process.genTauHadXSelectorSignal*process.tauHadSelector*
##     process.genTauHadXMatchedTauHadSelectorSignal*process.tauHadIsoAnalyzerSignal
##     )

## #sequence of modules to be run for the Z-->mumu tight muon isolation analysis
## process.tauHadIsolationZTauTauSequence = cms.Sequence(
##     process.genTauHadTauHadSelectorZTauTau*process.tauHadSelector*
##     process.genTauHadTauHadMatchedTauHadSelectorZTauTau*process.tauHadIsoAnalyzerZTauTau
##     )

#path
process.p = cms.Path(process.SEQUENCE)
